import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def clean_val(val_str):
    """Entfernt 'cm', 'm' und andere Reste und wandelt den Wert in Float um (NaN wird behandelt)."""
    v = val_str.lower().strip()
    if v == 'nan':
        return np.nan
    v = v.replace('cm', '').replace('m', '').replace('c', '')
    try:
        return float(v)
    except ValueError:
        return np.nan

def parse_muscle_log(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    num_steps = 0
    data = {} 
    mesh_names = {}
    solver_results = {}  # NEU
    max_node_idx = 0
    
    for line in lines:
        line = line.rstrip('\n')
        if not line or line.startswith('='):
            continue

        # NEU: SOLVER RESULTS EINLESEN
        if line.startswith("SolerResults:") or line.startswith("SolverResults:"):
            matches = re.findall(r'(\d+)="(.*?)"', line)
            for s_idx_str, s_msg in matches:
                solver_results[int(s_idx_str)] = s_msg
            continue
            
        if line.startswith("Meshes:"):
            matches = re.findall(r'(\d+)="(.*?)"', line)
            for m_idx_str, m_name in matches:
                mesh_names[int(m_idx_str)] = m_name
            continue
        
        if line.startswith("Parameter"):
            num_steps = line.count("Step")
            for s in range(num_steps):
                data[s] = {'eta': {}, 'phi': {}}
            continue
        
        if line.startswith("-") or line.startswith("MUSCLE LOG"):
            continue

        desc = line[:20].strip()
        vals_str = line[20:].split()
        
        if not desc.startswith("Node"):
            continue
            
        desc_parts = desc.split()
        if len(desc_parts) >= 4:
            node_idx = int(desc_parts[1])
            var_type = desc_parts[2]

            # NEU: max_node_idx tracken
            if var_type == "Pos":
                max_node_idx = max(max_node_idx, node_idx)
            
            if var_type in ["Eta", "Phi"]:
                mesh_idx = int(desc_parts[3])
                var_key = var_type.lower()
                
                for s in range(num_steps):
                    if mesh_idx not in data[s][var_key]:
                        data[s][var_key][mesh_idx] = {}
                        
                for s in range(num_steps):
                    if s < len(vals_str):
                        data[s][var_key][mesh_idx][node_idx] = clean_val(vals_str[s])
                    else:
                        data[s][var_key][mesh_idx][node_idx] = np.nan

    # NEU: num_nodes aus Pos-Zeilen
    num_nodes = max_node_idx + 1 if max_node_idx > 0 else 0

    for s in range(num_steps):
        for var_key in ['eta', 'phi']:
            for m_idx in data[s][var_key]:
                node_dict = data[s][var_key][m_idx]
                arr = [node_dict.get(n, np.nan) for n in range(num_nodes)]
                data[s][var_key][m_idx] = arr

    return num_steps, num_nodes, data, mesh_names, solver_results  # NEU: solver_results

def create_pdf_report(input_filepath, output_filepath):
    print(f"Verarbeite: {os.path.basename(input_filepath)}...")
    num_steps, num_nodes, data, mesh_names, solver_results = parse_muscle_log(input_filepath)
    
    if num_steps == 0 or num_nodes <= 2:
        print("  -> Keine oder zu wenige gültige Daten gefunden. Überspringe.")
        return

    inner_nodes = np.arange(1, num_nodes-1)
    
    mesh_colors = ['darkorange', 'royalblue', 'forestgreen', 'firebrick', 'mediumpurple', 
                   'sienna', 'hotpink', 'gray', 'olive', 'cyan']

    with PdfPages(output_filepath) as pdf:
        for step in range(num_steps):
            
            fig, (ax_eta, ax_phi, ax_comp) = plt.subplots(3, 1, figsize=(14, 15), sharex=False)
            
            mesh_indices = sorted(list(data[step]['eta'].keys()))
            num_meshes = len(mesh_indices)
            if num_meshes == 0:
                plt.close(fig)
                continue
            
            bar_width = 0.8 / num_meshes
            
            lines_eta, labels_eta = [], []
            lines_phi, labels_phi = [], []
            
            phi_valid_count = 0
            complementarity_sum = np.zeros(num_nodes)

            msg = solver_results.get(step, "Unknown")
            msg_color = 'green' if msg == "Solve_Succeeded" else 'red'
            
            
            neg_eta_lines = []
            neg_phi_lines = []

            for idx, m_idx in enumerate(mesh_indices):
                m_label = mesh_names.get(m_idx, f"M{m_idx}")
                
                etas = np.array(data[step]['eta'][m_idx])
                phis = np.array(data[step]['phi'][m_idx])

                for n in inner_nodes:
                    eta_val = etas[n]
                    phi_val = phis[n]
                    
                    if not np.isnan(eta_val) and eta_val < 0:
                        neg_eta_lines.append(f"Node {n:>2} | {m_label:<15} | η = {eta_val:.3e}")
                    if not np.isnan(phi_val) and phi_val < 0:
                        neg_phi_lines.append(f"Node {n:>2} | {m_label:<15} | φ = {phi_val:.3e}")
            # Negative Werte als Textbox in ax_eta
            if neg_eta_lines:
                neg_eta_text = "Negative η:\n" + "\n".join(neg_eta_lines)
                ax_eta.text(0.65, 0.95, neg_eta_text,
                            transform=ax_eta.transAxes,
                            fontsize=8, verticalalignment='center',
                            bbox=dict(boxstyle='round,pad=0.4', facecolor="#ededed", 
                                      edgecolor='black', alpha=0.9),
                            fontfamily='monospace', color='black')
            # Negative Werte als Textbox in ax_phi
            if neg_phi_lines:
                neg_phi_text = "Negative φ:\n" + "\n".join(neg_phi_lines)
                ax_phi.text(0.65, 0.95, neg_phi_text,
                            transform=ax_phi.transAxes,
                            fontsize=8, verticalalignment='center',
                            bbox=dict(boxstyle='round,pad=0.4', facecolor="#ededed",
                                      edgecolor='black', alpha=0.9),
                            fontfamily='monospace', color='black')



            for idx, m_idx in enumerate(mesh_indices):
                c = mesh_colors[idx % len(mesh_colors)]
                
                etas = np.array(data[step]['eta'][m_idx])
                phis = np.array(data[step]['phi'][m_idx])
                
                if not np.all(np.isnan(phis)):
                    phi_valid_count += 1

                products = etas * phis
                products = np.nan_to_num(products, nan=0.0)
                complementarity_sum += products
                
                offset = -0.4 + idx * bar_width + bar_width / 2.0
                
                # === NEU: Namen für die Legende auflösen ===
                m_label = mesh_names.get(m_idx, f"M{m_idx}") # Fallback auf M0, falls kein Name gefunden wurde
                
                b1 = ax_eta.bar(inner_nodes + offset, etas[1:-1], width=bar_width*0.9, color=c, 
                                label=rf'$\eta$ {m_label}')
                
                b2 = ax_phi.bar(inner_nodes + offset, phis[1:-1], width=bar_width*0.9, 
                                facecolor='none', edgecolor=c, linewidth=1.5, hatch='////',
                                label=rf'$\phi$ {m_label}')
                
                lines_eta.append(b1[0])
                labels_eta.append(b1.get_label())
                lines_phi.append(b2[0])
                labels_phi.append(b2.get_label())

            l_comp, = ax_comp.plot(inner_nodes, complementarity_sum[1:-1], color=msg_color, linewidth=2.5, 
                                   linestyle='-', marker='d', markersize=8, markerfacecolor=msg_color,
                                   zorder=10, label=r'Kompl. Fehler ($\sum \phi \cdot \eta$)')

            if phi_valid_count == 0 and step == 0:
                print(f"  [Warnung] Alle Phi-Werte im Step {step} sind NaN! Überprüfe den C++ Export.")

            for ax in [ax_eta, ax_phi, ax_comp]:
                ax.set_xticks(inner_nodes)
                ax.set_xticklabels([f"Node {n}" for n in inner_nodes], rotation=45, ha='center')
                # xlim von 0.5 bis 23.5 damit Node 1 und Node 23 gut sichtbar sind
                ax.set_xlim(0.5, num_nodes - 1.5) 
                ax.grid(True, axis='y', alpha=0.3)

            
            

            ax_eta.set_title(rf"Step {step} - Kraft ($\eta$)     [{msg}]",
                             fontsize=16, fontweight='bold', pad=10, color='black')
            ax_eta.set_title(rf"Step {step} - Kraft ($\eta$)", 
                 fontsize=16, fontweight='bold', pad=25)
            ax_eta.text(0.5, 1.01, f"[{msg}]", transform=ax_eta.transAxes,
                        ha='center', va='bottom', fontsize=13,
                        color=msg_color, fontweight='bold')
            ax_eta.set_ylabel(r"$\eta$ (Kräfte)", fontsize=12, color='black')
            ax_eta.axhline(0, color='black', linestyle='-', linewidth=1.2, alpha=0.8) 
            
            ax_phi.set_title(rf"Step {step} - Distanz ($\phi$)", 
                 fontsize=16, fontweight='bold', pad=25, color='dimgray')
            ax_phi.text(0.5, 1.01, f"[{msg}]", transform=ax_phi.transAxes,
                        ha='center', va='bottom', fontsize=13,
                        color=msg_color, fontweight='bold')
            ax_phi.set_ylabel(r"$\phi$ in m (Abstand)", fontsize=12, color='black')
            ax_phi.axhline(0, color='red', linestyle='--', linewidth=1.2, alpha=0.8) 
            
            ax_comp.set_title(rf"Step {step} - Komplementaritätsfehler (muss 0 sein!)", 
                  fontsize=16, fontweight='bold', pad=25, color='black')
            ax_comp.text(0.5, 1.01, f"[{msg}]", transform=ax_comp.transAxes,
                        ha='center', va='bottom', fontsize=13,
                        color=msg_color, fontweight='bold')
            ax_comp.set_ylabel(r"$\sum \phi \cdot \eta$", fontsize=12, color='black')
            ax_comp.axhline(0, color='black', linestyle='-', linewidth=1.5, alpha=0.8) 

            # Legenden
            ax_eta.legend(lines_eta, labels_eta, loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0.)
            ax_phi.legend(lines_phi, labels_phi, loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0.)
            ax_comp.legend([l_comp], [l_comp.get_label()], loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0.)

            plt.tight_layout()
            fig.subplots_adjust(right=0.72, hspace=0.4)# fig.subplots_adjust(right=0.82, hspace=0.4) 
            pdf.savefig(fig)
            plt.close(fig)
            
    print(f"  -> PDF gespeichert unter: {output_filepath}")

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.abspath(os.path.join(script_dir, "../../examples/results"))
    
    search_pattern = os.path.join(results_dir, "MuscleLog_*.txt")
    log_files = glob.glob(search_pattern)
    
    if not log_files:
        print(f"Keine MuscleLog-Dateien gefunden im Ordner: {results_dir}")
        return
        
    print(f"Gefundene Log-Dateien: {len(log_files)}")
    
    for filepath in log_files:
        filename = os.path.basename(filepath)
        pdf_filename = filename.replace('.txt', '.pdf')
        out_filepath = os.path.join(results_dir, pdf_filename)
        
        create_pdf_report(filepath, out_filepath)

if __name__ == "__main__":
    main()