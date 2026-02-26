import os
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

def parse_log_file(filepath):
    """
    Parst die exportierte TXT-Datei und liefert zwei Dictionaries zurück:
    joints = {'JointName': [val0, val1, ...], ...}
    muscles = {'MuscleName': [val0, val1, ...], ...}
    """
    joints = {}
    muscles = {}
    
    current_mode = None
    
    with open(filepath, 'r') as f:
        for line in f:
            # Modus umschalten anhand der Header
            if "JOINT ANGLES OVER TIME" in line:
                current_mode = "joints_wait"
                continue
            elif "MUSCLE LENGTHS OVER TIME" in line:
                current_mode = "muscles_wait"
                continue
                
            if current_mode == "joints_wait" and "Step / Joint Name" in line:
                current_mode = "joints"
                continue
            elif current_mode == "muscles_wait" and "Step / Muscle Name" in line:
                current_mode = "muscles"
                continue
                
            # Datenzeilen parsen
            if current_mode in ["joints", "muscles"]:
                # Separatoren und leere Zeilen überspringen
                if line.startswith("===") or line.startswith("---") or line.strip() == "":
                    continue
                
                # Wir wissen aus dem C++ Code: Die ersten 25 Zeichen sind der Name!
                name_part = line[:25].strip()
                data_part = line[25:].strip()
                
                if not name_part:
                    continue
                    
                # Werte parsen (behandelt auch 'NaN')
                values = []
                for val_str in data_part.split():
                    try:
                        values.append(float(val_str))
                    except ValueError:
                        values.append(np.nan)
                        
                if current_mode == "joints":
                    joints[name_part] = values
                elif current_mode == "muscles":
                    muscles[name_part] = values
                    
    return joints, muscles

def create_pdf_report(filepath, out_filepath):
    """
    Erstellt die PDF mit den gewünschten Plots (A4 Format).
    """
    joints, muscles = parse_log_file(filepath)
    
    if not joints and not muscles:
        print(f"WARNUNG: Keine gültigen Daten in {os.path.basename(filepath)} gefunden.")
        return

    # A4 Größe in Zoll (Portrait/Hochkant)
    a4_figsize = (8.27, 11.69)
    
    with PdfPages(out_filepath) as pdf:
        
        # =========================================================
        # SEITE 1: JOINT ANGLES (1 großes Diagramm)
        # =========================================================
        if joints:
            fig, ax = plt.subplots(figsize=a4_figsize)
            
            # Farbpalette für die Gelenke (sorgt für gute Kontraste)
            cmap = plt.get_cmap("tab20")
            colors = cmap(np.linspace(0, 1, len(joints)))
            
            for (j_name, j_vals), color in zip(joints.items(), colors):
                steps = range(len(j_vals))
                ax.plot(steps, j_vals, label=j_name, color=color, linewidth=2, marker='o', markersize=4)
                
            # Y-Achse hart auf +100 bis -100 festlegen
            ax.set_ylim(-100, 100)
            
            # Styling
            ax.set_title("Joint Angles Over Time", fontsize=16, fontweight='bold', pad=15)
            ax.set_xlabel("Simulation Steps", fontsize=12)
            ax.set_ylabel("Angle [°]", fontsize=12)
            ax.grid(True, linestyle='--', alpha=0.7)
            
            # Legende (außerhalb oder gut platziert)
            ax.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0), fontsize='small', ncol=2)
            
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
            
        # =========================================================
        # SEITE 2+: MUSCLE LENGTHS (Max 3 pro Seite)
        # =========================================================
        if muscles:
            muscle_items = list(muscles.items())
            num_muscles = len(muscle_items)
            
            # In 3er-Blöcke schneiden
            for i in range(0, num_muscles, 3):
                chunk = muscle_items[i:i+3]
                
                # 3 Subplots untereinander auf einer A4 Seite
                fig, axes = plt.subplots(nrows=3, ncols=1, figsize=a4_figsize)
                
                # Falls es unerwarteterweise nur 1 Achse ist (bei nrows=1), in Liste wandeln
                if not hasattr(axes, '__iter__'):
                    axes = [axes]
                
                # Die bis zu 3 Muskeln zeichnen
                for j, (m_name, m_vals) in enumerate(chunk):
                    ax = axes[j]
                    steps = range(len(m_vals))
                    
                    # Rote Farbe für Muskeln, um sich von Gelenken abzuheben
                    ax.plot(steps, m_vals, color='firebrick', linewidth=2.5, marker='s', markersize=5)
                    
                    ax.set_title(f"Muscle Length: {m_name}", fontsize=14, fontweight='bold')
                    ax.set_xlabel("Simulation Steps", fontsize=11)
                    ax.set_ylabel("Length", fontsize=11)
                    ax.grid(True, linestyle='--', alpha=0.7)
                    
                    # Optional: Y-Achse dynamisch aber mit etwas Puffer
                    if any(not np.isnan(v) for v in m_vals):
                        min_val = np.nanmin(m_vals)
                        max_val = np.nanmax(m_vals)
                        padding = (max_val - min_val) * 0.1 if max_val != min_val else 1.0
                        ax.set_ylim(min_val - padding, max_val + padding)
                
                # Wenn auf der letzten Seite weniger als 3 Muskeln sind, leere Plots verstecken
                for j in range(len(chunk), 3):
                    axes[j].set_visible(False)
                
                fig.tight_layout(pad=3.0) # Etwas mehr Platz zwischen den Plots
                pdf.savefig(fig)
                plt.close(fig)

    print(f"Erfolgreich gespeichert: {os.path.basename(out_filepath)}")


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Der Ordner, in den C++ die Ergebnisse schreibt
    results_dir = os.path.abspath(os.path.join(script_dir, "../../examples/results"))
    
    # Suche nach der Datei, die wir vorhin im C++ Code definiert haben (oder ähnlichen)
    search_pattern = os.path.join(results_dir, "*MuscleLength*.txt")
    log_files = glob.glob(search_pattern)
    
    if not log_files:
        print(f"Keine MuscleLengthLog-Dateien gefunden im Ordner: {results_dir}")
        return
        
    print(f"Gefundene Log-Dateien: {len(log_files)}")
    
    for filepath in log_files:
        filename = os.path.basename(filepath)
        pdf_filename = filename.replace('.txt', '.pdf')
        out_filepath = os.path.join(results_dir, pdf_filename)
        
        create_pdf_report(filepath, out_filepath)

if __name__ == "__main__":
    main()