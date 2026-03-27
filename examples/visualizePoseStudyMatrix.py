import os
import re
import pandas as pd
import matplotlib
matplotlib.use('Agg') # Zwingt Matplotlib, ohne grafische Oberfläche zu arbeiten
import matplotlib.pyplot as plt

# ==========================================
# EINSTELLUNGEN (Hier anpassen!)
# ==========================================
INPUT_FILE = "../examples/results/PoseStudy_Summary.txt" # Pfad zu deiner Textdatei
OUTPUT_DIR = "../examples/results/"                          # Speicherort der PNGs

ROWS_PER_TABLE = 4              # Nach wie vielen Zeilen (Muskeln) soll eine neue Tabelle beginnen?
SAVE_SEPARATE_PNGS = False      # True = Jede Tabelle als eigenes PNG / False = Alle Tabellen untereinander in einem PNG

# Farben für die Tabelle
COLOR_GREEN  = '#a8e6cf' # 10/10 Success
COLOR_YELLOW = '#ffd3b6' # Max Iterations (1)
COLOR_RED    = '#ff8b94' # Infeasible (2)
COLOR_GREY   = '#e0e0e0' # Fallback (z.B. andere Fehler)

def parse_pose_study(filepath):
    records = []
    current_record = ""

    # 1. Datei einlesen und Zeilenumbrüche reparieren
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('---') or line.startswith('Pose Name') or line.startswith('Score'):
                continue
            
            if 'CasSys_' in line and current_record:
                records.append(current_record)
                current_record = line
            else:
                if current_record:
                    current_record += " " + line
                else:
                    current_record = line
                    
        if current_record:
            records.append(current_record)

    data = []
    
    # 2. Daten extrahieren
    for rec in records:
        # Pose Name extrahieren (Alles vor CasSys_)
        pose_match = re.match(r'^(.*?)\s+CasSys_', rec)
        if not pose_match:
            continue
        pose = pose_match.group(1).strip()
        
        # Zeilenumbruch in die Pose einfügen (für schmalere Tabellenspalten)
        pose = pose.replace(" (", "\n(")

        # Muskel extrahieren
        muscle_match = re.search(r'(CasSys_\S+)', rec)
        muscle = muscle_match.group(1) if muscle_match else "Unknown"

        # Success extrahieren (z.B. "10/10" oder "9/10")
        succ_match = re.search(r'\s(\d+/\d+)\s', rec)
        success = succ_match.group(1) if succ_match else "0/10"

        # Schritte extrahieren
        steps = re.findall(r'\d\(\d+\)', rec)
        step_str = " ".join(steps)

        # 3. Logik für die Farbgebung
        if re.search(r'\b2\(\d+\)', step_str):
            color = COLOR_RED
        elif re.search(r'\b1\(\d+\)', step_str):
            color = COLOR_YELLOW
        elif success == '10/10':
            color = COLOR_GREEN
        else:
            color = COLOR_GREY

        data.append({
            'Pose': pose,
            'Muscle': muscle.replace("CasSys_", ""), 
            'Success': success,
            'Color': color
        })

    return pd.DataFrame(data)

def create_table_images(df, rows_per_table, separate_pngs, out_dir):
    if df.empty:
        print("Keine Daten zum Verarbeiten gefunden!")
        return

    # === DER FIX: Originale Reihenfolge merken ===
    ordered_poses = df['Pose'].unique()
    ordered_muscles = df['Muscle'].unique()
    # =============================================

    # Pivot-Tabellen generieren (Pandas sortiert hier leider intern alphabetisch)
    df_val = df.pivot(index='Muscle', columns='Pose', values='Success').fillna("")
    df_col = df.pivot(index='Muscle', columns='Pose', values='Color').fillna("#ffffff") 

    # === DER FIX TEIL 2: Die Tabelle wieder in die chronologische Reihenfolge zwingen ===
    df_val = df_val.reindex(index=ordered_muscles, columns=ordered_poses)
    df_col = df_col.reindex(index=ordered_muscles, columns=ordered_poses)
    # ====================================================================================

    num_rows = len(df_val)
    chunks_val = [df_val.iloc[i:i+rows_per_table] for i in range(0, num_rows, rows_per_table)]
    chunks_col = [df_col.iloc[i:i+rows_per_table] for i in range(0, num_rows, rows_per_table)]

    os.makedirs(out_dir, exist_ok=True)
    num_cols = df_val.shape[1]
    
    # Berechne optimale Breite
    fig_width = max(10, num_cols * 2.0 + 3)

    if separate_pngs:
        print(f"Erstelle {len(chunks_val)} separate Bilder...")
        for idx, (c_val, c_col) in enumerate(zip(chunks_val, chunks_col)):
            fig_height = len(c_val) * 0.5 + 2
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
            ax.axis('tight')
            ax.axis('off')
            
            table = ax.table(cellText=c_val.values, rowLabels=c_val.index, colLabels=c_val.columns, 
                             cellColours=c_col.values, loc='center', cellLoc='center')
            
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            table.scale(1, 2.5) 
            
            plt.tight_layout()
            out_file = os.path.join(out_dir, f"PoseStudy_Matrix_Part_{idx+1}.png")
            plt.savefig(out_file, dpi=200, bbox_inches='tight')
            plt.close()
            print(f"Gespeichert: {out_file}")
            
    else:
        print(f"Erstelle 1 kombiniertes Bild mit {len(chunks_val)} Tabellen...")
        fig_height = num_rows * 0.5 + (len(chunks_val) * 2.5)
        fig, axes = plt.subplots(len(chunks_val), 1, figsize=(fig_width, fig_height))
        
        if len(chunks_val) == 1:
            axes = [axes]
            
        for ax, c_val, c_col in zip(axes, chunks_val, chunks_col):
            ax.axis('tight')
            ax.axis('off')
            table = ax.table(cellText=c_val.values, rowLabels=c_val.index, colLabels=c_val.columns, 
                             cellColours=c_col.values, loc='center', cellLoc='center')
            
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            table.scale(1, 2.5)
            
        plt.tight_layout()
        out_file = os.path.join(out_dir, "PoseStudy_Matrix_Combined.png")
        plt.savefig(out_file, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"Gespeichert: {out_file}")


if __name__ == "__main__":
    if not os.path.exists(INPUT_FILE):
        print(f"FEHLER: Die Datei {INPUT_FILE} wurde nicht gefunden.")
    else:
        print("Lese Daten ein...")
        df = parse_pose_study(INPUT_FILE)
        create_table_images(df, ROWS_PER_TABLE, SAVE_SEPARATE_PNGS, OUTPUT_DIR)
        print("Fertig!")