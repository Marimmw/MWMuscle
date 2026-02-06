import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider
import re
import os

# ==============================================================================
# KONFIGURATION
# ==============================================================================

# Pfad zur Datei anpassen!
FILE_PATH = "./results/CasSys_PushedMuscleInputParametersELLIPSOID_PUSH_THROUGH.txt"

# Zuordnung von Mesh-Namen zu Dimensionen (Radien A, B, C)
# Falls ein Name im File auftaucht, der hier nicht steht, wird (0.1, 0.1, 0.1) genutzt.
MESH_CONFIG = {
    "Mesh_Left":   (0.2, 0.2, 0.2),    # Anker Links (Kugel)
    "Mesh_Right":  (0.2, 0.2, 0.2),    # Anker Rechts (Kugel)
    "Mesh_Pusher": (0.3, 0.21, 0.39),  # Der Pusher (Ellipsoid/Torus Hülle)
}

# Farben
COLOR_MESH = 'skyblue'
COLOR_ORIGIN = 'green'
COLOR_INSERTION = 'red'

# ==============================================================================
# PARSER LOGIK
# ==============================================================================

def parse_simulation_file(filepath):
    """Liest das txt file und extrahiert Daten pro Step."""
    if not os.path.exists(filepath):
        print(f"FEHLER: Datei nicht gefunden: {filepath}")
        return None

    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Finde Header Zeile
        header_row = 0
        for i, line in enumerate(lines):
            if "Step 0" in line:
                header_row = i
                break
        
        data_start = header_row + 2 
        parsed_data = {} 
        
        for line in lines[data_start:]:
            if len(line.strip()) < 5: continue
            
            # Name ist die ersten 60 Zeichen (oder bis zum ersten Leerzeichen danach)
            # In deinem File sind es keine Leerzeichen im Namen, also split() sicher.
            parts = line.split()
            if not parts: continue

            # Der erste Teil ist der Parametername (z.B. Mesh_Left_PosX)
            param_name = parts[0].strip()
            
            # Der Rest sind die Werte
            values_str = parts[1:]
            
            values = []
            for v in values_str:
                try:
                    values.append(float(v))
                except:
                    values.append(np.nan)
            
            if param_name:
                parsed_data[param_name] = values

        return parsed_data
        
    except Exception as e:
        print(f"Fehler beim Lesen der Datei: {e}")
        return None

def organize_data_per_step(raw_data):
    """Wandelt die flache Liste in strukturierte Objekte pro Step um."""
    
    # Anzahl Steps ermitteln
    if not raw_data: return []
    first_key = next(iter(raw_data))
    num_steps = len(raw_data[first_key])
    
    steps = []

    # 1. Mesh Namen finden (Suche nach Keys die auf "_PosX" enden)
    mesh_names = set()
    for key in raw_data.keys():
        if key.endswith("_PosX"):
            # "Mesh_Left_PosX" -> "Mesh_Left"
            name = key.replace("_PosX", "").strip()
            mesh_names.add(name)
            
    sorted_mesh_names = sorted(list(mesh_names))
    print(f"Gefundene Meshes: {sorted_mesh_names}")

    # 2. Origin/Insertion Prefix finden (z.B. "PushedMuscle")
    muscle_prefix = None
    for key in raw_data.keys():
        if "OriginX" in key:
            # "PushedMuscle_OriginX" -> "PushedMuscle"
            muscle_prefix = key.replace("_OriginX", "")
            break

    for t in range(num_steps):
        step_info = {
            'meshes': [],
            'origin': None,
            'insertion': None
        }
        
        # A. Meshes extrahieren
        for m_name in sorted_mesh_names:
            try:
                # Position (Suffix _PosX, _PosY, _PosZ)
                if f"{m_name}_PosX" not in raw_data: continue

                pos = np.array([
                    raw_data[f"{m_name}_PosX"][t],
                    raw_data[f"{m_name}_PosY"][t],
                    raw_data[f"{m_name}_PosZ"][t]
                ])
                
                # Rotation Matrix (Suffix _Ori_00, _Ori_01, etc.)
                rot = np.eye(3)
                for r in range(3):
                    for c in range(3):
                        key = f"{m_name}_Ori_{r}{c}" # z.B. Mesh_Left_Ori_00
                        if key in raw_data:
                            rot[r, c] = raw_data[key][t]
                
                step_info['meshes'].append({'name': m_name, 'pos': pos, 'rot': rot})
            except Exception as e:
                print(f"Warnung bei Mesh {m_name}: {e}")
                continue 

        # B. Origin / Insertion extrahieren
        if muscle_prefix:
            try:
                ox = raw_data[f"{muscle_prefix}_OriginX"][t]
                oy = raw_data[f"{muscle_prefix}_OriginY"][t]
                oz = raw_data[f"{muscle_prefix}_OriginZ"][t]
                step_info['origin'] = np.array([ox, oy, oz])

                ix = raw_data[f"{muscle_prefix}_InsertionX"][t]
                iy = raw_data[f"{muscle_prefix}_InsertionY"][t]
                iz = raw_data[f"{muscle_prefix}_InsertionZ"][t]
                step_info['insertion'] = np.array([ix, iy, iz])
            except:
                pass

        steps.append(step_info)
        
    return steps

# ==============================================================================
# 3D MATHEMATIK & PLOTTING
# ==============================================================================

def get_ellipsoid_data(center, rotation_matrix, radii, resolution=20):
    """Erzeugt Grid-Daten für ein rotiertes Ellipsoid."""
    u = np.linspace(0, 2 * np.pi, resolution)
    v = np.linspace(0, np.pi, resolution)
    
    x = radii[0] * np.outer(np.cos(u), np.sin(v))
    y = radii[1] * np.outer(np.sin(u), np.sin(v))
    z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
    
    points = np.stack([x.flatten(), y.flatten(), z.flatten()])
    transformed = rotation_matrix @ points
    transformed = transformed.T + center
    
    x_rot = transformed[:, 0].reshape(x.shape)
    y_rot = transformed[:, 1].reshape(y.shape)
    z_rot = transformed[:, 2].reshape(z.shape)
    
    return x_rot, y_rot, z_rot

def set_axes_equal(ax):
    """Erzwingt gleiche Skalierung der Achsen."""
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print(f"Lese Datei: {FILE_PATH} ...")
    raw_data = parse_simulation_file(FILE_PATH)
    
    if not raw_data:
        print("Keine Daten gefunden.")
        return

    sim_steps = organize_data_per_step(raw_data)
    num_steps = len(sim_steps)
    print(f"{num_steps} Zeitschritte verarbeitet.")
    
    if num_steps == 0:
        print("Keine Steps extrahiert. Überprüfe das Dateiformat.")
        return

    # --- PLOT SETUP ---
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')
    plt.subplots_adjust(bottom=0.15) 

    ax_slider = plt.axes([0.2, 0.05, 0.60, 0.03])
    slider = Slider(ax_slider, 'Step', 0, num_steps - 1, valinit=0, valstep=1)

    def update(val):
        step_idx = int(slider.val)
        step_data = sim_steps[step_idx]
        
        ax.clear()
        ax.set_title(f"Simulation Step {step_idx}")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")

        # 1. Meshes Plotten
        for mesh in step_data['meshes']:
            name = mesh['name']
            
            # Dimensionen aus Config holen oder Default
            dims = MESH_CONFIG.get(name, (0.1, 0.1, 0.1))

            X, Y, Z = get_ellipsoid_data(mesh['pos'], mesh['rot'], dims)
            
            # Surface Plot
            ax.plot_surface(X, Y, Z, color=COLOR_MESH, alpha=0.4, edgecolor='k', linewidth=0.1)
            
            # Koordinatenkreuz pro Mesh (Lokal)
            c = mesh['pos']
            R = mesh['rot']
            l = max(dims) * 1.2
            ax.quiver(c[0], c[1], c[2], R[0,0], R[1,0], R[2,0], length=l, color='r') # X
            ax.quiver(c[0], c[1], c[2], R[0,1], R[1,1], R[2,1], length=l, color='g') # Y
            ax.quiver(c[0], c[1], c[2], R[0,2], R[1,2], R[2,2], length=l, color='b') # Z
            
            # Text Label
            ax.text(c[0], c[1], c[2]+l, name, fontsize=8)

        # 2. Origin & Insertion
        if step_data['origin'] is not None:
            o = step_data['origin']
            ax.scatter(o[0], o[1], o[2], color=COLOR_ORIGIN, s=80, label='Origin', marker='o')
        
        if step_data['insertion'] is not None:
            ins = step_data['insertion']
            ax.scatter(ins[0], ins[1], ins[2], color=COLOR_INSERTION, s=80, label='Insertion', marker='^')

        # Legende (nur wenn Artists da sind)
        if step_data['origin'] is not None or step_data['insertion'] is not None:
            ax.legend()
            
        set_axes_equal(ax)
        fig.canvas.draw_idle()

    slider.on_changed(update)
    update(0)
    plt.show()

if __name__ == "__main__":
    main()