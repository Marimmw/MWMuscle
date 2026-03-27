import sys
import os
import requests

# === FAUbox Einstellungen ===


# ============================

def upload_to_faubox(file_path):
    if not os.path.exists(file_path):
        print(f"FEHLER: Datei {file_path} existiert nicht!")
        sys.exit(1)

    file_name = os.path.basename(file_path)
    upload_url = BASE_WEBDAV_URL + file_name

    print(f"Lade {file_name} in die FAUbox hoch...")

    try:
        with open(file_path, 'rb') as f:
            # PUT-Request mit Basis-Authentifizierung lädt die Datei hoch
            response = requests.put(upload_url, data=f, auth=(USERNAME, PASSWORD))
            
        if response.status_code in (200, 201, 204):
            print("Erfolgreich in die FAUbox hochgeladen!")
        else:
            print(f"Fehler beim Upload: HTTP {response.status_code} - {response.text}")
    except Exception as e:
        print(f"Ausnahmefehler beim Upload: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Verwendung: python sendReport.py <Dateipfad>")
        sys.exit(1)

    file_arg = sys.argv[1]
    upload_to_faubox(file_arg)