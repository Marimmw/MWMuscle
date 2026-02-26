#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>      // Für std::setw, std::setprecision
#include <filesystem>   // Braucht C++17 (in CMake sicherstellen!)


#include "simpleSimulation/SSMuscle.h"
#include "utils/MWMath.h"
#include "simpleSimulation/SSBody.h"

// In utility.h ganz oben (nach den #includes):

struct MeshParameter
{
    std::string Type;
    MWMath::Point3D PositionRelativeToBody;
    MWMath::Point3D RotationAxis2Parent;
    double RotationAngle2ParentDeg;
    std::vector<double> Params; // z.B. für Zylinder: [radius, height]
    MWMath::Point3D Color;
    int Segments = 16; // Für runde Meshes wie Zylinder oder Ellipsoid
};

struct MWJointParameter
{
    std::string TypeName; // z.B. "Revolute", "Universal", etc.
    std::string Name;
    MWMath::Point3D Position2Parent;
    MWMath::Point3D RotationAxis;
    double RotationAngleDeg;
    MWMath::Point3D Color;
    std::string ParentName;
    std::vector<MeshParameter> MeshParams;

    std::vector<double> JointStateDegrees; // Aktueller Zustand des Joints in Grad
    std::vector<MWMath::Point3D> RotationAxisList; // Liste der Rotationsachsen für mehrere DOFs
    std::vector<double> JointStartEndValues; // Start- und Endwerte für jeden DOF
};

// output structure for all muscle points over whole time casadi simulation for ONE muscle
struct MuscleAllResult {
    std::string muscleName;
    std::vector<std::vector<MWMath::Point3D>> allOptimizedPoints;
};

// structure for body/meshes geometry state
struct GeometryState {
    MWMath::Point3D Position;
    MWMath::RotMatrix3x3 Orientation;
};

// output structure for body states over whole time casadi simulation
struct BodyResult {
    std::string bodyName;
    std::vector<GeometryState> COMStates;
    std::vector<std::vector<GeometryState>> MeshesStates;
};




//namespace fs = std::filesystem;

inline void exportMuscleLog(const std::string& systemName, SSMuscle* muscle, std::vector<std::string> solverResults = {""}, std::string units="m") 
{
    namespace fs = std::filesystem;
    // 0. Sicherheitschecks
    if (!muscle) {
        std::cerr << "[Export] Fehler: Muskel-Pointer ist null." << std::endl;
        return;
    }
    if (muscle->MNodes.empty()) {
        std::cerr << "[Export] Warnung: Muskel " << muscle->Name << " hat keine Nodes." << std::endl;
        return;
    }

    // 1. Dateinamen und Pfad generieren
    // Format: MuscleLog_<SystemName>_<MuscleName>.txt
    std::string filename = "MuscleLog_" + systemName + "_" + muscle->Name + ".txt";
    
    fs::path targetDir = fs::path("..") / "examples" / "results";
    fs::path fullPath = targetDir / filename;

    // 2. Ordner erstellen
    try {
        if (!fs::exists(targetDir)) fs::create_directories(targetDir);
    } catch (const std::exception& e) {
        std::cerr << "[Export] Fehler beim Erstellen des Ordners: " << e.what() << std::endl;
        return;
    }

    // 3. Datei öffnen
    std::ofstream outFile(fullPath);
    if (!outFile.is_open()) {
        std::cerr << "[Export] Fehler: Konnte Datei nicht schreiben: " << fullPath << std::endl;
        return;
    }

    // 4. Simulations-Daten ermitteln (Anzahl Steps)
    // Wir nehmen die Step-Anzahl vom ersten Node (sollte überall gleich sein)
    size_t numSteps = muscle->MNodes[0].MNodeGlobalSteps.size();
    
    if (numSteps == 0) {
        outFile << "Keine Zeitschritte aufgezeichnet.\n";
        outFile.close();
        return;
    }

    // Formatierung
    const int widthDesc = 20; // Breite der Beschreibung
    const int widthVal  = 22; // Breite der Zahlenwerte

    // =========================================================
    // HEADER SCHREIBEN (Step 0   Step 1   Step 2 ...)
    // =========================================================
    outFile << "MUSCLE LOG: " << muscle->Name << " (System: " << systemName << ")\n";
    outFile << std::string(widthDesc + numSteps * widthVal, '=') << "\n";
    
    // --- NEU: MESH MAPPING IN DEN HEADER SCHREIBEN ---
    outFile << "SolerResults: ";
    for (size_t s = 0; s < solverResults.size(); ++s) {
        outFile << s << "=\"" << solverResults[s] << "\" ";
    }
    outFile << "\n";

    outFile << "Meshes: ";
    for (size_t m = 0; m < muscle->meshPtrs.size(); ++m) {
        if (muscle->meshPtrs[m]) {
            outFile << m << "=\"" << muscle->meshPtrs[m]->Name << "\" ";
        }
    }
    outFile << "\n";
    
    outFile << std::left << std::setw(widthDesc) << "Parameter";
    for (size_t s = 0; s < numSteps; ++s) {
        outFile << std::left << std::setw(widthVal) << ("Step " + std::to_string(s));
    }
    outFile << "\n";
    outFile << std::string(widthDesc + numSteps * widthVal, '-') << "\n";

    // =========================================================
    // DATEN SCHREIBEN (Zeilenweise: Node X, Y, Z, Etas...)
    // =========================================================

    for (size_t k = 0; k < muscle->MNodes.size(); ++k) {
        const auto& node = muscle->MNodes[k];
        std::string nodePrefix = "Node " + std::to_string(k);

        // --- A) POSITION (X, Y, Z) ---
        // Lambda zum Schreiben einer Zeile (spart Code-Duplizierung)
        auto writeRow = [&](std::string rowName, auto valueGetter, std::string unit) {
            outFile << std::left << std::setw(widthDesc) << rowName;
            for (size_t s = 0; s < numSteps; ++s) {
                if (s < node.MNodeGlobalSteps.size()) {
                    double val = valueGetter(s);
                    // Zusammenbauen in Stringstream für sauberes Layout
                    std::stringstream ss;
                    ss << std::fixed << std::setprecision(9) << val << unit;
                    outFile << std::left << std::setw(widthVal) << ss.str();
                } else {
                    outFile << std::left << std::setw(widthVal) << "NaN";
                }
            }
            outFile << "\n";
        };

        writeRow(nodePrefix + " Pos X", [&](size_t s){ return node.MNodeGlobalSteps[s].x; }, units);
        writeRow(nodePrefix + " Pos Y", [&](size_t s){ return node.MNodeGlobalSteps[s].y; }, units);
        writeRow(nodePrefix + " Pos Z", [&](size_t s){ return node.MNodeGlobalSteps[s].z; }, units);

        // --- B) ETAS (Falls vorhanden) ---
        // Wir prüfen im ersten Step, wie viele Etas dieser Node hat
        size_t numEtas = 0;
        if (!node.MNodeEtaSteps.empty() && !node.MNodeEtaSteps[0].empty()) {
            numEtas = node.MNodeEtaSteps[0].size();
        }

        for (size_t e = 0; e < numEtas; ++e) {
            
            // 1. Zeile schreiben: ETA
            std::string etaName = nodePrefix + " Eta " + std::to_string(e);
            outFile << std::left << std::setw(widthDesc) << etaName;
            for (size_t s = 0; s < numSteps; ++s) {
                if (s < node.MNodeEtaSteps.size() && e < node.MNodeEtaSteps[s].size()) {
                    double val = node.MNodeEtaSteps[s][e];
                    outFile << std::scientific << std::setprecision(9) 
                            << std::left << std::setw(widthVal) << val;
                } else {
                    outFile << std::setw(widthVal) << "0.0"; 
                }
            }
            outFile << "\n";

            // 2. Zeile schreiben: PHI (direkt unter dem passenden Eta)
            std::string phiName = nodePrefix + " Phi " + std::to_string(e);
            outFile << std::left << std::setw(widthDesc) << phiName;
            for (size_t s = 0; s < numSteps; ++s) {
                // Safety Check: Hat dieser Step Phi-Daten für dieses Mesh?
                if (s < node.MNodePhiSteps.size() && e < node.MNodePhiSteps[s].size()) {
                    double val = node.MNodePhiSteps[s][e];
                    std::stringstream ss;
                    ss << std::fixed << std::setprecision(9) << val << units; // "m" hier anhängen
                    outFile << std::left << std::setw(widthVal) << ss.str();
                } else {
                    outFile << std::setw(widthVal) << "NaN"; 
                }
            }
            outFile << "\n";
        }
        
        // Leerzeile zwischen Nodes für Lesbarkeit
        if (k < muscle->MNodes.size() - 1) outFile << "\n"; 
    }

    outFile.close();
    std::cout << "[Export] Muskel-Log gespeichert: " << fs::absolute(fullPath) << std::endl;
}

inline void exportMuscleLengthLog(const std::vector<SSMuscle*>& muscles, const std::vector<SSJoint*>& joints){
    namespace fs = std::filesystem;
    
    // 1. Pfad und Ordnerstruktur vorbereiten
    fs::path targetDir = fs::path("..") / "examples" / "results";
    std::string filename = "MuscleLength.txt";
    
    // Erstelle das Verzeichnis, falls es noch nicht existiert
    if (!fs::exists(targetDir)) {
        fs::create_directories(targetDir);
    }
    
    fs::path fullPath = targetDir / filename;
    std::ofstream outFile(fullPath);
    
    if (!outFile.is_open()) {
        std::cerr << "Fehler: Konnte Datei nicht oeffnen/erstellen: " << fullPath << std::endl;
        return;
    }

    // 2. Maximale Anzahl an Steps herausfinden (zur Sicherheit)
    size_t numSteps = 0;
    for (auto* j : joints) {
        if (j && j->DoneAngleSteps.size() > numSteps) numSteps = j->DoneAngleSteps.size();
    }
    for (auto* m : muscles) {
        if (m && m->MuscleLengthSteps.size() > numSteps) numSteps = m->MuscleLengthSteps.size();
    }

    if (numSteps == 0) {
        outFile << "Keine Daten vorhanden (0 Steps)." << std::endl;
        outFile.close();
        return;
    }

    // Präzision für alle kommenden Double-Werte auf 5 festlegen
    outFile << std::fixed << std::setprecision(5);

    // ==============================================================================
    // TABELLE 1: JOINT ANGLES
    // ==============================================================================
    outFile << "================================================================================\n";
    outFile << "                               JOINT ANGLES OVER TIME                           \n";
    outFile << "================================================================================\n";
    
    // Header (Steps)
    outFile << std::left << std::setw(25) << "Step / Joint Name";
    for (size_t t = 0; t < numSteps; ++t) {
        outFile << std::right << std::setw(12) << t;
    }
    outFile << "\n--------------------------------------------------------------------------------\n";

    // Gelenke-Daten eintragen
    for (auto* j : joints) {
        if (!j) continue;
        
        outFile << std::left << std::setw(25) << j->Name; // Name auf 25 Zeichen fixieren
        
        for (size_t t = 0; t < numSteps; ++t) {
            if (t < j->DoneAngleSteps.size()) {
                outFile << std::right << std::setw(12) << j->DoneAngleSteps[t];
            } else {
                outFile << std::right << std::setw(12) << "NaN"; // Falls Daten fehlen
            }
        }
        outFile << "\n";
    }
    
    outFile << "\n\n"; // Abstand zwischen den Tabellen

    // ==============================================================================
    // TABELLE 2: MUSCLE LENGTHS
    // ==============================================================================
    outFile << "================================================================================\n";
    outFile << "                              MUSCLE LENGTHS OVER TIME                          \n";
    outFile << "================================================================================\n";
    
    // Header (Steps) - Noch einmal zur besseren Übersicht wie gewünscht
    outFile << std::left << std::setw(25) << "Step / Muscle Name";
    for (size_t t = 0; t < numSteps; ++t) {
        outFile << std::right << std::setw(12) << t;
    }
    outFile << "\n--------------------------------------------------------------------------------\n";

    // Muskel-Daten eintragen
    for (auto* m : muscles) {
        if (!m) continue;
        
        outFile << std::left << std::setw(25) << m->Name; // Name auf 25 Zeichen fixieren
        
        for (size_t t = 0; t < numSteps; ++t) {
            if (t < m->MuscleLengthSteps.size()) {
                outFile << std::right << std::setw(12) << m->MuscleLengthSteps[t];
            } else {
                outFile << std::right << std::setw(12) << "NaN";
            }
        }
        outFile << "\n";
    }

    // 3. Abschluss
    outFile.close();
    std::cout << "Log erfolgreich exportiert nach: " << fullPath << std::endl;
}

void backupSourceCode() {
    namespace fs = std::filesystem;
    // 1. Pfad zur aktuellen Quelldatei (wird vom Compiler gesetzt)
    fs::path sourceFile = __FILE__; 
    
    // 2. Zielordner definieren (relativ zum Ausführungsort / build ordner)
    // Passe diesen Pfad an, falls dein build-Ordner woanders liegt
    fs::path targetDir = "../examples/results";

    // Ordner erstellen, falls nicht existent
    if (!fs::exists(targetDir)) {
        try {
            fs::create_directories(targetDir);
        } catch (const fs::filesystem_error& e) {
            std::cerr << "Fehler beim Erstellen des Backup-Ordners: " << e.what() << std::endl;
            return;
        }
    }

    // 3. Zeitstempel generieren (damit man die Dateien unterscheiden kann)
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d_%H-%M-%S");
    
    // Neuer Dateiname: main_2023-10-27_14-30-00.cpp
    std::string filename = "_main_.cpp";
    fs::path targetPath = targetDir / filename;

    // 4. Kopieren
    try {
        fs::copy_file(sourceFile, targetPath, fs::copy_options::overwrite_existing);
        std::cout << "--> Backup von main.cpp gespeichert unter: " << targetPath << std::endl;
    } catch (const fs::filesystem_error& e) {
        std::cerr << "--> Fehler beim Backup des Source Codes: " << e.what() << std::endl;
    }
}

inline void exportParameterLog(const std::vector<std::vector<double>>& values, 
                        const std::vector<std::vector<std::string>>& descriptions, 
                        const std::string& filename = "parameter_log.txt") 
{
    namespace fs = std::filesystem;
    // 1. Zielpfad definieren
    fs::path targetDir = fs::path("..") / "examples" / "results";
    fs::path fullPath = targetDir / filename;

    // 2. Ordner erstellen
    try {
        if (!fs::exists(targetDir)) fs::create_directories(targetDir);
    } catch (const std::exception& e) {
        std::cerr << "Fehler beim Erstellen des Ordners: " << e.what() << std::endl;
        return;
    }

    // 3. Datei öffnen
    std::ofstream outFile(fullPath);
    if (!outFile.is_open()) {
        std::cerr << "FEHLER: Konnte Datei nicht schreiben: " << fullPath << std::endl;
        return;
    }

    // Sicherheitscheck: Haben wir überhaupt Daten?
    if (values.empty()) {
        outFile << "Keine Daten vorhanden.\n";
        return;
    }

    // Anzahl der Steps und Parameter ermitteln
    size_t numSteps = values.size();
    size_t numParams = values[0].size(); // Annahme: Alle Steps haben gleich viele Params

    // Formatierungseinstellungen
    const int widthDesc = 25; // Breite der Beschreibungs-Spalte
    const int widthVal  = 25; // Breite der Wert-Spalten

    // =========================================================
    // HEADER (Step 0   Step 1   Step 2 ...)
    // =========================================================
    outFile << std::left << std::setw(widthDesc) << "Parameter Description";
    
    for (size_t s = 0; s < numSteps; ++s) {
        std::string header = "Step " + std::to_string(s);
        outFile << std::left << std::setw(widthVal) << header;
    }
    outFile << "\n";

    // Trennlinie
    outFile << std::string(widthDesc + numSteps * widthVal, '-') << "\n";

    // =========================================================
    // DATEN (Zeilenweise pro Parameter)
    // =========================================================
    
    // Wir iterieren über die Parameter-Indizes (Zeilen)
    for (size_t p = 0; p < numParams; ++p) {
        
        // 1. Beschreibung schreiben (aus dem ersten Step nehmen)
        std::string desc = "UNKNOWN";
        if (!descriptions.empty() && p < descriptions[0].size()) {
            desc = descriptions[0][p];
        }
        outFile << std::left << std::setw(widthDesc) << desc;

        // 2. Werte für alle Steps schreiben (Spalten)
        for (size_t s = 0; s < numSteps; ++s) {
            double val = 0.0;
            
            // Safety Check: Hat dieser Step auch diesen Parameter?
            if (p < values[s].size()) {
                val = values[s][p];
                outFile << std::fixed << std::setprecision(9) 
                        << std::left << std::setw(widthVal) << val;
            } else {
                outFile << std::setw(widthVal) << "NaN";
            }
        }
        outFile << "\n"; // Neue Zeile nach jedem Parameter
    }

    outFile.close();
    std::cout << "Parameter erfolgreich gespeichert unter: " << fs::absolute(fullPath) << std::endl;
}