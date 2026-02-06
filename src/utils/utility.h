#pragma once

#include <string>
#include <vector>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>      // Für std::setw, std::setprecision
#include <filesystem>   // Braucht C++17 (in CMake sicherstellen!)

#include "simpleSimulation/SSMuscle.h"
#include "utils/MWMath.h"


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







namespace fs = std::filesystem;

inline void exportMuscleLog(const std::string& systemName, SSMuscle* muscle) 
{
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
    const int widthVal  = 15; // Breite der Zahlenwerte

    // =========================================================
    // HEADER SCHREIBEN (Step 0   Step 1   Step 2 ...)
    // =========================================================
    outFile << "MUSCLE LOG: " << muscle->Name << " (System: " << systemName << ")\n";
    outFile << std::string(widthDesc + numSteps * widthVal, '=') << "\n";
    
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
        auto writeRow = [&](std::string rowName, auto valueGetter) {
            outFile << std::left << std::setw(widthDesc) << rowName;
            for (size_t s = 0; s < numSteps; ++s) {
                // Check, ob Step existiert
                if (s < node.MNodeGlobalSteps.size()) {
                    double val = valueGetter(s);
                    outFile << std::fixed << std::setprecision(6) 
                            << std::left << std::setw(widthVal) << val;
                } else {
                    outFile << std::setw(widthVal) << "NaN";
                }
            }
            outFile << "\n";
        };

        writeRow(nodePrefix + " Pos X", [&](size_t s){ return node.MNodeGlobalSteps[s].x; });
        writeRow(nodePrefix + " Pos Y", [&](size_t s){ return node.MNodeGlobalSteps[s].y; });
        writeRow(nodePrefix + " Pos Z", [&](size_t s){ return node.MNodeGlobalSteps[s].z; });

        // --- B) ETAS (Falls vorhanden) ---
        // Wir prüfen im ersten Step, wie viele Etas dieser Node hat
        size_t numEtas = 0;
        if (!node.MNodeEtaSteps.empty() && !node.MNodeEtaSteps[0].empty()) {
            numEtas = node.MNodeEtaSteps[0].size();
        }

        for (size_t e = 0; e < numEtas; ++e) {
            std::string etaName = nodePrefix + " Eta " + std::to_string(e);
            
            outFile << std::left << std::setw(widthDesc) << etaName;
            for (size_t s = 0; s < numSteps; ++s) {
                // Safety Check: Hat dieser Step auch Daten und ist der Eta-Index gültig?
                if (s < node.MNodeEtaSteps.size() && e < node.MNodeEtaSteps[s].size()) {
                    double val = node.MNodeEtaSteps[s][e];
                    outFile << std::scientific << std::setprecision(4) // Etas oft klein -> Scientific Notation
                            << std::left << std::setw(widthVal) << val;
                } else {
                    outFile << std::setw(widthVal) << "0.0"; // Fallback
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

inline void exportParameterLog(const std::vector<std::vector<double>>& values, 
                        const std::vector<std::vector<std::string>>& descriptions, 
                        const std::string& filename = "parameter_log.txt") 
{
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
    const int widthDesc = 30; // Breite der Beschreibungs-Spalte
    const int widthVal  = 15; // Breite der Wert-Spalten

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
                outFile << std::fixed << std::setprecision(6) 
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