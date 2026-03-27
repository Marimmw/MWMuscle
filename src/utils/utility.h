#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>      // Für std::setw, std::setprecision
#include <filesystem>   // Braucht C++17 (in CMake sicherstellen!)
#include <QDebug>
#include <QString>
#include <QProcess>
#include <QStringList>
#include <string>

#include <QDir>
#include <QDateTime>
#include <QFile>
#include <QJsonObject>
#include <QJsonArray>
#include <QJsonDocument>
#include <memory>

#include "simpleSimulation/SSMuscle.h"
#include "utils/MWMath.h"
#include "simpleSimulation/SSBody.h"
#include "simpleSimulation/SSMeshes.h"

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


inline std::vector<std::string> splitString(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    
    // std::getline liest aus dem Stream (ss) bis zum angegebenen Trennzeichen (delimiter)
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    
    return tokens;
}

//namespace fs = std::filesystem;
inline void printColorMsg(const std::string& message, int colorLevel) 
{
    std::string colorCode;
    const std::string resetCode = "\033[0m";

    switch(colorLevel) {
        case 0: 
            colorCode = "\033[31m";       // 0 = Rot (Fehler / Infeasible)
            break;
        case 1: 
            colorCode = "\033[38;5;208m"; // 1 = Orange (Max-Iter Warnung / Kritisch)
            break;
        case 2: 
            colorCode = "\033[33m";       // 2 = Gelb (Warnungen / Initialisierung)
            break;
        case 3: 
            colorCode = "\033[96m";       // 3 = Hellblau/Cyan (Infos / Status-Updates)
            break;
        case 4: 
            colorCode = "\033[32m";       // 4 = Grün (Erfolg / Solve_Succeeded)
            break;
        default: 
            colorCode = resetCode;        // Fallback: Standard Terminal-Farbe
            break;
    }

    // Alles zu einem einzigen String verschmelzen, damit qDebug keine Leerzeichen fälscht
    std::string fullMessage = colorCode + message + resetCode;
    
    // Sicher ausgeben
    qDebug().noquote() << QString::fromStdString(fullMessage);
}

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

inline void backupSourceCode(const std::string& callerFilePath) {
    namespace fs = std::filesystem;
    // 1. Pfad zur aktuellen Quelldatei (wird vom Compiler gesetzt)
    // fs::path sourceFile = __FILE__; 
    fs::path sourceFile = callerFilePath;

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




inline void uploadPoseStudyToFAUbox(const QString& txtFilePath) {
    qDebug() << "Starte FAUbox-Upload via Python...";

    QString pythonScript = "../examples/sendReport.py"; 

    QStringList arguments;
    arguments << pythonScript << txtFilePath;

    QProcess* process = new QProcess();
    
    QObject::connect(process, &QProcess::readyReadStandardOutput, [process]() {
        qDebug() << "FAUbox Upload Info:" << process->readAllStandardOutput().trimmed();
    });
    QObject::connect(process, &QProcess::readyReadStandardError, [process]() {
        qDebug() << "FAUbox Upload Error:" << process->readAllStandardError().trimmed();
    });

    process->start("python3", arguments);
    
    if (!process->waitForFinished(15000)) { // 15 Sekunden Timeout
        qDebug() << "FEHLER: Python-Skript hat nicht rechtzeitig geantwortet.";
    } else {
        qDebug() << "Upload-Skript beendet.";
    }
    
    process->deleteLater();
}



// export System
inline void addBodiesToJson(QJsonObject& rootJson, const std::vector<std::shared_ptr<SSTissue>>& bodies) {
    QJsonArray bodiesArray;

    for (const auto& body : bodies) {
        if (!body) continue; // Sicherheitshalber leere Pointer überspringen

        QJsonObject bodyObj;

        // 1. Name
        bodyObj["Name"] = QString::fromStdString(body->Name);

        // 2. Parent Name (Fallback zu "World")
        if (body->Parent != nullptr) {
            bodyObj["Parent"] = QString::fromStdString(body->Parent->Name);
        } else {
            bodyObj["Parent"] = "World";
        }

        // 3. MeshColor als Array [x, y, z]
        QJsonArray colorArray = {
            body->MeshColor.x, 
            body->MeshColor.y, 
            body->MeshColor.z
        };
        bodyObj["MeshColor"] = colorArray;

        // 4. Relative Position als Array [x, y, z]
        QJsonArray posArray = {
            body->Position2ParentRelInParentFrame.x, 
            body->Position2ParentRelInParentFrame.y, 
            body->Position2ParentRelInParentFrame.z
        };
        bodyObj["Position2ParentRelInParentFrame"] = posArray;

        // 5. Relative Orientierung als 3x3 Array (Array von Arrays)
        QJsonArray oriMatrix;
        for (int r = 0; r < 3; ++r) {
            QJsonArray rowArray = {
                body->Orientation2ParentRel.m[r][0],
                body->Orientation2ParentRel.m[r][1],
                body->Orientation2ParentRel.m[r][2]
            };
            oriMatrix.append(rowArray);
        }
        bodyObj["Orientation2ParentRel"] = oriMatrix;

        // 6. Children als Array von Strings
        QJsonArray childrenArray;
        for (const auto& child : body->Children) {
            if (child) {
                childrenArray.append(QString::fromStdString(child->Name));
            }
        }
        bodyObj["Children"] = childrenArray;

        QJsonArray meshNamesArray;
        for (const auto& mesh : body->Meshes) {
            if (mesh) {
                meshNamesArray.append(QString::fromStdString(mesh->Name));
            }
        }
        bodyObj["Meshes"] = meshNamesArray;

        // Das fertige Body-Objekt zum Haupt-Array hinzufügen
        bodiesArray.append(bodyObj);
    }

    // Das Array unter dem Schlüssel "Bodies" im Root-JSON speichern
    rootJson["Bodies"] = bodiesArray;
}

inline void addJointsToJson(QJsonObject& rootJson, const std::vector<std::shared_ptr<SSJoint>>& joints) {
    QJsonArray jointsArray;

    for (const auto& joint : joints) {
        if (!joint) continue;

        QJsonObject jointObj;

        // ==========================================
        // 1. BASIS-WERTE (von SSTissue geerbt)
        // ==========================================
        jointObj["Name"] = QString::fromStdString(joint->Name);

        if (joint->Parent != nullptr) {
            jointObj["Parent"] = QString::fromStdString(joint->Parent->Name);
        } else {
            jointObj["Parent"] = "World";
        }

        QJsonArray colorArray = {
            joint->MeshColor.x, 
            joint->MeshColor.y, 
            joint->MeshColor.z
        };
        jointObj["MeshColor"] = colorArray;

        QJsonArray posArray = {
            joint->Position2ParentRelInParentFrame.x, 
            joint->Position2ParentRelInParentFrame.y, 
            joint->Position2ParentRelInParentFrame.z
        };
        jointObj["Position2ParentRelInParentFrame"] = posArray;

        QJsonArray oriMatrix;
        for (int r = 0; r < 3; ++r) {
            QJsonArray rowArray = {
                joint->Orientation2ParentRel.m[r][0],
                joint->Orientation2ParentRel.m[r][1],
                joint->Orientation2ParentRel.m[r][2]
            };
            oriMatrix.append(rowArray);
        }
        jointObj["Orientation2ParentRel"] = oriMatrix;

        QJsonArray childrenArray;
        for (const auto& child : joint->Children) {
            if (child) {
                childrenArray.append(QString::fromStdString(child->Name));
            }
        }
        jointObj["Children"] = childrenArray;

        QJsonArray meshNamesArray;
        for (const auto& mesh : joint->Meshes) {
            if (mesh) {
                meshNamesArray.append(QString::fromStdString(mesh->Name));
            }
        }
        jointObj["Meshes"] = meshNamesArray;

        // ==========================================
        // 2. SSJOINT SPEZIFISCHE WERTE 
        // ==========================================
        
        // DoneAngleSteps -> AngleSteps
        QJsonArray angleStepsArray;
        for (double angle : joint->DoneAngleSteps) {
            angleStepsArray.append(angle);
        }
        jointObj["AngleSteps"] = angleStepsArray;

        // MaxAngles
        jointObj["MaxAngles"] = joint->MaxAngles;

        // RotationAxes als [x, y, z]
        QJsonArray rotAxisArray = {
            joint->RotationAxes.x,
            joint->RotationAxes.y,
            joint->RotationAxes.z
        };
        jointObj["RotationAxes"] = rotAxisArray;

        // TotalSteps
        jointObj["TotalSteps"] = joint->TotalSteps;

        // ==========================================

        jointsArray.append(jointObj);
    }

    // Unter dem Schlüssel "Joints" im JSON speichern
    rootJson["Joints"] = jointsArray; 
}

inline void addMeshesToJson(QJsonObject& rootJson, const std::vector<std::shared_ptr<SSMesh>>& meshes) {
    QJsonArray meshesArray;

    for (const auto& mesh : meshes) {
        if (!mesh) continue;

        QJsonObject meshObj;

        // ==========================================
        // 1. BASIS-WERTE (von SSTissue geerbt)
        // ==========================================
        meshObj["Name"] = QString::fromStdString(mesh->Name);

        if (mesh->Parent != nullptr) {
            meshObj["Parent"] = QString::fromStdString(mesh->Parent->Name);
        } else {
            meshObj["Parent"] = "World";
        }

        QJsonArray colorArray = {
            mesh->MeshColor.x, 
            mesh->MeshColor.y, 
            mesh->MeshColor.z
        };
        meshObj["MeshColor"] = colorArray;

        QJsonArray posArray = {
            mesh->Position2ParentRelInParentFrame.x, 
            mesh->Position2ParentRelInParentFrame.y, 
            mesh->Position2ParentRelInParentFrame.z
        };
        meshObj["Position2ParentRelInParentFrame"] = posArray;

        QJsonArray oriMatrix;
        for (int r = 0; r < 3; ++r) {
            QJsonArray rowArray = {
                mesh->Orientation2ParentRel.m[r][0],
                mesh->Orientation2ParentRel.m[r][1],
                mesh->Orientation2ParentRel.m[r][2]
            };
            oriMatrix.append(rowArray);
        }
        meshObj["Orientation2ParentRel"] = oriMatrix;

        QJsonArray childrenArray;
        for (const auto& child : mesh->Children) {
            if (child) {
                childrenArray.append(QString::fromStdString(child->Name));
            }
        }
        meshObj["Children"] = childrenArray;

        // ==========================================
        // 2. SSMESH SPEZIFISCHE WERTE
        // ==========================================
        
        // ViaPoint-Flag
        meshObj["bIsViaPoint"] = mesh->bIsViaPoint;

        // Typ und Dimensionen dynamisch ermitteln
        QString typeString = "UnknownMesh";
        QJsonArray dimensionsArray;

        // Prüfen auf Ellipsoid
        if (auto ell = dynamic_cast<SSEllipsoidMesh*>(mesh.get())) {
            typeString = "SSEllipsoidMesh";
            dimensionsArray = {ell->A, ell->B, ell->C};
        } 
        // Prüfen auf Zylinder
        else if (auto cyl = dynamic_cast<SSCylinderMesh*>(mesh.get())) {
            typeString = "SSCylinderMesh";
            dimensionsArray = {cyl->Radius, cyl->Height};
        } 
        // Prüfen auf Torus
        else if (auto tor = dynamic_cast<SSTorusMesh*>(mesh.get())) {
            typeString = "SSTorusMesh";
            // Da du im VTK-Code tor->A, tor->B, tor->C als Scale nutzt, 
            // packen wir sie der Vollständigkeit halber mit rein: [R, r, A, B, C]
            dimensionsArray = {tor->R, tor->r, tor->A, tor->B, tor->C};
        }

        meshObj["Type"] = typeString;
        meshObj["Dimensions"] = dimensionsArray;

        // ==========================================

        meshesArray.append(meshObj);
    }

    // Unter dem Schlüssel "Meshes" im JSON speichern
    rootJson["Meshes"] = meshesArray; 
}

inline void addMusclesToJson(QJsonObject& rootJson, const std::vector<SSMuscle*>& muscles) {
    QJsonArray musclesArray;

    for (const auto* mus : muscles) {
        if (!mus) continue;

        QJsonObject musObj;

        // Muskel Basis-Werte
        musObj["Name"] = QString::fromStdString(mus->Name);
        musObj["MNodesCount"] = mus->MNodesCount;
        musObj["MeshColor"] = QJsonArray{mus->MeshColor.x, mus->MeshColor.y, mus->MeshColor.z};
        
        musObj["OriginPointLocal"] = QJsonArray{mus->OriginPointLocal.x, mus->OriginPointLocal.y, mus->OriginPointLocal.z};
        musObj["parentMeshOrigin"] = (mus->parentMeshOrigin != nullptr) ? QString::fromStdString(mus->parentMeshOrigin->Name) : "World";
        musObj["InsertionPointLocal"] = QJsonArray{mus->InsertionPointLocal.x, mus->InsertionPointLocal.y, mus->InsertionPointLocal.z};
        musObj["parentMeshInsertion"] = (mus->parentMeshInsertion != nullptr) ? QString::fromStdString(mus->parentMeshInsertion->Name) : "World";

        QJsonArray meshPtrsArray;
        for (const auto* m : mus->meshPtrs) {
            if (m) meshPtrsArray.append(QString::fromStdString(m->Name));
        }
        musObj["meshPtrs"] = meshPtrsArray;

        // ==========================================
        // MUSCLE NODES & HISTORIE (Timesteps)
        // ==========================================
        QJsonArray nodesArray;
        for (const auto& node : mus->MNodes) {
            QJsonObject nodeObj;
            
            // Konstante Werte des Knotens
            nodeObj["bParentIsFixed"] = node.bParentIsFixed;

            // Zeitschritte durchlaufen
            QJsonArray timeStepsArray;
            size_t numSteps = node.MNodeGlobalSteps.size();

            for (size_t t = 0; t < numSteps; ++t) {
                QJsonObject stepObj; // Unser "Struct" für diesen Zeitschritt

                // 1. NEU: Parent-Mesh für EXAKT diesen Zeitschritt!
                // (Setzt voraus, dass du MNodeParentSteps in der MuscleNode Klasse hast)
                if (t < node.MNodeParentSteps.size() && node.MNodeParentSteps[t] != nullptr) {
                    stepObj["parentMesh"] = QString::fromStdString(node.MNodeParentSteps[t]->Name);
                } else {
                    stepObj["parentMesh"] = "None"; // Fallback
                }

                // 2. Globale Position
                if (t < node.MNodeGlobalSteps.size()) {
                    stepObj["GlobalPos"] = QJsonArray{
                        node.MNodeGlobalSteps[t].x, 
                        node.MNodeGlobalSteps[t].y, 
                        node.MNodeGlobalSteps[t].z
                    };
                }

                // 3. Initial Guess
                if (t < node.MNodeInitialGuessSteps.size()) {
                    stepObj["InitialGuess"] = QJsonArray{
                        node.MNodeInitialGuessSteps[t].x, 
                        node.MNodeInitialGuessSteps[t].y, 
                        node.MNodeInitialGuessSteps[t].z
                    };
                }

                // 4. Initial Guess Farbe
                if (t < node.MNodeInitialGuessColorSteps.size()) {
                    stepObj["InitialGuessColor"] = QJsonArray{
                        node.MNodeInitialGuessColorSteps[t].x, 
                        node.MNodeInitialGuessColorSteps[t].y, 
                        node.MNodeInitialGuessColorSteps[t].z
                    };
                }

                timeStepsArray.append(stepObj);
            }

            nodeObj["TimeSteps"] = timeStepsArray;
            
            // HIER wird der Node-Array befüllt (alles bleibt schön innerhalb des Muskels)
            nodesArray.append(nodeObj);
        }

        // HIER wird die Node-Liste dem Muskel zugewiesen
        musObj["MNodes"] = nodesArray;
        musclesArray.append(musObj);
    }

    rootJson["Muscles"] = musclesArray;
}

inline void addSettingsToJson(QJsonObject& rootJson, bool globalComputation, int objType, bool bSumPhiEta, bool bUseWarmstartEtas, double WarmstartEtaScaling, int maxIterations, double maxTol, bool bUseOwnGradient, int numTimeSteps) 
{   
    /* qDebug() << "Exportiere Simulationseinstellungen: " 
             << "\n  globalComputation:" << globalComputation
             << "\n  objType:" << objType
             << "\n  bSumPhiEta:" << bSumPhiEta
             << "\n  bUseWarmstartEtas:" << bUseWarmstartEtas
             << "\n  WarmstartEtaScaling:" << WarmstartEtaScaling
             << "\n  maxIterations:" << maxIterations
             << "\n  maxTol:" << maxTol
             << "\n  bUseOwnGradient:" << bUseOwnGradient
             << "\n  numTimeSteps:" << numTimeSteps; */
    QJsonObject cfgObj;
    
    // CasADi & Solver Settings
    cfgObj["globalComputation"] = globalComputation;
    cfgObj["objType"] = objType;
    cfgObj["bSumPhiEta"] = bSumPhiEta;
    cfgObj["bUseWarmstartEtas"] = bUseWarmstartEtas;
    cfgObj["WarmstartEtaScaling"] = WarmstartEtaScaling;
    cfgObj["maxIterations"] = maxIterations;
    cfgObj["maxTol"] = maxTol;
    cfgObj["bUseOwnGradient"] = bUseOwnGradient;
    
    // Globale Simulationseinstellungen
    cfgObj["numTimeSteps"] = numTimeSteps;

    // Alles unter dem Schlüssel "SimulationSettings" speichern
    rootJson["SimulationSettings"] = cfgObj;
}

inline void exportFullSceneToJson(const std::vector<std::shared_ptr<SSTissue>>& tissues, const std::vector<std::shared_ptr<SSMesh>>& meshes, const std::vector<SSMuscle*>& muscles,bool globalComputation, int objType, bool bSumPhiEta, bool bUseWarmstartEtas, double WarmstartEtaScaling, int maxIterations, double maxTol, bool bUseOwnGradient, int numTimeSteps) 
{
    std::vector<std::shared_ptr<SSTissue>> bodies;
    std::vector<std::shared_ptr<SSJoint>> joints;

    for (const auto& tissue : tissues) {
        // Versuche, das Tissue als SSJoint zu casten
        if (auto joint = std::dynamic_pointer_cast<SSJoint>(tissue)) {
            // Cast war erfolgreich -> Es ist ein Gelenk!
            joints.push_back(joint);
        } else {
            // Cast fehlgeschlagen -> Es ist ein normaler Body (Bone etc.)
            bodies.push_back(tissue);
        }
    }

    // 1. Ordnerpfad festlegen und erstellen, falls er nicht existiert
    QString outputFolder = "../examples/results/";
    QDir dir(outputFolder);
    if (!dir.exists()) {
        dir.mkpath("."); // Erstellt alle nötigen Unterordner
    }

    QString timestamp = QDateTime::currentDateTime().toString("yyyyMMdd_HHmmss");
    QString filename = "SceneExport_" + timestamp + ".json";
    QString fullPath = dir.filePath(filename);

    QJsonObject sceneJson;

    addSettingsToJson(sceneJson, globalComputation, objType, bSumPhiEta, 
                      bUseWarmstartEtas, WarmstartEtaScaling, maxIterations, 
                      maxTol, bUseOwnGradient, numTimeSteps);
    addBodiesToJson(sceneJson, bodies);
    addJointsToJson(sceneJson, joints);
    addMeshesToJson(sceneJson, meshes);
    addMusclesToJson(sceneJson, muscles);

    
    QJsonDocument doc(sceneJson);
    QFile jsonFile(fullPath);
    
    if (jsonFile.open(QIODevice::WriteOnly)) {
        // Indented sorgt für die schöne, gut lesbare Formatierung
        jsonFile.write(doc.toJson(QJsonDocument::Indented));
        jsonFile.close();
        qDebug() << "[JSON Export] Szene erfolgreich gespeichert unter:" << fullPath;
    } else {
        qDebug() << "[JSON Export] FEHLER: Konnte die Datei nicht erstellen:" << fullPath;
    }
}




inline bool loadSceneFromJson(const QString& filepath, 
                              std::vector<std::shared_ptr<SSTissue>>& outTissues, 
                              std::vector<std::shared_ptr<SSMesh>>& outMeshes, 
                              std::vector<SSMuscle*>& outMuscles,
                              std::shared_ptr<SSBody>& outRootSystem,
                              int& outNumTimeSteps) 
{
    QFile file(filepath);
    if (!file.open(QIODevice::ReadOnly)) {
        qDebug() << "[JSON Import] Fehler: Konnte Datei nicht oeffnen:" << filepath;
        return false;
    }

    QJsonDocument doc = QJsonDocument::fromJson(file.readAll());
    QJsonObject rootJson = doc.object();
    file.close();

    // LOADING SIM SETTINGS
    QJsonObject cfgObj = rootJson["SimulationSettings"].toObject();
    outNumTimeSteps = cfgObj["numTimeSteps"].toInt();
    // ...further

    // Hilfs-Maps zum schnellen Finden über den Namen
    std::map<std::string, std::shared_ptr<SSTissue>> tissueMap;
    std::map<std::string, std::shared_ptr<SSMesh>> meshMap;

    // Hilfs-Lambdas zum Lesen von Arrays
    auto parsePoint3D = [](const QJsonArray& arr) {
        return MWMath::Point3D(arr[0].toDouble(), arr[1].toDouble(), arr[2].toDouble());
    };
    auto parseRotMatrix = [](const QJsonArray& arr) {
        MWMath::RotMatrix3x3 mat;
        for (int i = 0; i < 3; ++i) {
            QJsonArray row = arr[i].toArray();
            mat.m[i][0] = row[0].toDouble();
            mat.m[i][1] = row[1].toDouble();
            mat.m[i][2] = row[2].toDouble();
        }
        return mat;
    };

    // =======================================================
    // PHASE 1: INSTANZIIERUNG (Bodies, Joints, Meshes)
    // =======================================================
    
    // 1.1 Bodies laden
    QJsonArray bodiesArray = rootJson["Bodies"].toArray();
    for (int i = 0; i < bodiesArray.size(); ++i) {
        QJsonObject bObj = bodiesArray[i].toObject();
        auto body = std::make_shared<SSBody>();
        body->Name = bObj["Name"].toString().toStdString();
        body->MeshColor = parsePoint3D(bObj["MeshColor"].toArray());
        body->Position2ParentRelInParentFrame = parsePoint3D(bObj["Position2ParentRelInParentFrame"].toArray());
        body->Orientation2ParentRel = parseRotMatrix(bObj["Orientation2ParentRel"].toArray());
        
        tissueMap[body->Name] = body;
        outTissues.push_back(body);
    }

    // 1.2 Joints laden
    QJsonArray jointsArray = rootJson["Joints"].toArray();
    for (int i = 0; i < jointsArray.size(); ++i) {
        QJsonObject jObj = jointsArray[i].toObject();
        auto joint = std::make_shared<SSJoint>();
        joint->Name = jObj["Name"].toString().toStdString();
        joint->MeshColor = parsePoint3D(jObj["MeshColor"].toArray());
        joint->Position2ParentRelInParentFrame = parsePoint3D(jObj["Position2ParentRelInParentFrame"].toArray());
        joint->Orientation2ParentRel = parseRotMatrix(jObj["Orientation2ParentRel"].toArray());
        
        // Joint-spezifisch
        joint->MaxAngles = jObj["MaxAngles"].toDouble();
        joint->TotalSteps = jObj["TotalSteps"].toInt();
        joint->RotationAxes = parsePoint3D(jObj["RotationAxes"].toArray());
        
        QJsonArray anglesArr = jObj["AngleSteps"].toArray();
        for (int a = 0; a < anglesArr.size(); ++a) joint->DoneAngleSteps.push_back(anglesArr[a].toDouble());
        
        outNumTimeSteps = std::max(outNumTimeSteps, joint->TotalSteps); // Ermittle numTimeSteps
        tissueMap[joint->Name] = joint;
        outTissues.push_back(joint);
    }

    // 1.3 Meshes laden
    QJsonArray meshesArray = rootJson["Meshes"].toArray();
    for (int i = 0; i < meshesArray.size(); ++i) {
        QJsonObject mObj = meshesArray[i].toObject();
        std::string type = mObj["Type"].toString().toStdString();
        std::string name = mObj["Name"].toString().toStdString();
        QJsonArray dims = mObj["Dimensions"].toArray();

        std::shared_ptr<SSMesh> mesh;
        
        // Konstruktoren OHNE den Namen aufrufen!
        if (type == "SSEllipsoidMesh") {
            mesh = std::make_shared<SSEllipsoidMesh>(dims[0].toDouble(), dims[1].toDouble(), dims[2].toDouble());
        } else if (type == "SSCylinderMesh") {
            mesh = std::make_shared<SSCylinderMesh>(dims[0].toDouble(), dims[1].toDouble());
        } else if (type == "SSTorusMesh") {
            mesh = std::make_shared<SSTorusMesh>(dims[0].toDouble(), dims[1].toDouble());
        } else {
            qDebug() << "[JSON Import] Unbekannter Mesh-Typ:" << QString::fromStdString(type) << "nutze Fallback.";
            // Fallback auf ein winziges Ellipsoid, da SSMesh abstrakt ist und nicht direkt erstellt werden darf
            mesh = std::make_shared<SSEllipsoidMesh>(0.01, 0.01, 0.01); 
        }

        // Namen NACHTRÄGLICH setzen
        mesh->Name = name;
        
        mesh->bIsViaPoint = mObj["bIsViaPoint"].toBool();
        mesh->MeshColor = parsePoint3D(mObj["MeshColor"].toArray());
        mesh->Position2ParentRelInParentFrame = parsePoint3D(mObj["Position2ParentRelInParentFrame"].toArray());
        mesh->Orientation2ParentRel = parseRotMatrix(mObj["Orientation2ParentRel"].toArray());

        meshMap[mesh->Name] = mesh;
        outMeshes.push_back(mesh);
    }

    // =======================================================
    // PHASE 2: BAUM VERKNÜPFEN (Parent, Children, Meshes)
    // =======================================================
    // Wir gehen durch die JSON nochmal durch und verknüpfen die Pointer
    
    auto linkTissue = [&](QJsonArray arr) {
        for (int i = 0; i < arr.size(); ++i) {
            QJsonObject obj = arr[i].toObject();
            auto tissue = tissueMap[obj["Name"].toString().toStdString()];
            
            // Parent setzen
            std::string parentName = obj["Parent"].toString().toStdString();
            if (parentName == "World") {
                outRootSystem = std::dynamic_pointer_cast<SSBody>(tissue);
            } else if (tissueMap.count(parentName)) {
                tissue->Parent = tissueMap[parentName];
            }

            // Children setzen
            QJsonArray childrenArr = obj["Children"].toArray();
            for (int c = 0; c < childrenArr.size(); ++c) {
                std::string childName = childrenArr[c].toString().toStdString();
                if (tissueMap.count(childName)) tissue->Children.push_back(tissueMap[childName]);
            }

            // Meshes setzen
            QJsonArray meshArr = obj["Meshes"].toArray();
            for (int m = 0; m < meshArr.size(); ++m) {
                std::string mName = meshArr[m].toString().toStdString();
                if (meshMap.count(mName)) {
                    tissue->Meshes.push_back(meshMap[mName]);
                    meshMap[mName]->Parent = tissue; // Mesh bekommt den Tissue als Parent
                }
            }
        }
    };
    linkTissue(bodiesArray);
    linkTissue(jointsArray);

    // =======================================================
    // PHASE 3: KINEMATISCHER PASS (Simulation abspielen)
    // =======================================================
    qDebug() << "[JSON Import] Rekonstruiere Kinematik für" << outNumTimeSteps << "Schritte...";
    if (outRootSystem) {
        for (int t = 0; t < outNumTimeSteps; ++t) {
            outRootSystem->update(t); // Update kaskadiert durch alle Children!
            
            // Mesh Historie füllen (analog zu deiner main.cpp)
            for (auto& m : outMeshes) {
                m->MeshPointsGlobal.push_back(m->PositionGlobal);
                m->allRMatrixGlobal.push_back(m->OrientationGlobal);
            }
            // Tissues Historie füllen
            for (auto& tObj : outTissues) {
                tObj->MeshPointsGlobal.push_back(tObj->PositionGlobal);
                tObj->allRMatrixGlobal.push_back(tObj->OrientationGlobal);
            }
        }
    } else {
        qDebug() << "[JSON Import] WARNUNG: Kein RootSystem gefunden!";
    }

    // =======================================================
    // PHASE 4: MUSKELN LADEN
    // =======================================================
    QJsonArray musclesArray = rootJson["Muscles"].toArray();
    for (int i = 0; i < musclesArray.size(); ++i) {
        QJsonObject musObj = musclesArray[i].toObject();
        
        // Parent Tissues ermitteln
        std::string oParentStr = musObj["parentMeshOrigin"].toString().toStdString();
        std::string iParentStr = musObj["parentMeshInsertion"].toString().toStdString();
        SSTissue* oriTissue = tissueMap.count(oParentStr) ? tissueMap[oParentStr].get() : outRootSystem.get();
        SSTissue* insTissue = tissueMap.count(iParentStr) ? tissueMap[iParentStr].get() : outRootSystem.get();

        MWMath::Point3D oriPos = parsePoint3D(musObj["OriginPointLocal"].toArray());
        MWMath::Point3D insPos = parsePoint3D(musObj["InsertionPointLocal"].toArray());
        int nodesCount = musObj["MNodesCount"].toInt();
        std::string name = musObj["Name"].toString().toStdString();

        SSMuscle* mus = new SSMuscle(name, nodesCount, oriTissue, oriPos, insTissue, insPos);
        mus->MeshColor = parsePoint3D(musObj["MeshColor"].toArray());

        // Verknüpfte Meshes
        QJsonArray mPtrsArr = musObj["meshPtrs"].toArray();
        for (int p = 0; p < mPtrsArr.size(); ++p) {
            std::string ptrName = mPtrsArr[p].toString().toStdString();
            if (meshMap.count(ptrName)) mus->meshPtrs.push_back(meshMap[ptrName].get());
        }

        // ===============================================
        // FEHLERBEHEBUNG: Vektor für Knoten reservieren!
        // ===============================================
        QJsonArray nodesArr = musObj["MNodes"].toArray();
        mus->MNodes.resize(nodesArr.size()); 

        for (int n = 0; n < nodesArr.size(); ++n) {
            QJsonObject nObj = nodesArr[n].toObject();
            mus->MNodes[n].bParentIsFixed = nObj["bParentIsFixed"].toBool();
            
            QJsonArray timeStepsArr = nObj["TimeSteps"].toArray();
            
            // Platz reservieren, um Neuallokationen zu sparen
            mus->MNodes[n].MNodeGlobalSteps.reserve(timeStepsArr.size());
            mus->MNodes[n].MNodeInitialGuessSteps.reserve(timeStepsArr.size());
            mus->MNodes[n].MNodeInitialGuessColorSteps.reserve(timeStepsArr.size());
            mus->MNodes[n].MNodeParentSteps.reserve(timeStepsArr.size());

            for (int t = 0; t < timeStepsArr.size(); ++t) {
                QJsonObject stepObj = timeStepsArr[t].toObject();
                
                // Globale Position laden
                mus->MNodes[n].MNodeGlobalSteps.push_back(parsePoint3D(stepObj["GlobalPos"].toArray()));
                
                // Initial Guess laden (falls vorhanden)
                if (stepObj.contains("InitialGuess")) {
                    mus->MNodes[n].MNodeInitialGuessSteps.push_back(parsePoint3D(stepObj["InitialGuess"].toArray()));
                }
                
                // Farben laden (falls vorhanden)
                if (stepObj.contains("InitialGuessColor")) {
                    mus->MNodes[n].MNodeInitialGuessColorSteps.push_back(parsePoint3D(stepObj["InitialGuessColor"].toArray()));
                }

                // Parent Mesh des Knotens laden (Wichtig für dynamische Reparametrisierung!)
                if (stepObj.contains("parentMesh")) {
                    std::string pName = stepObj["parentMesh"].toString().toStdString();
                    if (tissueMap.count(pName)) {
                        mus->MNodes[n].MNodeParentSteps.push_back(tissueMap[pName].get());
                    } else {
                        mus->MNodes[n].MNodeParentSteps.push_back(nullptr);
                    }
                }
            }
        }
        
        outMuscles.push_back(mus);
    }

    qDebug() << "[JSON Import] Erfolg! Geladen:" << outTissues.size() << "Tissues," 
             << outMeshes.size() << "Meshes," << outMuscles.size() << "Muscles.";
    
    return true;
}