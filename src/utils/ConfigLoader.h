#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <filesystem> // C++17 feature
#include <ctime>
#include <iomanip>
#include <QDebug>
#include <QVector>
#include <QString>

namespace fs = std::filesystem;

struct SimSettings {
    std::vector<double> MAXJOINTANGLES;
    bool bFingerIsEllipsoid = false;
    bool bCreateJointMeshes = true;
    std::vector<int> muscleNumPoints; // Liste für Muskeln
    std::string currentScene = "FT_translated";
    int numTimeSteps = 3;
    int objFunc = 0;
    int discretization = 7;
    std::string casadiParametrization = "local";
    std::string solverMethod = "VPP"; // "VPP", "VPPenalty", "Exp"
    bool bUseManualJacobian = false;
    bool dynamicReparametrization = true;
    bool bShowDiscretization = false;
    bool bShowSolverInVisualization = false;
    bool bSumPhiEta = false;
    bool bUseWarmstartEtas = true;
};

class ConfigManager {
public:
    static SimSettings loadConfig(const std::string& filepath) {
        SimSettings settings;
        std::ifstream file(filepath);
        
        if (!file.is_open()) {
            std::cerr << "ACHTUNG: Konnte config.txt nicht finden! Nutze Defaults." << std::endl;
            // Default Fallback
            settings.MAXJOINTANGLES = {90.0, 100.0, 90.0};
            settings.muscleNumPoints = {25};
            return settings;
        }

        std::string line;
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string key;
            ss >> key;

            if (key.empty() || key[0] == '#') continue; // Skip comments

            if (key == "MAXJOINTANGLES") {
                double val;
                while (ss >> val) settings.MAXJOINTANGLES.push_back(val);
            }
            else if (key == "bFingerIsEllipsoid") ss >> settings.bFingerIsEllipsoid;
            else if (key == "bCreateJointMeshes") ss >> settings.bCreateJointMeshes;
            else if (key == "muscleNumPoints") {
                int val;
                while (ss >> val) settings.muscleNumPoints.push_back(val);
            }
            else if (key == "currentScene") ss >> settings.currentScene;
            else if (key == "numTimeSteps") ss >> settings.numTimeSteps;
            else if (key == "objFunc") ss >> settings.objFunc;
            else if (key == "discretization") ss >> settings.discretization;
            else if (key == "casadiParametrization") ss >> settings.casadiParametrization;
            else if (key == "solverMethod") ss >> settings.solverMethod;
            else if (key == "bUseManualJacobian") ss >> settings.bUseManualJacobian;
            else if (key == "dynamicReparametrization") ss >> settings.dynamicReparametrization;
            else if (key == "bShowDiscretization") ss >> settings.bShowDiscretization;
            else if (key == "bShowSolverInVisualization") ss >> settings.bShowSolverInVisualization;
            else if (key == "bSumPhiEta") ss >> settings.bSumPhiEta;
            else if (key == "bUseWarmstartEtas") ss >> settings.bUseWarmstartEtas;
        }

        qDebug() << "-------------- CONFIG --------------";
        qDebug() << "   MAXJOINTANGLES:" << QVector<double>(settings.MAXJOINTANGLES.begin(), settings.MAXJOINTANGLES.end());
        qDebug() << "   bFingerIsEllipsoid:" << settings.bFingerIsEllipsoid;
        qDebug() << "   bCreateJointMeshes:" << settings.bCreateJointMeshes;
        qDebug() << "   muscleNumPoints:" << QVector<int>(settings.muscleNumPoints.begin(), settings.muscleNumPoints.end());
        qDebug() << "   currentScene:" << QString::fromStdString(settings.currentScene);
        qDebug() << "   numTimeSteps:" << settings.numTimeSteps;
        qDebug() << "   objFunc:" << settings.objFunc;
        qDebug() << "   discretization:" << settings.discretization;
        qDebug() << "   casadiParametrization:" << QString::fromStdString(settings.casadiParametrization);
        qDebug() << "   solverMethod:" << QString::fromStdString(settings.solverMethod);
        qDebug() << "   bUseManualJacobian:" << settings.bUseManualJacobian;
        qDebug() << "   dynamicReparametrization:" << settings.dynamicReparametrization;
        qDebug() << "   bShowDiscretization:" << settings.bShowDiscretization;
        qDebug() << "   bShowSolverInVisualization:" << settings.bShowSolverInVisualization;
        qDebug() << "   bSumPhiEta:" << settings.bSumPhiEta;
        qDebug() << "   bUseWarmstartEtas:" << settings.bUseWarmstartEtas;
        qDebug() << "------------------------------------";
        return settings;
    }

    static void saveCopy(const std::string& srcFile) {
        // Pfad zusammenbauen: ~/projects/MWmuscleWrapper/examples
        const char* homeDir = getenv("HOME");
        if (!homeDir) return;
        
        fs::path targetDir = fs::path(homeDir) / "projects/MWmuscleWrapper/examples/results";
        
        // Ordner erstellen falls nicht existent
        if (!fs::exists(targetDir)) {
            fs::create_directories(targetDir);
        }

        // Dateiname mit Zeitstempel: SystemParameters_YYYY-MM-DD_HH-MM-SS.txt
        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        std::ostringstream oss;
        oss << "SystemParameters_" ; //<< std::put_time(&tm, "%Y-%m-%d_%H-%M-%S") << ".txt";
        
        fs::path targetPath = targetDir / oss.str();

        // Kopieren
        try {
            fs::copy_file(srcFile, targetPath, fs::copy_options::overwrite_existing);
            std::cout << "Konfiguration gesichert unter: " << targetPath << std::endl;
        } catch (fs::filesystem_error& e) {
            std::cerr << "Fehler beim Sichern der Config: " << e.what() << std::endl;
        }
    }
};