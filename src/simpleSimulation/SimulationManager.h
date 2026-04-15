#pragma once

#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <iostream>
#include <chrono>
#include <mutex>   // NEU: Für Thread-Sicherheit beim Dateischreiben
#include <atomic>  // NEU: Für den globalen Counter

#include "utils/ConfigLoader.h"
#include "simpleSimulation/SSMuscle.h"
#include "simpleSimulation/SSMeshes.h"
#include "simpleSimulation/SSBody.h"
#include "simpleSimulation/casadiSystem.h"
#include "utils/ModelBuilder.h"

struct ParamDef {
    double start;
    double end;
    int steps; // Wie viele Zwischenschritte?
    std::string name; // Optional: Name des Parameters für die Zusammenfassung oder spätere Analyse
};

struct PoseDef {
    std::string PoseName;
    std::vector<double> SimulationJointAngles; // e.g. 10 values for 10 joints (max joint angles, sim is beginning from 0°)
    std::vector<std::string> JointNames; // for analyisis
};

class SimulationManager {
public:
    SimulationManager(const SimSettings& config);

    // Nimmt jetzt dynamisch eine Liste von Parameter-Konfigurationen auf
    void runParameterStudy(const std::vector<ParamDef>& paramDefs);

    void runPoseStudy(const std::vector<PoseDef>& poses);
    void runViaPointParamStudy(const std::vector<PoseDef>& poses);

    // Nimmt den fertigen Vektor mit den aktuellen Werten für diesen EINEN Run auf
    std::vector<std::string> runSingleSimulation(const std::vector<double>& params);
    std::vector<std::string> runSingleSimulationViaPointParam(const std::vector<double>& params);

    std::vector<PoseDef> createPoseDefs();
    std::vector<PoseDef> createParameterViaPointStudy();

    bool bExportLogSumValues = false; 

private:
    SimSettings m_cfg;

    // NEU: Rekursive Hilfsfunktion, die unsere dynamischen FOR-Schleifen baut
    void runCombinationsRecursive(
        const std::vector<ParamDef>& paramDefs,
        int currentDepth,
        std::vector<double>& currentValues,
        std::ofstream& outFile,
        int& currentRun,
        int totalRuns,
        std::chrono::time_point<std::chrono::high_resolution_clock> startTime); // Starttime

    void generateCombinations(
        const std::vector<ParamDef>& paramDefs,
        int currentDepth,
        std::vector<double>& currentValues,
        std::vector<std::vector<double>>& allJobs
    );
};