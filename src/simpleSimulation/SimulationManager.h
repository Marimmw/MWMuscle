#pragma once

#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <iostream>
#include <chrono>

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

class SimulationManager {
public:
    SimulationManager(const SimSettings& config);

    // Nimmt jetzt dynamisch eine Liste von Parameter-Konfigurationen auf
    void runParameterStudy(const std::vector<ParamDef>& paramDefs);

    // Nimmt den fertigen Vektor mit den aktuellen Werten für diesen EINEN Run auf
    std::string runSingleSimulation(const std::vector<double>& params);

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
};