#include "SimulationManager.h"
#include <iomanip>  // Für std::setw und std::left
#include <sstream>  // Für std::stringstream

SimulationManager::SimulationManager(const SimSettings& config) : m_cfg(config) {}

void SimulationManager::runParameterStudy(const std::vector<ParamDef>& paramDefs) {
    // 1. Gesamtzahl der Iterationen berechnen
    int totalRuns = 1;
    for (const auto& p : paramDefs) {
        totalRuns *= (p.steps > 0 ? p.steps : 1);
    }

    std::cout << "Starte Parameterstudie (" << totalRuns << " Iterationen, " 
              << paramDefs.size() << " Parameter)..." << std::endl;

    // Startzeit erfassen
    auto startTime = std::chrono::high_resolution_clock::now();

    // 2. Textdatei öffnen (Wir nutzen .txt für saubere Tabellen oder .csv)
    std::ofstream outFile("../examples/results/ParameterStudy_Summary.txt");
    
    int goodScore = m_cfg.numTimeSteps * 40;
    std::string scoreHeader = "Score(<" + std::to_string(goodScore) + ")";
    
    // HEADER SCHREIBEN (Sauber formatiert)
    outFile << std::left 
            << std::setw(15) << scoreHeader << "\t"
            << std::setw(15) << "Successes" << "\t";
    
    for (const auto& p : paramDefs) {
        outFile << std::setw(12) << p.name << "\t";
    }
    outFile << "Step_Details [Code(Iter)] (0=Succ, 1=Max, 2=Inf)\n";
    outFile << std::string(90 + paramDefs.size() * 15, '-') << "\n";

    // 3. Rekursion starten
    int currentRun = 0;
    std::vector<double> currentValues; 
    currentValues.reserve(paramDefs.size());

    runCombinationsRecursive(paramDefs, 0, currentValues, outFile, currentRun, totalRuns, startTime);

    outFile.close();
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto total_min = std::chrono::duration_cast<std::chrono::minutes>(endTime - startTime).count();
    std::cout << "Parameterstudie erfolgreich beendet in " << total_min << " Minuten!" << std::endl;
}

void SimulationManager::runCombinationsRecursive(
    const std::vector<ParamDef>& paramDefs,
    int currentDepth,
    std::vector<double>& currentValues,
    std::ofstream& outFile,
    int& currentRun,
    int totalRuns,
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime) 
{
    // ABBRUCHBEDINGUNG: Wir sind in der innersten Schleife angekommen!
    if (currentDepth == paramDefs.size()) {
        currentRun++; // Iteration hochzählen
        
        // 1. Simulation ausführen
        std::string resultLine = runSingleSimulation(currentValues);
        
        // 2. In Datei schreiben
        outFile << resultLine << "\n";
        outFile.flush(); 

        // 3. Konsolen-Output alle 5 Runs (oder beim letzten Run)
        if (currentRun % 1 == 0 || currentRun == totalRuns) {
            auto now = std::chrono::high_resolution_clock::now();
            auto elapsed_min = std::chrono::duration_cast<std::chrono::minutes>(now - startTime).count();
            std::cout << "[Parameter Study] Run " << currentRun << " / " << totalRuns 
                      << " | Verstrichene Zeit: " << elapsed_min << " min" << std::endl;
        }
        return; 
    }

    const ParamDef& currentDef = paramDefs[currentDepth];
    
    double stepSize = 0.0;
    if (currentDef.steps > 1) {
        stepSize = (currentDef.end - currentDef.start) / (currentDef.steps - 1);
    }

    for (int i = 0; i < currentDef.steps; ++i) {
        double val = currentDef.start + i * stepSize;
        currentValues.push_back(val); 
        
        runCombinationsRecursive(paramDefs, currentDepth + 1, currentValues, outFile, currentRun, totalRuns, startTime);
        
        currentValues.pop_back(); 
    }
}

std::string SimulationManager::runSingleSimulation(const std::vector<double>& params) {
    
    // 1. Container für diesen Durchlauf
    std::vector<std::shared_ptr<SSMesh>> meshes;
    std::vector<std::shared_ptr<SSTissue>> tissues;
    std::vector<SSMuscle*> musclePtrs;
    std::shared_ptr<SSBody> rootSystem;
    std::vector<CasadiSystem*> systems;

    // 2. Modell aufbauen (Hier gibst du deine Parameter p1-p4 an die Funktion weiter!)
    buildOHandModel(tissues, meshes, musclePtrs, rootSystem, m_cfg.numTimeSteps, m_cfg, 1.0, params);

    for (auto& m : meshes) {
        m->discretizeMesh(m_cfg.discretization);
    }

    // 3. Solver Setup
    for (auto* mus : musclePtrs) {
        std::vector<SSMuscle*> singleMuscleList = {mus};
        systems.push_back(new CasadiSystem(singleMuscleList, m_cfg.objFunc, m_cfg.solverMethod, m_cfg.casadiParametrization, m_cfg.bUseManualJacobian, m_cfg.bSumPhiEta, m_cfg.bUseWarmstartEtas, false));
        mus->initializeSimulationMuscle(m_cfg.numTimeSteps);
        mus->checkCollision();
    }

    // 4. SIMULATIONS-SCHLEIFE (Kein Viewer, keine unnötigen Prints)
    bool isSuccess = true;
    for(int t = 0; t < m_cfg.numTimeSteps; ++t) {
        if (rootSystem) {
            rootSystem->update(t); 
            for (auto& m : meshes) {
                m->MeshPointsGlobal.push_back(m->PositionGlobal);
                m->allRMatrixGlobal.push_back(m->OrientationGlobal);
                m->discretizeMesh(m_cfg.discretization);
            }
        }

        for(auto* mus : musclePtrs) {
            mus->OriginPointGlobal = mus->MNodes[0].predictNewGlobal();
            mus->InsertionPointGlobal = mus->MNodes.back().predictNewGlobal();
            mus->MNodes[0].PositionGlobal = mus->OriginPointGlobal;
            mus->MNodes.back().PositionGlobal = mus->InsertionPointGlobal;
            
            std::vector<MWMath::Point3D> guessPath;
            for(auto& node : mus->MNodes) guessPath.push_back(node.predictNewGlobal());
            mus->storeMNodesInitialGuess(t, guessPath);
        }

        // Casadi lösen
        for (auto* sys : systems) {
            sys->solveStepX(); 
            // Falls der Solver fehlschlägt, können wir isSuccess auf false setzen
            if (sys->SolverConvergenceMessages.back() != "Solve_Succeeded") {
                isSuccess = false;
            }
        }

        if (m_cfg.dynamicReparametrization || t < 1) {    
            for (auto* muscle : musclePtrs) muscle->updateMusclePointsParentsLocal();
        }

        for (auto* muscle : musclePtrs) {
            muscle->getViaPointNodeInfo();
            muscle->checkCollision();
            muscle->storeMNodesGlobalPositions();
            muscle->computeMuscleLength(true); // Länge berechnen
        }
    }

    // ==============================================================================
    // AUSWERTUNG DER CASADI LOGS
    // ==============================================================================
    int score = 0;
    int successCount = 0;
    std::string stepResults = "";

    // Wir werten das erste CasadiSystem aus (da dort die Logs deines Muskels liegen)
    if (!systems.empty()) {
        auto* sys = systems[0];
        
        for (size_t t = 0; t < sys->SolverConvergenceMessages.size(); ++t) {
            std::string msg = sys->SolverConvergenceMessages[t];
            int iters = sys->SolverConvergenceSteps[t];
            
            // Jeder Iterationsschritt kostet 1 Punkt
            score += iters;

            int statusCode = 3; // Fallback
            if (msg.find("Succeeded") != std::string::npos || msg.find("Success") != std::string::npos) {
                statusCode = 0;
                successCount++;
            } 
            else if (msg.find("Max") != std::string::npos || msg.find("Maximum") != std::string::npos) {
                statusCode = 1;
                // Optional: score += 5000; (Falls Max-Iter auch leicht bestraft werden soll)
            } 
            else if (msg.find("Infeasible") != std::string::npos) {
                statusCode = 2;
                score += 20000; // Harte Strafe für Infeasible
            }

            // String bauen: z.B. "0(45) "
            stepResults += std::to_string(statusCode) + "(" + std::to_string(iters) + ")  ";
        }
    }

    // MEMORY CLEANUP
    for (auto* sys : systems) delete sys;
    for (auto* mus : musclePtrs) delete mus;

    // ==============================================================================
    // TABELLEN-ZEILE GENERIEREN
    // ==============================================================================
    std::string succStr = std::to_string(successCount) + "/" + std::to_string(m_cfg.numTimeSteps);
    std::stringstream ss;
    ss << std::left 
       << std::setw(15) << score << "\t"
       << std::setw(15) << succStr << "\t";
       
    for (double p : params) {
        ss << std::setw(15) << p; // Parameter einfügen
    }
    ss << stepResults; // Am Ende die Historie anheften
                          
    return ss.str();
}