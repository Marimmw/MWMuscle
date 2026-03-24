#include "SimulationManager.h"
#include <iomanip>  // Für std::setw und std::left
#include <sstream>  // Für std::stringstream
#include <omp.h>

// ANSI Farbcodes
const std::string COLOR_RESET  = "\033[0m";
const std::string COLOR_GREEN  = "\033[32m";
const std::string COLOR_RED    = "\033[31m";
const std::string COLOR_YELLOW = "\033[33m";
const std::string COLOR_CYAN   = "\033[36m";

std::mutex casadiSetupMutex;

SimulationManager::SimulationManager(const SimSettings& config) : m_cfg(config) {}

void SimulationManager::runParameterStudy(const std::vector<ParamDef>& paramDefs) {
    /* // 1. Gesamtzahl der Iterationen berechnen
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
    std::cout << "Parameterstudie erfolgreich beendet in " << total_min << " Minuten!" << std::endl; */

    auto startTime = std::chrono::high_resolution_clock::now();

    // 1. JOBS GENERIEREN (Sekundensache)
    std::vector<std::vector<double>> allJobs;
    std::vector<double> currentValues;
    currentValues.reserve(paramDefs.size());
    generateCombinations(paramDefs, 0, currentValues, allJobs);
    
    int totalRuns = allJobs.size();
    std::cout << "Starte Parameterstudie (" << totalRuns << " Iterationen, auf mehreren Kernen)..." << std::endl;

    // 2. DATEI VORBEREITEN
    std::ofstream outFile("../examples/results/ParameterStudy_Summary.txt");
    int goodScore = m_cfg.numTimeSteps * 40;
    std::string scoreHeader = "Score(<" + std::to_string(goodScore) + ")";

    outFile << std::left << std::setw(35) << "Casadi System" << "\t" << std::setw(15) << scoreHeader << "\t" << std::setw(15) << "Successes" << "\t";
    for (const auto& p : paramDefs) outFile << std::setw(12) << p.name << "\t";
    outFile << "Step_Details [Code(Iter)] (0=Succ, 1=Max, 2=Inf)\n";
    outFile << std::string(90 + paramDefs.size() * 15, '-') << "\n";

    // 3. THREAD-SICHERE VARIABLEN
    std::mutex fileMutex;          // Sichert das Schreiben in die Datei
    std::atomic<int> currentRun{0}; // Zählt sicher hoch, egal welcher Thread gerade fertig wird

    // ==============================================================================
    // 4. PARALLELE SCHLEIFE (Hier passiert die Magie!)
    // schedule(dynamic) teilt die Jobs intelligent zu: Wenn ein Kern früher 
    // fertig ist (z.B. weil Casadi schnell infeasible war), bekommt er direkt den nächsten Job.
    // ==============================================================================
    int maxThreads = omp_get_max_threads();
    int useThreads = 1; // (maxThreads > 2) ? maxThreads - 2 : 1;
    std::cout << "Nutze " << useThreads << " von " << maxThreads << " verfuegbaren Threads." << std::endl;
    #pragma omp parallel for schedule(dynamic) num_threads(useThreads)

    for (int i = 0; i < totalRuns; ++i) {
        
        // Führe die Simulation lokal im Thread aus (NEU: fängt Vektor auf)
        std::vector<std::string> resultLines = runSingleSimulation(allJobs[i]);
        
        // --- KRITISCHER BEREICH ---
        {
            std::lock_guard<std::mutex> lock(fileMutex);
            
            // NEU: Für jeden Muskel eine eigene Zeile schreiben
            for (const auto& line : resultLines) {
                outFile << line << "\n";
            }
            outFile.flush();
            
            currentRun++;
            
            // Output alle 10 Runs (oder beim letzten)
            if (currentRun % 1 == 0 || currentRun == totalRuns) {
                auto now = std::chrono::high_resolution_clock::now();
                auto elapsed_min = std::chrono::duration_cast<std::chrono::minutes>(now - startTime).count();
                
                std::string printColor = COLOR_RESET;
                bool isSuccess = false;
                
                // Prüfe den Erfolg am ersten Muskel
                if (!resultLines.empty()) {
                    isSuccess = resultLines[0].find(std::to_string(m_cfg.numTimeSteps) + "/" + std::to_string(m_cfg.numTimeSteps)) != std::string::npos;
                }

                if (isSuccess) printColor = COLOR_GREEN;
                else printColor = COLOR_RED;

                std::cout << COLOR_CYAN << "[Study Progress] " << COLOR_RESET 
                          << "Run " << currentRun << " / " << totalRuns 
                          << " | Letzter Status: " << printColor << (isSuccess ? "SUCCESS" : "FAILED") << COLOR_RESET
                          << " | Zeit: " << COLOR_YELLOW << elapsed_min << " min" << COLOR_RESET 
                          << std::endl;
            }
        }
    }

    outFile.close();
    auto endTime = std::chrono::high_resolution_clock::now();
    auto total_min = std::chrono::duration_cast<std::chrono::minutes>(endTime - startTime).count();
    std::cout << "\n" << COLOR_GREEN << "Parameterstudie erfolgreich beendet in " << total_min << " Minuten!" << COLOR_RESET << std::endl;

}

void SimulationManager::runPoseStudy(const std::vector<PoseDef>& poses) {
    auto startTime = std::chrono::high_resolution_clock::now();
    int totalRuns = poses.size();
    
    std::cout << "Starte Posen-Studie (" << totalRuns << " Griffarten, auf mehreren Kernen)..." << std::endl;

    // 1. DATEI VORBEREITEN
    std::ofstream outFile("../examples/results/PoseStudy_Summary.txt");
    int goodScore = m_cfg.numTimeSteps * 40;
    std::string scoreHeader = "Score(<" + std::to_string(goodScore) + ")";

    // Header schreiben
    outFile << std::left << std::setw(35) << "Pose Name" << "\t" 
            << std::setw(35) << "Casadi System" << "\t"
            << std::setw(15) << scoreHeader << "\t" 
            << std::setw(15) << "Successes" << "\t";
            
    // Wir nehmen die Gelenk-Namen aus der ersten Pose für den Tabellenkopf
    if (!poses.empty() && !poses[0].JointNames.empty()) {
        for (const auto& jName : poses[0].JointNames) {
            outFile << std::setw(10) << jName << "\t";
        }
    }
    
    outFile << "Step_Details [Code(Iter)] (0=Succ, 1=Max, 2=Inf)\n";
    outFile << std::string(120, '-') << "\n";

    // 2. THREAD-SICHERE VARIABLEN
    std::mutex fileMutex;
    std::atomic<int> currentRun{0}; 

    int maxThreads = omp_get_max_threads();
    int useThreads = 1; // (maxThreads > 2) ? maxThreads - 2 : 1;
    std::cout << "Nutze " << useThreads << " von " << maxThreads << " verfuegbaren Threads." << std::endl;

    // 3. PARALLELE SCHLEIFE
    #pragma omp parallel for schedule(dynamic) num_threads(useThreads)
    for (int i = 0; i < totalRuns; ++i) {
        
        // Die Simulation wird mit den maximalen Gelenkwinkeln dieser Pose aufgerufen!
        // Der ModelBuilder in runSingleSimulation() schnappt sich params und steuert damit die Gelenke.
        /* std::string resultLine = runSingleSimulation(poses[i].SimulationJointAngles);
        
        // --- KRITISCHER BEREICH ---
        {
            std::lock_guard<std::mutex> lock(fileMutex);
            
            // Wir fügen den Posen-Namen GANZ VORNE in den Result-String ein, damit er in der Tabelle steht
            outFile << std::left << std::setw(20) << poses[i].PoseName << "\t" << resultLine << "\n";
            outFile.flush();
            
            currentRun++;
            
            auto now = std::chrono::high_resolution_clock::now();
            auto elapsed_min = std::chrono::duration_cast<std::chrono::minutes>(now - startTime).count();
            
            std::string printColor = COLOR_RESET;
            bool isSuccess = resultLine.find(std::to_string(m_cfg.numTimeSteps) + "/" + std::to_string(m_cfg.numTimeSteps)) != std::string::npos;
            if (isSuccess) printColor = COLOR_GREEN;
            else if (!isSuccess) printColor = COLOR_RED;

            std::cout << COLOR_CYAN << "[Pose Study] " << COLOR_RESET 
                      << "Finished " << currentRun << " / " << totalRuns 
                      << " | Pose: " << std::setw(15) << std::left << poses[i].PoseName
                      << " | Status: " << printColor << (isSuccess ? "SUCCESS" : "FAILED/MAX_ITER") << COLOR_RESET
                      << " | Zeit: " << COLOR_YELLOW << elapsed_min << " min" << COLOR_RESET 
                      << std::endl;
        } */
        std::vector<std::string> resultLines = runSingleSimulation(poses[i].SimulationJointAngles);
        
        // --- KRITISCHER BEREICH ---
        {
            std::lock_guard<std::mutex> lock(fileMutex);
            
            // Wir schreiben für jeden Muskel eine eigene Zeile (mit demselben Posen-Namen davor)
            for (const auto& line : resultLines) {
                outFile << std::left << std::setw(25) << poses[i].PoseName << "\t" << line << "\n";
            }
            outFile.flush();
            
            currentRun++;
            auto now = std::chrono::high_resolution_clock::now();
            auto elapsed_min = std::chrono::duration_cast<std::chrono::minutes>(now - startTime).count();
            
            std::string printColor = COLOR_RESET;
            bool isSuccess = false; // NEU: Standardmäßig auf false

            // NEU: Wir prüfen den ersten String im Vektor (falls vorhanden)
            if (!resultLines.empty()) {
                isSuccess = resultLines[0].find(std::to_string(m_cfg.numTimeSteps) + "/" + std::to_string(m_cfg.numTimeSteps)) != std::string::npos;
            }

            if (isSuccess) printColor = COLOR_GREEN;
            else printColor = COLOR_RED; // (Das else if (!isSuccess) ist redundant)

            std::cout << COLOR_CYAN << "[Pose Study] " << COLOR_RESET 
                      << "Finished " << currentRun << " / " << totalRuns 
                      << " | Pose: " << std::setw(15) << std::left << poses[i].PoseName
                      << " | Status: " << printColor << (isSuccess ? "SUCCESS" : "FAILED/MAX_ITER") << COLOR_RESET
                      << " | Zeit: " << COLOR_YELLOW << elapsed_min << " min" << COLOR_RESET 
                      << std::endl;
        }

    }

    outFile.close();
    auto endTime = std::chrono::high_resolution_clock::now();
    auto total_min = std::chrono::duration_cast<std::chrono::minutes>(endTime - startTime).count();
    std::cout << "\n" << COLOR_GREEN << "Posen-Studie erfolgreich beendet in " << total_min << " Minuten!" << COLOR_RESET << std::endl;
}

void SimulationManager::runCombinationsRecursive(const std::vector<ParamDef>& paramDefs, int currentDepth, std::vector<double>& currentValues, std::ofstream& outFile, int& currentRun, int totalRuns, std::chrono::time_point<std::chrono::high_resolution_clock> startTime) 
{
    // ABBRUCHBEDINGUNG: Wir sind in der innersten Schleife angekommen!
    if (currentDepth == paramDefs.size()) {
        currentRun++; // Iteration hochzählen
        
        // 1. Simulation ausführen
        std::vector<std::string> resultLines = runSingleSimulation(currentValues);
        for (const auto& line : resultLines) {
            outFile << line << "\n";
        }
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

std::vector<std::string> SimulationManager::runSingleSimulation(const std::vector<double>& params) {
    
    // 1. Container für diesen Durchlauf
    std::vector<std::shared_ptr<SSMesh>> meshes;
    std::vector<std::shared_ptr<SSTissue>> tissues;
    std::vector<SSMuscle*> musclePtrs;
    std::shared_ptr<SSBody> rootSystem;
    std::vector<CasadiSystem*> systems;

    // 2. Modell aufbauen (Hier gibst du deine Parameter p1-p4 an die Funktion weiter!)
    std::string buildResult = buildOHandModelOldExpandedViaX05(tissues, meshes, musclePtrs, rootSystem, m_cfg.numTimeSteps, m_cfg, 1.0, params);
    qDebug() << "        | ------------------------------------------------------- |";
    qDebug() << "        | Using System: " << QString::fromStdString(buildResult);
    qDebug() << "        | ------------------------------------------------------- |";

    for (auto& m : meshes) {
        m->discretizeMesh(m_cfg.discretization);
    }

    // 3. Solver Setup
    // --- KRITISCHER BEREICH: CasADi Setup ---
    {
        std::lock_guard<std::mutex> lock(casadiSetupMutex);
        for (auto* mus : musclePtrs) {
            std::vector<SSMuscle*> singleMuscleList = {mus};
            systems.push_back(new CasadiSystem(singleMuscleList, m_cfg.objFunc, m_cfg.solverMethod, m_cfg.casadiParametrization, m_cfg.bUseManualJacobian, m_cfg.bSumPhiEta, m_cfg.bUseWarmstartEtas, false, false));
            mus->initializeSimulationMuscle(m_cfg.numTimeSteps);
        }
    }
    // ------------------------------------------

    for (auto* mus : musclePtrs) {
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
            muscle->checkTorusSnapThrough();
        }
    }
    
    // ==============================================================================
    // AUSWERTUNG DER CASADI LOGS FÜR ALLE SYSTEME
    // ==============================================================================
    std::vector<std::string> resultLines;

    // Wir iterieren über ALLE Systeme (Reihenfolge ist identisch mit musclePtrs)
    for (size_t s = 0; s < systems.size(); ++s) {
        auto* sys = systems[s];
        auto* mus = musclePtrs[s]; 

        int score = 0;
        int successCount = 0;
        std::string stepResults = "";

        for (size_t t = 0; t < sys->SolverConvergenceMessages.size(); ++t) {
            std::string msg = sys->SolverConvergenceMessages[t];
            int iters = sys->SolverConvergenceSteps[t];
            
            score += iters;

            int statusCode = 3; // Fallback
            if (msg.find("Succeeded") != std::string::npos || msg.find("Success") != std::string::npos) {
                statusCode = 0;
                successCount++;
            } 
            else if (msg.find("Max") != std::string::npos || msg.find("Maximum") != std::string::npos) {
                statusCode = 1;
            } 
            else if (msg.find("Infeasible") != std::string::npos || msg.find("Invalid_Number") != std::string::npos) {
                statusCode = 2;
                score += 20000; 
            }

            stepResults += std::to_string(statusCode) + "(" + std::to_string(iters) + ")  ";
        }

        // ==============================================================================
        // TABELLEN-ZEILE GENERIEREN FÜR DIESEN MUSKEL
        // ==============================================================================
        std::string systemName = "CasSys_" + mus->Name;
        std::string succStr = std::to_string(successCount) + "/" + std::to_string(m_cfg.numTimeSteps);
        
        std::stringstream ss;
        ss << std::left 
           << std::setw(35) << systemName << "\t"
           << std::setw(15) << score
           << std::setw(15) << succStr;
           
        for (double p : params) {
            ss << std::setw(15) << p; 
        }
        ss << stepResults; 
        
        resultLines.push_back(ss.str());
    }

    // MEMORY CLEANUP
    for (auto* sys : systems) delete sys;
    for (auto* mus : musclePtrs) delete mus;

    return resultLines;
}

std::vector<PoseDef> SimulationManager::createPoseDefs()
{
    // 1. Erweitere die Header-Namen um die neue Spalte für die Knoten
    std::vector<std::string> jointNames = {"Wrist_F", "Wrist_A", "MCP_F", "MCP_A", "PIP", "DIP", "NumNodes", "VP_Size"};
    
    // 2. Deine Parameter
    std::vector<double> numNodes = {40,50,60,65,70,75,80,85,90,100,125};
    std::vector<double> vpSize = {0.05};

    std::vector<std::vector<double>> extraParameters = {numNodes, vpSize};
    // 3. Temporäre Struktur für die Basis-Posen (ohne numNodes)
    struct BasePose {
        std::string name;
        std::vector<double> angles;
    };
    
    std::vector<BasePose> basePoses = {
        {"Krampf Pose",    { 0.0,  80.0, 90.0,  0.0, 100.0, 80.0}},
        {"Full Fist",      { 0.0,   0.0, 90.0,  0.0, 100.0, 80.0}},
        {"Dach-Position",  { 0.0,   0.0, 90.0,  0.0,   0.0,  0.0}},
        {"Krallengriff",   { 0.0,   0.0,  0.0,  0.0, 100.0, 80.0}},
        {"Schraeger Griff",{ 0.0,   0.0, 45.0, 20.0,  45.0, 45.0}},
        {"Power Grip",     {20.0,   0.0, 90.0,  0.0, 100.0, 80.0}},
        {"Dart-Wurf",      { 0.0, -60.0, 90.0,  0.0, 100.0, 80.0}}
    };
    
    std::vector<PoseDef> variantList;
    
    // 4. Verschachtelte Schleife: Kombiniere jede Pose mit jeder Knotenanzahl
    bool bReverseOrder = true; // Wenn true, werden die Knoten von groß nach klein kombiniert
    if (bReverseOrder){
        for (double vp : vpSize) {
            for (double nodes : numNodes) {
                for (const auto& bp : basePoses) {

                    PoseDef p;

                    std::ostringstream vpStream;
                    vpStream << std::fixed << std::setprecision(4) << vp;

                    p.PoseName = bp.name + 
                        " (N=" + std::to_string(static_cast<int>(nodes)) +
                        ", vp=" + vpStream.str() + ")";

                    p.JointNames = jointNames;

                    // Reihenfolge bleibt: nodes, vp vorne
                    p.SimulationJointAngles = {nodes, vp};

                    p.SimulationJointAngles.insert(
                        p.SimulationJointAngles.end(),
                        bp.angles.begin(),
                        bp.angles.end()
                    );

                    variantList.push_back(p);
                }
            }
        }
    }
    else {
        for (const auto& bp : basePoses) {
            for (double nodes : numNodes) {
                for (double vp : vpSize) {

                    PoseDef p;

                    p.PoseName = bp.name + 
                        " (N=" + std::to_string(static_cast<int>(nodes)) +
                        ", vp=" + std::to_string(vp) + ")";

                    p.JointNames = jointNames;

                    p.SimulationJointAngles = bp.angles;

                    p.SimulationJointAngles.push_back(nodes); // NumNodes
                    p.SimulationJointAngles.push_back(vp);    // vpSize

                    variantList.push_back(p);
                }
            }
        }
    }
    
    return variantList;
}

// Rekursive Funktion füllt nur noch die Job-Liste auf
void SimulationManager::generateCombinations(const std::vector<ParamDef>& paramDefs,int currentDepth,std::vector<double>& currentValues,std::vector<std::vector<double>>& allJobs) 
{
    if (currentDepth == paramDefs.size()) {
        allJobs.push_back(currentValues);
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
        generateCombinations(paramDefs, currentDepth + 1, currentValues, allJobs);
        currentValues.pop_back(); 
    }
}


