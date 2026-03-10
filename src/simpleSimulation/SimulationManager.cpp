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

    outFile << std::left << std::setw(15) << scoreHeader << "\t" << std::setw(15) << "Successes" << "\t";
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
        
        // Führe die Simulation lokal im Thread aus
        std::string resultLine = runSingleSimulation(allJobs[i]);
        
        // --- KRITISCHER BEREICH (Immer nur 1 Thread darf hier gleichzeitig rein) ---
        {
            std::lock_guard<std::mutex> lock(fileMutex);
            outFile << resultLine << "\n";
            outFile.flush();
            
            currentRun++;
            
            // Output alle 10 Runs (oder beim letzten)
            if (currentRun % 1 == 0 || currentRun == totalRuns) {
                auto now = std::chrono::high_resolution_clock::now();
                auto elapsed_min = std::chrono::duration_cast<std::chrono::minutes>(now - startTime).count();
                
                std::string printColor = COLOR_RESET;
                bool isSuccess = resultLine.find(std::to_string(m_cfg.numTimeSteps) + "/" + std::to_string(m_cfg.numTimeSteps)) != std::string::npos;
                if (isSuccess) printColor = COLOR_GREEN;
                else if (!isSuccess) printColor = COLOR_RED;

                std::cout << COLOR_CYAN << "[Study Progress] " << COLOR_RESET 
                          << "Run " << currentRun << " / " << totalRuns 
                          << " | Letzter Status: " << printColor << (isSuccess ? "SUCCESS" : "FAILED") << COLOR_RESET
                          << " | Zeit: " << COLOR_YELLOW << elapsed_min << " min" << COLOR_RESET 
                          << std::endl;
            }
        }
        // --- KRITISCHER BEREICH ENDE ---
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
    outFile << std::left << std::setw(20) << "Pose Name" << "\t" 
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
        std::string resultLine = runSingleSimulation(poses[i].SimulationJointAngles);
        
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
    std::string buildResult = buildOHandModelOldExpandedViaX(tissues, meshes, musclePtrs, rootSystem, m_cfg.numTimeSteps, m_cfg, 1.0, params);
    qDebug() << "        | Using System: " << QString::fromStdString(buildResult);

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
       << std::setw(15) << score
       << std::setw(15) << succStr;
       
    for (double p : params) {
        ss << std::setw(15) << p; // Parameter einfügen
    }
    ss << stepResults; // Am Ende die Historie anheften
    qDebug() << "        | Result Line: " << QString::fromStdString(ss.str());
    return ss.str();
}

std::vector<PoseDef> SimulationManager::createPoseDefs()
{
    // 1. Erweitere die Header-Namen um die neue Spalte für die Knoten
    std::vector<std::string> jointNames = {"Wrist_F", "Wrist_A", "MCP_F", "MCP_A", "PIP", "DIP", "NumNodes"};
    
    // 2. Deine Parameter
    std::vector<double> numNodes = {20, 30, 50, 70};//{75, 100, 125, 150, 200};
    
    // 3. Temporäre Struktur für die Basis-Posen (ohne numNodes)
    struct BasePose {
        std::string name;
        std::vector<double> angles;
    };
    
    std::vector<BasePose> basePoses = {
        {"Full Fist",      { 0.0,   0.0, 90.0,  0.0, 100.0, 80.0}},
        {"Dach-Position",  { 0.0,   0.0, 90.0,  0.0,   0.0,  0.0}},
        {"Krallengriff",   { 0.0,   0.0,  0.0,  0.0, 100.0, 80.0}},
        {"Schraeger Griff",{ 0.0,   0.0, 45.0, 20.0,  45.0, 45.0}},
        {"Krampf Pose",    { 0.0,  80.0, 90.0,  0.0, 100.0, 80.0}},
        {"Power Grip",     {20.0,   0.0, 90.0,  0.0, 100.0, 80.0}},
        {"Dart-Wurf",      { 0.0, -60.0, 90.0,  0.0, 100.0, 80.0}}
    };
    
    std::vector<PoseDef> variantList;
    
    // 4. Verschachtelte Schleife: Kombiniere jede Pose mit jeder Knotenanzahl
    for (const auto& bp : basePoses) {
        for (double nodes : numNodes) {
            PoseDef p;
            
            // Name anpassen (z.B. "Full Fist (N=75)")
            p.PoseName = bp.name + " (N=" + std::to_string(static_cast<int>(nodes)) + ")";
            p.JointNames = jointNames;
            
            // Winkel übernehmen und am Ende die Knotenanzahl anhängen
            p.SimulationJointAngles = bp.angles;
            p.SimulationJointAngles.push_back(nodes); 
            
            variantList.push_back(p);
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


