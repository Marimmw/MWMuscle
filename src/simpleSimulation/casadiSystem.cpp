#include "casadiSystem.h"

#include <filesystem>


CasadiSystem::CasadiSystem(std::vector<SSMuscle*> muscles, int objType, std::string version, std::string parametrizationType, bool bUseCasGradient, bool bSumPhiEta, bool bUseWarmstartEtas, bool bDebug, bool bWriteFiles, double alpha)
    : m_muscles(muscles), objType(objType), Version(version), ParametrizationType(parametrizationType),
        bUseOwnGradient(bUseCasGradient), bSumPhiEta(bSumPhiEta), bUseWarmstartEtas(bUseWarmstartEtas), bDebug(bDebug), bWriteFiles(bWriteFiles), Alpha(alpha)
{
    CasadiSystemName = "CasSys_" + (m_muscles.empty() ? "Empty" : m_muscles[0]->Name);

    //setupCasadi();
    /* if (bSumPhiEta) {
        setupCasadiSum();
    }
    else {
        setupCasadi();
    } */
    if (bShowInfo) {
        qDebug() << "Initialized CasadiSystem: " << QString::fromStdString(CasadiSystemName);
        qDebug() << "  | Version: " << QString::fromStdString(Version);
        qDebug() << "  | ParamType: " << QString::fromStdString(ParametrizationType);
        qDebug() << "  | ObjType: " << objType;
        qDebug() << "  | bUseOwnGradient: " << bUseOwnGradient;
        qDebug() << "  | bSumPhiEta: " << bSumPhiEta;
        qDebug() << "  | bUseWarmstartEtas: " << bUseWarmstartEtas;
        qDebug() << "  | Alpha: " << Alpha;
    }

    setupCasadiViaSum();

    //setupCasadiToria();
    //setupCasadiSum();
    //setupCasadiViaSum_paramStudy();
    
}

void CasadiSystem::solveStepX(int step){
    /* if (bSumPhiEta) {
        solveStepSum();
    }
    else{solveStep();} */
    
    solveStepViaSum();

    //solveStepToria();
    //solveStepSum();
    //solveStepViaSum_paramStudy();
}

void CasadiSystem::setupCasadiSum()
{
    if (bDebug) qDebug() << "Setting up CasadiSystem with Sum Phi*Eta formulation...";
    using namespace casadi;
    MX all_x = MX::vertcat({});
    MX all_g = MX::vertcat({});
    MX all_p = MX::vertcat({});
    MX obj = 0;
    M = static_cast<int>(m_muscles.size());

    MX p_penalty = MX::sym("p_penalty", 1);  // <-- hier

    for (size_t m = 0; m < m_muscles.size(); ++m) {
        SSMuscle* mus = m_muscles[m];
        int num_inner = mus->MNodes.size() - 2;
        int K = mus->MNodes.size()-1;
        if (bDebug) qDebug() << "Muscle " << m << ": " << mus->Name.c_str() << ", K=" << K << ", num_inner=" << num_inner;
        int num_wrap = mus->meshPtrs.size();
        int num_eta_per_node = num_wrap;

        // --- VARIABLEN: [Punkte (3*n), Etas (n * Meshes)] ---
        MX x_mus = MX::sym("x_" + std::to_string(m), num_inner * 3 + num_inner * num_eta_per_node);
        all_x = MX::vertcat({all_x, x_mus});

        // --- PARAMETER: [Meshes (12*m), Origin (3), Insertion (3)] ---
        MX p_mus = MX::sym("p_" + std::to_string(m), 12 * num_wrap + 6);
        all_p = MX::vertcat({all_p, p_mus});

        MX P_orig = p_mus(Slice(12 * num_wrap, 12 * num_wrap + 3));
        MX P_ins  = p_mus(Slice(12 * num_wrap + 3, 12 * num_wrap + 6));
        int eta_offset = num_inner * 3;

        // --- CONSTRAINTS FÜR JEDEN KNOTEN ---
        for (int k = 0; k < num_inner; ++k) {
            MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
            MX g_prev = (k == 0) ? P_orig : x_mus(Slice((k - 1) * 3, k * 3));
            MX g_next = (k == num_inner - 1) ? P_ins : x_mus(Slice((k + 1) * 3, (k + 2) * 3));

            MX total_contact_force = 0;
            MX sum_h_eta = 0;  // Summierte Komplementarität

            // Über alle Wrapping-Meshes iterieren
            for (int j = 0; j < num_wrap; ++j) {
                MX q_j = p_mus(Slice(j * 12, (j + 1) * 12));
                MX eta_kj = x_mus(eta_offset + k * num_eta_per_node + j);
                

                // ============================================================
                // NEU: TORUS PENALTY (nur wenn Mesh ein Torus ist)
                // ============================================================
                if (bUseTorusPenalty) {
                    auto torus = dynamic_cast<SSTorusMesh*>(mus->meshPtrs[j]);
                    if (torus) {
                        MX C = q_j(Slice(0, 3));
                        MX N = q_j(Slice(9, 12)); // lokale Z-Achse = 3. Spalte der Rotationsmatrix

                        MX diff = g_k - C;
                        MX dot_diff_N = MX::dot(diff, N);
                        MX diff_projected = diff - dot_diff_N * N;
                        MX rho = MX::norm_2(diff_projected);

                        // Distanz-Gewichtung per Gauß
                        MX dist_to_center = MX::norm_2(diff);
                        double sigma = torus->R * 2.0;
                        MX weight = MX::exp(-(dist_to_center * dist_to_center) / (2.0 * sigma * sigma));

                        // Penalty wenn rho > R (Knoten außen)
                        MX violation = rho - torus->R; // R = torus->R
                        MX penalty = MX::fmax(MX(0.0), violation);
                        obj += TorusPenaltyWeight * p_penalty * weight * penalty * penalty;
                    }
                }
                // ============================================================

                // Distanz-Constraint berechnen
                MX h = mus->meshPtrs[j]->constraintDistance(g_k, q_j);

                // Gradient berechnen
                MX grad_h;
                if (bUseOwnGradient){
                    grad_h = mus->meshPtrs[j]->constraintJacobian(g_k, q_j);
                }
                else {
                    MX grad_h_full = MX::gradient(h, x_mus);
                    grad_h = grad_h_full(Slice(k * 3, (k + 1) * 3));
                }

                // Kontaktkraft akkumulieren
                total_contact_force += eta_kj * grad_h;
                
                // Komplementaritätsterme akkumulieren
                sum_h_eta += h * eta_kj;
                
                // EINZELNE CONSTRAINTS hinzufügen:
                // 1. Nicht-Eindringen: h >= 0
                all_g = MX::vertcat({all_g, h});
                
                // 2. Nur Druckkraft: eta >= 0 (wird durch Variablenbounds gehandhabt)
            }

            // SUMMIERTE KOMPLEMENTARITÄT: sum(h*eta) = 0
            all_g = MX::vertcat({all_g, sum_h_eta});

            // EULER-LAGRANGE GLEICHUNG (3 Komponenten für x, y, z)
            MX eq_el = (- g_prev + 2.0*g_k - g_next)*K - (1.0/K)*total_contact_force;
            all_g = MX::vertcat({all_g, eq_el});
        }
    }

    if (bUseTorusPenalty) {
        qDebug() << "Using torus penalty with weight: " << TorusPenaltyWeight;
        all_p = MX::vertcat({all_p, p_penalty});
    }


    opts["print_time"] = false;
    opts["record_time"] = false;
    opts["ipopt.print_level"] = 0;
    opts["ipopt.max_iter"] = maxIterations;
    opts["ipopt.tol"] = maxTol;
    opts["ipopt.linear_solver"] = "mumps";
    opts["ipopt.hessian_approximation"] = "limited-memory"; // "limited-memory" or "exact"
    if (bWriteFiles) {
        std::string dirPath = "../examples/results";
        if (!std::filesystem::exists(dirPath)) {
            std::filesystem::create_directories(dirPath);
        }
        // 3. Datei-Pfad zusammensetzen
        std::string filename = dirPath + "/solver_log_" + CasadiSystemName + ".txt";
        opts["ipopt.output_file"] = filename; 
        opts["ipopt.file_print_level"] = 5;
    }

    // extra
    /* opts["ipopt.max_soc"] = 4; // WICHTIG: In CasADi heißt die Option oft "max_soc", nicht "max_soc_iter"!
    opts["ipopt.alpha_red_factor"] = 0.5; // Line Search aggressiver abbremsen
    opts["ipopt.accept_every_trial_step"] = "no";    // Zwingt IPOPT, die Line Search genauer zu nehmen */

    MXDict nlp = {{"x", all_x}, {"f", obj}, {"g", all_g}, {"p", all_p}};
    solver = nlpsol("f", "ipopt", nlp, opts);

    if (bDebug) qDebug() << "CasadiSystem initialized with"
             << M << "muscles,"
             << "Objective Type:" << objType
             << "Total variables:" << all_x.size1()
             << "Total constraints:" << all_g.size1();
}

void CasadiSystem::solveStepSum() {
    if (bDebug) qDebug() << "Solving step with Sum Phi*Eta formulation...";
    using namespace casadi;

    std::vector<double> x0_all, p_all, lbg_all, ubg_all, lbx_all, ubx_all;
    std::vector<double> inputParamsForDebug;
    std::vector<std::string> inputParamDescriptionsForDebug;

    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;
        int num_wrap = mus->meshPtrs.size();
        int num_etas_total = num_inner * num_wrap;

        // --- INITIAL GUESS (x0) ---
        for (size_t k = 1; k < mus->MNodes.size() - 1; ++k) {
            MWMath::Point3D pWarm = mus->MNodes[k].predictNewGlobal();
            x0_all.push_back(pWarm.x);
            x0_all.push_back(pWarm.y);
            x0_all.push_back(pWarm.z);
        }
        
        if (mus->lastEtas.size() == (size_t)num_etas_total && bUseWarmstartEtas) {
            if (bDebug) qDebug() << "    Warmstart etas:" << QString::fromStdString(mus->Name);
            std::vector<double> scaledEtas;
            for (double eta : mus->lastEtas) {
                double scaledEta = eta * WarmstartEtaScaling;
                scaledEtas.push_back(scaledEta);
            }
            x0_all.insert(x0_all.end(), scaledEtas.begin(), scaledEtas.end());
        } else {
            if (bDebug) qDebug() << "    No warmstart etas:" << QString::fromStdString(mus->Name);
            for (int e = 0; e < num_etas_total; ++e) x0_all.push_back(0.0);
        }

        // --- PARAMETER (p) ---
        for (auto* mesh : mus->meshPtrs) {
            p_all.push_back(mesh->PositionGlobal.x);    inputParamsForDebug.push_back(mesh->PositionGlobal.x);  inputParamDescriptionsForDebug.push_back(mesh->Name + "_PosX");
            p_all.push_back(mesh->PositionGlobal.y);    inputParamsForDebug.push_back(mesh->PositionGlobal.y);  inputParamDescriptionsForDebug.push_back(mesh->Name + "_PosY");
            p_all.push_back(mesh->PositionGlobal.z);    inputParamsForDebug.push_back(mesh->PositionGlobal.z);  inputParamDescriptionsForDebug.push_back(mesh->Name + "_PosZ");
            
            MWMath::RotMatrix3x3 R = mesh->OrientationGlobal;
            for (int col = 0; col < 3; ++col) {
                for (int row = 0; row < 3; ++row) {
                    p_all.push_back(R.m[row][col]);
                    inputParamsForDebug.push_back(R.m[row][col]);
                    inputParamDescriptionsForDebug.push_back(mesh->Name + "_Ori_" + std::to_string(row) + std::to_string(col));
                }
            }
        }
        p_all.push_back(mus->OriginPointGlobal.x);    inputParamsForDebug.push_back(mus->OriginPointGlobal.x);  inputParamDescriptionsForDebug.push_back(mus->Name + "_OriginX");
        p_all.push_back(mus->OriginPointGlobal.y);    inputParamsForDebug.push_back(mus->OriginPointGlobal.y);  inputParamDescriptionsForDebug.push_back(mus->Name + "_OriginY");
        p_all.push_back(mus->OriginPointGlobal.z);    inputParamsForDebug.push_back(mus->OriginPointGlobal.z);  inputParamDescriptionsForDebug.push_back(mus->Name + "_OriginZ");
        p_all.push_back(mus->InsertionPointGlobal.x); inputParamsForDebug.push_back(mus->InsertionPointGlobal.x);  inputParamDescriptionsForDebug.push_back(mus->Name + "_InsertionX");
        p_all.push_back(mus->InsertionPointGlobal.y); inputParamsForDebug.push_back(mus->InsertionPointGlobal.y);  inputParamDescriptionsForDebug.push_back(mus->Name + "_InsertionY");
        p_all.push_back(mus->InsertionPointGlobal.z); inputParamsForDebug.push_back(mus->InsertionPointGlobal.z);  inputParamDescriptionsForDebug.push_back(mus->Name + "_InsertionZ");
        
        // --- CONSTRAINT BOUNDS (lbg/ubg) ---
        // Für jeden inneren Knoten:
        for (int k = 0; k < num_inner; ++k) {
            
            // 1. Nicht-Eindringen für jedes Mesh: h >= 0
            for (int j = 0; j < num_wrap; ++j) {
                lbg_all.push_back(0.0);    // h >= 0
                ubg_all.push_back(inf);
            }
            
            // 2. Eine summierte Komplementarität: sum(h*eta) = 0
            lbg_all.push_back(0.0);
            ubg_all.push_back(0.0);
            
            // 3. Euler-Lagrange Gleichung (3 Komponenten: x, y, z) = 0
            for (int d = 0; d < 3; ++d) { 
                lbg_all.push_back(0.0); 
                ubg_all.push_back(0.0); 
            }
        }

        // --- VARIABLE BOUNDS (lbx/ubx) ---
        // Gamma (Punkte) sind unbeschränkt
        for (int i = 0; i < num_inner * 3; ++i) {
            lbx_all.push_back(-inf); 
            ubx_all.push_back(inf);
        }
        // Eta (Kräfte) >= 0
        for (int i = 0; i < num_etas_total; ++i) {
            lbx_all.push_back(0.0); 
            ubx_all.push_back(inf);
        }
    }

    if (bUseTorusPenalty){
        double currentPenaltyWeight = (Step == 0) ? TorusPenaltyWeight : 0.0;
        qDebug() << "Torus penalty weight for this step: " << currentPenaltyWeight;
        p_all.push_back(currentPenaltyWeight);
    }

    allParameterInputsAllSteps.push_back(inputParamsForDebug);
    allParameterInputDescriptionsAllSteps.push_back(inputParamDescriptionsForDebug);

    // --- SOLVER CALL ---
    DMDict arg = {{"x0", x0_all}, {"p", p_all}, {"lbg", lbg_all}, {"ubg", ubg_all}, {"lbx", lbx_all}, {"ubx", ubx_all}};
    DMDict res = solver(arg);

    // --- EXTRACT RESULTS ---
    std::vector<double> res_x = std::vector<double>(res.at("x"));
    auto solverInfo = solver.stats();
    std::string status = solverInfo.at("return_status");

    const std::string C_RED   = "\033[31m";
    const std::string C_GREEN = "\033[32m";
    const std::string C_RESET = "\033[0m";
    int convSteps = int(solverInfo.at("iter_count"));
    if (true) { // Dein if(true) oder if(bDebug)
        std::string stepColor = (status == "Solve_Succeeded") ? C_GREEN : C_RED;
        std::string fullMessage = stepColor + "                Solver finished after " + std::to_string(convSteps) + " iterations" + "(" + status + ")" + C_RESET;
        qDebug().noquote() << QString::fromStdString(fullMessage);
    }
    SolverConvergenceMessages.push_back(QString::fromStdString(status).toStdString());
    SolverConvergenceSteps.push_back(convSteps);
    
    int current_x_offset = 0;
    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;
        int num_wrap = mus->meshPtrs.size();
        int num_etas_mus = num_inner * num_wrap;

        // Gamma (Punkte) zurückschreiben
        for (int i = 0; i < num_inner; ++i) {
            int mIdx = i + 1;
            MWMath::Point3D p;
            p.x = res_x[current_x_offset + i * 3 + 0];
            p.y = res_x[current_x_offset + i * 3 + 1];
            p.z = res_x[current_x_offset + i * 3 + 2];

            mus->MusclePointsGlobal[mIdx] = p;
            mus->MNodes[mIdx].PositionGlobal = p;
            mus->MNodes[mIdx].updateLocalFrame();
        }

        // Etas für Warmstart speichern
        int eta_start_idx = current_x_offset + (num_inner * 3);
        mus->lastEtas.assign(
            res_x.begin() + eta_start_idx, 
            res_x.begin() + eta_start_idx + num_etas_mus
        );

        // ETAS IN LISTE SCHREIBEN (für Export)
        int num_eta_per_node = num_wrap;
        for (int k = 0; k < num_inner; ++k) {
            int mIdx = k + 1;
            
            std::vector<double> currentStepEtas;
            
            for (int e = 0; e < num_eta_per_node; ++e) {
                int globalEtaIndex = k * num_eta_per_node + e;
                if (globalEtaIndex < mus->lastEtas.size()) {
                    currentStepEtas.push_back(mus->lastEtas[globalEtaIndex]);
                } else {
                    currentStepEtas.push_back(0.0);
                }
            }

            mus->MNodes[mIdx].MNodeEtaSteps.push_back(currentStepEtas);
            mus->MNodes[mIdx].lastEtas = mus->lastEtas; 
        }

        for (auto& node : mus->MNodes) {
            node.lastEtas = mus->lastEtas;
        }

        current_x_offset += (num_inner * 3 + num_etas_mus);
    }

    Step++;
}


// VIA POINTS
void CasadiSystem::solveStepVia()
{
    qDebug() << "          Solving step with Via Point formulation...";
    using namespace casadi;

    std::vector<double> x0_all, p_all, lbg_all, ubg_all, lbx_all, ubx_all;

    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;

        // ==========================================================
        // 0. MESHES SORTIEREN
        // ==========================================================
        std::vector<SSMesh*> sig_meshes;
        std::vector<SSMesh*> via_meshes;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) via_meshes.push_back(mesh);
            else sig_meshes.push_back(mesh);
        }
        
        int num_sig = sig_meshes.size();
        int num_via = via_meshes.size();
        int num_etas_sig = num_inner * num_sig;
        int num_etas_via = num_via;
        int num_etas_total = num_etas_sig + num_etas_via;

        // ==========================================================
        // 1. INITIAL GUESS (x0)
        // ==========================================================
        for (size_t k = 1; k < mus->MNodes.size() - 1; ++k) {
            MWMath::Point3D pWarm = mus->MNodes[k].predictNewGlobal();
            x0_all.push_back(pWarm.x);
            x0_all.push_back(pWarm.y);
            x0_all.push_back(pWarm.z);
        }
        
        // Warmstart für BEIDE Eta-Arten
        if (mus->lastEtas.size() == (size_t)num_etas_total && bUseWarmstartEtas) {
            for (double eta : mus->lastEtas) {
                x0_all.push_back(eta * WarmstartEtaScaling);
            }
        } else {
            for (int e = 0; e < num_etas_total; ++e) x0_all.push_back(0.0);
        }

        // ==========================================================
        // 2. PARAMETER (p)
        // ==========================================================
        // A. Signorini Meshes (12 Parameter: Pos + Rot Matrix)
        for (auto* mesh : sig_meshes) {
            p_all.push_back(mesh->PositionGlobal.x);
            p_all.push_back(mesh->PositionGlobal.y);
            p_all.push_back(mesh->PositionGlobal.z);
            
            MWMath::RotMatrix3x3 R = mesh->OrientationGlobal;
            for (int col = 0; col < 3; ++col) {
                for (int row = 0; row < 3; ++row) {
                    p_all.push_back(R.m[row][col]);
                }
            }
        }
        
        // B. Via-Point Meshes (Nur 3 Parameter: Position)
        for (auto* mesh : via_meshes) {
            p_all.push_back(mesh->PositionGlobal.x);
            p_all.push_back(mesh->PositionGlobal.y);
            p_all.push_back(mesh->PositionGlobal.z);
        }
        
        // C. Origin & Insertion
        p_all.push_back(mus->OriginPointGlobal.x);
        p_all.push_back(mus->OriginPointGlobal.y);
        p_all.push_back(mus->OriginPointGlobal.z);
        
        p_all.push_back(mus->InsertionPointGlobal.x);
        p_all.push_back(mus->InsertionPointGlobal.y);
        p_all.push_back(mus->InsertionPointGlobal.z);

        // ==========================================================
        // 3. CONSTRAINT BOUNDS (lbg/ubg)
        // ==========================================================
        // A. Globale Via-Points
        for (int v = 0; v < num_via; ++v) {
            // Komplementarität (Unterdruck-Trick)
            lbg_all.push_back(-inf); 
            ubg_all.push_back(0.0); // evtl Leichtes Spiel für Robustheit
            // Distanz-Constraint: h_v >= 0
            lbg_all.push_back(0.0); 
            ubg_all.push_back(inf);
        }

        // B. Knoten-Schleife (Signorini + Euler-Lagrange)
        for (int k = 0; k < num_inner; ++k) {
            // Non-penetration für jedes Signorini-Mesh (h_s >= 0)
            for (int s = 0; s < num_sig; ++s) {
                lbg_all.push_back(0.0); 
                ubg_all.push_back(inf);
            }
            // Summierte Signorini-Komplementarität (sum(h*eta) == 0)
            if (num_sig > 0) {
                lbg_all.push_back(0.0); 
                ubg_all.push_back(0.0);
            }
            // Euler-Lagrange Kraftbilanz (x, y, z == 0)
            for (int d = 0; d < 3; ++d) { 
                lbg_all.push_back(0.0); 
                ubg_all.push_back(0.0); 
            }
        }

        // ==========================================================
        // 4. VARIABLE BOUNDS (lbx/ubx)
        // ==========================================================
        // Gamma (Punkte)
        for (int i = 0; i < num_inner * 3; ++i) {
            lbx_all.push_back(-inf); ubx_all.push_back(inf);
        }
        // Alle Etas (Signorini + Via-Points) dürfen nur ziehen/drücken (>= 0)
        for (int i = 0; i < num_etas_total; ++i) {
            lbx_all.push_back(0.0); ubx_all.push_back(inf);
        }
    }

    // ==========================================================
    // SOLVER AUFRUF
    // ==========================================================
    DMDict arg = {{"x0", x0_all}, {"p", p_all}, {"lbg", lbg_all}, {"ubg", ubg_all}, {"lbx", lbx_all}, {"ubx", ubx_all}};
    DMDict res = solver(arg);

    // ==========================================================
    // ERGEBNISSE EXTRAHIEREN
    // ==========================================================
    std::vector<double> res_x = std::vector<double>(res.at("x"));
    auto solverInfo = solver.stats();
    std::string status = solverInfo.at("return_status");

    const std::string C_RED   = "\033[31m";
    const std::string C_GREEN = "\033[32m";
    const std::string C_RESET = "\033[0m";
    int convSteps = int(solverInfo.at("iter_count"));
    if (true) { // Dein if(true) oder if(bDebug)
        std::string stepColor = (status == "Solve_Succeeded") ? C_GREEN : C_RED;
        std::string fullMessage = stepColor + "    Solver finished after " + std::to_string(convSteps) + " iterations" + "(" + status + ")" + C_RESET;
        qDebug().noquote() << QString::fromStdString(fullMessage);
    }
    
    SolverConvergenceMessages.push_back(status);
    SolverConvergenceSteps.push_back(convSteps);
    
    int current_x_offset = 0;
    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;
        int num_sig = 0, num_via = 0;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) num_via++; else num_sig++;
        }
        
        int offset_nodes    = current_x_offset;
        int offset_sig_etas = offset_nodes + (num_inner * 3);
        int offset_via_etas = offset_sig_etas + (num_inner * num_sig);

        // 1. Positionen zurückschreiben
        for (int i = 0; i < num_inner; ++i) {
            int mIdx = i + 1;
            MWMath::Point3D p(
                res_x[offset_nodes + i * 3 + 0],
                res_x[offset_nodes + i * 3 + 1],
                res_x[offset_nodes + i * 3 + 2]
            );
            mus->MusclePointsGlobal[mIdx] = p;
            mus->MNodes[mIdx].PositionGlobal = p;
            mus->MNodes[mIdx].updateLocalFrame();
        }

        // 2. Globale Etas für Warmstart speichern
        int num_etas_total = (num_inner * num_sig) + num_via;
        mus->lastEtas.assign(res_x.begin() + offset_sig_etas, res_x.begin() + offset_sig_etas + num_etas_total);

        // ==========================================================
        // 3. MAPPING-TRICK FÜR EXPORT/VISUALISIERUNG
        // Wir dröseln die Ergebnisse für die alte Struktur wieder auf!
        // ==========================================================
        
        // Zuerst: Welcher Knoten ist jedem Via-Point am nächsten?
        std::vector<int> closest_node_for_via(num_via, -1);
        int via_counter = 0;
        for (auto* mesh : mus->meshPtrs) {
            if (!mesh->bIsViaPoint) continue;
            
            MWMath::Point3D v_pos = mesh->PositionGlobal;
            double min_dist = 1e9;
            int best_k = 0;
            
            for (int k = 0; k < num_inner; ++k) {
                double d = MWMath::distance(mus->MNodes[k + 1].PositionGlobal, v_pos);
                if (d < min_dist) { min_dist = d; best_k = k; }
            }
            closest_node_for_via[via_counter] = best_k;
            via_counter++;
        }

        // Nun füllen wir das Array für JEDEN Knoten
        for (int k = 0; k < num_inner; ++k) {
            int mIdx = k + 1;
            std::vector<double> currentStepEtas(mus->meshPtrs.size(), 0.0);

            int sig_idx = 0;
            int via_idx = 0;

            for (size_t m = 0; m < mus->meshPtrs.size(); ++m) {
                if (mus->meshPtrs[m]->bIsViaPoint) {
                    // Via-Point Kraft nur eintragen, wenn das der nächste Knoten ist
                    if (closest_node_for_via[via_idx] == k) {
                        currentStepEtas[m] = res_x[offset_via_etas + via_idx];
                    }
                    via_idx++;
                } else {
                    // Signorini Kraft ganz regulär auslesen
                    currentStepEtas[m] = res_x[offset_sig_etas + (k * num_sig) + sig_idx];
                    sig_idx++;
                }
            }
            
            mus->MNodes[mIdx].MNodeEtaSteps.push_back(currentStepEtas);
            mus->MNodes[mIdx].lastEtas = mus->lastEtas; 
        }

        current_x_offset += (num_inner * 3) + num_etas_total;
    }
}

void CasadiSystem::setupCasadiVia()
{
    qDebug() << "     Setting up CasadiSystem with Via Point formulation...";
    using namespace casadi;
    MX all_x = MX::vertcat({});
    MX all_g = MX::vertcat({});
    MX all_p = MX::vertcat({});
    MX obj = 0;
    
    

    for (size_t m = 0; m < m_muscles.size(); ++m) {
        SSMuscle* mus = m_muscles[m];
        int num_inner = mus->MNodes.size() - 2;
        int K = mus->MNodes.size() - 1;

        // ==========================================================
        // 0. MESHES SORTIEREN
        // ==========================================================
        std::vector<SSMesh*> sig_meshes;
        std::vector<SSMesh*> via_meshes;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) via_meshes.push_back(mesh);
            else sig_meshes.push_back(mesh);
        }
        int num_sig = sig_meshes.size();
        int num_via = via_meshes.size();

        // ==========================================================
        // 1. VARIABLEN DEFINIEREN (X-Vektor)
        // Aufbau: [ 3*N Knoten | N * M_sig Etas | M_via Etas ]
        // ==========================================================
        int num_etas_sig = num_inner * num_sig;
        int num_etas_via = num_via;
        
        MX x_mus = MX::sym("x_" + std::to_string(m), (num_inner * 3) + num_etas_sig + num_etas_via);
        all_x = MX::vertcat({all_x, x_mus});

        // Offsets zum Herausschneiden aus x_mus
        int offset_nodes    = 0;
        int offset_sig_etas = num_inner * 3;
        int offset_via_etas = offset_sig_etas + num_etas_sig;

        // ==========================================================
        // 2. PARAMETER DEFINIEREN (P-Vektor)
        // Aufbau: [ 12*M_sig Parameter | 3*M_via Parameter | 3 Orig | 3 Ins ]
        // ==========================================================
        MX p_mus = MX::sym("p_" + std::to_string(m), (12 * num_sig) + (3 * num_via) + 6);
        all_p = MX::vertcat({all_p, p_mus});

        // Offsets zum Herausschneiden aus p_mus
        int p_offset_sig  = 0;
        int p_offset_via  = num_sig * 12;
        int p_offset_ends = p_offset_via + (num_via * 3);

        MX P_orig = p_mus(Slice(p_offset_ends, p_offset_ends + 3));
        MX P_ins  = p_mus(Slice(p_offset_ends + 3, p_offset_ends + 6));

        // ==========================================================
        // 3. GLOBALE VIA-POINT BERECHNUNG (Wie bisher)
        // ==========================================================
        std::vector<MX> h_via_list; 
        // Parameter für Via-Points
        // double alpha = 100.0; // bei 50.0 ist der punkt weiter weg vom Mesh, sinnvoll ist 50-100 ab 300 "Invalid_Number_Detected" fehler
        double eps = 1e-8;
        for (int v = 0; v < num_via; ++v) {
            MX v_pos = p_mus(Slice(p_offset_via + v * 3, p_offset_via + (v + 1) * 3));
            MX eta_v = x_mus(offset_via_etas + v);
            
            // HIER HOLEN WIR UNS DEINE ECHTE C++ TOLERANZ!
            double r_tol = via_meshes[v]->MViaPointTolerance * 0.5; // ohne "0.5" irgendwie abstand zu groß 
            // Fallback
            if (r_tol <= 0.0001) r_tol = 0.005; 
            double r_sq = r_tol * r_tol;

            MX d_vector = MX::vertcat({});
            // for (int k = 0; k < num_inner; ++k) {
            //     MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
            //     d_vector = MX::vertcat({d_vector, sum1(sq(g_k - v_pos))});
            // }
            for (int k = 0; k < num_inner; ++k) {
                MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
                // LINEARE DISTANZ: sqrt( (dx)^2 + (dy)^2 + (dz)^2 + eps )
                MX dist = sqrt(sum1(sq(g_k - v_pos)) + eps);
                d_vector = MX::vertcat({d_vector, dist});
            }
            // geändert am 10.04.26
            /* d_vector = MX::vertcat({d_vector, sum1(sq(P_orig - v_pos))});
            d_vector = MX::vertcat({d_vector, sum1(sq(P_ins - v_pos))}); */
            d_vector = MX::vertcat({d_vector, sqrt(sum1(sq(P_orig - v_pos)) + eps)});
            d_vector = MX::vertcat({d_vector, sqrt(sum1(sq(P_ins - v_pos)) + eps)});

            MX D_v, h_v;
            if (bMaxLogSumTrick){
                // --- NUMERISCH STABILES LOGSUMEXP (MAX-TRICK) ---
                MX x = -Alpha * d_vector;
                // 1. Finde den größten Wert im Vektor (mmax)
                MX x_max = casadi::MX::mmax(x); 
                // 2. Ziehe das Maximum vor dem exp() ab. 
                MX sum_exp = sum1(exp(x - x_max));
                // 3. Logarithmus ziehen und das Maximum wieder addieren
                D_v = -(x_max + log(sum_exp)) / Alpha;
                h_v = r_tol - D_v; 
            }
            else{
                D_v = -casadi::MX::logsumexp(-Alpha * d_vector) / Alpha;
                h_v = r_sq - D_v; // hier eventuell linear statt quadratisch abziehen -> (quadratuisch - linear funkioniert aber irgendiwe besser?)
            }
            
            h_via_list.push_back(h_v);

            // Constraints für Via-Points
            all_g = MX::vertcat({all_g, h_v * eta_v}); 
            all_g = MX::vertcat({all_g, h_v});         
        }

        // #####################################

        // Gradienten für Via-Points vorab berechnen
        std::vector<MX> grad_h_via_full_list;
        for (int v = 0; v < num_via; ++v) {
            grad_h_via_full_list.push_back(MX::gradient(h_via_list[v], x_mus));
        }

        // ==========================================================
        // 4. KNOTEN-SCHLEIFE (Signorini & Euler-Lagrange)
        // ==========================================================
        for (int k = 0; k < num_inner; ++k) {
            MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
            MX g_prev = (k == 0) ? P_orig : x_mus(Slice((k - 1) * 3, k * 3));
            MX g_next = (k == num_inner - 1) ? P_ins : x_mus(Slice((k + 1) * 3, (k + 2) * 3));

            MX total_sig_force = 0;
            MX total_via_force = 0;
            MX sum_h_eta_sig = 0; 

            // --- A. SIGNORINI KRÄFTE (Klassisch aufsummiert) ---
            for (int s = 0; s < num_sig; ++s) {
                MX q_s = p_mus(Slice(p_offset_sig + s * 12, p_offset_sig + (s + 1) * 12));
                MX eta_ks = x_mus(offset_sig_etas + k * num_sig + s);
                
                MX h_s = sig_meshes[s]->constraintDistance(g_k, q_s);

                MX grad_h_s;
                if (bUseOwnGradient) grad_h_s = sig_meshes[s]->constraintJacobian(g_k, q_s);
                else {
                    MX grad_full = MX::gradient(h_s, x_mus);
                    grad_h_s = grad_full(Slice(k * 3, (k + 1) * 3));
                }

                total_sig_force += eta_ks * grad_h_s;
                sum_h_eta_sig += h_s * eta_ks;
                
                all_g = MX::vertcat({all_g, h_s}); // Non-penetration: h >= 0
            }

            // Summierte Signorini-Komplementarität hinzufügen
            if (num_sig > 0) {
                all_g = MX::vertcat({all_g, sum_h_eta_sig}); // sum(h*eta) == 0
            }

            // --- B. VIA-POINT KRÄFTE ---
            for (int v = 0; v < num_via; ++v) {
                MX eta_v = x_mus(offset_via_etas + v);
                MX grad_h_v = grad_h_via_full_list[v](Slice(k * 3, (k + 1) * 3)); 
                
                total_via_force += eta_v * grad_h_v;
            }

            // --- C. EULER-LAGRANGE KRAFTBILANZ ---
            // Summiere BEIDE Kontaktkraft-Arten in die Lagrange-Gleichung!
            MX total_contact_force = total_sig_force + total_via_force;
            MX eq_el = (-g_prev + 2.0 * g_k - g_next) * K - (1.0 / K) * total_contact_force;
            all_g = MX::vertcat({all_g, eq_el});
        }
    }

    // ==========================================================
    // SOLVER CONFIG (Bleibt unverändert)
    // ==========================================================
    opts["print_time"] = false;
    opts["record_time"] = false;
    opts["ipopt.print_level"] = 0;
    opts["ipopt.max_iter"] = maxIterations;
    opts["ipopt.tol"] = maxTol;
    opts["ipopt.linear_solver"] = "mumps";
    opts["ipopt.hessian_approximation"] = "limited-memory"; 
    
    // ... [Dein Code zum Schreiben der Datei hier] ...

    MXDict nlp = {{"x", all_x}, {"f", obj}, {"g", all_g}, {"p", all_p}};
    solver = nlpsol("f", "ipopt", nlp, opts);

    if (bDebug) qDebug() << "CasadiSystem initialized with"
             << M << "muscles,"
             << "Total variables:" << all_x.size1()
             << "Total constraints:" << all_g.size1();
}



// VIA POINTS SUM
void CasadiSystem::solveStepViaSum()
{
    qDebug() << "          Solving step with Via Point Sum formulation...";
    using namespace casadi;

    std::vector<double> x0_all, p_all, lbg_all, ubg_all, lbx_all, ubx_all;

    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;

        // ==========================================================
        // 0. MESHES SORTIEREN
        // ==========================================================
        std::vector<SSMesh*> sig_meshes;
        std::vector<SSMesh*> via_meshes;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) via_meshes.push_back(mesh);
            else sig_meshes.push_back(mesh);
        }
        
        int num_sig = sig_meshes.size();
        int num_via = via_meshes.size();
        int num_etas_sig = num_inner * num_sig;
        int num_etas_via = num_via;
        int num_etas_total = num_etas_sig + num_etas_via;

        // ==========================================================
        // 1. INITIAL GUESS (x0)
        // ==========================================================
        for (size_t k = 1; k < mus->MNodes.size() - 1; ++k) {
            MWMath::Point3D pWarm = mus->MNodes[k].predictNewGlobal();
            x0_all.push_back(pWarm.x);
            x0_all.push_back(pWarm.y);
            x0_all.push_back(pWarm.z);
        }
        
        // Warmstart für BEIDE Eta-Arten
        if (mus->lastEtas.size() == (size_t)num_etas_total && bUseWarmstartEtas) {
            for (double eta : mus->lastEtas) {
                x0_all.push_back(eta * WarmstartEtaScaling);
            }
        } else {
            for (int e = 0; e < num_etas_total; ++e) x0_all.push_back(0.0);
        }

        // ==========================================================
        // 2. PARAMETER (p)
        // ==========================================================
        // A. Signorini Meshes (12 Parameter: Pos + Rot Matrix)
        for (auto* mesh : sig_meshes) {
            p_all.push_back(mesh->PositionGlobal.x);
            p_all.push_back(mesh->PositionGlobal.y);
            p_all.push_back(mesh->PositionGlobal.z);
            
            MWMath::RotMatrix3x3 R = mesh->OrientationGlobal;
            for (int col = 0; col < 3; ++col) {
                for (int row = 0; row < 3; ++row) {
                    p_all.push_back(R.m[row][col]);
                }
            }
        }
        
        // B. Via-Point Meshes (Nur 3 Parameter: Position)
        for (auto* mesh : via_meshes) {
            p_all.push_back(mesh->PositionGlobal.x);
            p_all.push_back(mesh->PositionGlobal.y);
            p_all.push_back(mesh->PositionGlobal.z);
        }
        
        // C. Origin & Insertion
        p_all.push_back(mus->OriginPointGlobal.x);
        p_all.push_back(mus->OriginPointGlobal.y);
        p_all.push_back(mus->OriginPointGlobal.z);
        
        p_all.push_back(mus->InsertionPointGlobal.x);
        p_all.push_back(mus->InsertionPointGlobal.y);
        p_all.push_back(mus->InsertionPointGlobal.z);

        // ==========================================================
        // 3. CONSTRAINT BOUNDS (lbg/ubg)
        // ==========================================================
        // A. Globale Via-Points
        for (int v = 0; v < num_via; ++v) {
            // NUR NOCH Distanz-Constraint: h_v >= 0
            lbg_all.push_back(0.0); 
            ubg_all.push_back(inf);
        }
        // Summierte Via-Point-Komplementarität (sum(h_v * eta_v) == 0)
        if (num_via > 0) {
            lbg_all.push_back(-inf); // Exakt 0.0, wie von dir gewünscht!
            ubg_all.push_back(0.0);
        }


        // B. Knoten-Schleife (Signorini + Euler-Lagrange)
        for (int k = 0; k < num_inner; ++k) {
            // Non-penetration für jedes Signorini-Mesh (h_s >= 0)
            for (int s = 0; s < num_sig; ++s) {
                lbg_all.push_back(0.0); 
                ubg_all.push_back(inf);
            }
            // Summierte Signorini-Komplementarität (sum(h*eta) == 0)
            if (num_sig > 0) {
                lbg_all.push_back(0.0); 
                ubg_all.push_back(0.0);
            }
            // Euler-Lagrange Kraftbilanz (x, y, z == 0)
            for (int d = 0; d < 3; ++d) { 
                lbg_all.push_back(0.0); 
                ubg_all.push_back(0.0); 
            }
        }

        // ==========================================================
        // 4. VARIABLE BOUNDS (lbx/ubx)
        // ==========================================================
        // Gamma (Punkte)
        for (int i = 0; i < num_inner * 3; ++i) {
            lbx_all.push_back(-inf); ubx_all.push_back(inf);
        }
        // Alle Etas (Signorini + Via-Points) dürfen nur ziehen/drücken (>= 0)
        for (int i = 0; i < num_etas_total; ++i) {
            lbx_all.push_back(0.0); ubx_all.push_back(inf);
        }
    }

    // ==========================================================
    // SOLVER AUFRUF
    // ==========================================================
    DMDict arg = {{"x0", x0_all}, {"p", p_all}, {"lbg", lbg_all}, {"ubg", ubg_all}, {"lbx", lbx_all}, {"ubx", ubx_all}};
    DMDict res = solver(arg);

    // ==========================================================
    // ERGEBNISSE EXTRAHIEREN
    // ==========================================================
    std::vector<double> res_x = std::vector<double>(res.at("x"));
    auto solverInfo = solver.stats();
    std::string status = solverInfo.at("return_status");

    const std::string C_RED   = "\033[31m";
    const std::string C_GREEN = "\033[32m";
    const std::string C_RESET = "\033[0m";
    int convSteps = int(solverInfo.at("iter_count"));
    if (true) { // Dein if(true) oder if(bDebug)
        std::string stepColor = (status == "Solve_Succeeded") ? C_GREEN : C_RED;
        std::string fullMessage = stepColor + "    Solver finished after " + std::to_string(convSteps) + " iterations" + "(" + status + ")" + C_RESET;
        qDebug().noquote() << QString::fromStdString(fullMessage);
    }
    
    SolverConvergenceMessages.push_back(status);
    SolverConvergenceSteps.push_back(convSteps);
    
    int current_x_offset = 0;
    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;
        int num_sig = 0, num_via = 0;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) num_via++; else num_sig++;
        }
        
        int offset_nodes    = current_x_offset;
        int offset_sig_etas = offset_nodes + (num_inner * 3);
        int offset_via_etas = offset_sig_etas + (num_inner * num_sig);

        // 1. Positionen zurückschreiben
        for (int i = 0; i < num_inner; ++i) {
            int mIdx = i + 1;
            MWMath::Point3D p(
                res_x[offset_nodes + i * 3 + 0],
                res_x[offset_nodes + i * 3 + 1],
                res_x[offset_nodes + i * 3 + 2]
            );
            mus->MusclePointsGlobal[mIdx] = p;
            mus->MNodes[mIdx].PositionGlobal = p;
            mus->MNodes[mIdx].updateLocalFrame();
        }

        // 2. Globale Etas für Warmstart speichern
        int num_etas_total = (num_inner * num_sig) + num_via;
        mus->lastEtas.assign(res_x.begin() + offset_sig_etas, res_x.begin() + offset_sig_etas + num_etas_total);

        // ==========================================================
        // 3. MAPPING-TRICK FÜR EXPORT/VISUALISIERUNG
        // Wir dröseln die Ergebnisse für die alte Struktur wieder auf!
        // ==========================================================
        
        // Zuerst: Welcher Knoten ist jedem Via-Point am nächsten?
        std::vector<int> closest_node_for_via(num_via, -1);
        int via_counter = 0;
        for (auto* mesh : mus->meshPtrs) {
            if (!mesh->bIsViaPoint) continue;
            
            MWMath::Point3D v_pos = mesh->PositionGlobal;
            double min_dist = 1e9;
            int best_k = 0;
            
            for (int k = 0; k < num_inner; ++k) {
                double d = MWMath::distance(mus->MNodes[k + 1].PositionGlobal, v_pos);
                if (d < min_dist) { min_dist = d; best_k = k; }
            }
            closest_node_for_via[via_counter] = best_k;
            via_counter++;
        }

        // Nun füllen wir das Array für JEDEN Knoten
        for (int k = 0; k < num_inner; ++k) {
            int mIdx = k + 1;
            std::vector<double> currentStepEtas(mus->meshPtrs.size(), 0.0);

            int sig_idx = 0;
            int via_idx = 0;

            for (size_t m = 0; m < mus->meshPtrs.size(); ++m) {
                if (mus->meshPtrs[m]->bIsViaPoint) {
                    // Via-Point Kraft nur eintragen, wenn das der nächste Knoten ist
                    if (closest_node_for_via[via_idx] == k) {
                        currentStepEtas[m] = res_x[offset_via_etas + via_idx];
                    }
                    via_idx++;
                } else {
                    // Signorini Kraft ganz regulär auslesen
                    currentStepEtas[m] = res_x[offset_sig_etas + (k * num_sig) + sig_idx];
                    sig_idx++;
                }
            }
            
            mus->MNodes[mIdx].MNodeEtaSteps.push_back(currentStepEtas);
            mus->MNodes[mIdx].lastEtas = mus->lastEtas; 
        }

        current_x_offset += (num_inner * 3) + num_etas_total;
    }
}

void CasadiSystem::setupCasadiViaSum()
{
    qDebug() << "     Setting up CasadiSystem with Via Point formulation...";
    using namespace casadi;
    MX all_x = MX::vertcat({});
    MX all_g = MX::vertcat({});
    MX all_p = MX::vertcat({});
    MX obj = 0;
    
    

    for (size_t m = 0; m < m_muscles.size(); ++m) {
        SSMuscle* mus = m_muscles[m];
        int num_inner = mus->MNodes.size() - 2;
        int K = mus->MNodes.size() - 1;

        // ==========================================================
        // 0. MESHES SORTIEREN
        // ==========================================================
        std::vector<SSMesh*> sig_meshes;
        std::vector<SSMesh*> via_meshes;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) via_meshes.push_back(mesh);
            else sig_meshes.push_back(mesh);
        }
        int num_sig = sig_meshes.size();
        int num_via = via_meshes.size();

        // ==========================================================
        // 1. VARIABLEN DEFINIEREN (X-Vektor)
        // Aufbau: [ 3*N Knoten | N * M_sig Etas | M_via Etas ]
        // ==========================================================
        int num_etas_sig = num_inner * num_sig;
        int num_etas_via = num_via;
        
        MX x_mus = MX::sym("x_" + std::to_string(m), (num_inner * 3) + num_etas_sig + num_etas_via);
        all_x = MX::vertcat({all_x, x_mus});

        // Offsets zum Herausschneiden aus x_mus
        int offset_nodes    = 0;
        int offset_sig_etas = num_inner * 3;
        int offset_via_etas = offset_sig_etas + num_etas_sig;

        // ==========================================================
        // 2. PARAMETER DEFINIEREN (P-Vektor)
        // Aufbau: [ 12*M_sig Parameter | 3*M_via Parameter | 3 Orig | 3 Ins ]
        // ==========================================================
        MX p_mus = MX::sym("p_" + std::to_string(m), (12 * num_sig) + (3 * num_via) + 6);
        all_p = MX::vertcat({all_p, p_mus});

        // Offsets zum Herausschneiden aus p_mus
        int p_offset_sig  = 0;
        int p_offset_via  = num_sig * 12;
        int p_offset_ends = p_offset_via + (num_via * 3);

        MX P_orig = p_mus(Slice(p_offset_ends, p_offset_ends + 3));
        MX P_ins  = p_mus(Slice(p_offset_ends + 3, p_offset_ends + 6));

        // ==========================================================
        // 3. GLOBALE VIA-POINT BERECHNUNG (Wie bisher)
        // ==========================================================
        std::vector<MX> h_via_list; 
        MX sum_h_eta_via = 0;
        // Parameter für Via-Points
        // double alpha = 100.0; // bei 50.0 ist der punkt weiter weg vom Mesh, sinnvoll ist 50-100 ab 300 "Invalid_Number_Detected" fehler
        double eps = 1e-8;
        for (int v = 0; v < num_via; ++v) {
            MX v_pos = p_mus(Slice(p_offset_via + v * 3, p_offset_via + (v + 1) * 3));
            MX eta_v = x_mus(offset_via_etas + v);
            
            // HIER HOLEN WIR UNS DEINE ECHTE C++ TOLERANZ!
            double r_tol = via_meshes[v]->MViaPointTolerance * 0.5; // ohne "0.5" irgendwie abstand zu groß 
            // Fallback
            if (r_tol <= 0.0001) r_tol = 0.005; 
            double r_sq = r_tol * r_tol;

            MX d_vector = MX::vertcat({});
            // for (int k = 0; k < num_inner; ++k) {
            //     MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
            //     d_vector = MX::vertcat({d_vector, sum1(sq(g_k - v_pos))});
            // }
            for (int k = 0; k < num_inner; ++k) {
                MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
                // LINEARE DISTANZ: sqrt( (dx)^2 + (dy)^2 + (dz)^2 + eps )
                MX dist = sqrt(sum1(sq(g_k - v_pos)) + eps);
                d_vector = MX::vertcat({d_vector, dist});
            }
            // geändert am 10.04.26
            /* d_vector = MX::vertcat({d_vector, sum1(sq(P_orig - v_pos))});
            d_vector = MX::vertcat({d_vector, sum1(sq(P_ins - v_pos))}); */
            d_vector = MX::vertcat({d_vector, sqrt(sum1(sq(P_orig - v_pos)) + eps)});
            d_vector = MX::vertcat({d_vector, sqrt(sum1(sq(P_ins - v_pos)) + eps)});

            
            MX D_v, h_v;
            if (bMaxLogSumTrick){
                // --- NUMERISCH STABILES LOGSUMEXP (MAX-TRICK) ---
                MX x = -Alpha * d_vector;
                // 1. Finde den größten Wert im Vektor (mmax)
                MX x_max = casadi::MX::mmax(x); 
                // 2. Ziehe das Maximum vor dem exp() ab. 
                MX sum_exp = sum1(exp(x - x_max));
                // 3. Logarithmus ziehen und das Maximum wieder addieren
                D_v = -(x_max + log(sum_exp)) / Alpha;
                h_v = r_sq - D_v; 
            }
            else{
                D_v = -casadi::MX::logsumexp(-Alpha * d_vector) / Alpha;
                h_v = r_sq - D_v; // hier eventuell linear statt quadratisch abziehen -> (quadratuisch - linear funkioniert aber irgendiwe besser?)
            }
            
            
            h_via_list.push_back(h_v);

            sum_h_eta_via += h_v * eta_v;

            // Constraints für Via-Points: NUR NOCH DIE DISTANZ EINZELN HINZUFÜGEN
            all_g = MX::vertcat({all_g, h_v});     
        }
        // NEU: GLOBALES CONSTRAINT FÜR ALLE VIA-POINTS ANFÜGEN
        if (num_via > 0) {
            all_g = MX::vertcat({all_g, sum_h_eta_via});
        }
        // #####################################

        // Gradienten für Via-Points vorab berechnen
        std::vector<MX> grad_h_via_full_list;
        for (int v = 0; v < num_via; ++v) {
            grad_h_via_full_list.push_back(MX::gradient(h_via_list[v], x_mus));
        }

        // ==========================================================
        // 4. KNOTEN-SCHLEIFE (Signorini & Euler-Lagrange)
        // ==========================================================
        for (int k = 0; k < num_inner; ++k) {
            MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
            MX g_prev = (k == 0) ? P_orig : x_mus(Slice((k - 1) * 3, k * 3));
            MX g_next = (k == num_inner - 1) ? P_ins : x_mus(Slice((k + 1) * 3, (k + 2) * 3));

            MX total_sig_force = 0;
            MX total_via_force = 0;
            MX sum_h_eta_sig = 0; 

            // --- A. SIGNORINI KRÄFTE (Klassisch aufsummiert) ---
            for (int s = 0; s < num_sig; ++s) {
                MX q_s = p_mus(Slice(p_offset_sig + s * 12, p_offset_sig + (s + 1) * 12));
                MX eta_ks = x_mus(offset_sig_etas + k * num_sig + s);
                
                MX h_s = sig_meshes[s]->constraintDistance(g_k, q_s);

                MX grad_h_s;
                if (bUseOwnGradient) grad_h_s = sig_meshes[s]->constraintJacobian(g_k, q_s);
                else {
                    MX grad_full = MX::gradient(h_s, x_mus);
                    grad_h_s = grad_full(Slice(k * 3, (k + 1) * 3));
                }

                total_sig_force += eta_ks * grad_h_s;
                sum_h_eta_sig += h_s * eta_ks;
                
                all_g = MX::vertcat({all_g, h_s}); // Non-penetration: h >= 0
            }

            // Summierte Signorini-Komplementarität hinzufügen
            if (num_sig > 0) {
                all_g = MX::vertcat({all_g, sum_h_eta_sig}); // sum(h*eta) == 0
            }

            // --- B. VIA-POINT KRÄFTE ---
            for (int v = 0; v < num_via; ++v) {
                MX eta_v = x_mus(offset_via_etas + v);
                MX grad_h_v = grad_h_via_full_list[v](Slice(k * 3, (k + 1) * 3)); 
                
                total_via_force += eta_v * grad_h_v;
            }

            // --- C. EULER-LAGRANGE KRAFTBILANZ ---
            // Summiere BEIDE Kontaktkraft-Arten in die Lagrange-Gleichung!
            MX total_contact_force = total_sig_force + total_via_force;
            MX eq_el = (-g_prev + 2.0 * g_k - g_next) * K - (1.0 / K) * total_contact_force;
            all_g = MX::vertcat({all_g, eq_el});
        }
    }

    // ==========================================================
    // SOLVER CONFIG (Bleibt unverändert)
    // ==========================================================
    opts["print_time"] = false;
    opts["record_time"] = false;
    opts["ipopt.print_level"] = 0;
    opts["ipopt.max_iter"] = maxIterations;
    opts["ipopt.tol"] = maxTol;
    opts["ipopt.linear_solver"] = "mumps";
    opts["ipopt.hessian_approximation"] = "limited-memory"; 
    
    // ... [Dein Code zum Schreiben der Datei hier] ...

    MXDict nlp = {{"x", all_x}, {"f", obj}, {"g", all_g}, {"p", all_p}};
    solver = nlpsol("f", "ipopt", nlp, opts);

    if (bDebug) qDebug() << "CasadiSystem initialized with"
             << M << "muscles,"
             << "Total variables:" << all_x.size1()
             << "Total constraints:" << all_g.size1();
}


// TORUS-VIA-POINTS
/* 
// test gemini
void CasadiSystem::solveStepToria()
{
    qDebug() << "          Solving step with Via Point Sum formulation...";
    using namespace casadi;

    std::vector<double> x0_all, p_all, lbg_all, ubg_all, lbx_all, ubx_all;

    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;

        // ==========================================================
        // 0. MESHES SORTIEREN + TORUS-VIA-INFOS SAMMELN
        // ==========================================================
        std::vector<SSMesh*> sig_meshes;
        std::vector<SSMesh*> via_meshes;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) via_meshes.push_back(mesh);
            else                   sig_meshes.push_back(mesh);
        }

        struct TorusViaInfo { MWMath::Point3D position; double radius; };
        std::vector<TorusViaInfo> torus_via_infos;
        for (auto* mesh : mus->meshPtrs) {
            SSTorusMesh* torus = dynamic_cast<SSTorusMesh*>(mesh);
            if (torus) {
                torus_via_infos.push_back({ torus->PositionGlobal, torus->r });
            }
        }

        int num_sig        = (int)sig_meshes.size();
        int num_via        = (int)via_meshes.size();
        int num_torus_via  = (int)torus_via_infos.size();
        int num_etas_sig   = num_inner * num_sig;
        int num_etas_via   = num_via;
        int num_etas_torus = num_torus_via;
        int num_etas_total = num_etas_sig + num_etas_via + num_etas_torus;

        // ==========================================================
        // 1. INITIAL GUESS (x0)
        // ==========================================================
        for (size_t k = 1; k < mus->MNodes.size() - 1; ++k) {
            MWMath::Point3D pWarm = mus->MNodes[k].predictNewGlobal();
            x0_all.push_back(pWarm.x);
            x0_all.push_back(pWarm.y);
            x0_all.push_back(pWarm.z);
        }

        if (mus->lastEtas.size() == (size_t)num_etas_total && bUseWarmstartEtas) {
            for (double eta : mus->lastEtas)
                x0_all.push_back(eta * WarmstartEtaScaling);
        } else {
            for (int e = 0; e < num_etas_total; ++e) x0_all.push_back(0.0);
        }

        // ==========================================================
        // 2. PARAMETER (p)
        // ==========================================================
        // A. Signorini Meshes
        for (auto* mesh : sig_meshes) {
            p_all.push_back(mesh->PositionGlobal.x);
            p_all.push_back(mesh->PositionGlobal.y);
            p_all.push_back(mesh->PositionGlobal.z);
            MWMath::RotMatrix3x3 R = mesh->OrientationGlobal;
            for (int col = 0; col < 3; ++col)
                for (int row = 0; row < 3; ++row)
                    p_all.push_back(R.m[row][col]);
        }
        // B. Normale Via-Point Meshes
        for (auto* mesh : via_meshes) {
            p_all.push_back(mesh->PositionGlobal.x);
            p_all.push_back(mesh->PositionGlobal.y);
            p_all.push_back(mesh->PositionGlobal.z);
        }
        // C. NEU: Torus-Via-Positionen als Parameter
        for (auto& tvi : torus_via_infos) {
            p_all.push_back(tvi.position.x);
            p_all.push_back(tvi.position.y);
            p_all.push_back(tvi.position.z);
        }
        // D. Origin & Insertion
        p_all.push_back(mus->OriginPointGlobal.x);
        p_all.push_back(mus->OriginPointGlobal.y);
        p_all.push_back(mus->OriginPointGlobal.z);
        p_all.push_back(mus->InsertionPointGlobal.x);
        p_all.push_back(mus->InsertionPointGlobal.y);
        p_all.push_back(mus->InsertionPointGlobal.z);

        // ==========================================================
        // 3. CONSTRAINT BOUNDS (lbg/ubg)
        // ==========================================================
        // A. Normale Via-Points (h_v >= 0, sum == 0)
        for (int v = 0; v < num_via; ++v) {
            lbg_all.push_back(0.0); ubg_all.push_back(inf);
        }
        if (num_via > 0) {
            lbg_all.push_back(-inf); ubg_all.push_back(0.0);
        }

        // B. TORUS-VIA-POINTS DEAKTIVIEREN
        for (int t = 0; t < num_torus_via; ++t) {
            if (Step == 0) {
                // Schritt 0: Der Constraint gilt! (Abstand muss >= 0 sein)
                lbg_all.push_back(0.0); ubg_all.push_back(inf);
            } else {
                // Schritt > 0: Constraint komplett DEAKTIVIERT.
                // h_t darf jeden beliebigen Wert annehmen (von -inf bis +inf).
                // Der Solver muss dieses "Problem" nicht mehr lösen!
                lbg_all.push_back(-inf); ubg_all.push_back(inf);
            }
        }
        if (num_torus_via > 0) {
            if (Step == 0) {
                // Schritt 0: Summe(h*eta) <= 0 (Complementarity aktiv)
                lbg_all.push_back(-inf); ubg_all.push_back(0.0);
            } else {
                // Schritt > 0: Complementarity DEAKTIVIERT.
                lbg_all.push_back(-inf); ubg_all.push_back(inf);
            }
        }

        // C. Knoten-Schleife (Signorini + Euler-Lagrange)
        for (int k = 0; k < num_inner; ++k) {
            for (int s = 0; s < num_sig; ++s) {
                lbg_all.push_back(0.0); ubg_all.push_back(inf);
            }
            if (num_sig > 0) {
                lbg_all.push_back(0.0); ubg_all.push_back(0.0);
            }
            for (int d = 0; d < 3; ++d) {
                lbg_all.push_back(0.0); ubg_all.push_back(0.0);
            }
        }

        // ==========================================================
        // 4. VARIABLE BOUNDS (lbx/ubx)
        // ==========================================================
        // Knotenpositionen: unbegrenzt
        for (int i = 0; i < num_inner * 3; ++i) {
            lbx_all.push_back(-inf); ubx_all.push_back(inf);
        }
        // Signorini-Etas: >= 0
        for (int i = 0; i < num_etas_sig; ++i) {
            lbx_all.push_back(0.0); ubx_all.push_back(inf);
        }
        // Normale Via-Etas: >= 0
        for (int i = 0; i < num_etas_via; ++i) {
            lbx_all.push_back(0.0); ubx_all.push_back(inf);
        }
        // NEU: Torus-Etas (Variablen deaktivieren)
        for (int i = 0; i < num_etas_torus; ++i) {
            if (Step == 0) {
                // Schritt 0: Via-Point Eta darf frei arbeiten (>= 0)
                lbx_all.push_back(0.0); ubx_all.push_back(inf);
            } else {
                // Schritt > 0: Torus-Eta wird EXAKT AUF 0 GEZWUNGEN!
                // Wenn Eta = 0 ist, fällt die Via-Point Kraft (eta * grad_h) komplett weg.
                // Der Via-Point ist völlig unsichtbar für den Solver.
                lbx_all.push_back(0.0); ubx_all.push_back(0.0);
            }
        }
    }

    // ==========================================================
    // SOLVER AUFRUF
    // ==========================================================
    DMDict arg = {{"x0", x0_all}, {"p", p_all},
                  {"lbg", lbg_all}, {"ubg", ubg_all},
                  {"lbx", lbx_all}, {"ubx", ubx_all}};
    DMDict res = solver(arg);

    // ==========================================================
    // ERGEBNISSE EXTRAHIEREN
    // ==========================================================
    std::vector<double> res_x = std::vector<double>(res.at("x"));
    auto solverInfo = solver.stats();
    std::string status = solverInfo.at("return_status");

    const std::string C_RED   = "\033[31m";
    const std::string C_GREEN = "\033[32m";
    const std::string C_RESET = "\033[0m";
    int convSteps = int(solverInfo.at("iter_count"));
    {
        std::string stepColor   = (status == "Solve_Succeeded") ? C_GREEN : C_RED;
        std::string fullMessage = stepColor + "    Solver finished after "
                                + std::to_string(convSteps) + " iterations ("
                                + status + ")" + C_RESET;
        qDebug().noquote() << QString::fromStdString(fullMessage);
    }

    SolverConvergenceMessages.push_back(status);
    SolverConvergenceSteps.push_back(convSteps);

    // ==========================================================
    // ERGEBNISSE ZURÜCKSCHREIBEN
    // ==========================================================
    int current_x_offset = 0;
    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;

        std::vector<SSMesh*> sig_meshes, via_meshes;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) via_meshes.push_back(mesh);
            else                   sig_meshes.push_back(mesh);
        }

        struct TorusViaInfo { MWMath::Point3D position; double radius; };
        std::vector<TorusViaInfo> torus_via_infos;
        for (auto* mesh : mus->meshPtrs) {
            SSTorusMesh* torus = dynamic_cast<SSTorusMesh*>(mesh);
            if (torus) torus_via_infos.push_back({ torus->PositionGlobal, torus->r });
        }

        int num_sig        = (int)sig_meshes.size();
        int num_via        = (int)via_meshes.size();
        int num_torus_via  = (int)torus_via_infos.size();
        int num_etas_total = (num_inner * num_sig) + num_via + num_torus_via;

        int offset_nodes      = current_x_offset;
        int offset_sig_etas   = offset_nodes    + (num_inner * 3);
        int offset_via_etas   = offset_sig_etas + (num_inner * num_sig);
        int offset_torus_etas = offset_via_etas + num_via;

        // 1. Positionen zurückschreiben
        for (int i = 0; i < num_inner; ++i) {
            MWMath::Point3D p(
                res_x[offset_nodes + i * 3 + 0],
                res_x[offset_nodes + i * 3 + 1],
                res_x[offset_nodes + i * 3 + 2]);
            mus->MusclePointsGlobal[i + 1] = p;
            mus->MNodes[i + 1].PositionGlobal = p;
            mus->MNodes[i + 1].updateLocalFrame();
        }

        // 2. Alle Etas für Warmstart speichern
        mus->lastEtas.assign(
            res_x.begin() + offset_sig_etas,
            res_x.begin() + offset_sig_etas + num_etas_total);

        // 3. MAPPING FÜR EXPORT/VISUALISIERUNG
        // Via-Point → nächster Knoten
        std::vector<int> closest_node_for_via(num_via, -1);
        for (int v = 0; v < num_via; ++v) {
            MWMath::Point3D v_pos = via_meshes[v]->PositionGlobal;
            double min_dist = 1e9; int best_k = 0;
            for (int k = 0; k < num_inner; ++k) {
                double d = MWMath::distance(mus->MNodes[k + 1].PositionGlobal, v_pos);
                if (d < min_dist) { min_dist = d; best_k = k; }
            }
            closest_node_for_via[v] = best_k;
        }

        // Torus-Via → nächster Knoten (nur für Visualisierung relevant)
        std::vector<int> closest_node_for_torus(num_torus_via, -1);
        for (int t = 0; t < num_torus_via; ++t) {
            MWMath::Point3D t_pos = torus_via_infos[t].position;
            double min_dist = 1e9; int best_k = 0;
            for (int k = 0; k < num_inner; ++k) {
                double d = MWMath::distance(mus->MNodes[k + 1].PositionGlobal, t_pos);
                if (d < min_dist) { min_dist = d; best_k = k; }
            }
            closest_node_for_torus[t] = best_k;
        }

        // Etas pro Knoten füllen (für MNodeEtaSteps)
        for (int k = 0; k < num_inner; ++k) {
            // Größe: orig. meshPtrs + synthetische Torus-Via-Points
            std::vector<double> currentStepEtas(mus->meshPtrs.size() + num_torus_via, 0.0);

            int sig_idx = 0, via_idx = 0;
            for (size_t mi = 0; mi < mus->meshPtrs.size(); ++mi) {
                if (mus->meshPtrs[mi]->bIsViaPoint) {
                    if (closest_node_for_via[via_idx] == k)
                        currentStepEtas[mi] = res_x[offset_via_etas + via_idx];
                    via_idx++;
                } else {
                    currentStepEtas[mi] = res_x[offset_sig_etas + (k * num_sig) + sig_idx];
                    sig_idx++;
                }
            }
            // Torus-Etas anhängen
            for (int t = 0; t < num_torus_via; ++t) {
                if (closest_node_for_torus[t] == k)
                    currentStepEtas[mus->meshPtrs.size() + t] = res_x[offset_torus_etas + t];
            }

            mus->MNodes[k + 1].MNodeEtaSteps.push_back(currentStepEtas);
            mus->MNodes[k + 1].lastEtas = mus->lastEtas;
        }

        current_x_offset += (num_inner * 3) + num_etas_total;
    }

    // NEU: Schritt hochzählen — nach dem ersten Aufruf werden Torus-Etas deaktiviert
    Step++;
}

void CasadiSystem::setupCasadiToria()
{
    qDebug() << "     Setting up CasadiSystem with Via Point formulation...";
    using namespace casadi;
    MX all_x = MX::vertcat({});
    MX all_g = MX::vertcat({});
    MX all_p = MX::vertcat({});
    MX obj = 0;

    for (size_t m = 0; m < m_muscles.size(); ++m) {
        SSMuscle* mus = m_muscles[m];
        int num_inner = mus->MNodes.size() - 2;
        int K = mus->MNodes.size() - 1;

        // ==========================================================
        // 0. MESHES SORTIEREN + TORUS-FLAG MERKEN
        // ==========================================================
        std::vector<SSMesh*> sig_meshes;
        std::vector<SSMesh*> via_meshes;
        std::vector<bool>    via_is_torus; // NEU: merken welcher Via-Point ein Torus ist

        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) {
                via_meshes.push_back(mesh);
                // Torus-Erkennung per dynamic_cast
                via_is_torus.push_back(dynamic_cast<SSTorusMesh*>(mesh) != nullptr);
            } else {
                sig_meshes.push_back(mesh);
            }
        }
        
        // NEU: Automatisch Via-Points für jeden Torus anlegen (falls nicht schon vorhanden)
        // Diese werden als "synthetische" Via-Points am Torus-Zentrum eingefügt
        struct TorusViaInfo {
            MWMath::Point3D position;
            double          radius;
        };
        std::vector<TorusViaInfo> torus_via_infos;
        for (auto* mesh : mus->meshPtrs) {
            SSTorusMesh* torus = dynamic_cast<SSTorusMesh*>(mesh);
            if (torus) {
                // Torus-Zentrum als Via-Point-Position, Innenradius als Toleranz
                torus_via_infos.push_back({ torus->PositionGlobal, torus->r });
            }
        }
        int num_torus_via = (int)torus_via_infos.size(); // NEU

        int num_sig      = sig_meshes.size();
        int num_via      = via_meshes.size();
        int num_via_total = num_via + num_torus_via; // alle Via-Points inkl. Torus-synthetische

        int num_etas_sig      = num_inner * num_sig;
        int num_etas_via      = num_via;
        int num_etas_torus    = num_torus_via;                    // NEU
        int num_etas_total    = num_etas_sig + num_etas_via + num_etas_torus; // NEU

        // ==========================================================
        // 1. VARIABLEN DEFINIEREN (X-Vektor)
        // Aufbau: [ 3*N Knoten | N*M_sig Etas | M_via Etas | M_torus Etas ]
        // ==========================================================
        MX x_mus = MX::sym("x_" + std::to_string(m),
                           (num_inner * 3) + num_etas_sig + num_etas_via + num_etas_torus);
        all_x = MX::vertcat({all_x, x_mus});

        int offset_nodes     = 0;
        int offset_sig_etas  = num_inner * 3;
        int offset_via_etas  = offset_sig_etas + num_etas_sig;
        int offset_torus_etas = offset_via_etas + num_etas_via; // NEU

        // ==========================================================
        // 2. PARAMETER DEFINIEREN (P-Vektor)
        // Aufbau: [ 12*M_sig | 3*M_via | 3*M_torus | 3 Orig | 3 Ins ]
        // NEU: Torus-Via-Positionen als Parameter (damit Solver neu kompiliert werden kann)
        // ==========================================================
        MX p_mus = MX::sym("p_" + std::to_string(m),
                           (12 * num_sig) + (3 * num_via) + (3 * num_torus_via) + 6); // NEU
        all_p = MX::vertcat({all_p, p_mus});

        int p_offset_sig    = 0;
        int p_offset_via    = num_sig * 12;
        int p_offset_torus  = p_offset_via + (num_via * 3); // NEU
        int p_offset_ends   = p_offset_torus + (num_torus_via * 3); // NEU

        MX P_orig = p_mus(Slice(p_offset_ends,     p_offset_ends + 3));
        MX P_ins  = p_mus(Slice(p_offset_ends + 3, p_offset_ends + 6));

        // ==========================================================
        // 3a. NORMALE VIA-POINT BERECHNUNG (unverändert)
        // ==========================================================
        double eps = 1e-8;
        std::vector<MX> h_via_list;
        MX sum_h_eta_via = 0;

        for (int v = 0; v < num_via; ++v) {
            MX v_pos  = p_mus(Slice(p_offset_via + v * 3, p_offset_via + (v + 1) * 3));
            MX eta_v  = x_mus(offset_via_etas + v);

            double r_tol = via_meshes[v]->MViaPointTolerance * 0.5;
            if (r_tol <= 0.0001) r_tol = 0.005;
            double r_sq = r_tol * r_tol;

            MX d_vector = MX::vertcat({});
            for (int k = 0; k < num_inner; ++k) {
                MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
                MX dist = sqrt(sum1(sq(g_k - v_pos)) + eps);
                d_vector = MX::vertcat({d_vector, dist});
            }
            d_vector = MX::vertcat({d_vector, sqrt(sum1(sq(P_orig - v_pos)) + eps)});
            d_vector = MX::vertcat({d_vector, sqrt(sum1(sq(P_ins  - v_pos)) + eps)});

            MX D_v, h_v;
            if (bMaxLogSumTrick) {
                MX x    = -Alpha * d_vector;
                MX x_max = casadi::MX::mmax(x);
                MX sum_exp = sum1(exp(x - x_max));
                D_v = -(x_max + log(sum_exp)) / Alpha;
            } else {
                D_v = -casadi::MX::logsumexp(-Alpha * d_vector) / Alpha;
            }
            h_v = r_sq - D_v;

            h_via_list.push_back(h_v);
            sum_h_eta_via += h_v * eta_v;
            all_g = MX::vertcat({all_g, h_v}); // Distanz-Constraint
        }
        if (num_via > 0) {
            all_g = MX::vertcat({all_g, sum_h_eta_via});
        }

        // Gradienten für normale Via-Points
        std::vector<MX> grad_h_via_full_list;
        for (int v = 0; v < num_via; ++v) {
            grad_h_via_full_list.push_back(MX::gradient(h_via_list[v], x_mus));
        }

        // ==========================================================
        // 3b. TORUS VIA-POINT BERECHNUNG (NEU)
        // Schritt 0: echte h_v Constraints
        // Schritt > 0: h_v = r_sq (Konstante) → immer erfüllt, η=0 via ubx
        // ==========================================================
        std::vector<MX> h_torus_via_list;
        MX sum_h_eta_torus = 0;

        for (int t = 0; t < num_torus_via; ++t) {
            MX t_pos  = p_mus(Slice(p_offset_torus + t * 3, p_offset_torus + (t + 1) * 3));
            MX eta_t  = x_mus(offset_torus_etas + t);

            // Toleranz für den Via-Point (Wie beim normalen Via-Point)
            double r_tol = torus_via_infos[t].radius * 0.5;
            if (r_tol <= 0.0001) r_tol = 0.005;

            // Wir bauen IMMER die LogSumExp Funktion auf!
            MX d_vector = MX::vertcat({});
            for (int k = 0; k < num_inner; ++k) {
                MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
                MX dist = sqrt(sum1(sq(g_k - t_pos)) + eps);
                d_vector = MX::vertcat({d_vector, dist});
            }
            d_vector = MX::vertcat({d_vector, sqrt(sum1(sq(P_orig - t_pos)) + eps)});
            d_vector = MX::vertcat({d_vector, sqrt(sum1(sq(P_ins  - t_pos)) + eps)});

            MX D_t;
            if (bMaxLogSumTrick) {
                MX x     = -Alpha * d_vector;
                MX x_max = casadi::MX::mmax(x);
                MX sum_exp = sum1(exp(x - x_max));
                D_t = -(x_max + log(sum_exp)) / Alpha;
            } else {
                D_t = -casadi::MX::logsumexp(-Alpha * d_vector) / Alpha;
            }
            
            // h_t ist die Verletzung der Distanz (Radius - Approximierter Abstand)
            MX h_t = r_tol - D_t; 

            h_torus_via_list.push_back(h_t);
            sum_h_eta_torus += h_t * eta_t;
            all_g = MX::vertcat({all_g, h_t}); // Constraint hinzufügen
        }
        
        if (num_torus_via > 0) {
            all_g = MX::vertcat({all_g, sum_h_eta_torus});
        }

        // Gradienten für Torus-Via-Points berechnen
        std::vector<MX> grad_h_torus_full_list;
        for (int t = 0; t < num_torus_via; ++t) {
            grad_h_torus_full_list.push_back(MX::gradient(h_torus_via_list[t], x_mus));
        }

        // ==========================================================
        // 4. KNOTEN-SCHLEIFE (Signorini & Euler-Lagrange)
        // ==========================================================
        for (int k = 0; k < num_inner; ++k) {
            MX g_k    = x_mus(Slice(k * 3, (k + 1) * 3));
            MX g_prev = (k == 0)            ? P_orig : x_mus(Slice((k - 1) * 3, k * 3));
            MX g_next = (k == num_inner - 1) ? P_ins  : x_mus(Slice((k + 1) * 3, (k + 2) * 3));

            MX total_sig_force   = 0;
            MX total_via_force   = 0;
            MX total_torus_force = 0; // NEU
            MX sum_h_eta_sig     = 0;

            // --- A. SIGNORINI KRÄFTE ---
            for (int s = 0; s < num_sig; ++s) {
                MX q_s    = p_mus(Slice(p_offset_sig + s * 12, p_offset_sig + (s + 1) * 12));
                MX eta_ks = x_mus(offset_sig_etas + k * num_sig + s);

                MX h_s = sig_meshes[s]->constraintDistance(g_k, q_s);

                MX grad_h_s;
                if (bUseOwnGradient) grad_h_s = sig_meshes[s]->constraintJacobian(g_k, q_s);
                else {
                    MX grad_full = MX::gradient(h_s, x_mus);
                    grad_h_s = grad_full(Slice(k * 3, (k + 1) * 3));
                }

                total_sig_force += eta_ks * grad_h_s;
                sum_h_eta_sig   += h_s * eta_ks;
                all_g = MX::vertcat({all_g, h_s});
            }
            if (num_sig > 0) {
                all_g = MX::vertcat({all_g, sum_h_eta_sig});
            }

            // --- B. NORMALE VIA-POINT KRÄFTE ---
            for (int v = 0; v < num_via; ++v) {
                MX eta_v    = x_mus(offset_via_etas + v);
                MX grad_h_v = grad_h_via_full_list[v](Slice(k * 3, (k + 1) * 3));
                total_via_force += eta_v * grad_h_v;
            }

            // --- C. TORUS VIA-POINT KRÄFTE (NEU) ---
            // Schritt 0: echter Gradient → zieht Muskel durch Torus
            // Schritt > 0: Gradient von Konstante = 0 → kein Kraftbeitrag
            for (int t = 0; t < num_torus_via; ++t) {
                MX eta_t     = x_mus(offset_torus_etas + t);
                MX grad_h_t  = grad_h_torus_full_list[t](Slice(k * 3, (k + 1) * 3));
                total_torus_force += eta_t * grad_h_t;
            }

            // --- D. EULER-LAGRANGE KRAFTBILANZ ---
            MX total_contact_force = total_sig_force + total_via_force + total_torus_force;
            MX eq_el = (-g_prev + 2.0 * g_k - g_next) * K - (1.0 / K) * total_contact_force;
            all_g = MX::vertcat({all_g, eq_el});
        }
    }

    // ==========================================================
    // SOLVER CONFIG (unverändert)
    // ==========================================================
    opts["print_time"]                  = false;
    opts["record_time"]                 = false;
    opts["ipopt.print_level"]           = 0;
    opts["ipopt.max_iter"]              = maxIterations;
    opts["ipopt.tol"]                   = maxTol;
    opts["ipopt.linear_solver"]         = "mumps";
    opts["ipopt.hessian_approximation"] = "limited-memory";

    MXDict nlp = {{"x", all_x}, {"f", obj}, {"g", all_g}, {"p", all_p}};
    solver = nlpsol("f", "ipopt", nlp, opts);
}
 */


void CasadiSystem::setupCasadiToria()
{
    qDebug() << "     Setting up CasadiSystem Toria...";
    using namespace casadi;
    MX all_x = MX::vertcat({});
    MX all_g = MX::vertcat({});
    MX all_p = MX::vertcat({});
    MX obj   = 0;

    for (size_t m = 0; m < m_muscles.size(); ++m) {
        SSMuscle* mus = m_muscles[m];
        int num_inner = mus->MNodes.size() - 2;
        int K         = mus->MNodes.size() - 1;

        // ----------------------------------------------------------
        // 0. MESHES SORTIEREN
        //    Tori bleiben in sig_meshes (normales Signorini!)
        //    Zusätzlich sammeln wir ihre Positionen für Via-Points
        // ----------------------------------------------------------
        std::vector<SSMesh*> sig_meshes;
        std::vector<SSMesh*> via_meshes;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) via_meshes.push_back(mesh);
            else                   sig_meshes.push_back(mesh); // Tori landen hier!
        }

        // Tori extra merken für Via-Point-Zusatzterme
        struct TorusInfo { MWMath::Point3D pos; double r_via; };
        std::vector<TorusInfo> torus_infos;
        for (auto* mesh : mus->meshPtrs) {
            SSTorusMesh* torus = dynamic_cast<SSTorusMesh*>(mesh);
            if (torus) {
                double r_via = (torus->R - torus->r) * 1.1f ; //- torus->r ; 
                if (r_via < 0.005) r_via = 0.005;
                torus_infos.push_back({ torus->PositionGlobal, r_via });
            }
        }

        int num_sig       = (int)sig_meshes.size();
        int num_via       = (int)via_meshes.size();
        int num_torus_via = (int)torus_infos.size(); // extra Via-Point Etas für Tori

        // Eta-Layout:
        // [ N*num_sig (Signorini) | num_via (Via) | num_torus_via (Torus-Via, neu) ]
        int num_etas_sig      = num_inner * num_sig;
        int num_etas_via      = num_via;
        int num_etas_torus_via = num_torus_via;
        int num_etas_total    = num_etas_sig + num_etas_via + num_etas_torus_via;

        // ----------------------------------------------------------
        // 1. VARIABLEN
        // ----------------------------------------------------------
        MX x_mus = MX::sym("x_" + std::to_string(m),
                           (num_inner * 3) + num_etas_total);
        all_x = MX::vertcat({all_x, x_mus});

        int offset_nodes       = 0;
        int offset_sig_etas    = num_inner * 3;
        int offset_via_etas    = offset_sig_etas + num_etas_sig;
        int offset_tvia_etas   = offset_via_etas + num_etas_via;

        // ----------------------------------------------------------
        // 2. PARAMETER
        // [ 12*num_sig | 3*num_via | 3*num_torus_via | 3 Orig | 3 Ins ]
        // ----------------------------------------------------------
        MX p_mus = MX::sym("p_" + std::to_string(m),
                           (12 * num_sig) + (3 * num_via) + (3 * num_torus_via) + 6);
        all_p = MX::vertcat({all_p, p_mus});

        int p_off_sig   = 0;
        int p_off_via   = num_sig * 12;
        int p_off_tvia  = p_off_via  + num_via       * 3;
        int p_off_ends  = p_off_tvia + num_torus_via * 3;

        MX P_orig = p_mus(Slice(p_off_ends,     p_off_ends + 3));
        MX P_ins  = p_mus(Slice(p_off_ends + 3, p_off_ends + 6));

        constexpr double eps = 1e-8;

        // ----------------------------------------------------------
        // 3a. NORMALE VIA-POINTS (unverändert)
        // ----------------------------------------------------------
        std::vector<MX> h_via_list;
        MX sum_h_eta_via = 0;
        for (int v = 0; v < num_via; ++v) {
            MX v_pos = p_mus(Slice(p_off_via + v*3, p_off_via + (v+1)*3));
            MX eta_v = x_mus(offset_via_etas + v);

            double r_tol = via_meshes[v]->MViaPointTolerance * 1.5;
            if (r_tol <= 0.0001) r_tol = 0.005;
            double r_sq = r_tol * r_tol;

            MX d_vec = MX::vertcat({});
            for (int k = 0; k < num_inner; ++k) {
                MX g_k = x_mus(Slice(k*3, (k+1)*3));
                d_vec = MX::vertcat({d_vec, sqrt(sum1(sq(g_k - v_pos)) + eps)});
            }
            d_vec = MX::vertcat({d_vec, sqrt(sum1(sq(P_orig - v_pos)) + eps)});
            d_vec = MX::vertcat({d_vec, sqrt(sum1(sq(P_ins  - v_pos)) + eps)});

            MX D_v;
            if (bMaxLogSumTrick) {
                MX lx = -Alpha * d_vec;
                MX lx_max = casadi::MX::mmax(lx);
                D_v = -(lx_max + log(sum1(exp(lx - lx_max)))) / Alpha;
            } else {
                D_v = -casadi::MX::logsumexp(-Alpha * d_vec) / Alpha;
            }
            MX h_v = r_sq - D_v;

            h_via_list.push_back(h_v);
            sum_h_eta_via += h_v * eta_v;
            all_g = MX::vertcat({all_g, h_v});
        }
        if (num_via > 0) all_g = MX::vertcat({all_g, sum_h_eta_via});

        std::vector<MX> grad_h_via_list;
        for (int v = 0; v < num_via; ++v)
            grad_h_via_list.push_back(MX::gradient(h_via_list[v], x_mus));

        // ----------------------------------------------------------
        // 3b. TORUS VIA-POINT ZUSATZTERME (NUR FÜR STEP 0 aktiv)
        //     Zieht den Muskel beim ersten Schritt in den Torus rein.
        //     Ab Step 1: eta=0 (ubx=0) + lbg/ubg=[-inf,inf] → komplett inaktiv
        // ----------------------------------------------------------
        std::vector<MX> h_tvia_list;
        MX sum_h_eta_tvia = 0;

        for (int t = 0; t < num_torus_via; ++t) {
            MX t_pos  = p_mus(Slice(p_off_tvia + t*3, p_off_tvia + (t+1)*3));
            MX eta_tv = x_mus(offset_tvia_etas + t);

            double r_via = torus_infos[t].r_via;

            // Distanz jedes Knotens zum Torus-Zentrum via LogSumExp
            MX d_vec = MX::vertcat({});
            for (int k = 0; k < num_inner; ++k) {
                MX g_k = x_mus(Slice(k*3, (k+1)*3));
                d_vec = MX::vertcat({d_vec, sqrt(sum1(sq(g_k - t_pos)) + eps)});
            }
            d_vec = MX::vertcat({d_vec, sqrt(sum1(sq(P_orig - t_pos)) + eps)});
            d_vec = MX::vertcat({d_vec, sqrt(sum1(sq(P_ins  - t_pos)) + eps)});

            MX D_tv;
            if (bMaxLogSumTrick) {
                MX lx = -Alpha * d_vec;
                MX lx_max = casadi::MX::mmax(lx);
                D_tv = -(lx_max + log(sum1(exp(lx - lx_max)))) / Alpha;
            } else {
                D_tv = -casadi::MX::logsumexp(-Alpha * d_vec) / Alpha;
            }

            // h_tvia >= 0  ↔  mindestens ein Knoten ist näher als r_via am Torus-Zentrum
            MX h_tvia = r_via - D_tv;

            h_tvia_list.push_back(h_tvia);
            sum_h_eta_tvia += h_tvia * eta_tv;
            all_g = MX::vertcat({all_g, h_tvia}); // wird per lbg/ubg deaktiviert
        }
        if (num_torus_via > 0) all_g = MX::vertcat({all_g, sum_h_eta_tvia});

        std::vector<MX> grad_h_tvia_list;
        for (int t = 0; t < num_torus_via; ++t)
            grad_h_tvia_list.push_back(MX::gradient(h_tvia_list[t], x_mus));

        // ----------------------------------------------------------
        // 4. KNOTEN-SCHLEIFE: Signorini + Via + Torus-Via + EL
        // ----------------------------------------------------------
        for (int k = 0; k < num_inner; ++k) {
            MX g_k    = x_mus(Slice(k*3, (k+1)*3));
            MX g_prev = (k == 0)             ? P_orig : x_mus(Slice((k-1)*3, k*3));
            MX g_next = (k == num_inner - 1) ? P_ins  : x_mus(Slice((k+1)*3, (k+2)*3));

            MX F_sig  = 0, F_via = 0, F_tvia = 0;
            MX sum_h_eta_sig = 0;

            // A. Signorini (inkl. Tori — ganz normal!)
            for (int s = 0; s < num_sig; ++s) {
                MX q_s    = p_mus(Slice(p_off_sig + s*12, p_off_sig + (s+1)*12));
                MX eta_ks = x_mus(offset_sig_etas + k * num_sig + s);
                MX h_s    = sig_meshes[s]->constraintDistance(g_k, q_s);

                MX grad_h_s;
                if (bUseOwnGradient)
                    {grad_h_s = sig_meshes[s]->constraintJacobian(g_k, q_s);}
                else{
                    grad_h_s = MX::gradient(h_s, x_mus)(Slice(k*3, (k+1)*3));
                }

                F_sig         += eta_ks * grad_h_s;
                sum_h_eta_sig += h_s * eta_ks;
                all_g = MX::vertcat({all_g, h_s});
            }
            if (num_sig > 0) all_g = MX::vertcat({all_g, sum_h_eta_sig});

            // B. Normale Via-Points
            for (int v = 0; v < num_via; ++v) {
                MX eta_v = x_mus(offset_via_etas + v);
                F_via += eta_v * grad_h_via_list[v](Slice(k*3, (k+1)*3));
            }

            // C. Torus-Via-Point Kraft (Step 0: aktiv, Step>0: eta=0 → F=0)
            for (int t = 0; t < num_torus_via; ++t) {
                MX eta_tv = x_mus(offset_tvia_etas + t);
                F_tvia += eta_tv * grad_h_tvia_list[t](Slice(k*3, (k+1)*3));
            }

            // D. Euler-Lagrange
            MX eq_el = (-g_prev + 2.0*g_k - g_next) * K
                     - (1.0/K) * (F_sig + F_via + F_tvia);
            all_g = MX::vertcat({all_g, eq_el});
        }
    }

    opts["print_time"]                  = false;
    opts["record_time"]                 = false;
    opts["ipopt.print_level"]           = 0;
    opts["ipopt.max_iter"]              = maxIterations;
    opts["ipopt.tol"]                   = maxTol;
    opts["ipopt.linear_solver"]         = "mumps";
    opts["ipopt.hessian_approximation"] = "limited-memory";
    if (bWriteFiles) {
        std::string dirPath = "../examples/results";
        if (!std::filesystem::exists(dirPath)) {
            std::filesystem::create_directories(dirPath);
        }
        // 3. Datei-Pfad zusammensetzen
        std::string filename = dirPath + "/solver_log_" + CasadiSystemName + ".txt";
        opts["ipopt.output_file"] = filename; 
        opts["ipopt.file_print_level"] = 5;
    }

    MXDict nlp = {{"x", all_x}, {"f", obj}, {"g", all_g}, {"p", all_p}};
    solver = nlpsol("f", "ipopt", nlp, opts);
}

void CasadiSystem::solveStepToria()
{
    qDebug() << "          Solving Toria Step" << Step << "...";
    using namespace casadi;

    std::vector<double> x0_all, p_all, lbg_all, ubg_all, lbx_all, ubx_all;

    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;

        std::vector<SSMesh*> sig_meshes, via_meshes;
        // NEU: merken welche sig_meshes Tori sind
        std::vector<bool> sig_is_torus;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) {
                via_meshes.push_back(mesh);
            } else {
                sig_meshes.push_back(mesh);
                sig_is_torus.push_back(dynamic_cast<SSTorusMesh*>(mesh) != nullptr);
            }
        }

        struct TorusInfo { MWMath::Point3D pos; double r_via; };
        std::vector<TorusInfo> torus_infos;
        for (auto* mesh : mus->meshPtrs) {
            SSTorusMesh* torus = dynamic_cast<SSTorusMesh*>(mesh);
            if (torus) {
                double r_via = torus->R - torus->r;
                if (r_via < 0.005) r_via = 0.005;
                torus_infos.push_back({ torus->PositionGlobal, r_via });
            }
        }

        int num_sig            = (int)sig_meshes.size();
        int num_via            = (int)via_meshes.size();
        int num_torus_via      = (int)torus_infos.size();
        int num_etas_sig       = num_inner * num_sig;
        int num_etas_via       = num_via;
        int num_etas_torus_via = num_torus_via;
        int num_etas_total     = num_etas_sig + num_etas_via + num_etas_torus_via;

        // 1. INITIAL GUESS
        for (size_t k = 1; k < mus->MNodes.size() - 1; ++k) {
            MWMath::Point3D pw = mus->MNodes[k].predictNewGlobal();
            x0_all.push_back(pw.x);
            x0_all.push_back(pw.y);
            x0_all.push_back(pw.z);
        }
        if (mus->lastEtas.size() == (size_t)num_etas_total && bUseWarmstartEtas)
            for (double e : mus->lastEtas) x0_all.push_back(e * WarmstartEtaScaling);
        else
            for (int e = 0; e < num_etas_total; ++e) x0_all.push_back(0.0);

        // 2. PARAMETER (unverändert)
        for (auto* mesh : sig_meshes) {
            p_all.push_back(mesh->PositionGlobal.x);
            p_all.push_back(mesh->PositionGlobal.y);
            p_all.push_back(mesh->PositionGlobal.z);
            MWMath::RotMatrix3x3 R = mesh->OrientationGlobal;
            for (int c = 0; c < 3; ++c)
                for (int r = 0; r < 3; ++r)
                    p_all.push_back(R.m[r][c]);
        }
        for (auto* mesh : via_meshes) {
            p_all.push_back(mesh->PositionGlobal.x);
            p_all.push_back(mesh->PositionGlobal.y);
            p_all.push_back(mesh->PositionGlobal.z);
        }
        for (auto& ti : torus_infos) {
            p_all.push_back(ti.pos.x);
            p_all.push_back(ti.pos.y);
            p_all.push_back(ti.pos.z);
        }
        p_all.push_back(mus->OriginPointGlobal.x);
        p_all.push_back(mus->OriginPointGlobal.y);
        p_all.push_back(mus->OriginPointGlobal.z);
        p_all.push_back(mus->InsertionPointGlobal.x);
        p_all.push_back(mus->InsertionPointGlobal.y);
        p_all.push_back(mus->InsertionPointGlobal.z);

        // 3. CONSTRAINT BOUNDS
        // A. Normale Via-Points (immer aktiv)
        for (int v = 0; v < num_via; ++v) {
            lbg_all.push_back(0.0); ubg_all.push_back(inf);
        }
        if (num_via > 0) { lbg_all.push_back(-inf); ubg_all.push_back(0.0); }

        // B. Torus-Via-Point Constraints
        for (int t = 0; t < num_torus_via; ++t) {
            if (Step == 0) { lbg_all.push_back(0.0);  ubg_all.push_back(inf); }
            else           { lbg_all.push_back(-inf); ubg_all.push_back(inf); }
        }
        if (num_torus_via > 0) {
            if (Step == 0) { lbg_all.push_back(-inf); ubg_all.push_back(0.0); }
            else           { lbg_all.push_back(-inf); ubg_all.push_back(inf); }
        }

        // C. Signorini + EL — pro Knoten
        for (int k = 0; k < num_inner; ++k) {

            // h_s >= 0 pro Signorini-Mesh
            for (int s = 0; s < num_sig; ++s) {
                /* if (Step == 0 && sig_is_torus[s]) {
                    // Torus-Signorini bei Step 0 deaktivieren: h_s darf alles sein
                    lbg_all.push_back(-inf); ubg_all.push_back(inf);
                } else {
                    lbg_all.push_back(0.0); ubg_all.push_back(inf);
                } */
                lbg_all.push_back(0.0); ubg_all.push_back(inf);  // obstacle always active
            }

            // sum(h_s * eta_s) == 0
            if (num_sig > 0) {
                // Bei Step 0 mit Tori: die Summe enthält Torus-Terme mit eta=0,
                // also bleibt die Gleichung trotzdem konsistent → immer [0,0]
                lbg_all.push_back(0.0); ubg_all.push_back(0.0);
            }

            // Euler-Lagrange (immer Gleichung)
            for (int d = 0; d < 3; ++d) {
                lbg_all.push_back(0.0); ubg_all.push_back(0.0);
            }
        }

        // 4. VARIABLE BOUNDS
        // Knotenpositionen
        for (int i = 0; i < num_inner * 3; ++i) {
            lbx_all.push_back(-inf); ubx_all.push_back(inf);
        }

        // Signorini-Etas: Torus-Etas bei Step 0 auf 0 fixieren
        for (int k = 0; k < num_inner; ++k) {
            for (int s = 0; s < num_sig; ++s) {
                /* if (Step == 0 && sig_is_torus[s]) {
                    // Torus-Signorini-Eta = 0 → kein Kraftbeitrag, Constraint inaktiv
                    lbx_all.push_back(0.0); ubx_all.push_back(0.0);
                } else {
                    lbx_all.push_back(0.0); ubx_all.push_back(inf);
                } */
                lbx_all.push_back(0.0); ubx_all.push_back(inf); // obstacle always active
            }
        }

        // Via-Etas (immer aktiv)
        for (int i = 0; i < num_etas_via; ++i) {
            lbx_all.push_back(0.0); ubx_all.push_back(inf);
        }

        // Torus-Via-Etas
        for (int i = 0; i < num_etas_torus_via; ++i) {
            if (Step == 0) { lbx_all.push_back(0.0); ubx_all.push_back(inf); }
            else           { lbx_all.push_back(0.0); ubx_all.push_back(0.0); }
        }
    }

    // SOLVER
    DMDict arg = {{"x0", x0_all}, {"p", p_all},
                  {"lbg", lbg_all}, {"ubg", ubg_all},
                  {"lbx", lbx_all}, {"ubx", ubx_all}};
    DMDict res = solver(arg);

    std::vector<double> res_x = std::vector<double>(res.at("x"));
    auto solverInfo = solver.stats();
    std::string status    = solverInfo.at("return_status");
    int convSteps         = int(solverInfo.at("iter_count"));

    {
        const std::string C_RED = "\033[31m", C_GREEN = "\033[32m", C_RESET = "\033[0m";
        qDebug().noquote() << QString::fromStdString(
            ((status == "Solve_Succeeded") ? C_GREEN : C_RED)
            + "    Solver: " + std::to_string(convSteps)
            + " iters (" + status + ")" + C_RESET);
    }
    SolverConvergenceMessages.push_back(status);
    SolverConvergenceSteps.push_back(convSteps);

    // ERGEBNISSE ZURÜCKSCHREIBEN
    int cur_x = 0;
    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;

        std::vector<SSMesh*> sig_meshes, via_meshes;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) via_meshes.push_back(mesh);
            else                   sig_meshes.push_back(mesh);
        }

        struct TorusInfo { MWMath::Point3D pos; double r_via; };
        std::vector<TorusInfo> torus_infos;
        for (auto* mesh : mus->meshPtrs) {
            SSTorusMesh* torus = dynamic_cast<SSTorusMesh*>(mesh);
            if (torus) {
                double rv = torus->r * 0.5;
                if (rv < 0.005) rv = 0.005;
                torus_infos.push_back({ torus->PositionGlobal, rv });
            }
        }

        int num_sig        = (int)sig_meshes.size();
        int num_via        = (int)via_meshes.size();
        int num_torus_via  = (int)torus_infos.size();
        int num_etas_total = (num_inner * num_sig) + num_via + num_torus_via;

        int off_nodes  = cur_x;
        int off_sig    = off_nodes + num_inner * 3;
        int off_via    = off_sig   + num_inner * num_sig;
        int off_tvia   = off_via   + num_via;

        for (int i = 0; i < num_inner; ++i) {
            MWMath::Point3D p(res_x[off_nodes + i*3],
                              res_x[off_nodes + i*3 + 1],
                              res_x[off_nodes + i*3 + 2]);
            mus->MusclePointsGlobal[i+1]    = p;
            mus->MNodes[i+1].PositionGlobal = p;
            mus->MNodes[i+1].updateLocalFrame();
        }

        mus->lastEtas.assign(res_x.begin() + off_sig,
                             res_x.begin() + off_sig + num_etas_total);

        // Mapping nächster Knoten für Visualisierung
        std::vector<int> cv(num_via, -1);
        for (int v = 0; v < num_via; ++v) {
            MWMath::Point3D vp = via_meshes[v]->PositionGlobal;
            double md = 1e9; int bk = 0;
            for (int k = 0; k < num_inner; ++k) {
                double d = MWMath::distance(mus->MNodes[k+1].PositionGlobal, vp);
                if (d < md) { md = d; bk = k; }
            }
            cv[v] = bk;
        }
        std::vector<int> ct(num_torus_via, -1);
        for (int t = 0; t < num_torus_via; ++t) {
            MWMath::Point3D tp = torus_infos[t].pos;
            double md = 1e9; int bk = 0;
            for (int k = 0; k < num_inner; ++k) {
                double d = MWMath::distance(mus->MNodes[k+1].PositionGlobal, tp);
                if (d < md) { md = d; bk = k; }
            }
            ct[t] = bk;
        }

        for (int k = 0; k < num_inner; ++k) {
            std::vector<double> stepEtas(mus->meshPtrs.size() + num_torus_via, 0.0);
            int si = 0, vi = 0;
            for (size_t mi = 0; mi < mus->meshPtrs.size(); ++mi) {
                if (mus->meshPtrs[mi]->bIsViaPoint) {
                    if (cv[vi] == k) stepEtas[mi] = res_x[off_via + vi];
                    vi++;
                } else {
                    stepEtas[mi] = res_x[off_sig + k * num_sig + si++];
                }
            }
            for (int t = 0; t < num_torus_via; ++t)
                if (ct[t] == k) stepEtas[mus->meshPtrs.size() + t] = res_x[off_tvia + t];

            mus->MNodes[k+1].MNodeEtaSteps.push_back(stepEtas);
            mus->MNodes[k+1].lastEtas = mus->lastEtas;
        }

        cur_x += (num_inner * 3) + num_etas_total;
    }

    Step++;
}

/* 
void CasadiSystem::solveStepToria()
{
    qDebug() << "          Solving Toria Step" << Step << "...";
    using namespace casadi;

    std::vector<double> x0_all, p_all, lbg_all, ubg_all, lbx_all, ubx_all;

    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;

        std::vector<SSMesh*> sig_meshes, via_meshes;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) via_meshes.push_back(mesh);
            else                   sig_meshes.push_back(mesh);
        }

        struct TorusInfo { MWMath::Point3D pos; double r_via; };
        std::vector<TorusInfo> torus_infos;
        for (auto* mesh : mus->meshPtrs) {
            SSTorusMesh* torus = dynamic_cast<SSTorusMesh*>(mesh);
            if (torus) {
                double r_via = torus->r * 0.5;
                if (r_via < 0.005) r_via = 0.005;
                torus_infos.push_back({ torus->PositionGlobal, r_via });
            }
        }

        int num_sig           = (int)sig_meshes.size();
        int num_via           = (int)via_meshes.size();
        int num_torus_via     = (int)torus_infos.size();
        int num_etas_sig      = num_inner * num_sig;
        int num_etas_via      = num_via;
        int num_etas_torus_via = num_torus_via;
        int num_etas_total    = num_etas_sig + num_etas_via + num_etas_torus_via;

        // 1. INITIAL GUESS
        for (size_t k = 1; k < mus->MNodes.size() - 1; ++k) {
            MWMath::Point3D pw = mus->MNodes[k].predictNewGlobal();
            x0_all.push_back(pw.x);
            x0_all.push_back(pw.y);
            x0_all.push_back(pw.z);
        }
        if (mus->lastEtas.size() == (size_t)num_etas_total && bUseWarmstartEtas)
            for (double e : mus->lastEtas) x0_all.push_back(e * WarmstartEtaScaling);
        else
            for (int e = 0; e < num_etas_total; ++e) x0_all.push_back(0.0);

        // 2. PARAMETER
        for (auto* mesh : sig_meshes) {
            p_all.push_back(mesh->PositionGlobal.x);
            p_all.push_back(mesh->PositionGlobal.y);
            p_all.push_back(mesh->PositionGlobal.z);
            MWMath::RotMatrix3x3 R = mesh->OrientationGlobal;
            for (int c = 0; c < 3; ++c)
                for (int r = 0; r < 3; ++r)
                    p_all.push_back(R.m[r][c]);
        }
        for (auto* mesh : via_meshes) {
            p_all.push_back(mesh->PositionGlobal.x);
            p_all.push_back(mesh->PositionGlobal.y);
            p_all.push_back(mesh->PositionGlobal.z);
        }
        // Torus-Zentren (immer übergeben, Graph erwartet sie)
        for (auto& ti : torus_infos) {
            p_all.push_back(ti.pos.x);
            p_all.push_back(ti.pos.y);
            p_all.push_back(ti.pos.z);
        }
        p_all.push_back(mus->OriginPointGlobal.x);
        p_all.push_back(mus->OriginPointGlobal.y);
        p_all.push_back(mus->OriginPointGlobal.z);
        p_all.push_back(mus->InsertionPointGlobal.x);
        p_all.push_back(mus->InsertionPointGlobal.y);
        p_all.push_back(mus->InsertionPointGlobal.z);

        // 3. CONSTRAINT BOUNDS
        // A. Normale Via-Points (immer aktiv)
        for (int v = 0; v < num_via; ++v) {
            lbg_all.push_back(0.0); ubg_all.push_back(inf);
        }
        if (num_via > 0) { lbg_all.push_back(-inf); ubg_all.push_back(0.0); }

        // B. Torus-Via-Point Constraints
        for (int t = 0; t < num_torus_via; ++t) {
            if (Step == 0) { lbg_all.push_back(0.0);  ubg_all.push_back(inf); }
            else           { lbg_all.push_back(-inf); ubg_all.push_back(inf); } // deaktiviert
        }
        if (num_torus_via > 0) {
            if (Step == 0) { lbg_all.push_back(-inf); ubg_all.push_back(0.0); }
            else           { lbg_all.push_back(-inf); ubg_all.push_back(inf); } // deaktiviert
        }

        // C. Signorini + EL (immer aktiv)
        for (int k = 0; k < num_inner; ++k) {
            for (int s = 0; s < num_sig; ++s) {
                lbg_all.push_back(0.0); ubg_all.push_back(inf);
            }
            if (num_sig > 0) { lbg_all.push_back(0.0); ubg_all.push_back(0.0); }
            for (int d = 0; d < 3; ++d) {
                lbg_all.push_back(0.0); ubg_all.push_back(0.0);
            }
        }

        // 4. VARIABLE BOUNDS
        for (int i = 0; i < num_inner * 3; ++i) {
            lbx_all.push_back(-inf); ubx_all.push_back(inf);
        }
        for (int i = 0; i < num_etas_sig; ++i) {
            lbx_all.push_back(0.0); ubx_all.push_back(inf);
        }
        for (int i = 0; i < num_etas_via; ++i) {
            lbx_all.push_back(0.0); ubx_all.push_back(inf);
        }
        for (int i = 0; i < num_etas_torus_via; ++i) {
            if (Step == 0) { lbx_all.push_back(0.0); ubx_all.push_back(inf); }  // aktiv
            else           { lbx_all.push_back(0.0); ubx_all.push_back(0.0); }  // eta=0
        }
    }

    // SOLVER
    DMDict arg = {{"x0", x0_all}, {"p", p_all},
                  {"lbg", lbg_all}, {"ubg", ubg_all},
                  {"lbx", lbx_all}, {"ubx", ubx_all}};
    DMDict res = solver(arg);

    std::vector<double> res_x = std::vector<double>(res.at("x"));
    auto solverInfo = solver.stats();
    std::string status    = solverInfo.at("return_status");
    int convSteps         = int(solverInfo.at("iter_count"));

    {
        const std::string C_RED = "\033[31m", C_GREEN = "\033[32m", C_RESET = "\033[0m";
        qDebug().noquote() << QString::fromStdString(
            ((status == "Solve_Succeeded") ? C_GREEN : C_RED)
            + "    Solver: " + std::to_string(convSteps)
            + " iters (" + status + ")" + C_RESET);
    }
    SolverConvergenceMessages.push_back(status);
    SolverConvergenceSteps.push_back(convSteps);

    // ERGEBNISSE ZURÜCKSCHREIBEN
    int cur_x = 0;
    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;

        std::vector<SSMesh*> sig_meshes, via_meshes;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) via_meshes.push_back(mesh);
            else                   sig_meshes.push_back(mesh);
        }

        struct TorusInfo { MWMath::Point3D pos; double r_via; };
        std::vector<TorusInfo> torus_infos;
        for (auto* mesh : mus->meshPtrs) {
            SSTorusMesh* torus = dynamic_cast<SSTorusMesh*>(mesh);
            if (torus) {
                double rv = torus->r * 0.5;
                if (rv < 0.005) rv = 0.005;
                torus_infos.push_back({ torus->PositionGlobal, rv });
            }
        }

        int num_sig        = (int)sig_meshes.size();
        int num_via        = (int)via_meshes.size();
        int num_torus_via  = (int)torus_infos.size();
        int num_etas_total = (num_inner * num_sig) + num_via + num_torus_via;

        int off_nodes  = cur_x;
        int off_sig    = off_nodes + num_inner * 3;
        int off_via    = off_sig   + num_inner * num_sig;
        int off_tvia   = off_via   + num_via;

        for (int i = 0; i < num_inner; ++i) {
            MWMath::Point3D p(res_x[off_nodes + i*3],
                              res_x[off_nodes + i*3 + 1],
                              res_x[off_nodes + i*3 + 2]);
            mus->MusclePointsGlobal[i+1]    = p;
            mus->MNodes[i+1].PositionGlobal = p;
            mus->MNodes[i+1].updateLocalFrame();
        }

        mus->lastEtas.assign(res_x.begin() + off_sig,
                             res_x.begin() + off_sig + num_etas_total);

        // Mapping nächster Knoten für Visualisierung
        std::vector<int> cv(num_via, -1);
        for (int v = 0; v < num_via; ++v) {
            MWMath::Point3D vp = via_meshes[v]->PositionGlobal;
            double md = 1e9; int bk = 0;
            for (int k = 0; k < num_inner; ++k) {
                double d = MWMath::distance(mus->MNodes[k+1].PositionGlobal, vp);
                if (d < md) { md = d; bk = k; }
            }
            cv[v] = bk;
        }
        std::vector<int> ct(num_torus_via, -1);
        for (int t = 0; t < num_torus_via; ++t) {
            MWMath::Point3D tp = torus_infos[t].pos;
            double md = 1e9; int bk = 0;
            for (int k = 0; k < num_inner; ++k) {
                double d = MWMath::distance(mus->MNodes[k+1].PositionGlobal, tp);
                if (d < md) { md = d; bk = k; }
            }
            ct[t] = bk;
        }

        for (int k = 0; k < num_inner; ++k) {
            std::vector<double> stepEtas(mus->meshPtrs.size() + num_torus_via, 0.0);
            int si = 0, vi = 0;
            for (size_t mi = 0; mi < mus->meshPtrs.size(); ++mi) {
                if (mus->meshPtrs[mi]->bIsViaPoint) {
                    if (cv[vi] == k) stepEtas[mi] = res_x[off_via + vi];
                    vi++;
                } else {
                    stepEtas[mi] = res_x[off_sig + k * num_sig + si++];
                }
            }
            for (int t = 0; t < num_torus_via; ++t)
                if (ct[t] == k) stepEtas[mus->meshPtrs.size() + t] = res_x[off_tvia + t];

            mus->MNodes[k+1].MNodeEtaSteps.push_back(stepEtas);
            mus->MNodes[k+1].lastEtas = mus->lastEtas;
        }

        cur_x += (num_inner * 3) + num_etas_total;
    }

    Step++;
}
 */


// VIA POINTS SUM
void CasadiSystem::solveStepViaSum_paramStudy()
{
    //qDebug() << "          Solving step with Via Point Sum(paramStudy) formulation...";
    using namespace casadi;

    std::vector<double> x0_all, p_all, lbg_all, ubg_all, lbx_all, ubx_all;

    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;

        // ==========================================================
        // 0. MESHES SORTIEREN
        // ==========================================================
        std::vector<SSMesh*> sig_meshes;
        std::vector<SSMesh*> via_meshes;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) via_meshes.push_back(mesh);
            else sig_meshes.push_back(mesh);
        }
        
        int num_sig = sig_meshes.size();
        int num_via = via_meshes.size();
        int num_etas_sig = num_inner * num_sig;
        int num_etas_via = num_via;
        int num_etas_total = num_etas_sig + num_etas_via;

        // ==========================================================
        // 1. INITIAL GUESS (x0)
        // ==========================================================
        for (size_t k = 1; k < mus->MNodes.size() - 1; ++k) {
            MWMath::Point3D pWarm = mus->MNodes[k].predictNewGlobal();
            x0_all.push_back(pWarm.x);
            x0_all.push_back(pWarm.y);
            x0_all.push_back(pWarm.z);
        }
        
        // Warmstart für BEIDE Eta-Arten
        if (mus->lastEtas.size() == (size_t)num_etas_total && bUseWarmstartEtas) {
            for (double eta : mus->lastEtas) {
                x0_all.push_back(eta * WarmstartEtaScaling);
            }
        } else {
            for (int e = 0; e < num_etas_total; ++e) x0_all.push_back(0.0);
        }

        // ==========================================================
        // 2. PARAMETER (p)
        // ==========================================================
        // A. Signorini Meshes (12 Parameter: Pos + Rot Matrix)
        for (auto* mesh : sig_meshes) {
            p_all.push_back(mesh->PositionGlobal.x);
            p_all.push_back(mesh->PositionGlobal.y);
            p_all.push_back(mesh->PositionGlobal.z);
            
            MWMath::RotMatrix3x3 R = mesh->OrientationGlobal;
            for (int col = 0; col < 3; ++col) {
                for (int row = 0; row < 3; ++row) {
                    p_all.push_back(R.m[row][col]);
                }
            }
        }
        
        // B. Via-Point Meshes (Nur 3 Parameter: Position)
        for (auto* mesh : via_meshes) {
            p_all.push_back(mesh->PositionGlobal.x);
            p_all.push_back(mesh->PositionGlobal.y);
            p_all.push_back(mesh->PositionGlobal.z);
        }
        
        // C. Origin & Insertion
        p_all.push_back(mus->OriginPointGlobal.x);
        p_all.push_back(mus->OriginPointGlobal.y);
        p_all.push_back(mus->OriginPointGlobal.z);
        
        p_all.push_back(mus->InsertionPointGlobal.x);
        p_all.push_back(mus->InsertionPointGlobal.y);
        p_all.push_back(mus->InsertionPointGlobal.z);

        // ==========================================================
        // 3. CONSTRAINT BOUNDS (lbg/ubg)
        // ==========================================================
        // A. Globale Via-Points
        for (int v = 0; v < num_via; ++v) {
            // NUR NOCH Distanz-Constraint: h_v >= 0
            lbg_all.push_back(0.0); 
            ubg_all.push_back(inf);
        }
        // Summierte Via-Point-Komplementarität (sum(h_v * eta_v) <= oder == 0)
        if (num_via > 0) {
            lbg_all.push_back(-inf);
            ubg_all.push_back(0.0);
        }


        // B. Knoten-Schleife (Signorini + Euler-Lagrange)
        for (int k = 0; k < num_inner; ++k) {
            // Non-penetration für jedes Signorini-Mesh (h_s >= 0)
            for (int s = 0; s < num_sig; ++s) {
                lbg_all.push_back(0.0); 
                ubg_all.push_back(inf);
            }
            // Summierte Signorini-Komplementarität (sum(h*eta) == 0)
            if (num_sig > 0) {
                lbg_all.push_back(0.0); 
                ubg_all.push_back(0.0);
            }
            // Euler-Lagrange Kraftbilanz (x, y, z == 0)
            for (int d = 0; d < 3; ++d) { 
                lbg_all.push_back(0.0); 
                ubg_all.push_back(0.0); 
            }
        }

        // ==========================================================
        // 4. VARIABLE BOUNDS (lbx/ubx)
        // ==========================================================
        // Gamma (Punkte)
        for (int i = 0; i < num_inner * 3; ++i) {
            lbx_all.push_back(-inf); ubx_all.push_back(inf);
        }
        // Alle Etas (Signorini + Via-Points) dürfen nur ziehen/drücken (>= 0)
        for (int i = 0; i < num_etas_total; ++i) {
            lbx_all.push_back(0.0); ubx_all.push_back(inf);
        }
    }

    // ==========================================================
    // SOLVER AUFRUF
    // ==========================================================
    DMDict arg = {{"x0", x0_all}, {"p", p_all}, {"lbg", lbg_all}, {"ubg", ubg_all}, {"lbx", lbx_all}, {"ubx", ubx_all}};
    DMDict res = solver(arg);

    // ==========================================================
    // ERGEBNISSE EXTRAHIEREN
    // ==========================================================
    std::vector<double> res_x = std::vector<double>(res.at("x"));
    auto solverInfo = solver.stats();
    std::string status = solverInfo.at("return_status");

    const std::string C_RED   = "\033[31m";
    const std::string C_GREEN = "\033[32m";
    const std::string C_RESET = "\033[0m";
    int convSteps = int(solverInfo.at("iter_count"));
    if (bDebug) { // Dein if(true) oder if(bDebug)
        std::string stepColor = (status == "Solve_Succeeded") ? C_GREEN : C_RED;
        std::string fullMessage = stepColor + "    Solver finished after " + std::to_string(convSteps) + " iterations" + "(" + status + ")" + C_RESET;
        qDebug().noquote() << QString::fromStdString(fullMessage);
    }
    
    SolverConvergenceMessages.push_back(status);
    SolverConvergenceSteps.push_back(convSteps);
    
    int current_x_offset = 0;
    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;
        int num_sig = 0, num_via = 0;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) num_via++; else num_sig++;
        }
        
        int offset_nodes    = current_x_offset;
        int offset_sig_etas = offset_nodes + (num_inner * 3);
        int offset_via_etas = offset_sig_etas + (num_inner * num_sig);

        // 1. Positionen zurückschreiben
        for (int i = 0; i < num_inner; ++i) {
            int mIdx = i + 1;
            MWMath::Point3D p(
                res_x[offset_nodes + i * 3 + 0],
                res_x[offset_nodes + i * 3 + 1],
                res_x[offset_nodes + i * 3 + 2]
            );
            mus->MusclePointsGlobal[mIdx] = p;
            mus->MNodes[mIdx].PositionGlobal = p;
            mus->MNodes[mIdx].updateLocalFrame();
        }

        // 2. Globale Etas für Warmstart speichern
        int num_etas_total = (num_inner * num_sig) + num_via;
        mus->lastEtas.assign(res_x.begin() + offset_sig_etas, res_x.begin() + offset_sig_etas + num_etas_total);

        // ==========================================================
        // 3. MAPPING-TRICK FÜR EXPORT/VISUALISIERUNG
        // Wir dröseln die Ergebnisse für die alte Struktur wieder auf!
        // ==========================================================
        
        // Zuerst: Welcher Knoten ist jedem Via-Point am nächsten?
        std::vector<int> closest_node_for_via(num_via, -1);
        int via_counter = 0;
        for (auto* mesh : mus->meshPtrs) {
            if (!mesh->bIsViaPoint) continue;
            
            MWMath::Point3D v_pos = mesh->PositionGlobal;
            double min_dist = 1e9;
            int best_k = 0;
            
            for (int k = 0; k < num_inner; ++k) {
                double d = MWMath::distance(mus->MNodes[k + 1].PositionGlobal, v_pos);
                if (d < min_dist) { min_dist = d; best_k = k; }
            }
            closest_node_for_via[via_counter] = best_k;
            via_counter++;
        }

        // Nun füllen wir das Array für JEDEN Knoten
        for (int k = 0; k < num_inner; ++k) {
            int mIdx = k + 1;
            std::vector<double> currentStepEtas(mus->meshPtrs.size(), 0.0);

            int sig_idx = 0;
            int via_idx = 0;

            for (size_t m = 0; m < mus->meshPtrs.size(); ++m) {
                if (mus->meshPtrs[m]->bIsViaPoint) {
                    // Via-Point Kraft nur eintragen, wenn das der nächste Knoten ist
                    if (closest_node_for_via[via_idx] == k) {
                        currentStepEtas[m] = res_x[offset_via_etas + via_idx];
                    }
                    via_idx++;
                } else {
                    // Signorini Kraft ganz regulär auslesen
                    currentStepEtas[m] = res_x[offset_sig_etas + (k * num_sig) + sig_idx];
                    sig_idx++;
                }
            }
            
            mus->MNodes[mIdx].MNodeEtaSteps.push_back(currentStepEtas);
            mus->MNodes[mIdx].lastEtas = mus->lastEtas; 
        }

        current_x_offset += (num_inner * 3) + num_etas_total;
    }
}

void CasadiSystem::setupCasadiViaSum_paramStudy()
{
    qDebug() << "     Setting up CasadiSystem with Via Point(paramStudy) formulation...";
    using namespace casadi;
    MX all_x = MX::vertcat({});
    MX all_g = MX::vertcat({});
    MX all_p = MX::vertcat({});
    MX obj = 0;
    
    

    for (size_t m = 0; m < m_muscles.size(); ++m) {
        SSMuscle* mus = m_muscles[m];
        int num_inner = mus->MNodes.size() - 2;
        int K = mus->MNodes.size() - 1;

        // ==========================================================
        // 0. MESHES SORTIEREN
        // ==========================================================
        std::vector<SSMesh*> sig_meshes;
        std::vector<SSMesh*> via_meshes;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) via_meshes.push_back(mesh);
            else sig_meshes.push_back(mesh);
        }
        int num_sig = sig_meshes.size();
        int num_via = via_meshes.size();

        // ==========================================================
        // 1. VARIABLEN DEFINIEREN (X-Vektor)
        // Aufbau: [ 3*N Knoten | N * M_sig Etas | M_via Etas ]
        // ==========================================================
        int num_etas_sig = num_inner * num_sig;
        int num_etas_via = num_via;
        
        MX x_mus = MX::sym("x_" + std::to_string(m), (num_inner * 3) + num_etas_sig + num_etas_via);
        all_x = MX::vertcat({all_x, x_mus});

        // Offsets zum Herausschneiden aus x_mus
        int offset_nodes    = 0;
        int offset_sig_etas = num_inner * 3;
        int offset_via_etas = offset_sig_etas + num_etas_sig;

        // ==========================================================
        // 2. PARAMETER DEFINIEREN (P-Vektor)
        // Aufbau: [ 12*M_sig Parameter | 3*M_via Parameter | 3 Orig | 3 Ins ]
        // ==========================================================
        MX p_mus = MX::sym("p_" + std::to_string(m), (12 * num_sig) + (3 * num_via) + 6);
        all_p = MX::vertcat({all_p, p_mus});

        // Offsets zum Herausschneiden aus p_mus
        int p_offset_sig  = 0;
        int p_offset_via  = num_sig * 12;
        int p_offset_ends = p_offset_via + (num_via * 3);

        MX P_orig = p_mus(Slice(p_offset_ends, p_offset_ends + 3));
        MX P_ins  = p_mus(Slice(p_offset_ends + 3, p_offset_ends + 6));

        // ==========================================================
        // 3. GLOBALE VIA-POINT BERECHNUNG (Wie bisher)
        // ==========================================================
        std::vector<MX> h_via_list; 
        MX sum_h_eta_via = 0;
        // Parameter für Via-Points
        double eps = 1e-8;
        for (int v = 0; v < num_via; ++v) {
            MX v_pos = p_mus(Slice(p_offset_via + v * 3, p_offset_via + (v + 1) * 3));
            MX eta_v = x_mus(offset_via_etas + v);
            
            // HIER HOLEN WIR UNS DEINE ECHTE C++ TOLERANZ!
            double r_tol = via_meshes[v]->MViaPointTolerance; // ohne "0.5" irgendwie abstand zu groß 
            // Fallback
            if (r_tol <= 0.0001) r_tol = 0.005; 
            double r_sq = r_tol * r_tol;

            MX d_vector = MX::vertcat({});
            // for (int k = 0; k < num_inner; ++k) {
            //     MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
            //     d_vector = MX::vertcat({d_vector, sum1(sq(g_k - v_pos))});
            // }
            for (int k = 0; k < num_inner; ++k) {
                MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
                // LINEARE DISTANZ: sqrt( (dx)^2 + (dy)^2 + (dz)^2 + eps )
                MX dist = sqrt(sum1(sq(g_k - v_pos)) + eps);
                d_vector = MX::vertcat({d_vector, dist});
            }
            // geändert am 10.04.26
            /* d_vector = MX::vertcat({d_vector, sum1(sq(P_orig - v_pos))});
            d_vector = MX::vertcat({d_vector, sum1(sq(P_ins - v_pos))}); */
            d_vector = MX::vertcat({d_vector, sqrt(sum1(sq(P_orig - v_pos)) + eps)});
            d_vector = MX::vertcat({d_vector, sqrt(sum1(sq(P_ins - v_pos)) + eps)});

            MX D_v, h_v;
            if (bMaxLogSumTrick){
                // --- NUMERISCH STABILES LOGSUMEXP (MAX-TRICK) ---
                MX x = -Alpha * d_vector;
                // 1. Finde den größten Wert im Vektor (mmax)
                MX x_max = casadi::MX::mmax(x); 
                // 2. Ziehe das Maximum vor dem exp() ab. 
                MX sum_exp = sum1(exp(x - x_max));
                // 3. Logarithmus ziehen und das Maximum wieder addieren
                D_v = -(x_max + log(sum_exp)) / Alpha;
                h_v = r_tol - D_v; 
            }
            else{
                D_v = -casadi::MX::logsumexp(-Alpha * d_vector) / Alpha;
                h_v = r_tol - D_v; // hier eventuell linear statt quadratisch abziehen -> (quadratuisch - linear funkioniert aber irgendiwe besser?)
            }
        
            h_via_list.push_back(h_v);
            sum_h_eta_via += h_v * eta_v;

            // Constraints für Via-Points: NUR NOCH DIE DISTANZ EINZELN HINZUFÜGEN
            all_g = MX::vertcat({all_g, h_v});     
        }
        // NEU: GLOBALES CONSTRAINT FÜR ALLE VIA-POINTS ANFÜGEN
        if (num_via > 0) {
            all_g = MX::vertcat({all_g, sum_h_eta_via});
        }
        // #####################################

        // Gradienten für Via-Points vorab berechnen
        std::vector<MX> grad_h_via_full_list;
        for (int v = 0; v < num_via; ++v) {
            grad_h_via_full_list.push_back(MX::gradient(h_via_list[v], x_mus));
        }

        // ==========================================================
        // 4. KNOTEN-SCHLEIFE (Signorini & Euler-Lagrange)
        // ==========================================================
        for (int k = 0; k < num_inner; ++k) {
            MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
            MX g_prev = (k == 0) ? P_orig : x_mus(Slice((k - 1) * 3, k * 3));
            MX g_next = (k == num_inner - 1) ? P_ins : x_mus(Slice((k + 1) * 3, (k + 2) * 3));

            MX total_sig_force = 0;
            MX total_via_force = 0;
            MX sum_h_eta_sig = 0; 

            // --- A. SIGNORINI KRÄFTE (Klassisch aufsummiert) ---
            for (int s = 0; s < num_sig; ++s) {
                MX q_s = p_mus(Slice(p_offset_sig + s * 12, p_offset_sig + (s + 1) * 12));
                MX eta_ks = x_mus(offset_sig_etas + k * num_sig + s);
                
                MX h_s = sig_meshes[s]->constraintDistance(g_k, q_s);

                MX grad_h_s;
                if (bUseOwnGradient) grad_h_s = sig_meshes[s]->constraintJacobian(g_k, q_s);
                else {
                    MX grad_full = MX::gradient(h_s, x_mus);
                    grad_h_s = grad_full(Slice(k * 3, (k + 1) * 3));
                }

                total_sig_force += eta_ks * grad_h_s;
                sum_h_eta_sig += h_s * eta_ks;
                
                all_g = MX::vertcat({all_g, h_s}); // Non-penetration: h >= 0
            }

            // Summierte Signorini-Komplementarität hinzufügen
            if (num_sig > 0) {
                all_g = MX::vertcat({all_g, sum_h_eta_sig}); // sum(h*eta) == 0
            }

            // --- B. VIA-POINT KRÄFTE ---
            for (int v = 0; v < num_via; ++v) {
                MX eta_v = x_mus(offset_via_etas + v);
                MX grad_h_v = grad_h_via_full_list[v](Slice(k * 3, (k + 1) * 3)); 
                
                total_via_force += eta_v * grad_h_v;
            }

            // --- C. EULER-LAGRANGE KRAFTBILANZ ---
            // Summiere BEIDE Kontaktkraft-Arten in die Lagrange-Gleichung!
            MX total_contact_force = total_sig_force + total_via_force;
            MX eq_el = (-g_prev + 2.0 * g_k - g_next) * K - (1.0 / K) * total_contact_force;
            all_g = MX::vertcat({all_g, eq_el});
        }
    }

    // ==========================================================
    // SOLVER CONFIG (Bleibt unverändert)
    // ==========================================================
    opts["print_time"] = false;
    opts["record_time"] = false;
    opts["ipopt.print_level"] = 0;
    opts["ipopt.max_iter"] = maxIterations;
    opts["ipopt.tol"] = maxTol;
    opts["ipopt.linear_solver"] = "mumps";
    opts["ipopt.hessian_approximation"] = "limited-memory"; 
    
    // ... [Dein Code zum Schreiben der Datei hier] ...

    MXDict nlp = {{"x", all_x}, {"f", obj}, {"g", all_g}, {"p", all_p}};
    solver = nlpsol("f", "ipopt", nlp, opts);

    if (bDebug) qDebug() << "CasadiSystem initialized with"
             << M << "muscles,"
             << "Total variables:" << all_x.size1()
             << "Total constraints:" << all_g.size1();
}




// VIA LINE POINTS
void CasadiSystem::solveStepViaLine() {
    qDebug() << "Solving step with ViaLine Point formulation...";
    using namespace casadi;

    std::vector<double> x0_all, p_all, lbg_all, ubg_all, lbx_all, ubx_all;

    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;

        // ==========================================================
        // 0. MESHES SORTIEREN
        // ==========================================================
        std::vector<SSMesh*> sig_meshes;
        std::vector<SSMesh*> via_meshes;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) via_meshes.push_back(mesh);
            else sig_meshes.push_back(mesh);
        }
        
        int num_sig = sig_meshes.size();
        int num_via = via_meshes.size();
        int num_etas_sig = num_inner * num_sig;
        int num_etas_via = num_via;
        int num_etas_total = num_etas_sig + num_etas_via;

        // ==========================================================
        // 1. INITIAL GUESS (x0)
        // ==========================================================
        for (size_t k = 1; k < mus->MNodes.size() - 1; ++k) {
            MWMath::Point3D pWarm = mus->MNodes[k].predictNewGlobal();
            x0_all.push_back(pWarm.x);
            x0_all.push_back(pWarm.y);
            x0_all.push_back(pWarm.z);
        }
        
        // Warmstart für BEIDE Eta-Arten
        if (mus->lastEtas.size() == (size_t)num_etas_total && bUseWarmstartEtas) {
            for (double eta : mus->lastEtas) {
                x0_all.push_back(eta * WarmstartEtaScaling);
            }
        } else {
            for (int e = 0; e < num_etas_total; ++e) x0_all.push_back(0.0);
        }
        

        // ==========================================================
        // 2. PARAMETER (p)
        // ==========================================================
        // A. Signorini Meshes (12 Parameter: Pos + Rot Matrix)
        for (auto* mesh : sig_meshes) {
            p_all.push_back(mesh->PositionGlobal.x);
            p_all.push_back(mesh->PositionGlobal.y);
            p_all.push_back(mesh->PositionGlobal.z);
            
            MWMath::RotMatrix3x3 R = mesh->OrientationGlobal;
            for (int col = 0; col < 3; ++col) {
                for (int row = 0; row < 3; ++row) {
                    p_all.push_back(R.m[row][col]);
                }
            }
        }
        
        // B. Via-Point Meshes (Nur 3 Parameter: Position)
        for (auto* mesh : via_meshes) {
            p_all.push_back(mesh->PositionGlobal.x);
            p_all.push_back(mesh->PositionGlobal.y);
            p_all.push_back(mesh->PositionGlobal.z);
        }
        
        // C. Origin & Insertion
        p_all.push_back(mus->OriginPointGlobal.x);
        p_all.push_back(mus->OriginPointGlobal.y);
        p_all.push_back(mus->OriginPointGlobal.z);
        
        p_all.push_back(mus->InsertionPointGlobal.x);
        p_all.push_back(mus->InsertionPointGlobal.y);
        p_all.push_back(mus->InsertionPointGlobal.z);

        // ==========================================================
        // 3. CONSTRAINT BOUNDS (lbg/ubg)
        // ==========================================================
        // A. Globale Via-Points
        for (int v = 0; v < num_via; ++v) {
            // Komplementarität (Unterdruck-Trick)
            lbg_all.push_back(-inf); 
            ubg_all.push_back(0.0); // WICHTIG: 0.0001 Slack für lineare Distanzen!
            
            // Distanz-Constraint: h_v >= 0
            lbg_all.push_back(0.0); 
            ubg_all.push_back(inf);
        }

        // B. Knoten-Schleife (Signorini + Euler-Lagrange)
        for (int k = 0; k < num_inner; ++k) {
            // Non-penetration für jedes Signorini-Mesh (h_s >= 0)
            for (int s = 0; s < num_sig; ++s) {
                lbg_all.push_back(0.0); 
                ubg_all.push_back(inf);
            }
            // Summierte Signorini-Komplementarität (sum(h*eta) == 0)
            if (num_sig > 0) {
                lbg_all.push_back(0.0); 
                ubg_all.push_back(0.0);
            }
            // Euler-Lagrange Kraftbilanz (x, y, z == 0)
            for (int d = 0; d < 3; ++d) { 
                lbg_all.push_back(0.0); 
                ubg_all.push_back(0.0); 
            }
        }

        // ==========================================================
        // 4. VARIABLE BOUNDS (lbx/ubx)
        // ==========================================================
        // Gamma (Punkte)
        for (int i = 0; i < num_inner * 3; ++i) {
            lbx_all.push_back(-inf); ubx_all.push_back(inf);
        }
        // Alle Etas (Signorini + Via-Points) dürfen nur ziehen/drücken (>= 0)
        for (int i = 0; i < num_etas_total; ++i) {
            lbx_all.push_back(0.0); ubx_all.push_back(inf);
        }
    }

    // ==========================================================
    // SOLVER AUFRUF
    // ==========================================================
    DMDict arg = {{"x0", x0_all}, {"p", p_all}, {"lbg", lbg_all}, {"ubg", ubg_all}, {"lbx", lbx_all}, {"ubx", ubx_all}};
    DMDict res = solver(arg);

    // ==========================================================
    // ERGEBNISSE EXTRAHIEREN
    // ==========================================================
    std::vector<double> res_x = std::vector<double>(res.at("x"));
    auto solverInfo = solver.stats();
    std::string status = solverInfo.at("return_status");

    const std::string C_RED   = "\033[31m";
    const std::string C_GREEN = "\033[32m";
    const std::string C_RESET = "\033[0m";
    int convSteps = int(solverInfo.at("iter_count"));
    if (true) { // Dein if(true) oder if(bDebug)
        std::string stepColor = (status == "Solve_Succeeded") ? C_GREEN : C_RED;
        std::string fullMessage = stepColor + "    Solver finished after " + std::to_string(convSteps) + " iterations" + "(" + status + ")" + C_RESET;
        qDebug().noquote() << QString::fromStdString(fullMessage);
    }
    
    SolverConvergenceMessages.push_back(status);
    SolverConvergenceSteps.push_back(convSteps);
    
    int current_x_offset = 0;
    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;
        int num_sig = 0, num_via = 0;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) num_via++; else num_sig++;
        }
        
        int offset_nodes    = current_x_offset;
        int offset_sig_etas = offset_nodes + (num_inner * 3);
        int offset_via_etas = offset_sig_etas + (num_inner * num_sig);

        // 1. Positionen zurückschreiben
        for (int i = 0; i < num_inner; ++i) {
            int mIdx = i + 1;
            MWMath::Point3D p(
                res_x[offset_nodes + i * 3 + 0],
                res_x[offset_nodes + i * 3 + 1],
                res_x[offset_nodes + i * 3 + 2]
            );
            mus->MusclePointsGlobal[mIdx] = p;
            mus->MNodes[mIdx].PositionGlobal = p;
            mus->MNodes[mIdx].updateLocalFrame();
        }

        // 2. Globale Etas für Warmstart speichern
        int num_etas_total = (num_inner * num_sig) + num_via;
        mus->lastEtas.assign(res_x.begin() + offset_sig_etas, res_x.begin() + offset_sig_etas + num_etas_total);

        // ==========================================================
        // 3. MAPPING-TRICK FÜR EXPORT/VISUALISIERUNG
        // Wir dröseln die Ergebnisse für die alte Struktur wieder auf!
        // ==========================================================
        
        // Zuerst: Welcher Knoten ist jedem Via-Point am nächsten?
        std::vector<int> closest_node_for_via(num_via, -1);
        int via_counter = 0;
        for (auto* mesh : mus->meshPtrs) {
            if (!mesh->bIsViaPoint) continue;
            
            MWMath::Point3D v_pos = mesh->PositionGlobal;
            double min_dist = 1e9;
            int best_k = 0;
            
            for (int k = 0; k < num_inner; ++k) {
                double d = MWMath::distance(mus->MNodes[k + 1].PositionGlobal, v_pos);
                if (d < min_dist) { min_dist = d; best_k = k; }
            }
            closest_node_for_via[via_counter] = best_k;
            via_counter++;
        }

        // Nun füllen wir das Array für JEDEN Knoten
        for (int k = 0; k < num_inner; ++k) {
            int mIdx = k + 1;
            std::vector<double> currentStepEtas(mus->meshPtrs.size(), 0.0);

            int sig_idx = 0;
            int via_idx = 0;

            for (size_t m = 0; m < mus->meshPtrs.size(); ++m) {
                if (mus->meshPtrs[m]->bIsViaPoint) {
                    // Via-Point Kraft nur eintragen, wenn das der nächste Knoten ist
                    if (closest_node_for_via[via_idx] == k) {
                        currentStepEtas[m] = res_x[offset_via_etas + via_idx];
                    }
                    via_idx++;
                } else {
                    // Signorini Kraft ganz regulär auslesen
                    currentStepEtas[m] = res_x[offset_sig_etas + (k * num_sig) + sig_idx];
                    sig_idx++;
                }
            }
            
            mus->MNodes[mIdx].MNodeEtaSteps.push_back(currentStepEtas);
            mus->MNodes[mIdx].lastEtas = mus->lastEtas; 
        }

        current_x_offset += (num_inner * 3) + num_etas_total;
    }
}

void CasadiSystem::setupCasadiViaLine()
{
    qDebug() << "     Setting up CasadiSystem with ViaLine Point formulation...";
    using namespace casadi;
    MX all_x = MX::vertcat({});
    MX all_g = MX::vertcat({});
    MX all_p = MX::vertcat({});
    MX obj = 0;
    
    

    for (size_t m = 0; m < m_muscles.size(); ++m) {
        SSMuscle* mus = m_muscles[m];
        int num_inner = mus->MNodes.size() - 2;
        int K = mus->MNodes.size() - 1;

        // ==========================================================
        // 0. MESHES SORTIEREN
        // ==========================================================
        std::vector<SSMesh*> sig_meshes;
        std::vector<SSMesh*> via_meshes;
        for (auto* mesh : mus->meshPtrs) {
            if (mesh->bIsViaPoint) via_meshes.push_back(mesh);
            else sig_meshes.push_back(mesh);
        }
        int num_sig = sig_meshes.size();
        int num_via = via_meshes.size();

        // ==========================================================
        // 1. VARIABLEN DEFINIEREN (X-Vektor)
        // Aufbau: [ 3*N Knoten | N * M_sig Etas | M_via Etas ]
        // ==========================================================
        int num_etas_sig = num_inner * num_sig;
        int num_etas_via = num_via;
        
        MX x_mus = MX::sym("x_" + std::to_string(m), (num_inner * 3) + num_etas_sig + num_etas_via);
        all_x = MX::vertcat({all_x, x_mus});

        // Offsets zum Herausschneiden aus x_mus
        int offset_nodes    = 0;
        int offset_sig_etas = num_inner * 3;
        int offset_via_etas = offset_sig_etas + num_etas_sig;

        // ==========================================================
        // 2. PARAMETER DEFINIEREN (P-Vektor)
        // Aufbau: [ 12*M_sig Parameter | 3*M_via Parameter | 3 Orig | 3 Ins ]
        // ==========================================================
        MX p_mus = MX::sym("p_" + std::to_string(m), (12 * num_sig) + (3 * num_via) + 6);
        all_p = MX::vertcat({all_p, p_mus});

        // Offsets zum Herausschneiden aus p_mus
        int p_offset_sig  = 0;
        int p_offset_via  = num_sig * 12;
        int p_offset_ends = p_offset_via + (num_via * 3);

        MX P_orig = p_mus(Slice(p_offset_ends, p_offset_ends + 3));
        MX P_ins  = p_mus(Slice(p_offset_ends + 3, p_offset_ends + 6));

        // ==========================================================
        // 3. GLOBALE VIA-POINT BERECHNUNG (LINIEN-SEGMENTE!)
        // ==========================================================
        
        std::vector<MX> h_via_list; 
        Alpha = 150.0; // Besserer Wert für Linien
        double eps = 1e-8;

        std::vector<MX> path_nodes;
        path_nodes.push_back(P_orig);
        for (int k = 0; k < num_inner; ++k) {
            path_nodes.push_back(x_mus(Slice(k * 3, (k + 1) * 3)));
        }
        path_nodes.push_back(P_ins);

        for (int v = 0; v < num_via; ++v) {
            MX v_pos = p_mus(Slice(p_offset_via + v * 3, p_offset_via + (v + 1) * 3));
            MX eta_v = x_mus(offset_via_etas + v);
            
            double r_tol = via_meshes[v]->MViaPointTolerance * 0.5; 
            if (r_tol <= 0.0001) r_tol = 0.005; 

            MX d_vector = MX::vertcat({});
            
            for (size_t i = 0; i < path_nodes.size() - 1; ++i) {
                MX p1 = path_nodes[i];
                MX p2 = path_nodes[i+1];
                
                MX u = p2 - p1;          
                MX w = v_pos - p1;       
                
                MX c1 = sum1(w * u);
                MX c2 = sum1(u * u) + 1e-6; 
                MX t = c1 / c2;
                
                double eps_smooth = 1e-5; 
                MX t_min = 0.5 * (t + 1.0 - sqrt(sq(t - 1.0) + eps_smooth));
                MX t_clamped = 0.5 * (t_min + sqrt(sq(t_min) + eps_smooth));
                
                MX closest_point = p1 + t_clamped * u;
                MX dist = sqrt(sum1(sq(closest_point - v_pos)) + eps);
                d_vector = MX::vertcat({d_vector, dist});
            }

            // HIER HAT DIE KORREKTUR GEFEHLT! (+ log(n_segments)/Alpha)
            double n_segments = path_nodes.size() - 1.0;
            MX D_v = -casadi::MX::logsumexp(-Alpha * d_vector) / Alpha + (log(n_segments) / Alpha);
            
            MX h_v = r_tol - D_v; 
            h_via_list.push_back(h_v);

            // Regularisierung für eta_v
            obj += 1e-5 * sq(eta_v);

            all_g = MX::vertcat({all_g, h_v * eta_v}); 
            all_g = MX::vertcat({all_g, h_v});         
        }

        // HIER HAT DIE LÄNGEN-REGULARISIERUNG GEFEHLT!
        // Verhindert, dass Knoten unendlich weit auseinander rutschen
        for (size_t i = 0; i < path_nodes.size() - 1; ++i) {
            obj += 1e-5 * sum1(sq(path_nodes[i] - path_nodes[i+1]));
        }
        

        std::vector<MX> grad_h_via_full_list;
        for (int v = 0; v < num_via; ++v) {
            grad_h_via_full_list.push_back(MX::gradient(h_via_list[v], x_mus));
        }

        // ==========================================================
        // 4. KNOTEN-SCHLEIFE (Signorini & Euler-Lagrange)
        // ==========================================================
        for (int k = 0; k < num_inner; ++k) {
            MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
            MX g_prev = (k == 0) ? P_orig : x_mus(Slice((k - 1) * 3, k * 3));
            MX g_next = (k == num_inner - 1) ? P_ins : x_mus(Slice((k + 1) * 3, (k + 2) * 3));

            MX total_sig_force = 0;
            MX total_via_force = 0;
            MX sum_h_eta_sig = 0; 

            // --- A. SIGNORINI KRÄFTE (Klassisch aufsummiert) ---
            for (int s = 0; s < num_sig; ++s) {
                MX q_s = p_mus(Slice(p_offset_sig + s * 12, p_offset_sig + (s + 1) * 12));
                MX eta_ks = x_mus(offset_sig_etas + k * num_sig + s);
                
                MX h_s = sig_meshes[s]->constraintDistance(g_k, q_s);

                MX grad_h_s;
                if (bUseOwnGradient) grad_h_s = sig_meshes[s]->constraintJacobian(g_k, q_s);
                else {
                    MX grad_full = MX::gradient(h_s, x_mus);
                    grad_h_s = grad_full(Slice(k * 3, (k + 1) * 3));
                }

                total_sig_force += eta_ks * grad_h_s;
                sum_h_eta_sig += h_s * eta_ks;
                
                all_g = MX::vertcat({all_g, h_s}); // Non-penetration: h >= 0
            }

            // Summierte Signorini-Komplementarität hinzufügen
            if (num_sig > 0) {
                all_g = MX::vertcat({all_g, sum_h_eta_sig}); // sum(h*eta) == 0
            }

            // --- B. VIA-POINT KRÄFTE ---
            for (int v = 0; v < num_via; ++v) {
                MX eta_v = x_mus(offset_via_etas + v);
                MX grad_h_v = grad_h_via_full_list[v](Slice(k * 3, (k + 1) * 3)); 
                
                total_via_force += eta_v * grad_h_v;
            }

            // --- C. EULER-LAGRANGE KRAFTBILANZ ---
            // Summiere BEIDE Kontaktkraft-Arten in die Lagrange-Gleichung!
            MX total_contact_force = total_sig_force + total_via_force;
            MX eq_el = (-g_prev + 2.0 * g_k - g_next) * K - (1.0 / K) * total_contact_force;
            all_g = MX::vertcat({all_g, eq_el});
        }
    }

    // ==========================================================
    // SOLVER CONFIG (Bleibt unverändert)
    // ==========================================================
    opts["print_time"] = false;
    opts["record_time"] = false;
    opts["ipopt.print_level"] = 0;
    opts["ipopt.max_iter"] = maxIterations;
    opts["ipopt.tol"] = maxTol;
    opts["ipopt.linear_solver"] = "mumps";
    opts["ipopt.hessian_approximation"] = "limited-memory"; 
    
    // ... [Dein Code zum Schreiben der Datei hier] ...

    MXDict nlp = {{"x", all_x}, {"f", obj}, {"g", all_g}, {"p", all_p}};
    solver = nlpsol("f", "ipopt", nlp, opts);

    if (bDebug) qDebug() << "CasadiSystem initialized with"
             << M << "muscles,"
             << "Total variables:" << all_x.size1()
             << "Total constraints:" << all_g.size1();
}




// NORMAL
void CasadiSystem::setupCasadi()
{
    using namespace casadi;
    MX all_x = MX::vertcat({});
    MX all_g = MX::vertcat({});
    MX all_p = MX::vertcat({});
    MX obj = 0;
    M = static_cast<int>(m_muscles.size());

    for (size_t m = 0; m < m_muscles.size(); ++m) {
        SSMuscle* mus = m_muscles[m];
        int num_inner = mus->MNodes.size() - 2;
        int num_wrap = mus->meshPtrs.size();
        int num_eta_per_node = num_wrap;

        // --- VARIABLEN: [Punkte (3*n), Etas (n * Meshes)] ---
        MX x_mus = MX::sym("x_" + std::to_string(m), num_inner * 3 + num_inner * num_eta_per_node);
        all_x = MX::vertcat({all_x, x_mus});

        // --- PARAMETER: [Meshes (12*m), Origin (3), Insertion (3)] ---
        MX p_mus = MX::sym("p_" + std::to_string(m), 12 * num_wrap + 6);
        all_p = MX::vertcat({all_p, p_mus});

        MX P_orig = p_mus(Slice(12 * num_wrap, 12 * num_wrap + 3));
        MX P_ins  = p_mus(Slice(12 * num_wrap + 3, 12 * num_wrap + 6));
        int eta_offset = num_inner * 3;

        // --- GLEICHUNGEN (g) ---
        for (int k = 0; k < num_inner; ++k) {
            MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
            MX g_prev = (k == 0) ? P_orig : x_mus(Slice((k - 1) * 3, k * 3));
            MX g_next = (k == num_inner - 1) ? P_ins : x_mus(Slice((k + 1) * 3, (k + 2) * 3));

            MX total_contact_force = 0;
            int current_eta_idx = 0;

            for (int j = 0; j < num_wrap; ++j) {
                MX q_j = p_mus(Slice(j * 12, (j + 1) * 12));

                // Normale Signorini-Logik (Hindernisse)
                MX eta_kj = x_mus(eta_offset + k * num_eta_per_node + current_eta_idx);
                MX h = mus->meshPtrs[j]->constraintDistance(g_k, q_j);

                MX grad_h;
                if (bUseOwnGradient){
                    grad_h = mus->meshPtrs[j]->constraintJacobian(g_k, q_j);
                    //qDebug() << "using own gradient!";
                }
                else {
                    MX grad_h_full = MX::gradient(h, x_mus);
                    grad_h = grad_h_full(Slice(k * 3, (k + 1) * 3));
                    //qDebug() << "using auto-casadi gradient!";
                }

                total_contact_force += eta_kj * grad_h;
                
                // Signorini/Complementarity
                all_g = MX::vertcat({all_g, h * eta_kj, h, eta_kj}); 
                current_eta_idx++;
                
            }

            // Euler-Lagrange Kraftbilanz
            MX eq_el = (- g_prev + 2.0*g_k - g_next)*num_inner - (total_contact_force); // (g_next - 2.0 * g_k + g_prev) + total_contact_force;
            all_g = MX::vertcat({all_g, eq_el});
            
        }
    }

    opts["print_time"] = false;
    opts["record_time"] = false;
    opts["ipopt.print_level"] = 0;
    opts["ipopt.max_iter"] = maxIterations;
    opts["ipopt.tol"] = maxTol;
    opts["ipopt.linear_solver"] = "mumps";
    opts["ipopt.hessian_approximation"] = "limited-memory";
    std::string filename = "../examples/results/solver_log_" + CasadiSystemName + ".txt";
    opts["ipopt.output_file"] = filename; 
    opts["ipopt.file_print_level"] = 5;

    MXDict nlp = {{"x", all_x}, {"f", obj}, {"g", all_g}, {"p", all_p}};
    solver = nlpsol("f", "ipopt", nlp, opts);

    if (bDebug) qDebug() << "CasadiSystem initialized with"
             << M << "muscles,"
             << "Objective Type:" << objType
             << "Total variables:" << all_x.size1()
             << "Total constraints:" << all_g.size1();
}

void CasadiSystem::solveStep() {
    using namespace casadi;

    std::vector<double> x0_all, p_all, lbg_all, ubg_all, lbx_all, ubx_all;
    std::vector<double> inputParamsForDebug;
    std::vector<std::string> inputParamDescriptionsForDebug;

    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;
        int num_wrap = mus->meshPtrs.size();
        int num_etas_total = num_inner * num_wrap;

        // --- INITIAL GUESS (x0) ---
        for (size_t k = 1; k < mus->MNodes.size() - 1; ++k) {
            MWMath::Point3D pWarm = mus->MNodes[k].predictNewGlobal();
            x0_all.push_back(pWarm.x);
            x0_all.push_back(pWarm.y);
            x0_all.push_back(pWarm.z);
        }
        
        if (mus->lastEtas.size() == (size_t)num_etas_total && bUseWarmstartEtas) {
            if (bDebug) qDebug() << "    Warmstart etas:" << QString::fromStdString(mus->Name);
            x0_all.insert(x0_all.end(), mus->lastEtas.begin(), mus->lastEtas.end());
        } else {
            if (bDebug) qDebug() << "    No warmstart etas:" << QString::fromStdString(mus->Name);
            for (int e = 0; e < num_etas_total; ++e) x0_all.push_back(0.0);
        }

        // --- PARAMETER (p) ---
        for (auto* mesh : mus->meshPtrs) {
            p_all.push_back(mesh->PositionGlobal.x);    inputParamsForDebug.push_back(mesh->PositionGlobal.x);  inputParamDescriptionsForDebug.push_back(mesh->Name + "_PosX");
            p_all.push_back(mesh->PositionGlobal.y);    inputParamsForDebug.push_back(mesh->PositionGlobal.y);  inputParamDescriptionsForDebug.push_back(mesh->Name + "_PosY");
            p_all.push_back(mesh->PositionGlobal.z);    inputParamsForDebug.push_back(mesh->PositionGlobal.z);  inputParamDescriptionsForDebug.push_back(mesh->Name + "_PosZ");
            
           MWMath::RotMatrix3x3 R = mesh->OrientationGlobal;
           for (int col = 0; col < 3; ++col) {       // <--- Spalte ist außen!
                for (int row = 0; row < 3; ++row) {   // <--- Zeile ist innen!
                    // Wir greifen aber ganz normal auf m[row][col] zu
                    p_all.push_back(R.m[row][col]);
                    inputParamsForDebug.push_back(R.m[row][col]);
                    inputParamDescriptionsForDebug.push_back(mesh->Name + "_Ori_" + std::to_string(row) + std::to_string(col));
                }
            }
        }
        p_all.push_back(mus->OriginPointGlobal.x);    inputParamsForDebug.push_back(mus->OriginPointGlobal.x);  inputParamDescriptionsForDebug.push_back(mus->Name + "_OriginX");
        p_all.push_back(mus->OriginPointGlobal.y);    inputParamsForDebug.push_back(mus->OriginPointGlobal.y);  inputParamDescriptionsForDebug.push_back(mus->Name + "_OriginY");
        p_all.push_back(mus->OriginPointGlobal.z);    inputParamsForDebug.push_back(mus->OriginPointGlobal.z);  inputParamDescriptionsForDebug.push_back(mus->Name + "_OriginZ");
        p_all.push_back(mus->InsertionPointGlobal.x);     inputParamsForDebug.push_back(mus->InsertionPointGlobal.x);  inputParamDescriptionsForDebug.push_back(mus->Name + "_InsertionX");
        p_all.push_back(mus->InsertionPointGlobal.y);     inputParamsForDebug.push_back(mus->InsertionPointGlobal.y);  inputParamDescriptionsForDebug.push_back(mus->Name + "_InsertionY");
        p_all.push_back(mus->InsertionPointGlobal.z);     inputParamsForDebug.push_back(mus->InsertionPointGlobal.z);  inputParamDescriptionsForDebug.push_back(mus->Name + "_InsertionZ");
        // --- CONSTRAINTS SCHRANKEN (lbg/ubg) ---
        for (int i = 0; i < num_inner; ++i) {
            for (int j = 0; j < num_wrap; ++j) {
                // Komplementarität: h * eta = 0
                lbg_all.push_back(0.0); ubg_all.push_back(0.0); 
                // Nicht-Eindringen: h >= 0
                lbg_all.push_back(0.0); ubg_all.push_back(inf);
                // Nur Druckkraft: eta >= 0
                lbg_all.push_back(0.0); ubg_all.push_back(inf);
            }
            // EL = 0 (für x, y, z)
            for (int d = 0; d < 3; ++d) { 
                lbg_all.push_back(0.0); ubg_all.push_back(0.0); 
            }
        }

        // --- VARIABLEN SCHRANKEN (lbx/ubx) ---
        // Gamma (Punkte) sind frei
        for (int i = 0; i < num_inner * 3; ++i) {
            lbx_all.push_back(-inf); 
            ubx_all.push_back(inf);
        }
        // Eta (Kräfte) >= 0
        for (int i = 0; i < num_etas_total; ++i) {
            lbx_all.push_back(0.0); 
            ubx_all.push_back(inf);
        }
    }

    allParameterInputsAllSteps.push_back(inputParamsForDebug);
    allParameterInputDescriptionsAllSteps.push_back(inputParamDescriptionsForDebug);

    // --- SOLVER AUFRUF ---
    DMDict arg = {{"x0", x0_all}, {"p", p_all}, {"lbg", lbg_all}, {"ubg", ubg_all}, {"lbx", lbx_all}, {"ubx", ubx_all}};
    DMDict res = solver(arg);

    // --- ERGEBNISSE EXTRAHIEREN ---
    std::vector<double> res_x = std::vector<double>(res.at("x"));
    auto solverInfo = solver.stats();
    std::string status = solverInfo.at("return_status");

    const char* C_RED   = "\033[31m";
    const char* C_GREEN = "\033[32m";
    const char* C_RESET = "\033[0m";

    int convSteps = int(solverInfo.at("iter_count"));
    qDebug() << "Solver status:" << QString::fromStdString(status);
    if (bDebug) {
        // Du kannst auch die Anzahl der Iterationen passend einfärben:
        const char* stepColor = (status == "Solve_Succeeded") ? C_GREEN : C_RED;
        qDebug().noquote() << stepColor << "    Solver stopped after" << convSteps << "iterations" << C_RESET;
    }

    SolverConvergenceMessages.push_back(QString::fromStdString(status).toStdString());
    SolverConvergenceSteps.push_back(convSteps);
    
    int current_x_offset = 0;
    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;
        int num_wrap = mus->meshPtrs.size();
        int num_etas_mus = num_inner * (num_wrap);

        // Gamma (Punkte) zurückschreiben
        for (int i = 0; i < num_inner; ++i) {
            int mIdx = i + 1;
            MWMath::Point3D p;
            p.x = res_x[current_x_offset + i * 3 + 0];
            p.y = res_x[current_x_offset + i * 3 + 1];
            p.z = res_x[current_x_offset + i * 3 + 2];

            mus->MusclePointsGlobal[mIdx] = p;
            mus->MNodes[mIdx].PositionGlobal = p;
            mus->MNodes[mIdx].updateLocalFrame();
        }

        // Etas für Warmstart speichern
        int eta_start_idx = current_x_offset + (num_inner * 3);
        mus->lastEtas.assign(
            res_x.begin() + eta_start_idx, 
            res_x.begin() + eta_start_idx + num_etas_mus
        );

        // ETAS IN LISTE SCHRIBEN (für späteren Export)
        int num_eta_per_node = (num_inner > 0) ? (num_etas_mus / num_inner) : 0;
        for (int k = 0; k < num_inner; ++k) {
            int mIdx = k + 1; // Wir fangen bei Node 1 an (0 ist fix)
            
            std::vector<double> currentStepEtas;
            
            // Die Etas für diesen spezifischen Node aus dem großen Vektor holen
            for (int e = 0; e < num_eta_per_node; ++e) {
                // Index im 'lastEtas' Vektor berechnen
                int globalEtaIndex = k * num_eta_per_node + e;
                if (globalEtaIndex < mus->lastEtas.size()) {
                    currentStepEtas.push_back(mus->lastEtas[globalEtaIndex]);
                } else {
                    currentStepEtas.push_back(0.0);
                }
            }

            // Speichern für den späteren Export
            mus->MNodes[mIdx].MNodeEtaSteps.push_back(currentStepEtas);
            
            // Auch das 'lastEtas' Feld im Node updaten (für Debugging)
            mus->MNodes[mIdx].lastEtas = mus->lastEtas; 
        }

        for (auto& node : mus->MNodes) {
            node.lastEtas = mus->lastEtas;
        }

        current_x_offset += (num_inner * 3 + num_etas_mus);
    }
}


// alt test
void CasadiSystem::setupCasadiSumAlt()
{
    qDebug() << "Setting up CasadiSystem with Sum Phi*Eta formulation...";
    using namespace casadi;
    MX all_x = MX::vertcat({});
    MX all_g = MX::vertcat({});
    MX all_p = MX::vertcat({});
    MX obj = 0;
    M = static_cast<int>(m_muscles.size());

    for (size_t m = 0; m < m_muscles.size(); ++m) {
        SSMuscle* mus = m_muscles[m];
        int num_inner = mus->MNodes.size() - 2;
        int K = mus->MNodes.size() -1;
        int num_wrap = mus->meshPtrs.size();
        int num_eta_per_node = num_wrap;

        // --- VARIABLEN: [Punkte (3*n), Etas (n * Meshes)] ---
        MX x_mus = MX::sym("x_" + std::to_string(m), num_inner * 3 + num_inner * num_eta_per_node);
        all_x = MX::vertcat({all_x, x_mus});

        // --- PARAMETER: [Meshes (12*m), Origin (3), Insertion (3)] ---
        MX p_mus = MX::sym("p_" + std::to_string(m), 12 * num_wrap + 6);
        all_p = MX::vertcat({all_p, p_mus});

        MX P_orig = p_mus(Slice(12 * num_wrap, 12 * num_wrap + 3));
        MX P_ins  = p_mus(Slice(12 * num_wrap + 3, 12 * num_wrap + 6));
        int eta_offset = num_inner * 3;

        // --- CONSTRAINTS FÜR JEDEN KNOTEN ---
        for (int k = 0; k < num_inner; ++k) {
            MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
            MX g_prev = (k == 0) ? P_orig : x_mus(Slice((k - 1) * 3, k * 3));
            MX g_next = (k == num_inner - 1) ? P_ins : x_mus(Slice((k + 1) * 3, (k + 2) * 3));

            MX total_contact_force = 0;
            MX sum_h_eta = 0;  // Summierte Komplementarität

            // Über alle Wrapping-Meshes iterieren
            for (int j = 0; j < num_wrap; ++j) {
                MX q_j = p_mus(Slice(j * 12, (j + 1) * 12));
                MX eta_kj = x_mus(eta_offset + k * num_eta_per_node + j);
                
                // Distanz-Constraint berechnen
                MX h = mus->meshPtrs[j]->constraintDistance(g_k, q_j);

                // Gradient berechnen
                MX grad_h;
                if (bUseOwnGradient){
                    grad_h = mus->meshPtrs[j]->constraintJacobian(g_k, q_j);
                }
                else {
                    MX grad_h_full = MX::gradient(h, x_mus);
                    grad_h = grad_h_full(Slice(k * 3, (k + 1) * 3));
                }

                // Kontaktkraft akkumulieren
                total_contact_force += eta_kj * grad_h;
                
                // Komplementaritätsterme akkumulieren
                sum_h_eta += h * eta_kj;
                
                // EINZELNE CONSTRAINTS hinzufügen:
                // 1. Nicht-Eindringen: h >= 0
                all_g = MX::vertcat({all_g, h});
                
                // 2. Nur Druckkraft: eta >= 0 (wird durch Variablenbounds gehandhabt)
            }

            // SUMMIERTE KOMPLEMENTARITÄT: sum(h*eta) = 0
            all_g = MX::vertcat({all_g, sum_h_eta});

            // EULER-LAGRANGE GLEICHUNG (3 Komponenten für x, y, z)
            MX eq_el = (- g_prev + 2.0*g_k - g_next)*K - (1.0/K)*total_contact_force;
            all_g = MX::vertcat({all_g, eq_el});
        }
    }

    opts["print_time"] = false;
    opts["record_time"] = false;
    opts["ipopt.print_level"] = 0;
    opts["ipopt.max_iter"] = maxIterations;
    opts["ipopt.tol"] = maxTol;
    opts["ipopt.linear_solver"] = "mumps";
    opts["ipopt.hessian_approximation"] = "limited-memory"; // "limited-memory" or "exact"
    std::string filename = "../examples/results/solver_log_" + CasadiSystemName + ".txt";
    opts["ipopt.output_file"] = filename; 
    opts["ipopt.file_print_level"] = 5;

    MXDict nlp = {{"x", all_x}, {"f", obj}, {"g", all_g}, {"p", all_p}};
    solver = nlpsol("f", "ipopt", nlp, opts);

    qDebug() << "CasadiSystem initialized with"
             << M << "muscles,"
             << "Objective Type:" << objType
             << "Total variables:" << all_x.size1()
             << "Total constraints:" << all_g.size1();
}

void CasadiSystem::solveStepSumAlt() {
    qDebug() << "Solving step with Sum Phi*Eta formulation...";
    using namespace casadi;

    std::vector<double> x0_all, p_all, lbg_all, ubg_all, lbx_all, ubx_all;
    std::vector<double> inputParamsForDebug;
    std::vector<std::string> inputParamDescriptionsForDebug;

    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;
        int num_wrap = mus->meshPtrs.size();
        int num_etas_total = num_inner * num_wrap;

        // --- INITIAL GUESS (x0) ---
        for (size_t k = 1; k < mus->MNodes.size() - 1; ++k) {
            MWMath::Point3D pWarm = mus->MNodes[k].predictNewGlobal();
            x0_all.push_back(pWarm.x);
            x0_all.push_back(pWarm.y);
            x0_all.push_back(pWarm.z);
        }
        
        if (mus->lastEtas.size() == (size_t)num_etas_total && bUseWarmstartEtas) {
            qDebug() << "    Warmstart etas:" << QString::fromStdString(mus->Name);
            std::vector<double> scaledEtas;
            for (double eta : mus->lastEtas) {
                double scaledEta = eta * WarmstartEtaScaling;
                scaledEtas.push_back(scaledEta);
            }
            x0_all.insert(x0_all.end(), scaledEtas.begin(), scaledEtas.end());
        } else {
            qDebug() << "    No warmstart etas:" << QString::fromStdString(mus->Name);
            for (int e = 0; e < num_etas_total; ++e) x0_all.push_back(0.0);
        }

        // --- PARAMETER (p) ---
        for (auto* mesh : mus->meshPtrs) {
            p_all.push_back(mesh->PositionGlobal.x);    inputParamsForDebug.push_back(mesh->PositionGlobal.x);  inputParamDescriptionsForDebug.push_back(mesh->Name + "_PosX");
            p_all.push_back(mesh->PositionGlobal.y);    inputParamsForDebug.push_back(mesh->PositionGlobal.y);  inputParamDescriptionsForDebug.push_back(mesh->Name + "_PosY");
            p_all.push_back(mesh->PositionGlobal.z);    inputParamsForDebug.push_back(mesh->PositionGlobal.z);  inputParamDescriptionsForDebug.push_back(mesh->Name + "_PosZ");
            
            MWMath::RotMatrix3x3 R = mesh->OrientationGlobal;
            for (int col = 0; col < 3; ++col) {
                for (int row = 0; row < 3; ++row) {
                    p_all.push_back(R.m[row][col]);
                    inputParamsForDebug.push_back(R.m[row][col]);
                    inputParamDescriptionsForDebug.push_back(mesh->Name + "_Ori_" + std::to_string(row) + std::to_string(col));
                }
            }
        }
        p_all.push_back(mus->OriginPointGlobal.x);    inputParamsForDebug.push_back(mus->OriginPointGlobal.x);  inputParamDescriptionsForDebug.push_back(mus->Name + "_OriginX");
        p_all.push_back(mus->OriginPointGlobal.y);    inputParamsForDebug.push_back(mus->OriginPointGlobal.y);  inputParamDescriptionsForDebug.push_back(mus->Name + "_OriginY");
        p_all.push_back(mus->OriginPointGlobal.z);    inputParamsForDebug.push_back(mus->OriginPointGlobal.z);  inputParamDescriptionsForDebug.push_back(mus->Name + "_OriginZ");
        p_all.push_back(mus->InsertionPointGlobal.x); inputParamsForDebug.push_back(mus->InsertionPointGlobal.x);  inputParamDescriptionsForDebug.push_back(mus->Name + "_InsertionX");
        p_all.push_back(mus->InsertionPointGlobal.y); inputParamsForDebug.push_back(mus->InsertionPointGlobal.y);  inputParamDescriptionsForDebug.push_back(mus->Name + "_InsertionY");
        p_all.push_back(mus->InsertionPointGlobal.z); inputParamsForDebug.push_back(mus->InsertionPointGlobal.z);  inputParamDescriptionsForDebug.push_back(mus->Name + "_InsertionZ");
        
        // --- CONSTRAINT BOUNDS (lbg/ubg) ---
        // Für jeden inneren Knoten:
        for (int k = 0; k < num_inner; ++k) {
            
            // 1. Nicht-Eindringen für jedes Mesh: h >= 0
            for (int j = 0; j < num_wrap; ++j) {
                lbg_all.push_back(0.0);    // h >= 0
                ubg_all.push_back(inf);
            }
            
            // 2. Eine summierte Komplementarität: sum(h*eta) = 0
            lbg_all.push_back(0.0);
            ubg_all.push_back(0.0);
            
            // 3. Euler-Lagrange Gleichung (3 Komponenten: x, y, z) = 0
            for (int d = 0; d < 3; ++d) { 
                lbg_all.push_back(0.0); 
                ubg_all.push_back(0.0); 
            }
        }

        // --- VARIABLE BOUNDS (lbx/ubx) ---
        // Gamma (Punkte) sind unbeschränkt
        for (int i = 0; i < num_inner * 3; ++i) {
            lbx_all.push_back(-inf); 
            ubx_all.push_back(inf);
        }
        // Eta (Kräfte) >= 0
        for (int i = 0; i < num_etas_total; ++i) {
            lbx_all.push_back(0.0); 
            ubx_all.push_back(inf);
        }
    }

    allParameterInputsAllSteps.push_back(inputParamsForDebug);
    allParameterInputDescriptionsAllSteps.push_back(inputParamDescriptionsForDebug);

    // --- SOLVER CALL ---
    DMDict arg = {{"x0", x0_all}, {"p", p_all}, {"lbg", lbg_all}, {"ubg", ubg_all}, {"lbx", lbx_all}, {"ubx", ubx_all}};
    DMDict res = solver(arg);

    // --- EXTRACT RESULTS ---
    std::vector<double> res_x = std::vector<double>(res.at("x"));
    auto solverInfo = solver.stats();
    std::string status = solverInfo.at("return_status");

    if (status != "Solve_Succeeded") {
        qDebug() << "    Warning: Multi-Muscle Solver Status:" << QString::fromStdString(status);
    }
    int convSteps = int(solverInfo.at("iter_count"));
    qDebug() << "    Solver converged in" << convSteps << "iterations";
    SolverConvergenceMessages.push_back(QString::fromStdString(status).toStdString());
    SolverConvergenceSteps.push_back(convSteps);
    
    int current_x_offset = 0;
    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;
        int num_wrap = mus->meshPtrs.size();
        int num_etas_mus = num_inner * num_wrap;

        // Gamma (Punkte) zurückschreiben
        for (int i = 0; i < num_inner; ++i) {
            int mIdx = i + 1;
            MWMath::Point3D p;
            p.x = res_x[current_x_offset + i * 3 + 0];
            p.y = res_x[current_x_offset + i * 3 + 1];
            p.z = res_x[current_x_offset + i * 3 + 2];

            mus->MusclePointsGlobal[mIdx] = p;
            mus->MNodes[mIdx].PositionGlobal = p;
            mus->MNodes[mIdx].updateLocalFrame();
        }

        // Etas für Warmstart speichern
        int eta_start_idx = current_x_offset + (num_inner * 3);
        mus->lastEtas.assign(
            res_x.begin() + eta_start_idx, 
            res_x.begin() + eta_start_idx + num_etas_mus
        );

        // ETAS IN LISTE SCHREIBEN (für Export)
        int num_eta_per_node = num_wrap;
        for (int k = 0; k < num_inner; ++k) {
            int mIdx = k + 1;
            
            std::vector<double> currentStepEtas;
            
            for (int e = 0; e < num_eta_per_node; ++e) {
                int globalEtaIndex = k * num_eta_per_node + e;
                if (globalEtaIndex < mus->lastEtas.size()) {
                    currentStepEtas.push_back(mus->lastEtas[globalEtaIndex]);
                } else {
                    currentStepEtas.push_back(0.0);
                }
            }

            mus->MNodes[mIdx].MNodeEtaSteps.push_back(currentStepEtas);
            mus->MNodes[mIdx].lastEtas = mus->lastEtas; 
        }

        for (auto& node : mus->MNodes) {
            node.lastEtas = mus->lastEtas;
        }

        current_x_offset += (num_inner * 3 + num_etas_mus);
    }
}