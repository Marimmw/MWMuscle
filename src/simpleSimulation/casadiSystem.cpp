#include "casadiSystem.h"

#include <filesystem>


CasadiSystem::CasadiSystem(std::vector<SSMuscle*> muscles, int objType, std::string version, std::string parametrizationType, bool bUseCasGradient, bool bSumPhiEta, bool bUseWarmstartEtas, bool bDebug, bool bWriteFiles)
    : m_muscles(muscles), objType(objType), Version(version), ParametrizationType(parametrizationType),
        bUseOwnGradient(bUseCasGradient), bSumPhiEta(bSumPhiEta), bUseWarmstartEtas(bUseWarmstartEtas), bDebug(bDebug), bWriteFiles(bWriteFiles)
{
    CasadiSystemName = "CasSys_" + (m_muscles.empty() ? "Empty" : m_muscles[0]->Name);

    //setupCasadi();
    /* if (bSumPhiEta) {
        setupCasadiSum();
    }
    else {
        setupCasadi();
    } */

    setupCasadiViaSum_paramStudy();
    
    
    
}

void CasadiSystem::solveStepX(){
    /* if (bSumPhiEta) {
        solveStepSum();
    }
    else{solveStep();} */
    solveStepViaSum_paramStudy();
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
            d_vector = MX::vertcat({d_vector, sum1(sq(P_orig - v_pos))});
            d_vector = MX::vertcat({d_vector, sum1(sq(P_ins - v_pos))});

            MX D_v = -casadi::MX::logsumexp(-Alpha * d_vector) / Alpha;
            MX h_v = r_sq - D_v; // hier eventuell linear statt quadratisch abziehen -> (quadratuisch - linear funkioniert aber irgendiwe besser?)
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
            d_vector = MX::vertcat({d_vector, sum1(sq(P_orig - v_pos))});
            d_vector = MX::vertcat({d_vector, sum1(sq(P_ins - v_pos))});

            MX D_v = -casadi::MX::logsumexp(-Alpha * d_vector) / Alpha;
            MX h_v = r_sq - D_v; // hier eventuell linear statt quadratisch abziehen -> (quadratuisch - linear funkioniert aber irgendiwe besser?)
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



// VIA POINTS SUM
void CasadiSystem::solveStepViaSum_paramStudy()
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

void CasadiSystem::setupCasadiViaSum_paramStudy()
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
            d_vector = MX::vertcat({d_vector, sum1(sq(P_orig - v_pos))});
            d_vector = MX::vertcat({d_vector, sum1(sq(P_ins - v_pos))});

            MX D_v = -casadi::MX::logsumexp(-Alpha * d_vector) / Alpha;
            MX h_v = r_tol - D_v; // hier eventuell linear statt quadratisch abziehen -> (quadratuisch - linear funkioniert aber irgendiwe besser?)
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


