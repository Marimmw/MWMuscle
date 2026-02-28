#include "casadiSystem.h"

CasadiSystem::CasadiSystem(std::vector<SSMuscle*> muscles, int objType, std::string version, std::string parametrizationType, bool bUseCasGradient, bool bSumPhiEta, bool bUseWarmstartEtas, bool bDebug)
    : m_muscles(muscles), objType(objType), Version(version), ParametrizationType(parametrizationType),
        bUseOwnGradient(bUseCasGradient), bSumPhiEta(bSumPhiEta), bUseWarmstartEtas(bUseWarmstartEtas), bDebug(bDebug)
{
    CasadiSystemName = "CasSys_" + (m_muscles.empty() ? "Empty" : m_muscles[0]->Name);

   //setupCasadi();
   if (bSumPhiEta) {
    setupCasadiSum();
   }
   else {
    setupCasadi();
   }
    
}

void CasadiSystem::solveStepX(){
    if (bSumPhiEta) {
        solveStepSum();
    }
    else{solveStep();}
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
    std::string filename = "../examples/results/solver_log_" + CasadiSystemName + ".txt";
    opts["ipopt.output_file"] = filename; 
    opts["ipopt.file_print_level"] = 5;

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

    if (status != "Solve_Succeeded") {
        if (bDebug) qDebug() << "    Warning: Multi-Muscle Solver Status:" << QString::fromStdString(status);
    }
    int convSteps = int(solverInfo.at("iter_count"));
    if (bDebug) qDebug() << "    Solver converged in" << convSteps << "iterations";
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

    if (status != "Solve_Succeeded") {
        if (bDebug) qDebug() << "    Warnung: Multi-Muscle Solver Status:" << QString::fromStdString(status);
    }
    int convSteps = int(solverInfo.at("iter_count"));
    if (bDebug) qDebug() << "    Solver converged in" << convSteps << "iterations";
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


