#include "casadiSystem.h"

CasadiSystem::CasadiSystem(std::vector<SSMuscle*> muscles, int objType, std::string version, std::string parametrizationType, bool bUseCasGradient)
    : m_muscles(muscles), objType(objType), Version(version), ParametrizationType(parametrizationType),
        bUseOwnGradient(bUseCasGradient)
{
    CasadiSystemName = "CasSys_" + (m_muscles.empty() ? "Empty" : m_muscles[0]->Name);

    /* if (parametrizationType == "global"){
        qDebug() << "SET UP CASADI SYSTEM GLOBAL";
        if (version == "VPP"){
            setupCasadiVPP();
        }
        else if (version == "VPPenalty"){
            setupCasadiVPPenalty();
        }
        else if (version == "Exp"){
            setupCasadiEXP();
        }
        else{
            setupCasadi();
        }
    }
    else if (parametrizationType == "local"){
        qDebug() << "SET UP CASADI SYSTEM LOCAL";
        setupCasadiLocal();
    } */
   setupCasadi();
    
}

void CasadiSystem::solveStepX(){
    /* if (ParametrizationType == "global"){
        if (Version == "VPP"){
            solveStep();
        }
        else if (Version == "VPPenalty"){
            solveStepPenalty();
        }
        else if (Version == "Exp"){
            solveStepEXP();
        }
        else{
            solveStep();
        }
    }
    else if (ParametrizationType == "local"){
        qDebug() << "SOLVE STEP LOCAL";
        solveStepLocal();
    } */
    solveStep();
}


// LOCAL
void CasadiSystem::setupCasadiLocal()
{
    using namespace casadi;
    MX all_x = MX::vertcat({});
    MX all_g = MX::vertcat({});
    MX all_p = MX::vertcat({});
    MX obj = 1;
    M = static_cast<int>(m_muscles.size());

    for (size_t m = 0; m < m_muscles.size(); ++m) {
        SSMuscle* mus = m_muscles[m];
        int num_inner = mus->MNodes.size() - 2;
        int num_wrap = mus->meshPtrs.size();

        // 1. Via-Points identifizieren
        std::vector<int> via_indices;
        for (int j = 0; j < num_wrap; ++j) {
            if (mus->meshPtrs[j]->bIsViaPoint) via_indices.push_back(j);
        }
        int num_via = via_indices.size();
        int num_eta_per_node = num_wrap - num_via;

        // --- VARIABLEN (x) ---
        // Jetzt: LOKALE Punkte (3*n) + Etas
        // x_mus[0..2] ist jetzt l_k (lokal), nicht mehr g_k (global)!
        MX x_mus = MX::sym("x_" + std::to_string(m), num_inner * 3 + num_inner * num_eta_per_node);
        all_x = MX::vertcat({all_x, x_mus});

        // --- PARAMETER (p) ---
        // A. Meshes (Hindernisse): [Pos(3), Rot(9)] -> 12 pro Mesh
        // B. Origin & Insertion: [Pos(3)] -> 6 total (sind fix global oder werden global geupdatet)
        // C. NEU: Parent-Transformationen für jeden inneren Knoten! 
        //    Wir brauchen für jeden Knoten k wissen, an welchem Body er hängt, um l_k -> g_k zu rechnen.
        //    Format pro Knoten: [ParentPos(3), ParentRot(9)] -> 12 Parameter pro Knoten
        
        int params_meshes = 12 * num_wrap;
        int params_endpoints = 6;
        int params_node_parents = 12 * num_inner;

        MX p_mus = MX::sym("p_" + std::to_string(m), params_meshes + params_endpoints + params_node_parents);
        all_p = MX::vertcat({all_p, p_mus});

        // Slices für Parameter
        MX P_meshes = p_mus(Slice(0, params_meshes));
        MX P_orig   = p_mus(Slice(params_meshes, params_meshes + 3));
        MX P_ins    = p_mus(Slice(params_meshes + 3, params_meshes + 6));
        MX P_nodes  = p_mus(Slice(params_meshes + 6, p_mus.size1())); // Die neuen Parent-Trafos

        int eta_offset = num_inner * 3;

        // Container für Minimum-Abstände (Via-Points)
        std::vector<MX> min_dists_sq(num_via);

        // --- GLEICHUNGEN (g) ---
        for (int k = 0; k < num_inner; ++k) {
            
            // 1. LOKALE Koordinate aus Variable holen
            MX l_k = x_mus(Slice(k * 3, (k + 1) * 3));

            // 2. GLOBALE Koordinate berechnen (Forward Kinematics on-the-fly)
            // Parameter für diesen Knoten holen
            MX p_node_k = P_nodes(Slice(k * 12, (k + 1) * 12));
            MX parent_pos = p_node_k(Slice(0, 3));
            
            // Rotation rekonstruieren (Row-Major oder Col-Major? Konsistent bleiben!)
            // Hier Annahme: p enthält R00, R01, R02, R10... (Zeilenweise)
            MX parent_rot = MX::zeros(3, 3);
            parent_rot(0,0) = p_node_k(3); parent_rot(0,1) = p_node_k(4); parent_rot(0,2) = p_node_k(5);
            parent_rot(1,0) = p_node_k(6); parent_rot(1,1) = p_node_k(7); parent_rot(1,2) = p_node_k(8);
            parent_rot(2,0) = p_node_k(9); parent_rot(2,1) = p_node_k(10);parent_rot(2,2) = p_node_k(11);

            // Transformation: g_k = ParentPos + ParentRot * l_k
            MX g_k = parent_pos + mtimes(parent_rot, l_k);


            // 3. Nachbarn (g_prev, g_next) bestimmen
            MX g_prev;
            if (k == 0) {
                g_prev = P_orig; // Origin ist global
            } else {
                // Vorherigen Knoten auch erst globalisieren!
                MX l_prev = x_mus(Slice((k - 1) * 3, k * 3));
                
                MX p_node_prev = P_nodes(Slice((k - 1) * 12, k * 12));
                MX prev_parent_pos = p_node_prev(Slice(0, 3));
                MX prev_parent_rot = MX::zeros(3, 3);
                // ... Füllen (Copy-Paste oder Helper Funktion in CasADi wäre gut, hier inline):
                prev_parent_rot(0,0) = p_node_prev(3); prev_parent_rot(0,1) = p_node_prev(4); prev_parent_rot(0,2) = p_node_prev(5);
                prev_parent_rot(1,0) = p_node_prev(6); prev_parent_rot(1,1) = p_node_prev(7); prev_parent_rot(1,2) = p_node_prev(8);
                prev_parent_rot(2,0) = p_node_prev(9); prev_parent_rot(2,1) = p_node_prev(10);prev_parent_rot(2,2) = p_node_prev(11);

                g_prev = prev_parent_pos + mtimes(prev_parent_rot, l_prev);
            }

            MX g_next;
            if (k == num_inner - 1) {
                g_next = P_ins; // Insertion ist global
            } else {
                // Nächsten Knoten globalisieren
                MX l_next = x_mus(Slice((k + 1) * 3, (k + 2) * 3));
                
                MX p_node_next = P_nodes(Slice((k + 1) * 12, (k + 2) * 12));
                MX next_parent_pos = p_node_next(Slice(0, 3));
                MX next_parent_rot = MX::zeros(3, 3);
                next_parent_rot(0,0) = p_node_next(3); next_parent_rot(0,1) = p_node_next(4); next_parent_rot(0,2) = p_node_next(5);
                next_parent_rot(1,0) = p_node_next(6); next_parent_rot(1,1) = p_node_next(7); next_parent_rot(1,2) = p_node_next(8);
                next_parent_rot(2,0) = p_node_next(9); next_parent_rot(2,1) = p_node_next(10);next_parent_rot(2,2) = p_node_next(11);

                g_next = next_parent_pos + mtimes(next_parent_rot, l_next);
            }

            // --- Ab hier ist fast alles gleich, da wir g_k, g_prev, g_next haben ---

            MX total_contact_force = 0;
            int current_eta_idx = 0;
            int current_via_idx = 0;

            for (int j = 0; j < num_wrap; ++j) {
                MX q_j = P_meshes(Slice(j * 12, (j + 1) * 12)); // Mesh Parameter

                if (mus->meshPtrs[j]->bIsViaPoint) {
                    MX p_via = q_j(Slice(0, 3));
                    MX d2 = sumsqr(g_k - p_via);
                    if (k == 0) min_dists_sq[current_via_idx] = d2;
                    else min_dists_sq[current_via_idx] = fmin(min_dists_sq[current_via_idx], d2);
                    current_via_idx++;
                } else {
                    MX eta_kj = x_mus(eta_offset + k * num_eta_per_node + current_eta_idx);
                    
                    // Distanz & Gradient im GLOBALEN Raum
                    MX h = mus->meshPtrs[j]->constraintDistanceLocal(g_k, q_j);
                    MX grad_h_global = mus->meshPtrs[j]->constraintJacobianLocal(g_k, q_j);
                    
                    /* MX dummy_pos = MX::sym("dummy_pos", 3);
                    MX dummy_h = mus->meshPtrs[j]->constraintDistanceLocal(dummy_pos, q_j);
                    MX dummy_grad = MX::gradient(dummy_h, dummy_pos);
                    MX grad_h_global = MX::substitute(dummy_grad, dummy_pos, g_k); */
                    
                    total_contact_force += eta_kj * grad_h_global;
                    
                    all_g = MX::vertcat({all_g, h * eta_kj, h, eta_kj});
                    current_eta_idx++;
                }
            }

            // Euler-Lagrange (Kraftbilanz im GLOBALEN Raum)
            // Das ist physikalisch korrekt. Der Solver variiert l_k, das ändert g_k, das muss die Bilanz erfüllen.
            MX eq_el = (g_next - 2.0 * g_k + g_prev) + total_contact_force;
            all_g = MX::vertcat({all_g, eq_el});

            // Objective (Länge minimieren)
            if (objType == 1) obj = obj + norm_2(g_next - g_k);
            else if (objType == 2) { /*...*/ }
        }

        // Global Via-Point Constraints
        for (int v = 0; v < num_via; ++v) {
            double tol_sq = pow(mus->meshPtrs[via_indices[v]]->MViaPointTolerance, 2);
            all_g = MX::vertcat({all_g, min_dists_sq[v] - tol_sq});
        }
    }

    opts["print_time"] = false;
    opts["record_time"] = false;
    opts["ipopt.print_level"] = 0;
    opts["ipopt.max_iter"] = maxIterations;
    opts["ipopt.tol"] = maxTol;
    opts["ipopt.linear_solver"] = "mumps";

    std::string filename = "../examples/solverLogLocal_" + CasadiSystemName + ".txt";
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

void CasadiSystem::solveStepLocal() {
    using namespace casadi;
    std::vector<double> x0_all, p_all, lbg_all, ubg_all, lbx_all, ubx_all;

    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;
        int num_wrap = mus->meshPtrs.size();
        
        // Zähler
        int num_via = 0;
        for (auto* m : mus->meshPtrs) if (m->bIsViaPoint) num_via++;
        int num_eta_per_node = num_wrap - num_via;
        int num_etas_total = num_inner * num_eta_per_node;

        // --- 1. INITIAL GUESS (x0) - LOKAL! ---
        qDebug() << "INITIAL GEUESSES for solver";
        for (size_t k = 1; k < mus->MNodes.size() - 1; ++k) {
             SSTissue* parent = mus->MNodes[k].parentMesh->Parent.get();
            if (parent) {
                x0_all.push_back(mus->MNodes[k].PositionLocal.x);
                x0_all.push_back(mus->MNodes[k].PositionLocal.y);
                x0_all.push_back(mus->MNodes[k].PositionLocal.z);
                /* qDebug() << "Node" << k << " Local Pos: (" << mus->MNodes[k].PositionLocal.print().c_str() << ")"
                         << " Parent: " << parent->Name.c_str(); */
            } else {
                qDebug() << "Node" << k << "WARNING!!!!!!!!!!!!!!!!!! -> no parent Body";
            }  
           
        }
        
        // Etas (bleiben gleich)
        if (mus->lastEtas.size() == (size_t)num_etas_total && bUseWarmstartEtas) {
            x0_all.insert(x0_all.end(), mus->lastEtas.begin(), mus->lastEtas.end());
        } else {
            for (int e = 0; e < num_etas_total; ++e) x0_all.push_back(0.001);
        }

        // --- 2. PARAMETER (p) ---
        // A. Meshes
        for (auto* mesh : mus->meshPtrs) {
            p_all.push_back(mesh->PositionGlobal.x); p_all.push_back(mesh->PositionGlobal.y); p_all.push_back(mesh->PositionGlobal.z);
            for(int r=0;r<3;r++) for(int c=0;c<3;c++) p_all.push_back(mesh->OrientationGlobal.m[r][c]);
        }
        // B. Endpoints
        p_all.push_back(mus->OriginPointGlobal.x); p_all.push_back(mus->OriginPointGlobal.y); p_all.push_back(mus->OriginPointGlobal.z);
        p_all.push_back(mus->InsertionPointGlobal.x); p_all.push_back(mus->InsertionPointGlobal.y); p_all.push_back(mus->InsertionPointGlobal.z);

        // C. NODE PARENTS (NEU!)
        for (size_t k = 1; k < mus->MNodes.size() - 1; ++k) {
            SSTissue* refBody = mus->MNodes[k].parentMesh->Parent.get();
            if (refBody) {
                // refBody Pos
                //qDebug() << "Node" << k << " Parent for Casadi " << refBody->Name.c_str();
                //qDebug() << "   Pos: (" << refBody->PositionGlobal.x << refBody->PositionGlobal.y << refBody->PositionGlobal.z << ") +  ROT * ("<< mus->MNodes[k].PositionLocal.x << mus->MNodes[k].PositionLocal.y << mus->MNodes[k].PositionLocal.z<< ")";
                p_all.push_back(refBody->PositionGlobal.x);
                p_all.push_back(refBody->PositionGlobal.y);
                p_all.push_back(refBody->PositionGlobal.z);
                // Parent Rot (Zeilenweise für Casadi Setup)
                MWMath::RotMatrix3x3 R = refBody->OrientationGlobal;
                for(int r=0;r<3;r++) for(int c=0;c<3;c++) p_all.push_back(R.m[r][c]);

                /* qDebug() << "Node" << k << " Parent for Casadi " << refBody->Name.c_str();
                qDebug() << "   Pos:" << refBody->PositionGlobal.print().c_str(); 
                qDebug() << "   Rot:" << R.print().c_str(); */
                
            } /*else {
                // Identity Parent
                p_all.push_back(0); p_all.push_back(0); p_all.push_back(0); // Pos
                p_all.push_back(1); p_all.push_back(0); p_all.push_back(0); // Rot R0
                p_all.push_back(0); p_all.push_back(1); p_all.push_back(0); // Rot R1
                p_all.push_back(0); p_all.push_back(0); p_all.push_back(1); // Rot R2
            } */
        }

        // --- 3. BOUNDS (lbg, ubg, lbx, ubx) ---
        // (Constraints-Logik bleibt gleich, da g(x) immer noch Kräfte/Abstände sind)
        for (int i = 0; i < num_inner; ++i) {
            for (int j = 0; j < num_wrap; ++j) {
                if (!mus->meshPtrs[j]->bIsViaPoint) {
                    lbg_all.push_back(0.0); ubg_all.push_back(0.0);    // Komp
                    lbg_all.push_back(0.0); ubg_all.push_back(1e10);   // Dist
                    lbg_all.push_back(0.0); ubg_all.push_back(1e10);   // Kraft
                }
            }
            // EL = 0
            for (int d = 0; d < 3; ++d) { lbg_all.push_back(-ELTolerance); ubg_all.push_back(ELTolerance); }
        }
        // ViaPoints
        for (int v = 0; v < num_via; ++v) { lbg_all.push_back(-inf); ubg_all.push_back(0.0); }

        // Variablen Bounds
        for (int i = 0; i < num_inner * 3; ++i) { lbx_all.push_back(-inf); ubx_all.push_back(inf); } // Lokale Koordinaten sind frei
        for (int i = 0; i < num_etas_total; ++i) { lbx_all.push_back(0.0); ubx_all.push_back(inf); }
    }

    // --- SOLVER ---
    DMDict arg = {{"x0", x0_all}, {"p", p_all}, {"lbg", lbg_all}, {"ubg", ubg_all}, {"lbx", lbx_all}, {"ubx", ubx_all}};
    DMDict res = solver(arg);
    auto solverInfo = solver.stats();
    std::string status = solverInfo.at("return_status");

    if (status != "Solve_Succeeded") {
        qDebug() << "    Warnung: Multi-Muscle Solver Status:" << QString::fromStdString(status);
    }
    int convSteps = int(solverInfo.at("iter_count"));
    qDebug() << "    Solver converged in" << convSteps << "iterations";
    SolverConvergenceMessages.push_back(QString::fromStdString(status).toStdString());
    SolverConvergenceSteps.push_back(convSteps);
    
    // --- RESULTATE (Rückrechnung) ---
    std::vector<double> res_x = std::vector<double>(res.at("x"));
    
    int current_x_offset = 0;
    // Achtung: Wir müssen auch den Offset im p_all Vektor tracken, um die Parent-Trafos wiederzufinden!
    // Oder wir nutzen direkt die MNodes-Parent-Pointer, da die sich während dem Step nicht geändert haben.
    
    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;
        int num_wrap = mus->meshPtrs.size();
        int num_via = 0; for (auto* m : mus->meshPtrs) if (m->bIsViaPoint) num_via++;
        int num_etas_mus = num_inner * (num_wrap - num_via);

        for (int i = 0; i < num_inner; ++i) {
            int mIdx = i + 1;
            
            // Lokale Koordinaten aus Solver
            double lx = res_x[current_x_offset + i * 3 + 0];
            double ly = res_x[current_x_offset + i * 3 + 1];
            double lz = res_x[current_x_offset + i * 3 + 2];
            MWMath::Point3D l_p(lx, ly, lz);

            // Speichern als local
            mus->MNodes[mIdx].PositionLocal = l_p;
            // qDebug() << "Muscle" << mus->Name.c_str() << " Node" << mIdx << " Local Pos Solver: (" << l_p.print().c_str() << ")" << "parent:" << mus->MNodes[mIdx].parentMesh->Parent->Name.c_str();

            // In Global umrechnen für Visualisierung & Nächsten Schritt
            SSTissue* parent = mus->MNodes[mIdx].parentMesh->Parent.get();
            if (parent) {
                mus->MNodes[mIdx].PositionGlobal = parent->PositionGlobal + parent->OrientationGlobal.transform(l_p);
            } else {
                mus->MNodes[mIdx].PositionGlobal = l_p;
            }
            mus->MusclePointsGlobal[mIdx] = mus->MNodes[mIdx].PositionGlobal;
        }

        // Etas speichern
        int eta_start = current_x_offset + num_inner * 3;
        if(num_etas_mus > 0) 
            mus->lastEtas.assign(res_x.begin() + eta_start, res_x.begin() + eta_start + num_etas_mus);
        
        current_x_offset += (num_inner * 3 + num_etas_mus);
    }
}



void CasadiSystem::setupCasadiVPP(){
    qDebug() << "setupCasadiVPP";
    using namespace casadi;
    MX all_x = MX::vertcat({});
    MX all_g = MX::vertcat({});
    MX all_p = MX::vertcat({});
    MX obj = 1;
    M = static_cast<int>(m_muscles.size());

    for (size_t m = 0; m < m_muscles.size(); ++m) {
        SSMuscle* mus = m_muscles[m];
        int num_inner = mus->MNodes.size() - 2;
        int num_wrap = mus->meshPtrs.size();

        // 1. Via-Points identifizieren
        std::vector<int> via_indices;
        for (int j = 0; j < num_wrap; ++j) {
            if (mus->meshPtrs[j]->bIsViaPoint) via_indices.push_back(j);
        }
        int num_via = via_indices.size();
        int num_eta_per_node = num_wrap - num_via;

        // --- VARIABLEN: [Punkte (3*n), Etas (n * (Meshes-Vias))] ---
        MX x_mus = MX::sym("x_" + std::to_string(m), num_inner * 3 + num_inner * num_eta_per_node);
        all_x = MX::vertcat({all_x, x_mus});

        // --- PARAMETER: [Meshes (12*m), Origin (3), Insertion (3)] ---
        MX p_mus = MX::sym("p_" + std::to_string(m), 12 * num_wrap + 6);
        all_p = MX::vertcat({all_p, p_mus});

        MX P_orig = p_mus(Slice(12 * num_wrap, 12 * num_wrap + 3));
        MX P_ins  = p_mus(Slice(12 * num_wrap + 3, 12 * num_wrap + 6));
        int eta_offset = num_inner * 3;

        // Container für Minimum-Abstände (einer pro Via-Point über alle Knoten k)
        std::vector<MX> min_segment_dists_sq(num_via, 0.0);
        /* auto distSqSegmentPoint = [](MX A, MX B, MX P) { // HELPER
            MX ab = B - A;
            MX ap = P - A;
            // Projektion t auf die Linie (dot product)
            // t = (ap . ab) / (ab . ab)
            MX t = mtimes(ap.T(), ab) / mtimes(ab.T(), ab);
            // Clamp t zwischen 0.0 und 1.0 (wir wollen auf dem Segment bleiben!)
            // fmax(0, fmin(1, t))
            t = fmax(0.0, fmin(1.0, t));
            // Der nächste Punkt auf dem Segment
            MX closest = A + t * ab;
            // Quadratischer Abstand
            return sumsqr(closest - P);
        }; */
        auto distSqSegmentPoint = [](MX A, MX B, MX P) { 
            // Berechne den Mittelpunkt des Segments
            MX M = (A + B) / 2.0;
            // Gib den quadratischen Abstand zum ViaPoint zurück
            return sumsqr(M - P);
        };
        int current_via_idx = 0;

        // --- GLEICHUNGEN (g) ---
        for (int k = 0; k < num_inner; ++k) {
            MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
            MX g_prev = (k == 0) ? P_orig : x_mus(Slice((k - 1) * 3, k * 3));
            MX g_next = (k == num_inner - 1) ? P_ins : x_mus(Slice((k + 1) * 3, (k + 2) * 3));

            MX total_contact_force = 0;
            int current_eta_idx = 0;
            int current_via_idx = 0;

            for (int j = 0; j < num_wrap; ++j) {
                MX q_j = p_mus(Slice(j * 12, (j + 1) * 12));

                if (mus->meshPtrs[j]->bIsViaPoint) {
                    // Via-Point Logik: Minimum-Suche über alle Knoten k
                    MX q_j = p_mus(Slice(j * 12, (j + 1) * 12));
                    MX p_via = q_j(Slice(0, 3));
                    
                    // Berechne Abstand vom Segment (Prev -> Curr) zum ViaPoint
                    MX seg_dist_sq = distSqSegmentPoint(g_prev, g_k, p_via);
                    
                    // Sonderfall: Letztes Segment (NodeN -> Insertion)
                    // Das müssen wir nur einmal ganz am Ende (im letzten Loop) machen
                    MX last_seg_dist_sq = 1e10; // Default groß
                    if (k == num_inner - 1) {
                        last_seg_dist_sq = distSqSegmentPoint(g_k, P_ins, p_via);
                    }

                    if (k == 0) {
                        min_segment_dists_sq[current_via_idx] = seg_dist_sq;
                    } else {
                        // Minimum updaten
                        min_segment_dists_sq[current_via_idx] = fmin(min_segment_dists_sq[current_via_idx], seg_dist_sq);
                    }
                    
                    // Am Ende auch das allerletzte Segment prüfen
                    if (k == num_inner - 1) {
                        min_segment_dists_sq[current_via_idx] = fmin(min_segment_dists_sq[current_via_idx], last_seg_dist_sq);
                    }
                    current_via_idx++;
                } else {
                    // Normale Signorini-Logik (Hindernisse)
                    MX eta_kj = x_mus(eta_offset + k * num_eta_per_node + current_eta_idx);
                    MX h = mus->meshPtrs[j]->constraintDistance(g_k, q_j);

                    // MX grad_h = mus->meshPtrs[j]->constraintJacobian(g_k, q_j);
                    MX grad_h_full = MX::gradient(h, x_mus);
                    MX grad_h = grad_h_full(Slice(k * 3, (k + 1) * 3));

                    total_contact_force += eta_kj * grad_h;
                    
                    // Signorini/Complementarity
                    all_g = MX::vertcat({all_g, h * eta_kj, h, eta_kj}); 
                    current_eta_idx++;
                }
            }

            // Euler-Lagrange Kraftbilanz
            MX eq_el = (g_next - 2.0 * g_k + g_prev) + total_contact_force;
            all_g = MX::vertcat({all_g, eq_el});
            
            // Objective
            if (objType == 1) {
                obj = obj + norm_2(g_next - g_k);
            }
            else if (objType == 2) {
                MX prev_seg_sq = sumsqr(g_k - g_prev);
                MX current_seg_sq = sumsqr(g_next - g_k);
                double weight_equi = 10.0; 
                obj += weight_equi * sumsqr(prev_seg_sq - current_seg_sq);
            }
        }

        // Global Via-Point Constraints (one Constraint per Via-Point for whole muscle)
        for (int v = 0; v < num_via; ++v) {
            // Hole Toleranz aus dem Mesh Pointer
            double tol_sq = pow(mus->meshPtrs[via_indices[v]]->MViaPointTolerance, 2);
            
            // KORREKTUR 1: Richtiger Variablenname 'min_segment_dists_sq'
            MX dist_constraint = min_segment_dists_sq[v] - tol_sq; 
            
            // KORREKTUR 2: Expliziter Vektor für vertcat
            all_g = MX::vertcat({all_g, dist_constraint}); 
        }
    }

    opts["print_time"] = false;
    opts["record_time"] = false;
    opts["ipopt.print_level"] = 0;
    opts["ipopt.max_iter"] = maxIterations;
    opts["ipopt.tol"] = maxTol;
    opts["ipopt.linear_solver"] = "mumps";

    std::string filename = "../examples/solver_log_" + CasadiSystemName + ".txt";
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


void CasadiSystem::setupCasadiVPPenalty(){
    qDebug() << "setupCasadiVPP (Penalty Method)";
    using namespace casadi;
    MX all_x = MX::vertcat({});
    MX all_g = MX::vertcat({});
    MX all_p = MX::vertcat({});
    MX obj = 1;
    M = static_cast<int>(m_muscles.size());

    for (size_t m = 0; m < m_muscles.size(); ++m) {
        SSMuscle* mus = m_muscles[m];
        int num_inner = mus->MNodes.size() - 2;
        int num_wrap = mus->meshPtrs.size();

        // 1. Via-Points identifizieren
        std::vector<int> via_indices;
        for (int j = 0; j < num_wrap; ++j) {
            if (mus->meshPtrs[j]->bIsViaPoint) via_indices.push_back(j);
        }
        int num_via = via_indices.size();
        int num_eta_per_node = num_wrap - num_via;

        // --- VARIABLEN & PARAMETER ---
        MX x_mus = MX::sym("x_" + std::to_string(m), num_inner * 3 + num_inner * num_eta_per_node);
        all_x = MX::vertcat({all_x, x_mus});

        MX p_mus = MX::sym("p_" + std::to_string(m), 12 * num_wrap + 6);
        all_p = MX::vertcat({all_p, p_mus});

        MX P_orig = p_mus(Slice(12 * num_wrap, 12 * num_wrap + 3));
        MX P_ins  = p_mus(Slice(12 * num_wrap + 3, 12 * num_wrap + 6));
        int eta_offset = num_inner * 3;

        // Container für Minimum-Abstände
        std::vector<MX> min_segment_dists_sq(num_via, 0.0);
        
        // --- HELPER: ECHTE SEGMENT PROJEKTION (Besser für Penalty!) ---
        auto distSqSegmentPoint = [](MX A, MX B, MX P) { 
            MX ab = B - A;
            MX ap = P - A;
            // Projektion t auf die Linie
            MX t = mtimes(ap.T(), ab) / mtimes(ab.T(), ab);
            // Clamp t zwischen 0.0 und 1.0
            t = fmax(0.0, fmin(1.0, t));
            // Der nächste Punkt auf dem Segment
            MX closest = A + t * ab;
            return sumsqr(closest - P);
        };
        // -------------------------------------------------------------

        int current_via_idx = 0;

        // --- SCHLEIFE ÜBER KNOTEN ---
        for (int k = 0; k < num_inner; ++k) {
            MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
            MX g_prev = (k == 0) ? P_orig : x_mus(Slice((k - 1) * 3, k * 3));
            MX g_next = (k == num_inner - 1) ? P_ins : x_mus(Slice((k + 1) * 3, (k + 2) * 3));

            MX total_contact_force = 0;
            int current_eta_idx = 0;
            int current_via_idx = 0;

            for (int j = 0; j < num_wrap; ++j) {
                MX q_j = p_mus(Slice(j * 12, (j + 1) * 12));

                if (mus->meshPtrs[j]->bIsViaPoint) {
                    // --- VIA POINT LOGIK ---
                    MX p_via = q_j(Slice(0, 3));
                    MX seg_dist_sq = distSqSegmentPoint(g_prev, g_k, p_via);
                    
                    MX last_seg_dist_sq = 1e10; 
                    if (k == num_inner - 1) {
                        last_seg_dist_sq = distSqSegmentPoint(g_k, P_ins, p_via);
                    }

                    if (k == 0) {
                        min_segment_dists_sq[current_via_idx] = seg_dist_sq;
                    } else {
                        min_segment_dists_sq[current_via_idx] = fmin(min_segment_dists_sq[current_via_idx], seg_dist_sq);
                    }
                    
                    if (k == num_inner - 1) {
                        min_segment_dists_sq[current_via_idx] = fmin(min_segment_dists_sq[current_via_idx], last_seg_dist_sq);
                    }
                    current_via_idx++;

                } else {
                    // --- OBSTACLE LOGIK ---
                    MX eta_kj = x_mus(eta_offset + k * num_eta_per_node + current_eta_idx);
                    MX h = mus->meshPtrs[j]->constraintDistance(g_k, q_j);
                    MX grad_h_full = MX::gradient(h, x_mus);
                    MX grad_h = grad_h_full(Slice(k * 3, (k + 1) * 3));

                    total_contact_force += eta_kj * grad_h;
                    all_g = MX::vertcat({all_g, h * eta_kj, h, eta_kj}); 
                    current_eta_idx++;
                }
            }

            // Euler-Lagrange
            MX eq_el = (g_next - 2.0 * g_k + g_prev) + total_contact_force;
            all_g = MX::vertcat({all_g, eq_el});
            
            // Längen-Minimierung
            if (objType == 1) obj += sumsqr(g_next - g_k);
            else if (objType == 2) { /*...*/ }
        }

        // --- PENALTY STATT CONSTRAINT ---
        // Wir fügen KEINE Zeilen zu all_g hinzu!
        // Stattdessen erhöhen wir obj.
        double penalty_weight = 50000000000.0; // Hohes Gewicht!
        
        for (int v = 0; v < num_via; ++v) {
            // Wir bestrafen den Abstand zum Quadrat
            // obj = obj + weight * distance^2
            obj += penalty_weight * min_segment_dists_sq[v];
        }
    }

    opts["print_time"] = false;
    opts["record_time"] = false;
    opts["ipopt.print_level"] = 0;
    opts["ipopt.max_iter"] = maxIterations;
    opts["ipopt.tol"] = maxTol;
    opts["ipopt.linear_solver"] = "mumps";

    std::string filename = "../examples/solver_log_" + CasadiSystemName + ".txt";
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

void CasadiSystem::setupCasadiEXP()
{
    qDebug() << "setupCasadiEXP";
    using namespace casadi;
    MX all_x = MX::vertcat({});
    MX all_g = MX::vertcat({});
    MX all_p = MX::vertcat({});
    MX obj = 1;
    M = static_cast<int>(m_muscles.size());

    for (size_t m = 0; m < m_muscles.size(); ++m) {
        SSMuscle* mus = m_muscles[m];
        int num_inner = mus->MNodes.size() - 2;
        int num_wrap = mus->meshPtrs.size();

        // 1. Via-Points identifizieren
        std::vector<int> via_indices;
        for (int j = 0; j < num_wrap; ++j) {
            if (mus->meshPtrs[j]->bIsViaPoint) via_indices.push_back(j);
        }
        int num_via = via_indices.size();
        int num_eta_per_node = num_wrap - num_via;

        // --- VARIABLEN & PARAMETER ---
        MX x_mus = MX::sym("x_" + std::to_string(m), num_inner * 3 + num_inner * num_eta_per_node);
        all_x = MX::vertcat({all_x, x_mus});

        MX p_mus = MX::sym("p_" + std::to_string(m), 12 * num_wrap + 6);
        all_p = MX::vertcat({all_p, p_mus});

        MX P_orig = p_mus(Slice(12 * num_wrap, 12 * num_wrap + 3));
        MX P_ins  = p_mus(Slice(12 * num_wrap + 3, 12 * num_wrap + 6));
        int eta_offset = num_inner * 3;

        // Container für SoftMin: Sammle ALLE Abstände für jeden ViaPoint
        std::vector<std::vector<MX>> all_dists_per_vp(num_via);
        int current_via_idx = 0;

        // --- SCHLEIFE ÜBER KNOTEN ---
        for (int k = 0; k < num_inner; ++k) {
            MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
            MX g_prev = (k == 0) ? P_orig : x_mus(Slice((k - 1) * 3, k * 3));
            MX g_next = (k == num_inner - 1) ? P_ins : x_mus(Slice((k + 1) * 3, (k + 2) * 3));

            MX total_contact_force = 0;
            int current_eta_idx = 0;
            current_via_idx = 0;

            for (int j = 0; j < num_wrap; ++j) {
                MX q_j = p_mus(Slice(j * 12, (j + 1) * 12));

                if (mus->meshPtrs[j]->bIsViaPoint) {
                    // --- VIA POINT: Abstände sammeln ---
                    MX p_via = q_j(Slice(0, 3));
                    MX d2 = sumsqr(g_k - p_via); // Quadratischer Abstand
                    all_dists_per_vp[current_via_idx].push_back(d2);
                    current_via_idx++;
                } else {
                    // --- HINDERNISSE (Bleibt gleich) ---
                    MX eta_kj = x_mus(eta_offset + k * num_eta_per_node + current_eta_idx);
                    MX h = mus->meshPtrs[j]->constraintDistance(g_k, q_j);
                    MX grad_h_full = MX::gradient(h, x_mus);
                    MX grad_h = grad_h_full(Slice(k * 3, (k + 1) * 3));

                    total_contact_force += eta_kj * grad_h;
                    all_g = MX::vertcat({all_g, h * eta_kj, h, eta_kj}); 
                    current_eta_idx++;
                }
            }

            // Euler-Lagrange
            MX eq_el = (g_next - 2.0 * g_k + g_prev) + total_contact_force;
            all_g = MX::vertcat({all_g, eq_el});
            
            // Objective
            if (objType == 1) obj += norm_2(g_next - g_k);
            else if (objType == 2) { /*...*/ }
        }

        // --- GLOBAL VIA-POINT CONSTRAINTS (SOFTMIN) ---
        // Hier passiert die Magie: SoftMin Berechnung & Constraint
        
        double alpha = 50.0; // Härte des Minimums (Alpha=50 ist ein guter Startwert)

        for (int v = 0; v < num_via; ++v) {
            // 1. Alle Abstände als Vektor
            MX dists = MX::vertcat(all_dists_per_vp[v]);
            
            // 2. SoftMin (LogSumExp)
            // Dies ist der "glatte Ersatz" für min(dists)
            MX soft_min_d2 = -1.0/alpha * logsumexp(-alpha * dists);
            
            // 3. Constraint: SoftMin <= Toleranz^2
            double tol_sq = pow(mus->meshPtrs[via_indices[v]]->MViaPointTolerance, 2);
            
            // Wir fügen die Ungleichung hinzu: soft_min_d2 - tol_sq <= 0
            // In solveStep muss dafür ubg=0 gesetzt werden (lbg=-inf)
            all_g = MX::vertcat({all_g, soft_min_d2 - tol_sq}); 
        }
    }

    opts["ipopt.print_level"] = 0;
    opts["ipopt.max_iter"] = maxIterations;
    opts["ipopt.tol"] = maxTol;
    opts["ipopt.linear_solver"] = "mumps";

    std::string filename = "../examples/solver_log_" + CasadiSystemName + ".txt";
    opts["ipopt.output_file"] = filename; 
    opts["ipopt.file_print_level"] = 5;

    MXDict nlp = {{"x", all_x}, {"f", obj}, {"g", all_g}, {"p", all_p}};
    solver = nlpsol("f", "ipopt", nlp, opts);

    qDebug() << "CasadiSystem initialized with"
             << M << "muscles,"
             << "Total variables:" << all_x.size1()
             << "Total constraints:" << all_g.size1();
}

void CasadiSystem::solveStepPenalty(){
    using namespace casadi;

    std::vector<double> x0_all, p_all, lbg_all, ubg_all, lbx_all, ubx_all;

    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;
        int num_wrap = mus->meshPtrs.size();
        
        // Via-Points zählen
        int num_via = 0;
        for (auto* m : mus->meshPtrs) if (m->bIsViaPoint) num_via++;
        int num_eta_per_node = num_wrap - num_via;
        int num_etas_total = num_inner * num_eta_per_node;

        // --- INITIAL GUESS (x0) ---
        for (size_t k = 1; k < mus->MNodes.size() - 1; ++k) {
            MWMath::Point3D pWarm = mus->MNodes[k].predictNewGlobal();
            
            // OPTIONAL: Initial Guess Hilfe für Via-Points
            // Wenn der Punkt sehr nah an einem ViaPoint ist, setz ihn drauf.
            // Das hilft der Penalty-Methode enorm beim Start.
            /* for(auto* mesh : mus->meshPtrs) {
                if(mesh->bIsViaPoint) {
                    if(MWMath::distance(pWarm, mesh->PositionGlobal) < 0.2) { 
                         pWarm = mesh->PositionGlobal; 
                    }
                }
            }
            */
            
            x0_all.push_back(pWarm.x);
            x0_all.push_back(pWarm.y);
            x0_all.push_back(pWarm.z);
        }
        
        // Etas (nur für Hindernisse!)
        if (mus->lastEtas.size() == (size_t)num_etas_total && bUseWarmstartEtas) {
            x0_all.insert(x0_all.end(), mus->lastEtas.begin(), mus->lastEtas.end());
        } else {
            for (int e = 0; e < num_etas_total; ++e) x0_all.push_back(0.001);
        }

        // --- PARAMETER (p) ---
        for (auto* mesh : mus->meshPtrs) {
            // Position
            p_all.push_back(mesh->PositionGlobal.x);
            p_all.push_back(mesh->PositionGlobal.y);
            p_all.push_back(mesh->PositionGlobal.z);
            // Rotation
            for (int col = 0; col < 3; ++col) {
                for (int row = 0; row < 3; ++row) {
                    p_all.push_back(mesh->OrientationGlobal.m[row][col]);
                }
            }
        }
        p_all.push_back(mus->OriginPointGlobal.x); 
        p_all.push_back(mus->OriginPointGlobal.y); 
        p_all.push_back(mus->OriginPointGlobal.z);
        p_all.push_back(mus->InsertionPointGlobal.x); 
        p_all.push_back(mus->InsertionPointGlobal.y); 
        p_all.push_back(mus->InsertionPointGlobal.z);

        // --- CONSTRAINTS SCHRANKEN (lbg/ubg) ---
        for (int i = 0; i < num_inner; ++i) {
            for (int j = 0; j < num_wrap; ++j) {
                // WICHTIG: Bounds nur für Hindernisse hinzufügen, NICHT für ViaPoints
                if (!mus->meshPtrs[j]->bIsViaPoint) {
                    lbg_all.push_back(0.0); ubg_all.push_back(0.0);  // Komplementarität
                    lbg_all.push_back(0.0); ubg_all.push_back(1e10); // Abstand >= 0
                    lbg_all.push_back(0.0); ubg_all.push_back(1e10); // Kraft >= 0
                }
            }
            // EL = 0 (x,y,z)
            for (int d = 0; d < 3; ++d) { 
                lbg_all.push_back(-ELTolerance); ubg_all.push_back(ELTolerance); 
            }
        }

        // -----------------------------------------------------------------
        // WICHTIG: DIESEN LOOP LÖSCHEN !!!
        // (Er war verantwortlich für den Crash, weil die Constraint nicht mehr existiert)
        /*
        for (int v = 0; v < num_via; ++v) {
            lbg_all.push_back(-inf); 
            ubg_all.push_back(0.0);
        }
        */
        // -----------------------------------------------------------------

        // --- VARIABLEN SCHRANKEN (lbx/ubx) ---
        for (int i = 0; i < num_inner * 3; ++i) {
            lbx_all.push_back(-inf); ubx_all.push_back(inf);
        }
        for (int i = 0; i < num_etas_total; ++i) {
            lbx_all.push_back(0.0); ubx_all.push_back(inf);
        }
    }

    // --- SOLVER AUFRUF ---
    DMDict arg = {{"x0", x0_all}, {"p", p_all}, {"lbg", lbg_all}, {"ubg", ubg_all}, {"lbx", lbx_all}, {"ubx", ubx_all}};
    DMDict res = solver(arg);

    // --- ERGEBNISSE EXTRAHIEREN ---
    std::vector<double> res_x = std::vector<double>(res.at("x"));
    auto solverInfo = solver.stats();
    std::string status = solverInfo.at("return_status");

    if (status != "Solve_Succeeded") {
        qDebug() << "    Warnung: Solver Status:" << QString::fromStdString(status);
    }
    qDebug() << "    Solver converged in" << int(solverInfo.at("iter_count")) << "iterations";

    // (Rest der Funktion für Datenspeicherung bleibt identisch...)
    int current_x_offset = 0;
    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;
        int num_wrap = mus->meshPtrs.size();
        int num_via = 0;
        for (auto* m : mus->meshPtrs) if (m->bIsViaPoint) num_via++;
        int num_etas_mus = num_inner * (num_wrap - num_via);

        for (int i = 0; i < num_inner; ++i) {
            int mIdx = i + 1;
            MWMath::Point3D p;
            p.x = res_x[current_x_offset + i * 3 + 0];
            p.y = res_x[current_x_offset + i * 3 + 1];
            p.z = res_x[current_x_offset + i * 3 + 2];
            mus->MusclePointsGlobal[mIdx] = p;
            mus->MNodes[mIdx].PositionGlobal = p;
        }
        int eta_start_idx = current_x_offset + (num_inner * 3);
        if(num_etas_mus > 0) {
             mus->lastEtas.assign(res_x.begin() + eta_start_idx, res_x.begin() + eta_start_idx + num_etas_mus);
        }
        current_x_offset += (num_inner * 3 + num_etas_mus);
    }
}

void CasadiSystem::solveStepEXP()
{
    qDebug() << "solveStepEXP";
    using namespace casadi;

    std::vector<double> x0_all, p_all, lbg_all, ubg_all, lbx_all, ubx_all;

    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;
        int num_wrap = mus->meshPtrs.size();
        
        // Zählen der Variablen-Größen
        int num_via = 0;
        for (auto* m : mus->meshPtrs) if (m->bIsViaPoint) num_via++;
        int num_eta_per_node = num_wrap - num_via;
        int num_etas_total = num_inner * num_eta_per_node;

        // -------------------------------------------------------------------
        // 1. INITIAL GUESS (x0) MIT "SNAPPING"
        // -------------------------------------------------------------------
        // Wir berechnen den Vektor vom Start zum Ziel des Muskels.
        // Dann projizieren wir jeden ViaPoint auf diese Linie, um zu sehen,
        // bei welchem "Prozentsatz" (t) des Muskels er liegen sollte.
        // Den entsprechenden Knoten setzen wir hart auf die ViaPoint-Position.
        

        for (size_t k = 1; k < mus->MNodes.size() - 1; ++k) {
            // A) Standard: Vorhersage aus letztem Zeitschritt
            MWMath::Point3D pWarm = mus->MNodes[k].predictNewGlobal();

            x0_all.push_back(pWarm.x);
            x0_all.push_back(pWarm.y);
            x0_all.push_back(pWarm.z);
        }
        
        // --- ETAS (Kontaktkräfte) ---
        // Achtung: Nur für Hindernisse, nicht für ViaPoints!
        if (mus->lastEtas.size() == (size_t)num_etas_total && bUseWarmstartEtas) {
            x0_all.insert(x0_all.end(), mus->lastEtas.begin(), mus->lastEtas.end());
        } else {
            for (int e = 0; e < num_etas_total; ++e) x0_all.push_back(0.001);
        }

        // -------------------------------------------------------------------
        // 2. PARAMETER (p)
        // -------------------------------------------------------------------
        for (auto* mesh : mus->meshPtrs) {
            p_all.push_back(mesh->PositionGlobal.x);
            p_all.push_back(mesh->PositionGlobal.y);
            p_all.push_back(mesh->PositionGlobal.z);
            for (int col = 0; col < 3; ++col) {
                for (int row = 0; row < 3; ++row) {
                    p_all.push_back(mesh->OrientationGlobal.m[row][col]);
                }
            }
        }
        p_all.push_back(mus->OriginPointGlobal.x); 
        p_all.push_back(mus->OriginPointGlobal.y); 
        p_all.push_back(mus->OriginPointGlobal.z);
        p_all.push_back(mus->InsertionPointGlobal.x); 
        p_all.push_back(mus->InsertionPointGlobal.y); 
        p_all.push_back(mus->InsertionPointGlobal.z);

        // -------------------------------------------------------------------
        // 3. CONSTRAINTS SCHRANKEN (lbg/ubg)
        // -------------------------------------------------------------------
        // Hier muss die Struktur EXAKT zu setupCasadiSoftMin passen!
        
        for (int i = 0; i < num_inner; ++i) {
            for (int j = 0; j < num_wrap; ++j) {
                // WICHTIG: Bounds nur für echte Hindernisse hinzufügen!
                // In setupCasadiSoftMin wird der ViaPoint-Fall im 'if' behandelt (kein 'all_g' update),
                // und der Obstacle-Fall im 'else' (3x 'all_g' update).
                
                if (!mus->meshPtrs[j]->bIsViaPoint) {
                    lbg_all.push_back(0.0); ubg_all.push_back(0.0);  // 1. Komplementarität (h*eta = 0)
                    lbg_all.push_back(0.0); ubg_all.push_back(1e10); // 2. Nicht-Eindringen (h >= 0)
                    lbg_all.push_back(0.0); ubg_all.push_back(1e10); // 3. Nur Druckkraft (eta >= 0)
                }
            }
            
            // Euler-Lagrange Gleichgewicht (Fx, Fy, Fz = 0)
            for (int d = 0; d < 3; ++d) { 
                lbg_all.push_back(-ELTolerance); 
                ubg_all.push_back(ELTolerance); 
            }
        }

        // WICHTIG: KEIN Loop über Via-Points am Ende!
        // Da wir Penalty/SoftMin nutzen, gibt es keine globalen Constraints mehr in 'all_g'.

        // -------------------------------------------------------------------
        // 4. VARIABLEN SCHRANKEN (lbx/ubx)
        // -------------------------------------------------------------------
        // Knoten sind frei (-inf bis inf)
        for (int i = 0; i < num_inner * 3; ++i) {
            lbx_all.push_back(-inf); 
            ubx_all.push_back(inf);
        }
        // Kontaktkräfte sind positiv (0 bis inf)
        for (int i = 0; i < num_etas_total; ++i) {
            lbx_all.push_back(0.0); 
            ubx_all.push_back(inf);
        }

        // Global Via-Point Constraints
        for (int v = 0; v < num_via; ++v) {
            lbg_all.push_back(-inf); 
            ubg_all.push_back(0.0); // Constraint <= 0
        }
    }

    // -------------------------------------------------------------------
    // 5. SOLVER AUFRUF
    // -------------------------------------------------------------------
    DMDict arg = {{"x0", x0_all}, {"p", p_all}, {"lbg", lbg_all}, {"ubg", ubg_all}, {"lbx", lbx_all}, {"ubx", ubx_all}};
    
    // Debugging (falls es immer noch crasht, Dimensionen prüfen)
    /*
    qDebug() << "Solver Inputs:";
    qDebug() << "  x0 size:" << x0_all.size();
    qDebug() << "  lbg size:" << lbg_all.size();
    qDebug() << "  lbx size:" << lbx_all.size();
    */

    DMDict res = solver(arg);

    // -------------------------------------------------------------------
    // 6. ERGEBNISSE SPEICHERN
    // -------------------------------------------------------------------
    std::vector<double> res_x = std::vector<double>(res.at("x"));
    
    auto solverInfo = solver.stats();
    if (solverInfo.count("return_status")) {
        std::string status = solverInfo.at("return_status");
        if (status != "Solve_Succeeded") {
             qDebug() << "Warning: Solver Status:" << QString::fromStdString(status);
        }
    }

    int current_x_offset = 0;
    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;
        int num_wrap = mus->meshPtrs.size();
        int num_via = 0;
        for (auto* m : mus->meshPtrs) if (m->bIsViaPoint) num_via++;
        int num_etas_mus = num_inner * (num_wrap - num_via);

        // Punkte zurückschreiben
        for (int i = 0; i < num_inner; ++i) {
            int mIdx = i + 1;
            MWMath::Point3D p;
            p.x = res_x[current_x_offset + i * 3 + 0];
            p.y = res_x[current_x_offset + i * 3 + 1];
            p.z = res_x[current_x_offset + i * 3 + 2];

            mus->MusclePointsGlobal[mIdx] = p;
            mus->MNodes[mIdx].PositionGlobal = p;
        }

        // Etas speichern (für nächsten Warmstart)
        int eta_start_idx = current_x_offset + (num_inner * 3);
        if (num_etas_mus > 0) {
            mus->lastEtas.assign(
                res_x.begin() + eta_start_idx, 
                res_x.begin() + eta_start_idx + num_etas_mus
            );
            for (auto& node : mus->MNodes) node.lastEtas = mus->lastEtas;
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

        // 1. Via-Points identifizieren
        std::vector<int> via_indices;
        for (int j = 0; j < num_wrap; ++j) {
            if (mus->meshPtrs[j]->bIsViaPoint) via_indices.push_back(j);
        }
        int num_via = via_indices.size();
        int num_eta_per_node = num_wrap - num_via;

        // --- VARIABLEN: [Punkte (3*n), Etas (n * (Meshes-Vias))] ---
        MX x_mus = MX::sym("x_" + std::to_string(m), num_inner * 3 + num_inner * num_eta_per_node);
        all_x = MX::vertcat({all_x, x_mus});

        // --- PARAMETER: [Meshes (12*m), Origin (3), Insertion (3)] ---
        MX p_mus = MX::sym("p_" + std::to_string(m), 12 * num_wrap + 6);
        all_p = MX::vertcat({all_p, p_mus});

        MX P_orig = p_mus(Slice(12 * num_wrap, 12 * num_wrap + 3));
        MX P_ins  = p_mus(Slice(12 * num_wrap + 3, 12 * num_wrap + 6));
        int eta_offset = num_inner * 3;

        // Container für Minimum-Abstände (einer pro Via-Point über alle Knoten k)
        std::vector<MX> min_dists_sq(num_via);

        // --- GLEICHUNGEN (g) ---
        for (int k = 0; k < num_inner; ++k) {
            MX g_k = x_mus(Slice(k * 3, (k + 1) * 3));
            MX g_prev = (k == 0) ? P_orig : x_mus(Slice((k - 1) * 3, k * 3));
            MX g_next = (k == num_inner - 1) ? P_ins : x_mus(Slice((k + 1) * 3, (k + 2) * 3));

            MX total_contact_force = 0;
            int current_eta_idx = 0;
            int current_via_idx = 0;

            for (int j = 0; j < num_wrap; ++j) {
                MX q_j = p_mus(Slice(j * 12, (j + 1) * 12));

                if (mus->meshPtrs[j]->bIsViaPoint) {
                    // Via-Point Logik: Minimum-Suche über alle Knoten k
                    MX p_via = q_j(Slice(0, 3));
                    MX d2 = sumsqr(g_k - p_via); // Quadratischer Abstand Punkt k zu Via-Point j
                    
                    if (k == 0) {
                        // Beim ersten inneren Knoten: Initialisiere den Container
                        min_dists_sq[current_via_idx] = d2;
                    } else {
                        // Nutze casadi::fmin um das Minimum zu finden
                        min_dists_sq[current_via_idx] = fmin(min_dists_sq[current_via_idx], d2);
                    }
                    current_via_idx++;
                } else {
                    // Normale Signorini-Logik (Hindernisse)
                    MX eta_kj = x_mus(eta_offset + k * num_eta_per_node + current_eta_idx);
                    MX h = mus->meshPtrs[j]->constraintDistance(g_k, q_j);

                    MX grad_h;
                    if (bUseOwnGradient){
                        grad_h = mus->meshPtrs[j]->constraintJacobian(g_k, q_j);
                        qDebug() << "using own gradient!";
                    }
                    else {
                        MX grad_h_full = MX::gradient(h, x_mus);
                        grad_h = grad_h_full(Slice(k * 3, (k + 1) * 3));
                        qDebug() << "using auto-casadi gradient!";
                    }

                    total_contact_force -= eta_kj * grad_h;
                    
                    // Signorini/Complementarity
                    all_g = MX::vertcat({all_g, h * eta_kj, h, eta_kj}); 
                    current_eta_idx++;
                }
            }

            // Euler-Lagrange Kraftbilanz
            MX eq_el = (2.0 * g_k  - g_prev - g_next) - total_contact_force; // (g_next - 2.0 * g_k + g_prev) + total_contact_force;
            all_g = MX::vertcat({all_g, eq_el});
            
            // Objective
            if (objType == 1) {
                obj = obj + norm_2(g_next - g_k);
            }
            else if (objType == 2) {
                MX prev_seg_sq = sumsqr(g_k - g_prev);
                MX current_seg_sq = sumsqr(g_next - g_k);
                double weight_equi = 10.0; 
                obj += weight_equi * sumsqr(prev_seg_sq - current_seg_sq);
            }
        }

        // Global Via-Point Constraints (one Constraint per Via-Point for whole muscle)
        for (int v = 0; v < num_via; ++v) {
            double tol_sq = pow(mus->meshPtrs[via_indices[v]]->MViaPointTolerance, 2);
            all_g = MX::vertcat({all_g, min_dists_sq[v] - tol_sq}); // d_min^2 - tol^2 <= 0
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

    qDebug() << "CasadiSystem initialized with"
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
        
        // Via-Points zählen
        int num_via = 0;
        for (auto* m : mus->meshPtrs) if (m->bIsViaPoint) num_via++;
        int num_eta_per_node = num_wrap - num_via;
        int num_etas_total = num_inner * num_eta_per_node;

        // --- INITIAL GUESS (x0) ---
        for (size_t k = 1; k < mus->MNodes.size() - 1; ++k) {
            MWMath::Point3D pWarm = mus->MNodes[k].predictNewGlobal();
            x0_all.push_back(pWarm.x);
            x0_all.push_back(pWarm.y);
            x0_all.push_back(pWarm.z);
        }
        
        if (mus->lastEtas.size() == (size_t)num_etas_total && bUseWarmstartEtas) {
            qDebug() << "    Warmstart etas:" << QString::fromStdString(mus->Name);
            x0_all.insert(x0_all.end(), mus->lastEtas.begin(), mus->lastEtas.end());
        } else {
            qDebug() << "    No warmstart etas:" << QString::fromStdString(mus->Name);
            for (int e = 0; e < num_etas_total; ++e) x0_all.push_back(0.0);
        }

        // --- PARAMETER (p) ---
        for (auto* mesh : mus->meshPtrs) {
            p_all.push_back(mesh->PositionGlobal.x);    inputParamsForDebug.push_back(mesh->PositionGlobal.x);  inputParamDescriptionsForDebug.push_back(mesh->Name + "_PosX");
            p_all.push_back(mesh->PositionGlobal.y);    inputParamsForDebug.push_back(mesh->PositionGlobal.y);  inputParamDescriptionsForDebug.push_back(mesh->Name + "_PosY");
            p_all.push_back(mesh->PositionGlobal.z);    inputParamsForDebug.push_back(mesh->PositionGlobal.z);  inputParamDescriptionsForDebug.push_back(mesh->Name + "_PosZ");
            /* for (int col = 0; col < 3; ++col) {
                for (int row = 0; row < 3; ++row) {
                    p_all.push_back(mesh->OrientationGlobal.transposed().m[row][col]);
                }
            } */
           MWMath::RotMatrix3x3 R = mesh->OrientationGlobal;
           /* for (int row = 0; row < 3; ++row) {       // <--- Erst Zeile
                for (int col = 0; col < 3; ++col) {   // <--- Dann Spalte
                    p_all.push_back(R.m[row][col]);
                    inputParamsForDebug.push_back(R.m[row][col]);
                    inputParamDescriptionsForDebug.push_back(mesh->Name + "_Ori_" + std::to_string(row) + std::to_string(col));
                }
            } */
           for (int col = 0; col < 3; ++col) {       // <--- Spalte ist außen!
                for (int row = 0; row < 3; ++row) {   // <--- Zeile ist innen!
                    // Wir greifen aber ganz normal auf m[row][col] zu
                    p_all.push_back(R.m[row][col]);
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
                if (!mus->meshPtrs[j]->bIsViaPoint) {
                    // Komplementarität: h * eta = 0
                    lbg_all.push_back(0.0); ubg_all.push_back(0.0); 
                    // Nicht-Eindringen: h >= 0
                    lbg_all.push_back(0.0); ubg_all.push_back(inf);
                    // Nur Druckkraft: eta >= 0
                    lbg_all.push_back(0.0); ubg_all.push_back(inf);
                }
            }
            // EL = 0 (für x, y, z)
            for (int d = 0; d < 3; ++d) { 
                lbg_all.push_back(-ELTolerance); 
                ubg_all.push_back(ELTolerance); 
            }
        }

        /* // Via-Point Constraints: d_min^2 - tol^2 <= 0
        for (int v = 0; v < num_via; ++v) {
            lbg_all.push_back(-inf); 
            ubg_all.push_back(0.0);
        } */

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
        qDebug() << "    Warnung: Multi-Muscle Solver Status:" << QString::fromStdString(status);
    }
    int convSteps = int(solverInfo.at("iter_count"));
    qDebug() << "    Solver converged in" << convSteps << "iterations";
    SolverConvergenceMessages.push_back(QString::fromStdString(status).toStdString());
    SolverConvergenceSteps.push_back(convSteps);
    
    int current_x_offset = 0;
    for (auto* mus : m_muscles) {
        int num_inner = mus->MNodes.size() - 2;
        int num_wrap = mus->meshPtrs.size();
        int num_via = 0;
        for (auto* m : mus->meshPtrs) if (m->bIsViaPoint) num_via++;
        int num_etas_mus = num_inner * (num_wrap - num_via);

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


