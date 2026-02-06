/* 
#include "rawSimulation.h"
#include "model/MWTissueObject.h"
#include "model/MWTissueObjects/MWBody.h"
#include <vtkMatrix4x4.h>
#include <cmath>


void RawSimulation::run() {
    allMusclesResults.clear();
    allMeshesRMatrixResults.clear();

    auto cyl = std::make_shared<MWEllipsoidMesh>(0.8, 0.8, 1.6);
    cyl->PositionGlobalTest = MWMath::Point3D(0,0,-1.0);
    cyl->AttachedBody = new MWBody("pseudoBody", MWMath::Point3D(0,0,0), MWMath::Point3D(0,0,1), 0, MWMath::Point3D(0.6,0.1,0.1), nullptr, nullptr);
    
    // Initialisiere RMatrix Slot für das erste Mesh
    allMeshesRMatrixResults.push_back(std::vector<MWMath::RotMatrix3x3>());

    std::vector<double> angles;
    for(int i=0; i<=90; i+=10) angles.push_back((double)i);

    auto muscle1 = std::make_shared<MWMuscle>("M1", nullptr, MWMath::Point3D(1,0,0), nullptr, MWMath::Point3D(-1,0.5,0), 10);
    muscle1->updateIndiviual();
    muscle1->createMusclePoints();

    runMovementSequence(muscle1, {cyl}, angles);

    VTKSimViewer viewer(allMusclesResults, allMeshesRMatrixResults, {cyl}, angles);
    viewer.exec();

    delete cyl->AttachedBody;
}

void RawSimulation::runMovementSequence(std::shared_ptr<MWMuscle> muscle, std::vector<std::shared_ptr<MWMesh>> meshes, const std::vector<double> &zRotationSteps) {
    allMusclesResults.push_back({});
    size_t mIdx = allMusclesResults.size() - 1;

    // Wir speichern das letzte Ergebnis für den Warmstart
    std::vector<double> last_result;

    for (size_t i = 0; i < zRotationSteps.size(); ++i) {
        // 1. Rotation setzen
        MWMath::RotMatrix3x3 newRot = MWMath::axisAngle(MWMath::Point3D(0,1,0), zRotationSteps[i]); 
        static_cast<MWBody*>(meshes[0]->AttachedBody)->setRMatrixGlobal(newRot);
        allMeshesRMatrixResults[0].push_back(newRot);

        // 2. Solver aufrufen (mit Warmstart-Daten)
        std::vector<double> result = solveSingleMuscleGlobal(muscle, meshes, last_result);

        if (!result.empty()) {
            last_result = result; // Für den nächsten Schritt merken!
            
            int numInner = muscle->MusclePointsCount - 2;
            std::vector<MWMath::Point3D> inner = MWMath::MWDoubleVec2Point3DVec(
                std::vector<double>(result.begin(), result.begin() + numInner * 3)
            );

            std::vector<MWMath::Point3D> fullPath;
            fullPath.push_back(muscle->OriginPointGlobal);
            for(auto& p : inner) fullPath.push_back(p);
            fullPath.push_back(muscle->InsertionPointGlobal);

            allMusclesResults[mIdx].push_back(fullPath);
        }
    }
}

std::vector<double> RawSimulation::solveSingleMuscleGlobal(std::shared_ptr<MWMuscle> muscle, const std::vector<std::shared_ptr<MWMesh>>& activeMeshes,const std::vector<double>& warmstart) {
    using namespace casadi;
    int N = muscle->MusclePointsCount;
    int numInner = N - 2;
    double ds = 1.0 / (N - 1.0);

    std::vector<MX> y_syms;
    for (int k = 0; k < numInner; ++k) y_syms.push_back(MX::sym("y_" + std::to_string(k), 3));
    MX eta_vars = MX::sym("eta", numInner);
    
    // 1. KORREKTES MX::vertcat
    std::vector<MX> all_vars = y_syms;
    all_vars.push_back(eta_vars);
    MX x_nlp = MX::vertcat(all_vars); 

    MX y_start = MWMath::MWPoint2MXConst(muscle->OriginPointGlobal);
    MX y_end = MWMath::MWPoint2MXConst(muscle->InsertionPointGlobal);
    
    MX g = MX(); 
    std::vector<double> lbg, ubg;
    MX objective = 0; // Wir minimieren die Energie/Länge!

    for (int k = 0; k < numInner; ++k) {
        MX y_k = y_syms[k]; 
        MX y_prev = (k == 0) ? y_start : y_syms[k-1];
        MX y_next = (k == numInner - 1) ? y_end : y_syms[k+1];
        
        // Energie für das Objective (Längenminimierung)
        objective += sumsqr(y_k - y_prev) + sumsqr(y_next - y_k);

        MX F_eq = (1.0 / ds) * (2.0 * y_k - y_prev - y_next);
        
        MX phi_k = 0; 
        MX grad_phi_k = MX::zeros(3, 1);
        for (auto& mesh : activeMeshes) {
            phi_k += mesh->ConstraintMeshDistance(y_k); 
            grad_phi_k += mesh->ConstraintMeshJacobian(y_k);
        }

        MX Fk = F_eq - (ds / 2.0) * eta_vars(k) * grad_phi_k;

        // 2. KORREKTES MX::vertcat mit Initialisierungsliste
        g = MX::vertcat({g, Fk, phi_k, eta_vars(k), phi_k * eta_vars(k)});

        for(int i=0; i<3; ++i) { lbg.push_back(0); ubg.push_back(0); }
        lbg.push_back(0); ubg.push_back(inf);
        lbg.push_back(0); ubg.push_back(inf);
        lbg.push_back(0); ubg.push_back(0);
    }

    // Solver Setup (Objective f hinzufügen!)
    Dict opts;
    opts["ipopt.print_level"] = 5;
    opts["ipopt.tol"] = 1e-6;
    opts["ipopt.max_iter"] = 500;
    //opts["ipopt.mu_strategy"] = "adaptive"; // Hilft bei Infeasibility
    //opts["ipopt.acceptable_tol"] = 1e-2;
    Function solver = nlpsol("S", "ipopt", {{"x", x_nlp}, {"f", 0}, {"g", g}}, opts);
    
    // 3. WARMSTART LOGIK
    std::vector<double> x0;
    if (warmstart.empty()) {
        // Initialer Guess: Linie
        for (int k = 1; k <= numInner; ++k) {
            double a = (double)k / (N - 1);
            x0.push_back(muscle->OriginPointGlobal.x + a*(muscle->InsertionPointGlobal.x - muscle->OriginPointGlobal.x));
            x0.push_back(muscle->OriginPointGlobal.y + a*(muscle->InsertionPointGlobal.y - muscle->OriginPointGlobal.y));
            x0.push_back(muscle->OriginPointGlobal.z + a*(muscle->InsertionPointGlobal.z - muscle->OriginPointGlobal.z));
        }
        for (int k = 0; k < numInner; ++k) x0.push_back(0.01);
    } else {
        x0 = warmstart;
    }

    DMDict sol = solver(DMDict{{"x0", x0}, {"lbg", lbg}, {"ubg", ubg}});
    return std::vector<double>(sol.at("x"));
}

VTKSimViewer::VTKSimViewer(const std::vector<std::vector<std::vector<MWMath::Point3D>>>& muscleResults, 
                           const std::vector<std::vector<MWMath::RotMatrix3x3>>& meshResults,
                           std::vector<std::shared_ptr<MWMesh>> meshes, 
                           const std::vector<double>& angles)
    : m_muscleResults(muscleResults), m_meshResults(meshResults), m_meshes(meshes), m_angles(angles) 
{
    setWindowTitle("MW VTK 3D Analysis");
    resize(1200, 800);
    auto* layout = new QVBoxLayout(this);
    m_vtkWidget = new QVTKOpenGLNativeWidget(this);
    layout->addWidget(m_vtkWidget);
    m_slider = new QSlider(Qt::Horizontal);
    m_slider->setRange(0, (int)angles.size() - 1);
    layout->addWidget(m_slider);
    setupVtkPipeline();
    connect(m_slider, &QSlider::valueChanged, this, &VTKSimViewer::updateStep);
    updateStep(0);
}

void VTKSimViewer::updateMeshTransform(vtkActor* actor, const MWMath::RotMatrix3x3& R) {
    vtkSmartPointer<vtkMatrix4x4> vtkMat = vtkSmartPointer<vtkMatrix4x4>::New();
    vtkMat->Identity();
    // pos
    vtkMat->SetElement(0, 3, m_meshes[0]->PositionGlobalTest.x);
    vtkMat->SetElement(1, 3, m_meshes[0]->PositionGlobalTest.y);
    vtkMat->SetElement(2, 3, m_meshes[0]->PositionGlobalTest.z);
    for (int i = 0; i < 3; ++i) 
        for (int j = 0; j < 3; ++j) 
            vtkMat->SetElement(i, j, R.m[i][j]);
    actor->SetUserMatrix(vtkMat);
    qDebug() << "mesh... vtkMat:";
    vtkMat->Print(std::cout);
}

void VTKSimViewer::setupVtkPipeline() {
    m_renderer = vtkSmartPointer<vtkRenderer>::New();
    m_vtkWidget->renderWindow()->AddRenderer(m_renderer);
    m_renderer->SetBackground(0.1, 0.1, 0.1);

    for (auto& mesh : m_meshes) {
        auto sphere = vtkSmartPointer<vtkSphereSource>::New();
        sphere->SetThetaResolution(40); sphere->SetPhiResolution(40);
        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(sphere->GetOutputPort());
        auto actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        if (auto ell = std::dynamic_pointer_cast<MWEllipsoidMesh>(mesh)) actor->SetScale(ell->A, ell->B, ell->C);
        actor->GetProperty()->SetOpacity(0.5); actor->GetProperty()->SetColor(1, 0, 0);
        m_renderer->AddActor(actor);
        m_meshActors.push_back(actor);
    }

    for (size_t m = 0; m < m_muscleResults.size(); ++m) {
        auto points = vtkSmartPointer<vtkPoints>::New();
        auto line = vtkSmartPointer<vtkPolyData>::New();
        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(line);
        auto actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetLineWidth(4); actor->GetProperty()->SetColor(0, 0.7, 1.0);
        m_renderer->AddActor(actor);
        m_muscleActors.push_back(actor);
        m_muscleLines.push_back(line);
        m_musclePointsData.push_back(points);
    }
}

void VTKSimViewer::updateStep(int step) {
    for (size_t i = 0; i < m_meshActors.size(); ++i) {
        if (i < m_meshResults.size() && step < (int)m_meshResults[i].size()) 
            updateMeshTransform(m_meshActors[i], m_meshResults[i][step]);
    }
    for (size_t m = 0; m < m_muscleResults.size(); ++m) {
        if (step >= (int)m_muscleResults[m].size()) continue;
        const std::vector<MWMath::Point3D>& path = m_muscleResults[m][step];
        m_musclePointsData[m]->Reset();
        auto polyLine = vtkSmartPointer<vtkPolyLine>::New();
        polyLine->GetPointIds()->SetNumberOfIds(path.size());
        for (size_t i = 0; i < path.size(); ++i) {
            m_musclePointsData[m]->InsertNextPoint(path[i].x, path[i].y, path[i].z);
            polyLine->GetPointIds()->SetId(i, i);
            qDebug() << "muscle point <<" << i << ":" << path[i].x << path[i].y << path[i].z;
        }
        auto cells = vtkSmartPointer<vtkCellArray>::New();
        cells->InsertNextCell(polyLine);
        m_muscleLines[m]->SetPoints(m_musclePointsData[m]);
        m_muscleLines[m]->SetLines(cells);
    }
    m_vtkWidget->renderWindow()->Render();
}




 */