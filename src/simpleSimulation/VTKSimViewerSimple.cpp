#include "VTKSimViewerSimple.h"
#include <QKeyEvent>


VTKSimViewerSimple::VTKSimViewerSimple(const std::vector<std::vector<std::vector<MWMath::Point3D>>>& muscleResults, 
                           const std::vector<std::vector<MWMath::RotMatrix3x3>>& meshResults,
                           const std::vector<std::vector<std::vector<MWMath::Point3D>>>& initialGuesses,
                           const std::vector<std::vector<std::vector<MWMath::Point3D>>>& initialGuessColors,
                           std::vector<std::shared_ptr<SSMesh>> meshes, 
                            const std::vector<SSMuscle*>& muscles,
                           const std::vector<double>& angles,
                            const std::vector<std::vector<MWMath::Point3D>>& otherPointsWithColors,
                           const std::vector<std::vector<std::string>>& solverInfoText,
                           double scalerCM,
                            QWidget* parent)
    : QDialog(parent), m_muscleResults(muscleResults), m_meshResults(meshResults), 
      m_initialGuesses(initialGuesses), m_initialGuessColors(initialGuessColors), 
      m_meshes(meshes), m_muscles(muscles), m_angles(angles), m_otherPointsWithColors(otherPointsWithColors),
      m_solverInfoText(solverInfoText), m_scalerCM(scalerCM)
{
    setWindowTitle("MW CasADi 3D Analysis - Multi-Muscle System");
    resize(1200, 800);

    auto* layout = new QVBoxLayout(this);
    m_vtkWidget = new QVTKOpenGLNativeWidget(this);
    layout->addWidget(m_vtkWidget);

    m_slider = new QSlider(Qt::Horizontal);
    m_slider->setRange(0, angles.empty() ? 0 : (int)angles.size() - 1);
    auto* bottomLayout = new QHBoxLayout();
    bottomLayout->addWidget(m_slider);

    m_plotButton = new QPushButton("Generate Plots", this);
    bottomLayout->addWidget(m_plotButton);

    // Füge den kombinierten unteren Bereich zum Hauptlayout hinzu
    layout->addLayout(bottomLayout);

    // Verbinde den Klick-Event des Buttons mit unserer neuen Funktion
    connect(m_plotButton, &QPushButton::clicked, this, &VTKSimViewerSimple::generatePlots);

    setupVtkPipeline();

    connect(m_slider, &QSlider::valueChanged, this, &VTKSimViewerSimple::updateStep);
    
    // Initialer Render
    if (!angles.empty()) updateStep(0);
}

void VTKSimViewerSimple::setupVtkPipeline() {
    qDebug() << "Setting up VTK pipeline for" << m_muscleResults.size() << "muscles...";
    m_renderer = vtkSmartPointer<vtkRenderer>::New();
    m_vtkWidget->renderWindow()->AddRenderer(m_renderer);
    m_renderer->SetBackground(0.15, 0.15, 0.15);

    // Step-Text Actor (HUD)
    m_stepTextActor = vtkSmartPointer<vtkTextActor>::New();
    m_stepTextActor->SetInput(("Step: 0 / " + std::to_string(m_angles.size() - 1)).c_str()); // Initial
    m_stepTextActor->GetTextProperty()->SetFontSize(24);
    m_stepTextActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0); // Weiß
    m_stepTextActor->SetDisplayPosition(20, 30); // oben links (Pixel)
    m_renderer->AddActor2D(m_stepTextActor);

    // 1. Mesh Pipeline (Ellipsoide/Zylinder/Torus)
    for (auto& mesh : m_meshes) {
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();

        if (auto ell = std::dynamic_pointer_cast<SSEllipsoidMesh>(mesh)) {
            auto source = vtkSmartPointer<vtkSphereSource>::New();
            source->SetRadius(1.0);
            source->SetThetaResolution(40);
            source->SetPhiResolution(40);
            mapper->SetInputConnection(source->GetOutputPort());
            actor->SetScale(ell->A , ell->B, ell->C);
        } 
        else if (auto cyl = std::dynamic_pointer_cast<SSCylinderMesh>(mesh)) {
            auto source = vtkSmartPointer<vtkCylinderSource>::New();
            source->SetRadius(cyl->Radius);
            source->SetHeight(cyl->Height);
            source->SetResolution(40);

            auto transform = vtkSmartPointer<vtkTransform>::New();
            transform->RotateX(90); // VTK Y-Achse zu Z-Achse
            auto transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
            transformFilter->SetInputConnection(source->GetOutputPort());
            transformFilter->SetTransform(transform);
            transformFilter->Update();
            mapper->SetInputConnection(transformFilter->GetOutputPort());
        } 
        else if (auto tor = std::dynamic_pointer_cast<SSTorusMesh>(mesh)) {
            auto source = vtkSmartPointer<vtkParametricTorus>::New();
            source->SetRingRadius(tor->R);
            source->SetCrossSectionRadius(tor->r);
            
            auto pSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
            pSource->SetParametricFunction(source);
            pSource->SetScalarModeToNone();
            pSource->SetUResolution(60);
            pSource->SetVResolution(60);
            mapper->SetInputConnection(pSource->GetOutputPort());
            
            actor->SetScale(tor->A, tor->B , tor->C);
        }
        
        actor->SetMapper(mapper);
        actor->GetProperty()->SetOpacity(0.5);
        actor->GetProperty()->SetColor(mesh->MeshColor.x, mesh->MeshColor.y, mesh->MeshColor.z);
        m_renderer->AddActor(actor);
        m_meshActors.push_back(actor);
    }

    // 2. Muskel Pipeline (Polylines) - FÜR MEHRERE MUSKELN
    for (size_t m = 0; m < m_muscleResults.size(); ++m) {
        auto mPoints = vtkSmartPointer<vtkPoints>::New();
        auto mLineData = vtkSmartPointer<vtkPolyData>::New();
        auto mMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mMapper->SetInputData(mLineData);
        auto mActor = vtkSmartPointer<vtkActor>::New();
        mActor->SetMapper(mMapper);
        mActor->GetProperty()->SetLineWidth(4);
        
        // Verschiedene Farben für verschiedene Muskeln
        if (m == 0) {
            mActor->GetProperty()->SetColor(0.0, 0.8, 1.0); // Cyan für Muskel 1
        } else if (m == 1) {
            mActor->GetProperty()->SetColor(1.0, 0.5, 0.0); // Orange für Muskel 2
        } else {
            // Weitere Farben für zusätzliche Muskeln
            double hue = (m * 137.5) / 360.0; // Goldener Winkel für gute Farbverteilung
            mActor->GetProperty()->SetColor(
                0.5 + 0.5 * std::cos(hue * 2 * M_PI),
                0.5 + 0.5 * std::cos((hue + 0.33) * 2 * M_PI),
                0.5 + 0.5 * std::cos((hue + 0.67) * 2 * M_PI)
            );
        }
        
        m_renderer->AddActor(mActor);
        m_muscleActors.push_back(mActor);
        m_muscleLines.push_back(mLineData);
        m_musclePointsData.push_back(mPoints);
    }

    // 3. Initial Guess Pipeline - FÜR MEHRERE MUSKELN
    m_initialLines.clear();
    m_initialPointsData.clear();
    m_initialPointActors.clear();

    for (size_t m = 0; m < m_initialGuesses.size(); ++m) {
        // --- LINIE für Initial Guess ---
        auto igLinePoints = vtkSmartPointer<vtkPoints>::New();
        auto igLineData = vtkSmartPointer<vtkPolyData>::New();
        auto igLineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        igLineMapper->SetInputData(igLineData);
        auto igLineActor = vtkSmartPointer<vtkActor>::New();
        igLineActor->SetMapper(igLineMapper);
        igLineActor->GetProperty()->SetColor(0.8, 0.8, 0.0); // Gelb
        igLineActor->GetProperty()->SetLineStipplePattern(0xf0f0);
        igLineActor->GetProperty()->SetOpacity(0.5);
        m_renderer->AddActor(igLineActor);
        
        m_initialLines.push_back(igLineData);
        m_initialPointsData.push_back(igLinePoints);

        // --- KUGELN für Initial Guess Punkte ---
        std::vector<vtkSmartPointer<vtkActor>> currentMuscleSphereActors;
        
        // Anzahl Punkte aus dem ersten Timestep ermitteln
        int numPts = m_initialGuesses[m].empty() ? 0 : m_initialGuesses[m][0].size();

        for (int i = 0; i < numPts; ++i) {
            auto sphere = vtkSmartPointer<vtkSphereSource>::New();
            sphere->SetRadius(0.01*m_scalerCM);
            sphere->SetThetaResolution(16);
            sphere->SetPhiResolution(16);
            
            auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapper->SetInputConnection(sphere->GetOutputPort());
            
            auto actor = vtkSmartPointer<vtkActor>::New();
            actor->SetMapper(mapper);
            actor->GetProperty()->SetColor(1.0, 0.6, 0.0); // Default Orange
            
            m_renderer->AddActor(actor);
            currentMuscleSphereActors.push_back(actor);
        }
        m_initialPointActors.push_back(currentMuscleSphereActors);
    }

    // other points with colors pipeline -> as actors with spheres
    m_otherPointActors.clear();
    for (size_t t = 0; t < m_otherPointsWithColors.size(); ++t) {
        std::vector<vtkSmartPointer<vtkActor>> timeStepActors;
        const auto& ptsAndColors = m_otherPointsWithColors[t];
        for (size_t i = 0; i + 1 < ptsAndColors.size(); i += 2) {
            const auto& p = ptsAndColors[i];
            const auto& c = ptsAndColors[i + 1];
            auto sphere = vtkSmartPointer<vtkSphereSource>::New();
            sphere->SetRadius(0.008 * m_scalerCM);
            sphere->SetThetaResolution(12);
            sphere->SetPhiResolution(12);
            auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapper->SetInputConnection(sphere->GetOutputPort());
            auto actor = vtkSmartPointer<vtkActor>::New();
            actor->SetMapper(mapper);
            actor->GetProperty()->SetColor(c.x, c.y, c.z);
            actor->SetPosition(p.x, p.y, p.z);
            m_renderer->AddActor(actor);
            timeStepActors.push_back(actor);
        }
        m_otherPointActors.push_back(timeStepActors);
    }

    /* if (ShowCoordinateSystems > -1)
    {    // ---------------------------------------------------------
        // 7. Mesh Coordinate Systems (Pro Mesh)
        // ---------------------------------------------------------
        m_meshAxesActors.clear();
        for (size_t i = 0; i < m_meshes.size(); ++i) {
            auto axes = vtkSmartPointer<vtkAxesActor>::New();
            axes->SetTotalLength(0.5, 0.5, 0.5); // Etwas kleiner als Welt-Achsen
            axes->SetShaftTypeToLine(); // Dünnere Linien für Meshes
            // axes->SetShaftTypeToCylinder(); // Oder dicker, wenn gewünscht
            
            // Optional: Beschriftung ausschalten, damit es nicht zu voll wird
            axes->SetAxisLabels(0); 

            m_renderer->AddActor(axes);
            m_meshAxesActors.push_back(axes);
        }
        
        // ---------------------------------------------------------
        // 6. World Coordinate System (Ursprung)
        // ---------------------------------------------------------
        m_worldAxes = vtkSmartPointer<vtkAxesActor>::New();
        m_worldAxes->SetTotalLength(1.0, 1.0, 1.0); // Länge der Achsen
        m_worldAxes->SetShaftTypeToCylinder();
        m_worldAxes->SetCylinderRadius(0.02);
        m_worldAxes->SetAxisLabels(1); // X, Y, Z Labels anzeigen
        m_renderer->AddActor(m_worldAxes);
    } */
    m_meshAxesActors.clear();
    for (size_t i = 0; i < m_meshes.size(); ++i) {
        auto axes = vtkSmartPointer<vtkAxesActor>::New();
        axes->SetTotalLength(0.5 * m_scalerCM, 0.5 * m_scalerCM, 0.5 * m_scalerCM); 
        axes->SetShaftTypeToLine(); 
        axes->SetAxisLabels(0); 

        // Initiale Sichtbarkeit basierend auf ShowCoordinateSystems
        axes->SetVisibility(ShowCoordinateSystems == 2); 

        m_renderer->AddActor(axes);
        m_meshAxesActors.push_back(axes);
    }
    m_worldAxes = vtkSmartPointer<vtkAxesActor>::New();
    m_worldAxes->SetTotalLength(1.0 * m_scalerCM, 1.0 * m_scalerCM, 1.0 * m_scalerCM); // Länge der Achsen
    m_worldAxes->SetShaftTypeToCylinder();
    m_worldAxes->SetCylinderRadius(0.02 * m_scalerCM);
    m_worldAxes->SetAxisLabels(1); // X, Y, Z Labels anzeigen
    m_renderer->AddActor(m_worldAxes);
    
    // ---------------------------------------------------------
    // 6. World Coordinate System (Ursprung) - IMMER ERSTELLEN
    // ---------------------------------------------------------
    /* m_worldAxes = vtkSmartPointer<vtkAxesActor>::New();
    m_worldAxes->SetTotalLength(1.0, 1.0, 1.0); 
    m_worldAxes->SetShaftTypeToCylinder();
    m_worldAxes->SetCylinderRadius(0.02);
    m_worldAxes->SetAxisLabels(1);
    m_renderer->AddActor(m_worldAxes);
    
    // Initiale Sichtbarkeit basierend auf ShowCoordinateSystems
    m_worldAxes->SetVisibility(ShowCoordinateSystems >= 1); */
    

    // 1. Vector zu einem String mit Zeilenumbrüchen (\n) zusammenfügen
    std::string combinedText = "";
    for (const auto& line : m_solverInfoText[0]) {
        combinedText += line + "\n";
    }
    m_textActor = vtkSmartPointer<vtkTextActor>::New();
    m_textActor->SetInput(combinedText.c_str()); // Text setzen
    m_textActor->GetTextProperty()->SetFontSize(14);
    m_textActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0); // Weiß
    m_textActor->GetTextProperty()->SetFontFamilyToArial();
    m_textActor->GetTextProperty()->SetLineSpacing(1.2); 
    m_textActor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    m_textActor->SetPosition(0.98, 0.98); // Kleiner Abstand zum Rand
    m_textActor->GetTextProperty()->SetJustificationToRight();      // Rechtsbündig
    m_textActor->GetTextProperty()->SetVerticalJustificationToTop(); // Obenbündig
    m_renderer->AddActor2D(m_textActor);
    
    qDebug() << "Pipeline setup complete:" << m_muscleResults.size() << "muscles," 
             << m_meshes.size() << "meshes," << m_initialGuesses.size() << "initial guess paths";
}

void VTKSimViewerSimple::closeEvent(QCloseEvent *event)
{
    // Frage den User
    QMessageBox::StandardButton reply;
    reply = QMessageBox::question(this, "Exportieren?", 
                                  "Möchten Sie die Simulation vor dem Beenden als Bilder speichern?",
                                  QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel);

    if (reply == QMessageBox::Yes) {
        // 1. Export durchführen
        // (Achtung: Fenster muss sichtbar bleiben währenddessen)
        std::string exportPath = "../examples/results/resultPictures"; 
        this->exportAllStepsToImages(exportPath);
        
        // 2. Schließen akzeptieren
        event->accept();
    } 
    else if (reply == QMessageBox::No) {
        // Einfach schließen ohne Export
        event->accept();
    } 
    else {
        // "Abbrechen" gedrückt -> Fenster offen lassen
        event->ignore();
    }
}

void VTKSimViewerSimple::exportAllStepsToImages(std::string outputFolder)
{
    // 1. Ordner erstellen
    QDir dir(QString::fromStdString(outputFolder));
    if (!dir.exists()) {
        dir.mkpath(".");
    }

    qDebug() << "Starte automatischen Bild-Export nach:" << QString::fromStdString(outputFolder);

    int totalSteps = this->m_angles.size(); 

    for (int i = 0; i < totalSteps; ++i) {
        
        // A. Szene updaten (KORREKTUR: updateStep statt updateScene)
        this->updateStep(i); 

        // B. Render erzwingen (KORREKTUR: Zugriff über m_vtkWidget)
        this->m_vtkWidget->renderWindow()->Render(); 
        
        // C. Screenshot
        vtkNew<vtkWindowToImageFilter> windowToImageFilter;
        // KORREKTUR: Zugriff über m_vtkWidget
        windowToImageFilter->SetInput(this->m_vtkWidget->renderWindow());
        windowToImageFilter->ReadFrontBufferOff(); 
        windowToImageFilter->Update();

        // D. Speichern
        QString filename = QString("step%1.png").arg(i);
        QString fullPath = dir.filePath(filename);

        vtkNew<vtkPNGWriter> writer;
        writer->SetFileName(fullPath.toStdString().c_str());
        writer->SetInputConnection(windowToImageFilter->GetOutputPort());
        writer->Write();

        qDebug() << "Exportiert:" << filename;
    }
    qDebug() << "Export abgeschlossen.";
}

void VTKSimViewerSimple::updateMeshTransform(vtkProp3D* actor, const MWMath::RotMatrix3x3& R, const MWMath::Point3D& pos) {
    vtkSmartPointer<vtkMatrix4x4> vtkMat = vtkSmartPointer<vtkMatrix4x4>::New();
    vtkMat->Identity();

    // Translation
    vtkMat->SetElement(0, 3, pos.x);
    vtkMat->SetElement(1, 3, pos.y);
    vtkMat->SetElement(2, 3, pos.z);

    // Rotation
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            vtkMat->SetElement(i, j, R.m[i][j]);
        }
    }
    actor->SetUserMatrix(vtkMat);
}

void VTKSimViewerSimple::updateStep(int step) {

    if (m_stepTextActor) {
        std::string stepStr = "Step: " + std::to_string(step) + " / " + std::to_string(m_angles.size() - 1).c_str();
        m_stepTextActor->SetInput(stepStr.c_str());
    }

    // 1. Mesh Transformation
    for (size_t i = 0; i < m_meshActors.size(); ++i) {
        if (step < (int)m_meshes[i]->allRMatrixGlobal.size()) {
            updateMeshTransform(m_meshActors[i], 
                              m_meshes[i]->allRMatrixGlobal[step], 
                              m_meshes[i]->MeshPointsGlobal[step]);
        }
        if (!m_meshes[i]->AllScalerStepLists.empty()) {
            
            // Sicherstellen, dass wir nicht out-of-bounds lesen
            if (step < (int)m_meshes[i]->AllScalerStepLists.size()) {
                
                // Skalierungsvektor holen (x, y, z)
                MWMath::Point3D scale = m_meshes[i]->AllScalerStepLists[step];
                
                // Auf den Actor anwenden
                m_meshActors[i]->SetScale(scale.x, scale.y, scale.z);
            }
            else {
                // Fallback: Wenn Liste kürzer als Schritte -> Letzten Wert nehmen?
                // Oder einfach nichts tun (behält letzte Skalierung).
                // Hier Beispiel für letzten Wert:
                MWMath::Point3D scale = m_meshes[i]->AllScalerStepLists.back();
                m_meshActors[i]->SetScale(scale.x, scale.y, scale.z);
            }
        }

        // --- NEU: Achsenkreuz für dieses Mesh updaten ---
        /* if (ShowCoordinateSystems && i < m_meshAxesActors.size()) {
            if (step < (int)m_meshes[i]->allRMatrixGlobal.size()) {
                // Wir nutzen denselben Helper wie für das Mesh!
                // Das Achsenkreuz bewegt sich exakt mit dem Mesh mit.
                updateMeshTransform(m_meshAxesActors[i], 
                                  m_meshes[i]->allRMatrixGlobal[step], 
                                  m_meshes[i]->MeshPointsGlobal[step]);
            }
        } */
        if (i < m_meshAxesActors.size() && m_meshAxesActors[i]) {
            if (step < (int)m_meshes[i]->allRMatrixGlobal.size()) {
                updateMeshTransform(m_meshAxesActors[i], 
                                m_meshes[i]->allRMatrixGlobal[step], 
                                m_meshes[i]->MeshPointsGlobal[step]);
            }
        }
    }

    // 2. Muskel Update (Optimierter Pfad) - FÜR ALLE MUSKELN
    for (size_t m = 0; m < m_muscleResults.size(); ++m) {
        if (step >= (int)m_muscleResults[m].size()) continue;
        const auto& path = m_muscleResults[m][step];

        m_musclePointsData[m]->Reset();
        auto cells = vtkSmartPointer<vtkCellArray>::New();
        auto polyLine = vtkSmartPointer<vtkPolyLine>::New();
        polyLine->GetPointIds()->SetNumberOfIds(path.size());

        for (size_t i = 0; i < path.size(); ++i) {
            m_musclePointsData[m]->InsertNextPoint(path[i].x, path[i].y, path[i].z);
            polyLine->GetPointIds()->SetId(i, i);
        }
        cells->InsertNextCell(polyLine);
        m_muscleLines[m]->SetPoints(m_musclePointsData[m]);
        m_muscleLines[m]->SetLines(cells);
        m_muscleLines[m]->Modified();
    }

    // 3. Initial Guess Update (Linie) - FÜR ALLE MUSKELN
    for (size_t m = 0; m < m_initialGuesses.size(); ++m) {
        if (step >= (int)m_initialGuesses[m].size()) continue;
        
        const auto& igPath = m_initialGuesses[m][step];

        m_initialPointsData[m]->Reset();
        auto cells = vtkSmartPointer<vtkCellArray>::New();
        auto polyLine = vtkSmartPointer<vtkPolyLine>::New();
        polyLine->GetPointIds()->SetNumberOfIds(igPath.size());

        for (size_t i = 0; i < igPath.size(); ++i) {
            m_initialPointsData[m]->InsertNextPoint(igPath[i].x, igPath[i].y, igPath[i].z);
            polyLine->GetPointIds()->SetId(i, i);
        }
        
        cells->InsertNextCell(polyLine);
        m_initialLines[m]->SetPoints(m_initialPointsData[m]);
        m_initialLines[m]->SetLines(cells);
        m_initialLines[m]->Modified();
    }

    // 4. Initial Guess Punkte (Kugeln mit Farben) - FÜR ALLE MUSKELN
    for (size_t m = 0; m < m_initialGuesses.size(); ++m) {
        if (step >= (int)m_initialGuesses[m].size()) continue;

        const auto& path = m_initialGuesses[m][step];
        const auto& colors = m_initialGuessColors[m][step];
        auto& actors = m_initialPointActors[m];

        for (size_t i = 0; i < path.size() && i < actors.size(); ++i) {
            // Position setzen
            actors[i]->SetPosition(path[i].x, path[i].y, path[i].z);

            // Farbe setzen
            if (i < colors.size()) {
                actors[i]->GetProperty()->SetColor(colors[i].x, colors[i].y, colors[i].z);
            } else {
                actors[i]->GetProperty()->SetColor(1.0, 0.6, 0.0); // Fallback Orange
            }
            actors[i]->SetVisibility(1);
        }
        
        // Verstecke überschüssige Aktoren (falls path kürzer als actors)
        for (size_t i = path.size(); i < actors.size(); ++i) {
            actors[i]->SetVisibility(0);
        }
    } 
    // 4. Initial Guess / Berechnete Punkte (Kugeln mit Farben) - FÜR ALLE MUSKELN
    /* for (size_t m = 0; m < m_muscles.size(); ++m) {
        // Zugriff auf den aktuellen Muskel
        SSMuscle* mus = m_muscles[m]; 
        
        // Hole die Aktoren für diesen Muskel (Kugeln)
        // (Stelle sicher, dass m_initialPointActors in setupVtkPipeline groß genug initialisiert wurde!)
        if (m >= m_initialPointActors.size()) continue;
        auto& actors = m_initialPointActors[m];

        // Iteriere über alle Knoten des Muskels (MNodes)
        // Wir nehmen an: MNodes[0] bis MNodes[N-1]
        for (size_t i = 0; i < mus->MNodes.size(); ++i) {
            
            // Safety Check: Haben wir genug Aktoren?
            if (i >= actors.size()) break;

            MuscleNode& node = mus->MNodes[i];

            // Safety Check: Existieren Daten für diesen Step?
            bool hasData = (step < (int)node.MNodeLocalSteps.size()) && 
                           (step < (int)node.MNodeParentPosSteps.size()) && 
                           (step < (int)node.MNodeParentRotSteps.size());

            if (hasData) {
                // DATEN HOLEN
                MWMath::Point3D localPos = node.MNodeLocalSteps[step];
                MWMath::Point3D parentPos = node.MNodeParentPosSteps[step];
                MWMath::RotMatrix3x3 parentRot = node.MNodeParentRotSteps[step];

                // BERECHNUNG: Global = ParentPos + ParentRot * LocalPos
                MWMath::Point3D globalPos = parentPos + parentRot.transform(localPos);

                // POSITION SETZEN
                actors[i]->SetPosition(globalPos.x, globalPos.y, globalPos.z);
                
                // FARBE SETZEN (Falls vorhanden)
                if (step < (int)node.MNodeInitialGuessColorSteps.size()) {
                    MWMath::Point3D col = node.MNodeInitialGuessColorSteps[step];
                    actors[i]->GetProperty()->SetColor(col.x, col.y, col.z);
                } else {
                    actors[i]->GetProperty()->SetColor(1.0, 0.6, 0.0); // Fallback Orange
                }

                actors[i]->SetVisibility(1);
            } 
            else {
                // Keine Daten für diesen Step -> Verstecken
                actors[i]->SetVisibility(0);
            }
        }
        
        // Verstecke restliche Aktoren (falls der Muskel in diesem Frame weniger Knoten hätte, 
        // oder Actors-Pool größer ist als Node-Anzahl)
        for (size_t i = mus->MNodes.size(); i < actors.size(); ++i) {
            actors[i]->SetVisibility(0);
        }
    }
 */
    // 5. Other Points with Colors Update
    if (step < (int)m_otherPointActors.size()) {
        const auto& actors = m_otherPointActors[step];
        for (const auto& actor : actors) {
            actor->SetVisibility(1);
        }
    }
    for (size_t t = 0; t < m_otherPointActors.size(); ++t) {
        if (t != (size_t)step) {
            const auto& actors = m_otherPointActors[t];
            for (const auto& actor : actors) {
                actor->SetVisibility(0);
            }
        }
    }


    // 6. Solver Info Text Update
    if (m_textActor && step < (int)m_solverInfoText.size()) {
        
        std::string combinedText = "";
        double r = 1.0;
        double g = 1.0;
        double b = 1.0;
        for (const auto& line : m_solverInfoText[step]) {
            combinedText += line + "\n";
            if (line.find("Solve_Succeeded") != std::string::npos) {
                r = 0.6; g = 1.0; b = 0.6; 
            }
            else if (line.find("Infeasible") != std::string::npos || 
                     line.find("Failed") != std::string::npos) {
                // Leichtes Rot / Pastellrot
                r = 1.0; g = 0.5; b = 0.5; 
            }
        }

        m_textActor->SetInput(combinedText.c_str());
        m_textActor->GetTextProperty()->SetColor(r, g, b);
    }
    
    m_vtkWidget->renderWindow()->Render();
}

void VTKSimViewerSimple::keyPressEvent(QKeyEvent *event)
{
    if (event->key() == Qt::Key_C) {
        toggleCoordinateSystems();
    } 
    else if (event->key() == Qt::Key_Right) {
        // Einen Schritt vorwärts (Slider-Maximum beachten)
        if (m_slider->value() < m_slider->maximum()) {
            m_slider->setValue(m_slider->value() + 1);
        }
    } 
    else if (event->key() == Qt::Key_Left) {
        // Einen Schritt zurück (Slider-Minimum beachten)
        if (m_slider->value() > m_slider->minimum()) {
            m_slider->setValue(m_slider->value() - 1);
        }
    }
    else {
        // Andere Events an die Basisklasse weiterleiten
        QDialog::keyPressEvent(event);
    }
}

void VTKSimViewerSimple::toggleCoordinateSystems()
{
    // Modus weiterschalten: 0 -> 1 -> 2 -> 0
    ShowCoordinateSystems = (ShowCoordinateSystems + 1) % 3;

    // Feedback im Debug
    QString modeStr;
    if (ShowCoordinateSystems == 0) modeStr = "NONE";
    else if (ShowCoordinateSystems == 1) modeStr = "WORLD ONLY";
    else modeStr = "ALL (WORLD + MESHES)";
    //qDebug() << "Coordinate Systems Mode:" << modeStr;

    // --- 1. Welt-Achsenkreuz ---
    if (m_worldAxes) {
        // Sichtbar wenn Modus >= 1
        m_worldAxes->SetVisibility(ShowCoordinateSystems >= 1);
    }

    // --- 2. Mesh-Achsenkreuze ---
    for (auto& axes : m_meshAxesActors) {
        if (axes) {
            // Sichtbar nur wenn Modus == 2 (Alles)
            axes->SetVisibility(ShowCoordinateSystems == 2);
        }
    }

    // Render erzwingen, damit man es sofort sieht
    m_vtkWidget->renderWindow()->Render();
}

void VTKSimViewerSimple::generatePlots()
{
    std::cout << "\n[Python] Starte Plot-Generierung..." << std::endl;
    
    // Relativer Pfad zur Python-Executable in der venv
    // (Geht von 'build' einen Ordner hoch zu '.venv')
    std::string pythonExec = "../.venv/bin/python3";
    std::string pythonScript = "../src/pythonScripts/plot_muscle_logs.py";
    
    std::string pythonCommand = pythonExec + " " + pythonScript;
    
    int retCode = std::system(pythonCommand.c_str());
    
    if (retCode == 0) {
        std::cout << "[Python] PDFs erfolgreich erstellt!" << std::endl;
    } else {
        std::cerr << "[Python] Fehler beim Ausführen des Skripts (Code: " << retCode << ")." << std::endl;
        std::cerr << "         Stelle sicher, dass die .venv existiert und numpy/matplotlib installiert sind." << std::endl;
    }
}


