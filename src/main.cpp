#include <QApplication>
#include "utils/rawSimulation.h"
#include "utils/ConfigLoader.h"
#include "utils/utility.h"
#include <chrono>
// test
#include "simpleSimulation/casadiSystem.h"
#include "simpleSimulation/VTKSimViewerSimple.h"
#include "simpleSimulation/SSMuscle.h"
#include "simpleSimulation/SSMeshes.h"
#include "simpleSimulation/SSBody.h"



std::vector<double> MAXJOINTANGLES = {90.0, 100.0, 90.0}; // MCP, PIP, DIP
MWMath::RotMatrix3x3 FINGERSTARTORIENTATION = MWMath::axisAngle(MWMath::Point3D(0,1,0), 90.0); // Startorientierung des Fingermodells
MWMath::Point3D ROTAXIS = MWMath::Point3D(0,1,0).normed(); // Rotationsachse für Fingerbeugung
double ROTENDANBGLE = 90.0; // Endwinkel der Beugung

double PUSHA = 0.25, PUSHB = 0.4, PUSHC = 0.7; // Halbachsen des Pushers (Kugel)
MWMath::Point3D MOVEDIR = MWMath::Point3D(0, -1, 0).normed(); // Verschiebung des Fingermodells beim Drücken
MWMath::Point3D STARTP = MWMath::Point3D(0.0, 0.3, -0.0);//MWMath::Point3D(0.0, 0.5, -0.5); // Startpunkt des Drückers
double MOVEDIST = 1.0; // Distanz, die der Drücker bewegt wird


void updateSceneMovement(std::string sceneName, std::vector<std::shared_ptr<SSMesh>>& meshes, double progress){
    if (sceneName == "ELLIPSOID_PUSH_THROUGH" || sceneName == "TORUS_PUSH_THROUGH" || sceneName == "CYLINDER_PUSH_THROUGH") {
        
        if (meshes.size() < 3) return;

        //double startY = 0.5;
        //double endY = -1.5;
        double endRot = ROTENDANBGLE;
        //double currentY = startY + (endY - startY) * progress;
        double currentRot = 0.0 + (endRot - 0.0) * progress;

        // Wir nehmen an, der Pusher ist das letzte Mesh oder Index 2
        // (Reihenfolge in setupScene: Left(0), Right(1), Pusher(2))
        auto pusherMesh = meshes[2]; 
        
        if (pusherMesh->Parent) {
            pusherMesh->Parent->PositionGlobal = STARTP + MOVEDIR * MOVEDIST * progress;
            pusherMesh->PositionGlobal = pusherMesh->Parent->PositionGlobal; 
            pusherMesh->Parent->OrientationGlobal = FINGERSTARTORIENTATION * MWMath::axisAngle(ROTAXIS, currentRot);
            
            
            // WICHTIG beim Torus: Wenn er rotiert ist (Orientation2ParentRel),
            // muss OrientationGlobal auch geupdatet werden!
            // GlobalRot = ParentRot * LocalRot
            pusherMesh->OrientationGlobal = pusherMesh->Parent->OrientationGlobal * pusherMesh->Orientation2ParentRel;
            qDebug() << "Pusher Determinant: " << pusherMesh->OrientationGlobal.determinant();
        } else {
            pusherMesh->PositionGlobal = STARTP + MOVEDIR * MOVEDIST * progress;
        }
    }
    else if (sceneName == "TORUS_ROTATE_THROUGH") {
        
        if (meshes.size() < 3) return;

        double startY = 0.0;
        double startRot = 45.0;
        double endY = -180.0;
        double rotation = startRot + (endY - startRot) * progress;

        // Wir nehmen an, der Pusher ist das letzte Mesh oder Index 2
        // (Reihenfolge in setupScene: Left(0), Right(1), Pusher(2))
        auto pusherMesh = meshes[2]; 
        
        if (pusherMesh->Parent) {
            //pusherMesh->Parent->PositionGlobal.y = currentY;
            pusherMesh->PositionGlobal.y = 0.2; 
            
            // WICHTIG beim Torus: Wenn er rotiert ist (Orientation2ParentRel),
            // muss OrientationGlobal auch geupdatet werden!
            // GlobalRot = ParentRot * LocalRot
            pusherMesh->Parent->OrientationGlobal = MWMath::axisAngle({0, 1, 0}, rotation);
            pusherMesh->OrientationGlobal = pusherMesh->Parent->OrientationGlobal * pusherMesh->Orientation2ParentRel;
        } else {
            pusherMesh->OrientationGlobal = MWMath::axisAngle({0, 1, 0}, rotation) * pusherMesh->Orientation2ParentRel;
        }
    }
    else if (sceneName == "TWO_SEGMENTS_SIMPLE") {
        if (meshes.size() < 3) return; // Seg1, Seg2, Torus

        // --- Längen holen ---
        // Annahme: meshes[0] ist Seg1, meshes[1] ist Seg2
        double lenSeg1 = std::dynamic_pointer_cast<SSEllipsoidMesh>(meshes[0])->C;
        double lenSeg2 = std::dynamic_pointer_cast<SSEllipsoidMesh>(meshes[1])->C;

        // --- Winkel ---
        double angle = 90.0 * progress;
        MWMath::Point3D rotAxis = {0, 0, -1}; // Drehung um Z-Achse

        // --- Orientierung ---
        MWMath::RotMatrix3x3 rot1 = MWMath::RotMatrix3x3(); // Identity (X-Achse)
        MWMath::RotMatrix3x3 rot2 = rot1 * MWMath::axisAngle(rotAxis, angle); // Gedreht

        // ----------------------------------------------------------
        // 1. SEGMENT 1 (Proximal)
        // ----------------------------------------------------------
        // Wir wollen, dass das Segment bei STARTP *anfängt*.
        // Der Body ist aber in der *Mitte* (CoM).
        // Also: BodyPos = Start + (Länge/2 entlang X)
        MWMath::Point3D posBody1 = STARTP + rot1 * MWMath::Point3D(lenSeg1 / 2.0, 0, 0);

        auto body1 = meshes[0]->Parent;
        if(body1) {
            body1->PositionGlobal = posBody1;
            body1->OrientationGlobal = rot1;
        }

        // ----------------------------------------------------------
        // 2. GELENK PUNKT (Joint)
        // ----------------------------------------------------------
        // Das Gelenk ist am Ende von Seg1.
        // Vom Body1 (Mitte) aus gesehen ist das +Länge/2.
        MWMath::Point3D posJoint = posBody1 + rot1 * MWMath::Point3D(lenSeg1, 0, 0);

        // ----------------------------------------------------------
        // 3. SEGMENT 2 (Distal)
        // ----------------------------------------------------------
        // Der Body2 (Mitte) muss so liegen, dass das Mesh am Joint anfängt.
        // Also: Body2 = Joint + (Länge2/2 entlang der NEUEN Rotation)
        MWMath::Point3D posBody2 = posJoint + rot2 * MWMath::Point3D(lenSeg2, 0, 0);

        auto body2 = meshes[1]->Parent;
        if(body2) {
            body2->PositionGlobal = posBody2;
            body2->OrientationGlobal = rot2;
        }

        // ----------------------------------------------------------
        // SYNC MESHES
        // ----------------------------------------------------------
        for(auto& m : meshes) {
            if(m->Parent) {
                // WICHTIG: Das funktioniert nur, wenn im Setup der Offset (0,0,0) ist!
                // (Außer beim Torus, der hat einen Offset)
                m->OrientationGlobal = m->Parent->OrientationGlobal * m->Orientation2ParentRel;
                m->PositionGlobal = m->Parent->PositionGlobal + m->Parent->OrientationGlobal * m->Position2ParentRelInParentFrame;
            }
        }
    }
}


void setupSceneObjectOriented(std::string currentScene, std::vector<std::shared_ptr<SSTissue>>& tissues, std::vector<std::shared_ptr<SSMesh>>& meshes, std::vector<SSMuscle*>& muscles, std::shared_ptr<SSBody>& rootSystem,int numTimeSteps, const SimSettings& cfg) 
{
    meshes.clear();
    muscles.clear();

    if (currentScene == "NONE"){
        bool bFingerIsEllipsoid = cfg.bFingerIsEllipsoid;
        bool bCreateJointMeshes = cfg.bCreateJointMeshes;
        
        // --- PARAMETER ---
        double handLength = 1.8;
        // Längen der Fingerglieder (Halbachsen C): MC, PP, MP, DP
        std::vector<double> fLength = {0.3482 * handLength, 0.2027 * handLength, 0.1175 * handLength, 0.0882 * handLength};
        double dimFing = 0.2; // Dicken-Skalierung (A und B)
        double jRadFactor = 1.1; // Gelenk-Radius Faktor relativ zur Knochendicke
        
        // Hilfsvariablen für Geometrie
        double widthMC = dimFing * fLength[0];
        double widthPP = dimFing * fLength[1];
        double widthMP = dimFing * fLength[2];
        double widthDP = dimFing * fLength[3];

        // ==============================================================================
        // 1. SKELETT-STRUKTUR AUFBAUEN (Bodies & Joints)
        // ==============================================================================

        // Root
        rootSystem = std::make_shared<SSBody>("HandRoot", MWMath::Point3D(0,0,0));
        tissues.push_back(rootSystem);

        // --- MC (Metacarpal) ---
        // MC startet bei Root (0,0,0) oder verschoben? Sagen wir MC-Mitte ist bei C.
        auto bodyMC = std::make_shared<SSBody>("Body_MC", MWMath::Point3D(0,0, fLength[0]), MWMath::RotMatrix3x3(), rootSystem);
        tissues.push_back(bodyMC);
        
        // --- Joint MCP ---
        // Joint ist am Ende von MC. BodyMC ist in der Mitte (bei C). Also Joint ist bei +C relativ zu BodyMC.
        auto jointMCP = std::make_shared<SSJoint>("Joint_MCP", 
                                                MWMath::Point3D(0, 0, fLength[0]), // Relativ zu BodyMC (Mitte) -> Ende
                                                MWMath::RotMatrix3x3(), 
                                                bodyMC, 
                                                90.0, 
                                                MWMath::Point3D(1,0,0), 
                                                numTimeSteps);
        tissues.push_back(jointMCP);

        // --- PP (Proximal Phalanx) ---
        // BodyPP-Mitte soll Abstand C zum Joint haben (damit er danach kommt)
        auto bodyPP = std::make_shared<SSBody>("Body_PP", MWMath::Point3D(0,0, fLength[1]), MWMath::RotMatrix3x3(), jointMCP);
        tissues.push_back(bodyPP);

        // --- Joint PIP ---
        // Joint ist am Ende von PP. Relativ zur Mitte PP ist das +C.
        auto jointPIP = std::make_shared<SSJoint>("Joint_PIP", 
                                                MWMath::Point3D(0, 0, fLength[1]), 
                                                MWMath::RotMatrix3x3(), 
                                                bodyPP, 
                                                100.0, 
                                                MWMath::Point3D(1,0,0), 
                                                numTimeSteps);
        tissues.push_back(jointPIP);

        // --- MP (Middle Phalanx) ---
        // Mitte MP ist +C vom Joint entfernt
        auto bodyMP = std::make_shared<SSBody>("Body_MP", MWMath::Point3D(0,0, fLength[2]), MWMath::RotMatrix3x3(), jointPIP);
        tissues.push_back(bodyMP);

        // --- Joint DIP ---
        // Joint ist am Ende von MP (+C)
        auto jointDIP = std::make_shared<SSJoint>("Joint_DIP", 
                                                MWMath::Point3D(0, 0, fLength[2]), 
                                                MWMath::RotMatrix3x3(), 
                                                bodyMP, 
                                                90.0, 
                                                MWMath::Point3D(1,0,0), 
                                                numTimeSteps);
        tissues.push_back(jointDIP);

        // --- DP (Distal Phalanx) ---
        // Mitte DP ist +C vom Joint entfernt
        auto bodyDP = std::make_shared<SSBody>("Body_DP", MWMath::Point3D(0,0, fLength[3]), MWMath::RotMatrix3x3(), jointDIP);
        tissues.push_back(bodyDP);
        

        // ==============================================================================
        // 2. MESHES HINZUFÜGEN
        // ==============================================================================

        // Helper Lambda zum Hinzufügen (spart Code)
        auto addMesh = [&](std::shared_ptr<SSBody> body, std::shared_ptr<SSMesh> mesh) {
            body->Meshes.push_back(mesh);
            mesh->Parent = body;
            meshes.push_back(mesh);
        };
        auto addJointMesh = [&](std::shared_ptr<SSJoint> joint, std::shared_ptr<SSMesh> mesh) {
            joint->Meshes.push_back(mesh);
            mesh->Parent = joint;
            meshes.push_back(mesh);
        };

        // --- A) KNOCHEN (Ellipsoide) ---
        std::shared_ptr<SSEllipsoidMesh> meshMC, meshPP, meshMP, meshDP; // Pointer merken für Muskeln

        if (bFingerIsEllipsoid) {
            // MC
            meshMC = std::make_shared<SSEllipsoidMesh>(widthMC, widthMC, fLength[0]);
            meshMC->Name = "Mesh_MC"; meshMC->MeshColor = {1,0,0};
            // Mesh ist jetzt ZENTRIERT im Body (weil Body-Ursprung = Mitte)
            meshMC->Position2ParentRelInParentFrame = MWMath::Point3D(0, 0, 0); 
            addMesh(bodyMC, meshMC);

            // PP
            meshPP = std::make_shared<SSEllipsoidMesh>(widthPP, widthPP, fLength[1]);
            meshPP->Name = "Mesh_PP"; meshPP->MeshColor = {0,0,1};
            meshPP->Position2ParentRelInParentFrame = MWMath::Point3D(0, 0, 0); // Zentriert
            addMesh(bodyPP, meshPP);

            // MP
            meshMP = std::make_shared<SSEllipsoidMesh>(widthMP, widthMP, fLength[2]);
            meshMP->Name = "Mesh_MP"; meshMP->MeshColor = {0,1,0};
            meshMP->Position2ParentRelInParentFrame = MWMath::Point3D(0, 0, 0); // Zentriert
            addMesh(bodyMP, meshMP);

            // DP
            meshDP = std::make_shared<SSEllipsoidMesh>(widthDP, widthDP, fLength[3]);
            meshDP->Name = "Mesh_DP"; meshDP->MeshColor = {1,1,0};
            meshDP->Position2ParentRelInParentFrame = MWMath::Point3D(0, 0, 0); // Zentriert
            addMesh(bodyDP, meshDP);
        }
        else {
        }

        // --- B) RINGBÄNDER (Tori) ---
        auto t2 = std::make_shared<SSTorusMesh>(widthPP * 1.8, 0.10);
        if (currentScene == "NONE") {
            double torusMinor = 0.10;
            double tiltDeg = 0.0;

            // Torus 1 (an MC)
            /* auto t1 = std::make_shared<SSTorusMesh>(widthMC * 1.8, torusMinor);
            t1->Parent = bodyMC;
            t1->Name = "Torus_MC"; t1->MeshColor = {0,1,1};
            t1->Position2ParentRelInParentFrame = MWMath::Point3D(0, -widthMC * 1.5, fLength[0] * 0.3);
            t1->Orientation2ParentRel = MWMath::axisAngle({1,0,0}, tiltDeg);
            //t1->C = (meshMC ? meshMC->C : fLength[0]) * 7; t1->A = 0.5;
            addMesh(bodyMC, t1); */

            // Torus 2 (an PP)
            //auto t2 = std::make_shared<SSTorusMesh>(widthPP * 1.8, torusMinor);
            t2->Parent = bodyPP;
            t2->Name = "Torus_PP"; t2->MeshColor = {1,0,1};
            t2->Position2ParentRelInParentFrame = MWMath::Point3D(0, -widthPP * 1.3, fLength[1] * 0.3);
            t2->Orientation2ParentRel = MWMath::axisAngle({1,0,0}, tiltDeg);
            /* t2->C = (meshPP ? meshPP->C : fLength[1]) * 8; t2->A = 0.5; */
            addMesh(bodyPP, t2);

            // Torus 3 (an MP)
            /* auto t3 = std::make_shared<SSTorusMesh>(widthMP * 1.8, torusMinor);
            t3->Parent = bodyMP;
            t3->Name = "Torus_MP"; t3->MeshColor = {1,0.5,0};
            t3->Position2ParentRelInParentFrame = MWMath::Point3D(0, -widthMP * 1.5, fLength[2] * 0.3);
            t3->Orientation2ParentRel = MWMath::axisAngle({1,0,0}, tiltDeg);
            // t3->C = (meshMP ? meshMP->C : fLength[2]) * 12; t3->A = 0.5;
            addMesh(bodyMP, t3); */
        }
   
        // --- C) GELENK-MESHES (Kugeln) ---
        if (bCreateJointMeshes) {
            /* // MCP
            double r1 = widthMC * jRadFactor;
            auto mJ1 = std::make_shared<SSEllipsoidMesh>(r1, r1, r1);
            mJ1->Name = "JointMesh_MCP"; mJ1->MeshColor = {0.5, 1, 0.5}; mJ1->bIsJointMesh = true;
            addJointMesh(jointMCP, mJ1);

            // PIP
            double r2 = widthPP * jRadFactor;
            auto mJ2 = std::make_shared<SSEllipsoidMesh>(r2, r2, r2);
            mJ2->Name = "JointMesh_PIP"; mJ2->MeshColor = {0.5, 1, 0.5}; mJ2->bIsJointMesh = true;
            addJointMesh(jointPIP, mJ2);

            // DIP
            double r3 = widthMP * jRadFactor;
            auto mJ3 = std::make_shared<SSEllipsoidMesh>(r3, r3, r3);
            mJ3->Name = "JointMesh_DIP"; mJ3->MeshColor = {0.5, 1, 0.5}; mJ3->bIsJointMesh = true;
            addJointMesh(jointDIP, mJ3); */

            // MCP
            double r1 = widthMC * jRadFactor;
            auto mJ1 = std::make_shared<SSCylinderMesh>(r1, r1);
            mJ1->Name = "JointMesh_MCP"; mJ1->MeshColor = {0.5, 1, 0.5}; mJ1->bIsJointMesh = true;
            mJ1->Orientation2ParentRel = MWMath::axisAngle({0,1,0}, 90.0); // Zylinder um X-Achse drehen, damit er wie eine Kugel wirkt
            addJointMesh(jointMCP, mJ1);

            // PIP
            double r2 = widthPP * jRadFactor;
            auto mJ2 = std::make_shared<SSCylinderMesh>(r2, r2);
            mJ2->Name = "JointMesh_PIP"; mJ2->MeshColor = {0.5, 1, 0.5}; mJ2->bIsJointMesh = true;
            mJ2->Orientation2ParentRel = MWMath::axisAngle({0,1,0}, 90.0); // Zylinder um X-Achse drehen, damit er wie eine Kugel wirkt
            addJointMesh(jointPIP, mJ2);

            // DIP
            double r3 = widthMP * jRadFactor;
            auto mJ3 = std::make_shared<SSCylinderMesh>(r3, r3);
            mJ3->Name = "JointMesh_DIP"; mJ3->MeshColor = {0.5, 1, 0.5}; mJ3->bIsJointMesh = true;
            mJ3->Orientation2ParentRel = MWMath::axisAngle({0,1,0}, 90.0); // Zylinder um X-Achse drehen, damit er wie eine Kugel wirkt
            addJointMesh(jointDIP, mJ3);
        }


        // ==============================================================================
        // INITIALES UPDATE (Wichtig für Muskel-Punkte!)
        // ==============================================================================
        rootSystem->update(0); 


        // ==============================================================================
        // MUSKELN
        // ==============================================================================
        
        // Wir brauchen Pointer auf die Meshes für den Muskel.
        // Da wir shared_ptr haben, holen wir uns die Raw Pointers.
        std::vector<SSMesh*> muscleMeshPtrs;
        
        // NUR relevante Meshes hinzufügen (Knochen + Tori). Gelenke ignorieren?
        // In deinem alten Code hast du alles hinzugefügt.
        /* for(auto& m : meshes) {
            muscleMeshPtrs.push_back(m.get());
        }  */
       muscleMeshPtrs.push_back(t2.get());

        // FLEXOR
        int numPoints = (!cfg.muscleNumPoints.empty()) ? cfg.muscleNumPoints[0] : 25;
        static SSMuscle flexor("Flexor", numPoints, 
            meshMC.get(), {0.0, -(widthMC * 1.2), 0.0}, 
            meshDP.get(), {0.0, -(widthDP * 1.2), 0.0});
            // meshMP.get(), {0.0, -(widthMP * 1.2), 0.0}); // Skalierung angepasst
        
        flexor.meshPtrs = muscleMeshPtrs;
        flexor.createMusclePoints();
        muscles.push_back(&flexor);

        // EXTENSOR (Optional)
        /* static SSMuscle extensor("Extensor", numPoints, 
            meshMC.get(), {0.0, boneWidthMC * 1.1, 0.0}, 
            meshDP.get(), {0.0, boneWidthMP * 1.1, 0.0});
        
        // Extensor braucht keine Tori als Hindernis (evtl filtern?)
        std::vector<SSMesh*> extensorMeshPtrs;
        for(auto& m : meshes) {
            if(!dynamic_cast<SSTorusMesh*>(m.get())) extensorMeshPtrs.push_back(m.get());
        }
        extensor.meshPtrs = extensorMeshPtrs;
        extensor.createMusclePoints();
        muscles.push_back(&extensor); */
    }
    else if (currentScene == "ELLIPSOID_PUSH_THROUGH") {
        
        // --- PARAMETER ---
        double anchorRadius = 0.2;
        double pusherRadius = 0.3;
        double anchorDist = 2.0; // Abstand der Anker vom Zentrum (X-Achse)
        double pushStart = 0.5;  // Startposition des Pushers (Y-Achse)
        
        // ==============================================================================
        // 1. STRUKTUR AUFBAUEN (Root & Fixpunkte)
        // ==============================================================================

        // Root
        rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0,0,0));
        tissues.push_back(rootSystem);

        // --- Linker Anker (Fix) ---
        auto bodyLeft = std::make_shared<SSBody>("Body_Left", MWMath::Point3D(-anchorDist, 0, 0), MWMath::RotMatrix3x3(), rootSystem);
        rootSystem->Children.push_back(bodyLeft);
        tissues.push_back(bodyLeft);

        // --- Rechter Anker (Fix) ---
        auto bodyRight = std::make_shared<SSBody>("Body_Right", MWMath::Point3D(anchorDist, 0, 0), MWMath::RotMatrix3x3(), rootSystem);
        rootSystem->Children.push_back(bodyRight);
        tissues.push_back(bodyRight);

        // --- Pusher (Beweglich in der Mitte) ---
        // Startet bei Y = pushStart
        auto bodyPusher = std::make_shared<SSBody>("Body_Pusher", MWMath::Point3D(0, pushStart, 0), MWMath::RotMatrix3x3(), rootSystem);
        rootSystem->Children.push_back(bodyPusher);
        bodyPusher->PositionGlobal = STARTP;
        bodyPusher->OrientationGlobal = FINGERSTARTORIENTATION;
        tissues.push_back(bodyPusher);


        // ==============================================================================
        // 2. MESHES HINZUFÜGEN
        // ==============================================================================
        
        // Mesh Links
        auto meshLeft = std::make_shared<SSEllipsoidMesh>(anchorRadius, anchorRadius, anchorRadius);
        meshLeft->Name = "Mesh_Left"; meshLeft->MeshColor = {0,0,1}; // Blau
        bodyLeft->Meshes.push_back(meshLeft); meshLeft->Parent = bodyLeft; meshes.push_back(meshLeft);

        // Mesh Rechts
        auto meshRight = std::make_shared<SSEllipsoidMesh>(anchorRadius, anchorRadius, anchorRadius);
        meshRight->Name = "Mesh_Right"; meshRight->MeshColor = {0,1,0}; // Grün
        bodyRight->Meshes.push_back(meshRight); meshRight->Parent = bodyRight; meshes.push_back(meshRight);

        // Mesh Pusher (Rot)
        auto meshPusher = std::make_shared<SSEllipsoidMesh>(PUSHA, PUSHB, PUSHC); // Kugel
        // Oder flacheres Ellipsoid:
        // auto meshPusher = std::make_shared<SSEllipsoidMesh>(pusherRadius*1.5, pusherRadius, pusherRadius*0.5); 
        meshPusher->Name = "Mesh_Pusher"; meshPusher->MeshColor = {1,0,0}; 
        bodyPusher->Meshes.push_back(meshPusher); meshPusher->Parent = bodyPusher; meshes.push_back(meshPusher);


        // ==============================================================================
        // INITIALES UPDATE
        // ==============================================================================
        rootSystem->update(0); 
        updateSceneMovement(currentScene, meshes, 0);

        // ==============================================================================
        // 3. MUSKEL
        // ==============================================================================
        
        int numPoints = (!cfg.muscleNumPoints.empty()) ? cfg.muscleNumPoints[0] : 20;
        
        // Muskel spannt sich von Links nach Rechts
        // Start/Endpunkte relativ zu den Anker-Bodies (z.B. oben drauf bei +Y)
        static SSMuscle mus("PushedMuscle", numPoints, 
            meshLeft.get(), {anchorRadius*1.1, 0.0, 0.0}, 
            meshRight.get(), {-anchorRadius*1.1, 0.0, 0.0});
        
        // Alle Meshes als Hindernisse (Pusher + Anker, falls er sich um die Anker wickeln soll)
        for(auto& m : meshes) {
            mus.meshPtrs.push_back(m.get());
        }

        mus.createMusclePoints();
        muscles.push_back(&mus);
    }
    else if (currentScene == "TORUS_PUSH_THROUGH" || currentScene == "TORUS_ROTATE_THROUGH") {
        
        // --- PARAMETER ---
        double anchorRadius = 0.2;
        double anchorDist = 2.0; 
        double pushStart = 0.0;  
        
        // Torus Dimensionen
        double torusMajorR = 0.4; // Radius des Rings
        double torusMinorR = 0.25; // Dicke des Rohrs
        
        // ==============================================================================
        // 1. STRUKTUR (Root & Fixpunkte) - Identisch zum Ellipsoid-Case
        // ==============================================================================

        rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0,0,0));
        tissues.push_back(rootSystem);

        // Linker Anker
        auto bodyLeft = std::make_shared<SSBody>("Body_Left", MWMath::Point3D(-anchorDist, 0, 0), MWMath::RotMatrix3x3(), rootSystem);
        rootSystem->Children.push_back(bodyLeft);
        tissues.push_back(bodyLeft);

        // Rechter Anker
        auto bodyRight = std::make_shared<SSBody>("Body_Right", MWMath::Point3D(anchorDist, 0, 0), MWMath::RotMatrix3x3(), rootSystem);
        rootSystem->Children.push_back(bodyRight);
        tissues.push_back(bodyRight);

        // Pusher (Mitte)
        auto bodyPusher = std::make_shared<SSBody>("Body_Pusher", MWMath::Point3D(0, pushStart, 0), MWMath::axisAngle({0, 1, 0}, 45.0), rootSystem);
        rootSystem->Children.push_back(bodyPusher);
        tissues.push_back(bodyPusher);


        // ==============================================================================
        // 2. MESHES
        // ==============================================================================
        
        // Anker-Meshes (Blau)
        auto meshLeft = std::make_shared<SSEllipsoidMesh>(anchorRadius, anchorRadius, anchorRadius);
        meshLeft->Name = "Mesh_Left"; meshLeft->MeshColor = {0,0,1};
        bodyLeft->Meshes.push_back(meshLeft); meshLeft->Parent = bodyLeft; meshes.push_back(meshLeft);

        auto meshRight = std::make_shared<SSEllipsoidMesh>(anchorRadius, anchorRadius, anchorRadius);
        meshRight->Name = "Mesh_Right"; meshRight->MeshColor = {0,0,1};
        bodyRight->Meshes.push_back(meshRight); meshRight->Parent = bodyRight; meshes.push_back(meshRight);

        // --- TORUS PUSHER (Rot) ---
        auto meshTorus = std::make_shared<SSTorusMesh>(torusMajorR, torusMinorR);
        meshTorus->Name = "Mesh_Torus_Pusher"; 
        //meshTorus->B = 1.5; meshTorus->C = 2.0;
        meshTorus->MeshColor = {1, 0, 0}; // Rot
        
        // OPTIONAL: Rotation des Torus
        // 0 Grad: Liegt flach wie ein Donut auf dem Tisch (Muskel trifft auf die Seite)
        // 90 Grad um X: Steht wie ein Autoreifen (Muskel läuft über die Lauffläche) [Image of Torus Orientation 90deg]
        // Wir drehen ihn um 90 Grad um die X-Achse:
        meshTorus->Orientation2ParentRel = MWMath::axisAngle({0, 1, 0}, 0.0);
        // Oder leicht schräg für Chaos:
        // meshTorus->Orientation2ParentRel = MWMath::axisAngle({1, 0, 0}, 45.0);

        bodyPusher->Meshes.push_back(meshTorus); 
        meshTorus->Parent = bodyPusher; 
        meshes.push_back(meshTorus);


        // ==============================================================================
        // INITIALES UPDATE
        // ==============================================================================
        rootSystem->update(0); 
        // updateSceneMovement(currentScene, meshes, 0);

        // ==============================================================================
        // 3. MUSKEL
        // ==============================================================================
        
        int numPoints = (!cfg.muscleNumPoints.empty()) ? cfg.muscleNumPoints[0] : 25;
        
        static SSMuscle mus("TorusMuscle", numPoints, 
            meshLeft.get(), {0.0, anchorRadius*1.1, 0.0}, 
            meshRight.get(), {0.0, anchorRadius*1.1, 0.0});
        
        for(auto& m : meshes) mus.meshPtrs.push_back(m.get());

        mus.createMusclePoints();
        mus.updateMusclePointsParents();
        muscles.push_back(&mus);

        
    }
    else if (currentScene == "CYLINDER_PUSH_THROUGH") {
        
        // --- PARAMETER ---
        double anchorRadius = 0.2;
        double pusherRadius = 0.3; // Radius des Zylinders
        double pusherHeight = 0.8; // Länge des Zylinders
        double anchorDist = 2.0; 
        double pushStart = 0.5;  
        
        // ==============================================================================
        // 1. STRUKTUR (Identisch zu Ellipsoid)
        // ==============================================================================

        rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0,0,0));
        tissues.push_back(rootSystem);

        // Linker Anker
        auto bodyLeft = std::make_shared<SSBody>("Body_Left", MWMath::Point3D(-anchorDist, 0, 0), MWMath::RotMatrix3x3(), rootSystem);
        rootSystem->Children.push_back(bodyLeft);
        tissues.push_back(bodyLeft);

        // Rechter Anker
        auto bodyRight = std::make_shared<SSBody>("Body_Right", MWMath::Point3D(anchorDist, 0, 0), MWMath::RotMatrix3x3(), rootSystem);
        rootSystem->Children.push_back(bodyRight);
        tissues.push_back(bodyRight);

        // Pusher Body
        auto bodyPusher = std::make_shared<SSBody>("Body_Pusher", MWMath::Point3D(0, pushStart, 0), MWMath::RotMatrix3x3(), rootSystem);
        rootSystem->Children.push_back(bodyPusher);
        bodyPusher->PositionGlobal = STARTP;
        bodyPusher->OrientationGlobal = FINGERSTARTORIENTATION;
        tissues.push_back(bodyPusher);


        // ==============================================================================
        // 2. MESHES (Hier Zylinder statt Ellipsoid)
        // ==============================================================================
        
        // Links (Kugel/Ellipsoid)
        auto meshLeft = std::make_shared<SSEllipsoidMesh>(anchorRadius, anchorRadius, anchorRadius);
        meshLeft->Name = "Mesh_Left"; meshLeft->MeshColor = {0,0,1};
        bodyLeft->Meshes.push_back(meshLeft); meshLeft->Parent = bodyLeft; meshes.push_back(meshLeft);

        // Rechts (Kugel/Ellipsoid)
        auto meshRight = std::make_shared<SSEllipsoidMesh>(anchorRadius, anchorRadius, anchorRadius);
        meshRight->Name = "Mesh_Right"; meshRight->MeshColor = {0,1,0};
        bodyRight->Meshes.push_back(meshRight); meshRight->Parent = bodyRight; meshes.push_back(meshRight);

        // --- PUSHER (ZYLINDER) ---
        // Konstruktor Annahme: SSCylinderMesh(radius, height)
        auto meshPusher = std::make_shared<SSCylinderMesh>(pusherRadius, pusherHeight);
        
        meshPusher->Name = "Mesh_Pusher_Cyl"; 
        meshPusher->MeshColor = {1, 0, 0}; // Rot
        
        // OPTIONAL: Orientierung anpassen
        // Zylinder werden oft entlang der Y-Achse (Höhe) erstellt. 
        // Wenn der Muskel entlang der X-Achse läuft, wollen wir vielleicht, 
        // dass der Zylinder quer dazu liegt (z.B. entlang Z-Achse wie eine Walze).
        // Standardmäßig zeigt er oft nach oben (Y). Das ist okay, dann trifft der Muskel die runde Seite.
        // Falls du ihn drehen willst (z.B. liegend):
        // meshPusher->Orientation2ParentRel = MWMath::axisAngle({1, 0, 0}, 90.0);

        bodyPusher->Meshes.push_back(meshPusher); 
        meshPusher->Parent = bodyPusher; 
        meshes.push_back(meshPusher);


        // ==============================================================================
        // INITIALES UPDATE
        // ==============================================================================
        rootSystem->update(0); 
        updateSceneMovement(currentScene, meshes, 0);

        // ==============================================================================
        // 3. MUSKEL
        // ==============================================================================
        
        int numPoints = (!cfg.muscleNumPoints.empty()) ? cfg.muscleNumPoints[0] : 20;
        
        static SSMuscle mus("CylinderMuscle", numPoints, 
            meshLeft.get(), {anchorRadius*1.1, 0.0, 0.0}, 
            meshRight.get(), {-anchorRadius*1.1, 0.0, 0.0});
        
        for(auto& m : meshes) {
            mus.meshPtrs.push_back(m.get());
        }

        mus.createMusclePoints();
        muscles.push_back(&mus);
    }
    else if (currentScene == "TWO_SEGMENTS_SIMPLE") {
        
        double L1 = 0.6; 
        double L2 = 0.4; 
        double width = 0.15; 

        rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0,0,0));
        tissues.push_back(rootSystem);

        auto createSegment = [&](std::string name, double w, double length) -> std::shared_ptr<SSMesh> {
            auto body = std::make_shared<SSBody>("Body_" + name, MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), rootSystem);
            tissues.push_back(body);
            rootSystem->Children.push_back(body);

            auto mesh = std::make_shared<SSEllipsoidMesh>(w, w, length);
            mesh->Name = "Mesh_" + name;
            mesh->MeshColor = {0.8, 0.8, 0.8}; 
            mesh->C = length; 

            // Rotation: Z -> X
            mesh->Orientation2ParentRel = MWMath::axisAngle({0,1,0}, 90.0);
            
            // WICHTIG: Offset ist JETZT 0, weil der Body im Zentrum (CoM) ist!
            // In der vorherigen Version war es length/2.0. Jetzt ist es 0.
            mesh->Position2ParentRelInParentFrame = {0, 0, 0}; 
            

            body->Meshes.push_back(mesh);
            mesh->Parent = body;
            meshes.push_back(mesh);
            return mesh;
        };

        auto mSeg1 = createSegment("Seg1", width, L1); 
        auto mSeg2 = createSegment("Seg2", width, L2); 
        mSeg1->MeshColor = {1,1,0}; // Rot

        // --- TORUS ---
        // Der Torus soll am ENDE von Seg1 sitzen.
        // Seg1 Body ist in der Mitte. Das Ende ist bei +L1/2.
        auto body1 = mSeg1->Parent;
        auto mTorus = std::make_shared<SSTorusMesh>(width * 2.9, 0.2);
        mTorus->Name = "Torus_Joint";
        mTorus->Orientation2ParentRel = MWMath::axisAngle({0,1,0}, 90.0); // Um X-Achse
        mTorus->MeshColor = {1,0,0};

        // Offset vom Body1 (Mitte) zum Ende:
        mTorus->Position2ParentRelInParentFrame = {L1/3, 0, 0}; 

        body1->Meshes.push_back(mTorus);
        mTorus->Parent = body1; 
        meshes.push_back(mTorus); 

        updateSceneMovement(currentScene, meshes, 0.0);

        // --- MUSKEL ---
        int numPoints = 20;
        SSMuscle* flexor = new SSMuscle("Flexor", numPoints, 
            mSeg1.get(), {-L1/2.0, -width * 1.5, 0.0}, // Start am Anfang von Seg1 (Relativ zur Mitte = -L/2)
            mSeg2.get(), {L2/2.0, -width * 1.5, 0.0}); // Ende am Ende von Seg2 (Relativ zur Mitte = +L/2)

        for(auto& m : meshes) flexor->meshPtrs.push_back(m.get());
        flexor->createMusclePoints();
        flexor->updateMusclePointsParents();
        muscles.push_back(flexor);
        
    }
    else if(currentScene == "FINGER_SIMPLE2"){
        // --- PARAMETER ---
        double handLength = 1.8;
        std::vector<double> fLength = {0.3482 * handLength, 0.2027 * handLength, 0.1175 * handLength, 0.0882 * handLength};
        double L1 = fLength[0]; // Länge Segment 1
        double L2 = fLength[1]; // Länge Segment 2
        double width = 0.2; // Dicke

        // 1. ROOT
        rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), nullptr);
        tissues.push_back(rootSystem);

        // ==============================================================================
        // SEGMENT 1 (Proximal) - Fixiert an Root
        // ==============================================================================
        
        // Body 1 (Startpunkt/Gelenk 1)
        auto body1 = std::make_shared<SSBody>("Body_Seg1", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), rootSystem);
        tissues.push_back(body1);

        // Mesh 1
        auto mesh1 = std::make_shared<SSEllipsoidMesh>(width, width, L1);
        mesh1->Name = "Mesh_Seg1";
        mesh1->MeshColor = {0.8, 0.8, 0.0}; // Gelb

        // Ausrichtung: Z -> X (damit Länge entlang X ist)
        mesh1->Orientation2ParentRel = MWMath::axisAngle({0,1,0}, 90.0);
        // Position: Body ist am Anfang (0), Mesh-Mitte ist bei L/2
        mesh1->Position2ParentRelInParentFrame = {0, 0, 0};

        body1->Meshes.push_back(mesh1);
        mesh1->Parent = body1;
        meshes.push_back(mesh1);


        // ==============================================================================
        // GELENK (Joint) - Am Ende von Segment 1
        // ==============================================================================
        
        // Position relativ zu Body 1: Am Ende des Knochens (L1 entlang X)
        MWMath::Point3D jointPosRel = {L1, 0, 0};

        auto joint = std::make_shared<SSJoint>("Joint_1", 
                                             jointPosRel,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body1,                 // Parent Body
                                             90.0,                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        
        tissues.push_back(joint);


        // ==============================================================================
        // SEGMENT 2 (Distal) - Hängt am Joint
        // ==============================================================================

        // Body 2 (Startpunkt von Segment 2, sitzt direkt am Joint)
        // Relativ zum Joint ist die Position (0,0,0)
        auto body2 = std::make_shared<SSBody>("Body_Seg2", MWMath::Point3D(L2,0,0), MWMath::RotMatrix3x3(), joint);
        tissues.push_back(body2);

        // Mesh 2
        auto mesh2 = std::make_shared<SSEllipsoidMesh>(width, width, L2);
        mesh2->Name = "Mesh_Seg2";
        mesh2->MeshColor = {0.0, 0.8, 0.0}; // Grün
        mesh2->C = L2;

        // Ausrichtung: Z -> X
        mesh2->Orientation2ParentRel = MWMath::axisAngle({0,1,0}, 90.0);
        // Position: Body (Joint) ist am Anfang, Mesh-Mitte bei L/2
        mesh2->Position2ParentRelInParentFrame = {0, 0, 0};

        body2->Meshes.push_back(mesh2);
        mesh2->Parent = body2;
        meshes.push_back(mesh2);


        // ==============================================================================
        // TORUS (Am Ende von Segment 1 / auf dem Gelenk)
        // ==============================================================================
        
        // Wir hängen den Torus an Body 1 (er bewegt sich nicht mit dem Gelenk mit, sondern ist fest am 1. Knochen)
        auto mTorus = std::make_shared<SSTorusMesh>(mesh1->B * 1.5, mesh1->B * 0.5);
        mTorus->Name = "Torus_Joint";
        mTorus->MeshColor = {1, 0, 0}; // Rot

        // Ausrichtung: Loch entlang X (wie Knochen)
        mTorus->Orientation2ParentRel = MWMath::axisAngle({0,1,0}, 90.0);
        
        // Position: Am Ende von L1 (auf dem Gelenk)
        mTorus->Position2ParentRelInParentFrame = {0, 0, L1/2};

        body1->Meshes.push_back(mTorus);
        mTorus->Parent = body1;
        meshes.push_back(mTorus);


        // ==============================================================================
        // INITIALES UPDATE
        // ==============================================================================
        // Damit alle Matrizen einmal berechnet werden (wichtig für Muskel-Punkte)
        rootSystem->update(0); 


        // ==============================================================================
        // MUSKEL
        // ==============================================================================
        int numPoints = 25;
        MWMath::Point3D startOffset = {0.0, -mesh1->B * 1.1, 0.0};
        MWMath::Point3D endOffset = {0.0, -mesh2->B * 1.1, 0.0};
        SSMuscle* flexor = new SSMuscle("Flexor", numPoints, 
            mesh1.get(), startOffset, 
            mesh2.get(), endOffset);

        // Alle Hindernisse
        for(auto& m : meshes) {
             flexor->meshPtrs.push_back(m.get());
        }

        flexor->createMusclePoints();
        flexor->updateMusclePointsParents();
        muscles.push_back(flexor);
    }
    else if (currentScene == "FINGER_SIMPLE4"){
        // --- PARAMETER ---
        double handLength = 1.8;
        std::vector<double> fLength = {0.3482 * handLength, 0.2027 * handLength, 0.1175 * handLength, 0.0882 * handLength};
        double L1 = fLength[0]; // Länge Segment 1
        double L2 = fLength[1]; // Länge Segment 2
        double L3 = fLength[2]; // Länge Segment 3
        double L4 = fLength[3]; // Länge Segment 4
        std::vector<double> width = {0.15, 0.1, 0.07, 0.05}; // Dicke
        std::vector<double> relTorusPos = {0.6, 0.42, 0.35};
        std::vector<double> relTorusR = {1.85, 1.9, 2.3};
        std::vector<double> relTorusr = {1.0, 1.0, 1.2};

        // 1. ROOT
        rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), nullptr);
        tissues.push_back(rootSystem);

        // ==============================================================================
        // SEGMENT 1 (Proximal) - Fixiert an Root
        // ==============================================================================
        auto body1 = std::make_shared<SSBody>("Body_Seg1", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), rootSystem);
        tissues.push_back(body1);
        // Mesh 1
        auto mesh1 = std::make_shared<SSEllipsoidMesh>(width[0], width[0], L1, "Mesh_Seg1", body1, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(mesh1);
        // Torus 1
        auto mTorus1 = std::make_shared<SSTorusMesh>(mesh1->B * relTorusR[0], mesh1->B * relTorusr[0], 
             "Torus1", body1, MWMath::Point3D(0, 0, L1 * relTorusPos[0]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mTorus1);

        MWMath::Point3D jointPosRel1 = MWMath::Point3D(L1, 0, 0);
        auto joint = std::make_shared<SSJoint>("Joint_1", 
                                             jointPosRel1,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body1,                 // Parent Body
                                             45.0,                  // Max Winkel
                                             MWMath::Point3D(0,1,0),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint);

        auto jMesh1 = std::make_shared<SSEllipsoidMesh>(width[0]*1.1, width[0]*1.1, width[0]*1.1, "Mesh_Joint1", joint, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh1);
        


        // ==============================================================================
        // SEGMENT 2 (Distal) - Hängt am Joint
        // ==============================================================================
        auto body2 = std::make_shared<SSBody>("Body_Seg2", MWMath::Point3D(L2,0,0), MWMath::RotMatrix3x3(), joint);
        tissues.push_back(body2);
        auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2, 
             "Mesh_Seg2", body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.0, 0.8, 0.0));
        meshes.push_back(mesh2);
        // Torus 2
        auto mTorus2 = std::make_shared<SSTorusMesh>(mesh2->B * relTorusR[1], mesh2->B * relTorusr[1], 
             "Torus2", body2, MWMath::Point3D(0, 0, L2 * relTorusPos[1]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
        meshes.push_back(mTorus2);


        MWMath::Point3D jointPosRel2 = MWMath::Point3D(L2, 0, 0);
        auto joint2 = std::make_shared<SSJoint>("Joint_2", 
                                             jointPosRel2,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body2,                 // Parent Body
                                             100.0,                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint2);

        auto jMesh2 = std::make_shared<SSEllipsoidMesh>(width[1]*1.1, width[1]*1.1, width[1]*1.1, "Mesh_Joint2", joint2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh2);

        // SEGMENT 3
        auto body3 = std::make_shared<SSBody>("Body_Seg3", MWMath::Point3D(L3,0,0), MWMath::RotMatrix3x3(), joint2);
        tissues.push_back(body3);
        // Mesh 3
        auto mesh3 = std::make_shared<SSEllipsoidMesh>(width[2], width[2], L3, 
             "Mesh_Seg3", body3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), 
             MWMath::Point3D(0.0, 0.8, 0.8));
        meshes.push_back(mesh3);
        // Torus 3
        auto mTorus3 = std::make_shared<SSTorusMesh>(mesh3->B * relTorusR[2], mesh3->B * relTorusr[2], 
             "Torus3", body3, MWMath::Point3D(0, 0, L3 * relTorusPos[2]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
        meshes.push_back(mTorus3);


        MWMath::Point3D jointPosRel3 = MWMath::Point3D(L3, 0, 0);
        auto joint3 = std::make_shared<SSJoint>("Joint_3", 
                                             jointPosRel3,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body3,                 // Parent Body
                                             90.0,                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint3);

        auto jMesh3 = std::make_shared<SSEllipsoidMesh>(width[2]*1.1, width[2]*1.1, width[2]*1.1, "Mesh_Joint3", joint3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh3);

        // SEGMENT 4
        auto body4 = std::make_shared<SSBody>("Body_Seg4", MWMath::Point3D(L4,0,0), MWMath::RotMatrix3x3(), joint3);
        tissues.push_back(body4);
        // Mesh 4
        auto mesh4 = std::make_shared<SSEllipsoidMesh>(width[3], width[3], L4, 
             "Mesh_Seg4", body4, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), 
             MWMath::Point3D(0.9, 0.5, 0.3));
        meshes.push_back(mesh4);

        // ==============================================================================
        // INITIALES UPDATE
        // ==============================================================================
        // Damit alle Matrizen einmal berechnet werden (wichtig für Muskel-Punkte)
        for (auto& m : meshes) {m->InitializeMesh();}
        rootSystem->update(0); 


        // ==============================================================================
        // MUSKEL
        // ==============================================================================
        int numPoints = 25;
        MWMath::Point3D startOffset = MWMath::Point3D(0.0, -mesh1->B * 1.1, 0.0);
        MWMath::Point3D endOffset = MWMath::Point3D(0.0, -mesh4->B * 1.1, 0.0);
        SSMuscle* flexor = new SSMuscle("Flexor", numPoints, 
            mesh1.get(), startOffset, 
            mesh4.get(), endOffset);

        // Alle Hindernisse
        for(auto& m : meshes) {
             flexor->meshPtrs.push_back(m.get());
        }

        flexor->createMusclePoints();
        flexor->updateMusclePointsParents();
        muscles.push_back(flexor);
    }

}


int main(int argc, char** argv)
{

    // versuch simpel
    QApplication app(argc, argv); // VTK/Qt brauchen das QApplication Objekt
    std::string configPath = "../config.txt"; // Liegt im Build-Ordner
    SimSettings cfg = ConfigManager::loadConfig(configPath);
    ConfigManager::saveCopy(configPath);
    MAXJOINTANGLES = cfg.MAXJOINTANGLES; 
    std::string currentScene = cfg.currentScene;
    int numTimeSteps = cfg.numTimeSteps;
    int objFunc = cfg.objFunc;

    // DYNAMISCHES HINDERNIS-SZENARIO
    
    std::vector<std::shared_ptr<SSMesh>> meshes;
    std::vector<SSMuscle*> musclePtrs;
    std::shared_ptr<SSBody> rootSystem;
    std::vector<std::shared_ptr<SSTissue>> tissue;
    // hier szene objektorientert aufbauen per funktion
    //...setupScene(currentScene, meshes, musclePtrs);
    setupSceneObjectOriented(currentScene, tissue, meshes, musclePtrs, rootSystem, numTimeSteps, cfg);
    // re-dicretize meshes for new position
    for (auto& m : meshes) {
        m->discretizeMesh(cfg.discretization); // discretize for distance calculation 
    }

    // --- 3. SOLVER SETUP ---
    // CasadiSystem system(musclePtrs, objFunc);
    std::vector<CasadiSystem*> systems;
    for (auto* mus : musclePtrs) {
        std::vector<SSMuscle*> singleMuscleList = {mus};
        systems.push_back(new CasadiSystem(singleMuscleList, objFunc, cfg.solverMethod, cfg.casadiParametrization, cfg.bUseManualJacobian, cfg.bSumPhiEta));
        //systems.back()->CasadiSystemName = "CasSys_" + mus->Name;
        // qDebug() << "MUSCLE POINTS" << mus->MNodes.size();
    }

    // Container für n Muskeln
    int numMuscles = musclePtrs.size();
    std::vector<std::vector<std::vector<MWMath::Point3D>>> allMuscleResults(numMuscles);
    std::vector<std::vector<std::vector<MWMath::Point3D>>> allInitialGuesses(numMuscles);
    std::vector<std::vector<std::vector<MWMath::Point3D>>> allOptimizedPointColors(numMuscles);
    std::vector<std::vector<MWMath::RotMatrix3x3>> allMeshResults(meshes.size()); 
    std::vector<std::vector<MWMath::Point3D>> otherPointsWithColors(numTimeSteps); // {p1, c1, p2, c2, ... }
    std::vector<double> angles;

    // MUSCLE INITIALIZATION
    auto start = std::chrono::high_resolution_clock::now();
    qDebug() << "======================= MUSCLE INITIALIZATION ========================";
    for (auto* mus : musclePtrs) {
        mus->initializeSimulationMuscle(numTimeSteps);
        
        mus->getMuscleInfo();
        mus->checkCollision();
        qDebug() << "-----------------------------------------------------------";
    }

    // --- 4. SIMULATION LOOP ---
    for(int t = 0; t < numTimeSteps; ++t) {
        qDebug() << "====================== TIMESTEP " << t << " / " << numTimeSteps-1 << "============================================";
        double progress = (double)t / (double)numTimeSteps;
        angles.push_back((double)t);

        
        // MOVE MESHES / UPDATE SCENE
        if (currentScene == "NONE" || currentScene == "FINGER_SIMPLE2" || currentScene == "FINGER_SIMPLE4") {  
            if (rootSystem) {
                // Ein einziger Aufruf aktualisiert den gesamten Baum!
                rootSystem->update(t); 
                for (auto& m : meshes) {
                        m->MeshPointsGlobal.push_back(m->PositionGlobal);
                        m->allRMatrixGlobal.push_back(m->OrientationGlobal);
                }
                    
                // re-dicretize meshes for new position
                for (auto& m : meshes) {
                    m->discretizeMesh(cfg.discretization); // discretize for distance calculation 
                }

                /* for(auto tis : tissue){
                    tis->getInfo();
                } */
            }
        }
        else{
            updateSceneMovement(currentScene, meshes, progress);
            for (auto& m : meshes) {
                        m->MeshPointsGlobal.push_back(m->PositionGlobal);
                        m->allRMatrixGlobal.push_back(m->OrientationGlobal);
                }
                    
            // re-dicretize meshes for new position
            for (auto& m : meshes) {
                m->discretizeMesh(cfg.discretization); // discretize for distance calculation 
            }
        }

        // --- SCHRITT B & C: ENDPUNKTE & INITIAL GUESS PRO MUSKEL ---
        for(int m = 0; m < numMuscles; ++m) {
            auto* mus = musclePtrs[m];
            // Fixpunkte prädizieren
            mus->OriginPointGlobal = mus->MNodes[0].predictNewGlobal();
            mus->InsertionPointGlobal = mus->MNodes.back().predictNewGlobal();
            mus->MNodes[0].PositionGlobal = mus->OriginPointGlobal;
            mus->MNodes.back().PositionGlobal = mus->InsertionPointGlobal;
            // Guess Pfad sammeln
            std::vector<MWMath::Point3D> guessPath;
            for(auto& node : mus->MNodes) {
                guessPath.push_back(node.predictNewGlobal());
            }
            mus->storeMNodesInitialGuess(t, guessPath);
        }

        // --- SCHRITT D: MULTI-OPTIMIERUNG ---
        /* for (auto muscle : musclePtrs) {
            for (auto& node : muscle->MNodes) {
                qDebug() << "LocalPosition: " << node.PositionLocal.x << node.PositionLocal.y << node.PositionLocal.z;
            }
        } */

        for (auto* sys : systems) {
            qDebug() << "---------------------------" << " Solve CasadiSystem: " << QString::fromStdString(sys->CasadiSystemName) << " ---------------------------";
            sys->solveStepX(); // solves all systems
        }
        qDebug() << "-------------------------------------------------------------------------------------------------";
        
        /* for (auto muscle : musclePtrs) {
            for (auto& node : muscle->MNodes) {
                qDebug() << "LocalPosition: " << node.PositionLocal.x << node.PositionLocal.y << node.PositionLocal.z;
            }
        } */
        if (cfg.dynamicReparametrization)
        {    
            for (auto muscle : musclePtrs) {
                muscle->updateMusclePointsParentsLocal();// updateMusclePointsParents();
            }
        }

        for (auto muscle : musclePtrs) {
            muscle->getViaPointNodeInfo();
            muscle->checkCollision();
        }

        
        // --- SCHRITT E, F & G: PARENTS UPDATEN & ERGEBNISSE SPEICHERN ---
        for(int m = 0; m < numMuscles; ++m) {
            auto* mus = musclePtrs[m];

            // Farben sammeln
            std::vector<MWMath::Point3D> cols;
            for(auto& node : mus->MNodes) cols.push_back(node.getparentMeshColor());
            //GLOBALSAVE allOptimizedPointColors[m].push_back(cols);
            mus->storeMNodesInitialGuessColors(t, cols);

            // Pfad speichern
            mus->storeMNodesGlobalPositions();
            mus->computeMuscleLength(true);

            /* int im = 0;
            for (auto& node : mus->MNodes) {
                node.getMNodeInfoStep(im);
                im++;
            } 
            qDebug() << "Cfg: " << cfg.muscleNumPoints[0];
            qDebug() << "MuscleNodes" << im; */
        }

        // plot further points
        if (cfg.bShowDiscretization){    
            for (std::shared_ptr<SSMesh> mesh : meshes) {
                for (MWMath::Point3D& p : mesh->GlobalDiscreteMeshPoints) {
                    otherPointsWithColors[t].push_back(p);
                    otherPointsWithColors[t].push_back(mesh->MeshColor);
                }
            }
        }
    }
    // ------------------- END OF SIMULATION LOOP -------------------



    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // test zum neubefüllen
    allMuscleResults.clear();
    allInitialGuesses.clear();
    allOptimizedPointColors.clear();
    allMuscleResults.resize(numMuscles);
    allInitialGuesses.resize(numMuscles);
    allOptimizedPointColors.resize(numMuscles);

    for (int m = 0; m < numMuscles; ++m) {
        auto* mus = musclePtrs[m];
        
        // get all values per muscle
        std::vector<std::vector<MWMath::Point3D>> globalPoints;
        std::vector<std::vector<MWMath::Point3D>> localPoints;
        std::vector<std::vector<MWMath::Point3D>> initialGuessPoints;
        std::vector<std::vector<MWMath::Point3D>> initialGuessColors;
        std::vector<std::vector<std::vector<double>>> etaValues;
        mus->getAllMuscleMNodesStepValues(numTimeSteps, globalPoints, localPoints, initialGuessPoints, initialGuessColors, etaValues);
        allMuscleResults[m] = globalPoints;
        allInitialGuesses[m] = initialGuessPoints;
        allOptimizedPointColors[m] = initialGuessColors;
    }

    std::cout << "Dauer: " << duration.count()/1000 << " s" << std::endl;
    
    for (int m = 0; m < systems.size(); ++m) {
        auto* sys = systems[m];
        qDebug() << "=== Solver Convergence - " << QString::fromStdString(sys->CasadiSystemName) << "===";
        int i = 0;
        for (const auto& msg : sys->SolverConvergenceMessages) {
            qDebug() << "Step " << i << "/" << sys->SolverConvergenceMessages.size()-1 << ": " << QString::fromStdString(msg) << " (Iter: " << sys->SolverConvergenceSteps[i] << ")";
            i++;
            for (auto* mus : sys->m_muscles) {
                double len = mus->MuscleLengthSteps[i];
                qDebug() << "    " << QString::fromStdString(mus->Name) << " Length: " << len << " m";

            }
        }
    }

    // --- 5. VIEWER ---
    std::vector<std::vector<std::string>> viewerText;
    for (int t = 0; t < numTimeSteps; ++t) {
        std::vector<std::string> stepText;
        stepText.push_back("Scene: " + currentScene);
        stepText.push_back(std::string("OwnConstraintJacobian: ") + (cfg.bUseManualJacobian ? "Yes" : "No"));
        stepText.push_back("Time (s): " + QString::number(duration.count()/1000.0).toStdString());
        stepText.push_back("casadiParametrization: " + cfg.casadiParametrization);
        stepText.push_back(std::string("dynamicReparametrization: ") + (cfg.dynamicReparametrization ? "Yes" : "No"));
        stepText.push_back("Sucess: " + systems[0]->SolverConvergenceMessages[t]);
        viewerText.push_back(stepText);
    }
    VTKSimViewerSimple viewer(allMuscleResults, allMeshResults, allInitialGuesses, allOptimizedPointColors, meshes, musclePtrs, angles, otherPointsWithColors, viewerText);
    viewer.show();

    std::string outputName = systems[0]->CasadiSystemName + "InputParameters" + currentScene + ".txt";
    exportParameterLog(systems[0]->allParameterInputsAllSteps, systems[0]->allParameterInputDescriptionsAllSteps, outputName);

    for (auto* sys : systems) {
        for (auto* mus : sys->m_muscles) {
            exportMuscleLog(sys->CasadiSystemName, mus);
        }
    }

    return app.exec();

}