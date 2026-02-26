#include <QApplication>
#include "utils/rawSimulation.h"
#include "utils/ConfigLoader.h"
#include "utils/utility.h"
#include "utils/ModelBuilder.h"
#include <chrono>
#include <cstdlib>
#include <unordered_map>

// test
#include "simpleSimulation/casadiSystem.h"
#include "simpleSimulation/VTKSimViewerSimple.h"
#include "simpleSimulation/SSMuscle.h"
#include "simpleSimulation/SSMeshes.h"
#include "simpleSimulation/SSBody.h"

std::vector<double> FJAngles = {0.0, 0.0, 90.0, 90.0}; // MCP1, MCP2, PIP, DIP 
float VIEWSCALER = 1.0f; // VIEWER in [cm]
float GEOMETRYSCALER = 1.0f; // GEOMETRIE in [cm]
std::string UNITS = "cm"; // Einheiten für Export (z.B. "m" für Meter, "cm" für Zentimeter, etc.)



std::vector<JoingAngleDef> myMovement = {{5,  60.0},{100, 72.0},{2,  90.0}};






std::vector<double> MAXJOINTANGLES = {0.0, 0.0, 0.0}; // MCP, PIP, DIP
MWMath::RotMatrix3x3 FINGERSTARTORIENTATION = MWMath::axisAngle(MWMath::Point3D(0,0,1), 90.0); // Startorientierung des Fingermodells
MWMath::Point3D ROTAXIS = MWMath::Point3D(0,1,0).normed(); // Rotationsachse für Fingerbeugung
double ROTENDANBGLE = 90.0; // Endwinkel der Beugung

double PUSHA = 0.25, PUSHB = 0.4, PUSHC = 0.7; // Halbachsen des Pushers (Kugel)
MWMath::Point3D MOVEDIR = MWMath::Point3D(0, -1, 0).normed(); // Verschiebung des Fingermodells beim Drücken
MWMath::Point3D STARTP = MWMath::Point3D(0.0, 0.3, -0.0);//MWMath::Point3D(0.0, 0.5, -0.5); // Startpunkt des Drückers
double MOVEDIST = 1.0; // Distanz, die der Drücker bewegt wird

MWMath::Point3D COLORF1 = MWMath::Point3D(0.50, 0.00, 1.00); 
MWMath::Point3D COLORF2 = MWMath::Point3D(0.30, 0.30, 1.00); 
MWMath::Point3D COLORF3 = MWMath::Point3D(0.00, 0.50, 1.00); 
MWMath::Point3D COLORF4 = MWMath::Point3D(0.00, 0.90, 1.00);

MWMath::Point3D COLORJOINT = MWMath::Point3D(0.0, 0.9, 0.3);
MWMath::Point3D COLORFDEACTIVE = MWMath::Point3D(0.5, 0.5, 0.5);

void updateSceneMovement(std::string sceneName, std::vector<std::shared_ptr<SSMesh>>& meshes, double progress){
    if (sceneName == "ELLIPSOID_PUSH_THROUGH" || sceneName == "TORUS_PUSH_THROUGH" || sceneName == "CYLINDER_PUSH_THROUGH" || sceneName == "ELLTORUS") {
        
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
        double endY = -360.0;
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

    Hand::MMHandAngles MMHA;

    std::vector<double> jointAnglesPrescribed;
    createJointMovementVector(jointAnglesPrescribed, myMovement);

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
    else if (currentScene == "OAUDE08_UPPER_LIMB") {
        
        // --- PARAMETER ---
        double s = 0.1f; // Skalierung auf Meter (IPOPT bevorzugt Werte nahe 1)
        
        double rSphere = 80.0 * s;
        MWMath::Point3D posSphere(-3.92 * s, -0.46 * s, -0.28 * s);
        
        double rCyl = 40.0 * s;
        double hCyl = 800.0 * s; // Groß genug, um Z = -300 abzudecken
        MWMath::Point3D posCyl(0.0, 0.0, 0.0);
        
        // Da der Zylinder auf (0,0,0) liegt, sind diese globalen Punkte 
        // gleichzeitig die exakten lokalen Offsets zum Zylinder-Zentrum.
        MWMath::Point3D pOrigin(-80.0 * s, 20.0 * s, 80.0 * s); 
        MWMath::Point3D pInsertion(40.0 * s, 0.0, -300.0 * s); 

        // ==============================================================================
        // 1. STRUKTUR AUFBAUEN
        // ==============================================================================

        rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), nullptr);
        tissues.push_back(rootSystem);

        // --- Sphere Body ---
        auto bodySphere = std::make_shared<SSBody>("Body_Sphere", posSphere, MWMath::RotMatrix3x3(), rootSystem);
        rootSystem->Children.push_back(bodySphere);
        tissues.push_back(bodySphere);

        // --- Cylinder Body ---
        auto bodyCyl = std::make_shared<SSBody>("Body_Cylinder", posCyl, MWMath::RotMatrix3x3(), rootSystem);
        rootSystem->Children.push_back(bodyCyl);
        tissues.push_back(bodyCyl);

        // ==============================================================================
        // 2. MESHES HINZUFÜGEN
        // ==============================================================================
        
        auto meshSphere = std::make_shared<SSEllipsoidMesh>(rSphere, rSphere, rSphere);
        meshSphere->Name = "Mesh_Sphere"; meshSphere->MeshColor = {1.0, 0.3, 0.3}; // Rot
        bodySphere->Meshes.push_back(meshSphere); meshSphere->Parent = bodySphere; meshes.push_back(meshSphere);

        auto meshCyl = std::make_shared<SSCylinderMesh>(rCyl, hCyl);
        meshCyl->Name = "Mesh_Cylinder"; meshCyl->MeshColor = {0.3, 0.6, 1.0}; // Blau
        bodyCyl->Meshes.push_back(meshCyl); meshCyl->Parent = bodyCyl; meshes.push_back(meshCyl);

        // ==============================================================================
        // MANUAL INITIAL UPDATE
        // ==============================================================================
        
        bodySphere->PositionGlobal = posSphere;
        bodySphere->OrientationGlobal = MWMath::RotMatrix3x3();
        bodyCyl->PositionGlobal = posCyl;
        bodyCyl->OrientationGlobal = MWMath::RotMatrix3x3();

        meshSphere->PositionGlobal = posSphere;
        meshSphere->OrientationGlobal = MWMath::RotMatrix3x3();
        meshCyl->PositionGlobal = posCyl;
        meshCyl->OrientationGlobal = MWMath::RotMatrix3x3();

        for (auto& m : meshes) {
            m->InitializeMesh(); // Bereitet Buffer vor
            m->MeshPointsGlobal.push_back(m->PositionGlobal);
            m->allRMatrixGlobal.push_back(m->OrientationGlobal);
            
            // FIX: Sphären/Ellipsoide müssen ihre wahre Größe (A,B,C) als Skalierung pushen!
            if (auto ell = std::dynamic_pointer_cast<SSEllipsoidMesh>(m)) {
                m->AllScalerStepLists.push_back(MWMath::Point3D(ell->A, ell->B, ell->C));
            } else {
                // Zylinder etc. skalieren sich im Viewer selbst über ihre Parameter
                m->AllScalerStepLists.push_back(MWMath::Point3D(1.0, 1.0, 1.0));
            }
        }
        
        rootSystem->update(0);

        // ==============================================================================
        // 3. MUSKEL
        // ==============================================================================
        
        int numPoints = (!cfg.muscleNumPoints.empty()) ? cfg.muscleNumPoints[0] : 45; 
        
        // Wir hängen Start UND Ende an das Zylinder-Mesh. 
        // pOrigin und pInsertion sind die lokalen Offsets zum Zylinder.
        SSMuscle* mus = new SSMuscle("Aude08_UpperLimb", numPoints, 
            meshCyl.get(), pOrigin,     
            meshCyl.get(), pInsertion); 
        
        mus->meshPtrs.push_back(meshSphere.get());
        mus->meshPtrs.push_back(meshCyl.get());

        mus->createMusclePoints();
        mus->updateMusclePointsParents();
        muscles.push_back(mus);
    }
    else if (currentScene == "OKINEMATIC_ELBOW_MODEL") {
        
        // --- PARAMETER ---
        double s = 0.1; // Skalierung auf Meter (1 mm = 0.001 m)
        
        // Dimensionen (Halbachsen a, b, c) laut Paper
        double uaA = 45.0 * s, uaB = 45.0 * s, uaC = 174.0 * s; // Upper Arm
        double elA = 43.0 * s, elB = 45.0 * s, elC = 38.0 * s;  // Elbow Segment
        double laA = 38.0 * s, laB = 38.0 * s, laC = 142.0 * s; // Lower Arm

        // ==============================================================================
        // 1. STRUKTUR AUFBAUEN (Hierarchie & Gelenke)
        // ==============================================================================

        rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), nullptr);
        tissues.push_back(rootSystem);

        // --- Upper Arm (Fix im Ursprung) ---
        auto bodyUpper = std::make_shared<SSBody>("Body_UpperArm", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), rootSystem);
        rootSystem->Children.push_back(bodyUpper);
        tissues.push_back(bodyUpper);

        // --- 2-DOF GELENK (im Zentrum des Ellbogens bei Z = -174) ---
        
        // Gelenk 1: Flexion/Extension (Y-Achse)
        MWMath::Point3D jointPosRel1(0.0, 0.0, -174.0 * s); 
        auto jointFlex = std::make_shared<SSJoint>("Joint_ElbowFlexion", 
                                                   jointPosRel1, 
                                                   MWMath::RotMatrix3x3(),
                                                   bodyUpper, 
                                                   FJAngles.size() > 0 ? FJAngles[0] : 90.0, 
                                                   MWMath::Point3D(0,1,0), // Y-Achse
                                                   numTimeSteps);
        tissues.push_back(jointFlex);

        // Gelenk 2: Axiale Rotation (Z-Achse). Liegt exakt auf Gelenk 1 (Offset 0)
        MWMath::Point3D jointPosRel2(0.0, 0.0, 0.0); 
        auto jointAxial = std::make_shared<SSJoint>("Joint_ElbowAxial", 
                                                    jointPosRel2, 
                                                    MWMath::RotMatrix3x3(),
                                                    jointFlex, 
                                                    90.f, // FJAngles.size() > 1 ? FJAngles[1] : 0.0, 
                                                    MWMath::Point3D(0,1,0), // Z-Achse
                                                    numTimeSteps);
        tissues.push_back(jointAxial);

        // --- Elbow Segment ---
        // Hängt am axialen Gelenk. Da global bei -174 und Gelenk bei -174 -> lokaler Offset 0
        auto bodyElbow = std::make_shared<SSBody>("Body_Elbow", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), jointAxial);
        tissues.push_back(bodyElbow);

        // --- Lower Arm ---
        // Global bei -316. Ellbogen war bei -174. Differenz = -142.
        auto bodyLower = std::make_shared<SSBody>("Body_LowerArm", MWMath::Point3D(0.0, 0.0, -142.0 * s), MWMath::RotMatrix3x3(), bodyElbow);
        tissues.push_back(bodyLower);

        // ==============================================================================
        // 2. MESHES HINZUFÜGEN
        // ==============================================================================
        
        // Mesh Upper Arm (Rot)
        auto meshUpper = std::make_shared<SSEllipsoidMesh>(uaA, uaB, uaC);
        meshUpper->Name = "Mesh_UpperArm"; meshUpper->MeshColor = {1.0, 0.3, 0.3}; 
        meshUpper->Position2ParentRelInParentFrame = {0,0,0};
        bodyUpper->Meshes.push_back(meshUpper); meshUpper->Parent = bodyUpper; meshes.push_back(meshUpper);

        // Mesh Elbow (Grün)
        auto meshElbow = std::make_shared<SSEllipsoidMesh>(elA, elB, elC);
        meshElbow->Name = "Mesh_Elbow"; meshElbow->MeshColor = {0.3, 1.0, 0.3}; 
        meshElbow->Position2ParentRelInParentFrame = {0,0,0};
        bodyElbow->Meshes.push_back(meshElbow); meshElbow->Parent = bodyElbow; meshes.push_back(meshElbow);

        // Mesh Lower Arm (Blau)
        auto meshLower = std::make_shared<SSEllipsoidMesh>(laA, laB, laC);
        meshLower->Name = "Mesh_LowerArm"; meshLower->MeshColor = {0.3, 0.6, 1.0}; 
        meshLower->Position2ParentRelInParentFrame = {0,0,0};
        bodyLower->Meshes.push_back(meshLower); meshLower->Parent = bodyLower; meshes.push_back(meshLower);

        // ==============================================================================
        // INITIALES UPDATE & VTK-SCALING FIX
        // ==============================================================================
        for (auto& m : meshes) {m->InitializeMesh();}
        rootSystem->update(0); 

        // ==============================================================================
        // 3. MUSKELN (Triceps & Biceps)
        // ==============================================================================
        
        int ptsTri = (!cfg.muscleNumPoints.empty()) ? cfg.muscleNumPoints[0] : 25; 
        int ptsBic = (cfg.muscleNumPoints.size() > 1) ? cfg.muscleNumPoints[1] : ptsTri; 
        
        // --- TRICEPS ---
        MWMath::Point3D pOriginTri(45.0 * s, 0.0, 122.0 * s);
        MWMath::Point3D pInsertionTri(50.0 * s, 0.0, 100.0 * s);
        
        SSMuscle* triceps = new SSMuscle("Triceps", ptsTri, 
            meshUpper.get(), pOriginTri,     
            meshLower.get(), pInsertionTri); 
        
        // Alle 3 Meshes als potentielle Hindernisse
        triceps->meshPtrs.push_back(meshUpper.get());
        triceps->meshPtrs.push_back(meshElbow.get());
        triceps->meshPtrs.push_back(meshLower.get());

        triceps->createMusclePoints();
        triceps->updateMusclePointsParents();
        muscles.push_back(triceps);

        // --- BICEPS ---
        MWMath::Point3D pOriginBic(-45.0 * s, 0.0, 139.0 * s);
        MWMath::Point3D pInsertionBic(0.0, 38.0 * s, 85.0 * s);
        
        SSMuscle* biceps = new SSMuscle("Biceps", ptsBic, 
            meshUpper.get(), pOriginBic,     
            meshLower.get(), pInsertionBic); 
        
        // Alle 3 Meshes als potentielle Hindernisse
        biceps->meshPtrs.push_back(meshUpper.get());
        biceps->meshPtrs.push_back(meshElbow.get());
        biceps->meshPtrs.push_back(meshLower.get());

        biceps->createMusclePoints();
        biceps->updateMusclePointsParents();
        muscles.push_back(biceps);
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
    else if (currentScene == "ELLTORUS"){
        // --- PARAMETER ---
        double anchorRadius = 0.2;
        double anchorDist = 2.0; 
        double pushStart = 0.0;  
        
        // Elliptical Torus Dimensionen (Halbachsen)
        double RX = 0.4;  // Große Halbachse Pfad (Ring-Ausdehnung X)
        double RY = 0.3; // Kleine Halbachse Pfad (Ring-Ausdehnung Y)
        double R1 = 0.2;  // Horizontale Halbachse Rohr-Querschnitt
        double R2 = 0.5; // Vertikale Halbachse Rohr-Querschnitt
        
        // ==============================================================================
        // 1. STRUKTUR (Root & Fixpunkte)
        // ==============================================================================

        rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0,0,0));
        tissues.push_back(rootSystem);

        // Linker Anker (Fixpunkt Muskelursprung)
        auto bodyLeft = std::make_shared<SSBody>("Body_Left", MWMath::Point3D(-anchorDist, 0, 0), MWMath::RotMatrix3x3(), rootSystem);
        rootSystem->Children.push_back(bodyLeft);
        tissues.push_back(bodyLeft);

        // Rechter Anker (Fixpunkt Muskelansatz)
        auto bodyRight = std::make_shared<SSBody>("Body_Right", MWMath::Point3D(anchorDist, 0, 0), MWMath::RotMatrix3x3(), rootSystem);
        rootSystem->Children.push_back(bodyRight);
        tissues.push_back(bodyRight);

        // Pusher/Rotator Body (Trägt den elliptischen Torus)
        auto bodyPusher = std::make_shared<SSBody>("Body_Pusher_Ellip", 
                                                MWMath::Point3D(0, pushStart, 0), 
                                                MWMath::axisAngle({0, 0, 1}, 0.0), 
                                                rootSystem);
        rootSystem->Children.push_back(bodyPusher);
        tissues.push_back(bodyPusher);


        // ==============================================================================
        // 2. MESHES
        // ==============================================================================
        
        // Anker-Meshes (Blau)
        auto meshLeft = std::make_shared<SSEllipsoidMesh>(anchorRadius, anchorRadius, anchorRadius);
        meshLeft->Name = "Mesh_Left"; meshLeft->MeshColor = {0, 0, 1};
        bodyLeft->Meshes.push_back(meshLeft); meshLeft->Parent = bodyLeft; meshes.push_back(meshLeft);

        auto meshRight = std::make_shared<SSEllipsoidMesh>(anchorRadius, anchorRadius, anchorRadius);
        meshRight->Name = "Mesh_Right"; meshRight->MeshColor = {0, 0, 1};
        bodyRight->Meshes.push_back(meshRight); meshRight->Parent = bodyRight; meshes.push_back(meshRight);

        // --- ELLIPTICAL TORUS PUSHER (Rot) ---
        // Nutzt die neue Klasse mit elliptischen Radien
        auto meshTorus = std::make_shared<SSEllipticalTorusMesh>(R1, R2, RX, RY);
        meshTorus->Name = "Mesh_Elliptical_Torus_Pusher"; 
        meshTorus->MeshColor = {1, 0, 0}; // Rot
        
        // Orientierung: 90 Grad um X, damit der Torus "steht" und der Muskel durch das Loch gleiten kann
        meshTorus->Orientation2ParentRel = MWMath::axisAngle({1, 0, 0}, 90.0);

        bodyPusher->Meshes.push_back(meshTorus); 
        meshTorus->Parent = bodyPusher; 
        meshes.push_back(meshTorus);


        // ==============================================================================
        // INITIALES UPDATE & DISKRETISIERUNG (Wichtig für Viewer & Solver)
        // ==============================================================================
        rootSystem->update(0); 
        
        for (auto& m : meshes) {
            // Generiert die diskreten Punkte basierend auf der elliptischen Parametrisierung
            m->discretizeMesh(cfg.discretization); 
        }

        // ==============================================================================
        // 3. MUSKEL
        // ==============================================================================
        
        int numPoints = (!cfg.muscleNumPoints.empty()) ? cfg.muscleNumPoints[0] : 25;
        
        // Wir verwenden einen Pointer auf SSMuscle (Heapspeicher), um Segfaults beim Szenenwechsel zu vermeiden
        SSMuscle* mus = new SSMuscle("EllipTorusMuscle", numPoints, 
                                    meshLeft.get(), {0.0, anchorRadius * 1.1, 0.0}, 
                                    meshRight.get(), {0.0, anchorRadius * 1.1, 0.0});
        
        // Dem Muskel die Mesh-Liste für die Kollisionsprüfung/Wicklung übergeben
        for(auto& m : meshes) mus->meshPtrs.push_back(m.get());

        mus->createMusclePoints();
        mus->updateMusclePointsParents();
        muscles.push_back(mus); 
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
    
    else if (currentScene == "OFINGER_SIMPLE_BigTorus"){
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
                                             90.0,                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
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
    else if (currentScene == "OFINGER_SIMPLE_BigTorusBelow"){
        // --- PARAMETER ---
        double handLength = 1.8;
        std::vector<double> fLength = {0.3482 * handLength, 0.2027 * handLength, 0.1175 * handLength, 0.0882 * handLength};
        double L1 = fLength[0]; // Länge Segment 1
        double L2 = fLength[1]; // Länge Segment 2
        double L3 = fLength[2]; // Länge Segment 3
        double L4 = fLength[3]; // Länge Segment 4
        std::vector<double> width = {0.15, 0.1, 0.07, 0.05}; // Dicke
        std::vector<double> relTorusPos = {0.62, 0.45, 0.35};
        std::vector<double> relTorusR = {1.0, 1.0, 1.2};
        std::vector<double> relTorusr = {0.9, 0.9, 1.1};

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
             "Torus1", body1, MWMath::Point3D(0, -width[0]*0.8, L1 * relTorusPos[0]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mTorus1);

        MWMath::Point3D jointPosRel1 = MWMath::Point3D(L1, 0, 0);
        auto joint = std::make_shared<SSJoint>("Joint_1", 
                                                jointPosRel1,           // Wo sitzt das Gelenk relativ zum Parent?
                                                MWMath::RotMatrix3x3(),// Initiale Rotation
                                                body1,                 // Parent Body
                                                FJAngles[0],                  // Max Winkel
                                                MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                                numTimeSteps);
        tissues.push_back(joint);

        auto jMesh1 = std::make_shared<SSEllipsoidMesh>(width[0]*1.1, width[0]*1.1, width[0]*1.1, "Mesh_Joint1", joint, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh1);

        MWMath::Point3D jointPosRel11 = MWMath::Point3D(0, 0, 0);
        auto joint1 = std::make_shared<SSJoint>("Joint_11", 
                                                jointPosRel11,           // Wo sitzt das Gelenk relativ zum Parent?
                                                MWMath::axisAngle({0,1,0}, 0.0),// Initiale Rotation
                                                joint,                 // Parent Body
                                                FJAngles[1],                  // Max Winkel
                                                MWMath::Point3D(0,1,0),// Drehachse: Negative Z-Achse!
                                                numTimeSteps);
        tissues.push_back(joint1);
        


        // ==============================================================================
        // SEGMENT 2 (Distal) - Hängt am Joint
        // ==============================================================================
        auto body2 = std::make_shared<SSBody>("Body_Seg2", MWMath::Point3D(L2,0,0), MWMath::RotMatrix3x3(), joint1);
        tissues.push_back(body2);
        auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2, 
             "Mesh_Seg2", body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.0, 0.8, 0.0));
        meshes.push_back(mesh2);
        // Torus 2
        auto mTorus2 = std::make_shared<SSTorusMesh>(mesh2->B * relTorusR[1], mesh2->B * relTorusr[1], 
             "Torus2", body2, MWMath::Point3D(0, -width[1]*0.95, L2 * relTorusPos[1]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
        meshes.push_back(mTorus2);


        MWMath::Point3D jointPosRel2 = MWMath::Point3D(L2, 0, 0);
        auto joint2 = std::make_shared<SSJoint>("Joint_2", 
                                             jointPosRel2,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body2,                 // Parent Body
                                             FJAngles[2],                  // Max Winkel
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
             "Torus3", body3, MWMath::Point3D(0, -width[2], L3 * relTorusPos[2]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
        meshes.push_back(mTorus3);


        MWMath::Point3D jointPosRel3 = MWMath::Point3D(L3, 0, 0);
        auto joint3 = std::make_shared<SSJoint>("Joint_3", 
                                             jointPosRel3,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body3,                 // Parent Body
                                             FJAngles[3],                  // Max Winkel
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
        int numPoints = cfg.muscleNumPoints[0];
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
    else if (currentScene == "OFINGER_SIMPLE4BT"){
        // --- PARAMETER ---
        double handLength = 1.8;
        std::vector<double> fLength = {0.3482 * handLength, 0.2027 * handLength, 0.1175 * handLength, 0.0882 * handLength};
        double L1 = fLength[0]; // Länge Segment 1
        double L2 = fLength[1]; // Länge Segment 2
        double L3 = fLength[2]; // Länge Segment 3
        double L4 = fLength[3]; // Länge Segment 4
        std::vector<double> width = {0.15, 0.1, 0.07, 0.05}; // Dicke
        std::vector<double> relTorusPos = {0.62, 0.45, 0.35};
        std::vector<double> relTorusR = {1.2, 1.2, 1.3};
        std::vector<double> relTorusr = {0.7, 0.7, 0.8};

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
             "Torus1", body1, MWMath::Point3D(0, -width[0]*0.8, L1 * relTorusPos[0]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mTorus1);

        MWMath::Point3D jointPosRel1 = MWMath::Point3D(L1, 0, 0);
        auto joint = std::make_shared<SSJoint>("Joint_1", 
                                             jointPosRel1,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body1,                 // Parent Body
                                             90.0,                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
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
             "Torus2", body2, MWMath::Point3D(0, -width[1]*0.95, L2 * relTorusPos[1]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
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
             "Torus3", body3, MWMath::Point3D(0, -width[2], L3 * relTorusPos[2]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
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
        int numPoints = cfg.muscleNumPoints[0];
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
    else if (currentScene == "OFINGER_SIMPLESCALE"){
        // --- PARAMETER ---
        double handLength = 1.8;
        std::vector<double> fLength = {0.3482 * handLength, 0.2027 * handLength, 0.1175 * handLength, 0.0882 * handLength};
        double L1 = fLength[0]; // Länge Segment 1
        double L2 = fLength[1]; // Länge Segment 2
        double L3 = fLength[2]; // Länge Segment 3
        double L4 = fLength[3]; // Länge Segment 4
        std::vector<double> width = {0.15, 0.1, 0.07, 0.05}; // Dicke
        std::vector<double> relTorusPos = {0.6, 0.2, 0.35};
        std::vector<double> relTorusR = {0.9, 1.0, 1.3};
        std::vector<double> relTorusr = {0.5, 0.7, 0.8};

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
             "Torus1", body1, MWMath::Point3D(0, -width[0]*0.8, L1 * relTorusPos[0]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(1, 0, 0));
        mTorus1->C = 2.0;
        meshes.push_back(mTorus1);

        MWMath::Point3D jointPosRel1 = MWMath::Point3D(L1, 0, 0);
        auto joint = std::make_shared<SSJoint>("Joint_1", 
                                             jointPosRel1,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body1,                 // Parent Body
                                             90.0,                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
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
             "Torus2", body2, MWMath::Point3D(0, -width[1]*0.95, L2 * relTorusPos[1]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
        mTorus2->C = 2.4;
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
             "Torus3", body3, MWMath::Point3D(0, -width[2], L3 * relTorusPos[2]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
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
        int numPoints = cfg.muscleNumPoints[0];
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
    else if (currentScene == "OFINGER_SIMPLE4T_B+S"){
        // --- PARAMETER ---
        double handLength = 1.8;
        std::vector<double> fLength = {0.3482 * handLength, 0.2027 * handLength, 0.1175 * handLength, 0.0882 * handLength};
        double L1 = fLength[0]; // Länge Segment 1
        double L2 = fLength[1]; // Länge Segment 2
        double L3 = fLength[2]; // Länge Segment 3
        double L4 = fLength[3]; // Länge Segment 4
        std::vector<double> width = {0.15, 0.1, 0.07, 0.05}; // Dicke
        std::vector<double> relTorusPos = {0.62, 0.45, 0.35};
        std::vector<double> relTorusR = {1.0, 1.0, 1.2};
        std::vector<double> relTorusr = {0.9, 0.9, 1.1};

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
             "Torus1", body1, MWMath::Point3D(0, -width[0]*0.8, L1 * relTorusPos[0]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mTorus1);

        MWMath::Point3D jointPosRel1 = MWMath::Point3D(L1, 0, 0);
        auto joint = std::make_shared<SSJoint>("Joint_1", 
                                                jointPosRel1,           // Wo sitzt das Gelenk relativ zum Parent?
                                                MWMath::RotMatrix3x3(),// Initiale Rotation
                                                body1,                 // Parent Body
                                                FJAngles[0],                  // Max Winkel
                                                MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                                numTimeSteps);
        tissues.push_back(joint);

        auto jMesh1 = std::make_shared<SSEllipsoidMesh>(width[0]*1.1, width[0]*1.1, width[0]*1.1, "Mesh_Joint1", joint, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh1);

        MWMath::Point3D jointPosRel11 = MWMath::Point3D(0, 0, 0);
        auto joint1 = std::make_shared<SSJoint>("Joint_11", 
                                                jointPosRel11,           // Wo sitzt das Gelenk relativ zum Parent?
                                                MWMath::axisAngle({0,1,0}, 0.0),// Initiale Rotation
                                                joint,                 // Parent Body
                                                FJAngles[1],                  // Max Winkel
                                                MWMath::Point3D(0,1,0),// Drehachse: Negative Z-Achse!
                                                numTimeSteps);
        tissues.push_back(joint1);
        


        // ==============================================================================
        // SEGMENT 2 (Distal) - Hängt am Joint
        // ==============================================================================
        auto body2 = std::make_shared<SSBody>("Body_Seg2", MWMath::Point3D(L2,0,0), MWMath::RotMatrix3x3(), joint1);
        tissues.push_back(body2);
        auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2, 
             "Mesh_Seg2", body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.0, 0.8, 0.0));
        meshes.push_back(mesh2);
        // Torus 2
        auto mTorus2 = std::make_shared<SSTorusMesh>(mesh2->B * relTorusR[1], mesh2->B * relTorusr[1], 
             "Torus2", body2, MWMath::Point3D(0, -width[1]*0.95, L2 * relTorusPos[1]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
        meshes.push_back(mTorus2);
        auto mTorus21 = std::make_shared<SSTorusMesh>(mesh2->B * relTorusR[1]*1.4, mesh2->B * relTorusr[1]*1.3, 
             "Torus21", body2, MWMath::Point3D(0, -width[1]*0.95, -L2 * relTorusPos[1]*0.2), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
        meshes.push_back(mTorus21);


        MWMath::Point3D jointPosRel2 = MWMath::Point3D(L2, 0, 0);
        auto joint2 = std::make_shared<SSJoint>("Joint_2", 
                                             jointPosRel2,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body2,                 // Parent Body
                                             FJAngles[2],                  // Max Winkel
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
             "Torus3", body3, MWMath::Point3D(0, -width[2], L3 * relTorusPos[2]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
        meshes.push_back(mTorus3);


        MWMath::Point3D jointPosRel3 = MWMath::Point3D(L3, 0, 0);
        auto joint3 = std::make_shared<SSJoint>("Joint_3", 
                                             jointPosRel3,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body3,                 // Parent Body
                                             FJAngles[3],                  // Max Winkel
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
        int numPoints = cfg.muscleNumPoints[0];
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

    else if (currentScene == "O_1BTorNextBelow_wJ"){
        // --- PARAMETER ---
        double handLength = 1.8;
        std::vector<double> fLength = {0.3482 * handLength, 0.2027 * handLength, 0.1175 * handLength, 0.0882 * handLength};
        double L1 = fLength[0]; // Länge Segment 1
        double L2 = fLength[1]; // Länge Segment 2
        double L3 = fLength[2]; // Länge Segment 3
        double L4 = fLength[3]; // Länge Segment 4
        std::vector<double> width = {0.15, 0.1, 0.07, 0.05}; // Dicke
        std::vector<double> relTorusPos =   {0.62, 0.3, 0.0, 0.1};
        std::vector<double> relTorusR =     {1.2, 1.7, 1.7, 1.7};
        std::vector<double> relTorusr =     {0.7, 1.2, 1.1, 1.1};
        std::vector<double> relTorusY =     {0.8, 0.7, 0.7, 0.7};

        // 1. ROOT
        rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), nullptr);
        tissues.push_back(rootSystem);

        // ==============================================================================
        // SEGMENT 1 (Proximal) - Fixiert an Root
        // ==============================================================================
        auto body1 = std::make_shared<SSBody>("Body_Seg1", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), rootSystem);
        tissues.push_back(body1);
        // Mesh 1
        auto mesh1 = std::make_shared<SSEllipsoidMesh>(width[0], width[0], L1, "Mesh_Seg1", body1, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), COLORF1);
        meshes.push_back(mesh1);

        MWMath::Point3D jointPosRel1 = MWMath::Point3D(L1, 0, 0);
        auto joint = std::make_shared<SSJoint>("Joint_1", 
                                                jointPosRel1,           // Wo sitzt das Gelenk relativ zum Parent?
                                                MWMath::RotMatrix3x3(),// Initiale Rotation
                                                body1,                 // Parent Body
                                                FJAngles[0],                  // Max Winkel
                                                MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                                numTimeSteps);
        tissues.push_back(joint);

        auto jMesh1 = std::make_shared<SSEllipsoidMesh>(width[0]*1.1, width[0]*1.1, width[0]*1.1, "Mesh_Joint1", joint, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh1);

        MWMath::Point3D jointPosRel11 = MWMath::Point3D(0, 0, 0);
        auto joint1 = std::make_shared<SSJoint>("Joint_11", 
                                                jointPosRel11,           // Wo sitzt das Gelenk relativ zum Parent?
                                                MWMath::axisAngle({0,1,0}, 0.0),// Initiale Rotation
                                                joint,                 // Parent Body
                                                FJAngles[1],                  // Max Winkel
                                                MWMath::Point3D(0,1,0),// Drehachse: Negative Z-Achse!
                                                numTimeSteps);
        tissues.push_back(joint1);
        


        // ==============================================================================
        // SEGMENT 2 (Distal) - Hängt am Joint
        // ==============================================================================
        auto body2 = std::make_shared<SSBody>("Body_Seg2", MWMath::Point3D(L2,0,0), MWMath::RotMatrix3x3(), joint1);
        tissues.push_back(body2);
        auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2, 
             "Mesh_Seg2", body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), COLORF2);
        meshes.push_back(mesh2);
        // Torus 2
        auto mTorus2 = std::make_shared<SSTorusMesh>(mesh2->B * relTorusR[1], mesh2->B * relTorusr[1], 
             "Torus2", body2, MWMath::Point3D(0, -width[1]*relTorusY[1], -L2 * relTorusPos[1]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
        meshes.push_back(mTorus2);


        MWMath::Point3D jointPosRel2 = MWMath::Point3D(L2, 0, 0);
        auto joint2 = std::make_shared<SSJoint>("Joint_2", 
                                             jointPosRel2,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body2,                 // Parent Body
                                             FJAngles[2],                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint2);

        auto jMesh2 = std::make_shared<SSEllipsoidMesh>(width[1]*1.1, width[1]*1.1, width[1]*1.1, "Mesh_Joint2", joint2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), COLORJOINT);
        meshes.push_back(jMesh2);

        // SEGMENT 3
        auto body3 = std::make_shared<SSBody>("Body_Seg3", MWMath::Point3D(L3,0,0), MWMath::RotMatrix3x3(), joint2);
        tissues.push_back(body3);
        // Mesh 3
        auto mesh3 = std::make_shared<SSEllipsoidMesh>(width[2], width[2], L3, 
             "Mesh_Seg3", body3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), 
             COLORF3);
        meshes.push_back(mesh3);
        // Torus 3
        auto mTorus3 = std::make_shared<SSTorusMesh>(mesh3->B * relTorusR[2], mesh3->B * relTorusr[2], 
             "Torus3", body3, MWMath::Point3D(0, -width[2]*relTorusY[2], -L3 * relTorusPos[2]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
        meshes.push_back(mTorus3);


        MWMath::Point3D jointPosRel3 = MWMath::Point3D(L3, 0, 0);
        auto joint3 = std::make_shared<SSJoint>("Joint_3", 
                                             jointPosRel3,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body3,                 // Parent Body
                                             FJAngles[3],                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint3);

        auto jMesh3 = std::make_shared<SSEllipsoidMesh>(width[2]*1.1, width[2]*1.1, width[2]*1.1, "Mesh_Joint3", joint3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), COLORJOINT);
        meshes.push_back(jMesh3);

        // SEGMENT 4
        auto body4 = std::make_shared<SSBody>("Body_Seg4", MWMath::Point3D(L4,0,0), MWMath::RotMatrix3x3(), joint3);
        tissues.push_back(body4);
        // Mesh 4
        auto mesh4 = std::make_shared<SSEllipsoidMesh>(width[3], width[3], L4, 
             "Mesh_Seg4", body4, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), 
             COLORF4);
        meshes.push_back(mesh4);

        auto mTorus4 = std::make_shared<SSTorusMesh>(mesh4->B * relTorusR[3], mesh4->B * relTorusr[3], 
             "Torus4", body4, MWMath::Point3D(0, -width[3]*relTorusY[3], -L4 * relTorusPos[3]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
        meshes.push_back(mTorus4);

        // ==============================================================================
        // INITIALES UPDATE
        // ==============================================================================
        // Damit alle Matrizen einmal berechnet werden (wichtig für Muskel-Punkte)
        for (auto& m : meshes) {m->InitializeMesh();}
        rootSystem->update(0); 


        // ==============================================================================
        // MUSKEL
        // ==============================================================================
        int numPoints = cfg.muscleNumPoints[0];
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
    else if (currentScene == "O_SmultiTorBelow_wJ"){
        // --- PARAMETER ---
        double handLength = 1.8;
        std::vector<double> fLength = {0.3482 * handLength, 0.2027 * handLength, 0.1175 * handLength, 0.0882 * handLength};
        double L1 = fLength[0]; // Länge Segment 1
        double L2 = fLength[1]; // Länge Segment 2
        double L3 = fLength[2]; // Länge Segment 3
        double L4 = fLength[3]; // Länge Segment 4
        std::vector<double> width = {0.15, 0.1, 0.07, 0.05}; // Dicke
        std::vector<double> relTorusPos =   {0.62, 0.4, 0.15, 0.1};
        std::vector<double> relTorusR =     {1.2, 1.3, 1.3, 1.7};
        std::vector<double> relTorusr =     {0.7, 0.75, 0.8, 1.1};
        std::vector<double> relTorusY =     {0.6, 0.7, 0.7, 0.7};

        // 1. ROOT
        rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), nullptr);
        tissues.push_back(rootSystem);

        // ==============================================================================
        // SEGMENT 1 (Proximal) - Fixiert an Root
        // ==============================================================================
        auto body1 = std::make_shared<SSBody>("Body_Seg1", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), rootSystem);
        tissues.push_back(body1);
        // Mesh 1
        auto mesh1 = std::make_shared<SSEllipsoidMesh>(width[0], width[0], L1, "Mesh_Seg1", body1, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), COLORF1);
        meshes.push_back(mesh1);

        MWMath::Point3D jointPosRel1 = MWMath::Point3D(L1, 0, 0);
        auto joint = std::make_shared<SSJoint>("Joint_1", 
                                                jointPosRel1,           // Wo sitzt das Gelenk relativ zum Parent?
                                                MWMath::RotMatrix3x3(),// Initiale Rotation
                                                body1,                 // Parent Body
                                                FJAngles[0],                  // Max Winkel
                                                MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                                numTimeSteps);
        tissues.push_back(joint);

        auto jMesh1 = std::make_shared<SSEllipsoidMesh>(width[0]*1.1, width[0]*1.1, width[0]*1.1, "Mesh_Joint1", joint, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh1);

        MWMath::Point3D jointPosRel11 = MWMath::Point3D(0, 0, 0);
        auto joint1 = std::make_shared<SSJoint>("Joint_11", 
                                                jointPosRel11,           // Wo sitzt das Gelenk relativ zum Parent?
                                                MWMath::axisAngle({0,1,0}, 0.0),// Initiale Rotation
                                                joint,                 // Parent Body
                                                FJAngles[1],                  // Max Winkel
                                                MWMath::Point3D(0,1,0),// Drehachse: Negative Z-Achse!
                                                numTimeSteps);
        tissues.push_back(joint1);
        


        // ==============================================================================
        // SEGMENT 2 (Distal) - Hängt am Joint
        // ==============================================================================
        auto body2 = std::make_shared<SSBody>("Body_Seg2", MWMath::Point3D(L2,0,0), MWMath::RotMatrix3x3(), joint1);
        tissues.push_back(body2);
        auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2, 
             "Mesh_Seg2", body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), COLORF2);
        meshes.push_back(mesh2);
        // Torus 2
        for (int i = 0; i < 3; i++) {
             auto mTorus = std::make_shared<SSTorusMesh>(mesh2->B * relTorusR[1], mesh2->B * relTorusr[1], 
             "Torus2_" + std::to_string(i), body2, MWMath::Point3D(0, -width[1]*relTorusY[1], -L2 * relTorusPos[1]*(1.0 - 1.1*relTorusr[1]*i)), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
             meshes.push_back(mTorus);
        }


        MWMath::Point3D jointPosRel2 = MWMath::Point3D(L2, 0, 0);
        auto joint2 = std::make_shared<SSJoint>("Joint_2", 
                                             jointPosRel2,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body2,                 // Parent Body
                                             FJAngles[2],                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint2);

        auto jMesh2 = std::make_shared<SSEllipsoidMesh>(width[1]*1.1, width[1]*1.1, width[1]*1.1, "Mesh_Joint2", joint2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), COLORJOINT);
        meshes.push_back(jMesh2);

        // SEGMENT 3
        auto body3 = std::make_shared<SSBody>("Body_Seg3", MWMath::Point3D(L3,0,0), MWMath::RotMatrix3x3(), joint2);
        tissues.push_back(body3);
        // Mesh 3
        auto mesh3 = std::make_shared<SSEllipsoidMesh>(width[2], width[2], L3, 
             "Mesh_Seg3", body3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), 
             COLORF3);
        meshes.push_back(mesh3);
        // Torus 3
        for (int i = 0; i < 2; i++) {
             auto mTorus = std::make_shared<SSTorusMesh>(mesh3->B * relTorusR[2], mesh3->B * relTorusr[2], 
             "Torus3_" + std::to_string(i), body3, MWMath::Point3D(0, -width[2]*relTorusY[2], -L3 * relTorusPos[2]*(1.0 - 3.0*relTorusr[2]*i)), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
             meshes.push_back(mTorus);
        }


        MWMath::Point3D jointPosRel3 = MWMath::Point3D(L3, 0, 0);
        auto joint3 = std::make_shared<SSJoint>("Joint_3", 
                                             jointPosRel3,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body3,                 // Parent Body
                                             FJAngles[3],                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint3);

        auto jMesh3 = std::make_shared<SSEllipsoidMesh>(width[2]*1.1, width[2]*1.1, width[2]*1.1, "Mesh_Joint3", joint3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), COLORJOINT);
        meshes.push_back(jMesh3);

        // SEGMENT 4
        auto body4 = std::make_shared<SSBody>("Body_Seg4", MWMath::Point3D(L4,0,0), MWMath::RotMatrix3x3(), joint3);
        tissues.push_back(body4);
        // Mesh 4
        auto mesh4 = std::make_shared<SSEllipsoidMesh>(width[3], width[3], L4, 
             "Mesh_Seg4", body4, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), 
             COLORF4);
        meshes.push_back(mesh4);

        auto mTorus4 = std::make_shared<SSTorusMesh>(mesh4->B * relTorusR[3], mesh4->B * relTorusr[3], 
             "Torus4", body4, MWMath::Point3D(0, -width[3]*relTorusY[3], -L4 * relTorusPos[3]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
        meshes.push_back(mTorus4);

        // ==============================================================================
        // INITIALES UPDATE
        // ==============================================================================
        // Damit alle Matrizen einmal berechnet werden (wichtig für Muskel-Punkte)
        for (auto& m : meshes) {m->InitializeMesh();}
        rootSystem->update(0); 


        // ==============================================================================
        // MUSKEL
        // ==============================================================================
        int numPoints = cfg.muscleNumPoints[0];
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
    
    else if (currentScene == "OFINGER_SIMPLE_BigTorus_2DJoint"){
        // --- PARAMETER ---
        double handLength = 1.8*GEOMETRYSCALER;
        float GS = GEOMETRYSCALER;
        std::vector<double> fLength = {0.3482 * handLength, 0.2027 * handLength, 0.1175 * handLength, 0.0882 * handLength};
        double L1 = fLength[0]; // Länge Segment 1
        double L2 = fLength[1]; // Länge Segment 2
        double L3 = fLength[2]; // Länge Segment 3
        double L4 = fLength[3]; // Länge Segment 4
        double rWF = 0.7; // reduce width factor
        std::vector<double> width = {0.15*GS*rWF, 0.1*GS*rWF, 0.07*GS*rWF, 0.05*GS*rWF}; // Dicke
        std::vector<double> relTorusPos = {0.6, 0.42, 0.32};
        std::vector<double> relTorusR = {1.85/rWF, 2.0/rWF, 2.3/rWF}; // {1.85, 2.0, 2.3}
        std::vector<double> relTorusr = {1.0/rWF, 1.0/rWF, 1.2/rWF}; // {1.0, 1.0, 1.2}

        /* std::vector<double> width = {0.15*GS, 0.1*GS, 0.07*GS, 0.05*GS}; // Dicke
        std::vector<double> relTorusPos = {0.6, 0.42, 0.35};
        std::vector<double> relTorusR = {1.85, 2.0, 2.3};
        std::vector<double> relTorusr = {1.0, 1.0, 1.2}; */
        double off = 0.0; 

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
             "Torus1", body1, MWMath::Point3D(0, -width[0]*off, L1 * relTorusPos[0]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mTorus1);

        MWMath::Point3D jointPosRel1 = MWMath::Point3D(L1, 0, 0);
        auto joint = std::make_shared<SSJoint>("Joint_1", 
                                             jointPosRel1,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body1,                 // Parent Body
                                             FJAngles[0],                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint);
        
        double jointSize1 = width[0]*1.1/rWF;
        auto jMesh1 = std::make_shared<SSEllipsoidMesh>(jointSize1, jointSize1, jointSize1, "Mesh_Joint1", joint, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh1);
        
        MWMath::Point3D jointPosRel11 = MWMath::Point3D(0, 0, 0);
        auto joint1 = std::make_shared<SSJoint>("Joint_11", 
                                             jointPosRel11,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::axisAngle({0,1,0}, 0.0),// Initiale Rotation
                                             joint,                 // Parent Body
                                              FJAngles[1],                  // Max Winkel
                                             MWMath::Point3D(0,1,0),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint1);
        


        // ==============================================================================
        // SEGMENT 2 (Distal) - Hängt am Joint
        // ==============================================================================
        auto body2 = std::make_shared<SSBody>("Body_Seg2", MWMath::Point3D(L2,0,0), MWMath::RotMatrix3x3(), joint1);
        tissues.push_back(body2);
        auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2, 
             "Mesh_Seg2", body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.0, 0.8, 0.0));
        meshes.push_back(mesh2);
        // Torus 2
        auto mTorus2 = std::make_shared<SSTorusMesh>(mesh2->B * relTorusR[1], mesh2->B * relTorusr[1], 
             "Torus2", body2, MWMath::Point3D(0, -width[1]*off, L2 * relTorusPos[1]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
        meshes.push_back(mTorus2);


        MWMath::Point3D jointPosRel2 = MWMath::Point3D(L2, 0, 0);
        auto joint2 = std::make_shared<SSJoint>("Joint_2", 
                                             jointPosRel2,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body2,                 // Parent Body
                                              FJAngles[2],                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint2);
        
        double jointSize2 = width[1]*1.1/rWF;
        auto jMesh2 = std::make_shared<SSEllipsoidMesh>(jointSize2, jointSize2, jointSize2, "Mesh_Joint2", joint2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
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
        auto mTorus3 = std::make_shared<SSTorusMesh>(width[2] * relTorusR[2], width[2] * relTorusr[2], 
             "Torus3", body3, MWMath::Point3D(0, -width[2]*off, L3 * relTorusPos[2]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
        meshes.push_back(mTorus3);


        MWMath::Point3D jointPosRel3 = MWMath::Point3D(L3, 0, 0);
        auto joint3 = std::make_shared<SSJoint>("Joint_3", 
                                             jointPosRel3,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body3,                 // Parent Body
                                              FJAngles[3],                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        //joint3->AngleSteps = jointAnglesPrescribed;
        tissues.push_back(joint3);

        double jointSize3 = width[2]*1.1/rWF;
        auto jMesh3 = std::make_shared<SSEllipsoidMesh>(jointSize3, jointSize3, jointSize3, "Mesh_Joint3", joint3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
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
        int numPoints = cfg.muscleNumPoints[0];
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
    else if (currentScene == "OFINGER_SIMPLE_OnlyTorus"){
        // --- PARAMETER ---
        double handLength = 1.8*GEOMETRYSCALER;
        float GS = GEOMETRYSCALER;
        std::vector<double> fLength = {0.3482 * handLength, 0.2027 * handLength, 0.1175 * handLength, 0.0882 * handLength};
        double L1 = fLength[0]; // Länge Segment 1
        double L2 = fLength[1]; // Länge Segment 2
        double L3 = fLength[2]; // Länge Segment 3
        double L4 = fLength[3]; // Länge Segment 4
        double rWF = 0.7; // reduce width factor
        std::vector<double> width = {0.15*GS*rWF, 0.1*GS*rWF, 0.07*GS*rWF, 0.05*GS*rWF}; // Dicke
        std::vector<double> relTorusPos = {0.6, 0.42, 0.32};
        std::vector<double> relTorusR = {1.85/rWF, 2.0/rWF, 2.3/rWF}; // {1.85, 2.0, 2.3}
        std::vector<double> relTorusr = {1.0/rWF, 1.0/rWF, 1.2/rWF}; // {1.0, 1.0, 1.2}

        /* std::vector<double> width = {0.15*GS, 0.1*GS, 0.07*GS, 0.05*GS}; // Dicke
        std::vector<double> relTorusPos = {0.6, 0.42, 0.35};
        std::vector<double> relTorusR = {1.85, 2.0, 2.3};
        std::vector<double> relTorusr = {1.0, 1.0, 1.2}; */
        bool bShowBody = false;
        double off = 0.0; 

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
        if (bShowBody) meshes.push_back(mesh1);
        // Torus 1
        auto mTorus1 = std::make_shared<SSTorusMesh>(mesh1->B * relTorusR[0], mesh1->B * relTorusr[0], 
             "Torus1", body1, MWMath::Point3D(0, -width[0]*off, L1 * relTorusPos[0]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mTorus1);

        MWMath::Point3D jointPosRel1 = MWMath::Point3D(L1, 0, 0);
        auto joint = std::make_shared<SSJoint>("Joint_1", 
                                             jointPosRel1,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body1,                 // Parent Body
                                             FJAngles[0],                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint);
        
        double jointSize1 = width[0]*1.1/rWF;
        auto jMesh1 = std::make_shared<SSEllipsoidMesh>(jointSize1, jointSize1, jointSize1, "Mesh_Joint1", joint, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh1);
        
        MWMath::Point3D jointPosRel11 = MWMath::Point3D(0, 0, 0);
        auto joint1 = std::make_shared<SSJoint>("Joint_11", 
                                             jointPosRel11,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::axisAngle({0,1,0}, 0.0),// Initiale Rotation
                                             joint,                 // Parent Body
                                              FJAngles[1],                  // Max Winkel
                                             MWMath::Point3D(0,1,0),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint1);
        


        // ==============================================================================
        // SEGMENT 2 (Distal) - Hängt am Joint
        // ==============================================================================
        auto body2 = std::make_shared<SSBody>("Body_Seg2", MWMath::Point3D(L2,0,0), MWMath::RotMatrix3x3(), joint1);
        tissues.push_back(body2);
        auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2, 
             "Mesh_Seg2", body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.0, 0.8, 0.0));
        if (bShowBody) meshes.push_back(mesh2);
        // Torus 2
        auto mTorus2 = std::make_shared<SSTorusMesh>(mesh2->B * relTorusR[1], mesh2->B * relTorusr[1], 
             "Torus2", body2, MWMath::Point3D(0, -width[1]*off, L2 * relTorusPos[1]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
        meshes.push_back(mTorus2);


        MWMath::Point3D jointPosRel2 = MWMath::Point3D(L2, 0, 0);
        auto joint2 = std::make_shared<SSJoint>("Joint_2", 
                                             jointPosRel2,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body2,                 // Parent Body
                                              FJAngles[2],                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint2);
        
        double jointSize2 = width[1]*1.1/rWF;
        auto jMesh2 = std::make_shared<SSEllipsoidMesh>(jointSize2, jointSize2, jointSize2, "Mesh_Joint2", joint2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh2);

        // SEGMENT 3
        auto body3 = std::make_shared<SSBody>("Body_Seg3", MWMath::Point3D(L3,0,0), MWMath::RotMatrix3x3(), joint2);
        tissues.push_back(body3);
        // Mesh 3
        auto mesh3 = std::make_shared<SSEllipsoidMesh>(width[2], width[2], L3, 
             "Mesh_Seg3", body3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), 
             MWMath::Point3D(0.0, 0.8, 0.8));
        if (bShowBody) meshes.push_back(mesh3);
        // Torus 3
        auto mTorus3 = std::make_shared<SSTorusMesh>(width[2] * relTorusR[2], width[2] * relTorusr[2], 
             "Torus3", body3, MWMath::Point3D(0, -width[2]*off, L3 * relTorusPos[2]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
        meshes.push_back(mTorus3);


        MWMath::Point3D jointPosRel3 = MWMath::Point3D(L3, 0, 0);
        auto joint3 = std::make_shared<SSJoint>("Joint_3", 
                                             jointPosRel3,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body3,                 // Parent Body
                                              FJAngles[3],                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        //joint3->AngleSteps = jointAnglesPrescribed;
        tissues.push_back(joint3);

        double jointSize3 = width[2]*1.1/rWF;
        auto jMesh3 = std::make_shared<SSEllipsoidMesh>(jointSize3, jointSize3, jointSize3, "Mesh_Joint3", joint3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh3);

        // SEGMENT 4
        auto body4 = std::make_shared<SSBody>("Body_Seg4", MWMath::Point3D(L4,0,0), MWMath::RotMatrix3x3(), joint3);
        tissues.push_back(body4);
        // Mesh 4
        auto mesh4 = std::make_shared<SSEllipsoidMesh>(width[3], width[3], L4, 
             "Mesh_Seg4", body4, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), 
             MWMath::Point3D(0.9, 0.5, 0.3));
        if (bShowBody) meshes.push_back(mesh4);

        // ==============================================================================
        // INITIALES UPDATE
        // ==============================================================================
        // Damit alle Matrizen einmal berechnet werden (wichtig für Muskel-Punkte)
        for (auto& m : meshes) {m->InitializeMesh();}
        rootSystem->update(0); 


        // ==============================================================================
        // MUSKEL
        // ==============================================================================
        int numPoints = cfg.muscleNumPoints[0];
        MWMath::Point3D startOffset = MWMath::Point3D(0.0, -mesh1->B * 1.1, 0.0);
        MWMath::Point3D endOffset = MWMath::Point3D(mesh4->C, -mesh4->B * 1.1, 0.0);
        SSMuscle* flexor = new SSMuscle("Flexor", numPoints, 
            body1.get(), startOffset, 
            body4.get(), endOffset);

        // Alle Hindernisse
        for(auto& m : meshes) {
            if (dynamic_cast<SSTorusMesh*>(m.get()) || m->bIsJointMesh){
                flexor->meshPtrs.push_back(m.get());
            }
        }

        flexor->createMusclePoints();
        flexor->updateMusclePointsParents();
        muscles.push_back(flexor);
    }
    else if (currentScene == "OFINGER_SIMPLE_OnlyTorusSmall"){
        // --- PARAMETER ---
        double handLength = 1.8*GEOMETRYSCALER;
        float GS = GEOMETRYSCALER;
        std::vector<double> fLength = {0.3482 * handLength, 0.2027 * handLength, 0.1175 * handLength, 0.0882 * handLength};
        double L1 = fLength[0]; // Länge Segment 1
        double L2 = fLength[1]; // Länge Segment 2
        double L3 = fLength[2]; // Länge Segment 3
        double L4 = fLength[3]; // Länge Segment 4
        double rWF = 1.0; // reduce width factor
        std::vector<double> width = {0.15*GS*rWF, 0.1*GS*rWF, 0.07*GS*rWF, 0.05*GS*rWF}; // Dicke
        std::vector<double> relTorusPos = {0.6, 0.42, 0.32};
        std::vector<double> relTorusR = {1.3/rWF, 1.3/rWF, 1.6/rWF}; // {1.85, 2.0, 2.3}
        std::vector<double> relTorusr = {1.1/rWF, 1.1/rWF, 1.4/rWF}; // {1.0, 1.0, 1.2}

        /* std::vector<double> width = {0.15*GS, 0.1*GS, 0.07*GS, 0.05*GS}; // Dicke
        std::vector<double> relTorusPos = {0.6, 0.42, 0.35};
        std::vector<double> relTorusR = {1.85, 2.0, 2.3};
        std::vector<double> relTorusr = {1.0, 1.0, 1.2}; */
        bool bShowBody = true;
        double off = 0.8; 

        // 1. ROOT
        rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), nullptr);
        tissues.push_back(rootSystem);

        // ==============================================================================
        // SEGMENT 1 (Proximal) - Fixiert an Root
        // ==============================================================================
        auto body1 = std::make_shared<SSBody>("Body_Seg1", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), rootSystem);
        tissues.push_back(body1);
        // Mesh 1
        auto mesh1 = std::make_shared<SSEllipsoidMesh>(width[0], width[0], L1, "Mesh_Seg1", body1, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), COLORFDEACTIVE);
        if (bShowBody) meshes.push_back(mesh1);
        // Torus 1
        auto mTorus1 = std::make_shared<SSTorusMesh>(mesh1->B * relTorusR[0], mesh1->B * relTorusr[0], 
             "Torus1", body1, MWMath::Point3D(0, -width[0]*off, L1 * relTorusPos[0]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mTorus1);

        MWMath::Point3D jointPosRel1 = MWMath::Point3D(L1, 0, 0);
        auto joint = std::make_shared<SSJoint>("Joint_1", 
                                             jointPosRel1,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body1,                 // Parent Body
                                             FJAngles[0],                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint);
        
        double jointSize1 = width[0]*1.1/rWF;
        auto jMesh1 = std::make_shared<SSEllipsoidMesh>(jointSize1, jointSize1, jointSize1, "Mesh_Joint1", joint, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh1);
        
        MWMath::Point3D jointPosRel11 = MWMath::Point3D(0, 0, 0);
        auto joint1 = std::make_shared<SSJoint>("Joint_11", 
                                             jointPosRel11,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::axisAngle({0,1,0}, 0.0),// Initiale Rotation
                                             joint,                 // Parent Body
                                              FJAngles[1],                  // Max Winkel
                                             MWMath::Point3D(0,1,0),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint1);
        


        // ==============================================================================
        // SEGMENT 2 (Distal) - Hängt am Joint
        // ==============================================================================
        auto body2 = std::make_shared<SSBody>("Body_Seg2", MWMath::Point3D(L2,0,0), MWMath::RotMatrix3x3(), joint1);
        tissues.push_back(body2);
        auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2, 
             "Mesh_Seg2", body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), COLORFDEACTIVE);
        if (bShowBody) meshes.push_back(mesh2);
        // Torus 2
        auto mTorus2 = std::make_shared<SSTorusMesh>(mesh2->B * relTorusR[1], mesh2->B * relTorusr[1], 
             "Torus2", body2, MWMath::Point3D(0, -width[1]*off, L2 * relTorusPos[1]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
        meshes.push_back(mTorus2);


        MWMath::Point3D jointPosRel2 = MWMath::Point3D(L2, 0, 0);
        auto joint2 = std::make_shared<SSJoint>("Joint_2", 
                                             jointPosRel2,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body2,                 // Parent Body
                                              FJAngles[2],                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint2);
        
        double jointSize2 = width[1]*1.1/rWF;
        auto jMesh2 = std::make_shared<SSEllipsoidMesh>(jointSize2, jointSize2, jointSize2, "Mesh_Joint2", joint2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh2);

        // SEGMENT 3
        auto body3 = std::make_shared<SSBody>("Body_Seg3", MWMath::Point3D(L3,0,0), MWMath::RotMatrix3x3(), joint2);
        tissues.push_back(body3);
        // Mesh 3
        auto mesh3 = std::make_shared<SSEllipsoidMesh>(width[2], width[2], L3, 
             "Mesh_Seg3", body3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), 
             COLORFDEACTIVE);
        if (bShowBody) meshes.push_back(mesh3);
        // Torus 3
        auto mTorus3 = std::make_shared<SSTorusMesh>(width[2] * relTorusR[2], width[2] * relTorusr[2], 
             "Torus3", body3, MWMath::Point3D(0, -width[2]*off, L3 * relTorusPos[2]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
        meshes.push_back(mTorus3);


        MWMath::Point3D jointPosRel3 = MWMath::Point3D(L3, 0, 0);
        auto joint3 = std::make_shared<SSJoint>("Joint_3", 
                                             jointPosRel3,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body3,                 // Parent Body
                                              FJAngles[3],                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        //joint3->AngleSteps = jointAnglesPrescribed;
        tissues.push_back(joint3);

        double jointSize3 = width[2]*1.1/rWF;
        auto jMesh3 = std::make_shared<SSEllipsoidMesh>(jointSize3, jointSize3, jointSize3, "Mesh_Joint3", joint3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh3);

        // SEGMENT 4
        auto body4 = std::make_shared<SSBody>("Body_Seg4", MWMath::Point3D(L4,0,0), MWMath::RotMatrix3x3(), joint3);
        tissues.push_back(body4);
        // Mesh 4
        auto mesh4 = std::make_shared<SSEllipsoidMesh>(width[3], width[3], L4, 
             "Mesh_Seg4", body4, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), 
             COLORFDEACTIVE);
        if (bShowBody) meshes.push_back(mesh4);

        // ==============================================================================
        // INITIALES UPDATE
        // ==============================================================================
        // Damit alle Matrizen einmal berechnet werden (wichtig für Muskel-Punkte)
        for (auto& m : meshes) {m->InitializeMesh();}
        rootSystem->update(0); 


        // ==============================================================================
        // MUSKEL
        // ==============================================================================
        int numPoints = cfg.muscleNumPoints[0];
        MWMath::Point3D startOffset = MWMath::Point3D(0.0, -mesh1->B * 1.1, 0.0);
        MWMath::Point3D endOffset = MWMath::Point3D(mesh4->C, -mesh4->B * 1.1, 0.0);
        SSMuscle* flexor = new SSMuscle("Flexor", numPoints, 
            body1.get(), startOffset, 
            body4.get(), endOffset);

        // Alle Hindernisse
        for(auto& m : meshes) {
            if (dynamic_cast<SSTorusMesh*>(m.get()) || m->bIsJointMesh){
                flexor->meshPtrs.push_back(m.get());
            }
        }

        flexor->createMusclePoints();
        flexor->updateMusclePointsParents();
        muscles.push_back(flexor);
    }
    else if (currentScene == "OFINGER_SIMPLE_BigTorus_2DJointCylinder"){
        // --- PARAMETER ---
        double handLength = 1.8;
        std::vector<double> fLength = {0.3482 * handLength, 0.2027 * handLength, 0.1175 * handLength, 0.0882 * handLength};
        double L1 = fLength[0]; // Länge Segment 1
        double L2 = fLength[1]; // Länge Segment 2
        double L3 = fLength[2]; // Länge Segment 3
        double L4 = fLength[3]; // Länge Segment 4
        std::vector<double> width = {0.15, 0.1, 0.07, 0.05}; // Dicke
        std::vector<double> relTorusPos = {0.6, 0.42, 0.32};
        std::vector<double> relTorusR = {1.85, 2.0, 2.3};
        std::vector<double> relTorusr = {1.0, 1.0, 1.2};
        double off = 0.0; 

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
             "Torus1", body1, MWMath::Point3D(0, -width[0]*off, L1 * relTorusPos[0]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mTorus1);

        MWMath::Point3D jointPosRel1 = MWMath::Point3D(L1, 0, 0);
        auto joint = std::make_shared<SSJoint>("Joint_1", 
                                             jointPosRel1,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body1,                 // Parent Body
                                             FJAngles[0],                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint);

        auto jMesh1 = std::make_shared<SSCylinderMesh>(width[0]*1.1, width[0]*1.1, "Mesh_Joint1", joint, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh1);
        
        MWMath::Point3D jointPosRel11 = MWMath::Point3D(0, 0, 0);
        auto joint1 = std::make_shared<SSJoint>("Joint_11", 
                                             jointPosRel11,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::axisAngle({0,1,0}, 0.0),// Initiale Rotation
                                             joint,                 // Parent Body
                                             FJAngles[1],                  // Max Winkel
                                             MWMath::Point3D(0,1,0),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint1);
        


        // ==============================================================================
        // SEGMENT 2 (Distal) - Hängt am Joint
        // ==============================================================================
        auto body2 = std::make_shared<SSBody>("Body_Seg2", MWMath::Point3D(L2,0,0), MWMath::RotMatrix3x3(), joint1);
        tissues.push_back(body2);
        auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2, 
             "Mesh_Seg2", body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.0, 0.8, 0.0));
        meshes.push_back(mesh2);
        // Torus 2
        auto mTorus2 = std::make_shared<SSTorusMesh>(mesh2->B * relTorusR[1], mesh2->B * relTorusr[1], 
             "Torus2", body2, MWMath::Point3D(0, -width[1]*off, L2 * relTorusPos[1]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
        meshes.push_back(mTorus2);


        MWMath::Point3D jointPosRel2 = MWMath::Point3D(L2, 0, 0);
        auto joint2 = std::make_shared<SSJoint>("Joint_2", 
                                             jointPosRel2,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body2,                 // Parent Body
                                             FJAngles[2],                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint2);

        auto jMesh2 = std::make_shared<SSCylinderMesh>(width[1]*1.1, width[1]*1.1, "Mesh_Joint2", joint2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
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
             "Torus3", body3, MWMath::Point3D(0, -width[2]*off, L3 * relTorusPos[2]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
        meshes.push_back(mTorus3);


        MWMath::Point3D jointPosRel3 = MWMath::Point3D(L3, 0, 0);
        auto joint3 = std::make_shared<SSJoint>("Joint_3", 
                                             jointPosRel3,           // Wo sitzt das Gelenk relativ zum Parent?
                                             MWMath::RotMatrix3x3(),// Initiale Rotation
                                             body3,                 // Parent Body
                                             FJAngles[3],                  // Max Winkel
                                             MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse!
                                             numTimeSteps);
        tissues.push_back(joint3);

        auto jMesh3 = std::make_shared<SSCylinderMesh>(width[2]*1.1, width[2]*1.1, "Mesh_Joint3", joint3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
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
        int numPoints = cfg.muscleNumPoints[0];
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
    
    else if (currentScene == "OFINGER_2Segs"){
        // --- PARAMETER ---
        double handLength = 1.8 * GEOMETRYSCALER;
        float GS = GEOMETRYSCALER;
        
        // Wir nehmen die Werte der ehemals "hintersten" Segmente (L3 und L4)
        double L_Prox = 0.1175 * handLength; // Ursprünglich L3
        double L_Dist = 0.0882 * handLength; // Ursprünglich L4
        
        double w_Prox = 0.07 * GS;           // Ursprünglich width[2]
        double w_Dist = 0.05 * GS;           // Ursprünglich width[3]
        
        double torusRelPos = 0.32;           // Ursprünglich relTorusPos[2]
        double torusRelR = 2.3;              // Ursprünglich relTorusR[2]
        double torusRelr = 1.2;              // Ursprünglich relTorusr[2]
        double off = 0.0; 

        // 1. ROOT
        rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), nullptr);
        tissues.push_back(rootSystem);

        // ==============================================================================
        // SEGMENT 1 (Proximal) - Fixiert an Root
        // ==============================================================================
        auto bodyProx = std::make_shared<SSBody>("Body_Prox", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), rootSystem);
        tissues.push_back(bodyProx);
        
        // Mesh Proximal
        auto meshProx = std::make_shared<SSEllipsoidMesh>(w_Prox, w_Prox, L_Prox, 
             "Mesh_Prox", bodyProx, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), 
             MWMath::Point3D(0.0, 0.8, 0.8));
        meshes.push_back(meshProx);
        
        // Torus Proximal (führt den Muskel)
        auto mTorusProx = std::make_shared<SSTorusMesh>(meshProx->B * torusRelR, meshProx->B * torusRelr, 
             "Torus_Prox", bodyProx, MWMath::Point3D(0, -w_Prox*off, L_Prox * torusRelPos), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
        meshes.push_back(mTorusProx);

        // ==============================================================================
        // GELENK (Sitzt am Ende von Segment 1)
        // ==============================================================================
        MWMath::Point3D jointPosRel = MWMath::Point3D(L_Prox, 0, 0);
        auto joint = std::make_shared<SSJoint>("Joint_1", 
                                               jointPosRel,            // Position relativ zum Parent
                                               MWMath::RotMatrix3x3(), // Initiale Rotation
                                               bodyProx,               // Parent Body
                                               90.0, // Max Winkel aus Config
                                               MWMath::Point3D(0,0,-1),// Drehachse: Negative Z-Achse
                                               numTimeSteps);
        
        // Falls du vorgegebene Winkel-Kurven benutzt (auskommentieren falls nicht benötigt)
        // joint->AngleSteps = jointAnglesPrescribed; 
        
        tissues.push_back(joint);

        // Visuelle Kugel für das Gelenk
        auto jMesh = std::make_shared<SSEllipsoidMesh>(w_Prox*1.1, w_Prox*1.1, w_Prox*1.1, "Mesh_Joint", joint, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh);

        // ==============================================================================
        // SEGMENT 2 (Distal) - Hängt am Gelenk
        // ==============================================================================
        auto bodyDist = std::make_shared<SSBody>("Body_Dist", MWMath::Point3D(L_Dist,0,0), MWMath::RotMatrix3x3(), joint);
        tissues.push_back(bodyDist);
        
        // Mesh Distal
        auto meshDist = std::make_shared<SSEllipsoidMesh>(w_Dist, w_Dist, L_Dist, 
             "Mesh_Dist", bodyDist, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), 
             MWMath::Point3D(0.9, 0.5, 0.3));
        meshes.push_back(meshDist);

        // ==============================================================================
        // INITIALES UPDATE
        // ==============================================================================
        // Matrizen berechnen (wichtig für Muskel-Punkte Setup)
        for (auto& m : meshes) {m->InitializeMesh();}
        rootSystem->update(0); 

        // ==============================================================================
        // MUSKEL
        // ==============================================================================
        int numPoints = cfg.muscleNumPoints.empty() ? 25 : cfg.muscleNumPoints[0];
        
        // Ansatzpunkte unten am Knochen
        MWMath::Point3D startOffset = MWMath::Point3D(-handLength*0.7, -0.15*1.1, 0.0);
        MWMath::Point3D endOffset = MWMath::Point3D(meshDist->C, -meshDist->B * 1.1, 0.0);
        
        SSMuscle* flexor = new SSMuscle("Flexor", numPoints, 
            bodyProx.get(), startOffset, 
            bodyDist.get(), endOffset);

        // Alle vorhandenen Meshes als potentielle Hindernisse anmelden
        for(auto& m : meshes) {
             flexor->meshPtrs.push_back(m.get());
        }

        flexor->createMusclePoints();
        flexor->updateMusclePointsParents();
        muscles.push_back(flexor);
    }
    else if (currentScene == "OFINGER_3Segs"){
        // --- PARAMETER ---
        double handLength = 1.8 * GEOMETRYSCALER;
        float GS = GEOMETRYSCALER;
        
        // Die Längen für die letzten 3 Segmente
        std::vector<double> fLength = {0.3482 * handLength, 0.2027 * handLength, 0.1175 * handLength, 0.0882 * handLength};
        double L2 = fLength[1]; // Länge Segment 2 (Proximal in dieser Szene)
        double L3 = fLength[2]; // Länge Segment 3 (Mitte)
        double L4 = fLength[3]; // Länge Segment 4 (Distal)
        
        std::vector<double> width = {0.15*GS, 0.1*GS, 0.07*GS, 0.05*GS}; // Dicke
        std::vector<double> relTorusPos = {0.6, 0.42, 0.32};
        std::vector<double> relTorusR = {1.85, 2.0, 2.3};
        std::vector<double> relTorusr = {1.0, 1.0, 1.2};
        double off = 0.0; 

        // 1. ROOT
        rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), nullptr);
        tissues.push_back(rootSystem);

        // ==============================================================================
        // SEGMENT 2 (Proximal) - Fixiert an Root
        // ==============================================================================
        auto body2 = std::make_shared<SSBody>("Body_Seg2", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), rootSystem);
        tissues.push_back(body2);
        
        // Mesh 2
        auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2, 
             "Mesh_Seg2", body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.0, 0.8, 0.0));
        meshes.push_back(mesh2);
        
        // Torus 2
        /* auto mTorus2 = std::make_shared<SSTorusMesh>(mesh2->B * relTorusR[1], mesh2->B * relTorusr[1], 
             "Torus2", body2, MWMath::Point3D(0, -width[1]*off, L2 * relTorusPos[1]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
        meshes.push_back(mTorus2); */

        // ==============================================================================
        // GELENK 2 (Verbindet Seg 2 und Seg 3)
        // ==============================================================================
        MWMath::Point3D jointPosRel2 = MWMath::Point3D(L2, 0, 0);
        auto joint2 = std::make_shared<SSJoint>("Joint_2", 
                                               jointPosRel2,           
                                               MWMath::RotMatrix3x3(),
                                               body2,                  
                                               0.0, 
                                               MWMath::Point3D(0,0,-1), // Drehachse
                                               numTimeSteps);
        tissues.push_back(joint2);

        //auto jMesh2 = std::make_shared<SSEllipsoidMesh>(width[1]*1.1, width[1]*1.1, width[1]*1.1, "Mesh_Joint2", joint2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        //auto jMesh2 = std::make_shared<SSCylinderMesh>(width[1]*1.1, width[1]*1.1, "Mesh_Joint2", joint2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        //meshes.push_back(jMesh2);

        // ==============================================================================
        // SEGMENT 3 (Mitte) - Hängt am Joint 2
        // ==============================================================================
        auto body3 = std::make_shared<SSBody>("Body_Seg3", MWMath::Point3D(L3,0,0), MWMath::RotMatrix3x3(), joint2);
        tissues.push_back(body3);
        
        // Mesh 3
        /* auto mesh3 = std::make_shared<SSEllipsoidMesh>(width[2], width[2], L3, 
             "Mesh_Seg3", body3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), 
             MWMath::Point3D(0.0, 0.8, 0.8));
        meshes.push_back(mesh3); */
        
        // Torus 3
        auto mTorus3 = std::make_shared<SSTorusMesh>( width[2] * relTorusR[2], width[2] * relTorusr[2], 
             "Torus3", body3, MWMath::Point3D(0, -width[2]*off, L3 * relTorusPos[2]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
        meshes.push_back(mTorus3);

        // ==============================================================================
        // GELENK 3 (Verbindet Seg 3 und Seg 4)
        // ==============================================================================
        MWMath::Point3D jointPosRel3 = MWMath::Point3D(L3, 0, 0);
        auto joint3 = std::make_shared<SSJoint>("Joint_3", 
                                               jointPosRel3,           
                                               MWMath::RotMatrix3x3(),
                                               body3,                  
                                               90.0, 
                                               MWMath::Point3D(0,0,-1), // Drehachse
                                               numTimeSteps);
        
        // Falls du eine bestimmte Trajektorie für das vordere Gelenk hast:
        // joint3->AngleSteps = jointAnglesPrescribed; 
        
        tissues.push_back(joint3);

        //auto jMesh3 = std::make_shared<SSEllipsoidMesh>(width[2]*1.1, width[2]*1.1, width[2]*1.1, "Mesh_Joint3", joint3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        //auto jMesh3 = std::make_shared<SSCylinderMesh>(width[2]*1.1, width[2]*1.1, "Mesh_Joint3", joint3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        //meshes.push_back(jMesh3);

        // ==============================================================================
        // SEGMENT 4 (Distal) - Hängt am Joint 3
        // ==============================================================================
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
        // Matrizen berechnen (wichtig für Muskel-Punkte Setup)
        for (auto& m : meshes) {m->InitializeMesh();}
        rootSystem->update(0); 

        // ==============================================================================
        // MUSKEL (Spannung von Seg 2 bis Seg 4)
        // ==============================================================================
        int numPoints = cfg.muscleNumPoints.empty() ? 25 : cfg.muscleNumPoints[0];
        
        // Ansatzpunkte
        MWMath::Point3D startOffset = MWMath::Point3D(-handLength*0.7, -0.15*1.1, 0.0);
        MWMath::Point3D endOffset = MWMath::Point3D(mesh4->C, -mesh4->B * 1.1, 0.0);
        
        SSMuscle* flexor = new SSMuscle("Flexor", numPoints, 
            body2.get(), startOffset, 
            body4.get(), endOffset);

        // Alle vorhandenen Meshes als potentielle Hindernisse anmelden
        for(auto& m : meshes) {
             flexor->meshPtrs.push_back(m.get());
        }

        flexor->createMusclePoints();
        flexor->updateMusclePointsParents();
        muscles.push_back(flexor);
    }
    else if (currentScene == "OFINGER_Absolut"){
        // --- PARAMETER ---
        double handLength = 1.8 * GEOMETRYSCALER;
        float GS = GEOMETRYSCALER;
        std::vector<double> fLength = {0.3482 * handLength, 0.2027 * handLength, 0.1175 * handLength, 0.0882 * handLength};
        double L1 = fLength[0]; // Länge Segment 1
        double L2 = fLength[1]; // Länge Segment 2
        double L3 = fLength[2]; // Länge Segment 3
        double L4 = fLength[3]; // Länge Segment 4
        double rWF = 0.8; // reduce width factor
        std::vector<double> width = {0.15*GS*rWF, 0.1*GS*rWF, 0.07*GS*rWF, 0.05*GS*rWF}; // Dicke
        std::vector<double> relTorusPos = {0.6, 0.42, 0.32};
        std::vector<double> relTorusR = {1.85/rWF, 2.0/rWF, 2.3/rWF}; // {1.85, 2.0, 2.3}
        std::vector<double> relTorusr = {1.0/rWF, 1.0/rWF, 1.2/rWF}; // {1.0, 1.0, 1.2}
        double off = 0.0; 

        // bools
        int bm1 = 1; // Mesh Segment 1
        int bm2 = 1; // Mesh Segment 2
        int bm3 = 1; // Mesh Segment 3
        int bm4 = 1; // Mesh Segment 4

        int bT1 = 1;  // Torus Segment 1
        int bT2 = 1;  // Torus Segment 2
        int bT3 = 1;  // Torus Segment 3

        int bJ1 = 1; // Gelenk 1
        int bJ2 = 1; // Gelenk 2
        int bJ3 = 1; // Gelenk 3
        qDebug() << "Case:" << bm1 << bm2 << bm3 << bm4 << bT1 << bT2 << bT3 << bJ1 << bJ2 << bJ3;

        bool usePreScribed = 0;
        float boneWidthScaler = 1.0;
        

        // 1. ROOT
        rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), nullptr);
        tissues.push_back(rootSystem);

        // ==============================================================================
        // SEGMENT 1 (Proximal) - Referenziert auf Root
        // ==============================================================================
        MWMath::Point3D pBody1(0, 0, 0);

        auto body1 = std::make_shared<SSBody>("Body_Seg1", pBody1, MWMath::RotMatrix3x3(), rootSystem);
        tissues.push_back(body1);
        // Mesh 1 - Relativ zu Body 1
        auto mesh1 = std::make_shared<SSEllipsoidMesh>(width[0], width[0], L1, "Mesh_Seg1", body1, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0.0));
        if (bm1) meshes.push_back(mesh1);
        // Torus 1 - Relativ zu Body 1
        auto mTorus1 = std::make_shared<SSTorusMesh>(mesh1->B * relTorusR[0], mesh1->B * relTorusr[0], 
             "Torus1", body1, MWMath::Point3D(0, -width[0]*off, L1 * relTorusPos[0]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(1, 0, 0));
        if (bT1) meshes.push_back(mTorus1); 

        // ==============================================================================
        // JOINT 1 & 11 - Referenziert auf Root
        // ==============================================================================
        MWMath::Point3D pJoint1(L1, 0, 0);

        auto joint = std::make_shared<SSJoint>("Joint_1", pJoint1, MWMath::RotMatrix3x3(), rootSystem, FJAngles[0], MWMath::Point3D(0,0,-1), numTimeSteps);
        tissues.push_back(joint);
        auto jMesh1 = std::make_shared<SSEllipsoidMesh>(width[0]*1.1, width[0]*1.1, width[0]*1.1, "Mesh_Joint1", joint, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        if (bJ1) meshes.push_back(jMesh1);
        auto joint1 = std::make_shared<SSJoint>("Joint_11", pJoint1, MWMath::axisAngle({0,1,0}, 0.0), rootSystem, FJAngles[1], MWMath::Point3D(0,1,0), numTimeSteps);
        tissues.push_back(joint1);
       

        // ==============================================================================
        // SEGMENT 2 - Referenziert auf Root
        // ==============================================================================
        MWMath::Point3D pBody2 = pJoint1 + MWMath::Point3D(L2, 0, 0); // Pos: L1 + L2
        
        auto body2 = std::make_shared<SSBody>("Body_Seg2", pBody2, MWMath::RotMatrix3x3(), rootSystem);
        tissues.push_back(body2);
        // Mesh 2 - Relativ zu Body 2
        auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2, "Mesh_Seg2", body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.0, 0.8, 0.0));
        
        if (bm2) meshes.push_back(mesh2);
        // Torus 2 - Relativ zu Body 2
        auto mTorus2 = std::make_shared<SSTorusMesh>(mesh2->B * relTorusR[1], mesh2->B * relTorusr[1], 
             "Torus2", body2, MWMath::Point3D(0, -width[1]*off, L2 * relTorusPos[1]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
        if (bT2) meshes.push_back(mTorus2);

        // ==============================================================================
        // JOINT 2 - Referenziert auf Root
        // ==============================================================================
        MWMath::Point3D pJoint2 = pBody2 + MWMath::Point3D(L2, 0, 0); // Pos: L1 + 2*L2
        
        auto joint2 = std::make_shared<SSJoint>("Joint_2", pJoint2, MWMath::RotMatrix3x3(), rootSystem, FJAngles[2], MWMath::Point3D(0,0,-1), numTimeSteps);
        tissues.push_back(joint2);
        auto jMesh2 = std::make_shared<SSEllipsoidMesh>(width[1]*1.1, width[1]*1.1, width[1]*1.1, "Mesh_Joint2", joint2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        if (bJ2) meshes.push_back(jMesh2);

        // ==============================================================================
        // SEGMENT 3 - Referenziert auf Root
        // ==============================================================================
        MWMath::Point3D pBody3 = pJoint2 + MWMath::Point3D(L3, 0, 0); // Pos: L1 + 2*L2 + L3
        auto body3 = std::make_shared<SSBody>("Body_Seg3", pBody3, MWMath::RotMatrix3x3(), rootSystem);
        tissues.push_back(body3);
        
        // Mesh 3 - Relativ zu Body 3
        float mesh3B = width[2];
        auto mesh3 = std::make_shared<SSEllipsoidMesh>(width[2]*boneWidthScaler, width[2]*boneWidthScaler, L3, "Mesh_Seg3", body3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.0, 0.8, 0.8));
        //auto mesh3 = std::make_shared<SSCylinderMesh>(width[2]*boneWidthScaler, L3*2, "Mesh_Seg3", body3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.0, 0.8, 0.8));
        if (bm3) meshes.push_back(mesh3);
        
        // Torus 3 - Relativ zu Body 3
        auto mTorus3 = std::make_shared<SSTorusMesh>(mesh3B * relTorusR[2], mesh3B * relTorusr[2], 
             "Torus3", body3, MWMath::Point3D(0, -width[2]*off, L3 * relTorusPos[2]), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
        if (bT3) meshes.push_back(mTorus3); 

        // ==============================================================================
        // JOINT 3 - Referenziert auf Root
        // ==============================================================================
        MWMath::Point3D pJoint3 = pBody3 + MWMath::Point3D(L3, 0, 0); // Pos: L1 + 2*L2 + 2*L3
        auto joint3 = std::make_shared<SSJoint>("Joint_3", pJoint3, MWMath::RotMatrix3x3(), rootSystem, FJAngles[3], MWMath::Point3D(0,0,-1), numTimeSteps);
        if (usePreScribed) joint3->AngleSteps = jointAnglesPrescribed;
        tissues.push_back(joint3);

        auto jMesh3 = std::make_shared<SSEllipsoidMesh>(width[2]*1.1, width[2]*1.1, width[2]*1.1, "Mesh_Joint3", joint3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        if (bJ3) meshes.push_back(jMesh3);

        // ==============================================================================
        // SEGMENT 4 (Distal) - BEWEGLICH! Referenziert auf Joint 3!
        // ==============================================================================
        // Da das Gelenk bereits an Position pJoint3 sitzt, ist die lokale Position L4.
        auto body4 = std::make_shared<SSBody>("Body_Seg4", MWMath::Point3D(L4, 0, 0), MWMath::RotMatrix3x3(), joint3);
        tissues.push_back(body4);
        
        // Mesh 4 - Relativ zu Body 4
        auto mesh4 = std::make_shared<SSEllipsoidMesh>(width[3], width[3], L4, "Mesh_Seg4", body4, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(0.9, 0.5, 0.3));
        if (bm4) meshes.push_back(mesh4);

        // ==============================================================================
        // INITIALES UPDATE
        // ==============================================================================
        for (auto& m : meshes) {m->InitializeMesh();}
        rootSystem->update(0); 

        // ==============================================================================
        // MUSKEL
        // ==============================================================================
        int numPoints = cfg.muscleNumPoints[0];
        
        MWMath::Point3D startOffset = MWMath::Point3D(0.0, -width[0] * 1.1, 0.0);
        MWMath::Point3D endOffset = MWMath::Point3D(0.0, -mesh4->B * 1.1, 0.0);
        
        SSMuscle* flexor = new SSMuscle("Flexor", numPoints, 
            rootSystem.get(), startOffset, 
            mesh4.get(), endOffset);

        for(auto& m : meshes) {
             flexor->meshPtrs.push_back(m.get());
        }

        flexor->createMusclePoints();
        flexor->updateMusclePointsParents();
        muscles.push_back(flexor);
    }

    // elliptischer Torus tests
    else if (currentScene == "OFINGER_2SegsELLTORUS"){
        // --- PARAMETER ---
        double handLength = 1.8 * GEOMETRYSCALER;
        float GS = GEOMETRYSCALER;
        
        // Segmentlängen
        double L_Prox = 0.1175 * handLength; 
        double L_Dist = 0.0882 * handLength; 
        
        // Knochenbreiten
        double w_Prox = 0.07 * GS;           
        double w_Dist = 0.05 * GS;           
        
        // --- ELLIPTISCHER TORUS PARAMETER ---
        double torusRelPos = -0.2; // 0.32 
        
        // Pfad-Ellipse (Der Ring): 
        // RX ist der Radius in der Hauptebene, RY kann zur Abflachung des Rings genutzt werden
        double RX = w_Prox * 2.3; // entspricht torusRelR
        double RY = w_Prox * 2.1; // Beispiel für leichte Elliptizität des Pfades
        
        // Rohr-Ellipse (Der Querschnitt):
        // R1 ist die Breite des Schlauchs, R2 die Höhe
        double R1 = w_Prox * 1.0; // entspricht torusRelr
        double R2 = w_Prox * 3.0; // Beispiel für abgeflachten Schlauchquerschnitt
        
        double off = 0.0; 

        // 1. ROOT
        rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), nullptr);
        tissues.push_back(rootSystem);

        // ==============================================================================
        // SEGMENT 1 (Proximal) - Fixiert an Root
        // ==============================================================================
        auto bodyProx = std::make_shared<SSBody>("Body_Prox", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), rootSystem);
        tissues.push_back(bodyProx);
        
        // Mesh Proximal
        auto meshProx = std::make_shared<SSEllipsoidMesh>(w_Prox, w_Prox, L_Prox, 
            "Mesh_Prox", bodyProx, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), 
            MWMath::Point3D(0.0, 0.8, 0.8));
        meshes.push_back(meshProx);
        
        // --- ELLIPTISCHER TORUS Proximal ---
        // Nutzt den neuen Konstruktor: R1, R2, RX, RY
        auto mTorusProx = std::make_shared<SSEllipticalTorusMesh>(
            R1, R2, RX, RY,
            "Torus_Prox", 
            bodyProx, 
            MWMath::Point3D(0, -w_Prox * off, L_Prox * torusRelPos), 
            MWMath::axisAngle({0,1,0}, 90.0), 
            MWMath::Point3D(0, 0, 0.8)
        );
        meshes.push_back(mTorusProx);

        // ==============================================================================
        // GELENK (Sitzt am Ende von Segment 1)
        // ==============================================================================
        MWMath::Point3D jointPosRel = MWMath::Point3D(L_Prox, 0, 0);
        auto joint = std::make_shared<SSJoint>("Joint_1", 
                                            jointPosRel, 
                                            MWMath::RotMatrix3x3(), 
                                            bodyProx, 
                                            90.0, 
                                            MWMath::Point3D(0,0,-1),
                                            numTimeSteps);
        tissues.push_back(joint);

        // Visuelle Kugel für das Gelenk
        auto jMesh = std::make_shared<SSEllipsoidMesh>(w_Prox*1.1, w_Prox*1.1, w_Prox*1.1, "Mesh_Joint", joint, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh);

        // ==============================================================================
        // SEGMENT 2 (Distal) - Hängt am Gelenk
        // ==============================================================================
        auto bodyDist = std::make_shared<SSBody>("Body_Dist", MWMath::Point3D(L_Dist,0,0), MWMath::RotMatrix3x3(), joint);
        tissues.push_back(bodyDist);
        
        // Mesh Distal
        auto meshDist = std::make_shared<SSEllipsoidMesh>(w_Dist, w_Dist, L_Dist, 
            "Mesh_Dist", bodyDist, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({0,1,0}, 90.0), 
            MWMath::Point3D(0.9, 0.5, 0.3));
        meshes.push_back(meshDist);

        // ==============================================================================
        // INITIALES UPDATE & DISKRETISIERUNG
        // ==============================================================================
        for (auto& m : meshes) {
            m->InitializeMesh();
            // WICHTIG: Hier wird die neue discretizeMesh Methode für den elliptischen Torus aufgerufen
            m->discretizeMesh(cfg.discretization); 
        }
        rootSystem->update(0); 

        // ==============================================================================
        // MUSKEL
        // ==============================================================================
        int numPoints = cfg.muscleNumPoints.empty() ? 25 : cfg.muscleNumPoints[0];
        
        MWMath::Point3D startOffset = MWMath::Point3D(-handLength*0.7, -0.15*1.1, 0.0);
        MWMath::Point3D endOffset = MWMath::Point3D(meshDist->C, -meshDist->B * 1.1, 0.0);
        
        SSMuscle* flexor = new SSMuscle("Flexor", numPoints, 
            bodyProx.get(), startOffset, 
            bodyDist.get(), endOffset);

        for(auto& m : meshes) {
            flexor->meshPtrs.push_back(m.get());
        }

        flexor->createMusclePoints();
        flexor->updateMusclePointsParents();
        muscles.push_back(flexor);
    }
    else if (currentScene == "OFINGER_3SegsELLTORUS") {
        // --- PARAMETER ---
        double handLength = 1.8 * GEOMETRYSCALER;
        float GS = GEOMETRYSCALER;

        // Längen der Fingerglieder
        std::vector<double> fLength = { 0.3482 * handLength, 0.2027 * handLength, 0.1175 * handLength, 0.0882 * handLength };
        double L2 = fLength[1]; // Proximal
        double L3 = fLength[2]; // Mitte
        double L4 = fLength[3]; // Distal

        std::vector<double> width = { 0.15 * GS, 0.1 * GS, 0.07 * GS, 0.05 * GS };
        
        // Relative Radien für die Pfad-Ellipse (RX, RY) und die Rohr-Breite (R1)
        std::vector<double> relTorusRX = { 1.85, 2.0, 2.3 }; 
        std::vector<double> relTorusRY = { 2.0, 2.4, 2.6 }; // RY etwas kleiner für flachere Ringbänder
        std::vector<double> relTorusR1 = { 0.7, 0.7, 0.9 }; 
        
        double off = 0.0; 

        // 1. ROOT
        rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), nullptr);
        tissues.push_back(rootSystem);

        // ==============================================================================
        // SEGMENT 2 (Proximal)
        // ==============================================================================
        auto body2 = std::make_shared<SSBody>("Body_Seg2", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), rootSystem);
        tissues.push_back(body2);

        auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2, 
            "Mesh_Seg2", body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 0, 1, 0 }, 90.0), MWMath::Point3D(0.0, 0.8, 0.0));
        meshes.push_back(mesh2);

        // Elliptischer Torus für Segment 2
        // R1: Breite, R2: 0.6 * Mesh-Breite, RX/RY: Pfad
        auto mTorus2 = std::make_shared<SSEllipticalTorusMesh>(
            width[1] * relTorusR1[0], // R1
            L2 * 0.7,           // R2 (fest vorgegeben)
            width[1] * relTorusRX[0], // RX
            width[1] * relTorusRY[0], // RY
            "Torus_Seg2", body2, 
            MWMath::Point3D(0, 0, 0), // Relativ mittig auf dem Mesh
            MWMath::axisAngle({ 0, 1, 0 }, 90.0), 
            MWMath::Point3D(0.8, 0.8, 0)
        );
        meshes.push_back(mTorus2);

        // ==============================================================================
        // GELENK 2 ( MCP / PIP Übergang )
        // ==============================================================================
        MWMath::Point3D jointPosRel2 = MWMath::Point3D(L2, 0, 0);
        auto joint2 = std::make_shared<SSJoint>("Joint_2", jointPosRel2, MWMath::RotMatrix3x3(), body2, 90.0, MWMath::Point3D(0, 0, -1), numTimeSteps);
        tissues.push_back(joint2);

        // ==============================================================================
        // SEGMENT 3 (Mitte)
        // ==============================================================================
        auto body3 = std::make_shared<SSBody>("Body_Seg3", MWMath::Point3D(L3, 0, 0), MWMath::RotMatrix3x3(), joint2);
        tissues.push_back(body3);

        auto mesh3 = std::make_shared<SSEllipsoidMesh>(width[2], width[2], L3, 
            "Mesh_Seg3", body3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 0, 1, 0 }, 90.0), MWMath::Point3D(0.0, 0.8, 0.8));
        meshes.push_back(mesh3);

        // Elliptischer Torus für Segment 3
        // Positioniert in der Mitte des Segments (Point3D(0,0,0) da Body im Zentrum)
        auto mTorus3 = std::make_shared<SSEllipticalTorusMesh>(
            width[2] * relTorusR1[1], // R1
            L3 * 0.7,           // R2 (fest vorgegeben)
            width[2] * relTorusRX[1], // RX
            width[2] * relTorusRY[1], // RY
            "Torus_Seg3", body3, 
            MWMath::Point3D(0, 0, 0), // Zentriert auf L3
            MWMath::axisAngle({ 0, 1, 0 }, 90.0), 
            MWMath::Point3D(0, 0, 0.8)
        );
        meshes.push_back(mTorus3);

        // ==============================================================================
        // GELENK 3
        // ==============================================================================
        MWMath::Point3D jointPosRel3 = MWMath::Point3D(L3, 0, 0);
        auto joint3 = std::make_shared<SSJoint>("Joint_3", jointPosRel3, MWMath::RotMatrix3x3(), body3, 90.0, MWMath::Point3D(0, 0, -1), numTimeSteps);
        tissues.push_back(joint3);

        // ==============================================================================
        // SEGMENT 4 (Distal)
        // ==============================================================================
        auto body4 = std::make_shared<SSBody>("Body_Seg4", MWMath::Point3D(L4, 0, 0), MWMath::RotMatrix3x3(), joint3);
        tissues.push_back(body4);

        auto mesh4 = std::make_shared<SSEllipsoidMesh>(width[3], width[3], L4, 
            "Mesh_Seg4", body4, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 0, 1, 0 }, 90.0), MWMath::Point3D(0.9, 0.5, 0.3));
        meshes.push_back(mesh4);

        // ==============================================================================
        // INITIALES UPDATE & DISKRETISIERUNG
        // ==============================================================================
        for (auto& m : meshes) { 
            m->InitializeMesh(); 
            // Wichtig: Ruft die neue elliptische Diskretisierung auf
            m->discretizeMesh(cfg.discretization); 
        }
        rootSystem->update(0);

        // ==============================================================================
        // MUSKEL
        // ==============================================================================
        int numPoints = cfg.muscleNumPoints.empty() ? 25 : cfg.muscleNumPoints[0];

        MWMath::Point3D startOffset = MWMath::Point3D(-handLength * 0.7, -width[1] * 1.1, 0.0);
        MWMath::Point3D endOffset = MWMath::Point3D(mesh4->C, -width[3] * 1.1, 0.0);

        SSMuscle* flexor = new SSMuscle("Flexor", numPoints, body2.get(), startOffset, body4.get(), endOffset);

        for (auto& m : meshes) {
            flexor->meshPtrs.push_back(m.get());
        }

        flexor->createMusclePoints();
        flexor->updateMusclePointsParents();
        muscles.push_back(flexor);
    }

    else if (currentScene == "OHandTest"){
        // --- PARAMETER ---
        double HL = 1.8 * GEOMETRYSCALER;
        double HB = 0.8 * GEOMETRYSCALER;
        double scale = HL / Hand::RefHandLength;
        double Sign = 1.0; // 1.0 = Rechts
        float GS = GEOMETRYSCALER;
        std::unordered_map<std::string, SSMesh*> MeshMap; // Hilfreich für spätere Muskeldefinitionen

        std::vector<double> width = { 0.15*0.5 * GS, 0.1*0.5 * GS, 0.07*0.5 * GS, 0.05*0.5 * GS };
        
        // --- TORUS PARAMETER AUS DEINER VORLAGE ---
        double rWF = 1.0; 
        std::vector<double> relTorusPos = {0.28, 0.16, 0.06}; // {0.32, 0.2, 0.1};         // Position auf dem Knochen
        std::vector<double> relTorusR = {0.14, 0.19, 0.28}; //{1.6/rWF, 1.8/rWF, 1.3/rWF}; // Großer Radius
        std::vector<double> relTorusr = {0.13, 0.18, 0.26};// {1.4/rWF, 1.6/rWF, 1.0/rWF}; // Kleiner Radius (Dicke)
        std::vector<double> off = {0.06, 0.07, 0.18}; // Palmarer Offset (Verschiebung Richtung Handinnenfläche)

        double segLen[4];
        double segRadius[4];
        double jointRadius[3];
        std::vector<std::string> fingerNames = {"Thumb", "Index", "Middle", "Ring", "Little"};

        // 1. ROOT
        rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), nullptr);
        tissues.push_back(rootSystem);

        // Drehpunkt des Handgelenks
        MWMath::Point3D wristOffset(0.0, 0.003567 * scale, -0.003901 * scale * Sign);
        
        // Z-Achse (0,0,1) als Rotationsachse ist KORREKT für Handgelenks-Flexion/Extension!
        auto wristJointAbd = std::make_shared<SSJoint>("Wrist_JointAbd", wristOffset, MWMath::RotMatrix3x3(), rootSystem, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(wristJointAbd);
        auto jMeshWrist = std::make_shared<SSEllipsoidMesh>(0.01 * scale, 0.01 * scale, 0.01 * scale, "Mesh_Wrist_Joint", wristJointAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0,0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMeshWrist);
        MeshMap["Mesh_Wrist_Joint"] = jMeshWrist.get();

        auto wristJointFlex = std::make_shared<SSJoint>("Wrist_JointFlex", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointAbd, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(wristJointFlex);
        // ==========================================
        // 2. CARPALS (Die Handwurzel-Basisplatte)
        // ==========================================
        auto carpals = std::make_shared<SSBody>("Carpals", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointFlex);
        tissues.push_back(carpals);
        
        // KORREKTUR: 
        // 1. Offset auf +Y (z.B. +0.03), da die Handwurzel *vor* dem Handgelenk liegt.
        // 2. Wir definieren die Ellipsoid-Maße direkt passend zu den Achsen:
        //    X (Dicke), Y (Länge), Z (Breite). Keine zusätzliche Rotation nötig!
        double carpalThickness = width[0] * 1.5; // X
        double carpalLength = width[0] * 2.5;    // Y
        double carpalWidth = width[0] * 3.5;     // Z
        
        auto meshCarpals = std::make_shared<SSEllipsoidMesh>(carpalThickness, carpalLength, carpalWidth, "Mesh_Carpals", carpals, MWMath::Point3D(0, 0.0*scale, 0), MWMath::RotMatrix3x3(), MWMath::Point3D(0.5, 0.5, 0.5));
        meshes.push_back(meshCarpals);
        MeshMap["Mesh_Carpals"] = meshCarpals.get();

        // ==========================================
        // 3. WRAPPING SURFACE (Karpaltunnel-Zylinder)
        // ==========================================
        double flexCylRadius = 0.0115 * scale;
        double flexCylLength = 0.1400 * scale; // Notiz: 14cm ist sehr lang, in AnyBody oft so gemacht, damit Sehnen nicht abrutschen.
        
        // KORREKTUR:
        // Zylinder muss auf die Palmar-Seite (+X) geschoben werden und leicht nach vorne (+Y).
        MWMath::Point3D flexCylOffset(carpalThickness * 0.5 + flexCylRadius, 0.02 * scale, 0.0);

        // KORREKTUR ROTATION:
        // Der Zylinder muss entlang der Z-Achse liegen! 
        // Wenn dein SSCylinderMesh standardmäßig entlang der Y-Achse (Länge) generiert wird, 
        // kippt ihn eine Rotation um 90° um die X-Achse {1, 0, 0} exakt auf die Z-Achse.
        auto flexorCylMesh = std::make_shared<SSCylinderMesh>(flexCylRadius, flexCylLength, "Mesh_WristFlexorCyl", carpals, flexCylOffset, MWMath::axisAngle({ 1, 0, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));  
        meshes.push_back(flexorCylMesh);
        MeshMap["Mesh_WristFlexorCyl"] = flexorCylMesh.get();

        

        // INDEX // 
        for (int i = 1; i < 2; i++) {
            int fidx = i;
            std::string cFName = fingerNames[fidx];
            double jAngle;
            std::string prefN = std::to_string(fidx);

            // MCP
            segLen[0] = Hand::SegLenRatios[fidx][0] * HL;
            segLen[1] = Hand::SegLenRatios[fidx][1] * HL;
            segLen[2] = Hand::SegLenRatios[fidx][2] * HL;
            segLen[3] = Hand::SegLenRatios[fidx][3] * HL;

            segRadius[0] = Hand::FingerWrapData[fidx-1][0].radius * scale;
            segRadius[1] = Hand::FingerWrapData[fidx-1][1].radius * scale;
            segRadius[2] = Hand::FingerWrapData[fidx-1][2].radius * scale;
            segRadius[3] = Hand::FingerWrapData[fidx-1][3].radius * scale;

            jointRadius[0] = segRadius[0] * 1.1; /* Hand::JointCylRadiiRel[fidx].mcp * HB; */
            jointRadius[1] = segRadius[2] * 1.1; /* Hand::JointCylRadiiRel[fidx].pip * HB; */
            jointRadius[2] = segRadius[3] * 1.1; /* Hand::JointCylRadiiRel[fidx].dip * HB; */

            MWMath::Point3D posMCP(0.0, Hand::JointPosTable[fidx][0] * HL, Hand::JointPosTable[fidx][1] * HB * Sign);
            

            // CMC Joint + Mesh
            MWMath::Point3D posCMC = Hand::CMCOffsets[fidx - 1] * scale;
            posCMC.z *= Sign;
            MWMath::RotMatrix3x3 forientation = buildOrientation(posMCP, posCMC);
            jAngle = 0.0; // CMC_AA_M[fidx];
            std::string jName = "CMC" + prefN + "_JointAbd";
            auto jointCMCAbd = std::make_shared<SSJoint>(jName, posCMC, forientation, carpals, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
            tissues.push_back(jointCMCAbd);
            auto jMeshCMC = std::make_shared<SSEllipsoidMesh>(width[0], width[0], width[0], "Mesh_CMC"+prefN+"_Joint", jointCMCAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0 , 0 }, 90.0), MWMath::Point3D(0.8, 0.8, 0.0));
            meshes.push_back(jMeshCMC);
            MeshMap["Mesh_CMC"+prefN+"_Joint"] = jMeshCMC.get();

            jAngle = 0.0; // CMC_FE_M[fidx];
            jName = "CMC" + prefN + "_JointFlex";
            auto jointCMCFlex = std::make_shared<SSJoint>(jName, MWMath::Point3D(0.,0.,0.), MWMath::RotMatrix3x3(), jointCMCAbd, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
            tissues.push_back(jointCMCFlex);

            // ==========================================
            // MC SEGMENT
            // ==========================================
            double flength = MWMath::distance(posCMC, posMCP);
            // Body MC + Mesh
            jName = "MC" + prefN;
            auto bodyMC = std::make_shared<SSBody>(jName, MWMath::Point3D(0.,-flength*0.5,0.), MWMath::RotMatrix3x3(), jointCMCFlex);
            tissues.push_back(bodyMC);
            auto meshMC = std::make_shared<SSEllipsoidMesh>(segRadius[0], segRadius[0], flength*0.5, "Mesh_"+jName, bodyMC, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1, 0, 0 }, 90.0), COLORF1);
            meshes.push_back(meshMC);
            MeshMap["Mesh_" + jName] = meshMC.get();

            // --- NEU: TORUS MC ---
            double tR0 = segLen[0] * relTorusR[0];
            double tr0 = segLen[0] * relTorusr[0];
            double tPosY0 = -flength * relTorusPos[0]; // Position von Gelenk aus gerechnet
            auto mTorusMC = std::make_shared<SSTorusMesh>(tR0, tr0, "MeshTorus_"+jName, bodyMC, 
                 MWMath::Point3D(segLen[0] * off[0], tPosY0, 0), // +X ist Palmar!
                 MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(1, 0, 0));
            meshes.push_back(mTorusMC);
            MeshMap["MeshTorus_" + jName] = mTorusMC.get();

             
            // MCP Joint + Mesh
            MWMath::Point3D jointPosRelMCP = MWMath::Point3D(0, -flength*0.5, 0); 
            jAngle = 0.0; // MMHA.MCP_AA_M[fidx];
            jName = "MCP" + prefN + "_JointAbd";
            auto jointMCPAbd = std::make_shared<SSJoint>(jName, jointPosRelMCP, MWMath::RotMatrix3x3(), bodyMC, jAngle, MWMath::Point3D(1, 0, 0), numTimeSteps);
            tissues.push_back(jointMCPAbd);
            auto jMeshMCP = std::make_shared<SSEllipsoidMesh>(jointRadius[0], jointRadius[0], jointRadius[0], "Mesh_MCP"+prefN+"_Joint", jointMCPAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0 , 0 }, 90.0), MWMath::Point3D(0.8, 0.8, 0.0));
            meshes.push_back(jMeshMCP);
            MeshMap["Mesh_MCP"+prefN+"_Joint"] = jMeshMCP.get();

            jAngle = MMHA.MCP_FE_M[fidx];
            jName = "MCP" + prefN + "_JointFlex";
            auto jointMCPFlex = std::make_shared<SSJoint>(jName, MWMath::Point3D(0.,0.,0.), MWMath::RotMatrix3x3(), jointMCPAbd, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
            tissues.push_back(jointMCPFlex);
            

            // ==========================================
            // PP SEGMENT
            // ==========================================
            MWMath::Point3D posPP = MWMath::Point3D(0, -segLen[1]*0.5, 0);
            jName =  "PP" + prefN;
            auto bodyPP = std::make_shared<SSBody>(jName, posPP, MWMath::RotMatrix3x3(), jointMCPFlex);
            tissues.push_back(bodyPP);
            auto meshPP = std::make_shared<SSEllipsoidMesh>(segRadius[1], segRadius[1], segLen[1]*0.5, "Mesh_"+jName, bodyPP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1, 0, 0 }, 90.0), COLORF2);
            meshes.push_back(meshPP);
            MeshMap["Mesh_" + jName] = meshPP.get();

            // --- NEU: TORUS PP ---
            double tR1 = segLen[1] * relTorusR[1];
            double tr1 = segLen[1] * relTorusr[1];
            double tPosY1 = -segLen[1] *relTorusPos[1]; 
            auto mTorusPP = std::make_shared<SSTorusMesh>(tR1, tr1, "MeshTorus_"+jName, bodyPP, 
                 MWMath::Point3D(segLen[1] * off[1], tPosY1, 0), 
                 MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0.8, 0.3, 0));
            meshes.push_back(mTorusPP);
            MeshMap["MeshTorus_" + jName] = mTorusPP.get();

            
            // PIP Joint + Mesh
            MWMath::Point3D jointPosRelPP = MWMath::Point3D(0, -segLen[1]*0.5, 0); 
            jAngle = MMHA.PIP_M[fidx];
            jName = "PIP" + prefN + "_Joint";
            auto jointPIP = std::make_shared<SSJoint>(jName, jointPosRelPP, MWMath::RotMatrix3x3(), bodyPP, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
            tissues.push_back(jointPIP);
            //auto MeshPIP = std::make_shared<SSEllipsoidMesh>(jointRadius[1], jointRadius[1], jointRadius[1], cFName + "_Mesh_PIP_Joint", jointPIP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 0, 1, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
            auto MeshPIP = std::make_shared<SSCylinderMesh>(jointRadius[1], jointRadius[1]*3, "Mesh_"+jName, jointPIP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 0, 1, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
            meshes.push_back(MeshPIP);
            MeshMap["Mesh_" + jName] = MeshPIP.get();

            // ==========================================
            // MP SEGMENT
            // ==========================================
            MWMath::Point3D posMP = MWMath::Point3D(0., -segLen[2]*0.5, 0);
            jName =  "MP" + prefN;
            auto bodyMP = std::make_shared<SSBody>(jName, posMP, MWMath::RotMatrix3x3(), jointPIP);
            tissues.push_back(bodyMP);
            auto meshMP = std::make_shared<SSEllipsoidMesh>(segRadius[2], segRadius[2], segLen[2]*0.5, "Mesh_"+jName, bodyMP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0, 0 }, 90.0), COLORF3);
            meshes.push_back(meshMP);
            MeshMap["Mesh_" + jName] = meshMP.get();

            // --- NEU: TORUS MP ---
            double tR2 = segLen[2] * relTorusR[2];
            double tr2 = segLen[2] * relTorusr[2];
            double tPosY2 = -segLen[2] * relTorusPos[2]; 
            auto mTorusMP = std::make_shared<SSTorusMesh>(tR2, tr2, "MeshTorus_"+jName, bodyMP, 
                 MWMath::Point3D(segLen[2] * off[2], tPosY2, 0), 
                 MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0.8, 0.3, 0.2));
            meshes.push_back(mTorusMP);
            MeshMap["MeshTorus_" + jName] = mTorusMP.get();

            // DIP Joint + Mesh
            MWMath::Point3D jointPosRelMP = MWMath::Point3D(0,-segLen[2]*0.5, 0); 
            jAngle = MMHA.DIP_M[fidx];
            jName = "DIP" + prefN + "_Joint";
            auto jointDIP = std::make_shared<SSJoint>(jName, jointPosRelMP, MWMath::RotMatrix3x3(), bodyMP, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
            tissues.push_back(jointDIP);
            //auto MeshDIP = std::make_shared<SSEllipsoidMesh>(jointRadius[2], jointRadius[2], jointRadius[2], cFName + "_Mesh_DIP_Joint", jointDIP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 0, 1, 0 }, 90.0), MWMath::Point3D(0.8, 0.8, 0.0));
            auto MeshDIP = std::make_shared<SSCylinderMesh>(jointRadius[2], jointRadius[2]*3, "Mesh_"+jName, jointDIP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 0, 1, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
            meshes.push_back(MeshDIP);
            MeshMap["Mesh_" + jName] = MeshDIP.get();

            // ==========================================
            // DP SEGMENT
            // ==========================================
            MWMath::Point3D posDP = MWMath::Point3D(0., -segLen[3]*0.5, 0.);
            jName =  "DP" + prefN;
            auto bodyDP = std::make_shared<SSBody>(jName, posDP, MWMath::RotMatrix3x3(), jointDIP);
            tissues.push_back(bodyDP);
            auto meshDP = std::make_shared<SSEllipsoidMesh>(segRadius[3], segRadius[3], segLen[3]*0.5, "Mesh_"+jName, bodyDP, MWMath::Point3D(0., 0., 0.), MWMath::axisAngle({ 1,0,0}, 90.0), COLORF4);
            meshes.push_back(meshDP);
            MeshMap["Mesh_" + jName] = meshDP.get(); 
           
        }

        // print mesh map
        for (auto i : MeshMap){
            qDebug() << "MeshMap: " << QString::fromStdString(i.first) << " -> " << i.second->Name.c_str();
        }


        // ==============================================================================
        // INITIALES UPDATE & DISKRETISIERUNG
        // ==============================================================================
        for (auto& m : meshes) { 
            m->InitializeMesh(); 
            m->discretizeMesh(cfg.discretization); 
        }
        rootSystem->update(0);

        // ==============================================================================
        // MUSKEL
        // ==============================================================================
        int numPoints = cfg.muscleNumPoints.empty() ? 25 : cfg.muscleNumPoints[0];
        // "Mesh_WristFlexorCyl", "Mesh_CMC1_Joint"
        std::vector<std::string> refBodyNames = { "Mesh_WristFlexorCyl", "MeshTorus_MC1", "Mesh_MCP1_Joint", "MeshTorus_PP1", "Mesh_PIP1_Joint", "MeshTorus_MP1", "Mesh_DIP1_Joint"};
        MWMath::Point3D startOffset = MWMath::Point3D(0., 0.5, 0.0);
        MWMath::Point3D endOffset = MWMath::Point3D(0.05, 0.0, 0.0);
        SSMuscle* flexor = new SSMuscle("Flexor", numPoints, rootSystem.get(), startOffset, getRefMesh("Mesh_DP1", MeshMap), endOffset);

        for (auto& m : getRefMeshes(refBodyNames, MeshMap)) {
            flexor->meshPtrs.push_back(m);
            //if (dynamic_cast<SSTorusMesh*>(m.get())) flexor->meshPtrs.push_back(m.get());
        }

        //flexor->createMusclePoints();
        flexor->createMusclePointsComplexPath();
        flexor->updateMusclePointsParentsLocal();
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
        systems.push_back(new CasadiSystem(singleMuscleList, objFunc, cfg.solverMethod, cfg.casadiParametrization, cfg.bUseManualJacobian, cfg.bSumPhiEta, cfg.bUseWarmstartEtas));
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
        if (currentScene[0] == 'O' || currentScene == "NONE" || currentScene == "FINGER_SIMPLE2" || currentScene == "FINGER_SIMPLE4" || currentScene == "FINGER_SIMPLE4T"){ 
            if (rootSystem) {
                // Ein einziger Aufruf aktualisiert den gesamten Baum!
                /* rootSystem->Position2ParentRelInParentFrame += MWMath::Point3D(0.2,0.2,0.2);
                rootSystem->Orientation2ParentRel *= MWMath::axisAngle({1,1,0}, 10.0); */
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


        for (auto* sys : systems) {
            qDebug() << "---------------------------" << " Solve CasadiSystem: " << QString::fromStdString(sys->CasadiSystemName) << " ---------------------------";
            sys->solveStepX(); // solves all systems
        }
        qDebug() << "-------------------------------------------------------------------------------------------------";
        
        if (cfg.dynamicReparametrization || t < 1)
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
                qDebug() << "    " << QString::fromStdString(mus->Name) << " Length: " << len*0.01 << " m";

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
        stepText.push_back(std::string("useWarmstartEtas: ") + (cfg.bUseWarmstartEtas ? "Yes" : "No"));
        std::string numNodes;
        for (auto* mus : musclePtrs) {
            numNodes += mus->Name + ": " + std::to_string(mus->MNodes.size()) + "  ";
        }
        stepText.push_back("Muscle Nodes: " + numNodes);
        stepText.push_back("Sucess: " + systems[0]->SolverConvergenceMessages[t] + "(" + std::to_string(systems[0]->SolverConvergenceSteps[t]) + ")");
        std::string jointsangle;
        int tooManyJoints = 0;
        for (const auto& join : tissue) {
            if (auto jointPtr = std::dynamic_pointer_cast<SSJoint>(join)) {
                std::stringstream ss;
                if (t < jointPtr->DoneAngleSteps.size()) {
                    if (jointPtr->DoneAngleSteps[t+1] != 0.0)
                    {    ss << std::fixed << std::setprecision(2) << jointPtr->DoneAngleSteps[t+1];
                        jointsangle += jointPtr->Name + ": " + ss.str() + "°  ";
                        tooManyJoints++;
                    }
                }
            }
        }
        if (tooManyJoints < 8) stepText.push_back("Joint Angles: " + jointsangle);
        viewerText.push_back(stepText);
    }
    VTKSimViewerSimple viewer(allMuscleResults, allMeshResults, allInitialGuesses, allOptimizedPointColors, tissue, meshes, musclePtrs, angles, otherPointsWithColors, viewerText, 1.0);
    viewer.show();

    // ##############################################
    // EXPORTS 
    // ##############################################

    std::string outputName = systems[0]->CasadiSystemName + "InputParameters" + currentScene + ".txt";
    exportParameterLog(systems[0]->allParameterInputsAllSteps, systems[0]->allParameterInputDescriptionsAllSteps, outputName);
    
    for (auto* sys : systems) {
        for (auto* mus : sys->m_muscles) {
            exportMuscleLog(sys->CasadiSystemName, mus, sys->SolverConvergenceMessages, UNITS);
        }
    }
    backupSourceCode();

    std::vector<SSJoint*> joints;
    for (const auto& tis : tissue) {
        if (auto jointPtr = std::dynamic_pointer_cast<SSJoint>(tis)) {
            joints.push_back(jointPtr.get());
        }
    }
    exportMuscleLengthLog(musclePtrs, joints);

    return app.exec();

}