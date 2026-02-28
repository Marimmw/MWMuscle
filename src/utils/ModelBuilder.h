#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <cmath>

#include "utils/MWMath.h"

#include "simpleSimulation/SSMuscle.h"
#include "simpleSimulation/SSMeshes.h"
#include "simpleSimulation/SSBody.h"
#include "../utils/ConfigLoader.h" // Für SimSettings

// THIS FILE CONTAINS HELPER FUNCTIONS FOR CREATING MODELS

struct JoingAngleDef
{
    int steps;
    double finalAngleDeg;
};

struct FingerSegParams
{
    double length;
    double width;
    double jointRadius;
};





inline void createJointMovementVector(std::vector<double>& jointAngles, const std::vector<JoingAngleDef>& jointAngleDefs = {{10, 90.0}}) {
    double currentAngle;
    if (jointAngles.empty()){
        currentAngle = 0.0;
    }
    else{
        currentAngle = jointAngles.back();
    }
    for (const JoingAngleDef& def : jointAngleDefs) {
        double stepSize = (def.finalAngleDeg - currentAngle) / def.steps;
        for (int i = 0; i < def.steps; ++i) {
            currentAngle += stepSize;
            jointAngles.push_back(currentAngle);
        }
    }
}

inline MWMath::RotMatrix3x3 buildOrientation(const MWMath::Point3D& start, const MWMath::Point3D& end) {
        
    // 1. Die Y-Achse ist genau die Richtung des Knochens (Endpunkt minus Startpunkt)
    MWMath::Point3D yAxis = (end - start).normed();

    // 2. Die Handfläche liegt grob auf der globalen X-Achse
    MWMath::Point3D palmNormal(1.0, 0.0, 0.0);

    // 3. Z-Achse = Flexionsachse = Quer zum Finger und zur Handfläche
    // Nutzung DEINER cross() Funktion:
    MWMath::Point3D zAxis = palmNormal.cross(yAxis).normed();

    // 4. Echte lokale X-Achse korrigieren (damit alle 3 Achsen 100% senkrecht aufeinander stehen)
    MWMath::Point3D xAxis = yAxis.cross(zAxis).normed();

    // Matrix zusammenbauen (Die Vektoren bilden die Spalten der RotMatrix3x3)
    return MWMath::RotMatrix3x3(
        xAxis.x, yAxis.x, zAxis.x,
        xAxis.y, yAxis.y, zAxis.y,
        xAxis.z, yAxis.z, zAxis.z
    );
}

inline std::vector<SSMesh*> getRefMeshes(std::vector<std::string> meshNames, std::unordered_map<std::string, SSMesh*>& MeshMap) {
    std::vector<SSMesh*> meshPtrs;
    for (const std::string& name : meshNames) {
        if (MeshMap.find(name) != MeshMap.end()) {
            meshPtrs.push_back(MeshMap[name]);
        }
        else {
            std::cerr << "Warnung: Mesh mit Namen '" << name << "' nicht gefunden!" << std::endl;
        }
    }
    return meshPtrs;
}

inline SSMesh* getRefMesh(std::string meshName, std::unordered_map<std::string, SSMesh*>& MeshMap) {
    if (MeshMap.find(meshName) != MeshMap.end()) {
        return MeshMap[meshName];
    }
    else {
        std::cerr << "Warnung: Mesh mit Namen '" << meshName << "' nicht gefunden!" << std::endl;
        return nullptr;
    }
}
///// HAND MODEL BUILDING /////

namespace Hand {

    // Spalten-Index-Legende:
    //  0: CMC F/E Min (Ext)  |  1: CMC F/E Max (Flex)
    //  2: CMC A/A Min (Add)  |  3: CMC A/A Max (Abd)
    //  4: MCP F/E Min (Ext)  |  5: MCP F/E Max (Flex)
    //  6: MCP A/A Min (Add)  |  7: MCP A/A Max (Abd)
    //  8: PIP F/E Min (Ext)  |  9: PIP F/E Max (Flex)   (Reines Scharnier)
    // 10: DIP F/E Min (Ext)  | 11: DIP F/E Max (Flex)   (Reines Scharnier)
    
    constexpr double MinMaxAngles[6][12] = {
        // CMC_FE(0,1)    CMC_AA(2,3)    MCP_FE(4,5)    MCP_AA(6,7)    PIP_FE(8,9)   DIP_FE(10,11)
        {   0.0,   0.0,    0.0,   0.0,    0.0,   0.0,    0.0,   0.0,    0.0,   0.0,    0.0,   0.0}, // 0: Wrist/Dummy
        { -20.0,  45.0,  -20.0,  45.0,  -10.0,  60.0,  -10.0,  10.0,  -15.0,  90.0,    0.0,   0.0}, // 1: Thumb (DIP = 0, da nur ein IP-Gelenk)
        {   0.0,   5.0,    0.0,   0.0,  -30.0,  90.0,  -20.0,  20.0,    0.0, 110.0,    0.0,  90.0}, // 2: Index
        {   0.0,   5.0,    0.0,   0.0,  -30.0,  90.0,  -15.0,  15.0,    0.0, 110.0,    0.0,  90.0}, // 3: Middle
        {  -5.0,  10.0,    0.0,   0.0,  -30.0,  90.0,  -15.0,  15.0,    0.0, 110.0,    0.0,  90.0}, // 4: Ring
        { -10.0,  15.0,    0.0,   0.0,  -30.0,  90.0,  -25.0,  25.0,    0.0, 110.0,    0.0,  90.0}  // 5: Little
    };

    struct MMHandAngles {
        // --- Wrist ---
        double W_FE_m = -90; // Flexion/Extension Min (Extension)
        double W_FE_M = 90; // Flexion/Extension Max (Flexion)
        double W_AA_m = 45; // Abduction/Adduction Min (Adduction)
        double W_AA_M = -25; // Abduction/Adduction Max (Abduction)

        // --- CMC (Carpometacarpal) ---
        //                     [0]     [1]     [2]    [3]    [4]    [5]
        double CMC_FE_m[6] = { 0.0,  -20.0,   0.0,   0.0,  -5.0, -10.0 }; // Flexion/Extension Min
        double CMC_FE_M[6] = { 0.0,   45.0,   5.0,   5.0,  10.0,  15.0 }; // Flexion/Extension Max
        double CMC_AA_m[6] = { 0.0,  -20.0,   0.0,   0.0,   0.0,   0.0 }; // Abduction/Adduction Min
        double CMC_AA_M[6] = { 0.0,   45.0,   0.0,   0.0,   0.0,   0.0 }; // Abduction/Adduction Max

        // --- MCP (Metacarpophalangeal / Grundgelenk) ---
        //                     [0]     [1]     [2]    [3]    [4]    [5]
        double MCP_FE_m[6] = { 0.0,  -10.0, -30.0, -30.0, -30.0, -30.0 }; // Flex/Ext Min
        double MCP_FE_M[6] = { 0.0,   90.0,  90.0,  90.0,  90.0,  90.0 }; // Flex/Ext Max
        double MCP_AA_m[6] = { 0.0,  -10.0, -20.0, -15.0, -15.0, -25.0 }; // Abd/Add Min
        double MCP_AA_M[6] = { 0.0,   10.0,  20.0,  15.0,  15.0,  25.0 }; // Abd/Add Max

        // --- PIP (Proximal Interphalangeal / Mittelgelenk) ---
        //                     [0]     [1]     [2]    [3]    [4]    [5]
        double PIP_m[6] = { 0.0,  -15.0,   0.0,   0.0,   0.0,   0.0 }; // Reines Scharnier
        double PIP_M[6] = { 0.0,   100.0, 100.0, 100.0, 100.0, 100.0 }; 

        // --- DIP (Distal Interphalangeal / Endgelenk) ---
        //                     [0]     [1]     [2]    [3]    [4]    [5]
        double DIP_m[6] = { 0.0,    0.0,   0.0,   0.0,   0.0,   0.0 }; // Reines Scharnier
        double DIP_M[6] = { 0.0,   90.0,  90.0,  90.0,  90.0,  90.0 };
    };

    // Konstanten aus AnyScript
    constexpr double SegLenRatios[5][4] = {
        {0.118, 0.251, 0.196, 0.158}, 
        {0.463, 0.245, 0.143, 0.097}, 
        {0.446, 0.266, 0.170, 0.108}, 
        {0.421, 0.244, 0.165, 0.107}, 
        {0.414, 0.204, 0.117, 0.093}  
    };

    constexpr double JointPosTable[5][2] = {
        {-0.073,  0.196}, 
        {-0.447,  0.251}, 
        {-0.446,  0.000}, 
        {-0.409, -0.206}, 
        {-0.368, -0.402}  
    };

    const MWMath::Point3D CMCOffsets[4] = {
        {-0.004, -0.016,  0.006}, 
        {-0.003, -0.020, -0.008}, 
        { 0.000, -0.020, -0.020}, 
        { 0.004, -0.018, -0.025}  
    };

    constexpr double RefHandLength = 0.195;
    constexpr double ThumbRotX = -60.0; // Grad
    constexpr double ThumbRotY =  44.1; // Grad

    struct BoneWrapParams {
        double radius;
        double lengthCoef;
    };
    constexpr BoneWrapParams FingerWrapData[4][4] = {
        // {radius[mm], lengthRelativeToSegment}
        // Zeigefinger (F2)
        { {0.0070, 0.60}, {0.0045, 0.56}, {0.0060, 0.50}, { 0.0050, 0.50 }},
        // Mittelfinger (F3)
        { {0.0060, 0.50}, {0.0060, 0.50}, {0.0050, 0.50}, { 0.0050, 0.50 }},
        // Ringfinger (F4)
        { {0.0060, 0.50}, {0.0060, 0.50}, {0.0050, 0.50}, { 0.0050, 0.50 }},
        // Kleiner Finger (F5)
        { {0.0085, 0.70}, {0.0060, 0.60}, {0.0040, 0.60}, { 0.0050, 0.50 }}
    };

    
    struct JointRadiusRel {
        double mcp; // Grundgelenk (Knöchel)
        double pip; // Mittelgelenk (Beim Daumen ist das das IP-Gelenk)
        double dip; // Endgelenk (Beim Daumen nicht vorhanden)
    };

    // Die 14 Werte relativ zur Referenz-Handbreite (68 mm)
    // Der Compiler rechnet die Brüche zur Compile-Zeit exakt aus (z.B. 6.66 / 68.0 = 0.09794...)
    constexpr JointRadiusRel JointCylRadiiRel[5] = {
        // Daumen (F1): MCP = 6.66, IP = 6.55
        { 6.66 / 68.0,  6.55 / 68.0,  0.0 }, 
        // Zeigefinger (F2): MCP = 12.53, PIP = 2.96, DIP = 1.45
        { 12.53 / 68.0, 2.96 / 68.0,  1.45 / 68.0 },
        // Mittelfinger (F3): MCP = 8.27, PIP = 4.47, DIP = 2.19
        { 8.27 / 68.0,  4.47 / 68.0,  2.19 / 68.0 },
        // Ringfinger (F4): MCP = 5.63, PIP = 3.22, DIP = 1.54
        { 5.63 / 68.0,  3.22 / 68.0,  1.54 / 68.0 },
        // Kleiner Finger (F5): MCP = 7.73, PIP = 2.53, DIP = 1.00
        { 7.73 / 68.0,  2.53 / 68.0,  1.00 / 68.0 }
    };

    // ---------------------------------------------------------
    // Struktur für einen einzelnen Knochen (Bone)
    // ---------------------------------------------------------
    struct Bone {
        std::string name;
        MWMath::Point3D startPoint;  // Drehpunkt proximal
        MWMath::Point3D endPoint;    // Drehpunkt distal
        double length;
        MWMath::RotMatrix3x3 orientation; // Lokale X, Y (Knochenachse), Z (Drehachse)
    };

    // ---------------------------------------------------------
    // Struktur für einen kompletten Finger
    // ---------------------------------------------------------
    struct Finger {
        std::vector<Bone> bones; // Meta, Prox, Mid, Dist
    };

    class HandBuilder {
    public:
        double HL;   // HandLength
        double HB;   // HandBreadth
        double Sign; // 1.0 = Rechts

        HandBuilder(double handLength, double handBreadth, double sign = 1.0) 
            : HL(handLength), HB(handBreadth), Sign(sign) {}

        // Hilfsfunktion: Baut die "LookAt" Matrix
        /* MWMath::RotMatrix3x3 buildOrientation(const MWMath::Point3D& start, const MWMath::Point3D& end) {
        
            // 1. Die Y-Achse ist genau die Richtung des Knochens (Endpunkt minus Startpunkt)
            MWMath::Point3D yAxis = (end - start).normed();

            // 2. Die Handfläche liegt grob auf der globalen X-Achse
            MWMath::Point3D palmNormal(1.0, 0.0, 0.0);

            // 3. Z-Achse = Flexionsachse = Quer zum Finger und zur Handfläche
            // Nutzung DEINER cross() Funktion:
            MWMath::Point3D zAxis = palmNormal.cross(yAxis).normed();

            // 4. Echte lokale X-Achse korrigieren (damit alle 3 Achsen 100% senkrecht aufeinander stehen)
            MWMath::Point3D xAxis = yAxis.cross(zAxis).normed();

            // Matrix zusammenbauen (Die Vektoren bilden die Spalten der RotMatrix3x3)
            return MWMath::RotMatrix3x3(
                xAxis.x, yAxis.x, zAxis.x,
                xAxis.y, yAxis.y, zAxis.y,
                xAxis.z, yAxis.z, zAxis.z
            );
        } */


        // idx: 1 = Zeigefinger, 2 = Mittelfinger, 3 = Ringfinger, 4 = Kleiner Finger
        Finger buildFinger(int idx) {
            
            Finger f;
            
            // Dynamischer String für die Namen (idx 1 wird zu "2", etc.)
            std::string fNum = std::to_string(idx + 1); 

            // ==========================================
            // 1. BASISPUNKTE BERECHNEN
            // ==========================================
            
            MWMath::Point3D posMCP(
                0.0, 
                JointPosTable[idx][0] * HL, 
                JointPosTable[idx][1] * HB * Sign
            );

            // Beachte: CMCOffsets Array hat nur 4 Einträge (für Finger 2-5). 
            // Wenn idx = 1 (Zeigefinger), greifen wir auf CMCOffsets[0] zu.
            double scale = HL / RefHandLength;
            MWMath::Point3D posCMC = CMCOffsets[idx - 1] * scale;
            posCMC.z *= Sign;

            // ==========================================
            // 2. METACARPAL BONE BAUEN
            // ==========================================
            Bone meta;
            meta.name = "Metacarpal_" + fNum;
            meta.startPoint = posCMC;
            meta.endPoint = posMCP;
            meta.length = MWMath::distance(posCMC, posMCP);
            meta.orientation = buildOrientation(posCMC, posMCP);
            f.bones.push_back(meta);

            // ==========================================
            // 3. PHALANGES (Fingerglieder) BAUEN
            // ==========================================
            double lenProx = SegLenRatios[idx][1] * HL;
            double lenMid  = SegLenRatios[idx][2] * HL;
            double lenDist = SegLenRatios[idx][3] * HL;

            // Y-Richtung des Fingers (Zweite Spalte der Orientierungs-Matrix)
            MWMath::Point3D fingerDir(meta.orientation.m[0][1], meta.orientation.m[1][1], meta.orientation.m[2][1]);

            // Proximal
            Bone prox;
            prox.name = "Proximal_" + fNum;
            prox.startPoint = meta.endPoint;
            prox.endPoint = prox.startPoint + (fingerDir * lenProx);
            prox.length = lenProx;
            prox.orientation = meta.orientation;
            f.bones.push_back(prox);

            // Medial
            Bone mid;
            mid.name = "Medial_" + fNum;
            mid.startPoint = prox.endPoint;
            mid.endPoint = mid.startPoint + (fingerDir * lenMid);
            mid.length = lenMid;
            mid.orientation = meta.orientation;
            f.bones.push_back(mid);

            // Distal
            Bone dist;
            dist.name = "Distal_" + fNum;
            dist.startPoint = mid.endPoint;
            dist.endPoint = dist.startPoint + (fingerDir * lenDist);
            dist.length = lenDist;
            dist.orientation = meta.orientation;
            f.bones.push_back(dist);

            return f;
        }
    };

}

inline MWMath::Point3D getRotatedPosition(MWMath::Point3D center, double R, double angleDeg) {
    double angleRad = angleDeg * M_PI / 180.0;
    double newX = R * std::cos(angleRad);
    double newY = center.y + R * std::sin(angleRad);
    return MWMath::Point3D(newX, newY, 0.0);
}


//////////////////////////////////////////////////
// ------- DIFFERENT HAND MODEL BUILDER ------- //
inline std::string buildOHandModel(std::vector<std::shared_ptr<SSTissue>>& tissues,std::vector<std::shared_ptr<SSMesh>>& meshes,std::vector<SSMuscle*>& muscles,std::shared_ptr<SSBody>& rootSystem,int numTimeSteps,const SimSettings& cfg,float geometryScaler,const std::vector<double> processParams)
{
    meshes.clear();
    muscles.clear();

    // --- Initialisierung der Hand-Winkel (aus deiner Vorlage) ---
    MWMath::Point3D COLORF1 = MWMath::Point3D(0.50, 0.00, 1.00); 
    MWMath::Point3D COLORF2 = MWMath::Point3D(0.30, 0.30, 1.00); 
    MWMath::Point3D COLORF3 = MWMath::Point3D(0.00, 0.50, 1.00); 
    MWMath::Point3D COLORF4 = MWMath::Point3D(0.00, 0.90, 1.00);
    Hand::MMHandAngles MMHA;

    // --- PARAMETER ---
    double HL = 1.8 * geometryScaler;
    double HB = 0.8 * geometryScaler;
    double scale = HL / Hand::RefHandLength;
    double Sign = 1.0; // 1.0 = Rechts
    float GS = geometryScaler;
    
    std::unordered_map<std::string, SSMesh*> MeshMap;

    std::vector<double> width = { 0.15*0.5 * GS, 0.1*0.5 * GS, 0.07*0.5 * GS, 0.05*0.5 * GS };
    
    double rWF = 1.0; 
    double PF = processParams.size() > 0 ? processParams[0] : 1.0;
    double offF = processParams.size() > 1 ? processParams[1] : 1.0;
    double RF = processParams.size() > 2 ? processParams[2] : 1.0;
    double relT = processParams.size() > 3 ? processParams[3] : 1.0;
    
    std::vector<double> relTorusPos = {0.32*PF, 0.2*PF, 0.1*PF}; //{0.32*PF, 0.2*PF, 0.1*PF}; 
    std::vector<double> relTorusR = {0.14*RF, 0.19*RF, 0.28*RF}; 
    std::vector<double> relTorusr = {relT, relT, relT};
    std::vector<double> off = {0.06*offF, 0.07*offF, 0.18*offF}; 

    double segLen[4];
    double segRadius[4];
    double jointRadius[3];
    std::vector<std::string> fingerNames = {"Thumb", "Index", "Middle", "Ring", "Little"};

    // 1. ROOT
    MWMath::RotMatrix3x3 Ausrichtung = MWMath::axisAngle({0,0,1}, 90.0) * MWMath::axisAngle({0,1,0}, 180.0);
    rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0, 0, 0), Ausrichtung, nullptr);
    tissues.push_back(rootSystem);

    MWMath::Point3D wristOffset(0.0, 0.003567 * scale, -0.003901 * scale * Sign);
    
    auto wristJointAbd = std::make_shared<SSJoint>("Wrist_JointAbd", wristOffset, MWMath::RotMatrix3x3(), rootSystem, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointAbd);
    auto jMeshWrist = std::make_shared<SSEllipsoidMesh>(0.01 * scale, 0.01 * scale, 0.01 * scale, "Mesh_Wrist_Joint", wristJointAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0,0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMeshWrist);
    MeshMap["Mesh_Wrist_Joint"] = jMeshWrist.get();

    auto wristJointFlex = std::make_shared<SSJoint>("Wrist_JointFlex", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointAbd, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointFlex);

    // 2. CARPALS
    auto carpals = std::make_shared<SSBody>("Carpals", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointFlex);
    tissues.push_back(carpals);
    
    double carpalThickness = width[0] * 1.5; 
    double carpalLength = width[0] * 2.5;    
    double carpalWidth = width[0] * 3.5;     
    
    auto meshCarpals = std::make_shared<SSEllipsoidMesh>(carpalThickness, carpalLength, carpalWidth, "Mesh_Carpals", carpals, MWMath::Point3D(0, 0.0*scale, 0), MWMath::RotMatrix3x3(), MWMath::Point3D(0.5, 0.5, 0.5));
    meshes.push_back(meshCarpals);
    MeshMap["Mesh_Carpals"] = meshCarpals.get();

    // 3. WRAPPING SURFACE 
    double flexCylRadius = 0.0115 * scale;
    double flexCylLength = 0.1400 * scale; 
    MWMath::Point3D flexCylOffset(carpalThickness * 0.5 + flexCylRadius, 0.02 * scale, 0.0);

    auto flexorCylMesh = std::make_shared<SSCylinderMesh>(flexCylRadius, flexCylLength, "Mesh_WristFlexorCyl", carpals, flexCylOffset, MWMath::axisAngle({ 1, 0, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));  
    meshes.push_back(flexorCylMesh);
    MeshMap["Mesh_WristFlexorCyl"] = flexorCylMesh.get();

    // FINGER SCHLEIFE
    for (int i = 1; i < 2; i++) {
        int fidx = i;
        std::string cFName = fingerNames[fidx];
        double jAngle;
        std::string prefN = std::to_string(fidx);

        segLen[0] = Hand::SegLenRatios[fidx][0] * HL;
        segLen[1] = Hand::SegLenRatios[fidx][1] * HL;
        segLen[2] = Hand::SegLenRatios[fidx][2] * HL;
        segLen[3] = Hand::SegLenRatios[fidx][3] * HL;

        segRadius[0] = Hand::FingerWrapData[fidx-1][0].radius * scale;
        segRadius[1] = Hand::FingerWrapData[fidx-1][1].radius * scale;
        segRadius[2] = Hand::FingerWrapData[fidx-1][2].radius * scale;
        segRadius[3] = Hand::FingerWrapData[fidx-1][3].radius * scale;

        jointRadius[0] = segRadius[0] * 1.1; 
        jointRadius[1] = segRadius[2] * 1.1; 
        jointRadius[2] = segRadius[3] * 1.1; 

        MWMath::Point3D posMCP(0.0, Hand::JointPosTable[fidx][0] * HL, Hand::JointPosTable[fidx][1] * HB * Sign);
        
        // CMC
        MWMath::Point3D posCMC = Hand::CMCOffsets[fidx - 1] * scale;
        posCMC.z *= Sign;
        MWMath::RotMatrix3x3 forientation = buildOrientation(posMCP, posCMC);
        jAngle = 0.0; 
        std::string jName = "CMC" + prefN + "_JointAbd";
        auto jointCMCAbd = std::make_shared<SSJoint>(jName, posCMC, forientation, carpals, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointCMCAbd);
        auto jMeshCMC = std::make_shared<SSEllipsoidMesh>(width[0], width[0], width[0], "Mesh_CMC"+prefN+"_Joint", jointCMCAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0 , 0 }, 90.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMeshCMC);
        MeshMap["Mesh_CMC"+prefN+"_Joint"] = jMeshCMC.get();

        jAngle = 0.0;
        jName = "CMC" + prefN + "_JointFlex";
        auto jointCMCFlex = std::make_shared<SSJoint>(jName, MWMath::Point3D(0.,0.,0.), MWMath::RotMatrix3x3(), jointCMCAbd, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointCMCFlex);

        // MC SEGMENT
        double flength = MWMath::distance(posCMC, posMCP);
        jName = "MC" + prefN;
        auto bodyMC = std::make_shared<SSBody>(jName, MWMath::Point3D(0.,-flength*0.5,0.), MWMath::RotMatrix3x3(), jointCMCFlex);
        tissues.push_back(bodyMC);
        auto meshMC = std::make_shared<SSEllipsoidMesh>(segRadius[0], segRadius[0], flength*0.5, "Mesh_"+jName, bodyMC, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1, 0, 0 }, 90.0), COLORF1);
        meshes.push_back(meshMC);
        MeshMap["Mesh_" + jName] = meshMC.get();

        double tR0 = segLen[0] * relTorusR[0];
        double tr0 = tR0 * relTorusr[0]; 
        double tPosY0 = -flength * relTorusPos[0]; 
        auto mTorusMC = std::make_shared<SSTorusMesh>(tR0, tr0, "MeshTorus_"+jName, bodyMC, 
             MWMath::Point3D(segLen[0] * off[0], tPosY0, 0), 
             MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mTorusMC);
        MeshMap["MeshTorus_" + jName] = mTorusMC.get();

        // MCP 
        MWMath::Point3D jointPosRelMCP = MWMath::Point3D(0, -flength*0.5, 0); 
        jAngle = 0.0; 
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

        // PP SEGMENT
        MWMath::Point3D posPP = MWMath::Point3D(0, -segLen[1]*0.5, 0);
        jName =  "PP" + prefN;
        auto bodyPP = std::make_shared<SSBody>(jName, posPP, MWMath::RotMatrix3x3(), jointMCPFlex);
        tissues.push_back(bodyPP);
        auto meshPP = std::make_shared<SSEllipsoidMesh>(segRadius[1], segRadius[1], segLen[1]*0.5, "Mesh_"+jName, bodyPP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1, 0, 0 }, 90.0), COLORF2);
        meshes.push_back(meshPP);
        MeshMap["Mesh_" + jName] = meshPP.get();

        double tR1 = segLen[1] * relTorusR[1];
        double tr1 = tR1 * relTorusr[1]; 
        double tPosY1 = -segLen[1] *relTorusPos[1]; 
        auto mTorusPP = std::make_shared<SSTorusMesh>(tR1, tr1, "MeshTorus_"+jName, bodyPP, 
             MWMath::Point3D(segLen[1] * off[1], tPosY1, 0), 
             MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0.8, 0.3, 0));
        meshes.push_back(mTorusPP);
        MeshMap["MeshTorus_" + jName] = mTorusPP.get();

        // PIP 
        MWMath::Point3D jointPosRelPP = MWMath::Point3D(0, -segLen[1]*0.5, 0); 
        jAngle = MMHA.PIP_M[fidx];
        jName = "PIP" + prefN + "_Joint";
        auto jointPIP = std::make_shared<SSJoint>(jName, jointPosRelPP, MWMath::RotMatrix3x3(), bodyPP, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointPIP);
        auto MeshPIP = std::make_shared<SSCylinderMesh>(jointRadius[1], jointRadius[1]*3, "Mesh_"+jName, jointPIP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 0, 1, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(MeshPIP);
        MeshMap["Mesh_" + jName] = MeshPIP.get();

        // MP SEGMENT
        MWMath::Point3D posMP = MWMath::Point3D(0., -segLen[2]*0.5, 0);
        jName =  "MP" + prefN;
        auto bodyMP = std::make_shared<SSBody>(jName, posMP, MWMath::RotMatrix3x3(), jointPIP);
        tissues.push_back(bodyMP);
        auto meshMP = std::make_shared<SSEllipsoidMesh>(segRadius[2], segRadius[2], segLen[2]*0.5, "Mesh_"+jName, bodyMP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0, 0 }, 90.0), COLORF3);
        meshes.push_back(meshMP);
        MeshMap["Mesh_" + jName] = meshMP.get();

        double tR2 = segLen[2] * relTorusR[2];
        double tr2 = tR2 * relTorusr[2]; 
        double tPosY2 = -segLen[2] * relTorusPos[2]; 
        auto mTorusMP = std::make_shared<SSTorusMesh>(tR2, tr2, "MeshTorus_"+jName, bodyMP, 
             MWMath::Point3D(segLen[2] * off[2], tPosY2, 0), 
             MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0.8, 0.3, 0.2));
        meshes.push_back(mTorusMP);
        MeshMap["MeshTorus_" + jName] = mTorusMP.get();

        // DIP 
        MWMath::Point3D jointPosRelMP = MWMath::Point3D(0,-segLen[2]*0.5, 0); 
        jAngle = MMHA.DIP_M[fidx];
        jName = "DIP" + prefN + "_Joint";
        auto jointDIP = std::make_shared<SSJoint>(jName, jointPosRelMP, MWMath::RotMatrix3x3(), bodyMP, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointDIP);
        auto MeshDIP = std::make_shared<SSCylinderMesh>(jointRadius[2], jointRadius[2]*3, "Mesh_"+jName, jointDIP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 0, 1, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(MeshDIP);
        MeshMap["Mesh_" + jName] = MeshDIP.get();

        // DP SEGMENT
        MWMath::Point3D posDP = MWMath::Point3D(0., -segLen[3]*0.5, 0.);
        jName =  "DP" + prefN;
        auto bodyDP = std::make_shared<SSBody>(jName, posDP, MWMath::RotMatrix3x3(), jointDIP);
        tissues.push_back(bodyDP);
        auto meshDP = std::make_shared<SSEllipsoidMesh>(segRadius[3], segRadius[3], segLen[3]*0.5, "Mesh_"+jName, bodyDP, MWMath::Point3D(0., 0., 0.), MWMath::axisAngle({ 1,0,0}, 90.0), COLORF4);
        meshes.push_back(meshDP);
        MeshMap["Mesh_" + jName] = meshDP.get(); 
    }

    // INITIALES UPDATE & DISKRETISIERUNG
    for (auto& m : meshes) { 
        m->InitializeMesh(); 
        m->discretizeMesh(cfg.discretization); 
    }
    rootSystem->update(0);

    // MUSKEL
    int numPoints = cfg.muscleNumPoints.empty() ? 25 : cfg.muscleNumPoints[0];
    std::vector<std::string> refBodyNames = { "Mesh_WristFlexorCyl", "MeshTorus_MC1", "Mesh_MCP1_Joint", "MeshTorus_PP1", "Mesh_PIP1_Joint", "MeshTorus_MP1", "Mesh_DIP1_Joint"};
    MWMath::Point3D startOffset = MWMath::Point3D(0., 0.5, 0.0);
    MWMath::Point3D endOffset = MWMath::Point3D(0.05, 0.0, 0.0);
    
    SSMuscle* flexor = new SSMuscle("Flexor", numPoints, rootSystem.get(), startOffset, getRefMesh("Mesh_DP1", MeshMap), endOffset);

    for (auto& m : getRefMeshes(refBodyNames, MeshMap)) {
        flexor->meshPtrs.push_back(m);
    }

    flexor->createMusclePointsComplexPath();
    flexor->updateMusclePointsParentsLocal();
    muscles.push_back(flexor);


    std::string returnString = "";
    for (const auto param : processParams) {
        returnString += std::to_string(param) + "\t";
    }
    return returnString;
}

// --- Use only 1 torus and two limitting cylinders
inline std::string buildOHandModelCyl(std::vector<std::shared_ptr<SSTissue>>& tissues,std::vector<std::shared_ptr<SSMesh>>& meshes,std::vector<SSMuscle*>& muscles,std::shared_ptr<SSBody>& rootSystem,int numTimeSteps,const SimSettings& cfg,float geometryScaler,const std::vector<double> processParams)
{
    meshes.clear();
    muscles.clear();

    // --- Initialisierung der Hand-Winkel (aus deiner Vorlage) ---
    MWMath::Point3D COLORF1 = MWMath::Point3D(0.50, 0.00, 1.00); 
    MWMath::Point3D COLORF2 = MWMath::Point3D(0.30, 0.30, 1.00); 
    MWMath::Point3D COLORF3 = MWMath::Point3D(0.00, 0.50, 1.00); 
    MWMath::Point3D COLORF4 = MWMath::Point3D(0.00, 0.90, 1.00);
    Hand::MMHandAngles MMHA;

    // --- PARAMETER ---
    double HL = 1.8 * geometryScaler;
    double HB = 0.8 * geometryScaler;
    double scale = HL / Hand::RefHandLength;
    double Sign = 1.0; // 1.0 = Rechts
    float GS = geometryScaler;
    
    std::unordered_map<std::string, SSMesh*> MeshMap;

    std::vector<double> width = { 0.15*0.5 * GS, 0.1*0.5 * GS, 0.07*0.5 * GS, 0.05*0.5 * GS };
    
    double rWF = 1.0; 
    double PF = processParams.size() > 0 ? processParams[0] : 1.0;
    double offF = processParams.size() > 1 ? processParams[1] : 1.0;
    double RF = processParams.size() > 2 ? processParams[2] : 1.0;
    double relT = processParams.size() > 3 ? processParams[3] : 0.7;
    
    std::vector<double> relTorusPos = {0.25*PF, 0.1*PF, 0.0*PF}; //{0.32*PF, 0.2*PF, 0.1*PF}; 
    std::vector<double> relTorusR = {0.1*RF, 0.19*RF, 0.28*RF}; 
    std::vector<double> relTorusr = {relT, relT, relT};
    std::vector<double> off = {0.06*1.0, 0.07*4.0, 0.18*2.5}; 

    double segLen[4];
    double segRadius[4];
    double jointRadius[3];
    std::vector<std::string> fingerNames = {"Thumb", "Index", "Middle", "Ring", "Little"};

    // 1. ROOT
    MWMath::RotMatrix3x3 Ausrichtung = MWMath::axisAngle({0,0,1}, 90.0) * MWMath::axisAngle({0,1,0}, 180.0);
    rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0, 0, 0), Ausrichtung, nullptr);
    tissues.push_back(rootSystem);

    MWMath::Point3D wristOffset(0.0, 0.003567 * scale, -0.003901 * scale * Sign);
    
    auto wristJointAbd = std::make_shared<SSJoint>("Wrist_JointAbd", wristOffset, MWMath::RotMatrix3x3(), rootSystem, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointAbd);
    auto jMeshWrist = std::make_shared<SSEllipsoidMesh>(0.01 * scale, 0.01 * scale, 0.01 * scale, "Mesh_Wrist_Joint", wristJointAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0,0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMeshWrist);
    MeshMap["Mesh_Wrist_Joint"] = jMeshWrist.get();

    auto wristJointFlex = std::make_shared<SSJoint>("Wrist_JointFlex", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointAbd, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointFlex);

    // 2. CARPALS
    auto carpals = std::make_shared<SSBody>("Carpals", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointFlex);
    tissues.push_back(carpals);
    
    double carpalThickness = width[0] * 1.5; 
    double carpalLength = width[0] * 2.5;    
    double carpalWidth = width[0] * 3.5;     
    
    auto meshCarpals = std::make_shared<SSEllipsoidMesh>(carpalThickness, carpalLength, carpalWidth, "Mesh_Carpals", carpals, MWMath::Point3D(0, 0.0*scale, 0), MWMath::RotMatrix3x3(), MWMath::Point3D(0.5, 0.5, 0.5));
    meshes.push_back(meshCarpals);
    MeshMap["Mesh_Carpals"] = meshCarpals.get();

    // 3. WRAPPING SURFACE 
    double flexCylRadius = 0.0115 * scale;
    double flexCylLength = 0.1400 * scale; 
    MWMath::Point3D flexCylOffset(carpalThickness * 0.5 + flexCylRadius, 0.02 * scale, 0.0);

    auto flexorCylMesh = std::make_shared<SSCylinderMesh>(flexCylRadius, flexCylLength, "Mesh_WristFlexorCyl", carpals, flexCylOffset, MWMath::axisAngle({ 1, 0, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));  
    meshes.push_back(flexorCylMesh);
    MeshMap["Mesh_WristFlexorCyl"] = flexorCylMesh.get();

    // FINGER SCHLEIFE
    for (int i = 1; i < 2; i++) {
        int fidx = i;
        //std::string cFName = fingerNames[fidx];
        double jAngle;
        std::string prefN = std::to_string(fidx);

        segLen[0] = Hand::SegLenRatios[fidx][0] * HL;
        segLen[1] = Hand::SegLenRatios[fidx][1] * HL;
        segLen[2] = Hand::SegLenRatios[fidx][2] * HL;
        segLen[3] = Hand::SegLenRatios[fidx][3] * HL;

        segRadius[0] = Hand::FingerWrapData[fidx-1][0].radius * scale;
        segRadius[1] = Hand::FingerWrapData[fidx-1][1].radius * scale;
        segRadius[2] = Hand::FingerWrapData[fidx-1][2].radius * scale;
        segRadius[3] = Hand::FingerWrapData[fidx-1][3].radius * scale;

        jointRadius[0] = segRadius[0] * 1.1; 
        jointRadius[1] = segRadius[2] * 1.1; 
        jointRadius[2] = segRadius[3] * 1.1; 

        MWMath::Point3D posMCP(0.0, Hand::JointPosTable[fidx][0] * HL, Hand::JointPosTable[fidx][1] * HB * Sign);
        
        // CMC
        MWMath::Point3D posCMC = Hand::CMCOffsets[fidx - 1] * scale;
        posCMC.z *= Sign;
        MWMath::RotMatrix3x3 forientation = buildOrientation(posMCP, posCMC);
        jAngle = 0.0; 
        std::string jName = "CMC" + prefN + "_JointAbd";
        auto jointCMCAbd = std::make_shared<SSJoint>(jName, posCMC, forientation, carpals, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointCMCAbd);
        auto jMeshCMC = std::make_shared<SSEllipsoidMesh>(width[0], width[0], width[0], "Mesh_CMC"+prefN+"_Joint", jointCMCAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0 , 0 }, 90.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMeshCMC);
        MeshMap["Mesh_CMC"+prefN+"_Joint"] = jMeshCMC.get();

        jAngle = 0.0;
        jName = "CMC" + prefN + "_JointFlex";
        auto jointCMCFlex = std::make_shared<SSJoint>(jName, MWMath::Point3D(0.,0.,0.), MWMath::RotMatrix3x3(), jointCMCAbd, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointCMCFlex);

        // MC SEGMENT
        double flength = MWMath::distance(posCMC, posMCP);
        jName = "MC" + prefN;
        auto bodyMC = std::make_shared<SSBody>(jName, MWMath::Point3D(0.,-flength*0.5,0.), MWMath::RotMatrix3x3(), jointCMCFlex);
        tissues.push_back(bodyMC);
        auto meshMC = std::make_shared<SSEllipsoidMesh>(segRadius[0], segRadius[0], flength*0.5, "Mesh_"+jName, bodyMC, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1, 0, 0 }, 90.0), COLORF1);
        meshes.push_back(meshMC);
        MeshMap["Mesh_" + jName] = meshMC.get();

        double tR0 = segLen[0] * relTorusR[0];
        double tr0 = tR0 * relTorusr[0]; 
        double tPosY0 = -flength * relTorusPos[0]; 
        auto mTorusMC = std::make_shared<SSTorusMesh>(tR0, tr0, "MeshTorus_"+jName, bodyMC, 
             MWMath::Point3D(segLen[0] * off[0], tPosY0, 0), 
             MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mTorusMC);
        MeshMap["MeshTorus_" + jName] = mTorusMC.get();

        // MCP 
        MWMath::Point3D jointPosRelMCP = MWMath::Point3D(0, -flength*0.5, 0); 
        jAngle = 0.0; 
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

        // PP SEGMENT
        MWMath::Point3D posPP = MWMath::Point3D(0, -segLen[1]*0.5, 0);
        jName =  "PP" + prefN;
        auto bodyPP = std::make_shared<SSBody>(jName, posPP, MWMath::RotMatrix3x3(), jointMCPFlex);
        tissues.push_back(bodyPP);
        auto meshPP = std::make_shared<SSEllipsoidMesh>(segRadius[1], segRadius[1], segLen[1]*0.5, "Mesh_"+jName, bodyPP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1, 0, 0 }, 90.0), COLORF2);
        meshes.push_back(meshPP);
        MeshMap["Mesh_" + jName] = meshPP.get();

        double tR1 = segLen[1] * relTorusR[1];
        double tr1 = tR1 * relTorusr[1]; 
        double tPosY1 = -segLen[1] *relTorusPos[1]; 
        // SSCylinderMesh Parameter: (Radius, Höhe, Name, Parent, Pos, Rot, Color)
        // Wir nutzen tr1 als Radius und machen ihn etwas breiter als den Finger (width[1] * 1.5)
        auto mTorusPP = std::make_shared<SSCylinderMesh>(tr1, tr1*2.0, "MeshTorus_"+jName, bodyPP, 
                MWMath::Point3D(segLen[1] * off[1], tPosY1, 0), 
                MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.3, 0));
        meshes.push_back(mTorusPP);
        MeshMap["MeshTorus_" + jName] = mTorusPP.get();

        // PIP 
        MWMath::Point3D jointPosRelPP = MWMath::Point3D(0, -segLen[1]*0.5, 0); 
        jAngle = MMHA.PIP_M[fidx];
        jName = "PIP" + prefN + "_Joint";
        auto jointPIP = std::make_shared<SSJoint>(jName, jointPosRelPP, MWMath::RotMatrix3x3(), bodyPP, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointPIP);
        auto MeshPIP = std::make_shared<SSCylinderMesh>(jointRadius[1], jointRadius[1]*3, "Mesh_"+jName, jointPIP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 0, 1, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(MeshPIP);
        MeshMap["Mesh_" + jName] = MeshPIP.get();

        // MP SEGMENT
        MWMath::Point3D posMP = MWMath::Point3D(0., -segLen[2]*0.5, 0);
        jName =  "MP" + prefN;
        auto bodyMP = std::make_shared<SSBody>(jName, posMP, MWMath::RotMatrix3x3(), jointPIP);
        tissues.push_back(bodyMP);
        auto meshMP = std::make_shared<SSEllipsoidMesh>(segRadius[2], segRadius[2], segLen[2]*0.5, "Mesh_"+jName, bodyMP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0, 0 }, 90.0), COLORF3);
        meshes.push_back(meshMP);
        MeshMap["Mesh_" + jName] = meshMP.get();

        double tR2 = segLen[2] * relTorusR[2];
        double tr2 = tR2 * relTorusr[2]; 
        double tPosY2 = -segLen[2] * relTorusPos[2]; 
        auto mTorusMP = std::make_shared<SSCylinderMesh>(tr2, tr2*2.0, "MeshTorus_"+jName, bodyMP, 
                MWMath::Point3D(segLen[2] * off[2], tPosY2, 0), 
                MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.3, 0.2));
        meshes.push_back(mTorusMP);
        MeshMap["MeshTorus_" + jName] = mTorusMP.get();

        // DIP 
        MWMath::Point3D jointPosRelMP = MWMath::Point3D(0,-segLen[2]*0.5, 0); 
        jAngle = MMHA.DIP_M[fidx];
        jName = "DIP" + prefN + "_Joint";
        auto jointDIP = std::make_shared<SSJoint>(jName, jointPosRelMP, MWMath::RotMatrix3x3(), bodyMP, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointDIP);
        auto MeshDIP = std::make_shared<SSCylinderMesh>(jointRadius[2], jointRadius[2]*3, "Mesh_"+jName, jointDIP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 0, 1, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(MeshDIP);
        MeshMap["Mesh_" + jName] = MeshDIP.get();

        // DP SEGMENT
        MWMath::Point3D posDP = MWMath::Point3D(0., -segLen[3]*0.5, 0.);
        jName =  "DP" + prefN;
        auto bodyDP = std::make_shared<SSBody>(jName, posDP, MWMath::RotMatrix3x3(), jointDIP);
        tissues.push_back(bodyDP);
        auto meshDP = std::make_shared<SSEllipsoidMesh>(segRadius[3], segRadius[3], segLen[3]*0.5, "Mesh_"+jName, bodyDP, MWMath::Point3D(0., 0., 0.), MWMath::axisAngle({ 1,0,0}, 90.0), COLORF4);
        meshes.push_back(meshDP);
        MeshMap["Mesh_" + jName] = meshDP.get(); 
    }

    // INITIALES UPDATE & DISKRETISIERUNG
    for (auto& m : meshes) { 
        m->InitializeMesh(); 
        m->discretizeMesh(cfg.discretization); 
    }
    rootSystem->update(0);

    // MUSKEL
    int numPoints = cfg.muscleNumPoints.empty() ? 25 : cfg.muscleNumPoints[0];
    std::vector<std::string> refBodyNames = { "MeshTorus_MC1", "Mesh_MCP1_Joint", "MeshTorus_PP1", "Mesh_PIP1_Joint", "MeshTorus_MP1", "Mesh_DIP1_Joint"};
    MWMath::Point3D startOffset = MWMath::Point3D(0., 0.5, 0.0);
    MWMath::Point3D endOffset = MWMath::Point3D(0.05, 0.0, 0.0);
    
    SSMuscle* flexor = new SSMuscle("Flexor", numPoints, rootSystem.get(), startOffset, getRefMesh("Mesh_DP1", MeshMap), endOffset);

    for (auto& m : getRefMeshes(refBodyNames, MeshMap)) {
        flexor->meshPtrs.push_back(m);
    }

    flexor->createMusclePointsComplexPath();
    flexor->updateMusclePointsParentsLocal();
    muscles.push_back(flexor);


    std::string returnString = "";
    for (const auto param : processParams) {
        returnString += std::to_string(param) + "\t";
    }
    return returnString;
}

// --- based on the original OFingerSimple 
inline std::string buildOHandModelOld(std::vector<std::shared_ptr<SSTissue>>& tissues,std::vector<std::shared_ptr<SSMesh>>& meshes,std::vector<SSMuscle*>& muscles,std::shared_ptr<SSBody>& rootSystem,int numTimeSteps,const SimSettings& cfg,float geometryScaler,const std::vector<double> processParams)
{
    meshes.clear();
    muscles.clear();

    // --- Initialisierung der Hand-Winkel ---
    MWMath::Point3D COLORF1 = MWMath::Point3D(0.50, 0.00, 1.00); 
    MWMath::Point3D COLORF2 = MWMath::Point3D(0.30, 0.30, 1.00); 
    MWMath::Point3D COLORF3 = MWMath::Point3D(0.00, 0.50, 1.00); 
    MWMath::Point3D COLORF4 = MWMath::Point3D(0.00, 0.90, 1.00);
    MWMath::Point3D COLORFDEACTIVE = MWMath::Point3D(0.5, 0.5, 0.5); // Graue Farbe aus Vorlage
    Hand::MMHandAngles MMHA;

    // --- PARAMETER ---
    double HL = 1.8 * geometryScaler;
    double HB = 0.8 * geometryScaler;
    double scale = HL / Hand::RefHandLength;
    double Sign = 1.0; // 1.0 = Rechts
    float GS = geometryScaler;
    
    std::unordered_map<std::string, SSMesh*> MeshMap;

    // --- FINGER PARAMETER AUS "OFINGER_SIMPLE_OnlyTorusSmall" ---
    std::vector<double> fLength = {0.3482 * HL, 0.2027 * HL, 0.1175 * HL, 0.0882 * HL};
    double L1 = fLength[0]; // Länge Segment 1
    double L2 = fLength[1]; // Länge Segment 2
    double L3 = fLength[2]; // Länge Segment 3
    double L4 = fLength[3]; // Länge Segment 4
    
    double rWF = 0.5; 
    std::vector<double> width = {0.15*GS*rWF, 0.1*GS*rWF, 0.07*GS*rWF, 0.05*GS*rWF}; // Dicke

    // Parameter aus processParams ziehen (für die Torus-Logik)
    double PF = processParams.size() > 0 ? processParams[0] : 1.0;
    double offF = processParams.size() > 1 ? processParams[1] : 1.0;
    double RF = processParams.size() > 2 ? processParams[2] : 1.0;
    double relT = processParams.size() > 3 ? processParams[3] : 1.0;

    // ==============================================================================
    // KORREKTUR: KEINE DIVISION DURCH rWF MEHR HIER! 
    // Da width (und somit mesh->B) schon durch rWF skaliert ist, 
    // würde eine erneute Division den Torus wieder massiv/doppelt so groß machen!
    // ==============================================================================
    std::vector<double> relTorusPos = {0.1*PF, -0.08*PF, -0.18*PF}; 
    std::vector<double> relTorusR = {1.3*RF, 1.3*RF, 1.6*RF}; 
    std::vector<double> relTorusr = {1.1*relT, 1.1*relT, 1.4*relT}; 

    bool bShowBody = true;
    double off = 0.8 * offF; 

    // 1. ROOT
    MWMath::RotMatrix3x3 Ausrichtung = MWMath::axisAngle({0,0,1}, 90.0) * MWMath::axisAngle({0,1,0}, 180.0);
    rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0, 0, 0), Ausrichtung, nullptr);
    tissues.push_back(rootSystem);

    MWMath::Point3D wristOffset(0.0, 0.003567 * scale, -0.003901 * scale * Sign);
    
    auto wristJointAbd = std::make_shared<SSJoint>("Wrist_JointAbd", wristOffset, MWMath::RotMatrix3x3(), rootSystem, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointAbd);
    auto jMeshWrist = std::make_shared<SSEllipsoidMesh>(0.01 * scale, 0.01 * scale, 0.01 * scale, "Mesh_Wrist_Joint", wristJointAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0,0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMeshWrist);
    MeshMap["Mesh_Wrist_Joint"] = jMeshWrist.get();

    auto wristJointFlex = std::make_shared<SSJoint>("Wrist_JointFlex", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointAbd, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointFlex);

    // 2. CARPALS
    auto carpals = std::make_shared<SSBody>("Carpals", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointFlex);
    tissues.push_back(carpals);
    
    double carpalThickness = width[0] * 1.5; 
    double carpalLength = width[0] * 2.5;    
    double carpalWidth = width[0] * 3.5;     
    
    auto meshCarpals = std::make_shared<SSEllipsoidMesh>(carpalThickness, carpalLength, carpalWidth, "Mesh_Carpals", carpals, MWMath::Point3D(0, 0.0*scale, 0), MWMath::RotMatrix3x3(), MWMath::Point3D(0.5, 0.5, 0.5));
    meshes.push_back(meshCarpals);
    MeshMap["Mesh_Carpals"] = meshCarpals.get();

    // 3. WRAPPING SURFACE (Karpaltunnel)
    double flexCylRadius = 0.0115 * scale;
    double flexCylLength = 0.1400 * scale; 
    MWMath::Point3D flexCylOffset(carpalThickness * 0.5 + flexCylRadius, 0.02 * scale, 0.0);

    auto flexorCylMesh = std::make_shared<SSCylinderMesh>(flexCylRadius, flexCylLength, "Mesh_WristFlexorCyl", carpals, flexCylOffset, MWMath::axisAngle({ 1, 0, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));  
    meshes.push_back(flexorCylMesh);
    MeshMap["Mesh_WristFlexorCyl"] = flexorCylMesh.get();

    // ==============================================================================
    // FINGER AUFBAU (Index Finger = fidx 1)
    // ==============================================================================
    int fidx = 1; // Wir bauen den Zeigefinger
    std::string cFName = "Index";
    std::string prefN = std::to_string(fidx);
    double jAngle;

    // Basis-Position für CMC
    MWMath::Point3D posMCP(0.0, Hand::JointPosTable[fidx][0] * HL, Hand::JointPosTable[fidx][1] * HB * Sign);
    MWMath::Point3D posCMC = Hand::CMCOffsets[fidx - 1] * scale;
    posCMC.z *= Sign;
    MWMath::RotMatrix3x3 forientation = buildOrientation(posMCP, posCMC);

    // CMC Joint Abd
    jAngle = 0.0; 
    std::string jName = "CMC" + prefN + "_JointAbd";
    auto jointCMCAbd = std::make_shared<SSJoint>(jName, posCMC, forientation, carpals, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(jointCMCAbd);
    auto jMeshCMC = std::make_shared<SSEllipsoidMesh>(width[0], width[0], width[0], "Mesh_CMC"+prefN+"_Joint", jointCMCAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0 , 0 }, 90.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMeshCMC);
    MeshMap["Mesh_CMC"+prefN+"_Joint"] = jMeshCMC.get();

    // CMC Joint Flex
    jAngle = 0.0;
    jName = "CMC" + prefN + "_JointFlex";
    auto jointCMCFlex = std::make_shared<SSJoint>(jName, MWMath::Point3D(0.,0.,0.), MWMath::RotMatrix3x3(), jointCMCAbd, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(jointCMCFlex);

    // ==============================================================================
    // SEGMENT 1 (Proximal / MC)
    // ==============================================================================
    double flength = MWMath::distance(posCMC, posMCP);
    jName = "MC" + prefN;
    
    auto body1 = std::make_shared<SSBody>(jName, MWMath::Point3D(0.,-flength*0.5,0.), MWMath::RotMatrix3x3(), jointCMCFlex);
    tissues.push_back(body1);
    
    // Mesh 1
    auto mesh1 = std::make_shared<SSEllipsoidMesh>(width[0], width[0], flength*0.5, "Mesh_"+jName, body1, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), COLORFDEACTIVE);
    if (bShowBody) meshes.push_back(mesh1);
    MeshMap["Mesh_" + jName] = mesh1.get();
    
    // Torus 1
    auto mTorus1 = std::make_shared<SSTorusMesh>(mesh1->B * relTorusR[0], mesh1->B * relTorusr[0], 
         "MeshTorus_"+jName, body1, MWMath::Point3D(width[0]*off, -L1 * relTorusPos[0], 0), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(1, 0, 0));
    meshes.push_back(mTorus1);
    MeshMap["MeshTorus_" + jName] = mTorus1.get();

    // MCP Joint (Joint 1)
    MWMath::Point3D jointPosRel1 = MWMath::Point3D(0, -flength*0.5, 0); 
    jAngle = 0.0; 
    jName = "MCP" + prefN + "_JointAbd";
    auto jointMCPAbd = std::make_shared<SSJoint>(jName, jointPosRel1, MWMath::RotMatrix3x3(), body1, jAngle, MWMath::Point3D(1, 0, 0), numTimeSteps);
    tissues.push_back(jointMCPAbd);
    
    // KORREKTUR: jointSize OHNE /rWF !
    double jointSize1 = width[0]*1.1;
    auto jMeshMCP = std::make_shared<SSEllipsoidMesh>(jointSize1, jointSize1, jointSize1, "Mesh_MCP"+prefN+"_Joint", jointMCPAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0 , 0 }, 90.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMeshMCP);
    MeshMap["Mesh_MCP"+prefN+"_Joint"] = jMeshMCP.get();

    jAngle = MMHA.MCP_FE_M[fidx];
    jName = "MCP" + prefN + "_JointFlex";
    auto jointMCPFlex = std::make_shared<SSJoint>(jName, MWMath::Point3D(0.,0.,0.), MWMath::RotMatrix3x3(), jointMCPAbd, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(jointMCPFlex);

    // ==============================================================================
    // SEGMENT 2 (PP)
    // ==============================================================================
    jName =  "PP" + prefN;
    auto body2 = std::make_shared<SSBody>(jName, MWMath::Point3D(0, -L2*0.5, 0), MWMath::RotMatrix3x3(), jointMCPFlex);
    tissues.push_back(body2);
    
    auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2*0.5, "Mesh_"+jName, body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), COLORFDEACTIVE);
    if (bShowBody) meshes.push_back(mesh2);
    MeshMap["Mesh_" + jName] = mesh2.get();
    
    // Torus 2
    auto mTorus2 = std::make_shared<SSTorusMesh>(mesh2->B * relTorusR[1], mesh2->B * relTorusr[1], 
         "MeshTorus_"+jName, body2, MWMath::Point3D(width[1]*off, -L2 * relTorusPos[1], 0), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
    meshes.push_back(mTorus2);
    MeshMap["MeshTorus_" + jName] = mTorus2.get();

    // PIP Joint (Joint 2)
    MWMath::Point3D jointPosRel2 = MWMath::Point3D(0, -L2*0.5, 0); 
    jAngle = MMHA.PIP_M[fidx];
    jName = "PIP" + prefN + "_Joint";
    auto jointPIP = std::make_shared<SSJoint>(jName, jointPosRel2, MWMath::RotMatrix3x3(), body2, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(jointPIP);
    
    // KORREKTUR: jointSize OHNE /rWF !
    double jointSize2 = width[1]*1.1;
    auto MeshPIP = std::make_shared<SSEllipsoidMesh>(jointSize2, jointSize2, jointSize2, "Mesh_"+jName, jointPIP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 0, 1, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(MeshPIP);
    MeshMap["Mesh_" + jName] = MeshPIP.get();

    // ==============================================================================
    // SEGMENT 3 (MP)
    // ==============================================================================
    jName =  "MP" + prefN;
    auto body3 = std::make_shared<SSBody>(jName, MWMath::Point3D(0., -L3*0.5, 0), MWMath::RotMatrix3x3(), jointPIP);
    tissues.push_back(body3);
    
    auto mesh3 = std::make_shared<SSEllipsoidMesh>(width[2], width[2], L3*0.5, "Mesh_"+jName, body3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), COLORFDEACTIVE);
    if (bShowBody) meshes.push_back(mesh3);
    MeshMap["Mesh_" + jName] = mesh3.get();

    // Torus 3
    auto mTorus3 = std::make_shared<SSTorusMesh>(mesh3->B * relTorusR[2], mesh3->B * relTorusr[2], 
         "MeshTorus_"+jName, body3, MWMath::Point3D(width[2]*off, -L3 * relTorusPos[2], 0), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
    meshes.push_back(mTorus3);
    MeshMap["MeshTorus_" + jName] = mTorus3.get();

    // DIP Joint (Joint 3)
    MWMath::Point3D jointPosRel3 = MWMath::Point3D(0,-L3*0.5, 0); 
    jAngle = MMHA.DIP_M[fidx];
    jName = "DIP" + prefN + "_Joint";
    auto jointDIP = std::make_shared<SSJoint>(jName, jointPosRel3, MWMath::RotMatrix3x3(), body3, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(jointDIP);
    
    // KORREKTUR: jointSize OHNE /rWF !
    double jointSize3 = width[2]*1.1;
    auto MeshDIP = std::make_shared<SSEllipsoidMesh>(jointSize3, jointSize3, jointSize3, "Mesh_"+jName, jointDIP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 0, 1, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(MeshDIP);
    MeshMap["Mesh_" + jName] = MeshDIP.get();

    // ==============================================================================
    // SEGMENT 4 (DP)
    // ==============================================================================
    jName =  "DP" + prefN;
    auto body4 = std::make_shared<SSBody>(jName, MWMath::Point3D(0., -L4*0.5, 0.), MWMath::RotMatrix3x3(), jointDIP);
    tissues.push_back(body4);
    
    auto mesh4 = std::make_shared<SSEllipsoidMesh>(width[3], width[3], L4*0.5, "Mesh_"+jName, body4, MWMath::Point3D(0., 0., 0.), MWMath::axisAngle({1,0,0}, 90.0), COLORFDEACTIVE);
    if (bShowBody) meshes.push_back(mesh4);
    MeshMap["Mesh_" + jName] = mesh4.get(); 

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
    
    // Die Namen der Meshes für den Muskelpfad sammeln
    std::vector<std::string> refBodyNames = { 
        "Mesh_WristFlexorCyl", 
        "MeshTorus_MC1", 
        "Mesh_MCP1_Joint", 
        "MeshTorus_PP1", 
        "Mesh_PIP1_Joint", 
        "MeshTorus_MP1", 
        "Mesh_DIP1_Joint"
    };

    MWMath::Point3D startOffset = MWMath::Point3D(0., 0.5, 0.0);
    MWMath::Point3D endOffset = MWMath::Point3D(width[3], -width[3] * 1.1, 0.0);
    
    SSMuscle* flexor = new SSMuscle("Flexor", numPoints, rootSystem.get(), startOffset, getRefMesh("Mesh_DP1", MeshMap), endOffset);

    // Füge die Hindernisse aus der Liste zum Muskel hinzu
    for (auto& m : getRefMeshes(refBodyNames, MeshMap)) {
        flexor->meshPtrs.push_back(m);
    }

    flexor->createMusclePointsComplexPath();
    flexor->updateMusclePointsParentsLocal();
    muscles.push_back(flexor);

    std::string returnString = "";
    for (const auto param : processParams) {
        returnString += std::to_string(param) + "\t";
    }
    return returnString;
}


inline std::string buildOHandModelTorusAsJoint(std::vector<std::shared_ptr<SSTissue>>& tissues,std::vector<std::shared_ptr<SSMesh>>& meshes,std::vector<SSMuscle*>& muscles,std::shared_ptr<SSBody>& rootSystem,int numTimeSteps,const SimSettings& cfg,float geometryScaler,const std::vector<double> processParams)
{
    meshes.clear();
    muscles.clear();

    // --- Initialisierung der Hand-Winkel (aus deiner Vorlage) ---
    MWMath::Point3D COLORF1 = MWMath::Point3D(0.50, 0.00, 1.00); 
    MWMath::Point3D COLORF2 = MWMath::Point3D(0.30, 0.30, 1.00); 
    MWMath::Point3D COLORF3 = MWMath::Point3D(0.00, 0.50, 1.00); 
    MWMath::Point3D COLORF4 = MWMath::Point3D(0.00, 0.90, 1.00);
    Hand::MMHandAngles MMHA;

    // --- PARAMETER ---
    double HL = 1.8 * geometryScaler;
    double HB = 0.8 * geometryScaler;
    double scale = HL / Hand::RefHandLength;
    double Sign = 1.0; // 1.0 = Rechts
    float GS = geometryScaler;
    
    std::unordered_map<std::string, SSMesh*> MeshMap;

    std::vector<double> width = { 0.15*0.5 * GS, 0.1*0.5 * GS, 0.07*0.5 * GS, 0.05*0.5 * GS };
    std::vector<double> torusOffAngle = {50.0, 60.0, 50.0}; // {60.0, 60.0, 40.0}; // Winkel für die Torus-Offset-Rotation
    std::vector<double> torusRotR = {1.4, 1.6, 1.5}; // Rotation der Torus-Orientierung (für mehr "Drehung" bei den distalen Gelenken)
    
    double rWF = 1.0; 
    double PF = processParams.size() > 0 ? processParams[0] : 1.0;
    double offF = processParams.size() > 1 ? processParams[1] : 1.0;
    double RF = processParams.size() > 2 ? processParams[2] : 1.0;
    double relT = processParams.size() > 3 ? processParams[3] : 1.0;
    
    std::vector<double> relTorusPos = {0.32*PF, 0.2*PF, 0.1*PF}; 
    std::vector<double> relTorusR = {0.14*RF, 0.19*RF, 0.28*RF}; 
    std::vector<double> relTorusr = {relT, relT, relT};
    std::vector<double> off = {0.06*offF, 0.07*offF, 0.18*offF}; 

    double segLen[4];
    double segRadius[4];
    double jointRadius[3];
    std::vector<std::string> fingerNames = {"Thumb", "Index", "Middle", "Ring", "Little"};

    // 1. ROOT
    MWMath::RotMatrix3x3 Ausrichtung = MWMath::axisAngle({0,0,1}, 90.0) * MWMath::axisAngle({0,1,0}, 180.0);
    rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0, 0, 0), Ausrichtung, nullptr);
    tissues.push_back(rootSystem);

    MWMath::Point3D wristOffset(0.0, 0.003567 * scale, -0.003901 * scale * Sign);
    
    auto wristJointAbd = std::make_shared<SSJoint>("Wrist_JointAbd", wristOffset, MWMath::RotMatrix3x3(), rootSystem, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointAbd);
    auto jMeshWrist = std::make_shared<SSEllipsoidMesh>(0.01 * scale, 0.01 * scale, 0.01 * scale, "Mesh_Wrist_Joint", wristJointAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0,0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMeshWrist);
    MeshMap["Mesh_Wrist_Joint"] = jMeshWrist.get();

    auto wristJointFlex = std::make_shared<SSJoint>("Wrist_JointFlex", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointAbd, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointFlex);

    // 2. CARPALS
    auto carpals = std::make_shared<SSBody>("Carpals", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointFlex);
    tissues.push_back(carpals);
    
    double carpalThickness = width[0] * 1.5; 
    double carpalLength = width[0] * 2.5;    
    double carpalWidth = width[0] * 3.5;     
    
    auto meshCarpals = std::make_shared<SSEllipsoidMesh>(carpalThickness, carpalLength, carpalWidth, "Mesh_Carpals", carpals, MWMath::Point3D(0, 0.0*scale, 0), MWMath::RotMatrix3x3(), MWMath::Point3D(0.5, 0.5, 0.5));
    meshes.push_back(meshCarpals);
    MeshMap["Mesh_Carpals"] = meshCarpals.get();

    // 3. WRAPPING SURFACE 
    double flexCylRadius = 0.0115 * scale;
    double flexCylLength = 0.1400 * scale; 
    MWMath::Point3D flexCylOffset(carpalThickness * 0.5 + flexCylRadius, 0.02 * scale, 0.0);

    auto flexorCylMesh = std::make_shared<SSCylinderMesh>(flexCylRadius, flexCylLength, "Mesh_WristFlexorCyl", carpals, flexCylOffset, MWMath::axisAngle({ 1, 0, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));  
    meshes.push_back(flexorCylMesh);
    MeshMap["Mesh_WristFlexorCyl"] = flexorCylMesh.get();

    // FINGER SCHLEIFE
    for (int i = 1; i < 2; i++) {
        int fidx = i;
        std::string cFName = fingerNames[fidx];
        double jAngle;
        std::string prefN = std::to_string(fidx);

        segLen[0] = Hand::SegLenRatios[fidx][0] * HL;
        segLen[1] = Hand::SegLenRatios[fidx][1] * HL;
        segLen[2] = Hand::SegLenRatios[fidx][2] * HL;
        segLen[3] = Hand::SegLenRatios[fidx][3] * HL;

        segRadius[0] = Hand::FingerWrapData[fidx-1][0].radius * scale;
        segRadius[1] = Hand::FingerWrapData[fidx-1][1].radius * scale;
        segRadius[2] = Hand::FingerWrapData[fidx-1][2].radius * scale;
        segRadius[3] = Hand::FingerWrapData[fidx-1][3].radius * scale;

        jointRadius[0] = segRadius[0] * 1.1; 
        jointRadius[1] = segRadius[2] * 1.1; 
        jointRadius[2] = segRadius[3] * 1.1; 

        MWMath::Point3D posMCP(0.0, Hand::JointPosTable[fidx][0] * HL, Hand::JointPosTable[fidx][1] * HB * Sign);
        
        // CMC
        MWMath::Point3D posCMC = Hand::CMCOffsets[fidx - 1] * scale;
        posCMC.z *= Sign;
        MWMath::RotMatrix3x3 forientation = buildOrientation(posMCP, posCMC);
        jAngle = 0.0; 
        std::string jName = "CMC" + prefN + "_JointAbd";
        auto jointCMCAbd = std::make_shared<SSJoint>(jName, posCMC, forientation, carpals, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointCMCAbd);
        auto jMeshCMC = std::make_shared<SSEllipsoidMesh>(width[0], width[0], width[0], "Mesh_CMC"+prefN+"_Joint", jointCMCAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0 , 0 }, 90.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMeshCMC);
        MeshMap["Mesh_CMC"+prefN+"_Joint"] = jMeshCMC.get();

        jAngle = 0.0;
        jName = "CMC" + prefN + "_JointFlex";
        auto jointCMCFlex = std::make_shared<SSJoint>(jName, MWMath::Point3D(0.,0.,0.), MWMath::RotMatrix3x3(), jointCMCAbd, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointCMCFlex);

        // ==========================================
        // MC SEGMENT
        // ==========================================
        double flength = MWMath::distance(posCMC, posMCP);
        jName = "MC" + prefN;
        auto bodyMC = std::make_shared<SSBody>(jName, MWMath::Point3D(0.,-flength*0.5,0.), MWMath::RotMatrix3x3(), jointCMCFlex);
        tissues.push_back(bodyMC);
        auto meshMC = std::make_shared<SSEllipsoidMesh>(segRadius[0], segRadius[0], flength*0.5, "Mesh_"+jName, bodyMC, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1, 0, 0 }, 90.0), COLORF1);
        meshes.push_back(meshMC);
        MeshMap["Mesh_" + jName] = meshMC.get();

        // --- NEU: TORUS MC ---
        // Schlauch exakt so groß wie das Gelenk (tr0 = jointRadius)
        double tr0 = jointRadius[0] * relT; 
        double tR0 = jointRadius[0] * torusRotR[0]; 
        auto mTorusMC = std::make_shared<SSTorusMesh>(tR0, tr0, "MeshTorus_MC"+prefN, bodyMC, 
             getRotatedPosition(MWMath::Point3D(0., -flength*0.5, 0.), tR0, torusOffAngle[0]), // Sitzt exakt im Gelenkzentrum
             MWMath::axisAngle({0,0,1}, torusOffAngle[0]) * MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mTorusMC);
        MeshMap["MeshTorus_MC" + prefN] = mTorusMC.get();

         
        // MCP Joint + Mesh
        MWMath::Point3D jointPosRelMCP = MWMath::Point3D(0, -flength*0.5, 0); 
        jAngle = 0.0; 
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
        double tr1 = jointRadius[1] * relT;
        double tR1 = jointRadius[1] * torusRotR[1]; 
        auto mTorusPP = std::make_shared<SSTorusMesh>(tR1, tr1, "MeshTorus_PP"+prefN, bodyPP, 
             getRotatedPosition(MWMath::Point3D(0., -segLen[1]*0.5, 0.), tR1, torusOffAngle[1]), // Sitzt exakt im Gelenkzentrum
             MWMath::axisAngle({0,0,1}, torusOffAngle[1]) * MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0.8, 0.3, 0));
        meshes.push_back(mTorusPP);
        MeshMap["MeshTorus_PP" + prefN] = mTorusPP.get();

        
        // PIP Joint + Mesh
        MWMath::Point3D jointPosRelPP = MWMath::Point3D(0, -segLen[1]*0.5, 0); 
        jAngle = MMHA.PIP_M[fidx];
        jName = "PIP" + prefN + "_Joint";
        auto jointPIP = std::make_shared<SSJoint>(jName, jointPosRelPP, MWMath::RotMatrix3x3(), bodyPP, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointPIP);
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
        double tr2 = jointRadius[2] * relT; 
        double tR2 = jointRadius[2] * torusRotR[2]; 
        auto mTorusMP = std::make_shared<SSTorusMesh>(tR2, tr2, "MeshTorus_MP"+prefN, bodyMP, 
             getRotatedPosition(MWMath::Point3D(0., -segLen[2]*0.5, 0.), tR2, torusOffAngle[2]), // Sitzt exakt im Gelenkzentrum
             MWMath::axisAngle({0,0,1}, torusOffAngle[2]) * MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0.8, 0.3, 0.2));
        meshes.push_back(mTorusMP);
        MeshMap["MeshTorus_MP" + prefN] = mTorusMP.get();

        // DIP Joint + Mesh
        MWMath::Point3D jointPosRelMP = MWMath::Point3D(0,-segLen[2]*0.5, 0); 
        jAngle = MMHA.DIP_M[fidx];
        jName = "DIP" + prefN + "_Joint";
        auto jointDIP = std::make_shared<SSJoint>(jName, jointPosRelMP, MWMath::RotMatrix3x3(), bodyMP, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointDIP);
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

    // INITIALES UPDATE & DISKRETISIERUNG
    for (auto& m : meshes) { 
        m->InitializeMesh(); 
        m->discretizeMesh(cfg.discretization); 
    }
    rootSystem->update(0);

    // MUSKEL
    int numPoints = cfg.muscleNumPoints.empty() ? 25 : cfg.muscleNumPoints[0];
    //std::vector<std::string> refBodyNames = { "Mesh_MC1", "MeshTorus_MC1", "Mesh_PP1", "MeshTorus_PP1", "Mesh_MP1", "MeshTorus_MP1"};
    std::vector<std::string> refBodyNames = { "MeshTorus_MC1", "MeshTorus_PP1", "MeshTorus_MP1"};
    //std::vector<std::string> refBodyNames = { "Mesh_WristFlexorCyl", "MeshTorus_MC1", "Mesh_MCP1_Joint", "MeshTorus_PP1", "Mesh_PIP1_Joint", "MeshTorus_MP1", "Mesh_DIP1_Joint"};
    MWMath::Point3D startOffset = MWMath::Point3D(0., 0.5, 0.0);
    MWMath::Point3D endOffset = MWMath::Point3D(0.05, 0.0, 0.0);
    
    SSMuscle* flexor = new SSMuscle("Flexor", numPoints, rootSystem.get(), startOffset, getRefMesh("Mesh_DP1", MeshMap), endOffset);

    for (auto& m : getRefMeshes(refBodyNames, MeshMap)) {
        flexor->meshPtrs.push_back(m);
    }

    flexor->createMusclePointsComplexPath();
    flexor->updateMusclePointsParentsLocal();
    muscles.push_back(flexor);


    std::string returnString = "";
    for (const auto param : processParams) {
        returnString += std::to_string(param) + "\t";
    }
    return returnString;
}

inline std::string buildOHandModelTorusAsJointKreuzband(std::vector<std::shared_ptr<SSTissue>>& tissues,std::vector<std::shared_ptr<SSMesh>>& meshes,std::vector<SSMuscle*>& muscles,std::shared_ptr<SSBody>& rootSystem,int numTimeSteps,const SimSettings& cfg,float geometryScaler,const std::vector<double> processParams)
{
    meshes.clear();
    muscles.clear();

    // --- Initialisierung der Hand-Winkel (aus deiner Vorlage) ---
    MWMath::Point3D COLORF1 = MWMath::Point3D(0.50, 0.00, 1.00); 
    MWMath::Point3D COLORF2 = MWMath::Point3D(0.30, 0.30, 1.00); 
    MWMath::Point3D COLORF3 = MWMath::Point3D(0.00, 0.50, 1.00); 
    MWMath::Point3D COLORF4 = MWMath::Point3D(0.00, 0.90, 1.00);
    Hand::MMHandAngles MMHA;

    // --- PARAMETER ---
    double HL = 1.8 * geometryScaler;
    double HB = 0.8 * geometryScaler;
    double scale = HL / Hand::RefHandLength;
    double Sign = 1.0; // 1.0 = Rechts
    float GS = geometryScaler;
    
    std::unordered_map<std::string, SSMesh*> MeshMap;

    std::vector<double> width = { 0.15*0.5 * GS, 0.1*0.5 * GS, 0.07*0.5 * GS, 0.05*0.5 * GS };
    std::vector<double> torusOffAngle = {50.0, 60.0, 50.0}; // {60.0, 60.0, 40.0}; // Winkel für die Torus-Offset-Rotation
    std::vector<double> torusRotR = {1.4, 1.6, 1.5}; // Rotation der Torus-Orientierung (für mehr "Drehung" bei den distalen Gelenken)
    
    double rWF = 1.0; 
    double PF = processParams.size() > 0 ? processParams[0] : 1.0;
    double offF = processParams.size() > 1 ? processParams[1] : 1.0;
    double RF = processParams.size() > 2 ? processParams[2] : 1.0;
    double relT = processParams.size() > 3 ? processParams[3] : 1.0;
    
    std::vector<double> relTorusPos = {0.32*PF, 0.2*PF, 0.1*PF}; 
    std::vector<double> relTorusR = {0.14*RF, 0.19*RF, 0.28*RF}; 
    std::vector<double> relTorusr = {relT, relT, relT};
    std::vector<double> off = {0.06*offF, 0.07*offF, 0.18*offF}; 

    double segLen[4];
    double segRadius[4];
    double jointRadius[3];
    std::vector<std::string> fingerNames = {"Thumb", "Index", "Middle", "Ring", "Little"};

    // 1. ROOT
    MWMath::RotMatrix3x3 Ausrichtung = MWMath::axisAngle({0,0,1}, 90.0) * MWMath::axisAngle({0,1,0}, 180.0);
    rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0, 0, 0), Ausrichtung, nullptr);
    tissues.push_back(rootSystem);

    MWMath::Point3D wristOffset(0.0, 0.003567 * scale, -0.003901 * scale * Sign);
    
    auto wristJointAbd = std::make_shared<SSJoint>("Wrist_JointAbd", wristOffset, MWMath::RotMatrix3x3(), rootSystem, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointAbd);
    auto jMeshWrist = std::make_shared<SSEllipsoidMesh>(0.01 * scale, 0.01 * scale, 0.01 * scale, "Mesh_Wrist_Joint", wristJointAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0,0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMeshWrist);
    MeshMap["Mesh_Wrist_Joint"] = jMeshWrist.get();

    auto wristJointFlex = std::make_shared<SSJoint>("Wrist_JointFlex", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointAbd, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointFlex);

    // 2. CARPALS
    auto carpals = std::make_shared<SSBody>("Carpals", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointFlex);
    tissues.push_back(carpals);
    
    double carpalThickness = width[0] * 1.5; 
    double carpalLength = width[0] * 2.5;    
    double carpalWidth = width[0] * 3.5;     
    
    auto meshCarpals = std::make_shared<SSEllipsoidMesh>(carpalThickness, carpalLength, carpalWidth, "Mesh_Carpals", carpals, MWMath::Point3D(0, 0.0*scale, 0), MWMath::RotMatrix3x3(), MWMath::Point3D(0.5, 0.5, 0.5));
    meshes.push_back(meshCarpals);
    MeshMap["Mesh_Carpals"] = meshCarpals.get();

    // 3. WRAPPING SURFACE 
    double flexCylRadius = 0.0115 * scale;
    double flexCylLength = 0.1400 * scale; 
    MWMath::Point3D flexCylOffset(carpalThickness * 0.5 + flexCylRadius, 0.02 * scale, 0.0);

    auto flexorCylMesh = std::make_shared<SSCylinderMesh>(flexCylRadius, flexCylLength, "Mesh_WristFlexorCyl", carpals, flexCylOffset, MWMath::axisAngle({ 1, 0, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));  
    meshes.push_back(flexorCylMesh);
    MeshMap["Mesh_WristFlexorCyl"] = flexorCylMesh.get();

    // FINGER SCHLEIFE
     // FINGER SCHLEIFE
    for (int i = 1; i < 2; i++) {
        int fidx = i;
        std::string cFName = fingerNames[fidx];
        double jAngle;
        std::string prefN = std::to_string(fidx);

        segLen[0] = Hand::SegLenRatios[fidx][0] * HL;
        segLen[1] = Hand::SegLenRatios[fidx][1] * HL;
        segLen[2] = Hand::SegLenRatios[fidx][2] * HL;
        segLen[3] = Hand::SegLenRatios[fidx][3] * HL;

        segRadius[0] = Hand::FingerWrapData[fidx-1][0].radius * scale;
        segRadius[1] = Hand::FingerWrapData[fidx-1][1].radius * scale;
        segRadius[2] = Hand::FingerWrapData[fidx-1][2].radius * scale;
        segRadius[3] = Hand::FingerWrapData[fidx-1][3].radius * scale;

        jointRadius[0] = segRadius[0] * 1.1; 
        jointRadius[1] = segRadius[2] * 1.1; 
        jointRadius[2] = segRadius[3] * 1.1; 

        MWMath::Point3D posMCP(0.0, Hand::JointPosTable[fidx][0] * HL, Hand::JointPosTable[fidx][1] * HB * Sign);
        
        // ==========================================
        // CMC Joint
        // ==========================================
        MWMath::Point3D posCMC = Hand::CMCOffsets[fidx - 1] * scale;
        posCMC.z *= Sign;
        MWMath::RotMatrix3x3 forientation = buildOrientation(posMCP, posCMC);
        jAngle = 0.0; 
        std::string jName = "CMC" + prefN + "_JointAbd";
        auto jointCMCAbd = std::make_shared<SSJoint>(jName, posCMC, forientation, carpals, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointCMCAbd);
        auto jMeshCMC = std::make_shared<SSEllipsoidMesh>(width[0], width[0], width[0], "Mesh_CMC"+prefN+"_Joint", jointCMCAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0 , 0 }, 90.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMeshCMC);
        MeshMap["Mesh_CMC"+prefN+"_Joint"] = jMeshCMC.get();

        jAngle = 0.0;
        jName = "CMC" + prefN + "_JointFlex";
        auto jointCMCFlex = std::make_shared<SSJoint>(jName, MWMath::Point3D(0.,0.,0.), MWMath::RotMatrix3x3(), jointCMCAbd, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointCMCFlex);

        // ==========================================
        // MC SEGMENT (Proximal Phalanx / Handknochen)
        // ==========================================
        double flength = MWMath::distance(posCMC, posMCP);
        jName = "MC" + prefN;
        auto bodyMC = std::make_shared<SSBody>(jName, MWMath::Point3D(0.,-flength*0.5,0.), MWMath::RotMatrix3x3(), jointCMCFlex);
        tissues.push_back(bodyMC);
        auto meshMC = std::make_shared<SSEllipsoidMesh>(segRadius[0], segRadius[0], flength*0.5, "Mesh_"+jName, bodyMC, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1, 0, 0 }, 90.0), COLORF1);
        meshes.push_back(meshMC);
        MeshMap["Mesh_" + jName] = meshMC.get();

        // --- TORUS MC (Ringband) ---
        double tr0 = jointRadius[0] * relT; 
        double tR0 = jointRadius[0] * torusRotR[0]; 
        auto mTorusMC = std::make_shared<SSTorusMesh>(tR0, tr0, "MeshTorus_MC"+prefN, bodyMC, 
             getRotatedPosition(MWMath::Point3D(0., -flength*0.5, 0.), tR0, torusOffAngle[0]), 
             MWMath::axisAngle({0,0,1}, torusOffAngle[0]) * MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mTorusMC);
        MeshMap["MeshTorus_MC" + prefN] = mTorusMC.get();
         
        // ==========================================
        // MCP Joint
        // ==========================================
        MWMath::Point3D jointPosRelMCP = MWMath::Point3D(0, -flength*0.5, 0); 
        jAngle = 0.0; 
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
        // PP SEGMENT (Proximal Phalanx)
        // ==========================================
        MWMath::Point3D posPP = MWMath::Point3D(0, -segLen[1]*0.5, 0);
        jName =  "PP" + prefN;
        auto bodyPP = std::make_shared<SSBody>(jName, posPP, MWMath::RotMatrix3x3(), jointMCPFlex);
        tissues.push_back(bodyPP);
        auto meshPP = std::make_shared<SSEllipsoidMesh>(segRadius[1], segRadius[1], segLen[1]*0.5, "Mesh_"+jName, bodyPP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1, 0, 0 }, 90.0), COLORF2);
        meshes.push_back(meshPP);
        MeshMap["Mesh_" + jName] = meshPP.get();

        


        // --- TORUS PP (Normales Ringband unten am PIP Gelenk) ---
        double tr1 = jointRadius[1] * relT;
        double tR1 = jointRadius[1] * torusRotR[1]; 
        auto mTorusPP = std::make_shared<SSTorusMesh>(tR1, tr1, "MeshTorus_PP"+prefN, bodyPP, 
             getRotatedPosition(MWMath::Point3D(0., -segLen[1]*0.5, 0.), tR1, torusOffAngle[1]), // -L/2 (Distales Ende von PP)
             MWMath::axisAngle({0,0,1}, torusOffAngle[1]) * MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0.8, 0.3, 0));
        meshes.push_back(mTorusPP);
        MeshMap["MeshTorus_PP" + prefN] = mTorusPP.get();

        // --- KREUZBAND TEIL 1 (Distal am MCP Gelenk, gehört aber zu PP!) ---
        // Dieser Torus sitzt oben am PP-Knochen (pos: L/2) und ist negativ geneigt.
        auto mTorusPP_Cruciate = std::make_shared<SSTorusMesh>(tR1, tr1, "MeshTorus_Cruciate_MCP"+prefN, bodyPP, 
             getRotatedPosition(MWMath::Point3D(0., segLen[1]*0.5, 0.), tR1, -torusOffAngle[1]), // +L/2 (Proximales Ende von PP) und negativer Winkel
             MWMath::axisAngle({0,0,1}, -torusOffAngle[1]) * MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(1, 0.5, 0)); // Orange Farbe zur Unterscheidung
        meshes.push_back(mTorusPP_Cruciate);
        MeshMap["MeshTorus_Cruciate_MCP" + prefN] = mTorusPP_Cruciate.get();

        
        // ==========================================
        // PIP Joint
        // ==========================================
        MWMath::Point3D jointPosRelPP = MWMath::Point3D(0, -segLen[1]*0.5, 0); 
        jAngle = MMHA.PIP_M[fidx];
        jName = "PIP" + prefN + "_Joint";
        auto jointPIP = std::make_shared<SSJoint>(jName, jointPosRelPP, MWMath::RotMatrix3x3(), bodyPP, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointPIP);
        auto MeshPIP = std::make_shared<SSCylinderMesh>(jointRadius[1], jointRadius[1]*3, "Mesh_"+jName, jointPIP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 0, 1, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(MeshPIP);
        MeshMap["Mesh_" + jName] = MeshPIP.get();

        // ==========================================
        // MP SEGMENT (Middle Phalanx)
        // ==========================================
        MWMath::Point3D posMP = MWMath::Point3D(0., -segLen[2]*0.5, 0);
        jName =  "MP" + prefN;
        auto bodyMP = std::make_shared<SSBody>(jName, posMP, MWMath::RotMatrix3x3(), jointPIP);
        tissues.push_back(bodyMP);
        auto meshMP = std::make_shared<SSEllipsoidMesh>(segRadius[2], segRadius[2], segLen[2]*0.5, "Mesh_"+jName, bodyMP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0, 0 }, 90.0), COLORF3);
        meshes.push_back(meshMP);
        MeshMap["Mesh_" + jName] = meshMP.get();


        // --- TORUS MP (Normales Ringband unten am DIP Gelenk) ---
        double tr2 = jointRadius[2] * relT; 
        double tR2 = jointRadius[2] * torusRotR[2]; 
        auto mTorusMP = std::make_shared<SSTorusMesh>(tR2, tr2, "MeshTorus_MP"+prefN, bodyMP, 
             getRotatedPosition(MWMath::Point3D(0., -segLen[2]*0.5, 0.), tR2, torusOffAngle[2]), 
             MWMath::axisAngle({0,0,1}, torusOffAngle[2]) * MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0.8, 0.3, 0.2));
        meshes.push_back(mTorusMP);
        MeshMap["MeshTorus_MP" + prefN] = mTorusMP.get();

        // --- KREUZBAND TEIL 2 (Distal am PIP Gelenk, gehört zu MP!) ---
        auto mTorusMP_Cruciate = std::make_shared<SSTorusMesh>(tR2, tr2, "MeshTorus_Cruciate_PIP"+prefN, bodyMP, 
             getRotatedPosition(MWMath::Point3D(0., segLen[2]*0.5, 0.), tR2, -torusOffAngle[2]), // +L/2 und negativer Winkel
             MWMath::axisAngle({0,0,1}, -torusOffAngle[2]) * MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0.8, 0.5, 0.2));
        meshes.push_back(mTorusMP_Cruciate);
        MeshMap["MeshTorus_Cruciate_PIP" + prefN] = mTorusMP_Cruciate.get();

        // ==========================================
        // DIP Joint
        // ==========================================
        MWMath::Point3D jointPosRelMP = MWMath::Point3D(0,-segLen[2]*0.5, 0); 
        jAngle = MMHA.DIP_M[fidx];
        jName = "DIP" + prefN + "_Joint";
        auto jointDIP = std::make_shared<SSJoint>(jName, jointPosRelMP, MWMath::RotMatrix3x3(), bodyMP, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointDIP);
        auto MeshDIP = std::make_shared<SSCylinderMesh>(jointRadius[2], jointRadius[2]*3, "Mesh_"+jName, jointDIP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 0, 1, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(MeshDIP);
        MeshMap["Mesh_" + jName] = MeshDIP.get();

        // ==========================================
        // DP SEGMENT (Distal Phalanx)
        // ==========================================
        MWMath::Point3D posDP = MWMath::Point3D(0., -segLen[3]*0.5, 0.);
        jName =  "DP" + prefN;
        auto bodyDP = std::make_shared<SSBody>(jName, posDP, MWMath::RotMatrix3x3(), jointDIP);
        tissues.push_back(bodyDP);
        auto meshDP = std::make_shared<SSEllipsoidMesh>(segRadius[3], segRadius[3], segLen[3]*0.5, "Mesh_"+jName, bodyDP, MWMath::Point3D(0., 0., 0.), MWMath::axisAngle({ 1,0,0}, 90.0), COLORF4);
        meshes.push_back(meshDP);
        MeshMap["Mesh_" + jName] = meshDP.get(); 

        /* // --- KREUZBAND TEIL 3 (Distal am DIP Gelenk, gehört zu DP!) ---
        auto mTorusDP_Cruciate = std::make_shared<SSTorusMesh>(tR2, tr2, "MeshTorus_Cruciate_DIP"+prefN, bodyDP, 
             getRotatedPosition(MWMath::Point3D(0., segLen[3]*0.5, 0.), tR2, -torusOffAngle[2]), // +L/2 und negativer Winkel
             MWMath::axisAngle({0,0,1}, -torusOffAngle[2]) * MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0.8, 0.6, 0.4));
        meshes.push_back(mTorusDP_Cruciate);
        MeshMap["MeshTorus_Cruciate_DIP" + prefN] = mTorusDP_Cruciate.get(); */
    }
    

    // INITIALES UPDATE & DISKRETISIERUNG
    for (auto& m : meshes) { 
        m->InitializeMesh(); 
        m->discretizeMesh(cfg.discretization); 
    }
    rootSystem->update(0);

    // MUSKEL
    int numPoints = cfg.muscleNumPoints.empty() ? 25 : cfg.muscleNumPoints[0];
    //std::vector<std::string> refBodyNames = { "Mesh_MC1", "MeshTorus_MC1", "Mesh_PP1", "MeshTorus_PP1", "Mesh_MP1", "MeshTorus_MP1"};
    std::vector<std::string> refBodyNames = {
        "Mesh_WristFlexorCyl",
        //"Mesh_MC1", 
        "Mesh_WristFlexorCyl", 
        // --- MCP Bereich ---
        "MeshTorus_MC1",               // Normales Ringband VOR dem Gelenk
        "MeshTorus_Cruciate_MCP1",     // Kreuzband NACH dem Gelenk
        // --- PIP Bereich ---
        "MeshTorus_PP1",               // Normales Ringband VOR dem Gelenk
        "MeshTorus_Cruciate_PIP1",     // Kreuzband NACH dem Gelenk
        // --- DIP Bereich ---
        "MeshTorus_MP1",               // Normales Ringband VOR dem Gelenk
        };
    //std::vector<std::string> refBodyNames = { "Mesh_WristFlexorCyl", "MeshTorus_MC1", "Mesh_MCP1_Joint", "MeshTorus_PP1", "Mesh_PIP1_Joint", "MeshTorus_MP1", "Mesh_DIP1_Joint"};
    MWMath::Point3D startOffset = MWMath::Point3D(0., 0.5, 0.0);
    MWMath::Point3D endOffset = MWMath::Point3D(0.05, 0.0, 0.0);
    
    SSMuscle* flexor = new SSMuscle("Flexor", numPoints, rootSystem.get(), startOffset, getRefMesh("Mesh_DP1", MeshMap), endOffset);

    for (auto& m : getRefMeshes(refBodyNames, MeshMap)) {
        flexor->meshPtrs.push_back(m);
    }

    flexor->createMusclePointsComplexPath();
    flexor->updateMusclePointsParentsLocal();
    muscles.push_back(flexor);


    std::string returnString = "";
    for (const auto param : processParams) {
        returnString += std::to_string(param) + "\t";
    }
    return returnString;
}



inline std::string buildOHandModelCylEllHole(std::vector<std::shared_ptr<SSTissue>>& tissues,std::vector<std::shared_ptr<SSMesh>>& meshes,std::vector<SSMuscle*>& muscles,std::shared_ptr<SSBody>& rootSystem,int numTimeSteps,const SimSettings& cfg,float geometryScaler,const std::vector<double> processParams)
{
    meshes.clear();
    muscles.clear();

    // --- Initialisierung der Hand-Winkel (aus deiner Vorlage) ---
    MWMath::Point3D COLORF1 = MWMath::Point3D(0.50, 0.00, 1.00); 
    MWMath::Point3D COLORF2 = MWMath::Point3D(0.30, 0.30, 1.00); 
    MWMath::Point3D COLORF3 = MWMath::Point3D(0.00, 0.50, 1.00); 
    MWMath::Point3D COLORF4 = MWMath::Point3D(0.00, 0.90, 1.00);
    Hand::MMHandAngles MMHA;

    // --- PARAMETER ---
    double HL = 1.8 * geometryScaler;
    double HB = 0.8 * geometryScaler;
    double scale = HL / Hand::RefHandLength;
    double Sign = 1.0; // 1.0 = Rechts
    float GS = geometryScaler;
    
    std::unordered_map<std::string, SSMesh*> MeshMap;

    std::vector<double> width = { 0.15*0.5 * GS, 0.1*0.5 * GS, 0.07*0.5 * GS, 0.05*0.5 * GS };
    
    double rWF = 1.0; 
    double PF = processParams.size() > 0 ? processParams[0] : 1.0;
    double offF = processParams.size() > 1 ? processParams[1] : 1.0;
    double RF = processParams.size() > 2 ? processParams[2] : 1.0;
    double relT = processParams.size() > 3 ? processParams[3] : 0.7;
    
    std::vector<double> relTorusPos = {0.25*PF, 0.1*PF, 0.0*PF}; 
    std::vector<double> relTorusR = {0.1*RF, 0.19*RF, 0.28*RF}; 
    std::vector<double> relTorusr = {relT, relT, relT};
    std::vector<double> off = {0.06*1.0, 0.07*4.0, 0.18*2.5}; 

    double segLen[4];
    double segRadius[4];
    double jointRadius[3];
    std::vector<std::string> fingerNames = {"Thumb", "Index", "Middle", "Ring", "Little"};

    // 1. ROOT
    MWMath::RotMatrix3x3 Ausrichtung = MWMath::axisAngle({0,0,1}, 90.0) * MWMath::axisAngle({0,1,0}, 180.0);
    rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0, 0, 0), Ausrichtung, nullptr);
    tissues.push_back(rootSystem);

    MWMath::Point3D wristOffset(0.0, 0.003567 * scale, -0.003901 * scale * Sign);
    
    auto wristJointAbd = std::make_shared<SSJoint>("Wrist_JointAbd", wristOffset, MWMath::RotMatrix3x3(), rootSystem, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointAbd);
    auto jMeshWrist = std::make_shared<SSEllipsoidMesh>(0.01 * scale, 0.01 * scale, 0.01 * scale, "Mesh_Wrist_Joint", wristJointAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0,0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMeshWrist);
    MeshMap["Mesh_Wrist_Joint"] = jMeshWrist.get();

    auto wristJointFlex = std::make_shared<SSJoint>("Wrist_JointFlex", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointAbd, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointFlex);

    // 2. CARPALS
    auto carpals = std::make_shared<SSBody>("Carpals", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointFlex);
    tissues.push_back(carpals);
    
    double carpalThickness = width[0] * 1.5; 
    double carpalLength = width[0] * 2.5;    
    double carpalWidth = width[0] * 3.5;     
    
    auto meshCarpals = std::make_shared<SSEllipsoidMesh>(carpalThickness, carpalLength, carpalWidth, "Mesh_Carpals", carpals, MWMath::Point3D(0, 0.0*scale, 0), MWMath::RotMatrix3x3(), MWMath::Point3D(0.5, 0.5, 0.5));
    meshes.push_back(meshCarpals);
    MeshMap["Mesh_Carpals"] = meshCarpals.get();

    // 3. WRAPPING SURFACE 
    double flexCylRadius = 0.0115 * scale;
    double flexCylLength = 0.1400 * scale; 
    MWMath::Point3D flexCylOffset(carpalThickness * 0.5 + flexCylRadius, 0.02 * scale, 0.0);

    auto flexorCylMesh = std::make_shared<SSCylinderMesh>(flexCylRadius, flexCylLength, "Mesh_WristFlexorCyl", carpals, flexCylOffset, MWMath::axisAngle({ 1, 0, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));  
    meshes.push_back(flexorCylMesh);
    MeshMap["Mesh_WristFlexorCyl"] = flexorCylMesh.get();

    // FINGER SCHLEIFE
    for (int i = 1; i < 2; i++) {
        int fidx = i;
        double jAngle;
        std::string prefN = std::to_string(fidx);

        segLen[0] = Hand::SegLenRatios[fidx][0] * HL;
        segLen[1] = Hand::SegLenRatios[fidx][1] * HL;
        segLen[2] = Hand::SegLenRatios[fidx][2] * HL;
        segLen[3] = Hand::SegLenRatios[fidx][3] * HL;

        segRadius[0] = Hand::FingerWrapData[fidx-1][0].radius * scale;
        segRadius[1] = Hand::FingerWrapData[fidx-1][1].radius * scale;
        segRadius[2] = Hand::FingerWrapData[fidx-1][2].radius * scale;
        segRadius[3] = Hand::FingerWrapData[fidx-1][3].radius * scale;

        jointRadius[0] = segRadius[0] * 1.1; 
        jointRadius[1] = segRadius[2] * 1.1; 
        jointRadius[2] = segRadius[3] * 1.1; 

        MWMath::Point3D posMCP(0.0, Hand::JointPosTable[fidx][0] * HL, Hand::JointPosTable[fidx][1] * HB * Sign);
        
        // CMC
        MWMath::Point3D posCMC = Hand::CMCOffsets[fidx - 1] * scale;
        posCMC.z *= Sign;
        MWMath::RotMatrix3x3 forientation = buildOrientation(posMCP, posCMC);
        jAngle = 0.0; 
        std::string jName = "CMC" + prefN + "_JointAbd";
        auto jointCMCAbd = std::make_shared<SSJoint>(jName, posCMC, forientation, carpals, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointCMCAbd);
        auto jMeshCMC = std::make_shared<SSEllipsoidMesh>(width[0], width[0], width[0], "Mesh_CMC"+prefN+"_Joint", jointCMCAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0 , 0 }, 90.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMeshCMC);
        MeshMap["Mesh_CMC"+prefN+"_Joint"] = jMeshCMC.get();

        jAngle = 0.0;
        jName = "CMC" + prefN + "_JointFlex";
        auto jointCMCFlex = std::make_shared<SSJoint>(jName, MWMath::Point3D(0.,0.,0.), MWMath::RotMatrix3x3(), jointCMCAbd, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointCMCFlex);

        // ==========================================
        // MC SEGMENT
        // ==========================================
        double flength = MWMath::distance(posCMC, posMCP);
        jName = "MC" + prefN;
        auto bodyMC = std::make_shared<SSBody>(jName, MWMath::Point3D(0.,-flength*0.5,0.), MWMath::RotMatrix3x3(), jointCMCFlex);
        tissues.push_back(bodyMC);
        auto meshMC = std::make_shared<SSEllipsoidMesh>(segRadius[0], segRadius[0]*1.4, flength*0.5, "Mesh_"+jName, bodyMC, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1, 0, 0 }, 90.0), COLORF1);
        meshes.push_back(meshMC);
        MeshMap["Mesh_" + jName] = meshMC.get();

        // --- NEU: 3 ELLIPSOIDE ALS LOCH (Ersetzt den Torus am MC) ---
        // Basisposition für das "Loch"
        double tPosY0 = -flength * relTorusPos[0];
        
        // Parameter für die kleinen Ellipsoide (die Wände des Lochs)
        double wallRadius = segRadius[0] * relTorusr[0]; 
        double gapSize = 0.05 ;//segRadius[0] * relTorusR[0]; // Wie groß ist die Lücke für den Muskel
        double down = segLen[0] * 0.13;
        // 1. Linke Wand
        auto mWallLeft = std::make_shared<SSEllipsoidMesh>(wallRadius, wallRadius, wallRadius*2, "MeshHoleLeft_"+jName, bodyMC, 
             MWMath::Point3D(down, tPosY0, gapSize), // Verschoben in +Z
             MWMath::axisAngle({0,1,0}, -60.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mWallLeft);
        MeshMap["MeshHoleLeft_" + jName] = mWallLeft.get();

        // 2. Rechte Wand
        auto mWallRight = std::make_shared<SSEllipsoidMesh>(wallRadius, wallRadius, wallRadius*2, "MeshHoleRight_"+jName, bodyMC, 
             MWMath::Point3D(down, tPosY0, -gapSize), // Verschoben in -Z
             MWMath::axisAngle({0,1,0}, 60.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mWallRight);
        MeshMap["MeshHoleRight_" + jName] = mWallRight.get();
        
        /* // 3. Boden-Wand (unter dem Muskel, weiter palmar)
        auto mWallBottom = std::make_shared<SSEllipsoidMesh>(wallRadius*2, wallRadius, wallRadius, "MeshHoleBottom_"+jName, bodyMC, 
             MWMath::Point3D(segLen[0] * off[0] + gapSize, tPosY0, 0), // Verschoben in +X (palmar)
             MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mWallBottom);
        MeshMap["MeshHoleBottom_" + jName] = mWallBottom.get(); */

         
        // MCP Joint + Mesh
        MWMath::Point3D jointPosRelMCP = MWMath::Point3D(0, -flength*0.5, 0); 
        jAngle = 0.0; 
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

        // --- ZYLINDER PP ---
        double tr1 = segLen[1] * relTorusR[1];
        double tPosY1 = -segLen[1] *relTorusPos[1]; 
        auto mTorusPP = std::make_shared<SSCylinderMesh>(tr1, tr1*2.0, "MeshTorus_PP"+prefN, bodyPP, 
                MWMath::Point3D(segLen[1] * off[1], tPosY1, 0), 
                MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.3, 0));
        meshes.push_back(mTorusPP);
        MeshMap["MeshTorus_PP" + prefN] = mTorusPP.get();

        // PIP 
        MWMath::Point3D jointPosRelPP = MWMath::Point3D(0, -segLen[1]*0.5, 0); 
        jAngle = MMHA.PIP_M[fidx];
        jName = "PIP" + prefN + "_Joint";
        auto jointPIP = std::make_shared<SSJoint>(jName, jointPosRelPP, MWMath::RotMatrix3x3(), bodyPP, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointPIP);
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

        // --- ZYLINDER MP ---
        double tr2 = segLen[2] * relTorusR[2];
        double tPosY2 = -segLen[2] * relTorusPos[2]; 
        auto mTorusMP = std::make_shared<SSCylinderMesh>(tr2, tr2*2.0, "MeshTorus_MP"+prefN, bodyMP, 
                MWMath::Point3D(segLen[2] * off[2], tPosY2, 0), 
                MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.3, 0.2));
        meshes.push_back(mTorusMP);
        MeshMap["MeshTorus_MP" + prefN] = mTorusMP.get();

        // DIP 
        MWMath::Point3D jointPosRelMP = MWMath::Point3D(0,-segLen[2]*0.5, 0); 
        jAngle = MMHA.DIP_M[fidx];
        jName = "DIP" + prefN + "_Joint";
        auto jointDIP = std::make_shared<SSJoint>(jName, jointPosRelMP, MWMath::RotMatrix3x3(), bodyMP, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointDIP);
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

    // INITIALES UPDATE & DISKRETISIERUNG
    for (auto& m : meshes) { 
        m->InitializeMesh(); 
        m->discretizeMesh(cfg.discretization); 
    }
    rootSystem->update(0);

    // MUSKEL
    int numPoints = cfg.muscleNumPoints.empty() ? 25 : cfg.muscleNumPoints[0];
    
    // Die Namen der Meshes für den Muskelpfad. 
    // ACHTUNG: Der Muskel kann jetzt nicht mehr durch EINEN Torus, er muss die 3 Wände als Hindernis umgehen!
    std::vector<std::string> refBodyNames = { 
        "Mesh_WristFlexorCyl", 
        "Mesh_MC1",
        "MeshHoleLeft_MC1",
        "MeshHoleRight_MC1",
        "Mesh_MCP1_Joint", 
        "MeshTorus_PP1", 
        "Mesh_PIP1_Joint", 
        "MeshTorus_MP1", 
        "Mesh_DIP1_Joint"
    };

    MWMath::Point3D startOffset = MWMath::Point3D(0., 0.5, -0.1);
    MWMath::Point3D endOffset = MWMath::Point3D(0.05, 0.0, 0.0);
    
    SSMuscle* flexor = new SSMuscle("Flexor", numPoints, rootSystem.get(), startOffset, getRefMesh("Mesh_DP1", MeshMap), endOffset);

    for (auto& m : getRefMeshes(refBodyNames, MeshMap)) {
        flexor->meshPtrs.push_back(m);
    }

    flexor->createMusclePointsComplexPath();
    flexor->updateMusclePointsParentsLocal();
    muscles.push_back(flexor);


    std::string returnString = "";
    for (const auto param : processParams) {
        returnString += std::to_string(param) + "\t";
    }
    return returnString;
}

inline std::string buildOHandModelCylEll4Hole(std::vector<std::shared_ptr<SSTissue>>& tissues,std::vector<std::shared_ptr<SSMesh>>& meshes,std::vector<SSMuscle*>& muscles,std::shared_ptr<SSBody>& rootSystem,int numTimeSteps,const SimSettings& cfg,float geometryScaler,const std::vector<double> processParams)
{
    meshes.clear();
    muscles.clear();

    // --- Initialisierung der Hand-Winkel (aus deiner Vorlage) ---
    MWMath::Point3D COLORF1 = MWMath::Point3D(0.50, 0.00, 1.00); 
    MWMath::Point3D COLORF2 = MWMath::Point3D(0.30, 0.30, 1.00); 
    MWMath::Point3D COLORF3 = MWMath::Point3D(0.00, 0.50, 1.00); 
    MWMath::Point3D COLORF4 = MWMath::Point3D(0.00, 0.90, 1.00);
    Hand::MMHandAngles MMHA;

    // --- PARAMETER ---
    double HL = 1.8 * geometryScaler;
    double HB = 0.8 * geometryScaler;
    double scale = HL / Hand::RefHandLength;
    double Sign = 1.0; // 1.0 = Rechts
    float GS = geometryScaler;
    
    std::unordered_map<std::string, SSMesh*> MeshMap;

    std::vector<double> width = { 0.15*0.5 * GS, 0.1*0.5 * GS, 0.07*0.5 * GS, 0.05*0.5 * GS };
    
    double rWF = 1.0; 
    double PF = processParams.size() > 0 ? processParams[0] : 1.0;
    double offF = processParams.size() > 1 ? processParams[1] : 1.0;
    double RF = processParams.size() > 2 ? processParams[2] : 1.0;
    double relT = processParams.size() > 3 ? processParams[3] : 0.7;
    
    std::vector<double> relTorusPos = {0.35*PF, 0.1*PF, 0.0*PF}; 
    std::vector<double> relTorusR = {0.1*RF, 0.19*RF, 0.28*RF}; 
    std::vector<double> relTorusr = {relT, relT, relT};
    std::vector<double> off = {0.06*1.0, 0.07*4.0, 0.18*2.5}; 

    double segLen[4];
    double segRadius[4];
    double jointRadius[3];
    std::vector<std::string> fingerNames = {"Thumb", "Index", "Middle", "Ring", "Little"};

    // 1. ROOT
    MWMath::RotMatrix3x3 Ausrichtung = MWMath::axisAngle({0,0,1}, 90.0) * MWMath::axisAngle({0,1,0}, 180.0);
    rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0, 0, 0), Ausrichtung, nullptr);
    tissues.push_back(rootSystem);

    MWMath::Point3D wristOffset(0.0, 0.003567 * scale, -0.003901 * scale * Sign);
    
    auto wristJointAbd = std::make_shared<SSJoint>("Wrist_JointAbd", wristOffset, MWMath::RotMatrix3x3(), rootSystem, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointAbd);
    auto jMeshWrist = std::make_shared<SSEllipsoidMesh>(0.01 * scale, 0.01 * scale, 0.01 * scale, "Mesh_Wrist_Joint", wristJointAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0,0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMeshWrist);
    MeshMap["Mesh_Wrist_Joint"] = jMeshWrist.get();

    auto wristJointFlex = std::make_shared<SSJoint>("Wrist_JointFlex", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointAbd, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointFlex);

    // 2. CARPALS
    auto carpals = std::make_shared<SSBody>("Carpals", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointFlex);
    tissues.push_back(carpals);
    
    double carpalThickness = width[0] * 1.5; 
    double carpalLength = width[0] * 2.5;    
    double carpalWidth = width[0] * 3.5;     
    
    auto meshCarpals = std::make_shared<SSEllipsoidMesh>(carpalThickness, carpalLength, carpalWidth, "Mesh_Carpals", carpals, MWMath::Point3D(0, 0.0*scale, 0), MWMath::RotMatrix3x3(), MWMath::Point3D(0.5, 0.5, 0.5));
    meshes.push_back(meshCarpals);
    MeshMap["Mesh_Carpals"] = meshCarpals.get();

    // 3. WRAPPING SURFACE 
    double flexCylRadius = 0.0115 * scale;
    double flexCylLength = 0.1400 * scale; 
    MWMath::Point3D flexCylOffset(carpalThickness * 0.5 + flexCylRadius, 0.02 * scale, 0.0);

    auto flexorCylMesh = std::make_shared<SSCylinderMesh>(flexCylRadius, flexCylLength, "Mesh_WristFlexorCyl", carpals, flexCylOffset, MWMath::axisAngle({ 1, 0, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));  
    meshes.push_back(flexorCylMesh);
    MeshMap["Mesh_WristFlexorCyl"] = flexorCylMesh.get();

    // FINGER SCHLEIFE
    for (int i = 1; i < 2; i++) {
        int fidx = i;
        double jAngle;
        std::string prefN = std::to_string(fidx);

        segLen[0] = Hand::SegLenRatios[fidx][0] * HL;
        segLen[1] = Hand::SegLenRatios[fidx][1] * HL;
        segLen[2] = Hand::SegLenRatios[fidx][2] * HL;
        segLen[3] = Hand::SegLenRatios[fidx][3] * HL;

        segRadius[0] = Hand::FingerWrapData[fidx-1][0].radius * scale;
        segRadius[1] = Hand::FingerWrapData[fidx-1][1].radius * scale;
        segRadius[2] = Hand::FingerWrapData[fidx-1][2].radius * scale;
        segRadius[3] = Hand::FingerWrapData[fidx-1][3].radius * scale;

        jointRadius[0] = segRadius[0] * 1.1; 
        jointRadius[1] = segRadius[2] * 1.1; 
        jointRadius[2] = segRadius[3] * 1.1; 

        MWMath::Point3D posMCP(0.0, Hand::JointPosTable[fidx][0] * HL, Hand::JointPosTable[fidx][1] * HB * Sign);
        
        // CMC
        MWMath::Point3D posCMC = Hand::CMCOffsets[fidx - 1] * scale;
        posCMC.z *= Sign;
        MWMath::RotMatrix3x3 forientation = buildOrientation(posMCP, posCMC);
        jAngle = 0.0; 
        std::string jName = "CMC" + prefN + "_JointAbd";
        auto jointCMCAbd = std::make_shared<SSJoint>(jName, posCMC, forientation, carpals, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointCMCAbd);
        auto jMeshCMC = std::make_shared<SSEllipsoidMesh>(width[0], width[0], width[0], "Mesh_CMC"+prefN+"_Joint", jointCMCAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0 , 0 }, 90.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMeshCMC);
        MeshMap["Mesh_CMC"+prefN+"_Joint"] = jMeshCMC.get();

        jAngle = 0.0;
        jName = "CMC" + prefN + "_JointFlex";
        auto jointCMCFlex = std::make_shared<SSJoint>(jName, MWMath::Point3D(0.,0.,0.), MWMath::RotMatrix3x3(), jointCMCAbd, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointCMCFlex);

        // ==========================================
        // MC SEGMENT
        // ==========================================
        double flength = MWMath::distance(posCMC, posMCP);
        jName = "MC" + prefN;
        auto bodyMC = std::make_shared<SSBody>(jName, MWMath::Point3D(0.,-flength*0.5,0.), MWMath::RotMatrix3x3(), jointCMCFlex);
        tissues.push_back(bodyMC);
        auto meshMC = std::make_shared<SSEllipsoidMesh>(segRadius[0], segRadius[0]*1.4, flength*0.5, "Mesh_"+jName, bodyMC, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1, 0, 0 }, 90.0), COLORF1);
        meshes.push_back(meshMC);
        MeshMap["Mesh_" + jName] = meshMC.get();

        // --- NEU: 3 ELLIPSOIDE ALS LOCH (Ersetzt den Torus am MC) ---
        // Basisposition für das "Loch"
        double tPosY0 = -flength * relTorusPos[0];
        
        // Parameter für die kleinen Ellipsoide (die Wände des Lochs)
        double wallRadius = segRadius[0] * relTorusr[0]; 
        double gapSize = 0.1 ;//segRadius[0] * relTorusR[0]; // Wie groß ist die Lücke für den Muskel
        double down = segLen[0] * 0.13;
        // 1. Linke Wand
        auto mWallLeft = std::make_shared<SSEllipsoidMesh>(wallRadius*1.5, wallRadius*1.5, wallRadius*4, "MeshHoleLeft_"+jName, bodyMC, 
             MWMath::Point3D(down, tPosY0, gapSize), // Verschoben in +Z
             MWMath::axisAngle({0,1,0}, -90.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mWallLeft);
        MeshMap["MeshHoleLeft_" + jName] = mWallLeft.get();

        // 2. Rechte Wand
        auto mWallRight = std::make_shared<SSEllipsoidMesh>(wallRadius*1.5, wallRadius*1.5, wallRadius*4, "MeshHoleRight_"+jName, bodyMC, 
             MWMath::Point3D(down, tPosY0, -gapSize), // Verschoben in -Z
             MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mWallRight);
        MeshMap["MeshHoleRight_" + jName] = mWallRight.get();
        
        // 3. Boden-Wand (unter dem Muskel, weiter palmar)
        auto mWallBottom = std::make_shared<SSEllipsoidMesh>(wallRadius*2, wallRadius*2, wallRadius*6, "MeshHoleBottom_"+jName, bodyMC, 
             MWMath::Point3D(down*1.7, tPosY0, 0), // Verschoben in +X (palmar)
             MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mWallBottom);
        MeshMap["MeshHoleBottom_" + jName] = mWallBottom.get();

         
        // MCP Joint + Mesh
        MWMath::Point3D jointPosRelMCP = MWMath::Point3D(0, -flength*0.5, 0); 
        jAngle = 0.0; 
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

        // --- ZYLINDER PP ---
        double tr1 = segLen[1] * relTorusR[1];
        double tPosY1 = -segLen[1] *relTorusPos[1]; 
        auto mTorusPP = std::make_shared<SSCylinderMesh>(tr1, tr1*2.0, "MeshTorus_PP"+prefN, bodyPP, 
                MWMath::Point3D(segLen[1] * off[1], tPosY1, 0), 
                MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.3, 0));
        meshes.push_back(mTorusPP);
        MeshMap["MeshTorus_PP" + prefN] = mTorusPP.get();

        // PIP 
        MWMath::Point3D jointPosRelPP = MWMath::Point3D(0, -segLen[1]*0.5, 0); 
        jAngle = MMHA.PIP_M[fidx];
        jName = "PIP" + prefN + "_Joint";
        auto jointPIP = std::make_shared<SSJoint>(jName, jointPosRelPP, MWMath::RotMatrix3x3(), bodyPP, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointPIP);
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

        // --- ZYLINDER MP ---
        double tr2 = segLen[2] * relTorusR[2];
        double tPosY2 = -segLen[2] * relTorusPos[2]; 
        auto mTorusMP = std::make_shared<SSCylinderMesh>(tr2, tr2*2.0, "MeshTorus_MP"+prefN, bodyMP, 
                MWMath::Point3D(segLen[2] * off[2], tPosY2, 0), 
                MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.3, 0.2));
        meshes.push_back(mTorusMP);
        MeshMap["MeshTorus_MP" + prefN] = mTorusMP.get();

        // DIP 
        MWMath::Point3D jointPosRelMP = MWMath::Point3D(0,-segLen[2]*0.5, 0); 
        jAngle = MMHA.DIP_M[fidx];
        jName = "DIP" + prefN + "_Joint";
        auto jointDIP = std::make_shared<SSJoint>(jName, jointPosRelMP, MWMath::RotMatrix3x3(), bodyMP, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointDIP);
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

    // INITIALES UPDATE & DISKRETISIERUNG
    for (auto& m : meshes) { 
        m->InitializeMesh(); 
        m->discretizeMesh(cfg.discretization); 
    }
    rootSystem->update(0);

    // MUSKEL
    int numPoints = cfg.muscleNumPoints.empty() ? 25 : cfg.muscleNumPoints[0];
    
    // Die Namen der Meshes für den Muskelpfad. 
    // ACHTUNG: Der Muskel kann jetzt nicht mehr durch EINEN Torus, er muss die 3 Wände als Hindernis umgehen!
    std::vector<std::string> refBodyNames = { 
        "Mesh_WristFlexorCyl", 
        "Mesh_MC1",
        "MeshHoleBottom_MC1",
        "MeshHoleLeft_MC1",
        "MeshHoleRight_MC1",
        "Mesh_MCP1_Joint", 
        "MeshTorus_PP1", 
        "Mesh_PIP1_Joint", 
        "MeshTorus_MP1", 
        "Mesh_DIP1_Joint"
    };

    MWMath::Point3D startOffset = MWMath::Point3D(0., 0.5, -0.1);
    MWMath::Point3D endOffset = MWMath::Point3D(0.05, 0.0, 0.0);
    
    SSMuscle* flexor = new SSMuscle("Flexor", numPoints, rootSystem.get(), startOffset, getRefMesh("Mesh_DP1", MeshMap), endOffset);

    for (auto& m : getRefMeshes(refBodyNames, MeshMap)) {
        flexor->meshPtrs.push_back(m);
    }

    flexor->createMusclePointsComplexPath();
    flexor->updateMusclePointsParentsLocal();
    muscles.push_back(flexor);


    std::string returnString = "";
    for (const auto param : processParams) {
        returnString += std::to_string(param) + "\t";
    }
    return returnString;
}

inline std::string buildOHandModelCylCyl4Hole(std::vector<std::shared_ptr<SSTissue>>& tissues,std::vector<std::shared_ptr<SSMesh>>& meshes,std::vector<SSMuscle*>& muscles,std::shared_ptr<SSBody>& rootSystem,int numTimeSteps,const SimSettings& cfg,float geometryScaler,const std::vector<double> processParams)
{
    meshes.clear();
    muscles.clear();

    // --- Initialisierung der Hand-Winkel (aus deiner Vorlage) ---
    MWMath::Point3D COLORF1 = MWMath::Point3D(0.50, 0.00, 1.00); 
    MWMath::Point3D COLORF2 = MWMath::Point3D(0.30, 0.30, 1.00); 
    MWMath::Point3D COLORF3 = MWMath::Point3D(0.00, 0.50, 1.00); 
    MWMath::Point3D COLORF4 = MWMath::Point3D(0.00, 0.90, 1.00);
    Hand::MMHandAngles MMHA;

    // --- PARAMETER ---
    double HL = 1.8 * geometryScaler;
    double HB = 0.8 * geometryScaler;
    double scale = HL / Hand::RefHandLength;
    double Sign = 1.0; // 1.0 = Rechts
    float GS = geometryScaler;
    
    std::unordered_map<std::string, SSMesh*> MeshMap;

    std::vector<double> width = { 0.15*0.5 * GS, 0.1*0.5 * GS, 0.07*0.5 * GS, 0.05*0.5 * GS };
    
    double rWF = 1.0; 
    double PF = processParams.size() > 0 ? processParams[0] : 1.0;
    double offF = processParams.size() > 1 ? processParams[1] : 1.0;
    double RF = processParams.size() > 2 ? processParams[2] : 1.0;
    double relT = processParams.size() > 3 ? processParams[3] : 0.7;
    
    std::vector<double> relTorusPos = {0.25*PF, 0.1*PF, 0.0*PF}; 
    std::vector<double> relTorusR = {0.1*RF, 0.19*RF, 0.28*RF}; 
    std::vector<double> relTorusr = {relT, relT, relT};
    std::vector<double> off = {0.06*1.0, 0.07*4.0, 0.18*2.5}; 

    double segLen[4];
    double segRadius[4];
    double jointRadius[3];
    std::vector<std::string> fingerNames = {"Thumb", "Index", "Middle", "Ring", "Little"};

    // 1. ROOT
    MWMath::RotMatrix3x3 Ausrichtung = MWMath::axisAngle({0,0,1}, 90.0) * MWMath::axisAngle({0,1,0}, 180.0);
    rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0, 0, 0), Ausrichtung, nullptr);
    tissues.push_back(rootSystem);

    MWMath::Point3D wristOffset(0.0, 0.003567 * scale, -0.003901 * scale * Sign);
    
    auto wristJointAbd = std::make_shared<SSJoint>("Wrist_JointAbd", wristOffset, MWMath::RotMatrix3x3(), rootSystem, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointAbd);
    auto jMeshWrist = std::make_shared<SSEllipsoidMesh>(0.01 * scale, 0.01 * scale, 0.01 * scale, "Mesh_Wrist_Joint", wristJointAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0,0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMeshWrist);
    MeshMap["Mesh_Wrist_Joint"] = jMeshWrist.get();

    auto wristJointFlex = std::make_shared<SSJoint>("Wrist_JointFlex", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointAbd, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointFlex);

    // 2. CARPALS
    auto carpals = std::make_shared<SSBody>("Carpals", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointFlex);
    tissues.push_back(carpals);
    
    double carpalThickness = width[0] * 1.5; 
    double carpalLength = width[0] * 2.5;    
    double carpalWidth = width[0] * 3.5;     
    
    auto meshCarpals = std::make_shared<SSEllipsoidMesh>(carpalThickness, carpalLength, carpalWidth, "Mesh_Carpals", carpals, MWMath::Point3D(0, 0.0*scale, 0), MWMath::RotMatrix3x3(), MWMath::Point3D(0.5, 0.5, 0.5));
    meshes.push_back(meshCarpals);
    MeshMap["Mesh_Carpals"] = meshCarpals.get();

    // 3. WRAPPING SURFACE 
    double flexCylRadius = 0.0115 * scale;
    double flexCylLength = 0.1400 * scale; 
    MWMath::Point3D flexCylOffset(carpalThickness * 0.5 + flexCylRadius, 0.02 * scale, 0.0);

    auto flexorCylMesh = std::make_shared<SSCylinderMesh>(flexCylRadius, flexCylLength, "Mesh_WristFlexorCyl", carpals, flexCylOffset, MWMath::axisAngle({ 1, 0, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));  
    meshes.push_back(flexorCylMesh);
    MeshMap["Mesh_WristFlexorCyl"] = flexorCylMesh.get();

    // FINGER SCHLEIFE
    for (int i = 1; i < 2; i++) {
        int fidx = i;
        double jAngle;
        std::string prefN = std::to_string(fidx);

        segLen[0] = Hand::SegLenRatios[fidx][0] * HL;
        segLen[1] = Hand::SegLenRatios[fidx][1] * HL;
        segLen[2] = Hand::SegLenRatios[fidx][2] * HL;
        segLen[3] = Hand::SegLenRatios[fidx][3] * HL;

        segRadius[0] = Hand::FingerWrapData[fidx-1][0].radius * scale;
        segRadius[1] = Hand::FingerWrapData[fidx-1][1].radius * scale;
        segRadius[2] = Hand::FingerWrapData[fidx-1][2].radius * scale;
        segRadius[3] = Hand::FingerWrapData[fidx-1][3].radius * scale;

        jointRadius[0] = segRadius[0] * 1.1; 
        jointRadius[1] = segRadius[2] * 1.1; 
        jointRadius[2] = segRadius[3] * 1.1; 

        MWMath::Point3D posMCP(0.0, Hand::JointPosTable[fidx][0] * HL, Hand::JointPosTable[fidx][1] * HB * Sign);
        
        // CMC
        MWMath::Point3D posCMC = Hand::CMCOffsets[fidx - 1] * scale;
        posCMC.z *= Sign;
        MWMath::RotMatrix3x3 forientation = buildOrientation(posMCP, posCMC);
        jAngle = 0.0; 
        std::string jName = "CMC" + prefN + "_JointAbd";
        auto jointCMCAbd = std::make_shared<SSJoint>(jName, posCMC, forientation, carpals, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointCMCAbd);
        auto jMeshCMC = std::make_shared<SSEllipsoidMesh>(width[0], width[0], width[0], "Mesh_CMC"+prefN+"_Joint", jointCMCAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0 , 0 }, 90.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMeshCMC);
        MeshMap["Mesh_CMC"+prefN+"_Joint"] = jMeshCMC.get();

        jAngle = 0.0;
        jName = "CMC" + prefN + "_JointFlex";
        auto jointCMCFlex = std::make_shared<SSJoint>(jName, MWMath::Point3D(0.,0.,0.), MWMath::RotMatrix3x3(), jointCMCAbd, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointCMCFlex);

        // ==========================================
        // MC SEGMENT
        // ==========================================
        double flength = MWMath::distance(posCMC, posMCP);
        jName = "MC" + prefN;
        auto bodyMC = std::make_shared<SSBody>(jName, MWMath::Point3D(0.,-flength*0.5,0.), MWMath::RotMatrix3x3(), jointCMCFlex);
        tissues.push_back(bodyMC);
        auto meshMC = std::make_shared<SSEllipsoidMesh>(segRadius[0], segRadius[0]*1.4, flength*0.5, "Mesh_"+jName, bodyMC, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1, 0, 0 }, 90.0), COLORF1);
        meshes.push_back(meshMC);
        MeshMap["Mesh_" + jName] = meshMC.get();

        // --- NEU: 3 ELLIPSOIDE ALS LOCH (Ersetzt den Torus am MC) ---
        // Basisposition für das "Loch"
        double tPosY0 = -flength * relTorusPos[0];
        
        // Parameter für die kleinen Ellipsoide (die Wände des Lochs)
        double wallRadius = segRadius[0] * relTorusr[0]; 
        double gapSize = 0.1 ;//segRadius[0] * relTorusR[0]; // Wie groß ist die Lücke für den Muskel
        double down = segLen[0] * 0.13;
        // 1. Linke Wand
        auto mWallLeft = std::make_shared<SSCylinderMesh>(wallRadius*1.7, wallRadius*4, "MeshHoleLeft_"+jName, bodyMC, 
             MWMath::Point3D(down, tPosY0, gapSize), // Verschoben in +Z
             MWMath::axisAngle({0,1,0}, -90.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mWallLeft);
        MeshMap["MeshHoleLeft_" + jName] = mWallLeft.get();

        // 2. Rechte Wand
        auto mWallRight = std::make_shared<SSCylinderMesh>(wallRadius*1.7, wallRadius*4, "MeshHoleRight_"+jName, bodyMC, 
             MWMath::Point3D(down, tPosY0, -gapSize), // Verschoben in -Z
             MWMath::axisAngle({0,1,0}, 90.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mWallRight);
        MeshMap["MeshHoleRight_" + jName] = mWallRight.get();
        
        // 3. Boden-Wand (unter dem Muskel, weiter palmar)
        auto mWallBottom = std::make_shared<SSCylinderMesh>(wallRadius*1.8, wallRadius*6, "MeshHoleBottom_"+jName, bodyMC, 
             MWMath::Point3D(down*1.5, tPosY0, 0), // Verschoben in +X (palmar)
             MWMath::axisAngle({0,1,0}, 0.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mWallBottom);
        MeshMap["MeshHoleBottom_" + jName] = mWallBottom.get();

         
        // MCP Joint + Mesh
        MWMath::Point3D jointPosRelMCP = MWMath::Point3D(0, -flength*0.5, 0); 
        jAngle = 0.0; 
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

        // --- ZYLINDER PP ---
        double tr1 = segLen[1] * relTorusR[1];
        double tPosY1 = -segLen[1] *relTorusPos[1]; 
        auto mTorusPP = std::make_shared<SSCylinderMesh>(tr1, tr1*2.0, "MeshTorus_PP"+prefN, bodyPP, 
                MWMath::Point3D(segLen[1] * off[1], tPosY1, 0), 
                MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.3, 0));
        meshes.push_back(mTorusPP);
        MeshMap["MeshTorus_PP" + prefN] = mTorusPP.get();

        // PIP 
        MWMath::Point3D jointPosRelPP = MWMath::Point3D(0, -segLen[1]*0.5, 0); 
        jAngle = MMHA.PIP_M[fidx];
        jName = "PIP" + prefN + "_Joint";
        auto jointPIP = std::make_shared<SSJoint>(jName, jointPosRelPP, MWMath::RotMatrix3x3(), bodyPP, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointPIP);
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

        // --- ZYLINDER MP ---
        double tr2 = segLen[2] * relTorusR[2];
        double tPosY2 = -segLen[2] * relTorusPos[2]; 
        auto mTorusMP = std::make_shared<SSCylinderMesh>(tr2, tr2*2.0, "MeshTorus_MP"+prefN, bodyMP, 
                MWMath::Point3D(segLen[2] * off[2], tPosY2, 0), 
                MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.3, 0.2));
        meshes.push_back(mTorusMP);
        MeshMap["MeshTorus_MP" + prefN] = mTorusMP.get();

        // DIP 
        MWMath::Point3D jointPosRelMP = MWMath::Point3D(0,-segLen[2]*0.5, 0); 
        jAngle = MMHA.DIP_M[fidx];
        jName = "DIP" + prefN + "_Joint";
        auto jointDIP = std::make_shared<SSJoint>(jName, jointPosRelMP, MWMath::RotMatrix3x3(), bodyMP, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointDIP);
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

    // INITIALES UPDATE & DISKRETISIERUNG
    for (auto& m : meshes) { 
        m->InitializeMesh(); 
        m->discretizeMesh(cfg.discretization); 
    }
    rootSystem->update(0);

    // MUSKEL
    int numPoints = cfg.muscleNumPoints.empty() ? 25 : cfg.muscleNumPoints[0];
    
    // Die Namen der Meshes für den Muskelpfad. 
    // ACHTUNG: Der Muskel kann jetzt nicht mehr durch EINEN Torus, er muss die 3 Wände als Hindernis umgehen!
    std::vector<std::string> refBodyNames = { 
        "Mesh_WristFlexorCyl", 
        "Mesh_MC1",
        "MeshHoleBottom_MC1",
        "MeshHoleLeft_MC1",
        "MeshHoleRight_MC1",
        "Mesh_MCP1_Joint", 
        "MeshTorus_PP1", 
        "Mesh_PIP1_Joint", 
        "MeshTorus_MP1", 
        "Mesh_DIP1_Joint"
    };

    MWMath::Point3D startOffset = MWMath::Point3D(0., 0.5, -0.1);
    MWMath::Point3D endOffset = MWMath::Point3D(0.05, 0.0, 0.0);
    
    SSMuscle* flexor = new SSMuscle("Flexor", numPoints, rootSystem.get(), startOffset, getRefMesh("Mesh_DP1", MeshMap), endOffset);

    for (auto& m : getRefMeshes(refBodyNames, MeshMap)) {
        flexor->meshPtrs.push_back(m);
    }

    flexor->createMusclePointsComplexPath();
    flexor->updateMusclePointsParentsLocal();
    muscles.push_back(flexor);


    std::string returnString = "";
    for (const auto param : processParams) {
        returnString += std::to_string(param) + "\t";
    }
    return returnString;
}


