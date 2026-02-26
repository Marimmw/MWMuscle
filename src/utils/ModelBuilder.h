#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>

#include "utils/MWMath.h"

#include "simpleSimulation/SSMuscle.h"
#include "simpleSimulation/SSMeshes.h"
#include "simpleSimulation/SSBody.h"

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
