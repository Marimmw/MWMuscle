#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <cmath>

#include "utils/MWMath.h"
#include "utils/utility.h"

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

struct ViaPointDef
{
    std::string name;
    std::string refObj;
    MWMath::Point3D relPos;
    double tolerance = 0.01;
};

struct MuscleDef{
    std::string name;
    std::string originMeshName;
    MWMath::Point3D originRelPos;
    std::string insertionMeshName;
    MWMath::Point3D insertionRelPos;
    std::vector<std::string> meshes;
    std::vector<ViaPointDef> viaPoints;
    int numNodes = 2;
    MWMath::Point3D mcolor;
    double torusPathDirection = 1.0f;
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
            std::cerr << "     [getRefMeshes-WARNING] Didnt find Mesh name '" << name << "'!" << std::endl;
        }
    }
    return meshPtrs;
}

inline SSMesh* getRefMesh(std::string meshName, std::unordered_map<std::string, SSMesh*>& MeshMap) {
    if (MeshMap.find(meshName) != MeshMap.end()) {
        return MeshMap[meshName];
    }
    else {
        std::cerr << "     [getRefMesh-WARNING] Didnt find Mesh name '" << meshName << "' nicht gefunden!" << std::endl;
        return nullptr;
    }
}

inline SSTissue* getRefTissue(std::string tissueName, std::vector<std::shared_ptr<SSTissue>>& tissues) {
    for (const auto& tissue : tissues) {
        if (tissue->Name == tissueName) {
            return tissue.get();
        }
    }
    std::cerr << "Warnung: Tissue mit Namen '" << tissueName << "' nicht gefunden!" << std::endl;
    return nullptr;
}

inline SSMesh* getRefMeshFromList(std::string meshName, std::vector<std::shared_ptr<SSMesh>>& meshes) {
    for (const auto& mesh : meshes) {
        if (mesh->Name == meshName) {
            return mesh.get();
        }
    }
    std::cerr << "     [getRefMeshFromList-WARNING] Didnt find Mesh name '" << meshName << "' nicht gefunden!" << std::endl;
    return nullptr;
}

///// HAND MODEL BUILDING /////

namespace Hand {
    constexpr double RefHandLength = 0.195;
    constexpr double ThumbRotX = -60.0; // Grad
    constexpr double ThumbRotY =  44.1; // Grad


    const MWMath::Point3D FLEXORCOLOR(0.7,0.,0.);
    const MWMath::Point3D EXTENSORCOLOR(0.,0.,0.7);
    const MWMath::Point3D PINTEROSSEICCOLOR(0.6,0.3,0.);
    const MWMath::Point3D DINTEROSSEICCOLOR(0.6,0.6,0.);
    const MWMath::Point3D LUMBRICALCOLOR(0.,0.7,0.);
    const MWMath::Point3D OTHERMUSCLECOLOR(0.4,0.2,0.);

    inline std::vector<MuscleDef> muscleDefsAny = {
        // testmuscle
        {
            "TestMuscle", 
            "", MWMath::Point3D{1.8*0.1,1.8*0.3, 0.0}, 
            "Mesh_DP1", MWMath::Point3D{0.05*1.1, -1.8*0.097, 0.0}, 
            {"Mesh_CarpalTunnel","MeshVP_MC1","Mesh_MC1_Joint","Mesh_PP1","MeshVP_PP1","Mesh_PIP1_Joint","Mesh_MP1","MeshVP_MP1","Mesh_DIP1_Joint","Mesh_DP1"}, 
            {}, 
            2, MWMath::Point3D{1.0, 0.0, 0.0},
        },
        
        // ==========================================
        // --- DAUMEN (Finger 1) ---
        // ==========================================
        {
            "AbductorPollicisLongus", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_MC0", MWMath::Point3D{-0.01, 0.007, 0.0}, 
            {}, {}, 
            2, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "ExtensorPollicisBrevis", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_MC0", MWMath::Point3D{-0.01, -0.02, 0.0} , 
            {}, 
            { {"Via1_EPB", "Mesh_MC0", MWMath::Point3D{-0.01, 0.0, 0.0} } }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "ExtensorPollicisLongus", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_MC0", MWMath::Point3D{-0.01, -0.02, 0.0} , 
            {}, 
            { {"Via1_EPL", "Mesh_MC0", MWMath::Point3D{-0.01, 0.0, 0.0} } }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorPollicisLongus", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_MC0", MWMath::Point3D{0.008, -0.015, 0.0} , 
            {}, 
            { {"Via_FPL", "Mesh_MC0", MWMath::Point3D{0.008, 0.0, 0.0} } }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        /* 
        // ==========================================
        // --- ZEIGEFINGER (Finger 2) ---
        // ==========================================
        {
            "ExtensorIndicis", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal2", MWMath::Point3D{-0.01, -0.02, 0.0} , 
            {}, 
            { {"Via1_EI", "Mesh_Metacarpal2", MWMath::Point3D{-0.007, 0.02, 0.0} } }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "ExtensorCarpiRadialisLongus", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal2", MWMath::Point3D{-0.006, -0.01, 0.0} , 
            {}, {}, 
            2, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorCarpiRadialis", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal2", MWMath::Point3D{0.006, 0.03, 0.0} , 
            {}, {}, 
            2, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorDigitorumSuperficialis2", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal2", MWMath::Point3D{0.007, -0.025, 0.0} , 
            {}, 
            { {"Via_FDS2", "Mesh_Metacarpal2", MWMath::Point3D{0.008, 0.01, 0.0} } }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorDigitorumProfundus2", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal2", MWMath::Point3D{0.009, -0.02, 0.0} , 
            {}, 
            { {"Via_FDP2", "Mesh_Metacarpal2", MWMath::Point3D{0.009, 0.01, 0.0} } }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "ExtensorDigitorum2", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal2", MWMath::Point3D{-0.009, -0.02, 0.0} , 
            {}, 
            { {"Via_ED2", "Mesh_Metacarpal2", MWMath::Point3D{-0.009, 0.01, 0.0} } }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },
         */
        /* 
        // ==========================================
        // --- MITTELFINGER (Finger 3) ---
        // ==========================================
        {
            "ExtensorCarpiRadialisBrevis", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal3", MWMath::Point3D{-0.006, 0.025, 0.0} , 
            {}, 
            { {"Via_ECRB", "Mesh_Metacarpal3", MWMath::Point3D{-0.006, 0.027, 0.0} } }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorDigitorumSuperficialis3", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal3", MWMath::Point3D{0.007, -0.025, 0.0} , 
            {}, 
            { {"Via_FDS3", "Mesh_Metacarpal3", MWMath::Point3D{0.008, 0.01, 0.0} } }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorDigitorumProfundus3", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal3", MWMath::Point3D{0.009, -0.02, 0.0} , 
            {}, 
            { {"Via_FDP3", "Mesh_Metacarpal3", MWMath::Point3D{0.009, 0.01, 0.0} } }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "ExtensorDigitorum3", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal3", MWMath::Point3D{-0.009, -0.02, 0.0} , 
            {}, 
            { {"Via_ED3", "Mesh_Metacarpal3", MWMath::Point3D{-0.009, 0.01, 0.0} } }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },

        // ==========================================
        // --- RINGFINGER (Finger 4) ---
        // ==========================================
        {
            "PalmarisLongus", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal4", MWMath::Point3D{0.004, 0.02, 0.0} , 
            {}, {}, 
            2, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorDigitorumSuperficialis4", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal4", MWMath::Point3D{0.007, -0.025, 0.0} , 
            {}, 
            { {"Via_FDS4", "Mesh_Metacarpal4", MWMath::Point3D{0.008, 0.01, 0.0} } }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorDigitorumProfundus4", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal4", MWMath::Point3D{0.009, -0.02, 0.0} , 
            {}, 
            { {"Via_FDP4", "Mesh_Metacarpal4", MWMath::Point3D{0.009, 0.01, 0.0} } }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "ExtensorDigitorum4", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal4", MWMath::Point3D{-0.009, -0.02, 0.0} , 
            {}, 
            { {"Via_ED4", "Mesh_Metacarpal4", MWMath::Point3D{-0.009, 0.01, 0.0} } }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },

        // ==========================================
        // --- KLEINER FINGER (Finger 5) ---
        // ==========================================
        {
            "ExtensorCarpiUlnaris", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal5", MWMath::Point3D{-0.002, 0.02, -0.007} , 
            {}, {}, 
            2, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorCarpiUlnaris", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal5", MWMath::Point3D{0.009, 0.02, 0.007} , 
            {}, {}, 
            2, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorDigitorumSuperficialis5", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal5", MWMath::Point3D{0.009, -0.02, 0.0} , 
            {}, 
            { {"Via_FDS5", "Mesh_Metacarpal5", MWMath::Point3D{0.009, 0.01, 0.0} } }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorDigitorumProfundus5", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal5", MWMath::Point3D{0.009, -0.02, 0.0} , 
            {}, 
            { {"Via_FDP5", "Mesh_Metacarpal5", MWMath::Point3D{0.009, 0.01, 0.0} } }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "ExtensorDigitorum5", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal5", MWMath::Point3D{-0.009, -0.02, 0.0} , 
            {}, 
            { {"Via_ED5", "Mesh_Metacarpal5", MWMath::Point3D{-0.009, 0.01, 0.0} } }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "ExtensorDigitiMinimi", 
            "Mesh_Forearm_Placeholder", MWMath::Point3D{0.0, 0.0, 0.0}, 
            "Mesh_Metacarpal5", MWMath::Point3D{-0.009, -0.02, 0.0} , 
            {}, 
            { {"Via_EDM", "Mesh_Metacarpal5", MWMath::Point3D{-0.009, 0.01, 0.0} } }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        } */
    };

    inline MWMath::RotMatrix3x3 RotMiraToWorld = MWMath::RotMatrix3x3( 1, 0, 0, -1, 0, 0, 0, -1, 0);
    inline MWMath::RotMatrix3x3 RotMiraToMe = MWMath::RotMatrix3x3( 0, 0, 1, 0, 1, 0, -1, 0, 0);
    inline const MWMath::Point3D WristToUlnaStyloidMM = MWMath::Point3D(5.0, 30.0, 25.0);; // MWMath::Point3D(5.0, 30.0, 25.0); // selbst geschätzt // RefSys = carpals
    inline std::vector<MuscleDef> muscleDefsMira = {
        
        // ==========================================
        // --- UNTERARM / HANDGELENK (WRIST) ---
        // ==========================================
        {
            "ExtensorCarpiRadialisLongus", 
            "Mesh_Humerus", MWMath::Point3D{0.0194, 0.2427, 0.0136}  , // Op
            "Mesh_MC1", MWMath::Point3D{0.0505, -0.0368, 0.0028}  ,    // I
            {}, 
            { 
                {"Via_Od", "Mesh_Humerus", MWMath::Point3D{0.0313, 0.2135, 0.0125}  },
                {"Via_Rad", "Mesh_Radius", MWMath::Point3D{0.0439, 0.0074, 0.0064}  }
            }, 
            4, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "ExtensorCarpiRadialisBrevis", 
            "Mesh_Humerus", MWMath::Point3D{0.0306, 0.2072, 0.0172}  , 
            "Mesh_MC2", MWMath::Point3D{0.0351, -0.0392, 0.0068}  , 
            {}, 
            {
                {"Via_Od", "Mesh_Radius", MWMath::Point3D{0.0311, 0.1342, 0.0141}  },
                {"Via_Rad", "Mesh_Radius", MWMath::Point3D{0.0399, 0.0081, 0.0080}  }
            }, 
            4, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "ExtensorCarpiUlnaris", 
            "Mesh_Humerus", MWMath::Point3D{0.0176, 0.2034, 0.0217}  , 
            "Mesh_MC4", MWMath::Point3D{0.0001, -0.0384, -0.0076}  , 
            {}, 
            {
                {"Via_Od", "Mesh_Ulna", MWMath::Point3D{0.0048, 0.0929, 0.0124}  },
                {"Via_Rad", "Mesh_Radius", MWMath::Point3D{0.0078, 0.0014, -0.0004}  }
            }, 
            4, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "PalmarisLongus", 
            "Mesh_Humerus", MWMath::Point3D{-0.0056, 0.2017, -0.0201}  , 
            "Mesh_MC2", MWMath::Point3D{0.0363, -0.0823, -0.0103}  , 
            {}, 
            {
                {"Via_Od", "Mesh_Ulna", MWMath::Point3D{0.0332, 0.2259, -0.0017}  },
                {"Via_Rad", "Mesh_Radius", MWMath::Point3D{0.0275, 0.0005, -0.0266}  },
                {"Via_MC", "Mesh_MC2", MWMath::Point3D{0.0317, -0.0225, -0.0245}  }
            }, 
            5, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorCarpiUlnaris", 
            "Mesh_Humerus", MWMath::Point3D{-0.0147, 0.2200, -0.0154}  , 
            "Mesh_MC4", MWMath::Point3D{0.0091, -0.0246, -0.0242}  , 
            {}, 
            {
                {"Via_Op", "Mesh_Ulna", MWMath::Point3D{0.0004, 0.2123, 0.0099}  },
                {"Via_Od", "Mesh_Ulna", MWMath::Point3D{0.0027, 0.0926, 0.0068}  },
                {"Via_Rad", "Mesh_Radius", MWMath::Point3D{0.0080, 0.0127, -0.0200}  }
            }, 
            5, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorCarpiRadialis", 
            "Mesh_Humerus", MWMath::Point3D{-0.0053, 0.2243, -0.0252}  , 
            "Mesh_MC1", MWMath::Point3D{0.0459, -0.0475, -0.0058}  , 
            {}, 
            {
                {"Via_Rad", "Mesh_Radius", MWMath::Point3D{0.0441, -0.0122, -0.0284}  }
            }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },

        // ==========================================
        // --- DAUMEN (Thumb / Finger 0) ---
        // ==========================================
        {
            "ExtensorPollicisLongus", 
            "Mesh_Ulna", MWMath::Point3D{0.0017, 0.1271, 0.0182}  , 
            "Mesh_Ulna", MWMath::Point3D{0.0169, 0.0466, 0.0087}  , // Proximal to Wrap
            {}, {}, 
            2, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "ExtensorPollicisBrevis", 
            "Mesh_Ulna", MWMath::Point3D{0.0199, 0.0750, 0.0060}  , 
            "Mesh_Radius", MWMath::Point3D{0.0522, 0.0046, -0.0063}  , 
            {}, 
            {
                {"Via_Uln", "Mesh_Ulna", MWMath::Point3D{0.0346, 0.0337, 0.0152}  }
            }, 
            3, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "AbductorPollicisLongus", 
            "Mesh_Ulna", MWMath::Point3D{0.0110, 0.1436, 0.0159}  , 
            "Mesh_MC0", MWMath::Point3D{0.0659, -0.0297, -0.0176}  , 
            {}, 
            {
                {"Via_Uln", "Mesh_Ulna", MWMath::Point3D{0.0446, 0.0373, 0.0068}  },
                {"Via_Rad", "Mesh_Radius", MWMath::Point3D{0.0513, 0.0058, -0.0057}  }
            }, 
            4, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "AbductorPollicisBrevis", 
            "Mesh_MC2", MWMath::Point3D{0.0341, -0.0206, -0.0234}  , 
            "Mesh_PP0", MWMath::Point3D{0.0814, -0.0650, -0.0154}  , 
            {}, 
            {
                {"Via_Om", "Mesh_MC2", MWMath::Point3D{0.0283, -0.0308, -0.0218}  },
                {"Via_Od", "Mesh_MC2", MWMath::Point3D{0.0322, -0.0445, -0.0164}  }
            }, 
            4, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorPollicisBrevis_Sup", 
            "Mesh_MC2", MWMath::Point3D{0.0291, -0.0423, -0.0171}  , 
            "Mesh_PP0", MWMath::Point3D{0.0740, -0.0681, -0.0189}  , 
            {}, {}, 
            2, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "OpponensPollicis", 
            "Mesh_MC2", MWMath::Point3D{0.0297, -0.0293, -0.0207}  , 
            "Mesh_MC0", MWMath::Point3D{0.0766, -0.0609, -0.0162}  , 
            {}, 
            {
                {"Via_Od", "Mesh_MC2", MWMath::Point3D{0.0389, -0.0223, -0.0236}  },
                {"Via_Ip", "Mesh_MC0", MWMath::Point3D{0.0633, -0.0293, -0.0189}  }
            }, 
            4, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorPollicisLongus", 
            "Mesh_Radius", MWMath::Point3D{0.0194, 0.1541, 0.0005}  , 
            "Mesh_DP0", MWMath::Point3D{0.0875, -0.0969, -0.0010}  , 
            {}, 
            {
                {"Via_Od", "Mesh_Radius", MWMath::Point3D{0.0331, 0.0499, -0.0023}  },
                {"Via_MC3", "Mesh_MC2", MWMath::Point3D{0.0379, -0.0223, -0.0209}  },
                {"Via_MC1", "Mesh_MC0", MWMath::Point3D{0.0594, -0.0636, -0.0156}  },
                {"Via_PP1a", "Mesh_PP0", MWMath::Point3D{0.0741, -0.0785, -0.0091}  },
                {"Via_PP1b", "Mesh_PP0", MWMath::Point3D{0.0820, -0.0882, -0.0047}  }
            }, 
            7, MWMath::Point3D{1.0, 0.0, 0.0}
        },

        // ==========================================
        // --- ZEIGEFINGER (Index / Finger 1) ---
        // ==========================================
        {
            "ExtensorIndicis", 
            "Mesh_Ulna", MWMath::Point3D{0.0016, 0.0996, 0.0141}  , 
            "Mesh_Ulna", MWMath::Point3D{0.0156, 0.0438, 0.0027}  , 
            {}, {}, 
            2, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorDigitorumSuperficialis_1", 
            "Mesh_MC2", MWMath::Point3D{0.0313, -0.0307, -0.0202}  , // Originiert distal von FDS common
            "Mesh_MP1", MWMath::Point3D{0.0427, -0.1442, -0.0049}  , 
            {}, 
            {
                {"Via_MC", "Mesh_MC1", MWMath::Point3D{0.0435, -0.0848, -0.0091}  },
                {"Via_PPa", "Mesh_PP1", MWMath::Point3D{0.0465, -0.1044, -0.0082}  },
                {"Via_PPb", "Mesh_PP1", MWMath::Point3D{0.0472, -0.1143, -0.0037}  },
                {"Via_PPc", "Mesh_PP1", MWMath::Point3D{0.0465, -0.1297, -0.0043}  }
            }, 
            6, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorDigitorumProfundus_1", 
            "Mesh_MC1", MWMath::Point3D{0.0316, -0.0279, -0.0169}  , 
            "Mesh_DP1", MWMath::Point3D{0.0455, -0.1632, -0.0057}  , 
            {}, 
            {
                {"Via_MCb", "Mesh_MC1", MWMath::Point3D{0.0369, -0.0552, -0.0106}  },
                {"Via_MCc", "Mesh_MC1", MWMath::Point3D{0.0441, -0.0862, -0.0067}  },
                {"Via_PPa", "Mesh_PP1", MWMath::Point3D{0.0479, -0.1061, -0.0036}  },
                {"Via_PPb", "Mesh_PP1", MWMath::Point3D{0.0471, -0.1272, -0.0014}  },
                {"Via_MPa", "Mesh_MP1", MWMath::Point3D{0.0456, -0.1387, -0.0063}  },
                {"Via_MPb", "Mesh_MP1", MWMath::Point3D{0.0457, -0.1494, -0.0048}  }
            }, 
            8, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "Lumbricales_0", 
            "Mesh_MC1", MWMath::Point3D{0.0369, -0.0552, -0.0106}  , 
            "Mesh_PP1", MWMath::Point3D{0.0593, -0.1030, -0.0025}  , 
            {}, {}, 
            2, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "DorsalInterosseus_0", 
            "Mesh_MC0", MWMath::Point3D{0.0530, -0.0392, -0.0062}  , 
            "Mesh_PP1", MWMath::Point3D{0.0558, -0.1079, 0.0005}  , 
            {}, 
            {
                {"Via_MC0", "Mesh_MC0", MWMath::Point3D{0.0617, -0.0567, -0.0060}  },
                {"Via_MC1a", "Mesh_MC1", MWMath::Point3D{0.0464, -0.0497, -0.0053}  },
                {"Via_MC1b", "Mesh_MC1", MWMath::Point3D{0.0469, -0.0710, -0.0020}  },
                {"Via_MC1c", "Mesh_MC1", MWMath::Point3D{0.0579, -0.0883, -0.0047}  }
            }, 
            6, MWMath::Point3D{1.0, 0.0, 0.0}
        },

        // ==========================================
        // --- MITTELFINGER (Middle / Finger 2) ---
        // ==========================================
        {
            "FlexorDigitorumSuperficialis_2", 
            "Mesh_MC2", MWMath::Point3D{0.0282, -0.0321, -0.0202}  , 
            "Mesh_MP2", MWMath::Point3D{0.0240, -0.1532, -0.0157}  , 
            {}, 
            {
                {"Via_MC", "Mesh_MC2", MWMath::Point3D{0.0288, -0.0831, -0.0092}  },
                {"Via_PPa", "Mesh_PP2", MWMath::Point3D{0.0281, -0.1052, -0.0097}  },
                {"Via_PPb", "Mesh_PP2", MWMath::Point3D{0.0272, -0.1112, -0.0092}  },
                {"Via_PPc", "Mesh_PP2", MWMath::Point3D{0.0239, -0.1302, -0.0110}  }
            }, 
            6, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorDigitorumProfundus_2", 
            "Mesh_MC2", MWMath::Point3D{0.0267, -0.0257, -0.0139}  , 
            "Mesh_DP2", MWMath::Point3D{0.0246, -0.1698, -0.0291}  , 
            {}, 
            {
                {"Via_MCa", "Mesh_MC2", MWMath::Point3D{0.0274, -0.0571, -0.0083}  },
                {"Via_MCb", "Mesh_MC2", MWMath::Point3D{0.0290, -0.0850, -0.0066}  },
                {"Via_PPa", "Mesh_PP2", MWMath::Point3D{0.0271, -0.1042, -0.0051}  },
                {"Via_PPb", "Mesh_PP2", MWMath::Point3D{0.0259, -0.1293, -0.0090}  },
                {"Via_MPa", "Mesh_MP2", MWMath::Point3D{0.0253, -0.1396, -0.0143}  },
                {"Via_MPb", "Mesh_MP2", MWMath::Point3D{0.0278, -0.1539, -0.0158}  }
            }, 
            8, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "Lumbricales1", 
            "Mesh_MC2", MWMath::Point3D{0.0274, -0.0571, -0.0083}  , 
            "Mesh_PP2", MWMath::Point3D{0.0383, -0.1016, -0.0031}  , 
            {}, {}, 
            2, MWMath::Point3D{1.0, 0.0, 0.0}
        },

        // ==========================================
        // --- RINGFINGER (Ring / Finger 3) ---
        // ==========================================
        {
            "FlexorDigitorumSuperficialis_3", 
            "Mesh_MC2", MWMath::Point3D{0.0253, -0.0302, -0.0201}  , 
            "Mesh_MP3", MWMath::Point3D{-0.0041, -0.1421, -0.0084}  , 
            {}, 
            {
                {"Via_MC", "Mesh_MC3", MWMath::Point3D{0.0130, -0.0805, -0.0130}  },
                {"Via_PPa", "Mesh_PP3", MWMath::Point3D{0.0081, -0.0985, -0.0084}  },
                {"Via_PPb", "Mesh_PP3", MWMath::Point3D{0.0062, -0.1050, -0.0062}  },
                {"Via_PPc", "Mesh_PP3", MWMath::Point3D{0.0019, -0.1196, -0.0043}  }
            }, 
            6, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorDigitorumProfundus_3", 
            "Mesh_MC2", MWMath::Point3D{0.0244, -0.0252, -0.0145}  , 
            "Mesh_DP3", MWMath::Point3D{-0.0030, -0.1634, -0.0130}  , 
            {}, 
            {
                {"Via_MCa", "Mesh_MC3", MWMath::Point3D{0.0180, -0.0586, -0.0115}  },
                {"Via_MCb", "Mesh_MC3", MWMath::Point3D{0.0137, -0.0792, -0.0092}  },
                {"Via_PPa", "Mesh_PP3", MWMath::Point3D{0.0061, -0.0983, -0.0053}  },
                {"Via_PPb", "Mesh_PP3", MWMath::Point3D{0.0030, -0.1212, -0.0037}  },
                {"Via_MPa", "Mesh_MP3", MWMath::Point3D{-0.0007, -0.1341, -0.0087}  },
                {"Via_MPb", "Mesh_MP3", MWMath::Point3D{-0.0015, -0.1467, -0.0091}  }
            }, 
            8, MWMath::Point3D{1.0, 0.0, 0.0}
        },

        // ==========================================
        // --- KLEINER FINGER (Little / Finger 4) ---
        // ==========================================
        {
            "FlexorDigitorumSuperficialis_4", 
            "Mesh_MC2", MWMath::Point3D{0.0226, -0.0289, -0.0214}  , 
            "Mesh_MP4", MWMath::Point3D{-0.0194, -0.1215, -0.0104}  , 
            {}, 
            {
                {"Via_MC", "Mesh_MC4", MWMath::Point3D{0.0069, -0.0765, -0.0104}  },
                {"Via_PPa", "Mesh_PP4", MWMath::Point3D{-0.0049, -0.0884, -0.0101}  },
                {"Via_PPb", "Mesh_PP4", MWMath::Point3D{-0.0108, -0.0991, -0.0053}  },
                {"Via_PPc", "Mesh_PP4", MWMath::Point3D{-0.0133, -0.1063, -0.0041}  }
            }, 
            6, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorDigitorumProfundus_4", 
            "Mesh_MC2", MWMath::Point3D{0.0358, -0.0266, -0.0139}  , 
            "Mesh_DP4", MWMath::Point3D{-0.0151, -0.1379, -0.0115}  , 
            {}, 
            {
                {"Via_MCa", "Mesh_MC4", MWMath::Point3D{0.0151, -0.0602, -0.0100}  },
                {"Via_MCb", "Mesh_MC4", MWMath::Point3D{0.0021, -0.0809, -0.0086}  },
                {"Via_PPa", "Mesh_PP4", MWMath::Point3D{-0.0080, -0.0923, -0.0053}  },
                {"Via_PPb", "Mesh_PP4", MWMath::Point3D{-0.0208, -0.1027, 0.0001}  },
                {"Via_MPa", "Mesh_MP4", MWMath::Point3D{-0.0156, -0.1165, -0.0087}  },
                {"Via_MPb", "Mesh_MP4", MWMath::Point3D{-0.0161, -0.1247, -0.0106}  }
            }, 
            8, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "FlexorDigitiMinimi", 
            "Mesh_MC2", MWMath::Point3D{-0.0134, -0.0923, -0.0113}  , 
            "Mesh_PP4", MWMath::Point3D{0.0147, -0.0210, -0.0160}  , 
            {}, {}, 
            2, MWMath::Point3D{1.0, 0.0, 0.0}
        },
        {
            "AbductorDigitiMinimi", 
            "Mesh_MC2", MWMath::Point3D{0.0120, -0.0140, -0.0261}  , 
            "Mesh_PP4", MWMath::Point3D{-0.0152, -0.0863, -0.0091}  , 
            {}, {}, 
            2, MWMath::Point3D{1.0, 0.0, 0.0}
        }
    };
    
    inline const MWMath::RotMatrix3x3 RotWristToWorld = MWMath::RotMatrix3x3();
    inline const MWMath::RotMatrix3x3 RotCarpalsToHumerus(0.0,  0.0, -1.0,0.0,  1.0,  0.0,1.0,  0.0,  0.0); // (X=lateral, Y=proximal, Z=posterior)
    inline const MWMath::RotMatrix3x3 RotCarpalsToUlna(0.0,  0.0, -1.0,0.0, -1.0,  0.0,-1.0,  0.0,  0.0); // (X=medial, Y=distal, Z=posterior)
    inline const MWMath::RotMatrix3x3 RotCarpalsToRadius(0.0,  0.0, -1.0,1.0,  0.0,  0.0,0.0, -1.0,  0.0); // (X=proximal, Y=medial, Z=posterior)
    inline const MWMath::Point3D WristToRadiusMM(0.0, 240.0, 25.0);
    inline const MWMath::Point3D WristToUlnaMM(-15.0, 260.0, -0.0); //WristToUlnaMM(-15.0, 250.0, -25.0);
    inline const MWMath::Point3D WristToHumerusMM(0.0, 560.0, 0.0);
    /* inline std::vector<MuscleDef> muscleDefsExcel = {
        
        // ==========================================
        // --- THENAR (Daumen / Finger 0) ---
        // ==========================================
        { "AbductorPollicisBrevis_1", "Mesh_Carpals", MWMath::Point3D{8.9, 1.6, 14.0} , "Mesh_PP0", MWMath::Point3D{-2.9, 1.7, 3.5} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{4.0, 5.0, 18.0} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "AbductorPollicisBrevis_2", "Mesh_Carpals", MWMath::Point3D{6.7, 0.2, 8.1} , "Mesh_PP0", MWMath::Point3D{-2.7, -0.028, 3.5} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{11.0, 5.0, 12.2} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "AbductorPollicisBrevis_3", "Mesh_Carpals", MWMath::Point3D{7.4, 3.0, 14.0} , "Mesh_PP0", MWMath::Point3D{-3.4, 0.64, 3.3} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{9.0, 5.0, 13.0} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "AbductorPollicisBrevis_4", "Mesh_Carpals", MWMath::Point3D{4.3, 3.7, 12.0} , "Mesh_PP0", MWMath::Point3D{-2.5, 0.42, 3.5} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{11.0, 5.0, 14.6} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "AbductorPollicisBrevis_5", "Mesh_Carpals", MWMath::Point3D{7.1, 2.2, 10.0} , "Mesh_PP0", MWMath::Point3D{-2.3, 1.8, 2.8} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{6.6, 5.0, 15.0} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },

        { "FlexorPollicisBrevisSuperficial_1", "Mesh_Carpals", MWMath::Point3D{8.1, 0.7, 17.0} , "Mesh_PP0", MWMath::Point3D{2.0, 6.1, 6.3} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{5.0, 7.5, 10.0} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "FlexorPollicisBrevisSuperficial_2", "Mesh_Carpals", MWMath::Point3D{8.1, 0.6, 15.0} , "Mesh_PP0", MWMath::Point3D{2.0, 4.4, 5.3} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{8.0, 6.4, 8.0} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },

        { "FlexorPollicisBrevisDeep_1", "Mesh_Carpals", MWMath::Point3D{3.02, -4.0, -1.0} , "Mesh_PP0", MWMath::Point3D{1.3, 0.5, 2.6} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{11.0, 6.5, 4.9} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },

        { "OpponensPollicis_1", "Mesh_Carpals", MWMath::Point3D{9.1, -5.4, 19.0} , "Mesh_MC0", MWMath::Point3D{0.34, -5.6, 1.6} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{0.48, 8.8, 7.3} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "OpponensPollicis_2", "Mesh_Carpals", MWMath::Point3D{9.7, -4.9, 15.0} , "Mesh_MC0", MWMath::Point3D{3.3, -0.4, 1.4} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{5.5, 8.3, 3.0} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "OpponensPollicis_3", "Mesh_Carpals", MWMath::Point3D{8.3, -10.0, 18.0} , "Mesh_MC0", MWMath::Point3D{1.8, -1.1, 1.8} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{2.6, 8.5, 5.2} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "OpponensPollicis_4", "Mesh_Carpals", MWMath::Point3D{6.2, -8.8, 25.0} , "Mesh_MC0", MWMath::Point3D{1.8, -10.5, 2.1} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{-2.1, 9.2, 9.8} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },

        { "AdductorPollicisOblique_1", "Mesh_Carpals", MWMath::Point3D{6.7, -10.0, -0.4} , "Mesh_PP0", MWMath::Point3D{5.93, 6.8, -0.24} , {}, {}, 2, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "AdductorPollicisOblique_2", "Mesh_Carpals", MWMath::Point3D{6.3, -14.0, -5.0} , "Mesh_PP0", MWMath::Point3D{5.77, 8.4, -0.91} , {}, {}, 2, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "AdductorPollicisOblique_3", "Mesh_Carpals", MWMath::Point3D{6.82, -23.0, -4.8} , "Mesh_PP0", MWMath::Point3D{6.3, 10.0, -4.7} , {}, {}, 2, MWMath::Point3D{1.0, 0.0, 0.0} },

        { "AdductorPollicisTransverse_1", "Mesh_MC2", MWMath::Point3D{1.4, 23.0, -0.36} , "Mesh_PP0", MWMath::Point3D{3.7, 10.0, -7.5} , {}, {}, 2, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "AdductorPollicisTransverse_2", "Mesh_MC2", MWMath::Point3D{1.7, -2.4, -0.078} , "Mesh_PP0", MWMath::Point3D{3.1, 4.9, -6.1} , {}, {}, 2, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "AdductorPollicisTransverse_3", "Mesh_MC2", MWMath::Point3D{4.1, -10.0, -0.41} , "Mesh_PP0", MWMath::Point3D{1.2, 0.8, -5.0} , {}, {}, 2, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "AdductorPollicisTransverse_4", "Mesh_MC2", MWMath::Point3D{1.3, 10.0, -0.4} , "Mesh_PP0", MWMath::Point3D{3.5, 7.5, -7.2} , {}, {}, 2, MWMath::Point3D{1.0, 0.0, 0.0} },

        // ==========================================
        // --- HYPOTHENAR (Kleiner Finger / Finger 4) ---
        // ==========================================
        { "FlexorDigitiMinimiBrevis_1", "Mesh_Carpals", MWMath::Point3D{12.0, -15.0, -21.0} , "Mesh_PP4", MWMath::Point3D{6.4, 11.0, -2.3} , {}, { {"Via_1", "Mesh_MC4", MWMath::Point3D{6.4, 1.9, -3.9} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },

        { "AbductorDigitiMinimi_1", "Mesh_Carpals", MWMath::Point3D{12.0, -6.1, -25.0} , "Mesh_PP4", MWMath::Point3D{-0.631, 4.85, -4.56} , {}, { {"Via_1", "Mesh_MC4", MWMath::Point3D{2.1, 5.34, -8.5} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "AbductorDigitiMinimi_2", "Mesh_Carpals", MWMath::Point3D{8.1, -5.6, -30.0} , "Mesh_PP4", MWMath::Point3D{1.25, 4.1, -5.85} , {}, { {"Via_1", "Mesh_MC4", MWMath::Point3D{6.14, 4.745, -11.2} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },

        { "OpponensDigitiMinimi_1", "Mesh_Carpals", MWMath::Point3D{10.0, -22.0, -23.0} , "Mesh_MC4", MWMath::Point3D{2.3, -0.12, -4.0} , {}, { {"Via_1", "Mesh_MC4", MWMath::Point3D{3.0, 9.2, -3.2} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "OpponensDigitiMinimi_2", "Mesh_Carpals", MWMath::Point3D{10.0, -20.0, -21.0} , "Mesh_MC4", MWMath::Point3D{2.1, -7.2, -4.0} , {}, { {"Via_1", "Mesh_MC4", MWMath::Point3D{4.7, 8.7, 0.84} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },

        // ==========================================
        // --- LUMBRICALS (Binnenmuskeln) ---
        // ==========================================
        // { "Lumbricals1_1", "Mesh_Carpals", MWMath::Point3D{0.0, 0.0, 0.0} , "Mesh_PP1", MWMath::Point3D{-2.0, 5.1, 3.7} , {}, {}, 2, MWMath::Point3D{1.0, 0.0, 0.0} },
        // { "Lumbricals2_1", "Mesh_Carpals", MWMath::Point3D{0.0, 0.0, 0.0} , "Mesh_PP2", MWMath::Point3D{-0.36, 10.0, 4.1} , {}, {}, 2, MWMath::Point3D{1.0, 0.0, 0.0} },
        // { "Lumbricals3_1", "Mesh_Carpals", MWMath::Point3D{0.0, 0.0, 0.0} , "Mesh_PP3", MWMath::Point3D{-3.4, 5.9, 3.5} , {}, {}, 2, MWMath::Point3D{1.0, 0.0, 0.0} },
        // { "Lumbricals3_2", "Mesh_Carpals", MWMath::Point3D{0.0, 0.0, 0.0} , "Mesh_PP3", MWMath::Point3D{-2.1, 6.7, 4.1} , {}, {}, 2, MWMath::Point3D{1.0, 0.0, 0.0} },
        // { "Lumbricals4_1", "Mesh_Carpals", MWMath::Point3D{0.0, 0.0, 0.0} , "Mesh_PP4", MWMath::Point3D{-0.038, 4.8, 3.2} , {}, {}, 2, MWMath::Point3D{1.0, 0.0, 0.0} },
        // { "Lumbricals4_2", "Mesh_Carpals", MWMath::Point3D{0.0, 0.0, 0.0} , "Mesh_PP4", MWMath::Point3D{-0.6, 2.7, 2.7} , {}, {}, 2, MWMath::Point3D{1.0, 0.0, 0.0} },

        // ==========================================
        // --- PALMAR INTEROSSEI ---
        // ==========================================
        //{ "PalmarInterossei1_1", "Mesh_MC1", MWMath::Point3D{0.2, -7.0, -2.6} , "Mesh_PP1", MWMath::Point3D{2.9, 8.1, -5.1} , {}, { {"Via_1", "Mesh_MC1", MWMath::Point3D{2.2, -14.0, -6.1} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "PalmarInterossei1_2", "Mesh_MC1", MWMath::Point3D{2.1, 11.0, -2.1} , "Mesh_PP1", MWMath::Point3D{1.51, 9.3, -5.6} , {}, { {"Via_1", "Mesh_MC1", MWMath::Point3D{2.0, -3.2, -7.5} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //
        //{ "PalmarInterossei2_1", "Mesh_MC3", MWMath::Point3D{-0.53, 5.8, 1.6} , "Mesh_PP3", MWMath::Point3D{2.2, 12.0, 5.6} , {}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{3.1, -0.49, 5.3} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "PalmarInterossei2_2", "Mesh_MC3", MWMath::Point3D{-0.42, -5.3, 2.5} , "Mesh_PP3", MWMath::Point3D{1.2, 14.0, 6.5} , {}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{1.6, -10.7, 5.3} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //
        //{ "PalmarInterossei3_1", "Mesh_MC4", MWMath::Point3D{-2.0, -4.8, 3.9} , "Mesh_PP4", MWMath::Point3D{2.3, 7.4, 4.7} , {}, { {"Via_1", "Mesh_MC4", MWMath::Point3D{0.0, -10.2, 5.2} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "PalmarInterossei3_2", "Mesh_MC4", MWMath::Point3D{-2.0, 5.9, 2.1} , "Mesh_PP4", MWMath::Point3D{3.1, 7.1, 4.2} , {}, { {"Via_1", "Mesh_MC4", MWMath::Point3D{-0.5, -0.98, 3.6} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },

        // ==========================================
        // --- DORSAL INTEROSSEI ---
        // ==========================================
        //{ "DorsalInterossei1_1", "Mesh_MC1", MWMath::Point3D{-2.7, 14.0, 5.5} , "Mesh_PP1", MWMath::Point3D{-2.0, 12.0, 9.1} , {}, { {"Via_1", "Mesh_MC1", MWMath::Point3D{0.81, 1.9, 8.5} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei1_2", "Mesh_MC1", MWMath::Point3D{-2.5, 1.4, 4.3} , "Mesh_PP1", MWMath::Point3D{-4.1, 12.0, 8.4} , {}, { {"Via_1", "Mesh_MC1", MWMath::Point3D{-1.6, 1.1, 6.5} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei1_3", "Mesh_MC1", MWMath::Point3D{-2.1, -10.0, 4.5} , "Mesh_PP1", MWMath::Point3D{-2.8, 14.0, 9.3} , {}, { {"Via_1", "Mesh_MC1", MWMath::Point3D{-1.1, 0.93, 11.0} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei1_4", "Mesh_MC0", MWMath::Point3D{1.2, 5.7, -4.3} , "Mesh_PP1", MWMath::Point3D{-2.1, 11.0, 8.2} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{7.0, 8.0, -15.0} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei1_5", "Mesh_MC0", MWMath::Point3D{0.88, -1.3, -4.7} , "Mesh_PP1", MWMath::Point3D{-3.8, 11.0, 7.6} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{7.0, 0.0, -15.0} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },

        //{ "DorsalInterossei2_1", "Mesh_MC1", MWMath::Point3D{-3.3, 13.0, -5.4} , "Mesh_PP2", MWMath::Point3D{-2.0, 10.0, 5.6} , {}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-2.0, -9.0, 8.5} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei2_2", "Mesh_MC1", MWMath::Point3D{-3.7, 10.0, -4.7} , "Mesh_PP2", MWMath::Point3D{-2.6, 10.0, 5.6} , {}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-2.6, -20.6, 8.5} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei2_3", "Mesh_MC1", MWMath::Point3D{-3.2, -9.5, -4.6} , "Mesh_PP2", MWMath::Point3D{-3.5, 10.0, 5.1} , {}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-3.5, -29.6, 8.5} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei2_4", "Mesh_MC2", MWMath::Point3D{-5.6, -9.2, 3.4} , "Mesh_PP2", MWMath::Point3D{-3.3, 15.0, 4.6} , {}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-4.3, -26.0, 5.9} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei2_5", "Mesh_MC2", MWMath::Point3D{-5.9, 0.49, 3.4} , "Mesh_PP2", MWMath::Point3D{-0.96, 15.0, 4.6} , {}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-3.6, -26.0, 5.4} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei2_6", "Mesh_MC2", MWMath::Point3D{-5.6, 9.9, 2.5} , "Mesh_PP2", MWMath::Point3D{-5.6, 15.0, 4.6} , {}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-2.9, -26.0, 5.1} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        
        //{ "DorsalInterossei3_1", "Mesh_MC2", MWMath::Point3D{-4.7, 9.8, -3.1} , "Mesh_PP2", MWMath::Point3D{-0.21, 16.0, -6.5} , {}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-2.3, -21.0, -7.0} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei3_2", "Mesh_MC2", MWMath::Point3D{-4.2, 0.95, -3.7} , "Mesh_PP2", MWMath::Point3D{-2.9, 16.0, -5.0} , {}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-3.0, -20.0, -6.5} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei3_3", "Mesh_MC2", MWMath::Point3D{-5.5, -8.8, -4.5} , "Mesh_PP2", MWMath::Point3D{-1.6, 17.0, -5.9} , {}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-3.5, -21.0, -8.1} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei3_4", "Mesh_MC3", MWMath::Point3D{-3.3, 7.7, 2.6} , "Mesh_PP2", MWMath::Point3D{-1.8, 10.0, -5.6} , {}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-1.3, 20.0, -9.4} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei3_5", "Mesh_MC3", MWMath::Point3D{-2.8, -6.1, 3.2} , "Mesh_PP2", MWMath::Point3D{-2.7, 10.0, -4.9} , {}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-2.8, -19.0, -9.9} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei3_6", "Mesh_MC3", MWMath::Point3D{-2.2, 0.96, 2.7} , "Mesh_PP2", MWMath::Point3D{-0.18, 10.0, -6.2} , {}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-1.7, -20.0, -8.6} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },

        //{ "DorsalInterossei4_1", "Mesh_MC3", MWMath::Point3D{-2.4, -0.64, -4.4} , "Mesh_PP3", MWMath::Point3D{-2.2, 17.0, -6.4} , {}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{-1.9, -18.0, -7.1} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei4_2", "Mesh_MC3", MWMath::Point3D{-2.1, -6.7, -4.1} , "Mesh_PP3", MWMath::Point3D{-3.3, 17.0, -6.1} , {}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{-1.8, -18.0, -5.7} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei4_3", "Mesh_MC3", MWMath::Point3D{-2.7, 6.0, -4.1} , "Mesh_PP3", MWMath::Point3D{-1.3, 17.0, -6.6} , {}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{-0.56, -18.0, -6.0} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei4_4", "Mesh_MC4", MWMath::Point3D{-4.2, 7.1, 1.1} , "Mesh_PP3", MWMath::Point3D{-1.6, 6.0, -6.0} , {}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{-2.4, -16.0, -7.8} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei4_5", "Mesh_MC4", MWMath::Point3D{-4.8, 0.85, 1.8} , "Mesh_PP3", MWMath::Point3D{-1.1, 5.0, -6.0} , {}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{-2.1, -18.0, -6.3} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "DorsalInterossei4_6", "Mesh_MC4", MWMath::Point3D{-4.9, -6.4, 3.2} , "Mesh_PP3", MWMath::Point3D{-3.3, 5.0, -5.0} , {}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{-2.3, -16.0, -7.0} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },

        // ==========================================
        // --- EXTRINSIC FLEXORS (Von Unterarm) ---
        // ==========================================
        { "FlexorCarpiUlnaris_1", "Mesh_Humerus", MWMath::Point3D{-17.5, -355.6, -20.3} , "Mesh_Carpals", MWMath::Point3D{10.77, 2.472, -20.86} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{16.64, 58.93, 16.68} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "FlexorCarpiUlnaris_2", "Mesh_Humerus", MWMath::Point3D{-20.85, -359.3, -20.17} , "Mesh_Carpals", MWMath::Point3D{10.73, 3.284, -17.91} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{23.29, 58.17, 10.73} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        
        { "FlexorCarpiRadialis_1", "Mesh_Humerus", MWMath::Point3D{-12.92, -359.6, -16.45} , "Mesh_MC1", MWMath::Point3D{2.84, 24.64, -4.447} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{14.98, 67.29, -21.02} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "FlexorCarpiRadialis_2", "Mesh_Humerus", MWMath::Point3D{-13.64, -355.4, -13.84} , "Mesh_MC1", MWMath::Point3D{2.553, 28.05, -5.668} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{12.88, 67.91, -18.02} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },

        { "FlexorPollicisLongus_1", "Mesh_Radius", MWMath::Point3D{-89.66, -8.361, -4.372} , "Mesh_DP0", MWMath::Point3D{0.997, 5.414, 3.404} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-182.0, -13.33, -10.5} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "FlexorPollicisLongus_2", "Mesh_Radius", MWMath::Point3D{-77.92, -8.64, -2.982} , "Mesh_DP0", MWMath::Point3D{0.997, 5.414, 3.404} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-182.0, -8.82, -10.2} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        
        //{ "FlexorDigitorumSuperficialis_3_1", "Mesh_Humerus", MWMath::Point3D{-18.92, -330.2, -23.85} , "Mesh_MP3", MWMath::Point3D{2.652, 7.257, 2.47} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{20.03, 99.78, -20.78} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "FlexorDigitorumSuperficialisHumeroulnar_3_2", "Mesh_Humerus", MWMath::Point3D{-17.97, -338.3, -24.45} , "Mesh_MP3", MWMath::Point3D{2.03, 9.25, 1.878} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{18.18, 98.96, -17.18} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "FlexorDigitorumSuperficialis_4_1", "Mesh_Humerus", MWMath::Point3D{-21.94, -368.0, -24.17} , "Mesh_MP4", MWMath::Point3D{1.495, 5.604, -1.007} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{20.99, 99.75, -12.12} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "FlexorDigitorumSuperficialisHumeroulnar_4_2", "Mesh_Humerus", MWMath::Point3D{-21.41, -370.6, -23.6} , "Mesh_MP4", MWMath::Point3D{1.457, 5.692, -1.092} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{18.99, 99.53, -5.606} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        
        { "FlexorDigitorumSuperficialis_1_1", "Mesh_Radius", MWMath::Point3D{-61.7, -3.448, -4.993} , "Mesh_MP1", MWMath::Point3D{2.116, 6.523, 2.057} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-73.54, 2.824, -9.478} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "FlexorDigitorumSuperficialis_1_1", "Mesh_Radius", MWMath::Point3D{0.,0.,0.} , "Mesh_MP1", MWMath::Point3D{10.,0. ,0.} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-73.54, 2.824, -9.478} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        
        //{ "FlexorDigitorumSuperficialisRadial_1_2", "Mesh_Radius", MWMath::Point3D{-46.97, -0.906, -3.337} , "Mesh_MP1", MWMath::Point3D{2.393, 7.813, 2.173} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-74.02, 8.087, -13.73} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "FlexorDigitorumSuperficialis_2_1", "Mesh_Radius", MWMath::Point3D{-33.96, 1.06, -1.918} , "Mesh_MP2", MWMath::Point3D{2.047, 9.653, 0.5864} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{12.48, 96.83, -25.0} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "FlexorDigitorumSuperficialisRadial_2_2", "Mesh_Radius", MWMath::Point3D{-21.14, 3.75, -2.235} , "Mesh_MP2", MWMath::Point3D{2.234, 8.336, -0.335} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{16.33, 98.35, -24.18} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        
        { "FlexorDigitorumProfundus_1_1", "Mesh_Ulna", MWMath::Point3D{20.25, 100.3, -17.71} , "Mesh_DP1", MWMath::Point3D{0.918, 1.99, 1.897} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{19.01, 109.3, -16.66} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "FlexorDigitorumProfundus_1_2", "Mesh_Ulna", MWMath::Point3D{18.29, 75.68, -11.29} , "Mesh_DP1", MWMath::Point3D{1.54, 3.546, 2.038} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{26.45, 110.3, -14.92} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        { "FlexorDigitorumProfundus_2_1", "Mesh_Ulna", MWMath::Point3D{20.25, 100.3, -17.71} , "Mesh_DP2", MWMath::Point3D{0.918, 1.99, 1.897} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{19.01, 109.3, -16.66} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "FlexorDigitorumProfundus_2_2", "Mesh_Ulna", MWMath::Point3D{18.29, 75.68, -11.29} , "Mesh_DP2", MWMath::Point3D{1.54, 3.546, 2.038} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{26.45, 110.3, -14.92} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        { "FlexorDigitorumProfundus_3_1", "Mesh_Ulna", MWMath::Point3D{20.25, 100.3, -17.71} , "Mesh_DP3", MWMath::Point3D{0.918, 1.99, 1.897} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{19.01, 109.3, -16.66} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "FlexorDigitorumProfundus_3_2", "Mesh_Ulna", MWMath::Point3D{18.29, 75.68, -11.29} , "Mesh_DP3", MWMath::Point3D{1.54, 3.546, 2.038} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{26.45, 110.3, -14.92} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        { "FlexorDigitorumProfundus_4_1", "Mesh_Ulna", MWMath::Point3D{20.25, 100.3, -17.71} , "Mesh_DP4", MWMath::Point3D{0.918, 1.99, 1.897} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{19.01, 109.3, -16.66} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "FlexorDigitorumProfundus_4_2", "Mesh_Ulna", MWMath::Point3D{18.29, 75.68, -11.29} , "Mesh_DP4", MWMath::Point3D{1.54, 3.546, 2.038} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{26.45, 110.3, -14.92} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },

        { "PalmarisLongus_1", "Mesh_Humerus", MWMath::Point3D{-18.52, -363.9, -24.39} , "Mesh_Carpals", MWMath::Point3D{0.295, -4.826, 0.209} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{18.71, 75.68, -13.82} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },

        // ==========================================
        // --- EXTRINSIC EXTENSORS (Table 7) ---
        // ==========================================
        { "ExtensorCarpiRadialisLongus_1", "Mesh_Humerus", MWMath::Point3D{25.37, -320.2, 2.568} , "Mesh_MC1", MWMath::Point3D{-8.58, 28.16, 4.778} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{4.335, -17.12, 2.269} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "ExtensorCarpiRadialisLongus_2", "Mesh_Humerus", MWMath::Point3D{32.2, -342.3, 1.9} , "Mesh_MC1", MWMath::Point3D{-8.56, 25.58, 4.415} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-3.796, -20.34, 10.88} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        { "ExtensorCarpiRadialisBrevis_1", "Mesh_Humerus", MWMath::Point3D{39.0, -331.2, -1.775} , "Mesh_MC2", MWMath::Point3D{-12.26, 25.31, -0.861} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-33.14, -9.632, 16.33} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "ExtensorCarpiRadialisBrevis_2", "Mesh_Humerus", MWMath::Point3D{39.98, -336.8, -4.294} , "Mesh_MC2", MWMath::Point3D{-11.96, 27.73, -0.582} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-34.24, -16.86, 11.73} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        { "ExtensorCarpiUlnaris_1", "Mesh_Humerus", MWMath::Point3D{24.11, -346.4, 0.571} , "Mesh_MC4", MWMath::Point3D{-2.222, 17.97, -4.251} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{-4.28, 109.4, 4.36} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "ExtensorCarpiUlnaris_2", "Mesh_Humerus", MWMath::Point3D{30.96, -342.0, 2.139} , "Mesh_MC4", MWMath::Point3D{-2.792, 17.03, -4.285} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{-2.785, 109.5, 2.853} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },


        { "ExtensorDigitorum_1_1", "Mesh_Humerus", MWMath::Point3D{34.87, -339.1, -3.738} , "Mesh_DP1", MWMath::Point3D{-4.47, 4.714, -1.288} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.48, -11.22, 13.74} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "ExtensorDigitorum_1_2", "Mesh_Humerus", MWMath::Point3D{30.49, -340.9, -4.776} , "Mesh_DP1", MWMath::Point3D{-4.524, 5.019, 0.458} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.21, -4.973, 14.83} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        { "ExtensorDigitorum_2_1", "Mesh_Humerus", MWMath::Point3D{34.87, -339.1, -3.738} , "Mesh_DP2", MWMath::Point3D{-4.47, 4.714, -1.288} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.48, -11.22, 13.74} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "ExtensorDigitorum_2_2", "Mesh_Humerus", MWMath::Point3D{30.49, -340.9, -4.776} , "Mesh_DP2", MWMath::Point3D{-4.524, 5.019, 0.458} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.21, -4.973, 14.83} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        { "ExtensorDigitorum_3_1", "Mesh_Humerus", MWMath::Point3D{34.87, -339.1, -3.738} , "Mesh_DP3", MWMath::Point3D{-4.47, 4.714, -1.288} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.48, -11.22, 13.74} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "ExtensorDigitorum_3_2", "Mesh_Humerus", MWMath::Point3D{30.49, -340.9, -4.776} , "Mesh_DP3", MWMath::Point3D{-4.524, 5.019, 0.458} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.21, -4.973, 14.83} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        { "ExtensorDigitorum_4_1", "Mesh_Humerus", MWMath::Point3D{34.87, -339.1, -3.738} , "Mesh_DP4", MWMath::Point3D{-4.47, 4.714, -1.288} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.48, -11.22, 13.74} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "ExtensorDigitorum_4_2", "Mesh_Humerus", MWMath::Point3D{30.49, -340.9, -4.776} , "Mesh_DP4", MWMath::Point3D{-4.524, 5.019, 0.458} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.21, -4.973, 14.83} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },


        { "ExtensorDigitiMinimi_1", "Mesh_Humerus", MWMath::Point3D{26.15, -371.1, -3.214} , "Mesh_DP4", MWMath::Point3D{-3.714, 8.692, -0.0449} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{-15.26, 121.7, 7.688} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "ExtensorDigitiMinimi_2", "Mesh_Humerus", MWMath::Point3D{23.89, -371.0, -2.115} , "Mesh_DP4", MWMath::Point3D{-5.486, 4.598, -3.015} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{-0.121, 121.8, 6.031} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },

        { "ExtensorIndicis_1", "Mesh_Ulna", MWMath::Point3D{10.83, 226.3, -3.307} , "Mesh_DP1", MWMath::Point3D{-4.342, 6.562, 1.834} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-219.5, -10.35, 25.02} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} },
        //{ "ExtensorIndicis_2", "Mesh_Ulna", MWMath::Point3D{11.54, 210.5, -2.816} , "Mesh_DP1", MWMath::Point3D{-4.253, 5.784, 1.265} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-219.8, -5.307, 22.61} } }, 3, MWMath::Point3D{1.0, 0.0, 0.0} }
        
    };
     */
    
    inline std::vector<MuscleDef> muscleDefsExcel = {
        
        // ==========================================
        // --- THENAR (Daumen / Finger 0) ---
        // ==========================================
        /* 
        { "AbductorPollicisBrevis_1", "Mesh_Carpals", MWMath::Point3D{8.9, 1.6, 14.0} , "Mesh_PP0", MWMath::Point3D{-2.9, 1.7, 3.5} , {"Mesh_MC0", "Mesh_MCP0_Joint"}, {}, 10, OTHERMUSCLECOLOR, 1.0 },
        //{ "AbductorPollicisBrevis_2", "Mesh_Carpals", MWMath::Point3D{6.7, 0.2, 8.1} , "Mesh_PP0", MWMath::Point3D{-2.7, -0.028, 3.5} ,{"Mesh_MC0", "Mesh_MCP0_Joint"}, {}, 10, OTHERMUSCLECOLOR },
        //{ "AbductorPollicisBrevis_3", "Mesh_Carpals", MWMath::Point3D{7.4, 3.0, 14.0} , "Mesh_PP0", MWMath::Point3D{-3.4, 0.64, 3.3} , {"Mesh_MC0", "Mesh_MCP0_Joint"}, {}, 10, OTHERMUSCLECOLOR },
        //{ "AbductorPollicisBrevis_4", "Mesh_Carpals", MWMath::Point3D{4.3, 3.7, 12.0} , "Mesh_PP0", MWMath::Point3D{-2.5, 0.42, 3.5} , {"Mesh_MC0", "Mesh_MCP0_Joint"}, {}, 10, OTHERMUSCLECOLOR },
        //{ "AbductorPollicisBrevis_5", "Mesh_Carpals", MWMath::Point3D{7.1, 2.2, 10.0} , "Mesh_PP0", MWMath::Point3D{-2.3, 1.8, 2.8} ,  {"Mesh_MC0", "Mesh_MCP0_Joint"}, {}, 10, OTHERMUSCLECOLOR },

        { "FlexorPollicisBrevisSuperficial_1", "Mesh_Carpals", MWMath::Point3D{8.1, 0.7, 17.0} , "Mesh_PP0", MWMath::Point3D{2.0, 6.1, 6.3} , {"Mesh_MC0", "Mesh_MCP0_Joint",}, {}, 10, FLEXORCOLOR },
        //{ "FlexorPollicisBrevisSuperficial_2", "Mesh_Carpals", MWMath::Point3D{8.1, 0.6, 15.0} , "Mesh_PP0", MWMath::Point3D{2.0, 4.4, 5.3} , {"Mesh_MC0", "Mesh_MCP0_Joint",}, {}, 10, FLEXORCOLOR },

        { "FlexorPollicisBrevis_1", "Mesh_Carpals", MWMath::Point3D{3.02, -4.0, -1.0} , "Mesh_PP0", MWMath::Point3D{1.3, 0.5, 2.6} , {"Mesh_MC0", "Mesh_MCP0_Joint", "Mesh_PP0"}, {}, 10, FLEXORCOLOR },

        { "OpponensPollicis_1", "Mesh_Carpals", MWMath::Point3D{9.1, -5.4, 19.0} , "Mesh_MC0", MWMath::Point3D{0.34, -5.6, 1.6} ,{"Mesh_MC0", "Mesh_MCP0_Joint"}, {}, 10, OTHERMUSCLECOLOR },
        //{ "OpponensPollicis_2", "Mesh_Carpals", MWMath::Point3D{9.7, -4.9, 15.0} , "Mesh_MC0", MWMath::Point3D{3.3, -0.4, 1.4} , {"Mesh_MC0", "Mesh_MCP0_Joint"}, {}, 10, OTHERMUSCLECOLOR },
        //{ "OpponensPollicis_3", "Mesh_Carpals", MWMath::Point3D{8.3, -10.0, 18.0} , "Mesh_MC0", MWMath::Point3D{1.8, -1.1, 1.8} ,{"Mesh_MC0", "Mesh_MCP0_Joint"}, {}, 10, OTHERMUSCLECOLOR },
        //{ "OpponensPollicis_4", "Mesh_Carpals", MWMath::Point3D{6.2, -8.8, 25.0} , "Mesh_MC0", MWMath::Point3D{1.8, -10.5, 2.1} ,{"Mesh_MC0", "Mesh_MCP0_Joint"}, {}, 10, OTHERMUSCLECOLOR },

        { "AdductorPollicisOblique_1", "Mesh_Carpals", MWMath::Point3D{6.7, -10.0, -0.4} , "Mesh_PP0", MWMath::Point3D{5.93, 6.8, -0.24} , {"Mesh_MC0", "Mesh_MCP0_Joint"}, {}, 10, OTHERMUSCLECOLOR },
        //{ "AdductorPollicisOblique_2", "Mesh_Carpals", MWMath::Point3D{6.3, -14.0, -5.0} , "Mesh_PP0", MWMath::Point3D{5.77, 8.4, -0.91} , {"Mesh_MC0", "Mesh_MCP0_Joint"}, {}, 2, OTHERMUSCLECOLOR },
        //{ "AdductorPollicisOblique_3", "Mesh_Carpals", MWMath::Point3D{6.82, -23.0, -4.8} , "Mesh_PP0", MWMath::Point3D{6.3, 10.0, -4.7} , {"Mesh_MC0", "Mesh_MCP0_Joint"}, {}, 2, OTHERMUSCLECOLOR },
        { "AdductorPollicisTransverse_1", "Mesh_MC2", MWMath::Point3D{1.4, 23.0, -0.36} , "Mesh_PP0", MWMath::Point3D{3.7, 10.0, -7.5}, {"Mesh_MC0", "Mesh_MC1", "Mesh_MCP0_Joint"}, {}, 10, OTHERMUSCLECOLOR },
        //{ "AdductorPollicisTransverse_2", "Mesh_MC2", MWMath::Point3D{1.7, -2.4, -0.078} , "Mesh_PP0", MWMath::Point3D{3.1, 4.9, -6.1}, {"Mesh_MC0", "Mesh_MC1", "Mesh_MC2"}, {}, 2, OTHERMUSCLECOLOR },
        //{ "AdductorPollicisTransverse_3", "Mesh_MC2", MWMath::Point3D{4.1, -10.0, -0.41} , "Mesh_PP0", MWMath::Point3D{1.2, 0.8, -5.0}, {"Mesh_MC0", "Mesh_MC1", "Mesh_MC2"}, {}, 2, OTHERMUSCLECOLOR },
        //{ "AdductorPollicisTransverse_4", "Mesh_MC2", MWMath::Point3D{1.3, 10.0, -0.4} , "Mesh_PP0", MWMath::Point3D{3.5, 7.5, -7.2} ,  {"Mesh_MC0", "Mesh_MC1", "Mesh_MC2"}, {}, 2, OTHERMUSCLECOLOR },
         */

        //{ "AbductorPollicisLongus_1", "Mesh_Ulna", MWMath::Point3D{1.698, 100.1, 10.71} , "Mesh_MC0", MWMath::Point3D{-6.95, 12.86, -4.657} , {"MeshVP_APL", "Mesh_Carpals"}, {}, 30, EXTENSORCOLOR, 1.0 },
        //{ "AbductorPollicisLongus_2", "Mesh_Ulna", MWMath::Point3D{2.222, 71.91, 13.82} , "Mesh_MC0", MWMath::Point3D{-7.12, 10.74, -4.488} , {"MeshVP_APL", "Mesh_Carpals"}, {}, 50, EXTENSORCOLOR, 1.0 },
        //{ "AbductorPollicisLongus_3", "Mesh_Radius", MWMath::Point3D{-80.71, -1.14, 3.914} , "Mesh_MC0", MWMath::Point3D{-8.15, 11.34, -2.948} , {"MeshVP_APL", "Mesh_Carpals"}, {}, 50, EXTENSORCOLOR, 1.0 },
        //{ "AbductorPollicisLongus_4", "Mesh_Radius", MWMath::Point3D{-68.12, 1.48, 4.362} , "Mesh_MC0", MWMath::Point3D{-7.44, 10.35, -4.149} , {"MeshVP_APL", "Mesh_Carpals"}, {}, 50, EXTENSORCOLOR, 1.0 },

        {"ExtensorPollicisLongus_1", "Mesh_Ulna", MWMath::Point3D{11.53, 197.9, 3.395} , "Mesh_DP0", MWMath::Point3D{-5.998, 7.305, -1.757} , {"MeshVP_EPL", "Mesh_Carpals", "MeshVPE_MC0", "Mesh_MC0_Joint", "MeshVPE_PP0", "Mesh_PP0", "Mesh_IP0_Joint", "MeshVPE_DP0"}, {}, 30, OTHERMUSCLECOLOR },
        // {"ExtensorPollicisLongus_2", "Mesh_Ulna", MWMath::Point3D{9.356, 134.1, 1.838} , "Mesh_DP0", MWMath::Point3D{-6.048, 8.096, -0.257} , {"MeshVP_EPL", "Mesh_Carpals", "MeshVPE_MC0", "Mesh_MC0_Joint", "Mesh_PP0", "Mesh_IP0_Joint", "Mesh_PP0_Joint"}, {}, 40, OTHERMUSCLECOLOR },

        {"ExtensorPollicisBrevis_1", "Mesh_Radius", MWMath::Point3D{-167.6, 2.166, -0.1964} , "Mesh_PP0", MWMath::Point3D{-7.834, 13.34, -3.043} , {"MeshVP_EPB", "Mesh_Carpals", "MeshVPE_MC0", "Mesh_MCP0_Joint"}, {}, 30, OTHERMUSCLECOLOR },
        // {"ExtensorPollicisBrevis_2", "Mesh_Radius", MWMath::Point3D{-153.8,-0.05,0.0524} , "Mesh_PP0", MWMath::Point3D{-8.295, 11.88, -2.4} , {"MeshVP_EPB", "Mesh_Carpals", "MeshVPE_MC0", "Mesh_MCP0_Joint"}, {}, 30, OTHERMUSCLECOLOR },

        // ==========================================
        // --- HYPOTHENAR (Kleiner Finger / Finger 4) ---
        // ==========================================
        /* 
        { "FlexorDigitiMinimiBrevis_1", "Mesh_Carpals", MWMath::Point3D{12.0, -15.0, -21.0} , "Mesh_PP4", MWMath::Point3D{6.4, 11.0, -2.3} , {"Mesh_MC4", "Mesh_MCP4_Joint"}, {}, 60, FLEXORCOLOR },
        
        { "AbductorDigitiMinimi_1", "Mesh_Carpals", MWMath::Point3D{12.0, -6.1, -25.0} , "Mesh_PP4", MWMath::Point3D{-0.631, 4.85, -4.56} , {"Mesh_MC4", "Mesh_MCP4_Joint"}, {}, 10, OTHERMUSCLECOLOR },
        //{ "AbductorDigitiMinimi_2", "Mesh_Carpals", MWMath::Point3D{8.1, -5.6, -30.0} , "Mesh_PP4", MWMath::Point3D{1.25, 4.1, -5.85} ,   {"Mesh_MC4", "Mesh_MCP4_Joint"}, {} }, 10, OTHERMUSCLECOLOR },

        { "OpponensDigitiMinimi_1", "Mesh_Carpals", MWMath::Point3D{10.0, -22.0, -23.0} , "Mesh_MC4", MWMath::Point3D{2.3, -0.12, -4.0} ,   {"Mesh_MC4", "Mesh_MCP4_Joint"}, {}, 10, OTHERMUSCLECOLOR },
        //{ "OpponensDigitiMinimi_2", "Mesh_Carpals", MWMath::Point3D{10.0, -20.0, -21.0} , "Mesh_MC4", MWMath::Point3D{2.1, -7.2, -4.0} ,  {"Mesh_MC4", "Mesh_MCP4_Joint"}, {}, 10, OTHERMUSCLECOLOR },
       
        
        
        // ==========================================
        // --- LUMBRICALS (Binnenmuskeln) ---
        // ==========================================
        { "Lumbricals1_1", "Mesh_MC1", MWMath::Point3D{10.0, 0.0, 0.0} , "Mesh_PP1", MWMath::Point3D{-2.0, 5.1, 3.7} ,  {"Mesh_MC1", "Mesh_MCP1_Joint"}, {}, 10, LUMBRICALCOLOR },
        { "Lumbricals2_1", "Mesh_MC2", MWMath::Point3D{10.0, 0.0, 0.0} , "Mesh_PP2", MWMath::Point3D{-0.36, 10.0, 4.1}, {"Mesh_MC2", "Mesh_MCP2_Joint"}, {}, 10, LUMBRICALCOLOR },
        { "Lumbricals3_1", "Mesh_MC2", MWMath::Point3D{10.0, 0.0, 0.0} , "Mesh_PP3", MWMath::Point3D{-3.4, 5.9, 3.5} ,  {"Mesh_MC2", "Mesh_MCP3_Joint"}, {}, 10, LUMBRICALCOLOR },
        { "Lumbricals3_2", "Mesh_MC3", MWMath::Point3D{10.0, 0.0, 0.0} , "Mesh_PP3", MWMath::Point3D{-2.1, 6.7, 4.1} ,  {"Mesh_MC3", "Mesh_MCP3_Joint"}, {}, 10, LUMBRICALCOLOR },
        { "Lumbricals4_1", "Mesh_MC3", MWMath::Point3D{10.0, 0.0, 0.0} , "Mesh_PP4", MWMath::Point3D{-0.038, 4.8, 3.2}, {"Mesh_MC3", "Mesh_MCP4_Joint"}, {}, 10, LUMBRICALCOLOR },
        { "Lumbricals4_2", "Mesh_MC4", MWMath::Point3D{10.0, 0.0, 0.0} , "Mesh_PP4", MWMath::Point3D{-0.6, 2.7, 2.7} ,  {"Mesh_MC4", "Mesh_MCP4_Joint"}, {}, 10, LUMBRICALCOLOR },
        
        
        // ==========================================
        // --- PALMAR INTEROSSEI ---
        // ==========================================
        { "PalmarInterossei1_1", "Mesh_MC1", MWMath::Point3D{0.2, -7.0, -2.6} , "Mesh_PP1", MWMath::Point3D{2.9, 8.1, -5.1} ,  {"Mesh_MCP1_Joint", "Mesh_PP1"}, {}, 10, PINTEROSSEICCOLOR },
        //{ "PalmarInterossei1_2", "Mesh_MC1", MWMath::Point3D{2.1, 11.0, -2.1} , "Mesh_PP1", MWMath::Point3D{1.51, 9.3, -5.6},{"Mesh_MCP1_Joint"}, { {"Via_1", "Mesh_MC1", MWMath::Point3D{2.0, -3.2, -7.5} } }, 10, PINTEROSSEICCOLOR },
        { "PalmarInterossei2_1", "Mesh_MC3", MWMath::Point3D{-0.53, 5.8, 1.6} , "Mesh_PP3", MWMath::Point3D{2.2, 12.0, 5.6} , {"Mesh_MCP3_Joint", "Mesh_PP3"}, {}, 10, PINTEROSSEICCOLOR },
        //{ "PalmarInterossei2_2", "Mesh_MC3", MWMath::Point3D{-0.42, -5.3, 2.5} , "Mesh_PP3", MWMath::Point3D{1.2, 14.0, 6.5} , {"Mesh_MCP3_Joint"}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{1.6, -10.7, 5.3} } }, 10, PINTEROSSEICCOLOR },
        { "PalmarInterossei3_1", "Mesh_MC4", MWMath::Point3D{-2.0, -4.8, 3.9} , "Mesh_PP4", MWMath::Point3D{2.3, 7.4, 4.7} , {"Mesh_MCP4_Joint", "Mesh_PP4"}, {}, 10, PINTEROSSEICCOLOR },
        //{ "PalmarInterossei3_2", "Mesh_MC4", MWMath::Point3D{-2.0, 5.9, 2.1} , "Mesh_PP4", MWMath::Point3D{3.1, 7.1, 4.2} , {"Mesh_MCP4_Joint"}, { {"Via_1", "Mesh_MC4", MWMath::Point3D{-0.5, -0.98, 3.6} } }, 10, PINTEROSSEICCOLOR },

        // ==========================================
        // --- DORSAL INTEROSSEI ---
        // ==========================================
        { "DorsalInterossei1_1", "Mesh_MC1", MWMath::Point3D{-2.7, 14.0, 5.5} , "Mesh_PP1", MWMath::Point3D{-2.0, 12.0, 9.1} , {"Mesh_MCP1_Joint", "Mesh_PP1"}, {}, 10, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei1_2", "Mesh_MC1", MWMath::Point3D{-2.5, 1.4, 4.3} , "Mesh_PP1", MWMath::Point3D{-4.1, 12.0, 8.4} , {"Mesh_MCP1_Joint"}, {}, 10, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei1_3", "Mesh_MC1", MWMath::Point3D{-2.1, -10.0, 4.5} , "Mesh_PP1", MWMath::Point3D{-2.8, 14.0, 9.3} , {"Mesh_MCP1_Joint"}, {} }, 10, DINTEROSSEICCOLOR },
        { "DorsalInterossei1_4", "Mesh_MC0", MWMath::Point3D{1.2, 5.7, -4.3} , "Mesh_PP1", MWMath::Point3D{-2.1, 11.0, 8.2} , {"Mesh_MCP1_Joint", "Mesh_PP1"}, {}, 10, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei1_5", "Mesh_MC0", MWMath::Point3D{0.88, -1.3, -4.7} , "Mesh_PP1", MWMath::Point3D{-3.8, 11.0, 7.6} , {"Mesh_MCP1_Joint", "Mesh_PP1"}, {}, 10, DINTEROSSEICCOLOR },

        { "DorsalInterossei2_1", "Mesh_MC1", MWMath::Point3D{-3.3, 13.0, -5.4} , "Mesh_PP2", MWMath::Point3D{-2.0, 10.0, 5.6} , {"Mesh_MCP2_Joint", "Mesh_PP2"}, {}, 10, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei2_2", "Mesh_MC1", MWMath::Point3D{-3.7, 10.0, -4.7} , "Mesh_PP2", MWMath::Point3D{-2.6, 10.0, 5.6} , {"Mesh_MCP2_Joint", "Mesh_PP2"}, {}, 10, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei2_3", "Mesh_MC1", MWMath::Point3D{-3.2, -9.5, -4.6} , "Mesh_PP2", MWMath::Point3D{-3.5, 10.0, 5.1} , {"Mesh_MCP2_Joint", "Mesh_PP2"}, {}, 10, DINTEROSSEICCOLOR },
        { "DorsalInterossei2_4", "Mesh_MC2", MWMath::Point3D{-5.6, -9.2, 3.4} , "Mesh_PP2", MWMath::Point3D{-3.3, 15.0, 4.6} , {"Mesh_MCP2_Joint", "Mesh_PP2"}, {}, 10, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei2_5", "Mesh_MC2", MWMath::Point3D{-5.9, 0.49, 3.4} , "Mesh_PP2", MWMath::Point3D{-0.96, 15.0, 4.6} , {"Mesh_MCP2_Joint", "Mesh_PP2"}, {}, 10, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei2_6", "Mesh_MC2", MWMath::Point3D{-5.6, 9.9, 2.5} , "Mesh_PP2", MWMath::Point3D{-5.6, 15.0, 4.6} , {"Mesh_MCP2_Joint", "Mesh_PP2"}, {}, 10, DINTEROSSEICCOLOR },
    
        { "DorsalInterossei3_1", "Mesh_MC2", MWMath::Point3D{-4.7, 9.8, -3.1} , "Mesh_PP2", MWMath::Point3D{-0.21, 16.0, -6.5} , {"Mesh_MCP2_Joint", "Mesh_PP2"}, {}, 10, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei3_2", "Mesh_MC2", MWMath::Point3D{-4.2, 0.95, -3.7} , "Mesh_PP2", MWMath::Point3D{-2.9, 16.0, -5.0} , {"Mesh_MCP2_Joint", "Mesh_PP2"}, {}, 10, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei3_3", "Mesh_MC2", MWMath::Point3D{-5.5, -8.8, -4.5} , "Mesh_PP2", MWMath::Point3D{-1.6, 17.0, -5.9} , {"Mesh_MCP2_Joint", "Mesh_PP2"}, {}, 10, DINTEROSSEICCOLOR },
        { "DorsalInterossei3_4", "Mesh_MC3", MWMath::Point3D{-3.3, 7.7, 2.6} , "Mesh_PP2", MWMath::Point3D{-1.8, 10.0, -5.6} , {"Mesh_MCP2_Joint", "Mesh_PP2"}, {}, 10, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei3_5", "Mesh_MC3", MWMath::Point3D{-2.8, -6.1, 3.2} , "Mesh_PP2", MWMath::Point3D{-2.7, 10.0, -4.9} , {"Mesh_MCP2_Joint", "Mesh_PP2"}, {}, 10, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei3_6", "Mesh_MC3", MWMath::Point3D{-2.2, 0.96, 2.7} , "Mesh_PP2", MWMath::Point3D{-0.18, 10.0, -6.2} , {"Mesh_MCP2_Joint", "Mesh_PP2"}, {}, 10, DINTEROSSEICCOLOR },

        { "DorsalInterossei4_1", "Mesh_MC3", MWMath::Point3D{-2.4, -0.64, -4.4} , "Mesh_PP3", MWMath::Point3D{-2.2, 17.0, -6.4} , {"Mesh_MCP3_Joint", "Mesh_PP3"}, {}, 10, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei4_2", "Mesh_MC3", MWMath::Point3D{-2.1, -6.7, -4.1} , "Mesh_PP3", MWMath::Point3D{-3.3, 17.0, -6.1} , {"Mesh_MCP3_Joint", "Mesh_PP3"}, {}, 10, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei4_3", "Mesh_MC3", MWMath::Point3D{-2.7, 6.0, -4.1} , "Mesh_PP3", MWMath::Point3D{-1.3, 17.0, -6.6} , {"Mesh_MCP3_Joint", "Mesh_PP3"}, {}, 10, DINTEROSSEICCOLOR },
        { "DorsalInterossei4_4", "Mesh_MC4", MWMath::Point3D{-4.2, 7.1, 1.1} , "Mesh_PP3", MWMath::Point3D{-1.6, 6.0, -6.0} , {"Mesh_MCP3_Joint", "Mesh_PP3"}, {}, 10, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei4_5", "Mesh_MC4", MWMath::Point3D{-4.8, 0.85, 1.8} , "Mesh_PP3", MWMath::Point3D{-1.1, 5.0, -6.0} , {"Mesh_MCP3_Joint", "Mesh_PP3"}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{-2.1, -18.0, -6.3} } }, 10, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei4_6", "Mesh_MC4", MWMath::Point3D{-4.9, -6.4, 3.2} , "Mesh_PP3", MWMath::Point3D{-3.3, 5.0, -5.0} , {"Mesh_MCP3_Joint", "Mesh_PP3"}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{-2.3, -16.0, -7.0} } }, 10, DINTEROSSEICCOLOR },
        
       
        
        // ==========================================
        // --- EXTRINSIC FLEXORS (Von Unterarm) ---
        // ==========================================
        
        { "FlexorCarpiUlnaris_1", "Mesh_Humerus", MWMath::Point3D{-17.5, -355.6, -20.3} , "Mesh_Carpals", MWMath::Point3D{10.77, 2.472, -20.86} ,   {"MeshVP_FCU"}, {}, 30, FLEXORCOLOR, 1.0 },
        //{ "FlexorCarpiUlnaris_2", "Mesh_Humerus", MWMath::Point3D{-20.85, -359.3, -20.17} , "Mesh_Carpals", MWMath::Point3D{10.73, 3.284, -17.91},{"MeshVP_FCU"}, {}, 30, FLEXORCOLOR, 1.0 },
        
        { "FlexorCarpiRadialis_1", "Mesh_Humerus", MWMath::Point3D{-12.92, -359.6, -16.45} , "Mesh_MC1", MWMath::Point3D{2.84, 24.64, -4.447} ,  {"MeshVP_FCR"}, {}, 30, FLEXORCOLOR, 1.0 },
        //{ "FlexorCarpiRadialis_2", "Mesh_Humerus", MWMath::Point3D{-13.64, -355.4, -13.84} , "Mesh_MC1", MWMath::Point3D{2.553, 28.05, -5.668},{"MeshVP_FCR"}, {}, 30, FLEXORCOLOR, 1.0 },
       
        */
        { "FlexorPollicisLongus_1", "Mesh_Radius", MWMath::Point3D{-89.66, -8.361, -4.372} , "Mesh_DP0", MWMath::Point3D{0.997, 5.414, 3.404} , {"MeshVP_FPL", "MeshVPF_MC0", "Mesh_MC0", "Mesh_MCP0_Joint", "Mesh_PP0", "Mesh_IP0_Joint", "MeshVPF_PP0"}, {}, 50, FLEXORCOLOR, 1.0 },
        //{ "FlexorPollicisLongus_2", "Mesh_Radius", MWMath::Point3D{-77.92, -8.64, -2.982} , "Mesh_DP0", MWMath::Point3D{0.997, 5.414, 3.404} ,{"MeshVP_FPL", "MeshVP_MC0", "Mesh_MC0", "Mesh_MCP0_Joint", "Mesh_PP0", "Mesh_MCP0_Joint", "MeshVPF_PP0"}, {}, 10, FLEXORCOLOR, 1.0 },
        /*
        
        { "FlexorDigitorumSuperficialis_1_1", "Mesh_Radius", MWMath::Point3D{-61.7, -3.448, -4.993} , "Mesh_MP1", MWMath::Point3D{2.116, 6.523, 2.057} ,  {"MeshVP_FDS", "Mesh_MC1", "MeshVPF_MC1", "Mesh_MCP1_Joint", "Mesh_PP1" , "MeshVPF_PP1", "Mesh_PIP1_Joint", "Mesh_MP1"}, { }, 60, FLEXORCOLOR, 1.0},
        //{ "FlexorDigitorumSuperficialis_1_2", "Mesh_Radius", MWMath::Point3D{-46.97, -0.906, -3.337} , "Mesh_MP1", MWMath::Point3D{2.393, 7.813, 2.173},  {"MeshVP_FDS", "Mesh_MC1", "MeshVPF_MC1", "Mesh_MCP1_Joint", "Mesh_PP1" , "MeshVPF_PP1", "Mesh_PIP1_Joint"}, { } }, 60, FLEXORCOLOR, 1.0 },
        { "FlexorDigitorumSuperficialis_2_1", "Mesh_Radius", MWMath::Point3D{-33.96, 1.06, -1.918} , "Mesh_MP2", MWMath::Point3D{2.047, 9.653, 0.5864},     {"MeshVP_FDS", "Mesh_MC2", "MeshVPF_MC2", "Mesh_MCP2_Joint", "Mesh_PP2" , "MeshVPF_PP2", "Mesh_PIP2_Joint"},  {}, 60, FLEXORCOLOR, 1.0 },
        //{ "FlexorDigitorumSuperficialis_2_2", "Mesh_Radius", MWMath::Point3D{-21.14, 3.75, -2.235} , "Mesh_MP2", MWMath::Point3D{2.234, 8.336, -0.335},   {"MeshVP_FDS", "Mesh_MC2", "MeshVPF_MC2", "Mesh_MCP2_Joint", "Mesh_PP2" , "MeshVPF_PP2", "Mesh_PIP2_Joint"}, {}, 60, FLEXORCOLOR, 1.0 },
        { "FlexorDigitorumSuperficialis_3_1", "Mesh_Humerus", MWMath::Point3D{-18.92, -330.2, -23.85} , "Mesh_MP3", MWMath::Point3D{2.652, 7.257, 2.47} ,   {"MeshVP_FDS", "Mesh_MC3", "MeshVPF_MC3", "Mesh_MCP3_Joint", "Mesh_PP3" , "MeshVPF_PP3", "Mesh_PIP3_Joint"}, {}, 60, FLEXORCOLOR, 1.0 },
        //{ "FlexorDigitorumSuperficialis_3_2", "Mesh_Humerus", MWMath::Point3D{-17.97, -338.3, -24.45} , "Mesh_MP3", MWMath::Point3D{2.03, 9.25, 1.878} ,  {"MeshVP_FDS", "Mesh_MC3", "MeshVPF_MC3", "Mesh_MCP3_Joint", "Mesh_PP3" , "MeshVPF_PP3", "Mesh_PIP3_Joint"}, {}, 60, FLEXORCOLOR, 1.0 },
        { "FlexorDigitorumSuperficialis_4_1", "Mesh_Humerus", MWMath::Point3D{-21.94, -368.0, -24.17} , "Mesh_MP4", MWMath::Point3D{1.495, 5.604, -1.007},  {"MeshVP_FDS", "Mesh_MC4", "MeshVPF_MC4", "Mesh_MCP4_Joint", "Mesh_PP4" , "MeshVPF_PP4", "Mesh_PIP4_Joint"}, {}, 60, FLEXORCOLOR, 1.0 },
        //{ "FlexorDigitorumSuperficialis_4_2", "Mesh_Humerus", MWMath::Point3D{-21.41, -370.6, -23.6} , "Mesh_MP4", MWMath::Point3D{1.457, 5.692, -1.092}, {"MeshVP_FDS", "Mesh_MC4", "MeshVPF_MC4", "Mesh_MCP4_Joint", "Mesh_PP4" , "MeshVPF_PP4", "Mesh_PIP4_Joint"}, {}, 60, FLEXORCOLOR, 1.0 },
       
        
        { "FlexorDigitorumProfundus_1_1", "Mesh_Ulna", MWMath::Point3D{20.25, 100.3, -17.71} , "Mesh_DP1", MWMath::Point3D{0.918, 1.99, 1.897} ,            {"MeshVP_FDP", "Mesh_MC1", "MeshVPF_MC1", "Mesh_MCP1_Joint", "Mesh_PP1" , "MeshVPF_PP1", "Mesh_PIP1_Joint", "Mesh_MP1", "MeshVPF_MP1", "Mesh_DIP1_Joint"}, {}, 50, FLEXORCOLOR, 1.0 },
        //{ "FlexorDigitorumProfundus_1_2", "Mesh_Ulna", MWMath::Point3D{18.29, 75.68, -11.29} , "Mesh_DP1", MWMath::Point3D{1.54, 3.546, 2.038} ,          {"MeshVP_FDP", "Mesh_MC1", "MeshVPF_MC1", "Mesh_MCP1_Joint", "Mesh_PP1" , "MeshVPF_PP1", "Mesh_PIP1_Joint", "Mesh_MP1", "MeshVPF_MP1", "Mesh_DIP1_Joint"}, {} }, 50, FLEXORCOLOR, 1.0 },
        { "FlexorDigitorumProfundus_2_1", "Mesh_Ulna", MWMath::Point3D{20.25, 100.3, -17.71} , "Mesh_DP2", MWMath::Point3D{0.918, 1.99, 1.897} ,            {"MeshVP_FDP", "Mesh_MC2", "MeshVPF_MC2", "Mesh_MCP2_Joint", "Mesh_PP2" , "MeshVPF_PP2", "Mesh_PIP2_Joint", "Mesh_MP2", "MeshVPF_MP2", "Mesh_DIP2_Joint"}, {}, 50, FLEXORCOLOR, 1.0 },
        //{ "FlexorDigitorumProfundus_2_2", "Mesh_Ulna", MWMath::Point3D{18.29, 75.68, -11.29} , "Mesh_DP2", MWMath::Point3D{1.54, 3.546, 2.038} ,          {"MeshVP_FDP", "Mesh_MC2", "MeshVPF_MC2", "Mesh_MCP2_Joint", "Mesh_PP2" , "MeshVPF_PP2", "Mesh_PIP2_Joint", "Mesh_MP2", "MeshVPF_MP2", "Mesh_DIP2_Joint"}, {} }, 50, FLEXORCOLOR, 1.0 },
        { "FlexorDigitorumProfundus_3_1", "Mesh_Ulna", MWMath::Point3D{20.25, 100.3, -17.71} , "Mesh_DP3", MWMath::Point3D{0.918, 1.99, 1.897} ,            {"MeshVP_FDP", "Mesh_MC3", "MeshVPF_MC3", "Mesh_MCP3_Joint", "Mesh_PP3" , "MeshVPF_PP3", "Mesh_PIP3_Joint", "Mesh_MP3", "MeshVPF_MP3", "Mesh_DIP3_Joint"}, {}, 50, FLEXORCOLOR, 1.0 },
        //{ "FlexorDigitorumProfundus_3_2", "Mesh_Ulna", MWMath::Point3D{18.29, 75.68, -11.29} , "Mesh_DP3", MWMath::Point3D{1.54, 3.546, 2.038} ,          {"MeshVP_FDP", "Mesh_MC3", "MeshVPF_MC3", "Mesh_MCP3_Joint", "Mesh_PP3" , "MeshVPF_PP3", "Mesh_PIP3_Joint", "Mesh_MP3", "MeshVPF_MP3", "Mesh_DIP3_Joint"}, {} }, 50, FLEXORCOLOR, 1.0 },
       
        { "FlexorDigitorumProfundus_4_1", "Mesh_Ulna", MWMath::Point3D{20.25, 100.3, -17.71} , "Mesh_DP4", MWMath::Point3D{0.918, 1.99, 1.897} ,            {"Mesh_Carpals", "MeshVP_FDP", "Mesh_MC4", "MeshVPF_MC4", "Mesh_MCP4_Joint", "Mesh_PP4" , "MeshVPF_PP4", "Mesh_PIP4_Joint", "Mesh_MP4", "MeshVPF_MP4", "Mesh_DIP4_Joint"}, { }, 50, FLEXORCOLOR, 1.0 },
        //{ "FlexorDigitorumProfundus_4_2", "Mesh_Ulna", MWMath::Point3D{18.29, 75.68, -11.29} , "Mesh_DP4", MWMath::Point3D{1.54, 3.546, 2.038} ,          {"MeshVP_FDP", "Mesh_MC4", "MeshVPF_MC4", "Mesh_MCP4_Joint", "Mesh_PP4" , "MeshVPF_PP4", "Mesh_PIP4_Joint", "Mesh_MP4", "MeshVPF_MP4", "Mesh_DIP4_Joint"}, {}, 50, FLEXORCOLOR, 1.0 },
        
        { "PalmarisLongus_1", "Mesh_Humerus", MWMath::Point3D{-18.52, -363.9, -24.39} , "Mesh_Carpals", MWMath::Point3D{0.295, -4.826, 0.209} , {"Mesh_Carpals", "MeshVP_FPL"}, {}, 50, FLEXORCOLOR },
       
        // ==========================================
        // --- EXTRINSIC EXTENSORS (Table 7) ---
        // ==========================================
        { "ExtensorCarpiRadialisLongus_1", "Mesh_Humerus", MWMath::Point3D{25.37, -320.2, 2.568} , "Mesh_MC1", MWMath::Point3D{-8.58, 28.16, 4.778} , {"MeshVP_ECRL"}, {}, 10, EXTENSORCOLOR, 1.0 },
        //{ "ExtensorCarpiRadialisLongus_2", "Mesh_Humerus", MWMath::Point3D{32.2, -342.3, 1.9} , "Mesh_MC1", MWMath::Point3D{-8.56, 25.58, 4.415} ,     {"MeshVP_ECRL"}, {}, 10, EXTENSORCOLOR, 1.0 },
        { "ExtensorCarpiRadialisBrevis_1", "Mesh_Humerus", MWMath::Point3D{39.0, -331.2, -1.775} , "Mesh_MC2", MWMath::Point3D{-12.26, 25.31, -0.861} ,  {"MeshVP_ECRL"}, {}, 10, EXTENSORCOLOR, 1.0 },
        //{ "ExtensorCarpiRadialisBrevis_2", "Mesh_Humerus", MWMath::Point3D{39.98, -336.8, -4.294} , "Mesh_MC2", MWMath::Point3D{-11.96, 27.73, -0.582},{"MeshVP_ECRL"}, {}, 10, EXTENSORCOLOR, 1.0 },
        { "ExtensorCarpiUlnaris_1", "Mesh_Humerus", MWMath::Point3D{24.11, -346.4, 0.571} , "Mesh_MC4", MWMath::Point3D{-2.222, 17.97, -4.251} ,         {"MeshVP_ECRL"}, {}, 10, EXTENSORCOLOR, 1.0 },
        //{ "ExtensorCarpiUlnaris_2", "Mesh_Humerus", MWMath::Point3D{30.96, -342.0, 2.139} , "Mesh_MC4", MWMath::Point3D{-2.792, 17.03, -4.285} ,       {"MeshVP_ECRL"}, {}, 10, EXTENSORCOLOR, 1.0 },
        

        { "ExtensorDigitorum_1_1", "Mesh_Humerus", MWMath::Point3D{34.87, -339.1, -3.738} , "Mesh_DP1", MWMath::Point3D{-4.47, 4.714, -1.288} ,{"Mesh_UlnaPseudo", "Mesh_Carpals", "MeshVP_ED", "MeshVPE_MC1", "Mesh_MCP1_Joint", "Mesh_PP1" , "MeshVPE_PP1", "Mesh_PIP1_Joint", "Mesh_MP1", "MeshVPE_MP1", "Mesh_DIP1_Joint"}, {}, 80, EXTENSORCOLOR, 1.0},
        //{ "ExtensorDigitorum_1_2", "Mesh_Humerus", MWMath::Point3D{30.49, -340.9, -4.776} , "Mesh_DP1", MWMath::Point3D{-4.524, 5.019, 0.458} , {"Mesh_UlnaPseudo", "Mesh_Carpals", "MeshVP_ED", "MeshVPE_MC1", "Mesh_MCP1_Joint", "Mesh_PP1" , "MeshVPE_PP1", "Mesh_PIP1_Joint", "Mesh_MP1", "MeshVPE_MP1", "Mesh_DIP1_Joint"}, {}, 80, EXTENSORCOLOR },
        { "ExtensorDigitorum_2_1", "Mesh_Humerus", MWMath::Point3D{34.87, -339.1, -3.738} , "Mesh_DP2", MWMath::Point3D{-4.47, 4.714, -1.288} ,   {"Mesh_UlnaPseudo", "Mesh_Carpals", "MeshVP_ED", "MeshVPE_MC2", "Mesh_MCP2_Joint", "Mesh_PP2" , "MeshVPE_PP2", "Mesh_PIP2_Joint", "Mesh_MP2", "MeshVPE_MP2", "Mesh_DIP2_Joint"}, {}, 80, EXTENSORCOLOR, 1.0},
        //{ "ExtensorDigitorum_2_2", "Mesh_Humerus", MWMath::Point3D{30.49, -340.9, -4.776} , "Mesh_DP2", MWMath::Point3D{-4.524, 5.019, 0.458} , {"Mesh_UlnaPseudo", "Mesh_Carpals", "MeshVP_ED", "MeshVPE_MC1", "Mesh_MCP1_Joint", "Mesh_PP1" , "MeshVPE_PP1", "Mesh_PIP1_Joint", "Mesh_MP1", "MeshVPE_MP1", "Mesh_DIP1_Joint"}, {}, 80, EXTENSORCOLOR },
        { "ExtensorDigitorum_3_1", "Mesh_Humerus", MWMath::Point3D{34.87, -339.1, -3.738} , "Mesh_DP3", MWMath::Point3D{-4.47, 4.714, -1.288} ,   {"Mesh_UlnaPseudo", "Mesh_Carpals", "MeshVP_ED", "MeshVPE_MC3", "Mesh_MCP3_Joint", "Mesh_PP3" , "MeshVPE_PP3", "Mesh_PIP3_Joint", "Mesh_MP3", "MeshVPE_MP3", "Mesh_DIP3_Joint"}, {}, 80, EXTENSORCOLOR, 1.0},
        //{ "ExtensorDigitorum_3_2", "Mesh_Humerus", MWMath::Point3D{30.49, -340.9, -4.776} , "Mesh_DP3", MWMath::Point3D{-4.524, 5.019, 0.458} , {"Mesh_UlnaPseudo", "Mesh_Carpals", "MeshVP_ED", "MeshVPE_MC1", "Mesh_MCP1_Joint", "Mesh_PP1" , "MeshVPE_PP1", "Mesh_PIP1_Joint", "Mesh_MP1", "MeshVPE_MP1", "Mesh_DIP1_Joint"}, {} }, 80, EXTENSORCOLOR },
        { "ExtensorDigitorum_4_1", "Mesh_Humerus", MWMath::Point3D{34.87, -339.1, -3.738} , "Mesh_DP4", MWMath::Point3D{-4.47, 4.714, -1.288} ,   {"Mesh_UlnaPseudo", "Mesh_Carpals", "MeshVP_ED", "MeshVPE_MC4", "Mesh_MCP4_Joint", "Mesh_PP4" , "MeshVPE_PP4", "Mesh_PIP4_Joint", "Mesh_MP4", "MeshVPE_MP4", "Mesh_DIP4_Joint"}, {}, 80, EXTENSORCOLOR, 1.0},
        //{ "ExtensorDigitorum_4_2", "Mesh_Humerus", MWMath::Point3D{30.49, -340.9, -4.776} , "Mesh_DP4", MWMath::Point3D{-4.524, 5.019, 0.458} , {"Mesh_UlnaPseudo", "Mesh_Carpals", "MeshVP_ED", "MeshVPE_MC1", "Mesh_MCP1_Joint", "Mesh_PP1" , "MeshVPE_PP1", "Mesh_PIP1_Joint", "Mesh_MP1", "MeshVPE_MP1", "Mesh_DIP1_Joint"}, {}, 80, EXTENSORCOLOR },
        

        { "ExtensorDigitiMinimi_1", "Mesh_Humerus", MWMath::Point3D{26.15, -371.1, -3.214} , "Mesh_DP4", MWMath::Point3D{-3.714, 8.692, -0.0449} ,{"Mesh_UlnaPseudo", "Mesh_Carpals", "MeshVP_ED", "MeshVPE_MC4", "Mesh_MCP4_Joint", "Mesh_PP4" , "MeshVPE_PP4", "Mesh_PIP4_Joint", "Mesh_MP4", "MeshVPE_MP4", "Mesh_DIP4_Joint"}, {}, 80, EXTENSORCOLOR, 1.0},
        //{ "ExtensorDigitiMinimi_2", "Mesh_Humerus", MWMath::Point3D{23.89, -371.0, -2.115} , "Mesh_DP4", MWMath::Point3D{-5.486, 4.598, -3.015} , {"Mesh_UlnaPseudo", "Mesh_Carpals", "MeshVP_ED", "MeshVPE_MC4", "Mesh_MCP4_Joint", "Mesh_PP4" , "MeshVPE_PP4", "Mesh_PIP4_Joint", "Mesh_MP4", "MeshVPE_MP4", "Mesh_DIP4_Joint"}, {}, 80, EXTENSORCOLOR },

        { "ExtensorIndicis_1", "Mesh_Ulna", MWMath::Point3D{10.83, 226.3, -3.307} , "Mesh_DP1", MWMath::Point3D{-4.342, 6.562, 1.834} , {"Mesh_UlnaPseudo", "Mesh_Carpals", "MeshVP_EDM", "MeshVPE_MC1", "Mesh_MCP1_Joint", "Mesh_PP1" , "MeshVPE_PP1", "Mesh_PIP1_Joint", "Mesh_MP1", "MeshVPE_MP1", "Mesh_DIP1_Joint"}, {}, 50, EXTENSORCOLOR, 1.0},
        //{ "ExtensorIndicis_2", "Mesh_Ulna", MWMath::Point3D{11.54, 210.5, -2.816} , "Mesh_DP1", MWMath::Point3D{-4.253, 5.784, 1.265} , {"Mesh_UlnaPseudo", "Mesh_Carpals", "MeshVP_EDM", "MeshVPE_MC1", "Mesh_MCP1_Joint", "Mesh_PP1" , "MeshVPE_PP1", "Mesh_PIP1_Joint", "Mesh_MP1", "MeshVPE_MP1", "Mesh_DIP1_Joint"}, {} }, 50, EXTENSORCOLOR }
         */
    };
    
    inline std::vector<MuscleDef> muscleDefsExcelOld = {
        /* 
        // ==========================================
        // --- THENAR (Daumen / Finger 0) ---
        // ==========================================
        { "AbductorPollicisBrevis_1", "Mesh_Carpals", MWMath::Point3D{8.9, 1.6, 14.0} , "Mesh_PP0", MWMath::Point3D{-2.9, 1.7, 3.5} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{4.0, 5.0, 18.0} } }, 3, OTHERMUSCLECOLOR },
        //{ "AbductorPollicisBrevis_2", "Mesh_Carpals", MWMath::Point3D{6.7, 0.2, 8.1} , "Mesh_PP0", MWMath::Point3D{-2.7, -0.028, 3.5} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{11.0, 5.0, 12.2} } }, 3, OTHERMUSCLECOLOR },
        //{ "AbductorPollicisBrevis_3", "Mesh_Carpals", MWMath::Point3D{7.4, 3.0, 14.0} , "Mesh_PP0", MWMath::Point3D{-3.4, 0.64, 3.3} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{9.0, 5.0, 13.0} } }, 3, OTHERMUSCLECOLOR },
        //{ "AbductorPollicisBrevis_4", "Mesh_Carpals", MWMath::Point3D{4.3, 3.7, 12.0} , "Mesh_PP0", MWMath::Point3D{-2.5, 0.42, 3.5} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{11.0, 5.0, 14.6} } }, 3, OTHERMUSCLECOLOR },
        //{ "AbductorPollicisBrevis_5", "Mesh_Carpals", MWMath::Point3D{7.1, 2.2, 10.0} , "Mesh_PP0", MWMath::Point3D{-2.3, 1.8, 2.8} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{6.6, 5.0, 15.0} } }, 3, OTHERMUSCLECOLOR },

        { "FlexorPollicisBrevisSuperficial_1", "Mesh_Carpals", MWMath::Point3D{8.1, 0.7, 17.0} , "Mesh_PP0", MWMath::Point3D{2.0, 6.1, 6.3} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{5.0, 7.5, 10.0} } }, 3, FLEXORCOLOR },
        //{ "FlexorPollicisBrevisSuperficial_2", "Mesh_Carpals", MWMath::Point3D{8.1, 0.6, 15.0} , "Mesh_PP0", MWMath::Point3D{2.0, 4.4, 5.3} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{8.0, 6.4, 8.0} } }, 3, FLEXORCOLOR },

        { "FlexorPollicisBrevisDeep_1", "Mesh_Carpals", MWMath::Point3D{3.02, -4.0, -1.0} , "Mesh_PP0", MWMath::Point3D{1.3, 0.5, 2.6} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{11.0, 6.5, 4.9} } }, 3, FLEXORCOLOR },

        { "OpponensPollicis_1", "Mesh_Carpals", MWMath::Point3D{9.1, -5.4, 19.0} , "Mesh_MC0", MWMath::Point3D{0.34, -5.6, 1.6} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{0.48, 8.8, 7.3} } }, 3, OTHERMUSCLECOLOR },
        //{ "OpponensPollicis_2", "Mesh_Carpals", MWMath::Point3D{9.7, -4.9, 15.0} , "Mesh_MC0", MWMath::Point3D{3.3, -0.4, 1.4} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{5.5, 8.3, 3.0} } }, 3, OTHERMUSCLECOLOR },
        //{ "OpponensPollicis_3", "Mesh_Carpals", MWMath::Point3D{8.3, -10.0, 18.0} , "Mesh_MC0", MWMath::Point3D{1.8, -1.1, 1.8} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{2.6, 8.5, 5.2} } }, 3, OTHERMUSCLECOLOR },
        //{ "OpponensPollicis_4", "Mesh_Carpals", MWMath::Point3D{6.2, -8.8, 25.0} , "Mesh_MC0", MWMath::Point3D{1.8, -10.5, 2.1} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{-2.1, 9.2, 9.8} } }, 3, OTHERMUSCLECOLOR },

        { "AdductorPollicisOblique_1", "Mesh_Carpals", MWMath::Point3D{6.7, -10.0, -0.4} , "Mesh_PP0", MWMath::Point3D{5.93, 6.8, -0.24} , {}, {}, 2, OTHERMUSCLECOLOR },
        //{ "AdductorPollicisOblique_2", "Mesh_Carpals", MWMath::Point3D{6.3, -14.0, -5.0} , "Mesh_PP0", MWMath::Point3D{5.77, 8.4, -0.91} , {}, {}, 2, OTHERMUSCLECOLOR },
        //{ "AdductorPollicisOblique_3", "Mesh_Carpals", MWMath::Point3D{6.82, -23.0, -4.8} , "Mesh_PP0", MWMath::Point3D{6.3, 10.0, -4.7} , {}, {}, 2, OTHERMUSCLECOLOR },

        { "AdductorPollicisTransverse_1", "Mesh_MC2", MWMath::Point3D{1.4, 23.0, -0.36} , "Mesh_PP0", MWMath::Point3D{3.7, 10.0, -7.5} , {}, {}, 2, OTHERMUSCLECOLOR },
        //{ "AdductorPollicisTransverse_2", "Mesh_MC2", MWMath::Point3D{1.7, -2.4, -0.078} , "Mesh_PP0", MWMath::Point3D{3.1, 4.9, -6.1} , {}, {}, 2, OTHERMUSCLECOLOR },
        //{ "AdductorPollicisTransverse_3", "Mesh_MC2", MWMath::Point3D{4.1, -10.0, -0.41} , "Mesh_PP0", MWMath::Point3D{1.2, 0.8, -5.0} , {}, {}, 2, OTHERMUSCLECOLOR },
        //{ "AdductorPollicisTransverse_4", "Mesh_MC2", MWMath::Point3D{1.3, 10.0, -0.4} , "Mesh_PP0", MWMath::Point3D{3.5, 7.5, -7.2} , {}, {}, 2, OTHERMUSCLECOLOR },

        // ==========================================
        // --- HYPOTHENAR (Kleiner Finger / Finger 4) ---
        // ==========================================
        { "FlexorDigitiMinimiBrevis_1", "Mesh_Carpals", MWMath::Point3D{12.0, -15.0, -21.0} , "Mesh_PP4", MWMath::Point3D{6.4, 11.0, -2.3} , {}, { {"Via_1", "Mesh_MC4", MWMath::Point3D{6.4, 1.9, -3.9} } }, 3, FLEXORCOLOR },

        { "AbductorDigitiMinimi_1", "Mesh_Carpals", MWMath::Point3D{12.0, -6.1, -25.0} , "Mesh_PP4", MWMath::Point3D{-0.631, 4.85, -4.56} , {}, { {"Via_1", "Mesh_MC4", MWMath::Point3D{2.1, 5.34, -8.5} } }, 3, OTHERMUSCLECOLOR },
        //{ "AbductorDigitiMinimi_2", "Mesh_Carpals", MWMath::Point3D{8.1, -5.6, -30.0} , "Mesh_PP4", MWMath::Point3D{1.25, 4.1, -5.85} , {}, { {"Via_1", "Mesh_MC4", MWMath::Point3D{6.14, 4.745, -11.2} } }, 3, OTHERMUSCLECOLOR },

        { "OpponensDigitiMinimi_1", "Mesh_Carpals", MWMath::Point3D{10.0, -22.0, -23.0} , "Mesh_MC4", MWMath::Point3D{2.3, -0.12, -4.0} , {}, { {"Via_1", "Mesh_MC4", MWMath::Point3D{3.0, 9.2, -3.2} } }, 3, OTHERMUSCLECOLOR },
        //{ "OpponensDigitiMinimi_2", "Mesh_Carpals", MWMath::Point3D{10.0, -20.0, -21.0} , "Mesh_MC4", MWMath::Point3D{2.1, -7.2, -4.0} , {}, { {"Via_1", "Mesh_MC4", MWMath::Point3D{4.7, 8.7, 0.84} } }, 3, OTHERMUSCLECOLOR },
        */
        
        /*
        // ==========================================
        // --- LUMBRICALS (Binnenmuskeln) ---
        // ==========================================
        { "Lumbricals1_1", "Mesh_MC1", MWMath::Point3D{10.0, 0.0, 0.0} , "Mesh_PP1", MWMath::Point3D{-2.0, 5.1, 3.7} ,  {"Mesh_MC1", "MeshMCP1_Joint"}, {}, 2, LUMBRICALCOLOR },
        { "Lumbricals2_1", "Mesh_MC2", MWMath::Point3D{10.0, 0.0, 0.0} , "Mesh_PP2", MWMath::Point3D{-0.36, 10.0, 4.1}, {"Mesh_MC2", "MeshMCP2_Joint"}, {}, 2, LUMBRICALCOLOR },
        { "Lumbricals3_1", "Mesh_MC2", MWMath::Point3D{10.0, 0.0, 0.0} , "Mesh_PP3", MWMath::Point3D{-3.4, 5.9, 3.5} ,  {"Mesh_MC2", "MeshMCP3_Joint"}, {}, 2, LUMBRICALCOLOR },
        { "Lumbricals3_2", "Mesh_MC3", MWMath::Point3D{10.0, 0.0, 0.0} , "Mesh_PP3", MWMath::Point3D{-2.1, 6.7, 4.1} ,  {"Mesh_MC3", "MeshMCP3_Joint"}, {}, 2, LUMBRICALCOLOR },
        { "Lumbricals4_1", "Mesh_MC3", MWMath::Point3D{10.0, 0.0, 0.0} , "Mesh_PP4", MWMath::Point3D{-0.038, 4.8, 3.2}, {"Mesh_MC3", "MeshMCP4_Joint"}, {}, 2, LUMBRICALCOLOR },
        { "Lumbricals4_2", "Mesh_MC4", MWMath::Point3D{10.0, 0.0, 0.0} , "Mesh_PP4", MWMath::Point3D{-0.6, 2.7, 2.7} ,  {"Mesh_MC4", "MeshMCP4_Joint"}, {}, 2, LUMBRICALCOLOR },
         
        // ==========================================
        // --- PALMAR INTEROSSEI ---
        // ==========================================
        { "PalmarInterossei1_1", "Mesh_MC1", MWMath::Point3D{0.2, -7.0, -2.6} , "Mesh_PP1", MWMath::Point3D{2.9, 8.1, -5.1} ,  {"MeshMCP1_Joint"}, { {"Via_1", "Mesh_MC1", MWMath::Point3D{2.2, -14.0, -6.1} } }, 3, PINTEROSSEICCOLOR },
        //{ "PalmarInterossei1_2", "Mesh_MC1", MWMath::Point3D{2.1, 11.0, -2.1} , "Mesh_PP1", MWMath::Point3D{1.51, 9.3, -5.6},{"MeshMCP1_Joint"}, { {"Via_1", "Mesh_MC1", MWMath::Point3D{2.0, -3.2, -7.5} } }, 3, PINTEROSSEICCOLOR },
        
        { "PalmarInterossei2_1", "Mesh_MC3", MWMath::Point3D{-0.53, 5.8, 1.6} , "Mesh_PP3", MWMath::Point3D{2.2, 12.0, 5.6} , {"MeshMCP3_Joint"}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{3.1, -0.49, 5.3} } }, 3, PINTEROSSEICCOLOR },
        //{ "PalmarInterossei2_2", "Mesh_MC3", MWMath::Point3D{-0.42, -5.3, 2.5} , "Mesh_PP3", MWMath::Point3D{1.2, 14.0, 6.5} , {"MeshMCP3_Joint"}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{1.6, -10.7, 5.3} } }, 3, PINTEROSSEICCOLOR },
        
        { "PalmarInterossei3_1", "Mesh_MC4", MWMath::Point3D{-2.0, -4.8, 3.9} , "Mesh_PP4", MWMath::Point3D{2.3, 7.4, 4.7} , {"MeshMCP4_Joint"}, { {"Via_1", "Mesh_MC4", MWMath::Point3D{0.0, -10.2, 5.2} } }, 3, PINTEROSSEICCOLOR },
        //{ "PalmarInterossei3_2", "Mesh_MC4", MWMath::Point3D{-2.0, 5.9, 2.1} , "Mesh_PP4", MWMath::Point3D{3.1, 7.1, 4.2} , {"MeshMCP4_Joint"}, { {"Via_1", "Mesh_MC4", MWMath::Point3D{-0.5, -0.98, 3.6} } }, 3, PINTEROSSEICCOLOR },

        // ==========================================
        // --- DORSAL INTEROSSEI ---
        // ==========================================
        { "DorsalInterossei1_1", "Mesh_MC1", MWMath::Point3D{-2.7, 14.0, 5.5} , "Mesh_PP1", MWMath::Point3D{-2.0, 12.0, 9.1} , {"MeshMCP1_Joint"}, { {"Via_1", "Mesh_MC1", MWMath::Point3D{0.81, 1.9, 8.5} } }, 3, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei1_2", "Mesh_MC1", MWMath::Point3D{-2.5, 1.4, 4.3} , "Mesh_PP1", MWMath::Point3D{-4.1, 12.0, 8.4} , {"MeshMCP1_Joint"}, { {"Via_1", "Mesh_MC1", MWMath::Point3D{-1.6, 1.1, 6.5} } }, 3, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei1_3", "Mesh_MC1", MWMath::Point3D{-2.1, -10.0, 4.5} , "Mesh_PP1", MWMath::Point3D{-2.8, 14.0, 9.3} , {"MeshMCP1_Joint"}, { {"Via_1", "Mesh_MC1", MWMath::Point3D{-1.1, 0.93, 11.0} } }, 3, DINTEROSSEICCOLOR },
        { "DorsalInterossei1_4", "Mesh_MC0", MWMath::Point3D{1.2, 5.7, -4.3} , "Mesh_PP1", MWMath::Point3D{-2.1, 11.0, 8.2} , {"MeshMCP1_Joint"}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{7.0, 8.0, -15.0} } }, 3, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei1_5", "Mesh_MC0", MWMath::Point3D{0.88, -1.3, -4.7} , "Mesh_PP1", MWMath::Point3D{-3.8, 11.0, 7.6} , {"MeshMCP1_Joint"}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{7.0, 0.0, -15.0} } }, 3, DINTEROSSEICCOLOR },

        { "DorsalInterossei2_1", "Mesh_MC1", MWMath::Point3D{-3.3, 13.0, -5.4} , "Mesh_PP2", MWMath::Point3D{-2.0, 10.0, 5.6} , {"MeshMCP2_Joint"}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-2.0, -9.0, 8.5} } }, 3, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei2_2", "Mesh_MC1", MWMath::Point3D{-3.7, 10.0, -4.7} , "Mesh_PP2", MWMath::Point3D{-2.6, 10.0, 5.6} , {"MeshMCP2_Joint"}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-2.6, -20.6, 8.5} } }, 3, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei2_3", "Mesh_MC1", MWMath::Point3D{-3.2, -9.5, -4.6} , "Mesh_PP2", MWMath::Point3D{-3.5, 10.0, 5.1} , {"MeshMCP2_Joint"}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-3.5, -29.6, 8.5} } }, 3, DINTEROSSEICCOLOR },
        { "DorsalInterossei2_4", "Mesh_MC2", MWMath::Point3D{-5.6, -9.2, 3.4} , "Mesh_PP2", MWMath::Point3D{-3.3, 15.0, 4.6} , {"MeshMCP2_Joint"}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-4.3, -26.0, 5.9} } }, 3, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei2_5", "Mesh_MC2", MWMath::Point3D{-5.9, 0.49, 3.4} , "Mesh_PP2", MWMath::Point3D{-0.96, 15.0, 4.6} , {"MeshMCP2_Joint"}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-3.6, -26.0, 5.4} } }, 3, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei2_6", "Mesh_MC2", MWMath::Point3D{-5.6, 9.9, 2.5} , "Mesh_PP2", MWMath::Point3D{-5.6, 15.0, 4.6} , {"MeshMCP2_Joint"}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-2.9, -26.0, 5.1} } }, 3, DINTEROSSEICCOLOR },
    
        { "DorsalInterossei3_1", "Mesh_MC2", MWMath::Point3D{-4.7, 9.8, -3.1} , "Mesh_PP2", MWMath::Point3D{-0.21, 16.0, -6.5} , {"MeshMCP2_Joint"}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-2.3, -21.0, -7.0} } }, 3, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei3_2", "Mesh_MC2", MWMath::Point3D{-4.2, 0.95, -3.7} , "Mesh_PP2", MWMath::Point3D{-2.9, 16.0, -5.0} , {"MeshMCP2_Joint"}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-3.0, -20.0, -6.5} } }, 3, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei3_3", "Mesh_MC2", MWMath::Point3D{-5.5, -8.8, -4.5} , "Mesh_PP2", MWMath::Point3D{-1.6, 17.0, -5.9} , {"MeshMCP2_Joint"}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-3.5, -21.0, -8.1} } }, 3, DINTEROSSEICCOLOR },
        { "DorsalInterossei3_4", "Mesh_MC3", MWMath::Point3D{-3.3, 7.7, 2.6} , "Mesh_PP2", MWMath::Point3D{-1.8, 10.0, -5.6} , {"MeshMCP2_Joint"}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-1.3, 20.0, -9.4} } }, 3, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei3_5", "Mesh_MC3", MWMath::Point3D{-2.8, -6.1, 3.2} , "Mesh_PP2", MWMath::Point3D{-2.7, 10.0, -4.9} , {"MeshMCP2_Joint"}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-2.8, -19.0, -9.9} } }, 3, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei3_6", "Mesh_MC3", MWMath::Point3D{-2.2, 0.96, 2.7} , "Mesh_PP2", MWMath::Point3D{-0.18, 10.0, -6.2} , {"MeshMCP2_Joint"}, { {"Via_1", "Mesh_MC2", MWMath::Point3D{-1.7, -20.0, -8.6} } }, 3, DINTEROSSEICCOLOR },

        { "DorsalInterossei4_1", "Mesh_MC3", MWMath::Point3D{-2.4, -0.64, -4.4} , "Mesh_PP3", MWMath::Point3D{-2.2, 17.0, -6.4} , {"MeshMCP3_Joint"}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{-1.9, -18.0, -7.1} } }, 3, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei4_2", "Mesh_MC3", MWMath::Point3D{-2.1, -6.7, -4.1} , "Mesh_PP3", MWMath::Point3D{-3.3, 17.0, -6.1} , {"MeshMCP3_Joint"}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{-1.8, -18.0, -5.7} } }, 3, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei4_3", "Mesh_MC3", MWMath::Point3D{-2.7, 6.0, -4.1} , "Mesh_PP3", MWMath::Point3D{-1.3, 17.0, -6.6} , {"MeshMCP3_Joint"}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{-0.56, -18.0, -6.0} } }, 3, DINTEROSSEICCOLOR },
        { "DorsalInterossei4_4", "Mesh_MC4", MWMath::Point3D{-4.2, 7.1, 1.1} , "Mesh_PP3", MWMath::Point3D{-1.6, 6.0, -6.0} , {"MeshMCP3_Joint"}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{-2.4, -16.0, -7.8} } }, 3, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei4_5", "Mesh_MC4", MWMath::Point3D{-4.8, 0.85, 1.8} , "Mesh_PP3", MWMath::Point3D{-1.1, 5.0, -6.0} , {"MeshMCP3_Joint"}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{-2.1, -18.0, -6.3} } }, 3, DINTEROSSEICCOLOR },
        //{ "DorsalInterossei4_6", "Mesh_MC4", MWMath::Point3D{-4.9, -6.4, 3.2} , "Mesh_PP3", MWMath::Point3D{-3.3, 5.0, -5.0} , {"MeshMCP3_Joint"}, { {"Via_1", "Mesh_MC3", MWMath::Point3D{-2.3, -16.0, -7.0} } }, 3, DINTEROSSEICCOLOR },
        */
       
        /* 
        // ==========================================
        // --- EXTRINSIC FLEXORS (Von Unterarm) ---
        // ==========================================
        { "FlexorCarpiUlnaris_1", "Mesh_Humerus", MWMath::Point3D{-17.5, -355.6, -20.3} , "Mesh_Carpals", MWMath::Point3D{10.77, 2.472, -20.86} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{16.64, 58.93, 16.68} } }, 3, FLEXORCOLOR },
        //{ "FlexorCarpiUlnaris_2", "Mesh_Humerus", MWMath::Point3D{-20.85, -359.3, -20.17} , "Mesh_Carpals", MWMath::Point3D{10.73, 3.284, -17.91} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{23.29, 58.17, 10.73} } }, 3, FLEXORCOLOR },
        
        { "FlexorCarpiRadialis_1", "Mesh_Humerus", MWMath::Point3D{-12.92, -359.6, -16.45} , "Mesh_MC1", MWMath::Point3D{2.84, 24.64, -4.447} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{14.98, 67.29, -21.02} } }, 3, FLEXORCOLOR },
        //{ "FlexorCarpiRadialis_2", "Mesh_Humerus", MWMath::Point3D{-13.64, -355.4, -13.84} , "Mesh_MC1", MWMath::Point3D{2.553, 28.05, -5.668} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{12.88, 67.91, -18.02} } }, 3, FLEXORCOLOR },

        { "FlexorPollicisLongus_1", "Mesh_Radius", MWMath::Point3D{-89.66, -8.361, -4.372} , "Mesh_DP0", MWMath::Point3D{0.997, 5.414, 3.404} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-182.0, -13.33, -10.5} } }, 3, FLEXORCOLOR },
        //{ "FlexorPollicisLongus_2", "Mesh_Radius", MWMath::Point3D{-77.92, -8.64, -2.982} , "Mesh_DP0", MWMath::Point3D{0.997, 5.414, 3.404} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-182.0, -8.82, -10.2} } }, 3, FLEXORCOLOR },
        
        { "FlexorDigitorumSuperficialis_3_1", "Mesh_Humerus", MWMath::Point3D{-18.92, -330.2, -23.85} , "Mesh_MP3", MWMath::Point3D{2.652, 7.257, 2.47} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{20.03, 99.78, -20.78} } }, 3, FLEXORCOLOR },
        //{ "FlexorDigitorumSuperficialis_3_2", "Mesh_Humerus", MWMath::Point3D{-17.97, -338.3, -24.45} , "Mesh_MP3", MWMath::Point3D{2.03, 9.25, 1.878} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{18.18, 98.96, -17.18} } }, 3, FLEXORCOLOR },
        { "FlexorDigitorumSuperficialis_4_1", "Mesh_Humerus", MWMath::Point3D{-21.94, -368.0, -24.17} , "Mesh_MP4", MWMath::Point3D{1.495, 5.604, -1.007} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{20.99, 99.75, -12.12} } }, 3, FLEXORCOLOR },
        //{ "FlexorDigitorumSuperficialis_4_2", "Mesh_Humerus", MWMath::Point3D{-21.41, -370.6, -23.6} , "Mesh_MP4", MWMath::Point3D{1.457, 5.692, -1.092} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{18.99, 99.53, -5.606} } }, 3, FLEXORCOLOR },
        
        { "FlexorDigitorumSuperficialis_1_1", "Mesh_Radius", MWMath::Point3D{-61.7, -3.448, -4.993} , "Mesh_MP1", MWMath::Point3D{2.116, 6.523, 2.057} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-73.54, 2.824, -9.478} } }, 3, FLEXORCOLOR },
        //{ "FlexorDigitorumSuperficialis_1_1", "Mesh_Radius", MWMath::Point3D{0.,0.,0.} , "Mesh_MP1", MWMath::Point3D{10.,0. ,0.} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-73.54, 2.824, -9.478} } }, 3, FLEXORCOLOR },
        
        //{ "FlexorDigitorumSuperficialis_1_2", "Mesh_Radius", MWMath::Point3D{-46.97, -0.906, -3.337} , "Mesh_MP1", MWMath::Point3D{2.393, 7.813, 2.173} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-74.02, 8.087, -13.73} } }, 3, FLEXORCOLOR },
        { "FlexorDigitorumSuperficialis_2_1", "Mesh_Radius", MWMath::Point3D{-33.96, 1.06, -1.918} , "Mesh_MP2", MWMath::Point3D{2.047, 9.653, 0.5864} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{12.48, 96.83, -25.0} } }, 3, FLEXORCOLOR },
        //{ "FlexorDigitorumSuperficialis_2_2", "Mesh_Radius", MWMath::Point3D{-21.14, 3.75, -2.235} , "Mesh_MP2", MWMath::Point3D{2.234, 8.336, -0.335} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{16.33, 98.35, -24.18} } }, 3, FLEXORCOLOR },
        
        { "FlexorDigitorumProfundus_1_1", "Mesh_Ulna", MWMath::Point3D{20.25, 100.3, -17.71} , "Mesh_DP1", MWMath::Point3D{0.918, 1.99, 1.897} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{19.01, 109.3, -16.66} } }, 3, FLEXORCOLOR },
        //{ "FlexorDigitorumProfundus_1_2", "Mesh_Ulna", MWMath::Point3D{18.29, 75.68, -11.29} , "Mesh_DP1", MWMath::Point3D{1.54, 3.546, 2.038} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{26.45, 110.3, -14.92} } }, 3, FLEXORCOLOR },
        { "FlexorDigitorumProfundus_2_1", "Mesh_Ulna", MWMath::Point3D{20.25, 100.3, -17.71} , "Mesh_DP2", MWMath::Point3D{0.918, 1.99, 1.897} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{19.01, 109.3, -16.66} } }, 3, FLEXORCOLOR },
        //{ "FlexorDigitorumProfundus_2_2", "Mesh_Ulna", MWMath::Point3D{18.29, 75.68, -11.29} , "Mesh_DP2", MWMath::Point3D{1.54, 3.546, 2.038} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{26.45, 110.3, -14.92} } }, 3, FLEXORCOLOR },
        { "FlexorDigitorumProfundus_3_1", "Mesh_Ulna", MWMath::Point3D{20.25, 100.3, -17.71} , "Mesh_DP3", MWMath::Point3D{0.918, 1.99, 1.897} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{19.01, 109.3, -16.66} } }, 3, FLEXORCOLOR },
        //{ "FlexorDigitorumProfundus_3_2", "Mesh_Ulna", MWMath::Point3D{18.29, 75.68, -11.29} , "Mesh_DP3", MWMath::Point3D{1.54, 3.546, 2.038} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{26.45, 110.3, -14.92} } }, 3, FLEXORCOLOR },
        { "FlexorDigitorumProfundus_4_1", "Mesh_Ulna", MWMath::Point3D{20.25, 100.3, -17.71} , "Mesh_DP4", MWMath::Point3D{0.918, 1.99, 1.897} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{19.01, 109.3, -16.66} } }, 3, FLEXORCOLOR },
        //{ "FlexorDigitorumProfundus_4_2", "Mesh_Ulna", MWMath::Point3D{18.29, 75.68, -11.29} , "Mesh_DP4", MWMath::Point3D{1.54, 3.546, 2.038} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{26.45, 110.3, -14.92} } }, 3, FLEXORCOLOR },

        { "PalmarisLongus_1", "Mesh_Humerus", MWMath::Point3D{-18.52, -363.9, -24.39} , "Mesh_Carpals", MWMath::Point3D{0.295, -4.826, 0.209} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{18.71, 75.68, -13.82} } }, 3, FLEXORCOLOR },

        // ==========================================
        // --- EXTRINSIC EXTENSORS (Table 7) ---
        // ==========================================
        { "ExtensorCarpiRadialisLongus_1", "Mesh_Humerus", MWMath::Point3D{25.37, -320.2, 2.568} , "Mesh_MC1", MWMath::Point3D{-8.58, 28.16, 4.778} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{4.335, -17.12, 2.269} } }, 3, EXTENSORCOLOR },
        //{ "ExtensorCarpiRadialisLongus_2", "Mesh_Humerus", MWMath::Point3D{32.2, -342.3, 1.9} , "Mesh_MC1", MWMath::Point3D{-8.56, 25.58, 4.415} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-3.796, -20.34, 10.88} } }, 3, EXTENSORCOLOR },
        { "ExtensorCarpiRadialisBrevis_1", "Mesh_Humerus", MWMath::Point3D{39.0, -331.2, -1.775} , "Mesh_MC2", MWMath::Point3D{-12.26, 25.31, -0.861} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-33.14, -9.632, 16.33} } }, 3, EXTENSORCOLOR },
        //{ "ExtensorCarpiRadialisBrevis_2", "Mesh_Humerus", MWMath::Point3D{39.98, -336.8, -4.294} , "Mesh_MC2", MWMath::Point3D{-11.96, 27.73, -0.582} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-34.24, -16.86, 11.73} } }, 3, EXTENSORCOLOR },
        { "ExtensorCarpiUlnaris_1", "Mesh_Humerus", MWMath::Point3D{24.11, -346.4, 0.571} , "Mesh_MC4", MWMath::Point3D{-2.222, 17.97, -4.251} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{-4.28, 109.4, 4.36} } }, 3, EXTENSORCOLOR },
        //{ "ExtensorCarpiUlnaris_2", "Mesh_Humerus", MWMath::Point3D{30.96, -342.0, 2.139} , "Mesh_MC4", MWMath::Point3D{-2.792, 17.03, -4.285} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{-2.785, 109.5, 2.853} } }, 3, EXTENSORCOLOR },


        { "ExtensorDigitorum_1_1", "Mesh_Humerus", MWMath::Point3D{34.87, -339.1, -3.738} , "Mesh_DP1", MWMath::Point3D{-4.47, 4.714, -1.288} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.48, -11.22, 13.74} } }, 3, EXTENSORCOLOR },
        //{ "ExtensorDigitorum_1_2", "Mesh_Humerus", MWMath::Point3D{30.49, -340.9, -4.776} , "Mesh_DP1", MWMath::Point3D{-4.524, 5.019, 0.458} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.21, -4.973, 14.83} } }, 3, EXTENSORCOLOR },
        { "ExtensorDigitorum_2_1", "Mesh_Humerus", MWMath::Point3D{34.87, -339.1, -3.738} , "Mesh_DP2", MWMath::Point3D{-4.47, 4.714, -1.288} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.48, -11.22, 13.74} } }, 3, EXTENSORCOLOR },
        //{ "ExtensorDigitorum_2_2", "Mesh_Humerus", MWMath::Point3D{30.49, -340.9, -4.776} , "Mesh_DP2", MWMath::Point3D{-4.524, 5.019, 0.458} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.21, -4.973, 14.83} } }, 3, EXTENSORCOLOR },
        { "ExtensorDigitorum_3_1", "Mesh_Humerus", MWMath::Point3D{34.87, -339.1, -3.738} , "Mesh_DP3", MWMath::Point3D{-4.47, 4.714, -1.288} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.48, -11.22, 13.74} } }, 3, EXTENSORCOLOR },
        //{ "ExtensorDigitorum_3_2", "Mesh_Humerus", MWMath::Point3D{30.49, -340.9, -4.776} , "Mesh_DP3", MWMath::Point3D{-4.524, 5.019, 0.458} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.21, -4.973, 14.83} } }, 3, EXTENSORCOLOR },
        { "ExtensorDigitorum_4_1", "Mesh_Humerus", MWMath::Point3D{34.87, -339.1, -3.738} , "Mesh_DP4", MWMath::Point3D{-4.47, 4.714, -1.288} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.48, -11.22, 13.74} } }, 3, EXTENSORCOLOR },
        //{ "ExtensorDigitorum_4_2", "Mesh_Humerus", MWMath::Point3D{30.49, -340.9, -4.776} , "Mesh_DP4", MWMath::Point3D{-4.524, 5.019, 0.458} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.21, -4.973, 14.83} } }, 3, EXTENSORCOLOR },


        { "ExtensorDigitiMinimi_1", "Mesh_Humerus", MWMath::Point3D{26.15, -371.1, -3.214} , "Mesh_DP4", MWMath::Point3D{-3.714, 8.692, -0.0449} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{-15.26, 121.7, 7.688} } }, 3, EXTENSORCOLOR },
        //{ "ExtensorDigitiMinimi_2", "Mesh_Humerus", MWMath::Point3D{23.89, -371.0, -2.115} , "Mesh_DP4", MWMath::Point3D{-5.486, 4.598, -3.015} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{-0.121, 121.8, 6.031} } }, 3, EXTENSORCOLOR },

        { "ExtensorIndicis_1", "Mesh_Ulna", MWMath::Point3D{10.83, 226.3, -3.307} , "Mesh_DP1", MWMath::Point3D{-4.342, 6.562, 1.834} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-219.5, -10.35, 25.02} } }, 3, EXTENSORCOLOR },
        //{ "ExtensorIndicis_2", "Mesh_Ulna", MWMath::Point3D{11.54, 210.5, -2.816} , "Mesh_DP1", MWMath::Point3D{-4.253, 5.784, 1.265} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-219.8, -5.307, 22.61} } }, 3, EXTENSORCOLOR }
        */
    };
    

    inline std::vector<MuscleDef> muscleDefsExcelReduced = {
        
        // ==========================================
        // --- THENAR (Daumen / Finger 0) ---
        // ==========================================
        { "AbductorPollicisBrevis_1", "Mesh_Carpals", MWMath::Point3D{8.9, 1.6, 14.0} , "Mesh_PP0", MWMath::Point3D{-2.9, 1.7, 3.5} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{4.0, 5.0, 18.0} } }, 3, OTHERMUSCLECOLOR },

        { "FlexorPollicisBrevisSuperficial_1", "Mesh_Carpals", MWMath::Point3D{8.1, 0.7, 17.0} , "Mesh_PP0", MWMath::Point3D{2.0, 6.1, 6.3} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{5.0, 7.5, 10.0} } }, 3, FLEXORCOLOR },

        { "FlexorPollicisBrevisDeep_1", "Mesh_Carpals", MWMath::Point3D{3.02, -4.0, -1.0} , "Mesh_PP0", MWMath::Point3D{1.3, 0.5, 2.6} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{11.0, 6.5, 4.9} } }, 3, FLEXORCOLOR },

        { "OpponensPollicis_1", "Mesh_Carpals", MWMath::Point3D{9.1, -5.4, 19.0} , "Mesh_MC0", MWMath::Point3D{0.34, -5.6, 1.6} , {}, { {"Via_1", "Mesh_MC0", MWMath::Point3D{0.48, 8.8, 7.3} } }, 3, OTHERMUSCLECOLOR },

        { "AdductorPollicisOblique_1", "Mesh_Carpals", MWMath::Point3D{6.7, -10.0, -0.4} , "Mesh_PP0", MWMath::Point3D{5.93, 6.8, -0.24} , {}, {}, 2, OTHERMUSCLECOLOR },

        { "AdductorPollicisTransverse_1", "Mesh_MC2", MWMath::Point3D{1.4, 23.0, -0.36} , "Mesh_PP0", MWMath::Point3D{3.7, 10.0, -7.5} , {}, {}, 2, OTHERMUSCLECOLOR },
        
        // ==========================================
        // --- HYPOTHENAR (Kleiner Finger / Finger 4) ---
        // ==========================================
        { "FlexorDigitiMinimiBrevis_1", "Mesh_Carpals", MWMath::Point3D{12.0, -15.0, -21.0} , "Mesh_PP4", MWMath::Point3D{6.4, 11.0, -2.3} , {}, { {"Via_1", "Mesh_MC4", MWMath::Point3D{6.4, 1.9, -3.9} } }, 3, FLEXORCOLOR },

        { "AbductorDigitiMinimi_1", "Mesh_Carpals", MWMath::Point3D{12.0, -6.1, -25.0} , "Mesh_PP4", MWMath::Point3D{-0.631, 4.85, -4.56} , {}, { {"Via_1", "Mesh_MC4", MWMath::Point3D{2.1, 5.34, -8.5} } }, 3, OTHERMUSCLECOLOR },

        { "OpponensDigitiMinimi_1", "Mesh_Carpals", MWMath::Point3D{10.0, -22.0, -23.0} , "Mesh_MC4", MWMath::Point3D{2.3, -0.12, -4.0} , {}, { {"Via_1", "Mesh_MC4", MWMath::Point3D{3.0, 9.2, -3.2} } }, 3, OTHERMUSCLECOLOR },
        
        // ==========================================
        // --- EXTRINSIC FLEXORS (Von Unterarm) ---
        // ==========================================
        { "FlexorCarpiUlnaris_1", "Mesh_Humerus", MWMath::Point3D{-17.5, -355.6, -20.3} , "Mesh_Carpals", MWMath::Point3D{10.77, 2.472, -20.86} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{16.64, 58.93, 16.68} } }, 3, FLEXORCOLOR },
        
        { "FlexorCarpiRadialis_1", "Mesh_Humerus", MWMath::Point3D{-12.92, -359.6, -16.45} , "Mesh_MC1", MWMath::Point3D{2.84, 24.64, -4.447} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{14.98, 67.29, -21.02} } }, 3, FLEXORCOLOR },

        { "FlexorPollicisLongus_1", "Mesh_Radius", MWMath::Point3D{-89.66, -8.361, -4.372} , "Mesh_DP0", MWMath::Point3D{0.997, 5.414, 3.404} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-182.0, -13.33, -10.5} } }, 3, FLEXORCOLOR },
        
        { "FlexorDigitorumSuperficialis_1_1", "Mesh_Radius", MWMath::Point3D{-61.7, -3.448, -4.993} , "Mesh_MP1", MWMath::Point3D{2.116, 6.523, 2.057} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-73.54, 2.824, -9.478} } }, 3, FLEXORCOLOR },

        { "FlexorDigitorumProfundus_1_1", "Mesh_Ulna", MWMath::Point3D{20.25, 100.3, -17.71} , "Mesh_DP1", MWMath::Point3D{0.918, 1.99, 1.897} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{19.01, 109.3, -16.66} } }, 3, FLEXORCOLOR },
        
        { "FlexorDigitorumProfundus_2_1", "Mesh_Ulna", MWMath::Point3D{20.25, 100.3, -17.71} , "Mesh_DP2", MWMath::Point3D{0.918, 1.99, 1.897} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{19.01, 109.3, -16.66} } }, 3, FLEXORCOLOR },
        
        { "FlexorDigitorumProfundus_3_1", "Mesh_Ulna", MWMath::Point3D{20.25, 100.3, -17.71} , "Mesh_DP3", MWMath::Point3D{0.918, 1.99, 1.897} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{19.01, 109.3, -16.66} } }, 3, FLEXORCOLOR },
        
        { "FlexorDigitorumProfundus_4_1", "Mesh_Ulna", MWMath::Point3D{20.25, 100.3, -17.71} , "Mesh_DP4", MWMath::Point3D{0.918, 1.99, 1.897} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{19.01, 109.3, -16.66} } }, 3, FLEXORCOLOR },

        { "PalmarisLongus_1", "Mesh_Humerus", MWMath::Point3D{-18.52, -363.9, -24.39} , "Mesh_Carpals", MWMath::Point3D{0.295, -4.826, 0.209} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{18.71, 75.68, -13.82} } }, 3, FLEXORCOLOR },

        // ==========================================
        // --- EXTRINSIC EXTENSORS (Table 7) ---
        // ==========================================
        { "ExtensorCarpiRadialisLongus_1", "Mesh_Humerus", MWMath::Point3D{25.37, -320.2, 2.568} , "Mesh_MC1", MWMath::Point3D{-8.58, 28.16, 4.778} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{4.335, -17.12, 2.269} } }, 3, EXTENSORCOLOR },

        { "ExtensorCarpiRadialisBrevis_1", "Mesh_Humerus", MWMath::Point3D{39.0, -331.2, -1.775} , "Mesh_MC2", MWMath::Point3D{-12.26, 25.31, -0.861} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-33.14, -9.632, 16.33} } }, 3, EXTENSORCOLOR },

        { "ExtensorCarpiUlnaris_1", "Mesh_Humerus", MWMath::Point3D{24.11, -346.4, 0.571} , "Mesh_MC4", MWMath::Point3D{-2.222, 17.97, -4.251} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{-4.28, 109.4, 4.36} } }, 3, EXTENSORCOLOR },

        { "ExtensorDigitorum_1_1", "Mesh_Humerus", MWMath::Point3D{34.87, -339.1, -3.738} , "Mesh_DP1", MWMath::Point3D{-4.47, 4.714, -1.288} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.48, -11.22, 13.74} } }, 3, EXTENSORCOLOR },
        
        { "ExtensorDigitorum_2_1", "Mesh_Humerus", MWMath::Point3D{34.87, -339.1, -3.738} , "Mesh_DP2", MWMath::Point3D{-4.47, 4.714, -1.288} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.48, -11.22, 13.74} } }, 3, EXTENSORCOLOR },
        
        { "ExtensorDigitorum_3_1", "Mesh_Humerus", MWMath::Point3D{34.87, -339.1, -3.738} , "Mesh_DP3", MWMath::Point3D{-4.47, 4.714, -1.288} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.48, -11.22, 13.74} } }, 3, EXTENSORCOLOR },
        
        { "ExtensorDigitorum_4_1", "Mesh_Humerus", MWMath::Point3D{34.87, -339.1, -3.738} , "Mesh_DP4", MWMath::Point3D{-4.47, 4.714, -1.288} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-89.48, -11.22, 13.74} } }, 3, EXTENSORCOLOR },

        { "ExtensorDigitiMinimi_1", "Mesh_Humerus", MWMath::Point3D{26.15, -371.1, -3.214} , "Mesh_DP4", MWMath::Point3D{-3.714, 8.692, -0.0449} , {}, { {"Via_1", "Mesh_Ulna", MWMath::Point3D{-15.26, 121.7, 7.688} } }, 3, EXTENSORCOLOR },

        { "ExtensorIndicis_1", "Mesh_Ulna", MWMath::Point3D{10.83, 226.3, -3.307} , "Mesh_DP1", MWMath::Point3D{-4.342, 6.562, 1.834} , {}, { {"Via_1", "Mesh_Radius", MWMath::Point3D{-219.5, -10.35, 25.02} } }, 3, EXTENSORCOLOR }
    
    };
    
    
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
        // vgl Prometheus book
        // --- Wrist ---
        double W_FE_m = -60; // Flexion/Extension Min (Extension) -> Dorsalextension
        double W_FE_M = 80; // Flexion/Extension Max (Flexion) -> Palmarflexion
        double W_AA_m = -20; // Abduction/Adduction Min (Adduction) -> Radialabduktion
        double W_AA_M = 40; // Abduction/Adduction Max (Abduction) -> Ulnarabduktion

        // --- CMC (Carpometacarpal) ---
        //                     [0]     [1]     [2]    [3]    [4]
        double CMC_FE_m[5] = { -0.0,   0.0,   0.0,  -0.0, -0.0 }; // Flexion/Extension Min
        double CMC_FE_M[5] = {  0.0,   5.0,   5.0,  0.0,  0.0 }; // Flexion/Extension Max
        double CMC_AA_m[5] = { -0.0,   0.0,   0.0,   0.0,   0.0 }; // Abduction/Adduction Min
        double CMC_AA_M[5] = {  0.0,   0.0,   0.0,   0.0,   0.0 }; // Abduction/Adduction Max

        // --- MCP (Metacarpophalangeal / Grundgelenk) ---
        //                     [0]     [1]     [2]    [3]    [4] 
        double MCP_FE_m[5] = {-45.0, -45.0, -45.0, -45.0, -45.0 }; // Flex/Ext Min -> Extension/Dorsal
        double MCP_FE_M[5] = { 90.0, 90.0, 90.0, 90.0, 90.0 }; // Flex/Ext Max -> Flexion/Palmar
        double MCP_AA_m[5] = {-10.0, -20.0, -15.0, -15.0, -25.0 }; // Abd/Add Min -> direction away middlefinger
        double MCP_AA_M[5] = { 10.0,  20.0,  15.0,  15.0,  25.0 }; // Abd/Add Max -> direction to middlefinger

        // --- PIP (Proximal Interphalangeal / Mittelgelenk) ---
        //                     [0]     [1]     [2]    [3]    [4]  
        double PIP_m[5] = { -15.0,  15.0,  15.0,  15.0,  15.0 }; // Reines Scharnier
        double PIP_M[5] = { 100.0, 100.0, 100.0, 100.0, 100.0 }; 

        // --- DIP (Distal Interphalangeal / Endgelenk) ---
        //                     [0]     [1]     [2]    [3]    [4] 
        double DIP_m[5] = { 0.0,  -10.0, -10.0, -10.0, -10.0}; // Reines Scharnier
        double DIP_M[5] = { 0.0,   90.0,  90.0,  90.0,  90.0 };
    };

    // Konstanten aus AnyScript
    constexpr double SegLenRatios[5][4] = {
        {0.118, 0.251, 0.196, 0.158}, 
        {0.463, 0.245, 0.143, 0.097}, 
        {0.446, 0.266, 0.170, 0.108}, 
        {0.421, 0.244, 0.165, 0.107}, 
        {0.414, 0.204, 0.117, 0.093}  
    };

    constexpr double SegLenRatiosEngelhardt[5][4] = {
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
}

// wk10
inline std::string createMusclePathHand(std::shared_ptr<SSBody>& rootSystem, std::vector<std::shared_ptr<SSMesh>>& meshes, std::vector<SSMuscle*>& muscles, const SimSettings& cfg, int numNodes, double scale = 1.0f){
    qDebug() << "    [createMusclePathHand]...";
    std::vector<MuscleDef> allMusclesMira = Hand::muscleDefsMira;
    std::vector<MuscleDef> allMusclesExcel = Hand::muscleDefsExcel;
    
    double scaleFMira = 0.01;
    double scaleFExcel = 0.01;
    bool bUseViaPoints = false;

    std::vector<std::string> createdMuscleNames;

    for (const auto& mDef : allMusclesExcel) {
        qDebug() << "Processing muscle: " << mDef.name.c_str();
        // get meshes to consider/wrap
        std::vector<SSMesh*> pathMeshes;
        for (const auto& meshName : mDef.meshes) {
            if (getRefMeshFromList(meshName, meshes)) {
                pathMeshes.push_back(getRefMeshFromList(meshName, meshes));
            }
        }

        const MuscleDef muscleD = mDef;
        const MuscleDef insertM = mDef;
        SSMesh* oriMesh = getRefMeshFromList(muscleD.originMeshName, meshes);
        SSMesh* insMesh = getRefMeshFromList(muscleD.insertionMeshName, meshes);
        

        // check if muscle parents exists
        bool bSkipMuscle = true;
        if (insMesh != nullptr || insertM.insertionMeshName == "Mesh_Radius" || insertM.insertionMeshName == "Mesh_Ulna" || insertM.insertionMeshName == "Mesh_Humerus") {
            if (oriMesh != nullptr || muscleD.originMeshName == "Mesh_Radius" || muscleD.originMeshName == "Mesh_Ulna" || muscleD.originMeshName == "Mesh_Humerus") {
                bSkipMuscle = false;
            } else {
                std::cerr << "     [WARNING] Origin Mesh " << insertM.originMeshName << " nicht gefunden für Muskel " << mDef.name << ". Muskel wird übersprungen." << std::endl;
                bSkipMuscle = true;
            }
        }
        if (bSkipMuscle) continue;

        /* SSTissue* oriBody = oriMesh ? oriMesh->Parent.get() : rootSystem.get();
        MWMath::Point3D oriPos = oriBody->PositionGlobal + (Hand::WristToUlnaMM + Hand::RotMiraToWorld * originM.originRelPos) * scaleFMira;
        SSTissue* insBody = insMesh ? insMesh->Parent.get() : rootSystem.get();
        MWMath::Point3D insPos = insBody->PositionGlobal + insertM.insertionRelPos * scaleFExcel; */

        MWMath::Point3D boneOffset(0.f,0.f,0.f);
        MWMath::RotMatrix3x3 boneRotOffset;
        bool bIsExtrinsic = false;
        if (muscleD.originMeshName == "Mesh_Humerus") {
            boneOffset = Hand::WristToHumerusMM;
            boneRotOffset = Hand::RotCarpalsToHumerus;
            bIsExtrinsic = true;
        } 
        else if (muscleD.originMeshName == "Mesh_Radius") {
            boneOffset = Hand::WristToRadiusMM;
            boneRotOffset = Hand::RotCarpalsToRadius;
            bIsExtrinsic = true;
        } 
        else if (muscleD.originMeshName == "Mesh_Ulna") {
            boneOffset = Hand::WristToUlnaMM;
            boneRotOffset = Hand::RotCarpalsToUlna;
            bIsExtrinsic = true;
        }

        SSTissue* oriBody;
        MWMath::Point3D oriPos;
        SSTissue* insBody;
        MWMath::Point3D insPos;
        if (bIsExtrinsic){
            oriBody = (dynamic_cast<SSMesh*>(oriMesh)) ? oriMesh->Parent.get() : rootSystem.get();
            oriPos = oriBody->PositionGlobal + (Hand::RotWristToWorld * boneOffset * scaleFMira) + (Hand::RotWristToWorld * boneRotOffset * muscleD.originRelPos * scaleFMira);
            insBody = (dynamic_cast<SSMesh*>(insMesh)) ? insMesh->Parent.get() : rootSystem.get(); //dynamic_cast<SSBody*>(insMesh->Parent.get()); //s (dynamic_cast<SSMesh*>(insMesh)) ? insMesh->Parent.get() : rootSystem.get();
            insPos = muscleD.insertionRelPos * scaleFExcel + MWMath::Point3D(0., -dynamic_cast<SSEllipsoidMesh*>(insMesh)->C, 0.);  
        }
        else{
            oriBody = (oriMesh != nullptr) ? dynamic_cast<SSTissue*>(oriMesh) : rootSystem.get(); // (dynamic_cast<SSMesh*>(oriMesh)) ? oriMesh->Parent.get() : rootSystem.get();
            qDebug() << "  [ORIGIN] Parent Mesh: " << muscleD.originMeshName.c_str() << "->" << oriBody->Name.c_str() << " (Mesh Ptr: " << oriBody << ")";
            oriPos = muscleD.originRelPos * scaleFExcel; // + muscleD.originRelPos * scaleFExcel;// + (Hand::RotWristToWorld * boneOffset * scaleFMira) + (Hand::RotWristToWorld * boneRotOffset * muscleD.originRelPos * scaleFMira);
            insBody = (insMesh != nullptr) ? dynamic_cast<SSTissue*>(insMesh) : rootSystem.get(); //(dynamic_cast<SSMesh*>(insMesh)) ? insMesh->Parent.get() : rootSystem.get(); //dynamic_cast<SSBody*>(insMesh->Parent.get()); //s (dynamic_cast<SSMesh*>(insMesh)) ? insMesh->Parent.get() : rootSystem.get();
            qDebug() << "  [INSERT] Parent Mesh: " << muscleD.insertionMeshName.c_str() << "->" << insBody->Name.c_str() << " (Mesh Ptr: " << insBody << ")";
            insPos = muscleD.insertionRelPos * scaleFExcel;
            //MWMath::Point3D insPos = muscleD.insertionRelPos * scaleFExcel + MWMath::Point3D(0., -dynamic_cast<SSEllipsoidMesh*>(insMesh)->C, 0.);
        }
        

        
        // ---------------------------------------------------------
        // DEBUG LOGGING: Positionen überprüfen
        // ---------------------------------------------------------
        /* qDebug() << "==================================================";
        qDebug() << "MUSCLE: " << mDef.name.c_str();
        
        // Origin Logging
        std::string oName = oriBody ? oriBody->Name : "NULL";
        qDebug() << "  [ORIGIN] Parent Mesh: " << muscleD.originMeshName.c_str() << " (Body: " << oName.c_str() << ")";
        qDebug() << "  [ORIGIN] Parent PosGlobal:  X:" << oriBody->PositionGlobal.x << " Y:" << oriBody->PositionGlobal.y << " Z:" << oriBody->PositionGlobal.z;
        qDebug() << "  [ORIGIN] Calculated oriPos: X:" << oriPos.x << " Y:" << oriPos.y << " Z:" << oriPos.z;
        qDebug() << "  [Muscle] Calculated relInsert: X:" << oriPos.x << " Y:" << oriPos.y << " Z:" << oriPos.z;
        
        // Insertion Logging
        std::string iName = insBody ? insBody->Name : "NULL";
        qDebug() << "  [INSERT] Parent Mesh: " << mDef.insertionMeshName.c_str() << " (Body: " << iName.c_str() << ")";
        qDebug() << "  [INSERT] Parent PosGlobal:  X:" << insBody->PositionGlobal.x << " Y:" << insBody->PositionGlobal.y << " Z:" << insBody->PositionGlobal.z;
        qDebug() << "  [INSERT] Calculated insPos: X:" << insPos.x << " Y:"  << insPos.y << " Z:" << insPos.z;
        qDebug() << "==================================================";*/
        // ---------------------------------------------------------
        //int nodesCount = cfg.muscleNumPoints.empty() ? 2 : cfg.muscleNumPoints[0];
        int nodesCount = mDef.numNodes;
        SSMuscle* muscle = new SSMuscle(mDef.name, nodesCount, oriBody, oriPos, insBody, insPos);
        muscle->meshPtrs = pathMeshes;
        muscle->MeshColor = muscleD.mcolor;
        createdMuscleNames.push_back(mDef.name);

        muscle->TorusPathDirection = mDef.torusPathDirection;
        muscle->createMusclePointsComplexPath();
        muscle->updateMusclePointsParentsLocal();
        muscles.push_back(muscle);
    }

    return "Created muscles: " + std::to_string(createdMuscleNames.size());
}

inline MWMath::Point3D getRotatedPosition(MWMath::Point3D center, double R, double angleDeg) {
    double angleRad = angleDeg * M_PI / 180.0;
    double newX = R * std::cos(angleRad);
    double newY = center.y + R * std::sin(angleRad);
    return MWMath::Point3D(newX, newY, 0.0);
}



// thumb
inline void buildThumbModel( std::vector<std::shared_ptr<SSTissue>>& tissues, std::vector<std::shared_ptr<SSMesh>>& meshes, std::unordered_map<std::string, SSMesh*>& MeshMap, std::shared_ptr<SSBody> carpals, int numTimeSteps, const std::vector<double>& FJAngles, double scale, double Sign, double HL, double HB, const MWMath::Point3D& colorThumb, double rWF = 0.5, float GS = 1.)
{
    // ==============================================================================
    // 0. INIT & THUMB PARAMETERS
    // ==============================================================================
    std::string prefN = "0"; // Thumb = Finger 0
    double viaPointR = 0.05;
    double off = 1.0;
    
    // Daumen ist meist etwas dicker als der Zeigefinger
    std::vector<double> width = {0.16*GS*rWF, 0.12*GS*rWF, 0.08*GS*rWF}; 
    std::vector<double> relTorusPos = {0.7, 0.55, 0.45};

    // Knochenlängen: Der Daumen hat nur 3 Glieder! 
    // Nutze die Indizes 0 (MC), 1 (PP), und 3 (DP, in AnyBody oft so gemappt, da MP fehlt)
    double L1 = Hand::SegLenRatios[0][1] * HL *0.5; // Metacarpal
    double L2 = Hand::SegLenRatios[0][2] * HL *0.5; // Proximal
    double L3 = Hand::SegLenRatios[0][3] * HL *0.5; // Distal (IP)

    // ==============================================================================
    // 1. CMC GELENK & BASIS-ROTATION (Sattelgelenk)
    // ==============================================================================
    // 1. Position des Daumens: Wir nutzen Hand::JointPosTable[0]
    // In AnyBody ist JointPosTable[0] -> X=0, Y=-0.073*HL, Z=0.196*HB
    MWMath::Point3D posCMC(
        0.0, 
        Hand::JointPosTable[0][0] * HL, 
        Hand::JointPosTable[0][1] * HB * Sign
    );
    
    // 2. Rotation des Daumens: Wir nutzen deine definierten ThumbRotX und ThumbRotY Konstanten
    // Der Daumen wird erst um X gedreht, dann um Y in die Handfläche geklappt.
    MWMath::RotMatrix3x3 rotX = MWMath::axisAngle({1, 0, 0}, Hand::ThumbRotX * Sign);
    MWMath::RotMatrix3x3 rotY = MWMath::axisAngle({0, 1, 0}, Hand::ThumbRotY * Sign);
    
    MWMath::RotMatrix3x3 forientation = rotX * rotY; 

    // CMC Joint Abd (X-Achse)
    auto jointCMCAbd = std::make_shared<SSJoint>("CMC0_Joint", posCMC, forientation, carpals, 0.0, MWMath::Point3D(1, 0, 0), numTimeSteps);
    tissues.push_back(jointCMCAbd);
    
    auto jMeshCMC = std::make_shared<SSEllipsoidMesh>(width[0], width[0], width[0], "Mesh_CMC0", jointCMCAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0,0 }, 90.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMeshCMC);
    MeshMap["Mesh_CMC0_Joint"] = jMeshCMC.get();

    // CMC Joint Flex (Z-Achse)
    auto jointCMCFlex = std::make_shared<SSJoint>("CMC0_JointFlex", MWMath::Point3D(0,0,0), MWMath::RotMatrix3x3(), jointCMCAbd, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(jointCMCFlex);

    // ==============================================================================
    // 2. SEGMENT 1 (Metacarpal)
    // ==============================================================================
    auto body1 = std::make_shared<SSBody>("MC0", MWMath::Point3D(0, -L1, 0), MWMath::RotMatrix3x3(), jointCMCFlex);
    tissues.push_back(body1);
    
    auto mesh1 = std::make_shared<SSEllipsoidMesh>(width[0], width[0], L1, "Mesh_MC0", body1, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), colorThumb);
    meshes.push_back(mesh1);
    MeshMap["Mesh_MC0"] = mesh1.get();

    // ViaPoint Torus
    auto mVP1 = std::make_shared<SSEllipsoidMesh>(viaPointR, viaPointR, viaPointR, "MeshVPF_MC0", body1, MWMath::Point3D(width[0]*off, -L1 * relTorusPos[0], 0.), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(1, 0, 0));
    mVP1->bIsViaPoint = true;
    meshes.push_back(mVP1);
    MeshMap["MeshVPF_MC0"] = mVP1.get();

    // ==============================================================================
    // 3. MCP GELENK (Flexion & Abduktion nach AnyBody)
    // ==============================================================================
    MWMath::Point3D jointPosRel1(0, -L1, 0);
    
    // MCP Flexion (Z-Achse)
    auto jointMCPFlex = std::make_shared<SSJoint>("MCP0_JointFlex", jointPosRel1, MWMath::RotMatrix3x3(), body1, FJAngles[0], MWMath::Point3D(0,0,1), numTimeSteps);
    tissues.push_back(jointMCPFlex);
    
    double jointSize1 = width[0]*1.1;
    auto jMeshMCP = std::make_shared<SSEllipsoidMesh>(jointSize1, jointSize1, jointSize1, "Mesh_MCP0_Joint", jointMCPFlex, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMeshMCP);
    MeshMap["Mesh_MCP0_Joint"] = jMeshMCP.get();
    
    // MCP Abduktion (X-Achse)
    auto jointMCPAbd = std::make_shared<SSJoint>("MCP0_JointAbd", MWMath::Point3D(0,0,0), MWMath::axisAngle({1,0,0}, 0.0), jointMCPFlex, FJAngles[1], MWMath::Point3D(1,0,0), numTimeSteps);
    tissues.push_back(jointMCPAbd);

    // ==============================================================================
    // 4. SEGMENT 2 (Proximal Phalanx)
    // ==============================================================================
    auto body2 = std::make_shared<SSBody>("PP0", MWMath::Point3D(0, -L2, 0), MWMath::RotMatrix3x3(), jointMCPAbd);
    tissues.push_back(body2);
    
    auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2, "Mesh_PP0", body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), colorThumb);
    meshes.push_back(mesh2);
    MeshMap["Mesh_PP0"] = mesh2.get();

    auto mVP2 = std::make_shared<SSEllipsoidMesh>(viaPointR, viaPointR, viaPointR, "MeshVPF_PP0", body2, MWMath::Point3D(width[1]*off, -L2 * relTorusPos[1], 0.), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(1, 0, 0));
    mVP2->bIsViaPoint = true;
    meshes.push_back(mVP2);
    MeshMap["MeshVPF_PP0"] = mVP2.get();

    auto mVPE1 = std::make_shared<SSEllipsoidMesh>(viaPointR, viaPointR, viaPointR, "MeshVPE_PP0", jointMCPAbd, MWMath::Point3D(-jointSize1*1.2, 0., 0.), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0, 0, 1));
    mVPE1->bIsViaPoint = true;
    meshes.push_back(mVPE1);
    MeshMap["MeshVPE_PP0"] = mVPE1.get();

    // ==============================================================================
    // 5. IP GELENK (Interphalangeal - In AnyBody als DIP betitelt)
    // ==============================================================================
    MWMath::Point3D jointPosRel2(0, -L2, 0);
    auto jointIP = std::make_shared<SSJoint>("IP0_JointFlex", jointPosRel2, MWMath::RotMatrix3x3(), body2, FJAngles[2], MWMath::Point3D(0,0,1), numTimeSteps);
    tissues.push_back(jointIP);
    
    double jointSize2 = width[1]*1.1;
    auto jMeshIP = std::make_shared<SSEllipsoidMesh>(jointSize2, jointSize2, jointSize2, "Mesh_IP0_Joint", jointIP, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMeshIP);
    MeshMap["Mesh_IP0_Joint"] = jMeshIP.get();

    auto mVPE2 = std::make_shared<SSEllipsoidMesh>(viaPointR, viaPointR, viaPointR, "MeshVPE_DP0", jointIP, MWMath::Point3D(-jointSize2*1.2, 0., 0.), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0, 0, 1));
    mVPE2->bIsViaPoint = true;
    meshes.push_back(mVPE2);
    MeshMap["MeshVPE_DP0"] = mVPE2.get();

    // ==============================================================================
    // 6. SEGMENT 3 (Distal Phalanx)
    // ==============================================================================
    auto body3 = std::make_shared<SSBody>("DP0", MWMath::Point3D(0, -L3, 0), MWMath::RotMatrix3x3(), jointIP);
    tissues.push_back(body3);
    
    auto mesh3 = std::make_shared<SSEllipsoidMesh>(width[2], width[2], L3, "Mesh_DP0", body3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), colorThumb);
    meshes.push_back(mesh3);
    MeshMap["Mesh_DP0"] = mesh3.get();
}


inline void buildWristViaPoints(std::vector<std::shared_ptr<SSMesh>>& meshes, std::unordered_map<std::string, SSMesh*>& MeshMap, std::shared_ptr<SSBody> wristBody, double scale = 1.0, double wristWidth = 6.4, double viaPointR = 0.03) {

    if (!wristBody) return;
    double angleDeg = 5.0;
    MWMath::RotMatrix3x3 planeRotation = MWMath::axisAngle(MWMath::Point3D(1, 0, 0), angleDeg);

    
    struct VPDef {
        std::string name;
        double x;
        double z;
        int isThumb;
    };

    // Offsets (y = 0 in der gedrehten Ebene)
    std::vector<VPDef> vpDefs = {
        // in cm
        {"MeshVP_EDM",  -1.8,  1.8, 0},
        {"MeshVP_ED",  -1.8,  0.0, 0},
        {"MeshVP_ECRB", -1.8, -1.8, 0},
        {"MeshVP_EPL",  -1.8, -2.5, 2},
        {"MeshVP_ECRL", -1.0, -3.0, 0},
        {"MeshVP_EPB",   0.8, -4.0, 2},
        {"MeshVP_APL",   1.5, -4.0, 2},
        {"MeshVP_FCR",   2.5, -1.8, 1},
        {"MeshVP_FCU",   2.3,  2.5, 1},
        {"MeshVP_ADMB",  1.0,  3.3, 3},
        {"MeshVP_ECU",  -1.0,  2.8, 0},
        {"MeshVP_FPL",   1.5,  -1., 1},
        {"MeshVP_FDP",   1.5,  0.0, 1},
        {"MeshVP_FDS",   2.0,  0.0, 1}
    };
    

    for (const auto& def : vpDefs) {
        
        // Punkt in der ungedrehten lokalen Ebene (Y=0) konstruieren (mit Skalierung)
        MWMath::Point3D baseOffset(def.x*.01f * scale, 0.0, -def.z*.01f * scale);
        // Punkt in die um 10° gedrehte Ebene projizieren
        MWMath::Point3D rotatedOffset = planeRotation * baseOffset;
        
        MWMath::Point3D COLOR_VP;
        if (def.isThumb==1){ COLOR_VP = MWMath::Point3D(0.6, 0.1, 0.1);}
        else if (def.isThumb==0) {COLOR_VP = MWMath::Point3D(0.1, 0.1, 0.6);}
        else if (def.isThumb==2) {COLOR_VP = MWMath::Point3D(1.0, 0.5, 0.0);}
        else if (def.isThumb==3) {COLOR_VP = MWMath::Point3D(0.1, 1.0, 0.0);}

        // Mesh erstellen (an wristBody)
        auto vpMesh = std::make_shared<SSEllipsoidMesh>(
            viaPointR, viaPointR, viaPointR, 
            def.name, 
            wristBody, 
            rotatedOffset, 
            MWMath::RotMatrix3x3(), // Orientierung an die Ebene anpassen
            COLOR_VP
        );
        vpMesh->bIsViaPoint = true;
        vpMesh->MViaPointTolerance = viaPointR*1.4; 

        meshes.push_back(vpMesh);
        MeshMap[def.name] = vpMesh.get();
    }
}




// DONT TOUCH THIS FUNCTION, IT WORKS FOR NOW!!!!
inline std::string buildOHandModelOldExpandedDontTouch(std::vector<std::shared_ptr<SSTissue>>& tissues,std::vector<std::shared_ptr<SSMesh>>& meshes,std::vector<SSMuscle*>& muscles,std::shared_ptr<SSBody>& rootSystem,int numTimeSteps,const SimSettings& cfg,float geometryScaler,const std::vector<double> processParams)
{
    // -----------------------------------
    // This function represents the original "OFINGER_SIMPLE_OnlyTorusSmall" Case, where most of the 90°-Rotations 
    // worked, but now as a handmodel builder script.
    // -----------------------------------

    meshes.clear();
    muscles.clear();

    std::vector<double> FJAngles = {90.0, 0.0, 100.0, 90.0}; // MCP1, MCP2, PIP, DIP 
    MWMath::Point3D COLORF1 = MWMath::Point3D(0.50, 0.00, 1.00); 
    MWMath::Point3D COLORF2 = MWMath::Point3D(0.30, 0.30, 1.00); 
    MWMath::Point3D COLORF3 = MWMath::Point3D(0.00, 0.50, 1.00); 
    MWMath::Point3D COLORF4 = MWMath::Point3D(0.00, 0.90, 1.00);

    MWMath::Point3D COLORJOINT = MWMath::Point3D(0.0, 0.9, 0.3);
    MWMath::Point3D COLORFDEACTIVE = MWMath::Point3D(0.5, 0.5, 0.5);

    // --- PARAMETER ---
    double HL = 1.8;
    double HB = 0.8;
    double scale = HL / Hand::RefHandLength;
    std::unordered_map<std::string, SSMesh*> MeshMap;
    double Sign = 1.0; // Für Linke Hand -1, Rechte Hand +1
    float GS = 1.0;
    std::vector<double> fLength = {Hand::SegLenRatios[1][0] * HL, Hand::SegLenRatios[1][1] * HL, Hand::SegLenRatios[1][2] * HL, Hand::SegLenRatios[1][3] * HL};
    //std::vector<double> fLength = {0.3482 * HL, 0.2027 * HL, 0.1175 * HL, 0.0882 * HL};
    double L1 = fLength[0]; 
    double L2 = fLength[1]; 
    double L3 = fLength[2]; 
    double L4 = fLength[3]; 
    double rWF = 1.0; 
    std::vector<double> width = {0.15*GS*rWF, 0.1*GS*rWF, 0.07*GS*rWF, 0.05*GS*rWF}; 
    std::vector<double> relTorusPos = {0.6, 0.42, 0.32};
    std::vector<double> relTorusR = {1.3/rWF, 1.3/rWF, 1.6/rWF}; 
    std::vector<double> relTorusr = {1.1/rWF, 1.1/rWF, 1.4/rWF}; 

    bool bShowBody = true;
    double off = 1.0; 

    // 1. ROOT
    MWMath::RotMatrix3x3 Ausrichtung = MWMath::axisAngle({0,0,1}, 90.0) * MWMath::axisAngle({0,1,0}, 180.0);
    rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0, 0, 0), Ausrichtung, nullptr);
    tissues.push_back(rootSystem);

    MWMath::Point3D wristOffset(0.0, 0.003567 * scale, -0.003901 * scale * Sign);
    
    auto wristJointAbd = std::make_shared<SSJoint>("Wrist_JointAbd", wristOffset, MWMath::RotMatrix3x3(), rootSystem, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointAbd);
    /* auto jMeshWrist = std::make_shared<SSEllipsoidMesh>(0.01 * scale, 0.01 * scale, 0.01 * scale, "Mesh_Wrist_Joint", wristJointAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0,0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMeshWrist);
    MeshMap["Mesh_Wrist_Joint"] = jMeshWrist.get(); */

    auto wristJointFlex = std::make_shared<SSJoint>("Wrist_JointFlex", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointAbd, 0.0, MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointFlex);

    // 2. CARPALS
    auto carpals = std::make_shared<SSBody>("Carpals", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointFlex);
    tissues.push_back(carpals);
    
    double carpalThickness = width[0] * 0.5; 
    double carpalLength = width[0] * 1.25;    
    double carpalWidth = width[0] * 0.5 * 3.5;     
    
    auto meshCarpals = std::make_shared<SSEllipsoidMesh>(carpalThickness, carpalLength, carpalWidth, "Mesh_Carpals", carpals, MWMath::Point3D(0, 0.0*scale, 0), MWMath::RotMatrix3x3(), MWMath::Point3D(0.5, 0.5, 0.5));
    //meshes.push_back(meshCarpals);
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
    // SEGMENT 1 (Proximal)
    // ==============================================================================
    auto body1 = std::make_shared<SSBody>("Body_Seg1", MWMath::Point3D(0,-L1,0), MWMath::RotMatrix3x3(), jointCMCFlex);
    tissues.push_back(body1);
    
    // Mesh 1: Rotiert um X, damit es entlang -Y zeigt
    auto mesh1 = std::make_shared<SSEllipsoidMesh>(width[0], width[0], L1, "Mesh_Seg1", body1, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), COLORFDEACTIVE);
    if (bShowBody) meshes.push_back(mesh1);
    
    // Torus 1: Offset auf +X, Länge auf -Y
    auto mTorus1 = std::make_shared<SSTorusMesh>(mesh1->B * relTorusR[0], mesh1->B * relTorusr[0], 
            "Torus1", body1, MWMath::Point3D(width[0]*off, -L1 * relTorusPos[0], 0.), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(1, 0, 0));
    meshes.push_back(mTorus1);

    // Joint 1: Position am Ende von Segment 1 auf der Y-Achse
    MWMath::Point3D jointPosRel1 = MWMath::Point3D(0, -L1, 0);
    auto joint = std::make_shared<SSJoint>("Joint_1", 
                                            jointPosRel1,           // Position relativ zum Parent
                                            MWMath::RotMatrix3x3(), // Initiale Rotation
                                            body1,                  // Parent Body
                                            FJAngles[0],            // Max Winkel
                                            MWMath::Point3D(0,0,1),// Drehachse: Bleibt Z!
                                            numTimeSteps);
    tissues.push_back(joint);
    
    double jointSize1 = width[0]*1.1/rWF;
    // Gelenkkugel braucht keine spezifische Rotation
    auto jMesh1 = std::make_shared<SSEllipsoidMesh>(jointSize1, jointSize1, jointSize1, "Mesh_Joint1", joint, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMesh1);
    
    MWMath::Point3D jointPosRel11 = MWMath::Point3D(0, 0, 0);
    auto joint1 = std::make_shared<SSJoint>("Joint_11", 
                                            jointPosRel11, 
                                            MWMath::axisAngle({1,0,0}, 0.0),
                                            joint, 
                                            FJAngles[1], 
                                            MWMath::Point3D(1,0,0), // Abduktion wandert von Y auf X!
                                            numTimeSteps);
    tissues.push_back(joint1);
    

    // ==============================================================================
    // SEGMENT 2 (Distal)
    // ==============================================================================
    // Body 2 startet am Ende des vorherigen Segments
    auto body2 = std::make_shared<SSBody>("Body_Seg2", MWMath::Point3D(0,-L2,0), MWMath::RotMatrix3x3(), joint1);
    tissues.push_back(body2);
    
    auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2, 
            "Mesh_Seg2", body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), COLORFDEACTIVE);
    if (bShowBody) meshes.push_back(mesh2);
    
    // Torus 2
    auto mTorus2 = std::make_shared<SSTorusMesh>(mesh2->B * relTorusR[1], mesh2->B * relTorusr[1], 
            "Torus2", body2, MWMath::Point3D(width[1]*off, -L2 * relTorusPos[1], 0.), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
    meshes.push_back(mTorus2);

    MWMath::Point3D jointPosRel2 = MWMath::Point3D(0, -L2, 0);
    auto joint2 = std::make_shared<SSJoint>("Joint_2", 
                                            jointPosRel2, 
                                            MWMath::RotMatrix3x3(),
                                            body2, 
                                            FJAngles[2], 
                                            MWMath::Point3D(0,0,1), // Bleibt Z
                                            numTimeSteps);
    tissues.push_back(joint2);
    
    double jointSize2 = width[1]*1.1/rWF;
    auto jMesh2 = std::make_shared<SSEllipsoidMesh>(jointSize2, jointSize2, jointSize2, "Mesh_Joint2", joint2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMesh2);

    // ==============================================================================
    // SEGMENT 3
    // ==============================================================================
    auto body3 = std::make_shared<SSBody>("Body_Seg3", MWMath::Point3D(0,-L3,0), MWMath::RotMatrix3x3(), joint2);
    tissues.push_back(body3);
    
    auto mesh3 = std::make_shared<SSEllipsoidMesh>(width[2], width[2], L3, 
            "Mesh_Seg3", body3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), 
            COLORFDEACTIVE);
    if (bShowBody) meshes.push_back(mesh3);
    
    auto mTorus3 = std::make_shared<SSTorusMesh>(width[2] * relTorusR[2], width[2] * relTorusr[2], 
            "Torus3", body3, MWMath::Point3D(width[2]*off, -L3 * relTorusPos[2], 0.), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
    meshes.push_back(mTorus3);

    MWMath::Point3D jointPosRel3 = MWMath::Point3D(0, -L3, 0);
    auto joint3 = std::make_shared<SSJoint>("Joint_3", 
                                            jointPosRel3, 
                                            MWMath::RotMatrix3x3(),
                                            body3, 
                                            FJAngles[3], 
                                            MWMath::Point3D(0,0,1), // Bleibt Z
                                            numTimeSteps);
    tissues.push_back(joint3);

    double jointSize3 = width[2]*1.1/rWF;
    auto jMesh3 = std::make_shared<SSEllipsoidMesh>(jointSize3, jointSize3, jointSize3, "Mesh_Joint3", joint3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMesh3);

    // ==============================================================================
    // SEGMENT 4
    // ==============================================================================
    auto body4 = std::make_shared<SSBody>("Body_Seg4", MWMath::Point3D(0,-L4,0), MWMath::RotMatrix3x3(), joint3);
    tissues.push_back(body4);
    
    auto mesh4 = std::make_shared<SSEllipsoidMesh>(width[3], width[3], L4, 
            "Mesh_Seg4", body4, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), 
            COLORFDEACTIVE);
    if (bShowBody) meshes.push_back(mesh4);

    // ==============================================================================
    // INITIALES UPDATE
    // ==============================================================================
    for (auto& m : meshes) {m->InitializeMesh();}
    rootSystem->update(0); 

    // ==============================================================================
    // MUSKEL
    // ==============================================================================
    int numPoints = cfg.muscleNumPoints[0];
    
    // Muskelpunkte ebenfalls an neue Achsen angepasst: 
    // Offset auf der X-Achse, Länge (C) auf der negativen Y-Achse
    MWMath::Point3D startOffset = MWMath::Point3D(HL*0.03,HL*0.3, 0.0); //MWMath::Point3D(mesh1->B * 1.1, 0.0, 0.0);
    MWMath::Point3D endOffset = MWMath::Point3D(mesh4->B * 1.1, -mesh4->C, 0.0);
    
    SSMuscle* flexor = new SSMuscle("Flexor", numPoints, 
        rootSystem.get(), startOffset, 
        body4.get(), endOffset);

    // Alle Hindernisse
    for(auto& m : meshes) {
        /* if (dynamic_cast<SSTorusMesh*>(m.get()) || m->bIsJointMesh){
            flexor->meshPtrs.push_back(m.get());
        } */
        flexor->meshPtrs.push_back(m.get());
    }

    flexor->createMusclePointsComplexPath();
    flexor->updateMusclePointsParents();
    muscles.push_back(flexor);

    qDebug ()<< "GEHT BIS DA";
    return "buildOHandModelOldExpanded";
}

inline std::string buildOHandModelOldExpanded(std::vector<std::shared_ptr<SSTissue>>& tissues,std::vector<std::shared_ptr<SSMesh>>& meshes,std::vector<SSMuscle*>& muscles,std::shared_ptr<SSBody>& rootSystem,int numTimeSteps,const SimSettings& cfg,float geometryScaler,const std::vector<double> processParams)
{
    // -----------------------------------
    // This function represents the original "OFINGER_SIMPLE_OnlyTorusSmall" Case, where most of the 90°-Rotations 
    // worked, but now as a handmodel builder script.
    // -----------------------------------

    meshes.clear();
    muscles.clear();

    std::vector<double> FWAngles = {0.0, 0.0};
    std::vector<double> FJAngles = {90.0, 0.0,  100.0, 80.0}; // MCP1, MCP2, PIP, DIP 
    MWMath::Point3D COLORF1 = MWMath::Point3D(0.50, 0.00, 1.00); 
    MWMath::Point3D COLORF2 = MWMath::Point3D(0.30, 0.30, 1.00); 
    MWMath::Point3D COLORF3 = MWMath::Point3D(0.00, 0.50, 1.00); 
    MWMath::Point3D COLORF4 = MWMath::Point3D(0.00, 0.90, 1.00);

    MWMath::Point3D COLORJOINT = MWMath::Point3D(0.0, 0.9, 0.3);
    MWMath::Point3D COLORFDEACTIVE = MWMath::Point3D(0.5, 0.5, 0.5);

    // --- PARAMETER ---
    double HL = 1.8;
    double HB = 0.8;
    double scale = HL / Hand::RefHandLength;
    std::unordered_map<std::string, SSMesh*> MeshMap;
    double Sign = 1.0; // Für Linke Hand -1, Rechte Hand +1
    float GS = 1.0;
    std::vector<double> fLength = {Hand::SegLenRatios[1][0] * HL, Hand::SegLenRatios[1][1] * HL, Hand::SegLenRatios[1][2] * HL, Hand::SegLenRatios[1][3] * HL};
    //std::vector<double> fLength = {0.3482 * HL, 0.2027 * HL, 0.1175 * HL, 0.0882 * HL};
    double L1 = fLength[0]; 
    double L2 = fLength[1]; 
    double L3 = fLength[2]; 
    double L4 = fLength[3]; 
    double rWF = 1.0; 
    std::vector<double> width = {0.15*GS*rWF, 0.1*GS*rWF, 0.07*GS*rWF, 0.05*GS*rWF}; 
    std::vector<double> relTorusPos = {0.6, 0.42, 0.32};
    std::vector<double> relTorusR = {1.3/rWF, 1.3/rWF, 1.6/rWF}; 
    std::vector<double> relTorusr = {1.1/rWF, 1.1/rWF, 1.4/rWF}; 

    bool bShowBody = true;
    double off = 1.0; 

    // 1. ROOT
    MWMath::RotMatrix3x3 Ausrichtung = MWMath::axisAngle({0,0,1}, 90.0) * MWMath::axisAngle({0,1,0}, 180.0);
    rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0, 0, 0), Ausrichtung, nullptr);
    tissues.push_back(rootSystem);

    MWMath::Point3D wristOffset(0.0, 0.003567 * scale, -0.003901 * scale * Sign);
    
    auto wristJointAbd = std::make_shared<SSJoint>("Wrist_JointAbd", wristOffset, MWMath::RotMatrix3x3(), rootSystem, FWAngles[0], MWMath::Point3D(1, 0, 0), numTimeSteps);
    tissues.push_back(wristJointAbd);
    /* auto jMeshWrist = std::make_shared<SSEllipsoidMesh>(0.01 * scale, 0.01 * scale, 0.01 * scale, "Mesh_Wrist_Joint", wristJointAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0,0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMeshWrist);
    MeshMap["Mesh_Wrist_Joint"] = jMeshWrist.get(); */

    auto wristJointFlex = std::make_shared<SSJoint>("Wrist_JointFlex", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointAbd, FWAngles[1], MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointFlex);

    // 2. CARPALS
    auto carpals = std::make_shared<SSBody>("Carpals", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointFlex);
    tissues.push_back(carpals);
    
    double carpalThickness = width[0] * 0.5; 
    double carpalLength = width[0] * 1.25;    
    double carpalWidth = width[0] * 0.5 * 3.5;     
    
    auto meshCarpals = std::make_shared<SSEllipsoidMesh>(carpalThickness, carpalLength, carpalWidth, "Mesh_Carpals", carpals, MWMath::Point3D(0, 0.0*scale, 0), MWMath::RotMatrix3x3(), MWMath::Point3D(0.5, 0.5, 0.5));
    //meshes.push_back(meshCarpals);
    MeshMap["Mesh_Carpals"] = meshCarpals.get();

    // 3. WRAPPING SURFACE (Karpaltunnel)
    double flexCylRadius = 0.0115 * scale;
    double flexCylLength = 0.1400 * scale; 
    MWMath::Point3D flexCylOffset(carpalThickness * 0.5 + flexCylRadius, 0.02 * scale, 0.0);

    auto flexorCylMesh = std::make_shared<SSCylinderMesh>(flexCylRadius, flexCylLength, "Mesh_WristFlexorCyl", rootSystem, flexCylOffset, MWMath::axisAngle({ 1, 0, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));  
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
    // SEGMENT 1 (Proximal)
    // ==============================================================================
    auto body1 = std::make_shared<SSBody>("Body_Seg1", MWMath::Point3D(0,-L1,0), MWMath::RotMatrix3x3(), jointCMCFlex);
    tissues.push_back(body1);
    
    // Mesh 1: Rotiert um X, damit es entlang -Y zeigt
    auto mesh1 = std::make_shared<SSEllipsoidMesh>(width[0], width[0], L1, "Mesh_Seg1", body1, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), COLORFDEACTIVE);
    if (bShowBody) meshes.push_back(mesh1);
    
    // Torus 1: Offset auf +X, Länge auf -Y
    auto mTorus1 = std::make_shared<SSTorusMesh>(mesh1->B * relTorusR[0], mesh1->B * relTorusr[0], 
            "Torus1", body1, MWMath::Point3D(width[0]*off, -L1 * relTorusPos[0], 0.), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(1, 0, 0));
    meshes.push_back(mTorus1);

    // Joint 1: Position am Ende von Segment 1 auf der Y-Achse
    MWMath::Point3D jointPosRel1 = MWMath::Point3D(0, -L1, 0);
    auto joint = std::make_shared<SSJoint>("Joint_1", 
                                            jointPosRel1,           // Position relativ zum Parent
                                            MWMath::RotMatrix3x3(), // Initiale Rotation
                                            body1,                  // Parent Body
                                            FJAngles[0],            // Max Winkel
                                            MWMath::Point3D(0,0,1),// Drehachse: Bleibt Z!
                                            numTimeSteps);
    tissues.push_back(joint);
    
    double jointSize1 = width[0]*1.1/rWF;
    // Gelenkkugel braucht keine spezifische Rotation
    auto jMesh1 = std::make_shared<SSEllipsoidMesh>(jointSize1, jointSize1, jointSize1, "Mesh_Joint1", joint, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMesh1);
    
    MWMath::Point3D jointPosRel11 = MWMath::Point3D(0, 0, 0);
    auto joint1 = std::make_shared<SSJoint>("Joint_11", 
                                            jointPosRel11, 
                                            MWMath::axisAngle({1,0,0}, 0.0),
                                            joint, 
                                            FJAngles[1], 
                                            MWMath::Point3D(1,0,0), // Abduktion wandert von Y auf X!
                                            numTimeSteps);
    tissues.push_back(joint1);
    

    // ==============================================================================
    // SEGMENT 2 (Distal)
    // ==============================================================================
    // Body 2 startet am Ende des vorherigen Segments
    auto body2 = std::make_shared<SSBody>("Body_Seg2", MWMath::Point3D(0,-L2,0), MWMath::RotMatrix3x3(), joint1);
    tissues.push_back(body2);
    
    auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2, 
            "Mesh_Seg2", body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), COLORFDEACTIVE);
    if (bShowBody) meshes.push_back(mesh2);
    
    // Torus 2
    auto mTorus2 = std::make_shared<SSTorusMesh>(mesh2->B * relTorusR[1], mesh2->B * relTorusr[1], 
            "Torus2", body2, MWMath::Point3D(width[1]*off, -L2 * relTorusPos[1], 0.), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
    meshes.push_back(mTorus2);

    MWMath::Point3D jointPosRel2 = MWMath::Point3D(0, -L2, 0);
    auto joint2 = std::make_shared<SSJoint>("Joint_2", 
                                            jointPosRel2, 
                                            MWMath::RotMatrix3x3(),
                                            body2, 
                                            FJAngles[2], 
                                            MWMath::Point3D(0,0,1), // Bleibt Z
                                            numTimeSteps);
    tissues.push_back(joint2);
    
    double jointSize2 = width[1]*1.1/rWF;
    auto jMesh2 = std::make_shared<SSEllipsoidMesh>(jointSize2, jointSize2, jointSize2, "Mesh_Joint2", joint2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMesh2);

    // ==============================================================================
    // SEGMENT 3
    // ==============================================================================
    auto body3 = std::make_shared<SSBody>("Body_Seg3", MWMath::Point3D(0,-L3,0), MWMath::RotMatrix3x3(), joint2);
    tissues.push_back(body3);
    
    auto mesh3 = std::make_shared<SSEllipsoidMesh>(width[2], width[2], L3, 
            "Mesh_Seg3", body3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), 
            COLORFDEACTIVE);
    if (bShowBody) meshes.push_back(mesh3);
    
    auto mTorus3 = std::make_shared<SSTorusMesh>(width[2] * relTorusR[2], width[2] * relTorusr[2], 
            "Torus3", body3, MWMath::Point3D(width[2]*off, -L3 * relTorusPos[2], 0.), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
    meshes.push_back(mTorus3);

    MWMath::Point3D jointPosRel3 = MWMath::Point3D(0, -L3, 0);
    auto joint3 = std::make_shared<SSJoint>("Joint_3", 
                                            jointPosRel3, 
                                            MWMath::RotMatrix3x3(),
                                            body3, 
                                            FJAngles[3], 
                                            MWMath::Point3D(0,0,1), // Bleibt Z
                                            numTimeSteps);
    tissues.push_back(joint3);

    double jointSize3 = width[2]*1.1/rWF;
    auto jMesh3 = std::make_shared<SSEllipsoidMesh>(jointSize3, jointSize3, jointSize3, "Mesh_Joint3", joint3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMesh3);

    // ==============================================================================
    // SEGMENT 4
    // ==============================================================================
    auto body4 = std::make_shared<SSBody>("Body_Seg4", MWMath::Point3D(0,-L4,0), MWMath::RotMatrix3x3(), joint3);
    tissues.push_back(body4);
    
    auto mesh4 = std::make_shared<SSEllipsoidMesh>(width[3], width[3], L4, 
            "Mesh_Seg4", body4, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), 
            COLORFDEACTIVE);
    if (bShowBody) meshes.push_back(mesh4);

    // ==============================================================================
    // INITIALES UPDATE
    // ==============================================================================
    for (auto& m : meshes) {m->InitializeMesh();}
    rootSystem->update(0); 

    // ==============================================================================
    // MUSKEL
    // ==============================================================================
    int numPoints = cfg.muscleNumPoints[0];
    
    // Muskelpunkte ebenfalls an neue Achsen angepasst: 
    // Offset auf der X-Achse, Länge (C) auf der negativen Y-Achse
    MWMath::Point3D startOffset = MWMath::Point3D(HL*0.03,HL*0.3, 0.0); //MWMath::Point3D(mesh1->B * 1.1, 0.0, 0.0);
    MWMath::Point3D endOffset = MWMath::Point3D(mesh4->B * 1.1, -mesh4->C, 0.0);
    
    SSMuscle* flexor = new SSMuscle("Flexor", numPoints, 
        rootSystem.get(), startOffset, 
        body4.get(), endOffset);

    // Alle Hindernisse
    for(auto& m : meshes) {
        /* if (dynamic_cast<SSTorusMesh*>(m.get()) || m->bIsJointMesh){
            flexor->meshPtrs.push_back(m.get());
        } */
        flexor->meshPtrs.push_back(m.get());
    }

    flexor->createMusclePointsComplexPath();
    //flexor->createMusclePoints();
    flexor->updateMusclePointsParents();
    muscles.push_back(flexor);

    return "buildOHandModelOldExpanded";
}

inline std::string buildOHandModelOldExpanded_paramStudy(std::vector<std::shared_ptr<SSTissue>>& tissues,std::vector<std::shared_ptr<SSMesh>>& meshes,std::vector<SSMuscle*>& muscles,std::shared_ptr<SSBody>& rootSystem,int numTimeSteps,const SimSettings& cfg,float geometryScaler,const std::vector<double> processParams)
{
    // -----------------------------------
    // This function represents the original "OFINGER_SIMPLE_OnlyTorusSmall" Case, where most of the 90°-Rotations 
    // worked, but now as a handmodel builder script.
    // -----------------------------------
    qDebug ()<< "      Running 'buildOHandModelOldExpanded_paramStudy' PARAM STUDY: "<< processParams[0] << ", "<< processParams[1] << ", "<< processParams[2] << ", "<< processParams[3] << ", "<< processParams[4] << ", "<< processParams[5];
    meshes.clear();
    muscles.clear();

    std::vector<double> FWAngles(2); // = {45.0, 0.0};
    std::vector<double> FJAngles(4); // = {90.0, 0.0, 100.0, 90.0}; // MCP1, MCP2, PIP, DIP 
    FWAngles[0] = processParams[0]; // Wrist Abduction
    FWAngles[1] = processParams[1]; // Wrist Flexion
    FJAngles[0] = processParams[2]; // MCP Abd
    FJAngles[1] = processParams[3]; // MCP Flexion
    FJAngles[2] = processParams[4]; // PIP
    FJAngles[3] = processParams[5]; // DIP

    MWMath::Point3D COLORF1 = MWMath::Point3D(0.50, 0.00, 1.00); 
    MWMath::Point3D COLORF2 = MWMath::Point3D(0.30, 0.30, 1.00); 
    MWMath::Point3D COLORF3 = MWMath::Point3D(0.00, 0.50, 1.00); 
    MWMath::Point3D COLORF4 = MWMath::Point3D(0.00, 0.90, 1.00);

    MWMath::Point3D COLORJOINT = MWMath::Point3D(0.0, 0.9, 0.3);
    MWMath::Point3D COLORFDEACTIVE = MWMath::Point3D(0.5, 0.5, 0.5);

    // --- PARAMETER ---
    double HL = 1.8;
    double HB = 0.8;
    double scale = HL / Hand::RefHandLength;
    std::unordered_map<std::string, SSMesh*> MeshMap;
    double Sign = 1.0; // Für Linke Hand -1, Rechte Hand +1
    float GS = 1.0;
    std::vector<double> fLength = {Hand::SegLenRatios[1][0] * HL, Hand::SegLenRatios[1][1] * HL, Hand::SegLenRatios[1][2] * HL, Hand::SegLenRatios[1][3] * HL};
    //std::vector<double> fLength = {0.3482 * HL, 0.2027 * HL, 0.1175 * HL, 0.0882 * HL};
    double L1 = fLength[0]; 
    double L2 = fLength[1]; 
    double L3 = fLength[2]; 
    double L4 = fLength[3]; 
    double rWF = 1.0; 
    std::vector<double> width = {0.15*GS*rWF, 0.1*GS*rWF, 0.07*GS*rWF, 0.05*GS*rWF}; 
    std::vector<double> relTorusPos = {0.6, 0.42, 0.32};
    std::vector<double> relTorusR = {1.3/rWF, 1.3/rWF, 1.6/rWF}; 
    std::vector<double> relTorusr = {1.1/rWF, 1.1/rWF, 1.4/rWF}; 

    bool bShowBody = true;
    double off = 1.0; 

    // 1. ROOT
    MWMath::RotMatrix3x3 Ausrichtung = MWMath::axisAngle({0,0,1}, 90.0) * MWMath::axisAngle({0,1,0}, 180.0);
    rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0, 0, 0), Ausrichtung, nullptr);
    tissues.push_back(rootSystem);

    MWMath::Point3D wristOffset(0.0, 0.003567 * scale, -0.003901 * scale * Sign);
    
    auto wristJointAbd = std::make_shared<SSJoint>("Wrist_JointAbd", wristOffset, MWMath::RotMatrix3x3(), rootSystem, FWAngles[0], MWMath::Point3D(1, 0, 0), numTimeSteps);
    tissues.push_back(wristJointAbd);
    /* auto jMeshWrist = std::make_shared<SSEllipsoidMesh>(0.01 * scale, 0.01 * scale, 0.01 * scale, "Mesh_Wrist_Joint", wristJointAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0,0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMeshWrist);
    MeshMap["Mesh_Wrist_Joint"] = jMeshWrist.get(); */

    auto wristJointFlex = std::make_shared<SSJoint>("Wrist_JointFlex", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointAbd, FWAngles[1], MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointFlex);

    // 2. CARPALS
    auto carpals = std::make_shared<SSBody>("Carpals", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointFlex);
    tissues.push_back(carpals);
    
    double carpalThickness = width[0] * 0.5; 
    double carpalLength = width[0] * 1.25;    
    double carpalWidth = width[0] * 0.5 * 3.5;     
    
    auto meshCarpals = std::make_shared<SSEllipsoidMesh>(carpalThickness, carpalLength, carpalWidth, "Mesh_Carpals", carpals, MWMath::Point3D(0, 0.0*scale, 0), MWMath::RotMatrix3x3(), MWMath::Point3D(0.5, 0.5, 0.5));
    //meshes.push_back(meshCarpals);
    MeshMap["Mesh_Carpals"] = meshCarpals.get();

    // 3. WRAPPING SURFACE (Karpaltunnel)
    double flexCylRadius = 0.0115 * scale;
    double flexCylLength = 0.1400 * scale; 
    MWMath::Point3D flexCylOffset(carpalThickness * 0.5 + flexCylRadius, 0.02 * scale, 0.0);

    auto flexorCylMesh = std::make_shared<SSCylinderMesh>(flexCylRadius, flexCylLength, "Mesh_WristFlexorCyl", rootSystem, flexCylOffset, MWMath::axisAngle({ 1, 0, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));  
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
    // SEGMENT 1 (Proximal)
    // ==============================================================================
    auto body1 = std::make_shared<SSBody>("Body_Seg1", MWMath::Point3D(0,-L1,0), MWMath::RotMatrix3x3(), jointCMCFlex);
    tissues.push_back(body1);
    
    // Mesh 1: Rotiert um X, damit es entlang -Y zeigt
    auto mesh1 = std::make_shared<SSEllipsoidMesh>(width[0], width[0], L1, "Mesh_Seg1", body1, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), COLORFDEACTIVE);
    if (bShowBody) meshes.push_back(mesh1);
    
    // Torus 1: Offset auf +X, Länge auf -Y
    auto mTorus1 = std::make_shared<SSTorusMesh>(mesh1->B * relTorusR[0], mesh1->B * relTorusr[0], 
            "Torus1", body1, MWMath::Point3D(width[0]*off, -L1 * relTorusPos[0], 0.), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(1, 0, 0));
    meshes.push_back(mTorus1);

    // Joint 1: Position am Ende von Segment 1 auf der Y-Achse
    MWMath::Point3D jointPosRel1 = MWMath::Point3D(0, -L1, 0);
    auto joint = std::make_shared<SSJoint>("Joint_1", 
                                            jointPosRel1,           // Position relativ zum Parent
                                            MWMath::RotMatrix3x3(), // Initiale Rotation
                                            body1,                  // Parent Body
                                            FJAngles[0],            // Max Winkel
                                            MWMath::Point3D(0,0,1),// Drehachse: Bleibt Z!
                                            numTimeSteps);
    tissues.push_back(joint);
    
    double jointSize1 = width[0]*1.1/rWF;
    // Gelenkkugel braucht keine spezifische Rotation
    auto jMesh1 = std::make_shared<SSEllipsoidMesh>(jointSize1, jointSize1, jointSize1, "Mesh_Joint1", joint, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMesh1);
    
    MWMath::Point3D jointPosRel11 = MWMath::Point3D(0, 0, 0);
    auto joint1 = std::make_shared<SSJoint>("Joint_11", 
                                            jointPosRel11, 
                                            MWMath::axisAngle({1,0,0}, 0.0),
                                            joint, 
                                            FJAngles[1], 
                                            MWMath::Point3D(1,0,0), // Abduktion wandert von Y auf X!
                                            numTimeSteps);
    tissues.push_back(joint1);
    

    // ==============================================================================
    // SEGMENT 2 (Distal)
    // ==============================================================================
    // Body 2 startet am Ende des vorherigen Segments
    auto body2 = std::make_shared<SSBody>("Body_Seg2", MWMath::Point3D(0,-L2,0), MWMath::RotMatrix3x3(), joint1);
    tissues.push_back(body2);
    
    auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2, 
            "Mesh_Seg2", body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), COLORFDEACTIVE);
    if (bShowBody) meshes.push_back(mesh2);
    
    // Torus 2
    auto mTorus2 = std::make_shared<SSTorusMesh>(mesh2->B * relTorusR[1], mesh2->B * relTorusr[1], 
            "Torus2", body2, MWMath::Point3D(width[1]*off, -L2 * relTorusPos[1], 0.), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
    meshes.push_back(mTorus2);

    MWMath::Point3D jointPosRel2 = MWMath::Point3D(0, -L2, 0);
    auto joint2 = std::make_shared<SSJoint>("Joint_2", 
                                            jointPosRel2, 
                                            MWMath::RotMatrix3x3(),
                                            body2, 
                                            FJAngles[2], 
                                            MWMath::Point3D(0,0,1), // Bleibt Z
                                            numTimeSteps);
    tissues.push_back(joint2);
    
    double jointSize2 = width[1]*1.1/rWF;
    auto jMesh2 = std::make_shared<SSEllipsoidMesh>(jointSize2, jointSize2, jointSize2, "Mesh_Joint2", joint2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMesh2);

    // ==============================================================================
    // SEGMENT 3
    // ==============================================================================
    auto body3 = std::make_shared<SSBody>("Body_Seg3", MWMath::Point3D(0,-L3,0), MWMath::RotMatrix3x3(), joint2);
    tissues.push_back(body3);
    
    auto mesh3 = std::make_shared<SSEllipsoidMesh>(width[2], width[2], L3, 
            "Mesh_Seg3", body3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), 
            COLORFDEACTIVE);
    if (bShowBody) meshes.push_back(mesh3);
    
    auto mTorus3 = std::make_shared<SSTorusMesh>(width[2] * relTorusR[2], width[2] * relTorusr[2], 
            "Torus3", body3, MWMath::Point3D(width[2]*off, -L3 * relTorusPos[2], 0.), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
    meshes.push_back(mTorus3);

    MWMath::Point3D jointPosRel3 = MWMath::Point3D(0, -L3, 0);
    auto joint3 = std::make_shared<SSJoint>("Joint_3", 
                                            jointPosRel3, 
                                            MWMath::RotMatrix3x3(),
                                            body3, 
                                            FJAngles[3], 
                                            MWMath::Point3D(0,0,1), // Bleibt Z
                                            numTimeSteps);
    tissues.push_back(joint3);

    double jointSize3 = width[2]*1.1/rWF;
    auto jMesh3 = std::make_shared<SSEllipsoidMesh>(jointSize3, jointSize3, jointSize3, "Mesh_Joint3", joint3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMesh3);

    // ==============================================================================
    // SEGMENT 4
    // ==============================================================================
    auto body4 = std::make_shared<SSBody>("Body_Seg4", MWMath::Point3D(0,-L4,0), MWMath::RotMatrix3x3(), joint3);
    tissues.push_back(body4);
    
    auto mesh4 = std::make_shared<SSEllipsoidMesh>(width[3], width[3], L4, 
            "Mesh_Seg4", body4, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), 
            COLORFDEACTIVE);
    if (bShowBody) meshes.push_back(mesh4);

    // ==============================================================================
    // INITIALES UPDATE
    // ==============================================================================
    for (auto& m : meshes) {m->InitializeMesh();}
    rootSystem->update(0); 

    // ==============================================================================
    // MUSKEL
    // ==============================================================================
    int numPoints = cfg.muscleNumPoints[0];
    
    // Muskelpunkte ebenfalls an neue Achsen angepasst: 
    // Offset auf der X-Achse, Länge (C) auf der negativen Y-Achse
    MWMath::Point3D startOffset = MWMath::Point3D(HL*0.03,HL*0.3, 0.0); //MWMath::Point3D(mesh1->B * 1.1, 0.0, 0.0);
    MWMath::Point3D endOffset = MWMath::Point3D(mesh4->B * 1.1, -mesh4->C, 0.0);
    
    SSMuscle* flexor = new SSMuscle("Flexor", numPoints, 
        rootSystem.get(), startOffset, 
        body4.get(), endOffset);

    // Alle Hindernisse
    for(auto& m : meshes) {
        /* if (dynamic_cast<SSTorusMesh*>(m.get()) || m->bIsJointMesh){
            flexor->meshPtrs.push_back(m.get());
        } */
        flexor->meshPtrs.push_back(m.get());
    }

    flexor->createMusclePointsComplexPath();
    //flexor->createMusclePoints();
    flexor->updateMusclePointsParents();
    muscles.push_back(flexor);

    return "buildOHandModelOldExpanded_paramStudy";
}



inline std::string buildOHandModelOldExpandedLoop(std::vector<std::shared_ptr<SSTissue>>& tissues,std::vector<std::shared_ptr<SSMesh>>& meshes,std::vector<SSMuscle*>& muscles,std::shared_ptr<SSBody>& rootSystem,int numTimeSteps,const SimSettings& cfg,float geometryScaler,const std::vector<double> processParams)
{
    // -----------------------------------
    // This function represents the original "OFINGER_SIMPLE_OnlyTorusSmall" Case, where most of the 90°-Rotations 
    // worked, but now as a handmodel builder script.
    // -----------------------------------

    meshes.clear();
    muscles.clear();
    std::unordered_map<std::string, SSMesh*> MeshMap;
    std::vector<double> fLength;
    //std::vector<double> fLength = {0.3482 * HL, 0.2027 * HL, 0.1175 * HL, 0.0882 * HL};

    // --- STATIC ---
    std::vector<double> HAngles = processParams.empty() ? std::vector<double>{0.0,  0.0, 0.0,  0.0, 0.0, 0.0} : processParams;
    HAngles = { 0.0,  0.0, 90.0,  0.0, 100.0, 80.0};
    int numPoints = cfg.muscleNumPoints[0];
    if (processParams.size() > 6) numPoints = static_cast<int>(processParams[6]);
    std::vector<double> FWAngles = {HAngles[0], HAngles[1]}; // Wrist Abduction, Wrist Flexion
    std::vector<double> FJAngles = {HAngles[2], HAngles[3], HAngles[4], HAngles[5]}; // MCP1, MCP2, PIP, DIP 
    const MWMath::Point3D COLORF1 = MWMath::Point3D(0.50, 0.00, 1.00); 
    const MWMath::Point3D COLORF2 = MWMath::Point3D(0.30, 0.30, 1.00); 
    const MWMath::Point3D COLORF3 = MWMath::Point3D(0.00, 0.50, 1.00); 
    const MWMath::Point3D COLORF4 = MWMath::Point3D(0.00, 0.90, 1.00);
    const MWMath::Point3D COLORFLEXOR = MWMath::Point3D(03, 0.03, 0.03);
    const MWMath::Point3D COLORJOINT = MWMath::Point3D(0.0, 0.9, 0.3);
    const MWMath::Point3D COLORFDEACTIVE = MWMath::Point3D(0.5, 0.5, 0.5);
    const bool bShowBody = true;
    qDebug ()<< "      Running 'buildOHandModelOldExpandedLoop' with HAngles: "<< HAngles[0] << ", "<< HAngles[1] << ", "<< HAngles[2] << ", "<< HAngles[3] << ", "<< HAngles[4] << ", "<< HAngles[5] << " and numPoints: " << numPoints;
    

    // --- PARAMETER ---
    const double HL = 1.8;
    const double HB = 0.8;
    double scale = HL / Hand::RefHandLength;
    double Sign = 1.0; // Für Linke Hand -1, Rechte Hand +1
    float GS = 1.0;
    double rWF = 1.0;
    std::vector<double> width = {0.15*GS*rWF, 0.1*GS*rWF, 0.07*GS*rWF, 0.05*GS*rWF}; 
    
    // --- TORI ---
    std::vector<double> relTorusPos = {0.6, 0.42, 0.32};
    std::vector<double> relTorusR = {1.3/rWF, 1.3/rWF, 1.6/rWF}; 
    std::vector<double> relTorusr = {1.1/rWF, 1.1/rWF, 1.4/rWF}; 
    double off = 1.0; // vertical torus offset relative to width
    

    // --- ROOT ---
    MWMath::RotMatrix3x3 Ausrichtung = MWMath::axisAngle({0,0,1}, 90.0) * MWMath::axisAngle({0,1,0}, 180.0);
    rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0, 0, 0), Ausrichtung, nullptr);
    tissues.push_back(rootSystem);

    // --- WRIST
    MWMath::Point3D wristOffset(0.0, 0.003567 * scale, -0.003901 * scale * Sign);
    auto wristJointAbd = std::make_shared<SSJoint>("Wrist_JointAbd", wristOffset, MWMath::RotMatrix3x3(), rootSystem, FWAngles[0], MWMath::Point3D(1, 0, 0), numTimeSteps);
    tissues.push_back(wristJointAbd);
    /* auto jMeshWrist = std::make_shared<SSEllipsoidMesh>(0.01 * scale, 0.01 * scale, 0.01 * scale, "Mesh_Wrist_Joint", wristJointAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0,0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
    meshes.push_back(jMeshWrist);
    MeshMap["Mesh_Wrist_Joint"] = jMeshWrist.get(); */
    auto wristJointFlex = std::make_shared<SSJoint>("Wrist_JointFlex", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointAbd, FWAngles[1], MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointFlex);

    // --- CARPALS
    auto carpals = std::make_shared<SSBody>("Carpals", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointFlex);
    tissues.push_back(carpals);
    double carpalThickness = width[0] * 0.5; 
    double carpalLength = width[0] * 1.25;    
    double carpalWidth = width[0] * 0.5 * 3.5;     
    auto meshCarpals = std::make_shared<SSEllipsoidMesh>(carpalThickness, carpalLength, carpalWidth, "Mesh_Carpals", carpals, MWMath::Point3D(0, 0.0*scale, 0), MWMath::RotMatrix3x3(), MWMath::Point3D(0.5, 0.5, 0.5));
    //meshes.push_back(meshCarpals);
    MeshMap["Mesh_Carpals"] = meshCarpals.get();

    // --- WRAPPING SURFACE (Karpaltunnel)
    double flexCylRadius = 0.0115 * scale;
    double flexCylLength = 0.1400 * scale; 
    MWMath::Point3D flexCylOffset(carpalThickness * 0.5 + flexCylRadius, 0.02 * scale, 0.0);
    auto flexorCylMesh = std::make_shared<SSCylinderMesh>(flexCylRadius, flexCylLength, "Mesh_WristFlexorCyl", rootSystem, flexCylOffset, MWMath::axisAngle({ 1, 0, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));  
    meshes.push_back(flexorCylMesh);
    MeshMap["Mesh_WristFlexorCyl"] = flexorCylMesh.get();


    // ==============================================================================
    // FINGER AUFBAU (Index Finger = fidx 1)
    // ==============================================================================
    for (int i = 1; i < 2; i++){
        // ------------------------------------
        // INIT AND PARAMETERS
        // ------------------------------------
        int fidx = i;
        std::string cFName = "Index";
        std::string prefN = std::to_string(fidx);
        double jAngle;
        std::string meshName;
        fLength = {Hand::SegLenRatios[fidx][0] * HL, Hand::SegLenRatios[fidx][1] * HL, Hand::SegLenRatios[fidx][2] * HL, Hand::SegLenRatios[fidx][3] * HL};
        double L1 = fLength[0]; 
        double L2 = fLength[1]; 
        double L3 = fLength[2]; 
        double L4 = fLength[3]; 


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
        meshName = "Mesh_CMC" + prefN + "_Joint";
        auto jMeshCMC = std::make_shared<SSEllipsoidMesh>(width[0], width[0], width[0], meshName, jointCMCAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0 , 0 }, 90.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMeshCMC);
        MeshMap[meshName] = jMeshCMC.get();

        // CMC Joint Flex
        jAngle = 0.0;
        jName = "CMC" + prefN + "_JointFlex";
        auto jointCMCFlex = std::make_shared<SSJoint>(jName, MWMath::Point3D(0.,0.,0.), MWMath::RotMatrix3x3(), jointCMCAbd, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointCMCFlex);



        // ==============================================================================
        // SEGMENT 1 (Proximal)
        // ==============================================================================
        auto body1 = std::make_shared<SSBody>("MC"+prefN, MWMath::Point3D(0,-L1,0), MWMath::RotMatrix3x3(), jointCMCFlex);
        tissues.push_back(body1);
        
        // Mesh 1: Rotiert um X, damit es entlang -Y zeigt
        meshName = "Mesh_MC" + prefN;
        auto mesh1 = std::make_shared<SSEllipsoidMesh>(width[0], width[0], L1, meshName, body1, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), COLORF1);
        if (bShowBody) meshes.push_back(mesh1);
        MeshMap[meshName] = mesh1.get();
        
        // Torus 1: Offset auf +X, Länge auf -Y
        meshName = "MeshTorus_MC" + prefN;
        auto mTorus1 = std::make_shared<SSTorusMesh>(mesh1->B * relTorusR[0], mesh1->B * relTorusr[0], 
                "Torus1", body1, MWMath::Point3D(width[0]*off, -L1 * relTorusPos[0], 0.), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(1, 0, 0));
        meshes.push_back(mTorus1);
        MeshMap[meshName] = mTorus1.get();

        // Joint 1: Position am Ende von Segment 1 auf der Y-Achse
        MWMath::Point3D jointPosRel1 = MWMath::Point3D(0, -L1, 0);
        auto joint = std::make_shared<SSJoint>("Joint_1", 
                                                jointPosRel1,           // Position relativ zum Parent
                                                MWMath::RotMatrix3x3(), // Initiale Rotation
                                                body1,                  // Parent Body
                                                FJAngles[0],            // Max Winkel
                                                MWMath::Point3D(0,0,1),// Drehachse: Bleibt Z!
                                                numTimeSteps);
        tissues.push_back(joint);
        
        double jointSize1 = width[0]*1.1/rWF;
        // Gelenkkugel braucht keine spezifische Rotation
        meshName = "Mesh_MC" + prefN + "_Joint";
        auto jMesh1 = std::make_shared<SSEllipsoidMesh>(jointSize1, jointSize1, jointSize1, meshName, joint, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh1);
        MeshMap[meshName] = jMesh1.get();
        
        MWMath::Point3D jointPosRel11 = MWMath::Point3D(0, 0, 0);
        auto joint1 = std::make_shared<SSJoint>("Joint_11", 
                                                jointPosRel11, 
                                                MWMath::axisAngle({1,0,0}, 0.0),
                                                joint, 
                                                FJAngles[1], 
                                                MWMath::Point3D(1,0,0), // Abduktion wandert von Y auf X!
                                                numTimeSteps);
        tissues.push_back(joint1);
        

        // ==============================================================================
        // SEGMENT 2 (Distal)
        // ==============================================================================
        // Body 2 startet am Ende des vorherigen Segments
        auto body2 = std::make_shared<SSBody>("PP" + prefN, MWMath::Point3D(0,-L2,0), MWMath::RotMatrix3x3(), joint1);
        tissues.push_back(body2);
        
        meshName = "Mesh_PP" + prefN;
        auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2, 
                meshName, body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), COLORF2);
        if (bShowBody) meshes.push_back(mesh2);
        MeshMap[meshName] = mesh2.get();
        
        // Torus 2
        meshName = "MeshTorus_PP" + prefN;
        auto mTorus2 = std::make_shared<SSTorusMesh>(mesh2->B * relTorusR[1], mesh2->B * relTorusr[1], 
                meshName, body2, MWMath::Point3D(width[1]*off, -L2 * relTorusPos[1], 0.), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0.8, 0.8, 0));
        meshes.push_back(mTorus2);
        MeshMap[meshName] = mTorus2.get();

        MWMath::Point3D jointPosRel2 = MWMath::Point3D(0, -L2, 0);
        auto joint2 = std::make_shared<SSJoint>("Joint_2", 
                                                jointPosRel2, 
                                                MWMath::RotMatrix3x3(),
                                                body2, 
                                                FJAngles[2], 
                                                MWMath::Point3D(0,0,1), // Bleibt Z
                                                numTimeSteps);
        tissues.push_back(joint2);
        
        double jointSize2 = width[1]*1.1/rWF;
        meshName = "Mesh_PIP" + prefN + "_Joint";
        auto jMesh2 = std::make_shared<SSEllipsoidMesh>(jointSize2, jointSize2, jointSize2, meshName, joint2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh2);
        MeshMap[meshName] = jMesh2.get();

        // ==============================================================================
        // SEGMENT 3
        // ==============================================================================
        auto body3 = std::make_shared<SSBody>("MP" + prefN, MWMath::Point3D(0,-L3,0), MWMath::RotMatrix3x3(), joint2);
        tissues.push_back(body3);
        
        meshName = "Mesh_MP" + prefN;
        auto mesh3 = std::make_shared<SSEllipsoidMesh>(width[2], width[2], L3, 
                meshName, body3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), 
                COLORF3);
        if (bShowBody) meshes.push_back(mesh3);
        MeshMap[meshName] = mesh3.get();
        
        meshName = "MeshTorus_MP" + prefN;
        auto mTorus3 = std::make_shared<SSTorusMesh>(width[2] * relTorusR[2], width[2] * relTorusr[2], 
                meshName, body3, MWMath::Point3D(width[2]*off, -L3 * relTorusPos[2], 0.), MWMath::axisAngle({1,0,0}, 90.0), MWMath::Point3D(0, 0, 0.8));
        meshes.push_back(mTorus3);
        MeshMap[meshName] = mTorus3.get();

        MWMath::Point3D jointPosRel3 = MWMath::Point3D(0, -L3, 0);
        auto joint3 = std::make_shared<SSJoint>("Joint_3", 
                                                jointPosRel3, 
                                                MWMath::RotMatrix3x3(),
                                                body3, 
                                                FJAngles[3], 
                                                MWMath::Point3D(0,0,1), // Bleibt Z
                                                numTimeSteps);
        tissues.push_back(joint3);

        double jointSize3 = width[2]*1.1/rWF;
        meshName = "Mesh_PIP" + prefN + "_Joint";
        auto jMesh3 = std::make_shared<SSEllipsoidMesh>(jointSize3, jointSize3, jointSize3, meshName, joint3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));
        meshes.push_back(jMesh3);
        MeshMap[meshName] = jMesh3.get();

        // ==============================================================================
        // SEGMENT 4
        // ==============================================================================
        auto body4 = std::make_shared<SSBody>("DP" + prefN, MWMath::Point3D(0,-L4,0), MWMath::RotMatrix3x3(), joint3);
        tissues.push_back(body4);
        
        meshName = "Mesh_DP" + prefN;
        auto mesh4 = std::make_shared<SSEllipsoidMesh>(width[3], width[3], L4, 
                meshName, body4, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), 
                COLORF4);
        if (bShowBody) meshes.push_back(mesh4);
        MeshMap[meshName] = mesh4.get();
    }

    // ==============================================================================
    // INITIALES UPDATE
    // ==============================================================================
    for (auto& m : meshes) {m->InitializeMesh();}
    rootSystem->update(0); 

    // ==============================================================================
    // MUSKEL
    // ==============================================================================
    
    // Muskelpunkte ebenfalls an neue Achsen angepasst: 
    // Offset auf der X-Achse, Länge (C) auf der negativen Y-Achse
    MWMath::Point3D startOffset = MWMath::Point3D(HL*0.03,HL*0.3, 0.0); //MWMath::Point3D(mesh1->B * 1.1, 0.0, 0.0);
    auto originMesh = dynamic_cast<SSEllipsoidMesh*>(MeshMap["Mesh_MC1"]);
    auto insertionMesh = dynamic_cast<SSEllipsoidMesh*>(MeshMap["Mesh_DP1"]);
    // MWMath::Point3D endOffset = MWMath::Point3D(insertionMesh->B * 1.1, 0.0, 0.0);
    MWMath::Point3D endOffset = MWMath::Point3D(insertionMesh->B * 1.1, -insertionMesh->C, 0.0);
    
    SSMuscle* flexor = new SSMuscle("Flexor", numPoints, 
        rootSystem.get(), startOffset, 
        getRefTissue("DP1", tissues), endOffset);

    // Alle Hindernisse
    for(auto& m : meshes) {
        /* if (dynamic_cast<SSTorusMesh*>(m.get()) || m->bIsJointMesh){
            flexor->meshPtrs.push_back(m.get());
        } */
        flexor->meshPtrs.push_back(m.get());
    }

    flexor->createMusclePointsComplexPath();
    //flexor->createMusclePoints();
    flexor->updateMusclePointsParents();
    muscles.push_back(flexor);

    return "buildOHandModelOldExpandedLoop";
}



inline std::string buildOHandModelOldExpandedViaX05(std::vector<std::shared_ptr<SSTissue>>& tissues,std::vector<std::shared_ptr<SSMesh>>& meshes,std::vector<SSMuscle*>& muscles,std::shared_ptr<SSBody>& rootSystem,int numTimeSteps,const SimSettings& cfg,float geometryScaler,const std::vector<double> processParams) {
    // -----------------------------------
    // This function represents the original "OFINGER_SIMPLE_OnlyTorusSmall" Case, where most of the 90°-Rotations 
    // worked, but now as a handmodel builder script.
    // -----------------------------------

    meshes.clear();
    muscles.clear();
    std::unordered_map<std::string, SSMesh*> MeshMap;
    std::vector<double> fLength;
    //std::vector<double> fLength = {0.3482 * HL, 0.2027 * HL, 0.1175 * HL, 0.0882 * HL};

    // --- STATIC ---
    std::vector<double> HAngles = processParams.empty() ? std::vector<double>{0.0,  0.0, 0.0,  0.0, 0.0, 0.0} : processParams;
    // HAngles = { 0.0,  0.0, 90.0,  0.0, 100.0, 80.0}; qDebug() << "      Using default HAngles for 'buildOHandModelOldExpandedVia', since processParams is empty.";
    int numPoints = cfg.muscleNumPoints[0];
    if (processParams.size() > 6) numPoints = static_cast<int>(processParams[6]);
    std::vector<double> FWAngles = {HAngles[0], HAngles[1]}; // Wrist Abduction, Wrist Flexion
    std::vector<double> FJAngles = {HAngles[2], HAngles[3], HAngles[4], HAngles[5]}; // MCP1, MCP2, PIP, DIP 
    double BCol = 0.8;
    const MWMath::Point3D COLORF1 = MWMath::Point3D(BCol, BCol, BCol); // MWMath::Point3D(0.50, 0.00, 1.00); 
    const MWMath::Point3D COLORF2 = MWMath::Point3D(BCol, BCol, BCol); // MWMath::Point3D(0.30, 0.30, 1.00); 
    const MWMath::Point3D COLORF3 = MWMath::Point3D(BCol, BCol, BCol); // MWMath::Point3D(0.00, 0.50, 1.00); 
    const MWMath::Point3D COLORF4 = MWMath::Point3D(BCol, BCol, BCol); // MWMath::Point3D(0.00, 0.90, 1.00);
    const MWMath::Point3D COLORVPFLEXOR = MWMath::Point3D(0.6, 0.1, 0.1);
    const MWMath::Point3D COLORVPEXTENSOR = MWMath::Point3D(0.1, 0.1, 0.6);
    const MWMath::Point3D COLORJOINT = MWMath::Point3D(0.00, 0.90, 1.00); // MWMath::Point3D(0.0, 0.9, 0.3);
    const MWMath::Point3D COLORFDEACTIVE = MWMath::Point3D(0.5, 0.5, 0.5);
    const bool bShowBody = true;
    qDebug ()<< "      Running 'buildOHandModelOldExpandedViaX' with HAngles: "<< HAngles[0] << ", "<< HAngles[1] << ", "<< HAngles[2] << ", "<< HAngles[3] << ", "<< HAngles[4] << ", "<< HAngles[5] << " and numPoints: " << numPoints;
    

    // --- PARAMETER ---
    const double HL = 1.8;
    const double HB = 0.8;
    double scale = HL / Hand::RefHandLength;
    double Sign = 1.0; // Für Linke Hand -1, Rechte Hand +1
    float GS = 1.0;
    double rWF = 0.5;
    std::vector<double> width = {0.15*GS*rWF, 0.1*GS*rWF, 0.07*GS*rWF, 0.05*GS*rWF}; 
    double viaPointR = 0.05;
    
    // --- TORI ---
    std::vector<double> relTorusPos = {0.7, 0.55, 0.45}; //{0.6, 0.42, 0.32};
    std::vector<double> relTorusR = {1.3/rWF, 1.3/rWF, 1.6/rWF}; 
    std::vector<double> relTorusr = {1.1/rWF, 1.1/rWF, 1.4/rWF}; 
    double off = 1.0; // vertical torus offset relative to width
    

    // --- ROOT ---
    MWMath::RotMatrix3x3 Ausrichtung = MWMath::axisAngle({0,0,1}, 90.0) * MWMath::axisAngle({0,1,0}, 180.0);
    rootSystem = std::make_shared<SSBody>("Root", MWMath::Point3D(0, 0, 0), Ausrichtung, nullptr);
    tissues.push_back(rootSystem);

    // --- WRIST
    MWMath::Point3D wristOffset(0.0, 0.003567 * scale, -0.003901 * scale * Sign);
    auto wristJointAbd = std::make_shared<SSJoint>("Wrist_JointAbd", wristOffset, MWMath::RotMatrix3x3(), rootSystem, FWAngles[0], MWMath::Point3D(1, 0, 0), numTimeSteps);
    tissues.push_back(wristJointAbd);
    auto jMeshWrist = std::make_shared<SSEllipsoidMesh>(0.01 * scale, 0.01 * scale, 0.01 * scale, "Mesh_Wrist_Joint", wristJointAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0,0 }, 0.0), COLORJOINT);
    meshes.push_back(jMeshWrist);
    MeshMap["Mesh_Wrist_Joint"] = jMeshWrist.get();
    auto wristJointFlex = std::make_shared<SSJoint>("Wrist_JointFlex", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointAbd, FWAngles[1], MWMath::Point3D(0, 0, 1), numTimeSteps);
    tissues.push_back(wristJointFlex);

    // --- CARPALS
    auto carpals = std::make_shared<SSBody>("Carpals", MWMath::Point3D(0, 0, 0), MWMath::RotMatrix3x3(), wristJointFlex);
    tissues.push_back(carpals);
    double carpalThickness = 0.01*scale; // width[0] * 0.5; 
    double carpalLength = 0.174*HL*0.5; // relativ zur hand nach prometheus abbildung //0.02*scale; //width[0] * 1.25;    
    double carpalWidth = 1.6*carpalLength; // relativ zur hand nach prometheus abbildung //0.045*scale; // width[0] * 0.5 * 3.5;     
    auto meshCarpals = std::make_shared<SSEllipsoidMesh>(carpalThickness*1.7*0.7, carpalLength, carpalWidth, "Mesh_Carpals", carpals, MWMath::Point3D(-0.017, 0.0*scale, 0), MWMath::RotMatrix3x3(), COLORF1);
    //auto meshCarpals = std::make_shared<SSCylinderMesh>(carpalThickness*1.7, carpalWidth, "Mesh_Carpals", carpals, MWMath::Point3D(0, 0.0*scale, 0), MWMath::RotMatrix3x3(), COLORF1);
    meshes.push_back(meshCarpals);
    MeshMap["Mesh_Carpals"] = meshCarpals.get();


    // HELP MESHES...DEBUGGING PURPOSES
    double armScale = 0.01; //mm -> dm; 
    // --- HUMERUS ---
    auto bodyHumerus = std::make_shared<SSBody>("Body_Humerus", Hand::WristToHumerusMM * armScale, Hand::RotCarpalsToHumerus, rootSystem);
    tissues.push_back(bodyHumerus);
    // --- RADIUS ---
    auto bodyRadius = std::make_shared<SSBody>("Body_Radius", Hand::WristToRadiusMM * armScale, Hand::RotCarpalsToRadius, rootSystem);
    tissues.push_back(bodyRadius);
    // --- ULNA ---
    auto bodyUlna = std::make_shared<SSBody>("Body_Ulna", Hand::WristToUlnaMM * armScale, Hand::RotCarpalsToUlna, rootSystem);
    tissues.push_back(bodyUlna);
    MWMath::Point3D ulnaOffset = MWMath::Point3D(-Hand::WristToUlnaMM.x*armScale, Hand::WristToUlnaMM.y*armScale*0.5, -0.15);
    auto meshUlna = std::make_shared<SSEllipsoidMesh>(carpalThickness, Hand::WristToUlnaMM.y*armScale*0.4, carpalThickness, "Mesh_UlnaPseudo", bodyUlna, ulnaOffset, MWMath::RotMatrix3x3(), COLORJOINT);
    meshes.push_back(meshUlna);
    MeshMap["Mesh_UlnaPseudo"] = meshUlna.get();

    // --- WRAPPING SURFACE (Karpaltunnel)
    // "Mesh_BigCarpalTunnel", "Mesh_MC1"
    /* double flexCylRadius = 0.0115 * scale;
    double flexCylLength = 0.1400 * scale; 
    //MWMath::Point3D flexCylOffset(carpalThickness * 0.5 + flexCylRadius, 0.02 * scale, 0.0);
    MWMath::Point3D flexCylOffset(carpalThickness * 0.5 + flexCylRadius, 0.003567 * scale, 0.0);
    auto flexorCylMesh = std::make_shared<SSCylinderMesh>(flexCylRadius, flexCylLength, "Mesh_WristFlexorCyl", rootSystem, flexCylOffset, MWMath::axisAngle({ 1, 0, 0 }, 0.0), MWMath::Point3D(0.8, 0.8, 0.0));  
    meshes.push_back(flexorCylMesh);
    MeshMap["Mesh_WristFlexorCyl"] = flexorCylMesh.get(); */
    MWMath::RotMatrix3x3 torusRot = MWMath::axisAngle({1, 0, 0}, 90.0); 
    double holeRadius = 1.3*0.0025 * scale;   // Dicke des Knochen/Band-Rings (1 cm)
    double tubethickness = 1.3*0.01 * scale;
    // RETINACULUM FLEXORUM TORUS
    /* std::string meshRetFlexTorus = "MeshTorus_RetinaculumFlexorum"; 
    auto meshRetFlex = std::make_shared<SSTorusMesh>(tubethickness+holeRadius, tubethickness-holeRadius,  meshRetFlexTorus,  carpals,  
        MWMath::Point3D(carpalThickness+holeRadius*2, -carpalLength*0.5, 0.0), torusRot, COLORVPFLEXOR);
    meshes.push_back(meshRetFlex);
    MeshMap[meshRetFlexTorus] = meshRetFlex.get(); */
    // RETINACULUM EXTENSORUM TORUS
    /* std::string meshRetExtTorus = "MeshTorus_RetinaculumExtensorum"; 
    auto meshRetExt = std::make_shared<SSTorusMesh>(tubethickness+holeRadius, tubethickness-holeRadius,  meshRetExtTorus,  carpals,  
        MWMath::Point3D(-(carpalThickness+holeRadius*2), carpalLength*0.4, 0.0), torusRot,  COLORVPEXTENSOR);
    meshRetExt->B = 2.0;
    meshes.push_back(meshRetExt);
    MeshMap[meshRetExtTorus] = meshRetExt.get(); */
    // BIG CARPAL TUNNEL TORUS
    /* double bigHoleRadius = 0.0085 * scale;   // Dicke des Knochen/Band-Rings (1 cm)
    double bigTubethickness = 0.02 * scale;  // Abstand zur Schlauchmitte (1,5 cm)
    auto bigMeshCarpalTunnel = std::make_shared<SSTorusMesh>(bigTubethickness+bigHoleRadius, bigTubethickness-bigHoleRadius,  "Mesh_BigCarpalTunnel",  carpals,  
        MWMath::Point3D(0.0, -carpalLength*0.2, 0.0), torusRot,  MWMath::Point3D(0.6, 0.6, 0.5));
    meshes.push_back(bigMeshCarpalTunnel);
    MeshMap["Mesh_BigCarpalTunnel"] = bigMeshCarpalTunnel.get(); */

    buildWristViaPoints(meshes, MeshMap, carpals, scale=scale, viaPointR =viaPointR);


    // ==============================================================================
    // THUMB AUFBAU (Index 0)
    // ==============================================================================   
    buildThumbModel(tissues, meshes, MeshMap, carpals, numTimeSteps, FJAngles, scale, Sign, HL, HB, MWMath::Point3D(1.0, 0.5, 0.0));
    
    // ==============================================================================
    // FINGER AUFBAU (Index Finger = fidx 1)
    // ==============================================================================
    for (int i = 1; i < 5; i++){
        // ------------------------------------
        // INIT AND PARAMETERS
        // ------------------------------------
        int fidx = i;
        if (fidx > 4) {printColorMsg("Ddnt create more fingers than 5", 2); continue;}
        std::string cFName = "Index";
        std::string prefN = std::to_string(fidx);
        double jAngle;
        std::string meshName;
        fLength = {Hand::SegLenRatios[fidx][0] * HL, Hand::SegLenRatios[fidx][1] * HL, Hand::SegLenRatios[fidx][2] * HL, Hand::SegLenRatios[fidx][3] * HL};
        double L1 = fLength[0]; 
        double L2 = fLength[1]; 
        double L3 = fLength[2]; 
        double L4 = fLength[3]; 
        qDebug() << "      Building finger '" << cFName.c_str() << "' with segment lengths: " << L1 << ", " << L2 << ", " << L3 << ", " << L4;


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
        meshName = "Mesh_CMC" + prefN + "_Joint";
        auto jMeshCMC = std::make_shared<SSEllipsoidMesh>(width[0], width[0], width[0], meshName, jointCMCAbd, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({ 1,0 , 0 }, 90.0), COLORJOINT);
        meshes.push_back(jMeshCMC);
        MeshMap[meshName] = jMeshCMC.get();

        // CMC Joint Flex
        jAngle = 0.0;
        jName = "CMC" + prefN + "_JointFlex";
        auto jointCMCFlex = std::make_shared<SSJoint>(jName, MWMath::Point3D(0.,0.,0.), MWMath::RotMatrix3x3(), jointCMCAbd, jAngle, MWMath::Point3D(0, 0, 1), numTimeSteps);
        tissues.push_back(jointCMCFlex);



        // ==============================================================================
        // SEGMENT 1 (Proximal)
        // ==============================================================================
        auto body1 = std::make_shared<SSBody>("MC"+prefN, MWMath::Point3D(0,-L1*0.5,0), MWMath::RotMatrix3x3(), jointCMCFlex);
        tissues.push_back(body1);
        
        // Mesh 1: Rotiert um X, damit es entlang -Y zeigt
        meshName = "Mesh_MC" + prefN;
        auto mesh1 = std::make_shared<SSEllipsoidMesh>(width[0], width[0], L1*0.5, meshName, body1, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), COLORF1);
        if (bShowBody) meshes.push_back(mesh1);
        MeshMap[meshName] = mesh1.get();
        
        // VIA POINT MC
        meshName = "MeshVPF_MC" + prefN;
        auto mVPF1 = std::make_shared<SSEllipsoidMesh>(viaPointR, viaPointR, viaPointR, 
                meshName, body1, MWMath::Point3D(width[0]*off, -L1*0.5 * relTorusPos[0], 0.), MWMath::axisAngle({1,0,0}, 90.0), COLORVPFLEXOR);
        mVPF1->bIsViaPoint = true;
        meshes.push_back(mVPF1);
        MeshMap[meshName] = mVPF1.get();

        /* meshName = "MeshVPE_MC" + prefN;
        auto mVPE1 = std::make_shared<SSEllipsoidMesh>(viaPointR, viaPointR, viaPointR, 
                meshName, body1, MWMath::Point3D(-width[0]*off, -L1*0.5 * relTorusPos[0], 0.), MWMath::axisAngle({1,0,0}, 90.0), COLORVPEXTENSOR);
        mVPE1->bIsViaPoint = true;
        meshes.push_back(mVPE1);
        MeshMap[meshName] = mVPE1.get(); */

        // Joint 1: Position am Ende von Segment 1 auf der Y-Achse
        MWMath::Point3D jointPosRel1 = MWMath::Point3D(0, -L1*0.5, 0);
        auto joint = std::make_shared<SSJoint>("Joint_1", 
                                                jointPosRel1,           // Position relativ zum Parent
                                                MWMath::RotMatrix3x3(), // Initiale Rotation
                                                body1,                  // Parent Body
                                                FJAngles[0],            // Max Winkel
                                                MWMath::Point3D(0,0,1),// Drehachse: Bleibt Z!
                                                numTimeSteps);
        tissues.push_back(joint);
        
        double jointSize1 = width[0]*1.1;
        // Gelenkkugel braucht keine spezifische Rotation
        meshName = "Mesh_MCP" + prefN + "_Joint";
        auto jMesh1 = std::make_shared<SSEllipsoidMesh>(jointSize1, jointSize1, jointSize1, meshName, joint, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 0.0), COLORJOINT);
        meshes.push_back(jMesh1);
        MeshMap[meshName] = jMesh1.get();

        meshName = "MeshVPE_MC" + prefN;
        auto mVPE1 = std::make_shared<SSEllipsoidMesh>(viaPointR, viaPointR, viaPointR, 
                meshName, joint, MWMath::Point3D(-jointSize1*1.2, 0., 0.), MWMath::axisAngle({1,0,0}, 90.0), COLORVPEXTENSOR);
        mVPE1->bIsViaPoint = true;
        meshes.push_back(mVPE1);
        MeshMap[meshName] = mVPE1.get();
        
        MWMath::Point3D jointPosRel11 = MWMath::Point3D(0, 0, 0);
        auto joint1 = std::make_shared<SSJoint>("Joint_11", 
                                                jointPosRel11, 
                                                MWMath::axisAngle({1,0,0}, 0.0),
                                                joint, 
                                                FJAngles[1], 
                                                MWMath::Point3D(1,0,0), // Abduktion wandert von Y auf X!
                                                numTimeSteps);
        tissues.push_back(joint1);
        

        // ==============================================================================
        // SEGMENT 2 (Distal)
        // ==============================================================================
        // Body 2 startet am Ende des vorherigen Segments
        auto body2 = std::make_shared<SSBody>("PP" + prefN, MWMath::Point3D(0,-L2*0.5,0), MWMath::RotMatrix3x3(), joint1);
        tissues.push_back(body2);
        
        meshName = "Mesh_PP" + prefN;
        auto mesh2 = std::make_shared<SSEllipsoidMesh>(width[1], width[1], L2*0.5, 
                meshName, body2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), COLORF2);
        if (bShowBody) meshes.push_back(mesh2);
        MeshMap[meshName] = mesh2.get();
        
        // VIA POINTS PP
        meshName = "MeshVPF_PP" + prefN;
        auto mVPF2 = std::make_shared<SSEllipsoidMesh>(viaPointR, viaPointR, viaPointR, 
                meshName, body2, MWMath::Point3D(width[1]*off, -L2*0.5 * relTorusPos[1], 0.), MWMath::axisAngle({1,0,0}, 90.0), COLORVPFLEXOR);
        mVPF2->bIsViaPoint = true;
        meshes.push_back(mVPF2);
        MeshMap[meshName] = mVPF2.get();

        /* meshName = "MeshVPE_PP" + prefN;
        auto mVPE2 = std::make_shared<SSEllipsoidMesh>(viaPointR, viaPointR, viaPointR, 
                meshName, body2, MWMath::Point3D(-width[1]*off, -L2*0.5 * relTorusPos[1], 0.), MWMath::axisAngle({1,0,0}, 90.0), COLORVPEXTENSOR);
        mVPE2->bIsViaPoint = true;
        meshes.push_back(mVPE2);
        MeshMap[meshName] = mVPE2.get(); */

        MWMath::Point3D jointPosRel2 = MWMath::Point3D(0, -L2*0.5, 0);
        auto joint2 = std::make_shared<SSJoint>("Joint_2", 
                                                jointPosRel2, 
                                                MWMath::RotMatrix3x3(),
                                                body2, 
                                                FJAngles[2], 
                                                MWMath::Point3D(0,0,1), // Bleibt Z
                                                numTimeSteps);
        tissues.push_back(joint2);
        
        double jointSize2 = width[1]*1.1;
        meshName = "Mesh_PIP" + prefN + "_Joint";
        auto jMesh2 = std::make_shared<SSEllipsoidMesh>(jointSize2, jointSize2, jointSize2, meshName, joint2, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 0.0), COLORJOINT);
        meshes.push_back(jMesh2);
        MeshMap[meshName] = jMesh2.get();

        meshName = "MeshVPE_PP" + prefN;
        auto mVPE2 = std::make_shared<SSEllipsoidMesh>(viaPointR, viaPointR, viaPointR, 
                meshName, joint2, MWMath::Point3D(-jointSize2*1.2, 0., 0.), MWMath::axisAngle({1,0,0}, 90.0), COLORVPEXTENSOR);
        mVPE2->bIsViaPoint = true;
        meshes.push_back(mVPE2);
        MeshMap[meshName] = mVPE2.get();

        // ==============================================================================
        // SEGMENT 3
        // ==============================================================================
        auto body3 = std::make_shared<SSBody>("MP" + prefN, MWMath::Point3D(0,-L3*0.5,0), MWMath::RotMatrix3x3(), joint2);
        tissues.push_back(body3);
        
        meshName = "Mesh_MP" + prefN;
        auto mesh3 = std::make_shared<SSEllipsoidMesh>(width[2], width[2], L3*0.5, 
                meshName, body3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), 
                COLORF3);
        if (bShowBody) meshes.push_back(mesh3);
        MeshMap[meshName] = mesh3.get();
        
        meshName = "MeshVPF_MP" + prefN;
        auto mVPF3 = std::make_shared<SSEllipsoidMesh>(viaPointR, viaPointR, viaPointR, 
                meshName, body3, MWMath::Point3D(width[2]*off, -L3*0.5 * relTorusPos[2], 0.), MWMath::axisAngle({1,0,0}, 90.0), COLORVPFLEXOR);
        mVPF3->bIsViaPoint = true;
        meshes.push_back(mVPF3);
        MeshMap[meshName] = mVPF3.get();

        /* meshName = "MeshVPE_MP" + prefN;
        auto mVPE3 = std::make_shared<SSEllipsoidMesh>(viaPointR, viaPointR, viaPointR, 
                meshName, body3, MWMath::Point3D(-width[2]*off, -L3*0.5 * relTorusPos[2], 0.), MWMath::axisAngle({1,0,0}, 90.0), COLORVPEXTENSOR);
        mVPE3->bIsViaPoint = true;
        meshes.push_back(mVPE3);
        MeshMap[meshName] = mVPE3.get(); */

        MWMath::Point3D jointPosRel3 = MWMath::Point3D(0, -L3*0.5, 0);
        auto joint3 = std::make_shared<SSJoint>("Joint_3", 
                                                jointPosRel3, 
                                                MWMath::RotMatrix3x3(),
                                                body3, 
                                                FJAngles[3], 
                                                MWMath::Point3D(0,0,1), // Bleibt Z
                                                numTimeSteps);
        tissues.push_back(joint3);

        double jointSize3 = width[2]*1.1;
        meshName = "Mesh_DIP" + prefN + "_Joint";
        auto jMesh3 = std::make_shared<SSEllipsoidMesh>(jointSize3, jointSize3, jointSize3, meshName, joint3, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 0.0), COLORJOINT);
        meshes.push_back(jMesh3);
        MeshMap[meshName] = jMesh3.get();

        meshName = "MeshVPE_MP" + prefN;
        auto mVPE3 = std::make_shared<SSEllipsoidMesh>(viaPointR, viaPointR, viaPointR, 
                meshName, joint3, MWMath::Point3D(-jointSize3*1.2, 0., 0.), MWMath::axisAngle({1,0,0}, 90.0), COLORVPEXTENSOR);
        mVPE3->bIsViaPoint = true;
        meshes.push_back(mVPE3);
        MeshMap[meshName] = mVPE3.get();

        // ==============================================================================
        // SEGMENT 4
        // ==============================================================================
        auto body4 = std::make_shared<SSBody>("DP" + prefN, MWMath::Point3D(0,-L4*0.5,0), MWMath::RotMatrix3x3(), joint3);
        tissues.push_back(body4);
        
        meshName = "Mesh_DP" + prefN;
        auto mesh4 = std::make_shared<SSEllipsoidMesh>(width[3], width[3], L4*0.5, 
                meshName, body4, MWMath::Point3D(0, 0, 0), MWMath::axisAngle({1,0,0}, 90.0), 
                COLORF4);
        if (bShowBody) meshes.push_back(mesh4);
        MeshMap[meshName] = mesh4.get();
    }

    // ==============================================================================
    // INITIALES UPDATE
    // ==============================================================================
    for (auto& m : meshes) {m->InitializeMesh();}
    rootSystem->update(0); 

    // ==============================================================================
    // MUSKEL
    // ==============================================================================
    
    // Muskelpunkte ebenfalls an neue Achsen angepasst: 
    // Offset auf der X-Achse, Länge (C) auf der negativen Y-Achse
    MWMath::Point3D startOffset = MWMath::Point3D(HL*0.03,HL*0.3, 0.0); //MWMath::Point3D(mesh1->B * 1.1, 0.0, 0.0);
    auto originMesh = dynamic_cast<SSEllipsoidMesh*>(MeshMap["Mesh_MC1"]);
    auto insertionMesh = dynamic_cast<SSEllipsoidMesh*>(MeshMap["Mesh_DP1"]);
    std::vector<std::string> refBodyNames = {"Mesh_WristFlexorCyl","MeshVP_MC1","Mesh_MC1_Joint","Mesh_PP1","MeshVP_PP1","Mesh_PIP1_Joint","Mesh_MP1","MeshVP_MP1","Mesh_DIP1_Joint","Mesh_DP1"};
    // MWMath::Point3D endOffset = MWMath::Point3D(insertionMesh->B * 1.1, 0.0, 0.0);
    MWMath::Point3D endOffset = MWMath::Point3D(insertionMesh->B * 1.1, -insertionMesh->C, 0.0);
    
    /* SSMuscle* flexor = new SSMuscle("Flexor", numPoints, 
        rootSystem.get(), startOffset, 
        getRefTissue("DP1", tissues), endOffset); */

    // Alle Hindernisse
    /* for(auto& m : meshes) {
        if (dynamic_cast<SSTorusMesh*>(m.get()) || m->bIsJointMesh){
            flexor->meshPtrs.push_back(m.get());
        }
        
        qDebug() << "      Added mesh '" << m->Name.c_str() << "' as obstacle for muscle path.";
        flexor->meshPtrs.push_back(m.get());
    } */

    /* for (auto& m : getRefMeshes(refBodyNames, MeshMap)) {
        flexor->meshPtrs.push_back(m);
        qDebug() << "      Added mesh '" << m->Name.c_str() << "' as obstacle for muscle path.";
    }

    flexor->createMusclePointsComplexPath();
    flexor->updateMusclePointsParents();
    muscles.push_back(flexor); */

    qDebug() << createMusclePathHand(rootSystem, meshes, muscles, cfg, numPoints).c_str();

    return "buildOHandModelOldExpandedViaX05";
}

