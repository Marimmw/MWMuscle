#pragma once

#include <vector>
#include <string>
#include <casadi/casadi.hpp>
#include "utils/MWMath.h"
#include "simpleSimulation/SSMeshes.h"
#include "simpleSimulation/SSTissue.h"




class MuscleNode : public SSTissue {
public:
    MuscleNode() = default;
    virtual ~MuscleNode() = default;

    MWMath::Point3D PositionLocal;  // Koordinate relativ zum Parent
    SSTissue* parentMesh = nullptr;
    bool bParentIsFixed = false;
    std::vector<double> lastEtas;

    std::vector<MWMath::Point3D> MNodeLocalSteps;
    std::vector<MWMath::Point3D> MNodeGlobalSteps;
    std::vector<MWMath::Point3D> MNodeInitialGuessSteps;
    std::vector<MWMath::Point3D> MNodeInitialGuessColorSteps;

    // for visualizing/checking local parametrization
    std::vector<MWMath::Point3D> MNodeParentPosSteps;
    std::vector<MWMath::RotMatrix3x3> MNodeParentRotSteps;

    std::vector<std::vector<double>> MNodeEtaSteps;

    // Updates PositionLocal based on the parent reference Body and the computed PositionGlobal
    void updateLocalFrame() {
        if (!parentMesh) return;
        MWMath::Point3D diff = PositionGlobal - parentMesh->Parent->PositionGlobal;
        PositionLocal = parentMesh->Parent->OrientationGlobal.transposed().transform(diff);
    }

    // Berechnet die NEUE globale Position, nachdem sich das Mesh bewegt hat
    MWMath::Point3D predictNewGlobal() {
        if (!parentMesh) return PositionGlobal;
        if(!parentMesh->Parent) return PositionGlobal;
        return parentMesh->Parent->PositionGlobal + (parentMesh->Parent->OrientationGlobal.transform(PositionLocal));
    }

    MWMath::Point3D getparentMeshColor() {
        double colorDarkerFactor = 0.3;
        if (bParentIsFixed) {
            if (!parentMesh) return MWMath::Point3D(0.5, 0.5, 0.5); // Grau, wenn kein Parent
            MWMath::Point3D baseColor = parentMesh->MeshColor;
            return MWMath::Point3D(baseColor.x * colorDarkerFactor,
                                   baseColor.y * colorDarkerFactor,
                                   baseColor.z * colorDarkerFactor);
        } 
        else {
            if (!parentMesh) return MWMath::Point3D(1.0, 1.0, 1.0); // Weiß, wenn kein Parent
            return parentMesh->MeshColor;
        }
    }

    void getMNodeInfoStep(int idx);
private:
};

class SSMuscle : public SSTissue {
public:
    SSMuscle() = default;
    SSMuscle(std::string name, int nodesCount, SSTissue* originParent, MWMath::Point3D originPointLocal, SSTissue* insertionParent, MWMath::Point3D insertionPointLocal) 
        : Name(name), MNodesCount(nodesCount), parentMeshOrigin(originParent), OriginPointLocal(originPointLocal), parentMeshInsertion(insertionParent), InsertionPointLocal(insertionPointLocal) {}
    virtual ~SSMuscle() = default;

    std::vector<MWMath::Point3D> MusclePointsGlobal;
    std::vector<MuscleNode> MNodes; // all muscle nodes including origin and insertion
    std::map<SSMesh*, int> attractorToNodeIndex;
    std::vector<MWMath::Point3D> CurrentViaPoints;
    std::vector<double> ViaPointsTolerances;
    bool bHasViaPoints = false;
    int MNodesCount = 3;
    std::vector<double> lastEtas;
    std::string Name = "";
    std::vector<SSMesh*> meshPtrs;
    double MuscleLength;

    MWMath::Point3D OriginPointGlobal;
    MWMath::Point3D OriginPointLocal;
    SSTissue* parentMeshOrigin = nullptr;
    MWMath::Point3D InsertionPointGlobal;
    MWMath::Point3D InsertionPointLocal;
    SSTissue* parentMeshInsertion = nullptr;

    std::vector<double> MuscleLengthSteps;

    double tooFarAwayThreshold = 1000000.0; // [meter]

    // DEBUG prints
    bool bParentDebug = false;
    bool bMeshDistanceDebug = true;

    std::vector<std::vector<MWMath::Point3D>> allOptimizedPoints; // Alle optimierten Punkte für alle Simulationsschritte
    void createMusclePoints();
    void updateMusclePointsParents();
    void updateMusclePointsParentsLocal();
    void updateAttractorNodes();
    double computeMuscleLength(bool addToHistory = false, int stepIdx = -1);

    int checkCollision(std::vector<SSMesh*> allToCheckMeshes = {});

    void getAttractorNodeInfo();
    void getViaPointNodeInfo();
    void getMuscleInfo(std::string prefix="    ");

    // data management
    void storeMNodesInitialGuess(int stepIdx, std::vector<MWMath::Point3D>& guessPath);
    void storeMNodesInitialGuessColors(int stepIdx, std::vector<MWMath::Point3D>& colorPath);
    void storeMNodesGlobalPositions();
    void initializeSimulationMuscle(int& steps);
    std::vector<MWMath::Point3D> gatherMNodesGlobal(int& step);
    void getAllMuscleMNodesStepValues(int totalSteps,   
        std::vector<std::vector<MWMath::Point3D>>& globalPoints, 
                                      std::vector<std::vector<MWMath::Point3D>>& localPoints, 
                                      std::vector<std::vector<MWMath::Point3D>>& initialGuessPoints,
                                      std::vector<std::vector<MWMath::Point3D>>& initialGuessColors,
                                      std::vector<std::vector<std::vector<double>>>& etaValues);
    
    private:
};