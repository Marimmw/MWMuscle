#pragma once

#include <vector>
#include <string>
#include <casadi/casadi.hpp>
#include "utils/MWMath.h"
#include "simpleSimulation/SSMeshes.h"
#include "simpleSimulation/SSTissue.h"

class SSJoint; // Forward Declaration
struct MomentArmFoJoint {
    std::string JointName;
    SSJoint* Joint = nullptr;
    std::vector<double> MomentArmValues; // Momentenarm-Werte über die Zeit
    std::vector<double> JointAngleValues; // Gelenkwinkel-Werte über die Zeit
};


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

    std::vector<std::vector<double>> MNodeEtaSteps; // [ [eta_mesh1, eta_mesh2, ...], [eta_mesh1, eta_mesh2, ...], ... ]
    std::vector<std::vector<double>> MNodePhiSteps; // [ [phi_mesh1, phi_mesh2, ...], [phi_mesh1, phi_mesh2, ...], ... ] -> Phi-Werte für alle Meshes für jeden Schritt

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

    std::vector<double> computePhiForEachMeshAtStep(int stepIdx, const std::vector<SSMesh*>& meshPtrs) {
        // used after the simulaiton for checking and plotting constraints
        // computes a vector of phi values for this MuscleNode at a specific stepIdx for all referenced muscle-meshes
        
        int meshCount = meshPtrs.size();
        std::vector<double> phiValues(meshCount, 0.0);
        for (size_t i = 0; i < meshCount; ++i) {
            SSMesh* mesh = meshPtrs[i];
            double phi = 0.0;
            if (stepIdx < MNodeGlobalSteps.size() && stepIdx < mesh->MeshPointsGlobal.size() && stepIdx < mesh->allRMatrixGlobal.size()) {
                // 1. Daten holen
                MWMath::Point3D nodePos = MNodeGlobalSteps[stepIdx];
                MWMath::Point3D meshPosGlobal = mesh->MeshPointsGlobal[stepIdx];
                MWMath::RotMatrix3x3 meshRotGlobal = mesh->allRMatrixGlobal[stepIdx];

                // 2. Numerische Inputs für CasADi vorbereiten (Muss exakt deinem 'p'-Vektor aus solveStepSum entsprechen!)
                std::vector<double> gamma_val = {nodePos.x, nodePos.y, nodePos.z};
                
                std::vector<double> q_val;
                q_val.push_back(meshPosGlobal.x);
                q_val.push_back(meshPosGlobal.y);
                q_val.push_back(meshPosGlobal.z);
                
                // Rotation (Spaltenweise / Column-Major, wie in deiner solveStepSum)
                for (int col = 0; col < 3; ++col) {
                    for (int row = 0; row < 3; ++row) {
                        q_val.push_back(meshRotGlobal.m[row][col]);
                    }
                }

                // 3. CasADi Auswertung
                using namespace casadi;
                
                // Symbolische Variablen erstellen
                MX sym_gamma = MX::sym("gamma", 3);
                MX sym_q = MX::sym("q", 12);
                
                // Deine Mesh-Funktion aufrufen, um den MX-Graphen zu generieren
                MX phi_expr = mesh->constraintDistance(sym_gamma, sym_q);
                
                // Daraus eine kompilierte/auswertbare CasADi-Funktion machen
                Function phi_func("phi_eval", {sym_gamma, sym_q}, {phi_expr});
                
                // Numerische Werte übergeben (als casadi::DM)
                std::vector<DM> arg = {gamma_val, q_val};
                std::vector<DM> res = phi_func(arg); // Auswerten!
                
                // Das Ergebnis ist eine 1x1 Matrix, wir casten sie zu double
                phi = double(res.at(0));
            }
            phiValues[i] = phi;
        }
        return phiValues;
    }
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
    std::vector<MomentArmFoJoint> MuscleMomentArmResults;

    double tooFarAwayThreshold = 1000000.0; // [meter]

    // DEBUG prints
    bool bParentDebug = false;
    bool bMeshDistanceDebug = false;

    std::vector<std::vector<MWMath::Point3D>> allOptimizedPoints; // Alle optimierten Punkte für alle Simulationsschritte
    void createMusclePoints();
    void createMusclePointsComplexPath();
    void updateMusclePointsParents();
    void updateMusclePointsParentsLocal();
    void updateAttractorNodes();
    double computeMuscleLength(bool addToHistory = false, int stepIdx = -1);
    double computeMomentArm(int stepIdx = -1); // "-1" for all steps, otherwise specific step

    int checkCollision(std::vector<SSMesh*> allToCheckMeshes = {});

    // std::vector<SSMesh*> getReferencedMeshes(std::vector<std::string> meshNames);
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

    void exportMomentArms();
    
    private:
};