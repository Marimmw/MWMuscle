#pragma once
#include <casadi/casadi.hpp>
#include "SSMuscle.h"


#include <QDialog>
#include <QVBoxLayout>
#include <QSlider>
#include <QDebug>
#include <vector>

#include <memory>

#include "utils/MWMath.h"
#include "SSMeshes.h" // Für SSMesh / SSEllipsoidMesh
#include "SSTissue.h"
#include "SSBody.h"


class CasadiSystem {
public:
    // Konstruktor setzt jetzt direkt das System auf
    CasadiSystem(std::vector<SSMuscle*> muscles, int objType = 0, std::string version="new", std::string parametrizationType="global",
                    bool bUseCasGradient = false, bool bSumPhiEta = false);
    std::string CasadiSystemName;
    std::vector<std::string> SolverConvergenceMessages;
    std::vector<int> SolverConvergenceSteps;
    std::vector<SSMuscle*> m_muscles;

    std::vector<std::vector<double>> allParameterInputsAllSteps;
    std::vector<std::vector<std::string>> allParameterInputDescriptionsAllSteps;


    void setupCasadi();
    void solveStepX();

    void solveStepSum();
    void setupCasadiSum();
    
    void solveStep();



private:
    
    casadi::Function m_solver;
    std::string Version;
    std::string ParametrizationType; // global or local
    bool bUseOwnGradient;
    
    
    
    // Wir müssen speichern, welcher Muskel wo im großen x-Vektor steht
    struct MuscleOffsets {
        int x_start;      // Index im x-Vektor
        int p_start;      // Index im p-Vektor (Origin/Insertion)
        int num_nodes;    
        int num_etas;
    };
    std::vector<MuscleOffsets> m_offsets;
    int m_total_x_vars = 0;


    int M; // Anzahl Muskeln
    //int K; // unused currently // Anzahl Segmente (Punkte = K+1)
    SSMesh* generic_mesh; //SSEllipsoidMesh* ellipsoid;
    std::vector<SSMesh*> generic_meshes;
    bool globalComputation = true;
    int objType; // 0 = no(=0). min Length
    bool bSumPhiEta;
    bool bUseWarmstartEtas = false;
    double WarmstartEtaScaling = 1.0;
    int maxIterations = 2000;
    double maxTol = 1e-5;
    double ELTolerance = 0.0;
    //std::string hessianApproximation = "limited-memory"; // "limited-memory" or "exact"
    casadi::Dict opts;

    casadi::Function solver;
    
    // Dimensionen wie in MATLAB
    int num_musclenodes;
    int num_contactmultiplier;
    int num_total_vars;
};


