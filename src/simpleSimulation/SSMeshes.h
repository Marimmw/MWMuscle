#pragma once

#include <casadi/casadi.hpp>
#include "utils/MWMath.h"

#include "simpleSimulation/SSTissue.h"

#include <cmath>
#include <string>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


class SSMesh : public SSTissue, public std::enable_shared_from_this<SSMesh> {
public:
    SSMesh() = default;
    SSMesh(std::string name, 
           std::shared_ptr<SSTissue> parent, 
           MWMath::Point3D relPos, 
           MWMath::RotMatrix3x3 relRot, 
           MWMath::Point3D color) 
    {
        this->Name = name;
        this->Parent = parent; 
        this->Position2ParentRelInParentFrame = relPos;
        this->Orientation2ParentRel = relRot;
        this->MeshColor = color;
    }
    virtual ~SSMesh() {};
    bool bIsAttractor = false;
    bool bIsViaPoint = false;
    double MViaPointTolerance = 0.0; // Nur relevant, wenn bIsVia
    bool bIsJointMesh = false;
    

    // MWMath::Point3D PositionGlobal = MWMath::Point3D(0.0, 0.0, 0.0);
    // MWMath::RotMatrix3x3 OrientationGlobal = MWMath::RotMatrix3x3();
    // MWMath::Point3D MeshColor = MWMath::Point3D(0.7, 0.7, 0.7); // RGB Werte zwischen 0 und 1
    std::vector<MWMath::Point3D> MeshPointsGlobal; // Globale Punkte des Meshes
    std::vector<MWMath::RotMatrix3x3> allRMatrixGlobal; // Alle Rotationsmatrizen für alle Simulationsschritte
    std::vector<MWMath::Point3D> GlobalDiscreteMeshPoints;
    std::vector<MWMath::Point3D> AllScalerStepLists;

    void InitializeMesh();

    virtual casadi::MX constraintJacobian(casadi::MX gamma, casadi::MX q) = 0;
    virtual casadi::MX constraintDistance(casadi::MX gamma, casadi::MX q) = 0;

    virtual casadi::MX constraintJacobianLocal(casadi::MX gamma, casadi::MX q) {return casadi::MX();};
    virtual casadi::MX constraintDistanceLocal(casadi::MX gamma, casadi::MX q) {return casadi::MX();};

    virtual void getCasadiParentGPOs(casadi::MX& pos, casadi::MX& ori) override;
    virtual double getDistanceNumerically(MWMath::Point3D pGlobal, bool signedDistance=false) {return 0.0;};
    virtual void discretizeMesh(int discrCount) {return;};

    void updateMeshPosAndRot();

};

class SSEllipsoidMesh : public SSMesh {
public:
    SSEllipsoidMesh(double a, double b, double c) : A(a), B(b), C(c) {MViaPointTolerance = std::max({a, b, c}); }
    SSEllipsoidMesh(double a, double b, double c,
                    std::string name, 
                    std::shared_ptr<SSTissue> parent,
                    MWMath::Point3D relPos, 
                    MWMath::RotMatrix3x3 relRot, 
                    MWMath::Point3D color)
        : SSMesh(name, parent, relPos, relRot, color), // Ruft Basis auf
          A(a), B(b), C(c)                             // Initialisiert EIGENE Member
    {
        MViaPointTolerance = std::max({a, b, c});
    }
    ~SSEllipsoidMesh() override = default;

    casadi::MX constraintJacobian(casadi::MX gamma, casadi::MX q) override;
    casadi::MX constraintDistance(casadi::MX gamma, casadi::MX q) override;

    virtual casadi::MX constraintJacobianLocal(casadi::MX gamma, casadi::MX q) override;
    virtual casadi::MX constraintDistanceLocal(casadi::MX gamma, casadi::MX q) override;

    double getDistanceNumerically(MWMath::Point3D pGlobal, bool signedDistance=false) override;
    virtual void discretizeMesh(int discrCount) override;
    double A;
    double B;
    double C;
private:
};

class SSCylinderMesh : public SSMesh {
public:
    SSCylinderMesh(double r, double h) : Radius(r), Height(h) {MViaPointTolerance = r; }
    SSCylinderMesh(double r, double h,
                   std::string name, 
                   std::shared_ptr<SSTissue> parent,
                   MWMath::Point3D relPos, 
                   MWMath::RotMatrix3x3 relRot, 
                   MWMath::Point3D color)
        : SSMesh(name, parent, relPos, relRot, color), Radius(r), Height(h)
    {
        MViaPointTolerance = r;
    }
    ~SSCylinderMesh() override = default;   

    casadi::MX constraintJacobian(casadi::MX gamma, casadi::MX q) override;
    casadi::MX constraintDistance(casadi::MX gamma, casadi::MX q) override;
    double getDistanceNumerically(MWMath::Point3D pGlobal, bool signedDistance=false) override;
    double Radius;
    double Height;
private:
};

class SSTorusMesh : public SSMesh {
public:
    SSTorusMesh(double R, double r) : R(R), r(r) {MViaPointTolerance = R+r; }
    SSTorusMesh(double R, double r,
               std::string name, 
               std::shared_ptr<SSTissue> parent,
               MWMath::Point3D relPos, 
               MWMath::RotMatrix3x3 relRot, 
               MWMath::Point3D color)
        : SSMesh(name, parent, relPos, relRot, color), R(R), r(r)
    {
        MViaPointTolerance = R+r;
    }
    ~SSTorusMesh() override = default;   

    casadi::MX constraintJacobian(casadi::MX gamma, casadi::MX q) override;
    casadi::MX constraintDistance(casadi::MX gamma, casadi::MX q) override;
    double getDistanceNumerically(MWMath::Point3D pGlobal, bool signedDistance=false) override;
    double R;
    double r;

    virtual void discretizeMesh(int discrCount) override;

    double A = 1.0;
    double B = 1.0;
    double C = 1.0;
private:
};

class SSEllipticalTorusMesh : public SSMesh {
public:
    SSEllipticalTorusMesh(double R1, double R2, double RX, double RY) : R1(R1), R2(R2), RX(RX), RY(RY) {MViaPointTolerance = std::max({R1, R2, RX, RY}); }
    SSEllipticalTorusMesh(double R1, double R2, double RX, double RY,
                    std::string name, 
                    std::shared_ptr<SSTissue> parent,
                    MWMath::Point3D relPos, 
                    MWMath::RotMatrix3x3 relRot, 
                    MWMath::Point3D color)
        : SSMesh(name, parent, relPos, relRot, color), R1(R1), R2(R2), RX(RX), RY(RY)
    {
        MViaPointTolerance = std::max({R1, R2, RX, RY});
    }
    ~SSEllipticalTorusMesh() override = default;

    casadi::MX constraintJacobian(casadi::MX gamma, casadi::MX q) override;
    casadi::MX constraintDistance(casadi::MX gamma, casadi::MX q) override;

    double getDistanceNumerically(MWMath::Point3D pGlobal, bool signedDistance=false) override;
    virtual void discretizeMesh(int discrCount) override;

    // Extrusionsellipse 
    double R1; // ~x-achse
    double R2; // höhe=z-achse
    // pfadellipse (z-achse)
    double RX; // x-achse
    double RY; // y-achse

private:
};
