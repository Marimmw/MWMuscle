#pragma once
#include "utils/MWMath.h"
#include "simpleSimulation/SSTissue.h"
#include "simpleSimulation/SSMeshes.h"

class SSBody : public SSTissue {
public:
    SSBody() = default;
    SSBody(std::string name, 
           MWMath::Point3D relPos, 
           MWMath::RotMatrix3x3 relRot = MWMath::RotMatrix3x3(), // Identität als Default
           std::shared_ptr<SSTissue> parent = nullptr) 
    {
        Name = name;
        Position2ParentRelInParentFrame = relPos;
        Orientation2ParentRel = relRot;
        Parent = parent;
        
        if (Parent) {
            Parent->Children.push_back(std::shared_ptr<SSTissue>(this)); // Achtung: Shared Ptr Management hier ist tricky, besser extern machen!
            SystemLayer = Parent->SystemLayer + 1; // Child ist eine Schicht tiefer als Parent
        }
        else{
            SystemLayer = 0; // Root ist Schicht 0
        }
    }
    ~SSBody() override = default;

    
    virtual int update(int step=0) override;
private:
};

class SSJoint : public SSTissue {
public:
    SSJoint() = default;
    SSJoint(std::string name, 
            MWMath::Point3D relPos, 
            MWMath::RotMatrix3x3 relRot,
            std::shared_ptr<SSTissue> parent,
            double maxAngleDeg, 
            MWMath::Point3D axis, 
            int totalSteps)
    {
        Name = name;
        Position2ParentRelInParentFrame = relPos;
        Orientation2ParentRel = relRot;
        Parent = parent;
        
        TotalSteps = totalSteps;
        
        // Setup für EINE Achse (Index 0) - erweiterbar für mehr DoF
        RotationAxes = axis;
        MaxAngles = maxAngleDeg;
        CurrentAngles = 0.0; // Start bei 0
        if (Parent) {
            Parent->Children.push_back(std::shared_ptr<SSTissue>(this)); // Achtung: Shared Ptr Management hier ist tricky, besser extern machen!
            SystemLayer = Parent->SystemLayer + 1; // Child ist eine Schicht tiefer als Parent
        }
        else{
            SystemLayer = 0; // Root ist Schicht 0
        }
    }
    ~SSJoint() override = default;
    
    double CurrentAngles; // In Grad
    double MaxAngles; // In Grad
    MWMath::Point3D RotationAxes;
    int TotalSteps;

    MWMath::RotMatrix3x3 JointHalfRotation = MWMath::RotMatrix3x3();

    virtual int update(int step=0) override;
private:
};