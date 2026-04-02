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
    std::vector<double> AngleSteps; 
    std::vector<double> DoneAngleSteps;

    MWMath::RotMatrix3x3 JointHalfRotation = MWMath::RotMatrix3x3();

    virtual int update(int step=0) override;
private:
};


class SSTranslationalJoint : public SSTissue {
public:
    SSTranslationalJoint() = default;
    SSTranslationalJoint(std::string name, 
                         MWMath::Point3D relPos, 
                         MWMath::RotMatrix3x3 relRot,
                         std::shared_ptr<SSTissue> parent,
                         double maxTranslation, // Maximale Auslenkungsdistanz (statt Winkel)
                         MWMath::Point3D axis,  // Richtung (z.B. {1,0,0} für lokale X-Achse)
                         int totalSteps)
    {
        Name = name;
        Position2ParentRelInParentFrame = relPos;
        Orientation2ParentRel = relRot;
        Parent = parent;
        
        TotalSteps = totalSteps;
        
        TranslationAxis = axis.normed(); // Sicherstellen, dass es ein Einheitsvektor ist
        MaxTranslation = maxTranslation;
        CurrentTranslation = 0.0; 

        if (Parent) {
            Parent->Children.push_back(std::shared_ptr<SSTissue>(this));
            SystemLayer = Parent->SystemLayer + 1;
        }
        else{
            SystemLayer = 0; 
        }
    }
    ~SSTranslationalJoint() override = default;
    
    double CurrentTranslation; 
    double MaxTranslation; 
    MWMath::Point3D TranslationAxis;
    int TotalSteps;
    
    std::vector<double> TranslationSteps; // Optional: Vorgegebene Trajektorie (analog zu AngleSteps)
    std::vector<double> DoneTranslationSteps;

    // Für die Visualisierung: Der halbe Weg, z.B. für einen Teleskop-Zylinder
    MWMath::Point3D JointHalfTranslation = MWMath::Point3D(0,0,0);

    virtual int update(int step=0) override;
};


