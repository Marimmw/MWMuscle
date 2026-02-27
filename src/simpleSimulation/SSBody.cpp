
#include "simpleSimulation/SSBody.h"

int SSBody::update(int step)
{
    if (bDebug == true && step == 0) qDebug() << QString("   ").repeated(SystemLayer) << "|->" << " Updating Body: " <<  QString::fromStdString(Name);

    // 1. Parent-Transformation anwenden
    if (Parent) {
        // Position
        PositionGlobal = Parent->PositionGlobal + Parent->OrientationGlobal.transform(Position2ParentRelInParentFrame);
        // Orientierung
        OrientationGlobal = Parent->OrientationGlobal * Orientation2ParentRel;
    }
    else {
        PositionGlobal = Position2ParentRelInParentFrame;
        OrientationGlobal = Orientation2ParentRel;
    }
    for (auto& mesh : Meshes) {
        // Mesh-Position und -Orientierung aktualisieren
        /* mesh->OrientationGlobal = OrientationGlobal * mesh->Orientation2ParentRel;
        mesh->PositionGlobal = PositionGlobal + mesh->OrientationGlobal.transform(mesh->Position2ParentRelInParentFrame); */
        mesh->updateMeshPosAndRot();
        if (bDebug == true && step == 0) qDebug() << QString("   ").repeated(SystemLayer) << "|      " << " Mesh: " << QString::fromStdString(mesh->Name);
    }
    
    for (auto& child : Children) {
        child->update(step);
    }
    
    return 0;
}


// -----------------------------------------------------

int SSJoint::update(int step)
{  
    if (bDebug == true && step == 0) qDebug() << QString("   ").repeated(SystemLayer) << "|->" << " Updating Joint: " <<  QString::fromStdString(Name);
    MWMath::RotMatrix3x3 jointRotation = MWMath::RotMatrix3x3();
    MWMath::RotMatrix3x3 jointHalfRotation = MWMath::RotMatrix3x3();
    
            double newAngle = 0.0;
            MWMath::Point3D axis = RotationAxes;
            if (AngleSteps.empty()) {
                double denominator = (TotalSteps > 1) ? (double)(TotalSteps - 1) : 1.0;
                double progress = (double)step / denominator;
                newAngle = MaxAngles * progress;
            }
            else {
                if (step < AngleSteps.size()) {
                    newAngle = AngleSteps[step];
                }
                else {
                    newAngle = AngleSteps.back(); // Wenn mehr Schritte als definiert, letzten Winkel beibehalten
                }
            }
            jointRotation = jointRotation * MWMath::axisAngle(axis, newAngle);
            jointHalfRotation = jointHalfRotation * MWMath::axisAngle(axis, newAngle * 0.5);
            DoneAngleSteps.push_back(newAngle);
        
    // 1. Parent-Transformation anwenden
    if (Parent) {
        // Position
        PositionGlobal = Parent->PositionGlobal + Parent->OrientationGlobal.transform(Position2ParentRelInParentFrame);
        // Orientierung
        OrientationGlobal = Parent->OrientationGlobal * Orientation2ParentRel * jointRotation;
        JointHalfRotation = Parent->OrientationGlobal * Orientation2ParentRel * jointHalfRotation;
    }
    else {
        PositionGlobal = Position2ParentRelInParentFrame;
        OrientationGlobal = Orientation2ParentRel * jointRotation;
        JointHalfRotation = Orientation2ParentRel * jointHalfRotation;
    }

    for (auto& child : Children) {
        child->update(step);
    }
    for (auto& mesh : Meshes) {
        if (mesh->bIsJointMesh) {
            // Gelenk-Mesh um die halbe Rotation drehen
            
        }
        mesh->OrientationGlobal = JointHalfRotation * mesh->Orientation2ParentRel;
        mesh->PositionGlobal = PositionGlobal + mesh->OrientationGlobal.transform(mesh->Position2ParentRelInParentFrame);
        //mesh->updateMeshPosAndRot();
    }
    return 0;
}
