
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
    //qDebug() << "Global Position: " << Name.c_str() << "=" <<  PositionGlobal.print().c_str();
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
            //qDebug() << QString("   ").repeated(SystemLayer) << "|      " << QString::fromStdString(Name) <<  " Step: " << step << ", New Angle: " << newAngle;
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
    //qDebug() << "Global Position: " << Name.c_str() << "=" <<  PositionGlobal.print().c_str();
    return 0;
}


// -----------------------------------------------------

int SSTranslationalJoint::update(int step)
{  
    if (bDebug == true && step == 0) qDebug() << QString("   ").repeated(SystemLayer) << "|->" << " Updating Translational Joint: " <<  QString::fromStdString(Name);
    
    double newTranslation = 0.0;
    MWMath::Point3D axis = TranslationAxis;
    
    // 1. Berechne die aktuelle Auslenkung
    if (TranslationSteps.empty()) {
        double denominator = (TotalSteps > 1) ? (double)(TotalSteps - 1) : 1.0;
        double progress = (double)step / denominator;
        newTranslation = MaxTranslation * progress;
    }
    else {
        if (step < TranslationSteps.size()) {
            newTranslation = TranslationSteps[step];
        }
        else {
            newTranslation = TranslationSteps.back(); // Letzten Wert beibehalten
        }
    }

    CurrentTranslation = newTranslation;
    DoneTranslationSteps.push_back(newTranslation);

    // Vektor der Translation (Lokal)
    MWMath::Point3D localTranslation = axis * newTranslation;
    MWMath::Point3D localHalfTranslation = axis * (newTranslation * 0.5);

    // 2. Globale Transformationen berechnen
    if (Parent) {
        // Orientierung wird 1:1 vom Parent übernommen + Initiales Offset (keine Bewegung)
        OrientationGlobal = Parent->OrientationGlobal * Orientation2ParentRel;
        
        // Die Basis-Position des Gelenks an sich
        MWMath::Point3D basePos = Parent->PositionGlobal + Parent->OrientationGlobal.transform(Position2ParentRelInParentFrame);
        
        // Die neue Position wird entlang der lokalen Achse des Gelenks verschoben
        PositionGlobal = basePos + OrientationGlobal.transform(localTranslation);
        JointHalfTranslation = basePos + OrientationGlobal.transform(localHalfTranslation);
    }
    else {
        OrientationGlobal = Orientation2ParentRel;
        PositionGlobal = Position2ParentRelInParentFrame + OrientationGlobal.transform(localTranslation);
        JointHalfTranslation = Position2ParentRelInParentFrame + OrientationGlobal.transform(localHalfTranslation);
    }

    // 3. Update Children
    for (auto& child : Children) {
        child->update(step);
    }

    // 4. Update Meshes
    for (auto& mesh : Meshes) {
        mesh->OrientationGlobal = OrientationGlobal * mesh->Orientation2ParentRel;
        
        if (mesh->bIsJointMesh) {
            // Wenn das Mesh das Gelenk selbst repräsentiert, platzieren wir es in der Mitte der Auslenkung!
            mesh->PositionGlobal = JointHalfTranslation + mesh->OrientationGlobal.transform(mesh->Position2ParentRelInParentFrame);
            // Optional: Wenn du einen Zylinder hast, der "mitwachsen" soll, könntest du hier mesh->C = newTranslation setzen.
        }
        else {
            // Reguläre Meshes wandern voll mit der Spitze des Gelenks mit
            mesh->PositionGlobal = PositionGlobal + mesh->OrientationGlobal.transform(mesh->Position2ParentRelInParentFrame);
        }
    }

    return 0;
}



