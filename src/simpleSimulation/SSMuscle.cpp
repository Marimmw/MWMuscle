#include "SSMuscle.h"

#include <QDebug>
#include <QString>

void SSMuscle::createMusclePoints()
{
    MusclePointsGlobal.clear();
    MNodes.clear();
    MusclePointsGlobal.reserve(MNodesCount);
    MNodes.reserve(MNodesCount);

    InsertionPointGlobal = parentMeshInsertion ? parentMeshInsertion->PositionGlobal + parentMeshInsertion->OrientationGlobal.transform(InsertionPointLocal) : InsertionPointLocal;
    OriginPointGlobal = parentMeshOrigin ? parentMeshOrigin->PositionGlobal + parentMeshOrigin->OrientationGlobal.transform(OriginPointLocal) : OriginPointLocal;
    MWMath::Point3D direction = (InsertionPointGlobal - OriginPointGlobal) * (1/distance(InsertionPointGlobal, OriginPointGlobal));
    for (int i = 0; i < MNodesCount; ++i) {
        double t = static_cast<double>(i) / (MNodesCount - 1); // Richtiges t von 0 bis 1
        MWMath::Point3D point;
        point.x = OriginPointGlobal.x + t * (InsertionPointGlobal.x - OriginPointGlobal.x);
        point.y = OriginPointGlobal.y + t * (InsertionPointGlobal.y - OriginPointGlobal.y);
        point.z = OriginPointGlobal.z + t * (InsertionPointGlobal.z - OriginPointGlobal.z);
        MusclePointsGlobal.push_back(point);

        // Node initialisieren
        MuscleNode node;
        node.PositionGlobal = point;
        if (i == 0 || i == MNodesCount - 1)
        {
            node.bParentIsFixed = true;
            node.parentMesh = (i == 0) ? parentMeshOrigin : parentMeshInsertion; // Wird im ersten Simulationsschritt gesetzt
            MWMath::Point3D parentPosGlob = MWMath::Point3D(0,0,0);
            MWMath::RotMatrix3x3 parentRGlob = MWMath::RotMatrix3x3();
            if (i==0){
                parentPosGlob = parentMeshOrigin ? parentMeshOrigin->PositionGlobal : MWMath::Point3D(0,0,0);
                parentRGlob = parentMeshOrigin ? parentMeshOrigin->OrientationGlobal : MWMath::RotMatrix3x3();
                OriginPointGlobal = parentPosGlob + parentRGlob.transform(OriginPointLocal);
                node.PositionLocal = OriginPointLocal;
            }
            else{
                parentPosGlob = parentMeshInsertion ? parentMeshInsertion->PositionGlobal : MWMath::Point3D(0,0,0);
                parentRGlob = parentMeshInsertion ? parentMeshInsertion->OrientationGlobal : MWMath::RotMatrix3x3();
                InsertionPointGlobal = parentPosGlob + parentRGlob.transform(InsertionPointLocal);
                node.PositionLocal = InsertionPointLocal;
            }
        }
        else
        {
            node.PositionGlobal = point;
            node.PositionLocal = point;
            node.bParentIsFixed = false;
            node.parentMesh = nullptr; // Wird im ersten Simulationsschritt gesetzt 
        }
        
        MNodes.push_back(node);
    }
}

void SSMuscle::updateMusclePointsParents()
{
    for (size_t k = 1; k < MNodes.size()-1; ++k) {
            if (!MNodes[k].bParentIsFixed) {
                auto& node = MNodes[k];
                
                // Optimierte Position übernehmen
                node.PositionGlobal = MusclePointsGlobal[k];

                // Nächstes Mesh finden
                double minDist = 1e10;
                SSMesh* bestParent = nullptr;
                
                if (bMeshDistanceDebug) {
                    qDebug() << "   Node " << k << "(" << node.PositionGlobal.x << "," << node.PositionGlobal.y << "," << node.PositionGlobal.z << "):";
                }
                for (auto& m : meshPtrs) {
                    // avoid cylinder meshes (just because i think last time might didnt work yet)
                    if (dynamic_cast<SSCylinderMesh*>(m)) {
                        continue;
                    }
                    if (dynamic_cast<SSTorusMesh*>(m)) {
                        continue;
                    }
                    // avoid via point meshes -> because for parent choice the distance to the COM is currently used, and via points (if not moving but near to a joint tend to produce points inside meshes for the next step)
                    /* if (m->bIsViaPoint) {
                        continue;
                    } */
                    if(m->bIsJointMesh){
                        continue;
                    }

                    double d = m->getDistanceNumerically(node.PositionGlobal);
                    if (bMeshDistanceDebug) {qDebug() << "       Prüfe Mesh" << QString::fromStdString(m->Name) << "(" << m->MeshColor.x << "," << m->MeshColor.y << "," << m->MeshColor.z << ") Distanz =" << d << "m";}
                    if (d < minDist) {
                        minDist = d;
                        bestParent = m;
                    }
                }

                /* // Schwellwert: Wenn zu weit weg, kein Parent zuweisen
                if (minDist > tooFarAwayThreshold) {
                    node.parentMesh = nullptr;
                    node.PositionLocal = node.PositionGlobal;
                    if (bParentDebug) {qDebug() << "  Node" << k << "-> Kein Parent (Distanz:" << minDist << "m)";}
                } else {
                    node.parentMesh = bestParent;
                    node.PositionLocal = node.parentMesh->OrientationGlobal.transposed().transform(
                        node.PositionGlobal - node.parentMesh->PositionGlobal);
                    if (bParentDebug) {qDebug() << "  Node" << k << "-> Parent:"  << QString::fromStdString(bestParent->Name) << "(Distanz:" << minDist << "m)";}
                } */
                // Schwellwert: Wenn zu weit weg, kein Parent zuweisen
                if (minDist > tooFarAwayThreshold) {
                    node.parentMesh = nullptr;
                    node.PositionLocal = node.PositionGlobal;
                    if (bParentDebug) {qDebug() << "  Node" << k << "-> Kein Parent (Distanz:" << minDist << "m)";}
                } else {
                    node.parentMesh = bestParent;
                    MWMath::Point3D parentOfMeshPos = bestParent->Parent->PositionGlobal;
                    MWMath::RotMatrix3x3 parentOfMeshR = bestParent->Parent->OrientationGlobal;

                    node.PositionLocal = parentOfMeshR.transposed().transform(
                        node.PositionGlobal - parentOfMeshPos);
                    if (bParentDebug) {qDebug() << "  Node" << k << "-> Parent:"  << QString::fromStdString(bestParent->Name) << "(Distanz:" << minDist << "m)";}
                }

                bool bBigDebug = false;
                if (bBigDebug){
                    QString parentName = node.parentMesh ? QString::fromStdString(node.parentMesh->Name) : "Kein Parent";
                    QString possibleParents = "";
                    for (auto& m : meshPtrs) {
                        double d = std::abs(m->getDistanceNumerically(node.PositionGlobal));
                        possibleParents += QString::fromStdString(m->Name) + "(d=" + QString::number(d, 'f', 3) + "m), ";
                    }
                    qDebug() << "  Node" << k << "PositionGlobal: (" << node.PositionGlobal.x << "," << node.PositionGlobal.y << "," << node.PositionGlobal.z << ")"
                             << "       Parent:" << parentName
                             // << "       PossibleParents: " << possibleParents
                             << "       PositionLocal: (" << node.PositionLocal.x << "," << node.PositionLocal.y << "," << node.PositionLocal.z << ")";
                }
            }
        }
}

void SSMuscle::updateMusclePointsParentsLocal()
{
    for (size_t k = 1; k < MNodes.size()-1; ++k) {
            if (!MNodes[k].bParentIsFixed) {
                auto& node = MNodes[k];
                
                // Optimierte Position übernehmen
                node.PositionGlobal = MusclePointsGlobal[k];

                // Nächstes Mesh finden
                double minDist = 1e10;
                SSMesh* bestParent = nullptr;
                if (bMeshDistanceDebug) {
                    qDebug() << "   Node " << k << "(" << node.PositionGlobal.x << "," << node.PositionGlobal.y << "," << node.PositionGlobal.z << "):";
                }
                for (auto& m : meshPtrs) {
                    // avoid cylinder meshes (just because i think last time might didnt work yet)
                    if (dynamic_cast<SSCylinderMesh*>(m)) {
                        continue;
                    }
                    /* if (dynamic_cast<SSTorusMesh*>(m)) {
                        continue;
                    } */
                    // avoid via point meshes -> because for parent choice the distance to the COM is currently used, and via points (if not moving but near to a joint tend to produce points inside meshes for the next step)
                    /* if (m->bIsViaPoint) {
                        continue;
                    } */
                    if(m->bIsJointMesh){
                        continue;
                    }

                    double d = m->getDistanceNumerically(node.PositionGlobal);
                    if (bMeshDistanceDebug) {qDebug() << "       Prüfe Mesh" << QString::fromStdString(m->Name) << "(" << m->MeshColor.x << "," << m->MeshColor.y << "," << m->MeshColor.z << ") Distanz =" << d << "m";}
                    if (d < minDist) {
                        minDist = d;
                        bestParent = m;
                    }
                }


                
                node.parentMesh = bestParent;
                // WICHTIG: Wir wollen die lokale Position relativ zum BODY (Grandparent),
                // nicht zum Mesh, weil der Solver den Body nutzt!
                
                SSTissue* referenceBody = nullptr;
                
                // Hat das Mesh einen Parent (Body/Joint)?
                if (bestParent->Parent) {
                    referenceBody = bestParent->Parent.get();
                } else {
                    // Fallback, falls Mesh keinen Parent hat (sollte mit main-Fix nicht passieren)
                    //referenceBody = bestParent; 
                    qDebug() << "[FEHLER] Mesh" << QString::fromStdString(bestParent->Name) << "hat keinen Parent Body/Joint!";
                    referenceBody = nullptr;
                } 

                // Berechnung relativ zum Referenz-Body
                MWMath::Point3D diff = node.PositionGlobal - referenceBody->PositionGlobal;
                node.PositionLocal = referenceBody->OrientationGlobal.transposed().transform(diff); // node.PositionLocal = // referenceBody->OrientationGlobal.transform(diff);
                if (bParentDebug) {qDebug() << "  Node" << k << "-> Parent:"  << QString::fromStdString(bestParent->Name) << "(Distanz:" << minDist << "m)";}
                
            }
        }
}

void SSMuscle::updateAttractorNodes()
{
    // updates the map from attractor mesh to node index -> which node (index) is influenced by the attractor(s)
    for (SSMesh* mesh : meshPtrs){
        if (mesh->bIsAttractor){
            double minDist = 1e10;
            int bestNodeIdx = -1;
            for (int k = 0; k < MNodes.size(); ++k) {
                double d = MWMath::distance(MNodes[k].PositionGlobal, mesh->PositionGlobal);
                if (d < minDist) {
                    minDist = d;
                    bestNodeIdx = k;
                }
            }
            attractorToNodeIndex[mesh] = bestNodeIdx;
        }
    }
}

// computes the muscle length; default step=-1 uses current positions
double SSMuscle::computeMuscleLength(bool addToHistory, int stepIdx)
{
    double musclelength = 0.0;
    if (stepIdx < 0) {
        for (size_t i = 1; i < MNodes.size(); ++i) {
            musclelength += MWMath::distance(MNodes[i-1].PositionGlobal, MNodes[i].PositionGlobal);
        }
    }
    else{
        for (size_t i = 1; i < MNodes.size(); ++i) {
            if (stepIdx < MNodes[i].MNodeGlobalSteps.size() && stepIdx < MNodes[i-1].MNodeGlobalSteps.size()) {
                musclelength += MWMath::distance(MNodes[i-1].MNodeGlobalSteps[stepIdx], MNodes[i].MNodeGlobalSteps[stepIdx]);
            }
            else {
                qDebug() << "Fehler: Schritt" << stepIdx << "übersteigt gespeicherte globale Schritte für Node" << i << "oder" << i-1;
            }
        }
    }
    if (addToHistory){
        MuscleLengthSteps.push_back(musclelength);
    }
    MuscleLength = musclelength;
    return musclelength;
}

int SSMuscle::checkCollision(std::vector<SSMesh *> allToCheckMeshes)
{
    if (allToCheckMeshes.empty()) {
        allToCheckMeshes = meshPtrs; // if no meshes provided, check against own meshes
    }

    qDebug() << "      Checking Collisions - " << allToCheckMeshes.size() << " Meshes";
    // checks collision -> 0=no collision, 1=collision with "own muscle mesh", 2=collision with other mesh
    for (auto* m : allToCheckMeshes) {
        if (m->bIsViaPoint) {
            continue; // via point meshes are not considered for collision
        }
        for (size_t k = 0; k < MNodes.size(); ++k) {
            double d = m->getDistanceNumerically(MNodes[k].PositionGlobal, true);
            if (d < 0.0) {
                if (std::find(meshPtrs.begin(), meshPtrs.end(), m) != meshPtrs.end()) {
                    // collision with own muscle mesh
                    qDebug() << "        [WARNING] Collision detected with own muscle mesh:" << QString::fromStdString(m->Name) << "at Node Index:" << k;
                    return 1;
                }
                else {
                    // collision with other mesh
                    qDebug() << "        [WARNING] Collision detected with OTHER mesh:" << QString::fromStdString(m->Name) << "at Node Index:" << k;
                    return 2;
                }
            }
        }
    }
    return 0;
}

void SSMuscle::getAttractorNodeInfo()
{
    for (auto& pair : attractorToNodeIndex) {
        SSMesh* mesh = pair.first;
        int nodeIdx = pair.second;
        qDebug() << "   Attractor Mesh:" << QString::fromStdString(mesh->Name) << "-> Node Index:" << nodeIdx;
    }
}

void SSMuscle::getViaPointNodeInfo()
{
    if (bHasViaPoints == false){return;}

    qDebug() << "   Checking Muscle: " << QString::fromStdString(Name);
    for (auto* m : meshPtrs) {
        if (m->bIsViaPoint) {
            qDebug() << "       ViaPoint Mesh:" << QString::fromStdString(m->Name) << "(vptol=" << m->MViaPointTolerance << ")";
            double minDist = 1e10;
            int minDistIdx;
            for (size_t k = 0; k < MNodes.size(); ++k) {
                double d = MWMath::distance(MNodes[k].PositionGlobal, m->PositionGlobal);
                if (d < minDist) {
                    minDist = d;
                    minDistIdx = k;
                }
                // qDebug() << "       " << "Node Index:" << k << "; Distanz:" << d;
            }
            qDebug() << "           Nächstgelegener Node Index:" << minDistIdx << "; Distanz:" << minDist;
            if (minDist > m->MViaPointTolerance) {
                qDebug() << "           [WARNING] ViaPoint constraint violated at Node Index:" << minDistIdx << "-> " << minDist;
            }
        }
    } 

    /* for (auto* m : meshPtrs) {
        if (m->bIsViaPoint) {
            qDebug() << "       ViaPoint Mesh:" << QString::fromStdString(m->Name) << "(vptol=" << m->MViaPointTolerance << ")";
            double minDist = 1e10;
            int minDistIdx;
            for (size_t k = 0; k < MNodes.size()-1; ++k) {
                double lineLength = MWMath::distance(MNodes[k].PositionGlobal, MNodes[k+1].PositionGlobal);
                MWMath::Point3D dirNorm = (MNodes[k+1].PositionGlobal - MNodes[k].PositionGlobal) * (1/lineLength);
                MWMath::Point3D middleOfLine = MNodes[k].PositionGlobal + dirNorm * lineLength*0.5;
                double d = MWMath::distance(middleOfLine, m->PositionGlobal);
                if (d < minDist) {
                    minDist = d;
                    minDistIdx = k;
                }
                // qDebug() << "       " << "Node Index:" << k << "; Distanz:" << d;
            }
            qDebug() << "           Nächstgelegener LineSegment:" << minDistIdx << "->" << minDistIdx+1 << "; Distanz:" << minDist;
            if (minDist > m->MViaPointTolerance) {
                qDebug() << "           [WARNING] ViaPoint constraint violated at Node Index:" << minDistIdx << "-> " << minDist;
            }
        }
    } */
}

void SSMuscle::getMuscleInfo(std::string prefix)
{
    qDebug() << QString::fromStdString(prefix) + "Muscle Name:" << QString::fromStdString(Name);
    qDebug() << QString::fromStdString(prefix) + "   Origin Point Local (" << QString::fromStdString(parentMeshOrigin->Name) << "): (" << OriginPointLocal.x << "," << OriginPointLocal.y << "," << OriginPointLocal.z << ")";
    qDebug() << QString::fromStdString(prefix) + "   Insertion Point Local (" << QString::fromStdString(parentMeshInsertion->Name) << "): (" << InsertionPointLocal.x << "," << InsertionPointLocal.y << "," << InsertionPointLocal.z << ")";
    qDebug() << QString::fromStdString(prefix) + "   Number of Nodes:" << MNodes.size();
    qDebug() << QString::fromStdString(prefix) + "   Current Muscle Length:" << MuscleLength << "m";
    qDebug() << QString::fromStdString(prefix) + "   Has Via Points:" << (bHasViaPoints ? "Yes" : "No");
    qDebug() << QString::fromStdString(prefix) + "   Meshes:";
    for (auto& m : meshPtrs) {
        qDebug() << QString::fromStdString(prefix) + "    - " + QString::fromStdString(m->Name);
    }
    
}

void SSMuscle::storeMNodesInitialGuess(int stepIdx, std::vector<MWMath::Point3D>& guessPath)
{
    for (size_t k = 0; k < MNodes.size(); ++k) {
        if (stepIdx >= MNodes[k].MNodeInitialGuessSteps.size()) {
            MNodes[k].MNodeInitialGuessSteps.resize(stepIdx + 1);
        }
        MNodes[k].MNodeInitialGuessSteps[stepIdx] = guessPath[k];
    }
}

void SSMuscle::storeMNodesInitialGuessColors(int stepIdx, std::vector<MWMath::Point3D> &colorPath)
{
    for (size_t k = 0; k < MNodes.size(); ++k) {
        if (stepIdx >= MNodes[k].MNodeInitialGuessColorSteps.size()) {
            MNodes[k].MNodeInitialGuessColorSteps.resize(stepIdx + 1);
        }
        MNodes[k].MNodeInitialGuessColorSteps[stepIdx] = colorPath[k];
    }
}

void SSMuscle::storeMNodesGlobalPositions()
{
    // and local
    //MNodes[0].MNodeGlobalSteps.push_back(OriginPointGlobal);
    MNodes[0].MNodeLocalSteps.push_back(OriginPointLocal);
    for (size_t k = 0; k < MNodes.size(); ++k) {
        MNodes[k].MNodeGlobalSteps.push_back(MNodes[k].PositionGlobal);
        MNodes[k].MNodeLocalSteps.push_back(MNodes[k].PositionLocal);

        if (MNodes[k].parentMesh != nullptr){
            if (MNodes[k].parentMesh->Parent == nullptr){
                MNodes[k].MNodeParentPosSteps.push_back(MWMath::Point3D(0,0,0));
                MNodes[k].MNodeParentRotSteps.push_back(MWMath::RotMatrix3x3());
            }
            else{
                MNodes[k].MNodeParentPosSteps.push_back(MNodes[k].parentMesh->Parent->PositionGlobal);
                MNodes[k].MNodeParentRotSteps.push_back(MNodes[k].parentMesh->Parent->OrientationGlobal);
            }
        }
        else{
            MNodes[k].MNodeParentPosSteps.push_back(MWMath::Point3D(0,0,0));
            MNodes[k].MNodeParentRotSteps.push_back(MWMath::RotMatrix3x3());
            qDebug() << "Warning: Muscle" << QString::fromStdString(Name) << "Node" << k << "has no parent mesh assigned!";
        }
        
    }
    //MNodes[MNodes.size()-1].MNodeGlobalSteps.push_back(InsertionPointGlobal);
    MNodes[MNodes.size()-1].MNodeLocalSteps.push_back(InsertionPointLocal);

}

void SSMuscle::initializeSimulationMuscle(int &steps)
{
    MuscleLengthSteps.reserve(steps);
    for (auto& node : MNodes) {
        node.MNodeEtaSteps.reserve(steps);
        node.MNodeGlobalSteps.reserve(steps);
        node.MNodeLocalSteps.reserve(steps);
        node.MNodeInitialGuessSteps.reserve(steps);
        node.MNodeInitialGuessColorSteps.reserve(steps);
        node.MNodeParentPosSteps.reserve(steps);
        node.MNodeParentRotSteps.reserve(steps);
    }

    for (auto mesh : meshPtrs) {
        if (mesh->bIsViaPoint) {
            bHasViaPoints = true;
            break;
        }
    }
}

std::vector<MWMath::Point3D> SSMuscle::gatherMNodesGlobal(int &step)
{
    std::vector<MWMath::Point3D> musclePoints;
    musclePoints.resize(MNodes.size());
    for (size_t k = 0; k < MNodes.size(); ++k) {
        if (step < MNodes[k].MNodeGlobalSteps.size()) {
            musclePoints[k] = MNodes[k].MNodeGlobalSteps[step];
        }
        else {
            qDebug() << "Fehler: Schritt" << step << "übersteigt gespeicherte globale Schritte für Node" << k;
        }
    }
    return musclePoints;
}

void SSMuscle::getAllMuscleMNodesStepValues(int totalSteps, std::vector<std::vector<MWMath::Point3D>> &globalPoints, std::vector<std::vector<MWMath::Point3D>> &localPoints, std::vector<std::vector<MWMath::Point3D>> &initialGuessPoints, std::vector<std::vector<MWMath::Point3D>> &initialGuessColors, std::vector<std::vector<std::vector<double>>> &etaValues)
{
    globalPoints.resize(totalSteps);
    localPoints.resize(totalSteps);
    initialGuessPoints.resize(totalSteps);
    initialGuessColors.resize(totalSteps);
    etaValues.resize(totalSteps);

    for (int step = 0; step < totalSteps; ++step) {
        globalPoints[step].resize(MNodes.size());
        localPoints[step].resize(MNodes.size());
        initialGuessPoints[step].resize(MNodes.size());
        initialGuessColors[step].resize(MNodes.size());
        etaValues[step].resize(MNodes.size());

        for (size_t k = 0; k < MNodes.size(); ++k) {
            if (step < MNodes[k].MNodeGlobalSteps.size()) {
                globalPoints[step][k] = MNodes[k].MNodeGlobalSteps[step];
            }
            if (step < MNodes[k].MNodeLocalSteps.size()) {
                localPoints[step][k] = MNodes[k].MNodeLocalSteps[step];
            }
            if (step < MNodes[k].MNodeInitialGuessSteps.size()) {
                initialGuessPoints[step][k] = MNodes[k].MNodeInitialGuessSteps[step];
            }
            if (step < MNodes[k].MNodeInitialGuessColorSteps.size()) {
                initialGuessColors[step][k] = MNodes[k].MNodeInitialGuessColorSteps[step];
            }
            if (step < MNodes[k].MNodeEtaSteps.size()) {
                etaValues[step][k] = MNodes[k].MNodeEtaSteps[step];
            }
        }
    }
}



// #########################################
void MuscleNode::getMNodeInfoStep(int idx)
{
    qDebug() << " Muscle Node << " << idx;
    qDebug() << "   PositionGlobal: (" << PositionGlobal.x << "," << PositionGlobal.y << "," << PositionGlobal.z << ")";
    qDebug() << "   PositionLocal: (" << PositionLocal.x << "," << PositionLocal.y << "," << PositionLocal.z << ")";
    if (parentMesh){
        qDebug() << "   ParentMesh: " << QString::fromStdString(parentMesh->Name);
        if(parentMesh->Parent)
        {
            qDebug() << "   Parent: " << QString::fromStdString(parentMesh->Parent->Name);
            qDebug() << "   Parent Pos" << parentMesh->Parent->PositionGlobal.x << "," << parentMesh->Parent->PositionGlobal.y << "," << parentMesh->Parent->PositionGlobal.z << "";
            /* qDebug() << "   Parent Rot Matrix: ";
            for (int r=0;r<3;r++){
                qDebug() << "       " << parentMesh->Parent->OrientationGlobal.m[r][0] << "," << parentMesh->Parent->OrientationGlobal.m[r][1] << "," << parentMesh->Parent->OrientationGlobal.m[r][2];
                qDebug() << "";
            } */
        }
        else{
            qDebug() << "   Parent Mesh -> parent: nullptr";
        }
        
    }
    else{
        qDebug() << "   Parent: nullptr";
    }
}
