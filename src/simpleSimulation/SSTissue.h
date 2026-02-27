#pragma once

#include "utils/MWMath.h"

#include <casadi/casadi.hpp>
#include <string>
#include <vector>
#include <memory>   

#include <QDebug>
#include <QString>

class SSMesh;

class SSTissue{
public:
    SSTissue() = default;
    virtual ~SSTissue() = default;

    std::string Name = "";
    MWMath::Point3D PositionGlobal = MWMath::Point3D(0.0, 0.0, 0.0);
    MWMath::RotMatrix3x3 OrientationGlobal = MWMath::RotMatrix3x3();
    MWMath::Point3D MeshColor = MWMath::Point3D(0.7, 0.7, 0.7); // RGB Werte zwischen 0 und 1

    std::vector<std::shared_ptr<SSMesh>> Meshes; // for joints and bodies
    std::vector<std::shared_ptr<SSTissue>> Children;
    MWMath::Point3D Position2ParentRelInParentFrame = MWMath::Point3D(0.0, 0.0, 0.0);
    MWMath::RotMatrix3x3 Orientation2ParentRel = MWMath::RotMatrix3x3();
    std::shared_ptr<SSTissue> Parent = nullptr;
    int SystemLayer = 0; // Root = 0, Kinder von Root = 1, etc. -> für Update-Reihenfolge und Debug-Ausgaben
    bool bDebug = false;
    //void setParent(std::shared_ptr<SSTissue> parent) { Parent = parent; }
    //std::shared_ptr<SSTissue> getParent() const { return Parent; }
    virtual int update(int step=0) {return 0;};

    // get the global position and orientation of the parent for casadi -> local parametrization
    virtual void getCasadiParentGPOs(casadi::MX& pos, casadi::MX& ori) {
        pos = casadi::MX::vertcat({PositionGlobal.x, PositionGlobal.y, PositionGlobal.z});
        ori = casadi::MX::zeros(3,3);
        for (int i=0;i<3;i++){for (int j=0;j<3;j++){ori(i,j) = OrientationGlobal.m[i][j];}}
    };

    virtual double getDistanceNumerically(MWMath::Point3D pGlobal) {
        double diff = distance(pGlobal, PositionGlobal);
        return diff;
    }

    void getInfo(){
        qDebug() << "Tissue: " << QString::fromStdString(Name);
        qDebug() << "   PositionGlobal: (" << PositionGlobal.x << "," << PositionGlobal.y << "," << PositionGlobal.z << ")";
        qDebug() << "   OrientationGlobal: ";
        for (int r=0;r<3;r++){
            qDebug() << "       " << OrientationGlobal.m[r][0] << "," << OrientationGlobal.m[r][1] << "," << OrientationGlobal.m[r][2];
            qDebug() << "";
        }
        qDebug() << "   Relative Position to Parent: (" << Position2ParentRelInParentFrame.x << "," << Position2ParentRelInParentFrame.y << "," << Position2ParentRelInParentFrame.z << ")";
        qDebug() << "   Relative Orientation to Parent: ";
        for (int r=0;r<3;r++){
            qDebug() << "       " << Orientation2ParentRel.m[r][0] << "," << Orientation2ParentRel.m[r][1] << "," << Orientation2ParentRel.m[r][2];
        }
        if (Parent)
            qDebug() << "   Parent: " << QString::fromStdString(Parent->Name);
        else
            qDebug() << "   Parent: nullptr";
    }
private:
    
};