/* #pragma once
#include <QDialog>
#include <QSlider>
#include <QVBoxLayout>
#include <vector>
#include <memory>
#include <QVTKOpenGLNativeWidget.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <casadi/casadi.hpp>
#include "MWMath.h"
#include "../model/MWTissueObjects/MWMuscle.h"
#include "../model/MWMeshes/MWEllipsoidMesh.h"

class RawSimulation {
public:
    RawSimulation() = default;
    void run();
    std::vector<std::vector<std::vector<MWMath::Point3D>>> allMusclesResults;
    std::vector<std::vector<MWMath::RotMatrix3x3>> allMeshesRMatrixResults;

    void runMovementSequence(std::shared_ptr<MWMuscle> muscle, std::vector<std::shared_ptr<MWMesh>> meshes, const std::vector<double>& zRotationSteps);
    std::vector<double> solveSingleMuscleGlobal(
        std::shared_ptr<MWMuscle> muscle, 
        const std::vector<std::shared_ptr<MWMesh>>& activeMeshes,
        const std::vector<double>& warmstart // <--- DIESER PARAMETER MUSS HIER REIN
    );
};

class VTKSimViewer : public QDialog {
    Q_OBJECT
public:
    VTKSimViewer(const std::vector<std::vector<std::vector<MWMath::Point3D>>>& muscleResults, 
                 const std::vector<std::vector<MWMath::RotMatrix3x3>>& meshResults,
                 std::vector<std::shared_ptr<MWMesh>> meshes, 
                 const std::vector<double>& angles);
private slots:
    void updateStep(int step);
private:
    void setupVtkPipeline();
    void updateMeshTransform(vtkActor* actor, const MWMath::RotMatrix3x3& R);

    const std::vector<std::vector<std::vector<MWMath::Point3D>>>& m_muscleResults; // Hier umbenannt
    const std::vector<std::vector<MWMath::RotMatrix3x3>>& m_meshResults;
    std::vector<std::shared_ptr<MWMesh>> m_meshes;
    std::vector<double> m_angles;

    QVTKOpenGLNativeWidget* m_vtkWidget;
    QSlider* m_slider;
    vtkSmartPointer<vtkRenderer> m_renderer;
    std::vector<vtkSmartPointer<vtkActor>> m_meshActors;
    std::vector<vtkSmartPointer<vtkActor>> m_muscleActors;
    std::vector<vtkSmartPointer<vtkPolyData>> m_muscleLines;
    std::vector<vtkSmartPointer<vtkPoints>> m_musclePointsData;
}; */