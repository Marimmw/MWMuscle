
#include <QDialog>
#include <QVBoxLayout>
#include <QSlider>
#include <QDebug>
#include <vector>

#include <memory>

// VTK Includes
//#include <vtkGlyph3D.h>
#include <QVTKOpenGLNativeWidget.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkSphereSource.h>
#include <vtkCylinderSource.h>       
#include <vtkTransform.h>            
#include <vtkTransformPolyDataFilter.h>
#include <vtkParametricTorus.h>        
#include <vtkParametricFunctionSource.h>
#include <vtkProperty.h>
#include <vtkMatrix4x4.h>
#include <vtkPointData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkAxesActor.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>

#include <QMessageBox>
#include <QCloseEvent>
#include <QDebug>

#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <QDir>
#include <QString>

#include "utils/MWMath.h"
#include "SSMeshes.h" 
#include "SSMuscle.h"

class VTKSimViewerSimple : public QDialog {
    Q_OBJECT

public:
    VTKSimViewerSimple(const std::vector<std::vector<std::vector<MWMath::Point3D>>>& muscleResults, 
                       const std::vector<std::vector<MWMath::RotMatrix3x3>>& meshResults,
                       const std::vector<std::vector<std::vector<MWMath::Point3D>>>& initialGuesses,
                       const std::vector<std::vector<std::vector<MWMath::Point3D>>>& initialGuessColors,
                       std::vector<std::shared_ptr<SSMesh>> meshes,
                    const std::vector<SSMuscle*>& muscles, 
                       const std::vector<double>& angles,
                       const std::vector<std::vector<MWMath::Point3D>>& otherPointsWithColors,
                       const std::vector<std::vector<std::string>>& solverInfoText,
                       QWidget* parent = nullptr);
    void exportAllStepsToImages(std::string outputFolder);
    bool ShowCoordinateSystems = true;

protected:
    void closeEvent(QCloseEvent *event) override;
private slots:
    void updateStep(int step);
    
private:
    void setupVtkPipeline();
    void updateMeshTransform(vtkProp3D* actor, const MWMath::RotMatrix3x3& R, const MWMath::Point3D& pos);
    

    // Daten-Referenzen (für mehrere Muskeln)
    const std::vector<std::vector<std::vector<MWMath::Point3D>>>& m_muscleResults;      // [muskel][timestep][punkt]
    const std::vector<std::vector<std::vector<MWMath::Point3D>>>& m_initialGuesses;     // [muskel][timestep][punkt]
    const std::vector<std::vector<std::vector<MWMath::Point3D>>>& m_initialGuessColors; // [muskel][timestep][punkt]
    const std::vector<std::vector<MWMath::RotMatrix3x3>>& m_meshResults;                // [mesh][timestep]
    std::vector<std::shared_ptr<SSMesh>> m_meshes;
    std::vector<SSMuscle*> m_muscles;
    std::vector<double> m_angles;
    const std::vector<std::vector<MWMath::Point3D>>& m_otherPointsWithColors;    
    //std::vector<std::string> m_solverInfoText;       // [timestep][p1, c1, p2, c2, ...]
    const std::vector<std::vector<std::string>>& m_solverInfoText;       // [timestep][p1, c1, p2, c2, ...]
    std::vector<vtkSmartPointer<vtkAxesActor>> m_meshAxesActors; // Pro Mesh ein Achsenkreuz
    vtkSmartPointer<vtkAxesActor> m_worldAxes; // Ein Achsenkreuz für den Welt-Ursprung

    // UI & VTK
    QVTKOpenGLNativeWidget* m_vtkWidget;
    QSlider* m_slider;
    vtkSmartPointer<vtkRenderer> m_renderer;
    
    // VTK Actors und Daten (für mehrere Muskeln)
    std::vector<vtkSmartPointer<vtkActor>> m_meshActors;                               // [mesh]
    std::vector<vtkSmartPointer<vtkActor>> m_muscleActors;                             // [muskel]
    std::vector<vtkSmartPointer<vtkPolyData>> m_muscleLines;                           // [muskel]
    std::vector<vtkSmartPointer<vtkPoints>> m_musclePointsData;                        // [muskel]
    
    // Initial Guess Visualisierung (für mehrere Muskeln)
    std::vector<vtkSmartPointer<vtkPolyData>> m_initialLines;                          // [muskel] - Linien
    std::vector<vtkSmartPointer<vtkPoints>> m_initialPointsData;                       // [muskel] - Punkt-Daten für Linien
    std::vector<std::vector<vtkSmartPointer<vtkActor>>> m_initialPointActors;          // [muskel][punkt] - Kugel-Actors
    std::vector<std::vector<vtkSmartPointer<vtkActor>>> m_otherPointActors; 
    
    vtkSmartPointer<vtkTextActor> m_textActor;
    vtkSmartPointer<vtkTextActor> m_stepTextActor;// [timestep][punkt] - Kugel-Actors

};