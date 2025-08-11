/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

#ifndef __qSlicerKMAPModuleWidget_h
#define __qSlicerKMAPModuleWidget_h

// Slicer includes
#include "qSlicerAbstractModuleWidget.h"

#include "qSlicerKMAPModuleExport.h"
#include "qMRMLNodeComboBox.h"
#include <vtkCollection.h>
#include <vtkSmartPointer.h>
#include <vtkMRMLScene.h>
#include <vtkStringArray.h>
#include <vtkMRMLScalarVolumeNode.h>
#include <vtkMRMLSegmentationNode.h>
#include <vtkMRMLSubjectHierarchyNode.h>
#include <vector>
#include <QPointer>
#include <QCheckBox.h>
#include <vtkSegmentation.h>
#include <vtkSlicerKMAPLogic.h>
#include <QProgressBar.h>

#include <ctkDICOMDatabase.h>
#include <qSlicerApplication.h>
#include <dcmtk/dcmdata/dctk.h>
#include <PythonQt.h>
#include <qSlicerPythonManager.h>
#include <vtkMRMLPlotChartNode.h>
#include <vtkMRMLPlotSeriesNode.h>
#include <vtkMRMLLayoutNode.h>
#include <vtkMRMLPlotViewNode.h>

class qSlicerKMAPModuleWidgetPrivate;
class vtkMRMLNode;

class Q_SLICER_QTMODULES_KMAP_EXPORT qSlicerKMAPModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:

  typedef qSlicerAbstractModuleWidget Superclass;
  explicit qSlicerKMAPModuleWidget(QWidget* parent = nullptr);

  ~qSlicerKMAPModuleWidget() override;

  // void setNodeSelectorEnabled(qMRMLNodeComboBox* selector, bool enabled);
  // static std::map<std::string, vtkIdType> GetStudyAndPatientAncestors(vtkMRMLSubjectHierarchyNode* shNode, vtkIdType volumeNodeId);

  void enter() override;
  void exit() override;
  void setMRMLScene(vtkMRMLScene* scene) override;
  void getDurations();
  void clearTACdata();
  void clearFITdata();
  void enableTACbutton();
  void enableFITbutton();
  QVariantMap TACtoPythonDict();
  void RemoveExistingPlotChartAndTable();
  vtkMRMLPlotChartNode* GetOrCreatePlotChart();
  vtkMRMLTableNode* GetOrCreatePlotTable();

public slots:
  void onSubjectHierarchyChanged();
  void onPatChanged(int index);
  void onStuChanged(int index);
  void onCTChanged(int index);
  void onPETChanged(int index);
  void onSegChanged(int index);
  void onSegmentsChanged();
  void onTACbutton();
  void onSelectAllbutton();
  void onExcelPathChanged(const QString& path);
  void onSaveExcelbutton();
  void onPlotbutton();
  void onIFSelectionChanged(int index);
  void onVOISelectAllbutton();
  void onFITbutton();
  void onVOISegmentsChanged();
  void onModelsChanged();
  void onModelsAllbutton();
  void onVOISelectionChanged(int index);
  void onModelsTCMSelectAllbutton();
  void onPlotTCMbutton();

protected:
  QScopedPointer<qSlicerKMAPModuleWidgetPrivate> d_ptr;
  bool IsActive{false};

private:
  Q_DECLARE_PRIVATE(qSlicerKMAPModuleWidget);
  Q_DISABLE_COPY(qSlicerKMAPModuleWidget);
  vtkMRMLSubjectHierarchyNode* SubjectHierarchyNode;
  vtkIdType patID, stuID, ctID, petID, segID;
  std::vector<QString> segmentIDs;
  std::map<std::string, std::vector<VoxelStatistics>> segmentTACs;
  std::map<std::string, std::string> segmentTACsnames;
  std::map<std::string, std::map<std::string, TCMParameters>> segmentMTGA;
  std::map<std::string, std::map<std::string, TCMParameters>> segmentTCM;
  std::map<std::string, std::vector<std::vector<double>>> segmentTAC4TCMfits;
  std::map<std::string, std::map<std::string, double*>> segmentTCMfits;
  QProgressBar* ProgressBar;
  std::vector<double> timePoints, durations;
  QStringList checkboxNames, ModelsNames, StatsNames;
  std :: string IFID, plotTCMVOI;
  std :: vector < std :: string > VOIsegmentIDs, modelsID;
};

#endif
