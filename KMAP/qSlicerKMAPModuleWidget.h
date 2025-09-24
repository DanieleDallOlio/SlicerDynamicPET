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

#ifdef _WIN32
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#endif

// Slicer includes
#include "qSlicerAbstractModuleWidget.h"
#include "SegmentationChangeWatcher.h"
#include "KeyPressWatcher.h"

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
#include <QCheckBox>
#include <vtkSegmentation.h>
#include <vtkSlicerKMAPLogic.h>
#include <QProgressBar>
#include <QMessageBox>

#include <ctkDICOMDatabase.h>
#include <qSlicerApplication.h>
#include <dcmtk/dcmdata/dctk.h>
#include <PythonQt.h>
#include <qSlicerPythonManager.h>
#include <vtkMRMLPlotChartNode.h>
#include <vtkMRMLPlotSeriesNode.h>
#include <vtkMRMLLayoutNode.h>
#include <vtkMRMLPlotViewNode.h>
#include <vtkMRMLScalarVolumeDisplayNode.h>
#include <vtkMRMLProceduralColorNode.h>
#include <qMRMLPlotWidget.h>
#include <qSlicerLayoutManager.h>
#include <qMRMLPlotView.h>
#include <QKeyEvent>
#include <QMainWindow>
#include <vtkChartXY.h>
#include <vtkPlot.h>

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
  void clearFITMTGAdata();
  void enableTACbutton();
  void enableFITbutton();
  void enableFITMTGAbutton();
  void enableFITMTGAImgbutton();
  void enableFITTCMImgbutton();
  void runVuong(std::string sel1, std::string sel2, std::string segmentID);
  void runTCMstat(std::string sel1, std::string sel2, std::string segmentID);
  QVariantMap TACtoPythonDict();
  QVariantMap fittedTCMtoPythonDict();
  QVariantMap fittedMTGAtoPythonDict();
  QVariantMap TCMParamsToPythonDict();
  QVariantMap MTGAParamsToPythonDict();
  void RemoveExistingPlotChartAndTable();
  vtkMRMLPlotChartNode* GetOrCreatePlotChart();
  vtkMRMLTableNode* GetOrCreatePlotTable();
  bool checkdisplayedKMAP();

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
  void onExcelTCMfittedPathChanged(const QString& path);
  void onExcelMTGAfittedPathChanged(const QString& path);
  void onExcelTCMPathChanged(const QString& path);
  void onExcelMTGAPathChanged(const QString& path);
  void onSaveExcelbutton();
  void onSaveTCMExcelbutton();
  void onSaveMTGAExcelbutton();
  void onSaveTCMfittedExcelbutton();
  void onSaveMTGAfittedExcelbutton();
  void onPlotbutton();
  void onIFSelectionChanged(int index);
  void onIFMTGASelectionChanged(int index);
  void onIFImgSelectionChanged(int index);
  void onVOISelectAllbutton();
  void onVOIMTGASelectAllbutton();
  void onFITbutton();
  void onResetbutton();
  void onFITMTGAbutton();
  void onFITMTGAImgbutton();
  void onFITTCMImgbutton();
  void onVOISegmentsChanged();
  void onVOIMTGASegmentsChanged();
  void onModelsChanged();
  void onModelsMTGAChanged();
  void onModelsMTGAImgChanged();
  void onModelsTCMImgChanged();
  void onModelsAllbutton();
  void onModelsMTGAAllbutton();
  void onVOISelectionChanged(int index);
  void onVOIMTGASelectionChanged(int index);
  void onModelsTCMSelectAllbutton();
  void onModelsSelectAllMTGAImgbutton();
  void onModelsSelectAllTCMImgbutton();
  void onPlotTCMbutton();
  void onPlotMTGAbutton();
  void onOLSclicked();
  void onOLSImgclicked();
  void onWLSclicked();
  void onWLSImgclicked();
  void onRLSclicked();
  void onRLSImgclicked();
  void onSliderChanged(int index);
  void onSliderImgChanged(int index);
  void onStdFitclicked();
  void onStdFitImgclicked();
  void onWFitclicked();
  void onWFitImgclicked();
  void onMTGAModelBox(int index);
  void onTCMModelBox(int index);
  void onSelectedPoint(vtkStringArray* mrmlPlotSeriesIDs, vtkCollection* selectionCol);
  void onDeleteKeyPressed();

protected:
  QScopedPointer<qSlicerKMAPModuleWidgetPrivate> d_ptr;
  bool IsActive{false};
  vtkMRMLSequenceNode* sequencePETNode;
  vtkMRMLSequenceNode* segSequenceNode;
  vtkMRMLSequenceBrowserNode* sequenceBrowserPETNode;
  vtkSmartPointer<SegmentationChangeWatcher> SegWatcher;
private:
  Q_DECLARE_PRIVATE(qSlicerKMAPModuleWidget);
  Q_DISABLE_COPY(qSlicerKMAPModuleWidget);
  std::vector<double> extractColumn(const std::vector<std::vector<double>>& mat, const int index=0);
  vtkMRMLSubjectHierarchyNode* SubjectHierarchyNode;
  vtkIdType patID, stuID, ctID, petID, segID;
  vtkMRMLScalarVolumeNode* petNode;

  int PETdims[3];
  int numberOfTimepoints;
  std::vector<QString> segmentIDs;
  std::map<std::string, std::vector<VoxelStatistics>> segmentTACs;
  std::map<std::string, std::string> segmentTACsnames;
  std::map<std::string, std::map<std::string, MTGAParameters>> segmentMTGA;
  std::map<std::string, std::map<std::string, TCMParameters>> segmentTCM;
  std::map<std::string, std::vector<std::vector<double>>> segmentTAC4TCMfits;
  std::map<std::string, std::vector<bool>> segmentkeep4TCMfits;
  std::map<std::string, std::map<std::string, double*>> segmentTCMfits;
  std::map<std::string, std::vector<MTGAParameters>> MTGAImgOutcomes;
  std::map<std::string, std::vector<TCMParameters>> TCMImgOutcomes;

  QProgressBar* ProgressBar;
  QPushButton* stopButton;
  std::atomic<bool> stopRequested;
  std::vector<double> timePoints, durations;
  QStringList checkboxNames, ModelsNamesTCM, ModelsNamesMTGA, StatsNames;
  std :: string IFID, plotTCMVOI, plotMTGAVOI, plotMTGAModel;
  std :: vector < std :: string > VOIsegmentIDs, VOIMTGAsegmentIDs, modelsID, modelsMTGAID, modelsMTGAImgID, modelsTCMImgID;
  int PlotSelectedFrame;
  std :: string PlotSelectedVOI;
  KeyPressWatcher* keyWatcher{nullptr};
  QSet<QPair<QString, vtkIdType>> lastSelection;
  std::unordered_map<std::string, std::string> ColNameToSegmentID;
  QMap<QString, vtkPlot*> MapPlotSeriesNodeIDToPlot;
  std::vector<std::vector<double>> PET_flatten_values;
};

#endif
