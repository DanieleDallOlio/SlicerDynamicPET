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

#ifdef _WIN32
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#endif


// Qt includes
#include <QDebug>

// Slicer includes
#include "qSlicerKMAPModuleWidget.h"
#include "ui_qSlicerKMAPModuleWidget.h"
#include <QApplication>

#ifdef _WIN32
#include <sstream>
#include <iomanip>

// strptime replacement for Windows
inline char* strptime(const char* s, const char* f, struct tm* tm) {
    std::istringstream input(s);
    input >> std::get_time(tm, f);
    if (input.fail()) return nullptr;
    return (char*)(s + input.tellg());
}
#endif



//-----------------------------------------------------------------------------
class qSlicerKMAPModuleWidgetPrivate: public Ui_qSlicerKMAPModuleWidget
{
  Q_DECLARE_PUBLIC(qSlicerKMAPModuleWidget);

protected:
  qSlicerKMAPModuleWidget* const q_ptr;
  void setDoubleField(QLineEdit* le, double lo, double hi, int decimals);
  void setIntField(QLineEdit* le, int lo, int hi);
public:
  qSlicerKMAPModuleWidgetPrivate(qSlicerKMAPModuleWidget& object);
  ~qSlicerKMAPModuleWidgetPrivate()=default;
  void init();
  void populatePatientComboBox();
  void populateStudyComboBox(vtkIdType patientID);
  void populateNodeComboBox(QComboBox* comboBox, vtkIdType parentItemID, const char * requiredNodeType, const std :: string requiredModality);
  void populateSegmentCheckboxes(vtkIdType SegItemID);
  void populatePlotSegmentCheckboxes();
  void populateIF();
  void populateIFMTGA();
  void populateIFImg();
  void populateVOI(std :: string ifID);
  void populateVOIMTGA(std :: string ifID);
  void populateResultsVOI();
  void populateResultsVOIMTGA();
  void populateResultsTable(std :: string segmentID);
  void populateResultsMTGATable(std :: string segmentID);
  void populateModelsTCM(std :: string segmentID);
  void populateModelsMTGA(std :: string segmentID);
  void populateTimeBarMTGA();
  void populateTimeBarMTGAImg();
  void populateModelCombo(QComboBox* comboToFill,
                          const std::string& otherSelectedModel,
                          const std::string& currentSelectedModel,
                          const std::string& segmentID);
  void populateModelComboTCM(QComboBox* comboToFill,
                             const std::string& otherSelectedModel,
                             const std::string& currentSelectedModel,
                             const std::string& segmentID);
};

//-----------------------------------------------------------------------------
// qSlicerKMAPModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerKMAPModuleWidgetPrivate::qSlicerKMAPModuleWidgetPrivate(qSlicerKMAPModuleWidget& object): q_ptr(&object)
{
  Q_Q(qSlicerKMAPModuleWidget);
}


void qSlicerKMAPModuleWidgetPrivate::setDoubleField(QLineEdit* le, double lo, double hi, int decimals)
{

  QObject::connect(le, &QLineEdit::editingFinished, le, [le, lo, hi, decimals]() {
    const QLocale loc = le->locale(); // respect UI locale (comma/dot)
    bool ok = false;
    double x = loc.toDouble(le->text(), &ok);
    if (!ok) {
      le->setText(loc.toString(lo, 'f', decimals));
      return;
    }
    if (x < lo) x = lo;
    if (x > hi) x = hi;
    le->setText(loc.toString(x, 'f', decimals));
  });

  return ;
}

void qSlicerKMAPModuleWidgetPrivate::setIntField(QLineEdit* le, int lo, int hi)
{
  QObject::connect(le, &QLineEdit::editingFinished, le, [le, lo, hi]() {
    bool ok = false;
    int val = le->text().toInt(&ok);
    if (!ok) {
      le->setText(QString::number(lo));
      return;
    }
    if (val < lo) val = lo;
    if (val > hi) val = hi;
    le->setText(QString::number(val));
  });

  return ;
}

//-----------------------------------------------------------------------------
void qSlicerKMAPModuleWidgetPrivate::init()
{
  Q_Q(qSlicerKMAPModuleWidget);
  this->setupUi(q);


  this->StuSelector->setEnabled(false);
  this->CTSelector->setEnabled(false);
  this->PETSelector->setEnabled(false);
  this->SegSelector->setEnabled(false);
  this->segmentSelectAll->setEnabled(false);
  this->saveExcelButton->setEnabled(false);
  this->PlotsTabWidget->setTabEnabled(this->PlotsTabWidget->indexOf(this->MTGAWidget), false);
  this->MTGAWidget->setEnabled(false);
  this->PlotsTabWidget->setTabEnabled(this->PlotsTabWidget->indexOf(this->TCMWidget), false);
  this->TCMWidget->setEnabled(false);
  this->PlotsTabWidget->setTabEnabled(this->PlotsTabWidget->indexOf(this->ImagingWidget), false);
  this->ImagingWidget->setEnabled(false);

  this->VOIsegmentSelectAll->setEnabled(false);
  this->VOIMTGAsegmentSelectAll->setEnabled(false);
  this->saveTCMfittedExcelButton->setEnabled(false);
  this->saveMTGAfittedExcelButton->setEnabled(false);
  this->saveTCMExcelButton->setEnabled(false);
  this->saveMTGAExcelButton->setEnabled(false);

  this->TCMResultsTable->setSortingEnabled(true);
  this->MTGAResultsTable->setSortingEnabled(true);

  // Make connections
  QObject::connect( this->PatSelector, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onPatChanged(int)));
  QObject::connect( this->StuSelector, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onStuChanged(int)));
  QObject::connect( this->CTSelector, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onCTChanged(int)) );
  QObject::connect( this->PETSelector, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onPETChanged(int)) );
  QObject::connect( this->SegSelector, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onSegChanged(int)));
  QObject::connect( this->TACbutton, SIGNAL(clicked(bool)),
    q, SLOT(onTACbutton()));
  QObject::connect( this->segmentSelectAll, SIGNAL(clicked(bool)),
    q, SLOT(onSelectAllbutton()));
  QObject::connect( this->direxcel, SIGNAL(currentPathChanged(const QString&)),
    q, SLOT(onExcelPathChanged(const QString&)));
  QObject::connect( this->fileexcel, SIGNAL(textChanged(const QString&)),
    q, SLOT(onExcelPathChanged(const QString&)));
  QObject::connect( this->saveExcelButton, SIGNAL(clicked(bool)),
    q, SLOT(onSaveExcelbutton()));
  QObject::connect( this->plotButton, SIGNAL(clicked(bool)),
    q, SLOT(onPlotbutton()));
  QObject::connect( this->plotTCMButton, SIGNAL(clicked(bool)),
    q, SLOT(onPlotTCMbutton()));
  QObject::connect( this->plotMTGAButton, SIGNAL(clicked(bool)),
    q, SLOT(onPlotMTGAbutton()));
  // QObject::connect( this->PlotErrorCheckbox, SIGNAL(toggled(bool)),
  //   q, SLOT(onPlotbutton()));
  QObject::connect( this->IFSelector, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onIFSelectionChanged(int)));
  QObject::connect( this->IFSelectorMTGA, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onIFMTGASelectionChanged(int)));
  QObject::connect( this->IFSelectorImg, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onIFImgSelectionChanged(int)));
  QObject::connect( this->VOIsegmentSelectAll, SIGNAL(clicked(bool)),
    q, SLOT(onVOISelectAllbutton()));
  QObject::connect( this->VOIMTGAsegmentSelectAll, SIGNAL(clicked(bool)),
    q, SLOT(onVOIMTGASelectAllbutton()));
  QObject::connect( this->ModelsSelectAll, SIGNAL(clicked(bool)),
    q, SLOT(onModelsAllbutton()));
  QObject::connect( this->ModelsMTGASelectAll, SIGNAL(clicked(bool)),
      q, SLOT(onModelsMTGAAllbutton()));
  QObject::connect( this->ModelsTCMSelectAll, SIGNAL(clicked(bool)),
      q, SLOT(onModelsTCMSelectAllbutton()));
  QObject::connect( this->ModelsSelectAllMTGAImg, SIGNAL(clicked(bool)),
      q, SLOT(onModelsSelectAllMTGAImgbutton()));
  QObject::connect( this->ModelsSelectAllTCMImg, SIGNAL(clicked(bool)),
      q, SLOT(onModelsSelectAllTCMImgbutton()));
  QObject::connect( this->FITbutton, SIGNAL(clicked(bool)),
    q, SLOT(onFITbutton()));
  QObject::connect( this->FITbuttonTCMImg, SIGNAL(clicked(bool)),
    q, SLOT(onFITTCMImgbutton()));
  QObject::connect( this->RESETbutton, SIGNAL(clicked(bool)),
    q, SLOT(onResetbutton()));
  QObject::connect( this->FITMTGAbutton, SIGNAL(clicked(bool)),
    q, SLOT(onFITMTGAbutton()));
  QObject::connect( this->FITbuttonMTGAImg, SIGNAL(clicked(bool)),
    q, SLOT(onFITMTGAImgbutton()));
  QObject::connect( this->VOISelector, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onVOISelectionChanged(int)));
  QObject::connect( this->VOISelectorMTGA, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onVOIMTGASelectionChanged(int)));
  QObject::connect( this->direxceltcm, SIGNAL(currentPathChanged(const QString&)),
    q, SLOT(onExcelTCMPathChanged(const QString&)));
  QObject::connect( this->direxcelmtga, SIGNAL(currentPathChanged(const QString&)),
    q, SLOT(onExcelMTGAPathChanged(const QString&)));
  QObject::connect( this->fileexceltcm, SIGNAL(textChanged(const QString&)),
    q, SLOT(onExcelTCMPathChanged(const QString&)));
  QObject::connect( this->fileexcelmtga, SIGNAL(textChanged(const QString&)),
    q, SLOT(onExcelMTGAPathChanged(const QString&)));
  QObject::connect( this->direxceltcmfitted, SIGNAL(currentPathChanged(const QString&)),
    q, SLOT(onExcelTCMfittedPathChanged(const QString&)));
  QObject::connect( this->direxcelmtgafitted, SIGNAL(currentPathChanged(const QString&)),
    q, SLOT(onExcelMTGAfittedPathChanged(const QString&)));
  QObject::connect( this->fileexceltcmfitted, SIGNAL(textChanged(const QString&)),
    q, SLOT(onExcelTCMfittedPathChanged(const QString&)));
  QObject::connect( this->fileexcelmtgafitted, SIGNAL(textChanged(const QString&)),
    q, SLOT(onExcelMTGAfittedPathChanged(const QString&)));
  QObject::connect( this->saveTCMExcelButton, SIGNAL(clicked(bool)),
    q, SLOT(onSaveTCMExcelbutton()));
  QObject::connect( this->saveMTGAExcelButton, SIGNAL(clicked(bool)),
    q, SLOT(onSaveMTGAExcelbutton()));
  QObject::connect( this->saveTCMfittedExcelButton, SIGNAL(clicked(bool)),
    q, SLOT(onSaveTCMfittedExcelbutton()));
  QObject::connect( this->saveMTGAfittedExcelButton, SIGNAL(clicked(bool)),
    q, SLOT(onSaveMTGAfittedExcelbutton()));
  QObject::connect(this->olsFitCheckBox, SIGNAL(toggled(bool)),
    q, SLOT(onOLSclicked()));
  QObject::connect(this->olsFitCheckBoxImg, SIGNAL(toggled(bool)),
    q, SLOT(onOLSImgclicked()));
  QObject::connect(this->weightedFitCheckBox, SIGNAL(toggled(bool)),
    q, SLOT(onWLSclicked()));
  QObject::connect(this->weightedFitCheckBoxImg, SIGNAL(toggled(bool)),
    q, SLOT(onWLSImgclicked()));
  QObject::connect(this->robustFitCheckBox, SIGNAL(toggled(bool)),
    q, SLOT(onRLSclicked()));
  QObject::connect(this->robustFitCheckBoxImg, SIGNAL(toggled(bool)),
    q, SLOT(onRLSImgclicked()));
  QObject::connect(this->standardFitCheckBox, SIGNAL(toggled(bool)),
    q, SLOT(onStdFitclicked()));
  QObject::connect(this->standardFitCheckBoxImg, SIGNAL(toggled(bool)),
    q, SLOT(onStdFitImgclicked()));
  QObject::connect(this->weightFitCheckBox, SIGNAL(toggled(bool)),
    q, SLOT(onWFitclicked()));
  QObject::connect(this->weightFitCheckBoxImg, SIGNAL(toggled(bool)),
    q, SLOT(onWFitImgclicked()));
  QObject::connect(this->MTGAModel1, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onMTGAModelBox(int)));
  QObject::connect(this->MTGAModel2, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onMTGAModelBox(int)));
  QObject::connect(this->TCMModel1, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onTCMModelBox(int)));
  QObject::connect(this->TCMModel2, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onTCMModelBox(int)));


  // MTGA controls
  this->setDoubleField(this->framingNormEdit, 0.0, 3600.0, 2);
  this->setDoubleField(this->huberTuneEdit,  1e-3, 10.0,   6);
  this->setDoubleField(this->tolEdit,        1e-12, 1e-1, 12);
  this->setIntField   (this->maxIterEdit,    1,     100000);

  // MTGA imaging controls
  this->setDoubleField(this->framingNormEditImg, 0.0, 3600.0, 2);
  this->setDoubleField(this->huberTuneEditImg,  1e-3, 10.0,   6);
  this->setDoubleField(this->tolEditImg,        1e-12, 1e-1, 12);
  this->setIntField   (this->maxIterEditImg,    1,     100000);
  #ifdef HAVE_OPENMP
  this->setIntField(this->numThreadsMTGA, 1, omp_get_max_threads());
  this->numThreadsMTGA->setText(QString::number(omp_get_max_threads()));
  #else
  this->setIntField(this->numThreadsMTGA, 1, 1);
  #endif

  // TCM params
  this->setDoubleField(this->k1Initial, 0.0, 10.0, 6);
  this->setDoubleField(this->k1Lower,   0.0, 10.0, 6);
  this->setDoubleField(this->k1Upper,   0.0, 10.0, 6);

  this->setDoubleField(this->k2Initial, 0.0, 10.0, 6);
  this->setDoubleField(this->k2Lower,   0.0, 10.0, 6);
  this->setDoubleField(this->k2Upper,   0.0, 10.0, 6);

  this->setDoubleField(this->k3Initial, 0.0, 10.0, 6);
  this->setDoubleField(this->k3Lower,   0.0, 10.0, 6);
  this->setDoubleField(this->k3Upper,   0.0, 10.0, 6);

  this->setDoubleField(this->k4Initial, 0.0, 10.0, 6);
  this->setDoubleField(this->k4Lower,   0.0, 10.0, 6);
  this->setDoubleField(this->k4Upper,   0.0, 10.0, 6);

  this->setDoubleField(this->vbInitial, 0.0, 1.0, 6);
  this->setDoubleField(this->vbLower,   0.0, 1.0, 6);
  this->setDoubleField(this->vbUpper,   0.0, 1.0, 6);

  this->setDoubleField(this->tdInitial, -10.0, 600., 3);
  this->setDoubleField(this->tdLower,   -10.0, 600., 3);
  this->setDoubleField(this->tdUpper,   -10.0, 600., 3);

  this->setDoubleField(this->decayConstEdit, 1e-6, 10.0, 10);
  this->setDoubleField(this->timeStepEdit,   0.001, 60.0, 6);

  this->setDoubleField(this->pbrp1Edit, 0.0, 1.0, 6);
  this->setDoubleField(this->pbrp2Edit, 0.0, 1.0, 6);
  this->setDoubleField(this->pbrp3Edit, 0.0, 1.0, 6);

  this->setIntField(this->maxIterTCMEdit, 1, 100000);

  // TCM imaging params
  this->setDoubleField(this->k1InitialImg, 0.0, 10.0, 6);
  this->setDoubleField(this->k1LowerImg,   0.0, 10.0, 6);
  this->setDoubleField(this->k1UpperImg,   0.0, 10.0, 6);

  this->setDoubleField(this->k2InitialImg, 0.0, 10.0, 6);
  this->setDoubleField(this->k2LowerImg,   0.0, 10.0, 6);
  this->setDoubleField(this->k2UpperImg,   0.0, 10.0, 6);

  this->setDoubleField(this->k3InitialImg, 0.0, 10.0, 6);
  this->setDoubleField(this->k3LowerImg,   0.0, 10.0, 6);
  this->setDoubleField(this->k3UpperImg,   0.0, 10.0, 6);

  this->setDoubleField(this->k4InitialImg, 0.0, 10.0, 6);
  this->setDoubleField(this->k4LowerImg,   0.0, 10.0, 6);
  this->setDoubleField(this->k4UpperImg,   0.0, 10.0, 6);

  this->setDoubleField(this->vbInitialImg, 0.0, 1.0, 6);
  this->setDoubleField(this->vbLowerImg,   0.0, 1.0, 6);
  this->setDoubleField(this->vbUpperImg,   0.0, 1.0, 6);

  this->setDoubleField(this->tdInitialImg, -10.0, 600.0, 3);
  this->setDoubleField(this->tdLowerImg,   -10.0, 600.0, 3);
  this->setDoubleField(this->tdUpperImg,   -10.0, 600.0, 3);

  this->setDoubleField(this->decayConstEditImg, 1e-6, 10.0, 10);
  this->setDoubleField(this->timeStepEditImg,   0.001, 60.0, 6);

  this->setDoubleField(this->pbrp1EditImg, 0.0, 1.0, 6);
  this->setDoubleField(this->pbrp2EditImg, 0.0, 1.0, 6);
  this->setDoubleField(this->pbrp3EditImg, 0.0, 1.0, 6);

  this->setIntField(this->maxIterTCMEditImg, 1, 100000);
  #ifdef HAVE_OPENMP
  this->setIntField(this->numThreadsTCM, 1, omp_get_max_threads());
  this->numThreadsTCM->setText(QString::number(omp_get_max_threads()));
  #else
  this->setIntField(this->numThreadsTCM, 1, 1);
  #endif


  this->TACCollapsibleButton->setCollapsed(true);
  this->TCMCollapsibleButton->setCollapsed(true);
  this->MTGACollapsibleButton->setCollapsed(true);
  this->MTGAStatTestButton->setCollapsed(true);
  this->MTGAStatTestButton->setEnabled(false);
  this->TCMStatTestButton->setCollapsed(true);
  this->TCMStatTestButton->setEnabled(false);
  for (const QString& name : q->checkboxNames)
  {
    QCheckBox* cb = new QCheckBox(name, this->PlotStatsCheckContents);
    this->PlotStatsCheckLayout->addWidget(cb);
  }

  PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
  mainContext.evalScript(R"PYTHON(
try:
    import pandas as pd
except ImportError:
    import slicer
    slicer.util.pip_install("pandas")
    import pandas as pd

try:
    import xlsxwriter
except ImportError:
    import slicer
    slicer.util.pip_install("xlsxwriter")

try:
    import numpy as np
except ImportError:
    import slicer
    slicer.util.pip_install("numpy")

def DPE_save_multisheet_excel(filepath, sheet_data_dict):
    """
    filepath: str - full path to xlsx
    sheet_data_dict: dict[str, list[list[str]]] - sheet name to 2D table
    """
    with pd.ExcelWriter(filepath, engine="xlsxwriter") as writer:
        for sheet, data in sheet_data_dict.items():
            df = pd.DataFrame(data)[["Time(s)", "Duration", "Mean", "Median", "StDev","IQR","Min", "Max", "Q1", "Q3", "VoxelCount","Volume(mm3)","Volume(cm3)"]]
            df.to_excel(writer, sheet_name=sheet, index=False)

def DPE_saveTCM_multisheet_excel(filepath, sheet_data_dict):
    """
    filepath: str - full path to xlsx
    sheet_data_dict: dict[str, list[list[str]]] - sheet name to 2D table
    """
    with pd.ExcelWriter(filepath, engine="xlsxwriter") as writer:
        for sheet, data in sheet_data_dict.items():
            df = pd.DataFrame(data)[["Model", "K1", "k2", "k3", "k4", "vb", "td", "Ki", "DV", "AIC", "BIC", "MASE", "chi^2_nu"]]
            df.to_excel(writer, sheet_name=sheet, index=False)

def DPE_saveMTGA_multisheet_excel(filepath, sheet_data_dict):
    """
    filepath: str - full path to xlsx
    sheet_data_dict: dict[str, list[list[str]]] - sheet name to 2D table
    """
    with pd.ExcelWriter(filepath, engine="xlsxwriter") as writer:
        for sheet, data in sheet_data_dict.items():
            df = pd.DataFrame(data)[["Model", "Ki", "DV", "Intercept", "R2", "AIC", "MASE"]]
            df.to_excel(writer, sheet_name=sheet, index=False)

def DPE_generic_save_multisheet_excel(filepath, sheet_data_dict):
    """
    filepath: str - full path to xlsx
    sheet_data_dict: dict[str, list[dict]] - sheet name to list of row dicts
    """
    with pd.ExcelWriter(filepath, engine="xlsxwriter") as writer:
        for sheet, data in sheet_data_dict.items():
            if not data:
                continue  # skip empty sheets

            # Create DataFrame from list of dicts — columns inferred automatically
            df = pd.DataFrame(data)

            # Optional: ensure "Time(s)" is first column if present
            if "Time(s)" in df.columns:
                taccols = [x for x in df.columns if "TAC" in x]
                if len(taccols)>0:
                  cols = ["Time(s)"] + taccols + [c for c in df.columns if not np.isin(c, ["Time(s)"]+taccols)]
                else:
                  cols = ["Time(s)"] + [c for c in df.columns if c != "Time(s)"]
                df = df[cols]

            df.to_excel(writer, sheet_name=sheet, index=False)

def DPE_genericMTGA_save_multisheet_excel(filepath, sheet_data_dict):
    """
    filepath: str - full path to xlsx
    sheet_data_dict: dict[str, list[dict]] - sheet name to list of row dicts
    """
    with pd.ExcelWriter(filepath, engine="xlsxwriter") as writer:
        for sheet, data in sheet_data_dict.items():
            if not data:
                continue  # skip empty sheets

            # Create DataFrame from list of dicts — columns inferred automatically
            df = pd.DataFrame(data)

            # Optional: ensure "Time(s)" is first column if present
            cols = []
            if "Patlak_x" in df.columns:
              cols += ["Patlak_x", "Patlak_y", "Patlak_fitted"]
            if "Logan_x" in df.columns:
              cols += ["Logan_x", "Logan_y", "Logan_fitted"]
            if "RE_x" in df.columns:
              cols += ["RE_x", "RE_y", "RE_fitted"]
            if len(cols)>0:
              df = df[cols]

            df.to_excel(writer, sheet_name=sheet, index=False)
)PYTHON");

  for (const QString& name : q->ModelsNamesMTGA)
  {
    QCheckBox* cb = new QCheckBox(name, this->ModelsMTGACheckContents);
    this->ModelsMTGACheckLayout->addWidget(cb);
    QObject::connect(cb, SIGNAL(stateChanged(int)),
                q, SLOT(onModelsMTGAChanged()));
    QCheckBox* cb2 = new QCheckBox(name, this->ModelsMTGACheckContents);
    this->ModelsCheckLayoutMTGAImg->addWidget(cb2);
    QObject::connect(cb2, SIGNAL(stateChanged(int)),
                q, SLOT(onModelsMTGAImgChanged()));
  }

  for (const QString& name : q->ModelsNamesTCM)
  {
    QCheckBox* cb = new QCheckBox(name, this->ModelsCheckContents);
    this->ModelsCheckLayout->addWidget(cb);
    QObject::connect(cb, SIGNAL(stateChanged(int)),
                q, SLOT(onModelsChanged()));
    QCheckBox* cb2 = new QCheckBox(name, this->ModelsCheckContents);
    this->ModelsCheckLayoutTCMImg->addWidget(cb2);
    QObject::connect(cb2, SIGNAL(stateChanged(int)),
                q, SLOT(onModelsTCMImgChanged()));
  }

  for (int i = 0; i < q->StatsNames.size(); ++i)
  {
    this->StatSelector->addItem(q->StatsNames[i], q->StatsNames[i]);
    this->StatSelectorMTGA->addItem(q->StatsNames[i], q->StatsNames[i]);
    this->StatSelectorImg->addItem(q->StatsNames[i], q->StatsNames[i]);
  }

  this->StatSelector->setCurrentIndex(0);
  this->StatSelectorMTGA->setCurrentIndex(0);
  this->StatSelectorImg->setCurrentIndex(0);

}

void qSlicerKMAPModuleWidgetPrivate::populatePatientComboBox() {
  Q_Q(qSlicerKMAPModuleWidget);

  vtkIdType currentSelectedID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  int currentIndex = this->PatSelector->currentIndex();
  if (currentIndex >= 0)
  {
    currentSelectedID = this->PatSelector->itemData(currentIndex).value<vtkIdType>();
  }

  this->PatSelector->blockSignals(true);  // Optional: prevent signal emission
  this->PatSelector -> clear();
  this->PatSelector->addItem(QString::fromStdString("None"), QVariant::fromValue(vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID));


  vtkMRMLScene* scene = q->mrmlScene();

  if (!scene) {
    q->patID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
    this->populateStudyComboBox(q->patID);
    return;
  }

  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  vtkIdType rootID = shNode->GetSceneItemID();
  if (rootID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID){
    q->patID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
    this->populateStudyComboBox(q->patID);
    return;
  }

  std::function<void(vtkIdType)> visit;
  std::vector<vtkIdType> patients;
  visit = [&](vtkIdType itemID)
  {
    if (itemID != rootID && shNode->HasItemAttribute(itemID, "Level"))
    {
      std :: string level = shNode->GetItemAttribute(itemID, "Level");
      if (level == "Patient")
      {
        patients.push_back(itemID);
      }
    }

    std::vector<vtkIdType> children;
    shNode->GetItemChildren(itemID, children);
    for (vtkIdType childID : children)
    {
      visit(childID);
    }
  };

  visit(rootID);
  int restoredIndex = 0;
  for (vtkIdType id : patients)
  {
    std::string name = shNode->GetItemName(id);
    this->PatSelector->addItem(QString::fromStdString(name), QVariant::fromValue(id));

    if (id == currentSelectedID)
    {
      restoredIndex = this->PatSelector->count() - 1;
    }
  }

  // Restore previous selection if possible
  if (restoredIndex >= 0)
  {
    this->PatSelector->setCurrentIndex(restoredIndex);
  }

  this->PatSelector->blockSignals(false);  // Re-enable signals

  vtkIdType passonID = restoredIndex>0 ? currentSelectedID : vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  q->patID = passonID;
  this->populateStudyComboBox(q->patID);

  return;
}


void qSlicerKMAPModuleWidgetPrivate::populateStudyComboBox(vtkIdType patientID)
{
  Q_Q(qSlicerKMAPModuleWidget);
  // Save current selection by study ID
  vtkIdType currentSelectedStudyID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  int currentIndex = this->StuSelector->currentIndex();
  if (currentIndex >= 0)
  {
    currentSelectedStudyID = this->StuSelector->itemData(currentIndex).value<vtkIdType>();
  }

  this->StuSelector->blockSignals(true);
  this->StuSelector->clear();
  // Add a "None" option first
  this->StuSelector->addItem(QString::fromStdString("None"), QVariant::fromValue(vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID));
  if (patientID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
    this->StuSelector->setEnabled(false);
    // "None" selected — ignore or reset state
    q->stuID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
    this->populateNodeComboBox(this->CTSelector,
                               q->stuID,
                               "vtkMRMLScalarVolumeNode",
                               "CT"
                              );
    this->populateNodeComboBox(this->PETSelector,
                               q->stuID,
                               "vtkMRMLScalarVolumeNode",
                               "PT"
                              );
    this->populateNodeComboBox(this->SegSelector,
                               q->stuID,
                               "vtkMRMLSegmentationNode",
                               ""
                              );
    return;
  }
  this->StuSelector->setEnabled(true);


  vtkMRMLScene* scene = q->mrmlScene();
  if (!scene)
  {
    this->StuSelector->setEnabled(false);
    q->stuID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
    // "None" selected — ignore or reset state
    this->populateNodeComboBox(this->CTSelector,
                               q->stuID,
                               "vtkMRMLScalarVolumeNode",
                               "CT"
                              );
    this->populateNodeComboBox(this->PETSelector,
                               q->stuID,
                               "vtkMRMLScalarVolumeNode",
                               "PT"
                              );
    this->populateNodeComboBox(this->SegSelector,
                               q->stuID,
                               "vtkMRMLSegmentationNode",
                               ""
                              );
    return;
  }

  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (!shNode)
  {
    this->StuSelector->setEnabled(false);
    q->stuID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
    // "None" selected — ignore or reset state
    this->populateNodeComboBox(this->CTSelector,
                               q->stuID,
                               "vtkMRMLScalarVolumeNode",
                               "CT"
                              );
    this->populateNodeComboBox(this->PETSelector,
                               q->stuID,
                               "vtkMRMLScalarVolumeNode",
                               "PT"
                              );
    this->populateNodeComboBox(this->SegSelector,
                               q->stuID,
                               "vtkMRMLSegmentationNode",
                               ""
                              );
    return;
  }

  // Retrieve the children of the given patient
  std::vector<vtkIdType> children;
  shNode->GetItemChildren(patientID, children);

  // Index to restore, default to the "None" option, which is at index 0
  int restoredIndex = 0;
  for (vtkIdType childID : children)
  {
    if (shNode->HasItemAttribute(childID, "Level"))
    {
      std::string level = shNode->GetItemAttribute(childID, "Level");
      if (level == "Study")
      {
        std::string name = shNode->GetItemName(childID);
        this->StuSelector->addItem(QString::fromStdString(name), QVariant::fromValue(childID));

        // Check if this study was previously selected
        if (childID == currentSelectedStudyID)
        {
          restoredIndex = this->StuSelector->count() - 1;
        }
      }
    }
  }

  // Restore the previous selection (if still present)
  this->StuSelector->setCurrentIndex(restoredIndex);
  this->StuSelector->blockSignals(false);

  vtkIdType passonID = restoredIndex>0 ? currentSelectedStudyID : vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  q->stuID = passonID;
  this->populateNodeComboBox(this->CTSelector,
                             q->stuID,
                             "vtkMRMLScalarVolumeNode",
                             "CT"
                            );
  this->populateNodeComboBox(this->PETSelector,
                             q->stuID,
                             "vtkMRMLScalarVolumeNode",
                             "PT"
                            );
  this->populateNodeComboBox(this->SegSelector,
                             q->stuID,
                             "vtkMRMLSegmentationNode",
                             ""
                            );
  return;
}

void qSlicerKMAPModuleWidgetPrivate::populateNodeComboBox(
  QComboBox* comboBox,
  vtkIdType parentItemID,
  const char * requiredNodeType,
  const std :: string requiredModality = ""  // Optional: empty string disables filtering
)
{
  Q_Q(qSlicerKMAPModuleWidget);
  // Save current selection
  vtkIdType currentSelectedItemID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  int currentIndex = comboBox->currentIndex();
  if (currentIndex >= 0)
  {
    currentSelectedItemID = comboBox->itemData(currentIndex).value<vtkIdType>();
  }

  comboBox->blockSignals(true);
  comboBox->clear();
  comboBox->addItem("None", QVariant::fromValue(vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID));

  if (parentItemID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
    comboBox->setEnabled(false);
    comboBox->blockSignals(false);
    if (std::string(requiredNodeType)=="vtkMRMLSegmentationNode") {
      q->segID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
      this->populateSegmentCheckboxes(q->segID);
    } else {
      if (requiredModality=="CT")
        q->ctID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
      if (requiredModality=="PT")
        q->petID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
      q->enableTACbutton();
    }
    return;
  }

  vtkMRMLScene* scene = q->mrmlScene();
  if (!scene)
  {
    comboBox->setEnabled(false);
    comboBox->blockSignals(false);
    if (std::string(requiredNodeType)=="vtkMRMLSegmentationNode") {
      q->segID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
      this->populateSegmentCheckboxes(q->segID);
    } else {
      if (requiredModality=="CT")
        q->ctID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
      if (requiredModality=="PT")
        q->petID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
      q->enableTACbutton();
    }
    return;
  }

  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (!shNode)
  {
    comboBox->setEnabled(false);
    comboBox->blockSignals(false);
    if (std::string(requiredNodeType)=="vtkMRMLSegmentationNode") {
      q->segID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
      this->populateSegmentCheckboxes(q->segID);
    } else {
      if (requiredModality=="CT")
        q->ctID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
      if (requiredModality=="PT")
        q->petID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
      q->enableTACbutton();
    }
    return;
  }

  if (std::string(requiredNodeType)=="vtkMRMLSegmentationNode") {
    if (q->petID==vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID) {
      comboBox->setEnabled(false);
      comboBox->blockSignals(false);
      q->segID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
      this->populateSegmentCheckboxes(q->segID);
      return;
    }
  }
  comboBox->setEnabled(true);
  int restoredIndex = 0;

  std::function<void(vtkIdType)> collectItems;
  collectItems = [&](vtkIdType itemID)
  {
    vtkMRMLNode* dataNode = shNode->GetItemDataNode(itemID);
    if (dataNode && dataNode->IsA(requiredNodeType))
    {

      bool hasatt = shNode->HasItemAttribute(itemID, "DICOM.Modality");
      std :: string modalityAttr = hasatt ? shNode->GetItemAttribute(itemID, "DICOM.Modality") : "";
      bool modalityMatches = requiredModality=="" || requiredModality == modalityAttr;

      if (requiredModality=="PT") {
        std::vector< std::string > pt_attributes = dataNode->GetAttributeNames();
        auto hasseq = find(pt_attributes.begin(), pt_attributes.end(), "Sequences.BaseName");
        modalityMatches = modalityMatches && hasseq!=pt_attributes.end();
      }

      if (modalityMatches)
      {
        std::string name = shNode->GetItemName(itemID);
        comboBox->addItem(QString::fromStdString(name), QVariant::fromValue(itemID));

        if (itemID == currentSelectedItemID)
        {
          restoredIndex = comboBox->count() - 1;
        }
      }
    }

    std::vector<vtkIdType> children;
    shNode->GetItemChildren(itemID, children);
    for (vtkIdType childID : children)
    {
      collectItems(childID);
    }
  };

  collectItems(parentItemID);
  comboBox->setCurrentIndex(restoredIndex);
  comboBox->blockSignals(false);

  vtkIdType passonID = restoredIndex>0 ? currentSelectedItemID : vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  if (std::string(requiredNodeType)=="vtkMRMLSegmentationNode") {
    q->segID = passonID;
    this->populateSegmentCheckboxes(q->segID);
  } else {
    if (requiredModality=="CT")
      q->ctID = passonID;
    if (requiredModality=="PT")
      q->petID = passonID;
    q->enableTACbutton();
  }
}


void qSlicerKMAPModuleWidgetPrivate::populateSegmentCheckboxes(vtkIdType SegItemID)
{
  Q_Q(qSlicerKMAPModuleWidget);

  // Step 1: Save currently selected segment IDs
  QSet<QString> previouslySelectedIDs;
  for (int i = 0; i < this->segmentCheckLayout->count(); ++i)
  {
    QLayoutItem* item = this->segmentCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslySelectedIDs.insert(checkbox->property("SegmentID").toString());
    }
  }

  // Step 2: Clear existing checkboxes
  this->SegmentCheckContents->blockSignals(true);
  QLayoutItem* item;
  while ((item = this->segmentCheckLayout->takeAt(0)) != nullptr)
  {
    delete item->widget();
    delete item;
  }

  if (SegItemID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
    // this->segmentCheckLayout->setEnabled(false);
    this->SegmentCheckContents->blockSignals(false);
    q->segmentIDs.clear();
    q->enableTACbutton();
    this->segmentSelectAll->setEnabled(false);
    return;
  }

  vtkMRMLScene* scene = q->mrmlScene();
  if (!scene)
  {
    this->SegmentCheckContents->blockSignals(false);
    q->segmentIDs.clear();
    q->enableTACbutton();
    this->segmentSelectAll->setEnabled(false);
    return;
  }

  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (!shNode)
  {
    this->SegmentCheckContents->blockSignals(false);
    q->segmentIDs.clear();
    q->enableTACbutton();
    this->segmentSelectAll->setEnabled(false);
    return;
  }

  vtkMRMLSegmentationNode* segNode = vtkMRMLSegmentationNode::SafeDownCast(shNode->GetItemDataNode(SegItemID));
  // Step 3: Repopulate based on new segmentation node
  if (!segNode) {
    this->SegmentCheckContents->blockSignals(false);
    q->segmentIDs.clear();
    q->enableTACbutton();
    this->segmentSelectAll->setEnabled(false);
    return;
  }

  if (q->sequencePETNode == nullptr || q->segSequenceNode == nullptr) {
    this->SegmentCheckContents->blockSignals(false);
    q->segmentIDs.clear();
    q->enableTACbutton();
    this->segmentSelectAll->setEnabled(false);
    return;
  }

  vtkSegmentation* segmentation = segNode->GetSegmentation();
  if (!segmentation) {
    this->SegmentCheckContents->blockSignals(false);
    q->segmentIDs.clear();
    q->enableTACbutton();
    this->segmentSelectAll->setEnabled(false);
    return;
  }

  q->segmentIDs.clear();
  std :: vector<std :: string> segmentIDs = segmentation->GetSegmentIDs();
  for (vtkIdType i = 0; i < segmentIDs.size(); ++i)
  {
    std::string segmentID = segmentIDs[i];
    vtkSegment* vtksegment = segmentation->GetSegment(segmentID);
    std::string segmentName = vtksegment->GetName();

    QCheckBox* checkbox = new QCheckBox(QString::fromStdString(segmentName));
    checkbox->setProperty("SegmentID", QString::fromStdString(segmentID));

    bool wasSelected = previouslySelectedIDs.contains(QString::fromStdString(segmentID));
    checkbox->setChecked(wasSelected);
    this->segmentCheckLayout->addWidget(checkbox);
    QObject::connect(checkbox, SIGNAL(stateChanged(int)),
                 q, SLOT(onSegmentsChanged()));
    if (wasSelected) {
      q->segmentIDs.push_back(QString::fromStdString(segmentID));
    }
  }
  q->enableTACbutton();

  this->segmentCheckLayout->addStretch();
  this->segmentSelectAll->setEnabled(true);
  this->SegmentCheckContents->blockSignals(false);
}


void qSlicerKMAPModuleWidgetPrivate::populatePlotSegmentCheckboxes()
{
  Q_Q(qSlicerKMAPModuleWidget);

  // Step 1: Save currently selected segment IDs
  QSet<QString> previouslyPlotSelectedIDs;
  for (int i = 0; i < this->PlotsegmentCheckLayout->count(); ++i)
  {
    QLayoutItem* item = this->PlotsegmentCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslyPlotSelectedIDs.insert(checkbox->property("SegmentID").toString());
    }
  }

  // Step 2: Clear existing checkboxes
  this->PlotSegmentCheckContents->blockSignals(true);
  QLayoutItem* item;
  while ((item = this->PlotsegmentCheckLayout->takeAt(0)) != nullptr)
  {
    delete item->widget();
    delete item;
  }

  if (q->segmentTACsnames.empty() || q->segmentTACs.empty()) {
    this->TACCollapsibleButton->setEnabled(false);
    this->SegmentCheckContents->blockSignals(false);
    return;
  } else {
    this->TACCollapsibleButton->setEnabled(true);
  }

  for (auto it = q->segmentTACsnames.begin(); it!=q->segmentTACsnames.end(); ++it)
  {
    std::string segmentID = it->first;
    std::string segmentName = it->second;

    QCheckBox* checkbox = new QCheckBox(QString::fromStdString(segmentName));
    checkbox->setProperty("SegmentID", QString::fromStdString(segmentID));

    bool wasSelected = previouslyPlotSelectedIDs.contains(QString::fromStdString(segmentID));
    checkbox->setChecked(wasSelected);
    this->PlotsegmentCheckLayout->addWidget(checkbox);
    // QObject::connect(checkbox, SIGNAL(stateChanged(int)),
    //              q, SLOT(onPlotSegmentsChanged()));
  }

  this->PlotsegmentCheckLayout->addStretch();
  this->SegmentCheckContents->blockSignals(false);
}

void qSlicerKMAPModuleWidgetPrivate::populateTimeBarMTGA() {
  Q_Q(qSlicerKMAPModuleWidget);

  this->timeOffsetSlider->setMinimum(1);
  this->timeOffsetSlider->setMaximum(q->numberOfTimepoints);
  this->timeOffsetSlider->setValue(1);

  frameEdit->setReadOnly(true);
  timeSecEdit->setReadOnly(true);
  timeMinEdit->setReadOnly(true);
  q->onSliderChanged(1);

  QObject::connect( this->timeOffsetSlider, SIGNAL(valueChanged(int)),
    q, SLOT(onSliderChanged(int)));
}

void qSlicerKMAPModuleWidgetPrivate::populateTimeBarMTGAImg() {
  Q_Q(qSlicerKMAPModuleWidget);

  this->timeOffsetSliderImg->setMinimum(1);
  this->timeOffsetSliderImg->setMaximum(q->numberOfTimepoints);
  this->timeOffsetSliderImg->setValue(1);

  frameEdit->setReadOnly(true);
  timeSecEdit->setReadOnly(true);
  timeMinEdit->setReadOnly(true);
  q->onSliderImgChanged(1);

  QObject::connect( this->timeOffsetSliderImg, SIGNAL(valueChanged(int)),
    q, SLOT(onSliderImgChanged(int)));
}

void qSlicerKMAPModuleWidgetPrivate::populateIF()
{
  Q_Q(qSlicerKMAPModuleWidget);

  std :: string currentSelectedID = "";
  int currentIndex = this->IFSelector->currentIndex();
  if (currentIndex >= 0)
  {
    currentSelectedID = this->IFSelector->itemData(currentIndex).toString().toStdString();
  }

  this->IFSelector->blockSignals(true);  // Optional: prevent signal emission
  this->IFSelector->clear();
  this->IFSelector->addItem(QString::fromStdString("None"), QString::fromStdString(""));

  if (q->segmentTACsnames.empty() || q->segmentTACs.empty())
  {
    this->PlotsTabWidget->setTabEnabled(this->PlotsTabWidget->indexOf(this->TCMWidget), false);
    this->TCMWidget->setEnabled(false);
    this->IFSelector->blockSignals(false);
    return;
  }
  this->PlotsTabWidget->setTabEnabled(this->PlotsTabWidget->indexOf(this->TCMWidget), true);
  this->TCMWidget->setEnabled(true);

  int restoredIndex = 0;
  for (const auto& [segmentID, displayName] : q->segmentTACsnames)
  {
    this->IFSelector->addItem(QString::fromStdString(displayName), QString::fromStdString(segmentID));
    if (segmentID==currentSelectedID) {
      restoredIndex = this->IFSelector->count() - 1;
    }
  }

  // Restore previous selection if possible
  if (restoredIndex >= 0)
  {
    this->IFSelector->setCurrentIndex(restoredIndex);
  }

  this->IFSelector->blockSignals(false);

  std :: string passonID = restoredIndex>0 ? currentSelectedID : "";
  q-> IFID = passonID;
  this->populateVOI(q-> IFID);
  return;
}

void qSlicerKMAPModuleWidgetPrivate::populateIFMTGA()
{
  Q_Q(qSlicerKMAPModuleWidget);

  std :: string currentSelectedID = "";
  int currentIndex = this->IFSelectorMTGA->currentIndex();
  if (currentIndex >= 0)
  {
    currentSelectedID = this->IFSelectorMTGA->itemData(currentIndex).toString().toStdString();
  }

  this->IFSelectorMTGA->blockSignals(true);  // Optional: prevent signal emission
  this->IFSelectorMTGA->clear();
  this->IFSelectorMTGA->addItem(QString::fromStdString("None"), QString::fromStdString(""));

  if (q->segmentTACsnames.empty() || q->segmentTACs.empty())
  {
    this->PlotsTabWidget->setTabEnabled(this->PlotsTabWidget->indexOf(this->MTGAWidget), false);
    this->MTGAWidget->setEnabled(false);
    this->IFSelectorMTGA->blockSignals(false);
    return;
  }
  this->PlotsTabWidget->setTabEnabled(this->PlotsTabWidget->indexOf(this->MTGAWidget), true);
  this->MTGAWidget->setEnabled(true);

  int restoredIndex = 0;
  for (const auto& [segmentID, displayName] : q->segmentTACsnames)
  {
    this->IFSelectorMTGA->addItem(QString::fromStdString(displayName), QString::fromStdString(segmentID));
    if (segmentID==currentSelectedID) {
      restoredIndex = this->IFSelectorMTGA->count() - 1;
    }
  }

  // Restore previous selection if possible
  if (restoredIndex >= 0)
  {
    this->IFSelectorMTGA->setCurrentIndex(restoredIndex);
  }

  this->IFSelectorMTGA->blockSignals(false);

  std :: string passonID = restoredIndex>0 ? currentSelectedID : "";
  q-> IFID = passonID;
  this->populateVOIMTGA(q-> IFID);
  return;
}

void qSlicerKMAPModuleWidgetPrivate::populateIFImg()
{
  Q_Q(qSlicerKMAPModuleWidget);

  std :: string currentSelectedID = "";
  int currentIndex = this->IFSelectorImg->currentIndex();
  if (currentIndex >= 0)
  {
    currentSelectedID = this->IFSelectorImg->itemData(currentIndex).toString().toStdString();
  }

  this->IFSelectorImg->blockSignals(true);  // Optional: prevent signal emission
  this->IFSelectorImg->clear();
  this->IFSelectorImg->addItem(QString::fromStdString("None"), QString::fromStdString(""));

  if (q->segmentTACsnames.empty() || q->segmentTACs.empty())
  {
    this->PlotsTabWidget->setTabEnabled(this->PlotsTabWidget->indexOf(this->ImagingWidget), false);
    this->ImagingWidget->setEnabled(false);
    this->IFSelectorImg->blockSignals(false);
    return;
  }
  this->PlotsTabWidget->setTabEnabled(this->PlotsTabWidget->indexOf(this->ImagingWidget), true);
  this->ImagingWidget->setEnabled(true);

  int restoredIndex = 0;
  for (const auto& [segmentID, displayName] : q->segmentTACsnames)
  {
    this->IFSelectorImg->addItem(QString::fromStdString(displayName), QString::fromStdString(segmentID));
    if (segmentID==currentSelectedID) {
      restoredIndex = this->IFSelectorImg->count() - 1;
    }
  }

  // Restore previous selection if possible
  if (restoredIndex >= 0)
  {
    this->IFSelectorImg->setCurrentIndex(restoredIndex);
  }

  this->IFSelectorImg->blockSignals(false);

  return;
}


void qSlicerKMAPModuleWidgetPrivate::populateVOI(std :: string ifID)
{
  Q_Q(qSlicerKMAPModuleWidget);

  // Step 1: Save currently selected segment IDs
  QSet<QString> previouslySelectedIDs;
  for (int i = 0; i < this->VOICheckLayout->count(); ++i)
  {
    QLayoutItem* item = this->VOICheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslySelectedIDs.insert(checkbox->property("SegmentID").toString());
    }
  }

  // Clear previous VOI checkboxes
  this->VOICheckContents->blockSignals(true);
  QLayoutItem* child;
  while ((child = this->VOICheckLayout->takeAt(0)) != nullptr)
  {
    if (child->widget())
    {
      delete child->widget();
    }
    delete child;
  }

  if (ifID=="")
  {
    this->VOICheckContents->blockSignals(false);
    q->VOIsegmentIDs.clear();
    q->enableFITbutton();
    this->VOIsegmentSelectAll->setEnabled(false);
    return;
  }

  // Get the selected IF segment ID
  q->VOIsegmentIDs.clear();

  // Add checkboxes for all other segments
  for (const auto& [segmentID, displayName] : q->segmentTACsnames)
  {
    if (segmentID == ifID)
      continue;

    QCheckBox* cb = new QCheckBox(QString::fromStdString(displayName));
    cb->setProperty("SegmentID", QString::fromStdString(segmentID));
    bool wasSelected = previouslySelectedIDs.contains(QString::fromStdString(segmentID));
    cb->setChecked(wasSelected);
    this->VOICheckLayout->addWidget(cb);
    QObject::connect(cb, SIGNAL(stateChanged(int)),
                     q, SLOT(onVOISegmentsChanged()));
    if (wasSelected)
      q->VOIsegmentIDs.push_back(segmentID);
  }
  q->enableFITbutton();

  this->VOICheckLayout->addStretch();
  this->VOIsegmentSelectAll->setEnabled(true);
  this->VOICheckContents->blockSignals(false);

}

void qSlicerKMAPModuleWidgetPrivate::populateVOIMTGA(std :: string ifID)
{
  Q_Q(qSlicerKMAPModuleWidget);

  // Step 1: Save currently selected segment IDs
  QSet<QString> previouslySelectedIDs;
  for (int i = 0; i < this->VOIMTGACheckLayout->count(); ++i)
  {
    QLayoutItem* item = this->VOIMTGACheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslySelectedIDs.insert(checkbox->property("SegmentID").toString());
    }
  }

  // Clear previous VOI checkboxes
  this->VOIMTGACheckContents->blockSignals(true);
  QLayoutItem* child;
  while ((child = this->VOIMTGACheckLayout->takeAt(0)) != nullptr)
  {
    if (child->widget())
    {
      delete child->widget();
    }
    delete child;
  }

  if (ifID=="")
  {
    this->VOIMTGACheckContents->blockSignals(false);
    q->VOIMTGAsegmentIDs.clear();
    q->enableFITMTGAbutton();
    this->VOIMTGAsegmentSelectAll->setEnabled(false);
    return;
  }

  // Get the selected IF segment ID
  q->VOIMTGAsegmentIDs.clear();

  // Add checkboxes for all other segments
  for (const auto& [segmentID, displayName] : q->segmentTACsnames)
  {
    if (segmentID == ifID)
      continue;

    QCheckBox* cb = new QCheckBox(QString::fromStdString(displayName));
    cb->setProperty("SegmentID", QString::fromStdString(segmentID));
    bool wasSelected = previouslySelectedIDs.contains(QString::fromStdString(segmentID));
    cb->setChecked(wasSelected);
    this->VOIMTGACheckLayout->addWidget(cb);
    QObject::connect(cb, SIGNAL(stateChanged(int)),
                     q, SLOT(onVOIMTGASegmentsChanged()));
    if (wasSelected)
      q->VOIMTGAsegmentIDs.push_back(segmentID);
  }
  q->enableFITMTGAbutton();

  this->VOIMTGACheckLayout->addStretch();
  this->VOIMTGAsegmentSelectAll->setEnabled(true);
  this->VOIMTGACheckContents->blockSignals(false);

}


void qSlicerKMAPModuleWidgetPrivate::populateResultsVOI()
{
  Q_Q(qSlicerKMAPModuleWidget);

  std :: string currentSelectedID = "";
  int currentIndex = this->VOISelector->currentIndex();
  if (currentIndex >= 0)
  {
    currentSelectedID = this->VOISelector->itemData(currentIndex).toString().toStdString();
  }

  this->VOISelector->blockSignals(true);  // Optional: prevent signal emission
  this->VOISelector->clear();
  this->VOISelector->addItem(QString::fromStdString("None"), QString::fromStdString(""));

  if (q->segmentTCM.empty() || q->segmentTACsnames.empty() || q->segmentTACs.empty()) {
    this->TCMResultsButton->setEnabled(false);
    this->VOISelector->blockSignals(false);
    this->populateResultsTable("");
    return;
  }
  this->TCMResultsButton->setEnabled(true);

  int restoredIndex = 0;
  for (const auto& [segmentID, _] : q->segmentTCM)
  {
    this->VOISelector->addItem(QString::fromStdString(q->segmentTACsnames[segmentID]), QString::fromStdString(segmentID));
    if (segmentID==currentSelectedID) {
      restoredIndex = this->VOISelector->count() - 1;
    }
  }

  // Restore previous selection if possible
  if (restoredIndex >= 0)
  {
    this->VOISelector->setCurrentIndex(restoredIndex);
  }

  this->VOISelector->blockSignals(false);

  std :: string passonID = restoredIndex>0 ? currentSelectedID : "";
  q->plotTCMVOI = passonID;
  this->populateResultsTable(passonID);
  return;
}

void qSlicerKMAPModuleWidgetPrivate::populateResultsVOIMTGA()
{
  Q_Q(qSlicerKMAPModuleWidget);

  std :: string currentSelectedID = "";
  int currentIndex = this->VOISelectorMTGA->currentIndex();
  if (currentIndex >= 0)
  {
    currentSelectedID = this->VOISelectorMTGA->itemData(currentIndex).toString().toStdString();
  }

  this->VOISelectorMTGA->blockSignals(true);  // Optional: prevent signal emission
  this->VOISelectorMTGA->clear();
  this->VOISelectorMTGA->addItem(QString::fromStdString("None"), QString::fromStdString(""));

  if (q->segmentMTGA.empty() || q->segmentTACsnames.empty() || q->segmentTACs.empty()) {
    this->MTGAResultsButton->setEnabled(false);
    this->VOISelectorMTGA->blockSignals(false);
    this->populateResultsMTGATable("");
    return;
  }
  this->MTGAResultsButton->setEnabled(true);

  int restoredIndex = 0;
  for (const auto& [segmentID, _] : q->segmentMTGA)
  {
    this->VOISelectorMTGA->addItem(QString::fromStdString(q->segmentTACsnames[segmentID]), QString::fromStdString(segmentID));
    if (segmentID==currentSelectedID) {
      restoredIndex = this->VOISelectorMTGA->count() - 1;
    }
  }

  // Restore previous selection if possible
  if (restoredIndex >= 0)
  {
    this->VOISelectorMTGA->setCurrentIndex(restoredIndex);
  }

  this->VOISelectorMTGA->blockSignals(false);

  std :: string passonID = restoredIndex>0 ? currentSelectedID : "";
  q->plotMTGAVOI = passonID;
  this->populateResultsMTGATable(passonID);
  return;
}

auto makeNumericItem = [](double value, int precision = 6) {
    auto *item = new QTableWidgetItem(QString::number(value));
    item->setData(Qt::EditRole, value);  // ensures numeric sorting
    item->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
    return item;
};

void qSlicerKMAPModuleWidgetPrivate::populateResultsTable(std :: string segmentID)
{
  Q_Q(qSlicerKMAPModuleWidget);

  this->TCMResultsTable->clear();
  this->TCMResultsTable->setRowCount(0);
  this->TCMResultsTable->setColumnCount(0);

  if (segmentID.empty())
  {
    this->populateModelsTCM(segmentID);
    this->populateModelComboTCM(this->TCMModel1, "", "", segmentID);
    this->populateModelComboTCM(this->TCMModel2, "", "", segmentID);
    this->TCMLRTP->setText("");
    this->TCMVuongP->setText("");
    return;
  }

  // Define column headers
  QStringList headers = { "", "K1", "k2", "k3", "k4", "vb", "td", "Ki", "DV", "AIC", "BIC", "MASE", "chi^2_nu"};
  this->TCMResultsTable->setColumnCount(headers.size());
  this->TCMResultsTable->setHorizontalHeaderLabels(headers);

  // Get the parameter map for the selected segment
  const auto& labelMap = q->segmentTCM[segmentID];
  int row = 0;
  this->TCMResultsTable->setRowCount(labelMap.size());

  int totalRowHeight = 0;
  this->TCMResultsTable->resizeRowsToContents();
  for (const auto& [label, params] : labelMap)
  {
    this->TCMResultsTable->setItem(row, 0, new QTableWidgetItem(QString::fromStdString(label)));
    this->TCMResultsTable->setItem(row, 1, makeNumericItem(params.K1));
    this->TCMResultsTable->setItem(row, 2, makeNumericItem(params.k2));
    this->TCMResultsTable->setItem(row, 3, makeNumericItem(params.k3));
    this->TCMResultsTable->setItem(row, 4, makeNumericItem(params.k4));
    this->TCMResultsTable->setItem(row, 5, makeNumericItem(params.vb));
    this->TCMResultsTable->setItem(row, 6, makeNumericItem(params.td));
    this->TCMResultsTable->setItem(row, 7, makeNumericItem(params.Ki));
    this->TCMResultsTable->setItem(row, 8, makeNumericItem(params.DV));
    this->TCMResultsTable->setItem(row, 9, makeNumericItem(params.AIC));
    this->TCMResultsTable->setItem(row, 10, makeNumericItem(params.BIC));
    this->TCMResultsTable->setItem(row, 11, makeNumericItem(params.MASE));
    this->TCMResultsTable->setItem(row, 12, makeNumericItem(params.chi2));

    totalRowHeight += this->TCMResultsTable->rowHeight(row);
    ++row;
  }
  totalRowHeight += this->TCMResultsTable->horizontalHeader()->height();
  totalRowHeight += 2 * this->TCMResultsTable->frameWidth();
  this->TCMResultsTable->setMinimumHeight(totalRowHeight);
  this->TCMResultsTable->setMaximumHeight(totalRowHeight);
  // Make table read-only
  this->TCMResultsTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
  this->TCMResultsTable->resizeColumnsToContents();
  this->populateModelsTCM(segmentID);

  std::string sel1, sel2;
  int idx1 = this->TCMModel1->currentIndex();
  if (idx1 >= 0)
    sel1 = this->TCMModel1->itemData(idx1).toString().toStdString();
  int idx2 = this->TCMModel2->currentIndex();
  if (idx2 >= 0)
    sel2 = this->TCMModel2->itemData(idx2).toString().toStdString();
  this->populateModelComboTCM(this->TCMModel1, sel2, sel1, segmentID);
  this->populateModelComboTCM(this->TCMModel2, sel1, sel2, segmentID);
  if (idx1 > 0 & idx2 >0){
    q->runTCMstat(sel1, sel2, segmentID);
  } else {
    this->TCMLRTP->setText("");
    this->TCMVuongP->setText("");
  }
}

void qSlicerKMAPModuleWidgetPrivate::populateResultsMTGATable(std :: string segmentID)
{
  Q_Q(qSlicerKMAPModuleWidget);

  this->MTGAResultsTable->clear();
  this->MTGAResultsTable->setRowCount(0);
  this->MTGAResultsTable->setColumnCount(0);

  if (segmentID.empty())
  {
    this->populateModelsMTGA(segmentID);
    this->populateModelCombo(this->MTGAModel1, "", "", segmentID);
    this->populateModelCombo(this->MTGAModel2, "", "", segmentID);
    this->MTGAVuongP->setText("");
    return;
  }

  // Define column headers
  QStringList headers = { "", "Ki", "DV", "Intercept", "R2", "AIC", "MASE"};
  this->MTGAResultsTable->setColumnCount(headers.size());
  this->MTGAResultsTable->setHorizontalHeaderLabels(headers);

  // Get the parameter map for the selected segment
  const auto& labelMap = q->segmentMTGA[segmentID];
  int row = 0;
  this->MTGAResultsTable->setRowCount(labelMap.size());

  int totalRowHeight = 0;
  this->MTGAResultsTable->resizeRowsToContents();
  for (const auto& [label, params] : labelMap)
  {
    this->MTGAResultsTable->setItem(row, 0, new QTableWidgetItem(QString::fromStdString(label)));
    this->MTGAResultsTable->setItem(row, 1, makeNumericItem(params.Ki));
    this->MTGAResultsTable->setItem(row, 2, makeNumericItem(params.DV));
    this->MTGAResultsTable->setItem(row, 3, makeNumericItem(params.Intercept));
    this->MTGAResultsTable->setItem(row, 4, makeNumericItem(params.R2));
    this->MTGAResultsTable->setItem(row, 5, makeNumericItem(params.AIC));
    this->MTGAResultsTable->setItem(row, 6, makeNumericItem(params.MASE));
    // this->MTGAResultsTable->setItem(row, 7, makeNumericItem(params.chi2));

    totalRowHeight += this->MTGAResultsTable->rowHeight(row);
    ++row;
  }
  totalRowHeight += this->MTGAResultsTable->horizontalHeader()->height();
  totalRowHeight += 4 * this->MTGAResultsTable->frameWidth();
  this->MTGAResultsTable->setMinimumHeight(totalRowHeight);
  this->MTGAResultsTable->setMaximumHeight(totalRowHeight);
  // Make table read-only
  this->MTGAResultsTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
  this->MTGAResultsTable->resizeColumnsToContents();
  this->populateModelsMTGA(segmentID);

  std::string sel1, sel2;
  int idx1 = this->MTGAModel1->currentIndex();
  if (idx1 >= 0)
    sel1 = this->MTGAModel1->itemData(idx1).toString().toStdString();
  int idx2 = this->MTGAModel2->currentIndex();
  if (idx2 >= 0)
    sel2 = this->MTGAModel2->itemData(idx2).toString().toStdString();
  this->populateModelCombo(this->MTGAModel1, sel2, sel1, segmentID);
  this->populateModelCombo(this->MTGAModel2, sel1, sel2, segmentID);
  if (idx1 > 0 & idx2 >0) {
    q->runVuong(sel1, sel2, segmentID);
  } else {
    this->MTGAVuongP->setText("");
  }
}


void qSlicerKMAPModuleWidgetPrivate::populateModelsTCM(std :: string segmentID)
{
  Q_Q(qSlicerKMAPModuleWidget);
  // Step 1: Save currently selected segment IDs
  QSet<QString> previouslySelectedIDs;
  for (int i = 0; i < this->ModelsTCMCheckLayout->count(); ++i)
  {
    QLayoutItem* item = this->ModelsTCMCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslySelectedIDs.insert(checkbox->text());
    }
  }

  // Clear previous VOI checkboxes
  this->ModelsTCMCheckContents->blockSignals(true);
  QLayoutItem* child;
  while ((child = this->ModelsTCMCheckLayout->takeAt(0)) != nullptr)
  {
    if (child->widget())
    {
      delete child->widget();
    }
    delete child;
  }


  if (segmentID.empty())
  {
    this->ModelsTCMCheckContents->blockSignals(false);
    this->plotTCMButton->setEnabled(false);
    this->ModelsTCMSelectAll->setEnabled(false);
    return;
  }

  // Get the parameter map for the selected segment
  const auto& labelMap = q->segmentTCM[segmentID];
  for (const auto& [label, params] : labelMap) {
    QCheckBox* cb = new QCheckBox(QString::fromStdString(label));
    bool wasSelected = previouslySelectedIDs.contains(QString::fromStdString(label));
    cb->setChecked(wasSelected);
    this->ModelsTCMCheckLayout->addWidget(cb);
  }
  this->plotTCMButton->setEnabled(true);

  this->ModelsTCMCheckLayout->addStretch();
  this->ModelsTCMSelectAll->setEnabled(true);
  this->ModelsTCMCheckContents->blockSignals(false);
  q->onPlotTCMbutton();
}

void qSlicerKMAPModuleWidgetPrivate::populateModelsMTGA(std :: string segmentID)
{

  Q_Q(qSlicerKMAPModuleWidget);

  std :: string currentSelectedID = "";
  int currentIndex = this->MTGASelector->currentIndex();
  if (currentIndex >= 0)
  {
    currentSelectedID = this->MTGASelector->itemData(currentIndex).toString().toStdString();
  }

  this->MTGASelector->blockSignals(true);  // Optional: prevent signal emission
  this->MTGASelector->clear();
  this->MTGASelector->addItem(QString::fromStdString("None"), QString::fromStdString(""));

  if (segmentID.empty())
  {
    this->MTGASelector->blockSignals(false);
    this->plotMTGAButton->setEnabled(false);
    return;
  }

  const auto& labelMap = q->segmentMTGA[segmentID];
  int restoredIndex = 0;
  for (const auto& [label, params] : labelMap)
  {
    this->MTGASelector->addItem(QString::fromStdString(label), QString::fromStdString(label));
    if (label==currentSelectedID) {
      restoredIndex = this->MTGASelector->count() - 1;
    }
  }

  // Restore previous selection if possible
  if (restoredIndex >= 0)
  {
    this->MTGASelector->setCurrentIndex(restoredIndex);
  }
  this->plotMTGAButton->setEnabled(true);

  this->MTGASelector->blockSignals(false);

  std :: string passonID = restoredIndex>0 ? currentSelectedID : "";
  q->plotMTGAModel = passonID;
  q->onPlotMTGAbutton();
}

void qSlicerKMAPModuleWidgetPrivate::populateModelCombo(
    QComboBox* comboToFill,
    const std::string& otherSelectedModel,
    const std::string& currentSelectedModel,
    const std::string& segmentID)
{
  Q_Q(qSlicerKMAPModuleWidget);

  if (segmentID.empty()) {
    this->MTGAModel1->clear();
    this->MTGAModel2->clear();
    this->MTGAVuongP->setText("");
    this->MTGAStatTestButton->setCollapsed(true);
    this->MTGAStatTestButton->setEnabled(false);
    return;
  }

  auto it = q->segmentMTGA.find(segmentID);
  if (it == q->segmentMTGA.end()) {
    this->MTGAModel1->clear();
    this->MTGAModel2->clear();
    this->MTGAVuongP->setText("");
    this->MTGAStatTestButton->setCollapsed(true);
    this->MTGAStatTestButton->setEnabled(false);
    return;
  }
  this->MTGAStatTestButton->setEnabled(true);
  const auto& modelsForSegment = it->second;

  comboToFill->blockSignals(true);
  comboToFill->clear();
  comboToFill->addItem("", "");  // empty choice

  int restoredIndex = 0;
  for (const auto& [modelName, params] : modelsForSegment)
  {
    if (!otherSelectedModel.empty() && modelName == otherSelectedModel)
      continue;  // skip what’s selected in the other box

    comboToFill->addItem(QString::fromStdString(modelName), QString::fromStdString(modelName));

    if (modelName == currentSelectedModel)
    {
      restoredIndex = comboToFill->count() - 1;
    }
  }

  comboToFill->setCurrentIndex(restoredIndex);
  comboToFill->blockSignals(false);
}

void qSlicerKMAPModuleWidgetPrivate::populateModelComboTCM(
    QComboBox* comboToFill,
    const std::string& otherSelectedModel,
    const std::string& currentSelectedModel,
    const std::string& segmentID)
{
  Q_Q(qSlicerKMAPModuleWidget);

  if (segmentID.empty()) {
    this->TCMModel1->clear();
    this->TCMModel2->clear();
    this->TCMLRTP->setText("");
    this->TCMVuongP->setText("");
    this->TCMStatTestButton->setCollapsed(true);
    this->TCMStatTestButton->setEnabled(false);
    return;
  }

  auto it = q->segmentTCM.find(segmentID);
  if (it == q->segmentTCM.end()) {
    this->TCMModel1->clear();
    this->TCMModel2->clear();
    this->TCMLRTP->setText("");
    this->TCMVuongP->setText("");
    this->TCMStatTestButton->setCollapsed(true);
    this->TCMStatTestButton->setEnabled(false);
    return;
  }
  this->TCMStatTestButton->setEnabled(true);
  const auto& modelsForSegment = it->second;

  comboToFill->blockSignals(true);
  comboToFill->clear();
  comboToFill->addItem("", "");  // empty choice

  int restoredIndex = 0;
  for (const auto& [modelName, params] : modelsForSegment)
  {
    if (!otherSelectedModel.empty() && modelName == otherSelectedModel)
      continue;  // skip what’s selected in the other box

    comboToFill->addItem(QString::fromStdString(modelName), QString::fromStdString(modelName));

    if (modelName == currentSelectedModel)
    {
      restoredIndex = comboToFill->count() - 1;
    }
  }

  comboToFill->setCurrentIndex(restoredIndex);
  comboToFill->blockSignals(false);
}



//-----------------------------------------------------------------------------
// qSlicerKMAPModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerKMAPModuleWidget::qSlicerKMAPModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerKMAPModuleWidgetPrivate(*this) )
{
  Q_D(qSlicerKMAPModuleWidget);
  this->SubjectHierarchyNode = nullptr;
  this->patID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  this->stuID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  this->ctID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  this->petID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  this->segID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  this->sequencePETNode = nullptr;
  this->segSequenceNode = nullptr;
  this->sequenceBrowserPETNode = nullptr;
  this->numberOfTimepoints = 0;
  QMainWindow* mainWindow = qobject_cast<QMainWindow*>(qSlicerApplication::application()->mainWindow());
  this->ProgressBar = new QProgressBar(mainWindow);
  this->ProgressBar->setObjectName(QString::fromUtf8("ProgressBar"));
  this->ProgressBar->setMaximum(100);
  this->ProgressBar->setValue(0);
  this->ProgressBar->resize(300, 30);
  QRect parentGeometry = this->ProgressBar->parentWidget()->geometry();
  QRect barGeometry = this->ProgressBar->frameGeometry();
  this->ProgressBar->move(
      (parentGeometry.width() - barGeometry.width()) / 2,
      (parentGeometry.height() - barGeometry.height()) / 2
  );
  this->stopButton = new QPushButton("Stop", mainWindow);
  stopButton->move(this->ProgressBar->x(), this->ProgressBar->y() + 40);

  this->stopButton->setStyleSheet(
      "QPushButton {"
      "  background-color: red;"
      "  color: white;"        // Text color
      "  font-weight: bold;"
      "}"
      "QPushButton:pressed {"
      "  background-color: darkred;"
      "}"
  );
  this->stopRequested = false;
  QObject::connect(stopButton, &QPushButton::clicked, [this]() {
      this->stopRequested = true;
  });

  this->checkboxNames = QStringList{
    "Mean", "Median", "Min", "Max"//, "VoxelCount", "Volume(cc)"
  };
  this->ModelsNamesMTGA = QStringList{
    "Patlak", "Logan", "RE"
  };
  this->ModelsNamesTCM = QStringList{
    "1TCM", "1TdCM", "1TiCM", "1TidCM", "2TCM", "2dTCM", "2TiCM", "2TidCM"
  };
  this->StatsNames = QStringList{
    "Mean", "Median"
  };
  this->PlotSelectedFrame = -1;
  this->PlotSelectedVOI = "";
  // Install global key watcher
  if (mainWindow)
  {
    this->keyWatcher = new KeyPressWatcher(mainWindow);
    mainWindow->installEventFilter(this->keyWatcher);
    QObject::connect(this->keyWatcher, SIGNAL(deletePressed()), this, SLOT(onDeleteKeyPressed()));
  }

  // Optional: start with watcher disabled
  this->keyWatcher->setActive(true);
  this->PET_flatten_values.clear();
  this->PETdims[0] = 0;
  this->PETdims[1] = 0;
  this->PETdims[2] = 0;
  // qRegisterMetaType<std::vector<TCMParameters>>("std::vector<TCMParameters>");
  d->init();
}

//-----------------------------------------------------------------------------
qSlicerKMAPModuleWidget::~qSlicerKMAPModuleWidget()
{
}

// void qSlicerKMAPModuleWidget::setNodeSelectorEnabled(qMRMLNodeComboBox* selector, bool enabled)
// {
//   selector->setEnabled(enabled);
//   const auto children = selector->findChildren<QWidget*>();
//   for (QWidget* child : children)
//   {
//     child->setEnabled(enabled);
//   }
// }

// std::map<std::string, vtkIdType> qSlicerKMAPModuleWidget::GetStudyAndPatientAncestors(
//   vtkMRMLSubjectHierarchyNode* shNode,
//   vtkIdType itemID)
// {
//   std::map<std::string, vtkIdType> result;
//
//   if (shNode==nullptr)
//     return result;
//
//   if (itemID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
//     return result;
//
//   vtkIdType currentID = shNode->GetItemParent(itemID);  // skip self
//   while (currentID != vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
//   {
//     // Check and retrieve "Level" attribute
//     if (shNode->HasItemAttribute(currentID, "Level"))
//     {
//       std :: string levelStr = shNode->GetItemAttribute(currentID, "Level");
//       if (levelStr == "Study" || levelStr == "Patient")
//       {
//         // Insert if not already in map (closer ancestor takes precedence)
//         if (result.find(levelStr) == result.end())
//         {
//           result[levelStr] = currentID;
//         }
//
//         // Optional: stop once both found
//         if (result.size() == 2)
//           break;
//       }
//     }
//
//     currentID = shNode->GetItemParent(currentID);
//   }
//
//   return result;
// }


// void qSlicerKMAPModuleWidget::showProgressBar()
// {
//   this->ui->progressBar->setVisible(true);
//   this->ui->progressBar->setValue(0);
// }
//
// void qSlicerKMAPModuleWidget::updateProgress(double progress)
// {
//   this->ui->progressBar->setValue(static_cast<int>(progress * 100.0));
//   qApp->processEvents(); // To force GUI update
// }
//
// void qSlicerKMAPModuleWidget::hideProgressBar()
// {
//   this->ui->progressBar->setVisible(false);
// }

void qSlicerKMAPModuleWidget::enter()
{
  this->IsActive = true;
  this->Superclass::enter();

  // Optional: force refresh when the user enters the module
  this->onSubjectHierarchyChanged();
}

void qSlicerKMAPModuleWidget::exit()
{
  this->IsActive = false;
  this->Superclass::exit();
}


void qSlicerKMAPModuleWidget::onSubjectHierarchyChanged() {
  if (!this->IsActive)
  {
    return;  // Don't do anything if the module is not active
  }
  Q_D(qSlicerKMAPModuleWidget);
  d->populatePatientComboBox();
}

void qSlicerKMAPModuleWidget::setMRMLScene(vtkMRMLScene* scene) {
  this->Superclass::setMRMLScene(scene);
  Q_D(qSlicerKMAPModuleWidget);

  this->qvtkDisconnectAll();

  this->SubjectHierarchyNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (this->SubjectHierarchyNode)
  {
    this->qvtkConnect(this->SubjectHierarchyNode, vtkMRMLSubjectHierarchyNode::SubjectHierarchyItemAddedEvent,
                      this, SLOT(onSubjectHierarchyChanged()));
    this->qvtkConnect(this->SubjectHierarchyNode, vtkMRMLSubjectHierarchyNode::SubjectHierarchyItemRemovedEvent,
                      this, SLOT(onSubjectHierarchyChanged()));
    this->qvtkConnect(this->SubjectHierarchyNode, vtkMRMLSubjectHierarchyNode::SubjectHierarchyItemModifiedEvent,
                      this, SLOT(onSubjectHierarchyChanged()));
  }
  this->SegWatcher = vtkSmartPointer<SegmentationChangeWatcher>::New();
  this->SegWatcher->GetSequencePET = [this]() { return this->sequencePETNode; };
  this->SegWatcher->GetLogic = [this]() { return vtkSlicerKMAPLogic::SafeDownCast(this->logic()); };
  this->SegWatcher->GetsegmentTACs = [this]() { return &this->segmentTACs; };
  this->SegWatcher->GetSegEditCorr = [d]() { return d->PlotLiveSegEdit->isChecked(); };
  this->SegWatcher->RunPlot = [this]() { this->onPlotbutton(); };
  this->SegWatcher->GetCurrentSegID = [this]() { return this->SubjectHierarchyNode->GetItemName(this->segID); };

}

void qSlicerKMAPModuleWidget::getDurations()
{
  if (!this->durations.empty() || !this->timePoints.empty()) {
    return;
  }

  Q_D(qSlicerKMAPModuleWidget);
  vtkMRMLScene* scene = this->mrmlScene();
  if (!scene) {
    return;
  }
  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (!shNode)
  {
    return;
  }

  vtkMRMLScalarVolumeNode* petNode = vtkMRMLScalarVolumeNode::SafeDownCast(shNode->GetItemDataNode(this->petID));
  if (!petNode) {
    this->durations.clear();
    this->timePoints.clear();
    return;
  }
  this->petNode = petNode;

  ctkDICOMDatabase* db = qSlicerApplication::application()->dicomDatabase();

  vtkIdType shItemID = shNode->GetItemByDataNode(petNode);
  if (shItemID == vtkMRMLSubjectHierarchyNode::GetInvalidItemID())
  {
    qWarning() << "No SH item for PET node\n!";
  }
  std::string uidListString = shNode->GetItemUID(
    shItemID,
    vtkMRMLSubjectHierarchyConstants::GetDICOMInstanceUIDName()
  );
  std::vector<std::string> uidVector;
  vtkMRMLSubjectHierarchyNode::DeserializeUIDList(uidListString, uidVector);
  if (uidVector.empty())
  {
      qWarning() << "No DICOMInstanceUID found for PET node!";
      return;
  }
  QString instanceUID = QString::fromStdString(uidVector[0]);

  QString seriesUID = db->seriesForFile(db->fileForInstance(instanceUID));
  QStringList fileList = db->filesForSeries(seriesUID);

  if (fileList.empty())
    throw std::runtime_error(
        "No files found for " + seriesUID.toStdString()
    );

  std::map<double, double> timeSorteddurations;
  for (const QString& file : fileList)
  {
    DcmFileFormat fileformat;
    OFCondition status = fileformat.loadFile(file.toStdString().c_str(), EXS_Unknown, EGL_noChange, DCM_MaxReadLength, ERM_autoDetect);

    if (!status.good())
    {
      std::cerr << "Cannot read file: " << file.toStdString() << std::endl;
      continue;
    }

    DcmDataset* dataset = fileformat.getDataset();

    // 1. Get acquisition datetime
    OFString acqTimeStr;
    if (dataset->findAndGetOFString(DCM_AcquisitionDateTime, acqTimeStr).bad())
    {
      std::cerr << "Missing AcquisitionDateTime" << std::endl;
      continue;
    }

    // Convert to sortable float seconds since start
    struct tm tmTime = {};
    const char* dateTimeCStr = acqTimeStr.c_str();
    if (!strptime(dateTimeCStr, "%Y%m%d%H%M%S", &tmTime))
    {
      std::cerr << "Failed to parse time: " << dateTimeCStr << std::endl;
      continue;
    }
    time_t epochTime = mktime(&tmTime);
    double timeInSeconds = static_cast<double>(epochTime);

    // 2. Get frame duration
    double durationSec = -1.0;

    Float64 durVal;
    if (dataset->findAndGetFloat64(DcmTagKey(0x0018, 0x1242), durVal).good())
    {
      durationSec = durVal / 1000.0;
    }
    else if (dataset->findAndGetFloat64(DcmTagKey(0x0067, 0x1004), durVal).good())
    {
      durationSec = durVal;  // already in seconds
    }
    else
    {
      std::cerr << "No known duration tag in: " << file.toStdString() << std::endl;
      continue;
    }

    timeSorteddurations[timeInSeconds] = durationSec;
  }

  this->durations.clear();
  this->timePoints.clear();
  this->durations.reserve(timeSorteddurations.size());
  this->timePoints.reserve(timeSorteddurations.size());

  double timealong = 0.0;
  for (const auto& pair : timeSorteddurations)
  {
    timealong += pair.second;
    this->durations.push_back(pair.second);
    this->timePoints.push_back(timealong);
  }

}



void qSlicerKMAPModuleWidget::onPatChanged (int index) {
  Q_D(qSlicerKMAPModuleWidget);
  this->patID = d->PatSelector->itemData(index).value<vtkIdType>();
  d->populateStudyComboBox(this->patID);
  if (this->patID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
    // "None" selected — ignore or reset state
    d->fileexcel->setText(QString::fromStdString(".xlsx"));
    return;
  }
  std :: string name = this->SubjectHierarchyNode->GetItemName(this->patID);
  std :: string excelfile = name + "_TAC.xlsx";
  d->fileexcel->setText(QString::fromStdString(excelfile));
  std :: string excelfiletcm = name + "_TCMparameters.xlsx";
  d->fileexceltcm->setText(QString::fromStdString(excelfiletcm));
  std :: string excelfiletcmfitted = name + "_TCMfitted.xlsx";
  d->fileexceltcmfitted->setText(QString::fromStdString(excelfiletcmfitted));
  std :: string excelfilemtga = name + "_MTGAparameters.xlsx";
  d->fileexcelmtga->setText(QString::fromStdString(excelfilemtga));
  std :: string excelfilemtgafitted = name + "_MTGAfitted.xlsx";
  d->fileexcelmtgafitted->setText(QString::fromStdString(excelfilemtgafitted));
}


void qSlicerKMAPModuleWidget::onStuChanged (int index) {
  Q_D(qSlicerKMAPModuleWidget);
  this->stuID = d->StuSelector->itemData(index).value<vtkIdType>();
  d->populateNodeComboBox(d->CTSelector,
                          this->stuID,
                          "vtkMRMLScalarVolumeNode",
                          "CT"
                          );
  d->populateNodeComboBox(d->PETSelector,
                          this->stuID,
                          "vtkMRMLScalarVolumeNode",
                          "PT"
                          );
  d->populateNodeComboBox(d->SegSelector,
                          this->stuID,
                          "vtkMRMLSegmentationNode",
                          ""
                          );
}


void qSlicerKMAPModuleWidget::onCTChanged (int index) {
  Q_D(qSlicerKMAPModuleWidget);
  this->ctID = d->CTSelector->itemData(index).value<vtkIdType>();
  this->enableTACbutton();
}


void qSlicerKMAPModuleWidget::onPETChanged (int index) {
  Q_D(qSlicerKMAPModuleWidget);
  this->durations.clear();
  this->timePoints.clear();
  this->petID = d->PETSelector->itemData(index).value<vtkIdType>();

  this->sequencePETNode = nullptr;
  this->sequenceBrowserPETNode = nullptr;
  this->SegWatcher->browser = nullptr;

  vtkMRMLScene* scene = this->mrmlScene();
  if (scene==nullptr) {
    return;
  }
  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (!shNode) {
    return;
  }
  // Fetch PET
  if (this->petID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID) {
    d->populateNodeComboBox(d->SegSelector,
                            this->stuID,
                            "vtkMRMLSegmentationNode",
                            ""
                            );
    return;
  }
  vtkMRMLScalarVolumeNode* petNode = vtkMRMLScalarVolumeNode::SafeDownCast(shNode->GetItemDataNode(this->petID));
  if (!petNode) {
    d->populateNodeComboBox(d->SegSelector,
                            this->stuID,
                            "vtkMRMLSegmentationNode",
                            ""
                            );
    return;
  }

  // Collect the sequence for the dynamic PET
  for (int i = 0; i < scene->GetNumberOfNodesByClass("vtkMRMLSequenceBrowserNode"); ++i)
  {
    vtkMRMLSequenceBrowserNode* browser = vtkMRMLSequenceBrowserNode::SafeDownCast(scene->GetNthNodeByClass(i, "vtkMRMLSequenceBrowserNode"));
    if (!browser)
      continue;

    // Check if this browser is using our PET node as a proxy node
    vtkMRMLSequenceNode* seqNode = browser->GetSequenceNode(petNode);
    if (seqNode)
    {
      this->sequencePETNode = seqNode;
      this->sequenceBrowserPETNode = browser;
      this->SegWatcher->browser = browser;
      break;
    }
  }
  if (!this->sequencePETNode || !this->sequenceBrowserPETNode)
  {
    QMessageBox::warning(nullptr,
                         tr("Missing node"),
                         tr("Could not find sequence or browser node for PET."));
    return;
  }
  this->numberOfTimepoints = this->sequencePETNode->GetNumberOfDataNodes();

  d->populateNodeComboBox(d->SegSelector,
                          this->stuID,
                          "vtkMRMLSegmentationNode",
                          ""
                          );

  this->getDurations();
  this->enableTACbutton();

  // Get proxy node for current time/frame
  this->sequenceBrowserPETNode->SetSelectedItemNumber(this->numberOfTimepoints - 1);
  vtkMRMLScalarVolumeNode* proxyVolume =
      vtkMRMLScalarVolumeNode::SafeDownCast(this->sequenceBrowserPETNode->GetProxyNode(this->sequencePETNode));
  if (!proxyVolume)
  {
      qCritical() << "Cannot get proxy volume for PET frame";
      return;
  }

  vtkMRMLProceduralColorNode* petLUTNode =
      vtkMRMLProceduralColorNode::SafeDownCast(
          scene->GetFirstNodeByName("PET-DICOM"));
  if (!petLUTNode)
  {
      qCritical() << "Could not find PET-DICOM procedural color node in the scene";
      return;
  }

  // Apply LUT and auto window/level
  vtkMRMLScalarVolumeDisplayNode* displayNode =
      vtkMRMLScalarVolumeDisplayNode::SafeDownCast(proxyVolume->GetDisplayNode());

  if (displayNode)
  {
      displayNode->SetAndObserveColorNodeID(petLUTNode->GetID());
      displayNode->AutoWindowLevelOn();
      displayNode->SetAutoWindowLevel(1);
      displayNode->SetVisibility(true);
  }

  vtkSlicerApplicationLogic* appLogic = qSlicerApplication::application()->applicationLogic();
  if (appLogic)
  {
      vtkMRMLSliceCompositeNode* compositeNode;

      // Iterate over all slice views
      for (int i = 0; i < appLogic->GetMRMLScene()->GetNumberOfNodesByClass("vtkMRMLSliceCompositeNode"); ++i)
      {
          compositeNode = vtkMRMLSliceCompositeNode::SafeDownCast(
              appLogic->GetMRMLScene()->GetNthNodeByClass(i, "vtkMRMLSliceCompositeNode"));
          if (compositeNode)
          {
              compositeNode->SetBackgroundVolumeID(proxyVolume->GetID());
          }
      }
  }
}

void qSlicerKMAPModuleWidget::onSegChanged (int index) {
  Q_D(qSlicerKMAPModuleWidget);
  this->segID = d->SegSelector->itemData(index).value<vtkIdType>();

  this->segSequenceNode = nullptr;

  vtkMRMLScene* scene = this->mrmlScene();
  if (scene==nullptr) {
    return;
  }
  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (!shNode) {
    return;
  }
  // Fetch PET
  if (this->petID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID) {
    d->populateNodeComboBox(d->SegSelector,
                            this->stuID,
                            "vtkMRMLSegmentationNode",
                            ""
                            );
    return;
  }
  vtkMRMLScalarVolumeNode* petNode = vtkMRMLScalarVolumeNode::SafeDownCast(shNode->GetItemDataNode(this->petID));
  if (!petNode || this->sequencePETNode == nullptr || this->sequenceBrowserPETNode == nullptr) {
    d->populateNodeComboBox(d->SegSelector,
                            this->stuID,
                            "vtkMRMLSegmentationNode",
                            ""
                            );
    return;
  }
  // Fetch Segmentation
  if (this->segID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID) {
    d->populateSegmentCheckboxes(this->segID);
    return;
  }
  vtkMRMLSegmentationNode* segNode = vtkMRMLSegmentationNode::SafeDownCast(shNode->GetItemDataNode(this->segID));
  if (!segNode) {
    d->populateSegmentCheckboxes(this->segID);
    return;
  }
  // Make sure the source of the segmentation is a binary label map, alongside a created closed surface
  vtkSlicerKMAPLogic* logic = vtkSlicerKMAPLogic::SafeDownCast(this->logic());
  if (!logic) {
    return;
  }
  logic->setupSeg(segNode);


  vtkMRMLSequenceNode* seqNode =
    this->sequenceBrowserPETNode->GetSequenceNode(segNode);
  if (!seqNode)
  {
    vtkSmartPointer<vtkMRMLSequenceNode> newSeqNode =
      vtkSmartPointer<vtkMRMLSequenceNode>::New();
    newSeqNode->SetName(shNode->GetItemName(segID).c_str());
    scene->AddNode(newSeqNode);
    this->sequenceBrowserPETNode->AddProxyNode(segNode, newSeqNode, false);
    this->sequenceBrowserPETNode->SetSaveChanges(newSeqNode, true);
    std::string indexValue;
    for (int i = 0; i < this->numberOfTimepoints; ++i)
    {
      indexValue = this->sequencePETNode->GetNthIndexValue(i);
      if (!newSeqNode->GetDataNodeAtValue(indexValue))
      {
        newSeqNode->SetDataNodeAtValue(segNode, indexValue);
      }
    }
    this->SegWatcher->ObserveSegmentationNode(segNode);
    seqNode = newSeqNode;
  }
  this->segSequenceNode = seqNode;

  d->populateSegmentCheckboxes(this->segID);
  this->enableTACbutton();
}

void qSlicerKMAPModuleWidget::onSegmentsChanged()
{
  Q_D(qSlicerKMAPModuleWidget);
  // std::vector<QString> selectedSegmentIDs;
  //
  // for (int i = 0; i < d->segmentCheckLayout->count(); ++i)
  // {
  //   QLayoutItem* item = d->segmentCheckLayout->itemAt(i);
  //   QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
  //   if (checkbox && checkbox->isChecked())
  //   {
  //     QString segmentID = checkbox->property("SegmentID").toString();
  //     selectedSegmentIDs.push_back(segmentID);
  //   }
  // }
  // this->segmentIDs = selectedSegmentIDs;
  d->populateSegmentCheckboxes(this->segID);
}

void qSlicerKMAPModuleWidget::clearTACdata() {
  Q_D(qSlicerKMAPModuleWidget);
  this->RemoveExistingPlotChartAndTable();
  this->segmentTACs.clear();
  this->segmentTACsnames.clear();
  d->populatePlotSegmentCheckboxes();
  d->populateIF();
  d->populateIFMTGA();
  d->populateIFImg();
  d->TACCollapsibleButton->setCollapsed(true);
  for (int i = 0; i < d->PlotStatsCheckLayout->count(); ++i)
  {
    QWidget* widget = d->PlotStatsCheckLayout->itemAt(i)->widget();
    QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
    if (cb)
    {
      cb->setChecked(false);
    }
  }
  d->direxcel->setCurrentPath(QString::fromStdString(""));
  d->saveExcelButton->setEnabled(false);
  return;
}

void qSlicerKMAPModuleWidget::enableTACbutton() {
  Q_D(qSlicerKMAPModuleWidget);
  if (this->ctID==vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID) {
    d->TACbutton->setEnabled(false);
    this->clearTACdata();
    return;
  }
  if (this->petID==vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID) {
    d->TACbutton->setEnabled(false);
    this->clearTACdata();
    return;
  }
  if (this->segID==vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID) {
    d->TACbutton->setEnabled(false);
    this->clearTACdata();
    return;
  }
  if (this->segmentIDs.empty()) {
    d->TACbutton->setEnabled(false);
    this->clearTACdata();
    return;
  }
  if (this->sequencePETNode==nullptr || this->sequenceBrowserPETNode==nullptr) {
    d->TACbutton->setEnabled(false);
    this->clearTACdata();
    return;
  }
  if (this->segSequenceNode==nullptr) {
    d->TACbutton->setEnabled(false);
    this->clearTACdata();
    return;
  }
  d->TACbutton->setEnabled(true);
}


vtkMRMLTableNode* qSlicerKMAPModuleWidget::GetOrCreatePlotTable()
{
  vtkMRMLTableNode* tableNode = vtkMRMLTableNode::SafeDownCast(
    this->mrmlScene()->GetFirstNodeByName("KMAP.PlotTable"));
  if (!tableNode)
  {
    tableNode = vtkMRMLTableNode::New();
    tableNode->SetName("KMAP.PlotTable");
    this->mrmlScene()->AddNode(tableNode);
    tableNode->Delete();
  }
  else
  {
    tableNode->RemoveAllColumns();
  }
  return tableNode;
}

vtkMRMLPlotChartNode* qSlicerKMAPModuleWidget::GetOrCreatePlotChart()
{
  vtkMRMLPlotChartNode* chartNode = vtkMRMLPlotChartNode::SafeDownCast(
    this->mrmlScene()->GetFirstNodeByName("KMAP.PlotChart"));
  if (!chartNode)
  {
    chartNode = vtkMRMLPlotChartNode::New();
    chartNode->SetName("KMAP.PlotChart");
    this->mrmlScene()->AddNode(chartNode);
    chartNode->Delete();
  }
  else
  {
    chartNode->RemoveAllPlotSeriesNodeIDs();
  }
  return chartNode;
}

void qSlicerKMAPModuleWidget::RemoveExistingPlotChartAndTable()
{
  this->MapPlotSeriesNodeIDToPlot.clear();
  this->ColNameToSegmentID.clear();
  this->PlotSelectedFrame = -1;
  this->PlotSelectedVOI.clear();
  this->lastSelection.clear();

  vtkMRMLPlotChartNode* chartNode = vtkMRMLPlotChartNode::SafeDownCast(
    this->mrmlScene()->GetFirstNodeByName("KMAP.PlotChart"));
  if (chartNode)
  {
    vtkCollection* viewNodes = this->mrmlScene()->GetNodesByClass("vtkMRMLPlotViewNode");
    if (viewNodes)
    {
      for (int i = 0; i < viewNodes->GetNumberOfItems(); ++i)
      {
        vtkMRMLPlotViewNode* viewNode = vtkMRMLPlotViewNode::SafeDownCast(
          viewNodes->GetItemAsObject(i));
        if (viewNode && viewNode->GetPlotChartNodeID() &&
            std::string(viewNode->GetPlotChartNodeID()) == chartNode->GetID())
        {
          viewNode->SetPlotChartNodeID(nullptr);
        }
      }
      viewNodes->Delete();
    }

    std::vector<std::string> seriesIDs;
    chartNode->GetPlotSeriesNodeIDs(seriesIDs);
    this->mrmlScene()->RemoveNode(chartNode);
    for (const std::string& id : seriesIDs)
    {
      vtkMRMLNode* node = this->mrmlScene()->GetNodeByID(id);
      if (node)
        this->mrmlScene()->RemoveNode(node);
    }
  }

  vtkMRMLTableNode* tableNode = vtkMRMLTableNode::SafeDownCast(
    this->mrmlScene()->GetFirstNodeByName("KMAP.PlotTable"));
  if (tableNode)
    this->mrmlScene()->RemoveNode(tableNode);
}

void qSlicerKMAPModuleWidget::onTACbutton()
{
  Q_D(qSlicerKMAPModuleWidget);
  vtkMRMLScene* scene = this->mrmlScene();
  if (!scene) {
    return;
  }
  // Run TAC computation
  vtkSlicerKMAPLogic* logic = vtkSlicerKMAPLogic::SafeDownCast(this->logic());
  if (!logic) {
    return;
  }

  if (this->durations.empty() || this->timePoints.empty()) {
    std::cerr << "Missing frame time information!" << std::endl;
    return;
  }

  std::vector<QString> segmentsToCompute;

  for (const QString& segmentID_qt : this->segmentIDs)
  {
    std::string segmentID = segmentID_qt.toStdString();
    if (this->segmentTACsnames.find(segmentID) == this->segmentTACsnames.end())
    {
      segmentsToCompute.push_back(segmentID_qt);
    }
    // else if (this->segmentHasChanged(segmentID))  // Implement this!
    // {
    //   // Segment changed — needs recomputing
    //   segmentsToCompute.push_back(segmentID_qt);
    // }
  }

  this->ProgressBar->setVisible(true);
  this->ProgressBar->show();
  // logic->computeTAC(this->ctID, this->petID, this->segID, segmentsToCompute, this->segmentTACs, this->segmentTACsnames, this->ProgressBar);
  this->stopRequested = false;
  if (this->PET_flatten_values.empty()) {
    logic->Image2Flatten(this->petID, this->PET_flatten_values, this->PETdims, this->numberOfTimepoints, this->ProgressBar, this->stopButton, this->stopRequested);
  }
  if (this->stopRequested) {
    return;
  }
  logic->TAC(this->sequencePETNode, this->segSequenceNode, segmentsToCompute, this->segmentTACs, this->segmentTACsnames, this->ProgressBar, this->stopButton, this->stopRequested);
  this->ProgressBar->setValue(0);
  this->ProgressBar->setVisible(false);
  if (this->stopRequested) {
    return;
  }
  // qApp->processEvents();
  d->populatePlotSegmentCheckboxes();
  d->populateIF();
  d->populateTimeBarMTGA();
  d->populateTimeBarMTGAImg();
  d->populateIFMTGA();
  d->populateIFImg();
  return;
}

void qSlicerKMAPModuleWidget::onSelectAllbutton()
{
  Q_D(qSlicerKMAPModuleWidget);
  d->SegmentCheckContents->blockSignals(true);
  if (this->segmentIDs.size()==(d->segmentCheckLayout->count()-1)) {
    for (int i = 0; i < d->segmentCheckLayout->count(); ++i)
    {
      QWidget* widget = d->segmentCheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->blockSignals(true);
        cb->setChecked(false);
        cb->blockSignals(false);
      }
    }
  } else {
    for (int i = 0; i < d->segmentCheckLayout->count(); ++i)
    {
      QWidget* widget = d->segmentCheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->blockSignals(true);
        cb->setChecked(true);
        cb->blockSignals(false);
      }
    }
  }
  d->SegmentCheckContents->blockSignals(false);
  d->populateSegmentCheckboxes(this->segID);
}

void qSlicerKMAPModuleWidget::onExcelPathChanged(const QString& path)
{
  Q_D(qSlicerKMAPModuleWidget);
  d->saveExcelButton->setEnabled(!path.trimmed().isEmpty());
}

void qSlicerKMAPModuleWidget::onExcelTCMPathChanged(const QString& path)
{
  Q_D(qSlicerKMAPModuleWidget);
  d->saveTCMExcelButton->setEnabled(!path.trimmed().isEmpty());
}

void qSlicerKMAPModuleWidget::onExcelMTGAPathChanged(const QString& path)
{
  Q_D(qSlicerKMAPModuleWidget);
  d->saveMTGAExcelButton->setEnabled(!path.trimmed().isEmpty());
}

void qSlicerKMAPModuleWidget::onExcelTCMfittedPathChanged(const QString& path)
{
  Q_D(qSlicerKMAPModuleWidget);
  d->saveTCMfittedExcelButton->setEnabled(!path.trimmed().isEmpty());
}

void qSlicerKMAPModuleWidget::onExcelMTGAfittedPathChanged(const QString& path)
{
  Q_D(qSlicerKMAPModuleWidget);
  d->saveMTGAfittedExcelButton->setEnabled(!path.trimmed().isEmpty());
}

QVariantMap qSlicerKMAPModuleWidget::TACtoPythonDict()
{
  Q_D(qSlicerKMAPModuleWidget);
  QVariantMap out;

  for (const auto& [segmentName, statsVec] : this->segmentTACs)
  {
    QVariantList voxelStatsList;
    const size_t N = statsVec.size();
    if (N != this->durations.size() || N != this->timePoints.size())
    {
      std::cerr << "Mismatch in vector sizes: statsVec (" << N
                << "), durations (" << durations.size()
                << "), timePoints (" << timePoints.size() << ")" << std::endl;
    }
    else
    {
      for (size_t i = 0; i < N; ++i)
      {
        const auto& vs = statsVec[i];
        QVariantMap vsMap;
        // Add time and duration
        vsMap["Time(s)"] = timePoints[i];
        vsMap["Duration"] = durations[i];
        // Add stats
        vsMap["VoxelCount"] = vs.count;
        vsMap["Mean"] = vs.mean;
        vsMap["Median"] = vs.median;
        vsMap["Min"] = vs.min;
        vsMap["Max"] = vs.max;
        vsMap["StDev"] = vs.stddev;
        vsMap["Q1"] = vs.q1;
        vsMap["Q3"] = vs.q3;
        vsMap["IQR"] = vs.iqr;
        vsMap["Volume(mm3)"] = vs.volume_mm3;
        vsMap["Volume(cm3)"] = vs.volume_cm3;

        voxelStatsList.append(vsMap);
      }
    }

    std :: string sheetName = this->segmentTACsnames[segmentName];
    if (sheetName.length() > 30)
    {
      sheetName = sheetName.substr(0, 30);
    }
    out[QString::fromStdString(sheetName)] = voxelStatsList;
  }

  return out;
}

QVariantMap qSlicerKMAPModuleWidget::TCMParamsToPythonDict()
{
  QVariantMap out;

  for (const auto& [segmentName, modelParamsMap] : this->segmentTCM)
  {
    QVariantList rows;

    for (const auto& [modelName, params] : modelParamsMap)
    {
      QVariantMap row;
      row["Model"] = QString::fromStdString(modelName);
      row["K1"]    = params.K1;
      row["k2"]    = params.k2;
      row["k3"]    = params.k3;
      row["k4"]    = params.k4;
      row["vb"]    = params.vb;
      row["td"]    = params.td;
      row["Ki"]    = params.Ki;
      row["DV"]    = params.DV;
      row["AIC"]   = params.AIC;
      row["BIC"]   = params.BIC;
      row["MASE"]  = params.MASE;
      row["chi^2_nu"]  = params.chi2;

      rows.append(row);
    }

    // Ensure Excel-friendly sheet name
    std::string sheetName = this->segmentTACsnames[segmentName];
    if (sheetName.length() > 30)
    {
      sheetName = sheetName.substr(0, 30);
    }

    out[QString::fromStdString(sheetName)] = rows;
  }

  return out;
}

QVariantMap qSlicerKMAPModuleWidget::MTGAParamsToPythonDict()
{
  QVariantMap out;

  for (const auto& [segmentName, modelParamsMap] : this->segmentMTGA)
  {
    QVariantList rows;

    for (const auto& [modelName, params] : modelParamsMap)
    {
      QVariantMap row;
      row["Model"] = QString::fromStdString(modelName);
      row["Ki"]    = params.Ki;
      row["DV"]    = params.DV;
      row["Intercept"] = params.Intercept;
      row["R2"]   = params.R2;
      row["AIC"]   = params.AIC;
      row["MASE"]  = params.MASE;
      // row["chi^2_nu"]  = params.chi2;

      rows.append(row);
    }

    // Ensure Excel-friendly sheet name
    std::string sheetName = this->segmentTACsnames[segmentName];
    if (sheetName.length() > 30)
    {
      sheetName = sheetName.substr(0, 30);
    }

    out[QString::fromStdString(sheetName)] = rows;
  }

  return out;
}

QVariantMap qSlicerKMAPModuleWidget::fittedTCMtoPythonDict()
{
  Q_D(qSlicerKMAPModuleWidget);
  QVariantMap out;

  int statIDQString = d->StatSelector->currentIndex();
  std::string currentSelectedStatID = d->StatSelector->itemData(statIDQString).toString().toStdString();

  for (const auto& [segmentName, tacvoi] : this->segmentTAC4TCMfits)
  {
    QVariantList rowList; // Will hold rows for this VOI
    const size_t N = this->timePoints.size();
    auto fitMapIt = this->segmentTCMfits.find(segmentName);
    if (fitMapIt == this->segmentTCMfits.end())
      continue;
    const auto& fits = fitMapIt->second;
    for (size_t i = 0; i < N; ++i)
    {
      QVariantMap row;
      row["Time(s)"] = this->timePoints[i];

      // Add TAC VOI
      if (!tacvoi.empty() && tacvoi.size() == N)
      {
        row[QString::fromStdString("TAC (" + currentSelectedStatID + ")")] = tacvoi[i][0];
      }

      // Add each TCM fit value for this time point
      for (const auto& [modelName, fitPtr] : fits)
      {
        if (fitPtr != nullptr)
        {
          row[QString::fromStdString(modelName)] = fitPtr[i];
        }
        else
        {
          row[QString::fromStdString(modelName)] = QVariant(); // blank cell
        }
      }
      rowList.append(row);
    }

    // Shorten sheet name for Excel if needed
    std::string sheetName = this->segmentTACsnames[segmentName];
    if (sheetName.length() > 30)
    {
      sheetName = sheetName.substr(0, 30);
    }

    out[QString::fromStdString(sheetName)] = rowList;
  }

  return out;
}

QVariantMap qSlicerKMAPModuleWidget::fittedMTGAtoPythonDict()
{
  Q_D(qSlicerKMAPModuleWidget);
  QVariantMap out;

  int statIDQString = d->StatSelectorMTGA->currentIndex();
  std::string currentSelectedStatID = d->StatSelectorMTGA->itemData(statIDQString).toString().toStdString();

  for (const auto& [segmentName, modelParamsMap] : this->segmentMTGA)
  {
    QVariantList rowList;

    if (!modelParamsMap.empty())
    {
        const auto& firstParams = modelParamsMap.begin()->second;
        const size_t N = firstParams.x.size();

        for (size_t i = 0; i < N; ++i)
        {
          QVariantMap row;

          for (const auto& [modelName, params] : modelParamsMap)
          {

            row[QString::fromStdString(modelName + "_x")] = params.x[i];
            row[QString::fromStdString(modelName + "_y")] = params.y[i];
            row[QString::fromStdString(modelName + "_fitted")] = params.fitted[i];

          }
          rowList.append(row);
        }
    }

    // Shorten sheet name for Excel if needed
    std::string sheetName = this->segmentTACsnames[segmentName];
    if (sheetName.length() > 30)
    {
      sheetName = sheetName.substr(0, 30);
    }

    out[QString::fromStdString(sheetName)] = rowList;

  }

  return out;
}

void qSlicerKMAPModuleWidget::onSaveExcelbutton()
{
  Q_D(qSlicerKMAPModuleWidget);
  QString path = d->direxcel->currentPath();
  QString filename = d->fileexcel->text();
  QString fullPath = QDir(path).filePath(filename);

  QVariantMap segmentDict = this->TACtoPythonDict();
  PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
  PythonQtObjectPtr result = mainContext.call("DPE_save_multisheet_excel", QVariantList{ fullPath, segmentDict });
}

void qSlicerKMAPModuleWidget::onSaveTCMExcelbutton()
{
  Q_D(qSlicerKMAPModuleWidget);
  QString path = d->direxceltcm->currentPath();
  QString filename = d->fileexceltcm->text();
  QString fullPath = QDir(path).filePath(filename);

  QVariantMap segmentDict = this->TCMParamsToPythonDict();
  PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
  PythonQtObjectPtr result = mainContext.call("DPE_saveTCM_multisheet_excel", QVariantList{ fullPath, segmentDict });
}

void qSlicerKMAPModuleWidget::onSaveMTGAExcelbutton()
{
  Q_D(qSlicerKMAPModuleWidget);
  QString path = d->direxcelmtga->currentPath();
  QString filename = d->fileexcelmtga->text();
  QString fullPath = QDir(path).filePath(filename);

  QVariantMap segmentDict = this->MTGAParamsToPythonDict();
  PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
  PythonQtObjectPtr result = mainContext.call("DPE_saveMTGA_multisheet_excel", QVariantList{ fullPath, segmentDict });
}

void qSlicerKMAPModuleWidget::onSaveTCMfittedExcelbutton()
{
  Q_D(qSlicerKMAPModuleWidget);
  QString path = d->direxceltcmfitted->currentPath();
  QString filename = d->fileexceltcmfitted->text();
  QString fullPath = QDir(path).filePath(filename);

  QVariantMap segmentDict = this->fittedTCMtoPythonDict();
  PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
  PythonQtObjectPtr result = mainContext.call("DPE_generic_save_multisheet_excel", QVariantList{ fullPath, segmentDict });
}

void qSlicerKMAPModuleWidget::onSaveMTGAfittedExcelbutton()
{
  Q_D(qSlicerKMAPModuleWidget);
  QString path = d->direxcelmtgafitted->currentPath();
  QString filename = d->fileexcelmtgafitted->text();
  QString fullPath = QDir(path).filePath(filename);

  QVariantMap segmentDict = this->fittedMTGAtoPythonDict();
  PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
  PythonQtObjectPtr result = mainContext.call("DPE_genericMTGA_save_multisheet_excel", QVariantList{ fullPath, segmentDict });
}

void qSlicerKMAPModuleWidget::onSelectedPoint(vtkStringArray* mrmlPlotSeriesIDs, vtkCollection* selectionCol)
{
  if (!this->checkdisplayedKMAP())
    return;

  vtkMRMLScene* scene = this->mrmlScene();
  if (scene==nullptr)
    return;

  if (!mrmlPlotSeriesIDs || !selectionCol || this->MapPlotSeriesNodeIDToPlot.empty())
    return;

  QSet<QPair<QString, vtkIdType>> newSelection;

  // Loop over each series that has selected points
  int psf_value = -1;
  std :: string psv_value = "";
  std :: string lastseriesID = "";
  for (vtkIdType i = 0; i < mrmlPlotSeriesIDs->GetNumberOfValues(); ++i)
  {
    QString seriesID = QString::fromStdString(mrmlPlotSeriesIDs->GetValue(i));
    vtkIdTypeArray* selectedPoints = vtkIdTypeArray::SafeDownCast(
        selectionCol->GetItemAsObject(i));

    if (!selectedPoints)
      continue;

    vtkIdType pointIndex = selectedPoints->GetValue(selectedPoints->GetNumberOfTuples()-1);
    QPair<QString, vtkIdType> candidate(seriesID, pointIndex);
    vtkMRMLPlotSeriesNode* seriesNode = vtkMRMLPlotSeriesNode::SafeDownCast(
        scene->GetNodeByID(seriesID.toStdString()));
    if (!seriesNode)
        return;
    vtkMRMLTableNode* tableNode = seriesNode->GetTableNode();
    if (!tableNode)
        return;
    vtkTable* table = tableNode->GetTable();
    if (!table)
        return;
    std::string labelColName = seriesNode->GetLabelColumnName();
    vtkAbstractArray* labelArray = table->GetColumnByName(labelColName.c_str());
    if (!labelArray)
        return;

    vtkStringArray* strArray = vtkStringArray::SafeDownCast(labelArray);
    if (!strArray)
        return;

    QString labelValue = QString::fromStdString(strArray->GetValue(pointIndex));
    QStringList parts = labelValue.split(',');
    if (!parts.isEmpty())
    {
      QString framePart = parts[0].trimmed();
      framePart.remove("Frame:");
      psf_value = framePart.trimmed().toInt();
      psv_value = this->ColNameToSegmentID[seriesNode->GetName()];
    }
    if (!newSelection.contains(candidate))
      newSelection.insert(candidate);
    if (!this->lastSelection.contains(candidate))
    {
      this->PlotSelectedFrame = psf_value;
      this->PlotSelectedVOI = psv_value;
    } else {
      if (this->PlotSelectedFrame == psf_value && this->PlotSelectedVOI == psv_value) {
        lastseriesID = seriesID.toStdString();
        continue;
      }
      if (!this->MapPlotSeriesNodeIDToPlot.contains(seriesID)) {
        vtkGenericWarningMacro("MapPlotSeriesNodeIDToPlot does not contain seriesID=" << seriesID.toStdString());
      }
      vtkPlot* vtkplot = this->MapPlotSeriesNodeIDToPlot.value(seriesID);
      vtkSmartPointer<vtkIdTypeArray> emptySelection = vtkSmartPointer<vtkIdTypeArray>::New();
      vtkplot->SetSelection(emptySelection);
    }
  }
  if (newSelection.size()==0) {
    // qDebug() << "No point selected";
    this->PlotSelectedFrame = -1;
    this->PlotSelectedVOI = "";
  } else if (newSelection.size()==1 && lastseriesID!="") {
    // qDebug() << "Old point selected";
  } else {
    // qDebug() << "New point selected";
    if (!this->MapPlotSeriesNodeIDToPlot.contains(QString::fromStdString(lastseriesID))) {
        vtkGenericWarningMacro("MapPlotSeriesNodeIDToPlot does not contain seriesID=" << lastseriesID);
    }
    vtkPlot* vtkplot = this->MapPlotSeriesNodeIDToPlot.value(QString::fromStdString(lastseriesID));
    vtkSmartPointer<vtkIdTypeArray> emptySelection = vtkSmartPointer<vtkIdTypeArray>::New();
    vtkplot->SetSelection(emptySelection);
  }
  this->lastSelection = newSelection;
}

void qSlicerKMAPModuleWidget::onPlotbutton()
{
  Q_D(qSlicerKMAPModuleWidget);
  this->ColNameToSegmentID.clear();
  this->MapPlotSeriesNodeIDToPlot.clear();
  vtkMRMLScene* scene = this->mrmlScene();
  // Get selected segments
  std::vector<std::string> PlotSelectedIDs;
  for (int i = 0; i < d->PlotsegmentCheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->PlotsegmentCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      std::string segmentid = checkbox->property("SegmentID").toString().toStdString();
      PlotSelectedIDs.push_back(segmentid);
    }
  }

  // Get selected stats
  std::vector<std::string> PlotSelectedStats;
  for (int i = 0; i < d->PlotStatsCheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->PlotStatsCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      PlotSelectedStats.push_back(checkbox->text().toStdString());
    }
  }

  if (PlotSelectedIDs.empty() || PlotSelectedStats.empty())
    return;

  // Clear previous plot/chart/table
  this->RemoveExistingPlotChartAndTable();

  // Create or get table
  vtkSmartPointer<vtkMRMLTableNode> tableNode = this->GetOrCreatePlotTable();

  // Add time column
  vtkNew<vtkDoubleArray> timeArray;
  timeArray->SetName("Time (min)");
  vtkNew<vtkStringArray> labelArray;
  labelArray->SetName("ToolTipLabelTAC");
  for (int i=0; i < this->timePoints.size(); ++i) {
    double t = this->timePoints[i];
    timeArray->InsertNextValue(t/60.);
    std::ostringstream oss;
    oss << "Frame: " << i
        << ", Time(s): " << this->timePoints[i]
        << ", Time(min): " << this->timePoints[i]/60.0;
    labelArray->InsertNextValue(oss.str());
  }
  tableNode->AddColumn(timeArray);
  tableNode->AddColumn(labelArray);

  // Create plot chart
  vtkMRMLPlotChartNode* chartNode = this->GetOrCreatePlotChart();
  chartNode->SetTitle("Time Activity Curve");
  chartNode->SetXAxisTitle("Time (min)");
  chartNode->SetYAxisTitle("SUVbw (g/mL)");
  std::unordered_map<std::string, std::string> LabelToSeriesID;
  for (const std::string& segmentID : PlotSelectedIDs)
  {
    std :: string segmentName = this->segmentTACsnames[segmentID];
    for (const std::string& statName : PlotSelectedStats)
    {
      std::string colName = segmentName + " - " + statName;
      this->ColNameToSegmentID[colName] = segmentID;
      vtkNew<vtkDoubleArray> statArray;
      statArray->SetName(colName.c_str());

      vtkNew<vtkDoubleArray> statErrArray;
      std::string colErrName = colName + " Error";
      statErrArray->SetName(colErrName.c_str());

      vtkNew<vtkDoubleArray> statArrayLine;
      std::string statArrayLineName = colName + " Line";
      statArrayLine->SetName(statArrayLineName.c_str());

      for (int ivs=0; ivs<this->segmentTACs[segmentID].size(); ++ivs)
      {
        const VoxelStatistics& vs = this->segmentTACs[segmentID][ivs];
        double value = std::numeric_limits<double>::quiet_NaN();
        if (vs.keep) {
          if (statName == "Mean")
          {
            value = vs.mean;
            // if (d->PlotErrorCheckbox && d->PlotErrorCheckbox->isChecked())
            //   statErrArray->InsertNextValue(vs.stddev);
          }
          else if (statName == "Median")
          {
            value = vs.median;
            // if (d->PlotErrorCheckbox && d->PlotErrorCheckbox->isChecked())
            //   statErrArray->InsertNextValue(vs.iqr);
          }
          else if (statName == "VoxelCount") value = vs.count;
          else if (statName == "Min")        value = vs.min;
          else if (statName == "Max")        value = vs.max;
          else if (statName == "Volume(cc)") value = vs.volume_cm3;
          else vtkGenericWarningMacro("Unknown stat name: " << statName);
        }

        statArray->InsertNextValue(value);
        // Line points
        if (!std::isnan(value)) {
          statArrayLine->InsertNextValue(value);
        } else {
          double nextValue = std::numeric_limits<double>::quiet_NaN();
          double x1 = std::numeric_limits<double>::quiet_NaN();
          for (int next_ivs = ivs+1; next_ivs<this->timePoints.size(); ++next_ivs) {
            const VoxelStatistics& vsNext = this->segmentTACs[segmentID][next_ivs];
            if (vsNext.keep)
            {
                x1 = this->timePoints[next_ivs];
                if (statName == "Mean") nextValue = vsNext.mean;
                else if (statName == "Median") nextValue = vsNext.median;
                else if (statName == "VoxelCount") nextValue = vsNext.count;
                else if (statName == "Min") nextValue = vsNext.min;
                else if (statName == "Max") nextValue = vsNext.max;
                else if (statName == "Volume(cc)") nextValue = vsNext.volume_cm3;
                break;
            }
          }
          if (std::isnan(nextValue)) {
            statArrayLine->InsertNextValue(std::numeric_limits<double>::quiet_NaN());
            continue;
          }
          double prevValue = std::numeric_limits<double>::quiet_NaN();
          double x0 = std::numeric_limits<double>::quiet_NaN();
          for (int prev_ivs = ivs-1; prev_ivs>=0; --prev_ivs) {
            const VoxelStatistics& vsPrev = this->segmentTACs[segmentID][prev_ivs];
            if (vsPrev.keep)
            {
                x0 = this->timePoints[prev_ivs];
                if (statName == "Mean") prevValue = vsPrev.mean;
                else if (statName == "Median") prevValue = vsPrev.median;
                else if (statName == "VoxelCount") prevValue = vsPrev.count;
                else if (statName == "Min") prevValue = vsPrev.min;
                else if (statName == "Max") prevValue = vsPrev.max;
                else if (statName == "Volume(cc)") prevValue = vsPrev.volume_cm3;
                break;
            }
          }
          if (std::isnan(prevValue)) {
            statArrayLine->InsertNextValue(std::numeric_limits<double>::quiet_NaN());
            continue;
          }

          double x  = this->timePoints[ivs];
          // Proper linear interpolation
          value = prevValue + ((x - x0) / (x1 - x0)) * (nextValue - prevValue);
          statArrayLine->InsertNextValue(value);
        }

      }

      tableNode->AddColumn(statArray);
      tableNode->AddColumn(statArrayLine);
      // if (statErrArray->GetNumberOfTuples() > 0)
      //   tableNode->AddColumn(statErrArray);

      vtkSmartPointer<vtkMRMLPlotSeriesNode> lineSeries = vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
      scene->AddNode(lineSeries);
      lineSeries->SetName("");
      lineSeries->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);
      lineSeries->SetAndObserveTableNodeID(tableNode->GetID());
      lineSeries->SetXColumnName("Time (min)");
      lineSeries->SetYColumnName(statArrayLineName.c_str());
      lineSeries->SetLabelColumnName("ToolTipLabelTAC");
      lineSeries->SetUniqueColor();
      lineSeries->SetMarkerStyle(vtkMRMLPlotSeriesNode::MarkerStyleNone);
      chartNode->AddAndObservePlotSeriesNodeID(lineSeries->GetID());

      vtkSmartPointer<vtkMRMLPlotSeriesNode> series = vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
      scene->AddNode(series);
      series->SetName(colName.c_str());
      series->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);
      series->SetAndObserveTableNodeID(tableNode->GetID());
      series->SetXColumnName("Time (min)");
      series->SetYColumnName(colName.c_str());
      series->SetLabelColumnName("ToolTipLabelTAC");
      series->SetLineStyle(vtkMRMLPlotSeriesNode::LineStyleNone);
      series->SetColor(lineSeries->GetColor());
      chartNode->AddAndObservePlotSeriesNodeID(series->GetID());
      LabelToSeriesID[colName] = series->GetID();
    }
  }

  // Show plot view
  auto* layoutNode = vtkMRMLLayoutNode::SafeDownCast(scene->GetFirstNodeByClass("vtkMRMLLayoutNode"));
  if (layoutNode)
    layoutNode->SetViewArrangement(vtkMRMLLayoutNode::SlicerLayoutConventionalPlotView);

  vtkMRMLPlotViewNode* plotViewNode = vtkMRMLPlotViewNode::SafeDownCast(
    scene->GetFirstNodeByClass("vtkMRMLPlotViewNode"));
  if (plotViewNode)
  {
    plotViewNode->SetPlotChartNodeID(chartNode->GetID());
    qMRMLPlotWidget* plotWidget = nullptr;
    if (qSlicerApplication::application())
    {
      qSlicerLayoutManager* layoutManager =
          qSlicerApplication::application()->layoutManager();
      qMRMLPlotWidget* plotWidget = nullptr;
      plotWidget = layoutManager->plotWidget(0);
      qMRMLPlotView* plotView = plotWidget->plotView();
      if (plotView)
      {
        QObject::connect(plotView, SIGNAL(dataSelected(vtkStringArray*, vtkCollection*)),
                         this, SLOT(onSelectedPoint(vtkStringArray*, vtkCollection*)));

        vtkSmartPointer<vtkChartXY> chart = plotView->chart();
        for (int i = 0; i < chart->GetNumberOfPlots(); ++i)
        {
           vtkPlot* plot = chart->GetPlot(i);
           std :: string PlotLabel = plot->GetLabel();
           QString seriesNodeID = QString::fromStdString(LabelToSeriesID[PlotLabel]);
           this->MapPlotSeriesNodeIDToPlot[seriesNodeID] = plot;
        }
      }
    }
  }
}

void qSlicerKMAPModuleWidget::onIFSelectionChanged(int index)
{
  Q_D(qSlicerKMAPModuleWidget);
  this->IFID = d->IFSelector->itemData(index).toString().toStdString();
  d->populateVOI(this->IFID);
  d->IFSelectorMTGA->setCurrentIndex(index);
  d->populateVOIMTGA(this->IFID);
  d->IFSelectorImg->setCurrentIndex(index);
  if (this->IFID == "")
  {
    return;
  }

}

void qSlicerKMAPModuleWidget::onIFMTGASelectionChanged(int index)
{
  Q_D(qSlicerKMAPModuleWidget);
  this->IFID = d->IFSelectorMTGA->itemData(index).toString().toStdString();
  d->populateVOIMTGA(this->IFID);
  d->IFSelector->setCurrentIndex(index);
  d->populateVOI(this->IFID);
  d->IFSelectorImg->setCurrentIndex(index);
  if (this->IFID == "")
  {
    return;
  }

}

void qSlicerKMAPModuleWidget::onIFImgSelectionChanged(int index)
{
  Q_D(qSlicerKMAPModuleWidget);
  this->IFID = d->IFSelectorImg->itemData(index).toString().toStdString();
  d->IFSelector->setCurrentIndex(index);
  d->IFSelectorMTGA->setCurrentIndex(index);
  d->populateVOI(this->IFID);
  d->populateVOIMTGA(this->IFID);
  if (this->IFID == "")
  {
    return;
  }

}

void qSlicerKMAPModuleWidget::onVOISelectionChanged(int index)
{
  Q_D(qSlicerKMAPModuleWidget);
  std :: string segmentID = d->VOISelector->itemData(index).toString().toStdString();
  this->plotTCMVOI=segmentID;
  d->populateResultsTable(segmentID);
  if (segmentID == "")
  {
    return;
  }
}

void qSlicerKMAPModuleWidget::onVOIMTGASelectionChanged(int index)
{
  Q_D(qSlicerKMAPModuleWidget);
  std :: string segmentID = d->VOISelectorMTGA->itemData(index).toString().toStdString();
  this->plotMTGAVOI=segmentID;
  d->populateResultsMTGATable(segmentID);
  if (segmentID == "")
  {
    return;
  }
}

void qSlicerKMAPModuleWidget::onVOISelectAllbutton()
{
  Q_D(qSlicerKMAPModuleWidget);
  d->VOICheckContents->blockSignals(true);
  if (this->VOIsegmentIDs.size()==(d->VOICheckLayout->count()-1)) {
    for (int i = 0; i < d->VOICheckLayout->count(); ++i)
    {
      QWidget* widget = d->VOICheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->blockSignals(true);
        cb->setChecked(false);
        cb->blockSignals(false);
      }
    }
  } else {
    for (int i = 0; i < d->VOICheckLayout->count(); ++i)
    {
      QWidget* widget = d->VOICheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->blockSignals(true);
        cb->setChecked(true);
        cb->blockSignals(false);
      }
    }
  }
  d->VOICheckContents->blockSignals(false);
  d->populateVOI(this->IFID);
}

void qSlicerKMAPModuleWidget::onVOIMTGASelectAllbutton()
{
  Q_D(qSlicerKMAPModuleWidget);
  d->VOIMTGACheckContents->blockSignals(true);
  if (this->VOIMTGAsegmentIDs.size()==(d->VOIMTGACheckLayout->count()-1)) {
    for (int i = 0; i < d->VOIMTGACheckLayout->count(); ++i)
    {
      QWidget* widget = d->VOIMTGACheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->blockSignals(true);
        cb->setChecked(false);
        cb->blockSignals(false);
      }
    }
  } else {
    for (int i = 0; i < d->VOIMTGACheckLayout->count(); ++i)
    {
      QWidget* widget = d->VOIMTGACheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->blockSignals(true);
        cb->setChecked(true);
        cb->blockSignals(false);
      }
    }
  }
  d->VOIMTGACheckContents->blockSignals(false);
  d->populateVOIMTGA(this->IFID);
}

void qSlicerKMAPModuleWidget::onOLSclicked()
{
  Q_D(qSlicerKMAPModuleWidget);
  d->weightedFitCheckBox->blockSignals(true);
  d->robustFitCheckBox->blockSignals(true);
  d->weightedFitCheckBox->setChecked(false);
  d->robustFitCheckBox->setChecked(false);
  d->weightedFitCheckBox->blockSignals(false);
  d->robustFitCheckBox->blockSignals(false);
  d->robustParamsWidget->setVisible(false);
  if (!d->olsFitCheckBox->isChecked())
    d->olsFitCheckBox->setChecked(true);
}

void qSlicerKMAPModuleWidget::onOLSImgclicked()
{
  Q_D(qSlicerKMAPModuleWidget);
  d->weightedFitCheckBoxImg->blockSignals(true);
  d->robustFitCheckBoxImg->blockSignals(true);
  d->weightedFitCheckBoxImg->setChecked(false);
  d->robustFitCheckBoxImg->setChecked(false);
  d->weightedFitCheckBoxImg->blockSignals(false);
  d->robustFitCheckBoxImg->blockSignals(false);
  d->robustParamsWidgetImg->setVisible(false);
  if (!d->olsFitCheckBoxImg->isChecked())
    d->olsFitCheckBoxImg->setChecked(true);
}

void qSlicerKMAPModuleWidget::onWLSclicked()
{
  Q_D(qSlicerKMAPModuleWidget);
  d->olsFitCheckBox->blockSignals(true);
  d->robustFitCheckBox->blockSignals(true);
  d->olsFitCheckBox->setChecked(false);
  d->robustFitCheckBox->setChecked(false);
  d->olsFitCheckBox->blockSignals(false);
  d->robustFitCheckBox->blockSignals(false);
  d->robustParamsWidget->setVisible(false);
  if (!d->weightedFitCheckBox->isChecked())
    d->olsFitCheckBox->setChecked(true);
}

void qSlicerKMAPModuleWidget::onWLSImgclicked()
{
  Q_D(qSlicerKMAPModuleWidget);
  d->olsFitCheckBoxImg->blockSignals(true);
  d->robustFitCheckBoxImg->blockSignals(true);
  d->olsFitCheckBoxImg->setChecked(false);
  d->robustFitCheckBoxImg->setChecked(false);
  d->olsFitCheckBoxImg->blockSignals(false);
  d->robustFitCheckBoxImg->blockSignals(false);
  d->robustParamsWidgetImg->setVisible(false);
  if (!d->weightedFitCheckBoxImg->isChecked())
    d->olsFitCheckBoxImg->setChecked(true);
}

void qSlicerKMAPModuleWidget::onRLSclicked()
{
  Q_D(qSlicerKMAPModuleWidget);
  d->olsFitCheckBox->blockSignals(true);
  d->weightedFitCheckBox->blockSignals(true);
  d->olsFitCheckBox->setChecked(false);
  d->weightedFitCheckBox->setChecked(false);
  d->olsFitCheckBox->blockSignals(false);
  d->weightedFitCheckBox->blockSignals(false);
  if (!d->robustFitCheckBox->isChecked()) {
    d->olsFitCheckBox->setChecked(true);
    d->robustParamsWidget->setVisible(false);
  } else {
    d->robustParamsWidget->setVisible(true);
  }
}

void qSlicerKMAPModuleWidget::onRLSImgclicked()
{
  Q_D(qSlicerKMAPModuleWidget);
  d->olsFitCheckBoxImg->blockSignals(true);
  d->weightedFitCheckBoxImg->blockSignals(true);
  d->olsFitCheckBoxImg->setChecked(false);
  d->weightedFitCheckBoxImg->setChecked(false);
  d->olsFitCheckBoxImg->blockSignals(false);
  d->weightedFitCheckBoxImg->blockSignals(false);
  if (!d->robustFitCheckBoxImg->isChecked()) {
    d->olsFitCheckBoxImg->setChecked(true);
    d->robustParamsWidgetImg->setVisible(false);
  } else {
    d->robustParamsWidgetImg->setVisible(true);
  }
}

void qSlicerKMAPModuleWidget::onStdFitclicked()
{
  Q_D(qSlicerKMAPModuleWidget);
  d->weightFitCheckBox->blockSignals(true);
  d->weightFitCheckBox->setChecked(false);
  d->weightFitCheckBox->blockSignals(false);
  if (!d->standardFitCheckBox->isChecked()) {
    d->standardFitCheckBox->setChecked(true);
  }
}

void qSlicerKMAPModuleWidget::onStdFitImgclicked()
{
  Q_D(qSlicerKMAPModuleWidget);
  d->weightFitCheckBoxImg->blockSignals(true);
  d->weightFitCheckBoxImg->setChecked(false);
  d->weightFitCheckBoxImg->blockSignals(false);
  if (!d->standardFitCheckBoxImg->isChecked()) {
    d->standardFitCheckBoxImg->setChecked(true);
  }
}

void qSlicerKMAPModuleWidget::onWFitclicked()
{
  Q_D(qSlicerKMAPModuleWidget);
  d->standardFitCheckBox->blockSignals(true);
  d->standardFitCheckBox->setChecked(false);
  d->standardFitCheckBox->blockSignals(false);
  if (!d->weightFitCheckBox->isChecked()) {
    d->standardFitCheckBox->setChecked(true);
  }
}

void qSlicerKMAPModuleWidget::onWFitImgclicked()
{
  Q_D(qSlicerKMAPModuleWidget);
  d->standardFitCheckBoxImg->blockSignals(true);
  d->standardFitCheckBoxImg->setChecked(false);
  d->standardFitCheckBoxImg->blockSignals(false);
  if (!d->weightFitCheckBoxImg->isChecked()) {
    d->standardFitCheckBoxImg->setChecked(true);
  }
}

void qSlicerKMAPModuleWidget::onSliderChanged(int index)
{
  Q_D(qSlicerKMAPModuleWidget);
  double timeSec = this->timePoints[index-1];
  double timeMin = timeSec / 60.0;

  d->frameEdit->setText(QString::number(index));
  d->timeSecEdit->setText(QString::number(timeSec, 'f', 2));
  d->timeMinEdit->setText(QString::number(timeMin, 'f', 2));
}

void qSlicerKMAPModuleWidget::onSliderImgChanged(int index)
{
  Q_D(qSlicerKMAPModuleWidget);
  double timeSec = this->timePoints[index-1];
  double timeMin = timeSec / 60.0;

  d->frameEditImg->setText(QString::number(index));
  d->timeSecEditImg->setText(QString::number(timeSec, 'f', 2));
  d->timeMinEditImg->setText(QString::number(timeMin, 'f', 2));
}

void qSlicerKMAPModuleWidget::runVuong(std::string sel1,
                                       std::string sel2,
                                       std::string segmentID
                                     )
{
  Q_D(qSlicerKMAPModuleWidget);
  vtkSlicerKMAPLogic* logic = vtkSlicerKMAPLogic::SafeDownCast(this->logic());
  if (!logic) {
    std::cerr << "Missing Logic!" << std::endl;
    return;
  }

  if (segmentID.empty()) {
    d->MTGAModel1->clear();
    d->MTGAModel2->clear();
    d->MTGAVuongP->setText("");
    return;
  }

  auto it = this->segmentMTGA.find(segmentID);
  if (it == this->segmentMTGA.end()) {
    d->MTGAModel1->clear();
    d->MTGAModel2->clear();
    d->MTGAVuongP->setText("");
    return;
  }
  auto& modelsForSegment = it->second;

  const int N1 = modelsForSegment[sel1].y.size();
  const int N2 = modelsForSegment[sel2].y.size();
  if (N1 != N2) {
    throw std::runtime_error(
        sel1 + " has not been fitted with the same number of datapoints (" + std::to_string(N1) +
        ") of " + sel2 + " (" + std::to_string(N2) + ")."
    );
  }
  std::vector<double> w1 = modelsForSegment[sel1].weights;
  std::vector<double> w2 = modelsForSegment[sel2].weights;
  // Check they are the same length
  if (w1.size() != w2.size()) {
      throw std::invalid_argument("Weight vectors must have the same length");
  }
  // Compute average weights
  std::vector<double> wgt_avg(w1.size());
  for (size_t i = 0; i < w1.size(); ++i) {
      wgt_avg[i] = 0.5 * (w1[i] + w2[i]);
  }
  const std::vector<double>* wgt = &wgt_avg;
  double p = logic->computeVuongP(modelsForSegment[sel1].r,
                                  modelsForSegment[sel2].r,
                                  wgt,
                                  modelsForSegment[sel1].dof,
                                  modelsForSegment[sel2].dof,
                                  VuongCorrection::BIC,
                                  Tail::TwoSided
                                );
  d->MTGAVuongP->setText(QString::number(p, 'g', 4));
  // d->MTGAVuongP->adjustSize();
  return;
}

void qSlicerKMAPModuleWidget::runTCMstat(std::string sel1,
                                         std::string sel2,
                                         std::string segmentID
                                        )
{
  Q_D(qSlicerKMAPModuleWidget);
  vtkSlicerKMAPLogic* logic = vtkSlicerKMAPLogic::SafeDownCast(this->logic());
  if (!logic) {
    std::cerr << "Missing Logic!" << std::endl;
    return;
  }

  if (segmentID.empty()) {
    d->TCMModel1->clear();
    d->TCMModel2->clear();
    d->TCMLRTP->setText("");
    d->TCMVuongP->setText("");
    return;
  }

  auto it = this->segmentTCM.find(segmentID);
  if (it == this->segmentTCM.end()) {
    d->TCMModel1->clear();
    d->TCMModel2->clear();
    d->TCMLRTP->setText("");
    d->TCMVuongP->setText("");
    return;
  }
  auto& modelsForSegment = it->second;

  auto itTac = this->segmentTAC4TCMfits.find(segmentID);
  if (itTac == this->segmentTAC4TCMfits.end()) {
    d->TCMModel1->clear();
    d->TCMModel2->clear();
    d->TCMLRTP->setText("");
    d->TCMVuongP->setText("");
    return;
  }

  auto itFits = this->segmentTCMfits.find(segmentID);
  if (itFits == this->segmentTCMfits.end()) {
    d->TCMModel1->clear();
    d->TCMModel2->clear();
    d->TCMLRTP->setText("");
    d->TCMVuongP->setText("");
    return;
  }


  const int N1 = modelsForSegment[sel1].weights.size();
  const int N2 = modelsForSegment[sel2].weights.size();
  if (N1 != N2) {
    throw std::runtime_error(
        sel1 + " has not been fitted with the same number of datapoints (" + std::to_string(N1) +
        ") of " + sel2 + " (" + std::to_string(N2) + ")."
    );
  }
  std::vector<double> w1 = modelsForSegment[sel1].weights;
  std::vector<double> w2 = modelsForSegment[sel2].weights;
  // Check they are the same length
  if (w1.size() != w2.size()) {
      throw std::invalid_argument("Weight vectors must have the same length");
  }
  // Compute average weights
  std::vector<double> wgt_avg(w1.size());
  for (size_t i = 0; i < w1.size(); ++i) {
      wgt_avg[i] = 0.5 * (w1[i] + w2[i]);
  }
  const std::vector<double>* wgt = &wgt_avg;
  double p = logic->computeLRTP(modelsForSegment[sel1].loglik,
                                modelsForSegment[sel2].loglik,
                                modelsForSegment[sel1].dof,
                                modelsForSegment[sel2].dof
                                );
  d->TCMLRTP->setText(QString::number(p, 'g', 4));
  double pvuong = logic->computeVuongP(modelsForSegment[sel1].r,
                                       modelsForSegment[sel2].r,
                                       wgt,
                                       modelsForSegment[sel1].dof,
                                       modelsForSegment[sel2].dof,
                                       VuongCorrection::BIC,
                                       Tail::TwoSided
                                     );
  d->TCMVuongP->setText(QString::number(pvuong, 'g', 4));
  return;
}

void qSlicerKMAPModuleWidget::onMTGAModelBox(int index)
{
  Q_D(qSlicerKMAPModuleWidget);

  std::string selectedVOI = this->plotMTGAVOI;
  if (selectedVOI.empty()) {
    d->MTGAModel1->clear();
    d->MTGAModel2->clear();
    d->MTGAVuongP->setText("");
    return;
  }

  auto it = this->segmentMTGA.find(selectedVOI);
  if (it == this->segmentMTGA.end()) {
    d->MTGAModel1->clear();
    d->MTGAModel2->clear();
    d->MTGAVuongP->setText("");
    return;
  }
  const auto& modelsForSegment = it->second;
  std::string sel1, sel2;
  int idx1 = d->MTGAModel1->currentIndex();
  if (idx1 >= 0)
    sel1 = d->MTGAModel1->itemData(idx1).toString().toStdString();
  int idx2 = d->MTGAModel2->currentIndex();
  if (idx2 >= 0)
    sel2 = d->MTGAModel2->itemData(idx2).toString().toStdString();
  d->populateModelCombo(d->MTGAModel1, sel2, sel1, selectedVOI);
  d->populateModelCombo(d->MTGAModel2, sel1, sel2, selectedVOI);
  if (idx1>0 & idx2>0) {
    this->runVuong(sel1, sel2, selectedVOI);
  } else {
    d->MTGAVuongP->setText("");
  }
  return;
}

void qSlicerKMAPModuleWidget::onTCMModelBox(int index)
{
  Q_D(qSlicerKMAPModuleWidget);

  std::string selectedVOI = this->plotTCMVOI;
  if (selectedVOI.empty()) {
    d->TCMModel1->clear();
    d->TCMModel2->clear();
    d->TCMLRTP->setText("");
    d->TCMVuongP->setText("");
    return;
  }

  auto it = this->segmentTCM.find(selectedVOI);
  if (it == this->segmentTCM.end()) {
    d->TCMModel1->clear();
    d->TCMModel2->clear();
    d->TCMLRTP->setText("");
    d->TCMVuongP->setText("");
    return;
  }
  const auto& modelsForSegment = it->second;
  std::string sel1, sel2;
  int idx1 = d->TCMModel1->currentIndex();
  if (idx1 >= 0)
    sel1 = d->TCMModel1->itemData(idx1).toString().toStdString();
  int idx2 = d->TCMModel2->currentIndex();
  if (idx2 >= 0)
    sel2 = d->TCMModel2->itemData(idx2).toString().toStdString();
  d->populateModelComboTCM(d->TCMModel1, sel2, sel1, selectedVOI);
  d->populateModelComboTCM(d->TCMModel2, sel1, sel2, selectedVOI);
  if (idx1>0 & idx2>0) {
    this->runTCMstat(sel1, sel2, selectedVOI);
  } else {
    d->TCMLRTP->setText("");
    d->TCMVuongP->setText("");
  }
  return;
}

void qSlicerKMAPModuleWidget::clearFITdata() {
  Q_D(qSlicerKMAPModuleWidget);
  this->segmentTCM.clear();
  d->populateResultsVOI();
  d->TCMResultsButton->setCollapsed(true);
  return;

}

void qSlicerKMAPModuleWidget::clearFITMTGAdata() {
  Q_D(qSlicerKMAPModuleWidget);
  this->segmentMTGA.clear();
  d->populateResultsVOIMTGA();
  d->MTGAResultsButton->setCollapsed(true);
  return;
}



void qSlicerKMAPModuleWidget::enableFITbutton() {
  Q_D(qSlicerKMAPModuleWidget);
  if (this->IFID=="") {
    d->FITbutton->setEnabled(false);
    this->clearFITdata();
    return;
  }
  if (this->VOIsegmentIDs.empty()) {
    d->FITbutton->setEnabled(false);
    this->clearFITdata();
    return;
  }
  if (this->modelsID.empty()) {
    d->FITbutton->setEnabled(false);
    this->clearFITdata();
    return;
  }
  d->FITbutton->setEnabled(true);
}

void qSlicerKMAPModuleWidget::enableFITTCMImgbutton() {
  Q_D(qSlicerKMAPModuleWidget);
  if (this->IFID=="") {
    d->FITbuttonTCMImg->setEnabled(false);
    // this->clearFITTCMImgdata();
    return;
  }
  if (this->modelsTCMImgID.empty()) {
    d->FITbuttonTCMImg->setEnabled(false);
    // this->clearFITTCMImgdata();
    return;
  }
  d->FITbuttonTCMImg->setEnabled(true);
}

void qSlicerKMAPModuleWidget::enableFITMTGAbutton() {
  Q_D(qSlicerKMAPModuleWidget);
  if (this->IFID=="") {
    d->FITMTGAbutton->setEnabled(false);
    this->clearFITMTGAdata();
    return;
  }
  if (this->VOIMTGAsegmentIDs.empty()) {
    d->FITMTGAbutton->setEnabled(false);
    this->clearFITMTGAdata();
    return;
  }
  if (this->modelsMTGAID.empty()) {
    d->FITMTGAbutton->setEnabled(false);
    this->clearFITMTGAdata();
    return;
  }
  d->FITMTGAbutton->setEnabled(true);
}

void qSlicerKMAPModuleWidget::enableFITMTGAImgbutton() {
  Q_D(qSlicerKMAPModuleWidget);
  if (this->IFID=="") {
    d->FITbuttonMTGAImg->setEnabled(false);
    // this->clearFITMTGAImgdata();
    return;
  }
  if (this->modelsMTGAImgID.empty()) {
    d->FITbuttonMTGAImg->setEnabled(false);
    // this->clearFITMTGAImgdata();
    return;
  }
  d->FITbuttonMTGAImg->setEnabled(true);
}



void qSlicerKMAPModuleWidget::onVOISegmentsChanged()
{
  Q_D(qSlicerKMAPModuleWidget);
  std::vector<std::string> VOIselectedSegmentIDs;

  for (int i = 0; i < d->VOICheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->VOICheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      std :: string segmentID = checkbox->property("SegmentID").toString().toStdString();
      VOIselectedSegmentIDs.push_back(segmentID);
    }
  }
  this->VOIsegmentIDs = VOIselectedSegmentIDs;

  this->enableFITbutton();
}

void qSlicerKMAPModuleWidget::onVOIMTGASegmentsChanged()
{
  Q_D(qSlicerKMAPModuleWidget);
  std::vector<std::string> VOIselectedSegmentIDs;

  for (int i = 0; i < d->VOIMTGACheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->VOIMTGACheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      std :: string segmentID = checkbox->property("SegmentID").toString().toStdString();
      VOIselectedSegmentIDs.push_back(segmentID);
    }
  }
  this->VOIMTGAsegmentIDs = VOIselectedSegmentIDs;

  this->enableFITMTGAbutton();
}

std::vector<double> qSlicerKMAPModuleWidget::extractColumn(const std::vector<std::vector<double>>& mat, const int index)
{
    std::vector<double> col;
    col.reserve(mat.size());
    for (const auto& row : mat)
    {
        col.push_back(row[index]); // assumes at least one column
    }
    return col;
}

void qSlicerKMAPModuleWidget::onFITbutton()
{
  Q_D(qSlicerKMAPModuleWidget);

  if (segmentTACsnames.empty() || segmentTACs.empty()) {
    std::cerr << "Missing TACs!" << std::endl;
    return;
  }

  if (durations.empty() || timePoints.empty()) {
    std::cerr << "Missing frame time information!" << std::endl;
    return;
  }

  int statIDQString = d->StatSelector->currentIndex();
  if (statIDQString<0) {
    std::cerr << "Missing stat choice!" << std::endl;
    return;
  }
  std::string currentSelectedStatID = d->StatSelector->itemData(statIDQString).toString().toStdString();


  if (this->IFID.empty()) {
    std::cerr << "Missing input function!" << std::endl;
    return;
  }

  if (VOIsegmentIDs.empty()) {
    std::cerr << "Missing VOIs to fit!" << std::endl;
    return;
  }

  if (modelsID.empty()) {
    std::cerr << "Missing Models to fit!" << std::endl;
    return;
  }

  // Run TAC computation
  vtkSlicerKMAPLogic* logic = vtkSlicerKMAPLogic::SafeDownCast(this->logic());
  if (!logic) {
    std::cerr << "Missing Logic!" << std::endl;
    return;
  }

  // Collect parameters using a lambda for brevity
  auto getParamTriplet = [&](QLineEdit* init, QLineEdit* lb, QLineEdit* ub) {
    return std::tuple<double, double, double>{
      init->text().toDouble(), lb->text().toDouble(), ub->text().toDouble()
    };
  };

  auto [k1Init, k1Lower, k1Upper] = getParamTriplet(d->k1Initial, d->k1Lower, d->k1Upper);
  auto [k2Init, k2Lower, k2Upper] = getParamTriplet(d->k2Initial, d->k2Lower, d->k2Upper);
  auto [k3Init, k3Lower, k3Upper] = getParamTriplet(d->k3Initial, d->k3Lower, d->k3Upper);
  auto [k4Init, k4Lower, k4Upper] = getParamTriplet(d->k4Initial, d->k4Lower, d->k4Upper);
  auto [vbInit, vbLower, vbUpper] = getParamTriplet(d->vbInitial, d->vbLower, d->vbUpper);
  auto [tdInit, tdLower, tdUpper] = getParamTriplet(d->tdInitial, d->tdLower, d->tdUpper);

  const long Nframe = timePoints.size();
  const long Nvox = 1;

  std::vector<std::vector<double>> framing;
  framing.reserve(Nframe);
  for (double d : durations)
  {
    framing.emplace_back(1, d);  // Adds a vector with 1 element (column vector)
  }

  const std::vector<double>* wgt = nullptr;
  std::map< std::string, std::vector<double>> wgtVec;

  std::map<std::string, std::vector<std::vector<double>>> tac;
  std::map<std::string, std::vector<bool>> keeptacvec;
  for (const auto& [segmentName, statsVec] : segmentTACs)
  {
    if (statsVec.size() != static_cast<size_t>(Nframe))
    {
      std::cerr << "Mismatch in TAC frame size for segment " << segmentName << std::endl;
      return;
    }

    tac[segmentName].reserve(Nframe);
    keeptacvec[segmentName].reserve(Nframe);
    for (int ivs=0; ivs<statsVec.size(); ++ivs)
    {
      const auto& vs = statsVec[ivs];
      double value;
      if (currentSelectedStatID == "Mean") {
        value = vs.mean;
        if (d->weightFitCheckBox->isChecked()) {
          wgtVec[segmentName].push_back(1./(vs.stddev+1e-16));
        } else {
          wgtVec[segmentName].push_back(1.);
        }
      }
      else if (currentSelectedStatID == "Median") {
        value = vs.median;
        if (d->weightFitCheckBox->isChecked()) {
          wgtVec[segmentName].push_back(1./(vs.iqr+1e-16));
        } else {
          wgtVec[segmentName].push_back(1.);
        }
      }
      else
      {
        vtkGenericWarningMacro("Unknown stat: " << currentSelectedStatID);
        return;
      }
      if (!vs.keep) {
        wgtVec[segmentName][ivs] = 0.;
      }
      tac[segmentName].emplace_back(1, value);  // Adds one-element row (column vector)
      keeptacvec[segmentName].push_back(vs.keep);  // Adds one-element row (column vector)
    }
  }
  this->segmentTAC4TCMfits = tac;
  this->segmentkeep4TCMfits = keeptacvec;

  const double dk = d->decayConstEdit->text().toDouble();
  const double timestep = d->timeStepEdit->text().toDouble();

  auto [pbrp1, pbrp2, pbrp3] = getParamTriplet(d->pbrp1Edit, d->pbrp2Edit, d->pbrp3Edit);
  const double pbrp[] = {pbrp1, pbrp2, pbrp3};
  const int maxiter = d->maxIterTCMEdit->text().toInt();

  std::vector<std::vector<double>> Cp = tac[IFID];

  for (const std::string& segmentID : VOIsegmentIDs)
  {
    const auto& tacVOI = tac[segmentID];
    auto tac_flatten = extractColumn(tacVOI);
    wgt = &wgtVec[segmentID];

    for (const std::string& modelID : modelsID)
    {
      if (modelID == "1TCM") {
        bool sens[] = {true, true, true, false};
        double lb_1tcm[]   = {vbLower, k1Lower, k2Lower, 0.};
        double ub_1tcm[]   = {vbUpper, k1Upper, k2Upper, 0.};
        double init_1tcm[] = {vbInit,  k1Init,  k2Init,  0.};
        logic->callTCM(tacVOI, Cp, framing, Nframe, Nvox,
                       init_1tcm, lb_1tcm, ub_1tcm, sens,
                       dk, timestep, pbrp, maxiter, 1,
                       this->segmentTCM[segmentID]["1TCM"],
                       this->segmentTCMfits[segmentID]["1TCM"],
                       wgt
                   );
        // logic->getFittedTCM(this->segmentTCMfits[segmentID]["1TCM"],
        //                     Cp, framing, Nframe, Nvox, init_1tcm, lb_1tcm,
        //                     ub_1tcm, sens, dk, timestep, pbrp, maxiter,
        //                     1, this->segmentTCM[segmentID]["1TCM"]);
      }
      else if (modelID == "1TdCM") {
        bool sens[] = {true, true, true, true};
        double lb_1tcm[]   = {vbLower, k1Lower, k2Lower, tdLower};
        double ub_1tcm[]   = {vbUpper, k1Upper, k2Upper, tdUpper};
        double init_1tcm[] = {vbInit,  k1Init,  k2Init,  tdInit};
        logic->callTCM(tacVOI, Cp, framing, Nframe, Nvox,
                       init_1tcm, lb_1tcm, ub_1tcm, sens,
                       dk, timestep, pbrp, maxiter, 1,
                       this->segmentTCM[segmentID]["1TdCM"],
                       this->segmentTCMfits[segmentID]["1TdCM"],
                       wgt
                   );
        // logic->getFittedTCM(this->segmentTCMfits[segmentID]["1TdCM"],
        //                     Cp, framing, Nframe, Nvox, init_1tcm, lb_1tcm,
        //                     ub_1tcm, sens, dk, timestep, pbrp, maxiter,
        //                     1, this->segmentTCM[segmentID]["1TdCM"]);
      }
      else if (modelID == "1TiCM") {
        bool sens[] = {true, true, false, false};
        double lb_1tcm[]   = {vbLower, k1Lower, 0., 0.};
        double ub_1tcm[]   = {vbUpper, k1Upper, 0., 0.};
        double init_1tcm[] = {vbInit,  k1Init,  0.,  0.};
        logic->callTCM(tacVOI, Cp, framing, Nframe, Nvox,
                       init_1tcm, lb_1tcm, ub_1tcm, sens,
                       dk, timestep, pbrp, maxiter, 1,
                       this->segmentTCM[segmentID]["1TiCM"],
                       this->segmentTCMfits[segmentID]["1TiCM"],
                       wgt
                   );
        // logic->getFittedTCM(this->segmentTCMfits[segmentID]["1TiCM"],
        //                     Cp, framing, Nframe, Nvox, init_1tcm, lb_1tcm,
        //                     ub_1tcm, sens, dk, timestep, pbrp, maxiter,
        //                     1, this->segmentTCM[segmentID]["1TiCM"]);
      }
      else if (modelID == "1TidCM") {
        bool sens[] = {true, true, false, true};
        double lb_1tcm[]   = {vbLower, k1Lower, 0., tdLower};
        double ub_1tcm[]   = {vbUpper, k1Upper, 0., tdUpper};
        double init_1tcm[] = {vbInit,  k1Init,  0.,  tdInit};
        logic->callTCM(tacVOI, Cp, framing, Nframe, Nvox,
                       init_1tcm, lb_1tcm, ub_1tcm, sens,
                       dk, timestep, pbrp, maxiter, 1,
                       this->segmentTCM[segmentID]["1TidCM"],
                       this->segmentTCMfits[segmentID]["1TidCM"],
                       wgt
                   );
        // logic->getFittedTCM(this->segmentTCMfits[segmentID]["1TidCM"],
        //                     Cp, framing, Nframe, Nvox, init_1tcm, lb_1tcm,
        //                     ub_1tcm, sens, dk, timestep, pbrp, maxiter,
        //                     1, this->segmentTCM[segmentID]["1TidCM"]);
      }
      else if (modelID == "2TCM") {
        bool sens[] = {true, true, true, true, true, false};
        double lb_2tcm[]   = {vbLower, k1Lower, k2Lower, k3Lower, k4Lower, 0.};
        double ub_2tcm[]   = {vbUpper, k1Upper, k2Upper, k3Upper, k4Upper, 0.};
        double init_2tcm[] = {vbInit,  k1Init,  k2Init,  k3Init,  k4Init,  0.};
        logic->callTCM(tacVOI, Cp, framing, Nframe, Nvox,
                       init_2tcm, lb_2tcm, ub_2tcm, sens,
                       dk, timestep, pbrp, maxiter, 2,
                       this->segmentTCM[segmentID]["2TCM"],
                       this->segmentTCMfits[segmentID]["2TCM"],
                       wgt
                   );
        // logic->getFittedTCM(this->segmentTCMfits[segmentID]["2TCM"],
        //                     Cp, framing, Nframe, Nvox, init_2tcm, lb_2tcm,
        //                     ub_2tcm, sens, dk, timestep, pbrp, maxiter,
        //                     1, this->segmentTCM[segmentID]["2TCM"]);
      }
      else if (modelID == "2dTCM") {
        bool sens[] = {true, true, true, true, true, true};
        double lb_2tcm[]   = {vbLower, k1Lower, k2Lower, k3Lower, k4Lower, tdLower};
        double ub_2tcm[]   = {vbUpper, k1Upper, k2Upper, k3Upper, k4Upper, tdUpper};
        double init_2tcm[] = {vbInit,  k1Init,  k2Init,  k3Init,  k4Init,  tdInit};
        logic->callTCM(tacVOI, Cp, framing, Nframe, Nvox,
                       init_2tcm, lb_2tcm, ub_2tcm, sens,
                       dk, timestep, pbrp, maxiter, 2,
                       this->segmentTCM[segmentID]["2dTCM"],
                       this->segmentTCMfits[segmentID]["2dTCM"],
                       wgt
                     );
        // logic->getFittedTCM(this->segmentTCMfits[segmentID]["2dTCM"],
        //                     Cp, framing, Nframe, Nvox, init_2tcm, lb_2tcm,
        //                     ub_2tcm, sens, dk, timestep, pbrp, maxiter,
        //                     1, this->segmentTCM[segmentID]["2dTCM"]);
      }
      else if (modelID == "2TiCM") {
        bool sens[] = {true, true, true, true, false, false};
        double lb_2tcm[]   = {vbLower, k1Lower, k2Lower, k3Lower, 0., 0.};
        double ub_2tcm[]   = {vbUpper, k1Upper, k2Upper, k3Upper, 0., 0.};
        double init_2tcm[] = {vbInit,  k1Init,  k2Init,  k3Init,  0., 0.};
        logic->callTCM(tacVOI, Cp, framing, Nframe, Nvox,
                       init_2tcm, lb_2tcm, ub_2tcm, sens,
                       dk, timestep, pbrp, maxiter, 2,
                       this->segmentTCM[segmentID]["2TiCM"],
                       this->segmentTCMfits[segmentID]["2TiCM"],
                       wgt
                     );
        // logic->getFittedTCM(this->segmentTCMfits[segmentID]["2TiCM"],
        //                     Cp, framing, Nframe, Nvox, init_2tcm, lb_2tcm,
        //                     ub_2tcm, sens, dk, timestep, pbrp, maxiter,
        //                     1, this->segmentTCM[segmentID]["2TiCM"]);
        std::vector<double> fittedTCMvalues(this->segmentTCMfits[segmentID]["2TiCM"],
                                            this->segmentTCMfits[segmentID]["2TiCM"] + Nframe);
      }
      else if (modelID == "2TidCM") {
        bool sens[] = {true, true, true, true, false, true};
        double lb_2tcm[]   = {vbLower, k1Lower, k2Lower, k3Lower, 0., tdLower};
        double ub_2tcm[]   = {vbUpper, k1Upper, k2Upper, k3Upper, 0., tdUpper};
        double init_2tcm[] = {vbInit,  k1Init,  k2Init,  k3Init,  0.,  tdInit};
        logic->callTCM(tacVOI, Cp, framing, Nframe, Nvox,
                       init_2tcm, lb_2tcm, ub_2tcm, sens,
                       dk, timestep, pbrp, maxiter, 2,
                       this->segmentTCM[segmentID]["2TidCM"],
                       this->segmentTCMfits[segmentID]["2TidCM"],
                       wgt
                     );
        // logic->getFittedTCM(this->segmentTCMfits[segmentID]["2TidCM"],
        //                     Cp, framing, Nframe, Nvox, init_2tcm, lb_2tcm,
        //                     ub_2tcm, sens, dk, timestep, pbrp, maxiter,
        //                     1, this->segmentTCM[segmentID]["2TidCM"]);
      } else {
        std::cerr << "Unknown model ID: " << modelID << std::endl;
        return;
      }
    }
  }
  d->populateResultsVOI();
}

void qSlicerKMAPModuleWidget::onFITTCMImgbutton()
{
  Q_D(qSlicerKMAPModuleWidget);

  if ((this->PETdims[0]==0) | (this->PETdims[1]==0) | (this->PETdims[2]==0)) {
    QMessageBox::warning(nullptr,
                         tr("Missing PET dimensions"),
                         tr("Something went wrong to determine PET dimensions."));
  }
  if (this->numberOfTimepoints<1) {
    QMessageBox::warning(nullptr,
                         tr("Missing dynamic PET"),
                         tr("Dynamic PET has non-positive number of timepoints."));
  }

  if (this->PET_flatten_values.empty()) {
    QMessageBox::warning(nullptr,
                         tr("Missing dynamic PET values"),
                         tr("Dynamic PET values are empty."));
  }

  int numVoxels = this->PETdims[0]*this->PETdims[1]*this->PETdims[2];
  if (this->PET_flatten_values.size()!=numVoxels) {
    QMessageBox::warning(nullptr,
                         tr("Dynamic PET: mismatch values"),
                         tr("Number of dynamic PET is not as expected."));
  }

  if (this->PET_flatten_values[0].size()!=this->numberOfTimepoints) {
    QMessageBox::warning(nullptr,
                         tr("Dynamic PET: mismatch values"),
                         tr("Number of timepoints is not as expected."));
  }

  if (segmentTACsnames.empty() || segmentTACs.empty()) {
    std::cerr << "Missing TACs!" << std::endl;
    return;
  }

  if (durations.empty() || timePoints.empty()) {
    std::cerr << "Missing frame time information!" << std::endl;
    return;
  }

  int statIDQString = d->StatSelectorImg->currentIndex();
  if (statIDQString<0) {
    std::cerr << "Missing stat choice!" << std::endl;
    return;
  }
  std::string currentSelectedStatID = d->StatSelectorImg->itemData(statIDQString).toString().toStdString();


  if (this->IFID.empty()) {
    std::cerr << "Missing input function!" << std::endl;
    return;
  }

  if (modelsTCMImgID.empty()) {
    std::cerr << "Missing Models to fit!" << std::endl;
    return;
  }

  // Run TAC computation
  vtkSlicerKMAPLogic* logic = vtkSlicerKMAPLogic::SafeDownCast(this->logic());
  if (!logic) {
    std::cerr << "Missing Logic!" << std::endl;
    return;
  }

  // Collect parameters using a lambda for brevity
  auto getParamTriplet = [&](QLineEdit* init, QLineEdit* lb, QLineEdit* ub) {
    return std::tuple<double, double, double>{
      init->text().toDouble(), lb->text().toDouble(), ub->text().toDouble()
    };
  };

  auto [k1Init, k1Lower, k1Upper] = getParamTriplet(d->k1InitialImg, d->k1LowerImg, d->k1UpperImg);
  auto [k2Init, k2Lower, k2Upper] = getParamTriplet(d->k2InitialImg, d->k2LowerImg, d->k2UpperImg);
  auto [k3Init, k3Lower, k3Upper] = getParamTriplet(d->k3InitialImg, d->k3LowerImg, d->k3UpperImg);
  auto [k4Init, k4Lower, k4Upper] = getParamTriplet(d->k4InitialImg, d->k4LowerImg, d->k4UpperImg);
  auto [vbInit, vbLower, vbUpper] = getParamTriplet(d->vbInitialImg, d->vbLowerImg, d->vbUpperImg);
  auto [tdInit, tdLower, tdUpper] = getParamTriplet(d->tdInitialImg, d->tdLowerImg, d->tdUpperImg);

  const long Nframe = timePoints.size();
  const long Nvox = 1;

  std::vector<std::vector<double>> framing;
  framing.reserve(Nframe);
  for (double d : durations)
  {
    framing.emplace_back(1, d);  // Adds a vector with 1 element (column vector)
  }

  const std::vector<double>* wgt = nullptr;
  std::map< std::string, std::vector<double>> wgtVec;

  std::map<std::string, std::vector<std::vector<double>>> tac;
  std::map<std::string, std::vector<bool>> keeptacvec;

  std :: string& segmentName = this->IFID;
  const auto& statsVec = segmentTACs[segmentName];

  if (statsVec.size() != static_cast<size_t>(Nframe))
  {
    std::cerr << "Mismatch in TAC frame size for segment " << segmentName << std::endl;
    return;
  }

  tac[segmentName].reserve(Nframe);
  keeptacvec[segmentName].reserve(Nframe);
  for (int ivs=0; ivs<statsVec.size(); ++ivs)
  {
    const auto& vs = statsVec[ivs];
    double value;
    if (currentSelectedStatID == "Mean") {
      value = vs.mean;
      if (d->weightFitCheckBox->isChecked()) {
        wgtVec[segmentName].push_back(1./(vs.stddev+1e-16));
      } else {
        wgtVec[segmentName].push_back(1.);
      }
    }
    else if (currentSelectedStatID == "Median") {
      value = vs.median;
      if (d->weightFitCheckBox->isChecked()) {
        wgtVec[segmentName].push_back(1./(vs.iqr+1e-16));
      } else {
        wgtVec[segmentName].push_back(1.);
      }
    }
    else
    {
      vtkGenericWarningMacro("Unknown stat: " << currentSelectedStatID);
      return;
    }
    if (!vs.keep) {
      wgtVec[segmentName][ivs] = 0.;
    }
    tac[segmentName].emplace_back(1, value);  // Adds one-element row (column vector)
    keeptacvec[segmentName].push_back(vs.keep);  // Adds one-element row (column vector)
  }

  const double dk = d->decayConstEditImg->text().toDouble();
  const double timestep = d->timeStepEditImg->text().toDouble();

  auto [pbrp1, pbrp2, pbrp3] = getParamTriplet(d->pbrp1EditImg, d->pbrp2EditImg, d->pbrp3EditImg);
  const double pbrp[] = {pbrp1, pbrp2, pbrp3};
  const int maxiter = d->maxIterTCMEditImg->text().toInt();
  const int numThreads = d->numThreadsTCM->text().toInt();

  std::vector<std::vector<double>> Cp = tac[IFID];
  wgt = &wgtVec[segmentName];
  auto Cp_flatten = extractColumn(Cp);
  std::vector<double> framing_flatten = extractColumn(framing);

  this->stopRequested = false;
  TCMWorker* worker = new TCMWorker(
    logic,
    this->PET_flatten_values,  // voxels
    Cp_flatten,
    framing_flatten,
    modelsTCMImgID,            // vector<string> of model IDs
    vbInit, vbLower, vbUpper,
    k1Init, k1Lower, k1Upper,
    k2Init, k2Lower, k2Upper,
    k3Init, k3Lower, k3Upper,
    k4Init, k4Lower, k4Upper,
    tdInit, tdLower, tdUpper,
    dk,
    timestep,
    pbrp,
    maxiter,
    this->stopRequested,
    wgt,
    numThreads
  );
  this->ProgressBar->setMinimum(0);
  this->ProgressBar->setMaximum(100);
  this->ProgressBar->setValue(0);

  QObject::connect(worker, &TCMWorker::modelStarted, this, [this](const QString& modelID){
      this->ProgressBar->setFormat("Fitting " + modelID + " (%p%)");
      this->ProgressBar->setVisible(true);
      this->ProgressBar->setValue(0);
      this->stopButton->setVisible(true);
      this->stopButton->show();
  });

  QObject::connect(worker, &TCMWorker::progressChanged, this, [this](int value){
      this->ProgressBar->setValue(value);
  });

  QObject::connect(worker, &TCMWorker::canceled, this, [this, worker](const QString& modelID){
      this->ProgressBar->setVisible(false);
      this->stopButton->setVisible(false);
  });

  QObject::connect(worker, &TCMWorker::finishedProcessing,
    this, [this, logic, worker](const QString& modelID){//, const std::vector<TCMParameters>& results) {
      this->TCMImgOutcomes[modelID.toStdString()] = std::move(worker->results);

      auto getModelfields = [](const std::string& modelID) -> std::vector<std::string> {
          if (modelID == "1TCM")
              return {"K1", "k2", "vb", "Ki", "DV", "AIC", "MASE", "BIC", "chi2"};
          else if (modelID == "1TdCM")
              return {"K1", "k2", "vb", "td", "Ki", "DV", "AIC", "MASE", "BIC", "chi2"};
          else if (modelID == "1TiCM")
              return {"K1", "vb", "Ki", "AIC", "MASE", "BIC", "chi2"};
          else if (modelID == "1TidCM")
              return {"K1", "vb", "td", "Ki", "AIC", "MASE", "BIC", "chi2"};
          else if (modelID == "2TCM")
              return {"K1", "k2", "k3", "k4", "vb", "Ki", "DV", "AIC", "MASE", "BIC", "chi2"};
          else if (modelID == "2dTCM")
              return {"K1", "k2", "k3", "k4", "vb", "td", "Ki", "DV", "AIC", "MASE", "BIC", "chi2"};
          else if (modelID == "2TiCM")
              return {"K1", "k2", "k3", "vb", "Ki", "AIC", "MASE", "BIC", "chi2"};
          else if (modelID == "2TidCM")
              return {"K1", "k2", "k3", "vb", "td", "Ki", "AIC", "MASE", "BIC", "chi2"};
          else
              return {}; // unknown model
      };

      auto modelfields = getModelfields(modelID.toStdString());

      logic->CreateTCMParametricImages(
          this->TCMImgOutcomes[modelID.toStdString()],
          this->PETdims,
          modelfields,
          modelID.toStdString(),
          this->petNode,
          this->SubjectHierarchyNode,
          this->petID
      );
  });

  connect(this->stopButton, &QPushButton::clicked, this, [this](){
      this->stopRequested = true;
  });

  QObject::connect(worker, &TCMWorker::finishedAll, this, [this, worker](){
      this->ProgressBar->setVisible(false);
      this->stopButton->setVisible(false);
      worker->deleteLater();
  });

  worker->start();

}

void qSlicerKMAPModuleWidget::onFITMTGAbutton()
{
  Q_D(qSlicerKMAPModuleWidget);

  if (segmentTACsnames.empty() || segmentTACs.empty()) {
    std::cerr << "Missing TACs!" << std::endl;
    return;
  }

  if (durations.empty() || timePoints.empty()) {
    std::cerr << "Missing frame time information!" << std::endl;
    return;
  }

  int statIDQString = d->StatSelectorMTGA->currentIndex();
  if (statIDQString<0) {
    std::cerr << "Missing stat choice!" << std::endl;
    return;
  }
  std::string currentSelectedStatID = d->StatSelectorMTGA->itemData(statIDQString).toString().toStdString();


  if (this->IFID.empty()) {
    std::cerr << "Missing input function!" << std::endl;
    return;
  }

  if (this->VOIMTGAsegmentIDs.empty()) {
    std::cerr << "Missing VOIs to fit!" << std::endl;
    return;
  }

  if (this->modelsMTGAID.empty()) {
    std::cerr << "Missing Models to fit!" << std::endl;
    return;
  }

  // Run TAC computation
  vtkSlicerKMAPLogic* logic = vtkSlicerKMAPLogic::SafeDownCast(this->logic());
  if (!logic) {
    std::cerr << "Missing Logic!" << std::endl;
    return;
  }

  const long Nframe = timePoints.size();

  std::vector<std::vector<double>> framing;
  framing.reserve(Nframe);
  for (double d : durations)
  {
    framing.emplace_back(1, d);  // Adds a vector with 1 element (column vector)
  }
  std::vector<double> framing_flatten = extractColumn(framing);

  std::map<std::string, std::vector<std::vector<double>>> tac;
  const std::vector<double>* wgt = nullptr;
  std::map< std::string, std::vector<double>> wgtVec;

  for (const auto& [segmentName, statsVec] : segmentTACs)
  {
    if (statsVec.size() != static_cast<size_t>(Nframe))
    {
      std::cerr << "Mismatch in TAC frame size for segment " << segmentName << std::endl;
      return;
    }

    tac[segmentName].reserve(Nframe);
    for (int ivs=0; ivs<statsVec.size(); ++ivs)
    {
      const auto& vs = statsVec[ivs];
      double value;
      if (currentSelectedStatID == "Mean") {
        value = vs.mean;
        if (d->weightedFitCheckBox->isChecked()) {
          wgtVec[segmentName].push_back(1./(vs.stddev+1e-16));
        } else {
          wgtVec[segmentName].push_back(1.);
        }
      }
      else if (currentSelectedStatID == "Median") {
        value = vs.median;
        if (d->weightedFitCheckBox->isChecked()) {
          wgtVec[segmentName].push_back(1./(vs.iqr+1e-16));
        } else {
          wgtVec[segmentName].push_back(1.);
        }
      } else {
        std::cerr << "Unknown stat: " << currentSelectedStatID << std::endl;
        return;
      }
      if (!vs.keep) {
        wgtVec[segmentName][ivs] = 0.;
      }
      tac[segmentName].emplace_back(1, value);  // Adds one-element row (column vector)
    }
  }

  // const double timeOffset =  d->timeOffsetEdit->text().toDouble();
  const double framingNorm = d->framingNormEdit->text().toDouble();
  const double timeOffset = this->timePoints[d->timeOffsetSlider->value()-1] / framingNorm;
  const bool robust = d->robustFitCheckBox->isChecked();
  const bool std = d->standardizationCheckBox->isChecked();
  const double huber_tune = d->huberTuneEdit->text().toDouble();
  const double tol = d->tolEdit->text().toDouble();
  const int max_iter = d->maxIterEdit->text().toInt();


  std::vector<std::vector<double>> Cp = tac[IFID];
  auto Cp_flatten = extractColumn(Cp);
  for (const std::string& segmentID : this->VOIMTGAsegmentIDs)
  {
    const auto& tacVOI = tac[segmentID];
    auto tac_flatten = extractColumn(tacVOI);
    wgt = &wgtVec[segmentID];

    for (const std::string& modelID : this->modelsMTGAID)
    {
      if (modelID == "Patlak") {
        logic->Patlak(tac_flatten,
                      Cp_flatten,
                      framing_flatten,
                      this->segmentMTGA[segmentID]["Patlak"],
                      wgt,
                      timeOffset,
                      framingNorm,
                      robust,
                      std,
                      huber_tune,
                      tol,
                      max_iter
                      );
      }
      else if (modelID == "Logan") {
        logic->Logan(tac_flatten,
                     Cp_flatten,
                     framing_flatten,
                     this->segmentMTGA[segmentID]["Logan"],
                     wgt,
                     timeOffset,
                     framingNorm,
                     robust,
                     std,
                     huber_tune,
                     tol,
                     max_iter
                     );
      }
      else if (modelID == "RE") {
        logic->RE(tac_flatten,
                  Cp_flatten,
                  framing_flatten,
                  this->segmentMTGA[segmentID]["RE"],
                  wgt,
                  timeOffset,
                  framingNorm,
                  robust,
                  std,
                  huber_tune,
                  tol,
                  max_iter
                 );
      } else {
        std::cerr << "Unknown model ID: " << modelID << std::endl;
        return;
      }
    }
  }
  d->populateResultsVOIMTGA();
}

void qSlicerKMAPModuleWidget::onFITMTGAImgbutton()
{
  Q_D(qSlicerKMAPModuleWidget);

  if ((this->PETdims[0]==0) | (this->PETdims[1]==0) | (this->PETdims[2]==0)) {
    QMessageBox::warning(nullptr,
                         tr("Missing PET dimensions"),
                         tr("Something went wrong to determine PET dimensions."));
    return;
  }
  if (this->numberOfTimepoints<1) {
    QMessageBox::warning(nullptr,
                         tr("Missing dynamic PET"),
                         tr("Dynamic PET has non-positive number of timepoints."));
    return;
  }

  if (this->PET_flatten_values.empty()) {
    QMessageBox::warning(nullptr,
                         tr("Missing dynamic PET values"),
                         tr("Dynamic PET values are empty."));
    return;
  }

  int numVoxels = this->PETdims[0]*this->PETdims[1]*this->PETdims[2];
  if (this->PET_flatten_values.size()!=numVoxels) {
    QMessageBox::warning(nullptr,
                         tr("Dynamic PET: mismatch values"),
                         tr("Number of voxels is not as expected."));
    return;
  }

  if (this->PET_flatten_values[0].size()!=this->numberOfTimepoints) {
    QMessageBox::warning(nullptr,
                         tr("Dynamic PET: mismatch values"),
                         tr("Number of timepoints is not as expected."));
    return;
  }

  if (segmentTACsnames.empty() || segmentTACs.empty()) {
    std::cerr << "Missing TACs!" << std::endl;
    return;
  }

  if (durations.empty() || timePoints.empty()) {
    std::cerr << "Missing frame time information!" << std::endl;
    return;
  }

  int statIDQString = d->StatSelectorImg->currentIndex();
  if (statIDQString<0) {
    std::cerr << "Missing stat choice!" << std::endl;
    return;
  }
  std::string currentSelectedStatID = d->StatSelectorImg->itemData(statIDQString).toString().toStdString();

  if (this->IFID.empty()) {
    std::cerr << "Missing input function!" << std::endl;
    return;
  }

  if (this->modelsMTGAImgID.empty()) {
    std::cerr << "Missing Models to fit!" << std::endl;
    return;
  }

  // Run TAC computation
  vtkSlicerKMAPLogic* logic = vtkSlicerKMAPLogic::SafeDownCast(this->logic());
  if (!logic) {
    std::cerr << "Missing Logic!" << std::endl;
    return;
  }

  const long Nframe = timePoints.size();

  std::vector<std::vector<double>> framing;
  framing.reserve(Nframe);
  for (double d : durations)
  {
    framing.emplace_back(1, d);  // Adds a vector with 1 element (column vector)
  }
  std::vector<double> framing_flatten = extractColumn(framing);
  std::map<std::string, std::vector<std::vector<double>>> tac;
  const std::vector<double>* wgt = nullptr;
  std::map< std::string, std::vector<double>> wgtVec;

  std :: string& segmentName = this->IFID;
  const auto& statsVec = segmentTACs[segmentName];
  if (statsVec.size() != static_cast<size_t>(Nframe))
  {
    std::cerr << "Mismatch in TAC frame size for segment " << segmentName << std::endl;
    return;
  }

  tac[segmentName].reserve(Nframe);
  for (int ivs=0; ivs<statsVec.size(); ++ivs)
  {
    const auto& vs = statsVec[ivs];
    double value;
    if (currentSelectedStatID == "Mean") {
      value = vs.mean;
      if (d->weightedFitCheckBox->isChecked()) {
        wgtVec[segmentName].push_back(1./(vs.stddev+1e-16));
      } else {
        wgtVec[segmentName].push_back(1.);
      }
    }
    else if (currentSelectedStatID == "Median") {
      value = vs.median;
      if (d->weightedFitCheckBox->isChecked()) {
        wgtVec[segmentName].push_back(1./(vs.iqr+1e-16));
      } else {
        wgtVec[segmentName].push_back(1.);
      }
    } else {
      std::cerr << "Unknown stat: " << currentSelectedStatID << std::endl;
      return;
    }
    if (!vs.keep) {
      wgtVec[segmentName][ivs] = 0.;
    }
    tac[segmentName].emplace_back(1, value);  // Adds one-element row (column vector)
  }

  // const double timeOffset =  d->timeOffsetEdit->text().toDouble();
  const double framingNorm = d->framingNormEditImg->text().toDouble();
  const double timeOffset = this->timePoints[d->timeOffsetSliderImg->value()-1] / framingNorm;
  const bool robust = d->robustFitCheckBoxImg->isChecked();
  const bool std = d->standardizationCheckBoxImg->isChecked();
  const double huber_tune = d->huberTuneEditImg->text().toDouble();
  const double tol = d->tolEditImg->text().toDouble();
  const int max_iter = d->maxIterEditImg->text().toInt();
  const int numThreads = d->numThreadsMTGA->text().toInt();


  std::vector<std::vector<double>> Cp = tac[IFID];
  auto Cp_flatten = extractColumn(Cp);
  wgt = &wgtVec[segmentName];

  this->stopRequested = false;
  for (const std::string& modelID : this->modelsMTGAImgID) {
    if (modelID == "Patlak") {
      logic->Patlak4Img(
        this->PET_flatten_values,
        Cp_flatten,
        framing_flatten,
        wgt,
        timeOffset,
        framingNorm,
        robust,
        std,
        huber_tune,
        tol,
        max_iter,
        this->MTGAImgOutcomes[modelID],
        this->stopRequested,
        this->ProgressBar,
        numThreads,
        this->stopButton
      );
      if (this->stopRequested) {
        break;
      }
      logic->CreateMTGAParametricImages(
          this->MTGAImgOutcomes[modelID],
          this->PETdims,
          {"Ki", "Intercept", "AIC", "MASE", "R2", "chi2"},
          modelID,
          this->petNode,
          this->SubjectHierarchyNode,
          this->petID
        );
    } else if (modelID == "Logan") {
      logic->Logan4Img(
        this->PET_flatten_values,
        Cp_flatten,
        framing_flatten,
        wgt,
        timeOffset,
        framingNorm,
        robust,
        std,
        huber_tune,
        tol,
        max_iter,
        this->MTGAImgOutcomes[modelID],
        this->stopRequested,
        this->ProgressBar,
        numThreads,
        this->stopButton
      );
      if (this->stopRequested) {
        break;
      }
      logic->CreateMTGAParametricImages(
          this->MTGAImgOutcomes[modelID],
          this->PETdims,
          {"DV", "Intercept", "AIC", "MASE", "R2", "chi2"},
          modelID,
          this->petNode,
          this->SubjectHierarchyNode,
          this->petID
        );
    } else if (modelID == "RE") {
      logic->RE4Img(
        this->PET_flatten_values,
        Cp_flatten,
        framing_flatten,
        wgt,
        timeOffset,
        framingNorm,
        robust,
        std,
        huber_tune,
        tol,
        max_iter,
        this->MTGAImgOutcomes[modelID],
        this->stopRequested,
        this->ProgressBar,
        numThreads,
        this->stopButton
      );
      if (this->stopRequested) {
        break;
      }
      logic->CreateMTGAParametricImages(
          this->MTGAImgOutcomes[modelID],
          this->PETdims,
          {"DV", "Intercept", "AIC", "MASE", "R2", "chi2"},
          modelID,
          this->petNode,
          this->SubjectHierarchyNode,
          this->petID
        );
    } else {
      std::cerr << "Unknown model ID: " << modelID << std::endl;
      return;
    }
  }
}


void qSlicerKMAPModuleWidget::onModelsAllbutton()
{
  Q_D(qSlicerKMAPModuleWidget);

  d->ModelsCheckContents->blockSignals(true);
  std :: vector < std :: string > previouslySelectedModels;
  for (int i = 0; i < d->ModelsCheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->ModelsCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslySelectedModels.push_back(checkbox->text().toStdString());
    }
  }

  if (previouslySelectedModels.size()==(d->ModelsCheckLayout->count())) {
    this->modelsID.clear();
    for (int i = 0; i < d->ModelsCheckLayout->count(); ++i)
    {
      QLayoutItem* item = d->ModelsCheckLayout->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (checkbox)
      {
        checkbox->blockSignals(true);
        checkbox->setChecked(false);
        checkbox->blockSignals(false);
      }
    }
  } else {
    this->modelsID.clear();
    for (int i = 0; i < d->ModelsCheckLayout->count(); ++i)
    {
      QLayoutItem* item = d->ModelsCheckLayout->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (checkbox)
      {
        checkbox->blockSignals(true);
        checkbox->setChecked(true);
        checkbox->blockSignals(false);
        this->modelsID.push_back(checkbox->text().toStdString());
      }
    }
  }
  d->ModelsCheckContents->blockSignals(false);
  this->enableFITbutton();

}

void qSlicerKMAPModuleWidget::onModelsMTGAAllbutton()
{
  Q_D(qSlicerKMAPModuleWidget);

  d->ModelsMTGACheckContents->blockSignals(true);
  std :: vector < std :: string > previouslySelectedModels;
  for (int i = 0; i < d->ModelsMTGACheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->ModelsMTGACheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslySelectedModels.push_back(checkbox->text().toStdString());
    }
  }

  if (previouslySelectedModels.size()==(d->ModelsMTGACheckLayout->count())) {
    this->modelsMTGAID.clear();
    for (int i = 0; i < d->ModelsMTGACheckLayout->count(); ++i)
    {
      QLayoutItem* item = d->ModelsMTGACheckLayout->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (checkbox)
      {
        checkbox->blockSignals(true);
        checkbox->setChecked(false);
        checkbox->blockSignals(false);
      }
    }
  } else {
    this->modelsMTGAID.clear();
    for (int i = 0; i < d->ModelsMTGACheckLayout->count(); ++i)
    {
      QLayoutItem* item = d->ModelsMTGACheckLayout->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (checkbox)
      {
        checkbox->blockSignals(true);
        checkbox->setChecked(true);
        checkbox->blockSignals(false);
        this->modelsMTGAID.push_back(checkbox->text().toStdString());
      }
    }
  }
  d->ModelsMTGACheckContents->blockSignals(false);
  this->enableFITMTGAbutton();

}

void qSlicerKMAPModuleWidget::onModelsSelectAllMTGAImgbutton()
{
  Q_D(qSlicerKMAPModuleWidget);

  d->ModelsCheckContentsMTGAImg->blockSignals(true);
  std :: vector < std :: string > previouslySelectedModels;
  for (int i = 0; i < d->ModelsCheckLayoutMTGAImg->count(); ++i)
  {
    QLayoutItem* item = d->ModelsCheckLayoutMTGAImg->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslySelectedModels.push_back(checkbox->text().toStdString());
    }
  }

  if (previouslySelectedModels.size()==(d->ModelsCheckLayoutMTGAImg->count())) {
    this->modelsMTGAImgID.clear();
    for (int i = 0; i < d->ModelsCheckLayoutMTGAImg->count(); ++i)
    {
      QLayoutItem* item = d->ModelsCheckLayoutMTGAImg->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (checkbox)
      {
        checkbox->blockSignals(true);
        checkbox->setChecked(false);
        checkbox->blockSignals(false);
      }
    }
  } else {
    this->modelsMTGAImgID.clear();
    for (int i = 0; i < d->ModelsCheckLayoutMTGAImg->count(); ++i)
    {
      QLayoutItem* item = d->ModelsCheckLayoutMTGAImg->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (checkbox)
      {
        checkbox->blockSignals(true);
        checkbox->setChecked(true);
        checkbox->blockSignals(false);
        this->modelsMTGAImgID.push_back(checkbox->text().toStdString());
      }
    }
  }
  d->ModelsCheckContentsMTGAImg->blockSignals(false);
  this->enableFITMTGAImgbutton();

}

void qSlicerKMAPModuleWidget::onModelsChanged()
{
  Q_D(const qSlicerKMAPModuleWidget);

  this->modelsID.clear();
  for (int i = 0; i < d->ModelsCheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->ModelsCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      this->modelsID.push_back(checkbox->text().toStdString());
    }
  }
  this->enableFITbutton();
  return ;
}

void qSlicerKMAPModuleWidget::onModelsTCMImgChanged()
{
  Q_D(const qSlicerKMAPModuleWidget);

  this->modelsTCMImgID.clear();
  for (int i = 0; i < d->ModelsCheckLayoutTCMImg->count(); ++i)
  {
    QLayoutItem* item = d->ModelsCheckLayoutTCMImg->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      this->modelsTCMImgID.push_back(checkbox->text().toStdString());
    }
  }
  this->enableFITTCMImgbutton();
  return ;
}

void qSlicerKMAPModuleWidget::onModelsMTGAChanged()
{
  Q_D(const qSlicerKMAPModuleWidget);

  this->modelsMTGAID.clear();
  for (int i = 0; i < d->ModelsMTGACheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->ModelsMTGACheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      this->modelsMTGAID.push_back(checkbox->text().toStdString());
    }
  }
  this->enableFITMTGAbutton();
  return ;
}

void qSlicerKMAPModuleWidget::onModelsMTGAImgChanged()
{
  Q_D(const qSlicerKMAPModuleWidget);

  this->modelsMTGAImgID.clear();
  for (int i = 0; i < d->ModelsCheckLayoutMTGAImg->count(); ++i)
  {
    QLayoutItem* item = d->ModelsCheckLayoutMTGAImg->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      this->modelsMTGAImgID.push_back(checkbox->text().toStdString());
    }
  }
  this->enableFITMTGAImgbutton();
  return ;
}

void qSlicerKMAPModuleWidget::onModelsTCMSelectAllbutton()
{
  Q_D(qSlicerKMAPModuleWidget);
  d->ModelsTCMCheckContents->blockSignals(true);
  QSet<QString> previouslySelectedIDs;
  for (int i = 0; i < d->ModelsTCMCheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->ModelsTCMCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslySelectedIDs.insert(checkbox->text());
    }
  }

  if (previouslySelectedIDs.size()==(d->ModelsTCMCheckLayout->count()-1)) {
    for (int i = 0; i < d->ModelsTCMCheckLayout->count(); ++i)
    {
      QWidget* widget = d->ModelsTCMCheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->blockSignals(true);
        cb->setChecked(false);
        cb->blockSignals(false);
      }
    }
  } else {
    for (int i = 0; i < d->ModelsTCMCheckLayout->count(); ++i)
    {
      QWidget* widget = d->ModelsTCMCheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->blockSignals(true);
        cb->setChecked(true);
        cb->blockSignals(false);
      }
    }
  }
  d->ModelsTCMCheckContents->blockSignals(false);
  d->populateModelsTCM(this->plotTCMVOI);
}

void qSlicerKMAPModuleWidget::onModelsSelectAllTCMImgbutton()
{
  Q_D(qSlicerKMAPModuleWidget);

  d->ModelsCheckContentsTCMImg->blockSignals(true);
  std :: vector < std :: string > previouslySelectedModels;
  for (int i = 0; i < d->ModelsCheckLayoutTCMImg->count(); ++i)
  {
    QLayoutItem* item = d->ModelsCheckLayoutTCMImg->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslySelectedModels.push_back(checkbox->text().toStdString());
    }
  }

  if (previouslySelectedModels.size()==(d->ModelsCheckLayoutTCMImg->count())) {
    this->modelsTCMImgID.clear();
    for (int i = 0; i < d->ModelsCheckLayoutTCMImg->count(); ++i)
    {
      QLayoutItem* item = d->ModelsCheckLayoutTCMImg->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (checkbox)
      {
        checkbox->blockSignals(true);
        checkbox->setChecked(false);
        checkbox->blockSignals(false);
      }
    }
  } else {
    this->modelsTCMImgID.clear();
    for (int i = 0; i < d->ModelsCheckLayoutTCMImg->count(); ++i)
    {
      QLayoutItem* item = d->ModelsCheckLayoutTCMImg->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (checkbox)
      {
        checkbox->blockSignals(true);
        checkbox->setChecked(true);
        checkbox->blockSignals(false);
        this->modelsTCMImgID.push_back(checkbox->text().toStdString());
      }
    }
  }
  d->ModelsCheckContentsTCMImg->blockSignals(false);
  this->enableFITTCMImgbutton();

}


void qSlicerKMAPModuleWidget::onPlotTCMbutton()
{
  Q_D(qSlicerKMAPModuleWidget);
  this->ColNameToSegmentID.clear();
  this->MapPlotSeriesNodeIDToPlot.clear();
  if (this->plotTCMVOI.empty()) {
    return;
  }

  vtkMRMLScene* scene = this->mrmlScene();

  // Get selected VOI to plot TCM fit
  std :: string selectedVOI = this->plotTCMVOI;

  // Get TAC to plot
  std::vector<std::vector<double>> tacvoi = this->segmentTAC4TCMfits[selectedVOI];
  std::vector<bool> keepvoi = this->segmentkeep4TCMfits[selectedVOI];

  // Get TCM models to plot
  std::vector<std::string> PlotSelectedTCMs;
  for (int i = 0; i < d->ModelsTCMCheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->ModelsTCMCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      PlotSelectedTCMs.push_back(checkbox->text().toStdString());
    }
  }

  std::map<std::string, double*> TCMfits = this->segmentTCMfits[selectedVOI];

  if (selectedVOI.empty() || segmentTAC4TCMfits.empty() || PlotSelectedTCMs.empty() || TCMfits.empty())
    return;

  // Clear previous plot/chart/table
  this->RemoveExistingPlotChartAndTable();

  // Create or get table
  vtkSmartPointer<vtkMRMLTableNode> tableNode = this->GetOrCreatePlotTable();

  // Add time column (convert to minutes)
  vtkNew<vtkDoubleArray> timeArray;
  timeArray->SetName("Time (min)");
  for (double t : this->timePoints)
    timeArray->InsertNextValue(t / 60.0);
  tableNode->AddColumn(timeArray);

  // Create plot chart
  vtkMRMLPlotChartNode* chartNode = this->GetOrCreatePlotChart();

  // TAC as scatterpoints
  vtkNew<vtkDoubleArray> tacArray;
  tacArray->SetName("TAC");
  vtkNew<vtkStringArray> labelArray;
  labelArray->SetName("ToolTipLabelTAC");
  for (size_t i = 0; i < tacvoi.size(); ++i)
  {
      // Assuming tacvoi[0] holds measured values for VOI (adapt if structured differently)
      tacArray->InsertNextValue(keepvoi[i] ? tacvoi[i][0] : std::numeric_limits<double>::quiet_NaN());
      std::ostringstream oss;
      oss << "Frame: " << i
          << ", Time(s): " << this->timePoints[i]
          << ", Time(min): " << this->timePoints[i]/60.0;
      labelArray->InsertNextValue(oss.str());
  }
  tableNode->AddColumn(tacArray);
  tableNode->AddColumn(labelArray);

  // TCM fits as line plots
  for (const std::string& modelName : PlotSelectedTCMs)
  {
      auto it = TCMfits.find(modelName);
      if (it == TCMfits.end() || !it->second)
        continue;
      TCMParameters params = this->segmentTCM[selectedVOI][modelName];

      double* fitArrayPtr = it->second;

      vtkNew<vtkDoubleArray> fitArray;
      fitArray->SetName(modelName.c_str());
      for (size_t ivs = 0; ivs < this->timePoints.size(); ++ivs)
      {
        double fv = fitArrayPtr[ivs];
        fitArray->InsertNextValue(fv);
        // if (params.keep[ivs]) {
        //   double fv = fitArrayPtr[ivs];
        //   fitArray->InsertNextValue(fv);
        // } else {
        //   double nextValue = std::numeric_limits<double>::quiet_NaN();
        //   double x1 = std::numeric_limits<double>::quiet_NaN();
        //   for (int next_ivs = ivs+1; next_ivs<this->timePoints.size(); ++next_ivs) {
        //     if (params.keep[next_ivs])
        //     {
        //         x1 = this->timePoints[next_ivs];
        //         nextValue = fitArrayPtr[next_ivs];
        //     }
        //   }
        //   if (std::isnan(nextValue)) {
        //     fitArray->InsertNextValue(std::numeric_limits<double>::quiet_NaN());
        //     continue;
        //   }
        //   double prevValue = std::numeric_limits<double>::quiet_NaN();
        //   double x0 = std::numeric_limits<double>::quiet_NaN();
        //   for (int prev_ivs = ivs-1; prev_ivs>=0; --prev_ivs) {
        //     if (params.keep[prev_ivs])
        //     {
        //         x0 = this->timePoints[prev_ivs];
        //         prevValue = fitArrayPtr[prev_ivs];
        //     }
        //   }
        //   if (std::isnan(prevValue)) {
        //     fitArray->InsertNextValue(std::numeric_limits<double>::quiet_NaN());
        //     continue;
        //   }
        //
        //   double x  = this->timePoints[ivs];
        //   // Proper linear interpolation
        //   double value = prevValue + ((x - x0) / (x1 - x0)) * (nextValue - prevValue);
        //   fitArray->InsertNextValue(value);
        // }
      }

      tableNode->AddColumn(fitArray);

      vtkSmartPointer<vtkMRMLPlotSeriesNode> lineSeries = vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
      scene->AddNode(lineSeries);
      lineSeries->SetName(modelName.c_str());
      lineSeries->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);  // line plot
      lineSeries->SetAndObserveTableNodeID(tableNode->GetID());
      lineSeries->SetXColumnName("Time (min)");
      lineSeries->SetYColumnName(modelName.c_str());
      lineSeries->SetLabelColumnName("ToolTipLabelTAC");
      lineSeries->SetUniqueColor();
      lineSeries->SetMarkerStyle(vtkMRMLPlotSeriesNode::MarkerStyleNone);
      chartNode->AddAndObservePlotSeriesNodeID(lineSeries->GetID());
  }

  this->ColNameToSegmentID["TAC"]=selectedVOI;
  vtkSmartPointer<vtkMRMLPlotSeriesNode> scatterSeries = vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
  scene->AddNode(scatterSeries);
  scatterSeries->SetName("TAC");
  scatterSeries->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);  // scatter points
  scatterSeries->SetAndObserveTableNodeID(tableNode->GetID());
  scatterSeries->SetXColumnName("Time (min)");
  scatterSeries->SetYColumnName("TAC");
  scatterSeries->SetLabelColumnName("ToolTipLabelTAC");
  scatterSeries->SetUniqueColor();
  scatterSeries->SetLineStyle(vtkMRMLPlotSeriesNode::LineStyleNone);
  chartNode->AddAndObservePlotSeriesNodeID(scatterSeries->GetID());
  chartNode->SetTitle(this->segmentTACsnames[selectedVOI].c_str());
  chartNode->SetXAxisTitle("Time (min)");
  chartNode->SetYAxisTitle("SUVbw (g/mL)");

  // Show plot view
  auto* layoutNode = vtkMRMLLayoutNode::SafeDownCast(scene->GetFirstNodeByClass("vtkMRMLLayoutNode"));
  if (layoutNode)
    layoutNode->SetViewArrangement(vtkMRMLLayoutNode::SlicerLayoutConventionalPlotView);

  vtkMRMLPlotViewNode* plotViewNode = vtkMRMLPlotViewNode::SafeDownCast(
      scene->GetFirstNodeByClass("vtkMRMLPlotViewNode"));
  if (plotViewNode)
  {
      plotViewNode->SetPlotChartNodeID(chartNode->GetID());
      // qMRMLPlotWidget* plotWidget = nullptr;
      // if (qSlicerApplication::application())
      // {
      //     qSlicerLayoutManager* layoutManager =
      //         qSlicerApplication::application()->layoutManager();
      //     qMRMLPlotWidget* plotWidget = nullptr;
      //     plotWidget = layoutManager->plotWidget(0);
      //     qMRMLPlotView* plotView = plotWidget->plotView();
      //     if (plotView)
      //     {
      //       QObject::connect(plotView, SIGNAL(dataSelected(vtkStringArray*, vtkCollection*)),
      //                        this, SLOT(onSelectedPoint(vtkStringArray*, vtkCollection*)));
      //     }
      // }
  }

}

void qSlicerKMAPModuleWidget::onPlotMTGAbutton() {
  Q_D(qSlicerKMAPModuleWidget);
  this->ColNameToSegmentID.clear();
  this->MapPlotSeriesNodeIDToPlot.clear();

  if (this->plotMTGAVOI.empty()) {
    return;
  }

  vtkMRMLScene* scene = this->mrmlScene();

  // Get selected VOI for MTGA plot
  std::string selectedVOI = this->plotMTGAVOI;

  // Get selected MTGA model from combo box
  QString modelName_qstr = d->MTGASelector->currentText();
  if (modelName_qstr.toStdString().empty())
    return;
  std::string modelName = modelName_qstr.toStdString();

  // Check data availability
  auto voiIt = this->segmentMTGA.find(selectedVOI);
  if (voiIt == this->segmentMTGA.end())
    return;

  auto modelIt = voiIt->second.find(modelName);
  if (modelIt == voiIt->second.end())
    return;

  const MTGAParameters& params = modelIt->second;
  if (params.x.empty() || params.y.empty() || params.fitted.empty())
    return;

  // Clear previous plot/chart/table
  this->RemoveExistingPlotChartAndTable();

  // Create or get table
  vtkSmartPointer<vtkMRMLTableNode> tableNode = this->GetOrCreatePlotTable();

  vtkNew<vtkDoubleArray> xArray;
  vtkNew<vtkDoubleArray> yArray;
  vtkNew<vtkStringArray> labelArray;
  xArray->SetName("X");
  yArray->SetName("Data");
  labelArray->SetName("ToolTipData");
  for (int i=0;  i < params.x.size(); ++i) {
    double xv = params.x[i];
    double yv = params.y[i];
    xArray->InsertNextValue(xv);
    yArray->InsertNextValue(params.keep[i] ? yv : std::numeric_limits<double>::quiet_NaN());
    std::ostringstream oss;
    oss << "Frame: " << params.frame[i]-1
        << ", Time(s): " << this->timePoints[i]
        << ", Time(min): " << this->timePoints[i]/60.0;
    labelArray->InsertNextValue(oss.str());
  }
  tableNode->AddColumn(xArray);
  tableNode->AddColumn(yArray);
  tableNode->AddColumn(labelArray);


  // Create plot chart
  vtkMRMLPlotChartNode* chartNode = this->GetOrCreatePlotChart();

  // Fitted values as line plot
  vtkNew<vtkDoubleArray> fitArray;
  fitArray->SetName(modelName.c_str());
  for (int ivs=0; ivs<params.fitted.size(); ++ivs) {
    if (params.keep[ivs]) {
      double fv = params.fitted[ivs];
      fitArray->InsertNextValue(fv);
    } else {
      double nextValue = std::numeric_limits<double>::quiet_NaN();
      double x1 = std::numeric_limits<double>::quiet_NaN();
      for (int next_ivs = ivs+1; next_ivs<params.fitted.size(); ++next_ivs) {
        if (params.keep[next_ivs])
        {
            x1 = params.x[next_ivs];
            nextValue = params.fitted[next_ivs];
        }
      }
      if (std::isnan(nextValue)) {
        fitArray->InsertNextValue(std::numeric_limits<double>::quiet_NaN());
        continue;
      }
      double prevValue = std::numeric_limits<double>::quiet_NaN();
      double x0 = std::numeric_limits<double>::quiet_NaN();
      for (int prev_ivs = ivs-1; prev_ivs>=0; --prev_ivs) {
        if (params.keep[prev_ivs])
        {
            x0 = params.x[prev_ivs];
            prevValue = params.fitted[prev_ivs];
        }
      }
      if (std::isnan(prevValue)) {
        fitArray->InsertNextValue(std::numeric_limits<double>::quiet_NaN());
        continue;
      }

      double x  = params.x[ivs];
      // Proper linear interpolation
      double value = prevValue + ((x - x0) / (x1 - x0)) * (nextValue - prevValue);
      fitArray->InsertNextValue(value);
    }
  }
  tableNode->AddColumn(fitArray);

  vtkSmartPointer<vtkMRMLPlotSeriesNode> lineSeries = vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
  scene->AddNode(lineSeries);
  lineSeries->SetName(modelName.c_str());
  lineSeries->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);
  lineSeries->SetAndObserveTableNodeID(tableNode->GetID());
  lineSeries->SetXColumnName("X");
  lineSeries->SetYColumnName(modelName.c_str());
  lineSeries->SetLabelColumnName("ToolTipData");
  lineSeries->SetUniqueColor();
  lineSeries->SetMarkerStyle(vtkMRMLPlotSeriesNode::MarkerStyleNone);
  chartNode->AddAndObservePlotSeriesNodeID(lineSeries->GetID());

  // Scatter series for measured values
  this->ColNameToSegmentID["Data"] = selectedVOI;
  vtkSmartPointer<vtkMRMLPlotSeriesNode> scatterSeries = vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
  scene->AddNode(scatterSeries);
  scatterSeries->SetName("Data");
  scatterSeries->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);
  scatterSeries->SetAndObserveTableNodeID(tableNode->GetID());
  scatterSeries->SetXColumnName("X");
  scatterSeries->SetYColumnName("Data");
  scatterSeries->SetLabelColumnName("ToolTipData");
  scatterSeries->SetUniqueColor();
  scatterSeries->SetLineStyle(vtkMRMLPlotSeriesNode::LineStyleNone);
  chartNode->AddAndObservePlotSeriesNodeID(scatterSeries->GetID());
  chartNode->SetTitle((modelName + " - " + this->segmentTACsnames[selectedVOI]).c_str());
  if (modelName=="Patlak") {
    chartNode->SetXAxisTitle("intCp/Cp");
    chartNode->SetYAxisTitle("Ct/Cp");
  } else if (modelName=="Logan") {
    chartNode->SetXAxisTitle("intCp/Ct");
    chartNode->SetYAxisTitle("intCt/Ct");
  } else if (modelName=="RE") {
    chartNode->SetXAxisTitle("intCp/Cp");
    chartNode->SetYAxisTitle("intCt/Cp");
  } else {
    std::cerr << "Unknown model: " << modelName << std::endl;
  }

  // Show plot view
  auto* layoutNode = vtkMRMLLayoutNode::SafeDownCast(scene->GetFirstNodeByClass("vtkMRMLLayoutNode"));
  if (layoutNode)
    layoutNode->SetViewArrangement(vtkMRMLLayoutNode::SlicerLayoutConventionalPlotView);

  vtkMRMLPlotViewNode* plotViewNode = vtkMRMLPlotViewNode::SafeDownCast(
      scene->GetFirstNodeByClass("vtkMRMLPlotViewNode"));
  if (plotViewNode) {
    plotViewNode->SetPlotChartNodeID(chartNode->GetID());
    // qMRMLPlotWidget* plotWidget = nullptr;
    // if (qSlicerApplication::application())
    // {
    //     qSlicerLayoutManager* layoutManager =
    //         qSlicerApplication::application()->layoutManager();
    //     qMRMLPlotWidget* plotWidget = nullptr;
    //     plotWidget = layoutManager->plotWidget(0);
    //     qMRMLPlotView* plotView = plotWidget->plotView();
    //     if (plotView)
    //     {
    //       QObject::connect(plotView, SIGNAL(dataSelected(vtkStringArray*, vtkCollection*)),
    //                        this, SLOT(onSelectedPoint(vtkStringArray*, vtkCollection*)));
    //     }
    // }
  }
}

bool qSlicerKMAPModuleWidget::checkdisplayedKMAP() {
  vtkMRMLPlotViewNode* plotViewNode = nullptr;
  vtkCollection* viewNodes = this->mrmlScene()->GetNodesByClass("vtkMRMLPlotViewNode");
  if (!viewNodes || viewNodes->GetNumberOfItems() == 0)
  {
      return false; // no plot view node
  }

  plotViewNode = vtkMRMLPlotViewNode::SafeDownCast(viewNodes->GetItemAsObject(0));
  if (!plotViewNode)
  {
      return false;
  }

  // Get the currently displayed plot chart
  vtkMRMLPlotChartNode* currentPlot = vtkMRMLPlotChartNode::SafeDownCast(
      this->mrmlScene()->GetNodeByID(plotViewNode->GetPlotChartNodeID())
  );
  if (!currentPlot)
  {
      return false;
  }

  // Check if the currently displayed plot is your "KMAP.PlotChart"
  if (currentPlot->GetName() == nullptr || std::string(currentPlot->GetName()) != "KMAP.PlotChart")
  {
      return false; // not the KMAP plot, do nothing
  }
  return true;
}

void qSlicerKMAPModuleWidget::onDeleteKeyPressed() {
  Q_D(qSlicerKMAPModuleWidget);
  if (!this->checkdisplayedKMAP())
    return;

  if (this->PlotSelectedFrame == -1)
    return;

  if (this->PlotSelectedVOI == "")
    return;

  if (this->segmentTACs.empty())
    return;

  vtkGenericWarningMacro("Segment " << this->segmentTACsnames[this->PlotSelectedVOI] << " removed at frame " << this->PlotSelectedFrame);
  this->segmentTACs[this->PlotSelectedVOI][this->PlotSelectedFrame].keep = false;
  this->segmentTACs[this->PlotSelectedVOI][this->PlotSelectedFrame].empty = false;
  this->onPlotbutton();
  this->PlotSelectedFrame = -1;
  this->PlotSelectedVOI = "";
  return;
}

void qSlicerKMAPModuleWidget::onResetbutton()
{
  Q_D(qSlicerKMAPModuleWidget);
  vtkGenericWarningMacro("Restoring removed points ");

  if (this->segmentTACs.empty())
    return;

  for (auto& segmentPair : this->segmentTACs)
  {
    const std::string& segmentID = segmentPair.first;
    std::vector<VoxelStatistics>& tacVector = segmentPair.second;

    for (size_t frameID = 0; frameID < tacVector.size(); ++frameID)
    {
      VoxelStatistics& vs = tacVector[frameID];
      if (!vs.empty)
      {
          vs.keep = true;
      }
    }
  }
  if (this->checkdisplayedKMAP()) {
    this->onPlotbutton();
  }
  return;
}
