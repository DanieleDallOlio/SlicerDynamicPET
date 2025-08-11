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

// Qt includes
#include <QDebug>

// Slicer includes
#include "qSlicerKMAPModuleWidget.h"
#include "ui_qSlicerKMAPModuleWidget.h"
#include <QApplication>



//-----------------------------------------------------------------------------
class qSlicerKMAPModuleWidgetPrivate: public Ui_qSlicerKMAPModuleWidget
{
  Q_DECLARE_PUBLIC(qSlicerKMAPModuleWidget);

protected:
  qSlicerKMAPModuleWidget* const q_ptr;
public:
  qSlicerKMAPModuleWidgetPrivate(qSlicerKMAPModuleWidget& object);
  ~qSlicerKMAPModuleWidgetPrivate();
  void init();
  void populatePatientComboBox();
  void populateStudyComboBox(vtkIdType patientID);
  void populateNodeComboBox(QComboBox* comboBox, vtkIdType parentItemID, const char * requiredNodeType, const std :: string requiredModality);
  void populateSegmentCheckboxes(vtkIdType SegItemID);
  void populatePlotSegmentCheckboxes();
  void populateIF();
  void populateVOI(std :: string ifID);
  void populateResultsVOI();
  void populateResultsTable(std :: string segmentID);
  void populateModelsTCM(std :: string segmentID);
};

//-----------------------------------------------------------------------------
// qSlicerKMAPModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerKMAPModuleWidgetPrivate::qSlicerKMAPModuleWidgetPrivate(qSlicerKMAPModuleWidget& object): q_ptr(&object)
{
  Q_Q(qSlicerKMAPModuleWidget);
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
  this->SegSelector->setEnabled(false);
  this->segmentSelectAll->setEnabled(false);
  this->saveExcelButton->setEnabled(false);
  this->KineticsWidget->setEnabled(false);
  this->VOIsegmentSelectAll->setEnabled(false);

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
  // QObject::connect( this->PlotErrorCheckbox, SIGNAL(toggled(bool)),
  //   q, SLOT(onPlotbutton()));
  QObject::connect( this->IFSelector, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onIFSelectionChanged(int)));
  QObject::connect( this->VOIsegmentSelectAll, SIGNAL(clicked(bool)),
    q, SLOT(onVOISelectAllbutton()));
  QObject::connect( this->ModelsSelectAll, SIGNAL(clicked(bool)),
      q, SLOT(onModelsAllbutton()));
  QObject::connect( this->ModelsTCMSelectAll, SIGNAL(clicked(bool)),
      q, SLOT(onModelsTCMSelectAllbutton()));
  QObject::connect( this->FITbutton, SIGNAL(clicked(bool)),
    q, SLOT(onFITbutton()));
  QObject::connect( this->VOISelector, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onVOISelectionChanged(int)));

  this->TACCollapsibleButton->setCollapsed(true);
  this->TCMCollapsibleButton->setCollapsed(true);
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

def DPE_save_multisheet_excel(filepath, sheet_data_dict):
    """
    filepath: str - full path to xlsx
    sheet_data_dict: dict[str, list[list[str]]] - sheet name to 2D table
    """
    with pd.ExcelWriter(filepath, engine="xlsxwriter") as writer:
        for sheet, data in sheet_data_dict.items():
            df = pd.DataFrame(data)[["Time(s)", "Duration", "Mean", "Median", "StDev","IQR","Min", "Max", "Q1", "Q3", "VoxelCount","Volume(mm3)","Volume(cm3)"]]
            df.to_excel(writer, sheet_name=sheet, index=False)
)PYTHON");

  for (const QString& name : q->ModelsNames)
  {
    QCheckBox* cb = new QCheckBox(name, this->ModelsCheckContents);
    this->ModelsCheckLayout->addWidget(cb);
    QObject::connect(cb, SIGNAL(stateChanged(int)),
                q, SLOT(onModelsChanged()));
  }

  for (int i = 0; i < q->StatsNames.size(); ++i)
  {
    this->StatSelector->addItem(q->StatsNames[i], q->StatsNames[i]);
  }

  this->StatSelector->setCurrentIndex(0);

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
    // std::cout << "Patient: " << shNode->GetItemName(id)
    //           << " (ID: " << id << ")" << std::endl;
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
    std::string segmentName = segmentation->GetSegment(segmentID)->GetName();

    QCheckBox* checkbox = new QCheckBox(QString::fromStdString(segmentName));
    checkbox->setProperty("SegmentID", QString::fromStdString(segmentID));

    bool wasSelected = previouslySelectedIDs.contains(QString::fromStdString(segmentID));
    checkbox->setChecked(wasSelected);
    this->segmentCheckLayout->addWidget(checkbox);
    QObject::connect(checkbox, SIGNAL(stateChanged(int)),
                 q, SLOT(onSegmentsChanged()));
    if (wasSelected)
      q->segmentIDs.push_back(QString::fromStdString(segmentID));
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
    this->KineticsWidget->setEnabled(false);
    this->IFSelector->blockSignals(false);
    return;
  }
  this->KineticsWidget->setEnabled(true);

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

void qSlicerKMAPModuleWidgetPrivate::populateResultsTable(std :: string segmentID)
{
  Q_Q(qSlicerKMAPModuleWidget);

  this->TCMResultsTable->clear();
  this->TCMResultsTable->setRowCount(0);
  this->TCMResultsTable->setColumnCount(0);

  if (segmentID.empty())
  {
    this->populateModelsTCM(segmentID);
    return;
  }

  // Define column headers
  QStringList headers = { "", "K1", "k2", "k3", "k4", "vb", "td", "Ki", "DV", "AIC", "MASE" };
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
    this->TCMResultsTable->setItem(row, 1, new QTableWidgetItem(QString::number(params.K1)));
    this->TCMResultsTable->setItem(row, 2, new QTableWidgetItem(QString::number(params.k2)));
    this->TCMResultsTable->setItem(row, 3, new QTableWidgetItem(QString::number(params.k3)));
    this->TCMResultsTable->setItem(row, 4, new QTableWidgetItem(QString::number(params.k4)));
    this->TCMResultsTable->setItem(row, 5, new QTableWidgetItem(QString::number(params.vb)));
    this->TCMResultsTable->setItem(row, 6, new QTableWidgetItem(QString::number(params.td)));
    this->TCMResultsTable->setItem(row, 7, new QTableWidgetItem(QString::number(params.Ki)));
    this->TCMResultsTable->setItem(row, 8, new QTableWidgetItem(QString::number(params.DV)));
    this->TCMResultsTable->setItem(row, 9, new QTableWidgetItem(QString::number(params.AIC)));
    this->TCMResultsTable->setItem(row, 10, new QTableWidgetItem(QString::number(params.MASE)));

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
  this->ProgressBar = new QProgressBar();
  this->ProgressBar->setObjectName(QString::fromUtf8("ProgressBar"));
  this->ProgressBar->setMaximum(100);
  this->ProgressBar->setValue(0);
  this->checkboxNames = QStringList{
    "Mean", "Median", "Min", "Max", "VoxelCount", "Volume(cc)"
  };
  this->ModelsNames = QStringList{
    "Patlak", "Logan", "RE", "1TCM", "1TdCM", "1TiCM", "1TidCM", "2TCM", "2dTCM", "2TiCM", "2TidCM"
  };
  this->StatsNames = QStringList{
    "Mean", "Median"
  };
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

  ctkDICOMDatabase* db = qSlicerApplication::application()->dicomDatabase();
  QString instanceUID = petNode->GetAttribute("DICOM.instanceUIDs");
  QString seriesUID = db->seriesForFile(db->fileForInstance(instanceUID));
  QStringList fileList = db->filesForSeries(seriesUID);

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
  // std::cout << "Selected Patient ID: " << id << std::endl;
  std :: string name = this->SubjectHierarchyNode->GetItemName(this->patID);
  std :: string excelfile = name + "_TAC.xlsx";
  d->fileexcel->setText(QString::fromStdString(excelfile));
  //
  // std::cout << "Name: " << name
  //           << ", ID: " << id << std::endl;
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
  // if (this->stuID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  // {
  //   // "None" selected — ignore or reset state
  //   return;
  // }
  // std::cout << "Selected Study ID: " << id << std::endl;
  // std :: string name = this->SubjectHierarchyNode->GetItemName(id);
  //
  // std::cout << "Study: " << name
  //           << ", ID: " << id << std::endl;

}


void qSlicerKMAPModuleWidget::onCTChanged (int index) {
  Q_D(qSlicerKMAPModuleWidget);
  this->ctID = d->CTSelector->itemData(index).value<vtkIdType>();
  this->enableTACbutton();
  // if (this->ctID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  // {
  //   // "None" selected — ignore or reset state
  //   return;
  // }
  //
  // std :: string name = this->SubjectHierarchyNode->GetItemName(this->ctID);
  // std::cout << "CT: " << name
  //           << ", ID: " << this->ctID << std::endl;
}


void qSlicerKMAPModuleWidget::onPETChanged (int index) {
  Q_D(qSlicerKMAPModuleWidget);
  this->petID = d->PETSelector->itemData(index).value<vtkIdType>();
  this->getDurations();
  this->enableTACbutton();
  // if (this->petID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  // {
  //   // "None" selected — ignore or reset state
  //   return;
  // }
  // std :: string name = this->SubjectHierarchyNode->GetItemName(this->petID);
  // std::cout << "PET: " << name
  //           << ", ID: " << this->petID << std::endl;
}

void qSlicerKMAPModuleWidget::onSegChanged (int index) {
  Q_D(qSlicerKMAPModuleWidget);
  this->segID = d->SegSelector->itemData(index).value<vtkIdType>();
  d->populateSegmentCheckboxes(this->segID);
  this->enableTACbutton();
  // if (this->segID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  // {
  //   // "None" selected — ignore or reset state
  //   return;
  // }
  // std :: string name = this->SubjectHierarchyNode->GetItemName(this->segID);
  // std::cout << "Seg: " << name
  //           << ", ID: " << this->segID << std::endl;
}

void qSlicerKMAPModuleWidget::onSegmentsChanged()
{
  Q_D(qSlicerKMAPModuleWidget);
  std::vector<QString> selectedSegmentIDs;

  for (int i = 0; i < d->segmentCheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->segmentCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      QString segmentID = checkbox->property("SegmentID").toString();
      selectedSegmentIDs.push_back(segmentID);
    }
  }
  this->segmentIDs = selectedSegmentIDs;

  this->enableTACbutton();
}

void qSlicerKMAPModuleWidget::clearTACdata() {
  Q_D(qSlicerKMAPModuleWidget);
  this->segmentTACs.clear();
  this->segmentTACsnames.clear();
  d->populatePlotSegmentCheckboxes();
  d->populateIF();
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
  this->RemoveExistingPlotChartAndTable();
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
  vtkMRMLTableNode* tableNode = vtkMRMLTableNode::SafeDownCast(
    this->mrmlScene()->GetFirstNodeByName("KMAP.PlotTable"));
  if (tableNode)
    this->mrmlScene()->RemoveNode(tableNode);

  vtkMRMLPlotChartNode* chartNode = vtkMRMLPlotChartNode::SafeDownCast(
    this->mrmlScene()->GetFirstNodeByName("KMAP.PlotChart"));
  if (chartNode)
  {
    std::vector<std::string> seriesIDs;
    chartNode->GetPlotSeriesNodeIDs(seriesIDs);
    for (const std::string& id : seriesIDs)
    {
      vtkMRMLNode* node = this->mrmlScene()->GetNodeByID(id);
      if (node)
        this->mrmlScene()->RemoveNode(node);
    }
    this->mrmlScene()->RemoveNode(chartNode);
  }
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
  qApp->processEvents();
  logic->computeTAC(this->ctID, this->petID, this->segID, segmentsToCompute, this->segmentTACs, this->segmentTACsnames, this->ProgressBar);
  this->ProgressBar->setValue(0);
  this->ProgressBar->setVisible(false);
  qApp->processEvents();
  d->populatePlotSegmentCheckboxes();
  d->populateIF();
  return;
}

void qSlicerKMAPModuleWidget::onSelectAllbutton()
{
  Q_D(qSlicerKMAPModuleWidget);
  if (this->segmentIDs.size()==(d->segmentCheckLayout->count()-1)) {
    for (int i = 0; i < d->segmentCheckLayout->count(); ++i)
    {
      QWidget* widget = d->segmentCheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->setChecked(false);
      }
    }
  } else {
    for (int i = 0; i < d->segmentCheckLayout->count(); ++i)
    {
      QWidget* widget = d->segmentCheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->setChecked(true);
      }
    }
  }
  d->populateSegmentCheckboxes(this->segID);
}

void qSlicerKMAPModuleWidget::onExcelPathChanged(const QString& path)
{
  Q_D(qSlicerKMAPModuleWidget);
  d->saveExcelButton->setEnabled(!path.trimmed().isEmpty());
}

QVariantMap qSlicerKMAPModuleWidget::TACtoPythonDict()
{
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

void qSlicerKMAPModuleWidget::onPlotbutton()
{
  Q_D(qSlicerKMAPModuleWidget);

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
  for (double t : this->timePoints)
    timeArray->InsertNextValue(t/60.);
  tableNode->AddColumn(timeArray);

  // Create plot chart
  vtkMRMLPlotChartNode* chartNode = this->GetOrCreatePlotChart();

  for (const std::string& segmentID : PlotSelectedIDs)
  {
    std :: string segmentName = this->segmentTACsnames[segmentID];
    for (const std::string& statName : PlotSelectedStats)
    {
      std::string colName = segmentName + " - " + statName;
      vtkNew<vtkDoubleArray> statArray;
      statArray->SetName(colName.c_str());

      vtkNew<vtkDoubleArray> statErrArray;
      std::string colErrName = colName + " Error";
      statErrArray->SetName(colErrName.c_str());

      for (const VoxelStatistics& vs : this->segmentTACs[segmentID])
      {
        if (statName == "Mean")
        {
          statArray->InsertNextValue(vs.mean);
          if (d->PlotErrorCheckbox && d->PlotErrorCheckbox->isChecked())
            statErrArray->InsertNextValue(vs.stddev);
        }
        else if (statName == "Median")
        {
          statArray->InsertNextValue(vs.median);
          if (d->PlotErrorCheckbox && d->PlotErrorCheckbox->isChecked())
            statErrArray->InsertNextValue(vs.iqr);
        }
        else if (statName == "VoxelCount") statArray->InsertNextValue(vs.count);
        else if (statName == "Min")        statArray->InsertNextValue(vs.min);
        else if (statName == "Max")        statArray->InsertNextValue(vs.max);
        else if (statName == "Volume(cc)") statArray->InsertNextValue(vs.volume_cm3);
        else std::cerr << "Unknown stat name: " << statName << std::endl;
      }

      tableNode->AddColumn(statArray);
      // if (statErrArray->GetNumberOfTuples() > 0)
      //   tableNode->AddColumn(statErrArray);

      vtkSmartPointer<vtkMRMLPlotSeriesNode> series = vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
      scene->AddNode(series);
      series->SetName(colName.c_str());
      series->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);
      series->SetAndObserveTableNodeID(tableNode->GetID());
      series->SetXColumnName("Time (min)");
      series->SetYColumnName(colName.c_str());
      series->SetUniqueColor();


      // series->SetName(colName.c_str());
      // series->SetMarkerStyle(vtkMRMLPlotSeriesNode::MarkerStyleNone);
      // if (statErrArray->GetNumberOfTuples() > 0)
      //   series->SetAttribute("yErrorColumnName", colErrName.c_str());

      // chartNode->SetXAxisLabel("Time (s)");
      // chartNode->SetYAxisLabel("SUVbw (g/mL)");
      chartNode->AddAndObservePlotSeriesNodeID(series->GetID());
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
  }

}

void qSlicerKMAPModuleWidget::onIFSelectionChanged(int index)
{
  Q_D(qSlicerKMAPModuleWidget);
  this->IFID = d->IFSelector->itemData(index).toString().toStdString();
  d->populateVOI(this->IFID);
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
    // "None" selected — ignore or reset state
    // d->VOIfileexcel->setText(QString::fromStdString(".xlsx"));
    return;
  }
  // std::cout << "Selected Patient ID: " << id << std::endl;
  // std :: string name = this->SubjectHierarchyNode->GetItemName(this->patID);
  // std :: string excelfile = name + "_Kinetics.xlsx";
  // d->VOIfileexcel->setText(QString::fromStdString(excelfile));

}

void qSlicerKMAPModuleWidget::onVOISelectAllbutton()
{
  Q_D(qSlicerKMAPModuleWidget);
  if (this->VOIsegmentIDs.size()==(d->VOICheckLayout->count()-1)) {
    for (int i = 0; i < d->VOICheckLayout->count(); ++i)
    {
      QWidget* widget = d->VOICheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->setChecked(false);
      }
    }
  } else {
    for (int i = 0; i < d->VOICheckLayout->count(); ++i)
    {
      QWidget* widget = d->VOICheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->setChecked(true);
      }
    }
  }
  d->populateVOI(this->IFID);
}


void qSlicerKMAPModuleWidget::clearFITdata() {
  Q_D(qSlicerKMAPModuleWidget);
  this->segmentMTGA.clear();
  this->segmentTCM.clear();
  d->populateResultsVOI();
  d->TCMResultsButton->setCollapsed(true);
  // d->resultsdirexcel->setCurrentPath(QString::fromStdString(""));
  // d->saveResultsExcelButton->setEnabled(false);
  // this->RemoveExistingPlotChartAndTable();
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


  if (IFID.empty()) {
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

  std::map<std::string, std::vector<std::vector<double>>> tac;
  for (const auto& [segmentName, statsVec] : segmentTACs)
  {
    if (statsVec.size() != static_cast<size_t>(Nframe))
    {
      std::cerr << "Mismatch in TAC frame size for segment " << segmentName << std::endl;
      return;
    }

    tac[segmentName].reserve(Nframe);
    for (const auto& vs : statsVec)
    {
      double value;
      if (currentSelectedStatID == "Mean")
        value = vs.mean;
      else if (currentSelectedStatID == "Median")
        value = vs.median;
      else
      {
        std::cerr << "Unknown stat: " << currentSelectedStatID << std::endl;
        return;
      }
      tac[segmentName].emplace_back(1, value);  // Adds one-element row (column vector)
    }
  }
  this->segmentTAC4TCMfits = tac;

  const double dk = 0.;
  const double timestep = 1.;
  const double pbrp[] = {1.0, 0.0, 0.0};
  const int maxiter = 100;

  std::vector<std::vector<double>> Cp = tac[IFID];

  for (const std::string& segmentID : VOIsegmentIDs)
  {
    const auto& tacVOI = tac[segmentID];

    for (const std::string& modelID : modelsID)
    {
      if (modelID == "Patlak") {
        std::cout << "Running Patlak" << std::endl;
        // logic->Patlak(...);
      }
      else if (modelID == "Logan") {
        std::cout << "Running Logan" << std::endl;
        // logic->Logan(...);
      }
      else if (modelID == "RE") {
        std::cout << "Running RE" << std::endl;
        // logic->RE(...);
      }
      else if (modelID == "1TCM") {
        bool sens[] = {true, true, true, false};
        double lb_1tcm[]   = {vbLower, k1Lower, k2Lower, 0.};
        double ub_1tcm[]   = {vbUpper, k1Upper, k2Upper, 0.};
        double init_1tcm[] = {vbInit,  k1Init,  k2Init,  0.};
        // std::cout << "Running 1TCM" << std::endl;

        logic->callTCM(tacVOI, Cp, framing, Nframe, Nvox,
                       init_1tcm, lb_1tcm, ub_1tcm, sens,
                       dk, timestep, pbrp, maxiter, 1,
                       this->segmentTCM[segmentID]["1TCM"]);
        logic->getFittedTCM(this->segmentTCMfits[segmentID]["1TCM"],
                            Cp, framing, Nframe, Nvox, init_1tcm, lb_1tcm,
                            ub_1tcm, sens, dk, timestep, pbrp, maxiter,
                            1, this->segmentTCM[segmentID]["1TCM"]);
      }
      else if (modelID == "1TdCM") {
        bool sens[] = {true, true, true, true};
        double lb_1tcm[]   = {vbLower, k1Lower, k2Lower, tdLower};
        double ub_1tcm[]   = {vbUpper, k1Upper, k2Upper, tdUpper};
        double init_1tcm[] = {vbInit,  k1Init,  k2Init,  tdInit};
        // std::cout << "Running 1TdCM" << std::endl;
        logic->callTCM(tacVOI, Cp, framing, Nframe, Nvox,
                       init_1tcm, lb_1tcm, ub_1tcm, sens,
                       dk, timestep, pbrp, maxiter, 1,
                       this->segmentTCM[segmentID]["1TdCM"]);
        logic->getFittedTCM(this->segmentTCMfits[segmentID]["1TdCM"],
                            Cp, framing, Nframe, Nvox, init_1tcm, lb_1tcm,
                            ub_1tcm, sens, dk, timestep, pbrp, maxiter,
                            1, this->segmentTCM[segmentID]["1TdCM"]);
      }
      else if (modelID == "1TiCM") {
        bool sens[] = {true, true, false, false};
        double lb_1tcm[]   = {vbLower, k1Lower, 0., 0.};
        double ub_1tcm[]   = {vbUpper, k1Upper, 0., 0.};
        double init_1tcm[] = {vbInit,  k1Init,  0.,  0.};
        // std::cout << "Running 1TiCM" << std::endl;
        logic->callTCM(tacVOI, Cp, framing, Nframe, Nvox,
                       init_1tcm, lb_1tcm, ub_1tcm, sens,
                       dk, timestep, pbrp, maxiter, 1,
                       this->segmentTCM[segmentID]["1TiCM"]);
        logic->getFittedTCM(this->segmentTCMfits[segmentID]["1TiCM"],
                            Cp, framing, Nframe, Nvox, init_1tcm, lb_1tcm,
                            ub_1tcm, sens, dk, timestep, pbrp, maxiter,
                            1, this->segmentTCM[segmentID]["1TiCM"]);
      }
      else if (modelID == "1TidCM") {
        bool sens[] = {true, true, false, true};
        double lb_1tcm[]   = {vbLower, k1Lower, 0., tdLower};
        double ub_1tcm[]   = {vbUpper, k1Upper, 0., tdUpper};
        double init_1tcm[] = {vbInit,  k1Init,  0.,  tdInit};
        // std::cout << "Running 1TidCM" << std::endl;
        logic->callTCM(tacVOI, Cp, framing, Nframe, Nvox,
                       init_1tcm, lb_1tcm, ub_1tcm, sens,
                       dk, timestep, pbrp, maxiter, 1,
                       this->segmentTCM[segmentID]["1TidCM"]);
        logic->getFittedTCM(this->segmentTCMfits[segmentID]["1TidCM"],
                            Cp, framing, Nframe, Nvox, init_1tcm, lb_1tcm,
                            ub_1tcm, sens, dk, timestep, pbrp, maxiter,
                            1, this->segmentTCM[segmentID]["1TidCM"]);
      }
      else if (modelID == "2TCM") {
        bool sens[] = {true, true, true, true, true, false};
        double lb_2tcm[]   = {vbLower, k1Lower, k2Lower, k3Lower, k4Lower, 0.};
        double ub_2tcm[]   = {vbUpper, k1Upper, k2Upper, k3Upper, k4Upper, 0.};
        double init_2tcm[] = {vbInit,  k1Init,  k2Init,  k3Init,  k4Init,  0.};
        // std::cout << "Running 2TCM" << std::endl;
        logic->callTCM(tacVOI, Cp, framing, Nframe, Nvox,
                       init_2tcm, lb_2tcm, ub_2tcm, sens,
                       dk, timestep, pbrp, maxiter, 2,
                       this->segmentTCM[segmentID]["2TCM"]);
        logic->getFittedTCM(this->segmentTCMfits[segmentID]["2TCM"],
                            Cp, framing, Nframe, Nvox, init_2tcm, lb_2tcm,
                            ub_2tcm, sens, dk, timestep, pbrp, maxiter,
                            1, this->segmentTCM[segmentID]["2TCM"]);
      }
      else if (modelID == "2dTCM") {
        bool sens[] = {true, true, true, true, true, true};
        double lb_2tcm[]   = {vbLower, k1Lower, k2Lower, k3Lower, k4Lower, tdLower};
        double ub_2tcm[]   = {vbUpper, k1Upper, k2Upper, k3Upper, k4Upper, tdUpper};
        double init_2tcm[] = {vbInit,  k1Init,  k2Init,  k3Init,  k4Init,  tdInit};
        // std::cout << "Running 2dTCM" << std::endl;
        logic->callTCM(tacVOI, Cp, framing, Nframe, Nvox,
                       init_2tcm, lb_2tcm, ub_2tcm, sens,
                       dk, timestep, pbrp, maxiter, 2,
                       this->segmentTCM[segmentID]["2dTCM"]);
        logic->getFittedTCM(this->segmentTCMfits[segmentID]["2dTCM"],
                            Cp, framing, Nframe, Nvox, init_2tcm, lb_2tcm,
                            ub_2tcm, sens, dk, timestep, pbrp, maxiter,
                            1, this->segmentTCM[segmentID]["2dTCM"]);
      }
      else if (modelID == "2TiCM") {
        bool sens[] = {true, true, true, true, false, false};
        double lb_2tcm[]   = {vbLower, k1Lower, k2Lower, k3Lower, 0., 0.};
        double ub_2tcm[]   = {vbUpper, k1Upper, k2Upper, k3Upper, 0., 0.};
        double init_2tcm[] = {vbInit,  k1Init,  k2Init,  k3Init,  0., 0.};
        // std::cout << "Running 2TiCM" << std::endl;
        logic->callTCM(tacVOI, Cp, framing, Nframe, Nvox,
                       init_2tcm, lb_2tcm, ub_2tcm, sens,
                       dk, timestep, pbrp, maxiter, 2,
                       this->segmentTCM[segmentID]["2TiCM"]);
        logic->getFittedTCM(this->segmentTCMfits[segmentID]["2TiCM"],
                            Cp, framing, Nframe, Nvox, init_2tcm, lb_2tcm,
                            ub_2tcm, sens, dk, timestep, pbrp, maxiter,
                            1, this->segmentTCM[segmentID]["2TiCM"]);
      }
      else if (modelID == "2TidCM") {
        bool sens[] = {true, true, true, true, false, true};
        double lb_2tcm[]   = {vbLower, k1Lower, k2Lower, k3Lower, 0., tdLower};
        double ub_2tcm[]   = {vbUpper, k1Upper, k2Upper, k3Upper, 0., tdUpper};
        double init_2tcm[] = {vbInit,  k1Init,  k2Init,  k3Init,  0.,  tdInit};
        // std::cout << "Running 2TidCM" << std::endl;
        logic->callTCM(tacVOI, Cp, framing, Nframe, Nvox,
                       init_2tcm, lb_2tcm, ub_2tcm, sens,
                       dk, timestep, pbrp, maxiter, 2,
                       this->segmentTCM[segmentID]["2TidCM"]);
        logic->getFittedTCM(this->segmentTCMfits[segmentID]["2TidCM"],
                            Cp, framing, Nframe, Nvox, init_2tcm, lb_2tcm,
                            ub_2tcm, sens, dk, timestep, pbrp, maxiter,
                            1, this->segmentTCM[segmentID]["2TidCM"]);
      }
      else {
        std::cerr << "Unknown model ID: " << modelID << std::endl;
        return;
      }
    }
  }
  d->populateResultsVOI();
}

void qSlicerKMAPModuleWidget::onModelsAllbutton()
{
  Q_D(qSlicerKMAPModuleWidget);

  if (this->modelsID.size()==(d->ModelsCheckLayout->count())) {
    this->modelsID.clear();
    for (int i = 0; i < d->ModelsCheckLayout->count(); ++i)
    {
      QLayoutItem* item = d->ModelsCheckLayout->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (checkbox)
      {
        checkbox->setChecked(false);
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
        checkbox->setChecked(true);
        this->modelsID.push_back(checkbox->text().toStdString());
      }
    }
  }
  this->enableFITbutton();

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

void qSlicerKMAPModuleWidget::onModelsTCMSelectAllbutton()
{
  Q_D(qSlicerKMAPModuleWidget);
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
        cb->setChecked(false);
      }
    }
  } else {
    for (int i = 0; i < d->ModelsTCMCheckLayout->count(); ++i)
    {
      QWidget* widget = d->ModelsTCMCheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->setChecked(true);
      }
    }
  }
  d->populateModelsTCM(this->plotTCMVOI);
}


void qSlicerKMAPModuleWidget::onPlotTCMbutton()
{
  Q_D(qSlicerKMAPModuleWidget);

  vtkMRMLScene* scene = this->mrmlScene();

  // Get selected VOI to plot TCM fit
  std :: string selectedVOI = this->plotTCMVOI;

  // Get TAC to plot
  std::vector<std::vector<double>> tacvoi = this->segmentTAC4TCMfits[selectedVOI];

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
  for (size_t i = 0; i < tacvoi.size(); ++i)
  {
      // Assuming tacvoi[0] holds measured values for VOI (adapt if structured differently)
      tacArray->InsertNextValue(tacvoi[i][0]);
  }
  tableNode->AddColumn(tacArray);

  vtkSmartPointer<vtkMRMLPlotSeriesNode> scatterSeries = vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
  scene->AddNode(scatterSeries);
  scatterSeries->SetName("TAC");
  scatterSeries->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);  // scatter points
  scatterSeries->SetAndObserveTableNodeID(tableNode->GetID());
  scatterSeries->SetXColumnName("Time (min)");
  scatterSeries->SetYColumnName("TAC");
  scatterSeries->SetUniqueColor();
  scatterSeries->SetLineStyle(vtkMRMLPlotSeriesNode::LineStyleNone);
  chartNode->AddAndObservePlotSeriesNodeID(scatterSeries->GetID());
  chartNode->SetTitle(this->segmentTACsnames[selectedVOI].c_str());
  chartNode->SetXAxisTitle("Time (min)");
  chartNode->SetYAxisTitle("SUVbw (g/mL)");

  // TCM fits as line plots
  for (const std::string& modelName : PlotSelectedTCMs)
  {
      auto it = TCMfits.find(modelName);
      if (it == TCMfits.end() || !it->second)
        continue;


      double* fitArrayPtr = it->second;

      vtkNew<vtkDoubleArray> fitArray;
      fitArray->SetName(modelName.c_str());

      for (size_t i = 0; i < this->timePoints.size(); ++i)
      {
          fitArray->InsertNextValue(fitArrayPtr[i]);
      }

      tableNode->AddColumn(fitArray);

      vtkSmartPointer<vtkMRMLPlotSeriesNode> lineSeries = vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
      scene->AddNode(lineSeries);
      lineSeries->SetName(modelName.c_str());
      lineSeries->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);  // line plot
      lineSeries->SetAndObserveTableNodeID(tableNode->GetID());
      lineSeries->SetXColumnName("Time (min)");
      lineSeries->SetYColumnName(modelName.c_str());
      lineSeries->SetUniqueColor();
      lineSeries->SetMarkerStyle(vtkMRMLPlotSeriesNode::MarkerStyleNone);
      chartNode->AddAndObservePlotSeriesNodeID(lineSeries->GetID());
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
  }

}
