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
  // QObject::connect( this->PlotErrorCheckbox, SIGNAL(toggled(bool)),
  //   q, SLOT(onPlotbutton()));

  this->TACCollapsibleButton->setCollapsed(true);
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
  std :: string excelfile = name + ".xlsx";
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
  this->ProgressBar->setVisible(false);
  qApp->processEvents();
  d->populatePlotSegmentCheckboxes();
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
