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


}

void qSlicerKMAPModuleWidgetPrivate::populatePatientComboBox() {
  Q_Q(qSlicerKMAPModuleWidget);

  vtkIdType currentSelectedID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  int currentIndex = this->PatSelector->currentIndex();
  std :: cout << "currentIndex=" << currentIndex << std :: endl;
  if (currentIndex >= 0)
  {
    currentSelectedID = this->PatSelector->itemData(currentIndex).value<vtkIdType>();
  }

  this->PatSelector->blockSignals(true);  // Optional: prevent signal emission
  this->PatSelector -> clear();
  this->PatSelector->addItem(QString::fromStdString("None"), QVariant::fromValue(vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID));


  vtkMRMLScene* scene = q->mrmlScene();

  if (scene==nullptr)
    return ;

  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  vtkIdType rootID = shNode->GetSceneItemID();
  if (rootID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
    return ;

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
    std ::cout << "populatePatientComboBox" << std :: endl;
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

  std :: cout << "restoredIndex=" <<restoredIndex<<std::endl;
  this->populateStudyComboBox(restoredIndex);

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
    this->populateNodeComboBox(this->CTSelector,
                               vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID,
                               "vtkMRMLScalarVolumeNode",
                               "CT"
                              );
    this->populateNodeComboBox(this->PETSelector,
                               vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID,
                               "vtkMRMLScalarVolumeNode",
                               "PT"
                              );
    this->populateNodeComboBox(this->SegSelector,
                               vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID,
                               "vtkMRMLSegmentationNode",
                               ""
                              );
    return;
  }
  this->StuSelector->setEnabled(true);


  vtkMRMLScene* scene = q->mrmlScene();
  if (!scene)
  {
    this->StuSelector->blockSignals(false);
    return;
  }

  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (!shNode)
  {
    this->StuSelector->blockSignals(false);
    return;
  }

  // Retrieve the children of the given patient
  std::vector<vtkIdType> children;
  std ::cout << "populateStudyComboBox - ID:" << patientID << std :: endl;
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

  this->populateNodeComboBox(this->CTSelector,
                             restoredIndex,
                             "vtkMRMLScalarVolumeNode",
                             "CT"
                            );
  this->populateNodeComboBox(this->PETSelector,
                             restoredIndex,
                             "vtkMRMLScalarVolumeNode",
                             "PT"
                            );
  this->populateNodeComboBox(this->SegSelector,
                             restoredIndex,
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
      this->populateSegmentCheckboxes(vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID);
    }
    return;
  }

  vtkMRMLScene* scene = q->mrmlScene();
  if (!scene)
  {
    comboBox->blockSignals(false);
    return;
  }

  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (!shNode)
  {
    comboBox->blockSignals(false);
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
    std ::cout << "populateNodeComboBox" << std :: endl;
    shNode->GetItemChildren(itemID, children);
    for (vtkIdType childID : children)
    {
      collectItems(childID);
    }
  };

  collectItems(parentItemID);
  comboBox->setCurrentIndex(restoredIndex);
  comboBox->blockSignals(false);

  if (std::string(requiredNodeType)=="vtkMRMLSegmentationNode") {
    this->populateSegmentCheckboxes(restoredIndex);
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
    this->TACbutton->setEnabled(false);
    return;
  }

  vtkMRMLScene* scene = q->mrmlScene();
  if (!scene)
  {
    this->SegmentCheckContents->blockSignals(false);
    this->TACbutton->setEnabled(false);
    return;
  }

  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (!shNode)
  {
    this->SegmentCheckContents->blockSignals(false);
    this->TACbutton->setEnabled(false);
    return;
  }

  vtkMRMLSegmentationNode* segNode = vtkMRMLSegmentationNode::SafeDownCast(shNode->GetItemDataNode(SegItemID));
  // Step 3: Repopulate based on new segmentation node
  if (!segNode) {
    this->SegmentCheckContents->blockSignals(false);
    this->TACbutton->setEnabled(false);
    return;
  }

  vtkSegmentation* segmentation = segNode->GetSegmentation();
  if (!segmentation) {
    this->SegmentCheckContents->blockSignals(false);
    this->TACbutton->setEnabled(false);
    return;
  }

  // this->segmentCheckLayout->setEnabled(true);

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
  }
  if (previouslySelectedIDs.empty()) {
    this->TACbutton->setEnabled(false);
  } else {
    this->TACbutton->setEnabled(false);
  }

  this->segmentCheckLayout->addStretch();
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
  this -> SubjectHierarchyNode = nullptr;
  this -> patID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  this -> stuID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  this -> ctID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  this -> petID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  this -> segID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
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


void qSlicerKMAPModuleWidget::setMRMLScene(vtkMRMLScene* scene) {
  this->Superclass::setMRMLScene(scene);
  Q_D(qSlicerKMAPModuleWidget);

  this->qvtkDisconnectAll();

  this -> SubjectHierarchyNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (this -> SubjectHierarchyNode)
  {
    this->qvtkConnect(this -> SubjectHierarchyNode, vtkMRMLSubjectHierarchyNode::SubjectHierarchyItemAddedEvent,
                      this, SLOT(onSubjectHierarchyChanged()));
    this->qvtkConnect(this -> SubjectHierarchyNode, vtkMRMLSubjectHierarchyNode::SubjectHierarchyItemRemovedEvent,
                      this, SLOT(onSubjectHierarchyChanged()));
    this->qvtkConnect(this -> SubjectHierarchyNode, vtkMRMLSubjectHierarchyNode::SubjectHierarchyItemModifiedEvent,
                      this, SLOT(onSubjectHierarchyChanged()));
  }

  d->populatePatientComboBox();
}

void qSlicerKMAPModuleWidget::onSubjectHierarchyChanged() {
  Q_D(qSlicerKMAPModuleWidget);
  d->populatePatientComboBox();
}

void qSlicerKMAPModuleWidget::onPatChanged (int index) {
  Q_D(qSlicerKMAPModuleWidget);
  if (!this->mrmlScene() || this->mrmlScene()->IsBatchProcessing())
  {
    return;
  }
  this->SubjectHierarchyNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(this->mrmlScene());

  if (this->SubjectHierarchyNode== nullptr) {
    return;
  }


  this -> patID = d->PatSelector->itemData(index).value<vtkIdType>();
  std :: cout << "onPatChanged - ID:" << this -> patID << std :: endl;
  d->populateStudyComboBox(this -> patID);
  if (this -> patID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
    // "None" selected — ignore or reset state
    return;
  }

  // std::cout << "Selected Patient ID: " << id << std::endl;
  // std :: string name = this->SubjectHierarchyNode->GetItemName(id);
  //
  // std::cout << "Name: " << name
  //           << ", ID: " << id << std::endl;
}


void qSlicerKMAPModuleWidget::onStuChanged (int index) {
  Q_D(qSlicerKMAPModuleWidget);
  if (!this->mrmlScene() || this->mrmlScene()->IsBatchProcessing())
  {
    return;
  }
  this->SubjectHierarchyNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(this->mrmlScene());

  if (this->SubjectHierarchyNode== nullptr) {
    return;
  }

  this -> stuID = d->StuSelector->itemData(index).value<vtkIdType>();
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
  if (this -> stuID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
    // "None" selected — ignore or reset state
    return;
  }

  // std::cout << "Selected Study ID: " << id << std::endl;
  // std :: string name = this->SubjectHierarchyNode->GetItemName(id);
  //
  // std::cout << "Study: " << name
  //           << ", ID: " << id << std::endl;

}


void qSlicerKMAPModuleWidget::onCTChanged (int index) {
  Q_D(qSlicerKMAPModuleWidget);
  std::cout << "ehi1" << std ::endl;
  if (!this->mrmlScene() || this->mrmlScene()->IsBatchProcessing())
  {
    return;
  }
  this->SubjectHierarchyNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(this->mrmlScene());

  if (this->SubjectHierarchyNode== nullptr) {
    return;
  }

  this -> ctID = d->CTSelector->itemData(index).value<vtkIdType>();
  if (this -> ctID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
    // "None" selected — ignore or reset state
    return;
  }

  std :: string name = this->SubjectHierarchyNode->GetItemName(this -> ctID);
  std::cout << "CT: " << name
            << ", ID: " << this -> ctID << std::endl;
}


void qSlicerKMAPModuleWidget::onPETChanged (int index) {
  Q_D(qSlicerKMAPModuleWidget);
  std::cout << "ehi2" << std ::endl;
  if (!this->mrmlScene() || this->mrmlScene()->IsBatchProcessing())
  {
    return;
  }
  this->SubjectHierarchyNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(this->mrmlScene());

  if (this->SubjectHierarchyNode== nullptr) {
    return;
  }

  this -> petID = d->PETSelector->itemData(index).value<vtkIdType>();
  if (this -> petID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
    // "None" selected — ignore or reset state
    return;
  }

  std :: string name = this->SubjectHierarchyNode->GetItemName(this -> petID);
  std::cout << "PET: " << name
            << ", ID: " << this -> petID << std::endl;
}

void qSlicerKMAPModuleWidget::onSegChanged (int index) {
  Q_D(qSlicerKMAPModuleWidget);
  std::cout << "ehi3" << std ::endl;
  if (!this->mrmlScene() || this->mrmlScene()->IsBatchProcessing())
  {
    return;
  }
  this->SubjectHierarchyNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(this->mrmlScene());

  if (this->SubjectHierarchyNode== nullptr) {
    return;
  }

  this -> segID = d->SegSelector->itemData(index).value<vtkIdType>();
  d->populateSegmentCheckboxes(this -> segID);
  if (this -> segID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
    // "None" selected — ignore or reset state
    return;
  }

  std :: string name = this->SubjectHierarchyNode->GetItemName(this -> segID);
  std::cout << "Seg: " << name
            << ", ID: " << this -> segID << std::endl;
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
  this -> segmentIDs = selectedSegmentIDs;

  if (segmentIDs.empty()) {
    d->TACbutton->setEnabled(false);
    return;
  }

  vtkMRMLScene* scene = this->mrmlScene();
  if (!scene)
    return;
  this->SubjectHierarchyNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (!this->SubjectHierarchyNode)
    return;
  vtkMRMLSegmentationNode* segNode = vtkMRMLSegmentationNode::SafeDownCast(this->SubjectHierarchyNode->GetItemDataNode(this->segID));
  if (!segNode) {
    return;
  }

  vtkSegmentation* segmentation = segNode->GetSegmentation();
  if (!segmentation) {
    return;
  }
  // Print to console
  std::cout << "Selected segments:" << std::endl;

  for (const QString& id : selectedSegmentIDs)
  {
    std::string segmentName = segmentation->GetSegment(id.toStdString())->GetName();
    std::cout << segmentName << " - " << id.toStdString() << std::endl;
  }

  d->TACbutton->setEnabled(true);
}

void qSlicerKMAPModuleWidget::onTACbutton()
{
  std :: cout << "Get TAC!" << std :: endl;
  return;
}
