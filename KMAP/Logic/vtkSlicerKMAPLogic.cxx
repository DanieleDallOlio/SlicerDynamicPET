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

// KMAP Logic includes
#include "vtkSlicerKMAPLogic.h"

// MRML includes
#include <vtkMRMLScene.h>

// VTK includes
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>

// STD includes
#include <cassert>

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerKMAPLogic);

//----------------------------------------------------------------------------
vtkSlicerKMAPLogic::vtkSlicerKMAPLogic()
{
}

//----------------------------------------------------------------------------
vtkSlicerKMAPLogic::~vtkSlicerKMAPLogic()
{
}

//----------------------------------------------------------------------------
void vtkSlicerKMAPLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
void vtkSlicerKMAPLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//-----------------------------------------------------------------------------
void vtkSlicerKMAPLogic::RegisterNodes()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerKMAPLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerKMAPLogic
::OnMRMLSceneNodeAdded(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerKMAPLogic
::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
{
}


void vtkSlicerKMAPLogic::computeTAC(vtkIdType ctID, vtkIdType petID, vtkIdType segID, std::vector<QString> segmentsID)
{
  vtkMRMLScene* scene = this->GetMRMLScene();
  if (scene==nullptr) {
    return;
  }
  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (!shNode) {
    return;
  }
  // Fetch CT
  vtkMRMLScalarVolumeNode* ctNode = vtkMRMLScalarVolumeNode::SafeDownCast(shNode->GetItemDataNode(ctID));
  if (!ctNode) {
    return;
  }
  std :: cout << ctID << std :: endl;
  std :: string name = shNode->GetItemName(ctID);
  std::cout << "CT: " << name
            << ", ID: " << ctID << std::endl;
  // Fetch PET
  vtkMRMLScalarVolumeNode* petNode = vtkMRMLScalarVolumeNode::SafeDownCast(shNode->GetItemDataNode(petID));
  if (!petNode) {
    return;
  }
  name = shNode->GetItemName(petID);
  std::cout << "PET: " << name
            << ", ID: " << petID << std::endl;
  // Fetch Segmentation
  vtkMRMLSegmentationNode* segNode = vtkMRMLSegmentationNode::SafeDownCast(shNode->GetItemDataNode(segID));
  if (!segNode) {
    return;
  }
  name = shNode->GetItemName(segID);
  std::cout << "SEG: " << name
            << ", ID: " << segID << std::endl;
  // Fetch segments
  vtkSegmentation* seg = segNode->GetSegmentation();
  if (!seg) {
    return;
  }
  if (segmentsID.empty()) {
    return;
  }
  std::cout << "Selected segments:" << std::endl;
  for (const QString& id : segmentsID)
  {
    std::string segmentName = seg->GetSegment(id.toStdString())->GetName();
    std::cout << segmentName << " - " << id.toStdString() << std::endl;
  }
  // Temporary print-out
  std :: cout << "Ready for next step!" << std :: endl;
}
