// SegmentationChangeWatcher.cxx
#include "SegmentationChangeWatcher.h"

#include <vtkSegmentation.h>
#include <vtkMRMLScene.h>
#include <vtkObjectFactory.h>

vtkStandardNewMacro(SegmentationChangeWatcher);

SegmentationChangeWatcher::SegmentationChangeWatcher()
{
  this->Callback = vtkSmartPointer<vtkCallbackCommand>::New();
  this->Callback->SetClientData(this);
  this->Callback->SetCallback(SegmentationChangeWatcher::OnSegmentationChanged);
}

void SegmentationChangeWatcher::ObserveSegmentationNode(vtkMRMLSegmentationNode* segNode)
{
  if (!segNode) return;
  // Only add the observer once
  if (ObservedNodes.find(segNode) == ObservedNodes.end())
  {
      segNode->AddObserver(vtkMRMLSegmentationNode::SegmentationChangedEvent, this->Callback);
      ObservedNodes.insert(segNode);
  }

  // Skip if segmentation is null
  if (!segNode->GetSegmentation())
    return;

  vtkSegmentation* seg = segNode->GetSegmentation();
  this->SegToNode[seg] = segNode;

  // You can also optionally track per-segmentation to avoid double observers
  if (SegTags.find(seg) == SegTags.end())
  {
      unsigned long tag = seg->AddObserver(vtkSegmentation::RepresentationModified, this->Callback);
      SegTags[seg] = tag;
  }
}

void SegmentationChangeWatcher::OnSegmentationChanged(
  vtkObject* caller, unsigned long eid, void* clientData, void* callData)
{
  auto* self = static_cast<SegmentationChangeWatcher*>(clientData);
  if (!self->GetSegEditCorr()) {
    return;
  }

  auto* segmentTACs = self->GetsegmentTACs();
  if (segmentTACs == nullptr)
    return;
  if (segmentTACs->empty())
    return;

  vtkMRMLSegmentationNode* segNode = nullptr;
  std::string segmentId;

  if (auto* segmentation = vtkSegmentation::SafeDownCast(caller))
  {
      segNode = self->SegToNode[segmentation];
      if (callData)
          segmentId = static_cast<const char*>(callData);
  }
  else if (auto* node = vtkMRMLSegmentationNode::SafeDownCast(caller))
  {
      // node swapped segmentation, re-attach
      self->ObserveSegmentationNode(node);
      return;
  }

  if (!segNode || segmentId.empty())
      return;

  // Get the current frame index in the sequence
  int frameIndex = -1;
  if (self->browser)
  {
      frameIndex = self->browser->GetSelectedItemNumber();
  }
  if (frameIndex < 0)
      return;

  // Retrieve the corresponding PET volume
  vtkMRMLSequenceNode* PETSequenceNode = self->GetSequencePET();
  auto* PETVolume = vtkMRMLScalarVolumeNode::SafeDownCast(
      PETSequenceNode->GetDataNodeAtValue(
          PETSequenceNode->GetNthIndexValue(frameIndex)));

  if (!PETVolume)
      return;

  // Check if segment exists
  auto it = segmentTACs->find(segmentId);
  if (it == segmentTACs->end())
  {
    vtkGenericWarningMacro("Segment ID " << segmentId << " TAC yet to be available.");
    return;
  }

  // Check if frameIndex is in range
  if (frameIndex < 0 || frameIndex >= static_cast<int>(it->second.size()))
  {
    vtkGenericWarningMacro("Frame index " << frameIndex << " out of range for segment " << segmentId);
    return;
  }

  vtkNew<vtkStringArray> segmentArray;
  segmentArray->InsertNextValue(segmentId);

  auto* logic = self->GetLogic();
  vtkSmartPointer<vtkOrientedImageData> labelmap = vtkSmartPointer<vtkOrientedImageData>::New();
  vtkSlicerSegmentationsModuleLogic::GenerateMergedLabelmapInReferenceGeometry(segNode,
                                                                               PETVolume,
                                                                               segmentArray,
                                                                               vtkSegmentation::EXTENT_UNION_OF_EFFECTIVE_SEGMENTS,
                                                                               labelmap);

  VoxelStatistics stats;
  if (labelmap && labelmap->GetNumberOfPoints() > 0) {
    stats = logic->ComputeVoxelStatistics(PETVolume, labelmap, 1);
  } else {
    vtkGenericWarningMacro("Segment " << segmentId << " has been removed at frame " << frameIndex);
    stats = VoxelStatistics{};
  }

  (*segmentTACs)[segmentId][frameIndex] = stats;

  if (self->RunPlot)
  {
    self->RunPlot();  // calls your widget function safely
  }

}
