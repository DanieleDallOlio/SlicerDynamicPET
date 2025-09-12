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

  auto* segmentation = vtkSegmentation::SafeDownCast(caller);
  if (segmentation)
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

  if (self->GetCurrentSegID() != segNode->GetName())
    return;

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

  // Make sure the
  segmentation->SeparateSegmentLabelmap(segmentId);

  vtkSegment* segment = segmentation->GetSegment(segmentId);
  if (!segment || segment->GetRepresentation("Binary labelmap") == nullptr)
  {
    vtkGenericWarningMacro("Segment " << segmentId
                       << " is empty at frame " << frameIndex);
    (*segmentTACs)[segmentId][frameIndex] = VoxelStatistics{};
    (*segmentTACs)[segmentId][frameIndex].keep = false;
    (*segmentTACs)[segmentId][frameIndex].empty = true;
    if (self->RunPlot)
    {
      self->RunPlot();  // calls your widget function safely
    }
    return;
  }

  vtkOrientedImageData* segLabelmap = vtkOrientedImageData::SafeDownCast(
      segment->GetRepresentation(vtkSegmentationConverter::GetSegmentationBinaryLabelmapRepresentationName()));

  if (!segLabelmap)
  {
    vtkGenericWarningMacro("No binary labelmap for " << segmentId << " at frame " << frameIndex);
    (*segmentTACs)[segmentId][frameIndex] = VoxelStatistics{};
    (*segmentTACs)[segmentId][frameIndex].keep = false;
    (*segmentTACs)[segmentId][frameIndex].empty = true;
    if (self->RunPlot)
    {
      self->RunPlot();  // calls your widget function safely
    }
    return;
  }
  double range[2];
  segLabelmap->GetScalarRange(range);
  if (range[1] == 0.0) // max is zero → fully empty
  {
    vtkGenericWarningMacro("Segment " << segmentId << " is empty at frame " << frameIndex);
    (*segmentTACs)[segmentId][frameIndex] = VoxelStatistics{};
    (*segmentTACs)[segmentId][frameIndex].keep = false;
    (*segmentTACs)[segmentId][frameIndex].empty = true;
    if (self->RunPlot)
    {
      self->RunPlot();  // calls your widget function safely
    }
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

  int labelOrientedImageDataEffectiveExtent[6] = { 0, -1, 0, -1, 0, -1 };
  if (!vtkOrientedImageDataResample::CalculateEffectiveExtent(labelmap, labelOrientedImageDataEffectiveExtent))
  {
    vtkGenericWarningMacro("Segment " << segmentId << " has been removed at frame " << frameIndex);
    stats = VoxelStatistics{};
    stats.keep = false;
    stats.empty = true;
  } else {
    stats = logic->ComputeVoxelStatistics(PETVolume, labelmap, 1);
  }

  (*segmentTACs)[segmentId][frameIndex] = stats;

  if (self->RunPlot)
  {
    self->RunPlot();  // calls your widget function safely
  }

}
