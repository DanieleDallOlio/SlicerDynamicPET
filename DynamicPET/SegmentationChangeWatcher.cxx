// SegmentationChangeWatcher.cxx
#include "SegmentationChangeWatcher.h"

#include <vtkMRMLScene.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkOrientedImageData.h>
#include <vtkOrientedImageDataResample.h>
#include <vtkSegmentationConverter.h>
#include <vtkSlicerSegmentationsModuleLogic.h>
#include <vtkStringArray.h>

vtkStandardNewMacro(SegmentationChangeWatcher);

SegmentationChangeWatcher::~SegmentationChangeWatcher()
{
  this->Clear();
}

void SegmentationChangeWatcher::Clear()
{
  for (auto& kv : this->SegTags)
  {
    vtkSegmentation* seg = kv.first;
    if (!seg)
    {
      continue;
    }

    for (unsigned long tag : kv.second)
    {
      seg->RemoveObserver(tag);
    }
  }
  this->SegTags.clear();

  for (auto* node : this->ObservedNodes)
  {
    if (node)
    {
      node->RemoveObserver(this->Callback);
    }
  }
  this->ObservedNodes.clear();
  this->SegToNode.clear();
}

SegmentationChangeWatcher::SegmentationChangeWatcher()
{
  this->Callback = vtkSmartPointer<vtkCallbackCommand>::New();
  this->Callback->SetClientData(this);
  this->Callback->SetCallback(SegmentationChangeWatcher::OnSegmentationChanged);
}

void SegmentationChangeWatcher::ObserveSegmentationNode(vtkMRMLSegmentationNode* segNode)
{
  if (!segNode)
  {
    return;
  }

  if (this->ObservedNodes.find(segNode) == this->ObservedNodes.end())
  {
    segNode->AddObserver(
        vtkMRMLSegmentationNode::SegmentationChangedEvent,
        this->Callback);
    this->ObservedNodes.insert(segNode);
  }

  vtkSegmentation* seg = segNode->GetSegmentation();
  if (!seg)
  {
    return;
  }

  this->SegToNode[seg] = segNode;

  if (this->SegTags.find(seg) == this->SegTags.end())
  {
    std::vector<unsigned long> tags;
    tags.push_back(seg->AddObserver(
        vtkSegmentation::RepresentationModified,
        this->Callback));
    tags.push_back(seg->AddObserver(
        vtkSegmentation::SegmentAdded,
        this->Callback));
    tags.push_back(seg->AddObserver(
        vtkSegmentation::SegmentRemoved,
        this->Callback));
    this->SegTags[seg] = std::move(tags);
  }
}

void SegmentationChangeWatcher::SynchronizeSegmentAdded(
    vtkMRMLSegmentationNode* proxyNode,
    const std::string& segmentId)
{
  if (!proxyNode || segmentId.empty() || !this->GetSequenceSeg)
  {
    return;
  }

  vtkSegmentation* sourceSegmentation = proxyNode->GetSegmentation();
  if (!sourceSegmentation || !sourceSegmentation->GetSegment(segmentId))
  {
    return;
  }

  vtkMRMLSequenceNode* sequenceSeg = this->GetSequenceSeg();
  if (!sequenceSeg)
  {
    return;
  }

  const int numberOfFrames = sequenceSeg->GetNumberOfDataNodes();

  for (int i = 0; i < numberOfFrames; ++i)
  {
    vtkMRMLSegmentationNode* frameNode =
        vtkMRMLSegmentationNode::SafeDownCast(
            sequenceSeg->GetNthDataNode(i));

    if (!frameNode || !frameNode->GetSegmentation())
    {
      continue;
    }

    vtkSegmentation* frameSegmentation = frameNode->GetSegmentation();
    if (frameSegmentation->GetSegment(segmentId))
    {
      continue;
    }

    if (!frameSegmentation->CopySegmentFromSegmentation(
            sourceSegmentation,
            segmentId,
            false))
    {
      vtkGenericWarningMacro(
          "Could not propagate newly added segment "
          << segmentId << " to segmentation sequence frame " << i);
    }
  }
}

void SegmentationChangeWatcher::SynchronizeSegmentRemoved(
    const std::string& segmentId)
{
  if (segmentId.empty() || !this->GetSequenceSeg)
  {
    return;
  }

  vtkMRMLSequenceNode* sequenceSeg = this->GetSequenceSeg();
  if (!sequenceSeg)
  {
    return;
  }

  const int numberOfFrames = sequenceSeg->GetNumberOfDataNodes();

  for (int i = 0; i < numberOfFrames; ++i)
  {
    vtkMRMLSegmentationNode* frameNode =
        vtkMRMLSegmentationNode::SafeDownCast(
            sequenceSeg->GetNthDataNode(i));

    if (!frameNode || !frameNode->GetSegmentation())
    {
      continue;
    }

    vtkSegmentation* frameSegmentation = frameNode->GetSegmentation();
    if (frameSegmentation->GetSegment(segmentId))
    {
      frameSegmentation->RemoveSegment(segmentId);
    }
  }
}

void SegmentationChangeWatcher::OnSegmentationChanged(
    vtkObject* caller,
    unsigned long eid,
    void* clientData,
    void* callData)
{
  auto* self = static_cast<SegmentationChangeWatcher*>(clientData);
  if (!self)
  {
    return;
  }

  auto* segmentation = vtkSegmentation::SafeDownCast(caller);

  if (!segmentation)
  {
    if (auto* node = vtkMRMLSegmentationNode::SafeDownCast(caller))
    {
      // The proxy received a different vtkSegmentation object. Reattach
      // segmentation-level observers to the new object.
      self->ObserveSegmentationNode(node);
    }
    return;
  }

  auto ownerIt = self->SegToNode.find(segmentation);
  if (ownerIt == self->SegToNode.end())
  {
    return;
  }

  vtkMRMLSegmentationNode* segNode = ownerIt->second.GetPointer();
  if (!segNode)
  {
    return;
  }

  std::string segmentId;
  if (callData)
  {
    segmentId = static_cast<const char*>(callData);
  }

  if (segmentId.empty())
  {
    return;
  }

  // Segment membership is a sequence-wide structural property. When a segment
  // is moved/copied into the active proxy, Slicer only edits the current proxy
  // item. Propagate the segment structure across all dynamic segmentation
  // frames so later TAC extraction sees a consistent segment list.
  if (eid == vtkSegmentation::SegmentAdded)
  {
    self->SynchronizeSegmentAdded(segNode, segmentId);
    if (self->OnSegmentStructureChanged)
    {
      self->OnSegmentStructureChanged();
    }
    return;
  }

  if (eid == vtkSegmentation::SegmentRemoved)
  {
    self->SynchronizeSegmentRemoved(segmentId);
    if (self->OnSegmentStructureChanged)
    {
      self->OnSegmentStructureChanged();
    }
    return;
  }

  if (eid != vtkSegmentation::RepresentationModified)
  {
    return;
  }

  auto* segmentTACs = self->GetsegmentTACs ? self->GetsegmentTACs() : nullptr;
  if (!segmentTACs || segmentTACs->empty())
  {
    return;
  }

  // Only refresh a segment that already has a computed TAC. Newly added
  // segments are handled by the normal TAC button after structure sync.
  auto tacIt = segmentTACs->find(segmentId);
  if (tacIt == segmentTACs->end())
  {
    return;
  }

  int frameIndex = -1;
  if (self->browser)
  {
    frameIndex = self->browser->GetSelectedItemNumber();
  }
  if (frameIndex < 0 ||
      frameIndex >= static_cast<int>(tacIt->second.size()))
  {
    return;
  }

  vtkMRMLSequenceNode* PETSequenceNode =
      self->GetSequencePET ? self->GetSequencePET() : nullptr;
  if (!PETSequenceNode || frameIndex >= PETSequenceNode->GetNumberOfDataNodes())
  {
    return;
  }

  vtkMRMLScalarVolumeNode* PETVolume =
      vtkMRMLScalarVolumeNode::SafeDownCast(
          PETSequenceNode->GetNthDataNode(frameIndex));
  if (!PETVolume)
  {
    return;
  }

  vtkSegment* segment = segmentation->GetSegment(segmentId);
  if (!segment)
  {
    VoxelStatistics stats;
    stats.keep = false;
    stats.empty = true;
    tacIt->second[frameIndex] = stats;

    if (self->OnSegmentTACChanged)
    {
      self->OnSegmentTACChanged(segmentId);
    }
    return;
  }

  vtkNew<vtkStringArray> segmentArray;
  segmentArray->InsertNextValue(segmentId);

  vtkSmartPointer<vtkOrientedImageData> labelmap =
      vtkSmartPointer<vtkOrientedImageData>::New();

  vtkSlicerSegmentationsModuleLogic::GenerateMergedLabelmapInReferenceGeometry(
      segNode,
      PETVolume,
      segmentArray,
      vtkSegmentation::EXTENT_UNION_OF_EFFECTIVE_SEGMENTS,
      labelmap);

  VoxelStatistics stats;

  if (!labelmap ||
      !labelmap->GetPointData() ||
      !labelmap->GetPointData()->GetScalars())
  {
    stats.keep = false;
    stats.empty = true;
  }
  else
  {
    auto* logic = self->GetLogic ? self->GetLogic() : nullptr;
    if (!logic)
    {
      return;
    }
    stats = logic->ComputeVoxelStatistics(PETVolume, labelmap, 1);
  }

  tacIt->second[frameIndex] = stats;

  if (self->OnSegmentTACChanged)
  {
    self->OnSegmentTACChanged(segmentId);
  }

  // Plot refresh remains a user-selectable performance option. The TAC cache
  // itself is refreshed regardless so IF/ROI data cannot silently stay stale.
  if (self->GetSegEditCorr && self->GetSegEditCorr() && self->RunPlot)
  {
    self->RunPlot();
  }
}
