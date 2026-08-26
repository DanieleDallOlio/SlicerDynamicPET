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

#include "SegmentationChangeWatcher.h"

#include <vtkMRMLScene.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkOrientedImageData.h>
#include <vtkSegmentationConverter.h>
#include <vtkSlicerSegmentationsModuleLogic.h>
#include <vtkStringArray.h>
#include <vtkSegment.h>

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
  this->NodeToSeg.clear();
  this->NodeFrameIndex.clear();
  this->KnownSegmentNames.clear();
  this->SegmentationSequence = nullptr;
  this->StructureSyncInProgress = false;
}

SegmentationChangeWatcher::SegmentationChangeWatcher()
{
  this->Callback = vtkSmartPointer<vtkCallbackCommand>::New();
  this->Callback->SetClientData(this);
  this->Callback->SetCallback(SegmentationChangeWatcher::OnSegmentationChanged);
}

void SegmentationChangeWatcher::RefreshKnownSegmentNames(vtkSegmentation* segmentation)
{
  if (!segmentation)
  {
    return;
  }

  std::map<std::string, std::string> names;
  std::vector<std::string> segmentIDs;
  segmentation->GetSegmentIDs(segmentIDs);
  for (const std::string& segmentID : segmentIDs)
  {
    vtkSegment* segment = segmentation->GetSegment(segmentID);
    if (segment)
    {
      names[segmentID] = segment->GetName() ? segment->GetName() : "";
    }
  }
  this->KnownSegmentNames[segmentation] = std::move(names);
}

bool SegmentationChangeWatcher::SegmentNameChanged(
    vtkSegmentation* segmentation,
    const std::string& segmentId)
{
  if (!segmentation || segmentId.empty())
  {
    return false;
  }

  vtkSegment* segment = segmentation->GetSegment(segmentId);
  if (!segment)
  {
    return false;
  }

  const std::string currentName = segment->GetName() ? segment->GetName() : "";
  auto& names = this->KnownSegmentNames[segmentation];
  const auto oldIt = names.find(segmentId);
  const bool changed = oldIt != names.end() && oldIt->second != currentName;
  names[segmentId] = currentName;
  return changed;
}

void SegmentationChangeWatcher::ObserveSegmentationNode(
    vtkMRMLSegmentationNode* segNode,
    int frameIndex)
{
  if (!segNode)
  {
    return;
  }

  this->NodeFrameIndex[segNode] = frameIndex;

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

  // If Slicer replaced the vtkSegmentation object (for example during history
  // restoration), detach the old object-level observers before attaching the
  // new one. This avoids stale callbacks without treating simple browser
  // navigation as an edit.
  auto oldIt = this->NodeToSeg.find(segNode);
  if (oldIt != this->NodeToSeg.end() && oldIt->second.GetPointer() != seg)
  {
    vtkSegmentation* oldSeg = oldIt->second.GetPointer();
    if (oldSeg)
    {
      auto tagIt = this->SegTags.find(oldSeg);
      if (tagIt != this->SegTags.end())
      {
        for (unsigned long tag : tagIt->second)
        {
          oldSeg->RemoveObserver(tag);
        }
        this->SegTags.erase(tagIt);
      }
      this->SegToNode.erase(oldSeg);
      this->KnownSegmentNames.erase(oldSeg);
    }
  }

  this->NodeToSeg[segNode] = seg;
  this->SegToNode[seg] = segNode;

  if (this->SegTags.find(seg) == this->SegTags.end())
  {
    std::vector<unsigned long> tags;
    tags.push_back(seg->AddObserver(
        vtkSegmentation::RepresentationModified,
        this->Callback));
    tags.push_back(seg->AddObserver(
        vtkSegmentation::SourceRepresentationModified,
        this->Callback));
    tags.push_back(seg->AddObserver(
        vtkSegmentation::SegmentModified,
        this->Callback));
    tags.push_back(seg->AddObserver(
        vtkSegmentation::SegmentAdded,
        this->Callback));
    tags.push_back(seg->AddObserver(
        vtkSegmentation::SegmentRemoved,
        this->Callback));
    this->SegTags[seg] = std::move(tags);
    this->RefreshKnownSegmentNames(seg);
  }
}

int SegmentationChangeWatcher::ResolveFrameIndex(
    vtkMRMLSegmentationNode* node) const
{
  if (!node)
  {
    return StaticFrameIndex;
  }

  const auto it = this->NodeFrameIndex.find(node);
  if (it != this->NodeFrameIndex.end())
  {
    if (it->second != UseBrowserFrameIndex)
    {
      return it->second;
    }
  }

  if (this->Mode == AcquisitionMode::Static)
  {
    return StaticFrameIndex;
  }

  return this->browser ? this->browser->GetSelectedItemNumber() : -1;
}

void SegmentationChangeWatcher::SynchronizeSegmentAdded(
    vtkMRMLSegmentationNode* sourceNode,
    const std::string& segmentId)
{
  if (!sourceNode || segmentId.empty())
  {
    return;
  }

  vtkSegmentation* sourceSegmentation = sourceNode->GetSegmentation();
  if (!sourceSegmentation || !sourceSegmentation->GetSegment(segmentId))
  {
    return;
  }

  vtkMRMLSequenceNode* sequenceSeg = this->SegmentationSequence;
  if (!sequenceSeg && this->GetSequenceSeg)
  {
    sequenceSeg = this->GetSequenceSeg();
  }
  if (!sequenceSeg)
  {
    return;
  }

  const int numberOfFrames = sequenceSeg->GetNumberOfDataNodes();
  for (int i = 0; i < numberOfFrames; ++i)
  {
    vtkMRMLSegmentationNode* frameNode =
        vtkMRMLSegmentationNode::SafeDownCast(sequenceSeg->GetNthDataNode(i));
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
            sourceSegmentation, segmentId, false))
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
  if (segmentId.empty())
  {
    return;
  }

  vtkMRMLSequenceNode* sequenceSeg = this->SegmentationSequence;
  if (!sequenceSeg && this->GetSequenceSeg)
  {
    sequenceSeg = this->GetSequenceSeg();
  }
  if (!sequenceSeg)
  {
    return;
  }

  const int numberOfFrames = sequenceSeg->GetNumberOfDataNodes();
  for (int i = 0; i < numberOfFrames; ++i)
  {
    vtkMRMLSegmentationNode* frameNode =
        vtkMRMLSegmentationNode::SafeDownCast(sequenceSeg->GetNthDataNode(i));
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

void SegmentationChangeWatcher::SynchronizeSegmentRenamed(
    vtkMRMLSegmentationNode* sourceNode,
    const std::string& segmentId)
{
  if (!sourceNode || segmentId.empty() || !sourceNode->GetSegmentation())
  {
    return;
  }

  vtkSegment* sourceSegment = sourceNode->GetSegmentation()->GetSegment(segmentId);
  if (!sourceSegment)
  {
    return;
  }
  const std::string newName = sourceSegment->GetName() ? sourceSegment->GetName() : "";

  vtkMRMLSequenceNode* sequenceSeg = this->SegmentationSequence;
  if (!sequenceSeg && this->GetSequenceSeg)
  {
    sequenceSeg = this->GetSequenceSeg();
  }
  if (!sequenceSeg)
  {
    return;
  }

  const int numberOfFrames = sequenceSeg->GetNumberOfDataNodes();
  for (int i = 0; i < numberOfFrames; ++i)
  {
    vtkMRMLSegmentationNode* frameNode =
        vtkMRMLSegmentationNode::SafeDownCast(sequenceSeg->GetNthDataNode(i));
    vtkSegmentation* frameSegmentation = frameNode ? frameNode->GetSegmentation() : nullptr;
    vtkSegment* frameSegment = frameSegmentation
        ? frameSegmentation->GetSegment(segmentId) : nullptr;
    if (frameSegment &&
        std::string(frameSegment->GetName() ? frameSegment->GetName() : "") != newName)
    {
      frameSegment->SetName(newName.c_str());
    }
  }
}

void SegmentationChangeWatcher::ApplyDeferredStructureChange(
    vtkMRMLSegmentationNode* sourceNode,
    const std::string& segmentId,
    StructureChangeType changeType)
{
  if (this->StructureSyncInProgress)
  {
    return;
  }

  this->StructureSyncInProgress = true;
  if (this->Mode == AcquisitionMode::Dynamic)
  {
    if (changeType == StructureChangeType::Added)
    {
      this->SynchronizeSegmentAdded(sourceNode, segmentId);
    }
    else if (changeType == StructureChangeType::Removed)
    {
      this->SynchronizeSegmentRemoved(segmentId);
    }
    else if (changeType == StructureChangeType::Renamed)
    {
      this->SynchronizeSegmentRenamed(sourceNode, segmentId);
    }
  }

  for (auto& kv : this->KnownSegmentNames)
  {
    this->RefreshKnownSegmentNames(kv.first);
  }
  this->StructureSyncInProgress = false;
}

void SegmentationChangeWatcher::ApplyDeferredLegacyContentChange(
    vtkMRMLSegmentationNode* segNode,
    const std::string& segmentId,
    int frameIndex)
{
  this->RefreshLegacySegmentAtFrame(segNode, segmentId, frameIndex);
}

void SegmentationChangeWatcher::RefreshLegacySegmentAtCurrentFrame(
    vtkMRMLSegmentationNode* segNode,
    const std::string& segmentId)
{
  this->RefreshLegacySegmentAtFrame(
      segNode, segmentId, this->ResolveFrameIndex(segNode));
}

void SegmentationChangeWatcher::RefreshLegacySegmentAtFrame(
    vtkMRMLSegmentationNode* segNode,
    const std::string& segmentId,
    int frameIndex)
{
  if (!segNode || segmentId.empty())
  {
    return;
  }

  auto* segmentTACs = this->GetsegmentTACs ? this->GetsegmentTACs() : nullptr;
  if (!segmentTACs || segmentTACs->empty())
  {
    return;
  }

  auto tacIt = segmentTACs->find(segmentId);
  if (tacIt == segmentTACs->end())
  {
    return;
  }

  if (frameIndex < 0 || frameIndex >= static_cast<int>(tacIt->second.size()))
  {
    return;
  }

  vtkMRMLSequenceNode* PETSequenceNode =
      this->GetSequencePET ? this->GetSequencePET() : nullptr;
  if (!PETSequenceNode || frameIndex >= PETSequenceNode->GetNumberOfDataNodes())
  {
    return;
  }

  vtkMRMLScalarVolumeNode* PETVolume = vtkMRMLScalarVolumeNode::SafeDownCast(
      PETSequenceNode->GetNthDataNode(frameIndex));
  if (!PETVolume)
  {
    return;
  }

  // Once the VTK edit callback has returned, prefer the sequence item that
  // corresponds to the frame that actually emitted the edit. This prevents a
  // queued Single-mode refresh from accidentally following later browser
  // navigation. Fall back to the proxy for legacy scenes without a sequence.
  vtkMRMLSegmentationNode* statisticsSegNode = segNode;
  vtkMRMLSequenceNode* segmentationSequence =
      this->GetSequenceSeg ? this->GetSequenceSeg() : nullptr;
  if (segmentationSequence &&
      frameIndex < segmentationSequence->GetNumberOfDataNodes())
  {
    vtkMRMLSegmentationNode* frameNode = vtkMRMLSegmentationNode::SafeDownCast(
        segmentationSequence->GetNthDataNode(frameIndex));
    if (frameNode)
    {
      statisticsSegNode = frameNode;
    }
  }

  vtkSegmentation* segmentation = statisticsSegNode->GetSegmentation();
  vtkSegment* segment = segmentation ? segmentation->GetSegment(segmentId) : nullptr;
  if (!segment)
  {
    VoxelStatistics stats;
    stats.keep = false;
    stats.empty = true;
    tacIt->second[frameIndex] = stats;
    if (this->OnSegmentTACChanged)
    {
      this->OnSegmentTACChanged(segmentId);
    }
    return;
  }

  vtkNew<vtkStringArray> segmentArray;
  segmentArray->InsertNextValue(segmentId);
  vtkSmartPointer<vtkOrientedImageData> labelmap = vtkSmartPointer<vtkOrientedImageData>::New();
  vtkSlicerSegmentationsModuleLogic::GenerateMergedLabelmapInReferenceGeometry(
      statisticsSegNode,
      PETVolume,
      segmentArray,
      vtkSegmentation::EXTENT_UNION_OF_EFFECTIVE_SEGMENTS,
      labelmap);

  VoxelStatistics stats;
  if (!labelmap || !labelmap->GetPointData() || !labelmap->GetPointData()->GetScalars())
  {
    stats.keep = false;
    stats.empty = true;
  }
  else
  {
    auto* logic = this->GetLogic ? this->GetLogic() : nullptr;
    if (!logic)
    {
      return;
    }
    stats = logic->ComputeVoxelStatistics(PETVolume, labelmap, 1);
  }

  tacIt->second[frameIndex] = stats;
  if (this->OnSegmentTACChanged)
  {
    this->OnSegmentTACChanged(segmentId);
  }

  if (this->GetSegEditCorr && this->GetSegEditCorr() && this->RunPlot)
  {
    this->RunPlot();
  }
}

void SegmentationChangeWatcher::RefreshAllLegacySegmentsAtCurrentFrame(
    vtkMRMLSegmentationNode* segNode)
{
  auto* segmentTACs = this->GetsegmentTACs ? this->GetsegmentTACs() : nullptr;
  if (!segmentTACs)
  {
    return;
  }

  std::vector<std::string> ids;
  ids.reserve(segmentTACs->size());
  for (const auto& kv : *segmentTACs)
  {
    ids.push_back(kv.first);
  }
  for (const std::string& id : ids)
  {
    this->RefreshLegacySegmentAtCurrentFrame(segNode, id);
  }
}

void SegmentationChangeWatcher::OnSegmentationChanged(
    vtkObject* caller,
    unsigned long eid,
    void* clientData,
    void* callData)
{
  auto* self = static_cast<SegmentationChangeWatcher*>(clientData);
  if (!self || self->StructureSyncInProgress)
  {
    return;
  }

  auto* segmentation = vtkSegmentation::SafeDownCast(caller);
  if (!segmentation)
  {
    if (auto* node = vtkMRMLSegmentationNode::SafeDownCast(caller))
    {
      const auto frameIt = self->NodeFrameIndex.find(node);
      const int frameIndex = frameIt != self->NodeFrameIndex.end()
          ? frameIt->second : UseBrowserFrameIndex;
      // Slicer may replace a vtkSegmentation object during state restoration.
      // Reattach object-level observers. Actual content changes are handled by
      // SegmentModified/RepresentationModified; this avoids false TAC updates
      // when a dynamic proxy changes content simply because the browser moved.
      self->ObserveSegmentationNode(node, frameIndex);
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

  if (eid == vtkSegmentation::SegmentAdded || eid == vtkSegmentation::SegmentRemoved)
  {
    if (segmentId.empty())
    {
      return;
    }

    const StructureChangeType changeType = eid == vtkSegmentation::SegmentAdded
        ? StructureChangeType::Added : StructureChangeType::Removed;

    if (self->DetectionOnly)
    {
      if (self->OnSegmentStructureChangedDetailed)
      {
        self->OnSegmentStructureChangedDetailed(
            self->ContextIndex,
            self->ResolveFrameIndex(segNode),
            segNode,
            segmentId,
            changeType);
      }
      return;
    }

    if (changeType == StructureChangeType::Added)
    {
      self->SynchronizeSegmentAdded(segNode, segmentId);
    }
    else
    {
      self->SynchronizeSegmentRemoved(segmentId);
    }
    if (self->OnSegmentStructureChanged)
    {
      self->OnSegmentStructureChanged();
    }
    self->RefreshKnownSegmentNames(segmentation);
    return;
  }

  if (eid == vtkSegmentation::SegmentModified && !segmentId.empty() &&
      self->SegmentNameChanged(segmentation, segmentId))
  {
    if (self->DetectionOnly)
    {
      if (self->OnSegmentStructureChangedDetailed)
      {
        self->OnSegmentStructureChangedDetailed(
            self->ContextIndex,
            self->ResolveFrameIndex(segNode),
            segNode,
            segmentId,
            StructureChangeType::Renamed);
      }
    }
    else if (self->OnSegmentStructureChanged)
    {
      self->OnSegmentStructureChanged();
    }
    return;
  }

  if (eid != vtkSegmentation::RepresentationModified &&
      eid != vtkSegmentation::SourceRepresentationModified &&
      eid != vtkSegmentation::SegmentModified)
  {
    return;
  }

  if (segmentId.empty())
  {
    return;
  }

  if (self->DetectionOnly)
  {
    if (self->OnSegmentContentChanged)
    {
      self->OnSegmentContentChanged(
          self->ContextIndex,
          self->ResolveFrameIndex(segNode),
          segmentId);
    }
    return;
  }

  // Legacy Single mode used to recompute TAC synchronously from inside the
  // vtkSegmentation callback. Segment Editor undo/redo restores segmentation
  // state while those callbacks are still being emitted, so re-entering
  // segmentation conversion/statistics code here can crash. If the widget
  // supplies OnSegmentContentChanged then it owns a deferred Qt refresh.
  if (self->OnSegmentContentChanged)
  {
    self->OnSegmentContentChanged(
        self->ContextIndex,
        self->ResolveFrameIndex(segNode),
        segmentId);
    return;
  }

  self->RefreshLegacySegmentAtCurrentFrame(segNode, segmentId);
}
