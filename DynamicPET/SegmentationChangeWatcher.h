#ifndef __SegmentationChangeWatcher_h
#define __SegmentationChangeWatcher_h

#include <vtkCallbackCommand.h>
#include <vtkSegmentation.h>
#include <vtkWeakPointer.h>
#include <vtkSmartPointer.h>
#include <vtkMRMLSegmentationNode.h>
#include <vtkMRMLSequenceBrowserNode.h>
#include <vtkMRMLSequenceNode.h>
#include <vtkSlicerDynamicPETLogic.h>

#include <functional>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

class SegmentationChangeWatcher : public vtkObject
{
public:
  enum class AcquisitionMode
  {
    Dynamic,
    Static
  };

  enum class StructureChangeType
  {
    Added,
    Removed,
    Renamed
  };

  static constexpr int StaticFrameIndex = -1;
  static constexpr int UseBrowserFrameIndex = -2;

  vtkMRMLSequenceBrowserNode* browser{nullptr};
  std::function<vtkMRMLSequenceNode*()> GetSequencePET;
  std::function<vtkMRMLSequenceNode*()> GetSequenceSeg;
  std::function<vtkSlicerDynamicPETLogic*()> GetLogic;
  std::function<std::map<std::string, std::vector<VoxelStatistics>>*()> GetsegmentTACs;
  std::function<bool()> GetSegEditCorr;
  std::function<void()> RunPlot;
  std::function<std::string()> GetCurrentSegID;

  // Legacy Single-mode callbacks.
  std::function<void()> OnSegmentStructureChanged;
  std::function<void(const std::string&)> OnSegmentTACChanged;

  // Acquisition-aware notification callbacks. These are used by Multi mode.
  // The watcher detects which acquisition/frame/segment changed; the widget
  // owns provenance lookup and TAC recomputation.
  std::function<void(int, int, const std::string&)> OnSegmentContentChanged;
  std::function<void(
      int,
      int,
      vtkMRMLSegmentationNode*,
      const std::string&,
      StructureChangeType)> OnSegmentStructureChangedDetailed;

  AcquisitionMode Mode{AcquisitionMode::Dynamic};
  int ContextIndex{-1};
  bool DetectionOnly{false};

  // Optional direct sequence handle for acquisition-aware dynamic watchers.
  // Legacy Single mode may continue to use GetSequenceSeg.
  vtkWeakPointer<vtkMRMLSequenceNode> SegmentationSequence;

  static SegmentationChangeWatcher* New();
  vtkTypeMacro(SegmentationChangeWatcher, vtkObject);

  SegmentationChangeWatcher();
  ~SegmentationChangeWatcher() override;

  void Clear();

  // Attach to a segmentation node. Dynamic proxies use the default
  // UseBrowserFrameIndex; sequence data nodes use an explicit >=0 frame index;
  // static segmentations use StaticFrameIndex.
  void ObserveSegmentationNode(
      vtkMRMLSegmentationNode* segNode,
      int frameIndex = UseBrowserFrameIndex);

  // Apply structural synchronization after the original Slicer event callback
  // has returned. Multi mode intentionally defers these operations to Qt.
  void ApplyDeferredStructureChange(
      vtkMRMLSegmentationNode* sourceNode,
      const std::string& segmentId,
      StructureChangeType changeType);

  // Recompute one legacy Single-mode TAC observation after the originating
  // VTK callback has returned. This is intentionally public so the Qt widget
  // can defer undo/redo-driven updates to the next event-loop turn.
  void ApplyDeferredLegacyContentChange(
      vtkMRMLSegmentationNode* segNode,
      const std::string& segmentId,
      int frameIndex);

protected:
  static void OnSegmentationChanged(vtkObject* caller,
                                    unsigned long eid,
                                    void* clientData,
                                    void* callData);

  int ResolveFrameIndex(vtkMRMLSegmentationNode* node) const;
  void RefreshKnownSegmentNames(vtkSegmentation* segmentation);
  bool SegmentNameChanged(
      vtkSegmentation* segmentation,
      const std::string& segmentId);
  void RefreshLegacySegmentAtCurrentFrame(
      vtkMRMLSegmentationNode* segNode,
      const std::string& segmentId);
  void RefreshLegacySegmentAtFrame(
      vtkMRMLSegmentationNode* segNode,
      const std::string& segmentId,
      int frameIndex);
  void RefreshAllLegacySegmentsAtCurrentFrame(vtkMRMLSegmentationNode* segNode);

  void SynchronizeSegmentAdded(
      vtkMRMLSegmentationNode* sourceNode,
      const std::string& segmentId);
  void SynchronizeSegmentRemoved(
      const std::string& segmentId);
  void SynchronizeSegmentRenamed(
      vtkMRMLSegmentationNode* sourceNode,
      const std::string& segmentId);

  std::unordered_map<vtkSegmentation*, vtkWeakPointer<vtkMRMLSegmentationNode>> SegToNode;
  std::unordered_map<vtkMRMLSegmentationNode*, vtkWeakPointer<vtkSegmentation>> NodeToSeg;
  std::unordered_map<vtkMRMLSegmentationNode*, int> NodeFrameIndex;
  std::set<vtkMRMLSegmentationNode*> ObservedNodes;
  std::map<vtkSegmentation*, std::vector<unsigned long>> SegTags;
  std::map<vtkSegmentation*, std::map<std::string, std::string>> KnownSegmentNames;
  vtkSmartPointer<vtkCallbackCommand> Callback;
  bool StructureSyncInProgress{false};
};

#endif // SegmentationChangeWatcher_h
