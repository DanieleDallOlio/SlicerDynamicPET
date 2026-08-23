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
  vtkMRMLSequenceBrowserNode* browser{nullptr};
  std::function<vtkMRMLSequenceNode*()> GetSequencePET;
  std::function<vtkMRMLSequenceNode*()> GetSequenceSeg;
  std::function<vtkSlicerDynamicPETLogic*()> GetLogic;
  std::function<std::map<std::string, std::vector<VoxelStatistics>>*()> GetsegmentTACs;
  std::function<bool()> GetSegEditCorr;
  std::function<void()> RunPlot;
  std::function<std::string()> GetCurrentSegID;

  // Called after a segment is added/removed and its structure has been
  // synchronized across the segmentation sequence.
  std::function<void()> OnSegmentStructureChanged;

  // Called after an already-computed segment TAC has been refreshed for the
  // currently displayed dynamic frame.
  std::function<void(const std::string&)> OnSegmentTACChanged;

  static SegmentationChangeWatcher* New();
  vtkTypeMacro(SegmentationChangeWatcher, vtkObject);

  SegmentationChangeWatcher();
  ~SegmentationChangeWatcher() override;

  void Clear();

  // Attach to one segmentation proxy node.
  void ObserveSegmentationNode(vtkMRMLSegmentationNode* segNode);

protected:
  static void OnSegmentationChanged(vtkObject* caller,
                                    unsigned long eid,
                                    void* clientData,
                                    void* callData);

  void SynchronizeSegmentAdded(
      vtkMRMLSegmentationNode* proxyNode,
      const std::string& segmentId);
  void SynchronizeSegmentRemoved(
      const std::string& segmentId);

  std::unordered_map<vtkSegmentation*, vtkWeakPointer<vtkMRMLSegmentationNode>> SegToNode;
  std::set<vtkMRMLSegmentationNode*> ObservedNodes;
  std::map<vtkSegmentation*, std::vector<unsigned long>> SegTags;
  vtkSmartPointer<vtkCallbackCommand> Callback;
};

#endif // SegmentationChangeWatcher_h
