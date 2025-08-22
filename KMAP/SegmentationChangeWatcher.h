#pragma once

#include <vtkCallbackCommand.h>
#include <vtkSegmentation.h>
#include <vtkWeakPointer.h>
#include <vtkSmartPointer.h>
#include <vtkMRMLSegmentationNode.h>
#include <vtkMRMLSequenceBrowserNode.h>
#include <vtkSlicerKMAPLogic.h>

#include <unordered_map>

class SegmentationChangeWatcher : public vtkObject
{
public:
  vtkMRMLSequenceBrowserNode* browser;
  std::function<vtkMRMLSequenceNode*()> GetSequencePET;
  std::function<vtkSlicerKMAPLogic*()> GetLogic;
  std::function<std::map<std::string, std::vector<VoxelStatistics>>*()> GetsegmentTACs;
  std::function<bool()> GetSegEditCorr;
  std::function<void()> RunPlot;
  static SegmentationChangeWatcher* New();
  vtkTypeMacro(SegmentationChangeWatcher, vtkObject);

  // Attach to one segmentation node
  void ObserveSegmentationNode(vtkMRMLSegmentationNode* segNode);

protected:
  SegmentationChangeWatcher();

  // Called when a segment geometry is modified
  static void OnSegmentationChanged(vtkObject* caller,
                                    unsigned long eid,
                                    void* clientData,
                                    void* callData);

  // Map vtkSegmentation* → owning MRML node
  std::unordered_map<vtkSegmentation*, vtkWeakPointer<vtkMRMLSegmentationNode>> SegToNode;
  std::set<vtkMRMLSegmentationNode*> ObservedNodes;
  std::map<vtkSegmentation*, unsigned long> SegTags;
  vtkSmartPointer<vtkCallbackCommand> Callback;
};
