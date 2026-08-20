project(vtkSlicer${MODULE_NAME}ModuleLogic)

set(KIT ${PROJECT_NAME})

set(${KIT}_EXPORT_DIRECTIVE "VTK_SLICER_${MODULE_NAME_UPPER}_MODULE_LOGIC_EXPORT")


if (${OPENMP_FOUNDINT})
  if(MSVC)
    set(OPENMP_INCLUDE "")
    set(OPENMP_LIBRARIES "")
  endif()
  set(${KIT}_INCLUDE_DIRECTORIES
    PRIVATE
    ${EIGEN3_INCLUDE_DIRS}
    ${OPENMP_INCLUDE}
    )
  set(${KIT}_TARGET_LIBRARIES
    ${ITK_LIBRARIES}
    ${EIGEN3_LIBRARIES}
    ${OPENMP_LIBRARIES}
    vtkSlicerSegmentationsModuleLogic
    vtkSlicerSequencesModuleMRML
    kmap
    )
    add_definitions("-DHAVE_OPENMP")
else()
  set(${KIT}_INCLUDE_DIRECTORIES
    PRIVATE
    ${EIGEN3_INCLUDE_DIRS}
    )
  set(${KIT}_TARGET_LIBRARIES
    ${ITK_LIBRARIES}
    ${EIGEN3_LIBRARIES}
    vtkSlicerSegmentationsModuleLogic
    vtkSlicerSequencesModuleMRML
    kmap
    )
endif()

set(${KIT}_SRCS
  vtkSlicer${MODULE_NAME}Logic.cxx
  vtkSlicer${MODULE_NAME}Logic.h
  )


#-----------------------------------------------------------------------------
SlicerMacroBuildModuleLogic(
  NAME ${KIT}
  EXPORT_DIRECTIVE ${${KIT}_EXPORT_DIRECTIVE}
  INCLUDE_DIRECTORIES ${${KIT}_INCLUDE_DIRECTORIES}
  SRCS ${${KIT}_SRCS}
  TARGET_LIBRARIES ${${KIT}_TARGET_LIBRARIES}
  )
