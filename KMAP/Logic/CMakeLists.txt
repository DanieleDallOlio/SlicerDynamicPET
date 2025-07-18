project(vtkSlicer${MODULE_NAME}ModuleLogic)

set(KIT ${PROJECT_NAME})

set(${KIT}_EXPORT_DIRECTIVE "VTK_SLICER_${MODULE_NAME_UPPER}_MODULE_LOGIC_EXPORT")

set(${KIT}_INCLUDE_DIRECTORIES
  PRIVATE
  ${CMAKE_SOURCE_DIR}/KMAP-PY/include
  ${CMAKE_SOURCE_DIR}/KMAP-PY/hpp
  )

set(${KIT}_SRCS
  vtkSlicer${MODULE_NAME}Logic.cxx
  vtkSlicer${MODULE_NAME}Logic.h
  )

set(${KIT}_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  )

list(APPEND ${KIT}_TARGET_LIBRARIES
  vtkSlicerSegmentationsModuleLogic
  kmap::kmap
)


if(NOT TARGET kmap::kmap)
  add_library(kmap::kmap UNKNOWN IMPORTED GLOBAL)
  set_target_properties(kmap::kmap PROPERTIES
    IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/KMAP-PY/libkmap.dylib"
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/KMAP-PY/include"
  )
endif()
#-----------------------------------------------------------------------------
SlicerMacroBuildModuleLogic(
  NAME ${KIT}
  EXPORT_DIRECTIVE ${${KIT}_EXPORT_DIRECTIVE}
  INCLUDE_DIRECTORIES ${${KIT}_INCLUDE_DIRECTORIES}
  SRCS ${${KIT}_SRCS}
  TARGET_LIBRARIES ${${KIT}_TARGET_LIBRARIES}
  )
