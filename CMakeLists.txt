cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(DynamicPETExtension)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://www.slicer.org/wiki/Documentation/Nightly/Extensions/DynamicPETExtension")
set(EXTENSION_CONTRIBUTORS "Daniele Dall'Olio (University of Bologna)")
set(EXTENSION_DESCRIPTION "Tools to work with dynamic PET imaging")
set(EXTENSION_ICONURL "https://www.example.com/Slicer/Extensions/DynamicPETExtension.png")
set(EXTENSION_SCREENSHOTURLS "https://www.example.com/Slicer/Extensions/DynamicPETExtension/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

find_package(OpenMP REQUIRED COMPONENTS CXX)
set(OPENMP_FOUNDINT ${OPENMP_FOUND} CACHE INTERNAL "OPENMP_FOUND")
if (OPENMP_FOUND)
  message(STATUS "OpenMP found")
  if (OpenMP_CXX_VERSION_MAJOR LESS_EQUAL 2 AND OpenMP_CXX_VERSION_MINOR LESS 0)
    message(FATAL_ERROR " Your OpenMP is too old. Required OpenMP 2.0. Please upgrade.")
  endif()
  message(STATUS "OpenMP found")
  message(STATUS "C++ openmp flags: ${OpenMP_CXX_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  # if (APPLE)
  #   list(APPEND linked_libs OpenMP::OpenMP_CXX)
  # endif()
  if (APPLE)
    set(OPENMP_LIBRARIES OpenMP::OpenMP_CXX CACHE INTERNAL "OPENMP library target")
    get_target_property(_omp_inc OpenMP::OpenMP_CXX INTERFACE_INCLUDE_DIRECTORIES)
    set(OPENMP_INCLUDE ${_omp_inc} CACHE INTERNAL "OpenMP include directories")
  endif()
  add_definitions("-DHAVE_OPENMP")
endif()

#-----------------------------------------------------------------------------
# External library
add_subdirectory(KMAP-PY)
if(NOT APPLE)
  set_target_properties(kmap PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY $<TARGET_FILE_DIR:qSlicerKMAPModuleWidgets>
  )
  install(TARGETS kmap
    RUNTIME DESTINATION ${Slicer_INSTALL_QTLOADABLEMODULES_BIN_DIR}
    LIBRARY DESTINATION ${Slicer_INSTALL_QTLOADABLEMODULES_LIB_DIR}
    ARCHIVE DESTINATION ${Slicer_INSTALL_QTLOADABLEMODULES_LIB_DIR}
  )
endif()

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(KMAP)
## NEXT_MODULE
#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
