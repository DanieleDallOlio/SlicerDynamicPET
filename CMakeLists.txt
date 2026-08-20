cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(SlicerDynamicPET)

# -----------------------------------------------------------------------------
# CPU feature detection and compiler optimization flags
# -----------------------------------------------------------------------------
include(CheckCXXSourceCompiles)

# Check if AVX512 is supported by trying to compile a simple intrinsic
set(CMAKE_REQUIRED_FLAGS "-mavx512f")
check_cxx_source_compiles("
#include <immintrin.h>
int main() {
    __m512i a = _mm512_set1_epi32(1);
    __m512i b = _mm512_set1_epi32(2);
    __m512i c = _mm512_add_epi32(a,b);
    return _mm512_reduce_add_epi32(c);
}" HAVE_AVX512F)

if (MSVC)
  # Base optimization flags for MSVC
  set(EXTRA_OPT_FLAGS "/O2 /fp:fast /arch:AVX2")

  if (HAVE_AVX512F)
    message(STATUS "Enabling AVX512 for MSVC build")
    set(EXTRA_OPT_FLAGS "${EXTRA_OPT_FLAGS} /arch:AVX512")
  else()
    message(STATUS "Enabling AVX2 for MSVC build")
    set(EXTRA_OPT_FLAGS "${EXTRA_OPT_FLAGS} /arch:AVX2")
  endif()

  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${EXTRA_OPT_FLAGS}")

elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  # Base optimization flags for GCC/Clang
  set(EXTRA_OPT_FLAGS "-O3 -march=native -ffast-math -funroll-loops -flto")

  if (HAVE_AVX512F)
    message(STATUS "Enabling AVX512 for GCC/Clang build")
    set(EXTRA_OPT_FLAGS "${EXTRA_OPT_FLAGS} -mavx512f")
  endif()

  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${EXTRA_OPT_FLAGS}")
endif()

# -----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_CONTRIBUTORS "Daniele Dall'Olio (University of Bologna)")
set(EXTENSION_DESCRIPTION "Tools to work with dynamic PET imaging")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/DanieleDallOlio/SlicerDynamicPET/main/SlicerDynamicPET.png")
set(EXTENSION_SCREENSHOTURLS "")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

# -----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

# -----------------------------------------------------------------------------
# Optional OpenMP support
# -----------------------------------------------------------------------------

option(
  SLICERDYNAMICPET_ENABLE_OPENMP
  "Enable OpenMP acceleration when available"
  ON
)

set(SLICERDYNAMICPET_HAVE_OPENMP OFF)

if(SLICERDYNAMICPET_ENABLE_OPENMP)

  find_package(OpenMP QUIET COMPONENTS CXX)

  if(OpenMP_CXX_FOUND)

    set(SLICERDYNAMICPET_HAVE_OPENMP ON)

    message(STATUS "OpenMP found")
    message(STATUS "OpenMP version: ${OpenMP_CXX_VERSION}")
    message(STATUS "C++ OpenMP flags: ${OpenMP_CXX_FLAGS}")

    set(
      OPENMP_FOUNDINT
      TRUE
      CACHE INTERNAL "OPENMP_FOUND"
      FORCE
    )

    set(
      OPENMP_LIBRARIES
      OpenMP::OpenMP_CXX
      CACHE INTERNAL "OpenMP library target"
      FORCE
    )

    if(APPLE)
      get_target_property(
        _omp_inc
        OpenMP::OpenMP_CXX
        INTERFACE_INCLUDE_DIRECTORIES
      )

      set(
        OPENMP_INCLUDE
        ${_omp_inc}
        CACHE INTERNAL "OpenMP include directories"
        FORCE
      )
    endif()

    add_definitions(-DHAVE_OPENMP)

  else()

    message(
      STATUS
      "OpenMP not found - building SlicerDynamicPET without OpenMP acceleration"
    )

  endif()

else()

  message(
    STATUS
    "OpenMP explicitly disabled - building without OpenMP acceleration"
  )

endif()

# -----------------------------------------------------------------------------
# KMAP-CPP backend
# -----------------------------------------------------------------------------

set(
  KMAP_BUILD_EXAMPLES
  OFF
  CACHE BOOL
  "Build KMAP-CPP examples"
  FORCE
)

set(
  KMAP_INSTALL
  OFF
  CACHE BOOL
  "Install standalone KMAP-CPP package"
  FORCE
)

# KMAP-CPP should use OpenMP only if it was successfully found above.
set(
  OMP
  ${SLICERDYNAMICPET_HAVE_OPENMP}
  CACHE BOOL
  "Enable OpenMP in KMAP-CPP"
  FORCE
)

add_subdirectory(KMAP-CPP)
if(WIN32)

  set_target_properties(
    kmap
    PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY
        $<TARGET_FILE_DIR:qSlicerDynamicPETModuleWidgets>
  )

endif()


install(
  TARGETS kmap
  EXPORT KMAPTargets

  RUNTIME DESTINATION
    ${Slicer_INSTALL_QTLOADABLEMODULES_BIN_DIR}

  LIBRARY DESTINATION
    ${Slicer_INSTALL_QTLOADABLEMODULES_LIB_DIR}

  ARCHIVE DESTINATION
    ${Slicer_INSTALL_QTLOADABLEMODULES_LIB_DIR}
)

export(
  EXPORT KMAPTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/KMAPTargets.cmake"
  NAMESPACE KMAP::
)

#################################################################
#                          SUMMARY                              #
#################################################################

message(STATUS ""                                                                    )
message(STATUS "=================== SlicerDynamicPET configuration Summary =================="   )
message(STATUS "   SlicerDynamicPET "       )
message(STATUS ""                                                                    )
message(STATUS "   Build type : ${CMAKE_BUILD_TYPE}"                                 )
message(STATUS "   C++ :"                                                            )
message(STATUS "      C++ Compiler : ${CMAKE_CXX_COMPILER}"                          )
message(STATUS "      C++ flags    :"                                                )
foreach(FLAG ${CMAKE_CXX_FLAGS_LIST})
  message(STATUS "                    * ${FLAG}"                                     )
endforeach(FLAG)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "      C++ Debug flags:"                                            )
  foreach(FLAG ${CMAKE_CXX_FLAGS_DEBUG})
    message(STATUS "                    * ${FLAG}"                                   )
  endforeach(FLAG)
elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
  message(STATUS "      C++ Release flags  :"                                        )
  foreach(FLAG ${CMAKE_CXX_FLAGS_RELEASE})
    message(STATUS "                    * ${FLAG}"                                   )
  endforeach(FLAG)
endif()
message(STATUS ""                                                                    )
message(STATUS "   OpenMP support : ${OPENMP_FOUND}"                                 )
message(STATUS ""                                                                    )

# -----------------------------------------------------------------------------
# Extension modules
add_subdirectory(DynamicPET)
## NEXT_MODULE
# -----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
