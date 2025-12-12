cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(DynamicPETExtension)

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
set(EXTENSION_HOMEPAGE "https://www.slicer.org/wiki/Documentation/Nightly/Extensions/DynamicPETExtension")
set(EXTENSION_CONTRIBUTORS "Daniele Dall'Olio (University of Bologna)")
set(EXTENSION_DESCRIPTION "Tools to work with dynamic PET imaging")
set(EXTENSION_ICONURL "https://www.example.com/Slicer/Extensions/DynamicPETExtension.png")
set(EXTENSION_SCREENSHOTURLS "https://www.example.com/Slicer/Extensions/DynamicPETExtension/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

# -----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

find_package(OpenMP REQUIRED COMPONENTS CXX)
set(OPENMP_FOUNDINT ${OPENMP_FOUND} CACHE INTERNAL "OPENMP_FOUND")
if (OPENMP_FOUND)
  message(STATUS "OpenMP found")
  if (OpenMP_CXX_VERSION_MAJOR LESS_EQUAL 2 AND OpenMP_CXX_VERSION_MINOR LESS 0)
    message(FATAL_ERROR "Your OpenMP is too old. Required OpenMP 2.0. Please upgrade.")
  endif()
  message(STATUS "C++ openmp flags: ${OpenMP_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(OPENMP_LIBRARIES OpenMP::OpenMP_CXX CACHE INTERNAL "OPENMP library target")
  if (APPLE)
    get_target_property(_omp_inc OpenMP::OpenMP_CXX INTERFACE_INCLUDE_DIRECTORIES)
    set(OPENMP_INCLUDE ${_omp_inc} CACHE INTERNAL "OpenMP include directories")
  endif()
  add_definitions("-DHAVE_OPENMP")
endif()

# -----------------------------------------------------------------------------
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


#################################################################
#                          SUMMARY                              #
#################################################################

message(STATUS ""                                                                    )
message(STATUS "=================== DynamicPETExtension configuration Summary =================="   )
message(STATUS "   DynamicPETExtension "       )
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
add_subdirectory(KMAP)
## NEXT_MODULE
# -----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
