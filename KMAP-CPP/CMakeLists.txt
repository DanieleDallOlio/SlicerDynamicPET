cmake_minimum_required(VERSION 3.15...3.29)

# -----------------------------------------------------------------------------
# KMAP version
# -----------------------------------------------------------------------------

set(KMAP_MAJOR    1)
set(KMAP_MINOR    0)
set(KMAP_REVISION 0)
set(KMAP_VERSION ${KMAP_MAJOR}.${KMAP_MINOR}.${KMAP_REVISION})

project(
  KMAP
  LANGUAGES CXX
  VERSION ${KMAP_VERSION}
  DESCRIPTION "KMAP-CPP"
)

# -----------------------------------------------------------------------------
# Detect whether KMAP is built standalone or included by another project
# -----------------------------------------------------------------------------

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(KMAP_STANDALONE ON)
else()
  set(KMAP_STANDALONE OFF)
endif()

# -----------------------------------------------------------------------------
# Compile options
# -----------------------------------------------------------------------------

option(OMP
  "Enable OpenMP support"
  OFF
)

option(VERBOSE
  "Enable verbose stdout support"
  ON
)

option(KMAP_BUILD_EXAMPLES
  "Build example applications"
  ${KMAP_STANDALONE}
)

option(KMAP_INSTALL
  "Generate installation and CMake package files"
  ${KMAP_STANDALONE}
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")

# -----------------------------------------------------------------------------
# Default build type
# -----------------------------------------------------------------------------

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(default_build_type "Release")

  message(
    STATUS
    "Setting build type to '${default_build_type}' as none was specified."
  )

  set(
    CMAKE_BUILD_TYPE
    "${default_build_type}"
    CACHE STRING
    "Choose the type of build."
    FORCE
  )

  set_property(
    CACHE CMAKE_BUILD_TYPE
    PROPERTY STRINGS
      Debug
      Release
      MinSizeRel
      RelWithDebInfo
  )
endif()

# -----------------------------------------------------------------------------
# Directories
# -----------------------------------------------------------------------------

set(SRC_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/src
)

set(INC_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

set(HPP_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/hpp
)

set(EXAMPLE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/example
)

set(CMAKE_DEBUG_POSTFIX d)

set(kmaplib kmap)
set(demo1tcm 1tcm)
set(demo2tcm 2tcm)

# -----------------------------------------------------------------------------
# Sources
# -----------------------------------------------------------------------------

file(
  GLOB SRC
  CONFIGURE_DEPENDS
  "${SRC_DIR}/*.cpp"
)

file(
  GLOB HEADER
  CONFIGURE_DEPENDS
  "${INC_DIR}/*.h"
)

file(
  GLOB HPP
  CONFIGURE_DEPENDS
  "${HPP_DIR}/*.hpp"
)

# -----------------------------------------------------------------------------
# Main KMAP library
# -----------------------------------------------------------------------------

add_library(
  ${kmaplib}
  SHARED
  ${SRC}
  ${HEADER}
  ${HPP}
)

# Allows consumers in the same build tree to use KMAP::kmap
add_library(
  KMAP::kmap
  ALIAS
  ${kmaplib}
)

set_property(
  TARGET ${kmaplib}
  PROPERTY POSITION_INDEPENDENT_CODE ON
)

target_compile_definitions(
  ${kmaplib}
  PUBLIC
    KMAP_MAJOR=${KMAP_MAJOR}
    KMAP_MINOR=${KMAP_MINOR}
    KMAP_REVISION=${KMAP_REVISION}
)

target_include_directories(
  ${kmaplib}
  PUBLIC
    $<BUILD_INTERFACE:${INC_DIR}>
    $<BUILD_INTERFACE:${HPP_DIR}>
    $<INSTALL_INTERFACE:include/kmap>
)

# -----------------------------------------------------------------------------
# Compiler settings
# -----------------------------------------------------------------------------

if(
  CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
  OR CMAKE_CXX_COMPILER_ID MATCHES "Clang"
)

  target_compile_options(
    ${kmaplib}
    PRIVATE
      -Wall
      -Wextra
      -Wno-unused-result
      -Wno-unknown-pragmas
      -Wno-narrowing
      -Wpedantic
  )

  if(CMAKE_BUILD_TYPE STREQUAL "Debug")

    target_compile_options(
      ${kmaplib}
      PRIVATE
        -Og
    )

    target_compile_definitions(
      ${kmaplib}
      PRIVATE
        DEBUG
    )

  elseif(CMAKE_BUILD_TYPE STREQUAL "Release")

    # Preserve behavior close to the previous KMAP-PY build.
    target_compile_options(
      ${kmaplib}
      PRIVATE
        -Ofast
    )

  endif()

endif()


if(MSVC)

  target_compile_options(
    ${kmaplib}
    PRIVATE
      /wd4013
      /wd4018
      /wd4028
      /wd4047
      /wd4068
      /wd4090
      /wd4101
      /wd4113
      /wd4133
      /wd4190
      /wd4244
      /wd4267
      /wd4305
      /wd4477
      /wd4996
      /wd4819
      /fp:fast
  )

  target_compile_definitions(
    ${kmaplib}
    PRIVATE
      _CRT_RAND_S
      NOMINMAX
      _USE_MATH_DEFINES
      _CRT_SECURE_NO_DEPRECATE
      _SCL_SECURE_NO_WARNINGS
  )

  set_target_properties(
    ${kmaplib}
    PROPERTIES
      WINDOWS_EXPORT_ALL_SYMBOLS ON
  )

endif()

# -----------------------------------------------------------------------------
# OpenMP
# -----------------------------------------------------------------------------

if(OMP)

  find_package(
    OpenMP
    REQUIRED
    COMPONENTS CXX
  )

  if(
    DEFINED OpenMP_CXX_VERSION
    AND OpenMP_CXX_VERSION VERSION_LESS "4.5"
  )

    message(
      FATAL_ERROR
      "OpenMP is too old. Required OpenMP >= 4.5."
    )

  endif()

  target_link_libraries(
    ${kmaplib}
    PUBLIC
      OpenMP::OpenMP_CXX
  )

  target_compile_definitions(
    ${kmaplib}
    PUBLIC
      HAVE_OPENMP
  )

  message(
    STATUS
    "OpenMP found: ${OpenMP_CXX_VERSION}"
  )

else()

  message(STATUS "OpenMP disabled")

endif()

# -----------------------------------------------------------------------------
# Verbose output
# -----------------------------------------------------------------------------

if(VERBOSE)

  target_compile_definitions(
    ${kmaplib}
    PRIVATE
      VERBOSE
  )

  message(STATUS "Verbose state: ON")

else()

  message(STATUS "Verbose state: OFF")

endif()

# -----------------------------------------------------------------------------
# Examples
# -----------------------------------------------------------------------------

if(KMAP_BUILD_EXAMPLES)

  if(EXISTS "${EXAMPLE_DIR}/demo_1tcm_realdata.cpp")

    add_executable(
      ${demo1tcm}
      ${EXAMPLE_DIR}/demo_1tcm_realdata.cpp
    )

    target_link_libraries(
      ${demo1tcm}
      PRIVATE
        KMAP::kmap
    )

  endif()


  if(EXISTS "${EXAMPLE_DIR}/demo_2tcm_realdata.cpp")

    add_executable(
      ${demo2tcm}
      ${EXAMPLE_DIR}/demo_2tcm_realdata.cpp
    )

    target_link_libraries(
      ${demo2tcm}
      PRIVATE
        KMAP::kmap
    )

  endif()

endif()

# -----------------------------------------------------------------------------
# Apple RPATH
# -----------------------------------------------------------------------------

if(APPLE)

  set_target_properties(
    ${kmaplib}
    PROPERTIES
      INSTALL_RPATH "@loader_path/../lib"
      BUILD_WITH_INSTALL_RPATH TRUE
  )

  if(TARGET ${demo1tcm})
    set_target_properties(
      ${demo1tcm}
      PROPERTIES
        INSTALL_RPATH "@loader_path/../lib"
    )
  endif()

  if(TARGET ${demo2tcm})
    set_target_properties(
      ${demo2tcm}
      PROPERTIES
        INSTALL_RPATH "@loader_path/../lib"
    )
  endif()

endif()

# =============================================================================
# Installation / export
# =============================================================================

if(KMAP_INSTALL)

  include(GNUInstallDirs)
  include(CMakePackageConfigHelpers)

  set(
    KMAP_INSTALL_INCLUDE_DIR
    "${CMAKE_INSTALL_INCLUDEDIR}/kmap"
  )

  set(
    KMAP_INSTALL_CMAKE_DIR
    "${CMAKE_INSTALL_LIBDIR}/cmake/KMAP"
  )

  set(
    KMAP_INSTALL_PKGCONFIG_DIR
    "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
  )

  # ---------------------------------------------------------------------------
  # Install library
  # ---------------------------------------------------------------------------

  install(
    TARGETS ${kmaplib}

    EXPORT KMAPTargets

    RUNTIME DESTINATION
      ${CMAKE_INSTALL_BINDIR}

    LIBRARY DESTINATION
      ${CMAKE_INSTALL_LIBDIR}

    ARCHIVE DESTINATION
      ${CMAKE_INSTALL_LIBDIR}
  )

  # ---------------------------------------------------------------------------
  # Install public headers
  # ---------------------------------------------------------------------------

  install(
    DIRECTORY
      ${INC_DIR}/

    DESTINATION
      ${KMAP_INSTALL_INCLUDE_DIR}

    FILES_MATCHING
      PATTERN "*.h"
      PATTERN "*.hpp"
  )

  install(
    DIRECTORY
      ${HPP_DIR}/

    DESTINATION
      ${KMAP_INSTALL_INCLUDE_DIR}

    FILES_MATCHING
      PATTERN "*.h"
      PATTERN "*.hpp"
  )

  # ---------------------------------------------------------------------------
  # Install examples
  # ---------------------------------------------------------------------------

  if(KMAP_BUILD_EXAMPLES)

    if(TARGET ${demo1tcm})
      install(
        TARGETS ${demo1tcm}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      )
    endif()

    if(TARGET ${demo2tcm})
      install(
        TARGETS ${demo2tcm}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      )
    endif()

  endif()

  # ---------------------------------------------------------------------------
  # pkg-config
  # ---------------------------------------------------------------------------

  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/KMAP.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/KMAP.pc
    @ONLY
  )

  install(
    FILES
      ${CMAKE_CURRENT_BINARY_DIR}/KMAP.pc

    DESTINATION
      ${KMAP_INSTALL_PKGCONFIG_DIR}
  )

  # ---------------------------------------------------------------------------
  # CMake package configuration
  # ---------------------------------------------------------------------------

  configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/KMAPConfig.cmake.in

    ${CMAKE_CURRENT_BINARY_DIR}/KMAPConfig.cmake

    INSTALL_DESTINATION
      ${KMAP_INSTALL_CMAKE_DIR}
  )

  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/KMAPConfigVersion.cmake

    VERSION
      ${PROJECT_VERSION}

    COMPATIBILITY
      SameMajorVersion
  )

  # ---------------------------------------------------------------------------
  # Export installed target
  # ---------------------------------------------------------------------------

  install(
    EXPORT KMAPTargets

    FILE
      KMAPTargets.cmake

    NAMESPACE
      KMAP::

    DESTINATION
      ${KMAP_INSTALL_CMAKE_DIR}
  )

  # ---------------------------------------------------------------------------
  # Build-tree export
  # ---------------------------------------------------------------------------

  export(
    EXPORT KMAPTargets

    FILE
      ${CMAKE_CURRENT_BINARY_DIR}/KMAPTargets.cmake

    NAMESPACE
      KMAP::
  )

  # ---------------------------------------------------------------------------
  # Install package configuration files
  # ---------------------------------------------------------------------------

  install(
    FILES
      ${CMAKE_CURRENT_BINARY_DIR}/KMAPConfig.cmake
      ${CMAKE_CURRENT_BINARY_DIR}/KMAPConfigVersion.cmake

    DESTINATION
      ${KMAP_INSTALL_CMAKE_DIR}
  )

endif()

# =============================================================================
# Summary
# =============================================================================

message(STATUS "")
message(STATUS "=================== KMAP-CPP configuration ===================")
message(STATUS "KMAP version       : ${KMAP_VERSION}")
message(STATUS "Standalone build   : ${KMAP_STANDALONE}")
message(STATUS "Build type         : ${CMAKE_BUILD_TYPE}")
message(STATUS "OpenMP             : ${OMP}")
message(STATUS "Verbose            : ${VERBOSE}")
message(STATUS "Build examples     : ${KMAP_BUILD_EXAMPLES}")
message(STATUS "Install/export     : ${KMAP_INSTALL}")
message(STATUS "===============================================================")
message(STATUS "")
