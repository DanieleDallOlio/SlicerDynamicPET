# KMAP-CPP

## Overview

**KMAP-CPP** is a C++ tracer kinetic modeling library derived from the
[Open KMAP-C Toolkit](https://github.com/ShareKM/KMAP-C), developed as part
of the Open Kinetic Modeling Initiative.

The original **Kinetic Modeling and Analysis Package (KMAP)** is an open-source
software environment for implementing and applying tracer kinetic models to
dynamic positron emission tomography (PET) data. The project was initially
developed at the University of California, Davis, and subsequently released as
open-source software to support the
[Open Kinetic Modeling Initiative](https://www.openkmi.org/).

KMAP-CPP builds upon the original KMAP-C implementation while maintaining a
standalone C++ library that can also serve as the kinetic-modeling backend of
[SlicerDynamicPET](https://github.com/DanieleDallOlio/SlicerDynamicPET).

## Relationship to KMAP-C

KMAP-CPP is a fork and derivative of the
[ShareKM/KMAP-C](https://github.com/ShareKM/KMAP-C) project.

KMAP-C provides the original C/C++ implementations of core functionality for
dynamic PET kinetic modeling, including:

- input-function processing;
- tracer kinetic models;
- optimization algorithms;
- utility functions;
- voxel-wise kinetic modeling and OpenMP-based acceleration.

KMAP-CPP retains and extends this foundation with modifications developed for
standalone C++ use and integration with SlicerDynamicPET.

The upstream KMAP-C project remains the source and reference for the original
KMAP implementation. KMAP-CPP is independently maintained as a derivative of
that project.

## KMAP-CPP

The goal of KMAP-CPP is to provide a reusable C++ kinetic-modeling backend
without requiring MATLAB, Python, or 3D Slicer.

The library can be:

- built and used as a standalone C++ library;
- linked from other CMake projects;
- used for voxel-wise kinetic modeling with OpenMP acceleration;
- used as the computational kinetic-modeling backend of SlicerDynamicPET.

The current implementation includes extensions and modifications to the
upstream KMAP-C codebase, including changes to the build system, numerical and
optimization routines, utilities, and functionality required for efficient
voxel-wise parametric imaging.

## Building

KMAP-CPP requires a C++17-compatible compiler and CMake.

A standard Release build can be configured with:

```bash
cmake \
  -S . \
  -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DOMP=OFF
```

and built with:

```bash
cmake --build build --parallel
```

### OpenMP

OpenMP acceleration can be enabled with:

```bash
cmake \
  -S . \
  -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DOMP=ON
```

followed by:

```bash
cmake --build build --parallel
```

## Installation

KMAP-CPP provides CMake package configuration and pkg-config metadata.

For example:

```bash
cmake \
  -S . \
  -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DKMAP_INSTALL=ON \
  -DCMAKE_INSTALL_PREFIX=/path/to/install

cmake --build build --parallel
cmake --install build
```

An installed KMAP-CPP library can then be consumed from another CMake project
using:

```cmake
find_package(KMAP REQUIRED CONFIG)

target_link_libraries(
  my_application
  PRIVATE
    KMAP::kmap
)
```

## Examples and data

Example applications for the one-tissue and two-tissue compartment models are
provided in the `example/` directory.

Example/reference data used by the standalone implementations are available in
the `data/` directory.

Examples can be disabled when KMAP-CPP is embedded into another project:

```bash
-DKMAP_BUILD_EXAMPLES=OFF
```

## MATLAB/MEX sources

The repository retains MATLAB/MEX-related source material inherited from or
derived from the KMAP ecosystem in the `mex/` directory.

The core KMAP-CPP library itself is intended to remain independent of MATLAB.

## SlicerDynamicPET

KMAP-CPP is used as the kinetic-modeling backend of
[SlicerDynamicPET](https://github.com/DanieleDallOlio/SlicerDynamicPET), a
3D Slicer extension for dynamic PET analysis, ROI kinetic modeling, and
voxel-wise parametric imaging.

SlicerDynamicPET-specific functionality, including its graphical user
interface, MRML integration, image handling, and visualization, is maintained
separately from this library.

## Upstream project and attribution

KMAP-CPP is based on the open-source
[KMAP-C Toolkit](https://github.com/ShareKM/KMAP-C).

The original KMAP project was developed at the University of California,
Davis, and released to support the
[Open Kinetic Modeling Initiative](https://www.openkmi.org/).

Original KMAP-C authorship and contributions are retained and acknowledged.
See [CONTRIBUTORS](CONTRIBUTORS.md) for contributor information.

Users of KMAP-CPP are encouraged to also consult and acknowledge the original
KMAP-C project where appropriate.

## License

KMAP-CPP is distributed under the MIT License, consistent with the upstream
KMAP-C project.

The original copyright and license notices from KMAP-C are retained.

See [LICENSE](LICENSE) for details.
