# SlicerDynamicPET

<p align="center">
  <img src="SlicerDynamicPET.png" alt="SlicerDynamicPET" width="220">
</p>

**SlicerDynamicPET** is a 3D Slicer extension for importing, visualizing, and
quantitatively analyzing **dynamic positron emission tomography (PET)** data.

The extension combines a dedicated dynamic PET DICOM importer with tools for
time-activity curve extraction, kinetic modeling, and voxel-wise parametric
imaging in a single Slicer workflow.

---

## Included modules

### DynamicPET

**DynamicPET** is the main analysis module of SlicerDynamicPET.

It provides tools for:

- working with dynamic PET Volume Sequences;
- defining and extracting regional time-activity curves;
- graphical kinetic analysis;
- compartment-model fitting;
- voxel-wise kinetic modeling;
- generation and visualization of parametric images.

DynamicPET uses
[KMAP-CPP](https://github.com/DanieleDallOlio/KMAP-CPP)
as its C++ kinetic-modeling backend.

KMAP-CPP is derived from the open-source
[KMAP-C Toolkit](https://github.com/ShareKM/KMAP-C).

### dPETImporter

**dPETImporter** is the DICOM import component of the extension.

It recognizes dynamic PET acquisitions and reconstructs them as time-resolved
**Slicer Volume Sequences**.

The importer supports:

- 3D-per-frame dynamic PET datasets;
- 2D slice stacks grouped into temporal frames;
- extraction of frame timing information;
- optional SUVbw conversion when the required DICOM metadata are available;
- integration with Slicer's DICOM and Sequences infrastructure.

The importer is maintained independently at:

[dPETImporter](https://github.com/DanieleDallOlio/dPETImporter)

---

## Main features

- Dynamic PET DICOM import
- Reconstruction of time-resolved Volume Sequences
- Frame timing extraction and preservation
- Optional SUVbw conversion during import
- Regional time-activity curve extraction
- Graphical kinetic modeling
- Compartment-model fitting
- Voxel-wise parametric imaging
- Integration with Slicer segmentation, subject hierarchy, plotting, and visualization
- Optional OpenMP acceleration when available

---

## Typical workflow

1. Import a dynamic PET study using the **DICOM** module.
2. Load the acquisition using **dPETImporter** as a Volume Sequence.
3. Define or load anatomical or functional regions of interest.
4. Extract regional time-activity curves.
5. Select the desired kinetic-modeling approach.
6. Perform ROI-level graphical or compartmental modeling.
7. Generate voxel-wise parametric images when required.
8. Review the resulting curves, fitted parameters, and parametric volumes in Slicer.

---

## Installation

### Slicer Extensions Manager

SlicerDynamicPET is intended to be distributed through the
**3D Slicer Extensions Manager**.

After installation, restart Slicer if requested.

### Development build

SlicerDynamicPET contains a C++ loadable module and therefore must be built
against a compatible 3D Slicer source build.

The extension source tree contains:

```text
SlicerDynamicPET/
├── DynamicPET/
├── KMAP-CPP/
└── dPETImporter/
```

`KMAP-CPP` and `dPETImporter` are maintained as independent repositories and
are integrated into SlicerDynamicPET for packaging and distribution.

---

## OpenMP acceleration

SlicerDynamicPET can use **OpenMP** to accelerate computationally intensive
operations when OpenMP support is available.

The extension is intended to remain functional without OpenMP by falling back
to serial implementations.

---

## Screenshots

<p align="center">
  <img src="Screenshots/DynamicPET-main.png" alt="DynamicPET module" width="850">
</p>

Additional screenshots may include:

- dynamic PET DICOM import using dPETImporter;
- time-activity curve visualization and kinetic fitting;
- voxel-wise parametric imaging results.

---

## Architecture

SlicerDynamicPET separates user-interface and data-management functionality
from the kinetic-modeling backend:

```text
SlicerDynamicPET
├── dPETImporter
│   └── DICOM import and Volume Sequence creation
│
├── DynamicPET
│   └── Slicer user interface, TAC handling, kinetic analysis, and parametric imaging
│
└── KMAP-CPP
    └── standalone C++ kinetic-modeling backend
```

This separation allows KMAP-CPP and dPETImporter to remain independently
maintainable while SlicerDynamicPET provides the integrated end-user workflow.

---

## KMAP-CPP and upstream KMAP-C

The kinetic-modeling backend used by DynamicPET is
[KMAP-CPP](https://github.com/DanieleDallOlio/KMAP-CPP).

KMAP-CPP is a maintained C++ derivative of the
[Open KMAP-C Toolkit](https://github.com/ShareKM/KMAP-C), originally developed
within the Open Kinetic Modeling Initiative.

The original KMAP-C authorship, copyright notices, and MIT licensing are
preserved in KMAP-CPP.

---

## Research use

SlicerDynamicPET is research software under active development.

It is not a medical device and is not intended to replace validated clinical
software. Results should be independently verified for the datasets, tracers,
models, and workflows in which the extension is used.

---

## Citation

A publication specifically describing SlicerDynamicPET will be added here when
available.

When using the extension in scientific work, please also cite 3D Slicer and,
where appropriate, the original KMAP work and related methodological
publications used by the selected kinetic models.

---

## License

Original code developed for SlicerDynamicPET is distributed under the
**MIT License**.

Third-party and derived components retain their respective copyright and
license notices. In particular:

- **KMAP-CPP** retains the upstream KMAP-C MIT license and attribution;
- **3D Slicer** is distributed under its own BSD-style software license.

See [LICENSE](LICENSE) and the license files of bundled components for details.

---

## Author

**Daniele Dall'Olio**  
University of Bologna

- [GitHub](https://github.com/DanieleDallOlio)
- [University profile](https://www.unibo.it/sitoweb/daniele.dallolio)

<a href="https://github.com/UniboDIFABiophysics">
<div class="image">
<img src="https://cdn.rawgit.com/physycom/templates/697b327d/logo_unibo.png" width="45" height="45">
</div>
</a>

---

## Acknowledgments

SlicerDynamicPET is built on the 3D Slicer platform and its MRML, DICOM,
Sequences, segmentation, plotting, and visualization infrastructure.

The kinetic-modeling backend builds upon the Open KMAP-C Toolkit. Original
KMAP-C contributors and the Open Kinetic Modeling Initiative are acknowledged
for the underlying kinetic-modeling framework.
