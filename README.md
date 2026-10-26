# AMITEX phase-field / cohesive-zone solver

This repository contains an FFT-based solver for fracture in composite
materials, combining a bulk phase-field model with a frictional cohesive-zone
model at material interfaces. The implementation is built on [AMITEX](https://amitexfftp.github.io/AMITEX/index.html) and
uses composite voxels to represent the interfaces.

For the model formulation and validation, see the paper:

> “FFT phase-field model combined with cohesive composite voxels for fracture
> of composite materials with interfaces,” *Computational Mechanics*.
> [https://doi.org/10.1007/s00466-021-02041-1](https://doi.org/10.1007/s00466-021-02041-1)

## Repository contents

- `amitex_fftp/` — a copy of AMITEX-FFTP provided for completeness. Some
  folders have been removed to reduce the repository size, so this is not a
  complete AMITEX distribution.
- `amitex_fftp/libAmitex/src/user_pfczm/` — AMITEX user-model integration and
  the `amitex_fftp` executable target for this model.
- `examples/behaviours/` — Fortran material behaviours, including the
  phase-field bulk behaviour and the Alfano–Sacco cohesive-zone behaviour.
- `examples/` — material/loading/algorithm XML files and run scripts for three
  examples (`periodic2D`, `periodic3D`, `minicomp`). The microstructures are
  not included and must be generated first (see the README below).
- `examples/microstructures/code_python/` — Python equivalents of the MATLAB
  microstructure-generation scripts and helper functions. See its README for
  dependencies, usage, and how to generate the example microstructures.

## Requirements

A working build requires a Fortran compiler, a compatible MPI implementation
and Fortran MPI wrapper, FFTW, CMake, and `make`. The compiler used for AMITEX,
its MPI wrapper, and the material-behaviour library must be compatible. The
build script uses bundled AMITEX dependencies, but you may need to configure
paths and options for your system.

## Build

### Using the bundled AMITEX-FFTP copy

Because the bundled `amitex_fftp/` tree is trimmed, running `./install` may
print errors such as `No such file or directory` or `no makefile found` for
folders that are not included in this repository. These missing-folder errors
are expected in this PFCZM distribution and do not, by themselves, prevent
building or running this model.

1. From `amitex_fftp/`, edit the configuration near the top of `install` for
   your environment. In particular, select `FC` and `MPIFC`, and configure the
   FFT implementation and FFTW include/library paths. Ensure the selected MPI
   wrapper uses the same compiler as `FC`.
2. Run the installation script:

   ```sh
   cd amitex_fftp
   ./install
   ```

3. Build the PFCZM user executable from the repository root:

   ```sh
   make -C amitex_fftp/libAmitex/src/user_pfczm lAmitex=../..
   ```

   This creates `amitex_fftp/libAmitex/src/user_pfczm/amitex_fftp`.
4. Build the material-behaviour shared library if it is missing or needs to be
   rebuilt with your compiler. From the repository root:

   ```sh
   make -C examples/behaviours \
     lDecomp=../../amitex_fftp/lib_extern/2decomp_fft-1.5.847
   ```

   The example input files expect this library at
   `examples/behaviours/libUmatAmitex.so`. If you use a different compiler,
   update `examples/behaviours/Makefile` accordingly before building.

The installation script and example Makefiles contain platform-specific
settings; inspect them before building rather than assuming their defaults
match your system.

To clean the generated AMITEX, PFCZM, and material-behaviour binaries and
object/module files, run `./clean` from the repository root.

### Using a complete AMITEX-FFTP distribution

Alternatively, download the complete AMITEX-FFTP distribution from the
[AMITEX-FFTP homepage](https://amitexfftp.github.io/AMITEX/index.html). Copy
this repository's `amitex_fftp/libAmitex/src/user_pfczm/` folder into the
distribution's `libAmitex/src/` directory, then follow the same configuration
and installation steps above. Use the corresponding AMITEX source tree when
building the PFCZM executable.

## Run examples

### Step1: prepare input for AMITEX

Three examples are provided: `periodic2D`, `periodic3D` (periodic cell with a
matrix notch) and `minicomp` (finite mini-composite with a broken fibre). Their
XML files and run scripts are in `examples/`, but **the microstructures are not
included** and must be generated first with the Python tools:

```sh
python -m pip install -r examples/microstructures/code_python/requirements.txt
cd examples/microstructures
# run the periodic2D, periodic3D and minicomp commands listed under
# "Example usage" in code_python/README.md
```

This writes, for each case, the geometry (`.vtk`) and the composite-voxel files
(`rep_3_4_2/`) to `examples/microstructures/<case>/iUC<n>/vp<..>/h<..>/`, where
the scripts and `examples/materials/mat_*.xml` expect them. To use another
microstructure, change `MICR` in the script and the single
`Coeff_composite directory` line of the material file; the algorithm and
behaviour files need no change. See
[`examples/microstructures/code_python/README.md`](examples/microstructures/code_python/README.md)
for the generation options and the preparation of the material and load files.

### Step2: Run AMITEX

The example input paths and material configuration are relative to
`examples/`, so run the solver from that directory. After completing the build,
run:

```sh
cd examples
mpirun -np 1 ../amitex_fftp/libAmitex/src/user_pfczm/amitex_fftp \
  -nm microstructures/periodic2D/iUC1/vp0.1/h0.001/iUC1_h0.001.vtk \
  -m materials/mat_periodic2D.xml \
  -c loads/load_periodic2D.xml \
  -a algorithms/algo_Miehe2.xml \
  -s results/periodic2D/test
```

Change `-np` to suit your MPI environment and available resources. The output
prefix `results/periodic2D/test` can also be adjusted as needed.

The command-line arguments select the geometry (`-nm`), material definitions
(`-m`), loading and requested outputs (`-c`), algorithm parameters (`-a`), and
the output prefix (`-s`). The geometry and material XML refer to composite-voxel
data under `microstructures/` and to the behaviour library under
`behaviours/`; keep these paths available relative to the working directory.

`examples/script_periodic2D.sh`, `script_periodic3D.sh` and `script_minicomp.sh`
show the same invocation in script form for each example (run them from
`examples/`). Their environment setup and `mpirun -np` value may need
adjustment for your local checkout and installation.

## Reference

Please cite the paper above when using this solver or its model formulation.
For AMITEX-FFTP and bundled third-party components, consult their respective
documentation and license files in `amitex_fftp/`.
