# Python microstructure tools

The workflow (`workflows.py`) generates fibre-reinforced composite
microstructures for the AMITEX phase-field/cohesive-zone solver. Fibres and
matrix are voxelised on a regular grid, and the voxels cut by a fibre surface
are written as *composite voxels* (local fibre/matrix volume fractions, interface
area, normal and tangent) so that the fibre/matrix interface and its interphase
behaviour can be modelled without resolving it with the mesh. The matrix is grown
around the fibres until a target porosity is reached, which produces the
inter-fibre pores typical of chemical vapour infiltration (CVI) composites such
as SiC/SiC. Periodic cells and finite mini-composites with free boundaries are
supported.

This directory contains Python equivalents for the MATLAB source files in the
neighboring `code_matlab/` directory. The numerical helpers are collected in
`microstructure.py`; `workflows.py` provides the shared generation and export
pipeline, and the example files `SiCSiC_*.m` are provided as PRESET configurations 
in `workflow.py`.

## Requirements

Install the scientific Python dependencies with:

```sh
python -m pip install -r code_python/requirements.txt
```

`matplotlib` is only needed when using `--plot`.

## Run a workflow

### Command-line options (`python -m code_python.workflows --help`)

| Option | Description |
| --- | --- |
| `--preset {miniComp,periodic2D,periodic3D}` | Case to generate (default `periodic3D`); see the table below. |
| `--output DIR` | Output directory (default `python_generated/`). |
| `--layout FILE` | Reuse an existing fibre layout instead of generating one. Text file with three whitespace-separated columns in millimetres, `x_center y_center radius`; lines starting with `#` are ignored. With the mini-composite option, fibres touching or crossing the fibre-region boundary are dropped. |
| Preset parameters | `--length`, `--fibre-fraction`, `--pore-fraction`, `--voxel-size`, `--radius-*`, `--min-gap`, `--max-trials`, `--z-length`, `--radius-tolerance`, `--vibration-scale`, `--mini-region`, `--mini-margin`, `--thickness-gradient`, `--notch-length`, `--notch-thickness`, `--mini-composite`, `--notch` override the fields of the chosen preset; see *Configuration parameters* below. |
| `--seed N` | Seed of the random generator, for reproducible layouts. |
| `--iuc N` | Unit-cell number used in the output names (`iUC<N>`, `layouts/iUC_<N>.txt`). Defaults to the preset's value (1, or 3 for `periodic3D`). |
| `--matrix-growth {euclidean,cross}` | How the matrix is grown around the fibres to the target porosity in the periodic presets. `euclidean` (default): periodic Euclidean-distance transform, one pass, round matrix regions. `cross`: repeated dilation with the cross-shaped element of MATLAB `imdilate(..., strel('disk',1))`, giving more diamond/rectangular regions. Ignored by `miniComp`, which deposits matrix radially. |
| `--verbose` | Print progress of the iterative generation: radius-sampling attempts, layout retries, fibre-placement candidates and the overlap reached at each relaxation iteration. |
| `--plot` | Show the central z slice of the geometry (phases 0 pore, 1 matrix, 2 fibre) and save it as `layouts/iUC_<N>_slice.png`. Requires matplotlib. |
| `--verify-plot` | Save `layouts/iUC_<N>_composite_voxel_check.png`: local fibre fraction of every voxel with the composite-voxel normals (red) and in-plane tangents (white) overlaid, plus a zoom on the largest fibre. The stored tangent is the out-of-plane fibre axis (+z, cyan dots); the in-plane tangent is z × normal. A summary line checks that normals are unit length and point out of the fibres. |

### Configuration parameters

All parameters of a case live in the `Config` dataclass in `workflows.py`; each
preset (`PRESETS`) is a ready-made `Config`. Every field can be overridden
from the command line (lengths in mm), for example
`--preset periodic2D --fibre-fraction 0.4 --voxel-size 0.0005`.

| Field / option | Meaning | Default (`preset`) |
| --- | --- | --- |
| `length` / `--length` | Side length of the square cell (for the mini-composite: of the square region in which fibres are placed) | - |
| `fibre_fraction` / `--fibre-fraction` | Target fibre area fraction | - |
| `pore_fraction` / `--pore-fraction` | Target pore fraction after matrix growth | - |
| `voxel_size` / `--voxel-size` | Voxel edge length | - |
| `radius_mean` / `--radius-mean` | Mean fibre radius (Gaussian) | 7e-3 |
| `radius_min` / `--radius-min` | Smallest fibre radius | 2e-3 |
| `radius_max` / `--radius-max` | Largest fibre radius | 12e-3 |
| `radius_std` / `--radius-std` | Standard deviation of the fibre radius | 1e-3 |
| `min_gap` / `--min-gap` | Minimum gap between neighbouring fibres | 1.1e-3 |
| `max_trials` / `--max-trials` | Maximum relaxation iterations of the periodic layout generation (unused by `miniComp`) | 200 |
| `z_length` / `--z-length` | Extrusion length in z; unset gives a single z-plane | - |
| `radius_tolerance` / `--radius-tolerance` | Relative error allowed on `fibre_fraction` when sampling the radii | 0.01 |
| `vibration_scale` / `--vibration-scale` | Size of the random moves of the periodic layout relaxation, in units of `radius_mean` (unused by `miniComp`) | 0.2 |
| `mini_composite` / `--mini-composite`, `--no-mini-composite` | Finite, non-periodic mini-composite with free boundaries and radial matrix deposition instead of a periodic cell. Fibres and matrix lie fully inside the domain, which is sized automatically (see below). Off by default; on only in the `miniComp` preset | no |
| `mini_region` / `--mini-region` | `miniComp`: only fibres whose centre lies within this fraction of half the fibre region from its centre are kept and grown | 0.85 |
| `mini_margin` / `--mini-margin` | `miniComp`: pore border (mm) kept between the outermost fibre/matrix and the domain edge (rounded up to whole voxels) | 3e-3 |
| `thickness_gradient` / `--thickness-gradient` | `miniComp`: matrix thickness of the outermost fibre relative to the centre of the fibre region; each fibre's thickness is `ratio^(r/r_max)` with `r` its distance from the centre. `1` is uniform, `<1` thins the matrix outwards, `>1` thickens it. The MATLAB script corresponds to `e` | e (2.718) |
| `notch` / `--notch {none,matrix}` | Notch in the extruded volume: `matrix` removes matrix voxels over `notch_length` × `notch_thickness` in the mid-plane at the x = 0 edge | none |
| `notch_length` / `--notch-length` | Notch extent from the x = 0 edge (mm) | 5e-3 |
| `notch_thickness` / `--notch-thickness` | Notch opening in z, centred on the mid-plane (mm) | 2e-3 |
| `default_iuc` / `--iuc` | Unit-cell number used in output names | 1 |

The fibre radii are drawn from a Gaussian truncated to `radius_min`..`radius_max`
until the area fraction is within 1% of `fibre_fraction`.

### Preset cases

| Preset | Example and main distinction |
| --- | --- |
| `miniComp` | Mini-composite with free boundaries (the mini-composite option switched on): fibres are placed in a 0.2 × 0.2 mm region, 30% fibre fraction, 15% target porosity, 1 μm voxels, extruded through 0.1 mm. No periodic wrapping; matrix is deposited radially around the fibres in the central region. The domain is then sized so that all fibres and matrix lie fully inside it with a `mini_margin` pore border (see below). No notch. |
| `periodic2D` | Periodic 2-D fibre/matrix cell: 0.1 mm, 50% fibres, 4% target porosity, 0.2 μm voxels. Matrix is grown with periodic boundaries; exported as a single z-plane. |
| `periodic3D` | Periodic 3-D cell: 0.1 mm, 50% fibres, 15% porosity, 1 μm voxels, extruded through 0.2 mm, default unit-cell number 3. A matrix-only notch is cut in the mid-plane near the x = 0 edge. |

**Mini-composite domain.** With the mini-composite option the generated domain
is not `length`: fibres are placed in a square of side `length`, the matrix is
deposited, and the domain is then set to the smallest centred square that keeps
a pore border of `mini_margin` around all fibres and matrix (it can be larger
or smaller than `length`, e.g. a larger `thickness_gradient` gives a larger
domain). The size is printed (`miniComp domain: ...`). The saved layout text
file keeps the fibre-region coordinates; the geometry is offset by half the size
difference. Use `--mini-composite` on any preset, or `--no-mini-composite` to
switch it off.


### Example usage
Run commands from `examples/microstructures/`:

Generate a 2D periodic RVE
```sh
python -m code_python.workflows --preset periodic2D --seed 213 --plot --verbose --length 0.08 --fibre-fraction 0.2 --pore-fraction 0.1 --voxel-size 1e-3 --output periodic2D --min-gap 1.3e-3
```

Generate a 3D periodic RVE with notch in matrix
```sh
python -m code_python.workflows --seed 213 --length 0.05 --fibre-fraction 0.4 --pore-fraction 0.01 --voxel-size 1e-3 --min-gap 1.3e-3 --z-length 0.07 --notch matrix --notch-thickness 1.1e-3 --notch-length 0.005 --output periodic3D --plot --verbose
```

Generate a 3D mini-composite with notch in fibre (a broken fibre)
```sh
python -m code_python.workflows --seed 213 --length 0.05 --fibre-fraction 0.4 --pore-fraction 0.01 --voxel-size 1e-3 --min-gap 1.3e-3 --z-length 0.07 --notch fibre --notch-thickness 5e-3  --mini-composite --output minicomp --plot --verbose
```

###
Each run writes a fibre distribution layout to a `.txt` file, the AMITEX composite-voxel binary files,
and a binary VTK geometry file below the chosen output directory. The default
output is `python_generated/` under the current directory. Existing source
datasets and outputs are not overwritten.


###
The `example_dev` counterpart writes the small synthetic development case:

```sh
python -m code_python.example_dev --output generated/example_dev
```

## Microstructure generation in short

1. **Layout**: fibre radii are sampled from a truncated Gaussian and the fibres
   are placed and relaxed until the target fibre fraction is reached with a
   minimum gap (periodic cell, or finite mini-composite). `--seed` makes this
   reproducible; `--layout` reloads a saved layout.
2. **Voxelisation**: fibres and a matrix grown around them (up to the target
   `pore_fraction`) are rasterised on a regular grid. Voxels cut by a fibre
   surface become *composite voxels* with local fibre/matrix fractions, interface
   area, normal and tangent.
3. **Extrusion and notch** (3-D only): the 2-D section is extruded over
   `--z-length`. Optionally a matrix notch or a broken fibre is cut at mid-z
   (`--notch matrix|fibre`); removed composite voxels are also dropped from the
   binary files (the remaining linear indices are unchanged).
4. **Output** in `<output>/iUC<n>/vp<pore_fraction>/h<voxel_size>/`:
   `iUC<n>_h<voxel_size>[_notched].vtk` (geometry), `rep_3_4_2/` (composite-voxel
   binaries) and `list_composite_materials.txt`. The folder `iUC<n>/vp<..>/h<..>/`
   is what the material file points to.

## Preparing the AMITEX input files

Working examples are `examples/materials/mat_periodic2D.xml`, `mat_periodic3D.xml`
and `examples/loads/load_periodic2D.xml`, `load_periodic3D.xml`.

**Material file (`mat_*.xml`)**. The phase numbers are fixed: 1 pore, 2 matrix,
3 fibre, 4 interface (cohesive-zone interphase, `Interphase_material numM="4"`).
The elastic and fracture parameters of each phase are edited only when the
material properties change. For every new microstructure, **only one line must
be updated**: the path to the generated folder in `Material_composite`,

```xml
<Coeff_composite directory="microstructures/periodic3D/iUC3/vp0.01/h0.001/"/>
```

It must point to the `iUC<n>/vp<..>/h<..>/` directory written by the workflow,
i.e. the one containing `rep_3_4_2/` (path relative to where AMITEX is run).
The 2-D and 3-D material files differ only in this line.

**Load file (`load_*.xml`)**. It defines the output fields and the loading steps
(`Loading Tag="1"` elastic-range step, `Tag="2"` loading to failure). Each
component is `Driving="Strain"` (imposed, with `Value` the final strain) or
`Driving="Stress"` (`Constant`, i.e. zero stress). The 2-D example loads in
x (`xx` strain), the 3-D example loads along the fibres in z (`zz` strain);
all other components are stress-free. Adapt the loaded direction, strain
values, `Nincr`/`Tfinal` and output times to the case. The load file does not
depend on the microstructure itself. 

The algorithm and behaviour (UMAT library) files need no change.

More details, see [AMITEX website](https://amitexfftp.github.io/AMITEX/user_guide/to_begin.html)

## Running the three examples (periodic2D, periodic3D, minicomp)

The repository provides the run scripts (`examples/script_*.sh`) and the XML
inputs (`examples/materials/`, `examples/loads/`, `examples/algorithms/`) of
three examples. **The microstructures are not distributed**: they are written
to `examples/microstructures/<case>/` (git-ignored) and must be generated
first, from `examples/microstructures/`, with the commands in
[Example usage](#example-usage) (`--output periodic2D`, `periodic3D` and
`minicomp`, same `--seed`). The commands are chosen to give exactly the paths
used by the scripts and material files:

| Example | Microstructure expected by the script and `mat_*.xml` |
|---|---|
| `periodic2D` | `microstructures/periodic2D/iUC1/vp0.1/h0.001/iUC1_h0.001.vtk` |
| `periodic3D` | `microstructures/periodic3D/iUC3/vp0.01/h0.001/iUC3_h0.001_notched.vtk` |
| `minicomp` | `microstructures/minicomp/iUC3/vp0.01/h0.001/iUC3_h0.001_notched.vtk` |

Then, from `examples/`, run e.g. `./script_periodic3D.sh`. If you change
`--pore-fraction`, `--voxel-size`, `--iuc` or the notch option, the folder and
file names change (`vp<..>`, `h<..>`, `iUC<n>`, `_notched` only when a notch is
cut): update `MICR` in the script and the `Coeff_composite directory` line of
the material file accordingly.
