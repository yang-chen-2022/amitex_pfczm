"""Unidirectional (UD) fibre composite RVE with composite interface voxels.

Geometry (all parameters are defined in the PARAMETERS block below):
  * periodic along x; fibres are never cut by the y = 0 / y = y_length borders;
  * a void layer at the top (y = y_length) of the RVE;
  * a notch along x in the mid-plane z = z_length / 2, starting at y = 0;
  * voxels cut by a fibre surface are exported as AMITEX composite voxels;
  * the AMITEX inputs (mat.xml, load.xml, algo.xml, script.sh) are generated too.

Everything is configured in the PARAMETERS blocks below. Run with ``python ud_notched_rve.py`` (or ``python -m code_python.ud_notched_rve``
from the ``microstructures`` directory). Phase IDs in the VTK: 0 void/notch,
1 matrix, 2 fibre.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from code_python.microstructure import (  # noqa: E402
    CVparam_2,
    GenDataWithGaussianDistrib,
    generationFibres,
    getCompVox_circl_periodic,
    getOrientVec_periodic,
    saveBin,
    saveVolvtk_amitex,
    write_composite_voxel_data,
)
from code_python.workflows import _adjust_fractions  # noqa: E402

# ----------------------------------------------------------------------------
# PARAMETERS (lengths in mm)
# ----------------------------------------------------------------------------
X_LENGTH = 0.1            # RVE size along x (periodic)
Y_LENGTH = 0.12           # RVE size along y, including the void layer
Z_LENGTH = 0.1           # RVE size along z (fibre direction)
VOXEL_SIZE = 0.5e-3        # cubic voxel edge

FIBRE_FRACTION = 0.5      # target fibre volume fraction of the fibre-bearing region
RADIUS_MEAN = 9e-3        # fibre radius: truncated Gaussian (mean, std, min, max)
RADIUS_STD = 1e-3
RADIUS_MIN = 2e-3
RADIUS_MAX = 12e-3
MIN_GAP = 2*VOXEL_SIZE          # minimum surface-to-surface distance between fibres
EDGE_MARGIN = 2*VOXEL_SIZE     # minimum distance between fibres and the y borders / void layer
RADIUS_TOLERANCE = 0.01   # relative error allowed on the sampled fibre fraction

VOID_THICKNESS = VOXEL_SIZE*10    # void layer thickness at y = Y_LENGTH (0 disables it)
INTERPHASE_FRACTION = 1e-3  # fibre/matrix interphase fraction in composite voxels

NOTCH = True
NOTCH_LENGTH = 0.1*Y_LENGTH      # notch extent from the y = 0 edge
NOTCH_THICKNESS = 3*VOXEL_SIZE    # notch opening in z, centred on z = Z_LENGTH / 2
NOTCH_CUTS_FIBRES = True  # False: only matrix is removed; True: fibres too

MAX_TRIALS = 100          # relaxation iterations per attempt
MAX_ATTEMPTS = 100         # new radius samples if overlaps remain
VIBRATION_SCALE = 0.2     # random-move size, in units of RADIUS_MEAN
SEED = 0                  # None for a different layout at each run

# Directories. SIM_DIR is where AMITEX is launched from: paths written in the
# input files and in script.sh are relative to it (the `examples` folder).
SIM_DIR = Path(__file__).resolve().parents[2]
OUTPUT_DIR = SIM_DIR / "microstructures" / "ud_notched_rve"  # VTK + composite voxels
CASE = "ud_notched"       # names the input files: materials/mat_<CASE>.xml, etc.

# ----------------------------------------------------------------------------
# AMITEX SIMULATION PARAMETERS
# ----------------------------------------------------------------------------
def _compute_gc(lc: float, E: float, sigmac: float) -> float:
    """Compute the phase-field fracture toughness parameter g_c."""
    return 256 * sigmac**2 * lc / (27 * E)

# Elastic constants (E [MPa], nu) and phase-field data (lc [mm], gc [N/mm], sigmac [MPa]).
PORE = dict(E=1, nu=0.2, viscosity=0.0, lc=0.002, gc=0.01)
MATRIX = dict(E=3e3, nu=0.2, viscosity=0.0, lc=0.002, sigmac=75, gc=0.01)
FIBRE = dict(E=80e3, nu=0.2, viscosity=0.0, lc=0.002, sigmac=1.82e3, gc=0.01)
MATRIX["gc"] = _compute_gc(MATRIX["lc"], MATRIX["E"], MATRIX["sigmac"])
FIBRE["gc"] = _compute_gc(FIBRE["lc"], FIBRE["E"], FIBRE["sigmac"])

# Cohesive-zone law of the fibre/matrix interface (Alfano-Sacco).
# Each entry is either
#   * a scalar (same value everywhere), or
#   * a function f(d) -> array, where d is the array of distances [mm] between
#     each interface voxel centre and the top edge of the RVE (y = Y_LENGTH), or
#   * an array with one value per interface voxel.
# Per-voxel entries are written to binary files and referenced with
# Type="Constant_Zone"; zone2.bin then numbers the interface voxels 1..n.
# Write your own function, or use one of the three ready-made profiles:
def two_step(d_switch: float, near: float, far: float):
    """`near` for d < d_switch, `far` beyond."""
    return lambda d: np.where(d < d_switch, near, far)

def linear(d0: float, v0: float, d1: float, v1: float):
    """Linear from v0 at d = d0 to v1 at d = d1, constant outside [d0, d1]."""
    return lambda d: v0 + (v1 - v0) * np.clip((d - d0) / (d1 - d0), 0.0, 1.0)

def exponential(v_top: float, v_far: float, decay_length: float):
    """v_top at the top edge, tending to v_far with e-folding length `decay_length`."""
    return lambda d: v_far + (v_top - v_far) * np.exp(-d / decay_length)


# Examples (replace the scalars below):
#   sigma_n=two_step(0.03, near=50.0, far=100.0)
#   sigma_n=linear(0.0, 50.0, Y_LENGTH, 100.0)
#   sigma_n=exponential(50.0, 100.0, 0.02)
#   GcI=lambda d: 0.002 * (1 + d / Y_LENGTH)
INTERFACE = dict(
    sigma_n=75.0, GcI=1.6e-5,       # normal strength [MPa], fracture energy [N/mm]
    sigma_t=153.0, GcII=1.6e-5,      # tangential strength, fracture energy
    eta=0.999, mu=0.0,
)
# Reference medium: set both, or leave both None to use 0.5 * (min + max) of the
# Lame coefficients (lambda, mu separately) over the pore, matrix and fibre.
REFERENCE_E = None        # [MPa]
REFERENCE_NU = None
BEHAVIOUR_LIB = "behaviours/libUmatAmitex.so"  # relative to SIM_DIR

# Loading sequence, applied in order. Components: xx yy zz xy xz yz.
# ("Stress", None) -> constant stress; ("Strain", v) -> strain increasing linearly to v.
_FREE = {"xx": ("Stress", None), "yy": ("Stress", None), "xy": ("Stress", None),
         "xz": ("Stress", None), "yz": ("Stress", None)}
LOADINGS = [
    dict(nincr=10, tfinal=100, outputs=[10],
         load={**_FREE, "zz": ("Strain", 0.01)}),
    dict(nincr=500, tfinal=2100,
         outputs=[50, 100, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 300, 350, 400, 450, 500],
         load={**_FREE, "zz": ("Strain", 0.11)}),
]
TEMPERATURE = 20.0
VTK_STRESS, VTK_STRAIN = False, False
# Internal variables written to the VTK files, by material number
# (1 pore, 2 matrix, 3 fibre, 4 interface, 5 composite voxels).
VTK_INTVARS = {1: [], 2: [1], 3: [1], 4: [], 5: [28]}

# Algorithm
ALGO = dict(
    conv_criterion=1e-4, conv_acceleration=True, nitermax=3000,
    forced=True, nit_cvfor=1000, ncvfor=5,
    laminate_nmax_subdivision=8,
    laminate_accept_unconverged=True,
    laminate_jacobian="numerical",   # mfront_tangent, elastic, mfront_secant, numerical
    laminate_perturbation=1e-8,
    laminate_init="linear",          # proportional, linear, default
    laminate_conv_criterion=1e-2,
    filter="Default", small_perturbations=True,
    damage_tol=1e-6, damage_nitermax=3000, damage_max_composite=1e-3,
    damage_acceleration=True, laplacian="modifiedLap",  # defaultLap, modifiedLap
)

# Run
N_PROCS = 16
AMITEX_EXE = "../amitex_fftp/libAmitex/src/user_pfczm/amitex_fftp"   # relative to SIM_DIR
ENV_SCRIPT = "../amitex_fftp/env_amitex.sh"                        # relative to SIM_DIR
RESULTS_DIR = "results/ud_notched"                                  # relative to SIM_DIR
# ----------------------------------------------------------------------------


def generate_layout(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fibre centres/radii: periodic in x, fully inside the fibre band in y."""
    band = Y_LENGTH - VOID_THICKNESS  # matrix + fibres, below the void layer
    count = int(round(FIBRE_FRACTION * X_LENGTH * band / (np.pi * RADIUS_MEAN**2)))
    sample = {"nn": count, "AV": RADIUS_MEAN, "MIN": RADIUS_MIN, "MAX": RADIUS_MAX, "SIG": RADIUS_STD}
    param = {"eps": MIN_GAP, "trialsMAX": MAX_TRIALS, "ur_sc": VIBRATION_SCALE * RADIUS_MEAN}
    for _ in range(MAX_ATTEMPTS):
        for _ in range(20000):
            radii = GenDataWithGaussianDistrib(sample, rng)
            vf = np.sum(np.pi * radii**2) / (X_LENGTH * band)
            if abs(vf - FIBRE_FRACTION) / FIBRE_FRACTION <= RADIUS_TOLERANCE:
                break
        else:
            raise RuntimeError("Could not sample radii matching FIBRE_FRACTION")
        # Centres live in a y-range shrunk by the largest radius, so no fibre
        # can reach a y border (the shrunk range is treated as periodic during relaxation).
        pad = float(radii.max()) + EDGE_MARGIN
        if band - 2 * pad <= 0:
            raise ValueError("Y_LENGTH - VOID_THICKNESS is too small for the fibres")
        x, y, overlap = generationFibres(
            radii, param, [[0, X_LENGTH], [pad, band - pad]], rng
        )
        if overlap <= 0:
            return x, y, radii
    raise RuntimeError("Overlap-free layout not found: lower FIBRE_FRACTION or raise MAX_TRIALS")


def discretize(x_c, y_c, radii):
    """Voxelise: returns volume (ny, nx), composite-voxel data in 2D (F-order index)."""
    h = VOXEL_SIZE
    nx, ny = round(X_LENGTH / h), round(Y_LENGTH / h)
    cols, rows = np.meshgrid(np.arange(nx), np.arange(ny))
    xp, yp = (cols + 0.5) * h, (rows + 0.5) * h
    uc = {"xlim0": 0.0, "xlim1": X_LENGTH, "ylim0": 0.0, "ylim1": Y_LENGTH, "vxsiz": h}
    volume = np.ones((ny, nx), dtype=np.uint8)
    positions, fractions, areas, normals = [], [], [], []
    for cx, cy, r in zip(x_c, y_c, radii):
        dx = (xp - cx + X_LENGTH / 2) % X_LENGTH - X_LENGTH / 2
        volume[dx**2 + (yp - cy) ** 2 <= r**2] = 2
        # Fibres never touch the y borders, so the periodic distance is exact here.
        row, col = np.nonzero(getCompVox_circl_periodic(cx, cy, xp, yp, r, uc))
        if not len(row):
            continue
        normal, dist = getOrientVec_periodic(cx, cy, xp[row, col], yp[row, col], r, uc)
        local = np.column_stack((normal * ((r - dist) / h)[:, None], np.zeros(len(row))))
        local_n = np.column_stack((normal, np.zeros(len(row))))
        surface, frac = CVparam_2(local, local_n, 0)
        ghost = surface == 0
        frac[ghost] = [1, 0]
        for a, b in ((0, 1), (1, 0)):
            nan = np.isnan(frac[:, a])
            frac[nan, a] = 1 - frac[nan, b]
        positions.append(row + col * ny + 1)
        fractions.append(frac)
        areas.append(surface)
        normals.append(local_n)
    pos = np.concatenate(positions)
    # A voxel cut by two fibres would appear twice; AMITEX needs unique positions
    pos, keep = np.unique(pos, return_index=True)
    return volume, pos, keep, np.vstack(fractions), np.concatenate(areas), np.vstack(normals), nx, ny


def _lame(E: float, nu: float) -> tuple[float, float]:
    return E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))


def _rel(path: Path) -> str:
    return Path(os.path.relpath(path, SIM_DIR)).as_posix()


def _phase_field_material(num: int, name: str, p: dict) -> str:
    lam, mu = _lame(p["E"], p["nu"])
    return f"""<!-- MATERIAL {num}: {name} -->
<Material numM="{num}" Lib="{BEHAVIOUR_LIB}" Law="elasisoPF" >
    <Coeff Index="1" Type="Constant" Value="{lam:.10g}"/> <!--lambda [MPa]-->
    <Coeff Index="2" Type="Constant" Value="{mu:.10g}"/>  <!--mu [MPa]-->
    <Coeff Index="3" Type="Constant" Value="{p['viscosity']:.10g}"/>  <!--viscous parameter [MPa*s]-->
    <Coeff Index="4" Type="Constant" Value="{p['lc']:.10g}"/>   <!--characteristic length, lc [mm]-->
    <Coeff Index="5" Type="Constant" Value="{p['gc']:.10g}"/> <!--fracture energy, gc [N/mm]-->
    <IntVar Index="1" Type="Constant" Value="0."/> <!--damage variable-->
    <IntVar Index="2" Type="Constant" Value="0."/> <!--strain energy-->
    <Coeff_composite Index="1" Type="Constant" Value="1" /> <!--calculated in amitex-->
    <Coeff_composite Index="2" Type="Constant" Value="1" /> <!--calculated in amitex-->
</Material>
"""


def resolve_interface(distance: np.ndarray) -> dict:
    """Scalars stay floats; functions of the top-edge distance and arrays become float arrays."""
    n = len(distance)
    resolved = {}
    for key, value in INTERFACE.items():
        if callable(value):
            value = value(distance)
        if np.ndim(value) == 0:
            resolved[key] = float(value)
            continue
        value = np.asarray(value, dtype=float).reshape(-1)
        if value.size != n:
            raise ValueError(
                f"INTERFACE['{key}'] has {value.size} values, expected {n} interface voxels"
            )
        resolved[key] = value
    return resolved


def write_mat_xml(path: Path, interface_values: dict, n_zones: int) -> None:
    if (REFERENCE_E is None) != (REFERENCE_NU is None):
        raise ValueError("REFERENCE_E and REFERENCE_NU must both be set or both be None")
    if REFERENCE_E is None:
        lames = np.array([_lame(p["E"], p["nu"]) for p in (PORE, MATRIX, FIBRE)])
        lam0, mu0 = 0.5 * (lames.min(axis=0) + lames.max(axis=0))
    else:
        lam0, mu0 = _lame(REFERENCE_E, REFERENCE_NU)
    def coeff(index: int, key: str, comment: str) -> str:
        value = interface_values[key]
        if np.ndim(value) == 0:
            attrs = f'Type="Constant" Value="{value:.10g}"'
        else:
            file = OUTPUT_DIR / "interface_zones" / f"{key}.bin"
            saveBin(value, file, "double")
            attrs = f'Type="Constant_Zone" File="{_rel(file)}" Format="binary"'
        return f'    <Coeff Index="{index}" {attrs}/>  <!--{comment}-->\n'

    interface = f"""<!-- MATERIAL 4: INTERFACE -->
<Material numM="4" Lib="{BEHAVIOUR_LIB}" Law="CZM_AlfanoSacco_explicit" >
""" + (
        coeff(1, "sigma_n", "normal strength [MPa]")
        + coeff(2, "GcI", "GcI [N/mm]")
        + coeff(3, "sigma_t", "tangential strength [MPa]")
        + coeff(4, "GcII", "GcII [N/mm]")
        + coeff(5, "eta", "eta []")
        + coeff(6, "mu", "mu []")
    ) + """    <Coeff Index="7" Type="Constant" Value="0"/>  <!--interphase thickness [mm], recalculated in amitex-->
""" + "".join(
        f'    <Coeff Index="{k}" Type="Constant" Value="0"/>  <!--e{(k - 8) // 3 + 1}_{(k - 8) % 3 + 1}, allocated from the rep files-->\n'
        for k in range(8, 14)
    ) + """    <Coeff Index="14" Type="Constant" Value="0"/>  <!--viscous parameter for regularising interphase damage evolution-->
    <IntVar Index="1" Type="Constant" Value="0."/> <!--damage variable-->
    <IntVar Index="2" Type="Constant" Value="0."/> <!--residual inelastic shear displ.-->
    <IntVar Index="3" Type="Constant" Value="0."/> <!--dir. of shear at each step-->
    <IntVar Index="4" Type="Constant" Value="0."/> <!--dir. of shear at each step-->
    <IntVar Index="5" Type="Constant" Value="0."/> <!--damage variable at t-1-->
    <IntVar Index="6" Type="Constant" Value="0."/> <!--damage variable of the bulk-->
    <Coeff_composite Index="1" Type="Constant" Value="1" /> <!--calculated in amitex-->
    <Coeff_composite Index="2" Type="Constant" Value="1" /> <!--calculated in amitex-->
</Material>
"""
    text = f"""<?xml version="1.0" encoding="UTF-8"?>
<Materials>

<!-- REFERENCE MATERIAL -->
<Reference_Material Lambda0="{lam0:.10g}" Mu0="{mu0:.10g}"/>

{_phase_field_material(1, "PORE", PORE)}
{_phase_field_material(2, "MATRIX", MATRIX)}
{_phase_field_material(3, "FIBRE", FIBRE)}
{interface}
<!-- COMPOSITE MATERIALS -->
<Material_composite>
    <Coeff_composite directory="{_rel(OUTPUT_DIR)}/"/>
</Material_composite>

<!-- Specify which mat is the INTERPHASE -->
<Interphase>
     <Interphase_material  numM="4" Nzones="{n_zones}"/>
</Interphase>

</Materials>
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def write_load_xml(path: Path) -> None:
    flag = lambda b: "1" if b else "0"
    intvars = "\n".join(
        f'    <vtk_IntVarList numM="{m}"> {" ".join(map(str, v))} </vtk_IntVarList>'
        for m, v in sorted(VTK_INTVARS.items())
    )
    loadings = []
    for tag, step in enumerate(LOADINGS, 1):
        comps = []
        for c in ("xx", "yy", "zz", "xy", "xz", "yz"):
            driving, value = step["load"][c]
            if value is None:
                comps.append(f'    <{c} Driving="{driving}" Evolution="Constant"/>')
            else:
                comps.append(f'    <{c} Driving="{driving}" Evolution="Linear" Value="{value:.10g}"/>')
        loadings.append(f"""<Loading Tag="{tag}">
    <Time_Discretization Discretization="Linear" Nincr="{step['nincr']}" Tfinal="{step['tfinal']:.10g}" />
    <Output_vtkList> {" ".join(f"{t:.10g}" for t in step['outputs'])} </Output_vtkList>
{chr(10).join(comps)}
    <T Evolution="Constant"/>
    <Param Index="1" Evolution="Constant"/>
</Loading>
""")
    text = f"""<?xml version="1.0" encoding="UTF-8"?>
<Loading_Output>

<Output>
    <vtk_StressStrain Strain="{flag(VTK_STRAIN)}" Stress="{flag(VTK_STRESS)}"/>
{intvars}
</Output>

<InitLoadExt>
    <T Value="{TEMPERATURE:.10g}"/>
    <Param Index="1" Value="0."/>
</InitLoadExt>

{chr(10).join(loadings)}
</Loading_Output>
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def write_algo_xml(path: Path) -> None:
    a = ALGO
    b = lambda v: "true" if v else "false"
    text = f"""<?xml version="1.0" encoding="UTF-8"?>
<Algorithm_Parameters>
    <Algorithm Type="basic_scheme">
        <Convergence_Criterion Value="{a['conv_criterion']:g}"/>
        <Convergence_Acceleration Value="{b(a['conv_acceleration'])}"/>
        <Nitermax Value="{a['nitermax']}"/>
        <Convergence_Forced Value="{b(a['forced'])}" Nit_cvfor="{a['nit_cvfor']}" Ncvfor="{a['ncvfor']}"/>
    </Algorithm>

    <Algorithm_laminate>
        <Nmax_subdivision Value="{a['laminate_nmax_subdivision']}"/>
        <Accept_unconverged Value="{b(a['laminate_accept_unconverged'])}"/>
        <Jacobian_matrix Value="{a['laminate_jacobian']}"/>
        <Numerical_Jacobian_perturbation Value="{a['laminate_perturbation']:g}"/>
        <Initialisation_type Value="{a['laminate_init']}"/>
        <Convergence_Acceleration Value="true"/>
        <Convergence_Criterion Value="{a['laminate_conv_criterion']:g}"/>
    </Algorithm_laminate>

    <Mechanics>
        <Filter Type="{a['filter']}"/>
        <Small_Perturbations Value="{b(a['small_perturbations'])}"/>
    </Mechanics>

    <User Algo="standard">
         <P_real Index="1" Value="{a['damage_tol']:g}"/>  <!-- convergence criterion, damage field -->
         <P_real Index="2" Value="{a['damage_nitermax']}"/>  <!-- max iterations, damage field -->
         <P_real Index="3" Value="{a['damage_max_composite']:g}"/>  <!-- max damage in composite voxels -->
         <P_string Index="1" Value="{b(a['damage_acceleration'])}"/>  <!-- convergence acceleration, damage field -->
         <P_string Index="2" Value="{a['laplacian']}"/>  <!-- defaultLap, modifiedLap -->
    </User>

</Algorithm_Parameters>
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def write_script(path: Path, vtk: Path, mat: Path, load: Path, algo: Path) -> None:
    text = f"""#!/bin/bash
# Run from the simulation directory so that the relative paths below resolve.
cd "$(dirname "$0")"

source {ENV_SCRIPT}

AMITEX="{AMITEX_EXE}"
MPIRUNalias='mpirun -np {N_PROCS}'

MICR="{_rel(vtk)}"
MATE="{_rel(mat)}"
LOAD="{_rel(load)}"
ALGO="{_rel(algo)}"
RESU="{RESULTS_DIR}"

mkdir -p $RESU
$MPIRUNalias $AMITEX -nm $MICR -m $MATE -c $LOAD -a $ALGO -s $RESU"/test"
"""
    path.write_text(text)
    path.chmod(0o755)


def main() -> None:
    rng = np.random.default_rng(SEED)
    x_c, y_c, radii = generate_layout(rng)
    h = VOXEL_SIZE
    band = Y_LENGTH - VOID_THICKNESS
    assert np.all(y_c - radii >= EDGE_MARGIN - 1e-12) and np.all(y_c + radii <= band - EDGE_MARGIN + 1e-12)
    print(f"Fibres: {len(radii)}; fibre fraction in fibre band: "
          f"{np.sum(np.pi * radii**2) / (X_LENGTH * band):.4f}")

    volume2d, pos, keep, frac, area, normal, nx, ny = discretize(x_c, y_c, radii)
    print(pos)
    frac = _adjust_fractions(frac, INTERPHASE_FRACTION)[keep]
    area, normal = area[keep], normal[keep]

    # Void layer on top (y = Y_LENGTH): matrix voxels whose centre is above the band.
    y_centres = (np.arange(ny) + 0.5) * h
    volume2d[y_centres > band, :] = 0

    nz = round(Z_LENGTH / h)
    plane = nx * ny
    volume = np.repeat(volume2d[:, :, None], nz, axis=2)
    positions = (pos[None, :] + plane * np.arange(nz)[:, None]).ravel()
    fractions = np.tile(frac, (nz, 1))
    areas = np.tile(area, nz)
    normals = np.tile(normal, (nz, 1))

    if NOTCH:
        z = np.arange(nz)[None, None, :] + 0.5
        in_slab = np.abs(z * h - Z_LENGTH / 2) < NOTCH_THICKNESS / 2
        y = (np.arange(ny)[:, None, None] + 0.5) * h
        region = (y < NOTCH_LENGTH) & in_slab
        composite = np.zeros(volume.shape, dtype=bool)
        composite.ravel(order="F")[positions - 1] = True
        mask = region & (composite | (volume == 1) | (NOTCH_CUTS_FIBRES & (volume == 2)))
        volume[mask] = 0
        keep_cv = ~mask.ravel(order="F")[positions - 1]
        positions, fractions, areas, normals = (
            positions[keep_cv], fractions[keep_cv], areas[keep_cv], normals[keep_cv]
        )
        print(f"Notch: {np.count_nonzero(mask)} voxels removed")

    tangents = np.tile([0.0, 0.0, 1.0], (len(positions), 1))
    out = Path(OUTPUT_DIR)
    voxel_dir = out / "rep_3_4_2"  # phase IDs (fibre, interphase, matrix)
    voxel_dir.mkdir(parents=True, exist_ok=True)
    (out / "list_composite_materials.txt").write_text("3 4 2 laminate\n")
    n_cv = len(positions)
    # Distance of each interface voxel centre to the top edge (y = Y_LENGTH);
    # positions are 1-based Fortran-order indices over (y, x, z).
    rows = (positions.astype(np.int64) - 1) % ny
    distance = Y_LENGTH - (rows + 0.5) * h
    interface_values = resolve_interface(distance)
    per_voxel = any(np.ndim(v) > 0 for v in interface_values.values())
    # One zone per interface voxel only when some property varies voxel by voxel.
    zone2 = np.arange(1, n_cv + 1, dtype=np.uint64) if per_voxel else None
    write_composite_voxel_data(voxel_dir, positions, fractions, areas, normals, tangents, zone2)
    vtk_path = out / f"ud_notched_rve_h{h:g}.vtk"
    saveVolvtk_amitex(volume, vtk_path, "uint8", "cell", [h] * 3)
    np.savetxt(out / "fibre_layout.txt", np.column_stack((x_c, y_c, radii)), fmt="%.10g",
               header="x_center_mm y_center_mm radius_mm")
    mat, load, algo = (SIM_DIR / d / f"{n}_{CASE}.xml" for d, n in
                       (("materials", "mat"), ("loads", "load"), ("algorithms", "algo")))
    write_mat_xml(mat, interface_values, n_cv if per_voxel else 1)
    write_load_xml(load)
    write_algo_xml(algo)
    write_script(SIM_DIR / f"script_{CASE}.sh", vtk_path, mat, load, algo)
    print(f"AMITEX inputs written to {SIM_DIR} (script_{CASE}.sh)")
    print(f"Grid {nx} x {ny} x {nz}; composite voxels: {len(positions)}; output: {out}")


if __name__ == "__main__":
    main()
