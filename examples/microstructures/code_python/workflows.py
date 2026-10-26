"""Configured microstructure-generation workflows from the MATLAB examples."""

from __future__ import annotations

import argparse
from dataclasses import dataclass, fields, replace
from pathlib import Path

import numpy as np
from scipy.ndimage import binary_dilation, distance_transform_edt

from .microstructure import (
    CVparam_2,
    GenDataWithGaussianDistrib,
    generationFibres,
    getCompVox_circl_periodic,
    getOrientVec_periodic,
    saveVolvtk_amitex,
    write_composite_voxel_data,
)


@dataclass(frozen=True)
class Config:
    """All parameters of a microstructure case; ``PRESETS`` holds ready-made configurations."""

    length: float
    fibre_fraction: float
    pore_fraction: float
    voxel_size: float
    radius_mean: float = 7e-3
    radius_min: float = 2e-3
    radius_max: float = 12e-3
    radius_std: float = 1e-3
    min_gap: float = 1.1e-3
    max_trials: int = 200
    z_length: float | None = None
    mini_composite: bool = False
    notch: str | None = None
    default_iuc: int = 1
    radius_tolerance: float = 0.01  # relative error allowed on fibre_fraction
    vibration_scale: float = 0.2  # random-move size, in units of radius_mean
    mini_region: float = 0.85  # miniComp: radius of the kept region / half cell
    thickness_gradient: float = float(np.e)  # miniComp: outermost/central matrix thickness ratio
    mini_margin: float = 3e-3  # miniComp: pore border kept around fibres and matrix
    notch_length: float = 2e-2  # notch extent from the x = 0 edge
    notch_thickness: float = 4e-3  # notch opening in z (centred in the mid-plane)
    notch_fibre: int | None = None  # "fibre" notch: 0-based index of the fibre to break (random if unset)


NOTCH_TYPES = ("matrix", "fibre")


# Named, ready-made configurations selected with --preset.
PRESETS = {
    "miniComp": Config(
        0.2, 0.3, 0.15, 1e-3, min_gap=1.1e-3, z_length=0.1,
        mini_composite=True,
    ),
    "periodic2D": Config(0.1, 0.5, 0.04, 0.2e-3, min_gap=0.7e-3),
    "periodic3D": Config(
        0.1, 0.5, 0.15, 1e-3, min_gap=1e-3, z_length=0.2,
        default_iuc=3,
    ),
}


def _radius_sample(
    config: Config, rng: np.random.Generator, verbose: bool = False
) -> tuple[np.ndarray, float]:
    count = int(round(config.fibre_fraction * config.length**2 / (np.pi * config.radius_mean**2)))
    params = {
        "nn": count,
        "AV": config.radius_mean,
        "MIN": config.radius_min,
        "MAX": config.radius_max,
        "SIG": config.radius_std,
    }
    for attempt in range(1, 20001):
        radii = GenDataWithGaussianDistrib(params, rng)
        vf = float(np.sum(np.pi * radii**2) / config.length**2)
        if abs(vf - config.fibre_fraction) / config.fibre_fraction <= config.radius_tolerance:
            if verbose:
                print(
                    f"  Radius sampling accepted on attempt {attempt}: "
                    f"fibre fraction={vf:.6f}"
                )
            return radii, vf
        if verbose and attempt % 1000 == 0:
            print(f"  Radius sampling attempt {attempt}/20000; latest fraction={vf:.6f}")
    raise RuntimeError("Could not sample fibre radii within 1% of the target fibre fraction")


def _load_layout(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not path.is_file():
        raise FileNotFoundError(f"Fibre-layout text file does not exist: {path}")
    try:
        data = np.loadtxt(path, comments="#", ndmin=2)
    except ValueError as error:
        raise ValueError(
            f"Could not parse {path}; expected columns: x_center_mm y_center_mm radius_mm"
        ) from error
    if data.ndim != 2 or data.shape[1] != 3 or data.shape[0] == 0:
        raise ValueError(
            f"{path} must contain three numeric columns: x_center_mm y_center_mm radius_mm"
        )
    if not np.all(np.isfinite(data)) or np.any(data[:, 2] <= 0):
        raise ValueError(f"{path} contains non-finite values or non-positive radii")
    return data[:, 0], data[:, 1], data[:, 2]


def _generate_nonperiodic_layout(
    config: Config, rng: np.random.Generator, verbose: bool = False
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Place whole, non-overlapping fibres inside a finite cell without wrapping."""
    length = config.length
    for layout_attempt in range(1, 101):
        if verbose:
            print(f"  Finite placement attempt {layout_attempt}/100")
        radii, _ = _radius_sample(config, rng, verbose)
        accepted: list[tuple[float, float]] = []
        accepted_radii: list[float] = []
        for index in np.argsort(radii)[::-1]:
            radius = float(radii[index])
            if 2 * radius > length:
                break
            for _attempt in range(10000):
                x = float(rng.uniform(radius, length - radius))
                y = float(rng.uniform(radius, length - radius))
                if all(
                    np.hypot(x - old_x, y - old_y)
                    >= radius + old_radius + config.min_gap
                    for (old_x, old_y), old_radius in zip(accepted, accepted_radii)
                ):
                    accepted.append((x, y))
                    accepted_radii.append(radius)
                    if verbose:
                        print(
                            f"    Placed fibre {len(accepted)}/{len(radii)} "
                            f"after {_attempt + 1} candidate(s)"
                        )
                    break
            else:
                if verbose:
                    print(
                        f"    Could not place fibre {len(accepted) + 1}/{len(radii)} "
                        "within 10000 candidates; retrying layout"
                    )
                break
        if len(accepted) == len(radii):
            centers = np.asarray(accepted)
            placed_radii = np.asarray(accepted_radii)
            vf = float(np.sum(np.pi * placed_radii**2) / length**2)
            if verbose:
                print(f"  Finite placement succeeded; fibre fraction={vf:.6f}")
            return centers[:, 0], centers[:, 1], placed_radii, vf
    raise RuntimeError(
        "Could not place every miniComp fibre fully inside the non-periodic domain. "
        "Try a different --seed or lower the target fibre fraction."
    )


def _generate_layout(
    config: Config, rng: np.random.Generator, verbose: bool = False
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, int, float]:
    if config.mini_composite:
        x, y, radii, actual_vf = _generate_nonperiodic_layout(config, rng, verbose)
        return x, y, radii, actual_vf, 1, 0.0
    param = {
        "eps": config.min_gap,
        "trialsMAX": config.max_trials,
        "ur_sc": config.vibration_scale * config.radius_mean,
    }
    for retry in range(1, 21):
        if verbose:
            print(f"Layout generation retry {retry}/20")
        radii, actual_vf = _radius_sample(config, rng, verbose)
        x, y, max_overlap = generationFibres(
            radii, param, [[0, config.length], [0, config.length]], rng, verbose
        )
        if max_overlap <= 0 or retry == 20:
            if verbose:
                state = "without overlaps" if max_overlap <= 0 else "with unresolved overlaps"
                print(f"Layout generation finished {state}; max overlap={max_overlap:.6g}")
            return x, y, radii, actual_vf, retry, max_overlap
    raise RuntimeError("Unable to generate fibre layout")


def _periodic_dilate(mask: np.ndarray) -> np.ndarray:
    tiled = np.tile(mask, (3, 3))
    structure = np.asarray([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=bool)
    return binary_dilation(tiled, structure=structure)[
        mask.shape[0] : 2 * mask.shape[0], mask.shape[1] : 2 * mask.shape[1]
    ]


MATRIX_GROWTH_METHODS = ("euclidean", "cross")


def _grow_matrix(solid: np.ndarray, target_pore: float, method: str) -> np.ndarray:
    """Grow the fibre mask periodically until the pore fraction reaches the target.

    ``euclidean`` grows by a round (Euclidean-distance) radius in one pass;
    ``cross`` repeats the 4-connected dilation of the MATLAB ``strel('disk',1)``.
    """
    n_rows, n_cols = solid.shape
    if method == "euclidean":
        distance = distance_transform_edt(~np.tile(solid, (3, 3)))[
            n_rows : 2 * n_rows, n_cols : 2 * n_cols
        ]
        n_solid = int(np.ceil((1.0 - target_pore) * solid.size))
        if n_solid <= np.count_nonzero(solid):
            return solid
        threshold = np.partition(distance.ravel(), n_solid - 1)[n_solid - 1]
        return distance <= threshold
    if method != "cross":
        raise ValueError(f"Unknown matrix growth {method!r}; choose from {MATRIX_GROWTH_METHODS}")
    pore = 1.0 - np.count_nonzero(solid) / solid.size
    for _ in range(2 * max(solid.shape)):
        if pore - target_pore <= 1e-3:
            break
        solid = _periodic_dilate(solid)
        pore = 1.0 - np.count_nonzero(solid) / solid.size
    return solid


def _extrusion_depth(config: Config) -> int:
    if config.z_length is None:
        return 1
    depth = int(np.ceil(config.z_length / config.voxel_size))
    return depth + (depth % 2 == 0)


def _adjust_fractions(fractions: np.ndarray, interface_fraction: float) -> np.ndarray:
    result = np.asarray(fractions, dtype=float).copy()
    if len(result) == 0:
        return result
    for row in result:
        if np.any(row < 1e-4):
            largest, smallest = np.argmax(row), np.argmin(row)
            row[largest], row[smallest] = 1 - 1e-4, 1e-4
    rows = np.arange(len(result))
    result[rows, np.argmax(result, axis=1)] -= interface_fraction
    low = result[:, 0] < 1e-4
    delta = 1e-4 - result[low, 0]
    result[low, 1] -= delta
    result[low, 0] = 1e-4
    low = result[:, 1] < 1e-4
    delta = 1e-4 - result[low, 1]
    result[low, 0] -= delta
    result[low, 1] = 1e-4
    return result


def _grid_size(config: Config) -> int:
    n = int(np.ceil(round(config.length / config.voxel_size, 6)))
    return n + (n % 2 == 0)


def _thickness_factors(distances: np.ndarray, ratio: float) -> np.ndarray:
    """Relative matrix thickness of each fibre: ``ratio ** (distance / max distance)``.

    ``ratio`` is the thickness at the outermost fibre relative to the centre, so
    1 is uniform, values below 1 thin the matrix outwards and values above 1
    thicken it (the MATLAB script corresponds to ``ratio = e``).
    """
    if ratio <= 0:
        raise ValueError("thickness_gradient must be positive")
    max_distance = np.max(distances, initial=0.0)
    if max_distance <= 0:
        return np.ones_like(distances)
    return np.power(ratio, distances / max_distance)


def _mini_matrix(
    x_phys: np.ndarray,
    y_phys: np.ndarray,
    x_centers: np.ndarray,
    y_centers: np.ndarray,
    radii: np.ndarray,
    config: Config,
    region_length: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Radial matrix deposition of the miniComp; returns the solid mask and skipped fibres."""
    h = config.voxel_size
    center = config.length / 2
    mini_radius = region_length * 0.5 * config.mini_region
    distances = np.hypot(x_centers - center, y_centers - center)
    skipped = distances > mini_radius
    max_distance = np.max(distances, initial=0.0)
    gradient = _thickness_factors(distances, config.thickness_gradient)
    solid = np.zeros_like(x_phys, dtype=bool)
    last = int(np.ceil(region_length / h)) + 1
    for counter in range(1, last + 1):
        for i in np.flatnonzero(~skipped):
            dx, dy = x_phys - x_centers[i], y_phys - y_centers[i]
            solid |= dx**2 + dy**2 <= (radii[i] + gradient[i] * h * counter) ** 2
        pore = 1.0 - np.count_nonzero(solid) / (np.pi * (mini_radius / h) ** 2)
        if pore <= config.pore_fraction:
            break
    return solid, skipped


def _mini_padding(
    x_centers: np.ndarray, y_centers: np.ndarray, radii: np.ndarray, config: Config
) -> int:
    """Voxels to add (or remove, if negative) on every side of the fibre region.

    The result is the smallest symmetric domain in which fibres and matrix keep a
    pore border of ``mini_margin`` to the edge of the domain.
    """
    h = config.voxel_size
    margin = max(1, int(np.ceil(config.mini_margin / h)))
    pad = int(np.ceil(0.25 * config.length / h))
    for _ in range(12):
        trial = replace(config, length=config.length + 2 * pad * h)
        n = _grid_size(trial)
        grid = np.arange(1, n + 1) * h
        y_phys, x_phys = np.meshgrid(grid, grid, indexing="ij")
        solid, _ = _mini_matrix(
            x_phys, y_phys, x_centers + pad * h, y_centers + pad * h, radii,
            trial, config.length,
        )
        rows, cols = np.nonzero(solid)
        if len(rows) == 0:
            return 0
        free = min(rows.min(), cols.min(), n - 1 - rows.max(), n - 1 - cols.max())
        if free >= margin:
            return pad - (free - margin)
        pad *= 2
    raise RuntimeError("Could not fit the miniComp matrix inside a finite domain")


def _discretize(
    x_centers: np.ndarray,
    y_centers: np.ndarray,
    radii: np.ndarray,
    config: Config,
    matrix_growth: str = "euclidean",
    region_length: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Voxelize a layout; ``region_length`` is the side of the miniComp fibre region."""
    h, length = config.voxel_size, config.length
    region_length = region_length or length
    n = _grid_size(config)
    y_grid, x_grid = np.mgrid[1 : n + 1, 1 : n + 1]
    x_phys, y_phys = x_grid * h, y_grid * h
    uc = {"xlim0": 0.0, "xlim1": length, "ylim0": 0.0, "ylim1": length, "vxsiz": h}
    material = np.zeros((n, n), dtype=np.uint8) if config.mini_composite else np.full((n, n), 2, dtype=np.uint8)
    positions: list[np.ndarray] = []
    fractions: list[np.ndarray] = []
    areas: list[np.ndarray] = []
    normals: list[np.ndarray] = []
    tangents: list[np.ndarray] = []

    center = np.asarray([length / 2, length / 2])
    mini_radius = region_length * 0.5 * config.mini_region
    skipped = np.zeros(len(radii), dtype=bool)
    for i, (cx, cy, radius) in enumerate(zip(x_centers, y_centers, radii)):
        if config.mini_composite and np.linalg.norm(np.asarray([cx, cy]) - center) > mini_radius:
            skipped[i] = True
            continue
        if config.mini_composite:
            dx, dy = x_phys - cx, y_phys - cy
        else:
            dx = (x_phys - cx + length / 2) % length - length / 2
            dy = (y_phys - cy + length / 2) % length - length / 2
        inside = dx**2 + dy**2 <= radius**2
        material[inside] = 2 if config.mini_composite else 1
        interface = getCompVox_circl_periodic(
            cx, cy, x_phys, y_phys, radius, uc,
            periodic=not config.mini_composite,
        )
        row, col = np.nonzero(interface)
        if not len(row):
            continue
        position = row + col * n + 1
        normal, radial_distance = getOrientVec_periodic(
            cx, cy, x_phys[row, col], y_phys[row, col], radius, uc,
            periodic=not config.mini_composite,
        )
        local_points = np.column_stack(
            (normal * ((radius - radial_distance) / h)[:, None], np.zeros(len(row)))
        )
        local_normals = np.column_stack((normal, np.zeros(len(row))))
        surface, phase_fractions = CVparam_2(local_points, local_normals, 0)
        ghosts = surface == 0
        phase_fractions[ghosts, 0] = 1
        phase_fractions[ghosts, 1] = 0
        phase_fractions[np.isnan(phase_fractions[:, 0]), 0] = (
            1 - phase_fractions[np.isnan(phase_fractions[:, 0]), 1]
        )
        phase_fractions[np.isnan(phase_fractions[:, 1]), 1] = (
            1 - phase_fractions[np.isnan(phase_fractions[:, 1]), 0]
        )
        positions.append(position)
        fractions.append(phase_fractions)
        areas.append(surface)
        normals.append(local_normals)
        tangents.append(np.tile([0.0, 0.0, 1.0], (len(row), 1)))

    all_positions = np.concatenate(positions) if positions else np.empty(0, dtype=np.int64)
    all_fractions = np.vstack(fractions) if fractions else np.empty((0, 2))
    all_areas = np.concatenate(areas) if areas else np.empty(0)
    all_normals = np.vstack(normals) if normals else np.empty((0, 3))
    all_tangents = np.vstack(tangents) if tangents else np.empty((0, 3))
    if len(all_positions):
        all_fractions = _adjust_fractions(all_fractions, 1e-3)
        # A voxel cut by two fibres would appear twice; AMITEX needs unique positions
        _, order = np.unique(all_positions, return_index=True)
        all_positions = all_positions[order]
        all_fractions = all_fractions[order]
        all_areas = all_areas[order]
        all_normals = all_normals[order]
        all_tangents = all_tangents[order]

    if config.mini_composite:
        solid, _ = _mini_matrix(
            x_phys, y_phys, x_centers, y_centers, radii, config, region_length
        )
        core = material == 2
        material[solid] = 1
        material[core] = 2
    else:
        solid = material == 1
        initial_pore = 1 - np.count_nonzero(solid) / material.size
        if initial_pore < config.pore_fraction:
            raise ValueError("Requested porosity exceeds the initial available pore fraction")
        solid = _grow_matrix(solid, config.pore_fraction, matrix_growth)
        material = solid.astype(np.uint8)
        material[material == 1] = 1
        # Preserve original fibre cores as phase 2 after matrix dilation.
        core = np.zeros_like(material, dtype=bool)
        for cx, cy, radius in zip(x_centers, y_centers, radii):
            dx = (x_phys - cx + length / 2) % length - length / 2
            dy = (y_phys - cy + length / 2) % length - length / 2
            core |= dx**2 + dy**2 <= radius**2
        material[core] = 2
        material[~solid] = 0

    return material, all_positions, all_fractions, all_areas, all_normals, all_tangents


def _apply_notch(
    volume: np.ndarray,
    positions: np.ndarray,
    fractions: np.ndarray,
    surfaces: np.ndarray,
    normals: np.ndarray,
    tangents: np.ndarray,
    config: Config,
    fibre: tuple[float, float, float] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Cut a notch (pore) in the mid-plane.

    "matrix" removes matrix voxels over ``notch_length`` from the x = 0 edge.
    "fibre" breaks the single fibre ``fibre = (x, y, radius)`` (mm): its voxels,
    including the interphase layer, are removed within the notch_thickness slab.
    """
    if config.notch is None or volume.ndim != 3:
        return volume, positions, fractions, surfaces, normals, tangents
    h = config.voxel_size
    nx, ny, nz = volume.shape
    z = np.arange(1, nz + 1)[None, None, :]
    in_slab = np.abs(z - (nz + 1) / 2) < config.notch_thickness / 2 / h
    if config.notch == "fibre":
        cx, cy, radius = fibre
        x_phys = np.arange(1, ny + 1)[None, :] * h
        y_phys = np.arange(1, nx + 1)[:, None] * h
        dx, dy = x_phys - cx, y_phys - cy
        if not config.mini_composite:
            dx = (dx + config.length / 2) % config.length - config.length / 2
            dy = (dy + config.length / 2) % config.length - config.length / 2
        mask = (dx**2 + dy**2 <= (radius + h) ** 2)[:, :, None] & in_slab & (volume != 0)
    else:
        x = np.arange(1, ny + 1)[None, :, None]
        mask = (x < config.notch_length / h) & in_slab & (volume == 1)
    output = volume.copy()
    output[mask] = 0
    zero_indices = np.flatnonzero(mask.ravel(order="F")) + 1
    keep = ~np.isin(positions, zero_indices)
    return output, positions[keep], fractions[keep], surfaces[keep], normals[keep], tangents[keep]


def run_workflow(
    preset_name: str,
    output: Path,
    seed: int | None = None,
    layout_file: Path | None = None,
    iuc: int | None = None,
    show_plots: bool = False,
    verbose: bool = False,
    verify_plot: bool = False,
    matrix_growth: str = "euclidean",
    overrides: dict | None = None,
) -> tuple[Path, Path]:
    """Generate a layout, discretize it, and write AMITEX input files."""
    if preset_name not in PRESETS:
        raise ValueError(f"Unknown preset {preset_name!r}; choose from {', '.join(PRESETS)}")
    config = PRESETS[preset_name]
    overrides = {k: v for k, v in (overrides or {}).items() if v is not None}
    if overrides.get("notch") == "none":
        overrides["notch"] = None
    if overrides:
        config = replace(config, **overrides)
        if verbose:
            print(f"Preset overrides: {overrides}")
    iuc = iuc or config.default_iuc
    rng = np.random.default_rng(seed)
    if layout_file:
        if verbose:
            print(f"Loading fibre layout from {layout_file}")
        x_centers, y_centers, radii = _load_layout(layout_file)
        if config.mini_composite:
            length = config.length
            fits_domain = (
                (x_centers > radii)
                & (x_centers < length - radii)
                & (y_centers > radii)
                & (y_centers < length - radii)
            )
            removed = int(np.count_nonzero(~fits_domain))
            x_centers, y_centers, radii = (
                values[fits_domain] for values in (x_centers, y_centers, radii)
            )
            if removed:
                print(f"Removed {removed} miniComp fibre(s) that touch or cross the domain boundary.")
            if len(radii) == 0:
                raise ValueError("No fibres fit fully inside the miniComp domain.")
        actual_vf = float(np.sum(np.pi * radii**2) / config.length**2)
        retries, max_overlap = 0, 0.0
    else:
        if verbose:
            print(f"Generating layout for preset {preset_name!r}")
        x_centers, y_centers, radii, actual_vf, retries, max_overlap = _generate_layout(
            config, rng, verbose
        )

    layout_directory = output / "layouts"
    layout_directory.mkdir(parents=True, exist_ok=True)
    layout_path = layout_directory / f"iUC_{iuc}.txt"
    np.savetxt(
        layout_path,
        np.column_stack((x_centers, y_centers, radii)),
        fmt="%.10g",
        header=(
            f"Fibre layout: iUC={iuc}, fibre_fraction={actual_vf:.8g}, "
            f"retries={retries}, max_overlap={max_overlap:.8g}\n"
            "x_center_mm y_center_mm radius_mm"
        ),
    )

    region_length = config.length
    if config.mini_composite:
        pad = _mini_padding(x_centers, y_centers, radii, config)
        shift = pad * config.voxel_size
        x_centers, y_centers = x_centers + shift, y_centers + shift
        config = replace(config, length=config.length + 2 * shift)
        print(
            f"miniComp domain: {config.length:.6g} mm square "
            f"(fibre region {region_length:.6g} mm, border margin {config.mini_margin:g} mm)"
        )
    volume2d, positions, fractions, surfaces, normals, tangents = _discretize(
        x_centers, y_centers, radii, config, matrix_growth, region_length
    )
    if verify_plot:
        from .verification import plot_composite_voxel_vectors

        verification_path = layout_directory / f"iUC_{iuc}_composite_voxel_check.png"
        plot_composite_voxel_vectors(
            volume2d, positions, fractions, normals, tangents,
            x_centers, y_centers, radii, config.voxel_size,
            verification_path, show_plots,
        )
        print(f"Verification plot: {verification_path}")
    depth = _extrusion_depth(config)
    volume = np.repeat(volume2d[:, :, None], depth, axis=2)
    if depth > 1 and len(positions):
        plane_size = volume.shape[0] * volume.shape[1]
        layer_offsets = plane_size * np.arange(depth, dtype=np.uint64)
        positions = (positions.astype(np.uint64)[None, :] + layer_offsets[:, None]).ravel()
        fractions = np.tile(fractions, (depth, 1))
        surfaces = np.tile(surfaces, depth)
        normals = np.tile(normals, (depth, 1))
        tangents = np.tile(tangents, (depth, 1))
    notch_fibre = None
    if config.notch == "fibre":
        index = config.notch_fibre
        if index is None:
            index = int(rng.integers(len(radii)))
        elif not 0 <= index < len(radii):
            raise ValueError(f"notch_fibre {index} out of range (0-{len(radii) - 1})")
        notch_fibre = (float(x_centers[index]), float(y_centers[index]), float(radii[index]))
        print(f"Fibre notch: breaking fibre {index} at mid-z")
    solid_before = np.count_nonzero(volume)
    volume, positions, fractions, surfaces, normals, tangents = _apply_notch(
        volume, positions, fractions, surfaces, normals, tangents,
        config, notch_fibre,
    )
    if config.notch:
        print(f"Notch ({config.notch}): removed {solid_before - np.count_nonzero(volume)} voxels")

    data_directory = output / f"iUC{iuc}" / f"vp{config.pore_fraction:g}" / f"h{config.voxel_size:g}"
    # Phase IDs (fibre, interphase, matrix); the folder name encodes them
    phase_ids = "3 4 2"
    voxel_directory = data_directory / f"rep_{phase_ids.replace(' ', '_')}"
    voxel_directory.mkdir(parents=True, exist_ok=True)
    (data_directory / "list_composite_materials.txt").write_text(f"{phase_ids} laminate\n")
    write_composite_voxel_data(
        voxel_directory, positions, fractions, surfaces, normals, tangents
    )
    suffix = "_notched" if config.notch else ""
    vtk_path = data_directory / f"iUC{iuc}_h{config.voxel_size:g}{suffix}.vtk"
    saveVolvtk_amitex(volume, vtk_path, "uint8", "cell", [config.voxel_size] * 3)
    if show_plots:
        import matplotlib.pyplot as plt

        figure, axis = plt.subplots()
        image = axis.imshow(
            volume[:, :, depth // 2], origin="lower", vmin=0, vmax=2,
            extent=(0.5, volume.shape[1] + 0.5, 0.5, volume.shape[0] + 0.5),
        )
        axis.set_aspect("equal")
        axis.set_xlabel("x voxel")
        axis.set_ylabel("y voxel")
        figure.colorbar(image, ax=axis, ticks=[0, 1, 2], label="phase (0 pore, 1 matrix, 2 fibre)")
        slice_path = layout_directory / f"iUC_{iuc}_slice.png"
        figure.savefig(slice_path, dpi=150)
        print(f"Slice plot: {slice_path}")
        plt.show()
        plt.close(figure)
    print(f"Fibre fraction: {actual_vf:.6f}; fibres: {len(radii)}; retries: {retries}")
    print(f"Layout: {layout_path}")
    print(f"Geometry: {vtk_path}")
    print(f"Composite voxels: {len(positions)} ({voxel_directory})")
    return layout_path, vtk_path


def main(default_preset: str | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--preset",
        choices=sorted(PRESETS),
        default=default_preset or "periodic3D",
        help="Case to generate (default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("python_generated"),
        help="Output directory (default: %(default)s)",
    )
    parser.add_argument(
        "--layout",
        type=Path,
        help="Reuse a fibre layout (text file, columns: x_center_mm y_center_mm radius_mm) instead of generating one",
    )
    parser.add_argument("--seed", type=int, help="Seed for reproducible random generation")
    parser.add_argument("--iuc", type=int, help="Unit-cell number in output file names")
    parser.add_argument("--plot", action="store_true", help="Show and save the generated central z slice (layouts/iUC_<n>_slice.png)")
    parser.add_argument(
        "--verbose", action="store_true", help="Report progress during iterative generation"
    )
    parser.add_argument(
        "--verify-plot",
        action="store_true",
        help="Save a plot of composite-voxel normals/tangents over the local fibre fraction",
    )
    parser.add_argument(
        "--matrix-growth",
        choices=MATRIX_GROWTH_METHODS,
        default="euclidean",
        help="Periodic matrix growth: round Euclidean distance (default) or the "
        "MATLAB-equivalent cross-shaped repeated dilation (not used by miniComp)",
    )
    group = parser.add_argument_group(
        "configuration parameters",
        "Override individual fields of the chosen preset configuration (lengths in mm).",
    )
    for number in (
        ("length", float, "Side length of the square cell"),
        ("fibre_fraction", float, "Target fibre area fraction (0-1)"),
        ("pore_fraction", float, "Target pore fraction of the matrix growth (0-1)"),
        ("voxel_size", float, "Voxel edge length"),
        ("radius_mean", float, "Mean fibre radius of the Gaussian distribution"),
        ("radius_min", float, "Smallest allowed fibre radius"),
        ("radius_max", float, "Largest allowed fibre radius"),
        ("radius_std", float, "Standard deviation of the fibre radius"),
        ("min_gap", float, "Minimum gap between neighbouring fibres"),
        ("max_trials", int, "Maximum relaxation iterations of the periodic layout generation"),
        ("z_length", float, "Extrusion length in z (2-D single plane if unset)"),
        ("radius_tolerance", float, "Relative error allowed on the fibre fraction when sampling radii"),
        ("vibration_scale", float, "Size of the random moves of the layout relaxation, in units of radius_mean"),
        ("mini_region", float, "miniComp: radius of the retained region as a fraction of half the fibre region"),
        ("mini_margin", float, "miniComp: pore border kept between fibres/matrix and the domain edge"),
        ("thickness_gradient", float, "miniComp: matrix thickness of the outermost fibre relative to the centre (1 = uniform, <1 thinner, >1 thicker outwards)"),
        ("notch_length", float, "Notch extent from the x = 0 edge"),
        ("notch_thickness", float, "Notch opening in z, centred on the mid-plane"),
        ("notch_fibre", int, "Fibre notch: 0-based index of the fibre to break (default: random)"),
    ):
        name, kind, text = number
        group.add_argument(
            "--" + name.replace("_", "-"), dest=name, type=kind, metavar=kind.__name__.upper(),
            help=text + " (default: value of the preset)",
        )
    group.add_argument(
        "--mini-composite", action=argparse.BooleanOptionalAction, default=None,
        help="Finite, non-periodic mini-composite whose fibres and matrix lie fully inside the domain, instead of a periodic cell (default: off, except in the miniComp preset)",
    )
    group.add_argument(
        "--notch", choices=("none", *NOTCH_TYPES), default=None,
        help="Notch cut in the extruded 3-D volume: matrix = matrix only, over notch-length; fibre = break one fibre at mid-z (default: none)",
    )
    args = parser.parse_args()
    overrides = {
        f.name: getattr(args, f.name)
        for f in fields(Config)
        if f.name != "default_iuc" and hasattr(args, f.name)
    }
    run_workflow(
        args.preset, args.output, args.seed, args.layout, args.iuc, args.plot,
        args.verbose, args.verify_plot, args.matrix_growth, overrides,
    )

if __name__ == "__main__":
    main()
