"""Numerical helpers translated from ``code_matlab/fcts``.

Array ordering follows MATLAB/Fortran column-major order where it matters for
AMITEX voxel indices and binary output.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.spatial import ConvexHull, QhullError


_EDGES = np.asarray(
    [
        [0, 1], [1, 2], [3, 2], [0, 3],
        [4, 5], [5, 6], [7, 6], [4, 7],
        [0, 4], [1, 5], [2, 6], [3, 7],
    ],
    dtype=int,
)
_VERTICES = np.asarray(
    [
        [-0.5, -0.5, -0.5], [0.5, -0.5, -0.5],
        [0.5, 0.5, -0.5], [-0.5, 0.5, -0.5],
        [-0.5, -0.5, 0.5], [0.5, -0.5, 0.5],
        [0.5, 0.5, 0.5], [-0.5, 0.5, 0.5],
    ],
    dtype=float,
)


def _as_points(values: Any, columns: int = 3) -> np.ndarray:
    result = np.asarray(values, dtype=float)
    if result.ndim == 1:
        result = result.reshape(1, -1)
    if result.ndim != 2 or result.shape[1] != columns:
        raise ValueError(f"Expected an array with {columns} columns")
    return result


def _periodic_delta(delta: np.ndarray, length: float) -> np.ndarray:
    return (delta + length / 2.0) % length - length / 2.0


def dist2D_periodic_v2(
    A: Sequence[float], B: Sequence[float], box: Sequence[Sequence[float]]
) -> float:
    """Periodic distance between two 2-D points, matching MATLAB v2."""
    a = np.asarray(A, dtype=float)
    b = np.asarray(B, dtype=float)
    bounds = np.asarray(box, dtype=float)
    if a.shape != (2,) or b.shape != (2,) or bounds.shape != (2, 2):
        raise ValueError("A and B must be 2-vectors and box must have shape (2, 2)")
    lengths = bounds[:, 1] - bounds[:, 0]
    dx, dy = a - b
    distances = [np.hypot(dx, dy)]
    if a[0] - bounds[0, 0] <= lengths[0] / 2:
        distances.append(np.hypot((bounds[0, 1] + (a[0] - bounds[0, 0])) - b[0], dy))
    else:
        distances.append(np.hypot((bounds[0, 0] - (bounds[0, 1] - a[0])) - b[0], dy))
    if a[1] - bounds[1, 0] <= lengths[1] / 2:
        distances.append(np.hypot(dx, (bounds[1, 1] + (a[1] - bounds[1, 0])) - b[1]))
    else:
        distances.append(np.hypot(dx, (bounds[1, 0] - (bounds[1, 1] - a[1])) - b[1]))
    return float(min(distances))


def dist2D_periodic_v3(
    A: Sequence[float],
    B: Sequence[Sequence[float]],
    box: Sequence[Sequence[float]],
    return_index: bool = False,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """Vectorized minimum-image distance; optional index identifies chosen shift."""
    a = np.asarray(A, dtype=float)
    b = _as_points(B, 2)
    bounds = np.asarray(box, dtype=float)
    if a.shape != (2,) or bounds.shape != (2, 2):
        raise ValueError("A must be a 2-vector and box must have shape (2, 2)")
    lengths = bounds[:, 1] - bounds[:, 0]
    shifts = np.asarray(
        [
            [0.0, 0.0], [lengths[0], 0.0], [-lengths[0], 0.0],
            [0.0, lengths[1]], [0.0, -lengths[1]],
            [lengths[0], lengths[1]], [lengths[0], -lengths[1]],
            [-lengths[0], -lengths[1]], [-lengths[0], lengths[1]],
        ]
    )
    delta = a[None, None, :] + shifts[None, :, :] - b[:, None, :]
    candidates = np.linalg.norm(delta, axis=2)
    indices = np.argmin(candidates, axis=1)
    distances = candidates[np.arange(len(b)), indices]
    if return_index:
        return distances, indices
    return distances


def getCmax(
    x0_all: Sequence[float],
    y0_all: Sequence[float],
    r0_all: Sequence[float],
    eps: float,
    ROI: Sequence[Sequence[float]],
) -> np.ndarray:
    """Compute pairwise circle-overlap severity (rows are circles)."""
    centers = np.column_stack((x0_all, y0_all)).astype(float)
    radii = np.asarray(r0_all, dtype=float).reshape(-1)
    if len(centers) != len(radii):
        raise ValueError("Center and radius arrays must have equal lengths")
    result = np.zeros((len(radii), max(0, len(radii) - 1)), dtype=float)
    for i in range(len(radii)):
        others = np.arange(len(radii)) != i
        distances = dist2D_periodic_v3(centers[i], centers[others], ROI)
        required = radii[i] + radii[others] + eps
        result[i] = (required - distances) / required
    return result


def relaxCircle(
    x0: float,
    y0: float,
    Cij: Sequence[float],
    x1: Sequence[float],
    y1: Sequence[float],
    ind: Sequence[int],
    ROI: Sequence[Sequence[float]],
) -> tuple[float, float]:
    """Move one center in response to repulsive overlap forces."""
    bounds = np.asarray(ROI, dtype=float)
    overlap = np.asarray(Cij, dtype=float)
    indices = np.asarray(ind, dtype=int)
    x1a, y1a = np.asarray(x1, dtype=float), np.asarray(y1, dtype=float)
    nx = (x0 - x1a) / (1.0 - overlap) * overlap
    ny = (y0 - y1a) / (1.0 - overlap) * overlap
    wrap_x = (indices == 1) | (indices == 2) | (indices >= 5)
    wrap_y = (indices == 3) | (indices == 4) | (indices >= 5)
    nx[wrap_x] *= -1
    ny[wrap_y] *= -1
    x_new, y_new = x0 + np.sum(nx) / 3.0, y0 + np.sum(ny) / 3.0
    lx, ly = bounds[0, 1] - bounds[0, 0], bounds[1, 1] - bounds[1, 0]
    x_new = bounds[0, 0] + (x_new - bounds[0, 0]) % lx
    y_new = bounds[1, 0] + (y_new - bounds[1, 0]) % ly
    return float(x_new), float(y_new)


def relaxationCircles(
    x0_all: Sequence[float],
    y0_all: Sequence[float],
    r0_all: Sequence[float],
    id: Sequence[bool],
    eps: float,
    ROI: Sequence[Sequence[float]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply a relaxation step to overlapping periodic circles."""
    del id  # The MATLAB implementation intentionally tests every circle.
    x, y, radii = (np.asarray(v, dtype=float).copy() for v in (x0_all, y0_all, r0_all))
    relaxed = np.zeros(len(radii), dtype=bool)
    current_max = float(np.max(getCmax(x, y, radii, eps, ROI), initial=0.0))
    for i in range(len(radii)):
        other_ids = np.arange(len(radii)) != i
        distances, shifts = dist2D_periodic_v3(
            [x[i], y[i]], np.column_stack((x[other_ids], y[other_ids])), ROI, True
        )
        required = radii[i] + radii[other_ids] + eps
        overlaps = (required - distances) / required
        active = overlaps > 0
        if not np.any(active):
            continue
        new_x, new_y = relaxCircle(
            x[i], y[i], overlaps[active], x[other_ids][active],
            y[other_ids][active], shifts[active], ROI
        )
        after = dist2D_periodic_v3(
            [new_x, new_y], np.column_stack((x[other_ids], y[other_ids])), ROI
        )
        if np.max((required - after) / required, initial=-np.inf) < current_max:
            x[i], y[i], relaxed[i] = new_x, new_y, True
    return x, y, relaxed


def vibrationCircles(
    x0_all: Sequence[float],
    y0_all: Sequence[float],
    r0_all: Sequence[float],
    eps: float,
    ROI: Sequence[Sequence[float]],
    lstRelaxed: Sequence[bool],
    ur_sc: float,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Randomly perturb relaxed circles that remain overlapped."""
    rng = rng or np.random.default_rng()
    bounds = np.asarray(ROI, dtype=float)
    x, y, radii = (np.asarray(v, dtype=float).copy() for v in (x0_all, y0_all, r0_all))
    for i in np.flatnonzero(lstRelaxed):
        others = np.arange(len(radii)) != i
        distances = dist2D_periodic_v3([x[i], y[i]], np.column_stack((x[others], y[others])), ROI)
        if np.any(radii[i] + radii[others] + eps > distances):
            radius = rng.random() * ur_sc
            angle = rng.random() * 2 * np.pi
            x[i] += radius * np.cos(angle)
            y[i] += radius * np.sin(angle)
            x[i] = bounds[0, 0] + (x[i] - bounds[0, 0]) % (bounds[0, 1] - bounds[0, 0])
            y[i] = bounds[1, 0] + (y[i] - bounds[1, 0]) % (bounds[1, 1] - bounds[1, 0])
    return x, y


def generationFibres(
    r0_all: Sequence[float],
    param: Mapping[str, float],
    ROI: Sequence[Sequence[float]],
    rng: np.random.Generator | None = None,
    verbose: bool = False,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Place circles without overlap using relaxation and random vibration."""
    rng = rng or np.random.default_rng()
    radii = np.asarray(r0_all, dtype=float).reshape(-1)
    bounds = np.asarray(ROI, dtype=float)
    x = bounds[0, 0] + rng.random(len(radii)) * (bounds[0, 1] - bounds[0, 0])
    y = bounds[1, 0] + rng.random(len(radii)) * (bounds[1, 1] - bounds[1, 0])
    eps, max_trials = float(param["eps"]), int(param["trialsMAX"])
    overlap = getCmax(x, y, radii, eps, bounds)
    current = float(np.max(overlap, initial=0.0))
    if verbose:
        print(
            f"  Relaxation: initial maximum overlap={current:.6g}; "
            f"overlapping fibres={np.count_nonzero(np.max(overlap, axis=1) > 0)}"
        )
        if current <= 0:
            print("  Relaxation converged before the first iteration")
    for trial in range(max_trials + 1):
        if current <= 0:
            break
        x, y, relaxed = relaxationCircles(x, y, radii, np.max(overlap, axis=1) > 0, eps, bounds)
        x, y = vibrationCircles(x, y, radii, eps, bounds, relaxed, float(param["ur_sc"]), rng)
        overlap = getCmax(x, y, radii, eps, bounds)
        current = float(np.max(overlap, initial=0.0))
        if verbose:
            print(
                f"  Relaxation iteration {trial + 1}/{max_trials + 1}: "
                f"maximum overlap={current:.6g}; "
                f"overlapping fibres={np.count_nonzero(np.max(overlap, axis=1) > 0)}"
            )
    if verbose and current > 0:
        print(f"  Relaxation stopped with unresolved overlap={current:.6g}")
    return x, y, current


def GenDataWithGaussianDistrib(
    r: Mapping[str, float], rng: np.random.Generator | None = None
) -> np.ndarray:
    """Draw radii from a Gaussian distribution truncated to ``MIN..MAX``."""
    from scipy.special import ndtr

    rng = rng or np.random.default_rng()
    n = int(r["nn"])
    lower = ndtr((float(r["MIN"]) - float(r["AV"])) / float(r["SIG"]))
    upper = ndtr((float(r["MAX"]) - float(r["AV"])) / float(r["SIG"]))
    samples = rng.uniform(lower, upper, n)
    from scipy.special import ndtri

    return float(r["AV"]) + float(r["SIG"]) * ndtri(samples)


def getCompVox_circl_periodic(
    x0: float,
    y0: float,
    x: np.ndarray,
    y: np.ndarray,
    r0: float,
    UC: Mapping[str, float],
    periodic: bool = True,
) -> np.ndarray:
    """Find voxels whose four corners lie on both sides of a circular edge.

    Set ``periodic=False`` for a finite, non-periodic material domain.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.shape != y.shape:
        raise ValueError("x and y grids must have the same shape")
    signs = np.zeros((x.size, 4), dtype=np.int8)
    zero = np.zeros((x.size, 4), dtype=bool)
    corners = [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]
    h = float(UC["vxsiz"])
    box = [[UC["xlim0"], UC["xlim1"]], [UC["ylim0"], UC["ylim1"]]]
    for i, (cx, cy) in enumerate(corners):
        corners_x, corners_y = (x + cx * h).ravel(), (y + cy * h).ravel()
        if periodic:
            d = dist2D_periodic_v3(
                [x0, y0], np.column_stack((corners_x, corners_y)), box
            )
        else:
            d = np.hypot(corners_x - x0, corners_y - y0)
        signs[d >= r0, i] = 1
        signs[d < r0, i] = -1
        zero[:, i] = d == r0
    denom = np.sum(~zero, axis=1)
    mixed = np.divide(np.abs(np.sum(signs, axis=1)), denom, out=np.zeros_like(denom, dtype=float), where=denom != 0) < 1
    return mixed.reshape(x.shape)


def getOrientVec_periodic(
    x0: float,
    y0: float,
    iX: Sequence[float],
    iY: Sequence[float],
    r0: float,
    UC: Mapping[str, float],
    periodic: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Return outward circle normals and center-to-voxel distances.

    Set ``periodic=False`` for a finite, non-periodic material domain.
    """
    points = np.column_stack((iX, iY)).astype(float)
    center = np.asarray([x0, y0], dtype=float)
    bounds = np.asarray(
        [[UC["xlim0"], UC["xlim1"]], [UC["ylim0"], UC["ylim1"]]], dtype=float
    )
    lengths = bounds[:, 1] - bounds[:, 0]
    delta = points - center
    if periodic:
        outliers = np.linalg.norm(delta, axis=1) - r0 > float(UC["vxsiz"])
        delta[outliers] = _periodic_delta(delta[outliers], lengths)
    radii = np.linalg.norm(delta, axis=1)
    if np.any(radii == 0):
        raise ValueError("A voxel center coincides with a fibre center")
    return delta / radii[:, None], radii


def findInterSectVertices(
    P: Sequence[Sequence[float]],
    N: Sequence[Sequence[float]],
    xV: Sequence[Sequence[float]],
    yV: Sequence[Sequence[float]],
    zV: Sequence[Sequence[float]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Intersect plane(s) with the 12 edges of each cuboid."""
    p, n = _as_points(P), _as_points(N)
    vertices = np.stack((xV, yV, zV), axis=-1).astype(float)
    if vertices.ndim == 2:
        vertices = vertices[None, ...]
    if vertices.shape != (len(n), 8, 3) or len(p) != len(n):
        raise ValueError("Expected P/N shape (n,3) and vertex coordinate arrays (n,8)")
    result = np.full((len(n), 12, 3), np.nan, dtype=float)
    for edge_index, (i, j) in enumerate(_EDGES):
        v1, v2 = vertices[:, i], vertices[:, j]
        direction = v2 - v1
        numerator = np.sum(n * (p - v1), axis=1)
        denominator = np.sum(direction * n, axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            point = v1 + direction * (numerator / denominator)[:, None]
        valid = np.all((point - v1 >= -1e-12) & (point - v2 <= 1e-12), axis=1)
        valid &= np.isfinite(point).all(axis=1)
        result[valid, edge_index] = point[valid]
    return result[:, :, 0], result[:, :, 1], result[:, :, 2]


def polyarea3D(x: Sequence[float], y: Sequence[float], z: Sequence[float]) -> float:
    """Compute area of a planar 3-D polygon, tolerating repeated vertices."""
    points = np.unique(np.column_stack((x, y, z)).astype(float), axis=0)
    if len(points) < 3:
        return 0.0
    center = points.mean(axis=0)
    centered = points - center
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    projected = centered @ vh[:2].T
    angles = np.arctan2(projected[:, 1], projected[:, 0])
    polygon = projected[np.argsort(angles)]
    return float(
        0.5
        * abs(
            np.dot(polygon[:, 0], np.roll(polygon[:, 1], -1))
            - np.dot(polygon[:, 1], np.roll(polygon[:, 0], -1))
        )
    )


def area_element_2(
    N: Sequence[Sequence[float]],
    G: Sequence[Sequence[float]],
    xV: Sequence[Sequence[float]],
    yV: Sequence[Sequence[float]],
    zV: Sequence[Sequence[float]],
) -> np.ndarray:
    """Area where each plane intersects its cuboid."""
    n = _as_points(N)
    g = _as_points(G)
    xi, yi, zi = findInterSectVertices(g, n, xV, yV, zV)
    areas = np.zeros(len(n), dtype=float)
    for i in range(len(n)):
        valid = np.isfinite(xi[i])
        areas[i] = polyarea3D(xi[i, valid], yi[i, valid], zi[i, valid])
    return areas


def _plane_intersections_one(
    point: np.ndarray, normal: np.ndarray, vertices: np.ndarray
) -> np.ndarray:
    intersections: list[np.ndarray] = []
    for i, j in _EDGES:
        a, b = vertices[i], vertices[j]
        denom = float(np.dot(normal, b - a))
        if abs(denom) < 1e-14:
            continue
        t = float(np.dot(normal, point - a) / denom)
        if -1e-12 <= t <= 1.0 + 1e-12:
            intersections.append(a + np.clip(t, 0.0, 1.0) * (b - a))
    return np.unique(np.asarray(intersections), axis=0) if intersections else np.empty((0, 3))


def volume_element(
    N: Sequence[Sequence[float]],
    P: Sequence[Sequence[float]],
    xV: Sequence[Sequence[float]],
    yV: Sequence[Sequence[float]],
    zV: Sequence[Sequence[float]],
    eI: float,
) -> np.ndarray:
    """Volume fractions of the two cuboid regions outside a finite interface."""
    normals, points = _as_points(N), _as_points(P)
    verts = np.stack((xV, yV, zV), axis=-1).astype(float)
    if verts.ndim == 2:
        verts = verts[None, ...]
    if verts.shape != (len(normals), 8, 3):
        raise ValueError("Vertex coordinates must have shape (n, 8)")
    result = np.zeros((len(normals), 2), dtype=float)
    for i, (normal, point, vertices) in enumerate(zip(normals, points, verts)):
        normal_length = np.linalg.norm(normal)
        if normal_length == 0:
            raise ValueError("Plane normals must be non-zero")
        normal = normal / normal_length
        # The source MATLAB routine positions the two faces at P +/- eI*N.
        cuts = []
        for shifted_point in (point - eI * normal, point + eI * normal):
            cuts.append(_plane_intersections_one(shifted_point, normal, vertices))
        sides = (vertices - point) @ normal > 0
        for phase in range(2):
            face = cuts[phase]
            if len(face) < 3:
                result[i, phase] = 0.0
                continue
            cap_vertices = np.vstack((face, vertices[sides == (phase == 1)]))
            cap_vertices = np.unique(cap_vertices, axis=0)
            if len(cap_vertices) < 4:
                result[i, phase] = 0.0
                continue
            try:
                result[i, phase] = ConvexHull(cap_vertices).volume
            except QhullError:
                result[i, phase] = np.nan
    return result


def CVparam_2(
    P: Sequence[Sequence[float]], N: Sequence[Sequence[float]], eI: float
) -> tuple[np.ndarray, np.ndarray]:
    """Compute interface area and two bulk volume fractions for unit voxels."""
    points, normals = _as_points(P), _as_points(N)
    vertices = np.tile(_VERTICES, (len(points), 1, 1))
    x, y, z = (vertices[:, :, i] for i in range(3))
    areas = area_element_2(normals, points, x, y, z)
    fractions = volume_element(normals, points, x, y, z, float(eI))
    return areas, fractions


def CVparam(
    V: np.ndarray, theta: float, fvI: float, matID: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Legacy straight-interface voxel parameters (translated MATLAB CVparam)."""
    volume = np.asarray(V)
    if volume.ndim != 2:
        raise ValueError("V must be a 2-D material-ID array")
    rows, cols = volume.shape
    y, x = np.mgrid[1 : rows + 1, 1 : cols + 1]
    selected = volume.ravel(order="F") == matID
    positions = np.flatnonzero(selected)
    count = int(np.count_nonzero(selected))
    tangent = np.tile([np.cos(theta), np.sin(theta), 0.0], (count, 1))
    normal = np.tile([np.cos(theta + np.pi / 2), np.sin(theta + np.pi / 2), 0.0], (count, 1))
    phase1 = np.full(count, (1 - fvI) / 2.0)
    phase2 = np.full(count, (1 - fvI) / 2.0)
    phase3 = np.full(count, fvI)
    zone = np.ones(count, dtype=np.uint64)
    return tangent, normal, phase1, phase2, phase3, positions + 1, zone


def GetStiffIso(E: float, NU: float) -> np.ndarray:
    """Isotropic stiffness in the 9-entry engineering-Voigt layout."""
    lame_lambda = E * NU / ((1 + NU) * (1 - 2 * NU))
    shear = E / (2 * (1 + NU))
    return np.asarray(
        [
            lame_lambda + 2 * shear, lame_lambda, lame_lambda,
            lame_lambda + 2 * shear, lame_lambda, lame_lambda + 2 * shear,
            shear, shear, shear,
        ],
        dtype=float,
    )


def GetStress(C: Sequence[float], STRANT: Sequence[float]) -> np.ndarray:
    """Evaluate stress using the 9-entry stiffness format."""
    c = np.asarray(C, dtype=float)
    strain = np.asarray(STRANT, dtype=float)
    if c.size != 9 or strain.size != 6:
        raise ValueError("C must contain 9 values and STRANT must contain 6")
    return np.asarray(
        [
            c[0] * strain[0] + c[1] * strain[1] + c[2] * strain[2],
            c[1] * strain[0] + c[3] * strain[1] + c[4] * strain[2],
            c[2] * strain[0] + c[4] * strain[1] + c[5] * strain[2],
            c[6] * strain[3], c[7] * strain[4], c[8] * strain[5],
        ]
    )


def CDM_interphase_fct(
    lbda: float, mu: float, gc: float, lc: float, def_: Sequence[float]
) -> tuple[np.ndarray, float]:
    """Phase-field damaged isotropic stress and scalar damage."""
    strain = np.asarray(def_, dtype=float).reshape(6)
    positive = (strain + np.abs(strain)) / 2
    negative = (strain - np.abs(strain)) / 2
    trace = np.sum(strain[:3])
    trace_positive, trace_negative = max(trace, 0), min(trace, 0)
    energy = (
        0.5 * lbda * trace_positive**2
        + mu * (np.sum(positive[:3] ** 2) + np.sum(positive[3:] ** 2) / 2)
    )
    damage = 1 - gc / (gc + 2 * lc * energy)
    stress = np.zeros(6)
    stress[:3] = (1 - damage) ** 2 * (lbda * trace_positive + 2 * mu * positive[:3])
    stress[:3] += lbda * trace_negative + 2 * mu * negative[:3]
    stress[3:] = (1 - damage) ** 2 * mu * positive[3:] + mu * negative[3:]
    return stress, float(damage)


def plot_circle(x: float, y: float, r: float, c: str = "b", points: int = 200) -> Any:
    """Plot one circle using matplotlib and return the line handle."""
    import matplotlib.pyplot as plt

    theta = np.linspace(0, 2 * np.pi, points)
    return plt.plot(x + r * np.cos(theta), y + r * np.sin(theta), c)


def plotPeriodicCircles(
    x0_all: Sequence[float],
    y0_all: Sequence[float],
    r0_all: Sequence[float],
    UC: Mapping[str, float],
    *args: str,
) -> Any:
    """Plot periodic copies of circles within the configured unit cell."""
    import matplotlib.pyplot as plt

    offset = bool(args and args[0].lower() == "offsetim")
    xlo, xhi = float(UC["xlim0"]), float(UC["xlim1"])
    ylo, yhi = float(UC["ylim0"]), float(UC["ylim1"])
    h = float(UC.get("vxsiz", 1.0)) if offset else 1.0
    artists = []
    for x, y, r in zip(x0_all, y0_all, r0_all):
        for dx in (-xhi + xlo, 0, xhi - xlo):
            for dy in (-yhi + ylo, 0, yhi - ylo):
                artists.extend(plot_circle(float(x) + dx / h, float(y) + dy / h, float(r), "b"))
    plt.xlim(xlo, xhi)
    plt.ylim(ylo, yhi)
    plt.gca().set_aspect("equal")
    return artists


def saveBin(data: Sequence[Any], filename: str | Path, dtype: str | None = None) -> None:
    """Write the legacy MATLAB binary format (count/type header, big-endian data)."""
    array = np.asarray(data)
    type_map = {
        "double": ("double", ">f8"),
        "single": ("float", ">f4"),
        "uint8": ("unsigned_char", "u1"),
        "int8": ("char", "i1"),
        "uint16": ("unsigned_short", ">u2"),
        "int16": ("short", ">i2"),
        "uint32": ("unsigned_int", ">u4"),
        "int32": ("int", ">i4"),
        "uint64": ("unsigne_long", ">u8"),
        "int64": ("long", ">i8"),
    }
    kind = dtype or ("double" if array.dtype.kind == "f" else array.dtype.name)
    if kind not in type_map:
        raise ValueError(f"Unsupported saveBin dtype: {kind}")
    vtk_type, numpy_type = type_map[kind]
    flat = array.ravel(order="F").astype(numpy_type, copy=False)
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as stream:
        stream.write(f"{flat.size}\n{vtk_type}\n".encode("ascii"))
        flat.tofile(stream)


def saveVolvtk_amitex(
    volume: np.ndarray,
    filename: str | Path,
    dtype: str = "uint8",
    choice: str = "cell",
    spacing: Sequence[float] = (1.0, 1.0, 1.0),
) -> None:
    """Write a legacy binary STRUCTURED_POINTS VTK volume for AMITEX."""
    type_map = {
        "uint8": ("unsigned_char", ">u1"),
        "int8": ("char", ">i1"),
        "uint16": ("unsigned_short", ">u2"),
        "int16": ("short", ">i2"),
        "uint32": ("unsigned_int", ">u4"),
        "int32": ("int", ">i4"),
        "uint64": ("unsigned_long", ">u8"),
        "int64": ("long", ">i8"),
        "single": ("float", ">f4"),
        "double": ("double", ">f8"),
    }
    data = np.asarray(volume)
    if data.ndim == 2:
        data = data[:, :, None]
    if data.ndim != 3 or dtype not in type_map or choice not in {"cell", "point"}:
        raise ValueError("Expected a 3-D volume, supported VTK dtype, and cell or point data")
    vtk_type, numpy_type = type_map[dtype]
    nx, ny, nz = data.shape
    dimensions = (nx + 1, ny + 1, nz + 1) if choice == "cell" else (nx, ny, nz)
    data_count = nx * ny * nz
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as stream:
        header = (
            "# vtk DataFile Version 4.2\n"
            "tomo_Volume\n"
            "BINARY\n"
            "DATASET STRUCTURED_POINTS\n"
            f"DIMENSIONS {dimensions[0]} {dimensions[1]} {dimensions[2]}\n"
            "ORIGIN 0 0 0\n"
            f"SPACING {float(spacing[0]):e} {float(spacing[1]):e} {float(spacing[2]):e}\n"
            f"{'CELL_DATA' if choice == 'cell' else 'POINT_DATA'} {data_count}\n"
            f"SCALARS tomo_Volume {vtk_type}\n"
            "LOOKUP_TABLE default\n"
        )
        stream.write(header.encode("ascii"))
        data.ravel(order="F").astype(numpy_type, copy=False).tofile(stream)


def write_composite_voxel_data(
    directory: str | Path,
    positions: np.ndarray,
    fractions: np.ndarray,
    surfaces: np.ndarray,
    normals: np.ndarray,
    tangents: np.ndarray,
) -> None:
    """Write the AMITEX composite-voxel files produced by the MATLAB scripts."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    count = len(positions)
    if fractions.shape != (count, 2) or normals.shape != (count, 3) or tangents.shape != (count, 3):
        raise ValueError("Composite voxel arrays have inconsistent lengths")
    if surfaces.shape not in {(count,), (count, 1)}:
        raise ValueError("Surface array has inconsistent length")
    values = {
        "fv1.bin": (fractions[:, 0], "double"),
        "fv2.bin": (1.0 - np.sum(fractions, axis=1), "double"),
        "fv3.bin": (fractions[:, 1], "double"),
        "N12x.bin": (normals[:, 1], "double"),
        "N12y.bin": (normals[:, 0], "double"),
        "N12z.bin": (np.zeros(count), "double"),
        "Tx.bin": (tangents[:, 1], "double"),
        "Ty.bin": (tangents[:, 0], "double"),
        "Tz.bin": (tangents[:, 2], "double"),
        "S12.bin": (np.asarray(surfaces).reshape(-1), "double"),
        "pos.bin": (np.asarray(positions, dtype=np.uint64), "uint64"),
        "zone1.bin": (np.ones(count, dtype=np.uint64), "uint64"),
        "zone2.bin": (np.ones(count, dtype=np.uint64), "uint64"),
        "zone3.bin": (np.ones(count, dtype=np.uint64), "uint64"),
    }
    for name, (data, kind) in values.items():
        saveBin(data, directory / name, kind)
