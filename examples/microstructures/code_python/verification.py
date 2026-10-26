"""Visual check of composite-voxel normals and tangents."""

from __future__ import annotations

from pathlib import Path

import numpy as np

FIBRE_PHASE = 2  # material id of fibres in the voxel image
FIBRE_COLUMN = 0  # fibre column in the composite-voxel volume fractions


def fibre_fraction_image(
    material: np.ndarray, positions: np.ndarray, fractions: np.ndarray
) -> np.ndarray:
    """Fibre fraction per voxel: 0/1 for pure voxels, local value for composite voxels."""
    image = (material == FIBRE_PHASE).astype(float)
    rows, cols = _row_col(positions, material.shape[0])
    image[rows, cols] = fractions[:, FIBRE_COLUMN]
    return image


def _row_col(positions: np.ndarray, n_rows: int) -> tuple[np.ndarray, np.ndarray]:
    index = np.asarray(positions, dtype=np.int64) - 1  # one-based, column-major
    return index % n_rows, index // n_rows


def plot_composite_voxel_vectors(
    material: np.ndarray,
    positions: np.ndarray,
    fractions: np.ndarray,
    normals: np.ndarray,
    tangents: np.ndarray,
    x_centers: np.ndarray,
    y_centers: np.ndarray,
    radii: np.ndarray,
    voxel_size: float,
    output: Path | None = None,
    show: bool = False,
    max_arrows: int = 1500,
):
    """Overlay normals and tangents on the local fibre-fraction image.

    The stored tangent is the out-of-plane fibre axis, so it is drawn as dots;
    the in-plane tangent (fibre axis x normal) is drawn as arrows for checking
    that the normal is perpendicular to the interface.
    """
    import matplotlib

    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_rows, n_cols = material.shape
    rows, cols = _row_col(positions, n_rows)
    image = fibre_fraction_image(material, positions, fractions)
    normals = np.asarray(normals, dtype=float)[:, :2]
    tangents = np.asarray(tangents, dtype=float)
    in_plane_tangent = np.column_stack((-normals[:, 1], normals[:, 0]))

    # Normals should point from the fibre centre towards the matrix.
    centers = np.column_stack((x_centers, y_centers)) / voxel_size
    voxel_xy = np.column_stack((cols + 1, rows + 1)).astype(float)
    nearest = np.argmin(
        np.linalg.norm(voxel_xy[:, None, :] - centers[None, :, :], axis=2)
        - np.asarray(radii)[None, :] / voxel_size,
        axis=1,
    )
    radial = voxel_xy - centers[nearest]
    radial /= np.maximum(np.linalg.norm(radial, axis=1, keepdims=True), 1e-12)
    alignment = np.sum(radial * normals, axis=1)
    print(
        "Composite-voxel check: "
        f"{len(positions)} voxels; "
        f"|n| in [{np.linalg.norm(normals, axis=1).min():.3f}, "
        f"{np.linalg.norm(normals, axis=1).max():.3f}]; "
        f"mean n.radial={alignment.mean():.3f} (min {alignment.min():.3f}); "
        f"max |n.t|={np.abs(np.sum(normals * in_plane_tangent, axis=1)).max():.1e}; "
        f"stored tangent = {np.unique(tangents, axis=0).tolist()}"
    )

    stride = max(1, int(np.ceil(len(positions) / max_arrows)))
    pick = slice(None, None, stride)
    zoom = np.argmax(radii)
    zoom_half = 2.0 * radii[zoom] / voxel_size

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    for ax, dense in zip(axes, (False, True)):
        mesh = ax.imshow(
            image, origin="lower", extent=(0.5, n_cols + 0.5, 0.5, n_rows + 0.5),
            cmap="viridis", vmin=0, vmax=1, interpolation="nearest",
        )
        sel = slice(None) if dense else pick
        scale = 25 if dense else 40
        ax.quiver(
            cols[sel] + 1, rows[sel] + 1, normals[sel, 0], normals[sel, 1],
            color="red", angles="xy", scale_units="xy", scale=1 / (scale / 25),
            width=0.003, label="normal",
        )
        ax.quiver(
            cols[sel] + 1, rows[sel] + 1, in_plane_tangent[sel, 0],
            in_plane_tangent[sel, 1], color="white", angles="xy", scale_units="xy",
            scale=1 / (scale / 25), width=0.003, label="tangent (in-plane)",
        )
        ax.plot(
            cols[sel] + 1, rows[sel] + 1, ".", color="cyan", markersize=2,
            label="stored tangent (out of plane, +z)",
        )
        ax.set_xlabel("x voxel")
        ax.set_ylabel("y voxel")
        ax.set_aspect("equal")
        if dense:
            cx, cy = centers[zoom]
            ax.set_xlim(cx - zoom_half, cx + zoom_half)
            ax.set_ylim(cy - zoom_half, cy + zoom_half)
            ax.set_title("Zoom on the largest fibre (all composite voxels)")
        else:
            ax.set_title("Normals and tangents of composite voxels")
            ax.legend(loc="upper right", fontsize=7)
    fig.colorbar(mesh, ax=axes, label="local fibre fraction", shrink=0.8)
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return alignment
