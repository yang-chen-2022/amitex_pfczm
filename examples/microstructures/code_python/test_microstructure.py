import contextlib
from dataclasses import replace
import io
import tempfile
import unittest
from pathlib import Path

import numpy as np

from .microstructure import (
    CVparam_2,
    dist2D_periodic_v3,
    getCompVox_circl_periodic,
    getOrientVec_periodic,
    generationFibres,
    saveBin,
    saveVolvtk_amitex,
)
from .workflows import (
    PRESETS,
    _extrusion_depth,
    _generate_nonperiodic_layout,
    _apply_notch,
    _discretize,
    _grow_matrix,
    _load_layout,
    _mini_padding,
    _thickness_factors,
    run_workflow,
)
from .verification import plot_composite_voxel_vectors


class MicrostructureTests(unittest.TestCase):
    def test_periodic_distance_uses_minimum_image(self):
        distance = dist2D_periodic_v3(
            [0.02, 0.5], [[0.98, 0.5]], [[0, 1], [0, 1]]
        )
        np.testing.assert_allclose(distance, [0.04])

    def test_plane_through_unit_cube_bisects_volume(self):
        area, fractions = CVparam_2([[0, 0, 0]], [[1, 0, 0]], 0)
        np.testing.assert_allclose(area, [1.0])
        np.testing.assert_allclose(fractions, [[0.5, 0.5]])

    def test_plane_volume_fractions_sum_to_cube_volume(self):
        _, fractions = CVparam_2([[0.13, -0.1, 0.2]], [[1, 2, 3]], 0)
        self.assertAlmostEqual(float(np.sum(fractions)), 1.0)

    def test_legacy_binary_file_is_column_major_and_big_endian(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "pos.bin"
            saveBin(np.asarray([[1, 2], [3, 4]]), path, "uint64")
            raw = path.read_bytes()
            header, payload = raw.split(b"\n", 2)[:2], raw.split(b"\n", 2)[2]
            self.assertEqual(header, [b"4", b"unsigne_long"])
            np.testing.assert_array_equal(np.frombuffer(payload, dtype=">u8"), [1, 3, 2, 4])

    def test_vtk_writer_emits_valid_structured_points_header(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "volume.vtk"
            saveVolvtk_amitex(np.zeros((2, 3, 1), dtype=np.uint8), path, spacing=[0.1] * 3)
            contents = path.read_bytes()
            self.assertIn(b"DIMENSIONS 3 4 2\n", contents)
            self.assertIn(b"CELL_DATA 6\n", contents)

    def test_periodic2d_is_a_single_z_plane(self):
        self.assertEqual(_extrusion_depth(PRESETS["periodic2D"]), 1)
        self.assertEqual(_extrusion_depth(PRESETS["miniComp"]), 101)

    def test_text_layout_round_trip(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "layout.txt"
            np.savetxt(path, [[0.02, 0.03, 0.004], [0.06, 0.05, 0.003]], header="x y r")
            x, y, radius = _load_layout(path)
            np.testing.assert_allclose(x, [0.02, 0.06])
            np.testing.assert_allclose(y, [0.03, 0.05])
            np.testing.assert_allclose(radius, [0.004, 0.003])

    def test_mini_comp_layout_keeps_fibres_inside_nonperiodic_cell(self):
        preset = PRESETS["miniComp"]
        x, y, radii, _ = _generate_nonperiodic_layout(
            preset, np.random.default_rng(42)
        )
        self.assertTrue(np.all(x > radii))
        self.assertTrue(np.all(x < preset.length - radii))
        self.assertTrue(np.all(y > radii))
        self.assertTrue(np.all(y < preset.length - radii))
        for i in range(len(radii)):
            for j in range(i):
                self.assertGreaterEqual(
                    np.hypot(x[i] - x[j], y[i] - y[j]),
                    radii[i] + radii[j] + preset.min_gap,
                )

    def test_orientation_can_disable_periodic_wrapping(self):
        uc = {"xlim0": 0, "xlim1": 1, "ylim0": 0, "ylim1": 1, "vxsiz": 0.01}
        periodic_normal, _ = getOrientVec_periodic(
            0.01, 0.5, [0.99], [0.5], 0.1, uc, periodic=True
        )
        finite_normal, _ = getOrientVec_periodic(
            0.01, 0.5, [0.99], [0.5], 0.1, uc, periodic=False
        )
        self.assertLess(periodic_normal[0, 0], 0)
        self.assertGreater(finite_normal[0, 0], 0)

    def test_nonperiodic_interface_does_not_wrap_across_domain(self):
        coordinates = np.arange(1, 11, dtype=float) / 10
        y_grid, x_grid = np.meshgrid(coordinates, coordinates)
        uc = {"xlim0": 0, "xlim1": 1, "ylim0": 0, "ylim1": 1, "vxsiz": 0.1}
        finite = getCompVox_circl_periodic(
            0.1, 0.5, x_grid, y_grid, 0.15, uc, periodic=False
        )
        periodic = getCompVox_circl_periodic(
            0.1, 0.5, x_grid, y_grid, 0.15, uc, periodic=True
        )
        self.assertFalse(np.any(finite[-1, :]))
        self.assertTrue(np.any(periodic[-1, :]))

    def test_fibre_generation_can_report_iteration_progress(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            generationFibres(
                [0.45, 0.45],
                {"eps": 0.0, "trialsMAX": 0, "ur_sc": 0.01},
                [[0, 1], [0, 1]],
                np.random.default_rng(1),
                verbose=True,
            )
        self.assertIn("Relaxation: initial maximum overlap=", output.getvalue())
        self.assertIn("Relaxation iteration 1/1:", output.getvalue())

    def test_composite_voxel_normals_point_out_of_fibres(self):
        preset = PRESETS["miniComp"]
        x, y, radii, _ = _generate_nonperiodic_layout(preset, np.random.default_rng(3))
        material, positions, fractions, surfaces, normals, _ = _discretize(
            x, y, radii, preset
        )
        n_rows = material.shape[0]
        rows, cols = (positions - 1) % n_rows, (positions - 1) // n_rows
        centers = np.column_stack((x, y)) / preset.voxel_size
        voxels = np.column_stack((cols + 1, rows + 1))
        nearest = np.argmin(
            np.linalg.norm(voxels[:, None] - centers[None], axis=2)
            - radii[None] / preset.voxel_size,
            axis=1,
        )
        radial = voxels - centers[nearest]
        radial /= np.linalg.norm(radial, axis=1, keepdims=True)
        self.assertGreater(np.mean(np.sum(radial * normals[:, :2], axis=1)), 0.99)
        self.assertGreater(np.mean(surfaces > 0), 0.9)

    def test_verification_plot_is_written(self):
        preset = PRESETS["miniComp"]
        x, y, radii, _ = _generate_nonperiodic_layout(preset, np.random.default_rng(3))
        data = _discretize(x, y, radii, preset)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "check.png"
            with contextlib.redirect_stdout(io.StringIO()):
                plot_composite_voxel_vectors(
                    data[0], data[1], data[2], data[4], data[5],
                    x, y, radii, preset.voxel_size, path,
                )
            self.assertGreater(path.stat().st_size, 0)

    def test_euclidean_growth_is_round_periodic_and_hits_target(self):
        solid = np.zeros((101, 101), dtype=bool)
        solid[48:53, 48:53] = True
        grown = _grow_matrix(solid, 0.7, "euclidean")
        self.assertLessEqual(1 - grown.mean(), 0.7 + 1e-3)
        self.assertGreater(1 - grown.mean(), 0.7 - 0.02)
        cross = _grow_matrix(solid, 0.7, "cross")
        # for equal area, a diamond reaches further from the seed than a disc
        yy, xx = np.nonzero(grown)
        extent_grown = np.hypot(yy - 50, xx - 50).max()
        yy, xx = np.nonzero(cross)
        extent_cross = np.hypot(yy - 50, xx - 50).max()
        self.assertLess(extent_grown, extent_cross)
        corner = np.zeros((101, 101), dtype=bool)
        corner[0, 0] = True
        wrapped = _grow_matrix(corner, 0.9, "euclidean")
        self.assertTrue(wrapped[-1, -1] and wrapped[0, -1] and wrapped[-1, 0])

    def test_preset_overrides_are_applied(self):
        with tempfile.TemporaryDirectory() as temporary:
            with contextlib.redirect_stdout(io.StringIO()) as out:
                run_workflow(
                    "periodic2D", Path(temporary), seed=1, verbose=True,
                    overrides={"length": 0.04, "voxel_size": 0.001, "fibre_fraction": 0.4,
                               "notch": None, "z_length": None},
                )
            self.assertIn("Preset overrides", out.getvalue())
            self.assertTrue(any(Path(temporary).glob("iUC1/vp0.04/h0.001/*.vtk")))

    def test_notch_length_and_thickness_and_gradient_are_controllable(self):
        base = replace(
            PRESETS["periodic3D"], notch="matrix", notch_length=0.005, notch_thickness=0.002
        )
        zeros = {}
        for name, config in {
            "default": base,
            "long": replace(base, notch_length=0.01, notch_thickness=0.004),
        }.items():
            material, positions, fractions, surfaces, normals, tangents = _discretize(
                *_toy_layout(), config
            )
            volume = np.repeat(material[:, :, None], _extrusion_depth(config), axis=2)
            volume = _apply_notch(
                volume, positions, fractions, surfaces, normals, tangents, config
            )[0]
            zeros[name] = np.count_nonzero(volume == 0) - np.count_nonzero(material == 0) * volume.shape[2]
        self.assertGreater(zeros["long"], 2 * zeros["default"])
        mini = PRESETS["miniComp"]
        solid = {
            g: np.count_nonzero(
                _discretize(*_toy_layout(length=mini.length), replace(mini, thickness_gradient=g))[0] == 1
            )
            for g in (1.0, 3.0)
        }
        self.assertNotEqual(solid[1.0], solid[3.0])

    def test_mini_composite_domain_keeps_a_border_around_matrix(self):
        config = replace(PRESETS["miniComp"], voxel_size=2e-3, mini_margin=4e-3)
        x, y, radii = _toy_layout(length=config.length)
        pad = _mini_padding(x, y, radii, config)
        shift = pad * config.voxel_size
        domain = replace(config, length=config.length + 2 * shift)
        material = _discretize(
            x + shift, y + shift, radii, domain, region_length=config.length
        )[0]
        solid = np.nonzero(material > 0)
        margin = 2
        self.assertGreaterEqual(min(solid[0].min(), solid[1].min()), margin)
        self.assertGreaterEqual(
            min(material.shape[0] - 1 - solid[0].max(), material.shape[1] - 1 - solid[1].max()),
            margin,
        )

    def test_thickness_ratio_controls_direction_of_gradient(self):
        d = np.array([0.0, 0.5, 1.0])
        np.testing.assert_allclose(_thickness_factors(d, 1.0), [1, 1, 1])
        np.testing.assert_allclose(_thickness_factors(d, 0.5), [1, 0.5**0.5, 0.5])
        np.testing.assert_allclose(_thickness_factors(d, 4.0), [1, 2, 4])
        with self.assertRaises(ValueError):
            _thickness_factors(d, 0.0)


def _toy_layout(length=0.1):
    x = np.array([0.3, 0.5, 0.7]) * length
    y = np.array([0.3, 0.6, 0.4]) * length
    return x, y, np.array([0.006, 0.007, 0.006])


if __name__ == "__main__":
    unittest.main()
