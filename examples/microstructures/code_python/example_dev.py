"""Python counterpart of ``code_matlab/SiCSiCpore_dev.m`` test data."""

import argparse
from pathlib import Path

import numpy as np

from .microstructure import saveVolvtk_amitex, write_composite_voxel_data


def main(output: Path = Path("python_generated/SiCSiCpore_dev")) -> None:
    output.mkdir(parents=True, exist_ok=True)
    volume = np.asarray([1, 2, 2, 2, 1, 0], dtype=np.uint8).reshape((6, 1, 1))
    positions = np.asarray([2, 4], dtype=np.uint64)
    fractions = np.asarray([[0.5, 0.4], [0.5, 0.4]])
    normals = np.tile([1.0, 0.0, 0.0], (2, 1))
    tangents = np.tile([0.0, 0.0, 1.0], (2, 1))
    write_composite_voxel_data(
        output / "rep_3_4_2", positions, fractions, np.ones(2), normals, tangents
    )
    saveVolvtk_amitex(volume, output / "SiCSiCpore_dev.vtk", "uint8", "cell", [1, 1, 1])
    print(f"Example geometry written to {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=Path("python_generated/SiCSiCpore_dev")
    )
    main(parser.parse_args().output)
