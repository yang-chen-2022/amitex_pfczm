# -*- coding:utf-8 -*-
# Materials (and structure) definitions

from itertools import product
import os

from py4amitex.input import (
    Material,
    Materials,
    ReferenceMaterial,
    ReferenceMaterialD,
    Zone,
)
from py4amitex.input.materialbuilder import buildMaterialsFromVtk


def __makeRectangularZone(begin, end, dims):
    return Zone(
        dims,
        [
            (ix, iy, iz)
            for ix in range(begin[0], end[0])
            for iy in range(begin[1], end[1])
            for iz in range(begin[2], end[2])
        ],
    )


def m2mat_2zone_xml(grid):
    N = grid.dims()[0]
    LAMBDA = [[100.0e6, 200.0e6], [300.0e6, 400.0e6]]
    MU = [[100.0e6, 200.0e6], [300.0e6, 400.0e6]]

    materials = Materials()
    materials.referenceMaterial = ReferenceMaterial(250e6, 250e6)

    for m in range(2):
        mat = Material()
        mat.setLaw("elasiso", "")
        mat.setNumberCoeff(2)

        minY = N // 2 if m == 0 else 0
        maxY = N if m == 0 else N // 2

        for z in range(2):
            minZ = N // 2 if z == 0 else 0
            maxZ = N if z == 0 else N // 2
            zone = __makeRectangularZone((0, minY, minZ), (N, maxY, maxZ), (N, N, N))
            mat.addZone(zone, [LAMBDA[m][z], MU[m][z]])
        materials.add(mat)

    return materials


def __arlequin(material, nzones, nvox, coeffK=False):
    for izone in product(range(nzones), repeat=3):
        idx = izone[0] + nzones * (izone[1] + nzones * izone[2])
        zone = Zone([nzones * nvox for _ in range(3)])
        for ip in product(range(nvox), repeat=3):
            zone.add([izone[k] * nvox + ip[k] for k in range(3)])
        if coeffK:
            material.addZone(zone, coeffKs=[float(idx + 1)])
        else:
            material.addZone(zone)


def lin_diffusion():
    materials = Materials()
    materials.referenceMaterialD = ReferenceMaterialD(14.0)

    material = Material()
    material.setLawK("fourier_iso")
    material.setNumberCoeffK(1)
    __arlequin(material, 3, 5, True)
    materials.add(material)

    return materials


def __getCoeffFromTxt(path):
    with open(path) as file:
        return [float(c) for c in file.read().split()]


def polyx_GD():
    amitexPath = os.getenv("AMITEX_PATH", "..")

    materials = Materials()
    grid = buildMaterialsFromVtk(
        materials,
        zoneIdPath=f"{amitexPath}/cas_tests/microstructures/arlequin/arlequin_N2_r16.vtk",
    )

    materials.referenceMaterial = ReferenceMaterial(2.0431e5, 0.8756e5)

    material = materials.material(0)
    material.setLaw(
        "umatBCCGD",
        f"{amitexPath}/cas_tests/comportements/polyxGD/comportement_umat/libUmatAmitex.so",
    )
    material.setNumberCoeff(28)
    coeff1to25 = [
        236.412e3,
        0.35,
        275.2e3,
        112.4e3,
        87.56e3,
        363.0,
        410.0,
        480.0,
        0.0,
        1000.0,
        10e-6,
        1e11,
        2.481e-7,
        100,
        3.0,
        2.0e-6,
        0.84,
        5e4,
        5e4,
        1e-6,
        50,
        0,
        0.04,
        0.0,
        50,
    ]
    for idx, coeff in enumerate(coeff1to25):
        material.setCoeff(idx, coeff)
    for idx in range(25, 28):
        path = f"{amitexPath}/cas_tests/mate/coefficients/Coeff{idx}"
        material.setCoeffZone(idx, __getCoeffFromTxt(path))

    for idx in range(107):
        material.addIntVar(0.0)

    return (materials, grid)


def renfort2_i_m0():
    amitexPath = os.getenv("AMITEX_PATH", "..")

    materials = Materials()
    grid = buildMaterialsFromVtk(
        materials,
        materialIdPath=f"{amitexPath}/cas_tests/microstructures/beam/cube44_cylinder_mate.vtk",
        zoneIdPath=f"{amitexPath}/cas_tests/microstructures/beam/cube44_cylinder_zone.vtk",
        minId=0,
    )

    materials.referenceMaterial = ReferenceMaterial(
        2.884615384615384, 1.923076923076923
    )

    # MATERIAL 1 : surrounding void material
    material1 = materials.material(0)
    material1.setLaw("elasiso")
    material1.setCoeffs([0e3, 0e3])

    # MATERIAL 2 : matrix and inclusion
    material2 = materials.material(1)
    material2.setLaw("elasiso")
    material2.setNumberCoeff(2)
    #                          matrix              inclusion
    material2.setCoeffZone(0, [0.576923076923077e1, 0.0])
    material2.setCoeffZone(1, [0.384615384615385e1, 0.0])

    return (materials, grid)
