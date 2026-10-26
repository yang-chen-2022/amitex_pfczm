# -*- coding:utf-8 -*-
from py4amitex.input import (
    Input,
    Grid,
    Materials,
    ReferenceMaterial
)
from py4amitex.input.materialbuilder import buildMaterialsFromVtk

import algo
import char
import mate

def make_input():
    algoParams = algo.default()
    loadingOutput = char.traction_1()

    materials = Materials()
    materials.referenceMaterial = ReferenceMaterial(5.76923076923077e8, 3.84615384615385e8)
    grid = buildMaterialsFromVtk(materials, zoneIdPath="../cas_tests/microstructures/voxcomp/polyX_27G_R21.vtk")
    
    material = materials.material(0)
    material.setLaw("elasiso")
    material.setNumberCoeff(2)
    material.setNumberCoeffComposite(2)

    coeffdir = "../cas_tests/mate/coefficients"
    material.setCoeffZoneFromBin(0, f"{coeffdir}/Lambda1_polyx27G_R21.bin")
    material.setCoeffZoneFromBin(1, f"{coeffdir}/Mu1_polyx27G_R21.bin")
    material.setCoeffCompositeZoneFromBin(0, f"{coeffdir}/Lambda1_polyx27G_R21.bin")
    material.setCoeffCompositeZoneFromBin(1, f"{coeffdir}/Mu1_polyx27G_R21.bin")

    input = Input(grid, algoParams, materials, loadingOutput)
    return input
