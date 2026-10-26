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
    loadingOutput = char.def_imp()

    materials = Materials()    
    materials.referenceMaterial = ReferenceMaterial(22.21175e+10, 14.80765e+10)
    grid = buildMaterialsFromVtk(materials, zoneIdPath="../cas_tests/microstructures/al2p5_65/mat_al2p5_65.vtk")
    
    material = materials.material(0)
    material.setLaw("elasiso")
    material.setNumberCoeff(2)
    material.setCoeffZone(0, [4.0385e+11, 4.0385e+10])
    material.setCoeffZone(1, [2.6923e+11, 2.6923e+10])

    input = Input(grid, algoParams, materials, loadingOutput)
    return input
