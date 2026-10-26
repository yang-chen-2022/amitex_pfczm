# -*- coding:utf-8 -*-
from py4amitex.input import (
    Input,
    Grid,
)

import algo
import char
import mate


def make_input():
    grid = Grid([15, 15, 15], [1.0, 1.0, 1.0])

    algoParams = algo.default_diffusion()
    materials = mate.lin_diffusion()
    loadingOutput = char.diffusion_flux()

    return Input(grid, algoParams, materials, loadingOutput)
