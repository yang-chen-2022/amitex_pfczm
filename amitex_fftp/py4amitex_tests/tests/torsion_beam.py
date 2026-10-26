# -*- coding:utf-8 -*-
from py4amitex.input import (
    Input,
    Grid,
)

import algo
import char
import mate


def make_input():
    algoParams = algo.default()
    loadingOutput = char.defimp_torsion_beam()
    (materials, grid)  = mate.renfort2_i_m0()

    input = Input(grid, algoParams, materials, loadingOutput)
    return input
