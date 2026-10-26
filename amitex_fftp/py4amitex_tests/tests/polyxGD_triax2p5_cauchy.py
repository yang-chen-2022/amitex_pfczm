# -*- coding:utf-8 -*-
from py4amitex.input import (
    Input,
    Grid,
)

import algo
import char
import mate


def make_input():
    algoParams = algo.Smacro_GD()
    loadingOutput = char.triax2p5(cauchy=True)
    (materials, grid) = mate.polyx_GD()

    input = Input(grid, algoParams, materials, loadingOutput)
    return input
