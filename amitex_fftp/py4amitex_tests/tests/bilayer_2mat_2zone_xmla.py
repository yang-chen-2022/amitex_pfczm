# -*- coding:utf-8 -*-
from py4amitex.input import (
    Input,
    Grid,
)

import algo
import char
import mate


def make_input():
    N = 40
    DX = 0.025
    grid = Grid([N, N, N], [DX, DX, DX])
    algoParams = algo.default()
    loadingOutput = char.traction_1()
    materials = mate.m2mat_2zone_xml(grid)
    input = Input(grid, algoParams, materials, loadingOutput)
    return input
