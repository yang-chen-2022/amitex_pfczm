# -*- coding:utf-8 -*-
# Algorithm parameters definitions

from py4amitex.input import (
    Diffusion,
    Algorithm,
    AlgorithmParameters,
    Mechanics,
    Algorithm,
    AlgorithmParameters,
)


def default():
    mechanics = Mechanics(filter="Default", smallPerturbations=True)
    return AlgorithmParameters(
        Algorithm(type="Basic_Scheme", convergenceAcceleration=True),
        mechanics=mechanics,
    )


def default_diffusion():
    algorithm = Algorithm(type="Basic_Scheme", convergenceAcceleration=True)
    diffusion = Diffusion(filter="Default", stationary=True)
    return AlgorithmParameters(algorithm, diffusion=diffusion)


def Smacro_GD():
    algorithm = Algorithm(
        type="Basic_Scheme",
        convergenceCriterionSmacro=1.0e-3,
        convergenceAcceleration=True,
    )

    mechanics = Mechanics(filter="Default", smallPerturbations=False)

    return AlgorithmParameters(algorithm, mechanics=mechanics)
