# -*- coding:utf-8 -*-
# Loadings and Output definitions

from py4amitex.input import (
    Loading,
    MechanicDriving,
    Component,
    LoadingOutput,
    DiffusionDriving,
    Evolution,
    Output,
    DirStress,
)


def traction_1():
    loadingOutput = LoadingOutput()
    loading = Loading()
    loading.setTimeDiscretizationLinear(1, 100.0)
    loading.setLinearEvolution(Component.XX, MechanicDriving.Strain, 0.01)
    for i in range(3):
        for j in range(i, 3):
            if i != 0 or j != 0:
                loading.setLinearEvolution((i, j), MechanicDriving.Stress, 0.0)
    loadingOutput.add(loading)
    return loadingOutput


def diffusion_flux():
    loadingOutput = LoadingOutput()

    output = Output()
    output.setVtkFluxDGradD(1, 1)
    loadingOutput.output = output

    loading = Loading()
    loading.setTimeDiscretizationUser([1.0])
    loading.setOutputVtkList([1])
    for i in range(3):
        loading.setEvolution(i, DiffusionDriving.Flux, Evolution.Linear, i + 1)
    loadingOutput.add(loading)

    return loadingOutput


def def_imp():
    loadingOutput = LoadingOutput()

    output = Output()
    output.setVtkStressStrain(1, 1)
    loadingOutput.output = output

    loading = Loading()
    loading.setTimeDiscretizationUser([1.0])
    loading.setOutputVtkList([1])
    loading.setEvolution(Component.XX, MechanicDriving.Strain, Evolution.Linear, -5e-4)
    loading.setEvolution(Component.YY, MechanicDriving.Strain, Evolution.Linear, 1.5e-4)
    loading.setEvolution(Component.ZZ, MechanicDriving.Strain, Evolution.Linear, 1.5e-4)
    loading.setEvolution(Component.XY, MechanicDriving.Strain, Evolution.Linear, 0.0)
    loading.setEvolution(Component.XZ, MechanicDriving.Strain, Evolution.Linear, 0.0)
    loading.setEvolution(Component.YZ, MechanicDriving.Strain, Evolution.Linear, 0.0)
    loadingOutput.add(loading)

    return loadingOutput


def triax2p5(cauchy=False):
    loading = Loading()
    loading.setTimeDiscretizationLinear(40, 40.0)
    loading.setOutputVtkList([1])
    loading.setEvolution(
        Component.XX, MechanicDriving.Stress, Evolution.Constant, dirstress=1.0
    )
    loading.setEvolution(
        Component.YY, MechanicDriving.Stress, Evolution.Constant, dirstress=1.3
    )
    loading.setEvolution(
        Component.ZZ,
        MechanicDriving.Strain,
        Evolution.Linear,
        value=0.004,
        dirstress=1.6,
    )
    for i in range(3):
        for j in range(3):
            if i != j:
                loading.setEvolution(
                    (i, j), MechanicDriving.Stress, Evolution.Constant, dirstress=0
                )
    loading.setDirStress(DirStress.Cauchy if cauchy else DirStress.PK1)
    output = Output()
    output.setVtkStressStrain(0, 0)

    loadingOutput = LoadingOutput()
    loadingOutput.output = output
    loadingOutput.add(loading)
    return loadingOutput


def defimp_torsion_beam():
    loading = Loading()
    loading.setTimeDiscretizationLinear(1, 100.0)
    loading.setOutputVtkList([1])

    loading.setLinearEvolution(Component.XX, MechanicDriving.Strain, 0.01)
    for i in range(3):
        for j in range(3):
            if (i, j) != Component.XX:
                loading.setLinearEvolution((i, j), MechanicDriving.Strain, 0.0)
    #  gives epsxz,y =  1
    loading.setGradGradU(Component.ZX, Component.Y, Evolution.Linear, 0.002)
    # gives epsxy,z = -1
    loading.setGradGradU(Component.YX, Component.Z, Evolution.Linear, -0.002)

    output = Output()
    output.setVtkStressStrain(1, 1)

    loadingOutput = LoadingOutput()
    loadingOutput.add(loading)
    loadingOutput.output = output
    return loadingOutput
