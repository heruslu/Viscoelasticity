# Viscoelasticity
This repository explains how to build a simulation for viscoelastic material deformation using `team-pancho/FEM3D` MATLAB package. Note that this package has not yet released, but will be available under [team-pancho GitHub page](https://github.com/team-pancho) soon.

Here we will demonstrate how to build a simulation using this package. This example requires the following functions from `team-pancho/FEM3D`:
* `meshRectPrism`
* `window`
* `FEMViscoelastictySimulation3D`
* `interpolateFEM3D`
* `interpolateMesh`
* `animateMultSurfFEM3D`

Numerical solver is based on FEM for space and Convolution Quadrature (CQ).[^1] for time discretization. Viscoelastic model is a fractional isotropic Zener model. More information

# Time and space discretization
We start with setting the polinomial degree for the FEM solver. 

`` k = 3;``

We choose a final simulation time. It


[^1]: Footnote
