# Viscoelasticity
This repository explains how to build a simulation for viscoelastic material deformation using `team-pancho/FEM3D` MATLAB package. Note that this package has not yet released, but will be available under [team-pancho GitHub page](https://github.com/team-pancho) soon.

Here we will demonstrate how to build a simulation using this package. This example requires the following functions from `team-pancho/FEM3D`:
* `meshRectPrism`
* `window`
* `FEMViscoelastictySimulation3D`
* `interpolateFEM3D`
* `interpolateMesh`
* `animateMultSurfFEM3D`

## Goal
We want to compute a FEM approximation of the deformation of a viscoelastic material with isotropic fractional Zener model (See \[Section 4, 1\]).
Numerical solver is based on FEM for space and Convolution Quadrature (CQ) \[2\] for time discretization.

## Desiging the MATLAB script

### Space and time discretization
We start with setting the polinomial degree for the FEM solver. 

`` k = 3;``

We choose a final simulation time. It is know that CQ solver error accumulates in time \[3\], so for large simulation time one needs to choose small time stepping.

``
M = 40;
dt = 0.04;
``

### Generating the mesh for space and time
Because of the implementation of `meshRectPrism`, to mark the bottom faces of the prism we set
 
``DBC = 1``
 
We want to create the initial prism `[0,8] x [0,8] x [0,2]` with 2 refinements,
 
```
ref = 2;
ztop = 2; ytop = 8; xtop = 8;
T = meshRectPrism(xtop,ytop,ztop,ref,DBC);
```
and then scale the prism in `x` and `y` directions by `0.5`.
```
scale = 0.5; 
T.coordinates(1:2,:) = scale*T.coordinates(1:2,:);
```
Create enhanced properties of the mesh `T`,
```
T = edgesAndFaces(T); T = enhanceGrid3D(T);
```
and get the list of Neumann faces
```
neuList = {find(T.faces(4,:)==2)};
```

Generate the temporal mesh using the given time-step and final time
```
time = 0:dt:M; nt = length(time);
```

### Viscoelastic parameters
Define the mass density function,
```
rho_c = 0.01;
rho =@(x,y,z) rho_c+0*x;
```
and fractional zener parameters.
```
m_muc = 0.2; a_muc = 0.6; b_muc = 2.0; nu_muc = 0.7; 
m_mu = @(x,y,z) m_muc+0*x; a_mu = @(x,y,z) a_muc+0*x; b_mu = @(x,y,z) b_muc+0*x;
nu_mu = @(x,y,z) nu_muc+0*x;
m_lamc = 0.9; a_lamc = 0.3; b_lamc = 3.0; nu_lamc = 0.7;
m_lam = @(x,y,z) m_lamc+0*x; a_lam = @(x,y,z) a_lamc+0*x; b_lam = @(x,y,z) b_lamc+0*x;
nu_lam = @(x,y,z) nu_lamc+0*x;
```
See the documentation for what these parameters correspond in the matematical formulation of the viscoelastic problem. 




### References
\[1\]: Viscoelasticity paper
\[2\]Sayas, Hassell CQ
\[2\]: CQTR
