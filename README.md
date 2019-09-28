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
We want to compute a FEM approximation of the deformation of a viscoelastic material with isotropic fractional Zener model (See [Section 4, 1](https://github.com/heruslu/Viscoelasticity/blob/master/README.md#references)).
Numerical solver is based on FEM for space and Convolution Quadrature (CQ) \[2\] for time discretization.

## Desiging the MATLAB script

### Space and time discretization
We start with setting the polynomial degree for the FEM solver. 

`` k = 3;``

We choose a final simulation time. It is know that the error of the CQ solver accumulates in time \[3\]. Therefore, for large simulation time, one needs to choose small time stepping.

``
M = 40;
dt = 0.04;
``

### Generating the mesh for space and time
Because of the implementation of `meshRectPrism`, to mark the bottom faces of the prism, we set
 
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
and get the list of Neumann faces.
```
neuList = {find(T.faces(4,:)==2)};
```

Generate the temporal mesh using the given time-step and final time.
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
See Section 1.1 in the documentation for what these parameters correspond in the matematical formulation of the viscoelastic problem. 

### Boundary conditions and forcing term
We define a function for producing vectorized `0`-values.
```
zero = @(x,y,z,t) 0*x;
```
#### Dirichlet condition: Fixes bottom faces of the solid
```
uD{1} = zero; uD{2} = zero; uD{3} = zero;
```
#### Hand press as Neumann condition
Hand press will be activated in the time interval `[tA,tB]`. We define the time interval as a certain proportion of the final time (`M`) of the simulation. The rate of the proportion is controlled by `wP`.
```
wP = 10; tA = M/wP; tB = 2*M/wP;
```
Because of the stability issues, the hand press effect should be smooth in the interval `[tA,tB]`. We achieve this with the `window` function.
```
punchwf = window(tA,tA+2*dt,tB-2*dt,tB);
```
Import the function `hand` which is an indicator function for a human hand in `[0,0]x[1,1]`
```
human_hand_print
```
Set the hand press area and magnitude of this press.
```
punchArea = @(x,y,z) hand(0.5*(x-1),0.5*(y+0.5));
punchMagnitude = -2.8;
punch = @(x,y,z,t) 0*x + punchMagnitude.*punchwf(t).*punchArea(x,y,z);
```
Using `punch` function define the Neumann data.
```
sig{1,1} = zero; sig{1,2} = zero; sig{1,3} = zero;
sig{2,1} = zero; sig{2,2} = zero; sig{2,3} = zero; 
sig{3,1} = zero; sig{3,2} = zero; sig{3,3} = @(x,y,z,t) 0*x+punch(x,y,z,t);
Sig = {sig};
```

#### Gravity as a forcing term
```
gravitywf = window(3*dt,3.1*dt,M+1,M+2);
gravity = 9.8;
f{1} = zero; f{2} = zero; 
f{3} = @(x,y,z,t) 0*x - gravity.*gravitywf(t).*rho(x,y,z);
```

### Solver
We compute `uh, vh, wh` as `x,y,z` components of the displacement, and `sxxh,syyh,szzh,syzh,sxzh,sxyh` as the components of the symmetric stress tensor. These quantities are FEM coefficient matrices for each DOF and each time-step.
```
[uh,vh,wh,sxxh,syyh,szzh,syzh,sxzh,sxyh] = ...
    FEMViscoelasticitySimulation3D({m_mu,b_mu,a_mu,nu_mu},...
                    {m_lam,b_lam,a_lam,nu_lam},rho,f,uD,Sig,T,k,neuList,time,1);
```

### Creating an animation
Our FEM solver for this simualtion is based on `k=3` degree polynomials. To obtain a smooth simulation we interpolate our FEM solution at high order quadrature points. We only compute the deformation at the surface of the solid, and create a surface mesh at these quadrature points.
```
kinterp = 9;
uhtotSurf = interpolateFEM3D([uh;vh;wh],T,k,kinterp);
Tsurf = interpolateMesh(T,kinterp);
```
We also create a hand mesh and an approproate movement for it.
```
handmovewf = window(0.5*tA,tA+13*dt,tB-2*dt,tB);
hand_top_lev = 4; hand_bottom_lev = 1.9;
surfLevel = min(wh);
zmove = -0.2 + hand_top_lev - (hand_top_lev - hand_bottom_lev - 0.55*surfLevel).*handmovewf(time);
load('hand_mesh.mat');
Nvhand = length(Thand.Z);
whhand = zeros(Nvhand,nt) + zmove;
uhtotHand = [zeros(Nvhand,nt);zeros(Nvhand,nt)-0.9;whhand];
```
Finally we construct the animation with a choice of a pair of colors for the solid and the hand.
```
movie_fn =@() axis('off');
color = {[151,190,162]/256,[175, 110, 81]/256};
animateMultSurfFEM3D({uhtotSurf,uhtotHand},{Tsurf,Thand},...
    color,0,movie_fn);
```

### References
\[1\]: T.S. Brown, S. Du, H. Eruslu, and F.-J. Sayas. Analysis of models for viscoelastic wave propagation. Applied Mathematics and Nonlinear Sciences, 2018. ([DOI:10.21042/AMNS.2018.1.00006](https://doi.org/10.1093/imanum/drx079))

\[2\] M. Hassell, and F.-J. Sayas. Convolution Quadrature for Wave Simulations. Numerical simulation in physics and engineering. SEMA SIMAI Springer Ser.(9) pp 71-159, 2016. ([Chapter](https://link-springer-com.udel.idm.oclc.org/chapter/10.1007/978-3-319-32146-2_2))

\[3\]: H. Eruslu, and F.-J. Sayas. Brushing up a theorem by Lehel Banjai on the convergence of Trapezoidal Rule Convolution Quadrature, 2019 ([arxiv:1903.09031](https://arxiv.org/abs/1903.09031))
