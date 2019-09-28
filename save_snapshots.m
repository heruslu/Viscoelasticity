% Script to  produce the simulation of a hand press 
% Requires the data created by viscoelastic_deformation_hand_press.m

% Create the hand move
handmovewf = window(0.5*tA,tA+13*dt,tB-2*dt,tB);
hand_top_lev = 4; hand_bottom_lev = 1.9;
surfLevel = min(wh);
zmove = -0.2 + hand_top_lev - (hand_top_lev - hand_bottom_lev - 0.55*surfLevel).*handmovewf(time);
load('hand_mesh.mat'); % loads Thand, mesh for the hand
Nvhand = length(Thand.Z);
whhand = zeros(Nvhand,nt) + zmove;
% DOF for the FEM vector of the hand move
uhtotHand = [zeros(Nvhand,nt);zeros(Nvhand,nt)-0.9;whhand];

% Interpolate solution and mesh
kinterp = 9;
uhtotSurf = interpolateFEM3D([uh;vh;wh],T,k,kinterp);
Tsurf = interpolateMesh(T,kinterp);

% Create the animation
movie_fn =@() axis('off');
% Color for the material and the hand
color = {[151,190,162]/256,[175, 110, 81]/256};
animateMultSurfFEM3D({uhtotSurf,uhtotHand},{Tsurf,Thand},...
    color,0,movie_fn);