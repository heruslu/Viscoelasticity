% Script to compute the deformation of a viscoelastic material due to a 
% hand press

% Simulation design: Hasan Eruslu, University of Delaware
% Package: team-pancho/FEM3D
% Last modified: September 27, 2019

%% Time and space discretization
tic;
k = 3; % Polynomial degree for FEM approximation
M = 40; % Final simulation time
dt = 0.04; % Time steps 

%% Generating space and time mesh
% Rectangle Prism Mesh with approproate Neumann faces
ref = 2; % Refinement of the mesh
DBC = 1; % Bottom face is Dirichlet
% Create the initial prism as [0,8] x [0,8] x [0,2]
ztop = 2; ytop = 8; xtop = 8;
T = meshRectPrism(xtop,ytop,ztop,ref,DBC);
% Scale the prism in x and y directions by 0.5
scale = 0.5; 
T.coordinates(1:2,:) = scale*T.coordinates(1:2,:);
% Create enhanced properties of the mesh T
T = edgesAndFaces(T); % Information of edge and face indexing
T = enhanceGrid3D(T); % Face orientations, volume and area of elements
neuList = {find(T.faces(4,:)==2)}; % List of Neumann faces

% Temporal mesh generating
time = 0:dt:M; nt = length(time);

%% Viscoelastic parameters
% mass density
rho_c = 0.01;
rho =@(x,y,z) rho_c+0*x;
% fractional zener parameters
m_muc = 0.2; a_muc = 0.6; b_muc = 2.0; nu_muc = 0.7; 
m_mu = @(x,y,z) m_muc+0*x; a_mu = @(x,y,z) a_muc+0*x; b_mu = @(x,y,z) b_muc+0*x;
nu_mu = @(x,y,z) nu_muc+0*x;
m_lamc = 0.9; a_lamc = 0.3; b_lamc = 3.0; nu_lamc = 0.7;
m_lam = @(x,y,z) m_lamc+0*x; a_lam = @(x,y,z) a_lamc+0*x; b_lam = @(x,y,z) b_lamc+0*x;
nu_lam = @(x,y,z) nu_lamc+0*x;

%% Boundary conditions and forcing term
zero = @(x,y,z,t) 0*x;
% Dirichlet condition: Fixes bottom faces of the solid
uD{1} = zero; uD{2} = zero; uD{3} = zero;

% Hand press as Neumann condition
wP = 10; % Portion of the total time for hand press activation
% Hand press will be activated in the time interval [tA,tB]
tA = M/wP; tB = 2*M/wP;
punchwf = window(tA,tA+2*dt,tB-2*dt,tB); % Window function for hand press

human_hand_print % imports the function hand
punchArea = @(x,y,z) hand(0.5*(x-1),0.5*(y+0.5));
punchMagnitude = -2.8;
punch = @(x,y,z,t) 0*x + punchMagnitude.*punchwf(t).*punchArea(x,y,z);

sig{1,1} = zero; sig{1,2} = zero; sig{1,3} = zero;
sig{2,1} = zero; sig{2,2} = zero; sig{2,3} = zero; 
sig{3,1} = zero; sig{3,2} = zero; sig{3,3} = @(x,y,z,t) 0*x+punch(x,y,z,t);
Sig = {sig};

% Gravity as a forcing term
gravitywf = window(3*dt,3.1*dt,M+1,M+2);
gravity = 9.8;
f{1} = zero; f{2} = zero; 
f{3} = @(x,y,z,t) 0*x - gravity.*gravitywf(t).*rho(x,y,z);

%% Solver
% uh, vh, wh: DOF for the displacement
% sxxh,... : DOF for the stress tensor components
[uh,vh,wh,sxxh,syyh,szzh,syzh,sxzh,sxyh] = ...
    FEMViscoelasticitySimulation3D({m_mu,b_mu,a_mu,nu_mu},...
                    {m_lam,b_lam,a_lam,nu_lam},rho,f,uD,Sig,T,k,neuList,time,1);
%% Creating the animation
save_snapshots;
toc;