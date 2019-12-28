% example_leapfrog
% An example of incompressible Schroedinger flow producing leapfrogging
% vortex rings.
%
clear
clc
%% PARAMETERS
vol_size = {10,5,5};   % box size
vol_res = {128,64,64}; % volume resolution
hbar = 0.1;            % Planck constant
dt = 1/24;             % time step
tmax = 85;             % max time
background_vel = [-0.2,0,0]; % background velocity

r1 = 1.5;              % radius of 1st ring
r2 = 0.9;              % radius of 2nd ring
n1 = [-1,0,0];         % normal direction of 1st ring
n2 = [-1,0,0];         % normal direction of 2nd ring

cen1 = [vol_size{:}]/2; % center of 1st ring
cen2 = [vol_size{:}]/2; % center of 2nd ring

n_particles = 10000;   % number of particles



%% INITIALIZATION


isf = ISF(vol_size{:},vol_res{:});
isf.hbar = hbar;
isf.dt = dt;
isf.BuildSchroedinger;
% Set background velocity
kvec = background_vel/isf.hbar;
phase = kvec(1).*isf.px + kvec(2).*isf.py + kvec(3).*isf.pz;
psi1 = exp(1i*phase);
psi2 = 0.01*exp(1i*phase);
% Add vortex rings
d = isf.dx*5; % neighborhood thickness of seifert surface
psi1 = isf.AddCircle(psi1,cen1,n1,r1,d);
psi1 = isf.AddCircle(psi1,cen2,n2,r2,d);
[psi1,psi2] = isf.Normalize(psi1,psi2);
[psi1,psi2] = isf.PressureProject(psi1,psi2);

%% SET PARTICLES
uu = rand(n_particles,1);
vv = rand(n_particles,1);
party = 0.5 + 4*uu;
partz = 0.5 + 4*vv;
partx = 5*ones(size(party));
particle = Particles;
particle.x = partx;
particle.y = party;
particle.z = partz;

hpart = plot3(particle.x,particle.y,particle.z,'.','MarkerSize',1);
axis equal;
axis([0,vol_size{1},0,vol_size{2},0,vol_size{3}]);
cameratoolbar
%drawnow

%% MAIN ITERATION
itermax = ceil(tmax/dt);
for iter = 1:itermax
    t = iter*dt;
    % incompressible Schroedinger flow
    [psi1,psi2] = isf.SchroedingerFlow(psi1,psi2);
    [psi1,psi2] = isf.Normalize(psi1,psi2);
    [psi1,psi2] = isf.PressureProject(psi1,psi2);
    
    % particle visualization
    [vx,vy,vz] = isf.VelocityOneForm(psi1,psi2,isf.hbar);
    [vx,vy,vz] = isf.StaggeredSharp(vx,vy,vz);
    particle.StaggeredAdvect(isf,vx,vy,vz,isf.dt);
    particle.Keep(particle.x>0&particle.x<vol_size{1}&...
                  particle.y>0&particle.y<vol_size{2}&...
                  particle.z>0&particle.z<vol_size{3})
    set(hpart,'XData',particle.x,'YData',particle.y,'ZData',particle.z);
    title(['iter = ',num2str(iter)])
    drawnow
end