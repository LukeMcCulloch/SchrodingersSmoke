% example_jet
% An example of incompressible Schroedinger flow producing a jet.
%
clear
clc
%% PARAMETERS
vol_size = {4,2,2};   % box size
vol_res = {64,32,32}; % volume resolution
hbar = 0.1;            % Planck constant
dt = 1/48;             % time step
tmax = 50;             % max time

% another interesting set of parameter:
%vol_size = {4,2,2};   % box size
%vol_res = {128,64,64}; % volume resolution
%hbar = 0.02;            % Planck constant
%dt = 1/48;             % time step
%tmax = 10;             % max time

jet_velocity = [1,0,0]; % jet velocity

nozzle_cen = [2-1.7, 1-0.034, 1+0.066]; % nozzle center
nozzle_len = 0.5;                   % nozzle length
nozzle_rad = 0.5;                   % nozzle radius

n_particles = 50;   % number of particles



%% INITIALIZATION


isf = ISF(vol_size{:},vol_res{:});
isf.hbar = hbar;
isf.dt = dt;
isf.BuildSchroedinger;

% Set nozzle
isJet = abs(isf.px - nozzle_cen(1))<=nozzle_len/2 & ...
    (isf.py - nozzle_cen(2)).^2+(isf.pz - nozzle_cen(3)).^2<=nozzle_rad.^2;

% initialize psi
psi1 = ones(size(isf.px));
psi2 = 0.01*psi1;
[psi1,psi2] = isf.Normalize(psi1,psi2);

% constrain velocity
kvec = jet_velocity/isf.hbar;
omega = sum(jet_velocity.^2)/(2*isf.hbar);
phase = kvec(1).*isf.px + kvec(2).*isf.py + kvec(3).*isf.pz;
for iter = 1:10
    amp1 = abs(psi1);
    amp2 = abs(psi2);
    psi1(isJet) = amp1(isJet).*exp(1i*phase(isJet));
    psi2(isJet) = amp2(isJet).*exp(1i*phase(isJet));
    [psi1,psi2] = isf.PressureProject(psi1,psi2);
end


%% SET PARTICLES
particle = Particles;
particle.x = 0;
particle.y = 0;
particle.z = 0;

hpart = plot3(particle.x,particle.y,particle.z,'.','MarkerSize',1);
axis equal;
axis([0,vol_size{1},0,vol_size{2},0,vol_size{3}]);
cameratoolbar
drawnow

%% MAIN ITERATION
itermax = ceil(tmax/dt);
for iter = 1:itermax
    t = iter*dt;
    
    % incompressible Schroedinger flow
    [psi1,psi2] = isf.SchroedingerFlow(psi1,psi2);
    [psi1,psi2] = isf.Normalize(psi1,psi2);
    [psi1,psi2] = isf.PressureProject(psi1,psi2);
    
    % constrain velocity
    phase = kvec(1).*isf.px + kvec(2).*isf.py + kvec(3).*isf.pz - omega*t;
    amp1 = abs(psi1);
    amp2 = abs(psi2);
    psi1(isJet) = amp1(isJet).*exp(1i*phase(isJet));
    psi2(isJet) = amp2(isJet).*exp(1i*phase(isJet));
    [psi1,psi2] = isf.PressureProject(psi1,psi2);
    
    % particle birth
    rt = rand(n_particles,1)*2*pi;
    newx = nozzle_cen(1)*ones(size(rt));
    newy = nozzle_cen(2) + 0.9*nozzle_rad*cos(rt);
    newz = nozzle_cen(3) + 0.9*nozzle_rad*sin(rt);
    particle.x = [particle.x;newx];
    particle.y = [particle.y;newy];
    particle.z = [particle.z;newz];
    
    % advect and show particles
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