#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 00:30:07 2019

@author: lukemcculloch
"""


from __future__ import division

import numpy as np
from numpy import exp
fft = np.fft.fft
ifft = np.fft.ifft
fftn = np.fft.fftn
ifftn = np.fft.ifftn
fftshift = np.fft.fftshift

rand = np.random.rand
ones = np.ones
size = np.shape
ceil = np.ceil

sin = np.sin
pi = np.pi
imag = np.complex(0.,1.) #or just use, e.g. '1j' to matlab's '1i'

import matplotlib.pyplot as plt
from plotTrajectory import *

from ISF import ISF
from Particles import Particles



    
def plot3D(particle, canvas = None, color='black'):
    """
    """
    if canvas is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    else:
        ax = canvas
        
    ax.plot(particle.x,particle.y,particle.z)
    return ax


# example_leapfrog
# An example of incompressible Schroedinger flow producing leapfrogging
# vortex rings.


## PARAMETERS
vol_size = np.asarray([10.,5.,5.])   # box size
vol_res = [128,64,64] # volume resolution
hbar = 0.1            # Planck constant
dt = 1/24             # time step
tmax = 1             # max time
background_vel = np.asarray([-0.2,0,0]) # background velocity

r1 = 1.5              # radius of 1st ring
r2 = 0.9              # radius of 2nd ring
n1 = [-1,0,0]         # normal direction of 1st ring
n2 = [-1,0,0]         # normal direction of 2nd ring

cen1 = vol_size/2 # center of 1st ring
cen2 = vol_size/2 # center of 2nd ring

n_particles = 10000   # number of particles



## INITIALIZATION


isf = ISF(sizex = vol_size[0],
          sizey = vol_size[1],
          sizez = vol_size[2],
          resx  = vol_res[0],
          resy  = vol_res[1],
          resz  = vol_res[2],
          hbar = hbar)
self = isf

isf.hbar = hbar
isf.dt = dt
#"""
isf.BuildSchroedinger()
# Set background velocity
kvec = background_vel/isf.hbar
phase = kvec[0]*isf.px + kvec[1]*isf.py + kvec[2]*isf.pz
#"""
"""
python:
    isf.px[0:10,0:10,0].T
    phase[0:10,0:10,0].T
    
    psi1[0:10,0:10,0].T
    psi2[0:10,0:10,0].T
    
matlab:
    isf.px(1:10,1:10)
    phase(1:10,1:10)
    
    psi1(1:10,1:10)
    psi2(1:10,1:10)
#"""
#"""
psi1 = exp(1j*phase)
psi2 = 0.01*exp(1j*phase)
# Add vortex rings
d = isf.dx*5 # neighborhood thickness of seifert surface
#"""
"""
self = isf
psi=psi1
center=cen1
normal=n1
r=r1
d=d
#"""
#"""
psi1 = isf.AddCircle(psi1,cen1,n1,r1,d)
psi1 = isf.AddCircle(psi1,cen2,n2,r2,d)
[psi1,psi2] = isf.Normalize(psi1,psi2)
#"""
#"""
[psi1,psi2] = isf.PressureProject(psi1,psi2)

## SET PARTICLES
uu = rand(n_particles,1)
vv = rand(n_particles,1)
party = 0.5 + 4*uu
partz = 0.5 + 4*vv
partx = 5*ones(size(party))
particle = Particles(x=partx,
                     y=party,
                     z=partz)

#"""
"""
axis = plot3D(particle.x,particle.y,particle.z,'.','MarkerSize',1)
axis equal
axis([0,vol_size{1},0,vol_size{2},0,vol_size{3}])
cameratoolbar
drawnow
#"""

#"""
## MAIN ITERATION
itermax = int(ceil(tmax/dt))
#for iter = 1:itermax
iter = 0
while(iter < itermax):
    iter +=1
    t = iter*dt
    print t
    # incompressible Schroedinger flow
    [psi1,psi2] = isf.SchroedingerFlow(psi1,psi2)
    [psi1,psi2] = isf.Normalize(psi1,psi2)
    [psi1,psi2] = isf.PressureProject(psi1,psi2)
    
    # particle visualization
#    [vx,vy,vz] = isf.VelocityOneForm(psi1,psi2,isf.hbar)
#    [vx,vy,vz] = isf.StaggeredSharp(vx,vy,vz)
#    particle.StaggeredAdvect(isf,vx,vy,vz,isf.dt)
#    particle.Keep(particle.x>0&particle.x<vol_size{1}&...
#                  particle.y>0&particle.y<vol_size{2}&...
#                  particle.z>0&particle.z<vol_size{3})
#    set(hpart,'XData',particle.x,'YData',particle.y,'ZData',particle.z)
#    title(['iter = ',num2str(iter)])
#    drawnow

#"""