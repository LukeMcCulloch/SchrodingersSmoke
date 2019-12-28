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

sin = np.sin
pi = np.pi
imag = np.complex(0.,1.) #or just use, e.g. '1j' to matlab's '1i'

from ISF import ISF

# example_leapfrog
# An example of incompressible Schroedinger flow producing leapfrogging
# vortex rings.


## PARAMETERS
vol_size = np.asarray([[10.],[5.],[5.]])    # box size
vol_res = np.asarray([[128,64,64]]) # volume resolution
hbar = 0.1            # Planck constant
dt = 1/24             # time step
tmax = 85             # max time
background_vel = [-0.2,0,0] # background velocity

r1 = 1.5              # radius of 1st ring
r2 = 0.9              # radius of 2nd ring
n1 = [-1,0,0]         # normal direction of 1st ring
n2 = [-1,0,0]         # normal direction of 2nd ring

cen1 = vol_size/2 # center of 1st ring
cen2 = vol_size/2 # center of 2nd ring

n_particles = 10000   # number of particles



## INITIALIZATION