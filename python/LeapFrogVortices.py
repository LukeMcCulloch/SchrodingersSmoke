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


# example_leapfrog
# An example of incompressible Schroedinger flow producing leapfrogging
# vortex rings.