

import numpy as np

from TorusDEC import TorusDEC


class ISF(TorusDEC):
    """
    ISF a handle class for simulating incompressible Schroedinger flow.
        SYNTAX
          isf = ISF(sizex,sizey,sizez,resx,resy,resz)
      isf = ISF(sizex,sizey,sizez,res)
        DESCRIPTION
          ISF is a subclass of TorusDEC handle class. See TorusDEC for calling
      constructor.
          To setup, 
            isf = ISF(5,5,5,32,32,32);  call constructor
        isf.hbar = 0.1;                 specify Planck constant
        isf.dt   = 1/24;                specify time step
        isf.BuildSchroedinger;          this command builds coeff for solving
                                        Schroedinger equation in Fourier domain
          Useful functions:
            [psi1,psi2] = isf.SchroedingerFlow(psi1,psi2)
            solves Schroedinger equation for (psi1,psi2) for isf.dt time.
            [psi1,psi2] = isf.Normalize(psi1,psi2)
            normalizes (psi1,psi2).
            [psi1,psi2] = isf.PressureProject(psi1,psi2)
            projects (psi1,psi2) to satisfy divergence free constraint.
            [vx,vy,vz] = isf.VelocityOneForm(psi1,psi2,isf.hbar)
            extracts velocity 1-from from (psi1,psi2)
            [vx,vy,vz] = isf.VelocityOneForm(psi1,psi2)
            extracts velocity 1-form assuming hbar=1.
            See also TorusDEC
    """
    
    def __init__(self, sizex, sizey, sizez,
                 resx=None, resy=None, resz=None, res=None,
                 hbar=None, dt=None, 
                 SchroedingerMask=None):
        super().__init__(sizex, sizey, sizez,
                         resx, resy, resz, res)
        self.hbar = hbar             # reduced Planck constant
        self.dt   = dt             # time step
        self.SchroedingerMask = SchroedingerMask# Fourier coefficient for solving Schroedinger eq