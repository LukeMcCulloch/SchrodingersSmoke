
"""
Created on Dec... 25? 2019

@author: lukemcculloch, based on teh Caltech-TU Berlin Code
1st author Albert Chern.
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
        
        
    def BuildSchroedinger(self) :
        # builds coefficients in Fourier space.
        #
        nx=self.resx; ny=self.resy; nz=self.resz
        fac = -4*pi**2*self.hbar
        kx = (self.iix-1.-nx/2.)/(self.sizex)
        ky = (self.iiy-1.-ny/2.)/(self.sizey)
        kz = (self.iiz-1.-nz/2.)/(self.sizez)
        #lambda = fac*(kx**2+ky**2+kz**2)
        avar = fac*(kx**2+ky**2+kz**2)
        
        self.SchroedingerMask = exp(1.j*avar*self.dt/2.)
        return
        
    def SchroedingerFlow(self,psi1,psi2) :
        # solves Schroedinger equation for dt time.
        #
        psi1 = fftshift(fftn(psi1)); psi2 = fftshift(fftn(psi2))
        psi1 = np.multiply( psi1 , self.SchroedingerMask )
        psi2 = np.multiply( psi2 , self.SchroedingerMask )
        psi1 = ifftn(fftshift(psi1)); psi2 = ifftn(fftshift(psi2))
        return np.asarray([psi1,psi2])
        
    def PressureProject(self, psi1,psi2) :
        # Pressure projection of 2-component wave def.
        #
        [vx,vy,vz] = self.VelocityOneForm(psi1,psi2)
        div = self.Div(vx,vy,vz)
        q = self.PoissonSolve(div)
        #[psi1,psi2] = self.GaugeTransform(psi1,psi2,-q)
        #return np.asarray([psi1,psi2])
        return self.GaugeTransform(psi1,psi2,-q)
        
    def VelocityOneForm(self, psi1,psi2,hbar) :
        # extracts velocity 1-form from (psi1,psi2).
        # If hbar argument is empty, hbar=1 is assumed.
        ixp = mod(self.ix,self.resx) + 1.
        iyp = mod(self.iy,self.resy) + 1.
        izp = mod(self.iz,self.resz) + 1.
        vx = angle(conj(psi1) * psi1[ixp,:,:]  
                  +conj(psi2) * psi2[ixp,:,:] )
        vy = angle(conj(psi1) * psi1[:,iyp,:]
                  +conj(psi2) * psi2[:,iyp,:] )
        vz = angle(conj(psi1) * psi1[:,:,izp]
                  +conj(psi2) * psi2[:,:,izp] )
        if nargin<4:
            hbar = 1
            
        vx = vx*hbar
        vy = vy*hbar
        vz = vz*hbar
        return np.asarray([vx,vy,vz])
    
    def AddCircle(self, psi,center,normal,r,d) :
        # adds a vortex ring to a 1-component wave def psi.
        # Inputs center, normal, r specify the circle.
        # Input d specify the thickness around the disk to create a boost
        # in phase. Usually d = 5*dx where dx is grid edge length.
        rx = self.px - center[0]
        ry = self.py - center[1]
        rz = self.pz - center[2]
        #normal = normal/norm(normal,2) #matlab calls for the 2 norm
        normal = normal/np.linalg.norm(normal) #linalg gives the 2-norm by default
        alpha = np.zeros_like(rx)
        z = rx*normal[0] + ry*normal[1] + rz*normal[2]
        inCylinder = (rx**2+ry**2+rz**2 - z**2) < r**2
        inLayerP = (z> 0) & (z<= d/2) & inCylinder
        inLayerM = (z<=0) & (z>=-d/2) & inCylinder
        alpha[inLayerP] = -pi*(2*z(inLayerP)/d - 1.)
        alpha[inLayerM] = -pi*(2*z(inLayerM)/d + 1.)
        #psi = np.multiply( psi , exp(1j*alpha) )
        psi = psi * exp(1j*alpha)
        return psi


    @staticmethod
    def GaugeTransform(psi1,psi2,q) :
        # multiplies exp(i*q) to (psi1,psi2)
        #
        eiq = exp(1j*q)
        psi1 = np.multiply( psi1 , eiq )
        psi2 = np.multiply( psi2 , eiq )
        return np.asarray([psi1,psi2])
    
    
    @staticmethod
    def Hopf(psi1,psi2) :
        # extracts Clebsch variable s=(sx,sy,sz) from (psi1,psi2)
        #
        a = np.real(psi1)
        b = np.imag(psi1)
        c = np.real(psi2)
        d = np.imag(psi2)
        #sx = 2*( a.multiply( c ) + b.multiply( d ) )
        #sy = 2*( a.multiply( d ) - b.multiply( c ) )
        sx = 2*( a*c + b*d )
        sy = 2*( a*d - b*c )
        sz = a**2 + b**2 - c**2 - d**2
        return np.asarray([sx,sy,sz])
    
    
    @staticmethod
    def Normalize(psi1,psi2) :
        # normalizes (psi1,psi2).
        #
        psi_norm = np.sqrt(abs(psi1)**2 + abs(psi2)**2)
        psi1 = psi1/psi_norm
        psi2 = psi2/psi_norm
        return np.asarray([psi1,psi2])