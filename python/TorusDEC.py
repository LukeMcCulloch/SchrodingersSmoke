#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 22:39:20 2019

@author: lukemcculloch
"""

from __future__ import division
import numpy as np
#import scipy as sp
fft = np.fft.fft
ifft = np.fft.ifft
fftn = np.fft.fftn
ifftn = np.fft.ifftn

sin = np.sin
pi = np.pi

from utils import mod

def ndgrid(a,b,c):
    """
    matlab's ndgrid is similar but not quite the same as 
    numpy meshgrid
    
    code experimentations shows the
    need to transpose 0 and 1 dimensions
    
    see also:
    https://stackoverflow.com/questions/12402045/
    mesh-grid-functions-in-python-meshgrid-mgrid-ogrid-ndgrid
    """
    return np.meshgrid(b,a,c)

class TorusDEC(object):
    """
    TorusDEC a class simple 3D grid with basic exterior calculus operations.

    SYNTAX
    
    obj = TorusDEC(sizex,sizey,sizez,resx,resy,resz)
    obj = TorusDEC(sizex,sizey,sizez,res)
    classdef MyPDEProblem < TorusDEC
    
    DESCRIPTION
    TorusDEC is a handle class that an instance is a 3D grid with periodic
    boundaries in x,y,z direction, i.e. a 3-torus. DEC stands for "Discrete
    Exterior Calculus", a set of operations including exterior derivatives,
    codifferentials.
    obj = TorusDEC(sizex,sizey,sizez,resx,resy,resz) creates an instance
    obj, a 3D grid with size sizex, sizey, sizez and number of divisions
    resx,resy,resz in each dimension.
    obj = TorusDEC creates a default empty instance.
    
    obj = TorusDEC(obj2) copies obj2 to obj.

    obj = TorusDEC(sizex,sizey,sizez,res) creates a grid with size
    sizex,sizey,sizez so that edge lengths dx,dy,dz are equal (cubic 
    lattice).  Input res specifies the number of grid points in the longest
    dimension.
    classdef MyPDEProblem < TorusDEC defines MyPDEProblem as a subclass of
    TorusDEC. MyPDEProblem inherits all methods and members of TorusDEC.

    CLASS MEMBERS
    sizex,sizey,sizez - length in each dimension.
        resx, resy,resz - number of divisions in each dimension.
            px, py, pz - positions.  Each of px, py, pz is a 3D array 
                        carrying x,y,z coordinate of each grid vertex. 
                        px ranges from 0 to sizex, and similarly for py,pz.
            dx, dy, dz - edge lengths in each dimension.
            ix, iy, iz - convenient 1D arrays 1:resx, 1:resy, 1:resz.
        iix, iiy, iiz - convenient 3D arrays generated by ndgrid(ix,iy,iz).
    METHODS
    Suppose obj is an instance of TorusDEC.
    Exterior derivatives:
    [vx,vy,vz] = obj.DerivativeOfFunction(f)
        For a function f compute the 1-form df.
        f is a 3D array representing a scalar function on the grid. vx,vy,vz
        is the 1-form df integrated along edges. vx corresonds to edge 
        (i,j,k)->(i+1,j,k) and so on.
    [wx,wy,wz] = obj.DerivativeOfOneForm(vx,vy,vz)
        For a 1-form v compute the 2-form dv.
    f = obj.DerivativeOfTwoForm(wx,wy,wz)
        For a 2-form w compute the 3-form dw.
    Codifferentials:
    
    f = obj.Div(vx,vy,vz)
        For a 1-form v compute the function *d*v.
    Sharp Operator:
    [ux,uy,uz] = obj.Sharp(vx,vy,vz)
        For a 1-form v compute the corresponding vector field v^sharp by
        averaging to vertices
    [ux,uy,uz] = obj.StaggeredSharp(vx,vy,vz)
        For a 1-form v compute the corresponding vector field v^sharp as
        a staggered vector field living on edges
    Poisson Solve:
    u = obj.PoissonSolve(f)
        solves the Poisson equation L u = f, where L is the Laplacian on
        the 3-torus (negative semidefinite convension).  
        u and f has zero mean.
    """
    def __init__(self,
                 sizex, sizey, sizez, # size of grid
                 resx=None, resy=None, resz=None, res=None):   # number of grid points in each dimension)
        self.sizex = sizex
        self.sizey = sizey
        self.sizez = sizez
        if res is not None:
            self.resx = res
            self.resy = res
            self.resz = res
        else:
            self.resx = resx
            self.resy = resy
            self.resz = resz
        
        self.dx = self.sizex/self.resx #using python 3 division!
        self.dy = self.sizey/self.resy
        self.dz = self.sizez/self.resz
        
        self.ix = np.arange(0,self.resx)
        self.iy = np.arange(0,self.resy)
        self.iz = np.arange(0,self.resz)
        
        # numpy y coord == matlab z coord
        # numpy z coord == matlab y coord
        self.iix,self.iiy,self.iiz = ndgrid(self.ix,
                                               self.iy,
                                               self.iz)
        """
        self.iix,self.iiy,self.iiz = np.mgrid[self.ix,
                                                 self.iy,
                                                 self.iz]
        """
        
        self.px = (self.iix)*self.dx
        self.py = (self.iiy)*self.dy
        self.pz = (self.iiz)*self.dz
        
    
    def DerivativeOfFunction(self, f):
        """
        for the function f, compute the 1-form df
        
            f is a 3D array representing a scalar function on the grid.
            
            vx,vy,vz is the 1-form df integrated along edges. 
            
            vx corresonds to edge (i,j,k)->(i+1,j,k) and so on.
        """
        ixp = mod(self.ix,self.resx) + 1.
        iyp = mod(self.iy,self.resy) + 1.
        izp = mod(self.iz,self.resz) + 1.
    	#ixp = mod(self.ix+1,self.resx)
    	#iyp = mod(self.iy+1,self.resy)
    	#izp = mod(self.iz+1,self.resz)
        vx = f[ixp,:,:] - f
        vy = f[:,iyp,:] - f
        vz = f[:,:,izp] - f
        return np.asarray([vx,vy,vz])
        
    
    def DerivativeOfOneForm(self, vx,vy,vz):
        """
        For a 1-form v compute the 2-form dv
        
        f is a 3D array representing a scalar function on the grid.
        
        vx,vy,vz is the 1-form df integrated along edges. 
        
        vx corresonds to edge (i,j,k)->(i+1,j,k) and so on.
        """
        ixp = mod(self.ix,self.resx) + 1.
        iyp = mod(self.iy,self.resy) + 1.
        izp = mod(self.iz,self.resz) + 1.
        #ixp = mod(self.ix+1,self.resx)
        #iyp = mod(self.iy+1,self.resy)
        #izp = mod(self.iz+1,self.resz)
        wx = vy - vy[:,:,izp] + vz[:,iyp,:] - vz
        wy = vz - vz[ixp,:,:] + vx[:,:,izp] - vx
        wz = vx - vx[:,iyp,:] + vy[ixp,:,:] - vy
        return np.asarray([wx,wy,wz])
    
    
    def DerivativeOfTwoForm(self,wx,wy,wz):
        """
        #DerivativeOfTwoForm
        #For a 2-form w compute the 3-form dw
        """
        ixp = mod(self.ix,self.resx) + 1.
        iyp = mod(self.iy,self.resy) + 1.
        izp = mod(self.iz,self.resz) + 1.
        #ixp = mod(self.ix + 1,self.resx)
        #iyp = mod(self.iy + 1,self.resy)
        #izp = mod(self.iz + 1,self.resz)
        f =     wx[ixp,:,:] - wx
        f = f + wy[:,iyp,:] - wy
        f = f + wz[:,:,izp] - wz
        return f
        
    def Div(self,vx,vy,vz):
        """
        #Div
        #For a 1-form v compute the def *d*v
        """
        ixm = mod(self.ix-2,self.resx) #+ 1
        iym = mod(self.iy-2,self.resy) #+ 1
        izm = mod(self.iz-2,self.resz) #+ 1
        #ixm = mod(self.ix-2 + 1,self.resx)
        #iym = mod(self.iy-2 + 1,self.resy)
        #izm = mod(self.iz-2 + 1,self.resz)
        f =     [vx - vx[ixm,:,:]]/(self.dx**2)
        f = f + [vy - vy[:,iym,:]]/(self.dy**2)
        f = f + [vz - vz[:,:,izm]]/(self.dz**2)
        return f
        
    def Sharp(self,vx,vy,vz):
        """
        #Sharp
        #For a 1-form v compute the corresponding vector field v^sharp by
        #averaging to vertices
        """
        ixm = mod(self.ix-2.,self.resx) + 1.
        iym = mod(self.iy-2.,self.resy) + 1.
        izm = mod(self.iz-2.,self.resz) + 1.
        #ixm = mod(self.ix-2 + 1,self.resx)
        #iym = mod(self.iy-2 + 1,self.resy)
        #izm = mod(self.iz-2 + 1,self.resz)
        ux = 0.5*( vx[ixm,:,:] + vx )/self.dx
        uy = 0.5*( vy[:,iym,:] + vy )/self.dy
        uz = 0.5*( vz[:,:,izm] + vz )/self.dz
        return  np.asarray([ux,uy,uz])
        
    def StaggeredSharp(self,vx,vy,vz):
        """
        #StaggeredSharp
        #For a 1-form v compute the corresponding vector field v^sharp as
        #a staggered vector field living on edges
        """
        ux = vx/self.dx
        uy = vy/self.dy
        uz = vz/self.dz
        return  np.asarray([ux,uy,uz])
        
    def PoissonSolve(self,f):
        """
        #PoissonSolve by Spectral method
        """
        f = fftn(f)
        sx = sin(pi*(self.iix-1)/self.resx)/self.dx
        sy = sin(pi*(self.iiy-1)/self.resy)/self.dy
        sz = sin(pi*(self.iiz-1)/self.resz)/self.dz
        denom = sx**2 + sy**2 + sz**2
        fac = -0.25/denom
        fac[0,0,0] = 0
        f = f * fac
        #f = np.multiply(f,fac)
        f = ifftn(f)
        return f
        
        
        
if __name__ == '__main__':
    sizex = 10
    sizey = 10
    sizez = 10
    res = 10
    self = TorusDEC(sizex=sizex,
                   sizey=sizey,
                   sizez=sizez,
                   res=res)