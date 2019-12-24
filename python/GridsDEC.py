#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 05:18:31 2019

@author: lukemcculloch
"""


class GridDEC(object):
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
        pass
    