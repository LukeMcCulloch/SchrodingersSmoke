#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 23:53:24 2019

@author: lukemcculloch
"""
from __future__ import division

import numpy as np
floor = np.floor

from utils import mod, sub2ind

class Particles(object):
    # Particles - class of particle which can be advected by staggered velocity
    # field on a TorusDEC grid, using RK4 method.
    # Velocities are trilinearly interpolated.
    def __init__(self, x, y, z):
        self.x = x  #array of positions in x
        self.y = y  #array of positions in y
        self.z = z  #array of positions in z
        
    def StaggeredAdvect(self, particle,torus,vx,vy,vz,dt) :
        """ 
            advect particle positions using RK4 in a grid torus with
            staggered velocity vx,vy,vz, for dt period of time
        """
        [k1x,k1y,k1z] = particle.StaggeredVelocity(
                                particle.x,particle.y,particle.z,
                                torus,vx,vy,vz)
        
        [k2x,k2y,k2z] = particle.StaggeredVelocity(
                particle.x+k1x*dt/2,particle.y+k1y*dt/2,particle.z+k1z*dt/2,
                torus,vx,vy,vz)
        
        [k3x,k3y,k3z] = particle.StaggeredVelocity(
                particle.x+k2x*dt/2,particle.y+k2y*dt/2,particle.z+k2z*dt/2,
                torus,vx,vy,vz)
            
        [k4x,k4y,k4z] = particle.StaggeredVelocity(
                particle.x+k3x*dt,particle.y+k3y*dt,particle.z+k3z*dt,
                torus,vx,vy,vz)
            
        particle.x = particle.x + dt/6*(k1x+2*k2x+2*k3x+k4x)
        particle.y = particle.y + dt/6*(k1y+2*k2y+2*k3y+k4y)
        particle.z = particle.z + dt/6*(k1z+2*k2z+2*k3z+k4z)
        return
    
    
    def Keep(self, particle,ind) :
        """used for removing particles
        """
        particle.x = particle.x[ind]
        particle.y = particle.y[ind]
        particle.z = particle.z[ind]
        return
    
    
    @staticmethod
    def  StaggeredVelocity(px,py,pz,torus,vx,vy,vz) :
        """
        evaluates velocity at (px,py,pz) in the grid torus with staggered
        velocity vector field vx,vy,vz
        """
        px = mod(px,torus.sizex)
        py = mod(py,torus.sizey)
        pz = mod(pz,torus.sizez)
        
        ix = floor(px/torus.dx) #+ 1
        iy = floor(py/torus.dy) #+ 1
        iz = floor(pz/torus.dz) #+ 1
        ixp = mod(ix,torus.resx) #+ 1
        iyp = mod(iy,torus.resy) #+ 1
        izp = mod(iz,torus.resz) #+ 1
        ind0 = sub2ind(
                    np.asarray([torus.resx,torus.resy,torus.resz],ix,iy,iz)
                    )
        indxp = sub2ind(
                    np.asarray([torus.resx,torus.resy,torus.resz],ixp,iy,iz)
                    )
        indyp = sub2ind(
                    np.asarray([torus.resx,torus.resy,torus.resz],ix,iyp,iz)
                    )
        indzp = sub2ind(
                    np.asarray([torus.resx,torus.resy,torus.resz],ix,iy,izp)
                    )
        indxpyp = sub2ind(
                    np.asarray([torus.resx,torus.resy,torus.resz],ixp,iyp,iz)
                    )
        indypzp = sub2ind(
                    np.asarray([torus.resx,torus.resy,torus.resz],ix,iyp,izp)
                    )
        indxpzp = sub2ind(
                    np.asarray([torus.resx,torus.resy,torus.resz],ixp,iy,izp)
                    )
        
        wx = px - (ix-1)*torus.dx
        wy = py - (iy-1)*torus.dy
        wz = pz - (iz-1)*torus.dz
        ux = (1-wz)*((1-wy)*vx[ind0 ]+wy*vx[indyp  ]) + \
                wz *((1-wy)*vx[indzp]+wy*vx[indypzp])
        uy = (1-wz)*((1-wx)*vy[ind0 ]+wx*vy[indxp  ]) + \
                wz *((1-wx)*vy[indzp]+wx*vy[indxpzp])
        uz = (1-wy)*((1-wx)*vz[ind0 ]+wx*vz[indxp  ]) + \
                wy *((1-wx)*vz[indyp]+wx*vz[indxpyp])
                
        return np.asarray([ux,uy,uz])