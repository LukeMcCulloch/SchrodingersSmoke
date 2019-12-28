#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 23:53:24 2019

@author: lukemcculloch
"""
from __future__ import division

class Particles(object):
    # Particles - class of particle which can be advected by staggered velocity
    # field on a TorusDEC grid, using RK4 method.
    # Velocities are trilinearly interpolated.
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        
    def StaggeredAdvect(self, particle,torus,vx,vy,vz,dt) :
        """ 
            advect particle positions using RK4 in a grid torus with
            staggered velocity vx,vy,vz, for dt period of time
        """
        
        return
    
    
    def Keep(self, particle,ind) :
        """for removing particles
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
        ind0 = sub2ind([torus.resx,torus.resy,torus.resz],ix,iy,iz)
        indxp = sub2ind([torus.resx,torus.resy,torus.resz],ixp,iy,iz)
        indyp = sub2ind([torus.resx,torus.resy,torus.resz],ix,iyp,iz)
        indzp = sub2ind([torus.resx,torus.resy,torus.resz],ix,iy,izp)
        indxpyp = sub2ind([torus.resx,torus.resy,torus.resz],ixp,iyp,iz)
        indypzp = sub2ind([torus.resx,torus.resy,torus.resz],ix,iyp,izp)
        indxpzp = sub2ind([torus.resx,torus.resy,torus.resz],ixp,iy,izp)
        
        wx = px - (ix-1)*torus.dx
        wy = py - (iy-1)*torus.dy
        wz = pz - (iz-1)*torus.dz
        ux = (1-wz).*((1-wy).*vx(ind0 )+wy.*vx(indyp  )) + ...
                wz .*((1-wy).*vx(indzp)+wy.*vx(indypzp))
        uy = (1-wz).*((1-wx).*vy(ind0 )+wx.*vy(indxp  )) + ...
                wz .*((1-wx).*vy(indzp)+wx.*vy(indxpzp))
        uz = (1-wy).*((1-wx).*vz(ind0 )+wx.*vz(indxp  )) + ...
                wy .*((1-wx).*vz(indyp)+wx.*vz(indxpyp))
                
        return np.asarray([ux,uy,uz])