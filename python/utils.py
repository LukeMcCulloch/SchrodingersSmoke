#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 07:28:41 2019

@author: lukemcculloch
"""
import numpy as np


def mod(a,b) :
    """
    alternate syntax for modulo
    """
    return a%b



def sub2ind2D( size, rows, cols ) :
    return rows + (cols-1)*size[0];


def sub2ind( A, i1, i2, i3 ) :
    size = np.shape(A)
    i1 + (i2-1)*size[0] + (i3-1)*size[0]*size[1] 



def ndgrid(a,b,c) :
    """
    matlab's ndgrid is similar but not quite the same as 
    numpy meshgrid
    
    code experimentations shows the
    need to transpose 0 and 1 dimensions
    
    see also:
    https://stackoverflow.com/questions/12402045/
    mesh-grid-functions-in-python-meshgrid-mgrid-ogrid-ndgrid
    
    a = self.ix
    b = self.iy
    c = self.iz
    
    """
    l0 = np.shape(a)[0]
    l1 = np.shape(b)[0]
    l2 = np.shape(c)[0]
    
    e0 = np.zeros((l0,l1,l2),int)
    for i in range(l1):
        for j in range(l2):
            e0[:,i,j] = a
            
    e1 = np.zeros((l0,l1,l2),int)
    for i in range(l1):
        for j in range(l2):
            e1[i,:,j] = b
            
    e2 = np.zeros((l0,l1,l2),int)
    for i in range(l1):
        for j in range(l2):
            e2[i,j,:] = c
    
    #return np.meshgrid( (a,b,c), indexing='xy')
    return e0,e1,e2