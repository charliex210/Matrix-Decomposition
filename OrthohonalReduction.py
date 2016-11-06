# -*- coding: utf-8 -*-
"""
Created on Sat Nov 05 21:17:19 2016

@author: Charlie Young
"""

from __future__ import division
import numpy as np

def norm(x,p=2):
    """p norm of a vector"""
    return np.sum(np.array(x)**p)**(1./p)

def GivensReduction(A,precision=1e-8):
    """orthohonal reduction of A using Givens reduction (plane rotation)"""
    
    A = np.matrix(A)
    m,n = A.shape
    T = A.copy()
    P = []
    
    i = 0; j = 0
    while j < min(m,n):
        
        for i in range(j+1,m):
            Pij = np.identity(m)
    
            c = float(T[j,j]/np.sqrt(T[j,j]**2+T[i,j]**2))
            s = float(T[i,j]/np.sqrt(T[j,j]**2+T[i,j]**2))
            Pij[j,j] = c
            Pij[j,i] = s
            Pij[i,j] = -s
            Pij[i,i] = c
    
            T = np.dot(Pij,T)
            P.append(Pij)
    
        j += 1
    
    Q = np.identity(m)
    for k in range(len(P))[::-1]:
        Q = np.dot(Q,P[k])
    
    if precision:
        Q[abs(Q)<precision] = 0
        T[abs(T)<precision] = 0
        
    return Q,T

if __name__ == '__main__':
    A = np.matrix([[0,-20,-14],
                   [3,27,-4],
                   [4,11,-2]])
    
    Q,T=GivensReduction(A)
    print "Q: \n", Q
    print "T: \n", T
