# -*- coding: utf-8 -*-
"""
Created on Sat Nov 05 21:17:19 2016

@author: Charlie Young
"""


from __future__ import division
import numpy as np

A = np.matrix([[0,-20,-14],
              [3,27,-4],
              [4,11,-2]])


def norm(x,p=2):
    """p norm of a vector"""
    return np.sum(np.array(x)**p)**(1./p)
    
m,n = A.shape
P = {}

i = 0; j = 0
while j < min(m,n):
    
    P_ = np.identity(m)
    for i in range(j+1,m):
        Pij = np.identity(m)

        c = float(A[j,j]/np.sqrt(A[j,j]**2+A[i,j]**2))
        s = float(A[i,j]/np.sqrt(A[j,j]**2+A[i,j]**2))
        Pij[j,j] = c
        Pij[j,i] = s
        Pij[i,j] = -s
        Pij[i,i] = c

        A = np.dot(Pij,A)
        P_ = np.dot(P_,Pij).T

    P[j] = P_
    j += 1

Q = np.identity(m)
for k in range(len(P))[::-1]:
    Q = np.dot(Q,P[k])


