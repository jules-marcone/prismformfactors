import numpy as np

def SINCOS(x):
    """return sin(x), cos(x)"""
    return np.sin(x), np.cos(x)

def rot(phi,theta,psi): #Euler convention Z1Y2Z3 for rotation
    S1, C1 = SINCOS(phi)
    S2, C2 = SINCOS(theta)
    S3, C3 = SINCOS(psi)
    R=[[C1*C2*C3-S1*S3, -C3*S1-C1*C2*S3, C1*S2],
       [C1*S3+C2*C3*S1,  C1*C3-C2*S1*S3, S1*S2],
       [-C3*S2,                   S2*S3,    C2]]
    return (np.array(R))