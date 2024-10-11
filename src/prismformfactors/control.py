import numpy as np
from prismformfactors import functions
import scipy as sp
from prismformfactors import rotation as rot

def Amp_abc_disk_rod(q,R,L): # q with three components
    qa,qb,qc=q[0],q[1],q[2] # dimensions of the rod: disk of radius R and length L
    qpar=np.sqrt(qa*qa+qb*qb)
    qpar=np.asanyarray(qpar) #conversion to array
    cutoff=0.000000000000000000000000000001
    qpar_cutoff=np.where(qpar==0,cutoff,qpar) # replace by cutoff if equal to 0
    return 2*sp.special.jv(1, qpar_cutoff*R)/(qpar_cutoff*R)*functions.sinc(qc*L/2)

def I_disk_rod(q,R,L): # Intensity in 3D q with three components
    A=Amp_abc_disk_rod(q,R,L)
    I=(np.abs(A))**2
    I=I*functions.volume_disk_rod(R,L)**2 # multiplication by the particle volume at the power 2
    return I 

def Amp_abc_square_rod(q,a,L): # q with three components
    qa,qb,qc=q[0],q[1],q[2] # dimensions of the square rod are: a x a x L
    return functions.sinc(qa*a/2)*functions.sinc(qb*a/2)*functions.sinc(qc*L/2)

def volume_square_rod(a,L):
    return a*a*L

def I_square_rod(q,a,L): # Intensity in 3D q with three components
    A=Amp_abc_square_rod(q,a,L)
    I=(np.abs(A))**2
    I=I*(volume_square_rod(a,L))**2 # multiplication by the volume at the power 2
    return I 

def Ixy_disk_rod(qxy,phideg,thetadeg,psideg,R,L): # Intensity in xy plane with rotation 
    phi=phideg*np.pi/180 #conversion to radians
    theta=thetadeg*np.pi/180 #conversion to radians
    psi=psideg*np.pi/180 #conversion to radians
    # image is calculated in the (v1, v2) plane
    v1x,v1y,v1z =np.dot(rot.rot(phi,theta,psi),[1,0,0]) # v1 : rotation of the [1,0,0] direction
    v2x,v2y,v2z=np.dot(rot.rot(phi,theta,psi),[0,1,0])  # v2 : rotation of the [0,1,0] direction
    qrotx=v1x*qxy[0]+v2x*qxy[1] # x-component of the rotated q-vector
    qroty=v1y*qxy[0]+v2y*qxy[1] # y-component of the rotated q-vector
    qrotz=v1z*qxy[0]+v2z*qxy[1] # z-component of the rotated q-vector
    qrot=[qrotx,qroty,qrotz] # rotated q-vector
    return I_disk_rod(qrot,R,L)

def I_pave(q,edge,length):
    qa,qb,qc =q[0],q[1],q[2]
    A = np.sinc(qa*edge/2/np.pi)*np.sinc(qb*edge/2/np.pi)*np.sinc(qc*length/2/np.pi)
    I=(np.abs(A))**2
    I=I*(edge*edge*length)**2
    return I

