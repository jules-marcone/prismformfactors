from prismformfactors import wuttke2d
from prismformfactors import functions
import numpy as np
from prismformfactors import rotation as rot

def shape_generator(number_of_sides, radius):
    """Return the list of vertices of a n-sided regular polygon, of circumradius radius."""
    vertices = [[radius*np.cos(2*step*np.pi/number_of_sides), radius*np.sin(2*step*np.pi/number_of_sides)] for step in range(0,number_of_sides)]
    return vertices


def Amp_nanoprism(q,nsides,edge,L): # From factor in 3D of the nanoprism q with three components
    """Returns the form factor amplitude of a prism with a n-sided regular polygon cross-section and length L for a specific three dimensional q."""
    qa,qb,qc=q[0],q[1],q[2] 
    qab=[qa,qb]
    radius=functions.radius_from_edge(nsides,edge)
    vertices = shape_generator(nsides,radius)
#    A=formfactor(vertices,qab,0.)*sinc(qc*L/2)/surface_from_radius(nsides,radius)
    A=wuttke2d.formfactor(vertices,qab,0.)*functions.sinc(qc*L/2) # form factor is proportionnal to the area
    return A

def I_nanoprism(q,nsides,edge,L): # Intensity in 3D q with three components
    """Returns the scattered intensity of a prism with a n-sided regular polygon cross-section and length L for a specific three dimensional q."""
    A=Amp_nanoprism(q,nsides,edge,L)
    I=(np.abs(A))**2
#    I=I*(volume_nanoprism(nsides,edge,L))**2 # multiplication by the volume at the power 2
    I=I*(L)**2 # multiplication by the volume at the power 2
    return I  

def Ixy_nanoprism(qxy,phideg,thetadeg,psideg,nsides,edge,L): # Intensity in xy plane with rotation 
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
    return I_nanoprism(qrot,nsides,edge,L)