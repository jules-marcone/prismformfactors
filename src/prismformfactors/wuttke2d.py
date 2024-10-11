import numpy as np
from prismformfactors import functions
from prismformfactors import rotation as rot

def scalar_product(u,v):
    """Parameters : u and v two vectors as lists of their coordinates
    Preconditions : u and v are the same length
    Result : the scalar product of u and v"""
    somme = 0
    for i in range(len(u)):
        somme += u[i]*v[i]
    return somme

def edgecenters_generator(vertices:list):
    """Parameters : Vertices, list of the vertices of the 2D-shape each listed as a list of the 2D coordinates of the shape in the plane
    Preconditions : a loop, all in the plane
    Result : the list of the edge centers of the 2D shape"""
    extended_vertices = [] + vertices
    extended_vertices.append(vertices[0])
    edgecenter = []
    for i in range(len(vertices)):
        coordinates=[]
        for j in range(2):
            coordinates.append((extended_vertices[i+1][j]+extended_vertices[i][j])/2)
        edgecenter.append(coordinates)
    return edgecenter

def halfedges_generator(vertices:list):
    """Parameters : Vertices, list of the vertices of the 2D-shape each listed as a list of the 2D coordinates of the shape in the plane
    Preconditions : a loop, all in the plane
    Result : the list of the edge centers of the 2D shape"""
    extended_vertices = [] + vertices
    extended_vertices.append(vertices[0])
    halfedge = []
    for i in range(len(vertices)):
        coordinates=[]
        for j in range(2):
            coordinates.append((extended_vertices[i+1][j]-extended_vertices[i][j])/2)
        halfedge.append(coordinates)
    return halfedge

def formfactor(vertices:list, q, c): # gives the area of the polygon at q==[0,0]
    """Parameters : Vertices, list of the vertices of the 2D-shape each listed as a list of the 2D coordinates of the shape in the plane
    q (is q//), listed as a list of coordinates in the plane
    c an arbitrary constant
    Preconditions : vertice a loop, all in the plane
    Result : the complex form factor f(q) of the shape defined by the vertices."""
    qmodulus2=q[0]**2+q[1]**2
    qmodulus2=np.asanyarray(qmodulus2) #conversion to array
    cutoff=0.000000000000000000000000000001
    qmodulus2_cutoff=np.where(qmodulus2==0,cutoff,qmodulus2) # replace by cutoff if equal to 0
    #This case starts occuring for computations of the prism formfactor, where during the
    #calculation of the orientational average, q// may be equal to 0.
    #if qmodulus2==0:
    #    return 0
    edgecenters = edgecenters_generator(vertices)
    halfedges = halfedges_generator(vertices)
    sum = 0
    for i in range(len(vertices)):
        
        qEj = scalar_product(q, halfedges[i])

        triple_product=q[0]*halfedges[i][1]-q[1]*halfedges[i][0]

        #The exp(iqRj) is rewritten as a sum of cos+isin because of the way math.exp() works (not allowing complex as input)
        qRj = scalar_product(q, edgecenters[i])
        sum += triple_product * (functions.sinc(qEj)*(np.cos(qRj)+np.sin(qRj)*1J)-c)
    return(2/(1J*qmodulus2_cutoff)*sum)
#    return(2/(1J*qmodulus2)*sum)

def formfactor_intensity(vertices:list, q, c):
    """Parameters : vertices, list of the vertices of a 2D-shape, listed as a list of 2D-coordinates
    q is the planar component of the scattering vector, as a list of its coordinates
    c is an arbitrary constant
    Preconditions : vertices form a loop, and are all contained in the same 2D plane
    Result : the scattered intensity I(q) = f(q)f*(q)"""
    f = formfactor(vertices, q, c)
#    return f*np.conj(f)
    return (np.abs(f))**2

