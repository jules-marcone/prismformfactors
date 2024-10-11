from sasmodels import weights as wgh
import numpy as np
from prismformfactors import functions
from prismformfactors import average as ave

#Most fonctions here only generate a polydispersity on the size parameters, and then compute the functions from the average subpackage.

def poly_vals(param,sigma,npoints:int):
    """Generates different values of a polydispersity on a specified model of variable param including n points, with polydispersity parameter sigma.
    The polydispersity is taken from the sasmodel / SASView library
    More information available here: https://www.sasview.org/docs/user/qtgui/Perspectives/Fitting/pd/polydispersity.html"""
    values,weights = wgh.get_weights("gaussian", npoints, sigma, 1, param, [0,2*param], True)
    return values, weights

def Iiso_nanoprism_polydispersity_valreturn(nsides:int,edge,length,scale,background,sigma,npoints:int,qmin,qmax,nsteps:int,norder:int):
    """Linked polydispersity on the width and length.
    nsides : number of sides of the cross-section
    edge : edge length of the the cross-section edges
    length : prism length
    sigma : 1 polydispersity parameter identical for length and edge length
    npoints : number of points on which polydispersity will be calculated
    qmin, qmax : range of q values on which to compute the model
    nsteps : number of points on which the model will be calculated in the q-range
    norder : order of the lebedev quadrature"""
    # generation of the points of the polydispersity
    q_values=np.logspace(np.log(qmin)/np.log(10),np.log(qmax)/np.log(10),nsteps)
    values, weights = poly_vals(length,sigma,npoints)
    total_intensity = np.array([0. for i in range(nsteps)])
    for i in range(len(values)):
        li = values[i]
        ratio = li/length
        ei = functions.edge_from_gyration_radius(nsides, ratio*functions.gyration_radius_from_edge(nsides, edge))
        total_intensity += weights[i] * np.array(ave.Iiso_nanoprism_fixedq(nsides,ei,li,scale,background,q_values,norder)[1]) 
    return q_values, total_intensity

def Iiso_nanoprism_polydispersity_fixedq(nsides:int,edge,length,scale,background,sigma,npoints:int,q_values,norder:int):
    # generation of the points of the polydispersity
    values, weights = poly_vals(length,sigma,npoints)
    total_intensity = np.array([0. for i in range(len(q_values))])
    for i in range(len(values)):
        li = values[i]
        ratio = li/length
        ei = functions.edge_from_gyration_radius(nsides, ratio*functions.gyration_radius_from_edge(nsides, edge))
        total_intensity += weights[i] * np.array(ave.Iiso_nanoprism_fixedq(nsides,ei,li,scale,background,q_values,norder)[1]) 
    return q_values, total_intensity

def Iiso_nanoprism_polydispersity_length_valreturn(nsides:int,edge,length,scale,background,sigma,npoints:int,qmin,qmax,nsteps:int,norder:int):
    """Polydispersity only on the length of the prism"""
    q_values=np.logspace(np.log(qmin)/np.log(10),np.log(qmax)/np.log(10),nsteps)
    values, weights = poly_vals(length,sigma,npoints)
    total_intensity = np.array([0. for i in range(nsteps)])
    for i in range(len(values)):
        li = values[i]
        total_intensity += weights[i] * np.array(ave.Iiso_nanoprism_fixedq(nsides,edge,li,scale,background,q_values,norder)[1])
    return q_values, total_intensity

def Iiso_nanoprism_polydispersity_length_fixedq(nsides:int,edge,length,scale,background,sigma,npoints:int,q_values,norder:int):
    """Polydispersity only on the length of the prism"""
    values, weights = poly_vals(length,sigma,npoints)
    total_intensity = np.array([0. for i in range(len(q_values))])
    for i in range(len(values)):
        li = values[i]
        total_intensity += weights[i] * np.array(ave.Iiso_nanoprism_fixedq(nsides,edge,li,scale,background,q_values,norder)[1])
    return q_values, total_intensity

def Iiso_nanoprism_polydispersity_width_valreturn(nsides:int,edge,length,scale,background,sigma,npoints:int,qmin,qmax,nsteps:int,norder:int):
    """Polydispersity only on the width"""
    q_values=np.logspace(np.log(qmin)/np.log(10),np.log(qmax)/np.log(10),nsteps)
    values, weights = poly_vals(edge,sigma,npoints)
    total_intensity = np.array([0. for i in range(nsteps)])
    for i in range(len(values)):
        ei = values[i]
        total_intensity += weights[i] * np.array(ave.Iiso_nanoprism_fixedq(nsides,ei,length,scale,background,q_values,norder)[1])
    return q_values, total_intensity

def Iiso_nanoprism_polydispersity_width_fixedq(nsides:int,edge,length,scale,background,sigma,npoints:int,q_values,norder:int):
    """Polydispersity only on the width"""
    values, weights = poly_vals(edge,sigma,npoints)
    total_intensity = np.array([0. for i in range(len(q_values))])
    for i in range(len(values)):
        ei = values[i]
        total_intensity += weights[i] * np.array(ave.Iiso_nanoprism_fixedq(nsides,ei,length,scale,background,q_values,norder)[1])
    return q_values, total_intensity

def Iiso_nanoprism_polydispersity_full_valreturn(nsides:int,edge,length,scale,background,sigma_edge,sigma_length,npoints:int,qmin,qmax,nsteps:int,norder:int):
    """Two distinct polydispersity parameters on the edge and the length"""
    q_values=np.logspace(np.log(qmin)/np.log(10),np.log(qmax)/np.log(10),nsteps)
    values_edge, weights = poly_vals(edge,sigma_edge,npoints)
    values_length, weights = poly_vals(length,sigma_length,npoints)
    total_intensity = np.array([0. for i in range(nsteps)])
    for i in range(len(values_edge)):
        ei = values_edge[i]
        li = values_length[i]
        wi = weights[i]
        total_intensity += wi * np.array(ave.Iiso_nanoprism_fixedq(nsides,ei,li,scale,background,q_values,norder)[1])
    return q_values, total_intensity

def Iiso_nanoprism_polydispersity_full_fixedq(nsides:int,edge,length,scale,background,sigma_edge,sigma_length,npoints:int,q_values,norder:int):
    """Polydispersity on the width and length with two different polydispersity parameters for edge and length
    nsides : number of sides of the cross-section
    edge : edge length of the the cross-section edges
    length : prism length
    sigma_edge, sigma_length: polydispersity parameters the length and edge length
    npoints : number of points on which polydispersity will be calculated
    qmin, qmax : range of q values on which to compute the model
    nsteps : number of points on which the model will be calculated in the q-range
    norder : order of the lebedev quadrature"""
    values_edge, weights = poly_vals(edge,sigma_edge,npoints)
    values_length, weights = poly_vals(length,sigma_length,npoints)
    total_intensity = np.array([0. for i in range(len(q_values))])
    for i in range(len(values_edge)):
        ei = values_edge[i]
        li = values_length[i]
        wi = weights[i]
        total_intensity += wi * np.array(ave.Iiso_nanoprism_fixedq(nsides,ei,li,scale,background,q_values,norder)[1])
    return q_values, total_intensity