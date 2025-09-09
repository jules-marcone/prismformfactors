from pylebedev import PyLebedev
import numpy as np
from prismformfactors import controls as ctrl
from prismformfactors import formfactor
from matplotlib import pyplot as plt
from prismformfactors import functions

orderlist=[3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131]

leblib = PyLebedev()

def Iiso_disk_valreturn(R,L,qmin,qmax,nstep:int,norder:int):
    """Parameters : R, radius of the cylinder
    L : Length of the cylinder
    qmin, qmax : range of q values for which the model is computed
    nstep : number of points for which the model is calculated
    norder : precision of Lebedev quadrature used (ranging between 0 for order 3, and 31 for order 131)
    Returns the scattered intesity of a cylinder as two lists : q_values (values of q for which the model is computed, averaged in all directions) and I_values (computed intensities)."""
    q_values=np.logspace(np.log(qmin)/np.log(10),np.log(qmax)/np.log(10),nstep)
    I_values=[]
    # order can take 32 predefined values between 3 and 131 listed in orderlist
    order=orderlist[norder] #norder is between 0 and 31
    q_unit,w = leblib.get_points_and_weights(order) #integration points and associated weight on the unit sphere
    for q_modulus in q_values:
        integral = 1 * np.sum(w * ctrl.I_disk_rod([q_modulus*q_unit[:,0],q_modulus*q_unit[:,1],q_modulus*q_unit[:,2]],R,L))
        I_values.append(integral)
    #plt.loglog(q_values,I_values)
    return q_values, I_values 

def Iiso_disk_fixedq(R,L,scale,background,q_values,norder:int):
    """Parameters : R, radius of the cylinder
    L : Length of the cylinder
    q_values : values for which I(q) is computed
    nstep : number of points for which the model is calculated
    norder : precision of Lebedev quadrature used (ranging between 0 for order 3, and 31 for order 131)
    Returns the scattered intesity of a cylinder as two lists : q_values (values of q for which the model is computed, averaged in all directions) and I_values (computed intensities)."""
    I_values=[]
    # order can take 32 predefined values between 3 and 131 listed in orderlist
    order=orderlist[norder] #norder is between 0 and 31
    q_unit,w = leblib.get_points_and_weights(order) #integration points and associated weight on the unit sphere
    for q_modulus in q_values:
        integral = 1 * np.sum(w * ctrl.I_disk_rod([q_modulus*q_unit[:,0],q_modulus*q_unit[:,1],q_modulus*q_unit[:,2]],R,L)) / functions.volume_disk_rod(R,L)**2
        I_values.append(scale*integral+background)
    #plt.loglog(q_values,I_values)
    return q_values, I_values 

def Iiso_pave_valreturn(edge,length,qmin,qmax,nstep:int,norder:int):
    """Parameters : edge, edge length of the cross-section of the square prism
    length : length of the square prism
    qmin, qmax : range of q values for which the model is computed
    nstep : number of points for which the model is calculated
    norder : precision of Lebedev quadrature used (ranging between 0 for order 3, and 31 for order 131)
    Returns the scattered intesity of a square prism as two lists : q_values (values of q for which the model is computed, averaged in all directions) and I_values (computed intensities)."""
    q_values = np.logspace(np.log(qmin)/np.log(10),np.log(qmax)/np.log(10),nstep)
    I_values=[]
    order=orderlist[norder]
    q_unit,w=leblib.get_points_and_weights(order)
    for q_modulus in q_values:
        integral = 1 * np.sum(w*ctrl.I_pave([q_modulus*q_unit[:,0],q_modulus*q_unit[:,1],q_modulus*q_unit[:,2]],edge,length))
        I_values.append(integral)
    return q_values, I_values

def Iiso_nanoprism_valreturn(nsides:int,edge,L,scale,background,qmin,qmax,nstep:int,norder:int):
    """Parameters : nsides : number of sides of the regular polygon cross-section
    edge : length of the sides of the cross section
    L : Length of the prism
    qmin, qmax : range of q values for which the model is computed
    nstep : number of points for which the model is calculated
    norder : precision of Lebedev quadrature used (ranging between 0 for order 3, and 31 for order 131)
    Returns the scattered intesity of a prism of regular cross-section as two lists : q_values (values of q for which the model is computed, averaged in all directions) and I_values (computed intensities)."""
    q_values=np.logspace(np.log(qmin)/np.log(10),np.log(qmax)/np.log(10),nstep)
    I_values=[]
    # order can take 32 predefined values between 3 and 131 listed in orderlist
    order=orderlist[norder] #norder is between 0 and 31
    q_unit,w = leblib.get_points_and_weights(order) #integration points and associated weight on the unit sphere
    for q_modulus in q_values:
        integral = 1 * np.sum(w * formfactor.I_nanoprism([q_modulus*q_unit[:,0],q_modulus*q_unit[:,1],q_modulus*q_unit[:,2]],nsides,edge,L))/(functions.volume_nanoprism(nsides,edge,L))**2
        I_values.append(scale*integral+background)
    return q_values, I_values

def Iiso_nanoprism_fixedq(nsides:int,edge,L,scale,background,q_values,norder:int):
    "From a list of q_values, compute all related scattered intensities (with averageing) and return the list of q-values and the list of I-values."
    I_values=[]
    order=orderlist[norder]
    q_unit,w = leblib.get_points_and_weights(order)
    for q_modulus in q_values:
        integral = 1 * np.sum(w * formfactor.I_nanoprism([q_modulus*q_unit[:,0],q_modulus*q_unit[:,1],q_modulus*q_unit[:,2]],nsides,edge,L))/(functions.volume_nanoprism(nsides,edge,L))**2
        I_values.append(scale*integral+background)
    return q_values, I_values
