from prismformfactors import functions
from prismformfactors import wuttke2d
import numpy as np
from matplotlib import pyplot as plt
from prismformfactors import control as ctrl
from prismformfactors import formfactor
from prismformfactors import average as ave
from prismformfactors import dataload as data
from prismformfactors import polydispersity as pd

def formfactorgraph(n:int,radius):
    vertices = formfactor.shape_generator(n,radius)
    angles = [2*i*np.pi/(1000) for i in range(1000)] #calculations conducted for one fold of the symmetry
    q = [[np.cos(k),np.sin(k)] for k in angles]
    I = [wuttke2d.formfactor_intensity(vertices,k,0) for k in q]
    plt.plot(angles,I)
    plt.xlabel('angle (rad)',fontsize=16)
    plt.ylabel('I(q)',fontsize=16)
    plt.tick_params(pad=0,labelsize=14)
    plt.show()

def Ixy_disk_rod_2Dplot(phideg,thetadeg,psideg,R,L,qmax,nstep:int):
   
    Qx=np.linspace(-qmax,qmax,nstep)
    Qy=np.linspace(-qmax,qmax,nstep)

    qqx,qqy=np.meshgrid(Qx,Qy)
    
    plt.pcolormesh(Qx,Qy,np.log(ctrl.Ixy_disk_rod([qqx,qqy],phideg,thetadeg,psideg,R,L)),shading='auto')
    #plt.colorbar()
    plt.xlabel('qx',fontsize=16)
    plt.ylabel('qy',fontsize=16)
    plt.tick_params(pad=0,labelsize=14)
    plt.title("total Intensity (x,y) plane",fontsize=20)
    plt.show()

def Ixy_nanoprism_2Dplot(phideg,thetadeg,psideg,nsides:int,edge,L,qmax,nstep:int):
   
    Qx=np.linspace(-qmax,qmax,nstep)
    Qy=np.linspace(-qmax,qmax,nstep)

    qqx,qqy=np.meshgrid(Qx,Qy)
    
    plt.pcolormesh(Qx,Qy,np.log(formfactor.Ixy_nanoprism([qqx,qqy],phideg,thetadeg,psideg,nsides,edge,L)),shading='auto')
    #plt.colorbar()
    plt.xlabel('q$_{x}$ '+r'(Å$^{-1}$)',fontsize=16)
    plt.ylabel('q$_{y}$ '+r'(Å$^{-1}$)',fontsize=16)
    plt.title("n = "+str(nsides),fontsize=20)
    plt.tick_params(pad=0,labelsize=14)
    plt.show()

def Iiso_nanoprism(nsides:int,edge,L,qmin,qmax,nstep:int,norder:int):
    q_values, I_values = ave.Iiso_nanoprism_valreturn(nsides,edge,L,1,0,qmin,qmax,nstep,norder)
    plt.loglog(q_values,I_values)
    plt.xlabel('q '+r'(Å$^{-1}$)',fontsize=16)
    plt.ylabel('I(q)',fontsize=16)
    plt.tick_params(pad=0,labelsize=14)
    plt.show()

def viewdata(filepath:str,qmin=0.0001,qmax=0.2):
    q, Iq, Eq = data.load_data(filepath)
    plt.loglog(q,Iq,label='experimental')
    plt.xlabel('q '+r'(Å$^{-1}$)',fontsize=16)
    plt.ylabel('I(q)',fontsize=16)
    plt.tick_params(pad=0,labelsize=14)
    plt.legend()
    plt.xlim(left=qmin,right=qmax)
    plt.show()

def Iiso_nanoprism_plusdata(filepath:str,nsides:int,edge,L,scale,background,qmin,qmax,nstep:int,norder:int):
    q_values, I_values = ave.Iiso_nanoprism_valreturn(nsides,edge,L,scale,background,qmin,qmax,nstep,norder)
    q, Iq, Eq = data.load_data(filepath)
    #plt.loglog(q_values,I_values)
    plt.loglog(q_values,I_values,label='simulated')
    plt.loglog(q,Iq,label='experimental')
    plt.xlabel('q '+r'(Å$^{-1}$)',fontsize=16)
    plt.ylabel('I(q)',fontsize=16)
    plt.tick_params(pad=0,labelsize=14)
    plt.legend()
    plt.xlim(left=qmin,right=qmax)
    plt.show()

def Iiso_nanoprism_chisquare(filename:str,nsides:int,edge,L,scale,background,qmin,qmax,nqs:int,norder:int,chi_type:str="dI Data"):
    """Plot the data contained in filename, and the model with the following parameters."""
    q, Iq, Eq = data.load_data(filename)
    f, ax = plt.subplots()
    q_exp=[]
    Iq_exp = []
    for i in range(len(q)):
        if qmin <= q[i] <= qmax:
            q_exp.append(q[i])
            Iq_exp.append(Iq[i])
    q_values = [q_exp[int(len(q_exp)*i/nqs)] for i in range(int(nqs*(len(q_exp)-1)/len(q_exp))+1)]
    Iq_compare = [Iq_exp[int(len(Iq_exp)*i/nqs)] for i in range(int(nqs*(len(Iq_exp)-1)/len(Iq_exp))+1)]
    q_values, I_values = ave.Iiso_nanoprism_fixedq(nsides,edge,L,scale,background,q_values,norder)
    #plt.loglog(q_values,I_values)
    plt.loglog(q_values,I_values,label='simulated')
    plt.loglog(q,Iq,label='experimental')
    props = dict(boxstyle='square',facecolor='white')
    plt.text(0.25,0.25,r'$\chi^2$ = '+f"{functions.chi_square(Iq_compare,I_values,chi_type):.3f}",fontsize=14,horizontalalignment='center',transform=ax.transAxes,bbox=props)
    plt.xlabel('q '+r'(Å$^{-1}$)',fontsize=16)
    plt.ylabel('I(q)',fontsize=16)
    plt.tick_params(pad=0,labelsize=14)
    plt.legend()
    plt.xlim(left=qmin,right=qmax)
    plt.show()

def Iiso_nanoprism_final(filepath:str,nsides:int,edge,L,scale,background,sigma,sigma2,poly_points:int,qmin,qmax,nqs:int,norder:int,poly_type:str="None",chi_type:str="dI Data"):
    """Plot the data contained in filename, and the model with the following parameters.
    The different polydispersity types are : "None", "Length Only", "Width Only", "Linked Length and Width", "Distinct Length and Width"."""
    q, Iq, Eq = data.load_data(filepath)
    f, ax = plt.subplots()
    q_exp=[]
    Iq_exp = []
    for i in range(len(q)):
        if qmin <= q[i] <= qmax:
            q_exp.append(q[i])
            Iq_exp.append(Iq[i])
    q_values = [q_exp[int(len(q_exp)*i/nqs)] for i in range(int(nqs*(len(q_exp)-1)/len(q_exp))+1)]
    Iq_compare = [Iq_exp[int(len(Iq_exp)*i/nqs)] for i in range(int(nqs*(len(Iq_exp)-1)/len(Iq_exp))+1)]
    match poly_type:
        case "None" :
            q_values, I_values = ave.Iiso_nanoprism_fixedq(nsides,edge,L,scale,background,q_values,norder)
        case "Length Only" :
            q_values, I_values = pd.Iiso_nanoprism_polydispersity_length_fixedq(nsides,edge,L,scale,background,sigma,poly_points,q_values,norder)
        case "Width Only" :
            q_values, I_values = pd.Iiso_nanoprism_polydispersity_width_fixedq(nsides,edge,L,scale,background,sigma,poly_points,q_values,norder)
        case "Linked Length and Width" :
            q_values, I_values = pd.Iiso_nanoprism_polydispersity_fixedq(nsides,edge,L,scale,background,sigma,poly_points,q_values,norder)
        case "Distinct Length and Width" :
            if sigma==0 or sigma2==0:
                raise ValueError("One of the two polydispersities is 0. Please select the appropriate polydispersity mode instead.")
            q_values, I_values = pd.Iiso_nanoprism_polydispersity_full_fixedq(nsides,edge,L,scale,background,sigma,sigma2,poly_points,q_values,norder)
        case _:
            print('WARNING : poly_type not recognised. Check that the poly_type entered appears in the following list : "None", "Length Only", "Width Only", "Linked Length and Width", "Distinct Length and Width"')
            print("Now computing model with no polydispersity.")
            q_values, I_values = ave.Iiso_nanoprism_fixedq(nsides,edge,L,scale,background,q_values,norder)
    #plt.loglog(q_values,I_values)
    plt.loglog(q_values,I_values,label='simulated')
    plt.loglog(q,Iq,label='experimental')
    props = dict(boxstyle='square',facecolor='white')
    plt.text(0.25,0.25,r'$\chi^2$ = '+f"{functions.chi_square(Iq_compare,I_values,chi_type):.3f}",fontsize=14,horizontalalignment='center',transform=ax.transAxes,bbox=props)
    plt.xlabel('q '+r'(Å$^{-1}$)',fontsize=16)
    plt.ylabel('I(q)',fontsize=16)
    plt.tick_params(pad=0,labelsize=14)
    plt.legend()
    plt.xlim(left=qmin,right=qmax)
    plt.show()

def Iiso_nanoprism_chisquare_fixedgyr(filename:str,nsides:int,gyr,L,scale,background,qmin,qmax,nqs:int,norder:int,chi_type:str="dI Data"):
    """Plot the data contained in filename, and the model with the following parameters."""
    q, Iq, Eq = data.load_data(filename)
    f, ax = plt.subplots()
    q_exp=[]
    Iq_exp = []
    for i in range(len(q)):
        if qmin <= q[i] <= qmax:
            q_exp.append(q[i])
            Iq_exp.append(Iq[i])
    edge = functions.edge_from_gyration_radius(nsides,gyr)
    q_values = [q_exp[int(len(q_exp)*i/nqs)] for i in range(int(nqs*(len(q_exp)-1)/len(q_exp))+1)]
    Iq_compare = [Iq_exp[int(len(Iq_exp)*i/nqs)] for i in range(int(nqs*(len(Iq_exp)-1)/len(Iq_exp))+1)]
    q_values, I_values = ave.Iiso_nanoprism_fixedq(nsides,edge,L,scale,background,q_values,norder)
    #plt.loglog(q_values,I_values)
    plt.loglog(q_values,I_values,label='simulated')
    plt.loglog(q,Iq,label='experimental')
    props = dict(boxstyle='square',facecolor='white')
    plt.text(0.25,0.25,r'$\chi^2$ = '+f"{functions.chi_square(Iq_compare,I_values,chi_type):.3f}",fontsize=14,horizontalalignment='center',transform=ax.transAxes,bbox=props)
    plt.xlabel('q '+r'(Å$^{-1}$)',fontsize=16)
    plt.ylabel('I(q)',fontsize=16)
    plt.tick_params(pad=0,labelsize=14)
    plt.legend()
    plt.xlim(left=qmin,right=qmax)
    plt.show()

def Iiso_nanoprism_final_fixedgyr(filepath:str,nsides:int,gyr,L,scale,background,sigma,sigma2,poly_points:int,qmin,qmax,nqs:int,norder:int,poly_type:str="None",chi_type:str="dI Data"):
    """Plot the data contained in filename, and the model with the following parameters.
    The different polydispersity types are : "None", "Length Only", "Width Only", "Linked Length and Width", "Distinct Length and Width"."""
    q, Iq, Eq = data.load_data(filepath)
    f, ax = plt.subplots()
    q_exp=[]
    Iq_exp = []
    for i in range(len(q)):
        if qmin <= q[i] <= qmax:
            q_exp.append(q[i])
            Iq_exp.append(Iq[i])
    q_values = [q_exp[int(len(q_exp)*i/nqs)] for i in range(int(nqs*(len(q_exp)-1)/len(q_exp))+1)]
    Iq_compare = [Iq_exp[int(len(Iq_exp)*i/nqs)] for i in range(int(nqs*(len(Iq_exp)-1)/len(Iq_exp))+1)]
    edge = functions.edge_from_gyration_radius(nsides,gyr)
    match poly_type:
        case "None" :
            q_values, I_values = ave.Iiso_nanoprism_fixedq(nsides,edge,L,scale,background,q_values,norder)
        case "Length Only" :
            q_values, I_values = pd.Iiso_nanoprism_polydispersity_length_fixedq(nsides,edge,L,scale,background,sigma,poly_points,q_values,norder)
        case "Width Only" :
            q_values, I_values = pd.Iiso_nanoprism_polydispersity_width_fixedq(nsides,edge,L,scale,background,sigma,poly_points,q_values,norder)
        case "Linked Length and Width" :
            q_values, I_values = pd.Iiso_nanoprism_polydispersity_fixedq(nsides,edge,L,scale,background,sigma,poly_points,q_values,norder)
        case "Distinct Length and Width" :
            if sigma==0 or sigma2==0:
                raise ValueError("One of the two polydispersities is 0. Please select the appropriate polydispersity mode instead.")
            q_values, I_values = pd.Iiso_nanoprism_polydispersity_full_fixedq(nsides,edge,L,scale,background,sigma,sigma2,poly_points,q_values,norder)
        case _:
            print('WARNING : poly_type not recognised. Check that the poly_type entered appears in the following list : "None", "Length Only", "Width Only", "Linked Length and Width", "Distinct Length and Width"')
            print("Now computing model with no polydispersity.")
            q_values, I_values = ave.Iiso_nanoprism_fixedq(nsides,edge,L,scale,background,q_values,norder)
    #plt.loglog(q_values,I_values)
    plt.loglog(q_values,I_values,label='simulated')
    plt.loglog(q,Iq,label='experimental')
    props = dict(boxstyle='square',facecolor='white')
    plt.text(0.25,0.25,r'$\chi^2$ = '+f"{functions.chi_square(Iq_compare,I_values,chi_type):.3f}",fontsize=14,horizontalalignment='center',transform=ax.transAxes,bbox=props)
    plt.xlabel('q '+r'(Å$^{-1}$)',fontsize=16)
    plt.ylabel('I(q)',fontsize=16)
    plt.tick_params(pad=0,labelsize=14)
    plt.legend()
    plt.xlim(left=qmin,right=qmax)
    plt.show()