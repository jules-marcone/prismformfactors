import numpy as np
def sinc(x):
    """Returns the cardinal sinus of x."""
    return np.sinc(x/np.pi) #in numpy, sinc is defined with pi*x for the argument

def radius_from_surface(nsides:int, surface):
    """Returns the circumradius of an n-sided regular polygon of area surface."""
    return np.sqrt(surface/(nsides*np.sin(np.pi/nsides)*np.cos(np.pi/nsides)))

def surface_from_radius(nsides:int, radius):
    """Returns the area of an n-sided regular polygon of circumradius radius."""
    return nsides*(radius**2)*np.sin(np.pi/nsides)*np.cos(np.pi/nsides)

def radius_from_edge(nsides:int, edge):
    """Returns the circumradius of an n-sided regular polygon of edge length edge."""
    return edge/(2*np.sin(np.pi/nsides))

def edge_from_radius(nsides:int, radius):
    """Returns the edge length of an n-sided regular polygon of circumradius radius."""
    return 2*radius*np.sin(np.pi/nsides)

def edge_from_area(nsides:int,area):
    """Returns the edge length of an n-sided regular polygon of area area."""
    return 2*np.sqrt(area/nsides * np.tan(np.pi/nsides))

def radius_from_gyration_radius(nsides:int,gyrrad):
    """Returns the circumradius of an n-sided regular polygon of gyration radius gyrrad."""
    return gyrrad / np.sqrt(np.sinc(2/nsides))

def gyration_radius_from_radius(nsides:int,radius):
    """Returns the gyration radius of an n-sided regular polygon of circumradius radius."""
    return radius * np.sqrt(np.sinc(2/nsides)) #np.sinc(x) calcule sin(pi*x)/(pi*x)

def gyration_radius_from_edge(nsides:int, edge):
    """Returns the gyration radius of an n-sided regular polygon of edge length edge."""
    return (edge / (2*np.sin(np.pi/nsides))) * np.sqrt(np.sinc(2/nsides))

def edge_from_gyration_radius(nsides:int,gyr):
    """Returns the edge length of an n-sided regular polygon of gyration ratio gyr."""
    return (gyr * 2 * np.sin(np.pi/nsides)) / np.sqrt(np.sinc(2/nsides))

def volume_nanoprism(nsides:int,edge,L):
    """Returns the volume of a nanoprism with a cross-section in the shape of an n-sided regular polygon of edge length edge, and length L."""
    radius=radius_from_edge(nsides,edge)
    surface=surface_from_radius(nsides,radius)
    return surface*L

def volume_disk_rod(R,L):
    """Returns the volume of a cylinder of cross-section radius R, and of length L."""
    return np.pi*R*R*L

def chi_square(Lexp:list,Ltheo:list,weight_type:str="dI Data"):
    """Calculates the chi-square given two lists with the same number of arguments. 
    Lexp is the experimental data set, Ltheo is theoretical model data set.
    weight-type parameters can be:
    "None"
    "dI Data": artificially created by creating error bars as 10% of the value of the experimental data.
    "sqrt(I Data)"
    "I Data"
    """
    match weight_type:
        case "None":
            return np.sum((np.array(Lexp)-np.array(Ltheo))**2)/len(Lexp)
        case "dI Data":
            return np.sum(((np.array(Lexp)-np.array(Ltheo))/(0.1*np.array(Lexp)))**2)/len(Lexp)
        case "sqrt(I Data)":
            return np.sum((np.array(Lexp)-(np.array(Ltheo)))**2/(np.array(Lexp)))/len(Ltheo)
        case "I Data":
            return np.sum(((np.array(Lexp)-(np.array(Ltheo)))/(np.array(Lexp)))**2)/len(Ltheo)
        case _:
            print("WARNING : weight_type for chi² calculation not recognized. Check that the weight_type entered belongs to the following list : 'None', 'dI Data', 'sqrt(dI Data)', 'I Data'.")
            print("Now computing the chi-square using dI Data.")
            return np.sum(((np.array(Lexp)-np.array(Ltheo))/(0.1*np.array(Lexp)))**2)/len(Lexp)
