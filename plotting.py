import numpy as np 
from numba import jit
import matplotlib.pyplot as plt
import os

from force_energy import dist,sigma_h

def save_into_path(TITLE=None,save_path="./plot/"):
    """Save into a certain directory

    Args:
        TITLE (str, optional): filename of image. Defaults to None.
        save_path (str, optional): directory name. Defaults to "./plot/".
    """
    os.makedirs(name=save_path,exist_ok=True)
    plt.savefig(str(save_path)+str(TITLE)+".png")


def plot_block(coords,stopstep,T_init=None):
    """Plot the block of atoms with their coordinates.

    Args:
        coords (numpy.ndarray): A 2D array containing the x and y coordinates of atoms in the block.
        stopstep (int): The iteration number at which the simulation stopped.
        T_init(int or float): T_desire of the system. Just for plot's title. (Defaults to 'None')
    """
    n=int((len(coords)**0.5))
    if T_init == None:
        TITLE = "Initial " + str(n) + '*' + str(n) +" Block"
    elif stopstep==0:
        TITLE = "Initial " + str(n) + '*' + str(n) +" Block at " + str(T_init) + "K"
    else:
        TITLE = "MD-" + str(T_init) + "K relaxed " + str(n) + '*' + str(n) + " Block, iter=" +str(stopstep)
    plt.figure(figsize=(4,4))
    x = [i[0] for i in coords]
    y = [i[1] for i in coords]
    plt.scatter(x,y,s=8)
    plt.xlabel('x coordinate')
    plt.ylabel('y coordinate')
    plt.title(TITLE)
    save_into_path(TITLE,save_path="./plot/")
    plt.show()


@jit
def rdf(coords, RDF_MAX=12, n=400, RHO=400./((19*(3.405*(2**(1/6))))**2), PI=np.pi):
    """Calculate the radial distribution function (RDF) for a given set of atomic coordinates.

    Args:
        coords (numpy.ndarray): A 2D array containing the atomic coordinates.
        RDF_MAX (int, optional): Maximum distance to consider for the RDF calculation. Defaults to 12.
        n (int, optional): Number of atoms in the system. Defaults to 400.
        RHO (float, optional): Number density of the atomic system. Defaults to 400./((19*(3.405*(2**(1/6))))**2).
        PI (float, optional): The value of pi. Defaults to np.pi.

    Returns:
        tuple: A tuple includes a ndarray of radial distances (r) and the corresponding RDF values.
    """
    r = np.linspace(1,RDF_MAX,n)
    rdf_result = []
    dr = 0.1
    for r_value in r:
        count = 0   
        for i in range(n):
            for j in range(n):
                if i!=j:
                    if r_value <= dist(coords[i],coords[j]) <= r_value+dr :
                        count += 1
        rdf_result.append(count/(2*PI*r_value*dr*n*RHO))
    return r, rdf_result


def plot_rdf(coords,TITLE='Initial RDF'):
    """Plot the radial distribution function (RDF) of the system.

    Args:
        coords (numpy.ndarray): A 2D array containing the x and y coordinates of atoms in the system.
        title (str, optional): The title for the plot. Defaults to 'Initial RDF'.
    """
    r,rdf_result=rdf(coords)
    rdf_result = list(rdf_result)
    print('length of r array')
    print(len(r))
    print('begining and end of r array')
    print(r[0],r[-1])
    plt.figure(figsize=(8,5))
    plt.bar(r,rdf_result,width=0.2)
    plt.xlim((2,13))
    plt.title(TITLE)
    plt.xlabel('Distance (Å)')
    plt.ylabel('Radial distribution function [# atoms]')
    save_into_path(TITLE,save_path="./plot/")
    plt.show()


def Plot_Ep_Ek(E_p,E_k,stopstep,T_init=None):
    """Plot the potential, kinetic, and total energies of the system versus the number of iterations.

    Args:
        E_p (list): A list of potential energies at each iteration.
        E_k (list): A list of kinetic energies at each iteration.
        stopstep (int): The iteration number at which the simulation stopped.
        T_init (float, optional): The initial temperature of the system. Defaults to None.
    """
    xlist = []
    for i in range(stopstep+1):
        if i%50 ==0:
            xlist.append(i)
        if len(xlist)==len(E_p):
            break
    E_tot = np.array(E_p)+np.array(E_k)
    plt.figure(figsize=(7,5))
    plt.plot([i for i in xlist], [i for i in E_p],label='Potential Energy')
    plt.plot([i for i in xlist], [i for i in E_k],label='Kinetic Energy')
    plt.plot([i for i in xlist], [i for i in E_tot],label='Total Energy')
    if T_init !=None:
        TITLE = "Total energy versus iters, T=" +str(T_init)+"K"
    else:
        TITLE = "Total energy versus iters"
    plt.title(TITLE) 
    plt.xlabel('Iteration')
    plt.ylabel('Total Energy (eV)')
    plt.legend()
    save_into_path(TITLE,save_path="./plot/")
    plt.show()


def Plot_temperature(stopstep,T,T_init=None):
    """Plot the temperature of the system versus the number of iterations.

    Args:
        stopstep (int): The iteration number at which the simulation stopped.
        T (list): Temperature list appended and collected in main loop.
        T_init (float, optional): The initial temperature (T_desired) of the system. Defaults to None.
    """
    xlist = []
    for i in range(stopstep+1):
        if i%200 ==0:
            xlist.append(i)
        if len(xlist)==len(T):
            break
    plt.figure(figsize=(7,5))
    if T_init !=None:
        TITLE = "Temperature versus iters, T_init=" +str(T_init)+ 'K'
    else:
        TITLE = "Temperature versus iters"
    plt.title(TITLE)
    plt.xlabel('Iterations')
    plt.ylabel('Temperature (K)')
    plt.plot([i for i in xlist], [i for i in T])
    save_into_path(TITLE,save_path="./plot/")
    plt.show()


def Plot_HS_tot(stopstep,hydrostress,T_init=None):
    """Plot the total hydrostatic stress versus the number of iterations.

    Args:
        stopstep (int): The iteration number at which the simulation stopped.
        hydrostress (list): A list of total hydrostatic stress values at each iteration.
        T_init (float, optional): The initial temperature of the system. Defaults to None.
    """
    xlist = []
    for i in range(stopstep+1):
        if i%200 ==0:
            xlist.append(i)
        if len(xlist)==len(hydrostress):
            break
    plt.figure(figsize=(7,5))
    if T_init !=None:
        TITLE = "Hydrostatic Pressure versus Iters, T="+str(T_init)+'K'
        
    else:
        TITLE = "Hydrostatic Pressure versus Iters"
    plt.title(TITLE)
    plt.xlabel('Iteration')
    plt.ylabel('Stress (eV/Å²)')
    plt.plot([i for i in xlist], [i for i in hydrostress])
    save_into_path(TITLE,save_path="./plot/")
    plt.show()


def Plot_sigma_h(coords,T_init,Type='MD'):
    """Plot the hydrostatic stress distribution in the final atomic structure.

    Args:
        coords (numpy.ndarray): A 2D array containing the x and y coordinates of atoms in the system.
        T_init (float): The initial temperature of the system.
        Type (str, optional): The type of simulation used to generate the structure (e.g., 'MD' for molecular dynamics). Defaults to 'MD'.
    """
    coords_list = []
    for i in range(len(coords)):
        coords_list.append(list(coords[i]))
    plt.figure(figsize=(6.5,5))
    plt.xlabel('X coordinate')
    plt.ylabel('Y coordinate')
    TITLE = "Hydrostatic stress, "+str(Type)+"final structure at "+str(T_init)+'K'
    plt.title(TITLE)
    plt.scatter([i[0] for i in coords_list],[i[1] for i in coords_list],s=8,c=sigma_h(coords))
    plt.colorbar(label='eV/Å²')
    save_into_path(TITLE,save_path="./plot/")
    plt.show()


def Plot_results(coords,E_p,E_k,hydrostress,T,T_init,stopstep):
    """Plot various results of the atomic simulation, including relaxed block structure, energies, RDF, temperature, and hydrostatic stress.

    Args:
        coords (numpy.ndarray): A 2D array containing the x and y coordinates of atoms in the system.
        E_p (list): A list of potential energies at each iteration.
        E_k (list): A list of kinetic energies at each iteration.
        hydrostress (list): A list of total hydrostatic stress values at each iteration.
        T (list): A list of temperatures at each iteration.
        T_init (float): The initial temperature of the system.
        stopstep (int): The iteration number at which the simulation stopped.
    """
    
    plot_block(coords,stopstep,T_init)
    Plot_Ep_Ek(E_p,E_k,stopstep,T_init)
    plot_rdf(coords,title="MD-"+str(T_init)+"K relaxed RDF")
    Plot_temperature(stopstep,T,T_init)
    Plot_HS_tot(stopstep,hydrostress,T_init)
    Plot_sigma_h(coords,T_init)


def Plot_temperature_dependencies_scatter(E_x_finals,T,n=20,type='Potential',color='g'):
    """Plot the temperature dependencies of some physical parameters, 
        e.g.: potential energy (PE), kinetic energy (KE), total energy (E_t, PE+KE), hydrostatic energy (sigma_h)

    Args:
        E_x_finals (list): A list of any energy (Default: PE) value corresponds to T(T_init_list)
        T (list): Temperature list corresponds to E_x_finals.
        type (str, optional): What kind of data we want to plot versus temperature list. Defaults to 'Potential'.
        color (str, optional): Color of diagram. Defaults to 'g'(green).
    """
    plt.figure(figsize=(8,6))
    plt.scatter([i for i in T], [i for i in E_x_finals],label='Converged '+type)
    plt.plot([i for i in T], [i for i in E_x_finals],c=color,alpha=0.3)
    TITLE = "Temperature dependencies of " +type+'_'+str(n)+'*'+str(n)+" Block"
    plt.title()
    plt.xlabel('Temperature (K)')
    plt.legend()
    save_into_path(TITLE,save_path="./plot/")
    plt.show()


def C_v(E_list,T_list):
    """Relevant specific heat at constant volume.

    Args:
        E_list (list): Energy values corresponding to temperature list.
        T_list (list): Temperature list corresponding to energy values.

    Returns:
        Array: Relevant heat corresponding to temperature list. (len(T_list),)
    """
    C_v_list = []
    for i in range(len(E_list)-1):
        dE = E_list[i+1] - E_list[i]
        dT = T_list[i+1] - T_list[i]
        dE_dT = dE/dT
        T = T_list[i] + dT/2
        C_v_list.append(np.array([dE_dT,T]))
    return np.array(C_v_list).T


if __name__ == "__main__":
    from initialization import *
    c = init_coords(rmin=3.405*(2**(1/6)))
    plot_rdf(c,title='Initial RDF')
    print("hello")