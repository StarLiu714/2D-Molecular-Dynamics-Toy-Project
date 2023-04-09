from savedata import save_datapack
from initialization import *
from force_energy import *
from thermostat import *
from plotting import *
from update import *
from ACF import *

import numpy as np 
import matplotlib.pyplot as plt
from tqdm import trange


def main_loop(coords, iter_num=100000, delta_t=1e-14, m=6.63*1e-26, K_B=8.67*1e-5, E_kin=0, T_init=0, tolerance_kappa=2e-4, tolerance_T=2):
    """The mainloop follows Verlet Algorithm.

    Args:
        coords (numpy.ndarray): input(initial) atoms' positions, angstrom, (400,2)
        iter_num (int, optional): Number of iterations. Defaults to 100000.
        delta_t (float, optional): Time step, second(s). Defaults to 1e-14.
        m (float, optional): Mass of Argon atom, kg. Defaults to 6.63*1e-26.
        K_B (float, optional): Boltzmann Constant, eV/K. Defaults to 8.67*1e-5.
        E_kin (int, optional): System kinetic energy, eV. Defaults to 0.
        T_init (int, optional): Desired Temperature, K. Defaults to 0.
        tolerance_T (float, optional): Temperature threshold triggering rescaling,Kelvin(K). Defaults to 1.5.

    Returns:
        coords_J_1 (numpy.ndarray): Output(resulting) atoms' positions, angstrom, (400,2)
        E_p_verlet (list): Potential energy recorded every 50 timesteps, eV
        E_k_verlet (list): Kinetic energy recorded every 50 timesteps, eV
        hydrostress_verlet (list): Total hydrostatic stree of 400 atoms system, eV/(angstrom^2)
        v_J (numpy.ndarray): Output(resulting) velocities of of 400 atoms, angstrom/(s^2), (400,2)
        T (list): Real-time temperature recorded every 200 timesteps, kelvin
        iter (int): The resulting iteration number, scalar
    """
    
    atom_num=len(coords)
    T_initial = T_init
    # Initialize lists for plotting
    E_p_verlet, E_k_verlet, hydrostress_verlet, T =[],[],[],[]
    # Initialize variable for optimising
    # Initial velocity array, all zeros
    v = np.zeros(np.shape(coords))  
    # Initial position(400,2), vector(400,2) and systemforce(1,2)
    coords_J_1 = coords.copy()
    v_J_1 = v.copy()
    force_J_1 = force_arr(coords_J_1)

    for iter in trange(iter_num+1):

        coords_J = coords_J_1.copy()  # positions of 400 atoms at JΔt
        force_J = force_J_1.copy()  # forces of of 400 atoms at JΔt, eV/angstrom
        v_J = v_J_1  # velocities of of 400 atoms at JΔt
        # Plot intermedium relaxed structures
        if iter in [0,1000,2000,4000,8000,20000,50000,75000,100000,130000,200000]:
            plot_block(coords_J,iter,T_initial)
        # Update the position array, coords_J.shape=(400,2), force_J.shape=(1,2)
        # According to r_i(J+1) = r_i(J) + (dt/2m)*(F_i(J+1)+F_i(J))
        coords_J_1 = update_position_vectorize(coords_J,force_J,v_J,delta_t,m=6.63*1e-26)
        # Update the force array, force_J.shape=(1,2)
        force_J_1 = force_arr(coords_J_1)
        # Update the velocity array, v.shape=(400,2), v_i.shape=(1,2)
        # According to v_i(J+1) = v_i(J) + (dt/2m)*(F_i(J+1)+F_i(J))
        v_J_1 = update_velocity_vectorize(force_J,force_J_1,v_J,delta_t,m=6.63*1e-26)
        
        # Rescale the velocity
        # Check convergency on temperature
        if (iter>=5000) and iter%5000 == 0:
            # If kappa is around 1, stop the main loop
            T_now = temperature(v_J_1)
            T_previous = temperature(v_J)
            if T_init == 0:
                kk=0
            else:
                kk=kappa(T_initial,T_now)
            if np.absolute(kk-1) <= tolerance_kappa:
                return coords_J_1, E_p_verlet, E_k_verlet, hydrostress_verlet, v_J, T, iter
            # Check whether rescale is needed
            delta_T = np.absolute(T_now-T_previous)
            if delta_T <= tolerance_T:
                v_J_1 = v_new(v_J_1,kk)

        #Calculate Energy, hydrostatic pressure
        #You can do this at every tens or hundreds of iterations
        if (iter%50==0):
            E_p_i,E_k_i = 0,0
            E_p_i += E_p(coords_J_1)
            E_k_i += E_k(v_J_1)
            E_p_verlet.append(E_p_i)
            E_k_verlet.append(E_k_i)
        if (iter%200==0):
            hydrostress_verlet.append(total_hydrostatic_stress(coords_J_1))
            T.append(temperature(v_J))
    
    return coords_J_1, E_p_verlet, E_k_verlet, hydrostress_verlet, v_J, T, iter


def running_MD_for_different_temperatures(T_list, coords_init, iter_num, datapack_init={}, dt=5e-15, tolerance_T_run=1.5):
    """Auto running program, should given the wanted temperatures as Temperature List

    Args:
        T_list (list): the series of temperatures we want to run for MD, Kelvin
        coords_init (numpy.ndarray): initial coordinate array, (400,2), angstrom
        iters (int, optional): number of iteration
        datapack_init(dict): an empty dictionary
        dt (float, optional): time step, second. Defaults to 5e-15.
        tolerance_T_run (float, optional): temperature threshold triggering rescaling. Defaults to 1.5. 

    Returns:
        dict: a zipped datapack recording all the simulation results, for each temperature, the key is str(temperature)+'K',
            e.g.: datapack['10K'] = coords_10K, E_p_10K, E_k_10K, hydrostress_10K, v_J_10K, T_10K, stopstep_10K
    """

    len_T = len(T_list)
    n=int(len(coords_init)**0.5)
    datapack = datapack_init
    for i in range(len_T):
        T_initial = T_list[i]
        iters = iter_num
        datapack[str(T_list[i])+'K'] = np.array(main_loop(coords_init,iter_num=iters,delta_t=dt,T_init=T_initial,tolerance_T=tolerance_T_run))
        np.save('MD_'+str(n)+'Block_'+str(T_initial)+'K_datapack',datapack[str(T_initial)+'K'])
    return datapack



"""Constant
"""
K_B   = 8.617e-5   #eV/K
K_B_J = 1.38e-23   #J/K
m     = 6.63*1e-26 #kg

"""Initialize atomic positions for MD
"""
coords = init_coords(n=20)
coords_MD = coords.copy()
plot_block(coords_MD,0)

"""Initialize atomic velocities as 0
"""
v = np.zeros(np.shape(coords))  # initialize velocity array

"""Initialize datapack
"""
datapack = {}

"""Initialize Temperature List
"""
T_list_1 = [
    5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,100
]

"""Run MD at desire temperature list (5~100K)
"""
datapack_1 = running_MD_for_different_temperatures(
    T_list=T_list_1, 
    coords_init=coords_MD, 
    iter_num=150000, 
    dt=5e-15)

"""Save datapacks
"""
for i in range(len(T_list_1)):
    save_datapack(T_list_1[i],datapack_1)

"""Plot and output all the diagrams
"""
for i in range(len(T_list_1)):
    coords_1, E_p, E_k, hydrostress, v_J, T, stopstep = datapack_1[str(T_list_1[i])+'K']
    Plot_results(coords_1, E_p, E_k, hydrostress, T, T_list_1[i], stopstep)

"""Evaluate the total energy and total hydrostatic pressure and plot their temperature dependencies
"""
T_init_list = T_list_1
E_p_finals,E_k_finals,HS_finals,E_tot_finals,E_var_T=[],[],[],[],[]
for i in range(len(T_init_list)):
    T = T_init_list[i]
    E_p,E_k,_,HS = datapack[str(T)+'K'][1:5]
    E_tot = np.array(E_p[-1000:])+np.array(E_k[-1000:])
    E_var = np.var(E_tot)
    E_p_T,E_k_T,HS_T = E_p[-1],E_k[-1],HS[-1]
    E_tot_T = E_p_T+E_k_T
    E_p_finals.append(E_p_T)
    HS_finals.append(HS_T)
    E_tot_finals.append(E_tot_T)
    E_var_T.append(E_var)

"""Relavent Heat
"""
n1=20
Plot_temperature_dependencies_scatter(E_p_finals,T_init_list,n1)
Plot_temperature_dependencies_scatter(E_tot_finals,T_init_list,n1,type='Total Energy',color='r')
Plot_temperature_dependencies_scatter(HS_finals,T_init_list,n1,'Hydrostatic Pressure',color='b')
C_v_array = C_v(E_tot_finals,T_init_list)
c1,c2=C_v_array
Plot_temperature_dependencies_scatter(c1,c2,n=20,type='Potential',color='g')

"""Autocorrelation
"""
v_log_coords_dict = {}
for i in range(len(T_init_list)):
    v_log_coords_dict[str(T_init_list[i])+'K'] = collect_data_for_acf(T_init_list[i],datapack)
acf_results = {}
for i in range(len(T_init_list)):
    v_J_versus_t,coords = v_log_coords_dict[str(T_init_list[i])+'K']
    acf_results[str(T_init_list[i])+'K'] = ACF(v_J_versus_t)
dt = 1e-14
for i in range(len(T_init_list)):
    plot_acf(acf_results[str(T_init_list[i])+'K'],T_init_list[i],dt)

"""Self-diffusion coefficient
"""
D_list = D_T(T_init_list,acf_results)
# Discard the first 5 temperatures' data
# <25K the 2D argon is definitely solid with no linear relationship
plot_D_T(T_init_list[5:],D_list[5:])


if __name__ == "__main__":
    print("hello")
