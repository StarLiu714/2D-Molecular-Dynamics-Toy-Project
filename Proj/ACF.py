import numpy as np 
import matplotlib.pyplot as plt
from tqdm import trange

from update import *
from force_energy import force_arr

def further_loop(coords,v_J,delta_t=5e-15,iter_num=1000):
    """Run further iterations to get the Autocorrelation data

    Args:
        coords (numpy.ndarray): Final coordinate array from mainloop, initial coordinate array of further loop, (400,2)
        v_J (numpy.ndarray): Final velocity array from mainloop, initial velocity array of further loop, (400,2)
        delta_t (float, optional): Differential time quantity, second(s). Defaults to 5e-15.
        iter_num (int, optional): How many iters further carrying out. Defaults to 1000.

    Returns:
        v_log (numpy.ndarray): 3-dimensional velocity array, (iteration_num,num_of_atom,system_dimension)=(1000,20*20,2)
        coords_J_1 (numpy.ndarray): Final coordinate after further loop, (400,2)
    """
    
    # Initialize lists for plotting
    v_log = list()
    # Initialize variable for optimising
    # Initial position(400,2), vector(400,2) and systemforce(1,2)
    coords_J_1 = coords.copy()
    v_J_1 = v_J.copy()
    force_J_1 = force_arr(coords_J_1)

    for iter in trange(iter_num+1):

        coords_J = coords_J_1.copy()  # positions of 400 atoms at JΔt
        force_J = force_J_1.copy()  # forces of of 400 atoms at JΔt
        v_J = v_J_1  # velocities of of 400 atoms at JΔt

        # Update the position array, coords_J.shape=(400,2), force_J.shape=(1,2)
        # According to r_i(J+1) = r_i(J) + (dt/2m)*(F_i(J+1)+F_i(J))
        coords_J_1 = update_position_vectorize(coords_J,force_J,v_J,delta_t,m=6.63*1e-26)
        # Update the force array, force_J.shape=(1,2)
        force_J_1 = force_arr(coords_J_1)
        # Update the velocity array, v.shape=(400,2), v_i.shape=(1,2)
        # According to v_i(J+1) = v_i(J) + (dt/2m)*(F_i(J+1)+F_i(J))
        v_J_1 = update_velocity_vectorize(force_J,force_J_1,v_J,delta_t,m=6.63*1e-26)

        v_log.append(v_J_1.copy())  # v_J_1: (400,2), v_log: (1000,400,2)

    return np.array(v_log), coords_J_1


def ACF(v_J_versus_t):
    """Autocorrelation function. Calculating by using the further iterations after main loop.

    Args:
        v_J_versus_t (numpy.ndarray): Final 3-dimensional velocity array after further-loop, (iteration_num,num_of_atom,system_dimension)=(1000,20*20,2)

    Returns:
        tlist: continuous number of timestep, i.e.:1,2,3,...(int(iter_num/2))
        ac: autocorrelation value. i.e.:<v(0)v(t)>
    """

    # Figure out <v(t)>s, i.e. <v(0)>~<v(1000)>
    # v_t_i.shape=(1000,400)
    v_t_i = np.sqrt(np.sum(v_J_versus_t**2,axis=2))
    ttt = int(len(v_J_versus_t)/2)  # ttt=500
    M = len(v_J_versus_t)-ttt
    tlist,ac=[],[]
    for t in trange(1,ttt+1):
        # <v(0)v(t)> =c =1/M*[v(1)v(1+t)+v(2)v(2+t)+...+v(M)v(M+t)]
        c = 0
        for j in range(M):
            # c.shape=(400,)
            c += (v_t_i[j] * v_t_i[j+t])  # c_j= v(1)v(1+t),v(2)v(2+t),...,v(M)v(M+t)
        c = np.sum(c)/(M*400)  # c is scalar
        # ac =[<v(0)v(t)>,<v(1)v(t+1)>,<v(2)v(t+2)>,...,<v(1000-t)v(1000)>]
        ac.append(c)
        tlist.append(t)
    return tlist,ac


def collect_data_for_acf(T,datapack,dt=5e-15,iter_num=1000):
    """Choose the temperature we want to collect acf data.

    Args:
        T (int): Temperature we want to collect acf data.
        datapack (dict): a zipped datapack recording all the simulation results, for each temperature, the key is str(temperature)+'K'
        dt (float, optional): Differential value of time, second. Defaults to 5e-15.
        iter_num (int, optional): Number of iteration. Defaults to 1000.

    Returns:
        v_J_versus_t(numpy.ndarray): total velocity of the system for each time interval (dt).
        coords(numpy.ndarray): Final locations of each atom, (400,2)
    """
    # Collect coords&v_J arrays at every timestep
    # Read already-relaxed atom positions (400,2) and latest v_J (400,2)
    coords, _, _, _, v_J, _, _ = datapack[str(T)+'K']

    # Collect data from further loop, v_J_versus_t: (1000,400,2)
    v_J_versus_t,coords = further_loop(coords,v_J,dt,iter_num)
    
    return v_J_versus_t,coords


def plot_acf(acf_result,T,dt=1e-14):
    """Plot ACF versus a series of temperature.

    Args:
        acf_result (dict): A dictionary recording ACF results. e.g.: acf_result['10K']=(tlist,ac) 
        T (list): Temperature list, as the T_desire of certain ACF-iteration plot.
        dt (float, optional): Timestep, second(s). Defaults to 1e-14.
    """
    
    tlist,ac = acf_result
    # Plot the temperature dependence
    plt.figure(figsize=(8,6))
    plt.plot(tlist,ac,alpha=0.6)
    plt.title('Autocorrelation Function at '+str(T)+'K')
    plt.xlabel('Time ('+str(dt)+' second)')
    plt.ylabel('Velocity correlation function')
    plt.show()


def D_T(T_list,ACF_results,dt=1e-14):
    """Function of D(T)

    Args:
        T_list (int list): Series of ambient (desired) temperature, kelvin.
        acf_results (dict): recorded the acf results for each temperature, 
            i.e.: for 10K, ac_list = acf_results['10K']
        dt (float, optional): timestep. Defaults to 1e-14.

    Returns:
        list: self-diffusion coefficient list corresponds to the T_list
    """
    D = []
    for i in range(len(T_list)):
        _,ac_list = ACF_results[str(T_list[i])+'K']
        ac = np.array(ac_list)
        D.append(1/3*np.sum(ac*dt))
    return D


def plot_D_T(T_list,D_list):
    """Plot the relationship between log{self-diffusion coefficient (D)} and reciprocal{ambient temperature(T)}
        i.e.: log(D) versus (1/T)

    Args:
        T_list (int list): Series of ambient (desired) temperature, kelvin
        D_list (float list): Corresponding self-diffusion coefficient value
    """
    lnD = np.log(D_list[1:])
    T_reciprocal = 1/np.array(T_list[1:])
    plt.scatter(T_reciprocal,lnD)
    plt.plot(T_reciprocal,lnD,alpha=0.6)
    plt.xlabel('T⁻¹ (1/K)')
    plt.ylabel('ln(D)')
    plt.title('Diffusion Coefficient Function')
    plt.show()

