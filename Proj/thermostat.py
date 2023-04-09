import numpy as np 
from numba import jit
from force_energy import E_k

@jit
def temperature(v_J, N=400, K_B=8.67*1e-5):
    """Calculate the temperature of a system of atoms.

    Args:
        v_J (numpy.ndarray): A 2D array containing the velocities of atoms in the system. <---------------------
        N (int, optional): The total number of atoms in the system. Defaults to 400.
        K_B (float, optional): Boltzmann's constant in eV/K. Defaults to 8.67*1e-5.

    Returns:
        float: The temperature of the system.
    """
    E_kin = E_k(v_J)
    return E_kin/(N*K_B)


@jit
def kappa(T_wanted,T_now):
    """Scaling parameter kappa

    Args:
        T_wanted (float): The desired ambient temperature we set, scalar
        T_now (float): The real time temperature, scalar

    Returns:
        float: Rescaling parameter kappa, scalar
    """
    return (T_wanted/T_now)**0.5


@jit
def v_new(v_J,kk):
    """To carry out velocity scaling

    Args:
        v_J (numpy.ndarray): current velocity of 400 atoms, (400,2)
        kk (float): current kappa value

    Returns:
        array: the already rescaled velocity array, (400,2)
    """
    return v_J*kk

if __name__ == "__main__":
    print("hello")