from numba import jit

@jit
def update_velocity_vectorize(force_J,force_J_1,v_J,delta_t,m=6.63*1e-26):
    """The update function for each atom's velocity
    Formula:
    v_i(J+1) = v_i(J) + (dt/2m)*(F_i(J+1)+F_i(J))

    Args:
        force_J (numpy.ndarray): Force array of previous timestep, eV/angstrom, (400,2)
        force_J_1 (numpy.ndarray): Force array of current timestep, eV/angstrom, (400,2)
        v_J (numpy.ndarray): Velocity array of previous timestep, (400,2)
        delta_t (float): Timestep, second(s).
        m (float, optional): Atomic mass, kg. Defaults to 6.63*1e-26.

    Returns:
        numpy.ndarray: Atomic velocity of next timestep, angstrom/s, (400,2) 
    """
    delta_v = delta_t/(2*m) * (force_J + force_J_1)
    v_J += delta_v*16.02   # eV*s/(Å*kg) to Å/s, 1eV*s/(Å*kg) =16.02Å/s

    return v_J

@jit
def update_position_vectorize(coords_J,force_J,v_J,delta_t,m=6.63*1e-26):
    """The update function for each atom's coordinate

    Args:
        coords_J (numpy.ndarray): Atomic positions of current timestep, angstrom, (400,2)
        force_J (numpy.ndarray): Force array of previous timestep, eV/angstrom, (400,2)
        v_J (numpy.ndarray): Velocity array of previous timestep, angstrom/s, angstrom/s, (400,2)
        delta_t (float): timestep, second(s)
        m (float, optional): Atomic mass, kg. Defaults to 6.63*1e-26.

    Returns:
        numpy.ndarray: Atomic positions of next timestep, angstrom, (400,2) 
    """
    # eV*s²/(Å*kg) to Å, 1eV*s²/(Å*kg) =16.02Å
    delta_r = v_J *delta_t + (delta_t**2/(2*m))*force_J*16.02
    coords_J_1 = coords_J + delta_r
    return coords_J_1