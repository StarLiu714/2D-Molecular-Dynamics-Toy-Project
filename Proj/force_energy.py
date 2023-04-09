import numpy as np 
from numba import jit

@jit 
def dist(a, b):
    """Calculate the Euclidean distance between two points. 

    Args:
        a (numpy.ndarray): An array representing the first point (x1, y1).
        b (numpy.ndarray): An array representing the second point (x2, y2).

    Returns:
        numpy.ndarray: The Euclidean distance between the two points.
    """
    return ((a[0] - b[0])**2 + (a[1] - b[1])**2)**0.5


@jit
def LJ(r, A=-6.8102128e-3, B=-5.5640876e-3, R0=7.0, RCUT=7.5, EPSILON=0.010323, SIGMA = 3.405):  # r = distance from central atom
    """Calculate the LJ potential for a given distance.

    Args:
        r (numpy.ndarray): The distance between two atoms.
        A (float, optional): Constant parameter A of LJ potential function when (r_0 >r >rcut). Defaults to -6.8102128*10**-3.
        B (float, optional): Constant parameter B of LJ potential function when (r_0 >r >rcut). Defaults to -5.5640876*10**-3.
        R0 (float, optional): Distance threshold for the potential function. Defaults to 7.0.
        RCUT (float, optional): Cutoff distance for the potential function. Defaults to 7.5.
        EPSILON (float, optional): Depth of the potential well. Defaults to 0.010323.
        SIGMA (float, optional): Distance at which the potential is zero. Defaults to 3.405. 

    Returns:
        numpy.ndarray: The LJ potential for the given distance.
    """
    if r <= R0:
        return 4*EPSILON*( (SIGMA/r)**12 - (SIGMA/r)**6 )
    elif r<=RCUT:
        return A*(r-RCUT)**3+B*(r-RCUT)**2
    return 0


#we make a simple derivative of LJ 
@jit
def dLJ(r, A=-6.8102128e-3, B=-5.5640876e-3, R0=7.0, RCUT=7.5, EPSILON=0.010323, SIGMA = 3.405):  
    """Calculate the derivative of the LJ potential for a given distance.

    Args:
        r (numpy.ndarray): The distance between two atoms.
        A (float, optional): Constant parameter A of LJ potential function when (r_0 >r >rcut). Defaults to -6.8102128*10**-3.
        B (float, optional): Constant parameter B of LJ potential function when (r_0 >r >rcut). Defaults to -5.5640876*10**-3.
        R0 (float, optional): Distance threshold for the potential function. Defaults to 7.0.
        RCUT (float, optional): Cutoff distance for the potential function. Defaults to 7.5.
        EPSILON (float, optional): Depth of the potential well. Defaults to 0.010323.
        SIGMA (float, optional): Distance at which the potential is zero. Defaults to 3.405. 

    Returns:
        numpy.ndarray: The derivative of the LJ potential for the given distance.
    """
    if r <= R0:
        return 4*EPSILON*( -12.0/r *(SIGMA/r)**12 + 6.0/r *(SIGMA/r)**6 )
    if r <= RCUT:
        return 3*A*(r-RCUT)**2+2*B*(r-RCUT)
    return 0

@jit
def force(target,coords,RCUT=7.5,n=400):  
    """Computes the force acting on a target atom in a 2D system of atoms.

    The force is calculated by summing the LJ potential forces 
    between the target atom and all other atoms within the specified cutoff distance.

    Args:
        target (int): Index of the target atom, i.e., coords[target] would give the coordinates (x, y).
        coords (numpy.ndarray): ndarray of all atoms' positions, shape: (n, 2).
        RCUT (float, optional): Cutoff distance for the LJ potential calculation. Defaults to 7.5.
        n (int, optional): Total number of atoms in the system. Defaults to 400.

    Returns:
        numpy.ndarray: The net force acting on the target atom, shape: (2,), where the first element is the x-component and the second element is the y-component.
    """
    F = np.zeros(2)  # force_x,force_y = 0,0

    for i in range(n):  # sum forces over neighbors within potential cutoff
        d = dist(coords[target], coords[i])
        if i != target and d <= RCUT:
            r = coords[i] - coords[target]  # r[0],r[1] = x,y
            F += dLJ(d) / d * r

    return F



### get the central atom stress components
@jit
def stress(center, coords, n=20, rmin=3.405*(2**(1/6)), rcut=7.5):  
    """Calculate the stress components of a central atom in a system of atoms.

    Args:
        center (int): The index of the central atom in the coords array.
        coords (numpy.ndarray): A 2D array containing the atomic coordinates.
        n (int, optional): The number of atoms in each column. Defaults to 20.
        rmin (float, optional): Minimum distance between atoms. Defaults to 3.405*(2**(1/6)).<---------------------
        rcut (float, optional): Cutoff distance for stress calculation. Defaults to 7.5.

    Returns:
        numpy.ndarray: A 1D array containing the stress components (xx, xy, yx, yy) of the central atom.
    """
    AREA = (n-1)*3.405*rmin**2
    atom_num = n**2
    
    # stress xx,xy,yx,yy are s[0,0] [0,1] [1,0] [1,1] respectively
    s = np.zeros(4)

    for i in range(atom_num):
        d = dist(coords[i],coords[center])
        if i!=center and d<=rcut:

            #clarify what are these steps in relation to the formulae
            x = coords[i][0] - coords[center][0]
            y = coords[i][1] - coords[center][1]
            
            delta_matrix = np.array([x**2,x*y,y*x,y**2])
            s += delta_matrix *dLJ(d)/d

    s /= (AREA/atom_num)
    return s

@jit
def stress_total(coords,n=20):
    """Calculate the total stress components of a system of atoms.

    Args:
        coords (numpy.ndarray): A 2D array containing the atomic coordinates.
        n (int, optional): The number of atoms in each column. Defaults to 20.

    Returns:
        numpy.ndarray: A 1D array containing the total stress components (xx, xy, yx, yy) of the system.
    """
    atom_num = n**2
    stress_tot = np.array([0.,0.,0.,0.])
    
    # sum up 
    for i in range(atom_num):
        stress_tot += stress(i,coords)
    return stress_tot

# Hydrostatic stress for each single atom
@jit
def stress_h(center,coords,n=20,rmin=3.405*(2**(1/6)),rcut=7.5):
    """Calculate the hydrostatic stress of a single atom in a system of atoms.

    Args:
        center (int): The index of the center atom in the coords array.
        coords (numpy.ndarray): A 2D array containing the atomic coordinates.
        n (int, optional): The number of atoms in each column. Defaults to 20.
        rmin (float, optional): Minimum distance between atoms, essentially it is the distance of one atom's nearest neighbour. Defaults to 3.405*(2**(1/6)).<---------------------
        rcut (float, optional): Cutoff distance for stress calculation. Defaults to 7.5.

    Returns:
        float: The hydrostatic stress of the single atom.
    """
    AREA = ((n-1)*rmin)**2
    atom_num = n**2

    # stress xx,yy are s[0],s[1] respectively
    s = np.zeros(2)

    for i in range(atom_num):
        d = dist(coords[i],coords[center])
        if i!=center and d<=rcut:
            r = coords[i] - coords[center]
            delta_matrix = np.array([r[0]**2,r[1]**2])
            s += delta_matrix *dLJ(d)/d

    s /= (AREA/atom_num)

    stress_hs = - 0.5 * sum(s) # only the hydrostatic pressure is needed
    
    return stress_hs

# Total Hydrostatic stress for 20*20 block
@jit
def total_hydrostatic_stress(coords,n=20):
    """Calculate the total hydrostatic stress for a block of atoms.

    Args:
        coords (numpy.ndarray): A 2D array containing the x and y coordinates of atoms in the block.
        n (int, optional): The number of atoms for each column. Defaults to 20.

    Returns:
        numpy.ndarray: The total hydrostatic stress for the block.
    """
    sigma_h_tot = 0
    atom_num = n**2
    for i in range(atom_num):
        sigma_h_tot += stress_h(i,coords)
    return sigma_h_tot

@jit
def E_k(v_J, m=6.63*1e-26):
    """Calculate the kinetic energy of a system of atoms.

    Args:
        v_J (numpy.ndarray): Velocities of atoms at a certain timestep, angstrom/s, (400,2)
        m (float, optional): The mass of one single Argon atom. Defaults to 6.63*1e-26, kg

    Returns:
        float: The total kinetic energy of the system.
    """
    vv = v_J.flatten()
    # 1Å²*kg/s² =0.0624eV, 0.5*m*sum(vv**2) is of Å²*kg which should be turn into eV
    return 0.5*m*sum(vv**2)*0.0624

@jit
def E_p(coords,n=400):
    """Calculate the potential energy of the system for a given set of atomic coordinates.

    Args:
        coords (numpy.ndarray): A 2D array containing the atomic coordinates.
        n (int, optional): Number of atoms in the system. Defaults to 400.

    Returns:
        float: The potential energy of the system.
    """
    # initialze energy
    potential = 0
    for i in range(n):
        for j in range(i+1,n):
            d = dist(coords[i],coords[j])
            # equilibrium thus E(xy)=E(yx), and E(xx) and E(yy) did not be count on for now
            potential += LJ(d) 
    return potential

@jit
def force_arr(coords):
    """Calculate the force acting on each atom in the system.

    Args:
        coords (numpy.ndarray): A 2D array containing the x and y coordinates of atoms in the system.

    Returns:
        numpy.ndarray: A 2D array containing the force acting on each atom, eV/angstrom
    """
    f = np.zeros(np.shape(coords))
    for a in range(len(f)):
        f[a] = force(a,coords)
    return f

@jit
def sigma_h(coords):
    """Calculate the hydrostatic stress for each atom in the system.

    Args:
        coords (numpy.ndarray): A 2D array containing the x and y coordinates of atoms in the system.

    Returns:
        numpy.ndarray: A 1D array containing the hydrostatic stress values for each atom.
    """
    hs = np.zeros(len(coords))
    for i in range(400):
        hs[i] = stress_h(i,coords)
    return hs

if __name__ == "__main__":
    from initialization import *
    coords = init_coords(n=20)
    print(total_hydrostatic_stress(coords))