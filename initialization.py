import numpy as np 

# Create intial 20*20 blocks (with random displacement)
def init_coords(n=20,rmin=3.405*(2**(1/6))):
    """Initialize a 2D array of coordinates with random displacement.

    This function generates a 2D array of coordinates on a regular grid, and adds a small random
    displacement to each coordinate. The grid is created with a specified minimum distance between
    the coordinates (rmin) and a specified number of points in each dimension (n).

    Args:
        rmin (float, optional): Minimum distance between coordinates. Defaults to 3.405*(2**(1/6)).
        n (int, optional): Number of points in each dimension of the grid. Defaults to 20.

    Returns:
        numpy.ndarray: A 2D array of shape (n * n, 2) containing the coordinates with random displacement.
    """
    # 20*20 regular array
    coords_init = np.array([[rmin*float(i), rmin*float(j)] for i in range (0,n) for j in range(0,n)])
    # random displacement, range as (-0.5,0.5)*1e-2
    coords_init += (np.random.random(np.shape(coords_init))-0.5)*1e-2
    return coords_init


if __name__ == "__main__":
    # Let's create an initial regular 20*20 block.
    coords = init_coords(n=20)
    print(coords)