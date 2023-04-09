import numpy as np

def save_datapack(T_init,datapack):
    """Save simulation result data as .npy files.

    Args:
        T_init (int): The certain temperature we want to save its results.
        datapack (dict): a zipped datapack recording all the simulation results, for each temperature, the key is str(temperature)+'K'
    """
    np.save('MD_'+str(T_init)+'K_datapack',datapack[str(T_init)+'K'])