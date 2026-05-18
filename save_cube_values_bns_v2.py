'''
Created on February, 2023

@author: Christos Ntinos

This Module extracts the background noise values of a list of channels for a list of cube names.

Version: 2.0 (January 2026)

'''

from pyvims import VIMS
import numpy as np
import pandas as pd
import pickle
from tqdm import tqdm


def save_cube_values_bns(cube_name_all, output_name, root_path):
    """
    Parameters
    ----------
    cube_name_all: list of str
        names of the cubes to be used
    root_path: str
        path in repository to save cubes
    output_name: str
        name for binary output file

    Returns
    -------
    df0: pandas.Dataframe
        map of background noise (bns) for all bands

    This function calculates the map of the background noise for all bands.
    """
    columns = ['cube', 'pixel_index', 'exposure']
    for band_id in range(0, 256):
        columns.append(f'I/F_{band_id+97}_bns')  # start from 97 for IR

    # make dataframe for first cube
    cube_name = cube_name_all[0]
    cube = VIMS(cube_name, root=root_path)

    df0 = pd.DataFrame({'cube': cube.NL * cube.NS * [cube_name],
                        'pixel_index': range(0, cube.NL * cube.NS),
                        'exposure': cube.NL * cube.NS * [cube.expo]})

    data_array = np.empty((cube.NL * cube.NS, 256))

    for band_id in range(0, 256):

        data_array[:, band_id] = cube.background[band_id].flatten()

    df0[columns[3:]] = data_array

    counter = 1
    for counter in tqdm(range(1, len(cube_name_all)), desc="bns VIMS cubes (-1)"):

        cube_name = cube_name_all[counter]

        cube = VIMS(cube_name, root=root_path)

        df1 = pd.DataFrame({'cube': cube.NL * cube.NS * [cube_name],
                            'pixel_index': range(0, cube.NL * cube.NS),
                            'exposure': cube.NL * cube.NS * [cube.expo]})

        data_array = np.empty((cube.NL * cube.NS, 256))

        for band_id in range(0, 256):
            data_array[:, band_id] = cube.background[band_id].flatten()

        df1[columns[3:]] = data_array
        df0 = pd.concat([df0, df1], ignore_index=True)
        counter = counter + 1

    with open(f'_data/{output_name}.bin', 'wb') as handle:
        pickle.dump(df0, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return df0

