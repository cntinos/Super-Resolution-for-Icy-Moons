"""
Created on October, 2023

@author: Christos Ntinos

This Module creates a map from a list of cubes of the raw data.

Version: 2.0 (January 2026)
(+) raw cubes read with module read_qub.py (binary parser)
"""

import pandas as pd
import pickle
from tqdm import tqdm
from read_qub import read_qub


def save_map_raw(cube_name_all, output_name_raw, root_path_raw):
    """
    Parameters
    ----------
    cube_name_all: list of str
        names of the cubes to be used
    output_name_raw: str
        name for binary output file
    root_path_raw: str
        path in repository for raw cubes

    Returns
    -------
    df0: pandas.Dataframe
        map of raw data (DN) for all bands

    This function calculates the map of the raw data (DN) for all bands. The raw cubes (.qub) must be opened externally
    and provided as csv files in this module.
    """

    # make dataframe for first cube
    cube_name = cube_name_all[0]
    # read raw cube
    data_array, _, _ = read_qub(root_path_raw + f'/v{cube_name}.qub')

    df0 = pd.DataFrame({'cube': data_array.shape[1] * data_array.shape[2] * [cube_name],
                        'pixel_index': range(0, data_array.shape[1] * data_array.shape[2])})

    for i in range(97, 353):
        df0[f'I/F_{i}_raw'] = data_array[i-1, :, :].T.flatten()  # all samples of one line before moving to the next (time equation)

    for counter in tqdm(range(1, len(cube_name_all)), desc="raw VIMS cubes (-1)"):
        cube_name = cube_name_all[counter]
        data_array, _, _ = read_qub(root_path_raw + f'/v{cube_name}.qub')

        df1 = pd.DataFrame({'cube': data_array.shape[1] * data_array.shape[2] * [cube_name],
                            'pixel_index': range(0, data_array.shape[1] * data_array.shape[2])})

        for i in range(97, 353):
            df1[f'I/F_{i}_raw'] = data_array[i-1, :, :].T.flatten()  # all samples of one line before moving to the next (time equation)

        df0 = pd.concat([df0, df1], ignore_index=True)

    with open(f'_data/{output_name_raw}.bin', 'wb') as handle:
        pickle.dump(df0, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return df0

