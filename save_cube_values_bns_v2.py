"""
Created on February, 2023

@author: Christos Ntinos

This Module extracts the background noise values of a list of channels for a list of cube names.

Version: 2.0 (January 2026)
"""

import numpy as np
import pandas as pd
import pickle
from tqdm import tqdm
from read_qub import keyword, read_qub
from constants import BBLS_timecorrection


def save_cube_values_bns(cube_name_all, output_name, root_path_raw):
    """
    Parameters
    ----------
    cube_name_all: list of str
        names of the cubes to be used
    root_path_raw: str
        path in repository for raw cubes
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
        columns.append(f'I/F_{band_id+97}_bns')  # start from 97 for IR (numbering starts from 1, 96 if from zero)

    # make dataframe for first cube
    cube_name = cube_name_all[0]
    # read raw cube
    core, suffix, label = read_qub(root_path_raw + f'/v{cube_name}.qub')

    # read exposure from label
    exposure = [float(x) for x in keyword(label, "EXPOSURE_DURATION").strip('()').split(',')][0] * BBLS_timecorrection / 1000  # msec to sec

    df0 = pd.DataFrame({'cube': suffix.shape[1] * suffix.shape[2] * [cube_name],
                        'pixel_index': range(0, suffix.shape[1] * suffix.shape[2]),
                        'exposure': suffix.shape[1] * suffix.shape[2] * [exposure]})

    data_array = np.empty((suffix.shape[1] * suffix.shape[2], 256))

    for band_id in range(0, 256):

        data_array[:, band_id] = suffix[band_id + 96, :, :].flatten()

    df0[columns[3:]] = data_array

    # complete dataframe for all cubes
    for counter in tqdm(range(1, len(cube_name_all)), desc="bns VIMS cubes (-1)"):

        cube_name = cube_name_all[counter]

        core, suffix, label = read_qub(root_path_raw + f'/v{cube_name}.qub')
        exposure = [float(x) for x in keyword(label, "EXPOSURE_DURATION").strip('()').split(',')][0] * BBLS_timecorrection / 1000

        df1 = pd.DataFrame({'cube': suffix.shape[1] * suffix.shape[2] * [cube_name],
                            'pixel_index': range(0, suffix.shape[1] * suffix.shape[2]),
                            'exposure': suffix.shape[1] * suffix.shape[2] * [exposure]})

        data_array = np.empty((suffix.shape[1] * suffix.shape[2], 256))

        for band_id in range(0, 256):
            data_array[:, band_id] = suffix[band_id + 96, :, :].flatten()

        df1[columns[3:]] = data_array
        df0 = pd.concat([df0, df1], ignore_index=True)
        counter = counter + 1

    with open(f'_data/{output_name}.bin', 'wb') as handle:
        pickle.dump(df0, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return df0

