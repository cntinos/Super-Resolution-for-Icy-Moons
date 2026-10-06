'''
Created on May, 2022

@author: Christos Ntinos

This module creates and saves a map combining navigation and spectral information. List of cubes must be the same as the
one used in save_cube_values.py

Version: 2.0 (March 2026)
'''

import pickle
import pandas as pd
import read_cube_header
from planetaryimage import CubeFile
import numpy as np
from tqdm import tqdm


def save_map(cube_name_all, cube_flybys, num_points, method, method_sampling, root_path, output_name):
    """
    Parameters
    ----------
    cube_name_all: list of str
        names of the cubes to be used
    cube_flybys: list of str
        names of flybys corresponding to each cube to be used
    num_points: int
        number of points per edge of the polygon - pixel
    method: str
        "ELLIPSOID" or "DSK" for modeling the target body
    method_sampling: str
        "DISTANCE" or "DIAGONAL" for calculating pixel sampling
    root_path: str
        path in repository for saved calibrated cubes
    output_name: str
        name for binary output file

    Returns
    -------
    df0: pandas.Dataframe
        map of radiance factor for all bands
    polygon_coordinates: list of lists
        polygon_coordinates = [longitudes, latitudes, nav_cube_all_25, nav_cube_all_26]

    This function calculates the map for all bands.
    """

    # initialise outputs
    longitudes = []
    latitudes = []
    nav_cube_all_25 = []
    nav_cube_all_26 = []

    # initiate dataframe for first cube
    cube_name = cube_name_all[0]

    # load info from cube header
    r_size, c_size, target_cube, ntvtim_start, ntvtim_stop, mirtim, \
    R_offset, C_offset, modeid_IR, inttim = read_cube_header.main(cube_name, root_path)

    # assign time
    cube_year = ntvtim_start.split('-')[0]
    cube_day = ntvtim_start.split('-')[1].split('T')[0]
    if int(cube_day) < 182.5:
        cube_time = cube_year+'.0'
    else:
        cube_time = cube_year+'.5'  # ignore leap years

    # load navigation cube (deserialize)
    with open(f'_nav_cubes/num_points_{num_points}_{method[0:3]}_{method_sampling}/V{cube_name}.nav', 'rb') as handle:
        nav_cube, longitude_all, latitude_all = pickle.load(handle)

    pixel_cube_name = [cube_name] * nav_cube.shape[1]
    pixel_index = range(0, nav_cube.shape[1])
    pixel_cube_flybys = [cube_flybys[cube_name_all.index(cube_name)]] * nav_cube.shape[1]
    pixel_cube_time = [cube_time] * nav_cube.shape[1]

    longitudes.extend(longitude_all)
    latitudes.extend(latitude_all)

    sampling_all = nav_cube[31]  # pixel sampling
    nav_cube_all_24 = nav_cube[24]  # reflectance correction
    nav_cube_all_25.extend(nav_cube[25])  # N.P. in
    nav_cube_all_26.extend(nav_cube[26])  # S.P. in
    nav_cube_all_28 = nav_cube[28]  # incidence angle
    nav_cube_all_29 = nav_cube[29]  # emission angle
    nav_cube_all_30 = nav_cube[30]  # phase angle
    nav_cube_all_32 = nav_cube[32]  # azimuthal angle
    nav_cube_all_6 = nav_cube[6]
    nav_cube_all_7 = nav_cube[7]
    nav_cube_all_8 = nav_cube[8]
    nav_cube_all_9 = nav_cube[9]
    nav_cube_all_10 = nav_cube[10]
    nav_cube_all_15 = nav_cube[15]
    nav_cube_all_16 = nav_cube[16]
    nav_cube_all_17 = nav_cube[17]
    nav_cube_all_18 = nav_cube[18]
    nav_cube_all_19 = nav_cube[19]

    # read spectral info from cubes with planetaryimage module
    isisimage = CubeFile.open(root_path + '/C' + cube_name + '_ir.cub')
    cube_spectra = isisimage.data

    # create navigation cube dataframe to be pickled
    my_data = {'cube': pixel_cube_name,
               'pixel_index': pixel_index,
               'flyby': pixel_cube_flybys,
               'time': pixel_cube_time,
               'sampling': sampling_all,
               'refl_cor': nav_cube_all_24,
               'incidence': nav_cube_all_28,
               'emission': nav_cube_all_29,
               'phase': nav_cube_all_30,
               'azimuthal': nav_cube_all_32,
               'lon_0': nav_cube_all_6,
               'lon_1': nav_cube_all_7,
               'lon_2': nav_cube_all_8,
               'lon_3': nav_cube_all_9,
               'lon_4': nav_cube_all_10,
               'lat_0': nav_cube_all_15,
               'lat_1': nav_cube_all_16,
               'lat_2': nav_cube_all_17,
               'lat_3': nav_cube_all_18,
               'lat_4': nav_cube_all_19}

    df0 = pd.DataFrame(my_data)

    for band_id in range(97, 353):
        values = cube_spectra[band_id - 97, :].T.flatten()
        # replace -3.4028227e+38 with np.nan
        values[values < 0] = np.nan
        values[values > 1] = np.nan
        df0[f'I/F_{band_id}'] = values

    for counter in tqdm(range(1, len(cube_name_all)), desc="VIMS cubes"):

        cube_name = cube_name_all[counter]

        # load info from cube header
        r_size, c_size, target_cube, ntvtim_start, ntvtim_stop, mirtim, \
        R_offset, C_offset, modeid_IR, inttim = read_cube_header.main(cube_name, root_path)

        # assign time
        cube_year = ntvtim_start.split('-')[0]
        cube_day = ntvtim_start.split('-')[1].split('T')[0]
        if int(cube_day) < 182.5:
            cube_time = cube_year + '.0'
        else:
            cube_time = cube_year + '.5'  # ignore leap years

        # load navigation cube (deserialize)
        with open(f'_nav_cubes/num_points_{num_points}_{method[0:3]}_{method_sampling}/V{cube_name}.nav', 'rb') as handle:
            nav_cube, longitude_all, latitude_all = pickle.load(handle)

        pixel_cube_name = [cube_name] * nav_cube.shape[1]
        pixel_index = range(0, nav_cube.shape[1])
        pixel_cube_flybys = [cube_flybys[cube_name_all.index(cube_name)]] * nav_cube.shape[1]
        pixel_cube_time = [cube_time] * nav_cube.shape[1]

        longitudes.extend(longitude_all)
        latitudes.extend(latitude_all)

        sampling_all = nav_cube[31]  # pixel sampling
        nav_cube_all_24 = nav_cube[24]  # reflectance correction
        nav_cube_all_25.extend(nav_cube[25])  # N.P. in
        nav_cube_all_26.extend(nav_cube[26])  # S.P. in
        nav_cube_all_28 = nav_cube[28]  # incidence angle
        nav_cube_all_29 = nav_cube[29]  # emission angle
        nav_cube_all_30 = nav_cube[30]  # phase angle
        nav_cube_all_32 = nav_cube[32]  # azimuthal angle
        nav_cube_all_6 = nav_cube[6]
        nav_cube_all_7 = nav_cube[7]
        nav_cube_all_8 = nav_cube[8]
        nav_cube_all_9 = nav_cube[9]
        nav_cube_all_10 = nav_cube[10]
        nav_cube_all_15 = nav_cube[15]
        nav_cube_all_16 = nav_cube[16]
        nav_cube_all_17 = nav_cube[17]
        nav_cube_all_18 = nav_cube[18]
        nav_cube_all_19 = nav_cube[19]

        # read spectral info from cubes with planetaryimage module
        isisimage = CubeFile.open(root_path + '/C' + cube_name + '_ir.cub')
        cube_spectra = isisimage.data

        # create navigation cube dataframe to be pickled
        my_data = {'cube': pixel_cube_name,
                   'pixel_index': pixel_index,
                   'flyby': pixel_cube_flybys,
                   'time': pixel_cube_time,
                   'sampling': sampling_all,
                   'refl_cor': nav_cube_all_24,
                   'incidence': nav_cube_all_28,
                   'emission': nav_cube_all_29,
                   'phase': nav_cube_all_30,
                   'azimuthal': nav_cube_all_32,
                   'lon_0': nav_cube_all_6,
                   'lon_1': nav_cube_all_7,
                   'lon_2': nav_cube_all_8,
                   'lon_3': nav_cube_all_9,
                   'lon_4': nav_cube_all_10,
                   'lat_0': nav_cube_all_15,
                   'lat_1': nav_cube_all_16,
                   'lat_2': nav_cube_all_17,
                   'lat_3': nav_cube_all_18,
                   'lat_4': nav_cube_all_19}

        df1 = pd.DataFrame(my_data)

        for band_id in range(97, 353):
            values = cube_spectra[band_id - 97, :].T.flatten()
            # replace -3.4028227e+38 with np.nan
            values[values < 0] = np.nan
            values[values > 1] = np.nan
            df1[f'I/F_{band_id}'] = values

        df0 = pd.concat([df0, df1], ignore_index=True)

    polygon_coordinates = [longitudes, latitudes, nav_cube_all_25, nav_cube_all_26]

    # store map data (serialize)
    with open(f'_data/{output_name}.bin', 'wb') as handle:
        pickle.dump([df0, polygon_coordinates], handle, protocol=pickle.HIGHEST_PROTOCOL)

    return df0, polygon_coordinates

