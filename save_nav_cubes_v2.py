"""
Created on May, 2022

@author: Christos Ntinos

This module calculates and saves the navigation cubes specified in a list.

Version: 2.0 (January 2026)
(+) mt_version
"""

import spiceypy as spice
import read_cube_header
import createVIMSnav
import pickle
import auto_mk_selection
from tqdm import tqdm
import os


def save_nav_cubes(cube_name_all, root_path, num_points, target, method, method_sampling, mt_version):
    """
    Parameters
    ----------
    cube_name_all : list of str
        list of cube names to be calculated
    root_path : str
        path for .cub files downloaded from VIMS data portal
    num_points : int
        number of points per polygon edge
    target : str
        name of target
    method: str
        "ELLIPSOID" or "DSK" for modeling the target body
    method_sampling: str
        "DISTANCE" or "DIAGONAL" for calculating pixel sampling
    mt_version : str
        flag for using saved metakernel or search online for latest, 'saved' or 'latest'

    Returns
    -------
    return message : str
        'Calculation is completed.

    This function calculates the navigation cube for a provided cube name list.
    """
    for counter in tqdm(range(0, len(cube_name_all)), desc="Navigation cubes"):

        cube_name = cube_name_all[counter]

        # load info from cube header
        r_size, c_size, target_cube, ntvtim_start, ntvtim_stop, mirtim, \
        R_offset, C_offset, modeid_IR, inttim = read_cube_header.main(cube_name, root_path)

        target_year = ntvtim_start[0:4]

        # load SPICE kernels
        if mt_version == 'latest':
            metakernel = auto_mk_selection.main(target_year)
            spice.spiceypy.furnsh(metakernel)
        elif mt_version == 'saved':
            mt_saved = os.listdir('_Cassini_kernels/mk/')
            mt_target_year = [file for file in mt_saved if target_year in file]
            spice.spiceypy.furnsh('_Cassini_kernels/mk/' + mt_target_year[0])  # load oldest version
        spice.spiceypy.furnsh('_Cassini_kernels/mk/additional_dsk.tm')  # load DSK

        # create navigation cube
        nav_cube, longitude_all, latitude_all, rec_corners_all, area_edges_in_all, output, area_edges_all\
            = createVIMSnav.main(cube_name, num_points, target, root_path, method, method_sampling)

        # create navigation cube list to be pickled
        my_data = [nav_cube, longitude_all, latitude_all]

        # store navigation cube (serialize)
        if method == 'ELLIPSOID':
            with open(f'_nav_cubes/num_points_{num_points}_ELL_' + method_sampling + f'/V{cube_name}.nav', 'wb') as handle:
                pickle.dump(my_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
        elif method == 'DSK/UNPRIORITIZED':
            with open(f'_nav_cubes/num_points_{num_points}_DSK_' + method_sampling + f'/V{cube_name}.nav', 'wb') as handle:
                pickle.dump(my_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            raise Exception('Error: Unrecognised method.')

        spice.spiceypy.kclear()  # clean kernel pool to avoid overloading

    return 'Calculation is completed.'
