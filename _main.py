"""
Created on March, 2026

@author: Christos Ntinos

This Module runs the pipeline Super-Resolution for Icy Moons.

Version: 1.0

"""
from plot_map_new_v2 import plot_map
from constants import font_size
from save_nav_cubes_v2 import save_nav_cubes
from save_map_v2 import save_map
from save_cube_values_bns_v2 import save_cube_values_bns
from save_map_raw_v2 import save_map_raw
from radiometric_error_vims_v2 import radiometric_error_vims
from save_map_unsaturated_v2 import save_map_unsaturated
from wvln_shift_correction_v2 import wvln_shift_correction
from photometric_correction_coef_v2 import photometric_correction_coef
from photometric_correction_v2 import photometric_correction
from save_map_projected_v2 import save_map_projected
from save_map_sr_v2 import create_super_cells, load_map, calculate_overlapping_pixels, calculate_sr_map

import pickle
import matplotlib.pyplot as plt
import pandas as pd
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)  # supress df fragmentation warning

if __name__ == "__main__":

    # STEP 1: Calculate navigation cubes ###############################################################################

    root_path = '_Enceladus_cubes'  # for .cub files downloaded from VIMS data portal
    num_points = 5  # number of points per polygon-pixel edge
    target = 'ENCELADUS'  # impose target name because some cubes are targeted at SATURN
    method = 'DSK/UNPRIORITIZED'  # 'ELLIPSOID'  # 'DSK/UNPRIORITIZED'
    method_sampling = 'DIAGONAL'  # 'DISTANCE'  # 'DIAGONAL'
    mt_version = 'saved'

    # read cube names from text file
    f = open("Enceladus_cubes_test.txt", "r")
    cube_name_all = [x.replace("\n", "") for x in f.readlines()]
    f.close()

    exit_message = save_nav_cubes(cube_name_all, root_path, num_points, target, method, method_sampling, mt_version)
    print(exit_message)

    # STEP 2: Calculate map for all bands ##############################################################################

    # read cube names from text file
    f = open("Enceladus_cubes_ds4.txt", "r")
    cube_name_all = [x.replace("\n", "") for x in f.readlines()]
    f.close()

    # read cube flyby from text file
    f = open("Enceladus_cubes_flybys_ds4.txt", "r")
    cube_flybys = [x.replace("\n", "") for x in f.readlines()]
    f.close()

    # define output file (map) name
    output_name_bands = "global_map_all_bands_5_DSK_ds4"
    num_points = 5
    method = "DSK"
    method_sampling = "DIAGONAL"
    root_path = "_Enceladus_cubes"
    
    df, polygon_coordinates = save_map(cube_name_all, cube_flybys, num_points, method, method_sampling, root_path,
                                       output_name_bands)

    # STEP 3: VIMS error calculation ###################################################################################
    # define output file name for background noise (bns) extracted from cubes backplanes
    output_name_bns = output_name_bands + "_bns"
    df_bns = save_cube_values_bns(cube_name_all, output_name_bns, root_path)

    # define output file name for raw values (in DN) extracted from .qub files using IDL routine
    output_name_raw = output_name_bands + "_raw"
    # check if map is precalculated
    df_raw = save_map_raw(cube_name_all, output_name_raw, root_path)

    # define output file name for calculated sigma values (use of external files)
    output_name_sigma = output_name_bands + "_sigma"
    df_sigma = radiometric_error_vims(df_bns, df_raw, cube_name_all, output_name_sigma)

    # STEP 4: Wavelength shift #########################################################################################
    # read cube names from text file
    f = open("Enceladus_cubes_ds4.txt", "r")
    cube_name_all = [x.replace("\n", "") for x in f.readlines()]
    f.close()

    root_path = "_Enceladus_cubes"

    with open('_data/global_map_all_bands_5_DSK_ds4.bin', 'rb') as handle:
        df, polygon_coordinates = pickle.load(handle)

    with open('_data/df_all_bands_sigma_ds4.bin', 'rb') as handle:
        df_sigma = pickle.load(handle)

    with open('_data/df_all_bands_raw_ds4.bin', 'rb') as handle:
        df_raw = pickle.load(handle)

    with open('_data/df_all_bands_bns_ds4.bin', 'rb') as handle:
        df_bns = pickle.load(handle)

    output_name_wvln = 'global_map_all_wvlns_5_DSK_ds4'
    output_name_wvln_df = 'df_all_wvlns_[]_ds4'

    df_wvln_all, df_sigma, df_raw, df_bns_new = wvln_shift_correction(cube_name_all, root_path, df, polygon_coordinates,
                                                                      df_sigma, df_raw, df_bns, output_name_wvln,
                                                                      output_name_wvln_df)
    
    # STEP 5: Remove saturation ########################################################################################
    with open('_data/global_map_all_wvlns_5_DSK_ds4.bin', 'rb') as handle:
        df, polygon_coordinates = pickle.load(handle)

    with open('_data/df_all_wvlns_raw_ds4.bin', 'rb') as handle:
        df_raw = pickle.load(handle)

    with open('_data/df_all_wvlns_bns_ds4.bin', 'rb') as handle:
        df_bns = pickle.load(handle)

    margin = 0.97
    output_name_unsat = 'global_map_all_wvlns_5_DSK_ds4_unsaturated_with_bns_97'

    df_unsat = save_map_unsaturated(df, polygon_coordinates, df_raw, df_bns, margin, output_name_unsat)
    
    # STEP 6: Photometric correction coefficients ######################################################################
    # calculate coefficients a1, a0 of the linear fit of Aeq-g and their accompanying errors
    with open('_data/global_map_all_wvlns_5_DSK_ds4_unsaturated_with_bns_97.bin', 'rb') as handle:
        df_0, polygon_coordinates = pickle.load(handle)

    output_name_coef = 'photometric_coef_linear_all_wvlns_5_DSK_ds4_unsaturated_with_bns_97'
    df_coef = photometric_correction_coef(df_0, output_name_coef)

    # STEP 7: Photometric correction of map ############################################################################
    # define band (standard wavelength) at which to calculate the map
    band_id = 'I/F_3.09689'

    with open('_data/global_map_all_wvlns_5_DSK_ds4_unsaturated_with_bns_97.bin', 'rb') as handle:
        df_0, polygon_coordinates = pickle.load(handle)

    with open('_data/photometric_coef_linear_all_wvlns_5_DSK_ds4_unsaturated_with_bns_97.bin', 'rb') as handle:  # SAME for a0 and rf_cor??
        df_coef = pickle.load(handle)

    with open('_data/df_all_wvlns_sigma_ds4.bin', 'rb') as handle:
        df_sigma = pickle.load(handle)

    output_name_cor = 'global_map_wvln_3_09689_cor_0_0_0_5_DSK_ds4_unsaturated_with_bns_97'

    df_cor = photometric_correction(df_0, polygon_coordinates, df_coef, band_id, df_sigma, output_name_cor,
                                    method='rf_cor', target_geom=[0, 0, 0])

    fig, axs = plt.subplots(1, 1)
    plot_map(df_cor, polygon_coordinates, axs, band_id+'_cor', 'gray',
             f'cor I/F ({band_id.replace("I/F_", "").replace("_cor", "").replace("_sr", "")} '
             + r'$\mu m$)', vmin=0, vmax=1)
    axs.set_aspect('equal')
    axs.set_xlim(-180, 180)
    axs.set_ylim(-90, 90)
    axs.set_xticks(range(-180, 181, 60),
                   labels=[r'$180^{\circ}$W', r'$120^{\circ}$W', r'$60^{\circ}$W', r'$0^{\circ}$', r'60$^{\circ}$E',
                           r'120$^{\circ}$E', r'180$^{\circ}$E'], fontsize=font_size, ha='left')
    axs.set_yticks(range(-90, 91, 30),
                   labels=[r'90$^{\circ}$S', r'60$^{\circ}$S', r'30$^{\circ}$S', r'0$^{\circ}$', r'30$^{\circ}$N',
                           r'60$^{\circ}$N', r'90$^{\circ}$N'], fontsize=font_size, va='center')

    fig.suptitle(
        f'Enceladus at {band_id.replace("I/F_", "").replace("_sr", "")} μm corrected, 5 points/edge, DSK, data-set 4: 1566 cubes' + '\n' +
        r'target pixels, i, e $\leq$ 80$^{\circ}$, 0.001 $\leq$ I/F $\leq$ 1, sampling $\leq$ 78 km/pixel (156281 pixels)' +
        '\n saturation limit = (2**12-1-bns)*0.97, general correction')

    plt.show()

    # STEP 8: Save projected map #######################################################################################

    file_name = 'global_map_wvln_3_09689_cor_0_0_0_5_DSK_ds4_unsaturated_with_bns_97'

    df_proj, polygon_coordinates_proj = save_map_projected(file_name)

    # STEP 9: Calculate SR map #########################################################################################

    # define SR weight power (n=1, 2, 3,...)
    weight_power = 2

    # create new polygon pixels
    num_points_x = 1025  # num_cells+1
    num_points_y = 513

    # define the lower and upper limits for x and y
    minX, maxX, minY, maxY = -180, 180, -90, 90  # global

    # create regular rectangular grid of super-cells
    longitude_all_super, latitude_all_super = create_super_cells(minX, maxX, num_points_x, minY, maxY, num_points_y)

    # file name definition
    file_name = f'sr_global_map_ds4_wvln_3_09689_cor_0_0_0_5_DSK_unsaturated_with_bns_97_{num_points_x - 1}x{num_points_y - 1}'
    file_name_input = '_data/global_map_wvln_3_09689_cor_0_0_0_5_DSK_ds4_unsaturated_with_bns_97_projected.bin'

    # load photometrically corrected map of specific wavelength
    df, longitude_all, latitude_all = load_map(file_name_input)

    # calculate the polygons contained in each super-cell
    polygon_contains_all = calculate_overlapping_pixels(df, longitude_all, latitude_all, longitude_all_super,
                                                        latitude_all_super, file_name)

    # name of the band 'I/F_{standard_wavelength}'
    band_id = 'I/F_3.09689'

    # calculate the SR map
    df_sr = calculate_sr_map(df, band_id, weight_power, file_name, longitude_all, latitude_all, polygon_contains_all,
                             longitude_all_super, latitude_all_super)

