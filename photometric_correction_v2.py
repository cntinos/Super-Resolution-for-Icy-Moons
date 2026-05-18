'''
Created on June, 2022

@author: Christos Ntinos

This Module saves the photometrically corrected map for one band or for one wavelength.

Version: 2.0 (January 2026)

'''

import numpy as np
import pickle
from scipy.stats import gaussian_kde
from scipy import stats
from filter_map import filter_map_fun
from constants import disc_function


def photometric_correction(df_0, polygon_coordinates, df_coef, band_id, df_sigma, output_name_cor,
                           method='rf_cor', target_geom=None):
    """
    Calculates photometrically corrected map.

    Parameters
    ----------
    df_0 : pandas.DataFrame
        uncorrected map
    polygon_coordinates : list of lists
        polygon_coordinates = [longitudes, latitudes, nav_cube_all_25, nav_cube_all_26]
    df_coef : pandas.DataFrame
        coefficients of linear fit of equigonal albedo vs phase angle
    band_id : str
        standard wavelength at which to calculate map
    df_sigma : pandas.DataFrame
        standard deviation of radiance factor of uncorrected map
    output_name_cor : str
        name of output file containing the photometrically corrected map (.bin)
    method : str
        method for calculating the corrected radiance factor: 'rf_cor' for Xu et al. (2020) and 'a0' for equigonal
        albedo at zero phase angle (Robidel et al., 2020)
    target_geom : list of floats
        incidence, emission and phase angles in deg at which to normalise observations, applicable only for
        method = 'rf_cor'

    Returns
    -------
    df_cor : pandas.DataFrame
        photometrically corrected map
    """

    # filter pixels
    df = filter_map_fun(df_0, band_id)
    # define df with corrected values
    df_cor = df.iloc[:, :20]

    # calculate equigonal albedo
    a0_linear = df_coef[df_coef['band'] == band_id]['intercept'].values[0]
    a1_linear = df_coef[df_coef['band'] == band_id]['slope'].values[0]
    a0_linear_sigma = df_coef[df_coef['band'] == band_id]['intercept_stderr'].values[0]
    a1_linear_sigma = df_coef[df_coef['band'] == band_id]['slope_stderr'].values[0]

    inc = df['incidence'].values
    emi = df['emission'].values
    pha = df['phase'].values

    disc_fun = disc_function(inc, emi, pha)
    rf = df[band_id]
    phase_fun_all = rf / disc_fun
    df_cor['Aeq'] = phase_fun_all

    rf_sigma = df_sigma[f'{band_id}_sigma'].loc[df_cor.index]  # do I need .loc????
    phase_fun_sigma_all = rf_sigma / disc_fun
    df_cor['Aeq_sigma'] = phase_fun_sigma_all

    if method == 'a0':
        # calculate equigonal albedo at zero phase angle
        rf_cor_all = phase_fun_all - a1_linear * pha

        # calculate one sigma error of equigonal albedo at zero phase angle
        rf_cor_sigma_all = np.sqrt(phase_fun_sigma_all**2 + (pha*a1_linear_sigma)**2)

    elif method == 'rf_cor':
        # calculate corrected values
        inc_t, emi_t, pha_t = target_geom

        if target_geom == [0, 0, 0]:
            disc_fun_t = 1
        else:
            disc_fun_t = disc_function(inc_t, emi_t, pha_t)

        Aeq = a0_linear + a1_linear * pha
        Aeq_t = a0_linear + a1_linear * pha_t
        rf_cor_all = phase_fun_all * disc_fun_t * Aeq_t / Aeq

        # calculate one sigma error of corrected values
        Aeq_sigma = np.sqrt((a1_linear_sigma**2) * (pha**2) + a0_linear_sigma**2)
        Aeq_t_sigma = np.sqrt((a1_linear_sigma**2) * (pha_t**2) + a0_linear_sigma**2)
        rf_cor_sigma_all = np.abs(rf_cor_all) * np.sqrt((rf_sigma/rf)**2 + (Aeq_t_sigma/Aeq_t)**2 + (Aeq_sigma/Aeq)**2)
    else:
        raise Exception("Method not recognised; method='a0' OR 'rf_cor'.\n"
                        "For 'rf_cor' target geometry must be provided e.g. target_geom=[0, 0, 0]")

    df_cor[band_id + '_cor'] = rf_cor_all
    df_cor[band_id+'_cor_sigma'] = rf_cor_sigma_all

    # filter photometrically corrected I/F values before saving
    df_cor = df_cor[df_cor[band_id+'_cor'] <= 1]  # replace with nan instead?
    df_cor = df_cor[df_cor[band_id+'_cor'] >= 0]

    # store map data (serialize) ! polygon_coordinates saved without filtering to ensure index matching with df_cor
    with open(f'_data/{output_name_cor}.bin', 'wb') as handle:
        pickle.dump([df_cor, polygon_coordinates], handle, protocol=pickle.HIGHEST_PROTOCOL)

    return df_cor

