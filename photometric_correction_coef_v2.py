'''
Created on September, 2022

@author: Christos Ntinos

This Module calculates and saves the coefficients of the linear fit for the equigonal albedo - phase angle points.

Version: 2.0 (January 2026)

'''

import pickle
import pandas as pd
from scipy import stats
from filter_map import filter_map_fun
from constants import standard_wvln_2, disc_function
from tqdm import tqdm


def photometric_correction_coef(df_0, output_name_coef):
    """
    Parameters
    ----------
    df: pandas.DataFrame
        map of radiance factor for all wavelengths
    output_name_coef: str
        name for binary output file

    Returns
    -------
    df_coef: pandas.DataFrame
        coefficients of linear fit for all wavelengths

    This function calculates coefficients a0, a1 of the linear fit of the Aeq-g points with their accompanying errors.
    """

    band_id_all = []

    # for wavelengths
    for i in standard_wvln_2:
        band_id_all.append(f'I/F_{i}')
    wavelength_all = standard_wvln_2

    slope_all = []
    intercept_all = []
    slope_stderr = []
    intercept_stderr = []

    for counter in tqdm(range(0, len(band_id_all)), desc="bands"):
        band_id = band_id_all[counter]
        # filter pixels
        df = filter_map_fun(df_0, band_id)

        # calculate equigonal albedo (the linear fit of which provides the phase function)
        inc = df['incidence'].values
        emi = df['emission'].values
        pha = df['phase'].values

        disc_fun = disc_function(inc, emi, pha)
        cube_value = df[band_id]
        df['Aeq'] = cube_value / disc_fun

        x = df[df['Aeq'].notnull()]['phase']
        y = df[df['Aeq'].notnull()]['Aeq']

        # calculate linear regression of data-set and uncertainty
        res = stats.linregress(x.astype(float), y.astype(float))
        slope_all.append(res.slope)
        intercept_all.append(res.intercept)
        slope_stderr.append(res.stderr)
        intercept_stderr.append(res.intercept_stderr)

    # save data (serialise)
    my_data = {'band': band_id_all,
               'wavelength': wavelength_all,
               'slope': slope_all,
               'intercept': intercept_all,
               'slope_stderr': slope_stderr,
               'intercept_stderr': intercept_stderr}

    df_coef = pd.DataFrame(my_data)

    # save data (serialize)
    with open(f'_data/{output_name_coef}.bin', 'wb') as handle:
        pickle.dump(df_coef, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return df_coef

