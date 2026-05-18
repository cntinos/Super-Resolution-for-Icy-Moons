"""
Created on June, 2025

@author: Christos Ntinos (IPGP)

This Module identifies the parts of each calibrated spectrum of the global map that are saturated based on the
respective raw values in DN and replaces them with NaN.

Version: 2.0 (January 2026)
max DN = 2**12-1 NOT 2**12

"""


import pickle
from constants import standard_wvln_2
import numpy as np
from tqdm import tqdm


def save_map_unsaturated(df, polygon_coordinates, df_raw, df_bns, margin, output_name_unsat):
    """
    Parameters
    ----------
    df: pandas.DataFrame
        map of radiance factor for all wavelengths
    polygon_coordinates: list of lists
        polygon_coordinates = [longitudes, latitudes, nav_cube_all_25, nav_cube_all_26]
    df_raw: pandas.DataFrame
        map of raw DN values for all wavelengths
    df_bns: pandas.DataFrame
        map of background noise (bns) for all wavelengths
    margin: float
        margin below max DN value to consider saturation, margin=0.97 suggested value
    output_name_unsat: str
        name for binary output file

    Returns
    -------
    df_unsat: pandas.DataFrame
        map of radiance factor for all bands, with nan where saturated

    This function reads a map of all the wavelengths and replaces with nan the values corresponding to raw DN values
    above a certain threshold (saturated).
    """

    spectra_all = []
    for i in tqdm(range(0, df.shape[0]), desc="pixels"):
        spectrum = df.iloc[i][20:].values
        spectrum_raw = df_raw[(df_raw['cube'] == df.iloc[i]['cube']) & (
                    df_raw['pixel_index'] == df.iloc[i]['pixel_index'])].values.flatten()[20:]
        spectrum_bns = df_bns[(df_bns['cube'] == df.iloc[i]['cube']) & (
                df_bns['pixel_index'] == df.iloc[i]['pixel_index'])].values.flatten()[3:]
        spectrum = np.where(spectrum_raw >= ((2**12-1)-spectrum_bns)*margin, np.nan, spectrum)  # dynamic saturation limit
        spectra_all.append(spectrum)

    # pack data
    df_unsat = df.iloc[:, :20]
    df_unsat[[f'I/F_{standard_wvln_0}' for standard_wvln_0 in standard_wvln_2]] = spectra_all

    # store data (serialize) for polygons contained in each SR map
    with open(f'_data/{output_name_unsat}.bin', 'wb') as handle:
        pickle.dump([df_unsat, polygon_coordinates], handle, protocol=pickle.HIGHEST_PROTOCOL)

    return df_unsat

