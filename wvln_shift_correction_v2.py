'''
Created on May, 2024

@author: Christos Ntinos

This Module resamples shifted spectra to the standard wavelengths.

Version: 2.0 (January 2026)

'''

import pickle
import json
import pandas as pd
from constants import standard_wvln_2
from scipy.interpolate import interp1d
import numpy as np
from tqdm import tqdm


def wvln_shift_correction(cube_name_all, root_path, df, polygon_coordinates, df_sigma, df_raw, df_bns, output_name_wvln,
                          output_name_wvln_df):
    """
    Parameters
    ----------
    cube_name_all: list of str
        names of the cubes to be used
    root_path: str
        path in repository to save cubes
    df: pandas.DataFrame
        map of radiance factor for all wavelengths
    polygon_coordinates: list of lists
        polygon_coordinates = [longitudes, latitudes, nav_cube_all_25, nav_cube_all_26]
    df_sigma: pandas.DataFrame
        map of sigma values for all wavelengths
    df_raw: pandas.DataFrame
        map of raw DN values for all wavelengths
    df_bns: pandas.DataFrame
        map of background noise (bns) for all wavelengths
    output_name_wvln: str
        name for binary output file (map)
    output_name_wvln_df: str
        name for binary output file (df)

    Returns
    -------
    df_wvln_all: pandas.DataFrame
        map of radiance factor for all wavelengths
    df_sigma: pandas.DataFrame
        map of sigma values for all wavelengths
    df_raw: pandas.DataFrame
        map of raw values for all wavelengths
    df_bns_new: pandas.DataFrame
        map of background noise values for all wavelengths

    This function corrects the input map for the VIMS wavelngth shift and ouputs maps in standard wavelengths.
    """

    wvln_all = []
    for counter in tqdm(range(0, len(cube_name_all)), desc="VIMS cubes"):
        cube_name = cube_name_all[counter]

        # load info from cube header
        file_name = root_path + "/C" + cube_name + "_ir.cub"

        # read cube header up to '  Group = Archive'
        header = []
        with open(file_name, 'r') as f:
            for line in iter(lambda: f.readline().rstrip(), '  Group = Kernels'):
                header.append(line)

        for j in range(0, len(header)):
            if 'Center' in header[j]:
                wvln_0 = header[j].replace('    Center         = (', '[')  # first line
                wvln_0 = wvln_0.replace('-', '')
                for i in range(j+1, len(header)):
                    if 'MissionAverage' not in header[i]:
                        if '-' in header[i]:
                            wvln_1 = header[i].replace('                      ', '')
                            wvln_1 = wvln_1.replace('-', '')
                        else:
                            wvln_1 = header[i].replace('                      ', '')
                        wvln_0 = wvln_0 + wvln_1
                    else:
                        wvln = wvln_0.replace(')', ']')
                        wvln = json.loads(wvln)  # convert string of list to list of strings
                        wvln = [float(i) for i in wvln]  # convert string to float

            else:
                pass

        if len(wvln) != 256:
            raise Exception("Spectral bands missing")
        else:
            wvln_all.append(wvln)

    # calculate shifted spectra
    index = df.index.values.tolist()

    spectrum_wvln = []
    spectrum_wvln_sigma = []
    spectrum_wvln_raw = []
    longest_na_gap_all = []
    pixel_excluded = 0
    counter_sigma = 0
    counter_raw = 0
    for i in index:
        # check match between calibrated, sigma and raw values (cube and pixel_index is the same)
        if (df.loc[i]['cube'] == df_sigma.loc[i]['cube'] == df_raw.loc[i]['cube']) & \
                (df.loc[i]['pixel_index'] == df_sigma.loc[i]['pixel_index'] == df_raw.loc[i]['pixel_index']):

            cube_name = df.loc[i]['cube']
            index_wvln = cube_name_all.index(cube_name)
            spectrum_band = df.loc[i][20:].values
            spectrum_sigma = df_sigma.loc[i][2:].values
            spectrum_raw = df_raw.loc[i][2:].values

            df_wvln = pd.DataFrame({'wvln': wvln_all[index_wvln],
                                    'value': spectrum_band,
                                    'sigma': spectrum_sigma,
                                    'raw': spectrum_raw})

            # calculate continuously missing data
            na_groups = df_wvln['value'].notna().cumsum()[df_wvln['value'].isna()]
            lengths_consecutive_na = na_groups.groupby(na_groups).agg(len)
            longest_na_gap = lengths_consecutive_na.max()
            if np.isnan(longest_na_gap):
                longest_na_gap = 0

            if df_wvln['value'].isnull().sum() < 256:  # remove only pixels with no bands captured
                if df_wvln['value'].isnull().sum() >= 1:
                    print(df.loc[i]['cube'], df.loc[i]['pixel_index'])

                f_linear = interp1d(df_wvln['wvln'].values, df_wvln['value'].values, kind='linear', bounds_error=False,
                                    fill_value=np.nan)
                y_linear = f_linear(standard_wvln_2)
                spectrum_wvln.append(y_linear)

                # calculate continuously missing data - sigma
                na_groups = df_wvln['sigma'].notna().cumsum()[df_wvln['sigma'].isna()]
                lengths_consecutive_na = na_groups.groupby(na_groups).agg(len)
                longest_na_gap = lengths_consecutive_na.max()
                if np.isnan(longest_na_gap):
                    longest_na_gap = 0
                if longest_na_gap > 5:
                    counter_sigma = counter_sigma + 1
                    print(f"Sigma values not available: {counter_sigma}")

                # interpolate sigma (even with gaps)
                f_linear = interp1d(df_wvln['wvln'].values, df_wvln['sigma'].values, fill_value="extrapolate")
                y_linear_sigma = f_linear(standard_wvln_2)
                spectrum_wvln_sigma.append(y_linear_sigma)

                # calculate continuously missing data - raw
                na_groups = df_wvln['raw'].notna().cumsum()[df_wvln['raw'].isna()]
                lengths_consecutive_na = na_groups.groupby(na_groups).agg(len)
                longest_na_gap = lengths_consecutive_na.max()
                if np.isnan(longest_na_gap):
                    longest_na_gap = 0
                if longest_na_gap > 5:
                    counter_raw = counter_raw + 1
                    print(f"Raw values not available: {counter_raw}")

                # interpolate raw (no missing values in raw)
                f_linear = interp1d(df_wvln['wvln'].values, df_wvln['raw'].values, fill_value="extrapolate")
                y_linear_raw = f_linear(standard_wvln_2)
                spectrum_wvln_raw.append(y_linear_raw)

            else:
                spectrum_wvln.append(len(standard_wvln_2)*[np.nan])
                spectrum_wvln_sigma.append(len(standard_wvln_2) * [np.nan])
                spectrum_wvln_raw.append(len(standard_wvln_2) * [np.nan])
                pixel_excluded = pixel_excluded + 1
                longest_na_gap_all.append(longest_na_gap)

        else:
            raise Exception("Mismatch of data!")

    print('Pixels missing more than 5 continuous spectral values', pixel_excluded)

    # save new map
    # define df with shift corrected values
    df_wvln_all = df.iloc[:, :20]
    df_wvln_all[[f'I/F_{standard_wvln_0}' for standard_wvln_0 in standard_wvln_2]] = spectrum_wvln
    df_wvln_all = df_wvln_all[~df_wvln_all[f'I/F_{standard_wvln_2[0]}'].isna()]  # remove non-reconstructed spectra

    # define df with shift corrected sigma values
    df_sigma = df.iloc[:, :20]
    df_sigma[[f'I/F_{standard_wvln_0}_sigma' for standard_wvln_0 in standard_wvln_2]] = spectrum_wvln_sigma
    df_sigma = df_sigma.loc[df_wvln_all.index.values.tolist()]  # remove non-reconstructed spectra

    # define df with shift corrected raw values
    df_raw = df.iloc[:, :20]
    df_raw[[f'I/F_{standard_wvln_0}_raw' for standard_wvln_0 in standard_wvln_2]] = spectrum_wvln_raw
    df_raw = df_raw.loc[df_wvln_all.index.values.tolist()]  # remove non-reconstructed spectra

    # store map data (serialize) ! polygon_coordinates saved without filtering to ensure index matching with df_cor
    with open(f'_data/{output_name_wvln}.bin', 'wb') as handle:
        pickle.dump([df_wvln_all, polygon_coordinates], handle, protocol=pickle.HIGHEST_PROTOCOL)

    # store map data (serialize) for sigma values
    with open(f'_data/{output_name_wvln_df.replace("[]", "sigma")}.bin', 'wb') as handle:
        pickle.dump(df_sigma, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # store map data (serialize) for raw values
    with open(f'_data/{output_name_wvln_df.replace("[]", "raw")}.bin', 'wb') as handle:
        pickle.dump(df_raw, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # wavelength shift for background noise (needed for saturation filtering)
    spectrum_wvln_bns = []
    for i in df_bns.index.values.tolist():
        cube_name = df_bns.loc[i]['cube']
        index_wvln = cube_name_all.index(cube_name)
        spectrum_bns = df_bns.loc[i][3:].values

        df_wvln = pd.DataFrame({'wvln': wvln_all[index_wvln],
                                'value': spectrum_bns})

        f_linear = interp1d(df_wvln['wvln'].values, df_wvln['value'].values, kind='linear', bounds_error=False,
                            fill_value=np.nan)
        y_linear = f_linear(standard_wvln_2)
        spectrum_wvln_bns.append(y_linear)

    # define df with shift corrected background noise
    df_bns_new = df_bns.iloc[:, :3]
    df_bns_new[[f'I/F_{standard_wvln_0}_bns' for standard_wvln_0 in standard_wvln_2]] = spectrum_wvln_bns

    # store map data (serialize)
    with open(f'_data/{output_name_wvln_df.replace("[]", "bns")}_bns.bin', 'wb') as handle:
        pickle.dump(df_bns_new, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return df_wvln_all, df_sigma, df_raw, df_bns_new

