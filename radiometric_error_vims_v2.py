'''
Created on October, 2023

@author: Christos Ntinos

This Module calculates the standard deviation of the radiance factor I/F for the original (uncorrected) dataset based on
the analysis of [Es-sayeh et al., 2023].

Version: 2.0 (January 2026)

Adapted from Mael Es-sayer (IPGP/UPC, 2022) in IDL
Adapted from Thomas Cornet in IDL
Adapted from Sebastien Rodriguez in IDL

'''

import pickle
from statistics import mean, pstdev
import numpy as np
from scipy.interpolate import interp1d
from tqdm import tqdm

# read saved background noise
import pandas as pd


def radiometric_error_vims(df_bns, df_raw, cube_name_all, output_name_sigma):
    """
    Parameters
    ----------
    df_bns : pandas.DataFrame
        dataframe of the background noise (bns) for all pixels in map
    df_raw : pandas.DataFrame
        dataframe of the raw data (DN) for all pixels in map
    cube_name_all: list of str
        names of the cubes to be used
    output_name_sigma: str
        name for binary output file

    Returns
    -------
    df_sigma: pandas.Dataframe
        map of one sigma error of I/F for all bands

    This function calculates the map of the one sigma error of the radiance factor for all bands.
    """
    columns = []
    for band_id in range(97, 353):
        columns.append(f'I/F_{band_id}_bns_mean')

    bns_mean_array = np.empty((len(cube_name_all), 256))
    # calculate mean background noise for each cube
    i = 0  # counter

    for counter in tqdm(range(0, len(cube_name_all)), desc="VIMS cubes"):
        cube_name = cube_name_all[counter]
        df0 = df_bns[df_bns['cube'] == cube_name]
        for band_id in range(97, 353):

            bns_mean = mean(df0[f'I/F_{band_id}_bns'])
            bns_srd = pstdev(df0[f'I/F_{band_id}_bns'])

            bns_mean_new = mean(df0[df0[f'I/F_{band_id}_bns'] <= bns_mean+3*bns_srd][f'I/F_{band_id}_bns'])

            bns_mean_array[i, band_id-97] = bns_mean_new

        i = i + 1

    df_mean_bns = pd.DataFrame({'cube': cube_name_all})
    df_mean_bns[columns] = bns_mean_array

    # caclulate radiometric noise of VIMS
    gain = 300  # VIMS gain in DN^-1
    read_noise = 1  # in DN

    # read VIMS standard wavelength and Full Width Half Maximum (FWHM) for each channel [Clarke et al., 2018]
    with open('_shared_files_and_code/VIMS_IR_FWHM.txt') as f:
        lines = f.readlines()

    wavelength_all = []
    fwhm_all = []
    for line in lines:
        wavelength_all.append(float(line.split()[1]))  # standard wavelength in microns
        fwhm_all.append(float(line.split()[2]))  # FWHM in microns

    # read VIMS transfer function
    with open('_shared_files_and_code/ftransvis-IR_11-06-04.txt') as f:
        lines = f.readlines()

    transfer_fun = []
    for line in lines:
        transfer_fun.append(float(line.split()[0]))

    # read Sun radiance factor
    with open('_shared_files_and_code/vims_solar_func.txt') as f:
        lines = f.readlines()

    sun_wavelength = []
    sun_flux = []
    for line in lines:
        sun_wavelength.append(float(line.split()[1]))
        sun_flux.append(float(line.split()[0]))

    # interpolate solar spectrum to standard wavelength
    interp_func = interp1d(np.array(sun_wavelength), np.array(sun_flux)/(10**6))
    sun_flux_inter = interp_func(np.array(wavelength_all))

    # calculate standard deviation for I/F for each pixel and band
    df_sigma = df_raw.copy()[['cube', 'pixel_index']]

    c = 2.998*10**8  # speed of light in m/s
    h = 6.626*10**(-34)  # Planck's constant in Js
    alpha = 0.035  # IR telescope aperture in m^2 !!! different from report
    omega = 2.5*10**(-7)  # FOV solid angle subtended in srad
    distance = 9.5475  # mean Sun - Saturn distance in au

    for band_id in tqdm(range(97, 353), desc="bands"):

        sigma_all = []
        for i in df_raw.index.values.tolist():
            cube_name = df_raw['cube'][i]
            mean_bns = df_mean_bns[f'I/F_{band_id}_bns_mean'][df_mean_bns['cube'] == cube_name].values[0]
            exposure = df_bns['exposure'][i]  # exposure time in sec
            signal = df_raw[f'I/F_{band_id}_raw'][i]

            # shot noise of signal
            sn_signal = np.sqrt(signal * gain + mean_bns * gain)/gain

            # shot noise of background
            sn_back = np.sqrt(mean_bns * gain)/gain

            # total shot noise
            sn_total = np.sqrt(sn_signal**2 + sn_back**2)

            # total noise in DN
            total_noise = np.sqrt(sn_total**2 + read_noise**2)

            # standard deviation of I/F
            sigma = total_noise * (transfer_fun[band_id-97]/exposure) * (h*c/(wavelength_all[band_id-97]*10**(-6))) * \
                    (1/(alpha*omega*fwhm_all[band_id-97])) * (1/(sun_flux_inter[band_id-97]/(np.pi*distance**2)))  # !!! check devide by distance

            sigma_all.append(sigma)

        df_sigma[f'I/F_{band_id}_sigma'] = sigma_all

    # save dataframe
    with open(f'_data/{output_name_sigma}.bin', 'wb') as handle:
        pickle.dump(df_sigma, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return df_sigma

