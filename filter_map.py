'''
Created on August, 2023

@author: Christos Ntinos

This Module filters the dataframe containing the all the pixels in the data set, thus ensuring consistently filtering
across the pipeline.

Version: 1.0

'''


def filter_map_fun(df, band_id):
    """
    Parameters
    ----------
    df : pandas.DataFrame
        input dataframe
    band_id : str
        wavelength of the band to be filtered

    Returns
    -------
    df : pandas.DataFrame
        filtered dataframe

    This function filters the input df based on angles, sampling etc.
    """
    # filter data !! must use the same filters as in photometric_correction_coef_v2.py (recalculation)

    sampling_limit = 78  # km/pixel
    df = df[df['sampling'] <= sampling_limit]  # filter for sampling smaller than limit
    df = df[df['incidence'] <= 80]  # filter out night-side pixels
    df = df[df['emission'] <= 80]

    if band_id == 'I/F_3.09689':  # cut-off phase angle due to pole at g=-a0/a1
        df = df[df['phase'] <= 100]

    df = df[df[band_id] >= 0.001]
    df = df[df[band_id] <= 1]
    df = df[df['refl_cor'] == 1]

    return df

