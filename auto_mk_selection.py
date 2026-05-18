"""
Created on May, 2022

@author: Christos Ntinos

This Module checks in PDS for the latest metakernel of a certain target year; and if it does not match the one existing
locally it downloads it along with the kernels it lists.

Version: 1.0
"""

import requests
from bs4 import BeautifulSoup
import urllib.request
from pathlib import Path


def main(target_year):
    """
    Parameters
    ----------
    target_year : str
        year of cube acquisition

    Returns
    -------
    file_tm : str
        path to metakernel file

    This function checks for the latest version of the metakernel for a given year and downloads all kernels listed.
    """

    url_mk = 'https://naif.jpl.nasa.gov/pub/naif/pds/data/co-s_j_e_v-spice-6-v1.0/cosp_1000/extras/mk/'

    # check for latest version of metakernel
    response = requests.get(url_mk, params={})
    response_text = response.text
    soup = BeautifulSoup(response_text, 'html.parser')

    metakernels_in_PDS = [x for x in soup.text.split() if '.tm' in x]  # retrieve all mk names in PDS
    metakernels_target_year = [x for x in metakernels_in_PDS if f'{target_year}' in x]  # keep mk of the target year
    metakernels_target_year_last = max(metakernels_target_year)  # keep latest mk version

    mk_path_local = Path('_Cassini_kernels/mk/'+metakernels_target_year_last)
    if mk_path_local.is_file():
        pass
    else:  # mk not found locally
        urllib.request.urlretrieve(url_mk+metakernels_target_year_last, mk_path_local)  # download mk from PDS

        # change default path
        with open('_Cassini_kernels/mk/'+metakernels_target_year_last, 'r') as file:
            filedata = file.read()

        # Replace the target string
        filedata = filedata.replace("      PATH_VALUES     = ( './data' )",
                                    "      PATH_VALUES     = ( './_Cassini_kernels' )")

        # Write the file out again
        with open('_Cassini_kernels/mk/'+metakernels_target_year_last, 'w') as file:
            file.write(filedata)

    # check that all the kernels listed in the metakernel are available locally, if not download them from PDS
    url_kernels = 'https://naif.jpl.nasa.gov/pub/naif/CASSINI/kernels/'
    url_generic = 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/'

    mk_file = open('_Cassini_kernels/mk/'+metakernels_target_year_last, 'r')
    lines = mk_file.readlines()
    kernels_in_mk = [x for x in lines if '$KERNELS' in x]
    kernels = [x.replace("                          '$KERNELS/", "") for x in kernels_in_mk]
    kernels = [x.replace("'\n", "") for x in kernels]

    for kernel in kernels:
        kernel_path_local = Path('_Cassini_kernels/'+kernel)
        if kernel_path_local.is_file():  # kernle exist locally
            pass
        else:
            try:  # mission specific kernel
                urllib.request.urlretrieve(url_kernels + kernel, '_Cassini_kernels/'+kernel)  # download mk from PDS
            except:
                try:
                    urllib.request.urlretrieve(url_generic + kernel, '_Cassini_kernels/' + kernel)
                except:
                    url_generic_2 = 'https://naif.jpl.nasa.gov/pub/naif/pds/data/co-s_j_e_v-spice-6-v1.0/cosp_1000/data/'
                    urllib.request.urlretrieve(url_generic_2 + kernel, '_Cassini_kernels/' + kernel)

    file_tm = '_Cassini_kernels/mk/'+metakernels_target_year_last

    return file_tm

