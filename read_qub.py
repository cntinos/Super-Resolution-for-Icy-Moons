"""
Created on October 2026

@author: Christos Ntinos (IPGP)

This Module reads the raw VIMS data from a .qub file (binary parser for VIMS raw cubes).

Version: 1.0
"""

import numpy as np
import re


def read_label(raw):
    """
    Parameters
    ----------
    raw: bytes
        binary raw cube file

    Returns
    -------
    label: str
        raw cube header

    This function return the attached label text of the cube (everything up to the line that says END).
    """
    text = raw[:100000].decode('latin1')
    end = re.search(r'\r?\nEND\s*\r?\n', text)

    if end is None:
        raise ValueError('no END line found: not a PDS3 label?')

    label = text[:end.start()]
    return label


def keyword(text, name):
    """
    Parameters
    ----------
    text: str
        attached label to search
    name: str
        keyword string to look for

    Returns
    -------
    value: str
        value string

    This function isolates and returns the value of one label keyword from the cube's attached label.
    """
    match = re.search(r'^\s*%s\s*=\s*(.+?)\s*$' % re.escape(name), text, re.MULTILINE)
    if match is None:
        raise ValueError('keyword %s not found in the label' % name)

    value = match.group(1)
    return value


def read_qub(path):
    """
    Parameters
    ----------
    path: str
        path of raw cube (.qub) to read

    Returns
    -------
    core: np.array
        array of raw data (DN) for all bands (VIS + IR)
    bkg: np.array
        array of background noise (DN) for each line and for all bands (VIS + IR)
    header: str
        raw cube header

    This function opens a VIMS raw cube (.qub) and outputs the data in DN (spectral cube and background noise).
    """

    # open raw cube
    raw = open(path, 'rb').read()

    # read attached label and isolate parameters for parsing
    label = read_label(raw)
    rb = int(keyword(label, 'RECORD_BYTES'))
    start = (int(re.search(r'\^QUBE\s*=\s*(\d+)', label).group(1)) - 1)*rb
    ns, nb, nl = [int(x) for x in keyword(label, 'CORE_ITEMS').strip('()').split(',')]
    ss, bs, _ = [int(x) for x in keyword(label, 'SUFFIX_ITEMS').strip('()').split(',')]
    dtype_core = f">i{int(keyword(label, 'CORE_ITEM_BYTES'))}"
    dtype_suffix = f">i{int(keyword(label, 'SUFFIX_BYTES'))}"

    core = np.zeros((nb, nl, ns))
    bkg = np.zeros((nb, nl))
    p = start

    for l in range(nl):
        for b in range(nb):
            core[b, l] = np.frombuffer(raw, dtype_core, ns, p)
            p += 2*ns
            bkg[b, l] = np.frombuffer(raw, dtype_suffix, ss, p)[0]
            p += 4*ss
        p += bs*(ns + ss)*4  # band-suffix planes (housekeeping)

    return core.astype(float), bkg.astype(float)

