"""
This module extracts from the cube header the necessary information for the calculation of the navigation cube.

INPUT:
CUB_NAME: is the name of the .CUB file which we want to retrieve the header information

AUTHOR: Christos Ntinos (IPGP, November 2021)

VERSION: 1.0

"""

def main(cube_name, root_path):
    """
    Parameters
    ----------
    cube_name : str
        cube name to be calculated
    root_path :
        file path  for .cub files downloaded from VIMS data portal

    Returns
    -------
    output : list
        r_size : int
            number of samples
        c_size : int
            number of lines
        target : str
            name of target
        ntvtim_start : str
            start time of acquisition
        ntvtim_stop : str
            end time of acquisition
        mirtim : float
            inter line delay duration (sec)
        R_offset :  int
            X-axis offset
        C_offset : int
            Z-axis offset
        modeid_IR : str
            sampling mode ("NORMAL" or "HI-RES")
        inttim : float
            exposure duration (sec)

    This function reads the header of the cube and returns key parameters.
    """

    file_name = root_path + "/C" + cube_name + "_ir.cub"

    # set the time correction for the instrument's clock
    BBLS_timecorrection = 1.01725

    # read cube header up to '  Group = Archive'
    header = []
    with open(file_name, 'r') as f:
        for line in iter(lambda: f.readline().rstrip(), '  Group = Archive'):
            header.append(line)

    # retrieve information from the header of the calibrated cube
    r_size = int(header[8].replace('      Samples = ', ''))
    c_size = int(header[9].replace('      Lines   = ', ''))
    target = header[24].replace('    TargetName                = ', '')
    ntvtim_start = header[27].replace('    StartTime                 = ', '')
    ntvtim_stop = header[28].replace('    StopTime                  = ', '')
    mirtim = float(header[31].replace('    InterlineDelayDuration    = ', '')) * BBLS_timecorrection/1000
    R_offset = int(header[32].replace('    XOffset                   = ', ''))
    C_offset = int(header[33].replace('    ZOffset                   = ', ''))
    modeid_IR = header[36].replace('    SamplingMode              = ', '')
    inttim = float(header[38].replace('    ExposureDuration          = (', '').split(" <IR>", 1)[0]) * BBLS_timecorrection/1000

    output = [r_size, c_size, target, ntvtim_start, ntvtim_stop, mirtim, R_offset, C_offset, modeid_IR, inttim]

    return output

