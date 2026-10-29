"""
This module creates a nav_cube from kernel/SPICE functions.
The navigation cube produced refers only to the IR channel.

INPUT:
CUB_NAME: is the name of the .CUB file you want to create the nav
nav:
Layer #1     : Intersept of center of pixel not found/found [0, 1]
Layer #2-5   : Intersept of corners of pixel not found/found [0, 1]
Layer #6     : Number of corners of projected pixel [2, 3, 4, 5, 6, 7, 8]
Layer #7-15  : Longitude (E), center and up to 8 corners {CLW, TL=1}  [deg]
Layer #16-24 : Latitude (N),  center and up to 8 corners {CLW, TL=1}  [deg]
Layer #25    : Reflectance correction coefficient [>=1] / Ain/Atotal<=1
Layer #26    : North pole contained in the pixel [0, 1]
Layer #27    : South pole contained in the pixel [0, 1]
Layer #28    : Time et of pixel [...]
Layer #29    : Incidence angle (Ang B/W Illsrc & Srfnrm) [deg]
Layer #30    : Emission angle  (Ang B/W Obsrvr & Srfnrm) [deg]
Layer #31    : Phase angle     (Ang B/W Inc & Emi || to Obs&Src)[deg]
Layer #32    : Pixel sampling [km/pixel]
Layer #33    : Azimuthal angle between the planes of incidence and emergence [deg]

AUTHOR: Christos Ntinos (IPGP, October 2021) in Python
        Adapted from Nicholas Kutsop (Cornell/JPL, October  2018) in Matlab
        Adapted from Benoit Seignovert (ENSL/JPL,  May-July   2013)
        Adapted from Sebastien Rodriguez (LPGNantes, April 2008) in IDL

VERSION: 1.0

"""

import numpy as np
import spiceypy as spice
import read_cube_header
import CoverageTool_RectangularPIXEL
import refl_correction


def main(cube_name, num_points, target, root_path, method, method_sampling):
    """
    Parameters
    ----------
    cube_name : str
        cube name to be calculated
    num_points : int
        number of points per pixel-polygon edge
    target : str
        name of target
    root_path :
        file path  for .cub files downloaded from VIMS data portal
    method: str
        "ELLIPSOID" or "DSK" for modeling the target body
    method_sampling: str
        "DISTANCE" or "DIAGONAL" for calculating pixel sampling

    Returns
    -------
    nav_cube : ndarray
         navigation cube
    longitude_all : list of lists
        longitudes of all points in all pixels-polygons
    latitude_all : list of lists
        latitude of all points in all pixels-polygons
    rec_corners_all : list of lists
        corners of all pixels-polygons in rectangular coordinates (km)
    area_edges_in_all : list of ndarray
        for all pixels the points that define the edges of the part of the target inside the pixel (as projected on the
        limb plane) constructed py pixel corners and limb segments in rectangular coordinates (km)
    output : list of ndarray
        all points comprising the edges of all the pixels-polygons in rectangular coordinates (km)
    area_edges_all : list of ndarray
        for all pixels the points that define the edge of the projected pixel_polygon on the target in rectangular
        coordinates (km)

    This function creates the navigation cube for a given cube name.
    """

    # load info from cube header
    r_size, c_size, target_cube, ntvtim_start, ntvtim_stop, mirtim, \
    R_offset, C_offset, modeid_IR, inttim = read_cube_header.main(cube_name, root_path)

    decnames = 'IR'
    VIMSDN = 'IR'

    instnm = 'CASSINI_VIMS_IR'  # instrument name
    instid = spice.spiceypy.bodn2c(instnm)  # instrument ID

    abcorr = 'CN+S'
    room = 10  # max number of outputs for spice functions
    obsrvr = 'CASSINI'
    spc_id = spice.spiceypy.bodn2c('CASSINI')

    et_stop = spice.spiceypy.str2et(ntvtim_stop)
    et_start = spice.spiceypy.str2et(ntvtim_start)

    # body-fixed reference frame name
    fixref = 'IAU_' + target

    # retrieve target's radial dimensions
    dim, values = spice.spiceypy.bodvrd(target, 'RADII', 3)
    radii = values[0]

    rdecsz = 64  # r_size theorical max
    cdecsz = 64  # c_size theorical max

    IRexposec = inttim * r_size  # whiskbroom scanner (one pixel at a time)

    start_fstpix = et_start

    # calculate the time for each pixel
    IRpix_intstarts = np.zeros((r_size, c_size))

    for i in range(0, r_size):
        for j in range(0, c_size):
            IRrow_scantimes = i * inttim  # timing between rows
            IRlinedelay_times = j * (mirtim + r_size * inttim)  # timing between columns
            IRpix_intstarts[i][j] = IRlinedelay_times + IRrow_scantimes  # timing between each pixel

    start_arrpix = start_fstpix + IRpix_intstarts  # start time for each pixel in et
    et_arr = start_arrpix + 0.5 * inttim

    # retrieve instrument FOV parameters
    shape, frame, bsight, n, bounds = spice.spiceypy.getfov(instid, room)

    pxfovr = np.amax(np.abs(bounds[1] - bounds[2])) / rdecsz  # calculate dimensions of the FOV
    pxfovc = np.amax(np.abs(bounds[0] - bounds[1])) / cdecsz

    R_TL = bounds[1][0]
    C_TL = bounds[1][1]
    L_TL = bounds[1][2]

    shft = np.array([[0, 1, 1, -1, -1], [0, 1, -1, -1, 1], [0, 0, 0, 0, 0]])  # 3x5 array of corners and center of pixel
    j, h = np.meshgrid(np.linspace(0, c_size), np.linspace(0, r_size))

    if modeid_IR == 'NORMAL':
        ifovr = pxfovr
        ifovc = pxfovc
        rthrow = np.fix(R_offset - 1) * ifovr
        cthrow = np.fix(C_offset - 1) * ifovc
        sampl_coef = 0.0005  # used to calculate pixel sampling
    elif modeid_IR == 'HI-RES':
        ifovr = pxfovr / 2
        ifovc = pxfovc
        rthrow = np.fix(2 * (R_offset + (r_size - 1) / 4)) * ifovr
        cthrow = np.fix(C_offset - 1) * ifovc
        sampl_coef = 0.00025
    else:
        raise Exception(f"Unknown instrument mode {modeid_IR}")

    bsite = np.array([0.5 * ifovr, 0.5 * ifovc, 0])
    ifov_2 = np.array([0.5 * ifovr, 0.5 * ifovc, 0])  # ==bsite??

    # create
    tmpc, tmpr = np.meshgrid(
        np.arange(C_TL + cthrow, C_TL + cthrow + c_size * ifovc, ifovc),
        np.arange(R_TL + rthrow, R_TL + rthrow + r_size * ifovr, ifovr))

    tmpl = np.zeros(tmpc.shape) + L_TL
    tmpi = np.array((tmpr.flatten().T, tmpc.flatten().T, tmpl.flatten().T))

    arri = tmpi + np.tile(bsite, (tmpi.shape[1], 1)).T + np.tile(ifov_2, (tmpi.shape[1], 1)).T

    # retrieve ellipsoid semi-major axes
    dim, values = spice.spiceypy.bodvrd(target, 'RADII', 3)
    radii_a = values[0]
    radii_b = values[1]
    radii_c = values[2]

    output = []
    corners_limb = []
    rec_corners_all = []
    area_edges_all = []
    area_edges_in_all = []
    longitude_all = []
    latitude_all = []

    # create nav cube
    nav_cube = np.empty((33, r_size * c_size))
    nav_cube[:] = np.nan

    for i in range(r_size * c_size):
        et = et_arr.flatten()[i]

        # calcualte pixel sampling
        pos_spacecraft, litime = spice.spiceypy.spkpos(obsrvr, et, f'IAU_{target}', abcorr, target)

        for g in range(0, 5):
            try:
                dvec = arri[:, i] + shft[:, g] * ifov_2  # ray's direction vector
                spoint, trgepc, srfvec = spice.spiceypy.sincpt(method, target, et, fixref, abcorr, obsrvr, frame, dvec)
                nav_cube[0+g][i] = 1  # Layer 1-5, intercept of center/corners of pixel found

            except:
                nav_cube[0+g][i] = 0  # Layer 1-5, intercept NOT found

        if nav_cube[1][i] == 1 and nav_cube[2][i] == 1 and nav_cube[3][i] == 1 and nav_cube[4][i] == 1:
            # case of pixel completely inside the target
            for g in range(0, 5):
                dvec = arri[:, i] + shft[:, g] * ifov_2  # ray's direction vector
                spoint, trgepc, srfvec = spice.spiceypy.sincpt(method, target, et, fixref, abcorr, obsrvr, frame, dvec)

                radius, pclon, pclat = spice.spiceypy.reclat(spoint)

                nav_cube[6 + g][i] = np.rad2deg(pclon)  # Layer 7-11: longitude, center & corners [deg]
                nav_cube[15 + g][i] = np.rad2deg(pclat)  # Layer 16-20: latitude, center & corners [deg]

            nav_cube[11][i] = np.nan  # Layer 12-15: no 5th - 8th corners
            nav_cube[12][i] = np.nan
            nav_cube[13][i] = np.nan
            nav_cube[14][i] = np.nan
            nav_cube[20][i] = np.nan  # Layer 21-24: no 5th - 8th corners
            nav_cube[21][i] = np.nan
            nav_cube[22][i] = np.nan
            nav_cube[23][i] = np.nan

            # calculate pixel projection using the Coverage Tool
            dvec_all = []
            for j in range(1, 5):  # corners only
                dvec_all.append(
                    arri[:, i] + shft[:, j] * ifov_2)  # ray's direction vector for the 4 corners of the pixel

            area_corners, num_corners, output_all, area_edges, rec_corners, area_edges_in, limb_plane, north_in, south_in = \
                    CoverageTool_RectangularPIXEL.main(et, target, instnm, obsrvr, num_points, abcorr, dvec_all, method)
            output.extend(output_all)

            nav_cube[24][i] = 1  # Atotal/Ain
            nav_cube[25][i] = north_in
            nav_cube[26][i] = south_in

            # calculate photometric parameters
            pha_all = np.zeros(num_corners)
            inc_all = np.zeros(num_corners)
            emi_all = np.zeros(num_corners)
            for j in range(0, num_corners):
                trgepc, srfvec, pha_all[j], inc_all[j], emi_all[j] = spice.spiceypy.illumg(method, target, 'SUN',
                                                                  et, fixref, abcorr, obsrvr, area_corners[j])

            nav_cube[27][i] = et  # Layer 28
            nav_cube[28][i] = np.rad2deg(np.sum(inc_all)/num_corners)  # Layer 29: incidence angle [deg]
            nav_cube[29][i] = np.rad2deg(np.sum(emi_all)/num_corners)  # Layer 30: emission angle [deg]
            nav_cube[30][i] = np.rad2deg(np.sum(pha_all)/num_corners)  # Layer 31: phase angle [deg]

            # calculate pixel sampling
            if method_sampling == 'DIAGONAL':
            # by calculating the max length edge of the pixel
                edge_length_all = []
                for j in [*range(0, num_corners), 0]:
                    edge_length = 0
                    for k in range(0, num_points-1):
                        edge_length = edge_length + np.linalg.norm(area_edges[j][k+1]-area_edges[j][k])
                    edge_length_all.append(edge_length)

                # include diagonals
                for j in range(0, num_corners):
                    index = [*range(0, num_corners)]
                    for k in index:
                        point_1 = area_corners[j]
                        point_2 = area_corners[k]
                        diagonal_length = np.linalg.norm(point_2-point_1)
                        edge_length_all.append(diagonal_length)

                nav_cube[31][i] = max(edge_length_all)

            elif method_sampling == 'DISTANCE':
                # or more simply using the distance from the S/C
                vector_center = np.zeros(3)
                for j in range(0, num_corners):
                    vector_center = vector_center + area_corners[j]
                spoint = spice.spiceypy.surfpt([0, 0, 0], vector_center, radii_a, radii_b, radii_c)
                distance = np.linalg.norm(spoint - pos_spacecraft)
                nav_cube[31][i] = distance * sampl_coef

            else:
                raise Exception('Error: Unrecognised sampling method.')

            # calculate azimuthal angle (needed for the shadowing function in 6 parameter Hapke model)
            azccw = 1
            elplsz = 0
            obspos = spoint  # rectangular coordinates of pixel center
            azlsta_e, lt = spice.spiceypy.azlcpo('ELLIPSOID', obsrvr, et, abcorr, azccw, elplsz, obspos, target, fixref)
            azlsta_i, lt = spice.spiceypy.azlcpo('ELLIPSOID', 'SUN', et, abcorr, azccw, elplsz, obspos, target, fixref)

            if np.rad2deg(np.abs(azlsta_e[1]-azlsta_i[1]))<=180:
                nav_cube[32][i] = np.rad2deg(np.abs(azlsta_e[1]-azlsta_i[1]))
            else:
                nav_cube[32][i] = np.rad2deg(np.abs(azlsta_e[1] - azlsta_i[1]))-180

        elif nav_cube[1][i] == 0 and nav_cube[2][i] == 0 and nav_cube[3][i] == 0 and nav_cube[4][i] == 0:
            # case of pixel completely outside the target
            num_corners = 0
            nav_cube[24][i] = 0  # Atotal/Ain
            rec_corners = []
            area_edges_in = []
            area_edges = []

        else:
            # case of pixel crossing the limb
            nav_cube[24][i] = 0.5  # Atotal/Ain
            if method == 'DSK/UNPRIORITIZED':
                # for method='DSK/UNPRIORITIZED' we ignore limb-crossing pixels
                num_corners = 0
                rec_corners = []
                area_edges_in = []
                area_edges = []

            elif method == 'ELLIPSOID':
                dvec_all = []
                for g in range(1, 5):  # corners only
                    dvec_all.append(
                        arri[:, i] + shft[:, g] * ifov_2)  # ray's direction vector for the 4 corners of the pixel

                area_corners, num_corners, output_all, area_edges, rec_corners, area_edges_in, limb_plane, north_in, south_in = \
                    CoverageTool_RectangularPIXEL.main(et, target, instnm, obsrvr, num_points, abcorr, dvec_all, method)

                if num_corners != 0:  # agreement between nav_cub.T[i][1:5] and flag1
                    nav_cube[5][i] = num_corners  # Layer 6: number of corners of projected pixel
                    for j in range(0, num_corners):
                        radius, pclon, pclat = spice.spiceypy.reclat(area_corners[j])
                        nav_cube[7 + j][i] = np.rad2deg(pclon)  # Layer 8-15: longitude corners up to 8 [deg]
                        nav_cube[16 + j][i] = np.rad2deg(pclat)  # Layer 17-24: latitude corners up to 8 [deg]

                    try:
                        nav_cube[24][i] = refl_correction.main(area_corners, num_corners, limb_plane,
                                                               area_edges_in, rec_corners)
                    except:  # fully inside the target pixel originally recognised as limb-crossing ???
                        nav_cube[24][i] = 1

                    nav_cube[25][i] = north_in
                    nav_cube[26][i] = south_in

                    # calculate photometric parameters
                    pha_all = np.zeros(num_corners)
                    inc_all = np.zeros(num_corners)
                    emi_all = np.zeros(num_corners)
                    for j in range(0, num_corners):
                        trgepc, srfvec, pha_all[j], inc_all[j], emi_all[j] = spice.spiceypy.illumg('ELLIPSOID', target, 'SUN',
                                                                                                   et, fixref, abcorr, obsrvr,
                                                                                                   area_corners[j])

                    # nav_cube[0][i] = 1  # Layer: intercept found/not found [0/1]
                    nav_cube[27][i] = et  # Layer 28
                    nav_cube[28][i] = np.rad2deg(np.sum(inc_all) / num_corners)  # Layer 29: incidence angle [deg]
                    nav_cube[29][i] = np.rad2deg(np.sum(emi_all) / num_corners)  # Layer 30: emission angle [deg]
                    nav_cube[30][i] = np.rad2deg(np.sum(pha_all) / num_corners)  # Layer 31: phase angle [deg]

                    # calculate pixel sampling
                    if method_sampling == 'DIAGONAL':
                        # by calculating the max length edge of the pixel
                        edge_length_all = []
                        for j in [*range(0, num_corners), 0]:
                            edge_length = 0
                            for k in range(0, num_points-1):
                                edge_length = edge_length + np.linalg.norm(area_edges[j][k+1]-area_edges[j][k])
                            edge_length_all.append(edge_length)

                        # include diagonals
                        for j in range(0, num_corners):
                            index = [*range(0, num_corners)]
                            for k in index:
                                point_1 = area_corners[j]
                                point_2 = area_corners[k]
                                diagonal_length = np.linalg.norm(point_2 - point_1)
                                edge_length_all.append(diagonal_length)

                        nav_cube[31][i] = max(edge_length_all)

                    elif method_sampling == 'DISTANCE':
                        # or more simply using the distance from the S/C
                        vector_center = np.zeros(3)
                        for j in range(0, num_corners):
                            vector_center = vector_center + area_corners[j]
                        spoint = spice.spiceypy.surfpt([0, 0, 0], vector_center, radii_a, radii_b, radii_c)
                        distance = np.linalg.norm(spoint-pos_spacecraft)
                        nav_cube[31][i] = distance * sampl_coef
                    else:
                        raise Exception('Error: Unrecognised sampling method.')

                else:
                    # case of pixel completely outside the target originally recognised as limb-pixel
                    num_corners = 0
                    rec_corners = []
                    area_edges_in = []
                    area_edges = []
                    nav_cube[24][i] = 0  # Atotal/Ain

            else:
                raise Exception('Error: Invalid method.')

        rec_corners_all.append(rec_corners)
        area_edges_all.append(area_edges)
        area_edges_in_all.append(area_edges_in)

        longitude = []  # all the points of the edge of the pixel i
        latitude = []

        for i_new in range(0, num_corners):
            for j in range(0, num_points):
                radius, pclon, pclat = spice.spiceypy.reclat(area_edges[i_new][j])
                longitude.append(np.rad2deg(pclon))
                latitude.append(np.rad2deg(pclat))
        longitude_all.append(longitude)  # all the pixels
        latitude_all.append(latitude)

    return nav_cube, longitude_all, latitude_all, rec_corners_all, area_edges_in_all, output, area_edges_all