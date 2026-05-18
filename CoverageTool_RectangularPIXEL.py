"""
Created on November, 2021

@author: Christos Ntinos

This Module takes as input 1) the time in ET format, 2) the name of the target, 3) the name of the instrument, 4) the
name of the spacecraft, 5) the number of points for each of the edges of the FOV, 6) the aberration correction e.g.
"NONE", "LT+S", 7)  a matrix with the poistion of the corners of the FOV and 8) the computational method ("ELLIPSOID" or
"DSK/UNPRIORITIZED") and ouputs a three-column vector with rectangular coordinates of the projection of a rectangular FOV in
the IAU reference frame of the target.

Version: 1.0
"""

import spiceypy as spice
import numpy as np
import pixel_contains_poles


def main(time_et, target, instrument, spacecraft, num_points, aber_cor, dvec_all, method):
    """
    Parameters
    ----------
    time_et : float
        time of pixel acquisition
    target : str
        name of target
    instrument : str
        name of instrument
    spacecraft : str
        name of spacecraft
    num_points : int
        number of points per pixel-polygon edge
    aber_cor : str
        aberration correction flag
    dvec_all : ndarray
        ray's direction vector for center and corners
    method : str

    Returns
    -------
    area_corners : ndarray
        corners of pixel-polygon in rectangular coordinates (km)
    num_corners : int
        number of corners of pixel-polygon
    output_all : ndarray
        all points of the edges of pixel-polygon in rectangular coordinates (km)
    area_edges : ndarray
        points comprising the edges of the the pixel-polygon in rectangular coordinates (km)
    rec_corners : ndarray
        corners of pixel-polygon in rectangular coordinates (km)
    area_edges_in :
        points that define the edges of the part of the target inside the pixel (as projected on the limb plane)
        constructed py pixel corners and limb segments in rectangular coordinates (km)
    limb_plane : SPICE plane
        plane containing the limb
    area_edges_all : ndarray
        points that define the edge of the pixel-polygon in rectangular coordinates (km)
    north_in : bool
        North Pole contained in pixel (1) or not (0)
    south_in : bool
        South Pole contained in pixel (1) or not (0)

    This function calculates the projection on the target of a pixel for the different case of pixel-limb intersection.
    """
    start_index = 1  # the loops starts from 1 so as not to produce an error when aber_cor="NONE", else from 0

    # calculate position of spacecraft in the specified window
    pos_spacecraft, litime = spice.spiceypy.spkpos(spacecraft, time_et, f'IAU_{target}', aber_cor, target)

    # retrieve ellipsoid semi-major axes
    dim, values = spice.spiceypy.bodvrd(target, 'RADII', 3)
    radii_a = values[0]
    radii_b = values[1]
    radii_c = values[2]

    # create ellipsoid from parametric equation
    u, v = np.mgrid[0:2 * np.pi:100j, 0:np.pi:100j]
    xe = radii_a * np.cos(u) * np.sin(v)
    ye = radii_b * np.sin(u) * np.sin(v)
    ze = radii_c * np.cos(v)

    # calculate sub-observer point DO WE NEED THIS?
    subobsrv, trgepcsub, srfvecsub = spice.spiceypy.subpnt('INTERCEPT/ELLIPSOID', target, time_et,  # method=??
                                                           f'IAU_{target}', aber_cor, spacecraft)

    # create limb ellipse
    ellipse = spice.spiceypy.edlimb(radii_a, radii_b, radii_c, pos_spacecraft)
    limb_center, limb_major, limb_minor = spice.spiceypy.el2cgv(ellipse)
    ellipse_limb = np.zeros((3, 100))
    i = 0
    for t in np.linspace(0, 2 * np.pi, 100):
        ellipse_limb[0][i] = limb_center[0] + np.cos(t) * limb_major[0] + np.sin(t) * limb_minor[0]
        ellipse_limb[1][i] = limb_center[1] + np.cos(t) * limb_major[1] + np.sin(t) * limb_minor[1]
        ellipse_limb[2][i] = limb_center[2] + np.cos(t) * limb_major[2] + np.sin(t) * limb_minor[2]

        i = i + 1

    # calculate limb plane
    normal_vector_limb = np.cross(limb_major, limb_minor) / np.linalg.norm(np.cross(limb_major, limb_minor))
    constant = np.inner(ellipse_limb.T[0], normal_vector_limb)
    limb_plane = spice.spiceypy.nvc2pl(normal_vector_limb, constant)

    ####################################################################################################################

    # calculate limb points
    poscaswrttit, lt = spice.spiceypy.spkpos('CASSINI', time_et, f'IAU_{target}', aber_cor, target)
    possunwrttit, lt = spice.spiceypy.spkpos('SUN', time_et, f'IAU_{target}', aber_cor, target)
    vec1 = np.cross(poscaswrttit, possunwrttit)
    refvec = np.cross(poscaswrttit, vec1)
    ncuts = 50
    rolstp = 2 * np.pi / ncuts
    schstp = 4.0
    soltol = 1E-7
    maxn = ncuts

    if np.abs(np.arccos(
            np.dot(refvec, possunwrttit) / (np.linalg.norm(refvec) * np.linalg.norm(possunwrttit)))) > np.pi / 2:
        refvec = np.cross(vec1, poscaswrttit)
    else:
        refvec = np.cross(poscaswrttit, vec1)

    npts, limb_points, epochs, limb_tangts = spice.spiceypy.limbpt('TANGENT/ELLIPSOID', target, time_et,
                                                                   f'IAU_{target}', aber_cor, 'ELLIPSOID LIMB',
                                                                   spacecraft,
                                                                   refvec, rolstp, ncuts, schstp, soltol, maxn)
    # (0, 0, 1.0), np.pi, 2, 1.0, 1.0, 10)

    # calculate elliptical limb as the intersection of the ellipsoid and the limb plane
    ellipse = spice.spiceypy.inedpl(radii_a, radii_b, radii_c, limb_plane)

    # calculate elliptic limb parameters
    limb_center, limb_major, limb_minor = spice.spiceypy.el2cgv(ellipse)

    # create limb plane
    normal_vector_limb = np.cross(limb_major, limb_minor) / np.linalg.norm(np.cross(limb_major, limb_minor))
    constant = np.inner(limb_points[0], normal_vector_limb)
    limb_plane = spice.spiceypy.nvc2pl(normal_vector_limb, constant)

    #############################################################

    radii = np.linalg.norm(limb_major)
    radii_limb = np.linalg.norm(limb_major)
    # create circle in 3D representing the terminator
    trgeps, obspos, terminator_points = spice.spiceypy.edterm('UMBRAL', 'SUN', target, time_et,
                                                              f'IAU_{target}', aber_cor, spacecraft, 3)

    normal_vector_ter = np.cross(terminator_points[0], terminator_points[1]) / np.linalg.norm(
        np.cross(terminator_points[0], terminator_points[1]))
    unit_vector_ter = terminator_points[0] / np.linalg.norm(terminator_points[0])
    cross_vector_ter = np.cross(normal_vector_ter, unit_vector_ter)

    circle_terminator = np.zeros((3, num_points))
    i = 0
    for t in np.linspace(0, 2 * np.pi, num_points):
        circle_terminator[0][i] = radii * np.cos(t) * unit_vector_ter[0] + radii * np.sin(t) * cross_vector_ter[0]
        circle_terminator[1][i] = radii * np.cos(t) * unit_vector_ter[1] + radii * np.sin(t) * cross_vector_ter[1]
        circle_terminator[2][i] = radii * np.cos(t) * unit_vector_ter[2] + radii * np.sin(t) * cross_vector_ter[2]

        i = i + 1

    # retrieve instrument ID
    instrument_id = spice.spiceypy.bodn2c(instrument)

    # retrieve FOV characteristics (in the reference frame of the instrument)
    shape, frame, bsight, n, bounds = spice.spiceypy.getfov(instrument_id, 10)

    # calculate corners of fov DO I NEED THIS?
    width = np.maximum(spice.spiceypy.vsep(bounds[0], bounds[1]), spice.spiceypy.vsep(bounds[1], bounds[2])) / 2
    height = np.minimum(spice.spiceypy.vsep(bounds[0], bounds[1]), spice.spiceypy.vsep(bounds[1], bounds[2])) / 2

    fov_corner1 = np.matmul(spice.spiceypy.pxform(frame, f'IAU_{target}', time_et),
                            np.transpose(np.array(dvec_all[0])))
    fov_corner2 = np.matmul(spice.spiceypy.pxform(frame, f'IAU_{target}', time_et),
                            np.transpose(np.array(dvec_all[1])))
    fov_corner3 = np.matmul(spice.spiceypy.pxform(frame, f'IAU_{target}', time_et),
                            np.transpose(np.array(dvec_all[2])))
    fov_corner4 = np.matmul(spice.spiceypy.pxform(frame, f'IAU_{target}', time_et),
                            np.transpose(np.array(dvec_all[3])))

    # calculate intercepts of the corners of the fov with the limb plane
    npxpts1, xpt1 = spice.spiceypy.inrypl(pos_spacecraft, fov_corner1, limb_plane)
    npxpts2, xpt2 = spice.spiceypy.inrypl(pos_spacecraft, fov_corner2, limb_plane)
    npxpts3, xpt3 = spice.spiceypy.inrypl(pos_spacecraft, fov_corner3, limb_plane)
    npxpts4, xpt4 = spice.spiceypy.inrypl(pos_spacecraft, fov_corner4, limb_plane)

    # create one vector for the corner points for the fov on the limb plane
    rec_corners = np.array([xpt1, xpt2, xpt3, xpt4])

    # check that the fov intersects the limb plane
    if rec_corners.all() == 0:
        print('The FOV does not intersect the limb plane')
    else:
        # calculate the intersection point of the ellipsoid and the rays of the bounds of the rectangle
        rec_corners_e = np.zeros((4, 3))
        for i in range(0, 4):
            try:
                '''
                rec_corners_e[i] = spice.spiceypy.surfpt(pos_spacecraft, rec_corners[i] - pos_spacecraft,
                                                         radii_a, radii_b, radii_c)
                '''
                rec_corners_e[i], trgepc, srfvec = spice.spiceypy.sincpt(method, target, time_et,
                                                                            f'IAU_{target}', aber_cor, spacecraft,
                                                                            frame, dvec_all[i])
            except:
                pass

        flag1 = np.linalg.norm(rec_corners_e, axis=1)

        # calculate if the target is in the FOV of the instrument
        visible = spice.spiceypy.fovtrg(instrument, target, 'ELLIPSOID', f'IAU_{target}', aber_cor, spacecraft, time_et)

        # calculate the intersection point of the projection of the FOV on the limb plane with the limb ellipse
        rec_corners_e_limb = np.zeros((8, 3))

        if visible is True:
            for i in [0, 1, 2, -1]:
                try:
                    rec_corners_e_limb[2 * i][:] = spice.spiceypy.surfpt(rec_corners[i], rec_corners[i - 1] -
                                                                         rec_corners[i], radii_a, radii_b, radii_c)
                except:
                    pass

            for i in [0, 1, 2, -1]:
                try:
                    rec_corners_e_limb[2 * i + 1][:] = spice.spiceypy.surfpt(rec_corners[i],
                                                                             rec_corners[i + 1] - rec_corners[i],
                                                                             radii_a, radii_b, radii_c)
                except:
                    pass

        # remove intersection points that are not on the rectangle representing the FOV
        flag0 = np.zeros(8)
        flag0[1] = (rec_corners_e_limb[1][0] - rec_corners[0][0]) / (rec_corners[1][0] - rec_corners[0][0])
        flag0[2] = (rec_corners_e_limb[2][0] - rec_corners[0][0]) / (rec_corners[1][0] - rec_corners[0][0])
        flag0[3] = (rec_corners_e_limb[3][0] - rec_corners[1][0]) / (rec_corners[2][0] - rec_corners[1][0])
        flag0[4] = (rec_corners_e_limb[4][0] - rec_corners[1][0]) / (rec_corners[2][0] - rec_corners[1][0])
        flag0[5] = (rec_corners_e_limb[5][0] - rec_corners[2][0]) / (rec_corners[3][0] - rec_corners[2][0])
        flag0[6] = (rec_corners_e_limb[6][0] - rec_corners[2][0]) / (rec_corners[3][0] - rec_corners[2][0])
        flag0[7] = (rec_corners_e_limb[7][0] - rec_corners[3][0]) / (rec_corners[0][0] - rec_corners[3][0])
        flag0[0] = (rec_corners_e_limb[0][0] - rec_corners[3][0]) / (rec_corners[0][0] - rec_corners[3][0])

        for i in range(0, 8):
            if flag0[i] < 0 or flag0[i] > 1:
                rec_corners_e_limb[i][:] = np.zeros(3)

        # remove double intersection points
        for i in [0, 1, 2, 3, 4, 5, 6, -1]:
            if np.linalg.norm(rec_corners_e_limb[i]) != 0 and (
                    np.round(rec_corners_e_limb[i], 3) == np.round(rec_corners_e_limb[i + 1], 3)).all():
                rec_corners_e_limb[i + 1][:] = np.zeros(3)

        flag2 = np.linalg.norm(rec_corners_e_limb, axis=1)

        num_corners = 4  # DO I NEED THIS?
        area_corners = np.zeros((num_corners, 3))
        area_edges = np.zeros((num_corners, num_points, 3))
        # ----------------------------------------------------------------------------------------------------------------------
        # calculate the corners of the projection of the fov on the ellipsoid for all the different cases
        #print(flag1)
        #print(flag2)
        # case 1
        if flag1[0] == 0 and flag1[1] != 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] == 0 and flag2[1] != 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] == 0 and flag2[5] == 0 and flag2[6] == 0 and flag2[7] == 0:

            # print('case 1')
            num_corners = 3
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[1]
            area_corners[1] = rec_corners_e[1]
            area_corners[2] = rec_corners_e_limb[3]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[2] + i / num_points * (area_corners[0] - area_corners[2]) - limb_center
                area_edges[2][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

            # define edge of part of target inside the pixel (as projected on the limb plane)
            area_edges_in = np.concatenate((area_edges[2], [area_corners[0]], [rec_corners[1]], [area_edges[2][0]]), axis=0)

        # case 2
        elif flag1[0] == 0 and flag1[1] != 0 and flag1[2] != 0 and flag1[3] == 0 and flag2[0] == 0 and flag2[1] != 0 and \
                flag2[2] == 0 and flag2[3] == 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] == 0 and flag2[7] == 0:

            # print('case 2')
            num_corners = 4
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[1]
            area_corners[1] = rec_corners_e[1]
            area_corners[2] = rec_corners_e[2]
            area_corners[3] = rec_corners_e_limb[5]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[3] + i / num_points * (area_corners[0] - area_corners[3]) - limb_center
                area_edges[3][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

            area_edges_in = np.append(area_edges[3], [area_corners[0], rec_corners[1], rec_corners[2],
                                                      area_edges[3][0]], axis=0)

            area_edges_in = np.concatenate((area_edges[3], [area_corners[0]], [rec_corners[1]], [rec_corners[2]],
                                            [area_edges[3][0]]), axis=0)

        # case 3
        if flag1[0] == 0 and flag1[1] == 0 and flag1[2] != 0 and flag1[3] == 0 and flag2[0] == 0 and flag2[1] == 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] == 0 and flag2[7] == 0:

            # print('case 3')
            num_corners = 3
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[3]
            area_corners[1] = rec_corners_e[2]
            area_corners[2] = rec_corners_e_limb[5]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[2] + i / num_points * (area_corners[0] - area_corners[2]) - limb_center
                area_edges[2][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[2], [area_corners[0]], [rec_corners[2]],
                                            [area_edges[2][0]]), axis=0)

        # case 4
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] != 0 and flag1[3] != 0 and flag2[0] == 0 and flag2[1] == 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] == 0 and flag2[5] == 0 and flag2[6] == 0 and flag2[7] != 0:

            # print('case 4')
            num_corners = 4
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[7]
            area_corners[1] = rec_corners_e_limb[3]
            area_corners[2] = rec_corners_e[2]
            area_corners[3] = rec_corners_e[3]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]

            for i in range(start_index, num_points):
                ray = area_corners[0] + i / num_points * (area_corners[1] - area_corners[0]) - limb_center
                area_edges[0][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[0] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[0], [area_corners[1]], [rec_corners[2]], [rec_corners[3]],
                                            [area_edges[0][0]]), axis=0)

        # case 5
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] != 0 and flag2[0] == 0 and flag2[1] == 0 and \
                flag2[2] == 0 and flag2[3] == 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] == 0 and flag2[7] != 0:

            # print('case 5')
            num_corners = 3
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[7]
            area_corners[1] = rec_corners_e_limb[5]
            area_corners[2] = rec_corners_e[3]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]

            for i in range(start_index, num_points):
                ray = area_corners[0] + i / num_points * (area_corners[1] - area_corners[0]) - limb_center
                area_edges[0][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[0] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[0], [area_corners[1]], [rec_corners[3]],
                                            [area_edges[0][0]]), axis=0)

        # case 6
        elif flag1[0] != 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] != 0 and flag2[0] == 0 and flag2[1] != 0 and \
                flag2[2] == 0 and flag2[3] == 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] == 0 and flag2[7] == 0:

            # print('case 6')
            num_corners = 4
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e[0]
            area_corners[1] = rec_corners_e_limb[1]
            area_corners[2] = rec_corners_e_limb[5]
            area_corners[3] = rec_corners_e[3]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[1] + i / num_points * (area_corners[2] - area_corners[1]) - limb_center
                area_edges[1][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[0] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[1], [area_corners[2]], [rec_corners[3]], [rec_corners[0]],
                                            [area_edges[1][0]]), axis=0)

        # case 7
        elif flag1[0] != 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] == 0 and flag2[1] != 0 and \
                flag2[2] == 0 and flag2[3] == 0 and flag2[4] == 0 and flag2[5] == 0 and flag2[6] == 0 and flag2[7] != 0:

            # print('case 7')
            num_corners = 3
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e[0]
            area_corners[1] = rec_corners_e_limb[1]
            area_corners[2] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[1] + i / num_points * (area_corners[2] - area_corners[1]) - limb_center
                area_edges[1][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[0] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[1], [area_corners[2]], [rec_corners[0]],
                                            [area_edges[1][0]]), axis=0)

        # case 8
        elif flag1[0] != 0 and flag1[1] != 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] == 0 and flag2[1] == 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] == 0 and flag2[5] == 0 and flag2[6] == 0 and flag2[7] != 0:

            # print('case 8')
            num_corners = 4
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e[0]
            area_corners[1] = rec_corners_e[1]
            area_corners[2] = rec_corners_e_limb[3]
            area_corners[3] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[2] + i / num_points * (area_corners[3] - area_corners[2]) - limb_center
                area_edges[2][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[0] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[2], [area_corners[3]], [rec_corners[0]], [rec_corners[1]],
                                            [area_edges[2][0]]), axis=0)

        # case 9
        elif flag1[0] != 0 and flag1[1] != 0 and flag1[2] != 0 and flag1[3] != 0:

            # print('case 9')
            num_corners = 4
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e[0]
            area_corners[1] = rec_corners_e[1]
            area_corners[2] = rec_corners_e[2]
            area_corners[3] = rec_corners_e[3]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            '''
            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[0] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)
            '''

            for i in range(start_index, num_points):
                try:
                    area_edges[0][i][:], trgepc, srfvec = spice.spiceypy.sincpt(method, target, time_et,
                                                                            f'IAU_{target}', aber_cor, spacecraft,
                                                                            frame, dvec_all[0] + i / num_points * (
                                                                            dvec_all[1] - dvec_all[0]))
                except:
                    area_edges[0][i][0] = np.nan
                    area_edges[0][i][1] = np.nan
                    area_edges[0][i][2] = np.nan
                    print('not found')

                try:
                    area_edges[1][i][:], trgepc, srfvec = spice.spiceypy.sincpt(method, target, time_et,
                                                                            f'IAU_{target}', aber_cor, spacecraft,
                                                                            frame, dvec_all[1] + i / num_points * (
                                                                                    dvec_all[2] - dvec_all[1]))
                except:
                    area_edges[1][i][0] = np.nan
                    area_edges[1][i][1] = np.nan
                    area_edges[1][i][2] = np.nan
                    print('not found')

                try:
                    area_edges[2][i][:], trgepc, srfvec = spice.spiceypy.sincpt(method, target, time_et,
                                                                            f'IAU_{target}', aber_cor, spacecraft,
                                                                            frame, dvec_all[2] + i / num_points * (
                                                                                    dvec_all[3] - dvec_all[2]))
                except:
                    area_edges[2][i][0] = np.nan
                    area_edges[2][i][1] = np.nan
                    area_edges[2][i][2] = np.nan
                    print('not found')

                try:
                    area_edges[3][i][:], trgepc, srfvec = spice.spiceypy.sincpt(method, target, time_et,
                                                                            f'IAU_{target}', aber_cor, spacecraft,
                                                                            frame, dvec_all[3] + i / num_points * (
                                                                                    dvec_all[0] - dvec_all[3]))
                except:
                    area_edges[3][i][0] = np.nan
                    area_edges[3][i][1] = np.nan
                    area_edges[3][i][2] = np.nan
                    print('not found')

            area_edges_in = []  # Ain/Atotal=1

        '''
        # case 10
        if (flag1 == 0).all() and (flag2 == 0).all() and visible is True:
            num_corners = 2
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = ellipse_limb.T[0]
            area_corners[1] = ellipse_limb.T[int(num_points / 2)]

            unit_vector_limb = (area_corners[1] - limb_center) / np.linalg.norm((area_corners[1] - limb_center))
            cross_vector_limb = np.cross(normal_vector_limb, unit_vector_limb)
            ang_sep = spice.spiceypy.vsep(area_corners[1] - limb_center, area_corners[0] - limb_center)
            i = 0
            for t in np.linspace(0, ang_sep, num_points):
                area_edges[0][i][:] = radii_limb * np.cos(t) * unit_vector_limb + radii_limb * np.sin(
                    t) * cross_vector_limb + limb_center
                i = i + 1

            unit_vector_limb = (area_corners[0] - limb_center) / np.linalg.norm((area_corners[0] - limb_center))
            cross_vector_limb = np.cross(normal_vector_limb, unit_vector_limb)
            ang_sep = spice.spiceypy.vsep(area_corners[1] - limb_center, area_corners[0] - limb_center)
            i = 0
            for t in np.linspace(0, ang_sep, num_points):
                area_edges[1][i][:] = radii_limb * np.cos(t) * unit_vector_limb + radii_limb * np.sin(
                    t) * cross_vector_limb + limb_center
                i = i + 1

            area_edges_in = []
        '''

        # case 11
        if (flag1 == 0).all() and (flag2 == 0).all():  # and visible is False: !!! elif???

            # print('case 11')
            num_corners = 0
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            notfound_text = '{} is not in the FOV of {}'.format(target, instrument)

            area_edges_in = []

        # case 12
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] != 0 and flag2[1] != 0 and \
                flag2[2] != 0 and flag2[3] != 0 and flag2[4] != 0 and flag2[5] != 0 and flag2[6] != 0 and flag2[7] != 0:

            # print('case 12')
            num_corners = 8
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[0]
            area_corners[1] = rec_corners_e_limb[1]
            area_corners[2] = rec_corners_e_limb[2]
            area_corners[3] = rec_corners_e_limb[3]
            area_corners[4] = rec_corners_e_limb[4]
            area_corners[5] = rec_corners_e_limb[5]
            area_corners[6] = rec_corners_e_limb[6]
            area_corners[7] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]
            area_edges[5][0] = area_corners[5]
            area_edges[6][0] = area_corners[6]
            area_edges[7][0] = area_corners[7]

            for i in range(start_index, num_points):
                ray = area_corners[0] + i / num_points * (area_corners[1] - area_corners[0]) - limb_center
                area_edges[0][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[2] + i / num_points * (area_corners[3] - area_corners[2]) - limb_center
                area_edges[2][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[4] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[4] + i / num_points * (area_corners[5] - area_corners[4]) - limb_center
                area_edges[4][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[5][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[5] + i / num_points * (
                        area_corners[6] - area_corners[5]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[6] + i / num_points * (area_corners[7] - area_corners[6]) - limb_center
                area_edges[6][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[7][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[7] + i / num_points * (
                        area_corners[0] - area_corners[7]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[0], [area_corners[1]], area_edges[2], [area_corners[3]],
                                            area_edges[4], [area_corners[5]], area_edges[6], [area_corners[7]],
                                            [area_edges[0][0]]), axis=0)

        # case 13
        elif flag1[0] != 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] != 0 and flag2[0] == 0 and flag2[1] != 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] != 0 and flag2[5] != 0 and flag2[6] == 0 and flag2[7] == 0:

            # print('case 13')
            num_corners = 6
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e[0]
            area_corners[1] = rec_corners_e_limb[1]
            area_corners[2] = rec_corners_e_limb[3]
            area_corners[3] = rec_corners_e_limb[4]
            area_corners[4] = rec_corners_e_limb[5]
            area_corners[5] = rec_corners_e[3]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]
            area_edges[5][0] = area_corners[5]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[1] + i / num_points * (area_corners[2] - area_corners[1]) - limb_center
                area_edges[1][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[3] + i / num_points * (area_corners[4] - area_corners[3]) - limb_center
                area_edges[3][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[4][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[4] + i / num_points * (
                        area_corners[5] - area_corners[4]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[5][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[5] + i / num_points * (
                        area_corners[0] - area_corners[5]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate(([rec_corners[0]], area_edges[1], [area_corners[2]], area_edges[3],
                                            [area_corners[4]], [rec_corners[3]], [rec_corners[0]]), axis=0)

        # case 14
        elif flag1[0] != 0 and flag1[1] != 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] == 0 and flag2[1] == 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] != 0 and flag2[7] != 0:

            # print('case 14')
            num_corners = 6
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e[0]
            area_corners[1] = rec_corners_e[1]
            area_corners[2] = rec_corners_e_limb[3]
            area_corners[3] = rec_corners_e_limb[5]
            area_corners[4] = rec_corners_e_limb[6]
            area_corners[5] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]
            area_edges[5][0] = area_corners[5]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[2] + i / num_points * (area_corners[3] - area_corners[2]) - limb_center
                area_edges[2][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[4] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[4] + i / num_points * (area_corners[5] - area_corners[4]) - limb_center
                area_edges[4][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[5][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[5] + i / num_points * (
                        area_corners[0] - area_corners[5]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate(([rec_corners[0]], [rec_corners[1]], area_edges[2], [area_corners[3]],
                                            area_edges[4], [area_corners[5]], [rec_corners[0]]), axis=0)

        # case 15
        elif flag1[0] == 0 and flag1[1] != 0 and flag1[2] != 0 and flag1[3] == 0 and flag2[0] != 0 and flag2[1] != 0 and \
                flag2[2] == 0 and flag2[3] == 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] == 0 and flag2[7] != 0:

            # print('case 15')
            num_corners = 6
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[0]
            area_corners[1] = rec_corners_e_limb[1]
            area_corners[2] = rec_corners_e[1]
            area_corners[3] = rec_corners_e[2]
            area_corners[4] = rec_corners_e_limb[5]
            area_corners[5] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]
            area_edges[5][0] = area_corners[5]

            for i in range(start_index, num_points):
                ray = area_corners[0] + i / num_points * (area_corners[1] - area_corners[0]) - limb_center
                area_edges[0][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[4] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[4] + i / num_points * (area_corners[5] - area_corners[4]) - limb_center
                area_edges[4][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[5][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[5] + i / num_points * (
                        area_corners[0] - area_corners[5]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[0], [area_corners[1]], [rec_corners[1]], [rec_corners[2]],
                                            area_edges[4], [area_corners[5]], [area_edges[0][0]]), axis=0)

        # case 16
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] != 0 and flag1[3] != 0 and flag2[0] == 0 and flag2[1] != 0 and \
                flag2[2] != 0 and flag2[3] != 0 and flag2[4] == 0 and flag2[5] == 0 and flag2[6] == 0 and flag2[7] != 0:

            # print('case 16')
            num_corners = 6
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[7]
            area_corners[1] = rec_corners_e_limb[1]
            area_corners[2] = rec_corners_e_limb[2]
            area_corners[3] = rec_corners_e_limb[3]
            area_corners[4] = rec_corners_e[2]
            area_corners[5] = rec_corners_e[3]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]
            area_edges[5][0] = area_corners[5]

            for i in range(start_index, num_points):
                ray = area_corners[0] + i / num_points * (area_corners[1] - area_corners[0]) - limb_center
                area_edges[0][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[2] + i / num_points * (area_corners[3] - area_corners[2]) - limb_center
                area_edges[2][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[4] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[4][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[4] + i / num_points * (
                        area_corners[5] - area_corners[4]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[5][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[5] + i / num_points * (
                        area_corners[0] - area_corners[5]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[0], [area_corners[1]], area_edges[2], [area_corners[3]],
                                            [rec_corners[2]], [rec_corners[3]], [area_edges[0][0]]), axis=0)

        # case 17
        elif flag1[0] != 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] == 0 and flag2[1] != 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] != 0 and flag2[5] != 0 and flag2[6] != 0 and flag2[7] != 0:

            # print('case 17')
            num_corners = 7
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e[0]
            area_corners[1] = rec_corners_e_limb[1]
            area_corners[2] = rec_corners_e_limb[3]
            area_corners[3] = rec_corners_e_limb[4]
            area_corners[4] = rec_corners_e_limb[5]
            area_corners[5] = rec_corners_e_limb[6]
            area_corners[6] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]
            area_edges[5][0] = area_corners[5]
            area_edges[6][0] = area_corners[6]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[1] + i / num_points * (area_corners[2] - area_corners[1]) - limb_center
                area_edges[1][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[3] + i / num_points * (area_corners[4] - area_corners[3]) - limb_center
                area_edges[3][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[4][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[4] + i / num_points * (
                        area_corners[5] - area_corners[4]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[5] + i / num_points * (area_corners[6] - area_corners[5]) - limb_center
                area_edges[5][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[6][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[6] + i / num_points * (
                        area_corners[0] - area_corners[6]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate(([rec_corners[0]], area_edges[1], [area_corners[2]], area_edges[3],
                                            [area_corners[4]], area_edges[5], [area_corners[6]], [rec_corners[0]]), axis=0)

        # case 18
        elif flag1[0] == 0 and flag1[1] != 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] != 0 and flag2[1] != 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] != 0 and flag2[7] != 0:

            # print('case 18')
            num_corners = 7
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[1]
            area_corners[1] = rec_corners_e[1]
            area_corners[2] = rec_corners_e_limb[3]
            area_corners[3] = rec_corners_e_limb[5]
            area_corners[4] = rec_corners_e_limb[6]
            area_corners[5] = rec_corners_e_limb[7]
            area_corners[6] = rec_corners_e_limb[0]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]
            area_edges[5][0] = area_corners[5]
            area_edges[6][0] = area_corners[6]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[2] + i / num_points * (area_corners[3] - area_corners[2]) - limb_center
                area_edges[2][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[4] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[4] + i / num_points * (area_corners[5] - area_corners[4]) - limb_center
                area_edges[4][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[5][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[5] + i / num_points * (
                        area_corners[6] - area_corners[5]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[6] + i / num_points * (area_corners[0] - area_corners[6]) - limb_center
                area_edges[6][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate(([area_corners[0]], [rec_corners[1]], area_edges[2], [area_corners[3]],
                                            area_edges[4], [area_corners[5]], area_edges[6], [area_corners[0]]), axis=0)

        # case 19
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] != 0 and flag1[3] == 0 and flag2[0] != 0 and flag2[1] != 0 and \
                flag2[2] != 0 and flag2[3] != 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] == 0 and flag2[7] != 0:

            # print('case 19')
            num_corners = 7
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[1]
            area_corners[1] = rec_corners_e_limb[2]
            area_corners[2] = rec_corners_e_limb[3]
            area_corners[3] = rec_corners_e[2]
            area_corners[4] = rec_corners_e_limb[5]
            area_corners[5] = rec_corners_e_limb[7]
            area_corners[6] = rec_corners_e_limb[0]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]
            area_edges[5][0] = area_corners[5]
            area_edges[6][0] = area_corners[6]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[1] + i / num_points * (area_corners[2] - area_corners[1]) - limb_center
                area_edges[1][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[4] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[4] + i / num_points * (area_corners[5] - area_corners[4]) - limb_center
                area_edges[4][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[5][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[5] + i / num_points * (
                        area_corners[6] - area_corners[5]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[6] + i / num_points * (area_corners[0] - area_corners[6]) - limb_center
                area_edges[6][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate(([area_corners[0]], area_edges[1], [area_corners[2]], [rec_corners[2]],
                                            area_edges[4], [area_corners[5]], area_edges[6], [area_corners[0]]), axis=0)

        # case 20
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] != 0 and flag2[0] != 0 and flag2[1] != 0 and \
                flag2[2] != 0 and flag2[3] != 0 and flag2[4] != 0 and flag2[5] != 0 and flag2[6] == 0 and flag2[7] == 0:

            # print('case 20')
            num_corners = 7
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[1]
            area_corners[1] = rec_corners_e_limb[2]
            area_corners[2] = rec_corners_e_limb[3]
            area_corners[3] = rec_corners_e_limb[4]
            area_corners[4] = rec_corners_e_limb[5]
            area_corners[5] = rec_corners_e[3]
            area_corners[6] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]
            area_edges[5][0] = area_corners[5]
            area_edges[6][0] = area_corners[6]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[1] + i / num_points * (area_corners[2] - area_corners[1]) - limb_center
                area_edges[1][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[3] + i / num_points * (area_corners[4] - area_corners[3]) - limb_center
                area_edges[3][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[4][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[4] + i / num_points * (
                        area_corners[5] - area_corners[4]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[5][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[5] + i / num_points * (
                        area_corners[6] - area_corners[5]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[6] + i / num_points * (area_corners[0] - area_corners[6]) - limb_center
                area_edges[6][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate(([area_corners[0]], area_edges[1], [area_corners[2]], area_edges[3],
                                            [area_corners[4]], [rec_corners[3]], area_edges[6],
                                            [area_corners[0]]), axis=0)

        # case 21
        elif flag1[0] == 0 and flag1[1] != 0 and flag1[2] != 0 and flag1[3] != 0 and flag2[0] == 0 and flag2[1] != 0 and \
                flag2[2] == 0 and flag2[3] == 0 and flag2[4] == 0 and flag2[5] == 0 and flag2[6] == 0 and flag2[7] != 0:

            # print('case 21')
            num_corners = 5
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[1]
            area_corners[1] = rec_corners_e[1]
            area_corners[2] = rec_corners_e[2]
            area_corners[3] = rec_corners_e[3]
            area_corners[4] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[4] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[4] + i / num_points * (area_corners[0] - area_corners[4]) - limb_center
                area_edges[4][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate(([area_corners[0]], [rec_corners[1]], [rec_corners[2]], [rec_corners[3]],
                                            area_edges[4], [area_corners[3]]), axis=0)

        # case 22
        elif flag1[0] != 0 and flag1[1] == 0 and flag1[2] != 0 and flag1[3] != 0 and flag2[0] == 0 and flag2[1] != 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] == 0 and flag2[5] == 0 and flag2[6] == 0 and flag2[7] == 0:

            # print('case 22')
            num_corners = 5
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e[0]
            area_corners[1] = rec_corners_e_limb[1]
            area_corners[2] = rec_corners_e_limb[3]
            area_corners[3] = rec_corners_e[2]
            area_corners[4] = rec_corners_e[3]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[1] + i / num_points * (area_corners[2] - area_corners[1]) - limb_center
                area_edges[1][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[4] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[4][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[4] + i / num_points * (
                        area_corners[0] - area_corners[4]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate(([rec_corners[0]], area_edges[1], [area_corners[2]], [rec_corners[2]],
                                            [rec_corners[3]], [rec_corners[0]]), axis=0)

        # case 23
        elif flag1[0] != 0 and flag1[1] != 0 and flag1[2] == 0 and flag1[3] != 0 and flag2[0] == 0 and flag2[1] == 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] == 0 and flag2[7] == 0:

            # print('case 23')
            num_corners = 5
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e[0]
            area_corners[1] = rec_corners_e[1]
            area_corners[2] = rec_corners_e_limb[3]
            area_corners[3] = rec_corners_e_limb[5]
            area_corners[4] = rec_corners_e[3]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[2] + i / num_points * (area_corners[3] - area_corners[2]) - limb_center
                area_edges[2][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[4] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[4][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[4] + i / num_points * (
                        area_corners[0] - area_corners[4]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate(([rec_corners[0]], [rec_corners[1]], area_edges[2], [area_corners[3]],
                                            [rec_corners[3]], [rec_corners[0]]), axis=0)

        # case 24
        elif flag1[0] != 0 and flag1[1] != 0 and flag1[2] != 0 and flag1[3] == 0 and flag2[0] == 0 and flag2[1] == 0 and \
                flag2[2] == 0 and flag2[3] == 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] == 0 and flag2[7] != 0:

            # print('case 24')
            num_corners = 5
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e[0]
            area_corners[1] = rec_corners_e[1]
            area_corners[2] = rec_corners_e[2]
            area_corners[3] = rec_corners_e_limb[5]
            area_corners[4] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[3] + i / num_points * (area_corners[4] - area_corners[3]) - limb_center
                area_edges[3][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[4][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[4] + i / num_points * (
                        area_corners[0] - area_corners[4]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[0], area_edges[1], area_edges[2], [area_corners[3]],
                                            [area_corners[4]], [rec_corners[0]]), axis=0)

        # case 25
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] == 0 and flag2[1] != 0 and \
                flag2[2] != 0 and flag2[3] == 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] != 0 and flag2[7] == 0:

            # print('case 25')
            num_corners = 4
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[1]
            area_corners[1] = rec_corners_e_limb[2]
            area_corners[2] = rec_corners_e_limb[5]
            area_corners[3] = rec_corners_e_limb[6]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[1] + i / num_points * (area_corners[2] - area_corners[1]) - limb_center
                area_edges[1][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[3] + i / num_points * (area_corners[0] - area_corners[3]) - limb_center
                area_edges[3][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate(([area_corners[0]], area_edges[1], [area_corners[2]], area_edges[3],
                                            [area_corners[0]]), axis=0)

        # case 26
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] != 0 and flag2[1] == 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] != 0 and flag2[5] == 0 and flag2[6] == 0 and flag2[7] != 0:

            # print('case 26')
            num_corners = 4
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[0]
            area_corners[1] = rec_corners_e_limb[3]
            area_corners[2] = rec_corners_e_limb[4]
            area_corners[3] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]

            for i in range(start_index, num_points):
                ray = area_corners[0] + i / num_points * (area_corners[1] - area_corners[0]) - limb_center
                area_edges[0][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[2] + i / num_points * (area_corners[3] - area_corners[2]) - limb_center
                area_edges[2][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[0] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[0], [area_corners[1]], area_edges[2], [area_corners[3]],
                                            area_edges[0][0]), axis=0)

        # case 27
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] == 0 and flag2[1] != 0 and \
                flag2[2] != 0 and flag2[3] == 0 and flag2[4] == 0 and flag2[5] == 0 and flag2[6] == 0 and flag2[7] == 0:

            # print('case 27')
            num_corners = 2
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[1]
            area_corners[1] = rec_corners_e_limb[2]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[1] + i / num_points * (area_corners[0] - area_corners[1]) - limb_center
                area_edges[1][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[1], [area_corners[0]], [area_edges[1][0]]), axis=0)

        # case 28
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] == 0 and flag2[1] == 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] != 0 and flag2[5] == 0 and flag2[6] == 0 and flag2[7] == 0:

            # print('case 28')
            num_corners = 2
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[3]
            area_corners[1] = rec_corners_e_limb[4]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[1] + i / num_points * (area_corners[0] - area_corners[1]) - limb_center
                area_edges[1][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[1], [area_corners[0]], [area_edges[1][0]]), axis=0)

        # case 29
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] == 0 and flag2[1] == 0 and \
                flag2[2] == 0 and flag2[3] == 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] != 0 and flag2[7] == 0:

            # print('case 29')
            num_corners = 2
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[5]
            area_corners[1] = rec_corners_e_limb[6]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[1] + i / num_points * (area_corners[0] - area_corners[1]) - limb_center
                area_edges[1][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[1], [area_corners[0]], [area_edges[1][0]]), axis=0)

        # case 30
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] != 0 and flag2[1] == 0 and \
                flag2[2] == 0 and flag2[3] == 0 and flag2[4] == 0 and flag2[5] == 0 and flag2[6] == 0 and flag2[7] != 0:

            # print('case 30')
            num_corners = 2
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[7]
            area_corners[1] = rec_corners_e_limb[0]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[1] + i / num_points * (area_corners[0] - area_corners[1]) - limb_center
                area_edges[1][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[1], [area_corners[0]], [area_edges[1][0]]), axis=0)

        # case 31
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] != 0 and flag2[1] != 0 and \
                flag2[2] != 0 and flag2[3] == 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] != 0 and flag2[7] != 0:

            # print('case 31')
            num_corners = 6
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[0]
            area_corners[1] = rec_corners_e_limb[1]
            area_corners[2] = rec_corners_e_limb[2]
            area_corners[3] = rec_corners_e_limb[5]
            area_corners[4] = rec_corners_e_limb[6]
            area_corners[5] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]
            area_edges[5][0] = area_corners[5]

            for i in range(start_index, num_points):
                ray = area_corners[0] + i / num_points * (area_corners[1] - area_corners[0]) - limb_center
                area_edges[0][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[2] + i / num_points * (area_corners[3] - area_corners[2]) - limb_center
                area_edges[2][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[4] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[4] + i / num_points * (area_corners[5] - area_corners[4]) - limb_center
                area_edges[4][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[5][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[5] + i / num_points * (
                        area_corners[0] - area_corners[5]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[0], [area_corners[1]], area_edges[2], [area_corners[3]],
                                            area_edges[4], [area_corners[5]], [area_edges[0][0]]), axis=0)

        # case 32
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] != 0 and flag2[1] != 0 and \
                flag2[2] != 0 and flag2[3] != 0 and flag2[4] != 0 and flag2[5] == 0 and flag2[6] == 0 and flag2[7] != 0:

            # print('case 32')
            num_corners = 6
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[0]
            area_corners[1] = rec_corners_e_limb[1]
            area_corners[2] = rec_corners_e_limb[2]
            area_corners[3] = rec_corners_e_limb[3]
            area_corners[4] = rec_corners_e_limb[4]
            area_corners[5] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]
            area_edges[5][0] = area_corners[5]

            for i in range(start_index, num_points):
                ray = area_corners[0] + i / num_points * (area_corners[1] - area_corners[0]) - limb_center
                area_edges[0][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[2] + i / num_points * (area_corners[3] - area_corners[2]) - limb_center
                area_edges[2][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[4] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[4] + i / num_points * (area_corners[5] - area_corners[4]) - limb_center
                area_edges[4][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[5][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[5] + i / num_points * (
                        area_corners[0] - area_corners[5]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[0], [area_corners[1]], area_edges[2], [area_corners[3]],
                                            area_edges[2], [area_corners[5]], [area_edges[0][0]]), axis=0)

        # case 33
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] == 0 and flag2[1] != 0 and \
                flag2[2] != 0 and flag2[3] != 0 and flag2[4] != 0 and flag2[5] != 0 and flag2[6] != 0 and flag2[7] == 0:

            # print('case 33')
            num_corners = 6
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[6]
            area_corners[1] = rec_corners_e_limb[1]
            area_corners[2] = rec_corners_e_limb[2]
            area_corners[3] = rec_corners_e_limb[3]
            area_corners[4] = rec_corners_e_limb[4]
            area_corners[5] = rec_corners_e_limb[5]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]
            area_edges[5][0] = area_corners[5]

            for i in range(start_index, num_points):
                ray = area_corners[0] + i / num_points * (area_corners[1] - area_corners[0]) - limb_center
                area_edges[0][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[2] + i / num_points * (area_corners[3] - area_corners[2]) - limb_center
                area_edges[2][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[4] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[4] + i / num_points * (area_corners[5] - area_corners[4]) - limb_center
                area_edges[4][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[5][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[5] + i / num_points * (
                        area_corners[0] - area_corners[5]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[0], [area_corners[1]], area_edges[2], [area_corners[3]],
                                            area_edges[4], [area_corners[5]], area_edges[0][0]), axis=0)

        # case 34
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] != 0 and flag2[1] == 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] != 0 and flag2[5] != 0 and flag2[6] != 0 and flag2[7] != 0:

            # print('case 34')
            num_corners = 6
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[0]
            area_corners[1] = rec_corners_e_limb[3]
            area_corners[2] = rec_corners_e_limb[4]
            area_corners[3] = rec_corners_e_limb[5]
            area_corners[4] = rec_corners_e_limb[6]
            area_corners[5] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]
            area_edges[5][0] = area_corners[5]

            for i in range(start_index, num_points):
                ray = area_corners[0] + i / num_points * (area_corners[1] - area_corners[0]) - limb_center
                area_edges[0][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[2] + i / num_points * (area_corners[3] - area_corners[2]) - limb_center
                area_edges[2][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[4] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[4] + i / num_points * (area_corners[5] - area_corners[4]) - limb_center
                area_edges[4][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[5][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[5] + i / num_points * (
                        area_corners[0] - area_corners[5]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate((area_edges[0], [area_corners[1]], area_edges[2], [area_corners[3]],
                                            area_edges[4], [area_corners[5]], [area_edges[0][0]]), axis=0)

        # case 35
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] == 0 and flag2[1] != 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] == 0 and flag2[7] != 0:

            # print('case 35')
            num_corners = 2
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = ellipse_limb.T[0]
            area_corners[1] = ellipse_limb.T[int(num_points / 2)]

            area_edges[0][0] = area_corners[0]  # add first corner of edge because start_index=1
            area_edges[1][0] = area_corners[1]

            unit_vector_limb = (area_corners[1] - limb_center) / np.linalg.norm((area_corners[1] - limb_center))
            cross_vector_limb = np.cross(normal_vector_limb, unit_vector_limb)
            ang_sep = spice.spiceypy.vsep(area_corners[1] - limb_center, area_corners[0] - limb_center)
            i = 0
            for t in np.linspace(0, ang_sep, num_points):
                area_edges[0][i][:] = radii_limb * np.cos(t) * unit_vector_limb + radii_limb * np.sin(
                    t) * cross_vector_limb + limb_center
                i = i + 1

            unit_vector_limb = (area_corners[0] - limb_center) / np.linalg.norm((area_corners[0] - limb_center))
            cross_vector_limb = np.cross(normal_vector_limb, unit_vector_limb)
            ang_sep = spice.spiceypy.vsep(area_corners[1] - limb_center, area_corners[0] - limb_center)
            i = 0
            for t in np.linspace(0, ang_sep, num_points):
                area_edges[1][i][:] = radii_limb * np.cos(t) * unit_vector_limb + radii_limb * np.sin(
                    t) * cross_vector_limb + limb_center
                i = i + 1

            area_edges_in = ellipse_limb

        # case 36
        elif flag1[0] != 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] == 0 and flag2[1] != 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] != 0 and flag2[5] == 0 and flag2[6] == 0 and flag2[7] != 0:

            # print('case 36')
            num_corners = 5
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e[0]
            area_corners[1] = rec_corners_e_limb[1]
            area_corners[2] = rec_corners_e_limb[3]
            area_corners[3] = rec_corners_e_limb[4]
            area_corners[4] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[1] + i / num_points * (area_corners[2] - area_corners[1]) - limb_center
                area_edges[1][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[3] + i / num_points * (area_corners[4] - area_corners[3]) - limb_center
                area_edges[3][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[4][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[4] + i / num_points * (
                        area_corners[0] - area_corners[4]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate(([rec_corners[0]], area_edges[1], [area_corners[2]], area_edges[3],
                                            [area_corners[4]], [rec_corners[0]]), axis=0)

        # case 37
        elif flag1[0] != 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] == 0 and flag2[1] != 0 and \
                flag2[2] == 0 and flag2[3] == 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] != 0 and flag2[7] != 0:

            # print('case 37')
            num_corners = 5
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e[0]
            area_corners[1] = rec_corners_e_limb[1]
            area_corners[2] = rec_corners_e_limb[5]
            area_corners[3] = rec_corners_e_limb[6]
            area_corners[4] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[1] + i / num_points * (area_corners[2] - area_corners[1]) - limb_center
                area_edges[1][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[3] + i / num_points * (area_corners[4] - area_corners[3]) - limb_center
                area_edges[3][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[4][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[4] + i / num_points * (
                        area_corners[0] - area_corners[4]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate(([rec_corners[0]], area_edges[1], [area_corners[2]], area_edges[3],
                                            [area_corners[4]], [rec_corners[0]]), axis=0)

        # case 38
        elif flag1[0] == 0 and flag1[1] != 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] != 0 and flag2[1] != 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] == 0 and flag2[5] == 0 and flag2[6] == 0 and flag2[7] != 0:

            # print('case 38')
            num_corners = 5
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[0]
            area_corners[1] = rec_corners_e_limb[1]
            area_corners[2] = rec_corners_e[1]
            area_corners[3] = rec_corners_e_limb[3]
            area_corners[4] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]

            for i in range(start_index, num_points):
                ray = area_corners[0] + i / num_points * (area_corners[1] - area_corners[0]) - limb_center
                area_edges[0][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[3] + i / num_points * (area_corners[4] - area_corners[3]) - limb_center
                area_edges[3][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[4][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[4] + i / num_points * (
                        area_corners[0] - area_corners[4]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate(([rec_corners[1]], area_edges[3], [area_corners[4]], area_edges[0],
                                            [area_corners[1]], [rec_corners[1]]), axis=0)

        # case 39
        elif flag1[0] == 0 and flag1[1] != 0 and flag1[2] == 0 and flag1[3] == 0 and flag2[0] == 0 and flag2[1] != 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] != 0 and flag2[7] == 0:

            # print('case 39')
            num_corners = 5
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[1]
            area_corners[1] = rec_corners_e[1]
            area_corners[2] = rec_corners_e_limb[3]
            area_corners[3] = rec_corners_e_limb[5]
            area_corners[4] = rec_corners_e_limb[6]

            area_edges[0][0] = area_corners[0]
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[2] + i / num_points * (area_corners[3] - area_corners[2]) - limb_center
                area_edges[2][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[4] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[4] + i / num_points * (area_corners[0] - area_corners[4]) - limb_center
                area_edges[4][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate(([rec_corners[1]], area_edges[2], [area_corners[3]], area_edges[4],
                                            [area_corners[0]], [rec_corners[1]]), axis=0)

        # case 40
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] != 0 and flag1[3] == 0 and flag2[0] == 0 and flag2[1] != 0 and \
                flag2[2] != 0 and flag2[3] != 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] == 0 and flag2[7] == 0:

            # print('case 40')
            num_corners = 5
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[1]
            area_corners[1] = rec_corners_e_limb[2]
            area_corners[2] = rec_corners_e_limb[3]
            area_corners[3] = rec_corners_e[2]
            area_corners[4] = rec_corners_e_limb[5]

            area_edges[0][0] = area_corners[0]
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[1] + i / num_points * (area_corners[2] - area_corners[1]) - limb_center
                area_edges[1][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[4] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[4] + i / num_points * (area_corners[0] - area_corners[4]) - limb_center
                area_edges[4][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate(([rec_corners[2]], area_edges[4], [area_corners[0]], area_edges[1],
                                            [area_corners[2]], [rec_corners[2]]), axis=0)

        # case 41
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] != 0 and flag1[3] == 0 and flag2[0] != 0 and flag2[1] == 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] == 0 and flag2[7] != 0:

            # print('case 41')
            num_corners = 5
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[0]
            area_corners[1] = rec_corners_e_limb[3]
            area_corners[2] = rec_corners_e[2]
            area_corners[3] = rec_corners_e_limb[5]
            area_corners[4] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]

            for i in range(start_index, num_points):
                ray = area_corners[0] + i / num_points * (area_corners[1] - area_corners[0]) - limb_center
                area_edges[0][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[1][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[1] + i / num_points * (
                        area_corners[2] - area_corners[1]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[3] + i / num_points * (area_corners[4] - area_corners[3]) - limb_center
                area_edges[3][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[4][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[4] + i / num_points * (
                        area_corners[0] - area_corners[4]) - pos_spacecraft, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate(([rec_corners[2]], area_edges[3], [area_corners[4]], area_edges[0],
                                            [area_corners[1]], [rec_corners[2]]), axis=0)

        # case 42
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] != 0 and flag2[0] == 0 and flag2[1] != 0 and \
                flag2[2] != 0 and flag2[3] == 0 and flag2[4] == 0 and flag2[5] != 0 and flag2[6] == 0 and flag2[7] != 0:

            # print('case 42')
            num_corners = 5
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[1]
            area_corners[1] = rec_corners_e_limb[2]
            area_corners[2] = rec_corners_e_limb[5]
            area_corners[3] = rec_corners_e[3]
            area_corners[4] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[1] + i / num_points * (area_corners[2] - area_corners[1]) - limb_center
                area_edges[1][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[4] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[4] + i / num_points * (area_corners[0] - area_corners[4]) - limb_center
                area_edges[4][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate(([rec_corners[3]], area_edges[4], [area_corners[0]], area_edges[1],
                                            [area_corners[2]], [rec_corners[3]]), axis=0)

        # case 43
        elif flag1[0] == 0 and flag1[1] == 0 and flag1[2] == 0 and flag1[3] != 0 and flag2[0] == 0 and flag2[1] == 0 and \
                flag2[2] == 0 and flag2[3] != 0 and flag2[4] != 0 and flag2[5] != 0 and flag2[6] == 0 and flag2[7] != 0:

            # print('case 43')
            num_corners = 5
            area_corners = np.zeros((num_corners, 3))
            area_edges = np.zeros((num_corners, num_points, 3))

            area_corners[0] = rec_corners_e_limb[3]
            area_corners[1] = rec_corners_e_limb[4]
            area_corners[2] = rec_corners_e_limb[5]
            area_corners[3] = rec_corners_e[3]
            area_corners[4] = rec_corners_e_limb[7]

            area_edges[0][0] = area_corners[0]
            area_edges[1][0] = area_corners[1]
            area_edges[2][0] = area_corners[2]
            area_edges[3][0] = area_corners[3]
            area_edges[4][0] = area_corners[4]

            for i in range(start_index, num_points):
                area_edges[0][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[0] + i / num_points * (
                        area_corners[1] - area_corners[0]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[1] + i / num_points * (area_corners[2] - area_corners[1]) - limb_center
                area_edges[1][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

                area_edges[2][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[2] + i / num_points * (
                        area_corners[3] - area_corners[2]) - pos_spacecraft, radii_a, radii_b, radii_c)

                area_edges[3][i][:] = spice.spiceypy.surfpt(pos_spacecraft, area_corners[3] + i / num_points * (
                        area_corners[4] - area_corners[3]) - pos_spacecraft, radii_a, radii_b, radii_c)

                ray = area_corners[4] + i / num_points * (area_corners[0] - area_corners[4]) - limb_center
                area_edges[4][i][:] = spice.spiceypy.surfpt(limb_center, ray, radii_a, radii_b, radii_c)

            area_edges_in = np.concatenate(([rec_corners[3]], area_edges[4], [area_corners[0]], area_edges[1],
                                            [area_corners[2]], [rec_corners[3]]), axis=0)

        else:
            pass

    # check if any of the poles is contained in the pixel
    north_in, south_in = pixel_contains_poles.main(target, time_et, aber_cor, spacecraft, radii_c,
                                                   pos_spacecraft, limb_plane, xpt1, xpt2, xpt4)

    if 'area_edges_in' not in locals():
        print(f'pixel {i} case not found')
        area_edges_in = []
    else:
        pass

    output_all = np.resize(area_edges, (num_corners * num_points, 3))
    return area_corners, num_corners, output_all, area_edges, rec_corners, area_edges_in, limb_plane, north_in, south_in
