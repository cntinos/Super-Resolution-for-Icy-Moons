"""
Created on February, 2022

@author: Christos Ntinos

This Module calculates if the south or north pole is contained in the pixel.

Version: 1.0
"""

import numpy as np
import spiceypy as spice


def main(target, time_et, aber_cor, spacecraft, radii_c, pos_spacecraft, limb_plane, xpt1, xpt2, xpt4):
    """
    Parameters
    ----------
    target : str
        name of target
    time_et : float
        time of pixel acquisition
    aber_cor : str
        aberration correction flag
    spacecraft : str
        name of spacecraft
    radii_c : float
        radius of the triaxial ellipsoid in the Z direction
    pos_spacecraft : list of floats
        position of spacecraft in rectangular coordinates (km)
    limb_plane : SPICE plane object
        plane containing the limb
    xpt1 : list of floats
        intercepts of corner 1 of the FOV with the limb plane in rectangular coordinates (km)
    xpt2 : list of floats
        intercepts of corner 2 of the FOV with the limb plane in rectangular coordinates (km)
    xpt4 : list of floats
        intercepts of corner 4 of the FOV with the limb plane in rectangular coordinates (km)

    Returns
    -------
    north_in : bool
        North Pole contained in pixel (1) or not (0)
    south_in : bool
        South Pole contained in pixel (1) or not (0)

    This function checks if a pixel contains the North or South Pole and returns a flag 0/1.
    """

    # calculate if any of the poles are contained in the pixel
    trgepc, srfvec, phase, incdnc, emissn, north_visible, lit = spice.spiceypy.illumf(
        "ELLIPSOID", target, 'SUN', time_et, f'IAU_{target}', aber_cor, spacecraft, [0, 0, radii_c])
    trgepc, srfvec, phase, incdnc, emissn, south_visible, lit = spice.spiceypy.illumf(
        "ELLIPSOID", target, 'SUN', time_et, f'IAU_{target}', aber_cor, spacecraft, [0, 0, -radii_c])

    if north_visible is True and south_visible is False:
        south_in = 0
        # calculate projection of the North Pole on the limb plane
        nxpts, xpt = spice.spiceypy.inrypl(pos_spacecraft, [-pos_spacecraft[0], -pos_spacecraft[1], -pos_spacecraft[2]
                                                            + radii_c], limb_plane)
        # check if the projection of the North Pole is inside the pixel
        edge_1 = xpt2-xpt1
        edge_2 = xpt4-xpt1
        pole = xpt-xpt1
        print(xpt1, xpt2, xpt)
        if 0 <= np.inner(edge_1/np.linalg.norm(edge_1), pole)<=np.linalg.norm(edge_1) and 0 <= np.inner(edge_2 / np.linalg.norm(edge_2), pole) <= np.linalg.norm(edge_2):
            north_in = 1  # North Pole contained in the pixel
        else:
            north_in = 0

    elif north_visible is False and south_visible is True:
        north_in = 0
        # calculate projection of the South Pole on the limb plane
        nxpts, xpt = spice.spiceypy.inrypl(pos_spacecraft, [-pos_spacecraft[0], -pos_spacecraft[1], -pos_spacecraft[2]
                                                            - radii_c], limb_plane)
        # check if the projection of the North Pole is inside the pixel
        edge_1 = xpt2-xpt1
        edge_2 = xpt4-xpt1
        pole = xpt-xpt1
        if 0<=np.inner(edge_1/np.linalg.norm(edge_1), pole)<=np.linalg.norm(edge_1) and 0 <= np.inner(edge_2 / np.linalg.norm(edge_2), pole) <= np.linalg.norm(edge_2):
            south_in = 1  # South Pole contained in the pixel
        else:
            south_in = 0
    else:
        north_in = 0
        south_in = 0

    return north_in, south_in

