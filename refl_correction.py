"""
Created on February, 2022

@author: Christos Ntinos

This Module takes as input 1) the corners of the pixel, 2) the number of corners of the pixel, 3) the limb plane and
4) the vector containing all the points tracing the edge of the part of the pixel that corresponds to the target
and calculates the ration of the area inside the target of the limb-crossing pixels over the total area of the pixel.

Version: 1.0
"""
import numpy as np
import spiceypy as spice

def main(area_corners, num_corners, limb_plane, area_edges_in, rec_corners):

    vector_in = np.zeros(3)
    # calculate the vector passing through the area of the pixel corresponding to the target
    for i in range(0, len(area_edges_in)):
        vector_in = vector_in + area_edges_in[i]

    # calculate the intersection of the vector_in with the limb_plane
    nxpts, xpt = spice.spiceypy.inrypl(np.array([0, 0, 0]), vector_in, limb_plane)

    # calculate the area of the pixel corresponding to the targer
    area_in = 0
    for i in range(0, area_edges_in.shape[0]-1):
        area_in = area_in + 0.5 * np.linalg.norm(np.cross(area_edges_in[i]-xpt, area_edges_in[i+1]-xpt))

    # calculate the total area of the pixel on the limb plane
    area_all = np.linalg.norm(rec_corners[0]-rec_corners[1]) * np.linalg.norm(rec_corners[0] - rec_corners[3])

    # calculate the correction coefficient
    refl_correction_coef = area_in/area_all

    return refl_correction_coef