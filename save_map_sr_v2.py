"""
Created on December, 2025

@author: Christos Ntinos (IPGP)

This Module calculates the pixels from the dataset contained in each super-pixel of a regular rectangular grid.

Version: 2.0 (January 2026) (+) STRtree

"""

import pickle
import matplotlib
import numpy as np
import pandas as pd
from shapely.geometry import Polygon  # two different "Polygon" objects
from matplotlib.patches import Polygon as Polygon_patch
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
from filter_map import filter_map_fun
import os.path
from constants import font_size
import spiceypy as spice
import auto_mk_selection
from tqdm import tqdm
from shapely.strtree import STRtree
from pyproj import Geod


def create_super_cells(minX, maxX, num_points_x, minY, maxY, num_points_y):
    # create one-dimensional arrays for x and y
    x = np.linspace(minX, maxX, num_points_x)
    y = np.linspace(minY, maxY, num_points_y)

    # create super-cells
    longitude_all_super = []
    latitude_all_super = []
    for i in range(0, num_points_x-1):
        for j in range(0, num_points_y-1):
            longitude_all_super.append([x[i], x[i+1], x[i+1], x[i]])
            latitude_all_super.append([y[j+1], y[j+1], y[j], y[j]])
    return longitude_all_super, latitude_all_super


def load_map(file_name_input):
    # load data (deserialize)
    with open(file_name_input, 'rb') as handle:
        df, polygon_coordinates = pickle.load(handle)

    longitude_all, latitude_all, _, _ = polygon_coordinates  # all pixels (sky, limb, central)

    return df, longitude_all, latitude_all


def calculate_overlapping_pixels(df, longitude_all, latitude_all, longitude_all_super, latitude_all_super, file_name):

    # check if file with the overlapping pixels already exists, if True skip calculation (computationally expensive)
    if not os.path.exists('_data/polygon_contained_' + file_name + '.bin'):

        polygon_contains_all = []  # list of lists of indexes of pixels containing each super resolution pixel

        polygon_all = []
        polygon_indices = []
        for i in df.index.values.tolist():
            polygon = Polygon(list(zip(longitude_all[i], latitude_all[i])))
            polygon_all.append(polygon)
            polygon_indices.append(i)
        indices_all = np.array(polygon_indices)

        tree = STRtree(polygon_all)  # build spatial index using the Sort-Tile-Recursive (STR) algorithm

        for k in tqdm(range(0, len(longitude_all_super)), desc="Polygons contained"):
            polygon_super = Polygon(list(zip(longitude_all_super[k], latitude_all_super[k])))
            tree_indices = tree.query(polygon_super, predicate="intersects")  # fast bounding-box pruning
            polygon_contains_all.append(indices_all[tree_indices])

        # store data (serialize) for polygons contained in each SR pixel
        with open('_data/polygon_contained_' + file_name + '.bin', 'wb') as handle:
            pickle.dump([[longitude_all_super, latitude_all_super], polygon_contains_all], handle,
                        protocol=pickle.HIGHEST_PROTOCOL)

    else:
        # read data (deserialize) for polygons contained in each SR pixel
        print('Polygons contained pre-calculated.')
        with open('_data/polygon_contained_' + file_name + '.bin', 'rb') as handle:
            _, polygon_contains_all = pickle.load(handle)

    return polygon_contains_all


def calculate_sr_map(df, band_id, weight_power, file_name, longitude_all, latitude_all, polygon_contains_all,
                     longitude_all_super, latitude_all_super, weight_mode='2D'):


    # Use IAU mean radius = 252.1 km
    geod = Geod(a=256200, b=248600)  # metres

    # short pixels by sampling
    df_sorted = df.sort_values(by='sampling', ascending=False)
    index_sorted = df_sorted.index.values.tolist()

    # check if file with the super-resolution map already exists, if True skip calculation (computationally expensive)
    if not os.path.exists('_data/' + file_name + f'_{weight_power}_{weight_mode}_new.bin'):
        # calculate SR pixel value
        colormap_name = 'gray'
        cmap = matplotlib.colormaps[colormap_name]
        norm = Normalize(vmin=0.5, vmax=0.85)  # Normalize(vmin=np.min(polygons_all_value), vmax=np.max(polygons_all_value))
        plot_variable = band_id + '_cor'

        # calculate new pixel value
        super_value_all = []
        polygon_super_all = []
        super_value_all_avg = []
        super_value_num = []
        super_value_all_sigma = []
        sampling_all = []

        for i in tqdm(range(0, len(polygon_contains_all)), desc="SR map"):
            polygon_super = Polygon(list(zip(longitude_all_super[i], latitude_all_super[i])))
            polygon_patch = Polygon_patch(list(zip(longitude_all_super[i], latitude_all_super[i])))
            polygon_super_all.append(polygon_patch)

            index_contained = polygon_contains_all[i]
            super_value = df_sorted[df_sorted.index.isin(index_contained)].tail(1)[plot_variable]
            super_value_all.append(super_value)
            sampling_all.append(df_sorted[df_sorted.index.isin(index_contained)].tail(1)['sampling'].values)

            if weight_mode == '2D':
                index = polygon_contains_all[i]
                pixel_points_all = [list(zip(longitude_all[j], latitude_all[j])) for j in index]
                polygon_all = [Polygon(pixel_points) for pixel_points in pixel_points_all]

                weight_all_list = []
                for index in range(0, len(polygon_contains_all[i])):
                    print(df.loc[polygon_contains_all[i][index]])
                    polygon = polygon_all[index]
                    weight_all_list.append((polygon.intersection(polygon_super).area / polygon.area) ** weight_power)
                weight_all = np.array(weight_all_list)
                #weight_all = np.array([(polygon.intersection(polygon_super).area / polygon.area) ** weight_power for polygon in polygon_all])
                weight_sum = sum(weight_all)
                super_value_wavg = sum(weight_all * df.loc[index][band_id+'_cor'])
                super_value_sigma = sum((weight_all ** 2) * (df.loc[index][band_id + '_cor_sigma'] ** 2))

            elif weight_mode == '3D':
                weight_sum = 0
                super_value_wavg = 0
                super_value_sigma = 0
                for j in range(0, len(polygon_contains_all[i])):
                    index = polygon_contains_all[i][j]
                    pixel_points = list(zip(longitude_all[index], latitude_all[index]))
                    polygon = Polygon(pixel_points)
                    polygon_inter = polygon.intersection(polygon_super)

                    if polygon_inter.geom_type == 'MultiPolygon':
                        weight_multy = []
                        for n in range(0, len(polygon_inter.geoms)):
                            longitude_inter, latitude_inter = polygon_inter.geoms[n].exterior.coords.xy
                            area_inter, _ = geod.polygon_area_perimeter(longitude_inter, latitude_inter)
                            weight_multy.append(np.abs(area_inter))
                        area_inter = sum(weight_multy)
                    else:
                        longitude_inter, latitude_inter = polygon_inter.exterior.coords.xy
                        area_inter, _ = geod.polygon_area_perimeter(longitude_inter, latitude_inter)

                    area_pixel, _ = geod.polygon_area_perimeter(longitude_all[index], latitude_all[index])

                    weight = (np.abs(area_inter) / np.abs(area_pixel)) ** weight_power  # test different powers

                    weight_sum += weight
                    super_value_wavg += weight * df.loc[index][band_id + '_cor']
                    super_value_sigma += (weight * df.loc[index][band_id + '_cor_sigma'])**2

            try:
                super_value_num.append(super_value_wavg / weight_sum)
                super_value_all_sigma.append(np.sqrt(super_value_sigma/(weight_sum**2)))

            except:  # no pixels contained --> background pixel
                super_value_all_avg.append([1, 1, 1, 1])  # white background pixels
                super_value_num.append(np.nan)
                super_value_all_sigma.append(np.nan)

        d = {'lon_0': (np.array(longitude_all_super).T[0]+np.array(longitude_all_super).T[1])/2,
             'lon_1': np.array(longitude_all_super).T[0],
             'lon_2': np.array(longitude_all_super).T[1],
             'lon_3': np.array(longitude_all_super).T[2],
             'lon_4': np.array(longitude_all_super).T[3],
             'lat_0': (np.array(latitude_all_super).T[1]+np.array(latitude_all_super).T[2])/2,
             'lat_1': np.array(latitude_all_super).T[0],
             'lat_2': np.array(latitude_all_super).T[1],
             'lat_3': np.array(latitude_all_super).T[2],
             'lat_4': np.array(latitude_all_super).T[3],
             'sampling': sampling_all,
             f'{band_id}_sr': super_value_num,
             f'{band_id}_sr_sigma': super_value_all_sigma}

        df_sr = pd.DataFrame(data=d)

        nav_cube_all_25 = len(longitude_all_super) * [0]
        nav_cube_all_26 = len(longitude_all_super) * [0]
        polygon_coordinates_sr = [longitude_all_super, latitude_all_super, nav_cube_all_25, nav_cube_all_26]

        # store data (serialize) for polygons contained in each SR map
        with open('_data/' + file_name + f'_{weight_power}_{weight_mode}.bin', 'wb') as handle:
            pickle.dump([df_sr, polygon_coordinates_sr], handle, protocol=pickle.HIGHEST_PROTOCOL)

    else:
        # read data (deserialize) for polygons contained in each SR pixel
        print('Super-resolution map pre-calculated.')
        with open('_data/' + file_name + f'_{weight_power}.bin', 'rb') as handle:
            df_sr, polygon_coordinates_sr = pickle.load(handle)

        return df_sr


if __name__ == "__main__":
    # define SR weight power
    weight_power = 2

    # create new polygon pixels
    num_points_x = 1441  # num_cells+1
    num_points_y = 721

    # define the lower and upper limits for x and y
    minX, maxX, minY, maxY = -180, 180, -90, 90  # global

    longitude_all_super, latitude_all_super = create_super_cells(minX, maxX, num_points_x, minY, maxY, num_points_y)

    # file name definition
    band_id = 'I/F_1.80401'  # 'I/F_2.00167'  # 'I/F_1.65567'  #'I/F_3.09689'  # 'I/F_1.80401'
    file_name = f'sr_global_map_ds4_wvln_{band_id.replace("I/F_","").replace(".","_")}_cor_0_0_0_5_DSK_unsaturated_with_bns_97_{num_points_x - 1}x{num_points_y - 1}'
    file_name_input = f'_data/global_map_wvln_{band_id.replace("I/F_","").replace(".","_")}_cor_0_0_0_5_DSK_ds4_unsaturated_with_bns_97_projected.bin'

    df, longitude_all, latitude_all = load_map(file_name_input)
    #df = df[df['cube'] != '1487302144_1']  # remover HR reference

    polygon_contains_all = calculate_overlapping_pixels(df, longitude_all, latitude_all, longitude_all_super,
                                                        latitude_all_super, file_name)

    df_sr = calculate_sr_map(df, band_id, weight_power, file_name, longitude_all, latitude_all, polygon_contains_all,
                             longitude_all_super, latitude_all_super)

