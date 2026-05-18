"""
Created on September, 2022

@author: Christos Ntinos

This Module saves an equirectangular projection of the input map, by resolving the pixels containing the N or S pole and
the pixels that cross the antimeridian.

Version: 2.0 (April, 2026)
Fix bug caused by sort in pixels containing North/South Pole.

"""

import pickle
import pandas as pd
from shapely.geometry.polygon import LinearRing, LineString


def save_map_projected(file_name):
    """
    Calculates photometrically corrected map.

    Parameters
    ----------
    file_name : str
        name of output file containing the photometrically corrected map projected(.bin)

    Returns
    -------
    df_new : pandas.DataFrame
        photometrically corrected map projected
    polygon_coordinates: list of lists
        polygon_coordinates = [longitudes, latitudes, nav_cube_all_25, nav_cube_all_26]
    """

    # load data (deserialize)
    with open('_data/' + file_name + '.bin', 'rb') as handle:
        df, polygon_coordinates = pickle.load(handle)

    longitude_all, latitude_all, nav_cube_all_25, nav_cube_all_26 = polygon_coordinates  # all pixels (sky, limb, central)

    df_index_max = len(longitude_all)  # maximum index of the unfiltered map

    counter = 0
    rows_to_add = []

    for i in df.index.values.tolist():
        if len(longitude_all[i]) != 0:  # there are empty pixels [] corresponding to the sky

            pixel_points = list(zip(longitude_all[i], latitude_all[i]))
            ring = LinearRing(pixel_points)
            lin_string = LineString(pixel_points)

            if nav_cube_all_25[i] == 1:  # North Pole in the pixel
                index_max = longitude_all[i].index(max(longitude_all[i]))
                index_min = longitude_all[i].index(min(longitude_all[i]))
                if index_min < index_max:
                    x_pixel_pole = longitude_all[i][:index_min + 1] + [-180, -180, 180, 180] + longitude_all[i][index_min + 1:]
                    y_pixel_pole = latitude_all[i][:index_min + 1] + [latitude_all[i][index_min], 90, 90, latitude_all[i][index_max]] + latitude_all[i][index_min + 1:]
                elif index_min > index_max:
                    x_pixel_pole = longitude_all[i][:index_max + 1] + [180, 180, -180, -180] + longitude_all[i][index_max + 1:]
                    y_pixel_pole = latitude_all[i][:index_max + 1] + [latitude_all[i][index_max], 90, 90, latitude_all[i][index_min]] + latitude_all[i][index_max + 1:]
                else:
                    raise Exception("index_min=index_max in longitude!")

                longitude_all[i] = x_pixel_pole
                latitude_all[i] = y_pixel_pole

            elif nav_cube_all_26[i] == 1:  # South Pole in the pixel
                index_max = longitude_all[i].index(max(longitude_all[i]))
                index_min = longitude_all[i].index(min(longitude_all[i]))
                if index_min < index_max:
                    x_pixel_pole = longitude_all[i][:index_min + 1] + [-180, -180, 180, 180] + longitude_all[i][index_min + 1:]
                    y_pixel_pole = latitude_all[i][:index_min + 1] + [latitude_all[i][index_min], -90, -90, latitude_all[i][index_max]] + latitude_all[i][index_min + 1:]
                elif index_min > index_max:
                    x_pixel_pole = longitude_all[i][:index_max + 1] + [180, 180, -180, -180] + longitude_all[i][index_max + 1:]
                    y_pixel_pole = latitude_all[i][:index_max + 1] + [latitude_all[i][index_max], -90, -90, latitude_all[i][index_min]] + latitude_all[i][index_max + 1:]
                else:
                    raise Exception("index_min=index_max in longitude!")

                longitude_all[i] = x_pixel_pole
                latitude_all[i] = y_pixel_pole

            else:
                # extend NP and SP contained
                nav_cube_all_25.append(nav_cube_all_25[i])
                nav_cube_all_26.append(nav_cube_all_26[i])

                list_shift = [*longitude_all[i][1:], longitude_all[i][0]]
                if any(abs(lon1-lon2) > 180 for lon1, lon2 in zip(longitude_all[i], list_shift)):
                    print(i)

                    longitude_all_1 = []
                    latitude_all_1 = []
                    longitude_all_2 = []
                    latitude_all_2 = []

                    for j in range(0, len(longitude_all[i])):

                        latitude_all_1.append(latitude_all[i][j])
                        latitude_all_2.append(latitude_all[i][j])

                        if longitude_all[i][j] >= 0:
                            longitude_all_1.append(longitude_all[i][j])
                            longitude_all_2.append(longitude_all[i][j] - 360)
                        else:
                            longitude_all_1.append(longitude_all[i][j] + 360)
                            longitude_all_2.append(longitude_all[i][j])

                    longitude_all[i] = longitude_all_1
                    longitude_all.append(longitude_all_2)
                    latitude_all[i] = latitude_all_1
                    latitude_all.append(latitude_all_2)

                    # add line at the end of dataframe for duplicate pixel
                    rows_to_add.append(df.loc[i])

        else:
            pass
        counter = counter + 1

    # adjust index
    index_to_add = range(df_index_max, df_index_max+len(rows_to_add))
    df_to_add = pd.DataFrame(rows_to_add).set_axis(index_to_add)
    df_new = pd.concat([df, df_to_add])

    # pack data
    polygon_coordinates = longitude_all, latitude_all, nav_cube_all_25, nav_cube_all_26
    # store nmap (serialize)
    with open('_data/' + file_name + '_projected.bin', 'wb') as handle:
        pickle.dump([df_new, polygon_coordinates], handle, protocol=pickle.HIGHEST_PROTOCOL)

    return df_new, polygon_coordinates

