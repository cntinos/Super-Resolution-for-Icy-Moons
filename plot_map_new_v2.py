'''
Created on July, 2022

@author: Christos Ntinos

This Module plots maps of selected variable based on a dataframe input.
IAU standard reference frame is planetographic west positive (0-->360W).

Version: 2.0 (March, 2026)
vmin, vmax and method arguments (for West positive coordinates according to IAU) added
Fix bug caused by sort in pixels containing North/South Pole.

'''

import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib as mpl
from constants import  font_size


def plot_map(df, polygon_coordinates, ax, plot_variable, cmap_name, cbar_title, vmin=0, vmax=1, method='Zero_center',
             cbar_mode='on'):
    """
    Parameters
    ----------
    df : pandas.DataFrame
        input dataframe map
    polygon_coordinates : list of lists
        polygon_coordinates = [longitudes, latitudes, nav_cube_all_25, nav_cube_all_26]
    ax : Axes object
        axis on which to plot the map
    plot_variable : str
        variable to be plotted as the color scale (band, photometric angles etc.)
    cmap_name : str
        color map name to be used
    cbar_title : str
        color bar title
    vmin : float
        minimum value in the color scale, default value vmin=0
    vmax : float
        maximum value in the color scale, default value vmin=1
    method: str
        direction of the map: 'Zero_center' for positive east (180W-->180E)(default) ot 'West' for positive west
        (0-->360W)
    cbar_mode : str
        color bar added ('on') or not ('off'), default value cbar_mode='on'

    Returns
    -------
    out : Axis object
        map to be plotted

    This function plots a map from the input dataframe.
    """

    longitude_all, latitude_all, nav_cube_all_25, nav_cube_all_26 = polygon_coordinates  # all pixels (sky, limb, central)

    df_sorted = df.sort_values(by='sampling', ascending=False)
    index_sorted = df_sorted.index.values.tolist()

    # PLOT using matplotlib patches (static)

    polygons_all = []
    polygons_all_value = []

    if method == 'Zero_center':  # planetocentric centered at 0-0

        for i in index_sorted:
            if len(longitude_all[i]) != 0 and ~np.isnan(df[plot_variable][i]):  # there are empty pixels [] corresponding to the sky

                pixel_color = df[plot_variable][i]

                if nav_cube_all_25[i] == 1:  # north pole in the pixel
                    index_max = longitude_all[i].index(max(longitude_all[i]))
                    index_min = longitude_all[i].index(min(longitude_all[i]))
                    if index_min < index_max:
                        x_pixel_pole = longitude_all[i][:index_min + 1] + [-180, -180, 180, 180] + longitude_all[i][
                                                                                                   index_min + 1:]
                        y_pixel_pole = latitude_all[i][:index_min + 1] + [latitude_all[i][index_min], 90, 90,
                                                                          latitude_all[i][index_max]] + latitude_all[i][
                                                                                                        index_min + 1:]
                    elif index_min > index_max:
                        x_pixel_pole = longitude_all[i][:index_max + 1] + [180, 180, -180, -180] + longitude_all[i][
                                                                                                   index_max + 1:]
                        y_pixel_pole = latitude_all[i][:index_max + 1] + [latitude_all[i][index_max], 90, 90,
                                                                          latitude_all[i][index_min]] + latitude_all[i][
                                                                                                        index_max + 1:]
                    else:
                        raise Exception("index_min=index_max in longitude!")

                    polygon = Polygon(list(zip(x_pixel_pole, y_pixel_pole)))
                    polygons_all.append(polygon)
                    polygons_all_value.append(pixel_color)

                elif nav_cube_all_26[i] == 1:  # south pole in the pixel
                    index_max = longitude_all[i].index(max(longitude_all[i]))
                    index_min = longitude_all[i].index(min(longitude_all[i]))
                    if index_min < index_max:
                        x_pixel_pole = longitude_all[i][:index_min + 1] + [-180, -180, 180, 180] + longitude_all[i][
                                                                                                   index_min + 1:]
                        y_pixel_pole = latitude_all[i][:index_min + 1] + [latitude_all[i][index_min], -90, -90,
                                                                          latitude_all[i][index_max]] + latitude_all[i][
                                                                                                        index_min + 1:]
                    elif index_min > index_max:
                        x_pixel_pole = longitude_all[i][:index_max + 1] + [180, 180, -180, -180] + longitude_all[i][
                                                                                                   index_max + 1:]
                        y_pixel_pole = latitude_all[i][:index_max + 1] + [latitude_all[i][index_max], -90, -90,
                                                                          latitude_all[i][index_min]] + latitude_all[i][
                                                                                                        index_max + 1:]
                    else:
                        raise Exception("index_min=index_max in longitude!")

                    polygon = Polygon(list(zip(x_pixel_pole, y_pixel_pole)))
                    polygons_all.append(polygon)
                    polygons_all_value.append(pixel_color)

                else:
                    if len(np.unique(np.sign(longitude_all[i]))) == 1:
                        polygon = Polygon(list(zip(longitude_all[i], latitude_all[i])))
                        polygons_all.append(polygon)
                        polygons_all_value.append(pixel_color)
                    else:
                        list_shift = [*longitude_all[i][1:], longitude_all[i][0]]
                        if any(abs(lon1-lon2) > 180 for lon1, lon2 in zip(longitude_all[i], list_shift)):

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

                            polygon_1 = Polygon(list(zip(longitude_all_1, latitude_all_1)))
                            polygon_2 = Polygon(list(zip(longitude_all_2, latitude_all_2)))
                            polygons_all.extend([polygon_1, polygon_2])
                            polygons_all_value.extend([pixel_color, pixel_color])

                        else:
                            polygon = Polygon(list(zip(longitude_all[i], latitude_all[i])))
                            polygons_all.append(polygon)
                            polygons_all_value.append(pixel_color)

            else:
                pass

    elif method == 'West':  # planetocentric positive west (0-->360W)

        for i in index_sorted:
            if len(longitude_all[i]) != 0 and ~np.isnan(df[plot_variable][i]):  # there are empty pixels [] corresponding to the sky

                pixel_color = df[plot_variable][i]

                lon_west = np.where(np.array(longitude_all[i]) >= 0, 360-np.array(longitude_all[i]), -np.array(longitude_all[i]))

                if nav_cube_all_25[i] == 1:  # north pole in the pixel
                    x_pixel_pole = [x for x, _ in sorted(zip(lon_west, latitude_all[i]))]
                    y_pixel_pole = [y for _, y in sorted(zip(lon_west, latitude_all[i]))]

                    x_pixel_pole.extend([360, 360, 0, 0])
                    y_pixel_pole.extend([y_pixel_pole[-1], 90, 90, y_pixel_pole[0]])

                    polygon = Polygon(list(zip(x_pixel_pole, y_pixel_pole)))
                    polygons_all.append(polygon)
                    polygons_all_value.append(pixel_color)

                elif nav_cube_all_26[i] == 1:  # south pole in the pixel

                    x_pixel_pole = [x for x, _ in sorted(zip(lon_west, latitude_all[i]))]
                    y_pixel_pole = [y for _, y in sorted(zip(lon_west, latitude_all[i]))]

                    x_pixel_pole.extend([360, 360, -180, -180])
                    y_pixel_pole.extend([y_pixel_pole[-1], -90, -90, y_pixel_pole[0]])

                    polygon = Polygon(list(zip(x_pixel_pole, y_pixel_pole)))
                    polygons_all.append(polygon)
                    polygons_all_value.append(pixel_color)

                else:
                    if max(lon_west)-min(lon_west) > 180:

                        polygon_1 = Polygon(list(zip(-np.array(longitude_all[i]), latitude_all[i])))
                        polygon_2 = Polygon(list(zip(360-np.array(longitude_all[i]), latitude_all[i])))
                        polygons_all.extend([polygon_1, polygon_2])
                        polygons_all_value.extend([pixel_color, pixel_color])

                    else:
                        polygon = Polygon(list(zip(lon_west, latitude_all[i])))
                        polygons_all.append(polygon)
                        polygons_all_value.append(pixel_color)

            else:
                pass

    p = PatchCollection(polygons_all, cmap=cmap_name, alpha=1, norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax))
    p.set_array(polygons_all_value)
    out = ax.add_collection(p)
    if cbar_mode == 'on':
        cbar = plt.colorbar(p, orientation='vertical', shrink=0.6,
                            ticks=np.linspace(vmin, vmax, 5, endpoint=True), format="%.2f")
        cbar.ax.tick_params(labelsize=font_size)
        cbar.set_label(cbar_title, fontsize=font_size)
    elif cbar_mode == 'off':
        pass

    return out


if __name__ == "__main__":
    # run function
    # load data (deserialize) to map
    with open('_data/sr_global_map_ds4_wvln_1_80401_cor_0_0_0_5_DSK_unsaturated_with_bns_97_1440x720_2_2D.bin', 'rb') as handle:
        my_data = pickle.load(handle)

    # unpack my_data
    df, polygon_coordinates = my_data
    longitude_all, latitude_all, nav_cube_all_25, nav_cube_all_26 = polygon_coordinates

    # select band to plot
    band_id = 'I/F_1.80401_sr'

    range_min = df[band_id].min()
    range_max = df[band_id].max()

    fig, axs = plt.subplots(1, 1)
    plot_map(df, polygon_coordinates, axs, band_id, 'gray', 'SR I/F (1.80401 μm)',
             vmin=np.round(range_min), vmax=np.round(range_max), method='Zero_center')
    axs.set_aspect('equal')
    axs.set_xlim(-180, 180)
    axs.set_ylim(-90, 90)
    axs.set_xticks(range(-180, 181, 60),
               labels=[r'$180^{\circ}$W', r'$120^{\circ}$W', r'$60^{\circ}$W', r'$0^{\circ}$', r'60$^{\circ}$E',
                       r'120$^{\circ}$E', r'180$^{\circ}$E'], fontsize=font_size, ha='left')
    axs.set_yticks(range(-90, 91, 30),
               labels=[r'90$^{\circ}$S', r'60$^{\circ}$S', r'30$^{\circ}$S', r'0$^{\circ}$', r'30$^{\circ}$N',
                       r'60$^{\circ}$N', r'90$^{\circ}$N'], fontsize=font_size, va='center')

    fig.suptitle(
        f'Enceladus SR global map 1440X720 grid at {band_id.replace("I/F_", "").replace("_sr", "")} μm corrected, 5 points/edge, DSK, data-set 4: 1566 cubes' + '\n' +
        r'target pixels, i, e $\leq$ 80$^{\circ}$, 0.001 $\leq$ I/F $\leq$ 1, sampling $\leq$ 78 km/pixel (155720 pixels)' +
        '\n saturation limit = (2**12-bns)*0.97, general correction, wi = (Ai/Atotal)^4, 3D area')

    plt.show()

