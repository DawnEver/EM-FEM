from enum import Enum

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import LinearNDInterpolator

__all__ = ['plot_map', 'PlotMapType']


class PlotMapType(Enum):
    Scatter = 'scatter'
    Coutour = 'coutour'
    Coutourf = 'coutourf'
    Quiver = 'quiver'


def plot_map(
    title: str, vertInfoMat: np.ndarray, c_mat: np.ndarray, boundary, plot_type: PlotMapType = PlotMapType.Scatter
):
    plt.title(title)
    c_mat[c_mat < boundary[0]] = boundary[0]
    c_mat[c_mat > boundary[1]] = boundary[1]
    x_raw, y_raw = vertInfoMat[:, 0], vertInfoMat[:, 1]

    match plot_type:
        case PlotMapType.Scatter:
            cax = plt.scatter(x_raw, y_raw, c=c_mat)
        case PlotMapType.Quiver:
            cax = plt.quiver(x_raw, y_raw, c_mat[:, 0], c_mat[:, 1])
        case _:
            interp_func = LinearNDInterpolator(vertInfoMat, c_mat)
            levels = 20
            n_x_div = 100
            n_y_div = 100
            x_new = np.linspace(x_raw.min(), x_raw.max(), n_x_div)
            y_new = np.linspace(y_raw.min(), y_raw.max(), n_y_div)
            X, Y = np.meshgrid(x_new, y_new)
            sample_points = np.vstack([X.flatten(), Y.flatten()]).T
            interp_value = interp_func(sample_points).reshape([n_x_div, n_y_div])
            if plot_type == PlotMapType.Coutour:
                cax = plt.contour(X, Y, interp_value, levels=levels)
            elif plot_type == PlotMapType.Coutourf:
                cax = plt.contourf(X, Y, interp_value, levels=levels)
    if 0:
        from matplotlib import ticker

        plt.colorbar(cax, ticks=ticker.LogLocator(), format=ticker.LogFormatterSciNotation())
    else:
        plt.colorbar()
    plt.show()
