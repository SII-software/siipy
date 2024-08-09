import numpy as np

from matplotlib import pyplot as plt
from matplotlib.colors import Normalize

plt.rcParams["font.family"] = "serif"
plt.rcParams["text.usetex"] = True

plt.rcParams["figure.figsize"] = (8, 6)
plt.rcParams["figure.dpi"] = 300

plt.rcParams["axes.labelsize"] = 16
plt.rcParams["axes.labelweight"] = "bold"

plt.rcParams["xtick.labelsize"] = 14
plt.rcParams["ytick.labelsize"] = 14

from astropy import units as u
from astropy.visualization import quantity_support

from .utils import cen


class Plotter:

    def __init__(
            self,
    ):

        plt.style.use('dark_background')

    def plot_skymap(
            self,
            sky_x,
            sky_y,
            sky_map,
            zoom=4,
            level=40,
            ceiling=None,
            **ax_kwargs,
    ):
        cmap = ax_kwargs.pop("cmap", "gray")
        smap = cen(sky_map, zoom)
        map_max = smap.max()
        map_min = smap.min()

        mas = (1 * u.mas).to("", equivalencies=u.dimensionless_angles())
        sx = cen(sky_x, zoom) / mas
        sy = cen(sky_y, zoom) / mas

        if ceiling is not None:
            map_min, map_max = 0, np.max(ceiling, smap.max())

        levels = np.linspace(map_min, map_max, level)

        return self.plot_map(
            sx.to_value(),
            sy.to_value(),
            smap,
            levels,
            cmap=cmap,
            **ax_kwargs
        )

    def plot_groundmap(
            self,
            ground_x,
            ground_y,
            ground_map,
            zoom=8,
            level=20,
            ceiling=None,
            fceiling=None,
            **ax_kwargs,
    ):
        cmap = ax_kwargs.pop("cmap", "inferno")
        norm = ax_kwargs.get("norm", None)
        gx = cen(ground_x, zoom)
        gy = cen(ground_y, zoom)
        gmap = cen(ground_map, zoom)

        if ceiling is not None:
            map_min, map_max = 0, max(ceiling, ground_map.max())
        elif fceiling is not None:
            map_max = max(fceiling, ground_map.max())
            map_min = -map_max,
        else:
            map_max = ground_map.max()
            map_min = 0

        if norm is None:
            ax_kwargs["norm"] = Normalize(vmin=map_min, vmax=map_max)

        levels = np.linspace(map_min, map_max, level)
        return self.plot_map(
            gx.to_value("m"),
            gy.to_value("m"),
            gmap.to_value(),
            levels,
            cmap=cmap,
            **ax_kwargs
        )

    @staticmethod
    def plot_map(
            coord_x,
            coord_y,
            coord_map,
            levels,
            **ax_kwargs,

    ):
        ax = ax_kwargs.pop("ax", None)
        xlabel = ax_kwargs.pop("xlabel", None)
        ylabel = ax_kwargs.pop("ylabel", None)
        norm = ax_kwargs.pop("norm", None)

        if ax is None:
            fig, ax = plt.subplots()

        fig = plt.gcf()

        with (quantity_support()):
            cs = ax.contourf(coord_x, coord_y, coord_map, levels, norm=norm, **ax_kwargs)
            xmin = np.min(coord_x)
            xmax = np.abs(np.fmin(xmin, np.max(coord_x)))
            ax.set_xlim(xmin, xmax)
            ymin = np.min(coord_y)
            ymax = np.abs(np.fmin(ymin, np.max(coord_y)))
            ax.set_ylim(ymin, ymax)
            if xlabel is not None:
                ax.set_xlabel(xlabel)
            if ylabel is not None:
                ax.set_ylabel(ylabel)
            fig.colorbar(cs)
            fig.tight_layout()
            ax.set_aspect('equal')
