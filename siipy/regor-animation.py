#############################################################
##### Reimplementation of P. Saha and K. N. Rai code ########
#############################################################

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

plt.style.use('dark_background')

from astropy.time import Time
from astropy.constants import R_sun
from scipy.optimize import brentq

pi, cos, sin = np.pi, np.cos, np.sin


def calc_a1(a, q):
    return q / (1 + q) * a


def calc_a2(a, q):
    return -1 / (1 + q) * a


def pos(
        orbital_phase,
        e,
        omega,
        Omega,
        inclination
):
    ec = np.sqrt(1 - e * e)
    psi = brentq(lambda psi: psi - e * sin(psi) - orbital_phase, 0, 2 * pi)

    x, y = cos(psi) - e, ec * sin(psi)
    x, y = x * cos(omega) - y * sin(omega), x * sin(omega) + y * cos(omega)
    x, y = x, y * cos(inclination)
    x, y = x * cos(Omega) - y * sin(Omega), x * sin(Omega) + y * cos(Omega)
    return x, y


def ellipsis(sx, sy, rad, posx=0, posy=0, elong=1, pa=0):
    sx, sy = sx.to_value("m"), sy.to_value("m")
    dx, dy = sx - posx, sy - posy
    dx, dy = cos(pa) * dx + sin(pa) * dy, -sin(pa) * dx + cos(pa) * dy
    return np.where(dx ** 2 * elong + dy ** 2 / elong < rad.value ** 2, 1, 0)


def plot_skymap(x_coord, y_coord, smap, zoom=4,  ceil=None, **ax_kwargs):
    def cen(f):
        N = f.shape[0]
        M = N // (zoom * 2)
        return f[N // 2 - M:N // 2 + M, N // 2 - M:N // 2 + M]

    smap = cen(smap)

    mapmax = smap.max()
    mapmin = smap.min()

    mas = (1 * u.mas).to("", equivalencies=u.dimensionless_angles())
    sx, sy = cen(x_coord.to_value("m")) / mas, cen(y_coord.to_value("m")) / mas
    if ceil is not None:
        mapmin, mapmax = 0, max(ceil, smap.max())
    levels = np.linspace(mapmin, mapmax, 40)
    plt.clf()
    cs = plt.contourf(sx, sy, smap, levels, **ax_kwargs)
    plt.xlabel('mas')
    plt.colorbar(cs)
    plt.tight_layout()
    plt.gca().set_aspect('equal')


def run(
        start_time,
        end_time,
        time_steps,
        period,
        mass_ratio,
        semi_major_axis,
        eccentricity,
        T1,
        T2,
        r1,
        r2,
        omega,
        Omega,
        inclination,
        nbins,
        binsize,
        wavelength
):
    start = Time(start_time, format='jd')
    end = Time(end_time, format='jd')

    time_edges = np.arange(start, end, time_steps)
    phases = [2 * np.pi * ((time.jd / period.value) % 1.0) for time in time_edges]

    a1 = calc_a1(semi_major_axis.value, mass_ratio)
    a2 = calc_a2(semi_major_axis.value, mass_ratio)

    s = np.arange(-nbins // 2, nbins // 2) * binsize
    sky_x, sky_y = np.meshgrid(s, s)
    x = wavelength / (nbins * s)
    ground_x, ground_y = np.meshgrid(x, x)

    for phase in phases[1:]:
        posx, posy = pos(phase, eccentricity, omega, Omega, inclination)
        tmap1 = T1 * ellipsis(sky_x, sky_y, r1, a1 * posx, a1 * posy)
        tmap2 = T2 * ellipsis(sky_x, sky_y, r2, a2 * posx, a2 * posy)
        tmap = tmap1 + tmap2
        plot_skymap(sky_x, sky_y, tmap, zoom=4, cmap='gray')
        plt.pause(2)
    return sky_x, sky_y, ground_x, ground_y, tmap

dis = 380 * u.pc
a = (600 * u.lsec).to("pc") / dis
sky_x, sky_y, ground_x, ground_y, tmap = run(
    246e5,
    246e5 + 8,
    1 * u.h,
    4 * u.day,
    9 / 28,
    a,
    0.326,
    35e3,
    57e3,
    17 * R_sun / dis.to("m"),
    6 * R_sun / dis.to("m"),
    (140 * u.deg).to("rad"),
    (130 * u.deg).to("rad"),
    (117 * u.deg).to("rad"),
    1024,
    0.09 * u.nm,
    465 * u.nm,
)

print(tmap)
