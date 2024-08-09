#############################################################
##### Reimplementation of P. Saha and K. N. Rai code ########
#############################################################
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import matplotlib.pyplot as plt

from astropy import units as u
from astropy import constants
from astropy.time import Time
from astropy.modeling import models

from scipy.optimize import brentq

pi, cos, sin = np.pi, np.cos, np.sin

from numpy.fft import fft2, ifft2, fftshift

from plotting.plotting import Plotter

def fourier(f):
    return fftshift(fft2(fftshift(f)))

def ifourier(f):
    return fftshift(ifft2(fftshift(f)))

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
    #sx, sy = sx.to_value("m"), sy.to_value("m")
    dx, dy = sx - posx, sy - posy
    dx, dy = cos(pa) * dx + sin(pa) * dy, -sin(pa) * dx + cos(pa) * dy
    return np.where(dx ** 2 * elong + dy ** 2 / elong < rad.value ** 2, 1, 0)


def calc_sbright(T,lam,ds):
    z = (constants.h.si * constants.c.si) / (lam.to("m") * constants.k_B.si * (T + 2.725) * u.K)  # to avoid overflow in exp(z)
    z = np.clip(z.to_value(), 0, 20)
    return (ds / lam)**2 / (np.exp(z) - 1)
def calc_correlation_density(sbright,lam):
    N = sbright.shape[0]
    pcoh = np.sum(sbright)
    pflux = 2 * constants.c.si / lam.to("m") * pcoh
    mag = -2.5 * (np.log10(constants.h.si.value * pflux.value) + 22.44)
    print('AB = %5.2f' % mag)
    V = fourier(sbright)
    Vmax = V[N//2, N//2]
    V2 = np.abs(V / Vmax) ** 2
    f = pcoh * V2
    return f

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

    plotter = Plotter()

    fig = plt.figure()
    for phase in phases[1:]:
        ax = fig.subplots(1, 2)
        posx, posy = pos(phase, eccentricity, omega, Omega, inclination)
        tmap1 = T1 * ellipsis(sky_x, sky_y, r1, a1 * posx, a1 * posy)
        tmap2 = T2 * ellipsis(sky_x, sky_y, r2, a2 * posx, a2 * posy)
        tmap = tmap1 + tmap2
        plotter.plot_skymap(sky_x, sky_y, tmap, ax=ax[0],
                                  xlabel='mas', zoom=4, cmap='gray')
        s_bright = calc_sbright(tmap, wavelength, binsize)
        corr_density = calc_correlation_density(s_bright, wavelength)
        flux = 1 * np.sqrt(1e9 * 300) * corr_density
        normalized_flux = flux / np.max(flux)
        plotter.plot_groundmap(ground_x, ground_y, normalized_flux, ceiling=1,
                                     ax=ax[1], zoom=8, cmap='inferno')
        plt.pause(1)
        plt.clf()
        # plt.close()
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
    17 * constants.R_sun / dis.to("m"),
    6 * constants.R_sun / dis.to("m"),
    (140 * u.deg).to("rad"),
    (130 * u.deg).to("rad"),
    (117 * u.deg).to("rad"),
    1024,
    9.e-11,
    465 * u.nm,
)

print(tmap)
