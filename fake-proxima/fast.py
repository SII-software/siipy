import numpy as np
from hbtsimul import ampl, exposure, stats
from fakeprox import teff
from graphics import draw
from matplotlib.pyplot import pause, savefig, style

style.use('dark_background')

def illus(fname):
    #pause(0.1)
    print(fname)
    savefig(fname)

N = 2048
ds = 2e-11
lam = 800e-9

gz = 16

sx,sy,x,y,Tmap,vari = teff(N,ds,lam)


S = ampl(Tmap,lam,ds)

dt = 10
sdt = 0.01
F = 100

fr = 0
for t in np.linspace(0,dt,F):
    phi2 = exposure(S,t,sdt,lam)/sdt
    phi2 *= 1e6
    til = 'photons km$^{-2}$ s$^{-1}$ Hz$^{-1}$'
    draw(x,y,phi2,gz,'ground',ceil=0.48,title=til)
    fname = 'tmp/decoh'+('%i' % fr)
    illus(fname)
    fr += 1


X = 0*x
fr = 0
for t in np.linspace(10,10+dt,F):
    phi2 = exposure(S,t,sdt,lam)/sdt
    phi2 *= 1e6
    X += phi2/F * dt
    til = 'photons km$^{-2}$ s$^{-1}$ Hz$^{-1}$'
    draw(x,y,phi2,gz,'ground',ceil=.48,title=til)
    fname = 'tmp/devel'+('%i' % fr)
    illus(fname)
    fr += 1
til = 'photons km$^{-2}$'
draw(x,y,X,gz,'ground',title=til)
illus('figs/intens')

g = stats(X)
draw(x,y,g,gz,'ground',cmap='coolwarm',fceil=.001,title='correlation')
illus('figs/correl')


