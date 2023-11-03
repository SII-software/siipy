import numpy as np
import matplotlib.pyplot as pl
from visib import sbright, correldens
from graphics import draw

from fakeprox import teff
from layout import Layout

pl.style.use('dark_background')

source = Layout()

source.set_radec((14,29,42.94),(-62,-40,-46.13)) # Proxima Cen
source.set_latlon((-23,-16,-25.4856),(16,31,10.2432))  # HESS

def stretch(f):
    return np.sign(f) * np.abs(f)**(1/4)

def illus(fname):
    pl.pause(10)
    #print(fname)
    #pl.savefig(fname)

def run():
    F = 100
    N = 2048
    ds = 2e-11
    lam = 800e-9
    gz = 16
    sx,sy,x,y,Tmap,vari = teff(N,ds,lam)
    S = sbright(Tmap,lam,ds)
    f = correldens(S,lam)
    #fv = correldens(S*vari,lam)
    f /= np.max(f)
    #fv /= np.max(fv)
    masf = (np.pi/180/3.6e6/ds)**2
    smax = np.max(S*vari*masf)
    #draw(sx,sy,S*masf,8,'sky',ceil=smax,cmap='magma',title=til)
    #illus('figs/src_ud')
    til = ('photons / (m^2 s Hz mas^2) at %i nm' % round(1e9*lam))
    draw(sx,sy,S*masf,8,'sky',ceil=smax,cmap='magma',title=til)
    illus('figs/src_ud')
    draw(sx,sy,S*masf,8,'sky',ceil=smax,cmap='magma',title=til,lam=lam)
    illus('figs/src_udn')
    #til = ('|visib|^2 at %i nm' % (1e9*lam+0.5))
    #draw(x,y,f,gz,'ground',fceil=1,cmap='coolwarm',title=til)
    #illus('figs/visibq')
    #sig = stretch(f)
    #til = ('|visib| at %i nm' % (1e9*lam+0.5))
    #draw(x,y,sig,gz,'ground',fceil=1,cmap='coolwarm',title=til)
    #illus('figs/visib')
    sig = stretch(f)
    til = ('(corr)^(1/4) at %i nm' % (1e9*lam+0.5))
    
    for fr in range(48):
        jd = 246e5 + (fr-7)/48
        u,v,w = source.get_uvw(jd,dx,dy,0*dx)
        hx,hy,hz = source.get_xyz(jd,0,0,1)
        print('time %.3f altitude %.3f' % (jd-246e5,hz))
        draw(x,y,sig,gz,'ground',fceil=np.max(abs(sig)),cmap='coolwarm',title=til)
        if hz > 0:
            pl.plot(u,v,'+',color='white')
        print(fr)
        illus('tmp/corr%i'%fr)
        

x = []
y = []
fil = open('locs-hess.txt')
while True:
    s = fil.readline()
    if not s:
        break
    xs,ys,c = s.split()
    x.append(float(xs))
    y.append(float(ys))

dx = []
dy = []
for i in range(len(x)):
    for j in range(len(x)):
        if i != j:
            dx.append(x[i]-x[j])
            dy.append(y[i]-y[j])
dx = np.array(dx)
dy = np.array(dy)

run()

