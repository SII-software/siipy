import numpy as np
import matplotlib.pyplot as pl
from visib import sbright, correldens
from graphics import draw

from fakeprox import teff
from layout import Layout

pl.style.use('dark_background')

source = Layout()

source.set_radec((14,29,42.94),(-62,40,46.13)) # Proxima Cen
source.set_latlon((-23,16,25.4856),(16,31,10.2432))  # HESS

def stretch(f):
    return np.sign(f) * np.abs(f)**(1/4)

def run():
    F = 96
    N = 2048
    ds = 2.58e-11
    lam = 800e-9
    
    sx,sy,x,y,Tmap,vari = teff(N,ds,lam)
    S = sbright(Tmap,lam,ds)
    #pl.show()
    #f *= 1 * (1e9 * 300)**.5  # A = 1 m^2 and t_obs = 5 min
    f = correldens(S,lam)
    fv = correldens(S*vari,lam)
    f /= np.max(f)
    fv /= np.max(fv)
    til = ('photons / (m^2 s Hz sr)')
    draw(sx,sy,S*vari/ds**2,8,'sky',cmap='magma',title=til)
    #pl.savefig('tmp/src')
    pl.pause(2)
    sig = stretch(f)
    til = ('correlation^(1/4) at %i nm' % (1e9*lam+0.5))
    draw(x,y,sig,16,'ground',fceil=np.max(abs(sig)),cmap='coolwarm',title=til)
    #pl.savefig('tmp/corr')
    pl.pause(2)
    sig = stretch(fv-f)
    til = ('(corr diff)^(1/4) at %i nm' % (1e9*lam+0.5))
    
    draw(x,y,sig,16,'ground',fceil=np.max(abs(sig)),cmap='coolwarm',title=til)
    
    for fr in range(48):
        jd = 246e5 + (fr-7)/48
        u,v,w = source.get_uvw(jd,dx,dy,0*dx)
        hx,hy,hz = source.get_xyz(jd,0,0,1)
        print('time %.3f altitude %.3f' % (jd-246e5,hz))
        if hz > 0:
            pl.plot(u,v,'+',color='white')
        print(fr)
    pl.pause(10)

        

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

