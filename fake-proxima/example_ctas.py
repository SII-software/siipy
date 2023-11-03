import numpy as np
import matplotlib.pyplot as pl
from visib import sbright, correldens
from graphics import draw

from fakeprox import teff
from layout import Layout

pl.style.use('dark_background')

source = Layout()

source.set_radec((14,29,42.94),(-62,-40,-46.13)) # Proxima Cen
source.set_latlon((-24,-41,-0.34),(-70,-18,-58.84))  # CTA-S

def stretch(f):
    return np.sign(f) * np.abs(f)**(1/4)

def illus(fname):
    pl.pause(1)
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
    fv = correldens(S*vari,lam)
    f /= np.max(f)
    fv /= np.max(fv)
    til = ('photons / (m^2 s Hz mas^2)')
    masf = (np.pi/180/3.6e6/ds)**2
    smax = np.max(S*vari*masf)
    #draw(sx,sy,S*masf,8,'sky',ceil=smax,cmap='magma',title=til)
    #illus('figs/src_ud')
    draw(sx,sy,S*vari*masf,8,'sky',ceil=smax,cmap='magma',title=til)
    illus('figs/src')
    #til = ('|visib|^2 at %i nm' % (1e9*lam+0.5))
    #draw(x,y,f,gz,'ground',fceil=1,cmap='coolwarm',title=til)
    #illus('figs/visibq')
    #sig = stretch(f)
    #til = ('|visib| at %i nm' % (1e9*lam+0.5))
    #draw(x,y,sig,gz,'ground',fceil=1,cmap='coolwarm',title=til)
    #illus('figs/visib')
    sig = stretch(fv-f)
    til = ('(diff corr)^(1/4) at %i nm' % (1e9*lam+0.5))
    
    for fr in range(48):
        jd = 246e5 + fr/48
        u,v,w = source.get_uvw(jd,dx,dy,0*dx)
        hx,hy,hz = source.get_xyz(jd,0,0,1)
        print('time %.3f altitude %.3f' % (jd-246e5,hz))
        draw(x,y,sig,gz,'ground',fceil=np.max(abs(sig)),cmap='coolwarm',title=til)
        if hz > 0:
            pl.plot(u,v,'+',color='white')
        print(fr)
        illus('tmp/dcorr%i'%fr)
        

x = []
y = []
fil = open('locs-ctasold.txt')
while True:
    s = fil.readline()
    if not s:
        break
    xs,ys,c = s.split()
    if c == 'MST' or c == 'LST':
        x.append(float(xs))
        y.append(float(ys))

dx = []
dy = []
for i in range(len(x)):
    for j in range(len(x)):
        if i != j:
            bx = x[i]-x[j]
            by = y[i]-y[j]
            dx.append(bx)
            dy.append(by)
dx = np.array(dx)
dy = np.array(dy)



run()

