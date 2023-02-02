import numpy as np
import matplotlib.pyplot as pl
import matplotlib.colors as colors

mas = 1e-3/(180/np.pi*3600)

# Draws contour map of f (on sky or ground, directed)
# zoom factor (power of 2 preferred)
def draw(xc,yc,f,zoom,where,cmap='Greys_r',ceil=None,fceil=None,title=None):
    def cen(f):
        N = f.shape[0]
        M = N//(zoom*2)
        return f[N//2-M:N//2+M,N//2-M:N//2+M]
    f = cen(f)
    pl.clf()
    fmax = f.max()
    fmin = f.min()
    if where=='sky':
        sx,sy = cen(xc)/mas, cen(yc)/mas
        if ceil:
            fmin,fmax = 0,max(ceil,f.max())
        levs = np.linspace(fmin,fmax,40)
        cs = pl.contourf(sx,sy,f,levs,cmap=cmap)
        pl.xlabel('mas')
    if where=='ground':
        if xc[-1,-1] > 3e4:
            x,y = 1e-3*cen(xc), 1e-3*cen(yc)
            pl.xlabel('kilometres')
        elif xc[-1,-1] < 1:
            x,y = 1e3*cen(xc), 1e3*cen(yc)
            pl.xlabel('millimetres')
        else:
            x,y = cen(xc), cen(yc)
            pl.xlabel('metres')
        if ceil:
            fmin,fmax = 0,max(ceil,f.max())
            levs = np.linspace(fmin,fmax,20)
        elif fceil:
            fmax = max(fceil,f.max())
            fmin = -fmax
            levs = np.linspace(fmin,fmax,80)
        else:
            fmin,fmax = 0,f.max()
            levs = np.linspace(fmin,fmax,20)
        cs = pl.contourf(x,y,f,levs,norm=colors.Normalize(vmin=fmin,vmax=fmax),cmap=cmap)
    if fmax > 1000:
        fms = '%5.1e'
    elif 10 < fmax:
        fms = '%i'
    else:
        lgf = np.log10(fmax)
        ip = int(-lgf) + 2
        if lgf < -5:
            fms = '%7.1e'
        else:
            fms = '%' + '.%i' % ip + 'f'
    pl.colorbar(cs,format=fms)
    if title:
        pl.title(title)
    pl.gca().set_aspect('equal')
    pl.tight_layout()

  
