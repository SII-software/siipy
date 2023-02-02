import numpy as np
from numpy.random import seed,random
from fourier import fourier, ifourier

SI_c = 299792458
SI_h = 6.626070040e-34
SI_k = 1.38064852e-23

# Generates amplitude map from T map
def ampl(T,lam,ds):
    c = SI_c
    h = SI_h
    k = SI_k
    z = h*c/(lam*k*(T+2.725))  # to avoid overflow in exp(z)
    z = np.clip(z,0,20)
    return ds/lam * (np.exp(z)-1)**-.5

# Computes exposure X(t,dt) at wavelength lam
def exposure(S,t,dt,lam):
    N = S.shape[0]
    seed(0)
    rr = random(N**2).reshape((N,N)) - 0.5
    tpdnu = np.tan(rr*np.pi)
    X = 0*S
    Nav = 20
    Nsamp = int(dt*Nav)+1
#    print('Nsamp = ',Nsamp)
    for ts in np.linspace(t,t+dt,Nsamp):
#        print(ts,end=' ')
        ph = tpdnu*ts
        St = S*np.exp(1j*ph)
        V = fourier(St)
        X += abs(V*V)
#    print()
    X *= dt/Nsamp
    c = SI_c
    h = SI_h
    pflux = 2*c/lam * np.sum(X)/(N**2*dt)
    mag = -2.5*(np.log10(h*pflux)+22.44)
    print('%9.2e photons m^{-2} s^{-1} (ln Hz)^{-1} (unpolarized)  AB = %5.2f' % (pflux,mag))
    return X

# Computes maps of HBT correlation g and SNR per data point
def stats(X,msk=None):
    N = X.shape[0]
    Xav = np.sum(X)/N**2
    if 'NoneType' not in str(type(msk)):
        X *= msk
    X1 = X.reshape(N*N)
    N2 = len(X1[X1>0])
    C = abs(ifourier(abs(fourier(X)**2)))/N2
    g = C/Xav**2 - 1
    gmax = np.max(g)
    f = (C - Xav**2) / Xav
    snrmax = np.max(f)
    g[g<0] = 0
    print('Xav = %9.2e m^{-2}  Xmax/Xav = %5.2f' % (Xav,np.max(X)/Xav))
    print('max g = %9.2e  max SNR = %9.2e' % (gmax,snrmax))
    return g

