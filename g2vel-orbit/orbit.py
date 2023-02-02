import numpy as np
from scipy.optimize import brentq

class Orbit():

    def __init__(self,P,e,I,Omega,omega,jdperi,q):
        self.P, self.e, self.jdperi, self.q = P, e, jdperi, q
        pi, cos, sin = np.pi, np.cos, np.sin
        x = I*pi/180
        self.cosI = cos(x)
        x = Omega*pi/180
        self.cosOm, self.sinOm = cos(x), sin(x)
        x = omega*pi/180
        self.cosom, self.sinom = cos(x), sin(x)
        
    def pos(self,tph):
        ''' Keplerian x,y,z '''
        e = self.e
        # First solve Kepler's equation
        # tph is orbital phase (aka mean anomaly)
        psi = brentq(lambda psi: psi - e*np.sin(psi) - tph, 0, 2*np.pi)
        # then find x,y in orbital plane
        ec = (1-e*e)**.5
        x,y = np.cos(psi)-e, ec*np.sin(psi)
        # next rotate by omega
        cs, sn = self.cosom, self.sinom
        x,y = x*cs - y*sn, x*sn + y*cs
        # then inclination
        x,y = x, y*self.cosI
        # and finally rotate by Omega
        cs,sn = self.cosOm, self.sinOm
        x,y = x*cs - y*sn, x*sn + y*cs
        return x,y

    def binarypos(self,jd):
        q = self.q
        a1 = q/(1+q)
        a2 = -1/(1+q)
        tph = 2*np.pi*(jd-self.jdperi) / self.P
        tph = tph % (2*np.pi)
        xs,ys = self.pos(tph)
        return [a1*xs,a2*xs], [a1*ys,a2*ys]


