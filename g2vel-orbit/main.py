import numpy as np
import matplotlib.pyplot as plt

from orbit import Orbit

spica = Orbit(P=78.53, e=0.334, I=65.5,
              Omega=247.7-90, omega=67.7, jdperi=2400000.5+50120.4, q=9/28)


#jdp = 2400000.5+50120.4
#for jd in np.linspace(jdp,jdp+40,41):
for jd in np.linspace(2460000,2460080,81):
    x,y = spica.binarypos(jd)
    R = [0.065,.023]
    plt.clf()
    plt.title('%.2f' % jd)
    ax = plt.gca()
    plt.xlim((-1.2,1.2))
    plt.ylim((-1.2,1.2))
    ax.set_aspect('equal')
    for k in 0,1:
        circ = plt.Circle((-x[k],-y[k]), R[k], color='r')
        ax.add_patch(circ)
    #plt.plot(x[0]-x[1],y[0]-y[1],'.')
    plt.pause(0.04)
