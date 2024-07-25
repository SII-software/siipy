from dynesty.utils import resample_equal
from scipy.optimize import leastsq
import scipy.optimize as sciopt
import matplotlib.pyplot as plt
from obstime import observe
from orbit import orb_int
import numpy as np
import pickle

nt = observe()
o = orb_int() 

class fit_orb():
      """
      Create an object to observe the orbit on observational time of SII.
      Fitting of the complete close orbit can be done with Kepler's equation.
      Fitting of an arc of orbit can be done with a polynomial function. 
      Or the sum of polynomial fuction and oscillatory function.
      """      
      def orbit(self, start, ended, obj1, obj2, fname1, fname2, fname3):
          """
          The respective position (orbit) of stars on observational time for SII.
          
          Parameters :
          ----------
          
          start  : array
                   Starting time on each observational interval.
          ended  : array
                   Ending time on each observational interval.
          obj1   : int
                   The first object of binary in star system.
          obj2   : int
                   The second object of binary in star system.
          fname1 : str 
                   The name of .npy file containing X-coordinate of each object in system.
          fname2 : str 
                   The name of .npy file containing Y-coordinate of each object in system.
          fname3 : str 
                   The name of .npy file containing t-coordinate of each object in system.  
                 
          Returns :
          -------
                   Return the X-coordinate, Y-coordinate, and t-coordinate of orbit of binary on observational time of SII.
          """
          x = np.load(fname1)
          y = np.load(fname2)
          t = np.load(fname3)

          dx = x[obj1,::1] - x[obj2,::1] 
          dy = y[obj1,::1] - y[obj2,::1]

          Xth = []
          Yth = []
          th = []
          for i in range(0, len(start)):
              Xth.append(dx[start[i] : ended[i]])
              Yth.append(dy[start[i] : ended[i]])
              th.append(t[start[i] : ended[i], 0])

          return np.concatenate(Xth), np.concatenate(Yth), np.concatenate(th)

      # variable position of close orbit (X, Y)
      def close_xy(self, jul, P, a, e, inc, Omega, epoch, omega): 
          """
          """   
          x_t = []
          y_t = []
          for i in jul:
              ec = (1-e**2)**.5
              mean_an = 2*np.pi * (i-epoch)/P
              mean_an = mean_an % (2*np.pi)
              si = sciopt.brentq(lambda psi: psi - e*np.sin(psi) - mean_an, 0, 2*np.pi)
              x, y = np.cos(si)-e, ec*np.sin(si)
              pos = a * o.rotate(x, y, omega, inc, Omega)
              x_t.append(pos[0])
              y_t.append(pos[1])                                                     
          xt = np.array(x_t)
          yt = np.array(y_t)                                                 
          return xt, yt

      # the residuals
      def resi(self, zs, jul, X, Y):
          """
          """
          xt, yt = self.close_xy(jul, zs[0], zs[1], zs[2], zs[3], zs[4], zs[5], zs[6])
          return np.concatenate((X-xt, Y-yt))

      def corb_fit(self, zs_xy, jlday, Xcord, Ycord, fname): 
          """
          Fit the observed orbit of a close binary using Kepler's equation.
          
          Parameters :
          ----------
          
          zs_xy : list
                  The orbital parameters in order [P (days), a (au), e, I (deg), Omega (deg), epoch (days), omega (deg)].
          jlday : array
                  The array of julian day of all observation.
          Xcord : array
                  The array of X-coordinate of binary orbit on observational time.
          Ycord : array
                  The array of Y-coordinate of binary orbit on observational time. 
          fname : str
                  The output .png file.
              
          Returns :
          -------
                  It return the .png file which shows the fitted orbit with observed orbit and residuals.  
          """
          res = leastsq(self.resi, zs_xy, (jlday, Xcord, Ycord))                
          zs_xy = res[0]                                      
          z_res = self.resi(zs_xy, jlday, Xcord, Ycord)
          dx1 = z_res[:len(jlday)]                             
          dy1 = z_res[len(jlday):]       # find residuals dx1, dy1
          print('Period in days = ', zs_xy[0])
          print('Semi-major axis in light second = ', zs_xy[1])
          print('Eccentricity = ', zs_xy[2]) 
          print('Inclination in radian = ', zs_xy[3]) 
          print('Longitude of the ascending node in radian = ', zs_xy[4])
          print('Epoch in days = ', zs_xy[5])
          print('Argument of periapsis in radian = ', zs_xy[6])
          x_t, y_t = self.close_xy(jlday, zs_xy[0], zs_xy[1], zs_xy[2], zs_xy[3], zs_xy[4], zs_xy[5], zs_xy[6])
          # plot the residuals and fitted orbit with observed orbit
          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [12,12]
          plt.rcParams["font.weight"] = "bold"
          plt.rcParams["axes.labelweight"] = "bold"
          plt.plot(dx1, dy1, '.', label='residuals')
          plt.plot(Xcord, Ycord, '.', label='observed orbit')
          plt.plot(x_t, y_t, '.', label='fited orbit')
          plt.xlabel('X-coordinate in light second')
          plt.ylabel('Y-coordinate in light second')
          plt.title('Compare the fiting orbit together with residuals', fontweight='bold')
          plt.legend(loc="lower right")
          plt.savefig(fname + ".png")
          plt.show()
 

      # fitting the variable position of arc of an orbit (X, Y)
      # polynomial order function
      def polyfunc(self, zs, t):
          """
          """
          x_t = y_t = 0
          tp = 1
          for p in range(0, len(zs), 2):
              x_t += zs[p] * tp
              y_t += zs[p+1] * tp
              tp *= t
          return x_t, y_t

      def polynom(self, zs, t, X, Y):
          """
          """
          x_t, y_t = self.polyfunc(zs, t)
          return np.concatenate((X-x_t, Y-y_t))

      def resi_poly(self, zs_xy, tcord, Xcord, Ycord):
          """
          Fit the variable position of an arc of an orbit (X, Y).
          
          The orbit is fitted with Nth order polynomial function.
          
          Parameters :
          ----------
          
          zs_xy : list
                  The list of 2(N+1) zeros for Nth order polynomial.
          tcord : array
                  The array of t-coordinate of binary orbit on observational time.
          Xcord : array
                  The array of X-coordinate of binary orbit on observational time.
          Ycord : array
                  The array of Y-coordinate of binary orbit on observational time. 
              
          Returns :
          -------
          dx, dy  : array
                    The X and Y residuals of orbit.
                    It return the .png file which shows the fitted orbit with observed orbit and residuals.          
          """
          res = leastsq(self.polynom, zs_xy, (tcord, Xcord, Ycord))
          zs_xy = res[0]
          z_res = self.polynom(zs_xy, tcord, Xcord, Ycord)
          dx1 = z_res[:len(tcord)]
          dy1 = z_res[len(tcord):]
          x_t, y_t = self.polyfunc(zs_xy, tcord)
          print('The polynomial coefficients (ax[i], ay[i]) for the curvature of orbit', zs_xy)
          # plot the residuals and fitted orbit with observed orbit
          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [12, 12]
          plt.rcParams["font.weight"] = "bold"
          plt.rcParams["axes.labelweight"] = "bold"
          plt.plot(dx1, dy1, '.', label='residuals')
          plt.plot(Xcord, Ycord, '.', label='observed orbit')
          plt.plot(x_t, y_t, '.', label='fited orbit')
          plt.xlabel('X-coordinate in light second')
          plt.ylabel('Y-coordinate in light second')
          plt.title('Compare the fiting orbit together with residuals', fontweight='bold')
          plt.legend(loc="upper right")
          plt.savefig("fitorb.png")
          plt.show()
          return dx1, dy1

      # sinusoidal function to fit the planet orbit
      def sinusoid(self, os, t, X, Y):
          """
          """
          w = os[0]
          cs = np.cos(w*t)
          sn = np.sin(w*t)
          x_t = os[1]*cs + os[2]*sn
          y_t = os[3]*cs + os[4]*sn
          return np.concatenate((X-x_t, Y-y_t))

      # oscillatory residual 
      def resi_osl(self, os, t, dx, dy):
          """
          """
          ores = leastsq(self.sinusoid, os, (t, dx, dy))
          os = ores[0]
          yr = 86400*365.25                                                          
          per = 2*np.pi/os[0]/yr
          print('periodic part %6.3f yr' % per)
          z_ores = self.sinusoid(os, t, dx, dy)
          dx2 = z_ores[:len(t)]
          dy2 = z_ores[len(t):]
          return dx2, dy2

      def worb_fit(self, dfname, zs_xy, os_xy, tcord, Xcord, Ycord, title, fname):
          """
          Fit the variable position of an arc of an orbit (X, Y).
          
          The orbit is fitted with sum of Nth order polynomial and oscillatroy function.
           
          Parameters :
          ----------

          dfname : str
                   The input interferometric (the posterior distribution) .pkl file.
          zs_xy  : list
                   The list of 2(N+1) zeros for Nth order polynomial.
          os_xy  : list
                   The list of 5 initial parameters. First represents the angular frequency of oscillation of planet.
                   Other four represents the coefficient of cosine and sine function for X and Y orbit. 
          tcord  : array
                   The array of t-coordinate of binary orbit on observational time.
          Xcord  : array
                   The array of X-coordinate of binary orbit on observational time.
          Ycord  : array
                   The array of Y-coordinate of binary orbit on observational time. 
          title  : str
                   The title of output. It can be for different masses as well for different orbit.
          fname  : str
                   The output .png file.
              
          Returns :
          -------
                   Return .png file showing the residuals of binary orbit and planet orbit with interferometric error.          
          """
          infile = open(dfname, 'rb')
          results = pickle.load(infile)
          infile.close()
          weights = np.exp(results['logwt'] - results['logz'][-1]) # weighting the each nested 
          postsamples = resample_equal(results.samples, weights)   # posterior distribution after weighting each nested 

          # the 1-sigma standard deviation for x and y 
          p_x = np.percentile(postsamples[:,0],[16,50,84])
          p_y = np.percentile(postsamples[:,1],[16,50,84])

          err_x = p_x[2] - p_x[1]
          err_y = p_y[2] - p_y[1]

          dx1, dy1 = self.resi_poly(zs_xy, tcord, Xcord, Ycord)
          dx2, dy2 = self.resi_osl(os_xy, tcord, dx1, dy1)

          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [12,12]
          plt.rcParams["font.weight"] = "bold"
          plt.rcParams["axes.labelweight"] = "bold"
          plt.plot(dx1, dy1, '.-.', color = 'g', markersize='10', label='planet orbit')
          plt.plot(dx2, dy2, '--', color = 'r', markersize='10', label='signal residual')
          plt.errorbar(err_x, err_y, xerr = err_x, yerr = err_y, label='interferometric error', fmt='o', mfc='black', markersize='8', elinewidth=4, ecolor='black')
          plt.xlabel('X residual (light second)')
          plt.ylabel('Y residual (light second)')
          plt.title(title, fontweight='bold')
          plt.gca().set_aspect(aspect = 'equal', adjustable = 'datalim')
          plt.legend(loc="upper right")
          plt.savefig(fname + '.png')
          plt.show()


