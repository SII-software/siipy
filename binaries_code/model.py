from scipy.interpolate import RectBivariateSpline
from numpy.fft import fft2, ifft2, fftshift
from dynesty import DynamicNestedSampler
from multiprocessing import Pool
from orbitfit import fit_orb
import scipy.special as sp
from aperture import aper
import numpy as np
import pickle

ap = aper()
ob = fit_orb()

class modeling():
      """
      Creates an object to estimate the parameters of the binary star systems.
      The binaries are close binary (like Spica), Regor type binary, and wide binary.
      The model of stars can be uniform or linear limb-darkening depending on the linear limb-darkening coefficient.
      """
      # defined constant
      c = 299792458.0       # in m/s
      h = 6.62607015e-34    # J.s
      k = 1.380649e-23      # J/K

      def star_para(self, Tmp_a, Tmp_b, Tmp_r, dis_e, lam_c, lam_e):
          """
          The set of parameters of binary star system.
          
          Parameters :
          ----------
          Tmp_a : int
                  The temperature of source A.
          Tmp_b : int
                  The temperature of source B.
          Tmp_r : int
                  The temperature of Roche-lobe region.
          dis_e : float
                  The distance of star system from earth.
          lam_c : float
                  The wavelength of continuum line observation.
          lam_e : float
                  The wavelength of emission line observation.
                
          Returns :
          -------
                  None.
          """
          self.Tmp_a = Tmp_a
          self.Tmp_b = Tmp_b
          self.Tmp_r = Tmp_r
          self.dis_e = dis_e
          self.lam_c = lam_c
          self.lam_e = lam_e

      def sii_para(self, Area, delt, itta, loss, avgt, chnl):
          """
          The set of parameters of SII observatory.
          
          Parameters :
          ----------
          Area : float
                 Area of telescope aperture.
          delt : float
                 The resolution time of photon detector.
          itta : float
                 The quantum efficiency of  photon detector.
          loss : float
                 Loss in system.
          avgt : float
                 The average time of observation.
          chnl : int
                 The number of channel in observation.
             
          Returns :
          -------
                 None.
          """
          self.Area = Area
          self.delt = delt
          self.itta = itta
          self.loss = loss
          self.avgt = avgt
          self.chnl = chnl

      # photon flux from single source
      def pho_spec(self, R, T, lam):
          """
          """
          a1 = np.pi * R**2
          a1 /= (self.dis_e*lam)**2
          a2 = self.h*self.c/(lam*self.k*(T + 2.725))
          a2 = np.clip(a2,0,10)
          I = a1 / (np.exp(a2) - 1)
          return I

      def hbt(self, xcord, ycord, Xcord, Ycord, Rad_a, Rad_b, Rad_r, lin_c): 
          """
          The square visibility model of a binary star system with linear limb-darkeing coefficient.
          
          Parameters :
          ----------
          xcord : N x N 
                  For the baseline along x-direction.
          ycord : N x N 
                  For the baseline along y-direction.
          Xcord : int
                  X-coordinate of orbit of binary on observational time of SII. 
          Ycord : int
                  Y-coordinate of orbit of binary on observational time of SII. 
          Rad_a : float
                  The radius of source A.
          Rad_b : float
                  The radius of source B.
          Rad_r : float
                  The maximum radius of Roche lobe from the surface of star.
          lin_c : float
                  The linear limb-darkening coefficients.
          
          Returns :
          -------
                  The squared visibility data.
          """
          r = np.sqrt(xcord**2 + ycord**2)
          v = 2*np.pi/(self.lam_c*self.dis_e)

          b1 = v*r*Rad_a
          b2 = v*r*Rad_b
          b3 = v*r*Rad_r

          C = 1/2 - lin_c/6
          a = 1 - lin_c
          b = lin_c          
                       
          I_a = self.pho_spec(Rad_a, self.Tmp_a, self.lam_c)
          I_b = self.pho_spec(Rad_b, self.Tmp_b, self.lam_c)

          V_a = a*sp.jv(1, b1)/b1
          V_a += b*np.sqrt(np.pi/2)*sp.jv(3/2,b1)/b1**(3/2)
          V_a *= I_a
          
          V_1 = a*sp.jv(1, b2)/b2
          V_1 += b*np.sqrt(np.pi/2)*sp.jv(3/2,b2)/b2**(3/2)
          V_1 *= I_b

          if self.lam_e == self.lam_c:
             I_r = self.pho_spec(Rad_r, self.Tmp_r, self.lam_e)
             I_br = self.pho_spec(Rad_b, self.Tmp_r, self.lam_e)
             V_2 = 2 * I_br * sp.jv(1, b2)/b2                        # limb-darkening of Wolf Rayet star is not assumed
             V_3 = 2 * I_r * sp.jv(1, b3)/b3
          else:
             I_r = I_br = 0
             V_2 = V_3 = 0
          
          V_b = V_1 - V_2 + V_3

          V_ab = V_a * V_b * np.cos(v * (xcord*Xcord + ycord*Ycord))
          V = (V_a**2 + V_b**2 + 2 * V_ab)
          V /= (I_a*C + I_b*C - I_br + I_r)**2    
          return V

      def noise(self, Rad_a, Rad_b, Rad_r):
          """
          The noise calculation for a binary star system with SII.
          
          Parameters :
          ----------
          
          Rad_a : float
                  The radius of source A.
          Rad_b : float
                  The radius of source B.
          Rad_r : float
                  The maximum radius of Roche lobe from the surface of star.
              
          Returns :
          -------
                  float.
          """
          Phi = self.pho_spec(Rad_a, self.Tmp_a, self.lam_c)
          Phi += self.pho_spec(Rad_b, self.Tmp_b, self.lam_c)
          if self.lam_e == self.lam_c:
             Phi += self.pho_spec(Rad_r, self.Tmp_r, self.lam_e)
             Phi -= self.pho_spec(Rad_b, self.Tmp_r, self.lam_e)
          else:
             Phi += 0
          tspec = self.Area * self.itta * self.loss * Phi
          return np.sqrt(self.delt/self.avgt)/(tspec*self.chnl)

      def r_signal(self, xcord, ycord, zs_xy, tcord, Rad_a, Rad_b, Rad_r, lin_c):
          """
          """
          X, Y = ob.polyfunc(zs_xy, tcord)
          return self.hbt(xcord, ycord, X, Y, Rad_a, Rad_b, Rad_r, lin_c)

      def c_signal(self, xcord, ycord, Rad_a, Rad_b, Rad_r, lin_c, P, a, e, inc, Omega, epoch, omega, jlday):
          """
          """
          X, Y = ob.close_xy(jlday, P, a, e, inc, Omega, epoch, omega)
          return self.hbt(xcord, ycord, X, Y, Rad_a, Rad_b, Rad_r, lin_c)

      # Fourier transform of any function 
      def fourier(self, f):
          """
          """
          return fftshift(fft2(fftshift(f)))
      # Inverse Fourier transform of any function
      def ifourier(self, f):
          """
          """
          return fftshift(ifft2(fftshift(f)))

      # getting signal with w aperture which is masked with some function 
      def w_signal(self, x, y, X, Y, Rad_a, Rad_b, w):
          """   
                 
          """
          g = self.hbt(x, y, X, Y, Rad_a, Rad_b, 0, 0)
          fg = self.fourier(g)
          fw = self.fourier(w)
          dfun = self.fourier(0*x + 1)
          return self.ifourier(fg*fw*fw)/self.ifourier(dfun*fw*fw)

      def smdata(self, bs, uv, sf, jd, nn, XY=None, apt=None, pl=None, kp=None):
          """
          The noisy simulated observed data for a wide binary, close binary and Regor type object.
          In wide binary, both star is following the uniform model.
          In close binary, both star is following the same linear limb-darkening model and the variable orbit is following the Kepler's equation.          
          In Regor type binary, the Wolf Rayet star is following the uniform model. However the companion star is following the linear limb-darkening model         
          and the arc of variable orbit of this binary is following the polynomial function of order N.
          Parameters :
          ----------
          bs : list
               List of baseline on each observational time in xy direction i.e. bs = [baseX, baseY].
               baseX : float
                       The rotation of baseline along x-axis.
               baseY : float
                       The rotation of baseline along y-axis.
          uv : list
               List of variational xy plane i.e. uv = [xcord, ycord].
               1. xcord : N x N 
                          For the baseline along x-direction.
               2. ycord : N x N 
                          For the baseline along y-direction.
          sf : list
               List of radius of objects and linear limb-darkening coefficient i.e. sf = [Rad_a, Rad_b, Rad_r, lin_c].
               1. Rad_a : float
                          The radius of source A.
               2. Rad_b : float
                          The radius of source B.
               3. Rad_r : float
                          The maximum radius of Roche lobe from the surface of star.
               4. lin_c : float
                          The linear limb-darkening coefficients.
          jd : array
               The array of julian day of all observation. 
          nn : int
               nn = 1 for simulated noisy data and nn = 0 for model data.
          XY : list
               List of XY plane of orbit of binary i.e. [Xcord, Ycord].
               1. Xcord : int
                          X-coordinate of orbit of binary on observational time of SII. 
               2. Ycord : int
                          Y-coordinate of orbit of binary on observational time of SII.
               Note. XY = XY for wide binary, otherwise XY = None.
          apt: list
               List of parameters of masked aperture i.e. ap = [aperX, aperY, radii, width, orien].
               1. aperX : N x N 
                          For the aperture along x-direction.
               2. aperY : N x N 
                          For the aperture along y-direction.
               3. radii : int or float
                          Radius of telescopes.
               4. width : float
                          The width of masking strips.
               5. orien : float
                          The orientation of masking strips in radian.
               Note. apt = apt for wide binary, otherwise apt = None.
          pl : list
               List of polynomial coefficients and time coordinates of orbit i.e. pl = [zs_xy, tcord].
               1. zs_xy : list
                          The list of polynomial coefficient in x and y order (ax[i], ay[i]).
               2. tcord : array
                          t-coordinate of orbit of binary on observational time of SII.
               Note. pl = pl for Regor type binary, otherwise pl = None.
          kp : list
               List of orbital element of binary i.e. [cycle, se_ax, eccen, incli, Omega, epoch, omega].
               1. cycle : float
                          The period of an orbit.
               2. se_ax : float
                          The semi-major axis of an orbit.
               3. eccen : flaot
                          The eccentricity of an orbit.
               4. incli : float
                          The value of inclination of an orbit.
               5. Omega : float
                          The value of Longitude of the ascending node of an orbit.
               6. epoch : float
                          The value of epoch of an orbit
               7. omega : float
                          The value of Argument of Periapsis of an orbit.
               Note. kp = kp for close binary, otherwise kp = None.
          Returns :
          -------
               array : noisy simulated observed data. 
          """
          if apt is not None:
             w = ap.aper(apt[0], apt[1], apt[2], apt[3], apt[4])
             g = self.w_signal(uv[0], uv[1], XY[0], XY[1], sf[0], sf[1], w)
          elif kp is not None:
             g = self.c_signal(uv[0], uv[1], sf[0], sf[1], sf[2], sf[3], kp[0], kp[1], kp[2], kp[3], kp[4], kp[5], kp[6], jd)
          elif pl is not None:
             g = self.r_signal(uv[0], uv[1], pl[0], pl[1], sf[0], sf[1], sf[2], sf[3])
          sdata = []
          for i in range(len(bs[0])):
              gi = g[i,:,:]
              xgr = uv[0][i,0,:]
              ygr = uv[1][i,:,0]
              zi = np.transpose(abs(gi))
              spline = RectBivariateSpline(xgr,ygr,zi)
              sdata.append(spline.ev(bs[0][i,:], bs[1][i,:]) + nn * self.noise(sf[0], sf[1], sf[2]) * np.random.randn(len(jd)))
          return np.array(sdata)

      # define the transformation of prior for parameters to be estimated
      def prior(self, theta): 
          """
          """ 
          p_p = theta
          p = []
          for i in range(len(self.pnames)): 
              p.append(p_p[i] * (self.p_max[i] - self.p_min[i]) + self.p_min[i])
          return tuple(p)

      # define the liklihood for parameters to be estimated
      def likli(self, ln):
          """
          """
          p = ln
          sigma = self.noise(self.Rad_a, self.Rad_b, self.Rad_r)
          if self.Rad_r == 0:
             w = ap.aper(self.aperX, self.aperY, self.radii, self.width, self.orien)
             g = self.w_signal(self.xcord, self.ycord, p[0], p[1], self.Rad_a, self.Rad_b, w)

          elif self.Rad_r == 1e-25:
             if self.para == 2:
                g = self.c_signal(self.xcord, self.ycord, p[0], p[1], self.Rad_r, self.lin_c, self.cycle, self.se_ax, self.eccen, self.incli, self.Omega, self.epoch, self.omega, self.jlday)
             elif self.para == 3:
                g = self.c_signal(self.xcord, self.ycord, p[0], p[1], self.Rad_r, p[2], self.cycle, self.se_ax, self.eccen, self.incli, self.Omega, self.epoch, self.omega, self.jlday)
             elif self.para == 4:
                g = self.c_signal(self.xcord, self.ycord, p[0], p[1], self.Rad_r, p[2], self.cycle, p[3], self.eccen, self.incli, self.Omega, self.epoch, self.omega, self.jlday)
             elif self.para == 5:
                g = self.c_signal(self.xcord, self.ycord, p[0], p[1], self.Rad_r, p[2], self.cycle, p[3], p[4], self.incli, self.Omega, self.epoch, self.omega, self.jlday)
             elif self.para == 6:
                g = self.c_signal(self.xcord, self.ycord, p[0], p[1], self.Rad_r, p[2], self.cycle, p[3], p[4], p[5], self.Omega, self.epoch, self.omega, self.jlday)
             elif self.para == 7:
                g = self.c_signal(self.xcord, self.ycord, p[0], p[1], self.Rad_r, p[2], self.cycle, p[3], p[4], p[5], p[6], self.epoch, self.omega, self.jlday)
             elif self.para == 8:
                g = self.c_signal(self.xcord, self.ycord, p[0], p[1], self.Rad_r, p[2], self.cycle, p[3], p[4], p[5], p[6], self.epoch, p[7], self.jlday)
             elif self.para == 9:
                g = self.c_signal(self.xcord, self.ycord, p[0], p[1], self.Rad_r, p[2], p[8], p[3], p[4], p[5], p[6], self.epoch, p[7], self.jlday)
             else:
                g = self.c_signal(self.xcord, self.ycord, p[0], p[1], self.Rad_r, p[2], p[8], p[3], p[4], p[5], p[6], p[9], p[7], self.jlday)

          else:
             if self.lam_e == self.lam_c:
                if self.para == 2:
                   g = self.r_signal(self.xcord, self.ycord, self.zs_xy, self.tcord, p[0], p[1], self.Rad_r, self.lin_c) # only both main radii with emission line
                elif self.para == 3:
                   g = self.r_signal(self.xcord, self.ycord, self.zs_xy, self.tcord, p[0], p[1], p[2], self.lin_c) # both radii and roche lobe with emission line
                else:
                   zs = p[3:]
                   g = self.r_signal(self.xcord, self.ycord, zs, self.tcord, p[0], p[1], p[2], self.lin_c)   # all radii together with polynomial orbit with emission line
             else:
                if self.para == 2:
                   g = self.r_signal(self.xcord, self.ycord, self.zs_xy, self.tcord, p[0], p[1], self.Rad_r, self.lin_c) # only both main radii with continuum line
                elif self.para == 3:
                   g = self.r_signal(self.xcord, self.ycord, self.zs_xy, self.tcord, p[0], p[1], self.Rad_r, p[2]) # only both main radii and ld of BS with continuum line
                else:
                   zs = p[3:]
                   g = self.r_signal(self.xcord, self.ycord, zs, self.tcord, p[0], p[1], 0, p[2])   # only both main radii and ld of BS together with polynomial orbit with continuum line    
          model_data = []
          for i in range(len(self.baseX)):
             gi = g[i,:,:]
             xgr = self.xcord[i,0,:]
             ygr = self.ycord[i,:,0]
             zi = np.transpose(abs(gi))
             spline = RectBivariateSpline(xgr,ygr,zi)
             model_data.append(spline.ev(self.baseX[i,:], self.baseY[i,:]))       
          mdata = np.array(model_data)    
          sdata = np.load(self.dfname)              
          G = np.sum(sigma**(-2) * sdata * mdata)                                 
          W = np.sum(sigma**(-2) * mdata * mdata) 
          return 0.5 * (G*G/W - np.log(W)) 

      def estimate(self, fname, bs, uv, sf, p_min, p_max, para, apt=None, pl=None, kp=None, jlday=None):
          """
          Estimate the parameters of binary star system.
          For close binary, the order is [Rad_a, Rad_b, lin_c, se_ax, eccen, incli, Omega, omega, cycle, epoch].
          For Regor type binary, the order is [Rad_a, Rad_b, Rad_r/lin_c, zs_xy]. Rad_r with emission line and l with continuum line.
          For wide binary, the order is [Xcord, Ycord].
          
          Parameters :
          ----------
          fname : list
                  List of string i.e. fname = [fname2, dfname, pnames]
                  1. fname2 : str
                              The output .pkl file.
                  2. dfname : str
                              The input simulated or observed noisy data containing the HBT signal for a baseline.
                  3. pnames : tupple
                              Tupple of name of estimated parameters.
          bs    : list
                  List of baseline on each observational time in xy direction i.e. bs = [baseX, baseY].
                  baseX : float
                          The rotation of baseline along x-axis.
                  baseY : float
                          The rotation of baseline along y-axis.
          uv    : list
                  List of variational xy plane i.e. uv = [xcord, ycord].
                  1. xcord : N x N 
                             For the baseline along x-direction.
                  2. ycord : N x N 
                             For the baseline along y-direction.
          sf    : list
                  List of radius of objects i.e. sf = [Rad_a, Rad_b, Rad_r, lin_c].
                  1. Rad_a : float
                             The radius of source A.
                  2. Rad_b : float
                             The radius of source B.
                  3. Rad_r : float
                             The maximum radius of Roche lobe from the surface of star.
                  4. lin_c : float
                             The linear limb-darkening coefficients.
          p_min : list
                  The list of minimum value of all estimated parameters.
          p_max : list
                  The list of maximum value of all estimated parameters.
          para  : int
                  It decides the number of parameters to be estimated.
                  i.e. para = 2, only [Rad_a, Rad_b].
          apt   : list
                  List of parameters of masked aperture i.e. ap = [aperX, aperY, radii, width, orien].
                  1. aperX : N x N 
                             For the aperture along x-direction.
                  2. aperY : N x N 
                             For the aperture along y-direction.
                  3. radii : int or float
                             Radius of telescopes.
                  4. width : float
                             The width of masking strips.
                  5. orien : float
                             The orientation of masking strips in radian.
                  Note. It is used only for wide binary, otherwise None.
          pl    : list
                  List of polynomial coefficients, t-coordinates and wavelengths of orbit i.e. pl = [zs_xy, tcord, lam_c, lam_e].
                  1. zs_xy : list
                             The list of polynomial coefficient in x and y order (ax[i], ay[i]).
                  2. tcord : array
                             t-coordinate of orbit of binary on observational time of SII.
                  Note. It is used only for Regor type binary, otherwise None.
          kp    : list
                  List of orbital element of binary and jd i.e. [cycle, se_ax, eccen, incli, Omega, epoch, omega].
                  1. cycle : float
                             The period of an orbit.
                  2. se_ax : float
                             The semi-major axis of an orbit.
                  3. eccen : flaot
                             The eccentricity of an orbit.
                  4. incli : float
                             The value of inclination of an orbit.
                  5. Omega : float
                             The value of Longitude of the ascending node of an orbit.
                  6. epoch : float
                             The value of epoch of an orbit
                  7. omega : float
                             The value of Argument of Periapsis of an orbit.
                  Note. It is used only for close binary, otherwise None.
          jlday : array
                  The array of julian day of all observation.
          Returns :
          -------
                  Return the .pkl file which contains information about the posterior distribution of parameters.
          """
          self.baseX, self.baseY, self.xcord, self.ycord = bs[0], bs[1], uv[0], uv[1]
          self.Rad_a, self.Rad_b, self.Rad_r, self.lin_c = sf[0], sf[1], sf[2], sf[3]
          self.p_max, self.p_min, self.para = p_max, p_min, para
          self.dfname, self.pnames = fname[1], fname[2]
          if apt is not None:
             self.aperX, self.aperY, self.radii, self.width, self.orien = apt[0], apt[1], apt[2], apt[3], apt[4]
          elif pl is not None:
             self.zs_xy, self.tcord = pl[0], pl[1]
          elif kp is not None:
             self.cycle, self.se_ax, self.eccen, self.incli = kp[0], kp[1], kp[2], kp[3]
             self.Omega, self.epoch, self.omega, self.jlday = kp[4], kp[5], kp[6], jlday  
          ndim = len(self.pnames)
          pool = Pool()
          sampler = DynamicNestedSampler(self.likli, self.prior, ndim, pool=pool, queue_size=28)
          sampler.run_nested()
          dres = sampler.results
          out = open(fname[0],'wb')
          pickle.dump(dres,out)
          out.close()
          pool.close()

