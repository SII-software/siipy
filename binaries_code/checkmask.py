from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt
from model import modeling
from aperture import aper
import numpy as np

ap = aper()
md = modeling()

class chkmsk():
      """
      Create an object to check the mask's width and orientation of a circular 
      aperture to observe signals from a wide binary.
      """
      def chbt(self, fname, uv, pa, XY, sf):
          """
          The plot of 2-d hbt signal to check the approx width and orientation of fringes.
          
          Parameters :
          ----------
          
          fname : str
                  The output .png file.
          uv    : list
                  List of variational xy plane i.e. uv = [xcord, ycord].
                  1. xcord : N x N 
                             For the baseline along x-direction.
                  2. ycord : N x N 
                             For the baseline along y-direction.
          pa    : list
                  The list of parameters of binary star system i.e. [Tmp_a, Tmp_b, Tmp_r, dis_e, lam_c, lam_e].
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
          XY    : list
                  List of XY plane of orbit of binary i.e. [Xcord, Ycord].
                  1. Xcord : int
                             X-coordinate of orbit of binary on observational time of SII. 
                  2. Ycord : int
                             Y-coordinate of orbit of binary on observational time of SII.
          sf    : list
                  List of radius of objects and linear limb-darkening coefficient i.e. sf = [Rad_a, Rad_b, Rad_r, lin_c].
                  1. Rad_a : float
                             The radius of source A.
                  2. Rad_b : float
                             The radius of source B.
                  3. Rad_r : float
                             The maximum radius of Roche lobe from the surface of star.
                  4. lin_c : float
                             The linear limb-darkening coefficients. 
            
          Returns :
          -------
                  Return the .png file of 2-d hbt signal for a given baseline.   
          """
          md.Tmp_a = pa[0]
          md.Tmp_b = pa[1]
          md.Tmp_r = pa[2]
          md.dis_e = pa[3]
          md.lam_c = pa[4]
          md.lam_e = pa[5]
          g = md.hbt(uv[0], uv[1], XY[0], XY[1], sf[0], sf[1], sf[2], sf[3])
          plt.close()
          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [12,12]
          plt.rcParams["font.weight"] = "bold" 
          plt.rcParams["axes.labelweight"] = "bold"
          plt.contourf(uv[0], uv[1], g)
          plt.title('HBT signal on grids of (x, y) plane for a baseline rotated along earth', fontweight='bold')
          plt.xlabel("East x in (meter)")
          plt.ylabel("North y in (meter)")
          plt.savefig(fname)
          plt.show()

      def aplot(self, fname, apt, lim):
          """
          Plot the masked aperture with specific radius, mask's width and orientation.
          
          Parameters :
          ----------
          
          fname : str
                  The output .png file.
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
          lim   : list
                  List of limiting value for X and Y axis i.e. lim = [xlim0, xlim1, ylim0, ylim1].
                  1. xlim0 : float
                             The minimum value of x-axis.
                  2. xlim1 : float
                             The maximum value of x-axis.
                  3. ylim0 : float
                             The minimum value of y-axis.
                  4. ylim1 : float
                             The maximum value of y-axis.
              
          Returns :
          -------
                  It return the .png file of masked aperture.
          """
          w = ap.aper(apt[0], apt[1], apt[2], apt[3], apt[4])
          plt.close()
          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [12,12]
          plt.rcParams["font.weight"] = "bold" 
          plt.rcParams["axes.labelweight"] = "bold"
          plt.contourf(apt[0], apt[1], w) 
          plt.title('Mask width = ' + str(apt[3]) + ' m & orientation = ' + str(np.round(apt[4]*180/np.pi, decimals=1)) + '$^\circ$', fontweight='bold')
          plt.xlabel("East x in (meter)")
          plt.ylabel("North y in (meter)")
          plt.gca().set_aspect('equal')
          plt.xlim(lim[0], lim[1])
          plt.ylim(lim[2], lim[3])
          plt.savefig(fname)
          plt.show()

      def gplot(self, fname, apt, uv, pa, XY, sf): 
          """
          Plot the signal for a masked aperture.
          
          Parameters :
          ----------
          fname : str
                  The output .png file.
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
          uv    : list
                  List of variational xy plane i.e. uv = [xcord, ycord].
                  1. xcord : N x N 
                             For the baseline along x-direction.
                  2. ycord : N x N 
                             For the baseline along y-direction.
          pa    : list
                  The list of parameters of binary star system i.e. [Tmp_a, Tmp_b, Tmp_r, dis_e, lam_c, lam_e].
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
          XY    : list
                  List of XY plane of orbit of binary i.e. [Xcord, Ycord].
                  1. Xcord : int
                             X-coordinate of orbit of binary on observational time of SII. 
                  2. Ycord : int
                             Y-coordinate of orbit of binary on observational time of SII.
          sf    : list
                  List of radius of objects and linear limb-darkening coefficient i.e. sf = [Rad_a, Rad_b, Rad_r, lin_c].
                  1. Rad_a : float
                             The radius of source A.
                  2. Rad_b : float
                             The radius of source B.
                  3. Rad_r : float
                             The maximum radius of Roche lobe from the surface of star.
                  4. lin_c : float
                             The linear limb-darkening coefficients.
          
          Returns :
          -------
                  Return a .png file showing HBT Signal for a Particular Masked Aperture.
          """
          md.Tmp_a = pa[0]
          md.Tmp_b = pa[1]
          md.Tmp_r = pa[2]
          md.dis_e = pa[3]
          md.lam_c = pa[4]
          md.lam_e = pa[5]
          w = ap.aper(apt[0], apt[1], apt[2], apt[3], apt[4])
          g = md.w_signal(uv[0], uv[1], XY[0], XY[1], sf[0], sf[1], w)
          plt.close()
          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [12,12]
          plt.rcParams["font.weight"] = "bold" 
          plt.rcParams["axes.labelweight"] = "bold"
          plt.contourf(uv[0], uv[1], abs(g))
          plt.title('HBT Signal for a Particular Masked Aperture', fontweight='bold')
          plt.xlabel("East x in (meter)")
          plt.ylabel("North y in (meter)")
          plt.savefig(fname)
          plt.show()

      def cmask(self, fname, apt, uv, bs, pa, sf, XY, jlday, lim):
          """
          fname : str
                  The output .png file.
          apt   : list
                  List of parameters of masked aperture i.e. ap = [aperX, aperY, radii, width, orien].
                  1. aperX : N x N 
                             For the aperture along x-direction.
                  2. aperY : N x N 
                             For the aperture along y-direction.
                  3. radii : int or float
                             Radius of telescopes.
                  4. width : float
                             The range of width of masking strips [width1, width2].
                  5. orien : float
                             The range of orientation of masking strips in radian [orien1, orien2].
          uv    : list
                  List of variational xy plane i.e. uv = [xcord, ycord].
                  1. xcord : N x N 
                             For the baseline along x-direction.
                  2. ycord : N x N 
                             For the baseline along y-direction.
          bs    : list
                  List of baseline on each observational time in xy direction i.e. bs = [baseX, baseY].
                  baseX : float
                          The rotation of baseline along x-axis.
                  baseY : float
                          The rotation of baseline along y-axis.
          pa    : list
                  The list of parameters of binary star system i.e. [Tmp_a, Tmp_b, Tmp_r, dis_e, lam_c, lam_e].
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
          sf    : list
                  List of radius of objects and linear limb-darkening coefficient i.e. sf = [Rad_a, Rad_b, Rad_r, lin_c].
                  1. Rad_a : float
                             The radius of source A.
                  2. Rad_b : float
                             The radius of source B.
                  3. Rad_r : float
                             The maximum radius of Roche lobe from the surface of star.
                  4. lin_c : float
                             The linear limb-darkening coefficients.
          XY    : list
                  List of XY plane of orbit of binary i.e. [Xcord, Ycord].
                  1. Xcord : int
                             X-coordinate of orbit of binary on observational time of SII. 
                  2. Ycord : int
                             Y-coordinate of orbit of binary on observational time of SII.
          jlday : array
                  The array of julian day of all observation. 
          lim   : list
                  List of limiting value for Y axis i.e. lim = [ylim0, ylim1].
                  1. ylim0 : float
                             The minimum value of y-axis.
                  2. ylim1 : float
                             The maximum value of y-axis.

          Returns :
          -------
                  Return a .png file showing HBT Signal for different Masked Aperture of a baseline.
          """
          md.Tmp_a = pa[0]
          md.Tmp_b = pa[1]
          md.Tmp_r = pa[2]
          md.dis_e = pa[3]
          md.lam_c = pa[4]
          md.lam_e = pa[5]
          plt.close()
          plt.rcParams.update({'font.size': 10})
          plt.rcParams["figure.figsize"] = [16,8]
          plt.rcParams["font.weight"] = "bold" 
          plt.rcParams["axes.labelweight"] = "bold"
          text_kwargs = dict(ha='center', va='center', fontsize=10)          
          for n in np.linspace(apt[4][0], apt[4][1], 51):
              pic = 0
              for m in np.linspace(apt[3][0], apt[3][1], 51):
                  w = ap.aper(apt[0], apt[1], apt[2], m, n)
                  g = md.w_signal(uv[0], uv[1], XY[0], XY[1], sf[0], sf[1], w)
    
                  plt.subplot(1, 3, 2)
                  plt.contourf(apt[0], apt[1], w)
                  plt.colorbar(shrink=0.50)
                  plt.gca().set_aspect('equal') 
                  plt.title('Ori. =%5.1f' % (np.round(n*180/np.pi, decimals=1)) + '$^\circ$,', loc = "left", fontweight='bold')
                  plt.title('Width =%5.1f mm' % (1e3*m), loc = "right", fontweight='bold')

                  plt.subplot(1, 3, 3)
                  plt.contourf(uv[0], uv[1], abs(g))
                  plt.colorbar(shrink=0.50)
                  plt.plot(bs[0], bs[1], '.', color='red')
                  plt.gca().set_aspect('equal')
                  plt.title('Track of baseline', **text_kwargs, fontweight='bold')  

                  plt.subplot(1, 3, 1)
                  xgr = uv[0][0, :]
                  ygr = uv[1][:, 0]
                  z = np.transpose(abs(g))
                  spline = RectBivariateSpline(xgr, ygr, z)
                  gt = spline.ev(bs[0], bs[1])
                  plt.plot(jlday, gt)
                  plt.xlabel('Julian Day', ha='center')
                  plt.ylabel('correlated Signal $g_w(u)$')
                  plt.title('Accordingly mask orientation and width', **text_kwargs, fontweight='bold')
                  plt.ylim(lim[0], lim[1])
                  plt.savefig(fname + str(pic) + '.png')
                  pic += 1
                  plt.clf()

      def wmask(self, fname, apt, uv, bs, pa, sf, XY, jlday, lim):
          """
          fname : str
                  The output .png file.
          apt   : list
                  List of parameters of masked aperture i.e. ap = [aperX, aperY, radii, width, orien].
                  1. aperX : N x N 
                             For the aperture along x-direction.
                  2. aperY : N x N 
                             For the aperture along y-direction.
                  3. radii : int or float
                             Radius of telescopes.
                  4. width : float
                             The range of width of masking strips [width1, width2].
                  5. orien : float
                             The particular orientation of masking strips in radian.
          uv    : list
                  List of variational xy plane i.e. uv = [xcord, ycord].
                  1. xcord : N x N 
                             For the baseline along x-direction.
                  2. ycord : N x N 
                             For the baseline along y-direction.
          bs    : list
                  List of baseline on each observational time in xy direction i.e. bs = [baseX, baseY].
                  baseX : float
                          The rotation of baseline along x-axis.
                  baseY : float
                          The rotation of baseline along y-axis.
          pa    : list
                  The list of parameters of binary star system i.e. [Tmp_a, Tmp_b, Tmp_r, dis_e, lam_c, lam_e].
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
          sf    : list
                  List of radius of objects and linear limb-darkening coefficient i.e. sf = [Rad_a, Rad_b, Rad_r, lin_c].
                  1. Rad_a : float
                             The radius of source A.
                  2. Rad_b : float
                             The radius of source B.
                  3. Rad_r : float
                             The maximum radius of Roche lobe from the surface of star.
                  4. lin_c : float
                             The linear limb-darkening coefficients.
          XY    : list
                  List of XY plane of orbit of binary i.e. [Xcord, Ycord].
                  1. Xcord : int
                             X-coordinate of orbit of binary on observational time of SII. 
                  2. Ycord : int
                             Y-coordinate of orbit of binary on observational time of SII.
          jlday : array
                  The array of julian day of all observation. 
          lim   : list
                  List of limiting value for Y axis i.e. lim = [ylim0, ylim1].
                  1. ylim0 : float
                             The minimum value of y-axis.
                  2. ylim1 : float
                             The maximum value of y-axis.
          
          Returns :
          -------
                  Return a .png file showing HBT Signal for a different Masked Aperture with fixed orientation of a baseline.
          """
          md.Tmp_a = pa[0]
          md.Tmp_b = pa[1]
          md.Tmp_r = pa[2]
          md.dis_e = pa[3]
          md.lam_c = pa[4]
          md.lam_e = pa[5]
          plt.close()
          plt.rcParams.update({'font.size': 10})
          plt.rcParams["figure.figsize"] = [16,8]
          plt.rcParams["font.weight"] = "bold" 
          plt.rcParams["axes.labelweight"] = "bold"
          text_kwargs = dict(ha='center', va='center', fontsize=10)
          pic = 0
          for m in np.linspace(apt[3][0], apt[3][1], 100):  
              w = ap.aper(apt[0], apt[1], apt[2], m, apt[4])
              g = md.w_signal(uv[0], uv[1], XY[0], XY[1], sf[0], sf[1], w)
    
              plt.subplot(1, 3, 2)
              plt.contourf(apt[0], apt[1], w)
              plt.colorbar(shrink=0.50)
              plt.gca().set_aspect('equal')
              plt.title('Mask width %5.2f mm' % (1e3*m), **text_kwargs, fontweight='bold')

              plt.subplot(1, 3, 3)
              plt.contourf(uv[0], uv[1], abs(g))
              plt.colorbar(shrink=0.50)
              plt.plot(bs[0], bs[1], '.', color='red')
              plt.gca().set_aspect('equal')
              plt.title('Track of baseline',**text_kwargs, fontweight='bold')  

              plt.subplot(1, 3, 1)
              xgr = uv[0][0, :]
              ygr = uv[1][:, 0]
              z = np.transpose(abs(g))
              spline = RectBivariateSpline(xgr, ygr, z)
              gt = spline.ev(bs[0], bs[1])
              plt.plot(jlday, gt)
              plt.xlabel('Julian Day', ha='center')
              plt.ylabel('correlated Signal $g_w(u)$')
              plt.title('According to mask width', **text_kwargs, fontweight = 'bold')
              plt.ylim(lim[0], lim[1])
              a = 10000*m
              s = '%3.2f' %a
              plt.savefig(fname + str(pic) + '.png')
              pic += 1
              plt.clf()

      def omask(self, fname, apt, uv, bs, pa, sf, XY, jlday, lim):
          """
          fname : str
                  The output .png file.
          apt   : list
                  List of parameters of masked aperture i.e. ap = [aperX, aperY, radii, width, orien].
                  1. aperX : N x N 
                             For the aperture along x-direction.
                  2. aperY : N x N 
                             For the aperture along y-direction.
                  3. radii : int or float
                             Radius of telescopes.
                  4. width : float
                             The particular width of masking strips.
                  5. orien : float
                             The range of orientation of masking strips in radian [orien1, orien2].
          uv    : list
                  List of variational xy plane i.e. uv = [xcord, ycord].
                  1. xcord : N x N 
                             For the baseline along x-direction.
                  2. ycord : N x N 
                             For the baseline along y-direction.
          bs    : list
                  List of baseline on each observational time in xy direction i.e. bs = [baseX, baseY].
                  baseX : float
                          The rotation of baseline along x-axis.
                  baseY : float
                          The rotation of baseline along y-axis.
          pa    : list
                  The list of parameters of binary star system i.e. [Tmp_a, Tmp_b, Tmp_r, dis_e, lam_c, lam_e].
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
          sf    : list
                  List of radius of objects and linear limb-darkening coefficient i.e. sf = [Rad_a, Rad_b, Rad_r, lin_c].
                  1. Rad_a : float
                             The radius of source A.
                  2. Rad_b : float
                             The radius of source B.
                  3. Rad_r : float
                             The maximum radius of Roche lobe from the surface of star.
                  4. lin_c : float
                             The linear limb-darkening coefficients.
          XY    : list
                  List of XY plane of orbit of binary i.e. [Xcord, Ycord].
                  1. Xcord : int
                             X-coordinate of orbit of binary on observational time of SII. 
                  2. Ycord : int
                             Y-coordinate of orbit of binary on observational time of SII.
          jlday : array
                  The array of julian day of all observation. 
          lim   : list
                  List of limiting value for Y axis i.e. lim = [ylim0, ylim1].
                  1. ylim0 : float
                             The minimum value of y-axis.
                  2. ylim1 : float
                             The maximum value of y-axis.
          
          Returns :
          -------
                  Return a .png file showing HBT Signal for a different Masked Aperture with fixed width of a baseline.
          """
          md.Tmp_a = pa[0]
          md.Tmp_b = pa[1]
          md.Tmp_r = pa[2]
          md.dis_e = pa[3]
          md.lam_c = pa[4]
          md.lam_e = pa[5]
          plt.close()
          plt.rcParams.update({'font.size': 10})
          plt.rcParams["figure.figsize"] = [16,8]
          plt.rcParams["font.weight"] = "bold" 
          plt.rcParams["axes.labelweight"] = "bold"
          text_kwargs = dict(ha='center', va='center', fontsize=10)
          pic = 0
          for n in np.linspace(apt[4][0], apt[4][1], 100):  
              w = ap.aper(apt[0], apt[1], apt[2], apt[3], n)
              g = md.w_signal(uv[0], uv[1], XY[0], XY[1], sf[0], sf[1], w)
    
              plt.subplot(1, 3, 2)
              plt.contourf(apt[0], apt[1], w)
              plt.colorbar(shrink=0.50)
              plt.gca().set_aspect('equal')
              plt.title('Mask orientation %5.2f' % (np.round(n*180/np.pi, decimals=1))+'$^\circ$', **text_kwargs, fontweight='bold')

              plt.subplot(1, 3, 3)
              plt.contourf(uv[0], uv[1], abs(g))
              plt.colorbar(shrink=0.50)
              plt.plot(bs[0], bs[1], '.', color='red')
              plt.gca().set_aspect('equal')
              plt.title('Track of baseline', **text_kwargs, fontweight='bold')  

              plt.subplot(1, 3, 1)
              xgr = uv[0][0, :]
              ygr = uv[1][:, 0]
              z = np.transpose(abs(g))
              spline = RectBivariateSpline(xgr, ygr, z)
              gt = spline.ev(bs[0], bs[1])
              plt.plot(jlday, gt)
              plt.xlabel('Julian Day', ha='center')
              plt.ylabel('correlated Signal $g_w(u)$')
              plt.title('According to mask orientation', **text_kwargs, fontweight = 'bold')
              plt.ylim(lim[0], lim[1])
              a = n*1000
              s = '%3.2f' %a 
              plt.savefig(fname + str(pic) + '.png')
              pic += 1
              plt.clf()

      def var_w(self, fname, apt, uv, bs, pa, sf, XY, jlday, lim, color):
          """
          fname : str
                  The output .png file.
          apt   : list
                  List of parameters of masked aperture i.e. ap = [aperX, aperY, radii, width, orien].
                  1. aperX : N x N 
                             For the aperture along x-direction.
                  2. aperY : N x N 
                             For the aperture along y-direction.
                  3. radii : int or float
                             Radius of telescopes.
                  4. width : float
                             The list of width of masking strips.
                  5. orien : float
                             The particular orientation of masking strips in radian.
          uv    : list
                  List of variational xy plane i.e. uv = [xcord, ycord].
                  1. xcord : N x N 
                             For the baseline along x-direction.
                  2. ycord : N x N 
                             For the baseline along y-direction.
          bs    : list
                  List of baseline on each observational time in xy direction i.e. bs = [baseX, baseY].
                  baseX : float
                          The rotation of baseline along x-axis.
                  baseY : float
                          The rotation of baseline along y-axis.
          pa    : list
                  The list of parameters of binary star system i.e. [Tmp_a, Tmp_b, Tmp_r, dis_e, lam_c, lam_e].
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
          sf    : list
                  List of radius of objects and linear limb-darkening coefficient i.e. sf = [Rad_a, Rad_b, Rad_r, lin_c].
                  1. Rad_a : float
                             The radius of source A.
                  2. Rad_b : float
                             The radius of source B.
                  3. Rad_r : float
                             The maximum radius of Roche lobe from the surface of star.
                  4. lin_c : float
                             The linear limb-darkening coefficients.
          XY    : list
                  List of XY plane of orbit of binary i.e. [Xcord, Ycord].
                  1. Xcord : int
                             X-coordinate of orbit of binary on observational time of SII. 
                  2. Ycord : int
                             Y-coordinate of orbit of binary on observational time of SII.
          jlday : array
                  The array of julian day of all observation. 
          lim   : list
                  List of limiting value for Y axis i.e. lim = [ylim0, ylim1].
                  1. ylim0 : float
                             The minimum value of y-axis.
                  2. ylim1 : float
                             The maximum value of y-axis.
          color : list
                  The list of color for each baseline tracking.
          
          Returns :
          -------
                  Return a .png file showing HBT Signal for a different Masked Aperture's width and fixed orientation for a baseline.
          """
          md.Tmp_a = pa[0]
          md.Tmp_b = pa[1]
          md.Tmp_r = pa[2]
          md.dis_e = pa[3]
          md.lam_c = pa[4]
          md.lam_e = pa[5]
          plt.close()
          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [12,12]
          plt.rcParams["font.weight"] = "bold" 
          plt.rcParams["axes.labelweight"] = "bold"
          for i in range(len(apt[3])):
              w = ap.aper(apt[0], apt[1], apt[2], apt[3][i], apt[4])
              g = md.w_signal(uv[0], uv[1], XY[0], XY[1], sf[0], sf[1], w)
              xgr = uv[0][0, :]
              ygr = uv[1][:, 0]
              z = np.transpose(abs(g))
              spline = RectBivariateSpline(xgr, ygr, z)
              gt = spline.ev(bs[0], bs[1])
              plt.plot(jlday, gt, label=str(np.round(apt[3][i]*1000, decimals=1)) + 'mm', color = color[i])
          plt.xlabel('Julian Day', ha='center')
          plt.ylabel('Correlated signal $g_w(u)$')
          plt.title('According to variation of mask width', fontweight = 'bold')
          plt.ylim(lim[0], lim[1])
          plt.legend(loc="lower right")
          plt.gcf().set_size_inches(10, 10)
          plt.savefig(fname + '.png')
          plt.show()

      def var_o(self, fname, apt, uv, bs, pa, sf, XY, jlday, lim, color):
          """
          fname : str
                  The output .png file.
          apt   : list
                  List of parameters of masked aperture i.e. ap = [aperX, aperY, radii, width, orien].
                  1. aperX : N x N 
                             For the aperture along x-direction.
                  2. aperY : N x N 
                             For the aperture along y-direction.
                  3. radii : int or float
                             Radius of telescopes.
                  4. width : float
                             The list of width of masking strips.
                  5. orien : float
                             The particular orientation of masking strips in radian.
          uv    : list
                  List of variational xy plane i.e. uv = [xcord, ycord].
                  1. xcord : N x N 
                             For the baseline along x-direction.
                  2. ycord : N x N 
                             For the baseline along y-direction.
          bs    : list
                  List of baseline on each observational time in xy direction i.e. bs = [baseX, baseY].
                  baseX : float
                          The rotation of baseline along x-axis.
                  baseY : float
                          The rotation of baseline along y-axis.
          pa    : list
                  The list of parameters of binary star system i.e. [Tmp_a, Tmp_b, Tmp_r, dis_e, lam_c, lam_e].
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
          sf    : list
                  List of radius of objects and linear limb-darkening coefficient i.e. sf = [Rad_a, Rad_b, Rad_r, lin_c].
                  1. Rad_a : float
                             The radius of source A.
                  2. Rad_b : float
                             The radius of source B.
                  3. Rad_r : float
                             The maximum radius of Roche lobe from the surface of star.
                  4. lin_c : float
                             The linear limb-darkening coefficients.
          XY    : list
                  List of XY plane of orbit of binary i.e. [Xcord, Ycord].
                  1. Xcord : int
                             X-coordinate of orbit of binary on observational time of SII. 
                  2. Ycord : int
                             Y-coordinate of orbit of binary on observational time of SII.
          jlday : array
                  The array of julian day of all observation. 
          lim   : list
                  List of limiting value for Y axis i.e. lim = [ylim0, ylim1].
                  1. ylim0 : float
                             The minimum value of y-axis.
                  2. ylim1 : float
                             The maximum value of y-axis.
          color : list
                  The list of color for each baseline tracking.
          
          Returns :
          -------
          Return a .png file showing HBT Signal for a different Masked Aperture's orientation and fixed width for a baseline.
          """
          md.Tmp_a = pa[0]
          md.Tmp_b = pa[1]
          md.Tmp_r = pa[2]
          md.dis_e = pa[3]
          md.lam_c = pa[4]
          md.lam_e = pa[5]
          plt.close()
          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [12,12]
          plt.rcParams["font.weight"] = "bold" 
          plt.rcParams["axes.labelweight"] = "bold"
          for i in range(len(apt[4])):
              w = ap.aper(apt[0], apt[1], apt[2], apt[3], apt[4][i])
              g = md.w_signal(uv[0], uv[1], XY[0], XY[1], sf[0], sf[1], w)
              xgr = uv[0][0, :]
              ygr = uv[1][:, 0]
              z = np.transpose(abs(g))
              spline = RectBivariateSpline(xgr, ygr, z)
              gt = spline.ev(bs[0], bs[1])
              plt.plot(jlday, gt, label=str(np.round(apt[4][i]*180/np.pi, decimals=1))+'$^\circ$', color = color[i])
          plt.xlabel('Julian Day', ha='center')
          plt.ylabel('Correlated signal $g_w(u)$')
          plt.title('According to variation of mask orientation', fontweight = 'bold')
          plt.ylim(lim[0], lim[1])
          plt.legend(loc="lower right")
          plt.gcf().set_size_inches(10, 10)
          plt.savefig(fname + '.png')
          plt.show() 

