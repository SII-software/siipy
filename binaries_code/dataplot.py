from dynesty.utils import resample_equal
import matplotlib.pyplot as plt
from matplotlib import cm
from glob import glob
import numpy as np
import pickle
import corner
import math
import re

from model import modeling
from obspoint import vary_base
vb = vary_base()
cod = modeling()

class dplot():
      """
      Creates an object that plot the simulated, observed, binned or model data.
      It also plot all the interferometric data.
      """
      def cntrplot(self, pa, sf, XY, bs, step, nbase, color, fname):
          """
          The track of baselines on contour plot of Visibility signal of close or Regor type binary.
          This tracking is on each observational interval.
          
          Parameters :
          ----------
          
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
          bs    : list
                  List of baseline on each observational time in xy direction i.e. bs = [baseX, baseY].
                  baseX : float
                          The rotation of baseline along x-axis.
                  baseY : float
                          The rotation of baseline along y-axis.
          step  : array
                  The number of steps in each observational interval. 
          nbase : list
                  List of tupples, which are the name of baseline.
          color : list
                  The list of color for each baseline tracking.
          fname : str
                  The output .png file.
                  
          Returns :
          -------
                  Return the contour plot with tracking of baselines.
          """
          plt.close()
          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [14,12]
          plt.rcParams["font.weight"] = "bold"
          plt.rcParams["axes.labelweight"] = "bold"
          cod.Tmp_a = pa[0]
          cod.Tmp_b = pa[1]
          cod.Tmp_r = pa[2]
          cod.dis_e = pa[3]
          cod.lam_c = pa[4]
          cod.lam_e = pa[5]
          st = np.cumsum(step)
          x = np.linspace(bs[0].min(), bs[0].max(), 1000)
          y = np.linspace(bs[1].min(), bs[1].max(), 1000)
          gX, gY = np.meshgrid(x, y)
          for k in range(len(bs[0])):
              for i in range(len(step)):
                  pic = 0
                  for j in range(step[i]):
                      if i==0:
                         g = cod.hbt(gX, gY, XY[0][j], XY[1][j], sf[0], sf[1], sf[2], sf[3])
                         plt.plot(bs[0][k,:j+1], bs[1][k,:j+1], '.', markersize=12, c = color[k], label = nbase[k])
                      else:
                         g = cod.hbt(gX, gY, XY[0][j+st[i-1]], XY[1][j+st[i-1]], sf[0], sf[1], sf[2], sf[3])
                         plt.plot(bs[0][k,st[i-1]:j+1+st[i-1]], bs[1][k,st[i-1]:j+1+st[i-1]], '.', markersize=12, c = color[k], label = nbase[k])
                      plt.contourf(gX, gY, abs(g), cmap="Greens", vmin=0, vmax=1.1)
                      plt.title('On the observational interval ' + str(i + 1), fontweight='bold')       
                      plt.colorbar()               
                      plt.legend(loc="upper right")
                      plt.savefig(fname + str(pic) + '.png')
                      pic += 1
                      plt.clf() 

      def sigplot(self, sdata, step, nbase, jlday, color, fname):
          """
          The observaed signal in each observational interval for each baselines.
          
          Parameters :
          ----------
          
          sdata : str
                  The input .npy file containing HBT data for variational baseline.
          step  : array
                  The number of steps in each observational interval. 
          nbase : list
                  List of tupples, which are the name of baseline.
          jlday : array
                  The array of julian day of all observation.
          color : list
                  The list of color for each baseline tracking.
          fname : str
                  The output .png file.
                  
          Returns :
          -------
                  Return the .png file which shows the observaed signal in each observational interval for each baselines.
          """
          plt.close()
          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [12,12]
          plt.rcParams["font.weight"] = "bold"
          plt.rcParams["axes.labelweight"] = "bold"
          gtn = np.load(sdata)
          st = np.cumsum(step)
          for k in range(len(nbase)):
              for i in range(len(step)):
                  pic = 0
                  for j in range(step[i]):
                      if i==0:
                         plt.errorbar(jlday[:j+1], gtn[k,:j+1], yerr=np.std(gtn), color=color[k], fmt='o', markersize=8, label = nbase[k])
                         plt.xlim(jlday[0], jlday[st[i]-1])
                      else:
                         plt.errorbar(jlday[st[i-1]:j+1+st[i-1]], gtn[k,st[i-1]:j+1+st[i-1]], yerr=np.std(gtn), color=color[k], fmt='o', markersize=8, label = nbase[k])
                         plt.xlim(jlday[st[i-1]], jlday[st[i]-1])
                      plt.ylim(gtn.min()-np.std(gtn), gtn.max()+np.std(gtn))
                      plt.xlabel('Julian Day', ha='center')
                      plt.ylabel('correlated Signal $g_w(u)$ with noise')
                      plt.title('On the observational interval ' + str(i+1), fontweight='bold')
                      plt.legend(loc="upper right")
                      plt.savefig(fname + str(pic) + '.png')
                      pic += 1
                      plt.clf()

      def cmplot(self, sdata, pa, sf, XY, bs, step, nbase, jlday, color, fname):
          """
          The track of baselines on contour plot of Visibility signal of close or Regor type binary.
          This tracking is on each observational interval.
          
          Parameters :
          ----------
          
          sdata : str
                  The input .npy file containing HBT data for variational baseline.
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
          bs    : list
                  List of baseline on each observational time in xy direction i.e. bs = [baseX, baseY].
                  baseX : float
                          The rotation of baseline along x-axis.
                  baseY : float
                          The rotation of baseline along y-axis.
          step  : array
                  The number of steps in each observational interval. 
          nbase : list
                  List of tupples, which are the name of baseline.
          jlday : array
                  The array of julian day of all observation.
          color : list
                  The list of color for each baseline tracking.
          fname : str
                  The output .png file.
                  
          Returns :
          -------
                  Return the contour plot with tracking of baselines and the corresponding captured signal.
          """
          plt.close()
          plt.rcParams.update({'font.size': 10})
          plt.rcParams["figure.figsize"] = [18,8]
          text_kwargs = dict(ha='center', va='center', fontsize=12)
          text_kwarg = dict(ha='center', va='center', fontsize=18)
          plt.rcParams["font.weight"] = "bold"
          plt.rcParams["axes.labelweight"] = "bold"
          cod.Tmp_a = pa[0]
          cod.Tmp_b = pa[1]
          cod.Tmp_r = pa[2]
          cod.dis_e = pa[3]
          cod.lam_c = pa[4]
          cod.lam_e = pa[5]
          gtn = np.load(sdata)
          st = np.cumsum(step)
          x = np.linspace(bs[0].min(), bs[0].max(), 1000)
          y = np.linspace(bs[1].min(), bs[1].max(), 1000)
          gX, gY = np.meshgrid(x, y)
          for k in range(len(bs[0])):
              for i in range(len(step)):
                  pic = 0
                  for j in range(step[i]):     
                      plt.subplot(1,2,1)
                      if i==0:
                         g = cod.hbt(gX, gY, XY[0][j], XY[1][j], sf[0], sf[1], sf[2], sf[3])
                         plt.plot(bs[0][k,:j+1], bs[1][k,:j+1], '.', markersize=12, c = color[k], label = nbase[k])
                      else:
                         g = cod.hbt(gX, gY, XY[0][j+st[i-1]], XY[1][j+st[i-1]], sf[0], sf[1], sf[2], sf[3])
                         plt.plot(bs[0][k,st[i-1]:j+1+st[i-1]], bs[1][k,st[i-1]:j+1+st[i-1]], '.', markersize=12, c = color[k], label = nbase[k])
                      plt.contourf(gX, gY, abs(g), cmap="Greens", vmin=0, vmax=1.1)
                      plt.colorbar()  
                      plt.title('Track of baseline ' + str(nbase[k]) + ' on 2-d SII correlation', **text_kwargs, fontweight='bold')
                      plt.xlabel('Baseline in east (meter)', fontweight='bold')
                      plt.ylabel('Baseline in north (meter)', fontweight='bold')
                      
                      plt.subplot(1,2,2)
                      if i==0:
                         plt.errorbar(jlday[:j+1], gtn[k,:j+1], yerr=np.std(gtn), color=color[k], fmt='o', markersize=8, label = nbase[k])
                         plt.xlim(jlday[0], jlday[st[i]-1])
                      else:
                         plt.errorbar(jlday[st[i-1]:j+1+st[i-1]], gtn[k,st[i-1]:j+1+st[i-1]], yerr=np.std(gtn), color=color[k], fmt='o', markersize=8, label = nbase[k])
                         plt.xlim(jlday[st[i-1]], jlday[st[i]-1])
                      plt.ylim(gtn.min()-np.std(gtn), gtn.max()+np.std(gtn))
                      plt.title('The correlated signal $g_w(u)$ with noise for baseline ' + str(nbase[k]), **text_kwargs, fontweight='bold')
                      plt.xlabel('Julian Day', ha='center', fontweight='bold')
                      plt.suptitle('On the observational interval ' + str(i+1), color='darkred', **text_kwarg, fontweight='bold')
                      plt.savefig(fname + str(pic) + '.png')
                      pic += 1
                      plt.clf()
                      

      def allsigplot(self, sdata, fname, jlday, nbase, color):
          """
          The HBT data of all baseline together for one observational interval.
          
          Parameters :
          ----------
          
          sdata : str
                  The input .npy file containing HBT data for variational baseline.
          fname : str
                  The output .png file. 
          jlday : array
                  The array of julian day of all observation.
          nbase : list
                  List of tupples, which are the name of baseline.
          color : list
                  The list of color for each baseline tracking.
          
          Returns :
          -------
                  Return .png file showing the observed signal from each baseline.
          """
          data = np.load(sdata)
          plt.close()
          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [12,12]
          plt.rcParams["font.weight"] = "bold"
          plt.rcParams["axes.labelweight"] = "bold"
          fig, ax = plt.subplots(len(nbase), sharex=True)
          for i in range(len(nbase)):
              ax[i].plot(jlday, data[i], color = color[i], label = nbase[i])
              ax[i].legend(bbox_to_anchor=(0.90, 0.90), loc=3)
          fig.text(0.5, 0.04, 'Julian Observation Day', ha='center')
          fig.text(0.01, 0.5, 'Correlated Signal $g_w(u)$ According to baseline', va='center', rotation='vertical')
          fig.suptitle("Signal on observational interval", y=0.91, fontweight='bold')
          plt.subplots_adjust(left=0.1, bottom=None, right=None, top=None, wspace=None, hspace=None)
          plt.savefig(fname + ".png")
          plt.show()

      def spflux(self, fname, Tmper, Radii, dis_e, lam, del_lam):
          """
          Plot the Spectral photon flux for black body sources.
          
          Parameters :
          ----------
          fname : str
                  The output .png file.
          Tmper : int
                  List of teperature of objects.
          Radii : int
                  List of radius of objects.
          dis_e : float
                  The distance of star system from earth.
          freqn : float
                  List of observational frequencies of objects.
          del_freq: float
                  The value of bandwidth
            
          Returns :
          -------
                  Return the .png file showing the Spectral photon flux for black body sources.             
          
          """
          plt.close()
          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [12,12]
          plt.rcParams["font.weight"] = "bold"
          plt.rcParams["axes.labelweight"] = "bold"
          # defined constant
          c = 299792458.0       # in m/s
          h = 6.62607015e-34    # J.s
          k = 1.380649e-23      # J/K
          # Total flux
          for i in range(len(Radii)):
              sp = []
              for j in range(len(lam)):
                  a1 = np.pi * Radii[i]**2
                  a1 /= (lam[j] * dis_e)**2
                  a2 = h*c/(lam[j] * k * (Tmper[i] + 2.725))
                  a2 = np.clip(a2, 0, 1000)
                  sp.append(a1 / (np.exp(a2) - 1))
              fl = np.array(sp)
              plt.plot(lam*1e9, fl, linestyle='dashed', label = 'Object' + str(i+1))
          plt.legend(loc="upper right")
          plt.xlabel('$\\lambda$ (nm)')
          plt.ylabel('$\\Phi$ ($\\rm\\,photons\\;m^{-2}\;s^{-1}\;Hz^{-1}$)')        
          plt.tight_layout()
          plt.savefig(fname + '.png')
          plt.show()
          
          # Photon rate
          for i in range(len(Radii)):
              sp = []
              for j in range(len(lam)):
                  a1 = np.pi * Radii[i]**2
                  a1 /= (lam[j] * dis_e)**2
                  a2 = h*c/(lam[j] * k * (Tmper[i] + 2.725))
                  a2 = np.clip(a2, 0, 1000)
                  phi = a1 / (np.exp(a2) - 1)
                  sp.append(c*del_lam*phi/lam[j]**2)
              fl = np.array(sp)
              plt.plot(lam*1e9, fl, linestyle='dashed', label = 'Object' + str(i+1))
          plt.legend(loc="upper right")
          plt.xlabel('$\\lambda$ (nm)')
          plt.ylabel('$\\Phi \, \\Delta \\nu$ ($\\rm\\,photons\\;m^{-2}\;s^{-1}$)')
          plt.tight_layout()
          plt.savefig(fname + 'rate' + '.png')
          plt.show()

      def bindt(self, jlday, sdata, mdata, fname, Nbinn, int_b):
          """
          Plot the observaed data, model data and binned data in interval n.
          
          Parameters :
          ----------
          jlday : array
                  The array of julian day of all observation.
          sdata : str
                  The input .npy file containing observed HBT data for variational baseline.
          mdata : str
                  The input .npy file containing the model HBT data for variational baseline.
          fname : str
                  The output .png file. 
          Nbinn : int
                  The number of times to bin.
          int_b : int
                  The binning interval of time.
            
          Returns :
          -------
                  Return the .png file showing the observaed data, model data and binned data in interval n.
          """
          plt.close()
          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [12,12]
          plt.rcParams["font.weight"] = "bold"
          plt.rcParams["axes.labelweight"] = "bold"
          bj =[]
          data = []
          ebar = []
          for i in range(Nbinn):
              bj.append(jlday[i*int_b])
              bn = sdata[i*int_b : (i+1)*int_b]
              data.append(np.mean(bn))
              ebar.append(np.std(bn)/len(bn)**.5)
          jday = np.array(bj)
          bdata = np.array(data)
          sebar = np.array(ebar)
          plt.plot(jlday, sdata, '.', label='Observed data')
          plt.plot(jlday, mdata, label='Noiseless signal')
          plt.errorbar(jday, bdata, yerr=sebar, fmt='o', label='Binned data')
          plt.xlabel('Julian Day', ha='center')
          plt.ylabel('Correlated signal $g_w(u)$ with noise')
          plt.title('Observed signal for a baseline', fontweight='bold')
          plt.legend(loc="upper right")
          plt.savefig(fname + '.png')
          plt.show()

      def error_plot(self, nights, fnames, fname1, fname2):
          """
          To see the drop in interferometric error with number of observational interval.
          
          Parameters :
          ----------
          
          nights : int
                   The number of observational interval.
          fnames : str
                   The input files in .pkl format for different observational interval.
          fname1 : str
                   The output .png file for X-coordinate error.
          fname2 : str
                   The output .png file for Y-coordinate error.
                
          Returns :
          -------
                   The .png files showing the drop in X and Y coordinate interferometric error.
          """
          filnames = glob(fnames)
          filnames.sort(key=lambda x:[int(c) if c.isdigit() else c for c in re.split(r'(\d+)',x)])
          sample = []
          for f in filnames:
                infile = open(f,'rb')
                results = pickle.load(infile)
                infile.close()
                weights = np.exp(results['logwt'] - results['logz'][-1])
                sample.append(resample_equal(results.samples, weights))
          postsamples = np.array(sample, dtype=object)
          # error in x And y using postsamples
          err_x = []
          err_y = []
          for i in range(len(postsamples)):
              samplei = np.copy(postsamples[i])
              p_x = np.percentile(samplei[:,0],[16,50,84])
              p_y = np.percentile(samplei[:,1],[16,50,84])
              err_x.append(p_x[2] - p_x[1])
              err_y.append(p_x[2] - p_x[1])
          error_x = err_x
          error_y = err_y

          plt.close()
          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [12,12]
          plt.rcParams["font.weight"] = "bold"
          plt.rcParams["axes.labelweight"] = "bold"
          plt.plot(nights, error_x, '.', color = 'blue', markersize=15)
          plt.title('Drop in interferometric error with observational interval length', fontweight='bold')
          plt.xlabel('Number Of Night')
          plt.ylabel('Interferometric 1-sigma error in X (ls)')
          plt.savefig(fname1 + '.png')
          plt.show()

          plt.plot(nights, error_y, '.', color = 'green', markersize=15)
          plt.title('Drop in interferometric error with observational length', fontweight='bold')
          plt.xlabel('Number Of Night')
          plt.ylabel('Interferometric 1-sigma error in Y (ls)')
          plt.savefig(fname2 + '.png')
          plt.show()

      def para_plot(self, truth, pname, dfile, cfile):
          """
          Plot the corner plot for each parameters.
          
          Parameters:
          ----------
          
          truth : list
                  List of estimated parameters value to generate simulated noisy HBT data.
          pname : tupple
                  Tupple of name of estimated parameters.
          dfile : str
                  The simulated file of parameters in .pkl format.
          cfile : str
                  The output corner .png file.
                
          Returns :
          -------
                  Return the corner plot of each parameters in .png file.
          """
          ndim = len(pname)
          infile = open(dfile, 'rb')
          results = pickle.load(infile)
          infile.close()
          weights = np.exp(results['logwt'] - results['logz'][-1])    # weighting the each nested 
          postsamples = resample_equal(results.samples, weights)      # sample for posterior distribution after weighting the each nested    
          # print the percentile value      
          med = np.zeros(ndim)
          for i in range(ndim):
              p = np.percentile(postsamples[:, i], [16, 50, 84])      # the percentile parameters value 16%, 50%, 84% 
              med[i] = p[1]
              print('%s in range %7.16e %7.16e %7.16e' % (pname[i],p[0],p[1],p[2]))
          # plotting with posterior samples
          plt.close()
          plt.rcParams.update({'font.size': 7})
          plt.rcParams["figure.figsize"] = [11,11]
          plt.rcParams["font.weight"] = "bold"
          plt.rcParams["axes.labelweight"] = "bold" 
          tex = dict(ha='center', fontsize=10, fontweight='bold')
          text = dict(ha='center', fontsize=10, fontweight='bold')
          fig = corner.corner(postsamples, labels= pname, color='g', truths=truth, truth_color='r', plot_datapoints=False, label_kwargs=tex, title_kwargs= text, show_titles=True, title_fmt=".2e")
          plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.08, hspace=0.15)
          plt.savefig(cfile + ".png") 
          plt.show()

