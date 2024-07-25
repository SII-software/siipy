import numpy as np
import matplotlib.pyplot as plt


class orbplot():
      """
      Create an object to analyze the orbit's integration result.
      """
      def pos_plot(self, xdata, ydata, Nobj, pname, rname):
          """ 
          Plot the position of N objects of a star system.

          Plot the respective position of Nth objects (depends on Nobj) w.r.t. first object's position.

          Parameters :
          ----------

          xdata : str
                  The name of .npy file containing X-coordinate of objects.
          ydata : str
                  The name of .npy file containing Y-coordinate of objects. 
          Nobj  : int
                  The number of objects to plot. Also the respective position.
                  If Nobj = 3, the respective position of 3rd objects to the 1st object.
          pname : str
                  The name of output file in .png format for position of objects.
          rname : str
                  The name of output file in .png format for respective position of two objects.

          Returns :
          -------
                  Return the two-dimensional position of N objects in .png format.
                  Also, return the respective position of two objects in .png format.
          """

          plt.close()
          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [12,12]
          plt.rcParams["font.weight"] = "bold"
          plt.rcParams["axes.labelweight"] = "bold"
          x = np.load(xdata)
          y = np.load(ydata)
          # position of stars
          for i in range(Nobj):
              plt.plot(x[i,:], y[i,:], linestyle= 'dashdot', label='Object ' + str(i+1))
          plt.xlabel('X in light second')
          plt.ylabel('Y in light second')
          plt.title("Position of " + str(Nobj) + " objects in star system", fontweight='bold')
          plt.legend(loc="upper right")
          plt.savefig(pname + '.png')
          plt.show()
          # respective orbit of two body system from origin body 
          if Nobj==2:        
             plt.plot(x[0,:]-x[1,:], y[0,:]-y[1,:], color='cyan', linestyle= 'dashdot')
          elif Nobj==3:
             plt.plot(x[0,:]-x[2,:], y[0,:]-y[2,:], color='cyan', linestyle= 'dashdot')
          else:
             plt.plot(x[0,:]-x[3,:], y[0,:]-y[3,:], color='cyan', linestyle= 'dashdot')
          plt.xlabel('X in light second')
          plt.ylabel('Y in light second')
          plt.title("Respective position of object " + str(Nobj) + " with respect to the first object", fontweight='bold')
          plt.savefig(rname + '.png')
          plt.show() 
          
      def vel_plot(self, tname, vname, Nobj, steps, fname):
          """
          Plot the velocity of objects in star system.

          Parameters :
          ----------

          tname : str
                  The name of .npy file containing time coordinate of objects.
          vname : str
                  The name of .npy file containing Z-coordinate of velocity of objects.
          Nobj  : int
                  The number of objects to be ploted.
          steps : int
                  The number of steps to be ploted.
          fname : str
                  The name of output file in .png format.

          Returns :
          -------
                  Return the .png files of velocity distribution of objects
          """
          plt.close()
          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [12,12]
          plt.rcParams["font.weight"] = "bold"
          plt.rcParams["axes.labelweight"] = "bold"
          t = np.load(tname)
          vz = np.load(vname)
          for i in range(Nobj):
              plt.plot(t[0 : steps], vz[i, : steps], label='Object ' + str(i+1))
          plt.xlabel('Time in seconds')
          plt.ylabel('Velocity in light second')
          plt.title('Compare the velocities of objects in star system', fontweight='bold')
          plt.legend(loc="upper right")
          plt.savefig(fname + '.png')
          plt.show()             

      def c_plot(self, xdata, ydata, dstep, steps, Nn, fname):
          """
          Plot the position of close binary stars on observational time.
          Alternatively, plot the respective orbit of close binary stars on observational time.

          Parameters :
          ----------

          xdata : str
                  The name of .npy file containing X-coordinate of objects.
          ydata : str
                  The name of .npy file containing Y-coordinate of objects.
          dstep : array
                  The first step on each observational interval.
          steps : array
                  The number of steps in each observational interval.
          Nn    : int
                  If n = 0, the position of each objects.
                  If n = 1, the respective position of two objects.
          fname : str
                  The name of output file in .png format. 

          Returns :
          -------
                  Return the .png files for position of two objects on observation time.
                  Or the .png files for respective orbit of two objects on observation time.
          """
          plt.close()
          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [12,12]
          plt.rcParams["font.weight"] = "bold"
          plt.rcParams["axes.labelweight"] = "bold"
          x = np.load(xdata)
          y = np.load(ydata)
          for j in range(0, len(dstep), 1):
              pic = 0
              for i in range(0, steps[j], 1):
                  if Nn==0:
                     plt.plot(x[0, :], y[0, :], color = 'red', linestyle= 'dotted', linewidth = 2)
                     plt.plot(x[1, :], y[1, :], color = 'green', linestyle= 'dotted', linewidth = 2)
                     plt.plot(x[0, dstep[j] : dstep[j]+i], y[0, dstep[j] : dstep[j]+i], color = 'black', linewidth = 3)
                     plt.plot(x[1, dstep[j] : dstep[j]+i], y[1, dstep[j] : dstep[j]+i], color = 'blue', linewidth = 3)
                     plt.plot(x[0, dstep[j]+i], y[0, dstep[j]+i], marker='*', color='red', markersize=15, label = 'Object 1')
                     plt.plot(x[1, dstep[j]+i], y[1, dstep[j]+i], marker='*', color='green', markersize=15, label = 'Object 2')
                     plt.xlabel('X in light second')
                     plt.ylabel('Y in light second')
                     plt.title("Objects 1 and 2 on " + 'observational interval ' + str(j+1), fontweight='bold')
                     plt.legend(loc="upper right")
                     plt.savefig(fname + str(pic) + '.png')
                     pic += 1
                     plt.clf()
                  elif Nn==1:
                     plt.plot(x[0, :] - x[1, :], y[0, :] - y[1, :], color = 'cyan', linestyle= 'dotted', linewidth = 2)
                     plt.plot(x[0, dstep[j] : dstep[j]+i] - x[1, dstep[j] : dstep[j]+i], y[0, dstep[j] : dstep[j]+i] - y[1, dstep[j] : dstep[j]+i], color = 'black', linewidth = 3)
                     plt.plot(x[0, dstep[j]+i] - x[1, dstep[j]+i], y[0, dstep[j]+i] - y[1, dstep[j]+i], marker='*', color='red', markersize=15, label = 'Orbit of Objects')
                     plt.xlabel('X in light second')
                     plt.ylabel('Y in light second')
                     plt.title("Orbit of objects on " + 'observational interval ' + str(j+1), fontweight='bold')
                     plt.legend(loc="upper right")
                     plt.savefig(fname + str(pic) + '.png')
                     pic += 1
                     plt.clf()

      def w_plot(self, xdata, ydata, Nobj, steps, slices, fname):
          """
          Plot the variable position of wide binary with a planet around first objects.

          Parameters :
          ----------

          xdata : str
                  The name of .npy file containing X-coordinate of objects.
          ydata : str
                  The name of .npy file containing Y-coordinate of objects.
          Nobj  : int
                  The object to be plotted. Nobj = 1 for first star only.
          steps : int
                  The number of steps of observation to be ploted.
          slices: int
                  The steps to be ignored between two observation.
          fname : str
                  The name of output file in .png format.

          Returns :
          -------
                  Return the .png files of position of a wide binary star system.   
          """
          plt.close()
          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [12,12]
          plt.rcParams["font.weight"] = "bold"
          plt.rcParams["axes.labelweight"] = "bold"
          x = np.load(xdata)
          y = np.load(ydata)
          color = ['red', 'green', 'blue']
          line = ['dashdot', 'dotted', 'dashed']
          for i in range(0, steps, slices):
              pic = 0
              for n in range(Nobj):
                  if n == 1:
                     marker = 'o'
                  else:
                     marker = '*'
                  plt.plot(x[n, : i], y[n, : i], color = color[n], linestyle= line[n], linewidth = 3, label = 'object ' + str(n+1))
                  plt.plot(x[n, i], y[n, i], marker = marker, color = color[n], markersize = 12)
              plt.title('Orbit plots of three objects', fontweight='bold')
              plt.xlabel('x-Coordinates (light second)')
              plt.ylabel('y-coordinate (light second)')
              plt.legend(loc="upper right")
              plt.savefig(fname + str(pic) + '.png')
              pic += 1
              plt.clf()

      def wob_plot(self, xdata, ydata, Nobj, steps, slices, lim, fname):
          """
          To show the wobbling nature of objects.

          Parameters :
          ----------

          xdata : str
                  The name of .npy file containing X-coordinate of objects.
          ydata : str
                  The name of .npy file containing Y-coordinate of objects. 
          Nobj  : int
                  The object to be plotted. obj=(1, 2, 3) for (first star, planet, second star)
          steps : int
                  The number of steps of observation to be ploted.
          slices: int
                  The steps to be ignored between two observation.
          lim   : list
                  The element number of minimum and maximum value of X and Y coordinate of given Nobj.
                  i.e. [xlim1, xlim2, ylim1, ylim2]
          fname : str
                  The name of output file in .png format.
                  
          Returns :
          -------
                  Return the .png file of orbit of objects according to obj.
          """
          plt.close()
          plt.rcParams.update({'font.size': 14})
          plt.rcParams["figure.figsize"] = [12,12]
          plt.rcParams["font.weight"] = "bold"
          plt.rcParams["axes.labelweight"] = "bold"
          x = np.load(xdata)
          y = np.load(ydata)
          color = ['red', 'green', 'blue']
          pic = 0
          for i in range(0, steps, slices):
              plt.plot(x[Nobj-1, : i], y[Nobj-1, : i], color = color[Nobj-1], linestyle= 'solid', linewidth = 1.5)
              if Nobj == 2:
                 marker = 'o'
              else:
                 marker = '*'
              plt.plot(x[Nobj-1, i], y[Nobj-1, i], marker=marker, color = color[Nobj-1], markersize=15)
              plt.title('Wobble nature of taken object', fontweight='bold')
              plt.xlabel('x-Coordinates (light second)')
              plt.ylabel('y-coordinate (light second)')
              plt.xlim(x[Nobj-1, lim[0]], x[Nobj-1, lim[1]])
              plt.ylim(y[Nobj-1, lim[2]], y[Nobj-1, lim[3]])
              plt.savefig(fname + str(pic) + '.png')
              pic += 1
              plt.clf()


