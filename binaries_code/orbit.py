import scipy.optimize as sciopt
import numpy as np

class orb_int():
      """
      Create an object to do the integration for an orbit.
      It uses Kepler's equation to get the initial position and velocity of the system.
      For the integration of orbit, the leap-frog method is used.
      """
 
      def obs_para(self, avgt, step):
           """
           Set the parameters for integration of orbit.
          
           Parameters :
           ----------

           avgt : int
                  The time interval between two step.
           step : int
                  Number of steps for integration.

           Returns :
           -------
                  None.
           """
           self.avgt = avgt
           self.step = step

      def parameter(self, fname, Nobj, *obj):
           """
           Save all the orbital parameters of star system with their mass in .txt format.

           According to the number of objects in system, call obj1, obj2,....
           
           Parameters :
           ----------

           fname : str
                   The name of output file.
           Nobj  : int
                   Number of objects in system.
           *obj  : str
                   The number of objects in system i.e. *obj = obj1, obj2,......, objN.
                   The parameters are in order 
                   objN = ['M (Msun)', 'P (days)', 'a (au)', 'e', 'I (deg)', 'Omega (deg)', 'epoch (days)', 'omega (deg)'].
                   

           Returns :
           -------
                   fname.txt.
           """

           par =  "#M (Msun), P (days), a (au), e, I (deg), Omega (deg), epoch (days), omega (deg)"
           f = open(fname, "w")
           f.write(par)                                           
           f.write('\n')     
           for i in range(Nobj):
               for j in obj[i]:
                   f.write(j + str('\t')) 
               f.write('\n')              
           f.close() 
           
      # barycentric position of N objects.
      def bary(self, m, x, N):
           """       
           """
           cm = 0*m
           cm[0] = m[0]                                                            
           for n in range(1,N):
               cm[n] = cm[n-1] + m[n]                                              
           cx = 0*x                                                                
           for n in range(1,N):
               x[n] += cx[n-1]                                                     
               cx[n] = (cm[n-1]*cx[n-1] + m[n]*x[n])/cm[n];                        
           b = cx[N-1]                                                            
           x[:] -= b                                                               
           return x
           
      # The orbit in terms of node (Omega), inclination (I) and periapsis (omega).
      def rotate(self, x, y, omega, I, Omega):   
           """
           """              
           cs, sn = np.cos(omega), np.sin(omega)
           x, y, z = x*cs - y*sn, x*sn + y*cs, 0
           cs, sn = np.cos(I), np.sin(I)
           x,y,z = x, y*cs - z*sn, y*sn + z*cs
           Omega += np.pi/2                                                           # rotate (north, east)
           cs,sn = np.cos(Omega), np.sin(Omega)
           x,y,z = x*cs - y*sn, x*sn + y*cs, z
           return np.array((x,y,z))

      # The position and velocity of N objects according to Kepler's equation. 
      def posvel(self, P, a, e, omega, I, Omega, mean_an):
           """
           """
           ec = (1-e*e)**.5
           mean_an = mean_an % (2*np.pi)                                              
           psi = sciopt.brentq(lambda psi: psi - e*np.sin(psi) - mean_an, 0, 2*np.pi)    
           cs,sn = np.cos(psi), np.sin(psi)
           x,y = cs-e, ec*sn                                                       
           vx,vy = -sn/(1-e*cs), ec*cs/(1-e*cs)                                   
           pos = a * self.rotate(x,y,omega,I,Omega)                                     
           vel = 2*np.pi * a/P * self.rotate(vx,vy,omega,I,Omega)                          
           return pos, vel
           
      # Initial position and velocity of objects on observational time (now).
      def getposvel(self, fname, now, N):    
           """
           """                                
           fil = open(fname)
           fil.readline()
           star = []
           for n in range(N):
               pars = fil.readline().split()
               pars = [float(s) for s in pars]                                    
               star.append(pars)                                                   
           c = 299792458
           days = 86400
           au = 149597870700 / c
           deg = np.pi/180
           GMsun = 4.9254909e-6            
           mass = np.zeros(N)
           pos = np.zeros((N,3))      
           vel = np.zeros((N,3))           
           for n in range(N):
               mass[n] = float(star[n][0]) * GMsun                                 
           for n in range(1,N): 
               P,a,e,I,Omega,ep,omega = star[n][1:]                                
               mean_an = 2*np.pi * (now-ep)/P                                         
               P *= days                                                           
               a *= au                                                             
               I *= deg                                                            
               Omega *= deg                                                        
               omega *= deg                                                        
               pos[n],vel[n] = self.posvel(P,a,e,omega,I,Omega,mean_an)                 
           return mass, self.bary(mass,pos,N), self.bary(mass,vel,N)

      def initial(self, fname, fname1, now, Nobj):  
           """
           Calculate the initial position and velocity of N number of objects in system.

           It uses Kepler's equation.

           Parameters :
           ----------

           fname  : str
                    Input as fname.txt file with mass and orbital parameter.
           fname1 : str
                    Output as fname1.txt with initial position and velocity.
           now    : float
                    The julian day on the first observation.
           Nobj   : int
                    Number of objects in system.

           Returns :
           -------
                    fname1.txt.
           """                        
           mass, pos, vel = self.getposvel(fname, now, Nobj)                                
           fil = open(fname1, 'w')                             
           for l in range(Nobj):
               fil.write("%22.15e" %mass[l])                                     
               fil.write('\t')                       
           fil.write("\n")                         
           for l in range(Nobj):
               for k in range(3):                                                
                   fil.write("%22.15e" %pos[l,k])        
                   fil.write('\t')
               fil.write("\n")                           
           for l in range(Nobj):
               for k in range(3):                                                     
                   fil.write("%22.15e" %vel[l,k])        
                   fil.write('\t')
               fil.write("\n")                           
           fil.close()    
           print('center of mass position in light second = ',np.sum(mass*pos[:,0])/np.sum(mass))                                                        
           
      # Calculate total energy and ratio of Kinetic and potential energy of N object's system.
      def checks(self, mass, pos, vel, N):
           """
           """
           # center of mass for position and velocity
           cmx = 0*pos[0]
           cmv = 0*vel[0]
           for i in range(N):
               cmx += mass[i]*pos[i]
               cmv += mass[i]*vel[i]
           # kinetic energy and potential energy
           kin = pot = 0
           for i in range(N):
                kin += mass[i]*np.sum(vel[i]*vel[i])/2           
                for j in range(i):
                    dr = pos[i]-pos[j]
                    r = np.sum(dr*dr)**.5                        
                    pot -= mass[i]*mass[j]/r                     
           print('Total energy of the star system (kinetic + potential)', kin+pot)
           print('Ratio of kinetic and potential energy (Virial theorem)', kin/pot)                  
  
      def intgration(self, Iname, Nobj, *fname):  
           """ 
           Integrate the orbit of N objects using leap-frog method.

           Return five .npy format files (fname1, fname2, fname3, .......).

           Parameters :
           ----------

           Iname  : str
                    Input as Iname.txt with initial position and velocity.
           Nobj   : int
                    Number of objects in system.
           *fname : str
                    Output files in .npy format for X, Y, Z, and t coordinates also z-velocity of N objects.
                    *fname = fname1, fname2, fname3, fname4, fname5

           Returns :
           -------
                    X, Y, Z and time coordinates of position of N objects also the line of sight velocity.
                             
           """
           fil = open(Iname, "r")
           l = fil.readline().split() 
           mass = np.zeros(Nobj)                    
           for i in range(len(l)):                             
               mass[i] = float(l[i])                                
           pos = np.zeros(shape=(Nobj,3))
           for i in range(Nobj):
               l = fil.readline().split()
               for j in range(len(l)):
                   pos[i][j] = float(l[j])                           
           vel = np.zeros(shape=(Nobj,3))
           for i in range(Nobj):
               l = fil.readline().split()
               for j in range(len(l)):
                   vel[i][j] = float(l[j]) 

           # calculate the center of mass of position and velocity
           cmx = 0*pos[0]
           cmv = 0*vel[0]
           for i in range(Nobj):
               cmx += mass[i]*pos[i]
               cmv += mass[i]*vel[i]
           cmx /= np.sum(mass)                                      
           cmv /= np.sum(mass)   
           print('The center of mass of position :', cmx)
           print('The center of mass of velocity :', cmv)                                   
           
           # Barycentric position and velocity of each objects (pos and vel of Nth body w.r.t. center of mass)
           for i in range(Nobj):
               pos[i] -= cmx                                        
               vel[i] -= cmv                                        
                         
           # calculate the next position and velocity (z component) using Leapfrog algorithm
           x = np.zeros((Nobj, self.step), dtype = float)
           y = np.zeros((Nobj, self.step), dtype = float)
           z = np.zeros((Nobj, self.step), dtype = float)
           t_stp = np.zeros((self.step,1), dtype = float)
           vz = np.zeros((Nobj, self.step), dtype = float)
           for t in range(self.step):
               if t % self.avgt == 0:
                    self.checks(mass, pos, vel, Nobj)
               pos += vel*self.avgt/2                                      
               acc = 0*pos                                          
               t_stp[t] = t_stp[t-1] + self.avgt                            
               for i in range(Nobj):
                   for j in range(i):
                       dr = pos[i]-pos[j]                           
                       r2 = dr[0]**2 + dr[1]**2 + dr[2]**2          
                       r3 = r2**1.5                                  
                       acc[j,:] += mass[i]*dr/r3                    
                       acc[i,:] -= mass[j]*dr/r3                                
               vel += acc*self.avgt                                        
               pos += vel*self.avgt/2
               x[:,t] = pos[:,0]
               y[:,t] = pos[:,1]
               z[:,t] = pos[:,2]
               vz[:,t] = vel[:,2]                                   

           # save the position of N objects for n-coordinates also and time and z-velocity of N objects
           np.save(fname[0], x)           
           np.save(fname[1], y)           
           np.save(fname[2], z)           
           np.save(fname[3], t_stp)  
           np.save(fname[4], vz)                             

