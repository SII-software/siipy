import numpy as np
from orbit import orb_int
from orbitplot import orbplot
from obstime import observe
from orbitfit import fit_orb
from aperture import aper
from obspoint import vary_base
from model import modeling
from dataplot import dplot

################## The orbital integration ##########
orbint = orb_int()
orbint.obs_para(avgt = 900, step = 7540)
obj1 = ['28.5']                                                                  
obj2 = ['9.0', '78.53', '1.2', '0.334', '65.5', '247.7', '2450120.9', '67.4']    
obj3 = ['14', '1.6e8', '13791', '0.326', '65', '120', '2450120.9', '6']           
N = 3
now = 2460311.23264          # 1st Jan 2024, 7 pm
fname = "../html/html/regor/orbit/para.txt" 
fname1 = "../html/html/regor/orbit/ini.txt" 
xname = '../html/html/regor/orbit/x_coord.npy'
yname = '../html/html/regor/orbit/y_coord.npy'
zname = '../html/html/regor/orbit/z_coord.npy'
tname = '../html/html/regor/orbit/t_coord.npy'
vname = '../html/html/regor/orbit/v_coord.npy'         
#a = orbint.parameter(fname, N, obj1, obj2, obj3)
#b = orbint.initial(fname, fname1, now, N) 
#c = orbint.intgration(fname1, N, xname, yname, zname, tname, vname)

##################### the orbit plot #################
ob = orbplot()
obs = observe()
vo = fit_orb()
start = [2460311.23264, 2460312.23264, 2460313.23264, 2460314.23264, 2460315.23264, 
         2460316.23264, 2460317.23264, 2460318.23264, 2460319.23264, 2460320.23264, 
         2460321.23264, 2460322.24097, 2460323.27986, 2460324.32292, 2460325.36875, 
         2460326.41458, 2460327.45972]
end = [2460311.40208, 2460312.44028, 2460313.47708, 2460314.51389, 2460315.55000, 
       2460316.58611, 2460317.62222, 2460318.65972, 2460319.67292, 2460320.67292, 
       2460321.67292, 2460322.67292, 2460323.67292, 2460324.67292, 2460325.67292, 
       2460326.67292, 2460327.67292]
pname = '../html/html/regor/orbit/binary2'
rname = '../html/html/regor/orbit/orbit2'
vel = '../html/html/regor/orbit/vel'
pos = '../html/html/regor/orbit/pos/pos'
orb = '../html/html/regor/orbit/orb/orb'
#ob.pos_plot(xname, yname, 2, pname, rname)  
#ob.vel_plot(tname, vname, 2, orbint.step, vel) 
step = obs.obslen(start, end, 0.0104167)     
dst, den = obs.nobs(start, step, 0.0104167)  
jd = obs.julday(start, end, step)
#ob.c_plot(xname, yname, dst, step, 0, pos)
X, Y, t = vo.orbit(dst, den, 0, 1, xname, yname, tname) 
zs = np.zeros(8)
#vo.resi_poly(zs, t, X, Y)

#################### The Intensity Interferometry integration ##########
# telescope position, name, baseline name and baseline length
ap = aper()
x = np.linspace(-100, 100, 500)
y = np.linspace(-100, 100, 500)
Tname = ['T1', 'T2', 'T3', 'T4'] 
Tpos = [-0.16-85.04j, 85.07-0.37j, 0.24+85.04j, -85.04+0.28j]
fname = 'regor/intfer/telescope.png'
#ap.telescope(Tpos, Tname, fname, x, y, 6) 
tel, base = ap.baseline(Tpos, Tname) 
x = np.real(base) 
y = np.imag(base)
z = 1e-6

# variational baseline according to the earth rotation
vb = vary_base()
la = vb.raddeg(-23, 16, 17)
lo = vb.raddeg(16, 30, 00)
r = vb.radhr(8, 9, 31.95013)
de = vb.raddeg(-47, 20, 11.7108)
vb.ps_para(lat=la, lon=lo, rac=r, dec=de)
dist = []                  # the baseline position according to earth rotation
for i in range(len(x)):
    dist.append(vb.rotbase(x[i], y[i], z, jd))

distance = np.array(dist)
xt = distance[:,0]         # baselines in east direction
yt = distance[:,1]         # baselines in north direction
zt = distance[:,2]         # baseline in up

grid = []                  # grids for each baselines as earth rotate
for i in range(len(x)):
    grid.append(vb.grids(xt[i,:], yt[i,:], 6, len(jd)))

tgrid = np.array(grid)
gX = tgrid[:,0]            # X axis grids of every baseline
gY = tgrid[:,1]            # Y axis grids of every baseline
wX = tgrid[:,2]            # X axis grids of telescope
wY = tgrid[:,3]            # Y axis grids of telescope

md = modeling()
dfname = '../html/html/regor/intfer/sdata_c.npy'
dfname1 = '../html/html/regor/intfer/mdata_c.npy'
md.sii_para(Area = 100, delt = 5e-9, itta = 0.30, loss = 0.50, avgt = orbint.avgt, chnl=1)
md.star_para(Tmp_a = 35000, Tmp_b = 57000, Tmp_r = 13880, dis_e = 3.458e10, lam_c = 465e-9, lam_e = 465e-9)
# create the simulated noisy and model data
XY = [X, Y]
bs = [xt, yt]
uv = [gX, gY]
pa = [md.Tmp_a, md.Tmp_b, md.Tmp_r, md.dis_e, md.lam_c, md.lam_e]
sf = [39.4503, 13.9236, 82, 0.6]
zs = [-2.66085443e+02, 3.42304944e+02, -3.31245787e-04, -2.70726990e-05, 5.25464975e-11, -1.70730610e-10, 1.16312155e-16, 9.31503649e-18]
pl = [zs, t]

# create the simulated noisy and model data
#data = md.smdata(bs, uv, sf, jd, 1, pl=pl)
#data1 = md.smdata(bs, uv, sf, jd, 0, pl=pl)
#np.save(dfname, data)  
#np.save(dfname1, data1)

# plot the simulated noisy and model data
dtp = dplot()
fname = '../html/html/regor/intfer/sdata/cntr/cntr'
color = ['blue', 'orange', 'black', 'cyan', 'red', 'purple']
#dtp.cntrplot(pa, sf, XY, bs, step, tel, color, fname)
fname = '../html/html/regor/intfer/sdata/sig/sig'
#dtp.sigplot(dfname, step, tel, jd, color, fname) 
fname = '../html/html/regor/intfer/sdata/comp/comp'
#dtp.cmplot(dfname, pa, sf, XY, bs, step, tel, jd, color, fname)

# estimate the parameter
fname2 = '../html/html/regor/intfer/dres_c_all.pkl'
pnames = ('R_a (ls)', 'R_b (ls)', 'R_r (ls)', 'x0', 'y0', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3')
fname = [fname2, dfname, pnames]
p_min = [0.5*sf[0], 0.5*sf[1], 0.5*sf[2], 0.5*zs[0], 0.5*zs[1], 0.5*zs[2], 0.5*zs[3], 0.5*zs[4], 0.5*zs[5], 0.5*zs[6], 0.5*zs[7]]
p_max = [1.5*sf[0], 1.5*sf[1], 1.5*sf[2], 1.5*zs[0], 1.5*zs[1], 1.5*zs[2], 1.5*zs[3], 1.5*zs[4], 1.5*zs[5], 1.5*zs[6], 1.5*zs[7]]
para = 11
md.estimate(fname, bs, uv, sf, p_min, p_max, para, pl=pl)  # estimate the parameters


# Corner plot for each parameter
truth = [sf[0], sf[1], sf[2], zs[0], zs[1], zs[2], zs[3], zs[4], zs[5], zs[6], zs[7]]
cfile = '../html/html/regor/intfer/para_c_all'
dtp.para_plot(truth, pnames, fname2, cfile)

