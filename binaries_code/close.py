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
orbint.obs_para(avgt = 900, step = 800)
obj1 = ['11.43']                                                                     
obj2 = ['7.21', '4.0145', '0.131', '0.133', '63.1', '309.938', '2440678.008', '255.6'] 
N = 2
now = 2460311.29167            # 1st Jan 2024, 7 pm
fname = "../draft_material/html/close/orbit/para.txt"       
fname1 = "../draft_material/html/close/orbit/ini.txt"
xname = '../draft_material/html/close/orbit/x_coord.npy'
yname = '../draft_material/html/close/orbit/y_coord.npy'
zname = '../draft_material/html/close/orbit/z_coord.npy'
tname = '../draft_material/html/close/orbit/t_coord.npy'
vname = '../draft_material/html/close/orbit/v_coord.npy' 
#a = orbint.parameter(fname, N, obj1, obj2)                                                                                                                  
#b = orbint.initial(fname, fname1, now, N)                                               
#c = orbint.intgration(fname1, N, xname, yname, zname, tname, vname)                     

##################### the orbit plot #################
ob = orbplot()
obs = observe()
vo = fit_orb()
start = [2460311.29167, 2460312.29167, 2460313.29167, 2460314.29167]   # (from Jan 1, 2024 to Jan 8, 2024, 7pm)
end = [2460311.70833, 2460312.70833, 2460313.70833, 2460314.70833]     # (from Jan 1, 2023 to Jan 8, 2024, 5am)
pname = '../draft_material/html/close/orbit/binary2'
rname = '../draft_material/html/close/orbit/orbit2'
vel = '../draft_material/html/close/orbit/vel'
pos = '../draft_material/html/close/orbit/pos/pos'
orb = '../draft_material/html/close/orbit/orb/orb'
fname = '../draft_material/html/close/orbit/fitorb'
#ob.pos_plot(xname, yname, 2, pname, rname) 
#ob.vel_plot(tname, vname, 2, orbint.step, vel)                               
step = obs.obslen(start, end, 0.0104167) 
dst, den = obs.nobs(start, step, 0.0104167)
jd = obs.julday(start, end, step)
#ob.c_plot(xname, yname, dst, step, 1, orb) 
X, Y, t = vo.orbit(dst, den, 0, 1, xname, yname, tname)
zs = [4, 66, 0.133, 1.1, 5, 2440678, 4] 
#vo.corb_fit(zs, jd, X, Y, fname)

##### The Intensity Interferometry integration ######
# telescope position, name, baseline name and baseline length
ap = aper()
x = np.linspace(-60, 150, 500)
y = np.linspace(-60, 80, 500)
Tname = ['T1', 'T2', 'T3', 'T4'] 
Tpos = [135 - 15j, 40 - 50j, 30 + 60j, -40 + 10j]
fname = '../draft_material/html/close/intfer/telescope.png'
#ap.telescope(Tpos, Tname, fname, x, y, 6)
tel, base = ap.baseline(Tpos, Tname)
x = np.real(base) 
y = np.imag(base)
z = 1e-6

# variational baseline according to the earth rotation
vb = vary_base()
la = vb.raddeg(31, 40, 30)
lo = vb.raddeg(110, 57, 7)
r = vb.radhr(13, 25, 11.579)
de = vb.raddeg(-11, 9, 40.75)
vb.ps_para(lat=la, lon=lo, rac=r, dec=de)
dist = []                    # the baseline position according to earth rotation
for i in range(len(x)):
    dist.append(vb.rotbase(x[i], y[i], z, jd))

distance = np.array(dist)
xt = distance[:,0]           # baselines in east direction
yt = distance[:,1]           # baselines in north direction
zt = distance[:,2]           # baseline in up

grid = []                    # grids for each baselines as earth rotate
for i in range(len(x)):
    grid.append(vb.grids(xt[i,:], yt[i,:], 6, len(jd)))

tgrid = np.array(grid)
gX = tgrid[:,0]              # X axis grids of every baseline
gY = tgrid[:,1]              # Y axis grids of every baseline
wX = tgrid[:,2]              # X axis grids of telescope
wY = tgrid[:,3]              # Y axis grids of telescope

md = modeling()
dfname = '../draft_material/html/close/intfer/sdata.npy' 
dfname1 = '../draft_material/html/close/intfer/mdata.npy'
c = 299792458
au = 149597870700 / c
deg = np.pi/180
md.sii_para(Area = 100, delt = 5e-9, itta = 0.30, loss = 0.50, avgt = orbint.avgt, chnl=1) 
md.star_para(Tmp_a = 25300, Tmp_b = 20900, Tmp_r = 0, dis_e = 7.89e9, lam_c = 570e-9, lam_e = 470e-9)

# create the simulated noisy and model data
XY = [X, Y]
bs = [xt, yt]
uv = [gX, gY]
pa = [md.Tmp_a, md.Tmp_b, md.Tmp_r, md.dis_e, md.lam_c, md.lam_e]
sf = [17.33, 8.68, 1e-25, 0.6]
kp = [4.0145, 0.131*au, 0.133, 63.1*deg, 309.938*deg, 2440678.008, 255.6*deg]
#data = md.smdata(bs, uv, sf, jd, 1, kp=kp)
#data1 = md.smdata(bs, uv, sf, jd, 0, kp=kp)
#np.save(dfname, data)  
#np.save(dfname1, data1)

# plot the simulated noisy and model data
dtp = dplot()
fname = '../draft_material/html/close/intfer/sdata/cntr/cntr'
color = ['blue', 'orange', 'cyan', 'olive', 'red', 'purple']
#dtp.cntrplot(pa, sf, XY, bs, step, tel, color, fname)
fname = '../draft_material/html/close/intfer/sdata/sig/sig'
#dtp.sigplot(dfname, step, tel, jd, color, fname)
fname = '../draft_material/html/close/intfer/sdata/comp/comp'
#dtp.cmplot(dfname, pa, sf, XY, bs, step, tel, jd, color, fname)

# estimate the parameter
pnames = ('R_a (ls)', 'R_b (ls)', 'ld', 'a (ls)', 'e', 'i (rad)', r'$\Omega (rad)$', r'$\omega (rad)$')
fname2 = '../draft_material/html/close/intfer/dres.pkl'
fname = [fname2, dfname, pnames]
p_min = [0.5*sf[0], 0.5*sf[1], 0.5*sf[3], 0.5*kp[1], 0.5*kp[2], 0.1*kp[3], 0.58*kp[4], 0.58*kp[6]]
p_max = [1.8*sf[0], 1.8*sf[1], 1.7*sf[3], 1.5*kp[1], 1.5*kp[2], 1.5*kp[3], 1.16*kp[4], 1.16*kp[6]]
para = 8 # it define the number of parameters to be estimated
#md.estimate(fname, bs, uv, sf, p_min, p_max, para, kp=kp, jlday=jd)  

# Corner plot for each parameter
truth = [sf[0], sf[1], sf[3], kp[1], kp[2], kp[3], kp[4], kp[6]]

cfile = '../draft_material/html/close/intfer/para'
dtp.para_plot(truth, pnames, fname2, cfile)

