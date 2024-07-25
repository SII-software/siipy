import numpy as np
from orbit import orb_int
from orbitplot import orbplot
from aperture import aper
from obstime import observe
from obspoint import vary_base
from orbitfit import fit_orb
from model import modeling
from checkmask import chkmsk
from dataplot import dplot

################### The orbital integration ############
orbint = orb_int()
orbint.obs_para(avgt = 1000000, step = 18000)
obj1 = ['1.1']          
obj2 = ['1e-2', '365.25', '1', '0', '80', '0', '0', '0'] 
obj3 = ['0.9', '29183.475', '23.5', '0.518', '79', '204.85', '2405889.5', '231.65']
obj4 = ['0.1221', '199791750', '8700', '0.5', '107', '126', '105815172.5', '72'] 
N = 4
now = 2460311.29167         # 1st Jan 2024, 7 pm
fname = "../draft_material/html/wide/orbit/para.txt"   
fname1 = "../draft_material/html/wide/orbit/ini.txt"  
xname = '../draft_material/html/wide/orbit/x_coord.npy'
yname = '../draft_material/html/wide/orbit/y_coord.npy'
zname = '../draft_material/html/wide/orbit/z_coord.npy'
tname = '../draft_material/html/wide/orbit/t_coord.npy'
vname = '../draft_material/html/wide/orbit/v_coord.npy'     
#a = orbint.parameter(fname, N, obj1, obj2, obj3, obj4) 
#b = orbint.initial(fname, fname1, now, N) 
#c = orbint.intgration(fname1, N, xname, yname, zname, tname, vname) 

###################### the orbit plot ###################
ob = orbplot()
pname = '../draft_material/html/wide/orbit/binary4' 
rname = '../draft_material/html/wide/orbit/orbit4' 
vel = '../draft_material/html/wide/orbit/vel4'  
pos = '../draft_material/html/wide/orbit/pos/pos' 
wob = '../draft_material/html/wide/orbit/wobP/wobP'
#ob.pos_plot(xname, yname, 4, pname, rname) 
#ob.vel_plot(tname, vname, 4, 100, vel)
ob.w_plot(xname, yname, 3, 4000, 20, pos)
x = np.load(xname)
y = np.load(yname)
lim = [x[0, : 100].argmin(), x[0, : 100].argmax(), y[0, : 100].argmin(), y[0, : 100].argmax()]
ob.wob_plot(xname, yname, 1, 100, 1, lim, wob)

####### The Intensity Interferometry integration ##########
# the telescopes and baseline
ap = aper()
x = np.linspace(0.5, 3.5, 1000)
y = np.linspace(1.5, 5.5, 1000)
Tname = ['T1', 'T2', 'T3', 'T4']
Tpos = [2+2j, 2+5j, 1+4j, 3+4j]
fname = '../draft_material/html/wide/intfer/telescope.png'
#ap.telescope(Tpos, Tname, fname, x, y, 0.2, width=0.08, orien=1.82)
tel, base = ap.baseline(Tpos, Tname)
x = np.real(base)
y = np.imag(base)
z = 1e-6

# variational baseline
vb = vary_base()
obs = observe()
vo = fit_orb()
la = vb.raddeg(-43,59,12.0) 
lo = vb.raddeg(170,27,54.0)
r = vb.radhr(14,39,35.06311)
de = vb.raddeg(-60,50,15.0992)
vb.ps_para(lat=la, lon=lo, rac=r, dec=de)
start = [2460311.29167]  
end = [2460311.70833] # .70833
step = obs.obslen(start, end, 1.15741e-5)
dst, den = obs.nobs(start, step, 1.15741e-5) 
jd = obs.julday(start, end, step)
X, Y, t = vo.orbit(dst, den, 0, 2, xname, yname, tname)

dist = []                           # the baseline position according to earth rotation
for i in range(len(x)):
    dist.append(vb.rotbase(x[i], y[i], z, jd))

distance = np.array(dist)
xt = distance[:,0]                 # baselines in east direction
yt = distance[:,1]                 # baselines in north direction
zt = distance[:,2]                 # baseline in up

grid = []                          # grids for each baselines as earth rotate
for i in range(len(x)):
    grid.append(vb.grids(xt[i,:], yt[i,:], 0.2, 1024))  # radius of aperture 0.2 m

tgrid = np.array(grid)
gX = tgrid[:,0]                    # X axis grids of every baseline
gY = tgrid[:,1]                    # Y axis grids of every baseline
wX = tgrid[:,2]                    # X axis grids of telescope
wY = tgrid[:,3]                    # Y axis grids of telescope

# the SII and wide binary parameters
md = modeling()
wm = chkmsk()
md.sii_para(Area = 0.063, delt = 4e-10, itta = 0.30, loss = 0.50, avgt = 1, chnl=100)
md.star_para(Tmp_a = 5790, Tmp_b = 5260, Tmp_r = 0, dis_e = 0.14e9, lam_c = 6e-7, lam_e = 5e-7)
pa = [md.Tmp_a, md.Tmp_b, md.Tmp_r, md.dis_e, md.lam_c, md.lam_e]
sf = [2.94, 2.07, 0, 0]
XY = [X[0], Y[0]]

# Check the mask's width and orientation and accordingly the simulated data
# to check the mask's parameter, keep observational time for some minutes only
bs = [xt[0], yt[0]]
uv = [gX[0], gY[0]]
fname = '../draft_material/html/wide/intfer/hbt.png'
#wm.chbt(fname, uv, pa, XY, sf)

apt = [wX[0], wY[0], 0.2, 0.0146, 1.456]
lim = [-0.4, 0.4, -0.4, 0.4]
fname = '../draft_material/html/wide/intfer/aper.png'
wm.aplot(fname, apt, lim)
fname = '../draft_material/html/wide/intfer/sig.png'
#wm.gplot(fname, apt, uv, pa, XY, sf)

fname = '../draft_material/html/wide/intfer/signal/sig'
lim = [0.4, 0.8]
apt = [wX[0], wY[0], 0.2, [0.0140, 0.0150], [1.450, 1.460]]
#wm.cmask(fname, apt, uv, bs, pa, sf, XY, jd, lim)

fname = '../draft_material/html/wide/intfer/gif/sig'
lim = [0.4, 0.8]
apt = [wX[0], wY[0], 0.2, [0.0140, 0.0150], 1.456]
#wm.wmask(fname, apt, uv, bs, pa, sf, XY, jd, lim)

fname = '../draft_material/html/wide/intfer/gif1/sig'
lim = [0.4, 0.8]
apt = [wX[0], wY[0], 0.2, 0.0146, [1.450, 1.460]]
#wm.omask(fname, apt, uv, bs, pa, sf, XY, jd, lim)

fname = '../draft_material/html/wide/intfer/varwdth'
b = [0.0140, 0.0142, 0.0144, 0.0146, 0.0148, 0.0150, 0.0152]
lim = [0.4, 0.8]
apt = [wX[0], wY[0], 0.2, b, 1.456]
color = ['black', 'green', 'cyan', 'blue', 'red', 'deeppink', 'lime']
#wm.var_w(fname, apt, uv, bs, pa, sf, XY, jd, lim, color)

fname = '../draft_material/html/wide/intfer/varorn'
psi = [1.40, 1.41, 1.42, 1.436, 1.44, 1.45, 1.46]
lim = [0.4, 0.8]
apt = [wX[0], wY[0], 0.2, 0.0146, psi]
color = ['black', 'green', 'cyan', 'blue', 'red', 'deeppink', 'lime']
#wm.var_o(fname, apt, uv, bs, pa, sf, XY, jd, lim, color)

# save the simulated noisy and model data 
dfname = '../draft_material/html/wide/intfer/sdata.npy'                
dfname1 = '../draft_material/html/wide/intfer/mdata.npy'
bs = [xt, yt]
uv = [gX, gY]
apt = [wX, wY, 0.2, 0.0146, 1.456]
#data = md.smdata(bs, uv, sf, jd, 1, XY=XY, apt=apt)
#data1 = md.smdata(bs, uv, sf, jd, 0, XY=XY, apt=apt)
#np.save(dfname, data)  
#np.save(dfname1, data1) 

# plot the data for given telescopes
dp = dplot()
fname = '../draft_material/html/wide/intfer/mdata'
color = ['blue', 'orange', 'cyan', 'olive', 'red', 'purple']
#dp.allsigplot(dfname1, fname, jd, tel, color)

fname = '../draft_material/html/wide/intfer/spflux'
T = [5790, 5260]
R = [2.94, 2.07]
D = 0.14e9
nu = np.linspace(4.2e14,7.8e14,100)
#dp.spflux(fname, T, R, D, nu)

fname = '../draft_material/html/wide/intfer/bin'
sdata = np.load(dfname)
mdata = np.load(dfname1)
#dp.bindt(jd[0:1000], sdata[0, :1000], mdata[0, :1000], fname, 45, 20) 

fname2 = '../draft_material/html/wide/intfer/dres.pkl'
pnames = ('X (ls)', 'Y (ls)')
fname = [fname2, dfname, pnames]
p_min = [0.5 * XY[0], 0.5 * XY[1]]
p_max = [1.5 * XY[0], 1.5 * XY[1]]
para = 2
#md.estimate(fname, bs, uv, sf, p_min, p_max, para, apt=apt)

# Corner plot for each parameter
truth = [X[0], Y[0]]
cfile = '../draft_material/html/wide/intfer/para'
#dp.para_plot(truth, pnames, fname2, cfile)

fname = '../draft_material/html/wide/intfer/diff_len/dres/*dres*'
fname1 = 'wide/intfer/errorx'
fname2 = 'wide/intfer/errory'
night = np.arange(1,36,1)
#dp.error_plot(night, fname, fname1, fname2)

dfname = '../draft_material/html/wide/intfer/diff_mass/dres/dres1.pkl'
fname = '../draft_material/html/wide/intfer/lsq1'
title = 'Mass of Planet =  3330.3 $M_\oplus$'
zs = np.zeros(14)
os = [2.5e-7,1,1,1,1]
ts = np.load(tname)
x = np.load(xname)
y = np.load(yname)
X = x[0,:200] - x[2,:200]
Y = y[0,:200] - y[2,:200]
t = ts[:200,0]
#vo.worb_fit(dfname, fname, title, zs, os, t, X, Y)

