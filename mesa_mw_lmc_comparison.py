import numpy as np
import matplotlib.pyplot as plt
import os
home = os.path.expanduser('~')

dirpath = '/agb_star_and_dust/data/'
dirpath_mw = '/agb_star_and_dust/data/agb_mw/'
plotpath = '/agb_star_and_dust/plots/'

fig1 = plt.figure(figsize=(8,6))   #Triple alpha comparison
ax1 = fig1.add_subplot(111)
fig2 = plt.figure(figsize=(8,6))   #Star mass comparison
ax2 = fig2.add_subplot(111)
fig3 = plt.figure(figsize=(8,6))   #log Center T comparison
ax3 = fig3.add_subplot(111)
fig4 = plt.figure(figsize=(8,6))   #HR Diagram
ax4 = fig4.add_subplot(111)
fig5 = plt.figure(figsize=(8,6))   #star age vs log R
ax5 = fig5.add_subplot(111)

filename = home + dirpath + '0.5m_agb.data'
data05 = np.genfromtxt(filename, skip_header=5, dtype='str')
filename = home + dirpath + '0.6m_agb.data'
data06 = np.genfromtxt(filename, skip_header=5, dtype='str')
filename = home + dirpath + '0.7m_agb.data'
data07 = np.genfromtxt(filename, skip_header=5, dtype='str')
filename = home + dirpath + '0.8m_agb.data'
data08 = np.genfromtxt(filename, skip_header=5, dtype='str')
filename = home + dirpath + '0.9m_agb.data'
data09 = np.genfromtxt(filename, skip_header=5, dtype='str')
filename = home + dirpath + '1m_agb.data'
data10 = np.genfromtxt(filename, skip_header=5, dtype='str')
filename = home + dirpath + '2m_agb.data'
data20 = np.genfromtxt(filename, skip_header=5, dtype='str')
filename = home + dirpath + '3m_agb.data'
data30 = np.genfromtxt(filename, skip_header=5, dtype='str')
filename = home + dirpath + '4m_agb.data'
data40 = np.genfromtxt(filename, skip_header=5, dtype='str')
filename = home + dirpath + '5m_agb.data'
data50 = np.genfromtxt(filename, skip_header=5, dtype='str')

filename = home + dirpath_mw + '0.5m_agb.data'
data05_mw = np.genfromtxt(filename, skip_header=5, dtype='str')
filename = home + dirpath_mw + '0.6m_agb.data'
data06_mw = np.genfromtxt(filename, skip_header=5, dtype='str')
filename = home + dirpath_mw + '0.7m_agb.data'
data07_mw = np.genfromtxt(filename, skip_header=5, dtype='str')
filename = home + dirpath_mw + '0.8m_agb.data'
data08_mw = np.genfromtxt(filename, skip_header=5, dtype='str')
filename = home + dirpath_mw + '0.9m_agb.data'
data09_mw = np.genfromtxt(filename, skip_header=5, dtype='str')
filename = home + dirpath_mw + '1m_agb.data'
data10_mw = np.genfromtxt(filename, skip_header=5, dtype='str')
filename = home + dirpath_mw + '2m_agb.data'
data20_mw = np.genfromtxt(filename, skip_header=5, dtype='str')
filename = home + dirpath_mw + '3m_agb.data'
data30_mw = np.genfromtxt(filename, skip_header=5, dtype='str')
filename = home + dirpath_mw + '4m_agb.data'
data40_mw = np.genfromtxt(filename, skip_header=5, dtype='str')
filename = home + dirpath_mw + '5m_agb.data'
data50_mw = np.genfromtxt(filename, skip_header=5, dtype='str')

linecolor = ['Purple','Indigo','Blue','Orange','OrangeRed','Red',]
label = [r'$\mathrm{0.5 M_{\odot}~(LMC)}$',r'$\mathrm{1 M_{\odot}~(LMC)}$',r'$\mathrm{5 M_{\odot}~(LMC)}$',r'$\mathrm{0.5 M_{\odot}~(MW)}$',r'$\mathrm{1 M_{\odot}~(MW)}$',r'$\mathrm{5 M_{\odot}~(MW)}$']

name = data05[0,:]
#Plot the triple alpha reaction rate
x = np.where(name == 'star_age')
y = np.where(name == 'tri_alfa')
ax1.plot(data05[1:,x].astype('float').flat,data05[1:,y].astype('float').flat,'-',color=linecolor[0],label=label[0])
ax1.plot(data10[1:,x].astype('float').flat,data10[1:,y].astype('float').flat,'-',color=linecolor[1],label=label[1])
ax1.plot(data50[1:,x].astype('float').flat,data50[1:,y].astype('float').flat,'-',color=linecolor[2],label=label[2])
ax1.plot(data05_mw[1:,x].astype('float').flat,data05_mw[1:,y].astype('float').flat,'-',color=linecolor[3],label=label[3])
ax1.plot(data10_mw[1:,x].astype('float').flat,data10_mw[1:,y].astype('float').flat,'-',color=linecolor[4],label=label[4])
ax1.plot(data50_mw[1:,x].astype('float').flat,data50_mw[1:,y].astype('float').flat,'-',color=linecolor[5],label=label[5])
ax1.set_xlabel(r'$\mathrm{Star~Age~(yr)}$',fontsize=14)
ax1.set_ylabel(r'$\mathrm{Triple-\alpha~Reaction~Rate}$',fontsize=14)
ax1.set_xscale('log')
ax1.set_xlim([1e5,1e11])
ax1.legend(loc='upper left',numpoints=1)
fig1.savefig(home+plotpath+'triple_alpha_mw_lmc_comparison.eps',format='eps',dpi=300)
ax1.cla()
fig1.clf()
#Plot the star mass comparison
x = np.where(name == 'star_age')
y = np.where(name == 'star_mass')
ax2.plot(data05[1:,x].astype('float').flat,data05[1:,y].astype('float').flat,'-',color=linecolor[0],label=label[0])
ax2.plot(data10[1:,x].astype('float').flat,data10[1:,y].astype('float').flat,'-',color=linecolor[1],label=label[1])
ax2.plot(data50[1:,x].astype('float').flat,data50[1:,y].astype('float').flat,'-',color=linecolor[2],label=label[2])
ax2.plot(data05_mw[1:,x].astype('float').flat,data05_mw[1:,y].astype('float').flat,'-',color=linecolor[3],label=label[3])
ax2.plot(data10_mw[1:,x].astype('float').flat,data10_mw[1:,y].astype('float').flat,'-',color=linecolor[4],label=label[4])
ax2.plot(data50_mw[1:,x].astype('float').flat,data50_mw[1:,y].astype('float').flat,'-',color=linecolor[5],label=label[5])
ax2.set_xlabel(r'$\mathrm{Star~Age~(yr)}$',fontsize=14)
ax2.set_ylabel(r'$\mathrm{Star~Mass (M_{\odot})}$',fontsize=14)
ax2.set_xscale('log')
ax2.set_ylim([0,6])
ax2.legend(loc='lower left',numpoints=1)
fig2.savefig(home+plotpath+'star_mass_mw_lmc_comparison.eps',format='eps',dpi=300)
ax2.cla()
fig2.clf()
#Plot the star age vs log center T
x = np.where(name == 'star_age')
y = np.where(name == 'log_center_T')
ax3.plot(data05[1:,x].astype('float').flat,data05[1:,y].astype('float').flat,'-',color=linecolor[0],label=label[0])
ax3.plot(data10[1:,x].astype('float').flat,data10[1:,y].astype('float').flat,'-',color=linecolor[1],label=label[1])
ax3.plot(data50[1:,x].astype('float').flat,data50[1:,y].astype('float').flat,'-',color=linecolor[2],label=label[2])
ax3.plot(data05_mw[1:,x].astype('float').flat,data05_mw[1:,y].astype('float').flat,'-',color=linecolor[3],label=label[3])
ax3.plot(data10_mw[1:,x].astype('float').flat,data10_mw[1:,y].astype('float').flat,'-',color=linecolor[4],label=label[4])
ax3.plot(data50_mw[1:,x].astype('float').flat,data50_mw[1:,y].astype('float').flat,'-',color=linecolor[5],label=label[5])
ax3.set_xlabel(r'$\mathrm{Star~Age~(yr)}$',fontsize=14)
ax3.set_ylabel(r'$\mathrm{log~T_{center}}$',fontsize=14)
ax3.set_xscale('log')
ax3.legend(loc='lower left',numpoints=1)
fig3.savefig(home+plotpath+'log_center_T_mw_lmc_comparison.eps',format='eps',dpi=300)
ax3.cla()
fig3.clf()
#Plot the HR diagram
x = np.where(name == 'star_age')
y = np.where(name == 'tri_alfa')
lmc = [data05,data06,data07,data08,data09,data10,data20,data30,data40,data50]
mw = [data05_mw,data06_mw,data07_mw,data08_mw,data09_mw,data10_mw,data20_mw,data30_mw,data40_mw,data50_mw]
mass = [0.5,0.6,0.7,0.8,0.9,1.0,2.0,3.0,4.0,5.0]
hf_time_lmc = []
hf_mass_lmc = []
hf_time_mw = []
hf_mass_mw = []
for i in range(0,len(mass)):
	data = lmc[i]
	helium_abundance = data[1:,y].astype('float')
	ind = np.where(helium_abundance == max(helium_abundance))
	hf_time_lmc.append((data[1:,x].astype('float'))[ind])
	hf_mass_lmc.append(mass[i])
for i in range(0,len(mass)):
	data = mw[i]
	helium_abundance = data[1:,y].astype('float')
	ind = np.where(helium_abundance == max(helium_abundance))
	hf_time_mw.append((data[1:,x].astype('float'))[ind])
	hf_mass_mw.append(mass[i])
ax4.plot(hf_mass_lmc,hf_time_lmc,'-',color='Blue',label=r'$\mathrm{Z_{LMC}}$')
ax4.plot(hf_mass_mw,hf_time_mw,'-',color='Red',label=r'$\mathrm{Z_{MW}}$')
ax4.set_xlabel(r'$\mathrm{Mass~(M_{\odot})}$',fontsize=14)
ax4.set_ylabel(r'$\mathrm{Time~(yr)}$',fontsize=14)
ax4.set_yscale('log')
ax4.set_xscale('log')
ax4.legend(loc='upper right',numpoints=1)
ax4.set_xlim([0.3,6])
fig4.savefig(home+plotpath+'helium_flash_time_mw_lmc_comparison.eps',format='eps',dpi=300)
ax4.cla()
fig4.clf()
#Plot the star age vs log R
x = np.where(name == 'star_age')
y = np.where(name == 'log_R')
ax5.plot(data05[1:,x].astype('float').flat,data05[1:,y].astype('float').flat,'-',color=linecolor[0],label=label[0])
ax5.plot(data10[1:,x].astype('float').flat,data10[1:,y].astype('float').flat,'-',color=linecolor[1],label=label[1])
ax5.plot(data50[1:,x].astype('float').flat,data50[1:,y].astype('float').flat,'-',color=linecolor[2],label=label[2])
ax5.plot(data05_mw[1:,x].astype('float').flat,data05_mw[1:,y].astype('float').flat,'-',color=linecolor[3],label=label[3])
ax5.plot(data10_mw[1:,x].astype('float').flat,data10_mw[1:,y].astype('float').flat,'-',color=linecolor[4],label=label[4])
ax5.plot(data50_mw[1:,x].astype('float').flat,data50_mw[1:,y].astype('float').flat,'-',color=linecolor[5],label=label[5])
ax5.set_xlabel(r'$\mathrm{Star~Age~(yr)}$',fontsize=14)
ax5.set_ylabel(r'$\mathrm{log~R}$',fontsize=14)
ax5.set_xscale('log')
ax5.legend(loc='lower left',numpoints=1)
fig5.savefig(home+plotpath+'log_r_mw_lmc_comparison.eps',format='eps',dpi=300)
ax5.cla()
fig5.clf()