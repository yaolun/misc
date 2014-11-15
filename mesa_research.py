import numpy as np
import matplotlib.pyplot as plt
import os
home = os.path.expanduser('~')

dirpath = ['/agb_star_and_dust/data/','/agb_star_and_dust/data/agb_mw/']
plotpath = '/agb_star_and_dust/plots/'

filepath = ['0.5m_agb','0.6m_agb','0.7m_agb','0.8m_agb','0.9m_agb','1m_agb','2m_agb','3m_agb','4m_agb','5m_agb']#,'6m_agb','7m_agb','8m_agb'
linecolor = ['Purple','DarkViolet','DarkSlateBlue','Blue','ForestGreen','DarkOrange','GoldenRod','IndianRed','OrangeRed','Red',]
label = [r'$\mathrm{0.5 M_{\odot}}$',r'$\mathrm{0.6 M_{\odot}}$',r'$\mathrm{0.7 M_{\odot}}$',r'$\mathrm{0.8 M_{\odot}}$',r'$\mathrm{0.9 M_{\odot}}$',r'$\mathrm{1 M_{\odot}}$',r'$\mathrm{2 M_{\odot}}$',r'$\mathrm{3 M_{\odot}}$',r'$\mathrm{4 M_{\odot}}$',r'$\mathrm{5 M_{\odot}}$',r'$\mathrm{6 M_{\odot}}$',r'$\mathrm{7 M_{\odot}}$',r'$\mathrm{8 M_{\odot}}$']
mass = np.array([0.5,0.6,0.7,0.8,0.9,1.0,2.0,3.0,4.0,5.0])
mass_density = 0.10498*mass**(-2.35)
time_bin = np.array([1e7,5e7,1e8,5e8,1e9,5e9,1e10,5e10,1e11])
suffix = ['lmc','mw']

for datapath in dirpath:
	fig1 = plt.figure(figsize=(8,6))   #HR diagram
	ax1 = fig1.add_subplot(111)
	fig2 = plt.figure(figsize=(8,6))   #star age vs star mass
	ax2 = fig2.add_subplot(111)
	fig3 = plt.figure(figsize=(8,6))   #Equation of states
	ax3 = fig3.add_subplot(111)
	fig4 = plt.figure(figsize=(8,6))   #star age vs center H
	ax4 = fig4.add_subplot(111)
	fig5 = plt.figure(figsize=(8,6))   #star age vs center He
	ax5 = fig5.add_subplot(111)
	fig6 = plt.figure(figsize=(8,6))   #star age vs triple alpha
	ax6 = fig6.add_subplot(111)
	fig7 = plt.figure(figsize=(8,6))   #time of helium flash of different mass stars
	ax7 = fig7.add_subplot(111)
	fig8 = plt.figure(figsize=(8,6))
	ax8 = fig8.add_subplot(111)

	hf_mass = []
	hf_time = []
	for path in filepath:
		filename = home+datapath+path+'.data'
		data = np.genfromtxt(filename, skip_header=5, dtype='str')
		print 'Load %s Data' % path
		name = data[0,:]
		#Plot the HR diagram
		x = np.where(name == 'log_Teff')
		y = np.where(name == 'log_L')
		ax1.plot(data[1:,x].astype('float').flat,data[1:,y].astype('float').flat,'-',color=linecolor[filepath.index(path)],label=label[filepath.index(path)])
		#Plot the star age versus star mass
		x = np.where(name == 'star_age')
		y = np.where(name == 'star_mass')
		ax2.plot((data[1:,x].astype('float')).flat,data[1:,y].astype('float').flat,'-',color=linecolor[filepath.index(path)],label=label[filepath.index(path)])
		#Plot the equation of states
		x = np.where(name == 'log_center_Rho')
		y = np.where(name == 'log_center_T')
		ax3.plot((data[1:,x].astype('float')).flat,data[1:,y].astype('float').flat,'-',color=linecolor[filepath.index(path)],label=label[filepath.index(path)])
		#Plot of star age and center hydrogen abundance
		x = np.where(name == 'star_age')
		y = np.where(name == 'center_h1')
		ind = np.nonzero(data[1:,y].astype('float').flat > 1e-10)
		ax4.plot(((data[1:,x].astype('float')))[ind].flat,(data[1:,y].astype('float'))[ind].flat,'-',color=linecolor[filepath.index(path)],label=label[filepath.index(path)])
		#Plot of star age and center helium abundance
		x = np.where(name == 'star_age')
		y = np.where(name == 'center_he4')
		ind = np.nonzero(data[1:,y].astype('float').flat > 0)
		ax5.plot(((data[1:,x].astype('float')))[ind].flat,(data[1:,y].astype('float'))[ind].flat,'-',color=linecolor[filepath.index(path)],label=label[filepath.index(path)])	
		#Plot of star age and the triple alpha reaction rate
		x = np.where(name == 'star_age')
		y = np.where(name == 'tri_alfa')
		ax6.plot(((data[1:,x].astype('float')))[ind].flat,(data[1:,y].astype('float'))[ind].flat,'-',color=linecolor[filepath.index(path)],label=label[filepath.index(path)])	
		#Plot the time of helium flash of different mass stars
		x = np.where(name == 'star_age')
		y = np.where(name == 'tri_alfa')
		helium_abundance = data[1:,y].astype('float')
		ind = np.where(helium_abundance == max(helium_abundance))
		hf_time.append((data[1:,x].astype('float'))[ind])
		hf_mass.append(mass[filepath.index(path)])
		#Plot the luminosity of H and triple alpha
		x = np.where(name == 'star_age')
		y = np.where(name == 'tri_alfa')
		y_lh = np.where(name == 'log_LH')
		#ax8.plot((data[1:,x].astype('float')).flat,data[1:,y].astype('float').flat,'-',color=linecolor[filepath.index(path)],label=label[filepath.index(path)])
		ax8.plot((data[1:,x].astype('float')).flat,data[1:,y_lh].astype('float').flat,'-',color=linecolor[filepath.index(path)],label=label[filepath.index(path)])
	ax1.set_xlabel(r'$\mathrm{log~T_{eff}}$',fontsize=14)
	ax1.set_ylabel(r'$\mathrm{log~L}$',fontsize=14)
	ax1.set_xlim([6,3])
	ax1.legend(loc='upper left',numpoints=1)

	ax2.set_xlabel(r'$\mathrm{Star~Age~(yr)}$',fontsize=14)
	ax2.set_ylabel(r'$\mathrm{Star~Mass~(M_{\odot})}$',fontsize=14)
	ax2.set_xlim([1e5,1e11])
	ax2.set_xscale('log')
	ax2.legend(loc='upper right',numpoints=1)

	ax3.set_xlabel(r'$\mathrm{log~\rho}$',fontsize=14)
	ax3.set_ylabel(r'$\mathrm{log~T}$',fontsize=14)
	ax3.legend(loc='upper left',numpoints=1)

	ax4.set_xlabel(r'$\mathrm{Star~Age~(yr)}$',fontsize=14)
	ax4.set_ylabel(r'$\mathrm{Center~Hydrogen~Mass~Fraction~(\%)}$',fontsize=14)
	ax4.set_xlim([0,1e11])
	ax4.legend(loc='upper right',numpoints=1)

	ax5.set_xlabel(r'$\mathrm{Star~Age~(yr)}$',fontsize=14)
	ax5.set_ylabel(r'$\mathrm{Center~Helium~Mass~Fraction~(\%)}$',fontsize=14)
	ax5.set_xlim([0,1e11])
	ax5.legend(loc='upper right',numpoints=1)

	ax6.set_xlabel(r'$\mathrm{Star~Age~(yr)}$',fontsize=14)
	ax6.set_ylabel(r'$\mathrm{Triple-\alpha~Reaction~Rate}$',fontsize=14)
	ax6.set_xscale('log')
	ax6.set_xlim([1e5,1e11])
	ax6.legend(loc='upper left',numpoints=1)

	#hf_time = np.log10(hf_time)
	#hf_mass = np.log10(hf_mass)
	ax7.plot(hf_mass,hf_time,'o')
	ax7.set_xlabel(r'$\mathrm{Mass~(M_{\odot})}$',fontsize=14)
	ax7.set_ylabel(r'$\mathrm{Time~(yr)}$',fontsize=14)
	ax7.set_yscale('log')
	ax7.set_xscale('log')
	ax7.set_xlim([0.3,10])

	ax8.set_xlabel(r'$\mathrm{Star~Age~(yr)}$',fontsize=14)
	ax8.set_ylabel(r'$\mathrm{Triple-\alpha~Rate~and~log~L_{H}}$',fontsize=14)
	#ax8.set_xlim([0,1e11])
	ax8.set_xscale('log')
	ax8.set_ylim([-5,10])
	#ax8.legend(loc='upper right',numpoints=1)
	print 'Change the axis'
	fig1.savefig(home+plotpath+suffix[dirpath.index(datapath)]+'_HR_diagram_all.eps',format='eps',dpi=300)
	fig2.savefig(home+plotpath+suffix[dirpath.index(datapath)]+'_star_age_star_mass_all.eps',format='eps',dpi=300)
	fig3.savefig(home+plotpath+suffix[dirpath.index(datapath)]+'_EoS_all.eps',format='eps',dpi=300)
	fig4.savefig(home+plotpath+suffix[dirpath.index(datapath)]+'_center_h_star_age.eps',format='eps',dpi=300)
	fig5.savefig(home+plotpath+suffix[dirpath.index(datapath)]+'_center_he_star_age.eps',format='eps',dpi=300)
	fig6.savefig(home+plotpath+suffix[dirpath.index(datapath)]+'_triple_alpha_star_age.eps',format='eps',dpi=300)
	fig7.savefig(home+plotpath+suffix[dirpath.index(datapath)]+'_helium_flash_time.eps',format='eps',dpi=300)
	fig8.savefig(home+plotpath+suffix[dirpath.index(datapath)]+'_triple_alpha_log_LH_star_age.eps',format='eps',dpi=300)
	ax1.cla()
	fig1.clf()
	ax2.cla()
	fig2.clf()
	ax3.cla()
	fig3.clf()
	ax4.cla()
	fig4.clf()
	ax5.cla()
	fig5.clf()
	ax6.cla()
	fig6.clf()
	ax7.cla()
	fig7.clf()
	ax8.cla()
	fig8.clf()
	print 'Finish %s' % suffix[dirpath.index(datapath)]
print 'Done!'