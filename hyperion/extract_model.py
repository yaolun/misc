def extract_hyperion(filename,indir=None,outdir=None,dstar=178.0,wl_aper=None,save=True):
	def l_bol(wl,fv,dist=178.0):
		import numpy as np
		import astropy.constants as const
		# wavelength unit: um
		# Flux density unit: Jy
		#
		# constants setup
		#
		c = const.c.cgs.value
		pc = const.pc.cgs.value
		PI = np.pi
		SL = const.L_sun.cgs.value
		# Convert the unit from Jy to erg s-1 cm-2 Hz-1
		fv = np.array(fv)*1e-23
		freq = c/(1e-4*np.array(wl))
		
		diff_dum = freq[1:]-freq[0:-1]
		freq_interpol = np.hstack((freq[0:-1]+diff_dum/2.0,freq[0:-1]+diff_dum/2.0,freq[0],freq[-1]))
		freq_interpol = freq_interpol[np.argsort(freq_interpol)[::-1]]
		fv_interpol = np.empty(len(freq_interpol))
		# calculate the histogram style of spectrum
		#
		for i in range(0,len(fv)):
			if i == 0:
				fv_interpol[i] = fv[i]
			else:
				fv_interpol[2*i-1] = fv[i-1]
				fv_interpol[2*i] = fv[i]
		fv_interpol[-1] = fv[-1]
		
		dv = freq_interpol[0:-1]-freq_interpol[1:]
		dv = np.delete(dv,np.where(dv==0))

		fv = fv[np.argsort(freq)]
		freq = freq[np.argsort(freq)]

		return (np.trapz(fv,freq)*4.*PI*(dist*pc)**2)/SL


	import matplotlib.pyplot as plt
	import numpy as np
	import os
	from hyperion.model import ModelOutput
	from hyperion.model import Model
	from scipy.interpolate import interp1d
	from hyperion.util.constants import pc, c, lsun

	# Read in the observation data and calculate the noise & variance
	if indir == None:
		indir = '/Users/yaolun/bhr71/'
	if outdir == None:
		outdir = '/Users/yaolun/bhr71/hyperion/'

	# assign the file name from the input file
	print_name = os.path.splitext(os.path.basename(filename))[0]
	#
	[wl_pacs,flux_pacs,unc_pacs] = np.genfromtxt(indir+'BHR71_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_continuum.txt',\
										dtype='float',skip_header=1).T
	# Convert the unit from Jy to erg cm-2 Hz-1
	flux_pacs = flux_pacs*1e-23
	[wl_spire,flux_spire] = np.genfromtxt(indir+'BHR71_spire_corrected_continuum.txt',dtype='float',skip_header=1).T
	flux_spire = flux_spire*1e-23 
	wl_obs = np.hstack((wl_pacs,wl_spire))
	flux_obs = np.hstack((flux_pacs,flux_spire))

	[wl_pacs_data,flux_pacs_data,unc_pacs_data] = np.genfromtxt(indir+'BHR71_centralSpaxel_PointSourceCorrected_CorrectedYES_trim.txt',\
												  dtype='float').T
	[wl_spire_data,flux_spire_data] = np.genfromtxt(indir+'BHR71_spire_corrected.txt',\
													dtype='float').T

	[wl_pacs_flat,flux_pacs_flat,unc_pacs_flat] = np.genfromtxt(indir+'BHR71_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_flat_spectrum.txt',\
										dtype='float',skip_header=1).T
	[wl_spire_flat,flux_spire_flat] = np.genfromtxt(indir+'BHR71_spire_corrected_flat_spectrum.txt',dtype='float',skip_header=1).T

	# Convert the unit from Jy to erg cm-2 Hz-1
	flux_pacs_flat = flux_pacs_flat*1e-23 
	flux_spire_flat = flux_spire_flat*1e-23
	flux_pacs_data = flux_pacs_data*1e-23
	flux_spire_data = flux_spire_data*1e-23


	wl_pacs_noise = wl_pacs_data
	flux_pacs_noise = flux_pacs_data-flux_pacs-flux_pacs_flat
	wl_spire_noise = wl_spire_data
	flux_spire_noise = flux_spire_data-flux_spire-flux_spire_flat

	# Read in the Spitzer IRS spectrum
	[wl_irs, flux_irs]= (np.genfromtxt(indir+'bhr71_spitzer_irs.txt',skip_header=2,dtype='float').T)[0:2]
	# Convert the unit from Jy to erg cm-2 Hz-1
	flux_irs = flux_irs*1e-23
	# Remove points with zero or negative flux 
	ind = flux_irs > 0
	wl_irs = wl_irs[ind]
	flux_irs = flux_irs[ind]
	# Calculate the local variance (for spire), use the instrument uncertainty for pacs
	#
	wl_noise_5 = wl_spire_noise[(wl_spire_noise > 194)*(wl_spire_noise <= 304)]
	flux_noise_5 = flux_spire_noise[(wl_spire_noise > 194)*(wl_spire_noise <= 304)]
	wl_noise_6 = wl_spire_noise[wl_spire_noise > 304]
	flux_noise_6 = flux_spire_noise[wl_spire_noise > 304]
	wl_noise = [wl_pacs_data[wl_pacs_data<=190.31],wl_noise_5,wl_noise_6]
	flux_noise = [unc_pacs[wl_pacs_data<=190.31],flux_noise_5,flux_noise_6]
	sig_num = 20
	sigma_noise = []
	for i in range(0,len(wl_noise)):
		sigma_dum = np.zeros([len(wl_noise[i])])
		for iwl in range(0,len(wl_noise[i])):
			if iwl < sig_num/2:
				sigma_dum[iwl] = np.std(np.hstack((flux_noise[i][0:sig_num/2],flux_noise[i][0:sig_num/2-iwl])))
			elif len(wl_noise[i])-iwl < sig_num/2:
				sigma_dum[iwl] = np.std(np.hstack((flux_noise[i][iwl:],flux_noise[i][len(wl_noise[i])-sig_num/2:])))
			else:
				sigma_dum[iwl] = np.std(flux_noise[i][iwl-sig_num/2:iwl+sig_num/2])
		sigma_noise = np.hstack((sigma_noise,sigma_dum))
	sigma_noise = np.array(sigma_noise)

	# Read in the photometry data
	phot = np.genfromtxt(indir+'bhr71.txt',dtype=None,skip_header=1,comments='%')
	wl_phot = []
	flux_phot = []
	flux_sig_phot = []
	note = []
	for i in range(0,len(phot)):
		wl_phot.append(phot[i][0])
		flux_phot.append(phot[i][1])
		flux_sig_phot.append(phot[i][2])
		note.append(phot[i][4])
	wl_phot = np.array(wl_phot)
	# Convert the unit from Jy to erg cm-2 Hz-1
	flux_phot = np.array(flux_phot)*1e-23
	flux_sig_phot = np.array(flux_sig_phot)*1e-23

	# Print the observed L_bol
	wl_tot = np.hstack((wl_irs,wl_obs,wl_phot))
	flux_tot = np.hstack((flux_irs,flux_obs,flux_phot))
	flux_tot = flux_tot[np.argsort(wl_tot)]
	wl_tot = wl_tot[np.argsort(wl_tot)]
	l_bol_obs = l_bol(wl_tot,flux_tot*1e23)             


	# Open the model
	m = ModelOutput(filename)

	if wl_aper == None:
		wl_aper = [3.6, 4.5, 5.8, 8.0, 10, 16, 20, 24, 35, 70, 100, 160, 250, 350, 500, 850]

	# Create the plot
	mag = 1.5
	fig = plt.figure(figsize=(8*mag,6*mag))
	ax_sed = fig.add_subplot(1, 1, 1)

	# Plot the observed SED
	# plot the observed spectra
	pacs, = ax_sed.plot(np.log10(wl_pacs),np.log10(c/(wl_pacs*1e-4)*flux_pacs),'-',color='Green',linewidth=1.5*mag)
	spire, = ax_sed.plot(np.log10(wl_spire),np.log10(c/(wl_spire*1e-4)*flux_spire),'r-',linewidth=1.5*mag)
	irs, = ax_sed.plot(np.log10(wl_irs),np.log10(c/(wl_irs*1e-4)*flux_irs),'-',color='Blue',linewidth=1.5*mag)
	ax_sed.text(0.75,0.9,r'$\mathrm{L_{bol}= %5.2f L_{\odot}}$' % l_bol_obs,fontsize=mag*16,transform=ax_sed.transAxes) 

	# plot the observed photometry data
	photometry, = ax_sed.plot(np.log10(wl_phot),np.log10(c/(wl_phot*1e-4)*flux_phot),'s',mfc='DimGray',mec='k',markersize=8)
	ax_sed.errorbar(np.log10(wl_phot),np.log10(c/(wl_phot*1e-4)*flux_phot),\
		yerr=[np.log10(c/(wl_phot*1e-4)*flux_phot)-np.log10(c/(wl_phot*1e-4)*(flux_phot-flux_sig_phot)),\
			  np.log10(c/(wl_phot*1e-4)*(flux_phot+flux_sig_phot))-np.log10(c/(wl_phot*1e-4)*flux_phot)],\
		fmt='s',mfc='DimGray',mec='k',markersize=8)

	# Extract the SED for the smallest inclination and largest aperture, and
	# scale to 300pc. In Python, negative indices can be used for lists and
	# arrays, and indicate the position from the end. So to get the SED in the
	# largest aperture, we set aperture=-1.
	# aperture group is aranged from smallest to infinite
	sed_inf = m.get_sed(group=0, inclination=0, aperture=-1, distance=dstar * pc)

	# l_bol_sim = l_bol(sed_inf.wav, sed_inf.val/(c/sed_inf.wav*1e4)*1e23)
	# print sed.wav, sed.val
	# print 'Bolometric luminosity of simulated spectrum: %5.2f lsun' % l_bol_sim


	# plot the simulated SED
	sim, = ax_sed.plot(np.log10(sed_inf.wav), np.log10(sed_inf.val), '-', color='GoldenRod', linewidth=1.5*mag)
	# get flux at different apertures
	flux_aper = np.empty_like(wl_aper)
	unc_aper = np.empty_like(wl_aper)
	for i in range(0, len(wl_aper)):
		sed_dum = m.get_sed(group=i+1, inclination=0, aperture=-1, distance=dstar * pc)
		f = interp1d(sed_dum.wav, sed_dum.val)
		flux_aper[i] = f(wl_aper[i])
		# # interpolate the uncertainty (maybe not the best way to do this)
		# print sed_dum.unc
		# f = interp1d(sed_dum.wav, sed_dum.unc)
		# unc_aper[i] = f(wl_aper[i])
		# ax_sed.plot(np.log10(sed_dum.wav), np.log10(sed_dum.val), '-', linewidth=1.5*mag)
		# print l_bol(sed_dum.wav, sed_dum.val/(c/sed_dum.wav*1e4)*1e23)
	aper, = ax_sed.plot(np.log10(wl_aper),np.log10(flux_aper),'o',mfc='None',mec='k',markersize=12,markeredgewidth=3)
	# calculate the bolometric luminosity of the aperture 
	l_bol_sim = l_bol(wl_aper, flux_aper/(c/np.array(wl_aper)*1e4)*1e23)
	print 'Bolometric luminosity of simulated spectrum: %5.2f lsun' % l_bol_sim

	# print out the sed into ascii file for reading in later
	if save == True:
		# unapertured SED
		foo = open(outdir+print_name+'_sed_inf.txt','w')
		foo.write('%12s \t %12s \n' % ('wave','vSv'))
		for i in range(0, len(sed_inf.wav)):
			foo.write('%12g \t %12g \n' % (sed_inf.wav[i], sed_inf.val[i]))
		foo.close()
		# SED with convolution of aperture sizes
		foo = open(outdir+print_name+'_sed_w_aperture.txt','w')
		foo.write('%12s \t %12s \n' % ('wave','vSv'))
		for i in range(0, len(wl_aper)):
			foo.write('%12g \t %12g \n' % (wl_aper[i], flux_aper[i]))
		foo.close()

	# Read in and plot the simulated SED produced by RADMC-3D using the same parameters
	# [wl,fit] = np.genfromtxt(indir+'hyperion/radmc_comparison/spectrum.out',dtype='float',skip_header=3).T
	# l_bol_radmc = l_bol(wl,fit*1e23/dstar**2)
	# radmc, = ax_sed.plot(np.log10(wl),np.log10(c/(wl*1e-4)*fit/dstar**2),'-',color='DimGray', linewidth=1.5*mag, alpha=0.5)

	# print the L bol of the simulated SED (both Hyperion and RADMC-3D)
	# lg_sim = ax_sed.legend([sim,radmc],[r'$\mathrm{L_{bol,sim}=%5.2f~L_{\odot},~L_{center}=9.18~L_{\odot}}$' % l_bol_sim, \
	# 	r'$\mathrm{L_{bol,radmc3d}=%5.2f~L_{\odot},~L_{center}=9.18~L_{\odot}}$' % l_bol_radmc],\
	# 	loc='lower right',fontsize=mag*16)

	# read the input central luminosity by reading in the source information from output file
	dum = Model()
	dum.use_sources(filename)
	L_cen = dum.sources[0].luminosity/lsun

	lg_sim = ax_sed.legend([sim],[r'$\mathrm{L_{bol,sim}=%5.2f~L_{\odot},~L_{center}=%5.2f~L_{\odot}}$' % (l_bol_sim, L_cen)], \
        loc='lower right',fontsize=mag*16)
	# plot setting
	ax_sed.set_xlabel(r'$\mathrm{log~\lambda~({\mu}m)}$',fontsize=mag*20)
	ax_sed.set_ylabel(r'$\mathrm{log~\nu S_{\nu}~(erg/cm^{2}/s)}$',fontsize=mag*20)
	[ax_sed.spines[axis].set_linewidth(1.5*mag) for axis in ['top','bottom','left','right']]
	ax_sed.minorticks_on()
	ax_sed.tick_params('both',labelsize=mag*18,width=1.5*mag,which='major',pad=15,length=5*mag)
	ax_sed.tick_params('both',labelsize=mag*18,width=1.5*mag,which='minor',pad=15,length=2.5*mag)

	ax_sed.set_ylim([-14,-7])
	ax_sed.set_xlim([0,3])

	lg_data = ax_sed.legend([irs, pacs, spire,photometry],[r'$\mathrm{{\it Spitzer}-IRS}$',r'$\mathrm{{\it Herschel}-PACS}$',r'$\mathrm{{\it Herschel}-SPIRE}$',r'$\mathrm{Photometry}$'],\
							loc='best',fontsize=14*mag,numpoints=1,framealpha=0.3)
	plt.gca().add_artist(lg_sim)

	# Write out the plot
	fig.savefig(outdir+print_name+'_sed.pdf',format='pdf',dpi=300,bbox_inches='tight')
	fig.clf()

	# Package for matching the colorbar
	from mpl_toolkits.axes_grid1 import make_axes_locatable

	# Extract the image for the first inclination, and scale to 300pc. We
	# have to specify group=1 as there is no image in group 0.
	image = m.get_image(group=len(wl_aper)+1, inclination=0, distance=dstar * pc, units='MJy/sr')
	# image = m.get_image(group=14, inclination=0, distance=dstar * pc, units='MJy/sr')
	# Open figure and create axes
	# fig = plt.figure(figsize=(8, 8))
	fig, axarr = plt.subplots(3, 3, sharex='col', sharey='row',figsize=(13.5,12))

	# Pre-set maximum for colorscales
	VMAX = {}
	# VMAX[3.6] = 10.
	# VMAX[24] = 100.
	# VMAX[160] = 2000.
	# VMAX[500] = 2000.
	VMAX[100] = 10.
	VMAX[250] = 100.
	VMAX[500] = 2000.
	VMAX[1000] = 2000.

	# We will now show four sub-plots, each one for a different wavelength
	# for i, wav in enumerate([3.6, 24, 160, 500]):
	# for i, wav in enumerate([100, 250, 500, 1000]):
	# for i, wav in enumerate([4.5, 9.7, 24, 40, 70, 100, 250, 500, 1000]):
	for i, wav in enumerate([250, 500, 1000, 250, 500, 1000, 250, 500, 1000]):


		# ax = fig.add_subplot(3, 3, i + 1)
		ax = axarr[i/3, i%3]

		# Find the closest wavelength
		iwav = np.argmin(np.abs(wav - image.wav))

		# Calculate the image width in arcseconds given the distance used above
		w = np.degrees((1.5 * pc) / image.distance) * 60.

		# Image in the unit of MJy/sr
		# Change it into erg/s/cm2/Hz/sr
		factor = 1e-23*1e6
		# avoid zero in log
		val = image.val[:, :, iwav] * factor + 1e-30

		# This is the command to show the image. The parameters vmin and vmax are
		# the min and max levels for the colorscale (remove for default values).
		im = ax.imshow(np.log10(val), vmin= -22, vmax= -12,
				  cmap=plt.cm.jet, origin='lower', extent=[-w, w, -w, w], aspect=1)

		# Colorbar setting
		# create an axes on the right side of ax. The width of cax will be 5%
		# of ax and the padding between cax and ax will be fixed at 0.05 inch.
		if (i+1) % 3 == 0:
			divider = make_axes_locatable(ax)
			cax = divider.append_axes("right", size="5%", pad=0.05)
			cb = fig.colorbar(im, cax=cax)
			cb.solids.set_edgecolor("face")
			cb.ax.minorticks_on()
			cb.ax.set_ylabel(r'$\mathrm{log(I_{\nu})~[erg/s/cm^{2}/Hz/sr]}$',fontsize=12)
			cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
			plt.setp(cb_obj,fontsize=12)

		if (i+1) == 7:
			# Finalize the plot
			ax.set_xlabel('RA Offset (arcsec)')
			ax.set_ylabel('Dec Offset (arcsec)')

		ax.tick_params(axis='both', which='major', labelsize=16)
		ax.set_adjustable('box-forced')
		ax.text(0.7,0.88,str(wav) + r'$\mathrm{~\mu m}$',fontsize=16,color='white', transform=ax.transAxes)

	fig.subplots_adjust(hspace=0,wspace=-0.2)

	# Adjust the spaces between the subplots 
	# plt.tight_layout()
	fig.savefig(outdir+print_name+'_cube_plot.png', format='png', dpi=300, bbox_inches='tight')
	fig.clf()

# indir = '/Users/yaolun/bhr71/obs_for_radmc/'
# outdir = '/Users/yaolun/bhr71/hyperion/'
# extract_hyperion('/Users/yaolun/test/model125/model125.rtout',indir=indir,outdir='/Users/yaolun/test/')
# extract_hyperion('/hyperion/best_model_bettyjo.rtout',indir=indir,outdir=outdir+'bettyjo/')
# extract_hyperion('/hyperion/old_setup2.rtout',indir=indir,outdir=outdir)
