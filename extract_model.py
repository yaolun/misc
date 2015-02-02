def extract_hyperion(filename,indir=None,outdir=None,dstar=178.0):
	def l_bol(wl,fv,dist=178.0):
		import numpy as np
		import astropy.constants as const
		import astropy.units as uni
		import matplotlib.pyplot as plt
		# wavelength unit: um
		# Flux density unit: Jy
		#
		# constants setup
		#
		c = const.c.cgs.value
		pc = 3.086e+18
		PI = np.pi
		SL = 3.846e+33
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

	from hyperion.model import ModelOutput
	from hyperion.util.constants import pc, c

	# Read in the observation data and calculate the noise & variance
	if indir == None:
		indir = '/Users/yaolun/bhr71/'
	if outdir == None:
		outdir = '/Users/yaolun/bhr71/hyperion/'
	[wl_pacs,flux_pacs,unc_pacs] = np.genfromtxt(indir+'obs_for_radmc/BHR71_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_continuum.txt',\
										dtype='float',skip_header=1).T
	# Convert the unit from Jy to erg cm-2 Hz-1
	flux_pacs = flux_pacs*1e-23
	[wl_spire,flux_spire] = np.genfromtxt(indir+'obs_for_radmc/BHR71_spire_corrected_continuum.txt',dtype='float',skip_header=1).T
	flux_spire = flux_spire*1e-23 
	wl_obs = np.hstack((wl_pacs,wl_spire))
	flux_obs = np.hstack((flux_pacs,flux_spire))

	[wl_pacs_data,flux_pacs_data,unc_pacs_data] = np.genfromtxt(indir+'obs_for_radmc/BHR71_centralSpaxel_PointSourceCorrected_CorrectedYES_trim.txt',\
												  dtype='float').T
	[wl_spire_data,flux_spire_data] = np.genfromtxt(indir+'obs_for_radmc/BHR71_spire_corrected.txt',\
													dtype='float').T

	[wl_pacs_flat,flux_pacs_flat,unc_pacs_flat] = np.genfromtxt(indir+'obs_for_radmc/BHR71_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_flat_spectrum.txt',\
										dtype='float',skip_header=1).T
	[wl_spire_flat,flux_spire_flat] = np.genfromtxt(indir+'obs_for_radmc/BHR71_spire_corrected_flat_spectrum.txt',dtype='float',skip_header=1).T

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
	[wl_irs, flux_irs]= (np.genfromtxt(indir+'obs_for_radmc/bhr71_spitzer_irs.txt',skip_header=2,dtype='float').T)[0:2]
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
	phot = np.genfromtxt(indir+'obs_for_radmc/bhr71.txt',dtype=None,skip_header=1,comments='%')
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
	m = ModelOutput(indir+filename)

	# Create the plot
	mag = 1.5
	fig = plt.figure(figsize=(8*mag,6*mag))
	ax_sed = fig.add_subplot(1, 1, 1)

	# Extract the SED for the smallest inclination and largest aperture, and
	# scale to 300pc. In Python, negative indices can be used for lists and
	# arrays, and indicate the position from the end. So to get the SED in the
	# largest aperture, we set aperture=-1.
	sed = m.get_sed(inclination=0, aperture=-1, distance=dstar * pc)
	l_bol_sim = l_bol(sed.wav, sed.val/(c/sed.wav*1e4)*1e23)
	# print sed.wav, sed.val
	print 'Bolometric luminosity of simulated spectrum: %5.2f' % l_bol_sim

	# Plot the observed SED
	# plot the observed spectra
	pacs, = ax_sed.plot(np.log10(wl_pacs),np.log10(c/(wl_pacs*1e-4)*flux_pacs),'-',color='Green',linewidth=1.5*mag)
	spire, = ax_sed.plot(np.log10(wl_spire),np.log10(c/(wl_spire*1e-4)*flux_spire),'r-',linewidth=1.5*mag)
	irs, = ax_sed.plot(np.log10(wl_irs),np.log10(c/(wl_irs*1e-4)*flux_irs),'-',color='Blue',linewidth=1.5*mag)
	ax_sed.text(0.75,0.9,r'$\mathrm{L_{bol}= %5.2f L_{\odot}}$' % l_bol_obs,fontsize=mag*16,transform=ax_sed.transAxes) 

	# plot the observed photometry data
	photometry, = ax_sed.plot(np.log10(wl_phot),np.log10(c/(wl_phot*1e-4)*flux_phot),'s',mfc='DimGray',mec='k',markersize=8)
	ax_sed.errorbar(np.log10(wl_phot),np.log10(c/(wl_phot*1e-4)*flux_phot),yerr=flux_sig_phot,fmt='s',mfc='DimGray',mec='k',markersize=8)

	# plot the simulated SED
	sim, = ax_sed.plot(np.log10(sed.wav), np.log10(sed.val), '-', color='GoldenRod', linewidth=1.5*mag)

	# Read in and plot the simulated SED produced by RADMC-3D using the same parameters
	[wl,fit] = np.genfromtxt(indir+'hyperion/radmc_comparison/spectrum.out',dtype='float',skip_header=3).T
	l_bol_radmc = l_bol(wl,fit*1e23/dstar**2)
	radmc, = ax_sed.plot(np.log10(wl),np.log10(c/(wl*1e-4)*fit/dstar**2),'-',color='DimGray', linewidth=1.5*mag, alpha=0.5)

	# print the L bol of the simulated SED (both Hyperion and RADMC-3D)
	lg_sim = ax_sed.legend([sim,radmc],[r'$\mathrm{L_{bol,sim}=%5.2f~L_{\odot},~L_{center}=9.18~L_{\odot}}$' % l_bol_sim, \
		r'$\mathrm{L_{bol,radmc3d}=%5.2f~L_{\odot},~L_{center}=9.18~L_{\odot}}$' % l_bol_radmc],\
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
	fig.savefig(outdir+'best_model_sed.pdf',format='pdf',dpi=300,bbox_inches='tight')
	fig.clf()

	# Package for matching the colorbar
	from mpl_toolkits.axes_grid1 import make_axes_locatable

	# Extract the image for the first inclination, and scale to 300pc. We
	# have to specify group=1 as there is no image in group 0.
	image = m.get_image(inclination=0, distance=dstar * pc, units='MJy/sr')

	# Open figure and create axes
	fig = plt.figure(figsize=(8, 8))

	# Pre-set maximum for colorscales
	VMAX = {}
	VMAX[3.6] = 10.
	VMAX[24] = 100.
	VMAX[160] = 2000.
	VMAX[500] = 2000.

	# We will now show four sub-plots, each one for a different wavelength
	for i, wav in enumerate([3.6, 24, 160, 500]):

		ax = fig.add_subplot(2, 2, i + 1)

		# Find the closest wavelength
		iwav = np.argmin(np.abs(wav - image.wav))

		# Calculate the image width in arcseconds given the distance used above
		w = np.degrees((1.5 * pc) / image.distance) * 60.

		# Image in the unit of MJy/sr
		# Change it into erg/s/cm2/Hz/sr
		factor = 1e-23*1e6
		# avoid zero in log
		image.val[:, :, iwav] = image.val[:, :, iwav] * factor + 1e-30

		# This is the command to show the image. The parameters vmin and vmax are
		# the min and max levels for the colorscale (remove for default values).
		im = ax.imshow(np.log10(image.val[:, :, iwav]).T, vmin= -22, vmax= -12,
				  cmap=plt.cm.jet, origin='lower', extent=[-w, w, -w, w])

		# Colorbar setting
		# create an axes on the right side of ax. The width of cax will be 5%
		# of ax and the padding between cax and ax will be fixed at 0.05 inch.
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		cb = fig.colorbar(im, cax=cax)
		cb.solids.set_edgecolor("face")
		cb.ax.minorticks_on()
		cb.ax.set_ylabel(r'$\mathrm{log(I_{\nu})~[erg/s/cm^{2}/Hz/sr]}$',fontsize=12)
		cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
		plt.setp(cb_obj,fontsize=12)
		
		# Finalize the plot
		ax.tick_params(axis='both', which='major', labelsize=10)
		ax.set_xlabel('x (arcsec)')
		ax.set_ylabel('y (arcsec)')
		ax.text(0.5,0.88,str(wav) + r'$\mathrm{~\mu m}$',fontsize=16,color='white', transform=ax.transAxes)
	# Adjust the spaces between the subplots 
	plt.tight_layout()
	fig.savefig(outdir+'simple_cube_plot.pdf', format='pdf', dpi=300, bbox_inches='tight') 

indir = '/Users/yaolun/bhr71/'
outdir = '/Users/yaolun/bhr71/hyperion/'
# extract_hyperion('/hyperion/best_model.rtout',indir=indir)
# extract_hyperion('/hyperion/best_model_bettyjo.rtout',indir=indir,outdir=outdir+'bettyjo/')
extract_hyperion('/hyperion/best_model.rtout',indir=indir,outdir=outdir+'bettyjo/newton_')
