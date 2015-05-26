def get_bhr71_obs(indir):
	import numpy as np
	from spitzer_unc import spitzer_unc

	# Read in Herschel data
	# continuum
	[wl_pacs,flux_pacs,unc_pacs] = np.genfromtxt(indir+'BHR71_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_continuum.txt',dtype='float',skip_header=1).T
	[wl_spire,flux_spire] = np.genfromtxt(indir+'BHR71_spire_corrected_continuum.txt',dtype='float',skip_header=1).T
	# original spectra
	[wl_pacs_data,flux_pacs_data,unc_pacs_data] = np.genfromtxt(indir+'BHR71_centralSpaxel_PointSourceCorrected_CorrectedYES_trim.txt',dtype='float').T
	[wl_spire_data,flux_spire_data] = np.genfromtxt(indir+'BHR71_spire_corrected.txt',dtype='float').T
	# flat spectra
	[wl_pacs_flat,flux_pacs_flat,unc_pacs_flat] = np.genfromtxt(indir+'BHR71_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_flat_spectrum.txt',dtype='float',skip_header=1).T
	[wl_spire_flat,flux_spire_flat] = np.genfromtxt(indir+'BHR71_spire_corrected_flat_spectrum.txt',dtype='float',skip_header=1).T

	flux_pacs_noise = flux_pacs_data-flux_pacs-flux_pacs_flat
	flux_spire_noise = flux_spire_data-flux_spire-flux_spire_flat

	# Read in the Spitzer IRS spectrum
	# [wl_irs, flux_irs]= (np.genfromtxt(indir+'bhr71_spitzer_irs.txt',skip_header=2,dtype='float').T)[0:2]
	# # Remove points with zero or negative flux 
	# ind = flux_irs > 0
	# wl_irs = wl_irs[ind]
	# flux_irs = flux_irs[ind]
	wl_irs, flux_irs, unc_irs = spitzer_unc(indir+'bhr71_spitzer_irs.txt')
	# Calculate the local variance (for spire), use the instrument uncertainty for pacs
	#
	# Spitzer noise is not considered now
	wl_noise = [wl_pacs_data[wl_pacs_data <= 190.31],wl_spire[(wl_spire > 194) & (wl_spire <= 304)],wl_spire[wl_spire > 304]]
	flux_noise = [unc_pacs[wl_pacs_data <= 190.31],flux_spire_noise[(wl_spire > 194) & (wl_spire <= 304)],flux_spire_noise[wl_spire > 304]]
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
	sigma_noise = np.hstack((unc_irs, sigma_noise))
	# print len(wl_spec), len(sigma_noise)
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

	# Print the observed L_bol
	wl_spec = np.hstack((wl_irs,wl_pacs,wl_spire))
	flux_spec = np.hstack((flux_irs,flux_pacs,flux_spire))
	flux_spec = flux_spec[np.argsort(wl_spec)]
	wl_spec = wl_spec[np.argsort(wl_spec)]

	return {'spec': (wl_spec, flux_spec, sigma_noise), 'phot': (wl_phot, flux_phot, flux_sig_phot)}
