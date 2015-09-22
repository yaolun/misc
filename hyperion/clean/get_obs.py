def get_obs(indir, obj='BHR71'):
	"""
	obj input in uppercase.  But check the path to make sure.
	"""
	import numpy as np
	from spitzer_unc import spitzer_unc

	# Read in Herschel data
	# continuum
	[wl_pacs,flux_pacs,unc_pacs] = np.genfromtxt(indir+obj+'_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_continuum.txt',dtype='float',skip_header=1).T
	[wl_spire,flux_spire] = np.genfromtxt(indir+obj+'_spire_corrected_continuum.txt',dtype='float',skip_header=1).T
	# noise spectra
	[wl_pacs_noise, flux_pacs_noise,unc_pacs_noise] = np.genfromtxt(indir+obj+'_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_noise_spectrum.txt',dtype='float',skip_header=1).T
	[wl_spire_noise,flux_spire_noise] = np.genfromtxt(indir+obj+'_spire_corrected_noise_spectrum.txt',dtype='float',skip_header=1).T

	# Read in Spitzer data
	wl_irs, flux_irs, unc_irs = spitzer_unc(indir+obj.lower()+'_spitzer_irs.txt')
	
	# Calculate the local variance (for spire), use the instrument uncertainty for pacs
	#
	wl_noise = [wl_pacs_noise, wl_spire_noise]
	flux_noise = [flux_pacs_noise, flux_spire_noise]
	sig_num = 20
	sigma_noise = []
	for i in range(0,len(wl_noise)):
		sigma_dum = np.zeros_like(wl_noise[i])
		for iwl in range(0,len(wl_noise[i])):
			if iwl < sig_num/2:
				sigma_dum[iwl] = np.std(np.hstack((flux_noise[i][0:sig_num/2],flux_noise[i][0:sig_num/2-iwl])))
			elif len(wl_noise[i])-iwl < sig_num/2:
				sigma_dum[iwl] = np.std(np.hstack((flux_noise[i][iwl:],flux_noise[i][len(wl_noise[i])-sig_num/2:])))
			else:
				sigma_dum[iwl] = np.std(flux_noise[i][iwl-sig_num/2:iwl+sig_num/2])
		sigma_noise = np.hstack((sigma_noise, sigma_dum))
	sigma_noise = np.hstack((unc_irs, sigma_noise))

	# Read in the photometry data
	phot = np.genfromtxt(indir+obj.lower()+'.txt',dtype=None,skip_header=1,comments='%')
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
	flux_phot = np.array(flux_phot)
	flux_sig_phot = np.array(flux_sig_phot)

	# Print the observed L_bol
	wl_spec = np.hstack((wl_irs,wl_pacs,wl_spire))
	flux_spec = np.hstack((flux_irs,flux_pacs,flux_spire))
	flux_spec = flux_spec[np.argsort(wl_spec)]
	wl_spec = wl_spec[np.argsort(wl_spec)]

	return {'spec': (wl_spec, flux_spec, sigma_noise), 'phot': (wl_phot, flux_phot, flux_sig_phot)}