def t_bol(wave, flux, freq=False):
	"""
	wave: wavelength in um unless freq is set to True.  If freq is True, then wave is frequency in Hz
	flux: flux in Jy as always.
	"""
	import numpy as np
	import scipy.special as spec
	import astropy.constants as const

	# constants setup
	c = const.c.cgs.value
	h = const.h.cgs.value
	k = const.k_B.cgs.value

	# convert unit from (um, Jy) -> (Hz, erg s-1 cm-2 Hz-1)
	fv = np.array(flux) * 1e-23
	if freq == False:
		freq = c/(1e-4*np.array(wave))
	else:
		freq = wave

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

	# calculate the zeroth and first moment
	I1 = np.trapz(fv*freq, freq)
	I0 = np.trapz(fv, freq)

	# T_bol equation from Myers & Ladd 1993
	t_bol = (spec.zetac(4)+1)/(4*(spec.zetac(5)+1))*h/k * (I1/I0)

	return t_bol


# import numpy as np
# wave, flux, unc = np.genfromtxt('/Users/yaolun/bhr71/fitting/latest/pacs/data/BHR71_centralSpaxel_PointSourceCorrected_CorrectedYES_trim.txt',skip_header=1).T
# print t_bol(wave, flux)
