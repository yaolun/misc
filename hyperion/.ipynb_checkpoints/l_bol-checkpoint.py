def l_bol(wl,fv,dist):
	import numpy as np
	import astropy.constants as const
	import astropy.units as uni
	import matplotlib.pyplot as plt
	# wavelength unit: um
	# Flux density unit: Jy
	# distance unit: pc
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
