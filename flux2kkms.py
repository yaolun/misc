def flux2kkms(wl, flux, pix_size):
	"""
	wl: in um
	flux: in erg s-1 cm-2
	pix_size: in arcsec
	"""
	import numpy as np
	import astropy.constants as const

	k = const.k_B.cgs.value

	omega_b = np.pi/4/np.log(2)*pix_size**2*(1/3600.0/180*np.pi)**2		# in sr
	kkms = (wl*1e-4)**3/2/k * flux / omega_b / 1e5

	return kkms

def coord_offset(center, offset, unit='deg'):
	"""
	center = (RA, Dec) tuple in degree
	offset = (delta_RA, delta_Dec) tuple in arcsec
	"""
	import numpy as np

	RA_cen, Dec_cen = center
	d_RA, d_Dec = offset

	RA_off = d_RA / 3600. / np.cos(np.radians(Dec_cen)) + RA_cen
	Dec_off = d_Dec / 3600. + Dec_cen

	return RA_off, Dec_off

import numpy as np
import astropy.constants as const
c = const.c.cgs.value
pix_size = 9.4
assumed_velo = 1.

# [OI] 63 um
print 'L1551-IRS5 1-D'
wl = 63.174267
flux = 5.3951846e-19 * 1e7
print flux2kkms(wl, flux, pix_size) / assumed_velo

print 'L1551-IRS5 cube at center'
wl = 63.174233
flux = 3.3077937e-19 * 1e7
print flux2kkms(wl, flux, pix_size) / assumed_velo

print 'estimate red-shifted peak'
# offset to the center: (d_RA, d_Dec) = (6, 4).  velocity offset ~ +40 km/s.
# together with v_lsr = 6.5 km/s and Earth v_lsr ~ -10 km/s, the line will center around 32.5 km/s
center = (67.89196135, 18.13468391)
offset = (6,4)
print coord_offset(center, offset)
print 'estimated flux at the coordinates above'
rest_wl = 63.18367004
wl = rest_wl + (40.)/c*rest_wl
flux = 2.4e-14*(9.4/2)**2*np.pi
print flux2kkms(wl, flux, pix_size) / assumed_velo

print 'estimate blue-shifted peak'
# offset to the center: (d_RA, d_Dec) = (-6, 0), velocity offset ~ -80 km/s.
# together with v_lsr = 6.5 km/s (assume Earth v_lsr within +/- 20 km/s), the line will center around -85 km/s +/- 20 km/s
offset = (-6, 0)
print coord_offset(center, offset)
print 'estimated flux at the coordinates above'
wl = rest_wl + (-80.)/c*rest_wl
flux = 3.6e-14*(9.4/2)**2*np.pi
print flux2kkms(wl, flux, pix_size) / assumed_velo
