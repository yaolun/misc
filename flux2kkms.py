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
from astropy.coordinates import SkyCoord
from astropy import units as u
c = const.c.cgs.value
pix_size = 9.4
assumed_velo = 1.

# [OI] 63 um & [CII] 158 um
print 'L1551-IRS5 1-D'
wl = np.array([63.174267, 157.76002])
flux = np.array([5.3951846e-19, 1.3200978e-20]) * 1e7
print flux2kkms(wl, flux, pix_size) / assumed_velo

print 'L1551-IRS5 cube at center'
wl = np.array([63.174233, 157.75969])
flux = np.array([3.3077937e-19, 4.3660825e-21]) * 1e7
print flux2kkms(wl, flux, pix_size) / assumed_velo

# offset to the center: (d_RA, d_Dec) = (6, 4).  velocity offset ~ +40 km/s.
# together with v_lsr = 6.5 km/s and Earth v_lsr ~ -10 km/s, the line will center around 32.5 km/s
center = (67.89196135, 18.13468391)
# coordinates transformation
cen = SkyCoord(ra=center[0]*u.degree, dec=center[1]*u.degree, frame='icrs')
print cen.ra.hms, cen.dec.dms

print 'estimate red-shifted peak'
offset = (6,4)

blue_coord = coord_offset(center, offset)
blue = SkyCoord(ra=blue_coord[0]*u.degree, dec=blue_coord[1]*u.degree, frame='icrs')
print blue.ra.hms, blue.dec.dms

print 'estimated flux at the coordinates above'
rest_oi_wl = 63.18367004
rest_cii_wl = 157.6922760
print (63.174233-rest_oi_wl)/rest_oi_wl*c/1e5
wl = np.array([rest_oi_wl + (50.)/c*rest_oi_wl, rest_cii_wl + (150.)/c*rest_cii_wl])
flux = np.array([3.6e-14, 8e-16])*(9.4/2)**2*np.pi
print flux2kkms(wl, flux, pix_size) / assumed_velo

print 'estimate blue-shifted peak'
# offset to the center: (d_RA, d_Dec) = (-6, 0), velocity offset ~ -80 km/s.
# together with v_lsr = 6.5 km/s (assume Earth v_lsr within +/- 20 km/s), the line will center around -85 km/s +/- 20 km/s
offset = (-6, -4)

red_coord = coord_offset(center, offset)
red = SkyCoord(ra=red_coord[0]*u.degree, dec=red_coord[1]*u.degree, frame='icrs')
print red.ra.hms, red.dec.dms

print 'estimated flux at the coordinates above'
wl = np.array([rest_oi_wl + (-90.)/c*rest_oi_wl, rest_cii_wl + (100.)/c*rest_cii_wl])
flux = np.array([3.6e-14, 6e-16])*(9.4/2)**2*np.pi
print flux2kkms(wl, flux, pix_size) / assumed_velo

print 'off-source coordinates'
offset = (-22.18, 33.28)
off_coord = coord_offset(center, offset)
off = SkyCoord(ra=off_coord[0]*u.degree, dec=off_coord[1]*u.degree, frame='icrs')
print off.ra.hms, off.dec.dms
