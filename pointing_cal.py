import numpy as np
from astropy.coordinates import ICRS
from astropy import units as u
# Pointing Position Calculation

# Pixel		RA			Dec			Velocity
# Pixel 13 	67.891847	18.134617	-44.502739
# Pxiel  7	67.889155	18.133292	 83.823706

ra13  = 67.891847
dec13 = 18.134617
ra7   = 67.889155
dec7  = 18.133292
v13   = -44.502739
v7    = 83.823706
# Convert to projected coordinate ( RA = RA * cos(Dec) )

ra13  = ra13*np.cos(dec13*np.pi/180.)
ra7   = ra7 *np.cos(dec7 *np.pi/180.)

# Pixel		RA			Dec
# Pixel 13 	61.3146902	18.134617
# Pixel  7  61.3131884	18.133292

# Interpolation

ra  = ra13 + (ra7 - ra13)  * (0-v13)/(v7-v13)
dec = dec13+ (dec7- dec13) * (0-v13)/(v7-v13)

# Calculate the off-source position
#
# Outflow direction: 50 degree respect to the North
# Off by 100 arcsec distance

ra_off1  = ra  + 100./3600*np.cos(50*np.pi/180.)
dec_off1 = dec - 100./3600*np.sin(50*np.pi/180.)

ra_off2  = ra  - 100./3600*np.cos(50*np.pi/180.)
dec_off2 = dec  + 100./3600*np.sin(50*np.pi/180.)

# Convert them back to observation coordinates

ra      = ra / np.cos(dec*np.pi/180.)
ra_off1 = ra_off1 / np.cos(dec_off1*np.pi/180.)
ra_off2 = ra_off2 / np.cos(dec_off2*np.pi/180.)

# Convert from degree to hh:mm:ss dd:mm:ss

on   = ICRS(ra=ra,dec=dec,unit=(u.degree,u.degree))
off1 = ICRS(ra=ra_off1,dec=dec_off1,unit=(u.degree,u.degree))
off2 = ICRS(ra=ra_off2,dec=dec_off2,unit=(u.degree,u.degree))
JL   = ICRS('4h31m34.1s +18d08m4.9s')
print JL

print 'On   > RA: %s, Dec: %s' % (on.ra.to_string(u.hour),on.dec.to_string(u.degree))
print 'Off1 > RA: %s, Dec: %s' % (off1.ra.to_string(u.hour),off1.dec.to_string(u.degree))
print 'Off2 > RA: %s, Dec: %s' % (off2.ra.to_string(u.hour),off2.dec.to_string(u.degree))

