##############################################
#                                            #
#  Short script to convert CASA viewer fits  #
#    output to CLASS readable fits file.     #
#    for data from CASA 3.4, and CLASS       #
#               June 2014                    #
#                                            #
#              Magnus Persson                #
#        magnusp@strw.leidenuniv.nl          #
##############################################
"""
Change filenames after you preferences, then
import and view the new FITS file/spectrum in CLASS:

fits read spw0_newhdr.fits /check
go browse 

One key difference is that CLASS uses the RESTFREQ
keyword together with the CRVAL1 keyword to create the
frequency array.
"""
# YOUR definitions
filename = '/Volumes/SD-Mac/Google Drive/research/bhr71_infall/analysis/cs-line-0p2.fits'
outfile = '/Volumes/SD-Mac/Google Drive/research/bhr71_infall/analysis/cs-line-0p2-converted.fits'
vlsr = 8.3e3  # in m/s
objstr = 'B335'
unit = 'Jy'
#############################################
# imports
try:
    from astropy.io import fits
except:
    try:
        import pyfits as fits
    except:
        raise(ImportError('Need either AstroPy or PyFits installed'))
#~ import astropy.wcs
from scipy import array, zeros, sign
#############################################
# other constant(s)
# speed of light in m/s (from cm/s)!
CC = 29979245800.0e-2 # m/s
#############################################
# Keywords to DELET
delwords = ['SPECSYS', 'ALTRVAL' , 'ALTRPIX' , 'VELREF']
#############################################
# get the header and data to work with
hdulist = fits.open(filename)
# and from that, the primary header
hdr = hdulist[0].header
data = hdulist[0].data
hdulist.close()
#############################################
# calculate the RA and DEC position in degrees
pos_str = hdr['POSITION']
ra_str, dec_str = pos_str[:12], pos_str[12:]
# now as e.g. '16:34:23.66' and '24d20m30.0s'
# we need to extract all of it and convert
# to degrees for both RA and Dec
# make translation table, i.e. what to
# to replace and with what
# T =  maketrans('md:', ':::')                  
# ra = array([float(i) for i in translate(ra_str, T).split(':')])
# dec = array([float(i) for i in translate(dec_str, T).split(':')])
# modify by YLY
T =  str.maketrans('md:', ':::')
ra = array([float(i) for i in ra_str.translate(T).split(':')])
dec = array([float(i) for i in dec_str.translate(T).split(':')])

ra_tran = array([1,sign(ra[0])*60,sign(ra[0])*3600])
dec_tran = array([1,sign(dec[0])*60,sign(dec[0])*3600])
# sum to get the float values
ra_deg = (ra/ra_tran).sum()
# multiply with 15 to get degrees
ra_deg *= 15.0                                
dec_deg = (dec/dec_tran).sum()
# BOOM, coordinates done, not really necessary though.
#############################################
# extended data with 2 axis
tmp_data = zeros((1,1,len(data)))
tmp_data[0,0] = data.copy()
data = tmp_data
#############################################
# calculate the channel width in m/s
fdelt = hdr['CDELT1']
fref = hdr['RESTFRQ']
# add RESTFREQ if it doesn't exist
if 'RESTFREQ' not in hdr.keys():
    print('No RESTFREQ keyword, creating')
    hdr.set('RESTFREQ', fref)
# doppler formula to get the delta-V
deltav = fdelt * CC / fref
print('DELTFREQ : ' + str(fdelt))
print('DELTAV : ' + str(deltav))
#############################################
# define keywords to change/add in a list
#kw = ['RESTFREQ',      'VELO-LSR', 'DELTAV', 'OBJECT', 'BUNIT']
kw = ['VELO-LSR', 'DELTAV', 'OBJECT', 'BUNIT']
#val = [hdr['RESTFRQ '],   vlsr,     deltav,  objstr,    unit  ]
val = [ vlsr,     deltav,  objstr,    unit  ]
# for RA
kw.append('CTYPE2')
val.append('RA---SIN')
kw.append('CRVAL2')
val.append(ra_deg)
kw.append('CDELT2')
val.append(0.0)
kw.append('CRPIX2')
val.append(0.0)
# for Dec
kw.append('CTYPE3')
val.append('DEC--SIN')
kw.append('CRVAL3')
val.append(dec_deg)
kw.append('CDELT3')
val.append(0.0)
kw.append('CRPIX3')
val.append(0.0)

#############################################
# do the change
for keyword,value in zip(kw, val):
    hdr.set(str(keyword), value=value)
for keyword in delwords:
    hdr.remove(keyword)
print('Min: {0:.4f}, Max: {1:.4f}'.format(data.min(), data.max()))
hdr.set('DATAMIN', data.min())
hdr.set('DATAMAX', data.max())

# >>>> UGLY^3 "FIX" <<<<
# first copy the CRVAL1 and RESTFREQ keywords to backup keywords
hdr.set('CR1OLD', hdr['CRVAL1'], 'BACKUP : old CRVAL1 value')
hdr.set('RFQOLD', hdr['RESTFREQ'], 'BACKUP : old RESTFREQ value')
# now define RESTFREQ to CRVAL1
# and CRVAL1 to 0.0
# Add one CDELT1, because? I dunno, it seemed to make it work,
# at least for this data (Tim's ALMA band 9 data on IRAS 16293+2422)
hdr['RESTFREQ'] = hdr['CRVAL1'] + hdr['CDELT1']
#~ hdr.remove('RESTFREQ') # test
hdr['CRVAL1'] = 0.0

fits.writeto(outfile, data.astype('float32'), header=hdr, output_verify="fix", overwrite=True)
