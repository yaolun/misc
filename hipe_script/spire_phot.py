
##############################################################################
# Example for performing Aperture Photometry for Extended Emission within HIPE
##############################################################################

############################### Setup #####################################

from herschel.share.unit import Constant
c  = Constant.SPEED_OF_LIGHT.value
from java.lang.Math import PI

############################### Setup #####################################


# Loading an observation (Gamma Draconis) from the HSA
obsidFTS = 1342248249 # main obsid
obs = getObservation(obsid=obsidFTS, useHsa=True)
spec = obs.refs['level2'].product.refs['HR_spectrum_point_apod'].product
# spec = fitsReader(file = '/Users/yaolun/bhr71/calibration_testing/data/1342248249/level2/HR_spectrum_point_apod/hspirespectrometer1342248249_a1060001_spgApod_HR_20spss_1454323899152.fits')

outdir = '/Users/yaolun/bhr71/best_calibrated/'
outdir = '/Users/yaolun/test/B335/'

ra  = spec.meta["ra"].value
dec = spec.meta["dec"].value
# For a source with spectrum S(nu) proportional to S^alpha
# SPIRE default pipeline assumes alpha = -1
# For BHR71
alpha = [2.04, 2.35, 2.59]

# beta =1.0
# tempK =16.0

# Run for an individual SPIRE band
array = ["PSW", "PMW","PLW"]  # SPIRE Array Bands: "PSW", "PMW", "PLW"

############################ Import Data ##################################
# Loading an observation of Gamma Dra from the HSA
obsList = getPhotObsidsForFts(spec.obsid)
obsid = obsList[0]
obs     = getObservation(obsid, useHsa=True, instrument='SPIRE')
# For observation from your own local pool use the following line instead
# obs     = getObservation(obsid, poolName='mypool', instrument='SPIRE')


# Extract the Point Source (Jy/beam) and Extended (MJy/sr) calibrated maps
# from the Observation Context
mapPsrcPSW = obs.level2.refs["psrc"+array[0]].product
mapExtdPSW = obs.level2.refs["extd"+array[0]].product
mapPsrcPMW = obs.level2.refs["psrc"+array[1]].product
mapExtdPMW = obs.level2.refs["extd"+array[1]].product
mapPsrcPLW = obs.level2.refs["psrc"+array[2]].product
mapExtdPLW = obs.level2.refs["extd"+array[2]].product

######################## Correction Parameters ############################
# Values are obtained from the SPIRE calibration tree assuming
# a point source with the spectral index alpha specified above
# N.B.: spire_cal_12_02 or later version is needed
# Beam Area for pipeline (alpha=-1)
# Beam Corrections: Beam(alpha=-1)/Beam(alpha)
# Aperture Corrections
# Colour Corrections for point sources (pipeline assumes alpha=-1)


beamCorrTable  = cal.phot.refs["ColorCorrBeam"].product
aperCorrTable  = cal.phot.colorCorrApertureList.refs[0].product
kCorrPsrcTable = cal.phot.colorCorrKList.refs[0].product

beamAreaPSW  = beamCorrTable.meta["beamPipeline"+array[0].title()+"Arc"].double
beamCorrPSW  = beamCorrTable.getAlphaCorrection(alpha[0], array[0])
kCorrPsrcPSW = kCorrPsrcTable.getAlphaCorrection(alpha[0], array[0])
aperCorrPSW  = aperCorrTable.getApertColorCorrection(alpha[0], array[0])

beamAreaPMW  = beamCorrTable.meta["beamPipeline"+array[1].title()+"Arc"].double
beamCorrPMW  = beamCorrTable.getAlphaCorrection(alpha[1], array[1])
kCorrPsrcPMW = kCorrPsrcTable.getAlphaCorrection(alpha[1], array[1])
aperCorrPMW  = aperCorrTable.getApertColorCorrection(alpha[1], array[1])

beamAreaPLW  = beamCorrTable.meta["beamPipeline"+array[2].title()+"Arc"].double
beamCorrPLW  = beamCorrTable.getAlphaCorrection(alpha[2], array[2])
kCorrPsrcPLW = kCorrPsrcTable.getAlphaCorrection(alpha[2], array[2])
aperCorrPLW  = aperCorrTable.getApertColorCorrection(alpha[2], array[2])

print "beamCorrPSW = ",beamCorrPSW
print "beamCorrPMW = ",beamCorrPMW
print "beamCorrPLW = ",beamCorrPLW

# K4p is the conversation factor (colour correction) applied in the standard
# pipeline to give flux for nu S_nu= constant point source
# (converts from alpha=0 to alpha=-1)
# Since this factor is automatically applied in the standard pipeline,
# we need it for removal and replacement by that for extended sources.
k4PPSW = cal.phot.fluxConvList.refs[2].product[array[0]].meta['k4P_'+array[0]].double
k4PPMW = cal.phot.fluxConvList.refs[2].product[array[1]].meta['k4P_'+array[1]].double
k4PPLW = cal.phot.fluxConvList.refs[2].product[array[2]].meta['k4P_'+array[2]].double

# K4e is similar to K4p but includes the frequency-dependent beam,
# making it applicable to fully extended sources with nu F nu=constant spectrum
k4EPSW = cal.phot.fluxConvList.refs[2].product[array[0]].meta['k4E_'+array[0]].double
k4EPMW = cal.phot.fluxConvList.refs[2].product[array[1]].meta['k4E_'+array[1]].double
k4EPLW = cal.phot.fluxConvList.refs[2].product[array[2]].meta['k4E_'+array[2]].double

# Therefore the conversion between exteneded and point source calibration is
K4EdivK4PPSW = k4EPSW / k4PPSW
K4EdivK4PPMW = k4EPMW / k4PPMW
K4EdivK4PPLW = k4EPLW / k4PPLW

# FWHM of point source beam profile taken from Calibration Tree
fwhmPSW = cal.phot.getBeamProf(array[0]).meta['FWHM_gMean%s'%array[0].title()].value
fwhmPMW = cal.phot.getBeamProf(array[1]).meta['FWHM_gMean%s'%array[1].title()].value
fwhmPLW = cal.phot.getBeamProf(array[2]).meta['FWHM_gMean%s'%array[2].title()].value

print "fwhmPSW", fwhmPSW
print "fwhmPMW", fwhmPMW
print "fwhmPLW", fwhmPLW

rPSW = SQRT(fitted_size**2 + fwhmPSW**2)/2.0
rPMW = SQRT(fitted_size**2 + fwhmPMW**2)/2.0
rPLW = SQRT(fitted_size**2 + fwhmPLW**2)/2.0

print "rPSW", rPSW
print "rPMW", rPMW
print "rPLW", rPLW

# radius (aperture radius needed by timeline fitter and aperture photometry)
# Default map pixel size
#c
radiusValues        = {"PSW":rPSW,  "PMW":rPMW,   "PLW":rPLW}
pixelSizeValues     = {"PSW":6.0,   "PMW":10.0,   "PLW":14.0}

radiusPSW = radiusValues["PSW"]
radiusPMW = radiusValues["PMW"]
radiusPLW = radiusValues["PLW"]
# Radii of inner and outer annulus for background estimation
# These radii are not important. I do not want to subtract any background from the total flux in the aperture used.
innerArcsec=60.0
outerArcsec=90.0

format = '%30s = %2.3f Jy'            #Format string for printout
################################################################################

################################################################################
######################## Aperture Photometry ##################################
# For optimal results begin with MJy/sr calibrated maps
# ExtdPxW maps include the relative gain correction which flatfield the
# detectors of the arrays for their integral beam profiles rather than their peaks

# Remove colour correction for extended sources and apply that for point sources
mapExtdPSW = imageDivide(image1= mapExtdPSW, scalar=K4EdivK4PPSW)
mapExtdPMW = imageDivide(image1= mapExtdPMW, scalar=K4EdivK4PPMW)
mapExtdPLW = imageDivide(image1= mapExtdPLW, scalar=K4EdivK4PPLW)

# Convert map units to [Jy/pixel] for Aperture Photometry algorithm.
mapExtdPSW = convertImageUnit(mapExtdPSW, newUnit="Jy/pixel", beamArea=beamAreaPSW)
mapExtdPMW = convertImageUnit(mapExtdPMW, newUnit="Jy/pixel", beamArea=beamAreaPMW)
mapExtdPLW = convertImageUnit(mapExtdPLW, newUnit="Jy/pixel", beamArea=beamAreaPLW)

################################################################################
#################### Running Annular Aperture Photometry #######################
# We use standard radii for aperture and background annulus
resultAperPSW = annularSkyAperturePhotometry(image=mapExtdPSW, centerRA="%s"%(ra), centerDec="%s"%(dec),\
fractional=1, centroid=False, radiusArcsec=radiusPSW,innerArcsec=innerArcsec, outerArcsec=outerArcsec)

resultAperPMW = annularSkyAperturePhotometry(image=mapExtdPMW, centerRA="%s"%(ra), centerDec="%s"%(dec),\
fractional=1, centroid=False, radiusArcsec=radiusPMW,innerArcsec=innerArcsec, outerArcsec=outerArcsec)

resultAperPLW = annularSkyAperturePhotometry(image=mapExtdPLW, centerRA="%s"%(ra), centerDec="%s"%(dec),\
fractional=1, centroid=False, radiusArcsec=radiusPLW,innerArcsec=innerArcsec, outerArcsec=outerArcsec)

# Apply beam correction, colour correction, aperture correction for given alpha
# Using "getTargetPlusSkyTotal()" means that there is not background subtraction

fluxAperPSW = resultAperPSW.getTargetPlusSkyTotal() * beamCorrPSW * aperCorrPSW * kCorrPsrcPSW
fluxAperPMW = resultAperPMW.getTargetPlusSkyTotal() * beamCorrPMW * aperCorrPMW * kCorrPsrcPMW
fluxAperPLW = resultAperPLW.getTargetPlusSkyTotal() * beamCorrPLW * aperCorrPLW * kCorrPsrcPLW

errorAperPSW = resultAperPSW.getTargetPlusSkyError() * beamCorrPSW * aperCorrPSW * kCorrPsrcPSW
errorAperPMW = resultAperPMW.getTargetPlusSkyError() * beamCorrPMW * aperCorrPMW * kCorrPsrcPMW
errorAperPLW = resultAperPLW.getTargetPlusSkyError() * beamCorrPLW * aperCorrPLW * kCorrPsrcPLW

print 'kCorrPsrcPSW =', beamCorrPSW * aperCorrPSW * kCorrPsrcPSW
print 'kCorrPsrcPMW =', beamCorrPMW * aperCorrPMW * kCorrPsrcPMW
print 'kCorrPsrcPLW =', beamCorrPLW * aperCorrPLW * kCorrPsrcPLW

print format%('Flux from Aperture Photometry: PSW', fluxAperPSW)
print format%('Flux from Aperture Photometry: PMW', fluxAperPMW)
print format%('Flux from Aperture Photometry: PLW', fluxAperPLW)

print format%('Flux Error from Aperture Photometry: PSW', errorAperPSW)
print format%('Flux Error from Aperture Photometry: PMW', errorAperPMW)
print format%('Flux Error from Aperture Photometry: PLW', errorAperPLW)

# Write to ASCII table
tds = TableDataset()
wave = Float1d([250, 350, 500])
flux = Float1d([fluxAperPSW, fluxAperPMW, fluxAperPLW])
err =  Float1d([errorAperPSW, errorAperPMW, errorAperPLW])
aperture = Float1d([2*rPSW, 2*rPMW, 2*rPLW])
tds.addColumn("wavelength(um)",Column(wave))
tds.addColumn("flux(Jy)",Column(flux))
tds.addColumn("uncertainty(Jy)",Column(err))
tds.addColumn("aperture(arcsec)",Column(aperture))
asciiTableWriter(file=outdir+"phot_sect.txt",table=tds, writeMetadata=False)
############################### End of script ##################################
################################################################################
