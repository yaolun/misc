 
#  This file is part of Herschel Common Science System (HCSS).
#  Copyright 2001-2014 Herschel Science Ground Segment Consortium
# 
#  HCSS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
# 
#  HCSS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
# 
#  You should have received a copy of the GNU Lesser General
#  Public License along with HCSS.
#  If not, see <http://www.gnu.org/licenses/>.
# 
"""    
A CHANGE
  $Id: ChopNodRangeScanPointingCorrection.py,v 1.15.2.1 2016/11/28 08:19:28 rik Exp $      

----------------------------USER INSTRUCTIONS--------------------------------
-----------------------PLEASE READ BEFORE YOU BEGIN--------------------------
-----------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!--------------------------

Available via:
 Pipelines -> PACS -> Spectrometer -> Chopped line scan -> Point Source Background Normalization
 Also to be found in your HIPE installation directory in: scripts/pacs/scripts/ipipe/spec

Purpose:
  PACS Spectrometer interactive pipeline processing from level 0 to level 2 for 
  chopped Line or short range Spectroscopy observations. Other scripts for 
  reducing PACS spectroscopy are available via the Pipeline menu in HIPE.
  This script is intended for point sources only, since it introduces a 
  correction to the telescope pointing by comparing your data to the telescope 
  beams. It uses the telescope background normalisation method for doing the 
  flux calibration.
 
  This script does work and has been tested on some targets, but it is still in 
  the process of full performance and limit testing. It is therefore recommended 
  that you check the effect of the pointing correction on your target, e.g. by comparing the 
  results of this script with those of the background normalisation pipeline script for range scans.
  Please also note that this pipeline is currently not completely documented in the PACS Data Reduction Guide.

Note: 
  This file is the reference copy, taken from your current installation of HIPE.
  If you need to edit it, we recommend to copy it to a different location first.
   
Description:
- This script starts with various setup options. In particular the OBSID
  and camera ('red' or 'blue') should be tuned everytime you run it.
- This script is supposed to be run line-by-line.
- At every moment, you can save your intermediate products to a pool with
  saveSlicedCopy(mySlicedPacsProduct, poolName[, poolLocation]), or load a 
  product previously saved with readSliced(poolName [,poolLocation])
  Note that you can't overwrite an existing pool with the same name
- Comments are included in this script, but for a full explanation see the 
  PACS Data Reduction Guide (PDRG): 
    -> Help Menu -> Show Contents -> PACS Data Reduction Guide
  The pipeline tasks themselves are explained in their "URM" entries, which you
  can find in the HIPE help page, accessed via the HIPE GUI menu
     Help -> Help Contents
  under:
     References -> PACS User's Reference Manual

Inspection:
- At each and every moment between Level 0 and Level 2, 
  slicedSummary(mySlicedPacsProduct) will return information about the contents 
  and layout of your data. If you aim at a more visual representation of your
  data structure, use slicedSummaryPlot(mySlicedPacsProduct).
- Similarly, you can call maskSummary(mySlicedPacsProduct[, slice=<number>]) to 
  get an overview of the masks that exist and/or are taken into account 
  (i.e. 'active') at that moment.
- There exist a number of additional quick and easy visualisation tools, like
  MaskViewer, plotCubes, plotCubesMasks, plotPixel, slicedPlotPointing, etc. 
  Some of them appear below, others not. See the documentation for a complete list.
- Any individual PACS' products (i.e. any slice) can be inspected in detail in 
  the "Spectrum Explorer". For instance, you may do 
	oneFrame = slicedFrames.get(slice)
	openVariable("oneFrame","Spectrum Explorer")
- The source code of the pipeline tasks can be accessed the following way:
  In the 'Tasks' panel (Hipe Menu - Window - Show View - Workbench - Tasks),
  open the 'All', right click on the Task and 'view source'.

Author: PACS ICC
"""

import os,sys
# ------------------------------------------------------------------------------
# 
# Preparation
# 
# ------------------------------------------------------------------------------
# verbose: False - silent, execute the pipeline only
#          True  - will trigger diagnostic output on the screen, plots, and displays
verbose = True

# ------------------------------------------------------------------------------
# GET THE DATA
# ------------------------------------------------------------------------------
#    
# First, set the OBSID of the observation to process.
# CHANGE THE OBSID here to your own.
#
# As this script is also run as part of the ChopNod multiObs script(s), the
# following "if" tests for the existence of a variable called multiObs, which
# will be present if you are running the multiObs script. If multiObs is
# present, the obsid will have been set already, and if not then the obsid is set
# here. (If you get a NameError, then the obsid had not been set.)
if ((not locals().has_key('multiObs')) or (not multiObs)):
	obsid = 1342229701

# Next, get the data
useHsa = 1
obs    = getObservation(obsid, verbose=True, useHsa=useHsa, poolLocation=None, poolName=None)
#if useHsa: saveObservation(obs, poolLocation=None, poolName=None)

# Show an overview of the uplink parameters of this observation 
if verbose: obsSummary(obs)

# Extract the level-0 products from the ObservationContext
pacsPropagateMetaKeywords(obs,'0', obs.level0)
level0 = PacsContext(obs.level0)

# ------------------------------------------------------------------------------
# COMPUTE A NEW POINTING PRODUCT
# ------------------------------------------------------------------------------
# If level0 was generated with SPG 11 or 12 (check in obs.meta["creator"]),
# run this block. For SPG 13 and later, don't
"""
acmsProduct = obs.auxiliary.acms
tcHistoryProduct = obs.auxiliary.teleCommHistory
oldpp = obs.auxiliary.pointing
newPP = calcAttitude(oldpp, acmsProduct, tcHistoryProduct)
obs.auxiliary.pointing = newPP
"""

# ------------------------------------------------------------------------------
# GET THE CALIBRATION FILES
# ------------------------------------------------------------------------------
# Set up the calibration tree. We take the most recent calibration files, 
# for the specific time of your observation (obs=obs)
calTree = getCalTree(obs=obs)
if verbose:
	print calTree
	print calTree.common
	print calTree.spectrometer

# ------------------------------------------------------------------------------
# SELECT DATA FROM ONE CAMERA 
# ------------------------------------------------------------------------------

# Red or blue camera ? 
if ((not locals().has_key('multiObs')) or (not multiObs)):
	camera    = 'blue'

# For your camera, extract the Frames (scientific data), the rawramps (raw data 
# for one pixel), and the DMC header (the mechanisms' status information, 
# sampled at a high frequency) 
slicedFrames  = SlicedFrames(level0.fitted.getCamera(camera).product)
slicedRawRamp = level0.raw.getCamera(camera).product
slicedDmcHead = level0.dmc.getCamera(camera).product    

# Get an overview of the basic structure of the data, prior to any processing
if verbose: slicedSummary(slicedFrames)

# ------------------------------------------------------------------------------
# SETUP -> OUTPUT
# ------------------------------------------------------------------------------

# Extract basic information
calSet = str(calTree.version)
try:
	target      = str(obs.meta["object"].value).replace(" ","").replace(".","").replace("+","plus")
	od          = "OD"+str(int(obs.meta["odNumber"].value)).zfill(4)
	hipeVersion = (str(Configuration.getProjectInfo().track) + '_' + str(Configuration.getProjectInfo().build)).replace(".","_")
except:
	target,od,hipeVersion = "","",""

# saveOutput: False - nothing is saved
#             True  - the output directory 'outputDir' will be used to store the 
#	              products of this pipeline (intermediate and final). 
# Example: outputDir = "/home/me/Herschel/PacsSpectro/pipelineOutput/"
# When saveOutput is True, nameBasis will be used as basis for the filenames of all outputs
saveOutput = True
outputDir  = str(Configuration.getWorkDir())+"/pacsSpecOut/"
if (not os.path.exists(outputDir)): os.mkdir(outputDir)
if verbose and saveOutput: print "The products of this pipeline will be saved in ",outputDir

nameBasis  = str(obsid)+"_"+target+"_"+od+"_Hipe_"+hipeVersion+"_calSet_"+calSet+"_"+camera+"_norm_poc"

# ------------------------------------------------------------------------------
#        Processing      Level 0 -> Level 0.5
# ------------------------------------------------------------------------------

# flag the saturated data in a mask "SATURATION" (and "RAWSATURATION": this 
# uses the raw data we get for some pixels)
# used cal files: RampSatLimits and SignalSatLimits
slicedFrames = specFlagSaturationFrames(slicedFrames, rawRamp = slicedRawRamp, calTree=calTree, copy=1)

# Convert digital units to Volts, used cal file: Readouts2Volts
slicedFrames = specConvDigit2VoltsPerSecFrames(slicedFrames, calTree=calTree)

# Identify the calibration blocks and fill the CALSOURCE Status entry
slicedFrames = detectCalibrationBlock(slicedFrames)

# Add the time information in UTC to the Status
slicedFrames = addUtc(slicedFrames, obs.auxiliary.timeCorrelation)

# Add the pointing information of the central spaxel to the Status
#   Uses the pointing, horizons product (solar system object ephemeries), 
#   orbitEphemeris products, and the SIAM cal file.
slicedFrames = specAddInstantPointing(slicedFrames, obs.auxiliary.pointing, calTree = calTree, orbitEphem = obs.auxiliary.orbitEphemeris, horizonsProduct = obs.auxiliary.horizons)    

# If SSO, move SSO target to a fixed position in sky. This is needed for mapping SSOs.
if (isSolarSystemObject(obs)):
  slicedFrames = correctRaDec4Sso (slicedFrames, timeOffset=0, orbitEphem=obs.auxiliary.orbitEphemeris, horizonsProduct=obs.auxiliary.horizons, linear=0)

# Extend the Status of Frames with the grating and chopper parameters GRATSCAN, CHOPPER, CHOPPOS
# used cal file: ChopperThrowDescription
slicedFrames = specExtendStatus(slicedFrames, calTree=calTree)

# Convert the chopper readouts to an angle wrt. focal plane unit and the sky
# and add this to the Status, used cal files: ChopperAngle and ChopperSkyAngle
slicedFrames = convertChopper2Angle(slicedFrames, calTree=calTree)

# Add the positions for each pixel (ra and dec datasets)
# used cal files: ArrayInstrument and ModuleArray
slicedFrames = specAssignRaDec(slicedFrames, calTree=calTree)

# Show spatial footprint
if verbose:  ppoint = slicedPlotPointing(slicedFrames, plotBoresight=False)

# Add the wavelength for each pixel (wave dataset), used cal file: WavePolynomes
slicedFrames = waveCalc(slicedFrames, calTree=calTree)

# Correct the wavelength for the spacecraft velocity. 
# Uses the pointing, orbitEphemeris and timeCorrelation product.
slicedFrames = specCorrectHerschelVelocity(slicedFrames, obs.auxiliary.orbitEphemeris, obs.auxiliary.pointing, obs.auxiliary.timeCorrelation, obs.auxiliary.horizons, force=True)

# Find the major logical blocks of this observation and organise them in the 
# BlockTable attached to the Frames; used cal file: ObcpDescription
slicedFrames = findBlocks(slicedFrames, calTree = calTree)

# Flag the known bad or noisy pixels in the masks "BADPIXELS" and "NOISYPIXELS"
# used cal files: BadPixelMask and NoisyPixelMask
slicedFrames = specFlagBadPixelsFrames(slicedFrames, calTree=calTree)

# Summary of the slices 
if verbose: slicedSummary(slicedFrames)

# Slice the data by Line/Range, Raster Point, nod position, nod cycle, on/off position and per band. 
# The parameters removeUndefined and removeMasked are for cleaning purposes
# If you want to keep the OUTOFBAND data, set removeMasked to False
slicedFrames, additionalOutContexts = pacsSliceContext(slicedFrames,[slicedDmcHead],removeUndefined=True, removeMasked=True)
slicedDmcHead = additionalOutContexts[0]

# The internal structure of your data has changed
if verbose: slicedSummary(slicedFrames)

# Flag the data affected by the chopper movement in the mask "UNCLEANCHOP"
# Uses the high resolution Dec/Mec header and the cal files ChopperAngle and ChopperJitterThreshold 
slicedFrames = flagChopMoveFrames(slicedFrames, dmcHead=slicedDmcHead, calTree=calTree)

# Flag the data effected by the grating movement in the mask "GRATMOVE"
# Uses the high resolution Dec/Mec header and the cal file GratingJitterThreshold 
slicedFrames = flagGratMoveFrames(slicedFrames, dmcHead=slicedDmcHead, calTree=calTree)

# ------------------------------------------------------------------------------
#         Processing      Level 0.5 -> Level 1
# ------------------------------------------------------------------------------

# De-activate all masks before running the glitch flagging
slicedFrames = activateMasks(slicedFrames, String1d([" "]), exclusive = True)

# Detect and flag glitches ("GLITCH" mask)
slicedFrames = specFlagGlitchFramesQTest(slicedFrames,copy=1)

# Activate all masks, for all slices
slicedFrames = activateMasks(slicedFrames, slicedFrames.get(0).getMaskTypes())

# Convert the signal level to the reference integration capacitance (which is the smallest one)
# used cal file: capacitanceRatios
slicedFrames = convertSignal2StandardCap(slicedFrames, calTree=calTree)

# Compute the differential signal of each on-off pair of datapoints, for each chopper cycle
#  -->computes the normalised difference: 2 * (on - off)/(on + off)
# The calibration block is cut out of the slicedFrames, so only the scientific slices remain.
slicedFrames = specDiffChop(slicedFrames, scical = "sci", keepall = False, normalize=True)

if verbose:
	# Data Structure. Only science blocks are left
	slicedSummary(slicedFrames)
	# Detector signal: signal is now differential
	pbasic = plotSignalBasic(slicedFrames, slice=0, titleText="plotSignalBasic")

if saveOutput:
	name=nameBasis+"_slicedFrames_B4FF"
	try:
		saveSlicedCopy(slicedFrames,name, poolLocation=outputDir)
	except:
		print "Exception raised: ",sys.exc_info()
		print "You may have to remove: ", outputDir+'/'+name
	# To restore the data:
	# slicedCubes = readSliced(name, poolLocation=outputDir)

# Compute the telescope background flux from the off-source data, and normalise the on-source data by this
# For more detail on the parameters, see the URM for this task
# The task uses the offRatio calibration product
slicedFrames, background = specRespCalToTelescope(slicedFrames, obs.auxiliary.hk, calTree = calTree, reduceNoise = 1)

if verbose:
    # Visualise the telescope background spectrum in Jy
    pBack = PlotXY(background.refs[0].product.wavelengths, background.refs[0].product.fluxJy, xtitle="Wavelength [$\mu$m]", ytitle="Flux [Jy]", titleText="Telescope Background")
    for slice in range(1,len(background.refs)):
        pBack.addLayer(LayerXY(background.refs[slice].product.wavelengths, background.refs[slice].product.fluxJy))
    #ptel = plotTelescopeModel(telescope=background, stroke=2, color=Color.black)

slicedFrames = activateMasks(slicedFrames, String1d(["SATURATION","RAWSATURATION","NOISYPIXELS","BADPIXELS","UNCLEANCHOP","GRATMOVE","GLITCH"]), exclusive = True)

# --- POINTING OFFSET CORRECTION

# The next tasks will determine the pointing, and compute a flux correction wrt the pointing error.
# The correction is based on the assumption of (hence only valid for) a point source.
# The instantaneous pointing can be estimated from the flux distribution between the central spaxels,
# compared to the flux distribution of a perfectly pointed point source (derived from the beam calfiles).
# This is the path followed when usePointingProduct = False.
# This relies entirely on the flux measured in the central 9 spaxels, hence it requires sufficient S/N.

# The pointing error can also be seen as made of two components: an absolute offset, and a jitter around it.
# The absolute offset can only be determined thanks to the tasks described above,
# but one can determine an median offset per pointing, i.e. per nod. That allows
# to determine the absolute offset from coadded data, and thus have a better result for weaker sources.
# This is the path followed when usePointingProduct = True
# The jitter is then derived from the gyro-propagated pointing products that come with the observation.
# NB: these products have been improved by new calculations (gyro-propagation), so make sure
# your level0 products have been processed with SPG 11 or higher (see the PACS and Herschel wiki)

usePointingProduct = 1

# STEP 0 
# The point source must be located on the central 3x3 spaxels for this to work
# => we must only consider those spaxels
spaxels = Int1d([6,7,8,11,12,13,16,17,18])

# Slice per grating scan to save memory when performing the next tasks
rules = [SlicingRule("ScanDir",1)]
slicedFrames = pacsSliceContext(slicedFrames, slicingRules = rules, removeUndefined=1)[0]

if verbose: slicedSummary(slicedFrames)

if saveOutput: saveSlicedCopy(slicedFrames, nameBasis+"_slicedFrames_prePOC", poolLocation=outputDir)

# STEP 1
# Determine a map that is the chi-squared of the differences between the flux of 
# your point source and the flux of the centred beam
# Task uses the various beams calibration products 
# Some parameters:
#   oversample: oversampling factor to interpolate the beams over
#       The higher the oversample, the more precise the position determination will be
#       The requested computer resources are proportional to oversample^2.
#   smoothFactor: smoothing factor along the spectral domain (working in grating positions rather than wavelength)
#       Every grating plateau (set of data-points taken at the same grating position) is medianed by default. 
#       Or you can median every "smoothFactor" grating plateau instead, rather than just one: leading to a higher S/N
#       When perNod is True, smoothFactor is ignored and a smoothFactor is calculated in a way that around 5 grating plateaux will be left. 
oversampleBeam = 5
smoothFactor1 = 10
chiSqr = specDetermineChiSquare(slicedFrames, calTree = calTree, spaxels = spaxels, oversample = oversampleBeam, smoothFactor = smoothFactor1, perNod = usePointingProduct)

# STEP 2
# Using the chiSqr, now determine the pointing offsets (i.e. those which caused the measured flux offsets)
# Some parameters:
#     chiThreshold: defines the minimum area of the chi square map: MIN(Chi2) < chi2 < chi_threshold*MIN(Chi2)
#     smoothWidth: median box width for smoothing the signal along the timeline: default is 0, i.e. no smoothing
# The output product contains the derived source position, stored in the array coordinates 
# of the beam calibration products (25x25 images with a step of 2.5 arcsec)
pointingOffset = specDeterminePointingOffset(slicedFrames, chiSqr, spaxels=spaxels, chiThreshold = 1.2, smoothWidth = 0, sigmaClip = 2.5,searchMax=1,useMedianPointingInLeak=1)

# STEP 3
# For slicedFrames that have Ra/Dec arrays calculated using the gyro-propagated pointing product the following task creates a pointing offset
# By default, we get the median pointingOffset from steps 1 & 2 above, and the jitter from the spacecraft's pointing product (here)
# Alternatively, you may also skip steps 1 & 2, and directly plug any pre-determined pointingOffset here (e.g. derived with another method, or from the other camera)
if usePointingProduct:
    smoothFactor2 = 3
    pointingOffset = specDeterminePreCalculatedPointing(slicedFrames, pointingOffset, oversample =  oversampleBeam, smoothFactor=smoothFactor2)

# STEP 4 
# Calculate the flux correction based on the pointing offset determined in the previous step 
# and add the flux correction factors to the pointingOffset product.
if usePointingProduct:
	pointingOffset = specDeterminePointingOffsetFromPreCalculatedPointing(slicedFrames, pointingOffset, calTree=calTree, spaxels=spaxels)

# Plot the pointing jitter 
if verbose:
    pjitter = plotPointingOffset(pointingOffset, calTree=calTree)

# STEP 5
# Apply the pointing corrections factors to the flux of the 9 central spaxels.
slicedFrames = specApplyPointingCorrection(slicedFrames, pointingOffset, spaxels = spaxels)

# Apply a small time- and wavelength-dependent correction to the telescope background to recover the abs. flux calibration
if calTree.version < 61 or calTree.version > 64:
    slicedFrames, telBackCor = specCorrectTelescopeBackground(slicedFrames, calTree = calTree)

# STEP 6 Re-slice the sliced products according to the standard rules
slicedFrames = pacsSliceContext(slicedFrames)[0]
point, chi = reSlicePointingProducts(slicedFrames, pointingOffset)
del(chiSqr,chi)

if verbose: slicedSummary(slicedFrames)

# --- END OF POINTING OFFSET CORRECTION

if saveOutput: 
    saveSlicedCopy(slicedFrames, nameBasis+"_slicedFrames_preFF", poolLocation=outputDir)
    saveSlicedCopy(point, nameBasis+ "_pointingOffsetProduct", poolLocation=outputDir) 

# Refine the spectral flatfield.
# For very short ranges (<~5 um, when the continuum can be approximated by a straight line)
# consider using the line-scan pipeline script.
# . excludeLeaks=True will remove wavelength ranges affected by spectral leakages.
#   It is always recommended, and mandatory when using the present script
# . selectedRange allows one to restrict the flat-field to a user-defined range
#   The spectrum outside of the selected range must be discarded from any scientific analysis.
# . useSplinesModel branches between fitting splines or polynomials.
#	- Parameter "polyOrder" is overridden when using splines.
#	- when using polynomials : polyOrder = 3 will work for most SEDs scans.

useSplinesModel = True
excludeLeaks = True
slicedFrames = specFlatFieldRange(slicedFrames,polyOrder=3, verbose=verbose, excludeLeaks=excludeLeaks, selectedRange=None, useSplinesModel=useSplinesModel)

# Convert the Frames to PacsCubes
slicedCubes = specFrames2PacsCube(slicedFrames)

# compute ra/dec meta keywords
slicedCubes = centerRaDecMetaData(slicedCubes)

# ------------------------------------------------------------------------------
#         Processing      Level 1 -> Level 2
# ------------------------------------------------------------------------------
#

# Building a wavelength grid. 
# Used cal file: wavelengthGrid
# The wavelength grid is used to perform a final sigma clipping on the spectrum (specFlagOutliers)
#     and then to rebin the spectrum (specWaveRebin)
# Note that for the final cube rebinning it is recommended that: 
#    - For single scan SED or Nyquist sampled range scans, it is recommended to perform the final 
#        rebinning with oversample=2, upsample>=1, which corresponds to the native resolution of 
#        the instrument
#    - For line scan, deep range scans or Nyquist sampled range scans with repetition factor > 1, 
#        oversampling > 2 is made possible by the high degree of redundancy provided by the observation

if ((not locals().has_key('multiObs')) or (not multiObs)):
	oversample = 2
	upsample   = 2

waveGrid=wavelengthGrid(slicedCubes, oversample=oversample, upsample=upsample, calTree = calTree)

# Activate all masks 
slicedCubes = activateMasks(slicedCubes, slicedCubes.get(0).maskTypes, exclusive = True)

# Flag the remaining outliers (sigma-clipping in wavelength domain)
slicedCubes = specFlagOutliers(slicedCubes, waveGrid, nSigma=5, nIter=1)

# Activate all masks 
slicedCubes = activateMasks(slicedCubes, slicedCubes.get(0).maskTypes, exclusive = True)

# Rebin all selected cubes on the consistent wavelength grids (mandatory for the next task specAddNod)
# To compare the rebinned spectrum to the dot cloud (possibly with various rebinning parameters), see the PDRG (Chap 3)
slicedRebinnedCubes = specWaveRebin(slicedCubes, waveGrid)

# compute ra/dec meta keywords
slicedRebinnedCubes = centerRaDecMetaData(slicedRebinnedCubes)

if verbose:
    slicedSummary(slicedRebinnedCubes)
    # Sky footprint
    # overlay parameter: Plot the footprint on an image.
    # values: [None, 'WISE', 'MSX', obsid, image]
    overlay = None 
    p9 = plotCubesRaDec(slicedRebinnedCubes, overlay = overlay)

# Average the nod-A & nod-B rebinned cubes.
# All cubes at the same raster position are averaged.
# This is the final science-grade product currently possible, for all observations except 
# the spatially oversampled rasters
slicedFinalCubes = specAddNodCubes(slicedRebinnedCubes)

# compute ra/dec meta keywords
slicedFinalCubes = centerRaDecMetaData(slicedFinalCubes)

if verbose: 
	slicedSummary(slicedFinalCubes)
	# plot all slices for a given spaxel (at this stage, normally only 1 / range and / raster position)
	x,y = 2,2
	p11 = plotCubes(slicedFinalCubes,x=x,y=y)
	#5x5 plot of one of the rebinned cubes (one line/range at one raster position)
	slice = 0
	p10 = plotCube5x5(slicedFinalCubes.refs[slice].product)


if saveOutput:
    name = nameBasis+"_slicedCubes"
    saveSlicedCopy(slicedCubes, name, poolLocation=outputDir)
    name = nameBasis+"_slicedRebinnedCubes"
    saveSlicedCopy(slicedRebinnedCubes, name, poolLocation=outputDir)
    name = nameBasis+"_slicedFinalCubes"
    saveSlicedCopy(slicedFinalCubes, name, poolLocation=outputDir)
    
# ------------------------------------------------------------------------------
#         Post-Processing
# ------------------------------------------------------------------------------

# Before you proceed with the scientific analysis, you will need to make a choice
# based on the nature of your object (extended or not), the type of observation
# that was performed (pointed or raster mapping), and its spatial sampling
# (oversampled or not). 
# - Spatially-oversampled raster maps: are projected to a single data cube, 
# thereby recovering the fullest possible spatial information. The final 
# science-grade product is a "projected" or "drizzled" cube.
# - Single pointing observations of point sources: to recover the point-source 
# calibrated spectrum, run extractCentralSpectrum. 
# - Single pointing of extended source or spatially undersampled rasters: 
# the final science-grade product is the "rebinned" cubes, or "interpolated" cubes

# ------- POINT SOURCES --------------------------------------------------------

# For single pointing observations of point sources, the calibrated spectrum is
# obtained by correcting for the fraction of the PSF that falls out of the (central) spaxel.
# That is what we call the point source correction. This is done in extractCentralSpectrum.

# extractCentralSpectrum produces three extra spectra.
# - "c1" is the central spaxel's spectrum with the point source correction applied 
# - "c9" is based on the integral of the central 9 spaxels, with the appropriate
#    point source correction applied. This spectrum, though usually more noisy, 
#    is more robust against slight pointing offsets (< 1/2 spaxel) and jitter (in range-scans). 
# - "c129" is the central spaxel, scaled to "c9". This combines the advantages 
#    of the best S/N (c1) with the robustness of c9.
#    c129 and c9 should not always be used. For low flux sources (~<10Jy 
#    continuum in central spaxel), the non-central spaxels contain essentially
#    noise resulting in uncertain 'scaling'.

# If your source is not in the central spaxel, you cannot use "extractCentralSpectrum". 
# Instead, see the dedicated script in menu 'Scripts'.

# The pointing correction applied below can be computed at every 
# wavelength, but must be smoothed in order to avoid introducing noise 
# in the resulting spectrum. There are two ways to smooth the correction:
# A. With classical filters (smoothing=='filter').
#    Two parameters control the smoothing process:
#    1.  "width" (=gaussianFilterWidth) controls the smoothing itself. If width = 0, a 
#    wavelength-independent correction (i.e. a scalar) is applied (one per range)
#    2.  "preFilterWidth" (=medianFilterWidth) is used to avoid any influence of / on the 
#    spectral lines in the correction
# B. Via wavelet decomposition. The intensity of the smoothing is controlled by
#    the parameter nLowFreq, which describes the number of low-frequency layers to 
#    keep from the wavelet decomposition to build the smoothed correction.
#    The lower nLowFreq, the smoother the result.
# NB: as written below, the call to extractCentralSpectrum demands a value
# for gaussianFilterWidth, medianFilterWidth & nLowFreq, but the non-relevant
# ones are automatically ignored by the function.

# If you wish for a correction of the shape of the continuum, you should find 
# the optimum value for these parameters, so that the correction 
# that is applied is smooth, but respects its own global shape. 
# The 'verbose plot' should help you to make a reasonable choice. Smoothing 
# too aggressively, or not enough will under/over correct the effect of the pointing jitter on the spectral shape.
# Typical values of 'width' should be of a few tens (e.g. 50), and will depend on the sampling 
# of the observation (Nyquist/SED or deep scan), over- and up-sampling adopted for the rebinning,
# the influence of the pointing jitter on your peculiar observation, etc.

smoothing = 'wavelet'
nLowFreq            = 4
# or
smoothing = 'filter'
gaussianFilterWidth = 50
medianFilterWidth   = 15

for slice in range(len(slicedFinalCubes.refs)):
    if verbose: print
    # a. Extract central spectrum, incl. point source correction (c1) 
    # b. Extract SUM(central_3x3 spaxels), incl. point source correction (c9)
    # c. Scale central spaxel to the level of the central_3x3 (c129 -> See PDRG & print extractCentralSpectrum.__doc__)
    c1,c9,c129 = extractCentralSpectrum(slicedFinalCubes.get(slice), smoothing=smoothing, width=gaussianFilterWidth, preFilterWidth=medianFilterWidth, nLowFreq=nLowFreq, calTree=calTree, verbose=verbose)
    #
    # Save to Fits
    if saveOutput:
        # ! The Pointing Offset Correction is computed for the SUM(9 central spaxels)
        # ! After it, this 'c1' spectrum should ideally not be used => we won't save it
        #name = nameBasis + "_centralSpaxel_PointSourceCorrected_Corrected3x3NO_slice_"
        #simpleFitsWriter(product=c1,file = outputDir + name + str(slice).zfill(2)+".fits")
        name = nameBasis + "_central9Spaxels_PointSourceCorrected_slice_"
        simpleFitsWriter(product=c9,file = outputDir + name + str(slice).zfill(2)+".fits")
        name = nameBasis + "_centralSpaxel_PointSourceCorrected_Corrected3x3YES_slice_"
        simpleFitsWriter(product=c129,file = outputDir + name + str(slice).zfill(2)+".fits")
