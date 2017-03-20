# 
#  This file is part of Herschel Common Science System (HCSS).
#  Copyright 2001-2011 Herschel Science Ground Segment Consortium
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
  $Id: ChopNodBackgroundNormalizationRange.py,v 1.40.2.1 2016/11/28 08:19:28 rik Exp $      

----------------------------USER INSTRUCTIONS--------------------------------
-----------------------PLEASE READ BEFORE YOU BEGIN--------------------------
-----------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!--------------------------

Available via:
 Pipeline -> PACS -> Spectrometer -> Chopped Large Range -> Background normalization

Purpose:
  PACS Spectrometer interactive pipeline processing from level 0 to level 2 for 
  chopped Line or short range spectroscopy observations. Other scripts for 
  reducing PACS spectroscopy are available via the Pipeline menu in HIPE.
  The differences between working with point and extended sources occur
  after Level 2 (formal end of the pipeline), in the post-processing part.
 
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
print 'The output pool directory is ', outputDir
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
	obsid = 1342269932

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
#outputDir  = str(Configuration.getWorkDir())+"/pacsSpecOut/"
#outputDir = '/home/bettyjo/yaolun/L1251B/ov8_up7/'
if (not os.path.exists(outputDir)): os.mkdir(outputDir)
if verbose and saveOutput: print "The products of this pipeline will be saved in ",outputDir

nameBasis  = str(obsid)+"_"+target+"_"+od+"_Hipe_"+hipeVersion+"_calSet_"+calSet+"_"+camera+"_telescope"

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
# 
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
slicedFrames = specCorrectHerschelVelocity(slicedFrames, obs.auxiliary.orbitEphemeris, obs.auxiliary.pointing, obs.auxiliary.timeCorrelation, obs.auxiliary.horizons)

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

# Refine the spectral flatfield.
# For very short ranges (<~5 um, when the continuum can be approximated by a straight line)
# consider using the line-scan pipeline script.
# . excludeLeaks=True will remove wavelength ranges affected by spectral leakages.
#   It is always recommended, and is mandatory if you use the present pipeline script
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

if verbose: slicedSummary(slicedCubes)

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

# Activate all masks - but see the PDRG Chap. 3 for more information about activating the GLITCH mask at this point in the pipeline
slicedCubes = activateMasks(slicedCubes, slicedCubes.get(0).maskTypes, exclusive = True)

# Flag the remaining outliers (sigma-clipping in wavelength domain)
slicedCubes = specFlagOutliers(slicedCubes, waveGrid, nSigma=5, nIter=1)

# Activate all masks
slicedCubes = activateMasks(slicedCubes, String1d([str(i) for i in slicedCubes.get(0).maskTypes if i not in ["INLINE", "OUTLIERS_B4FF"]]), exclusive = True)

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
# This is the final science-grade product, for all observations except the spatially 
# oversampled rasters (for point sources, see the mandatory post-processing below)
slicedFinalCubes = specAddNodCubes(slicedRebinnedCubes)

# compute ra/dec meta keywords
slicedFinalCubes = centerRaDecMetaData(slicedFinalCubes)

# Computes the telescope background flux and scales the normalized signal with the telescope background flux
slicedFinalCubes, background = specRespCalToTelescope(slicedFinalCubes, obs.auxiliary.hk, calTree = calTree)

if verbose:
    slicedSummary(slicedFinalCubes)
    # Visualise the telescope background spectrum in Jy
    pBack = PlotXY(background.refs[0].product.wavelengths, background.refs[0].product.fluxJy, xtitle="Wavelength [$\mu$m]", ytitle="Flux [Jy]", titleText="Telescope Background")
    for slice in range(1,len(background.refs)):
	pBack.addLayer(LayerXY(background.refs[slice].product.wavelengths, background.refs[slice].product.fluxJy))
    #ptel = plotTelescopeModel(telescope=background, stroke=2, color=Color.black)
    #
    # Central Spaxel for all raster positions (at this stage, normally only 1 / range and / raster position)
    x,y = 2,2
    pfinal = plotCubes(slicedFinalCubes, x=x, y=y,stroke=1,title="plotCubes - "+str(obsid)+" "+camera,subtitle="Final Rebinned Spectrum. Spaxel ["+str(x)+","+str(y)+"].\n No point source correction applied")
    #
    # Same as above, now overplotting an estimate of the continuum RMS (biased 'inside' the spectral lines).
    # For line scans, the continuum level is estimated from a polynomial of order 1 fitted outside the spectral line.
    pstd  = plotCubesStddev(slicedFinalCubes, plotLowExp=1, plotDev=0, nsigma=3, isLineScan=-1, filterWidth=20, smoothLevel=6, spaxelX=x, spaxelY=y, verbose=verbose, calTree=calTree)
    pstd.titleText, pstd.subtitleText="plotCubesStddev - "+str(obsid)+" "+camera,"Final Rebinned Spectrum & uncertainties. Spaxel ["+str(x)+","+str(y)+"].\n No point source correction applied"
    #
    # 5x5 plot of one of the rebinned cubes (one line/range at one raster position)
    slice = 0
    p55 = plotCube5x5(slicedFinalCubes.get(slice), frameTitle="plotCube5x5 - "+str(obsid)+" "+camera+" slice "+str(slice))

if saveOutput:
	name = nameBasis+"_slicedCubes"
	saveSlicedCopy(slicedCubes, name, poolLocation=outputDir)
	name = nameBasis+"_slicedRebinnedCubes"
	saveSlicedCopy(slicedRebinnedCubes, name, poolLocation=outputDir)
	name = nameBasis+"_slicedFinalCubes"
	saveSlicedCopy(slicedFinalCubes, name, poolLocation=outputDir)

# ------------------------------------------------------------------------------
#         Processing Level 2.0
# ------------------------------------------------------------------------------
#
#
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

# ------- EXTENDED SOURCES - SPATIALLY OVERSAMPLED RASTERS ---------------------

# Spatially oversampled rasters can be combined in a single (projected)cube via 
# spectral projection.
# See the PDRG Chap. 3 and the URM entry for a longer explanation of this 
# spectral projection and the parameters of the task

# SPECPROJECT
# specProject works on the final rebinned cubes.
slicedProjectedCubes = specProject(slicedFinalCubes,outputPixelsize=3.0)

# compute ra/dec meta keywords
slicedProjectedCubes = centerRaDecMetaData(slicedProjectedCubes)

if saveOutput:
  name = nameBasis+"_slicedProjectedCubes"
  saveSlicedCopy(slicedProjectedCubes, name, poolLocation=outputDir)

# Display the projected cube in the Spectrum Explorer
if verbose and (obs.obsMode=="Mapping"):
	slice = 0
	oneProjectedSlice = slicedProjectedCubes.refs[slice].product
	openVariable("oneProjectedSlice", "Spectrum Explorer")

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
        name = nameBasis + "_centralSpaxel_PointSourceCorrected_Corrected3x3NO_slice_"
        simpleFitsWriter(product=c1,file = outputDir + name + str(slice).zfill(2)+".fits")
        name = nameBasis + "_central9Spaxels_PointSourceCorrected_slice_"
        simpleFitsWriter(product=c9,file = outputDir + name + str(slice).zfill(2)+".fits")
        name = nameBasis + "_centralSpaxel_PointSourceCorrected_Corrected3x3YES_slice_"
        simpleFitsWriter(product=c129,file = outputDir + name + str(slice).zfill(2)+".fits")

# ------- UNDERSAMPLED MAPS & SINGLE POINTING OF EXTENDED SOURCES --------------

# The final science grade products are the previously produced "slicedFinalCubes"
# You can also interpolate these cubes on a regular spatial grid using specInterpolate

outputPixelsize = 3.0
slicedInterpolatedCubes = specInterpolate(slicedFinalCubes, outputPixelsize=outputPixelsize, conserveFlux=True)

# compute ra/dec meta keywords
slicedInterpolatedCubes = centerRaDecMetaData(slicedInterpolatedCubes)

if saveOutput:
	name = nameBasis+"_slicedInterpolatedCubes"
	saveSlicedCopy(slicedInterpolatedCubes, name, poolLocation=outputDir)
