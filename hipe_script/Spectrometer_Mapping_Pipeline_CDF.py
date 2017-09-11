# 
#  This file is part of Herschel Common Science System (HCSS).
#  Copyright 2001-2015 Herschel Science Ground Segment Consortium
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
###########################################################################
###          SPIRE Spectrometer Mapping User Reprocessing Script        ###
###########################################################################
#  Purpose:  A simplified version of the SPIRE Mapping SPG pipeline script. 
#            This script will reprocess a specific observation using the 
#            latest SPIRE calibration.          
#
#  Usage:    Options for reprocessing can be set in the 
#            "User Selectable Options" section at the beginning 
#            of the script. 
#            (A): Specify the observation to be processed with an 
#                 observation ID
#            (C): Specify an output directory, where the the preCubes 
#                 (Spectrum2ds that lists the individual spectra) and 
#                 the spectral cubes will be saved as FITS files for 
#                 each array.
#            (D): ***For H+LR observations only*** 
#                 Choose whether to process the LR or HR data. 
#                 The Default is LR.
#            (E): The spatial pixel size of the spectral cubes can be 
#                 changed
#            (F): Loading an Observation Context into HIPE 
#                 There are three options provided:
#                   1. (default) will search your local store
#                   2. will load data from a specified local pool
#                   3. will download directly from the HSA
#
#  Updated:  08/Mar/2016
#
###########################################################################
#  Changes for v10:
#     * Apodization moved into the spectral domain at the end of the pipeline
#  Changes for v11:
#     * Bright mode now uses same pipeline steps and cal files as nominal mode
#     * Temperature Drift Correction task removed
#     * FITS header items now sorted into more logical order
#     * ssds products now averaged before making the cube
#     * Small change to algorithm to calculate optimum WCS
#  Changes for v12:
#     * Include vignetted detectors in maps
#     * Correct the frequency axis to the Local Standard of Rest
#  Changes for v13:
#     * set newThresholdCoef=True in waveletDeglitcher
#     * use new spectral apodization task
#     * write to obs context in new HIPE v13 format
#     * any overlapping SPIRE photometer observations are added to the 
#       Spectrometer metaData
#     * the individual extended calibrated spectra are now stored in the 
#       preCube variable, and moved from Level-1 to Level-2 in the 
#       observation context
#  Changes for v14:
#     * New task added to correct the extended source calibration
#     * Spectral cubes gridded with the convolution algorithm are now produced in 
#       addition to those gridded with the naive projection algorithm
from herschel.spire.ia.util import MetaDataDictionary 
from herschel.ia.toolbox.util import RestoreTask
restore = RestoreTask()
def addMultiObsMeta(meta, obsList):
    """
    Insert metadata for related observations.

    Input:
        meta  metadata to be added to
        obsList  list of observation ids
    """
    dict = MetaDataDictionary.getInstance()
    myBaseParam = 'photObsid'
    for i in range(0,len(obsList)):
        param = dict.newParameter(myBaseParam, Long(obsList[i]))
        param.setDescription("%s %03d"%(param.getDescription(), i))
        meta.set("%s%03d"%(myBaseParam, i), param)
pass

def getPhotObsidsForFts(id):
    """
    Get a list of photometer obsids that a given FTS obsid is contained in.
    Input:
        id FTS obsid as a number
    Output:
        list of photometer obsids that overlap with the FTS one
    """
    if not globals().has_key("photObsidForFts"):
        hcssDir = Configuration.getProperty("var.hcss.dir")
        restore(hcssDir+"/data/spire/ia/pipeline/scripts/merging/photObsidsForFts.ser")
    key = "0x%X"%id
    if globals().has_key("photObsidForFts"):
        if photObsidForFts.has_key(key):
            return photObsidForFts[key]
        else:
            print "No photometer observations found for %i"%id
            return []
    print "Cannot find photObsidForFts"
    return []
pass

lrSlwChannels = [ \
    "SLWB2",
    "SLWB3",
    "SLWC2",
    "SLWC3",
    "SLWC4",
    "SLWD2",
    "SLWD3",
]

lrSlwChannelsToRemove = [ \
    "SLWA1",
    "SLWA2",
    "SLWA3",
    "SLWB1",
    "SLWB4",
    "SLWC1",
    "SLWC5",
    "SLWD1",
    "SLWD4",
    "SLWE1",
    "SLWE2",
    "SLWE3",
]
###########################################################################
###                     User Selectable Options                         ###
###########################################################################
#
# (A) Specify OBSID:
myObsid = 1342251289
# L1251B
myObsid = 1342268303
#
# (B) N/A
#
# (C) Specify the output directory for writing the resulting spectra and 
#     cubes into FITS files
#     Apodized spectra will only be saved to FITS files if apodize=1
outDir = "/home/bettyjo/yaolun/L1251B/spire/hipe14/"
apodize = 1
#
# (D) For H+L observations only - changing this option has no effect for
#     observations that were not observed in "High+Low" mode:
#     Choose whether to process LR or HR data (from a HR+LR observation)
processRes="LR"
#
# (E) Specify the map pixel size for the final data cubes (SSW and SLW)
#     depending on the spatial sampling in units of degree:
sswFwhm = 19.0 / 3600.0
slwFwhm = 35.0 / 3600.0
gridSpacing={"full":        {"SSW": 0.5 * sswFwhm, "SLW": 0.5 * slwFwhm}, \
             "intermediate":{"SSW": 1.0 * sswFwhm, "SLW": 1.0 * slwFwhm}, \
             "sparse":      {"SSW": 2.0 * sswFwhm, "SLW": 2.0 * slwFwhm}}
#
# (F) Load an observation context into HIPE:
#     1. To search for the observation in your local store use:
obs = getObservation(myObsid, useHsa=True, save=True)
#     2. To specify a pool name in your local store use:
#obs = getObservation(myObsid, poolName="poolName") 
#     3. To load data directly from the HSA use:
#obs = getObservation(myObsid, useHsa=True)
#
# (G) Calibration Context and Calibration Files 
#     Load the latest calibration tree from the HSA
cal = spireCal(calTree = "spire_cal_14_3")
#
#     To load the calibration tree from the HSA 
#     and save it as a pool use:
#cal = spireCal(calTree="spire_cal", saveTree=True)
#
#     If the latest calibration tree is already saved
#     locally use:
#cal = spireCal(pool="poolName")
#
#     For more details on loading and saving the SPIRE
#     calibration Context, see the "Calibration" chapter 
#     in the SPIRE Data Reduction Guide
#
############################################################################

print "Processing observation %i"%(myObsid)

# Check that the script is appropriate for the observation
if (obs.meta['instMode'].value == 'SOF1' and obs.meta['obsMode'].value != 'Raster'):
   check=raw_input("Warning! You are trying to process a point observation with the Mapping Pipeline\nDo you wish to continue? (y/n)")
   if String(check).toLowerCase() == 'y':
       pass
   else:
      raise Exception("Incorrect script being used. Stopping.")
pass 

# Attach the updated calibration tree to the observation context
obs.calibration.update(cal)

# Find out the bias mode of the observation (nominal/bright)
biasMode = obs.meta["biasMode"].value

# Extract necessary Calibration Products from the Observation Context
nonLinCorr     = obs.calibration.spec.nonLinCorr
chanNum        = obs.calibration.spec.chanNum
bolPar         = obs.calibration.spec.bolPar
lpfPar         = obs.calibration.spec.lpfPar
phaseCorrLim   = obs.calibration.spec.phaseCorrLim
chanTimeConst  = obs.calibration.spec.chanTimeConst
bsmPos         = obs.calibration.spec.bsmPos
detAngOff      = obs.calibration.spec.detAngOff
smecZpd        = obs.calibration.spec.smecZpd
chanTimeOff    = obs.calibration.spec.chanTimeOff
smecStepFactor = obs.calibration.spec.smecStepFactor
opdLimits      = obs.calibration.spec.opdLimits
bandEdge       = obs.calibration.spec.bandEdge
brightGain     = obs.calibration.spec.brightGain
extCorr        = obs.calibration.spec.extCorr
# teleModel contains the OPD-dependent emissivity correction
# factors that are applied to the telescope model.
teleModel      = obs.calibration.spec.teleModel

# Extract necessary Auxiliary Products from the Observation Context
hpp  = obs.auxiliary.pointing
siam = obs.auxiliary.siam
hk   = obs.auxiliary.hk

# Get the Level-0.5 data directly from the observation
level0_5 = obs.level0_5

# Start to process the observation from Level 0.5
# Process each SMEC scan building block (0xa106) individually, append to a list 
sdsList = SpireMapContext()
sdsList_apod = SpireMapContext()

# Propagate the metadata from the Level 0.5 context to the lists 
for key in level0_5.meta.keySet():
    if key != "creator" and key != "creationDate" and key != "fileName" and \
    key != "type" and key != "description": 
        sdsList.meta[key]=level0_5.meta[key].copy()
        sdsList_apod.meta[key]=level0_5.meta[key].copy()

# Each building block is a different jiggle position
for bbid in level0_5.getBbids(0xa106):
    print"Processing building block 0x%X (%i/%i)"%(bbid, bbid-0xa1060000L, \
                           len(obs.level0_5.getBbids(0xa106)))
    sdt   = level0_5.get(bbid).sdt
    # Record the calibration tree version used by the pipeline:
    sdt.calVersion = obs.calibration.version
    nhkt  = level0_5.get(bbid).nhkt
    smect = level0_5.get(bbid).smect
    bsmt  = level0_5.get(bbid).bsmt
    # Extract the jiggle ID from the metadata:
    jiggId = sdt.meta["jiggId"].value
    # Extract raster ID from the metadata:
    raster = sdt.meta["pointNum"].value
    # -----------------------------------------------------------
    # 1st level deglitching:    
    # Consult the Pipeline Specification Manual for more options 
    # of the waveletDeglitcher
    sdt = waveletDeglitcher(sdt, optionReconstruction="polynomialAdaptive10",\
                            newThresholdCoef=True)
    # -----------------------------------------------------------
    # Run the Non-linearity correction:
    sdt = specNonLinearityCorrection(sdt, nonLinCorr=nonLinCorr)
    # -----------------------------------------------------------
    # Repair clipped samples where needed:
    sdt = clippingCorrection(sdt)
    # -----------------------------------------------------------
    # Time domain phase correction:
    sdt = timeDomainPhaseCorrection(sdt, nhkt, lpfPar=lpfPar, \
               phaseCorrLim=phaseCorrLim, chanTimeConst=chanTimeConst)        
    # -----------------------------------------------------------
    # Add pointing info:
    bat = calcBsmAngles(bsmt, bsmPos=bsmPos)
    spp = createSpirePointing(hpp=hpp, siam=siam, \
                            detAngOff=detAngOff, bat=bat)
    # -----------------------------------------------------------
    # Create interferogram:
    sdi = createIfgm(sdt, smect=smect, nhkt=nhkt, spp=spp, \
                     smecZpd=smecZpd,\
                     chanTimeOff=chanTimeOff,\
                     smecStepFactor=smecStepFactor)
    # -----------------------------------------------------------
    # Update the resolution if processing a H+L observation as LR
    if obs.meta["commandedResolution"].value == "H+LR" and processRes == "LR":
        sdi.processResolution = "LR"
    # Adjust OPD ranges to ensure that they are the same for all scans
    sdi = makeSameOpds(sdi, opdLimits=opdLimits)
    # -----------------------------------------------------------
    # Baseline correction:
    sdi = baselineCorrection(sdi, type="fourier", threshold=4)
    # -----------------------------------------------------------
    # 2nd level deglitching:
    sdi = deglitchIfgm(sdi, deglitchType="MAD")
    # -----------------------------------------------------------
    # Phase correction
    # The phase correction is calculated from an averaged LR interferogram:
    avgSdiFull = averageInterferogram(sdi)
    lowResSdi  = avgSdiFull.copy()
    lowResSdi.processResolution = "LR"
    lowResSdi  = makeSameOpds(lowResSdi, opdLimits=opdLimits)
    # Apply the phase correction:
    sdi = phaseCorrection(sdi, avgSdi=lowResSdi, avgSdiFull=avgSdiFull, \
                            spectralUnit="GHz")
    # -----------------------------------------------------------
    # Fourier transform:
    ssds = fourierTransform(sdi, ftType="postPhaseCorr", zeroPad="standard", \
                            spectralUnit="GHz")
    # -----------------------------------------------------------
    # Get the RSRF calibration products
    # Note: this will only work if the raw data was processed with 
    # HIPE v7 and above. If you get an error here, re-downloading 
    # the observation from the HSA should fix it
    instRsrf = obs.calibration.spec.instRsrfList.getProduct(ssds)
    teleRsrf = obs.calibration.spec.teleRsrfList.getProduct(ssds)
    # -----------------------------------------------------------
    # Remove out of band data:
    ssds = removeOutOfBand(ssds, bandEdge=bandEdge)
    # -----------------------------------------------------------
    # Apply bright gain correction (bright source setting only):
    if biasMode == "bright":
        ssds = specApplyBrightGain(ssds, brightGain=brightGain)
    # -----------------------------------------------------------
    # Correction for instrument emission:
    ssds = instCorrection(ssds, nhkt=nhkt, instRsrf=instRsrf)
    # -----------------------------------------------------------
    # Get the flux conversion calibration products and apply to spectra:
    ssds = specExtendedFluxConversion(ssds, teleRsrf=teleRsrf)
    # -----------------------------------------------------------
    # Correction for telescope emission:
    ssds = telescopeCorrection(ssds, hk=hk, teleModel=teleModel)
    # -----------------------------------------------------------
    # If LR, apply the LR correction
    if ssds.processResolution=="LR":
        lrCorr = obs.calibration.spec.lrCorr
        beamParam = obs.calibration.spec.beamParamList.getProduct(ssds)
        ssdsFull = ssds.copy()
        ssdsFull = filterChannels(ssdsFull, removeChannels=String1d(lrSlwChannels))
        ssds = filterChannels(ssds, keepChannels=String1d(lrSlwChannels))
        ssds = specPointFluxConversion(ssds, beamParam=beamParam)
        ssds = applyLrCorr(ssds, lrCorr=lrCorr)
        ssds = specLrExtFluxConv(ssds, ssdsFull, beamParam=beamParam)
    # -----------------------------------------------------------
    # Average across all scans:
    ssds = averageSpectra(ssds)
    # ----------------------------------------------------------- 
    # Apodization:
    ssds_apod = apodizeSpectra(ssds.copy(), apodName="aNB_15")
    # -----------------------------------------------------------
    # Correct the frequency scale to be in the Local Standard of Rest
    ssds = applyRadialVelocity(ssds, targetFrame="lsr")
    ssds_apod = applyRadialVelocity(ssds_apod, targetFrame="lsr")
    # -----------------------------------------------------------
    # Sort the metadata into a logical order
    ssds = metaDataSorter(ssds)
    ssds_apod = metaDataSorter(ssds_apod)
    # ----------------------------------------------------------- 
    # Append this scan to the list, taking account whether the
    # resolution was H+L or not
    if obs.meta["commandedResolution"].value == "H+LR":
        # for processing all scans as LR
        if processRes == "LR":
            sdsList.setProduct("%d"%(sdsList.refs.size()), ssds)
            sdsList_apod.setProduct("%d"%(sdsList_apod.refs.size()), ssds_apod)
        # for processing the HR scans
        elif processRes == ssds.processResolution:
            sdsList.setProduct("%d"%(sdsList.refs.size()), ssds)
            sdsList_apod.setProduct("%d"%(sdsList_apod.refs.size()), ssds_apod)
    # or, otherwise
    else:
        sdsList.setProduct("%d"%(sdsList.refs.size()),ssds)
        sdsList_apod.setProduct("%d"%(sdsList_apod.refs.size()),ssds_apod)
        #
    # -----------------------------------------------------------
    # Save the interferogram data back into the Observation Context:
    if obs.level1:
        res = ssds.processResolution
        obs.level1.getProduct("Point_%i_Jiggle_%i_%s"%(raster,\
            jiggId, res)).setProduct("interferogram", sdi)
        # Remove the old style products (pre-HIPE13) if they exist
        obs.level1.getProduct("Point_%i_Jiggle_%i_%s"%(raster,\
            jiggId, res)).refs.remove("apodized_spectrum")
        obs.level1.getProduct("Point_%i_Jiggle_%i_%s"%(raster,\
            jiggId, res)).refs.remove("unapodized_spectrum")

# ---------------------------------------------------------------
# Carry out regridding
mapSampling = obs.meta['mapSampling'].value
for array in ["SSW", "SLW"]:
    if array=="SLW":
        beamDiamDet="SLWC3"
    if array=="SSW":
        beamDiamDet="SSWD4"
    # -----------------------------------------------------------
    # Create a listing of all spectra and positions in a spectrum2d
    preCube = spirePreprocessCube(context=sdsList, arrayType=array, \
                             unvignetted=False)
    preCube_apod = spirePreprocessCube(context=sdsList_apod, arrayType=array, \
                             unvignetted=False)
    # -----------------------------------------------------------
    # Apply the extended calibration correction to extended spectra
    preCube = applyExtCalCorr(preCube, specExtCorr=extCorr)
    preCube_apod = applyExtCalCorr(preCube_apod, specExtCorr=extCorr)
    # -----------------------------------------------------------
    # Set up the grid - covering the RA and Dec of observed points 
    # using specified gridSpacing:
    wcs = SpecWcsCreator.createWcs(preCube, gridSpacing[mapSampling][array], \
                              gridSpacing[mapSampling][array])
    # -----------------------------------------------------------
    # If LR, remove the SLW vignetted detectors from cube
    if preCube.meta["processResolution"].value=="LR" and array=="SLW":
        preCubeCopy = preCube.copy()
        preCube.deleteSpectra(lrSlwChannelsToRemove)
        preCube_apodCopy = preCube_apod.copy()
        preCube_apod.deleteSpectra(lrSlwChannelsToRemove)
    pass
    # -----------------------------------------------------------
    # Regrid the data using the Naive Projection algorithm:
    cube = spireProjection(spc=preCube, wcs=wcs, projectionType="naive")
    cube_apod = spireProjection(spc=preCube_apod, wcs=wcs, projectionType="naive")
    # -----------------------------------------------------------
    # Regrid the data using the Convolution algorithm:
    beamParam = obs.calibration.spec.beamParamList.getProduct(preCube.meta["processResolution"].value, \
                               preCube.startDate)
    beamDiam = beamParam[0][beamDiamDet]["beamDiam"].data
    cube_convol = spireProjection(spc=preCube, wcs=wcs, projectionType="convolution", \
                               beamWidthArray=beamDiam)
    cube_convol_apod = spireProjection(spc=preCube_apod, wcs=wcs, \
                               projectionType="convolution", beamWidthArray=beamDiam)
    # -----------------------------------------------------------
    # Sort the metadata into a logical order
    cube = metaDataSorter(cube)
    cube_apod = metaDataSorter(cube_apod)
    cube_convol = metaDataSorter(cube_convol)
    cube_convol_apod = metaDataSorter(cube_convol_apod)
    # -----------------------------------------------------------
    # Check to see if there are any photometer observations that overlap
    # with this observation and if so, append their obsids to the metadata
    obsList = getPhotObsidsForFts(obs.obsid)
    # -----------------------------------------------------------
    # If LR, save the full preCubes
    if preCube.meta["processResolution"].value=="LR" and array=="SLW":
        preCube = preCubeCopy.copy()
        preCube_apod = preCube_apodCopy.copy()
    if obsList!=None:
        #
        for l2ProdRef in obs.level2.refs:
            addMultiObsMeta(preCube.meta, obsList)
            addMultiObsMeta(cube.meta, obsList)
            addMultiObsMeta(cube_convol.meta, obsList)
            addMultiObsMeta(preCube_apod.meta, obsList)
            addMultiObsMeta(cube_apod.meta, obsList)
            addMultiObsMeta(cube_convol_apod.meta, obsList)
        pass    
        addMultiObsMeta(obs.meta, obsList)
    pass
    # -----------------------------------------------------------
    # Tweak the cube variable name to specify the array:
    exec("cube%s = cube"%array)
    exec("cube%s_apod = cube_apod"%array)
    exec("cube%s_convol = cube_convol"%array)
    exec("cube%s_convol_apod = cube_convol_apod"%array)
    # -----------------------------------------------------------
    # Save the data back to the Observation Context, making sure
    # that it is saved in the new (from HIPE v13 onwards) format.
    if obs.level2:
        res = preCube.meta["processResolution"].value
        # Save the preCube
        obs.level2.setProduct("%s_%s_spectrum2d"%(res, array), preCube)
        obs.level2.setProduct("%s_%s_spectrum2d_apod"%(res, array), preCube_apod)
        # Save the cube (removing old style products if they exist)
        obs.level2.refs.remove("%s_%s_unapodized_spectrum"%(res, array))
        obs.level2.refs.remove("%s_%s_apodized_spectrum"%(res, array))
        obs.level2.setProduct("%s_%s_cube"%(res, array), cube)
        obs.level2.setProduct("%s_%s_cube_convol"%(res, array), cube_convol)
        obs.level2.setProduct("%s_%s_cube_apod"%(res, array), cube_apod)
        obs.level2.setProduct("%s_%s_cube_convol_apod"%(res, array), cube_convol_apod)
    # -----------------------------------------------------------
    # Save the full preCubes and cubes to FITS:
    simpleFitsWriter(preCube, "%s%i_%s_%s_spectrum2d.fits"%(outDir, myObsid, \
                     res, array))
    simpleFitsWriter(cube, "%s%i_%s_%s_cube.fits"%(outDir, myObsid, \
                     res, array))
    simpleFitsWriter(cube_convol, "%s%i_%s_%s_cube_convol.fits"%(outDir, myObsid, \
                     res, array))
    # and if required, save the apodized data
    if apodize:
        simpleFitsWriter(preCube_apod, "%s%i_%s_%s_spectrum2d_apod.fits"%(outDir, myObsid, \
                     res, array))
        simpleFitsWriter(cube_apod, "%s%i_%s_%s_cube_apod.fits"%(outDir, myObsid, \
                     res, array))
        simpleFitsWriter(cube_convol_apod, "%s%i_%s_%s_cube_convol_apod.fits"%(outDir, myObsid, \
                     res, array))

# Finally we can save the new reprocessed observation back to your hard disk.
# Uncomment the next line and choose a poolName, either the existing one or a new one
#saveObservation(obs, poolName="enter-a-poolname", saveCalTree=True)

print "Processing of observation %i complete"%(myObsid)
#### End of the script ####
