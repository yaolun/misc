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
###    SPIRE Spectrometer Single Pointing User Reprocessing Script      ###
###########################################################################
#  Purpose:  A simplified version of the SPIRE Point Source SPG pipeline
#            script. This script will reprocess a specific observation
#            using the latest SPIRE calibration.
#
#  Usage:    Options for reprocessing can be set in the
#            "User Selectable Options" section at the beginning
#            of the script.
#            (A): Specify the observation to be processed with an
#                 observation ID
#            (B): By default the script only processes the centre
#                 detectors SLWC3 and SSWD4.
#                 Set processOnlyCenterDetectors = 0 to process all
#                 detectors. Note that processing all detectors will
#                 take significantly more memory and time.
#            (C): Specify an output directory, where the extended-source
#                 and point-source flux calibrated spectra will be saved
#                 as FITS files.
#            (D): ***For H+LR observations only***
#                 Choose whether to process the LR or HR data.
#                 The Default is LR.
#            (F): Loading an Observation Context into HIPE
#                 There are three options provided:
#                   1. (default) will search your local store
#                   2. will load data from a specified local pool
#                   3. will download directly from the HSA
#
#  Updated:  23/Nov/2015
#
###########################################################################
#  Changes for v10:
#     * Apodization moved into the spectral domain at the end of the pipeline
#  Changes for v11:
#     * Bright mode now uses same pipeline steps and cal files as nominal mode
#     * pointSourceSds product now contains all unvignetted detectors
#     * Temperature Drift Correction task removed
#     * FITS header items now sorted into more logical order
#  Changes for v12:
#     * Correct the frequency axis to the Local Standard of Rest
#  Changes for v13:
#     * set newThresholdCoef=True in waveletDeglitcher
#     * use new spectral apodization task
#     * write to obs context in new HIPE v13 format
#     * any overlapping SPIRE photometer observations are added to the
#       Spectrometer metaData
#  Changes for v14:
#     * New task added to correct the extended source calibration
#     * New task added to correct low resolution observations
from herschel.spire.ia.util import MetaDataDictionary
from herschel.ia.toolbox.util import RestoreTask
import os
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
            print "No photometer observations found for 0x%X"%id
            return []
    print "Cannot find photObsidForFts"
    return []
pass

#####################################
# Full AOR List for SPIRE-FTS (JDG) #
#####################################

#Object    OBSID
#DIGIT DISKS
#HD 100546 1342202273
#ABAur     1342216887
#HD 169142 1342216904
#HD 144432 1342214830
#HD 100453 1342224748
#HD 36112  1342216886
#HD 97048  1342216877
#HD 104237 1342216876
#HD 104237 1342208382 -- failed
#HD 142527 1342214821
#HD 163296 1342216906
#HD 179218 1342243607
#RY Tau    1342214857

#Klaus P.
#RNO90     1342228704


#FOOSH
#HBC722R1  1342210857
#V1735Cyg  1342219560
#V1515Cyg  1342221685
#V1331Cyg  1342221694
#V1057Cyg  1342221695
#HBC722R2  1342221696
#FUOrionis 1342230412

#COPS-SPIRE
#RCRA-IRS7B 1342242620
#RCRA-IRS7C 1342242621
#HH46       1342245084
#L723-MM    1342245094
#L1014      1342245857
#L1157      1342247625
#Ced110     1342248246
#BHR71      1342248249
#IRAS03245  1342249053
#L1551-IRS5 1342249470
#L1489      1342249473
#L1455-IRS3 1342249474
#B1-a       1342249475
#B1-c       1342249476
#IRAS03301  1342249477
#TMR1       1342250509
#TMC1A      1342250510
#L1527      1342250511
#TMC1       1342250512
#IRAS15398  1342250515
#RNO91      1342251285
#GSS30      1342251286
#VLA1623    1342251287
#IRS44/46   1342251289
#WL12       1342251290
#HH100      1342252897
#RCrA-IRS5A 1342253646
#L483       1342253649
#B335       1342253652
#DKCha      1342254037
#L1489Deep  1342214854
#L1489Off1  1342214855
#L1489Off2  1342214856
#B335- GT2  1342243609

#SPIRE HD Targets
#TWHya      1342210862

#Jeong-Eun L1251B
#L1251B     1342268303

###################################

#B335 - GT2
#Obsid = [1342243609]

#Steve Skinner - IRAS 20198+3716
#Obsid = [1342255815]

#Herbigs
#Obsid = [1342202273,1342216887,1342216904,1342214830,1342224748,1342216886,1342216877,1342216876,1342214821,1342216906,1342243607,1342214857]
#Obsid = [1342216876,1342243607]

#Herbig Failed
#Obsid = [1342208382]

#Klaus T Tauri
#Obsid = [1342228704]

#TedHD
#Obsid = [1342210862]

#Jeong-Eun
# Obsid = [1342268303]

#FOOSH list
#Obsid =[1342210857]
#Obsid =[1342219560,1342221685,1342221694,1342221695,1342221696,1342230412]

#Full COPS list
Obsid=[1342242620,1342242621,1342245084,1342245094,1342245857,
       1342247625,1342248246,1342248249,1342249053,1342249470,
       1342249474,1342249475,1342249476,1342249477,1342250509,
       1342250510,1342250512,1342250515,1342251285,1342251286,
       1342251287,1342251290,1342253646,1342253649,1342253652,
       1342254037,1342252897,1342251289]

# SECT cannot converage at L1489 1342249473, L1527 1342250511, HH100 1342252897
# mapping observation IRS44/46 1342251289

obj_list = ['RCrA-IRS7B','RCrA-IRS7C','HH46','L723-MM','L1014',
            'L1157','Ced110','BHR71','IRAS03245','L1551-IRS5',
            'L1455-IRS3','B1-a','B1-c','IRAS03301','TMR1',
            'TMC1A','TMC1','IRAS15398','RNO91','GSS30-IRS1',
            'VLA1623','WL12','RCrA-IRS5A','L483','B335',
            'DKCha','HH100','IRS44/46']

#L1489 Deep Scan
#Obsid=[1342214855,1342214856,1342214854]

#Elias29 Deep Scan
#Obsid=[1342204896,1342204897,1342204895]

# BHR71
# Obsid=[1342248249]
#Obsid=[1342249464]
useHsa=1

# collection of fitted size
size_list = []
# collection of photometry
phot_list = []

# first obsid
start_obsid = 1342251289

# the overall output directory
outdir = "/home/bettyjo/yaolun/CDF_archive_v2/"

# for reducing L1251B
Obsid = [1342268303]
obj_list = ['L1251B']
outdir = '/home/bettyjo/yaolun/L1251B/spire/'

for i in range(len(Obsid)):
    start_ind = Obsid.index(start_obsid)
    if i < start_ind:
        continue
    ###########################################################################
    ###                     User Selectable Options                         ###
    ###########################################################################
    #
    # (A) Specify OBSID:
    myObsid = Obsid[i]
    #
    # (B) Only process the centre detector if processOnlyCenterDetectors = 1.
    #     Process all detector channels if processOnlyCenterDetectors = 0:
    processOnlyCenterDetectors = 0
    #
    # (C) Specify the output directory for writing the resulting spectra into
    #     FITS files. Apodized spectra will only be saved to FITS files if apodize=1
    outDir = outdir + obj_list[i]+"/spire/data/"
    # create directory
    if not os.path.exists(outDir+'fits/'):
        os.makedirs(outDir+'fits/')
    apodize = 1
    if apodize:
        apodName = "aNB_15"
        if processOnlyCenterDetectors:
            apodName="apodized_centralpix"
        else:
            apodName = "unapod"
    #
    # (D) For H+L observations only - changing this option has no effect for
    #     observations that were not observed in "High+Low" mode:
    #     Choose whether to process LR or HR data (from a HR+LR observation)
    processRes = "HR"
    #
    # (E) N/A
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
    #     Read the latest calibration tree relevant to HCSS v14
    #     from the local disc:
    # cal = spireCal(pool="spire_cal_14_3")
    cal = spireCal(calTree="spire_cal_14_3")
    #
    #     The following command can be used to load and save the calibration tree
    #     from the HSA into a pool called "spire_cal_14_1":
    #cal = spireCal(calTree="spire_cal_14_1", saveTree=True)
    #     For more details, see the "Calibration" chapter in the SPIRE
    #     Data Reduction Guide
    #
    ###########################################################################

    # skip the cube reduction part if the FITS files already exist
    if not os.path.exists(outDir+'fits/'+str(myObsid)+'_'+processRes+'_spectrum_point_apod.fits'):

        print "Processing observation %s"%(myObsid)

        # Check that the script is appropriate for the observation
        if (obs.meta['instMode'].value  != 'SOF1' or obs.meta['obsMode'].value == 'Raster'):
           check=raw_input("Warning! You are trying to process a mapping observation with the Point Pipeline\nDo you wish to continue? (y/n)")
           if String(check).toLowerCase() == 'y':
               pass
           else:
              raise Exception("Incorrect script being used. Stopping.")
        pass

        # Attach the updated calibration tree to the observation context:
        obs.calibration.update(cal)

        # Define the central detectors:
        centreDetectors = ["SLWC3","SSWD4"]


        # Find out the bias mode of the observation (nominal/bright):
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
        lrCorr         = obs.calibration.spec.lrCorr
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
        # context, and then merge.
        bbList = SpireListContext()
        for bbid in level0_5.getBbids(0xa106):
            print "Processing building block 0x%X (%i/%i)"%(bbid, bbid-0xa1060000L, len(obs.level0_5.getBbids(0xa106)))
            sdt   = level0_5.get(bbid).sdt
            # Record the calibration tree version used by the pipeline:
            sdt.calVersion = obs.calibration.version
            nhkt  = level0_5.get(bbid).nhkt
            smect = level0_5.get(bbid).smect
            bsmt  = level0_5.get(bbid).bsmt
            # Remove all detectors except the centre ones:
            if processOnlyCenterDetectors:
                sdt = filterChannels(sdt, keepChannels=centreDetectors)
            # -----------------------------------------------------------
            # Consult the Pipeline Specification Manual for more options of the waveletDeglitcher
            # 1st level deglitching:
            sdt = waveletDeglitcher(sdt, optionReconstruction="polynomialAdaptive10", newThresholdCoef=True)
            # -----------------------------------------------------------
            # Run the Non-linearity correction:
            sdt = specNonLinearityCorrection(sdt, nonLinCorr=nonLinCorr)
            # -----------------------------------------------------------
            # Correct clipping where needed:
            sdt = clippingCorrection(sdt)
            # -----------------------------------------------------------
            # Time domain phase correction:
            sdt = timeDomainPhaseCorrection(sdt, nhkt, lpfPar=lpfPar, \
                       phaseCorrLim=phaseCorrLim, chanTimeConst=chanTimeConst)
            # -----------------------------------------------------------
            # Get pointing info:
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
            # Append this building block to the list, taking account whether the resolution
            # was H+L or not
            bbMap = SpireMapContext()
            if obs.meta["commandedResolution"].value == "H+LR":
                # for processing all scans as LR
                if processRes == "LR":
                    sdi.processResolution = "LR"
                    bbMap.setProduct("ifgm", sdi)
                    bbMap.setProduct("nhkt", nhkt)
                    bbList.addProduct(bbMap)
                # for processing the HR scans
                elif processRes == sdi.processResolution:
                    bbMap.setProduct("ifgm", sdi)
                    bbMap.setProduct("nhkt", nhkt)
                    bbList.addProduct(bbMap)
            else:
                bbMap.setProduct("ifgm", sdi)
                bbMap.setProduct("nhkt", nhkt)
                bbList.addProduct(bbMap)
        # Loop over building blocks ends here
        # -----------------------------------------------------------
        # Merge all the building blocks into one:
        merged = mergeFtsBuildingBlocks(bbList)
        sdi    = merged.getProduct("ifgm")
        nhkt   = merged.getProduct("nhkt")

        # -----------------------------------------------------------
        # Make sure all interferograms are the same length:
        sdi = makeSameOpds(sdi, opdLimits=opdLimits)

        # -----------------------------------------------------------
        # Baseline correction:
        sdi = baselineCorrection(sdi, type="fourier", threshold=4)

        # -----------------------------------------------------------
        # 2nd-level deglitching:
        sdi = deglitchIfgm(sdi, deglitchType="MAD")

        # -----------------------------------------------------------
        # Phase correction
        # The phase correction is calculated from an averaged LR interferogram:
        avgSdiFull = averageInterferogram(sdi)
        lowResSdi  = avgSdiFull.copy()
        lowResSdi.processResolution = "LR"
        lowResSdi  = makeSameOpds(lowResSdi, opdLimits=opdLimits)
        # Apply the phase correction:
        sdi = phaseCorrection(sdi, avgSdi=lowResSdi, avgSdiFull=avgSdiFull, spectralUnit="GHz")

        # -----------------------------------------------------------
        # Fourier transform:
        ssds = fourierTransform(sdi, ftType="postPhaseCorr", zeroPad="standard", \
                                spectralUnit="GHz")

        # -----------------------------------------------------------
        # Get the RSRF calibration products
        # Note: this will only work if the raw data was processed with HIPE v7 and above
        # If you get an error here, re-downloading the observation from the HSA should fix it
        instRsrf  = obs.calibration.spec.instRsrfList.getProduct(ssds)
        teleRsrf  = obs.calibration.spec.teleRsrfList.getProduct(ssds)
        beamParam = obs.calibration.spec.beamParamList.getProduct(ssds)

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
        # Apply the extended flux calibration:
        extended = specExtendedFluxConversion(ssds, teleRsrf=teleRsrf)

        # -----------------------------------------------------------
        # Correction for telescope emission:
        extended = telescopeCorrection(extended, hk=hk, teleModel=teleModel)

        # -----------------------------------------------------------
        # Apply point-source flux calibration (*copying*, rather than replacing
        # the variable "extended").
        # Keep all detectors for which calibration data exists
        pointSourceSds = filterChannels(extended.copy(), keepChannels=beamParam.uniqueDetectors)
        pointSourceSds = specPointFluxConversion(pointSourceSds, beamParam=beamParam)

        # -----------------------------------------------------------
        # If LR, apply the LR correction to the point-source calibrated data
        if pointSourceSds.processResolution=="LR":
            pointSourceSds = applyLrCorr(pointSourceSds, lrCorr=lrCorr)

        # -----------------------------------------------------------
        # Average across all scans:
        extended = averageSpectra(extended)
        pointSourceSds = averageSpectra(pointSourceSds)

        # -----------------------------------------------------------
        # Apodization:
        # For the extended calibrated data
        extended_apod = apodizeSpectra(extended.copy(), apodName="aNB_15")
        # For the point-source calibrated data
        pointSourceSds_apod = apodizeSpectra(pointSourceSds.copy(), apodName="aNB_15")

        # -----------------------------------------------------------
        # Correct the frequency scale to be in the Local Standard of Rest
        extended       = applyRadialVelocity(extended, targetFrame="lsr")
        pointSourceSds = applyRadialVelocity(pointSourceSds, targetFrame="lsr")
        # and for apodized data
        extended_apod = applyRadialVelocity(extended_apod, targetFrame="lsr")
        pointSourceSds_apod = applyRadialVelocity(pointSourceSds_apod, targetFrame="lsr")

        # -----------------------------------------------------------
        # Sort the metadata into a logical order
        extended = metaDataSorter(extended)
        pointSourceSds = metaDataSorter(pointSourceSds)
        # and for apodized data
        extended_apod = metaDataSorter(extended_apod)
        pointSourceSds_apod = metaDataSorter(pointSourceSds_apod)

        # -----------------------------------------------------------
        # Apply the extended calibration correction to the extended spectra
        extended = applyExtCalCorr(extended, specExtCorr=extCorr)
        extended_apod = applyExtCalCorr(extended_apod, specExtCorr=extCorr)

        # -----------------------------------------------------------
        # Check to see if there are any photometer observations that overlap
        # with this observation and if so, append their obsids to the metadata
        obsList = getPhotObsidsForFts(obs.obsid)
        if obsList!=None:
            #
            for l2ProdRef in obs.level2.refs:
                addMultiObsMeta(extended.meta, obsList)
                addMultiObsMeta(pointSourceSds.meta, obsList)
                addMultiObsMeta(extended_apod.meta, obsList)
                addMultiObsMeta(pointSourceSds_apod.meta, obsList)
            pass
            addMultiObsMeta(obs.meta, obsList)
        pass

        # -----------------------------------------------------------
        # Save the final spectra to FITS (both extended and point-source calibrated):
        simpleFitsWriter(extended, "%s%i_%s_spectrum_extended.fits"\
                 %(outDir+'fits/', myObsid, extended.processResolution))
        simpleFitsWriter(pointSourceSds, "%s%i_%s_spectrum_point.fits"\
                 %(outDir+'fits/', myObsid, pointSourceSds.processResolution))
        # and if required, save the apodized data
        if apodize:
            simpleFitsWriter(extended_apod, "%s%i_%s_spectrum_extended_apod.fits"\
                 %(outDir+'fits/', myObsid, extended_apod.processResolution))
            simpleFitsWriter(pointSourceSds_apod, "%s%i_%s_spectrum_point_apod.fits"\
                 %(outDir+'fits/', myObsid, pointSourceSds_apod.processResolution))

        # -----------------------------------------------------------
        # Update the Observation context format if necessary:
        obs = updateObsContext(obs)

        # -----------------------------------------------------------
        # Save the processed products back into the Observation Context.
        # This is only possible if both Level-1 and Level-2 contexts are present.
        if obs.level1 and obs.level2:
            res = pointSourceSds.processResolution
            # Save the products back into the right places inside the Observation Context
            obs.level1.getProduct("Point_0_Jiggle_0_%s"%res).setProduct("interferogram", sdi)
            obs.level2.setProduct("%s_spectrum_ext"%res, extended)
            obs.level2.setProduct("%s_spectrum_point"%res, pointSourceSds)
            obs.level2.setProduct("%s_spectrum_ext_apod"%res, extended_apod)
            obs.level2.setProduct("%s_spectrum_point_apod"%res, pointSourceSds_apod)
        # Finally we can save the new reprocessed observation back to your hard disk.
        # Uncomment the next line and choose a poolName, either the existing one or a new one
        #saveObservation(obs, poolName="enter-a-poolname", saveCalTree=True)

    ############################################
    #  SECT Corrections along with Photometry  #
    ############################################

    ############################### Setup #####################################

    from herschel.share.unit import Constant
    c  = Constant.SPEED_OF_LIGHT.value
    from java.lang.Math import PI

    ############################### Setup #####################################



    if obj_list[i] == 'HH100':
        continue

    spec = obs.refs['level2'].product.refs['HR_spectrum_point_apod'].product

    sect_spec = semiExtendedCorrector(spectrum=spec, doPlots=False, calibration=cal,
    	optimiseDiameter=True, sourceModel=SemiExtendedSourceModelShape("gaussian", 30.0, 0.0, 0.0, 0.0, 0.0, 1.0, 257, 257))

    fitted_size = sect_spec.meta['sourceDiameter'].double

    correctedSpectrum = semiExtendedCorrector(spectrum=spec, calibration=cal, optimiseDiameter=True, doPlots=False,\
     	gaussRefBeamDiam=fitted_size, sourceModel=SemiExtendedSourceModelShape("gaussian", fitted_size, 0.0, 0.0, 0.0, 0.0, 1.0, 257, 257))
    ra  = spec.meta["ra"].value
    dec = spec.meta["dec"].value
    simpleFitsWriter(correctedSpectrum, "%s%i_%s_spectrum_SECT.fits" \
    	%(outDir+'fits/', myObsid, processRes))
    spire_spec = convertWavescale(ds=spireProduct2SimpleSpectrum(input=correctedSpectrum),
                                    to='micrometer', referenceUnit='kHz')
    exportSpectrumToAscii(ds=spire_spec,file=outDir+str(myObsid)+'_spire_sect.txt',meta=False)
    # check whether it has photometry
    obsList = getPhotObsidsForFts(spec.obsid)

    # Write out the fitted size and the OBSID of the photometry data.
    # if len(obsList) == 0:
    #     misc_data = Float1d(fitted_size)
    # else:
    #     misc_data = Float1d([fitted_size]).append(Float1d(obsList))
    tds = TableDataset()
    tds.addColumn("size", Column(Float1d([fitted_size])))
    tds.addColumn("phot_obsid", Column(obsList))
    asciiTableWriter(file=outDir+str(myObsid)+"_spire_sect_meta.txt",table=tds, writeMetadata=False)

    print "Processed of observation %s complete :-)"%(myObsid)
