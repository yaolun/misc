
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
# Obsid=[1342242620,1342242621,1342245084,1342245094,1342245857,1342247625,1342248246,1342248249,\
# 1342249053,1342249470,1342249473,1342249474,1342249475,1342249476,1342249477,1342250509,\
# 1342250510,1342250511,1342250512,1342250515,1342251285,1342251286,1342251287,1342251289,\
# 1342251290,1342252897,1342253646,1342253649,1342253652,1342254037]

# COPS list selected from successful SECT Reduction
Obsid=[1342242620,1342242621,1342245084,1342245094,1342245857,
       1342247625,1342248246,1342248249,1342249053,1342249470,
       1342249473,1342249474,1342249475,1342249476,1342249477,
       1342250509,1342250510,1342250512,1342250515,1342251285,
       1342251286,1342251287,1342251290,1342253646,1342253649,
       1342253652,1342254037]

obj_list = ['RCrA-IRS7B','RCrA-IRS7C','HH46','L723-MM','L1014',
            'L1157','Ced110','BHR71','IRAS03245','L1551-IRS5',
            'L1489','L1455-IRS3','B1-a','B1-c','IRAS03301',
            'TMR1','TMC1A','TMC1','IRAS15398','RNO91',
            'GSS30-IRS1','VLA1623','WL12','RCrA-IRS5A','L483',
            'B335','DKCha']

phot_list = [[1342216002,1342206678,1342206677], [1342216002,1342206678,1342206677], 0, [1342229605], [1342220631,1342219974],
             [1342189844,1342189843],[1342213179,1342213178], [1342226633], [1342190327,1342190326], [1342202251,1342202250],
             [1342190327,1342190326], [1342190327,1342190326], [1342190327,1342190326], [1342190327,1342190326], [1342202253,1342202252],
             [1342202253,1342202252], [1342202253,1342202252], [1342213183,1342213182], [1342263845,1342263844], [1342205094,1342205093,1342203074],
             [1342205094,1342205093,1342203074], [1342205094,1342205093], [1342205094,1342205093], [1342216002,1342206678,1342206677], [1342229186],
             [1342192685], [1342213181,1342213180]]

size_list =[37.0, 38.0, 16.0, 12.75, 38.0,
            11.0, 30.5, 15.5, 14.5, 14.25,
            54.0, 38.0, 14.0, 49.0, 26.5,
            18.25, 30.5, 35.0, 21.0, 38.5,
            25.75, 41.0, 41.0, 40.0, 22.0,
            13.5, 13.25]


############################### Input #####################################
indir = '/home/bettyjo/yaolun/CDF_SPIRE_reduction/'
outdir = '/home/bettyjo/yaolun/CDF_SPIRE_reduction/photometry/'
fitted_size = 0
# Radii of inner and outer annulus for background estimation
# These radii are not important. I do not want to subtract any background from the total flux in the aperture used.
innerArcsec=60.0
outerArcsec=90.0

############################### Setup #####################################

from herschel.share.unit import Constant
c  = Constant.SPEED_OF_LIGHT.value
from java.lang.Math import PI

for obsidFTS in Obsid:
    print 'processing ', obsidFTS, obj_list[Obsid.index(obsidFTS)]
    print 'fitted size is ', size_list[Obsid.index(obsidFTS)]
    fitted_size = size_list[Obsid.index(obsidFTS)]
    obs = getObservation(obsid=obsidFTS, useHsa=True)
    spec = obs.refs['level2'].product.refs['HR_spectrum_point_apod'].product

    ra  = spec.meta["ra"].value
    dec = spec.meta["dec"].value

    # Run for an individual SPIRE band
    array = ["PSW", "PMW","PLW"]  # SPIRE Array Bands: "PSW", "PMW", "PLW"

    ############################ Import Data ##################################
    # Loading an observation of Gamma Dra from the HSA
    phot_obs = phot_list[Obsid.index(obsidFTS)]
    wave = []
    phot = []
    error = []
    phot_aper = []
    alpha_data = asciiTableReader(file=indir+str(obsidFTS)+'_alpha.txt', tableType='SPACES')
    alpha = [float(alpha_data[0].data[1]), float(alpha_data[1].data[1]), float(alpha_data[2].data[1])]

    if phot_obs == 0:
        print 'No photometry data found for ', obj_list[Obsid.index(obsidFTS)]
        continue

    for phot_o in phot_obs:
        obs     = getObservation(phot_o, useHsa=True, instrument='SPIRE')
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

        # print "beamCorrPSW = ",beamCorrPSW
        # print "beamCorrPMW = ",beamCorrPMW
        # print "beamCorrPLW = ",beamCorrPLW

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

        # print "fwhmPSW", fwhmPSW
        # print "fwhmPMW", fwhmPMW
        # print "fwhmPLW", fwhmPLW

        rPSW = SQRT(fitted_size**2 + fwhmPSW**2)/2.0
        rPMW = SQRT(fitted_size**2 + fwhmPMW**2)/2.0
        rPLW = SQRT(fitted_size**2 + fwhmPLW**2)/2.0

        # print "rPSW", rPSW
        # print "rPMW", rPMW
        # print "rPLW", rPLW

        # radius (aperture radius needed by timeline fitter and aperture photometry)
        # Default map pixel size
        #c
        radiusValues        = {"PSW":rPSW,  "PMW":rPMW,   "PLW":rPLW}
        pixelSizeValues     = {"PSW":6.0,   "PMW":10.0,   "PLW":14.0}

        radiusPSW = radiusValues["PSW"]
        radiusPMW = radiusValues["PMW"]
        radiusPLW = radiusValues["PLW"]

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

        # print 'kCorrPsrcPSW =', beamCorrPSW * aperCorrPSW * kCorrPsrcPSW
        # print 'kCorrPsrcPMW =', beamCorrPMW * aperCorrPMW * kCorrPsrcPMW
        # print 'kCorrPsrcPLW =', beamCorrPLW * aperCorrPLW * kCorrPsrcPLW

        print format%('Flux from Aperture Photometry: PSW', fluxAperPSW)
        print format%('Flux from Aperture Photometry: PMW', fluxAperPMW)
        print format%('Flux from Aperture Photometry: PLW', fluxAperPLW)

        print format%('Flux Error from Aperture Photometry: PSW', errorAperPSW)
        print format%('Flux Error from Aperture Photometry: PMW', errorAperPMW)
        print format%('Flux Error from Aperture Photometry: PLW', errorAperPLW)

        wave.extend([250,350,500])
        phot.extend([fluxAperPSW, fluxAperPMW, fluxAperPLW])
        error.extend([errorAperPSW, errorAperPMW, errorAperPLW])
        phot_aper.extend([2*rPSW, 2*rPMW, 2*rPLW])

    # Write to ASCII table
    wave = Float1d(wave)
    phot = Float1d(phot)
    error = Float1d(error)
    phot_aper = Float1d(phot_aper)
    phot_obsidForwrite = []
    for o in phot_obs:
        phot_obsidForwrite.extend([o,o,o])
    phot_obsidForwrite = String1d(phot_obsidForwrite)
    #
    tds = TableDataset()
    tds.addColumn("OBSID", Column(phot_obsidForwrite))
    tds.addColumn("wavelength(um)",Column(wave))
    tds.addColumn("flux(Jy)",Column(phot))
    tds.addColumn("uncertainty(Jy)",Column(error))
    tds.addColumn("aperture(arcsec)",Column(phot_aper))
    asciiTableWriter(file=outdir+obj_list[Obsid.index(obsidFTS)]+"_spire_phot.txt",table=tds, writeMetadata=False)
    ############################### End of script ##################################
    ################################################################################
