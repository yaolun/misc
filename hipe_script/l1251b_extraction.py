slwSpectrum2d = obs.level2.getProduct('HR_SLW_spectrum2d_apod')
sswSpectrum2d = obs.level2.getProduct('HR_SSW_spectrum2d_apod')
# RA and Dec of L1251B
coord = [339.6949583,75.1923889]
pointSpectrum = specPointSourceExtractor(slwSpectra=slwSpectrum2d,
                                         sswSpectra=sswSpectrum2d,
                                         coords=coord, calibration=obs.calibration,
                                         maxDistArcsec=31.8/2, averageSpectra=False)
correctedSpectrum = semiExtendedCorrector(spectrum=pointSpectrum,
                                          sourceModel=SemiExtendedSourceModelShape("gaussian", 18.0, 0.0, 0.0, 0.0, 0.0, 1.0, 257, 257),
                                          calibration=obs.calibration, doPlots=False,
                                          optimiseDiameter=True)
fitted_size = correctedSpectrum.meta['sourceDiameter'].double
correctedSpectrum = semiExtendedCorrector(spectrum=pointSpectrum,
                                          sourceModel=SemiExtendedSourceModelShape("gaussian", 18.0, 0.0, 0.0, 0.0, 0.0, 1.0, 257, 257),
                                          calibration=obs.calibration, optimiseDiameter=True,
                                          gaussRefBeamDiam=fitted_size)
simpleFitsWriter(pointSpectrum, "L1251B_PointSourceExtraction.fits")
simpleFitsWriter(correctedSpectrum, "L1251B_PointSourceExtraction_SECT.fits")
