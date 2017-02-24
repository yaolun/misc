obsid_phot = {'Source': ['L1157', 'L1014', 'IRAS03301+3111', 'B1-c', 'B1-a',
                         'IRAS03245+3002', 'L1455-IRS3', 'TMR1', 'TMC1', 'TMC1A',
                         'L1551-IRS5', 'B335', 'GSS30-IRS1', 'VLA1623-243', 'WL12',
                         'IRS46', 'RCrA-IRS7C', 'RCrA-IRS5A', 'RCrA-IRS7B', 'BHR71',
                         'DKCha', 'Elias29', 'IRAM04191+1522', 'IRS63', 'L1448-MM', 'L1489',
                         'L1527', 'Serpens-SMM3', 'Serpens-SMM4', 'Ced110-IRS4', 'HH100',
                         'HH46', 'IRAS15398-3359', 'L483', 'L723-MM', 'RNO91'],
              'obsid': [[1342224778, 1342224779, 1342189845],
                        [1342225450, 1342225449],
                        [1342227103, 1342227104],
                        [1342267246, 1342267247],
                        [1342227103, 1342227104],
                        [1342227103, 1342227104],
                        [1342227103, 1342227104],
                        [1342228175, 1342228174],
                        [1342202252],
                        [1342202252],
                        [1342202251],
                        [1342196030, 1342196031],
                        [1342227148, 1342227149, 1342205093, 1342205094],
                        [1342205093, 1342205094],
                        [1342238817, 1342238816],
                        [1342205093, 1342205094],
                        [1342184510, 1342184511],
                        [1342267429, 1342267427, 1342267428, 1342267430, 1342242076, 1342241402,
                         1342241519, 1342241403, 1342242555, 1342242077, 1342241314, 1342241520,
                         1342241313, 1342242554],
                        [1342218806, 1342218807],
                        [1342224922, 1342224925, 1342224924, 1342224923],
                        [1342212709, 1342212708, 1342213180],
                        [1342205093, 1342205094],
                        [1342190941, 1342190942, 1342241875, 1342241874],
                        [1342227037, 1342227038],
                        [],
                        [1342216036, 1342216037],
                        [1342243453, 1342243455, 1342243454],
                        [],
                        [1342229079, 1342229080],
                        [1342223480, 1342223481, 1342224782, 1342224783],
                        [1342267769, 1342267768],
                        [],
                        [1342226706, 1342226705],
                        [1342228398, 1342228397, 1342228395, 1342228396],
                        [1342231917, 1342231918],
                        [1342263844, 1342263845]],
                 'ra': ['20:39:06.3','21:24:07.5','03:33:12.8','03:33:17.9','03:33:16.7',
                        '03:27:39.1','03:28:00.4','04:39:13.9','04:41:12.7','04:39:35',
                        '04:31:34.1','19:37:00.9','16:26:21.4','16:26:26.4','16:26:44.2',
                        '16:27:29.4','19:01:55.3','19:01:48.1','19:01:56.4','12:01:36.3',
                        '12:53:17.2','16:27:09.4','04:21:56.9','16:31:35.6','03:25:38.9',
                        '04:04:42.9','04:39:53.9','18:29:59.3','18:29:56.7','11:06:47',
                        '19:01:49.1','08:25:43.9','15:42:01.3','18:17:29.9','19:17:53.7',
                        '16:34:29.3'],
                'dec': ['+68:02:16','+49:59:09','+31:21:24.2','+31:09:31.9','+31:07:55.2',
                        '+30:13:03.1','+30:08:01.3','+25:53:20.6','+25:46:35.9','+25:41:45.5',
                        '+18:08:04.9','+07:34:09.7','-24:23:04.3','-24:24:30','-24:34:48.4',
                        '-24:39:16.1','-36:57:17','-36:57:22.7','-36:57:28.3','-65:08:53',
                        '-77:07:10.7','-24:37:18.6','+15:29:45.9','-24:01:29.3','+30:44:05.4',
                        '+26:18:56.3','+26:03:09.8','+01:14:01.7','+01:13:17.2','-77:22:32.4',
                        '-36:58:16','-51:00:36','-34:09:15','-04:39:39.5','+19:12:20',
                        '-15:47:01.4']
}

outdir = '/home/bettyjo/yaolun/CDF_archive_v2/'
aper_data = asciiTableReader(file=outdir+'pacs_1d_apertures.txt', tableType='SPACES')

import sys

# aperture size in radius
# aper = 31.8/2

start_from = 'L1157'
skip = True

# PACS aperture photometry
for i in range(len(obsid_phot['Source'])):
    if obsid_phot['Source'][i] == start_from:
        skip = False
    if skip:
        print obsid_phot['Source'][i], ' is skipped.'
        continue
    obsid_dum = obsid_phot['obsid'][i]
    if len(obsid_dum) == 0:
        print 'No photometry obsid found for ', obsid_phot['Source'][i]
        continue

    if not obsid_phot['Source'][i] in aper_data['Object'].data:
        print obsid_phot['Source'][i], ' not found in the processed object list.'
        continue

    radius = aper_data['aperture'].data[list(aper_data['Object'].data).index(obsid_phot['Source'][i])] / 2

    if radius > 26.5:
        radius = 26.5

    wave = []
    flux = []
    error = []

    # get aperture fitted by matching PACS and SPIRE spectra


    for obsid in obsid_dum:
    	obs = getObservation(obsid=obsid,useHsa=True)
    	caltree=getCalTree(obs=obs)

        innerArcsec = 150.0
        outerArcsec = 200.0
        if obsid_phot['Source'][i] == 'RCrA-IRS7C':
            innerArcsec = 40.0
            outerArcsec = 60.0

	if obs.meta['instrument'].string == 'SPIRE':
	    print obsid, ' is SPIRE'
	    continue

    	blue_cam = obs.refs['level2'].product.refs['HPPPMAPB'].product
    	red_cam = obs.refs['level2'].product.refs['HPPPMAPR'].product

        blue_wave = blue_cam.meta['wavelength'].double
        red_wave =  red_cam.meta['wavelength'].double

    	blue_phot = annularSkyAperturePhotometry(image=blue_cam, fractional=1,\
    				centerRa=obsid_phot['ra'][i], centerDec=obsid_phot['dec'][i], radiusArcsec=radius,
    				innerArcsec=innerArcsec, outerArcsec=outerArcsec)
    	red_phot = annularSkyAperturePhotometry(image=red_cam, fractional=1,\
    				centerRa=obsid_phot['ra'][i], centerDec=obsid_phot['dec'][i], radiusArcsec=radius,
    				innerArcsec=innerArcsec, outerArcsec=outerArcsec)

    	blue_phot_cor = photApertureCorrectionPointSource(apphot=blue_phot, band="blue", \
    				calTree=caltree, responsivityVersion=6)

     	red_phot_cor = photApertureCorrectionPointSource(apphot=red_phot, band="red", \
    				calTree=caltree, responsivityVersion=6)

        print obsid_phot['Source'][i]
        print blue_phot_cor.total[2], blue_phot_cor.error[2]
    	print red_phot_cor.total[2], red_phot_cor.error[2]

        wave.append(blue_wave)
        wave.append(red_wave)
        flux.append(blue_phot_cor.total[2])
        flux.append(red_phot_cor.total[2])
        error.append(blue_phot_cor.error[2])
        error.append(red_phot_cor.error[2])

    # Write to ASCII table
    tds = TableDataset()
    wave = Float1d(wave)
    flux = Float1d(flux)
    err =  Float1d(error)
    tds.addColumn("wavelength(um)",Column(wave))
    tds.addColumn("flux(Jy)",Column(flux))
    tds.addColumn("uncertainty(Jy)",Column(err))
    asciiTableWriter(file=outdir+obsid_phot['Source'][i]+'/pacs/data/'+obsid_phot['Source'][i]+'_pacs_phot.txt',
                     table=tds, writeMetadata=False)
