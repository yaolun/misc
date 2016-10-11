def azimuthal_avg_radial_intensity(wave, rtout, plotname, dstar,
                                   annulus_width=10, rrange=[10,200], group=8, obs=None,
                                   other_obs=None):

    """
    The 'obs' option only works for Herschel PACS/SPIRE image.
    The 'obs' option now accept
    """

    import numpy as np
    import matplotlib as mpl
    # to avoid X server error
    mpl.use('Agg')
    from astropy.io import ascii, fits
    import matplotlib.pyplot as plt
    from photutils import aperture_photometry as ap
    from photutils import CircularAperture, CircularAnnulus
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astropy import wcs
    from hyperion.model import ModelOutput
    import astropy.constants as const
    import os

    pc = const.pc.cgs.value
    AU = const.au.cgs.value

    # radial grid in arcsec
    # make the annulus center on
    r = np.arange(rrange[0], rrange[1], annulus_width, dtype=float) - annulus_width*0.5
    # r = np.arange(rrange[0], rrange[1], annulus_width, dtype=float) - annulus_width*0.55

    # source_center = '12 01 36.3 -65 08 53.0'

    def ExtractIntensityObs(rrange, annulus_width, obs):
        import numpy as np
        from astropy.io import fits
        from astropy.coordinates import SkyCoord
        from astropy import wcs
        from photutils import aperture_photometry as ap
        from photutils import CircularAperture, CircularAnnulus

        r = np.arange(rrange[0], rrange[1], annulus_width, dtype=float) - annulus_width*0.5
        # r = np.arange(rrange[0], rrange[1], annulus_width, dtype=float) - annulus_width*0.55

        imgpath = obs['imgpath']
        source_center = obs['source_center']
        # Read in data and set up coversions
        im_hdu = fits.open(imgpath)
        im = im_hdu[1].data
        wave = im_hdu[0].header['WAVELNTH']
        # error
        if (wave < 200.0) & (wave > 70.0):
            im_err = im_hdu[5].data
        elif (wave > 200.0) & (wave < 670.0):
            im_err = im_hdu[2].data
        else:
            im_err_exten = raw_input('The extension that includes the image error: ')
            im_err = im_hdu[int(im_err_exten)].data

        w = wcs.WCS(im_hdu[1].header)

        coord = SkyCoord(source_center, unit=(u.hourangle, u.deg))
        pixcoord = w.wcs_world2pix(coord.ra.degree, coord.dec.degree, 1)
        pix2arcsec = abs(im_hdu[1].header['CDELT1'])*3600.

        # determine whether need to convert the unit
        factor = 1
        print 'Image unit is ', im_hdu[1].header['BUNIT']
        if im_hdu[1].header['BUNIT'] != 'Jy/pixel':
            print 'Image unit is ', im_hdu[1].header['BUNIT']

            if im_hdu[1].header['BUNIT'] == 'MJy/sr':
                # convert intensity unit from MJy/sr to Jy/pixel
                factor = 1e6/4.25e10*abs(im_hdu[1].header['CDELT1']*im_hdu[1].header['CDELT2'])*3600**2
            else:
                factor = raw_input('What is the conversion factor to Jy/pixel?')

        I = np.empty_like(r[:-1])
        I_low = np.empty_like(r[:-1])
        I_hi = np.empty_like(r[:-1])
        I_err = np.empty_like(r[:-1])

        # for calculating the uncertainty from the variation within each annulus
        # construct the x- and y-matrix
        grid_x, grid_y = np.meshgrid(np.linspace(0,len(im[0,:])-1,len(im[0,:])),
                                     np.linspace(0,len(im[:,0])-1,len(im[:,0])))

        grid_dist = ((grid_x-pixcoord[0])**2+(grid_y-pixcoord[1])**2)**0.5

        # iteration
        for ir in range(len(r)-1):
            aperture = CircularAnnulus((pixcoord[0],pixcoord[1]), r_in=r[ir]/pix2arcsec, r_out=r[ir+1]/pix2arcsec)
            phot = ap(im, aperture, error=im_err)
            I[ir] = phot['aperture_sum'].data * factor / aperture.area()

            # uncertainty
            im_dum = np.where((grid_dist < r[ir+1]/pix2arcsec) & (grid_dist >= r[ir]/pix2arcsec), im, np.nan)

            # estimate the uncertainty by offsetting the annulus by +/- 1 pixel
            offset = -1
            if r[ir]/pix2arcsec + offset < 0:
                offset = -r[ir]/pix2arcsec
            aperture = CircularAnnulus((pixcoord[0],pixcoord[1]),
                                r_in=r[ir]/pix2arcsec + offset, r_out=r[ir+1]/pix2arcsec + offset)
            phot = ap(im, aperture, error=im_err)
            I_low[ir] = phot['aperture_sum'].data * factor / aperture.area()

            offset = 1
            aperture = CircularAnnulus((pixcoord[0],pixcoord[1]),
                                r_in=r[ir]/pix2arcsec + offset, r_out=r[ir+1]/pix2arcsec + offset)
            phot = ap(im, aperture, error=im_err)
            I_hi[ir] = phot['aperture_sum'].data * factor / aperture.area()

        I_err = (abs(I_low - I) + abs(I_hi - I))/2.

        return r, I, I_err

    if obs != None:
        I_obs = []
        for o in obs:
            if 'label' not in o.keys():
                label_dum = r'$\rm{observation}$'
                color_dum = 'g'
                linestyle_dum = '-'
                rrange_dum = rrange
                annulus_width_dum = annulus_width
            else:
                label_dum = o['label']
                color_dum = o['plot_color']
                linestyle_dum = o['plot_linestyle']
                rrange_dum = o['rrange']
                annulus_width_dum = o['annulus_width']

            r_dum, I_dum, I_err_dum = ExtractIntensityObs(rrange_dum, annulus_width_dum, o)
            # determine the label
            I_obs.append({'imgpath':o['imgpath'], 'r':r_dum, 'I':I_dum, 'I_err':I_err_dum, 'label': label_dum,
                          'plot_color':color_dum, 'plot_linestyle':linestyle_dum})

        # The first image should be the one to be compared primarily, and written out
        I = I_obs[0]['I']
        I_err = I_obs[0]['I_err']
        imgpath = I_obs[0]['imgpath']
        #

    # read in from RTout
    rtout = ModelOutput(rtout)

    im = rtout.get_image(group=group, inclination=0, distance=dstar*pc, units='Jy', uncertainties=True)
    factor = 1

    # Find the closest wavelength
    iwav = np.argmin(np.abs(wave - im.wav))
    # avoid zero when log, and flip the image
    val = im.val[::-1, :, iwav]
    unc = im.unc[::-1, :, iwav]

    w = np.degrees(max(rtout.get_quantities().r_wall) / im.distance) * 3600
    npix = len(val[:,0])
    pix2arcsec = 2*w/npix

    I_sim = np.empty_like(r[:-1])
    I_sim_hi = np.empty_like(r[:-1])
    I_sim_low = np.empty_like(r[:-1])
    I_sim_err = np.empty_like(r[:-1])

    # for calculating the uncertainty from the variation within each annulus
    # construct the x- and y-matrix
    grid_x, grid_y = np.meshgrid(np.linspace(0,npix-1,npix),
                                 np.linspace(0,npix-1,npix))

    dist_x = abs(grid_x - ((npix-1)/2.))
    dist_y = abs(grid_y - ((npix-1)/2.))

    grid_dist = (dist_x**2+dist_y**2)**0.5

    # iteration
    for ir in range(len(r)-1):
        aperture = CircularAnnulus((npix/2.+0.5, npix/2.+0.5),
                            r_in=r[ir]/pix2arcsec, r_out=r[ir+1]/pix2arcsec)
        phot = ap(val, aperture, error=unc)
        I_sim[ir] = phot['aperture_sum'].data / aperture.area()

        # uncertainty
        im_dum = np.where((grid_dist < r[ir+1]/pix2arcsec) & (grid_dist >= r[ir]/pix2arcsec), val, np.nan)
        # I_sim_err[ir] = phot['aperture_sum_err'].data / aperture.area()
        # I_sim_err[ir] = (np.nanstd(im_dum)**2+phot['aperture_sum_err'].data**2)**0.5 * factor / aperture.area()

        offset = -1
        aperture = CircularAnnulus((npix/2.+0.5, npix/2.+0.5),
                            r_in=r[ir]/pix2arcsec + offset, r_out=r[ir+1]/pix2arcsec + offset)
        phot = ap(val, aperture, error=unc)
        I_sim_low[ir] = phot['aperture_sum'].data * factor / aperture.area()

        offset = 1
        aperture = CircularAnnulus((npix/2.+0.5, npix/2.+0.5),
                            r_in=r[ir]/pix2arcsec + offset, r_out=r[ir+1]/pix2arcsec + offset)
        phot = ap(val, aperture, error=unc)
        I_sim_hi[ir] = phot['aperture_sum'].data * factor / aperture.area()

    I_sim_err = (abs(I_sim_low - I_sim)+ abs(I_sim_hi - I_sim))/2.


    if obs != None:
        # write the numbers into file
        foo = open(plotname+'_radial_profile_'+str(wave)+'um.txt', 'w')
        # print some header info
        foo.write('# wavelength '+str(wave)+' um \n')
        foo.write('# image file '+os.path.basename(imgpath)+' \n')
        foo.write('# annulus width '+str(annulus_width)+' arcsec \n')
        # write profiles
        foo.write('r_in \t I \t I_err \t I_sim \t I_sim_err \n')
        foo.write('# [arcsec] \t [Jy/pixel] \t [Jy/pixel] \t [Jy/pixel] \t [Jy/pixel] \n')
        for i in range(len(I)):
            foo.write('%f \t %e \t %e \t %e \t %e \n' % (r[i]+annulus_width/2., I[i], I_err[i], I_sim[i], I_sim_err[i]))
        foo.close()
    else:
        # write the numbers into file
        foo = open(plotname+'_radial_profile_'+str(wave)+'um.txt', 'w')
        # print some header info
        foo.write('# wavelength '+str(wave)+' um \n')
        foo.write('# annulus width '+str(annulus_width)+' arcsec \n')
        # write profiles
        foo.write('r_in \t I_sim \t I_sim_err \n')
        foo.write('# [arcsec] \t [Jy/pixel] \t [Jy/pixel] \n')
        for i in range(len(I_sim)):
            foo.write('%f \t %e \t %e \n' % (r[i]+annulus_width/2., I_sim[i], I_sim_err[i]))
        foo.close()

    # plot
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)

    I_sim_hi = np.log10((I_sim+I_sim_err)/I_sim.max())-np.log10(I_sim/I_sim.max())
    I_sim_low = np.log10(I_sim/I_sim.max())-np.log10((I_sim-I_sim_err)/I_sim.max())
    i_sim = ax.errorbar(np.log10(r[:-1]*dstar), np.log10(I_sim/I_sim.max()), color='b',
                    yerr=(I_sim_low, I_sim_hi), marker='o', linestyle='-', mec='None', markersize=5,
                    ecolor='b', elinewidth=1.5, capthick=1.5, barsabove=True)

    if obs != None:
        plot_profile = []
        plot_label = []
        for o in I_obs:
            I_hi = np.log10((o['I']+o['I_err'])/o['I'].max())-np.log10(o['I']/o['I'].max())
            I_low = np.log10(o['I']/o['I'].max())-np.log10((o['I']-o['I_err'])/o['I'].max())
            i = ax.errorbar(np.log10(o['r'][:-1]*dstar), np.log10(o['I']/o['I'].max()), color=o['plot_color'],
                            yerr=(I_low, I_hi), marker='o', linestyle=o['plot_linestyle'], mec='None', markersize=5,
                            ecolor=o['plot_color'], elinewidth=1.5, capthick=1.5, barsabove=True)
            plot_profile.append(i)
            plot_label.append(o['label'])

        plot_profile.append(i_sim)
        plot_label.append(r'$\rm{simulation}$')
        ax.legend(plot_profile, plot_label,
                  fontsize=16, numpoints=1, loc='best')
    else:
        ax.legend([i_sim], [r'$\rm{simulation}$'], fontsize=16, numpoints=1, loc='best')

    # limit radius
    ax.axvline([np.log10(100*dstar)], color='k', linestyle='--', linewidth=1)
    #
    [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
    ax.minorticks_on()
    ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=10,length=5)
    ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=10,length=2.5)
    ax.set_xlabel(r'$\rm{log(\it{b})\,[\rm{AU}]}$', fontsize=18)
    ax.set_ylabel(r'$\rm{log(I\,/\,I_{max})}$', fontsize=18)

    # fix the tick label font
    ticks_font = mpl.font_manager.FontProperties(family='STIXGeneral',size=18)
    for label in ax.get_xticklabels():
        label.set_fontproperties(ticks_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(ticks_font)

    fig.savefig(plotname+'_radial_profile_'+str(wave)+'um.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()

# obs_azi = [{'imgpath': '/Users/yaolun/test/hpacs1342224922_20hpppmapr_00_1431606963820.fits',
#           'source_center': '12:01:36.81 -65:08:49.22'}]#,
#         #   {'imgpath': '/Users/yaolun/test/extdPLW_jypx.fits',
        #    'source_center': '12:01:36.81 -65:08:49.22',
        #    'label': r'$\rm{SPIRE\,500\,\mu m}$', 'plot_color':'k', 'plot_linestyle':'-',
        #    'rrange': [10,200], 'annulus_width': 10}]
# azimuthal_avg_radial_intensity(160.0, '/Volumes/SD-Mac/model32_nontsc.rtout',
#                                '/Users/yaolun/bhr71/hyperion/model32_nontsc', 200.0, obs=obs_azi)
