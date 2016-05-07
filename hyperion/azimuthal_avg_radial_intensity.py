def azimuthal_avg_radial_intensity(wave, imgpath, source_center, rtout, plotname,
                                   annulus_width=10, group=8, dstar=200.):

    import numpy as np
    import matplotlib as mpl
    from astropy.io import ascii, fits
    import matplotlib.pyplot as plt
    from photutils import aperture_photometry as ap
    from photutils import CircularAperture, CircularAnnulus
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astropy import wcs
    from hyperion.model import ModelOutput
    import astropy.constants as const

    pc = const.pc.cgs.value
    AU = const.au.cgs.value

    # source_center = '12 01 36.3 -65 08 53.0'

    # Read in data and set up coversions
    im_hdu = fits.open(imgpath)
    im = im_hdu[1].data
    w = wcs.WCS(im_hdu[1].header)

    coord = SkyCoord(source_center, unit=(u.hourangle, u.deg))
    pixcoord = w.wcs_world2pix(coord.ra.degree, coord.dec.degree, 1)
    pix2arcsec = abs(im_hdu[1].header['CDELT1'])*3600.
    # convert intensity unit from MJy/sr to Jy/pixel
    factor = 1e6/4.25e10*abs(im_hdu[1].header['CDELT1']*im_hdu[1].header['CDELT2'])*3600**2

    # radial grid in arcsec
    # annulus_width = 10
    r = np.arange(10, 200, annulus_width, dtype=float)
    I = np.empty_like(r[:-1])

    # iteration
    for ir in range(len(r)-1):
        aperture = CircularAnnulus((pixcoord[0],pixcoord[1]), r_in=r[ir]/pix2arcsec, r_out=r[ir+1]/pix2arcsec)
    #     print aperture.r_in
        phot = ap(im, aperture)
        I[ir] = phot['aperture_sum'].data * factor / aperture.area()
        # print r[ir], I[ir]

    # read in from RTout
    rtout = ModelOutput(rtout)
    # setting up parameters
    # dstar = 200.
    # group = 8
    # wave = 500.0

    im = rtout.get_image(group=group, inclination=0, distance=dstar*pc, units='Jy')

    # Find the closest wavelength
    iwav = np.argmin(np.abs(wave - im.wav))
    # avoid zero when log, and flip the image
    val = im.val[::-1, :, iwav]

    w = np.degrees(max(rtout.get_quantities().r_wall) / im.distance) * 3600
    npix = len(val[:,0])
    pix2arcsec = 2*w/npix

    # radial grid in arcsec
    # annulus_width = 10
    r = np.arange(10, 200, annulus_width, dtype=float)
    I_sim = np.empty_like(r[:-1])

    # iteration
    for ir in range(len(r)-1):
        aperture = CircularAnnulus((npix/2.+0.5, npix/2.+0.5), r_in=r[ir]/pix2arcsec, r_out=r[ir+1]/pix2arcsec)
    #     print aperture.r_in
        phot = ap(val, aperture)
        I_sim[ir] = phot['aperture_sum'].data / aperture.area()
        # print r[ir], I_sim[ir]

    # plot
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)

    i_sim, = ax.plot(np.log10(r[:-1]*dstar), np.log10(I_sim/I_sim.max()), 'o-', mec='None', markersize=10)
    i, = ax.plot(np.log10(r[:-1]*dstar), np.log10(I/I.max()), 'o-', mec='None', markersize=10)

    ax.legend([i, i_sim], [r'$\rm{observation}$', r'$\rm{simulation}$'], fontsize=16, numpoints=1, loc='upper right')
    [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
    ax.minorticks_on()
    ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=10,length=5)
    ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=10,length=2.5)
    ax.set_xlabel('log(Radius) [AU]', fontsize=18)
    ax.set_ylabel(r'$\rm{log(I\,/\,I_{max})}$', fontsize=18)

    # fix the tick label font
    ticks_font = mpl.font_manager.FontProperties(family='STIXGeneral',size=18)
    for label in ax.get_xticklabels():
        label.set_fontproperties(ticks_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(ticks_font)


    fig.savefig(plotname+'_radial_profile_'+str(wave)+'um.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()
