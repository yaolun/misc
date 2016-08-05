def hyperion_image(rtout, wave, plotdir, printname, dstar=200., group=0, marker=0,
                    size='full', convolve=False, unit=None, scalebar=None):
    # to avoid X server error
    import matplotlib as mpl
    mpl.use('Agg')
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import astropy.constants as const
    from hyperion.model import ModelOutput
    # Package for matching the colorbar
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    pc = const.pc.cgs.value

    if unit == None:
        unit = 'erg\,s^{-1}\,cm^{-2}\,Hz^{-1}\,sr^{-1}'

    m = ModelOutput(rtout)

    # Extract the image.
    image = m.get_image(group=group, inclination=0, distance=dstar * pc, units='MJy/sr')

    # print np.shape(image.val)
    # Open figure and create axes
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    # Find the closest wavelength
    iwav = np.argmin(np.abs(wave - image.wav))

    # Calculate the image width in arcseconds given the distance used above
    # get the max radius
    rmax = max(m.get_quantities().r_wall)
    w = np.degrees(rmax / image.distance) * 3600.

    # Image in the unit of MJy/sr
    # Change it into erg/s/cm2/Hz/sr
    # factor = 1e-23*1e6
    factor = 1
    # avoid zero in log
    # flip the image, because the setup of inclination is upside down
    val = image.val[::-1, :, iwav] * factor + 1e-30

    if convolve:
        from astropy.convolution import convolve, Gaussian2DKernel
        img_res = 2*w/len(val[:,0])
        kernel = Gaussian2DKernel(0.27/2.354/img_res)
        val = convolve(val, kernel)

    if size != 'full':
        pix_e2c = (w-size/2.)/w * len(val[:,0])/2
        val = val[pix_e2c:-pix_e2c, pix_e2c:-pix_e2c]
        w = size/2.

    # This is the command to show the image. The parameters vmin and vmax are
    # the min and max levels for the colorscale (remove for default values).
    # cmap = sns.cubehelix_palette(start=0.1, rot=-0.7, gamma=0.2, as_cmap=True)
    cmap = plt.cm.CMRmap
    im = ax.imshow(val,
            norm=mpl.colors.LogNorm(vmin=1.515e-01, vmax=4.118e+01),
            cmap=cmap, origin='lower', extent=[-w, w, -w, w], aspect=1)

    # draw the flux extraction regions
    # x = 100
    # y = 100
    # area = x*y / 4.25e10
    # offset = 50
    #
    # pos_n = (len(val[0,:])/2.-1,len(val[0,:])/2.-1 + offset*len(val[0,:])/2/w)
    # pos_s = (len(val[0,:])/2.-1,len(val[0,:])/2.-1 - offset*len(val[0,:])/2/w)
    #
    # import matplotlib.patches as patches
    # ax.add_patch(patches.Rectangle((-x/2, -y), x, y, fill=False, edgecolor='lime'))
    # ax.add_patch(patches.Rectangle((-x/2, 0), x, y, fill=False, edgecolor='lime'))

    # plot the marker for center position by default or user input offset
    ax.plot([0],[-marker], '+', color='lime', markersize=10, mew=2)
    ax.set_xlim([-w,w])
    ax.set_ylim([-w,w])
    # ax.plot([0],[-10], '+', color='m', markersize=10, mew=2)

    # plot scalebar
    if scalebar != None:
        ax.plot([0.85*w-scalebar, 0.85*w], [-0.8*w, -0.8*w], color='w', linewidth=3)
        # add text
        ax.text(0.85*w-scalebar/2, -0.9*w, r'$\rm{'+str(scalebar)+"\,arcsec}$",
                color='w', fontsize=18, fontweight='bold', ha='center')

    # fix the tick label font
    ticks_font = mpl.font_manager.FontProperties(family='STIXGeneral',size=16)
    for label in ax.get_xticklabels():
        label.set_fontproperties(ticks_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(ticks_font)

    # Colorbar setting
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = fig.colorbar(im, cax=cax)
    cb.solids.set_edgecolor("face")
    cb.ax.minorticks_on()
    cb.ax.set_ylabel(r'$\rm{Intensity\,['+unit+']}$',fontsize=16)
    cb.ax.tick_params('both', width=1.5, which='major', length=3)
    cb.ax.tick_params('both', width=1.5, which='minor', length=2)
    cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
    plt.setp(cb_obj,fontsize=18)
    # fix the tick label font
    ticks_font = mpl.font_manager.FontProperties(family='STIXGeneral',size=18)
    for label in cb.ax.get_yticklabels():
        label.set_fontproperties(ticks_font)

    ax.set_xlabel(r'$\rm{RA\,Offset\,[arcsec]}$', fontsize=16)
    ax.set_ylabel(r'$\rm{Dec\,Offset\,[arcsec]}$', fontsize=16)

    # set the frame color
    ax.spines['bottom'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.spines['right'].set_color('white')

    ax.tick_params(axis='both', which='major', width=1.5, labelsize=18, color='white', length=5)
    ax.text(0.7,0.88,str(wave) + r'$\rm{\,\mu m}$',fontsize=20,color='white', transform=ax.transAxes)

    fig.savefig(plotdir+printname+'_image_'+str(wave)+'.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()

rtout = '/Volumes/SD-Mac/model2.rtout'
wave = 3.6
plotdir = '/Users/yaolun/bhr71/hyperion/controlled/'
hyperion_image(rtout, wave, plotdir, 'model2', group=0, dstar=200., marker=0, size=100, convolve=False,
        unit='MJy\,sr^{-1}', scalebar=30)
