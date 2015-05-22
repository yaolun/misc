def alma_cavity(freq, outdir, vlim, units='MJy/sr', pix=300, filename=None, label=None):
    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.constants as const
    from hyperion.model import ModelOutput
    from matplotlib.ticker import MaxNLocator

    # constants setup
    c = const.c.cgs.value
    pc = const.pc.cgs.value
    au = const.au.cgs.value
    # Image in the unit of MJy/sr
    # Change it into erg/s/cm2/Hz/sr
    if units == 'erg/s/cm2/Hz/sr':
        factor = 1e-23*1e6
        cb_label = r'$\rm{I_{\nu}\,(erg\,s^{-1}\,cm^{-2}\,Hz^{-1}\,sr^{-1})}$'
    elif units == 'MJy/sr':
        factor = 1
        cb_label = r'$\rm{I_{\nu}\,(MJy\,sr^{-1})}$'

    if filename == None:
        # input files setup
        filename_reg = '/Users/yaolun/test/model12.rtout'
        filename_r2  = '/Users/yaolun/test/model13.rtout'
        filename_r15 = '/Users/yaolun/test/model17.rtout'
        filename_uni = '/Users/yaolun/test/model62.rtout'
    else:
        filename_reg = filename['reg']
        filename_r2  = filename['r2']
        filename_r15 = filename['r15']
        filename_uni = filename['uni']

    if label == None:
        label_reg = r'$\rm{const.+r^{-2}}$'
        label_r2  = r'$\rm{r^{-2}}$'
        label_r15 = r'$\rm{r^{-1.5}}$'
        label_uni = r'$\rm{uniform}$'
    else:
        label_reg = label['reg']
        label_r2  = label['r2']
        label_r15 = label['r15']
        label_uni = label['uni']

    wl_aper = [3.6, 4.5, 5.8, 8.0, 8.5, 9, 9.7, 10, 10.5, 11, 16, 20, 24, 35, 70, 100, 160, 250, 350, 500, 850]
    wav = c/freq/1e9*1e4
    # wav = 40

    # read in
    # regular cavity setting
    m_reg = ModelOutput(filename_reg)
    image_reg = m_reg.get_image(group=len(wl_aper)+1, inclination=0, distance=178.0*pc, units='MJy/sr')
    # Calculate the image width in arcseconds given the distance used above
    w = np.degrees((1.5 * pc) / image_reg.distance) * 60.
    pix_num = len(image_reg.val[:,0,0])
    pix2arcsec = 2*w/pix_num
    pix2au = np.radians(2*w/pix_num/3600.)*image_reg.distance/au

    iwav = np.argmin(np.abs(wav - image_reg.wav))
    # avoid zero in log
    val_reg = image_reg.val[:, :, iwav] * factor + 1e-30

    # r^-2 cavity setting
    m_r2 = ModelOutput(filename_r2)
    image_r2 = m_r2.get_image(group=len(wl_aper)+1, inclination=0, distance=178.0*pc, units='MJy/sr')
    # Calculate the image width in arcseconds given the distance used above
    w = np.degrees((1.5 * pc) / image_r2.distance) * 60.
    pix_num = len(image_reg.val[:,0,0])
    pix2arcsec = 2*w/pix_num
    pix2au = np.radians(2*w/pix_num/3600.)*image_reg.distance/au
    iwav = np.argmin(np.abs(wav - image_r2.wav))
    # avoid zero in log
    val_r2 = image_r2.val[:, :, iwav] * factor + 1e-30

    # r^-1.5 cavity setting
    m_r15 = ModelOutput(filename_r15)
    image_r15 = m_r15.get_image(group=len(wl_aper)+1, inclination=0, distance=178.0*pc, units='MJy/sr')
    # Calculate the image width in arcseconds given the distance used above
    w = np.degrees((1.5 * pc) / image_r15.distance) * 60.
    pix_num = len(image_reg.val[:,0,0])
    pix2arcsec = 2*w/pix_num
    pix2au = np.radians(2*w/pix_num/3600.)*image_reg.distance/au
    iwav = np.argmin(np.abs(wav - image_r15.wav))
    # avoid zero in log
    val_r15 = image_r15.val[:, :, iwav] * factor + 1e-30

    # uniform cavity setting
    m_uni = ModelOutput(filename_uni)
    image_uni = m_uni.get_image(group=len(wl_aper)+1, inclination=0, distance=178.0*pc, units='MJy/sr')
    # Calculate the image width in arcseconds given the distance used above
    w = np.degrees((1.5 * pc) / image_uni.distance) * 60.
    pix_num = len(image_reg.val[:,0,0])
    pix2arcsec = 2*w/pix_num
    pix2au = np.radians(2*w/pix_num/3600.)*image_reg.distance/au
    iwav = np.argmin(np.abs(wav - image_uni.wav))
    # avoid zero in log
    val_uni = image_uni.val[:, :, iwav] * factor + 1e-30

    # 1-D radial intensity profile
    # get y=0 plane, and plot it
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)

    reg, = ax.plot(np.linspace(-pix/2,pix/2,num=pix)*pix2arcsec, val_reg[:,pix/2-1], color='b', linewidth=2)
    r2,  = ax.plot(np.linspace(-pix/2,pix/2,num=pix)*pix2arcsec, val_r2[:,pix/2-1], color='r', linewidth=1.5)
    r15, = ax.plot(np.linspace(-pix/2,pix/2,num=pix)*pix2arcsec, val_r15[:,pix/2-1], '--', color='r', linewidth=1.5)
    uni, = ax.plot(np.linspace(-pix/2,pix/2,num=pix)*pix2arcsec, val_uni[:,pix/2-1], color='k', linewidth=1.5)
    ax.legend([reg, r2, r15, uni], [label_reg, label_r2, label_r15, label_uni],\
              numpoints=1, loc='lower center', fontsize=18)

    ax.set_xlim([-1,1])
    ax.set_xlabel(r'$\rm{offset\,(arcsec)}$', fontsize=24)
    # ax.set_ylabel(r'$\rm{I_{\nu}~(erg~s^{-1}~cm^{-2}~Hz^{-1}~sr^{-1})}$', fontsize=16)
    ax.set_ylabel(cb_label, fontsize=24)

    [ax.spines[axis].set_linewidth(2) for axis in ['top','bottom','left','right']]
    ax.minorticks_on()
    ax.tick_params('both',labelsize=16,width=2,which='major',pad=10,length=5)
    ax.tick_params('both',labelsize=16,width=2,which='minor',pad=10,length=2.5)

    fig.savefig(outdir+'cavity_intensity_'+str(freq)+'.pdf', format='pdf', dpi=300, bbox_inches='tight')

    # 2-D intensity map
    from mpl_toolkits.axes_grid1 import AxesGrid
    image_grid = [val_reg, val_uni,val_r2, val_r15]
    label_grid = [label_reg, label_r2, label_r15, label_uni]
    fig = plt.figure(figsize=(30,30))
    grid = AxesGrid(fig, 142, # similar to subplot(142)
                        nrows_ncols = (2, 2),
                        axes_pad = 0,
                        share_all=True,
                        label_mode = "L",
                        cbar_location = "right",
                        cbar_mode="single",
                        )
    for i in range(4):
        offset = np.linspace(-pix/2,pix/2,num=pix)*pix2arcsec
        trim = np.where(abs(offset)<=2)
        im = grid[i].pcolor(np.linspace(-pix/2,pix/2,num=pix)*pix2arcsec, np.linspace(-pix/2,pix/2,num=pix)*pix2arcsec,\
            image_grid[i], cmap=plt.cm.jet, vmin=vlim[0], vmax=vlim[1])#vmin=(image_grid[i][trim,trim]).min(), vmax=(image_grid[i][trim,trim]).max())
        grid[i].set_xlim([-4,4])
        grid[i].set_ylim([-4,4])
        grid[i].set_xlabel(r'$\rm{RA\,offset\,(arcsec)}$', fontsize=14)
        grid[i].set_ylabel(r'$\rm{Dec\,offset\,(arcsec)}$', fontsize=14)
        # lg = grid[i].legend([label_grid[i]], loc='upper center', numpoints=1, fontsize=16)
        # for text in lg.get_texts():
        #     text.set_color('w')
        grid[i].text(0.5,0.8, label_grid[i], color='w', weight='heavy', fontsize=18, transform=grid[i].transAxes, ha='center')
        grid[i].locator_params(axis='x', nbins=5)
        grid[i].locator_params(axis='y', nbins=5)
        [grid[i].spines[axis].set_linewidth(1.2) for axis in ['top','bottom','left','right']]
        grid[i].tick_params('both',labelsize=12,width=1.2,which='major',pad=10,color='white',length=5)
        grid[i].tick_params('both',labelsize=12,width=1.2,which='minor',pad=10,color='white',length=2.5)

        # fix the overlap tick labels
        if i != 0:
            x_nbins = len(grid[i].get_xticklabels())
            y_nbins = len(grid[i].get_yticklabels())
            grid[i].yaxis.set_major_locator(MaxNLocator(nbins=5, prune='upper'))
            if i != 2:
                grid[i].xaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))

    [grid[0].spines[axis].set_color('white') for axis in ['bottom','right']]
    [grid[1].spines[axis].set_color('white') for axis in ['bottom','left']]
    [grid[2].spines[axis].set_color('white') for axis in ['top','right']]
    [grid[3].spines[axis].set_color('white') for axis in ['top','left']]

    #     ax.set_aspect('equal')
    cb = grid.cbar_axes[0].colorbar(im)
    cb.solids.set_edgecolor("face")
    cb.ax.minorticks_on()
    cb.ax.set_ylabel(cb_label,fontsize=12)
    cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
    plt.setp(cb_obj,fontsize=12)

    # fig.text(0.5, -0.05 , r'$\rm{RA~offset~(arcsec)}$', fontsize=12, ha='center')
    # fig.text(0, 0.5, r'$\rm{Dec~offset~(arcsec)}$', fontsize=12, va='center', rotation='vertical')

    fig.savefig(outdir+'cavity_2d_intensity_'+str(freq)+'.png',format='png',dpi=300,bbox_inches='tight')

filename = {'reg': '/Users/yaolun/bhr71/hyperion/alma/model2.rtout', \
            'r2': '/Users/yaolun/bhr71/hyperion/alma/model10.rtout', \
            'r15': '/Users/yaolun/bhr71/hyperion/alma/model13.rtout', \
            'uni': '/Users/yaolun/bhr71/hyperion/alma/model16.rtout'}
filename = {'reg': '/Users/yaolun/test/model91_old_aper.rtout',\
            'r2': '/Users/yaolun/test/model91_old_aper.rtout',\
            'r15': '/Users/yaolun/test/model91.rtout',\
            'uni': '/Users/yaolun/test/model91.rtout'}
# filename = {'reg': '/home/bettyjo/yaolun/hyperion/bhr71/alma/model2/model2.rtout', \
#             'r2': '/home/bettyjo/yaolun/hyperion/bhr71/alma/model10/model10.rtout', \
#             'r15': '/home/bettyjo/yaolun/hyperion/bhr71/alma/model13/model13.rtout', \
#             'uni': '/home/bettyjo/yaolun/hyperion/bhr71/alma/model16/model16.rtout'}
label = {'reg': r'$\rm{const.+r^{-2}}$', 'r2': r'$\rm{r^{-2}}$', 'r15': r'$\rm{r^{-1.5}}$', 'uni': r'$\rm{uniform}$'}
freq = [230,345,460]
freq = [83275.682]
# freq = [345]
vlim = [[150,800],[600,3200],[1500,9000]]
vlim = [[0,100]]
# vlim = [[600,2600]]
outdir = '/Users/yaolun/test/'
# outdir = '/home/bettyjo/yaolun/test/'
for f in freq:
    print 'Frequency: %5f GHz' % f
    # if f == 345:
    alma_cavity(f,outdir, vlim[freq.index(f)], filename=filename, label=label, pix=300)