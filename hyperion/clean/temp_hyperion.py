def temp_hyperion(rtout,outdir, bb_dust=False):
    import numpy as np
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import os
    from hyperion.model import ModelOutput
    import astropy.constants as const
    from matplotlib.colors import LogNorm

    # seaborn colormap
    import seaborn.apionly as sns

    # constants setup
    AU = const.au.cgs.value

    # misc variable setup
    print_name = os.path.splitext(os.path.basename(rtout))[0]

    m = ModelOutput(rtout)
    q = m.get_quantities()

    # get the grid info
    ri, thetai = q.r_wall, q.t_wall
    rc     = 0.5*( ri[0:len(ri)-1]     + ri[1:len(ri)] )
    thetac = 0.5*( thetai[0:len(thetai)-1] + thetai[1:len(thetai)] )

    # get the temperature profile
    # and average across azimuthal angle
    # temperature array in [phi, theta, r]
    temp = q['temperature'][0].array.T
    temp2d = np.sum(temp**2, axis=2)/np.sum(temp, axis=2)
    temp2d_exp = np.hstack((temp2d,temp2d,temp2d[:,0:1]))
    thetac_exp = np.hstack((thetac-np.pi/2, thetac+np.pi/2, thetac[0]-np.pi/2))

    mag = 1
    fig = plt.figure(figsize=(mag*8,mag*6))
    ax = fig.add_subplot(111, projection='polar')

    # cmap = sns.cubehelix_palette(light=1, as_cmap=True)
    cmap = plt.cm.CMRmap
    im = ax.pcolormesh(thetac_exp, rc/AU, temp2d_exp, cmap=cmap, norm=LogNorm(vmin=5, vmax=100))
    im.set_edgecolor('face')

    ax.set_xlabel(r'$\rm{Polar\,angle\,(Degree)}$',fontsize=20)
    ax.set_ylabel(r'$\rm{Radius\,(AU)}$',fontsize=20, labelpad=-140, color='grey')
    ax.tick_params(labelsize=16)
    ax.tick_params(axis='y', colors='grey')
    ax.set_yticks(np.arange(0,int(R_env_max/AU/10000.)*10000, 10000))

    ax.set_xticklabels([r'$\rm{90^{\circ}}$',r'$\rm{45^{\circ}}$',r'$\rm{0^{\circ}}$',r'$\rm{-45^{\circ}}$',\
                            r'$\rm{-90^{\circ}}$',r'$\rm{-135^{\circ}}$',r'$\rm{180^{\circ}}$',r'$\rm{135^{\circ}}$'])
    cb = fig.colorbar(im, pad=0.1)
    cb.ax.set_ylabel(r'$\rm{Averaged\,Temperature\,(K)}$',fontsize=20)
    cb.set_ticks([5,10,20,30,40,50,60,70,80,90,100])
    cb.set_ticklabels([r'$\rm{5}$',r'$\rm{10}$',r'$\rm{20}$',r'$\rm{30}$',r'$\rm{40}$',r'$\rm{50}$',r'$\rm{60}$',r'$\rm{70}$',r'$\rm{80}$',r'$\rm{90}$',r'$\rm{>100}$'])
    cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
    plt.setp(cb_obj,fontsize=20)

    # fix the tick label font
    ticks_font = mpl.font_manager.FontProperties(family='STIXGeneral',size=20)
    for label in ax.get_yticklabels():
        label.set_fontproperties(ticks_font)

    fig.savefig(outdir+print_name+'_temperature.png', format='png', dpi=300, bbox_inches='tight')
    fig.clf()

    # Plot the radial temperature profile
    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(111)

    plot_grid = [0,99,199]
    label_grid = [r'$\rm{outflow}$', r'$\rm{45^{\circ}}$', r'$\rm{midplane}$']
    alpha = np.linspace(0.3,1.0,len(plot_grid))
    color_list = [[0.8507598215729224, 0.6322174528970308, 0.6702243543099417],\
                  [0.5687505862870377, 0.3322661256969763, 0.516976691731939],\
                  [0.1750865648952205, 0.11840023306916837, 0.24215989137836502]]

    for i in plot_grid:
        temp_rad, = ax.plot(np.log10(rc/AU), np.log10(temp2d[:,i]),'-',color=color_list[plot_grid.index(i)],\
                            linewidth=2, markersize=3,label=label_grid[plot_grid.index(i)])

    # plot the theoretical prediction for black body dust without considering the extinction
    if bb_dust == True:
        from hyperion.model import Model
        sigma = const.sigma_sb.cgs.value
        lsun = const.L_sun.cgs.value

        dum = Model()
        dum.use_sources(rtout)
        L_cen = dum.sources[0].luminosity/lsun

        t_bbdust = (L_cen*lsun/(16*np.pi*sigma*rc**2))**(0.25)
        temp_bbdust, = ax.plot(np.log10(rc/AU), np.log10(t_bbdust), '--', color='r', linewidth=2.5,label=r'$\rm{blackbody\,dust}$')

    ax.legend(loc='upper right', numpoints=1, fontsize=24)
    ax.set_xlabel(r'$\rm{log\,R\,(AU)}$',fontsize=24)
    ax.set_ylabel(r'$\rm{log\,T\,(K)}$',fontsize=24)
    [ax.spines[axis].set_linewidth(2) for axis in ['top','bottom','left','right']]
    ax.minorticks_on()
    ax.tick_params('both',labelsize=24,width=2,which='major',pad=15,length=5)
    ax.tick_params('both',labelsize=24,width=2,which='minor',pad=15,length=2.5)

    # fix the tick label font
    ticks_font = mpl.font_manager.FontProperties(family='STIXGeneral',size=24)
    for label in ax.get_xticklabels():
        label.set_fontproperties(ticks_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(ticks_font)

    ax.set_ylim([0,4])
    fig.gca().set_xlim(left=np.log10(0.05))
    # ax.set_xlim([np.log10(0.8),np.log10(10000)])

    fig.savefig(outdir+print_name+'_temp_radial.pdf',format='pdf',dpi=300,bbox_inches='tight')
    fig.clf()

# rtout = '/Users/yaolun/bhr71/hyperion/controlled/model224.rtout'
# outdir = '/Users/yaolun/bhr71/hyperion/controlled/'
# temp_hyperion(rtout, outdir, bb_dust=True)
