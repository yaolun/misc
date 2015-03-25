def temp_hyperion(rtout,outdir):
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    from hyperion.model import ModelOutput
    import astropy.constants as const
    from matplotlib.colors import LogNorm

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

    mag = 1.5
    fig = plt.figure(figsize=(mag*8,mag*6))
    ax = fig.add_subplot(111, projection='polar')

    cmap = 'jet'
    im = ax.pcolormesh(thetac_exp, rc/AU, temp2d_exp, cmap=cmap, norm=LogNorm(vmin=temp2d.min(), vmax=100))
    im.set_edgecolor('face')

    ax.set_xlabel(r'$\mathrm{Polar~angle~(Degree)}$',fontsize=20)
    ax.set_ylabel(r'$\mathrm{Radius~(AU)}$',fontsize=20)
    ax.tick_params(labelsize=20)
    ax.set_yticks(np.arange(0,max(ri)/AU,max(ri)/AU/5))

    ax.set_xticklabels([r'$\mathrm{90^{\circ}}$',r'$\mathrm{45^{\circ}}$',r'$\mathrm{0^{\circ}}$',r'$\mathrm{-45^{\circ}}$',\
                            r'$\mathrm{-90^{\circ}}$',r'$\mathrm{-135^{\circ}}$',r'$\mathrm{180^{\circ}}$',r'$\mathrm{135^{\circ}}$'])
    ax.grid(True)
    cb = fig.colorbar(im, pad=0.1)
    cb.ax.set_ylabel(r'$\mathrm{Averaged~Temperature~(K)}$',fontsize=20)
    cb.set_ticks([5,10,20,30,40,50,60,70,80,90,100])
    cb.set_ticklabels([5,10,20,30,40,50,60,70,80,90,'$>$100'])
    cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
    plt.setp(cb_obj,fontsize=20)
    fig.savefig(outdir+print_name+'_temperature.png', format='png', dpi=300, bbox_inches='tight')
    fig.clf()

    # Plot the radial density profile
    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(111)

    plot_grid = [0,99,199]
    alpha = np.linspace(0.3,1.0,len(plot_grid))
    for i in plot_grid:
        rho_rad, = ax.plot(np.log10(rc/AU), np.log10(temp2d[:,i]),'-',color='b',linewidth=2, markersize=3,alpha=alpha[plot_grid.index(i)])

    ax.set_xlabel(r'$\mathrm{log~R~(AU)}$',fontsize=20)
    ax.set_ylabel(r'$\mathrm{log~T~(K)}$',fontsize=20)
    [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
    ax.minorticks_on()
    ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
    ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)
    ax.set_ylim([0,4])
    fig.gca().set_xlim(left=np.log10(0.05))
    # ax.set_xlim([np.log10(0.8),np.log10(10000)])

    fig.savefig(outdir+print_name+'_temp_radial.pdf',format='pdf',dpi=300,bbox_inches='tight')
    fig.clf()

# rtout = '/Users/yaolun/bhr71/hyperion/cycle3/model10/model10.rtout'
# outdir = '/Users/yaolun/test/'
# temp_hyperion(rtout, outdir)