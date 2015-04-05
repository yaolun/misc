def sed_grid_cs_age(indir, array, outdir, cslist, agelist):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from hyperion.model import ModelOutput
    import astropy.constants as const
    # constants setup
    AU = const.au.cgs.value

    row, col = np.shape(array)

    # col = col+1
    # row = row+1

    fig, axarr = plt.subplots(row, col, sharex='col', sharey='row', figsize=(16,12))

    for rr in range(0, row):
        for cc in range(0, col):
            ax = axarr[rr,cc]
            # sed part
            # if rr+1 != row:
            # infinite aperture
            (wave_inf, sed_inf) = np.genfromtxt(indir+'/model'+str(array[rr,cc])+'_sed_inf.txt', skip_header=1).T
            # sed with apertures
            (wave, sed) = np.genfromtxt(indir+'/model'+str(array[rr,cc])+'_sed_w_aperture.txt', skip_header=1).T

            ax.plot(np.log10(wave_inf), np.log10(sed_inf), color='k', linewidth=1)
            ax.plot(np.log10(wave), np.log10(sed), 'o-',mfc='None',mec='b',markersize=6,markeredgewidth=2)

            ax.set_ylim([-15,-8])

            [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
            ax.minorticks_on() 
            ax.tick_params('both',labelsize=14,width=1.5,which='major',pad=15,length=5)
            ax.tick_params('both',labelsize=14,width=1.5,which='minor',pad=15,length=2.5)

            if rr+1 == row:
                if cc == 0:
                    ax.set_xlabel(r'$\mathrm{log(wavelength)~(\mu m)}$', fontsize=14)
                    ax.set_ylabel(r'$\mathrm{log~\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$', fontsize=14)
            # # density part
            # else:
            #     m = ModelOutput(indir+'/model'+str(array[1, cc])+'.rtout')
            #     q = m.get_quantities()

            #     # get the grid info
            #     ri, thetai = q.r_wall, q.t_wall
            #     rc     = 0.5*( ri[0:len(ri)-1]     + ri[1:len(ri)] )
            #     thetac = 0.5*( thetai[0:len(thetai)-1] + thetai[1:len(thetai)] )

            #     # get the density profile (dust)
            #     rho = q['density'][0].array.T
            #     rho2d = np.sum(rho**2, axis=2)/np.sum(rho, axis=2)
            #     # rho2d_exp = np.hstack((rho2d,rho2d,rho2d[:,0:1]))
            #     # thetac_exp = np.hstack((thetac-np.pi/2, thetac+np.pi/2, thetac[0]-np.pi/2))
            #     # midplane
            #     ax.plot(np.log10(rc/AU), np.log10(rho2d[:,0]),'-',color='k',linewidth=1)
            #     # pole
            #     ax.plot(np.log10(rc/AU), np.log10(rho2d[:,199]),'--',color='k',linewidth=1)

            #     [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
            #     ax.minorticks_on() 
            #     ax.tick_params('both',labelsize=14,width=1.5,which='major',pad=15,length=5)
            #     ax.tick_params('both',labelsize=14,width=1.5,which='minor',pad=15,length=2.5)

            #     ax.yaxis.tick_right()
            #     ax.yaxis.set_label_position("right")
            #     if rr+1 == row:
            #         ax.set_xlabel(r'$\mathrm{radius~[AU]}$', fontsize=14)
            #         ax.set_ylabel(r'$\mathrm{dust~density~[g~cm^{-3}]}$', fontsize=14)

            # fix the overlap tick labels
            x_nbins = len(ax.get_xticklabels())
            y_nbins = len(ax.get_yticklabels())
            if (rr != 0) & (cc != 0):
                ax.xaxis.set_major_locator(MaxNLocator(nbins=x_nbins, prune='lower'))
                ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune='upper'))

    fig.text(0.5, 0.03 , r'$\mathrm{age~[yr]~(1\times 10^{4},~2.5\times 10^{4},~5\times 10^{4},~7.5\times 10^{4},~1\times 10^{5})}$', fontsize=20, ha='center')
    fig.text(0.05, 0.5, r'$\mathrm{sound~speed~[km~s^{-1}]~(0.1,~0.2,~0.3)}$', fontsize=20, va='center', rotation='vertical')

    fig.subplots_adjust(hspace=0,wspace=0)
    fig.savefig(outdir+'sed_cs_age.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()

import numpy as np
indir = '/Users/yaolun/bhr71/hyperion/controlled/'
outdir = '/Users/yaolun/test/'
array = np.array([[1,6,11],[2,7,12],[3,8,13],[4,9,14],[5,10,15]])
array = np.array([[1,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15]])
cslist = [0.1,0.2,0.3]
agelist = [1e4,2.5e4,5e4,7.5e4,1e5]
sed_grid_cs_age(indir, array, outdir, cslist, agelist)