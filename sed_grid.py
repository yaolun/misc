def sed_grid_cs_age(indir, array, outdir, cslist, agelist, obs=None):
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

    fig, axarr = plt.subplots(row, col, sharex='col', sharey='row', figsize=(12,7.1))

    for rr in range(0, row):
        for cc in range(0, col):
            ax = axarr[rr,cc]
            # sed part
            # if rr+1 != row:
            # infinite aperture
            (wave_inf, sed_inf) = np.genfromtxt(indir+'/model'+str(array[rr,cc])+'_sed_inf.txt', skip_header=1).T
            # sed with apertures
            (wave, sed) = np.genfromtxt(indir+'/model'+str(array[rr,cc])+'_sed_w_aperture.txt', skip_header=1).T

            ax.plot(np.log10(wave_inf), np.log10(sed_inf), color='k', linewidth=0.7)
            if obs != None:  
                import sys
                sys.path.append('/Users/yaolun/programs/misc/hyperion')
                from get_bhr71_obs import get_bhr71_obs
                c = const.c.cgs.value

                bhr71 = get_bhr71_obs(obs)  # in um and Jy
                wave_obs, flux_obs, noise_obs = bhr71['spec']
                ax.plot(np.log10(wave_obs[wave_obs<50]), np.log10(c/(wave_obs[wave_obs<50]*1e-4)*flux_obs[wave_obs<50]*1e-23), color='r', alpha=0.7, linewidth=1)
                ax.plot(np.log10(wave_obs[(wave_obs>50)&(wave_obs<190.31)]), np.log10(c/(wave_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-4)*flux_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-23), color='r', alpha=0.7, linewidth=1)
                ax.plot(np.log10(wave_obs[wave_obs>194]), np.log10(c/(wave_obs[wave_obs>194]*1e-4)*flux_obs[wave_obs>194]*1e-23), color='r', alpha=0.7, linewidth=1)

            ax.plot(np.log10(wave), np.log10(sed), 'o-',mfc='b',mec='b',markersize=4,markeredgewidth=1,linewidth=1.2)

            ax.set_ylim([-15,-8])

            [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
            ax.minorticks_on() 
            ax.tick_params('both',labelsize=14,width=1,which='major',pad=15,length=5)
            ax.tick_params('both',labelsize=14,width=1,which='minor',pad=15,length=2.5)

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

    fig.text(0.5, -0.05 , r'$\mathrm{age~[yr]~(1\times 10^{4},~2.5\times 10^{4},~5\times 10^{4},~7.5\times 10^{4},~1\times 10^{5})}$', fontsize=20, ha='center')
    fig.text(0, 0.5, r'$\mathrm{sound~speed~[km~s^{-1}]~(0.3,~0.2,~0.1)}$', fontsize=20, va='center', rotation='vertical')

    fig.subplots_adjust(hspace=0,wspace=0)
    fig.savefig(outdir+'sed_cs_age.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()

def sed_omega(indir, array, outdir, obs=None):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from hyperion.model import ModelOutput
    import astropy.constants as const

    # constants setup
    AU = const.au.cgs.value

    col, = np.shape(array)
    row = 1

    fig, axarr = plt.subplots(row, col, sharex='col', sharey='row', figsize=(12,2.7))

    for cc in range(0, col):
        ax = axarr[cc]
        # sed part
        # if rr+1 != row:
        # infinite aperture
        (wave_inf, sed_inf) = np.genfromtxt(indir+'/model'+str(array[cc])+'_sed_inf.txt', skip_header=1).T
        # sed with apertures
        (wave, sed) = np.genfromtxt(indir+'/model'+str(array[cc])+'_sed_w_aperture.txt', skip_header=1).T

        ax.plot(np.log10(wave_inf), np.log10(sed_inf), color='k', linewidth=1, alpha=0.7)
        if obs != None:  
            import sys
            sys.path.append('/Users/yaolun/programs/misc/hyperion')
            from get_bhr71_obs import get_bhr71_obs
            c = const.c.cgs.value

            bhr71 = get_bhr71_obs(obs)  # in um and Jy
            wave_obs, flux_obs, noise_obs = bhr71['spec']
            ax.plot(np.log10(wave_obs[wave_obs<50]), np.log10(c/(wave_obs[wave_obs<50]*1e-4)*flux_obs[wave_obs<50]*1e-23), color='r', alpha=0.7, linewidth=1)
            ax.plot(np.log10(wave_obs[(wave_obs>50)&(wave_obs<190.31)]), np.log10(c/(wave_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-4)*flux_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-23), color='r', alpha=0.7, linewidth=1)
            ax.plot(np.log10(wave_obs[wave_obs>194]), np.log10(c/(wave_obs[wave_obs>194]*1e-4)*flux_obs[wave_obs>194]*1e-23), color='r', alpha=0.7, linewidth=1)
 
        ax.plot(np.log10(wave), np.log10(sed), 'o-',mfc='b',mec='b',markersize=4,markeredgewidth=1,linewidth=1.2)

        ax.set_ylim([-14,-8])

        [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
        ax.minorticks_on() 
        ax.tick_params('both',labelsize=12,width=1.2,which='major',pad=15,length=5)
        ax.tick_params('both',labelsize=12,width=1.2,which='minor',pad=15,length=2.5)

        if cc == 0:
            ax.set_xlabel(r'$\mathrm{log(wavelength)~(\mu m)}$', fontsize=14)
            ax.set_ylabel(r'$\mathrm{log~\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$', fontsize=14)

        # fix the overlap tick labels
        x_nbins = len(ax.get_xticklabels())
        y_nbins = len(ax.get_yticklabels())
        if (cc != 0):
            ax.xaxis.set_major_locator(MaxNLocator(nbins=x_nbins, prune='lower'))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune='upper'))

    fig.text(0.5, -0.13 , r'$\mathrm{\Omega_{\circ}~[s^{-1}]~(1\times 10^{-13},~5\times 10^{-14},~1\times 10^{-14})}$', fontsize=14, ha='center')

    fig.subplots_adjust(hspace=0,wspace=0)
    fig.savefig(outdir+'sed_omega0.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()

def sed_disk(indir, array, outdir, xlabel, plotname, obs=None):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from hyperion.model import ModelOutput
    import astropy.constants as const
    # constants setup
    AU = const.au.cgs.value

    col, = np.shape(array)
    row = 1

    fig, axarr = plt.subplots(row, col, sharex='col', sharey='row', figsize=(12,2.7))

    for cc in range(0, col):
        ax = axarr[cc]
        # sed part
        # if rr+1 != row:
        # infinite aperture
        (wave_inf, sed_inf) = np.genfromtxt(indir+'/model'+str(array[cc])+'_sed_inf.txt', skip_header=1).T
        # sed with apertures
        (wave, sed) = np.genfromtxt(indir+'/model'+str(array[cc])+'_sed_w_aperture.txt', skip_header=1).T

        ax.plot(np.log10(wave_inf), np.log10(sed_inf), color='k', linewidth=1)
        if obs != None:  
            import sys
            sys.path.append('/Users/yaolun/programs/misc/hyperion')
            from get_bhr71_obs import get_bhr71_obs
            c = const.c.cgs.value

            bhr71 = get_bhr71_obs(obs)  # in um and Jy
            wave_obs, flux_obs, noise_obs = bhr71['spec']
            ax.plot(np.log10(wave_obs[wave_obs<50]), np.log10(c/(wave_obs[wave_obs<50]*1e-4)*flux_obs[wave_obs<50]*1e-23), color='r', alpha=0.7, linewidth=1)
            ax.plot(np.log10(wave_obs[(wave_obs>50)&(wave_obs<190.31)]), np.log10(c/(wave_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-4)*flux_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-23), color='r', alpha=0.7, linewidth=1)
            ax.plot(np.log10(wave_obs[wave_obs>194]), np.log10(c/(wave_obs[wave_obs>194]*1e-4)*flux_obs[wave_obs>194]*1e-23), color='r', alpha=0.7, linewidth=1)

        ax.plot(np.log10(wave), np.log10(sed), 'o-',mfc='b',mec='b',markersize=4,markeredgewidth=1,linewidth=1.2)

        ax.set_ylim([-14,-8])

        [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
        ax.minorticks_on() 
        ax.tick_params('both',labelsize=12,width=1.2,which='major',pad=15,length=5)
        ax.tick_params('both',labelsize=12,width=1.2,which='minor',pad=15,length=2.5)

        if cc == 0:
            ax.set_xlabel(r'$\mathrm{log(wavelength)~(\mu m)}$', fontsize=14)
            ax.set_ylabel(r'$\mathrm{log~\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$', fontsize=14)

        # fix the overlap tick labels
        x_nbins = len(ax.get_xticklabels())
        y_nbins = len(ax.get_yticklabels())
        if (cc != 0):
            ax.xaxis.set_major_locator(MaxNLocator(nbins=x_nbins, prune='lower'))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune='upper'))

    fig.text(0.5, -0.13 , xlabel, fontsize=14, ha='center')

    fig.subplots_adjust(hspace=0,wspace=0)
    fig.savefig(outdir+'sed_disk_'+plotname+'.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()

def sed_disk_exist_com(indir, array, outdir,obs=None):
    import numpy as np
    import matplotlib.pyplot as plt
    from hyperion.model import ModelOutput
    import astropy.constants as const
    # constant setup
    AU = const.cgs.value

    # disk part - model

    fig = plt.figure(figsize=(8,6))

def sed_grid_theta_cav_incl(indir, array, outdir, obs=None):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from hyperion.model import ModelOutput
    import astropy.constants as const
    # constants setup
    AU = const.au.cgs.value

    row, col = np.shape(array)

    fig, axarr = plt.subplots(row, col, sharex='col', sharey='row', figsize=(12,7.1))

    for rr in range(0, row):
        for cc in range(0, col):
            ax = axarr[rr,cc]
            # sed part
            # if rr+1 != row:
            # infinite aperture
            (wave_inf, sed_inf) = np.genfromtxt(indir+'/model'+str(array[rr,cc])+'_sed_inf.txt', skip_header=1).T
            # sed with apertures
            (wave, sed) = np.genfromtxt(indir+'/model'+str(array[rr,cc])+'_sed_w_aperture.txt', skip_header=1).T

            ax.plot(np.log10(wave_inf), np.log10(sed_inf), color='k', linewidth=0.7)
            if obs != None:  
                import sys
                sys.path.append('/Users/yaolun/programs/misc/hyperion')
                from get_bhr71_obs import get_bhr71_obs
                c = const.c.cgs.value

                bhr71 = get_bhr71_obs(obs)  # in um and Jy
                wave_obs, flux_obs, noise_obs = bhr71['spec']
                ax.plot(np.log10(wave_obs[wave_obs<50]), np.log10(c/(wave_obs[wave_obs<50]*1e-4)*flux_obs[wave_obs<50]*1e-23), color='r', alpha=0.7, linewidth=1)
                ax.plot(np.log10(wave_obs[(wave_obs>50)&(wave_obs<190.31)]), np.log10(c/(wave_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-4)*flux_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-23), color='r', alpha=0.7, linewidth=1)
                ax.plot(np.log10(wave_obs[wave_obs>194]), np.log10(c/(wave_obs[wave_obs>194]*1e-4)*flux_obs[wave_obs>194]*1e-23), color='r', alpha=0.7, linewidth=1)

            ax.plot(np.log10(wave), np.log10(sed), 'o-',mfc='b',mec='b',markersize=4,markeredgewidth=1,linewidth=1.2)

            ax.set_ylim([-15,-8])

            [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
            ax.minorticks_on() 
            ax.tick_params('both',labelsize=14,width=1,which='major',pad=15,length=5)
            ax.tick_params('both',labelsize=14,width=1,which='minor',pad=15,length=2.5)

            if rr+1 == row:
                if cc == 0:
                    ax.set_xlabel(r'$\mathrm{log(wavelength)~(\mu m)}$', fontsize=14)
                    ax.set_ylabel(r'$\mathrm{log~\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$', fontsize=14)
        
            # fix the overlap tick labels
            x_nbins = len(ax.get_xticklabels())
            y_nbins = len(ax.get_yticklabels())
            if (rr != 0) & (cc != 0):
                ax.xaxis.set_major_locator(MaxNLocator(nbins=x_nbins, prune='lower'))
                ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune='upper'))

    fig.text(0.5, -0.05 , r'$\mathrm{\theta_{\rm cav}~[deg.]~(15^{\circ},~20^{\circ},~25^{\circ},~30^{\circ},~35^{\circ})}$', fontsize=20, ha='center')
    fig.text(0, 0.5, r'$\mathrm{\theta_{\rm incl.}~[deg.]~(60^{\circ},~70^{\circ},~8 0^{\circ})}$', fontsize=20, va='center', rotation='vertical')

    fig.subplots_adjust(hspace=0,wspace=0)
    fig.savefig(outdir+'sed_theta_cav_incl.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()
def sed_grid_rho_cav_centeredge(indir, array, outdir, obs=None):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from hyperion.model import ModelOutput
    import astropy.constants as const
    # constants setup
    AU = const.au.cgs.value

    row, col = np.shape(array)

    fig, axarr = plt.subplots(row, col, sharex='col', sharey='row', figsize=(9.6,7.1))

    for rr in range(0, row):
        for cc in range(0, col):
            ax = axarr[rr,cc]
            # sed part
            # if rr+1 != row:
            # infinite aperture
            (wave_inf, sed_inf) = np.genfromtxt(indir+'/model'+str(array[rr,cc])+'_sed_inf.txt', skip_header=1).T
            # sed with apertures
            (wave, sed) = np.genfromtxt(indir+'/model'+str(array[rr,cc])+'_sed_w_aperture.txt', skip_header=1).T

            ax.plot(np.log10(wave_inf), np.log10(sed_inf), color='k', linewidth=0.7)
            if obs != None:  
                import sys
                sys.path.append('/Users/yaolun/programs/misc/hyperion')
                from get_bhr71_obs import get_bhr71_obs
                c = const.c.cgs.value

                bhr71 = get_bhr71_obs(obs)  # in um and Jy
                wave_obs, flux_obs, noise_obs = bhr71['spec']
                ax.plot(np.log10(wave_obs[wave_obs<50]), np.log10(c/(wave_obs[wave_obs<50]*1e-4)*flux_obs[wave_obs<50]*1e-23), color='r', alpha=0.7, linewidth=1)
                ax.plot(np.log10(wave_obs[(wave_obs>50)&(wave_obs<190.31)]), np.log10(c/(wave_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-4)*flux_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-23), color='r', alpha=0.7, linewidth=1)
                ax.plot(np.log10(wave_obs[wave_obs>194]), np.log10(c/(wave_obs[wave_obs>194]*1e-4)*flux_obs[wave_obs>194]*1e-23), color='r', alpha=0.7, linewidth=1)

            ax.plot(np.log10(wave), np.log10(sed), 'o-',mfc='b',mec='b',markersize=4,markeredgewidth=1,linewidth=1.2)

            ax.set_ylim([-15,-8])

            [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
            ax.minorticks_on() 
            ax.tick_params('both',labelsize=14,width=1,which='major',pad=15,length=5)
            ax.tick_params('both',labelsize=14,width=1,which='minor',pad=15,length=2.5)

            if rr+1 == row:
                if cc == 0:
                    ax.set_xlabel(r'$\mathrm{log(wavelength)~(\mu m)}$', fontsize=14)
                    ax.set_ylabel(r'$\mathrm{log~\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$', fontsize=14)
        
            # fix the overlap tick labels
            x_nbins = len(ax.get_xticklabels())
            y_nbins = len(ax.get_yticklabels())
            if (rr != 0) & (cc != 0):
                ax.xaxis.set_major_locator(MaxNLocator(nbins=x_nbins, prune='lower'))
                ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune='upper'))

    fig.text(0.5, -0.05 , r'$\mathrm{\rho_{cav,\circ}~[g~cm^{-3}]~(1\times 10^{-18},~5\times 10^{-18},~1\times 10^{-17},~5\times 10^{-17})}$', fontsize=20, ha='center')
    fig.text(0, 0.5, r'$\mathrm{R_{cav,\circ}~[AU]~(40,~30,~20)}$', fontsize=20, va='center', rotation='vertical')

    fig.subplots_adjust(hspace=0,wspace=0)
    fig.savefig(outdir+'sed_rho_cav_centeredge.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()


import numpy as np
indir = '/Users/yaolun/bhr71/hyperion/controlled/'
outdir = '/Users/yaolun/Copy/Papers/yaolun/bhr71/figures/'
obs='/Users/yaolun/bhr71/obs_for_radmc/'

# grid of cs and age
array = np.array([[1,6,11],[2,7,12],[3,8,13],[4,9,14],[5,10,15]])
array = np.array([[1,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15]])
cslist = [0.1,0.2,0.3]
agelist = [1e4,2.5e4,5e4,7.5e4,1e5]
sed_grid_cs_age(indir, array, outdir, cslist, agelist, obs=obs)

# grid of Omega0
array = np.array([16,17,18])
sed_omega(indir, array, outdir, obs=obs)

# grid of disk parameters
# disk mass
array = np.array([19,20,21,22,23])
xlabel = r'$\mathrm{M_{disk}~[M_{\odot}]~(0.1,~0.3,~0.5,~0.7,~1.0)}$'
plotname = 'mdisk'
sed_disk(indir, array, outdir, xlabel, plotname, obs=obs)
# flare power
array = np.array([29,30,31,32,33])
xlabel = r'$\mathrm{\beta~(1.0,~1.2,~1.4,~1.6,~1.8)}$'
plotname = 'beta'
sed_disk(indir, array, outdir, xlabel, plotname, obs=obs)

# grid of theta_cav and incl.
array = np.array([[35,36,37,38,39],[40,41,42,43,44],[45,46,47,48,49]])
sed_grid_theta_cav_incl(indir, array, outdir, obs=obs)

# grid of rho_cav_center and sed_rho_cav_edge
array = np.array([[49,50,51,52],[53,54,55,56],[57,58,59,60]])
sed_grid_rho_cav_centeredge(indir, array, outdir, obs=obs)