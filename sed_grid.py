def sed_grid_cs_age(indir, array, outdir, cslist, agelist, obs=None):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from hyperion.model import ModelOutput
    import astropy.constants as const
    import sys
    sys.path.append('/Users/yaolun/programs/misc/hyperion')
    from t_bol import t_bol

    # constants setup
    AU = const.au.cgs.value
    c = const.c.cgs.value

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

            # ax.plot(np.log10(wave_inf), np.log10(sed_inf), color='k', linewidth=0.7)
            if obs != None:  
                from get_bhr71_obs import get_bhr71_obs

                bhr71 = get_bhr71_obs(obs)  # in um and Jy
                wave_obs, flux_obs, noise_obs = bhr71['spec']
                ax.plot(np.log10(wave_obs[wave_obs<50]), np.log10(c/(wave_obs[wave_obs<50]*1e-4)*flux_obs[wave_obs<50]*1e-23), color='r', alpha=0.7, linewidth=1)
                ax.plot(np.log10(wave_obs[(wave_obs>50)&(wave_obs<190.31)]), np.log10(c/(wave_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-4)*flux_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-23), color='r', alpha=0.7, linewidth=1)
                ax.plot(np.log10(wave_obs[wave_obs>194]), np.log10(c/(wave_obs[wave_obs>194]*1e-4)*flux_obs[wave_obs>194]*1e-23), color='r', alpha=0.7, linewidth=1)

            ax.plot(np.log10(wave), np.log10(sed), 'o-',mfc='b',mec='b',markersize=4,markeredgewidth=1,linewidth=1.2)
            # print the bolometric temperature
            ax.text(0.4, 0.1, r'$\mathrm{T_{bol}= %4.1f~K}$' % t_bol(wave,sed*wave*1e-4/c), fontsize=12, transform=ax.transAxes)

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

def sed_omega(indir, array, outdir, obs=None, compact=False):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from hyperion.model import ModelOutput
    import astropy.constants as const

    # constants setup
    AU = const.au.cgs.value

    if compact == False:
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

            # ax.plot(np.log10(wave_inf), np.log10(sed_inf), color='k', linewidth=1, alpha=0.7)
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
            ax.set_xlim([0.4,2])

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
    else:
        color_list = ['#ffeda0','#feb24c','#f03b20']
        label = [r'$\mathrm{\Omega_{\circ}=1\times 10^{-14}~s^{-1}}$',r'$\mathrm{\Omega_{\circ}=5\times 10^{-14}~s^{-1}}$',r'$\mathrm{\Omega_{\circ}=1\times 10^{-13}~s^{-1}}$']
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)

        for i in range(0, len(array)):
            # sed with infinite aperture
            (wave_inf, sed_inf) = np.genfromtxt(indir+'/model'+str(array[i])+'_sed_inf.txt', skip_header=1).T
            # sed with apertures
            (wave, sed) = np.genfromtxt(indir+'/model'+str(array[i])+'_sed_w_aperture.txt', skip_header=1).T

            ax.plot(np.log10(wave), np.log10(sed), 'o-',mfc=color_list[i],mec=color_list[i],color=color_list[i],markersize=4,markeredgewidth=1,linewidth=1.2,label=label[i])
        ax.legend(loc='lower right', numpoints=1, framealpha=0.3, fontsize=16)

        [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
        ax.minorticks_on() 
        ax.tick_params('both',labelsize=14,width=1.5,which='major',pad=15,length=5)
        ax.tick_params('both',labelsize=14,width=1.5,which='minor',pad=15,length=2.5)

        ax.set_xlabel(r'$\mathrm{log(wavelength)~(\mu m)}$', fontsize=18)
        ax.set_ylabel(r'$\mathrm{log~\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$', fontsize=18)
        ax.set_ylim([-14,-8])
        # ax.set_xlim([0.4,2])
        # ax.locator_params(axis='x', nbins=5)

        fig.savefig(outdir+'sed_omega0.pdf', format='pdf', dpi=300, bbox_inches='tight')
        fig.clf()

def sed_five(indir, array, outdir, xlabel, plotname, obs=None, zoom=False, tbol=False, compact=None, yrange=None):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from hyperion.model import ModelOutput
    import astropy.constants as const
    import sys
    sys.path.append('/Users/yaolun/programs/misc/hyperion')
    from t_bol import t_bol

    # constants setup
    AU = const.au.cgs.value
    c = const.c.cgs.value

    if compact == None:
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

            # ax.plot(np.log10(wave_inf), np.log10(sed_inf), color='k', linewidth=1)
            if obs != None:  
                
                from get_bhr71_obs import get_bhr71_obs

                bhr71 = get_bhr71_obs(obs)  # in um and Jy
                wave_obs, flux_obs, noise_obs = bhr71['spec']
                ax.plot(np.log10(wave_obs[wave_obs<50]), np.log10(c/(wave_obs[wave_obs<50]*1e-4)*flux_obs[wave_obs<50]*1e-23), color='r', alpha=0.7, linewidth=1)
                ax.plot(np.log10(wave_obs[(wave_obs>50)&(wave_obs<190.31)]), np.log10(c/(wave_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-4)*flux_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-23), color='r', alpha=0.7, linewidth=1)
                ax.plot(np.log10(wave_obs[wave_obs>194]), np.log10(c/(wave_obs[wave_obs>194]*1e-4)*flux_obs[wave_obs>194]*1e-23), color='r', alpha=0.7, linewidth=1)

            ax.plot(np.log10(wave), np.log10(sed), 'o-',mfc='b',mec='b',markersize=4,markeredgewidth=1,linewidth=1.2)

            if tbol == True:
                 ax.text(0.4, 0.1, r'$\mathrm{T_{bol}= %4.1f~K}$' % t_bol(wave, sed*wave*1e-4/c), fontsize=12, transform=ax.transAxes)

            if zoom == True:
                ax.set_xlim([0.4,2])

            [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
            ax.minorticks_on() 
            ax.tick_params('both',labelsize=10,width=1.2,which='major',pad=15,length=5)
            ax.tick_params('both',labelsize=10,width=1.2,which='minor',pad=15,length=2.5)

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
        fig.savefig(outdir+'sed_'+plotname+'.pdf', format='pdf', dpi=300, bbox_inches='tight')
        fig.clf()
    else:
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)

        # color map setting
        # color_list = plt.cm.Set1(np.linspace(0, 1, len(array)+1))
        color_list = ['#fed976','#feb24c','#fd8d3c','#f03b20','#bd0026' ]
        # color_list = ['#feebe2','#fbb4b9','#f768a1','#c51b8a','#7a0177']

        if obs != None:  
            from get_bhr71_obs import get_bhr71_obs
            bhr71 = get_bhr71_obs(obs)  # in um and Jy
            wave_obs, flux_obs, noise_obs = bhr71['spec']
            ax.plot(np.log10(wave_obs[wave_obs<50]), np.log10(c/(wave_obs[wave_obs<50]*1e-4)*flux_obs[wave_obs<50]*1e-23), color='r', alpha=0.7, linewidth=1)
            ax.plot(np.log10(wave_obs[(wave_obs>50)&(wave_obs<190.31)]), np.log10(c/(wave_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-4)*flux_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-23), color='r', alpha=0.7, linewidth=1)
            ax.plot(np.log10(wave_obs[wave_obs>194]), np.log10(c/(wave_obs[wave_obs>194]*1e-4)*flux_obs[wave_obs>194]*1e-23), color='r', alpha=0.7, linewidth=1)

        for i in range(0, len(array)):
            # infinite aperture
            (wave_inf, sed_inf) = np.genfromtxt(indir+'/model'+str(array[i])+'_sed_inf.txt', skip_header=1).T
            # sed with apertures
            (wave, sed) = np.genfromtxt(indir+'/model'+str(array[i])+'_sed_w_aperture.txt', skip_header=1).T

            ax.plot(np.log10(wave), np.log10(sed), 'o-',mfc=color_list[i],mec=color_list[i],color=color_list[i],markersize=4,markeredgewidth=1,linewidth=1.2,label=compact[i])
            ax.legend(loc='lower right', numpoints=1, framealpha=0.3, fontsize=16)
        ax.set_xlabel(r'$\mathrm{log(wavelength)~(\mu m)}$', fontsize=18)
        ax.set_ylabel(r'$\mathrm{log~\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$', fontsize=18)
        ax.set_ylim([-14,-8])
        if zoom == True:
            ax.set_xlim([0.4,2])
        if yrange == None:
            ax.set_ylim([-14,-8])
        else:
            ax.set_ylim(yrange)

        [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
        ax.minorticks_on() 
        ax.tick_params('both',labelsize=16,width=1.2,which='major',pad=15,length=5)
        ax.tick_params('both',labelsize=16,width=1.2,which='minor',pad=15,length=2.5)

        fig.savefig(outdir+'sed_'+plotname+'.pdf', format='pdf', dpi=300, bbox_inches='tight')


def sed_grid_theta_cav_incl(indir, array, outdir, obs=None, compact=False):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from hyperion.model import ModelOutput
    import astropy.constants as const
    # constants setup
    AU = const.au.cgs.value

    if compact == False:
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

                # ax.plot(np.log10(wave_inf), np.log10(sed_inf), color='k', linewidth=0.7)
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
                ax.set_xlim([0.4,2])
                ax.locator_params(axis='x', nbins=5)

                [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
                ax.minorticks_on() 
                ax.tick_params('both',labelsize=10,width=1,which='major',pad=15,length=5)
                ax.tick_params('both',labelsize=10,width=1,which='minor',pad=15,length=2.5)

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
    else:
        # plot multiple SEDs in same panel
        num, col = np.shape(array)
        row = 1

        # set up the color map
        # color_list = ['#f03b20','#feb24c','#ffeda0']
        color_list = ['#ffeda0','#feb24c','#f03b20']
        # label setup
        label = [ r'$\mathrm{\theta_{incl.}=80^{\circ}}$', r'$\mathrm{\theta_{incl.}=70^{\circ}}$',\
                  r'$\mathrm{\theta_{incl.}=60^{\circ}}$']

        fig, axarr = plt.subplots(row, col, sharex='col', sharey='row', figsize=(12,2.7))

        for cc in range(0, col):
            ax = axarr[cc]
            # obs part
            if obs != None:  
                
                from get_bhr71_obs import get_bhr71_obs

                bhr71 = get_bhr71_obs(obs)  # in um and Jy
                wave_obs, flux_obs, noise_obs = bhr71['spec']
                ax.plot(np.log10(wave_obs[wave_obs<50]), np.log10(c/(wave_obs[wave_obs<50]*1e-4)*flux_obs[wave_obs<50]*1e-23), color='r', alpha=0.7, linewidth=1)
                ax.plot(np.log10(wave_obs[(wave_obs>50)&(wave_obs<190.31)]), np.log10(c/(wave_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-4)*flux_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-23), color='r', alpha=0.7, linewidth=1)
                ax.plot(np.log10(wave_obs[wave_obs>194]), np.log10(c/(wave_obs[wave_obs>194]*1e-4)*flux_obs[wave_obs>194]*1e-23), color='r', alpha=0.7, linewidth=1)

            for i in range(0, num):
                # infinite aperture
                (wave_inf, sed_inf) = np.genfromtxt(indir+'/model'+str(array[i,cc])+'_sed_inf.txt', skip_header=1).T
                # sed with apertures
                (wave, sed) = np.genfromtxt(indir+'/model'+str(array[i,cc])+'_sed_w_aperture.txt', skip_header=1).T

                ax.plot(np.log10(wave), np.log10(sed), 'o-',mfc=color_list[i],mec=color_list[i], color=color_list[i], markersize=4,markeredgewidth=1,linewidth=1.2, label=label[i])
                if cc == 0:
                    ax.legend(loc='upper left', numpoints=1,framealpha=0.3,fontsize=10)

            ax.set_ylim([-14,-8])
            ax.set_xlim([0.4,2])
            ax.locator_params(axis='x', nbins=5)

            [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
            ax.minorticks_on() 
            ax.tick_params('both',labelsize=12,width=1.2,which='major',pad=10,length=5)
            ax.tick_params('both',labelsize=12,width=1.2,which='minor',pad=10,length=2.5)

            if cc == 0:
                ax.set_xlabel(r'$\mathrm{log(wavelength)~(\mu m)}$', fontsize=16)
                ax.set_ylabel(r'$\mathrm{log~\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$', fontsize=16)

            # fix the overlap tick labels
            x_nbins = len(ax.get_xticklabels())
            y_nbins = len(ax.get_yticklabels())
            if (cc != 0):
                ax.xaxis.set_major_locator(MaxNLocator(nbins=x_nbins, prune='lower'))
                ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune='upper'))

        # fig.text(0.5, -0.05 , r'$\mathrm{\rho_{cav,\circ}~[g~cm^{-3}]~(1\times 10^{-20},~5\times 10^{-20},~1\times 10^{-19},~5\times 10^{-19})}$', fontsize=20, ha='center')
        fig.text(0.5, -0.15, r'$\mathrm{\theta_{cav}~[deg.]~(15^{\circ},~20^{\circ},~25^{\circ},~30^{\circ},~35^{\circ})}$', fontsize=20, ha='center' )

        fig.subplots_adjust(hspace=0,wspace=0)
        fig.savefig(outdir+'sed_theta_cav_incl.pdf', format='pdf', dpi=300, bbox_inches='tight')
        fig.clf()
def sed_grid_rho_cav_centeredge(indir, array, outdir, obs=None, compact=False):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from hyperion.model import ModelOutput
    import astropy.constants as const
    # constants setup
    AU = const.au.cgs.value

    if compact == False:
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

                # ax.plot(np.log10(wave_inf), np.log10(sed_inf), color='k', linewidth=0.7)
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
                ax.set_xlim([0.4,2])

                [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
                ax.minorticks_on() 
                ax.tick_params('both',labelsize=10,width=1,which='major',pad=15,length=5)
                ax.tick_params('both',labelsize=10,width=1,which='minor',pad=15,length=2.5)

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

        fig.text(0.5, -0.05 , r'$\mathrm{\rho_{cav,\circ}~[g~cm^{-3}]~(1\times 10^{-20},~5\times 10^{-20},~1\times 10^{-19},~5\times 10^{-19})}$', fontsize=20, ha='center')
        fig.text(0, 0.5, r'$\mathrm{R_{cav,\circ}~[AU]~(40,~30,~20)}$', fontsize=20, va='center', rotation='vertical')

        fig.subplots_adjust(hspace=0,wspace=0)
        fig.savefig(outdir+'sed_rho_cav_centeredge.pdf', format='pdf', dpi=300, bbox_inches='tight')
        fig.clf()
    else:
        # plot multiple SEDs in same panel
        col, num = np.shape(array)
        row = 1

        # set up the color map
        color_list = plt.cm.gnuplot2(np.linspace(0, 1, num+1))
        # label setup
        label = [ r'$\mathrm{1\times 10^{-20}}$', r'$\mathrm{5\times 10^{-20}}$',\
                  r'$\mathrm{1\times 10^{-19}}$', r'$\mathrm{5\times 10^{-19}}$']

        fig, axarr = plt.subplots(row, col, sharex='col', sharey='row', figsize=(12,2.7))

        for cc in range(0, col):
            ax = axarr[cc]
            # obs part
            if obs != None:  
                
                from get_bhr71_obs import get_bhr71_obs

                bhr71 = get_bhr71_obs(obs)  # in um and Jy
                wave_obs, flux_obs, noise_obs = bhr71['spec']
                ax.plot(np.log10(wave_obs[wave_obs<50]), np.log10(c/(wave_obs[wave_obs<50]*1e-4)*flux_obs[wave_obs<50]*1e-23), color='r', alpha=0.7, linewidth=1)
                ax.plot(np.log10(wave_obs[(wave_obs>50)&(wave_obs<190.31)]), np.log10(c/(wave_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-4)*flux_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-23), color='r', alpha=0.7, linewidth=1)
                ax.plot(np.log10(wave_obs[wave_obs>194]), np.log10(c/(wave_obs[wave_obs>194]*1e-4)*flux_obs[wave_obs>194]*1e-23), color='r', alpha=0.7, linewidth=1)

            for i in range(0, num):
                # infinite aperture
                (wave_inf, sed_inf) = np.genfromtxt(indir+'/model'+str(array[cc,i])+'_sed_inf.txt', skip_header=1).T
                # sed with apertures
                (wave, sed) = np.genfromtxt(indir+'/model'+str(array[cc,i])+'_sed_w_aperture.txt', skip_header=1).T

                ax.plot(np.log10(wave), np.log10(sed), 'o-',mfc=color_list[i],mec=color_list[i], color=color_list[i], markersize=4,markeredgewidth=1,linewidth=1.2, label=label[i])
                if cc == col-1:
                    ax.legend(loc='lower right', numpoints=1,framealpha=0.3,fontsize=12)

            ax.set_ylim([-14,-8])
            ax.set_xlim([0.4,2])
            ax.locator_params(axis='x', nbins=5)

            [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
            ax.minorticks_on() 
            ax.tick_params('both',labelsize=12,width=1.2,which='major',pad=10,length=5)
            ax.tick_params('both',labelsize=12,width=1.2,which='minor',pad=10,length=2.5)

            if cc == 0:
                ax.set_xlabel(r'$\mathrm{log(wavelength)~(\mu m)}$', fontsize=16)
                ax.set_ylabel(r'$\mathrm{log~\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$', fontsize=16)

            # fix the overlap tick labels
            x_nbins = len(ax.get_xticklabels())
            y_nbins = len(ax.get_yticklabels())
            if (cc != 0):
                ax.xaxis.set_major_locator(MaxNLocator(nbins=x_nbins, prune='lower'))
                ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune='upper'))

        # fig.text(0.5, -0.05 , r'$\mathrm{\rho_{cav,\circ}~[g~cm^{-3}]~(1\times 10^{-20},~5\times 10^{-20},~1\times 10^{-19},~5\times 10^{-19})}$', fontsize=20, ha='center')
        fig.text(0.5, -0.15, r'$\mathrm{R_{cav,\circ}~[AU]~(20,~30,~40)}$', fontsize=20, ha='center' )

        fig.subplots_adjust(hspace=0,wspace=0)
        fig.savefig(outdir+'sed_rho_cav_centeredge.pdf', format='pdf', dpi=300, bbox_inches='tight')
        fig.clf()



def disk_exist_com(indir, array, outdir, obs=None):
    import numpy as np
    import matplotlib.pyplot as plt
    from hyperion.model import ModelOutput
    import astropy.constants as const
    # constants setup
    AU = const.au.cgs.value

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)

    # get data
    # disk part
    (d_wave_inf, d_sed_inf) = np.genfromtxt(indir+'/model'+str(array[0])+'_sed_inf.txt', skip_header=1).T
    (d_wave, d_sed) = np.genfromtxt(indir+'/model'+str(array[0])+'_sed_w_aperture.txt', skip_header=1).T

    # no disk part
    (nd_wave_inf, nd_sed_inf) = np.genfromtxt(indir+'/model'+str(array[1])+'_sed_inf.txt', skip_header=1).T
    (nd_wave, nd_sed) = np.genfromtxt(indir+'/model'+str(array[1])+'_sed_w_aperture.txt', skip_header=1).T

    disk, = ax.plot(np.log10(nd_wave), np.log10(nd_sed), 'o-',mfc='b',mec='b',markersize=5,markeredgewidth=1,color='b', linewidth=1.5)
    nodisk, = ax.plot(np.log10(d_wave), np.log10(d_sed), 'o-',mfc='k',mec='k',markersize=5,markeredgewidth=1,color='k', linewidth=1.5)

    if obs != None:  
        import sys
        sys.path.append('/Users/yaolun/programs/misc/hyperion')
        from get_bhr71_obs import get_bhr71_obs
        c = const.c.cgs.value

        bhr71 = get_bhr71_obs(obs)  # in um and Jy
        wave_obs, flux_obs, noise_obs = bhr71['spec']
        obs_data, = ax.plot(np.log10(wave_obs[wave_obs<50]), np.log10(c/(wave_obs[wave_obs<50]*1e-4)*flux_obs[wave_obs<50]*1e-23), color='r', alpha=0.7, linewidth=1)
        ax.plot(np.log10(wave_obs[(wave_obs>50)&(wave_obs<190.31)]), np.log10(c/(wave_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-4)*flux_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-23), color='r', alpha=0.7, linewidth=1)
        ax.plot(np.log10(wave_obs[wave_obs>194]), np.log10(c/(wave_obs[wave_obs>194]*1e-4)*flux_obs[wave_obs>194]*1e-23), color='r', alpha=0.7, linewidth=1)


    [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
    ax.minorticks_on() 
    ax.tick_params('both',labelsize=16,width=1.5,which='major',pad=15,length=5)
    ax.tick_params('both',labelsize=16,width=1.5,which='minor',pad=15,length=2.5)

    ax.set_xlabel(r'$\mathrm{log(wavelength)~(\mu m)}$', fontsize=18)
    ax.set_ylabel(r'$\mathrm{log~\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$', fontsize=18)

    plt.legend([disk, nodisk, obs_data], [r'$\mathrm{w/~disk}$', r'$\mathrm{w/o~disk}$',r'$\mathrm{observation}$'], numpoints=1, loc='lower right', fontsize=16)

    fig.savefig(outdir+'sed_disk_com.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()


def sed_tstar(indir, array, outdir, obs=None):
    import numpy as np
    import matplotlib.pyplot as plt
    from hyperion.model import ModelOutput
    import astropy.constants as const
    # constants setup
    AU = const.au.cgs.value

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)

    # get data
    # tstar = 4500 K
    (t1_wave_inf, t1_sed_inf) = np.genfromtxt(indir+'/model'+str(array[0])+'_sed_inf.txt', skip_header=1).T
    (t1_wave, t1_sed) = np.genfromtxt(indir+'/model'+str(array[0])+'_sed_w_aperture.txt', skip_header=1).T

    # tstar = 5000 K
    (t2_wave_inf, t2_sed_inf) = np.genfromtxt(indir+'/model'+str(array[1])+'_sed_inf.txt', skip_header=1).T
    (t2_wave, t2_sed) = np.genfromtxt(indir+'/model'+str(array[1])+'_sed_w_aperture.txt', skip_header=1).T

    # tstar = 5500 K
    (t3_wave_inf, t3_sed_inf) = np.genfromtxt(indir+'/model'+str(array[2])+'_sed_inf.txt', skip_header=1).T
    (t3_wave, t3_sed) = np.genfromtxt(indir+'/model'+str(array[2])+'_sed_w_aperture.txt', skip_header=1).T

    t1, = ax.plot(np.log10(t1_wave), np.log10(t1_sed), 'o-',mfc='Magenta',mec='Magenta',markersize=5,markeredgewidth=1,color='Magenta', linewidth=1.5)
    t2, = ax.plot(np.log10(t2_wave), np.log10(t2_sed), 'o-',mfc='r',mec='r',markersize=5,markeredgewidth=1,color='r', linewidth=1.5)
    t3, = ax.plot(np.log10(t3_wave), np.log10(t3_sed), 'o-',mfc='b',mec='b',markersize=5,markeredgewidth=1,color='b', linewidth=1.5)

    if obs != None:  
        import sys
        sys.path.append('/Users/yaolun/programs/misc/hyperion')
        from get_bhr71_obs import get_bhr71_obs
        c = const.c.cgs.value

        bhr71 = get_bhr71_obs(obs)  # in um and Jy
        wave_obs, flux_obs, noise_obs = bhr71['spec']
        obs_data, = ax.plot(np.log10(wave_obs[wave_obs<50]), np.log10(c/(wave_obs[wave_obs<50]*1e-4)*flux_obs[wave_obs<50]*1e-23), color='k', alpha=0.7, linewidth=1)
        ax.plot(np.log10(wave_obs[(wave_obs>50)&(wave_obs<190.31)]), np.log10(c/(wave_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-4)*flux_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-23), color='k', alpha=0.7, linewidth=1)
        ax.plot(np.log10(wave_obs[wave_obs>194]), np.log10(c/(wave_obs[wave_obs>194]*1e-4)*flux_obs[wave_obs>194]*1e-23), color='k', alpha=0.7, linewidth=1)


    [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
    ax.minorticks_on() 
    ax.tick_params('both',labelsize=14,width=1.5,which='major',pad=15,length=5)
    ax.tick_params('both',labelsize=14,width=1.5,which='minor',pad=15,length=2.5)

    ax.set_xlabel(r'$\mathrm{log(wavelength)~(\mu m)}$', fontsize=16)
    ax.set_ylabel(r'$\mathrm{log~\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$', fontsize=16)
    ax.set_ylim([-13,-7])

    plt.legend([t1, t2, t3, obs_data], [r'$\mathrm{T_{\star}=4500~K}$', r'$\mathrm{T_{\star}=5000~K}$',r'$\mathrm{T_{\star}=5500~K}$',r'$\mathrm{observation}$'], numpoints=1, loc='lower right', fontsize=16)

    fig.savefig(outdir+'sed_tstar.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()

def sed_rstar(indir, array, outdir, obs=None):
    import numpy as np
    import matplotlib.pyplot as plt
    from hyperion.model import ModelOutput
    import astropy.constants as const
    # constants setup
    AU = const.au.cgs.value

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)

    # get data
    # tstar = 4500 K
    (r1_wave_inf, r1_sed_inf) = np.genfromtxt(indir+'/model'+str(array[0])+'_sed_inf.txt', skip_header=1).T
    (r1_wave, r1_sed) = np.genfromtxt(indir+'/model'+str(array[0])+'_sed_w_aperture.txt', skip_header=1).T

    # tstar = 5000 K
    (r2_wave_inf, r2_sed_inf) = np.genfromtxt(indir+'/model'+str(array[1])+'_sed_inf.txt', skip_header=1).T
    (r2_wave, r2_sed) = np.genfromtxt(indir+'/model'+str(array[1])+'_sed_w_aperture.txt', skip_header=1).T

    # tstar = 5500 K
    (r3_wave_inf, r3_sed_inf) = np.genfromtxt(indir+'/model'+str(array[2])+'_sed_inf.txt', skip_header=1).T
    (r3_wave, r3_sed) = np.genfromtxt(indir+'/model'+str(array[2])+'_sed_w_aperture.txt', skip_header=1).T

    r1, = ax.plot(np.log10(r1_wave), np.log10(r1_sed), 'o-',mfc='Magenta',mec='Magenta',markersize=5,markeredgewidth=1,color='Magenta', linewidth=1.5)
    r2, = ax.plot(np.log10(r2_wave), np.log10(r2_sed), 'o-',mfc='r',mec='r',markersize=5,markeredgewidth=1,color='r', linewidth=1.5)
    r3, = ax.plot(np.log10(r3_wave), np.log10(r3_sed), 'o-',mfc='b',mec='b',markersize=5,markeredgewidth=1,color='b', linewidth=1.5)

    if obs != None:  
        import sys
        sys.path.append('/Users/yaolun/programs/misc/hyperion')
        from get_bhr71_obs import get_bhr71_obs
        c = const.c.cgs.value

        bhr71 = get_bhr71_obs(obs)  # in um and Jy
        wave_obs, flux_obs, noise_obs = bhr71['spec']
        obs_data, = ax.plot(np.log10(wave_obs[wave_obs<50]), np.log10(c/(wave_obs[wave_obs<50]*1e-4)*flux_obs[wave_obs<50]*1e-23), color='k', alpha=0.7, linewidth=1)
        ax.plot(np.log10(wave_obs[(wave_obs>50)&(wave_obs<190.31)]), np.log10(c/(wave_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-4)*flux_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-23), color='k', alpha=0.7, linewidth=1)
        ax.plot(np.log10(wave_obs[wave_obs>194]), np.log10(c/(wave_obs[wave_obs>194]*1e-4)*flux_obs[wave_obs>194]*1e-23), color='k', alpha=0.7, linewidth=1)


    [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
    ax.minorticks_on() 
    ax.tick_params('both',labelsize=14,width=1.5,which='major',pad=15,length=5)
    ax.tick_params('both',labelsize=14,width=1.5,which='minor',pad=15,length=2.5)

    ax.set_xlabel(r'$\mathrm{log(wavelength)~(\mu m)}$', fontsize=16)
    ax.set_ylabel(r'$\mathrm{log~\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$', fontsize=16)
    ax.set_ylim([-13,-7])

    plt.legend([r1, r2, r3, obs_data], [r'$\mathrm{R_{\star}=0.3~R_{\odot}}$', r'$\mathrm{R_{\star}=0.5~R_{\odot}}$',r'$\mathrm{R_{\star}=0.7~R_{\odot}}$',r'$\mathrm{observation}$'], numpoints=1, loc='lower right', fontsize=16)

    fig.savefig(outdir+'sed_rstar.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()

def sed_cav_powerlaw(indir, array, outdir, obs=None):
    import numpy as np
    import matplotlib.pyplot as plt
    from hyperion.model import ModelOutput
    import astropy.constants as const
    # constants setup
    AU = const.au.cgs.value

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)

    # get data
    # tstar = 4500 K
    (r1_wave_inf, r1_sed_inf) = np.genfromtxt(indir+'/model'+str(array[0])+'_sed_inf.txt', skip_header=1).T
    (r1_wave, r1_sed) = np.genfromtxt(indir+'/model'+str(array[0])+'_sed_w_aperture.txt', skip_header=1).T

    # tstar = 5000 K
    (r2_wave_inf, r2_sed_inf) = np.genfromtxt(indir+'/model'+str(array[1])+'_sed_inf.txt', skip_header=1).T
    (r2_wave, r2_sed) = np.genfromtxt(indir+'/model'+str(array[1])+'_sed_w_aperture.txt', skip_header=1).T

    # tstar = 5500 K
    (r3_wave_inf, r3_sed_inf) = np.genfromtxt(indir+'/model'+str(array[2])+'_sed_inf.txt', skip_header=1).T
    (r3_wave, r3_sed) = np.genfromtxt(indir+'/model'+str(array[2])+'_sed_w_aperture.txt', skip_header=1).T

    r1, = ax.plot(np.log10(r1_wave), np.log10(r1_sed), 'o-',mfc='Magenta',mec='Magenta',markersize=5,markeredgewidth=1,color='Magenta', linewidth=1.5)
    r2, = ax.plot(np.log10(r2_wave), np.log10(r2_sed), 'o-',mfc='r',mec='r',markersize=5,markeredgewidth=1,color='r', linewidth=1.5)
    r3, = ax.plot(np.log10(r3_wave), np.log10(r3_sed), 'o-',mfc='b',mec='b',markersize=5,markeredgewidth=1,color='b', linewidth=1.5)

    if obs != None:  
        import sys
        sys.path.append('/Users/yaolun/programs/misc/hyperion')
        from get_bhr71_obs import get_bhr71_obs
        c = const.c.cgs.value

        bhr71 = get_bhr71_obs(obs)  # in um and Jy
        wave_obs, flux_obs, noise_obs = bhr71['spec']
        obs_data, = ax.plot(np.log10(wave_obs[wave_obs<50]), np.log10(c/(wave_obs[wave_obs<50]*1e-4)*flux_obs[wave_obs<50]*1e-23), color='k', alpha=0.7, linewidth=1)
        ax.plot(np.log10(wave_obs[(wave_obs>50)&(wave_obs<190.31)]), np.log10(c/(wave_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-4)*flux_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-23), color='k', alpha=0.7, linewidth=1)
        ax.plot(np.log10(wave_obs[wave_obs>194]), np.log10(c/(wave_obs[wave_obs>194]*1e-4)*flux_obs[wave_obs>194]*1e-23), color='k', alpha=0.7, linewidth=1)


    [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
    ax.minorticks_on() 
    ax.tick_params('both',labelsize=14,width=1.5,which='major',pad=15,length=5)
    ax.tick_params('both',labelsize=14,width=1.5,which='minor',pad=15,length=2.5)

    ax.set_xlabel(r'$\mathrm{log(wavelength)~(\mu m)}$', fontsize=16)
    ax.set_ylabel(r'$\mathrm{log~\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$', fontsize=16)
    ax.set_ylim([-13,-7])

    plt.legend([r1, r2, r3, obs_data], [r'$\mathrm{\rho_{cav,\circ}=5\times 10^{-15}~g~cm^{-3}}$', r'$\mathrm{\rho_{cav,\circ}=5\times 10^{-16}~g~cm^{-3}}$',r'$\mathrm{\rho_{cav,\circ}=5\times 10^{-17}~g~cm^{-3}}$',r'$\mathrm{observation}$'], numpoints=1, loc='lower right', fontsize=16)

    fig.savefig(outdir+'sed_cav_cont_powerlaw.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()

def sed_cav_struc_com(indir, array, outdir, obs=None):
    import numpy as np
    import matplotlib.pyplot as plt
    from hyperion.model import ModelOutput
    import astropy.constants as const
    # constants setup
    AU = const.au.cgs.value

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)

    # get data
    # tstar = 4500 K
    (r1_wave_inf, r1_sed_inf) = np.genfromtxt(indir+'/model'+str(array[0])+'_sed_inf.txt', skip_header=1).T
    (r1_wave, r1_sed) = np.genfromtxt(indir+'/model'+str(array[0])+'_sed_w_aperture.txt', skip_header=1).T

    # tstar = 5000 K
    (r2_wave_inf, r2_sed_inf) = np.genfromtxt(indir+'/model'+str(array[1])+'_sed_inf.txt', skip_header=1).T
    (r2_wave, r2_sed) = np.genfromtxt(indir+'/model'+str(array[1])+'_sed_w_aperture.txt', skip_header=1).T

    # tstar = 5500 K
    (r3_wave_inf, r3_sed_inf) = np.genfromtxt(indir+'/model'+str(array[2])+'_sed_inf.txt', skip_header=1).T
    (r3_wave, r3_sed) = np.genfromtxt(indir+'/model'+str(array[2])+'_sed_w_aperture.txt', skip_header=1).T

    r1, = ax.plot(np.log10(r1_wave), np.log10(r1_sed), 'o-',mfc='Magenta',mec='Magenta',markersize=5,markeredgewidth=1,color='Magenta', linewidth=1.5)
    r2, = ax.plot(np.log10(r2_wave), np.log10(r2_sed), 'o-',mfc='r',mec='r',markersize=5,markeredgewidth=1,color='r', linewidth=1.5)
    r3, = ax.plot(np.log10(r3_wave), np.log10(r3_sed), 'o-',mfc='b',mec='b',markersize=5,markeredgewidth=1,color='b', linewidth=1.5)

    if obs != None:  
        import sys
        sys.path.append('/Users/yaolun/programs/misc/hyperion')
        from get_bhr71_obs import get_bhr71_obs
        c = const.c.cgs.value

        bhr71 = get_bhr71_obs(obs)  # in um and Jy
        wave_obs, flux_obs, noise_obs = bhr71['spec']
        obs_data, = ax.plot(np.log10(wave_obs[wave_obs<50]), np.log10(c/(wave_obs[wave_obs<50]*1e-4)*flux_obs[wave_obs<50]*1e-23), color='k', alpha=0.7, linewidth=1)
        ax.plot(np.log10(wave_obs[(wave_obs>50)&(wave_obs<190.31)]), np.log10(c/(wave_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-4)*flux_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-23), color='k', alpha=0.7, linewidth=1)
        ax.plot(np.log10(wave_obs[wave_obs>194]), np.log10(c/(wave_obs[wave_obs>194]*1e-4)*flux_obs[wave_obs>194]*1e-23), color='k', alpha=0.7, linewidth=1)


    [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
    ax.minorticks_on() 
    ax.tick_params('both',labelsize=14,width=1.5,which='major',pad=15,length=5)
    ax.tick_params('both',labelsize=14,width=1.5,which='minor',pad=15,length=2.5)

    ax.set_xlabel(r'$\mathrm{log(wavelength)~(\mu m)}$', fontsize=16)
    ax.set_ylabel(r'$\mathrm{log~\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$', fontsize=16)
    ax.set_ylim([-13,-7])

    plt.legend([r1, r2, r3, obs_data], [r'$\mathrm{\rho(r)\propto r^{-2}}$', r'$\mathrm{\rho(r)\propto r^{-1.5}}$',r'$\mathrm{const.+r^{-2}}$',r'$\mathrm{observation}$'], numpoints=1, loc='lower right', fontsize=16)

    fig.savefig(outdir+'sed_cav_struc_com.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()


import numpy as np
indir = '/Users/yaolun/bhr71/hyperion/controlled/'
outdir = '/Users/yaolun/Copy/Papers/yaolun/bhr71/figures/'
obs = '/Users/yaolun/bhr71/obs_for_radmc/'
# obs = None

# # grid of cs and age
# array = np.array([[1,6,11],[2,7,12],[3,8,13],[4,9,14],[5,10,15]])
# array = np.array([[1,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15]])
# cslist = [0.1,0.2,0.3]
# agelist = [1e4,2.5e4,5e4,7.5e4,1e5]
# sed_grid_cs_age(indir, array, outdir, cslist, agelist, obs= None)

# grid of Omega0
array = np.array([18,17,16])
sed_omega(indir, array, outdir, obs= None, compact=True)

# # grid of disk parameters
# # disk mass
# array = np.array([19,20,21,22,23])
# xlabel = r'$\mathrm{M_{disk}~[M_{\odot}]~(0.1,~0.3,~0.5,~0.7,~1.0)}$'
# compact = [r'$\mathrm{M_{disk}=0.1~M_{\odot}}$',r'$\mathrm{M_{disk}=0.3~M_{\odot}}$',r'$\mathrm{M_{disk}=0.5~M_{\odot}}$',r'$\mathrm{M_{disk}=0.7~M_{\odot}}$',r'$\mathrm{M_{disk}=1.0~M_{\odot}}$']
# plotname = 'disk_mdisk'
# sed_five(indir, array, outdir, xlabel, plotname, obs= None, zoom=True, compact=compact, yrange=[-13,-8])
# # flare power
# array = np.array([29,30,31,32,33])
# xlabel = r'$\mathrm{\beta~(1.0,~1.2,~1.4,~1.6,~1.8)}$'
# compact = [r'$\mathrm{\beta=1.0}$',r'$\mathrm{\beta=1.2}$',r'$\mathrm{\beta=1.4}$',r'$\mathrm{\beta=1.6}$',r'$\mathrm{\beta=1.8}$']
# plotname = 'disk_beta'
# sed_five(indir, array, outdir, xlabel, plotname, obs= None, zoom=True, compact=compact, yrange=[-13,-8])

# # grid of theta_cav and incl.
# array = np.array([[35,36,37,38,39],[40,41,42,43,44],[45,46,47,48,49]])
# # sed_grid_theta_cav_incl(indir, array, outdir, obs= None)
# sed_grid_theta_cav_incl(indir, array, outdir, obs= None, compact=True)

# # grid of rho_cav_center and sed_rho_cav_edge
# array = np.array([[49,50,51,52],[53,54,55,56],[57,58,59,60]])
# # sed_grid_rho_cav_centeredge(indir, array, outdir, obs= None)
# sed_grid_rho_cav_centeredge(indir, array, outdir, obs= None, compact=True)

# # # disk & no dis comparison
# array = np.array([16,61])
# disk_exist_com(indir, array, outdir, obs=obs)

# # grid of tstar
# array = np.array([69,70,71])
# sed_tstar(indir, array, outdir, obs=obs)

# # grid of rstar
# array = np.array([72,73,74])
# sed_rstar(indir, array, outdir, obs=obs)

# grid of R_env_max
array = np.array([63,64,65,66,67])
xlabel = r'$\mathrm{R_{env,max}~[AU]~(7.5\times 10^{3},~1\times 10^{4},~2.5\times 10^{4},~5\times 10^{4},~7.5\times 10^{4})}$'
compact = [r'$\mathrm{R_{env,max}=7.5\times 10^{3}~AU}$',r'$\mathrm{R_{env,max}=1.0\times 10^{4}~AU}$',r'$\mathrm{R_{env,max}=2.5\times 10^{4}~AU}$',r'$\mathrm{R_{env,max}=5.0\times 10^{4}~AU}$',r'$\mathrm{R_{env,max}=7.5\times 10^{4}~AU}$']
plotname = 'r_max'
sed_five(indir, array, outdir, xlabel, plotname, obs= None, tbol=True, compact=compact)

# # grid of continuous cavity power law
# array = np.array([13,14,15])
# sed_cav_powerlaw('/Users/yaolun/bhr71/hyperion/cycle5', array, outdir, obs=obs)

# # grid of cavity structure comparison
# array = np.array([14,17,1])
# sed_cav_struc_com('/Users/yaolun/bhr71/hyperion/cycle5', array, outdir, obs=obs)
