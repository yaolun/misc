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
            ax.text(0.4, 0.1, r'$T_{bol}= %4.1f\/K$' % t_bol(wave,sed*wave*1e-4/c), fontsize=12, transform=ax.transAxes)

            ax.set_ylim([-15,-8])

            [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
            ax.minorticks_on() 
            ax.tick_params('both',labelsize=14,width=1,which='major',pad=15,length=5)
            ax.tick_params('both',labelsize=14,width=1,which='minor',pad=15,length=2.5)

            if rr+1 == row:
                if cc == 0:
                    ax.set_xlabel(r'$log(wavelength)\/(\mu m)$', fontsize=14)
                    ax.set_ylabel(r'$log\/\nu S_{\nu}\/(erg\/s^{-1}\/cm^{-2})$', fontsize=14)
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
            #         ax.set_xlabel(r'$\rm{radius\/[AU]}$', fontsize=14)
            #         ax.set_ylabel(r'$\rm{dust\/density\/[g\/cm^{-3}]}$', fontsize=14)

            # fix the overlap tick labels
            x_nbins = len(ax.get_xticklabels())
            y_nbins = len(ax.get_yticklabels())
            if (rr != 0) & (cc != 0):
                ax.xaxis.set_major_locator(MaxNLocator(nbins=x_nbins, prune='lower'))
                ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune='upper'))

    fig.text(0.5, -0.05 , r'$age\/[yr]\/(1\times 10^{4},\/2.5\times 10^{4},\/5\times 10^{4},\/7.5\times 10^{4},\/1\times 10^{5})$', fontsize=20, ha='center')
    fig.text(0, 0.5, r'$sound\/speed\/[km\/s^{-1}]\/(0.38,\/0.25,\/0.2)$', fontsize=20, va='center', rotation='vertical')

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
                ax.set_xlabel(r'$log(wavelength)\/(\mu m)}$', fontsize=14)
                ax.set_ylabel(r'$log\/\nu S_{\nu}\/(erg\/s^{-1}\/cm^{-2})$', fontsize=14)

            # fix the overlap tick labels
            x_nbins = len(ax.get_xticklabels())
            y_nbins = len(ax.get_yticklabels())
            if (cc != 0):
                ax.xaxis.set_major_locator(MaxNLocator(nbins=x_nbins, prune='lower'))
                ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune='upper'))

        fig.text(0.5, -0.13 , r'$\Omega_{\circ}\/[s^{-1}]\/(1\times 10^{-13},\/5\times 10^{-14},\/1\times 10^{-14})$', fontsize=14, ha='center')

        fig.subplots_adjust(hspace=0,wspace=0)
        fig.savefig(outdir+'sed_omega0.pdf', format='pdf', dpi=300, bbox_inches='tight')
        fig.clf()
    else:
        # color_list = ['#ffeda0','#feb24c','#f03b20']
        # color_list = ['#fd8d3c','#f03b20','#bd0026']
        color_list = [[0.8507598215729224, 0.6322174528970308, 0.6702243543099417],\
                     [0.5687505862870377, 0.3322661256969763, 0.516976691731939],\
                     [0.1750865648952205, 0.11840023306916837, 0.24215989137836502]]

        label = [r'$\Omega_{\circ}=1\times 10^{-14}\/s^{-1}$',r'$\Omega_{\circ}=5\times 10^{-14}\/s^{-1}$',r'$\Omega_{\circ}=1\times 10^{-13}\/s^{-1}$']
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)

        for i in range(0, len(array)):
            # sed with infinite aperture
            (wave_inf, sed_inf) = np.genfromtxt(indir+'/model'+str(array[i])+'_sed_inf.txt', skip_header=1).T
            # sed with apertures
            (wave, sed) = np.genfromtxt(indir+'/model'+str(array[i])+'_sed_w_aperture.txt', skip_header=1).T

            ax.plot(np.log10(wave), np.log10(sed), 'o-',mfc=color_list[i],mec=color_list[i],color=color_list[i],markersize=7,markeredgewidth=1,linewidth=2,label=label[i])
        ax.legend(loc='lower right', numpoints=1, framealpha=0.3, fontsize=16)

        [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
        ax.minorticks_on() 
        ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
        ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)

        ax.set_xlabel(r'$log(wavelength)\/(\mu m)$', fontsize=20)
        ax.set_ylabel(r'$log\/\nu S_{\nu}\/(erg\/s^{-1}\/cm^{-2})$', fontsize=20)
        ax.set_ylim([-13,-8])
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

            ax.plot(np.log10(wave), np.log10(sed), 'o-',mfc='b',mec='b',markersize=3,markeredgewidth=1,linewidth=1.2)

            if tbol == True:
                 ax.text(0.4, 0.1, r'$T_{bol}= %4.1f\/K$' % t_bol(wave, sed*wave*1e-4/c), fontsize=12, transform=ax.transAxes)

            if zoom == True:
                ax.set_xlim([0.4,2])

            [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
            ax.minorticks_on() 
            ax.tick_params('both',labelsize=10,width=1.2,which='major',pad=15,length=5)
            ax.tick_params('both',labelsize=10,width=1.2,which='minor',pad=15,length=2.5)

            if cc == 0:
                ax.set_xlabel(r'$log(wavelength)\/(\mu m)$', fontsize=14)
                ax.set_ylabel(r'$log\/\nu S_{\nu}\/(erg\/s^{-1}\/cm^{-2})$', fontsize=14)

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
        color_list = ['#ecd078','#d95b43','#c02942','#542437','#53777a']
        color_list = [[0.8507598215729224, 0.6322174528970308, 0.6702243543099417],\
                     [0.7357826498167007, 0.4722583075098643, 0.5939898816486836],\
                     [0.5687505862870377, 0.3322661256969763, 0.516976691731939],\
                     [0.37315740206144277, 0.21948554297767336, 0.40755444345087455],\
                     [0.1750865648952205, 0.11840023306916837, 0.24215989137836502]]
        if len(array) == 3:
            color_list = [color_list[0],color_list[2],color_list[4]]

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

            ax.plot(np.log10(wave), np.log10(sed), 'o-',mfc=color_list[i],mec=color_list[i],color=color_list[i],markersize=7,markeredgewidth=1,linewidth=2,label=compact[i])
            ax.legend(loc='lower right', numpoints=1, framealpha=0.3, fontsize=16)
        ax.set_xlabel(r'$log(wavelength)\/(\mu m)$', fontsize=18)
        ax.set_ylabel(r'$log\/\nu S_{\nu}\/(erg\/s^{-1}\/cm^{-2})$', fontsize=18)
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
                        ax.set_xlabel(r'$log(wavelength)\/(\mu m)$', fontsize=14)
                        ax.set_ylabel(r'$log\/\nu S_{\nu}\/(erg\/s^{-1}\/cm^{-2})$', fontsize=14)
            
                # fix the overlap tick labels
                x_nbins = len(ax.get_xticklabels())
                y_nbins = len(ax.get_yticklabels())
                if (rr != 0) & (cc != 0):
                    ax.xaxis.set_major_locator(MaxNLocator(nbins=x_nbins, prune='lower'))
                    ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune='upper'))

        fig.text(0.5, -0.05 , r'$\theta_{\rm cav}\/[deg.]\/(15^{\circ},\/20^{\circ},\/25^{\circ},\/30^{\circ},\/35^{\circ})$', fontsize=20, ha='center')
        fig.text(0, 0.5, r'$\theta_{\rm incl.}\/[deg.]\/(60^{\circ},\/70^{\circ},\/8 0^{\circ})$', fontsize=20, va='center', rotation='vertical')

        fig.subplots_adjust(hspace=0,wspace=0)
        fig.savefig(outdir+'sed_theta_cav_incl.pdf', format='pdf', dpi=300, bbox_inches='tight')
        fig.clf()
    else:
        # plot multiple SEDs in same panel
        num, col = np.shape(array)
        row = 1

        # set up the color map
        # color_list = ['#f03b20','#feb24c','#ffeda0']
        # color_list = ['#ffeda0','#feb24c','#f03b20']
        color_list = ['#fd8d3c','#f03b20','#bd0026']
        color_list = [[0.8507598215729224, 0.6322174528970308, 0.6702243543099417],\
                     [0.5687505862870377, 0.3322661256969763, 0.516976691731939],\
                     [0.1750865648952205, 0.11840023306916837, 0.24215989137836502]]

        # label setup
        label = [ r'$\theta_{incl.}=80^{\circ}$', r'$\theta_{incl.}=65^{\circ}$',\
                  r'$\theta_{incl.}=50^{\circ}$']

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

                ax.plot(np.log10(wave), np.log10(sed), 'o-',mfc=color_list[i],mec=color_list[i], color=color_list[i], markersize=3,markeredgewidth=1,linewidth=1.2, label=label[i])
                if cc == 0:
                    ax.legend(loc='upper left', numpoints=1,framealpha=0.3,fontsize=10)

            ax.set_ylim([-14,-8])
            ax.set_xlim([0.4,2])
            ax.locator_params(axis='x', nbins=5)

            [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
            ax.minorticks_on() 
            ax.tick_params('both',labelsize=14,width=1.2,which='major',pad=10,length=5)
            ax.tick_params('both',labelsize=14,width=1.2,which='minor',pad=10,length=2.5)

            if cc == 0:
                ax.set_xlabel(r'$log(wavelength)\/(\mu m)$', fontsize=16)
                ax.set_ylabel(r'$log\/\nu S_{\nu}\/(erg\/s^{-1}\/cm^{-2})$', fontsize=16)

            # fix the overlap tick labels
            x_nbins = len(ax.get_xticklabels())
            y_nbins = len(ax.get_yticklabels())
            if (cc != 0):
                ax.xaxis.set_major_locator(MaxNLocator(nbins=x_nbins, prune='lower'))
                ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune='upper'))

        # fig.text(0.5, -0.05 , r'$\rho_{cav,\circ}\/[g\/cm^{-3}]\/(1\times 10^{-20},\/5\times 10^{-20},\/1\times 10^{-19},\/5\times 10^{-19})$', fontsize=20, ha='center')
        fig.text(0.5, -0.15, r'$\theta_{cav}\/[deg.]\/(15^{\circ},\/20^{\circ},\/25^{\circ},\/30^{\circ},\/35^{\circ})$', fontsize=20, ha='center' )

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
                        ax.set_xlabel(r'$log(wavelength)\/(\mu m)$', fontsize=14)
                        ax.set_ylabel(r'$log\/\nu S_{\nu}\/(erg\/s^{-1}\/cm^{-2})$', fontsize=14)
            
                # fix the overlap tick labels
                x_nbins = len(ax.get_xticklabels())
                y_nbins = len(ax.get_yticklabels())
                if (rr != 0) & (cc != 0):
                    ax.xaxis.set_major_locator(MaxNLocator(nbins=x_nbins, prune='lower'))
                    ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune='upper'))

        fig.text(0.5, -0.05 , r'$\rho_{cav,\circ}\/[g\/cm^{-3}]\/(1\times 10^{-20},\/5\times 10^{-20},\/1\times 10^{-19},\/5\times 10^{-19})$', fontsize=20, ha='center')
        fig.text(0, 0.5, r'$R_{cav,\circ}\/[AU]\/(40,\/30,\/20)$', fontsize=20, va='center', rotation='vertical')

        fig.subplots_adjust(hspace=0,wspace=0)
        fig.savefig(outdir+'sed_rho_cav_centeredge.pdf', format='pdf', dpi=300, bbox_inches='tight')
        fig.clf()
    else:
        # plot multiple SEDs in same panel
        col, num = np.shape(array)
        row = 1

        # set up the color map
        color_list = plt.cm.gnuplot2(np.linspace(0, 1, num+1))
        color_list = [[0.8507598215729224, 0.6322174528970308, 0.6702243543099417],\
                     [0.7357826498167007, 0.4722583075098643, 0.5939898816486836],\
                     [0.5687505862870377, 0.3322661256969763, 0.516976691731939],\
                     [0.37315740206144277, 0.21948554297767336, 0.40755444345087455],\
                     [0.1750865648952205, 0.11840023306916837, 0.24215989137836502]]
        color_list = color_list[len(color_list)-num:]

        # label setup
        label = [ r'$5\times 10^{-20}$', r'$1\times 10^{-19}$',\
                  r'$5\times 10^{-19}$', r'$1\times 10^{-18}$']

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

                ax.plot(np.log10(wave), np.log10(sed), 'o-',mfc=color_list[i],mec=color_list[i], color=color_list[i], markersize=4,markeredgewidth=1,linewidth=1.5, label=label[i])
                if cc == col-1:
                    ax.legend(loc='lower right', numpoints=1,framealpha=0.3,fontsize=12)

            ax.set_ylim([-14,-8])
            ax.set_xlim([0.4,2])
            ax.locator_params(axis='x', nbins=5)

            [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
            ax.minorticks_on() 
            ax.tick_params('both',labelsize=16,width=1.2,which='major',pad=10,length=5)
            ax.tick_params('both',labelsize=16,width=1.2,which='minor',pad=10,length=2.5)

            if cc == 0:
                ax.set_xlabel(r'$log(wavelength)\/(\mu m)$', fontsize=18)
                ax.set_ylabel(r'$log\/\nu S_{\nu}\/(erg\/s^{-1}\/cm^{-2})$', fontsize=18)

            # fix the overlap tick labels
            x_nbins = len(ax.get_xticklabels())
            y_nbins = len(ax.get_yticklabels())
            if (cc != 0):
                ax.xaxis.set_major_locator(MaxNLocator(nbins=x_nbins, prune='lower'))
                ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune='upper'))

        # fig.text(0.5, -0.05 , r'$\rho_{cav,\circ}\/[g\/cm^{-3}]\/(1\times 10^{-20},\/5\times 10^{-20},\/1\times 10^{-19},\/5\times 10^{-19})$', fontsize=20, ha='center')
        fig.text(0.5, -0.15, r'$R_{cav,\circ}\/[AU]\/(20,\/30,\/40)$', fontsize=20, ha='center' )

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

    disk, = ax.plot(np.log10(d_wave), np.log10(d_sed), 'o-',mfc='b',mec='b',markersize=7,markeredgewidth=1,color='b', linewidth=2)
    nodisk, = ax.plot(np.log10(nd_wave), np.log10(nd_sed), 'o-',mfc='k',mec='k',markersize=7,markeredgewidth=1,color='k', linewidth=2)

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
    ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
    ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)

    ax.set_xlabel(r'$log(wavelength)\/(\mu m)$', fontsize=20)
    ax.set_ylabel(r'$log\/\nu S_{\nu}\/(erg\/s^{-1}\/cm^{-2})$', fontsize=20)

    plt.legend([disk, nodisk, obs_data], [r'$w/\/disk$', r'$w/o\/disk$',r'$observation$'], numpoints=1, loc='lower right', fontsize=16)

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

    ax.set_xlabel(r'$log(wavelength)\/(\mu m)$', fontsize=16)
    ax.set_ylabel(r'$log\/\nu S_{\nu}\/(erg\/s^{-1}\/cm^{-2})$', fontsize=16)
    ax.set_ylim([-13,-7])

    plt.legend([t1, t2, t3, obs_data], [r'$T_{\star}=4500\/K$', r'$T_{\star}=5000\/K$',r'$T_{\star}=5500\/K$',r'$observation$'], numpoints=1, loc='lower right', fontsize=16)

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

    ax.set_xlabel(r'$log(wavelength)\/(\mu m)$', fontsize=16)
    ax.set_ylabel(r'$log\/\nu S_{\nu}\/(erg\/s^{-1}\/cm^{-2})$', fontsize=16)
    ax.set_ylim([-13,-7])

    plt.legend([r1, r2, r3, obs_data], [r'$R_{\star}=0.3\/R_{\odot}$', r'$R_{\star}=0.5\/R_{\odot}$',r'$R_{\star}=0.7\/R_{\odot}$',r'$observation$'], numpoints=1, loc='lower right', fontsize=16)

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

    r1, = ax.plot(np.log10(r1_wave), np.log10(r1_sed), 'o-',mfc='Magenta',mec='Magenta',markersize=7,markeredgewidth=1,color='Magenta', linewidth=2)
    r2, = ax.plot(np.log10(r2_wave), np.log10(r2_sed), 'o-',mfc='r',mec='r',markersize=7,markeredgewidth=1,color='r', linewidth=2)
    r3, = ax.plot(np.log10(r3_wave), np.log10(r3_sed), 'o-',mfc='b',mec='b',markersize=7,markeredgewidth=1,color='b', linewidth=2)

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
    ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
    ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)

    ax.set_xlabel(r'$log(wavelength)\/(\mu m)$', fontsize=20)
    ax.set_ylabel(r'$log\/\nu S_{\nu}\/(erg\/s^{-1}\/cm^{-2})$', fontsize=20)
    ax.set_ylim([-13,-7])

    plt.legend([r1, r2, r3, obs_data], [r'$\rho_{cav,\circ}=5\times 10^{-15}\/g\/cm^{-3}$', r'$\rho_{cav,\circ}=5\times 10^{-16}\/g\/cm^{-3}$',r'$\rho_{cav,\circ}=5\times 10^{-17}\/g\/cm^{-3}$',r'$observation$'], numpoints=1, loc='lower right', fontsize=16)

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
    # rho(r) \/ r^-2
    (r1_wave_inf, r1_sed_inf) = np.genfromtxt(indir+'/model'+str(array[0])+'_sed_inf.txt', skip_header=1).T
    (r1_wave, r1_sed) = np.genfromtxt(indir+'/model'+str(array[0])+'_sed_w_aperture.txt', skip_header=1).T

    # rho(r) \/ r^-1.5
    (r2_wave_inf, r2_sed_inf) = np.genfromtxt(indir+'/model'+str(array[1])+'_sed_inf.txt', skip_header=1).T
    (r2_wave, r2_sed) = np.genfromtxt(indir+'/model'+str(array[1])+'_sed_w_aperture.txt', skip_header=1).T

    # rho(r) \/ cont.s + r^-2
    (r3_wave_inf, r3_sed_inf) = np.genfromtxt(indir+'/model'+str(array[2])+'_sed_inf.txt', skip_header=1).T
    (r3_wave, r3_sed) = np.genfromtxt(indir+'/model'+str(array[2])+'_sed_w_aperture.txt', skip_header=1).T

    # rho(r) \/ uniform
    (r4_wave_inf, r4_sed_inf) = np.genfromtxt(indir+'/model'+str(array[3])+'_sed_inf.txt', skip_header=1).T
    (r4_wave, r4_sed) = np.genfromtxt(indir+'/model'+str(array[3])+'_sed_w_aperture.txt', skip_header=1).T

    r1, = ax.plot(np.log10(r1_wave), np.log10(r1_sed), 'o-',mfc='Magenta',mec='Magenta',markersize=7,markeredgewidth=1,color='Magenta', linewidth=2)
    r2, = ax.plot(np.log10(r2_wave), np.log10(r2_sed), 'o-',mfc='r',mec='r',markersize=7,markeredgewidth=1,color='r', linewidth=2)
    r3, = ax.plot(np.log10(r3_wave), np.log10(r3_sed), 'o-',mfc='b',mec='b',markersize=7,markeredgewidth=1,color='b', linewidth=2)
    r4, = ax.plot(np.log10(r4_wave), np.log10(r4_sed), 'o-',mfc='k',mec='k',markersize=7,markeredgewidth=1,color='k', linewidth=2)
 
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
    ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
    ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)

    ax.set_xlabel(r'$log(wavelength)\/(\mu m)$', fontsize=20)
    ax.set_ylabel(r'$log\/\nu S_{\nu}\/(erg\/s^{-1}\/cm^{-2})$', fontsize=20)
    ax.set_ylim([-13,-7])

    plt.legend([r1, r2, r3, r4, obs_data], [r'$\rho(r)\propto r^{-2}$', r'$\rho(r)\propto r^{-1.5}$',\
                r'$const.+r^{-2}$',r'$uniform$',r'$observation$'],\
                numpoints=1, loc='lower right', fontsize=16)
    # plt.legend([r1, r2, r4, obs_data], [r'$\rho(r)\propto r^{-2}$', r'$\rho(r)\propto r^{-1.5}$',\
    #             r'$uniform$',r'$observation$'],\
    #             numpoints=1, loc='lower right', fontsize=16)

    fig.savefig(outdir+'sed_cav_struc_com.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()

def sed_lum(indir, array, outdir, obs=None):
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

    ax.set_xlabel(r'$log(wavelength)\/(\mu m)$', fontsize=16)
    ax.set_ylabel(r'$log\/\nu S_{\nu}\/(erg\/s^{-1}\/cm^{-2})$', fontsize=16)
    ax.set_ylim([-13,-7])

    plt.legend([r1, r2, r3], [r'$4500\/K$', r'$5000\/K$',r'$5500\/K$'], numpoints=1, loc='lower right', fontsize=16)

    fig.savefig(outdir+'sed_lstar.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()

def model_vs_obs(modelname,indir,outdir,obs=None,dstar=178.0,wl_aper=None,rtout=False):
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    from hyperion.model import ModelOutput
    from hyperion.model import Model
    from scipy.interpolate import interp1d
    from hyperion.util.constants import pc, c, lsun    
    import sys
    sys.path.append('/Users/yaolun/programs/misc/hyperion')
    from get_bhr71_obs import get_bhr71_obs
    from l_bol import l_bol    

    # Read in the observation data and calculate the noise & variance
    if indir == None:
        indir = '/Users/yaolun/bhr71/'
    if outdir == None:
        outdir = '/Users/yaolun/bhr71/hyperion/'

    bhr71 = get_bhr71_obs(obs)  # in um and Jy
    wave_obs, flux_obs, noise_obs = bhr71['spec']
    wave_phot, flux_phot, flux_sig_phot = bhr71['phot']

    # Convert the unit from Jy to erg cm-2 Hz-1
    flux_obs = flux_obs*1e-23
    noise_obs = noise_obs*1e-23
    flux_phot = flux_phot*1e-23
    flux_sig_phot = flux_sig_phot*1e-23

    # Print the observed L_bol
    wl_tot = np.hstack((wave_obs,wave_phot))
    flux_tot = np.hstack((flux_obs,flux_phot))
    flux_tot = flux_tot[np.argsort(wl_tot)]
    wl_tot = wl_tot[np.argsort(wl_tot)]
    l_bol_obs = l_bol(wl_tot,flux_tot*1e23, dstar)             

    # Open the model
    m = ModelOutput(indir+modelname+'.rtout')

    if wl_aper == None:
        wl_aper = [3.6, 4.5, 5.8, 8.0, 8.5, 9, 9.7, 10, 10.5, 11, 16, 20, 24, 35, 70, 100, 160, 250, 350, 500, 850]

    # Create the plot
    mag = 1.5
    fig = plt.figure(figsize=(8*mag,6*mag))
    ax_sed = fig.add_subplot(1, 1, 1)

    # Plot the observed SED
    pacs, = ax_sed.plot(np.log10(wave_obs[wave_obs<50]), np.log10(c/(wave_obs[wave_obs<50]*1e-4)*flux_obs[wave_obs<50]),'-',color='DimGray', alpha=0.7, linewidth=1.5*mag)
    spire, = ax_sed.plot(np.log10(wave_obs[(wave_obs>50)&(wave_obs<190.31)]), np.log10(c/(wave_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-4)*flux_obs[(wave_obs>50)&(wave_obs<190.31)]),'-',color='DimGray', alpha=0.7, linewidth=1.5*mag)
    irs, = ax_sed.plot(np.log10(wave_obs[wave_obs>194]), np.log10(c/(wave_obs[wave_obs>194]*1e-4)*flux_obs[wave_obs>194]),'-',color='DimGray', alpha=0.7, linewidth=1.5*mag)


    # plot the observed photometry data
    photometry, = ax_sed.plot(np.log10(wave_phot),np.log10(c/(wave_phot*1e-4)*flux_phot),'s',mfc='DimGray',mec='k',markersize=8)
    ax_sed.errorbar(np.log10(wave_phot),np.log10(c/(wave_phot*1e-4)*flux_phot),\
        yerr=[np.log10(c/(wave_phot*1e-4)*flux_phot)-np.log10(c/(wave_phot*1e-4)*(flux_phot-flux_sig_phot)),\
              np.log10(c/(wave_phot*1e-4)*(flux_phot+flux_sig_phot))-np.log10(c/(wave_phot*1e-4)*flux_phot)],\
        fmt='s',mfc='DimGray',mec='k',markersize=8)

    # perform the same procedure of flux extraction of aperture flux with observed spectra
    wl_aper = np.array(wl_aper)
    obs_aper_wl = wl_aper[(wl_aper >= min(wave_obs)) & (wl_aper <= max(wave_obs))]
    obs_aper_sed = np.empty_like(obs_aper_wl)
    sed_tot = c/(wl_tot*1e-4)*flux_tot
    # wl_tot and flux_tot are already hstacked and sorted by wavelength
    for i in range(0, len(obs_aper_wl)):
        if (obs_aper_wl[i] < 50.) & (obs_aper_wl[i] >= 5):
            res = 60.
        elif obs_aper_wl[i] < 5:
            res = 10.
        else:
            res = 1000.
        ind = np.where((wl_tot < obs_aper_wl[i]*(1+1./res)) & (wl_tot > obs_aper_wl[i]*(1-1./res)))
        if len(ind[0]) != 0:
            obs_aper_sed[i] = np.mean(sed_tot[ind])
        else:
            f = interp1d(wl_tot, sed_tot)
            obs_aper_sed[i] = f(wl_aper[i])

    # read in the raw output from hyperion and then save the extraction results in text file, otherwise read the text files instead.
    if rtout == True:
        sed_inf = m.get_sed(group=0, inclination=0, aperture=-1, distance=dstar * pc)
        flux_aper = np.empty_like(wl_aper)
        unc_aper = np.empty_like(wl_aper)
        for i in range(0, len(wl_aper)):
            sed_dum = m.get_sed(group=i+1, inclination=0, aperture=-1, distance=dstar * pc)
            # use a rectangle function the average the simulated SED
            # apply the spectral resolution
            if (wl_aper[i] < 50.) & (wl_aper[i] >= 5):
                res = 60.
            elif wl_aper[i] < 5:
                res = 10.
            else:
                res = 1000.
            ind = np.where((sed_dum.wav < wl_aper[i]*(1+1./res)) & (sed_dum.wav > wl_aper[i]*(1-1./res)))
            if len(ind[0]) != 0:
                flux_aper[i] = np.mean(sed_dum.val[ind])
            else:
                f = interp1d(sed_dum.wav, sed_dum.val)
                flux_aper[i] = f(wl_aper[i])
        
        # make the variables consistent with others
        sim_inf = sed_inf.wav
        sim_sed_inf = sed_inf.val

        # save the results in text files
        # unapertured SED
        foo = open(outdir+modelname+'_sed_inf.txt','w')
        foo.write('%12s \t %12s \n' % ('wave','vSv'))
        for i in range(0, len(sed_inf.wav)):
            foo.write('%12g \t %12g \n' % (sed_inf.wav[i], sed_inf.val[i]))
        foo.close()
        # SED with convolution of aperture sizes
        foo = open(outdir+modelname+'_sed_w_aperture.txt','w')
        foo.write('%12s \t %12s \n' % ('wave','vSv'))
        for i in range(0, len(wl_aper)):
            foo.write('%12g \t %12g \n' % (wl_aper[i], flux_aper[i]))
        foo.close()
    else:
        # read in the extracted text files
        (sim_inf, sim_sed_inf) = np.genfromtxt(indir+modelname+'_sed_inf.txt', skip_header=1).T
        (wl_aper, flux_aper) = np.genfromtxt(indir+modelname+'_sed_w_aperture.txt', skip_header=1).T

    aper_obs, = ax_sed.plot(np.log10(obs_aper_wl),np.log10(obs_aper_sed), 's-', mec='None', mfc='r', color='r',markersize=10, linewidth=1.5)
    aper, = ax_sed.plot(np.log10(wl_aper),np.log10(flux_aper),'o-', mec='Blue', mfc='None', color='b',markersize=12, markeredgewidth=3, linewidth=1.7)
    # calculate the bolometric luminosity of the aperture 
    l_bol_sim = l_bol(wl_aper, flux_aper/(c/np.array(wl_aper)*1e4)*1e23, dstar)
    print 'Bolometric luminosity of simulated spectrum: %5.2f lsun' % l_bol_sim

    # plot setting
    ax_sed.set_xlabel(r'$\rm{log\,\lambda\,({\mu}m)}$',fontsize=mag*20)
    ax_sed.set_ylabel(r'$\rm{log\,\nu S_{\nu}\,(erg\,cm^{-2}\,s^{-1})}$',fontsize=mag*20)
    [ax_sed.spines[axis].set_linewidth(1.5*mag) for axis in ['top','bottom','left','right']]
    ax_sed.minorticks_on()
    ax_sed.tick_params('both',labelsize=mag*18,width=1.5*mag,which='major',pad=15,length=5*mag)
    ax_sed.tick_params('both',labelsize=mag*18,width=1.5*mag,which='minor',pad=15,length=2.5*mag)

    ax_sed.set_ylim([-13,-7.5])
    ax_sed.set_xlim([0,3])

    lg_data = ax_sed.legend([irs, photometry, aper, aper_obs],\
        [r'$\rm{observation}$',\
        r'$\rm{photometry}$',r'$\rm{F_{aper,sim}}$',r'$\rm{F_{aper,obs}}$'],\
        loc='upper left',fontsize=14*mag,numpoints=1,framealpha=0.3)

    # Write out the plot
    fig.savefig(outdir+modelname+'_sed.pdf',format='pdf',dpi=300,bbox_inches='tight')
    fig.clf()

def three_model_vs_obs(modelname,indir,outdir,label,obs=None,dstar=178.0,wl_aper=None,rtout=False):
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    from hyperion.model import ModelOutput
    from hyperion.model import Model
    from scipy.interpolate import interp1d
    from hyperion.util.constants import pc, c, lsun    
    import sys
    sys.path.append('/Users/yaolun/programs/misc/hyperion')
    from get_bhr71_obs import get_bhr71_obs
    from l_bol import l_bol    

    color_list = ['CornflowerBlue','SkyBlue','Blue']
    style = ['-','--','-.']
    style = ['-','-','-']

    # Read in the observation data and calculate the noise & variance
    if indir == None:
        indir = '/Users/yaolun/bhr71/'
    if outdir == None:
        outdir = '/Users/yaolun/bhr71/hyperion/'

    bhr71 = get_bhr71_obs(obs)  # in um and Jy
    wave_obs, flux_obs, noise_obs = bhr71['spec']
    wave_phot, flux_phot, flux_sig_phot = bhr71['phot']

    # Convert the unit from Jy to erg cm-2 Hz-1
    flux_obs = flux_obs*1e-23
    noise_obs = noise_obs*1e-23
    flux_phot = flux_phot*1e-23
    flux_sig_phot = flux_sig_phot*1e-23

    # Print the observed L_bol
    wl_tot = np.hstack((wave_obs,wave_phot))
    flux_tot = np.hstack((flux_obs,flux_phot))
    flux_tot = flux_tot[np.argsort(wl_tot)]
    wl_tot = wl_tot[np.argsort(wl_tot)]
    l_bol_obs = l_bol(wl_tot,flux_tot*1e23, dstar)             


    if wl_aper == None:
        wl_aper = [3.6, 4.5, 5.8, 8.0, 8.5, 9, 9.7, 10, 10.5, 11, 16, 20, 24, 35, 70, 100, 160, 250, 350, 500, 850]

    # Create the plot
    mag = 1.5
    fig = plt.figure(figsize=(8*mag,6*mag))
    ax_sed = fig.add_subplot(1, 1, 1)

    # Plot the observed SED
    ax_sed.plot(np.log10(wave_obs[wave_obs<50]), np.log10(c/(wave_obs[wave_obs<50]*1e-4)*flux_obs[wave_obs<50]),'-',color='DimGray', alpha=0.7, linewidth=1.5*mag)
    ax_sed.plot(np.log10(wave_obs[(wave_obs>50)&(wave_obs<190.31)]), np.log10(c/(wave_obs[(wave_obs>50)&(wave_obs<190.31)]*1e-4)*flux_obs[(wave_obs>50)&(wave_obs<190.31)]),'-',color='DimGray', alpha=0.7, linewidth=1.5*mag)
    ax_sed.plot(np.log10(wave_obs[wave_obs>194]), np.log10(c/(wave_obs[wave_obs>194]*1e-4)*flux_obs[wave_obs>194]),'-',color='DimGray', alpha=0.7, linewidth=1.5*mag, label=r'$\rm{observations}$')


    # plot the observed photometry data
    ax_sed.plot(np.log10(wave_phot),np.log10(c/(wave_phot*1e-4)*flux_phot),'s',mfc='DimGray',mec='k',markersize=8,label=r'$\rm{photometry}$')
    ax_sed.errorbar(np.log10(wave_phot),np.log10(c/(wave_phot*1e-4)*flux_phot),\
        yerr=[np.log10(c/(wave_phot*1e-4)*flux_phot)-np.log10(c/(wave_phot*1e-4)*(flux_phot-flux_sig_phot)),\
              np.log10(c/(wave_phot*1e-4)*(flux_phot+flux_sig_phot))-np.log10(c/(wave_phot*1e-4)*flux_phot)],\
        fmt='s',mfc='DimGray',mec='k',markersize=8)

    # perform the same procedure of flux extraction of aperture flux with observed spectra
    wl_aper = np.array(wl_aper)
    obs_aper_wl = wl_aper[(wl_aper >= min(wave_obs)) & (wl_aper <= max(wave_obs))]
    obs_aper_sed = np.empty_like(obs_aper_wl)
    sed_tot = c/(wl_tot*1e-4)*flux_tot
    # wl_tot and flux_tot are already hstacked and sorted by wavelength
    for i in range(0, len(obs_aper_wl)):
        if (obs_aper_wl[i] < 50.) & (obs_aper_wl[i] >= 5):
            res = 60.
        elif obs_aper_wl[i] < 5:
            res = 10.
        else:
            res = 1000.
        ind = np.where((wl_tot < obs_aper_wl[i]*(1+1./res)) & (wl_tot > obs_aper_wl[i]*(1-1./res)))
        if len(ind[0]) != 0:
            obs_aper_sed[i] = np.mean(sed_tot[ind])
        else:
            f = interp1d(wl_tot, sed_tot)
            obs_aper_sed[i] = f(wl_aper[i])
    aper_obs, = ax_sed.plot(np.log10(obs_aper_wl),np.log10(obs_aper_sed), 's-', mec='None', mfc='r', color='r',markersize=10, linewidth=1.5, label=r'$\rm{F_{aper,obs}}$')

    # iterate through three models, best-fit, Kristensen 2012, and Bourke 1997 geometry
    for mod in modelname:

        # read in the raw output from hyperion and then save the extraction results in text file, otherwise read the text files instead.
        if rtout == True:
            # Open the model
            m = ModelOutput(mod+'.rtout')

            sed_inf = m.get_sed(group=0, inclination=0, aperture=-1, distance=dstar * pc)
            flux_aper = np.empty_like(wl_aper)
            unc_aper = np.empty_like(wl_aper)
            for i in range(0, len(wl_aper)):
                sed_dum = m.get_sed(group=i+1, inclination=0, aperture=-1, distance=dstar * pc)
                # use a rectangle function the average the simulated SED
                # apply the spectral resolution
                if (wl_aper[i] < 50.) & (wl_aper[i] >= 5):
                    res = 60.
                elif wl_aper[i] < 5:
                    res = 10.
                else:
                    res = 1000.
                ind = np.where((sed_dum.wav < wl_aper[i]*(1+1./res)) & (sed_dum.wav > wl_aper[i]*(1-1./res)))
                if len(ind[0]) != 0:
                    flux_aper[i] = np.mean(sed_dum.val[ind])
                else:
                    f = interp1d(sed_dum.wav, sed_dum.val)
                    flux_aper[i] = f(wl_aper[i])
            
            # make the variables consistent with others
            sim_inf = sed_inf.wav
            sim_sed_inf = sed_inf.val

            # save the results in text files
            # unapertured SED
            foo = open(outdir+mod+'_sed_inf.txt','w')
            foo.write('%12s \t %12s \n' % ('wave','vSv'))
            for i in range(0, len(sed_inf.wav)):
                foo.write('%12g \t %12g \n' % (sed_inf.wav[i], sed_inf.val[i]))
            foo.close()
            # SED with convolution of aperture sizes
            foo = open(outdir+mod+'_sed_w_aperture.txt','w')
            foo.write('%12s \t %12s \n' % ('wave','vSv'))
            for i in range(0, len(wl_aper)):
                foo.write('%12g \t %12g \n' % (wl_aper[i], flux_aper[i]))
            foo.close()
        else:
            # read in the extracted text files
            (sim_inf, sim_sed_inf) = np.genfromtxt(mod+'_sed_inf.txt', skip_header=1).T
            (wl_aper, flux_aper) = np.genfromtxt(mod+'_sed_w_aperture.txt', skip_header=1).T

        aper, = ax_sed.plot(np.log10(wl_aper),np.log10(flux_aper),'o', linestyle=style[modelname.index(mod)], mec=color_list[modelname.index(mod)], mfc='None', color=color_list[modelname.index(mod)],markersize=12, markeredgewidth=3, linewidth=1.7, label=label[modelname.index(mod)])
        # calculate the bolometric luminosity of the aperture 
        l_bol_sim = l_bol(wl_aper, flux_aper/(c/np.array(wl_aper)*1e4)*1e23, dstar)
        print 'Bolometric luminosity of simulated spectrum: %5.2f lsun' % l_bol_sim

    # plot setting
    ax_sed.set_xlabel(r'$\rm{log\,\lambda\,({\mu}m)}$',fontsize=mag*20)
    ax_sed.set_ylabel(r'$\rm{log\,\nu S_{\nu}\,(erg\,cm^{-2}\,s^{-1})}$',fontsize=mag*20)
    [ax_sed.spines[axis].set_linewidth(1.5*mag) for axis in ['top','bottom','left','right']]
    ax_sed.minorticks_on()
    ax_sed.tick_params('both',labelsize=mag*18,width=1.5*mag,which='major',pad=15,length=5*mag)
    ax_sed.tick_params('both',labelsize=mag*18,width=1.5*mag,which='minor',pad=15,length=2.5*mag)

    ax_sed.set_ylim([-13,-7.5])
    ax_sed.set_xlim([0.3,3])

    lg_data = ax_sed.legend(loc='lower right',fontsize=12*mag,numpoints=1,framealpha=0.3)

    # Write out the plot
    fig.savefig(outdir+'three_models_vs_obs_sed.pdf',format='pdf',dpi=300,bbox_inches='tight')
    fig.clf()

import numpy as np
indir = '/Users/yaolun/bhr71/hyperion/controlled/'
outdir = '/Users/yaolun/Copy/Papers/yaolun/bhr71/figures/'
obs = '/Users/yaolun/bhr71/obs_for_radmc/'
# obs = None

# # grid of cs and age
# array = np.array([[7,8,9,10,11],[12,13,14,15,16],[17,18,19,20,21]])
# cslist = [0.1,0.2,0.3]
# agelist = [1e4,2.5e4,5e4,7.5e4,1e5]
# sed_grid_cs_age(indir, array, outdir, cslist, agelist, obs= None)

# # grid of Omega0
# array = np.array([24,23,22])
# sed_omega(indir, array, outdir, obs= None, compact=True)

# # # grid of disk parameters
# # disk mass
# array = np.array([71,72,73,74,75])
# xlabel = r'$M_{disk}\/[M_{\odot}]\/(0.01,\/0.03,\/0.05,\/0.07,\/0.1)$'
# compact = [r'$M_{disk}=0.01\/M_{\odot}$',r'$M_{disk}=0.03\/M_{\odot}$',r'$M_{disk}=0.05\/M_{\odot}$',r'$M_{disk}=0.07\/M_{\odot}$',r'$M_{disk}=0.1\/M_{\odot}$']
# plotname = 'disk_mdisk'
# sed_five(indir, array, outdir, xlabel, plotname, obs= None, zoom=True, compact=compact, yrange=[-13,-8])
# # flare power
# array = np.array([31,32,33,34,35])
# xlabel = r'$\beta\/(1.0,\/1.2,\/1.4,\/1.6,\/1.8)$'
# compact = [r'$\beta=1.0$',r'$\beta=1.2$',r'$\beta=1.4$',r'$\beta=1.6$',r'$\beta=1.8$']
# plotname = 'disk_beta'
# sed_five(indir, array, outdir, xlabel, plotname, obs= None, zoom=True, compact=compact, yrange=[-13,-8])
# # scale height
# array = np.array([36,37,38,39,40])
# xlabel = r'$h_{100}\/[AU]\/(6,\/8,\/10\/,12,\/14)$'
# compact = [r'$h_{100}=6\/AU$',r'$h_{100}=8\/AU$',r'$h_{100}=10\/AU$',r'$h_{100}=12\/AU$',r'$h_{100}=14\/AU$']
# plotname = 'disk_h100'
# sed_five(indir, array, outdir, xlabel, plotname, obs=None, zoom=True, compact=compact, yrange=[-13,-8])

# # grid of theta_cav and incl.
# array = np.array([[44,45,46,47,48],[49,50,51,52,53],[54,55,56,57,58]])
# # sed_grid_theta_cav_incl(indir, array, outdir, obs= None)
# sed_grid_theta_cav_incl(indir, array, outdir, obs= None, compact=True)
# # only for incl. = 50
# array = np.array([54,55,56,57,58])
# xlabel = r'$\theta_{cav}\/[deg.]\/(15^{\circ}, 20^{\circ}, 25^{\circ}, 30^{\circ}, 35^{\circ})$'
# plotname = 'theta_cav_incl50'
# compact = [r'$\theta_{cav}=15^{\circ}$',r'$\theta_{cav}=20^{\circ}$',r'$\theta_{cav}=25^{\circ}$',\
#            r'$\theta_{cav}=30^{\circ}$',r'$\theta_{cav}=35^{\circ}$']
# sed_five(indir, array, outdir, xlabel, plotname, obs=None, compact=compact)

# # grid of rho_cav_center and sed_rho_cav_edge
# array = np.array([[59,60,61,62],[63,64,65,66],[67,68,69,70]])
# # sed_grid_rho_cav_centeredge(indir, array, outdir, obs= None)
# sed_grid_rho_cav_centeredge(indir, array, outdir, obs= None, compact=True)

# # disk & no dis comparison
# disk & no disk
# array = np.array([71,76])
# disk_exist_com(indir, array, outdir, obs=obs)

# # grid of tstar
# array = np.array([69,70,71])
# sed_tstar(indir, array, outdir, obs=obs)

# # grid of rstar
# array = np.array([72,73,74])
# sed_rstar(indir, array, outdir, obs=obs)

# grid of R_env_max
# array = np.array([63,64,65,66,67])
# xlabel = r'$R_{env,max}\/[AU]\/(7.5\times 10^{3},\/1\times 10^{4},\/2.5\times 10^{4},\/5\times 10^{4},\/7.5\times 10^{4})}$'
# compact = [r'$R_{env,max}=7.5\times 10^{3}\/AU}$',r'$R_{env,max}=1.0\times 10^{4}\/AU}$',r'$R_{env,max}=2.5\times 10^{4}\/AU}$',r'$R_{env,max}=5.0\times 10^{4}\/AU}$',r'$R_{env,max}=7.5\times 10^{4}\/AU}$']
# array = np.array([4,5,6])
# xlabel = r'$R_{env,max}\/[AU]\/(7.5\times 10^{3},\/1\times 10^{4},\/2.5\times 10^{4})$'
# compact = [r'$R_{env,max}=7.5\times 10^{3}\/AU$',r'$R_{env,max}=1.0\times 10^{4}\/AU$',r'$R_{env,max}=2.5\times 10^{4}\/AU$']
# plotname = 'r_max'
# sed_five(indir, array, outdir, xlabel, plotname, obs= None, tbol=True, compact=compact)

# # grid of continuous cavity power law
# array = np.array([13,14,15])
# sed_cav_powerlaw(indir, array, outdir, obs=obs)

# grid of cavity structure comparison
# array = np.array([42,43,71,41])  # r-2, r-1.5, const.+r-2, uniform
# sed_cav_struc_com(indir, array, outdir, obs=obs)

# grid of tstar with the same lstar
# array = np.array([1,2,3])
# sed_lum(indir, array, outdir)

# model_vs_obs('model46', '/Users/yaolun/bhr71/hyperion/cycle7/', '/Users/yaolun/test/', obs=obs)

three_model_vs_obs(['/Users/yaolun/bhr71/hyperion/cycle5/model16','/Users/yaolun/bhr71/hyperion/cycle7/model53','/Users/yaolun/bhr71/hyperion/cycle7/model46'],\
    '/Users/yaolun/bhr71/hyperion/cycle7/', '/Users/yaolun/test/',\
    [r'$\rm{Kristensen\,et.\,al.\,2012}$', r'$\rm{geometry\,from\,Bourke\,et.\,al.\,1997}$',r'$\rm{This\,Study}$'], obs)