def fir_chi2_2d(array_list, keywords, obs, wl_aper=None, fixed=False, ref=None, spitzer_only=False, \
    herschel_only=False, plot_model=True, zoom_1d=None, shade=True):
    """
    array_list: contains dictionaries, each dictionary represents a location of 'model_list.txt', and the model numbers within.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import ascii
    import sys
    import os
    sys.path.append(os.path.expanduser('~')+'/programs/spectra_analysis/')
    from phot_filter import phot_filter
    from get_obs import get_obs
    import astropy.constants as const
    from scipy.interpolate import interp1d
    from scipy.interpolate import griddata
    from hyperion.model import ModelOutput
    from matplotlib.colors import LogNorm
    import copy
    import collections
    import seaborn.apionly as sns
    import matplotlib as mpl

    # function for checking duplicate in lists
    compare = lambda x, y: collections.Counter(x) == collections.Counter(y)

    # constant setup
    c = const.c.cgs.value
    LS = const.L_sun.cgs.value
    pc = const.pc.cgs.value
    dstar = 200.0

    # chi2 function
    def fir_chi2(obs, sim, wave=[70.,100.,160.,250.,350.,500.], log=False, robust_est=False):

        # experimental procedure
        # introduce a systematic error estimated from the gap between PACS and SPIRE.
        # Take the ratio of this gap to the peak value as the fraction in fractional uncertainty.
        # Then add two uncertainties in quadrature

        chi2 = 0
        if log == False:
            for w in wave:
                val = (sim['sed'][sim['wave'] == w] - obs['sed'][obs['wave'] == w]) / obs['sed'][obs['wave'] == w]
                unc_2 = (sim['sed'][sim['wave'] == w]/obs['sed'][obs['wave'] == w])**2 *\
                        ( (sim['sigma'][sim['wave'] == w]/sim['sed'][sim['wave'] == w])**2 + (obs['sigma'][obs['wave'] == w]/obs['sed'][obs['wave'] == w])**2 )
                # unc_2 = 1
                # print w, val, unc_2**0.5
                chi2 = chi2 + val**2 / unc_2
        else:
            # not proper functioning at this moment
            for w in wave:
                # print w, (sim['sed'][sim['wave'] == w]-obs['sed'][obs['wave'] == w])**2
                # print 'print'
                # print np.log10(obs['sed'][obs['wave'] == w])-np.log10(obs['sed'][obs['wave'] == w]-obs['sigma'][obs['wave'] == w])
                # print np.log10(obs['sed'][obs['wave'] == w]+obs['sigma'][obs['wave'] == w])-np.log10(obs['sed'][obs['wave'] == w])
                chi2 = chi2 + ((np.log10(sim['sed'][sim['wave'] == w])-np.log10(obs['sed'][obs['wave'] == w]))**2) /\
                         (max(np.log10(obs['sed'][obs['wave'] == w])-np.log10(obs['sed'][obs['wave'] == w]-obs['sigma'][obs['wave'] == w]), np.log10(obs['sed'][obs['wave'] == w]+obs['sigma'][obs['wave'] == w])-np.log10(obs['sed'][obs['wave'] == w])))**2
        return chi2, len(wave)

    # setup the aperture size
    if wl_aper == None:
        # wl_aper = [3.6, 4.5, 5.8, 8.0, 8.5, 9, 9.7, 10, 10.5, 11, 16, 20, 24, 35, 70, 100, 160, 250, 350, 500, 850]
        # wl_aper = [3.6, 4.5, 8.5, 9, 9.7, 10, 16, 20, 24, 35, 70, 100, 160, 250, 350, 500]
        wl_aper = [3.6, 4.5, 8.5, 9, 9.7, 10, 16, 20, 24, 30, 70, 100, 160, 250, 350, 500]

        # wl_aper = [70., 100., 160., 250., 350., 500.]
    if spitzer_only:
        wl_aper = [5.8, 8.0, 8.5, 9, 9.7, 10, 10.5, 11, 16, 20, 24, 35]
    if herschel_only:
        wl_aper = [70, 100, 160, 250, 350, 500, 1300]
        # wl_aper = [35, 70, 85, 100, 120, 140, 160, 200, 250, 300, 350, 400, 500, 600]
    # test version:
    # the current wavelength channals: [35, 70, 85, 100, 120, 140, 160, 200, 250, 300, 350, 400, 500, 600, 850]
    wl_aper = np.array(wl_aper, dtype=float)

    # read the observed SED and extract with apertures
    bhr71 = get_obs(obs)
    wave_obs, flux_obs, sed_obs_noise = bhr71['spec']
    # add IRAC1, IRAC2, and SEST 1.3 mm photometry data
    wave_obs = np.hstack((wave_obs, bhr71['phot'][0][0:2], bhr71['phot'][0][7]))
    flux_obs = np.hstack((flux_obs, bhr71['phot'][1][0:2], bhr71['phot'][1][7]))
    sed_obs_noise = np.hstack((sed_obs_noise, bhr71['phot'][2][0:2], bhr71['phot'][2][7]))

    sed_obs = c/(wave_obs*1e-4)*flux_obs*1e-23
    sed_obs_noise = c/(wave_obs*1e-4)*sed_obs_noise*1e-23
    # print sed_obs_noise
    obs_aper_sed = np.empty_like(wl_aper)
    obs_aper_sed_noise = np.empty_like(wl_aper)

    # sort the wavelength
    sed_obs = sed_obs[np.argsort(wave_obs)]
    sed_obs_noise = sed_obs_noise[np.argsort(wave_obs)]
    wave_obs = wave_obs[np.argsort(wave_obs)]

    # calculate the spectrophotometry
    for i in range(0, len(wl_aper)):
        # apply the filter function
        # decide the filter name
        if wl_aper[i] == 70:
            fil_name = 'Herschel PACS 70um'
        elif wl_aper[i] == 100:
            fil_name = 'Herschel PACS 100um'
        elif wl_aper[i] == 160:
            fil_name = 'Herschel PACS 160um'
        elif wl_aper[i] == 250:
            fil_name = 'Herschel SPIRE 250um'
        elif wl_aper[i] == 350:
            fil_name = 'Herschel SPIRE 350um'
        elif wl_aper[i] == 500:
            fil_name = 'Herschel SPIRE 500um'
        elif wl_aper[i] == 3.6:
            fil_name = 'IRAC Channel 1'
        elif wl_aper[i] == 4.5:
            fil_name = 'IRAC Channel 2'
        elif wl_aper[i] == 5.8:
            fil_name = 'IRAC Channel 3'
        elif wl_aper[i] == 8.0:
            fil_name = 'IRAC Channel 4'
        elif wl_aper[i] == 24:
            fil_name = 'MIPS 24um'
        # elif wl_aper[i] == 850:
        #     fil_name = 'SCUBA 850WB'
        # do not have SCUBA spectra
        else:
            fil_name = None

        if fil_name != None:
            filter_func = phot_filter(fil_name, '/Users/yaolun/programs/misc/hyperion/')
            # Observed SED needs to be trimmed before applying photometry filters
            filter_func = filter_func[(filter_func['wave']/1e4 >= min(wave_obs))*\
                                      ((filter_func['wave']/1e4 >= 54.8)+(filter_func['wave']/1e4 <= 36.0853))*\
                                      ((filter_func['wave']/1e4 <= 95.05)+(filter_func['wave']/1e4 >=103))*\
                                      ((filter_func['wave']/1e4 <= 190.31)+(filter_func['wave']/1e4 >= 195))*\
                                      (filter_func['wave']/1e4 <= max(wave_obs))]
            f = interp1d(wave_obs, sed_obs)
            f_unc = interp1d(wave_obs, sed_obs_noise)
            obs_aper_sed[i] = np.trapz(filter_func['wave']/1e4, f(filter_func['wave']/1e4)*filter_func['transmission'])/np.trapz(filter_func['wave']/1e4, filter_func['transmission'])
            obs_aper_sed_noise[i] = abs(np.trapz((filter_func['wave']/1e4)**2, (f_unc(filter_func['wave']/1e4)*filter_func['transmission'])**2))**0.5 / abs(np.trapz(filter_func['wave']/1e4, filter_func['transmission']))
        else:
            # use a rectangle function the average the simulated SED
            # apply the spectral resolution
            if (wl_aper[i] < 50.) & (wl_aper[i] >= 5):
                res = 60.
            elif wl_aper[i] < 5:
                res = 10.
            else:
                res = 1000.
            ind = np.where((wave_obs < wl_aper[i]*(1+1./res)) & (wave_obs > wl_aper[i]*(1-1./res)))
            if len(ind[0]) != 0:
                obs_aper_sed[i] = np.mean(sed_obs[ind])
                obs_aper_sed_noise[i] = np.mean(sed_obs_noise[ind])
            else:
                f = interp1d(wave_obs, sed_obs)
                f_unc = interp1d(wave_obs, sed_obs_noise)
                obs_aper_sed[i] = f(wl_aper[i])
                obs_aper_sed_noise[i] = f_unc(wl_aper[i])
    if 1300. in wl_aper:
        np.append(obs_aper_sed, 3.8)
        np.append(obs_aper_sed_noise, 0.57)
    # calculate Chi2 from simulated SED
    p1 = []
    p2 = []
    model_label = []
    chi2 = []
    total_chi2 = []

    # plot the simulation on top of the observation
    fig = plt.figure(figsize=(8,6))
    ax_sim = fig.add_subplot(111)

    num_file = len(array_list)
    for ifile in range(0,num_file):
        listpath = array_list[ifile]['listpath']
        datapath = array_list[ifile]['datapath']
        model_num = array_list[ifile]['model_num']

        model_list = ascii.read(listpath)

        if ref != None:
            ignore_col = ['d_sub', 'M_env_dot', 'R_inf', 'R_cen', 'mstar', 'M_tot_gas', 'M_disk']
            ref_params = copy.copy(model_list)
            if fixed == True:
                ref_params.remove_column(keywords['col'][0])
            else:
                ref_params.remove_columns(keywords['col'])
            ref_params.remove_columns(ignore_col)
            ref_params.remove_column('Model#')
            ref_params = (ref_params[:][model_list['Model#'] == 'Model'+str(ref)])[0].data

        yerr = np.hstack(( obs_aper_sed_noise[wl_aper < 50.],\
              ( obs_aper_sed_noise[(wl_aper >= 50.) & (wl_aper < 70.)]**2 + (obs_aper_sed[(wl_aper >= 50.) & (wl_aper < 70)]*0.11)**2 )**0.5,\
              ( obs_aper_sed_noise[(wl_aper >= 70.) & (wl_aper < 200.)]**2 + (obs_aper_sed[(wl_aper >= 70.) & (wl_aper < 200)]*0.12)**2 )**0.5,\
              ( obs_aper_sed_noise[wl_aper >= 200.]**2 + (obs_aper_sed[wl_aper >= 200.]*0.07)**2 )**0.5))

        ax_sim.errorbar(wl_aper, obs_aper_sed, yerr=yerr, fmt='s', color='DimGray')
        ax_sim.plot(bhr71['spec'][0], c/(bhr71['spec'][0]*1e-4)*bhr71['spec'][1]*1e-23, '-', color='DimGray')

        for i in range(0, len(model_num)):
            imod = model_num[i]
            # if (model_list[keywords['col'][0]][model_list['Model#'] == 'Model'+str(imod)]).data not in [100., 1000., 10000., 15000.]:
            #     continue
            model_dum = ascii.read(datapath+'/model'+str(imod)+'_sed_w_aperture.txt')
            obs = {'wave': np.array(wl_aper), 'sed': obs_aper_sed, 'sigma': yerr}
            sim = {'wave': model_dum['wave'].data, 'sed': model_dum['vSv'].data, 'sigma': model_dum['sigma_vSv'].data}
            chi2_dum, n = fir_chi2(obs, sim, wave=wl_aper)

            # make a test plot
            # ax_sim.errorbar(sim['wave'], sim['sed'], yerr=sim['sigma'], fmt='o', mec='None',alpha=0.5,\
                # label=r'$\rm{'+str(int((model_list[keywords['col'][0]][model_list['Model#'] == 'Model'+str(imod)]).data))+'\,yr}$')

            reduced_chi2_dum = chi2_dum/(n-2-1)
            total_chi2.extend(reduced_chi2_dum)

            if ref == None:
                # read the parameter values
                p1.extend((model_list[keywords['col'][0]][model_list['Model#'] == 'Model'+str(imod)]).data)
                if fixed == False:
                    p2.extend((model_list[keywords['col'][1]][model_list['Model#'] == 'Model'+str(imod)]).data)
                model_label.append(str(imod))
                # print reduced_chi2_dum
                chi2.extend(reduced_chi2_dum)
            else:
                # get other parameters of model i
                dum_params = copy.copy(model_list)
                # clean up test parameters and irrelvant parameters
                if fixed == True:
                    dum_params.remove_column(keywords['col'][0])
                else:
                    dum_params.remove_columns(keywords['col'])
                dum_params.remove_columns(ignore_col)
                dum_params.remove_column('Model#')
                dum_params = (dum_params[:][model_list['Model#'] == 'Model'+str(imod)])[0].data
                # compare the controlled parameters with the reference
                if compare(ref_params, dum_params) == False:
                    continue
                else:
                    # print (model_list[keywords['col'][0]][model_list['Model#'] == 'Model'+str(imod)]).data, imod
                    p1.extend((model_list[keywords['col'][0]][model_list['Model#'] == 'Model'+str(imod)]).data)
                    if imod == ref:
                        ref_p1 = float((model_list[keywords['col'][0]][model_list['Model#'] == 'Model'+str(imod)]).data)
                    if fixed == False:
                        p2.extend((model_list[keywords['col'][1]][model_list['Model#'] == 'Model'+str(imod)]).data)
                        if imod == ref:
                            ref_p2 = float((model_list[keywords['col'][1]][model_list['Model#'] == 'Model'+str(imod)]).data)
                    model_label.append(str(imod))

                    # print reduced_chi2_dum
                    chi2.extend(reduced_chi2_dum)


    # finish the plot
    ax_sim.legend(loc='best', numpoints=1)
    ax_sim.set_xlim([30,700])
    fig.savefig('/Users/yaolun/test/chi2_simulation_summary.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()

    # convert to array
    p1 = np.array(p1)
    p2 = np.array(p2)
    chi2 = np.array(chi2)

    # if dealing with age convert to better unit for plotting
    for i in range(len(keywords['col'])):
        if keywords['col'][i] == 'age':
            if i == 0:
                p1 = p1/1e4
            elif (i ==1) and (fixed == False):
                p2 = p2/1e4

    # plot 1d relation
    if fixed == True:
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)

        p1 = np.squeeze(p1)
        chi2 = np.squeeze(chi2)

        p1 = np.array(p1); chi2 = np.array(chi2)
        ax.plot(p1[np.argsort(p1)], chi2[np.argsort(p1)], 'o-', mec='None', color='Green', linewidth=1, markersize=4)
        ax.set_xlabel(keywords['label'][0], fontsize=18)
        ax.set_ylabel(r'$\rm{\chi^{2}_{reduced}}$', fontsize=18)

        ax.axvline(2.3059, color='k', linestyle='--', linewidth=1)

        if shade:
            if ref != None:
                # mark the region where the chi-squared ranging from lowest possible value to double
                ax.axhspan(chi2[p1*1e4 == ref_p1], 2*chi2[p1*1e4 == ref_p1], color='grey', alpha=0.3)
                # make vertical lines to show the uncertainty with a certain chi-square criteria
                ax.axvspan(min(p1[chi2 <= min(chi2)*2]), max(p1[chi2 <= min(chi2)*2]),
                           color='b', alpha=0.3)

        # ax.set_yscale('log')

        if zoom_1d != None:
            ax.set_xlim(zoom_1d)
        # else:
        #     # fig.gca().set_xlim(left=0)
        #     ax.set_xlim([0,10])
        # ax.set_ylim([0, 50])

        [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
        ax.minorticks_on()
        ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
        ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)

        fig.savefig('/Users/yaolun/test/chi2_'+str(keywords['col'][0])+'_1d.pdf', format='pdf', dpi=300, bbox_inches='tight')
        fig.clf()

    else:

        # write out the chi2
        foo = open('/Users/yaolun/bhr71/hyperion/chi2_grid/chi2-2d.txt', 'w')
        foo.write('{} \t {} \t {} \n'.format('chisq', keywords['col'][0], keywords['col'][1]))
        for i in range(len(chi2)):
            foo.write('{} \t {} \t {} \n'.format(chi2[i], p1[i], p2[i]))
        foo.close()

        # rebin the data and plot 2D contour
        from scipy.interpolate import griddata
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        # need to normalize the array to [0,1]
        p1_norm = (p1-min(p1))/(max(p1)-min(p1))
        p2_norm = (p2-min(p2))/(max(p2)-min(p2))

        # special treatment for rho_cav_center parameter
        if keywords['col'][0] == 'rho_cav_center':
            p1_norm = (np.log10(p1)-min(np.log10(p1)))/(max(np.log10(p1))-min(np.log10(p1)))
            p1_label = []
            for ip1 in p1_norm:
                if ip1 not in p1_label:
                    p1_label.append(ip1)
        if keywords['col'][1] == 'rho_cav_center':
            p2_norm = (np.log10(p2)-min(np.log10(p2)))/(max(np.log10(p2))-min(np.log10(p2)))
            p2_label = []
            for ip2 in p2_norm:
                if ip2 not in p2_label:
                    p2_label.append(ip2)

        x = np.linspace(0, 1, 50)
        y = np.linspace(0, 1, 50)

        z = griddata((p1_norm, p2_norm), chi2, (x[None,:], y[:,None]), method='nearest')
        if z[np.isnan(z) != True].min() < 0:
            print 'Minimum of z is below zero.  Change the interpolation method to linear'
            z = griddata((p1_norm, p2_norm), chi2, (x[None,:], y[:,None]), method='linear')
        masked_array = np.ma.array (z, mask=np.isnan(z))

        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111)

        # plot the contour with color and lines
        # ax.contour(x, y, z, 10, linewidths=0.5,colors='k', norm=LogNorm(vmin=chi2.min(), vmax=chi2.max()))
        cmap = plt.cm.CMRmap
        # import custom colormap
        from custom_colormap import custom_colormap
        cmap = mpl.colors.ListedColormap(custom_colormap())
        cmap.set_bad('Red',1.)
        im = ax.imshow(masked_array, cmap=cmap, origin='lower', extent=[0,1,0,1],\
            norm=LogNorm(vmin=chi2.min(), vmax=chi2.max()))  # chi2.max()  # 1e7
        # Blues_r
        ax.set_xticks(np.linspace(0, 1, 5))
        ax.set_xticklabels(np.linspace(min(p1), max(p1), 5))
        ax.set_yticks(np.linspace(0, 1, 5))
        #
        # p2 = 0.3 *pc/10e5 / np.tan(np.radians(p2)) / (3600.*24*365) / 1e4
        #
        ax.set_yticklabels(np.linspace(min(p2), max(p2), 5))

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = fig.colorbar(im, cax=cax)
        cb.solids.set_edgecolor("face")
        cb.ax.minorticks_on()
        cb.ax.set_ylabel(r'$\rm{\chi^{2}_{reduce}}$', fontsize=20)
        cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
        plt.setp(cb_obj,fontsize=12)
        if plot_model:
            # plot the original data points
            ori_data = ax.scatter(p1_norm,p2_norm, marker='o',c='b',s=5)

        # print the model number near the points
        for i in range(len(model_label)):
            ax.annotate(model_label[i], (p1_norm[i], p2_norm[i]))

        ax.set_xlabel(keywords['label'][0], fontsize=20)
        ax.set_ylabel(keywords['label'][1], fontsize=20)
        [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
        ax.minorticks_on()
        ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
        ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)

        # special treatment fo rho_cav_center
        if keywords['col'][0] == 'rho_cav_center':
            ax.set_xticks(p1_label)
            ax.set_xticklabels([r'$\rm{5\times10^{-20}}$',r'$\rm{10^{-19}}$',r'$\rm{5\times10^{-19}}$',r'$\rm{10^{-18}}$',r'$\rm{5\times10^{-18}}$'])
            ax.tick_params(axis='x', which='minor', bottom='off', top='off')
        if keywords['col'][1] == 'rho_cav_center':
            ax.set_yticks(p2_label)
            ax.set_yticklabels([r'$\rm{5\times10^{-20}}$',r'$\rm{10^{-19}}$',r'$\rm{5\times10^{-19}}$',r'$\rm{10^{-18}}$',r'$\rm{5\times10^{-18}}$'])
            ax.tick_params(axis='y', which='minor', left='off', right='off')

        fig.savefig('/Users/yaolun/test/chi2_%s_%s.pdf' % (keywords['col'][0], keywords['col'][1]), format='pdf', dpi=300, bbox_inches='tight')

    return p1, p2, chi2

import numpy as np
# array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/cycle5/model_list.txt',
#                'datapath': '/Users/yaolun/bhr71/hyperion/cycle5',
#                'model_num': np.arange(38,46)},
#               {'listpath': '/Users/yaolun/bhr71/hyperion/cycle6/model_list.txt',
#                'datapath': '/Users/yaolun/bhr71/hyperion/cycle6',
#                'model_num': np.arange(39,49)}]
array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/chi2_grid/model_list.txt',
               'datapath': '/Users/yaolun/bhr71/hyperion/chi2_grid',
               'model_num': np.arange(1,26)}]

# array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/controlled/model_list.txt',
#                'datapath': '/Users/yaolun/bhr71/hyperion/controlled',
#                'model_num': np.arange(1,77)}]
# keywords = {'col':['age','view_angle'], 'label': [r'$\rm{age\,[10^{4}\,yr]}$', r'$\rm{viewing\,angle}$']}
# keywords = {'col':['age','theta_cav'], 'label': [r'$\rm{age\,[10^{4}\,yr]}$', r'$\rm{\theta_{cav}}$']}
# keywords = {'col':['age','Cs'], 'label': [r'$\rm{age\,[10^{4}\,yr]}$', r'$\rm{c_{s}\,[km\,s^{-1}]}$']}
# keywords = {'col':['view_angle','theta_cav'], 'label': [r'$\rm{viewing\,angle}$', r'$\rm{\theta_{cav}}$']}
# keywords = {'col':['age','rho_cav_edge'], 'label': [r'$\rm{age\,[10^{4}\,yr]}$', r'$\rm{R_{cav,\circ}\,[AU]}$']}
# keywords = {'col':['age','rho_cav_center'], 'label': [r'$\rm{age\,[10^{4}\,yr]}$', r'$\rm{\rho_{cav,\circ}\,[g\,cm^{-3}]}$']}
# keywords = {'col':['rho_cav_center','rho_cav_edge'], 'label': [r'$\rm{\rho_{cav,\circ}\,[g\,cm^{-3}]}$', r'$\rm{R_{cav,\circ}\,[AU]}$']}
# keywords_list = [{'col':['age','theta_cav'], 'label': [r'$\rm{age\,[10^{4}\,yr]}$', r'$\rm{\theta_{cav}\,[deg.]}$']},\
#                  {'col':['age','view_angle'], 'label': [r'$\rm{age\,[10^{4}\,yr]}$', r'$\rm{\theta_{incl}\,[deg.]}$']},\
#                  {'col':['view_angle','theta_cav'], 'label': [r'$\rm{\theta_{incl}\,[deg.]}$', r'$\rm{\theta_{cav}\,[deg.]}$']}]
keywords_list = [{'col':['age','view_angle'], 'label': [r'$\rm{age\,[10^{4}\,yr]}$', r'$\rm{\theta_{incl}\,[deg.]}$']}]
# keywords_list = [{'col':['theta_cav','view_angle'], 'label': [r'$\rm{\theta_{cav}\,[deg.]}$', r'$\rm{\theta_{incl}\,[deg.]}$']}]
obs = '/Users/yaolun/bhr71/best_calibrated/'

# for tstar & age
# keywords_list = [{'col':['tstar','age'], 'label': [r'$\rm{T_{\star}\,[K]}$', r'$\rm{age\,[10^{4}\,yr]}$']}]
#
for keywords in keywords_list:
    p1, p2, chi2 = fir_chi2_2d(array_list, keywords, obs, ref=19, plot_model=False, herschel_only=False)
    for i in range(len(p1)):
        print p1[i], p2[i], chi2[i]
    fir_chi2_2d(array_list, keywords, obs)

# # 1-D rho_cav_center
# array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/cycle9/model_list.txt',
#                'datapath': '/Users/yaolun/bhr71/hyperion/cycle9',
#                'model_num': np.arange(34,43)}]
# keywords = {'col':['rho_cav_center'], 'label': [r'$\rm{\rho_{cav}\,[g\,cm^{-3}]}$']}
# fir_chi2_2d(array_list, keywords, obs, fixed=True, ref=34)

# # 1-D inclination
# array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/cycle9/model_list.txt',
#                'datapath': '/Users/yaolun/bhr71/hyperion/cycle9',
#                'model_num': np.hstack((np.arange(44,50), 34))}]
# keywords = {'col':['view_angle'], 'label': [r'$\rm{\theta_{incl.}\,[deg.]}$']}
# fir_chi2_2d(array_list, keywords, obs, fixed=True, ref=34)

# # 1-D rho_cav_edge
# array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/cycle9/model_list.txt',
#                'datapath': '/Users/yaolun/bhr71/hyperion/cycle9',
#                'model_num': np.hstack((np.arange(50,54), 34, np.arange(68,71)))}]
# keywords = {'col':['rho_cav_edge'], 'label': [r'$\rm{R_{cav,\circ}\,[AU]}$']}
# fir_chi2_2d(array_list, keywords, obs, fixed=True, ref=34)

# 1-D age
# array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/cycle9/model_list.txt',
#                'datapath': '/Users/yaolun/bhr71/hyperion/cycle9',
#                # 'model_num': np.hstack((34,71))}]
#                'model_num': np.hstack((np.arange(54,68), 34, 71))}]

# array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/controlled/model_list.txt',
#                'datapath': '/Users/yaolun/bhr71/hyperion/controlled',
#                # 'model_num': np.hstack((np.arange(2,8),np.arange(14,47)))}]
#                'model_num': np.arange(86,120)}]
# array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/controlled/model_list.txt',
#                'datapath': '/Users/yaolun/bhr71/hyperion/controlled',
#                # 'model_num': np.hstack((np.arange(2,8),np.arange(14,47)))}]
#                'model_num': np.arange(155,223)}]
# keywords = {'col':['age'], 'label': [r'$\rm{t_{col}\,[10^{4}\,year]}$']}
# fir_chi2_2d(array_list, keywords, obs, fixed=True, ref=181, herschel_only=True, zoom_1d=[0,5])

# For estimating the effect of the uncertainty in dust opacity at submm wavelength.
# array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/test/model_list.txt',
#                'datapath': '/Users/yaolun/bhr71/hyperion/test/',
#                'model_num':np.array([17,18])}]
# keywords = {'col':['age'], 'label': [r'$\rm{t_{col}\,[10^{4}\,year]}$']}
# fir_chi2_2d(array_list, keywords, obs, fixed=True, herschel_only=True, zoom_1d=[0,5])

# For fitting the best age for p25 dust opactity
array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/controlled/model_list.txt',
               'datapath': '/Users/yaolun/bhr71/hyperion/controlled/',
               'model_num': np.arange(268,303)}]
            #    'model_num':np.hstack((17,np.arange(19,34)))}]
keywords = {'col':['age'], 'label': [r'$\rm{t_{col}\,[10^{4}\,year]}$']}
# keywords = {'col':['view_angle'], 'label':[r'$\rm{\theta_{incl.}\,[deg.]}$']}
# keywords = {'col':['tstar'], 'label': [r'$\rm{T_{\star}\,[K]}$']}
fir_chi2_2d(array_list, keywords, obs, fixed=True, herschel_only=True, ref=285, shade=False)
