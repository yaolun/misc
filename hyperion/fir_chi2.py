def fir_chi2_2d(array_list, keywords, obs, wl_aper=None, fixed=False, ref=None, spitzer_only=False, herschel_only=False):
    """
    array_list: contains dictionaries, each dictionary represents a location of 'model_list.txt', and the model numbers within.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import ascii
    from get_bhr71_obs import get_bhr71_obs
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
    dstar = 178.0

    # chi2 function
    def fir_chi2(obs, sim, wave=[70.,100.,160.,250.,350.,500.], log=False):

        # experimental procedure
        # introduce a systematic error estimated from the gap between PACS and SPIRE.
        # Take the ratio of this gap to the peak value as the fraction in fractional uncertainty.
        # Then add two uncertainties in quadrature

        chi2 = 0
        if log == False:
            for w in wave:
                # print w, (sim['sed'][sim['wave'] == w]-obs['sed'][obs['wave'] == w])**2
                val = (sim['sed'][sim['wave'] == w] - obs['sed'][obs['wave'] == w]) / obs['sed'][obs['wave'] == w]
                unc_2 = (sim['sed'][sim['wave'] == w]/obs['sed'][obs['wave'] == w])**2 *\
                        ( (sim['sigma'][sim['wave'] == w]/sim['sed'][sim['wave'] == w])**2 + (obs['sigma'][obs['wave'] == w]/obs['sed'][obs['wave'] == w])**2 ) + \
                        2 * (obs['sigma'][obs['wave'] == w]/obs['sed'][obs['wave'] == w])**2
                # unc = unc_2**0.5
                # print val**2, unc_2
                chi2 = chi2 + val**2 / unc_2
        else:
            # not proper functioning at this moment
            for w in wave:
                # print w, (sim['sed'][sim['wave'] == w]-obs['sed'][obs['wave'] == w])**2
                print 'print'
                print np.log10(obs['sed'][obs['wave'] == w])-np.log10(obs['sed'][obs['wave'] == w]-obs['sigma'][obs['wave'] == w])
                print np.log10(obs['sed'][obs['wave'] == w]+obs['sigma'][obs['wave'] == w])-np.log10(obs['sed'][obs['wave'] == w])
                chi2 = chi2 + ((np.log10(sim['sed'][sim['wave'] == w])-np.log10(obs['sed'][obs['wave'] == w]))**2) /\
                         (max(np.log10(obs['sed'][obs['wave'] == w])-np.log10(obs['sed'][obs['wave'] == w]-obs['sigma'][obs['wave'] == w]), np.log10(obs['sed'][obs['wave'] == w]+obs['sigma'][obs['wave'] == w])-np.log10(obs['sed'][obs['wave'] == w])))**2
        return chi2, len(wave)

    # setup the aperture size
    if wl_aper == None:
        # wl_aper = [3.6, 4.5, 5.8, 8.0, 8.5, 9, 9.7, 10, 10.5, 11, 16, 20, 24, 35, 70, 100, 160, 250, 350, 500, 850]
        wl_aper = [5.8, 8.0, 8.5, 9, 9.7, 10, 10.5, 11, 16, 20, 24, 35, 70, 100, 160, 250, 350, 500]
        # wl_aper = [70., 100., 160., 250., 350., 500.]
    if spitzer_only:
        wl_aper = [5.8, 8.0, 8.5, 9, 9.7, 10, 10.5, 11, 16, 20, 24, 35]
    if herschel_only:
        wl_aper = [70., 100., 160., 250., 350., 500.]

    # read the observed SED and extract with apertures
    bhr71 = get_bhr71_obs(obs)
    wave_obs, flux_obs, sed_obs_noise = bhr71['spec']
    # add IRAC1 and IRAC2 photometry data
    wave_obs = np.hstack((wave_obs, bhr71['phot'][0][0:2]))
    flux_obs = np.hstack((flux_obs, bhr71['phot'][1][0:2]))
    sed_obs_noise = np.hstack((sed_obs_noise, bhr71['phot'][2][0:2]))
    # wave_obs = wave_obs[wave_obs > 50]
    # flux_obs = flux_obs[wave_obs > 50]
    sed_obs = c/(wave_obs*1e-4)*flux_obs*1e-23
    sed_obs_noise = c/(wave_obs*1e-4)*sed_obs_noise*1e-23
    # print sed_obs_noise
    obs_aper_sed = np.empty_like(wl_aper)
    obs_aper_sed_noise = np.empty_like(wl_aper)

    # sort the wavelength
    sed_obs = sed_obs[np.argsort(wave_obs)]
    sed_obs_noise = sed_obs_noise[np.argsort(wave_obs)]
    wave_obs = wave_obs[np.argsort(wave_obs)]

    # setup resolution
    for i in range(0, len(wl_aper)):
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
            f_noise = interp1d(wave_obs, sed_obs_noise)
            obs_aper_sed[i] = f(wl_aper[i])
            obs_aper_sed_noise[i] = f_noise(wl_aper[i])

    # calculate Chi2 from simulated SED
    p1 = []
    p2 = []
    model_label = []
    chi2 = []
    total_chi2 = []

    num_file = len(array_list)
    for ifile in range(0,num_file):
        listpath = array_list[ifile]['listpath']
        datapath = array_list[ifile]['datapath']
        model_num = array_list[ifile]['model_num']

        model_list = ascii.read(listpath)

        if ref != None:
            ignore_col = ['d_sub', 'M_env_dot', 'R_inf', 'R_cen', 'mstar', 'total_mass', 'M_disk']
            ref_params = copy.copy(model_list)
            if fixed == True:
                ref_params.remove_column(keywords['col'][0])
            else:
                ref_params.remove_columns(keywords['col'])
            ref_params.remove_columns(ignore_col)
            ref_params.remove_column('Model#')
            ref_params = (ref_params[:][model_list['Model#'] == 'Model'+str(ref)])[0].data

        # find the model has minimun chi2 first
        # for i in range(0, len(model_num)):
        #     imod = model_num[i]
        #     model_dum = ascii.read(datapath+'/model'+str(imod)+'_sed_w_aperture.txt')
        #     chi2_dum, n = fir_chi2({'wave': np.array(wl_aper), 'sed': obs_aper_sed, 'sigma': sed_obs_noise}, 
        #         {'wave': model_dum['wave'].data, 'sed': model_dum['vSv'].data, 'sigma': model_dum['sigma_vSv'].data}, wave=wl_aper)
        #     reduced_chi2_dum = chi2_dum/(n-2-1)
        #     total_chi2.extend(reduced_chi2_dum)
        # print total_chi2
        # print min(total_chi2), np.where(total_chi2 == min(total_chi2))

        for i in range(0, len(model_num)):
            imod = model_num[i]
            model_dum = ascii.read(datapath+'/model'+str(imod)+'_sed_w_aperture.txt')
            chi2_dum, n = fir_chi2({'wave': np.array(wl_aper), 'sed': obs_aper_sed, 'sigma': sed_obs_noise}, 
                {'wave': model_dum['wave'].data, 'sed': model_dum['vSv'].data, 'sigma': model_dum['sigma_vSv'].data}, wave=wl_aper)
            reduced_chi2_dum = chi2_dum/(n-2-1)
            total_chi2.extend(reduced_chi2_dum)
            # manually exclude much older age
            # if keywords['col'][0] == 'age':
            #     if (model_list[keywords['col'][0]][model_list['Model#'] == 'Model'+str(imod)]).data >= 5e5:
            #         continue
            # if keywords['col'][1] == 'age':
            #     if (model_list[keywords['col'][1]][model_list['Model#'] == 'Model'+str(imod)]).data >= 5e5:
            #         continue                    
            if ref == None:
                # read the parameter values
                p1.extend((model_list[keywords['col'][0]][model_list['Model#'] == 'Model'+str(imod)]).data)
                if fixed == False:
                    p2.extend((model_list[keywords['col'][1]][model_list['Model#'] == 'Model'+str(imod)]).data)
                model_label.append(str(imod))
            else:
                # get other parameters of model i
                dum_params = copy.copy(model_list)
                if fixed == True:
                    dum_params.remove_column(keywords['col'][0])
                else:
                    dum_params.remove_columns(keywords['col'])
                dum_params.remove_columns(ignore_col)
                dum_params.remove_column('Model#')
                dum_params = (dum_params[:][model_list['Model#'] == 'Model'+str(imod)])[0].data

                if compare(ref_params, dum_params) == False:
                    # print dum_params
                    continue
                else:
                    p1.extend((model_list[keywords['col'][0]][model_list['Model#'] == 'Model'+str(imod)]).data)
                    if fixed == False:
                        p2.extend((model_list[keywords['col'][1]][model_list['Model#'] == 'Model'+str(imod)]).data)
                    model_label.append(str(imod))

            chi2.extend(reduced_chi2_dum)


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

        ax.plot(p1[np.argsort(p1)], chi2[np.argsort(p1)], 'o-', mec='None', color='Green', linewidth=1.5)
        ax.set_xlabel(keywords['label'][0], fontsize=18)
        # ax.set_ylabel(r'$\rm{\Sigma(sim.-obs.)^{2}/(\sigma_{data}^{2}+\sigma_{sys}^{2})}$', fontsize=18)
        ax.set_ylabel(r'$\rm{\chi^{2}_{reduced}}$', fontsize=18)

        ax.set_yscale('log')

        [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
        ax.minorticks_on() 
        ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
        ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)

        fig.savefig('/Users/yaolun/test/chi2_'+str(keywords['col'][0])+'_1d.pdf', format='pdf', dpi=300, bbox_inches='tight')
        fig.clf()

    else:
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

        z = griddata((p1_norm, p2_norm), chi2, (x[None,:], y[:,None]), method='cubic')
        if z.min() < 0:
            z = griddata((p1_norm, p2_norm), chi2, (x[None,:], y[:,None]), method='linear')

        # print z.min()

        # z = np.log10(z)

        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111)

        # plot the contour with color and lines
        ax.contour(x, y, z, 10, linewidths=0.5,colors='k')
        # cs = ax.contourf(x,y,z,15,cmap=plt.cm.jet)
        # cmap = sns.cubehelix_palette(light=1, as_cmap=True, reverse=True)
        cmap = plt.cm.CMRmap
        # import custom colormap
        from custom_colormap import custom_colormap
        cmap = mpl.colors.ListedColormap(custom_colormap())
        im = ax.imshow(z, cmap=cmap, origin='lower', extent=[0,1,0,1],\
            norm=LogNorm(vmin=chi2.min(), vmax=chi2.max()))
        # Blues_r
        ax.set_xticks(np.linspace(0, 1, 5))
        ax.set_xticklabels(np.linspace(min(p1), max(p1), 5))
        ax.set_yticks(np.linspace(0, 1, 5))
        ax.set_yticklabels(np.linspace(min(p2), max(p2), 5))

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = fig.colorbar(im, cax=cax)
        cb.solids.set_edgecolor("face")
        cb.ax.minorticks_on()
        # cb.ax.set_ylabel(r'$\rm{\Sigma(sim./obs.-1)^{2}/(\sigma_{combine}^{2})}$',fontsize=16)
        cb.ax.set_ylabel(r'$\rm{\chi^{2}_{reduce}}$', fontsize=20)
        cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
        plt.setp(cb_obj,fontsize=12)
        # plot the original data points
        ori_data = ax.scatter(p1_norm,p2_norm, marker='o',c='b',s=5)

        # print the model number near the points
        # for i in range(len(model_label)):
            # print model_label[i]
            # ax.annotate(model_label[i], (p1_norm[i], p2_norm[i]))

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




    # # plot the contour
    # p1 = np.squeeze(p1)
    # p2 = np.squeeze(p2)
    # chi2 = np.squeeze(chi2)
    # X, Y =  np.meshgrid(p1,p2)
    # Z = np.empty_like(X)
    # for i in range(0, len(Z[0])):
    #     for j in range(0, len(Z[1])):
    #         Z[i,j] = chi2[(p1 == X[i,j])*(p2 == Y[i,j]).T]
    #         # print X[i,j], Y[i,j], Z[i,j]

    # fig = plt.figure(figsize=(8,6))
    # ax = fig.add_subplot(111)

    # CS = ax.contourf(X/1e4, Y, Z, 10, # [-1, -0.1, 0, 0.1],
    #                     #alpha=0.5,
    #                     cmap=plt.cm.bone,
    #                     origin='lower')
    # # Make a colorbar for the ContourSet returned by the contourf call.
    # cbar = plt.colorbar(CS)
    # cbar.ax.set_ylabel(r'$\mathrm{Reduced~\chi^{2}}$', fontsize=16)
    # # Add the contour line levels to the colorbar
    # # cbar.add_lines(CS2)

    # ax.set_xlabel(keywords['label'][0], fontsize=20)
    # ax.set_ylabel(keywords['label'][1], fontsize=20)
    # [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
    # ax.minorticks_on() 
    # ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
    # ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)

    # fig.savefig('/Users/yaolun/test/chi2_agscs.pdf', format='pdf', dpi=300, bbox_inches='tight')

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
               'model_num': np.arange(1,126)}]
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
keywords_list = [{'col':['age','theta_cav'], 'label': [r'$\rm{age\,[10^{4}\,yr]}$', r'$\rm{\theta_{cav}\,[deg.]}$']},\
                 {'col':['age','view_angle'], 'label': [r'$\rm{age\,[10^{4}\,yr]}$', r'$\rm{\theta_{incl}\,[deg.]}$']},\
                 {'col':['view_angle','theta_cav'], 'label': [r'$\rm{\theta_{incl}\,[deg.]}$', r'$\rm{\theta_{cav}\,[deg.]}$']}]

obs = '/Users/yaolun/bhr71/obs_for_radmc/'

# for keywords in keywords_list:
#     p1, p2, chi2 = fir_chi2_2d(array_list, keywords, obs, ref=32)
#     fir_chi2_2d(array_list, keywords, obs)

# 1-D rho_cav_center
array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/cycle9/model_list.txt',
               'datapath': '/Users/yaolun/bhr71/hyperion/cycle9',
               'model_num': np.arange(34,43)}]
keywords = {'col':['rho_cav_center'], 'label': [r'$\rm{\rho_{cav}\,[g\,cm^{-3}]}$']}
fir_chi2_2d(array_list, keywords, obs, fixed=True, ref=34)

# 1-D inclination
array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/cycle9/model_list.txt',
               'datapath': '/Users/yaolun/bhr71/hyperion/cycle9',
               'model_num': np.hstack((np.arange(44,50), 34))}]
keywords = {'col':['view_angle'], 'label': [r'$\rm{\theta_{incl.}\,[deg.]}$']}
fir_chi2_2d(array_list, keywords, obs, fixed=True, ref=34)

# 1-D rho_cav_edge
array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/cycle9/model_list.txt',
               'datapath': '/Users/yaolun/bhr71/hyperion/cycle9',
               'model_num': np.hstack((np.arange(30,34), 14))}]
keywords = {'col':['rho_cav_edge'], 'label': [r'$\rm{R_{cav,\circ}\,[AU]}$']}
fir_chi2_2d(array_list, keywords, obs, fixed=True, ref=14)
