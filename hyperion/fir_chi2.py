def fir_chi2_2d(array_list, keywords, obs, wl_aper=None, fixed_cs=False, ref=None):
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
    import copy
    import collections

    # function for checking duplicate in lists
    compare = lambda x, y: collections.Counter(x) == collections.Counter(y)

    # constant setup
    c = const.c.cgs.value

    # chi2 function
    def fir_chi2(obs, sim, wave=[70.,100.,160.,250.,350.,500.]):

        chi2 = 0
        for w in wave:
            # print w, (sim['sed'][sim['wave'] == w]-obs['sed'][obs['wave'] == w])**2
            chi2 = chi2 + ((sim['sed'][sim['wave'] == w]-obs['sed'][obs['wave'] == w])**2)# / (obs['sigma'][obs['wave'] == w])**2
        return chi2, len(wave)

    # setup the aperture size
    if wl_aper == None:
        wl_aper = [3.6, 4.5, 5.8, 8.0, 8.5, 9, 9.7, 10, 10.5, 11, 16, 20, 24, 35, 70, 100, 160, 250, 350, 500, 850]
        wl_aper = [5.8, 8.0, 8.5, 9, 9.7, 10, 10.5, 11, 16, 20, 24, 35, 70, 100, 160, 250, 350, 500]
        # wl_aper = [70., 100., 160., 250., 350., 500.]

    # read the observed SED and extract with apertures
    bhr71 = get_bhr71_obs(obs)
    wave_obs, flux_obs, sed_obs_noise = bhr71['spec']
    # wave_obs = wave_obs[wave_obs > 50]
    # flux_obs = flux_obs[wave_obs > 50]
    sed_obs = c/(wave_obs*1e-4)*flux_obs*1e-23
    # sed_obs_noise = c/(wave_obs*1e-4)*sed_obs_noise*1e-23
    obs_aper_sed = np.empty_like(wl_aper)

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
        else:
            f = interp1d(wave_obs, sed_obs)
            obs_aper_sed[i] = f(wl_aper[i])

    # print obs_aper_sed

    # calculate Chi2 from simulated SED
    p1 = []
    p2 = []
    chi2 = []

    num_file = len(array_list)
    for ifile in range(0,num_file):
        listpath = array_list[ifile]['listpath']
        datapath = array_list[ifile]['datapath']
        model_num = array_list[ifile]['model_num']

        model_list = ascii.read(listpath)

        if ref != None:
            ignore_col = ['d_sub', 'M_env_dot', 'R_inf', 'R_cen', 'mstar', 'total_mass']
            ref_params = copy.copy(model_list)
            ref_params.remove_columns(keywords['col'])
            ref_params.remove_columns(ignore_col)
            ref_params.remove_column('Model#')
            ref_params = (ref_params[:][model_list['Model#'] == 'Model'+str(ref)])[0].data
            # print ref_params

        for i in range(0, len(model_num)):
            imod = model_num[i]
            # manually exclude much older age
            if keywords['col'][0] == 'age':
                if (model_list[keywords['col'][0]][model_list['Model#'] == 'Model'+str(imod)]).data >= 5e5:
                    continue
            if keywords['col'][1] == 'age':
                if (model_list[keywords['col'][1]][model_list['Model#'] == 'Model'+str(imod)]).data >= 5e5:
                    continue                    
            if ref == None:
                # read the parameter values
                p1.extend((model_list[keywords['col'][0]][model_list['Model#'] == 'Model'+str(imod)]).data)
                p2.extend((model_list[keywords['col'][1]][model_list['Model#'] == 'Model'+str(imod)]).data)
            else:
                # get other parameters of model i
                dum_params = copy.copy(model_list)
                dum_params.remove_columns(keywords['col'])
                dum_params.remove_columns(ignore_col)
                dum_params.remove_column('Model#')
                dum_params = (dum_params[:][model_list['Model#'] == 'Model'+str(imod)])[0].data

                if compare(ref_params, dum_params) == False:
                    # print dum_params
                    continue

            # read the simulated SED
            model_dum = ascii.read(datapath+'/model'+str(imod)+'_sed_w_aperture.txt')
            # print p1[-1], p2[-1]
            # print imod
            # print model_dum
            # plug them into the chi2 function
            chi2_dum, n = fir_chi2({'wave': np.array(wl_aper), 'sed': obs_aper_sed, 'sigma': sed_obs_noise}, {'wave': model_dum['wave'].data, 'sed': model_dum['vSv'].data}, wave=wl_aper)
            reduced_chi2_dum = chi2_dum/(n-2-1)
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
            elif i ==1:
                p2 = p2/1e4

    # plot 1d relation
    if fixed_cs == True:
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)

        p1 = np.squeeze(p1)
        p2 = np.squeeze(p2)
        chi2 = np.squeeze(chi2)

        p1 = np.array(p1); chi2 = np.array(chi2)

        ax.plot(p1[np.argsort(p1)]/1e4, chi2[np.argsort(p1)], 'o-', mec='None', color='Green', linewidth=1.5)
        ax.set_xlabel(keywords['label'][0], fontsize=18)
        ax.set_ylabel(r'$\rm{\Sigma(sim.-obs.)^{2}}$', fontsize=18)

        [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
        ax.minorticks_on() 
        ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
        ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)

        fig.savefig('/Users/yaolun/test/chi2_agscs_1d.pdf', format='pdf', dpi=300, bbox_inches='tight')
        fig.clf()

    # rebin the data and plot 2D contour
    from scipy.interpolate import griddata
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    x = np.linspace(min(p1), max(p1), 50)
    y = np.linspace(min(p2), max(p2), 50)
    z = griddata((p1, p2), chi2, (x[None,:], y[:,None]), method='cubic')

    # z = np.log10(z)

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    # plot the contour with color and lines
    ax.contour((x-min(x))* 50/(max(x)-min(x)),\
                (y-min(y))* 50/(max(y)-min(y)),z,15,linewidths=0.5,colors='k', vmin=chi2.min())
    # cs = ax.contourf(x,y,z,15,cmap=plt.cm.jet)
    im = ax.imshow(z, cmap='Blues_r', origin='lower', vmin=chi2.min())
    ax.set_xticks(np.linspace(0, 50, 5))
    ax.set_xticklabels(np.linspace(min(p1), max(p1), 5))
    ax.set_yticks(np.linspace(0, 50, 5))
    ax.set_yticklabels(np.linspace(min(p2), max(p2), 5))

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = fig.colorbar(im, cax=cax)
    cb.solids.set_edgecolor("face")
    cb.ax.minorticks_on()
    cb.ax.set_ylabel(r'$\rm{\Sigma(sim.-obs.)^{2}}$',fontsize=16)
    cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
    plt.setp(cb_obj,fontsize=12)
    # plot the original data points
    ori_data = ax.scatter((p1-min(p1))* 50/(max(p1)-min(p1)),\
                          (p2-min(p2))* 50/(max(p2)-min(p2)), marker='o',c='b',s=5)

    ax.set_xlabel(keywords['label'][0], fontsize=20)
    ax.set_ylabel(keywords['label'][1], fontsize=20)
    [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
    ax.minorticks_on() 
    ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
    ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)

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
# array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/cycle6/model_list.txt',
#                'datapath': '/Users/yaolun/bhr71/hyperion/cycle6',
#                'model_num': np.arange(1,91)}]
array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/controlled/model_list.txt',
               'datapath': '/Users/yaolun/bhr71/hyperion/controlled',
               'model_num': np.arange(1,77)}]
# keywords = {'col':['age','view_angle'], 'label': [r'$\rm{age\,[10^{4}\,yr]}$', r'$\rm{viewing\,angle}$']}
# keywords = {'col':['age','theta_cav'], 'label': [r'$\rm{age\,[10^{4}\,yr]}$', r'$\rm{\theta_{cav}}$']}
keywords = {'col':['age','Cs'], 'label': [r'$\rm{age\,[10^{4}\,yr]}$', r'$\rm{c_{s}\,[km\,s^{-1}]}$']}
# keywords = {'col':['view_angle','theta_cav'], 'label': [r'$\rm{viewing\,angle}$', r'$\rm{\theta_{cav}}$']}


obs = '/Users/yaolun/bhr71/obs_for_radmc/'

# keywords = {'col':['age','Cs'], 'label': [r'$\mathrm{age~[10^{4}~yr]}$', r'$\mathrm{c_{s}~[km~s^{-1}]}$']}
# obs = '/Users/yaolun/bhr71/obs_for_radmc/'
p1, p2, chi2 = fir_chi2_2d(array_list, keywords, obs, ref=22)
for i in range(0, len(p2)):
    print p1[i], p2[i], chi2[i]
