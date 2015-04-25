def fir_chi2_2d(array_list, keywords, obs, wl_aper=None):
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

    # constant setup
    c = const.c.cgs.value

    # chi2 function
    def fir_chi2(obs, sim, wave=[70.,100.,160.,250.,350.,500.]):

        chi2 = 0
        for w in wave:
            print w, (sim['sed'][sim['wave'] == w]-obs['sed'][obs['wave'] == w])**2
            chi2 = chi2 + ((sim['sed'][sim['wave'] == w]-obs['sed'][obs['wave'] == w])**2)# / (obs['sigma'][obs['wave'] == w])**2
        print chi2
        return chi2, len(wave)

    # setup the aperture size
    if wl_aper == None:
        wl_aper = [3.6, 4.5, 5.8, 8.0, 8.5, 9, 9.7, 10, 10.5, 11, 16, 20, 24, 35, 70, 100, 160, 250, 350, 500, 850]
        wl_aper = [70., 100., 160., 250., 350., 500.]

    # read the observed SED and extract with apertures
    bhr71 = get_bhr71_obs(obs)
    wave_obs, flux_obs, noise_obs = bhr71['spec']
    wave_obs = wave_obs[wave_obs > 50]
    flux_obs = flux_obs[wave_obs > 50]
    sed_obs = c/(wave_obs*1e-4)*flux_obs*1e-23
    sed_obs_noise = c/(wave_obs*1e-4)*noise_obs*1e-23
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

        for i in range(0, len(model_num)):
            imod = model_num[i]
            # read the parameter values
            p1.append(model_list[keywords['col'][0]][model_list['Model#'] == 'Model'+str(imod)])
            p2.append(model_list[keywords['col'][1]][model_list['Model#'] == 'Model'+str(imod)])

            # read the simulated SED
            model_dum = ascii.read(datapath+'/model'+str(imod)+'_sed_w_aperture.txt')
            print p1[-1], p2[-1]
            # plug them into the chi2 function
            chi2_dum, n = fir_chi2({'wave': np.array(wl_aper), 'sed': obs_aper_sed, 'sigma': sed_obs_noise}, {'wave': model_dum['wave'], 'sed': model_dum['vSv']}, wave=wl_aper)
            reduced_chi2_dum = chi2_dum/(n-2-1)
            chi2.append(reduced_chi2_dum)

    # plot the contour
    p1 = np.squeeze(p1)
    p2 = np.squeeze(p2)
    chi2 = np.squeeze(chi2)
    X, Y =  np.meshgrid(p1,p2)
    Z = np.empty_like(X)
    for i in range(0, len(Z[0])):
        for j in range(0, len(Z[1])):
            Z[i,j] = chi2[(p1 == X[i,j])*(p2 == Y[i,j]).T]
            # print X[i,j], Y[i,j], Z[i,j]

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)

    CS = ax.contourf(X/1e4, Y, Z, 10, # [-1, -0.1, 0, 0.1],
                        #alpha=0.5,
                        cmap=plt.cm.bone,
                        origin='lower')
    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = plt.colorbar(CS)
    cbar.ax.set_ylabel(r'$\mathrm{Reduced~\chi^{2}}$', fontsize=16)
    # Add the contour line levels to the colorbar
    # cbar.add_lines(CS2)

    ax.set_xlabel(keywords['label'][0], fontsize=20)
    ax.set_ylabel(keywords['label'][1], fontsize=20)
    [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
    ax.minorticks_on() 
    ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
    ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)

    fig.savefig('/Users/yaolun/test/chi2_agscs.pdf', format='pdf', dpi=300, bbox_inches='tight')

    return p1, p2, chi2

import numpy as np
array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/controlled/model_list.txt',
               'datapath': '/Users/yaolun/bhr71/hyperion/controlled',
               'model_num': np.arange(7,22)}]
               #,
              # {'listpath': '/Users/yaolun/bhr71/hyperion/cycle5/model_list.txt',
              #  'datapath': '/Users/yaolun/bhr71/hyperion/cycle5',
              #  'model_num': np.arange(38,46)}
keywords = {'col':['age','Cs'], 'label': [r'$\mathrm{age~[10^{4}~yr]}$', r'$\mathrm{c_{s}~[km~s^{-1}]}$']}
obs = '/Users/yaolun/bhr71/obs_for_radmc/'
p1, p2, chi2 = fir_chi2_2d(array_list, keywords, obs)
# for i in range(0, len(p2)):
    # print p1[i], p2[i], chi2[i]
