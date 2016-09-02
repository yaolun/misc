def radial_chisq(array_list, keywords, filename_ext, plotpath, rmax=None, ref=None, zoom_1d=None):
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import ascii
    import copy
    import collections

    # function for calculating chisq
    def chisq(obs, sim, r_in=None, rmax=None):
        """
        r_in is the radius array corresponds to each intensity value.
        It is the inner radius of the annulus used in photometry.

        rmax is the maximum inner radius that is allowed for chisq calculation.

        The first elements should always be the same
        """
        chi2 = 0

        # trim the intensity list if specified
        if (r_in != None) and (rmax != None):
            ind = r_in <= rmax

        val2 = (obs['I'][ind]-sim['I'][ind])**2
        unc2 = obs['I_err'][ind]**2 + sim['I_err'][ind]**2

        # log version
        # import numpy as np
        # unc2_obs = (np.log10(obs['I'][ind]+obs['I_err'][ind])-np.log10(obs['I'][ind])) * (np.log10(obs['I'][ind])-np.log10(obs['I'][ind]-obs['I_err'][ind]))
        # unc2_sim = (np.log10(sim['I'][ind]+sim['I_err'][ind])-np.log10(sim['I'][ind])) * (np.log10(sim['I'][ind])-np.log10(sim['I'][ind]-sim['I_err'][ind]))
        # unc2 = unc2_obs + unc2_sim
        #
        # val2 = (np.log10(obs['I'][ind])-np.log10(sim['I'][ind]))**2
        #

        chi2 = np.nansum(val2[1:]/unc2[1:])
        n = len(val2)

        return chi2, n

    compare = lambda x, y: collections.Counter(x) == collections.Counter(y)


    p1 = []
    model_label = []
    chi2 = []
    total_chi2 = []

    num_file = len(array_list)
    for ifile in range(0,num_file):
        listpath = array_list[ifile]['listpath']
        datapath = array_list[ifile]['datapath']
        model_num = array_list[ifile]['model_num']

        model_list = ascii.read(listpath)

        # get the model parameters of the reference model
        if ref != None:
            ignore_col = ['d_sub', 'M_env_dot', 'R_inf', 'R_cen', 'mstar', 'M_tot_gas', 'M_disk']
            ref_params = copy.copy(model_list)
            ref_params.remove_column(keywords['col'][0])
            ref_params.remove_columns(ignore_col)
            ref_params.remove_column('Model#')
            ref_params = (ref_params[:][model_list['Model#'] == 'Model'+str(ref)])[0].data

        # iterate through the models to calculate chisq
        for i in range(0, len(model_num)):
            imod = model_num[i]
            model_dum = ascii.read(datapath+'/model'+str(imod)+filename_ext+'.txt')
            # exclude the first data point (smallest radius), because it is where two profiles
            # normalized.
            obs = {'I': model_dum['I'].data/model_dum['I'].data.max(),
                   'I_err': model_dum['I_err'].data/model_dum['I'].data.max()}
            sim = {'I': model_dum['I_sim'].data/model_dum['I_sim'].data.max(),
                   'I_err': model_dum['I_sim_err'].data/model_dum['I_sim'].data.max()}
            chi2_dum, n = chisq(obs, sim, r_in=model_dum['r_in'].data, rmax=rmax)

            # reduced_chi2_dum = chi2_dum/(n-2-1)
            reduced_chi2_dum = chi2_dum / (n-1)
            total_chi2.append(reduced_chi2_dum)

            if ref == None:
                # read the parameter values
                p1.extend((model_list[keywords['col'][0]][model_list['Model#'] == 'Model'+str(imod)]).data)
                model_label.append(str(imod))
                # print reduced_chi2_dum
                chi2.append(reduced_chi2_dum)
            else:
                # get other parameters of model i
                dum_params = copy.copy(model_list)
                # clean up test parameters and irrelvant parameters
                dum_params.remove_column(keywords['col'][0])
                dum_params.remove_columns(ignore_col)
                dum_params.remove_column('Model#')
                dum_params = (dum_params[:][model_list['Model#'] == 'Model'+str(imod)])[0].data
                # compare the controlled parameters with the reference
                if compare(ref_params, dum_params) == False:
                    continue
                else:
                    p1.extend((model_list[keywords['col'][0]][model_list['Model#'] == 'Model'+str(imod)]).data)
                    if imod == ref:
                        ref_p1 = float((model_list[keywords['col'][0]][model_list['Model#'] == 'Model'+str(imod)]).data)
                    model_label.append(str(imod))

                    print p1[-1], reduced_chi2_dum
                    chi2.append(reduced_chi2_dum)

    # plot the simulation on top of the observation
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)

    if keywords['col'][0] == 'age':
        p1 = np.array(np.squeeze(p1))/1e4
    else:
        p1 = np.array(np.squeeze(p1))

    chi2 = np.array(np.squeeze(chi2))

    print p1[chi2 == chi2.min()]

    ax.plot(p1[np.argsort(p1)], chi2[np.argsort(p1)], 'o-', mec='None', color='Green', linewidth=1, markersize=4)
    ax.set_xlabel(keywords['label'][0], fontsize=18)
    ax.set_ylabel(r'$\rm{\chi^{2}_{reduced}}$', fontsize=18)

    ax.axvline(2.2362, color='k', linestyle='--', linewidth=1)
    ax.axhline(1, color='k')
    ax.axhline(2, color='k', linestyle=':')
    # ax.axvspan(min(p1[chi2 <= min(chi2)*2]), max(p1[chi2 <= min(chi2)*2]),
    #            color='b', alpha=0.3)

    # ax.set_yscale('log')

    if zoom_1d != None:
        ax.set_xlim(zoom_1d)
    # else:
    #     # fig.gca().set_xlim(left=0)
    #     ax.set_xlim([0,10])
    ax.set_ylim([0, 10])

    [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
    ax.minorticks_on()
    ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
    ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)

    fig.savefig('/Users/yaolun/test/radial_chi2_'+str(keywords['col'][0])+'_1d.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()

import numpy as np
array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/controlled/model_list.txt',
               'datapath': '/Users/yaolun/bhr71/hyperion/controlled/',
               'model_num': np.arange(99,133)}]
keywords = {'col':['age'], 'label': [r'$\rm{t_{col}\,[10^{4}\,year]}$']}
# array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/model_list.txt',
#                'datapath': '/Users/yaolun/bhr71/hyperion/',
#                'model_num': np.arange(96,103)}]
# keywords = {'col':['view_angle'], 'label': [r'$\rm{\theta_{incl.}\,[deg.]}$']}
filename_ext = '_radial_profile_160.0um'
radial_chisq(array_list, keywords, filename_ext, '/Users/yaolun/test/', rmax=99.0, ref=115, zoom_1d=[0,7])
