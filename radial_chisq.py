def radial_chisq(array_list, keywords, filename_ext, plotpath, rmax=None, ref=None):
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import ascii
    import copy

    # function for calculating chisq
    def chisq(obs, sim, r_in=None, rmax=None):
        """
        r_in is the radius array corresponds to each intensity value.
        It is the inner radius of the annulus used in photometry.

        rmax is the maximum inner radius that is allowed for chisq calculation.
        """
        chi2 = 0

        # trim the intensity list if specified
        if (r_in != None) and (rmax != None):
            obs = obs[r_in <= rmax]
            sim = sim[r_in <= rmax]

        val2 = (obs['I']-sim['I'])**2
        unc2 = obs['I_err']**2 + sim['I_err']**2

        chi2 = np.sum(val2/unc2)
        n = len(obs)

        return chi2, n

    p1 = []
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
            print model_dum
            obs = {'I': model_dum['I'].data, 'I_err': model_dum['I_err'].data}
            sim = {'I': model_dum['I_sim'].data, 'I_err': model_dum['I_sim_err'].data}
            chi2_dum, n = chisq(obs, sim, r_in=model_dum['r_in[arcsec]'], rmax=rmax)

            reduced_chi2_dum = chi2_dum/(n-2-1)
            total_chi2.extend(reduced_chi2_dum)

            if ref == None:
                # read the parameter values
                p1.extend((model_list[keywords['col'][0]][model_list['Model#'] == 'Model'+str(imod)]).data)
                model_label.append(str(imod))
                # print reduced_chi2_dum
                chi2.extend(reduced_chi2_dum)
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

                    # print reduced_chi2_dum
                    chi2.extend(reduced_chi2_dum)

import numpy as np
array_list = [{'listpath': '/Users/yaolun/bhr71/hyperion/model_list.txt',
               'datapath': '/Users/yaolun/bhr71/hyperion/',
               'model_num': np.arange(55,67)}]
keywords = {'col':['age'], 'label': [r'$\rm{t_{col}\,[10^{4}\,year]}$']}
filename_ext = '_radial_profile_160.0um'
radial_chisq(array_list, keywords, filename_ext, '/Users/yaolun/test/', rmax=140.0, ref=61)
