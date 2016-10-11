def ConvolveSimulation(spec, conv_wl, filter_dir, filter_func=True):
    """
    spec = [wave_grid, vSv, unc_vSv]
    conv_wl = a list of wavelength channels for convolution
    output: wl, vSv, unc_vSv
    """
    import numpy as np
    from scipy.interpolate import interp1d
    from phot_filter import phot_filter

    # function for properly calculating uncertainty of spectrophotometry value
    def unc_spectrophoto(wl, unc, trans):
        # adopting smiliar procedure as Trapezoidal rule
        # (b-a) * [ f(a) + f(b) ] / 2
        #
        return ( np.sum( trans[:-1]**2 * unc[:-1]**2 * (wl[1:]-wl[:-1])**2 ) / np.trapz(trans, x=wl)**2 )**0.5

    #######################################
    # get fluxes with different apertures #
    #######################################
    # this is non-reduced wavelength array because this is for printing out fluxes at all channels specified by users
    conv_flux = np.zeros_like(conv_wl, dtype=float)
    conv_unc = np.zeros_like(conv_wl, dtype=float)
    a = np.zeros_like(conv_wl) + 1

    # sort by wavelength first.
    sort_wl = np.argsort(spec[0])
    val_sort = spec[1][sort_wl]
    unc_sort = spec[2][sort_wl]
    wav_sort = spec[0][sort_wl]

    # Before doing that, convert vSv to F_lambda
    flux_dum = val_sort / wav_sort
    unc_dum  = unc_sort / wav_sort

    for i in range(0, len(conv_wl)):

        # If no using filter function to extract the spectrophotometry,
        # then use the spectral resolution.
        if filter_func == False:
            # use a rectangle function the average the simulated SED
            # apply the spectral resolution
            if (conv_wl[i] < 50.) & (conv_wl[i] >= 5):
                res = 60.
            elif conv_wl[i] < 5:
                res = 10.
            else:
                res = 1000.
            ind = np.where((wav_sort < conv_wl[i]*(1+1./res)) & (wav_sort > conv_wl[i]*(1-1./res)))
            if len(ind[0]) != 0:
                conv_flux[i] = np.mean(flux_dum[ind])
                conv_unc[i]  = np.mean(unc_dum[ind])
            else:
                f = interp1d(wav_sort, flux_dum)
                f_unc = interp1d(wav_sort, unc_dum)
                conv_flux[i] = f(conv_wl[i])
                conv_unc[i]  = f_unc(conv_wl[i])
        # Using photometry filter function to extract spectrophotometry values
        else:
            # apply the filter function
            # decide the filter name
            if conv_wl[i] == 70:
                fil_name = 'Herschel PACS 70um'
            elif conv_wl[i] == 100:
                fil_name = 'Herschel PACS 100um'
            elif conv_wl[i] == 160:
                fil_name = 'Herschel PACS 160um'
            elif conv_wl[i] == 250:
                fil_name = 'Herschel SPIRE 250um'
            elif conv_wl[i] == 350:
                fil_name = 'Herschel SPIRE 350um'
            elif conv_wl[i] == 500:
                fil_name = 'Herschel SPIRE 500um'
            elif conv_wl[i] == 3.6:
                fil_name = 'IRAC Channel 1'
            elif conv_wl[i] == 4.5:
                fil_name = 'IRAC Channel 2'
            elif conv_wl[i] == 5.8:
                fil_name = 'IRAC Channel 3'
            elif conv_wl[i] == 8.0:
                fil_name = 'IRAC Channel 4'
            elif conv_wl[i] == 24:
                fil_name = 'MIPS 24um'
            elif conv_wl[i] == 850:
                fil_name = 'SCUBA 850WB'
            else:
                fil_name = None

            if fil_name != None:
                filter_func = phot_filter(fil_name, filter_dir)
                # Simulated SED should have enough wavelength coverage for applying photometry filters.
                f = interp1d(wav_sort, flux_dum)
                f_unc = interp1d(wav_sort, unc_dum)

                if (filter_func['wave'].min() < wav_sort.min()*1e4):
                    print fil_name, 'has wavelength smaller than input spectrum.  The filter function is trimmed.'
                    filter_func = filter_func[filter_func['wave']/1e4 >= wav_sort.min()]
                if (filter_func['wave'].max() > wav_sort.max()*1e4):
                    print fil_name, 'has wavelength greater than input spectrum.  The filter function is trimmed.'
                    filter_func = filter_func[filter_func['wave']/1e4 <= wav_sort.max()]

                conv_flux[i] = np.trapz(f(filter_func['wave']/1e4)*\
                                          filter_func['transmission'],x=filter_func['wave']/1e4 )/\
                               np.trapz(filter_func['transmission'], x=filter_func['wave']/1e4)
                conv_unc[i] = unc_spectrophoto(filter_func['wave']/1e4,
                                    f_unc(filter_func['wave']/1e4), filter_func['transmission'])
            else:
                # use a rectangle function the average the simulated SED
                # apply the spectral resolution
                if (conv_wl[i] < 50.) & (conv_wl[i] >= 5):
                    res = 60.
                elif conv_wl[i] < 5:
                    res = 10.
                else:
                    res = 1000.
                ind = np.where((wav_sort < conv_wl[i]*(1+1./res)) & (wav_sort > conv_wl[i]*(1-1./res)))
                if len(ind[0]) != 0:
                    conv_flux[i] = np.mean(flux_dum[ind])
                    conv_unc[i]  = np.mean(unc_dum[ind])
                else:
                    f = interp1d(wav_sort, flux_dum)
                    f_unc = interp1d(wav_sort, unc_dum)
                    conv_flux[i] = f(conv_wl[i])
                    conv_unc[i]  = f_unc(conv_wl[i])
    # temperory step: solve issue of uncertainty greater than the value
    for i in range(len(conv_wl)):
        if conv_unc[i] >= conv_flux[i]:
            conv_unc[i] = conv_flux[i] - 1e-20

    return (conv_wl, conv_flux*conv_wl, conv_unc*conv_wl)

# import numpy as np
# wl, flux, unc = np.genfromtxt('/Users/yaolun/bhr71/hyperion/controlled/model80_sed_inf.txt', skip_header=1).T
# conv_wl = [3.6, 4.5, 5.8, 8.0, 8.5, 9, 9.7, 10, 10.5, 11, 16, 20, 24, 30, 70, 100, 160, 250, 350, 500, 1300]
# filter_dir = '/Users/yaolun/programs/misc/hyperion/'
# conv = ConvolveSimulation((wl, flux, unc), conv_wl, filter_dir, filter_func=True)
# print conv[0]
# print conv[1]
