def extract_hyperion(filename,indir=None,outdir=None,dstar=178.0,aperture=None,save=True,filter_func=False,
    plot_all=False,clean=False,exclude_wl=[],log=True,image=True,obj='BHR71',print_data_w_aper=False):
    """
    filename: The path to Hyperion output file
    indir: The path to the directory which contains observations data
    outdir: The path to the directory for storing extracted plots and ASCII files
    """
    def l_bol(wl,fv,dstar=178.0):
        import numpy as np
        import astropy.constants as const
        # wavelength unit: um
        # Flux density unit: Jy
        # constants setup
        #
        c = const.c.cgs.value
        pc = const.pc.cgs.value
        PI = np.pi
        SL = const.L_sun.cgs.value
        # Convert the unit from Jy to erg s-1 cm-2 Hz-1
        fv = np.array(fv)*1e-23
        freq = c/(1e-4*np.array(wl))

        diff_dum = freq[1:]-freq[0:-1]
        freq_interpol = np.hstack((freq[0:-1]+diff_dum/2.0,freq[0:-1]+diff_dum/2.0,freq[0],freq[-1]))
        freq_interpol = freq_interpol[np.argsort(freq_interpol)[::-1]]
        fv_interpol = np.empty(len(freq_interpol))
        # calculate the histogram style of spectrum
        #
        for i in range(0,len(fv)):
            if i == 0:
                fv_interpol[i] = fv[i]
            else:
                fv_interpol[2*i-1] = fv[i-1]
                fv_interpol[2*i] = fv[i]
        fv_interpol[-1] = fv[-1]

        dv = freq_interpol[0:-1]-freq_interpol[1:]
        dv = np.delete(dv,np.where(dv==0))

        fv = fv[np.argsort(freq)]
        freq = freq[np.argsort(freq)]

        return (np.trapz(fv,freq)*4.*PI*(dstar*pc)**2)/SL

    # function for properly calculating uncertainty of spectrophotometry value
    def unc_spectrophoto(wl, unc, trans):
        # adopting smiliar procedure as Trapezoidal rule
        # (b-a) * [ f(a) + f(b) ] / 2
        #
        return ( np.sum( trans[:-1]**2 * unc[:-1]**2 * (wl[1:]-wl[:-1])**2 ) / np.trapz(trans, x=wl)**2 )**0.5

    # to avoid X server error
    import matplotlib as mpl
    mpl.use('Agg')
    #
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    from hyperion.model import ModelOutput, Model
    from scipy.interpolate import interp1d
    from hyperion.util.constants import pc, c, lsun, au
    from astropy.io import ascii
    import sys
    from phot_filter import phot_filter
    from get_obs import get_obs

    # Open the model
    m = ModelOutput(filename)

    # Read in the observation data and calculate the noise & variance
    if indir == None:
        indir = '/Users/yaolun/bhr71/'
    if outdir == None:
        outdir = '/Users/yaolun/bhr71/hyperion/'

    # assign the file name from the input file
    print_name = os.path.splitext(os.path.basename(filename))[0]

    # use a canned function to extract observational data
    obs_data = get_obs(indir, obj=obj)        # unit in um, Jy
    wl_tot, flux_tot, unc_tot = obs_data['spec']
    flux_tot = flux_tot*1e-23    # convert unit from Jy to erg s-1 cm-2 Hz-1
    unc_tot = unc_tot*1e-23
    l_bol_obs = l_bol(wl_tot, flux_tot*1e23, dstar=dstar)

    wl_phot, flux_phot, flux_sig_phot = obs_data['phot']
    flux_phot = flux_phot*1e-23   # convert unit from Jy to erg s-1 cm-2 Hz-1
    flux_sig_phot = flux_sig_phot*1e-23

    if aperture == None:
        aperture = {'wave': [3.6, 4.5, 5.8, 8.0, 8.5, 9, 9.7, 10, 10.5, 11, 16, 20, 24, 35, 70, 100, 160, 250, 350, 500, 850],\
                    'aperture': [7.2, 7.2, 7.2, 7.2, 7.2, 7.2, 7.2, 7.2, 7.2, 7.2, 20.4, 20.4, 20.4, 20.4, 24.5, 24.5, 24.5, 24.5, 24.5, 24.5, 24.5]}
    # assign wl_aper and aper from dictionary of aperture
    wl_aper = aperture['wave']
    aper    = aperture['aperture']
    # create the non-repetitive aperture list and index array
    aper_reduced = list(set(aper))
    index_reduced = np.arange(1, len(aper_reduced)+1)  # '+1': the zeroth slice corresponds to infinite aperture

    # Create the plot
    mag = 1.5
    fig = plt.figure(figsize=(8*mag,6*mag))
    ax_sed = fig.add_subplot(1, 1, 1)

    # Plot the observed SED
    if not clean:
        color_seq = ['Green','Red','Black']
    else:
        color_seq = ['DimGray','DimGray','DimGray']
    # plot the observations
    # plot in log scale
    if log:
        pacs, = ax_sed.plot(np.log10(wl_tot[(wl_tot>40) & (wl_tot<190.31)]),\
                            np.log10(c/(wl_tot[(wl_tot>40) & (wl_tot<190.31)]*1e-4)*flux_tot[(wl_tot>40) & (wl_tot<190.31)]),\
                            '-',color=color_seq[0],linewidth=1.5*mag, alpha=0.7)
        spire, = ax_sed.plot(np.log10(wl_tot[wl_tot > 194]),np.log10(c/(wl_tot[wl_tot > 194]*1e-4)*flux_tot[wl_tot > 194]),\
                            '-',color=color_seq[1],linewidth=1.5*mag, alpha=0.7)
        irs, = ax_sed.plot(np.log10(wl_tot[wl_tot < 40]),np.log10(c/(wl_tot[wl_tot < 40]*1e-4)*flux_tot[wl_tot < 40]),\
                            '-',color=color_seq[2],linewidth=1.5*mag, alpha=0.7)
        photometry, = ax_sed.plot(np.log10(wl_phot),np.log10(c/(wl_phot*1e-4)*flux_phot),'s',mfc='DimGray',mec='k',markersize=8)
        # plot the observed photometry data
        ax_sed.errorbar(np.log10(wl_phot),np.log10(c/(wl_phot*1e-4)*flux_phot),\
            yerr=[np.log10(c/(wl_phot*1e-4)*flux_phot)-np.log10(c/(wl_phot*1e-4)*(flux_phot-flux_sig_phot)),\
                  np.log10(c/(wl_phot*1e-4)*(flux_phot+flux_sig_phot))-np.log10(c/(wl_phot*1e-4)*flux_phot)],\
            fmt='s',mfc='DimGray',mec='k',markersize=8)
    # plot in normal scale
    else:
        pacs, = ax_sed.plot(np.log10(wl_tot[(wl_tot>40) & (wl_tot<190.31)]),\
                            c/(wl_tot[(wl_tot>40) & (wl_tot<190.31)]*1e-4)*flux_tot[(wl_tot>40) & (wl_tot<190.31)],\
                            '-',color=color_seq[0],linewidth=1.5*mag, alpha=0.7)
        spire, = ax_sed.plot(np.log10(wl_tot[wl_tot > 194]),c/(wl_tot[wl_tot > 194]*1e-4)*flux_tot[wl_tot > 194],\
                            '-',color=color_seq[1],linewidth=1.5*mag, alpha=0.7)
        irs, = ax_sed.plot(np.log10(wl_tot[wl_tot < 40]),c/(wl_tot[wl_tot < 40]*1e-4)*flux_tot[wl_tot < 40],\
                            '-',color=color_seq[2],linewidth=1.5*mag, alpha=0.7)
        photometry, = ax_sed.plot(wl_phot,c/(wl_phot*1e-4)*flux_phot,'s',mfc='DimGray',mec='k',markersize=8)
        # plot the observed photometry data
        ax_sed.errorbar(np.log10(wl_phot),c/(wl_phot*1e-4)*flux_phot,\
            yerr=[c/(wl_phot*1e-4)*flux_phot-c/(wl_phot*1e-4)*(flux_phot-flux_sig_phot),\
                  c/(wl_phot*1e-4)*(flux_phot+flux_sig_phot)-c/(wl_phot*1e-4)*flux_phot],\
            fmt='s',mfc='DimGray',mec='k',markersize=8)

    # if keyword 'clean' is not set, print L_bol derived from observations at upper right corner.
    if not clean:
        ax_sed.text(0.75,0.9,r'$\rm{L_{bol}= %5.2f L_{\odot}}$' % l_bol_obs,fontsize=mag*16,transform=ax_sed.transAxes)

    # getting SED with infinite aperture
    sed_inf = m.get_sed(group=0, inclination=0, aperture=-1, distance=dstar * pc, uncertainties=True)

    # plot the simulated SED with infinite aperture
    if clean == False:
        sim, = ax_sed.plot(np.log10(sed_inf.wav), np.log10(sed_inf.val), '-', color='GoldenRod', linewidth=0.5*mag)
        ax_sed.fill_between(np.log10(sed_inf.wav), np.log10(sed_inf.val-sed_inf.unc), np.log10(sed_inf.val+sed_inf.unc),\
            color='GoldenRod', alpha=0.5)

    # get fluxes with different apertures
    # this is non-reduced wavelength array because this is for printing out fluxes at all channels specified by users
    flux_aper = np.zeros_like(wl_aper, dtype=float)
    unc_aper = np.zeros_like(wl_aper, dtype=float)
    a = np.zeros_like(wl_aper) + 1
    color_list = plt.cm.jet(np.linspace(0, 1, len(wl_aper)+1))
    for i in range(0, len(wl_aper)):
        # occasionally users might want not to report some wavelength channels
        if wl_aper[i] in exclude_wl:
            continue
        # getting simulated SED from Hyperion output. (have to match with the reduced index)
        sed_dum = m.get_sed(group=index_reduced[np.where(aper_reduced == aper[i])], \
                            inclination=0, aperture=-1, distance=dstar * pc, uncertainties=True)
        # plot the whole SED from this aperture (optional)
        if plot_all == True:
            ax_sed.plot(np.log10(sed_dum.wav), np.log10(sed_dum.val),'-', color=color_list[i])
            ax_sed.fill_between(np.log10(sed_dum.wav), np.log10(sed_dum.val-sed_dum.unc), np.log10(sed_dum.val+sed_dum.unc),\
                color=color_list[i], alpha=0.5)
        # Extracting spectrophotometry values from simulated SED
        # Not using the photometry filer function to extract spectrophotometry values
        #
        # sort by wavelength first.  It seems needed with different python or numpy version
        sort_wl = np.argsort(sed_dum.wav)
        val_sort = sed_dum.val[sort_wl]
        unc_sort = sed_dum.unc[sort_wl]
        wav_sort = sed_dum.wav[sort_wl]
        # Before doing that, convert vSv to F_lambda
        flux_dum = val_sort / wav_sort
        unc_dum  = unc_sort / wav_sort
        if filter_func == False:
            # use a rectangle function the average the simulated SED
            # apply the spectral resolution
            if (wl_aper[i] < 50.) & (wl_aper[i] >= 5):
                res = 60.
            elif wl_aper[i] < 5:
                res = 10.
            else:
                res = 1000.
            ind = np.where((wav_sort < wl_aper[i]*(1+1./res)) & (wav_sort > wl_aper[i]*(1-1./res)))
            if len(ind[0]) != 0:
                flux_aper[i] = np.mean(flux_dum[ind])
                unc_aper[i]  = np.mean(unc_dum[ind])
            else:
                f = interp1d(wav_sort, flux_dum)
                f_unc = interp1d(wav_sort, unc_dum)
                flux_aper[i] = f(wl_aper[i])
                unc_aper[i]  = f_unc(wl_aper[i])
        # Using photometry filter function to extract spectrophotometry values
        else:
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
            elif wl_aper[i] == 850:
                fil_name = 'SCUBA 850WB'
            else:
                fil_name = None

            if fil_name != None:
                filter_func = phot_filter(fil_name)
                # Simulated SED should have enough wavelength coverage for applying photometry filters.
                f = interp1d(wav_sort, flux_dum)
                f_unc = interp1d(wav_sort, unc_dum)
                flux_aper[i] = np.trapz(f(filter_func['wave']/1e4)*filter_func['transmission'], x=filter_func['wave']/1e4 )/\
                                np.trapz(filter_func['transmission'], x=filter_func['wave']/1e4)
                # fix a bug
                unc_aper[i] = unc_spectrophoto(filter_func['wave']/1e4, f_unc(filter_func['wave']/1e4), filter_func['transmission'])
            else:
                # use a rectangle function the average the simulated SED
                # apply the spectral resolution
                if (wl_aper[i] < 50.) & (wl_aper[i] >= 5):
                    res = 60.
                elif wl_aper[i] < 5:
                    res = 10.
                else:
                    res = 1000.
                ind = np.where((wav_sort < wl_aper[i]*(1+1./res)) & (wav_sort > wl_aper[i]*(1-1./res)))
                if len(ind[0]) != 0:
                    flux_aper[i] = np.mean(flux_dum[ind])
                    unc_aper[i]  = np.mean(unc_dum[ind])
                else:
                    f = interp1d(wav_sort, flux_dum)
                    f_unc = interp1d(wav_sort, unc_dum)
                    flux_aper[i] = f(wl_aper[i])
                    unc_aper[i]  = f_unc(wl_aper[i])
    # temperory step: solve issue of uncertainty greater than the value
    for i in range(len(wl_aper)):
        if unc_aper[i] >= flux_aper[i]:
            unc_aper[i] = flux_aper[i] - 1e-20

    # perform the same procedure of flux extraction of aperture flux with observed spectra
    # wl_aper = np.array(wl_aper, dtype=float)
    obs_aper_wl = wl_aper[(wl_aper >= min(wl_tot)) & (wl_aper <= max(wl_tot))]
    obs_aper_flux = np.zeros_like(obs_aper_wl)
    obs_aper_unc = np.zeros_like(obs_aper_wl)
    # have change the simulation part to work in F_lambda for fliter convolution
    # flux_tot and unc_tot have units of erg/s/cm2/Hz.  Need to convert it to F_lambda (erg/s/cm2/um)
    fnu2fl = c/(wl_tot*1e-4)/wl_tot
    #
    # sed_tot = c/(wl_tot*1e-4)*flux_tot
    # sed_unc_tot = c/(wl_tot*1e-4)*unc_tot
    #
    # wl_tot and flux_tot are already hstacked and sorted by wavelength
    for i in range(0, len(obs_aper_wl)):
        # sometime users want not report some wavelength channels
        if obs_aper_wl[i] in exclude_wl:
            continue
        if filter_func == False:
            # use a rectangle function the average the simulated SED
            # apply the spectral resolution
            if (obs_aper_wl[i] < 50.) & (obs_aper_wl[i] >= 5):
                res = 60.
            elif obs_aper_wl[i] < 5:
                res = 10.
            else:
                res = 1000.
            ind = np.where((wl_tot < obs_aper_wl[i]*(1+1./res)) & (wl_tot > obs_aper_wl[i]*(1-1./res)))
            if len(ind[0]) != 0:
                obs_aper_flux[i] = np.mean(fnu2fl[ind]*flux_tot[ind])
                obs_aper_unc[i] = np.mean(fnu2fl[ind]*unc_tot[ind])
            else:
                f = interp1d(wl_tot, fnu2fl*flux_tot)
                f_unc = interp1d(wl_tot, fnu2fl*unc_tot)
                obs_aper_flux[i] = f(obs_aper_wl[i])
                obs_aper_unc[i] = f_unc(obs_aper_wl[i])
        else:
            # apply the filter function
            # decide the filter name
            if obs_aper_wl[i] == 70:
                fil_name = 'Herschel PACS 70um'
            elif obs_aper_wl[i] == 100:
                fil_name = 'Herschel PACS 100um'
            elif obs_aper_wl[i] == 160:
                fil_name = 'Herschel PACS 160um'
            elif obs_aper_wl[i] == 250:
                fil_name = 'Herschel SPIRE 250um'
            elif obs_aper_wl[i] == 350:
                fil_name = 'Herschel SPIRE 350um'
            elif obs_aper_wl[i] == 500:
                fil_name = 'Herschel SPIRE 500um'
            elif obs_aper_wl[i] == 3.6:
                fil_name = 'IRAC Channel 1'
            elif obs_aper_wl[i] == 4.5:
                fil_name = 'IRAC Channel 2'
            elif obs_aper_wl[i] == 5.8:
                fil_name = 'IRAC Channel 3'
            elif obs_aper_wl[i] == 8.0:
                fil_name = 'IRAC Channel 4'
            elif obs_aper_wl[i] == 24:
                fil_name = 'MIPS 24um'
            # elif obs_aper_wl[i] == 850:
            #     fil_name = 'SCUBA 850WB'
            # do not have SCUBA spectra
            else:
                fil_name = None

            if fil_name != None:
                filter_func = phot_filter(fil_name)
                # Observed SED needs to be trimmed before applying photometry filters
                filter_func = filter_func[(filter_func['wave']/1e4 >= min(wl_tot))*\
                                          ((filter_func['wave']/1e4 >= 54.8)+(filter_func['wave']/1e4 <= 36.0853))*\
                                          ((filter_func['wave']/1e4 <= 95.05)+(filter_func['wave']/1e4 >=103))*\
                                          ((filter_func['wave']/1e4 <= 190.31)+(filter_func['wave']/1e4 >= 195))*\
                                          (filter_func['wave']/1e4 <= max(wl_tot))]
                f = interp1d(wl_tot, fnu2fl*flux_tot)
                f_unc = interp1d(wl_tot, fnu2fl*unc_tot)
                obs_aper_flux[i] = np.trapz(f(filter_func['wave']/1e4)*filter_func['transmission'], x=filter_func['wave']/1e4)/\
                                   np.trapz(filter_func['transmission'], x=filter_func['wave']/1e4)
                # obs_aper_flux_unc[i] = abs(np.trapz((filter_func['wave']/1e4)**2, (f_unc(filter_func['wave']/1e4)*filter_func['transmission'])**2))**0.5 / abs(np.trapz(filter_func['wave']/1e4, filter_func['transmission']))
                # fix a bug
                obs_aper_unc[i] = unc_spectrophoto(filter_func['wave']/1e4, f_unc(filter_func['wave']/1e4), filter_func['transmission'])
            else:
                # use a rectangle function the average the simulated SED
                # apply the spectral resolution
                if (obs_aper_wl[i] < 50.) & (obs_aper_wl[i] >= 5):
                    res = 60.
                elif obs_aper_wl[i] < 5:
                    res = 10.
                else:
                    res = 1000.
                ind = np.where((wl_tot < obs_aper_wl[i]*(1+1./res)) & (wl_tot > obs_aper_wl[i]*(1-1./res)))
                if len(ind[0]) != 0:
                    obs_aper_flux[i] = np.mean(fnu2fl[ind]*flux_tot[ind])
                    obs_aper_unc[i] = np.mean(fnu2fl[ind]*unc_tot[ind])
                else:
                    f = interp1d(wl_tot, fnu2fl*flux_tot)
                    f_unc = interp1d(wl_tot, fnu2fl*unc_tot)
                    obs_aper_flux[i] = f(obs_aper_wl[i])
                    obs_aper_unc[i] = f_unc(obs_aper_wl[i])

    # plot the aperture-extracted spectrophotometry fluxes from observed spectra and simulations
    # in log-scale
    if log:
        aper_obs = ax_sed.errorbar(np.log10(obs_aper_wl), np.log10(obs_aper_flux * obs_aper_wl ),\
            yerr=[np.log10(obs_aper_flux*obs_aper_wl)-np.log10(obs_aper_flux*obs_aper_wl-obs_aper_unc*obs_aper_wl), np.log10(obs_aper_flux*obs_aper_wl+obs_aper_unc*obs_aper_wl)-np.log10(obs_aper_flux*obs_aper_wl)],\
            fmt='s', mec='None', mfc='r', markersize=10, linewidth=1.5, ecolor='Red', elinewidth=3, capthick=3, barsabove=True)
        aper = ax_sed.errorbar(np.log10(wl_aper),np.log10(flux_aper*wl_aper),\
            yerr=[np.log10(flux_aper*wl_aper)-np.log10(flux_aper*wl_aper-unc_aper*wl_aper), np.log10(flux_aper*wl_aper+unc_aper*wl_aper)-np.log10(flux_aper*wl_aper)],\
            fmt='o', mec='Blue', mfc='None', color='b',markersize=12, markeredgewidth=2.5, linewidth=1.7, ecolor='Blue', elinewidth=3, barsabove=True)
        ax_sed.set_ylim([-14,-7])
        ax_sed.set_xlim([0,3.2])
    # in normal scale
    else:
        aper_obs = ax_sed.errorbar(np.log10(obs_aper_wl), obs_aper_flux*obs_aper_wl, yerr=obs_aper_unc*obs_aper_wl,\
            fmt='s', mec='None', mfc='r', markersize=10, linewidth=1.5, ecolor='Red', elinewidth=3, capthick=3, barsabove=True)
        aper = ax_sed.errorbar(np.log10(wl_aper),flux_aper*wl_aper, yerr=unc_aper*wl_aper,\
            fmt='o', mec='Blue', mfc='None', color='b',markersize=12, markeredgewidth=2.5, linewidth=1.7, ecolor='Blue', elinewidth=3, barsabove=True)
        # ax_sed.set_xlim([1, 1000])
        ax_sed.set_xlim([0,3.2])
        # ax_sed.set_ylim([1e-14, 1e-8])
    # calculate the bolometric luminosity of the aperture
    # print flux_aper
    l_bol_sim = l_bol(wl_aper, flux_aper*wl_aper/(c/np.array(wl_aper)*1e4)*1e23, dstar=dstar)
    print 'Bolometric luminosity of simulated spectrum: %5.2f lsun' % l_bol_sim

    # print out the sed into ascii file for reading in later
    if save == True:
        # unapertured SED
        foo = open(outdir+print_name+'_sed_inf.txt','w')
        foo.write('%12s \t %12s \t %12s \n' % ('wave','vSv','sigma_vSv'))
        for i in range(0, len(sed_inf.wav)):
            foo.write('%12g \t %12g \t %12g \n' % (sed_inf.wav[i], sed_inf.val[i], sed_inf.unc[i]))
        foo.close()
        # SED with convolution of aperture sizes
        foo = open(outdir+print_name+'_sed_w_aperture.txt','w')
        foo.write('%12s \t %12s \t %12s \n' % ('wave','vSv','sigma_vSv'))
        for i in range(0, len(wl_aper)):
            foo.write('%12g \t %12g \t %12g \n' % (wl_aper[i], flux_aper[i]*wl_aper[i], unc_aper[i]*wl_aper[i]))
        foo.close()
        # print out the aperture-convolved fluxex from observations
        if print_data_w_aper:
            foo = open(outdir+print_name+'_obs_w_aperture.txt','w')
            foo.write('%12s \t %12s \t %12s \n' % ('wave','Jy','sigma_Jy'))
            for i in range(0, len(obs_aper_wl)):
                foo.write('%12g \t %12g \t %12g \n' % (obs_aper_wl[i], obs_aper_flux[i]*obs_aper_wl[i]/(c/obs_aper_wl[i]*1e4)*1e23, obs_aper_unc[i]*obs_aper_wl[i]/(c/obs_aper_wl[i]*1e4)*1e23))
            foo.close()
    # Read in and plot the simulated SED produced by RADMC-3D using the same parameters
    # [wl,fit] = np.genfromtxt(indir+'hyperion/radmc_comparison/spectrum.out',dtype='float',skip_header=3).T
    # l_bol_radmc = l_bol(wl,fit*1e23/dstar**2)
    # radmc, = ax_sed.plot(np.log10(wl),np.log10(c/(wl*1e-4)*fit/dstar**2),'-',color='DimGray', linewidth=1.5*mag, alpha=0.5)

    # print the L bol of the simulated SED (both Hyperion and RADMC-3D)
    # lg_sim = ax_sed.legend([sim,radmc],[r'$\rm{L_{bol,sim}=%5.2f\,L_{\odot},\,L_{center}=9.18\,L_{\odot}}$' % l_bol_sim, \
    #   r'$\rm{L_{bol,radmc3d}=%5.2f\,L_{\odot},\,L_{center}=9.18\,L_{\odot}}$' % l_bol_radmc],\
    #   loc='lower right',fontsize=mag*16)

    # read the input central luminosity by reading in the source information from output file
    dum = Model()
    dum.use_sources(filename)
    L_cen = dum.sources[0].luminosity/lsun

    # legend
    lg_data = ax_sed.legend([irs, photometry, aper, aper_obs],\
                            [r'$\rm{observation}$',\
                            r'$\rm{photometry}$',r'$\rm{F_{aper,sim}}$',r'$\rm{F_{aper,obs}}$'],\
                            loc='upper left',fontsize=14*mag,numpoints=1,framealpha=0.3)
    if clean == False:
        lg_sim = ax_sed.legend([sim],[r'$\rm{L_{bol,sim}=%5.2f\,L_{\odot},\,L_{center}=%5.2f\,L_{\odot}}$' % (l_bol_sim, L_cen)], \
                               loc='lower right',fontsize=mag*16)
        plt.gca().add_artist(lg_data)

    # plot setting
    ax_sed.set_xlabel(r'$\rm{log\,\lambda\,[{\mu}m]}$',fontsize=mag*20)
    ax_sed.set_ylabel(r'$\rm{log\,\nu S_{\nu}\,[erg\,cm^{-2}\,s^{-1}]}$',fontsize=mag*20)
    [ax_sed.spines[axis].set_linewidth(1.5*mag) for axis in ['top','bottom','left','right']]
    ax_sed.minorticks_on()
    ax_sed.tick_params('both',labelsize=mag*18,width=1.5*mag,which='major',pad=15,length=5*mag)
    ax_sed.tick_params('both',labelsize=mag*18,width=1.5*mag,which='minor',pad=15,length=2.5*mag)

    # fix the tick label font
    ticks_font = mpl.font_manager.FontProperties(family='STIXGeneral',size=mag*18)
    for label in ax_sed.get_xticklabels():
        label.set_fontproperties(ticks_font)
    for label in ax_sed.get_yticklabels():
        label.set_fontproperties(ticks_font)

    # Write out the plot
    fig.savefig(outdir+print_name+'_sed.pdf',format='pdf',dpi=300,bbox_inches='tight')
    fig.clf()

    # option for suppress image plotting (for speed)
    if image:
        # Package for matching the colorbar
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        # Extract the image for the first inclination, and scale to 300pc. We
        # have to specify group=1 as there is no image in group 0.
        image = m.get_image(group=len(aper_reduced)+1, inclination=0, distance=dstar * pc, units='MJy/sr')
        # image = m.get_image(group=14, inclination=0, distance=dstar * pc, units='MJy/sr')
        # Open figure and create axes
        # fig = plt.figure(figsize=(8, 8))
        fig, axarr = plt.subplots(3, 3, sharex='col', sharey='row',figsize=(13.5,12))

        # Pre-set maximum for colorscales
        VMAX = {}
        VMAX[100] = 10.
        VMAX[250] = 100.
        VMAX[500] = 2000.
        VMAX[1000] = 2000.

        # We will now show four sub-plots, each one for a different wavelength
        for i, wav in enumerate([3.6, 8.0, 9.7, 24, 40, 100, 250, 500, 1000]):

            # ax = fig.add_subplot(3, 3, i + 1)
            ax = axarr[i/3, i%3]

            # Find the closest wavelength
            iwav = np.argmin(np.abs(wav - image.wav))

            # Calculate the image width in arcseconds given the distance used above
            # get the max radius
            rmax = max(m.get_quantities().r_wall)
            w = np.degrees(rmax / image.distance) * 3600.

            # Image in the unit of MJy/sr
            # Change it into erg/s/cm2/Hz/sr
            factor = 1e-23*1e6
            # avoid zero in log
            # flip the image, because the setup of inclination is upside down
            val = image.val[::-1, :, iwav] * factor + 1e-30

            # This is the command to show the image. The parameters vmin and vmax are
            # the min and max levels for the colorscale (remove for default values).
            # cmap = sns.cubehelix_palette(start=0.1, rot=-0.7, gamma=0.2, as_cmap=True)
            cmap = plt.cm.CMRmap
            im = ax.imshow(np.log10(val), vmin= -22, vmax= -12,
                      cmap=cmap, origin='lower', extent=[-w, w, -w, w], aspect=1)

            # fix the tick label font
            ticks_font = mpl.font_manager.FontProperties(family='STIXGeneral',size=14)
            for label in ax.get_xticklabels():
                label.set_fontproperties(ticks_font)
            for label in ax.get_yticklabels():
                label.set_fontproperties(ticks_font)

            # Colorbar setting
            # create an axes on the right side of ax. The width of cax will be 5%
            # of ax and the padding between cax and ax will be fixed at 0.05 inch.
            if (i+1) % 3 == 0:
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cb = fig.colorbar(im, cax=cax)
                cb.solids.set_edgecolor("face")
                cb.ax.minorticks_on()
                cb.ax.set_ylabel(r'$\rm{log(I_{\nu})\,[erg\,s^{-1}\,cm^{-2}\,Hz^{-1}\,sr^{-1}]}$',fontsize=12)
                cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
                plt.setp(cb_obj,fontsize=12)
                # fix the tick label font
                ticks_font = mpl.font_manager.FontProperties(family='STIXGeneral',size=12)
                for label in cb.ax.get_yticklabels():
                    label.set_fontproperties(ticks_font)

            if (i+1) == 7:
                # Finalize the plot
                ax.set_xlabel(r'$\rm{RA\,Offset\,[arcsec]}$', fontsize=14)
                ax.set_ylabel(r'$\rm{Dec\,Offset\,[arcsec]}$', fontsize=14)

            ax.tick_params(axis='both', which='major', labelsize=16)
            ax.set_adjustable('box-forced')
            ax.text(0.7,0.88,str(wav) + r'$\rm{\,\mu m}$',fontsize=16,color='white', transform=ax.transAxes)

        fig.subplots_adjust(hspace=0,wspace=-0.2)

        # Adjust the spaces between the subplots
        # plt.tight_layout()
        fig.savefig(outdir+print_name+'_image_gridplot.png', format='png', dpi=300, bbox_inches='tight')
        fig.clf()

# indir = '/Users/yaolun/bhr71/obs_for_radmc/'
# outdir = '/Users/yaolun/bhr71/hyperion/'
# import numpy as np
# wl_aper, aper_arcsec = np.genfromtxt(indir+'aperture.txt', skip_header=1, dtype=float).T
# aperture = {'wave': wl_aper, 'aperture': aper_arcsec}
# extract_hyperion('/Users/yaolun/test/model224.rtout',indir=indir,outdir='/Users/yaolun/test/',\
#                  aperture=aperture,filter_func=True,plot_all=False,clean=True,image=False,print_data_w_aper=True)
# extract_hyperion('/Users/yaolun/test/model_test_1e4_ics_gra3opc.rtout',indir=indir,outdir='/Users/yaolun/test/',\
#                  wl_aper=wl_aper,filter_func=True,plot_all=False,clean=True)
# extract_hyperion('/Users/yaolun/test/model_test_1e4_ics_gra2opc.rtout',indir=indir,outdir='/Users/yaolun/test/',\
#                  wl_aper=wl_aper,filter_func=True,plot_all=False,clean=True)
