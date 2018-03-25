def get_obs_CDF(cdfdir, obj, spitzer_file=None, photfile=None):
    """
    obj input in uppercase.  But check the path to make sure.
    """
    import numpy as np
    from astropy.io import ascii
    
    def spitzer_unc(filename, R=60., width=2.5):
        """
        R is the resolving power (lambda/delta_lambda)
        width = number of resolution elements
        """

        irs = ascii.read(filename, data_start=2, header_start=None, comment='%')
        wl_irs, flux_irs = irs['col1'], irs['col2']
        # [wl_irs, flux_irs]= (np.genfromtxt(filename,skip_header=2,dtype='float').T)[0:2]
        # Remove points with zero or negative flux 
        ind = (flux_irs > 0) & (np.isnan(flux_irs) == False)
        wl_irs = wl_irs[ind]
        flux_irs = flux_irs[ind]
        unc_irs = np.empty_like(flux_irs)

        oversample = (wl_irs[1]-wl_irs[0] + wl_irs[2]-wl_irs[1])/2 / (wl_irs[1]/R)

        j = 0
        edge = []
        for i in range(len(wl_irs)):
            if (wl_irs[i]-width/2 * wl_irs[i]/R >= min(wl_irs)) and (wl_irs[i]+width/2 * wl_irs[i]/R <= max(wl_irs)):
                wl_dum = wl_irs[(wl_irs >= wl_irs[i]-width/2*wl_irs[i]/R) & (wl_irs <= wl_irs[i]+width/2*wl_irs[i]/R)]
                flux_dum = flux_irs[(wl_irs >= wl_irs[i]-width/2*wl_irs[i]/R) & (wl_irs <= wl_irs[i]+width/2*wl_irs[i]/R)]
                # return the coefficient, highest power first.
                fit_dum = np.polyfit(wl_dum, flux_dum, 3)
                base_dum = fit_dum[0]*wl_dum**3 + fit_dum[1]*wl_dum**2 + fit_dum[2]*wl_dum + fit_dum[3]


                unc_irs[i] = np.std(flux_dum-base_dum) / np.sqrt(oversample)
                if j == 0:
                    edge.append(unc_irs[i])
                j += 1
                edge_dum = unc_irs[i]
        edge.append(edge_dum)
        # print edge
        for i in range(len(wl_irs)):
            if wl_irs[i]-width/2 * wl_irs[i]/R < min(wl_irs):
                unc_irs[i] = edge[0]
            if wl_irs[i]+width/2 * wl_irs[i]/R > max(wl_irs):
                unc_irs[i] = edge[1]
            if flux_irs[i] - unc_irs[i] < 0:
                unc_irs[i] = 1/3. * flux_irs[i]

        return wl_irs, flux_irs, unc_irs


    output = {}


    # Read in Herschel data
    # TODO: case for the sources without advanced products.
    # continuum
    [wl_pacs,flux_pacs] = np.genfromtxt(cdfdir+obj+'/pacs/advanced_products/'+obj+'_pacs_weighted_continuum.txt',dtype='float',skip_header=1).T
    [wl_spire,flux_spire] = np.genfromtxt(cdfdir+obj+'/spire/advanced_products/'+obj+'_spire_corrected_continuum.txt',dtype='float',skip_header=1).T
    # noise spectra
    [wl_pacs_noise, flux_pacs_noise] = np.genfromtxt(cdfdir+obj+'/pacs/advanced_products/'+obj+'_pacs_weighted_residual_spectrum.txt',dtype='float',skip_header=1).T
    [wl_spire_noise,flux_spire_noise] = np.genfromtxt(cdfdir+obj+'/spire/advanced_products/'+obj+'_spire_corrected_residual_spectrum.txt',dtype='float',skip_header=1).T

    # Calculate the local variance (for spire), use the instrument uncertainty for pacs
    #
    wl_noise = [wl_pacs_noise, wl_spire_noise]
    flux_noise = [flux_pacs_noise, flux_spire_noise]
    sig_num = 20
    sigma_noise = []
    for i in range(0, len(wl_noise)):
        sigma_dum = np.zeros_like(wl_noise[i])
        for iwl in range(0, len(wl_noise[i])):
            if iwl < sig_num/2:
                sigma_dum[iwl] = np.std(np.hstack((flux_noise[i][0:int(sig_num/2)], flux_noise[i][0:int(sig_num/2)-iwl])))
            elif len(wl_noise[i])-iwl < sig_num/2:
                sigma_dum[iwl] = np.std(np.hstack((flux_noise[i][iwl:], flux_noise[i][len(wl_noise[i])-int(sig_num/2):])))
            else:
                sigma_dum[iwl] = np.std(flux_noise[i][iwl-int(sig_num/2):iwl+int(sig_num/2)])
        sigma_noise = np.hstack((sigma_noise, sigma_dum))


    # Read in Spitzer data
    if spitzer_file != None:
        wl_irs, flux_irs, unc_irs = spitzer_unc(spitzer_file)
        wl_spec = np.hstack((wl_irs, wl_pacs, wl_spire))
        flux_spec = np.hstack((flux_irs, flux_pacs, flux_spire))
        sigma_noise = np.hstack((unc_irs, sigma_noise))
    else:
        wl_spec = np.hstack((wl_pacs,wl_spire))
        flux_spec = np.hstack((flux_pacs,flux_spire))

    flux_spec = flux_spec[np.argsort(wl_spec)]
    sigma_noise = sigma_noise[np.argsort(wl_spec)]
    wl_spec = wl_spec[np.argsort(wl_spec)]
    
    # filter NaN value
    wl_spec = wl_spec[np.isnan(flux_spec) == False]
    sigma_noise = sigma_noise[np.isnan(flux_spec) == False]
    flux_spec = flux_spec[np.isnan(flux_spec) == False]
    
    output['spec'] = (wl_spec, flux_spec, sigma_noise)

    if photfile!= None:
        # Read in the photometry data
        phot = ascii.read(photfile, comment='%')
        # phot = np.genfromtxt(photfile, dtype=None, skip_header=1, comments='%')
        # wl_phot = []
        # flux_phot = []
        # flux_sig_phot = []
        # # note = []
        # for i in range(0,len(phot)):
        #     wl_phot.append(phot[i][0])
        #     flux_phot.append(phot[i][1])
        #     flux_sig_phot.append(phot[i][2])
        #     # note.append(phot[i][4])
        # wl_phot = np.array(wl_phot)
        # flux_phot = np.array(flux_phot)
        # flux_sig_phot = np.array(flux_sig_phot)
        wl_phot = phot['wavelength']
        flux_phot = phot['flux(Jy)']
        flux_sig_phot = phot['error(Jy)']
        selector = (wl_phot != 70) & (wl_phot != 100) & (wl_phot != 160) & (wl_phot != 250) & (wl_phot != 350) & (wl_phot != 500)
        wl_phot = wl_phot[selector]
        flux_phot = flux_phot[selector]
        flux_sig_phot = flux_sig_phot[selector]
        
        # Read in CDF photometry
        phot_pacs = ascii.read(cdfdir+obj+'/pacs/data/'+obj+'_pacs_phot.txt', data_start=4)
        phot_spire = ascii.read(cdfdir+obj+'/spire/data/'+obj+'_spire_phot.txt', data_start=4)
        # average the photometry
        phot_cdf = {'wave': [], 'flux': [], 'unc':[]}
        # PACS
        for i, w in enumerate(set(phot_pacs['wavelength(um)'])):
            phot_cdf['wave'].append(w)
            phot_cdf['flux'].append(np.mean(phot_pacs['flux(Jy)'][phot_pacs['wavelength(um)'] == w]))
            phot_cdf['unc'].append((np.sum(phot_pacs['uncertainty(Jy)'][phot_pacs['wavelength(um)'] == w]**2)/len(phot_pacs['uncertainty(Jy)'][phot_pacs['wavelength(um)'] == w]))**0.5)
        # SPIRE
        for i, w in enumerate(set(phot_spire['wavelength(um)'])):
            phot_cdf['wave'].append(w)
            phot_cdf['flux'].append(np.mean(phot_spire['flux(Jy)'][phot_spire['wavelength(um)'] == w]))
            phot_cdf['unc'].append((np.sum(phot_spire['uncertainty(Jy)'][phot_spire['wavelength(um)'] == w]**2)/len(phot_spire['uncertainty(Jy)'][phot_spire['wavelength(um)'] == w]))**0.5)

        # combine photoemtry
        wl_phot = np.hstack((wl_phot, np.array(phot_cdf['wave'])))
        flux_phot = np.hstack((flux_phot, np.array(phot_cdf['flux'])))
        flux_sig_phot = np.hstack((flux_sig_phot, np.array(phot_cdf['unc'])))
        
        # filter NaN values
        wl_phot = wl_phot[np.isnan(flux_phot) == False]
        flux_sig_phot = flux_sig_phot[np.isnan(flux_phot) == False]
        flux_phot = flux_phot[np.isnan(flux_phot) == False]
        
        output['phot'] = (wl_phot, flux_phot, flux_sig_phot)

    return output
