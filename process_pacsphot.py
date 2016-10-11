def proc_PACSphot(filename):
    from astropy.io import ascii
    import numpy as np
    data = ascii.read(filename, data_start=4)

    # construct non-redundant wavelength points
    wave = np.sort(list(set(data['wavelength(um)'].data)))
    phot = np.empty_like(wave)
    phot_err = np.empty_like(wave)

    # calculate the mean flux and uncertainty for each wavelength point
    for i in range(len(wave)):
        phot[i] = np.mean(data['flux(Jy)'][data['wavelength(um)'] == wave[i]])
        phot_err[i] = np.sqrt(np.sum(data['uncertainty(Jy)'][data['wavelength(um)'] == wave[i]]**2))/\
                      len(data['wavelength(um)'][data['wavelength(um)'] == wave[i]])

    return {'wavelength(um)': wave,
            'flux(Jy)': phot,
            'uncertainty(Jy)': phot_err}

# filename = '/Volumes/SD-Mac/pacsphot_cdf/BHR71_phot_sect.txt'
# proc_PACSphot(filename)
