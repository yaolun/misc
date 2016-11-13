def PacsSpire_SpecMatch(pacs, spire, threshold):
    """
    Compare the median flux at the last/first 5 um of PACS and SPIRE spectra.
    If both of them agree with each other within a certain threshold, return 0.
    Threshold is specified in percentage.
    If PACS flux is lower, return -1; if PACS flux is higher, return 1.
    Use the spectrum that has lower flux level as the base for threshold calculation.
    """
    import numpy as np

    pacs_level = np.nanmedian(pacs['Flux_Density(Jy)'][pacs['Wavelength(um)'] >= (pacs['Wavelength(um)'].max()-5)])
    spire_level = np.nanmedian(spire['Flux_Density(Jy)'][spire['Wavelength(um)'] <= (spire['Wavelength(um)'].min()+5)])

    base = min(pacs_level, spire_level)
    if abs(spire_level-pacs_level) <= threshold*base:
        return 0
    else:
        if spire_level > pacs_level:
            return 1
        else:
            return -1


# test script
# from astropy.io import ascii
# pacs = ascii.read('/Users/yaolun/bhr71/best_calibrated/BHR71_pacs_weighted.txt')
# spire = ascii.read('/Users/yaolun/bhr71/best_calibrated/BHR71_spire_corrected.txt')
# ans = PacsSpire_SpecMatch(pacs, spire, 0.1)
# print ans
