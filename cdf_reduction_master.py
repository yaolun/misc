# Master script for CDF reduction with the method used in BHR71 paper.

# library import
import os
import pidly
idl = pidly.IDL('/opt/local/exelis/idl83/bin/idl')
from astropy.io import ascii
import sys
sys.path.append(os.path.expanduser('~')+'/programs/line_fitting/')
sys.path.append(os.path.expanduser('~')+'/programs/spectra_analysis/')
from PreFittedModify import *
from spire_spectral_index import spire_spectral_index
from cdfPacs1d import cdfPacs1d

# some info about the observational programs

# COPS-SPIRE source list
obsid_spire = [1342242620,1342242621,1342245084,1342245094,1342245857,
               1342247625,1342248246,1342248249,1342249053,1342249470,
               1342249474,1342249475,1342249476,1342249477,1342250509,
               1342250510,1342250512,1342250515,1342251285,1342251286,
               1342251287,1342251290,1342253646,1342253649,1342253652,
               1342254037]

# SECT cannot converage at L1489 1342249473, L1527 1342250511, HH100 1342252897
# mapping observation IRS44/46 1342251289

obj_list_spire = ['RCrA-IRS7B','RCrA-IRS7C','HH46','L723-MM','L1014',
                  'L1157','Ced110','BHR71','IRAS03245','L1551-IRS5',
                  'L1455-IRS3','B1-a','B1-c','IRAS03301','TMR1',
                  'TMC1A','TMC1','IRAS15398','RNO91','GSS30-IRS1',
                  'VLA1623','WL12','RCrA-IRS5A','L483','B335',
                  'DKCha']

# all source info
obsid = [['AB_Aur','1342217842','1342217843','0'],\
         ['AS205','1342215737','1342215738','0'],\
         ['B1-a','1342216182','1342216183','1342249475'],\
         ['B1-c','1342216213','1342216214','1342249476'],\
         ['B335','1342208889','1342208888','1342253652'],\
         ['BHR71','1342212230','1342212231','1342248249'],\
         ['Ced110','0','0','1342248246'],\
         ['DG_Tau','1342225730','1342225731','0'],\
         ['EC82','1342192975','1342219435','0'],\
         ['Elias29','1342228519','1342228520','0'],\
         ['FUOri','1342250907','1342250908','1342230412'],\
         ['GSS30-IRS1','1342215678','1342215679','0'],\
         ['GSS30','0','0','1342251286'],\
         ['HD100453','1342211695','1342211696','0'],\
         ['HD100546','1342188037','1342188038','0'],\
         ['HD104237','1342207819','1342207820','0'],\
         ['HD135344B-1','1342213921','1342213922','0'],\
         ['HD139614','1342215683','1342215684','0'],\
         ['HD141569','1342213913','0','0'],\
         ['HD142527','1342216174','1342216175','0'],\
         ['HD142666','1342213916','0','0'],\
         ['HD144432','1342213919','0','0'],\
         ['HD144668','1342215641','1342215642','0'],\
         ['HD150193','1342227068','0','0'],\
         ['HD163296','1342217819','1342217820','0'],\
         ['HD169142','1342206987','1342206988','0'],\
         ['HD179218','1342208884','1342208885','0'],\
         ['HD203024','1342206975','0','0'],\
         ['HD245906','1342228528','0','0'],\
         ['HD35187','1342217846','0','0'],\
         ['HD36112','1342228247','1342228248','0'],\
         ['HD38120','1342226212','1342226213','0'],\
         ['HD50138','1342206991','1342206992','0'],\
         ['HD97048','1342199412','1342199413','0'],\
         ['HD98922','1342210385','0','0'],\
         ['HH46','0','0','1342245084'],\
         ['HH100','0','0','1342252897'],\
         ['HT_Lup','1342213920','0','0'],\
         ['IRAM04191','1342216654','1342216655','0'],\
         ['IRAS03245','1342214677','1342214676','1342249053'],\
         ['IRAS03301','1342215668','1342216181','1342249477'],\
         ['DKCha','1342188039','1342188040','1342254037'],\
         ['IRAS15398','0','0','1342250515'],\
         ['IRS46','1342228474','1342228475','1342251289'],\
         ['IRS48','1342227069','1342227070','0'],\
         ['IRS63','1342228473','1342228472','0'],\
         ['L1014','1342208911','1342208912','1342245857'],\
         ['L1157','1342208909','1342208908','1342247625'],\
         ['L1448-MM','1342213683','1342214675','0'],\
         ['L1455-IRS3','1342204122','1342204123','1342249474'],\
         ['L1489','1342216216','1342216215','0'],\
         ['L1527','1342192981','1342192982','0'],\
         ['L1551-IRS5','1342192805','1342229711','1342249470'],\
         ['L483','0','0','1342253649'],\
         ['L723-MM','0','0','1342245094'],\
         ['RCrA-IRS5A','1342207806','1342207805','1342253646'],\
         ['RCrA-IRS7B','1342207807','1342207808','1342242620'],\
         ['RCrA-IRS7C','1342206990','1342206989','1342242621'],\
         ['RNO90','1342228206','0','0'],\
         ['RNO91','0','0','1342251285'],\
         ['RU_Lup','1342215682','0','0'],\
         ['RY_Lup','1342216171','0','0'],\
         ['S_Cra','1342207809','1342207810','0'],\
         ['SR21','1342227209','1342227210','0'],\
         ['Serpens-SMM3','1342193216','1342193214','0'],\
         ['Serpens-SMM4','1342193217','1342193215','0'],\
         ['TMC1','1342225803','1342225804','1342250512'],\
         ['TMC1A','1342192987','1342192988','1342250510'],\
         ['TMR1','1342192985','1342192986','1342250509'],\
         ['V1057_Cyg','1342235853','1342235852','1342221695'],\
         ['V1331_Cyg','1342233446','1342233445','1342221694'],\
         ['V1515_Cyg','1342235691','1342235690','1342221685'],\
         ['V1735_Cyg','1342235849','1342235848','1342219560'],\
         ['VLA1623','1342213918','1342213917','1342251287'],\
         ['WL12','1342228187','1342228188','1342251290']]

# reduction parameters
outdir = '/home/bettyjo/yaolun/CDF_archive_v2/'
pacsdatadir = '/scratch/CDF_PACS_HSA/'
reduction_name = 'CDF_archive_v2'
if not os.path.exists(outdir):
    os.mkdir(outdir)

# create the header in the output file for line fitting
# line strength reported in flux unit.
header = ['Object','Line','LabWL(um)','ObsWL(um)','Sig_Cen(um)','Str(W/cm2)',
          'Sig_str(W/cm2)','FWHM(um)','Sig_FWHM(um)','Base(W/cm2/um)',
          'Noise(W/cm2/um)','SNR','E_u(K)','A(s-1)','g','RA(deg)','Dec(deg)',
          'Pixel_No.','Blend','Validity']
for element in header:
    header[header.index(element)] = element + '  '
foo = open(outdir+reduction_name+'_lines.txt', 'w')
foo.write('{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s} \n'.format(*header))
foo.close()


# SPIRE reduction first for matching the PACS 1-D spectra with SECT-corrected spectra

###################### STEP 1 ######################
# HIPE execution:
#   "Spectrometer_Point_Pipeline_CDF.py" in hipe_script folder
#   Make sure the source list provided in above is the same as the source list defined in the script used here.
#   Check "outdir" to be the same as the outdir for the whole reduction.
#   This script will output cube FITS files reduced by 4 standard SPIRE options, SECT-reduced FITS file, SECT-reduced ASCII file, and fitted size/phot_obsid for further use.

###################### STEP 2 ######################
# Parse the cube FITS file into individual ASCII files.
# Use the "HR_spectrum_extened" product by default
spirecubever = 'HR_spectrum_extened'
idl('.r /home/bettyjo/yaolun/programs/line_fitting/get_spire.pro')
idl('.r /home/bettyjo/yaolun/programs/line_fitting/get_radec_spire_py.pro')

for o in obsid_spire:
    obj = obj_list_spire[obsid_spire.index(o)]
    print 'Step 2 - ', obj
    idl.pro('get_spire', outdir=outdir+obj+'/spire/data/cube/', object=obj,
            filename=outdir+obj+'/spire/data/fits/'+o+'_HR_spectrum_extended_apod.fits',
            brightness=1)
    # get RA and Dec
    idl.pro('get_radec_spire', filename=outdir+obj+'/spire/data/fits/'+o+'_HR_spectrum_extended_apod.fits',
            slw=1, write=outdir+obj+'/spire/data/cube/')
    idl.pro('get_radec_spire', filename=outdir+obj+'/spire/data/fits/'+o+'_HR_spectrum_extended_apod.fits',
            ssw=1, write=outdir+obj+'/spire/data/cube/'+obj)
    # the output RA/Dec files are named as obj+radec_slw[ssw].txt

###################### STEP 3 ######################
# line fitting on cube products
idl('.r /home/bettyjo/yaolun/programs/line_fitting/extract_spire.pro')
idl('.r /home/bettyjo/yaolun/programs/gauss.pro')

for o in obsid_spire:
    obj = obj_list_spire[obsid_spire.index(o)]
    print 'Step 3 - ', obj
    # read in RA/Dec
    radec_slw = ascii.read(outdir+obj+'/spire/data/cube/'+obj+'_radec_slw.txt')
    radec_ssw = ascii.read(outdir+obj+'/spire/data/cube/'+obj+'_radec_slw.txt')
    # SLW
    idl.pro('extract_spire', indir=outdir+obj+'/spire/data/cube/', outdir=outdir+obj+'/spire/advanced_products/cube/',
            plotdir=outdir+obj+'/spire/advanced_products/cube/plots/', localbaseline=10, global_noise=20,
            ra=radec_slw['RA(deg)'].data, dec=radec_slw['Dec(deg)'].data, coordpix=radec_slw['Pixel'].data,
            slw=1, noiselevel=3, brightness=1, object=obj, flat=1, continuum_sub=1, current_pix=1, double_gauss=1,
            print_all=outdir+reduction_name+'_lines.txt')
    # SSW
    idl.pro('extract_spire', indir=outdir+obj+'/spire/data/cube/', outdir=outdir+obj+'/spire/advanced_products/cube/',
            plotdir=outdir+obj+'/spire/advanced_products/cube/plots/', localbaseline=10, global_noise=20,
            ra=radec_ssw['RA(deg)'].data, dec=radec_ssw['Dec(deg)'].data, coordpix=radec_ssw['Pixel'].data,
            ssw=1, noiselevel=3, brightness=1, object=obj, flat=1, continuum_sub=1, current_pix=1, double_gauss=1,
            print_all=outdir+reduction_name+'_lines.txt')

###################### STEP 4 ######################
# re-format the SECT-reduced 1-D product and perform fitting
print 'Step 4'
SPIRE1D_run(obsid=obsid, indir=outdir+obj+'/spire/', outdir=outdir+obj+'/spire/', global_dir=outdir+reduction_name)

###################### STEP 5 ######################
# Fit alpha for three photometric bands for later measuring photometric fluxes.
for o in obsid_spire:
    obj = obj_list_spire[obsid_spire.index(o)]
    print 'Step 5 - ', obj
    spire_spectral_index(outdir+obj+'/spire/data/', o, obj)

###################### STEP 6 ######################
# [HIPE] execution
#   "spire_phot.py" in hipe_script folder
#   measure the photometry fluxes at 250 um, 350 um, and 500 um with the convolved apeture sizes

sys.exit("SPIRE part successfully finished!")

# PACS reduction
###################### STEP 1 ######################
# [HIPE] execution
#   getPacsDIGIT.py in the home directory of bettyjo.
#   It will download the latest data and save at 'scratch'

###################### STEP 2 ######################
# Parse the cube file into individual ASCII files.
# Calculate the 1-D PACS spectrum with a given aperture size
skip = True
for o in obsid:
    if o[3] == '0':
        continue
    if o[1] == '0':
        continue
    if o[0] == 'IRS46':
        # for skipping processed objects
        skip = False
        continue
    if skip:
        continue
    # load aperture from SPIRE SECT reduction
    if os.path.exists(outdir+str(o[0])+'/spire/data/'+str(o[0])+'_spire_phot.txt'):
        spire_phot = ascii.read(outdir+str(o[0])+'/spire/data/'+str(o[0])+'_spire_phot.txt', data_start=4)
        aper_size = spire_phot['aperture(arcsec)'][spire_phot['wavelength(um)'] == spire_phot['wavelength(um)'].min()][0]
    else:
        aper_size = 31.8
    aper_size_fitted = cdfPacs1d(o[1:3], pacsdatadir, outdir+o[0], o[0])
    print o[0], aper_size_fitted

# need to modify this part for it to automatically figure out the aperture size that will result in a well-matched spectrum

###################### STEP 3 ######################
# Line fitting on both cube and 1-D products

idl('.r /home/bettyjo/yaolun/programs/line_fitting/extract_pacs.pro')
idl('.r /home/bettyjo/yaolun/programs/gauss.pro')

for o in obsid:
    if o[3] == '0':
        continue
    if o[1] == '0':
        continue

    # cube
    for ip in range(1,26):
        idl.pro('extract_pacs', indir=outdir+str(o[0])+'/pacs/data/cube/', filename=str(o[0])+'_pacs_pixel'+str(ip)+'_hsa',
                outdir=outdir+str(o[0])+'/pacs/advanced_products/cube/', plotdir=outdir+str(o[0])+'/pacs/advanced_products/cube/plots/',
                noiselevel=3, localbaseline=10, global_noise=20, fixed_width=1, opt_width=1, continuum=1, flat=1, object=str(o[0]),
                current_pix=str(ip), double_gauss=1, print_all=outdir+reduction_name+'_lines')
    # 1-D
    # read in RA/Dec
    radec_pacs = ascii.read(outdir+str(o[0])+'/pacs/data/cube/'+str(o[0])+'_pacs_pixel13_hsa_coord.txt')
    idl.pro('extract_pacs', indir=outdir+str(o[0])+'/pacs/data/cube/', filename=str(o[0])+'_pacs_weighted',
            outdir=outdir+str(o[0])+'/pacs/advanced_products/', plotdir=outdir+str(o[0])+'/pacs/advanced_products/plots/',
            noiselevel=3, localbaseline=10, global_noise=20, fixed_width=1, opt_width=1, continuum=1, flat=1, object=str(o[0]),
            ra=np.mean(radec_pacs['RA(deg)']), dec=np.mean(radec_pacs['Dec(deg)']), current_pix=str(ip),
            double_gauss=1, print_all=outdir+reduction_name+'_lines')

###################### STEP 4 ######################
# extract photometry with the aperture size specified or derived in above
# [HIPE] execution
#   - "pacs_phot_cdf.py" in hipe_script folder

# the default aperture size now is 31.8".
