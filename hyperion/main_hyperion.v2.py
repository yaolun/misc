# Parse the command-line options
import argparse

# get command-line arguments (ver 2.0)
parser = argparse.ArgumentParser(description='Options for manually control the modeling')
# parser.add_argument('--norun', action='store_false',
                    # help='disable the Hyperion run and only run the extraction part')
parser.add_argument('--norecord', action='store_true',
                    help='do not keep the record of this run')
parser.add_argument('--mono_wave',
                    help='for monochromatic Hyperion run (i.e. produce images at given wavelengths)')
parser.add_argument('--extract_only', action='store_true',
                    help='only run the extraction part and skip the Hyperion run')
parser.add_argument('--temp', action='store_true',
                    help='extract the temperature profile')
parser.add_argument('--skip_regular', action='store_true',
                    help='skip both image and SED extractions')
parser.add_argument('--core_num', default='20',
                    help='numbers of core used for Hyperion (default: 20)')
parser.add_argument('--mc_photons', default='1e6', type=float,
                    help='numbers of photons for calculating specfic energy (default: 1e6)')
parser.add_argument('--im_photons', default='1e6', type=float,
                    help='numbers of photons for the SED/image calculation (default: 1e6)')
parser.add_argument('--test', action='store_true')
parser.add_argument('--fast_plot', action='store_true',
                    help='do not make plots of gas density')
parser.add_argument('--low_res', action='store_true',
                    help='reduce the number of cell in r-axis (from 300 to 100)')
parser.add_argument('--prefix', default='',
                    help='prefix for the input table and output model list')
parser.add_argument('--ellipsoid', action='store_true',
                    help='use ellipsoid cavity instead of a power-law cavity (not fully tested in v2)')
parser.add_argument('--no_azimuthal', action='store_true',
                    help='disable the extraction of the radial profile of brightness')

args = vars(parser.parse_args())

# option for high resolution r-grid
# the angular range at which the azimuthal averaged radial intensity will perform.
args['rrange'] = [10, 200]
# the range of wavelength.  [wav_min, wav_max, wav_num]
args['wav_range'] = (1.0, 25., 1400)

def run_hyperion(**args):

    import numpy as np
    import os
    import sys
    from subprocess import Popen, call
    from pprint import pprint
    from setup_model import setup_model
    from input_reader import input_reader_table
    from extract_model import extract_hyperion
    from temp_hyperion import temp_hyperion
    from hyperion_image import hyperion_image
    from azimuthal_avg_radial_intensity import azimuthal_avg_radial_intensity
    import time

    # Default parameter
    fix_params = None

    print('wavelength ranging from ', wav_range[0], ' to ', wav_range[1],
          ' with ', wav_num, ' channels')
    # require additional input for monochromatic radiative transfer
    if mono_wave != None:
        image_only = True
        azimuthal = False
        print('Monochromatic RT now force "image_only" simulations.')

        # Here are some pre-set wavelengths
        if mono_wave == 'NIR':
            mono_wave = [1.25, 1.53]
        elif mono_wave == 'IRAC':
            mono_wave = [3.6, 4.5, 5.8, 8.0]
        elif mono_wave == 'MIPS':
            mono_wave = [24., 70., 160.]
        elif mono_wave == 'PACS':
            mono_wave = [70., 100., 160.]
        elif mono_wave == 'SPIRE':
            mono_wave = [250., 350., 500.]
        elif mono_wave == 'H':
            mono_wave = [1.6]
        else:
            # if none of those presets is specified, the program will prompt for input
            mono_wave = input('What are the bands for monochromatic RT? (separated by space)')
            mono_wave = [float(i) for i in mono_wave.split()]

        print('Simulations will be performed at the following wavelengths: ', mono_wave)

        unit = 'MJy\,sr^{-1}'
        print('Image unit is set to be MJy/sr')

    # path setting version 1.1
    # The path file "run_hyperion_path.txt" has to be placed at the same directory as main_hyperion.py
    home = os.path.expanduser('~')
    path_list = np.genfromtxt('run_hyperion_path.txt', dtype=str).T
    dict_path = {}
    for name, val in zip(path_list[0],path_list[1]):
        dict_path[name] = val
    obj = dict_path['object']  # object name
    dstar = float(dict_path['dstar'].data)  # distance in parsec

    # TODO: delete
    # if 'ext_source' not in path_list[0]:
    #     ext_source_path = None
    # else:
    #     ext_source_path = home + dict_path['ext_source']

    print('Current path setting --')
    pprint(dict_path)
    #
    # Read in aperture info - under obs_dir with filename "aperture.txt"
    wl_aper, aper_arcsec = np.genfromtxt(home+dict_path['obs_dir']+'aperture.txt', skip_header=1, dtype=float).T
    aperture = {'wave': wl_aper, 'aperture': aper_arcsec}

    if test:
        print('testing mode...')
        prefix = 'test'

    if not 'params_table' in globals():
        print('using the default path to input table')
        params_table = home+dict_path['input_table']+'input_table'+prefix+'.txt'
        # check existence of parameter table
        if not os.path.exists(params_table):
            raise RuntimeError('parameter table path invalid.')
        outdir = home + dict_path['outdir']+prefix+'/'

    params = input_reader_table(params_table)

    # Get current model number
    if not os.path.exists(outdir+'model_list.txt'):
        last_model_num = '0'
    else:
        last = open(outdir+'model_list.txt','r').readlines()[-1]
        last_model_num = (last.split('M_env_dot')[0]).split('Model')[1].split()[0]

    # setup the model
    if not extract_only:
        # determine starting model number for the current batch
        model_num = str(int(last_model_num)+1)
        #
        # iterate through the sets of parameters
        for i in range(0, len(params)):
            # get the current set of parameters
            params_dict = params[i]
            if 'cav_power' in params_dict.keys():
                power = params_dict['cav_power']
            else:
                power = 2.0
                print('there is no cav_power column in input_table, setting it to 2')
            # create model directory
            if not os.path.exists(outdir+'model'+str(int(model_num)+i)+'/'):
                os.makedirs(outdir+'model'+str(int(model_num)+i)+'/')
            outdir_dum = outdir+'model'+str(int(model_num)+i)+'/'

            # print out some information about the current calculating model
            print('Model'+str(int(model_num)+i))
            pprint(params_dict)
            # calculate the initial dust profile
            # option to fix some parameter
            # fix_params = {'R_min': 0.14}
            m = setup_model(outdir_dum,outdir, 'model'+str(int(model_num)+i),
                            params_dict, home+dict_path['dust_file'], wav_range, aperture,
                            plot=True, fast_plot=fast_plot, idl=True, norecord=norecord,
                            mono_wave=mono_wave, fix_params=fix_params,
                            low_res=low_res, power=power, im_photons=im_photons,
                            mc_photons=mc_photons, ellipsoid=ellipsoid, dstar=dstar,
                            TSC_dir=home+dict_path['TSC_dir'],
                            IDL_path=dict_path['IDL_path'], image_only=image_only)

            # Run hyperion
            print('Running with Hyperion')
            hyp_foo = open(outdir_dum+'hyperion.log','w')
            hyp_err = open(outdir_dum+'hyperion.err','w')
            run = Popen(['mpirun', '-n', str(core_num),
                         'hyperion_sph_mpi', '-f',
                         outdir_dum+'model'+str(int(model_num)+i)+'.rtin',
                         outdir_dum+'model'+str(int(model_num)+i)+'.rtout'],
                         stdout=hyp_foo, stderr=hyp_err)
            run.communicate()
            # Extract the results
            # the indir here is the dir that contains the observed spectra.
            print('Seems finish, lets check out the results')
            # extract SEDs and images
            if mono_wave == None:
                extract_hyperion(outdir_dum+'model'+str(int(model_num)+i)+'.rtout',
                                 indir=home+dict_path['obs_dir'], outdir=outdir_dum,
                                 aperture=aperture, filter_func=True, obj=obj, dstar=dstar)
            else:
                for w in mono_wave:
                    hyperion_image(outdir_dum+'model'+str(int(model_num)+i)+'.rtout', w, outdir_dum,
                        'model'+str(int(model_num)+i),dstar=dstar, unit=unit)
            # extract the temperature profile
            if temp:
                temp_hyperion(outdir_dum+'model'+str(int(model_num)+i)+'.rtout', outdir=outdir_dum)

            # extract the radial profile
            if not no_azimuthal:
                # the wavelength for plotting azimuthal-averaged radial intensity
                azi_wave = dict_path['azi_wave'].split('_')

                # whether a reference image file is specified
                if 'img_name' not in dict_path.keys():
                    obs_azi = None
                else:
                    imgpath = home+dict_path['obs_dir']+dict_path['img_name']+'.fits'
                    source_center = dict_path['source_ra']+' '+dict_path['source_dec']
                    obs_azi = [{'imgpath': imgpath, 'source_center': source_center}]

                aper_reduced = list(set(aperture['aperture']))

                for azi_w in azi_wave:
                    azimuthal_avg_radial_intensity(float(azi_w),
                            outdir_dum+'model'+str(int(model_num)+i)+'.rtout',
                            outdir_dum+'model'+str(int(model_num)+i), dstar,
                            annulus_width=10, group=len(aper_reduced)+1,
                            obs=obs_azi, rrange=rrange)
    else:
        print 'You have entered the extract-only mode...'
        num_min = input('What is the number of the first model?')
        num_max = input('What is the number of the last model?')
        num_min = int(num_min)
        num_max = int(num_max)
        # iterate the models between the min and max model numbers given
        for i in range(num_min, num_max+1):
            if not os.path.exists(outdir+'model'+str(i)+'/'):
                os.makedirs(outdir+'model'+str(i)+'/')
            outdir_dum = outdir+'model'+str(i)+'/'
            # print out some information about the current calculating model
            print 'Extracting Model'+str(i)
            # Extract the results
            # the indir here is the dir that contains the observed spectra.
            if not skip_regular:
                if mono_wave != None:
                    extract_hyperion(outdir_dum+'model'+str(i)+'.rtout',
                                     indir=home+dict_path['obs_dir'],
                                     outdir=outdir_dum, aperture=aperture,
                                     filter_func=True, obj=obj, dstar=dstar)
                else:
                    for w in mono_wave:
                        hyperion_image(outdir_dum+'model'+str(i)+'.rtout', w,
                            outdir_dum, 'model'+str(i),dstar=dstar, unit=unit)
            if temp:
                temp_hyperion(outdir_dum+'model'+str(i)+'.rtout',outdir=outdir_dum)

            if not no_azimuthal:
                # the wavelength for plotting azimuthal-averaged radial intensity
                azi_wave = dict_path['azi_wave'].split('_')

                # whether a reference image file is specified
                if 'img_name' not in dict_path.keys():
                    obs_azi = None
                else:
                    imgpath = home+dict_path['obs_dir']+dict_path['img_name']+'.fits'
                    source_center = dict_path['source_ra']+' '+dict_path['source_dec']
                    obs_azi = [{'imgpath': imgpath, 'source_center': source_center}]

                aper_reduced = list(set(aperture['aperture']))

                for azi_w in azi_wave:
                    azimuthal_avg_radial_intensity(float(azi_w),
                            outdir_dum+'model'+str(i)+'.rtout',
                            outdir_dum+'model'+str(i), dstar,
                            annulus_width=10, group=len(aper_reduced)+1,
                            obs=obs_azi, rrange=rrange)
