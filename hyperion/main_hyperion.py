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

# option for high resolution r-grid
# !!!
low_res = True

# Default setting
run = True
record = True
mono = False
mono_wave = None
control = False
extract_only = False
temp = False
alma=False
core_num = 20
better_im = False
chi2 = False
test = False
ellipsoid = False
fast_plot = False
image_only=False
azimuthal=True
fix_params = {}

# Get command-line arguments
if 'norun' in sys.argv:
    run = False
if 'norecord' in sys.argv:
    record = False
if 'mono' in sys.argv:
    mono = True
    image_only = True
    print 'Monochromatic RT now force "image_only" simulations.'
    print 'Need to go the code for more options.'
if 'control' in sys.argv:
    control = True
if 'extract_only' in sys.argv:
    extract_only = True
if 'temp' in sys.argv:
    temp = True
if 'alma' in sys.argv:
    alma = True
if '18' in sys.argv:
    core_num = 18
if 'better_im' in sys.argv:
    better_im = True
if 'chi2' in sys.argv:
    chi2 = True
if 'test' in sys.argv:
    test = True
if 'ellipsoid' in sys.argv:
    ellipsoid = True
if 'fast_plot' in sys.argv:
    fast_plot = True

print 'Setting - run: %s, record: %s, mono: %s' % (run,record,mono)

# require additional input for monochromatic radiative transfer
if mono:
    mono_wave = raw_input('What are the bands for monochromatic RT?')
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

    print 'Simulations will be performed at the following wavelengths: ', mono_wave
# path setting version 1.1
# The path file "run_hyperion_path.txt" has to be placed at the same directory as main_hyperion.py
home = os.path.expanduser('~')
path_list = np.genfromtxt('run_hyperion_path.txt', dtype=str).T
dict_path = {}
for name, val in zip(path_list[0],path_list[1]):
    dict_path[name] = val
obj = dict_path['object']
dstar = float(dict_path['dstar'].data)
print 'Current path setting --'
pprint(dict_path)
#
# Read in aperture info - under obs_dir with filename "aperture.txt"
wl_aper, aper_arcsec = np.genfromtxt(home+dict_path['obs_dir']+'aperture.txt', skip_header=1, dtype=float).T
aperture = {'wave': wl_aper, 'aperture': aper_arcsec}

if control == True:
    print 'Running the controlled grids for paper...'
    params_table = home + dict_path['input_table']+'input_table_control.txt'
    outdir = home + dict_path['outdir']+'controlled/'
if alma == True:
    print 'Running for ALMA proposal...'
    params_table = home + dict_path['input_table']+'input_table_alma.txt'
    outdir = home + dict_path['outdir']+'alma/'
if chi2 == True:
    print 'Running for chi2 grid...'
    params_table = home + dict_path['input_table']+'input_table_chi2.txt'
    outdir = home + dict_path['outdir']+'chi2_grid/'
if test == True:
    print 'testing mode...'
    params_table = home + dict_path['input_table']+'test_input.txt'
    outdir = home + dict_path['outdir']+'test/'
if ellipsoid == True:
    print 'Running with ellipsoid cavities...'
    params_table = home + dict_path['input_table']+'input_table_ellipsoid.txt'
    outdir = home + dict_path['outdir']+'ellipsoid/'
if not 'params_table' in globals():
    print 'using the default path to input table'
    params_table = home+dict_path['input_table']+'input_table.txt'
    outdir = home + dict_path['outdir']

params = input_reader_table(params_table)

# Get current model number
if not os.path.exists(outdir+'model_list.txt'):
    last_model_num = '0'
else:
    foo = open(outdir+'model_list.txt','r')
    for line in foo.readlines():
        pass
    last = line
    foo.close()
    last_model_num = (last.split('M_env_dot')[0]).split('Model')[1].split()[0]

if extract_only == False:
    model_num = str(int(last_model_num)+1)
    #
    for i in range(0, len(params)):
        params_dict = params[i]
        if 'cav_power' in params_dict.keys():
            power = params_dict['cav_power']
        else:
            power = 2.0
        if not os.path.exists(outdir+'model'+str(int(model_num)+i)+'/'):
            os.makedirs(outdir+'model'+str(int(model_num)+i)+'/')
        outdir_dum = outdir+'model'+str(int(model_num)+i)+'/'
        # print out some information about the current calculating model
        print 'Model'+str(int(model_num)+i)
        pprint(params_dict)
        # calculate the initial dust profile
        # option to fix some parameter
        # fix_params = {'R_min': 0.14}
        m = setup_model(outdir_dum,outdir,'model'+str(int(model_num)+i),params_dict,
                        home+dict_path['dust_file'],plot=True,fast_plot=fast_plot,
                        idl=True,record=record,mono=mono,mono_wave=mono_wave,
                        aperture=aperture,fix_params=fix_params, low_res=low_res,
                        power=power,better_im=better_im,ellipsoid=ellipsoid,
                        dstar=dstar,TSC_dir=home+dict_path['TSC_dir'],
                        IDL_path=dict_path['IDL_path'], image_only=image_only)
        if run == False:
            print 'Hyperion run is skipped. Make sure you have run this model before'
        else:
            # Run hyperion
            print 'Running with Hyperion'
            hyp_foo = open(outdir_dum+'hyperion.log','w')
            hyp_err = open(outdir_dum+'hyperion.err','w')
            run = Popen(['mpirun','-n',str(core_num),'hyperion_sph_mpi','-f',outdir_dum+'model'+str(int(model_num)+i)+'.rtin',outdir_dum+'model'+str(int(model_num)+i)+'.rtout'], stdout=hyp_foo, stderr=hyp_err)
            run.communicate()
        # Extract the results
        # the indir here is the dir that contains the observed spectra.
        print 'Seems finish, lets check out the results'
        if not mono:
            extract_hyperion(outdir_dum+'model'+str(int(model_num)+i)+'.rtout',
                             indir=home+dict_path['obs_dir'],outdir=outdir_dum,
                             aperture=aperture,filter_func=True,obj=obj,dstar=dstar)
        else:
            if type(mono_wave) is str:
                hyperion_image(outdir_dum+'model'+str(int(model_num)+i)+'.rtout',
                        float(mono_wave), outdir_dum, 'model'+str(int(model_num)+i),dstar=dstar)
            else:
                for w in mono_wave:
                    hyperion_image(outdir_dum+'model'+str(int(model_num)+i)+'.rtout', w, outdir_dum,
                        'model'+str(int(model_num)+i),dstar=dstar)
        if temp:
            temp_hyperion(outdir_dum+'model'+str(int(model_num)+i)+'.rtout',outdir=outdir_dum)

        if azimuthal:
            imgpath = home+dict_path['obs_dir']+'hspireplw1342226633_20pxmp_1431669349619.fits'
            source_center = dict_path['source_ra']+' '+dict_path['source_dec']
            aper_reduced = list(set(aperture['aperture']))
            azimuthal_avg_radial_intensity(500.0, imgpath, source_center,
                    outdir_dum+'model'+str(int(model_num)+i)+'.rtout',
                    'model'+str(int(model_num)+i),
                    annulus_width=10, group=len(aper_reduced)+1, dstar=dstar)
else:
    print 'You have entered the extract-only mode...'
    num_min = raw_input('What is the number of the first model?')
    num_max = raw_input('What is the number of the last model?')
    num_min = int(num_min)
    num_max = int(num_max)
    for i in range(num_min, num_max+1):
        if not os.path.exists(outdir+'model'+str(i)+'/'):
            os.makedirs(outdir+'model'+str(i)+'/')
        outdir_dum = outdir+'model'+str(i)+'/'
        # print out some information about the current calculating model
        print 'Extracting Model'+str(i)
        # Extract the results
        # the indir here is the dir that contains the observed spectra.
        if not mono:
            extract_hyperion(outdir_dum+'model'+str(i)+'.rtout',indir=home+dict_path['obs_dir'],
                             outdir=outdir_dum,aperture=aperture,
                             filter_func=True,obj=obj,dstar=dstar)
        else:
            if type(mono_wave) is str:
                hyperion_image(outdir_dum+'model'+str(i)+'.rtout',
                        float(mono_wave), outdir_dum, 'model'+str(i),dstar=dstar)
            else:
                for w in mono_wave:
                    hyperion_image(outdir_dum+'model'+str(i)+'.rtout', w,
                        outdir_dum, 'model'+str(i),dstar=dstar)
        if temp:
            temp_hyperion(outdir_dum+'model'+str(i)+'.rtout',outdir=outdir_dum)

        if azimuthal:
            imgpath = home+dict_path['obs_dir']+'hspireplw1342226633_20pxmp_1431669349619.fits'
            source_center = dict_path['source_ra']+' '+dict_path['source_dec']
            aper_reduced = list(set(aperture['aperture']))
            azimuthal_avg_radial_intensity(500.0, imgpath, source_center,
                    outdir_dum+'model'+str(i)+'.rtout',
                    'model'+str(i),
                    annulus_width=10, group=len(aper_reduced)+1, dstar=dstar)
