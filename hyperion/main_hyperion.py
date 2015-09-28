import numpy as np
import os
import sys
from subprocess import Popen, call
from pprint import pprint
from setup_model import setup_model
from input_reader import input_reader_table
from extract_model import extract_hyperion
from temp_hyperion import temp_hyperion
import time

# Default setting
run = True
record = True
mono = False
control = False
extract_only = False
temp = True
alma=False
core_num = 20
better_im = False
chi2 = False
test = False
ellipsoid = False
fix_params = {}

# Get command-line arguments
if 'norun' in sys.argv:
    run = False
if 'norecord' in sys.argv:
    record = False
if 'mono' in sys.argv:
    mono = True
if 'control' in sys.argv:
    control = True
if 'extract_only' in sys.argv:
    extract_only = True
if 'no_temp' in sys.argv:
    temp = False
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


print 'Setting - run: %s, record: %s, mono: %s' % (run,record,mono)

# path setting
# home = os.path.expanduser('~')
# outdir = home + '/hyperion/bhr71/'
# dust_file = home + '/programs/misc/oh5_hyperion.txt'
# params_table = home + '/programs/misc/hyperion/input_table.txt'
# obs_dir = home + '/radmc_simulation/bhr71/cycle1/observations/'
#
# path setting version 1.1
# The path file "run_hyperion_path.txt" has to be placed at the same directory as main_hyperion.py
home = os.path.expanduser('~')
path_list = np.genfromtxt('run_hyperion_path.txt', dtype=str).T
dict_path = {}
for name, val in zip(path_list[0],path_list[1]):
    dict_path[name] = val
print 'Current path setting --'
pprint(dict_path)
#
# Read in aperture info - under obs_dir with filename "aperture.txt"
wl_aper, aper_arcsec = np.genfromtxt(home+dict_path['obs_dir']+'aperture.txt', skip_header=1, dtype=float).T
aperture = {'wave': wl_aper, 'aperture': aper_arcsec}
# wl_aper = [3.6, 4.5, 5.8, 8.0, 8.5, 9, 9.7, 10, 10.5, 11, 16, 20, 24, 35, 70, 100, 160, 250, 350, 500, 850]
#
# if control == True:
#     print 'Running the controlled grids for paper...'
#     params_table = home + '/programs/misc/hyperion/input_table_control.txt'
#     outdir = home + '/hyperion/bhr71/controlled/'
# if alma == True:
#     print 'Running for ALMA proposal...'
#     params_table = home + '/programs/misc/hyperion/input_table_alma.txt'
#     outdir = home + '/hyperion/bhr71/alma/'
# if chi2 == True:
#     print 'Running for chi2 grid...'
#     params_table = home + '/programs/misc/hyperion/input_table_chi2.txt'
#     outdir = home + '/hyperion/bhr71/chi2_grid/'
# if test == True:
#     print 'testing mode...'
#     params_table = home + '/programs/misc/hyperion/test_input.txt'
#     outdir = home + '/hyperion/bhr71/test/'
# if ellipsoid == True:
#     print 'Running with ellipsoid cavities...'
#     params_table = home + '/programs/misc/hyperion/input_table_ellipsoid.txt'
#     outdir = home + '/hyperion/bhr71/ellipsoid/'
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
if params_table in locals() == False:
    params_table = home+dict_path['input_table']+'input_table.txt'

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
        fix_params = {'R_min': 0.14}
        m = setup_model(outdir_dum,outdir,'model'+str(int(model_num)+i),params_dict,home+dict_path['dust_file'],\
            plot=True,idl=True,record=record,mono=mono,aperture=aperture,fix_params=fix_params,alma=alma,\
            power=power,better_im=better_im,ellipsoid=ellipsoid,TSC_dir=home+dict_path['TSC_dir'],\
            IDL_path=dict_path['IDL_path'])
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
        extract_hyperion(outdir_dum+'model'+str(int(model_num)+i)+'.rtout',indir=dict_path['obs_dir'],outdir=outdir_dum,aperture=aperture,filter_func=True)
        temp_hyperion(outdir_dum+'model'+str(int(model_num)+i)+'.rtout',outdir=outdir_dum)
else:
    print 'You are entering the extract-only mode...'
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
        extract_hyperion(outdir_dum+'model'+str(i)+'.rtout',indir=dict_path['obs_dir'],outdir=outdir_dum,aperture=aperture,filter_func=True)
        if temp == True:
            temp_hyperion(outdir_dum+'model'+str(i)+'.rtout',outdir=outdir_dum)
