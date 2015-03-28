import numpy as np
import os
import sys
from subprocess import Popen
from pprint import pprint
from setup_model import setup_model
from input_reader import input_reader_table
from extract_model import extract_hyperion
from temp_hyperion import temp_hyperion

# Default setting
run = True
record = True
mono = False

# Get command-line arguments
if 'norun' in sys.argv:
    run = False
if 'norecord' in sys.argv:
    record = False
if 'mono' in sys.argv:
    mono = True

print 'Setting - run: %s, record: %s, mono: %s' % (run,record,mono)

# path setting
home = os.path.expanduser('~')
outdir = home + '/hyperion/bhr71/'
dust_file = home + '/programs/misc/dustkappa_oh5_extended.inp'
params_table = home + '/programs/misc/hyperion/input_table.txt'
obs_dir = home + '/radmc_simulation/bhr71/observations/'
# temp fix for the broken /opt/local/ of bettyjo
# outdir = home+'/test/hyperion/'
# obs_dir = home+'/bhr71/obs_for_radmc/'

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

model_num = str(int(last_model_num)+1)
# 
model_num = 47
#
# for i in range(0, len(params)):
for i in range(0,9):
    # params_dict = params[i]
    if not os.path.exists(outdir+'model'+str(int(model_num)+i)+'/'):
        os.makedirs(outdir+'model'+str(int(model_num)+i)+'/')
    outdir_dum = outdir+'model'+str(int(model_num)+i)+'/'
    # print out some information about the current calculating model
    print 'Model'+str(int(model_num)+i)
    # pprint(params_dict)
    # calculate the initial dust profile
    wl_aper = [3.6, 4.5, 5.8, 8.0, 10, 16, 20, 24, 35, 70, 100, 160, 250, 350, 500, 850]
    # option to fix some parameter
    fix_params = {'R_min': 0.14}
    # m = setup_model(outdir_dum,outdir,'model'+str(int(model_num)+i),params_dict,dust_file,plot=True,idl=True,record=record,mono=mono,wl_aper=wl_aper,fix_params=fix_params)
    if run == False:
        print 'Hyperion run is skipped. Make sure you have run this model before'
    else:
        # Run hyperion
        print 'Running with Hyperion'
        hyp_foo = open(outdir_dum+'hyperion.log','w')
        hyp_err = open(outdir_dum+'hyperion.err','w')
        run = Popen(['mpirun','-n','20','hyperion_sph_mpi','-f',outdir_dum+'model'+str(int(model_num)+i)+'.rtin',outdir_dum+'model'+str(int(model_num)+i)+'.rtout'], stdout=hyp_foo, stderr=hyp_err)
        run.communicate()
    # Extract the results
    # the indir here is the dir that contains the observed spectra.
    print 'Seems finish, lets check out the results'
    extract_hyperion(outdir_dum+'model'+str(int(model_num)+i)+'.rtout',indir=obs_dir,outdir=outdir_dum,wl_aper=wl_aper)
    temp_hyperion(outdir_dum+'model'+str(int(model_num)+i)+'.rtout',outdir=outdir_dum)