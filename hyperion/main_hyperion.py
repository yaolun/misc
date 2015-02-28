import numpy as np
import os
from subprocess import Popen
from pprint import pprint
from setup_model import setup_model
from input_reader import input_reader_table

# path setting
home = os.path.expanduser('~')
outdir = home + '/hyperion/bhr71/'
dust_file = home + '/programs/misc/dustkappa_oh5_extended.inp'
params_table = home + '/programs/misc/hyperion/input_table.txt'
obs_dir = home + '/radmc_simulation/bhr71/observations/'

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

for i in range(0, len(params)):
	params_dict = params[i]
	if not os.path.exists(outdir+'model'+str(int(model_num)+i)+'/'):
		os.makedirs(outdir+'model'+str(int(model_num)+i)+'/')
	outdir_dum = outdir+'model'+str(int(model_num)+i)+'/'
	# print out some information about the current calculating model
	print 'Model'+str(int(model_num)+i)
	pprint(params_dict)
	# calculate the initial dust profile
	m = setup_model(outdir_dum,outdir,'model'+str(int(model_num)+i),params_dict,dust_file,plot=True,idl=True)
	# Run hyperion
	print 'Running with Hyperion'
	# m.run(outdir_dum+'model'+str(int(model_num)+i)+'.rtout', mpi=True, n_processes=22)
	run = Popen('mpirun -n 22 hyperion_sph_mpi -f '+outdir_dum+'model'+str(int(model_num)+i)+'.rtin '+outdir_dum+'model'+str(int(model_num)+i)+'.rtout')
	run.communicate()
	# Extract the results
	# the indir here is the dir that contains the observed spectra.
	print 'Seems finish, lets check out the results'
	extract_hyperion(outdir_dum+'model'+str(int(model_num)+i)+'.rtout',indir=obs_dir,outdir=outdir_dum)
