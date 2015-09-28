def record_hyperion(dict_params, outdir):
	import numpy as np
	import os
	from pprint import pprint

	# the column head of the output record table
	# If this is changed, make sure change the output tuple at the end of this function as well.
	colhead = ('Model#','age','Cs','Omega0','tstar','R_env_max','theta_cav','rho_cav_center','rho_cav_edge','rstar','M_disk_dot',\
			   'M_disk','beta','h100','rho_cav','d_sub','M_env_dot','R_inf','R_cen','mstar','percentile','absolute','relative','view_angle','total_mass')

	if not os.path.exists(outdir+'model_list.txt'):
		foo = open(outdir+'model_list.txt','w')
		foo.write('%10s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s \n' % colhead) 
		foo.close()

	foo = open(outdir+'model_list.txt','r')
	for line in foo.readlines():
		pass
	last = line
	foo.close()
	last_model_num = (last.split('M_env_dot')[0]).split('Model')[1].split()[0]
	if last_model_num == '#':
		model_num = '1'
	else:
		model_num = str(int(last_model_num)+1)
	foo = open(outdir+'model_list.txt','a')
	output = 'Model'+model_num
	output = (output,dict_params['age'],dict_params['Cs'],dict_params['Omega0'],dict_params['tstar'],\
			  dict_params['R_env_max'],dict_params['theta_cav'],dict_params['rho_cav_center'],dict_params['rho_cav_edge'],\
			  dict_params['rstar'],dict_params['M_disk_dot'],dict_params['M_disk'],dict_params['beta'],dict_params['h100'],\
			  dict_params['rho_cav'],dict_params['d_sub'],dict_params['M_env_dot'],dict_params['R_inf'],dict_params['R_cen'],\
			  dict_params['mstar'],dict_params['percentile'],dict_params['absolute'],dict_params['relative'],\
			  dict_params['view_angle'],dict_params['total_mass'])
	foo.write('%10s  %14e  %14e  %14e  %14e  %14e  %14e  %14e  %14e  %14e  %14e  %14e  %14e  %14e  %14e  %14e  %14e  %14e  %14e  %14e  %14e  %14e  %14e  %14e  %14e \n' % output)

# from input_reader import input_reader
# filename = '/Users/yaolun/programs/misc/hyperion/tsc_params.dat'
# dict_params = input_reader(filename)
# outdir = '/Users/yaolun/programs/misc/hyperion/'
# record_hyperion(dict_params,outdir)