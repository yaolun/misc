def grid_create(list_params, outdir):

	import numpy as np
	import itertools as iter
	import copy
	from pprint import pprint
 	
	# default parameter setup
	ref = {'age': 3e4, 'Cs': 0.38, 'Omega0': 1e-13, 'tstar': 5267., 'R_env_max': 1e4, 'theta_cav': 20., 'rho_cav_center': 5e-19,\
		   'rho_cav_edge': 40., 'rstar': 5., 'M_disk_dot': 6.5e-7, 'M_disk': 0.01, 'beta': 1.093, 'rho_cav':1e-21, 'h100': 8.123,\
		   'percentile': 95., 'absolute': 2.0, 'relative': 1.02, 'view_angle': 40., 'cav_power': 2.0}
	# This is the right order
	colhead = ('age','Cs','Omega0','tstar','R_env_max','theta_cav','rho_cav_center','rho_cav_edge','rstar','M_disk_dot',\
				'M_disk','beta','h100','rho_cav','percentile','absolute','relative','view_angle','cav_power')
	# cartiesian product of lists
 	product = [x for x in apply(iter.product, list_params.values())]
 	# iterate through models and print as an input table file

 	# print the column names
 	foo = open(outdir+'input_table_chi2.txt', 'w')
 	foo.write('%14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s \n' % colhead) 

 	for mod in product:
 		params_dum = copy.copy(ref)
 		for col in list_params.keys():
 			params_dum[col] = mod[list_params.keys().index(col)]
 		# pprint(params_dum)
 		output = (params_dum['age'],params_dum['Cs'],params_dum['Omega0'],params_dum['tstar'],\
			  params_dum['R_env_max'],params_dum['theta_cav'],params_dum['rho_cav_center'],params_dum['rho_cav_edge'],\
			  params_dum['rstar'],params_dum['M_disk_dot'],params_dum['M_disk'],params_dum['beta'],params_dum['h100'],\
			  params_dum['rho_cav'],params_dum['percentile'],params_dum['absolute'],params_dum['relative'],\
			  params_dum['view_angle'],params_dum['cav_power'])
		foo.write('%14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s \n' % output)
	foo.close()
	return list_params.keys

import os
home = os.path.expanduser('~')

list_params = {'age': [1e4, 2e4, 3e4, 4e4, 5e4],\
			   'rho_cav_center': [5e-20, 1e-19, 5e-19, 1e-18, 5e-18],\
			   'rho_cav_edge': [20,30,40,50,60]}

# outdir = '/Users/yaolun/test/'
outdir = home + '/programs/misc/hyperion/'
grid_create(list_params,outdir)