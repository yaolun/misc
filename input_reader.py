def input_reader(filename, dict_params=None, default=True):
	"""
	Read the text file version of parameters into a dictionary
	- minor: range one parameter for given start, end, and number of separations
	"""

	import numpy as np
	import astropy.constants as const
	import pprint as pp

	# Constants setup
	c         = const.c.cgs.value
	AU        = const.au.cgs.value                         # Astronomical Unit       [cm]
	pc        = const.pc.cgs.value                         # Parsec                  [cm]
	MS        = const.M_sun.cgs.value                      # Solar mass              [g]
	LS        = const.L_sun.cgs.value                      # Solar luminosity        [erg/s]
	RS        = const.R_sun.cgs.value                      # Solar radius            [cm]
	G         = const.G.cgs.value                          # Gravitational constant  [cm^3/g/s^2]
	yr        = 60*60*24*365.                              # Years in seconds        [s]
	PI        = np.pi                                      # PI constant
	sigma     = const.sigma_sb.cgs.value                   # Stefan-Boltzmann constant 
	mh        = const.m_p.cgs.value + const.m_e.cgs.value  # Mass of Hydrogen atom   [g]
	# Maybe don't need this
	u = dict([('c',c), ('AU',AU), ('pc',pc), ('MS',MS), ('LS',LS), ('RS',RS), ('G',G), ('yr',yr), ('sigma',sigma), ('PI', np.pi), ('mh',mh)])

	# if there is no dict_params provided, then I must get the parameters from somewhere
	if dict_params == None:
		params = np.genfromtxt(filename, dtype=None)
		dict_params = {}
		for var, val in params:
			dict_params[var] = val
		# Use the default parameters for mostly fixed parameter (e.g. disk)
		if default == True:
			dict_params['M_disk'] = 5e-1
			dict_params['M_disk_dot'] = 6.5e-7
			dict_params['beta'] = 1.093
			dict_params['h100'] = 8.123
			dict_params['rho_cav'] = 1e-21
		pp.pprint(dict_params)
	else:
		# the dict_params can be provided by reading in .prev_params
		print 'The parameters are provided by user'
		pp.pprint(dict_params)

	return dict_params
filename = '/Users/yaolun/bhr71/radmc3d_params/tsc_params.dat'
input_reader(filename)