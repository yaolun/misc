def inf_rad_check(gridi, gridc, rhoenv, R_inf):
	import numpy as np
	import astropy.constants as const

	# constant setup
	AU = const.au.cgs.value

	(ri, thetai, phii) = gridi
	(rc, thetac, phic) = gridc

	if R_inf > max(ri):
		print 'Infall radius is greater than the cloud size'