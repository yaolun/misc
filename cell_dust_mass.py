def cell_dust_mass(number,ri,thetai,phii, rho):
	"""
	number: [irc, ithetac, iphic]. This indice are the indice for the grid of cell center.
	[ri, thetai, phii] -> wall grid
	"""
	import numpy as np

	r     = [ri[number[0]], ri[number[0]+1]]
	theta = [thetai[number[1]], thetai[number[1]+1]]
	phi   = [phii[number[2]], phii[number[2]+1]]
	print (1/3.)*(r[1]**3 - r[0]**3) * (phi[1]-phi[0]) * -(np.cos(theta[1])-np.cos(theta[0]))
	print rho[number[0],number[1],number[2]]
	mass = rho[number[0],number[1],number[2]] * (1/3.)*(r[1]**3 - r[0]**3) * (phi[1]-phi[0]) * -(np.cos(theta[1])-np.cos(theta[0]))

	return mass

def los_index(img_coord, incl, ri, thetai, phii):
	"""
	return the cell index along a line of sight
	img_coord: (x,y) in cm indicates the point that the line of sight intersect with the plane through the origin. 
			   Note that (x,y) are the coordinates after the inclination angle rotation.
	incl: The inclination angle in degree with respect to the polar axis (z)
	"""
	import numpy as np

	coord = np.empty((len(ri)*len(thetai)*len(phii),3),dtype='float')
	icell = 0
	for ir in range(0, len(ri)):
		for itheta in range(0, len(thetai)):
			for iphi in range(0, len(phii)):
				coord[icell] = np.array([ri[ir], thetai[itheta], phii[iphi]])
				icell += 1
	# Deal with the inclination angle
	# always rotate the model with y-axis as the rotational axis
	incl = incl * np.pi/180.0
	rot_matrix = np.matrix([[np.cos(incl), 0, np.sin(incl)],\
							[0,            1,            0],\
							[-np.sin(incl),0, np.cos(incl)]])
	coord_rot = np.empty_like(coord)
	for ic in range(0, len(coord[:,0])):
		coord_rot[ic] = (rot_matrix * coord[ic][:, np.newaxis]).T
	x = img_coord[0]
	y = img_coord[1]
	# coord_rot has the same index as coord
	theta_b = np.arctan(x/y)
	# correct for third and fourth quarants
	if y < 0:
		theta_b = theta_b + np.pi/2
	# calculate the impact parameter
	b = np.sqrt(x**2+y**2)

	# Roughly narrow down the number of cells
	r_lim = [np.sort(b-ri)[0], max(ri)]
	theta_lim = [np.sort(theta_b-thetai)[0], np.sort(thetai-theta_b)]
	if x > 0:
		phi_lim = [0, np.pi]
	else:
		phi_lim = [np.pi, 2*np.pi]
	sub_coord_rot = 
