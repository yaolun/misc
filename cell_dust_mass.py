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
	