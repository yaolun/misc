def co_A(J_low, J_hi):
	# Calculate the Einstein-A, upper energy and g-value for 12C16O
	import numpy as np
	import astropy.constants as const

	# Constants setup
	mu = 0.10984*1e-18
	# mu = 0.110*1e-18
	c = const.c.cgs.value
	h = const.h.cgs.value
	k = const.k_B.cgs.value
	h_bar = h/2/np.pi
	pi = np.pi

	J = np.arange(J_low, J_hi+1).astype('float')
	# Calculate the frequency using co_energy.so
	def E(Ji, vi = 0, delta_J = 1, delta_v = 0, massC = 12.000000, massO = 15.994915, u_only=False):
		import co_energy as co
		# Energy in unit of MHz
		Ei = co.co_energy(vi, Ji, massC, massO)
		Ef = co.co_energy(vi-delta_v, Ji-delta_J, massC, massO)
		if u_only == False:
			return (Ei-Ef)*1e6
		else:
			return (Ei-co.co_energy(vi, 0, massC, massO))*1e6

	v = J*0
	E_u = J*0
	for i in range(0,len(J)):
		v[i] = E(J[i])
		E_u[i] = E(J[i], u_only=True)*h/k
	# frequency from J=40-39 to J=48-47
	# v = np.array([4.5640055992E+12,4.6756792650E+12,4.7871737770E+12,4.8984848364E+12,5.0096081499E+12,5.1205394291E+12,5.2312743909E+12,\
	# 	 5.3418087577E+12,5.4521382575E+12])
	# J = np.array([1.0])
	# v = np.array([1.1527120118E+11])

	# Calculate the Einstein-A
	# Consider the distortion of the rotation stretching the bond resulting in the change of the dipole moment.
	# Goorvitch 1994 ApJS 95 535G
	M0 = -1.1013e-1
	def F(Ji,delta_J):
		Jf = Ji - delta_J
		m = (Ji*(Ji+1)-Jf*(Jf+1))/2.
		# The coefficients of the polynomial fit for v = 0
		b0 = 9.9985e-1
		b1 = 3.0160e-7
		b2 = -1.0290e-4
		b3 = -2.8660e-11
		b4 = -3.4601e-10
		return b0*m**0 + b1*m**1 + b2*m**2 + b3*m**3 + b4*m**4
	mu = (M0 * F(J,1.0)) * 1e-18
	A = (8*pi*h*v**3/c**3)*((2*J-1)/(2*J+1))*(8*pi**3/3/h**2)*mu**2*(J)/(2*J-1)
	# upper energy
	# CO39-38: 4293.64 K

	# CO J=1-0
	# A = 7.203e-8
	# E_tran = h*v/k
	# E_u = E_tran*0
	# for i in range(0,len(E_tran)):
	# 	E_u[i] = 4293.64 + np.sum(E_tran[0:i+1])

	for i in range(0,len(J)):
		print 'CO J=%2d-%2d: Wave= %.8f um, Eu = %.8f K, A = %.8e s-1, g = %2d' % (J[i], J[i]-1, c/v[i]*1e4, E_u[i], A[i], 2*J[i]+1)

co_A(1,49)