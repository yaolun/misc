def ch_plus(J_low, J_hi):
	# 	 FREQ	  ERR   LGINT DR	  ELO GUP	    TAG QNFMT, QN0, QN0
 #  835078.9500  0.0750 -0.0958 2    0.0000  3 -130031303 1 0 0       0 0 0       
 # 1669173.3384  0.5485  0.7208 2   27.8552  5  130031303 2 0 0       1 0 0       
 # 2501299.9662  2.0000  1.1049 2   83.5329  7  130031303 3 0 0       2 0 0       
 # 3330478.3679  4.5342  1.2782 2  166.9673  9  130031303 4 0 0       3 0 0       
 # 4155732.1742  7.9429  1.3102 2  278.0601 11  130031303 5 0 0       4 0 0       
 # 4976090.4780 11.6492  1.2320 2  416.6804 13  130031303 6 0 0       5 0 0       
 # 5790589.2001 14.9534  1.0607 2  582.6649 15  130031303 7 0 0       6 0 0       
 # 6598272.4546 18.2536  0.8063 2  775.8181 17  130031303 8 0 0       7 0 0       
 # 7398193.9149 26.2137  0.4757 2  995.9128 19  130031303 9 0 0       8 0 0  
 	import numpy as np
 	import astropy.constants as const

 	c = const.c.cgs.value
 	h = const.h.cgs.value
 	k = const.k_B.cgs.value
 	h_bar = h/2/np.pi
	pi = np.pi

	J = np.arange(J_low, J_hi+1).astype('float')
 	v = 1e6*np.array([ 835078.9500,1669173.3384,2501299.9662,3330478.3679,4155732.1742,4976090.4780,5790589.2001,6598272.4546,7398193.9149])
 	wave = c/v*1e4
 	E_u = np.array([ 27.8552,83.5329,166.9673,278.0601,416.6804,582.6649,775.8181,995.9128,995.9128+h*7398193.9149*1e6/k])

 	mu = 1.683 * 1e-18 # reference CDMS catalog
 	A = (8*pi*h*v**3/c**3)*((2*J-1)/(2*J+1))*(8*pi**3/3/h**2)*mu**2*(J)/(2*J-1)

 	for i in range(len(J)):
 		print 'CH+ J=%2d-%2d: Wave= %.8f um, Eu = %.8f K, A = %.8e s-1, g = %2d' % (J[i], J[i]-1, c/v[i]*1e4, E_u[i], A[i], 2*J[i]+1)

# ch_plus(1,9)

def count_ch_plus(filepath):
	import numpy as np
	import astropy.io.ascii as ascii
	chline = ['CH+1-0','CH+2-1','CH+3-2','CH+4-3','CH+5-4','CH+6-5']
	line = 0
	for path in filepath:
		data = ascii.read(path)
		for ch in chline:
			print ch, len(data[(data['Validity'] == 1) & (data['SNR'] >= 5) & (data['Line'] == ch) & (data['Str(W/cm2)'] > 0)])
			print data[(data['Validity'] == 1) & (data['SNR'] >= 5) & (data['Line'] == ch) & (data['Str(W/cm2)'] > 0)]
			line += len(data[(data['Validity'] == 1) & (data['SNR'] >= 5) & (data['Line'] == ch) & (data['Str(W/cm2)'] > 0)])
	print line
filepath = ['/Users/yaolun/test/CDF_archive_pacs_cube_lines.txt','/Users/yaolun/test/CDF_archive_spire_cube_lines.txt']
count_ch_plus(filepath)
