def phot_reader(datadir, objname, pacs=True, spire=True):
	import numpy as np
	from astropy.io import ascii

	filepath = datadir+str(objname).lower()+'.txt'
	data = ascii.read(filepath)
	wave = []
	flux = [] 	# in Jy

	if pacs == True:
		wave.extend((data['wavelength'][data['Name'] == 'PACS']).data.tolist())
		flux.extend((data['flux(Jy)'][data['Name'] == 'PACS']).data.tolist())
	if spire == True:
		wave.extend((data['wavelength'][data['Name'] == 'SPIRE']).data.tolist())
		flux.extend((data['flux(Jy)'][data['Name'] == 'SPIRE']).data.tolist())
	# check duplicate wavelength
	import collections
	for x, y in collections.Counter(wave).items():
		if y > 1:
			print 'Duplicate Photometry found'
			for i in range(len(wave)):
				print wave[i], flux[i]
	return {'wave': np.array(wave), 'flux': np.array(flux)}

# from pprint import pprint
# datadir = '/Users/yaolun/data/herschel_phot/'
# objname = 'TMC1'
# pprint(phot_reader(datadir, objname))