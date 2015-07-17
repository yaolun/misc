def spitzer_unc(filename, R=60., width=2.5):
	"""
	R is the resolving power (lambda/delta_lambda)
	width = number of resolution elements
	"""
	import numpy as np
	import matplotlib.pyplot as plt

	[wl_irs, flux_irs]= (np.genfromtxt(filename,skip_header=2,dtype='float').T)[0:2]
	# Remove points with zero or negative flux 
	ind = (flux_irs > 0) & (np.isnan(flux_irs) == False)
	wl_irs = wl_irs[ind]
	flux_irs = flux_irs[ind]
	unc_irs = np.empty_like(flux_irs)

	oversample = (wl_irs[1]-wl_irs[0] + wl_irs[2]-wl_irs[1])/2 / (wl_irs[1]/R)
	# print oversample

	# fig = plt.figure(figsize=(8,6))
	# ax = fig.add_subplot(111)
	# ax.plot(wl_irs, flux_irs, color='k')
	j = 0
	edge = []
	for i in range(len(wl_irs)):
		if (wl_irs[i]-width/2 * wl_irs[i]/R >= min(wl_irs)) and (wl_irs[i]+width/2 * wl_irs[i]/R <= max(wl_irs)):
			wl_dum = wl_irs[(wl_irs >= wl_irs[i]-width/2*wl_irs[i]/R) & (wl_irs <= wl_irs[i]+width/2*wl_irs[i]/R)]
			flux_dum = flux_irs[(wl_irs >= wl_irs[i]-width/2*wl_irs[i]/R) & (wl_irs <= wl_irs[i]+width/2*wl_irs[i]/R)]
			# return the coefficient, highest power first.
			fit_dum = np.polyfit(wl_dum, flux_dum, 3)
			# base_dum = fit_dum[0]*wl_dum**2 + fit_dum[1]*wl_dum + fit_dum[2]
			base_dum = fit_dum[0]*wl_dum**3 + fit_dum[1]*wl_dum**2 + fit_dum[2]*wl_dum + fit_dum[3]
			# base_dum[base_dum <= 0] = base_dum[base_dum <= 0]
			# try fit a straight line to prevent a negative baseline value
			# fit_dum = np.polyfit(wl_dum, flux_dum, 1)
			# base_dum = fit_dum[0]*wl_dum + fit_dum[1]
			# base_dum = fit_dum[0]*wl_dum + fit_dum[1]
			# ax.plot(wl_dum, base_dum, color='b')

			unc_irs[i] = np.std(flux_dum-base_dum) / np.sqrt(oversample)
			if j == 0:
				edge.append(unc_irs[i])
			j += 1
			edge_dum = unc_irs[i]
	edge.append(edge_dum)
	# print edge
	for i in range(len(wl_irs)):
		if wl_irs[i]-width/2 * wl_irs[i]/R < min(wl_irs):
			unc_irs[i] = edge[0]
		if wl_irs[i]+width/2 * wl_irs[i]/R > max(wl_irs):
			unc_irs[i] = edge[1]
		if flux_irs[i] - unc_irs[i] < 0:
			unc_irs[i] = 1/3. * flux_irs[i]

	# fig.savefig('/Users/yaolun/test/spitzer_unc_plot.pdf', format='pdf', dpi=300, bbox_inches='tight')
	# fig.clf()

	# fig = plt.figure(figsize=(8,6))
	# ax = fig.add_subplot(111)

	# ax.plot(wl_irs, flux_irs, color='k')
	# ax.errorbar(wl_irs, flux_irs, unc_irs, color='k')

	# fig.savefig('/Users/yaolun/test/spitzer_unc_plot_spec.pdf', format='pdf', dpi=300, bbox_inches='tight')

	return wl_irs, flux_irs, unc_irs

# import numpy as np
# filename = '/Users/yaolun/bhr71/obs_for_radmc/bhr71_spitzer_irs.txt'
# wl_irs, flux_irs, unc_irs = spitzer_unc(filename)
# print np.log10(flux_irs-unc_irs)