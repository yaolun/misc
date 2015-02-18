def plot_density(rho, rc, thetac, plotdir, plotname=None, zmin=1e-22, cmap='jet'):
	import numpy as np
	import matplotlib.pyplot as plt
	import astropy.constants as const
	from matplotlib.colors import LogNorm
	import matplotlib as mpl
	mpl.rcParams['text.usetex']=True

	AU = const.au.cgs.value
	PI = np.pi
	# Do the azimuthal average of the 3-D density profile
	print 'Averaging the density profile azimuthally.'
	rho2d = np.sum(rho**2,axis=2)/np.sum(rho,axis=2)

	fig = plt.figure(figsize=(8,6))
	ax_env  = fig.add_subplot(111,projection='polar')

	img_env = ax_env.pcolormesh(thetac,rc/AU,rho2d,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=np.nanmax(rho2d)))
	ax_env.pcolormesh(thetac-PI,rc/AU,rho2d,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=np.nanmax(rho2d)))

	ax_env.set_xlabel(r'$\mathrm{Polar~angle~(Degree)}$',fontsize=20)
	ax_env.set_ylabel(r'$\mathrm{Radius~(AU)}$',fontsize=20)
	ax_env.tick_params(labelsize=20)
	# ax_env.set_yticks(np.arange(0,R_env_max/AU,R_env_max/AU/5))
	# ax_env.set_ylim([0,10])
	ax_env.grid(True)
	cb = fig.colorbar(img_env)
	cb.ax.set_ylabel(r'$\mathrm{Surface~Density~(g/cm^{2})}$',fontsize=20)
	cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
	plt.setp(cb_obj,fontsize=20)

	if plotname == None:
		plotname = 'density_profile'
	print 'Saving the plot...'
	fig.savefig(plotdir+plotname+'.png', format='png', dpi=300, bbox_inches='tight')
	fig.clf()