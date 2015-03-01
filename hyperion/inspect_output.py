def inspect_output(rtout,plotdir,quantities=None):
	import numpy as np
	import hyperion as hp
	import matplotlib.pyplot as plt
	from hyperion.model import *
	import astropy.constants as const
	import matplotlib as mat
	from matplotlib.colors import LogNorm


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

	m = ModelOutput(rtout)
	grid = m.get_quantities()
	rc     = 0.5*(grid.r_wall[0:-1]+grid.r_wall[1:])
	thetac = 0.5*(grid.t_wall[0:-1]+grid.t_wall[1:])
	phic   = 0.5*(grid.p_wall[0:-1]+grid.p_wall[1:])

	print 'Only works for density now'
	if quantities == None:
		quantities = input('What quantity you want to take a look at? ')
	elif quantities == 'density':
		rho = grid[quantities][0].array.T
		rho2d = np.sum(rho**2,axis=2)/np.sum(rho,axis=2)

        rho2d_exp = np.hstack((rho2d,rho2d,rho2d[:,0:1]))
        thetac_exp = np.hstack((thetac-PI/2, thetac+PI/2, thetac[0]-PI/2))

		# Make the plot
		fig = plt.figure(figsize=(8,6))
		ax_env  = fig.add_subplot(111,projection='polar')

		zmin = 1e-22/mh
		cmap = 'jet'
		img_env = ax_env.pcolormesh(thetac_exp,rc/AU,rho2d_exp/mh,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=np.nanmax(rho2d_exp/mh)))
		ax_env.pcolormesh(thetac_exp-PI,rc/AU,rho2d_exp/mh,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=np.nanmax(rho2d_exp/mh)))

		ax_env.set_xlabel(r'$\mathrm{Polar~angle~(Degree)}$',fontsize=20)
		ax_env.set_ylabel(r'$\mathrm{Radius~(AU)}$',fontsize=20)
		ax_env.tick_params(labelsize=20)
		# ax_env.set_yticks(np.arange(0,R_env_max/AU,R_env_max/AU/5))
		ax_env.set_xticklabels([r'$\mathrm{90^{\circ}}$',r'$\mathrm{45^{\circ}}$',r'$\mathrm{0^{\circ}}$',r'$\mathrm{-45^{\circ}}$',\
								r'$\mathrm{-90^{\circ}}$',r'$\mathrm{-135^{\circ}}$',r'$\mathrm{180^{\circ}}$',r'$\mathrm{135^{\circ}}$'])
		ax_env.set_ylim([0,10000])
		ax_env.grid(True)
		cb = fig.colorbar(img_env, pad=0.1)
		cb.ax.set_ylabel(r'$\mathrm{Averaged~Density~(cm^{-3})}$',fontsize=20)
		cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
		plt.setp(cb_obj,fontsize=20)
		fig.savefig(plotdir+'dust_density.png',format='png',dpi=300,bbox_inches='tight')
		fig.clf()

		# Radial density plot
		fig = plt.figure(figsize=(12,9))
		ax = fig.add_subplot(111)

		plot_grid = [0,39,79,119,159,199]

		for i in plot_grid:
			rho_plot,  = ax.plot(np.log10(rc/AU), np.log10(rho2d[:,i]),'-',color='b',linewidth=1.5, markersize=3)

		# lg = plt.legend([wrong, wrong2, wrong_mid, wrong2_mid],\
		#                 [r'$\mathrm{Before~fixing~\theta~(pole)}$',r'$\mathrm{After~fixing~\theta~(pole)}$',r'$\mathrm{Before~fixing~\theta~(midplane)}$',r'$\mathrm{After~fixing~\theta~(midplane)}$'],\
		#                 fontsize=20, numpoints=1)
		ax.set_xlabel(r'$\mathrm{log(Radius)~(AU)}$',fontsize=20)
		ax.set_ylabel(r'$\mathrm{log(Density)~(g~cm^{-3})}$',fontsize=20)
		[ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
		ax.minorticks_on()
		ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
		ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)
		ax.set_ylim([-23,-11])
		ax.set_xlim([np.log10(0.8),np.log10(10000)])

		fig.savefig(plotdir+'radial_density.pdf',format='pdf',dpi=300,bbox_inches='tight')
		fig.clf()