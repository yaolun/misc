def inspect_output(rtout,plotdir,quantities=None):
	import numpy as np
	import hyperion as hp
	import os
	import matplotlib.pyplot as plt
	from hyperion.model import ModelOutput
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

	# Get the dir path of rtout file
	indir = os.path.dirname(rtout)

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

		# Read in TSC-only envelope
		rho_tsc = np.genfromtxt(indir+'/rhoenv.dat').T
		# extrapolate
		def poly(x, y, x0, deg=1):
		    import numpy as np
		    p = np.polyfit(x, y, deg)
		    y0 = 0
		    for i in range(0, len(p)):
		        y0 = y0 + p[i]*x0**(len(p)-i-1)
		    return y0

		print 'Warning: hard coded infall radius (3500 AU) is used for extrapolating TSC envelope'
		r_inf = 3500 * AU
		rhoenv = rho_tsc.copy()
		for ithetac in range(0, len(thetac)):
		    rho_dum = np.log10(rhoenv[(rc > 1.1*r_inf) & (np.isnan(rhoenv[:,ithetac]) == False),ithetac])
		    rc_dum = np.log10(rc[(rc > 1.1*r_inf) & (np.isnan(rhoenv[:,ithetac]) == False)])
		#     rho_dum_nan = np.log10(rhoenv[(rc > 1.1*r_inf) & (np.isnan(rhoenv[:,ithetac]) == True),ithetac])
		    rc_dum_nan = np.log10(rc[(rc > 1.1*r_inf) & (np.isnan(rhoenv[:,ithetac]) == True)])
		    for i in range(0, len(rc_dum_nan)):
		        rho_extrapol = poly(rc_dum, rho_dum, rc_dum_nan[i])
		        rhoenv[(np.log10(rc) == rc_dum_nan[i]),ithetac] = 10**rho_extrapol
		rho_tsc = rhoenv
		rho_tsc3d = np.empty_like(rho)
		for i in range(0, len(rho[0,0,:])):
		    rho_tsc3d[:,:,i] = rho_tsc
		rho_tsc2d = np.sum(rho_tsc3d**2,axis=2)/np.sum(rho_tsc3d,axis=2)

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
		ax_env.set_ylim([0,100])
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

		plot_grid = [0,19,39,59,79,99,119,139,159,179,199]
		c_range = range(len(plot_grid))
		cNorm  = mat.colors.Normalize(vmin=0, vmax=c_range[-1])
		# color map 1
		cm1 = plt.get_cmap('Blues') 
		scalarMap1 = mat.cm.ScalarMappable(norm=cNorm, cmap=cm1)
		# color map 2
		cm2 = plt.get_cmap('Reds') 
		scalarMap2 = mat.cm.ScalarMappable(norm=cNorm, cmap=cm2)

		for i in plot_grid:
			colorVal1 = scalarMap1.to_rgba(c_range[plot_grid.index(i)])
			colorVal2 = scalarMap2.to_rgba(c_range[plot_grid.index(i)])			
			rho_plot, = ax.plot(np.log10(rc/AU), np.log10(rho2d[:,i]),'o-',color=colorVal1,linewidth=1.5, markersize=3)
			tsc_only, = ax.plot(np.log10(rc/AU), np.log10(rho_tsc2d[:,i]),'o-',color=colorVal2,linewidth=1.5, markersize=3)
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
