def hyperion_sedcom(modellist, outdir, plotname, obs_data=None, labellist=None, lbol=False, legend=True, mag=1.5,\
					obs_preset='sh', dstar=1, aper=[3.6, 4.5, 5.8, 8.0, 10, 20, 24, 70, 160, 250, 350, 500, 850]):
	"""
	obs_data: dictionary which obs_data['spec'] is spectrum and obs_data['phot'] is photometry
			  obs_data['label'] = (wave, Fv, err) in um and Jy by default
	"""

	import numpy as np
	import os
	import matplotlib.pyplot as plt
	import astropy.constants as const
	from hyperion.model import ModelOutput
	from scipy.interpolate import interp1d
	from l_bol import l_bol
	import seaborn as sb
	# from seaborn import color_palette
	# from seaborn_color import seaborn_color

	# constant setup
	c = const.c.cgs.value
	pc = const.pc.cgs.value

	if labellist == None:
		if legend == True:
			print 'Model labels are not provided.  Use their filename instead.'
		labellist = []
		for i in range(0, len(modellist)):
			labellist.append(r'$\mathrm{'+os.path.splitext(os.path.basename(modellist[i]))[0]+'}$')

	# cm = seaborn_color('colorblind',len(modellist))
	sb.set(style="white")
	cm = sb.color_palette('husl', len(modellist))

	# create figure object
	fig = plt.figure(figsize=(8*mag,6*mag))
	ax = fig.add_subplot(111)
	# sb.set_style('ticks')

	print 'plotting with aperture at ', aper, 'um'

	# if the obs_data is provided than plot the observation first.  In this way, models won't be blocked by data
	if obs_data != None:
		if 'spec' in obs_data.keys():
			(wave, fv, err) = obs_data['spec']
			vfv = c/(wave*1e-4)*fv*1e-23
			l_bol_obs = l_bol(wave, fv, dstar)
			if legend == True:
				ax.text(0.75,0.9,r'$\mathrm{L_{bol}= %5.2f L_{\odot}}$' % l_bol_obs,fontsize=mag*16,transform=ax.transAxes)

			# general plotting scheme
			if obs_preset == None:
				spec, = ax.plot(np.log10(wave),np.log10(vfv),'-',color='k',linewidth=1.5*mag, label=r'$\mathrm{observations}$')
			# plot spitzer, Herschel pacs and spire in different colors
			elif obs_preset == 'sh':
				# spitzer
				spitz, = ax.plot(np.log10(wave[wave < 50]),np.log10(vfv[wave < 50]),'-',color='b',linewidth=1*mag,\
									label=r'$\mathrm{\it Spitzer}$')
				# herschel
				pacs, = ax.plot(np.log10(wave[(wave < 190.31) & (wave > 50)]),np.log10(vfv[(wave < 190.31) & (wave > 50)]),'-',\
									color='Green',linewidth=1*mag, label=r'$\mathrm{{\it Herschel}-PACS}$')
				spire, = ax.plot(np.log10(wave[wave >= 190.31]),np.log10(vfv[wave >= 190.31]),'-',color='k',linewidth=1*mag,\
									label=r'$\mathrm{{\it Herschel}-SPIRE}$')
				spec = [spitz, pacs, spire]

		if 'phot' in obs_data.keys():
			(wave_p, fv_p, err_p) = obs_data['phot']
			vfv_p = c/(wave_p*1e-4)*fv_p*1e-23
			vfv_p_err = c/(wave_p*1e-4)*err_p*1e-23
			phot, = ax.plot(np.log10(wave_p),np.log10(vfv_p),'s',mfc='DimGray',mec='k',markersize=8)
			ax.errorbar(np.log10(wave_p),np.log10(vfv_p),yerr=[np.log10(vfv_p)-np.log10(vfv_p-vfv_p_err), np.log10(vfv_p+vfv_p_err)-np.log10(vfv_p)],\
						fmt='s',mfc='DimGray',mec='k',markersize=8)

	modplot = dict()
	for imod in range(0, len(modellist)):
		m = ModelOutput(modellist[imod])
		# if not specified, distance of the star will be taken as 1 pc. 
		if aper == None:
			sed_dum = m.get_sed(group=0, inclination=0, aperture=-1, distance=dstar * pc)
			modplot['mod'+str(imod+1)], = ax_sed.plot(np.log10(sed_dum.wav), np.log10(sed_dum.val), '-', color='GoldenRod', linewidth=1.5*mag)
		else:
			vfv_aper = np.empty_like(aper)
			for i in range(0, len(aper)):
				sed_dum = m.get_sed(group=i+1, inclination=0, aperture=-1, distance=dstar * pc)
				f = interp1d(sed_dum.wav, sed_dum.val)
				vfv_aper[i] = f(aper[i])
			modplot['mod'+str(imod+1)], = ax.plot(np.log10(aper),np.log10(vfv_aper),'o',mfc='None',mec=cm[imod],markersize=12,\
													markeredgewidth=3, label=labellist[imod], linestyle='-',color=cm[imod],linewidth=1.5*mag)

	# plot fine tune
	ax.set_xlabel(r'$\mathrm{log~\lambda~({\mu}m)}$',fontsize=mag*20)
	ax.set_ylabel(r'$\mathrm{log~\nu S_{\nu}~(erg/cm^{2}/s)}$',fontsize=mag*20)
	[ax.spines[axis].set_linewidth(1.5*mag) for axis in ['top','bottom','left','right']]
	ax.minorticks_on()
	ax.tick_params('both',labelsize=mag*18,width=1.5*mag,which='major',pad=15,length=5*mag)
	ax.tick_params('both',labelsize=mag*18,width=1.5*mag,which='minor',pad=15,length=2.5*mag)

	if obs_preset == 'sh':
		ax.set_ylim([-14,-7])
		ax.set_xlim([0,3])

	if legend == True:
		lg = ax.legend(loc='best',fontsize=14*mag,numpoints=1,framealpha=0.3)

	# Write out the plot
	fig.savefig(outdir+plotname+'.pdf',format='pdf',dpi=300,bbox_inches='tight')
	fig.clf()



# import numpy as np
# from get_bhr71_obs import get_bhr71_obs

# obs_data = get_bhr71_obs('/Users/yaolun/bhr71/obs_for_radmc/')
# mod_num = [32,56]
# modellist = []
# modir = '/Users/yaolun/test/model'
# for mod in mod_num:
# 	modellist.append(modir+str(mod)+'/model'+str(mod)+'.rtout')
# outdir = '/Users/yaolun/test/'
# hyperion_sedcom(modellist, outdir, 'test', obs_data=obs_data, lbol=True, dstar=178)

