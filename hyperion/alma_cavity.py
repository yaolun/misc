def alma_cavity(freq, outdir, vlim, units='MJy/sr'):
	import numpy as np
	import matplotlib.pyplot as plt
	import astropy.constants as const
	from hyperion.model import ModelOutput
	# constants setup
	c = const.c.cgs.value
	pc = const.pc.cgs.value
	au = const.au.cgs.value
	# Image in the unit of MJy/sr
	# Change it into erg/s/cm2/Hz/sr
	if units == 'erg/s/cm2/Hz/sr':
		factor = 1e-23*1e6
	elif units == 'MJy/sr':
		factor = 1

	# input files setup
	filename_reg = '/Users/yaolun/test/model12.rtout'
	filename_r2 = '/Users/yaolun/test/model13.rtout'
	filename_r15 = '/Users/yaolun/test/model17.rtout'
	filename_uni = '/Users/yaolun/test/model62.rtout'
	wl_aper = [3.6, 4.5, 5.8, 8.0, 8.5, 9, 9.7, 10, 10.5, 11, 16, 20, 24, 35, 70, 100, 160, 250, 350, 500, 850]
	wav = c/freq/1e9*1e4
	# wav = 40

	# read in
	# regular cavity setting
	m_reg = ModelOutput(filename_reg)
	image_reg = m_reg.get_image(group=len(wl_aper)+1, inclination=0, distance=178.0*pc, units='MJy/sr')
	# Calculate the image width in arcseconds given the distance used above
	w = np.degrees((1.5 * pc) / image_reg.distance) * 60.
	pix_num = len(image_reg.val[:,0,0])
	pix2arcsec = 2*w/pix_num
	pix2au = np.radians(2*w/pix_num/3600.)*image_reg.distance/au

	iwav = np.argmin(np.abs(wav - image_reg.wav))
	# avoid zero in log
	val_reg = image_reg.val[:, :, iwav] * factor + 1e-30

	# r^-2 cavity setting
	m_r2 = ModelOutput(filename_r2)
	image_r2 = m_r2.get_image(group=len(wl_aper)+1, inclination=0, distance=178.0*pc, units='MJy/sr')
	# Calculate the image width in arcseconds given the distance used above
	w = np.degrees((1.5 * pc) / image_r2.distance) * 60.
	pix_num = len(image_reg.val[:,0,0])
	pix2arcsec = 2*w/pix_num
	pix2au = np.radians(2*w/pix_num/3600.)*image_reg.distance/au
	iwav = np.argmin(np.abs(wav - image_r2.wav))
	# avoid zero in log
	val_r2 = image_r2.val[:, :, iwav] * factor + 1e-30

	# r^-1.5 cavity setting
	m_r15 = ModelOutput(filename_r15)
	image_r15 = m_r15.get_image(group=len(wl_aper)+1, inclination=0, distance=178.0*pc, units='MJy/sr')
	# Calculate the image width in arcseconds given the distance used above
	w = np.degrees((1.5 * pc) / image_r15.distance) * 60.
	pix_num = len(image_reg.val[:,0,0])
	pix2arcsec = 2*w/pix_num
	pix2au = np.radians(2*w/pix_num/3600.)*image_reg.distance/au
	iwav = np.argmin(np.abs(wav - image_r15.wav))
	# avoid zero in log
	val_r15 = image_r15.val[:, :, iwav] * factor + 1e-30

	# uniform cavity setting
	m_uni = ModelOutput(filename_uni)
	image_uni = m_uni.get_image(group=len(wl_aper)+1, inclination=0, distance=178.0*pc, units='MJy/sr')
	# Calculate the image width in arcseconds given the distance used above
	w = np.degrees((1.5 * pc) / image_uni.distance) * 60.
	pix_num = len(image_reg.val[:,0,0])
	pix2arcsec = 2*w/pix_num
	pix2au = np.radians(2*w/pix_num/3600.)*image_reg.distance/au
	iwav = np.argmin(np.abs(wav - image_uni.wav))
	# avoid zero in log
	val_uni = image_uni.val[:, :, iwav] * factor + 1e-30

	# 1-D radial intensity profile
	# get y=0 plane, and plot it
	fig = plt.figure(figsize=(8,6))
	ax = fig.add_subplot(111)

	reg, = ax.plot(np.linspace(-150,150,num=300)*pix2arcsec, val_reg[:,149], color='k', linewidth=1.5)
	r2,  = ax.plot(np.linspace(-150,150,num=300)*pix2arcsec, val_r2[:,149], color='b', linewidth=1.5)
	r15, = ax.plot(np.linspace(-150,150,num=300)*pix2arcsec, val_r15[:,149], color='r', linewidth=1.5)
	uni, = ax.plot(np.linspace(-150,150,num=300)*pix2arcsec, val_uni[:,149], color='g', linewidth=1.5)
	ax.legend([reg, r2, r15, uni], [r'$\mathrm{const.+r^{-2}}$', r'$\mathrm{r^{-2}}$', r'$\mathrm{r^{-1.5}}$', r'$\mathrm{uniform~density}$'],\
	          numpoints=1, loc='lower center', fontsize=16)

	ax.set_xlim([-1,1])
	ax.set_xlabel(r'$\mathrm{angular~distance~(arcsec)}$', fontsize=16)
	# ax.set_ylabel(r'$\mathrm{I_{\nu}~(erg~s^{-1}~cm^{-2}~Hz^{-1}~sr^{-1})}$', fontsize=16)
	ax.set_ylabel(r'$\mathrm{I_{\nu}~(MJy~sr^{-1})}$', fontsize=16)

	[ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
	ax.minorticks_on()
	ax.tick_params('both',labelsize=16,width=1.5,which='major',pad=15,length=5)
	ax.tick_params('both',labelsize=16,width=1.5,which='minor',pad=15,length=2.5)

	fig.savefig(outdir+'cavity_intensity_'+str(freq)+'.pdf', format='pdf', dpi=300, bbox_inches='tight')

	# 2-D intensity map
	from mpl_toolkits.axes_grid1 import AxesGrid
	image_grid = [val_reg, val_uni,val_r2, val_r15]
	fig = plt.figure(figsize=(30,30))
	grid = AxesGrid(fig, 142, # similar to subplot(142)
	                    nrows_ncols = (2, 2),
	                    axes_pad = 0.15,
	                    share_all=True,
	                    label_mode = "L",
	                    cbar_location = "right",
	                    cbar_mode="single",
	                    )
	for i in range(4):
	    im = grid[i].pcolor(np.linspace(-150,150,num=300)*pix2arcsec, np.linspace(-150,150,num=300)*pix2arcsec,\
	                       image_grid[i], cmap=plt.cm.jet, vmin=vlim[0], vmax=vlim[1])
	    grid[i].set_xlim([-2,2])
	    grid[i].set_ylim([-2,2])
	    grid[i].set_xlabel(r'$\mathrm{RA~offset~(arcsec)}$', fontsize=12)
	    grid[i].set_ylabel(r'$\mathrm{Dec~offset~(arcsec)}$', fontsize=12)
	    grid[i].tick_params('both',labelsize=10,which='major')
	    grid[i].tick_params('both',labelsize=10,which='minor')
	#     ax.set_aspect('equal')
	cb = grid.cbar_axes[0].colorbar(im)
	cb.solids.set_edgecolor("face")
	cb.ax.minorticks_on()
	cb.ax.set_ylabel(r'$\mathrm{log(I_{\nu})~[MJy/sr]}$',fontsize=12)
	cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
	plt.setp(cb_obj,fontsize=12)

	# fig.text(0.5, -0.05 , r'$\mathrm{RA~offset~(arcsec)}$', fontsize=12, ha='center')
	# fig.text(0, 0.5, r'$\mathrm{Dec~offset~(arcsec)}$', fontsize=12, va='center', rotation='vertical')

	fig.savefig(outdir+'cavity_2d_intensity_'+str(freq)+'.png',format='png',dpi=300,bbox_inches='tight')

freq = [230,345,460]
freq = [230]
vlim = [[400,2000],[600,3200],[1500,9000]]
outdir = '/Users/yaolun/test/'
for f in freq:
	alma_cavity(f,outdir, vlim[freq.index(f)])