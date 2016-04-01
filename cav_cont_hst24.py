import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from hyperion.model import ModelOutput
import astropy.constants as const
from mpl_toolkits.axes_grid1 import make_axes_locatable

# constant
pc = const.pc.cgs.value

foo_reg = '/Users/yaolun/bhr71/hyperion/model12.rtout'  # it has three bands images
foo_const = '/Users/yaolun/bhr71/hyperion/model13.rtout'
foo_r2 = '/Users/yaolun/bhr71/hyperion/model14.rtout'
foo_incl = '/Users/yaolun/bhr71/hyperion/model15.rtout'

filename = [foo_reg, foo_const, foo_r2, foo_incl]

wave = 1.6

# construct the cavity density profile
# set up the density profile in the outflow cavity
rc = np.sort(np.hstack((np.arange(0.14, 41253, 41252/100),40.)))
# d2 and d15 have not yet scaled by the density offset
d2 = (rc/0.14)**-2*5e-16
dconst = np.zeros_like(rc)+1e-20
db = np.empty_like(rc)
for i in range(len(rc)):
    if rc[i] <= 60:
        db[i] = 1.0
    else:
        db[i] = (rc[i]/60.)**-2
db = db*5e-19

density = [db, dconst, d2]

# plot
# fig, axarr = plt.subplots(2, 4, sharex='col', sharey='row',figsize=(12,6))
fig = plt.figure(figsize=(12,8))
grid = ImageGrid(fig, 111,  # similar to subplot(111)
                 nrows_ncols=(1,4),
                 direction="row",
                 add_all=True,
                 label_mode="1",
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="10%",
                 cbar_pad=0.05,
                 )
for i in range(4):
    # get the H-band simulated image
    m = ModelOutput(filename[i])
    group = 0

    image = m.get_image(group=group, inclination=0, distance=178 * pc, units='MJy/sr')

    # Find the closest wavelength
    iwav = np.argmin(np.abs(wave - image.wav))

    # Calculate the image width in arcseconds given the distance used above
    # get the max radius
    rmax = max(m.get_quantities().r_wall)
    w = np.degrees(rmax / image.distance) * 3600.

    # Image in the unit of MJy/sr
    # Change it into erg/s/cm2/Hz/sr
    factor = 1e-23*1e6
    # avoid zero in log
    # flip the image, because the setup of inclination is upside down
    val = image.val[:, :, iwav] * factor + 1e-30

    # This is the command to show the image. The parameters vmin and vmax are
    # the min and max levels for the colorscale (remove for default values).
    cmap = plt.cm.CMRmap
    im = grid[i].imshow(np.log10(val), vmin= -20, vmax= -17,
              cmap=cmap, origin='lower', extent=[-w, w, -w, w], aspect=1)

    cb = grid[i].cax.colorbar(im)
    cb.solids.set_edgecolor("face")
    cb.ax.minorticks_on()
    cb.ax.set_ylabel(r'$\rm{log(I_{\nu})\,[erg\,s^{-1}\,cm^{-2}\,Hz^{-1}\,sr^{-1}]}$',fontsize=12)
    cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
    plt.setp(cb_obj,fontsize=12)
    ticks_font = mpl.font_manager.FontProperties(family='STIXGeneral',size=12)
    for label in cb.ax.get_yticklabels():
        label.set_fontproperties(ticks_font)

    grid[i].set_xlabel(r'$\rm{RA\,Offset\,[arcsec]}$', fontsize=14)
    grid[i].set_ylabel(r'$\rm{Dec\,Offset\,[arcsec]}$', fontsize=14)

fig.savefig('/Users/yaolun/Dropbox/HST_cycle24/cav_struc_H.pdf', format='pdf', dpi=300, bbox_inches='tight')
fig.clf()

def adjustFigAspect(fig,aspect=1):
    '''
    Adjust the subplot parameters so that the figure has the correct
    aspect ratio.
    '''
    xsize,ysize = fig.get_size_inches()
    minsize = min(xsize,ysize)
    xlim = .4*minsize/xsize
    ylim = .4*minsize/ysize
    if aspect < 1:
        xlim *= aspect
    else:
        ylim /= aspect
    fig.subplots_adjust(left=.5-xlim,
                        right=.5+xlim,
                        bottom=.5-ylim,
                        top=.5+ylim)

fig = plt.figure(figsize=(12,2))
grid = ImageGrid(fig, 111,  # similar to subplot(111)
                 nrows_ncols=(1,4),
                 direction="row",
                 add_all=True,
                 label_mode="1",
                 share_all=True)
for j in range(4):
    if j != 3:
        d = density[j]
        grid[j].plot(np.log10(rc), np.log10(d))

    # [grid[j].spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
    # grid[j].minorticks_on()
    grid[j].tick_params('both',which='both',bottom='off',top='off',left='off',right='off')
    # grid[j].tick_params('both',labelsize=12,width=1.5,which='minor',pad=15,length=2.5)

    grid[j].set_xlabel(r'$\rm{log(Radius)\,[AU]}$')
    grid[j].set_ylabel(r'$\rm{log(dust\,density)\,[g\,cm^{-3}]}$')


fig.savefig('/Users/yaolun/Dropbox/HST_cycle24/cav_struc_density.pdf', format='pdf', dpi=300, bbox_inches='tight')
fig.clf()
