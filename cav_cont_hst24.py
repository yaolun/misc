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
size = 100.
convolve = True

# construct the cavity density profile
# set up the density profile in the outflow cavity
rc = np.sort(np.hstack((np.arange(0.14, 41253, 41252/10),60.)))
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

density = [db, dconst, d2, db]

# function for stretch
def scale(array, (data_start, data_end), (start,end)):
    # print (array.max() - array.min())
    array = (array-data_start)*(end-start)/(data_end-data_start) + start

    return array

# plot
fig = plt.figure(figsize=(12,8))
grid = ImageGrid(fig, 111,  # similar to subplot(111)
                 nrows_ncols=(2,4),
                 direction="row",
                 add_all=True,
                 label_mode="1",
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="7%",
                 cbar_pad=0.05,
                 )
w = size/2.
for i in range(4):
    d = density[i]
    r_s = scale(np.log10(rc), (np.log10(0.14), np.log10(41253)), (-w,w))
    d_s = scale(np.log10(d), (-25,-15), (-w,w))
    grid[i].plot(r_s, d_s, linewidth=1.5)

    if i == 0:
        d_ticks = scale(np.array([-24,-22,-20,-18,-16,-14]), (-25,-14), (-w,w))
        d_tick_labels = [-24,-22,-20,-18,-16,-14]
        grid[i].set_ylabel(r'$\rm{log(dust\,density)\,[g\,cm^{-3}]}$', fontsize=14)
        grid[i].set_yticks(d_ticks)
        grid[i].set_yticklabels(d_tick_labels)

        grid[i].text(-0.15, 0.1, '-24', transform=grid[i].transAxes, horizontalalignment='left', fontsize=14)
        grid[i].text(-0.15, 0.3, '-22', transform=grid[i].transAxes, horizontalalignment='left', fontsize=14)
        grid[i].text(-0.15, 0.5, '-20', transform=grid[i].transAxes, horizontalalignment='left', fontsize=14)
        grid[i].text(-0.15, 0.7, '-18', transform=grid[i].transAxes, horizontalalignment='left', fontsize=14)
        grid[i].text(-0.15, 0.9, '-16', transform=grid[i].transAxes, horizontalalignment='left', fontsize=14)
        grid[i].text(-0.28, 0.5, r'$\rm{log(dust\,density)\,[g\,cm^{-3}]}$',
                    transform=grid[i].transAxes, verticalalignment='center',
                    rotation='vertical', fontsize=16)

        # second x-label
        r_ticks = scale(np.array([0,1,2,3,4]), (np.log10(0.14), np.log10(41253)), (-w,w))
        r_tick_labels = [0,1,2,3,4]
        ax_top = grid[i].twiny()
        ax_top.set_xlabel(r'$\rm{log(radius)\,[AU]}$', fontsize=16)
        ax_top.set_xticks(r_ticks)
        ax_top.set_xticklabels(r_tick_labels)
        ax_top.tick_params('x', labelsize=14)

    else:
        r_ticks = scale(np.array([0,1,2,3,4]), (np.log10(0.14), np.log10(41253)), (-w,w))
        ax_top = grid[i].twiny()
        ax_top.set_xticks(r_ticks)
        ax_top.set_xticklabels([])
    grid[i].tick_params('both',labelsize=14)

for i in range(4,8):
    # get the H-band simulated image
    m = ModelOutput(filename[i-4])
    image = m.get_image(group=0, inclination=0, distance=178 * pc, units='MJy/sr')

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
    val = image.val[::-1, :, iwav] * factor + 1e-30

    if size != 'full':
        pix_e2c = (w-size/2.)/w * len(val[:,0])/2
        val = val[pix_e2c:-pix_e2c, pix_e2c:-pix_e2c]
        w = size/2.

    if convolve:
        from astropy.convolution import convolve, Gaussian2DKernel

        kernel = Gaussian2DKernel(0.16/(len(val[:,0])/2/w)/2.354)
        val = convolve(val, kernel)

    # This is the command to show the image. The parameters vmin and vmax are
    # the min and max levels for the colorscale (remove for default values).
    cmap = plt.cm.CMRmap
    im = grid[i].imshow(np.log10(val), vmin= -21, vmax= -12,
              cmap=cmap, origin='lower', extent=[-w, w, -w, w], aspect=1)

    grid[i].plot([0],[0], '+', color='ForestGreen', markersize=10, mew=1.5)

    cb = grid[i].cax.colorbar(im)
    cb.solids.set_edgecolor("face")
    cb.ax.minorticks_on()
    cb.ax.set_ylabel(r'$\rm{log(I_{\nu})\,[erg\,s^{-1}\,cm^{-2}\,Hz^{-1}\,sr^{-1}]}$',fontsize=16)
    cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
    plt.setp(cb_obj,fontsize=14)
    ticks_font = mpl.font_manager.FontProperties(family='STIXGeneral',size=14)
    for label in cb.ax.get_yticklabels():
        label.set_fontproperties(ticks_font)

    off_tick = [-40,-20,0,20,40]
    off_tick_label = [-40,-20,0,20,40]
    if i == 4:
        grid[i].set_xlabel(r'$\rm{RA\,Offset\,[arcsec]}$', fontsize=16)
        grid[i].set_ylabel(r'$\rm{Dec\,Offset\,[arcsec]}$', fontsize=16)
        grid[i].set_xticks(off_tick)
        grid[i].set_yticks(off_tick)
        grid[i].set_xticklabels(off_tick_label)
        grid[i].set_yticklabels(off_tick_label)

    grid[i].set_xlim([-w,w])
    grid[i].set_ylim([-w,w])
    grid[i].tick_params('both',labelsize=14)


fig.savefig('/Users/yaolun/Dropbox/HST_cycle24/cav_struc_H.pdf',
            format='pdf', dpi=300, bbox_inches='tight')
fig.clf()
