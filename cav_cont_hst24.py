import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from hyperion.model import ModelOutput
import astropy.constants as const
from mpl_toolkits.axes_grid1 import make_axes_locatable

# constant
pc = const.pc.cgs.value

foo_reg = '/Users/yaolun/bhr71/hyperion/model11.rtout'  # it has three bands images
foo_const = '/Users/yaolun/bhr71/hyperion/model11.rtout'
foo_r2 = '/Users/yaolun/bhr71/hyperion/model11.rtout'
foo_incl = '/Users/yaolun/bhr71/hyperion/model11.rtout'

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
                 nrows_ncols=(2,4),
                 direction="row",
                 axes_pad=0.05,
                 add_all=True,
                 label_mode="1",
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="10%",
                 cbar_pad=0.05,
                 )
for i in range(4,8):
    # get the H-band simulated image
    m = ModelOutput(filename[i-4])
    group = 1
    if i == 0:
        group = 1

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

# for cavity density profiles
def scale(rc):
    r_scale = rc * (2*w)/(max(rc)-min(rc)) - w

    return r_scale

for j in range(3):
    d = density[j]
    grid[j].plot(scale(np.log10(rc)), scale(d))

fig.savefig('/Users/yaolun/Dropbox/HST_cycle24/cav_struc_H.pdf', format='pdf', dpi=300, bbox_inches='tight')
fig.clf()
