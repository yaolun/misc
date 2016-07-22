def OutflowIntensity(model, dstar, group, wave):
    import numpy as np
    import matplotlib.pyplot as plt
    from hyperion.model import ModelOutput
    from photutils import aperture_photometry as ap
    from photutils import CircularAperture, RectangularAperture
    import astropy.constants as const

    pc = const.pc.cgs.value

    # m = ModelOutput('/Volumes/SD-Mac/model131.rtout')
    image = m.get_image(group=group, inclination=0, distance=dstar * pc, units='MJy/sr')

    # The radius of the aperture in arcsec
    # radius = 10
    x = 100
    y = 100
    # area = np.pi*radius**2 / 4.25e10 # in sr
    area = x*y / 4.25e10
    # The offset to the center
    offset = 60

    if wave != list:
        wave = [wave]

    iwav = np.argmin(np.abs(wave - image.wav))
    # Image in the unit of MJy/sr, change it into erg/s/cm2/Hz/sr
    factor = 1e-23*1e6
    val = image.val[::-1, :, iwav] * factor + 1e-30

    # Calculate the image width in arcseconds given the distance used above
    # get the max radius
    rmax = max(m.get_quantities().r_wall)
    w = np.degrees(rmax / image.distance) * 3600.

    pos_n = (len(val[0,:])/2.-1,len(val[0,:])/2.-1 + offset*len(val[0,:])/2/w)
    pos_s = (len(val[0,:])/2.-1,len(val[0,:])/2.-1 - offset*len(val[0,:])/2/w)
    # aper_n = CircularAperture(pos_n, r=radius * len(val[0,:])/2/w )
    # aper_s = CircularAperture(pos_s, r=radius * len(val[0,:])/2/w )

    # # plot to make sure the selection is correct
    # from astropy.convolution import Gaussian1DKernel, convolve
    # g = Gaussian1DKernel(stddev=20)
    #
    # fig = plt.figure(figsize=(8,6))
    # ax = fig.add_subplot(111)
    #
    # ax.plot(np.arange(-w, w, 2*w/len(val[0,:])), convolve(np.sum(val[len(val[0,:])/2.-11:len(val[0,:])/2.+9,:], axis=0), g, boundary='extend'), color='b')
    # ax.plot(np.arange(-w, w, 2*w/len(val[0,:])), convolve(np.sum(val[:,len(val[0,:])/2.-11:len(val[0,:])/2.+9], axis=1), g, boundary='extend'), color='r')
    # # ax.set_xscale('log')
    #
    # fig.savefig('/Users/yaolun/test/im_test.pdf', format='pdf', dpi=300, bbox_inches='tight')
    # fig.clf()

    aper_n = RectangularAperture(pos_n, w=x*len(val[0,:])/2/w, h=y*len(val[0,:])/2/w, theta=0)
    aper_s = RectangularAperture(pos_s, w=x*len(val[0,:])/2/w, h=y*len(val[0,:])/2/w, theta=0)

    # multiply the aperture size in sr and convert to Jy
    phot_n = ap(val, aper_n)['aperture_sum'].data * area * 1e23
    phot_s = ap(val, aper_s)['aperture_sum'].data * area * 1e23

    return phot_n, phot_s

import numpy as np
from hyperion.model import ModelOutput
# to avoid X server error
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import astropy.constants as const
pc = const.pc.cgs.value

indir = '/home/bettyjo/yaolun/hyperion/bhr71/controlled/'
# array = np.array([70,71,72,73,74])
# view_angle = np.array([90,80,70,60,50])
# indir = '/Volumes/SD-Mac/'
array = np.array([70,71,72,73,74,75,76,77])
view_angle = np.array([90,80,70,60,50,40,30,20])

wave = 3.6
# ref = 14.0
ref = 1.97

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)

for imod in range(len(array)):
    m = ModelOutput(indir+'model'+str(array[imod])+'/'+'model'+str(array[imod])+'.rtout')
    # m = ModelOutput(indir+'model'+str(array[imod])+'.rtout')

    n, s = OutflowIntensity(m, 200.0, 0, wave)
    ax.plot(view_angle[imod], s/n, 'bo', mec='None')

ax.axhline(ref, linestyle='--', color='k', linewidth=1.2)

# fix the tick label font
ticks_font = mpl.font_manager.FontProperties(family='STIXGeneral',size=18)
for label in ax.get_xticklabels():
    label.set_fontproperties(ticks_font)
for label in ax.get_yticklabels():
    label.set_fontproperties(ticks_font)

ax.set_xlim([15,95])
# ax.set_ylim([0, max([ref+0.5])])
ax.set_xlabel(r'$\rm{Inclination\,Angle\,[deg.]}$',fontsize=20)
ax.set_ylabel(r'$\rm{F_{south}/F_{north}}$',fontsize=20)
[ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
ax.minorticks_on()
ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)

fig.savefig('/home/bettyjo/yaolun/NS_comparison.pdf', format='pdf', dpi=300, bbox_inches='tight')
# fig.savefig('/Users/yaolun/test/NS_comparison.pdf', format='pdf', dpi=300, bbox_inches='tight')
