# to avoid X server error
import matplotlib as mpl
mpl.use('Agg')
#
import numpy as np
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
import astropy.constants as const

g2d = 100
mmw = 2.37
mh = const.m_p.cgs.value + const.m_e.cgs.value
AU = const.au.cgs.value

model = np.arange(99,133).astype('str')
# color map
cmap = plt.cm.viridis
color_array = [cmap(np.linspace(0.9, 0, len(model))[i]) for i in range(len(model))]

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)

for i in range(len(model)):
    m = ModelOutput('/home/bettyjo/yaolun/hyperion/bhr71/controlled/model'+model[i]+'/model'+model[i]+'.rtout')
    q = m.get_quantities()
    r = q.r_wall
    rc = 0.5*(r[0:len(r)-1]+r[1:len(r)])
    rho = q['density'][0].array
    rho2d = np.sum(rho**2,axis=0)/np.sum(rho,axis=0)
    plt.plot(np.log10(rc[rc > 0.14*AU]/AU), np.log10(rho2d[199,rc > 0.14*AU]/g2d/mmw/mh)-0.1*i, '-',
             color=color_array[i], linewidth=1)
ax.set_ylim([-2,9])
ax.set_xlabel(r'$\rm{log(Radius)\,(AU)}$',fontsize=20)
ax.set_ylabel(r'$\rm{log(Dust\,Density)\,(cm^{-3})}$',fontsize=20)
[ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
ax.minorticks_on()
ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)

# fix the tick label font
ticks_font = mpl.font_manager.FontProperties(family='STIXGeneral',size=18)
for label in ax.get_xticklabels():
    label.set_fontproperties(ticks_font)
for label in ax.get_yticklabels():
    label.set_fontproperties(ticks_font)

ax.set_ylim([0,11])
fig.gca().set_xlim(left=np.log10(0.05))
fig.savefig('/home/bettyjo/yaolun/gas_radial_age.pdf',format='pdf',dpi=300,bbox_inches='tight')
fig.clf()
