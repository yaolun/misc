def getRadialDensity(rtout, angle, plotdir):
    """
    """
    import numpy as np
    from hyperion.model import ModelOutput


    m = ModelOutput(rtout)
    q = m.get_quantities()
    r_wall = q.r_wall; theta_wall = q.t_wall; phi_wall = q.p_wall
    # get the cell coordinates
    rc = r_wall[0:len(r_wall)-1] + 0.5*(r_wall[1:len(r_wall)]-r_wall[0:len(r_wall)-1])
    thetac = theta_wall[0:len(theta_wall)-1] + \
             0.5*(theta_wall[1:len(theta_wall)]-theta_wall[0:len(theta_wall)-1])
    phic = phi_wall[0:len(phi_wall)-1] + \
           0.5*(phi_wall[1:len(phi_wall)]-phi_wall[0:len(phi_wall)-1])
    #
    rho = q['density'].array[0]

    # find the closest angle in the thetac grid
    ind = np.argsort(abs(thetac-angle*np.pi/180.))[0]

    return rc, rho[0,ind,:]

import matplotlib.pyplot as plt
import astropy.constants as const
import numpy as np

# constants
AU = const.au.cgs.value
pc = const.pc.cgs.value
g2d = 100
mmw = 2.37
mh = const.m_p.cgs.value+const.m_e.cgs.value

rtout = '/Volumes/SD-Mac/model65.rtout'
plotdir = '/Users/yaolun/test/'
angle = 90
rc, rho = getRadialDensity(rtout, angle, plotdir)


fig  = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)

ax.plot(np.log10(rc[rc > 0.15*AU]/AU), np.log10(rho[rc > 0.15*AU]/mmw/mh))

# for comparing with K12
ax.plot(np.log10([24.8, 1000.]), np.log10([9.4e8,1.8e6]), 'o', color='r', mec='None', markersize=10)
ax.plot(np.log10([24.8, 1000.]), np.log10([9.4e6,1.8e4]), 's', color='r', mec='None', markersize=10, alpha=0.7)

ax.set_xlabel(r'$\rm{log(radius)\,[AU]}$', fontsize=18)
ax.set_ylabel(r'$\rm{log(dust\,desnity)\,[cm^{-3}]}$', fontsize=18)
[ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
ax.minorticks_on()
ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)
ax.set_ylim([0,11])
fig.gca().set_xlim(left=np.log10(0.05))

fig.savefig(plotdir+'radial_density_profile.pdf', format='pdf', dpi=300, bbox_inches='tight')
