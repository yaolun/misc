def tsc_com(plot=True, disk=False):
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    import astropy.constants as const
    import os
    import scipy as sci
    from scipy.optimize import fsolve
    from scipy.optimize import newton
    from scipy.integrate import nquad
    import sys
    sys.path.append("/Users/yaolun/programs/misc/hyperion/")
    from input_reader import input_reader_table
    home = os.path.expanduser('~')

    # read the TSC model
    # Constant Setup
    #
    AU        = 1.49598e13     # Astronomical Unit       [cm]
    pc        = 3.08572e18     # Parsec                  [cm]
    MS        = 1.98892e33     # Solar mass              [g]
    LS        = 3.8525e33      # Solar luminosity        [erg/s]
    RS        = 6.96e10        # Solar radius            [cm]
    G         = 6.67259e-8     # Gravitational constant  [cm3/g/s^2]
    yr        = 60*60*24*365   # Years in seconds
    PI        = np.pi          # PI constant
    sigma     = const.sigma_sb.cgs.value  # Stefan-Boltzmann constant
    mh        = const.m_p.cgs.value + const.m_e.cgs.value

    # read parameter from input_table
    params_table = '/Users/yaolun/programs/misc/hyperion/test_input.txt'
    params = input_reader_table(params_table)[0]
    # force omega = 4.1e-13 to emphasize the difference
    params['Omega0'] = 4.1e-13
    #
    rstar     = params['rstar'] * RS
    tstar     = params['tstar']
    R_env_max = params['R_env_max'] * AU
    T_sub = 1600
    a     = 1   #in micron
    d_sub = (LS/16./np.pi/sigma/AU**2*(4*np.pi*rstar**2*sigma*tstar**4/LS)/T_sub**4)**0.5 *AU
    R_env_min = d_sub
    R_cen     = params['Omega0']**2 * G**3 * (0.975*(params['Cs']*1e5)**3/G*params['age']*yr)**3 /(16*(params['Cs']*1e5)**8)
    R_inf     = params['Cs']*1e5*params['age']*yr
    R_disk_min= d_sub
    R_disk_max= R_cen
    theta_cav = params['theta_cav']
    beta      = params['beta']
    h100      = params['h100'] * AU
    M_env_dot = 0.975*(params['Cs']*1e5)**3/G
    M_disk    = params['M_disk'] * MS
    mstar     = M_env_dot * params['age']*yr
    rin       = rstar
    rout      = R_env_max
    rho_cav_center = params['rho_cav_center']
    rho_cav_edge = params['rho_cav_edge'] * AU

    # Grid Parameters
    nx        = 100L
    ny        = 400L
    nz        = 50L

    # Make the Coordinates
    #
    ri           = rin * (rout/rin)**(np.arange(nx+1).astype(dtype='float')/float(nx))
    ri           = np.hstack((0.0, ri))
    thetai       = PI*np.arange(ny+1).astype(dtype='float')/float(ny)
    phii         = PI*2.0*np.arange(nz+1).astype(dtype='float')/float(nz)

    # Keep the constant cell size in r-direction
    #

    ri_cellsize = ri[1:-1]-ri[0:-2]
    ind = np.where(ri_cellsize/AU > 100.0)[0][0]       # The largest cell size is 100 AU
    ri = np.hstack((ri[0:ind],ri[ind]+np.arange(np.ceil((rout-ri[ind])/100/AU))*100*AU))
    nx = len(ri)-1

    # Assign the coordinates of the center of cell as its coordinates.
    #
    rc           = 0.5*( ri[0:nx]     + ri[1:nx+1] )
    thetac       = 0.5*( thetai[0:ny] + thetai[1:ny+1] )
    phic         = 0.5*( phii[0:nz]   + phii[1:nz+1] )

    if disk == False:
        rho_env_tsc_idl = np.genfromtxt('/Users/yaolun/test/rhoenv.dat').T
    else:
        rho_env_tsc_idl = np.genfromtxt('/Users/yaolun/bhr71/hyperion/cycle9/rhoenv_disk.dat').T

    rc_idl = rc[(rc < min([R_inf,max(ri)]))]

    # because only region within infall radius is calculated by IDL program, need to project it to the original grid
    rho_env_tsc = np.zeros([len(rc), len(thetac)])
    for irc in range(len(rc)):
        if rc[irc] in rc_idl:
            rho_env_tsc[irc,:] = rho_env_tsc_idl[np.where(rc_idl == rc[irc]),:]

    # extrapolate for the NaN values at the outer radius, usually at radius beyond the infall radius
    # using r^-2 profile at radius greater than infall radius
    # and map the 2d strcuture onto 3d grid
    def poly(x, y, x0, deg=2):
        import numpy as np
        p = np.polyfit(x, y, deg)
        y0 = 0
        for i in range(0, len(p)):
            y0 = y0 + p[i]*x0**(len(p)-i-1)
        return y0
    # map TSC solution from IDL to actual 2-D grid
    rho_env_tsc2d = np.empty((nx,ny))
    if max(ri) > R_inf:
        ind_infall = np.where(rc <= R_inf)[0][-1]
        for i in range(0, len(rc)):
            if i <= ind_infall:
                rho_env_tsc2d[i,:] = rho_env_tsc[i,:]
            else:
                rho_env_tsc2d[i,:] =  10**(np.log10(rho_env_tsc[ind_infall,:]) - 2*(np.log10(rc[i]/rc[ind_infall])))
    else:
        rho_env_tsc2d = rho_env_tsc
    # map it to 3-D grid
    rho_env_tsc = np.empty((nx,ny,nz))
    for i in range(0, nz):
        rho_env_tsc[:,:,i] = rho_env_tsc2d

    # calculate the infall-only solution

    import hyperion as hp
    from hyperion.model import Model
    from hyperion.model import AnalyticalYSOModel

    m = AnalyticalYSOModel()

    # Define the luminsoity source
    source = m.add_spherical_source()
    source.luminosity = (4*PI*rstar**2)*sigma*(tstar**4)  # [ergs/s]
    source.radius = rstar  # [cm]
    source.temperature = tstar  # [K]
    source.position = (0., 0., 0.)
    source.mass = mstar
    print 'L_center =  % 5.2f L_sun' % ((4*PI*rstar**2)*sigma*(tstar**4)/LS)

    # Envelope structure
    #
    envelope = m.add_ulrich_envelope()
    envelope.mdot = M_env_dot    # Infall rate
    envelope.rmin = rin          # Inner radius
    envelope.rc   = R_cen        # Centrifugal radius
    envelope.rmax = R_env_max    # Outer radius
    envelope.star = source

    grid = hp.grid.SphericalPolarGrid(ri, thetai, phii)

    rho_env_ulrich = envelope.density(grid).T
    rho_env_ulrich2d = np.sum(rho_env_ulrich**2,axis=2)/np.sum(rho_env_ulrich,axis=2)

    # calculate the full density field

    # Grids and Density
    # Calculation inherited from the script used for RADMC-3D

    # Make the dust density model
    # Make the density profile of the envelope
    #
    print 'Calculating the dust density profile...'
    if theta_cav != 0:
        c0 = (1e4*AU)**(-0.5)*np.sqrt(1/np.sin(np.radians(theta_cav))**3-1/np.sin(np.radians(theta_cav)))
    else:
        c0 = 0
        print 'No cavity is applied'
    rho_tsc  = np.zeros([len(rc),len(thetac),len(phic)])
    rho_ulrich = np.zeros([len(rc),len(thetac),len(phic)])
    rho_disk = np.zeros([len(rc), len(thetac), len(phic)])

    # function for normalizing the disk mass
    def f(w,z,beta,rstar,h100):
        f = 2*PI*w*(1-np.sqrt(rstar/w))*(rstar/w)**(beta+1)*np.exp(-0.5*(z/(w**beta*h100/100**beta))**2)
        return f
    rho_0 = M_disk/(nquad(f,[[R_disk_min,R_disk_max],[-R_env_max,R_env_max]], args=(beta,rstar,h100)))[0]
    #
    total_mass_tsc = 0
    total_mass_ulrich = 0
    #
    for ir in range(0,len(rc)):
        for itheta in range(0,len(thetac)):
            for iphi in range(0,len(phic)):
                if rc[ir] > R_env_min:
                    # Envelope profile
                    w = abs(rc[ir]*np.cos(np.pi/2 - thetac[itheta]))
                    z = rc[ir]*np.sin(np.pi/2 - thetac[itheta])
                    z_cav = c0*abs(w)**1.5
                    if z_cav == 0:
                        z_cav = R_env_max
                    # Cavity
                    if abs(z) > abs(z_cav):
                        # Modification for using density gradient in the cavity
                        # option for using a power law profile without constant region
                        if rho_cav_edge == 0:
                            rho_cav_edge = R_env_min
                        # the rho_cav_center is the dust density calculated from mass loss rate
                        # gas-to-dust ratio of 100 is applied after the whole calculation, therefore need to time 100 now
                        if (rc[ir] <= rho_cav_edge) & (rc[ir] >= R_env_min):
                            rho_env_tsc[ir,itheta,iphi] = 100 * rho_cav_center#*((rc[ir]/AU)**2)
                            rho_env_ulrich[ir,itheta,iphi] = 100 * rho_cav_center
                        else:
                            rho_env_tsc[ir,itheta,iphi] = 100 * rho_cav_center*(rho_cav_edge/rc[ir])**2
                            rho_env_ulrich[ir,itheta,iphi] = 100 * rho_cav_center*(rho_cav_edge/rc[ir])**2

                    # manually calculate the infall-only solution
                    else:
                        mu = abs(np.cos(thetac[itheta]))
                        # Implement new root finding algorithm
                        roots = np.roots(np.array([1.0, 0.0, rc[ir]/R_cen-1.0, -mu*rc[ir]/R_cen]))
                        if len(roots[roots.imag == 0]) == 1:
                            if (abs(roots[roots.imag == 0]) - 1.0) <= 0.0:
                                mu_o_dum = roots[roots.imag == 0]
                            else:
                                mu_o_dum = -0.5
                                print 'Problem with cubic solving, cos(theta) = ', mu_o_dum
                                print 'parameters are ', np.array([1.0, 0.0, rc[ir]/R_cen-1.0, -mu*rc[ir]/R_cen])
                        else:
                            mu_o_dum = -0.5
                            for imu in range(0, len(roots)):
                                if roots[imu]*mu >= 0.0:
                                    if (abs((abs(roots[imu]) - 1.0)) <= 1e-5):
                                        mu_o_dum = 1.0 * np.sign(mu)
                                    else:
                                        mu_o_dum = roots[imu]
                            if mu_o_dum == -0.5:
                                print 'Problem with cubic solving, roots are: ', roots
                        mu_o = mu_o_dum.real
                        rho_env_ulrich[ir,itheta,iphi] = M_env_dot/(4*PI*(G*mstar*R_cen**3)**0.5)*(rc[ir]/R_cen)**(-3./2)*(1+mu/mu_o)**(-0.5)*(mu/mu_o+2*mu_o**2*R_cen/rc[ir])**(-1)

                    # Disk profile
                    if ((w >= R_disk_min) and (w <= R_disk_max)) == True:
                        h = ((w/(100*AU))**beta)*h100
                        rho_disk[ir,itheta,iphi] = rho_0*(1-np.sqrt(rstar/w))*(rstar/w)**(beta+1)*np.exp(-0.5*(z/h)**2)

                    # Combine envelope and disk
                    rho_tsc[ir,itheta,iphi] = rho_disk[ir,itheta,iphi] + rho_env_tsc[ir,itheta,iphi]
                    rho_ulrich[ir,itheta,iphi] = rho_disk[ir,itheta,iphi] + rho_env_ulrich[ir,itheta,iphi]# rho_env_ulrich[ir,itheta,iphi]
                else:
                    rho_tsc[ir,itheta,iphi] = 1e-40
                    rho_ulrich[ir,itheta,iphi] = 1e-40
                # add the dust mass into the total count
                cell_mass_tsc = rho_tsc[ir, itheta, iphi] * (1/3.)*(ri[ir+1]**3 - ri[ir]**3) * (phii[iphi+1]-phii[iphi]) * -(np.cos(thetai[itheta+1])-np.cos(thetai[itheta]))
                total_mass_tsc = total_mass_tsc + cell_mass_tsc

                cell_mass_ulrich = rho_ulrich[ir, itheta, iphi] * (1/3.)*(ri[ir+1]**3 - ri[ir]**3) * (phii[iphi+1]-phii[iphi]) * -(np.cos(thetai[itheta+1])-np.cos(thetai[itheta]))
                total_mass_ulrich = total_mass_ulrich + cell_mass_ulrich

    print total_mass_tsc, total_mass_ulrich
    # create 2d projection
    rho_tsc2d = np.sum(rho_tsc**2,axis=2)/np.sum(rho_tsc,axis=2)
    rho_ulrich2d = np.sum(rho_ulrich**2,axis=2)/np.sum(rho_ulrich,axis=2)

    # print min(rc)/AU, max(rc)/AU

    if plot == True:

        # make plots

        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)

        plot_grid = [199]
        # alpha = np.linspace(0.3,1.0,len(plot_grid))
        alpha = [1]
        for i in plot_grid:
            tsc, = ax.plot(np.log10(rc/AU), np.log10(rho_env_tsc2d[:,i]/mh), alpha=alpha[plot_grid.index(i)], color='b', linewidth=2)
            ulrich, = ax.plot(np.log10(rc/AU), np.log10(rho_env_ulrich2d[:,i]/mh), alpha=alpha[plot_grid.index(i)], color='r', linewidth=2)

        rinf = ax.axvline(np.log10(R_inf/AU), linestyle='--', color='k', linewidth=1.5)
        cen_r = ax.axvline(np.log10(R_cen/AU), linestyle=':', color='k', linewidth=1.5)

        ax.legend([tsc, ulrich, rinf, cen_r], [r'$\rm{full\,TSC}$', r'$\rm{infall-only\,TSC}$', r'$\rm{infall\,radius}$', r'$\rm{centrifugal\,radius}$'],\
                  fontsize=16, numpoints=1, loc='lower center')

        ax.set_ylim([0, 15])
        ax.set_xlim(left=np.log10(d_sub/AU))
        ax.set_xlabel(r'$\rm{log(radius)\,[AU]}$', fontsize=18)
        ax.set_ylabel(r'$\rm{log(gas\,density)\,[g\,cm^{-3}]}$', fontsize=18)
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

        fig.savefig('/Users/yaolun/test/tsc_comparison.pdf', format='pdf', dpi=300, bbox_inches='tight')

    return rho_tsc/100, rho_ulrich/100
tsc_com(disk=False)
