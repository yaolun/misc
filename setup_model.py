def setup_model(indir,outdir,model=False,denser_wall=False,plot=False,low_res=False,flat=True,scale=1.0):
    import numpy as np
    import astropy.constants as const
    import scipy as sci
    import matplotlib.pyplot as plt
    import matplotlib as mat
    import os
    from matplotlib.colors import LogNorm
    from scipy.optimize import fsolve
    from scipy.integrate import nquad
    from envelope_func import func
    from hyperion.model import Model

    # Constants setup
    c         = const.c.cgs.value
    AU        = 1.49598e13     # Astronomical Unit       [cm]
    pc        = 3.08572e18     # Parsec                  [cm]
    MS        = 1.98892e33     # Solar mass              [g]
    LS        = 3.8525e33      # Solar luminosity        [erg/s]
    RS        = 6.96e10        # Solar radius            [cm]
    G         = 6.67259e-8     # Gravitational constant  [cm3/g/s^2]
    yr        = 60*60*24*365   # Years in seconds
    PI        = np.pi          # PI constant
    sigma     = const.sigma_sb.cgs.value  # Stefan-Boltzmann constant 


    m = Model()

    # Create dust properties

    # Hyperion needs nu, albedo, chi, g, p_lin_max
    from hyperion.dust import HenyeyGreensteinDust
    # Read in the dust opacity table used by RADMC-3D
    dust_radmc = dict()
    [dust_radmc['wl'], dust_radmc['abs'], dust_radmc['scat'], dust_radmc['g']] = np.genfromtxt('dustkappa_oh5_extended.inp',skip_header=2).T
    # opacity per mass of dust?
    dust_hy = dict()
    dust_hy['nu'] = c/dust_radmc['wl']*1e4
    ind = np.argsort(dust_hy['nu'])
    dust_hy['nu'] = dust_hy['nu'][ind]
    dust_hy['albedo'] = (dust_radmc['scat']/(dust_radmc['abs']+dust_radmc['scat']))[ind]
    dust_hy['chi'] = (dust_radmc['abs']+dust_radmc['scat'])[ind]
    dust_hy['g'] = dust_radmc['g'][ind]
    dust_hy['p_lin_max'] = 0*dust_radmc['wl'][ind]     # assume no polarization

    d = HenyeyGreensteinDust(dust_hy['nu'], dust_hy['albedo'], dust_hy['chi'], dust_hy['g'], dust_hy['p_lin_max'])
    # dust sublimation does not occur
    # d.set_sublimation_temperature(None)
    d.write(outdir+'oh5.hdf5')
    d.plot(outdir+'oh5.png')

    # Grids and Density
    # Calculation inherited from the script used for RADMC-3D

    # Parameters setup
    # Import the model parameters from another file 
    #
    params     = np.genfromtxt(indir+'/params.dat',dtype=None)
    tstar      = params[0][1]
    mstar      = params[1][1]*MS
    rstar      = params[2][1]*RS
    M_env_dot  = params[3][1]*MS/yr
    M_disk_dot = params[4][1]*MS/yr
    R_env_max  = params[5][1]*AU
    R_env_min  = params[6][1]*AU
    theta_cav  = params[7][1]
    R_disk_max = params[8][1]*AU
    R_disk_min = params[9][1]*AU
    R_cen      = R_disk_max
    M_disk     = params[10][1]*MS
    beta       = params[11][1]
    h100       = params[12][1]*AU
    rho_cav    = params[13][1]
    if denser_wall == True:
        wall       = params[14][1]*AU
        rho_wall   = params[15][1]
    rho_cav_center = params[16][1]
    rho_cav_edge   = params[17][1]*AU
    
    # Grid Parameters
    nx        = 300L
    if low_res == True:
        nx    = 100L
    ny        = 400L
    nz        = 50L
    [nx, ny, nz] = [scale*nx, scale*ny, scale*nz]
    # nx        = 20
    # ny        = 40
    # nz        = 5

    
    # Model Parameters
    #
    rin       = rstar
    rout      = R_env_max
    rcen      = R_cen

    # Star Parameters
    #
    mstar    = mstar
    rstar    = rstar*0.9999
    tstar    = tstar
    pstar    = [0.,0.,0.]

    # Make the Coordinates
    #
    ri           = rin * (rout/rin)**(np.arange(nx+1).astype(dtype='float')/float(nx))
    ri           = np.hstack((0.0, ri))
    thetai       = PI*np.arange(ny+1).astype(dtype='float')/float(ny)
    phii         = PI*2.0*np.arange(nz+1).astype(dtype='float')/float(nz)
    
    # Keep the constant cell size in r-direction
    #
    if flat == True:
        ri_cellsize = ri[1:-1]-ri[0:-2]
        ind = np.where(ri_cellsize/AU > 100.0)[0][0]       # The largest cell size is 100 AU
        ri = np.hstack((ri[0:ind],ri[ind]+np.arange(np.ceil((rout-ri[ind])/100/AU))*100*AU))
        nx = len(ri)-1    

    # Assign the coordinates of the center of cell as its coordinates.
    #
    rc           = 0.5*( ri[0:nx]     + ri[1:nx+1] )
    thetac       = 0.5*( thetai[0:ny] + thetai[1:ny+1] )
    phic         = 0.5*( phii[0:nz]   + phii[1:nz+1] )
    # phic         = 0.5*( phii[0:nz-1]   + phii[1:nz] )

    # Make the dust density model
    # Make the density profile of the envelope
    #
    print 'Calculating the dust density profile...'
    if theta_cav != 0:
        c = R_env_max**(-0.5)*np.sqrt(1/np.sin(np.radians(theta_cav))**3-1/np.sin(np.radians(theta_cav)))
    else:
        c = 0
    rho_env  = np.zeros([len(rc),len(thetac),len(phic)])
    rho_disk = np.zeros([len(rc),len(thetac),len(phic)])
    rho      = np.zeros([len(rc),len(thetac),len(phic)])
    def f(w,z,beta,rstar,h100):
        f = 2*PI*w*(1-np.sqrt(rstar/w))*(rstar/w)**(beta+1)*np.exp(-0.5*(z/(w**beta*h100/100**beta))**2)
        return f
    rho_0 = M_disk/(nquad(f,[[R_disk_min,R_disk_max],[-R_env_max,R_env_max]], args=(beta,rstar,h100)))[0]
    i = 0
    j = 0
    if 'rho_cav_center' in locals() == False:
        rho_cav_center = 5.27e-18 # 1.6e-17  # 5.27e-18
        print 'Use 5.27e-18 as the default value for cavity center'
    if 'rho_cav_edge' in locals() == False:
        rho_cav_edge = 40*AU
        print 'Use 40 AU as the default value for size of the inner region'
    discont = 1
    if denser_wall == False:
        for ir in range(0,len(rc)):
            for itheta in range(0,len(thetac)):
                for iphi in range(0,len(phic)):
                    if rc[ir] > R_env_min:
                        # Envelope profile
                        w = abs(rc[ir]*np.cos(thetac[itheta]))
                        z = rc[ir]*np.sin(thetac[itheta])
                        z_cav = c*abs(w)**1.5
                        if z_cav == 0:
                            z_cav = R_env_max
                        if abs(z) > abs(z_cav):
                            # rho_env[ir,itheta,iphi] = rho_cav
                            # Modification for using density gradient in the cavity
                            if rc[ir] <= rho_cav_edge:
                                rho_env[ir,itheta,iphi] = rho_cav_center#*((rc[ir]/AU)**2)
                            else:
                                rho_env[ir,itheta,iphi] = rho_cav_center*discont*(rho_cav_edge/rc[ir])**2
                            i += 1
                        else:
                            j += 1
                            mu = abs(np.cos(thetac[itheta]))
                            mu_o = np.abs(fsolve(func,[0.5,0.5,0.5],args=(rc[ir],rcen,mu))[0])
                            rho_env[ir,itheta,iphi] = M_env_dot/(4*PI*(G*mstar*rcen**3)**0.5)*(rc[ir]/rcen)**(-3./2)*(1+mu/mu_o)**(-0.5)*(mu/mu_o+2*mu_o**2*rcen/rc[ir])**(-1)
                        # Disk profile
                        if ((w >= R_disk_min) and (w <= R_disk_max)) == True:
                            h = ((w/(100*AU))**beta)*h100
                            rho_disk[ir,itheta,iphi] = rho_0*(1-np.sqrt(rstar/w))*(rstar/w)**(beta+1)*np.exp(-0.5*(z/h)**2)
                        # Combine envelope and disk
                        rho[ir,itheta,iphi] = rho_disk[ir,itheta,iphi] + rho_env[ir,itheta,iphi]
                    else:
                        rho[ir,itheta,iphi] = 1e-30
        rho_env  = rho_env  + 1e-40
        rho_disk = rho_disk + 1e-40
        rho      = rho      + 1e-40
    else:
        for ir in range(0,len(rc)):
            for itheta in range(0,len(thetac)):
                for iphi in range(0,len(phic)):
                    # Envelope profile
                    w = abs(rc[ir]*np.cos(thetac[itheta]))
                    z = rc[ir]*np.sin(thetac[itheta])
                    z_cav = c*abs(w)**1.5
                    z_cav_wall = c*abs(w-wall)**1.5
                    if z_cav == 0:
                        z_cav = R_env_max
                    if abs(z) > abs(z_cav):
                        # rho_env[ir,itheta,iphi] = rho_cav
                        # Modification for using density gradient in the cavity
                        if rc[ir] <= 20*AU:
                            rho_env[ir,itheta,iphi] = rho_cav_center*((rc[ir]/AU)**2)
                        else:
                            rho_env[ir,itheta,iphi] = rho_cav_center*discont*(20*AU/rc[ir])**2
                        i += 1
                    elif (abs(z) > abs(z_cav_wall)) and (abs(z) < abs(z_cav)):
                        rho_env[ir,itheta,iphi] = rho_wall
                    else:
                        j += 1
                        mu = abs(np.cos(thetac[itheta]))
                        mu_o = np.abs(fsolve(func,[0.5,0.5,0.5],args=(rc[ir],rcen,mu))[0])
                        rho_env[ir,itheta,iphi] = M_env_dot/(4*PI*(G*mstar*rcen**3)**0.5)*(rc[ir]/rcen)**(-3./2)*(1+mu/mu_o)**(-0.5)*(mu/mu_o+2*mu_o**2*rcen/rc[ir])**(-1)
                    # Disk profile
                    if ((w >= R_disk_min) and (w <= R_disk_max)) == True:
                        h = ((w/(100*AU))**beta)*h100
                        rho_disk[ir,itheta,iphi] = rho_0*(1-np.sqrt(rstar/w))*(rstar/w)**(beta+1)*np.exp(-0.5*(z/h)**2)
                    # Combine envelope and disk
                    rho[ir,itheta,iphi] = rho_disk[ir,itheta,iphi] + rho_env[ir,itheta,iphi]
        rho_env  = rho_env  + 1e-40
        rho_disk = rho_disk + 1e-40
        rho      = rho      + 1e-40 

    # Insert the calculated grid and dust density profile into hyperion
    m.set_spherical_polar_grid(ri, thetai, phii)
    m.add_density_grid(rho.T, outdir+'oh5.hdf5')    # numpy read the array in reverse order

    # Define the luminsoity source
    source = m.add_spherical_source()
    source.luminosity = (4*PI*rstar**2)*sigma*(tstar**4)  # [ergs/s]
    source.radius = rstar  # [cm]
    source.temperature = tstar  # [K]
    source.position = (0., 0., 0.)
    print 'L_center =  % 5.2f L_sun' % ((4*PI*rstar**2)*sigma*(tstar**4)/LS)

    # Setting up images and SEDs
    image = m.add_peeled_images()
    image.set_wavelength_range(300, 2.0, 670.0)
    # pixel number
    image.set_image_size(300, 300)
    image.set_image_limits(-R_env_max, R_env_max, -R_env_max, R_env_max)
    image.set_viewing_angles([82.0], [0.0])
    image.set_uncertainties(True)
    # output as 64-bit
    image.set_output_bytes(8)

    # Radiative transfer setting

    # number of photons for temp and image
    m.set_raytracing(True)
    m.set_n_photons(initial=1000000, imaging=1000000, raytracing_sources=1000000, raytracing_dust=1000000)
    # number of iteration to compute dust specific energy (temperature)
    m.set_n_initial_iterations(20)
    m.set_convergence(True, percentile=99., absolute=1.5, relative=1.02)
    m.set_mrw(True)   # Gamma = 1 by default

    # Output setting
    # Density
    m.conf.output.output_density = 'last'

    # Density difference (shows where dust was destroyed)
    m.conf.output.output_density_diff = 'none'

    # Energy absorbed (using pathlengths)
    m.conf.output.output_specific_energy = 'last'

    # Number of unique photons that passed through the cell
    m.conf.output.output_n_photons = 'last'

    m.write(outdir+'best_model.rtin')






indir = '/Users/yaolun/bhr71/radmc3d_params'
outdir = '/Users/yaolun/bhr71/hyperion/'
setup_model(indir,outdir,low_res=True,scale=1.0)