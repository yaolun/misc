def setup_model_shell(indir,outdir,outname,rin_shell=None,denser_wall=False,tsc=True,idl=False,plot=False,low_res=False,flat=True,scale=1.0):
    import numpy as np
    import astropy.constants as const
    import scipy as sci
    import matplotlib.pyplot as plt
    import matplotlib as mat
    import os
    from matplotlib.colors import LogNorm
    from scipy.optimize import fsolve
    from scipy.optimize import newton
    from scipy.integrate import nquad
    from envelope_func import func
    import hyperion as hp
    from hyperion.model import Model
    from plot_density import plot_density

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

    # Grid Parameters
    nx        = 300L
    if low_res == True:
        nx    = 100L
    ny        = 400L
    nz        = 50L
    [nx, ny, nz] = [scale*nx, scale*ny, scale*nz]

    if tsc == False:
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
    else:
        # TSC model input setting
        params    = np.genfromtxt(indir+'/tsc_params.dat', dtype=None)
        # TSC model parameter
        M_env_dot = params[0][1]*MS/yr
        R_cen     = params[1][1]*AU
        R_inf     = params[2][1]*AU
        # protostar parameter
        tstar     = params[3][1]
        R_env_max = params[4][1]*AU
        theta_cav = params[5][1]
        rho_cav_center = params[6][1]
        rho_cav_edge   = params[7][1]*AU
        rstar     = params[8][1]*RS
        # Calculate the dust sublimation radius
        T_sub = 2000
        a     = 1   #in micron
        d_sub = (306.86*(a/0.1)**-0.4 * (4*np.pi*rstar**2*sigma*tstar**4/LS) / T_sub)**0.5 *AU
        # use the dust sublimation radius as the inner radius of disk and envelope
        R_disk_min = d_sub
        R_env_min  = d_sub
        rin        = rstar
        rout       = R_env_max
        R_disk_max = R_cen
        # mostly fixed parameter
        M_disk    = 0.5*MS
        beta      = 1.093
        h100      = 8.123*AU
        rho_cav   = 1e-21

        # Do the variable conversion
        cs = (G * M_env_dot / 0.975)**(1/3.)  # cm/s
        t = R_inf / cs / yr   # in year
        mstar = M_env_dot * t * yr
        omega = (R_cen * 16*cs**8 / (G**3 * mstar**3))**0.5

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
    if rin_shell == None:
        rin_shell = 0.3*R_env_max
    if tsc == False:
        print 'Calculating the dust density profile with infall solution...'

        rho_env  = np.zeros([len(rc),len(thetac),len(phic)])
        rho      = np.zeros([len(rc),len(thetac),len(phic)])

        for ir in range(0,len(rc)):
            for itheta in range(0,len(thetac)):
                for iphi in range(0,len(phic)):
                    if rc[ir] > rin_shell:
                        # Envelope profile
                        mu = abs(np.cos(PI/2 - thetac[itheta]))
                        # Implement new root finding algorithm
                        roots = np.roots(np.array([1.0, 0.0, rc[ir]/rcen-1.0, -mu*rc[ir]/rcen]))
                        if len(roots[roots.imag == 0]) == 1:
                            if (abs(roots[roots.imag == 0]) - 1.0) <= 0.0:
                                mu_o_dum = roots[roots.imag == 0]
                            else:
                                mu_o_dum = -0.5
                                print 'Problem with cubic solving, cos(theta) = ', mu_o_dum
                                print 'parameters are ', np.array([1.0, 0.0, rc[ir]/rcen-1.0, -mu*rc[ir]/rcen])
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
                        rho_env[ir,itheta,iphi] = M_env_dot/(4*PI*(G*mstar*rcen**3)**0.5)*(rc[ir]/rcen)**(-3./2)*(1+mu/mu_o)**(-0.5)*(mu/mu_o+2*mu_o**2*rcen/rc[ir])**(-1)
                        rho[ir,itheta,iphi] = rho_env[ir,itheta,iphi]
                    else:
                        rho[ir,itheta,iphi] = 1e-25
        rho_env  = rho_env  + 1e-40
        rho      = rho      + 1e-40
    # TSC model
    else:
        print 'Calculating the dust density profile with TSC solution...'
        # If needed, calculate the TSC model via IDL
        #
        if idl == True:
            print 'Using IDL to calculate the TSC model.  Make sure you are running this on mechine with IDL.'
            import pidly
            idl = pidly.IDL('/Applications/exelis/idl82/bin/idl')
            idl('.r ~/programs/misc/TSC/tsc.pro')
            idl.pro('tsc_run', outdir=outdir, grid=[nx,ny,nz], time=t, c_s=cs, omega=omega, rstar=rstar, renv_min=R_env_min, renv_max=R_env_max)
        else:
            print 'Read the pre-computed TSC model.'
        # read in the exist file
        rho_env_tsc = np.genfromtxt(outdir+'rhoenv.dat').T
        print np.shape(rho_env_tsc)
        # extrapolate for the NaN values at the outer radius, usually at radius beyond the infall radius
        # map the 2d strcuture onto 3d grid
        def poly(x, y, x0, deg=1):
            import numpy as np
            p = np.polyfit(x, y, deg)
            y0 = 0
            for i in range(0, len(p)):
                y0 = y0 + p[i]*x0**(len(p)-i-1)
            return y0
        rho_env_copy = np.array(rho_env_tsc)
        for ithetac in range(0, len(thetac)):
            rho_dum = np.log10(rho_env_copy[(rc > 1.1*R_inf) & (np.isnan(rho_env_copy[:,ithetac]) == False),ithetac])
            rc_dum = np.log10(rc[(rc > 1.1*R_inf) & (np.isnan(rho_env_copy[:,ithetac]) == False)])
            rc_dum_nan = np.log10(rc[(rc > 1.1*R_inf) & (np.isnan(rho_env_copy[:,ithetac]) == True)])
            for i in range(0, len(rc_dum_nan)):
                rho_extrapol = poly(rc_dum, rho_dum, rc_dum_nan[i])
                rho_env_copy[(np.log10(rc) == rc_dum_nan[i]),ithetac] = 10**rho_extrapol
        rho_env2d = rho_env_copy
        rho_env = np.empty((nx,ny,nz))
        for i in range(0, nz):
            rho_env[:,:,i] = rho_env2d
        # create the array of density of disk and the whole structure
        #
        rho      = np.zeros([len(rc),len(thetac),len(phic)])
        # The function for calculating the normalization of disk using the total disk mass
        #
        for ir in range(0,len(rc)):
            for itheta in range(0,len(thetac)):
                for iphi in range(0,len(phic)):
                    if rc[ir] > rin_shell:
                        # Envelope profile
                        rho[ir,itheta,iphi] = rho_env[ir,itheta,iphi]
                    else:
                        rho[ir,itheta,iphi] = 1e-25
        rho_env  = rho_env  + 1e-40
        rho      = rho      + 1e-40

    # Call function to plot the density
    plot_density(rho, rc, thetac,'/Users/yaolun/bhr71/hyperion/', plotname='shell')
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

    # Setting up the wavelength for monochromatic radiative transfer
    lambda0 = 0.1
    lambda1 = 2.0
    lambda2 = 50.0
    lambda3 = 95.0
    lambda4 = 200.0
    lambda5 = 314.0
    lambda6 = 670.0
    n01     = 10.0
    n12     = 20.0
    n23     = (lambda3-lambda2)/0.02
    n34     = (lambda4-lambda3)/0.03
    n45     = (lambda5-lambda4)/0.1
    n56     = (lambda6-lambda5)/0.1

    lam01   = lambda0 * (lambda1/lambda0)**(np.arange(n01)/n01)
    lam12   = lambda1 * (lambda2/lambda1)**(np.arange(n12)/n12)
    lam23   = lambda2 * (lambda3/lambda2)**(np.arange(n23)/n23)
    lam34   = lambda3 * (lambda4/lambda3)**(np.arange(n34)/n34)
    lam45   = lambda4 * (lambda5/lambda4)**(np.arange(n45)/n45)
    lam56   = lambda5 * (lambda6/lambda5)**(np.arange(n56+1)/n56)

    lam     = np.concatenate([lam01,lam12,lam23,lam34,lam45,lam56])
    nlam    = len(lam)

    # Create camera wavelength points
    n12     = 70.0
    n23     = 70.0
    n34     = 70.0
    n45     = 50.0
    n56     = 50.0
    
    lam12   = lambda1 * (lambda2/lambda1)**(np.arange(n12)/n12)
    lam23   = lambda2 * (lambda3/lambda2)**(np.arange(n23)/n23)
    lam34   = lambda3 * (lambda4/lambda3)**(np.arange(n34)/n34)
    lam45   = lambda4 * (lambda5/lambda4)**(np.arange(n45)/n45)
    lam56   = lambda5 * (lambda6/lambda5)**(np.arange(n56+1)/n56)

    lam_cam = np.concatenate([lam12,lam23,lam34,lam45,lam56])
    n_lam_cam = len(lam_cam)

    # Radiative transfer setting

    # number of photons for temp and image
    m.set_raytracing(True)
    m.set_monochromatic(True, wavelengths=[3.6, 4.5, 5.8, 8.0, 24, 70, 100, 160, 250, 350, 500])
    m.set_n_photons(initial=1000000, imaging_sources=1000000, imaging_dust=1000000,raytracing_sources=1000000, raytracing_dust=1000000)
    # imaging=100000, raytracing_sources=100000, raytracing_dust=100000
    # number of iteration to compute dust specific energy (temperature)
    m.set_n_initial_iterations(5)
    m.set_convergence(True, percentile=99., absolute=1.5, relative=1.02)
    m.set_mrw(True)   # Gamma = 1 by default
    # m.set_forced_first_scattering(forced_first_scattering=True)
    # Setting up images and SEDs
    image = m.add_peeled_images()
    # image.set_wavelength_range(300, 2.0, 670.0)
    # use the index of wavelength array used by the monochromatic radiative transfer
    image.set_wavelength_index_range(2,12)
    # pixel number
    image.set_image_size(300, 300)
    image.set_image_limits(-R_env_max, R_env_max, -R_env_max, R_env_max)
    image.set_viewing_angles([82.0], [0.0])
    image.set_uncertainties(True)
    # output as 64-bit
    image.set_output_bytes(8)

    # Output setting
    # Density
    m.conf.output.output_density = 'last'

    # Density difference (shows where dust was destroyed)
    m.conf.output.output_density_diff = 'none'

    # Energy absorbed (using pathlengths)
    m.conf.output.output_specific_energy = 'last'

    # Number of unique photons that passed through the cell
    m.conf.output.output_n_photons = 'last'

    m.write(outdir+outname+'.rtin')






indir = '/Users/yaolun/bhr71/radmc3d_params'
outdir = '/Users/yaolun/bhr71/hyperion/'
setup_model_shell(indir,outdir,'tsc_shell',low_res=True,scale=1,idl=True)