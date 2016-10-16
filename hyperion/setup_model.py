def setup_model(outdir,record_dir,outname,params,dust_file,tsc=True,idl=False,plot=False,\
                low_res=True,flat=True,scale=1,radmc=False,mono=False,mono_wave=None,
                record=True,dstar=200.,aperture=None,dyn_cav=False,fix_params=None,
                power=2,better_im=False,ellipsoid=False,TSC_dir='~/programs/misc/TSC/',
                IDL_path='/Applications/exelis/idl83/bin/idl',auto_disk=0.25,fast_plot=False,
                image_only=False, tsc_com=False, ext_source=None):
    """
    params = dictionary of the model parameters
    'alma' keyword is obsoleted
    outdir: The directory for storing Hyperion input files
    record_dir: The directory contains "model_list.txt" for recording parameters
    TSC_dir: Path the TSC-related IDL routines
    IDL_path: The IDL executable
    fast_plot: Do not plot the polar plot of the density because the rendering
               takes quite a lot of time.
    mono: monochromatic radiative transfer mode (need to specify the wavelength
          or a list of wavelength with 'mono_wave')
    image_only: only run for images
    """
    import numpy as np
    import astropy.constants as const
    from astropy.io import ascii
    import scipy as sci
    # to avoid X server error
    import matplotlib as mpl
    mpl.use('Agg')
    #
    import matplotlib.pyplot as plt
    import os
    from matplotlib.colors import LogNorm
    from scipy.integrate import nquad
    from hyperion.model import Model
    from record_hyperion import record_hyperion
    from outflow_inner_edge import outflow_inner_edge
    from pprint import pprint

    # Constants setup
    c         = const.c.cgs.value
    AU        = const.au.cgs.value     # Astronomical Unit       [cm]
    pc        = const.pc.cgs.value     # Parsec                  [cm]
    MS        = const.M_sun.cgs.value  # Solar mass              [g]
    LS        = const.L_sun.cgs.value  # Solar luminosity        [erg/s]
    RS        = const.R_sun.cgs.value  # Solar radius            [cm]
    G         = const.G.cgs.value      # Gravitational constant  [cm3/g/s^2]
    yr        = 60*60*24*365           # Years in seconds
    PI        = np.pi                  # PI constant
    sigma     = const.sigma_sb.cgs.value  # Stefan-Boltzmann constant
    mh        = const.m_p.cgs.value + const.m_e.cgs.value
    g2d       = 100.
    mmw       = 2.37                   # Kauffmann 2008

    m = Model()

    # min and max wavelength to compute (need to define them first for checking dust properties)
    # !!!
    wav_min = 2.0
    wav_max = 1400.
    wav_num = 1400

    # Create dust properties
    # Hyperion needs nu, albedo, chi, g, p_lin_max
    from hyperion.dust import HenyeyGreensteinDust
    # Read in the dust opacity table used by RADMC-3D
    dust = dict()
    [dust['nu'], dust['albedo'], dust['chi'], dust['g']] = np.genfromtxt(dust_file).T
    d = HenyeyGreensteinDust(dust['nu'], dust['albedo'], dust['chi'], dust['g'], dust['g']*0)
    # dust sublimation option
    d.set_sublimation_temperature('slow', temperature=1600.0)
    d.set_lte_emissivities(n_temp=3000,
                           temp_min=0.1,
                           temp_max=2000.)
    # if the min and/or max wavelength fall out of range
    if c/wav_min/1e-4 > dust['nu'].max():
        d.optical_properties.extrapolate_nu(dust['nu'].min(), c/wav_min/1e-4)
        print 'minimum wavelength is out of dust model.  The dust model is extrapolated.'
    if c/wav_max/1e-4 < dust['nu'].min():
        d.optical_properties.extrapolate_nu(c/wav_max/1e-4, dust['nu'].max())
        print 'maximum wavelength is out of dust model.  The dust model is extrapolated.'

    # try to solve the freq. problem
    d.optical_properties.extrapolate_nu(3.28e15, 5e15)
    #
    d.write(outdir+os.path.basename(dust_file).split('.')[0]+'.hdf5')
    d.plot(outdir+os.path.basename(dust_file).split('.')[0]+'.png')
    plt.clf()

    # Grids and Density

    # Grid Parameters
    nx        = 300L
    if low_res == True:
        nx    = 100L
    ny        = 400L
    nz        = 50L
    [nx, ny, nz] = [int(scale*nx), int(scale*ny), int(scale*nz)]

    # TSC model input setting
    dict_params = params
    # TSC model parameter
    cs        = dict_params['Cs']*1e5
    t         = dict_params['age']  # year
    omega     = dict_params['Omega0']
    # calculate related parameters
    M_env_dot = 0.975*cs**3/G
    mstar     = M_env_dot * t * yr
    R_cen     = omega**2 * G**3 * mstar**3 /(16*cs**8)
    R_inf     = cs * t * yr
    # protostar parameter
    tstar     = dict_params['tstar']
    R_env_max = dict_params['R_env_max']*AU
    theta_cav = dict_params['theta_cav']
    rho_cav_center = dict_params['rho_cav_center']
    rho_cav_edge   = dict_params['rho_cav_edge']*AU
    rstar     = dict_params['rstar']*RS
    # Mostly fixed parameter
    M_disk    = dict_params['M_disk']*MS
    beta      = dict_params['beta']
    h100      = dict_params['h100']*AU
    rho_cav   = dict_params['rho_cav']
    # make M_disk varies with mstar, which is the mass of star+disk
    if auto_disk != None:
        if M_disk != 0:
            print 'M_disk is reset to %4f of mstar (star+disk)' % auto_disk
            M_disk = mstar * auto_disk
        else:
            print 'M_disk = 0 is found.  M_disk is set to 0.'

    # ellipsoid cavity parameter
    if ellipsoid == True:
        # the numbers are given in arcsec
        a_out = 130 * dstar * AU
        b_out = 50  * dstar * AU
        z_out = a_out
        a_in  = dict_params['a_in'] * dstar * AU
        b_in  = a_in/a_out*b_out
        z_in  = a_in
        rho_cav_out = dict_params['rho_cav_out'] * mh
        rho_cav_in  = dict_params['rho_cav_in']  * mh
    # Calculate the dust sublimation radius
    T_sub = 1600
    a     = 1   # in micron
    # realistic dust
    # d_sub = 2.9388e7*(a/0.1)**-0.2 * (4*np.pi*rstar**2*sigma*tstar**4/LS)**0.5 / T_sub**3 *AU
    # black body dust
    d_sub = (LS/16./np.pi/sigma/AU**2*(4*np.pi*rstar**2*sigma*tstar**4/LS)/T_sub**4)**0.5 *AU
    # use the dust sublimation radius as the inner radius of disk and envelope
    R_disk_min = d_sub
    R_env_min  = d_sub
    rin        = rstar
    rout       = R_env_max
    R_disk_max = R_cen

    # print the variables
    print 'Dust sublimation radius %6f AU' % (d_sub/AU)
    print 'M_star %4f Solar mass' % (mstar/MS)
    print 'Infall radius %4f AU' % (R_inf / AU)

    # if there is any parameter found in fix_params, then fix them
    if fix_params != None:
        if 'R_min' in fix_params.keys():
            R_disk_min = fix_params['R_min']*AU
            R_env_min  = fix_params['R_min']*AU

    # Make the Coordinates
    #
    ri           = rin * (rout/rin)**(np.arange(nx+1).astype(dtype='float')/float(nx))
    ri           = np.hstack((0.0, ri))
    if ext_source:
        ri = np.hstack((ri, ri.max()+AU))
    thetai       = PI*np.arange(ny+1).astype(dtype='float')/float(ny)
    phii         = PI*2.0*np.arange(nz+1).astype(dtype='float')/float(nz)

    # Keep the constant cell size in r-direction at large radii
    #
    if flat == True:
        ri_cellsize = ri[1:-1]-ri[0:-2]
        ind = np.where(ri_cellsize/AU > 100.0)[0][0]       # The largest cell size is 100 AU
        ri = np.hstack((ri[0:ind],ri[ind]+np.arange(np.ceil((rout-ri[ind])/100/AU))*100*AU))
        nxx = nx
        nx = len(ri)-1
    # Assign the coordinates of the center of cell as its coordinates.
    #
    rc           = 0.5*( ri[0:nx]     + ri[1:nx+1] )
    thetac       = 0.5*( thetai[0:ny] + thetai[1:ny+1] )
    phic         = 0.5*( phii[0:nz]   + phii[1:nz+1] )

    # for non-TSC model
    if tsc_com:
        import hyperion as hp
        from hyperion.model import AnalyticalYSOModel

        non_tsc = AnalyticalYSOModel()

        # Define the luminsoity source
        nt_source = non_tsc.add_spherical_source()
        nt_source.luminosity = (4*PI*rstar**2)*sigma*(tstar**4)  # [ergs/s]
        nt_source.radius = rstar  # [cm]
        nt_source.temperature = tstar  # [K]
        nt_source.position = (0., 0., 0.)
        nt_source.mass = mstar

        # Envelope structure
        #
        nt_envelope = non_tsc.add_ulrich_envelope()
        nt_envelope.mdot = M_env_dot    # Infall rate
        nt_envelope.rmin = rin          # Inner radius
        nt_envelope.rc   = R_cen        # Centrifugal radius
        nt_envelope.rmax = R_env_max    # Outer radius
        nt_envelope.star = nt_source

        nt_grid = hp.grid.SphericalPolarGrid(ri, thetai, phii)

        rho_env_ulrich = nt_envelope.density(nt_grid).T
        rho_env_ulrich2d = np.sum(rho_env_ulrich**2,axis=2)/np.sum(rho_env_ulrich,axis=2)

    # Make the dust density model
    # Make the density profile of the envelope
    #
    total_mass = 0
    if tsc == False:
        print 'Calculating the dust density profile with infall solution...'
        if theta_cav != 0:
            # using R = 10000 AU as the reference point
            c0 = (10000.*AU)**(-0.5)*np.sqrt(1/np.sin(np.radians(theta_cav))**3-1/np.sin(np.radians(theta_cav)))
        else:
            c0 = 0
        rho_env  = np.zeros([len(rc),len(thetac),len(phic)])
        rho_disk = np.zeros([len(rc),len(thetac),len(phic)])
        rho      = np.zeros([len(rc),len(thetac),len(phic)])

        if dyn_cav == True:
            print 'WARNING: Calculation of interdependent cavity property has not implemented in infall-only solution!'
        # Normalization for the total disk mass
        def f(w,z,beta,rstar,h100):
            f = 2*PI*w*(1-np.sqrt(rstar/w))*(rstar/w)**(beta+1)*np.exp(-0.5*(z/(w**beta*h100/100**beta))**2)
            return f

        rho_0 = M_disk/(nquad(f,[[R_disk_min,R_disk_max],[-R_env_max,R_env_max]], args=(beta,rstar,h100)))[0]
        i = 0
        j = 0
        if 'rho_cav_center' in locals() == False:
            rho_cav_center = 5e-19
            print 'Use 5.27e-18 as the default value for cavity center'
        if 'rho_cav_edge' in locals() == False:
            rho_cav_edge = 40*AU
            print 'Use 40 AU as the default value for size of the inner region'
        discont = 1
        for ir in range(0,len(rc)):
            for itheta in range(0,len(thetac)):
                for iphi in range(0,len(phic)):
                    if rc[ir] > R_env_min:
                        # Envelope profile
                        w = abs(rc[ir]*np.cos(np.pi/2 - thetac[itheta]))
                        z = rc[ir]*np.sin(np.pi/2 - thetac[itheta])

                        if ellipsoid == False:
                            z_cav = c0*abs(w)**1.5
                            if z_cav == 0:
                                z_cav = R_env_max
                            cav_con = abs(z) > abs(z_cav)
                            if theta_cav == 90:
                                cav_con = True
                        else:
                            # condition for the outer ellipsoid
                            cav_con = (2*(w/b_out)**2 + ((abs(z)-z_out)/a_out)**2) < 1
                        if cav_con:
                            # open cavity
                            if ellipsoid == False:
                                if rho_cav_edge == 0:
                                    rho_cav_edge = R_env_min
                                if (rc[ir] <= rho_cav_edge) & (rc[ir] >= R_env_min):
                                    rho_env[ir,itheta,iphi] = g2d * rho_cav_center
                                else:
                                    rho_env[ir,itheta,iphi] = g2d * rho_cav_center*discont*(rho_cav_edge/rc[ir])**power
                                i += 1
                            else:
                                # condition for the inner ellipsoid
                                if (2*(w/b_in)**2 + ((abs(z)-z_in)/a_in)**2) > 1:
                                    rho_env[ir,itheta,iphi] = rho_cav_out
                                else:
                                    rho_env[ir,itheta,iphi] = rho_cav_in
                                i +=1
                        else:
                            j += 1
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
                            rho_env[ir,itheta,iphi] = M_env_dot/(4*PI*(G*mstar*R_cen**3)**0.5)*(rc[ir]/R_cen)**(-3./2)*(1+mu/mu_o)**(-0.5)*(mu/mu_o+2*mu_o**2*R_cen/rc[ir])**(-1)
                        # Disk profile
                        if ((w >= R_disk_min) and (w <= R_disk_max)) == True:
                            h = ((w/(100*AU))**beta)*h100
                            rho_disk[ir,itheta,iphi] = rho_0*(1-np.sqrt(rstar/w))*(rstar/w)**(beta+1)*np.exp(-0.5*(z/h)**2)
                        # Combine envelope and disk
                        rho[ir,itheta,iphi] = rho_disk[ir,itheta,iphi] + rho_env[ir,itheta,iphi]
                    else:
                        rho[ir,itheta,iphi] = 1e-30
                    # add the dust mass into the total count
                    cell_mass = rho[ir, itheta, iphi] * (1/3.)*(ri[ir+1]**3 - ri[ir]**3) * (phii[iphi+1]-phii[iphi]) * -(np.cos(thetai[itheta+1])-np.cos(thetai[itheta]))
                    total_mass = total_mass + cell_mass

        rho_env  = rho_env  + 1e-40
        rho_disk = rho_disk + 1e-40
        rho      = rho      + 1e-40
    # TSC model
    else:
        print 'Calculating the dust density profile with TSC solution...'
        if theta_cav != 0:
            c0 = (1e4*AU)**(-0.5)*np.sqrt(1/np.sin(np.radians(theta_cav))**3-1/np.sin(np.radians(theta_cav)))
        else:
            c0 = 0
        # If needed, calculate the TSC model via IDL
        #
        if idl == True:
            print 'Using IDL to calculate the TSC model.  Make sure you are running this on mechine with IDL.'
            import pidly
            idl = pidly.IDL(IDL_path)
            idl('.r '+TSC_dir+'tsc.pro')
            idl('.r '+TSC_dir+'tsc_run.pro')
            #
            # only run TSC calculation within infall radius
            # modify the rc array
            ind_infall = np.where(rc >= R_inf)[0][0]
            if max(ri) > R_inf:
                rc_idl = rc[0:ind_infall+1]
            else:
                rc_idl = rc[rc < max(ri)]
            idl.pro('tsc_run', indir=TSC_dir, outdir=outdir, rc=rc_idl, thetac=thetac, time=t,
                    c_s=cs, omega=omega, renv_min=R_env_min)
            file_idl = 'rhoenv.dat'
        else:
            print 'Read the pre-computed TSC model.'
            ind_infall = np.where(rc >= R_inf)[0][0]
            if max(ri) > R_inf:
                rc_idl = rc[0:ind_infall+1]
            else:
                rc_idl = rc[rc < max(ri)]
            if idl != False:
                file_idl = idl

        # read in the exist file
        rho_env_tsc_idl = np.genfromtxt(outdir+file_idl).T
        # because only region within infall radius is calculated by IDL program,
        # need to project it to the original grid
        rho_env_tsc = np.zeros([len(rc), len(thetac)])
        for irc in range(len(rc)):
            if rc[irc] in rc_idl:
                rho_env_tsc[irc,:] = rho_env_tsc_idl[np.squeeze(np.where(rc_idl == rc[irc])),:]

        # extrapolate for the NaN values at the outer radius, usually at radius beyond the infall radius
        # using r^-2 profile at radius greater than infall radius
        # and map the 2d strcuture onto 3-D grid
        # map TSC solution from IDL to actual 2-D grid
        rho_env_tsc2d = np.empty((nx,ny))
        if max(ri) > R_inf:
            for i in range(0, len(rc)):
                if i <= ind_infall:
                    rho_env_tsc2d[i,:] = rho_env_tsc[i,:]
                else:
                    rho_env_tsc2d[i,:] =  10**(np.log10(rho_env_tsc[ind_infall,:]) - 2*(np.log10(rc[i]/rc[ind_infall])))
        else:
            rho_env_tsc2d = rho_env_tsc
        # map it to 3-D grid
        rho_env = np.empty((nx,ny,nz))
        for i in range(0, nz):
            rho_env[:,:,i] = rho_env_tsc2d

        # typical no used.  Just an approach I tried to make the size of the
        # constant desnity region self-consistent with the outflow cavity.
        if dyn_cav == True:
            print 'Calculate the cavity properties using the criteria that swept-up mass = outflowed mass'
            # using swept-up mass = flow mass to derive the edge of the extended flat density region
            v_outflow = 1e2 * 1e5
            rho_cav_edge = outflow_inner_edge(np.copy(rho_env), (ri,thetai,phii),M_env_dot,v_outflow,theta_cav, R_env_min)
            dict_params['rho_cav_edge'] = rho_cav_edge
            # assume gas-to-dust ratio = 100
            rho_cav_center = 0.01 * 0.1*M_env_dot*rho_cav_edge/v_outflow/2 / (2*np.pi/3*rho_cav_edge**3*(1-np.cos(np.radians(theta_cav))))
            dict_params['rho_cav_center'] = rho_cav_center
            print 'inner edge is %5f AU and density is %e g/cm3' % (rho_cav_edge/AU, rho_cav_center)

        # create the array of density of disk and the whole structure
        #
        rho_disk = np.zeros([len(rc),len(thetac),len(phic)])
        rho      = np.zeros([len(rc),len(thetac),len(phic)])

        # non-TSC option
        if tsc_com:
            rho_ulrich = np.zeros([len(rc),len(thetac),len(phic)])

        # Calculate the disk scale height by the normalization of h100
        def f(w,z,beta,rstar,h100):
            f = 2*PI*w*(1-np.sqrt(rstar/w))*(rstar/w)**(beta+1)*np.exp(-0.5*(z/(w**beta*h100/100**beta))**2)
            return f
        # The function for calculating the normalization of disk using the total disk mass
        #
        rho_0 = M_disk/(nquad(f,[[R_disk_min,R_disk_max],[-R_env_max,R_env_max]], args=(beta,rstar,h100)))[0]
        i = 0
        j = 0

        # put in default outflow cavity setting if nothing is specified
        if 'rho_cav_center' in locals() == False:
            rho_cav_center = 5e-19
            print 'Use 5e-19 as the default value for cavity center'
        if 'rho_cav_edge' in locals() == False:
            rho_cav_edge = 40*AU
            print 'Use 40 AU as the default value for size of the inner region'
        discont = 1
        for ir in range(0,len(rc)):
            for itheta in range(0,len(thetac)):
                for iphi in range(0,len(phic)):
                    if rc[ir] > R_env_min:
                        # Envelope profile
                        w = abs(rc[ir]*np.cos(np.pi/2 - thetac[itheta]))
                        z = rc[ir]*np.sin(np.pi/2 - thetac[itheta])

                        if ellipsoid == False:
                            z_cav = c0*abs(w)**1.5
                            if z_cav == 0:
                                z_cav = R_env_max
                            cav_con = abs(z) > abs(z_cav)
                        else:
                            # condition for the outer ellipsoid
                            cav_con = (2*(w/b_out)**2 + ((abs(z)-z_out)/a_out)**2) < 1
                        if cav_con:
                            # open cavity
                            if ellipsoid == False:
                                if rho_cav_edge == 0:
                                    rho_cav_edge = R_env_min
                                if (rc[ir] <= rho_cav_edge) & (rc[ir] >= R_env_min):
                                    rho_env[ir,itheta,iphi] = g2d * rho_cav_center#*((rc[ir]/AU)**2)
                                    if tsc_com:
                                        rho_env_ulrich[ir,itheta,iphi] = rho_env[ir,itheta,iphi]
                                else:
                                    rho_env[ir,itheta,iphi] = g2d * rho_cav_center*discont*(rho_cav_edge/rc[ir])**power
                                    if tsc_com:
                                        rho_env_ulrich[ir,itheta,iphi] = rho_env[ir,itheta,iphi]
                                i += 1
                            else:
                                # condition for the inner ellipsoid
                                if (2*(w/b_in)**2 + ((abs(z)-z_in)/a_in)**2) > 1:
                                    rho_env[ir,itheta,iphi] = rho_cav_out
                                    if tsc_com:
                                        rho_env_ulrich[ir,itheta,iphi] = rho_env[ir,itheta,iphi]
                                else:
                                    rho_env[ir,itheta,iphi] = rho_cav_in
                                    if tsc_com:
                                        rho_env_ulrich[ir,itheta,iphi] = rho_env[ir,itheta,iphi]
                                i +=1

                        # Disk profile
                        if ((w >= R_disk_min) and (w <= R_disk_max)) == True:
                            h = ((w/(100*AU))**beta)*h100
                            rho_disk[ir,itheta,iphi] = rho_0*(1-np.sqrt(rstar/w))*(rstar/w)**(beta+1)*np.exp(-0.5*(z/h)**2)
                        # Combine envelope and disk
                        rho[ir,itheta,iphi] = rho_disk[ir,itheta,iphi] + rho_env[ir,itheta,iphi]
                        if tsc_com:
                            rho_ulrich[ir,itheta,iphi] = rho_disk[ir,itheta,iphi] + rho_env_ulrich[ir,itheta,iphi]
                    else:
                        rho[ir,itheta,iphi] = 1e-40
                        if tsc_com:
                            rho[ir,itheta,iphi] = 1e-40
                    # add the dust mass into the total count
                    cell_mass = rho[ir, itheta, iphi] * (1/3.)*(ri[ir+1]**3 - ri[ir]**3) * (phii[iphi+1]-phii[iphi]) * -(np.cos(thetai[itheta+1])-np.cos(thetai[itheta]))
                    total_mass = total_mass + cell_mass
    # apply gas-to-dust ratio of 100
    rho_dust = rho/g2d
    if tsc_com:
        rho_ulrich_dust = rho_ulrich/g2d
    total_mass_dust = total_mass/MS/g2d
    print 'Total dust mass = %f Solar mass' % total_mass_dust

    if record == True:
        # Record the input and calculated parameters
        params = dict_params.copy()
        params.update({'d_sub': d_sub/AU, 'M_env_dot': M_env_dot/MS*yr, 'R_inf': R_inf/AU, 'R_cen': R_cen/AU, 'mstar': mstar/MS, 'M_tot_gas': total_mass/MS})
        record_hyperion(params,record_dir)

    if plot == True:
        # rho2d is the 2-D projection of gas density
        # take the weighted average
        rho2d = np.sum(rho**2,axis=2)/np.sum(rho,axis=2)

        if tsc_com:
            rho2d = np.sum(rho_ulrich**2,axis=2)/np.sum(rho_ulrich,axis=2)

        if fast_plot == False:
            # Plot the azimuthal averaged density
            fig = plt.figure(figsize=(8,6))
            ax_env  = fig.add_subplot(111,projection='polar')

            zmin = 1e-22/mmw/mh
            zmin = 1e-1
            cmap = plt.cm.CMRmap
            rho2d_exp = np.hstack((rho2d,rho2d,rho2d[:,0:1]))
            thetac_exp = np.hstack((thetac-PI/2, thetac+PI/2, thetac[0]-PI/2))
            # plot the gas density
            img_env = ax_env.pcolormesh(thetac_exp,rc/AU,rho2d_exp/mmw/mh,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=1e6)) # np.nanmax(rho2d_exp/mmw/mh)

            ax_env.set_xlabel(r'$\rm{Polar\,angle\,(Degree)}$',fontsize=20)
            ax_env.set_ylabel('',fontsize=20, labelpad=-140)
            ax_env.tick_params(labelsize=18)
            ax_env.set_yticks(np.hstack((np.arange(0,(int(R_env_max/AU/10000.)+1)*10000, 10000),R_env_max/AU)))
            ax_env.set_xticklabels([r'$\rm{90^{\circ}}$',r'$\rm{45^{\circ}}$',r'$\rm{0^{\circ}}$',r'$\rm{-45^{\circ}}$',\
                                    r'$\rm{-90^{\circ}}$',r'$\rm{-135^{\circ}}$',r'$\rm{180^{\circ}}$',r'$\rm{135^{\circ}}$'])
            ax_env.set_yticklabels([])
            # fix the tick label font
            ticks_font = mpl.font_manager.FontProperties(family='STIXGeneral',size=20)
            for label in ax_env.get_yticklabels():
                label.set_fontproperties(ticks_font)

            ax_env.grid(True, color='LightGray', linewidth=1)
            cb = fig.colorbar(img_env, pad=0.1)
            cb.ax.set_ylabel(r'$\rm{Averaged\,Gas\,Density\,(cm^{-3})}$',fontsize=20)
            # cb.set_ticks([1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9])
            # cb.set_ticklabels([r'$\rm{10^{2}}$',r'$\rm{10^{3}}$',r'$\rm{10^{4}}$',r'$\rm{10^{5}}$',r'$\rm{10^{6}}$',\
            #                    r'$\rm{10^{7}}$',r'$\rm{10^{8}}$',r'$\rm{\geq 10^{9}}$'])
            # lower density ticks
            cb.set_ticks([1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6])
            cb.set_ticklabels([r'$\rm{10^{-1}}$',r'$\rm{10^{0}}$',r'$\rm{10^{1}}$',r'$\rm{10^{2}}$',r'$\rm{10^{3}}$',
                               r'$\rm{10^{4}}$',r'$\rm{10^{5}}$',r'$\rm{\geq 10^{6}}$'])

            cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
            plt.setp(cb_obj,fontsize=20)
            fig.savefig(outdir+outname+'_gas_density.png', format='png', dpi=300, bbox_inches='tight')
            fig.clf()

        # Plot the radial density profile
        fig = plt.figure(figsize=(12,9))
        ax = fig.add_subplot(111)

        plot_grid = [0,49,99,149,199]
        color_grid = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']
        label = [r'$\rm{\theta='+str(int(np.degrees(thetai[plot_grid[0]])))+'^{\circ}}$',\
                 r'$\rm{\theta='+str(int(np.degrees(thetai[plot_grid[1]])))+'^{\circ}}$',\
                 r'$\rm{\theta='+str(1+int(np.degrees(thetai[plot_grid[2]])))+'^{\circ}}$',\
                 r'$\rm{\theta='+str(int(np.degrees(thetai[plot_grid[3]])))+'^{\circ}}$',\
                 r'$\rm{\theta='+str(1+int(np.degrees(thetai[plot_grid[4]])))+'^{\circ}}$']
        alpha = np.linspace(0.3,1.0,len(plot_grid))
        for i in plot_grid:
            ax.plot(np.log10(rc[rc > 0.14*AU]/AU), np.log10(rho2d[rc > 0.14*AU,i]/g2d/mmw/mh)+plot_grid[::-1].index(i)*-0.2,'-',color=color_grid[plot_grid.index(i)],mec='None',linewidth=2.5, \
                    markersize=3, label=label[plot_grid.index(i)]) # alpha=alpha[plot_grid.index(i)],
        ax.axvline(np.log10(R_inf/AU), linestyle='--', color='k', linewidth=1.5, label=r'$\rm{infall\,radius}$')
        ax.axvline(np.log10(R_cen/AU), linestyle=':', color='k', linewidth=1.5, label=r'$\rm{centrifugal\,radius}$')

        lg = plt.legend(fontsize=20, numpoints=1, ncol=2, framealpha=0.7, loc='upper right')

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
        fig.savefig(outdir+outname+'_gas_radial.pdf',format='pdf',dpi=300,bbox_inches='tight')
        fig.clf()

    # Insert the calculated grid and dust density profile into hyperion
    m.set_spherical_polar_grid(ri, thetai, phii)

    m.add_density_grid(rho_dust.T, d)
    # for non-TSC option
    if tsc_com:
        m.add_density_grid(rho_ulrich_dust.T, d)

    # Define the luminsoity source
    source = m.add_spherical_source()
    source.luminosity = (4*PI*rstar**2)*sigma*(tstar**4)  # [ergs/s]
    source.radius = rstar  # [cm]
    source.temperature = tstar  # [K]
    source.position = (0., 0., 0.)
    print 'L_center =  % 5.2f L_sun' % ((4*PI*rstar**2)*sigma*(tstar**4)/LS)

    if ext_source != None:
        # add external heating - ISRF
        # use standard receipe from Hyperion doc
        isrf = ascii.read(ext_source, names=['wavelength', 'J_lambda'])
        isrf_nu = c/(isrf['wavelength']*1e-4)
        isrf_jnu = isrf['J_lambda']*isrf['wavelength']/isrf_nu

        if 'mmp83' in ext_source:
            FOUR_PI_JNU = 0.0217
        else:
            FOUR_PI_JNU = raw_input('What is the FOUR_PI_JNU value?')

        s_isrf = m.add_external_spherical_source()
        s_isrf.radius = ri.max()
        s_isrf.spectrum = (isrf_nu, isrf_jnu)
        s_isrf.luminosity = PI * s_isrf.radius**2 * FOUR_PI_JNU

    m.set_raytracing(True)
    # option of using more photons for imaging
    if better_im == False:
        im_photon = 1e6
    else:
        im_photon = 5e7

    if mono == True:
        if (type(mono_wave) == int) or (type(mono_wave) == float) or (type(mono_wave) == str):
            mono_wave = float(mono_wave)
            mono_wave = [mono_wave]

        # Monochromatic radiative transfer setting
        m.set_monochromatic(True, wavelengths=mono_wave)
        m.set_n_photons(initial=1e6, imaging_sources=im_photon,
            imaging_dust=im_photon,raytracing_sources=1e6, raytracing_dust=1e6)
    else:
        # regular wavelength grid setting
        m.set_n_photons(initial=1e6, imaging=im_photon,raytracing_sources=1e6,
                        raytracing_dust=1e6)
    # number of iteration to compute dust specific energy (temperature)
    m.set_n_initial_iterations(20)
    m.set_convergence(True, percentile=dict_params['percentile'],
                            absolute=dict_params['absolute'],
                            relative=dict_params['relative'])
    m.set_mrw(True)   # Gamma = 1 by default

    # Setting up images and SEDs
    if not image_only:
        # SED setting
        # Infinite aperture
        syn_inf = m.add_peeled_images(image=False)
        # use the index of wavelength array used by the monochromatic radiative transfer
        if mono == False:
            syn_inf.set_wavelength_range(wav_num, wav_min, wav_max)
        syn_inf.set_viewing_angles([dict_params['view_angle']], [0.0])
        syn_inf.set_uncertainties(True)
        syn_inf.set_output_bytes(8)

        # aperture
        # 7.2 in 10 um scaled by lambda / 10
        # flatten beyond 20 um
        # default aperture (should always specify a set of apertures)
        if aperture == None:
            aperture = {'wave': [3.6, 4.5, 5.8, 8.0, 8.5, 9, 9.7, 10, 10.5, 11, 16, 20, 24, 30, 70, 100, 160, 250, 350, 500, 1300],\
                        'aperture': [7.2, 7.2, 7.2, 7.2, 7.2, 7.2, 7.2, 7.2, 7.2, 7.2, 20.4, 20.4, 20.4, 20.4, 24.5, 24.5, 24.5, 24.5, 24.5, 24.5, 101]}
        # assign wl_aper and aper from dictionary of aperture
        wl_aper = aperture['wave']
        aper    = aperture['aperture']
        # create the non-repetitive aperture list and index array
        aper_reduced = list(set(aper))
        index_reduced = np.arange(1, len(aper_reduced)+1)

        dict_peel_sed = {}
        for i in range(0, len(aper_reduced)):
            aper_dum = aper_reduced[i]/2 * (1/3600.*np.pi/180.)*dstar*pc
            dict_peel_sed[str(index_reduced[i])] = m.add_peeled_images(image=False)
            # use the index of wavelength array used by the monochromatic radiative transfer
            if mono == False:
                dict_peel_sed[str(index_reduced[i])].set_wavelength_range(wav_num, wav_min, wav_max)
            dict_peel_sed[str(index_reduced[i])].set_viewing_angles([dict_params['view_angle']], [0.0])
            # aperture should be given in cm and its the radius of the aperture
            dict_peel_sed[str(index_reduced[i])].set_aperture_range(1, aper_dum, aper_dum)
            dict_peel_sed[str(index_reduced[i])].set_uncertainties(True)
            dict_peel_sed[str(index_reduced[i])].set_output_bytes(8)

    # image setting
    syn_im = m.add_peeled_images(sed=False)
    # use the index of wavelength array used by the monochromatic radiative transfer
    if mono == False:
        syn_im.set_wavelength_range(wav_num, wav_min, wav_max)
    # pixel number
    # !!!
    if not mono:
        pix_num = 300
    else:
        pix_num = 8000
    #
    syn_im.set_image_size(pix_num, pix_num)
    syn_im.set_image_limits(-R_env_max, R_env_max, -R_env_max, R_env_max)
    syn_im.set_viewing_angles([dict_params['view_angle']], [0.0])
    syn_im.set_uncertainties(True)
    syn_im.set_output_bytes(8)

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

    if radmc == True:
        # RADMC-3D still use a pre-defined aperture with lazy for-loop
        aper = np.zeros([len(lam)])
        ind = 0
        for wl in lam:
            if wl < 5:
                aper[ind] = 8.4
            elif wl >= 5 and wl < 14:
                aper[ind] = 1.8 * 4
            elif wl >= 14 and wl < 40:
                aper[ind] = 5.1 * 4
            else:
                aper[ind] = 24.5
            ind += 1

        # Write the wavelength_micron.inp file
        #
        f_wave = open(outdir+'wavelength_micron.inp','w')
        f_wave.write('%d \n' % int(nlam))
        for ilam in range(0,nlam):
            f_wave.write('%f \n' % lam[ilam])
        f_wave.close()

        # Write the camera_wavelength_micron.inp file
        #
        f_wave_cam = open(outdir+'camera_wavelength_micron.inp','w')
        f_wave_cam.write('%d \n' % int(nlam))
        for ilam in range(0,nlam):
            f_wave_cam.write('%f \n' % lam[ilam])
        f_wave_cam.close()

        # Write the aperture_info.inp
        #
        f_aper = open(outdir+'aperture_info.inp','w')
        f_aper.write('1 \n')
        f_aper.write('%d \n' % int(nlam))
        for iaper in range(0, len(aper)):
            f_aper.write('%f \t %f \n' % (lam[iaper],aper[iaper]/2))
        f_aper.close()

        # Write the stars.inp file
        #
        f_star = open(outdir+'stars.inp','w')
        f_star.write('2\n')
        f_star.write('1 \t %d \n' % int(nlam))
        f_star.write('\n')
        f_star.write('%e \t %e \t %e \t %e \t %e \n' % (rstar*0.9999,mstar,0,0,0))
        f_star.write('\n')
        for ilam in range(0,nlam):
            f_star.write('%f \n' % lam[ilam])
        f_star.write('\n')
        f_star.write('%f \n' % -tstar)
        f_star.close()

        # Write the grid file
        #
        f_grid = open(outdir+'amr_grid.inp','w')
        f_grid.write('1\n')                               # iformat
        f_grid.write('0\n')                               # AMR grid style  (0=regular grid, no AMR)
        f_grid.write('150\n')                             # Coordinate system  coordsystem<100: Cartisian; 100<=coordsystem<200: Spherical; 200<=coordsystem<300: Cylindrical
        f_grid.write('0\n')                               # gridinfo
        f_grid.write('1 \t 1 \t 1 \n')                    # Include x,y,z coordinate
        f_grid.write('%d \t %d \t %d \n' % (int(nx)-1,int(ny),int(nz)))    # Size of the grid
        [f_grid.write('%e \n' % ri[ir]) for ir in range(1,len(ri))]
        [f_grid.write('%f \n' % thetai[itheta]) for itheta in range(0,len(thetai))]
        [f_grid.write('%f \n' % phii[iphi]) for iphi in range(0,len(phii))]
        f_grid.close()

        # Write the density file
        #
        f_dust = open(outdir+'dust_density.inp','w')
        f_dust.write('1 \n')                      # format number
        f_dust.write('%d \n' % int((nx-1)*ny*nz)) # Nr of cells
        f_dust.write('1 \n')                      # Nr of dust species
        for iphi in range(0,len(phic)):
            for itheta in range(0,len(thetac)):
                for ir in range(1,len(rc)):
                    f_dust.write('%e \n' % rho_dust[ir,itheta,iphi])
        f_dust.close()

        # Write the dust opacity table
        f_dustkappa = open(outdir+'dustkappa_oh5_extended.inp','w')
        f_dustkappa.write('3 \n')                       # format index for including g-factor
        f_dustkappa.write('%d \n' % len(dust['nu']))    # number of wavlength/frequency in the table
        for i in range(len(dust['nu'])):
            f_dustkappa.write('%f \t %f \t %f \t %f \n' % (c/dust['nu'][i]*1e4, dust['chi'][i], dust['chi'][i]*dust['albedo'][i]/(1-dust['albedo'][i]), dust['g'][i]))
        f_dustkappa.close()

        # Write the Dust opacity control file
        #
        f_opac = open(outdir+'dustopac.inp','w')
        f_opac.write('2               Format number of this file\n')
        f_opac.write('1               Nr of dust species\n')
        f_opac.write('============================================================================\n')
        f_opac.write('1               Way in which this dust species is read\n')
        f_opac.write('0               0=Thermal grain\n')
        # f_opac.write('klaus           Extension of name of dustkappa_***.inp file\n')
        f_opac.write('oh5_extended    Extension of name of dustkappa_***.inp file\n')
        f_opac.write('----------------------------------------------------------------------------\n')
        f_opac.close()

        # Write the radmc3d.inp control file
        #
        f_control = open(outdir+'radmc3d.inp','w')
        f_control.write('nphot = %d \n' % 100000)
        f_control.write('scattering_mode_max = 2\n')
        f_control.write('camera_min_drr = 0.1\n')
        f_control.write('camera_min_dangle = 0.1\n')
        f_control.write('camera_spher_cavity_relres = 0.1\n')
        f_control.write('istar_sphere = 1\n')
        f_control.write('modified_random_walk = 1\n')
        f_control.close()

    return m

# from input_reader import input_reader_table
# from pprint import pprint
# import numpy as np
# filename = '/Users/yaolun/programs/misc/hyperion/test_input.txt'
# params = input_reader_table(filename)
# pprint(params[0])
# outdir = '/Users/yaolun/test/'
# record_dir = '/Users/yaolun/test/'
# dust_file = '/Users/yaolun/programs/misc/oh5_hyperion.txt'
# wl_aper, aper_arcsec = np.genfromtxt('/Users/yaolun/bhr71/best_calibrated/aperture.txt',
#                                      skip_header=1, dtype=float).T
# aperture = {'wave': wl_aper, 'aperture': aper_arcsec}
# setup_model(outdir,record_dir,'model_powerlaw',params[0],dust_file,plot=True,record=False,\
#     aperture=aperture,idl='model32_nontsc.dat',radmc=False,tsc=False)
