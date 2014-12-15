def problem_setup(indir,outdir,model=False,denser_wall=False,plot=False,low_res=False,flat=True):
    # coding: utf-8

    # In[ ]:

    # def a function with input of outdir / parameters set
    # Apply OH5 dust opacity table
    # Test for using sloppy mide within radmc3d.inp control file


    # In[145]:

    # get_ipython().magic(u'matplotlib inline')
    # %pylab inline --no-import-all
    import numpy as np
    import scipy as sci
    import matplotlib.pyplot as plt
    import matplotlib as mat
    import os
    from matplotlib.colors import LogNorm
    from scipy.optimize import fsolve
    from scipy.integrate import nquad
    from envelope_func import func


    # In[146]:

    # outdir = os.path.expanduser('~')+'/radmc_simulation/bhr71/'

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outdir = outdir+'/'

    # In[147]:

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


    # In[148]:

    # Parameters setup
    # Import the model parameters from another file 
    #
    # params     = np.genfromtxt(outdir+'params.dat',dtype=None)
    # params     = np.genfromtxt(indir+model+'.dat',dtype=None)
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

    # In[149]:

    # # Parameters Setup
    # #
    # tstar       = 4.218e3          # Stellar temperature         [K]
    # mstar       = 1.375            # Stellar mass                [Solar Mass]
    # rstar       = 8.426            # Stellar radius              [Solar radius]
    # M_env_dot   = MS*4.674e-4/yr   # Envelope accretion rate     [g/s]
    # M_disk_dot  = 1.040*1e-7*MS/yr # Disk accretion rate         [g/s]
    # R_env_max   = 17860*AU         # Envelope outer radius       [cm]
    # R_env_min   = 6.268*AU         # Envelope inner radius       [cm]
    # theta_cav   = 31.34            # Outflow cavity opening angle[deg]
    # R_disk_max  = 45.08*AU         # Disk outer radius           [cm]
    # R_disk_min  = 1.897*AU         # Disk inner radius           [cm]
    # R_cen       = R_disk_max       # Centrifugal radius          [cm]
    # M_disk      = 1.120*1e-2*MS    # Disk mass                   [g]
    # beta        = 1.154            # Disk flare factor           []
    # h100        = 6.854*AU         # Disk scale height at 100 AU [cm]
    # rho_cav     = 1.516e-20        # Outflow cavity density      [g/cm3]]


    # In[150]:

    # Monte Carlo Parameter
    #
    nphot     = 1000000

    # Grid Parameters
    #
    nx        = 300L
    if low_res == True:
        nx    = 100L
    ny        = 400L
    nz        = 50L

    # nx        = 20
    # ny        = 10
    # nz        = 5
    
    # Model Parameters
    #
    # rin       = R_env_min
    rin       = rstar
    rout      = R_env_max
    rcen      = R_cen

    # Star Parameters
    #
    mstar    = mstar
    rstar    = rstar*0.9999
    tstar    = tstar
    pstar    = [0.,0.,0.]


    # In[151]:

    # Make the Coordinates
    #
    ri           = rin * (rout/rin)**(np.arange(nx+1).astype(dtype='float')/float(nx))
    thetai       = PI*np.arange(ny+1).astype(dtype='float')/float(ny)
    phii         = PI*2.0*np.arange(nz+1).astype(dtype='float')/float(nz)
    # phii         = PI*2.0*np.arange(nz+1).astype(dtype='float')/float(nz+1)
    
    # Keep the constant cell size in r-direction
    #
    if flat == True:
        ri_cellsize = ri[1:-1]-ri[0:-2]
        ind = np.where(ri_cellsize/AU > 100.0)[0][0]       # The largest cell size is 100 AU
        ri = np.hstack((ri[0:ind],ri[ind]+np.arange(np.ceil((rout-ri[ind])/100/AU))*100*AU))
        nx = len(ri)-1    

    # In[152]:

    # print (ri[20]-ri[19])/ri[19]


    # In[153]:

    # Assign the coordinates of the center of cell as its coordinates.
    #
    rc           = 0.5*( ri[0:nx]     + ri[1:nx+1] )
    thetac       = 0.5*( thetai[0:ny] + thetai[1:ny+1] )
    phic         = 0.5*( phii[0:nz]   + phii[1:nz+1] )
    # phic         = 0.5*( phii[0:nz-1]   + phii[1:nz] )


    # In[155]:

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


    # In[156]:

    # print rho_disk.min(),rho_disk.max()

    if plot == True:

        # In[157]:

        # Plot the cavity wall profile
        # 
        print 'Plotting...'
        # print outdir
        # c = R_env_max**(-0.5)*np.sqrt(1/np.sin(np.radians(theta_cav))**3-1/np.sin(np.radians(theta_cav)))
        # cav_wall = []
        # for ir in range(0,len(rc)):
        #     cav_wall.append(c*rc[ir]**1.5)

        # cav_wall = np.array(cav_wall)
        # fig_cavwall = plt.figure(figsize=(9,6))
        # ax_cavwall  = fig_cavwall.add_subplot(111)
        # ax_cavwall.plot(rc/AU,cav_wall/AU,'x')
        # ax_cavwall.set_xlabel(r'$\mathrm{Radius~(AU)}$',fontsize=14)
        # ax_cavwall.set_ylabel(r'$\mathrm{Height~(AU)}$',fontsize=14)
        # ax_cavwall.set_xlim([R_env_min/AU,R_env_max/AU])
        # ax_cavwall.set_ylim([R_env_min/AU,R_env_max/AU])
        # fig_cavwall.savefig(outdir+'cavity_wall_profile.eps',format='eps',dpi=300)
        # ax_cavwall.cla()
        # fig_cavwall.clf()


        # In[162]:

        # fig_env = plt.figure()
        # # Plot the envelope density profile in r-theta coordinates
        
        # ax_env  = fig_env.add_subplot(122,projection='polar')

        # z = np.sum(rho_env**2,axis=2)/np.sum(rho_env,axis=2)

        # zmin = 1e-21
        # cmap = 'jet'
        # img_env = ax_env.pcolormesh(thetai,ri/AU,z,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=z.max()))
        # ax_env.pcolormesh(PI-thetai,ri/AU,z,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=z.max()))
        # ax_env.pcolormesh(-thetai,ri/AU,z,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=z.max()))
        # ax_env.pcolormesh(PI+thetai,ri/AU,z,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=z.max()))

        # ax_env.set_xlabel(r'$\mathrm{Polar~angle~(Degree)}$',fontsize=14)
        # ax_env.set_ylabel(r'$\mathrm{Radius~(AU)}$',fontsize=14)
        # ax_env.tick_params(labelsize=14)
        # ax_env.set_yticks(np.arange(0,R_env_max/AU,R_env_max/AU/5))
        # ax_env.grid(True)
        # cb = fig_env.colorbar(img_env)
        # cb.ax.set_ylabel(r'$\mathrm{Surface~Density~(g/cm^{2})}$',fontsize=14)
        # cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
        # plt.setp(cb_obj,fontsize=14)
        # Plot the mid plane slice of the envelope density profile in r-phi coordinates
        # 
        # ax_env_midplane = fig_env.add_subplot(111,projection='polar')

        # z = np.sum(rho_env**2,axis=1)/np.sum(rho_env,axis=1)

        # img_env = ax_env_midplane.pcolormesh(phii,ri/AU,z,cmap=cmap,  norm=LogNorm(vmin=zmin,vmax=z.max()))

        # ax_env_midplane.set_xlabel(r'$\mathrm{Azimuthal~angle~(Degree)}$',fontsize=14)
        # ax_env_midplane.set_ylabel(r'$\mathrm{Radius~(AU)}$',fontsize=14)
        # ax_env_midplane.tick_params(labelsize=14)
        # ax_env_midplane.set_yticks(np.arange(0,R_env_max/AU,R_env_max/AU/5))
        # ax_env_midplane.grid(True)
        # fig_env.savefig(outdir+'env_density_profile.pdf',format='pdf',dpi=300)
        # ax_env.cla()
        # ax_env_midplane.cla()
        # fig_env.clf()

        # # In[163]:

        # # Disk density profile
        # #
        # fig_disk = plt.figure(figsize=(20,7))
        # # Plot the disk density profile in r-theta coordinates
        # # 
        # ax_disk  = fig_disk.add_subplot(122,projection='polar')

        # z = np.sum(rho_disk**2,axis=2)/np.sum(rho_disk,axis=2)
        # print z.min(),z.max()
        # zmin = 1e-21
        # cmap = 'jet'
        # img_disk = ax_disk.pcolormesh(thetai,ri/AU,z,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=z.max()))
        # ax_disk.pcolormesh(PI-thetai,ri/AU,z,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=z.max()))
        # ax_disk.pcolormesh(-thetai,ri/AU,z,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=z.max()))
        # ax_disk.pcolormesh(PI+thetai,ri/AU,z,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=z.max()))

        # ax_disk.set_xlabel(r'$\mathrm{Polar~angle~(Degree)}$',fontsize=14)
        # ax_disk.set_ylabel(r'$\mathrm{Radius~(AU)}$',fontsize=14)
        # ax_disk.set_ylim([R_disk_min/AU,R_disk_max/AU])
        # ax_disk.tick_params(labelsize=14)
        # ax_disk.grid(True)
        # cb = fig_disk.colorbar(img_disk)
        # cb.ax.set_ylabel(r'$\mathrm{Surface~Density~(g/cm^{2})}$',fontsize=14)
        # cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
        # plt.setp(cb_obj,fontsize=14)
        # # Plot the mid plane slice of the disk density profile in r-phi coordinates
        # # 
        # ax_disk_midplane = fig_disk.add_subplot(121,projection='polar')

        # z = np.sum(rho_disk**2,axis=1)/np.sum(rho_disk,axis=1)

        # img_disk = ax_disk_midplane.pcolormesh(phii,ri/AU,z,cmap=cmap,  norm=LogNorm(vmin=zmin,vmax=z.max()))

        # ax_disk_midplane.set_xlabel(r'$\mathrm{Azimuthal~angle~(Degree)}$',fontsize=14)
        # ax_disk_midplane.set_ylabel(r'$\mathrm{Radius~(AU)}$',fontsize=14)
        # ax_disk_midplane.set_ylim([R_disk_min/AU,R_disk_max/AU])
        # ax_disk_midplane.tick_params(labelsize=14)
        # ax_disk_midplane.grid(True)
        # fig_disk.savefig(outdir+'disk_density_profile.eps',format='eps',dpi=300)
        # ax_disk.cla()
        # ax_disk_midplane.cla()
        # fig_disk.clf()

        # # In[164]:

        fig_all = plt.figure(figsize=(16,12))
        # Plot the density profile in r-theta coordinates
        # 
        ax_all  = fig_all.add_subplot(111,projection='polar')

        z = np.sum(rho**2,axis=2)/np.sum(rho,axis=2)

        zmin = 1e-22
        zmax = 1e-17# z.max()
        cmap = 'jet'
        img_1 = ax_all.pcolormesh(thetai,ri/AU,z,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=zmax))
        # img_2 = ax_all.pcolormesh(PI-thetai,ri/AU,z,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=z.max()))
        # img_3 = ax_all.pcolormesh(-thetai,ri/AU,z,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=z.max()))
        img_4 = ax_all.pcolormesh(PI+thetai,ri/AU,z,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=zmax))

        img_1.set_edgecolor('face')
        # img_2.set_edgecolor('face')
        # img_3.set_edgecolor('face')
        img_4.set_edgecolor('face')


        # Rasterize
        # img_1.set_rasterized(True)
        # img_2.set_rasterized(True)
        # img_3.set_rasterized(True)
        # img_4.set_rasterized(True)

        ax_all.set_xlabel(r'$\mathrm{Polar~angle~(Degree)}$',fontsize=32)
        ax_all.set_ylabel(r'$\mathrm{Radius~(AU)}$',fontsize=32)
        ax_all.tick_params(labelsize=28)
        ax_all.set_yticks(np.arange(0,R_env_max/AU,R_env_max/AU/5))
        ax_all.grid(True)

        cb = fig_all.colorbar(img_1)
        cb.solids.set_edgecolor('face')
        cb.ax.set_ylabel(r'$\mathrm{Averaged~Density~(g/cm^{3})}$',fontsize=32)
        cb_obj = plt.getp(cb.ax.axes,'yticklabels')
        plt.setp(cb_obj,fontsize=28)
        
        fig_all.savefig(outdir+'density_profile_edge.pdf',format='pdf',dpi=300,bbox_inches='tight')
        ax_all.cla()
        fig_all.clf()

        # # Plot the mid plane slice of the density profile in r-phi coordinates
        # # 
        # ax_all_midplane = fig_all.add_subplot(121,projection='polar')

        # z = np.sum(rho**2,axis=1)/np.sum(rho,axis=1)

        # img_all = ax_all_midplane.pcolormesh(phii,ri/AU,z,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=z.max()))

        # ax_all_midplane.set_xlabel(r'$\mathrm{Azimuthal~angle~(Degree)}$',fontsize=14)
        # ax_all_midplane.set_ylabel(r'$\mathrm{Radius~(AU)}$',fontsize=14)
        # ax_all_midplane.tick_params(labelsize=14)
        # ax_all_midplane.set_yticks(np.arange(0,R_env_max/AU,R_env_max/AU/5))
        # ax_all_midplane.grid(True)
        # fig_all.savefig(outdir+'density_profile.eps',format='eps',dpi=300)
        # ax_all.cla()
        # ax_all_midplane.cla()
        # fig_all.clf()

        # # In[165]:

        # # Plot the center structure of disk + envelope density profile
        # # 
        # fig_center = plt.figure(figsize=(20,7))
        # # Plot the center structure of disk + envelope density profile in r-theta coordinates
        # # 
        # ax_center  = fig_center.add_subplot(122,projection='polar')

        # z = np.sum(rho**2,axis=2)/np.sum(rho,axis=2)

        # zmin = 1e-21
        # cmap = 'jet'
        # img_center = ax_center.pcolormesh(thetai,ri/AU,z,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=z.max()))
        # ax_center.pcolormesh(PI-thetai,ri/AU,z,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=z.max()))
        # ax_center.pcolormesh(-thetai,ri/AU,z,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=z.max()))
        # ax_center.pcolormesh(PI+thetai,ri/AU,z,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=z.max()))

        # ax_center.set_xlabel(r'$\mathrm{Polar~angle~(Degree)}$',fontsize=14)
        # ax_center.set_ylabel(r'$\mathrm{Radius~(AU)}$',fontsize=14)
        # ax_center.tick_params(labelsize=14)
        # ax_center.set_ylim([0,200])
        # ax_center.set_yticks(np.arange(0,200,50))
        # ax_center.grid(True)
        # cb = fig_center.colorbar(img_center)
        # cb.ax.set_ylabel(r'$\mathrm{Surface~Density~(g/cm^{2})}$',fontsize=14)
        # cb_obj = plt.getp(cb.ax.axes,'yticklabels')
        # plt.setp(cb_obj,fontsize=14)

        # # Plot the slice at the mid plane of the center structure of disk + envelope density profile in r-theta coordinates
        # # 
        # ax_center_midplane  = fig_center.add_subplot(121,projection='polar')

        # z = np.sum(rho**2,axis=1)/np.sum(rho,axis=1)

        # img_center_midplane = ax_center_midplane.pcolormesh(phii,ri/AU,z,cmap=cmap,norm=LogNorm(vmin=zmin,vmax=z.max()))

        # ax_center_midplane.set_xlabel(r'$\mathrm{Polar~angle~(Degree)}$',fontsize=14)
        # ax_center_midplane.set_ylabel(r'$\mathrm{Radius~(AU)}$',fontsize=14)
        # ax_center_midplane.tick_params(labelsize=14)
        # ax_center_midplane.set_ylim([0,200])
        # ax_center_midplane.set_yticks(np.arange(0,200,50))
        # ax_center_midplane.grid(True)
        # fig_center.savefig(outdir+'density_profile_center.eps',format='eps',dpi=300)
        # ax_center.cla()
        # ax_center_midplane.cla()
        # fig_center.clf()

    # In[166]:

    # Write the wavelength_micron.inp file
    #

    print 'Printing the output files...'
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

    # Create the user-defined circular aperture size for given wavelength
    
    aper = np.zeros([len(lam_cam)])
    ind = 0
    for wl in lam_cam:
        if wl < 5:
            aper[ind] = 8
        elif wl >= 5 and wl < 10:
            aper[ind] = 18
        elif wl >= 10 and wl < 50:
            aper[ind] = 20
        else:
            aper[ind] = 24.5
        ind += 1

    # In[107]:

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
    f_wave_cam.write('%d \n' % int(n_lam_cam))
    for ilam in range(0,n_lam_cam):
        f_wave_cam.write('%f \n' % lam_cam[ilam])
    f_wave_cam.close()
    # In[108]:

    # Write the aperture_info.inp
    #
    f_aper = open(outdir+'aperture_info.inp','w')
    f_aper.write('1 \n')
    f_aper.write('%d \n' % int(n_lam_cam))
    for iaper in range(0, len(aper)):
        f_aper.write('%f \t %f \n' % (lam_cam[iaper],aper[iaper]))
    f_aper.close()

    # Write the stars.inp file
    #
    f_star = open(outdir+'stars.inp','w')
    f_star.write('2\n')
    f_star.write('1 \t %d \n' % int(nlam))
    f_star.write('\n')
    f_star.write('%e \t %e \t %e \t %e \t %e \n' % (rstar,mstar,pstar[0],pstar[1],pstar[2]))
    f_star.write('\n')
    for ilam in range(0,nlam):
        f_star.write('%f \n' % lam[ilam])
    f_star.write('\n')
    f_star.write('%f \n' % -tstar)
    f_star.close()


    # In[109]:

    # Write the grid file
    #
    f_grid = open(outdir+'amr_grid.inp','w')
    f_grid.write('1\n')                               # iformat
    f_grid.write('0\n')                               # AMR grid style  (0=regular grid, no AMR)
    f_grid.write('150\n')                             # Coordinate system  coordsystem<100: Cartisian; 100<=coordsystem<200: Spherical; 200<=coordsystem<300: Cylindrical
    f_grid.write('0\n')                               # gridinfo
    f_grid.write('1 \t 1 \t 1 \n')                    # Include x,y,z coordinate
    f_grid.write('%d \t %d \t %d \n' % (int(nx),int(ny),int(nz)))    # Size of the grid
    [f_grid.write('%e \n' % ri[ir]) for ir in range(0,len(ri))]
    [f_grid.write('%f \n' % thetai[itheta]) for itheta in range(0,len(thetai))]
    [f_grid.write('%f \n' % phii[iphi]) for iphi in range(0,len(phii))]
    f_grid.close()


    # In[110]:

    # Write the density file
    #
    f_dust = open(outdir+'dust_density.inp','w')
    f_dust.write('1 \n')                      # format number
    f_dust.write('%d \n' % int(nx*ny*nz))         # Nr of cells
    f_dust.write('1 \n')                      # Nr of dust species
    for iphi in range(0,len(phic)):
        for itheta in range(0,len(thetac)):
            for ir in range(0,len(rc)):
                f_dust.write('%e \n' % rho_env[ir,itheta,iphi])
    f_dust.close()


    # In[111]:

    # Write the Dust opacity control file
    # 
    f_opac = open(outdir+'dustopac.inp','w')
    f_opac.write('2               Format number of this file\n')
    f_opac.write('1               Nr of dust species\n')
    f_opac.write('============================================================================\n')
    f_opac.write('1               Way in which this dust species is read\n')
    f_opac.write('0               0=Thermal grain\n')
    # f_opac.write('klaus           Extension of name of dustkappa_***.inp file\n')
    f_opac.write('oh5_extended    Extension of name of dustkappe_***.inp file\n')
    f_opac.write('----------------------------------------------------------------------------\n')
    f_opac.close()
    
    # Copy the dust opacity table to the folder
    #
    import shutil
    shutil.copy2('/home/bettyjo/yaolun/radmc_simulation/bhr71/dustkappa_klaus.inp', outdir)
    shutil.copy2('/home/bettyjo/yaolun/radmc_simulation/bhr71/dustkappa_oh5_extended.inp', outdir)
    

    # In[112]:

    # Write the radmc3d.inp control file
    #
    f_control = open(outdir+'radmc3d.inp','w')
    f_control.write('nphot = %d \n' % nphot)
    f_control.write('scattering_mode_max = 2\n')
    f_control.write('camera_min_drr = 0.1\n')
    f_control.write('camera_min_dangle = 0.1\n')
    f_control.write('camera_spher_cavity_relres = 0.1\n')
    f_control.write('istar_sphere = 1\n')
    f_control.write('modified_random_walk = 1\n')
    f_control.close()


