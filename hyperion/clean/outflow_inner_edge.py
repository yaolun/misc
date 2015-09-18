def outflow_inner_edge(rho, grid, M_env_dot0, v, theta_cav, R_min, R_max=None):
    """
    grid = (ri, thetai, phii), walls
    M_env_dot0 in g/s
    v in cm/s
    theta_cav in degree
    R_min in cm
    """
    import numpy as np
    (ri, thetai, phii) = grid
    rc           = 0.5*( ri[0:len(ri)-1]     + ri[1:len(ri)] )
    thetac       = 0.5*( thetai[0:len(thetai)-1] + thetai[1:len(thetai)] )
    phic         = 0.5*( phii[0:len(phii)-1]   + phii[1:len(phii)] )

    # calculate the enclosed mass in the outflow cavity at different radius
    c0 = max(ri)**(-0.5)*np.sqrt(1/np.sin(np.radians(theta_cav))**3-1/np.sin(np.radians(theta_cav)))

    cavity_mass = np.empty(len(ri)-1)
    # cavity_mass_tot = 0
    for ir in range(1, len(ri)):
        if ri[ir] < R_min:
            cavity_mass[ir-1] = 0
            continue
        if R_max != None:
            if ri[ir] > R_max:
                cavity_mass[ir-1] = 0
                continue
        mass_dum = 0
        for itheta in range(1, len(thetai)):
            for iphi in range(1, len(phii)):
                w = abs(rc[ir-1]*np.cos(np.pi/2 - thetac[itheta-1]))
                z = rc[ir-1]*np.sin(np.pi/2 - thetac[itheta-1])
                z_cav = c0*abs(w)**1.5
                # truncate non-outflow region
                if abs(z) < abs(z_cav):
                    rho[ir-1, itheta-1, iphi-1] = 0.0
                V = (1/3.)*(ri[ir]**3 - ri[ir-1]**3) * (phii[iphi]-phii[iphi-1]) * -(np.cos(thetai[itheta])-np.cos(thetai[itheta-1]))
                if V < 0:
                    print ri[ir], thetai[itheta], phii[iphi]
                    print (1/3.)*(ri[ir]**3 - ri[ir-1]**3) * (phii[iphi]-phii[iphi-1]) * -(np.cos(thetai[itheta])-np.cos(thetai[itheta-1]))
                cell_mass = rho[ir-1, itheta-1, iphi-1] * (1/3.)*(ri[ir]**3 - ri[ir-1]**3) * (phii[iphi]-phii[iphi-1]) * -(np.cos(thetai[itheta])-np.cos(thetai[itheta-1]))
                mass_dum = mass_dum + cell_mass
        cavity_mass[ir-1] = mass_dum
        # cavity_mass_tot = cavity_mass_tot + mass_dum

    # flow mass derived from mass accretion rate
    flow_mass = M_env_dot0 * ri[1:]/v

    cavity_mass = cavity_mass[cavity_mass != 0]
    flow_mass = flow_mass[cavity_mass != 0]

    diff = cavity_mass - flow_mass
    ind = np.argsort(abs(diff))[1]

    return ri[ri > R_min][ind]#, cavity_mass, cavity_mass_tot
