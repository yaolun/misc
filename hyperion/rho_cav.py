def rho_cav(t):
    # t in years
    import numpy as np
    import astropy.constants as const

    # constants setup
    RS = const.R_sun.cgs.value
    LS = const.L_sun.cgs.value
    MS = const.M_sun.cgs.value
    G  = const.G.cgs.value
    AU = const.au.cgs.value
    yr = 3600*24*365.

    # variable setup
    mstar = 0.105 * MS
    rstar = 3 * RS
    L_bol = 13.92 * LS
    # theta_cav = np.radians(30)
    theta_cav = np.radians(20)
    v = 3e7
    t = 3600.*24*365 * t
    g2d = 100

    M_dot_o = rstar * L_bol / (10.*G*mstar) / g2d
    # derive from B97 CO momentum
    p = 11 * MS*1e5 # M_sun km/s to g cm/s
    t_outflow = 1e4 *yr
    M_dot_o = p/t_outflow/v / g2d # g/s in dust mass

    c0 = (v*t)**-0.5 * np.cos(theta_cav)/np.sin(theta_cav)**1.5

    V = (v*t)**3 * (4/3.*np.pi - 2*np.pi*np.cos(theta_cav) + 2*np.pi*np.cos(theta_cav)**2 - 2*np.pi/3.*np.cos(theta_cav)**3) + 6/7.*np.pi/c0**2*np.cos(theta_cav)**(7/3.)*(v*t)**(-2/3.)

    return M_dot_o*t/V, v*t/AU, M_dot_o/MS*yr, V, M_dot_o /MS*yr

print rho_cav(1.)
