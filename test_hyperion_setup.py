import numpy as np
import astropy.constants as const

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

indir = '/Users/yaolun/bhr71/radmc3d_params'
outdir = '/Users/yaolun/programs/misc/TSC/'

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
# mostly fixed parameter
M_disk    = 0.5*MS
beta      = 1.093
h100      = 8.123*AU
rho_cav   = 1e-21

nx        = 100L
ny        = 400L
nz        = 50L
grid = np.array([nx, ny, nz])

# Do the variable conversion
cs = (G * M_env_dot / 0.975)**(1/3.)  # cm/s
t = R_inf / cs / yr   # in year
mstar = M_env_dot * t * yr
omega = (R_cen * 16*cs**8 / (G**3 * mstar**3))**0.5

print d_sub/AU

import pidly
idl = pidly.IDL('/Applications/exelis/idl82/bin/idl')
idl('.r ~/programs/misc/TSC/tsc.pro')
idl.pro('tsc_run', outdir=outdir, grid=[nx,ny,nz], time=t, c_s=cs, omega=omega, rstar=rstar, renv_min=R_env_min, renv_max=R_env_max)
