# Test script for calculating Ulrich envelope using hyperion build-in function
import numpy as np
import astropy.constants as const
from hyperion.model import Model
from hyperion.model import AnalyticalYSOModel

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


# Variable setup
#
rstar     = 5 * RS
tstar     = 4500
R_env_max = 1.000000e+04 * AU
R_env_min = 8.000000e-01 * AU
R_cen     = 1.500000e+01 * AU
theta_cav = 0.0            # no cavity
M_env_dot = 3.000000e-05 * MS/yr
mstar     = 9.844506e-01 * MS
rin       = rstar
rout      = R_env_max
rcen      = R_cen

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
envelope.rc   = rcen         # Centrifugal radius
envelope.rmax = R_env_max    # Outer radius
envelope.star = source

# grid
# Make the Coordinates
# Grid Parameters
nx        = 100L
ny        = 400L
nz        = 50L
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

m.set_spherical_polar_grid(ri, thetai, phii)

# dust
# Hyperion needs nu, albedo, chi, g, p_lin_max
# from hyperion.dust import HenyeyGreensteinDust
# # Read in the dust opacity table used by RADMC-3D
# dust_radmc = dict()
# [dust_radmc['wl'], dust_radmc['abs'], dust_radmc['scat'], dust_radmc['g']] = np.genfromtxt('dustkappa_oh5_extended.inp',skip_header=2).T
# # opacity per mass of dust?
# dust_hy = dict()
# dust_hy['nu'] = c/dust_radmc['wl']*1e4
# ind = np.argsort(dust_hy['nu'])
# dust_hy['nu'] = dust_hy['nu'][ind]
# dust_hy['albedo'] = (dust_radmc['scat']/(dust_radmc['abs']+dust_radmc['scat']))[ind]
# dust_hy['chi'] = (dust_radmc['abs']+dust_radmc['scat'])[ind]
# dust_hy['g'] = dust_radmc['g'][ind]
# dust_hy['p_lin_max'] = 0*dust_radmc['wl'][ind]     # assume no polarization

# d = HenyeyGreensteinDust(dust_hy['nu'], dust_hy['albedo'], dust_hy['chi'], dust_hy['g'], dust_hy['p_lin_max'])

envelope.dust = 'oh5.hdf5'

# Setting up images and SEDs
image = m.add_peeled_images()
image.set_wavelength_range(100, 2.0, 670.0)
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
m.set_n_photons(initial=100000, imaging=100000, raytracing_sources=100000, raytracing_dust=100000)
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

m.write('ulrich_env.rtin')

