def gauss2d(x, y, sigmax, sigmay=None):
    # sigmax needs to be in pixel coordinates
    if sigmay == None:
        sigmay = sigmax
    return 1/(2*np.pi*sigmax*sigmay) * np.exp( -(x**2/2./sigmax**2 + y**2/2./sigmay**2) )

def gridding(gauss_beam, size=1201., phys_size=120.):
    # gauss_beam is the FWHM of the Gaussian
    # convert from aperture (diameter) to radius

    grid_x, grid_y = np.meshgrid(np.linspace(0,size-1,size), np.linspace(0,size-1,size))
    grid_x = grid_x - (size-1)/2.
    grid_y = grid_y - (size-1)/2.
    grid_gauss2d = gauss2d(grid_x,grid_y, sigmax=(size-1)/phys_size*(gauss_beam/2.354))
    dA = ((1/((size-1)/2.))*phys_size)**2

    return (grid_x, grid_y, grid_gauss2d)

# Numerically Gaussian integral within a certain aperture at a given location
# map a 120" x 120" grid with 1200 x 1200 pixels
def Gaussian_anywhere(grid, ra_offset, dec_offset, pixel_aper, size=1201., phys_size=120.):
    grid_x, grid_y, grid_gauss2d = grid

    # convert from physcial coordinates to pixel coordinates
    x = (ra_offset-phys_size/2.) * (size-1)/2./(phys_size/2.) + (size-1)/2.
    y = (dec_offset-phys_size/2.) * (size-1)/2./(phys_size/2.) + (size-1)/2.
    r_pix = (pixel_aper/2.) * (size-1)/phys_size
    grid_dist = ((grid_x-x)**2+(grid_y-y)**2)**0.5
    gauss2d_mask = np.where(grid_dist<=r_pix, grid_gauss2d,0)

    return np.sum(gauss2d_mask)

def spire_extration(spec, gauss_beam, modules='both'):
    
