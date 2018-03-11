def make_image(filepath, outdir=None, coord=None, coord_unit='hourangle',
       size=0.03, plotname=None, stretch='log', vmin=None, vmax=None,
       bar_size=30, aper=None, int_unit='MJy/sr', obj_text=None,
       scalebar=True, framecolor='white', no_marker=False, marker_list=[]):
    """
    size in degree
    """
    import numpy as np
    import aplpy as apl
    import matplotlib
    import matplotlib.pyplot as plt
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    import os

    def coord2SkyCoord(coord, coord_unit):
        if coord_unit == 'hourangle':
            c = SkyCoord(coord, unit=(u.hourangle, u.deg))
        elif coord_unit == 'deg':
            c = SkyCoord(coord, unit=(u.deg, u.deg))

        return c

    im = apl.FITSFigure(filepath, north=True)
    im.frame.set_linewidth(2)
    im.frame.set_color(framecolor)
    cmap = plt.cm.viridis
    im.show_colorscale(cmap=cmap, stretch=stretch, vmin=vmin, vmax=vmax)

    im.set_nan_color('white')
    im.tick_labels.set_font(size=18)

    im.add_colorbar()
    im.colorbar.set_font(size=18)
    im.colorbar.set_axis_label_text('Intensity ('+int_unit+')')
    im.colorbar.set_axis_label_font(size=18)
    im.colorbar.set_frame_linewidth(2)
    im.tick_labels.set_xformat('hh:mm:ss')
    im.tick_labels.set_yformat('dd:mm:ss')
    im.ticks.set_xspacing('auto')
    im.ticks.set_yspacing('auto')
    im.axis_labels.set_font(size=18)
    im.ticks.set_linewidth(2)
    im.ticks.set_color(framecolor)

    if coord != None:
        center = coord2SkyCoord(coord, coord_unit)

    # re-center the image based on the input central coordinates
        im.recenter(center.ra.deg, center.dec.deg, radius=size)

    if not no_marker:
        im.show_markers(center.ra.degree,
                        center.dec.degree, marker='+', c='red', s=120, linewidth=2)

    if len(marker_list) > 0:
        for marker in marker_list:
            c = coord2SkyCoord(marker, coord_unit)
            im.show_markers(c.ra.degree, c.dec.degree,
                            marker='+', c='white', s=120, linewidth=2)

    if scalebar:
        im.add_scalebar(bar_size/3600.)
        im.scalebar.set_length(bar_size * u.arcsecond)
        im.scalebar.set_label(r'$\rm{'+str(bar_size)+'\,arcsec}$')
        im.scalebar.set_font(size=22, weight='bold')
        im.scalebar.set(linestyle='solid', color='w', linewidth=3, corner='bottom')

    # print image information
    if obj_text != None:
        im.add_label(0.85, 0.92, obj_text, size='xx-large', weight='bold', relative=True, color='w')

    # plot a circular region for aperture
    if aper != None:
        im.show_circles([center.ra.deg], [center.dec.deg], [aper/3600.], color='lime', linewidth=3)

    if outdir == None:
        outdir = os.path.dirname(filepath)

    if plotname == None:
        plotname = 'image'
    im.save(outdir+'/'+plotname+'.pdf', format='pdf', transparent=True, dpi=300)

    return im
