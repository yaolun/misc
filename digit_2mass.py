from astropy import units as u
from subprocess import Popen, call
from astropy.coordinates import SkyCoord
from astropy.io import fits, ascii
import os
from astropy.io.votable import parse
import urllib

def make_image(filepath, coord, size=0.03, plotname=None, plotdir=None, stretch='log', vmin=None,
               bar_size=30, aper=None, int_unit='MJy/sr', text=None):
    import numpy as np
    import aplpy as apl
    import matplotlib
    import matplotlib.pyplot as plt
    from astropy.io import fits
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    import os

    mag = 1.5

    print coord

    im = apl.FITSFigure(filepath,north=True)
    cmap = plt.cm.viridis
    im.show_colorscale(cmap=cmap,stretch=stretch, vmin=vmin)

    im.recenter(coord.ra.degree,coord.dec.degree,radius=size)
    im.show_markers(coord.ra.degree, coord.dec.degree, marker='+', c='r', s=120, linewidth=2)

    im.set_nan_color('black')
    im.tick_labels.set_font(size=18)

    im.add_colorbar()
    im.colorbar.set_font(size=20)
    im.colorbar.set_axis_label_text('Intensity ('+int_unit+')')
    im.colorbar.set_axis_label_font(size=20)
    im.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss')
    im.ticks.set_xspacing(2/60.)
    im.ticks.set_yspacing(1/60.)
    im.axis_labels.set_font(size=20)
    im.ticks.set_linewidth(2)

    im.add_scalebar(bar_size/3600.)
    im.scalebar.set_length(bar_size * u.arcsecond)
    im.scalebar.set_label(r'$\rm{'+str(bar_size)+'\,arcsec}$')
    im.scalebar.set_font(size=20, weight='bold')
    im.scalebar.set(linestyle='solid', color='white', linewidth=3)

    # print image information
    if text != None:
        im.add_label(0.8, 0.9, text, size='xx-large', weight='bold', relative=True, color='white')

    # plot a circular region for aperture
    if aper != None:
        im.show_circles([coord.ra.degree], [coord.dec.degree], [aper/3600.], color='lime', linewidth=3)

    if plotdir == None:
        plotdir = os.path.dirname(filepath)
    if plotname == None:
        plotname = raw_input('What is the plot name: ')
        print plotname

    im.save(plotdir+'/'+plotname+'.pdf',format='pdf',transparent=True, dpi=300)

# read in DIGIT source list
filename = '/Users/yaolun/data/digit_source'
digit = ascii.read(filename)
digit_coord = []
for i in range(len(digit)):
    c = SkyCoord(digit['RA'][i]+' '+digit['DEC'][i], unit=(u.hourangle, u.deg))
    digit_coord.append(c)

objdir = '/Users/yaolun/data/digit_hst/2MASS/'

for coord in digit_coord:
    # search on 2MASS image server
    obj = digit['Name'][digit_coord.index(coord)]
    print obj

    if not os.path.exists(objdir+obj):
        os.mkdir(objdir+obj)

    size = '0.01'
    band = 'H'
    run = Popen(['curl','-o', objdir+obj+'/'+obj+'.xml',
                 'http://irsa.ipac.caltech.edu/cgi-bin/2MASS/IM/nph-im_sia?'+\
                 'POS='+str(coord.ra.degree)+','+str(coord.dec.degree)+'&'+\
                 'SIZE='+size+'&band='+band+'&FORMAT=image/fits'])
    run.communicate()

    # parse xml
    votable = parse(objdir+obj+'/'+obj+'.xml')
    table = votable.get_first_table()
    # this url list may contain more than one url.
    url = table.array['download'].data

    if len(url) == 0:
        print 'Image not found, skipping '+obj
        continue

    print 'Downloading '+str(len(url))+' file(s)'
    filepath = []
    for link in url:
        if not os.path.exists(objdir+obj+'/data'):
            os.mkdir(objdir+obj+'/data')
        os.chdir(objdir+obj+'/data')
        urllib.urlretrieve(link, link.split('name=')[1])
        filepath.append(objdir+obj+'/data/'+link.split('name=')[1])

    for foo in filepath:
    #     int_unit = fits.open(foo)[1].header
        make_image(foo, coord, size=0.03, plotname=os.path.basename(foo).split('.')[0],
                   plotdir=objdir+obj+'/data/', stretch='arcsinh', vmin=None,
                   bar_size=30, aper=None, int_unit='A.U.', text=obj+' 2MASS H band')