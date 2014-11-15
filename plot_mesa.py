import numpy as np
import matplotlib.pyplot as plt
import os
home = os.path.expanduser('~')

dirpath = '/agb_star_and_dust/data/agb_mw/'
#dirpath = '/hw3/'
plotpath = '/agb_star_and_dust/plots/'
#plotpath = dirpath

filepath = str(raw_input('Filename:'))
filename = home+dirpath+filepath+'.data'
data = np.genfromtxt(filename, skip_header=5, dtype='str')
name = data[0,:]
x_name = str(raw_input('Name of X axis:'))
y_name = str(raw_input('Name of Y axis:'))
x = np.where(name == x_name)
y = np.where(name == y_name)
x_name = x_name.replace('_',' ')
y_name = y_name.replace('_',' ')
plt.plot(data[1:,x].astype('float').flat,data[1:,y].astype('float').flat,'-')
plt.xlabel(x_name)
plt.ylabel(y_name)
if x_name == 'log Teff':
    plt.xlim([max(data[1:,x].astype('float').flat),min(data[1:,x].astype('float').flat)])
    print 'Change the axis'
#if x_name == 'star age':
#	plt.xlim([1.2e10,1.25e10])
plt.savefig(home+plotpath+filepath+'_'+x_name+'_'+y_name+'_mw.eps',format='eps',dpi=300)
plt.cla()
plt.clf()
print 'Done!'
