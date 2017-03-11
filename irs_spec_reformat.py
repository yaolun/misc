filename = 'wl12-j'

data = open('/Users/yaolun/Dropbox/cops-spire/IRS_spec/'+filename+'.txt','r')

obj = 'wl12'
data_start = 2
# columns for wavelength, flux, error
selected_cols = [0,1,2]
error_NaN = False

# new file to write
foo = open('/Users/yaolun/Dropbox/cops-spire/IRS_spec/reformatted/'+obj+'.txt','w')
foo.write('{:>18s}{:>18s}{:>18s}\n'.format('Wavelength(um)','Flux_Density(Jy)','Uncertainty(Jy)'))

if not error_NaN:
    for i in data.readlines()[data_start:]:
        foo.write('{:>18s}{:>18s}{:>18s}\n'.format(i.split()[selected_cols[0]],
                                                   i.split()[selected_cols[1]],
                                                   i.split()[selected_cols[2]],))
else:
    for i in data.readlines()[data_start:]:
        foo.write('{:>18s}{:>18s}{:>18s}\n'.format(i.split()[selected_cols[0]],
                                                   i.split()[selected_cols[1]],
                                                   'NaN',))

foo.close()
data.close()
