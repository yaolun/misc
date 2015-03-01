def read_params_table(filename,outdir,num=None):
    import numpy as np
    import os
    import shutil as sh
    data = np.genfromtxt(filename,skip_header=1,dtype=None)
    if num == None:
        num = len(data[:])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for i in range(0,num):
        real_outdir = os.path.dirname(outdir+str(data[i][0])+'/params.dat')
        model = str(data[i][0])
        if not os.path.exists(real_outdir):
            os.makedirs(real_outdir)
        print data[i][0]
        foo = open(outdir+str(data[i][0])+'/params.dat','w')
        foo.write('# Parameters Setup \n# \n')
        foo.write('tstar \t \t %e \t # Stellar temperature         [K]\n' % data[i][4])
        foo.write('mstar \t \t %e \t # Stellar mass                [Solar Mass]\n' % data[i][2])
        foo.write('rstar \t \t %e \t # Stellar radius              [Solar radius]\n' % data[i][3])
        foo.write('M_env_dot \t %e \t # Envelope accretion rate     [Solar mass/yr]\n' % data[i][5])
        foo.write('M_disk_dot \t %e \t # Disk accretion rate         [Solar mass/yr]\n' % data[i][15])
        foo.write('R_env_max \t %e \t # Envelope outer radius       [AU]\n' % data[i][6])
        foo.write('R_env_min \t %e \t # Envelope inner radius       [AU]\n' % data[i][7])
        foo.write('theta_cav \t %f \t \t # Outflow cavity opening angle[deg]\n' % data[i][8])
        foo.write('R_disk_max \t %e \t # Disk outer radius           [AU]\n' % data[i][10])
        foo.write('R_disk_min \t %e \t # Disk inner radius           [AU]\n' % data[i][11])
        foo.write('M_disk \t \t %e \t # Disk mass                   [Solar mass]\n' % data[i][9])
        foo.write('beta \t \t %f \t \t # Disk flare factor           []\n' % data[i][13])
        foo.write('h100 \t \t %f \t \t # Disk scale height at 100 AU [AU]\n' % data[i][16])
        foo.write('rho_cav \t %e \t # Outflow cavity density      [g/cm3]\n' % data[i][14])
        foo.write('wall \t \t %e \t \t # Thickness of denser cavity wall [AU]\n' % data[i][17])
        foo.write('rho_wall \t %e \t # Wall density          [g/cm3]\n' % data[i][18])
        if len(data[i]) >= 18:
            foo.write('rho_cav_center \t %e \t # Density of Cavity center [g/cm3]\n' % data[i][19])
            foo.write('rho_cav_edge \t %e \t # Size of the inner cavity region [AU]\n' % data[i][20])
        foo.write('\n')
        foo.close()
        # Save a copy in the 'radmc3d_params' folder
        if not os.path.exists(os.path.dirname(copy_dir+model+'.dat')):
            os.makedirs(os.path.dirname(copy_dir+model+'.dat'))
        sh.copyfile(outdir+model+'/params.dat',copy_dir+model+'.dat')
        
        # Flags check
        # denser cavity wall
        if data[i][17] == 0:
            denser_wall = False
        else:
            denser_wall = True
        
        # Initial inputs calculation
        problem_setup(outdir+model,outdir+model,denser_wall=denser_wall,plot=model_plot,low_res=low_res)
        
        
        