def lin_leastsqfit(xdata,ydata,sig_y,nofit=False):#,xlabel,ylabel,outdir,plotname
	import numpy as np
	import matplotlib
	import matplotlib.pyplot as plt
	import os
	from numpy.linalg import inv
	#home = os.path.expanduser('~')
	#Read the data	
	#data = np.loadtxt(home+'/bhr71/data/'+filename+'.txt')
	#data = [x,y,sig_y]
	#Construct the matrix
	n = len(xdata)#len(data[:,0])
	m = 1
	x = np.zeros((n,2))
	x = np.matrix(x)
	for i in range(0,n):#(len(data[:,0])):
		x[i,1] = xdata[i]#data[i,0]
		x[i,0] = 1
	w = np.zeros((n,n))
	w = np.matrix(w)
	for i in range(0,n):
		w[i,i] = 1/sig_y[i]**2
	y = np.zeros((n,1))
	for i in range(0,n):
		y[i,0] = ydata[i]#data[i,1]
	y = np.matrix(y)
	#Least square fitting part
	N = x.transpose()*w*x
	a_hat = inv(N)*x.transpose()*w*y
	cov = inv(N)
	s_min = (y-x*a_hat).transpose()*w*(y-x*a_hat)
	cov_hat = float(s_min)/(n-m-1)*cov
	t_rot = -1/a_hat[1]*np.log10(np.e)
	sig_t_rot = -t_rot*cov_hat[1,1]**0.5/a_hat[1]*np.log10(np.e)
	#print str(t_rot)+'+/-'+str(sig_t_rot)+'K'
	#ploting
	if nofit == True:
		#only data
		plt.plot(xdata,ydata+36,'go',linestyle='None')
		plt.errorbar(xdata,ydata+36,yerr=sig_y,linestyle='None',color='g')
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		plt.savefig(outdir+plotname+'_nofit.eps',format='eps',dpi=300)
		plt.cla()
		plt.clf()
		return np.squeeze(xdata),np.squeeze(ydata),np.squeeze(sig_y)
	else:
	#data w/ 1d fit w/ yerr
		yerr = (cov_hat[0,0]+2*cov_hat[0,1]*xdata+cov_hat[1,1]*xdata**2)**0.5
		yerr.shape = (n,1)
		yfit = a_hat[0]+a_hat[1]*xdata
		yfit.shape = (n,1)
		xi = xdata
		xi.shape = (n,1)
		#data, = plt.plot(xdata,ydata+36,'o',color='DarkGreen',linestyle='None')
		#fit, = plt.plot(xdata,yfit+36)
		#plt.plot(xi,yfit+yerr+36,'m--')
		#plt.plot(xi,yfit-yerr+36,'m--')
		#plt.fill_between(xi,yfit+yerr,yfit-yerr,facecolor='green',alpha=0.5)
		#plt.errorbar(xdata,ydata+36,yerr=sig_y,linestyle='None',color='DarkGreen')
        #plt.xlabel(xlabel,fontsize=14)
        #plt.ylabel(ylabel,fontsize=14)
        #plt.xlim([500,3500])
        #plt.ylim([42,48])
        #plt.legend([data,fit],['Data','Fit'],numpoints=1,loc='upper right')
        #plt.title('%8.4f +/- %8.6f K' % (t_rot, sig_t_rot),fontsize=14)
        #plt.savefig(outdir+plotname+'.eps',format='eps',dpi=300)
        #plt.cla()
        #plt.clf()
        return np.squeeze(np.array(yfit)),np.squeeze(np.array(yerr)),t_rot,sig_t_rot,float(s_min)/(n-m-1),a_hat[0]
# pix_name = [3,4,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
# for pix in pix_name:
# 	import numpy as np
# 	import os

# 	home = os.path.expanduser('~')
# 	filename = home+'/bhr71/pacs_co_rot_pixel'+str(pix)+'.txt'
# 	data = np.loadtxt(filename)
# 	wl = data[:,0]
# 	flux = data[:,1]*1e7
# 	sig_flux = data[:,2]*1e7
# 	A = data[:,3]
# 	Eu = data[:,4]
# 	g = data[:,5]
# 	c = 3e10
# 	h = 6.626e-27
# 	p=3.086
# 	v = c/(wl*1e-4)
# 	N = 4*np.pi*flux*(178*p)**2/(A*h*v)
# 	N_sigma = 4*np.pi*sig_flux*(178*p)**2/(A*h*v)
# 	x = Eu
# 	y = np.log10(N/g)
# 	yerr_hi = np.log10((N+N_sigma)/g)-np.log10(N/g)
# 	yerr_low = np.log10(N/g)-np.log10((N-N_sigma)/g)
# 	y_sig = y*0
# 	for i in range(0,len(y)):
# 		y_sig[i] = max(yerr_hi[i], yerr_low[i])
# 	lin_leastsqfit(x, y, y_sig,'$E_{u} (K)$','$log(\mathcal{N}/g)$',home+'/bhr71/plots/rotational_diagram/','bhr71_pacs_pixel'+str(pix)+'_co')
