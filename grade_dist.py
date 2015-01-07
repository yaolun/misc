def grade_dist(filepath,title,outdir):
	import numpy as np
	import matplotlib.pyplot as plt


	data = np.genfromtxt(filepath,dtype=None,delimiter=',',skip_header=1,skip_footer=1).T
	score = np.empty((len(data)))
	for i in range(0,len(data)):
		score[i] = float(data[i][5])

	# Histogram plot
	mag = 1.5
	fig = plt.figure(figsize=(mag*8,mag*6))
	ax = fig.add_subplot(111)

	n, bins, patches = ax.hist(score,bins=10,range=(5,105),histtype='bar',align='mid',linewidth=1.5*mag)
	plt.setp(patches, 'facecolor', 'MediumAquaMarine', 'alpha', 0.75)
	ax.set_title(title,fontsize=16*mag)
	ax.set_xlabel('Score',fontsize=16*mag)
	ax.set_ylabel('No. of students',fontsize=16*mag)
	ax.set_xlim([0,105])
	ax.set_ylim([0,30])
	[ax.spines[axis].set_linewidth(1.5*mag) for axis in ['top','bottom','left','right']]
	ax.tick_params('both',labelsize=mag*16,width=1.5*mag,which='major',length=5*mag)
	fig.savefig(filename=outdir+'score_distribution.eps',format='eps',dpi=300,bbox_inches='tight')

#import matplotlib
#print matplotlib

filepath = '/Users/yaolun/Dropbox/UT-Austin/TA/AST301_Robinson/Final_exam/12151059C.al.csv'
title = r'$\mathrm{Final~Exam}$'
outdir = '/Users/yaolun/Dropbox/UT-Austin/TA/AST301_Robinson/Final_exam/'

grade_dist(filepath,title,outdir)
