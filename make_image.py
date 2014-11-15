def make_image(filepath,scale=None,output=False,mag=1.5):
	import numpy as np
	import matplotlib.pyplot as plt
	from astropy.io import fits
	import os

	data = fits.open(filepath).data

	fig = plt.figure(figsize=(6*mag,8*mag))
	ax_im = fig.add_subplot(111)

	ax_im.imshow(data)