pro aper_flux, indir=indir,filename=filename,aperture=aperture,object=object

; Constants setup
; 
c = 2.99792458e10

; Read the aperture profile
; Unit: um & arcsec (radius)
readcol,aperture,format='D,D', wl_a, r, /silent

if not keyword_set(filename) then filename='image.fits'
; Image has the unit of Jy/pixel
imag = readfits(indir+filename, hdr)
freq = readfits(indir+filename,exten=1)
wl = c/freq*1e4
pixsize_x = sxpar(hdr,'CDELT1')*3600
pixsize_y = sxpar(hdr,'CDELT2')*3600
pixsize = abs(pixsize_x)
xc = sxpar(hdr,'CRPIX1')
yc = sxpar(hdr,'CRPIX2')
ra = sxpar(hdr,'CRVAL1')
dec = sxpar(hdr,'CRVAL2')
; IDL takes array index starting from 0
xc = xc-1
yc = yc-1
sky_inner = 20./pixsize
sky_outer = 30./pixsize
num_slice = n_elements(wl)
total_flux = []
total_flux_err = []
print, format='(a7,2x,a8,2x,a8,2x,a5)','Wave','Flux','Flux_sig','rad.'
for i = 0, num_slice-1 do begin
	r_dum = interpol(r, wl_a, wl[i])
	r_pix = r_dum/pixsize
	aper,imag[*,*,i],xc,yc,/exact,/flux,flux,errap,sky,skyerr,1,[r_pix],-1,setskyval=0,/nan,/silent;[sky_inner,sky_outer]
	; The output 'flux' from aper is the direct summation within the givern aperture
	; pixel size is diameter or radius?
	flux = flux*(!PI*r_pix^2);*(!pi/4/alog(2)*r_dum^2*(!PI/180/3600)^2)           				    ; unit: Jy
	errap = errap*(!PI*r_pix^2);*(!pi/4/alog(2)*r_dum^2*(!PI/180/3600)^2)           				; unit: Jy
	total_flux = [total_flux,flux]
	total_flux_err = [total_flux_err,errap]
	print, format='(g7.4,2x,g8.3,2x,g8.3,2x,g5.4)',wl[i],flux,errap,r_dum
endfor

; Plotting
;
; Create filled circle symbol
; Make a vector of 16 points, A[i] = 2pi/16:
A = FINDGEN(17) * (!PI*2/16.)
; Define the symbol to be a unit circle with 16 points, 
; and set the filled flag:
size = 0.5
USERSYM, size*COS(A), size*SIN(A), /FILL

;
set_plot,'ps'
!p.font=0
loadct,13,/silent
device, filename=indir+'photo_sed.eps',/helvetica,/portrait,/encapsulated,isolatin=1,font_size=16,decomposed=0,/color
!p.thick=4 & !x.thick=4 & !y.thick=4	
wl = alog10(wl)
vsv = alog10(total_flux*1e-23*c/(wl*1e-4))
plot, wl, vsv, psym=8, $
	xtitle = 'log '+greek('lambda',/append_font)+' ['+greek('mu',/append_font)+'m]', ytitle = 'log '+greek('nu',/append_font)+'S!d'+greek('nu',/append_font)+'!n [erg/s/cm!u2!n]'
al_legend, [object],/left,box=0
device, /close_file, decomposed=1
!p.multi=0

end
