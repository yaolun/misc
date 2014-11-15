pro sofia_plot
; Constants setup
;
lab = 63.1836709
c   = 2.99792458e10
k   = 1.380658e-16
pix_size = 9.4^2

; um & Jy
readcol, '~/data/digit_v65_slim/L1551-IRS5/L1551-IRS5_centralSpaxel_PointSourceCorrected_CorrectedYES_trim.txt',format='D,D',wl1,flux1
readcol, '~/data/digit_v65_slim/L1551-IRS5/data/L1551-IRS5_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_smooth_sed.txt',format='D,D',wl,cont_dum
cont = flux1*0
for i =0, n_elements(wl1)-1 do begin
	for j = 0, n_elements(wl)-1 do begin
		if wl1[i] eq wl[j] then begin
			cont[i] = cont_dum[j]
			break
		endif
	endfor
endfor
flux1 = (flux1-cont)*1e-23
; Jy to erg/s/cm2/Hz

; km/s offset & K
readcol, '~/l1551-irs5_co21.txt',format='D,D,D,D',velo2,t2,dum,dumm
; For [OI] observation proposal
velo1 = c*(wl1-lab)/lab/1e5
omega_b = !pi/4/alog(2)*pix_size*(1/3600.0/180*!pi)^2		; sr
t1 = (wl1*1e-4)^2/2/k*flux1/omega_b
ind = where(velo1 lt 400 and velo1 gt -400)
velo1 = velo1[ind]
t1 = t1[ind]
print, min(t2)
t1 = t1-min(t2)+0.0000001
t2 = t2-min(t2)+0.0000001
print, min(t2)

set_plot,'ps'
!p.font=0
loadct,13,/silent
device, filename='~/Dropbox/proposal/sofia/oi_co2-1.eps',/helvetica,/portrait,/encapsulated,isolatin=1,font_size=12,decomposed=0,/color
!p.thick=4 & !x.thick=4 & !y.thick=4	
plot, velo2,alog10(t2), xrange=[-300,300], psym=10, xtitle = 'Velocity (km/s)', ytitle = 'log T (K)',yrange=[-1,1]
oplot, velo1, alog10(t1),color=70;,psym=10
al_legend,['[OI] 63 '+greek('mu',/append_font)+'m ','CO J=2-1'],textcolor=[70,0],/right
al_legend, ['L1551-IRS5'],/left,box=0
device, /close_file, decomposed=1
!p.multi=0

end