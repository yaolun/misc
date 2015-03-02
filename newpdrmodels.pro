pro newpdrmodels, object, molecule, indir=indir, outdir=outdir, plotdir=plotdir, moddir=moddir

;  if not keyword_set(pixels) then pixels = slw
  if not keyword_set(indir) then indir = '~/' + object + '/lines/'
  if not keyword_set(outdir) then outdir = '~/' + object + '/lines/'
  if not keyword_set(plotdir) then plotdir = '~/' + object + '/plots/'
  if not keyword_set(moddir) then moddir = '~/models/'

  line_name_co = ['CO15-14','CO14-13','CO13-12','CO12-11','CO11-10','CO10-9','CO9-8','CO8-7','CO7-6','CO6-5','CO5-4','CO4-3','CO3-2','CO2-1','CO1-0']
  line_center_co = [173.6314272,185.9992957,200.27751475,216.93275100,236.61923625,260.24634206,289.12760810,325.23334516,371.65973939,433.56713410,520.24411585,650.26787364,866.96,1300.40,2600.76]
  ylabel = 'log Intensity (erg s!u-1!n cm!u-2!n arcsec!u-2!n)'  ;brightness switch
  xlabel = 'log Wavelength (' + Greek('mu',/append_font) + 'm)'
  legend = ['Molecular Cloud', 'Shock Model', 'PDR Model']
  colors = [31,96,192];center 95 - lt blue, outflows 112 - purple, nonoutflows 31 - green, shock model 96 - blue, pdr model 192 - red

  k = 1.3806488d-16 ; boltzmann constant in erg/K
  
;read in new pdr model in K km/s units
readcol, moddir + 'pdrintensitieskkms_new.txt', format='A,A,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', $
  model,g0, co10, co21, co32, co43, co54, co65, co76, co87, co98, co109, co1110, co1211, co1312, co1413, co1514, co1615, co1716
  
;convert from K km/s to erg/s/cm^2/arcsec^2 to compare to data
  ;(((7.6373435*1e5)*2*k)/((2600.76/1e4)^3))/4.25e10
  co10 = (((co10*1e5)*2*k)/((line_center_co[0]/1e4)^3))/4.25e10
  co21 = (((co21*1e5)*2*k)/((line_center_co[1]/1e4)^3))/4.25e10
  co32 = (((co32*1e5)*2*k)/((line_center_co[2]/1e4)^3))/4.25e10
  co43 = (((co43*1e5)*2*k)/((line_center_co[3]/1e4)^3))/4.25e10
  co54 = (((co54*1e5)*2*k)/((line_center_co[4]/1e4)^3))/4.25e10
  co65 = (((co65*1e5)*2*k)/((line_center_co[5]/1e4)^3))/4.25e10
  co76 = (((co76*1e5)*2*k)/((line_center_co[6]/1e4)^3))/4.25e10
  co87 = (((co87*1e5)*2*k)/((line_center_co[7]/1e4)^3))/4.25e10
  co98 = (((co98*1e5)*2*k)/((line_center_co[8]/1e4)^3))/4.25e10
  co109 = (((co109*1e5)*2*k)/((line_center_co[9]/1e4)^3))/4.25e10
  co1110 = (((co1110*1e5)*2*k)/((line_center_co[10]/1e4)^3))/4.25e10
  co1211 = (((co1211*1e5)*2*k)/((line_center_co[11]/1e4)^3))/4.25e10
  co1312 = (((co1312*1e5)*2*k)/((line_center_co[12]/1e4)^3))/4.25e10
  co1413 = (((co1413*1e5)*2*k)/((line_center_co[13]/1e4)^3))/4.25e10
  co1514 = (((co1514*1e5)*2*k)/((line_center_co[14]/1e4)^3))/4.25e10
 ; co1615 = (((co1615*1e5)*2*k)/((line_center_co[15]/1e4)^3))/4.25e10
 ; co1716 = (((co1716*1e5)*2*k)/((line_center_co[16]/1e4)^3))/4.25e10

;make column labels with model names
leg_model1 = 'PDR: '+ strcompress(model[0]) +' ' + strcompress(g0[0])
leg_model2 = 'PDR: '+ strcompress(model[1]) +' ' +  strcompress(g0[1])
leg_model3 = 'PDR: '+ strcompress(model[2]) +' ' +  strcompress(g0[2])
leg_model4 = 'PDR: '+ strcompress(model[3]) +' ' +  strcompress(g0[3])
leg_model5 = 'PDR: '+ strcompress(model[4]) +' ' +  strcompress(g0[4])
leg_model6 = 'PDR: '+ strcompress(model[5]) +' ' +  strcompress(g0[5])


;make new file to plot models in converted units
;switch from row based to column based
;col 1 = model 1, col 2 = model 2, etc
  openw, lun, moddir + 'pdrintensitiescgs_new.txt', /get_lun
 ; printf, format='A,A,A,A,A,A', model1, model2, model3, model4, model5, model6
    for i=0, n_elements(model) -1 do begin
      printf, lun, format='(2(a16,2X),15(g16.10,2X))', $
        model[i],g0[i], co10[i], co21[i], co32[i], co43[i], co54[i], co65[i], co76[i], co87[i], $
          co98[i], co109[i], co1110[i], co1211[i], co1312[i], co1413[i], co1514[i];, co1615[i], co1716[i]
    endfor
 free_lun, lun
 close, lun
 
;read in models in correct units then arrange in format where col1 = model1, col2 = model2, etc 
 readcol, moddir + 'pdrintensitiescgs_new.txt', format='A,A,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', $;,D,D', $
   model,g0, co10, co21, co32, co43, co54, co65, co76, co87, co98, co109, co1110, co1211, co1312, co1413, co1514;, co1615, co1716
  
;large array of all models 
  set = [ [co10], [co21], [co32], [co43], [co54], [co65], [co76], [co87], [co98], [co109], [co1110],$
           [co1211], [co1312], [co1413], [co1514] ] ;, [co1615], [co1716] ]

; make array of each model   
  model1 = set[0,*]
    model1 = reverse(model1)
  model2 = set[1,*]
    model2 = reverse(model2)
  model3 = set[2,*]
    model2 = reverse(model2)  
  model4 = set[3,*]
    model2 = reverse(model2)
  model5 = set[4,*]
    model2 = reverse(model2)
  model6 = set[5,*] 
    model2 = reverse(model2)

;get shock model information
  readcol, '~/models/intensitiesergscmarcsec.txt', format = 'D,D,D,D,D,D,D,D,D', $
    n25v2b1, n25v3b1, n25v3b3, n30v2b1, n30v3b1, n30v3b3, n35v2b1, n35v3b1, n35v3b3
     
      n25v2b1 = reverse(n25v2b1) ;reverse so that the line_center_co and line_name_co are in right order
      n25v3b1 = reverse(n25v3b1) ;goes from CO15-14 to CO1-0
      n25v3b3 = reverse(n25v3b3)
      n30v2b1 = reverse(n30v2b1)
      n30v3b1 = reverse(n30v3b1)
      n30v3b3 = reverse(n30v3b3)
      n35v2b1 = reverse(n35v2b1)
      n35v3b1 = reverse(n35v3b1)
      n35v3b3 = reverse(n35v3b3)
  
;; read in L1527 to plot against new pdrmodels
    readcol, indir + object + '_continuumsubtracted_'+ molecule + '_lines.txt', format = 'A, F, D, D, D, D, D, D, D, D, D, D, D, F, D, D, A, A, A', $
      line_name_n, lab_wl_n, cen_wl_n, sig_cen_wl_n, str_n, sig_str_n, fwhm_n, sig_fwhm_n, base_str_n, noise_n, snr_n, $
      E_u_n, A_n, g_n, ra_n, dec_n, pix_n, blend_msg_all, lowest
    
      str_n=str_n*1e7  ;convert units to erg s^-1 cm^-2 arcsec ^-2
      sig_str_n = sig_str_n*1e7
      base_str_n = base_str_n*1e7
      noise_n = noise_n*1e7
    ;labels
      nonoutflows = where(pix_n eq 'averagenooutflow')
    ;for nonoutflow
      wl_non = cen_wl_n[nonoutflows]
      str_non = str_n[nonoutflows]
      str_non = str_non[1:3]
      nondetect = where(pix_n eq 'averagenooutflow' and snr_n le 3.0) ; should be [0], [4]
      upperlimits = noise_n[nondetect]*(fwhm_n[nondetect]/2.354)*sqrt(2*!pi)*3 ;to get the 3sigma upperlimits take baseline*(fwhm/2.354)*sqrt2pi*3
      wl_nondetect = cen_wl_n[nondetect]
    ;error bars
      nondy = sig_str_n[nonoutflows]
      nondy = nondy[1:3]
    print, 'averagenooutflows ', wl_non, str_non, nondy ;check

;Plot of best fit model: n30v3b1 vs each new pdr model
    set_plot, 'ps'
    loadct,12,/silent
    device, filename= plotdir + object + '_pdrcompare_n30v3b1.eps', encapsulated=eps, $
      /helvetica, /isolatin1, /landscape, /color, font_size=10, decomposed=0 ;sizing to match model plots
    ;Make plot symbols circles
      plotsym,0,1, /fill ;circle
      x = alog10(wl_non[1:3])
      y = str_non
    ;plotmodels
      plot, x, y,xrange=[2.2,3.5], yrange=[1e-22,1e-16], xthick=4, ythick=4, thick=4, charsize=1.1, charthick=4, XTITLE =xlabel, YTITLE = ylabel, /ylog, ystyle=1, xstyle=1, /nodata
      oplot, alog10(line_center_co), n30v3b1, psym=8, color=colors[1] ;points
      oplot, alog10(line_center_co), n30v3b1, color=colors[1], thick=4 ;line
      oplot, reverse(alog10(line_center_co)), model4, psym=8, color=colors[2] ;points
      oplot, reverse(alog10(line_center_co)), model4, color=colors[2], thick=4 ;line
      oplot, x, y, PSYM=8, symsize=1, color=colors[0] ;green
      ;nondetection points
      plotsym,5,1, /fill ;down triangle
      oplot, alog10(wl_nondetect), upperlimits, PSYM=8, symsize=1,color=colors[0]
      nondyup = nondy
      nondydown = nondy
      oploterror, x, y, nondyup, psym=3, ERRcolor=colors[0], /hibar
      oploterror, x, y, nondydown, psym=3, ERRcolor=colors[0],  /lobar
      al_legend, ['n = 10!u3.0!n cm!u-3!n','v = 3 km s!u-1!n', 'B = 4 ' + Greek('mu', /append_font) + 'G', leg_model4], charthick=4, /right, box=0
      al_legend, legend ,textcolors = colors,/bottom, charthick=4, /right
    device, /close_file,decomposed=1
    set_plot, 'x'

stop

end