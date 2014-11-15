PRO evol_indicators_norb,pacs_file,irs_file, photom_file,distance,ttl,spire_file,indir=indir,outdir=outdir,outeps=outeps,usertbol=usertbol,userlbol=userlbol,em=em,$
  nosh=nosh,yrange=yrange,xrange=xrange,leg_right=leg_right,nolegend=nolegend,nobb=nobb

;CALLING SEQUENCE
;IDL> evol_indicators,'path_to_sed.txt',distance_in_pc,
;'file_for_plot.ps','title_of_source'
;
;Program to calculate the following evolutionary indicators:
;
;     Lbol
;     Tbol
;     Lbol/Lsmm
;
;The program will look for a 2-column file giving wavelength in 
;microns and flux density in Jy. A third column for uncertainty is currently ignored, and the program will work whether or not the column is present
;in the input file.
;
;All integrations are done using TSUM.pro (former
;integrate_trapezoid.pro is the same).  This requires
;the x-array (nu for all integrations) to be in increasing order.
;Thus the program will check that wavelength is in decreasing order.
;If not, it will re-order.
;
;
;Flux <= 0 points are removed.
;Distance is input by the user in pc
;No errors in returned quantities
;A line containing the source (title), Tbol, and Lbol will be printed.
;A plot of the SED will be displayed and saved under 'plotfile'
;
;Code last modified 09/19/2011 by MRR
;JDG - 10/02/2012; added ability to remove PACS data from fit, added SPIRE data keyword/optional
;JDG heavily modified from Michelle's original code -- can remove IRS from fit, added WISE, FORCAST keywords

;---Read in contents of files.  If there are lines that should be
;omitted, begin them with a %.

;if no indir set, set it
if not keyword_set(indir) then indir='~/jeong-eun/michelle/SED_final_cal/'

;PACS
if pacs_file ne "none" then begin
  readcol,indir+pacs_file,wvp,sjyp,format='D,D', comment='%',/silent
endif

;IRS
if irs_file ne "none" then begin
    readcol,indir+irs_file, wvi,sjyi, format='D,D',comment='%',/silent
    irs_sort=sort(wvi) ;Sort for plotting
    wvi=wvi(irs_sort)
    sjyi=sjyi(irs_sort)
endif

;SPIRE
if spire_file ne "none" then begin
  readcol,indir+spire_file,wvs,sjys,format='D,D',comment='%',/silent
endif

;Photometry
readcol, indir+photom_file, wvph, sjyph, esjyph, inst, format='D,D,D,A', comment='%', /silent

;Combine IRS+Photometry
if irs_file eq "none" then begin
    wvnoh=wvph                    ;noh = no herschel
    sjynoh=sjyph
endif else begin
    wvnoh=[wvi, wvph]
    sjynoh=[sjyi, sjyph]
endelse

;Add PACS
if pacs_file ne "none" then begin
  wv=[wvnoh, wvp]
  sjy=[sjynoh, sjyp]
endif else begin
  wv=wvnoh
  sjy=sjynoh
endelse

;Add SPIRE
if spire_file ne "none" then begin
wv=[wv,wvs]
sjy=[sjy,sjys]
endif

c=2.99792d+14
pi=double(!pi)
lsun=3.826d33
;lsun=3.86d33
dcgs=distance*3.0857d18                   ;Convert distance from pc to cm
h=6.626*1.0d-27
k=1.38*1.0d-16
ccm=c*1.0d-4

;convert flux density to Jy
;sjy=smJy*1.0e-3

;check that wavelength is in decreasing order
;if not, re-order (just reverses order)
;if(wv[0] lt wv[sz-1]) then begin
;    temp_wv=wv
;    temp_sjy=sjy
;    for i=0,sz-1 do begin
;        wv[i]=temp_wv[sz-1-i]
;        sjy[i]=temp_sjy[sz-1-i]
;    endfor
;endif

;re-order array that may be in random order to be in decreasing order
temp_wv=wv
temp_sjy=sjy
slist=sort(-1*temp_wv)
wv=temp_wv(slist)
sjy=temp_sjy(slist)

slistnoh=sort(-1*wvnoh)
wvnoh=wvnoh(slistnoh)
sjynoh=sjynoh(slistnoh)

;remove zero (and negative) fluxes
w=where(sjy gt 0.0)
wv=wv[w]
sjy=sjy[w]

wnoh=where(sjynoh gt 0.0)
wvnoh=wvnoh(wnoh)
sjynoh=sjynoh(wnoh)

wvcgs=wv*(1.0e-4)                         ;Convert wavelength to cgs
nu=c/wv                                   ;Convert wavelength to frequency
scgs=sjy*(1.0e-23)                        ;Convert Jy to cgs

wvcgsnoh=wvnoh*(1.0e-4)
nunoh=c/wvnoh
scgsnoh=sjynoh*(1.0e-23)

if pacs_file ne "none" then begin
  wvcgsp=wvp*(1.0e-4) ;pacs         
  nup=c/wvp                                  
  scgsp=sjyp*(1.0e-23)
endif

if spire_file ne "none" then begin
  wvcgss=wvs*(1.0e-4) ;pacs         
  nus=c/wvs                                  
  scgss=sjys*(1.0e-23)
endif

if irs_file ne "none" then begin
    wvcgsi=wvi*(1.0e-4)           
    nui=c/wvi                                  
    scgsi=sjyi*(1.0e-23)
endif

wvcgsph=wvph*(1.0e-4)           
nuph=c/wvph                                  
scgsph=sjyph*(1.0e-23)

;------------------------------------------------------------------
;-------------------------Lbol calculation-------------------------
;------------------------------------------------------------------


temp=tsum(nu,scgs)
lbol=(double(4.)*pi*dcgs*dcgs*temp)/lsun

if n_elements(scgsnoh) GT 1 then begin
  tempnoh=tsum(nunoh, scgsnoh)
  lbol_noh=(double(4.)*pi*dcgs*dcgs*tempnoh)/lsun
endif else begin
  tempnoh=temp
  lbol_noh=lbol
  print,'No non-Herschel data, so just avoiding errors...'
endelse

;------------------------------------------------------------------
;-------------------------Lbol/Lsmm calculation--------------------
;------------------------------------------------------------------

w=where(wv ge 350.0)
if(w[0] eq -1) or n_elements(w) EQ 1 then begin
   lbolsmm=-999
endif else begin
   nu_lsmm=nu[w]
   scgs_lsmm=scgs[w]
   sz_lsmm=n_elements(nu_lsmm)
   temp=tsum(nu_lsmm,scgs_lsmm)
   lsmm=(double(4.)*pi*dcgs*dcgs*temp)/lsun
   lbolsmm=lbol/lsmm
   ;lbolsmm_noh=lbol_noh/lsmm
endelse
;
print,lbolsmm


;------------------------------------------------------------------
;-------------------------Tbol calculation-------------------------
;------------------------------------------------------------------

moment_first=tsum(nu,nu*scgs)
moment_zeroth=tsum(nu,scgs)
nubar=moment_first/moment_zeroth
tbol=1.25e-11*nubar

mf_noh=tsum(nunoh, nunoh*scgsnoh) ; Same Tbol calc without herschel
mz_noh=tsum(nunoh, scgsnoh)
nubar_noh=mf_noh/mz_noh
tbol_noh=1.25e-11*nubar_noh


;print results
;results=[lbol,tbol,lbolsmm]
;print, ttl, '  Tbol:',tbol,'  Lbol:',lbol, '  Lbol/Lsubmm:',lbolsmm
;print, 'W/O Herschel    Tbol:',tbol_noh,'  Lbol:',lbol_noh, '    Lbol/Lsubmm:', lbolsmm_noh


;------------------------------------------------------------------
;-------------------------------Plot-------------------------------
;------------------------------------------------------------------


if keyword_set(usertbol) then tbol=usertbol
if keyword_set(userlbol) then lbol=userlbol


tbolsf=strnsignif(tbol,3)
lbolsf=strnsignif(lbol,3)

tbolsf_noh=strnsignif(tbol_noh,3)
lbolsf_noh=strnsignif(lbol_noh,3)

;nutex=TextoIDL('log \nuS_{\nu}(erg/s/cm^2)')
nutex='Log('+Greek('nu')+'S!D'+Greek('nu')+'!Nerg s!U-1!N cm!U2!N)'

;wavtex=TextoIDL('log \lambda(\mum)')
wavtex='Log '+Greek('lambda')+' ('+Greek('mu')+'m)'

;print, textoidl('\mu')

; Values for ploting
nufnu=alog10(nu*sjy*1.0d-23) ;Y values, nu*F(nu)
wvlog=alog10(wv) ; X values, Log(wavelength)
maxy=max(nufnu)+0.2 ;Max Y value for plotting
miny=min(nufnu)-0.2 ;Min Y value for plotting
nufnu_noh=alog10(nunoh*sjynoh*1.0d-23)
wvlog_noh=alog10(wvnoh)

if pacs_file ne "none" then begin
  nufnup=alog10(nup*sjyp*1.0d-23)
  wvlogp=alog10(wvp)
endif

if spire_file ne "none" then begin
  nufnus=alog10(nus*sjys*1.0d-23)
  wvlogs=alog10(wvs)
endif

if irs_file ne "none" then begin
    nufnui=alog10(nui*sjyi*1.0d-23)
    wvlogi=alog10(wvi)
endif
nufnuph=alog10(nuph*sjyph*1.0d-23)
wvlogph=alog10(wvph)

;maxx=max(wvlog)+0.1 ;Max X value for plotting
if min(wvlog) ge 0.5 then begin ; Min X value for plotting. 0 or less.
  minx=0.5
  maxx=4.0
endif else begin
    if min(wvlog) ge 0 then begin
        minx=0.0
        if max(wvlog) ge 3 then maxx=4.0
        if max(wvlog) lt 3 then maxx=3.0
    endif else begin
        if min(wvlog) le 0 then begin
            minx=-0.5
            maxx=4.0
        endif
    endelse
endelse

;plot, [minx, maxx],[miny,maxy],title=ttl, xtitle=wavtex, ytitle=nutex,/nodata
;oplot, wvlog, nufnu, psym=1
;oplot, wvlog_noh, nufnu_noh, psym=6, color=cgcolor('blue')

xwv=(findgen(1.0d4)+1) ;Make 10000 points for BB plot, start at 1.

const=(lbol*lsun/(4*!pi*dcgs^2))*(15*ccm^2*h^3/(2*!pi^4*k^4*tbol^4))
ynufnu=const*(2*h*ccm*ccm)/(((xwv*1.0d-4)^4)*(exp(h*ccm/(xwv*1.0d-4*k*tbol))-1))
if keyword_set(em) then ynufnu_em=ynufnu*(250/xwv)^em
xlogwv=alog10(xwv)
ylognufnu=alog10(ynufnu)
if keyword_set(em) then ylognufnu=alog10(ynufnu_em)

const_noh=(lbol_noh*lsun/(4*!pi*dcgs^2))*(15*ccm^2*h^3/(2*!pi^4*k^4*tbol_noh^4))
ynufnu_noh=const_noh*(2*h*ccm*ccm)/(((xwv*1.0d-4)^4)*(exp(h*ccm/(xwv*1.0d-4*k*tbol_noh))-1))
ylognufnu_noh=alog10(ynufnu_noh)



;oplot, xlogwv, ylognufnu, linestyle=0
;oplot, xlogwv, ylognufnu_noh, linestyle=0, color=cgcolor('blue')

;JDG removed this line
;!PATH = Expand_Path('+~/coyoteprograms') + ':' + !PATH

;al_legend,['T!Dbol!N='+tbolsf+' K','('+tbolsf_noh+' K)','L!Dbol!N='+lbolsf+' L!I!9n!3!N','('+lbolsf_noh+' K)'],/right_legend,/top_legend, charsize=1.5, charthick=2, box=0,textcolors=[cgcolor('white'),cgcolor('blue'),cgcolor('white'),cgcolor('blue')]
;al_legend,'hi',/right_legend,/top_legend, charsize=1.5, charthick=2, box=0
;legend,'L!Dbol!N='+lbolsf+' L!I!9n!3!N',position=[0.7,0.85],/norm, charsize=1.5, charthick=2, box=0



;;; Separate 2MASS, IRAS, IRAC, MIPS, ISO, SPIRE, PACS

phot_len=n_elements(inst)-1

twomass_wv=[0]
iras_wv=[0]
irac_wv=[0]
mips_wv=[0]
iso_wv=[0]
spire_wv=[0]
pacs_wv=[0]
forcast_wv=[0]
wise_wv=[0]
other_wv=[0]

twomass_fx=[0]
iras_fx=[0]
irac_fx=[0]
mips_fx=[0]
iso_fx=[0]
spire_fx=[0]
pacs_fx=[0]
forcast_fx=[0]
wise_fx=[0]
other_fx=[0]


for i=0, phot_len DO BEGIN
    if inst(i) eq '2MASS' THEN BEGIN
        twomass_wv=[twomass_wv, wvlogph(i)]
        twomass_fx=[twomass_fx, nufnuph(i)]
    endif else begin
        if inst(i) eq 'IRAS' then begin
            iras_wv=[iras_wv, wvlogph(i)]
            iras_fx=[iras_fx, nufnuph(i)]
        endif else begin
            if inst(i) eq 'IRAC' then begin
                irac_wv=[irac_wv, wvlogph(i)]
                irac_fx=[irac_fx, nufnuph(i)]
            endif else begin
                if inst(i) eq 'MIPS' then begin
                    mips_wv=[mips_wv, wvlogph(i)]
                    mips_fx=[mips_fx, nufnuph(i)]
                endif else begin
                    if inst(i) eq 'ISO' then begin
                        iso_wv=[iso_wv, wvlogph(i)]
                        iso_fx=[iso_fx, nufnuph(i)]
                    endif else begin
                      if inst(i) eq 'SPIRE' then begin
                          spire_wv=[spire_wv, wvlogph(i)]
                          spire_fx=[spire_fx, nufnuph(i)]
                      endif else begin
                        if inst(i) eq 'PACS' or inst(i) eq 'PACSLine' then begin
                          pacs_wv=[pacs_wv, wvlogph(i)]
                          pacs_fx=[pacs_fx, nufnuph(i)]
                        endif else begin
                          if inst(i) eq 'FORCAST' then begin
                            forcast_wv=[forcast_wv, wvlogph(i)]
                            forcast_fx=[forcast_fx, nufnuph(i)]
                          endif else begin
                              if inst(i) eq 'WISE' then begin
                                wise_wv=[wise_wv, wvlogph(i)]
                                wise_fx=[wise_fx, nufnuph(i)]
                              endif else begin                             
                                other_wv=[other_wv, wvlogph(i)]
                                other_fx=[other_fx, nufnuph(i)]
                            endelse
                          endelse
                        endelse
                      endelse
                    endelse
                endelse
            endelse
         endelse
     endelse    
 endfor

;set outdir -- JDG
if not keyword_set(outdir) then outdir=indir
if not keyword_set(outeps) then begin
  plotfile=outdir+'SED_plots/no_rb_factor/with_blackbody/'+ttl+'.eps'
  plotfile2=outdir+'SED_plots/no_rb_factor/'+ttl+'.eps'
endif else begin
  plotfile=outdir+'SED_plots/no_rb_factor/with_blackbody/'+outeps
  plotfile2=outdir+'SED_plots/no_rb_factor/'+outeps
endelse

set_plot, 'ps'
device, bits_per_pixel=8, color=1,/portrait, /encapsul,filename=plotfile
loadct,0
;miny=-14 ;special for L1448-C

plot, [minx, maxx],[miny,maxy],title=ttl, xtitle=wavtex, ytitle=nutex,/nodata, xthick=8, ythick=8, thick=8,charsize=1.25,charthick=5, xtickinterval=0.5, ytickinterval=1, xstyle=1, ystyle=1,xcharsize=1.25,ycharsize=1.25,yr=yrange,xr=xrange
leg_inst=['Start'] ;Start the legend instrument array
leg_sym=[0]      ;Start the legend symbol array
leg_size=[0]     ;Start the legend symbol size
leg_thick=[0]      ;Start the legend thickness



; Make filled square
A = [-1,-1,1,1]
B = [-1,1,1,-1]

USERSYM, A, B, /Fill

;
;if n_elements(twomass_wv) ne 1 then begin
;    twomass_wv=twomass_wv[1:(n_elements(twomass_wv)-1)]
;    twomass_fx=twomass_fx[1:(n_elements(twomass_fx)-1)]
;    oplot, twomass_wv, twomass_fx, psym=8, symsize=1.2
;    leg_inst=[leg_inst, '2MASS']
;    leg_sym=[leg_sym, 15]
;    leg_size=[leg_size, 1.2]
;    leg_thick=[leg_thick, 3]
;endif
;
;; Make filled diamond
X = [-1, 0, 1, 0, -1]  
Y = [0, 1, 0, -1, 0]  
USERSYM, X, Y ,/Fill
;
;if n_elements(irac_wv) ne 1 then begin
;    irac_wv=irac_wv[1:(n_elements(irac_wv)-1)]
;    irac_fx=irac_fx[1:(n_elements(irac_fx)-1)]
;    oplot, irac_wv, irac_fx, psym=8, symsize=1.2
;    leg_inst=[leg_inst, 'IRAC']
;    leg_sym=[leg_sym, 14]
;    leg_size=[leg_size, 1.2]
;    leg_thick=[leg_thick, 3]
;endif

; Make a vector of 16 points, A[i] = 2pi/16:  
A = FINDGEN(17) * (!PI*2/16.)  
; Define the symbol to be a unit circle with 16 points,   
; and set the filled flag:  
USERSYM, COS(A), SIN(A), /FILL  

if n_elements(other_wv) ne 1 then begin
    other_wv=other_wv[1:(n_elements(other_wv)-1)]
    other_fx=other_fx[1:(n_elements(other_fx)-1)]
    oplot, other_wv, other_fx, psym=8, symsize=1.2
endif

if keyword_set(em) then begin
  ;emissivity normalized
  oplot, xlogwv, ylognufnu
  print, 'emissivity law = lambda^'+trim(em,2)
endif else begin
  ;regular (default, actual blackbody) -- disable for NJE/BOV plot
  if not keyword_set(nobb) then oplot, xlogwv, ylognufnu, linestyle=2, thick=3
endelse


if pacs_file ne "none" then begin
  ;;PACS spec
  ;;--split at 100um (log(100)=2)--
  lop=where(wvlogp lt 2.0)
  hip=where(wvlogp gt 2.0)
  wvlogp_lo=wvlogp(lop)
  nufnup_lo=nufnup(lop)
  wvlogp_hi=wvlogp(hip)
  nufnup_hi=nufnup(hip)
endif

;oplot, wvlogp,    nufnup,    thick=3,  color=cgcolor('pink')
if spire_file ne "none" then begin
  if pacs_file ne "none" then begin
    oplot, wvlogp_lo, nufnup_lo, thick=2,  color=cgcolor('green')
    oplot, wvlogp_hi, nufnup_hi, thick=2,  color=cgcolor('green')
  endif else oplot,wvlogs,nufnus,thick=2,  color=cgcolor('red')
endif else begin
  if pacs_file ne "none" then begin
    oplot, wvlogp_lo, nufnup_lo, thick=2,  color=cgcolor('red')
    oplot, wvlogp_hi, nufnup_hi, thick=2,  color=cgcolor('red')
  endif
endelse
  
;;SPIRE spec (overplot the SPIRE data in a different color, for log(195)=2.29
if spire_file ne "none" then begin
;  vhip=where(wvlogp gt 2.29)
;  wvlogp_vhi=wvlogp(vhip)
;  nufnup_vhi=nufnup(vhip)
  oplot,wvlogs,nufnus, color=cgcolor('red')
endif

;;IRS spec
if irs_file ne "none" then begin
    if keyword_set(nosh) then begin
      oplot,wvlogi[where(wvlogi LT alog10(14))], nufnui[where(wvlogi LT alog10(14))],thick=2,  color=cgcolor('blue')
      oplot,wvlogi[where(wvlogi GT alog10(18))], nufnui[where(wvlogi GT alog10(18))],thick=2,  color=cgcolor('blue')
    endif else  oplot, wvlogi, nufnui, thick=3,  color=cgcolor('blue')
endif


; Make filled square
A = [-1,-1,1,1]
B = [-1,1,1,-1]
USERSYM, A, B, /Fill

if n_elements(twomass_wv) ne 1 then begin
    twomass_wv=twomass_wv[1:(n_elements(twomass_wv)-1)]
    twomass_fx=twomass_fx[1:(n_elements(twomass_fx)-1)]
    oplot, twomass_wv, twomass_fx, psym=8, symsize=1.2
    leg_inst=[leg_inst, '2MASS']
    leg_sym=[leg_sym, 15]
    leg_size=[leg_size, 1.2]
    leg_thick=[leg_thick, 3]
endif

;; Make filled diamond
X = [-1, 0, 1, 0, -1]  
Y = [0, 1, 0, -1, 0]  
USERSYM, X, Y ,/Fill

if n_elements(irac_wv) ne 1 then begin
    irac_wv=irac_wv[1:(n_elements(irac_wv)-1)]
    irac_fx=irac_fx[1:(n_elements(irac_fx)-1)]
    oplot, irac_wv, irac_fx, psym=8, symsize=1.2
    leg_inst=[leg_inst, 'IRAC']
    leg_sym=[leg_sym, 14]
    leg_size=[leg_size, 1.2]
    leg_thick=[leg_thick, 3]
endif

if n_elements(mips_wv) ne 1 then begin
    mips_wv=mips_wv[1:(n_elements(mips_wv)-1)]
    mips_fx=mips_fx[1:(n_elements(mips_fx)-1)]
    oplot, mips_wv, mips_fx, psym=6, symsize=1.2, thick=3
    leg_inst=[leg_inst, 'MIPS']
    leg_sym=[leg_sym, 6]
    leg_size=[leg_size, 1.2]
    leg_thick=[leg_thick, 3]
endif

if n_elements(forcast_wv) ne 1 then begin
    forcast_wv=forcast_wv[1:(n_elements(forcast_wv)-1)]
    forcast_fx=forcast_fx[1:(n_elements(forcast_fx)-1)]
    oplot, forcast_wv, forcast_fx, psym=2, symsize=1.2, thick=3
    leg_inst=[leg_inst, 'FORCAST']
    leg_sym=[leg_sym, 2]
    leg_size=[leg_size, 1.2]
    leg_thick=[leg_thick, 3]
endif


if n_elements(wise_wv) ne 1 then begin
    wise_wv=wise_wv[1:(n_elements(wise_wv)-1)]
    wise_fx=wise_fx[1:(n_elements(wise_fx)-1)]
    oplot, wise_wv, wise_fx, psym=1, symsize=1.2, thick=3
    leg_inst=[leg_inst, 'WISE']
    leg_sym=[leg_sym, 1]
    leg_size=[leg_size, 1.2]
    leg_thick=[leg_thick, 3]
endif


if n_elements(iras_wv) ne 1 then begin
    iras_wv=iras_wv[1:(n_elements(iras_wv)-1)]
    iras_fx=iras_fx[1:(n_elements(iras_fx)-1)]
    oplot, iras_wv, iras_fx, psym=5, symsize=1.2, thick=3
    leg_inst=[leg_inst, 'IRAS']
    leg_sym=[leg_sym, 5]
    leg_size=[leg_size, 1.2]
    leg_thick=[leg_thick, 3]
endif

if n_elements(iso_wv) ne 1 then begin
    iso_wv=iso_wv[1:n_elements(iso_wv)-1]
    iso_fx=iso_fx[1:n_elements(iso_fx)-1]
    oplot, iso_wv, iso_fx, psym=4, symsize=1.2, thick=3
    leg_inst=[leg_inst, 'ISO']
    leg_sym=[leg_sym, 4]
    leg_size=[leg_size, 1.2]
    leg_thick=[leg_thick, 3]
endif

if n_elements(spire_wv) ne 1 then begin
    spire_wv=spire_wv[1:n_elements(spire_wv)-1]
    spire_fx=spire_fx[1:n_elements(spire_fx)-1]
    oplot, spire_wv, spire_fx, psym=7, symsize=1.2, thick=3
    leg_inst=[leg_inst, 'SPIRE']
    leg_sym=[leg_sym, 7]
    leg_size=[leg_size, 1.2]
    leg_thick=[leg_thick, 3]
endif

; Make a vector of 16 points, A[i] = 2pi/16:  
A = FINDGEN(17) * (!PI*2/16.)  
; Define the symbol to be a unit circle with 16 points,   
; and set the filled flag:  
USERSYM, COS(A), SIN(A), /FILL   
if n_elements(pacs_wv) ne 1 then begin
    pacs_wv=pacs_wv[1:n_elements(pacs_wv)-1]
    pacs_fx=pacs_fx[1:n_elements(pacs_fx)-1]
    oplot, pacs_wv, pacs_fx, psym=7, symsize=1.2, thick=3
    leg_inst=[leg_inst, 'PACS']
    leg_sym=[leg_sym, 7]
    leg_size=[leg_size, 1.2]
    leg_thick=[leg_thick, 3]
endif


len_leg=n_elements(leg_inst)

if len_leg gt 1 then begin
    leg_inst=[leg_inst[1:len_leg-1], 'Other']
    leg_sym=[leg_sym[1:len_leg-1], 16]
    leg_size=[leg_size[1:len_leg-1], 1.2]
    leg_thick=[leg_thick[1:len_leg-1], 3]
    leg_line=[indgen(len_leg)*0, -99, -99]
endif else begin
    leg_inst=['Other']
    leg_sym=[16]
    leg_size=[1.2]
    leg_thick=[3]
    leg_line=[0,-99,-99]
endelse

;empty strings in line below for NJE/BOV plot
;leg_inst=[leg_inst, '','']
leg_inst=[leg_inst,'T!Dbol!N='+tbolsf+' K','L!Dbol!N='+lbolsf+' L!I!9n!3!N']
leg_sym=[leg_sym, 0, 0]
leg_size=[leg_size, 1.2, 1.2]
leg_thick=[leg_thick, 3, 3]

if keyword_set(leg_right) then begin
  if not keyword_set(nolegend) then al_legend,leg_inst,psym=leg_sym, /bottom, /right, charthick=3, symsize=leg_size, thick=leg_thick, linestyle=leg_line, bthick=5,charsize=1.25
endif else if not keyword_set(nolegend) then al_legend,leg_inst,psym=leg_sym, /bottom, /center, charthick=3, symsize=leg_size, thick=leg_thick, linestyle=leg_line, bthick=5,charsize=1.25
;al_legend,['2MASS','IRAC','MIPS','IRAS','ISO','Other'],psym=[15,14,6,5,4,16], /bottom, /center, charthick=3, symsize=[1.2,1.2,1.2,1.2,1.2,1.2], thick=[3,3,3,3,3,3]

;al_legend,['T!Dbol!N='+tbolsf+' K','L!Dbol!N='+lbolsf+' L!I!9n!3!N'],/right_legend,/top_legend, charsize=1.5, charthick=4, box=1

;print, minx, maxx, miny, maxy

set_plot, 'ps'
device, bits_per_pixel=8, color=1,/portrait, /encapsul,filename=plotfile2
plot, [minx, maxx],[miny,maxy],title=ttl, xtitle=wavtex, ytitle=nutex,/nodata, charthick=5, xthick=8, ythick=8, xtickinterval=0.5, ytickinterval=1, xstyle=1, ystyle=1,thick=8,charsize=1.25,xcharsize=1.25,ycharsize=1.25,yr=yrange,xr=xrange

; Make filled square
A = [-1,-1,1,1]
B = [-1,1,1,-1]

USERSYM, A, B, /Fill
;
;
;if n_elements(twomass_wv) ne 1 then begin
;    twomass_wv=twomass_wv[1:(n_elements(twomass_wv)-1)]
;    twomass_fx=twomass_fx[1:(n_elements(twomass_fx)-1)]
;    oplot, twomass_wv, twomass_fx, psym=8, symsize=1.2
;    leg_inst=[leg_inst, '2MASS']
;    leg_sym=[leg_sym, 15]
;    leg_size=[leg_size, 1.2]
;    leg_thick=[leg_thick, 3]
;    leg_line=[leg_line,-99]
;endif
;
;;; Make filled diamond
;X = [-1, 0, 1, 0, -1]  
;Y = [0, 1, 0, -1, 0]  
;USERSYM, X, Y ,/Fill
;
;if n_elements(irac_wv) ne 1 then begin
;    irac_wv=irac_wv[1:(n_elements(irac_wv)-1)]
;    irac_fx=irac_fx[1:(n_elements(irac_fx)-1)]
;    oplot, irac_wv, irac_fx, psym=8, symsize=1.2
;    leg_inst=[leg_inst, 'IRAC']
;    leg_sym=[leg_sym, 14]
;    leg_size=[leg_size, 1.2]
;    leg_thick=[leg_thick, 3]
;    leg_line=[leg_line,-99];,-99]
;endif

; Make a vector of 16 points, A[i] = 2pi/16:  
A = FINDGEN(17) * (!PI*2/16.)  
; Define the symbol to be a unit circle with 16 points,   
; and set the filled flag:  
USERSYM, COS(A), SIN(A), /FILL  


if n_elements(other_wv) ne 1 then begin
    oplot, other_wv, other_fx, psym=8, symsize=1.2
endif

if spire_file ne "none" then begin
  if pacs_file ne "none" then begin
    oplot, wvlogp_lo, nufnup_lo, thick=3,  color=cgcolor('green')
    oplot, wvlogp_hi, nufnup_hi, thick=3,  color=cgcolor('green')
  endif else oplot,wvlogs,nufnus, thick=3, color=cgcolor('red')
endif else begin
  if pacs_file ne "none" then begin
    oplot, wvlogp_lo, nufnup_lo, thick=3,  color=cgcolor('red')
    oplot, wvlogp_hi, nufnup_hi, thick=3,  color=cgcolor('red')
  endif
endelse
  
if irs_file ne "none" then begin
    oplot, wvlogi, nufnui,  thick=3, color=cgcolor('blue')
endif



;;SPIRE spec (overplot the SPIRE data in a different color, for log(195)=2.29
if spire_file ne "none" then begin
;  vhip=where(wvlogp gt 2.29)
;  wvlogp_vhi=wvlogp(vhip)
;  nufnup_vhi=nufnup(vhip)
  oplot,wvlogs,nufnus, color=cgcolor('red')
endif
;
;
;if n_elements(mips_wv) ne 1 then begin
;    mips_wv=mips_wv[1:(n_elements(mips_wv)-1)]
;    mips_fx=mips_fx[1:(n_elements(mips_fx)-1)]
;    oplot, mips_wv, mips_fx, psym=6, symsize=1.2, thick=3
;    leg_inst=[leg_inst, 'MIPS']
;    leg_sym=[leg_sym, 6]
;    leg_size=[leg_size, 1.2]
;    leg_thick=[leg_thick, 3]
;endif
;
;if n_elements(iras_wv) ne 1 then begin
;    iras_wv=iras_wv[1:(n_elements(iras_wv)-1)]
;    iras_fx=iras_fx[1:(n_elements(iras_fx)-1)]
;    oplot, iras_wv, iras_fx, psym=5, symsize=1.2, thick=3
;    leg_inst=[leg_inst, 'IRAS']
;    leg_sym=[leg_sym, 5]
;    leg_size=[leg_size, 1.2]
;    leg_thick=[leg_thick, 3]
;endif
;
;if n_elements(iso_wv) ne 1 then begin
;    iso_wv=iso_wv[1:n_elements(iso_wv)-1]
;    iso_fx=iso_fx[1:n_elements(iso_fx)-1]
;    oplot, iso_wv, iso_fx, psym=4, symsize=1.2, thick=3
;    leg_inst=[leg_inst, 'ISO']
;    leg_sym=[leg_sym, 4]
;    leg_size=[leg_size, 1.2]
;    leg_thick=[leg_thick, 3]
;endif
;
;if n_elements(spire_wv) ne 1 then begin
;    spire_wv=spire_wv[1:n_elements(spire_wv)-1]
;    spire_fx=spire_fx[1:n_elements(spire_fx)-1]
;    oplot, spire_wv, spire_fx, psym=7, symsize=1.2, thick=3
;    leg_inst=[leg_inst, 'SPIRE']
;    leg_sym=[leg_sym, 7]
;    leg_size=[leg_size, 1.2]
;    leg_thick=[leg_thick, 3]
;endif
;
;if n_elements(pacs_wv) ne 1 then begin
;    pacs_wv=pacs_wv[1:n_elements(pacs_wv)-1]
;    pacs_fx=pacs_fx[1:n_elements(pacs_fx)-1]
;    oplot, pacs_wv, pacs_fx, psym=8, symsize=1.2, thick=3
;    leg_inst=[leg_inst, 'PACS']
;    leg_sym=[leg_sym, 8]
;    leg_size=[leg_size, 1.2]
;    leg_thick=[leg_thick, 3]
;endif
;
;
;len_leg=n_elements(leg_inst)
;
;if len_leg gt 1 then begin
;    leg_inst=[leg_inst[1:len_leg-1], 'Other']
;    leg_sym=[leg_sym[1:len_leg-1], 16]
;    leg_size=[leg_size[1:len_leg-1], 1.2]
;    leg_thick=[leg_thick[1:len_leg-1], 3]
;    leg_line=[indgen(len_leg)*0, -99, -99]
;    ;leg_line=[indgen(len_leg)*0]
;endif else begin
;    leg_inst=['Other']
;    leg_sym=[16]
;    leg_size=[1.2]
;    leg_thick=[3]
;    leg_line=[0,-99,-99]
;    ;leg_line=[0]
;endelse



if n_elements(iso_wv) ne 1 then begin
    oplot, iso_wv, iso_fx, psym=4, symsize=1.2, thick=5
endif

if n_elements(mips_wv) ne 1 then begin
    oplot, mips_wv, mips_fx, psym=6, symsize=1.2, thick=5
endif

if n_elements(iras_wv) ne 1 then begin
    oplot, iras_wv, iras_fx, psym=5, symsize=1.2, thick=5
endif

    
oplot, wvlogph, nufnuph, psym=4, thick=3, symsize=1.2


if not keyword_set(nolegend) then al_legend,leg_inst,psym=leg_sym, /bottom, /center, charthick=3, symsize=leg_size, thick=leg_thick, linestyle=leg_line, bthick=5,charsize=1.25

;al_legend,['T!Dbol!N='+tbolsf+' K','L!Dbol!N='+lbolsf+' L!I!9n!3!N'],/right_legend,/top_legend, charsize=1.5, charthick=4, box=1

device, /close_file

end

pro fix_irsspec,scalefactor,indir=indir,infile=infile,outdir=outdir,outfile=outfile

if not keyword_set(indir) then indir='~/Desktop/jeong-eun/michelle/l1448mm/IRS/'
if not keyword_set(outdir) then outdir=indir
if not keyword_set(infile) then infile='L1448-MM_SH_avg.txt'
if not keyword_set(outfile) then outfile=strmid(infile,0,strlen(infile)-4)+'_scaledtoLH.txt'

readcol,indir+infile,wave,flux,comment='%',format='F,F'
writecol,outdir+outfile,wave,flux*scalefactor

end
