pro all_plot_v3
readcol,'~/programs/all_data.txt',format='A,D,D,D,D,D,D,D,D,D,D,D,D',name,a,b,c,d,e,f,g,h,i,j,k,l
data = dblarr(12,n_elements(a))
data[0,*]=a & data[1,*]=b & data[2,*]=c & data[3,*]=d & data[4,*]=e & data[5,*]=f & data[6,*]=g & data[7,*]=h & data[8,*]=i & data[9,*]=j
data[10,*]=k & data[11,*]=l
readcol,'~/programs/all_err.txt',format='A,D,D,D,D,D,D,D,D,D,D,D,D',name1,a,b,c,d,e,f,g,h,i,j,k,l
err = dblarr(12,n_elements(a))
err[0,*]=a & err[1,*]=b & err[2,*]=c & err[3,*]=d & err[4,*]=e & err[5,*]=f & err[6,*]=g & err[7,*]=h & err[8,*]=i & err[9,*]=j
err[10,*]=k & err[11,*]=l
;------------------------------------------------------------------------------------
i=[-0.5,0,0.5,0,-0.5]
j=[0,-0.5,0,0.5,0]
usersym,i,j,color=60,/fill
;--------------------------------------------------------------------------------------
s0 = data[*,15] & s0_err = err[*,15]
s1 = data[*,11] & s1_err = err[*,11]
s2 = data[*,8] & s2_err = err[*,8]
s3 = data[*,6] & s3_err = err[*,6]
s4 = data[*,4] & s4_err = err[*,4]
s5 = data[*,2] & s5_err = err[*,2]
s6 = data[*,1] & s6_err = err[*,1]
s7 = data[*,0] & s7_err = err[*,0]
;---------------------------------------------------------------------------------------
HI = data[*,30] & HI_err = err[*,30]
irac1 = data[*,18] & irac1_err = err[*,18]
irac2 = data[*,19] & irac2_err = err[*,19]
irac3 = data[*,20] & irac3_err = err[*,20]
irac4 = data[*,21] & irac4_err = err[*,21]
mips24 = data[*,22] & mips24_err = err[*,22]  
mips70 = data[*,23] & mips70_err = err[*,23]
mips160 = data[*,24] & mips160_err = err[*,24]
pacs100 = data[*,25] & pacs100_err = err[*,25]
spire250 = data[*,27] & spire250_err = err[*,27]
spire500 = data[*,29] & spire500_err = err[*,29]
Neii = data[*,9] & Neii_err = err[*,9]
Neiii = data[*,10] & Neiii_err = err[*,10]
;---------------------------------------------------------------------------------------
;You can change the species you want to plot.
a = mips24 & b = mips70 & c = mips160 & d = spire250
aa = mips24_err & bb = mips70_err & cc = mips160_err & dd = spire250_err
;---------------------------------------------------------------------------------------
set_plot,'PS'
!p.font = 0
device, filename = 'color-color.ps', /helvetica, /portrait, /encapsulated, isolatin = 1, font_size = 16, decomposed = 0, /color
!p.thick = 4 & !x.thick = 3 & !y.thick = 3
plot,a/b,c/d,psym=8,xtitle='F!d24!n/F!d70!n',ytitle='F!d160!n/F!d250!n',xstyle=2,ystyle=2
oploterror,a/b,c/d,a/b*((aa/a)^2+(bb/b)^2)^0.5,c/d*((cc/c)^2+(dd/d)^2)^0.5,psym=8
device, /close_file, decomposed = 1
!p.multi = 0
;---------------------------------------------------------------------------------------
;Plot of all line emission verse HI gas
set_plot,'PS'
!p.font = 0
device, filename='All_line_strength_verse_HI.ps', /helvetica, /portrait, /encapsulated, isolatin = 1, font_size = 14, decomposed = 0, /color
!p.multi=[0,4,2]
!p.thick = 4 & !x.thick = 3 & !y.thick = 3
plot, HI,s0*1e10,psym=8,title = 'H!d2!n S(0)',xtitle='K km/s',ytitle='Flux density (10 !u-10!nJy)'
oploterror,HI,s0*1e10,HI_err,s0_err*1e10,psym=8
plot, HI,s1*1e10,psym=8,title = 'H!d2!n S(1)',xtitle='K km/s',ytitle='Flux density (10 !u-10!nJy)'
oploterror,HI,s1*1e10,HI_err,s1_err*1e10,psym=8
plot, HI,s2*1e10,psym=8,title = 'H!d2!n S(2)',xtitle='K km/s',ytitle='Flux density (10 !u-10!nJy)'
oploterror,HI,s2*1e10,HI_err,s2_err*1e10,psym=8
plot, HI,s3*1e10,psym=8,title = 'H!d2!n S(3)',xtitle='K km/s',ytitle='Flux density (10 !u-10!nJy)'
oploterror,HI,s3*1e10,HI_err,s3_err*1e10,psym=8
plot, HI,s4*1e10,psym=8,title = 'H!d2!n S(4)',xtitle='K km/s',ytitle='Flux density (10 !u-10!nJy)'
oploterror,HI,s4*1e10,HI_err,s4_err*1e10,psym=8
plot, HI,s5*1e10,psym=8,title = 'H!d2!n S(5)',xtitle='K km/s',ytitle='Flux density (10 !u-10!nJy)'
oploterror,HI,s5*1e10,HI_err,s5_err*1e10,psym=8
plot, HI,s6*1e10,psym=8,title = 'H!d2!n S(6)',xtitle='K km/s',ytitle='Flux density (10 !u-10!nJy)'
oploterror,HI,s6*1e10,HI_err,s6_err*1e10,psym=8
plot, HI,s7*1e10,psym=8,title = 'H!d2!n S(7)',xtitle='K km/s',ytitle='Flux density (10 !u-10!nJy)'
oploterror,HI,s7*1e10,HI_err,s7_err*1e10,psym=8
device, /close_file, decomposed = 1
!p.multi = 0
end
