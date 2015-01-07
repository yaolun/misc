function cubesolve,ab
;adapted from TSC ncofrac.f
;returns 1d array with 4 elements: first three are the roots
;               the 4th gives info about the roots
;               mtype= 1    one real, two complex conjugate roots
;               mtype= 0    3 real roots
;               mtype=-1    3 real distinct roots

x=dblarr(4)
x[0]=-999.0
x[1]=-999.0
x[2]=-999.0
x[3]=-999.0
pi=!pi

if(abs(ab[0]) lt 1.0e-10) then begin
    print,'Coefficient of cubic term in CUBESOLVE < 1.0e-10'
    print,'Aborting!!!'
    stop
endif

b=ab

;preliminary definitions
p=b[1]/b[0]
q=b[2]/b[0]
r=b[3]/b[0]
ba=q-(p*p/3.)
bb=((2.*(p^3))-(9.*p*q)+(27.*r))/27.0
det=((bb*bb)/4.)+((ba^3)/27.)

;find solution, check sign of 'det' to find root types

;one real, two complex
if(det gt 0.0) then begin
    ca=((-bb/2.)+sqrt(det))
    if(ca ge 0.0) then sign=1.0
    if(ca lt 0.0) then sign=-1.0
    capa=1.0*sign*(abs(ca)^(1./3.))
    cb=((-bb/2.)-sqrt(det))
    if(cb ge 0.0) then sign=1.0
    if(cb lt 0.0) then sign=-1.0
    capb=1.0*sign*(abs(cb)^(1./3.))
    root1=capa+capb
    x[0]=root1-(p/3.)
    x[3]=1.0
endif

;three real distinct roots
if(det lt 0.0) then begin
    rg=(-bb/2.)/sqrt(-ba^(3./27.))
    phi=acos(rg)
    cons=2.0*sqrt(-ba/3.)
    x[0]=cons*cos(phi/3.)-(p/3.)
    x[1]=cons*cos((phi/3.)+(2.*pi/3.))-(p/3.)
    x[2]=cons*cos((phi/3.)+(4.*pi/3.))-(p/3.)
    x[3]=-1.0
endif

;three real roots
if(det eq 0.0) then begin
    bbtemp=-1.0*bb
    if(bbtemp ge 0.0) then sign=1.0
    if(bbtemp lt 0.0) then sign=-1.0
    capa=1.0*sign*(abs(-bb/2.)^(1./3.))
    x[0]=2.*capa-(p/3.)
    x[1]=(-1.0*capa)-(p/3.)
    x[3]=x[2]
    x[4]=0.0
endif

return,x

end
