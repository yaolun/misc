pro tsc, modeltime, c_s, Omega_0, r, theta, rhoenv, rhofloor, renv_in, indir=indir, outdir=outdir
;
; Set up the envelope
; Modified by MMD - This version calculates the TSC profile
;                   Read in the dd data from Susan's data file.
;                   Calculate x, tau, y for each point in the
;                   model (r,theta) grid.
;                   Put y into a 1-D array.
;                   Calculate nearest point in dd array for each value
;                   of y.
;                   Calculate the alphas using the IDL code I
;                   translated from Susan's fortran code.
;                   Also calculate the V's
;                   Put alphas and V's into 2D (r,theta) array.
;                   Write out velocity if last point is infalling.
;                   Calculate rho(r,theta) from the Cassen 
;                   & Moosman solution if x<tau^2                  
;                   Calculate rho(r,theta) from the equation in
;                   TSC84 Figure 3 caption, if x>=tau^2
;
;                   Set rho=rhofloor for all radii less than renv_in
;                   Set rho=rhofloor for all r>rout_new if t>time_past
;                   Calculate and print out envelope mass
;
; Modified by YLY
; INPUTS
; indir           - The input directory of the grid data from Susan. If not 
;                   assigned, the code will take its directory instead.
; outdir          - Output directory path
; modeltime       - The age of the protostellar system in years. 
; c_s             - Sound speed
; Omega_0         - The rotating speed as the perturbation to the collapse model
; r               - R-coordinate of the model
; theta           - Theta-coordinate of the model
; timestep        - (Not sure about the function of this right now)
; rhofloor        - The lowest density allowed in the model. Primarily designed for preventing zero in calculation.
; renv_in         - The inner radius of the envelope in cm.
;
; OUTPUTS
; rhoenv          - The density profile of envelope returned from the code.
; 



; Modified by YLY
; input setup
if not keyword_set(indir) then indir = ''
GG = 6.67259e-8        ; Gravitational constant
nr = n_elements(r)
nt = n_elements(theta)

print,''
print,'Using the TSC84 density profile'
print,''
mcol=12
lcol=987
Delta_Q=-3.52e-4
readcol,indir+'grid.plt34.mod',dd1,dd2,dd3,dd4,dd5,dd6,dd7,dd8,dd9,dd10,dd11,dd12,format='D,D,D,D,D,D,D,D,D,D,D,D',skipline=1
dd=[dd1,dd2,dd3,dd4,dd5,dd6,dd7,dd8,dd9,dd10,dd11,dd12]
sz_dd=n_elements(dd)

x_tsc=dblarr(nr,nt)
tau_tsc=dblarr(nr,nt)
y_tsc=dblarr(nr,nt)
y_tsc_1d=dblarr(nr*nt)
step=0
for i=0,nr-1 do begin
    for j=0,nt-1 do begin
        x_tsc[i,j]=r[i]/(c_s*modeltime*365.0*24.0*3600.0)
        tau_tsc[i,j]=Omega_0*(modeltime*365.0*24.0*3600.0)
        p2=(1./2.)*((3.0*cos(theta[j])*cos(theta[j]))-1.0)
        y_tsc[i,j]=x_tsc[i,j]*(1.+((tau_tsc[i,j]^2.0)*Delta_Q*p2))
        y_tsc_1d[step]=x_tsc[i,j]*(1.+((tau_tsc[i,j]^2.0)*Delta_Q*p2))
        step=step+1
    endfor
endfor
sz_y_tsc_1d=n_elements(y_tsc_1d)

lp=intarr(sz_y_tsc_1d)
for i=0,sz_y_tsc_1d-1 do begin
    for j=0,sz_dd-2 do begin
        if(j gt (lcol-1)) then begin
            lp[i]=999.0
            break
        endif
        if((y_tsc_1d[i]-dd[j])*(dd[j+1]-y_tsc_1d[i]) ge 0.0) then begin
            lp[i]=j
            break
        endif
    endfor
endfor

alpha_0_tsc_1d=dblarr(sz_y_tsc_1d)
alpha_M_tsc_1d=dblarr(sz_y_tsc_1d)
alpha_Q_tsc_1d=dblarr(sz_y_tsc_1d)
V_0_tsc_1d=dblarr(sz_y_tsc_1d)
V_M_tsc_1d=dblarr(sz_y_tsc_1d)
V_Q_tsc_1d=dblarr(sz_y_tsc_1d)
for i=0,sz_y_tsc_1d-1 do begin
    if(lp[i] ne 999.0) then begin
        alpha_0_tsc_1d[i]=dd[(lcol*1)+lp[i]]
        alpha_M_tsc_1d[i]=dd[(lcol*2)+lp[i]]
        alpha_Q_tsc_1d[i]=(-4./3.)*dd[(lcol*3)+lp[i]]
        V_0_tsc_1d[i]=dd[(lcol*4)+lp[i]]
        V_M_tsc_1d[i]=dd[(lcol*5)+lp[i]]
        V_Q_tsc_1d[i]=(-2./3.)*dd[(lcol*6)+lp[i]]
    endif
    if(lp[i] eq 999.0) then begin
        winterp=where(lp lt lp[i-1])
        temp=max(winterp)
        winterp=winterp[temp]
        alpha_0_tsc_1d[i]=2.0*(y_tsc_1d[i]^(-2.0))
        alpha_M_tsc_1d[i]=loglin_interp2pt(y_tsc_1d[winterp],alpha_M_tsc_1d[winterp],y_tsc_1d[i-1],alpha_M_tsc_1d[i-1],y_tsc_1d[i])
        alpha_Q_tsc_1d[i]=-1.0*loglin_interp2pt(y_tsc_1d[winterp],-1.0*alpha_Q_tsc_1d[winterp],y_tsc_1d[i-1],-1.0*alpha_Q_tsc_1d[i-1],y_tsc_1d[i])
        V_0_tsc_1d[i]=loglin_interp2pt(y_tsc_1d[winterp],V_0_tsc_1d[winterp],y_tsc_1d[i-1],V_0_tsc_1d[i-1],y_tsc_1d[i])
        V_M_tsc_1d[i]=loglin_interp2pt(y_tsc_1d[winterp],V_M_tsc_1d[winterp],y_tsc_1d[i-1],V_M_tsc_1d[i-1],y_tsc_1d[i])
        if(V_Q_tsc_1d[i-1] lt 0) then begin
            V_Q_tsc_1d[i]=-1.0*loglin_interp2pt(y_tsc_1d[winterp],-1.0*V_Q_tsc_1d[winterp],y_tsc_1d[i-1],-1.0*V_Q_tsc_1d[i-1],y_tsc_1d[i])
        endif
        if(V_Q_tsc_1d[i-1] gt 0) then begin
            V_Q_tsc_1d[i]=loglin_interp2pt(y_tsc_1d[winterp],V_Q_tsc_1d[winterp],y_tsc_1d[i-1],V_Q_tsc_1d[i-1],y_tsc_1d[i])
        endif
    endif
endfor

nan='nonan'
alpha_0_tsc=dblarr(nr,nt)
alpha_M_tsc=dblarr(nr,nt)
alpha_Q_tsc=dblarr(nr,nt)
rhoenv=dblarr(nr,nt)
V_0_tsc=dblarr(nr,nt)
V_M_tsc=dblarr(nr,nt)
V_Q_tsc=dblarr(nr,nt)
u_r=dblarr(nr,nt)
step=0
common share1, _a,_b
for i=0,nr-1 do begin
    for j=0,nt-1 do begin
        alpha_0_tsc[i,j]=alpha_0_tsc_1d[step]
        alpha_M_tsc[i,j]=alpha_M_tsc_1d[step]
        alpha_Q_tsc[i,j]=alpha_Q_tsc_1d[step]
        p2=(1./2.)*((3.0*cos(theta[j])*cos(theta[j]))-1.0)
        V_0_tsc[i,j]=V_0_tsc_1d[step]
        V_M_tsc[i,j]=V_M_tsc_1d[step]
        V_Q_tsc[i,j]=V_Q_tsc_1d[step]

        ;calculate radial velocity at each grid point
        ;if alog10(y)>=0 at (largest radius,pole)
        ;then write out this velocity to a data file
        ;if not, write out 0 to a data file
        u_r[i,j]=c_s*(V_0_tsc[i,j]+((tau_tsc[i,j]^2.0)*(V_M_tsc[i,j]+(V_Q_tsc[i,j]*p2))))
        if(alog10(y_tsc[nr-1,nt-1]) gt 0.0) then begin
            openw,lun,'velocity.dat',/get_lun, /append
            ; openw,2,'RESULTS/velocity.dat.'+strcompress(string(long(timestep)),/remove_all)
            printf,lun,0.0
            ; printf,2,0.0
            free_lun,lun
            close,lun
            ; close,2
        endif
        if(alog10(y_tsc[nr-1,nt-1]) le 0.0) then begin
            openw,lun,'velocity.dat',/get_lun, /append
            ; openw,2,'RESULTS/velocity.dat.'+strcompress(string(long(timestep)),/remove_all)
            printf,lun,(-1.0)*u_r[nr-1,nt-1]
            ; printf,2,(-1.0)*u_r[nr-1,nt-1]
            free_lun,lun
            close,lun
            ; close,2
        endif

        if(x_tsc[i,j] lt (tau_tsc[i,j]^2.0)) then begin    ;inner solution
            xmo=0.973546863
            if(x_tsc[i,j] eq 0.0) then rhoenv[i,j]=0.0
            if(x_tsc[i,j] gt 0.0) then begin
                fv=(-1.0)*sqrt(xmo/x_tsc[i,j])
                tsq=tau_tsc[i,j]*tau_tsc[i,j]           
                zeta=(xmo^3)*tsq/(16.*x_tsc[i,j])       ; TSC84, eq. 83 and 85 plus some information from Shu77
                ao=xmo/((x_tsc[i,j])^2.0)

                ;solve transcendental equation for theta_0
                ;first take care of special cases
                ;then, if theta_0 not yet assigned, call cubesolve
                theta_0=-999.0
                pi2=!pi/2.0
                if(zeta eq 0.0) then theta_0=theta[j]
                if(theta[j] eq 0.0) then theta_0=0.0
                if(abs((theta[j]/!pi)-1) le 1.0d-5) then theta_0=!pi
                if(abs((theta[j]/pi2)-1) le 1.0d-5) then begin
                    if(zeta le 1.0) then begin
                        theta_0=pi2
                    endif else begin
                        theta_0=-100.0
                    endelse
                endif
                if(theta_0 eq -999.0) then begin ;continue on if needed
                    temp1=zeta
                    temp2=0.0
                    temp3=(1.0-zeta)
                    temp4=(-1.0)*cos(theta[j])
                    tsc_b=[temp1,temp2,temp3,temp4]
                    xone=1.0
                    ; result_cubesolve=cubesolve(tsc_b)
                    ; if(result_cubesolve[3] gt 0.0) then begin
                    ;     if((result_cubesolve[0] ge ((-1.0)*xone)) and (result_cubesolve[0] le xone)) then begin
                    ;         costheta_0=result_cubesolve[0]
                    ;     endif else begin
                    ;         costheta_0=-0.5
                    ;         print,'Problem with cube root, costheta_0=',result_cubesolve[0]
                    ;     endelse
                    ; endif else begin
                    ;     costheta_0=-0.5
                    ;     iroot=0
                    ;     for z=0,2 do begin
                    ;         if(((cos(theta[j])*result_cubesolve[z]) ge 0.0) and (result_cubesolve[z] ge ((-1.0)*xone)) and (result_cubesolve[z] le xone)) then begin
                    ;             costheta_0=p[z]
                    ;             iroot=iroot+1
                    ;         endif
                    ;     endfor
                    ; endelse
                    ; Modified by YLY.  Implement a new cubic solver
                    result_cubesolve = cuberoot([temp4,temp3,temp2,temp1])
                    if n_elements(where(result_cubesolve le -1e20)) eq 2 then begin
                        if (result_cubesolve[0] ge -1d0) and (result_cubesolve[0] le 1d0) then begin
                            costheta_0=result_cubesolve[0]
                        endif else begin
                            costheta_0 = -0.5
                            print,'Problem with cube root, costheta_0=',result_cubesolve[0]
                            print, '[x^3, x^2, x^1, x^0]', [temp1,temp2,temp3,temp4]
                            stop
                        endelse
                    endif else begin
                        costheta_0 = -0.5
                        iroot = 0
                        for z = 0, 2 do begin
                            if(((cos(theta[j])*result_cubesolve[z]) ge 0.0) and (result_cubesolve[z] ge -1d0) and (result_cubesolve[z] le 1d0)) then begin
                                costheta_0=result_cubesolve[z]
                                iroot=iroot+1
                            endif
                        endfor
                    endelse

                    theta_0=acos(costheta_0)
                endif

                ;now i have theta_0
                ;special case 1:  in disk, on equator
                rhoenv[i,j]=-999.0
                if(theta_0 lt -1.0) then begin
                    rhoenv[i,j]=0.0
                endif
                ;special case 2: outside disk, on equator
                if(abs((theta[j]/pi2)-1.0) le 1.0e-5) then begin
                    f4=1.+(2.*zeta*((1.-(1.5*sin(theta_0)))^(2.0)))
                    rhoenv[i,j]=(1./(4.*!pi*GG*((modeltime*365.0*24.0*3600.0)^2.0)))*((-1.0)*ao)/(fv*sqrt(2.)*f4)
                endif
                ;general case if special cases not used
                if(rhoenv[i,j] eq -999.0) then begin
                    f1=sqrt(1.+(cos(theta[j])/cos(theta_0)))
                    f4=1.+(2.*zeta*((1.-(1.5*sin(theta_0)))^(2.0)))
                    rhoenv[i,j]=(1./(4.*!pi*GG*((modeltime*365.0*24.0*3600.0)^2.0)))*((-1.0)*ao)/(fv*f1*f4)
                endif
            endif
            if finite(rhoenv[i,j]) eq 0 then stop
            ; if (i eq 0) and (j eq 67) then stop 
            ;rhoenv[i,j]=1.0e-4*(muh2*1.67e-24) ;switch off inner solution
        endif
        if(x_tsc[i,j] ge (tau_tsc[i,j]^2.0)) then begin    ;TSC solution
            rhoenv[i,j]=(1./(4.*!pi*GG*((modeltime*365.0*24.0*3600.0)^2.0)))*(alpha_0_tsc[i,j]+((tau_tsc[i,j]^2.0)*(alpha_M_tsc[i,j]+(alpha_Q_tsc[i,j]*p2))))
        endif
        if(r[i] lt renv_in) then rhoenv[i,j] = rhofloor
        step=step+1

        ; if (rhoend[i,j] lt rhofloor) then rhoenv[i,j] = rhofloor
        ; if(finite(rhoenv[i,j]) eq 0) then begin
        ;    rhoenv[i,j]=0.0
        ;    nan='nan'
        ; endif

        ;rhoenv[i,j]=0.0 ;switch off the envelope

    endfor
endfor

; mass_envelope_temp=mass(rhoenv/gastodust,r,theta,gastodust=gastodust)
; ratio=menv_input/mass_envelope_temp
; rhoenv=rhoenv*ratio
; mass_envelope_final=mass(rhoenv/gastodust,r,theta,gastodust=gastodust)
; print,''
; print,'The envelope mass should be:  ',menv_input/MS,' Msun'
; print,'The envelope mass was originally:  ',mass_envelope_temp/MS,'  Msun'
; print,'The ratio of the two is:  ',ratio
; print,'The modified, final envelope mass is:  ',mass_envelope_final/MS,' Msun'
; print,''
; openw,1,'RESULTS_MENV/menv_'+strcompress(string(long(timestep)),/remove_all)+'.dat'
; printf,1,modeltime,mass_envelope_final/MS,'    ',nan
; close,1

;
; Mutually exclude the densities
; 
; ii = where( rhodisk gt rhoenv ) 
; if ii[0] ge 0 then rhoenv[ii] = rhofloor
; ii = where( rhodisk le rhoenv ) 
; if ii[0] ge 0 then rhodisk[ii] = rhofloor

; rho=rhodisk+rhoenv

end

