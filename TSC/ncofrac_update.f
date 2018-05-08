        program ncofrac
!---------------------------------------------------
! calculates spectral line profile for collapsing cloud
! input: ND-time tau, inclination of rotation axis w.r.t. l.o.s.
!       range in ND-distance x, number of theta,phi grid points
!       output file name, input data file 'grid.dta:?'
! output: ND-l.o.s. velocity vz, ND-mass element, position,velocity
!----------------------------------------------------
        parameter(nangle=10,npoints=1000)
        character*20 infile,outfile,txtfile,oldfile,velfile
        dimension aid(nangle),sum(nangle),vlos(npoints)
        logical error,test
        common /climits/ xrmin,xrmax,delxr,ntmax,npmax
        common /collpar/ d2,time,pi,pi2
        common /slimits/ nvmax,vmin,vmax,delv,ninc,ai(nangle)
        common /sprofile/ data(npoints,nangle)
!
        print *,' program ncofrac:'
        print *,' calculates spectral line profile data for collapsing cloud'
! input file names
10      continue
        infile = 'grid.plt34.mod'
        call oneddata(infile,error)
        if(error)goto10
! output file names
        txtfile = 'output_summary'
        open(unit=7, file=txtfile, status='replace')
        outfile = 'output_data'
        velfile = 'rho_v_env'
        test = .FALSE.

C         read the parameters from a file
          open(unit=13, file='tsc.par')
          read(13,*) vmin, vmax, delv, ninc, (aid(ni),ni=1,ninc), tau
          read(13,*) xrmin, xrmax, delxr, ntmax, npmax

          print *, vmin, vmax, delv, ninc, (aid(ni),ni=1,ninc), tau
          print *, xrmin, xrmax, delxr, ntmax, npmax
          nvmax= (vmax-vmin)/delv +1.0

! calculate solution
        do 40 ni=1,ninc
          ai(ni)=aid(ni)*pi/180.
40      continue

!       YLY: for debug
C        open(unit=11,file='ncollapse_debug',status='new')
!       YLY: for writing out density and velocities for given (r, theta)
        open(unit=4,file=velfile,status='replace')
        write(4,*) 'lp xr theta ro ur utheta uphi'
!
        call ncollapse(tau)
!
        open(unit=9,file=outfile,status='replace')
        write(9,802) nvmax,vmin,vmax,delv,tau,ninc,(aid(ni),ni=1,ninc)
        do 80 nv=1,nvmax
            v= vmin +float(nv-1)*delv
            vlos(nv)=v
            do 70 ni=1,ninc
              data(nv,ni)=data(nv,ni)/delv
70          continue
            write(9,*) v,(data(nv,ni),ni=1,ninc)
80      continue
          close(9)
            do 90 ni=1,ninc
              sum(ni)=avint(vlos,data(1,ni),nvmax,vlos(1),vlos(nvmax),ind)
              if(ind.eq.0) write(7,*)' avint problems, ni=',ni
90          continue
        write(7,*)' sum over profiles for each inclination'
        write(7,*) (sum(ni),ni=1,ninc)
        print *,' program nco sucessfully completed'
!       write(0,'("program nco sucessfully completed")')
        write(7,*)' program nco sucessfully completed'
        close(7)
!       YLY: close files
        close(4)
!
999     continue
        stop
802     format(1xi6,2f10.1,f6.1,f6.3,i3,16f5.1)
        end
        subroutine oneddata(datafile,error)
!----------------------------------------------------
!  read 1d data file 'datafile', store in array 'dd'
!----------------------------------------------------
        character*(*) datafile
        logical error
        parameter (lmsize=500000)
        common /pldata/ lcol,mcol,dd(lmsize)
!
        lun=8
        open(lun,file=datafile,status='old',err=999)
!
        read(lun,*) mcol,lcol
        print *,mcol,lcol
        do 10 l=1,lcol
          print *,l,m
          read(lun,*,end=900) (dd((m-1)*lcol+l),m=1,mcol)
!         add by YLY
          error=.false.
!         print *,dd((m-1)*lcol+l),m,mcol
10      continue
!
        close(lun)
!
        return
!
900     continue
        write(6,*) ' file smaller than assumed bounds'
        error=.true.
        return
999     continue
        write(6,*) ' error in opening file',datafile
        error=.true.
        return
        end
        subroutine ncollapse(tau)
!----------------------------------------------------------
! modified from thesis subroutine COLLAPSE
! calculates non-dim-density and velocity of collapse model as a function
! of position at a given non-dimensional time tau
!----------------------------------------------------------
        parameter(nangle=10,npoints=1000)
        common /pldata/ lc,mc,dd(1)
        common /climits/ xrmin,xrmax,delxr,ntmax,npmax
        common /collpar/ d2,time,pi,pi2
        common /slimits/ nvmax,vmin,vmax,delv,ninc,ai(nangle)
        common /sprofile/ data(npoints,nangle)
        dimension vz(20),vzold(20)
        data  pi2/1.5707963/,pi/3.1415927/
        data d2/-3.52E-4/

!
! YLY: Legendre polynomial - P2(cos(theta)) = 1-3/2*sin(theta)
        p2(ang)= 1. -1.5*(sin(ang)**2)
!
        time=tau
        tsq= tau*tau
! set up
        nrmax=(xrmax-xrmin)/delxr +1.0
        delr=delxr
        thetamin=0.
        thetamax=pi
        phimin=0.
        phimax=2.*pi
        delth=(thetamax-thetamin)/float(ntmax-1)
! set so program won't do both zero and 2pi in angle phi
        delph=(phimax-phimin)/float(npmax)
!
! initialize lp
! YLY: don't know what "lp" stands for
        lp=0
! ang1=theta, ang2=phi
!  step through volume, increment in xr,theta,phi
! if xr .eq. zero, then skip iteration
        nrmin=1
        if(xrmin.eq.0.) nrmin=2
        nskip=0
        do 600 nr=nrmin,nrmax
          xr= xrmin +float(nr-1)*delxr
! related to volume element
! YLY: write out msg if running the last element in r-grid
          if(nr.eq.nrmax) then
            delr=abs(xrmax-xr)
            write(7,*)' last step nr,xr,xrmax,delr',nr,xr,xrmax,delr
          end if
! YLY: looks like a volume integral element
          const= xr**2 *delth*delph*delr
! loop over angle theta
          do 500 nt= 1,ntmax
            theta=thetamin +float(nt-1)*delth
! YLY: theta_i in theta-grid
            ang1=theta
! YLY: dV (volume)
            dvol=abs(sin(theta)*const)
! if theta= pi set dvol=0 to avoid roundoff
            if(nt.eq.ntmax) dvol=0.
!
!  calculate nearest point in dd array
! YLY: Not sure the "p2" is properly defined
            yr=xr*(1. +d2*p2(ang1)*tau**2)
C            write(11,*), "xr, theta", xr, theta
C            write(11,*), "ang1, p2(ang1)", ang1, p2(ang1)
C            write(11,*), "yr", yr
! YLY: Where is "lc" defined?
! YLY: "lp" turns out always be 0.
            do 75 l=lp-1,lc-1
              if(l.le.0) then
                if(yr.lt.dd(1).or.yr.gt.dd(lc)) then
                        goto100
                else
                        goto50
                end if
              end if
              if( (yr-dd(l))*(dd(l+1)-yr).ge.0)then
                lp=l
                        goto100
              end if
50          continue
75          continue
            lp=0
!
100         continue
! YLY: Looks like the above do loop is to figure which value of lp
!      should be used
!
!  calculate density,velocity at position xr,theta
!       calculate inner solution
                if (xr.lt.tsq) then
                  call cassen(tau,ang1,xr,ro,ur,utheta,uphi)
                else
                  if (lp.eq.0) then
                      ro=0.
                      ur=0.
                      utheta=0.
                      uphi=0.
                  else
! YLY: This part looks like linear interpolation for density and velocities.
                      r2=dd(lp+1)
                      r1=dd(lp)
                      call outersoln(tau,ang1,lp,ro1,ur1,utheta1,uphi1)
                      call outersoln(tau,ang1,lp+1,ro2,ur2,utheta2,uphi2)
                      delta=(yr-r1)/(r2-r1)
! YLY: sound speed information seems missing
                      ro= ro1 +(ro2-ro1)*delta
                      ur= ur1 +(ur2-ur1)*delta
                      utheta= utheta1 +(utheta2-utheta1)*delta
                      uphi= uphi1 +(uphi2-uphi1)*delta
                  end if
                end if
! YLY: Try to output the density and velocity in spherical coordinates for
!      given (r, theta)
          write(4,*) lp, xr, theta, ro, ur, utheta, uphi
!
! calculate mass element, l.o.s. velocity for all inclination angles
          delm= ro*dvol
! initialize vzold for loop over angle phi
! YLY: not sure why need to do this
! YLY: In axisymmetric model, phi won't be matter for getting los velocity.
!      This line could be for preventing error if phimin = 0.
          phi= phimin -delph
          do 200 ni=1,ninc
           call vzpvect(ai(ni),theta,phi,ur,utheta,uphi,vzold(ni))
200       continue
! loop over angle phi
          do 400 np= 1,npmax
                phi= phimin +float(np-1)*delph
! YLY: Don't see this variable "ang2" is used in the following lines of
!      ncollapse.  Not sure it is for sharing outside of the subroutine.
                ang2=phi
! YLY: loop over inclination angles
          do 300 ni=1,ninc
! calculate l.o.s. velocity
            call vzpvect(ai(ni),theta,phi,ur,utheta,uphi,vz(ni))
! check if vz-vzold < delv
! YLY: The only difference between these two "vzpvect" is the input phi.
!      One is phimin-delph, and another one is looping over phi-grid.
!      But I don't know what's the reason behind this.
            if( abs((vz(ni)-vzold(ni))/delv) .gt.0.5) then
!              print *,' velocity channel skip',vz(ni),vzold(ni),xr,theta,phi
               write(7,*)' vel. channel skip',vz(ni),vzold(ni),xr,theta,phi
               nskip=nskip +1
! YLY: write out the reason for stop
!              if(nskip .gt. 30) stop
!
               if(nskip .gt. 30) then
                   write(7,*)'skipped channels exceeded max value', nskip
!                   stop
               end if
            end if
!           YLY: debug
            write(11,*), "nr, nt, np", nr, nt, np
            write(11,*), "vz(ni), vzold(ni)",vz(ni), vzold(ni)
            write(11,*), "vmin", vmin
            write(11,*), "delv", delv
            write(11,*), "nskip", nskip

!
! calculate velocity channel appropriate to velocity vz
            nv=nint( (vz(ni)-vmin)/delv +1.0 )
! calculate velocity at left hand side of channel
            vl= vmin +(float(nv-1)-0.5)*delv
! calculate velocity at right hand side of channel
            vu= vl+delv
! YLY: make sure the index of velocity not exceeding the size of velocity-grid
!      in the input file.  max size = 1000
            if (nv.ge. 1 .and. nv.le. npoints) then
! put mass element into velocity profile for inclination ai(ni)
             if( vl.le.vzold(ni).and.vzold(ni).le.vu ) then
!             write(7,*)' ni,nv,vz,vzold,v',ni,nv,vz(ni),vzold(ni),v
              temp= data(nv,ni)
              data(nv,ni)= data(nv,ni) +delm
!             write(7,*)' d(nv,ni)new,old',data(nv,ni),temp
             elseif( vzold(ni).lt.vl) then
              frac=abs((vz(ni)-vl)/(vz(ni)-vzold(ni)))
!             write(7,*)' ni,nv,vz,vzold,v,frac',ni,nv,vz(ni),vzold(ni),vl,frac
!             write(7,*)' d(nv,ni),d(nv-1,ni)old',data(nv,ni),data(nv-1,ni)
              data(nv,ni)=data(nv,ni) +delm*frac
              if(nv.gt.1) data(nv-1,ni)=data(nv-1,ni)+delm*(1.0-frac)
!             write(7,*)' d(nv,ni),d(nv-1,ni)new',data(nv,ni),data(nv-1,ni)
             else
              frac=abs((vz(ni)-vu)/(vz(ni)-vzold(ni)))
!             write(7,*)' ni,nv,vz,vzold,v,frac',ni,nv,vz(ni),vzold(ni),vu,frac
!             write(7,*)' d(nv,ni),d(nv+1,ni)old',data(nv,ni),data(nv+1,ni)
              data(nv,ni)=data(nv,ni)+delm*frac
              if(nv.lt.npoints) data(nv+1,ni)=data(nv+1,ni)+delm*(1.0-frac)
!             write(7,*)' d(nv,ni),d(nv+1,ni)new',data(nv,ni),data(nv+1,ni)
              end if
            else
              write(7,*)' vel. outside bounds',vz(ni)
!             print *,' vel. outside bounds',vz(ni)
            end if
! save velocity for checking velocity channel skip
            vzold(ni)=vz(ni)
!         write(7,800) vz(ni),delm,xr,theta,phi,ro,ur,utheta,uphi
300       continue
! if theta=0 then only do phi=0 iteration
          if(nt.eq.1 .and. thetamin.eq.0.) goto500
! if theta=2pi then only do phi=0 iteration
          if(nt.eq.ntmax) goto500
! end loop through volume
!         write(7,800) vz(ni),delm,xr,theta,phi,ro,ur,utheta,uphi,dvol
400         continue
500       continue
          write(7,800) vz(ni),delm,xr,theta,phi,ro,ur,utheta,uphi,dvol
600     continue
!
        return
800     format(1x9g12.4)
        end
        subroutine outersoln(tau,ang,lp,ro,ur,utheta,uphi)
!-------------------------------------------------------
! calculates physical variables density,velocities
! for collapse solution
!-------------------------------------------------------
        common /pldata/ lc,nc,dd(1)
!
        p2(a)= 1. -1.5*(sin(a))**2
        dp2(a)= -3.*sin(a)*cos(a)
!
        tsq=tau*tau
! density-----------------------------------------
10      continue
        alpha=dd(lc*1 +lp)
        alpham=dd(lc*2 +lp)
        alphaq= -4./3. *dd(lc*3+lp)
        ro= alpha +tsq*(alpham +p2(ang)*alphaq)
! Ur----------------------------------------------
20      continue
        v=dd(lc*4+lp)
        vm=dd(lc*5+lp)
        vq=-2./3.*dd(lc*6+lp)
        ur= v +tsq*( vm+ p2(ang)*vq)
! Utheta------------------------------------------
30      continue
        wq=-1./3.*dd(lc*7+lp)
        utheta=tsq*dp2(ang)*wq
! Uphi--------------------------------------------
40      continue
        f=dd(lc*8+lp)
        x=dd(lp)
        if(x.eq.0) then
                uphi=0.
        else
                uphi= tau*sin(ang)*f/x
        end if
        return
        entry outercart(tau,ang,lp,ro,ux,uy,uz)
!-------------------------------------------------
! calculate collapse solution, velocity in cartesian coordinates
! YLY: This entry is not used by the code here.
!-------------------------------------------------
        tsq=tau*tau
! density-----------------------------------------
45      continue
        alpha=dd(lc*1 +lp)
        alpham=dd(lc*2 +lp)
        alphaq= -4./3. *dd(lc*3+lp)
        ro= alpha +tsq*(alpham +p2(ang)*alphaq)
! Uw----------------------------------------------
!       ur
50      continue
        v=dd(lc*4+lp)
        vm=dd(lc*5+lp)
        vq=-2./3.*dd(lc*6+lp)
        ur= v +tsq*( vm+ p2(ang)*vq)
!       utheta
        wq=-1./3.*dd(lc*7+lp)
        utheta=tsq*dp2(ang)*wq
!
        uw= ur*sin(ang) +utheta*cos(ang)
! Uz----------------------------------------------
!
        uz= ur*cos(ang) -utheta*sin(ang)
! Ux----------------------------------------------
!   uphi
70      continue
        f=dd(lc*8+lp)
        x=dd(lp)
        if(x.eq.0) then
                value=0.
        else
                value= tau*sin(ang)*f/x
        end if
!
!
        ux= uw*cos(ang2) -uphi*sin(ang2)
! Uy----------------------------------------------
!
        uy= uw*sin(ang2) +uphi*cos(ang2)
!
        return
        end
        subroutine cassen(tau,theta,xr,ro,ur,utheta,uphi)
!--------------------------------------------------------------
! calculate Cassen and Moosman inner solution
!--------------------------------------------------------------
        common /collpar/ d2,time,pi,pi2
        data xmo/ .97546863/
        f1(a1,a2)= sqrt( 1.+cos(a1)/cos(a2))
        f2(a1,a2)= sqrt( 1.-cos(a1)/cos(a2))
        f3(a1,a2)= (cos(a1)-cos(a2))/sin(a1)
        f4(z1,a1)= 1. +2.*z1*(1.-1.5*(sin(a1))**2)
!
!       print *,' tau,theta,xr',tau,theta,xr
        save
        if(xr.eq.0.) then
                ur=0.
                utheta=0.
                uphi=0.
                ro=0.
                return
        end if
!
        fv= - sqrt(xmo/xr)
!
        tsq= tau*tau
        zeta= (xmo**3)*tsq/(16.*xr)
        ao= xmo/(xr)**2
!       print *,' zeta,ao,fv',zeta,ao,fv
!
        call method1(zeta,theta,thetao)
        if( thetao.lt.-1.) then
!               print *,' no thetao: x,zeta,theta=',xr,zeta,theta
! in disk, on equator
                ur=0.
                utheta=0.
                uphi=0.
                ro=0.
                return
        end if
        if( abs(theta/pi2-1.) .le. 1.e-5) then
! outside disk, on equator
                ur= fv*sqrt(2.)
                utheta=0.
                uphi=0.
                ro= -ao/(fv*sqrt(2.)*f4(zeta,thetao))
                return
        end if
!       print *,' theta,thetao',theta,thetao
!
! density-----------------------------------------
10      continue
        ro= -ao/(fv*f1(theta,thetao)*f4(zeta,thetao))
!       print *,' ro',ro
! Ur----------------------------------------------
20      continue
        ur= fv*f1(theta,thetao)
!       print *,' ur',ur
! Utheta------------------------------------------
30      continue
        if(theta.eq.0. .or. abs(theta/pi-1.).le.1e-5)then
         utheta=0.
        else
         utheta= f3(theta,thetao)*fv*f1(theta,thetao)
        end if
!       print *,' utheta',utheta
! Uphi--------------------------------------------
40      continue
        if(theta.eq.0. .or. abs(theta/pi-1.).le.1e-5)then
         uphi=0.
        else
         uphi= -fv*sin(thetao)/sin(theta) * f2(theta,thetao)
        end if
!       print *,' uphi',uphi
! Uw----------------------------------------------
!       ur
50      continue
!       ur= fv*f1(theta,thetao)
!       utheta
!       utheta= f3(theta,thetao)*ur
!
!       uw= ur*sin(theta) +utheta*cos(theta)
!       v=uw
!       return
! Uz----------------------------------------------
!       ur
60      continue
!       ur= fv*f1(theta,thetao)
!       utheta
!       utheta= f3(theta,thetao)*ur
!
!       uz= ur*cos(theta) -utheta*sin(theta)
!       v= uz
!
        return
        end
        subroutine method1(azeta,atheta,athetao)
!----------------------------------------------------------
!  solve transcendental equation
!  derived by Cassen and Moosman for parabolic trajectories
!  zeta * sin**2(thetao) = 1 - cos(theta)/cos(thetao)
!  YLY: Eq. 84 in TSC84
!  given zeta, theta solve for thetao
!----------------------------------------------------------
        implicit real*8 (b-h,o-z)
        dimension b(4),x(3)
        common /cparam/ pi
!
        if (azeta.eq.0.) then
                athetao= atheta
!               type 10, azeta,atheta,athetao
10      format(' special case: zeta=',f8.2,' theta=',f8.2,' thetao=',f8.2)
         return
        end if
        if (atheta.eq.0.) then
                athetao= 0.
!               type 10, azeta,atheta,athetao
                return
        end if
        if ( abs(atheta/pi-1.).le. 1.d-5 ) then
                athetao= pi
!               type 10, azeta,atheta,athetao
                return
        end if
        pi2=pi/2.0d0
        if ( abs(atheta/pi2-1.).le. 1.d-5 ) then
                if(azeta.le.1.) then
                 athetao= pi2
                else
                 athetao= -100.
                end if
!               type 10, azeta,atheta,athetao
                return
        end if
! single precision to double
        zeta= azeta
        theta= atheta
!
        costheta = dcos(theta)
        b(1)= zeta
        b(2)= 0.
        b(3)= 1.- zeta
        b(4)= - costheta
        xone= 1.
!
! find roots of cubic equation, put result into array x
        call cubesolv(b,x,mtype)
!
        if (mtype.gt.0) then
                if (x(1).ge.-xone .and. x(1).le.xone) then
                        costhetao= x(1)
                else
                        costhetao= -.5
!                       print *,' problem with cube root, costhetao=',x(1)
                        write(7,*)' problem with cube root, costhetao=',x(1)
                end if
!               type *,'  x=',x(1),' re(xc)=',x(2),' im(xc)=',x(3)
!               print *,'  x=',x(1)
        else
                costhetao= -.5
                iroot=0
                do 50 i=1,3
                  if(costheta*x(i) .ge. 0.
     &                  .and. x(i).ge.-xone .and. x(i).le.xone) then
                                costhetao= x(i)
                                iroot=iroot+1
!                               type *, ' ***x=', x(i)
                        end if
!                       print *,'  x=',x(i)
50              continue
                if(iroot.eq.0 .or. iroot .gt. 1) then
!                 print *,' problem with cube root, costhetao=,iroot=',x,iroot
                  write(7,*)' problem with cube root, costhetao,iroot=',x,iroot
                end if
        end if
!
55      continue
        thetao= dacos(costhetao)
!       write (*,5)mtype,zeta,theta,thetao
5       format(' mt=',i5,' zeta=',f8.2,' theta=',f8.2,' thetao=',f8.2)
!
! double precison to single precision
        athetao= thetao
        return
        end
        subroutine vzpvect(ai,theta,phi,vr,vtheta,vphi,vzp)
!-------------------------------------------
! given spherical velocity vector vr,vtheta,vph
! spherical coordinates theta,phi and inclination ai (in radians)
! then find the line of sight velocity component vzp
!-------------------------------------------
        sint=sin(theta)
        cost=cos(theta)
        sinp=sin(phi)
        cosp=cos(phi)
        sini=sin(ai)
        cosi=cos(ai)
!       print *,' cos(theta,phi,ai),sin( )',cost,cosp,cosi,sint,sinp,sini
!
        terma=vr*(sint*sinp*sini +cost*cosi)
        termb=vtheta*(cost*sinp*sini -sint*cosi)
        termc=vphi*(cosp*sini)
        vzp=terma +termb +termc
!       print *,' vz,terma,b,c',vz,terma,termb,termc
!
        return
        end
        subroutine cubesolv(ab,ax,mtype)
!--------------------------------------------------
! solves for roots of cubic equation of form
!   b(1)x**3 +b(2)x**2 +b(3)x +b(4) =0
! mtype = 1   one real, two complex conjugate roots
!       = 0   3 real roots
!       =-1   3 real distinct roots
! solution from CRC Handbook p. A-245
!---------------------------------------------------
        implicit real*8(a-h,o-z)
        dimension b(4),x(3)
        dimension ab(4),ax(3)
        common /cparam/ pi
        data pi/3.141592653589793238/
!
        if( abs(ab(1)).lt. 1e-10) then
          print *, ' coefficient of cubic term in CUBESOLV .lt. 1e-10'
!         write(0,'("coefficient of cubic term in CUBESOLV .lt. 1e-10")')
          stop
        end if
! real to double precision
        do 10 i=1,4
          b(i)=ab(i)
10      continue
! preliminary definitions
        p= b(2)/b(1)
        q= b(3)/b(1)
        r= b(4)/b(1)
!       print *,' p,q,r',p,q,r
        ba= q - p*p/3.
        bb= (2.*p**3 -9.*p*q +27.*r)/27.
        det= bb*bb/4. + (ba**3)/27.
!       print *,' ba,bb,det',ba,bb,det
! find solution, check sign of 'det' to find root types
! one real, two complex
        if (det.gt.0.) then
          ca=(-bb/2. +dsqrt(det))
          capa=dsign(1.0_8,ca)*(abs(ca))**(1.0_8/3.0_8)
          cb=(-bb/2. -dsqrt(det))
          capb=dsign(1.0_8,cb)*(abs(cb))**(1.0_8/3.0_8)
          root1= capa + capb
!         print *,' capa,capb,root1',capa,capb,root1
          x(1)=root1- p/3.
          mtype=1
          goto20
        else if (det.lt.0.) then
! three real distinct roots
          rg= (-bb/2.)/sqrt(-ba**3/27.)
          phi= dacos(rg)
          cons= 2.*dsqrt(-ba/3.)
!         print *,' rg,phi,cons',rg,phi,cons
          x(1)= cons*dcos(phi/3.) -p/3.
          x(2)= cons*dcos(phi/3. + 2.*pi/3.) -p/3.
          x(3)= cons*dcos(phi/3. + 4.*pi/3.) -p/3.
          mtype= -1
          goto20
        else
! three real roots
          capa=dsign(1.0_8,-bb)*(abs(-bb/2.0_8))**(1.0_8/3.0_8)
!         print *,' capa',capa
          x(1)=2.*capa -p/3.
          x(2)=-capa -p/3.
          x(3)= x(2)
          mtype=0
          goto20
        end if
!
20      continue
! complex to real conversion
        do 30 i=1,3
         ax(i)=x(i)
30      continue
        return
        end
        function avint(x,y,n,xlo,xup,ind)
!
!-------------------------------------------------
!      integrates y on nonuniform grid
!      input independent value array x and dependent
!      value array y of size n
!      xlo and xup are desired bounds of integration
!      ind is an error flag, ind =0 for error
!      x in ascending order; n .ge. 3
!      xlo,xup in interval [x(1),x(n)]
!
!--------------------------------------------------
!
!       implicit real*8 (a-h,o-z)
        dimension x(n),y(n)
        ind=0
        if(n.lt.3) return
        do 10 i=2,n
        if(x(i).le.x(i-1)) return
10      continue
        sum=0
        if(xlo.le.xup) goto5
        syl=xup
        xup=xlo
        xlo=syl
        ind=-1
        goto6
5       ind=1
        syl=xlo
6       ib=1
        j=n
        do 1 i=1,n
        if(x(i).ge.xlo) goto7
1       ib=ib+1
7       ib=max(2,ib)
        ib=min(ib,n-1)
        do 2 i=1,n
        if(xup.ge.x(j)) goto8
2       j=j-1
8       j=min(j,n-1)
        j=max(ib,j-1)
        do 3 jm=ib,j
        x1=x(jm-1)
        x2=x(jm)
        x3=x(jm+1)
        term1=y(jm-1)/((x1-x2)*(x1-x3))
        term2=y(jm)/((x2-x1)*(x2-x3))
        term3=y(jm+1)/((x3-x1)*(x3-x2))
        a=term1+term2+term3
        b=-(x2+x3)*term1-(x1+x3)*term2-(x1+x2)*term3
        c=x2*x3*term1+x1*x3*term2+x1*x2*term3
        if(jm.gt.ib) goto14
        ca=a
        cb=b
        cc=c
        goto15
14      ca=.5*(a+ca)
        cb=.5*(b+cb)
        cc=.5*(c+cc)
15      sum=sum+ca*(x2**3-syl**3)/3.+cb*.5*(x2**2-syl**2)
     .  +cc*(x2-syl)
        ca=a
        cb=b
        cc=c
3       syl=x2
        avint=sum+ca*(xup**3-syl**3)/3.+cb*.5*(xup**2-syl**2)
     .  +cc*(xup-syl)
        if(ind.eq.1) return
        ind=1
        syl=xup
        xup=xlo
        xlo=syl
        avint=-avint
        return
        end
