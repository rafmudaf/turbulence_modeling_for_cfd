      program general
c=======================================================================
c         G E N E R A T E   A   M E S H   A N D   I N I T I A L
c   C O N D I T I O N S   F O R   A   G E N E R A L   G E O M E T R Y
c-----------------------------------------------------------------------
c    Copyright (c) 2006 DCW Industries, Inc.  All rights reserved.
c=======================================================================
      include 'come2c.'
      include 'cpuid.'
c=======================================================================
c     C O N S T A N T S   A N D   I N I T I A L   P A R A M E T E R S
c=======================================================================
      call namin(icpu)
c=======================================
c           Gas Properties
c=======================================
      gm1=gam-1.
      gp1=gam+1.
      ggm1=gam*gm1
      cv=rgas/gm1
c=======================================
c           Flow Properties
c=======================================
      if(iadbwl.eq.1) twall=tinf*(1.+.5*gm1*prt*xma**2)
      eiwall=cv*twall
      rhoinf=pinf/(rgas*tinf)
      cinf=sqrt(gam*rgas*tinf)
      uinf=cinf*xma
      if(abs(vispow-1.5).lt..01) then
        muflag=0
      else
        muflag=1
      endif
      rmuinf=viscos(tinf)
      rhouin=rhoinf*uinf
      eiinf=cv*tinf
      einf=rhoinf*(eiinf+.5*uinf**2)
      if(rl.gt.0.) then
        rel=rhouin*rl/rmuinf
      else
        rl=rel*rmuinf/(uinf*rhoinf)
      endif
c=======================================
c         Counters and Flags
c=======================================
      ilm1=il-1
      jlm1=jl-1
      nadv=0
c=======================================
c     Turbulence Model Parameters
c=======================================
      epinf=retinf*rmuinf
      rhokin=rhoinf*tkeinf
      rhomin=rhoinf*rhokin/epinf
      rhomw=2500./rough**2
      omegin=rhomin/rhoinf
c=======================================
c      Miscellaneous Parameters
c=======================================
      pi=355./113.
      theta=pi*theta/180.
      time=0.
      faxi=.1
      n=1
      m=1
      nend=0
      do 20 i=1,il
        do 10 j=1,jl
          dvol(i,j)=0.
          ep(i,j)=0.
 10     continue
 20   continue
c=======================================
c     Initial Flowfield and Mesh
c=======================================
      write(*,100)
      call mesh
      write(*,110)
      call initl(omegin)
      call gmtry
      call bc
      write(*,120)
      call saveit
      call prntxy
      call prntff
c=======================================
c   Close general.prt and input2c.dat
c=======================================
      close(iunit1)
      close(iunit4)
c=======================================
c    Pause to keep window visible after
c    the run completes...
c=======================================
      write(*,130)
      pause
c=======================================================================
 100  format(/1x,'*** Generating the finite-difference mesh...')
 110  format(/1x,'*** Setting initial conditions...')
 120  format(/1x,'*** Saving EDDY2C input data file to disk...')
 130  format(///1x,'*** Mesh and initial conditions successfully '
     *             'created for Program EDDY2C ***')
c=======================================================================
      end
      subroutine bc
c=======================================================================
c                B O U N D A R Y   C O N D I T I O N S
c=======================================================================
      include 'come2c.'
c=======================================================================
      do 10 i=2,ilm1
        v(i,1)=-v(i,2)
        p(i,1)=p(i,2)
        rho(i,1,m)=rho(i,2,m)
        e(i,1,m)=e(i,2,m)
c===================
c  Symmetry plane
c===================
        if(i.lt.ile) then
          u(i,1)=u(i,2)
          ei(i,1)=ei(i,2)
          rhou(i,1,m)=rhou(i,2,m)
          rhov(i,1,m)=-rhov(i,2,m)
          if(model.eq.0.or.model.eq.1) then
            tke(i,1)=tke(i,2)
            omeg(i,1)=omeg(i,2)
            if(model.eq.1) then
              rhotxx(i,1,m)=rhotxx(i,2,m)
              rhotxy(i,1,m)=rhotxy(i,2,m)
              rhotyy(i,1,m)=rhotyy(i,2,m)
              rhokme(i,1,m)=rhokme(i,2,m)
            endif
          endif
c===================
c  Solid boundary
c===================
        else
          u(i,1)=-u(i,2)
          if(iadbwl.eq.0) then
            ei(i,1)=2.*eiwall-ei(i,2)
          else
            ei(i,1)=ei(i,2)
          endif
          sx=y(i,2)-y(i+1,2)
          sy=x(i+1,2)-x(i,2)
          sj=sqrt(sx*sx+sy*sy)
          sxp=sx/sj
          syp=sy/sj
          rup=-sxp*rhov(i,2,m)+syp*rhou(i,2,m)
          rvp=-syp*rhov(i,2,m)-sxp*rhou(i,2,m)
          rhou(i,1,m)=sxp*rvp+syp*rup
          rhov(i,1,m)=syp*rvp-sxp*rup
          if(model.eq.0.or.model.eq.1) then
            tw=.5*(ei(i,1)+ei(i,2))/cv
            rmu=viscos(tw)
            omegw=2.*rhomw*rmu/(rho(i,1,m)+rho(i,2,m))
            tke(i,1)=-tke(i,2)
            omeg(i,1)=2.*omegw-omeg(i,2)
            if(model.eq.1) then
              rhotxx(i,1,m)=-rhotxx(i,2,m)
              rhotxy(i,1,m)=-rhotxy(i,2,m)
              rhotyy(i,1,m)=-rhotyy(i,2,m)
              rhokme(i,1,m)=-rhokme(i,2,m)
            endif
          endif
        endif
        if(model.eq.0.or.model.eq.1) then
          rhok(i,1,m)=rho(i,1,m)*tke(i,1)
          rhomeg(i,1,m)=rho(i,1,m)*omeg(i,1)
        endif
 10   continue
      do 20 j=1,jl
c===================
c  Symmetry plane
c===================
        if(icase.eq.0) then
          rho(1,j,m)=rho(2,j,m)
          rhou(1,j,m)=rhou(2,j,m)
          rhov(1,j,m)=-rhov(2,j,m)
          e(1,j,m)=e(2,j,m)
          u(1,j)=u(2,j)
          v(1,j)=-v(2,j)
          ei(1,j)=ei(2,j)
          p(1,j)=p(2,j)
          if(model.eq.0.or.model.eq.1) then
            rhok(1,j,m)=rhok(2,j,m)
            rhomeg(1,j,m)=rhomeg(2,j,m)
            tke(1,j)=tke(2,j)
            omeg(1,j)=omeg(2,j)
            if(model.eq.1) then
              rhotxx(1,j,m)=rhotxx(2,j,m)
              rhotxy(1,j,m)=rhotxy(2,j,m)
              rhotyy(1,j,m)=rhotyy(2,j,m)
              rhokme(1,j,m)=rhokme(2,j,m)
            endif
          endif
        endif
c===================
c    Exit plane
c===================
        rho(il,j,m)=rho(ilm1,j,m)
        rhou(il,j,m)=rhou(ilm1,j,m)
        rhov(il,j,m)=rhov(ilm1,j,m)
        e(il,j,m)=e(ilm1,j,m)
        u(il,j)=u(ilm1,j)
        v(il,j)=v(ilm1,j)
        ei(il,j)=ei(ilm1,j)
        p(il,j)=p(ilm1,j)
        if(model.eq.0.or.model.eq.1) then
          rhok(il,j,m)=rhok(ilm1,j,m)
          rhomeg(il,j,m)=rhomeg(ilm1,j,m)
          tke(il,j)=tke(ilm1,j)
          omeg(il,j)=omeg(ilm1,j)
          if(model.eq.1) then
            rhotxx(il,j,m)=rhotxx(ilm1,j,m)
            rhotxy(il,j,m)=rhotxy(ilm1,j,m)
            rhotyy(il,j,m)=rhotyy(ilm1,j,m)
            rhokme(il,j,m)=rhokme(ilm1,j,m)
          endif
        endif
 20   continue
      return
      end
      subroutine gmtry
c=======================================================================
c       C A L C U L A T E   C E L L   C O R N E R   P O I N T S
c                        A N D   V O L U M E S
c=======================================================================
      include 'come2c.'
c=======================================================================
      do 20 j=1,jl
        do 10 i=1,il
          if(iaxi.eq.0) then
            vol(i,j)=(abs((x(i,j)-x(i+1,j))*y(i+1,j+1)
     *              +(x(i+1,j)-x(i+1,j+1))*y(i,j)
     *              +(x(i+1,j+1)-x(i,j))*y(i+1,j))
     *              +abs((x(i,j)-x(i+1,j+1))*y(i,j+1)
     *              +(x(i+1,j+1)-x(i,j+1))*y(i,j)
     *              +(x(i,j+1)-x(i,j))*y(i+1,j+1)))/2.
          else
            x1=x(i,j)
            y1=y(i,j)
            z1=y1*faxi
            x2=x(i,j+1)
            y2=y(i,j+1)
            z2=y2*faxi
            x3=x(i+1,j+1)
            y3=y(i+1,j+1)
            z3=y3*faxi
            x4=x(i+1,j)
            y4=y(i+1,j)
            z4=y4*faxi
            vol1=abs(z1*(x1*(y2-y4)-x2*(y1-y4)+x4*(y1-y2)))
            vol2=abs(z2*(x1*(y2-y3)-x2*(y1-y3)+x3*(y1-y2)))
            vol3=abs(z3*(x2*(y3-y4)-x3*(y2-y4)+x4*(y2-y3)))
            vol4=abs(z4*(x1*(y3-y4)-x3*(y1-y4)+x4*(y1-y3)))
            a234=x2*(y3-y4)-x3*(y2-y4)+x4*(y2-y3)
            a134=x1*(y3-y4)-x3*(y1-y4)+x4*(y1-y3)
            a124=x1*(y2-y4)-x2*(y1-y4)+x4*(y1-y2)
            a123=x1*(y2-y3)-x2*(y1-y3)+x3*(y1-y2)
            vol5=abs(z1*a234+z2*a134+z3*a124+z4*a123)/2.
            vol(i,j)=(vol1+vol2+vol3+vol4+vol5)/3.
          endif
          dvol(i,j)=0.
 10     continue
 20   continue
      return
      end
      subroutine initl(omegin)
c=======================================================================
c                I N I T I A L   P R O F I L E S
c=======================================================================
      include 'come2c.'
      parameter    (jebl=301)
      common/bl/   edim(jebl),edimp(jebl),eldm(jebl),eldmp(jebl),
     *             rdim(jebl),rdimp(jebl),tdim(jebl),tdimp(jebl),
     *             udim(jebl),udimp(jebl),uudm(jebl),uudmp(jebl),
     *             uvdm(jebl),uvdmp(jebl),vvdm(jebl),vvdmp(jebl),
     *             wdim(jebl),wdimp(jebl),ydim(jebl),yplus(jebl),
     *             uplus(jebl)
      common/atmp/ eib(jmx),elb(jmx),omegb(jmx),rhob(jmx),tkeb(jmx),
     *             txxb(jmx),txyb(jmx),tyyb(jmx),ub(jmx),vb(jmx)
      common/extra/yy(jmxp1)
c=======================================================================
c                   I N I T I A L   F L O W F I E L D
c=======================================================================
c Read boundary-layer profiles and then
c   convert to conservative stresses
c=======================================
      ipoint=1
      read(iunit2,140) iedge
      do 10 i=1,iedge
        read(iunit2,150) ydim(i),udim(i),tdim(i),rdim(i),edim(i),
     *                   wdim(i),eldm(i),uudm(i),uvdm(i),vvdm(i)
 10   continue
      if(iaxi.ne.0) ydim(2)=.5*ydim(3)
      do 20 i=1,iedge
        tdim(i)=cv*tdim(i)
        wdim(i)=wdim(i)*ydim(i)**2
        eldm(i)=rdim(i)*eldm(i)
        uudm(i)=rdim(i)*uudm(i)+.66666667*eldm(i)
        uvdm(i)=rdim(i)*uvdm(i)
        vvdm(i)=rdim(i)*vvdm(i)+.66666667*eldm(i)
 20   continue
c=======================================
c        Compute surface slopes
c=======================================
      call slope(iedge)
c=======================================
c    Compute spline-fit coefficients
c=======================================
      call spline(iedge,ydim,udim,udimp,1)
      call spline(iedge,ydim,tdim,tdimp,2)
      call spline(iedge,ydim,rdim,rdimp,3)
      call spline(iedge,ydim,edim,edimp,4)
      call spline(iedge,ydim,wdim,wdimp,5)
      call spline(iedge,ydim,eldm,eldmp,6)
      call spline(iedge,ydim,uudm,uudmp,7)
      call spline(iedge,ydim,uvdm,uvdmp,8)
      call spline(iedge,ydim,vvdm,vvdmp,9)
c=======================================
c   Interpolate profiles onto the mesh
c=======================================
      write(iunit1,100)
      yedge=ydim(iedge)
      do 30 j=2,jl
        if(yy(j).le.yedge) then
          call interp(yy(j),tkeb(j),rhob(j),eib(j),ub(j),omegb(j),
     *                elb(j),txxb(j),txyb(j),tyyb(j),ipoint,iedge)
          if(tkeb(j).le.0.) tkeb(j)=tkeinf
          omegb(j)=omegb(j)/yy(j)**2
        else
          rhob(j)=rhoinf
          ub(j)=uinf
          eib(j)=eiinf
          tkeb(j)=tkeinf
          omegb(j)=omegin
          elb(j)=0.
          txxb(j)=0.
          txyb(j)=0.
          tyyb(j)=0.
        endif
 30     continue
c=======================================
c     Set values at reflected point
c=======================================
        rhob(1)=rhob(2)
        eib(1)=eib(2)
        omegb(1)=omegb(2)
        ub(1)=-ub(2)
        tkeb(1)=-tkeb(2)
        elb(1)=-elb(2)
        txxb(1)=-txxb(2)
        txyb(1)=-txyb(2)
        tyyb(1)=-tyyb(2)
c=======================================
c     Print the upstream profiles
c=======================================
        do 40 j=1,jl
          write(iunit1,110) j,yy(j),ub(j),rhob(j),eib(j),tkeb(j),
     *                   omegb(j),elb(j),txxb(j),txyb(j),tyyb(j)
 40     continue
      do 50 i=1,il
        do 50 j=1,jl
          p(i,j)=pinf
          rho(i,j,1)=rhob(j)
          u(i,j)=ub(j)
          v(i,j)=0.
          ei(i,j)=eib(j)
          tke(i,j)=tkeb(j)
          omeg(i,j)=omegb(j)
          rhokme(i,j,1)=elb(j)
          rhotxx(i,j,1)=txxb(j)
          rhotxy(i,j,1)=txyb(j)
          rhotyy(i,j,1)=tyyb(j)
 50   continue
c=======================================================================
c                B O U N D A R Y   C O N D I T I O N S
c             A T  T H E   T O P   O F   T H E   M E S H
c=======================================================================
      do 80 i=1,il
          p(i,jl)=pinf
          rho(i,jl,1)=rhoinf
          u(i,jl)=uinf
          v(i,jl)=0.
          ei(i,jl)=eiinf
          tke(i,jl)=tkeinf
          omeg(i,jl)=omegin
          rhokme(i,jl,1)=0.
          rhotxx(i,jl,1)=0.
          rhotxy(i,jl,1)=0.
          rhotyy(i,jl,1)=0.
 80   continue
      return
c=======================================================================
 100  format (//40x,'*   I N L E T   P R O F I L E S   *'//1x,
     *      '   j','    y     ','    u     ','   rho    ',
     *      '    ei    ','   tke    ','   omeg   ','  rhok-e  ',
     *      '  rhotxx  ','  rhotxy  ','  rhotyy  ',
     *        /1x,'   -','    -     ','    -     ','   ---    ',
     *      '    --    ','   ---    ','   ----   ','  ------  ',
     *      '  ------  ','  ------  ','  ------  ')
 110  format(1x,i4,1p10e10.2)
 140  format(i6)
 150  format(10e12.4)
c=======================================================================
      end
      subroutine interp(yi,tkei,rhoi,eii,ui,omi,eli,txxi,txyi,tyyi,i,n)
c=======================================================================
c  I N T E R P O L A T E   B O U N D A R Y   L A Y E R   P R O F I L E S
c=======================================================================
      parameter    (jebl=301)
      common/bl/   edim(jebl),edimp(jebl),eldm(jebl),eldmp(jebl),
     *             rdim(jebl),rdimp(jebl),tdim(jebl),tdimp(jebl),
     *             udim(jebl),udimp(jebl),uudm(jebl),uudmp(jebl),
     *             uvdm(jebl),uvdmp(jebl),vvdm(jebl),vvdmp(jebl),
     *             wdim(jebl),wdimp(jebl),ydim(jebl),yplus(jebl),
     *             uplus(jebl)
c=======================================================================
      if(i.le.0) i=1
      if(i.ge.n) i=n-1
      if((ydim(i)-yi)*(ydim(i+1)-yi).le.0.) go to 20
      in=sign(1.,(ydim(i+1)-ydim(i))*(yi-ydim(i)))
 10   if((i+in).le.0) go to 20
      if((i+in).ge.n) go to 20
      i=i+in
      if((ydim(i)-yi)*(ydim(i+1)-yi).gt.0.) go to 10
 20   continue
      hh=ydim(i+1)-ydim(i)
      dels=(yi-ydim(i))/hh
      delsp=(ydim(i+1)-yi)/hh
      phi1=dels**2*(3.-2.*dels)
      phi0=1.-phi1
      psi1=-hh*delsp*dels**2
      psi0=hh*dels*delsp**2
      tkei=edim(i)*phi0+edim(i+1)*phi1+edimp(i)*psi0+edimp(i+1)*psi1
      rhoi=rdim(i)*phi0+rdim(i+1)*phi1+rdimp(i)*psi0+rdimp(i+1)*psi1
      eii =tdim(i)*phi0+tdim(i+1)*phi1+tdimp(i)*psi0+tdimp(i+1)*psi1
      ui  =udim(i)*phi0+udim(i+1)*phi1+udimp(i)*psi0+udimp(i+1)*psi1
      omi =wdim(i)*phi0+wdim(i+1)*phi1+wdimp(i)*psi0+wdimp(i+1)*psi1
      eli =eldm(i)*phi0+eldm(i+1)*phi1+eldmp(i)*psi0+eldmp(i+1)*psi1
      txxi=uudm(i)*phi0+uudm(i+1)*phi1+uudmp(i)*psi0+uudmp(i+1)*psi1
      txyi=uvdm(i)*phi0+uvdm(i+1)*phi1+uvdmp(i)*psi0+uvdmp(i+1)*psi1
      tyyi=vvdm(i)*phi0+vvdm(i+1)*phi1+vvdmp(i)*psi0+vvdmp(i+1)*psi1
      return
      end
      subroutine mesh
c=======================================================================
c       S E T U P   M E S H   C R E A T E D   E X T E R N A L L Y
c=======================================================================
      include 'come2c.'
      common/extra/yy(jmxp1)
      double precision xdp(imxp1,jmx-1),ydp(imxp1,jmx-1)
c=======================================================================
      ilp1=il+1
      jlp1=jl+1
c=======================================
c  Read Mesh Coordinates from grid.dat
c=======================================
      read(iunit3,*) i0,j0
      if(i0.ne.il.or.j0.ne.jl-2) then
        write(*,*) 'Mesh dimensions do not match input values'
        write(*,*) 'i0 = ',i0,'   il   = ',il
        write(*,*) 'j0 = ',j0,'   jl-2 = ',jl-2
        stop
      endif
      read(iunit3,*) ((xdp(i,j),ydp(i,j),j=1,jl-2),i=1,il)
      do 20 i=1,il
        do 10 j=2,jlm1
          x(i,j)=xdp(i,j-1)
          y(i,j)=ydp(i,j-1)
 10     continue
 20   continue
      do 30 j=2,jlm1
        x(ilm1,j)=2.*x(ilm1-1,j)-x(ilm1-2,j)
        y(ilm1,j)=2.*y(ilm1-1,j)-y(ilm1-2,j)
        x(il,j)=2.*x(ilm1,j)-x(ilm1-1,j)
        y(il,j)=2.*y(ilm1,j)-y(ilm1-1,j)
        x(ilp1,j)=2.*x(il,j)-x(ilm1,j)
        y(ilp1,j)=2.*y(il,j)-y(ilm1,j)
 30   continue
c=======================================
c     Compute coordinates of the
c     upper and lower boundaries
c=======================================
      do 40 i=1,ilp1
        x(i,1)=2.*x(i,2)-x(i,3)
        y(i,1)=2.*y(i,2)-y(i,3)
        x(i,jl)=2.*x(i,jlm1)-x(i,jl-2)
        y(i,jl)=2.*y(i,jlm1)-y(i,jl-2)
        x(i,jlp1)=2.*x(i,jl)-x(i,jlm1)
        y(i,jlp1)=2.*y(i,jl)-y(i,jlm1)
 40   continue
c======================================
c  Temporary fix for axisymmetric flows
c======================================
      if(iaxi.ne.0) then
        do 45 i=1,ilp1
          do 42 j=1,jlp1
            y(i,j)=radius+y(i,j)
 42       continue
 45     continue
      endif
c=======================================
c   Compute j-line midpoints at inlet
c=======================================
      do 50 j=1,jl
        yy(j)=.5*(y(1,j)+y(1,j+1))-y(1,2)
 50   continue
c=======================================================================
 100  format(i4,1x,i4)
 110  format(e23.14,2x,e23.14)
c=======================================================================
      end
      subroutine namin(icpu)
c=======================================================================
c           COORDINATE READING INPUT DATA FILE 'EDDY2C.DAT'
c=======================================================================
      parameter (niu=6,nn=22,nxc=6,nxcm=6,nxd=24,nxg=3,nng=2)
      character a*12,newfil*7,pform*9
      common/close/ xc(nxc)
      common/closem/xcm(nxcm)
      common/data/  xd(nxd)
      common/iodef/ iu(niu)
      common/ndata/ n(nn)
      common/gridr/ xg(nxg)
      common/gridi/ ng(nng)
c=======================================================================
c----------------------------------
c         Read units type
c----------------------------------
      call namsys(icpu,iin,iv,msdos,newfil,pform)
      open(iin,file='eddy2c.dat',status='old')
      iu(1)=2
      do 10 i=2,niu-2
        iu(i)=i+5
 10   continue
      iu(6)=10
      io=iu(1)
      open(io,file='general.prt',form=pform,status=newfil)
      write(io,200)
      read(iin,100) a,iu(5)
      write(io,100) a,iu(5)
c----------------------------------
c    Read general input (real)
c----------------------------------
      do 20 i=1,nxd
        read(iin,300) a,xd(i)
        write(io,300) a,xd(i) 
 20   continue
c----------------------------------
c  Read general input (integer)
c----------------------------------
      do 30 i=1,nn
        read(iin,100) a,n(i)
        write(io,100) a,n(i) 
 30   continue
c----------------------------------
c Read k-omega closure coefficients
c----------------------------------
      do 40 i=1,nxc
        read(iin,300) a,xc(i)
        write(io,300) a,xc(i) 
 40   continue
c----------------------------------
c Read Stress-omega closure coeffs
c----------------------------------
      do 50 i=1,nxcm
        read(iin,300) a,xcm(i)
        write(io,300) a,xcm(i) 
 50   continue
c----------------------------------
c    Read mesh parameters (real)
c----------------------------------
      do 60 i=1,nxg
        read(iin,300) a,xg(i)
        write(io,300) a,xg(i) 
 60   continue
c----------------------------------
c  Read mesh parameters (integer)
c----------------------------------
      do 70 i=1,nng
        read(iin,100) a,ng(i)
        write(io,100) a,ng(i) 
 70   continue
      write(io,400)
c----------------------------------
c        Open disk files
c----------------------------------
      open(iu(2),file='profil2c.dat',status='old')
      open(iu(3),file='grid.dat',status='old')
      open(iu(4),file='input2c.dat',form='unformatted',status=newfil)
      return
c=======================================================================
 100  format(1x,a12,i6)
 200  format(1x,'Program GENERAL input data:'/,
     *       1x,'___________________________'/)
 300  format(1x,a12,1pe13.6)
 400  format(1x,'___________________________'//)
c=======================================================================
      end
      subroutine namsys(icpu,iin,iv,msdos,newfil,pform)
c=======================================================================
c  S E T   U P   S Y S T E M - D E P E N D E N T   P A R A M E T E R S
c=======================================================================
      character newfil*7,pform*9
c-----------------------------------------------------------------------
c  Must define the following...
c   iv     = Standard console unit number
c   msdos  = 0...Printer opened as 'prn'
c            1...Don't open 'prn'
c   pform  = 'printer' for SVS Fortran
c          = 'print' for SUN Fortran
c          = 'formatted' for all others
c   newfil = 'new' if compiler writes over an existing file
c          = 'unknown' for Ansi-77 standard operation
c-----------------------------------------------------------------------
c=======================================================================
c----------------------------------
c  Define main input data file LUN
c----------------------------------
      iin=15
c----------------------------------
c       Define console LUN
c----------------------------------
      iv=5
c----------------------------------
c  SVS Fortran (680x0 and 80x86)
c----------------------------------
      if(icpu.eq.0) then
        msdos=0
        pform='printer'
        newfil='new'
c----------------------------------
c  Lahey/Microsoft Fortran (80x86)
c----------------------------------
      elseif(icpu.eq.1) then
        msdos=0
        pform='formatted'
        newfil='unknown'
c----------------------------------
c             VAX/VMS
c----------------------------------
      elseif(icpu.eq.2) then
        msdos=1
        pform='formatted'
        newfil='new'
c----------------------------------
c   SUN Fortran...SUN Workstation
c       and Definicon SPARC
c----------------------------------
      elseif(icpu.eq.3) then
        msdos=0
        pform='print'
        newfil='unknown'
c----------------------------------
c     Cray Fortran...Unicos
c----------------------------------
      elseif(icpu.eq.4) then
        msdos=1
        pform='formatted'
        newfil='unknown'
c----------------------------------
c      Silicon Graphics Iris
c----------------------------------
      elseif(icpu.eq.5) then
        msdos=0
        pform='formatted'
        newfil='unknown'
c----------------------------------
c      Error...say so and quit
c----------------------------------
      else
        write(*,*) 'icpu = ',icpu,' is not supported!!!'
        stop
      endif
      return
      end
      subroutine prntff
c=======================================================================
c                     P R I N T   P R O F I L E S
c=======================================================================
      include 'come2c.'
c=======================================================================
      write(iunit1,100)
      write(iunit1,110) time,nadv
      do 20 i=1,il
        write(iunit1,120) i
        write(iunit1,130)
        do 10 j=1,jl
          cij=sqrt(ggm1*abs(ei(i,j)))
          write(iunit1,140) j,rho(i,j,m),v(i,j),u(i,j),p(i,j),ei(i,j),
     *                      cij
 10     continue
 20   continue
      return
c=======================================================================
 100  format(/28x,'*   F L O W F I E L D   *'//)
 110  format(21x,'Time = ',1pe14.7,2x,'Nadv = ',i9)
 120  format(/8x,'Column  ',i5)
 130  format(5x,'j',5x,'rho',9x,'v',10x,'u',10x,'p',9x,'ei',10x,'c')
 140  format(1x,i5,1p6e11.4)
c=======================================================================
      end
      subroutine prntxy
c=======================================================================
c             P R I N T   M E S H   P A R A M E T E R S
c=======================================================================
      include 'come2c.'
c=======================================================================
      ilp1=il+1
      jlp1=jl+1
      write(iunit1,100)
      do 30 i=1,ilp1,2
        ip1=i+1
        if(ip1.le.ilp1) then
          write(iunit1,110) i,ip1
          do 10 j=1,jlp1
            write(iunit1,120) j,x(i,j),y(i,j),vol(i,j),x(ip1,j),
     *                        y(ip1,j),vol(ip1,j)
 10       continue
        else
          write(iunit1,130) i
          do 20 j=1,jlp1
            write(iunit1,140) j,x(i,j),y(i,j),vol(i,j)
 20       continue
        endif
 30   continue
      return
c=======================================================================
 100  format(/21x,'*   M E S H   C O O R I N A T E S    *'//)
 110  format(/8x,'Column  ',i5,27x,'Column  ',i5/,
     *   5x,'j',5x,'x',10x,'y',9x,'Vol',16x,'x',10x,'y',9x,'Vol')
 120  format(1x,i5,1p3e11.4,7x,3e11.4)
 130  format(/8x,'Column  ',i5/,5x,'j',5x,'x',10x,'y',9x,'Vol')
 140  format(1x,i5,1p3e11.4)
c=======================================================================
      end
      subroutine saveit
c=======================================================================
c      S A V E   F L O W F I E L D   A R R A Y S   O N   D I S K
c=======================================================================
      include 'come2c.'
c=======================================================================
      write(iunit4) nend,time,il,jl,cfl,ntau,iadbwl,eiwall
      write(iunit4) ((rho(i,j,1),u(i,j),v(i,j),ei(i,j),p(i,j),
     *              ep(i,j),vol(i,j),dvol(i,j),i=1,il),j=1,jl)
      ilp1=il+1
      jlp1=jl+1
      write(iunit4) ((x(i,j),y(i,j),i=1,ilp1),j=1,jlp1)
      write(iunit4) (eta(j),j=2,jlp1)
      if(model.eq.0.or.model.eq.1)
     *  write(iunit4) ((tke(i,j),omeg(i,j),i=1,il),j=1,jl)
      if(model.eq.1)
     *  write(iunit4) ((rhotxx(i,j,1),rhotxy(i,j,1),rhotyy(i,j,1),
     *                                rhokme(i,j,1),i=1,il),j=1,jl)
      close(iunit4)
      return
      end
      subroutine slope(iedge)
c=======================================================================
c           C O M P U T E   S U R F A C E   S L O P E S
c=======================================================================
      parameter    (jebl=301)
      common/bl/   edim(jebl),edimp(jebl),eldm(jebl),eldmp(jebl),
     *             rdim(jebl),rdimp(jebl),tdim(jebl),tdimp(jebl),
     *             udim(jebl),udimp(jebl),uudm(jebl),uudmp(jebl),
     *             uvdm(jebl),uvdmp(jebl),vvdm(jebl),vvdmp(jebl),
     *             wdim(jebl),wdimp(jebl),ydim(jebl),yplus(jebl),
     *             uplus(jebl)
c=======================================================================
      dy1=ydim(2)-ydim(1)
      dy2=ydim(3)-ydim(2)
      aa=-(2.*dy1+dy2)/(dy1*(dy1+dy2))
      bb=(dy1+dy2)/(dy1*dy2)
      cc=-dy1/(dy2*(dy1+dy2))
      edimp(1)=aa*edim(1)+bb*edim(2)+cc*edim(3)
      rdimp(1)=aa*rdim(1)+bb*rdim(2)+cc*rdim(3)
      tdimp(1)=aa*tdim(1)+bb*tdim(2)+cc*tdim(3)
      udimp(1)=aa*udim(1)+bb*udim(2)+cc*udim(3)
      wdimp(1)=aa*wdim(1)+bb*wdim(2)+cc*wdim(3)
      eldmp(1)=aa*eldm(1)+bb*eldm(2)+cc*eldm(3)
      uudmp(1)=aa*uudm(1)+bb*uudm(2)+cc*uudm(3)
      uvdmp(1)=aa*uvdm(1)+bb*uvdm(2)+cc*uvdm(3)
      vvdmp(1)=aa*vvdm(1)+bb*vvdm(2)+cc*vvdm(3)
      edimp(iedge)=0.
      rdimp(iedge)=0.
      tdimp(iedge)=0.
      udimp(iedge)=0.
      wdimp(iedge)=0.
      eldmp(iedge)=0.
      uudmp(iedge)=0.
      uvdmp(iedge)=0.
      vvdmp(iedge)=0.
      return
      end
      subroutine spline(n,vari,vard,vardp,li)
c=======================================================================
c      C O M P U T E   S P L I N E   F I T   C O E F F I C I E N T S
c=======================================================================
      parameter    (jebl=301)
      dimension    vari(1),vard(1),vardp(1),p1(jebl),p2(jebl)
      common/iodef/    iunit1 ,  iunit2 ,  iunit3 ,  iunit4 ,  iutype ,
     *                 iuplot
c=======================================================================
c    Make sure the table is in order
c=======================================
      if(vari(2)-vari(1)) 10,10,20
 10   k=li
      write(iunit1,100) j,k,(vari(j),j=1,n),(vard(j),j=1,n)
      stop
 20   do 30 j=2,n
        if(vari(j)-vari(j-1))10,10,30
 30   continue
c=======================================
c    Set up the tridiagonal matrix
c=======================================
      nm=n-1
      p1(1)=vardp(1)
      p2(1)=0.
      do 40 j=2,nm
        h=vari(j)-vari(j-1)
        hp=vari(j+1)-vari(j)
        g=3.*((vard(j+1)-vard(j))/hp**2+(vard(j)-vard(j-1))/h**2)
        a=1./h
        c=1./hp
        b=2.*(a+c)
        bs=b+a*p2(j-1)
        gs=g-a*p1(j-1)
        p1(j)=gs/bs
        p2(j)=-c/bs
 40   continue
c=======================================
c         Invert the matrix
c=======================================
      k=nm
      do 50 j=1,nm
        vardp(k)=p1(k)+p2(k)*vardp(k+1)
        k=k-1
 50   continue
      return
c=======================================================================
 100  format(1x,49hTable below out of order for SPLINE at position  ,
     *       i5/16hx Table number  ,i6,//(8g15.8))
c=======================================================================
      end
      function viscos(temp)
c=======================================================================
c                M O L E C U L A R   V I S C O S I T Y
c=======================================================================
      include 'come2c.'
c=======================================================================
      if(muflag.eq.0) then
        viscos=viscon*sqrt(temp**3)/(temp+su)
      else
        viscos=viscon*temp**vispow/(temp+su)
      endif
      return
      end
