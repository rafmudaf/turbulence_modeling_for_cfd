      program step
c=======================================================================
c      G E N E R A T E   A   M E S H   F O R   A   B A C K W A R D
c                       F A C I N G   S T E P
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
      call mesh
      call initl(omegin)
      call gmtry
      call bc
      call saveit
      call prntxy
      call prntff
c=======================================
c    Close step.prt and input2c.dat
c=======================================
      close(iunit1)
      close(iunit4)
c=======================================
c    Pause to keep window visible after
c    the run completes...
c=======================================
      write(*,100)
      pause
c=======================================================================
 100  format(///1x,'*** Mesh and initial conditions successfully '
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
c           A D J U S T   M E S H   A T   B O U N D A R I E S
c=======================================================================
      include 'come2c.'
      common/extra/yy(jmxp1)
      common/gridr/ xin,xout,ratey
      common/gridi/ iup,istep,ics
      double precision xdp(imxp1,jmx-1),ydp(imxp1,jmx-1)
c=======================================================================
      ilp1=il+1
      jlp1=jl+1
      xmin=xin/rl
      xmax=xout/rl
      ymax=hmax/rl
      etabl=hfm/hmax
c=======================================
c       Compute Mesh Coordinates
c=======================================
      call backstep(il,jl-2,iup,istep,xmin,xmax,ymax,ratey,etabl,jlfm,
     *          xdp,ydp)
      do 20 i=1,il
        do 10 j=2,jlm1
          x(i,j)=rl*sngl(xdp(i,j-1))
          y(i,j)=rl*sngl(ydp(i,j-1))
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
      parameter (niu=6,nn=22,nxc=6,nxcm=6,nxd=24,nxg=3,nng=3)
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
      open(io,file='step.prt',form=pform,status=newfil)
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
      open(iu(4),file='input2c.dat',form='unformatted',status=newfil)
      return
c=======================================================================
 100  format(1x,a12,i6)
 200  format(1x,'Program STEP input data:'/,
     *       1x,'________________________'/)
 300  format(1x,a12,1pe13.6)
 400  format(1x,'________________________'//)
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

c=======================================================================
c  The following routines generate the backward-facing step mesh in
c  double precision...the routines above are the generic routines that
c  adjust the mesh at the boundaries and set initial conditions
c=======================================================================

      subroutine backstep(imax,jmax,iup,istep,xmin,xmax,ymax,ratey,
     *                    etabl,jbl,x,y)
c=======================================================================
c G E N E R A T E   G R I D   F O R   B A C K S T E P   G E O M E T R Y
c=======================================================================
      double precision eta(201),xi(501)
      double precision x(501,201),y(501,201)
      double precision pi,etabl,ratey,ximax
      double precision ymax,ymax0,ymax1,ymax2
      double precision etamax,etamax0,etamax1,etamax2
      double precision xid,xid0,xid1,xid2
      double precision xiu,xiu0,xiu1,xiu2
      double precision xmin,xmin0,xmin1,xmin2
      double precision xmax,xmax0,xmax1,xmax2
      common/const/ pi
c=======================================================================
c=======================================
c    Generate grid in xi,eta space
c=======================================
      pi=4.d0*datan(1.d0)
c-----------------------------------
c First pre-Newton's iteration pass
c-----------------------------------
      xid0=2.*xmax
      xiu0=2.*dabs(xmin)
      etamax0=2.*ymax
      call grid(imax,jmax,iup,istep,xiu0,xid0,etamax0,ratey,
     *          etabl,jbl,xi,eta,1)
      call map(xi(imax),eta(1),xmax0,ydummy)
      call map(xi(1),eta(1),xmin0,ydummy)
      call map(xi(1),eta(jmax),xdummy,ymax0)
c-----------------------------------
c Second pre-Newton's iteration pass
c-----------------------------------
      factor=xmax/xmax0
      xid1=xid0*factor
      factor=dabs(xmin/xmin0)
      xiu1=xiu0*factor
      factor=ymax/ymax0
      etamax1=etamax0*factor
      call grid(imax,jmax,iup,istep,xiu1,xid1,etamax1,ratey,
     *          etabl,jbl,xi,eta,1)
      call map(xi(imax),eta(1),xmax1,ydummy)
      call map(xi(1),eta(1),xmin1,ydummy)
      call map(xi(1),eta(jmax),xdummy,ymax1)
c-----------------------------------
c      Newton's iterations to
c    determine xi and eta that
c   yield the precise grid size
c-----------------------------------
      maxit=20
      do 10 i=1,20
        xid2=xid0+(xmax-xmax0)*(xid1-xid0)/(xmax1-xmax0)
        xiu2=xiu0+(xmin-xmin0)*(xiu1-xiu0)/(xmin1-xmin0)
        etamax2=etamax0+(ymax-ymax0)*(etamax1-etamax0)/(ymax1-ymax0)
        call grid(imax,jmax,iup,istep,xiu2,xid2,etamax2,ratey,
     *            etabl,jbl,xi,eta,1)
        call map(xi(imax),eta(1),xmax2,ydummy)
        call map(xi(1),eta(1),xmin2,ydummy)
        call map(xi(1),eta(jmax),xdummy,ymax2)
        dxmax=xmax2-xmax1
        dxmin=xmin2-xmin1
        detamax=etamax2-etamax1
        if(abs(dxmin).lt.1.e-6*abs(xmin).and.
     *     abs(dxmax).lt.1.e-6*abs(xmax).and.
     *     abs(detamax).lt.1.e-6*ymax) go to 20
        xid0=xid1
        xiu0=xiu1
        etamax0=etamax1
        xmax0=xmax1
        xmin0=xmin1
        ymax0=ymax1
        xid1=xid2
        xiu1=xiu2
        etamax1=etamax2
        xmax1=xmax2
        xmin1=xmin2
        ymax1=ymax2
 10   continue
 20   xid=xid2
      xiu=xiu2
      etamax=etamax2
      if(i.lt.10) then
        write(*,100) i
      elseif(i.lt.20) then
        write(*,110) i
      else
        write(*,120) maxit
      endif
      call grid(imax,jmax,iup,istep,xiu,xid,etamax,ratey,
     *          etabl,jbl,xi,eta,0)
c=======================================
c     Apply the Schwarz-Christoffel
c  transformation from: xi,eta --> x,y
c=======================================
      do 40 i=1,imax
        do 30 j=1,jmax
          call map(xi(i),eta(j),x(i,j),y(i,j))
 30     continue
 40   continue
      ybar=y(1,1)+2.*etabl*ymax
      do 50 j=jbl+1,jmax
        if(y(1,j).gt.ybar) go to 60
 50   continue
 60   jbar=j-1
c=======================================
c        Adjust upper boundary
c=======================================
      call adjustu(imax,jmax,iup,istep,ymax,jbar,x,y)
      write(*,130) x(1,1),x(imax,1),y(1,jmax),y(1,2)-y(1,1)
c=======================================
c    Make sure base lies at x = 0
c=======================================
      do 70 i=iup,iup+istep-1
        x(i,1)=0.d0
 70   continue
      return
c=======================================================================
 100  format(' Grid generation required ',i1,' iterations to converge')
 110  format(' Grid generation required ',i2,' iterations to converge')
 120  format(' Grid generation failed to converge after ',i2,
     *       ' iterations')
 130  format(//' Backstep grid dimensions...'
     *        /' Left  boundary is located at  x/h = ',1pe10.3,
     *        /' Right Boundary is located at  x/h = ',1pe10.3,
     *        /' Upper boundary is located at  y/h = ',1pe10.3,
     *        /' Cell nearest the surface has dy/h = ',1pe10.3)
 200  format(i4,2x,i4)
c=======================================================================
      end
      subroutine adjustu(imax,jmax,iup,istep,ymax,jbar,x,y)
c=======================================================================
c   M A K E   U P P E R   G R I D   B O U N D A R Y   S T R A I G H T
c=======================================================================
      double precision x(501,201),y(501,201)
      double precision dyn,dyo,ratio,ratmax,ratmin,ymax
c=======================================================================
      ratmin=100.
      ratmax=0.
      do 20 i=1,imax
        dyo=y(i,jmax)-y(i,jbar)
        dyn=ymax-y(i,jbar)
        ratio=dyn/dyo
        ratmax=dmax1(ratmax,ratio)
        ratmin=dmin1(ratmin,ratio)
        do 10 j=jbar+1,jmax
          y(i,j)=y(i,jbar)+ratio*(y(i,j)-y(i,jbar))
 10     continue
 20   continue
      return
      end
      subroutine drtni(x,f,fp,fct,xst,eps,iend,ier)
c=======================================================================
c  Solve f(x) = 0 using Newton's method...double precision
c=======================================================================
      double precision a,dx,eps,f,fp,tol,tolf,x,xst
c=======================================================================
c  Prepare iteration
c===================
      ier=0
      x=xst
      tol=x
      call fct(tol,f,fp)
      tolf=100.*eps
c======================
c  Start iteration loop
c======================
      do 10 i=1,iend
c=============
c  Converged ?
c=============
        if(f.eq.0.) go to 20
c=============
c  df/dx = 0 ?
c=============
        if(fp.eq.0.) then
          ier=2
          go to 20
        endif
c=======================
c  Iteration is possible
c=======================
        dx=-f/fp
        x=x+dx
        tol=x
        call fct(tol,f,fp)
c==================
c  Convergence test
c==================
        tol=eps
        a=dabs(x)
        if(a.gt.1.) tol=tol*a
        if(dabs(dx).le.tol) then
          if(dabs(f).le.tolf) go to 20
        endif
 10   continue
c============================================
c  If no convergence after iend iterations...
c  return with ier set to 1
c============================================
      ier=1
 20   return
      end
      subroutine fct(x,f,fp)
c=======================================================================
c     F I N I T E   D I F F E R E N C E   G R I D   F U N C T I O N
c=======================================================================
      common/fctcom/   aa,aam,jend,jendm
      double precision aa,aam
      double precision f,fp,r1,r2,x
c=======================================================================
      f=-aam
      fp=aa
      r1=aa
      do 10 j=1,jendm-1
        r2=r1*x
        f=f+r2
        fp=fp+(j+1)*r2
        r1=r2
 10   continue
      f=f+r1*x
      return
      end
      subroutine grid(imax,jmax,iup,istep,xiu,xid,etamax,ratey,
     *                etabl,jbl,xi,eta,iflag)
c=======================================================================
c  G E N E R A T E   T H E   F I N I T E   D I F F E R E N C E   G R I D
c=======================================================================
      external fct
      common/fctcom/   aa,aam,jend,jendm
      double precision eta(500),xi(500)
      double precision aa,aam,dxi,dxi0,etabl,etamax,rate,ratei,ximax,
     *                 xiu,xid,ratey
      double precision denom,etatop,f,fp
c=======================================================================
      idown=imax-istep-iup
      ximax=xiu
      dxi0=1.d0/(istep-1)
      aa=dxi0
      aam=ximax-aa
      jend=iup
      jendm=jend-1
c=================================
c  First solve for upstream RATE
c=================================
      ratei=1.d0+5.d0/jmax
      call drtni(rate,f,fp,fct,ratei,1.d-6,500,ier)
      if(ier.eq.1) then
        write(*,100)
        pause
        stop
      elseif(ier.eq.2) then
        write(*,110)
        pause
        stop
      endif
      if(iflag.eq.0) write(*,120) rate
      dxi=dxi0
      xi(iup-1)=-dxi
      do 10 i=iup-2,1,-1
        dxi=dxi*rate
        xi(i)=xi(i+1)-dxi
 10   continue
      xi(iup)=0.d0
      do 20 i=iup+1,iup+istep
        xi(i)=xi(i-1)+dxi0
 20   continue
c=================================
c  Now solve for downstream RATE
c=================================
      ximax=xid
      dxi0=1.d0/(istep-1)
      aa=dxi0
      aam=ximax-aa
      jend=idown
      jendm=jend-1
      ratei=1.d0+5.d0/jmax
      call drtni(rate,f,fp,fct,ratei,1.d-6,500,ier)
      if(ier.eq.1) then
        write(*,130)
        pause
        stop
      elseif(ier.eq.2) then
        write(*,140)
        pause
        stop
      endif
      if(iflag.eq.0) write(*,150) rate
      dxi=dxi0
      xi(iup+istep+1)=xi(iup+istep)+dxi
      do 30 i=iup+istep+2,imax
        dxi=dxi*rate
        xi(i)=xi(i-1)+dxi
 30   continue
c=================================
c    Set up the vertical grid
c=================================
      jtop=jmax-1
      etatop=etamax
      rate=ratey
      denom=1.
      aam=1.      
      do 40 j=1,jbl-2
        aam=aam*rate
        denom=denom+aam
 40   continue
      deta=etabl*etamax/denom
      eta(1)=0.d0
      eta(2)=deta
      aa=deta
      do 50 j=3,jbl
        aa=aa*rate
        eta(j)=eta(j-1)+aa
 50   continue
      etabar=etatop-eta(jbl-1)
      deta=eta(jbl)-eta(jbl-1)
      aa=deta
      aam=etabar-aa
      jend=jtop-jbl+2
      jendm=jend-1
      ratei=1.d0+5.d0/(jend+1)
      call drtni(rate,f,fp,fct,ratei,1.d-6,500,ier)
      if(ier.eq.1) then
        write(*,160)
        pause
        stop
      elseif(ier.eq.2) then
        write(*,170)
        pause
        stop
      endif
      if(iflag.eq.0) then
        write(*,180) ratey
        write(*,190) rate
      endif
      do 60 j=jbl+1,jtop
        deta=deta*rate
        eta(j)=eta(j-1)+deta
 60   continue
      eta(jtop+1)=etatop
      return
c=======================================================================
 100  format(//' Newton iterations failed:'/'Convergence not attained',
     *         ' in computing stretching rate upstream of the step')
 110  format(//' Newton iterations failed:'/'df/dx = 0 encountered',
     *         ' in computing stretching rate upstream of the step')
 120  format(//' Stretching rate for points upstream of the corner',
     *         '    = ',1pe9.3)
 130  format(//' Newton iterations failed:'/'Convergence not attained',
     *         ' in computing stretching rate upstream of the step')
 140  format(//' Newton iterations failed:'/'df/dx = 0 encountered',
     *         ' in computing stretching rate upstream of the step')
 150  format(' Stretching rate for points downstream of the corner',
     *       '  = ',1pe9.3)
 160  format(//' Newton iterations failed:'/'Convergence not attained',
     *         ' in computing stretching rate above boundary layer')
 170  format(//' Newton iterations failed:'/'df/dx = 0 encountered',
     *         ' in computing stretching rate above boundary layer')
 180  format(' Stretching rate for points within the boundary layer',
     *       ' = ',1pe9.3)
 190  format(' Stretching rate for points above the boundary layer',
     *       '  = ',1pe9.3)
c=======================================================================
      end
      subroutine map(xi,eta,x,y)
c=======================================================================
c     M A P   F R O M   Z E T A   P L A N E   T O   Z   P L A N E
c=======================================================================
      complex*16 zeta,root1,root2,zbar1,zbar2,z
      double precision eta,xi
      double precision pi,x,y
      common/const/ pi
c=======================================================================
      zeta=dcmplx(xi,eta)
      root1=cdsqrt(zeta)
      root2=cdsqrt(zeta-1.)
      zbar1=2./pi*root1*root2
      zbar2=2./pi*cdlog(root1-root2)
      z=zbar1-zbar2
      x=dimag((0,1)*z)
      y=dimag(z)
      return
      end
