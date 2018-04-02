      program defect
c=======================================================================
c        I N C O M P R E S S I B L E   D E F E C T   L A Y E R
c
c    Copyright (c) 2006 DCW Industries, Inc.  All rights reserved.
c-----------------------------------------------------------------------
c            MODEL = 0.........................k-omega Model
c            MODEL = 2, KERNG = 0...Standard k-epsilon Model
c            MODEL = 2, KERNG > 0........RNG k-epsilon Model
c            MODEL = 3...................Baldwin-Barth Model
c            MODEL = 4................Spalart-Allmaras Model
c            MODEL = 99...................User-Defined Model
c=======================================================================
      parameter    (jmx=201)
      common/bcom/ climit(jmx),deta(jmx),dudxi(jmx),dxi(jmx),en(jmx),
     *             eps(jmx),et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),
     *             v(jmx),vextra(jmx),wt(jmx),wtm(jmx),xi(jmx),y1(jmx),
     *             y2(jmx),y3(jmx)
c-----------------------------------------------------------------------
      common/data/     betat
      common/extra/    acon   ,  akappa ,  cflmax ,  dti    ,  eone   ,
     *                 epsin  ,  epsone ,  error  ,  etamax ,  etamin ,
     *                 etin   ,  jend   ,  nstep  ,  pitld  ,  umax   ,
     *                 umin   ,  u1     ,  uone   ,  vcon   ,  wcon   ,
     *                 wtin   ,  wtmax  ,  wtone  ,  ximax  ,  ximin  ,
     *                 xk
      common/iodef/    iunit1 ,  iunit2 ,  model  ,  kerng
      common/ndata/    jmax   ,  maxn   ,  nedit  ,  nfreq
      common/oneeqn/   acmax  ,  biga   ,  bigb   ,  bigc   ,  cw1    ,
     *                 cw2    ,  cw3    ,  en1    ,  fwo
      common/rng/      ceps2p ,  phi    ,  rbeta  ,  rlamda
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c-----------------------------------------------------------------------
      include 'cpuid.'
c=======================================================================
c           Read input data
c=======================================
      call namin(icpu)
c=======================================
c            Initialization
c=======================================
      write(iunit1,100)
c=======================================
c    Calculate auxilliary constants
c=======================================
      jend=jmax-1
      vcon=1.+2.*betat
      cflmax=10.+jend/10.
      if(betat.lt.0.) cflmax=.4*cflmax
      error=1.e-5
      etamin=.0001
      acmax=50.
      umax=100.
c=======================================
c         Generate the grid
c=======================================
      call start(0)
      ximin=dlog(etamin)/akappa
      ximax=40.*exp(-betat/3.)-10.
      xk=1.
      call grid
c=======================================
c        Set initial conditions
c=======================================
      call start(1)
      wtmax=100.*wtone
      nstep=0
c=======================================
c   Identify flow and turbulence model
c=======================================
      write(*,110) betat
      if(model.eq.0) then
        umin=.05/jend
        write(*,120)
      elseif(model.eq.2) then
        umin=.001
        if(kerng.eq.0.) then
          write(*,130)
        else
          write(*,135)
        endif
      elseif(model.eq.3) then
        umin=.05/jend
        write(*,140)
      elseif(model.eq.4) then
        umin=.05/jend
        cflmax=5.+3.*(jend-100.)/20.
        if(betat.lt.0.) cflmax=.3*cflmax
        write(*,150)
      elseif(model.eq.99) then
        umin=.001
        write(*,160)
      endif
      write(*,100)
c=======================================================================
c                 BEGINNING OF MAIN INTEGRATION LOOP
c=======================================================================
 10   nstep=nstep+1
c=======================================
c Calculate quantities needed in EDDY
c=======================================
      call geteta
      call calcs
      call tmestp
c=======================================
c Calculate velocity, turbulence energy
c        and dissipation rate
c=======================================
      call eddy
c=======================================
c          Check convergence
c=======================================
      erro=0.
      do 20 j=1,jmax
        err=abs(u(j)-um(j))
        if(err.gt.erro) erro=err
 20   continue
c=======================================
c          Converged solution
c=======================================
      if(erro.le.error) then
        call geteta
        write(iunit1,170) nstep,erro,acon,pitld
        write(*,170) nstep,erro,acon,pitld
        write(*,180)
c=======================================
c       Un-converged solution
c=======================================
      else
        if(mod(nstep,nfreq).eq.0) then
          write(iunit1,170) nstep,erro,acon,pitld
          write(*,170) nstep,erro,acon,pitld
        endif
        if(nstep.ne.maxn) then
c=======================================
c          Print profiles
c=======================================
          if(mod(nstep,nedit).eq.0) call edit
c=======================================
c      Prepare for next timestep
c=======================================
          do 30 j=1,jmax
            etm(j)=et(j)
            um(j)=u(j)
            wtm(j)=wt(j)
            if(model.eq.3.or.model.eq.4) eps(j)=en(j)
 30       continue
          go to 10
c=======================================================================
c                   END OF MAIN INTEGRATION LOOP
c=======================================================================
        else
          write(iunit1,190) nstep
          write(*,190) nstep
        endif
      endif
      call edit
      close(iunit1)
c=======================================
c    Write plotting data disk file
c=======================================
      write(iunit2,200) jmax,betat
      write(iunit2,210) (eta(j),u(j),j=1,jmax)
      endfile(iunit2)
      close(iunit2)
c=======================================
c    Pause to keep window visible after
c    the run completes...
c=======================================
      pause
c=======================================================================
 100  format(1x,'Convergence history:'/)
 110  format(' DEFECT LAYER...Betat = ',1pe10.3)
 120  format(' k-omega Model')
 130  format(' Standard k-epsilon Model')
 135  format(' RNG k-epsilon Model')
 140  format(' Baldwin-Barth Model')
 150  format(' Spalart-Allmaras Model')
 160  format(' User-Defined Model')
 170  format(1x,'Iteration ',i5,'...Max Error =',1pe10.3,
     *          ':  A =',e10.3,', Pi =',e10.3)
 180  format(/19x,' ***** SOLUTION CONVERGED *****'/)
 190  format(/8x,'Computation terminated for max iteration',i6/)
 200  format(1x,i4,1pe14.6)
 210  format(1p2e14.6)
c=======================================================================
      end
      subroutine calcs
c=======================================================================
c C A L C U L A T E   Q U A N T I T I E S   N E E D E D   I N   E D D Y
c=======================================================================
      parameter    (jmx=201)
      common/bcom/ climit(jmx),deta(jmx),dudxi(jmx),dxi(jmx),en(jmx),
     *             eps(jmx),et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),
     *             v(jmx),vextra(jmx),wt(jmx),wtm(jmx),xi(jmx),y1(jmx),
     *             y2(jmx),y3(jmx)
c-----------------------------------------------------------------------
      common/extra/    acon   ,  akappa ,  cflmax ,  dti    ,  eone   ,
     *                 epsin  ,  epsone ,  error  ,  etamax ,  etamin ,
     *                 etin   ,  jend   ,  nstep  ,  pitld  ,  umax   ,
     *                 umin   ,  u1     ,  uone   ,  vcon   ,  wcon   ,
     *                 wtin   ,  wtmax  ,  wtone  ,  ximax  ,  ximin  ,
     *                 xk
      common/iodef/    iunit1 ,  iunit2 ,  model  ,  kerng
      common/ndata/    jmax   ,  maxn   ,  nedit  ,  nfreq
      common/oneeqn/   acmax  ,  biga   ,  bigb   ,  bigc   ,  cw1    ,
     *                 cw2    ,  cw3    ,  en1    ,  fwo
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c=======================================================================
c   Impose the integral constraint
c=======================================
      elog=dlog(etamin)
      acon=um(1)+(1.+u1*etamin)*dlog(etamin)/akappa
      uint=0.
      do 10 j=1,jend
        uint=uint+.5*(eps(j)*um(j)+eps(j+1)*um(j+1))*dxi(j)
 10   continue
      uintp=uint+(acon-(elog-1.+.25*u1*etamin*(2.*elog-1.))/akappa)
     *          *etamin
      renorm=1./uintp
      do 20 j=1,jend
        um(j)=renorm*um(j)
 20   continue
      uint=uint/uintp
      acon=(1.-uint)/etamin+(elog-1.+.25*u1*etamin*(2.*elog-1.))/akappa
      if(acon.lt.-5.) acon=-5.
      if(acon.gt.acmax) acon=acmax
      uone=acon-elog*(1.+u1*etamin)/akappa
c=======================================
c     Compute Coles' wake strength
c=======================================
      do 25 j=2,jmax
        if(u(j).le..01) then
          etaedg=eta(j)
          go to 28
        endif
 25   continue
 28   continue
      pitld=.5*(akappa*acon-dlog(etaedg))
c=======================================
c     Compute eddy viscosity, etc.
c=======================================
      eps(jmax)=epsin
      dudxi(jmax)=0.
      v(jmax)=-vcon*etamax
      vextra(jmax)=0.
      eps(1)=epsone
      dudxi(1)=-1.
      v(1)=-vcon*etamin
      vextra(1)=0.
      do 30 j=2,jend
        vextra(j)=0.
        climit(j)=1.
        if(model.eq.3.or.model.eq.4) then
          eps(j)=en(j)
        else
          eps(j)=etm(j)/wtm(j)
          if(model.eq.2) eps(j)=etm(j)*eps(j)
        endif
        dudxi(j)=(um(j+1)-um(j-1))/(dxi(j)+dxi(j-1))
        v(j)=-vcon*eta(j)
        if(model.eq.0.and.j.gt.5) then
          if(abs(clim*dudxi(j)).gt.etm(j)) then
            climit(j)=etm(j)/abs(clim*dudxi(j))
          endif
        endif
 30   continue
      climit(1)=climit(2)
      climit(jmax)=climit(jend)
c=======================================
c      k-omega model computations
c=======================================
      if(model.eq.0) then
        do 40 j=2,jend
          dkdxi=(1./et(j+1)-1./et(j-1))/(dxi(j-1)+dxi(j))
          dwdxi=dlog(wt(j+1)/wt(j-1))/(dxi(j-1)+dxi(j))
          chik=-alphas*dkdxi*dwdxi
          if(chik.gt.0.) then
            vextra(j)=sigmad*dlog(et(j+1)/et(j-1))/(dxi(j-1)+dxi(j))
          endif
 40     continue
      endif
      return
      end
      subroutine eddy
c=======================================================================
c      A D V A N C E   T H E   S O L U T I O N   I N   T I M E
c=======================================================================
      parameter    (jmx=201)
      dimension    a(jmx),b(jmx),c(jmx)
      common/bcom/ climit(jmx),deta(jmx),dudxi(jmx),dxi(jmx),en(jmx),
     *             eps(jmx),et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),
     *             v(jmx),vextra(jmx),wt(jmx),wtm(jmx),xi(jmx),y1(jmx),
     *             y2(jmx),y3(jmx)
c-----------------------------------------------------------------------
      common/data/     betat
      common/extra/    acon   ,  akappa ,  cflmax ,  dti    ,  eone   ,
     *                 epsin  ,  epsone ,  error  ,  etamax ,  etamin ,
     *                 etin   ,  jend   ,  nstep  ,  pitld  ,  umax   ,
     *                 umin   ,  u1     ,  uone   ,  vcon   ,  wcon   ,
     *                 wtin   ,  wtmax  ,  wtone  ,  ximax  ,  ximin  ,
     *                 xk
      common/iodef/    iunit1 ,  iunit2 ,  model  ,  kerng
      common/ndata/    jmax   ,  maxn   ,  nedit  ,  nfreq
      common/oneeqn/   acmax  ,  biga   ,  bigb   ,  bigc   ,  cw1    ,
     *                 cw2    ,  cw3    ,  en1    ,  fwo
      common/rng/      ceps2p ,  phi    ,  rbeta  ,  rlamda
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c=======================================================================
      if(model.eq.3.or.model.eq.4) go to 50
c=======================================
c          Dissipation rate
c=======================================
      do 10 j=2,jend
        prodw=wcon*eps(j)-beta*etm(j)/alphas
c-----------------------------
c     RNG k-epsilon term          
c-----------------------------
        if(model.eq.2.and.kerng.ne.0) then
          rdenom=rbeta*abs(dudxi(j))**3/betas+alphas*etm(j)**3
          rrng=(alphas*etm(j)-abs(dudxi(j))/rlamda)
     *        *abs(dudxi(j))**3/rdenom
          prodw=prodw-rrng
        endif
c-----------------------------
        if(etm(j).ne.0.) prodw=prodw+alpha*alphas*climit(j)
     *                                  *dudxi(j)**2/etm(j)
        a(j)=sigma*y3(j)+(v(j)-vextra(j))*y2(j)
        b(j)=-sigma*(y1(j)+y3(j))-2.*dti
        c(j)=sigma*y1(j)-(v(j)-vextra(j))*y2(j)
        wt(j)=-a(j)*wtm(j-1)-(b(j)+4.*dti)*wtm(j)-c(j)*wtm(j+1)
        if(prodw.gt.0.) then
          wt(j)=wt(j)-2.*prodw*wtm(j)
        else
          b(j)=b(j)+2.*prodw
        endif
 10   continue
      call tri(a,b,c,wt,1,jmax,0,wtone,0,wtin)
      do 20 j=2,jend
        if(wt(j).lt.0.) wt(j)=.5*(wtm(j-1)+wtm(j+1))
        if(wt(j).gt.wtmax) wt(j)=wtmax
 20   continue
c=======================================
c      Turbulence kinetic energy
c=======================================
      do 30 j=2,jend
        prode=2.*betat*eps(j)-alphas*etm(j)
        if(etm(j).ne.0.) prode=prode+alphas*climit(j)
     *                            *dudxi(j)**2/etm(j)
        a(j)=sigmas*y3(j)+v(j)*y2(j)
        b(j)=-sigmas*(y1(j)+y3(j))-2.*dti
        c(j)=sigmas*y1(j)-v(j)*y2(j)
        et(j)=-a(j)*etm(j-1)-(b(j)+4.*dti)*etm(j)-c(j)*etm(j+1)
        if(prode.gt.0.) then
          et(j)=et(j)-2.*prode*etm(j)
        else
          b(j)=b(j)+2.*prode
        endif
 30   continue
      call tri(a,b,c,et,1,jmax,0,eone,0,etin)
      do 40 j=2,jend
        if(et(j).lt.0.) et(j)=.5*(etm(j-1)+etm(j+1))
 40   continue
      go to 80
c=======================================
c         One-equation models
c=======================================
 50   do 60 j=2,jend
        sigmam=biga+.5*bigb*(1.+eps(j-1)/eps(j))
        sigmap=biga+.5*bigb*(1.+eps(j+1)/eps(j))
        if(model.eq.3) then
          if(sigmam.lt.0.) sigmam=biga/3.
          if(sigmap.lt.0.) sigmap=biga/3.
        endif
        prodn=-wcon*eps(j)-bigc*dudxi(j)
c-----------------------------
c    Spalart-Allmaras term
c-----------------------------
        if(model.eq.4) then
          if(abs(dudxi(j)).lt.1.e-10) then
            fw=fwo
          else
            rr=-(eps(j)/(akappa*eta(j)))**2/dudxi(j)
            ggocw3=(rr+cw2*(rr**6-rr))/cw3
            if(ggocw3.lt..25) then
              fw=fwo*ggocw3*(1.-ggocw3**6/6.)
            elseif(ggocw3.gt.4.) then
              fw=fwo*(1.-1./(6.*ggocw3**6))
            else
              fw=fwo/(1.+1./ggocw3**6)**(1./6.)
            endif
          endif
          prodn=prodn-cw1*fw*(eps(j)/eta(j))**2
        endif
c-----------------------------
        a(j)=sigmam*y3(j)+v(j)*y2(j)
        b(j)=-sigmap*y1(j)-sigmam*y3(j)-2.*dti
        c(j)=sigmap*y1(j)-v(j)*y2(j)
        en(j)=-a(j)*eps(j-1)-(b(j)+4.*dti)*eps(j)-c(j)*eps(j+1)
        if(prodn.gt.0.) then
          en(j)=en(j)-2.*prodn*eps(j)
        else
          b(j)=b(j)+2.*prodn
        endif
 60   continue
      call tri(a,b,c,en,1,jmax,0,epsone,0,epsin)
      do 70 j=2,jend
        if(en(j).lt.0.) en(j)=.5*(eps(j-1)+eps(j+1))
 70   continue
c=======================================
c              Momentum
c=======================================
 80   do 90 j=2,jend
        produ=2.*betat*eps(j)
        dalfdx=y2(j)*(climit(j+1)-climit(j-1))
        a(j)=climit(j)*y3(j)+(v(j)-dalfdx)*y2(j)
        b(j)=-climit(j)*y1(j)-climit(j)*y3(j)-2.*dti
        c(j)=climit(j)*y1(j)-(v(j)-dalfdx)*y2(j)
        u(j)=-a(j)*um(j-1)-(b(j)+4.*dti)*um(j)-c(j)*um(j+1)
        if(produ.gt.0.) then
          u(j)=u(j)-2.*produ*um(j)
        else
          b(j)=b(j)+2.*produ
        endif
 90   continue
      call tri(a,b,c,u,1,jmax,0,uone,0,0.)
c=======================================
c      Reduce ximax as required
c=======================================
      jm=jend/10
      if(u(1).gt.umax.or.u(jend-jm).lt.umin) then
        ximax=ximax-.01*(ximax-ximin)
        call grid
      endif
c=======================================
c     Increase ximax as required
c=======================================
      if(u(jend).gt.umin) then
        ximax=ximax+.01*(ximax-ximin)
        call grid
      endif
      return
      end
      subroutine edit
c=======================================================================
c                   P R I N T I N G   R O U T I N E
c=======================================================================
      parameter    (jmx=201)
      common/bcom/ climit(jmx),deta(jmx),dudxi(jmx),dxi(jmx),en(jmx),
     *             eps(jmx),et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),
     *             v(jmx),vextra(jmx),wt(jmx),wtm(jmx),xi(jmx),y1(jmx),
     *             y2(jmx),y3(jmx)
      common/fcom/ el(jmx),wtprnt(jmx)
c-----------------------------------------------------------------------
      common/extra/    acon   ,  akappa ,  cflmax ,  dti    ,  eone   ,
     *                 epsin  ,  epsone ,  error  ,  etamax ,  etamin ,
     *                 etin   ,  jend   ,  nstep  ,  pitld  ,  umax   ,
     *                 umin   ,  u1     ,  uone   ,  vcon   ,  wcon   ,
     *                 wtin   ,  wtmax  ,  wtone  ,  ximax  ,  ximin  ,
     *                 xk
      common/iodef/    iunit1 ,  iunit2 ,  model  ,  kerng
      common/ndata/    jmax   ,  maxn   ,  nedit  ,  nfreq
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c=======================================================================
      write(iunit1,100) nstep
c=======================================
c        One-equation models
c=======================================
      if(model.eq.3.or.model.eq.4) then
        write(iunit1,110)
        do 10 j=1,jmax
          write(iunit1,120) j,eta(j),u(j),eps(j)
 10     continue
c=======================================
c        Two-equation models
c=======================================
      else
        do 20 j=1,jmax
          el(j)=sqrt(abs(et(j)))/wt(j)
          if(model.eq.0) then
            wtprnt(j)=wt(j)
          elseif(model.eq.2) then
            wtprnt(j)=wt(j)/(betas*et(j))
            el(j)=et(j)*el(j)
          endif
 20     continue
        write(iunit1,130)
        do 30 j=1,jmax
          write(iunit1,140) j,eta(j),u(j),et(j),wtprnt(j),eps(j),el(j)
 30     continue
      endif
      return
c=======================================================================
 100  format(1h1,4x,'Timestep  ',i6)
 110  format(/5x,'j',5x,'eta',10x,'u',10x,'nut')
 120  format(i6,1p3e12.5)
 130  format(/5x,'j',5x,'eta',10x,'u',11x,'k',9x,'omega',8x,'nut',
     *       10x,'L')
 140  format(i6,1p6e12.5)
c=======================================================================
      end
      subroutine geteta
c=======================================================================
c                      C O M P U T E   E T A
c=======================================================================
      parameter    (jmx=201)
      common/bcom/ climit(jmx),deta(jmx),dudxi(jmx),dxi(jmx),en(jmx),
     *             eps(jmx),et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),
     *             v(jmx),vextra(jmx),wt(jmx),wtm(jmx),xi(jmx),y1(jmx),
     *             y2(jmx),y3(jmx)
c-----------------------------------------------------------------------
      common/extra/    acon   ,  akappa ,  cflmax ,  dti    ,  eone   ,
     *                 epsin  ,  epsone ,  error  ,  etamax ,  etamin ,
     *                 etin   ,  jend   ,  nstep  ,  pitld  ,  umax   ,
     *                 umin   ,  u1     ,  uone   ,  vcon   ,  wcon   ,
     *                 wtin   ,  wtmax  ,  wtone  ,  ximax  ,  ximin  ,
     *                 xk
      common/ndata/    jmax   ,  maxn   ,  nedit  ,  nfreq
c=======================================================================
      eta(1)=etamin
      eps1=eps(1)
      do 10 j=2,jend
        eps2=eps(j)
        eta(j)=eta(j-1)+.5*(eps1+eps2)*dxi(j-1)
        deta(j-1)=eta(j)-eta(j-1)
        eps1=eps2
 10   continue
      deta(jend)=deta(jend-1)
      eta(jmax)=eta(jend)+deta(jend)
      etamax=eta(jmax)
      return
      end
      subroutine grid
c=======================================================================
c  G E N E R A T E   T H E   F I N I T E   D I F F E R E N C E   G R I D
c=======================================================================
      parameter    (jmx=201)
      common/bcom/ climit(jmx),deta(jmx),dudxi(jmx),dxi(jmx),en(jmx),
     *             eps(jmx),et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),
     *             v(jmx),vextra(jmx),wt(jmx),wtm(jmx),xi(jmx),y1(jmx),
     *             y2(jmx),y3(jmx)
c-----------------------------------------------------------------------
      common/extra/    acon   ,  akappa ,  cflmax ,  dti    ,  eone   ,
     *                 epsin  ,  epsone ,  error  ,  etamax ,  etamin ,
     *                 etin   ,  jend   ,  nstep  ,  pitld  ,  umax   ,
     *                 umin   ,  u1     ,  uone   ,  vcon   ,  wcon   ,
     *                 wtin   ,  wtmax  ,  wtone  ,  ximax  ,  ximin  ,
     *                 xk
      common/ndata/    jmax   ,  maxn   ,  nedit  ,  nfreq
c=======================================================================
      if(abs(xk-1.).lt.1.e-6) then
        delta=(ximax-ximin)/jend
      else
        delta=(ximax-ximin)*(xk-1.)/(xk**jend-1.)
      endif
      xi(1)=ximin
      dxi(1)=delta
      do 10 j=2,jmax
        dxi(j)=xk*dxi(j-1)
        xi(j)=xi(j-1)+dxi(j-1)
 10   continue
c=======================================
c   Compute differencing coefficients
c=======================================
      do 20 j=2,jend
        denom=dxi(j)+dxi(j-1)
        y1(j)=2./(dxi(j)*denom)
        y2(j)=1./denom
        y3(j)=2./(dxi(j-1)*denom)
 20   continue
      return
      end
      subroutine namin(icpu)
c=======================================================================
c            COORDINATE READING INPUT DATA FILE 'DEFECT.DAT'
c=======================================================================
      parameter (niu=4,nn=4,nxd=1)
      character a*12,newfil*7,pform*9
      common/data/  xd(nxd)
      common/iodef/ iu(niu)
      common/ndata/ n(nn)
c=======================================================================
c----------------------------------
c         Read I/O units
c----------------------------------
      call namsys(icpu,iin,iv,msdos,newfil,pform)
      open(iin,file='defect.dat',status='old')
      iu(1)=2
      iu(2)=7
      io=iu(1)
      open(io,file='defect.prt',form=pform,status=newfil)
      write(io,200)
      do 10 i=3,niu
        read(iin,100) a,iu(i)
        write(io,100) a,iu(i) 
 10   continue
c----------------------------------
c    Read floating-point input
c----------------------------------
      do 20 i=1,nxd
        read(iin,300) a,xd(i)
        write(io,300) a,xd(i) 
 20   continue
c----------------------------------
c      Read fixed-point input
c----------------------------------
      do 30 i=1,nn
        read(iin,100) a,n(i)
        write(io,100) a,n(i) 
 30   continue
      write(io,400)
c----------------------------------
c   Open plotting-data disk files
c----------------------------------
      open(iu(2),file='plotd.dat',form=pform,status=newfil)
      return
c=======================================================================
 100  format(1x,a11,i5)
 200  format(1x,'Program DEFECT input data:'/,
     *       1x,'__________________________'/)
 300  format(1x,a12,1pe13.6)
 400  format(1x,'__________________________'//)
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
c  SVS Fortran (680x0 and 80386)
c----------------------------------
      if(icpu.eq.0) then
        msdos=0
        pform='printer'
        newfil='new'
c----------------------------------
c  Microsoft Fortran (8088/80x86)
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
      subroutine start(iflag)
c=======================================================================
c          S E T   U P   I N I T I A L   P R O F I L E S
c=======================================================================
c    NOTE:  We call START with IFLAG = 0 to determine AKAPPA for
c           Subroutine GRID...otherwise do all computations.
c-----------------------------------------------------------------------
      parameter    (jmx=201)
      common/bcom/ climit(jmx),deta(jmx),dudxi(jmx),dxi(jmx),en(jmx),
     *             eps(jmx),et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),
     *             v(jmx),vextra(jmx),wt(jmx),wtm(jmx),xi(jmx),y1(jmx),
     *             y2(jmx),y3(jmx)
c-----------------------------------------------------------------------
      common/data/     betat
      common/extra/    acon   ,  akappa ,  cflmax ,  dti    ,  eone   ,
     *                 epsin  ,  epsone ,  error  ,  etamax ,  etamin ,
     *                 etin   ,  jend   ,  nstep  ,  pitld  ,  umax   ,
     *                 umin   ,  u1     ,  uone   ,  vcon   ,  wcon   ,
     *                 wtin   ,  wtmax  ,  wtone  ,  ximax  ,  ximin  ,
     *                 xk
      common/iodef/    iunit1 ,  iunit2 ,  model  ,  kerng
      common/ndata/    jmax   ,  maxn   ,  nedit  ,  nfreq
      common/oneeqn/   acmax  ,  biga   ,  bigb   ,  bigc   ,  cw1    ,
     *                 cw2    ,  cw3    ,  en1    ,  fwo
      common/rng/      ceps2p ,  phi    ,  rbeta  ,  rlamda
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c=======================================================================
      if(model.eq.2.and.kerng.ne.0) then
        betas=.085
      else
        betas=.09
      endif
      alphas=sqrt(betas)
      etin=1.e-8
c=======================================
c     k-omega model coefficients
c=======================================
      if(model.eq.0) then
        alpha=.52
        beta=.0708
        sigma=.5
        sigmad=.125
        sigmas=.6
        akappa=sqrt(alphas*(beta/betas-alpha)/sigma)
        clim=7./8.
        wcon=1.+4.*betat
        c1=.5*sigmas*akappa**2/alphas
        c2=beta/(alpha*betas)
        e1=2.*betat/(akappa*(c1-1.))
        w1=(1.+.5*sigmad/sigmas)*c1*e1/(1.-c2)
        u1=(c2+.5*sigmad/sigmas)*c1*e1/(1.-c2)
        if(betat.ge.0.) then
          wtin=alphas*wcon/beta
        else
          wtin=alphas/beta
        endif
        epsin=1.e-4
        etin=epsin*wtin
c=======================================
c Standard k-epsilon model coefficients
c=======================================
      elseif(model.eq.2) then
        if(kerng.eq.0) then
          ceps1=1.44
          ceps2=1.92
          alpha=ceps1
          beta=betas*ceps2
          sigma=10./13.
          sigmas=1.
          akappa=sqrt(alphas*(ceps2-ceps1)/sigma)
          wcon=1.+6.*betat
          c1=.5*sigmas*akappa**2/alphas
          c2=1.+2.*c1
          e1=2.*betat/(akappa*(c1-1.))
          com=.5*e1/(ceps1-ceps2)
          w1=com*(c2*ceps1-ceps2)
          u1=com*(c2*ceps2-ceps1)
          if(betat.ge.0.) then
            wtin=.01*wcon/(alphas*ceps2)
          else
            wtin=.01/(alphas*ceps2)
          endif
          wtin=wtin**2
          epsin=1.e-4
          etin=sqrt(epsin*wtin)
c=======================================
c   RNG k-epsilon model coefficients
c=======================================
        else
          ceps1=1.42
          ceps2=1.68
          alpha=ceps1
          beta=betas*ceps2
          sigma=1.39
          sigmas=1.39
          rbeta=.012
          rlamda=4.38
          ceps2p=(betas-alphas/rlamda)/(alphas*betas+rbeta)
          phi=3.-1./(rlamda*alphas-1.)-3.*rbeta/(alphas*betas+rbeta)
          akappa=sqrt(alphas*(ceps2+ceps2p-ceps1)/sigma)
          wcon=1.+6.*betat
          c1=.5*sigmas*akappa**2/alphas
          e1=2.*betat/(akappa*(c1-1.))
          w1=.5*(sigma-ceps1*sigmas+.5*ceps2p*phi*sigmas)*e1/sigma
          u1=w1-(1.+c1)*e1
          if(betat.ge.0.) then
            wtin=.01*wcon/(alphas*ceps2)
          else
            wtin=.01/(alphas*ceps2)
          endif
          wtin=wtin**2
          epsin=1.e-4
          etin=sqrt(epsin*wtin)
        endif
c=======================================
c   Baldwin-Barth model coefficients
c=======================================
      elseif(model.eq.3) then
        akappa=.41
        ceps1=1.2
        ceps2=2.
        sigma=(ceps2-ceps1)*alphas/akappa**2
        wcon=1.+2.*betat
        u1=-betat/akappa
        en1=u1
        biga=3.*sigma
        bigb=-2.*sigma
        bigc=sigma*akappa**2
        epsin=.002
c=======================================
c  Spalart-Allmaras model coefficients
c=======================================
      elseif(model.eq.4) then
        akappa=.41
        cb1=.1355
        cb2=.622
        sigma=1.5
        cw1=cb1/akappa**2+(1.+cb2)*sigma
        cw2=.3
        cw3=2.
        fwo=(1.+cw3**6)**(1./6.)
        wcon=1.+2.*betat
        c1=cw1*(1.+5.*cw2)*cw3**6/(1.+cw3**6)*akappa**2
        c2=2.*sigma*akappa**2
        u1=2.*(-cb1+c2*(2.+cb2)-c1)/(4.*cb1-c2+2.*c1)*betat/akappa
        en1=-u1-2.*betat/akappa
        biga=(1.-cb2)*sigma
        bigb=cb2*sigma
        bigc=cb1
        epsin=0.
      endif
c=======================================================================
c          C O M P U T E   I N I T I A L   P R O F I L E S
c-----------------------------------------------------------------------
c The following profiles are correlations of converged k-omega solutions
c=======================================================================
      if(iflag.ne.0) then
c-----------------------
c       Velocity
c-----------------------
        u0=14.+16.*(betat/9.)
        ulam=.0083+.003*betat
        do 10 j=1,jend
          if(xi(j).le.-16.) then
            u(j)=u0-(xi(j)+16.)
          else
            u(j)=u0*exp(-ulam*(xi(j)+16.))
          endif
 10     continue
        u(jmax)=0.
c-----------------------
c    Eddy viscosity
c-----------------------
        epsmax=.015+.01*exp(-betat)
        epsmin=akappa*etamin
        xi0=8.*exp(-betat/3.)-11.5
        xi1=2.*xi0-ximin
        do 20 j=1,jend
          if(xi(j).le.xi0) then
            eps(j)=epsmin+(epsmax-epsmin)*(xi(j)-ximin)/(xi0-ximin)
          elseif(xi(j).le.xi1) then
            eps(j)=epsmax+(epsmax-epsmin)*(xi0-xi(j))/(xi0-ximin)
          else
            eps(j)=epsmin+(epsin-epsmin)*(xi(j)-xi1)/(ximax-xi1)
          endif
 20     continue
        eps(jmax)=epsmax
c-----------------------
c  Integral constraint
c-----------------------
        uint=0.
        do 30 j=1,jend
          uint=uint+.5*(eps(j)*u(j)+eps(j+1)*u(j+1))*dxi(j)
 30     continue
        renorm=1./uint
        do 40 j=1,jmax
          u(j)=renorm*u(j)
 40     continue
        uone=u(1)
c-----------------------
c  Turbulence energy
c-----------------------
        elog=dlog(etamin)
        eloge=etamin*elog
        eone=1.+e1*eloge
        etmax=1+.5*betat
        elam=.18+.57*(betat/9.)
        do 50 j=1,jend
          if(xi(j).le.-10.) then
            et(j)=eone+(etmax-eone)*(xi(j)-ximin)/(10.+ximin)
          else
            et(j)=etmax*exp(-elam*(xi(j)+10.))
          endif
          if(et(j).lt.etin) et(j)=etin
 50     continue
        et(jmax)=etin
c-----------------------
c   Dissipation rate
c-----------------------
        wtone=(1.+w1*eloge)/(akappa*etamin)
        do 60 j=1,jend
          wt(j)=et(j)/eps(j)
          if(model.eq.2) wt(j)=et(j)*wt(j)
          if(wt(j).lt.wtin) wt(j)=wtin
 60     continue
        wt(jmax)=wtin
        do 70 j=1,jmax
          um(j)=u(j)
          etm(j)=et(j)
          wtm(j)=wt(j)
          en(j)=eps(j)
 70     continue
c-----------------------
c    Eddy viscosity
c       at etamin
c-----------------------
        if(model.eq.3.or.model.eq.4) then
          epsone=akappa*etamin*(1.+en1*eloge)
        else
          epsone=eone/wtone
          if(model.eq.2) epsone=eone*epsone
        endif
      endif
      return
      end
      subroutine tmestp
c=======================================================================
c                 C O M P U T E   T I M E S T E P
c=======================================================================
      parameter    (jmx=201)
      common/bcom/ climit(jmx),deta(jmx),dudxi(jmx),dxi(jmx),en(jmx),
     *             eps(jmx),et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),
     *             v(jmx),vextra(jmx),wt(jmx),wtm(jmx),xi(jmx),y1(jmx),
     *             y2(jmx),y3(jmx)
c-----------------------------------------------------------------------
      common/extra/    acon   ,  akappa ,  cflmax ,  dti    ,  eone   ,
     *                 epsin  ,  epsone ,  error  ,  etamax ,  etamin ,
     *                 etin   ,  jend   ,  nstep  ,  pitld  ,  umax   ,
     *                 umin   ,  u1     ,  uone   ,  vcon   ,  wcon   ,
     *                 wtin   ,  wtmax  ,  wtone  ,  ximax  ,  ximin  ,
     *                 xk
c=======================================================================
c    Ramp the CFL number up during 
c       the first 500 timesteps
c=======================================
      if(nstep.lt.100) then
        cfl=.1
      elseif(nstep.lt.500) then
        cfl=.1+(cflmax-.1)*(nstep-100.)/400.
      else
        cfl=cflmax
      endif
c=======================================
c     Now determine the timestep
c=======================================
      dtimax=0.
      do 10 j=2,jend
        dti=abs(v(j))/dxi(j-1)
        dtimax=dmax1(dtimax,dti)
        dti=2./dxi(j-1)**2
        dtimax=dmax1(dtimax,dti)
 10   continue
      dtimax=dtimax/cfl
      dti=dtimax
      return
      end
      subroutine tri(a,b,c,d,il,iu,ibcl,dl,ibcu,du)
c=======================================================================
c         T R I - D I A G O N A L   M A T R I X   S O L V E R
c=======================================================================
      dimension a(*),b(*),c(*),d(*)
c=======================================================================
c  Set the lower boundary condition
c==================================
      if(ibcl.eq.1) then
c=======================================
c  Neumann boundary condition for ibcl=1
c=======================================
        a(il)=-dl
        b(il)=1.
c===============================================
c  Otherwise assume Dirichlet boundary condition
c===============================================
      else
        a(il)=dl
        b(il)=0.
      endif
c===============
c  Forward sweep
c===============
      ilp=il+1
      ium=iu-1
      do 10 i=ilp,ium
        bs=b(i)+a(i)*b(i-1)
        a(i)=(d(i)-a(i)*a(i-1))/bs
        b(i)=-c(i)/bs
 10   continue
c==================================
c  Set the upper boundary condition
c==================================
      if(ibcu.eq.1) then
c=======================================
c  Neumann boundary condition for ibcu=1
c=======================================
        d(iu)=(du+a(ium))/(1.-b(ium))
c===============================================
c  Otherwise assume Dirichlet boundary condition
c===============================================
      else
        d(iu)=du
      endif
c===================================================
c  Back substitution to determine the final solution
c===================================================
      it=il+iu
      do 20 ii=ilp,iu
        i=it-ii
        d(i)=a(i)+b(i)*d(i+1)
 20   continue
      return
      end
