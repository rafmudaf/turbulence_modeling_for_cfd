      program jet
c=======================================================================
c                 I N C O M P R E S S I B L E   J E T
c
c    Copyright (c) 2006 DCW Industries, Inc.  All rights reserved.
c-----------------------------------------------------------------------
c            MODEL = 0.........................k-omega Model
c            MODEL = 1....................Stress-omega Model
c            MODEL = 2, KERNG = 0...Standard k-epsilon Model
c            MODEL = 2, KERNG > 0........RNG k-epsilon Model
c            MODEL = 4................Spalart-Allmaras Model
c            MODEL = 9...................Mixing-Length Model
c            MODEL = 99...................User-Defined Model
c-----------------------------------------------------------------------
c                        JAXI = -1...Radial Jet
c                        JAXI =  0...Plane  Jet
c                        JAXI =  1...Round  Jet
c-----------------------------------------------------------------------
c                IPOPE = 0.....No Round Jet Correction
c                IPOPE = 1...Pope Round Jet Correction
c=======================================================================
      parameter (jmx=201)
      common/bcom/ deta(jmx),dudxi(jmx),dxi(jmx),ef(jmx),eps(jmx),
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),v(jmx),
     *             wt(jmx),wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx),
     *             z1(jmx),z3(jmx),zdiss(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx)
      common/ecom/ climit(jmx),fbeta(jmx),uoum(jmx),vextra(jmx)
c-----------------------------------------------------------------------
      common/data/     etin   ,  wtin
      common/extra/    ceps3  ,  cfl    ,  dti    ,  encon  ,  error  ,
     *                 fcon   ,  gcon   ,  jbar   ,  jend   ,  nstep  ,
     *                 pow    ,  scon   ,  tcon   ,  umax   ,  umin1  ,
     *                 umin2  ,  wtmax  ,  ximax  ,  xk     ,  zcon
      common/iodef/    iunit1 ,  iunit2 ,  model  ,  kerng
      common/ndata/    ipope  ,  jaxi   ,  jmax   ,  maxn   ,  nedit  ,
     *                 nfreq
      common/mixlen/   amix
      common/oneeqn/   aa1    ,  aa2    ,  aa3
      common/stress/   alphah ,  betah  ,  c1     ,  gammah ,  onemab
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c-----------------------------------------------------------------------
      include 'cpuid.'
      data pi/3.14159265/
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
      error=1.e-6
      jbar=iabs(jaxi)
      encon=.5*(jbar-1.)
      fcon=1./(2.-jbar)
      scon=(3.+jbar)/2.
      tcon=1.+jbar
      zcon=(5.+3.*jbar)/2.
      xk=1.
      umin1=1.e-2
      umin2=1.e-4
      wtmax=1000.
      if(jaxi.lt.0) then
        gcon=4.*pi
        ximax=50.
        umax=3.
        if(model.eq.1) then
          cfl=1.5
        elseif(model.eq.2.and.kerng.eq.0) then
          cfl=1.
        else
          cfl=5.
        endif
      elseif(jaxi.eq.0) then
        gcon=2.
        ximax=60.
        umax=3.
        if(model.eq.1) then
          cfl=2.
        elseif(model.eq.2.and.kerng.eq.0) then
          cfl=4.
        else
          cfl=5.
        endif
      else
        gcon=2.*pi
        ximax=25.
        umax=12.
        if(model.eq.1) then
          cfl=.5
        else
          cfl=2.
        endif
      endif
c=======================================
c         Generate the grid
c=======================================
      call grid
c=======================================
c        Set initial conditions
c=======================================
      call start
      nstep=0
c=======================================
c   Identify flow and turbulence model
c=======================================
      if(jaxi.lt.0) then
        write(*,110)
      elseif(jaxi.eq.0) then
        write(*,120)
      else
        write(*,130)
      endif
      if(model.eq.0) then
        write(*,140)
      elseif(model.eq.1) then
        write(*,150)
      elseif(model.eq.2) then
        if(kerng.eq.0.) then
          write(*,160)
        else
          write(*,165)
        endif
      elseif(model.eq.4) then
        write(*,170)
      elseif(model.eq.9) then
        write(*,180)
      elseif(model.eq.99) then
        write(*,190)
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
      do 30 j=1,jmax
        err=abs(u(j)-um(j))
        if(err.ge.erro) then
          erro=err
        endif
 30   continue
c=======================================
c          Converged solution
c=======================================
      if(erro.le.error) then
        call geteta
        call grow(dldx,jmax)
        write(iunit1,200) nstep,erro,dldx
        write(*,200) nstep,erro,dldx
        write(*,210)
c=======================================
c       Un-converged solution
c=======================================
      else
        if(mod(nstep,nfreq).eq.0) then
          call geteta
          call grow(dldx,jmax)
          write(iunit1,200) nstep,erro,dldx
          write(*,200) nstep,erro,dldx
        endif
        if(nstep.ne.maxn) then
c=======================================
c          Print profiles
c=======================================
          if(mod(nstep,nedit).eq.0) call edit
c=======================================
c      Prepare for next timestep
c=======================================
          do 40 j=1,jmax
            etm(j)=et(j)
            um(j)=u(j)
            wtm(j)=wt(j)
            if(model.eq.1) then
              sigxm(j)=sigx(j)
              sigym(j)=sigy(j)
              taum(j)=tau(j)
            endif
 40       continue
          go to 10
c=======================================================================
c                   END OF MAIN INTEGRATION LOOP
c=======================================================================
        else
          write(iunit1,220) nstep
          write(*,220) nstep
        endif
      endif
      call edit
      close(iunit1)
c=======================================
c    Write plotting data disk files
c=======================================
      write(iunit2,230) jmax,jaxi
      write(iunit2,240) (eta(j),uoum(j),j=1,jmax)
      endfile(iunit2)
      close(iunit2)
c=======================================
c    Pause to keep window visible after
c    the run completes...
c=======================================
      pause
c=======================================================================
 100  format(' Convergence history:'/)
 110  format(' RADIAL JET')
 120  format(' PLANE JET')
 130  format(' ROUND JET')
 140  format(' k-omega Model')
 150  format(' Stress-omega Model')
 160  format(' Standard k-epsilon Model')
 165  format(' RNG k-epsilon Model')
 170  format(' Spalart-Allmaras Model')
 180  format(' Mixing-Length Model')
 190  format(' User-Defined Model')
 200  format(1x,'At Iteration ',i5,'...Max Error = ',1pe12.5,
     *          '  and dldx =',e10.3)
 210  format(/19x,' ***** SOLUTION CONVERGED *****'/)
 220  format(/8x,'Computation terminated for max iteration',i5/)
 230  format(1x,2i4)
 240  format(1p2e14.6)
c=======================================================================
      end
      subroutine calcs
c=======================================================================
c C A L C U L A T E   Q U A N T I T I E S   N E E D E D   I N   E D D Y
c=======================================================================
      parameter (jmx=201)
      common/bcom/ deta(jmx),dudxi(jmx),dxi(jmx),ef(jmx),eps(jmx),
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),v(jmx),
     *             wt(jmx),wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx),
     *             z1(jmx),z3(jmx),zdiss(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx)
      common/ecom/ climit(jmx),fbeta(jmx),uoum(jmx),vextra(jmx)
c-----------------------------------------------------------------------
      common/data/     etin   ,  wtin
      common/extra/    ceps3  ,  cfl    ,  dti    ,  encon  ,  error  ,
     *                 fcon   ,  gcon   ,  jbar   ,  jend   ,  nstep  ,
     *                 pow    ,  scon   ,  tcon   ,  umax   ,  umin1  ,
     *                 umin2  ,  wtmax  ,  ximax  ,  xk     ,  zcon
      common/iodef/    iunit1 ,  iunit2 ,  model  ,  kerng
      common/ndata/    ipope  ,  jaxi   ,  jmax   ,  maxn   ,  nedit  ,
     *                 nfreq
      common/mixlen/   amix
      common/stress/   alphah ,  betah  ,  c1     ,  gammah ,  onemab
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c=======================================================================
c  Streamfunction & 'vertical' velocity
c=======================================
      if(jaxi.ne.1) then
        dum=u(1)*eps(1)
      else
        dum=0.
      endif
      ef(1)=0.
      do 10 j=2,jmax
        dup=u(j)*eps(j)
        if(jaxi.eq.1) dup=dup*eta(j)
        ef(j)=ef(j-1)+.5*(dum+dup)*dxi(j-1)
        dum=dup
 10   continue
      do 20 j=2,jmax
        v(j)=-fcon*ef(j)
        vextra(j)=0.
        if(jaxi.eq.1) v(j)=v(j)/eta(j)
 20   continue
      v(1)=0.
      vextra(1)=0.
c=======================================
c  Eddy viscosity and velocity gradient
c=======================================
      dudxi(1)=0.
      dudxi(jmax)=0.
      do 30 j=1,jmax
        climit(j)=1.
        if(model.eq.0) then
          eps(j)=et(j)/wt(j)
          if(j.ne.1.and.j.ne.jmax) dudxi(j)=y2(j)*(u(j+1)-u(j-1))
          if(abs(clim*dudxi(j)/alphas).gt.et(j)) then
            climit(j)=alphas*et(j)/abs(clim*dudxi(j))
          endif
        elseif(model.eq.1) then
          eps(j)=et(j)/wt(j)
          if(j.ne.1.and.j.ne.jmax) dudxi(j)=y2(j)*(u(j+1)-u(j-1))
        elseif(model.eq.2) then
          eps(j)=betas*et(j)**2/wt(j)
          dudxi(j)=v(j)*u(j)
        elseif(model.eq.4) then
          eps(j)=et(j)
          dudxi(j)=v(j)*u(j)
        elseif(model.eq.9) then
          if(j.ne.1.and.j.ne.jmax) dudxi(j)=y2(j)*(u(j+1)-u(j-1))
          eps(j)=amix*sqrt(abs(dudxi(j)))
        endif
 30   continue
c=======================================
c  Pope's round-jet correction term
c---------------------------------------
c  NOTE: The variable ef is eta*F(eta)
c        for the round jet, so that
c        ef/eta**2 is F(eta)/eta
c=======================================
      if(model.eq.2.and.jaxi.ne.0.and.ipope.ne.0) then
        do 50 j=2,jend
          if(jaxi.eq.1) then
            ubar=um(j)-ef(j)/eta(j)**2
          else
            ubar=um(j)
          endif
          chi=.25*dudxi(j)**2*ubar
          zdiss(j)=ceps3*chi/(betas*wtm(j))
 50     continue
      endif
c=======================================
c      k-omega model computations
c=======================================
      if(model.eq.0.or.model.eq.1) then
        fbeta(1)=1.
        do 60 j=2,jend
          vextra(j)=0.
          dkdxi=(1./et(j+1)-1./et(j-1))/(dxi(j-1)+dxi(j))
          dwdxi=dlog(wt(j+1)/wt(j-1))/(dxi(j-1)+dxi(j))
          chik=-dkdxi*dwdxi
          if(chik.gt.0.) vextra(j)=1000.
          if(jaxi.ne.0) then
            if(jaxi.eq.1) then
              ubar=um(j)-ef(j)/eta(j)**2
            else
              ubar=um(j)
            endif
            chiw=.25*dudxi(j)**2*ubar/(betas**3*etm(j)**2*wtm(j))
            fbeta(j)=(1.+85.*abs(chiw))/(1.+100.*abs(chiw))
          endif
 60     continue
        fbeta(jmax)=fbeta(jend)
      endif
c=======================================
c     k-epsilon model computations
c=======================================
      if(model.eq.2.and.kerng.ne.0) then
        fbeta(1)=1.
        do 70 j=2,jend
          alamm=abs(dudxi(j))/(betas*et(j))
          factor=alamm**3*(1.-alamm/4.38)/(1.+.012*alamm**3)
          fbeta(j)=1.+betas/beta*factor
 70     continue
        fbeta(jmax)=fbeta(jend)
      endif
      return
      end
      subroutine eddy
c=======================================================================
c      A D V A N C E   T H E   S O L U T I O N   I N   T I M E
c=======================================================================
      parameter (jmx=201)
      dimension    a(jmx),b(jmx),c(jmx)
      common/bcom/ deta(jmx),dudxi(jmx),dxi(jmx),ef(jmx),eps(jmx),
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),v(jmx),
     *             wt(jmx),wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx),
     *             z1(jmx),z3(jmx),zdiss(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx)
      common/ecom/ climit(jmx),fbeta(jmx),uoum(jmx),vextra(jmx)
c-----------------------------------------------------------------------
      common/data/     etin   ,  wtin
      common/extra/    ceps3  ,  cfl    ,  dti    ,  encon  ,  error  ,
     *                 fcon   ,  gcon   ,  jbar   ,  jend   ,  nstep  ,
     *                 pow    ,  scon   ,  tcon   ,  umax   ,  umin1  ,
     *                 umin2  ,  wtmax  ,  ximax  ,  xk     ,  zcon
      common/iodef/    iunit1 ,  iunit2 ,  model  ,  kerng
      common/mixlen/   amix
      common/ndata/    ipope  ,  jaxi   ,  jmax   ,  maxn   ,  nedit  ,
     *                 nfreq
      common/oneeqn/   aa1    ,  aa2    ,  aa3
      common/stress/   alphah ,  betah  ,  c1     ,  gammah ,  onemab
      common/stres2/   cf1    ,  cf2    ,  cf3    ,  cfx    ,  cfy
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c=======================================================================
      if(model.eq.4) then
        go to 50
      elseif(model.eq.9) then
        go to 80
      endif
c=======================================
c          Dissipation rate
c=======================================
      do 10 j=2,jend
        if(model.eq.0.or.model.eq.1) then
          prodw=scon*eps(j)*um(j)-beta*fbeta(j)*etm(j)
          if(etm(j).ne.0.) prodw=prodw+alpha*v(j)*um(j)*dudxi(j)/etm(j)
          cross=etm(j)**pow
          prodw=prodw*cross
          if(vextra(j).gt.1.) then
            etmh=(.5*(etm(j-1)+etm(j)))**pow
            etph=(.5*(etm(j+1)+etm(j)))**pow
            a(j)=sigma*etmh*y3(j)+v(j)*y2(j)*cross
            b(j)=-sigma*(etph*y1(j)+etmh*y3(j))-2.*dti*cross
            c(j)=sigma*etph*y1(j)-v(j)*y2(j)*cross
          else
            a(j)=(sigma*y3(j)+v(j)*y2(j))*cross
            b(j)=(-sigma*(y1(j)+y3(j))-2.*dti)*cross
            c(j)=(sigma*y1(j)-v(j)*y2(j))*cross
          endif
          wt(j)=-a(j)*wtm(j-1)-(b(j)+4.*dti*cross)*wtm(j)-c(j)*wtm(j+1)
        elseif(model.eq.2) then
          prodw=zcon*eps(j)*um(j)-beta*betas*fbeta(j)*etm(j)+zdiss(j)
          if(etm(j).ne.0.) prodw=prodw+alpha*v(j)*um(j)*dudxi(j)/etm(j)
          a(j)=sigma*y3(j)+v(j)*y2(j)
          b(j)=-sigma*(y1(j)+y3(j))-2.*dti
          c(j)=sigma*y1(j)-v(j)*y2(j)
          wt(j)=-a(j)*wtm(j-1)-(b(j)+4.*dti)*wtm(j)-c(j)*wtm(j+1)
        endif
        if(prodw.gt.0.) then
          wt(j)=wt(j)-2.*prodw*wtm(j)
        else
          b(j)=b(j)+2.*prodw
        endif
 10   continue
      call tri(a,b,c,wt,1,jmax,1,0.,0,wtin)
      do 20 j=2,jend
        if(wt(j).lt.0.) wt(j)=.5*(wtm(j-1)+wtm(j+1))
        if(wt(j).gt.wtmax) wt(j)=wtmax
 20   continue
      if(wt(1).lt.0.) wt(1)=wt(2)
c=======================================
c      Turbulence kinetic energy
c=======================================
      do 30 j=2,jend
        prode=tcon*eps(j)*um(j)-betas*etm(j)
        if(etm(j).ne.0.) then
          prode=prode+v(j)*um(j)*dudxi(j)/etm(j)
        endif
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
      call tri(a,b,c,et,1,jmax,1,0.,0,etin)
      do 40 j=2,jend
        if(et(j).lt.0.) et(j)=.5*(etm(j-1)+etm(j+1))
 40   continue
      if(et(1).lt.0.) et(1)=et(2)
      go to 80
c=======================================
c  One-equation model: eddy viscosity
c=======================================
 50   epsm=.5*aa3*(etm(1)+etm(2))
      do 60 j=2,jend
        epsp=.5*aa3*(etm(j+1)+etm(j))
        adiff=aa2-epsm/etm(j)
        cdiff=aa2-epsp/etm(j)
        prode=encon*etm(j)+aa1*abs(dudxi(j))
        a(j)=adiff*y3(j)+v(j)*y2(j)
        b(j)=-cdiff*y1(j)-adiff*y3(j)-2.*dti
        c(j)=cdiff*y1(j)-v(j)*y2(j)
        et(j)=-a(j)*etm(j-1)-(b(j)+4.*dti)*etm(j)-c(j)*etm(j+1)
        et(j)=et(j)-2.*prode*etm(j)
        epsm=epsp
 60   continue
      call tri(a,b,c,et,1,jmax,1,0.,0,etin)
      do 70 j=2,jend
        if(et(j).lt.0.) et(j)=.5*(etm(j-1)+etm(j+1))
 70   continue
      if(et(1).lt.0.) et(1)=et(2)
c=======================================
c              Momentum
c=======================================
 80   call calcs
      if(model.eq.1) then
        do 85 j=2,jend
          a(j)=y3(j)+v(j)*y2(j)
          b(j)=-y1(j)-y3(j)-2.*dti
          c(j)=y1(j)-v(j)*y2(j)
          u(j)=-a(j)*um(j-1)-(b(j)+4.*dti)*um(j)-c(j)*um(j+1)
     *         +2.*(eta(j+1)**jbar*dtau(j+1)-eta(j-1)**jbar*dtau(j-1))
     *         /(eta(j)**jbar*(dxi(j-1)+dxi(j)))
     *         -2.*fcon*eps(j)*um(j)**2
 85     continue
        call tri(a,b,c,u,1,jmax,1,0.,0,0.)
      else
        uint=0.
        dum=0.
        do 90 j=2,jmax
          dup=v(j)/climit(j)
          uint=uint+.5*(dum+dup)*dxi(j-1)
          if(uint.gt.-20.) then
            u(j)=u(1)*exp(uint)
          else
            u(j)=0.
          endif
          dum=dup
 90     continue
      endif
c=======================================
c      Impose integral constraint
c=======================================
      uint=0.
      if(jaxi.ne.1) then
        dum=eps(1)*u(1)**2
      else
        dum=0.
      endif
      do 100 j=2,jmax
        dup=eps(j)*u(j)**2
        if(jaxi.eq.1) dup=dup*eta(j)
        uint=uint+.5*(dum+dup)*dxi(j-1)
        dum=dup
 100  continue
      uint=gcon*uint
      uint=sqrt(uint)
      do 110 j=1,jmax
        u(j)=u(j)/uint
        taum(j)=u(j)*v(j)
 110  continue
c-----------------------------------------------------------------------
c                    Reynolds-stress equations
c-----------------------------------------------------------------------
      if(model.eq.1) then
        call calcs
c=======================================
c             Shear Stress
c=======================================
        do 120 j=2,jend
          a(j)=sigmas*y3(j)+v(j)*y2(j)
          b(j)=-sigmas*(y1(j)+y3(j))-betas*c1*et(j)-2.*dti
          c(j)=sigmas*y1(j)-v(j)*y2(j)
          tau(j)=-a(j)*taum(j-1)-(b(j)+4.*dti)*taum(j)-c(j)*taum(j+1)
     *           -2.*(cf1*sigxm(j)+cf2*sigym(j)+cf3*et(j))*dudxi(j)
 120    continue
        call tri(a,b,c,tau,1,jmax,1,0.,0,0.)
c=======================================
c          Normal Stresses
c=======================================
        do 140 j=2,jend
          a(j)=sigmas*y3(j)+v(j)*y2(j)
          b(j)=-sigmas*(y1(j)+y3(j))-betas*c1*et(j)-2.*dti
          c(j)=sigmas*y1(j)-v(j)*y2(j)
          sigx(j)=-a(j)*sigxm(j-1)-(b(j)+4.*dti)*sigxm(j)
     *            -c(j)*sigxm(j+1)-2.*cfx*tau(j)*dudxi(j)
 140    continue
        call tri(a,b,c,sigx,1,jmax,1,0.,0,0.)
        do 150 j=2,jend
          a(j)=sigmas*y3(j)+v(j)*y2(j)
          b(j)=-sigmas*(y1(j)+y3(j))-betas*c1*et(j)-2.*dti
          c(j)=sigmas*y1(j)-v(j)*y2(j)
          sigy(j)=-a(j)*sigym(j-1)-(b(j)+4.*dti)*sigym(j)
     *            -c(j)*sigym(j+1)-2.*cfy*tau(j)*dudxi(j)
 150    continue
        call tri(a,b,c,sigy,1,jmax,1,0.,0,0.)
c=======================================
c        Realizability checks
c=======================================
        relax=.05
        relaxm=1.-relax
        do 160 j=2,jend
          dtau0=tau(j)-dudxi(j)
          dtau(j)=relax*dtau0+relaxm*dtau(j)
          term=.66666667*et(j)
          upup=sigx(j)+term
          if(upup.lt.0.) sigx(j)=-term
          if(upup.gt.et(j)) sigx(j)=.5*term
          vpvp=sigy(j)+term
          if(vpvp.lt.0.) sigy(j)=-term
          if(vpvp.gt.et(j)) sigy(j)=.5*term
 160    continue
      endif
c-----------------------------------------------------------------------
      return
      end
      subroutine edit
c=======================================================================
c                   P R I N T I N G   R O U T I N E
c=======================================================================
      parameter (jmx=201)
      dimension    upup(jmx),vpvp(jmx)
      common/bcom/ deta(jmx),dudxi(jmx),dxi(jmx),ef(jmx),eps(jmx),
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),v(jmx),
     *             wt(jmx),wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx),
     *             z1(jmx),z3(jmx),zdiss(jmx)
      common/bdkw/ cnvk(jmx),cnvw(jmx),dflk(jmx),dflw(jmx),
     *             dftk(jmx),dftw(jmx),disk(jmx),disw(jmx),
     *             esk(jmx),esw(jmx),prdk(jmx),prdw(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx)
      common/ecom/ climit(jmx),fbeta(jmx),uoum(jmx),vextra(jmx)
      common/fcom/ dissk(jmx),el(jmx),prodk(jmx),wtprnt(jmx)
c-----------------------------------------------------------------------
      common/data/     etin   ,  wtin
      common/extra/    ceps3  ,  cfl    ,  dti    ,  encon  ,  error  ,
     *                 fcon   ,  gcon   ,  jbar   ,  jend   ,  nstep  ,
     *                 pow    ,  scon   ,  tcon   ,  umax   ,  umin1  ,
     *                 umin2  ,  wtmax  ,  ximax  ,  xk     ,  zcon
      common/iodef/    iunit1 ,  iunit2 ,  model  ,  kerng
      common/mixlen/   amix
      common/ndata/    ipope  ,  jaxi   ,  jmax   ,  maxn   ,  nedit  ,
     *                 nfreq
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c=======================================================================
      call grow(dldx,jmax)
      write(iunit1,100) nstep,dldx
      do 10 j=1,jmax
        uoum(j)=u(j)/u(1)
 10   continue
c-----------------------------
c      One-equation and
c    Mixing-length models
c-----------------------------
      if(model.eq.4.or.model.eq.9) then
        do 20 j=1,jmax
          tau(j)=dudxi(j)
 20     continue
        write(iunit1,110)
        write(iunit1,120) (j,eta(j),uoum(j),u(j),v(j),eps(j),tau(j),
     *                     j=1,jmax)
c-----------------------------
c k-omega and k-epsilon models
c-----------------------------
      elseif(model.eq.0.or.model.eq.2) then
        do 30 j=1,jmax
          tau(j)=dudxi(j)
          if(model.eq.0) then
            wtprnt(j)=wt(j)
            dissk(j)=betas*et(j)*wt(j)
          else
            wtprnt(j)=wt(j)/(betas*et(j))
            dissk(j)=wt(j)
          endif
          prodk(j)=tau(j)*dudxi(j)/eps(j)
          el(j)=sqrt(abs(et(j)))/wtprnt(j)
 30     continue
        write(iunit1,130)
        write(iunit1,140) (j,eta(j),uoum(j),u(j),v(j),et(j),wtprnt(j),
     *                 eps(j),tau(j),prodk(j),dissk(j),el(j),j=1,jmax)
c-----------------------------
c     Stress-omega model
c-----------------------------
      elseif(model.eq.1) then
        do 40 j=1,jmax
          term=.666666667*et(j)
          upup(j)=term+sigx(j)
          vpvp(j)=term+sigy(j)
          el(j)=sqrt(abs(et(j)))/wt(j)
 40     continue
        write(iunit1,150)
        write(iunit1,140) (j,eta(j),uoum(j),u(j),v(j),et(j),wt(j),
     *                eps(j),tau(j),upup(j),vpvp(j),el(j),j=1,jmax)
      endif
      return
c=======================================================================
 100  format(/1x,'nstep =',i6,5x,'dl/dx =',1pe11.4)
 110  format(/5x,'j',4x,'eta',6x,'u/Um',7x,'u',9x,'v',9x,'nut',7x,'tau')
 120  format(i6,1p6e10.3)
 130  format(/5x,'j',4x,'y/x',6x,'u/Um',7x,'u',9x,'v',9x,'k',8x,'omega',
     *        6x,'nut',7x,'tau',6x,'prod',6x,'diss',7x,'el')
 140  format(i6,1p11e10.3)
 150  format(/5x,'j',4x,'eta',6x,'u/Um',7x,'u',9x,'v',9x,'k',8x,'omega',
     *        6x,'nut',7x,'tau',7x,'uu',8x,'vv',8x,'el')
c=======================================================================
      end
      subroutine geteta
c=======================================================================
c                      C O M P U T E   E T A
c=======================================================================
      parameter (jmx=201)
      common/bcom/ deta(jmx),dudxi(jmx),dxi(jmx),ef(jmx),eps(jmx),
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),v(jmx),
     *             wt(jmx),wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx),
     *             z1(jmx),z3(jmx),zdiss(jmx)
      common/xxx/  climit(jmx)
c-----------------------------------------------------------------------
      common/data/     etin   ,  wtin
      common/extra/    ceps3  ,  cfl    ,  dti    ,  encon  ,  error  ,
     *                 fcon   ,  gcon   ,  jbar   ,  jend   ,  nstep  ,
     *                 pow    ,  scon   ,  tcon   ,  umax   ,  umin1  ,
     *                 umin2  ,  wtmax  ,  ximax  ,  xk     ,  zcon
      common/ndata/    ipope  ,  jaxi   ,  jmax   ,  maxn   ,  nedit  ,
     *                 nfreq
c=======================================================================
      eta(1)=0.
      eps1=eps(1)
      do 10 j=2,jend
        eps2=eps(j)
        eta(j)=eta(j-1)+.5*(eps1+eps2)*dxi(j-1)
        deta(j-1)=eta(j)-eta(j-1)
        eps1=eps2
 10   continue
      deta(jend)=deta(jend-1)
      eta(jmax)=eta(jend)+deta(jend)
      if(jaxi.eq.1) then
        do 20 j=2,jend
          y1(j)=.5*z1(j)*(eta(j+1)+eta(j))/eta(j)
          y3(j)=.5*z3(j)*(eta(j-1)+eta(j))/eta(j)
 20     continue
      endif
      return
      end
      subroutine grid
c=======================================================================
c  G E N E R A T E   T H E   F I N I T E   D I F F E R E N C E   G R I D
c=======================================================================
      parameter (jmx=201)
      common/bcom/ deta(jmx),dudxi(jmx),dxi(jmx),ef(jmx),eps(jmx),
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),v(jmx),
     *             wt(jmx),wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx),
     *             z1(jmx),z3(jmx),zdiss(jmx)
c-----------------------------------------------------------------------
      common/data/     etin   ,  wtin
      common/extra/    ceps3  ,  cfl    ,  dti    ,  encon  ,  error  ,
     *                 fcon   ,  gcon   ,  jbar   ,  jend   ,  nstep  ,
     *                 pow    ,  scon   ,  tcon   ,  umax   ,  umin1  ,
     *                 umin2  ,  wtmax  ,  ximax  ,  xk     ,  zcon
      common/ndata/    ipope  ,  jaxi   ,  jmax   ,  maxn   ,  nedit  ,
     *                 nfreq
c=======================================================================
c           Set up the grid
c=======================================
      if(abs(xk-1.).lt.1.e-6) then
        delta=ximax/jend
      else
        delta=ximax*(xk-1.)/(xk**jend-1.)
      endif
      xi(1)=0.
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
        if(jaxi.eq.1) then
          z1(j)=y1(j)
          z3(j)=y3(j)
        endif
 20   continue
      return
      end
      subroutine grow(dldx,jmax)
c=======================================================================
      parameter (jmx=201)
      common/bcom/ deta(jmx),dudxi(jmx),dxi(jmx),ef(jmx),eps(jmx),
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),v(jmx),
     *             wt(jmx),wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx),
     *             z1(jmx),z3(jmx),zdiss(jmx)
c-----------------------------------------------------------------------
      common/extra/    ceps3  ,  cfl    ,  dti    ,  encon  ,  error  ,
     *                 fcon   ,  gcon   ,  jbar   ,  jend   ,  nstep  ,
     *                 pow    ,  scon   ,  tcon   ,  umax   ,  umin1  ,
     *                 umin2  ,  wtmax  ,  ximax  ,  xk     ,  zcon
c=======================================================================
      uhalf=.5*u(1)
      do 10 j=2,jmax
        if(u(j).lt.uhalf.and.u(j-1).ge.uhalf) then
          jedgem=j
          go to 20
        endif
 10   continue
 20   deluu=u(jedgem)-u(jedgem-1)
      if(deluu.eq.0.) then
        deltau=0.
      else
        deltau=(u(jedgem)-uhalf)*deta(jedgem-1)/deluu
      endif
      dldx=eta(jedgem)-deltau
c=======================================
c      Reduce ximax as required
c=======================================
      jmid=jmax/2
      if(u(1).gt.umax.or.u(jmid).lt.umin1.or.u(jmax-10).lt.umin2) then
        ximax=.95*ximax
        call grid
      endif
c=======================================
c     Increase ximax as required
c=======================================
      if(u(jmax-10).gt.umin1) then
        ximax=1.05*ximax
        call grid
      endif
      return
      end
      subroutine namin(icpu)
c=======================================================================
c            COORDINATE READING INPUT DATA FILE 'JET.DAT'
c=======================================================================
      parameter (niu=4,nn=6,nxd=2)
      character a*12,newfil*7,pform*9
      common/data/  xd(nxd)
      common/iodef/ iu(niu)
      common/ndata/ n(nn)
c=======================================================================
c----------------------------------
c         Read I/O units
c----------------------------------
      call namsys(icpu,iin,iv,msdos,newfil,pform)
      open(iin,file='jet.dat',status='old')
      iu(1)=2
      iu(2)=7
      io=iu(1)
      open(io,file='jet.prt',form=pform,status=newfil)
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
c  Open plotting-data disk file
c----------------------------------
      open(iu(2),file='plotj.dat',form=pform,status=newfil)
      return
c=======================================================================
 100  format(1x,a12,i4)
 200  format(1x,'Program  JET  input data:'/,
     *       1x,'_________________________'/)
 300  format(1x,a12,1pe13.6)
 400  format(1x,'_________________________'/)
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
      subroutine start
c=======================================================================
c          S E T   U P   I N I T I A L   P R O F I L E S
c=======================================================================
      parameter (jmx=201)
      common/bcom/ deta(jmx),dudxi(jmx),dxi(jmx),ef(jmx),eps(jmx),
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),v(jmx),
     *             wt(jmx),wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx),
     *             z1(jmx),z3(jmx),zdiss(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx)
      common/ecom/ climit(jmx),fbeta(jmx),uoum(jmx),vextra(jmx)
c-----------------------------------------------------------------------
      common/data/     etin   ,  wtin
      common/extra/    ceps3  ,  cfl    ,  dti    ,  encon  ,  error  ,
     *                 fcon   ,  gcon   ,  jbar   ,  jend   ,  nstep  ,
     *                 pow    ,  scon   ,  tcon   ,  umax   ,  umin1  ,
     *                 umin2  ,  wtmax  ,  ximax  ,  xk     ,  zcon
      common/iodef/    iunit1 ,  iunit2 ,  model  ,  kerng
      common/mixlen/   amix
      common/ndata/    ipope  ,  jaxi   ,  jmax   ,  maxn   ,  nedit  ,
     *                 nfreq
      common/oneeqn/   aa1    ,  aa2    ,  aa3
      common/stress/   alphah ,  betah  ,  c1     ,  gammah ,  onemab
      common/stres2/   cf1    ,  cf2    ,  cf3    ,  cfx    ,  cfy
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c=======================================================================
      betas=.09
      alphas=sqrt(betas)
c=======================================
c     k-omega model coefficients
c=======================================
      if(model.eq.0) then
        alpha=.52
        beta=.0708
        sigma=.5
        sigmad=.125
        sigmas=.6
        clim=7./8.
        pow=sigmad/sigma
c=======================================
c   Stress-omega model coefficients
c=======================================
      elseif(model.eq.1) then
        alpha=.52
        beta=.0708
        sigma=.5
        sigmad=.125
        sigmas=.6
        pow=sigmad/sigma
        c1=1.8
        c2=10./19.
        alphah=(8.+c2)/11.
        betah=(8.*c2-2.)/11.
        gammah=(60.*c2-4.)/55.
        onemal=1.-alphah
        onemab=onemal-betah
        cf1=-betah
        cf2=onemal
        cf3=.66666667*(onemab+.75*gammah)
        cfx=.66666667*(2.*onemal+betah)
        cfy=-.66666667*(onemal+2.*betah)
c=======================================
c    k-epsilon model coefficients
c=======================================
      elseif(model.eq.2) then
        ceps3=.79
        if(kerng.eq.0) then
          alpha=1.44
          beta=1.92
          sigma=10./13.
          sigmas=1.
        else
          alpha=1.42
          beta=1.68
          betas=.085
          sigma=1./1.39
          sigmas=1./1.39
        endif
c=======================================
c  Spalart-Allmaras model coefficients
c=======================================
      elseif(model.eq.4) then
        cb1=.1355
        cb2=.622
        sigma=2./3.
        aa1=cb1
        aa2=(1.-cb2)/sigma
        aa3=-cb2/sigma
c=======================================
c   Mixing length model coefficients
c=======================================
      elseif(model.eq.9) then
        if(jaxi.eq.-1) then
          amix=.037
        elseif(jaxi.eq.0) then
          amix=.024
        else
          amix=.0186
        endif
      endif
c=======================================
c    Compute initial profiles from
c     the mixing-length solution
c=======================================
      if(jaxi.eq.1) then
        et0=2.
        wt0=260.
      else
        et0=.5
        wt0=50.
      endif
      do 30 j=1,jmax
        eta(j)=xi(j)*et0/wt0
        deta(j)=dxi(j)*et0/wt0
        if(jaxi.ne.1) then
          y=eta(j)/.13
          if(y.lt.1.) then
            u(j)=1.-tanh(y)**2
          else
            u(j)=0.
          endif
        else
          u(j)=1./(1.+57.53*eta(j)**2)**2
        endif
        y=4.*eta(j)
        if(y.lt.1.) then
          factor=1.-y**2
          et(j)=et0*factor
          wt(j)=wt0*factor
          if(et(j).le.etin) et(j)=etin
          if(wt(j).le.wtin) wt(j)=wtin
        else
          et(j)=etin
          wt(j)=wtin
        endif
        eps(j)=et0/wt0
 30   continue
c=======================================
c    Impose the integral constraint
c=======================================
      uint=0.
      if(jaxi.ne.1) then
        dum=1.
      else
        dum=0.
      endif
      do 40 j=2,jmax
        dup=u(j)**2
        if(jaxi.eq.1) dup=dup*eta(j)
        uint=uint+.5*(dum+dup)*deta(j-1)
        dum=dup
 40   continue
      uint=gcon*uint
      uint=sqrt(uint)
      do 50 j=1,jmax
        u(j)=u(j)/uint
 50   continue
c=======================================
c         Final computations
c=======================================
      if(model.eq.2) wtin=betas*etin*wtin
      do 60 j=1,jmax
        um(j)=u(j)
        if(model.eq.4) et(j)=eps(j)
        etm(j)=et(j)
        if(model.eq.2) wt(j)=betas*et(j)*wt(j)
        wtm(j)=wt(j)
        zdiss(j)=0.
        fbeta(j)=1.
        dtau(j)=0.
        if(model.eq.1) then
          sigx(j)=.1*et(j)
          sigy(j)=-sigx(j)
          tau(j)=-.3*et(j)
          sigxm(j)=sigx(j)
          sigym(j)=sigy(j)
          taum(j)=tau(j)
        endif
 60   continue
      return
      end
      subroutine tmestp
c=======================================================================
c                 C O M P U T E   T I M E S T E P
c=======================================================================
      parameter (jmx=201)
      common/bcom/ deta(jmx),dudxi(jmx),dxi(jmx),ef(jmx),eps(jmx),
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),v(jmx),
     *             wt(jmx),wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx),
     *             z1(jmx),z3(jmx),zdiss(jmx)
c-----------------------------------------------------------------------
      common/data/     etin   ,  wtin
      common/extra/    ceps3  ,  cfl    ,  dti    ,  encon  ,  error  ,
     *                 fcon   ,  gcon   ,  jbar   ,  jend   ,  nstep  ,
     *                 pow    ,  scon   ,  tcon   ,  umax   ,  umin1  ,
     *                 umin2  ,  wtmax  ,  ximax  ,  xk     ,  zcon
      common/ndata/    ipope  ,  jaxi   ,  jmax   ,  maxn   ,  nedit  ,
     *                 nfreq
c=======================================================================
      dtimax=0.
      do 10 j=2,jend
        dti=abs(v(j))/dxi(j-1)
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
