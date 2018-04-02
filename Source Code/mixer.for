      program mixer
c=======================================================================
c         C O M P R E S S I B L E   M I X I N G   L A Y E R
c
c    Copyright (c) 2006 DCW Industries, Inc.  All rights reserved.
c-----------------------------------------------------------------------
c            MODEL = 0.........................k-omega Model
c            MODEL = 1....................Stress-omega Model
c            MODEL = 2, KERNG = 0...Standard k-epsilon Model
c            MODEL = 2, KERNG > 0........RNG k-epsilon Model
c            MODEL = 4................Spalart-Allmaras Model
c            MODEL = 9...................Mixing Length Model
c            MODEL = 99...................User-Defined Model
c-----------------------------------------------------------------------
c            IMACH = 0.........Sarkar's Compressibility Term
c            IMACH = 1.........Wilcox's Compressibility Term
c            Otherwise..........Zeman's Compressibility Term
c=======================================================================
      parameter (jmx=201)
      common/bcom/ deta(jmx),dudxi(jmx),dxi(jmx),ef(jmx),eps(jmx),
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),wt(jmx),
     *             wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx)
      common/dcom/ h(jmx),hm(jmx),rho(jmx),yox(jmx)
      common/ecom/ climit(jmx),fbeta(jmx),betasm(jmx),vextra(jmx)
c-----------------------------------------------------------------------
      common/data/     etin   ,  gam    ,  prt    ,  rho2   ,  u2ou1  ,
     *                 wtin   ,  xis    ,  xma1   ,  xma2   ,  xmt0
      common/extra/    cflmax ,  dti    ,  error  ,  et2    ,  hcon   ,
     *                 hconbs ,  h2oh1  ,  nstep  ,  jend   ,  jzerom ,
     *                 jzerop ,  pow    ,  pphi   ,  pplo   ,  umax1  ,
     *                 umax2  ,  umin1  ,  umin2  ,  wt2    ,  wtmax  ,
     *                 ximax  ,  xk     ,  xmac   ,  xmsq   ,  xmtsq
      common/iodef/    iunit1 ,  iunit2 ,  model  ,  kerng
      common/mixlen/   amix
      common/ndata/    imach  ,  jmax   ,  jzero  ,  maxn   ,  nedit  ,
     *                 nfreq  ,  nthick
      common/oneeqn/   aa1    ,  aa2    ,  aa3
      common/stress/   alphah ,  betah  ,  c1     ,  gammah ,  onemab
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
      if(xma2.ge.xma1) xma2=0.
      write(iunit1,100)
c=======================================
c    Calculate auxilliary constants
c=======================================
      jzerom=jzero-1
      jzerop=jzero+1
      jend=jmax-1
      et2=etin
      wt2=wtin
      xmtsq=xmt0**2
      xmsq=xma1**2
      hcon=(gam-1.)*xmsq
      h2oh1=1.+.5*hcon
      h2oh1=h2oh1/(1.+.5*(gam-1.)*xma2**2)
      if(xma1.ne.0.) then
        u2ou1=xma2*sqrt(h2oh1)/xma1
      else
        h2oh1=1./rho2
      endif
      xmac=xma1*(1.-u2ou1)/(1.+sqrt(h2oh1))
      error=1.e-5
      if(xma1.le.1.75) then
        ximax=300.
      else
        ximax=400.*(xma1-1.)
      endif
      xk=1.
      umin1=1.e-5
      umax1=1.-umin1
      umin2=1.e-8
      umax2=1.-umin2
      wtmax=100.
      call pitot(xmsq,pphi)
      call pitot(xma2**2,pplo)
c=======================================
c         Generate the grid
c=======================================
      call grid
c=======================================
c        Set initial conditions
c=======================================
      call start
      hconbs=hcon*betas
      nstep=0
c=======================================
c   Identify flow and turbulence model
c=======================================
      write(*,110)
      if(model.eq.0) then
        cflmax=1.5
        write(*,120)
      elseif(model.eq.1) then
        error=1.e-6
        cflmax=2.
        write(*,130)
      elseif(model.eq.2) then
        cflmax=1.
        if(kerng.eq.0) then
          write(*,140)
        else
          write(*,145)
        endif
      elseif(model.eq.4) then
        cflmax=5.
        write(*,150)
      elseif(model.eq.9) then
        cflmax=10.
        write(*,160)
      elseif(model.eq.99) then
        cflmax=1.
        write(*,170)
      endif
      write(*,100)
c=======================================================================
c                 BEGINNING OF MAIN INTEGRATION LOOP
c=======================================================================
 10   nstep=nstep+1
c=======================================
c Calculate quantities needed in EDDY
c=======================================
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
        if(err.ge.erro) then
          erro=err
        endif
 20   continue
c=======================================
c          Converged solution
c=======================================
      if(erro.le.error) then
        call geteta
        call grow(dldx)
        write(iunit1,180) nstep,erro,dldx
        write(*,180) nstep,erro,dldx
        write(*,190)
c=======================================
c       Un-converged solution
c=======================================
      else
        if(mod(nstep,nfreq).eq.0) then
          call geteta
          call grow(dldx)
          write(iunit1,180) nstep,erro,dldx
          write(*,180) nstep,erro,dldx
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
            um(j)=u(j)
            hm(j)=h(j)
            etm(j)=et(j)
            wtm(j)=wt(j)
            if(model.eq.1) then
              sigxm(j)=sigx(j)
              sigym(j)=sigy(j)
              taum(j)=tau(j)
            endif
 30       continue
          go to 10
c=======================================================================
c                   END OF MAIN INTEGRATION LOOP
c=======================================================================
        else
          write(iunit1,200) nstep
          write(*,200) nstep
        endif
      endif
      call edit
      close(iunit1)
c=======================================
c    Write plotting data disk file
c=======================================
      call yoverx
      write(iunit2,210) jmax
      write(iunit2,220) (yox(j),u(j),j=1,jmax)
      endfile(iunit2)
      close(iunit2)
      call geteta
      write(*,230) xmac
      nthick=-1
      call grow(dldx)
      write(*,240) dldx
      nthick=0
      call grow(dldx)
      write(*,250) dldx
      nthick=1
      call grow(dldx)
      write(*,260) dldx
      nthick=2
      call grow(dldx)
      write(*,270) dldx
c=======================================
c    Pause to keep window visible after
c    the run completes...
c=======================================
      pause
c=======================================================================
 100  format(1x,'Convergence history:'/)
 110  format(' MIXING LAYER')
 120  format(' k-omega Model')
 130  format(' Stress-omega Model')
 140  format(' Standard k-epsilon Model')
 145  format(' RNG k-epsilon Model')
 150  format(' Spalart-Allmaras Model')
 160  format(' Mixing Length Model')
 170  format(' User-Defined Model')
 180  format(1x,'At Iteration ',i5,'...Max Error = ',1pe12.5,
     *          '  and dl/dx = ',e10.3)
 190  format(/19x,' ***** SOLUTION CONVERGED *****'/)
 200  format(/8x,'Computation terminated for max iteration',i5/)
 210  format(1x,i4)
 220  format(1p2e14.6)
 230  format(1x,'Convective Mach number     = ',1pe10.3)
 240  format(1x,'Vorticity thickness        = ',1pe10.3)
 250  format(1x,'Kinetic energy thickness   = ',1pe10.3)
 260  format(1x,'Pitot pressure thickness   = ',1pe10.3)
 270  format(1x,'Sullins momentum thickness = ',1pe10.3)
c=======================================================================
      end
      subroutine calcs
c=======================================================================
c C A L C U L A T E   Q U A N T I T I E S   N E E D E D   I N   E D D Y
c=======================================================================
      parameter (jmx=201)
      common/bcom/ deta(jmx),dudxi(jmx),dxi(jmx),ef(jmx),eps(jmx),
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),wt(jmx),
     *             wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx)
      common/dcom/ h(jmx),hm(jmx),rho(jmx),yox(jmx)
      common/ecom/ climit(jmx),fbeta(jmx),betasm(jmx),vextra(jmx)
      common/xcom/ chiwx(jmx)
c-----------------------------------------------------------------------
      common/data/     etin   ,  gam    ,  prt    ,  rho2   ,  u2ou1  ,
     *                 wtin   ,  xis    ,  xma1   ,  xma2   ,  xmt0
      common/extra/    cflmax ,  dti    ,  error  ,  et2    ,  hcon   ,
     *                 hconbs ,  h2oh1  ,  nstep  ,  jend   ,  jzerom ,
     *                 jzerop ,  pow    ,  pphi   ,  pplo   ,  umax1  ,
     *                 umax2  ,  umin1  ,  umin2  ,  wt2    ,  wtmax  ,
     *                 ximax  ,  xk     ,  xmac   ,  xmsq   ,  xmtsq
      common/iodef/    iunit1 ,  iunit2 ,  model  ,  kerng
      common/mixlen/   amix
      common/ndata/    imach  ,  jmax   ,  jzero  ,  maxn   ,  nedit  ,
     *                 nfreq  ,  nthick
      common/stress/   alphah ,  betah  ,  c1     ,  gammah ,  onemab
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c=======================================================================
c   C A L C U L A T E   V A L U E S   F O R   J = 1   A N D   J = JMAX
c-----------------------------------------------------------------------
c  Eddy viscosity and velocity gradient
c=======================================
      rho(1)=1./h2oh1
      rho(jmax)=1.
      if(model.eq.9) then
        eps(jmax)=eps(jend)
        eps(1)=eps(2)
      elseif(model.eq.4) then
        eps(jmax)=0.
        eps(1)=0.
      else
        eps(jmax)=etin/wtin
        eps(1)=et2/(wt2*h2oh1**2)
        if(model.eq.2) then
          eps(1)=betas*et2*eps(1)
          eps(jmax)=betas*etin*eps(jmax)
        endif
        xmtone=2.*xmsq*et2/h2oh1
        xmtinf=2.*xmsq*etin
c=======================================
c   Sarkar's Compressibility Term
c=======================================
        if(imach.eq.0) then
          betasm(1)=1.+xis*xmtone
          betasm(jmax)=1.+xis*xmtinf
c=======================================
c    Wilcox's Compressibility Term
c=======================================
        elseif(imach.eq.1) then
          if(xmtone.le.xmtsq) then
            betasm(1)=1.
          else
            betasm(1)=1.+xis*(xmtone-xmtsq)
          endif
          if(xmtinf.le.xmtsq) then
            betasm(jmax)=1.
          else
            betasm(jmax)=1.+xis*(xmtinf-xmtsq)
          endif
c=======================================
c    Zeman's Compressibility Term
c=======================================
        else
          if(xmtone.le..02/(gam+1.)) then
            betasm(1)=1.
          else
            xmt=sqrt(.5*(gam+1.)*xmtone)
            betasm(1)=1.+xis*(1.-exp(-((xmt-.1)/.6)**2))
          endif
          if(xmtinf.le..02/(gam+1.)) then
            betasm(jmax)=1.
          else
            xmt=sqrt(.5*(gam+1.)*xmtinf)
            betasm(jmax)=1.+xis*(1.-exp(-((xmt-.1)/.6)**2))
          endif
        endif
      endif
c-----------------------------------------------------------------------
c          M A I N   C O M P U T A T I O N A L   L O O P
c-----------------------------------------------------------------------
      do 10 j=2,jend
c=======================================
c  Eddy viscosity and velocity gradient
c=======================================
        climit(j)=1.
        rho(j)=1./h(j)
        dudxi(j)=y2(j)*(u(j+1)-u(j-1))
        if(model.eq.9) then
          eps(j)=amix*rho(j)*sqrt(abs(dudxi(j)))
        elseif(model.eq.4) then
          eps(j)=rho(j)**2*et(j)
        else
          eps(j)=rho(j)**2*et(j)/wt(j)
          if(model.eq.0) then
            if(abs(clim*dudxi(j)/alphas).gt.rho(j)*et(j)) then
              climit(j)=alphas*rho(j)*et(j)/abs(clim*dudxi(j))
            endif
          endif
          if(model.eq.2) eps(j)=betas*et(j)*eps(j)
          xmt2=2.*xmsq*rho(j)*et(j)
c=======================================
c   Sarkar's Compressibility Term
c=======================================
          if(imach.eq.0) then
            betasm(j)=1.+xis*xmt2
c=======================================
c     Wilcox's Compressibility Term
c=======================================
          elseif(imach.eq.1) then
            if(xmt2.le.xmtsq) then
              betasm(j)=1.
            else
              betasm(j)=1.+xis*(xmt2-xmtsq)
            endif
c=======================================
c    Zeman's Compressibility Term
c=======================================
          else
            if(xmt2.le..02/(gam+1.)) then
              betasm(j)=1.
            else
              xmt=sqrt(.5*(gam+1.)*xmt2)
              betasm(j)=1.+xis*(1.-exp(-((xmt-.1)/.6)**2))
            endif
          endif
        endif
 10   continue
      climit(1)=1.
      climit(jmax)=climit(jend)
c=======================================
c            Streamfunction
c=======================================
      ef(jzero)=0.
      ueps1=u(jzero)*eps(jzero)
      do 20 j=jzerop,jmax
        ueps2=u(j)*eps(j)
        ef(j)=ef(j-1)+.5*(ueps1+ueps2)*dxi(j-1)
        ueps1=ueps2
 20   continue
      ueps1=u(jzero)*eps(jzero)
      do 30 j=1,jzerom
        k=jzero-j
        ueps2=u(k)*eps(k)
        ef(k)=ef(k+1)-.5*(ueps1+ueps2)*dxi(k)
        ueps1=ueps2
 30   continue
c=======================================
c    k-omega and Stress-omega models
c=======================================
      if(model.eq.0.or.model.eq.1) then
        do 50 j=2,jend
          vextra(j)=0.
          dkdxi=(1./et(j+1)-1./et(j-1))/(dxi(j-1)+dxi(j))
          dwdxi=dlog(wt(j+1)/wt(j-1))/(dxi(j-1)+dxi(j))
          chik=-dkdxi*dwdxi/rho(j)**2
          if(chik.gt.0.) vextra(j)=1000.
 50     continue
      endif
c=======================================
c     k-epsilon model computations
c=======================================
      if(model.eq.2.and.kerng.ne.0) then
        fbeta(1)=1.
        do 60 j=2,jend
          alamm=abs(dudxi(j))/(betas*et(j))
          factor=alamm**3*(1.-alamm/4.38)/(1.+.012*alamm**3)
          fbeta(j)=1.+betas/beta*factor
 60     continue
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
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),wt(jmx),
     *             wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx)
      common/dcom/ h(jmx),hm(jmx),rho(jmx),yox(jmx)
      common/ecom/ climit(jmx),fbeta(jmx),betasm(jmx),vextra(jmx)
c-----------------------------------------------------------------------
      common/data/     etin   ,  gam    ,  prt    ,  rho2   ,  u2ou1  ,
     *                 wtin   ,  xis    ,  xma1   ,  xma2   ,  xmt0
      common/extra/    cflmax ,  dti    ,  error  ,  et2    ,  hcon   ,
     *                 hconbs ,  h2oh1  ,  nstep  ,  jend   ,  jzerom ,
     *                 jzerop ,  pow    ,  pphi   ,  pplo   ,  umax1  ,
     *                 umax2  ,  umin1  ,  umin2  ,  wt2    ,  wtmax  ,
     *                 ximax  ,  xk     ,  xmac   ,  xmsq   ,  xmtsq
      common/iodef/    iunit1 ,  iunit2 ,  model  ,  kerng
      common/ndata/    imach  ,  jmax   ,  jzero  ,  maxn   ,  nedit  ,
     *                 nfreq  ,  nthick
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
          if(model.eq.0) then
            prodw=-(beta+betas*(1.-betasm(j)))*etm(j)*rho(j)**2
     *            +eps(j)*um(j)
            if(etm(j).ne.0.) prodw=prodw+alpha*climit(j)
     *                               *dudxi(j)**2/etm(j)
          else
            prodw=-(beta+betas*(1.-betasm(j)))*etm(j)*rho(j)**2
     *            +eps(j)*um(j)
            if(etm(j).ne.0.) prodw=prodw+alpha*rho(j)*taum(j)*dudxi(j)
     *                            /etm(j)
          endif
          cross=etm(j)**pow
          prodw=prodw*cross
          if(vextra(j).gt.1.) then
            etmh=(.5*(etm(j-1)+etm(j)))**pow
            etph=(.5*(etm(j+1)+etm(j)))**pow
            a(j)=sigma*etmh*y3(j)-ef(j)*y2(j)*cross
            b(j)=-sigma*(etph*y1(j)+etmh*y3(j))-2.*dti*cross
            c(j)=sigma*etph*y1(j)+ef(j)*y2(j)*cross
          else
            a(j)=(sigma*y3(j)-ef(j)*y2(j))*cross
            b(j)=(-sigma*(y1(j)+y3(j))-2.*dti)*cross
            c(j)=(sigma*y1(j)+ef(j)*y2(j))*cross
          endif
          wt(j)=-a(j)*wtm(j-1)-(b(j)+4.*dti*cross)*wtm(j)-c(j)*wtm(j+1)
        elseif(model.eq.2) then
          prodw=-betas*beta*fbeta(j)*etm(j)*rho(j)**2+eps(j)*um(j)
          if(etm(j).ne.0.) prodw=prodw+alpha*dudxi(j)**2/etm(j)
          a(j)=sigma*y3(j)-ef(j)*y2(j)
          b(j)=-sigma*(y1(j)+y3(j))-2.*dti
          c(j)=sigma*y1(j)+ef(j)*y2(j)
          wt(j)=-a(j)*wtm(j-1)-(b(j)+4.*dti)*wtm(j)-c(j)*wtm(j+1)
        endif
        if(prodw.gt.0.) then
          wt(j)=wt(j)-2.*prodw*wtm(j)
        else
          b(j)=b(j)+2.*prodw
        endif
 10   continue
      call tri(a,b,c,wt,1,jmax,0,wt2,0,wtin)
      do 20 j=2,jend
        if(wt(j).lt.0.) wt(j)=.5*(wtm(j-1)+wtm(j+1))
        if(wt(j).gt.wtmax) wt(j)=wtmax
 20   continue
c=======================================
c      Turbulence kinetic energy
c=======================================
      do 30 j=2,jend
        prode=-betas*betasm(j)*etm(j)*rho(j)**2
        if(etm(j).ne.0.) then
          if(model.eq.0.or.model.eq.2) then
            prode=prode+climit(j)*dudxi(j)**2/etm(j)
          elseif(model.eq.1) then
            prode=prode+rho(j)*taum(j)*dudxi(j)/etm(j)
          endif
        endif
        a(j)=sigmas*y3(j)-ef(j)*y2(j)
        b(j)=-sigmas*(y1(j)+y3(j))-2.*dti
        c(j)=sigmas*y1(j)+ef(j)*y2(j)
        et(j)=-a(j)*etm(j-1)-(b(j)+4.*dti)*etm(j)-c(j)*etm(j+1)
        if(prode.gt.0.) then
          et(j)=et(j)-2.*prode*etm(j)
        else
          b(j)=b(j)+2.*prode
        endif
 30   continue
      call tri(a,b,c,et,1,jmax,0,et2,0,etin)
      do 40 j=2,jend
        if(et(j).lt.0.) et(j)=.5*(etm(j-1)+etm(j+1))
 40   continue
      go to 80
c=======================================
c  One-equation model: eddy viscosity
c=======================================
 50   epsm=.5*aa3*(rho(1)*etm(1)+rho(2)*etm(2))
      rhom=.5*(rho(1)+rho(2))
      do 60 j=2,jend
        epsp=.5*aa3*(etm(j+1)+etm(j))
        rhop=.5*(rho(j+1)+rho(j))
        adiff=rho(j)*(aa2*rhom-epsm/etm(j))
        cdiff=rho(j)*(aa2*rhop-epsp/etm(j))
        prode=aa1*rho(j)*abs(dudxi(j))-u(j)*etm(j)
        a(j)=adiff*y3(j)-ef(j)*y2(j)
        b(j)=-cdiff*y1(j)-adiff*y3(j)-2.*dti
        c(j)=cdiff*y1(j)+ef(j)*y2(j)
        et(j)=-a(j)*etm(j-1)-(b(j)+4.*dti)*etm(j)-c(j)*etm(j+1)
        et(j)=et(j)-2.*prode*etm(j)
        epsm=epsp
        rhom=rhop
 60   continue
      call tri(a,b,c,et,1,jmax,0,0.,0,0.)
      do 70 j=2,jend
        if(et(j).lt.0.) et(j)=.5*(etm(j-1)+etm(j+1))
 70   continue
      if(et(1).lt.0.) et(1)=et(2)
c=======================================
c              Momentum
c=======================================
 80   call calcs
      do 90 j=2,jend
        if(model.ne.1) then
          dalfdx=y2(j)*(climit(j+1)-climit(j-1))
          a(j)=climit(j)*y3(j)-(ef(j)+dalfdx)*y2(j)
          b(j)=-climit(j)*y1(j)-climit(j)*y3(j)-2.*dti
          c(j)=climit(j)*y1(j)+(ef(j)+dalfdx)*y2(j)
        else
          a(j)=y3(j)-ef(j)*y2(j)
          b(j)=-y1(j)-y3(j)-2.*dti
          c(j)=y1(j)+ef(j)*y2(j)
        endif
        u(j)=-a(j)*um(j-1)-(b(j)+4.*dti)*um(j)-c(j)*um(j+1)
        if(model.eq.1) u(j)=u(j)+2.*y2(j)*(dtau(j+1)-dtau(j-1))
 90   continue
      call tri(a,b,c,u,1,jmax,0,u2ou1,0,1.)
c=======================================
c             Mean Energy
c=======================================
      do 100 j=2,jend
        a(j)=y3(j)/prt-ef(j)*y2(j)
        b(j)=-(y1(j)+y3(j))/prt-2.*dti
        c(j)=y1(j)/prt+ef(j)*y2(j)
        h(j)=-a(j)*hm(j-1)-(b(j)+4.*dti)*hm(j)-c(j)*hm(j+1)
        if(model.eq.4.or.model.eq.9) then
          h(j)=h(j)-2.*hcon*dudxi(j)**2
        else
          h(j)=h(j)-2.*hconbs*betasm(j)*(et(j)*rho(j))**2
        endif
 100  continue
      call tri(a,b,c,h,1,jmax,0,h2oh1,0,1.)
c-----------------------------------------------------------------------
c                    Reynolds-stress equations
c-----------------------------------------------------------------------
      if(model.eq.1) then
        call calcs
c=======================================
c             Shear Stress
c=======================================
        do 110 j=2,jend
          a(j)=sigmas*y3(j)-ef(j)*y2(j)
          b(j)=-sigmas*(y1(j)+y3(j))-2.*dti
     *         -betas*betasm(j)*c1*et(j)*rho(j)**2
          c(j)=sigmas*y1(j)+ef(j)*y2(j)
          tau(j)=-a(j)*taum(j-1)-(b(j)+4.*dti)*taum(j)-c(j)*taum(j+1)
     *           -2.*(cf1*sigxm(j)+cf2*sigym(j)
     *           +cf3*et(j))*rho(j)*dudxi(j)
 110    continue
        call tri(a,b,c,tau,1,jmax,0,0.,0,0.)
c=======================================
c          Normal Stresses
c=======================================
        do 130 j=2,jend
          a(j)=sigmas*y3(j)-ef(j)*y2(j)
          b(j)=-sigmas*(y1(j)+y3(j))-2.*dti
     *         -betas*betasm(j)*c1*et(j)*rho(j)**2
          c(j)=sigmas*y1(j)+ef(j)*y2(j)
          sigx(j)=-a(j)*sigxm(j-1)-(b(j)+4.*dti)*sigxm(j)
     *            -c(j)*sigxm(j+1)-2.*cfx*rho(j)*tau(j)*dudxi(j)
 130    continue
        call tri(a,b,c,sigx,1,jmax,0,0.,0,0.)
        do 140 j=2,jend
          a(j)=sigmas*y3(j)-ef(j)*y2(j)
          b(j)=-sigmas*(y1(j)+y3(j))-2.*dti
     *         -betas*betasm(j)*c1*et(j)*rho(j)**2
          c(j)=sigmas*y1(j)+ef(j)*y2(j)
          sigy(j)=-a(j)*sigym(j-1)-(b(j)+4.*dti)*sigym(j)
     *            -c(j)*sigym(j+1)-2.*cfy*rho(j)*tau(j)*dudxi(j)
 140    continue
        call tri(a,b,c,sigy,1,jmax,0,0.,0,0.)
c=======================================
c        Realizability checks
c=======================================
        relax=.05
        relaxm=1.-relax
        do 150 j=2,jend
          dtau0=rho(j)*tau(j)-dudxi(j)
          dtau(j)=relax*dtau0+relaxm*dtau(j)
          term=.66666667*et(j)
          upup=sigx(j)+term
          if(upup.lt.0.) sigx(j)=-term
          if(upup.gt.et(j)) sigx(j)=.5*term
          vpvp=sigy(j)+term
          if(vpvp.lt.0.) sigy(j)=-term
          if(vpvp.gt.et(j)) sigy(j)=.5*term
 150    continue
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
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),wt(jmx),
     *             wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx)
      common/bdkw/ cnvk(jmx),cnvw(jmx),dflk(jmx),dflw(jmx),
     *             dftk(jmx),dftw(jmx),disk(jmx),disw(jmx),
     *             esk(jmx),esw(jmx),prdk(jmx),prdw(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx)
      common/dcom/ h(jmx),hm(jmx),rho(jmx),yox(jmx)
      common/ecom/ climit(jmx),fbeta(jmx),betasm(jmx),vextra(jmx)
      common/xcom/ chiwx(jmx)
c-----------------------------------------------------------------------
      common/data/     etin   ,  gam    ,  prt    ,  rho2   ,  u2ou1  ,
     *                 wtin   ,  xis    ,  xma1   ,  xma2   ,  xmt0
      common/extra/    cflmax ,  dti    ,  error  ,  et2    ,  hcon   ,
     *                 hconbs ,  h2oh1  ,  nstep  ,  jend   ,  jzerom ,
     *                 jzerop ,  pow    ,  pphi   ,  pplo   ,  umax1  ,
     *                 umax2  ,  umin1  ,  umin2  ,  wt2    ,  wtmax  ,
     *                 ximax  ,  xk     ,  xmac   ,  xmsq   ,  xmtsq
      common/iodef/    iunit1 ,  iunit2 ,  model  ,  kerng
      common/ndata/    imach  ,  jmax   ,  jzero  ,  maxn   ,  nedit  ,
     *                 nfreq  ,  nthick
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c=======================================================================
      call geteta
      call grow(dldx)
      call yoverx
      write(iunit1,100) nstep,dldx,xmac
      if(model.ne.1) then
        do 10 j=1,jmax
          tau(j)=h(j)*dudxi(j)
 10     continue
        write(iunit1,110)
        write(iunit1,120) (j,yox(j),u(j),ef(j),h(j),et(j),wt(j),
     *                                  eps(j),tau(j),j=1,jmax)
      else
        do 20 j=1,jmax
          term=.666666667*et(j)
          upup(j)=term+sigx(j)
          vpvp(j)=term+sigy(j)
 20     continue
        write(iunit1,130)
        write(iunit1,140) (j,yox(j),u(j),ef(j),h(j),et(j),wt(j),eps(j),
     *                                tau(j),upup(j),vpvp(j),j=1,jmax)
      endif
      return
c=======================================================================
 100  format(/1x,'nstep =',i6,5x,'dl/dx =',1pe11.4,5x,'Mc =',e10.3)
 110  format(/5x,'j',4x,'y/x',7x,'u',9x,'f',9x,'h',9x,'k',8x,'omega',
     *        6x,'nut',7x,'tau')
 120  format(i6,1p8e10.3)
 130  format(/5x,'j',4x,'y/x',7x,'u',9x,'f',9x,'h',9x,'k',8x,'omega',
     *        6x,'nut',7x,'tau',7x,'uu',8x,'vv')
 140  format(i6,1p10e10.3)
c=======================================================================
      end
      subroutine geteta
c=======================================================================
c                      C O M P U T E   E T A
c=======================================================================
      parameter (jmx=201)
      common/bcom/ deta(jmx),dudxi(jmx),dxi(jmx),ef(jmx),eps(jmx),
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),wt(jmx),
     *             wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx)
      common/dcom/ h(jmx),hm(jmx),rho(jmx),yox(jmx)
c-----------------------------------------------------------------------
      common/data/     etin   ,  gam    ,  prt    ,  rho2   ,  u2ou1  ,
     *                 wtin   ,  xis    ,  xma1   ,  xma2   ,  xmt0
      common/extra/    cflmax ,  dti    ,  error  ,  et2    ,  hcon   ,
     *                 hconbs ,  h2oh1  ,  nstep  ,  jend   ,  jzerom ,
     *                 jzerop ,  pow    ,  pphi   ,  pplo   ,  umax1  ,
     *                 umax2  ,  umin1  ,  umin2  ,  wt2    ,  wtmax  ,
     *                 ximax  ,  xk     ,  xmac   ,  xmsq   ,  xmtsq
      common/ndata/    imach  ,  jmax   ,  jzero  ,  maxn   ,  nedit  ,
     *                 nfreq  ,  nthick
c=======================================================================
c    Above the dividing streamline
c=======================================
      eta(jzero)=0.
      eps1=eps(jzero)
      do 10 j=jzerop,jend
        eps2=eps(j)
        eta(j)=eta(j-1)+.5*(eps1+eps2)*dxi(j-1)
        deta(j-1)=eta(j)-eta(j-1)
        eps1=eps2
 10   continue
      deta(jend)=deta(jend-1)
      eta(jmax)=eta(jend)+deta(jend)
c=======================================
c    Below the dividing streamline
c=======================================
      eps1=eps(jzero)
      do 20 j=jzerom,2,-1
        eps2=eps(j)
        eta(j)=eta(j+1)-.5*(eps1+eps2)*dxi(j)
        deta(j)=eta(j+1)-eta(j)
        eps1=eps2
 20   continue
      deta(1)=deta(2)
      eta(1)=eta(2)-deta(1)
      return
      end
      subroutine grid
c=======================================================================
c  G E N E R A T E   T H E   F I N I T E   D I F F E R E N C E   G R I D
c=======================================================================
      parameter (jmx=201)
      common/bcom/ deta(jmx),dudxi(jmx),dxi(jmx),ef(jmx),eps(jmx),
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),wt(jmx),
     *             wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx)
c-----------------------------------------------------------------------
      common/data/     etin   ,  gam    ,  prt    ,  rho2   ,  u2ou1  ,
     *                 wtin   ,  xis    ,  xma1   ,  xma2   ,  xmt0
      common/extra/    cflmax ,  dti    ,  error  ,  et2    ,  hcon   ,
     *                 hconbs ,  h2oh1  ,  nstep  ,  jend   ,  jzerom ,
     *                 jzerop ,  pow    ,  pphi   ,  pplo   ,  umax1  ,
     *                 umax2  ,  umin1  ,  umin2  ,  wt2    ,  wtmax  ,
     *                 ximax  ,  xk     ,  xmac   ,  xmsq   ,  xmtsq
      common/ndata/    imach  ,  jmax   ,  jzero  ,  maxn   ,  nedit  ,
     *                 nfreq  ,  nthick
c=======================================================================
      nn=jmax-jzero
      if(abs(xk-1.).lt.1.e-6) then
        delta=ximax/nn
      else
        delta=ximax*(xk-1.)/(xk**nn-1.)
      endif
c=======================================
c    Above the dividing streamline
c=======================================
      xi(jzero)=0.
      dxi(jzero)=delta
      do 10 j=jzerop,jend
        dxi(j)=xk*dxi(j-1)
        xi(j)=xi(j-1)+dxi(j-1)
 10   continue
      dxi(jmax)=0.
      xi(jmax)=ximax
c=======================================
c    Below the dividing streamline
c=======================================
      do 20 j=jzerom,1,-1
        dxi(j)=xk*dxi(j+1)
        xi(j)=xi(j+1)-dxi(j)
 20   continue
c=======================================
c   Compute differencing coefficients
c=======================================
      do 30 j=2,jend
        denom=dxi(j)+dxi(j-1)
        y1(j)=2./(dxi(j)*denom)
        y2(j)=1./denom
        y3(j)=2./(dxi(j-1)*denom)
 30   continue
      return
      end
      subroutine grow(dldx)
c=======================================================================
      parameter (jmx=201)
      common/bcom/ deta(jmx),dudxi(jmx),dxi(jmx),ef(jmx),eps(jmx),
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),wt(jmx),
     *             wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx)
      common/dcom/ h(jmx),hm(jmx),rho(jmx),yox(jmx)
c-----------------------------------------------------------------------
      common/data/     etin   ,  gam    ,  prt    ,  rho2   ,  u2ou1  ,
     *                 wtin   ,  xis    ,  xma1   ,  xma2   ,  xmt0
      common/extra/    cflmax ,  dti    ,  error  ,  et2    ,  hcon   ,
     *                 hconbs ,  h2oh1  ,  nstep  ,  jend   ,  jzerom ,
     *                 jzerop ,  pow    ,  pphi   ,  pplo   ,  umax1  ,
     *                 umax2  ,  umin1  ,  umin2  ,  wt2    ,  wtmax  ,
     *                 ximax  ,  xk     ,  xmac   ,  xmsq   ,  xmtsq
      common/ndata/    imach  ,  jmax   ,  jzero  ,  maxn   ,  nedit  ,
     *                 nfreq  ,  nthick
c=======================================================================
      dldx=0.
c=======================================
c   Bogdanoff's vorticity thickness
c=======================================
      if(nthick.eq.-1) then
        vortmx=0.
        do 10 j=1,jmax
          if(eps(j).ne.0.) then
            vort=abs(rho(j)*dudxi(j)/eps(j))
          else
            vort=0.
          endif
          if(vort.gt.vortmx) vortmx=vort
 10     continue
        if(vortmx.ne.0.) then
          dldx=1./vortmx
        else
          dldx=0.
        endif
c=======================================
c      Birch's energy thickness
c=======================================
      elseif(nthick.eq.0) then
        ulow=u2ou1+.316*(1.-u2ou1)
        uhi=u2ou1+.949*(1.-u2ou1)
        do 20 j=1,jmax
          if(u(j).ge.ulow.and.u(j-1).lt.ulow) jl=j
          if(u(j).ge.uhi.and.u(j-1).lt.uhi) ju=j
          if(u(j).ge.ulow.and.u(j).le.uhi) then
            dldx=dldx+.5*(h(j)+h(j+1))*deta(j)
          endif
 20     continue
        delul=u(jl)-u(jl-1)
        if(delul.eq.0.) then
          deltal=0.
        else
          hl=.5*(h(jl)+h(jl-1))
          deltal=(u(jl)-ulow)*hl*deta(jl-1)/delul
        endif
        deluu=u(ju)-u(ju-1)
        if(deluu.eq.0.) then
          deltau=0.
        else
          hu=.5*(h(ju)+h(ju-1))
          deltau=(u(ju)-uhi)*hu*deta(ju-1)/deluu
        endif
        dldx=dldx+deltal-deltau
c=======================================
c      Roshko's pitot thickness
c=======================================
      elseif(nthick.eq.1) then
        fractu=0.
        fractl=0.
        if(xma1.ge..3) then
          pprato=1.
          do 30 j=1,jmax
            em2=xmsq*rho(j)*u(j)**2
            call pitot(em2,pp)
            pprat=(pp-pplo)/(pphi-pplo)
            if(pprato.ge..95.and.pprat.lt..95) then
              jl=j
              if(pprato.ne.pprat) fractl=(pprat-.95)/(pprat-pprato)
            endif
            if(pprato.ge..05.and.pprat.lt..05) then
              ju=j
              if(pprato.ne.pprat) fractu=(pprat-.05)/(pprat-pprato)
            endif
            if(pprat.ge..05.and.pprat.le..95) then
              dldx=dldx+.5*(h(j)+h(j+1))*deta(j)
            endif
            pprato=pprat
 30       continue
        else
          ulow=u2ou1+.2236*(1.-u2ou1)
          uhi=u2ou1+.9747*(1.-u2ou1)
          do 40 j=1,jmax
            if(u(j).ge.ulow.and.u(j-1).lt.ulow) then
              jl=j
              if(u(jl).ne.u(jl-1)) fractl=(u(jl)-ulow)/(u(jl)-u(jl-1))
            endif
            if(u(j).ge.uhi.and.u(j-1).lt.uhi) then
              ju=j
              if(u(ju).ne.u(ju-1)) fractu=(u(ju)-uhi)/(u(ju)-u(ju-1))
            endif
            if(u(j).ge.ulow.and.u(j).le.uhi) then
              dldx=dldx+.5*(h(j)+h(j+1))*deta(j)
            endif
 40       continue
        endif
        hl=.5*(h(jl)+h(jl-1))
        hu=.5*(h(ju)+h(ju-1))
        deltal=fractl*hl*deta(jl-1)
        deltau=fractu*hu*deta(ju-1)
        dldx=dldx+deltal-deltau
c=======================================
c      Sullin's momentum thickness
c=======================================
      elseif(nthick.eq.2) then
        ulow=u2ou1+.1*(1.-u2ou1)
        uhi=u2ou1+.9*(1.-u2ou1)
        do 50 j=1,jmax
          if(u(j).ge.ulow.and.u(j-1).lt.ulow) jl=j
          if(u(j).ge.uhi.and.u(j-1).lt.uhi) ju=j
          if(u(j).ge.ulow.and.u(j).le.uhi) then
            dldx=dldx+.5*(h(j)+h(j+1))*deta(j)
          endif
 50     continue
        delul=u(jl)-u(jl-1)
        if(delul.eq.0.) then
          deltal=0.
        else
          hl=.5*(h(jl)+h(jl-1))
          deltal=(u(jl)-ulow)*hl*deta(jl-1)/delul
        endif
        deluu=u(ju)-u(ju-1)
        if(deluu.eq.0.) then
          deltau=0.
        else
          hu=.5*(h(ju)+h(ju-1))
          deltau=(u(ju)-uhi)*hu*deta(ju-1)/deluu
        endif
        dldx=dldx+deltal-deltau
      endif
c=======================================
c      Reduce ximax as required
c=======================================
      if(u(3).lt.umin2.and.u(jmax-2).gt.umax2) then
        ximax=.95*ximax
        call grid
      endif
c=======================================
c     Increase ximax as required
c=======================================
      if(u(3).gt.umin1.and.u(jmax-2).lt.umax1) then
        ximax=1.05*ximax
        call grid
      endif
      return
      end
      subroutine namin(icpu)
c=======================================================================
c            COORDINATE READING INPUT DATA FILE 'MIXER.DAT'
c=======================================================================
      parameter (niu=4,nn=7,nxd=10)
      character a*12,newfil*7,pform*9
      common/data/  xd(nxd)
      common/iodef/ iu(niu)
      common/ndata/ n(nn)
c=======================================================================
c----------------------------------
c         Read I/O units
c----------------------------------
      call namsys(icpu,iin,iv,msdos,newfil,pform)
      open(iin,file='mixer.dat',status='old')
      iu(1)=2
      iu(2)=7
      io=iu(1)
      open(io,file='mixer.prt',form=pform,status=newfil)
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
c   Open plotting-data disk file
c----------------------------------
      open(iu(2),file='plotm.dat',form=pform,status=newfil)
      return
c=======================================================================
 100  format(1x,a12,i4)
 200  format(1x,'Program MIXER input data:'/,
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
      subroutine pitot(em2,pp)
c=======================================================================
c            C O M P U T E   P I T O T   P R E S S U R E
c=======================================================================
      common/data/     etin   ,  gam    ,  prt    ,  rho2   ,  u2ou1  ,
     *                 wtin   ,  xis    ,  xma1   ,  xma2   ,  xmt0
c=======================================================================
      gm1=gam-1.
      gogm1=gam/gm1
      if(em2.le.1.) then
        pp=(1.+.5*gm1*em2)**gogm1
      else
        gp1=gam+1.
        pp=(.5*gp1*em2)**gogm1*(2.*gam*em2/gp1-gm1/gp1)**(-1./gm1)
      endif
      return
      end
      subroutine start
c=======================================================================
c          S E T   U P   I N I T I A L   P R O F I L E S
c=======================================================================
      parameter (jmx=201)
      common/bcom/ deta(jmx),dudxi(jmx),dxi(jmx),ef(jmx),eps(jmx),
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),wt(jmx),
     *             wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx)
      common/dcom/ h(jmx),hm(jmx),rho(jmx),yox(jmx)
      common/ecom/ climit(jmx),fbeta(jmx),betasm(jmx),vextra(jmx)
c-----------------------------------------------------------------------
      common/data/     etin   ,  gam    ,  prt    ,  rho2   ,  u2ou1  ,
     *                 wtin   ,  xis    ,  xma1   ,  xma2   ,  xmt0
      common/extra/    cflmax ,  dti    ,  error  ,  et2    ,  hcon   ,
     *                 hconbs ,  h2oh1  ,  nstep  ,  jend   ,  jzerom ,
     *                 jzerop ,  pow    ,  pphi   ,  pplo   ,  umax1  ,
     *                 umax2  ,  umin1  ,  umin2  ,  wt2    ,  wtmax  ,
     *                 ximax  ,  xk     ,  xmac   ,  xmsq   ,  xmtsq
      common/iodef/    iunit1 ,  iunit2 ,  model  ,  kerng
      common/mixlen/   amix
      common/ndata/    imach  ,  jmax   ,  jzero  ,  maxn   ,  nedit  ,
     *                 nfreq  ,  nthick
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
        amix=.0175
      endif
c=======================================
c  Mean velocity and enthalpy profiles
c=======================================
      ucon1=.5*(1.+u2ou1)
      ucon2=.5*(1.-u2ou1)
      hcon1=.5*(1.+h2oh1)
      hcon2=.5*(1.-h2oh1)
      do 10 j=2,jend
        erfeta=derf(xi(j)/54.)
        u(j)=ucon1+ucon2*erfeta
        h(j)=hcon1+hcon2*erfeta
 10   continue
      u(1)=u2ou1
      h(1)=h2oh1
      u(jmax)=1.
      h(jmax)=1.
c=======================================
c   Now compute the turbulence kinetic
c      energy and dissipation rate
c=======================================
      do 20 j=2,jend
        climit(j)=1.
        y=abs(xi(j)/108.)
        if(y.le.1.) then
          et(j)=.003*(cos(1.5707963*y))**2
          if(j.gt.jzero.and.et(j).lt.etin) et(j)=etin
          if(j.lt.jzero.and.et(j).lt.et2) et(j)=et2
          wt(j)=730.*et(j)
          if(j.gt.jzero.and.wt(j).lt.wtin) wt(j)=wtin
          if(j.lt.jzero.and.wt(j).lt.wt2) wt(j)=wt2
        elseif(j.gt.jzero) then
          et(j)=etin
          wt(j)=wtin
        else
          et(j)=et2
          wt(j)=wt2
        endif
 20   continue
      et(1)=et2
      wt(1)=wt2
      et(jmax)=etin
      wt(jmax)=wtin
      climit(1)=1.
      climit(jmax)=1.
c=======================================
c         Final computations
c=======================================
      if(model.eq.2) then
        wtin=betas*etin*wtin
        wt2=betas*et2*wt2
      endif
      do 30 j=1,jmax
        um(j)=u(j)
        hm(j)=h(j)
        rho(j)=1./h(j)
        etm(j)=et(j)
        if(model.eq.4) then
          et(j)=et(j)/wt(j)
        else
          eps(j)=rho(j)**2*et(j)/wt(j)
          if(model.eq.2) wt(j)=betas*et(j)*wt(j)
        endif
        wtm(j)=wt(j)
        fbeta(j)=1.
        dtau(j)=0.
        if(model.eq.1) then
          sigx(j)=.1*et(j)
          sigy(j)=-sigx(j)
          tau(j)=.3*et(j)
          sigxm(j)=sigx(j)
          sigym(j)=sigy(j)
          taum(j)=tau(j)
        endif
 30   continue
      return
      end
      subroutine tmestp
c=======================================================================
c             S E T   U P   T I M E S T E P   A R R A Y
c=======================================================================
      parameter (jmx=201)
      common/bcom/ deta(jmx),dudxi(jmx),dxi(jmx),ef(jmx),eps(jmx),
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),wt(jmx),
     *             wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx)
c-----------------------------------------------------------------------
      common/data/     etin   ,  gam    ,  prt    ,  rho2   ,  u2ou1  ,
     *                 wtin   ,  xis    ,  xma1   ,  xma2   ,  xmt0
      common/extra/    cflmax ,  dti    ,  error  ,  et2    ,  hcon   ,
     *                 hconbs ,  h2oh1  ,  nstep  ,  jend   ,  jzerom ,
     *                 jzerop ,  pow    ,  pphi   ,  pplo   ,  umax1  ,
     *                 umax2  ,  umin1  ,  umin2  ,  wt2    ,  wtmax  ,
     *                 ximax  ,  xk     ,  xmac   ,  xmsq   ,  xmtsq
      common/ndata/    imach  ,  jmax   ,  jzero  ,  maxn   ,  nedit  ,
     *                 nfreq  ,  nthick
c=======================================================================
c    Ramp the CFL number up during 
c       the first 200 timesteps
c=======================================
      if(nstep.lt.100) then
        cfl=.25*cflmax
      elseif(nstep.lt.200) then
        cfl=cflmax*(.25+.75*(nstep-100)/100.)
      else
        cfl=cflmax
      endif
c=======================================
c     Now determine the timestep
c=======================================
      dtimax=0.
      do 10 j=2,jend
        dti=abs(ef(j))/dxi(j-1)
        dtimax=dmax1(dtimax,dti)
 10   continue
      dti=dtimax/cfl
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
      subroutine yoverx
c=======================================================================
      parameter (jmx=201)
      common/bcom/ deta(jmx),dudxi(jmx),dxi(jmx),ef(jmx),eps(jmx),
     *             et(jmx),eta(jmx),etm(jmx),u(jmx),um(jmx),wt(jmx),
     *             wtm(jmx),xi(jmx),y1(jmx),y2(jmx),y3(jmx)
      common/dcom/ h(jmx),hm(jmx),rho(jmx),yox(jmx)
c-----------------------------------------------------------------------
      common/extra/    cflmax ,  dti    ,  error  ,  et2    ,  hcon   ,
     *                 hconbs ,  h2oh1  ,  nstep  ,  jend   ,  jzerom ,
     *                 jzerop ,  pow    ,  pphi   ,  pplo   ,  umax1  ,
     *                 umax2  ,  umin1  ,  umin2  ,  wt2    ,  wtmax  ,
     *                 ximax  ,  xk     ,  xmac   ,  xmsq   ,  xmtsq
      common/ndata/    imach  ,  jmax   ,  jzero  ,  maxn   ,  nedit  ,
     *                 nfreq  ,  nthick
c=======================================================================
      yox(jzero)=0.
      do 10 j=jzerop,jmax
        yox(j)=yox(j-1)+.5*(h(j)+h(j-1))*deta(j-1)
 10   continue
      do 20 j=jzerom,1,-1
        yox(j)=yox(j+1)-.5*(h(j)+h(j+1))*deta(j)
 20   continue
      return
      end
