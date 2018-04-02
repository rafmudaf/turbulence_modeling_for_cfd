      program pipe
c=======================================================================
c I N C O M P R E S S I B L E   C H A N N E L   O R   P I P E   F L O W
c
c    Copyright (c) 2006 DCW Industries, Inc.  All rights reserved.
c-----------------------------------------------------------------------
c                 MODEL = 0.............k-omega Model
c                 MODEL = 1........Stress-omega Model
c                 MODEL = 3.......Baldwin-Barth Model
c                 MODEL = 4....Spalart-Allmaras Model
c                 MODEL = 5........Cebeci-Smith Model
c                 MODEL = 6.......Baldwin-Lomax Model
c                 MODEL = 7........Johnson-King Model
c                 MODEL = 99.......User-Defined Model
c-----------------------------------------------------------------------
c                       JAXI = 0...Channel Flow
c                       JAXI = 1......Pipe Flow
c=======================================================================
      parameter (jmx=201)
      common/acom/ elmix(jmx),f2(jmx),fw(jmx),stilde(jmx)
      common/bcom/ dlti(jmx),dudy(jmx),dyplus(jmx),eps(jmx),et(jmx),
     *             etm(jmx),ret(jmx),uplus(jmx),uplusm(jmx),wt(jmx),
     *             wtm(jmx),y1(jmx),y2(jmx),y3(jmx),yoh(jmx),yplus(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx),vextra(jmx)
      common/pcom/ diss(jmx),prod(jmx)
      common/ucom/ uoum(jmx),upup(jmx),vpvp(jmx)
      common/vcom/ bbeta(jmx),climit(jmx),fc1(jmx),fgamma(jmx),
     *             fgammas(jmx),gamma(jmx),gammas(jmx)
c-----------------------------------------------------------------------
      common/data/     omegw  ,  retau  ,  yone
      common/extra/    akappa ,  cf     ,  cfl    ,  erro   ,  error  ,
     *                 jskipm ,  jskipp ,  nstep  ,  reh    ,  wn
      common/fctcom/   aa     ,  aam    ,  jend   ,  jendm
      common/iodef/    iunit1 ,  iunit2 ,  model
      common/mixlen/   aplus  ,  c1ml   ,  ccp    ,  ckleb  ,  clausr ,  
     *                 cwk
      common/ndata/    iruff  ,  jaxi   ,  jmax   ,  jskip  ,  maxn   ,
     *                 nedit  ,  nfreq  ,  nvisc
      common/stress/   alphah ,  betah  ,  c1     ,  gammah ,  onemab
      common/stres2/   cf1    ,  cf2    ,  cf3    ,  cfx    ,  cfy
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
      common/viscus/   alf0   ,  alfs0  ,  bbeta0 ,  gamh0  ,  rb     ,
     *                 rk     ,  rw
c-----------------------------------------------------------------------
      include 'cpuid.'
c=======================================================================
c           Read input data
c=======================================
      call namin(icpu)
c=======================================
c    Calculate auxilliary constants
c=======================================
      jend=jmax-1
      jendm=jend-1
      jskipm=jskip-1
      jskipp=jskip+1
      error=1.e-5
      if(model.ne.4) then
        cfl=4.
      else
        if(retau.gt.150.) then
          cfl=.5
        else
          cfl=.08
        error=5.e-4
        endif
      endif
c=======================================
c  Generate the grid and timestep array
c=======================================
      aa=yone/retau
      aam=1.-aa
      call grid
      do 10 j=2,jend
        dlti(j)=1./(cfl*dyplus(j-1))
 10   continue
      dlti(1)=dlti(2)
c=======================================
c       Set initial conditions
c=======================================
      call start
      wo=1.-wn
      nstep=0
c=======================================
c   Identify flow and turbulence model
c=======================================
      if(jaxi.eq.0) then
        write(iunit1,100)
        write(*,100)
      else
        write(iunit1,110)
        write(*,110)
      endif
      if(model.eq.0) then
        write(iunit1,120)
        write(*,120)
      elseif(model.eq.1) then
        write(iunit1,130)
        write(*,130)
      elseif(model.eq.3) then
        write(iunit1,140)
        write(*,140)
      elseif(model.eq.4) then
        write(iunit1,150)
        write(*,150)
      elseif(model.eq.5) then
        write(iunit1,160)
        write(*,160)
      elseif(model.eq.6) then
        write(iunit1,170)
        write(*,170)
      elseif(model.eq.7) then
        write(iunit1,180)
        write(*,180)
      elseif(model.eq.99) then
        write(iunit1,190)
        write(*,190)
      endif
      write(iunit1,200)
      write(*,200)
c=======================================================================
c                  BEGINNING OF MAIN INTEGRATION LOOP
c=======================================================================
 20   nstep=nstep+1
c=======================================
c  Calculate quantities needed in EDDY
c=======================================
      if(model.eq.0.or.model.eq.1) then
        call calcs
      elseif(model.eq.3) then
        call bbarth
      elseif(model.eq.4) then
        call spalart
      elseif(model.eq.5) then
        call csmith
      elseif(model.eq.6) then
        call blomax
      elseif(model.eq.7) then
        call jking
      endif
c=======================================
c Calculate velocity, turbulence energy
c        and dissipation rate
c=======================================
      call eddy
c=======================================
c           Check convergence
c=======================================
      erro=0.
      do 30 j=1,jmax
        err=abs(uplus(j)-uplusm(j))
        if(err.ge.erro) then
          erro=err
        endif
 30   continue
c=======================================
c         Converged solution
c=======================================
      if(erro.le.error) then
        write(iunit1,210) nstep,erro,reh,cf
        write(*,210) nstep,erro,reh,cf
        write(*,220)
c=======================================
c       Un-converged solution
c=======================================
      else
        if(mod(nstep,nfreq).eq.0) then
          write(iunit1,210) nstep,erro,reh,cf
          write(*,210) nstep,erro,reh,cf
        endif
        if(nstep.ne.maxn) then
c=======================================
c            Print profiles
c=======================================
          if(mod(nstep,nedit).eq.0) call edit
c=======================================
c      Prepare for next timestep
c    for all but algebraic models
c=======================================
          if(model.ge.0.and.model.le.4) then
            do 40 j=1,jmax
              etm(j)=et(j)
              uplusm(j)=uplus(j)
              wtm(j)=wt(j)
 40         continue
c=======================================
c  Use relaxation for algebraic models
c=======================================
          else
            do 50 j=1,jmax
              uplus(j)=wo*uplusm(j)+wn*uplus(j)
              uplusm(j)=uplus(j)
 50         continue
          endif
          go to 20
c=======================================================================
c  	       	     END OF MAIN INTEGRATION LOOP
c=======================================================================
        else
          write(iunit1,230) nstep
          write(*,230) nstep
        endif
      endif
      call edit
      close(iunit1)
c=======================================
c    Write plotting data disk file
c=======================================
      write(iunit2,240) jmax,jaxi,retau
      if(model.ge.3.and.model.le.7) then
        ncols=5
        write(iunit2,250) ncols
        write(iunit2,260) (yoh(j),uoum(j),yplus(j),uplus(j),tau(j),
     *                     j=1,jmax)
      else
        ncols=8
        write(iunit2,250) ncols
        do 60 j=2,jmax
          diss(j)=-betas*bbeta(j)*et(j)*wt(j)
          prod(j)=tau(j)*dudy(j)
 60     continue
        if(nvisc.ne.0) diss(1)=2.*diss(2)-diss(3)
        prod(1)=0.
        write(iunit2,270) (yoh(j),uoum(j),yplus(j),uplus(j),tau(j),
     *                     et(j),diss(j),prod(j),j=1,jmax)
      endif
      endfile(iunit2)
      close(iunit2)
c=======================================
c    Pause to keep window visible after
c    the run completes...
c=======================================
      pause
c=======================================================================
 100  format(' CHANNEL FLOW')
 110  format(' PIPE FLOW')
 120  format(' k-omega Model')
 130  format(' Stress-omega Model')
 140  format(' Baldwin-Barth Model')
 150  format(' Spalart-Allmaras Model')
 160  format(' Cebeci-Smith Model')
 170  format(' Baldwin-Lomax Model')
 180  format(' Johnson-King Model')
 190  format(' User-Defined Model')
 200  format(1x,'Convergence history:'/)
 210  format(1x,'At Iteration',i5,'...Max Error =',1pe10.3,
     *          '  Re =',e10.3,'  Cf =',e10.3)
 220  format(/21x,'***** SOLUTION CONVERGED *****'/)
 230  format(/13x,'Computation terminated for max iteration',i5/)
 240  format(1x,2i4,1pe14.6)
 250  format(1x,i2)
 260  format(1p5e14.6)
 270  format(1p8e14.6)
c=======================================================================
      end
      subroutine bbarth
c=======================================================================
c                B A L D W I N - B A R T H   M O D E L
c=======================================================================
      parameter (jmx=201)
      dimension    a(jmx),b(jmx),c(jmx)
      common/acom/ elmix(jmx),f2(jmx),fw(jmx),stilde(jmx)
      common/bcom/ dlti(jmx),dudy(jmx),dyplus(jmx),eps(jmx),et(jmx),
     *             etm(jmx),ret(jmx),uplus(jmx),uplusm(jmx),wt(jmx),
     *             wtm(jmx),y1(jmx),y2(jmx),y3(jmx),yoh(jmx),yplus(jmx)
c-----------------------------------------------------------------------
      common/data/     omegw  ,  retau  ,  yone
      common/extra/    akappa ,  cf     ,  cfl    ,  erro   ,  error  ,
     *                 jskipm ,  jskipp ,  nstep  ,  reh    ,  wn
      common/fctcom/   aa     ,  aam    ,  jend   ,  jendm
      common/mixlen/   aplus  ,  c1ml   ,  ccp    ,  ckleb  ,  clausr ,  
     *                 cwk
      common/ndata/    iruff  ,  jaxi   ,  jmax   ,  jskip  ,  maxn   ,
     *                 nedit  ,  nfreq  ,  nvisc
      common/oneeqn/   aplus0 ,  aplus2 ,  cb1    ,  cb2    ,  ceps1  ,  
     *                 ceps2  ,  cmu    ,  cv1    ,  cw1    ,  cw2    ,
     *                 cw3    ,  sigmae
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c=======================================================================
c  NOTE: The ET and ETM arrays are used to store R_T
c=======================================================================
      do 10 j=2,jend
        dudy(j)=(uplus(j+1)-uplus(j-1))/(dyplus(j)+dyplus(j-1))
        d1=1.-exp(-yplus(j)/aplus0)
        d2=1.-exp(-yplus(j)/aplus2)
        f2(j)=ceps1/ceps2+(1.-ceps1/ceps2)*(1./(akappa*yplus(j))
     *       +d1*d2)*(sqrt(d1*d2)+yplus(j)/sqrt(d1*d2)
     *       *(d2/aplus0*exp(-yplus(j)/aplus0)
     *        +d1/aplus2*exp(-yplus(j)/aplus2)))
        eps(j)=cmu*et(j)*d1*d2
 10   continue
      dudy(1)=1.
      dudy(jmax)=0.
      eps(1)=0.
      eps(jmax)=eps(jend)
      epsm=1.+.5*(eps(1)+eps(2))/sigmae
      do 40 j=2,jend
        eps0=1.+eps(j)/sigmae
        epsp=1.+.5*(eps(j)+eps(j+1))/sigmae
        prode=(ceps2*f2(j)-ceps1)*dudy(j)*sqrt(eps(j)/etm(j))
        dti=dlti(j)
        a(j)=(2.*eps0-epsm)*y3(j)
        b(j)=-(2.*eps0-epsp)*y1(j)-(2.*eps0-epsm)*y3(j)+prode-2.*dti
        c(j)=(2.*eps0-epsp)*y1(j)
        et(j)=-a(j)*etm(j-1)-(b(j)+4.*dti)*etm(j)-c(j)*etm(j+1)
        epsm=epsp
 40   continue
      etl=etm(1)
      call tri(a,b,c,et,1,jmax,0,etl,1,0.)
      do 50 j=2,jmax
        if(et(j).lt.0.) et(j)=.1*etm(j)
 50   continue
      do 60 j=2,jend
        d1=1.-exp(-yplus(j)/aplus0)
        d2=1.-exp(-yplus(j)/aplus2)
        eps(j)=cmu*et(j)*d1*d2
 60   continue
      eps(1)=0.
      eps(jmax)=eps(jend)
      return
      end
      subroutine blomax
c=======================================================================
c      B A L D W I N - L O M A X   E D D Y   V I S C O S I T Y
c=======================================================================
      parameter (jmx=201)
      common/acom/ elmix(jmx),f2(jmx),fw(jmx),stilde(jmx)
      common/bcom/ dlti(jmx),dudy(jmx),dyplus(jmx),eps(jmx),et(jmx),
     *             etm(jmx),ret(jmx),uplus(jmx),uplusm(jmx),wt(jmx),
     *             wtm(jmx),y1(jmx),y2(jmx),y3(jmx),yoh(jmx),yplus(jmx)
c-----------------------------------------------------------------------
      common/data/     omegw  ,  retau  ,  yone
      common/extra/    akappa ,  cf     ,  cfl    ,  erro   ,  error  ,
     *                 jskipm ,  jskipp ,  nstep  ,  reh    ,  wn
      common/fctcom/   aa     ,  aam    ,  jend   ,  jendm
      common/mixlen/   aplus  ,  c1ml   ,  ccp    ,  ckleb  ,  clausr ,  
     *                 cwk
      common/ndata/    iruff  ,  jaxi   ,  jmax   ,  jskip  ,  maxn   ,
     *                 nedit  ,  nfreq  ,  nvisc
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c=======================================================================
      fmax=0.
      ymax=0.
      do 10 j=2,jend
        dudy(j)=(uplus(j+1)-uplus(j-1))/(dyplus(j)+dyplus(j-1))
        elmix(j)=akappa*yplus(j)*(1.-exp(-yplus(j)/aplus))
        f=elmix(j)*dudy(j)
        if(f.gt.fmax) then
          fmax=f
          ymax=yplus(j)
        endif
 10   continue
      dudy(1)=1.
      dudy(jmax)=0.
      fmax=fmax/akappa
      fwake1=ymax*fmax
      fwake2=cwk*ymax*uplus(jmax)**2/fmax
      fwake=dmin1(fwake1,fwake2)
      fwake=clausr*ccp*fwake
      iflag=0
      eps(1)=0.
      ypluss=ymax/ckleb
      do 20 j=2,jmax
        fkleb=1./(1.+5.5*(yplus(j)/ypluss)**6)
        eps(j)=fwake*fkleb
        if(iflag.eq.0) then
          epsi=elmix(j)**2*dudy(j)
          if(epsi.lt.eps(j)) then
            eps(j)=epsi
          else
            iflag=1
          endif
        endif
 20   continue
      return
      end
      subroutine calcs
c=======================================================================
c   C O M P U T E   Q U A N T I T I E S   N E E D E D   I N   E D D Y
c     K - O M E G A   A N D   S T R E S S - O M E G A   M O D E L S
c=======================================================================
      parameter (jmx=201)
      common/acom/ elmix(jmx),f2(jmx),fw(jmx),stilde(jmx)
      common/bcom/ dlti(jmx),dudy(jmx),dyplus(jmx),eps(jmx),et(jmx),
     *             etm(jmx),ret(jmx),uplus(jmx),uplusm(jmx),wt(jmx),
     *             wtm(jmx),y1(jmx),y2(jmx),y3(jmx),yoh(jmx),yplus(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx),vextra(jmx)
      common/vcom/ bbeta(jmx),climit(jmx),fc1(jmx),fgamma(jmx),
     *             fgammas(jmx),gamma(jmx),gammas(jmx)
c-----------------------------------------------------------------------
      common/data/     omegw  ,  retau  ,  yone
      common/extra/    akappa ,  cf     ,  cfl    ,  erro   ,  error  ,
     *                 jskipm ,  jskipp ,  nstep  ,  reh    ,  wn
      common/fctcom/   aa     ,  aam    ,  jend   ,  jendm
      common/iodef/    iunit1 ,  iunit2 ,  model
      common/ndata/    iruff  ,  jaxi   ,  jmax   ,  jskip  ,  maxn   ,
     *                 nedit  ,  nfreq  ,  nvisc
      common/stress/   alphah ,  betah  ,  c1     ,  gammah ,  onemab
      common/stres2/   cf1    ,  cf2    ,  cf3    ,  cfx    ,  cfy
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
      common/viscus/   alf0   ,  alfs0  ,  bbeta0 ,  gamh0  ,  rb     ,
     *                 rk     ,  rw
c=======================================================================
c    Stress-omega model
c=============================
      if(model.eq.1) then
        eps(1)=0.
        dudy(1)=1.
        climit(1)=1.
        dudy(jmax)=0.
        do 10 j=2,jmax
          climit(j)=1.
          ret(j)=etm(j)/wtm(j)
          if(j.ne.jmax) then
            dudy(j)=(uplusm(j+1)-uplusm(j-1))/(dyplus(j)+dyplus(j-1))
          endif
 10     continue
        if(nvisc.eq.0) then
          bbeta(1)=1.
          fc1(1)=1.
          gammas(1)=1.
          gamma(1)=1.
          fgamma(1)=1.
          fgammas(1)=1.
        else
          bbeta(1)=bbeta0
          fc1(1)=5./3.
          gammas(1)=alfs0
          gamma(1)=alf0/alfs0
          fgamma(1)=0.
          fgammas(1)=gamh0
        endif
        do 20 j=2,jmax
          gammas(j)=(rk*alfs0+ret(j))/(rk+ret(j))
          gamma(j)=(rw*alf0+ret(j))/(rw+ret(j))
     *            /(3.*rw*alfs0+ret(j))*(3.*rw+ret(j))
          ret4=ret(j)**4
          rb4=rb**4
          fc1(j)=(5./3.*rb4+ret4)/(rb4+ret4)
          fgamma(j)=ret4/(rb4+ret4)
          fgammas(j)=(gamh0*rb4+ret4)/(rb4+ret4)
          bbeta(j)=(bbeta0*rb4+ret4)/(rb4+ret4)
          eps(j)=gammas(j)*ret(j)
 20     continue
c=============================
c       k-omega model
c=============================
      else
        if(nvisc.eq.0) then
          bbeta(1)=1.
          gammas(1)=1.
          gamma(1)=1.
        else
          bbeta(1)=bbeta0
          gammas(1)=alfs0
          gamma(1)=alf0
        endif
        do 30 j=1,jmax
          ret(j)=etm(j)/wtm(j)
          if(j.ne.1) then
            gammas(j)=(alfs0*rk+ret(j))/(rk+ret(j))
            gamma(j)=(alf0*rw+ret(j))/(rw+ret(j))
            if(nvisc.ne.0) then
              bbeta(j)=(bbeta0+(ret(j)/rb)**4)/(1.+(ret(j)/rb)**4)
            else
              bbeta(j)=1.
            endif
          endif
          eps(j)=gammas(j)*ret(j)
 30     continue
c=============================
c  Stress-limiter correction
c=============================
        climit(1)=1.
        dudy(1)=(1.-yoh(1))/(1.+climit(1)*eps(1))
        do 35 j=2,jend
          pod=clim*sqrt(gammas(j))*dudy(j)/(alphas*wt(j))
          if(pod.gt.1.) then
            climit(j)=1./pod
            gamma(j)=gamma(j)*climit(j)
            gammas(j)=gammas(j)*climit(j)
          else
            climit(j)=1.
          endif
          dudy(j)=(1.-yoh(j))/(1.+climit(j)*eps(j))
 35     continue
        climit(jmax)=climit(jend)
        dudy(jmax)=(1.-yoh(jmax))/(1.+climit(jmax)*eps(jmax))
      endif
c=============================
c  Cross-diffusion correction
c=============================
      do 40 j=2,jend
        dkdy=(et(j+1)-et(j-1))/(dyplus(j-1)+dyplus(j))
        dwdy=(wt(j+1)-wt(j-1))/(dyplus(j-1)+dyplus(j))
        chik=dkdy*dwdy/wt(j)**3
        if(chik.gt.0.) then
          vextra(j)=sigmad*dkdy/wt(j)
        else
          vextra(j)=0.
        endif
 40   continue
      return
      end
      subroutine csmith
c=======================================================================
c      C E B E C I - S M I T H   E D D Y   V I S C O S I T Y
c=======================================================================
      parameter (jmx=201)
      common/acom/ elmix(jmx),f2(jmx),fw(jmx),stilde(jmx)
      common/bcom/ dlti(jmx),dudy(jmx),dyplus(jmx),eps(jmx),et(jmx),
     *             etm(jmx),ret(jmx),uplus(jmx),uplusm(jmx),wt(jmx),
     *             wtm(jmx),y1(jmx),y2(jmx),y3(jmx),yoh(jmx),yplus(jmx)
c-----------------------------------------------------------------------
      common/data/     omegw  ,  retau  ,  yone
      common/extra/    akappa ,  cf     ,  cfl    ,  erro   ,  error  ,
     *                 jskipm ,  jskipp ,  nstep  ,  reh    ,  wn
      common/fctcom/   aa     ,  aam    ,  jend   ,  jendm
      common/mixlen/   aplus  ,  c1ml   ,  ccp    ,  ckleb  ,  clausr ,  
     *                 cwk
      common/ndata/    iruff  ,  jaxi   ,  jmax   ,  jskip  ,  maxn   ,
     *                 nedit  ,  nfreq  ,  nvisc
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c=======================================================================
      duo=1.
      deltas=0.
      do 10 j=2,jmax
        dun=1.-uplus(j)/uplus(jmax)
        if(jaxi.ne.0) dun=dun*(1.-yoh(j))
        deltas=deltas+.5*(duo+dun)*dyplus(j)
        duo=dun
        if(j.ne.jmax) then
          dudy(j)=(uplus(j+1)-uplus(j-1))/(dyplus(j)+dyplus(j-1))
        endif
 10   continue
      if(jaxi.ne.0) deltas=sqrt(2.*retau*deltas)
      dudy(1)=1.
      dudy(jmax)=0.
      fwake=clausr*uplus(jmax)*deltas
      iflag=0
      eps(1)=0.
      do 20 j=2,jmax
        fkleb=1./(1.+5.5*yoh(j)**6)
        eps(j)=fwake*fkleb
        if(iflag.eq.0) then
          elmix0=akappa*yplus(j)*(1.-exp(-yplus(j)/aplus))
          epsi=elmix0**2*dudy(j)
          if(epsi.lt.eps(j)) then
            eps(j)=epsi
          else
            iflag=1
          endif
        endif
 20   continue
      return
      end
      subroutine eddy
c=======================================================================
c      A D V A N C E   T H E   S O L U T I O N   I N   T I M E
c=======================================================================
      parameter (jmx=201)
      dimension    a(jmx),b(jmx),c(jmx)
      dimension    axx(jmx),bxx(jmx),cxx(jmx)
      dimension    ayy(jmx),byy(jmx),cyy(jmx)
      common/acom/ elmix(jmx),f2(jmx),fw(jmx),stilde(jmx)
      common/bcom/ dlti(jmx),dudy(jmx),dyplus(jmx),eps(jmx),et(jmx),
     *             etm(jmx),ret(jmx),uplus(jmx),uplusm(jmx),wt(jmx),
     *             wtm(jmx),y1(jmx),y2(jmx),y3(jmx),yoh(jmx),yplus(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx),vextra(jmx)
      common/vcom/ bbeta(jmx),climit(jmx),fc1(jmx),fgamma(jmx),
     *             fgammas(jmx),gamma(jmx),gammas(jmx)
c-----------------------------------------------------------------------
      common/data/     omegw  ,  retau  ,  yone
      common/extra/    akappa ,  cf     ,  cfl    ,  erro   ,  error  ,
     *                 jskipm ,  jskipp ,  nstep  ,  reh    ,  wn
      common/fctcom/   aa     ,  aam    ,  jend   ,  jendm
      common/iodef/    iunit1 ,  iunit2 ,  model
      common/ndata/    iruff  ,  jaxi   ,  jmax   ,  jskip  ,  maxn   ,
     *                 nedit  ,  nfreq  ,  nvisc
      common/stress/   alphah ,  betah  ,  c1     ,  gammah ,  onemab
      common/stres2/   cf1    ,  cf2    ,  cf3    ,  cfx    ,  cfy
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
      common/viscus/   alf0   ,  alfs0  ,  bbeta0 ,  gamh0  ,  rb     ,
     *                 rk     ,  rw
      data nnmax/4/
c=======================================================================
c       Dissipation-rate equation
c=======================================
      if(model.ge.3.and.model.le.7) go to 110
      epsm=1.+.5*sigma*(eps(jskip)+eps(jskipp))
      do 10 j=jskipp,jend
        epsp=1.+.5*sigma*(eps(j)+eps(j+1))
        if(model.eq.0) then
          prodw=alpha*gamma(j)*dudy(j)**2/wtm(j)-beta*wtm(j)
        else
          prodw=alpha*gamma(j)*taum(j)*dudy(j)/etm(j)-beta*wtm(j)
        endif
        dti=dlti(j)
        if(abs(prodw).gt.1.) dti=dti*abs(prodw)
        a(j)=epsm*y3(j)-vextra(j)*y2(j)
        b(j)=-epsp*y1(j)-epsm*y3(j)+prodw-2.*dti
        c(j)=epsp*y1(j)+vextra(j)*y2(j)
        wt(j)=-a(j)*wtm(j-1)-(b(j)+4.*dti)*wtm(j)-c(j)*wtm(j+1)
        epsm=epsp
 10   continue
      wtl=wtm(jskip)
      call tri(a,b,c,wt,jskip,jmax,0,wtl,1,0.)
      do 20 j=1,jskipm
        wt(j)=wtm(j)
 20   continue
      do 30 j=1,jmax
        if(wt(j).lt.0.) wt(j)=.1*wtm(j)
 30   continue
c=======================================
c      Turbulence-energy equation
c=======================================
      epsm=1.+.5*sigmas*(eps(1)+eps(2))
      do 40 j=2,jend
        epsp=1.+.5*sigmas*(eps(j)+eps(j+1))
        if(model.eq.0) then
          prode=gammas(j)*dudy(j)**2/wt(j)-betas*bbeta(j)*wt(j)
        else
          prode=taum(j)*dudy(j)/etm(j)-betas*bbeta(j)*wt(j)
        endif
        dti=dlti(j)
        if(abs(prode).gt.1.) dti=dti*abs(prode)
        a(j)=epsm*y3(j)
        b(j)=-epsp*y1(j)-epsm*y3(j)+prode-2.*dti
        c(j)=epsp*y1(j)
        et(j)=-a(j)*etm(j-1)-(b(j)+4.*dti)*etm(j)-c(j)*etm(j+1)
        epsm=epsp
 40   continue
      call tri(a,b,c,et,1,jmax,0,0.,1,0.)
      if(model.eq.0) then
        eps(1)=0.
        do 50 j=2,jmax
          if(et(j).lt.0.) et(j)=.1*etm(j)
          ret(j)=et(j)/wt(j)
          gammas(j)=(alfs0*rk+ret(j))/(rk+ret(j))
          eps(j)=gammas(j)*ret(j)
 50     continue
c=============================
c  Stress-limiter correction
c=============================
        climit(1)=1.
        do 55 j=2,jend
          pod=clim*sqrt(gammas(j))*dudy(j)/(alphas*wt(j))
          if(pod.gt.1.) then
            climit(j)=1./pod
            gamma(j)=gamma(j)*climit(j)
            gammas(j)=gammas(j)*climit(j)
          else
            climit(j)=1.
          endif
 55     continue
        climit(jmax)=climit(jend)
        go to 110
      endif
c=======================================
c     Reynolds-stress equations
c=======================================
      do 90 nn=1,nnmax
        epsm=1.+.5*sigmas*(eps(1)+eps(2))
        do 60 j=2,jend
          epsp=1.+.5*sigmas*(eps(j)+eps(j+1))
          prodt=-c1*fc1(j)*betas*bbeta(j)*wt(j)
          dti=dlti(j)/nnmax
          if(abs(prodt).gt.1.) dti=dti*abs(prodt)
          a(j)=epsm*y3(j)
          axx(j)=a(j)
          ayy(j)=a(j)
          b(j)=-epsp*y1(j)-epsm*y3(j)+prodt-2.*dti
          bxx(j)=b(j)
          byy(j)=b(j)
          c(j)=epsp*y1(j)
          cxx(j)=c(j)
          cyy(j)=c(j)
          if(jaxi.eq.1) then
            b(j)=b(j)-(1.+sigmas*eps(j))/(retau-yplus(j))**2
            byy(j)=byy(j)-4.*(1.+sigmas*eps(j))/(retau-yplus(j))**2
          endif
          sigx(j)=-axx(j)*sigxm(j-1)-(bxx(j)+4.*dti)*sigxm(j)
     *            -cxx(j)*sigxm(j+1)-2.*cfx*fgamma(j)*taum(j)*dudy(j)
          tau(j)=-a(j)*taum(j-1)-(b(j)+4.*dti)*taum(j)-c(j)*taum(j+1)
     *           -2.*((cf1*sigxm(j)+cf2*sigym(j))*fgamma(j)
     *           +(cf3*fgamma(j)+.5*gammah*fgammas(j))*et(j))*dudy(j)
          sigy(j)=-ayy(j)*sigym(j-1)-(byy(j)+4.*dti)*sigym(j)
     *            -cyy(j)*sigym(j+1)-2.*cfy*fgamma(j)*taum(j)*dudy(j)
          if(jaxi.eq.1) then
            sigy(j)=sigy(j)+4.*(1.+sigmas*eps(j))*sigxm(j)
     *              /(retau-yplus(j))**2
          endif
          epsm=epsp
 60     continue
        call tri(a,b,c,tau,1,jmax,0,0.,0,0.)
        call tri(axx,bxx,cxx,sigx,1,jmax,0,0.,1,0.)
        call tri(ayy,byy,cyy,sigy,1,jmax,0,0.,1,0.)
        do 70 j=2,jend
          if(tau(j).lt.0.) tau(j)=.5*(tau(j-1)+taum(j+1))
          upup=sigx(j)+.6666667*et(j)
          if(upup.lt.0.) sigx(j)=-.6666667*et(j)
          if(upup.gt.2.*et(j)) sigx(j)=.6666667*et(j)
          vpvp=sigy(j)+.6666667*et(j)
          if(vpvp.lt.0.) sigy(j)=-.6666667*et(j)
          if(vpvp.gt.2.*et(j)) sigy(j)=.6666667*et(j)
 70     continue
        do 80 j=1,jmax
          sigxm(j)=sigx(j)
          sigym(j)=sigy(j)
          taum(j)=tau(j)
 80     continue
 90   continue
      do 100 j=1,jmax
        dtau(j)=tau(j)-eps(j)*dudy(j)
 100  continue
c=======================================
c             Momentum
c=======================================
 110  uplus(1)=0.
      dtauo=0.
      dtauh=0.
      dtaun=0.
      epso=1.+climit(1)*eps(1)
      epsio=1./epso
      do 120 j=2,jmax
        epsn=1.+climit(j)*eps(j)
        epsin=1./epsn
        epsion=2./(epso+epsn)
        if(model.ne.0) then
          dtaun=dtau(j)
          dtauh=.5*(dtauo+dtaun)
        endif
        ddy=dyplus(j-1)
        fk1=ddy*(1.-yoh(j)-dtauo)*epsio
        fk2=ddy*(1.-yoh(j)-dtauh)*epsion
        fk3=ddy*(1.-yoh(j)-dtauh)*epsion
        fk4=ddy*(1.-yoh(j)-dtaun)*epsin
        uplus(j)=uplus(j-1)+.16666667*(fk1+2.*(fk2+fk3)+fk4)
        epso=epsn
        epsio=epsin
        if(model.ne.0) dtauo=dtaun
 120  continue
      uavg=0.
      uo=0.
      do 130 j=2,jmax
        un=uplus(j)
        if(jaxi.ne.0) un=un*(1.-yoh(j))
        uavg=uavg+.5*(uo+un)*dyplus(j-1)
        uo=un
 130  continue
      uavg=(1.+jaxi)*uavg/retau
      reh=2.*retau*uavg
      cf=2./uavg**2
      return
      end
      subroutine edit
c=======================================================================
c                   P R I N T I N G   R O U T I N E
c=======================================================================
      parameter (jmx=201)
      character*1  width(2)
      common/acom/ elmix(jmx),f2(jmx),fw(jmx),stilde(jmx)
      common/bcom/ dlti(jmx),dudy(jmx),dyplus(jmx),eps(jmx),et(jmx),
     *             etm(jmx),ret(jmx),uplus(jmx),uplusm(jmx),wt(jmx),
     *             wtm(jmx),y1(jmx),y2(jmx),y3(jmx),yoh(jmx),yplus(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx),vextra(jmx)
      common/ucom/ uoum(jmx),upup(jmx),vpvp(jmx)
      common/vcom/ bbeta(jmx),climit(jmx),fc1(jmx),fgamma(jmx),
     *             fgammas(jmx),gamma(jmx),gammas(jmx)
c-----------------------------------------------------------------------
      common/data/     omegw  ,  retau  ,  yone
      common/extra/    akappa ,  cf     ,  cfl    ,  erro   ,  error  ,
     *                 jskipm ,  jskipp ,  nstep  ,  reh    ,  wn
      common/fctcom/   aa     ,  aam    ,  jend   ,  jendm
      common/iodef/    iunit1 ,  iunit2 ,  model
      common/ndata/    iruff  ,  jaxi   ,  jmax   ,  jskip  ,  maxn   ,
     *                 nedit  ,  nfreq  ,  nvisc
      common/stress/   alphah ,  betah  ,  c1     ,  gammah ,  onemab
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
      data width/'h','d'/
c=======================================================================
      write(iunit1,100) nstep
      write(iunit1,110) width(jaxi+1),reh,cf
      do 10 j=1,jmax
        uoum(j)=uplus(j)/uplus(jmax)
        if(model.ne.1) tau(j)=eps(j)*dudy(j)
 10   continue
c-----------------------------
c      One-equation and
c      algebraic models
c-----------------------------
      if(model.ge.3.and.model.le.7) then
        write(iunit1,120) width(jaxi+1)
        write(iunit1,130) (j,yoh(j),uoum(j),yplus(j),uplus(j),eps(j),
     *                     tau(j),j=1,jmax)
c-----------------------------
c        k-omega model
c-----------------------------
      elseif(model.eq.0) then
        write(iunit1,140) width(jaxi+1)
        write(iunit1,150) (j,yoh(j),uoum(j),yplus(j),uplus(j),
     *        climit(j)*eps(j),climit(j)*tau(j),et(j),wt(j),j=1,jmax)
c-----------------------------
c      Stress-omega model
c-----------------------------
      else
        do 20 j=1,jmax
          e23=.66666667*et(j)
          upup(j)=e23+sigx(j)
          vpvp(j)=e23+sigy(j)
 20     continue
        write(iunit1,160) width(jaxi+1)
        write(iunit1,170) (j,yoh(j),uoum(j),yplus(j),uplus(j),eps(j),
     *                  tau(j),upup(j),vpvp(j),et(j),wt(j),j=1,jmax)
      endif
c=======================================================================
 100  format(//,4x,'nstep =',i5)
 110  format(//4x,'***** Re',a1,' =',1pe11.3/4x,'***** Cf  =',e11.3//)
 120  format(/5x,'j',2x,'y/(',a1,'/2)',4x,'u/Um',6x,'yplus',5x,'uplus',
     *        6x,'nut',6x,5h-u'v')
 130  format(i6,1p6e10.3)
 140  format(/5x,'j',2x,'y/(',a1,'/2)',4x,'u/Um',6x,'yplus',5x,'uplus',
     *        6x,'nut',6x,5h-u'v',3x,'k/utau**2',3x,'omega+')
 150  format(i6,1p8e10.3)
 160  format(/5x,'j',2x,'y/(',a1,'/2)',4x,'u/Um',6x,'yplus',5x,'uplus',
     *        6x,'nut',6x,5h-u'v',6x,4hu'u',6x,4hv'v',3x,'k/utau**2',
     *        3x,'omega+')
 170  format(i6,1p10e10.3)
c=======================================================================
      return
      end
      subroutine fct(x,f,fp)
c=======================================================================
c     F I N I T E   D I F F E R E N C E   G R I D   F U N C T I O N
c=======================================================================
      common/fctcom/   aa     ,  aam    ,  jend   ,  jendm
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
      subroutine grid
c=======================================================================
c  G E N E R A T E   T H E   F I N I T E   D I F F E R E N C E   G R I D
c=======================================================================
      external fct
      parameter (jmx=201)
      common/acom/ elmix(jmx),f2(jmx),fw(jmx),stilde(jmx)
      common/bcom/ dlti(jmx),dudy(jmx),dyplus(jmx),eps(jmx),et(jmx),
     *             etm(jmx),ret(jmx),uplus(jmx),uplusm(jmx),wt(jmx),
     *             wtm(jmx),y1(jmx),y2(jmx),y3(jmx),yoh(jmx),yplus(jmx)
c-----------------------------------------------------------------------
      common/data/     omegw  ,  retau  ,  yone
      common/extra/    akappa ,  cf     ,  cfl    ,  erro   ,  error  ,
     *                 jskipm ,  jskipp ,  nstep  ,  reh    ,  wn
      common/fctcom/   aa     ,  aam    ,  jend   ,  jendm
      common/iodef/    iunit1 ,  iunit2 ,  model
      common/ndata/    iruff  ,  jaxi   ,  jmax   ,  jskip  ,  maxn   ,
     *                 nedit  ,  nfreq  ,  nvisc
      common/stress/   alphah ,  betah  ,  c1     ,  gammah ,  onemab
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c=======================================================================
c First solve for RATE
c=====================
      ratei=1.+10./jmax
      call rtni(rate,f,fp,fct,ratei,1.e-6,500,ier)
      if(ier.eq.1) then
        write(*,*) 'Convergence not attained in computing RATE'
        stop
      elseif(ier.eq.2) then
        write(*,*) 'df/dx = 0 encountered in computing RATE'
        stop
      endif
c====================
c Now set up the grid
c====================
      yplus(1)=0.
      yoh(1)=0.
      dyplus(1)=yone
      do 10 j=2,jend
        dyplus(j)=rate*dyplus(j-1)
        yplus(j)=yplus(j-1)+dyplus(j-1)
        yoh(j)=yplus(j)/retau
 10   continue
      dyplus(jmax)=dyplus(jend)
      yplus(jmax)=retau
      yoh(jmax)=1.
      do 20 j=2,jend
        denom=dyplus(j)+dyplus(j-1)
        y1(j)=2./(dyplus(j)*denom)
        y2(j)=1./denom
        y3(j)=2./(dyplus(j-1)*denom)
        if(jaxi.ne.0) then
          y1(j)=.5*y1(j)*(2.-yoh(j+1)-yoh(j))/(1.-yoh(j))
          y3(j)=.5*y3(j)*(2.-yoh(j-1)-yoh(j))/(1.-yoh(j))
        endif
 20   continue
      return
      end
      subroutine jking
c=======================================================================
c     J O H N S O N - K I N G   1 / 2 - E Q U A T I O N   M O D E L
c=======================================================================
      parameter (jmx=201)
      common/acom/ elmix(jmx),f2(jmx),fw(jmx),stilde(jmx)
      common/bcom/ dlti(jmx),dudy(jmx),dyplus(jmx),eps(jmx),et(jmx),
     *             etm(jmx),ret(jmx),uplus(jmx),uplusm(jmx),wt(jmx),
     *             wtm(jmx),y1(jmx),y2(jmx),y3(jmx),yoh(jmx),yplus(jmx)
c-----------------------------------------------------------------------
      common/data/     omegw  ,  retau  ,  yone
      common/extra/    akappa ,  cf     ,  cfl    ,  erro   ,  error  ,
     *                 jskipm ,  jskipp ,  nstep  ,  reh    ,  wn
      common/fctcom/   aa     ,  aam    ,  jend   ,  jendm
      common/mixlen/   aplus  ,  c1ml   ,  ccp    ,  ckleb  ,  clausr ,  
     *                 cwk
      common/ndata/    iruff  ,  jaxi   ,  jmax   ,  jskip  ,  maxn   ,
     *                 nedit  ,  nfreq  ,  nvisc
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c=======================================================================
      duo=1.
      deltas=0.
      do 10 j=2,jmax
        dun=1.-uplus(j)/uplus(jmax)
        if(jaxi.ne.0) dun=dun*(1.-yoh(j))
        deltas=deltas+.5*(duo+dun)*dyplus(j)
        duo=dun
        if(j.ne.jmax) then
          dudy(j)=(uplus(j+1)-uplus(j-1))/(dyplus(j)+dyplus(j-1))
        endif
 10   continue
      if(jaxi.ne.0) deltas=sqrt(2.*retau*deltas)
      dudy(1)=1.
      dudy(jmax)=0.
      taumax=0.
      do 20 j=2,jmax
        tauloc=eps(j)*dudy(j)
        if(tauloc.gt.taumax) then
          taumax=tauloc
          yplusm=yplus(j)
        else
          go to 30
        endif
 20   continue
 30   if(taumax.ne.0.) then
        umplus=sqrt(taumax)
      else
        umplus=1.
        yplusm=10.
      endif
      rlen=c1ml*retau/akappa
      if(yplusm.le.rlen) then
        elm=akappa*yplusm
      else
        elm=c1ml*retau
      endif
      elc=elm/(1.+umplus)
      fwake=clausr*uplus(jmax)*deltas
      eps(1)=0.
      do 40 j=2,jmax
        fkleb=1./(1.+5.5*yoh(j)**6)
        epso=fwake*fkleb
        gamma2=xtanh(yplus(j)/elc)
        usplus=1.-gamma2*(1.-umplus)
        epsi=akappa*usplus*yplus(j)*(1.-exp(-yplus(j)/aplus))**2
        eps(j)=epso*xtanh(epsi/epso)
 40   continue
      return
      end
      subroutine namin(icpu)
c=======================================================================
c            COORDINATE READING INPUT DATA FILE 'PIPE.DAT'
c=======================================================================
      parameter (niu=3,nn=8,nxd=3)
      character a*12,newfil*7,pform*9
      common/data/  xd(nxd)
      common/iodef/ iu(niu)
      common/ndata/ n(nn)
c=======================================================================
c----------------------------------
c         Read I/O units
c----------------------------------
      call namsys(icpu,iin,iv,msdos,newfil,pform)
      open(iin,file='pipe.dat',status='old')
      iu(1)=2
      iu(2)=7
      io=iu(1)
      open(io,file='pipe.prt',form=pform,status=newfil)
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
      open(iu(2),file='plotp.dat',form=pform,status=newfil)
      return
c=======================================================================
 100  format(1x,a12,i4)
 200  format(1x,'Program PIPE input data:'/,
     *       1x,'________________________'/)
 300  format(1x,a12,1pe13.6)
 400  format(1x,'________________________'/)
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
      subroutine rtni(x,f,fp,fct,xst,eps,iend,ier)
c=======================================================================
c      SOLVE  f(x) = 0  USING NEWTON'S METHOD.  SINGLE PRECISION.
c-----------------------------------------------------------------------
      external fct
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
        a=abs(x)
        if(a.gt.1.) tol=tol*a
        if(abs(dx).le.tol) then
          if(abs(f).le.tolf) go to 20
        endif
 10   continue
c============================================
c  If no convergence after iend iterations...
c  return with ier set to 1
c============================================
      ier=1
 20   return
      end
      subroutine spalart
c=======================================================================
c            S P A L A R T - A L L M A R A S   M O D E L
c=======================================================================
      parameter (jmx=201)
      dimension    a(jmx),b(jmx),c(jmx)
      common/acom/ elmix(jmx),f2(jmx),fw(jmx),stilde(jmx)
      common/bcom/ dlti(jmx),dudy(jmx),dyplus(jmx),eps(jmx),et(jmx),
     *             etm(jmx),ret(jmx),uplus(jmx),uplusm(jmx),wt(jmx),
     *             wtm(jmx),y1(jmx),y2(jmx),y3(jmx),yoh(jmx),yplus(jmx)
c-----------------------------------------------------------------------
      common/data/     omegw  ,  retau  ,  yone
      common/extra/    akappa ,  cf     ,  cfl    ,  erro   ,  error  ,
     *                 jskipm ,  jskipp ,  nstep  ,  reh    ,  wn
      common/fctcom/   aa     ,  aam    ,  jend   ,  jendm
      common/mixlen/   aplus  ,  c1ml   ,  ccp    ,  ckleb  ,  clausr ,  
     *                 cwk
      common/ndata/    iruff  ,  jaxi   ,  jmax   ,  jskip  ,  maxn   ,
     *                 nedit  ,  nfreq  ,  nvisc
      common/oneeqn/   aplus0 ,  aplus2 ,  cb1    ,  cb2    ,  ceps1  ,  
     *                 ceps2  ,  cmu    ,  cv1    ,  cw1    ,  cw2    ,
     *                 cw3    ,  sigmae
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c=======================================================================
c  NOTE: The ET and ETM arrays are used to store \tilde{nu}
c=======================================================================
      do 10 j=2,jend
        dudy(j)=(uplus(j+1)-uplus(j-1))/(dyplus(j)+dyplus(j-1))
        fv1=et(j)**3/(et(j)**3+cv1**3)
        fv2=1.-et(j)/(1.+et(j)*fv1)
        stilde(j)=dudy(j)+et(j)*fv2/(akappa*yplus(j))**2
        rr=et(j)/(stilde(j)*(akappa*yplus(j))**2)
        gg=rr+cw2*(rr**6-rr)
        fw(j)=gg*((1.+cw3**6)/(gg**6+cw3**6))**(1/6)
 10   continue
      dudy(1)=1.
      dudy(jmax)=0.
      eps(1)=0.
      eps(jmax)=eps(jend)
      epsm=(1.+.5*(1.+cb2)*(eps(1)+eps(2)))/sigmae
      do 40 j=2,jend
        eps0=cb2*eps(j)/sigmae
        epsp=(1.+.5*(1.+cb2)*(eps(j)+eps(j+1)))/sigmae
        prode=cb1*stilde(j)-cw1*fw(j)*etm(j)/yplus(j)**2
        dti=dlti(j)
        a(j)=(epsm-eps0)*y3(j)
        b(j)=-(epsp-eps0)*y1(j)-(epsm-eps0)*y3(j)-2.*dti
        c(j)=(epsp-eps0)*y1(j)
        et(j)=-a(j)*etm(j-1)-(b(j)+4.*dti)*etm(j)-c(j)*etm(j+1)
        if(prode.gt.0.) then
          et(j)=et(j)-prode*etm(j)
        else
          b(j)=b(j)+prode
        endif
        epsm=epsp
 40   continue
      etl=etm(1)
      call tri(a,b,c,et,1,jmax,0,etl,1,0.)
      do 50 j=2,jmax
        if(et(j).lt.0.) et(j)=.1*etm(j)
 50   continue
      do 60 j=2,jend
        fv1=et(j)**3/(et(j)**3+cv1**3)
        eps(j)=fv1*et(j)
 60   continue
      eps(1)=0.
      eps(jmax)=eps(jend)
      return
      end
      subroutine start
c=======================================================================
c            S E T   U P   I N I T I A L   P R O F I L E S
c=======================================================================
      parameter (jmx=201)
      common/acom/ elmix(jmx),f2(jmx),fw(jmx),stilde(jmx)
      common/bcom/ dlti(jmx),dudy(jmx),dyplus(jmx),eps(jmx),et(jmx),
     *             etm(jmx),ret(jmx),uplus(jmx),uplusm(jmx),wt(jmx),
     *             wtm(jmx),y1(jmx),y2(jmx),y3(jmx),yoh(jmx),yplus(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx),vextra(jmx)
      common/vcom/ bbeta(jmx),climit(jmx),fc1(jmx),fgamma(jmx),
     *             fgammas(jmx),gamma(jmx),gammas(jmx)
c-----------------------------------------------------------------------
      common/data/     omegw  ,  retau  ,  yone
      common/extra/    akappa ,  cf     ,  cfl    ,  erro   ,  error  ,
     *                 jskipm ,  jskipp ,  nstep  ,  reh    ,  wn
      common/fctcom/   aa     ,  aam    ,  jend   ,  jendm
      common/iodef/    iunit1 ,  iunit2 ,  model
      common/mixlen/   aplus  ,  c1ml   ,  ccp    ,  ckleb  ,  clausr ,  
     *                 cwk
      common/ndata/    iruff  ,  jaxi   ,  jmax   ,  jskip  ,  maxn   ,
     *                 nedit  ,  nfreq  ,  nvisc
      common/oneeqn/   aplus0 ,  aplus2 ,  cb1    ,  cb2    ,  ceps1  ,  
     *                 ceps2  ,  cmu    ,  cv1    ,  cw1    ,  cw2    ,
     *                 cw3    ,  sigmae
      common/stress/   alphah ,  betah  ,  c1     ,  gammah ,  onemab
      common/stres2/   cf1    ,  cf2    ,  cf3    ,  cfx    ,  cfy
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
      common/viscus/   alf0   ,  alfs0  ,  bbeta0 ,  gamh0  ,  rb     ,
     *                 rk     ,  rw
c=======================================================================
c Set common k-omega/stress-omega and
c  high-Reynolds-number coefficients
c=======================================
      akappa=.41
      alpha=.52
      beta=.0708
      betas=.09
      sigma=.5
      sigmad=.125
      sigmas=.6
      clim=7./8.
      alfs0=1.
      alf0=1.
      gamh0=1.
      rb=0.
      rk=0.
      rw=0.
      alphas=sqrt(betas)
c=======================================
c     k-omega model coefficients
c=======================================
      if(model.eq.0.and.nvisc.ne.0) then
        alfs0=beta/3.
        alf0=1./9.
        bbeta0=beta/(3.*betas)
        rb=8.
        rk=6.
        rw=2.61
c=======================================
c   Stress-omega model coefficients
c=======================================
      elseif(model.eq.1) then
        c1=1.8
        c2=10./19.
        alphah=(8.+c2)/11.
        betah=(8.*c2-2.)/11.
        gammah=(60.*c2-4.)/55.
        if(nvisc.ne.0) then
          alfs0=beta/3.
          alf0=1./9.
          gamh0=21./2000.
          bbeta0=beta/(3.*betas)
          rb=8.
          rk=6.
          rw=22./9.
        endif
        onemal=1.-alphah
        onemab=onemal-betah
        cf1=-betah
        cf2=onemal
        cf3=.66666667*onemab
        cfx=.66666667*(2.*onemal+betah)
        cfy=-.66666667*(onemal+2.*betah)
c=======================================
c   Baldwin-Barth model coefficients
c=======================================
      elseif(model.eq.3) then
        aplus0=26.
        aplus2=10.
        ceps1=1.2
        ceps2=2.0
        cmu=.09
        sigmae=akappa**2/((ceps2-ceps1)*sqrt(cmu))
c=======================================
c  Spalart-Allmaras model coefficients
c=======================================
      elseif(model.eq.4) then
        cb1=.1355
        cb2=.622
        cv1=7.1
        sigmae=2./3.
        cw1=cb1/akappa**2+(1.+cb2)/sigmae
        cw2=.3
        cw3=2.
c=======================================
c   Cebeci-Smith model coefficients
c=======================================
      elseif(model.eq.5) then
        akappa=.4
        aplus=26.
        clausr=.0168
        wn=.65
c=======================================
c   Baldwin-Lomax model coefficients
c=======================================
      elseif(model.eq.6) then
        akappa=.4
        aplus=26.
        clausr=.0168
        ccp=1.6
        ckleb=.3
        cwk=1.
        wn=.4
c=======================================
c   Johnson-King model coefficients
c=======================================
      elseif(model.eq.7) then
        akappa=.4
        aplus=17.
        c1ml=.09
        clausr=.0168
        wn=.6
      endif
c=============================
c Initial profile coefficients
c=============================
      ymatch=10.
      c=5.
      if(iruff.eq.0) then
        omegi=1./omegw
      else
        omegi=0.
      endif
      do 10 j=1,jmax
        if(yplus(j).le.ymatch) then
          et(j)=(yplus(j)/ymatch)**2/alphas
          uplus(j)=yplus(j)
        else
          et(j)=1./alphas
          uplus(j)=dlog(yplus(j))/akappa+c
        endif
 10   continue
      con=sqrt(beta/6.)
      ymatch=6.*alphas*akappa/beta
      ww=sqrt(omegi)
      do 20 j=1,jmax
        if(yplus(j).le.ymatch) then
          yhat=con*yplus(j)+ww
          if(iruff.eq.1.and.j.eq.1) yhat=con*yplus(2)
          wt(j)=1./yhat**2
        else
          wt(j)=1./(alphas*akappa*yplus(j))
        endif
        if(wt(j).gt.omegw) wt(j)=omegw
 20   continue
      do 30 j=1,jmax
        bbeta(j)=1.
        uplusm(j)=uplus(j)
        etm(j)=et(j)
        wtm(j)=wt(j)
        eps(j)=etm(j)/wtm(j)
        dtau(j)=0.
        if(model.eq.1) then
          sigx(j)=.243*et(j)
          sigy(j)=-.239*et(j)
          tau(j)=.3*et(j)
          sigxm(j)=sigx(j)
          sigym(j)=sigy(j)
          taum(j)=tau(j)
        endif
        if(j.ne.1.and.j.ne.jmax)
     *     dudy(j)=(uplus(j+1)-uplus(j-1))/(dyplus(j)+dyplus(j-1))
        if(model.eq.3.or.model.eq.4) then
          et(j)=eps(j)
          etm(j)=eps(j)
        endif
 30   continue
      uavg=0.
      uo=0.
      dudy(1)=1.
      dudy(jmax)=0.
      climit(1)=1.
      do 40 j=2,jmax
        un=uplus(j)
        if(jaxi.ne.0) un=un*(1.-yoh(j))
        uavg=uavg+.5*(uo+un)*dyplus(j-1)
        uo=un
        climit(j)=1.
 40   continue
      uavg=(1.+jaxi)*uavg/retau
      reh=2.*retau*uavg
      cf=2./uavg**2
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
      function xtanh(x)
c=======================================================================
c              H Y P E R B O L I C   T A N G E N T
c=======================================================================
      if(x.eq.0.) then
        xtanh=0.
      else
        y=exp(-2.*x)
        xtanh=(1.-y)/(1.+y)
      endif
      return
      end
