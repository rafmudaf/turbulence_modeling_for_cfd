      program sublay
c=======================================================================
c            I N C O M P R E S S I B L E   S U B L A Y E R
c
c    Copyright (c) 2006 DCW Industries, Inc.  All rights reserved.
c-----------------------------------------------------------------------
c                 MODEL = 0.............k-omega Model
c                 MODEL = 1........Stress-omega Model
c                 MODEL = 99.......User-Defined Model
c=======================================================================
      parameter (jmx=401)
      common/bcom/ dlti(jmx),dudy(jmx),dyplus(jmx),eps(jmx),et(jmx),
     *             etm(jmx),ret(jmx),uplus(jmx),uplusm(jmx),wt(jmx),
     *             wtm(jmx),y1(jmx),y2(jmx),y3(jmx),yplus(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx),vextra(jmx)
      common/vcom/ bbeta(jmx),climit(jmx),fc1(jmx),fgamma(jmx),
     *             fgammas(jmx),gamma(jmx),gammas(jmx)
c-----------------------------------------------------------------------
      common/data/     omegw  ,  vwplus ,  ymax   ,  yone
      common/extra/    cfl    ,  ctilda ,  cwall  ,  cc0    ,  cc1    ,
     *                 erro   ,  error  ,  jskipm ,  jskipp ,  nstep  ,
     *                 rate   ,  umax
      common/fctcom/   aa     ,  aam    ,  jend   ,  jendm
      common/iodef/    iunit1 ,  iunit2 ,  model
      common/ndata/    iruff  ,  jmax   ,  jskip  ,  maxn   ,  nedit  ,
     *                 nfreq  ,  nvisc
      common/stress/   alphah ,  betah  ,  c1     ,  gammah ,  onemab
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
c     Calculate auxiliary constants
c=======================================
      jend=jmax-1
      jendm=jend-1
      jskipm=jskip-1
      jskipp=jskip+1
      if(abs(vwplus).gt.0.5) then
        cfl=2./(1.+9.*nvisc)
      elseif(abs(vwplus).gt.0.1) then
        cfl=2.
      else
        cfl=5.
      endif
      error=1.e-5
c=======================================
c  Generate the grid and timestep array
c=======================================
      aa=yone/ymax
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
      nstep=0
c=======================================
c   Identify flow and turbulence model
c=======================================
      write(iunit1,100)
      write(*,100)
      if(model.eq.0) then
        write(iunit1,110)
        write(*,110)
      elseif(model.eq.1) then
        write(iunit1,120)
        write(*,120)
      elseif(model.eq.99) then
        write(iunit1,130)
        write(*,130)
      endif
      if(iruff.eq.0) then
        write(*,140) omegw
      else
        write(*,150)
      endif
      write(iunit1,160)
      write(*,160)
c=======================================================================
c  	       	    BEGINNING OF MAIN INTEGRATION LOOP
c=======================================================================
 20   nstep=nstep+1
c=======================================
c  Calculate quantities needed in EDDY
c=======================================
      call calcs
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
        err=abs(uplus(j)-uplusm(j))
        if(err.ge.erro) then
          erro=err
        endif
 30   continue
c=======================================
c         Converged solution
c=======================================
      if(erro.le.error) then
        const=uplus(jmax)-umax
        write(iunit1,170) nstep,erro,const
        write(*,170) nstep,erro,const
        write(*,180)
c=======================================
c       Un-converged solution
c=======================================
      else
        if(mod(nstep,nfreq).eq.0) then
          const=uplus(jmax)-umax
          write(iunit1,170) nstep,erro,const
          write(*,170) nstep,erro,const
        endif
        if(nstep.lt.maxn) then
c=======================================
c            Print profiles
c=======================================
          if(mod(nstep,nedit).eq.0) call edit
c=======================================
c      Prepare for next timestep
c=======================================
          do 40 j=1,jmax
            etm(j)=et(j)
            uplusm(j)=uplus(j)
            wtm(j)=wt(j)
 40       continue
          go to 20
c=======================================================================
c  	       	     END OF MAIN INTEGRATION LOOP
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
      do 50 j=2,jmax
        wtm(j)=-betas*bbeta(j)*et(j)*wt(j)
        if(model.ne.1) tau(j)=eps(j)*dudy(j)
        etm(j)=tau(j)*dudy(j)
 50   continue
      if(nvisc.ne.0) then
        wtm(1)=2.*wtm(2)-wtm(3)
      else
        wtm(1)=0.
      endif
      etm(1)=0.
      write(iunit2,200) jmax,vwplus
      write(iunit2,210) (yplus(j),uplus(j),wtm(j),etm(j),j=1,jmax)
      endfile(iunit2)
      close(iunit2)
c=======================================
c    Pause to keep window visible after
c    the run completes...
c=======================================
      pause
c=======================================================================
 100  format(' VISCOUS SUBLAYER')
 110  format(' k-omega Model')
 120  format(' Stress-omega Model')
 130  format(' User-Defined Model')
 140  format(' omegw =',1pe10.2)
 150  format(' smooth wall')
 160  format(1x,'Convergence history:'/)
 170  format(1x,'At Iteration',i5,'...Max Error =',1pe10.3,
     *          '  C =',e10.3)
 180  format(/14x,'***** SOLUTION CONVERGED *****'/)
 190  format(/13x,'Computation terminated for max iteration',i5/)
 200  format(1x,i4,1pe14.6)
 210  format(1p4e14.6)
c=======================================================================
      end
      subroutine calcs
c=======================================================================
c   C O M P U T E   Q U A N T I T I E S   N E E D E D   I N   E D D Y
c=======================================================================
      parameter (jmx=401)
      common/bcom/ dlti(jmx),dudy(jmx),dyplus(jmx),eps(jmx),et(jmx),
     *             etm(jmx),ret(jmx),uplus(jmx),uplusm(jmx),wt(jmx),
     *             wtm(jmx),y1(jmx),y2(jmx),y3(jmx),yplus(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx),vextra(jmx)
      common/vcom/ bbeta(jmx),climit(jmx),fc1(jmx),fgamma(jmx),
     *             fgammas(jmx),gamma(jmx),gammas(jmx)
c-----------------------------------------------------------------------
      common/data/     omegw  ,  vwplus ,  ymax   ,  yone
      common/extra/    cfl    ,  ctilda ,  cwall  ,  cc0    ,  cc1    ,
     *                 erro   ,  error  ,  jskipm ,  jskipp ,  nstep  ,
     *                 rate   ,  umax
      common/fctcom/   aa     ,  aam    ,  jend   ,  jendm
      common/iodef/    iunit1 ,  iunit2 ,  model
      common/ndata/    iruff  ,  jmax   ,  jskip  ,  maxn   ,  nedit  ,
     *                 nfreq  ,  nvisc
      common/stress/   alphah ,  betah  ,  c1     ,  gammah ,  onemab
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
      common/viscus/   alf0   ,  alfs0  ,  bbeta0 ,  gamh0  ,  rb     ,
     *                 rk     ,  rw
c=======================================================================
c    Stress-omega model
c=============================
      if(model.eq.1) then
        ret(1)=0.
        eps(1)=0.
        dudy(1)=1.
        climit(1)=1.
        do 10 j=2,jmax
          climit(j)=1.
          ret(j)=etm(j)/wtm(j)
          if(j.ne.jmax) then
            dudy(j)=(uplusm(j+1)-uplusm(j-1))/(dyplus(j)+dyplus(j-1))
          else
            dudy(j)=2.*dudy(jend)-dudy(jend-1)
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
              retorb=(ret(j)/rb)**4
              bbeta(j)=(bbeta0+retorb)/(1.+retorb)
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
        dudy(1)=1.
c----------------------------
c For numerical stability...
c----------------------------
        if(abs(vwplus).gt.0.5.and.nvisc.ne.0.and.nstep.lt.200) then
          dlim=1./6.
        else
          dlim=clim
        endif
c----------------------------
        do 40 j=2,jend
          pod=dlim*sqrt(gammas(j))*dudy(j)/(alphas*wt(j))
          if(pod.gt.1.) then
            climit(j)=1./pod
            gamma(j)=gamma(j)*climit(j)
            gammas(j)=gammas(j)*climit(j)
          else
            climit(j)=1.
          endif
          dudy(j)=(1.+vwplus*uplusm(j))/(1.+climit(j)*eps(j))
 40     continue
        climit(jmax)=climit(jend)
        dudy(jmax)=dudy(jend)
        endif
c=============================
c  Cross-diffusion correction
c=============================
      do 50 j=2,jend
        dkdy=(et(j+1)-et(j-1))/(dyplus(j-1)+dyplus(j))
        dwdy=(wt(j+1)-wt(j-1))/(dyplus(j-1)+dyplus(j))
        chik=dkdy*dwdy/wt(j)**3
        if(chik.gt.0.) then
          vextra(j)=sigmad*dkdy/wt(j)
        else
          vextra(j)=0.
        endif
 50   continue
      return
      end
      subroutine eddy
c=======================================================================
c        A D V A N C E   T H E   S O L U T I O N   I N   T I M E
c=======================================================================
      parameter (jmx=401)
      dimension    a(jmx),b(jmx),c(jmx)
      dimension    axx(jmx),bxx(jmx),cxx(jmx)
      dimension    ayy(jmx),byy(jmx),cyy(jmx)
      common/bcom/ dlti(jmx),dudy(jmx),dyplus(jmx),eps(jmx),et(jmx),
     *             etm(jmx),ret(jmx),uplus(jmx),uplusm(jmx),wt(jmx),
     *             wtm(jmx),y1(jmx),y2(jmx),y3(jmx),yplus(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx),vextra(jmx)
      common/vcom/ bbeta(jmx),climit(jmx),fc1(jmx),fgamma(jmx),
     *             fgammas(jmx),gamma(jmx),gammas(jmx)
c-----------------------------------------------------------------------
      common/data/     omegw  ,  vwplus ,  ymax   ,  yone
      common/extra/    cfl    ,  ctilda ,  cwall  ,  cc0    ,  cc1    ,
     *                 erro   ,  error  ,  jskipm ,  jskipp ,  nstep  ,
     *                 rate   ,  umax
      common/fctcom/   aa     ,  aam    ,  jend   ,  jendm
      common/iodef/    iunit1 ,  iunit2 ,  model
      common/ndata/    iruff  ,  jmax   ,  jskip  ,  maxn   ,  nedit  ,
     *                 nfreq  ,  nvisc
      common/stress/   alphah ,  betah  ,  c1     ,  gammah ,  onemab
      common/stres2/   cf1    ,  cf2    ,  cf3    ,  cfx    ,  cfy
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
      common/viscus/   alf0   ,  alfs0  ,  bbeta0 ,  gamh0  ,  rb     ,
     *                 rk     ,  rw
      data nnmax/4/
c=======================================================================
c Dissipation-rate equation
c==========================
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
        a(j)=epsm*y3(j)+y2(j)-vextra(j)/(dyplus(j)+dyplus(j-1))
        b(j)=-epsp*y1(j)-epsm*y3(j)+prodw-2.*dti
        c(j)=epsp*y1(j)-y2(j)+vextra(j)/(dyplus(j)+dyplus(j-1))
        wt(j)=-a(j)*wtm(j-1)-(b(j)+4.*dti)*wtm(j)-c(j)*wtm(j+1)
        epsm=epsp
 10   continue
      wl=wtm(jskip)
      dlyi=1./((dlog(yplus(jmax))-dlog(yplus(jend)))*yplus(jmax))
      wu=(uplusm(jmax)-uplusm(jend))*dlyi/alphas
      if(abs(vwplus).lt..00001) wu=wtm(jmax)
      if(abs(wu).lt.1.e-10) wu=wtm(jmax)
      call tri(a,b,c,wt,jskip,jmax,0,wl,0,wu)
      do 20 j=1,jskipm
        wt(j)=wtm(j)
 20   continue
      do 30 j=1,jmax
        if(wt(j).lt.0.) wt(j)=.1*wtm(j)
 30   continue
c=============================
c  Turbulence-energy equation
c=============================
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
        a(j)=epsm*y3(j)+y2(j)
        b(j)=-epsp*y1(j)-epsm*y3(j)+prode-2.*dti
        c(j)=epsp*y1(j)-y2(j)
        et(j)=-a(j)*etm(j-1)-(b(j)+4.*dti)*etm(j)-c(j)*etm(j+1)
        epsm=epsp
 40   continue
      etl=etm(1)
      etu=(1.+vwplus*uplusm(jmax))/alphas
      call tri(a,b,c,et,1,jmax,0,etl,0,etu)
      if(model.eq.0) then
c----------------------------
c For numerical stability...
c----------------------------
        if(abs(vwplus).gt.0.5.and.nvisc.ne.0.and.nstep.lt.200) then
          dlim=1./6.
        else
          dlim=clim
        endif
c----------------------------
        eps(1)=0.
        do 50 j=2,jmax
          if(et(j).lt.0.) et(j)=.1*etm(j)
          ret(j)=et(j)/wt(j)
          gammas(j)=(alfs0*rk+ret(j))/(rk+ret(j))
          eps(j)=gammas(j)*ret(j)
          pod=dlim*sqrt(gammas(j))*dudy(j)/(alphas*wt(j))
          if(pod.gt.1.) then
            climit(j)=1./pod
          else
            climit(j)=1.
          endif
 50     continue
        climit(jmax)=climit(jend)
      endif
      if(model.eq.0) go to 110
c=============================
c  Reynolds-stress equations
c=============================
      do 90 nn=1,nnmax
        epsm=1.+.5*sigmas*(eps(1)+eps(2))
        do 60 j=2,jend
          epsp=1.+.5*sigmas*(eps(j)+eps(j+1))
          prodt=-c1*fc1(j)*betas*bbeta(j)*wt(j)
          dti=dlti(j)/nnmax
          if(abs(prodt).gt.1.) dti=dti*abs(prodt)
          a(j)=epsm*y3(j)+y2(j)
          axx(j)=a(j)
          ayy(j)=a(j)
          b(j)=-epsp*y1(j)-epsm*y3(j)+prodt-2.*dti
          bxx(j)=b(j)
          byy(j)=b(j)
          c(j)=epsp*y1(j)-y2(j)
          cxx(j)=c(j)
          cyy(j)=c(j)
          sigx(j)=-a(j)*sigxm(j-1)-(b(j)+4.*dti)*sigxm(j)
     *            -c(j)*sigxm(j+1)-2.*cfx*fgamma(j)*taum(j)*dudy(j)
          tau(j)=-a(j)*taum(j-1)-(b(j)+4.*dti)*taum(j)-c(j)*taum(j+1)
     *           -2.*((cf1*sigxm(j)+cf2*sigym(j))*fgamma(j)
     *           +(cf3*fgamma(j)+.5*gammah*fgammas(j))*et(j))*dudy(j)
          sigy(j)=-a(j)*sigym(j-1)-(b(j)+4.*dti)*sigym(j)
     *            -c(j)*sigym(j+1)-2.*fgamma(j)*cfy*taum(j)*dudy(j)
          epsm=epsp
 60     continue
        tauu=taum(jmax)
        sigxu=sigxm(jmax)
        sigyu=sigym(jmax)
        call tri(a,b,c,tau,1,jmax,0,0.,0,tauu)
        call tri(axx,bxx,cxx,sigx,1,jmax,0,0.,0,sigxu)
        call tri(ayy,byy,cyy,sigy,1,jmax,0,0.,0,sigyu)
        do 70 j=2,jend
          if(tau(j).lt.0.) tau(j)=.5*(tau(j-1)+taum(j+1))
          upup=sigx(j)+.6666667*et(j)
          if(upup.lt.0.) sigx(j)=-.6666667*et(j)
          if(upup.gt.2.*et(j)) sigx(j)=.6666667*et(j)
          vpvp=sigy(j)+.6666667*et(j)
          if(vpvp.lt.0.) sigy(j)=-.6666667*et(j)
          if(vpvp.gt.2.*et(j)) sigy(j)=.6666667*et(j)
 70     continue
        tau(jmax)=2.*tau(jend)-tau(jend-1)
        sigx(jmax)=2.*sigx(jend)-sigx(jend-1)
        sigy(jmax)=2.*sigy(jend)-sigy(jend-1)
        do 80 j=1,jmax
          sigxm(j)=sigx(j)
          sigym(j)=sigy(j)
          taum(j)=tau(j)
 80     continue
 90   continue
      if(vwplus.eq.0.) then
        do 100 j=1,jmax
          dtau(j)=tau(j)-eps(j)*dudy(j)
 100    continue
      endif
c=============================
c          Momentum
c=============================
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
        fk1=ddy*(vwplus*uplus(j-1)+1.-dtauo)*epsio
        fk2=ddy*(vwplus*(uplus(j-1)+.5*fk1)+1.-dtauh)*epsion
        fk3=ddy*(vwplus*(uplus(j-1)+.5*fk2)+1.-dtauh)*epsion
        fk4=ddy*(vwplus*(uplus(j-1)+fk3)+1.-dtaun)*epsin
        uplus(j)=uplus(j-1)+.16666667*(fk1+2.*(fk2+fk3)+fk4)
        epso=epsn
        epsio=epsin
        if(model.ne.0) dtauo=dtaun
 120  continue
      return
      end
      subroutine edit
c=======================================================================
c              D I S P L A Y   P R O G R A M   O U T P U T
c=======================================================================
      parameter (jmx=401)
      dimension ell(jmx)
      common/bcom/ dlti(jmx),dudy(jmx),dyplus(jmx),eps(jmx),et(jmx),
     *             etm(jmx),ret(jmx),uplus(jmx),uplusm(jmx),wt(jmx),
     *             wtm(jmx),y1(jmx),y2(jmx),y3(jmx),yplus(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx),vextra(jmx)
      common/vcom/ bbeta(jmx),climit(jmx),fc1(jmx),fgamma(jmx),
     *             fgammas(jmx),gamma(jmx),gammas(jmx)
      common/prnt/ diss(jmx),enut(jmx),prod(jmx)
c-----------------------------------------------------------------------
      common/data/     omegw  ,  vwplus ,  ymax   ,  yone
      common/extra/    cfl    ,  ctilda ,  cwall  ,  cc0    ,  cc1    ,
     *                 erro   ,  error  ,  jskipm ,  jskipp ,  nstep  ,
     *                 rate   ,  umax
      common/fctcom/   aa     ,  aam    ,  jend   ,  jendm
      common/iodef/    iunit1 ,  iunit2 ,  model
      common/ndata/    iruff  ,  jmax   ,  jskip  ,  maxn   ,  nedit  ,
     *                 nfreq  ,  nvisc
      common/stress/   alphah ,  betah  ,  c1     ,  gammah ,  onemab
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
c=======================================================================
      cc1=0.
      ctilda=uplus(jmax)-umax
      if(abs(vwplus).gt.1.e-6) cc1=(ctilda-cc0)/vwplus
      write(iunit1,100) nstep
      write(iunit1,110) ctilda,cc1
      con=1./sqrt(alphas)
      do 10 j=1,jmax
        ell(j)=con*sqrt(et(j))/wt(j)
        enut(j)=climit(j)*eps(j)
        if(model.eq.0) then
          prod(j)=enut(j)*dudy(j)**2
        else
          prod(j)=tau(j)*dudy(j)
        endif
        diss(j)=betas*bbeta(j)*et(j)*wt(j)
 10   continue
      if(nvisc.ne.0) diss(1)=2.*diss(2)-diss(3)
      write(iunit1,120)
      write(iunit1,130) (j,yplus(j),uplus(j),et(j),wt(j),enut(j),
     *                           ell(j),prod(j),diss(j), j=1,jmax)
      return
c=======================================================================
 100  format(//,'Nstep =',i5)
 110  format(/5x,'Constant in law of wall = ',1pe12.5,' = CO + ',
     *                                        e12.5,' * Vw/Utau')
 120  format(/5x,'j',4x,'yplus',7x,'uplus',7x,'kplus',6x,'omega+',
     *        7x,'nut',9x,'el+',8x,'prodk',7x,'dissk')
 130  format(i6,1p8e12.5)
c=======================================================================
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
      parameter (jmx=401)
      common/bcom/ dlti(jmx),dudy(jmx),dyplus(jmx),eps(jmx),et(jmx),
     *             etm(jmx),ret(jmx),uplus(jmx),uplusm(jmx),wt(jmx),
     *             wtm(jmx),y1(jmx),y2(jmx),y3(jmx),yplus(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx),vextra(jmx)
c-----------------------------------------------------------------------
      common/data/     omegw  ,  vwplus ,  ymax   ,  yone
      common/extra/    cfl    ,  ctilda ,  cwall  ,  cc0    ,  cc1    ,
     *                 erro   ,  error  ,  jskipm ,  jskipp ,  nstep  ,
     *                 rate   ,  umax
      common/fctcom/   aa     ,  aam    ,  jend   ,  jendm
      common/iodef/    iunit1 ,  iunit2 ,  model
      common/ndata/    iruff  ,  jmax   ,  jskip  ,  maxn   ,  nedit  ,
     *                 nfreq  ,  nvisc
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
      dyplus(1)=yone
      do 10 j=2,jend
        dyplus(j)=rate*dyplus(j-1)
        yplus(j)=yplus(j-1)+dyplus(j-1)
 10   continue
      dyplus(jmax)=dyplus(jend)
      yplus(jmax)=ymax
      do 20 j=2,jend
        denom=dyplus(j)+dyplus(j-1)
        y1(j)=2./(dyplus(j)*denom)
        y2(j)=vwplus/denom
        y3(j)=2./(dyplus(j-1)*denom)
 20   continue
      return
      end
      subroutine namin(icpu)
c=======================================================================
c            COORDINATE READING INPUT DATA FILE 'SUBLAY.DAT'
c=======================================================================
      parameter (niu=3,nn=7,nxd=4)
      character a*12,newfil*7,pform*9
      common/data/  xd(nxd)
      common/iodef/ iu(niu)
      common/ndata/ n(nn)
c=======================================================================
c----------------------------------
c         Read I/O units
c----------------------------------
      call namsys(icpu,iin,iv,msdos,newfil,pform)
      open(iin,file='sublay.dat',status='old')
      iu(1)=2
      iu(2)=7
      io=iu(1)
      open(io,file='sublay.prt',form=pform,status=newfil)
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
      open(iu(2),file='plots.dat',form=pform,status=newfil)
      return
c=======================================================================
 100  format(1x,a12,i4)
 200  format(1x,'Program SUBLAY input data:'/,
     *       1x,'__________________________'/)
 300  format(1x,a12,1pe13.6)
 400  format(1x,'__________________________'/)
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
      subroutine start
c=======================================================================
c            S E T   U P   I N I T I A L   P R O F I L E S
c=======================================================================
      parameter (jmx=401)
      common/bcom/ dlti(jmx),dudy(jmx),dyplus(jmx),eps(jmx),et(jmx),
     *             etm(jmx),ret(jmx),uplus(jmx),uplusm(jmx),wt(jmx),
     *             wtm(jmx),y1(jmx),y2(jmx),y3(jmx),yplus(jmx)
      common/ccom/ dtau(jmx),sigx(jmx),sigxm(jmx),sigy(jmx),sigym(jmx),
     *             tau(jmx),taum(jmx),vextra(jmx)
      common/vcom/ bbeta(jmx),climit(jmx),fc1(jmx),fgamma(jmx),
     *             fgammas(jmx),gamma(jmx),gammas(jmx)
c-----------------------------------------------------------------------
      common/data/     omegw  ,  vwplus ,  ymax   ,  yone
      common/extra/    cfl    ,  ctilda ,  cwall  ,  cc0    ,  cc1    ,
     *                 erro   ,  error  ,  jskipm ,  jskipp ,  nstep  ,
     *                 rate   ,  umax
      common/fctcom/   aa     ,  aam    ,  jend   ,  jendm
      common/iodef/    iunit1 ,  iunit2 ,  model
      common/ndata/    iruff  ,  jmax   ,  jskip  ,  maxn   ,  nedit  ,
     *                 nfreq  ,  nvisc
      common/stress/   alphah ,  betah  ,  c1     ,  gammah ,  onemab
      common/stres2/   cf1    ,  cf2    ,  cf3    ,  cfx    ,  cfy
      common/twoeqn/   alpha  ,  alphas ,  beta   ,  betas  ,  clim   ,
     *                 sigma  ,  sigmad ,  sigmas
      common/viscus/   alf0   ,  alfs0  ,  bbeta0 ,  gamh0  ,  rb     ,
     *                 rk     ,  rw
c=======================================================================
c Set viscous modification coefficients
c=======================================
      alfs0=1.
      alf0=1.
      gamh0=1.
      rb=0.
      rk=0.
      rw=0.
c=======================================
c     k-omega model coefficients
c=======================================
      akappa=.41
      betas=.09
      alphas=sqrt(betas)
      if(model.eq.0) then
        alpha=.52
        beta=.0708
        sigma=.5
        sigmad=.125
        sigmas=.6
        clim=7./8.
        if(nvisc.ne.0) then
          alfs0=beta/3.
          alf0=1./9.
          bbeta0=beta/(3.*betas)
          rb=8.
          rk=6.
          rw=2.61
        endif
c==============================
c  Stress-omega coefficients
c==============================
      elseif(model.eq.1) then
        alpha=.52
        beta=.0708
        sigma=.5
        sigmad=.125
        sigmas=.6
        c1=1.8
        c2=10./19.
        alphah=(8.+c2)/11.
        betah=(8.*c2-2.)/11.
        gammah=(60.*c2-4.)/55.
        if(nvisc.ne.0) then
          alfs0=beta/3.
          alf0=1./9.
          bbeta0=beta/(3.*betas)
          gamh0=21./2000.
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
      endif
c=============================
c Initial profile coefficients
c=============================
      ymatch=10.
      cc0=5.
      cwall=2.*(1.-sigma)/sigma+2.*akappa*cc0-1.
      umax=(1.+.25*(cwall+dlog(ymax))*vwplus/akappa)*dlog(ymax)/akappa
      if(iruff.eq.0) then
        omegi=1./omegw
      else
        omegi=0.
      endif
      con=sqrt(beta/6.)
      ww=(1.+.25*(cwall+2.*dlog(ymatch))*vwplus/akappa)
     *                          /(alphas*akappa*ymatch)
      uu=(1.+.25*(cwall+dlog(ymatch))*vwplus/akappa)*dlog(ymatch)/akappa
      ctilda=ymatch-uu
      uu=ymatch
      ee=(1.+uu*vwplus)/alphas
c===========================================
c Main loops for setting up initial profiles
c===========================================
      do 10 j=1,jmax
        y=yplus(j)
        if(y.le.ymatch) then
          yhat=con*yplus(j)+sqrt(omegi)
          if(iruff.ne.0.and.j.eq.1) yhat=con*yplus(2)
          wt(j)=1./yhat**2
          et(j)=ee*(y/ymatch)**2
          uplus(j)=y
        else
          wt(j)=(1.+.25*(cwall+2.*dlog(y))*vwplus/akappa)
     *                                 /(alphas*akappa*y)
          if(wt(j).gt.omegw) wt(j)=omegw
          uplus(j)=(1.+.25*(cwall+dlog(y))*vwplus/akappa)*dlog(y)/akappa
     *            +ctilda
          et(j)=(1.+vwplus*uplus(j))/alphas
        endif
 10   continue
      do 20 j=1,jmax
        bbeta(j)=1.
        uplusm(j)=uplus(j)
        etm(j)=et(j)
        wtm(j)=wt(j)
        dtau(j)=0.
        if(model.ne.0) then
          sigx(j)=.243*et(j)
          sigy(j)=-.239*et(j)
          tau(j)=.3*et(j)
          sigxm(j)=sigx(j)
          sigym(j)=sigy(j)
          taum(j)=tau(j)
        endif
 20   continue
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
