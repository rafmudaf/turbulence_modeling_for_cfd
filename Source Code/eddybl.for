      program eddybl
c=======================================================================
c                     P R O G R A M   E D D Y B L
c-----------------------------------------------------------------------
c  A two-dimensional/axisymmetric compressible boundary layer program
c                             Created by
c                  D. C. Wilcox, DCW Industries, Inc.
c                         La Canada, California
c-----------------------------------------------------------------------
c    Copyright (c) 2006 DCW Industries, Inc.  All rights reserved.
c-----------------------------------------------------------------------
c
c      ----------------------------------------------------------
c        GENERAL PURPOSE ROUTINES - AFFECTED BY TURBULENCE MODEL
c      ----------------------------------------------------------
c
c      EDDYBL...The main program - it coordinates the computation.
c      CHECK....Check for solution convergence (called by EDDYBL).
c      CLOSUR...Set closure coefficients and other turbulence
c               model parameters (called by EDDYBL).
c      DIMSAV...Write profiles on disk file IUPLOT in terms of 
c               physical (dimensional) quantities.  Useful for
c               plotting software (called by PRINT).
c      MUT......Select and call appropriate eddy-viscosity
c               subroutine (called by EDDYBL).
c      PRINT....Print profiles, call DUMPIT to save data for a 
c               restart run, and call DIMSAV to save data for
c               plotting (called by EDDYBL).
c      SOLVFT...Solve the mean momentum and energy equations
c               (called by EDDYBL).
c      STRTME...Read initial profiles from disk file IUNIT4 and
c               edge and surface conditions from disk file IUNIT2
c               (called by EDDYBL).
c      WALL.....Print integral parameters, edge and surface
c               conditions, convergence data, etc. (called by EDDYBL).
c
c      ----------------------------------------------------------
c      GENERAL PURPOSE ROUTINES - INDEPENDENT OF TURBULENCE MODEL
c      ----------------------------------------------------------
c
c      ADLINE...Adds mesh points as needed to accomodate boundary
c               layer growth (called by WIDTH).
c      DUMPIT...Writes information to disk file IUNIT3 sufficient 
c               for a restart run (called by PRINT).
c      EDGCON...Compute boundary-layer-edge and surface conditions
c               using cubic splines (called by EDDYBL).
c      EDGEX....Locate the boundary-layer edge (called by EDDYBL).
c      INTGRL...Compute integral parameters (called by EDDYBL).
c      INVERT...Invert a 3x3 matrix (called by SOLVSW).
c      NAMIN....Read input data file EDDYBL.DAT (called by EDDYBL).
c      NAMSYS...Set up system-dependent parameters dictated by 
c               parameter ICPU in include file CPUID (called by NAMIN).
c      PROPM....Compute centered and extrapolated mean flow properties
c               (called by SOLVBB, SOLVFT, SOLVKE, SOLVKW, SOLVSA,
c               and SOLVSW).
c      RANKIN...Compute freestream conditions from the Rankine-Hugoniot
c               relations (called by EDDYBL).
c      SOLVV....Solve for the vertical velocity (called by EDDYBL).
c      STABLE...Compute diagonal-dominance weighting factor (called 
c               by SOLVKE and SOLVKW).
c      UPDATE...Compute molecular viscosity, velocity gradient and 
c               temperature gradient (called by EDDYBL).
c      WIDTH....Check for sufficient grid width (called by EDDYBL).
c      XGRID....Compute difference coefficients for streamwise 
c               derivatives and extrapolation (called by EDDYBL).
c      YGRID....Compute normal-grid difference coefficients (called
c               by EDDYBL).
c
c      ----------------------------------------------------------
c                 TURBULENCE-MODEL SPECIFIC ROUTINES
c      ----------------------------------------------------------
c
c      EFMU12...Compute low-Reynolds-number terms for k-epsilon
c               models (called by SOLVKE).
c      SOLVBB...Solve the Baldwin-Barth model equation and
c               compute the eddy viscosity (called by MUT).
c      SOLVBL...Compute the Baldwin-Lomax model algebraic eddy
c               viscosity (called by MUT).
c      SOLVCS...Compute the Cebeci-Smith model algebraic eddy
c               viscosity (called by MUT).
c      SOLVJK...Solve the Johnson-King model equation and
c               compute the eddy viscosity (called by MUT).
c      SOLVKE...Solve the k-epsilon model equations and
c               compute the eddy viscosity (called by MUT).
c      SOLVKW...Solve the k-omega model equations and
c               compute the eddy viscosity (called by MUT).
c      SOLVSA...Solve the Spalart-Allmaras model equation and
c               compute the eddy viscosity (called by MUT).
c      SOLVSW...Solve the Stress-omega model equations (called by MUT).
c      WBCKE....Set surface boundary conditions on k and epsilon
c               for the k-epsilon model (called by SOLVKE).
c      WBCKW....Set surface boundary conditions on k and omega
c               for the k-omega model (called by SOLVKW).
c=======================================================================
      include 'comeb.'
      include 'cpuid.'
c=======================================================================
c          Read input data
c=======================================
      call namin(icpu)
c=======================================
c          Initialization
c=======================================
      if(iutype.eq.0) then
        heateq=778.26
      else
        heateq=1.
      endif
c-------------------
c   Stress limiter
c-------------------
      clim=7./8.
      clim2=clim**2
      do 5 n=1,jk
        climit(n)=1.
 5    continue
c-------------------
c    Laminar run
c-------------------
      if(model.lt.0) then
        lamnar=-1
        ibound=0
c-------------------
c   Turbulent run
c-------------------
      else
        lamnar=0
      endif
      ipoint=0
      isep=0
      iset=0
      mestop=0
      nskip=2
      emax=ziotae
      dsmin=.1*ds
c=======================================
c      Set closure coefficients
c=======================================
      call closur
c=======================================
c  Compute freestream conditions from
c    the Rankine-Hugoniot relations
c=======================================
      call rankin
c=======================================
c   Read starting profiles from disk
c=======================================
      call strtme
      do 10 numbr=1,ipro
        if(si.ge.proval(numbr)) proval(numbr)=1.e10
 10   continue
c=======================================
c     Set up grid normal to wall
c=======================================
      iedgem=iedge-1
      call ygrid
      ete=eteo
      wte=wteo
      eteo=eten
      wteo=wten
      ue=ustrt
      re=restrt
      te=testrt
      xnue=te**vispow*(1.+su/tref)/(te+su/tref)
      rmi=rstrt
      x=xistrt
      oy=2.*x
      dx1=done
      s=si
      sm1=si1
      dx2=dtwo
      dxds=dhelp
      dx1ds=d1help
c========================================
c    Compute initial viscosity and 
c   velocity gradient distributions
c========================================
      tr=t10*tc/te
      do 20 n=1,iedge
        xlm11(n)=(1.+tr)*(1.+to(n))**vispom/(1.+to(n)+tr)
        if(n.ne.1.and.n.ne.iedge) then
          fz(n)=y4(n)*fp(n+1)-y5(n)*fp(n)-y6(n)*fp(n-1)
        endif
 20   continue
      xlm1p=xlm11(1)
      fz(1)=(-ww1*fp(1)+ww2*fp(2)-ww3*fp(3)+ww4*fp(4))/(ww5*dy(1))
      fz(iedge)=0.
c=======================================================================
c   THIS IS THE BEGINNING OF THE MAIN LOOP FOR MARCHING ALONG SURFACE
c=======================================================================
      do 100 m=mstart,iend1
c=======================================
c    Increment s and x, set flags,
c         and save rmi and ue
c=======================================
        xistrt=x
        dhelp=dxds
        d1help=dx1ds
        done=dx1
        si1=sm1
        si=s
        sm2=sm1
        sm1=s
        ds=sm1-sm2
        if(iset.ne.0) then
          ds=dssav
          iset=0
c=============================
c  Change stepsize depending
c   on how well solution is
c         converging
c=============================
        elseif(m.ne.2) then
          if(itstop.lt.11) then
            dsmult=1.1-.02*(itstop-1)
          else
            dsmult=.9
          endif
          ds=dsmult*ds
          if(ds.lt.dsmin) ds=dsmin
        endif
        s=sm1+ds
c=============================
c   Check for profile print
c        and/or sstop
c=============================
        iset=0
        do 30 numbr=1,ipro
          if(s.ge.proval(numbr)) go to 40
 30     continue
        go to 50
 40     s=proval(numbr)
        dssav=ds
        ds=s-sm1
        proval(numbr)=1.e10
        iset=1
 50     if(s.ge.sstop) then
          s=sstop
          ds=s-sm1
          mestop=1
        endif
c=============================
c  Shift edge values for k,
c  omega, rmi, ue.....also,
c  reset iteration counters
c=============================
        eten=eteo
        wten=wteo
        eteo=ete
        wteo=wte
        rstrtm=rmi
        ustrtm=ue
        itstop=1
        itro=0 
c=======================================
c Read B.L. edge and surface conditions
c=======================================
        call edgcon(dpeds,drdz,curv,ipoint)
        rstrt=rmi
        if(kodwal.ne.1) qw=qw*visref*u1*u1/heateq
        rmi1j=rmi**j
        rmi2j=rmi1j*rmi1j
        phi=atan(drdz)
        costh=cos(phi)
        pp=dpeds
        if(m.eq.mstart.and.m.ne.2) then
          ue=ustrt
        else
          ue=sqrt(2.*t10*(1.-(pe/p10)**gm1og))
        endif
        ustrt=ue
        te=t10-.5*ue*ue
        re=pe/(gm1og*te)
        restrt=re
        testrt=te
        xnue=te**vispow*(1.+t10*tc)/(te+t10*tc)
c=======================================
c   Compute incremental changes in x
c=======================================
        dx2ds=dx1ds
        dx1ds=dxds
        dxds=re*ue*xnue*rmi2j
        if(m.ne.2) then
          ckk=(s-sm1)/(sm1-sm2)
          if(abs(ckk-1.).gt.1.e-7) then
            dx2=(s-sm1)*(dx1ds+dxds)/2.
          else
            dx2=(s-sm1)*(3.*dxds+6.*dx1ds-dx2ds)/8.
          endif
        endif
        dtwo=dx2
        x=x+dx2
        oy=2.*x
        oz=sqrt(oy)
c=======================================
c Compute quantities used repeatedly in
c  the main loop which remain constant
c         for all iterations
c=======================================
        duedx=-pp/(re*ue*dxds)
        dtedx=-ue*duedx
c=============================
c Uncouple momentum and energy
c    equations for M -> 0
c=============================
        if(xma.gt..2) then
          xal=ue*ue/te
        else
          xal=0.
        endif
c=============================
c    Useful coefficients
c=============================
        xbe=oy*duedx/ue
        faa=oz/(re*ue*rmi1j)
        tzmu=-heateq*gm1og*qw*eps*pr*oz/(visref*tref*r*te*ue*xnue*rmi1j)
c=============================
c   Set vertical velocity 
c      at the surface
c=============================
        vwall=faa*rvwald/(r1*u1*eps*xnue)
c=============================
c   Coefficients specific to
c      turbulence models
c=============================
        deltd=eps*curv
c-------------------
c   k-omega model
c-------------------
        coef2=oy/(ue*eps)
        coef7=2.*(xbe-1.)
        coef9=1./(re*xnue*rmi2j)
        coef11=oz*deltd/(re*ue*rmi1j)
        coef15=4.5*coef11
        coef69=re*oy/(xnue*eps**2)
c-------------------
c Stress-omega model
c-------------------
        coef3=oz/(eps*xnue*rmi1j)
        coef5=re*rmi1j*oz/eps
        coef13=coef3*coef11
c-------------------
c  k-epsilon model
c-------------------
        coef8=2.*(2.*xbe-1.)
        coef68=2.*oy*xnue/(re*ue**2)
c-------------------
c   Cebeci-Smith
c-------------------
        coef1=(re*ue)**2*rmi1j/(oz*xnue*eps**3)
        coef4=eps*xnue*ue**2*rmi1j/oz
        coef6=xbe*coef4**2/ue
c-------------------
c   Baldwin-Lomax
c-------------------
        coef10=re/(xnue*eps**2)
        coef12=re*ue**2*rmi1j/(eps*oz)
c-------------------
c   Johnson-King
c-------------------
        coef14=oy/(eps*re*ue*xnue*rmi2j)
c-------------------
c Spalart-Allmaras
c-------------------
        coef16=oy**2/(re*ue**2*xnue*rmi2j)
        coef17=coef3*(eps*akappa)**2/coef16
c==============================
c    Turbulence-model edge 
c         conditions
c==============================
        tr=t10*tc/te
        call xgrid(m)
c-------------------
c   Baldwin-Barth
c-------------------
        if(model.eq.3) then
          ete=30./(betas*coef69)
          wte=wteo
c-------------------
c Spalart-Allmaras
c-------------------
        elseif(model.eq.4) then
          ete=eteo*(1.-2.*dx2/oy)
          wte=wteo
c-------------------
c k-omega/Stress-omega
c-------------------
        else
          a=2.*dx2*beta*coef9/oy
          b=z1+2.*dx2*coef7/oy
          c=z3*wten-z2*wteo
          d=b*b-4.*a*c
          wte=wteo
          if(d.gt.0.) wte=.5*(sqrt(d)-b)/a
          ete=(z2*eteo-z3*eten)/(z1+2.*dx2*(2.*xbe+betas*coef9*wte)/oy)
c-------------------
c  k-epsilon model
c-------------------
          epe=betas*ete*wte
        endif
c=======================================
c  Calculate distance normal to surface
c          and eddy viscosity
c=======================================
        if(lamnar.ge.0.and.m.eq.mstart) then
          xn(1)=0.
          tpk=0.
          do 60 n=2,iedge
            tpk=tpk+.5*(2.+tp(n-1)+tp(n))*dy(n-1)
            xn(n)=faa*tpk
 60       continue
          call edgex(ierr)
          xnden=xn(netaeg)
          call mut(1)
        endif
c=======================================================================
c  Beginning of the main program loop.   The  equations  of  motion are
c  integrated from eta=0 to the boundary-layer edge in this loop.  Note
c  that the momentum and energy equations  are solved  first  (iturb=0)
c  and then the model equations are solved  (iturb=1).   This procedure
c  continues iteratively (itro=1,2,3,...) until convergence is attained.
c  The various subroutines called and their functions are as follows...
c-----------------------------------------------------------------------
c  SOLVFT..Solve the Mean-Flow Equations for velocity and temperature
c  SOLVV...Solve for vertical velocity
c  WIDTH...Check for sufficient grid width and growth
c  EDGEX...Determine location of boundary-layer edge
c  UPDATE..Update molecular viscosity, dU/dy, dT/dy, etc.
c  MUT.....Update eddy viscosity
c  CHECK...Check for convergence
c=======================================================================
 70     iturb=0
c-------------------
c Momentum & energy
c-------------------
        call solvft(isep)
        if(isep.ne.0) then
          write(*,230) s
          write(iunit1,230) s
          s=s-ds
          go to 110
        endif
c-------------------
c   Continuity
c-------------------
        call solvv
c-------------------
c Check grid width
c-------------------
        call width(iedgem,ierr)
        if(ierr.ne.0) then
          write(iunit1,200) jk
          go to 110
        endif
c-------------------
c  Find B.L. edge
c-------------------
        call edgex(ierr)
        if(ierr.ne.0) then
          write(iunit1,210)
          go to 110
        endif
c-------------------
c Update Uy, Ty, etc
c-------------------
        call update
        if(lamnar.ge.0) then
          iturb=1
c-------------------
c  Eddy viscosity
c-------------------
          call mut(0)
c-------------------
c Convergence check
c-------------------
          call check
c=======================================
c   Not converged...Store arrays and 
c     return for another iteration
c=======================================
          if(itro.ne.itstop) then
            itro=itro+1
            do 80 n=2,iedge
              fpp(n)=fp(n)
              etpp(n)=etp(n)
c-------------------
c  k-epsilon model
c-------------------
              if(model.eq.2) then
                eppp(n)=epp(n)
                wtp(n)=epp(n)/(betas*etp(n))
c-------------------
c k-omega/Stress-omega
c-------------------
              else
                wtpp(n)=wtp(n)
              endif
 80         continue
c-------------------
c Johnson-King model
c-------------------
            uvmxpp=uvmaxp
            go to 70
          endif
        endif
c=======================================================================
c          THIS IS THE END OF THE MAIN INTEGRATION LOOP
c=======================================================================
        if((jprint+kprint).ne.0) mestop=99
        if(m.lt.(mstart+8).and.mestop.ne.1) mestop=0
c=======================================
c     Compute integral parameters
c=======================================
        call intgrl
c=======================================
c        Normal wall printout
c=======================================
        if(ishort.ne.0) then
          call wall(m)
c=======================================
c        Short wall printout
c=======================================
        else
          write(iunit1,220) m,s,res,rethet,cfe,thadis
        endif
c=======================================
c          Print to console
c=======================================
        write(*,220) m,s,res,rethet,cfe,thadis
c=======================================
c           Print profiles
c=======================================
        if(iset.ne.0) call print(m)
c=======================================
c  Check for premature run termination
c=======================================
        if(mestop.ne.0.or.m.eq.iend1) go to 110
c=======================================
c    Update variables for marching
c=======================================
        do 90 n=1,iedge
c-------------------
c     Mean-flow
c-------------------
          fn(n)=fo(n)
          fo(n)=fp(n)
          tn(n)=to(n)
          to(n)=tp(n)
          vn(n)=vo(n)
          vo(n)=vp(n)
          raton(n)=ratoo(n)
          ratoo(n)=ratop(n)
c-------------------
c k-omega/Stress-omega
c-------------------
          etn(n)=eto(n)
          eto(n)=etp(n)
          wtn(n)=wto(n)
          wto(n)=wtp(n)
c-------------------
c  k-epsilon model
c-------------------
          epn(n)=epo(n)
          epo(n)=epp(n)
c-------------------
c Stress-omega model
c-------------------
          uun(n)=uuo(n)
          uuo(n)=uup(n)
          uvn(n)=uvo(n)
          uvo(n)=uvp(n)
          vvn(n)=vvo(n)
          vvo(n)=vvp(n)
          wwn(n)=wwo(n)
          wwo(n)=wwp(n)
 90     continue
        dx1=dx2
c-------------------
c Johnson-King model
c-------------------
        if(model.eq.7) then
          uvmaxn=uvmaxo
          uvmaxo=uvmaxp
        endif
 100  continue
c=======================================================================
c            END OF DO LOOP FOR MARCHING ALONG SURFACE
c=======================================================================
c    Print profiles for last station
c=======================================
 110  if(iset.ne.1) call print(m)
c=======================================
c   Terminate wall-property plot file
c=======================================
      dum=-999.
      write(iunit5,260) dum,dum
      close(iunit5)
      if(mestop.eq.99) then
        if(model.eq.2) then
          write(*,240) 
          write(iunit1,240)
        else
          write(*,250) 
          write(iunit1,250)
        endif
      endif
c=======================================
c     Close eddybl.prt, output.dat,
c      profil.dat and profil2c.dat
c=======================================
      close(iunit1)
      close(iunit3)
      close(iunit6)
      close(iuplot)
c=======================================
c    Pause to keep window visible after
c    the run completes...
c=======================================
      pause
c=======================================================================
 200  format(/1x,'*** Grid width too small ***'/
     *        1x,'Maximum limit of',i4,' for IEDGE has been exceeded')
 210  format(1x,'*** fp(i) never = or > .9999 ***')
 220  format(1x,'m =',i4,': s =',1pe9.2,1x,'Res =',e9.2,1x,'Re0 =',
     *                            e9.2,1x,'Cf =',e9.2,1x'H =',e9.2)
 230  format(/1x,'*** The boundary layer has separated at s = ',
     *                                            1pe9.2,' ***')
 240  format(/1x,'*** Negative k and/or epsilon...run terminated ***')
 250  format(/1x,'*** Negative k and/or omega...run terminated ***')
 260  format(1p2e12.4)
c=======================================================================
      end
      subroutine adline
c=======================================================================
c     A D D   M E S H   P O I N T S   T O   T H E   G R I D   A S
c     N E E D E D   T O   A C C O M O D A T E   G R O W T H   O F
c               T H E   B O U N D A R Y   L A Y E R
c-----------------------------------------------------------------------
c  NOTE: fp(n), tp(n), uup(n), uvp(n), vvp(n), etc. are set to 
c        constant values in Subroutine STRTME for iedge < n < jk.
c=======================================================================
      include 'comeb.'
c=======================================================================
      iedgem=iedge-1
      dym=xk*dy(iedgem)
      dy(iedge)=dym
c-------------------
c     Mean-flow
c-------------------
      raton(iedge)=raton(iedgem)
      ratoo(iedge)=ratoo(iedgem)
      ratop(iedge)=ratop(iedgem)
      vn(iedge)=vn(iedgem)-dym
      vo(iedge)=vo(iedgem)-dym
      vp(iedge)=vp(iedgem)-dym
c-------------------
c k-omega/Stress-omega
c-------------------
      etn(iedge)=etn(iedgem)
      eto(iedge)=eto(iedgem)
      etp(iedge)=etp(iedgem)
      wtn(iedge)=wtn(iedgem)
      wto(iedge)=wto(iedgem)
      wtp(iedge)=wtp(iedgem)
c-------------------
c  k-epsilon model
c-------------------
      epn(iedge)=epn(iedgem)
      epo(iedge)=epo(iedgem)
      epp(iedge)=epp(iedgem)
      return
      end
      subroutine check
c=======================================================================
c             C H E C K   F O R   C O N V E R G E N C E
c-----------------------------------------------------------------------
c       NOTE: The tolerance 'errmax' is set in subroutine closur.
c=======================================================================
      include 'comeb.'
c=======================================================================
c      First, determine maximum values
c============================================
      emax=0.
      epsmax=0.
      nepmax=1
      do 10 n=2,iedge
        rhoe=etp(n)/(1.+tp(n))
        if(rhoe.gt.emax) emax=rhoe
        if(ep(n).gt.epsmax) then
          epsmax=ep(n)
          nepmax=n
        endif
 10   continue
      lamnar=0 
      if(epsmax.gt..001) lamnar=1
      if(epsmax.ge.10.) ibound=1
      epsmax=epsmax*xlm11(nepmax)*(1.+tp(nepmax))**2
      nerror=0
c============================================
c  One-, two-equation & Stress-omega models
c============================================
      if(model.ge.0.and.model.le.4) then
        do 20 n=2,iedge
          errorf=abs((fp(n)-fpp(n))/fpp(n))
          errork=abs((etp(n)-etpp(n))/etpp(n))
          if(errorf.gt.errmax) nerror=nerror+1
          if(errork.gt.errmax) nerror=nerror+1
 20     continue
c============================================
c       Laminar and all other models
c============================================
      else
        lamnar=1
        do 30 n=2,iedge
          errorf=abs((fp(n)-fpp(n))/fpp(n))
          if(errorf.gt.errmax) nerror=nerror+1
 30     continue
        if(model.eq.7) then
          errors=abs((uvmaxp-uvmxpp)/uvmxpp)
          if(errors.gt.errmax) nerror=nerror+3
        endif
      endif
c============================================
c        Increment iteration counter
c============================================
      if(nerror.gt.2) itstop=itstop+1
      if(itstop.gt.20) itstop=20
      return
      end
      subroutine closur
c=======================================================================
c           S E T   C L O S U R E   C O E F F I C I E N T S
c=======================================================================
      include 'comeb.'
c=======================================================================
c     Miscellaneous coefficients
c=======================================
      akappa=.41
      beta50=sqrt(betas)
      beta25=sqrt(beta50)
      errmax=.0005
      rghcon=sqrt(beta/6.)
      xmtsq=xmt0**2
c-------------------
c Default values...
c    Viscous mods
c      excluded
c-------------------
      alf0=1.
      alfs0=1.
      bbeta0=1.
      gamh0=1.
      rb=0.
      rk=0.
      rw=0.
c------------------------
c     k-omega model
c------------------------
      if(model.eq.0) then
        sigmad=psieps
        if(nvisc.ne.0) then
          alfs0=beta/3.
          alf0=1./9.
          bbeta0=beta/(3.*betas)
          rb=8.
          rk=6.
          rw=2.61
        endif
c------------------------
c   Stress-omega model
c------------------------
      elseif(model.eq.1) then
        sigmad=psieps
        cc2=10./19.
        alphah=(8.+cc2)/11.
        betah=(8.*cc2-2.)/11.
        gammah=(60.*cc2-4.)/55.
        csw1=1.-alphah
        csw2=csw1-betah
        csw3=.6666667*csw2
        csw4=.66666667*(2.*csw1+betah)
        csw5=-.66666667*(csw1+2.*betah)
        csw6=1.-alphah+betah
        if(nvisc.ne.0) then
          alfs0=beta/3.
          alf0=1./9.
          bbeta0=beta/(3.*betas)
          gamh0=21./2000.
          rb=8.
          rk=6.
          rw=22./9.
        endif
c------------------------
c    k-epsilon model
c------------------------
      elseif(model.eq.2) then
        beta=betas*(ceps2-1.)
c------------------------
c  Baldwin-Barth model
c------------------------
      elseif(model.eq.3) then
        a0plus=psieps
        a2plus=sigmak
        sigmae=akappa**2/(sqrt(cmu)*(ceps2-ceps1))
        cbb1=ceps1/ceps2
        cbb2=1.-cbb1
        cbb3=1./(a0plus*a2plus)
        cbb4=.5*(a0plus+a2plus)*cbb3
        cbb5=3.*cbb2*sqrt(cbb3)/akappa
        cbb6=5.*cbb4/6.
        errmax=.0001
c------------------------
c Spalart-Allmaras model
c------------------------
      elseif(model.eq.4) then
        cb1=ceps1
        cb2=ceps2
        cv1=cmu
        cw1=cb1/akappa**2+(1.+cb2)/sigmae
        cw2=sigmak
        cw3=psieps
        csa1=cv1**3
        csa2=cw3**6
        errmax=.0001
        relax1=.85
        relax2=1.-relax1
c------------------------
c   Cebeci-Smith model
c------------------------
      elseif(model.eq.5) then
        akappa=.4
        alphac=sigmae
        a0plus=psieps
        errmax=.0001
c------------------------
c  Baldwin-Lomax model
c------------------------
      elseif(model.eq.6) then
        akappa=.4
        ccp=cmu
        ckleb=ceps1
        cwk=ceps2
        alphac=sigmae
        a0plus=psieps
        errmax=.0001
c------------------------
c   Johnson-King model
c------------------------
      elseif(model.eq.7) then
        akappa=.4
        ajk1=ceps1
        cdif=ceps2
        cjk1=cmu
        cjk2=sigmak
        alphac=sigmae
        a0plus=psieps
        errmax=.0001
        relax1=.5
        relax2=1.-relax1
      endif
      return
      end
      subroutine dimsav(m)
c=======================================================================
c  W R I T E   P R O F I L E S   O N   D I S K   I N   T E R M S   O F
c   P H Y S I C A L   ( D I M E N S I O N A L )   Q U A N T I T I E S
c=======================================================================
      include 'comeb.'
c=======================================================================
      ydd=eps
      udd=u1*ue
      tdd=tref*te
      rdd=r1*re
      edd=udd**2
      wdd=u1*ue**2/oy
      write(iunit6,100) iedge
      write(iuplot,100) iedge
      do 10 n=1,iedge
        tp1=1.+tp(n)
        ydim=ydd*xn(n)
        udim=udd*fp(n)
        tdim=tdd*tp1
        rdim=rdd/tp1
c=======================================
c            Laminar flow
c=======================================
        if(model.eq.-1) then
          edim=0.
          wdim=0.
c=======================================
c  Two-equation and Stress-omega models
c=======================================
        elseif(model.ge.0.and.model.le.2) then
          edim=edd*etp(n)
          wdim=wdd*wtp(n)
c=======================================
c        One-equation models
c=======================================
        elseif(model.eq.3.or.model.eq.4) then
          edim=coef69*etp(n)/(xlm11(n)*(1.+tp(n))**2)
          wdim=0.
c=======================================
c  Algebraic and 1/2-equation models
c=======================================
        else
          edim=0.
          wdim=0.
        endif
c=======================================
c          Stress-omega model
c=======================================
        if(model.eq.1) then
          uvdim=edd*uvp(n)
          uudim=-edd*(uup(n)+.66666667*etp(n))
          vvdim=-edd*(vvp(n)+.66666667*etp(n))
          wwdim=-edd*(wwp(n)+.66666667*etp(n))
c=======================================
c          All other models
c=======================================
        else
          eldim=0.
          uvdim=edd*xlm11(n)*ratop(n)**j*tp1*ep(n)*fz(n)/coef3
          uudim=0.
          vvdim=0.
          wwdim=0.
        endif
        write(iunit6,110) ydim,udim,tdim,rdim,edim,wdim,wwdim,
     *                    uudim,uvdim,vvdim
        write(iuplot,120) yplus(n),uplus(n)
 10   continue
      endfile iunit6
      endfile iuplot
      rewind iunit6
      rewind iuplot
      write(iunit1,130) m
      return
c=======================================================================
 100  format(i6)
 110  format(1p10e12.4)
 120  format(1p2e12.4)
 130  format(//2x,'Disk files PROFIL.DAT, WALL.DAT and PROFIL2C.DAT '
     *            'written at MSTART =',i4)
c=======================================================================
      end
      subroutine dumpit(m)
c=======================================================================
c          W R I T E   I N F O R M A T I O N   O N   D I S K
c          S U F F I C I E N T   F O R   R E S T A R T I N G
c=======================================================================
      include 'comeb.'
c=======================================================================
      write(iunit3) m,ibound,iedge
      do 10 n=1,iedge
        write(iunit3) dy(n),etn(n),eto(n),fn(n),fo(n),tn(n),to(n),
     *                vn(n),vo(n),wtn(n),wto(n),raton(n),ratoo(n),
     *                wwn(n),wwo(n),uun(n),uuo(n),uvn(n),uvo(n),
     *                vvn(n),vvo(n),deluv(n)
 10   continue
      write(iunit3) dhelp,done,dtwo,d1help,eten,eteo,itstop,
     *              restrt,rstrt,rstrtm,si,si1,testrt,ustrt,
     *              ustrtm,wten,wteo,xistrt
      endfile iunit3
      rewind iunit3
      write(iunit1,100) m
      return
c=======================================================================
 100  format(//2x,'Disk file "OUTPUT.DAT" written at MSTART =',i4)
c=======================================================================
      end
      subroutine edgcon(dpeds,drdz,curv,i)
c=======================================================================
c       C O M P U T E   B O U N D A R Y - L A Y E R - E D G E
c            A N D   S U R F A C E   C O N D I T I O N S
c=======================================================================
      include 'comeb.'
c=======================================================================
c    Reset pointer if out of range
c=======================================
      if(i.le.0) i=1
      if(i.ge.number) i=number-1
      ip=i+1
c=======================================
c   Point below minimum s for which
c    edge conditions are specified
c=======================================
      if(s.le.as(1)) then
        i=1
        pe=ape(1)
        rmi=armi(1)
        tw=atw(1)
        z=az(1)
        dpeds=0.
        rvwald=arhov(1)
        drdz=0.
        qw=aqw(1)
        curv=azcrv(1)
c=======================================
c   Point above maximum s for which
c    edge conditions are specified
c=======================================
      elseif(s.ge.as(number)) then
        i=number-1
        pe=ape(number)
        rmi=armi(number)
        tw=atw(number)
        z=az(number)
        dpeds=0.
        rvwald=arhov(number)
        drdz=0.
        qw=aqw(number)
        curv=azcrv(number)
c=======================================
c     Point  within range for which
c     edge conditions are specified
c=======================================
      else
        if((as(i)-s)*(as(ip)-s).gt.0.) then
          in=sign(1.,(as(ip)-as(i))*(s-as(i)))
 10       if((i+in).gt.0.and.(i+in).lt.number) then
            i=i+in
            ip=i+1
            if((as(i)-s)*(as(ip)-s).gt.0.) go to 10
          endif
        endif
        hh=as(ip)-as(i)
        dels=(s-as(i))/hh
        delsp=(as(ip)-s)/hh
        phi1=dels**2*(3.-2.*dels)
        phi0=1.-phi1
        psi1=-hh*delsp*dels**2
        psi0=hh*dels*delsp**2
        phi1p=6.*dels*(1.-dels)/hh
        phi0p=-phi1p
        psi1p=-dels*(2.*delsp-dels)
        psi0p=delsp*(delsp-2.*dels)
        pe=ape(i)*phi0+ape(ip)*phi1+apep(i)*psi0+apep(ip)*psi1
        rmi=armi(i)*phi0+armi(ip)*phi1+armip(i)*psi0+armip(ip)*psi1
        tw=atw(i)*phi0+atw(ip)*phi1+atwp(i)*psi0+atwp(ip)*psi1
        z=az(i)*phi0+az(ip)*phi1+azp(i)*psi0+azp(ip)*psi1
        rvwald=arhov(i)*phi0+arhov(ip)*phi1+arhovp(i)*psi0
     *                                    +arhovp(ip)*psi1
        qw=aqw(i)*phi0+aqw(ip)*phi1+aqwp(i)*psi0+aqwp(ip)*psi1
        curv=azcrv(i)*phi0+azcrv(ip)*phi1+azcrvp(i)*psi0+azcrvp(ip)*psi1
        dpeds=ape(i)*phi0p+ape(ip)*phi1p+apep(i)*psi0p+apep(ip)*psi1p
        rprime=armi(i)*phi0p+armi(ip)*phi1p+armip(i)*psi0p
     *                                    +armip(ip)*psi1p
        zprime=az(i)*phi0p+az(ip)*phi1p+azp(i)*psi0p+azp(ip)*psi1p
        if(abs(zprime).gt.1.e-10) then
          drdz=rprime/zprime
        else
          drdz=0.
        endif
      endif
      return
      end
      subroutine edgex(ierr)
c=======================================================================
c        D E T E R M I N E   T H E   L O C A T I O N   O F
c          T H E   B O U N D A R Y   L A Y E R   E D G E
c=======================================================================
      include 'comeb.'
c=======================================================================
      ierr=0
      ncount=0
      do 10 i=1,iedge
        ncount=ncount+1
        if(fp(i).ge..9999) then
          netaeg=ncount
          go to 20
        endif
 10   continue
      ierr=1
 20   return
      end
      subroutine efmu12(n)
c=======================================================================
c     K - E P S I L O N   M O D E L   V I S C O U S   T E R M S
c-----------------------------------------------------------------------
c    NOTE:   Some of the very near wall asymptotic forms below use the
c    exact theoretical relationship between Ry and Ret, Ry=sqrt(2*Ret)
c=======================================================================
      include 'comeb.'
c=======================================================================
c         Jones-Launder Model
c=======================================
      if(nvisc.eq.0) then
        efmu(n)=exp(-2.5/(1.+ret(n)/50.))
        ef1(n)=1.
        ef2(n)=1.-.3*exp(-ret(n)**2)
        dkhdy=y2(n)*(sqrt(etp(n+1))-sqrt(etp(n-1)))
        eps0(n)=2./coef9*xlm11(n)*(dkhdy*ratop(n)**j)**2
        toth1=ratop(n+1)**j/(1.+tp(n+1))
        toth2=ratop(n)**j/(1.+tp(n))
        toth3=ratop(n-1)**j/(1.+tp(n-1))
        cxp=.5*(toth1+toth2)
        cxm=.5*(toth2+toth3)
        fyy=y1(n)*cxp*(fp(n+1)-fp(n))-y3(n)*cxm*(fp(n)-fp(n-1))
        ehat(n)=2./coef9*cmu*efmu(n)*ret(n)*(xlm11(n)*toth2*fyy)**2
c=======================================
c         Launder-Sharma Model
c=======================================
      elseif(nvisc.eq.1) then
        efmu(n)=exp(-3.4/(1.+ret(n)/50.)**2)
        ef1(n)=1.
        ef2(n)=1.-.3*exp(-ret(n)**2)
        dkhdy=y2(n)*(sqrt(etp(n+1))-sqrt(etp(n-1)))
        eps0(n)=2./coef9*xlm11(n)*(dkhdy*ratop(n)**j)**2
        toth1=ratop(n+1)**j/(1.+tp(n+1))
        toth2=ratop(n)**j/(1.+tp(n))
        toth3=ratop(n-1)**j/(1.+tp(n-1))
        cxp=.5*(toth1+toth2)
        cxm=.5*(toth2+toth3)
        fyy=y1(n)*cxp*(fp(n+1)-fp(n))-y3(n)*cxm*(fp(n)-fp(n-1))
        ehat(n)=2./coef9*cmu*efmu(n)*ret(n)*(xlm11(n)*toth2*fyy)**2
c=======================================
c         Lam-Bremhorst model
c=======================================
      elseif(nvisc.eq.2) then
        if(ry(n).gt.1.) then
          efmu(n)=(1.-exp(-.0165*ry(n)))**2*(1.+20.5/ret(n))
        else
          efmu(n)=(.0165*ry(n))**2*(1.-.0165*ry(n))*(1.+20.5/ret(n))
        endif
        if(ret(n).lt.1..and.efmu(n).gt.1.) efmu(n)=.01116225
        if(efmu(n).lt..01116225) efmu(n)=.01116225
        ef1(n)=1.+(.05/efmu(n))**3
        ef2(n)=1.-exp(-ret(n)**2)
        eps0(n)=0.
        ehat(n)=0.
c=======================================
c             Chien Model
c=======================================
      elseif(nvisc.eq.3) then
        if(yplus(n).gt..01) then
          efmu(n)=1.-exp(-.0115*yplus(n))
        else
          efmu(n)=.0115*yplus(n)
        endif
        ef1(n)=1.
        ef2(n)=1.-2./9.*exp(-(ret(n)/6.)**2)
        eps0(n)=coef68*xlm11(n)*etp(n)*((1.+tp(n))/xn(n))**2
        ehat(n)=-coef68*xlm11(n)*epp(n)*((1.+tp(n))/xn(n))**2
     *          *coef9*exp(-.5*yplus(n))
c=======================================
c         Yang and Shih Model
c=======================================
      elseif(nvisc.eq.4) then
        if(ry(n).gt.1.) then
          reth=sqrt(ret(n))
          fact=1.-exp(-1.5e-4*ry(n)-5.e-7*ry(n)**3-1.e-10*ry(n)**5)
        else
          reth=.5*ry(n)
          fact=1.5e-4*ry(n)*(1.+.75e-4*ry(n))
        endif
        ef1(n)=reth/(1.+reth)
        ef2(n)=ef1(n)
        efmu(n)=sqrt(fact)/ef1(n)
        eps0(n)=0.
        toth1=ratop(n+1)**j/(1.+tp(n+1))
        toth2=ratop(n)**j/(1.+tp(n))
        toth3=ratop(n-1)**j/(1.+tp(n-1))
        cxp=.5*(toth1+toth2)
        cxm=.5*(toth2+toth3)
        fyy=y1(n)*cxp*(fp(n+1)-fp(n))-y3(n)*cxm*(fp(n)-fp(n-1))
        ehat(n)=1./coef9*cmu*efmu(n)*ret(n)*(xlm11(n)*toth2*fyy)**2
c=======================================
c  Fan, Lakshminarayana, Barnett Model
c=======================================
      elseif(nvisc.eq.5) then
        ryh=sqrt(ry(n))
        if(ryh.gt..1) then
          efw=1.-exp(-ryh/2.3+(ryh/2.3-ry(n)/8.89)
     *                       *(1.-exp(-ry(n)/20.))**3)
          fbar=efw/sqrt(ret(n))
        else
          efw=ryh/2.3
          fbar=.614875/ryh
        endif
        efmu(n)=.4*fbar+(1.-.4*fbar)*(1.-exp(-ry(n)/42.63))**3
        ef1(n)=1.
        ef2(n)=(1.-2./9.*exp(-(ret(n)/6.)**2))*efw**2
        eps0(n)=0.
        ehat(n)=0.
      endif
      return
      end
      subroutine intgrl
c=======================================================================
c        C O M P U T E   I N T E G R A L   P A R A M E T E R S
c=======================================================================
      include 'comeb.'
c=======================================================================
c  Momentum and displacement thickness
c=======================================
      tpp=1.+tp(1)
      co=0.
      cvo=1.
      disp=0.
      theta=0.
      deltav=0.
      do 10 n=2,iedge
        c=fp(n)*(1.-fp(n))
        cv=(1.+tp(n))*(1.-fp(n))/ratop(n)**j
        disp=disp+.5*(tp(n-1)+tp(n)+2.-fp(n-1)-fp(n))*dy(n-1)
        theta=theta+.5*(co+c)*dy(n-1)
        deltav=deltav+.5*(cvo+cv)*dy(n-1)
        co=c
        cvo=cv
 10   continue
c=======================================
c   Put things in desired dimensional
c        or nondimensional form
c=======================================
      deltav=deltav*eps*faa
      theta=theta*eps*faa
      disp=disp*eps*faa
      disp=abs(disp)
      thadis=disp/theta
      redelt=(re*ue*disp/xnue)*reyref
      redelv=(re*ue*deltav/xnue)*reyref
      rethet=redelt/thadis
      res=(re*ue*s/xnue)*reyref
      xmae=ue/sqrt(gm1*te)
      taud=visref*u1*xlm1p*re*xnue*ue*ue*rmi1j*fz(1)/(oz*eps)
      cfe=taud/(.5*r1*u1*u1*re*ue*ue)
      cfw=cfe*tw/te
      if(kodwal.eq.1) then
        qx=-xlm1p*re*ue*te*xnue*rmi1j*tz(1)/(pr*oz*eps)
        qsd=qx*r*visref*tref/(heateq*gm1og)
      else
        qsd=qw
        recov=tp(1)/(tt1/(te*tref)-1.)
      endif
      twd=tref*tw
      if(lamnar.eq.0) then
        rf=sqrt(pr)
      else
        rf=.89
      endif
      tawd=rf*(tt1-te*tref)+te*tref
      if(abs(twd-tawd).gt.1.e-10) then
        hd=qsd/(twd-tawd)
      else
        hd=0.
      endif
      che=heateq*hd*gm1og/(r*r1*u1*re*ue)
      chw=che*tw/te
      hsd=hd*s
      aked=r*visref*xnue/(heateq*pr*gm1og)
      akwd=aked*xlm1p*tpp
      anue=hsd/aked
      anuw=hsd/akwd
c=======================================
c           Save plotting data
c=======================================
      pestar=pe*r1*u1*u1
      testar=te*tref
      twd=testar*(1.+tp(1))
      write(iunit5,100) s,cfe
      return
c=======================================================================
 100  format(1p2e12.4)
c=======================================================================
      end
      subroutine invert(a,ierr)
c=======================================================================
c              I N V E R T   A   3 X 3   M A T R I X
c=======================================================================
      dimension a(3,3),b(3,3)
c=======================================================================
c      Compute the determinant
c=======================================
      b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
      b(2,1)=a(2,3)*a(3,1)-a(2,1)*a(3,3)
      b(3,1)=a(2,1)*a(3,2)-a(2,2)*a(3,1)
      det=b(1,1)*a(1,1)+b(2,1)*a(1,2)+b(3,1)*a(1,3)
c=======================================
c     Set ierr=1 if it's singular
c=======================================
      if(det.eq.0.) then
        ierr=1
c=======================================
c        Compute the inverse
c=======================================
      else
        ierr=0
        b(1,2)=a(1,3)*a(3,2)-a(1,2)*a(3,3)
        b(1,3)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
        b(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1)
        b(2,3)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
        b(3,2)=a(1,2)*a(3,1)-a(1,1)*a(3,2)
        b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
        deti=1./det
        do 10 i=1,3
          do 10 k=1,3
            a(i,k)=deti*b(i,k)
 10     continue
      endif
      return
      end
      subroutine mut(justmu)
c=======================================================================
c     S E L E C T   E D D Y   V I S C O S I T Y   R O U T I N E
c-----------------------------------------------------------------------
c  model = 0: SOLVKW...k-omega and Stress-omega models
c          1: SOLVSW...Stress-omega model
c          2: SOLVKE...k-epsilon models
c          3: SOLVBB...Baldwin-Barth model
c          4: SOLVSA...Spalart-Allmaras model
c          5: SOLVCS...Cebeci-Smith model
c          6: SOLVBL...Baldwin-Lomax model
c          7: SOLVJK...Johnson-King model
c=======================================================================
      include 'comeb.'
c=======================================================================
      if(model.eq.0.or.model.eq.1) then
        call solvkw(justmu)
        if(model.eq.1.and.justmu.eq.0) call solvsw
      elseif(model.eq.2) then
        call solvke(justmu)
      elseif(model.eq.3) then
        call solvbb(justmu)
      elseif(model.eq.4) then
        call solvsa(justmu)
      elseif(model.eq.5) then
        call solvcs
      elseif(model.eq.6) then
        call solvbl
      elseif(model.eq.7) then
        call solvjk
      endif
      return
      end
      subroutine namin(icpu)
c=======================================================================
c           COORDINATE READING INPUT DATA FILE 'EDDYBL.DAT'
c=======================================================================
      parameter (ipro=10,niu=8,nke=6,nnid=2,nn=11,nxc=6,nxcm=6,nxd=20)
      character a*12,newfil*7,pform*9
      common/close/ xc(nxc)
      common/closem/xcm(nxcm)
      common/data/  xd(nxd)
      common/iodef/ iu(niu)
      common/keps/  xke(nke)
      common/moddef/nid(nnid)
      common/ndata/ n(nn)
      common/pprint/proval(ipro)
c=======================================================================
c----------------------------------
c         Read units type
c----------------------------------
      call namsys(icpu,iin,iv,msdos,newfil,pform)
      open(iin,file='eddybl.dat',status='old')
      iu(1)=2
      do 10 i=2,niu-1
        iu(i)=i+5
 10   continue
      io=iu(1)
      open(io,file='eddybl.prt',form=pform,status=newfil)
      write(io,200)
      read(iin,100) a,iu(8)
      write(io,100) a,iu(8) 
c----------------------------------
c  Read model identification flags
c----------------------------------
      do 20 i=1,nnid
        read(iin,100) a,nid(i)
        write(io,100) a,nid(i) 
 20   continue
c----------------------------------
c Read k-omega closure coefficients
c----------------------------------
      do 30 i=1,nxc
        read(iin,300) a,xc(i)
        write(io,300) a,xc(i) 
 30   continue
c----------------------------------
c  Read k-epsilon closure coeffs
c----------------------------------
      do 40 i=1,nke
        read(iin,300) a,xke(i)
        write(io,300) a,xke(i) 
 40   continue
c----------------------------------
c Read Stress-omega closure coeffs
c----------------------------------
      do 50 i=1,nxcm
        read(iin,300) a,xcm(i)
        write(io,300) a,xcm(i) 
 50   continue
c----------------------------------
c    Read general input (real)
c----------------------------------
      do 60 i=1,nxd
        read(iin,300) a,xd(i)
        write(io,300) a,xd(i) 
 60   continue
c----------------------------------
c  Read general input (integer)
c----------------------------------
      do 70 i=1,nn
        read(iin,100) a,n(i)
        write(io,100) a,n(i) 
 70   continue
c----------------------------------
c   Read profile-print stations
c----------------------------------
      do 80 i=1,ipro
        read(iin,300) a,proval(i)
        write(io,300) a,proval(i) 
 80   continue
c----------------------------------
c Skip over, but print, quantities 
c    used for initial profiles
c----------------------------------
      do 90 i=1,7
        read(iin,300) a,dummy
        write(io,300) a,dummy
 90   continue
      write(io,400)
c----------------------------------
c       Open disk files
c----------------------------------
      open(iu(2),file='table.dat',form='unformatted',status='old')
      open(iu(3),file='output.dat',form='unformatted',status=newfil)
      open(iu(4),file='input.dat',form='unformatted',status='old')
      open(iu(5),file='wall.dat',status=newfil)
      open(iu(6),file='profil2c.dat',status=newfil)
      open(iu(7),file='profil.dat',status=newfil)
      return
c=======================================================================
 100  format(1x,a12,i4)
 200  format(1x,'Program EDDYBL input data:'/,
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
        pause
        stop
      endif
      return
      end
      subroutine print(m)
c=======================================================================
c    P R I N T   P R O F I L E S   A N D   S A V E   D A T A   F O R
c      A   R E S T A R T   R U N   A N D   F O R   P L O T T I N G
c=======================================================================
      include 'comeb.'
c=======================================================================
      wscale=xnue*eps**2/(re*oy)
      elscal=oy/(ue*eps*xnden*beta25)
      xmf4=eps*xnue*rmi**j*xlm11(1)*fz(1)/oz
      uplusb=ue*sqrt(xmf4*(1.+tp(1)))
      tscale=re*uplusb/(eps*xnue*xlm11(1)*(1.+tp(1))**2)
      uvscal=1./(xlm11(1)*fz(1))
      if(model.eq.1) uvscal=coef3*uvscal
      asqrd=.5*(g-1.)*prt*xmae**2*te/tw
      aaa=sqrt(asqrd)
      bnew=te/tw-1.+asqrd
      denom=sqrt(bnew**2+4.*asqrd)
      arcit=asin(bnew/denom)
      uplus(1)=0.
      do 10 n=1,iedge
        if(n.ne.1) uplus(n)=ue*(asin((2.*asqrd*fp(n)-bnew)/denom)+arcit)
     *                     /(aaa*uplusb)
        yplus(n)=tscale*xn(n)
        yodelt=xn(n)/xnden
        temp=1.+tp(n)
c=======================================
c            Laminar flow
c=======================================
        if(model.eq.-1) then
          if(n.eq.1) write(iunit1,100) s
          write(iunit1,110) n,yodelt,fp(n),yplus(n),uplus(n),temp
c=======================================
c  Two-equation and Stress-omega models
c=======================================
        elseif(model.ge.0.and.model.le.2) then
          omeg=wscale*wtp(n)
          esize=elscal*sqrt(etp(n))/wtp(n)
          if(model.eq.1) then
            uvbar=uvscal*uvp(n)/temp
          else
            uvbar=uvscal*ep(n)*ratop(n)**j*xlm11(n)*fz(n)
          endif
          if(n.eq.1) write(iunit1,120) s
          write(iunit1,130) n,yodelt,fp(n),yplus(n),uplus(n),etp(n),
     *                      omeg,ep(n),esize,uvbar,temp
c=======================================
c        One-equation models
c=======================================
        elseif(model.eq.3.or.model.eq.4) then
          chi=coef69*etp(n)/(xlm11(n)*temp**2)
          uvbar=uvscal*ep(n)*ratop(n)**j*xlm11(n)*fz(n)
          if(n.eq.1) write(iunit1,140) s
          write(iunit1,150) n,yodelt,fp(n),yplus(n),uplus(n),
     *                      chi,ep(n),uvbar,temp
c=======================================
c  Algebraic and 1/2-equation models
c=======================================
        else
          uvbar=uvscal*ep(n)*ratop(n)**j*xlm11(n)*fz(n)
          if(n.eq.1) write(iunit1,160) s
          write(iunit1,170) n,yodelt,fp(n),yplus(n),uplus(n),
     *                      ep(n),uvbar,temp
        endif
 10   continue
      call dumpit(m)
      call dimsav(m)
      return
c=======================================================================
 100  format(//1x,'x  ='1pe14.4,'  Profile'//5x,'i',3x,'y/delta',
     *       5x,'u/Ue',7x,'yplus',6x,'uplus',5x,'T/Te')
 110  format(i6,1p5e11.3)
 120  format(//1x,'x  ='1pe14.4,'  Profile'//5x,'i',3x,'y/delta',
     *       5x,'u/Ue',7x,'yplus',6x,'uplus',5x,'k/Ue**2',5x,'omega',
     *       5x,'eps/mu',5x,'L/delta',4x,'uv/tauw',5x,'T/Te')
 130  format(i6,1p10e11.3)
 140  format(//1x,'x  ='1pe14.4,'  Profile'//5x,'i',3x,'y/delta',
     *       5x,'u/Ue',7x,'yplus',6x,'uplus',6x,'chi',7x,'eps/mu',
     *       5x,'uv/tauw',5x,'T/Te')
 150  format(i6,1p8e11.3)
 160  format(//1x,'x  ='1pe14.4,'  Profile'//5x,'i',3x,'y/delta',
     *       5x,'u/Ue',7x,'yplus',6x,'uplus',5x,'eps/mu',5x,'uv/tauw',
     *       5x,'T/Te')
 170  format(i6,1p7e11.3)
c=======================================================================
      end
      subroutine propm(n)
c=======================================================================
c      S E T   O R   E X T R A P O L A T E   Q U A N T I T I E S
c                 F O R   E A C H   I T E R A T I O N
c=======================================================================
      include 'comeb.'
c=======================================================================
      rat2jm=1.
      rat2jp=1.
      if((itro+iturb).eq.0) then
c=======================================
c   Special computations for itro=0
c=======================================
        fm1=z4*fo(n)-z5*fn(n)
        tm1=z4*to(n)-z5*tn(n)
        if(tm1.le.-1.) tm1=.5*(to(n)+tn(n))
        vm1=z4*vo(n)-z5*vn(n)
        rato=z4*ratoo(n)-z5*raton(n)
        if(j.ne.0) then
          rat2jm=(.5*(z4*ratoo(n-1)-z5*raton(n-1)+rato))**2
          rat2jp=(.5*(z4*ratoo(n+1)-z5*raton(n+1)+rato))**2
        endif
        fy=y2(n)*(fo(n+1)-fo(n-1))
        ty=y2(n)*(to(n+1)-to(n-1))
c=======================================
c    Normal computations for itro > 0
c=======================================
      else
        fm1=fp(n)
        tm1=tp(n)
        vm1=vp(n)
        rato=ratop(n)
        if(j.ne.0) then
          rat2jm=(.5*(ratop(n-1)+rato))**2
          rat2jp=(.5*(ratop(n+1)+rato))**2
        endif
        fy=fz(n)
        ty=tz(n)
      endif
      tm1p=1.+tm1
      xlm1=xlm11(n)
      xlp=.5*(xlm1+xlm11(n+1))
      rxm=rat2jm*xlm
      rxp=rat2jp*xlp
      rato1j=rato**j
      rato2j=rato1j*rato1j
      fm2=z2*fo(n)-z3*fn(n)
      tm2=z2*to(n)-z3*tn(n)
      fym=fy-coef11*tm1p*fm1/rato1j
      fyp=fy+coef11*tm1p*fm1/rato1j
      return
      end
      subroutine rankin
c=======================================================================
c  C O M P U T E   F R E E S T R E A M   C O N D I T I O N S   F R O M
c      T H E   R A N K I N E - H U G O N I O T   R E L A T I O N S
c=======================================================================
      include 'comeb.'
c=======================================================================
      gm1=g-1.
      gm1og=gm1/g
      xmac=1.+.5*gm1*xma**2
      rt1=pt1/(r*tt1)
      p1=pt1/(xmac)**(1./gm1og)
      r1=rt1/(xmac)**(1./gm1)
      t1=tt1/xmac
      aa1=sqrt(g*p1/r1)
      u1=xma*aa1
      tref=gm1og*u1**2/r
      vis1=viscon*t1**vispow/(t1+su)
      visref=viscon*tref**vispow/(tref+su)
      rey=r1*u1/vis1
      reyref=rey*vis1/visref
      eps=1./sqrt(reyref)
      xwave=.0174533*wave
      abc=(xma*sin(xwave))**2
      if(abs(xwave).le..0000001) then
        p10=pt1/(r1*u1*u1)
      else
        p10=(1./(g*xma*xma))*((xmac*abc*(g+1.))
     *     /(abc*gm1+2.))**(1./gm1og)*((g+1.)/(2.*g*abc-gm1))**(1./gm1)
      endif
      t10=.5+1./(xma*xma*gm1)
      tc=su/(t1*xmac)
      prr=1./pr
      prtr=1./prt
      vispom=vispow-1.
      return
      end
      subroutine solvbb(justmu)
c=======================================================================
c       I N V E R T   T R I D I A G O N A L   M A T R I C E S :
c       B A L D W I N - B A R T H   M O D E L   S O L U T I O N
c-----------------------------------------------------------------------
c  Subroutines referenced:
c            PROPM........Mean-flow properties used repeatedly
c=======================================================================
      include 'comeb.'
c=======================================================================
      if(justmu.ne.0) go to 40
      etp(1)=0.
      pp1(1)=etp(1)
      pp2(1)=0.
      etpp(1)=etp(1)
c-----------------------------------------------------------------------
c                      F O R W A R D   S W E E P
c-----------------------------------------------------------------------
      iedgem=iedge-1
      tot1m=.5*(ry(1)+ry(2))
      tot2m=.5*(eps0(1)+eps0(2))
      do 10 n=2,iedgem
c=======================================
c  Compute quantities used repeatedly
c=======================================
        call propm(n)
        etm1=etp(n)
        etm2=z2*eto(n)-z3*etn(n)
c-----------------------------
c    Special computations: 
c       itro=0,iturb=1
c-----------------------------
        if(itro.eq.0) then
          etm1=z4*eto(n)-z5*etn(n)
          if(etm1.lt.0.) etm1=.5*(eto(n)+etn(n))
          etpp(n)=etm1
        endif
c=======================================
c     E D D Y   V I S C O S I T Y
c---------------------------------------
c   Compute centered viscosity terms
c        and the h coefficients
c=======================================
        tot1p=.5*(ry(n)+ry(n+1))
        tot2p=.5*(eps0(n)+eps0(n+1))
        epcon1=xlm1*tm1p*(1.+2.*ep(n)/sigmae)
        epcon2=1./tm1p
        epnm=epcon1*tot1m-epcon2*tot2m
        epnp=epcon1*tot1p-epcon2*tot2p
        epnh=xlm1*rato2j*(1.+ep(n)/sigmae)
        if(epnm.lt.0.) epnm=epnh
        if(epnp.lt.0.) epnp=epnh
        h1=.5*oy*fm1/dx2
        h2=y2(n)*vm1
        h7=2.*fm1
c=======================================
c        Compute source term
c=======================================
        prodn=(ceps2*ef2(n)-ceps1)*coef3*sqrt(cmu*efmu(n)/xlm1)
     *       *rato1j*fy/tm1p**2
c=======================================
c     Compute matrix coefficients
c=======================================
        a1n=-y3(n)*epnm-h2
        b1n= y3(n)*epnm+y1(n)*epnp+h7+z1*h1
        c1n=-y1(n)*epnp+h2
        g1n=h1*etm2
        if(prodn.gt.0.) then
          g1n=g1n+prodn*etm1
        else
          b1n=b1n-prodn
        endif
c=======================================
c  Do first step of inversion procedure
c=======================================
        b1ns=b1n+a1n*pp2(n-1)
        g1ns=g1n-a1n*pp1(n-1)
        pp1(n)=g1ns/b1ns
        pp2(n)=-c1n/b1ns
        tot1m=tot1p
        tot2m=tot2p
 10   continue
c-----------------------------------------------------------------------
c                  B A C K   S U B S T I T U T I O N
c-----------------------------------------------------------------------
c=======================================
c  Set condition at B.L. edge and solve
c=======================================
      etp(iedge)=ete
      jprint=0
      kon=iedgem
      do 30 n=2,iedgem
        etp(kon)=pp1(kon)+pp2(kon)*etp(kon+1)
        kon=kon-1
 30   continue
c-----------------------------------------------------------------------
c           C O M P U T E   E D D Y   V I S C O S I T Y
c-----------------------------------------------------------------------
 40   twall=1.+tp(1)
      fz(1)=(-ww1*fp(1)+ww2*fp(2)-ww3*fp(3)+ww4*fp(4))/(ww5*dy(1))
      dudyw=abs(fz(1))
      ypscal=eps*sqrt(coef1*xlm11(1)*twall*dudyw)
      do 50 n=2,iedge
        if(etp(n).lt.0.) etp(n)=.5*(etpp(n-1)+etpp(n+1))
        eltsq=xlm11(n)*(1.+tp(n))**2
        yplus(n)=ypscal*xn(n)/eltsq
        if(yplus(n).lt..2) then
          efmu(n)=cbb3*(1.-cbb4*yplus(n))*yplus(n)**2
          ef2(n)=cbb1+cbb5*(1.-cbb6*yplus(n))
        else
          exp1=exp(-yplus(n)/a0plus)
          exp2=exp(-yplus(n)/a2plus)
          damp1=1.-exp1
          damp2=1.-exp2
          efmu(n)=damp1*damp2
          rootef=sqrt(efmu(n))
          ef2(n)=cbb1+cbb2*(1./(akappa*yplus(n))+efmu(n))*(rootef
     *     +yplus(n)*(damp2*exp1/a0plus+damp1*exp2/a2plus)/rootef)
        endif
        ef1(n)=cmu*coef69*efmu(n)*etp(n)
        ep(n)=ef1(n)/eltsq
        ry(n)=ratop(n)**(2*j)/(1.+tp(n))
        eps0(n)=eltsq*ry(n)*ep(n)/sigmae
 50   continue
      ef1(1)=0.
      ef2(1)=cbb1+cbb5
      ep(1)=0.
      eps0(1)=0.
      ry(1)=1./twall
      ef1(iedge+1)=ef1(iedge)
      ef2(iedge+1)=ef2(iedge)
      ep(iedge+1)=ep(iedge)
      eps0(iedge+1)=eps0(iedge)
      ry(iedge+1)=ry(iedge)
      return
      end
      subroutine solvbl
c=======================================================================
c  B A L D W I N - L O M A X   M O D E L   E D D Y   V I S C O S I T Y
c=======================================================================
      include 'comeb.'
c=======================================================================
c      Compute scaling factors
c=======================================
      twall=1.+tp(1)
      fz(1)=(-ww1*fp(1)+ww2*fp(2)-ww3*fp(3)+ww4*fp(4))/(ww5*dy(1))
      dudyw=abs(fz(1))
      ypscal=eps*sqrt(coef1*dudyw/(xlm11(1)*twall**3))
c=======================================
c        Compute mixing length
c        and velocity gradient
c=======================================
      fmax=0.
      ymax=0.
      do 10 n=2,iedge
        yplus(n)=ypscal*xn(n)
        ef1(n)=akappa*eps*xn(n)*(1.-exp(-yplus(n)/a0plus))
        fbar=ef1(n)*ratop(n)**j*fz(n)/(1.+tp(n))
        if(fbar.gt.fmax) then
          fmax=fbar
          ymax=eps*xn(n)
        endif
 10   continue
      fmax=coef12*fmax/akappa
      fwake=ymax*fmax
      fwake2=cwk*ymax*ue**2/fmax
      if(fwake2.lt.fwake) fwake=fwake2
      outer=alphac*ccp*coef10*fwake
      xndens=ymax/(ckleb*eps)
c=======================================
c       Compute eddy viscosity
c=======================================
      muflag=0
      do 20 n=1,iedge
c-----------------------------
c   Outer layer viscosity
c-----------------------------
        ep(n)=outer/((1.+5.5*(xn(n)/xndens)**6)*xlm11(n)*(1.+tp(n))**2)
c-----------------------------
c   Inner layer viscosity
c-----------------------------
        if(muflag.eq.0) then
          emuti=coef1*ratop(n)**j*ef1(n)**2*fz(n)
     *                    /(xlm11(n)*(1.+tp(n))**3)
c-----------------------------
c  Select smaller of the two
c-----------------------------
          if(emuti.lt.ep(n)) then
            ep(n)=emuti
          else
            muflag=1
          endif
        endif
 20   continue
      return
      end
      subroutine solvcs
c=======================================================================
c   C E B E C I - S M I T H   M O D E L   E D D Y   V I S C O S I T Y
c=======================================================================
      include 'comeb.'
c=======================================================================
c      Compute velocity thickness
c        and velocity gradient
c=======================================
      twall=1.+tp(1)
      co=twall
      deltav=0.
      fz(1)=(-ww1*fp(1)+ww2*fp(2)-ww3*fp(3)+ww4*fp(4))/(ww5*dy(1))
      do 10 n=2,iedge
        c=(1.+tp(n))*(1.-fp(n))/ratop(n)**j
        deltav=deltav+.5*(co+c)*dy(n-1)
        co=c
        fz(n)=y4(n)*fp(n+1)-y5(n)*fp(n)-y6(n)*fp(n-1)
 10   continue
c=======================================
c      Compute scaling factors
c=======================================
      dudyw=abs(fz(1))
      ypscal=eps*sqrt(coef1*dudyw/(xlm11(1)*twall**3))
      outer=alphac*coef3*deltav
      utaub=sqrt(coef4*xlm11(1)*twall*dudyw)
      pplus=coef6/utaub**3
      vwplus=rvwald*twall/(utaub*re*r1*u1)
c=======================================
c       Compute eddy viscosity
c=======================================
      muflag=0
      do 20 n=1,iedge
c-----------------------------
c   Outer layer viscosity
c-----------------------------
        ep(n)=outer/((1.+5.5*(xn(n)/xnden)**6)*xlm11(n)*(1.+tp(n))**2)
c-----------------------------
c   Inner layer viscosity
c-----------------------------
        if(muflag.eq.0) then
          yplus(n)=ypscal*xn(n)
          aplus=a0plus*xlm11(n)/xlm11(1)*((1.+tp(n))/twall)**1.5
          if(vwplus.ne.0.) then
            expon=exp(11.8*xlm11(1)*twall*vwplus/(xlm11(n)*(1.+tp(n))))
            ensqrd=xlm11(n)*twall**2*(1.-expon)*pplus/vwplus+expon
          else
            ensqrd=1.-11.8*xlm11(1)*twall**3*pplus
          endif
          aplus=aplus/sqrt(ensqrd)
          elmix=akappa*eps*xn(n)*(1.-exp(-yplus(n)/aplus))
          emuti=coef1*ratop(n)**j*elmix**2*fz(n)
     *                 /(xlm11(n)*(1.+tp(n))**3)
c-----------------------------
c  Select smaller of the two
c-----------------------------
          if(emuti.lt.ep(n)) then
            ep(n)=emuti
          else
            muflag=1
          endif
        endif
 20   continue
      return
      end
      subroutine solvft(isep)
c=======================================================================
c  M O M E N T U M  /  E N E R G Y   E Q U A T I O N   S O L U T I O N
c=======================================================================
      include 'comeb.'
c=======================================================================
c    Set surface boundary conditions
c=======================================
c          No-slip
c-----------------------------
      fp(1)=0.
      pp1(1)=fp(1)
      pp2(1)=0.
      pp3(1)=0.
c-----------------------------
c  Wall temperature specified
c-----------------------------
      if(kodwal.eq.1) then
        tp(1)=tw/te-1.
c-----------------------------
c    Heat flux specified
c-----------------------------
      else
        tz(1)=tzmu/xlm1p
      endif
      qq1(1)=tp(1)
      qq2(1)=0.
      qq3(1)=0.
c-----------------------------------------------------------------------
c                      F O R W A R D   S W E E P
c-----------------------------------------------------------------------
      iedgem=iedge-1
      epsm=.5*(ep(1)+ep(2))
      xlm=.5*(xlm11(1)+xlm11(2))
      do 10 n=2,iedgem
        epsp=.5*(ep(n)+ep(n+1))
c=======================================
c  Compute quantities used repeatedly
c=======================================
        call propm(n)
c=======================================
c   Compute centered viscosity terms
c        and the h coefficients
c=======================================
        epfm=rxm*(1.+epsm)
        epfp=rxp*(1.+epsp)
        eptm=rxm*(prr+prtr*epsm)
        eptp=rxp*(prr+prtr*epsp)
        h1=.5*oy*fm1/dx2
        h2=vm1*y2(n)
        h3=xal*rato2j*xlm1*fy
        if(model.gt.2) h3=h3*(1.+ep(n))
        h4=2.*y2(n)*h3
        h5=xbe*fm1
        h6=xal*dissk(n)
        if(model.eq.1) then
          delg1=2.*y2(n)*(deluv(n+1)-deluv(n))
        else
          delg1=0.
        endif
c=======================================
c     Compute matrix coefficients
c=======================================
        a1=-y3(n)*epfm-h2
        b1= y3(n)*epfm+y1(n)*epfp+h5+z1*h1
        c1=-y1(n)*epfp+h2
        d1=0.
        e1=-xbe
        f1=0.
        g1=h1*fm2+xbe+delg1
        a2= h4
        b2=0.
        c2=-h4
        d2=-y3(n)*eptm-h2
        e2= y3(n)*eptm+y1(n)*eptp+z1*h1
        f2=-y1(n)*eptp+h2
        g2=h1*tm2-h3*fy+h6
c=======================================
c  Do first step of inversion procedure
c=======================================
c  Special computations for
c surface heat flux specified
c-----------------------------
        if(kodwal.ne.1.and.n.eq.2) then
          did=(c2*d1-c1*d2)-((c2*f1-c1*f2)*(((1.+xk)**2)-1.))
          qq1(1)=((c2*g1-c1*g2)
     *           +(c2*f1-c1*f2)*(xk*(1+xk)*dy(1))*tz(1))/did
          qq2(1)=-(c2*b1-c1*b2)/did
          qq3(1)=-((c2*e1-c1*e2)+((c2*f1-c1*f2)*((1.+xk)**2)))/did
        endif
c-----------------------------
c    Normal computations
c-----------------------------
        b1s=b1+a1*pp2(n-1)+d1*qq2(n-1)
        b2s=b2+a2*pp2(n-1)+d2*qq2(n-1)
        e1s=e1+a1*pp3(n-1)+d1*qq3(n-1)
        e2s=e2+a2*pp3(n-1)+d2*qq3(n-1)
        g1s=g1-a1*pp1(n-1)-d1*qq1(n-1)
        g2s=g2-a2*pp1(n-1)-d2*qq1(n-1)
        d=1./(b1s*e2s-e1s*b2s)
        pp1(n)=d*(g1s*e2s-g2s*e1s)
        pp2(n)=d*(e1s*c2-c1*e2s)
        pp3(n)=d*(e1s*f2-f1*e2s)
        qq1(n)=d*(b1s*g2s-b2s*g1s)
        qq2(n)=d*(c1*b2s-b1s*c2)
        qq3(n)=d*(f1*b2s-b1s*f2)
        epsm=epsp
        xlm=xlp
 10   continue
c-----------------------------------------------------------------------
c                  B A C K   S U B S T I T U T I O N
c-----------------------------------------------------------------------
c=======================================
c     Set conditions at B.L. edge
c=======================================
      fp(iedge)=1.
      tp(iedge)=0.
      kon=iedgem
      do 20 n=2,iedgem
c=======================================
c        Update the velocity
c=======================================
        fp(kon)=pp1(kon)+pp2(kon)*fp(kon+1)+pp3(kon)*tp(kon+1)
c-----------------------------
c    Check for separation
c-----------------------------
        if(fp(kon).lt.0.) then
          isep=1
          go to 30
c-----------------------------
c    Do a little smoothing
c-----------------------------
        elseif(fp(kon).gt.1.) then
          fp(kon)=1.
        endif
c=======================================
c       Update the temperature
c=======================================
        tp(kon)=qq1(kon)+qq2(kon)*fp(kon+1)+qq3(kon)*tp(kon+1)
        kon=kon-1
 20   continue
      if(kodwal.ne.1) then
        tp(1)=(xk*(1.+xk)*dy(1)*tz(1)-(1.+xk)**2*tp(2)+tp(3))
     *       /(1.-(1.+xk)**2)
        tw=(1.+tp(1))*te
      endif
c=======================================
c   Reset velocity and temperature if
c    the boundary layer separates
c=======================================
 30   if(isep.ne.0) then
        do 40 n=1,iedge
          fp(n)=fo(n)
          tp(n)=to(n)
 40     continue
      endif
      return
      end
      subroutine solvjk
c=======================================================================
c       I N V E R T   T R I D I A G O N A L   M A T R I C E S :
c        J O H N S O N - K I N G   M O D E L   S O L U T I O N
c=======================================================================
      include 'comeb.'
c=======================================================================
c      Compute velocity thickness
c        and velocity gradient
c=======================================
      twall=1.+tp(1)
      co=twall
      deltav=0.
      fz(1)=(-ww1*fp(1)+ww2*fp(2)-ww3*fp(3)+ww4*fp(4))/(ww5*dy(1))
      do 10 n=2,iedge
        c=(1.+tp(n))*(1.-fp(n))/ratop(n)**j
        deltav=deltav+.5*(co+c)*dy(n-1)
        co=c
        fz(n)=y4(n)*fp(n+1)-y5(n)*fp(n)-y6(n)*fp(n-1)
 10   continue
c=======================================
c      Compute scaling factors
c=======================================
      dudyw=abs(fz(1))
      utaub=sqrt(coef4*xlm11(1)*twall*dudyw)
      ypscal=eps*sqrt(coef1*dudyw/(xlm11(1)*twall**3))
      ym=xn(nmaxjk)
      tm=1.+tp(nmaxjk)
      umax=sqrt(uvmaxp)
      umplus=ue*umax/utaub
      outer=alphac*coef3*deltav
      udplus=1.
      if(umplus.gt.1.) udplus=umplus
      elm=akappa*ym
      elm2=cjk1*xnden
      if(elm2.lt.elm) elm=elm2
      elc=elm/(1.+umplus*sqrt(twall/tm))
c=======================================
c  Compute equilibrium eddy viscosity
c=======================================
      do 20 n=1,iedge
c-----------------------------
c   Outer layer viscosity
c-----------------------------
        ef1(n)=outer/((1.+5.5*(xn(n)/xnden)**6)*xlm11(n)*(1.+tp(n))**2)
c-----------------------------
c   Inner layer viscosity
c-----------------------------
        gamma2=tanh(xn(n)/elc)
        usplus=(1.-gamma2)*sqrt((1.+tp(n))/twall)
     *         +gamma2*umplus*sqrt((1.+tp(n))/tm)
        yplus(n)=ypscal*xn(n)
        fact1=xlm11(1)/xlm11(n)*(twall/(1.+tp(n)))**2
        damp=1.-exp(-fact1*udplus*yplus(n)/a0plus)
        ef2(n)=akappa*usplus*fact1*damp**2*yplus(n)
 20   continue
c=======================================
c  Determine nonequilibrium parameter
c=======================================
      uveq=0.
      do 30 n=2,iedge
        ep(n)=ef1(n)*tanh(ef2(n)/ef1(n))
        uvbar=xlm11(n)*(1.+tp(n))*ratop(n)**j*ep(n)*fz(n)/coef3
        if(uvbar.gt.uveq) then
          uveq=uvbar
        else
          go to 40
        endif
 30   continue
 40   nmaxjk=n-1
      sigmaj=uvmaxp/uveq
c=======================================
c       Recompute eddy viscosity
c           and shear stress
c=======================================
      do 50 n=1,iedge
        ep(n)=sigmaj*ef1(n)*tanh(ef2(n)/(sigmaj*ef1(n)))
        uvp(n)=xlm11(n)*(1.+tp(n))*ratop(n)**j*ep(n)*fz(n)/coef3
 50   continue
c=======================================
c  Solve equation for maximum stress
c=======================================
      fm1=fp(nmaxjk)
      h1=.5*oy*fm1/dx2
      h6=2.*xbe*fm1
      umaxeq=sqrt(uveq)
      uvm2=z2*uvmaxo-z3*uvmaxn
      b1=h1*z1+h6
      g1=h1*uvm2
      prodjk=ajk1*coef14*(umaxeq-umax)/elm
      if(sigmaj.gt.1.) prodjk=prodjk-cdif*coef14*umax
     *                       *(1.-sqrt(sigmaj))/(cjk2*xnden-ym)
      if(prodjk.gt.0.) then
        g1=g1+prodjk*uvmaxp
      else
        b1=b1-prodjk
      endif
      uvmaxp=g1/b1
      uvmaxp=relax1*uvmxpp+relax2*uvmaxp
      return
      end
      subroutine solvke(justmu)
c=======================================================================
c       I N V E R T   T R I D I A G O N A L   M A T R I C E S :
c   K - E P S I L O N   M O D E L   E Q U A T I O N   S O L U T I O N
c-----------------------------------------------------------------------
c  Subroutines referenced:
c            EFMU12.......Viscous damping terms
c            PROPM........Mean-flow properties used repeatedly
c            STABLE.......Matrix diagonal-dominance routine
c            WBCKE........Surface boundary condition for k and epsilon
c=======================================================================
      include 'comeb.'
c=======================================================================
      if(justmu.ne.0) go to 50
      call wbcke(etp(1),epp(1))
      pp1(1)=etp(1)
      pp2(1)=0.
      etpp(1)=etp(1)
      qq1(1)=epp(1)
      qq2(1)=0.
      eppp(1)=epp(1)
c-----------------------------------------------------------------------
c                      F O R W A R D   S W E E P
c-----------------------------------------------------------------------
      netegm=netaeg-1
      epsm=.5*(ep(1)+ep(2))
      xlm=.5*(xlm11(1)+xlm11(2))
      do 10 n=2,netegm
        epsp=.5*(ep(n)+ep(n+1))
c=======================================
c  Compute quantities used repeatedly
c=======================================
        call propm(n)
        etm1=etp(n)
        epm1=epp(n)
        etm2=z2*eto(n)-z3*etn(n)
        epm2=z2*epo(n)-z3*epn(n)
c-----------------------------
c    Special computations: 
c       itro=0,iturb=1
c-----------------------------
        if(itro.eq.0) then
          etm1=z4*eto(n)-z5*etn(n)
          if(etm1.lt.0.) etm1=.5*(eto(n)+etn(n))
          epm1=z4*epo(n)-z5*epn(n)
          if(epm1.lt.0.) epm1=.5*(epo(n)+epn(n))
          etpp(n)=etm1
          eppp(n)=epm1
        endif
c=======================================
c  T U R B U L E N C E   K I N E T I C
c             E N E R G Y
c---------------------------------------
c   Compute centered viscosity terms
c        and the h coefficients
c=======================================
        epkm=rxm*(1.+epsm/sigmak)
        epkp=rxp*(1.+epsp/sigmak)
        h1=.5*oy*fm1/dx2
        h2=vm1*y2(n)
        h6=2.*xbe*fm1
c=======================================
c Compute production & dissipation terms
c=======================================
        prodk=coef69*cmu*efmu(n)*rato2j*etm1/epm1*(fym/tm1p)**2
        dissk(n)=coef9*(epm1+eps0(n))/etm1
c=======================================
c    Compressibility modification
c=======================================
        xmt2=2.*xmae**2*etm1/tm1p-xmtsq
        if(xmt2.gt.0.) dissk(n)=dissk(n)*(1.+xis*xmt2)
c=======================================
c       Make sure the matrix is
c         diagonally dominant
c=======================================
        call stable(prodk,dissk(n),psik,0)
c=======================================
c     Compute matrix coefficients
c=======================================
        a1k=-y3(n)*epkm-h2
        b1k= y3(n)*epkm+y1(n)*epkp+h6+z1*h1-prodk+dissk(n)+psik
        c1k=-y1(n)*epkp+h2
        g1k=h1*etm2+psik*etm1
c=======================================
c  Do first step of inversion procedure
c=======================================
        b1ks=b1k+a1k*pp2(n-1)
        g1ks=g1k-a1k*pp1(n-1)
        pp1(n)=g1ks/b1ks
        pp2(n)=-c1k/b1ks
c=======================================
c   D I S S I P A T I O N   R A T E
c---------------------------------------
c   Compute centered viscosity terms
c        and the h coefficients
c=======================================
        epem=rxm*(1.+epsm/sigmae)
        epep=rxp*(1.+epsp/sigmae)
        h7=coef8*fm1
c=======================================
c Compute production & dissipation terms
c=======================================
        prode=(cmu*ceps1*ef1(n)*efmu(n)*coef69*rato2j*etm1
     *       *(fym/tm1p)**2+ehat(n))/epm1
        disse=ceps2*ef2(n)*coef9*epm1/etm1
c=======================================
c     Compute matrix coefficients
c=======================================
        a1e=-y3(n)*epem-h2
        b1e=y3(n)*epem+y1(n)*epep+h7+z1*h1+disse+psieps*prode
        c1e=-y1(n)*epep+h2
        g1e=h1*epm2+(1.+psieps)*prode*epm1
c=======================================
c  Do first step of inversion procedure
c=======================================
        b1es=b1e+a1e*qq2(n-1)
        g1es=g1e-a1e*qq1(n-1)
        qq1(n)=g1es/b1es
        qq2(n)=-c1e/b1es
        epsm=epsp
        xlm=xlp
 10   continue
c-----------------------------------------------------------------------
c                  B A C K   S U B S T I T U T I O N
c-----------------------------------------------------------------------
c=======================================
c Set conditions at B.L. edge and solve
c=======================================
      etp(netaeg)=ete
      epp(netaeg)=epe
      jprint=0
      kprint=0
      kon=netegm
      ktol=kon-4
      nedge=netaeg
      do 30 n=2,netegm
        etp(kon)=pp1(kon)+pp2(kon)*etp(kon+1)
        epp(kon)=qq1(kon)+qq2(kon)*epp(kon+1)
c=======================================
c Check for sharp turbulent-nonturbulent
c  interface unless the flow is laminar
c=======================================
        if(kon.ge.ktol) then
          if(epp(kon).lt.0..or.etp(kon).lt.0.) then
            if(ibound.ne.0) nedge=kon
            do 20 ki=kon,netegm
              epp(ki)=epe
              etp(ki)=ete
 20         continue
          endif
          if(epp(kon).lt.0.) jprint=jprint+1
          if(etp(kon).lt.0.) kprint=kprint+1
        endif
        kon=kon-1
 30   continue
c=======================================
c     Recompute netaeg if necessary
c=======================================
      if(nedge.ne.netaeg) netaeg=nedge
      do 40 n=netaeg,iedge
        etpp(n)=ete
        eppp(n)=epe
        etp(n)=ete
        epp(n)=epe
 40   continue
c-----------------------------------------------------------------------
c           C O M P U T E   E D D Y   V I S C O S I T Y
c-----------------------------------------------------------------------
 50   con=oy*re/(xnue*eps**2)
      conry=re*ue/(xnue*eps)
      fz(1)=(-ww1*fp(1)+ww2*fp(2)-ww3*fp(3)+ww4*fp(4))/(ww5*dy(1))
      xmf4=eps*xnue*rmi**j*xlm11(1)*fz(1)/oz
      uplusb=ue*sqrt(xmf4*(1.+tp(1)))
      tscale=re*uplusb/(eps*xnue*xlm11(1)*(1.+tp(1))**2)
      do 60 n=2,iedge
        if(epp(n).ne.0.) then
          if(etp(n).lt.0.) then
            if(n.lt.4) then
              etp(n)=etpp(4)*(xn(n)/xn(4))**2
            else
              etp(n)=.5*(etpp(n-1)+etpp(n+1))
            endif
          endif
          if(epp(n).lt.0.) epp(n)=.5*(eppp(n-1)+eppp(n+1))
          ret(n)=con*(etp(n)/(1.+tp(n)))**2/(epp(n)*xlm11(n))
          ry(n)=conry*xn(n)*sqrt(etp(n))/(xlm11(n)*(1.+tp(n))**2)
        else
          ret(n)=0.
          ry(n)=0.
        endif
        yplus(n)=tscale*xn(n)
        if(justmu.eq.0) then
          yplus(n)=.5*(yplus(n)+tscale*xn(n))
        else
          yplus(n)=tscale*xn(n)
        endif
        call efmu12(n)
        ep(n)=cmu*efmu(n)*ret(n)
        dissk(n)=etp(n)*dissk(n)
 60   continue
      ep(1)=0.
      ry(1)=0.
      ret(1)=0.
      ep(iedge+1)=ep(iedge)
      ry(iedge+1)=ry(iedge)
      ret(iedge+1)=ret(iedge)
      return
      end
      subroutine solvkw(justmu)
c=======================================================================
c  T U R B U L E N C E - M O D E L   E Q U A T I O N   S O L U T I O N
c         K - O M E G A / S T R E S S - O M E G A   M O D E L S
c-----------------------------------------------------------------------
c  Subroutines referenced:
c            PROPM........Mean-flow properties used repeatedly
c            STABLE.......Matrix diagonal-dominance routine
c            WBCKW........Surface boundary condition for k and omega
c=======================================================================
      include 'comeb.'
c=======================================================================
      if(justmu.ne.0) go to 50
      call wbckw(etp(1),wtp(1))
      pp1(1)=etp(1)
      pp2(1)=0.
      etpp(1)=etp(1)
      qq1(1)=wtp(1)
      qq2(1)=0.
      wtpp(1)=wtp(1)
c-----------------------------------------------------------------------
c                      F O R W A R D   S W E E P
c-----------------------------------------------------------------------
      netegm=netaeg-1
      epsm=.5*(ep(1)+ep(2))
      xlm=.5*(xlm11(1)+xlm11(2))
      do 10 n=2,netegm
        epsp=.5*(ep(n)+ep(n+1))
c=======================================
c      Prescribe dissipation rate
c      analytically near the wall
c=======================================
        if(n.lt.nskip) then
          yplus(n)=yscale*xn(n)
          yr2=(epsr+rghcon*yplus(n))**2
          wtp(n)=wwall/yr2
          qq1(n)=wtp(n)
          qq2(n)=0.
        endif
c=======================================
c  Compute quantities used repeatedly
c=======================================
        call propm(n)
        etm1=etp(n)
        wtm1=wtp(n)
        etm2=z2*eto(n)-z3*etn(n)
        wtm2=z2*wto(n)-z3*wtn(n)
c-----------------------------
c Special computations, itro=0
c-----------------------------
        if(itro.eq.0) then
          etm1=z4*eto(n)-z5*etn(n)
          if(etm1.lt.0.) etm1=.5*(eto(n)+etn(n))
          if(n.ge.nskip) then
            wtm1=z4*wto(n)-z5*wtn(n)
            if(wtm1.le.0.) wtm1=.5*(wto(n)+wtn(n))
          endif
          etpp(n)=etm1
          wtpp(n)=wtm1
        endif
c=======================================
c  T U R B U L E N C E   K I N E T I C
c             E N E R G Y
c---------------------------------------
c   Compute centered viscosity terms
c        and the h coefficients
c=======================================
        epkm=rxm*(1.+sigmas*epsm)
        epkp=rxp*(1.+sigmas*epsp)
        h1=.5*oy*fm1/dx2
        h2=vm1*y2(n)
        h6=2.*xbe*fm1
c=======================================
c Compute production & dissipation terms
c=======================================
        dissk(n)=betas*coef9*wtm1
        if(model.eq.0) then
          prodk=coef69*(rato2j*fym/tm1p-coef15*rato1j*fm1)*fym
     *         /(wtm1*tm1p)
          prodk=prodk*(alfs0*rk+ret(n))/(rk+ret(n))
          pod=clim2*prodk/dissk(n)
          climit(n)=1.
          if(pod.gt.1.) then
            climit(n)=1./sqrt(pod)
            prodk=climit(n)*prodk
          endif
          epkm=rxm*(1.+sigmas*epsm/climit(n))
          epkp=rxp*(1.+sigmas*epsp/climit(n))
        else
          prodk=coef3*rato1j*uvp(n)*fym/(etm1*tm1p)
        endif
        if(nvisc.ne.0) dissk(n)=dissk(n)*(bbeta0+(ret(n)/rb)**4)
     *                                      /(1.+(ret(n)/rb)**4)
c=======================================
c     Cross-diffusion modification
c=======================================
        dkdeta=y4(n)*etp(n+1)-y5(n)*etp(n)-y6(n)*etp(n-1)
        dwdeta=y4(n)*wtp(n+1)-y5(n)*wtp(n)-y6(n)*wtp(n-1)
        dkdw=dkdeta*dwdeta
        vextra=0.
        if(dkdw.gt.0.) then
          vextra=sigmad*coef9*(coef5/tm1p)**2*rato2j*dkdeta/wtm1
        endif
c=======================================
c    Compressibility modification
c=======================================
        xmt2=2.*xmae**2*etm1/tm1p-xmtsq
        if(xmt2.gt.0.) dissk(n)=dissk(n)*(1.+xis*xmt2)
c=======================================
c       Make sure the matrix is
c         diagonally dominant
c=======================================
        call stable(prodk,dissk(n),psik,0)
c=======================================
c     Compute matrix coefficients
c=======================================
        a1k=-y3(n)*epkm-h2
        b1k= y3(n)*epkm+y1(n)*epkp+h6+z1*h1-prodk+dissk(n)+psik
        c1k=-y1(n)*epkp+h2
        g1k=h1*etm2+psik*etm1
c=======================================
c  Do first step of inversion procedure
c=======================================
        b1ks=b1k+a1k*pp2(n-1)
        g1ks=g1k-a1k*pp1(n-1)
        pp1(n)=g1ks/b1ks
        pp2(n)=-c1k/b1ks
c=======================================
c           S P E C I F I C   
c   D I S S I P A T I O N   R A T E
c---------------------------------------
c Skip these computations if dissipation
c     rate is specified analytically
c=======================================
        if(n.ge.nskip) then
          epwm=rxm*(1.+sigma*epsm/climit(n))
          epwp=rxp*(1.+sigma*epsp/climit(n))
          h7=coef7*fm1
c=======================================
c Compute production & dissipation terms
c=======================================
          if(model.eq.0) then
            prodw=alpha*coef69*rato2j*(fym/tm1p)**2/wtm1
            prodw=climit(n)*prodw
            prodw=prodw*(rw*alf0+ret(n))/(rw+ret(n))
          else
            prodw=alpha*prodk
            prodw=prodw*(rw*alf0+ret(n))/(rw+ret(n))
     *          /(3.*rw*alfs0+ret(n))*(3.*rw+ret(n))
          endif
          dissw=beta*coef9*wtm1
c=======================================
c    Compressibility modification
c=======================================
          xmt2=2.*xmae**2*etm1/tm1p-xmtsq
          if(xmt2.gt.0.) dissw=dissw*(1.-xis*betas*xmt2/beta)
c=======================================
c       Make sure the matrix is
c         diagonally dominant
c=======================================
          call stable(prodw,dissw,psiw,1)
c=======================================
c     Compute matrix coefficients
c=======================================
          h2=(vm1-vextra)*y2(n)
          a1w=-y3(n)*epwm-h2
          b1w=y3(n)*epwm+y1(n)*epwp+h7+z1*h1-prodw+dissw+psiw
          c1w=-y1(n)*epwp+h2
          g1w=h1*wtm2+psiw*wtm1
c=======================================
c  Do first step of inversion procedure
c=======================================
          b1ws=b1w+a1w*qq2(n-1)
          g1ws=g1w-a1w*qq1(n-1)
          qq1(n)=g1ws/b1ws
          qq2(n)=-c1w/b1ws
        endif
        epsm=epsp
        xlm=xlp
 10   continue
c-----------------------------------------------------------------------
c                  B A C K   S U B S T I T U T I O N
c-----------------------------------------------------------------------
c=======================================
c Set conditions at B.L. edge and solve
c=======================================
      etp(netaeg)=ete
      wtp(netaeg)=wte
      jprint=0
      kprint=0
      kon=netegm
      ktol=kon-4
      nedge=netaeg
      do 30 n=2,netegm
        etp(kon)=pp1(kon)+pp2(kon)*etp(kon+1)
        if(kon.ge.nskip) then
          wtp(kon)=qq1(kon)+qq2(kon)*wtp(kon+1)
c=======================================
c Check for sharp turbulent-nonturbulent
c  interface unless the flow is laminar
c=======================================
          if(kon.ge.ktol) then
            if(wtp(kon).lt.0..or.etp(kon).lt.0.) then
            if(ibound.ne.0) nedge=kon
              do 20 ki=kon,netegm
                wtp(ki)=wte
                etp(ki)=ete
 20           continue
            endif
            if(wtp(kon).lt.0.) jprint=jprint+1
            if(etp(kon).lt.0.) kprint=kprint+1
          endif
        endif
        kon=kon-1
 30   continue
c=======================================
c     Recompute netaeg if necessary
c=======================================
      if(nedge.ne.netaeg) netaeg=nedge
      do 40 n=netaeg,iedge
        etpp(n)=ete
        wtpp(n)=wte
        etp(n)=ete
        wtp(n)=wte
 40   continue
c-----------------------------------------------------------------------
c           C O M P U T E   E D D Y   V I S C O S I T Y
c-----------------------------------------------------------------------
 50   con=oy*re/(xnue*eps**2)
      do 60 n=2,iedge
        if(wtp(n).ne.0.) then
          if(etp(n).lt.0.) etp(n)=.5*(etpp(n-1)+etpp(n+1))
          if(wtp(n).lt.0.) wtp(n)=.5*(wtpp(n-1)+wtpp(n+1))
          ret(n)=con*etp(n)/(wtp(n)*xlm11(n)*(1.+tp(n))**2)
        else
          ret(n)=0.
        endif
        ep(n)=climit(n)*ret(n)*(alfs0*rk+ret(n))/(rk+ret(n))
        dissks(n)=dissk(n)
        dissk(n)=etp(n)*dissk(n)
 60   continue
      ep(1)=0.
      ret(1)=0.
      ep(iedge+1)=ep(iedge)
      ret(iedge+1)=ret(iedge)
      return
      end
      subroutine solvsa(justmu)
c=======================================================================
c       I N V E R T   T R I D I A G O N A L   M A T R I C E S :
c    S P A L A R T - A L L M A R A S   M O D E L   S O L U T I O N
c-----------------------------------------------------------------------
c  Subroutines referenced:
c            PROPM........Mean-flow properties used repeatedly
c=======================================================================
      include 'comeb.'
c=======================================================================
      if(justmu.ne.0) go to 40
      etp(1)=0.
      pp1(1)=etp(1)
      pp2(1)=0.
      etpp(1)=etp(1)
c-----------------------------------------------------------------------
c                      F O R W A R D   S W E E P
c-----------------------------------------------------------------------
      iedgem=iedge-1
      tot1m=.5*(ry(1)+ry(2))
      tot2m=.5*(eps0(1)+eps0(2))
      do 10 n=2,iedgem
c=======================================
c  Compute quantities used repeatedly
c=======================================
        call propm(n)
        etm1=etp(n)
        etm2=z2*eto(n)-z3*etn(n)
c-----------------------------
c    Special computations: 
c       itro=0,iturb=1
c-----------------------------
        if(itro.eq.0) then
          etm1=z4*eto(n)-z5*etn(n)
          if(etm1.lt.0.) etm1=.5*(eto(n)+etn(n))
          etpp(n)=etm1
        endif
c=======================================
c     E D D Y   V I S C O S I T Y
c---------------------------------------
c   Compute centered viscosity terms
c        and the h coefficients
c=======================================
        tot1p=.5*(ry(n)+ry(n+1))
        tot2p=.5*(eps0(n)+eps0(n+1))
        epcon1=1./(sigmae*tm1p)
        epcon2=cb2*coef69*etm1/(sigmae*tm1p)
        epnm=epcon1*tot1m-epcon2*tot2m
        epnp=epcon1*tot1p-epcon2*tot2p
        epnh=xlm1*rato2j*(1.+coef69*etm1/(xlm1*tm1p**2))/sigmae
        if(epnm.lt.0.) epnm=epnh
        if(epnp.lt.0.) epnp=epnh
        h1=.5*oy*fm1/dx2
        h2=y2(n)*vm1
        h7=2.*fm1
c=======================================
c        Compute source term
c=======================================
        prodn=cb1*(coef3*rato1j*fy/tm1p
     *       +ef2(n)/akappa**2*coef16*etm1/(eps*xn(n))**2)
        dissn=cw1*ef1(n)*coef16*etm1/(eps*xn(n))**2
c=======================================
c     Compute matrix coefficients
c=======================================
        a1n=-y3(n)*epnm-h2
        b1n= y3(n)*epnm+y1(n)*epnp+h7+z1*h1+dissn
        c1n=-y1(n)*epnp+h2
        g1n=h1*etm2
        if(prodn.gt.0.) then
          g1n=g1n+prodn*etm1
        else
          b1n=b1n-prodn
        endif
c=======================================
c  Do first step of inversion procedure
c=======================================
        b1ns=b1n+a1n*pp2(n-1)
        g1ns=g1n-a1n*pp1(n-1)
        pp1(n)=g1ns/b1ns
        pp2(n)=-c1n/b1ns
        tot1m=tot1p
        tot2m=tot2p
 10   continue
c-----------------------------------------------------------------------
c                  B A C K   S U B S T I T U T I O N
c-----------------------------------------------------------------------
c=======================================
c  Set condition at B.L. edge and solve
c=======================================
      etp(iedge)=ete
      kon=iedgem
      do 30 n=2,iedgem
        etp(kon)=pp1(kon)+pp2(kon)*etp(kon+1)
        kon=kon-1
 30   continue
c-----------------------------------------------------------------------
c           C O M P U T E   E D D Y   V I S C O S I T Y
c-----------------------------------------------------------------------
 40   twall=1.+tp(1)
      do 50 n=2,iedge
        if(etp(n).lt.0.) etp(n)=.5*(etpp(n-1)+etpp(n+1))
        etp(n)=relax1*etpp(n)+relax2*etp(n)
        fz(n)=y4(n)*fp(n+1)-y5(n)*fp(n)-y6(n)*fp(n-1)
        chi=coef69*etp(n)/(xlm11(n)*(1.+tp(n))**2)
        fv1=chi**3/(chi**3+csa1)
        fv2=1.-chi/(1.+fv1*chi)
        rbar=etp(n)/(fv2*etp(n)+coef17*ratop(n)**j*xn(n)**2
     *                                  *fz(n)/(1.+tp(n)))
        gbar=rbar+cw2*(rbar**6-rbar)
        fw=gbar*((1.+csa2)/(gbar**6+csa2))**(1./6.)
        ep(n)=fv1*chi
        ef1(n)=fw
        ef2(n)=fv2
        eps0(n)=ratop(n)**(2*j)/(1.+tp(n))
        ry(n)=xlm11(n)*ratop(n)**(2*j)*(1.+tp(n))*(1.+(1.+cb2)*chi)
 50   continue
      ep(1)=0.
      eps0(1)=1./twall
      ef1(1)=1.
      ef2(1)=1.
      ry(1)=xlm11(1)*twall
      ep(iedge+1)=ep(iedge)
      eps0(iedge+1)=eps0(iedge)
      ef1(iedge+1)=ef1(iedge)
      ef2(iedge+1)=ef2(iedge)
      ry(iedge+1)=ry(iedge)
      return
      end
      subroutine solvsw
c=======================================================================
c  S O L V E   T H E   R E Y N O L D S - S T R E S S   E Q U A T I O N S
c=======================================================================
      dimension aa(3,3)
      include 'comeb.'
c=======================================================================
      pp1(1)=0.
      pp2(1)=0.
      pp3(1)=0.
      do 10 ii=1,3
        do 10 jj=1,3
          qm(ii,jj,1)=0.
 10   continue
      iedgem=iedge-1
      epsm=.5*(ep(1)+ep(2))
      do 30 n=2,iedgem
        epsp=.5*(ep(n)+ep(n+1))
c=======================================
c  Compute quantities used repeatedly
c=======================================
        itrosv=itro
        itro=1
        call propm(n)
        itro=itrosv
        uum2=z2*uuo(n)-z3*uun(n)
        uvm2=z2*uvo(n)-z3*uvn(n)
        vvm2=z2*vvo(n)-z3*vvn(n)
c=======================================
c   Compute centered viscosity terms
c        and the h coefficients
c=======================================
        epkm=rxm*(1.+sigmas*epsm)
        epkp=rxp*(1.+sigmas*epsp)
        h1=.5*oy*fm1/dx2
        h2=2.*xbe*fm1
        h3=h1*z1+h2
        h4=dissks(n)
        h5=coef3*rato1j/tm1p
        h6=h1*uvm2
        h7=h1*uum2
        h8=h1*vvm2
        h9=vm1*y2(n)+epkm*y3(n)
        ha=vm1*y2(n)-epkp*y1(n)
c=======================================
c  Compute viscous modification terms
c=======================================
        if(nvisc.eq.0) then
          fvisc1=1.
          fvisc2=1.
          fvisc3=1.
        else
          retorb=(ret(n)/rb)**4
          fvisc1=(5./3.+retorb)/(1.+retorb)
          fvisc2=retorb/(1.+retorb)
          fvisc3=(gamh0+retorb)/(1.+retorb)
        endif
c=======================================
c   Compute additional coefficients
c=======================================
        a1=h3+h4*cc1*fvisc1+epkp*y1(n)+epkm*y3(n)
        b1=h5*fym
        c1=(1.+csw6*fvisc2)*coef13*fm1
        c2=2.*c1
c=======================================
c     Compute production terms
c=======================================
        pxx=csw4*fvisc2*b1+c2
        pyy=csw5*fvisc2*b1-c2
        ptx=-betah*fvisc2*b1-c1
        pty=csw1*fvisc2*b1+c1
        ptk=(csw3*fvisc2+.5*gammah*fvisc3)*b1
c=======================================
c       Set up the 3x3 matrix
c=======================================
        aa(1,1)=-h9*qm(1,1,n-1)+a1
        aa(1,2)=-h9*qm(1,2,n-1)-ptx
        aa(1,3)=-h9*qm(1,3,n-1)-pty
        aa(2,1)=-h9*qm(2,1,n-1)-pxx
        aa(2,2)=-h9*qm(2,2,n-1)+a1
        aa(2,3)=-h9*qm(2,3,n-1)
        aa(3,1)=-h9*qm(3,1,n-1)-pyy
        aa(3,2)=-h9*qm(3,2,n-1)
        aa(3,3)=-h9*qm(3,3,n-1)+a1
c=======================================
c      Adjust the r.h.s. vector
c=======================================
        pp1n=h9*pp1(n-1)+h6+ptk*etp(n)
        pp2n=h9*pp2(n-1)+h7
        pp3n=h9*pp3(n-1)+h8
c=======================================
c   Invert the matrix and compute the
c      new influence coefficients
c=======================================
        call invert(aa,ierr)
        if(ierr.ne.0) then
          write(*,*) 'Singular matrix...run terminated'
          pause
          stop
        else
          pp1(n)=aa(1,1)*pp1n+aa(1,2)*pp2n+aa(1,3)*pp3n
          pp2(n)=aa(2,1)*pp1n+aa(2,2)*pp2n+aa(2,3)*pp3n
          pp3(n)=aa(3,1)*pp1n+aa(3,2)*pp2n+aa(3,3)*pp3n
          do 20 ii=1,3
            do 20 jj=1,3
              qm(ii,jj,n)=-ha*aa(ii,jj)
 20       continue
        endif
        epsm=epsp
 30   continue
      uvp(iedge)=0.
      uup(iedge)=0.
      vvp(iedge)=0.
      do 40 n=iedgem,2,-1
        np=n+1
        uvp(n)=qm(1,1,n)*uvp(np)+qm(1,2,n)*uup(np)+qm(1,3,n)*vvp(np)
     *        +pp1(n)
        uup(n)=qm(2,1,n)*uvp(np)+qm(2,2,n)*uup(np)+qm(2,3,n)*vvp(np)
     *        +pp2(n)
        vvp(n)=qm(3,1,n)*uvp(np)+qm(3,2,n)*uup(np)+qm(3,3,n)*vvp(np)
     *        +pp3(n)
 40   continue
c=======================================
c   Compute the deluv and wwp arrays
c=======================================
      do 50 n=2,iedgem
        itrosv=itro
        itro=1
        call propm(n)
        itro=itrosv
        deluv(n)=.5*(coef3*rato1j*(uvp(n)+uvp(n-1))/tm1p
     *          -rato2j*xlm1*(ep(n)+ep(n-1))*(fp(n)-fp(n-1))/dy(n-1))
        wwp(n)=-uup(n)-vvp(n)
 50   continue
      return
      end
      subroutine solvv
c=======================================================================
c       T R A P E Z O I D A L   R U L E   Q U A D R A T U R E :
c       C O N T I N U I T Y   E Q U A T I O N   S O L U T I O N
c=======================================================================
      include 'comeb.'
c=======================================================================
      bo=0.
      vp(1)=vwall
      do 10 n=2,iedge
        b=fp(n)+x*(z1*fp(n)-z2*fo(n)+z3*fn(n))/dx2
        vp(n)=vp(n-1)-.5*(b+bo)*dy(n-1)
        bo=b
 10   continue
      return
      end
      subroutine stable(prod,diss,psifix,ii)
c=======================================================================
c         C O M P U T E   D I A G O N A L   D O M I N A N C E
c                  W E I G H T I N G   F A C T O R
c=======================================================================
c  This holds for either model equation
c=======================================
      if(prod.gt.diss) then
        psifix=prod+.3*diss
      else
        psifix=1.3*diss
      endif
c=======================================
c     Can't make it too big in the
c           omega equation
c=======================================
      if(ii.ne.0) then
        psistp=.5*prod-.25*diss
        if(psistp.gt.psifix) psifix=psistp
      endif
      return
      end
      subroutine strtme
c=======================================================================
c       R E A D   S T A R T I N G   D A T A   F R O M   D I S K
c=======================================================================
      include 'comeb.'
c=======================================================================
c    Read profiles from 'input.dat'
c=======================================
      read(iunit4) mstart,ibound,iedge
      do 10 n=1,iedge
        read(iunit4) dy(n),etn(n),eto(n),fn(n),fo(n),tn(n),to(n),
     *               vn(n),vo(n),wtn(n),wto(n),raton(n),ratoo(n),
     *               wwn(n),wwo(n),uun(n),uuo(n),uvn(n),uvo(n),
     *               vvn(n),vvo(n),deluv(n)
 10   continue
      read(iunit4) dhelp,done,dtwo,d1help,eten,eteo,itstop,
     *             restrt,rstrt,rstrtm,si,si1,testrt,ustrt,
     *             ustrtm,wten,wteo,xistrt
      write(iunit1,200) mstart
c=======================================
c        Initialize arrays
c=======================================
      do 20 n=1,iedge
        fp(n)=fo(n)
        fpp(n)=fo(n)
        ratop(n)=ratoo(n)
        tp(n)=to(n)
        vp(n)=vo(n)
        etp(n)=eto(n)
        wtp(n)=wto(n)
        etpp(n)=eto(n)
        wtpp(n)=wto(n)
        uup(n)=uuo(n)
        uvp(n)=uvo(n)
        vvp(n)=vvo(n)
        wwp(n)=wwo(n)
 20   continue
      iedgep=iedge+1
      if(iedgep.le.jk) then
        do 30 i=iedgep,jk
c-------------------
c     Mean flow
c-------------------
          fn(i)=1.
          fo(i)=1.
          fp(i)=1.
          fpp(i)=1.
          raton(i)=1.
          ratoo(i)=1.
          ratop(i)=1.
          tn(i)=0.
          to(i)=0.
          tp(i)=0.
c-------------------
c k-omega/Stress-omega
c-------------------
          etn(i)=eten
          wtn(i)=wten
          eto(i)=eteo
          wto(i)=wteo
          etpp(i)=eteo
          wtpp(i)=wteo
c-------------------
c Stress-omega model
c-------------------
          uun(i)=0.
          uuo(i)=0.
          uup(i)=0.
          uvn(i)=0.
          uvo(i)=0.
          uvp(i)=0.
          vvn(i)=0.
          vvo(i)=0.
          vvp(i)=0.
          wwn(i)=0.
          wwo(i)=0.
          wwp(i)=0.
          deluv(i)=0.
 30     continue
      endif
c-------------------
c  k-epsilon model
c-------------------
      if(model.eq.2) then
        do 40 n=1,jk
          epn(n)=betas*etn(n)*wtn(n)
          epo(n)=betas*eto(n)*wto(n)
          epp(n)=epn(n)
          eppp(n)=epn(n)
 40     continue
c-------------------
c Baldwin-Barth and
c Spalart-Allmaras
c-------------------
      elseif(model.eq.3.or.model.eq.4) then
        wtbar=1./betas
        do 45 n=1,jk
          etn(n)=etn(n)/(betas*wtn(n))
          eto(n)=eto(n)/(betas*wto(n))
          etp(n)=etn(n)
          etpp(n)=etn(n)
          wtn(n)=wtbar
          wto(n)=wtbar
          wtp(n)=wtbar
          wtpp(n)=wtbar
 45     continue
        eten=eten/(betas*wten)
        eteo=eteo/(betas*wteo)
        wten=wtbar
        wteo=wtbar
      endif
c=======================================
c  Compute dissipation for models that
c  solve for turbulence energy....it's
c  needed for the mean energy equation
c=======================================
      do 50 n=1,jk
        if(model.ge.0.and.model.le.2) then
          dissk(n)=betas*etp(n)*wtp(n)
        else
          dissk(n)=0.
        endif
 50   continue
c=======================================
c    Determine maximum shear stress 
c      for the Johnson-King model
c=======================================
      if(model.eq.7) then
        uvmaxn=0.
        do 60 n=2,jk
          if(uvn(n).gt.uvmaxn) then
            uvmaxn=uvn(n)
          else
            go to 70
          endif
 60     continue
 70     uvmaxo=0.
        do 80 n=2,jk
          if(uvo(n).gt.uvmaxo) then
            uvmaxo=uvo(n)
          else
            go to 90
          endif
 80     continue
 90     nmaxjk=n-1
        uvmaxp=uvmaxo
        uvmxpp=uvmaxp
      endif
c=======================================
c   Read edge and surface conditions
c=======================================
      do 100 i=1,number
        read(iunit2) as(i),ape(i),apep(i),aqw(i),aqwp(i),
     *               arhov(i),arhovp(i),armi(i),armip(i),atw(i),
     *               atwp(i),az(i),azp(i),azcrv(i),azcrvp(i)
 100  continue
      close(iunit2)
      close(iunit4)
      return
c=======================================================================
 200  format(2x,'COMPUTATION INITIATED AT MSTART =',i4)
c=======================================================================
      end
      subroutine update
c=======================================================================
c       C O M P U T E   V I S C O S I T Y ,   V E L O C I T Y
c     A N D   T E M P E R A T U R E   G R A D I E N T ,   E T C .
c=======================================================================
      include 'comeb.'
c=======================================================================
c   Calculate wall and initial values
c   required for basic B.L. parameters
c=======================================
      xlm11(1)=(1.+tr)*(1.+tp(1))**vispom/(1.+tp(1)+tr)
      xlm1p=xlm11(1)
      if(kodwal.eq.1)
     *tz(1)=(-ww1*tp(1)+ww2*tp(2)-ww3*tp(3)+ww4*tp(4))/(ww5*dy(1))
      fz(1)=(-ww1*fp(1)+ww2*fp(2)-ww3*fp(3)+ww4*fp(4))/(ww5*dy(1))
      xmf4=eps*xnue*rmi1j*xlm1p*fz(1)/oz
      tpp=1.+tp(1)
      uplusb=ue*sqrt(xmf4*tpp)
      tscale=re*uplusb/(eps*xnue*xlm1p*tpp*tpp)
      tpk=0.
      bex1=sgn*rmi1j/(eps*costh)
      bex3=2.*sgn*eps*costh*oz/(re*ue*rmi2j)
      nuflag=0
      nskip=2
      xn(1)=0.
      ratop(1)=1.
      do 10 n=2,iedge
c=======================================
c  Calculate distance normal to surface
c       and transverse curvature
c=======================================
        tpk=tpk+.5*(2.+tp(n-1)+tp(n))*dy(n-1)
        xn(n)=faa*tpk
        if(j.ne.0) then
          arg=bex3*tpk
          if(arg.lt.-1.) go to 20
          ratop(n)=sqrt(1.+arg)
          xn(n)=bex1*(ratop(n)-1.)
        endif
c=======================================
c    Determine point to which omega
c      is prescribed analytically
c=======================================
        if(nuflag.ne.1) then
          yplus(n)=tscale*xn(n)
          if(yplus(n).le.ustop) then
            nskip=n
          else
            nuflag=1
          endif
        endif
c=======================================
c  Velocity and temperature derivatives
c=======================================
        if(n.lt.iedge) then
          fz(n)=y4(n)*fp(n+1)-y5(n)*fp(n)-y6(n)*fp(n-1)
          tz(n)=y4(n)*tp(n+1)-y5(n)*tp(n)-y6(n)*tp(n-1)
        endif
c=======================================
c         Molecular viscosity
c=======================================
        xlm11(n)=(1.+tr)*(1.+tp(n))**vispom/(1.+tp(n)+tr)
 10   continue
c=======================================
c    Compute edge value of coordinate
c    normal to surface and velocity &
c         temperature gradients
c=======================================
 20   fz(iedge)=(1.-fp(iedge-1))/(dy(iedge-1)+dy(iedge))
      tz(iedge)=-tp(iedge-1)/(dy(iedge-1)+dy(iedge))
      xnden=(xn(netaeg-2)*(1.-fp(netaeg-1))
     *      -xn(netaeg-1)*(1.-fp(netaeg-2)))
     *         /(fp(netaeg-2)-fp(netaeg-1))
      return
      end
      subroutine wall(m)
c=======================================================================
c         P R I N T   I N T E G R A L   P A R A M E T E R S ,
c             E D G E   C O N D I T I O N S ,   E T C .
c=======================================================================
      include 'comeb.'
      character*7 label
c=======================================================================
      negtiv=jprint+kprint
      pestar=pe*r1*u1*u1
      testar=te*tref
      restar=r1*re
      uestar=u1*ue
      upstar=u1*uplusb
      amustr=xnue*visref
      yestar=xnden*eps
      dueds=-pp/(re*ue)
      dteds=-ue*dueds
      if(model.eq.7) then
        dsmxo=sigmaj
        label='sigma ='
      elseif(model.ge.0.and.model.le.2) then
        dsmxo=2.*beta50*emax/cfe
        label='kmax  ='
      else
        dsmxo=epsmax/redelv
        label='epmax ='
      endif
      write(iunit1,100) res,s,uestar,dueds,qsd,m,cfe,taud,pestar,
     *                  pp,hd,netaeg,cfw,theta,testar,dteds,che,
     *                  iedge,rethet,disp,restar,xbe,chw,itro,
     *                  redelt,yestar,xmae,x,anue,nerror,thadis,
     *                  rvwald,amustr,rmi,anuw,negtiv,label,dsmxo,
     *                  upstar,yplus(2),z,recov,nskip
      return
c=======================================================================
 100  format(
     *      /2x,'Res   =',1pe12.5,2x,'s     =',e12.5,2x,'Ue    =',e12.5,
     *       2x,'dUe/ds=',  e12.5,2x,'qw    =',e12.5,2x,'M     =',i5   ,
     *      /2x,'Cfe   =',  e12.5,2x,'tauw  =',e12.5,2x,'Pe    =',e12.5,
     *       2x,'dPe/ds=',  e12.5,2x,'hdot  =',e12.5,2x,'Ne    =',i5   ,
     *      /2x,'Cfw   =',  e12.5,2x,'theta =',e12.5,2x,'Te    =',e12.5,
     *       2x,'dTe/ds=',  e12.5,2x,'Nste  =',e12.5,2x,'Iedge =',i5   ,
     *      /2x,'Rethet=',  e12.5,2x,'delta*=',e12.5,2x,'Rhoe  =',e12.5,
     *       2x,'beta  =',  e12.5,2x,'Nstw  =',e12.5,2x,'Itro  =',i5   ,
     *      /2x,'Redel*=',  e12.5,2x,'delta =',e12.5,2x,'Me    =',e12.5,
     *       2x,'xi    =',  e12.5,2x,'Nue   =',e12.5,2x,'Nerror=',i5   ,
     *      /2x,'H     =',  e12.5,2x,'rho*vw=',e12.5,2x,'Mue   =',e12.5,
     *       2x,'radius=',  e12.5,2x,'Nuw   =',e12.5,2x,'Negtiv=',i5   ,
     *      /2x,   a7    ,  e12.5,2x,'utau  =',e12.5,2x,'yplus =',e12.5,
     *       2x,'z     =',  e12.5,2x,'Recov =',e12.5,2x,'Nskip =',i5)
      end
      subroutine wbcke(etp1,epp1)
c=======================================================================
c      C O M P U T E   Q U A N T I T I E S   N E E D E D   T O
c      S P E C I F Y   T U R B U L E N T   E N E R G Y   A N D
c  D I S S I P A T I O N   A T   S U R F A C E : K - E P S I L O N
c-----------------------------------------------------------------------
c  NOTE:  Using conventional difference formulas for d2kdy2 and dkhdy
c         causes numerical problems with the k-epsilon model.
c=======================================================================
      include 'comeb.'
c=======================================================================
      etp1=0.
      if(kebc.eq.0) then
        epp1=0.
      elseif(kebc.eq.1) then
        d2kdy2=2.*etp(3)/((1.+xk)*dy(1))**2
        epp1=xlm11(1)*d2kdy2/coef9
      elseif(kebc.eq.2) then
        dkhdy=sqrt(etp(3))/((1.+xk)*dy(1))
        epp1=2.*xlm11(1)*dkhdy**2/coef9
      elseif(kebc.eq.3) then
        epp1=eppp(2)
      endif
      return
      end
      subroutine wbckw(etp1,wtp1)
c=======================================================================
c      C O M P U T E   Q U A N T I T I E S   N E E D E D   T O
c      S P E C I F Y   T U R B U L E N T   E N E R G Y   A N D
c    D I S S I P A T I O N   A T   S U R F A C E : K - O M E G A
c=======================================================================
      include 'comeb.'
c=======================================================================
c        Roughness contribution
c=======================================
      xnuwal=viscon*(tw*tref)**vispow/(su+tw*tref)
      rekedg=rek*re*ue*vis1/(visref*xnue)
      rekwal=rekedg*te*xnue*visref/(tw*xnuwal)
      utkonu=rekwal*uplusb/ue
      srin=(utkonu/200.)**2
      if(utkonu.gt.5.) then
        if(model.eq.0) then
          if(nvisc.eq.0) then
            srinf=100./utkonu
          else
            srinf=60./utkonu
          endif
        elseif(model.eq.1) then
          if(nvisc.eq.0) then
            srinf=70./utkonu
          else
            srinf=50./utkonu
          endif
        endif
        srin=1./(srinf+(1./srin-srinf)*exp(5.-utkonu))
      endif
c=======================================
c        Blowing contribution
c=======================================
      vwouta=rvwald*tw/(te*uplusb*re*r1*u1)
      if(vwouta.lt.0.) vwouta=0.
      if(model.eq.0) then
        if(nvisc.eq.0) then
          sbin=(1.+5.*vwouta)*vwouta/25.
        else
          sbin=(1.+5.*vwouta)*vwouta/14.
        endif
      elseif(model.eq.1) then
        if(nvisc.eq.0) then
          sbin=(1.+5.*vwouta)*vwouta/24.
        else
          sbin=(1.+5.*vwouta)*vwouta/11.
        endif
      endif
c=======================================
c       Combined contribution
c=======================================
      esinv=srin+sbin
      epsr=sqrt(esinv)
      wwall=oy*uplusb*uplusb*u1*r1*re*te/(ue*ue*xnuwal*tw)
      yscale=r1*u1*eps*re*te*uplusb/(xnuwal*tw)
c=======================================
c     Surface values of etp and wtp
c=======================================
      etp1=0.
      if(esinv.lt.1.e-10) then
        es=1.e15
        wtp1=1.e15
      else
        es=1./esinv
        wtp1=wwall/esinv
      endif
      return
      end
      subroutine width(iedgem,ierr)
c=======================================================================
c    C H E C K   F O R   S U F F I C I E N T   G R I D   W I D T H
c=======================================================================
      include 'comeb.'
c=======================================================================
      ierr=0
      delfp=abs(fp(iedgem)-fp(iedge-2))
      deltp=abs(tp(iedgem)-tp(iedge-2))
      if(delfp.gt..0001.or.deltp.gt..0001) then
        if(iedgem.gt.jk) then
          ierr=1
        else
          iedgem=iedge
          iedge=iedge+1
          call adline
        endif
      endif
      return
      end
      subroutine xgrid(m)
c=======================================================================
c C O M P U T E   T H E   D I F F E R E N C E   C O E F F I C I E N T S
c          F O R   S T R E A M W I S E   D E R I V A T I V E S
c=======================================================================
      include 'comeb.'
c=======================================================================
c      Extrapolation coefficients
c=======================================
      z4=(dx1+dx2)/dx1
      z5=dx2/dx1
c=======================================
c Zero-gradient difference coefficients
c=======================================
      if(m.eq.2) then
        z1=1.
        z2=1.
        z3=0.
c=======================================
c  General-case difference coefficients
c=======================================
      else
        z1=2.*((dx1+2.*dx2)/(dx1+dx2))
        z2=2.*(dx1+dx2)/dx1
        z3=2.*((dx2*dx2)/(dx1*(dx1+dx2)))
      endif
      return
      end
      subroutine ygrid
c=======================================================================
c                 C O M P U T E   N O R M A L - G R I D
c             D I F F E R E N C E   C O E F F I C I E N T S
c=======================================================================
      include 'comeb.'
c=======================================================================
      jkm=jk-1
      dy1=dy(1)
      do 10 n=2,jkm
        dy2=xk*dy1
        y2(n)=1./(dy1+dy2)
        y1(n)=2.*y2(n)/dy2
        y3(n)=2.*y2(n)/dy1
        y4(n)=dy1/(dy2*(dy1+dy2))
        y5(n)=(dy1-dy2)/(dy1*dy2)
        y6(n)=dy2/(dy1*(dy1+dy2))
        dy1=dy2
 10   continue
      w1=xk
      w2=1.+w1
      w3=1.+w1+w1*w1
      ww1=w3*w3*(w1*w2-1.)+w2
      ww2=w1*w2*w3*w3
      ww3=w3*w3
      ww4=w2
      ww5=w1*w1*w1*w2*w3
      return
      end
