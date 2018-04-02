      program eddy2c
c=======================================================================
c        TWO-DIMENSIONAL NAVIER-STOKES PROGRAM...COMPRESSIBLE FLOW
c                  MacCormack 1985 Implicit Algorithm
c-----------------------------------------------------------------------
c                   MODEL =-1...........Laminar Flow
c                   MODEL = 0..........k-omega Model
c                   MODEL = 1.....Stress-omega Model
c                   MODEL = 2....Baldwin-Lomax Model
c                   --------------------------------
c                   NVISC = 0...Exclude Viscous Mods
c                   NVISC = 1...Include Viscous Mods
c-----------------------------------------------------------------------
c    Copyright (c) 2006 DCW Industries, Inc.  All rights reserved.
c=======================================================================
      include 'come2c.'
      include 'cpuid.'
c==============================
c  Lahey system-time utility
c==============================
      call cpu_time(start_time)
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
      gmoprl=gam/prl
      gmoprt=gam/prt
      cv=rgas/gm1
      betam=gm1
c=======================================
c           Flow Properties
c=======================================
      if(iadbwl.eq.1) twall=tinf*(1.+.5*gm1*prt*xma**2)
      tt1=tinf*(1.+.5*gm1*xma**2)
      pt1=pinf*(1.+.5*gm1*xma**2)**(gam/gm1)
      rt1=pt1/(rgas*tt1)
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
      dynp=.5*rhoinf*uinf**2
      qdm=gam*((1.+.5*gm1*xma**2)*eiinf-eiwall)*rhouin
c=======================================
c         Counters and Flags
c=======================================
      if(ile.lt.1) ile=1
      ilm1=il-1
      jlm1=jl-1
      iskip=il/16
      if(16*iskip.ne.il) iskip=iskip+1
c=======================================
c     Turbulence-Model Parameters
c=======================================
      call vcoefs
      epinf=retinf*rmuinf
      rhokin=rhoinf*tkeinf
      rhomin=rhoinf*rhokin/epinf
      rhomw=40000./rough**2
      rhokmn=1.e-9*rhokin
      rhowmn=1.e-2*rhomin
      xmtsq=xmt0**2
c=======================================
c      Miscellaneous Parameters
c=======================================
      resid=0.
      faxi=.1
      implt=1
      n=1
      m=1
c=======================================
c     Initial Flowfield and Mesh
c=======================================
      call readit
c=======================================
c    More Miscellaneous Parameters
c=======================================
      nstart=nadv+1
      call bc
      if(model.eq.0.or.model.eq.1) then
        call epskw
        call vgrads
        call crdiff
        call strech
        if(model.eq.1.and.ntau.ne.0) call dtau
      endif
      k1=1
      k2=2
      imesh=0
c=======================================================================
c             * * * * *   M A I N   L O O P   * * * * *
c=======================================================================
      do 10 nadv=nstart,nend
        call tmestp(dt)
        if(model.eq.2) call turbfl
        call l(dt)
        write(*,100) nadv,resid
c====================================
c  Realizability violation checking
c====================================
c        if(nadv.eq.1) then
c          write(*,100) nadv,resid
c        else
c          realiz=100.*float(il*jl-irealz)/float(il*jl)
c          write(*,110) nadv,resid,realiz
c        endif
c====================================
        if(mod(nadv,nprnt).eq.0.or.nadv.eq.nend) call prtwal
 10   continue
 20   nend=nadv-1
      call saveit
c=======================================
c  Print profiles and save plotting data
c=======================================
      if(ipro.gt.0) call prntff
      close(iunit1)
      call plotme
      close(iuplot)
c=======================================
c    Pause to keep window visible after
c    the run completes...
c=======================================
c  Lahey system-time utility
c==============================
      call cpu_time(end_time)
      cputime=end_time-start_time
      write(*,*) 'CPU Time = ',cputime
      pause
c=======================================================================
 100  format(1x,'Nadv = ',i6,3x,'Residual = ',1pe9.3)
c====================================
c  Realizability violation checking
c====================================
c 110  format(1x,'Nadv = ',i6,3x,'Residual = ',1pe9.3,3x,
c     *          'Realizability Index = ',0pf5.1,'%')
c=======================================================================
      end
      subroutine bc
c=======================================================================
c                B O U N D A R Y   C O N D I T I O N S
c=======================================================================
      include 'come2c.'
c=======================================================================
      pe=pinf+dynp*pexit
c==========================================
c         Lower surface (j=1)
c==========================================
      do 30 i=2,ilm1
        v(i,1)=-v(i,2)
        p(i,1)=p(i,2)
        rho(i,1,m)=rho(i,2,m)
        e(i,1,m)=e(i,2,m)
c===================
c  Symmetry plane
c===================
        if(i.lt.ile.or.i.gt.ite) then
          u(i,1)=u(i,2)
          ei(i,1)=ei(i,2)
          rhou(i,1,m)=rhou(i,2,m)
          rhov(i,1,m)=-rhov(i,2,m)
          if(model.eq.0.or.model.eq.1) then
            tke(i,1)=tke(i,2)
            omeg(i,1)=omeg(i,2)
            if(model.eq.1) then
              rhotxx(i,1,m)=rhotxx(i,2,m)
              rhotxy(i,1,m)=-rhotxy(i,2,m)
              rhotyy(i,1,m)=-rhotyy(i,2,m)+1.333333*rho(i,2,m)*tke(i,2)
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
            endif
          endif
        endif
        if(model.eq.0.or.model.eq.1) then
          rhok(i,1,m)=rho(i,1,m)*tke(i,1)
          rhomeg(i,1,m)=rho(i,1,m)*omeg(i,1)
        endif
 30   continue
c==========================================
c         Symmetry plane (i=1)
c==========================================
      do 40 j=1,jl
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
              rhotxy(1,j,m)=-rhotxy(2,j,m)
              rhotyy(1,j,m)=-rhotyy(2,j,m)+1.333333*rhok(2,j,m)
            endif
          endif
        endif
c==========================================
c                Exit plane
c==========================================
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
          tke(il,j)=rhok(il,j,m)/rho(il,j,m)
          omeg(il,j)=rhomeg(il,j,m)/rho(il,j,m)
          if(model.eq.1) then
            rhotxx(il,j,m)=rhotxx(ilm1,j,m)
            rhotxy(il,j,m)=rhotxy(ilm1,j,m)
            rhotyy(il,j,m)=rhotyy(ilm1,j,m)
          endif
        endif
        sx=y(il,j+1)-y(il,j)
        sy=-(x(il,j+1)-x(il,j))
        si=sqrt(sx*sx+sy*sy)
        sxp=sx/si
        syp=sy/si
        up=sxp*u(ilm1,j)+syp*v(ilm1,j)
        cc=sqrt(ggm1*abs(ei(il,j)))
c==========================================
c Set exit pressure if the flow is subsonic
c     iexit = 0...Supersonic
c     Otherwise...Subsonic
c==========================================
        if(iexit.ne.0) then
          if(up.lt.cc) then
            dr=(pe-p(il,j))/(cc*cc)
            p(il,j)=pe
            rho(il,j,m)=rho(il,j,m)+dr
            rhou(il,j,m)=rho(il,j,m)*u(il,j)
            rhov(il,j,m)=rho(il,j,m)*v(il,j)
            ei(il,j)=p(il,j)/(gm1*rho(il,j,m))
            e(il,j,m)=rho(il,j,m)*(ei(il,j)+.5*(u(il,j)**2+v(il,j)**2)
     *               +tke(il,j))
          endif
c======================
c Prevent inflow across
c   the exit boundary
c======================
          if(j.ne.1) then
            up=abs(up)
            vp=-syp*u(il,j)+sxp*v(il,j)
            uu=sxp*up-syp*vp
            vv=syp*up+sxp*vp
            rhou(il,j,m)=rho(il,j,m)*uu
            rhov(il,j,m)=rho(il,j,m)*vv
            u(il,j)=uu
            v(il,j)=vv
            p(il,j)=p(il,j)-rho(il,j,m)*cc*(sxp*(uu-u(ilm1,j))
     *                                     +syp*(vv-v(ilm1,j)))
            ei(il,j)=p(il,j)/(gm1*rho(il,j,m))
            e(il,j,m)=rho(il,j,m)*(ei(il,j)+.5*(u(il,j)**2+v(il,j)**2)
     *               +tke(il,j))
          endif
        endif
 40   continue
c==========================================
c   Symmetry plane/Freestream (j = jl)
c==========================================
      if(jsym.ne.0) then
        do 50 i=1,il
          rho(i,jl,m)=rho(i,jlm1,m)
          rhou(i,jl,m)=rhou(i,jlm1,m)
          u(i,jl)=u(i,jlm1)
          ei(i,jl)=ei(i,jlm1)
          p(i,jl)=p(i,jlm1)
          v(i,jl)=2.*(y(i+1,jl)-y(i,jl))/(x(i+1,jl)-x(i,jl))*u(i,jlm1)
     *           -v(i,jlm1)
          rhov(i,jl,m)=rho(i,jl,m)*v(i,jl)
          e(i,jl,m)=rho(i,jl,m)*(ei(i,jl)+.5*(u(i,jl)**2+v(i,jl)**2)
     *             +tke(i,jl))
 50     continue
        if(model.eq.0.or.model.eq.1) then
          do 60 i=1,il
            rhok(i,jl,m)=rhok(i,jlm1,m)
            rhomeg(i,jl,m)=rhomeg(i,jlm1,m)
            tke(i,jl)=tke(i,jlm1)
            omeg(i,jl)=omeg(i,jlm1)
 60       continue
        endif
        if(model.eq.1) then
          do 70 i=1,il
            rhotxx(i,jl,m)=rhotxx(i,jlm1,m)
            rhotxy(i,jl,m)=rhotxy(i,jlm1,m)
            rhotyy(i,jl,m)=rhotyy(i,jlm1,m)
 70       continue
        endif
      endif
      return
      end
      block data blkdat
c=======================================================================
c             A R R A Y   I N I T I A L I Z A T I O N
c=======================================================================
      include 'come2c.'
      parameter (jmx21=21*jmx)
c=======================================================================
      data     fs1,fs2,fs3,fs4,fs5,fs6,fs7,fs8,fs9/18*0./
      data     g11,g12,g13,g14,g21,g22,g23,g24,g31,g32,g33,g34,
     *         g41,g42,g43,g44,g55,g66,g77,g88,g99/jmx21*0./
      data     a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,
     *         a41,a42,a43,a44,a55,a66,a77,a88,a99/21*0./
      data     b11,b12,b13,b14,b21,b22,b23,b24,b31,b32,b33,b34,
     *         b41,b42,b43,b44,b55,b66,b77,b88,b99/21*0./
      data     c11,c12,c13,c14,c21,c22,c23,c24,c31,c32,c33,c34,
     *         c41,c42,c43,c44,c55,c66,c77,c88,c99/21*0./
      end
      subroutine btdj(dt,i,iadd,jadd)
c=======================================================================
c           B L O C K   T R I D I A G O N A L   M A T R I X
c  I N V E R S I O N   P R O C E D U R E   I N   J - D I R E C T I O N
c=======================================================================
      include 'come2c.'
c=======================================================================
      g11(jl)=0.
      g12(jl)=0.
      g13(jl)=0.
      g14(jl)=0.
      g21(jl)=0.
      g22(jl)=0.
      g23(jl)=0.
      g24(jl)=0.
      g31(jl)=0.
      g32(jl)=0.
      g33(jl)=0.
      g34(jl)=0.
      g41(jl)=0.
      g42(jl)=0.
      g43(jl)=0.
      g44(jl)=0.
      if(model.eq.0.or.model.eq.1) then
        g55(jl)=0.
        g66(jl)=0.
        if(model.eq.1) then
          g77(jl)=0.
          g88(jl)=0.
          g99(jl)=0.
        endif
      endif
c=======================================
c Sweep from top (j=jl) to bottom (j=2)
c=======================================
      call implj(dt,i,jl,jadd)
      do 10 jx=2,jlm1
        j=jlm1+2-jx
        vr=vol(i,j+1)/vol(i,j)
c=======================================
c     {Du}j = {Du}j - [b]j {D}j+1
c=======================================
        du1=du(i,j+1,1)
        du2=du(i,j+1,2)
        du3=du(i,j+1,3)
        du4=du(i,j+1,4)
        du(i,j,1)=du(i,j,1)-vr*(b11*du1+b12*du2+b13*du3+b14*du4)
        du(i,j,2)=du(i,j,2)-vr*(b21*du1+b22*du2+b23*du3+b24*du4)
        du(i,j,3)=du(i,j,3)-vr*(b31*du1+b32*du2+b33*du3+b34*du4)
        du(i,j,4)=du(i,j,4)-vr*(b41*du1+b42*du2+b43*du3+b44*du4)
        if(model.eq.0.or.model.eq.1) then
          du(i,j,5)=du(i,j,5)-vr*b55*du(i,j+1,5)
          du(i,j,6)=du(i,j,6)-vr*b66*du(i,j+1,6)
          if(model.eq.1) then
            du(i,j,7)=du(i,j,7)-vr*b77*du(i,j+1,7)
            du(i,j,8)=du(i,j,8)-vr*b88*du(i,j+1,8)
            du(i,j,9)=du(i,j,9)-vr*b99*du(i,j+1,9)
          endif
        endif
c=======================================
c        [A]j = [b]j [g]j+1
c=======================================
        a11=-vr*(c11+b11*g11(j+1)+b12*g21(j+1)
     *              +b13*g31(j+1)+b14*g41(j+1))+1.
        a12=-vr*(c12+b11*g12(j+1)+b12*g22(j+1)
     *              +b13*g32(j+1)+b14*g42(j+1))
        a13=-vr*(c13+b11*g13(j+1)+b12*g23(j+1)
     *              +b13*g33(j+1)+b14*g43(j+1))
        a14=-vr*(c14+b11*g14(j+1)+b12*g24(j+1)
     *              +b13*g34(j+1)+b14*g44(j+1))
        a21=-vr*(c21+b21*g11(j+1)+b22*g21(j+1)
     *              +b23*g31(j+1)+b24*g41(j+1))
        a22=-vr*(c22+b21*g12(j+1)+b22*g22(j+1)
     *              +b23*g32(j+1)+b24*g42(j+1))+1.
        a23=-vr*(c23+b21*g13(j+1)+b22*g23(j+1)
     *              +b23*g33(j+1)+b24*g43(j+1))
        a24=-vr*(c24+b21*g14(j+1)+b22*g24(j+1)
     *              +b23*g34(j+1)+b24*g44(j+1))
        a31=-vr*(c31+b31*g11(j+1)+b32*g21(j+1)
     *              +b33*g31(j+1)+b34*g41(j+1))
        a32=-vr*(c32+b31*g12(j+1)+b32*g22(j+1)
     *              +b33*g32(j+1)+b34*g42(j+1))
        a33=-vr*(c33+b31*g13(j+1)+b32*g23(j+1)
     *              +b33*g33(j+1)+b34*g43(j+1))+1.
        a34=-vr*(c34+b31*g14(j+1)+b32*g24(j+1)
     *              +b33*g34(j+1)+b34*g44(j+1))
        a41=-vr*(c41+b41*g11(j+1)+b42*g21(j+1)
     *              +b43*g31(j+1)+b44*g41(j+1))
        a42=-vr*(c42+b41*g12(j+1)+b42*g22(j+1)
     *              +b43*g32(j+1)+b44*g42(j+1))
        a43=-vr*(c43+b41*g13(j+1)+b42*g23(j+1)
     *              +b43*g33(j+1)+b44*g43(j+1))
        a44=-vr*(c44+b41*g14(j+1)+b42*g24(j+1)
     *              +b43*g34(j+1)+b44*g44(j+1))+1.
        if(model.eq.0.or.model.eq.1) then
          a55=-vr*(c55+b55*g55(j+1))+1.+dt*(pk(i,j)
     *        +2.*betas*fbet(i,j)*omeg(i,j))
          a66=-vr*(c66+b66*g66(j+1))+1.+dt*(alpha*pk(i,j)
     *        +2.*beta*fbeta(i,j)*omeg(i,j))
          if(model.eq.1) then
            a77=-vr*(c77+b77*g77(j+1))+1.+dt*ptij(i,j)
            a88=-vr*(c88+b88*g88(j+1))+1.+dt*ptij(i,j)
            a99=-vr*(c99+b99*g99(j+1))+1.+dt*ptij(i,j)
          endif
        endif
c=======================================
c    Compute j-flux contributions
c=======================================
        call implj(dt,i,j,jadd)
c=======================================
c    Compute i-flux contributions
c=======================================
        call impli(dt,i,j,iadd)
        r11=a11
        r21=a21
        r31=a31
        r41=a41
        s12=a12/r11
        s13=a13/r11
        s14=a14/r11
        r22=a22-r21*s12
        r32=a32-r31*s12
        r42=a42-r41*s12
        s23=(a23-r21*s13)/r22
        s24=(a24-r21*s14)/r22
        r33=a33-r31*s13-r32*s23
        r43=a43-r41*s13-r42*s23
        s34=(a34-r31*s14-r32*s24)/r33
        r44=a44-r41*s14-r42*s24-r43*s34
        yy1=du(i,j,1)/r11
        yy2=(du(i,j,2)-r21*yy1)/r22
        yy3=(du(i,j,3)-r31*yy1-r32*yy2)/r33
        yy4=(du(i,j,4)-r41*yy1-r42*yy2-r43*yy3)/r44
        yy3=yy3-s34*yy4
        yy2=yy2-s23*yy3-s24*yy4
        yy1=yy1-s12*yy2-s13*yy3-s14*yy4
        du(i,j,1)=yy1
        du(i,j,2)=yy2
        du(i,j,3)=yy3
        du(i,j,4)=yy4
        if(model.eq.0.or.model.eq.1) then
          du(i,j,5)=du(i,j,5)/a55
          du(i,j,6)=du(i,j,6)/a66
          if(model.eq.1) then
            du(i,j,7)=du(i,j,7)/a77
            du(i,j,8)=du(i,j,8)/a88
            du(i,j,9)=du(i,j,9)/a99
          endif
        endif
c=======================================
c             Compute [g]
c=======================================
        if(j.ne.2) then
          yy1=c11/r11
          yy2=(c21-r21*yy1)/r22
          yy3=(c31-r31*yy1-r32*yy2)/r33
          g41(j)=(c41-r41*yy1-r42*yy2-r43*yy3)/r44
          g31(j)=yy3-s34*g41(j)
          g21(j)=yy2-s23*g31(j)-s24*g41(j)
          g11(j)=yy1-s12*g21(j)-s13*g31(j)-s14*g41(j)
          yy1=c12/r11
          yy2=(c22-r21*yy1)/r22
          yy3=(c32-r31*yy1-r32*yy2)/r33
          g42(j)=(c42-r41*yy1-r42*yy2-r43*yy3)/r44
          g32(j)=yy3-s34*g42(j)
          g22(j)=yy2-s23*g32(j)-s24*g42(j)
          g12(j)=yy1-s12*g22(j)-s13*g32(j)-s14*g42(j)
          yy1=c13/r11
          yy2=(c23-r21*yy1)/r22
          yy3=(c33-r31*yy1-r32*yy2)/r33
          g43(j)=(c43-r41*yy1-r42*yy2-r43*yy3)/r44
          g33(j)=yy3-s34*g43(j)
          g23(j)=yy2-s23*g33(j)-s24*g43(j)
          g13(j)=yy1-s12*g23(j)-s13*g33(j)-s14*g43(j)
          yy1=c14/r11
          yy2=(c24-r21*yy1)/r22
          yy3=(c34-r31*yy1-r32*yy2)/r33
          g44(j)=(c44-r41*yy1-r42*yy2-r43*yy3)/r44
          g34(j)=yy3-s34*g44(j)
          g24(j)=yy2-s23*g34(j)-s24*g44(j)
          g14(j)=yy1-s12*g24(j)-s13*g34(j)-s14*g44(j)
          if(model.eq.0.or.model.eq.1) then
            g55(j)=c55/a55
            g66(j)=c66/a66
            if(model.eq.1) then
              g77(j)=c77/a77
              g88(j)=c88/a88
              g99(j)=c99/a99
            endif
          endif
        endif
 10   continue
c=======================================
c            Update {du}
c=======================================
      do 20 j=3,jlm1
        du1=du(i,j-1,1)
        du2=du(i,j-1,2)
        du3=du(i,j-1,3)
        du4=du(i,j-1,4)
        du(i,j,1)=du(i,j,1)-g11(j)*du1-g12(j)*du2-g13(j)*du3
     *                     -g14(j)*du4
        du(i,j,2)=du(i,j,2)-g21(j)*du1-g22(j)*du2-g23(j)*du3
     *                     -g24(j)*du4
        du(i,j,3)=du(i,j,3)-g31(j)*du1-g32(j)*du2-g33(j)*du3
     *                     -g34(j)*du4
        du(i,j,4)=du(i,j,4)-g41(j)*du1-g42(j)*du2-g43(j)*du3
     *                     -g44(j)*du4
        if(model.eq.0.or.model.eq.1) then
          du(i,j,5)=du(i,j,5)-g55(j)*du(i,j-1,5)
          du(i,j,6)=du(i,j,6)-g66(j)*du(i,j-1,6)
          if(model.eq.1) then
            du(i,j,7)=du(i,j,7)-g77(j)*du(i,j-1,7)
            du(i,j,8)=du(i,j,8)-g88(j)*du(i,j-1,8)
            du(i,j,9)=du(i,j,9)-g99(j)*du(i,j-1,9)
          endif
        endif
 20   continue
      return
      end
      subroutine crdiff
c=======================================================================
c         C R O S S - D I F F U S I O N   F U N C T I O N
c=======================================================================
      include 'come2c.'
c=======================================================================
      do 20 i=2,ilm1
        do 10 j=2,jlm1
          dke=.5*(tke(i+1,j)-tke(i-1,j))
          dwe=.5*(omeg(i+1,j)-omeg(i-1,j))
          dxe=.5*(x(i+1,j+1)+x(i+1,j)-x(i,j+1)-x(i,j))
          dye=.5*(y(i+1,j+1)+y(i+1,j)-y(i,j+1)-y(i,j))
          dkn=.5*(tke(i,j+1)-tke(i,j-1))
          dwn=.5*(omeg(i,j+1)-omeg(i,j-1))
          dxn=.5*(x(i+1,j+1)+x(i,j+1)-x(i+1,j)-x(i,j))
          dyn=.5*(y(i+1,j+1)+y(i,j+1)-y(i+1,j)-y(i,j))
          dxy=1./(dxe*dyn-dxn*dye)
          dkdx=(dke*dyn-dkn*dye)*dxy
          dwdx=(dwe*dyn-dwn*dye)*dxy
          dkdy=(dkn*dxe-dke*dxn)*dxy
          dwdy=(dwn*dxe-dwe*dxn)*dxy
          chik=(dkdx*dwdx+dkdy*dwdy)/omeg(i,j)**3
          if(chik.le.0.) then
            cross(i,j)=0.
          else
            cross(i,j)=sigmad*chik*omeg(i,j)
          endif
 10     continue
 20   continue
      return
      end
      subroutine dtau
c=======================================================================
c              S T R E S S   D I F F E R E N T I A L S
c=======================================================================
      include 'come2c.'
c=======================================================================
      do 20 j=2,jlm1
        do 10 i=2,ilm1
          rlmbda=-.66666667*ep(i,j)
          divu=dudx(i,j)+dvdy(i,j)
          if(iaxi.eq.1) divu=divu+vor(i,j)
          dtxx(i,j)=rhotxx(i,j,m)-rlmbda*divu-2.*ep(i,j)*dudx(i,j)
          dtxy(i,j)=rhotxy(i,j,m)-ep(i,j)*(dudy(i,j)+dvdx(i,j))
          dtyy(i,j)=rhotyy(i,j,m)-rlmbda*divu-2.*ep(i,j)*dvdy(i,j)
 10     continue
 20   continue
      return
      end
      subroutine epskw
c=======================================================================
c        E D D Y   V I S C O S I T Y / K - O M E G A   M O D E L
c=======================================================================
      include 'come2c.'
c=======================================================================
      do 20 i=1,il
        ep(i,1)=0.
        ep(i,jl)=epinf
        do 10 j=2,jlm1
          ep(i,j)=rhok(i,j,m)/omeg(i,j)
 10     continue
 20   continue
      if(nvisc.ne.0.and.model.eq.0) then
          call viskw
      elseif(nvisc.ne.0.and.model.eq.1) then
          call vissw
      endif
      return
      end
      subroutine fsi(i,j,ii)
c=======================================================================
c    T R A N S P O R T   A N D   S T R E S S...X - D I R E C T I O N
c=======================================================================
      include 'come2c.'
c=======================================================================
c          I N V I S C I D   T R A N S P O R T   T E R M S
c=======================================================================
c   Cell geometry
c-------------------
      sx=(y(i+1,j+1)-y(i+1,j))
      sy=(x(i+1,j)-x(i+1,j+1))
      si=sqrt(sx*sx+sy*sy)
      sxp=sx/si
      syp=sy/si
c-------------------
c   Axisymmetric
c-------------------
      if(iaxi.eq.1) then
        factor=(y(i+1,j+1)+y(i+1,j))*faxi
        sx=sx*factor
        sy=sy*factor
        si=si*factor
      endif
c-----------------------------------------------------------------------
c          N O N - F L U X   S P L I T T I N G   O P T I O N
c-----------------------------------------------------------------------
      if(i.eq.1.and.icase.eq.0) then
c-------------------
c Transport velocity
c-------------------
        qsi=sx*u(i,j)+sy*v(i,j)
        qsip1=sx*u(i+1,j)+sy*v(i+1,j)
        piij=p(ii,j)
        qs=0.
c-------------------
c Flux vector, {F}i
c-------------------
        fs1(k2)=rho(ii,j,n)*qs
        fs2(k2)=rhou(ii,j,n)*qs+piij*sx
        fs3(k2)=rhov(ii,j,n)*qs+piij*sy
        fs4(k2)=(e(ii,j,n)+piij)*qs
        if(model.eq.0.or.model.eq.1) then
          fs5(k2)=rhok(ii,j,n)*qs
          fs6(k2)=rhomeg(ii,j,n)*qs
          if(model.eq.1) then
            fs7(k2)=rhotxx(ii,j,n)*qs
            fs8(k2)=rhotxy(ii,j,n)*qs
            fs9(k2)=rhotyy(ii,j,n)*qs
          endif
        endif
c-----------------------------------------------------------------------
c            F L U X   S P L I T T I N G   O P T I O N
c-----------------------------------------------------------------------
      else
        call fsifs(i,j,ii)
      endif
c=======================================================================
c           V I S C O U S   T R A N S P O R T   T E R M S
c=======================================================================
c-------------------
c    Gradients
c-------------------
      uiij=u(ii,j)
      viij=v(ii,j)
      due=u(i+1,j)-u(i,j)
      dve=v(i+1,j)-v(i,j)
      dee=ei(i+1,j)-ei(i,j)
      dxe=.25*(x(i+2,j+1)+x(i+2,j)-x(i,j+1)-x(i,j))
      dye=.25*(y(i+2,j+1)+y(i+2,j)-y(i,j+1)-y(i,j))
      dun=.25*(u(i,j+1)+u(i+1,j+1)-u(i,j-1)-u(i+1,j-1)) 
      dvn=.25*(v(i,j+1)+v(i+1,j+1)-v(i,j-1)-v(i+1,j-1)) 
      den=.25*(ei(i,j+1)+ei(i+1,j+1)-ei(i,j-1)-ei(i+1,j-1))
      dxn=x(i+1,j+1)-x(i+1,j)
      dyn=y(i+1,j+1)-y(i+1,j)
      dxy=1./(dxe*dyn-dxn*dye)
      dux=(due*dyn-dun*dye)*dxy
      dvx=(dve*dyn-dvn*dye)*dxy
      dex=(dee*dyn-den*dye)*dxy
      duy=(dun*dxe-due*dxn)*dxy
      dvy=(dvn*dxe-dve*dxn)*dxy
      dey=(den*dxe-dee*dxn)*dxy
c-------------------
c    Viscosities
c    -----------
c Suppress turbulent
c transport  for the
c implicit algorithm
c-------------------
      if(implt.eq.1) then
        thin=0.
      else
        thin=1.
      endif
      temp=abs(.5*(ei(i,j)+ei(i+1,j))/cv)
      rmu=viscos(temp)
      rkap=gmoprl*rmu
      rhok23=0.
      qx=0.
      qy=0.
      if(model.lt.0) then
        epi=0.
      elseif(model.eq.2) then
        epi=ep(i,j)*thin
      else
        epi=.5*(ep(i,j)+ep(i+1,j))*thin
        rhok23=.3333333*(rhok(i,j,n)+rhok(i+1,j,n))*thin
        dke=tke(i+1,j)-tke(i,j)
        doe=omeg(i+1,j)-omeg(i,j)
        dkn=.25*(tke(i,j+1)+tke(i+1,j+1)-tke(i,j-1)-tke(i+1,j-1))
        don=.25*(omeg(i,j+1)+omeg(i+1,j+1)-omeg(i,j-1)-omeg(i+1,j-1))
        dkx=(dke*dyn-dkn*dye)*dxy
        dox=(doe*dyn-don*dye)*dxy
        dky=(dkn*dxe-dke*dxn)*dxy
        doy=(don*dxe-doe*dxn)*dxy
        rmuk=rmu+sigmas*epi
        rmuo=rmu+sigma*epi
        diffk=rmuk*(dkx*sx+dky*sy)
        diffo=rmuo*(dox*sx+doy*sy)
        qx=rmuk*dkx
        qy=rmuk*dky
        if(model.eq.1) then
          dtxxe=txx(i+1,j)-txx(i,j)
          dtxye=txy(i+1,j)-txy(i,j)
          dtyye=tyy(i+1,j)-tyy(i,j)
          dtxxn=.25*(txx(i,j+1)+txx(i+1,j+1)-txx(i,j-1)-txx(i+1,j-1))
          dtxyn=.25*(txy(i,j+1)+txy(i+1,j+1)-txy(i,j-1)-txy(i+1,j-1))
          dtyyn=.25*(tyy(i,j+1)+tyy(i+1,j+1)-tyy(i,j-1)-tyy(i+1,j-1))
          dtxxx=(dtxxe*dyn-dtxxn*dye)*dxy
          dtxyx=(dtxye*dyn-dtxyn*dye)*dxy
          dtyyx=(dtyye*dyn-dtyyn*dye)*dxy
          dtxxy=(dtxxn*dxe-dtxxe*dxn)*dxy
          dtxyy=(dtxyn*dxe-dtxye*dxn)*dxy
          dtyyy=(dtyyn*dxe-dtyye*dxn)*dxy
          diffxx=rmuk*(dtxxx*sx+dtxxy*sy)
          diffxy=rmuk*(dtxyx*sx+dtxyy*sy)
          diffyy=rmuk*(dtyyx*sx+dtyyy*sy)
        endif
      endif
      rmu=rmu+fmodel*epi*climit(i,j)
      rkap=rkap+gmoprt*epi*climit(i,j)
      rlmbda=-.66666667*rmu
c-------------------
c     Stresses
c-------------------
      stress=rhok23-rlmbda*(dux+dvy)
      str(k2)=0.
      if(iaxi.eq.1.and.i.gt.1) then
        stress=stress-rlmbda*viij/(.5*(y(i+1,j)+y(i+1,j+1)))
        str(k2)=stress
      endif
      sigx=stress-2.*rmu*dux
      sigy=stress-2.*rmu*dvy
      tau=-rmu*(duy+dvx)
      if(model.eq.1) then
        if(ntau.eq.0) then
          sigx=sigx-.5*(rhotxx(i,j,n)+rhotxx(i+1,j,n))*thin
          sigy=sigy-.5*(rhotyy(i,j,n)+rhotyy(i+1,j,n))*thin
          tau=tau-.5*(rhotxy(i,j,n)+rhotxy(i+1,j,n))*thin
          if(iaxi.eq.1.and.i.ne.1)
     *       str(k2)=str(k2)+(rhotxx(ii,j,n)+rhotyy(ii,j,n))*thin
        else
          sigx=sigx-.5*(dtxx(i,j)+dtxx(i+1,j))*thin
          sigy=sigy-.5*(dtyy(i,j)+dtyy(i+1,j))*thin
          tau=tau-.5*(dtxy(i,j)+dtxy(i+1,j))*thin
          if(iaxi.eq.1.and.i.ne.1)
     *       str(k2)=str(k2)+(dtxx(ii,j)+dtyy(ii,j))*thin
        endif
      endif
      qx=qx+rkap*dex
      qy=qy+rkap*dey
c-------------------
c   Flux vector
c  {F}={F}i+{F}v
c-------------------
      fs2(k2)=fs2(k2)+sigx*sx+tau*sy
      fs3(k2)=fs3(k2)+tau*sx+sigy*sy
      fs4(k2)=fs4(k2)+(sigx*uiij+tau*viij-qx)*sx
     *               +(tau*uiij+sigy*viij-qy)*sy
      if(model.eq.0.or.model.eq.1) then
        fs5(k2)=fs5(k2)-diffk
        fs6(k2)=fs6(k2)-diffo
        if(model.eq.1) then
          fs7(k2)=fs7(k2)-diffxx
          fs8(k2)=fs8(k2)-diffxy
          fs9(k2)=fs9(k2)-diffyy
        endif
      endif
      return
      end
      subroutine fsifs(i,j,ii)
c=======================================================================
c         F L U X   S P L I T T I N G...X - D I R E C T I O N
c=======================================================================
      include 'come2c.'
c=======================================================================
      do 10 nn=1,2
c=======================================
c      Alternating-index method
c=======================================
        if(iave.eq.0) then
          rr=rho(ii,j,n)
          uu=u(ii,j)
          vv=v(ii,j)
          cc=sqrt(ggm1*abs(ei(ii,j)))
c=======================================
c          Averaging method
c=======================================
        else
          rr=.5*(rho(i,j,n)+rho(i+1,j,n))
          uu=.5*(rhou(i,j,n)+rhou(i+1,j,n))/rr
          vv=.5*(rhov(i,j,n)+rhov(i+1,j,n))/rr
          cc=sqrt(ggm1*abs(.5*(e(i,j,n)+e(i+1,j,n))/rr
     *                    -.5*(uu*uu+vv*vv)))
        endif
c=======================================
c   Flux splitting procedure...do left
c   moving contribution for  nn=1  and
c   right moving contribution for nn=2
c=======================================
        i1=i+2-nn
        i2=i1+3-2*nn
c=======================================
c     Check for boundary or shock
c=======================================
        if(i.eq.1.or.i.eq.ilm1) i2=i1
c=======================================
c Steger-Warming option (for stability)
c=======================================
        pmin=abs(p(i,j))
        pmin2=abs(p(i+1,j))
        if(pmin2.lt.pmin) pmin=pmin2
        pg=abs(p(i+1,j)-p(i,j))/pmin
        wt1=1./(1.+pg*pg)
        wt2=1.-wt1
        rr=wt1*rr+wt2*rho(i1,j,n)
        uu=wt1*uu+wt2*u(i1,j)
        vv=wt1*vv+wt2*v(i1,j)
        cc=wt1*cc+wt2*sqrt(ggm1*abs(ei(i1,j)))
        up=sxp*uu+syp*vv
        vp=sxp*vv-syp*uu
        alp=.5*(up*up+vp*vp)
c=======================================
c    Extrapolate for (i+1/2) plane
c=======================================
        pmin=abs(p(i1,j))
        pmin2=abs(p(i2,j))
        if(pmin2.lt.pmin) pmin=pmin2
        pg=pg+abs(p(i2,j)-p(i1,j))/pmin
        wt1=1./(1.+pg*pg)
        w1=1.+.5*wt1
        w2=-.5*wt1
        dr=w1*rho(i1,j,n)+w2*rho(i2,j,n)
        dru=w1*rhou(i1,j,n)+w2*rhou(i2,j,n)
        drv=w1*rhov(i1,j,n)+w2*rhov(i2,j,n)
        de=w1*e(i1,j,n)+w2*e(i2,j,n)
        dp=w1*p(i1,j)+w2*p(i2,j)
        if(model.eq.0.or.model.eq.1) then
          drhok=w1*rhok(i1,j,n)+w2*rhok(i2,j,n)
          dromeg=w1*rhomeg(i1,j,n)+w2*rhomeg(i2,j,n)
          if(model.eq.1) then
            drhtxx=w1*rhotxx(i1,j,n)+w2*rhotxx(i2,j,n)
            drhtxy=w1*rhotxy(i1,j,n)+w2*rhotxy(i2,j,n)
            drhtyy=w1*rhotyy(i1,j,n)+w2*rhotyy(i2,j,n)
          endif
        endif
c=======================================
c           {xx} = [S] {u}
c=======================================
        xx1=dr-dp/(cc*cc)
        xx2=(-up*dr+sxp*dru+syp*drv)*cc+dp
        xx3=(-vp*dr-syp*dru+sxp*drv)/rr
        xx4=-xx2+2.*dp
c=======================================
c           [dd] = [Lambda]
c=======================================
        dd1=si*up
        dd2=si*(up+cc)
        dd4=si*(up-cc)
c=======================================
c [yy] = [Lambda]- (Left  moving); nn=1
c [yy] = [Lambda]+ (Right moving); nn=2
c=======================================
        sgns=2*nn-3
        fact1=.5*(dd1+sgns*abs(dd1))
        fact2=.5*(dd2+sgns*abs(dd2))
        fact4=.5*(dd4+sgns*abs(dd4))
        yy1=fact1*xx1
        yy2=fact2*xx2
        yy3=fact1*xx3
        yy4=fact4*xx4
c=======================================
c        Auxiliary parameters
c=======================================
        t0=yy1+.5*(yy2+yy4)/(cc*cc)
        t1=up*t0+.5*(yy2-yy4)/cc
        t2=vp*t0+rr*yy3
c=======================================
c   Flux vector...
c     {F} = [S]inv [Lambda] [S] {u}
c=======================================
        fs1(k2)=(nn-1)*fs1(k2)+t0
        fs2(k2)=(nn-1)*fs2(k2)+sxp*t1-syp*t2
        fs3(k2)=(nn-1)*fs3(k2)+syp*t1+sxp*t2
        fs4(k2)=(nn-1)*fs4(k2)+alp*t0+rr*vp*yy3
     *         +.5*(up*(yy2-yy4)/cc+(yy2+yy4)/gm1)
        if(model.eq.0.or.model.eq.1) then
          fs5(k2)=(nn-1)*fs5(k2)+fact1*drhok
          fs6(k2)=(nn-1)*fs6(k2)+fact1*dromeg
          if(model.eq.1) then
            fs7(k2)=(nn-1)*fs7(k2)+fact1*drhtxx
            fs8(k2)=(nn-1)*fs8(k2)+fact1*drhtxy
            fs9(k2)=(nn-1)*fs9(k2)+fact1*drhtyy
          endif
        endif
 10   continue
c=======================================
c        Transport velocity
c=======================================
      rhomn=rho(i,j,n)
      rhomn2=rho(i+1,j,n)
      if(rhomn2.lt.rhomn) rhomn=rhomn2
      ut(i,j)=fs1(k2)/(si*rhomn)
      return
      end
      subroutine fsj(i,j,jj,dt)
c=======================================================================
c    T R A N S P O R T   A N D   S T R E S S...Y - D I R E C T I O N
c=======================================================================
      include 'come2c.'
c=======================================================================
c          I N V I S C I D   T R A N S P O R T   T E R M S
c=======================================================================
c   Cell geometry
c-------------------
      sx=(y(i,j+1)-y(i+1,j+1))
      sy=(x(i+1,j+1)-x(i,j+1))
      sj=sqrt(sx*sx+sy*sy)
      sxp=sx/sj
      syp=sy/sj
c-------------------
c   Axisymmetric
c-------------------
      if(iaxi.eq.1) then
        factor=(y(i+1,j+1)+y(i,j+1))*faxi
        sx=sx*factor
        sy=sy*factor
        sj=sj*factor
      endif
c-----------------------------------------------------------------------
c          N O N - F L U X   S P L I T T I N G   O P T I O N
c-----------------------------------------------------------------------
      if(j.eq.1) then
c-------------------
c Transport velocity
c-------------------
        qsj=sx*u(i,j)+sy*v(i,j)
        qsjp1=sx*u(i,j+1)+sy*v(i,j+1)
        pijj=p(i,jj)
        qs=0.
        pijj=pijj*(1.-(sxp*u(i,2)+syp*v(i,2))/sqrt(ggm1*ei(i,2)))**2
c-------------------
c Flux vector, {G}i
c-------------------
        fs1(k2)=rho(i,jj,n)*qs
        fs2(k2)=rhou(i,jj,n)*qs+pijj*sx
        fs3(k2)=rhov(i,jj,n)*qs+pijj*sy
        fs4(k2)=(e(i,jj,n)+pijj)*qs
        if(model.eq.0.or.model.eq.1) then
          fs5(k2)=rhok(i,jj,n)*qs
          fs6(k2)=rhomeg(i,jj,n)*qs
          if(model.eq.1) then
            fs7(k2)=rhotxx(i,jj,n)*qs
            fs8(k2)=rhotxy(i,jj,n)*qs
            fs9(k2)=rhotyy(i,jj,n)*qs
          endif
        endif
c-----------------------------------------------------------------------
c            F L U X   S P L I T T I N G   O P T I O N
c-----------------------------------------------------------------------
      else
        call fsjfs(i,j,jj,dt)
      endif
c=======================================================================
c           V I S C O U S   T R A N S P O R T   T E R M S
c=======================================================================
c-------------------
c    Gradients
c-------------------
      uijj=.5*(u(i,j)+u(i,j+1))
      vijj=.5*(v(i,j)+v(i,j+1))
      due=.25*(u(i+1,j)+u(i+1,j+1)-u(i-1,j)-u(i-1,j+1))
      dve=.25*(v(i+1,j)+v(i+1,j+1)-v(i-1,j)-v(i-1,j+1))
      dee=.25*(ei(i+1,j)+ei(i+1,j+1)-ei(i-1,j)-ei(i-1,j+1))
      dxe=x(i+1,j+1)-x(i,j+1)
      dye=y(i+1,j+1)-y(i,j+1)
      dun=u(i,j+1)-u(i,j)
      dvn=v(i,j+1)-v(i,j)
      den=ei(i,j+1)-ei(i,j)
      dxn=.25*(x(i+1,j+2)+x(i,j+2)-x(i+1,j)-x(i,j))
      dyn=.25*(y(i+1,j+2)+y(i,j+2)-y(i+1,j)-y(i,j))
      dxy=1./(dxe*dyn-dxn*dye)
      dux=(due*dyn-dun*dye)*dxy
      dvx=(dve*dyn-dvn*dye)*dxy
      dex=(dee*dyn-den*dye)*dxy
      duy=(dun*dxe-due*dxn)*dxy
      dvy=(dvn*dxe-dve*dxn)*dxy
      dey=(den*dxe-dee*dxn)*dxy
c-------------------
c   Viscosities
c-------------------
      temp=abs(.5*(ei(i,j)+ei(i,j+1))/cv)
      rmu=viscos(temp)
      rkap=gmoprl*rmu
      rhok23=0.
      qx=0.
      qy=0.
      if(model.lt.0) then
        epj=0.
      elseif(model.eq.2) then
        epj=ep(i,j)
      else
        epj=.5*(ep(i,j)+ep(i,j+1))
        if(j.ne.1) rhok23=.3333333*(rhok(i,j,n)+rhok(i,j+1,n))
        dke=.25*(tke(i+1,j)+tke(i+1,j+1)-tke(i-1,j)-tke(i-1,j+1))
        doe=.25*(omeg(i+1,j)+omeg(i+1,j+1)-omeg(i-1,j)-omeg(i-1,j+1))
        dkn=tke(i,j+1)-tke(i,j)
        don=omeg(i,j+1)-omeg(i,j)
        dkx=(dke*dyn-dkn*dye)*dxy
        dox=(doe*dyn-don*dye)*dxy
        dky=(dkn*dxe-dke*dxn)*dxy
        doy=(don*dxe-doe*dxn)*dxy
        rmuk=rmu+sigmas*epj
        rmuo=rmu+sigma*epj
        diffk=rmuk*(dkx*sx+dky*sy)
        diffo=rmuo*(dox*sx+doy*sy)
        qx=rmuk*dkx
        qy=rmuk*dky
        if(model.eq.1) then
          dtxxe=.25*(txx(i+1,j)+txx(i+1,j+1)-txx(i-1,j)-txx(i-1,j+1))
          dtxye=.25*(txy(i+1,j)+txy(i+1,j+1)-txy(i-1,j)-txy(i-1,j+1))
          dtyye=.25*(tyy(i+1,j)+tyy(i+1,j+1)-tyy(i-1,j)-tyy(i-1,j+1))
          dtxxn=txx(i,j+1)-txx(i,j)
          dtxyn=txy(i,j+1)-txy(i,j)
          dtyyn=tyy(i,j+1)-tyy(i,j)
          dtxxx=(dtxxe*dyn-dtxxn*dye)*dxy
          dtxyx=(dtxye*dyn-dtxyn*dye)*dxy
          dtyyx=(dtyye*dyn-dtyyn*dye)*dxy
          dtxxy=(dtxxn*dxe-dtxxe*dxn)*dxy
          dtxyy=(dtxyn*dxe-dtxye*dxn)*dxy
          dtyyy=(dtyyn*dxe-dtyye*dxn)*dxy
          diffxx=rmuk*(dtxxx*sx+dtxxy*sy)
          diffxy=rmuk*(dtxyx*sx+dtxyy*sy)
          diffyy=rmuk*(dtyyx*sx+dtyyy*sy)
        endif
      endif
      rmu=rmu+fmodel*epj*climit(i,j)
      rkap=rkap+gmoprt*epj*climit(i,j)
      rlmbda=-.66666667*rmu
c-------------------
c     Stresses
c-------------------
      stress=rhok23-rlmbda*(dux+dvy)
      if(iaxi.eq.1) then
        rad=.5*(y(i,j+1)+y(i+1,j+1))
        if(rad.gt.0.) stress=stress-rlmbda*vijj/rad
        str(k2)=stress
      endif
      sigx=stress-2.*rmu*dux
      sigy=stress-2.*rmu*dvy
      tau=-rmu*(duy+dvx)
      if(model.eq.1.and.j.ne.1) then
        if(ntau.eq.0) then
          sigx=sigx-.5*(rhotxx(i,j,n)+rhotxx(i,j+1,n))
          sigy=sigy-.5*(rhotyy(i,j,n)+rhotyy(i,j+1,n))
          tau=tau-.5*(rhotxy(i,j,n)+rhotxy(i,j+1,n))
          if(iaxi.eq.1) str(k2)=str(k2)+rhotxx(i,jj,n)+rhotyy(i,jj,n)
        else
          sigx=sigx-.5*(dtxx(i,j)+dtxx(i,j+1))
          sigy=sigy-.5*(dtyy(i,j)+dtyy(i,j+1))
          tau=tau-.5*(dtxy(i,j)+dtxy(i,j+1))
          if(iaxi.eq.1) str(k2)=str(k2)+dtxx(i,jj)+dtyy(i,jj)
        endif
      endif
      qx=qx+rkap*dex
      qy=qy+rkap*dey
c-------------------
c   Flux vector
c  {G}={G}i+{G}v
c-------------------
      fs2(k2)=fs2(k2)+sigx*sx+tau*sy
      fs3(k2)=fs3(k2)+tau*sx+sigy*sy
      fs4(k2)=fs4(k2)+(sigx*uijj+tau*vijj-qx)*sx
     *               +(tau*uijj+sigy*vijj-qy)*sy
      if(model.eq.0.or.model.eq.1) then
        fs5(k2)=fs5(k2)-diffk
        fs6(k2)=fs6(k2)-diffo
        if(model.eq.1) then
          fs7(k2)=fs7(k2)-diffxx
          fs8(k2)=fs8(k2)-diffxy
          fs9(k2)=fs9(k2)-diffyy
        endif
      endif
      return
      end
      subroutine fsjfs(i,j,jx,dt)
c=======================================================================
c         F L U X   S P L I T T I N G...Y - D I R E C T I O N
c=======================================================================
      include 'come2c.'
c=======================================================================
      vt(i,j)=0.
      jj=jx
      do 10 nn=1,2
c=======================================
c      Alternating-index method
c=======================================
        if(iave.eq.0) then
          rr=rho(i,jj,n)
          uu=u(i,jj)
          vv=v(i,jj)
          cc=sqrt(ggm1*abs(ei(i,jj)))
c=======================================
c          Averaging method
c=======================================
        else
          rr=.5*(rho(i,j,n)+rho(i,j+1,n))
          uu=.5*(rhou(i,j,n)+rhou(i,j+1,n))/rr
          vv=.5*(rhov(i,j,n)+rhov(i,j+1,n))/rr
          cc=sqrt(ggm1*abs(.5*(e(i,j,n)+e(i,j+1,n))/rr
     *                    -.5*(uu*uu+vv*vv)))
        endif
c=======================================
c   Flux  splitting  procedure...do up
c   moving contribution for  nn=1  and
c   down moving  contribution for nn=2
c=======================================
        j1=j+2-nn
        j2=j1+3-2*nn
c=======================================
c     Check for boundary or shock
c=======================================
        if(j.ge.jl+nn-3) j2=j1
c=======================================
c Steger-Warming option (for stability)
c=======================================
        pmin=abs(p(i,j))
        pmin2=abs(p(i,j+1))
        if(pmin2.lt.pmin) pmin=pmin2
        pg=abs(p(i,j+1)-p(i,j))/pmin
        wt1=1./(1.+pg*pg)
        wt2=1.-wt1
        rr=wt1*rr+wt2*rho(i,j1,n)
        uu=wt1*uu+wt2*u(i,j1)
        vv=wt1*vv+wt2*v(i,j1)
        cc=wt1*cc+wt2*sqrt(ggm1*abs(ei(i,j1)))
        wsj=0.
        vsj=sxp*u(i,j+1)+syp*v(i,j+1)
        csj=sqrt(ggm1*ei(i,j+1))
        if(-vsj.gt.csj) then
          wsh=-vsj-sqrt((gp1*p(i,j)/p(i,j+1)+gm1)/(2.*gam))*csj
          wsj=.5*(wsh-abs(wsh))*sj
          rr=rho(i,j+1,n)
          uu=u(i,j+1)
          vv=v(i,j+1)
          cc=csj
          pg=0.
        endif
        up=syp*uu-sxp*vv
        vp=sxp*uu+syp*vv
        alp=.5*(up*up+vp*vp)
c=======================================
c    Extrapolate for (j+1/2) plane
c=======================================
        pmin=abs(p(i,j1))
        pmin2=abs(p(i,j2))
        if(pmin2.lt.pmin) pmin=pmin2
        pg=pg+abs(p(i,j2)-p(i,j1))/pmin
        wt1=1./(1.+pg*pg)
        w1=1.+.5*wt1
        w2=-.5*wt1
        dr=w1*rho(i,j1,n)+w2*rho(i,j2,n)
        dru=w1*rhou(i,j1,n)+w2*rhou(i,j2,n)
        drv=w1*rhov(i,j1,n)+w2*rhov(i,j2,n)
        de=w1*e(i,j1,n)+w2*e(i,j2,n)
        dp=w1*p(i,j1)+w2*p(i,j2)
        if(model.eq.0.or.model.eq.1) then
          drhok=w1*rhok(i,j1,n)+w2*rhok(i,j2,n)
          dromeg=w1*rhomeg(i,j1,n)+w2*rhomeg(i,j2,n)
          if(model.eq.1) then
            drhtxx=w1*rhotxx(i,j1,n)+w2*rhotxx(i,j2,n)
            drhtxy=w1*rhotxy(i,j1,n)+w2*rhotxy(i,j2,n)
            drhtyy=w1*rhotyy(i,j1,n)+w2*rhotyy(i,j2,n)
          endif
        endif
c=======================================
c           {xx} = [S] {u}
c=======================================
        xx1=dr-dp/(cc*cc)
        xx2=(-up*dr+syp*dru-sxp*drv)/rr
        xx3=(-vp*dr+sxp*dru+syp*drv)*cc+dp
        xx4=-xx3+2.*dp
c=======================================
c           [dd] = [Lambda]
c=======================================
        dd1=sj*vp
        dd3=sj*(vp+cc)
        dd4=sj*(vp-cc)
c=======================================
c [yy] = [Lambda]- (Up   moving); nn=1
c [yy] = [Lambda]+ (Down moving); nn=2
c=======================================
        sgns=2*nn-3
        fact1=.5*(dd1+sgns*abs(dd1))-sgns*wsj
        fact3=.5*(dd3+sgns*abs(dd3))-sgns*wsj
        fact4=.5*(dd4+sgns*abs(dd4))-sgns*wsj
        yy1=fact1*xx1
        yy2=fact1*xx2
        yy3=fact3*xx3
        yy4=fact4*xx4
c=======================================
c        Auxiliary parameters
c=======================================
        t0=yy1+.5*(yy3+yy4)/(cc*cc)
        t1=up*t0+rr*yy2
        t2=vp*t0+.5*(yy3-yy4)/cc
c=======================================
c   Flux vector...
c     {G} = [S]inv [Lambda] [S] {u}
c=======================================
        fs1(k2)=(nn-1)*fs1(k2)+t0
        fs2(k2)=(nn-1)*fs2(k2)+syp*t1+sxp*t2
        fs3(k2)=(nn-1)*fs3(k2)-sxp*t1+syp*t2
        fs4(k2)=(nn-1)*fs4(k2)
     *         +alp*t0+rr*up*yy2+.5*(vp*(yy3-yy4)/cc+(yy3+yy4)/gm1)
        if(model.eq.0.or.model.eq.1) then
          fs5(k2)=(nn-1)*fs5(k2)+fact1*drhok
          fs6(k2)=(nn-1)*fs6(k2)+fact1*dromeg
          if(model.eq.1) then
            fs7(k2)=(nn-1)*fs7(k2)+fact1*drhtxx
            fs8(k2)=(nn-1)*fs8(k2)+fact1*drhtxy
            fs9(k2)=(nn-1)*fs9(k2)+fact1*drhtyy
          endif
        endif
 10   continue
c=======================================
c        Transport velocity
c=======================================
      rhomn=rho(i,j,n)
      rhomn2=rho(i,j+1,n)
      if(rhomn2.lt.rhomn) rhomn=rhomn2
      vt(i,j)=fs1(k2)/(sj*rhomn)
      return
      end
      subroutine impli(dt,i,j,iadd)
c=======================================================================
c   I M P L I C I T   P R O C E D U R E   I N   I - D I R E C T I O N
c=======================================================================
      include 'come2c.'
c=======================================================================
      sgn=1.
      dtvol=dt/vol(i,j)
      do 40 nn=1,2
c=======================================
c           Cell geometry
c=======================================
        ii=i+2-nn
        si=sqrt((x(ii,j+1)-x(ii,j))**2+(y(ii,j+1)-y(ii,j))**2)
        sxp=(y(ii,j+1)-y(ii,j))/si
        syp=(x(ii,j)-x(ii,j+1))/si
        if(iaxi.eq.1) si=si*(y(ii,j+1)+y(ii,j))*faxi
        if(ni.eq.2.and.nn.eq.1) then
          iflag=2
        else
          iflag=1
        endif
c=======================================
c      Alternating-index method
c=======================================
 10     if(iave.eq.0) then
          ii=i+1-nn+iadd
          rr=rho(ii,j,n)
          uu=u(ii,j)
          vv=v(ii,j)
          cc=sqrt(ggm1*abs(ei(ii,j)))
c=======================================
c          Averaging method
c=======================================
        else
          ii=i+1-nn
          rr=.5*(rho(ii,j,n)+rho(ii+1,j,n))
          uu=.5*(rhou(ii,j,n)+rhou(ii+1,j,n))/rr
          vv=.5*(rhov(ii,j,n)+rhov(ii+1,j,n))/rr
          cc=sqrt(ggm1*abs(.5*(e(ii,j,n)+e(ii+1,j,n))/rr
     *                    -.5*(uu*uu+vv*vv)))
        endif
c=======================================
c         Auxiliary parameters
c=======================================
        rc=rr*cc
        i1=i+2-nn+1-iflag
        rri=rho(i1,j,n)
        uui=u(i1,j)
        vvi=v(i1,j)
        cci=sqrt(ggm1*abs(ei(i1,j)))
        upi=syp*vvi+sxp*uui
        vpi=sxp*vvi-syp*uui
        alpi=.5*(upi*upi+vpi*vpi)
c=======================================
c Steger-Warming option (for stability)
c=======================================
        i1=i+2-nn
        i2=i1-1
        pmin=abs(p(i1,j))
        pmin2=abs(p(i2,j))
        if(pmin2.lt.pmin) pmin=pmin2
        pg=abs(p(i1,j)-p(i2,j))/pmin
        wt1=1./(1.+pg*pg)
        wt2=1.-wt1
        rr=wt1*rr+wt2*rri
        uu=wt1*uu+wt2*uui
        vv=wt1*vv+wt2*vvi
        cc=wt1*cc+wt2*cci
        up=syp*vv+sxp*uu
        vp=sxp*vv-syp*uu
        alp=.5*(up*up+vp*vp)
        rccsq=1./(cc*cc)
        dd=-1.5*sgn*(3-2*iflag)*si*dtvol
c=======================================
c          [dd] = [Lambda]
c=======================================
        dd1=up*dd
        dd2=(up+cc)*dd
        dd4=(up-cc)*dd
        coef=1.
        upi1=sxp*u(i1,j)+syp*v(i1,j)
        upi2=sxp*u(i2,j)+syp*v(i2,j)
        cci1=sqrt(ggm1*ei(i1,j))
        cci2=sqrt(ggm1*ei(i2,j))
        ep0=coef*abs(dd*(up-ut(i1-1,j)))
     *     +coef*abs(dd*(p(i1,j)-p(i2,j)))/rc
        ep1=coef*abs(dd*(upi1-upi2))+ep0
        ep2=coef*abs(dd*((upi1+cci1)-(upi2+cci2)))+ep0
        ep4=coef*abs(dd*((upi1-cci1)-(upi2-cci2)))+ep0
        dd1=.5*(dd1+sqrt(dd1*dd1+ep1*ep1))
        dd2=.5*(dd2+sqrt(dd2*dd2+ep2*ep2))
        dd4=.5*(dd4+sqrt(dd4*dd4+ep4*ep4))
c=======================================
c        [t] = [Lambda] [S]
c=======================================
        uu=uui
        vv=vvi
        alp=alpi
        t11=dd1*(1.-alp*betam*rccsq)
        t12=dd1*uu*betam*rccsq
        t13=dd1*vv*betam*rccsq
        t14=-dd1*betam*rccsq
        t21=dd2*(alp*betam-up*cc)
        t22=dd2*(sxp*cc-uu*betam)
        t23=dd2*(syp*cc-vv*betam)
        t24=dd2*betam
        t31=-dd1*vp/rr
        t32=-dd1*syp/rr
        t33=dd1*sxp/rr
        t34=0.
        t41=dd4*(alp*betam+up*cc)
        t42=dd4*(-sxp*cc-uu*betam)
        t43=dd4*(-syp*cc-vv*betam)
        t44=dd4*betam
        uu=sxp*up-syp*vp
        vv=syp*up+sxp*vp
        alp=.5*(up*up+vp*vp)
        if(iflag.eq.1) then
c=======================================
c          [bp] = [S]inv
c=======================================
          bp11=1.
          bp12=.5*rccsq
          bp13=0.
          bp14=.5*rccsq
          bp21=uu
          bp22=.5*(uu+sxp*cc)*rccsq
          bp23=-syp*rr
          bp24=.5*(uu-sxp*cc)*rccsq
          bp31=vv
          bp32=.5*(vv+syp*cc)*rccsq
          bp33=sxp*rr
          bp34=.5*(vv-syp*cc)*rccsq
          bp41=alp
          bp42=.5*((alp+up*cc)*rccsq+1./betam)
          bp43=rr*vp
          bp44=.5*((alp-up*cc)*rccsq+1./betam)
c=======================================
c       [b] = [S]inv [Lambda] [S]
c=======================================
          b11x=-(bp11*t11+bp12*t21+bp13*t31+bp14*t41)
          b12x=-(bp11*t12+bp12*t22+bp13*t32+bp14*t42)
          b13x=-(bp11*t13+bp12*t23+bp13*t33+bp14*t43)
          b14x=-(bp11*t14+bp12*t24+bp13*t34+bp14*t44)
          b21x=-(bp21*t11+bp22*t21+bp23*t31+bp24*t41)
          b22x=-(bp21*t12+bp22*t22+bp23*t32+bp24*t42)
          b23x=-(bp21*t13+bp22*t23+bp23*t33+bp24*t43)
          b24x=-(bp21*t14+bp22*t24+bp23*t34+bp24*t44)
          b31x=-(bp31*t11+bp32*t21+bp33*t31+bp34*t41)
          b32x=-(bp31*t12+bp32*t22+bp33*t32+bp34*t42)
          b33x=-(bp31*t13+bp32*t23+bp33*t33+bp34*t43)
          b34x=-(bp31*t14+bp32*t24+bp33*t34+bp34*t44)
          b41x=-(bp41*t11+bp42*t21+bp43*t31+bp44*t41)
          b42x=-(bp41*t12+bp42*t22+bp43*t32+bp44*t42)
          b43x=-(bp41*t13+bp42*t23+bp43*t33+bp44*t43)
          b44x=-(bp41*t14+bp42*t24+bp43*t34+bp44*t44)
          iflag=2
          ab=1.
          ac=0.
          if(nn.eq.1) go to 20
          go to 30
        else
c=======================================
c           [cp] = [S]inv
c=======================================
          cp11=1.
          cp12=.5*rccsq
          cp13=0.
          cp14=.5*rccsq
          cp21=uu
          cp22=.5*(uu+sxp*cc)*rccsq
          cp23=-syp*rr
          cp24=.5*(uu-sxp*cc)*rccsq
          cp31=vv
          cp32=.5*(vv+syp*cc)*rccsq
          cp33=sxp*rr
          cp34=.5*(vv-syp*cc)*rccsq
          cp41=alp
          cp42=.5*((alp+up*cc)*rccsq+1./betam)
          cp43=rr*vp
          cp44=.5*((alp-up*cc)*rccsq+1./betam)
c=======================================
c       [c] = [S]inv [Lambda] [S]
c=======================================
          c11x=-(cp11*t11+cp12*t21+cp13*t31+cp14*t41)
          c12x=-(cp11*t12+cp12*t22+cp13*t32+cp14*t42)
          c13x=-(cp11*t13+cp12*t23+cp13*t33+cp14*t43)
          c14x=-(cp11*t14+cp12*t24+cp13*t34+cp14*t44)
          c21x=-(cp21*t11+cp22*t21+cp23*t31+cp24*t41)
          c22x=-(cp21*t12+cp22*t22+cp23*t32+cp24*t42)
          c23x=-(cp21*t13+cp22*t23+cp23*t33+cp24*t43)
          c24x=-(cp21*t14+cp22*t24+cp23*t34+cp24*t44)
          c31x=-(cp31*t11+cp32*t21+cp33*t31+cp34*t41)
          c32x=-(cp31*t12+cp32*t22+cp33*t32+cp34*t42)
          c33x=-(cp31*t13+cp32*t23+cp33*t33+cp34*t43)
          c34x=-(cp31*t14+cp32*t24+cp33*t34+cp34*t44)
          c41x=-(cp41*t11+cp42*t21+cp43*t31+cp44*t41)
          c42x=-(cp41*t12+cp42*t22+cp43*t32+cp44*t42)
          c43x=-(cp41*t13+cp42*t23+cp43*t33+cp44*t43)
          c44x=-(cp41*t14+cp42*t24+cp43*t34+cp44*t44)
          ab=0.
          ac=1.
          if(nn.eq.1) go to 30
          if(ni.ne.1) then
c=======================================
c        {du} = {du} - [c] {du}
c=======================================
            i1=i-1
            du1=du(i1,j,1)
            du2=du(i1,j,2)
            du3=du(i1,j,3)
            du4=du(i1,j,4)
            du(i,j,1)=du(i,j,1)-(c11x*du1+c12x*du2+c13x*du3+c14x*du4)
            du(i,j,2)=du(i,j,2)-(c21x*du1+c22x*du2+c23x*du3+c24x*du4)
            du(i,j,3)=du(i,j,3)-(c31x*du1+c32x*du2+c33x*du3+c34x*du4)
            du(i,j,4)=du(i,j,4)-(c41x*du1+c42x*du2+c43x*du3+c44x*du4)
            if(model.eq.0.or.model.eq.1) then
              du(i,j,5)=du(i,j,5)+dd1*du(i1,j,5)
              du(i,j,6)=du(i,j,6)+dd1*du(i1,j,6)
              if(model.eq.1) then
                du(i,j,7)=du(i,j,7)+dd1*du(i1,j,7)
                du(i,j,8)=du(i,j,8)+dd1*du(i1,j,8)
                du(i,j,9)=du(i,j,9)+dd1*du(i1,j,9)
              endif
            endif
          endif
          go to 40
        endif
c=======================================
c        {du} = {du} - [b] {du}
c=======================================
 20     i1=i+1
        du1=du(i1,j,1)
        du2=du(i1,j,2)
        du3=du(i1,j,3)
        du4=du(i1,j,4)
        bdu1=b11x*du1+b12x*du2+b13x*du3+b14x*du4
        bdu2=b21x*du1+b22x*du2+b23x*du3+b24x*du4
        bdu3=b31x*du1+b32x*du2+b33x*du3+b34x*du4
        bdu4=b41x*du1+b42x*du2+b43x*du3+b44x*du4
        du(i,j,1)=du(i,j,1)-bdu1
        du(i,j,2)=du(i,j,2)-bdu2
        du(i,j,3)=du(i,j,3)-bdu3
        du(i,j,4)=du(i,j,4)-bdu4
        if(model.eq.0.or.model.eq.1) then
          bdu5=-dd1*du(i1,j,5)
          bdu6=-dd1*du(i1,j,6)
          du(i,j,5)=du(i,j,5)-bdu5
          du(i,j,6)=du(i,j,6)-bdu6
          if(model.eq.1) then
            bdu7=-dd1*du(i1,j,7)
            bdu8=-dd1*du(i1,j,8)
            bdu9=-dd1*du(i1,j,9)
            du(i,j,7)=du(i,j,7)-bdu7
            du(i,j,8)=du(i,j,8)-bdu8
            du(i,j,9)=du(i,j,9)-bdu9
          endif
        endif
c=======================================
c       {dus} = {dus} - [b] {du}
c=======================================
        dus(i,j,1)=dus(i,j,1)-bdu1
        dus(i,j,2)=dus(i,j,2)-bdu2
        dus(i,j,3)=dus(i,j,3)-bdu3
        dus(i,j,4)=dus(i,j,4)-bdu4
        if(model.eq.0.or.model.eq.1) then
          dus(i,j,5)=dus(i,j,5)-bdu5
          dus(i,j,6)=dus(i,j,6)-bdu6
          if(model.eq.1) then
            dus(i,j,7)=dus(i,j,7)-bdu7
            dus(i,j,8)=dus(i,j,8)-bdu8
            dus(i,j,9)=dus(i,j,9)-bdu9
          endif
        endif
        go to 10
c=======================================
c  [a] = [a] - [b]  or  [a] = [a] - [c]
c=======================================
 30     a11=a11-ab*b11x-ac*c11x
        a12=a12-ab*b12x-ac*c12x
        a13=a13-ab*b13x-ac*c13x
        a14=a14-ab*b14x-ac*c14x
        a21=a21-ab*b21x-ac*c21x
        a22=a22-ab*b22x-ac*c22x
        a23=a23-ab*b23x-ac*c23x
        a24=a24-ab*b24x-ac*c24x
        a31=a31-ab*b31x-ac*c31x
        a32=a32-ab*b32x-ac*c32x
        a33=a33-ab*b33x-ac*c33x
        a34=a34-ab*b34x-ac*c34x
        a41=a41-ab*b41x-ac*c41x
        a42=a42-ab*b42x-ac*c42x
        a43=a43-ab*b43x-ac*c43x
        a44=a44-ab*b44x-ac*c44x
        if(model.eq.0.or.model.eq.1) then
          a55=a55+dd1
          a66=a66+dd1
          if(model.eq.1) then
            a77=a77+dd1
            a88=a88+dd1
            a99=a99+dd1
          endif
        endif
c=======================================
c           Axisymmetric
c=======================================
        if(iaxi.eq.1.and.nn.eq.1) then
          add=-1.5*sgn*(-(x(i+1,j+1)-x(i+1,j))*(y(i+1,j+1)+y(i+1,j))
     *       +(x(i,j+1)-x(i,j))*(y(i,j+1)+y(i,j)))*faxi*betam*dtvol
          a31=a31+add*(.5*(u(i,j)**2+v(i,j)**2))
          a32=a32+add*(-u(i,j))
          a33=a33+add*(-v(i,j))
          a34=a34+add
          rmu=0.
          rad=.25*(y(i,j)+y(i+1,j)+y(i,j+1)+y(i+1,j+1))
          a33=a33-add*(rmu/rad)/betam
        endif
        if(nn.eq.2) go to 10
 40   continue
      return
      end
      subroutine implj(dt,i,j,jadd)
c=======================================================================
c   I M P L I C I T   P R O C E D U R E   I N   J - D I R E C T I O N
c=======================================================================
      include 'come2c.'
c=======================================================================
      sgn=1.
c=======================================
c            Cell geometry
c=======================================
      sj=sqrt((x(i+1,j)-x(i,j))**2+(y(i+1,j)-y(i,j))**2)
      syp=(x(i+1,j)-x(i,j))/sj
      sxp=(y(i,j)-y(i+1,j))/sj
      if(iaxi.eq.1) sj=sj*(y(i+1,j)+y(i,j))*faxi
c=======================================
c        Viscosity coefficients
c=======================================
      temp=abs(.5*(ei(i,j)+ei(i,j-1))/cv)
      rmu=viscos(temp)
      rkap=gmoprl*rmu
      if(model.lt.0) then
        epj=0.
      elseif(model.eq.2) then
        epj=ep(i,j-1)
      else
        epj=.5*(ep(i,j)+ep(i,j-1))
        rmuk=rmu+sigmas*epj
        rmuo=rmu+sigma*epj
      endif
      rmu=rmu+epj*climit(i,j)
      rkap=rkap+gmoprt*epj*climit(i,j)
      rlmbda=-.66666667*rmu
      if(j.eq.2) then
        jj=j
      else
        jj=j-1+jadd
      endif
      dtvol=dt/vol(i,j)
      upbar=.5*(syp*(u(i,j)+u(i,j-1))-sxp*(v(i,j)+v(i,j-1)))
      vpbar=.5*(sxp*(u(i,j)+u(i,j-1))+syp*(v(i,j)+v(i,j-1)))
      iflag=1
c=======================================
c      Alternating-index method
c=======================================
 10   if(iave.eq.0) then
        rr=rho(i,jj,n)
        uu=u(i,jj)
        vv=v(i,jj)
        cc=sqrt(ggm1*abs(ei(i,jj)))
c=======================================
c          Averaging method
c=======================================
      else
        rr=.5*(rho(i,j,n)+rho(i,j-1,n))
        uu=.5*(rhou(i,j,n)+rhou(i,j-1,n))/rr
        vv=.5*(rhov(i,j,n)+rhov(i,j-1,n))/rr
        cc=sqrt(ggm1*abs(.5*(e(i,j,n)+e(i,j-1,n))/rr
     *                  -.5*(uu*uu+vv*vv)))
      endif
c=======================================
c         Auxiliary parameters
c=======================================
      rc=rr*cc
      j1=j+1-iflag
      rrj=p(i,j1)/(gm1*ei(i,j1))
      uuj=u(i,j1)
      vvj=v(i,j1)
      upj=syp*uuj-sxp*vvj
      vpj=sxp*uuj+syp*vvj
      ccj=sqrt(ggm1*abs(ei(i,j1)))
      alpj=.5*(uuj*uuj+vvj*vvj)
c=======================================
c          Check for shock
c=======================================
      j1=j
      j2=j-1
      pmin=abs(p(i,j1))
      pmin2=abs(p(i,j2))
      if(pmin2.lt.pmin) pmin=pmin2
      pg=abs(p(i,j1)-p(i,j2))/pmin
      wt1=1./(1.+pg*pg)
      wt2=1.-wt1
      rr=wt1*rr+wt2*rrj
      uu=wt1*uu+wt2*uuj
      vv=wt1*vv+wt2*vvj
      cc=wt1*cc+wt2*ccj
      coef=1.
      wsj=0.
      vsj=sxp*u(i,j)+syp*v(i,j)
      csj=sqrt(ggm1*ei(i,j))
      if(-vsj.gt.csj) then
        wsj=-vsj-sqrt((gp1*p(i,j-1)/p(i,j)+gm1)/(2.*gam))*csj
        wsj=.5*(wsj-abs(wsj))
        rr=rho(i,j,n)
        uu=u(i,j)
        vv=v(i,j)
        cc=csj
        coef=0.
      endif
      up=syp*uu-sxp*vv
      vp=sxp*uu+syp*vv
      alp=.5*(up*up+vp*vp)
      rccsq=1./(cc*cc)
      dd=-1.5*sgn*(3-2*(iflag))*sj*dtvol
c=======================================
c    Boundary condition parameters
c=======================================
      ttb=1.
      ttc=1.
      tt=0.
      if(j.eq.2) then
        if(iflag.eq.2) go to 20
        dd=0.
        ttc=0.
        ttu=1.
        tte=1-2*iadbwl
        tt=1.
        if(model.eq.0.or.model.eq.1) then
          ttk=1.
          tto=1.
        endif
      endif
c=======================================
c           [dd] = [Lambda]
c=======================================
      dd1=vp*dd
      dd3=(vp+cc)*dd
      dd4=(vp-cc)*dd
      j1=j+1-iflag
      j2=j-2+iflag
      vpj1=sxp*u(i,j1)+syp*v(i,j1)
      vpj2=sxp*u(i,j2)+syp*v(i,j2)
      ccj1=sqrt(ggm1*ei(j1,j))
      ccj2=sqrt(ggm1*ei(j2,j))
      ep0=coef*abs(dd*(vp-vt(i,j-1)))
     *   +coef*abs(dd*(p(i,j1)-p(i,j2)))/rc
      ep1=coef*abs(dd*(vpj1-vpj2))+ep0
      ep3=coef*abs(dd*((vpj1+ccj1)-(vpj2+ccj2)))+ep0
      ep4=coef*abs(dd*((vpj1-ccj1)-(vpj2-ccj2)))+ep0
      dd1=.5*(dd1+sqrt(dd1*dd1+ep1*ep1))+abs(dd*wsj)
      dd3=.5*(dd3+sqrt(dd3*dd3+ep3*ep3))+abs(dd*wsj)
      dd4=.5*(dd4+sqrt(dd4*dd4+ep4*ep4))+abs(dd*wsj)
c=======================================
c        [t] = [Lambda] [S]
c=======================================
      uu=uuj
      vv=vvj
      alp=alpj
      t11=dd1*(1.-alp*betam*rccsq)
      t12=dd1*uu*betam*rccsq
      t13=dd1*vv*betam*rccsq
      t14=-dd1*betam*rccsq
      t21=-dd1*up/rr
      t22=dd1*syp/rr
      t23=-dd1*sxp/rr
      t24=0.
      t31=dd3*(alp*betam-vp*cc)
      t32=dd3*( sxp*cc-uu*betam)
      t33=dd3*( syp*cc-vv*betam)
      t34=dd3*betam
      t41=dd4*(alp*betam+vp*cc)
      t42=dd4*(-sxp*cc-uu*betam)
      t43=dd4*(-syp*cc-vv*betam)
      t44=dd4*betam
      uu=sxp*vp+syp*up
      vv=syp*vp-sxp*up
      alp=.5*(up*up+vp*vp)
      if(iflag.eq.1) then
c=======================================
c          [bp] = [S]inv
c=======================================
        bp11=1.
        bp12=0.
        bp13=.5*rccsq
        bp14=.5*rccsq
        bp21=uu
        bp22=syp*rr
        bp23=(uu+sxp*cc)*.5*rccsq
        bp24=(uu-sxp*cc)*.5*rccsq
        bp31=vv
        bp32=-sxp*rr
        bp33=(vv+syp*cc)*.5*rccsq
        bp34=(vv-syp*cc)*.5*rccsq
        bp41=alp
        bp42=rr*up
        bp43=.5*((alp+vp*cc)*rccsq+1./betam)
        bp44=.5*((alp-vp*cc)*rccsq+1./betam)
c=======================================
c  [b] = [S]inv [Lambda] [S] + [M] [N]
c         (Centered at j+1/2)
c=======================================
        vdd=sj*dtvol*sj/(.5*rrj*(vol(i,j)+vol(i,j-1)))
        rnu1=ttb*(1.+tt*ttu)*rmu
        rnu2=ttb*(1.+tt)*(rlmbda+2.*rmu)
        rnu3=ttb*(1.+tt*tte)*rkap
        b11=-(bp11*t11+bp12*t21+bp13*t31+bp14*t41)
        b12=-(bp11*t12+bp12*t22+bp13*t32+bp14*t42)
        b13=-(bp11*t13+bp12*t23+bp13*t33+bp14*t43)
        b14=-(bp11*t14+bp12*t24+bp13*t34+bp14*t44)
        b21=-(bp21*t11+bp22*t21+bp23*t31+bp24*t41)
     *      +vdd*( syp*upj*rnu1+sxp*vpj*rnu2)
        b22=-(bp21*t12+bp22*t22+bp23*t32+bp24*t42)
     *      -vdd*(syp*syp*rnu1+sxp*sxp*rnu2)
        b23=-(bp21*t13+bp22*t23+bp23*t33+bp24*t43)
     *      +vdd*sxp*syp*(rnu1-rnu2)
        b24=-(bp21*t14+bp22*t24+bp23*t34+bp24*t44)
        b31=-(bp31*t11+bp32*t21+bp33*t31+bp34*t41)
     *      +vdd*(-sxp*upj*rnu1+syp*vpj*rnu2)
        b32=-(bp31*t12+bp32*t22+bp33*t32+bp34*t42)
     *      +vdd*sxp*syp*(rnu1-rnu2)
        b33=-(bp31*t13+bp32*t23+bp33*t33+bp34*t43)
     *      -vdd*(sxp*sxp*rnu1+syp*syp*rnu2)
        b34=-(bp31*t14+bp32*t24+bp33*t34+bp34*t44)
        b41=-(bp41*t11+bp42*t21+bp43*t31+bp44*t41)
     *      +vdd*(upbar*upj*rnu1+vpbar*vpj*rnu2
     *      -(alpj-ei(i,j))*rnu3)
        b42=-(bp41*t12+bp42*t22+bp43*t32+bp44*t42)
     *      -vdd*( upbar*syp*rnu1+vpbar*sxp*rnu2-uuj*rnu3)
        b43=-(bp41*t13+bp42*t23+bp43*t33+bp44*t43)
     *      -vdd*(-upbar*sxp*rnu1+vpbar*syp*rnu2-vvj*rnu3)
        b44=-(bp41*t14+bp42*t24+bp43*t34+bp44*t44)
     *      -vdd*rnu3
        if(model.eq.0.or.model.eq.1) then
          rnu5=ttb*(1.+tt*ttk)*rmuk
          rnu6=ttb*(1.+tt*tto)*rmuo
          b55=-dd1-vdd*rnu5
          b66=-dd1-vdd*rnu6
          if(model.eq.1) then
            b77=-dd1
            b88=-dd1
            b99=-dd1
          endif
        endif
c=======================================
c        Boundary conditions
c=======================================
        if(j.eq.2) then
          d1=1.5*sj*dt/vol(i,j)
          bp1=alpj*betam*d1
          bp2=(-uuj*betam-sxp*ccj)*d1
          bp3=(-vvj*betam-syp*ccj)*d1
          bp4=betam*d1
          b21=b21+sxp*bp1
          b22=b22+sxp*bp2
          b23=b23+sxp*bp3
          b24=b24+sxp*bp4
          b31=b31+syp*bp1
          b32=b32+syp*bp2
          b33=b33+syp*bp3
          b34=b34+syp*bp4
        endif
c=======================================
c   Compute diagonal block matrix [a]
c=======================================
        a11=a11-b11
        a12=a12-b12
        a13=a13-b13
        a14=a14-b14
        a21=a21-b21
        a22=a22-b22
        a23=a23-b23
        a24=a24-b24
        a31=a31-b31
        a32=a32-b32
        a33=a33-b33
        a34=a34-b34
        a41=a41-b41
        a42=a42-b42
        a43=a43-b43
        a44=a44-b44
        if(model.eq.0.or.model.eq.1) then
          a55=a55-b55
          a66=a66-b66
          if(model.eq.1) then
            a77=a77-b77
            a88=a88-b88
            a99=a99-b99
          endif
        endif
        if(iaxi.eq.1) then
          d1=-1.5*sgn*((x(i+1,j+1)-x(i,j+1))*(y(i+1,j+1)+y(i,j+1))
     *       -(x(i+1,j)-x(i,j))*(y(i+1,j)+y(i,j)))*faxi*dtvol
          a31=a31+d1*(.5*(u(i,j)**2+v(i,j)**2)*betam)
          a32=a32+d1*(-u(i,j)*betam)
          a33=a33+d1*(-v(i,j)*betam)
          a34=a34+d1*betam
          temp=abs(ei(i,j))/cv
          rmu=viscos(temp)
          if(model.lt.0) then
            epj=0.
          elseif(model.eq.2) then
            epj=.5*(ep(i,j)+ep(i,j-1))
          else
            epj=climit(i,j)*ep(i,j)
          endif
          rmu=rmu+epj
          rad=.25*(y(i,j)+y(i+1,j)+y(i,j+1)+y(i+1,j+1))
          a31=a31+d1*2.*rmu*v(i,j)/(rad*rho(i,j,n))
          a33=a33-d1*2.*rmu/(rad*rho(i,j,n))
        endif
        iflag=2
        go to 10
      else
c=======================================
c          [cp] = [S]inv
c=======================================
        cp11=1.
        cp12=0.
        cp13=.5*rccsq
        cp14=.5*rccsq
        cp21=uu
        cp22=syp*rr
        cp23=(uu+sxp*cc)*.5*rccsq
        cp24=(uu-sxp*cc)*.5*rccsq
        cp31=vv
        cp32=-sxp*rr
        cp33=(vv+syp*cc)*.5*rccsq
        cp34=(vv-syp*cc)*.5*rccsq
        cp41=alp
        cp42=rr*up
        cp43=.5*((alp+vp*cc)*rccsq+1./betam)
        cp44=.5*((alp-vp*cc)*rccsq+1./betam)
c=======================================
c  [c] = [S]inv [Lambda] [S] + [M] [N]
c         (Centered at j-1/2)
c=======================================
        vdd=vdd*rho(i,j,n)/rrj
        rnu1=ttc*(1.+tt*ttu)*rmu
        rnu2=ttc*(1.+tt)*(rlmbda+2.*rmu)
        rnu3=ttc*(1.+tt*tte)*rkap
        c11=-(cp11*t11+cp12*t21+cp13*t31+cp14*t41)
        c12=-(cp11*t12+cp12*t22+cp13*t32+cp14*t42)
        c13=-(cp11*t13+cp12*t23+cp13*t33+cp14*t43)
        c14=-(cp11*t14+cp12*t24+cp13*t34+cp14*t44)
        c21=-(cp21*t11+cp22*t21+cp23*t31+cp24*t41)
     *      +vdd*( syp*upj*rnu1+sxp*vpj*rnu2)
        c22=-(cp21*t12+cp22*t22+cp23*t32+cp24*t42)
     *      -vdd*(syp*syp*rnu1+sxp*sxp*rnu2)
        c23=-(cp21*t13+cp22*t23+cp23*t33+cp24*t43)
     *      +vdd*sxp*syp*(rnu1-rnu2)
        c24=-(cp21*t14+cp22*t24+cp23*t34+cp24*t44)
        c31=-(cp31*t11+cp32*t21+cp33*t31+cp34*t41)
     *      +vdd*(-sxp*upj*rnu1+syp*vpj*rnu2)
        c32=-(cp31*t12+cp32*t22+cp33*t32+cp34*t42)
     *      +vdd*sxp*syp*(rnu1-rnu2)
        c33=-(cp31*t13+cp32*t23+cp33*t33+cp34*t43)
     *      -vdd*(sxp*sxp*rnu1+syp*syp*rnu2)
        c34=-(cp31*t14+cp32*t24+cp33*t34+cp34*t44)
        c41=-(cp41*t11+cp42*t21+cp43*t31+cp44*t41)
     *      +vdd*(upbar*upj*rnu1+vpbar*vpj*rnu2
     *      -(alpj-ei(i,j-1))*rnu3)
        c42=-(cp41*t12+cp42*t22+cp43*t32+cp44*t42)
     *      -vdd*( upbar*syp*rnu1+vpbar*sxp*rnu2-uuj*rnu3)
        c43=-(cp41*t13+cp42*t23+cp43*t33+cp44*t43)
     *      -vdd*(-upbar*sxp*rnu1+vpbar*syp*rnu2-vvj*rnu3)
        c44=-(cp41*t14+cp42*t24+cp43*t34+cp44*t44)
     *      -vdd*rnu3
        if(model.eq.0.or.model.eq.1) then
          rnu5=ttc*(1.+tt*ttk)*rmuk
          rnu6=ttc*(1.+tt*tto)*rmuo
          c55=-dd1-vdd*rnu5
          c66=-dd1-vdd*rnu6
          if(model.eq.1) then
            c77=-dd1-vdd*rnu5
            c88=-dd1-vdd*rnu5
            c99=-dd1-vdd*rnu5
          endif
        endif
      endif
 20   return
      end
      subroutine l(dt)
c=======================================================================
c             U P D A T E   F L O W   V A R I A B L E S
c=======================================================================
      include 'come2c.'
c=======================================================================
      str(k1)=0.
      str(k2)=0.
      iadd=mod(nadv,2)
      jadd=mod(((iadd+nadv)/2),2)
c-----------------------------------------------------------------------
c                 X - T R A N S P O R T   T E R M S
c-----------------------------------------------------------------------
      do 100 n=1,2
        m=3-n
        nm1=n-1
        do 20 j=2,jlm1
          call fsi(1,j,1+iadd)
          do 10 i=2,ilm1
            k=k1
            k1=k2
            k2=k
            call fsi(i,j,i+iadd)
            dtv=dt/vol(i,j)
            du(i,j,1)=-dtv*(fs1(k2)-fs1(k1))
            du(i,j,2)=-dtv*(fs2(k2)-fs2(k1))
            du(i,j,3)=-dtv*(fs3(k2)-fs3(k1))
            du(i,j,4)=-dtv*(fs4(k2)-fs4(k1))
            if(model.eq.0.or.model.eq.1) then
              du(i,j,5)=-dtv*(fs5(k2)-fs5(k1))
              du(i,j,6)=-dtv*(fs6(k2)-fs6(k1))
              if(model.eq.1) then
                du(i,j,7)=-dtv*(fs7(k2)-fs7(k1))
                du(i,j,8)=-dtv*(fs8(k2)-fs8(k1))
                du(i,j,9)=-dtv*(fs9(k2)-fs9(k1))
              endif
            endif
c===================
c    Axisymmetric
c===================
            if(iaxi.eq.1) then
              ss=(-(x(i+1,j+1)-x(i+1,j))*(y(i+1,j+1)+y(i+1,j))
     *            +(x(i,j+1)-x(i,j))*(y(i,j+1)+y(i,j)))*faxi
              rad=.25*(y(i,j)+y(i+1,j)+y(i,j+1)+y(i+1,j+1))
              temp=ei(i,j)/cv
              rmu=viscos(temp)
              if(model.lt.0) then
                if(model.eq.2) then
                  rmu=rmu+.5*(ep(i,j)+ep(i-1,j))
                else
                  rmu=rmu+fmodel*climit(i,j)*ep(i,j)
                endif
              endif
              sigtt=-p(i,j)+2.*rmu*v(i,j)/rad-.5*(str(k1)+str(k2))
              du(i,j,3)=du(i,j,3)-dtv*ss*sigtt
            endif
 10       continue
 20     continue
c-----------------------------------------------------------------------
c                 Y - T R A N S P O R T   T E R M S
c-----------------------------------------------------------------------
        do 40 i=2,ilm1
          call fsj(i,1,1+jadd,dt)
          do 30 j=2,jlm1
            k=k1
            k1=k2
            k2=k
            call fsj(i,j,j+jadd,dt)
            dtv=dt/vol(i,j)
            du(i,j,1)=du(i,j,1)-dtv*(fs1(k2)-fs1(k1))
            du(i,j,2)=du(i,j,2)-dtv*(fs2(k2)-fs2(k1))
            du(i,j,3)=du(i,j,3)-dtv*(fs3(k2)-fs3(k1))
            du(i,j,4)=du(i,j,4)-dtv*(fs4(k2)-fs4(k1))
            if(model.eq.0.or.model.eq.1) then
              call prodk(i,j,dt)
              du(i,j,5)=du(i,j,5)-dtv*(fs5(k2)-fs5(k1))+srk
              du(i,j,6)=du(i,j,6)-dtv*(fs6(k2)-fs6(k1))+sro
              if(model.eq.1) then
                du(i,j,7)=du(i,j,7)-dtv*(fs7(k2)-fs7(k1))+srtxx
                du(i,j,8)=du(i,j,8)-dtv*(fs8(k2)-fs8(k1))+srtxy
                du(i,j,9)=du(i,j,9)-dtv*(fs9(k2)-fs9(k1))+srtyy
              endif
            endif
c===================
c    Axisymmetric
c===================
            if(iaxi.eq.1) then
              ss=(+(x(i+1,j+1)-x(i,j+1))*(y(i+1,j+1)+y(i,j+1))
     *            -(x(i+1,j)-x(i,j))*(y(i+1,j)+y(i,j)))*faxi
              rad=.25*(y(i,j)+y(i+1,j)+y(i,j+1)+y(i+1,j+1))
              temp=ei(i,j)/cv
              rmu=viscos(temp)
              if(model.lt.0) then
                if(model.eq.2) then
                  rmu=rmu+.5*(ep(i,j)+ep(i,j-1))
                else
                  rmu=rmu+fmodel*climit(i,j)*ep(i,j)
                endif
              endif
              sigtt=-p(i,j)+2.*rmu*v(i,j)/rad-.5*(str(k1)+str(k2))
              du(i,j,3)=du(i,j,3)-dtv*ss*sigtt
            endif
 30       continue
 40     continue
c-----------------------------------------------------------------------
c               C O M P U T E   R E S I D U A L   S U M
c-----------------------------------------------------------------------
        resid=0.
        do 60 i=2,ilm1
          do 50 j=2,jlm1
            resid=resid+du(i,j,1)**2*vol(i,j)
  50      continue
  60    continue
        resid=resid/((il-2)*(jl-2)*rhoinf*dt)
c-----------------------------------------------------------------------
c                  C O M P U T E   I N C R E M E N T S
c-----------------------------------------------------------------------
        if(implt.gt.0) call limp(dt,iadd,jadd)
c-----------------------------------------------------------------------
c              U P D A T E   F L O W   P R O P E R T I E S
c-----------------------------------------------------------------------
        rn=n
        rnm1=rn-1.
        irealz=0
        do 90 i=2,ilm1
          do 70 jx=2,jlm1
            j=jl+1-jx
            rho(i,j,m)=(rnm1*rho(i,j,m)+rho(i,j,n)+du(i,j,1))/rn
            rhou(i,j,m)=(rnm1*rhou(i,j,m)+rhou(i,j,n)+du(i,j,2))/rn
            rhov(i,j,m)=(rnm1*rhov(i,j,m)+rhov(i,j,n)+du(i,j,3))/rn
            e(i,j,m)=(rnm1*e(i,j,m)+e(i,j,n)+du(i,j,4))/rn
            if(model.eq.0.or.model.eq.1) then
              rhok(i,j,m)=(rnm1*rhok(i,j,m)+rhok(i,j,n)+du(i,j,5))/rn
              rhomeg(i,j,m)=(rnm1*rhomeg(i,j,m)+rhomeg(i,j,n)
     *                     +du(i,j,6))/rn
c==================================
c   Check for negative turbulent
c    energy or dissipation rate
c==================================
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              if(rhok(i,j,m).le.rhokmn) then
                rhok(i,j,m)=.5*abs((rhok(i,j-1,n)+rhok(i,j+1,n)))
                irealz=irealz+1
              endif
              if(rhomeg(i,j,m).le.rhowmn) then
                rhomeg(i,j,m)=.5*abs((rhomeg(i,j-1,n)+rhomeg(i,j+1,n)))
                irealz=irealz+1
              endif
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              tke(i,j)=rhok(i,j,m)/rho(i,j,m)
              omeg(i,j)=rhomeg(i,j,m)/rho(i,j,m)
              if(model.eq.1) then
                rhotxx(i,j,m)=(rnm1*rhotxx(i,j,m)+rhotxx(i,j,n)
     *                       +du(i,j,7))/rn
                rhotxy(i,j,m)=(rnm1*rhotxy(i,j,m)+rhotxy(i,j,n)
     *                       +du(i,j,8))/rn
                rhotyy(i,j,m)=(rnm1*rhotyy(i,j,m)+rhotyy(i,j,n)
     *                       +du(i,j,9))/rn
              endif
            endif
            u(i,j)=rhou(i,j,m)/rho(i,j,m)
            v(i,j)=rhov(i,j,m)/rho(i,j,m)
            ei(i,j)=e(i,j,m)/rho(i,j,m)-tke(i,j)
     *             -.5*(u(i,j)*u(i,j)+v(i,j)*v(i,j))
            if(ei(i,j).lt.0.) then
              write(*,*) 'Negative energy at i,j = ',i,j
              ei(i,j)=abs(ei(i,j))
            endif
            p(i,j)=gm1*rho(i,j,m)*ei(i,j)
 70       continue
          if(model.eq.1) then
c==================================
c    Check u'u', v'v' and u'v'
c        for realizability
c==================================
            do 80 j=2,jlm1
              if(abs(rhotxx(i,j,m)).lt.rhokmn) then
                rhotxx(i,j,m)=.5*(rhotxx(i,j-1,n)+rhotxx(i,j+1,n))
              endif
              if(abs(rhotxy(i,j,m)).lt.rhokmn) then
                rhotxy(i,j,m)=.5*(rhotxy(i,j-1,n)+rhotxy(i,j+1,n))
              endif
              if(abs(rhotyy(i,j,m)).lt.rhokmn) then
                rhotyy(i,j,m)=.5*(rhotyy(i,j-1,n)+rhotyy(i,j+1,n))
              endif
              boundu=1.6*rhok(i,j,m)
              sumrhk=-(rhotxx(i,j,m)+rhotyy(i,j,m))
              if(sumrhk.gt.boundu) then
                rhotxx(i,j,m)=rhotxx(i,j,m)*boundu/sumrhk
                rhotyy(i,j,m)=rhotyy(i,j,m)*boundu/sumrhk
                rhotmx=.3*rhok(i,j,m)
                if(abs(rhotxy(i,j,m)).gt.rhotmx) then
                  if(rhotxy(i,j,m).lt.0.) then
                    rhotxy(i,j,m)=-rhotmx
                  else
                    rhotxy(i,j,m)=rhotmx 
                  endif
                endif
c====================================
c  Realizability violation checking
c====================================
c                irealz=irealz+1
c====================================
              endif
              txx(i,j)=rhotxx(i,j,m)/rho(i,j,m)
              txy(i,j)=rhotxy(i,j,m)/rho(i,j,m)
              tyy(i,j)=rhotyy(i,j,m)/rho(i,j,m)
 80         continue
          endif
 90     continue
c==================================
c     Set boundary conditions
c==================================
        call bc
        if(model.eq.0.or.model.eq.1) then
          call vgrads
          call epskw
          call crdiff
          call strech
          if(model.eq.1.and.ntau.ne.0) call dtau
        endif
        iadd=mod((iadd+1),2)
        jadd=mod((jadd+1),2)
 100  continue
      write(iunit1,110) nadv,time,cfl,resid
      return
c=======================================================================
 110  format(3x,'Nadv=',i5,3x,'Time=',1pe12.4,3x,'cfl=',e12.4,
     *                                   3x,'Residual=',e12.4)
c=======================================================================
      end
      subroutine limp(dt,iadd,jadd)
c=======================================================================
c    C O O R D I N A T E   G A U S S - S E I D E L   S O L U T I O N
c                      ni = 1...Backward sweep
c                      ni = 2...Forward  sweep
c=======================================================================
      include 'come2c.'
c=======================================================================
c  Set number of equations to be solved
c=======================================
      if(model.eq.1) then
        lmax=9
      elseif(model.eq.0) then
        lmax=6
      else
        lmax=4
      endif
c=======================================
c  Save a copy of the explicit {du} to
c  impose  b.c.'s on the implicit {du}
c=======================================
      do 30 l=1,lmax
        do 20 i=1,il
          do 10 j=1,jl
            if(i.eq.il.or.j.eq.jl) du(i,j,l)=0.
            dus(i,j,l)=du(i,j,l)
 10       continue
 20     continue
 30   continue
      do 70 ni=1,2
        do 60 ix=2,ilm1
c=======================================
c     Prepare for backward sweep
c=======================================
          if(ni.eq.1) then
            i=ilm1+2-ix
c=======================================
c    Set b.c.'s for forward sweep
c=======================================
          else
            i=ix
            do 50 l=1,lmax
              do 40 j=2,jl
                if(i.eq.2) then
                  du(il,j,l)=du(ilm1,j,l)
                  if(icase.gt.0) then
                    du(1,j,l)=0.
                  else
                    du(1,j,l)=du(2,j,l)
                    if(l.eq.3) du(1,j,l)=-du(1,j,l)
                  endif
                endif
                du(i,j,l)=dus(i,j,l)
 40           continue
 50         continue
          endif
c=======================================
c    Do block-tridiagonal inversion
c=======================================
          call btdj(dt,i,iadd,jadd)
 60     continue
 70   continue
      return
      end
      subroutine namin(icpu)
c=======================================================================
c           COORDINATE READING INPUT DATA FILE 'EDDY2C.DAT'
c=======================================================================
      parameter (niu=6,nn=22,nxc=6,nxcm=6,nxd=24)
      character a*12,newfil*7,pform*9
      common/close/ xc(nxc)
      common/closem/xcm(nxcm)
      common/data/  xd(nxd)
      common/iodef/ iu(niu)
      common/ndata/ n(nn)
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
      open(io,file='eddy2c.prt',form=pform,status=newfil)
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
      write(io,400)
c----------------------------------
c       Open disk files
c----------------------------------
      open(iu(3),file='input2c.dat',form='unformatted',status='old')
      open(iu(4),file='output2c.dat',form='unformatted',status=newfil)
      open(iu(6),file='wall2c.dat',status=newfil)
      return
c=======================================================================
 100  format(1x,a12,i6)
 200  format(1x,'Program EDDY2C input data:'/,
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
        stop
      endif
      return
      end
      subroutine plotme
c=======================================================================
c        S A V E   S K I N   F R I C T I O N   A N D   W A L L
c              P R E S S U R E   F O R   P L O T T I N G
c=======================================================================
      include 'come2c.'
c=======================================================================
      write(iuplot,100) ilm1-1,ipwall
      if(inletp.eq.0) then
        pref=pinf
      else
        pref=p(3,jlm1)
      endif
      do 10 i=2,ilm1
        xin=.5*zscale*(x(i,2)+x(i+1,2))
        if(ipwall.eq.0) then
          pwop0=.5*(p(i,1)+p(i,2))/pinf
        else
          pwop0=(.5*(p(i,1)+p(i,2))-pref)/dynp
        endif
        if(i.lt.ile) then
          cf=0.
        else
          sx=x(i+1,2)-x(i,2)
          sy=y(i+1,2)-y(i,2)
          sj=sqrt(sx*sx+sy*sy)
          sxp=sx/sj
          syp=sy/sj
          dxe=x(i+1,2)-x(i,2)
          dye=y(i+1,2)-y(i,2)
          dxn=.25*(x(i+1,3)+x(i,3)-x(i+1,1)-x(i,1))
          dyn=.25*(y(i+1,3)+y(i,3)-y(i+1,1)-y(i,1))
          dxy=1./(dxe*dyn-dxn*dye)
          dun=sxp*(u(i,2)-u(i,1))+syp*(v(i,2)-v(i,1))
          dudn=(sxp*dxe+syp*dye)*dun*dxy
          temp=abs(.5*(ei(i,1)+ei(i,2))/cv)
          rmu=viscos(temp)
          cf=rmu*dudn/dynp
        endif
        write(iuplot,110) xin,pwop0,cf
 10   continue
      return
c=======================================================================
 100  format(2i6)
 110  format(1p3e12.4)
c=======================================================================
      end
      subroutine prntff
c=======================================================================
c                     P R I N T   P R O F I L E S
c=======================================================================
      include 'come2c.'
c=======================================================================
c          Flowfield arrays
c=======================================
      uinf2=uinf**2
      oscale=rmuinf/(rhouin*uinf)
      tscale=1./(rhoinf*uinf2)
      write(iunit1,100) time,nadv-1
      do 20 i=1,il,ipro
        xbar=zscale*x(i,2)
        write(iunit1,110) i,xbar
        if(model.eq.-1.or.model.eq.2) then
          write(iunit1,120)
        elseif(model.eq.0) then
          write(iunit1,130)
        else
          write(iunit1,140)
        endif
        do 10 j=1,jl
          uoui=u(i,j)/uinf
          voui=v(i,j)/uinf
          rhorhi=rho(i,j,1)/rhoinf
          popi=p(i,j)/pinf
          toti=popi/rhorhi
          temp=abs(ei(i,j)/cv)
          rmu=viscos(temp)
          ret=ep(i,j)/rmu
          if(model.eq.-1.or.model.eq.2) then
            write(iunit1,150) j,uoui,voui,rhorhi,popi,toti,ret
          else
            tkbar=tke(i,j)/uinf2
            ombar=oscale*omeg(i,j)
            rhtxy=tscale*climit(i,j)*ep(i,j)*(dudy(i,j)+dvdx(i,j))
            if(model.eq.0) then
              write(iunit1,160) j,uoui,voui,rhorhi,popi,toti,ret,
     *                          tkbar,omeg(i,j),rhtxy
            else
              term=.6666667*rhok(i,j,1)
              rhtxx=tscale*(rhotxx(i,j,1)-term)
              rhtxy=tscale*rhotxy(i,j,1)
              rhtyy=tscale*(rhotyy(i,j,1)-term)
              write(iunit1,170) j,uoui,voui,rhorhi,popi,toti,ret,
     *                          tkbar,ombar,rhtxx,rhtxy,rhtyy
            endif
          endif
 10     continue
 20   continue
      return
c=======================================================================
 100  format(21x,'Time = ',1pe14.7,2x,'Nadv = ',i9)
 110  format(/8x,'Column  ',i3,5x,'x = ',1pe12.3)
 120  format(3x,'j',7x,'u/U',8x,'v/U',8x,'rho',9x,'p',10x,'T',9x,'Ret')
 130  format(3x,'j',7x,'u/U',8x,'v/U',8x,'rho',9x,'p',10x,'T',9x,'Ret',
     *       7x,'k/U**2',5x,'omega',5x,'rhotxy')
 140  format(3x,'j',6x,'u/U',7x,'v/U',7x,'rho',8x,'p',9x,'T',8x,'Ret',
     *       6x,'k/U**2',4x,'omega',5x,'rhotxx',4x,'rhotxy',4x,'rhotyy')
 150  format(1x,i3,3x,1p6e11.4)
 160  format(1x,i3,3x,1p9e11.4)
 170  format(1x,i3,3x,1p11e10.3)
c=======================================================================
      end
      subroutine prodk(i,j,dt)
c=======================================================================
c      T U R B U L E N C E   M O D E L   S O U R C E   T E R M S
c=======================================================================
      include 'come2c.'
c=======================================================================
      divu=dudx(i,j)+dvdy(i,j)
      emt2=2.*rhok(i,j,m)/(gam*p(i,j))-xmtsq
      if(emt2.le.0.) emt2=0.
c-----------------------------------------------------------------------
c                      K - O M E G A   M O D E L
c-----------------------------------------------------------------------
c=======================================
c     Left-hand-side source terms
c=======================================
      if(model.eq.0) then
        pk1=(dudy(i,j)+dvdx(i,j))**2+2.*(dudx(i,j)**2+dvdy(i,j)**2)
        if(iaxi.eq.1) then
          divu=divu+vor(i,j)
          pk1=pk1+2.*vor(i,j)**2
        endif
        pk1=pk1-.6666667*divu**2
        pk(i,j)=pk1*falfs(i,j)/omeg(i,j)
        climit(i,j)=1.
        pod=abs(clim2*pk(i,j))/(betas*omeg(i,j))
        if(pod.gt.1) then
          climit(i,j)=1./sqrt(pod)
        endif
        pk(i,j)=climit(i,j)*pk(i,j)-.6666667*divu
c=======================================
c     Right-hand-side source terms
c=======================================
        srk=(pk(i,j)-betas*fbet(i,j)*(1.+xis*emt2)*omeg(i,j))
     *     *rhok(i,j,m)*dt
        sro=(alpha*falf(i,j)*pk(i,j)
     *     -(beta*fbeta(i,j)-betas*fbet(i,j)*xis*emt2)*omeg(i,j))
     *     *rhomeg(i,j,m)*dt
c------------------------
c  Add cross diffusion
c------------------------
        sro=sro+cross(i,j)*rhomeg(i,j,m)*dt
        pk(i,j)=abs(pk(i,j))
c-----------------------------------------------------------------------
c                 S T R E S S - O M E G A   M O D E L
c-----------------------------------------------------------------------
      else
        climit(i,j)=1.
c------------------------
c   P-production tensor
c------------------------
        pxx=2.*(dudx(i,j)*rhotxx(i,j,m)+dudy(i,j)*rhotxy(i,j,m))
        pxy=dvdx(i,j)*rhotxx(i,j,m)+dudy(i,j)*rhotyy(i,j,m)
     *          +divu*rhotxy(i,j,m)
        pyy=2.*(dvdx(i,j)*rhotxy(i,j,m)+dvdy(i,j)*rhotyy(i,j,m))
c------------------------
c   D-production tensor
c------------------------
        dxx=2.*(dudx(i,j)*rhotxx(i,j,m)+dvdx(i,j)*rhotxy(i,j,m))
        dxy=dudy(i,j)*rhotxx(i,j,m)+dvdx(i,j)*rhotyy(i,j,m)
     *          +divu*rhotxy(i,j,m)
        dyy=2.*(dudy(i,j)*rhotxy(i,j,m)+dvdy(i,j)*rhotyy(i,j,m))
c------------------------
c    k-Production term
c------------------------
        prod=.5*(pxx+pyy)
c------------------------
c      Axisymmetric
c------------------------
        if(iaxi.eq.1) then
          divu=divu+vor(i,j)
          prod=prod-vor(i,j)*(rhotxx(i,j,m)+rhotyy(i,j,m))
        endif
        divu33=divu/3.
c------------------------
c   Strain-rate tensor
c------------------------
        sxx=dudx(i,j)-divu33
        sxy=.5*(dudy(i,j)+dvdx(i,j))
        syy=dvdy(i,j)-divu33
c------------------------
c   Dissipation term
c------------------------
        diss=betas*fbet(i,j)*(1.+xis*emt2)*fcc1(i,j)*cc1*omeg(i,j)
c=======================================
c     Right-hand-side source terms
c=======================================
        srtxx=(falbet(i,j)*(-csw1*pxx+betah*dxx+csw2*prod)
     *        -diss*rhotxx(i,j,m)
     *       +(csw3*falbet(i,j)+gammah*fgamh(i,j))*sxx*rhok(i,j,m))*dt
        srtxy=(falbet(i,j)*(-csw1*pxy+betah*dxy)
     *        -diss*rhotxy(i,j,m)
     *       +(csw3*falbet(i,j)+gammah*fgamh(i,j))*sxy*rhok(i,j,m))*dt
        srtyy=(falbet(i,j)*(-csw1*pyy+betah*dyy+csw2*prod)
     *        -diss*rhotyy(i,j,m)
     *       +(csw3*falbet(i,j)+gammah*fgamh(i,j))*syy*rhok(i,j,m))*dt
        prod=prod-.66666667*divu*rhok(i,j,m)
        pk(i,j)=prod/rhok(i,j,m)
        srk=(prod-betas*fbet(i,j)*(1.+xis*emt2)*omeg(i,j)
     *     *rhok(i,j,m))*dt
        sro=(alpha*falf(i,j)*pk(i,j)
     *     -(beta*fbeta(i,j)-betas*fbet(i,j)*xis*emt2)*omeg(i,j))
     *     *rhomeg(i,j,m)*dt
c------------------------
c  Add cross diffusion
c------------------------
        sro=sro+cross(i,j)*rhomeg(i,j,m)*dt
c=======================================
c   Extra terms for axisymmetric flow
c=======================================
        if(iaxi.eq.1) then
          temp=abs(ei(i,j)/cv)
          rmu=viscos(temp)
          rmuk=rmu+sigmas*ep(i,j)
          tzz=-txx(i,j)-tyy(i,j)
          radius=.25*(y(i,j+1)+y(i,j)+y(i+1,j+1)+y(i+1,j))
          srtxy=srtxy-dt*rmuk*txy(i,j)/radius**2
          srtyy=srtyy-2.*dt*rmuk*(tyy(i,j)-tzz)/radius**2
        endif
c=======================================
c     Left-hand-side source terms
c=======================================
        pk(i,j)=abs(pk(i,j))
        ptij(i,j)=4.*diss
      endif
      return
      end
      subroutine prtwal
c=======================================================================
c          P R I N T   S U R F A C E   P R O P E R T I E S
c=======================================================================
      include 'come2c.'
c=======================================================================
c     Display time and step number
c=======================================
      write(iunit1,100) time,nadv,cfl
      write(*,100) time,nadv,cfl
      if(ipwall.eq.0) then
        write(iunit1,110)
        write(*,110)
      else
        write(iunit1,120)
        write(*,120)
      endif
      isep=0
      if(inletp.eq.0) then
        pref=pinf
      else
        pref=p(3,jlm1)
      endif
      do 10 i=1,il
        xin=.5*zscale*(x(i,2)+x(i+1,2))
        pwall=.5*(p(i,1)+p(i,2))
        if(ipwall.eq.0) then
          pwop0=pwall/pinf
        else
          pwop0=(pwall-pref)/dynp
        endif
        if(i.lt.ile) then
          cf=0.
          q=0.
          yplus=0.
        else
          sx=x(i+1,2)-x(i,2)
          sy=y(i+1,2)-y(i,2)
          sj=sqrt(sx*sx+sy*sy)
          sxp=sx/sj
          syp=sy/sj
          dxe=x(i+1,2)-x(i,2)
          dye=y(i+1,2)-y(i,2)
          dxn=.25*(x(i+1,3)+x(i,3)-x(i+1,1)-x(i,1))
          dyn=.25*(y(i+1,3)+y(i,3)-y(i+1,1)-y(i,1))
          dxy=1./(dxe*dyn-dxn*dye)
          dun=sxp*(u(i,2)-u(i,1))+syp*(v(i,2)-v(i,1))
          dudn=(sxp*dxe+syp*dye)*dun*dxy
          temp=abs(.5*(ei(i,1)+ei(i,2))/cv)
          rmu=viscos(temp)
          cf=rmu*dudn/dynp
          rhow=pwall/(rgas*temp)
          dyp=sxp*dyn-syp*dxn
          uplus=sqrt(abs(cf*dynp/rhow))
          yplus=.5*rhow*uplus*dyp/rmu
          if(iadbwl.eq.0) then
            dee=ei(i+1,2)-ei(i,2)
            den=.5*(ei(i,2)-ei(i,1))
            dex=(dee*dyn-den*dye)*dxy
            dey=(den*dxe-dee*dxn)*dxy
            dedn=sxp*dey-syp*dex
            q=gmoprl*rmu*dedn/qdm
          else
            q=0.
          endif
        endif
        write(iunit1,130) i,xin,pwop0,cf,q,yplus
        if(mod(i,iskip).eq.0) write(*,130) i,xin,pwop0,cf,q,yplus
        if(i.eq.ile) then
          cfsav=cf
          xsav=xin
          xsep=xin
          xret=xin
        elseif(i.gt.ile) then
          if(cf*cfsav.lt.0.) then
            if(cf.lt.0.) then
              xsep=xsav+cfsav*(xin-xsav)/(cfsav-cf)
              isep=1
            else
              xret=xsav+cfsav*(xin-xsav)/(cfsav-cf)
            endif
          endif
          cfsav=cf
          xsav=xin
        endif
 10   continue
      if(isep.ne.0) then
        xbub=xret-xsep
        write(iunit1,140) xsep,xret,xbub
        write(*,140) xsep,xret,xbub
      endif
      return
c=======================================================================
 100  format(//4x,'Time = ',1pe13.6,5x,'Nadv = ',i6,5x,'Cfl = ',e10.2)
 110  format(/4x,'i',5x,'x',9x,'P/Po',9x,'Cf',11x,'q',9x,'yplus')
 120  format(/4x,'i',5x,'x',10x,'Cp',10x,'Cf',11x,'q',9x,'yplus')
 130  format(i5,1p5e12.5)
 140  format(/1x,'Separation bubble data: Xs =',1pe10.3,'  Xr =',e10.3,
     *           '  Xr-Xs =',e10.3/)
c=======================================================================
      end
      subroutine readit
c=======================================================================
c     R E A D   F L O W F I E L D   A R R A Y S   F R O M   D I S K
c=======================================================================
      include 'come2c.'
c=======================================================================
c              M E A N - F L O W   P R O P E R T I E S
c----------------------------------------------------------------------
c      Read time, indices, etc.
c=======================================
      read(iunit3) nadv,time,ilold,jlold,cfl,nvsold,iadold,eiwold
c=======================================
c   Compare with current input data
c=======================================
      if(ilold.ne.il.or.jlold.ne.jl) then
        write(iunit1,200) il,jl,ilold,jlold
        write(*,200) il,jl,ilold,jlold
        stop
      elseif(nvsold.ne.ntau.or.iadold.ne.iadbwl) then
        write(iunit1,210) nvsold,iadold,ntau,iadbwl
        write(*,210) nvsold,iadold,ntau,iadbwl
      endif
c=======================================
c       Read mean-flow arrays
c=======================================
      read(iunit3) ((rho(i,j,1),u(i,j),v(i,j),ei(i,j),p(i,j),
     *             ep(i,j),vol(i,j),dvol(i,j),i=1,il),j=1,jl)
c=======================================
c        Read mesh coordinates
c=======================================
      ilp1=il+1
      jlp1=jl+1
      read(iunit3) ((x(i,j),y(i,j),i=1,ilp1),j=1,jlp1)
      read(iunit3) (eta(j),j=2,jlp1)
c=======================================
c     Initialize flowfield arrays
c=======================================
      do 30 j=1,jl
        do 20 i=1,il
          rhou(i,j,1)=rho(i,j,1)*u(i,j)
          rhov(i,j,1)=rho(i,j,1)*v(i,j)
          e(i,j,1)=rho(i,j,1)*(ei(i,j)+.5*(u(i,j)*u(i,j)
     *                                    +v(i,j)*v(i,j)))
          rho(i,j,2)=rho(i,j,1)
          rhou(i,j,2)=rhou(i,j,1)
          rhov(i,j,2)=rhov(i,j,1)
          e(i,j,2)=e(i,j,1)
          ut(i,j)=0.
          vt(i,j)=0.
          tke(i,j)=0.
          climit(i,j)=1.
          do 10 k=1,9
            du(i,j,k)=0.
            dus(i,j,k)=0.
 10       continue
 20     continue
 30   continue
c-----------------------------------------------------------------------
c       T U R B U L E N C E - F I E L D   P R O P E R T I E S
c-----------------------------------------------------------------------
      if(model.eq.0.or.model.eq.1) then
c=======================================
c       Read k and omega arrays
c=======================================
        read(iunit3) ((tke(i,j),omeg(i,j),i=1,il),j=1,jl)
        do 50 j=1,jl
          do 40 i=1,il
            rhok(i,j,1)=rho(i,j,1)*tke(i,j)
            rhomeg(i,j,1)=rho(i,j,1)*omeg(i,j)
            rhok(i,j,2)=rhok(i,j,1)
            rhomeg(i,j,2)=rhomeg(i,j,1)
            e(i,j,1)=e(i,j,1)+rhok(i,j,1)
            e(i,j,2)=e(i,j,1)
 40       continue
 50     continue
c=======================================
c     Read Reynolds-stress arrays
c=======================================
        if(model.eq.1) then
          read(iunit3) ((rhotxx(i,j,1),rhotxy(i,j,1),rhotyy(i,j,1),
     *                                 rhokme(i,j,1),i=1,il),j=1,jl)
          do 70 j=1,jl
            do 60 i=1,il
              rhotxx(i,j,2)=rhotxx(i,j,1)
              rhotxy(i,j,2)=rhotxy(i,j,1)
              rhotyy(i,j,2)=rhotyy(i,j,1)
              txx(i,j)=rhotxx(i,j,1)/rho(i,j,1)
              txy(i,j)=rhotxy(i,j,1)/rho(i,j,1)
              tyy(i,j)=rhotyy(i,j,1)/rho(i,j,1)
              dtxx(i,j)=0.
              dtxy(i,j)=0.
              dtyy(i,j)=0.
 60         continue
 70       continue
        endif
c=======================================
c    Initialize viscous-mod arrays
c=======================================
        do 90 i=1,il
          do 80 j=1,jl
            ep(i,j)=abs(rhok(i,j,1)/omeg(i,j))
 80       continue
 90     continue
        if(model.eq.0) then
          if(nvisc.ne.0) then
            call viskw
          else
            do 110 i=1,il
              do 100 j=1,jl
                falfs(i,j)=1.
                falf(i,j)=1.
                fbet(i,j)=1.
 100          continue
 110        continue
          endif
        elseif(model.eq.1) then
          if(nvisc.ne.0) then
            call vissw
          else
            do 130 i=1,il
              do 120 j=1,jl
                falfs(i,j)=1.
                falf(i,j)=1.
                fbet(i,j)=1.
                falbet(i,j)=1.
                fgamh(i,j)=1.
                fcc1(i,j)=1.
 120          continue
 130        continue
          endif
        endif
      endif
      close(iunit3)
      return
c=======================================================================
 200  format(//1x,'-------INPUT GRID SIZE ERROR-------',
     *        /1x,'New grid size: il = ',i4,', jl = ',i4
     *        /1x,'Old grid size: il = ',i4,', jl = ',i4
     *        /1x,'----------RUN TERMINATED-----------')
 210  format(//1x,'------------INPUT DATA CHANGES-----------',
     *        /1x,'Saved on disk: ntau = ',i4,', iadbwl = ',i4,
     *        /1x,'Input values:  ntau = ',i4,', iadbwl = ',i4,
     *        /1x,'-----------------------------------------'///)
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
      if(model.eq.0.or.model.eq.1) then
        write(iunit4) ((tke(i,j),omeg(i,j),i=1,il),j=1,jl)
        if(model.eq.1)
     *    write(iunit4) ((rhotxx(i,j,1),rhotxy(i,j,1),rhotyy(i,j,1),
     *                                  rhokme(i,j,1),i=1,il),j=1,jl)
      endif
      close(iunit4)
      return
      end
      subroutine strech
c=======================================================================
c         V O R T E X - S T R E T C H I N G   F U N C T I O N
c=======================================================================
      include 'come2c.'
c=======================================================================
      do 20 i=2,ilm1
        do 10 j=2,jlm1
          if(iaxi.eq.0) then
            fbeta(i,j)=1.
          else
            omegxy=.5*(dvdx(i,j)-dudy(i,j))
            chiw=abs(vor(i,j)*omegxy*omegxy/(betas*omeg(i,j)**3))
            fbeta(i,j)=(1.+85.*chiw)/(1.+100.*chiw)
          endif
 10     continue
 20   continue
      return
      end
      subroutine tmestp(dt)
c=======================================================================
c                           T I M E S T E P
c=======================================================================
      include 'come2c.'
c=======================================================================
      if(nadv.le.ncfl) cfl=cflmul*cfl
      do 20 i=2,il
        do 10 j=2,jl
          si=sqrt((x(i,j+1)-x(i,j))**2+(y(i,j+1)-y(i,j))**2)
          sxp=(y(i,j+1)-y(i,j))/si
          syp=(x(i,j)-x(i,j+1))/si
          uu=sxp*u(i,j)+syp*v(i,j)
          sj=sqrt((x(i+1,j)-x(i,j))**2+(y(i+1,j)-y(i,j))**2)
          sxp=(y(i,j)-y(i+1,j))/sj
          syp=(x(i+1,j)-x(i,j))/sj
          vv=sxp*u(i,j)+syp*v(i,j)
          cc=sqrt(ggm1*abs(ei(i,j)))
          if(iaxi.eq.1) then
            factor=.5*(y(i,j)+y(i+1,j)+y(i,j+1)+y(i+1,j+1))*faxi
            si=si*factor
            sj=sj*factor
          endif
          dti=.5*vol(i,j)/(si*(abs(uu)+cc))
          dtj=.5*vol(i,j)/(sj*(abs(vv)+cc))
          if(i.eq.2.and.j.eq.2) then
            dte=dti
            dt=cfl*dti
          elseif(dti.lt.dte) then
            dte=dti
          endif
          if(dtj.lt.dte) dte=dtj
 10     continue
 20   continue
      time=time+dt
      return
      end
      subroutine turbfl
c=======================================================================
c               B A L D W I N - L O M A X   M O D E L
c=======================================================================
      include 'come2c.'
      dimension yi(jmx),aomega(jmx)
c=======================================================================
      aplus=26.
      ccp=1.6
      ckleb=.3
      cwk=1.
      smlk=.4
      bigk=.0168
      prt=.9
c=======================================
c     Compute umax and (du/dy)max
c=======================================
      do 40 i=1,il
        sj=sqrt((x(i+1,2)-x(i,2))**2+(y(i+1,2)-y(i,2))**2)
        sxp=(y(i,2)-y(i+1,2))/sj
        syp=(x(i+1,2)-x(i,2))/sj
        if(iaxi.eq.1) sj=sj*(y(i+1,2)+y(i,2))*faxi
        dudymx=2.*abs(syp*u(i,2)-sxp*v(i,2))*sj/vol(i,2)
        u2max=u(i,2)*u(i,2)+v(i,2)*v(i,2)
        do 10 j=2,jlfm
          u2=u(i,j+1)*u(i,j+1)+v(i,j+1)*v(i,j+1)
          if(u2.gt.u2max) u2max=u2
 10     continue
c=======================================
c   Compute vorticity, fmax and ymax
c=======================================
        temp=.5*(ei(i,1)+ei(i,2))/cv
        rmu=viscos(temp)
        ra=sqrt(rho(i,2,1)*rmu*dudymx)/(aplus*rmu)
        fmax=0.
        ymax=1.
        iflag=0
        do 20 j=2,jlm1
          if(i.eq.1) then
            ii=2
          elseif(ii.eq.il) then
            ii=ilm1
          else
            ii=i
          endif
          due=.5*(u(ii+1,j)-u(ii-1,j))
          dve=.5*(v(ii+1,j)-v(ii-1,j))
          dxe=.5*(x(i+1,j+1)+x(i+1,j)-x(i,j+1)-x(i,j))
          dye=.5*(y(i+1,j+1)+y(i+1,j)-y(i,j+1)-y(i,j))
          if(j.eq.2) then
            dun=u(i,j+1)-u(i,j)
            dvn=v(i,j+1)-v(i,j)
            dxn=.25*(x(i+1,j+2)+x(i,j+2)-x(i+1,j)-x(i,j))
            dyn=.25*(y(i+1,j+2)+y(i,j+2)-y(i+1,j)-y(i,j))
          else
            dun=.5*(u(i,j+1)-u(i,j-1))
            dvn=.5*(v(i,j+1)-v(i,j-1))
            dxn=.5*(x(i+1,j+1)+x(i,j+1)-x(i+1,j)-x(i,j))
            dyn=.5*(y(i+1,j+1)+y(i,j+1)-y(i+1,j)-y(i,j))
          endif
          dxy=1./(dxe*dyn-dxn*dye)
          dvx=(dve*dyn-dvn*dye)*dxy
          duy=(dun*dxe-due*dxn)*dxy
          yj=abs(sxp*(x(i,j+1)-x(i,2))+syp*(y(i,j+1)-y(i,2)))
          yi(j)=yj
          aomega(j)=abs(duy-dvx)
          fj=yj*aomega(j)*(1.-exp(-yj*ra))
          if(fj.gt.fmax.and.j.le.jlfm) then
            fmax=fj
            ymax=yj
            iflag=1
          endif
 20     continue
c=======================================
c        Compute wake factor
c=======================================
        fwake=ymax*fmax
        if(iflag.eq.1) then
          fwakep=cwk*ymax*u2max/fmax
          if(fwakep.lt.fwake) fwake=fwakep
        endif
        jep=jl
c=======================================
c       Compute eddy viscosity
c=======================================
        do 30 j=2,jlm1
          yj=yi(j)
          epij=ep(i,j)
          epi=rho(i,j,1)*((smlk*yj)*(1.-exp(-yj*ra)))**2*aomega(j)
          fkleb=1./(1.+5.5*(ckleb*yj/ymax)**6)
          epo=bigk*ccp*rho(i,j,1)*fwake*fkleb
          ep(i,j)=epi
          if(epo.lt.epi) ep(i,j)=epo
          if(j.lt.jep.and.epi.gt.epo) jep=j
          if(j.gt.jep) ep(i,j)=epo
          ep(i,j)=.5*(ep(i,j)+epij)
 30     continue
        ep(i,1)=0.
        ep(i,jl)=ep(i,jlm1)
 40   continue
      return
      end
      subroutine vcoefs
c=======================================================================
c         T U R B U L E N C E - M O D E L   P A R A M E T E R S
c=======================================================================
      include 'come2c.'
c=======================================================================
      clim=7./8.
      clim2=clim**2
      sigmad=.125
      alf0=1.
      alfs0=1.
      bbeta0=1.
      gamh0=1.
      rb=0.
      rk=0.
      rw=0.
c=============================
c    Stress-omega model
c=============================
      if(model.eq.1) then
        if(ntau.eq.0) then
          fmodel=0.
        else
          fmodel=1.
        endif
        csw1=1.-alphah
        csw2=2./3.*(1.-alphah-betah)
        csw3=2.*csw2
        if(nvisc.ne.0) then
          alf0=1./9.
          alfs0=beta/3.
          bbeta0=beta/(3.*betas)
          gamh0=21./2000.
          rb=8.
          rk=6.
          rw=22./9.
        endif
c=============================
c       k-omega model
c=============================
      else
        fmodel=1.
        if(nvisc.ne.0) then
          alf0=1./9.
          alfs0=beta/3.
          bbeta0=beta/(3.*betas)
          rb=8.
          rk=6.
          rw=2.61
        endif
      endif
      return
      end
      subroutine vgrads
c=======================================================================
c                V E L O C I T Y   G R A D I E N T S
c=======================================================================
      include 'come2c.'
c=======================================================================
      do 20 i=2,ilm1
        do 10 j=2,jlm1
          due=.5*(u(i+1,j)-u(i-1,j))
          dve=.5*(v(i+1,j)-v(i-1,j))
          dxe=.5*(x(i+1,j+1)+x(i+1,j)-x(i,j+1)-x(i,j))
          dye=.5*(y(i+1,j+1)+y(i+1,j)-y(i,j+1)-y(i,j))
          dun=.5*(u(i,j+1)-u(i,j-1))
          dvn=.5*(v(i,j+1)-v(i,j-1))
          dxn=.5*(x(i+1,j+1)+x(i,j+1)-x(i+1,j)-x(i,j))
          dyn=.5*(y(i+1,j+1)+y(i,j+1)-y(i+1,j)-y(i,j))
          dxy=1./(dxe*dyn-dxn*dye)
          dudx(i,j)=(due*dyn-dun*dye)*dxy
          dvdx(i,j)=(dve*dyn-dvn*dye)*dxy
          dudy(i,j)=(dun*dxe-due*dxn)*dxy
          dvdy(i,j)=(dvn*dxe-dve*dxn)*dxy
          vor(i,j)=4.*v(i,j)/(y(i,j+1)+y(i,j)+y(i+1,j+1)+y(i+1,j))
 10     continue
 20   continue
      return
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
      subroutine viskw
c=======================================================================
c        V I S C O U S   M O D S :  K - O M E G A   M O D E L
c=======================================================================
      include 'come2c.'
c=======================================================================
      rb4=rb**4
      do 20 i=1,il
        do 10 j=1,jl
          temp=abs(ei(i,j)/cv)
          rmu=viscos(temp)
          ret=ep(i,j)/rmu
          ret4=ret**4
          falfs(i,j)=(alfs0*rk+ret)/(rk+ret)
          falf(i,j)=(alf0*rw+ret)/((rw+ret)*falfs(i,j))
          fbet(i,j)=(bbeta0*rb4+ret4)/(rb4+ret4)
          ep(i,j)=falfs(i,j)*ep(i,j)
 10     continue
 20   continue
      return
      end
      subroutine vissw
c=======================================================================
c    V I S C O U S   M O D S :  S T R E S S - O M E G A   M O D E L
c=======================================================================
      include 'come2c.'
c=======================================================================
      do 20 i=1,il
        do 10 j=1,jl
          temp=abs(ei(i,j)/cv)
          rmu=viscos(temp)
          ret=ep(i,j)/rmu
          falfs(i,j)=(alfs0*rk+ret)/(rk+ret)
          falf(i,j)=(alf0*rw+ret)/(rw+ret)
     *             *(3.*rw+ret)/(3.*alfs0*rw+ret)
          retorb=(ret/rb)**4
          fbet(i,j)=(bbeta0+retorb)/(1.+retorb)
          falbet(i,j)=retorb/(1.+retorb)
          fgamh(i,j)=(gamh0+retorb)/(1.+retorb)
          fcc1(i,j)=(5./3.+retorb)/(1.+retorb)
          ep(i,j)=falfs(i,j)*ep(i,j)
 10     continue
 20   continue
      return
      end

