      program eddybl_start
c=======================================================================
c    W R I T E   E D G E   C O N D I T I O N   A N D   I N I T I A L
c                  P R O F I L E   D A T A   F I L E S
c-----------------------------------------------------------------------
c    Copyright (c) 2006 DCW Industries, Inc.  All rights reserved.
c=======================================================================
      common/close/    alpha  ,  beta   ,  betas  ,  sigma  ,  sigmas ,
     *                 prt
      common/closem/   alphah ,  betah  ,  gammah ,  cc1    ,  xis    ,
     *                 xmt0
      common/data/     cone   ,  ds     ,  g      ,  pr     ,  pt1    ,
     *                 r      ,  rek    ,  rmi0   ,  sgn    ,  sstop  ,
     *                 su     ,  tt1    ,  ustop  ,  viscon ,  vispow ,
     *                 wave   ,  xk     ,  xma    ,  ziotae ,  ziotal
      common/iodef/    iunit1 ,  iunit2 ,  iunit3 ,  iunit4 ,  iunit5 ,
     *                 iunit6 ,  iuplot ,  iutype
      common/keps/     ceps1  ,  ceps2  ,  cmu    ,  sigmae ,  sigmak ,
     *                 psieps
      common/moddef/   model  ,  nvisc
      common/ndata/    ibody  ,  ibound ,  iedge  ,  iend1  ,  ishort ,
     *                 j      ,  kebc   ,  kodwal ,  mstart ,  nflag  ,  
     *                 number
      common/pprint/   proval(10)
      common/pro1/     alamm  ,  cf     ,  delta  ,  h      ,  r2d    ,
     *                 rethet ,  si
      common/pro2/     betab  ,  drdz   ,  pe     ,  tw     ,  ue     ,
     *                 vw
      common/unityp/   heateq
c=======================================================================
c       Set logical unit numbers
c         and read input data
c=======================================
      iu15=15
      open(iu15,file='eddybl.dat',status='old',err=90)
      call namin(iu15)
c=======================================
c          Initialization
c=======================================
      if(iutype.eq.0) then
        heateq=778.26
      else
        heateq=1.
      endif
c=======================================
c      Generate edge conditions
c=======================================
      call table(nerror)
c=======================================
c  If successful...
c           generate initial profiles
c=======================================
      if(nerror.eq.0) then
        call interp(ierr)
        if(ierr.ne.0) then
          if(ierr.gt.0) then
            write(*,100) ierr
          else
            write(*,110)
          endif
          write(*,120)
          read(*,200) ians
          go to 10
        endif
        write(*,130)
c=======================================
c  Initial-profile display options...
c   idisk.ne.0  ==> save profiles in
c                   eddybl_start.prt
c   ivideo.ne.0 ==> display profiles
c                   on the console
c=======================================
        idisk=1
        ivideo=1
      endif
c=======================================
c        Now...generate them
c=======================================
      call start(idisk,ivideo,merror)
      if(merror.eq.0.and.idisk.eq.1) write(*,140)
      write(*,120)
      read(*,200) ians
      go to 10
 90   write(*,150)
 10   stop
c=======================================================================
 100  format(///,10x,' ERROR...Edge pressure exceeds total',
     *               ' pressure near point',i3)
 110  format(///,21x,' Diskfile TABLE.DAT does not exist ',
     *         //24x,'   If you get this message...   ',
     *          /24x,'  either your disk is full, or  ',
     *          /24x,'your computer is malfunctioning!')
 120  format(//26x,'Press ENTER to continue...',$)
 130  format(/20x,'Diskfile TABLE.DAT successfully written')
 140  format(//6x,'Initial profiles written in Ascii form to diskfile ',
     *            'EDDYBL_START.PRT')
 150  format(///19x,'ERROR: Diskfile EDDYBL.DAT is not present.')
 200  format(a1)
c=======================================================================
      end
      subroutine dissip(s,wscale,rghcon,beta50,akappa,deltap)
c=======================================================================
c    C O M P U T E   D I S S I P A T I O N - R A T E   P R O F I L E
c=======================================================================
      parameter (jk=301)
      common/array/    al(jk),dy(jk),et(jk),f(jk),rato(jk),t(jk),
     *                 uplus(jk),v(jk),w(jk),xl(jk),y(jk)
      common/extra/    tem(jk),w1(jk),w2(jk),w3(jk),y1(jk),y2(jk),
     *                 y3(jk),yodelt(jk),yplus(jk)
      common/lamina/   pp1(jk),pp2(jk),pp3(jk),qq1(jk),qq2(jk),qq3(jk),
     *                 rato2(jk)
      common/stressw/  uu(jk),uv(jk),vv(jk),ww(jk)
      common/mz/       deluv(jk)
      common/namo/     jsign  ,  re     ,  recov  ,  rw     ,  t1     ,
     *                 te     ,  u1     ,  vis1   ,  xend   ,  xi     ,
     *                 xnue
      common/close/    alpha  ,  beta   ,  betas  ,  sigma  ,  sigmas ,
     *                 prt
      common/closem/   alphah ,  betah  ,  gammah ,  cc1    ,  xis    ,
     *                 xmt0
      common/data/     cone   ,  ds     ,  g      ,  pr     ,  pt1    ,
     *                 r      ,  rek    ,  rmi0   ,  sgn    ,  sstop  ,
     *                 su     ,  tt1    ,  ustop  ,  viscon ,  vispow ,
     *                 wave   ,  xk     ,  xma    ,  ziotae ,  ziotal
      common/iodef/    iunit1 ,  iunit2 ,  iunit3 ,  iunit4 ,  iunit5 ,
     *                 iunit6 ,  iuplot ,  iutype
      common/keps/     ceps1  ,  ceps2  ,  cmu    ,  sigmae ,  sigmak ,
     *                 psieps
      common/moddef/   model  ,  nvisc
      common/ndata/    ibody  ,  ibound ,  iedge  ,  iend1  ,  ishort ,
     *                 j      ,  kebc   ,  kodwal ,  mstart ,  nflag  ,  
     *                 number
      common/pprint/   proval(10)
      common/pro1/     alamm  ,  cf     ,  delta  ,  h      ,  r2d    ,
     *                 rethet ,  si
      common/pro2/     betab  ,  drdz   ,  pe     ,  tw     ,  ue     ,
     *                 vw
      common/unityp/   heateq
c=======================================================================
      w(1)=wscale*s
      roottw=sqrt(tw/t1)
      coef3=sqrt(2./cf)/(roottw*ziotal*deltap)
      do 10 i=2,iedge
        w1(i)=s/(1.+rghcon*yplus(i))**2
        w3(i)=coef3*sqrt(et(i))
 10   continue
c=======================================
c  Turbulent case:  Find lower match
c    point and set lower w profile
c=======================================
      if(model.ge.0.and.ibound.ne.0) then
        do 20 i=2,iedge
          w2(i)=alamm/(beta50*akappa*yplus(i)*roottw)*sqrt(1.+t(i))
          if(w2(i).gt.w1(i)) go to 30
          w(i)=wscale*w1(i)
 20     continue
 30     iml=i
c=======================================
c   Turbulent case:  Find upper match
c  point and set intermediate w profile
c=======================================
        do 40 i=iml,iedge
          w2(i)=alamm/(beta50*akappa*yplus(i)*roottw)*sqrt(1.+t(i))
          if(w3(i).gt.w2(i)) go to 50
          w(i)=wscale*w2(i)
 40     continue
 50     imu=i
c=======================================
c    Laminar case:  Find lower match
c     point and set lower w profile
c=======================================
      else
        do 60 i=2,iedge
          if(w3(i).gt.w1(i)) go to 70
          w(i)=wscale*w1(i)
 60     continue
 70     iml=i
        imu=i
      endif
c=======================================
c        Compute outer profile
c=======================================
      do 80 i=imu,iedge
        w(i)=wscale*w3(i)
 80   continue
      return
      end
      subroutine energy(beta50,pih,deltap)
c=======================================================================
c   C O M P U T E   T U R B U L E N C E - E N E R G Y   P R O F I L E
c=======================================================================
      parameter (jk=301)
      common/array/    al(jk),dy(jk),et(jk),f(jk),rato(jk),t(jk),
     *                 uplus(jk),v(jk),w(jk),xl(jk),y(jk)
      common/extra/    tem(jk),w1(jk),w2(jk),w3(jk),y1(jk),y2(jk),
     *                 y3(jk),yodelt(jk),yplus(jk)
      common/lamina/   pp1(jk),pp2(jk),pp3(jk),qq1(jk),qq2(jk),qq3(jk),
     *                 rato2(jk)
      common/stressw/  uu(jk),uv(jk),vv(jk),ww(jk)
      common/mz/       deluv(jk)
      common/namo/     jsign  ,  re     ,  recov  ,  rw     ,  t1     ,
     *                 te     ,  u1     ,  vis1   ,  xend   ,  xi     ,
     *                 xnue
      common/close/    alpha  ,  beta   ,  betas  ,  sigma  ,  sigmas ,
     *                 prt
      common/closem/   alphah ,  betah  ,  gammah ,  cc1    ,  xis    ,
     *                 xmt0
      common/data/     cone   ,  ds     ,  g      ,  pr     ,  pt1    ,
     *                 r      ,  rek    ,  rmi0   ,  sgn    ,  sstop  ,
     *                 su     ,  tt1    ,  ustop  ,  viscon ,  vispow ,
     *                 wave   ,  xk     ,  xma    ,  ziotae ,  ziotal
      common/iodef/    iunit1 ,  iunit2 ,  iunit3 ,  iunit4 ,  iunit5 ,
     *                 iunit6 ,  iuplot ,  iutype
      common/keps/     ceps1  ,  ceps2  ,  cmu    ,  sigmae ,  sigmak ,
     *                 psieps
      common/moddef/   model  ,  nvisc
      common/ndata/    ibody  ,  ibound ,  iedge  ,  iend1  ,  ishort ,
     *                 j      ,  kebc   ,  kodwal ,  mstart ,  nflag  ,  
     *                 number
      common/pprint/   proval(10)
      common/pro1/     alamm  ,  cf     ,  delta  ,  h      ,  r2d    ,
     *                 rethet ,  si
      common/pro2/     betab  ,  drdz   ,  pe     ,  tw     ,  ue     ,
     *                 vw
      common/unityp/   heateq
c=======================================================================
c     Compute match point values
c=======================================
      escale=.5*alamm*cf/beta50
      et0=escale*tw*(deltap/20.)**2/t1
      do 10 i=1,iedge
        if(yplus(i).le.20.) then
          et(i)=et0*yodelt(i)**2
        else
          et(i)=escale*(1.+t(i))*cos(pih*yodelt(i))**2
          if(et(i).lt.ziotae.or.yodelt(i).gt.1.) go to 20
        endif
 10   continue
 20   ilow=i
      do 30 i=ilow,iedge
        et(i)=ziotae
 30   continue
      return
      end
      subroutine geom(xk,iedge)
c=======================================================================
c     G E N E R A T E   F I N I T E - D I F F E R E N C E   G R I D
c=======================================================================
      parameter (jk=301)
      common/array/    al(jk),dy(jk),et(jk),f(jk),rato(jk),t(jk),
     *                 uplus(jk),v(jk),w(jk),xl(jk),y(jk)
      common/extra/    tem(jk),w1(jk),w2(jk),w3(jk),y1(jk),y2(jk),
     *                 y3(jk),yodelt(jk),yplus(jk)
      common/lamina/   pp1(jk),pp2(jk),pp3(jk),qq1(jk),qq2(jk),qq3(jk),
     *                 rato2(jk)
      common/stressw/  uu(jk),uv(jk),vv(jk),ww(jk)
      common/mz/       deluv(jk)
      common/namo/     jsign  ,  re     ,  recov  ,  rw     ,  t1     ,
     *                 te     ,  u1     ,  vis1   ,  xend   ,  xi     ,
     *                 xnue
c=======================================================================
      y(1)=0.
      dy(1)=xend*(1.-xk)/(1.-xk**(iedge-1))
      do 10 n=2,iedge
        dy(n)=xk*dy(n-1)
        y(n)=y(n-1)+dy(n-1)
        denom=dy(n)+dy(n-1)
        y1(n)=2./(denom*dy(n))
        y2(n)=1./denom
        y3(n)=2./(denom*dy(n-1))
 10   continue
      return
      end
      subroutine interp(ierr)
c=======================================================================
c  C O M P U T E   I N I T I A L   C O N D I T I O N S   R E Q U I R E D
c               F O R   S U B R O U T I N E   S T A R T
c=======================================================================
      common/close/    alpha  ,  beta   ,  betas  ,  sigma  ,  sigmas ,
     *                 prt
      common/closem/   alphah ,  betah  ,  gammah ,  cc1    ,  xis    ,
     *                 xmt0
      common/data/     cone   ,  ds     ,  g      ,  pr     ,  pt1    ,
     *                 r      ,  rek    ,  rmi0   ,  sgn    ,  sstop  ,
     *                 su     ,  tt1    ,  ustop  ,  viscon ,  vispow ,
     *                 wave   ,  xk     ,  xma    ,  ziotae ,  ziotal
      common/iodef/    iunit1 ,  iunit2 ,  iunit3 ,  iunit4 ,  iunit5 ,
     *                 iunit6 ,  iuplot ,  iutype
      common/keps/     ceps1  ,  ceps2  ,  cmu    ,  sigmae ,  sigmak ,
     *                 psieps
      common/moddef/   model  ,  nvisc
      common/ndata/    ibody  ,  ibound ,  iedge  ,  iend1  ,  ishort ,
     *                 j      ,  kebc   ,  kodwal ,  mstart ,  nflag  ,  
     *                 number
      common/pprint/   proval(10)
      common/pro1/     alamm  ,  cf     ,  delta  ,  h      ,  r2d    ,
     *                 rethet ,  si
      common/pro2/     betab  ,  drdz   ,  pe     ,  tw     ,  ue     ,
     *                 vw
      common/unityp/   heateq
      dimension  as(50),ape(50),apep(50),aqw(50),aqwp(50),
     *           arhov(50),arhovp(50),armi(50),armip(50),atw(50),
     *           atwp(50),az(50),azp(50),azcrv(50),azcrvp(50)
c=======================================================================
c Open data file 'TABLE.DAT' or return
c with error flag set if it isn't there
c=======================================
      ierr=0
      open(iunit5,file='table.dat',form='unformatted',status='old',
     *                                                      err=10)
      go to 20
 10   ierr=-1
      go to 70
c=======================================
c  Compute freestream conditions from
c    the Rankine-Hugoniot relations
c=======================================
 20   gm1=g-1.
      gm1og=gm1/g
      xmac=1.+.5*gm1*xma**2
      rt1=pt1/(r*tt1)
      p1=pt1/xmac**(1./gm1og)
      r1=rt1/xmac**(1./gm1)
      aa1=sqrt(g*p1/r1)
      u1=xma*aa1
      tref=gm1og*u1**2/r
      xwave=.0174533*wave
      abc=(xma*sin(xwave))**2
      p10=pt1/(r1*u1*u1)
      if(abs(xwave).gt..0000001)p10=(1./(g*xma*xma))*((xmac*abc*(g+1.))
     *     /(abc*gm1+2.))**(1./gm1og)*((g+1.)/(2.*g*abc-gm1))**(1./gm1)
      t10=.5+1./(xma*xma*gm1)
c=======================================
c      Read data file 'TABLE.DAT'
c=======================================
      do 30 i=1,number
        read(iunit5) as(i),ape(i),apep(i),aqw(i),aqwp(i),
     *               arhov(i),arhovp(i),armi(i),armip(i),atw(i),
     *               atwp(i),az(i),azp(i),azcrv(i),azcrvp(i)
 30   continue
      close(iunit5)
c=======================================
c     Compute boundary-layer-edge
c       and surface conditions
c=======================================
      i=1
      if((as(i)-si)*(as(i+1)-si).le.0.) go to 60
      in=sign(1.,(as(i+1)-as(i))*(si-as(i)))
 50   if((i+in).le.0) go to 60
      if((i+in).ge.number) go to 60
      i=i+in
      if((as(i)-si)*(as(i+1)-si).gt.0.) go to 50
 60   continue
      hh=as(i+1)-as(i)
      dels=(si-as(i))/hh
      delsp=(as(i+1)-si)/hh
      phi1=dels**2*(3.-2.*dels)
      phi0=1.-phi1
      psi1=-hh*delsp*dels**2
      psi0=hh*dels*delsp**2
      phi1p=6.*dels*(1.-dels)/hh
      phi0p=-phi1p
      psi1p=-dels*(2.*delsp-dels)
      psi0p=delsp*(delsp-2.*dels)
      pe=ape(i)*phi0+ape(i+1)*phi1+apep(i)*psi0+apep(i+1)*psi1
      rmi=armi(i)*phi0+armi(i+1)*phi1+armip(i)*psi0+armip(i+1)*psi1
      tw=atw(i)*phi0+atw(i+1)*phi1+atwp(i)*psi0+atwp(i+1)*psi1
      rvwald=arhov(i)*phi0+arhov(i+1)*phi1+arhovp(i)*psi0
     *                                  +arhovp(i+1)*psi1
      dpeds=ape(i)*phi0p+ape(i+1)*phi1p+apep(i)*psi0p+apep(i+1)*psi1p
      rprime=armi(i)*phi0p+armi(i+1)*phi1p+armip(i)*psi0p
     *                                  +armip(i+1)*psi1p
      zprime=az(i)*phi0p+az(i+1)*phi1p+azp(i)*psi0p+azp(i+1)*psi1p
      drdz=0.
      if(abs(zprime).gt.1.e-10) drdz=rprime/zprime
c=======================================
c         Compute ue, te and re
c=======================================
      uesq=2.*t10*(1.-(pe/p10)**gm1og)
      if(uesq.lt.0.) then
        ierr=i
        go to 70
      endif
      ue=sqrt(uesq)
      te=t10-.5*ue*ue
      re=pe/(gm1og*te)
c=======================================
c  Compute the four quantities needed
c         for subroutine START
c=======================================
      rmi0=rmi
      betab=-2.*si*dpeds/(re*ue*ue)
      ue=u1*ue
      vw=rvwald*tw/(te*r1*re)
      pe=pe*r1*u1**2
      tw=tw*tref
 70   return
      end
      subroutine intgrl(disp,thet,h,iedge)
c=======================================================================
c       T R A P E Z O I D A L - R U L E   I N T E G R A T I O N
c=======================================================================
      parameter (jk=301)
      common/array/    al(jk),dy(jk),et(jk),f(jk),rato(jk),t(jk),
     *                 uplus(jk),v(jk),w(jk),xl(jk),y(jk)
      common/extra/    tem(jk),w1(jk),w2(jk),w3(jk),y1(jk),y2(jk),
     *                 y3(jk),yodelt(jk),yplus(jk)
      common/lamina/   pp1(jk),pp2(jk),pp3(jk),qq1(jk),qq2(jk),qq3(jk),
     *                 rato2(jk)
      common/stressw/  uu(jk),uv(jk),vv(jk),ww(jk)
      common/mz/       deluv(jk)
      common/namo/     jsign  ,  re     ,  recov  ,  rw     ,  t1     ,
     *                 te     ,  u1     ,  vis1   ,  xend   ,  xi     ,
     *                 xnue
c=======================================================================
      thet=0.
      disp=0.
      co=0.
      do 10 i=2,iedge
        disp=disp+.5*(2.+t(i)+t(i-1)-f(i)-f(i-1))*dy(i-1)
        c=f(i)*(1.-f(i))
        thet=thet+.5*(co+c)*dy(i-1)
        co=c
 10   continue
      h=disp/thet
      return
      end
      subroutine namin(iu15)
c=======================================================================
c          COORDINATE READING INPUT DATA FILE 'EDDYBL.DAT'
c=======================================================================
      parameter (ipro=10,niu=8,nke=6,nnid=2,nn=11,nxc=6,nxcm=6,nxd=20)
      character*12  a
      common/close/ xc(nxc)
      common/closem/xcm(nxcm)
      common/data/  xd(nxd)
      common/iodef/ iu(niu)
      common/keps/  xke(nke)
      common/moddef/nid(nnid)
      common/ndata/ n(nn)
      common/pprint/proval(ipro)
      common/pro1/  pr1(7)
c=======================================================================
      iu(1)=iu15
      do 10 i=2,niu-1
        iu(i)=i+5
 10   continue
      read(iu15,100) a,iu(8)
      do 20 i=1,nnid
        read(iu15,100) a,nid(i)
 20   continue
      do 30 i=1,nxc
        read(iu15,200) a,xc(i)
 30   continue
      do 40 i=1,nke
        read(iu15,200) a,xke(i)
 40   continue
      do 50 i=1,nxcm
        read(iu15,200) a,xcm(i)
 50   continue
      do 60 i=1,nxd
        read(iu15,200) a,xd(i)
 60   continue
      do 70 i=1,nn
        read(iu15,100) a,n(i)
 70   continue
      do 80 i=1,ipro
        read(iu15,200) a,proval(i)
 80   continue
      do 90 i=1,7
        read(iu15,200) a,pr1(i)
 90   continue
c=======================================
c  The rest of  eddybl.dat is read in
c  subroutine TABLE, so leave it open
c=======================================
      return
c=======================================================================
 100  format(1x,a12,i4)
 200  format(1x,a12,e13.6)
c=======================================================================
      end
      subroutine newton(xx,aa,bb)
c=======================================================================
c     S O L V E   A   T R A N S C E N D E N T A L   E Q U A T I O N
c              U S I N G   N E W T O N ' S   M E T H O D
c=======================================================================
c  Compatibility with Lahey "-dbl" compiler switch
c-----------------------------------------------------------------------
      double precision aa
c=======================================================================
      icount=0
      elog=dlog(aa)
 10   icount=icount+1
      ax=aa**(-xx)
      fx=ax+bb*xx-1.
      fp=-elog*ax+bb
      dx=-fx/fp
      xx=xx+dx
      if(abs(dx).gt.1.e-6*xx.and.icount.lt.10) go to 10
      return
      end
      subroutine noturb(alphab,tcoef,nedge,etamax)
c=======================================================================
c          C O M P U T E   L A M I N A R   P R O F I L E S
c=======================================================================
      parameter (jk=301)
      common/array/    al(jk),dy(jk),et(jk),f(jk),rato(jk),t(jk),
     *                 uplus(jk),v(jk),w(jk),xl(jk),y(jk)
      common/extra/    tem(jk),w1(jk),w2(jk),w3(jk),y1(jk),y2(jk),
     *                 y3(jk),yodelt(jk),yplus(jk)
      common/lamina/   pp1(jk),pp2(jk),pp3(jk),qq1(jk),qq2(jk),qq3(jk),
     *                 rato2(jk)
      common/stressw/  uu(jk),uv(jk),vv(jk),ww(jk)
      common/mz/       deluv(jk)
      common/namo/     jsign  ,  re     ,  recov  ,  rw     ,  t1     ,
     *                 te     ,  u1     ,  vis1   ,  xend   ,  xi     ,
     *                 xnue
      common/close/    alpha  ,  beta   ,  betas  ,  sigma  ,  sigmas ,
     *                 prt
      common/closem/   alphah ,  betah  ,  gammah ,  cc1    ,  xis    ,
     *                 xmt0
      common/data/     cone   ,  ds     ,  g      ,  pr     ,  pt1    ,
     *                 r      ,  rek    ,  rmi0   ,  sgn    ,  sstop  ,
     *                 su     ,  tt1    ,  ustop  ,  viscon ,  vispow ,
     *                 wave   ,  xk     ,  xma    ,  ziotae ,  ziotal
      common/iodef/    iunit1 ,  iunit2 ,  iunit3 ,  iunit4 ,  iunit5 ,
     *                 iunit6 ,  iuplot ,  iutype
      common/keps/     ceps1  ,  ceps2  ,  cmu    ,  sigmae ,  sigmak ,
     *                 psieps
      common/moddef/   model  ,  nvisc
      common/ndata/    ibody  ,  ibound ,  iedge  ,  iend1  ,  ishort ,
     *                 j      ,  kebc   ,  kodwal ,  mstart ,  nflag  ,  
     *                 number
      common/pprint/   proval(10)
      common/pro1/     alamm  ,  cf     ,  delta  ,  h      ,  r2d    ,
     *                 rethet ,  si
      common/pro2/     betab  ,  drdz   ,  pe     ,  tw     ,  ue     ,
     *                 vw
      common/unityp/   heateq
c=======================================================================
c           Initialization
c=======================================
      ipass=0
      xend=20.
      sum=su/te
      iedgem=iedge-1
      vispom=vispow-1.
c=======================================
c           Set up the grid
c=======================================
 10   if(jsign.lt.0.and.xend.gt.etamax) xend=etamax
      call geom(xk,iedge)
c=======================================
c       Initialize variables
c=======================================
      f(1)=0.
      if(rw.lt.1.) then
        t(1)=1./rw-1.
      else
        t(1)=0.
      endif
      ggo=0.
      tp1o=1.+t(1)
      xl(1)=(1.+sum)*tp1o**vispom/(tp1o+sum)
      rato2(1)=1.
      tint=0.
      icount=1
      do 20 n=2,iedge
        if(y(n).lt.5.) then
          f(n)=(y(n)/5.)**2
          t(n)=t(1)*(1.-f(n))
        else
          f(n)=1.
          t(n)=0.
        endif
        ggn=f(n) 
        v(n)=v(n-1)-.5*(ggn+ggo)*dy(n-1)
        ggo=ggn
        tp1n=1.+t(n) 
        xl(n)=(1.+sum)*tp1n**vispom/(tp1n+sum)
        tint=tint+.5*(tp1o+tp1n)*dy(n-1)
        rato2(n)=1.+tcoef*tint
        tp1o=tp1n
 20   continue
c=======================================
c  Solve momentum and energy equations
c      using Newton's iterations
c=======================================
 30   f(1)=0.
      t(1)=1./rw-1.
      pp1(1)=f(1)
      pp2(1)=0.
      pp3(1)=0.
      qq1(1)=t(1)
      qq2(1)=0.
      qq3(1)=0.
      epsm=.5*(xl(1)+xl(2)*rato2(2))
      do 40 n=2,iedgem
        epsp=.5*(xl(n)*rato2(n)+xl(n+1)*rato2(n+1))
        fy=(f(n+1)-f(n-1))*y2(n)
        disst=alphab*xl(n)*rato2(n)*fy
        a1=epsm*y3(n)+v(n)*y2(n)
        b1=-epsp*y1(n)-epsm*y3(n)-2.*betab*f(n)
        c1=epsp*y1(n)-v(n)*y2(n)
        d1=0.
        e1=betab
        f1=0.
        g1=-betab*(1.+f(n)**2)
        a2=-2.*disst*y2(n)
        b2=0.
        c2=-a2
        d2=epsm*y3(n)/pr+v(n)*y2(n)
        e2=-(epsp*y1(n)+epsm*y3(n))/pr
        f2=epsp*y1(n)/pr-v(n)*y2(n)
        g2=disst*fy
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
 40   continue
      kon=iedgem
      f(iedge)=1.
      t(iedge)=0.
      do 50 n=2,iedgem
        f(kon)=pp1(kon)+pp2(kon)*f(kon+1)+pp3(kon)*t(kon+1)
        if(f(kon).gt.1.) f(kon)=1.
        t(kon)=qq1(kon)+qq2(kon)*f(kon+1)+qq3(kon)*t(kon+1)
        kon=kon-1
 50   continue
c=======================================
c   Compute normal velocity profile
c=======================================
      ggo=0.
      do 60 n=2,iedge
        ggn=f(n) 
        v(n)=v(n-1)-.5*(ggn+ggo)*dy(n-1)
        ggo=ggn
 60   continue
c=======================================
c   Compute transverse curvature and
c          viscosity profiles
c=======================================
      tp1o=1.+t(1)
      xl(1)=(1.+sum)*tp1o**vispom/(tp1o+sum)
      rato(1)=1.
      rato2(1)=1.
      tint=0.
      do 70 n=2,iedge
        tp1n=1.+t(n) 
        xl(n)=(1.+sum)*tp1n**vispom/(tp1n+sum)
        tint=tint+.5*(tp1o+tp1n)*dy(n-1)
        rato2(n)=1.+tcoef*tint
        rato(n)=sqrt(rato2(n))
        tp1o=tp1n
 70   continue
      icount=icount+1
      if(icount.lt.10) go to 30
c=======================================
c  Locate boundary-layer edge & rescale
c=======================================
      do 80 n=1,iedge
        if(abs(1.-f(n)).lt..001.and.abs(t(n)).lt..001) go to 90
 80   continue
 90   nedge=n
      if(ipass.eq.0) then
        xend=1.2*y(nedge)
        ipass=1
        go to 10
      endif
      return
      end
      subroutine spline(n,vari,vard,vardp,li,ierr)
c=======================================================================
c    G E N E R A T E   S P L I N E   F I T   C O E F F I C I E N T S
c=======================================================================
      dimension vari(1),vard(1),vardp(1),p1(50),p2(50)
c=======================================================================
c    Make sure the table is in order
c=======================================
      ierr=0
      do 10 j=2,n
        if(vari(j).le.vari(j-1)) then
          ierr=1
          go to 20
        endif
 10   continue
 20   if(ierr.eq.1) then
        k=li
        write(*,100) j,k,(vari(j),j=1,n),(vard(j),j=1,n)
c=======================================
c    Set up the tridiagonal matrix
c=======================================
      else
        nm=n-1
        p1(1)=vardp(1)
        p2(1)=0.
        do 30 j=2,nm
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
 30     continue
c=======================================
c         Invert the matrix
c=======================================
        k=nm
        do 40 j=1,nm
          vardp(k)=p1(k)+p2(k)*vardp(k+1)
          k=k-1
 40     continue
      endif
      return
c=======================================================================
 100  format(//' Table below out of order at position  ',
     *             i5/1x,' Table Number  ',i6,//(4g15.8))
c=======================================================================
      end
      subroutine start(idisk,ivideo,merror)
c=======================================================================
c        G E N E R A T E   S T A R T I N G   P R O F I L E S
c                 F O R   P R O G R A M   E D D Y B L
c=======================================================================
c  Compatibility with Lahey "-dbl" compiler switch
c-----------------------------------------------------------------------
      double precision aa
c-----------------------------------------------------------------------
      parameter (jk=301)
      common/array/    al(jk),dy(jk),et(jk),f(jk),rato(jk),t(jk),
     *                 uplus(jk),v(jk),w(jk),xl(jk),y(jk)
      common/extra/    tem(jk),w1(jk),w2(jk),w3(jk),y1(jk),y2(jk),
     *                 y3(jk),yodelt(jk),yplus(jk)
      common/lamina/   pp1(jk),pp2(jk),pp3(jk),qq1(jk),qq2(jk),qq3(jk),
     *                 rato2(jk)
      common/stressw/  uu(jk),uv(jk),vv(jk),ww(jk)
      common/mz/       deluv(jk)
      common/namo/     jsign  ,  re     ,  recov  ,  rw     ,  t1     ,
     *                 te     ,  u1     ,  vis1   ,  xend   ,  xi     ,
     *                 xnue
      common/close/    alpha  ,  beta   ,  betas  ,  sigma  ,  sigmas ,
     *                 prt
      common/closem/   alphah ,  betah  ,  gammah ,  cc1    ,  xis    ,
     *                 xmt0
      common/data/     cone   ,  ds     ,  g      ,  pr     ,  pt1    ,
     *                 r      ,  rek    ,  rmi0   ,  sgn    ,  sstop  ,
     *                 su     ,  tt1    ,  ustop  ,  viscon ,  vispow ,
     *                 wave   ,  xk     ,  xma    ,  ziotae ,  ziotal
      common/iodef/    iunit1 ,  iunit2 ,  iunit3 ,  iunit4 ,  iunit5 ,
     *                 iunit6 ,  iuplot ,  iutype
      common/keps/     ceps1  ,  ceps2  ,  cmu    ,  sigmae ,  sigmak ,
     *                 psieps
      common/moddef/   model  ,  nvisc
      common/ndata/    ibody  ,  ibound ,  iedge  ,  iend1  ,  ishort ,
     *                 j      ,  kebc   ,  kodwal ,  mstart ,  nflag  ,  
     *                 number
      common/pprint/   proval(10)
      common/pro1/     alamm  ,  cf     ,  delta  ,  h      ,  r2d    ,
     *                 rethet ,  si
      common/pro2/     betab  ,  drdz   ,  pe     ,  tw     ,  ue     ,
     *                 vw
      common/unityp/   heateq
      data akappa/.41/,pih/1.5707963/
c=======================================================================
c       Freestream conditions
c=======================================
      merror=1
      jsign=int(sgn)
      beta50=sqrt(betas)
      recov=prt
      gm1=g-1.
      gm1og=gm1/g
      xmac=1.+.5*gm1*xma**2
      rt1=pt1/(r*tt1)
      p1=pt1/xmac**(1./gm1og)
      r1=rt1/xmac**(1./gm1)
      t1=tt1/xmac
      aa1=sqrt(g*p1/r1)
      u1=xma*aa1
      tref=gm1og*u1**2/r
      xwave=.0174533*wave
      abc=(xma*sin(xwave))**2
      if(abs(xwave).le..0000001) then
        p10=pt1/(r1*u1*u1)
      else
        p10=(1./(g*xma*xma))*((xmac*abc*(g+1.))
     *     /(abc*gm1+2.))**(1./gm1og)*((g+1.)/(2.*g*abc-gm1))**(1./gm1)
      endif
      t10=.5+1./(xma*xma*gm1)
      vis1=viscon*t1**vispow/(t1+su)
      visref=viscon*tref**vispow/(tref+su)
      rey=r1*u1/vis1
      reyref=rey*vis1/visref
      tc=su/(t1*xmac)
c=======================================
c  Provision for stagnation-point flow
c=======================================
      if((model.lt.0.or.ibound.eq.0).and.ibody.eq.1) then
        si=ds
        rmi0=si
      elseif(si.eq.0.) then
        si=ds
      endif
c=======================================
c    Boundary-layer-edge conditions
c=======================================
      phi=atan(drdz)
      costh=cos(phi)
      ue=ue/u1
      pe=pe/(r1*u1**2)
      te=t10-.5*ue*ue
      re=pe/(gm1og*te)
      xnue=te**vispow*(1.+t10*tc)/(te+t10*tc)
      ziotas=ziotal
      rw=tref*te/tw
      vw=vw/u1
      eps2=1./reyref
      eps=sqrt(eps2)
      if(j.eq.1) then
        r2d=1.
        rmi1j=rmi0
        rmi2j=rmi1j**2
      else
        rmi1j=1.
        rmi2j=1.
      endif
c=======================================
c           Turbulent flow
c=======================================
      if(model.ge.0.and.ibound.gt.0) then
        ziotal=.049295
        visw=viscon*tw**vispow/(tw+su)
        scf=sqrt(.5*cf/rw)      
        cq=rw*vw/ue
        xi=(eps*xnue*rmi1j)**2*abs(rethet-cq)/cf
        oy=2.*xi
        oz=sqrt(oy)
        yscale=(1.+j)*eps*oz/(re*ue*rmi2j)
        deltap=scf*rw*ue*delta*visref/(eps2*visw)
        xl(1)=rw*visw/(visref*xnue)
        elam=eps*xnue*rmi1j/oz
        redelb=re*ue*delta*(1.+.5*j*sgn*delta*costh/rmi1j)/(xnue*eps2)
        etarat=redelb-.5*(g-1.)*recov*xma*xma*rw*rethet-(1.-rw)*h*rethet
        etadel=elam*etarat
        xend=1.2*etadel
        if(jsign.lt.0) then
          if(j.eq.0) then
            etamax=ue*r2d/(eps*oz)
          else
            etamax=ue*rmi2j/(2.*eps*oz)
          endif
          if(xend.gt.etamax) then
            write(*,180)
            go to 90
          endif
        endif
c-----------------------------
c     Construct the mesh
c-----------------------------
        aa=.5*etarat*cf/xl(1)
        bb=.5*(3.*(redelb-(1.+h)*rethet)/etarat-1.)
        xx=1./bb
        call newton(xx,aa,bb)
        call geom(xk,iedge)
c-----------------------------
c Streamwise velocity profile
c-----------------------------
        call veloc(oz,eps,rmi1j,ye,ym,xx)
c-----------------------------
c  Temperature and viscosity
c-----------------------------
        call temper
c=======================================
c            Laminar flow
c=======================================
      else
        if(ibody.eq.2) then
          xi=re*ue*xnue*rmi2j*si
        else
          xi=.5*re*ue*xnue*rmi2j*ds/(j+1.)
          betab=1./(1.+j)
        endif
        oy=2.*xi
        oz=sqrt(oy)
        alphab=ue**2/te
        rw=tref*te/tw
        tcoef=2.*sgn*eps*oz*costh/(re*ue*rmi2j)
        v(1)=oz*rw*vw/(eps*re*ue*xnue*rmi1j)
        if(jsign.lt.0) then
          if(j.eq.0) then
            etamax=ue*r2d/(eps*oz)
          else
            etamax=ue*rmi2j/(2.*eps*oz)
          endif
        endif
        call noturb(alphab,tcoef,nedge,etamax)
        visw=viscon*tw**vispow/(tw+su)
        cf=2.*xl(1)*eps*xnue*rmi1j*(f(2)-f(1))/(oz*dy(1))
        scf=sqrt(.5*cf/rw)      
        yscale=(1.+j)*eps*oz/(re*ue*rmi2j)
        deltap=scf*rw*ue*visref/(eps2*visw)
      endif
c=======================================
c      Compute y/delta and yplus
c=======================================
      ybar=0.
      rato(1)=1.
      yodelt(1)=0.
      yplus(1)=0.
      do 10 i=2,iedge
        ybar=ybar+.5*(2.+t(i-1)+t(i))*yscale*dy(i-1)
        if(j.eq.0) then
          rato(i)=1.
          yodelt(i)=ybar
        else
          rato(i)=sqrt(1.+sgn*ybar)
          yodelt(i)=sgn*rmi1j*(rato(i)-1.)
        endif
        if(model.ge.0.and.ibound.ne.0) yodelt(i)=yodelt(i)/delta
        yplus(i)=deltap*yodelt(i)
 10   continue
      if(model.lt.0.or.ibound.eq.0) then
        delta=yodelt(nedge)
        deltap=delta*deltap
        do 15 i=2,iedge
          yodelt(i)=yodelt(i)/delta
 15     continue
      endif
c=======================================
c     Rescale the velocity profile
c=======================================
      if(model.ge.0.and.ibound.ne.0) then
        call veloc2(ym,iadj)
        call intgrl(disp,thet,hb,iedge)
        tscale=yscale*rmi1j/(1.+j)
        thet=tscale*thet
        disp=tscale*disp
        rthetb=re*ue*thet/(eps2*xnue)
        xi=xi*(rethet/rthetb)**2
        rthetb=rethet
        oy=2.*xi
        oz=sqrt(oy)
        yscale=(1.+j)*eps*oz/(re*ue*rmi2j)
        elam=eps*xnue*rmi1j/oz
        etadel=y(iadj)
        ctemp=6.*(h-hb)*elam*rethet/etadel**3
        vispom=vispow-1.
        sum=su/t1
        iadjm=iadj-1
        do 20 i=2,iadjm
          t(i)=t(i)+ctemp*y(i)*(etadel-y(i))
          xl(i)=(1.+t(i))**vispom*(1.+sum)/(1.+t(i)+sum)
 20     continue
c=======================================
c   Compute normal velocity profile
c=======================================
        v(1)=oz*rw*vw/(eps*re*ue*xnue*rmi1j)
        ggo=0.
        do 30 k=2,iedge
          gg=f(k) 
          v(k)=v(k-1)-.5*(gg+ggo)*dy(k-1)
          ggo=gg
 30     continue
      endif
c=======================================
c  Make sure yplus(2) is small enough
c=======================================
      if(yplus(2).gt.1.) then
        write(*,100) yplus(2)
c=======================================
c   Compute turbulence energy profile
c=======================================
      else
        if(model.ge.0.and.ibound.ne.0) then
          call energy(beta50,pih,deltap)
        else
          do 40 k=2,iedge
            et(k)=ziotae
 40       continue
          et(1)=0.
        endif
c=======================================
c    Compute Stress-omega quantities
c=======================================
        uu(1)=0.
        uv(1)=0.
        vv(1)=0.
        ww(1)=0.
        deluv(1)=0.
        do 50 k=2,iedge
          uu(k)=.243*et(k)
          uv(k)=.3*et(k)
          vv(k)=-.239*et(k)
          ww(k)=-.004*et(k)
          deluv(k)=0.
 50     continue
c=======================================
c   Compute dissipation-rate profile
c=======================================
        vwp=vw/(ue*scf)
        if(vwp.lt.0.) vwp=0.
        akplus=scf*rek*(t1/tw)**2/xl(1)
        sri=.01*akplus
        if(akplus.lt.25.) sri=4.*sri**2
        sbi=.05*vwp*(1.+5.*vwp)
        s=1./(sri+sbi)
        wscale=.5*oy*cf*rw/(eps2*xnue*xl(1))
        rghcon=sqrt(beta*s/6.)
        call dissip(s,wscale,rghcon,beta50,akappa,deltap)
        w(iedge)=w(iedge)*ziotal/ziotas
c=======================================
c         Prepare for printing
c=======================================
        call intgrl(disp,thet,hb,iedge)
        tscale=yscale*rmi1j/(1.+j)
        thet=tscale*thet
        disp=tscale*disp
        rthetb=re*ue*thet/(eps2*xnue)
        if(idisk.eq.0.and.ivideo.eq.0) go to 70
        elscal=oy/(ue*delta)
        do 60 i=1,iedge
          uplus(i)=f(i)/scf
          tem(i)=1.+t(i)
          al(i)=elscal*sqrt(et(i)/beta50)/w(i)
 60     continue
c==============================================
c   Write the profiles to 'eddybl_start.prt'
c==============================================
        if(idisk.ne.0) then
          open(iunit2,file='eddybl_start.prt',status='unknown')
          write(iunit2,110)
          write(iunit2,120) (n,y(n),yodelt(n),f(n),tem(n),
     *                       yplus(n),uplus(n),et(n),al(n),n=1,iedge)
          write(iunit2,130) thet,disp,hb,rthetb,cf,yplus(2)
          close(iunit2)
        endif
c=======================================
c     Display profiles on video
c=======================================
        if(ivideo.ne.0) then
          write(*,110)
          do 65 n=1,iedge
            write(*,120) n,y(n),yodelt(n),f(n),tem(n),
     *                   yplus(n),uplus(n),et(n),al(n)
 65       continue
        endif
 70     write(*,130) thet,disp,hb,rthetb,cf,yplus(2)
c-----------------------------
c  Warning for k-epsilon if 
c   yplus(2) is too small
c-----------------------------
        if(model.eq.2.and.yplus(2).lt..1) then
          write(*,140) yplus(2)
        endif
        dhelp=re*ue*xnue*rmi2j
        if((model.lt.0.or.ibound.eq.0).and.ibody.eq.1) then
          done=.5*dhelp*ds/(j+1.)
          d1help=0.
          rstrt=(2.*ds)**j
          si1=0.
          ustrt=2.*ue
          wten=0.
        else
          done=dhelp*ds
          d1help=dhelp
          rstrt=rmi1j
          si1=si-ds
          ustrt=ue
          wten=w(iedge)
          ibody=2
        endif
        dtwo=done
        eten=et(iedge)
        eteo=eten
        itstop=1
        restrt=re
        rstrtm=rmi1j
        testrt=te
        ustrtm=ue
        wteo=w(iedge)
        xistrt=xi
        mstart=2
c=======================================
c    Write file 'INPUT.DAT' on disk
c=======================================
        open(iunit5,file='input.dat',form='unformatted',
     *                                 status='unknown')
        write(iunit5) mstart,ibound,iedge
        do 80 n=1,iedge
          if((model.lt.0.or.ibound.eq.0).and.ibody.eq.1) then
            write(iunit5) dy(n),et(n),et(n),f(n),f(n),t(n),t(n),
     *                    v(n),v(n),0.,w(n),rato(n),rato(n),
     *                    ww(n),ww(n),uu(n),uu(n),uv(n),uv(n),
     *                    vv(n),vv(n),deluv(n)
          else
            write(iunit5) dy(n),et(n),et(n),f(n),f(n),t(n),t(n),
     *                    v(n),v(n),w(n),w(n),rato(n),rato(n),
     *                    ww(n),ww(n),uu(n),uu(n),uv(n),uv(n),
     *                    vv(n),vv(n),deluv(n)
          endif
 80     continue
        write(iunit5) dhelp,done,dtwo,d1help,eten,eteo,itstop,
     *                restrt,rstrt,rstrtm,si,si1,testrt,ustrt,
     *                ustrtm,wten,wteo,xistrt
        endfile iunit5
        close(iunit5)
        write(*,150)
      endif
      merror=0
      ziotal=ziotas
 90   return
c=======================================================================
 100  format(//20x,'YPLUS(2) = ',1pe13.6,' IS TOO LARGE ')
 110  format(/3x,1hn,3x,3heta,4x,7hy/delta,4x,4hu/Ue,5x,4hT/Te,
     *        4x,5hyplus,4x,5huplus,3x,7hk/Ue**2,2x,7hL/delta)
 120  format(1x,i3,1p8e9.2)
 130  format(//27x'** INTEGRAL PARAMETERS **'/29x,'THETA   = ',1pe11.4/
     *                   29x,'DELTA*  = ',e11.4/29x,'H       = ',e11.4/
     *                   29x,'RETHETA = ',e11.4/29x,'CF      = ',e11.4/
     *                   29x,'YPLUS(2)= ',e11.4)
 140  format(//10x,'WARNING: This run will start with ',
     *             'YPLUS(2) = ',1pe11.4,'.'/
     *         16x,'Having YPLUS(2) less than 0.1 sometimes inhibits'/
     *         22x,'convergence for the k-epsilon model.')
 150  format(/19x,'Diskfile INPUT.DAT successfully written')
 180  format(//22x,'MESH EXTENDS BEYOND THE CENTERLINE')
c=======================================================================
      end
      subroutine table(nerror)
c=======================================================================
c  S E T   U P   T H E   B O U N D A R Y - L A Y E R - E D G E   A N D
c S U R F A C E   C O N D I T I O N S   N E E D E D   B Y   E D D Y B L 
c=======================================================================
      common/close/    alpha  ,  beta   ,  betas  ,  sigma  ,  sigmas ,
     *                 prt
      common/closem/   alphah ,  betah  ,  gammah ,  cc1    ,  xis    ,
     *                 xmt0
      common/data/     cone   ,  ds     ,  g      ,  pr     ,  pt1    ,
     *                 r      ,  rek    ,  rmi0   ,  sgn    ,  sstop  ,
     *                 su     ,  tt1    ,  ustop  ,  viscon ,  vispow ,
     *                 wave   ,  xk     ,  xma    ,  ziotae ,  ziotal
      common/iodef/    iunit1 ,  iunit2 ,  iunit3 ,  iunit4 ,  iunit5 ,
     *                 iunit6 ,  iuplot ,  iutype
      common/keps/     ceps1  ,  ceps2  ,  cmu    ,  sigmae ,  sigmak ,
     *                 psieps
      common/moddef/   model  ,  nvisc
      common/ndata/    ibody  ,  ibound ,  iedge  ,  iend1  ,  ishort ,
     *                 j      ,  kebc   ,  kodwal ,  mstart ,  nflag  ,  
     *                 number
      common/pprint/   proval(10)
      common/pro1/     alamm  ,  cf     ,  delta  ,  h      ,  r2d    ,
     *                 rethet ,  si
      common/pro2/     betab  ,  drdz   ,  pe     ,  tw     ,  ue     ,
     *                 vw
      common/unityp/   heateq
      dimension peb(50),pebp(50),rvwald(50),rwp(50),s(50),
     *          zcurv(50),zcurvp(50),rmi(50),rp(50),z(50),
     *          zp(50),qw(50),qwp(50),twb(50),twbp(50)
      data rmi/50*1./
      data pebp,rwp,zcurvp,rvwald,zcurv,s,rp,z,zp,qw,qwp,twbp/600*0./
c========================================================================
c      Display console message
c=======================================
      nerror=0
      write(*,200)
c=======================================
c    Compute dimensional quantities
c=======================================
      gm1=g-1.
      gm1og=gm1/g
      xmac=1.+.5*gm1*xma**2
      rt1=pt1/(r*tt1)
      r1=rt1/xmac**(1./gm1)
      t1=tt1/xmac
      aa1=sqrt(g*r*t1)
      u1=xma*aa1
      tref=gm1og*u1**2/r
      visref=viscon*tref**vispow/(tref+su)
c=======================================
c   Read edge and wall distributions
c   from input data file EDDYBL.DAT
c=======================================
      read(iunit1,120) (s(i),peb(i),twb(i),qw(i),rvwald(i),rmi(i),
     *                  z(i),zcurv(i), i=1,number)
      read(iunit1,130) pebp(1),pebp(number)
      read(iunit1,130) twbp(1),twbp(number)
      read(iunit1,130) qwp(1),qwp(number)
      read(iunit1,130) rwp(1),rwp(number)
      read(iunit1,130) rp(1),rp(number)
      read(iunit1,130) zp(1),zp(number)
      read(iunit1,130) zcurvp(1),zcurvp(number)
      close(iunit1)
c=======================================
c        Nondimensionalization
c=======================================
      conp=1./(r1*u1*u1)
      cont=1./tref
      conq=heateq/(visref*u1*u1)
      do 20 i=1,number
        peb(i)=conp*peb(i)
        twb(i)=cont*twb(i)
        qw(i)=conq*qw(i)
        pebp(i)=conp*pebp(i)
        twbp(i)=cont*twbp(i)
        qwp(i)=conq*qwp(i)
 20   continue
c=======================================
c   Generate spline-fit coefficients
c=======================================
      call spline(number,s,peb,pebp,1,ierr)
      if(ierr.ne.0) go to 70
      call spline(number,s,rvwald,rwp,2,ierr)
      call spline(number,s,zcurv,zcurvp,3,ierr)
      call spline(number,s,rmi,rp,4,ierr)
      call spline(number,s,z,zp,5,ierr)
      call spline(number,s,qw,qwp,6,ierr)
      call spline(number,s,twb,twbp,7,ierr)
c=======================================
c  Write the coefficients to 'TABLE.DAT'
c=======================================
      open(iunit5,file='table.dat',form='unformatted',
     *                               status='unknown')
      do 30 i=1,number
        write(iunit5) s(i),peb(i),pebp(i),qw(i),qwp(i),rvwald(i),rwp(i),
     *                rmi(i),rp(i),twb(i),twbp(i),z(i),zp(i),zcurv(i),
     *                zcurvp(i)
 30   continue
      endfile iunit5
      close(iunit5)
      go to 80
c=======================================
c          Error handling
c=======================================
 70   close(iunit5)
      nerror=1
 80   return
c=======================================================================
 120  format(8e14.6)
 130  format(2e14.6)
 200  format(1x,'Generating edge conditions...')
c=======================================================================
      end
      subroutine temper
c=======================================================================
c       C O M P U T E   T E M P E R A T U R E   P R O F I L E
c=======================================================================
      parameter (jk=301)
      common/array/    al(jk),dy(jk),et(jk),f(jk),rato(jk),t(jk),
     *                 uplus(jk),v(jk),w(jk),xl(jk),y(jk)
      common/extra/    tem(jk),w1(jk),w2(jk),w3(jk),y1(jk),y2(jk),
     *                 y3(jk),yodelt(jk),yplus(jk)
      common/lamina/   pp1(jk),pp2(jk),pp3(jk),qq1(jk),qq2(jk),qq3(jk),
     *                 rato2(jk)
      common/stressw/  uu(jk),uv(jk),vv(jk),ww(jk)
      common/mz/       deluv(jk)
      common/namo/     jsign  ,  re     ,  recov  ,  rw     ,  t1     ,
     *                 te     ,  u1     ,  vis1   ,  xend   ,  xi     ,
     *                 xnue
      common/close/    alpha  ,  beta   ,  betas  ,  sigma  ,  sigmas ,
     *                 prt
      common/closem/   alphah ,  betah  ,  gammah ,  cc1    ,  xis    ,
     *                 xmt0
      common/data/     cone   ,  ds     ,  g      ,  pr     ,  pt1    ,
     *                 r      ,  rek    ,  rmi0   ,  sgn    ,  sstop  ,
     *                 su     ,  tt1    ,  ustop  ,  viscon ,  vispow ,
     *                 wave   ,  xk     ,  xma    ,  ziotae ,  ziotal
      common/iodef/    iunit1 ,  iunit2 ,  iunit3 ,  iunit4 ,  iunit5 ,
     *                 iunit6 ,  iuplot ,  iutype
      common/keps/     ceps1  ,  ceps2  ,  cmu    ,  sigmae ,  sigmak ,
     *                 psieps
      common/moddef/   model  ,  nvisc
      common/ndata/    ibody  ,  ibound ,  iedge  ,  iend1  ,  ishort ,
     *                 j      ,  kebc   ,  kodwal ,  mstart ,  nflag  ,  
     *                 number
      common/pprint/   proval(10)
      common/pro1/     alamm  ,  cf     ,  delta  ,  h      ,  r2d    ,
     *                 rethet ,  si
      common/pro2/     betab  ,  drdz   ,  pe     ,  tw     ,  ue     ,
     *                 vw
      common/unityp/   heateq
c=======================================================================
      thw=tw/t1-1.
      coef2=.5*(g-1.)*recov*xma*xma
      vispom=vispow-1.
      sum=su/t1
      do 10 i=1,iedge
        t(i)=(thw+coef2*f(i))*(1.-f(i))
        xl(i)=(1.+t(i))**vispom*(1.+sum)/(1.+t(i)+sum)
 10   continue
      return
      end
      subroutine upchar(s)
c=======================================================================
c    C O N V E R T   C H A R A C T E R   T O   U P P E R   C A S E
c=======================================================================
      character*1 s
c=======================================================================
      if(s.ge.'a'.and.s.le.'z') s=char(ichar(s)-32)
      return
      end
      subroutine veloc(oz,eps,rmi1j,ye,ym,xx)
c=======================================================================
c         C O M P U T E   V E L O C I T Y   P R O F I L E
c=======================================================================
      parameter (jk=301)
      common/array/    al(jk),dy(jk),et(jk),f(jk),rato(jk),t(jk),
     *                 uplus(jk),v(jk),w(jk),xl(jk),y(jk)
      common/extra/    tem(jk),w1(jk),w2(jk),w3(jk),y1(jk),y2(jk),
     *                 y3(jk),yodelt(jk),yplus(jk)
      common/lamina/   pp1(jk),pp2(jk),pp3(jk),qq1(jk),qq2(jk),qq3(jk),
     *                 rato2(jk)
      common/stressw/  uu(jk),uv(jk),vv(jk),ww(jk)
      common/mz/       deluv(jk)
      common/namo/     jsign  ,  re     ,  recov  ,  rw     ,  t1     ,
     *                 te     ,  u1     ,  vis1   ,  xend   ,  xi     ,
     *                 xnue
      common/close/    alpha  ,  beta   ,  betas  ,  sigma  ,  sigmas ,
     *                 prt
      common/closem/   alphah ,  betah  ,  gammah ,  cc1    ,  xis    ,
     *                 xmt0
      common/data/     cone   ,  ds     ,  g      ,  pr     ,  pt1    ,
     *                 r      ,  rek    ,  rmi0   ,  sgn    ,  sstop  ,
     *                 su     ,  tt1    ,  ustop  ,  viscon ,  vispow ,
     *                 wave   ,  xk     ,  xma    ,  ziotae ,  ziotal
      common/iodef/    iunit1 ,  iunit2 ,  iunit3 ,  iunit4 ,  iunit5 ,
     *                 iunit6 ,  iuplot ,  iutype
      common/keps/     ceps1  ,  ceps2  ,  cmu    ,  sigmae ,  sigmak ,
     *                 psieps
      common/moddef/   model  ,  nvisc
      common/ndata/    ibody  ,  ibound ,  iedge  ,  iend1  ,  ishort ,
     *                 j      ,  kebc   ,  kodwal ,  mstart ,  nflag  ,  
     *                 number
      common/pprint/   proval(10)
      common/pro1/     alamm  ,  cf     ,  delta  ,  h      ,  r2d    ,
     *                 rethet ,  si
      common/pro2/     betab  ,  drdz   ,  pe     ,  tw     ,  ue     ,
     *                 vw
      common/unityp/   heateq
c=======================================================================
c     Compute profile parameters
c=======================================
      xlw=xl(1)
      coef1=oz/(eps*xnue*rmi1j)
      coef3=2.*xlw/(coef1*cf)
      redelt=re*ue*delta*(1.+.5*delta/rmi1j)**j/(eps*eps*xnue)
c=======================================
c    Compute exponent reciprocal, enr
c=======================================
      enr=(xx+2.)/(xx-1.)
c=======================================
c         Compute match point
c=======================================
      ye=redelt/(coef1*(1.+(tw/t1-1.)/(enr+1.)
     *  +.5*(g-1.)*recov*xma*xma*enr/((enr+1.)*(enr+2.))))
      ym=(coef3**enr/ye)**(1./(enr-1.))
      coef3=1./coef3
      expon=1./enr
c=======================================
c     Compute analytical profile
c=======================================
      do 10 i=1,iedge
        if(y(i).ge.ye) then
          f(i)=1.
        elseif(y(i).le.ym) then
          f(i)=coef3*y(i)
        else
          f(i)=(y(i)/ye)**expon
        endif
 10   continue
      return
      end
      subroutine veloc2(ym,iadj)
c=======================================================================
c          R E S C A L E   V E L O C I T Y   P R O F I L E
c=======================================================================
      parameter (jk=301)
      common/array/    al(jk),dy(jk),et(jk),f(jk),rato(jk),t(jk),
     *                 uplus(jk),v(jk),w(jk),xl(jk),y(jk)
      common/extra/    tem(jk),w1(jk),w2(jk),w3(jk),y1(jk),y2(jk),
     *                 y3(jk),yodelt(jk),yplus(jk)
      common/lamina/   pp1(jk),pp2(jk),pp3(jk),qq1(jk),qq2(jk),qq3(jk),
     *                 rato2(jk)
      common/stressw/  uu(jk),uv(jk),vv(jk),ww(jk)
      common/mz/       deluv(jk)
      common/namo/     jsign  ,  re     ,  recov  ,  rw     ,  t1     ,
     *                 te     ,  u1     ,  vis1   ,  xend   ,  xi     ,
     *                 xnue
      common/close/    alpha  ,  beta   ,  betas  ,  sigma  ,  sigmas ,
     *                 prt
      common/closem/   alphah ,  betah  ,  gammah ,  cc1    ,  xis    ,
     *                 xmt0
      common/data/     cone   ,  ds     ,  g      ,  pr     ,  pt1    ,
     *                 r      ,  rek    ,  rmi0   ,  sgn    ,  sstop  ,
     *                 su     ,  tt1    ,  ustop  ,  viscon ,  vispow ,
     *                 wave   ,  xk     ,  xma    ,  ziotae ,  ziotal
      common/iodef/    iunit1 ,  iunit2 ,  iunit3 ,  iunit4 ,  iunit5 ,
     *                 iunit6 ,  iuplot ,  iutype
      common/keps/     ceps1  ,  ceps2  ,  cmu    ,  sigmae ,  sigmak ,
     *                 psieps
      common/moddef/   model  ,  nvisc
      common/ndata/    ibody  ,  ibound ,  iedge  ,  iend1  ,  ishort ,
     *                 j      ,  kebc   ,  kodwal ,  mstart ,  nflag  ,  
     *                 number
      common/pprint/   proval(10)
      common/pro1/     alamm  ,  cf     ,  delta  ,  h      ,  r2d    ,
     *                 rethet ,  si
      common/pro2/     betab  ,  drdz   ,  pe     ,  tw     ,  ue     ,
     *                 vw
      common/unityp/   heateq
c=======================================================================
      scf=sqrt(.5*cf/rw)
      do 10 i=1,iedge
        if(yodelt(i).ge.1.) go to 20
 10   continue
 20   iadj=i
      if(iadj.gt.iedge) iadj=iedge
      fscale=1./f(iadj)
      do 30 i=1,iedge
        if(i.ge.iadj) then
          f(i)=1.
        else
          f(i)=fscale*f(i)
        endif
 30   continue
      do 40 i=1,iedge
        if(y(i).ge.ym) go to 50
 40   continue
 50   istar=i+2
      fbar=scf*yplus(istar)
      fstar=f(istar)
      ystar=y(istar)
      rm1=fstar/fbar-1.
      do 60 i=1,istar
        phi=1.+rm1*(y(i)/ystar)**2
        f(i)=scf*yplus(i)*phi
 60   continue
      return
      end
