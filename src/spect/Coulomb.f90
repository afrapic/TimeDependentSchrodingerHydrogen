	subroutine scoul(z,e,l,r,f,fp,g,gp,err)      
!**************************************************************************       
!	This subroutine computes the sine (regular) and cosine (irregular) 
!	Coulomb wave functions .                                                            
!                                                                       
!  	input arguments:                                                     
!     	z ........ field strength, i.e. value of r*v(r) (assumed          
!                constant).                                             
!     	e ........ particle kinetic energy (positive).                    
!     	l ........ angular momentum quantum number kappa (.ne.0).         
!     	r ........ radial distance (positive).                            
!                                                                       
!  	output arguments:                                                    
!     	f, fp .... regular radial schrodinger-coulomb function and        
!                its derivative.                                        
!     	g, gp .... irregular radial schrodinger-coulomb function          
!                and its derivative.                                    
!     	err ...... accuracy of the computed functions (relative un-       
!                certainty).                                                            
!     radial functions are normalized so that, for large r, they        
!  	oscillate with unit amplitude.     
!	From package RADIAL, by F Salvat, J M Fernández-Varea and W Williamson, Jr
!****************************************************************************************                                                   
      implicit double precision (a-b,d-h,o-z), complex*16 (c)           
      common/ofcoul/delta0                                              
      common/ocoul/wavnum,eta,delta                                     
!                                                                       
      if(e.lt.0.0001d0.or.l.lt.0) then                                  
        f=0.0d0                                                         
        fp=0.0d0                                                        
        g=1.0d35                                                        
        gp=-1.0d35                                                      
        err=1.0d0                                                       
        if(e.lt.0.0001d0) write(6,2101)
 2101   format(1x,'*** error in scoul: e is too small.')
        if(l.lt.0) write(6,2102)                                        
 2102   format(1x,'*** error in scoul: l.lt.0.')                        
        return                                                          
      endif                                                             
!                                                                       
!  ****  parameters.                                                    
!                                                                       
      wavnum=dsqrt(e+e)                                                 
      if(dabs(z).gt.0.00001d0) then                                     
        eta=z/wavnum                                                    
        ical=0                                                          
      else                                                              
        eta=0.0d0                                                       
        ical=1                                                          
      endif                                                             
      rlamb=l                                                           
      x=wavnum*r                                                        
      if(ical.eq.1) goto 1                                              
!                                                                       
!  ************  coulomb functions.                                     
!                                                                       
      delta0=1.0d30                                                     
      call fcoul(eta,rlamb,x,f,fp,g,gp,err)                             
      fp=fp*wavnum                                                      
      gp=gp*wavnum                                                      
      if(delta0.lt.1.0d29) then                                         
        delta=delta0                                                    
      else                                                              
        delta=deltac(eta,rlamb)                                         
      endif                                                             
      return                                                            
!                                                                       
!  ************  z=0. spherical bessel functions.                       
!                                                                       
    1 continue                                                          
      f=x*besjn(1,l,x)                                                  
      g=-x*besjn(2,l,x)                                                 
      fp=((l+1)*besjn(1,l,x)-x*besjn(1,l+1,x))*wavnum                   
      gp=-((l+1)*besjn(2,l,x)-x*besjn(2,l+1,x))*wavnum                  
      delta=0.0d0                                                       
      err=0.0d0                                                         
      return                                                            
      end     
!-------------------------------------------------------------------------------------      
      	subroutine fcoul(eta,rlamb,x,f,fp,g,gp,err)                       
!***************************************************************************                                                                       
!	Calculation of coulomb functions for real eta, rlamb.gt.-1        
!	and x larger than, or of the order of xtp0 (the turning point        
!	for rlamb=0). steed's continued fraction method is combined          
!	with recursion relations and an asymptotic expansion. the            
!	output value err=1.0d0 indicates that the adopted evaluation         
!	algorithm is not applicable (x is too small).                        
!***************************************************************************                                                                     
      implicit double precision (a-b,d-h,o-z), complex*16 (c)           
      parameter (pi=3.1415926535897932d0,pih=0.5d0*pi,tpi=pi+pi,eps=1.0d-16,top=1.0d5,nterm=1000)                               
      common/ofcoul/delta                                               
!                                                                       
      if(rlamb.lt.-0.999d0) then                                        
        write(6,'(1x,''*** error in rcoul: rlamb.lt.-0.999'')')         
        stop                                                            
      endif                                                             
      if(x.lt.eps) goto 10                                              
!                                                                       
!  ****  numerical constants.                                           
!                                                                       
      ci=dcmplx(0.0d0,1.0d0)                                            
      ci2=2.0d0*ci                                                      
      cieta=ci*eta                                                      
      x2=x*x                                                            
      eta2=eta*eta                                                      
!                                                                       
!  ****  turning point (xtp). (44)                                      
!                                                                       
      if(rlamb.ge.0.0d0) then                                           
        xtp=eta+dsqrt(eta2+rlamb*(rlamb+1.0d0))                         
      else                                                              
        xtp=eps                                                         
      endif                                                             
      errs=10.0d0                                                       
      if(x.lt.xtp) goto 1                                               
!                                                                       
!  ************  asymptotic expansion. (71-75)                          
!                                                                       
!  ****  coulomb phase-shift.                                           
      delta=deltac(eta,rlamb)                                           
!                                                                       
      cpa=cieta-rlamb                                                   
      cpb=cieta+rlamb+1.0d0                                             
      cpz=ci2*x                                                         
      call sum2f0(cpa,cpb,cpz,c2f0,err1)                                
      cqa=cpa+1.0d0                                                     
      cqb=cpb+1.0d0                                                     
      call sum2f0(cqa,cqb,cpz,c2f0p,err2)                               
      c2f0p=ci*c2f0p*cpa*cpb/(2.0d0*x2)                                 
!  ****  functions.                                                     
      theta=x-eta*dlog(2.0d0*x)-rlamb*pih+delta                         
      if(theta.gt.1.0d4) theta=dmod(theta,tpi)                          
      ceith=cdexp(ci*theta)                                             
      cgif=c2f0*ceith                                                   
      g=cgif                                                            
      f=-ci*cgif                                                        
!  ****  derivatives.                                                   
      cgifp=(c2f0p+ci*(1.0d0-eta/x)*c2f0)*ceith                         
      gp=cgifp                                                          
      fp=-ci*cgifp                                                      
!  ****  global uncertainty. the wronskian may differ from 1 due        
!        to truncation and roundoff errors.                             
      err=dmax1(err1,err2,dabs(g*fp-f*gp-1.0d0))                        
      if(err.le.eps) return                                             
      errs=err                                                          
!                                                                       
!  ************  steed's continued fraction method.                     
!                                                                       
    1 continue                                                          
      cieta2=cieta+cieta                                                
      etax=eta*x                                                        
!                                                                       
!  ****  continued fraction for f. (60-70)                              
!                                                                       
      inull=0                                                           
      rlambn=rlamb+1.0d0                                                
      a1=-(rlambn+1.0d0)*(rlambn**2+eta2)*x/rlambn                      
      b0=(rlambn/x)+(eta/rlambn)                                        
      b1=(2.0d0*rlambn+1.0d0)*(rlambn*(rlambn+1.0d0)+etax)              
      fa3=b0                                                            
      fa2=b0*b1+a1                                                      
      fb3=1.0d0                                                         
      fb2=b1                                                            
      rf=fa3                                                            
!                                                                       
      do 2 n=2,nterm                                                    
      rfo=rf                                                            
      daf=dabs(rf)                                                      
      rlambn=rlamb+n                                                    
      an=-(rlambn**2-1.0d0)*(rlambn**2+eta2)*x2                         
      bn=(2.0d0*rlambn+1.0d0)*(rlambn*(rlambn+1.0d0)+etax)              
      fa1=fa2*bn+fa3*an                                                 
      fb1=fb2*bn+fb3*an                                                 
      tst=dabs(fb1)                                                     
!                                                                       
      if(tst.lt.1.0d-25) then                                           
        if(inull.gt.0) stop                                             
        inull=1                                                         
        fa3=fa2                                                         
        fa2=fa1                                                         
        fb3=fb2                                                         
        fb2=fb1                                                         
        rf=rfo                                                          
      else                                                              
        fa3=fa2/tst                                                     
        fa2=fa1/tst                                                     
        fb3=fb2/tst                                                     
        fb2=fb1/tst                                                     
        rf=fa2/fb2                                                      
        if(dabs(rf-rfo).lt.eps*daf) goto 3                              
      endif                                                             
    2 continue                                                          
    3 continue                                                          
      if(daf.gt.1.0d-25) then                                           
        errf=dabs(rf-rfo)/daf                                           
      else                                                              
        errf=eps                                                        
      endif                                                             
      if(errf.gt.errs) then                                             
        err=errs                                                        
        return                                                          
      endif                                                             
!                                                                       
!  ****  downward recursion for f and fp. only if rlamb.gt.1 and        
!        x.lt.xtp. (48,49)                                              
!                                                                       
      rlamb0=rlamb                                                      
      if(x.ge.xtp.or.rlamb0.lt.1.0d0) then                              
        ishift=0                                                        
      else                                                              
        ft=1.0d0                                                        
        ftp=rf                                                          
        is0=rlamb0+1.0d-6                                               
        tst=x*(x-2.0d0*eta)                                             
        do 4 i=1,is0                                                    
        etarl0=eta/rlamb0                                               
        rl=dsqrt(1.0d0+etarl0**2)                                       
        sl=(rlamb0/x)+etarl0                                            
        rlamb0=rlamb0-1.0d0                                             
        fto=ft                                                          
        ft=(sl*ft+ftp)/rl                                               
        ftp=sl*ft-rl*fto                                                
        if(ft.gt.1.0d10) then                                           
          ftp=ftp/ft                                                    
          ft=1.0d0                                                      
        endif                                                           
        rl1t=rlamb0*(rlamb0+1.0d0)                                      
        if(tst.gt.rl1t) then                                            
          ishift=i                                                      
          goto 5                                                        
        endif                                                           
    4   continue                                                        
        ishift=is0                                                      
    5   continue                                                        
        xtpc=eta+dsqrt(eta2+rl1t)                                       
        rfm=ftp/ft                                                      
      endif                                                             
!                                                                       
!  ****  continued fraction for p+ci*q with rlamb0. (76-79)             
!                                                                       
      inull=0                                                           
      can=cieta-eta2-rlamb0*(rlamb0+1.0d0)                              
      cb0=x-eta                                                         
      cbn=2.0d0*(x-eta+ci)                                              
      cfa3=cb0                                                          
      cfa2=cb0*cbn+can                                                  
      cfb3=1.0d0                                                        
      cfb2=cbn                                                          
      cpiq=cfa3                                                         
!                                                                       
      do 6 n=2,nterm                                                    
      cpiqo=cpiq                                                        
      dapiq=cdabs(cpiq)                                                 
      can=can+cieta2+(n+n-2)                                            
      cbn=cbn+ci2                                                       
      cfa1=cfa2*cbn+cfa3*can                                            
      cfb1=cfb2*cbn+cfb3*can                                            
      tst=cdabs(cfb1)                                                   
!                                                                       
      if(tst.lt.1.0d-25) then                                           
        if(inull.gt.0) stop                                             
        inull=1                                                         
        cfa3=cfa2                                                       
        cfa2=cfa1                                                       
        cfb3=cfb2                                                       
        cfb2=cfb1                                                       
        cpiq=cpiqo                                                      
      else                                                              
        cfa3=cfa2/tst                                                   
        cfa2=cfa1/tst                                                   
        cfb3=cfb2/tst                                                   
        cfb2=cfb1/tst                                                   
        cpiq=cfa2/cfb2                                                  
        if(cdabs(cpiq-cpiqo).lt.eps*dapiq) goto 7                       
      endif                                                             
    6 continue                                                          
    7 continue                                                          
      if(dapiq.gt.1.0d-25) then                                         
        errpiq=cdabs(cpiq-cpiqo)/dapiq                                  
      else                                                              
        errpiq=eps                                                      
      endif                                                             
      if(errpiq.gt.errs) then                                           
        err=errs                                                        
        return                                                          
      endif                                                             
      cpiq=ci*cpiq/x                                                    
!                                                                       
      rp=cpiq                                                           
      rq=-ci*cpiq                                                       
      if(rq.le.1.0d-25) goto 10                                         
      err=dmax1(errf,errpiq)                                            
!                                                                       
!  ****  inverting steed's transformation. (57,58)                      
!                                                                       
      if(ishift.lt.1) then                                              
        rfp=rf-rp                                                       
        f=dsqrt(rq/(rfp**2+rq**2))                                      
        if(fb2.lt.0.0d0) f=-f                                           
        fp=rf*f                                                         
        g=rfp*f/rq                                                      
        gp=(rp*rfp-rq**2)*f/rq                                          
        if(x.lt.xtp.and.g.gt.top*f) goto 10                             
      else                                                              
        rfp=rfm-rp                                                      
        fm=dsqrt(rq/(rfp**2+rq**2))                                     
        g=rfp*fm/rq                                                     
        gp=(rp*rfp-rq**2)*fm/rq                                         
        if(x.lt.xtpc.and.g.gt.top*fm) goto 10                           
!  ****  upward recursion for g and gp (if ishift.gt.0). (50,51)        
        do 8 i=1,ishift                                                 
        rlamb0=rlamb0+1.0d0                                             
        etarl0=eta/rlamb0                                               
        rl=dsqrt(1.0d0+etarl0**2)                                       
        sl=(rlamb0/x)+etarl0                                            
        go=g                                                            
        g=(sl*go-gp)/rl                                                 
        gp=rl*go-sl*g                                                   
        if(g.gt.1.0d35) goto 10                                         
    8   continue                                                        
    9   w=rf*g-gp                                                       
        f=1.0d0/w                                                       
        fp=rf/w                                                         
      endif                                                             
!  ****  the wronskian may differ from 1 due to roundoff errors.        
      err=dmax1(err,dabs(fp*g-f*gp-1.0d0))                              
      return                                                            
!                                                                       
   10 f=0.0d0                                                           
      fp=0.0d0                                                          
      g=1.0d35                                                          
      gp=-1.0d35                                                        
      err=1.0d0                                                         
      return                                                            
      end   
!-------------------------------------------------------------------------------  
	subroutine sum2f0(ca,cb,cz,cf,err)                                
!***********************************************************************                                                                 
!     	summation of the 2f0(ca,cb;cs) hypergeometric asymptotic          
!  	series. the positive and negative contributions to the real          
!  	and imaginary parts are added separately to obtain an estimate       
!  	of rounding errors.                                                  
!***********************************************************************                                                                       
      implicit double precision (a-b,d-h,o-z), complex*16 (c)           
      parameter (eps=1.0d-16,accur=0.5d-15,nterm=75)                    
      rrp=1.0d0                                                         
      rrn=0.0d0                                                         
      rip=0.0d0                                                         
      rin=0.0d0                                                         
      cdf=1.0d0                                                         
      err2=0.0d0                                                        
      err3=1.0d0                                                        
      do 1 i=1,nterm                                                    
      j=i-1                                                             
      cdf=cdf*(ca+j)*(cb+j)/(i*cz)                                      
      err1=err2                                                         
      err2=err3                                                         
      err3=cdabs(cdf)                                                   
      if(err1.gt.err2.and.err2.lt.err3) goto 2                          
      ar=cdf                                                            
      if(ar.gt.0.0d0) then                                              
        rrp=rrp+ar                                                      
      else                                                              
        rrn=rrn+ar                                                      
      endif                                                             
      ai=dcmplx(0.0d0,-1.0d0)*cdf                                       
      if(ai.gt.0.0d0) then                                              
        rip=rip+ai                                                      
      else                                                              
        rin=rin+ai                                                      
      endif                                                             
      cf=dcmplx(rrp+rrn,rip+rin)                                        
      af=cdabs(cf)                                                      
      if(af.gt.1.0d25) then                                             
        cf=0.0d0                                                        
        err=1.0d0                                                       
        return                                                          
      endif                                                             
      if(err3.lt.1.0d-25*af.or.err3.lt.eps) then                        
         err=eps                                                        
         return                                                         
      endif                                                             
    1 continue                                                          
!  ****  roundoff error.                                                
    2 continue                                                          
      tr=dabs(rrp+rrn)                                                  
      if(tr.gt.1.0d-25) then                                            
        errr=(rrp-rrn)*accur/tr                                         
      else                                                              
        errr=1.0d0                                                      
      endif                                                             
      ti=dabs(rip+rin)                                                  
      if(ti.gt.1.0d-25) then                                            
        erri=(rip-rin)*accur/ti                                         
      else                                                              
        erri=1.0d0                                                      
      endif                                                             
!  ****  ... and truncation error.                                      
      if(ar.gt.1.0d-25) then                                            
      err=dmax1(errr,erri)+err2/af                                      
      else                                                              
      err=dmax1(errr,erri)                                              
      endif                                                             
      return                                                            
      end                                         
!---------------------------------------------------------------------------------------                            
      function deltac(eta,rlamb)                                        
!***********************************************************************                                                               
!     calculation of coulomb phase shift (modulus 2*pi). (47)           
!***********************************************************************                                                                       
      implicit double precision (a-b,d-h,o-z), complex*16 (c)           
      parameter (pi=3.1415926535897932d0,tpi=pi+pi)                     
      ci=dcmplx(0.0d0,1.0d0)                                            
!  ****  coulomb phase-shift.                                           
      deltac=-ci*clgam(rlamb+1.0d0+ci*eta)                              
      if(deltac.ge.0.0d0) then                                          
        deltac=dmod(deltac,tpi)                                         
      else                                                              
        deltac=-dmod(-deltac,tpi)                                       
      endif                                                             
      return                                                            
      end                                                                      
!--------------------------------------------------------------------------
      	function clgam(cz)                                                
!***********************************************************************                                                                       
!     	this function gives log(gamma(cz)) for complex arguments.         
!                                                                       
!   	ref.: m. abramowitz and i.a. stegun, 'handbook of mathemati-        
!         cal functions'. dover, new york (1974). pp 255-257.           
!***********************************************************************                                                                       
      implicit double precision (a-b,d-h,o-z), complex*16 (c)           
      parameter (pi=3.1415926535897932d0)                               
      cza=cz                                                            
      iconj=0                                                           
      ar=cza                                                            
      clgam=36.84136149d0                                               
      if(cdabs(cza).lt.1.0d-16) return                                  
!                                                                       
      ai=cza*dcmplx(0.0d0,-1.0d0)                                       
      if(ai.gt.0.0d0) then                                              
        iconj=0                                                         
      else                                                              
        iconj=1                                                         
        cza=dconjg(cza)                                                 
      endif                                                             
!                                                                       
      czfac=1.0d0                                                       
      czfl=0.0d0                                                        
    1 czfac=czfac/cza                                                   
      if(cdabs(czfac).gt.1.0d8) then                                    
        czfl=czfl+cdlog(czfac)                                          
        czfac=1.0d0                                                     
      endif                                                             
      cza=cza+1.0d0                                                     
      ar=cza                                                            
      if(cdabs(cza).lt.1.0d-16) return                                  
      if(cdabs(cza).gt.15.0d0.and.ar.gt.0.0d0) goto 2                   
      goto 1                                                            
!  ****  stirling's expansion of cdlog(gamma(cza)).                     
    2 czi2=1.0d0/(cza*cza)                                              
      czs=(43867.0d0/244188.0d0)*czi2                                   
      czs=(czs-3617.0d0/122400.0d0)*czi2                                
      czs=(czs+1.0d0/156.0d0)*czi2                                      
      czs=(czs-691.0d0/360360.0d0)*czi2                                 
      czs=(czs+1.0d0/1188.0d0)*czi2                                     
      czs=(czs-1.0d0/1680.0d0)*czi2                                     
      czs=(czs+1.0d0/1260.0d0)*czi2                                     
      czs=(czs-1.0d0/360.0d0)*czi2                                      
      czs=(czs+1.0d0/12.0d0)/cza                                        
      clgam=(cza-0.5d0)*cdlog(cza)-cza+9.1893853320467274d-1+czs+czfl+cdlog(czfac)                                           
      if(iconj.eq.1) clgam=dconjg(clgam)                                
      return                                                            
      end                                                               
!--------------------------------------------------------------------------------------      
      function besjn(jy,n,x)                                            
!***********************************************************************                                                                     
!      this function computes the spherical bessel functions of         
!   the first kind and spherical bessel functions of the second         
!   kind (also known as spherical neumann functions) for real           
!   positive arguments.                                                 
!***********************************************************************                                
      implicit double precision (a-h,o-z)                               
      parameter (pi=3.1415926535897932d0,tpi=pi+pi)                     
      if(x.lt.0) then                                                   
        write(6,1000)                                                   
 1000   format(1x,'*** negative argument in function besjn.')           
        stop                                                            
      endif                                                             
!  ****  order and phase correction for neumann functions.              
!        abramowitz and stegun, eq. 10.1.15.                            
      if(jy.eq.2) then                                                  
        nl=-n-1                                                         
        iph=2*mod(iabs(n),2)-1                                          
      else                                                              
        nl=n                                                            
        iph=1                                                           
      endif                                                             
!  ****  selection of calculation mode.                                 
      if(nl.lt.0) goto 10                                               
      if(x.gt.1.0d0*nl) goto 7                                          
      xi=x*x                                                            
      if(xi.gt.nl+nl+3.0d0) goto 4                                      
!  ****  power series for small arguments and positive orders.          
!        abramowitz and stegun, eq. 10.1.2.                             
      f1=1.0d0                                                          
      ip=1                                                              
      if(nl.ne.0) then                                                  
        do 1 i=1,nl                                                     
        ip=ip+2                                                         
    1   f1=f1*x/ip                                                      
      endif                                                             
      xi=0.5d0*xi                                                       
      besjn=1.0d0                                                       
      ps=1.0d0                                                          
      do 2 i=1,500                                                      
      ip=ip+2                                                           
      ps=-ps*xi/(i*ip)                                                  
      besjn=besjn+ps                                                    
      if(dabs(ps).lt.1.0d-18*dabs(besjn)) goto 3                        
    2 continue                                                          
    3 besjn=iph*f1*besjn                                                
      return                                                            
!  ****  miller's method for positive orders and intermediate           
!        arguments. abramowitz and stegun, eq. 10.1.19.                 
    4 xi=1.0d0/x                                                        
      f2=0.0d0                                                          
      f3=1.0d-35                                                        
      ip=2*(nl+31)+3                                                    
      do 5 i=1,31                                                       
      f1=f2                                                             
      f2=f3                                                             
      ip=ip-2                                                           
      f3=ip*xi*f2-f1                                                    
      if(dabs(f3).gt.1.0d30) then                                       
        f2=f2/f3                                                        
        f3=1.0d0                                                        
      endif                                                             
    5 continue                                                          
      besjn=1.0d0                                                       
      f2=f2/f3                                                          
      f3=1.0d0                                                          
      do 6 i=1,nl                                                       
      f1=f2                                                             
      f2=f3                                                             
      ip=ip-2                                                           
      f3=ip*xi*f2-f1                                                    
      if(dabs(f3).gt.1.0d30) then                                       
        besjn=besjn/f3                                                  
        f2=f2/f3                                                        
        f3=1.0d0                                                        
      endif                                                             
    6 continue                                                          
      besjn=iph*xi*dsin(x)*besjn/f3                                     
      return                                                            
!  ****  recurrence relation for arguments greater than order.          
!        abramowitz and stegun, eq. 10.1.19.                            
    7 xi=1.0d0/x                                                        
      f3=xi*dsin(x)                                                     
      if(nl.eq.0) goto 9                                                
      f2=f3                                                             
      f3=xi*(f2-dcos(x))                                                
      if(nl.eq.1) goto 9                                                
      ip=1                                                              
      do 8 i=2,nl                                                       
      f1=f2                                                             
      f2=f3                                                             
      ip=ip+2                                                           
    8 f3=ip*xi*f2-f1                                                    
    9 besjn=iph*f3                                                      
      return                                                            
!  ****  recurrence relation for negative orders.                       
!        abramowitz and stegun, eq. 10.1.19.                            
   10 nl=iabs(nl)                                                       
      if(x.lt.7.36d-1*(nl+1)*1.0d-35**(1.0d0/(nl+1))) then              
        besjn=-1.0d35                                                   
        return                                                          
      endif                                                             
      xi=1.0d0/x                                                        
      f3=xi*dsin(x)                                                     
      f2=xi*(f3-dcos(x))                                                
      ip=3                                                              
      do 11 i=1,nl                                                      
      f1=f2                                                             
      f2=f3                                                             
      ip=ip-2                                                           
      f3=ip*xi*f2-f1                                                    
      if(dabs(f3).gt.1.0d35) then                                       
        besjn=-1.0d35                                                   
        return                                                          
      endif                                                             
   11 continue                                                          
      besjn=iph*f3                                                      
      return                                                            
      end                                                                                                                    
