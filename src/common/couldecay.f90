	subroutine couldecay(charge,energ,r,ibound,psi,psid)
!**********************************************************************
!	Asymptote of the Coulomb wave for complex energy
!	if real(energ)>0
!	psi = exp(+/- (i k r+sigma+lq pi/2)) * (2 k r)**(-i Z /k )
!	if real(energ)<0
!	psi = exp(-k r)
!	with psid the derivative
!	Original subroutine by D. Mitnik, modified by A L Frapiccini
!**********************************************************************
	use complex_module   
	integer ibound
	real*8 r,charge
	complex*16 energ,psi,psid,kr

	!charge=(-1d0)*charge

	if(real(energ).gt.0d0)then
	kr=sqrt(2d0*energ)
	psi=ibound*ci*(kr*r+charge*log(2d0*kr*r)/kr) 
	psi=exp(psi)
	psid=ci*(kr+charge/kr/r)*psi
	endif
	if(real(energ).lt.0d0)then
	kr=sqrt(-2d0*energ)
	psi=-(kr*r)
	psi=exp(psi)
	psid=-kr*psi
	endif

	return
	end
!----------------------------------------------------------------------------
      subroutine gammaln(xxx,gammln)
!*********************************************************************************
!     Gives ln(gamma), with gamma the regular gamma function and complex argument
!*********************************************************************************
      integer s
      real*8 stp,cof(6)
      complex*16 im,gammln,xxx,xx,yy,ser,tmp
      parameter(im=(0.d0,1.d0))

      save cof,stp
      data cof,stp/76.18009172947146d0,-86.50532032941677d0,&
     	24.01409824083091d0,-1.231739572450155d0,&
       0.1208650973866179d-2,-0.5395239384953d-5,2.5066282746310005d0/
      
        xx=xxx
        yy=xx
        tmp=xx+5.5d0
        tmp=(xx+0.5d0)*log(tmp)-tmp
        ser=1.000000000190015d0
        do s=1,6
           yy=yy+1.d0
           ser=ser+cof(s)/yy
        end do
        gammln=tmp+log(stp*ser/xx)
      
        end
