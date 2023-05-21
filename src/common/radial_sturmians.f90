	subroutine radial_sturmians(l,kp,nb,tx,rmax,norder,r,coef,nderiv,func)
!*********************************************************************************************
!	Sturmian function from 1 to norder evaluated in r normalized as <S(ord)|S(ord)>=1
!	Output in func
!	subroutine by A L Frapiccini
!*********************************************************************************************
	use complex_module
	integer l,kp,nb,norder,nd
	real*8 tx(l+2*kp-1),rmax
	real*8 r,bb1(nderiv,nb)
	complex*16 func(nderiv,1:norder),coef(1:nb,1:norder)
	integer left,mflag,leftmk,leftmk1,nderiv
	real*8,allocatable :: fn(:,:),a(:,:)

	if(r.gt.rmax)then
	func=c0
	return
	endif

	!Evaluate the Bsplines in r
	bb1=0d0
	allocate(fn(kp,nderiv),a(kp,kp))
	call interv(tx,l+2*kp-1,r,left,mflag)
	leftmk=left-kp
	call bsplvd(tx,kp,r,left,a,fn,nderiv)
	leftmk1=leftmk+1
	do nd=1,nderiv
          	do j=leftmk1,left
           	bb1(nd,j)=fn(j-leftmk,nd)
		end do
    	end do
	deallocate(fn,a)

	do nd=1,nderiv
	do i=1,norder
	func(nd,i)=c0	
		do j=1,nb
		func(nd,i)=func(nd,i)+coef(j,i)*bb1(nd,j)
		enddo
	enddo
	enddo

	return
	end subroutine 
