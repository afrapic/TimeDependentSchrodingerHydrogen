	subroutine evaluate_quad_sturmians(nquad1,rmax1,xq1,wq1,rc1,delta1,xs1,rc2,delta2,xs2,rc3,delta3,xs3)
!***********************************************************************************
!	This subroutine evaluates the Legendre-Gauss quadrature weights wq1 and 
!	points xq1 to integrate the Sturmians in (0,rmax1) 
!	subroutine by Dr. A L Frapiccini
!***********************************************************************************
	implicit none
	real*8 endpts(2),rmax1,rc1,delta1,rc2,delta2,rc3,delta3
	real*8,allocatable :: work1(:),xd1(:),wd1(:)
	real*8 xq1(nquad1),wq1(nquad1),xs1(nquad1),xs2(nquad1),xs3(nquad1)
	real*8 xb,xa,xbpa,xbma
	integer j,nquad1

	!Legendre quadrature in (0,rmax1)
 	!allocate(xd1(1000),wd1(1000),work1(1000))
 	!call gaussq(1,nquad1,0.0d0,0.0d0,0,endpts,work1,xd1,wd1)
 	!deallocate(work1)
	allocate(xd1(nquad1),wd1(nquad1))
	call gauleg(-1d0,1d0,xd1,wd1,nquad1)
	xb = rmax1
	xa = 0d0
	xbpa=(xb+xa)*0.5d0
	xbma=(xb-xa)*0.5d0
	do j=1,nquad1
	xq1(j) = xbpa+xbma*xd1(j)
	enddo
 	wq1(1:nquad1)=wd1(1:nquad1)
 	deallocate(wd1) 

	xb = rc1+delta1
	xa = rc1
	xbpa=(xb+xa)*0.5d0
	xbma=(xb-xa)*0.5d0
	do j=1,nquad1
	xs1(j) = xbpa+xbma*xd1(j)
	enddo

	xb = rc2+delta2
	xa = rc2
	xbpa=(xb+xa)*0.5d0
	xbma=(xb-xa)*0.5d0
	do j=1,nquad1
	xs2(j) = xbpa+xbma*xd1(j)
	enddo

	xb = rc3+delta3
	xa = rc3
	xbpa=(xb+xa)*0.5d0
	xbma=(xb-xa)*0.5d0
	do j=1,nquad1
	xs3(j) = xbpa+xbma*xd1(j)
	enddo
	deallocate(xd1)

	return
	end subroutine
!--------------------------------------------------------------------------------------
	subroutine evaluate_quad_sturmians1(nquad1,rmin1,rmax1,xq1,wq1)
!***********************************************************************************
!	This subroutine evaluates the Legendre-Gauss quadrature weights wq1 and 
!	points xq1 to integrate the Sturmians in (rmin1,rmax1) 
!	subroutine by Dr. A L Frapiccini
!***********************************************************************************
	implicit none
	real*8 endpts(2),rmax1,rmin1
	real*8,allocatable :: work1(:),xd1(:),wd1(:)
	real*8 xq1(nquad1),wq1(nquad1)
	real*8 xb,xa,xbpa,xbma
	integer j,nquad1

	!Legendre quadrature in (0,rmax1)
 	!allocate(xd1(1000),wd1(1000),work1(1000))
 	!call gaussq(1,nquad1,0.0d0,0.0d0,0,endpts,work1,xd1,wd1)
 	!deallocate(work1)
	allocate(xd1(nquad1),wd1(nquad1))
	call gauleg(-1d0,1d0,xd1,wd1,nquad1)
	xb = rmax1
	xa = rmin1
	xbpa=(xb+xa)*0.5d0
	xbma=(xb-xa)*0.5d0
	do j=1,nquad1
	xq1(j) = xbpa+xbma*xd1(j)
	enddo
 	wq1(1:nquad1)=wd1(1:nquad1)
 	deallocate(wd1,xd1) 

	return
	end subroutine
