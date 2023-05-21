	subroutine bspline_data(i,l,tx1,kp1,rmax1,nquad_bsp1,nb1)
	use matrix_module
	integer i,l,kp1,nquad_bsp1,nb1
	real*8 tx1(l+2*kp1-1),rmax1

	call knots(l,kp1,rmax1,tx1)
	call evaluate_quad_bsplines(nquad_bsp1)
	call evaluate_bsplines(l,nb1,kp1,nquad_bsp1,tx1)

	return
	end
!----------------------------------------------------------------------------------
	subroutine evaluate_quad_bsplines(nquad_bsp)
!***********************************************************************************
!	This subroutine evaluates the Legendre-Gauss quadrature weights wd and 
!	points xd to integrate the Bsplines in the knot sequence.
!	subroutine by Dr. A L Frapiccini
!***********************************************************************************
	use matrix_module
	implicit none
	integer nquad_bsp
	real*8 endpts(2)
	real*8,allocatable :: work1(:),xd1(:),wd1(:)

 	allocate(xd1(1000),wd1(1000),work1(1000))
 	call gaussq(1,nquad_bsp,0.0d0,0.0d0,0,endpts,work1,xd1,wd1)
 	deallocate(work1)
 	allocate(xd(nquad_bsp),wd(nquad_bsp))
 	xd(1:nquad_bsp)=xd1(1:nquad_bsp)
 	wd(1:nquad_bsp)=wd1(1:nquad_bsp)	
 	deallocate(xd1,wd1)

	return
	end subroutine
!--------------------------------------------------------------------------------------------------
	subroutine evaluate_bsplines(l,nb1,kp1,nquad_bsp1,tx1)
!***********************************************************************************
!	This subroutine evaluates the Bsplines and their derivative in the quadrature
!	points xd.
!	subroutine by Dr. A L Frapiccini
!***********************************************************************************
	use matrix_module
	implicit none
	integer l,i,j,nb1,kp1,nquad_bsp1,nq
	integer left,mflag,leftmk,leftmk1,nderiv
	real*8 xb,xa,xbpa,xbma,tx1(l+2*kp1-1)
	real*8,allocatable :: fn(:,:),a(:,:)
	real*8,allocatable :: tmpx(:)

	nq=l*nquad_bsp1

	!First calculates the total grid in x for the quadrature and stores in tmpx
	allocate(tmpx(nq))
	do i=1,l
	xb = tx1(kp1+i)
	xa = tx1(kp1+i-1)
	xbpa=(xb+xa)*0.5d0
	xbma=(xb-xa)*0.5d0
	do j=1,nquad_bsp1
	tmpx((i-1)*nquad_bsp1+j) = xbpa+xbma*xd(j)
	enddo
	enddo
	deallocate(xd)
	!Saves the new quadrature in xd
	allocate(xd(nq))
	xd=tmpx
	deallocate(tmpx)

	!Evaluate the Bsplines in the quadrature points
	allocate(bb(nb1,nq),dbb(nb1,nq))
	bb=0d0
	dbb=0d0
	nderiv=2
	allocate(fn(kp1,nderiv),a(kp1,kp1))

  	do i=1,nq
        call interv(tx1,l+2*kp1-1,xd(i),left,mflag)
        leftmk=left-kp1
        call bsplvd(tx1,kp1,xd(i),left,a,fn,nderiv)
        leftmk1=leftmk+1
          	do j=leftmk1,left
           	bb(j,i)=fn(j-leftmk,1)
           	dbb(j,i)=fn(j-leftmk,2)
		end do
    	end do

	deallocate(fn,a)

	return
	end subroutine

