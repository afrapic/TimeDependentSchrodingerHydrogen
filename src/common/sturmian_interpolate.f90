	subroutine interpolation(nr,rgrid,xss,st1,st_int)
	use matrix_module
	use complex_module
	use datainput_module
	implicit none	
	complex*16 a(nquad),b(nquad),c(nquad),d(nquad)
	real*8 rgrid(nr),xss(nquad)
	integer l,i,j,it,ne,nr
	complex*16 st_int(0:lq,norder,nr),st1(0:lq,norder,nquad)

	st_int=c0
	do l=0,lq
		do i=1,norder
		call spline_complex(xss,st1(l,i,1:nquad),a,b,c,d,c0,c0,nquad)
			do j=1,nr
			call findi(xss,rgrid(j),nquad,it)
			if(it.gt.0)then
			st_int(l,i,j)=a(it)+rgrid(j)*(b(it)+rgrid(j)*(c(it)+rgrid(j)*d(it)))
			endif
			enddo
		enddo
	enddo

	return
	end
!----------------------------------------------------------------------------------------------------------------
	subroutine spline_complex(x,y,a,b,c,d,s1,sn,n)
!***************************************************************************
!	Calculates the Cubic-Spline interpolating coefficients
!	for a complex y(n) function evaluated in the x(n) points.
!	subroutine by A L frapiccini
!***************************************************************************
	use complex_module
	implicit none
	integer n
	real*8 x(n)
	complex*16 y(n),a(n),b(n),c(n),d(n),s1,sn
	real*8 ry(n),ra(n),rb(n),rc(n),rd(n),rs1,rsn
	real*8 zy(n),za(n),zb(n),zc(n),zd(n),zs1,zsn

	ry=real(y)
	rs1=real(s1)
	rsn=real(sn)
	call spline(x,ry,ra,rb,rc,rd,rs1,rsn,n)

	zy=aimag(y)
	zs1=aimag(s1)
	zsn=aimag(sn)
	call spline(x,zy,za,zb,zc,zd,zs1,zsn,n)

	a=ra+ci*za
	b=rb+ci*zb
	c=rc+ci*zc
	d=rd+ci*zd

	return
	end
!----------------------------------------------------------------------------------------------------------------	
      subroutine spline(x,y,a,b,c,d,s1,sn,n)
!********************************************************************************
!      cubic spline interpolation between tabulated data.
!   input:
!     x(i) (i=1, ...,n) ........ grid points.
!                    (the x values must be in increasing order).
!     y(i) (i=1, ...,n) ........ corresponding function values.
!     s1,sn ..... second derivatives at x(1) and x(n).
!            (the natural spline corresponds to taking s1=sn=0).
!     n ........................ number of grid points.
!      the interpolating polynomial in the i-th interval, from
!   x(i) to x(i+1), is
!            pi(x) = a(i)+x*(b(i)+x*(c(i)+x*d(i)))
!   output:
!     a(i),b(i),c(i),d(i) ...... spline coefficients.
!
!      ref.: m.j. maron, 'numerical analysis: a practical
!            approach', macmillan publ. co., new york 1982.
!********************************************************************************
      implicit double precision (a-h,o-z)
      dimension x(n),y(n),a(n),b(n),c(n),d(n)
      if(n.lt.4) then
      write(6,10) n
   10 format(5x,'spline interpolation cannot be performed with',i4,' points. stop.')
      stop
      endif
      n1=n-1
      n2=n-2
!  ****  auxiliary arrays h(=a) and delta(=d).
      do 1 i=1,n1
      if(x(i+1)-x(i).lt.1.0d-10) then
      write(6,11)
   11 format(5x,'spline x values not in increasing order. stop.')
      stop
      endif
      a(i)=x(i+1)-x(i)
    1 d(i)=(y(i+1)-y(i))/a(i)
!  ****  symmetric coefficient matrix (augmented).
      do 2 i=1,n2
      b(i)=2.0d0*(a(i)+a(i+1))
      k=n1-i+1
    2 d(k)=6.0d0*(d(k)-d(k-1))
      d(2)=d(2)-a(1)*s1
      d(n1)=d(n1)-a(n1)*sn
!  ****  gauss solution of the tridiagonal system.
      do 3 i=2,n2
      r=a(i)/b(i-1)
      b(i)=b(i)-r*a(i)
    3 d(i+1)=d(i+1)-r*d(i)
!  ****  the sigma coefficients are stored in array d.
      d(n1)=d(n1)/b(n2)
      do 4 i=2,n2
      k=n1-i+1
    4 d(k)=(d(k)-a(k)*d(k+1))/b(k-1)
      d(n)=sn
!  ****  spline coefficients.
      si1=s1
      do 5 i=1,n1
      si=si1
      si1=d(i+1)
      h=a(i)
      hi=1.0d0/h
      a(i)=(hi/6.0d0)*(si*x(i+1)**3-si1*x(i)**3)+hi*(y(i)*x(i+1)-y(i+1)*x(i))+(h/6.0d0)*(si1*x(i)-si*x(i+1))
      b(i)=(hi/2.0d0)*(si1*x(i)**2-si*x(i+1)**2)+hi*(y(i+1)-y(i))+(h/6.0d0)*(si-si1)
      c(i)=(hi/2.0d0)*(si*x(i+1)-si1*x(i))
    5 d(i)=(hi/6.0d0)*(si1-si)
      return
      end
!----------------------------------------------------------------------------------
      subroutine findi(x,xc,n,i)
!****************************************************************
!      finds the interval (x(i),x(i+1)) containing the value xc.
!   input:
!     x(i) (i=1, ...,n) ........ grid points.
!                    (the x values must be in increasing order).
!     xc ....................... point to be located.
!     n ........................ number of grid points.
!   output:
!     i ........................ interval index.
!****************************************************************
      implicit double precision (a-h,o-z)
     dimension x(n)
      if(xc.gt.x(n)) then
      i=-1
      return
      endif
      if(xc.lt.x(1)) then
      i=-2
      return
      endif
      i=1
      i1=n
    1 it=(i+i1)/2
      if(xc.gt.x(it)) i=it
      if(xc.le.x(it)) i1=it
      if(i1-i.gt.1) goto 1
      return
      end
