	subroutine evaluate_bsplines1(l,nb1,kp1,r,tx1,bb1)
!***********************************************************************************
!	This subroutine evaluates the Bsplines and their derivative in the quadrature
!	points xd.
!	subroutine by Dr. A L Frapiccini
!***********************************************************************************
	implicit none
	integer l,nb1,kp1,j
	integer left,mflag,leftmk,leftmk1,nderiv
	real*8 tx1(l+2*kp1-1),r,bb1(nb1)
	real*8,allocatable :: fn(:,:),a(:,:)
	real*8,allocatable :: tmpx(:)

	!Evaluate the Bsplines in r
	nderiv=1
	allocate(fn(kp1,nderiv),a(kp1,kp1))

      call interv(tx1,l+2*kp1-1,r,left,mflag)
      leftmk=left-kp1
      call bsplvd(tx1,kp1,r,left,a,fn,nderiv)
      leftmk1=leftmk+1
          	do j=leftmk1,left
           	bb1(j)=fn(j-leftmk,1)
		end do
	deallocate(fn,a)

	return
	end subroutine
!-----------------------------------------------------------------------
	Subroutine knots(l,kp,rmax,v)
!***********************************************************************************
!	ALL BSPLINE CALCULATION ARE FROM "A PRACTICAL GUIDE TO SPLINES" BY C. DE BOOR
!	Gives the knots {vi}, i=1,...,n+k with n=l+k-1 in the interval [0,rmax]
!	with l subdivisions from the breakpoints calculated with bps
!	Endpoints with multiplicity k
!***********************************************************************************
	implicit none
	Integer n,m,l,kp
	Real*8 v(1:l+2*kp-1),rmax
	Real*8 x(1:l+1)

	do m=1,l+1
	x(m)=rmax*((m-1d0)/(l*1d0))
	enddo
	n=l+kp-1

	do m=1,kp
	v(m)=x(1)
	enddo
	do m=1,n-kp
	v(m+kp)=x(m+1)
	enddo
	do m=1,kp
	v(n+m)=x(l+1)
	enddo

	return
	end subroutine knots

!---------------------------------------------------------------------------------------
      subroutine bsplvd ( t, k, x, left, a, dbiatx, nderiv )
!********************************************************************************************
!  	From  * a practical guide to splines *  by c. de Boor (7 may 92)    
!	calls bsplvb
!	calculates value and deriv.s of all b-splines which do not vanish at x
!
!	******  i n p u t  ******
!  		t     the knot array, of length left+k (at least)
!  		k     the order of the b-splines to be evaluated
!  		x     the point at which these values are sought
!  		left  an integer indicating the left endpoint of the interval of
!        		interest. the  k  b-splines whose support contains the interval
!               		(t(left), t(left+1))
!        		are to be considered.
!  		a s s u m p t i o n  - - -  it is assumed that
!               		t(left) .lt. t(left+1)
!        		division by zero will result otherwise (in  b s p l v b ).
!        		also, the output is as advertised only if
!               		t(left) .le. x .le. t(left+1) .
!  		nderiv   an integer indicating that values of b-splines and their
!        		derivatives up to but not including the  nderiv-th  are asked
!        		for. ( nderiv  is replaced internally by the integer  m h i g h
!        		in  (1,k)  closest to it.)
!
!		******  w o r k   a r e a  ******
!  		a     an array of order (k,k), to contain b-coeff.s of the derivat-
!        		ives of a certain order of the  k  b-splines of interest.
!
!		******  o u t p u t  ******
!  		dbiatx   an array of order (k,nderiv). its entry  (i,m)  contains
!        		value of  (m-1)st  derivative of  (left-k+i)-th  b-spline of
!        		order  k  for knot sequence  t , i=1,...,k, m=1,...,nderiv.
!************************************************************************************************
      integer k,left,nderiv,   i,ideriv,il,j,jlow,jp1mid,kp1,kp1mm ,ldummy,m,mhigh                            
      real*8 a(k,k),dbiatx(k,nderiv),t(1),x,   factor,fkp1mm,sum
      mhigh = max0(min0(nderiv,k),1)
!     mhigh is usually equal to nderiv.
      kp1 = k+1
      call bsplvb(t,kp1-mhigh,1,x,left,dbiatx)
      if (mhigh .eq. 1)                 go to 99
!     the first column of  dbiatx  always contains the b-spline values
!     for the current order. these are stored in column k+1-current
!     order  before  bsplvb  is called to put values for the next
!     higher order on top of it.
      ideriv = mhigh
      do 15 m=2,mhigh
         jp1mid = 1
         do 11 j=ideriv,k
            dbiatx(j,ideriv) = dbiatx(jp1mid,1)
   11       jp1mid = jp1mid + 1
         ideriv = ideriv - 1
         call bsplvb(t,kp1-ideriv,2,x,left,dbiatx)
   15    continue
!
!     at this point,  b(left-k+i, k+1-j)(x) is in  dbiatx(i,j) for
!     i=j,...,k and j=1,...,mhigh ('=' nderiv). in particular, the
!     first column of  dbiatx  is already in final form. to obtain cor-
!     responding derivatives of b-splines in subsequent columns, gene-
!     rate their b-repr. by differencing, then evaluate at  x.
!
      jlow = 1
      do 20 i=1,k
         do 19 j=jlow,k
   19       a(j,i) = 0.
         jlow = i
   20    a(i,i) = 1.
!     at this point, a(.,j) contains the b-coeffs for the j-th of the
!     k  b-splines of interest here.
!
      do 40 m=2,mhigh
         kp1mm = kp1 - m
         fkp1mm = float(kp1mm)
         il = left
         i = k
!
!        for j=1,...,k, construct b-coeffs of  (m-1)st  derivative of
!        b-splines from those for preceding derivative by differencing
!        and store again in  a(.,j) . the fact that  a(i,j) = 0  for
!        i .lt. j  is used.
         do 25 ldummy=1,kp1mm
            factor = fkp1mm/(t(il+kp1mm) - t(il))
!           the assumption that t(left).lt.t(left+1) makes denominator
!           in  factor  nonzero.
            do 24 j=1,i
   24          a(i,j) = (a(i,j) - a(i-1,j))*factor
            il = il - 1
   25       i = i - 1
!
!        for i=1,...,k, combine b-coeffs a(.,i) with b-spline values
!        stored in dbiatx(.,m) to get value of  (m-1)st  derivative of
!        i-th b-spline (of interest here) at  x , and store in
!        dbiatx(i,m). storage of this value over the value of a b-spline
!        of order m there is safe since the remaining b-spline derivat-
!        ives of the same order do not use this value due to the fact
!        that  a(j,i) = 0  for j .lt. i .
         do 40 i=1,k
            sum = 0.
            jlow = max0(i,m)
            do 35 j=jlow,k
   35          sum = a(j,i)*dbiatx(j,m) + sum
   40       dbiatx(i,m) = sum
   99                                   return
      end
!----------------------------------------------------------------------------------------------------
      subroutine bsplvb ( t, jhigh, index, x, left, biatx )
!**********************************************************************************
!  	From  * a practical guide to splines *  by c. de boor    
!	calculates the value of all possibly nonzero b-splines at  x  of order
!
!       	        jout  =  max( jhigh , (j+1)*(index-1) )
!
! 	 with knot sequence  t .
!
!	******  i n p u t  ******
!  	t.....knot sequence, of length  left + jout  , assumed to be nonde-
!       	 creasing.  a s s u m p t i o n . . . .
!       	                t(left)  .lt.  t(left + 1)   .
!   	d i v i s i o n  b y  z e r o  will result if  t(left) = t(left+1)
!  	jhigh,
!  	index.....integers which determine the order  jout = max(jhigh,
!       	 (j+1)*(index-1))  of the b-splines whose values at  x  are to
!       	 be returned.  index  is used to avoid recalculations when seve-
!       	 ral columns of the triangular array of b-spline values are nee-
!       	 ded (e.g., in  bsplpp  or in  bsplvd ). precisely,
!       	              if  index = 1 ,
!       	 the calculation starts from scratch and the entire triangular
!       	 array of b-spline values of orders 1,2,...,jhigh  is generated
!       	 order by order , i.e., column by column .
!                     if  index = 2 ,
!       	 only the b-spline values of order  j+1, j+2, ..., jout  are ge-
!       	 nerated, the assumption being that  biatx , j , deltal , deltar
!       	 are, on entry, as they were on exit at the previous call.
!       	    in particular, if  jhigh = 0, then  jout = j+1, i.e., just
!       	 the next column of b-spline values is generated.
!
!  	w a r n i n g . . .  the restriction   jout .le. jmax (= 20)  is im-
!       	 posed arbitrarily by the dimension statement for  deltal  and
!       	 deltar  below, but is  n o w h e r e  c h e c k e d  for .
!
!  	x.....the point at which the b-splines are to be evaluated.
!  	left.....an integer chosen (usually) so that
!                  t(left) .le. x .le. t(left+1)  .
!
!	******  o u t p u t  ******
!  	biatx.....array of length  jout , with  biatx(i)  containing the val-
!       	 ue at  x  of the polynomial of order  jout  which agrees with
!       	 the b-spline  b(left-jout+i,jout,t)  on the interval (t(left),
!       	 t(left+1)) .
!********************************************************************************************
      integer index,jhigh,left, i,j,jmax,jp1
      parameter (jmax = 20)
      real*8 biatx(jhigh),t(1),x, deltal(jmax),deltar(jmax),saved,term
!     real biatx(jhigh),t(1),x,   deltal(20),deltar(20),saved,term
!     dimension biatx(jout), t(left+jout)
!     current fortran standard makes it impossible to specify the length of
!  t  and of  biatx  precisely without the introduction of otherwise
!  superfluous additional arguments.
      data j/1/
      save j,deltal,deltar 
!
                                        go to (10,20), index
   10 j = 1
      biatx(1) = 1.
      if (j .ge. jhigh)                 go to 99
!
   20    jp1 = j + 1
         deltar(j) = t(left+j) - x
         deltal(j) = x - t(left+1-j)
         saved = 0.
         do 26 i=1,j
            term = biatx(i)/(deltar(i) + deltal(jp1-i))
            biatx(i) = saved + deltar(i)*term
   26       saved = deltal(jp1-i)*term
         biatx(jp1) = saved
         j = jp1
         if (j .lt. jhigh)              go to 20
!
!         Write(6,*)'x in sub',x
   99                                   return
      end
!--------------------------------------------------------------------------------------------
      subroutine interv ( xt, lxt, x, left, mflag )
!********************************************************************************************
!  	from  * a practical guide to splines *  by c. de Boor    
!  	computes  left = max( i :  xt(i) .lt. xt(lxt) .and.  xt(i) .le. x )  .
!
!	******  i n p u t  ******
!  	xt.....a real sequence, of length  lxt , assumed to be nondecreasing
!  	lxt.....number of terms in the sequence  xt .
!  	x.....the point whose location with respect to the sequence  xt  is
!        	to be determined.
!
!	******  o u t p u t  ******
!  	left, mflag.....both integers, whose value is
!
!   		1     -1      if               x .lt.  xt(1)
!   		i      0      if   xt(i)  .le. x .lt. xt(i+1)
!   		i      0      if   xt(i)  .lt. x .eq. xt(i+1) .eq. xt(lxt)
!   		i      1      if   xt(i)  .lt.        xt(i+1) .eq. xt(lxt) .lt. x
!
!        	In particular,  mflag = 0  is the 'usual' case.  mflag .ne. 0
!        	indicates that  x  lies outside the CLOSED interval
!       	xt(1) .le. y .le. xt(lxt) . The asymmetric treatment of the
!        	intervals is due to the decision to make all pp functions cont-
!        	inuous from the right, but, by returning  mflag = 0  even if
!        	x = xt(lxt), there is the option of having the computed pp function
!        	continuous from the left at  xt(lxt) .
!**********************************************************************************************
      integer left,lxt,mflag,ihi,ilo,istep,middle
      real*8 x,xt(lxt)
      data ilo /1/
      save ilo 
!     Write(6,*)'t in interv',(xt(i),i=1,lxt) 
      ihi = ilo + 1
      if (ihi .lt. lxt)                 go to 20
         if (x .ge. xt(lxt))            go to 110
         if (lxt .le. 1)                go to 90
         ilo = lxt - 1
         ihi = lxt
!
   20 if (x .ge. xt(ihi))               go to 40
      if (x .ge. xt(ilo))               go to 100
!
!              **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
      istep = 1
   31    ihi = ilo
         ilo = ihi - istep
         if (ilo .le. 1)                go to 35
         if (x .ge. xt(ilo))            go to 50
         istep = istep*2
                                        go to 31
   35 ilo = 1
      if (x .lt. xt(1))                 go to 90
                                        go to 50
!              **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
   40 istep = 1
   41    ilo = ihi
         ihi = ilo + istep
         if (ihi .ge. lxt)              go to 45
         if (x .lt. xt(ihi))            go to 50
         istep = istep*2
                                        go to 41
   45 if (x .ge. xt(lxt))               go to 110
      ihi = lxt
!
!           **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval.
   50 middle = (ilo + ihi)/2
      if (middle .eq. ilo)              go to 100
!     note. it is assumed that middle = ilo in case ihi = ilo+1 .
      if (x .lt. xt(middle))            go to 53
         ilo = middle
                                        go to 50
   53    ihi = middle
                                        go to 50
!**** set output and return.
   90 mflag = -1
      left = 1
                                        return
  100 mflag = 0
      left = ilo
                                        return
  110 mflag = 1
	  if (x .eq. xt(lxt)) mflag = 0
      left = lxt
  111 if (left .eq. 1)                  return
	  left = left - 1
	  if (xt(left) .lt. xt(lxt))        return
										go to 111
      end
