	real*8 function potential(r,itype,charge,alpha,u0,rc,delta,rmax)
!************************************************************************************
!     Short range potential as a function of r
!     itype	--  0: Coulomb Potential ::::    -1 / r 
!          	--  1: Yukawa Potential  ::::    -exp(-alpha*r)/r
!          	--  2: Hulthen Potential  :::  -exp(-r/alpha)/(1-exp(-r/alpha) )
!   	        --  3: Exponential Potential  ::: -exp(-r*alpha)
!          	--  4: Box Potential :::  -1
!	All if r<rmax, 0 otherwise 
!	Function provided by D Mitnik, additions by A L Frapiccini
!*************************************************************************************
    	implicit none
    	real*8 charge,alpha,rmax,u0,rc,delta
    	real*8 r
    	integer itype
    	real*8 rexpo

	potential = 0d0

	if (itype.eq.0) then
	!.......    Coulombic Potential
 	if(r.lt.rmax) potential = -1d0 / r
	return
	endif

	if (itype.eq.1) then
	!....... Yukawa Potential
	if(r.lt.rmax) potential =  -exp(-alpha*r)/r
	return
	endif

	if (itype.eq.2) then
	!....... Hulthen Potential
	rexpo = exp(-r/alpha)
	if(r.lt.rmax) potential = -rexpo/(1d0 - rexpo)
	return
	endif

	if (itype.eq.3) then
	!....... Exponential Potential
	rexpo = exp(-r*alpha)
	if(r.lt.rmax) potential = -rexpo
	return
	endif

	if (itype.eq.4) then
	!....... Well Potential
	if (r.lt.rmax)  potential = -1d0
	return
	endif

	if (itype.eq.5) then
	!....... Fullerene Potential
	if(r.lt.rmax)then
		potential=-1d0/r
	endif
	return	
	endif

	if (itype.eq.6) then
	!....... Fullerene Potential
	if(r.lt.rmax)then
		potential=-1d0/r
	endif
	return	
	endif

	if (itype.eq.7) then
	!....... Fullerene Potential
	if(r.lt.rmax)then
		potential=-1d0/r
	endif
	return	
	endif

    	return
	end function potential
