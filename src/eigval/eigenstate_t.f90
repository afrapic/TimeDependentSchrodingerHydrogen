	subroutine eigenstate_t(h_atomic,m_over,h_int_up,h_int_low,tsave,isave)
	use datainput_module
	use file_io_module
	use scratchdir_module
	use matrix_module
	use complex_module
	implicit none
	integer nn,info,i,j,l,ir,isave
	real*8 h_atomic(ntot1,norder),m_over(ntot1,norder),h_int_low(ntot2,norder),h_int_up(ntot2,norder) 
	real*8 tsave
	complex*16 gt,vl,wp,wpl(0:1)
	complex*16,allocatable :: htot(:,:),mtot(:,:),bet(:),vr(:,:),work(:)
	real*8,allocatable :: rwork(:),eigvals(:),eigvalp(:),alph(:)
	integer,allocatable :: indx(:)

	call pulse_time(tsave,gt)
	print*,gt
	write(6,*)

	nn=2*norder
	allocate(htot(nn,nn),mtot(nn,nn))
	htot=c0
	mtot=c0
	htot(1:norder,1:norder)=c1*h_atomic(1:norder,1:norder)
	htot(norder+1:nn,norder+1:nn)=c1*h_atomic(norder+1:nn,1:norder)
	htot(norder+1:nn,1:norder)=gt*h_int_up(1:norder,1:norder)
	htot(1:norder,norder+1:nn)=conjg(transpose(gt*h_int_up(1:norder,1:norder)))
	
	mtot(1:norder,1:norder)=c1*m_over(1:norder,1:norder)
	mtot(norder+1:nn,norder+1:nn)=c1*m_over(norder+1:nn,1:norder)

	allocate(alph(nn),work(2*nn),rwork(3*nn-2))
	call zhegv(1,'V','U',nn,htot,nn,mtot,nn,alph,work,2*nn,rwork,info)
	deallocate(mtot,work,rwork)
	if(info.eq.0)then
      write(99,*) alph(1)
	deallocate(alph)
	write(6,*)
	endif

	j=1
	do ir=1,nquad
	wp=c0
		do l=0,1
		wpl(l)=c0
			do i=1,norder
			wpl(l)=wpl(l)+htot(l*norder+i,j)*st(l,i,ir)
			enddo
		wp=wp+abs(wpl(l))**2
		enddo
	write(100+isave,'(f20.10,2es30.15)') xq(ir),real(wp),aimag(wp)
	enddo
	deallocate(htot)	

	return
	end
!--------------------------------------------------------------------------------------------------------
!	subroutine pulse_time(t,gt0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!	The time dependency of the interaction with the pulse is written as
!!!!	H(int)=g(t)*mat_int*constants
!!!!	with g(t)=	A(t) for velocity gauge
!!!!			E(t)=-dA(t)/dt for length gauge
!!!!	Then constants are added in the multiplication routines.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	use datainput_module
!	use pi_module
!	use complex_module
!	implicit none   
!	real*8 t,t1,t2
!	complex*16 gt0,cte

!	gt0=c0
!	if(gauge.eq.1)then
!	cte=cmplx(amax)
!	endif
!	if(gauge.eq.2)then
!	cte=-ci*amax
!	endif

!	if(env.eq.1)then
!        gt0=cte*sin(omega*t)*sin(pi*t/tau)**2
!	endif

!	if(env.eq.2)then
!		if(t.ge.tini1.and.t.le.tfin1)then
!		gt0=cte*sin(omega*t)*sin(pi*t/tau)**2
!		endif
!		if(t.ge.tini2.and.t.le.tfin2)then
!		gt0=cte*sin(omega1*(t-tfin1-tsep))*sin(pi*(t-tfin1-tsep)/tau1)**2
!		endif
!	endif

!	if(env.eq.3)then
!		if(t.ge.tini1.and.t.le.tau)then
!		gt0=cte*sin(omega*t)*sin(pi*t/2d0/tau)**2
!		endif
!		if(t.gt.tau.and.t.le.tau+tau1)then
!		gt0=cte*sin(omega*t)
!		endif
!		if(t.gt.tau+tau1.and.t.le.tau1+2*tau)then
!		gt0=cte*sin(omega*t)*sin(pi*(t-tau1)/2d0/tau)**2
!		endif
!	endif

!      	return
!      	end
