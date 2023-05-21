!---------------------------------------------------------------------------------
	subroutine integ_st_1s(stl,integ)
	use datainput_module
	use matrix_module
	use complex_module
	implicit none
	integer i,j,jp
	real*8 integ(norder)
	real*8 a0,vr
	complex*16 stl(norder,nquad),res

	a0=rmax/2d0
	integ=0d0
	do j=1,norder
		res=c0
		do i=1,nquad
		vr=2d0*xq(i)*exp(-xq(i))
		res=res+wq(i)*stl(j,i)*vr*a0
		enddo
	integ(j)=real(res)
	enddo

	return
	end


!---------------------------------------------------------------------------------
	subroutine integ_st_pot(stl,integ)
	use datainput_module
	use matrix_module
	use complex_module
	implicit none
	integer i,j,jp
	complex*16 integ(norder,norder)
	real*8 a0,vr
	real*8,external :: potential
	complex*16 stl(norder,nquad),res

	a0=rmax/2d0
	integ=c0
	do j=1,norder
		do jp=j,norder
		res=c0
			do i=1,nquad
			vr=potential(xq(i),ishort,charge,alpha,1d0,1d0,1d0,rmax)
			res=res+wq(i)*stl(j,i)*stl(jp,i)*vr*a0
			enddo
		integ(jp,j)=res
		integ(j,jp)=integ(jp,j)
		enddo
	enddo

	return
	end
!-----------------------------------------------------------------------------------------
	subroutine integ_st_over(stl,integ)
	use datainput_module
	use matrix_module
	use complex_module
	use Complex_Incomplete_Gamma
	implicit none
	integer i,j,jp
	complex*16 integ(norder,norder)
	real*8 a0
	complex*16 stl(norder,nquad),res,asymp,eta,nu

	a0=rmax/2d0
	integ=c0
	do j=1,norder
		do jp=j,norder
		res=c0
			do i=1,nquad
			res=res+wq(i)*stl(j,i)*stl(jp,i)*a0
			enddo
		integ(jp,j)=res
		integ(j,jp)=integ(jp,j)
		enddo
	enddo

	return
	end
!---------------------------------------------------------------------------------
	subroutine integ_st_vw0(stl,integ)
	use datainput_module
	use matrix_module
	use complex_module
	implicit none
	integer i,j,jp
	complex*16 integ(norder,norder)
	real*8 a0,ra,rb,vr
	complex*16 stl(norder,nquad),res

	rb=rc0+delta0
	ra=rc0
	a0=(rb-ra)/2d0
	integ=c0
	do j=1,norder
		do jp=j,norder
		res=c0
			do i=1,nquad
			vr=-u0
			res=res+wq(i)*stl(j,i)*stl(jp,i)*a0*vr
			enddo
		integ(jp,j)=res
		integ(j,jp)=integ(jp,j)
		enddo
	enddo

	return
	end
!---------------------------------------------------------------------------------
	subroutine integ_st_vw1(stl,integ)
	use datainput_module
	use matrix_module
	use complex_module
	implicit none
	integer i,j,jp
	complex*16 integ(norder,norder)
	real*8 a0,ra,rb,vr
	complex*16 stl(norder,nquad),res

	rb=rc1+delta1
	ra=rc1
	a0=(rb-ra)/2d0
	integ=c0
	do j=1,norder
		do jp=j,norder
		res=c0
			do i=1,nquad
			vr=-u1
			res=res+wq(i)*stl(j,i)*stl(jp,i)*a0*vr
			enddo
		integ(jp,j)=res
		integ(j,jp)=integ(jp,j)
		enddo
	enddo

	return
	end
!---------------------------------------------------------------------------------
	subroutine integ_st_vw2(stl,integ)
	use datainput_module
	use matrix_module
	use complex_module
	implicit none
	integer i,j,jp
	complex*16 integ(norder,norder)
	real*8 a0,ra,rb,vr
	complex*16 stl(norder,nquad),res

	rb=rc2+delta2
	ra=rc2
	a0=(rb-ra)/2d0
	integ=c0
	do j=1,norder
		do jp=j,norder
		res=c0
			do i=1,nquad
			vr=-u2
			res=res+wq(i)*stl(j,i)*stl(jp,i)*a0*vr
			enddo
		integ(jp,j)=res
		integ(j,jp)=integ(jp,j)
		enddo
	enddo

	return
	end
!-----------------------------------------------------------------------------------------
	subroutine integ_st_dr(dstl,stlp,integ)
	use datainput_module
	use matrix_module
	use complex_module
	implicit none
	integer i,j,jp
	complex*16 integ(norder,norder)
	real*8 a0
	complex*16 dstl(norder,nquad),stlp(norder,nquad),res

	a0=rmax/2d0
	integ=c0
	do j=1,norder
		do jp=1,norder
		res=c0
			do i=1,nquad
			res=res+wq(i)*dstl(j,i)*stlp(jp,i)*a0
			enddo
		integ(jp,j)=res
		enddo
	enddo

	return
	end
!-----------------------------------------------------------------------------------------
	subroutine integ_st_rq(q,stl,stlp,integ)
	use datainput_module
	use matrix_module
	use complex_module
	use Complex_Incomplete_Gamma
	implicit none
	integer i,j,jp,q
	complex*16 integ(norder,norder)
	real*8 a0,r
	complex*16 stl(norder,nquad),stlp(norder,nquad),res,asymp,eta,nu

	a0=rmax/2d0
	integ=c0
	do j=1,norder
		do jp=1,norder
		res=c0
			do i=1,nquad
			r=xq(i)
			res=res+wq(i)*stl(j,i)*stlp(jp,i)*r**(q)*a0
			enddo
		integ(jp,j)=res
		enddo
	enddo

	return
	end
!-----------------------------------------------------------------------------------------------
	subroutine integ_st_long(stl,integ)
	use datainput_module
	use matrix_module
	use complex_module
	implicit none
	integer i,j,jp
	complex*16 integ(norder,norder)
	real*8 a0,vr
	real*8,external :: potential
	complex*16 stl(norder,nquad),res

	a0=rmax/2d0
	integ=c0
	do j=1,norder
		do jp=j,norder
		res=c0
			do i=1,nquad
			vr=potential(xq(i),ilong,charge,alpha,u0,rc0,delta0,rmax)
			res=res+wq(i)*stl(j,i)*stl(jp,i)*vr*a0
			enddo
		integ(jp,j)=res
		integ(j,jp)=integ(jp,j)
		enddo
	enddo

	return
	end
!-------------------------------------------------------------------------------------------------------
	subroutine groundstate(h_atomic,m_over,yt)
	use datainput_module
	use complex_module
      	use scratchdir_module
      	use file_io_module
	use matrix_module
	implicit none
	real*8 h_atomic(ntot1,norder),m_over(ntot1,norder)
	complex*16 yt(ntot1),wp
	real*8 cwork
	real*8 vr,nrm
	integer info,lwork
	real*8,allocatable :: work(:),ma(:,:),mb(:,:),eigenvalues(:),rhs(:)
	integer nn,ll,ele,i,j,ir
	character*2 cl,cn
	real*8,external :: potential
	integer,allocatable :: ipiv(:)
	complex*16,allocatable :: stl(:,:)
	
	if(qnumber(1).ne.0)then
	ll=qnumber(2)
	nn=qnumber(1)-ll
	ele=ll

	allocate(ma(norder,norder),mb(norder,norder),eigenvalues(norder))
	ma(:,:)=h_atomic(ele*norder+1:(ele+1)*norder,:)
	mb(:,:)=m_over(ele*norder+1:(ele+1)*norder,:)

	lwork=-1
	call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,cwork,lwork,info)

	lwork=int(cwork)
	allocate(work(lwork))
	call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,work,lwork,info)
		if(info.ne.0)then
		write(6,*) 'info:',info,'in dsygv in groundstate'
		endif
        deallocate(work)

	write(6,*) 'Bound state energy (nl):',nn,ll,eigenvalues(nn)
	call write_matrix(eigenvalues,trim(dir_output)//'eigval-'//trim(file_output)//'.dat')

	yt=c0
	yt(ele*norder+1:(ele+1)*norder)=cmplx(ma(:,nn))
	call norm2(m_over,yt,nrm)
	write(6,*) 'Norm:',nrm
	deallocate(ma,mb,eigenvalues)

!!
	else
	allocate(mb(norder,norder))
	ele=0
	mb(:,:)=m_over(ele*norder+1:(ele+1)*norder,:)
	allocate(rhs(norder),stl(norder,nquad))
	stl(:,:)=st(0,:,:)
	call integ_st_1s(stl,rhs)
	deallocate(stl)

	allocate(ipiv(norder))
	lwork=-1
	call dsysv('U',norder,1,mb,norder,ipiv,rhs,norder,cwork,lwork,info)

	lwork=int(cwork)
	allocate(work(lwork))
	call dsysv('U',norder,1,mb,norder,ipiv,rhs,norder,work,lwork,info)
		if(info.ne.0)then
		write(6,*) 'info:',info,'in dsysv in groundstate'
		endif
      deallocate(work,mb,ipiv)

	yt=c0
	yt(ele*norder+1:(ele+1)*norder)=cmplx(rhs(:))
	call norm2(m_over,yt,nrm)
	write(6,*) 'Norm:',nrm
	deallocate(rhs)
	endif

	return
	end
!-------------------------------------------------------------------------------------------------------
	subroutine boundstate(h_atomic,m_over)
	use datainput_module
	use complex_module
      	use scratchdir_module
      	use file_io_module
	use matrix_module
	implicit none
	real*8 h_atomic(ntot1,norder),m_over(ntot1,norder)
	complex*16 yt(ntot1),wp
	real*8 ma(norder,norder),cwork,mb(norder,norder)
	real*8 eigenvalues(norder),vr,nrm
	integer info,lwork
	real*8,allocatable :: work(:)
	integer nn,ll,ele,i,j,ir
	character*2 cl,cn
	real*8,external :: potential

	do ele=0,1
	write(cl,'(i2.2)') ele

	ma(:,:)=h_atomic(ele*norder+1:(ele+1)*norder,:)
	mb(:,:)=m_over(ele*norder+1:(ele+1)*norder,:)

	lwork=-1
	call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,cwork,lwork,info)

	lwork=int(cwork)
	allocate(work(lwork))
	call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,work,lwork,info)
		if(info.ne.0)then
		write(6,*) 'info:',info,'in dsygv in groundstate'
		endif
        	deallocate(work)

		do nn=1,3
		write(cn,'(i2.2)') nn
		yt=c0
		yt(ele*norder+1:(ele+1)*norder)=cmplx(ma(:,nn))
		call write_matrix_raw(yt,trim(dir_scratch)//'yt-n'//cn//'l'//cl//'-'//trim(file_output)//'.dat')
		call plot_bstate(nn,ele,yt)
		enddo
	enddo


	return
	end
!-------------------------------------------------------------------------------------------------------
	subroutine tdep_eigenstate(h_atomic,m_over)
	use datainput_module
	use complex_module
      use scratchdir_module
      use file_io_module
	use matrix_module
	implicit none
	real*8 h_atomic(ntot1,norder),m_over(ntot1,norder)
	complex*16,allocatable :: ma(:,:),mb(:,:),work(:)
	real*8,allocatable :: h_int_up(:,:)
	complex*16 gt0,cwork
	real*8 t,dt
	real*8,allocatable :: rwork(:),eigval(:)
	integer lwork,info,i

	dt=tfin_pulse/dble(100-1)
	open(unit=5000,file=trim(dir_output)//'eigval_tdep-'//trim(file_output)//'.dat')

	do i=1,100
	t=(i-1d0)*dt
	print*,t
	call pulse_time(t,gt0)

	allocate(h_int_up(ntot2,norder))
	call read_matrix_raw(h_int_up,trim(dir_scratch)//'hintup-'//trim(file_output)//'.dat')
	close(1000)
	
	allocate(ma(2*norder,2*norder),mb(2*norder,2*norder))
	ma=c0
	mb=c0
	ma(1:norder,1:norder)=c1*h_atomic(1:norder,:)
	ma(norder+1:2*norder,norder+1:2*norder)=c1*h_atomic(norder+1:2*norder,:)
	ma(norder+1:2*norder,1:norder)=gt0*h_int_up(1:norder,1:norder)
	mb(1:norder,1:norder)=c1*m_over(1:norder,:)
	mb(norder+1:2*norder,norder+1:2*norder)=c1*m_over(norder+1:2*norder,:)
	deallocate(h_int_up)

	lwork=-1
	allocate(eigval(2*norder),rwork(3*2*norder-2))
	call zhegv(1,'N','L',norder*2,ma,norder*2,mb,norder*2,eigval,cwork,lwork,rwork,info)

	lwork=int(cwork)
	allocate(work(lwork))
	call zhegv(1,'N','L',norder*2,ma,norder*2,mb,norder*2,eigval,work,lwork,rwork,info)
		if(info.ne.0)then
		write(6,*) 'info:',info,'in zhegv in tdep_eigensate'
		endif
      deallocate(work,rwork,mb)

	write(5000,'(51f20.10)') t,eigval(1:50)

	deallocate(ma,eigval)
	enddo
	close(5000)

	return
	end
!-------------------------------------------------------------------------------------------
	complex*16 function scalar_conjg(x,y)
	use complex_module
	use datainput_module
	implicit none
	complex*16 x(ntot1),y(ntot1)
	complex*16 sum,sump
	integer l,j

	sum=c0
	do l=0,lq
	sump=dot_product(x(l*norder+1:(l+1)*norder),y(l*norder+1:(l+1)*norder))
	sum=sum+sump
	enddo

	scalar_conjg=sum

	return
	end
!---------------------------------------------------------------
	subroutine hatomic_vec(h_atomic,vec,x)
	use complex_module
	use datainput_module
	implicit none
	real*8 h_atomic(ntot1,norder)
	complex*16 vec(ntot1),x(ntot1)
	integer l,j,jp

	x=c0
	do l=0,lq
	call my_zgemv('N',norder,norder,c1,h_atomic(l*norder+1:(l+1)*norder,:),norder,&
			 vec(l*norder+1:(l+1)*norder),1,c0,x(l*norder+1:(l+1)*norder),1)
	enddo

	return
	end
!-----------------------------------------------------------------------------------------------
	subroutine hint_vec(gt,h_int_up,vec,x)
	use complex_module
	use datainput_module
	implicit none
	real*8 h_int_up(ntot2,norder)
	complex*16 vec(ntot1),x(ntot1),gt
	complex*16,allocatable :: tmp(:)
	integer l,j,jp,nil,nfl,nil1,nfl1

	x=c0
	allocate(tmp(ntot1))
	tmp=c0

	do l=0,lq
		if(l.eq.0)then
		nil=l*norder+1
		nfl=(l+1)*norder
		nil1=(l+1)*norder+1
		nfl1=(l+2)*norder
		call my_zgemv('N',norder,norder,gt,h_int_up(nil:nfl,:),norder,vec(nil1:nfl1),1,c0,x(nil:nfl),1)
		endif

		if(l.eq.lq)then
		nil=(l-1)*norder+1
		nfl=(l)*norder
		nil1=l*norder+1
		nfl1=(l+1)*norder
		call my_zgemv('C',norder,norder,gt,h_int_up(nil:nfl,:),norder,vec(nil:nfl),1,c0,x(nil1:nfl1),1)
		endif

		if(l.gt.0.and.l.lt.lq)then
		nil=l*norder+1
		nfl=(l+1)*norder
		nil1=(l+1)*norder+1
		nfl1=(l+2)*norder
		call my_zgemv('N',norder,norder,gt,h_int_up(nil:nfl,:),norder,vec(nil1:nfl1),1,c0,x(nil:nfl),1)
		nil=(l-1)*norder+1
		nfl=(l)*norder
		nil1=l*norder+1
		nfl1=(l+1)*norder
		call my_zgemv('C',norder,norder,gt,h_int_up(nil:nfl,:),norder,vec(nil:nfl),1,c0,tmp(nil1:nfl1),1)
		endif
	enddo

	x=x+tmp
	deallocate(tmp)

	return
	end
!--------------------------------------------------------------------------------
    subroutine norm2(m_over,v,nrm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!	Norm of the wave-packet in the orthogonal basis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use datainput_module
	use complex_module
	implicit none
	integer l,i,j
	real*8 nrm,m_over(ntot1,norder)
	complex*16  v(ntot1),norm
	complex*16,allocatable :: x(:)
	complex*16,external :: scalar_conjg

	norm=c0
	allocate(x(ntot1))
	call hatomic_vec(m_over,v,x)
	norm=scalar_conjg(v,x)
	deallocate(x)
	
	nrm=real(norm)
	nrm=sqrt(nrm)

    return
    end 
!---------------------------------------------------------------
	subroutine mvec(job,m_chol,vec,x)
	use complex_module
	use datainput_module
	implicit none
	real*8 m_chol(ntot1,norder)
	complex*16 vec(ntot1),x(ntot1)
	integer l,j,jp
	character*1 job

	x=vec
	do l=0,lq
	call my_ztrmv('U',job,'N',norder,m_chol(l*norder+1:(l+1)*norder,:),norder,x(l*norder+1:(l+1)*norder),1)
	enddo

	return
	end
!---------------------------------------------------------
	subroutine mdecomp(m_over,m_chol)
!cholesky decomposition of overlap matrix (real symmetric)
	use datainput_module
	use complex_module
	use scratchdir_module
	use file_io_module
	implicit none
	integer info,l,j,jp
	real*8 m_over(ntot1,norder)
	real*8 m_chol(ntot1,norder)
	real*8 matrix1(norder,norder)
	character*2 cl

	m_chol=0d0
	do l=0,lq
	matrix1=0d0
	matrix1(:,:)=m_over(l*norder+1:(l+1)*norder,:)
	call dpotrf('U',norder,matrix1,norder,info)
		if(info.ne.0)then
		write(6,*) 'info:',info,'in dpotrf in mdecomp'
		stop
		endif
	m_chol(l*norder+1:(l+1)*norder,:)=matrix1(:,:)
	enddo

	return
	end
!-------------------------------------------------------------------------------
	subroutine msolve(job,m_chol,x,y)
	!solves U*y=x or U**T*y=x
	use datainput_module
	use complex_module
	implicit none
	real*8 m_chol(ntot1,norder)
	complex*16 y(ntot1),x(ntot1)
	integer j,jp,l
	character*1 job
	
	y=x
	do l=0,lq
	call my_ztrsv('U',job,'N',norder,m_chol(l*norder+1:(l+1)*norder,:),norder,y(l*norder+1:(l+1)*norder),1)
	enddo

	return
	end
!**********************************************************************************
	subroutine minverse(m_chol,x,y)
	use datainput_module
	use complex_module
	use scratchdir_module
	use file_io_module
	implicit none
	integer info,l,j,jp
	real*8 m_chol(ntot1,norder),xr(ntot1),xi(ntot1)
	complex*16 y(ntot1),x(ntot1)
	character*2 cl

	y=c0
	xr=real(x)
	xr=aimag(x)
	do l=0,lq
	call dpotrs('U',norder,1,m_chol(l*norder+1:(l+1)*norder,:),norder,xr(l*norder+1:(l+1)*norder),norder,info )
	call dpotrs('U',norder,1,m_chol(l*norder+1:(l+1)*norder,:),norder,xi(l*norder+1:(l+1)*norder),norder,info )
	y(l*norder+1:(l+1)*norder)=xr(l*norder+1:(l+1)*norder)+ci*xi(l*norder+1:(l+1)*norder)
	enddo

	return
	end
!--------------------------------------------------------------------------------------
	subroutine pulse_time(t,gt0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!	The time dependency of the interaction with the pulse is written as
!!!	H(int)=g(t)*mat_int*constants
!!!	with g(t)=	A(t) for velocity gauge
!!!			E(t)=-dA(t)/dt for length gauge
!!!	Then constants are added in the multiplication routines.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use datainput_module
	use pi_module
	use complex_module
	implicit none   
	real*8 t,t1,t2
	complex*16 gt0,cte

	gt0=c0
	if(gauge.eq.1)then
	cte=cmplx(amax)
	endif
	if(gauge.eq.2)then
	cte=-ci*amax
	endif

	if(env.eq.1)then
        gt0=cte*sin(omega*t)*sin(pi*t/tau)**2
	endif

	if(env.eq.2)then
		if(t.ge.tini1.and.t.le.tfin1)then
		gt0=cte*sin(omega*t)*sin(pi*t/tau)**2
		endif
		if(t.ge.tini2.and.t.le.tfin2)then
		gt0=cte*sin(omega1*(t-tfin1-tsep))*sin(pi*(t-tfin1-tsep)/tau1)**2
		endif
	endif

	if(env.eq.3)then
		if(t.ge.tini1.and.t.le.tau)then
		gt0=cte*sin(omega*t)*sin(pi*t/2d0/tau)**2
		endif
		if(t.gt.tau.and.t.le.tau+tau1)then
		gt0=cte*sin(omega*t)
		endif
		if(t.gt.tau+tau1.and.t.le.tau1+2*tau)then
		gt0=cte*sin(omega*t)*sin(pi*(t-tau1)/2d0/tau)**2
		endif
	endif

      	return
      	end



