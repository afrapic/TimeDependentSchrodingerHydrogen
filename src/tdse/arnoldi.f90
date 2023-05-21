	subroutine propag_arnoldi(h_atomic,m_over,m_chol,h_int_up,yt)
	use datainput_module
	use file_io_module
	use scratchdir_module
	use complex_module
    	implicit none
	real*8 h_atomic(ntot1,norder)
	real*8 m_over(ntot1,norder)
	real*8 m_chol(ntot1,norder)
	real*8 h_int_up(ntot2,norder)
    	integer it,isave,l,j
	real*8 nrm,timei,timef,time1,tsave,dt_save
	real*8 t,t_new,h0,tfinal
	complex*16 yt(ntot1),gt
	real*8,allocatable :: erq(:)
	complex*16,allocatable :: y_save(:),rndm(:),vrq(:,:)
	character*2 cisave
	real*8,external :: dsecnd
	integer,parameter :: maxstep=2000000000

	open(unit=55,file=trim(dir_scratch)//'tsave-'//trim(file_output)//'.dat')
	!open(unit=85,file=trim(dir_output)//'arnoldi-'//trim(file_output)//'.dat')


	allocate(y_save(ntot1),rndm(ntot1))
	do j=1,ntot1
	call random_number(nrm)
	rndm(j)=cmplx(nrm)/tolrndm
	enddo
	allocate(vrq(ntot1,0:nkry),erq(0:nkry))
	vrq=c0
	erq=0d0

	!Initialize variables
	t=tini_pulse
	tfinal=tfin_pulse
	tfin=tfin_pulse
	call norm2(m_over,yt,nrm)
	yt=yt/nrm
	y_save=yt
	h0=h
	dt_save=(tfin-tini_pulse)/dble(num_save_vec)

	!prints initial vector
	call write_coef(0,tini_pulse,yt)

	write(6,*) 'Pulse propagation'
	isave=1
	tsave=tini_pulse+isave*dt_save

	timei=dsecnd()
	!Iteration to perform the propagation
	do it=1,maxstep
	y_save=yt

	!Arnoldi scheme
	call arnoldi(t,h_atomic,h_int_up,m_chol,rndm,yt,vrq,erq)
	t_new=t+h0
	y_save=yt
	
	call norm2(m_over,yt,nrm)
	write(6,*) t+h0,h0,nrm

		!Exit if tfin reached
		if(t_new.eq.tfin)then
		timef=dsecnd()
		call write_coef(-1,tfin_pulse,yt)
		call write_coef(num_save_vec,tfin_pulse,yt)
		write(cisave,'(i2.2)') num_save_vec
		!call write_matrix_raw(vrq,trim(dir_scratch)//'eigv_propag-'//trim(file_output)//'-'//cisave//'.dat')
		!write(85,'(i2,3es30.15)') isave,erq(0),erq(1),erq(2)
		exit
		endif

		!Saves vector if reached t_save
		if(isave.ne.num_save_vec)then
		if(t_new.ge.tsave)then
		call write_coef(isave,t_new,yt)
		!write(cisave,'(i2.2)') isave
		!call write_matrix_raw(vrq,trim(dir_scratch)//'eigv_propag-'//trim(file_output)//'-'//cisave//'.dat')
		!print*,vrq(:,0)
		!write(85,'(i2,3es30.15)') isave,erq(0),erq(1),erq(2)
		isave=isave+1
		tsave=tini_pulse+isave*dt_save
		endif
		endif
		
		!Adjust time step at the end of the pulse and at the end of free propagation
		if(t_new.gt.tfinal) then
		h0=tfinal-t
		t_new=t+h0
        	t=t_new
		yt=y_save
		else
         	t=t_new
         	endif

		!Starts free propagation at the end of the pulse in tfin_pulse
!		if(t_new.eq.tfin_pulse)then
!		write(6,*)
!		write(6,*) 'Free propagation'
!		tfinal=tfin
!		h0=h
!		endif

	enddo

	call print_wtime(timef-timei)
	deallocate(rndm,y_save)

    	return
    	end
!-------------------------------------------------------------------------------------------------
	subroutine arnoldi(t,h_atomic,h_int_up,m_chol,rndm,ct,vq,eq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!	Arnoldi scheme to solve
!!!	c'(t)=-i*inv(B)*H(t)*c(t) from t to t+h
!!!	by transforming first 
!!!	B=U**T*U and d=U*c so
!!!	d'(t)=-i*Ĥ(t)*d with Ĥ(t)=inv(U**T)*H(t)*inv(U)
!!!	Then you have 
!!!	d(t+h)=exp(-i*Ĥ(t)*h)*d(t)
!!!	and this is solve using a Krylov subspace of Ĥ(t)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use datainput_module
	use complex_module
	implicit none 
	real*8 t,tolq
	real*8 h_atomic(ntot1,norder),m_chol(ntot1,norder)
	real*8 h_int_up(ntot2,norder),eq(0:nkry)
	complex*16 ct(ntot1),rndm(ntot1),vq(ntot1,0:nkry)
	complex*16,allocatable :: q(:,:),v(:,:),hij(:,:),tmp(:,:)
	complex*16,allocatable :: qj(:),x(:),temp(:),temp1(:),yt(:),eig(:),vrq(:,:)
	complex*16 scalp,scalp1
	complex*16,external :: scalar_conjg
	integer i,j
	complex*16,allocatable :: vr(:,:)

	!Tolerance for the orthogonality
	tolq=1d-12

	!Initializes Krylov space and matrix coefficients
	allocate(q(ntot1,0:nkry+1),v(ntot1,0:nkry+1))
	q=c0
	v=c0

	!Transforms yt=U*ct
	allocate(yt(ntot1))
	call mvec('N',m_chol,ct,yt)

	!It adds a random vector to the initial state because otherwise is not linearly independent
	!if(t.lt.tini_pulse+0.1d0)then 
	q(:,0)=rndm(:)
	!endif
	q(:,0)=q(:,0)+yt(:)
	scalp=scalar_conjg(q(:,0),q(:,0))
	q(:,0)=q(:,0)/sqrt(scalp)

	!hij is a Hessenber matrix.
	allocate(hij(0:nkry+1,0:nkry+1))
	hij=c0

	!Arnoldi-Lanczos algorithm
	do j=0,nkry
	!vec(j)=Ĥ*q(j)
	allocate(qj(ntot1),x(ntot1))
	qj(:)=q(:,j)
	call hamilt_vec(t,h_atomic,h_int_up,m_chol,qj,x)
	v(:,j)=x(:)
		do i=0,j
		scalp=scalar_conjg(q(:,i),v(:,j))
		hij(i,j)=scalp
		v(:,j)=v(:,j)-hij(i,j)*q(:,i)
		enddo
		scalp=scalar_conjg(v(:,j),v(:,j))
		hij(j+1,j)=sqrt(scalp)
		if(abs(hij(j+1,j)).lt.1d-8)then
		print*,hij(j+1,j)
		print*,'Arnoldi break'
		stop
		endif
	deallocate(qj,x)
	q(:,j+1)=v(:,j)/hij(j+1,j)
	scalp=scalar_conjg(q(:,j+1),q(:,j+1))
	q(:,j+1)=q(:,j+1)/sqrt(scalp)
			!Checks orthogonality after renormalizing
			do i=0,j
			scalp=scalar_conjg(q(:,i),q(:,j+1))
				if(abs(scalp).gt.tolq)then
				scalp1=scalar_conjg(q(:,i),q(:,i))
				q(:,j+1)=q(:,j+1)-scalp*q(:,i)/scalp1
				endif
			scalp1=scalar_conjg(q(:,j+1),q(:,j+1))
			q(:,j+1)=q(:,j+1)/sqrt(scalp1)
			enddo
	enddo
	deallocate(v)

	!Diagonalization of hij that is upper Hessenberg
	allocate(vr(0:nkry,0:nkry),eig(0:nkry))
	call diagonalizacion_krylov(hij,eig,vr)
	eq=real(eig)
	deallocate(hij)
	!write(85,'(6es15.5)') t,real(eig(0:4))

	!qj=vr*q**('T)*dv=d"(t)
	allocate(tmp(ntot1,0:nkry),temp(0:nkry))
	tmp(:,0:nkry)=q(:,0:nkry)
	deallocate(q)
	allocate(vrq(ntot1,0:nkry))
	call transf_vec_kry(tmp,vr,yt,temp,vrq)
	
	!diag(exp(-i*h*eig(0)),....,exp(-i*h*eig(n)))*qj=x=d"(t+h)
	allocate(temp1(0:nkry))
	temp1=c0
	do j=0,nkry
	temp1(j)=exp(-ci*h*eig(j))*temp(j)
	enddo
	deallocate(temp)

	!dv=q*inv(vr)*x=d(t+h)
	call backtransf_vec_kry(tmp,vr,temp1,yt)

	!ct=inv(U)*yt
	call msolve('N',m_chol,yt,ct)

	do j=0,nkry
	vq(:,j)=c0
	call msolve('N',m_chol,vrq(:,j),vq(:,j))
	enddo
	
	if(allocated(q)) deallocate(q)
	if(allocated(v)) deallocate(v)
	if(allocated(hij)) deallocate(hij)
	if(allocated(qj)) deallocate(qj)
	if(allocated(x)) deallocate(x)	
	if(allocated(vr)) deallocate(vr)
	if(allocated(temp)) deallocate(temp)
	if(allocated(temp1)) deallocate(temp1)
	if(allocated(eig)) deallocate(eig)
	if(allocated(tmp)) deallocate(tmp)
	if(allocated(yt)) deallocate(yt)
	if(allocated(vrq)) deallocate(vrq)

	return
	end

!-----------------------------------------------------------
	subroutine diagonalizacion_krylov(mh,eigval,eigvec)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!	Diagonalization of the upper Hessenber matrix, but in this case,
!!!	since the hamiltonian is hermitian, the matrix should be tridiagonal and
!!!	symmetric, so a general subroutine for hermitian matrix is used to
!!!	obtain the eigenvectors and eigenvalues.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use datainput_module
	use complex_module
	implicit none
	integer info,j
	complex*16 mh(0:nkry+1,0:nkry+1),eigval(0:nkry),eigvec(0:nkry,0:nkry)
	real*8,allocatable :: rwork(:),tmp1(:) 
	complex*16,allocatable :: tmp(:,:),work(:)

	allocate(tmp(0:nkry,0:nkry),tmp1(0:nkry),work(2*(nkry+1)-1),rwork(3*(nkry+1)-2))
	tmp(0:nkry,0:nkry)=mh(0:nkry,0:nkry)
	call ZHEEV('v','u',nkry+1,tmp,nkry+1,tmp1,work,2*(nkry+1)-1,rwork,info)
	deallocate(work,rwork)
		if(info.ne.0)then
		print*,'info arnoldi ZHEEV:',info
		stop
		endif	
	eigval=tmp1+ci*0d0
	eigvec=transpose(tmp)
	deallocate(tmp,tmp1)

	return
	end
!-------------------------------------------------------------------
	subroutine transf_vec_kry(mq,eigvec,vec,x,eigvecq)
	!Transforms in the Krylov space the vector vec=yt(t) in yt"(t)=eigvec*mq'*yt(t). Output in x
	use datainput_module
	use complex_module
	integer i
	complex*16 mq(ntot1,0:nkry),eigvec(0:nkry,0:nkry),vec(ntot1),x(0:nkry),eigvecq(ntot1,0:nkry)
	complex*16,allocatable :: tmp(:)
	complex*16,external :: scalar_conjg

	!mq'*vec=tmp
	allocate(tmp(0:nkry))
	do i=0,nkry
	tmp(i)=scalar_conjg(mq(:,i),vec(:))
	enddo

	!eigvec*tmp=x 
	call ZGEMV('n',nkry+1,nkry+1,c1,eigvec,nkry+1,tmp,1,c0,x,1)
	deallocate(tmp)	

	!mq*eigvec=eigvecq 
	do i=0,nkry
	call ZGEMV('n',ntot1,nkry+1,c1,mq(:,:),ntot1,eigvec(:,i),1,c0,eigvecq(:,i),1)
	enddo
	
	return
	end subroutine
!-------------------------------------------------------------------
	subroutine backtransf_vec_kry(mq,eigvec,vec,x)
	!Transforms to regular space the vector vec=yt"(t+h) in yt(t+h)=mq*inv(eigvec)*yt(t+h). Output in x
	use datainput_module
	use complex_module
      	integer info,i
	integer,allocatable :: ipiv(:)
	complex*16 mq(ntot1,0:nkry),eigvec(0:nkry,0:nkry),vec(0:nkry),x(ntot1)

	!inv(eigvec)*vec=qj solving eigvec*qj=vec, stores in vec
	allocate(ipiv(0:nkry))
	call ZGESV(nkry+1,1,eigvec,nkry+1,ipiv,vec,nkry+1,info)
	deallocate(ipiv)
	if(info.ne.0)then
	print*,'info arnoldi ZGESV:',info
	stop
	endif

	!mq*vec=x
	do i=0,lq
	call ZGEMV('N',norder,nkry+1,c1,mq(i*norder+1:(i+1)*norder,:),norder,vec,1,c0,x(i*norder+1:(i+1)*norder),1)
	enddo

	return
	end subroutine
!---------------------------------------------------------------------------------------
	subroutine hamilt_vec(t,h_atomic,h_int_up,m_chol,vec,x0)
	!Calculates H*vec=x. mh is the full hamiltonian
	!ouput in x
	use complex_module
	use datainput_module
	implicit none
	real*8 t
	complex*16 vec(ntot1),x0(ntot1)
	real*8 h_atomic(ntot1,norder)
	real*8 m_chol(ntot1,norder)
	real*8 h_int_up(ntot2,norder)
	complex*16 gt
	integer l,j
	complex*16,allocatable :: x(:),xi(:)

	!Hamiltonian for t
	call pulse_time(t,gt)

	!x0=inv(U)*vec
	call msolve('N',m_chol,vec,x0)

	!Atomic hamiltonian
	allocate(x(ntot1))
	call hatomic_vec(h_atomic,x0,x) 

	!Dipole interaction
	allocate(xi(ntot1))
	xi=c0
	if(t.le.tfin_pulse)then
		call hint_vec(gt,h_int_up,x0,xi)
	endif
	x=x+xi
	deallocate(xi)

	!x0=inv(U**T)*x
      x0=c0
	call msolve('T',m_chol,x,x0)
	deallocate(x)

	return
	end subroutine
!--------------------------------------------------------------------------------------------------
	subroutine write_coef(isave,tsave,yt)
	use datainput_module
	use file_io_module
	use complex_module
	use scratchdir_module
	implicit none
	integer isave
	real*8 tsave
	complex*16 yt(ntot1)
	character*2 cisave

	if(isave.ge.0)then
	write(cisave,'(i2.2)') isave
	call write_matrix_raw(yt,trim(dir_scratch)//'yt-'//trim(file_output)//'-'//cisave//'.dat')
	write(55,*) isave,tsave
	else
	call write_matrix_raw(yt,trim(dir_scratch)//'yt-'//trim(file_output)//'-tfin_pulse.dat')
	endif
	

	return
	end


