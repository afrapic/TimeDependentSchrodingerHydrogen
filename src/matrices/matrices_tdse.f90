	subroutine evaluate_matrix_tdse(h_atomic,h_int_low,h_int_up,m_over,m_chol)
	use complex_module
	use datainput_module
	use file_io_module
	use scratchdir_module
	use matrix_module
	implicit none
	real*8 h_atomic(ntot1,norder)
	real*8 h_int_low(ntot2,norder)
	real*8 h_int_up(ntot2,norder) 
	real*8 m_over(ntot1,norder)
	real*8 m_chol(ntot1,norder)
	integer l,j,jp
	real*8,allocatable :: matrix1(:,:),matrix2(:,:),matrix3(:,:),matrix4(:,:),matrix5(:,:)
	character*2 cl,clp

	h_atomic=0d0
	h_int_low=0d0
	h_int_up=0d0
	m_over=0d0
	m_chol=0d0

	!Atomic hamiltonian, interaction and overlap matrices
	do l=0,lq
	allocate(matrix1(norder,norder),matrix2(norder,norder))
	write(cl,'(i2.2)') l
	call read_matrix_raw(matrix1,trim(dir_scratch)//'integ_pot-l'//cl//'-'//trim(file_output)//'.dat')
	call read_matrix_raw(matrix2,trim(dir_scratch)//'integ_over-l'//cl//'-'//trim(file_output)//'.dat')
	if(ilong.eq.5)then
	allocate(matrix3(norder,norder))
	call read_matrix_raw(matrix3,trim(dir_scratch)//'integ_vw0-l'//cl//'-'//trim(file_output)//'.dat')
	endif
	if(ilong.eq.6)then
	allocate(matrix3(norder,norder))
	call read_matrix_raw(matrix3,trim(dir_scratch)//'integ_vw0-l'//cl//'-'//trim(file_output)//'.dat')
	allocate(matrix4(norder,norder))
	call read_matrix_raw(matrix4,trim(dir_scratch)//'integ_vw1-l'//cl//'-'//trim(file_output)//'.dat')
	endif
	if(ilong.eq.7)then
	allocate(matrix3(norder,norder))
	call read_matrix_raw(matrix3,trim(dir_scratch)//'integ_vw0-l'//cl//'-'//trim(file_output)//'.dat')
	allocate(matrix4(norder,norder))
	call read_matrix_raw(matrix4,trim(dir_scratch)//'integ_vw1-l'//cl//'-'//trim(file_output)//'.dat')
	allocate(matrix5(norder,norder))
	call read_matrix_raw(matrix5,trim(dir_scratch)//'integ_vw2-l'//cl//'-'//trim(file_output)//'.dat')
	endif
		do j=1,norder
			do jp=1,norder
			h_atomic(l*norder+j,jp)=real(matrix2(j,jp)*energy-beta(l,jp)*matrix1(j,jp))
			if(ilong.eq.5)then
			h_atomic(l*norder+j,jp)=h_atomic(l*norder+j,jp)+real(matrix3(j,jp))
			endif
			if(ilong.eq.6)then
			h_atomic(l*norder+j,jp)=h_atomic(l*norder+j,jp)+real(matrix3(j,jp))+real(matrix4(j,jp))
			endif
			if(ilong.eq.7)then
			h_atomic(l*norder+j,jp)=h_atomic(l*norder+j,jp)+real(matrix3(j,jp))+real(matrix4(j,jp))+real(matrix5(j,jp))
			endif
			m_over(l*norder+j,jp)=real(matrix2(j,jp))
			enddo
		enddo
	if(allocated(matrix1)) deallocate(matrix1)
	if(allocated(matrix2)) deallocate(matrix2)
	if(allocated(matrix3)) deallocate(matrix3)
	if(allocated(matrix4)) deallocate(matrix4)
	if(allocated(matrix5)) deallocate(matrix5)
	if(l.lt.lq)then
	write(clp,'(i2.2)') l+1
		if(gauge.eq.1)then
		allocate(matrix1(norder,norder))
		call read_matrix_raw(matrix1,trim(dir_scratch)//'integ_r-l'//cl//'-lp'//clp//'-'//trim(file_output)//'.dat')
		do j=1,norder
			do jp=1,norder
			h_int_low(l*norder+j,jp)=real(matrix1(j,jp))*sqrt((dble(l+1d0)**2)/(4d0*dble(l+1d0)**2-1d0))
			enddo
		enddo
		if(allocated(matrix1)) deallocate(matrix1)
		endif
		if(gauge.eq.2)then
		allocate(matrix1(norder,norder),matrix2(norder,norder))
		call read_matrix_raw(matrix1,trim(dir_scratch)//'integ_cou-l'//cl//'-lp'//clp//'-'//trim(file_output)//'.dat')
		call read_matrix_raw(matrix2,trim(dir_scratch)//'integ_dr-l'//cl//'-lp'//clp//'-'//trim(file_output)//'.dat')
		do j=1,norder
			do jp=1,norder
			h_int_low(l*norder+j,jp)=real(matrix2(j,jp)-(l+1d0)*matrix1(j,jp))*&
									sqrt((dble(l+1d0)**2)/(4d0*dble(l+1d0)**2-1d0))
			enddo
		enddo
		if(allocated(matrix1)) deallocate(matrix1)
		if(allocated(matrix2)) deallocate(matrix2)
		endif
	endif
	if(l.gt.0)then
	write(clp,'(i2.2)') l-1
		if(gauge.eq.1)then
		allocate(matrix1(norder,norder))
		call read_matrix_raw(matrix1,trim(dir_scratch)//'integ_r-l'//cl//'-lp'//clp//'-'//trim(file_output)//'.dat')
		do j=1,norder
			do jp=1,norder
			h_int_up((l-1)*norder+j,jp)=real(matrix1(j,jp))*sqrt((dble(l)**2)/(4d0*dble(l)**2-1d0))
			enddo
		enddo
		if(allocated(matrix1)) deallocate(matrix1)
		endif
		if(gauge.eq.2)then
		allocate(matrix1(norder,norder),matrix2(norder,norder))
		call read_matrix_raw(matrix1,trim(dir_scratch)//'integ_cou-l'//cl//'-lp'//clp//'-'//trim(file_output)//'.dat')
		call read_matrix_raw(matrix2,trim(dir_scratch)//'integ_dr-l'//cl//'-lp'//clp//'-'//trim(file_output)//'.dat')
		do j=1,norder
			do jp=1,norder
			h_int_up((l-1)*norder+j,jp)=real(matrix2(j,jp)+(l*1d0)*matrix1(j,jp))*sqrt((dble(l)**2)/(4d0*dble(l)**2-1d0))
			enddo
		enddo
		if(allocated(matrix1)) deallocate(matrix1)
		if(allocated(matrix2)) deallocate(matrix2)
		endif
	endif
	enddo

	!Cholesky decomposition of overlap overlap=L*L**T
	call mdecomp(m_over,m_chol)

	return
	end
!---------------------------------------------------------------------------------------------------------
	subroutine evaluate_integrals_tdse
	use matrix_module
	use datainput_module
	use file_io_module
	use scratchdir_module
	integer ele,j,jp
	complex*16,allocatable :: stl(:,:),dstl(:,:)
	complex*16,allocatable :: stlp(:,:)
	complex*16,allocatable :: matrix(:,:)
	character*2 cl,clp

	allocate(matrix(norder,norder),stl(norder,nquad))
	do ele=0,lq
	write(cl,'(i2.2)') ele
	stl(:,:)=st(ele,:,:)
	matrix=c0
	call integ_st_pot(stl,matrix)
	call write_matrix_raw(real(matrix),trim(dir_scratch)//'integ_pot-l'//cl//'-'//trim(file_output)//'.dat')
	matrix=c0
	call integ_st_over(stl,matrix)
	call write_matrix_raw(real(matrix),trim(dir_scratch)//'integ_over-l'//cl//'-'//trim(file_output)//'.dat')
	if(ilong.eq.5)then
	stl(:,:)=ss0(ele,:,:)
	call integ_st_vw0(stl,matrix)
	call write_matrix_raw(real(matrix),trim(dir_scratch)//'integ_vw0-l'//cl//'-'//trim(file_output)//'.dat')
	endif
	if(ilong.eq.6)then
	stl(:,:)=ss0(ele,:,:)
	call integ_st_vw0(stl,matrix)
	call write_matrix_raw(real(matrix),trim(dir_scratch)//'integ_vw0-l'//cl//'-'//trim(file_output)//'.dat')
	stl(:,:)=ss1(ele,:,:)
	call integ_st_vw1(stl,matrix)
	call write_matrix_raw(real(matrix),trim(dir_scratch)//'integ_vw1-l'//cl//'-'//trim(file_output)//'.dat')
	endif
	if(ilong.eq.7)then
	stl(:,:)=ss0(ele,:,:)
	call integ_st_vw0(stl,matrix)
	call write_matrix_raw(real(matrix),trim(dir_scratch)//'integ_vw0-l'//cl//'-'//trim(file_output)//'.dat')
	stl(:,:)=ss1(ele,:,:)
	call integ_st_vw1(stl,matrix)
	call write_matrix_raw(real(matrix),trim(dir_scratch)//'integ_vw1-l'//cl//'-'//trim(file_output)//'.dat')
	stl(:,:)=ss2(ele,:,:)
	call integ_st_vw2(stl,matrix)
	call write_matrix_raw(real(matrix),trim(dir_scratch)//'integ_vw2-l'//cl//'-'//trim(file_output)//'.dat')
	endif
	stl(:,:)=st(ele,:,:)
		if(ele.lt.lq)then
		write(clp,'(i2.2)') ele+1
		allocate(stlp(norder,nquad),dstl(norder,nquad))
		stlp(:,:)=st(ele+1,:,:)
		matrix=c0
		call integ_st_rq(1,stl,stlp,matrix)
		call write_matrix_raw(real(matrix),trim(dir_scratch)//'integ_r-l'//cl//'-lp'//clp//'-'//trim(file_output)//'.dat')
		matrix=c0
		call integ_st_rq(-1,stl,stlp,matrix)
		call write_matrix_raw(real(matrix),trim(dir_scratch)//'integ_cou-l'//cl//'-lp'//clp//'-'//trim(file_output)//'.dat')
		dstl(:,:)=dst(ele,:,:)
		matrix=c0
		call integ_st_dr(dstl,stlp,matrix)
		call write_matrix_raw(real(matrix),trim(dir_scratch)//'integ_dr-l'//cl//'-lp'//clp//'-'//trim(file_output)//'.dat')
		deallocate(stlp,dstl)
		endif
		if(ele.gt.0)then
		write(clp,'(i2.2)') ele-1
		allocate(stlp(norder,nquad),dstl(norder,nquad))
		stlp(:,:)=st(ele-1,:,:)
		matrix=c0
		call integ_st_rq(1,stl,stlp,matrix)
		call write_matrix_raw(real(matrix),trim(dir_scratch)//'integ_r-l'//cl//'-lp'//clp//'-'//trim(file_output)//'.dat')
		matrix=c0
		call integ_st_rq(-1,stl,stlp,matrix)
		call write_matrix_raw(real(matrix),trim(dir_scratch)//'integ_cou-l'//cl//'-lp'//clp//'-'//trim(file_output)//'.dat')
		dstl(:,:)=dst(ele,:,:)
		matrix=c0
		call integ_st_dr(dstl,stlp,matrix)
		call write_matrix_raw(real(matrix),trim(dir_scratch)//'integ_dr-l'//cl//'-lp'//clp//'-'//trim(file_output)//'.dat')
		deallocate(stlp,dstl)
		endif
	enddo
	deallocate(matrix,stl)

	return
	end
