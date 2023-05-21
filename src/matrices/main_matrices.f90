	PROGRAM MATRICES
	use datainput_module
	use file_io_module
	use scratchdir_module
	use matrix_module
	use complex_module
	implicit none
	integer ncom
	real*8,external :: dsecnd
	real*8 timei,timef,time1,time2
	integer l,ele
	character*2 cl
	complex*16,allocatable :: eigvalues(:),sturmq(:,:),dsturmq(:,:)
	real*8,allocatable :: h_atomic(:,:),h_int_low(:,:),h_int_up(:,:),m_over(:,:),m_chol(:,:)
	complex*16,allocatable :: yt(:)

	timei=dsecnd()

	!Check problem specific input file
	ncom = COMMAND_ARGUMENT_COUNT()
   	select case(ncom)
   	case(0)
		write(*,*) "Insuficient number of arguments, stop"
		stop
   	case(1)
       	call get_command_argument(1,file_input)
   	end select

	write(6,*) '****************************************'
	write(6,*) '   MATRICES FOR TDSE PROPAGATION        '
	write(6,*) '****************************************'
	write(6,*)

	!Reads input for problem
	call read_input

	!Reads quadrature weights and points
	allocate(xq(nquad),wq(nquad))
	call read_matrix_raw(xq,trim(dir_scratch)//'xq-'//trim(file_output)//'.dat')
	call read_matrix_raw(wq,trim(dir_scratch)//'wq-'//trim(file_output)//'.dat')
	if(ilong.eq.5)then
	allocate(xs0(nquad))
	call read_matrix_raw(xs0,trim(dir_scratch)//'xs0-'//trim(file_output)//'.dat')
	endif
	if(ilong.eq.6)then
	allocate(xs0(nquad))
	call read_matrix_raw(xs0,trim(dir_scratch)//'xs0-'//trim(file_output)//'.dat')
	allocate(xs1(nquad))
	call read_matrix_raw(xs1,trim(dir_scratch)//'xs1-'//trim(file_output)//'.dat')
	endif
	if(ilong.eq.7)then
	allocate(xs0(nquad))
	call read_matrix_raw(xs0,trim(dir_scratch)//'xs0-'//trim(file_output)//'.dat')
	allocate(xs1(nquad))
	call read_matrix_raw(xs1,trim(dir_scratch)//'xs1-'//trim(file_output)//'.dat')
	allocate(xs2(nquad))
	call read_matrix_raw(xs2,trim(dir_scratch)//'xs2-'//trim(file_output)//'.dat')
	endif

	!Reads the Sturmian basis in quadrature points for each angular momentum
	write(6,*)
      	write(6,*) 'Sturmian basis read'
	allocate(st(0:lq,norder,nquad),dst(0:lq,norder,nquad),beta(0:lq,norder))
	if(ilong.eq.5)then
	allocate(ss0(0:lq,norder,nquad))
	endif
	if(ilong.eq.6)then
	allocate(ss0(0:lq,norder,nquad))
	allocate(ss1(0:lq,norder,nquad))
	endif
	if(ilong.eq.7)then
	allocate(ss0(0:lq,norder,nquad))
	allocate(ss1(0:lq,norder,nquad))
	allocate(ss2(0:lq,norder,nquad))
	endif
	do ele=0,lq
      	time1=dsecnd()
	write(6,*) 'ele:',ele 
	write(cl,'(i2.2)') ele
	allocate(eigvalues(norder))
	allocate(sturmq(norder,nquad))
	allocate(dsturmq(norder,nquad))
	call read_matrix_raw(eigvalues,trim(dir_scratch)//'eigenvalues-l'//cl//'-'//trim(file_output)//'.dat')
	call read_matrix_raw(sturmq,trim(dir_scratch)//'st-l'//cl//'-'//trim(file_output)//'.dat')
	call read_matrix_raw(dsturmq,trim(dir_scratch)//'dst-l'//cl//'-'//trim(file_output)//'.dat')
	beta(ele,:)=eigvalues(:)
	st(ele,:,:)=sturmq(:,:)
	dst(ele,:,:)=dsturmq(:,:)
		if(ilong.eq.5)then
		call read_matrix_raw(sturmq,trim(dir_scratch)//'ss0-l'//cl//'-'//trim(file_output)//'.dat')
		ss0(ele,:,:)=sturmq(:,:)
		endif
		if(ilong.eq.6)then
		call read_matrix_raw(sturmq,trim(dir_scratch)//'ss0-l'//cl//'-'//trim(file_output)//'.dat')
		ss0(ele,:,:)=sturmq(:,:)
		call read_matrix_raw(sturmq,trim(dir_scratch)//'ss1-l'//cl//'-'//trim(file_output)//'.dat')
		ss1(ele,:,:)=sturmq(:,:)
		endif
		if(ilong.eq.7)then
		call read_matrix_raw(sturmq,trim(dir_scratch)//'ss0-l'//cl//'-'//trim(file_output)//'.dat')
		ss0(ele,:,:)=sturmq(:,:)
		call read_matrix_raw(sturmq,trim(dir_scratch)//'ss1-l'//cl//'-'//trim(file_output)//'.dat')
		ss1(ele,:,:)=sturmq(:,:)
		call read_matrix_raw(sturmq,trim(dir_scratch)//'ss2-l'//cl//'-'//trim(file_output)//'.dat')
		ss2(ele,:,:)=sturmq(:,:)
		endif
	if(allocated(eigvalues)) deallocate(eigvalues)
	if(allocated(sturmq)) deallocate(sturmq)
	if(allocated(dsturmq)) deallocate(dsturmq)
	time2=dsecnd()
      	call print_wtime(time2-time1)
	write(6,*)
      	enddo

	write(6,*) 'Evaluation of one dimensional integrals'
	time1=dsecnd()
	call evaluate_integrals_tdse
	time2=dsecnd()
	call print_wtime(time2-time1)
	deallocate(st,dst,xq,wq)
	if(allocated(ss0)) deallocate(ss0)
	if(allocated(xs0)) deallocate(xs0)
	if(allocated(ss1)) deallocate(ss1)
	if(allocated(xs1)) deallocate(xs1)
	if(allocated(ss2)) deallocate(ss2)
	if(allocated(xs2)) deallocate(xs2)
	write(6,*)

	write(6,*) 'Evaluation of matrices'
	time1=dsecnd()
	allocate(h_atomic(ntot1,norder))
	allocate(m_over(ntot1,norder))
	allocate(m_chol(ntot1,norder))
	allocate(h_int_low(ntot2,norder))
	allocate(h_int_up(ntot2,norder))
	call evaluate_matrix_tdse(h_atomic,h_int_low,h_int_up,m_over,m_chol)
	call write_matrix_raw(h_atomic,trim(dir_scratch)//'hatomic-'//trim(file_output)//'.dat')
	call write_matrix_raw(m_over,trim(dir_scratch)//'mover-'//trim(file_output)//'.dat')
	call write_matrix_raw(m_chol,trim(dir_scratch)//'mchol-'//trim(file_output)//'.dat')
	call write_matrix_raw(h_int_low,trim(dir_scratch)//'hintlow-'//trim(file_output)//'.dat')
	call write_matrix_raw(h_int_up,trim(dir_scratch)//'hintup-'//trim(file_output)//'.dat')
	time2=dsecnd()
	call print_wtime(time2-time1)
	write(6,*)

	deallocate(h_atomic,m_over,m_chol,h_int_low,h_int_up)
	if(allocated(beta)) deallocate(beta)
	if(allocated(st)) deallocate(st)
	if(allocated(dst)) deallocate(dst)

	write(6,*) 'END STAGE'
	timef=dsecnd()
	call print_wtime(timef-timei)

	end
