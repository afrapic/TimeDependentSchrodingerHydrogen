	PROGRAM EIGVAL
	use datainput_module
	use file_io_module
	use scratchdir_module
	use matrix_module
	use complex_module
	implicit none
	integer ncom
	real*8,external :: dsecnd
	real*8 timei,timef,time1,time2,tsave
	integer l,ele,i,nn
	character*2 cl
	complex*16,allocatable :: sturmq(:,:)
	real*8,allocatable :: h_atomic(:,:),m_over(:,:),m_chol(:,:),h_int_low(:,:),h_int_up(:,:)

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
	write(6,*) '   CALCULATION OF BOUND STATES          '
	write(6,*) '****************************************'
	write(6,*)

	!Reads input for problem
	call read_input

	!Reads quadrature weights and points
	allocate(xq(nquad),wq(nquad))
	call read_matrix_raw(xq,trim(dir_scratch)//'xq-'//trim(file_output)//'.dat')
	call read_matrix_raw(wq,trim(dir_scratch)//'wq-'//trim(file_output)//'.dat')
	!Reads the Sturmian basis in quadrature points for each angular momentum
	write(6,*)
	write(6,*) 'Sturmian basis read'
	allocate(st(0:lq,norder,nquad))
	do ele=0,lq
	time1=dsecnd()
	write(6,*) 'ele:',ele 
	write(cl,'(i2.2)') ele
	allocate(sturmq(norder,nquad))
	call read_matrix_raw(sturmq,trim(dir_scratch)//'st-l'//cl//'-'//trim(file_output)//'.dat')
	st(ele,:,:)=sturmq(:,:)
	if(allocated(sturmq)) deallocate(sturmq)
	time2=dsecnd()
	call print_wtime(time2-time1)
	write(6,*)
	enddo

	write(6,*) 'Reading matrices'
	time1=dsecnd()
	allocate(h_atomic(ntot1,norder))
	allocate(m_over(ntot1,norder))
	allocate(h_int_low(ntot2,norder))
	allocate(h_int_up(ntot2,norder)) 
	call read_matrix_raw(h_atomic,trim(dir_scratch)//'hatomic-'//trim(file_output)//'.dat')
	call read_matrix_raw(m_over,trim(dir_scratch)//'mover-'//trim(file_output)//'.dat')
	call read_matrix_raw(h_int_low,trim(dir_scratch)//'hintlow-'//trim(file_output)//'.dat')
	call read_matrix_raw(h_int_up,trim(dir_scratch)//'hintup-'//trim(file_output)//'.dat')
	
	time2=dsecnd()
	call print_wtime(time2-time1)
	write(6,*)

	write(6,*) 'Eigenvalues calculation'
	open(unit=55,file=trim(dir_scratch)//'tsave-'//trim(file_output)//'.dat')
	do i=1,num_save_vec+1
	read(55,*) nn,tsave
	call eigenstate_t(h_atomic,m_over,h_int_up,h_int_low,tsave,i)
	enddo

	deallocate(h_atomic,m_over,h_int_up,h_int_low)
	if(allocated(xq)) deallocate(xq)
	if(allocated(st)) deallocate(st)
	if(allocated(wq)) deallocate(wq)

	write(6,*) 'END STAGE'
	timef=dsecnd()
	call print_wtime(timef-timei)

	end
