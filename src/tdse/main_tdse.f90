	PROGRAM TDSE
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
	character*2 cl,cn0
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
	write(6,*) '   TDSE PROPAGATION WITH PULSE STAGE    '
	write(6,*) '****************************************'
	write(6,*)

	!Reads input for problem
	call read_input

	write(6,*) 'Reading matrices'
	time1=dsecnd()
	allocate(h_atomic(ntot1,norder))
	allocate(m_over(ntot1,norder))
	allocate(m_chol(ntot1,norder))
	allocate(h_int_up(ntot2,norder))
	call read_matrix_raw(h_atomic,trim(dir_scratch)//'hatomic-'//trim(file_output)//'.dat')
	call read_matrix_raw(m_over,trim(dir_scratch)//'mover-'//trim(file_output)//'.dat')
	call read_matrix_raw(m_chol,trim(dir_scratch)//'mchol-'//trim(file_output)//'.dat')
	call read_matrix_raw(h_int_up,trim(dir_scratch)//'hintup-'//trim(file_output)//'.dat')
	time2=dsecnd()
	call print_wtime(time2-time1)
	write(6,*)

	write(6,*) 'Bound state reading'
	allocate(yt(ntot1))
	yt=c0
	call read_matrix_raw(yt,trim(dir_scratch)//'yt0-'//trim(file_output)//'.dat')
	write(6,*)
	
	write(6,*) 'Time propagation'
	!call propag_fatunla(h_atomic,m_over,m_chol,h_int_up,yt)
	call propag_arnoldi(h_atomic,m_over,m_chol,h_int_up,yt)
	write(6,*)

	deallocate(h_atomic,m_over,m_chol,h_int_up,yt)


	write(6,*) 'END STAGE'
	timef=dsecnd()
	call print_wtime(timef-timei)

	end
