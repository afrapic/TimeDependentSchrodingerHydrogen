	PROGRAM SPECTRUM
	use datainput_module
	use file_io_module
	use scratchdir_module
	use matrix_module
	use complex_module
	use pi_module
	implicit none
	integer ncom,ik,ir,iq,nk
	real*8,external :: dsecnd
	real*8 timei,timef,time1,time2,kmax,dk,spct,ang
	real*8,allocatable :: spctl(:),xq2(:),wq2(:)
	complex*16 en
	complex*16,allocatable :: integ_cou1(:,:,:),integ_cou(:,:,:)
	integer j
	character*6 cr

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
	write(6,*) '         SPECTRUM PLOT STAGE            '
	write(6,*) '****************************************'
	write(6,*)

	!Reads input for problem
	call read_input

	call plot_wp
	!call plot_flux
	!
	if(plot_pop.eq.1)then
	write(6,*) 'Population'
	!call population_cont_p
	call population
	endif
	if(allocated(dst)) deallocate(dst)

	if(plot_spect_td.eq.1)then
	write(6,*) 'Energy spectrum from td'
	nk=ngridsk
	allocate(integ_cou1(nk,0:lq,norder))
	print*,'integral'
	allocate(xq2(nquad),wq2(nquad))
	call read_matrix_raw(xq2,trim(dir_scratch)//'xq-'//trim(file_output)//'.dat')
	call read_matrix_raw(wq2,trim(dir_scratch)//'wq-'//trim(file_output)//'.dat')
	call integrals_spect(nk,integ_cou1,xq2,wq2)
	call spect1(integ_cou1)
	endif

	if(plot_rm.eq.1)then
	write(6,*) 'Mean radius'
	call plot_rmean
	endif


	if(allocated(xq)) deallocate(xq)
	if(allocated(wq)) deallocate(wq)
	if(allocated(st)) deallocate(st)
	if(allocated(dst)) deallocate(dst)

	write(6,*)
	write(6,*) 'END STAGE'
	timef=dsecnd()
	call print_wtime(timef-timei)

	end
