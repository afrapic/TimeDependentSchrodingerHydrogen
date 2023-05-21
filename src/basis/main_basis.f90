	PROGRAM BUILDBASIS
!####################################################################################################
!	This stage computes the Sturmian basis set using Bsplines to expand the radial part of the
!	sturmian function. 
!	The output is in the /scratch folder, it prints the eigenvalues, the quadrature weights and
!	points and the sturmian functions evaluated in the quadrature for each electron and each
!	angular momentum asked.
!	subroutine by Dr. A L Frapiccini
!	report bugs to afrapic@uns.edu.ar
!####################################################################################################
	use datainput_module
	use file_io_module
	use scratchdir_module
	use matrix_module
	implicit none
	integer ncom
	real*8 timei,timef,time1,time2
	integer l
	real*8,external :: dsecnd

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
	write(6,*) '     STURMIAN BASIS SET CALCULATION     '
	write(6,*) '****************************************'
	write(6,*)

	!Reads input for problem
	call read_input

	write(6,*) '-----Bound basis for TDSE-----'
	time1=dsecnd()
	!Bsplines construction
	write(6,*) 'Bspline construction'
	l=nb-kp+1
	allocate(tx(l+2*kp-1))
	call bspline_data(1,l,tx,kp,rmax,nquad_bsp,nb)

	!Evaluates and prints quadrature weights and points for printing sturmians
	allocate(xq(nquad),wq(nquad),xs0(nquad),xs1(nquad),xs2(nquad))
	call evaluate_quad_sturmians(nquad,rmax,xq,wq,rc0,delta0,xs0,rc1,delta1,xs1,rc2,delta2,xs2)
	call write_matrix_raw(xq,trim(dir_scratch)//'xq-'//trim(file_output)//'.dat')
	call write_matrix_raw(wq,trim(dir_scratch)//'wq-'//trim(file_output)//'.dat')
	if(ilong.eq.5)then
	call write_matrix_raw(xs0,trim(dir_scratch)//'xs0-'//trim(file_output)//'.dat')
	endif
	if(ilong.eq.6)then
	call write_matrix_raw(xs0,trim(dir_scratch)//'xs0-'//trim(file_output)//'.dat')
	call write_matrix_raw(xs1,trim(dir_scratch)//'xs1-'//trim(file_output)//'.dat')
	endif
	if(ilong.eq.7)then
	call write_matrix_raw(xs0,trim(dir_scratch)//'xs0-'//trim(file_output)//'.dat')
	call write_matrix_raw(xs1,trim(dir_scratch)//'xs1-'//trim(file_output)//'.dat')
	call write_matrix_raw(xs2,trim(dir_scratch)//'xs2-'//trim(file_output)//'.dat')
	endif

	!Build and prints the Sturmian basis in quadrature points for each angular momentum
      write(6,*) 'Sturmian basis build'
	call build_basis
	deallocate(tx,xq,wq,wd,xd,bb,dbb,xs0,xs1,xs2)
	time2=dsecnd()
	call print_wtime(time2-time1)
	write(6,*)
	
	write(6,*) 'END STAGE'
	timef=dsecnd()
	call print_wtime(timef-timei)

	end
