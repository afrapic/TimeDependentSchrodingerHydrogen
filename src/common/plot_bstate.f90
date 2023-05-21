	subroutine plot_bstate(nn,l,yt)
	use datainput_module
	use complex_module
	use matrix_module
	use file_io_module
	use pi_module
	use scratchdir_module
	integer l,ir,i,ele,j,nn
	complex*16 yt(ntot1),wp
	character*2 cl,cn

	write(cn,'(i2.2)') nn
	write(cl,'(i2.2)') l
	open(unit=100,file=trim(dir_output)//'bstate-n'//cn//'l'//cl//'-'//trim(file_output)//'.dat')
	
	do ir=1,nquad
	wp=c0
		do i=1,norder
		wp=wp+yt(l*norder+i)*st(l,i,ir)
		enddo
	write(100,'(f20.10,2es30.15)') xq(ir),real(wp),aimag(wp)
	enddo
	close(100)	

	return
	end
