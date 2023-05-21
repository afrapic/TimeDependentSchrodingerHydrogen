!---------------------------------------------------------------------------------
	subroutine integrals_spect(ngridsk1,integ_cou,xq2,wq2)
	use matrix_module
	use datainput_module
	use complex_module
	use scratchdir_module
	use file_io_module
	integer l,ik
	complex*16,allocatable :: stl(:,:),vector(:)
	complex*16 integ_cou(ngridsk1,0:lq,norder)
	character*2 cl,clp
	real*8 kmax,dk,en,den,xq1(nquad),wq1(nquad)


	den=emax/dble(ngridsk1+1)

	do l=0,lq
	write(cl,'(i2.2)') l
	allocate(stl(norder,nquad))
	call read_matrix_raw(stl,trim(dir_scratch)//'st-l'//cl//'-'//trim(file_output)//'.dat')
		do ik=1,ngridsk1
		allocate(vector(norder))
		vector=c0
		en=ik*den 
		if(en.gt.0.0001d0)then
		call integ_st_cou(en,l,stl,vector,xq2,wq2)
		endif
		integ_cou(ik,l,:)=vector(:)
		deallocate(vector)
		enddo
	deallocate(stl)
	enddo

	return
	end
!-----------------------------------------------------------------------------------------
	subroutine integ_st_cou(en1,l1,stl,integ,xq2,wq2)
	use datainput_module
	use complex_module
	use pi_module
	use scratchdir_module
	use file_io_module
	implicit none
	integer i,j,l1
	complex*16 integ(norder)
	real*8 a0,en1,k,f,fp,g,gp,err,xq2(nquad),wq2(nquad)
	complex*16 stl(norder,nquad),gam,sigma
	real*8,allocatable :: fc(:)


 	allocate(fc(nquad))
	fc=0d0
	if(ilong.eq.0)then
      do i=1,nquad
	call scoul(-charge,en1,l1,xq2(i),f,fp,g,gp,err) 
	fc(i)=f
	enddo
	endif
	if(ilong.eq.5)then
	call svw(en1,l1,fc,xq2)
	endif
	if(ilong.eq.6)then
	call svw(en1,l1,fc,xq2)
	endif

	a0=rmax/2d0
	k=sqrt(2d0*en1)
	call gammaln(l1+1d0-ci/k,gam)
	sigma=atan2(aimag(exp(gam)),real(exp(gam)))
	do j=1,norder
	integ(j)=c0
		do i=1,nquad
		integ(j)=integ(j)+wq2(i)*stl(j,i)*fc(i)*a0*exp(ci*sigma)
		enddo
	enddo

	if(allocated(fc)) deallocate(fc)

	return
	end
