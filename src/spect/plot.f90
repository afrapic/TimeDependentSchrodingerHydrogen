	subroutine population_cont_p
	use datainput_module
	use complex_module
        use scratchdir_module
        use file_io_module
	use matrix_module
	implicit none
	real*8,allocatable :: h_atomic(:,:),m_over(:,:)
	complex*16 cnp(200)
	complex*16,allocatable :: yt(:),ytemp(:),ynp(:),integ_cou(:,:),vector(:),stl(:,:)
	real*8,allocatable :: ma(:,:),mb(:,:),eigenvalues(:)
	integer info,lwork
	real*8 tsave,en,den,cwork(1)
	real*8,allocatable :: work(:),xq2(:),wq2(:)
	integer nn,ll,isave,i,j,ik
        character*2 cisave,cl
	complex*16,external :: scalar_conjg

	allocate(h_atomic(ntot1,norder))
	allocate(m_over(ntot1,norder))
	call read_matrix_raw(h_atomic,trim(dir_scratch)//'hatomic-'//trim(file_output)//'.dat')
	call read_matrix_raw(m_over,trim(dir_scratch)//'mover-'//trim(file_output)//'.dat')

	ll=1
	allocate(ma(norder,norder),mb(norder,norder))
	ma=0d0
	mb=0d0
	ma(:,:)=h_atomic(ll*norder+1:(ll+1)*norder,:)
	mb(:,:)=m_over(ll*norder+1:(ll+1)*norder,:)

	lwork=-1
	allocate(eigenvalues(norder))
	call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,cwork,lwork,info)

	lwork=int(cwork(1))
	allocate(work(lwork))
	call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,work,lwork,info)
		if(info.ne.0)then
		write(6,*) 'info:',info,'in dsygv'
		endif
       deallocate(work,h_atomic,mb)
	
!	open(unit=80,file=trim(dir_output)//'population-eigval-cont-'//trim(file_output)//'.dat')
!	ll=1
!	write(cl,'(i2.2)') ll
!	allocate(xq2(nquad),wq2(nquad))
!	call read_matrix_raw(xq2,trim(dir_scratch)//'xq-'//trim(file_output)//'.dat')
!	call read_matrix_raw(wq2,trim(dir_scratch)//'wq-'//trim(file_output)//'.dat')
!	allocate(stl(norder,nquad))
!	call read_matrix_raw(stl,trim(dir_scratch)//'st-l'//cl//'-'//trim(file_output)//'.dat')
!	den=emax/30000d0
!	allocate(integ_cou(100,norder))
!	integ_cou=c0
!		do ik=1,100
!		allocate(vector(norder))
!		vector=c0
!		en=ik*den 
!		write(80,'(f20.10)') ik*den
!		call integ_st_cou(en,ll,stl,vector,xq2,wq2)
!		integ_cou(ik,:)=vector(:)
!		deallocate(vector)
!		enddo
!	deallocate(stl,xq2,wq2)
	

	open(unit=100,file=trim(dir_output)//'population_np-'//trim(file_output)//'.dat')
	open(unit=55,file=trim(dir_scratch)//'tsave-'//trim(file_output)//'.dat')
	open(unit=80,file=trim(dir_output)//'population-eigval_np-'//trim(file_output)//'.dat')
	cnp=c0
	do isave=0,num_save_vec
		allocate(yt(ntot1))
		write(cisave,'(i2.2)') isave
		call read_matrix_raw(yt,trim(dir_scratch)//'yt-'//trim(file_output)//'-'//cisave//'.dat')
	
		allocate(ytemp(ntot1))
		call hatomic_vec(m_over,yt,ytemp)
		deallocate(yt)
		ik=1
			do j=1,200
!			call zgemv('N',100,norder,c1,integ_cou,100,yt(norder+1:2*norder),1,c0,cnp,1)
				!if(eigenvalues(j).gt.-0.5d0)then
				if(eigenvalues(j).lt.1d0)then
				allocate(ynp(ntot1))
				ynp=c0
				ynp(ll*norder+1:(ll+1)*norder)=cmplx(ma(:,j))
				cnp(ik)=scalar_conjg(ynp,ytemp)
				ik=ik+1
				deallocate(ynp)	
					if(isave.eq.0)then
					write(80,'(f20.10)') eigenvalues(j)
					endif			
				endif
				!endif
			enddo
	read(55,*) nn,tsave
	write(100,'(f20.10,200es20.10)') tsave,abs(cnp)**2
	deallocate(ytemp)
	enddo	

	close(55)
	close(100)
	close(80)

	if(allocated(m_over)) deallocate(m_over)
	if(allocated(ytemp)) deallocate(ytemp)
	if(allocated(integ_cou)) deallocate(integ_cou)
	if(allocated(ma)) deallocate(ma)
	if(allocated(eigenvalues)) deallocate(eigenvalues)

	return
	end

!---------------------------------------------------------------------------------------------------------------------------------------
	subroutine population
	use datainput_module
	use complex_module
        use scratchdir_module
        use file_io_module
	use matrix_module
	implicit none
	real*8,allocatable :: h_atomic(:,:),m_over(:,:)
	complex*16 c1s,c2s,c2p,c3p,cnp,c3s,c3d,c4d
	complex*16,allocatable :: yt(:),ytemp(:),y1s(:),y2p(:),y3p(:),y2s(:),ynp(:,:),y3s(:),y3d(:),y4d(:)
	real*8 ma(norder,norder),cwork(1),mb(norder,norder)
	real*8 eigenvalues(norder),tsave
	integer info,lwork
	real*8,allocatable :: work(:)
	integer nn,ll,isave,i,j
      	character*2 cisave
	complex*16,external :: scalar_conjg

	allocate(h_atomic(ntot1,norder))
	allocate(m_over(ntot1,norder))
	call read_matrix_raw(h_atomic,trim(dir_scratch)//'hatomic-'//trim(file_output)//'.dat')
	call read_matrix_raw(m_over,trim(dir_scratch)//'mover-'//trim(file_output)//'.dat')

	open(unit=80,file=trim(dir_output)//'population-eigval-'//trim(file_output)//'.dat')

	!Population 1s
	ll=0
	nn=1

	ma(:,:)=h_atomic(ll*norder+1:(ll+1)*norder,:)
	mb(:,:)=m_over(ll*norder+1:(ll+1)*norder,:)

	lwork=-1
	call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,cwork,lwork,info)

	lwork=int(cwork(1))
	allocate(work(lwork))
	call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,work,lwork,info)
		if(info.ne.0)then
		write(6,*) 'info:',info,'in dsygv in groundstate'
		endif
        deallocate(work)

	allocate(y1s(ntot1))
	y1s=c0
	y1s(ll*norder+1:(ll+1)*norder)=cmplx(ma(:,nn))
	write(80,*) nn,ll,eigenvalues(nn)

	!Population 2s
	!ll=0
	!nn=2

	!ma(:,:)=h_atomic(ll*norder+1:(ll+1)*norder,:)
	!mb(:,:)=m_over(ll*norder+1:(ll+1)*norder,:)

	!lwork=-1
	!call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,cwork,lwork,info)

	!lwork=int(cwork(1))
	!allocate(work(lwork))
	!call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,work,lwork,info)
	!	if(info.ne.0)then
	!	write(6,*) 'info:',info,'in dsygv in groundstate'
	!	endif     
	!  deallocate(work)

	!allocate(y2s(ntot1))
	!y2s=c0
	!y2s(ll*norder+1:(ll+1)*norder)=cmplx(ma(:,nn))
	!write(80,*) nn,ll,eigenvalues(nn)

	!Population 3s
	!ll=0
	!nn=3

	!ma(:,:)=h_atomic(ll*norder+1:(ll+1)*norder,:)
	!mb(:,:)=m_over(ll*norder+1:(ll+1)*norder,:)

	!lwork=-1
	!call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,cwork,lwork,info)

	!lwork=int(cwork(1))
	!allocate(work(lwork))
	!call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,work,lwork,info)
	!	if(info.ne.0)then
	!	write(6,*) 'info:',info,'in dsygv in groundstate'
	!	endif
      !  deallocate(work)

	!allocate(y3s(ntot1))
	!y3s=c0
	!y3s(ll*norder+1:(ll+1)*norder)=cmplx(ma(:,nn))
	!write(80,*) nn,ll,eigenvalues(nn)

	!Population 2p
!	ll=1
!	nn=1

!	ma(:,:)=h_atomic(ll*norder+1:(ll+1)*norder,:)
!	mb(:,:)=m_over(ll*norder+1:(ll+1)*norder,:)

!	lwork=-1
!	call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,cwork,lwork,info)

!	lwork=int(cwork(1))
!	allocate(work(lwork))
!	call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,work,lwork,info)
!		if(info.ne.0)then
!		write(6,*) 'info:',info,'in dsygv in groundstate'
!		endif
!        deallocate(work)

!	allocate(y2p(ntot1))
!	y2p=c0
!	y2p(ll*norder+1:(ll+1)*norder)=cmplx(ma(:,nn))
!	write(80,*) nn,ll,eigenvalues(nn)

!	!Population 3p
!	ll=1
!	nn=2

!	ma(:,:)=h_atomic(ll*norder+1:(ll+1)*norder,:)
!	mb(:,:)=m_over(ll*norder+1:(ll+1)*norder,:)

!	lwork=-1
!	call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,cwork,lwork,info)

!	lwork=int(cwork(1))
!	allocate(work(lwork))
!	call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,work,lwork,info)
!		if(info.ne.0)then
!		write(6,*) 'info:',info,'in dsygv in groundstate'
!		endif
!        deallocate(work)

!	allocate(y3p(ntot1))
!	y3p=c0
!	y3p(ll*norder+1:(ll+1)*norder)=cmplx(ma(:,nn))
!	write(80,*) nn,ll,eigenvalues(nn)

!	!Population 3d
!	ll=2
!	nn=1

!	ma(:,:)=h_atomic(ll*norder+1:(ll+1)*norder,:)
!	mb(:,:)=m_over(ll*norder+1:(ll+1)*norder,:)

!	lwork=-1
!	call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,cwork,lwork,info)

!	lwork=int(cwork(1))
!	allocate(work(lwork))
!	call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,work,lwork,info)
!		if(info.ne.0)then
!		write(6,*) 'info:',info,'in dsygv in groundstate'
!		endif
!        deallocate(work)

!	allocate(y3d(ntot1))
!	y3d=c0
!	y3d(ll*norder+1:(ll+1)*norder)=cmplx(ma(:,nn))
!	write(80,*) nn,ll,eigenvalues(nn)

!	!Population 4d
!	ll=2
!	nn=2

!	ma(:,:)=h_atomic(ll*norder+1:(ll+1)*norder,:)
!	mb(:,:)=m_over(ll*norder+1:(ll+1)*norder,:)

!	lwork=-1
!	call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,cwork,lwork,info)

!	lwork=int(cwork(1))
!	allocate(work(lwork))
!	call dsygv(1,'V','L',norder,ma,norder,mb,norder,eigenvalues,work,lwork,info)
!		if(info.ne.0)then
!		write(6,*) 'info:',info,'in dsygv in groundstate'
!		endif
!        deallocate(work)

!	allocate(y4d(ntot1))
!	y4d=c0
!	y4d(ll*norder+1:(ll+1)*norder)=cmplx(ma(:,nn))
!	write(80,*) nn,ll,eigenvalues(nn)


	deallocate(h_atomic)

	open(unit=100,file=trim(dir_output)//'population1s-'//trim(file_output)//'.dat')
!	open(unit=200,file=trim(dir_output)//'population2p-'//trim(file_output)//'.dat')
!	open(unit=400,file=trim(dir_output)//'population3p-'//trim(file_output)//'.dat')
!	open(unit=500,file=trim(dir_output)//'population2s-'//trim(file_output)//'.dat')
!	open(unit=600,file=trim(dir_output)//'population3s-'//trim(file_output)//'.dat')
!	open(unit=700,file=trim(dir_output)//'population3d-'//trim(file_output)//'.dat')
!	open(unit=800,file=trim(dir_output)//'population4d-'//trim(file_output)//'.dat')
	open(unit=55,file=trim(dir_scratch)//'tsave-'//trim(file_output)//'.dat')
	do isave=0,num_save_vec
		allocate(yt(ntot1))
		write(cisave,'(i2.2)') isave
		call read_matrix_raw(yt,trim(dir_scratch)//'yt-'//trim(file_output)//'-'//cisave//'.dat')
	
		allocate(ytemp(ntot1))
		call hatomic_vec(m_over,yt,ytemp)
		deallocate(yt)
		c1s=scalar_conjg(y1s,ytemp)
!		c2s=scalar_conjg(y2s,ytemp)
!		c2p=scalar_conjg(y2p,ytemp)
!		c3p=scalar_conjg(y3p,ytemp)
!		c3s=scalar_conjg(y3s,ytemp)
!		c3d=scalar_conjg(y3d,ytemp)
!		c4d=scalar_conjg(y4d,ytemp)
		read(55,*) nn,tsave
		write(100,'(f20.10,f30.20)') tsave,abs(c1s)**2
!		write(200,'(f20.10,f30.20)') tsave,abs(c2p)**2 
!		write(400,'(f20.10,f30.20)') tsave,abs(c3p)**2 
!		write(500,'(f20.10,f30.20)') tsave,abs(c2s)**2
!		write(600,'(f20.10,f30.20)') tsave,abs(c3s)**2 
!		write(700,'(f20.10,f30.20)') tsave,abs(c3d)**2 
!		write(800,'(f20.10,f30.20)') tsave,abs(c4d)**2		

	deallocate(ytemp)
	enddo	

	deallocate(y1s) !,y2p,y3p,y2s,y3s,y3d,y4d)

	close(55)
	close(100)
!	close(200)
!	close(400)
!	close(500)
!	close(600)
!	close(700)
!	close(800)

	if(allocated(m_over)) deallocate(m_over)

	return
	end
!--------------------------------------------------------------------------------------------------------------------------
	subroutine plot_wp
	use datainput_module
	use complex_module
	use matrix_module
	use file_io_module
	use pi_module
	use scratchdir_module
	integer l,ir,i,ele,j,isave,iq,nn
	real*8 ang,vr,tsave
	complex*16 wp,ylm,c_i
	character*2 cisave,cl
	character*2 cn
	complex*16,allocatable :: sturmq(:,:),yt(:),st1(:,:,:),wpl(:),y0(:),vec(:),vrq(:,:)
	real*8,allocatable :: xq1(:),rgrid(:),m_over(:,:)
	real*8,external :: potential1
	complex*16,external :: scalar_conjg

	if(plot_vq.eq.1)then
	write(6,*) 'Plotting eigenvectors'
	allocate(dst(0:lq,norder,nquad))
	do ele=0,lq
		write(cl,'(i2.2)') ele
		allocate(sturmq(norder,nquad))
		call read_matrix_raw(sturmq,trim(dir_scratch)//'st-l'//cl//'-'//trim(file_output)//'.dat')
		dst(ele,:,:)=sturmq(:,:)
		if(allocated(sturmq)) deallocate(sturmq)
      	enddo

	allocate(xq1(nquad),st1(0:lq,norder,ngrids),rgrid(ngrids))
	call read_matrix_raw(xq1,trim(dir_scratch)//'xq-'//trim(file_output)//'.dat')
	do j=1,ngrids
		rgrid(j)=(j-1)*(rmax/dble(ngrids-1))
	enddo
	call interpolation(ngrids,rgrid,xq1,dst,st1)
	deallocate(dst,xq1)

	do isave=1,num_save_vec
	allocate(vrq(ntot1,0:nkry),wpl(0:lq))
	write(cisave,'(i2.2)') isave
	call read_matrix_raw(vrq,trim(dir_scratch)//'eigv_propag-'//trim(file_output)//'-'//cisave//'.dat')

	do j=1,3
	write(cn,'(i1.1)') j
	open(unit=100,file=trim(dir_output)//'eigvr-'//trim(file_output)//'-n'//cn//'-'//cisave//'.dat')
	
		do ir=1,ngrids
		wp=c0
			do l=0,lq
			wpl(l)=c0
				do i=1,norder
				wpl(l)=wpl(l)+vrq(l*norder+i,j)*st1(l,i,ir)
				enddo
			wp=wp+abs(wpl(l))**2
			enddo
		write(100,'(f20.10,es30.15)') rgrid(ir),abs(wp)
		enddo
	close(100)
	enddo	
	deallocate(vrq,wpl)
	enddo
	endif

	if(plot_rdens.eq.1)then
	write(6,*) 'Plotting radial density'
	allocate(dst(0:lq,norder,nquad))
	do ele=0,lq
		write(cl,'(i2.2)') ele
		allocate(sturmq(norder,nquad))
		call read_matrix_raw(sturmq,trim(dir_scratch)//'st-l'//cl//'-'//trim(file_output)//'.dat')
		dst(ele,:,:)=sturmq(:,:)
		if(allocated(sturmq)) deallocate(sturmq)
      	enddo

	allocate(xq1(nquad),st1(0:lq,norder,ngrids),rgrid(ngrids))
	call read_matrix_raw(xq1,trim(dir_scratch)//'xq-'//trim(file_output)//'.dat')
	do j=1,ngrids
		rgrid(j)=(j-1)*(rmax/dble(ngrids-1))
	enddo
	call interpolation(ngrids,rgrid,xq1,dst,st1)
	deallocate(dst,xq1)

	do isave=0,num_save_vec
	allocate(yt(ntot1),wpl(0:lq))
	write(cisave,'(i2.2)') isave
	call read_matrix_raw(yt,trim(dir_scratch)//'yt-'//trim(file_output)//'-'//cisave//'.dat')

	if(plot_extract.eq.1)then
	allocate(m_over(ntot1,norder))
	call read_matrix_raw(m_over,trim(dir_scratch)//'mover-'//trim(file_output)//'.dat')

		!do l=0,
		l=qnumber(2)
		write(cl,'(i2.2)') l
			!do i=1,2
			i=qnumber(1)
			write(cn,'(i2.2)') i
			allocate(y0(ntot1))
			call read_matrix_raw(y0,trim(dir_scratch)//'yt-n'//cn//'l'//cl//'-'//trim(file_output)//'.dat')
			allocate(vec(ntot1))
			call hatomic_vec(m_over,yt,vec)
			c_i=scalar_conjg(y0,vec)
			deallocate(vec)
			yt=yt-c_i*y0
			deallocate(y0)
			!enddo
		!enddo
	deallocate(m_over)
	endif		

	open(unit=100,file=trim(dir_output)//'rdens-'//trim(file_output)//'-'//cisave//'.dat')
	
		do ir=1,ngrids
		wp=c0
			do l=0,lq
			wpl(l)=c0
				do i=1,norder
				wpl(l)=wpl(l)+yt(l*norder+i)*st1(l,i,ir)
				enddo
			wp=wp+abs(wpl(l))**2d0
			enddo
		write(100,'(f20.10,5es30.15)') rgrid(ir),real(wp),real(wpl(1)),aimag(wpl(1)),real(wpl(0)),aimag(wpl(0))
		enddo
	close(100)
	deallocate(yt,wpl)
	enddo	
	endif

	if(plot_wavepkt.eq.1)then
	write(6,*) 'Plotting wave packet'
	if(allocated(dst)) deallocate(dst)
	if(allocated(st1)) deallocate(st1)
	if(allocated(xq1)) deallocate(xq1)
	if(allocated(rgrid)) deallocate(rgrid)
	allocate(dst(0:lq,norder,nquad))
	do ele=0,lq
		write(cl,'(i2.2)') ele
		allocate(sturmq(norder,nquad))
		call read_matrix_raw(sturmq,trim(dir_scratch)//'st-l'//cl//'-'//trim(file_output)//'.dat')
		dst(ele,:,:)=sturmq(:,:)
		if(allocated(sturmq)) deallocate(sturmq)
      	enddo

	allocate(xq1(nquad),st1(0:lq,norder,ngrids),rgrid(ngrids))
	call read_matrix_raw(xq1,trim(dir_scratch)//'xq-'//trim(file_output)//'.dat')
	do j=1,ngrids
		rgrid(j)=j*(rmax/dble(ngrids+1))
	enddo
	call interpolation(ngrids,rgrid,xq1,dst,st1)
	deallocate(dst,xq1)

	ang=2d0*pi/(100d0+1d0)
	do isave=0,num_save_vec
	allocate(yt(ntot1))
	write(cisave,'(i2.2)') isave
	call read_matrix_raw(yt,trim(dir_scratch)//'yt-'//trim(file_output)//'-'//cisave//'.dat')

	if(plot_extract.eq.1)then
	allocate(m_over(ntot1,norder))
	call read_matrix_raw(m_over,trim(dir_scratch)//'mover-'//trim(file_output)//'.dat')

		do l=0,1
		write(cl,'(i2.2)') l
			do i=1,2
			write(cn,'(i2.2)') i
			allocate(y0(ntot1))
			call read_matrix_raw(y0,trim(dir_scratch)//'yt-n'//cn//'l'//cl//'-'//trim(file_output)//'.dat')
			allocate(vec(ntot1))
			call hatomic_vec(m_over,yt,vec)
			c_i=scalar_conjg(y0,vec)
			deallocate(vec)
			yt=yt-c_i*y0
			deallocate(y0)
			enddo
		enddo
	deallocate(m_over)
	endif	

	open(unit=100,file=trim(dir_output)//'wp-'//trim(file_output)//'-'//cisave//'.dat')
		
		do ir=1,ngrids
		do iq=0,101
		wp=c0
			do l=0,lq
			call spherical_harmonic(l,0,iq*ang,0d0,ylm)
				do i=1,norder
				wp=wp+yt(l*norder+i)*st1(l,i,ir)*ylm
				enddo
			enddo
		write(100,'(f20.10,f20.10,es30.15)') rgrid(ir)*cos(iq*ang),rgrid(ir)*sin(iq*ang),abs(wp)**2
		enddo
		write(100,*)
		enddo
	close(100)
	deallocate(yt)
	enddo	
	endif

	if(allocated(dst)) deallocate(dst)
	if(allocated(st1)) deallocate(st1)
	if(allocated(xq1)) deallocate(xq1)
	if(allocated(rgrid)) deallocate(rgrid)

	return
	end
!---------------------------------------------------------------------------
	subroutine plot_flux
	use datainput_module
	use complex_module
	use matrix_module
	use file_io_module
	use pi_module
	use scratchdir_module
	integer l,ir,i,ele,j,isave,iq,nn
	real*8 ang,vr,tsave
	complex*16 flx
	character*2 cisave,cl,cn
	complex*16,allocatable :: sturmq(:,:),st1(:,:,:),dst1(:,:,:),yt(:),wpl(:),dwpl(:),flxl(:)
	real*8,allocatable :: xq1(:),rgrid(:)
	complex*16,external :: scalar_conjg

	if(plot_rdens.eq.2)then
	write(6,*) 'Plotting radial flux'
	allocate(st(0:lq,norder,nquad),dst(0:lq,norder,nquad))
	do ele=0,lq
		write(cl,'(i2.2)') ele
		allocate(sturmq(norder,nquad))
		call read_matrix_raw(sturmq,trim(dir_scratch)//'st-l'//cl//'-'//trim(file_output)//'.dat')
		st(ele,:,:)=sturmq(:,:)
		call read_matrix_raw(sturmq,trim(dir_scratch)//'dst-l'//cl//'-'//trim(file_output)//'.dat')
		dst(ele,:,:)=sturmq(:,:)
		if(allocated(sturmq)) deallocate(sturmq)
      	enddo

	allocate(xq1(nquad),st1(0:lq,norder,ngrids),dst1(0:lq,norder,ngrids),rgrid(ngrids))
	call read_matrix_raw(xq1,trim(dir_scratch)//'xq-'//trim(file_output)//'.dat')
	do j=1,ngrids
		rgrid(j)=j*(rmax/dble(ngrids+1))
	enddo
	call interpolation(ngrids,rgrid,xq1,st,st1)
	call interpolation(ngrids,rgrid,xq1,dst,dst1)
	deallocate(st,dst,xq1)


	do isave=0,num_save_vec
	allocate(yt(ntot1),wpl(0:lq),dwpl(0:lq),flxl(0:lq))
	write(cisave,'(i2.2)') isave
	call read_matrix_raw(yt,trim(dir_scratch)//'yt-'//trim(file_output)//'-'//cisave//'.dat')

	open(unit=100,file=trim(dir_output)//'flux-'//trim(file_output)//'-'//cisave//'.dat')
	
		do ir=1,ngrids
		flxl=c0
		flx=c0
			do l=0,lq
			wpl(l)=c0
			dwpl(l)=c0
				do i=1,norder
				wpl(l)=wpl(l)+yt(l*norder+i)*st1(l,i,ir)/rgrid(ir)
				dwpl(l)=dwpl(l)+yt(l*norder+i)*(-1d0*st1(l,i,ir)/rgrid(ir)**2+dst1(l,i,ir)/rgrid(ir))
				enddo
			flxl(l)=(1d0/(2d0*ci))*rgrid(ir)**2*(conjg(wpl(l))*dwpl(l)-conjg(dwpl(l))*wpl(l))
			flx=flx+flxl(l)
			enddo
		!write(*,*) rgrid(ir),flxl(1)
		write(100,'(f20.10,3es30.15)') rgrid(ir),real(flx),real(flxl(0)),real(flxl(1))
		enddo
	close(100)
	deallocate(yt,wpl,dwpl,flxl)
	enddo	
	endif

	if(allocated(st)) deallocate(st)
	if(allocated(dst1)) deallocate(dst1)
	if(allocated(dst)) deallocate(dst)
	if(allocated(st1)) deallocate(st1)
	if(allocated(xq1)) deallocate(xq1)
	if(allocated(rgrid)) deallocate(rgrid)

	return
	end
!---------------------------------------------------------------------------
	subroutine spect1(integ_cou)
	use datainput_module
	use complex_module
	use matrix_module
	use file_io_module
	use scratchdir_module
	use file_io_module
	use pi_module
	integer j,l,i,ik
	complex*16 integ_cou(ngridsk,0:lq,norder),ck1
	complex*16,allocatable :: yt(:),ckl(:),cki(:)
	character*2 csave
	real*8 kmax,dk,ck,en,den,kk,integ
!	real*8,allocatable :: h_atomic(:,:),m_over(:,:),ma(:,:),mb(:,:),work(:),temp(:,:)
	integer info,lwork
	real*8 cwork

	den=emax/dble(ngridsk+1)

	!do j=0,num_save_vec
	!write(csave,'(i2.2)') j

	allocate(yt(ntot1))
	call read_matrix_raw(yt,trim(dir_scratch)//'yt-'//trim(file_output)//'-tfin_pulse.dat')
	open(unit=200,file=trim(dir_output)//'spect-'//trim(file_output)//'-tfin_pulse.dat')
	!call read_matrix_raw(yt,trim(dir_scratch)//'yt-'//trim(file_output)//'-'//csave//'.dat')
	!open(unit=200,file=trim(dir_output)//'spect-'//trim(file_output)//'-'//csave//'.dat')
!	open(unit=300,file=trim(dir_output)//'integrated_spect-'//trim(file_output)//'-tfin_pulse.dat')


	allocate(ckl(0:lq),cki(1:ngridsk))
	cki=0d0
	ckl=c0
	do ik=1,ngridsk
	ck=0d0
	en=ik*den
	kk=sqrt(2d0*en)
		do l=0,lq
		ckl(l)=c0
			do i=1,norder
			ckl(l)=ckl(l)+yt(l*norder+i)*integ_cou(ik,l,i)
			enddo
		ck=ck+abs(ckl(l))**2
		enddo
	cki(ik)=ck
	if(en.gt.0.0001d0) write(200,'(f20.10,6es30.15)') en,ck,abs(ckl(0))**2,abs(ckl(1))**2,abs(ckl(2))**2,abs(ckl(3))**2,abs(ckl(4))**2
	enddo	

	close(200)	

!	integ=0d0
!	do ik=2,ngridsk
!	integ=integ+(cki(ik)+cki(i-1))*den/2d0
!	enddo

!	write(300,'(2f20.10)') omega,integ
!	close(300)

	if(allocated(ckl)) deallocate(ckl)
	if(allocated(cki)) deallocate(cki)
	if(allocated(yt)) deallocate(yt)
!	enddo

	return
	end
!--------------------------------------------------------------------------------------------------------------------------
	subroutine plot_rmean
	use datainput_module
	use complex_module
	use matrix_module
	use file_io_module
	use pi_module
	use scratchdir_module
	integer l,ir,i,ele,j,isave,iq,nquad2,nn
	complex*16 wp,ylm,rm,rml(0:lq),wpl(0:lq)
	character*2 cisave,cl
	complex*16,allocatable :: sturmq(:,:),yt(:),dsti(:,:,:)
	real*8,allocatable :: xq1(:),wq1(:),wq2(:),xq2(:),xsi(:),xs3(:),xs4(:)
	real*8 a0,a1,tsave

	if(allocated(dst)) deallocate(dst)
	allocate(xq1(nquad),wq1(nquad))
	call read_matrix_raw(xq1,trim(dir_scratch)//'xq-'//trim(file_output)//'.dat')
	!call read_matrix_raw(xq1,trim(dir_scratch)//'xsmid-'//trim(file_output)//'.dat')
	call read_matrix_raw(wq1,trim(dir_scratch)//'wq-'//trim(file_output)//'.dat')
	nquad2=100
	allocate(xq2(nquad2),wq2(nquad2),xsi(nquad2),xs3(nquad2),xs4(nquad2))
	call evaluate_quad_sturmians(nquad2,pi,xq2,wq2,rc0,delta0,xsi,rc1,delta1,xs3,rc2,delta2,xs4)
	deallocate(xsi,xs3,xs4)

	allocate(dst(0:lq,norder,nquad))
	do ele=0,lq
		write(cl,'(i2.2)') ele
		allocate(sturmq(norder,nquad))
		call read_matrix_raw(sturmq,trim(dir_scratch)//'st-l'//cl//'-'//trim(file_output)//'.dat')
		!call read_matrix_raw(sturmq,trim(dir_scratch)//'ssmid-l'//cl//'-'//trim(file_output)//'.dat')
		dst(ele,:,:)=sturmq(:,:)
		if(allocated(sturmq)) deallocate(sturmq)
      enddo

	open(unit=55,file=trim(dir_scratch)//'tsave-'//trim(file_output)//'.dat')
	open(unit=100,file=trim(dir_output)//'rmeanmid-'//trim(file_output)//'.dat')

	do isave=0,num_save_vec
	allocate(yt(ntot1))
	write(cisave,'(i2.2)') isave
	call read_matrix_raw(yt,trim(dir_scratch)//'yt-'//trim(file_output)//'-'//cisave//'.dat')
	
	a0=rmax/2d0
	!a0=(rc1-rc0-delta0)/2d0
	a1=pi/2d0
	rm=c0
	rml=c0
		do ir=1,nquad
		do iq=1,nquad2
		wp=c0
			do l=0,lq
			wpl(l)=c0
			call spherical_harmonic(l,0,xq2(iq),0d0,ylm)
				do i=1,norder
				wp=wp+yt(l*norder+i)*dst(l,i,ir)*ylm
				wpl(l)=wpl(l)+yt(l*norder+i)*dst(l,i,ir)*ylm
				enddo
			rml(l)=rml(l)+2d0*pi*wq1(ir)*wq2(iq)*sin(xq2(iq))*xq1(ir)*abs(wpl(l))**2*a0*a1
			!rml(l)=rml(l)+2d0*pi*wq1(ir)*wq2(iq)*sin(xq2(iq))*abs(wpl(l))**2*a0*a1
			enddo
		rm=rm+2d0*pi*wq1(ir)*wq2(iq)*sin(xq2(iq))*xq1(ir)*abs(wp)**2*a0*a1
		!rm=rm+2d0*pi*wq1(ir)*wq2(iq)*sin(xq2(iq))*abs(wp)**2*a0*a1
		enddo
		enddo
	read(55,*) nn,tsave
	write(100,'(f20.10,4es30.15)') tsave,real(rm),real(rml(0)),real(rml(1)),real(rml(2))

	deallocate(yt)
	enddo	

	close(55)
	close(100)

	return
	end
!--------------------------------------------------------------------
	subroutine spherical_harmonic(l,m,theta,phi,ylm)
!***********************************************************************************
!	This subroutine evaluates the regular spherical harmonics
!	y(l,m) in angles (theta,phi) with the ANGLES IN RADIANS
!	subroutine by Dr. A L Frapiccini
!***********************************************************************************
	use complex_module
	use pi_module
	implicit none
	integer l,m
	real*8 theta,phi,plm
	real*8,external :: plgndr
	complex*16 ylm,gam1,gam2

	call gammaln(c1*(l-m+1d0),gam1)
	call gammaln(c1*(l+m+1d0),gam2)
	plm=plgndr(l,m,cos(theta))

	ylm=plm*exp(ci*m*phi)
	ylm=ylm*sqrt((2d0*l+1d0)/4d0/pi)*sqrt(exp(gam1-gam2))

	return
	end
!-------------------------------------------------------------------------------------
      real*8 FUNCTION plgndr(l,m,x)
!	Computes the associated Legendre polynomial P(l,m)(x). 
!	Here m and l are integers satisfying 0 <=m<=l, while x lies in the range −1<=x<=1.
!	From Numerical Recipes in Fortran 77
      INTEGER l,m
      REAL*8 x
      INTEGER i,ll
      REAL*8 fact,pll,pmm,pmmp1,somx2

      pmm=1d0
      if(m.gt.0) then
        somx2=sqrt((1d0-x)*(1d0+x))
        fact=1d0
        do i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2d0
       enddo
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2d0*m+1d0)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do ll=m+2,l
            pll=(x*(2d0*ll-1d0)*pmmp1-(ll+m-1d0)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
	  enddo
          plgndr=pll
        endif
      endif
      return
      END
!--------------------------------------------------------------------------------
        subroutine factt
!********************************************************************************
!     Calculate the log of factorials of all integers from 1 to nfact
!     and store them in the array gamma.
!     Factorials up to 120 are stored in fact.
!********************************************************************************
      	implicit real*8(a-h,o-z)
	parameter(nfact=800)
      	common/fact/gamma(nfact)
      	common/bkfactorials/fact(0:121),fact2(0:121)
	dimension gamma2(nfact)
	data rzero,one/0.0d0,1.0d0/

      	gamma=rzero
      	gamma(2)=rzero
       	fact(0)=one
       	fact=one
      	gamma2=rzero
      	gamma2(2)=rzero
       	fact2(0)=one
       	fact2=one

      	do 10 k=3,nfact
       		j=k-1
       		a=dfloat(j)
       		gamma(k)=gamma(j)+dlog(a)
		iodd = mod(j,2)
		if (iodd.eq.1) gamma2(k) = gamma2(k-2)+dlog(a)
       		if (k.le.121) fact(k-1)=exp(gamma(k))
       		if (k.le.121) fact2(k-1)=exp(gamma2(k))
10    	continue
     
      	return
      	end

