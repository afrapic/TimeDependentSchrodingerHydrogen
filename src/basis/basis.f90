	subroutine build_basis
	use matrix_module
	use datainput_module
	use scratchdir_module
	use file_io_module
	integer ele,ik,l,i,j
	complex*16,allocatable :: eigvalues(:),sturmq(:,:),dsturmq(:,:),func(:,:),coef(:,:),sturmq1(:,:),sturms0(:,:),sturms1(:,:),sturms2(:,:)
	complex*16 en
	character*3 ck
	character*1 ce
	character*2 cl

	l=nb-kp+1
	do ele=0,lq
	write(cl,'(i2.2)') ele
	allocate(eigvalues(norder))
	allocate(sturmq(norder,nquad))
	allocate(dsturmq(norder,nquad),coef(nb,norder))
	allocate(sturms0(norder,nquad))
	allocate(sturms1(norder,nquad))
	allocate(sturms2(norder,nquad))
	call basis(ele,norder,nquad,nb,l,kp,tx,ilong,ishort,alpha,charge,u0,rc0,delta0,u1,rc1,delta1,&
		u2,rc2,delta2,iboundary,rmax,nquad_bsp,energy,eigvalues,coef,xq,sturmq,dsturmq,xs0,sturms0,xs1,sturms1,xs2,sturms2)
	call write_matrix_raw(eigvalues,trim(dir_scratch)//'eigenvalues-l'//cl//'-'//trim(file_output)//'.dat')
	call write_matrix_raw(sturmq,trim(dir_scratch)//'st-l'//cl//'-'//trim(file_output)//'.dat')
	call write_matrix_raw(dsturmq,trim(dir_scratch)//'dst-l'//cl//'-'//trim(file_output)//'.dat')
		if(ilong.eq.5)then
		call write_matrix_raw(sturms0,trim(dir_scratch)//'ss0-l'//cl//'-'//trim(file_output)//'.dat')
		endif
		if(ilong.eq.6)then
		call write_matrix_raw(sturms0,trim(dir_scratch)//'ss0-l'//cl//'-'//trim(file_output)//'.dat')
		call write_matrix_raw(sturms1,trim(dir_scratch)//'ss1-l'//cl//'-'//trim(file_output)//'.dat')
		endif
		if(ilong.eq.7)then
		call write_matrix_raw(sturms0,trim(dir_scratch)//'ss0-l'//cl//'-'//trim(file_output)//'.dat')
		call write_matrix_raw(sturms1,trim(dir_scratch)//'ss1-l'//cl//'-'//trim(file_output)//'.dat')
		call write_matrix_raw(sturms2,trim(dir_scratch)//'ss2-l'//cl//'-'//trim(file_output)//'.dat')
		endif
	if(allocated(eigvalues)) deallocate(eigvalues)
	if(allocated(sturmq)) deallocate(sturmq)
	if(allocated(sturms0)) deallocate(sturms0)
	if(allocated(sturms1)) deallocate(sturms1)
	if(allocated(sturms2)) deallocate(sturms2)
	if(allocated(dsturmq)) deallocate(dsturmq)
	if(allocated(coef)) deallocate(coef)
     	enddo

	return
	end
!------------------------------------------------------------------------------------------------------------------------------
	subroutine basis(ele,norder,nquad,nb,l,kp,tx,ilong,ishort,alpha,charge,u0,rc0,delta0,u1,rc1,delta1,u2,rc2,delta2,iboundary,rmax,nquad_bsp,energy,eigval &
 			,coef,xq,sturmq,dsturmq,xs,sturms,xs1,sturms1,xs2,sturms2)
!***********************************************************************************
!	Solves the generalized eigenvalue problem obtained by using the Bspline
!	expansion of the Sturmians and projecting onto the Bsplines for a given
!	electron nel and angular momentum ele
!	The output are the eigenvalues (eigval), the eigenvectors (coef),
!	the Sturmians evaluated in xq and their derivatives (sturmq and dsturmq),
!	and the Sturmians evalueated in xs.
!***********************************************************************************
	use complex_module
	use file_io_module
	use scratchdir_module
	implicit none
	integer ele,ncase,norder,nb,l,kp,ilong,ishort,iboundary,nquad_bsp,nquad,ierr
	real*8 alpha,charge,rmax,xq(nquad),u0,rc0,delta0,u1,rc1,delta1,u2,rc2,delta2,xs(nquad),xs1(nquad),xs2(nquad)
	complex*16 eigval(norder),coef(nb,norder),energy
	complex*16,allocatable :: mh(:,:),mv(:,:),mc(:,:)
	complex*16,allocatable :: alp(:),bet(:),vl(:,:),vr(:,:),work(:)
	real*8,allocatable :: dmh(:,:),dmv(:,:),w(:),dvr(:,:),dvl(:,:),dalpr(:),dalpi(:),dbet(:)
	real*8,allocatable :: rwork(:)
	integer info,lwork
	complex*16 cwork(1),yn,ynp,ratio
	integer lda,j,i
	complex*16,allocatable :: eigenvalues(:),eigenvectors(:,:),nm(:),func(:,:)
	complex*16 sturmq(norder,nquad),sturms(norder,nquad),sturms1(norder,nquad),sturms2(norder,nquad)
	complex*16 dsturmq(norder,nquad)
	real*8 crwork(1),tx(l+2*kp-1)
	integer,allocatable :: ipiv(:)
	character*2 cl
	
	!Matrices h and v with the boundary conditions included
	allocate(mh(1:nb-2,1:nb-2),mv(1:nb-2,1:nb-2))
	call hamiltonian(ele,kp,nb,l,tx,rmax,ilong,u0,rc0,delta0,u1,rc1,delta1,u2,rc2,delta2,charge,energy,nquad_bsp,iboundary,mh)
	call srpotential(kp,nb,l,tx,rmax,energy,ishort,alpha,charge,u0,rc0,delta0,u1,rc1,delta1,u2,rc2,delta2,nquad_bsp,iboundary,mv)

	!Solves the eigenproblem
	!Chooses the lapack subroutine suitable for the energy+boundary case
	if(real(energy).gt.0d0)then
	if(iboundary.ne.0)then
	ncase=0
	else
	ncase=1
	endif
	endif
	if(real(energy).lt.0d0)then
	if(iboundary.ne.0)then
	ncase=2
	else
	ncase=1
	endif
	endif

	select case(ncase)
	!General complex
	case(0)
	allocate(ipiv(nb-2))
	call zgetrf(nb-2,nb-2,mv,nb-2,ipiv,info)
	lwork=nb-2
	allocate(work(lwork))
	call zgetri(nb-2,mv,nb-2,ipiv,work,lwork,info)
	deallocate(work,ipiv)
	allocate(mc(nb-2,nb-2))
	call zgemm('N','N',nb-2,nb-2,nb-2,c1,mv,nb-2,mh,nb-2,c0,mc,nb-2)
	deallocate(mv,mh)
	allocate(vl(nb-2,nb-2),vr(nb-2,nb-2),eigenvalues(nb-2),rwork(2*(nb-2)))
	lwork=-1
	call zgeev('N','V',nb-2,mc,nb-2,eigenvalues,vl,nb-2,vr,nb-2,cwork,lwork,rwork,info)
	lwork=int(cwork(1))
	allocate(work(lwork))
	call zgeev('N','V',nb-2,mc,nb-2,eigenvalues,vl,nb-2,vr,nb-2,work,lwork,rwork,info)
	if(info.ne.0)then
	write(6,*) 'info:',info,'in zgeev in basis.f90'
	stop
	endif
	deallocate(vl,rwork,work,mc)
	!Real symmetric
	case(1)
	allocate(dmh(nb-2,nb-2),dmv(nb-2,nb-2))
	dmh=real(mh)
	dmv=real(mv)
	deallocate(mh,mv)
	allocate(w(nb-2))
	lwork=-1
	call dsygv(1,'V','U',nb-2,dmh,nb-2,dmv,nb-2,w,crwork,lwork,info)
	lwork=int(crwork(1))
	allocate(rwork(lwork))
	call dsygv(1,'V','U',nb-2,dmh,nb-2,dmv,nb-2,w,rwork,lwork,info)
	if(info.ne.0)then
	write(6,*) 'info:',info,'in dsygv in basis.f90'
	stop
	endif
	deallocate(rwork,dmv)
	allocate(eigenvalues(nb-2))
	eigenvalues=cmplx(w)
	deallocate(w)
	allocate(vr(nb-2,nb-2))
	vr=cmplx(dmh)
	deallocate(dmh)
	!General real
	case(2)
	allocate(dmh(nb-2,nb-2),dmv(nb-2,nb-2))
	dmh=real(mh)
	dmv=real(mv)
	deallocate(mh,mv)
	allocate(dalpr(nb-2),dalpi(nb-2),dbet(nb-2),dvr(nb-2,nb-2),dvl(nb-2,nb-2))
	lwork=-1
	call dggev('N','V',nb-2,dmh,nb-2,dmv,nb-2,dalpr,dalpi,dbet,dvl,nb-2,dvr,nb-2,crwork,lwork,info)
	lwork=int(crwork(1))
	allocate(rwork(lwork))
	call dggev('N','V',nb-2,dmh,nb-2,dmv,nb-2,dalpr,dalpi,dbet,dvl,nb-2,dvr,nb-2,rwork,lwork,info)
	if(info.ne.0)then
	write(6,*) 'info:',info,'in dggev in basis.f90'
	stop
	endif
	deallocate(dmh,dmv,dvl,rwork)
	allocate(eigenvalues(nb-2))
	eigenvalues=(dalpr+ci*dalpi)/dbet
	deallocate(dalpr,dalpi,dbet)
	allocate(vr(nb-2,nb-2))
	vr=cmplx(dvr)
	deallocate(dvr)
	end select
	
	!Sorts eigenvalues and coefficients in increasing absolute value
	lda=size(vr,2)
	call sort(lda,eigenvalues,vr)

	!Coefficient eigenvectors(nb,:) according to the boundary condition
	!eigenvectors(j,#eigenvalue)
	allocate(eigenvectors(1:nb,1:nb-2))
	eigenvectors=c0
	if(iboundary.ne.0)then
	eigenvectors(1,1:nb-2)=c0
	call couldecay(charge,energy,rmax,iboundary,yn,ynp)
	ratio=ynp/yn
		do j=1,nb-2
		eigenvectors(2:nb-1,j)=yn*(kp*1d0-1d0-ratio*(rmax-tx(nb)))*vr(1:nb-2,j)/vr(nb-2,j)/(kp*1d0-1d0)
		enddo
		eigenvectors(nb,1:nb-2)=yn		
	else
	eigenvectors(1,1:nb-2)=c0
	eigenvectors(2:nb-1,1:nb-2)=vr(1:nb-2,1:nb-2)
	eigenvectors(nb,1:nb-2)=c0
	endif
	deallocate(vr)

	allocate(nm(1:nb-2))
	call norm(rmax,charge,energy,l,kp,nb,tx,nquad_bsp,eigenvectors,nm)

	!Saves the norder(nel) asked eigenvalues and normalized eigenvectors
	eigval=c0
	coef=c0
	do j=1,norder
	eigval(j)=eigenvalues(j)
	coef(:,j)=eigenvectors(:,j)*nm(j)
	print*,j,eigval(j)
	enddo
	stop
	deallocate(nm,eigenvectors,eigenvalues)


	!Calculates Sturmians in the Legendre quadrature and in between the quadrature points
	do i=1,nquad
	allocate(func(2,1:norder))
	call radial_sturmians(l,kp,nb,tx,rmax,norder,xq(i),coef,2,func)
	sturmq(1:norder,i)=func(1,1:norder)
	dsturmq(1:norder,i)=func(2,1:norder)
	deallocate(func)
	enddo

	sturms=c0
	sturms1=c0
	if(ilong.eq.5)then
	!Calculates Sturmians in the Legendre quadrature and in between the quadrature points
	do i=1,nquad
	allocate(func(2,1:norder))
	call radial_sturmians(l,kp,nb,tx,rmax,norder,xs(i),coef,2,func)
	sturms(1:norder,i)=func(1,1:norder)
	deallocate(func)
	enddo
	endif

	if(ilong.eq.6)then
	!Calculates Sturmians in the Legendre quadrature and in between the quadrature points
	do i=1,nquad
	allocate(func(2,1:norder))
	call radial_sturmians(l,kp,nb,tx,rmax,norder,xs(i),coef,2,func)
	sturms(1:norder,i)=func(1,1:norder)
	deallocate(func)
	enddo
	do i=1,nquad
	allocate(func(2,1:norder))
	call radial_sturmians(l,kp,nb,tx,rmax,norder,xs1(i),coef,2,func)
	sturms1(1:norder,i)=func(1,1:norder)
	deallocate(func)
	enddo
!!!!!!!
	allocate(w(nquad),rwork(nquad))
	call evaluate_quad_sturmians1(nquad,0d0,rc0,rwork,w)
	call write_matrix_raw(rwork,trim(dir_scratch)//'xscou-'//trim(file_output)//'.dat')
	allocate(mv(norder,nquad))
	do i=1,nquad
	allocate(func(2,1:norder))
	call radial_sturmians(l,kp,nb,tx,rmax,norder,rwork(i),coef,2,func)
	mv(1:norder,i)=func(1,1:norder)
	deallocate(func)
	enddo
	write(cl,'(i2.2)') ele
	call write_matrix_raw(mv,trim(dir_scratch)//'sscou-l'//cl//'-'//trim(file_output)//'.dat')
	deallocate(w,rwork,mv)
	allocate(w(nquad),rwork(nquad))
	call evaluate_quad_sturmians1(nquad,rc0+delta0,rc1,rwork,w)
	call write_matrix_raw(rwork,trim(dir_scratch)//'xsmid-'//trim(file_output)//'.dat')
	allocate(mv(norder,nquad))
	do i=1,nquad
	allocate(func(2,1:norder))
	call radial_sturmians(l,kp,nb,tx,rmax,norder,rwork(i),coef,2,func)
	mv(1:norder,i)=func(1,1:norder)
	deallocate(func)
	enddo
	write(cl,'(i2.2)') ele
	call write_matrix_raw(mv,trim(dir_scratch)//'ssmid-l'//cl//'-'//trim(file_output)//'.dat')
	deallocate(w,rwork,mv)
	endif

	if(ilong.eq.7)then
	!Calculates Sturmians in the Legendre quadrature and in between the quadrature points
	do i=1,nquad
	allocate(func(2,1:norder))
	call radial_sturmians(l,kp,nb,tx,rmax,norder,xs(i),coef,2,func)
	sturms(1:norder,i)=func(1,1:norder)
	deallocate(func)
	enddo
	do i=1,nquad
	allocate(func(2,1:norder))
	call radial_sturmians(l,kp,nb,tx,rmax,norder,xs1(i),coef,2,func)
	sturms1(1:norder,i)=func(1,1:norder)
	deallocate(func)
	enddo
	do i=1,nquad
	allocate(func(2,1:norder))
	call radial_sturmians(l,kp,nb,tx,rmax,norder,xs2(i),coef,2,func)
	sturms2(1:norder,i)=func(1,1:norder)
	deallocate(func)
	enddo
	endif

	return
	end subroutine 
!-----------------------------------------------------------------------------------------------	
	subroutine hamiltonian(ele,kp,nb,l,tx,rmax,ilong,u0,rc0,delta0,u1,rc1,delta1,u2,rc2,delta2,charge,energy,nquad_bsp,iboundary,mh)
!********************************************************************************
!	Builds the matrix <B(j,k)|T+V0-E|B(i,k)> where T is the kynetic energy
!       and V0=-charge/r. 
!********************************************************************************
	use complex_module
	implicit none
	integer ele,kp,nb,l,nquad_bsp,ilong
	integer k,j,i,iboundary
	complex*16 mh(1:nb-2,1:nb-2)
	complex*16,allocatable :: du_kin(:,:),du_l(:,:),du_cou(:,:),du_ov(:,:),du_h(:,:)
	complex*16 yn,ynp,ratio,energy
	real*8 rmax,tx(l+2*kp-1),charge,u0,rc0,delta0,u1,rc1,delta1,u2,rc2,delta2

	allocate(du_kin(2*kp-1,nb-1),du_l(2*kp-1,nb-1),du_cou(2*kp-1,nb-1),du_ov(2*kp-1,nb-1))
	call integ_kin(l,kp,nb,tx,nquad_bsp,du_kin)
	call integ_l(ele,l,kp,nb,tx,nquad_bsp,du_l)
	call integ_cou(ilong,u0,rc0,delta0,u1,rc1,delta1,u2,rc2,delta2,charge,rmax,l,kp,nb,tx,nquad_bsp,du_cou)
	call integ_over(l,kp,nb,tx,nquad_bsp,du_ov)

	allocate(du_h(2*kp-1,nb-1))
	du_h=c0
	du_h=du_kin+du_l+charge*du_cou-energy*du_ov
	deallocate(du_kin,du_l,du_cou,du_ov)	

	if(iboundary.ne.0)then
	call couldecay(charge,energy,rmax,iboundary,yn,ynp)
	ratio=ynp/yn
	do i=2,kp
	du_h(i,nb-2)=du_h(i,nb-2)+du_h(i-1,nb-1)*(kp-1d0)/(kp-1d0-ratio*(rmax-tx(nb)))
	enddo
	endif

	mh=c0
	do j=1,nb-2
	k=kp-j
		do i=max(1,j-kp+1),min(nb-2,j+kp-1)
		mh(i,j)=du_h(k+i,j)
		enddo
	enddo
	deallocate(du_h)

	if(allocated(du_kin)) deallocate(du_kin)  
	if(allocated(du_l)) deallocate(du_l)
	if(allocated(du_cou)) deallocate(du_cou)
	if(allocated(du_ov)) deallocate(du_ov)

	return
	end subroutine 
!----------------------------------------------------------------------------------------------------	
	subroutine srpotential(kp,nb,l,tx,rmax,energy,ishort,alpha,charge,u0,rc0,delta0,u1,rc1,delta1,u2,rc2,delta2,nquad_bsp,iboundary,v)
!***************************************************************************************
!	Builds the matrix <B(j,k)|V|B(i,k)> where V is the short range potential 
!***************************************************************************************
	use complex_module
	implicit none
	integer kp,nb,l,ishort,nquad_bsp
	integer i,j,k,iboundary
	complex*16 v(1:nb-2,1:nb-2)
	complex*16 ratio,yn,ynp,energy
	complex*16,allocatable :: du_pot(:,:)
	real*8 rmax,tx(l+2*kp-1),alpha,charge
	real*8 u0,rc0,delta0,u1,rc1,delta1,u2,rc2,delta2

	allocate(du_pot(2*kp-1,nb-1))
	call integ_pot(l,kp,nb,tx,rmax,ishort,alpha,charge,u0,rc0,delta0,u1,rc1,delta1,u2,rc2,delta2,nquad_bsp,du_pot)
	
	if(iboundary.ne.0)then
	call couldecay(charge,energy,rmax,iboundary,yn,ynp)
	ratio=ynp/yn
	do i=2,kp 
	du_pot(i,nb-2)=du_pot(i,nb-2)+du_pot(i-1,nb-1)*(kp-1d0)/(kp-1d0-ratio*(rmax-tx(nb)))
	enddo
	endif

	v=c0
	do j=1,nb-2
	k=kp-j
		do i=max(1,j-kp+1),min(nb-2,j+kp-1)
		v(i,j)=du_pot(k+i,j)
		enddo
	enddo
	deallocate(du_pot)

	return
	end subroutine 
!--------------------------------------------------------------------------------------------	
	subroutine norm(rmax,charge,energy,l,kp,nb,tx,nquad_bsp,x,nm)
!*****************************************************************************
!	Normalization of the coefficients such that <Sn|Sn>=1
!       Output is the normalization vector nm
!*****************************************************************************
	use complex_module
	implicit none
	integer l,kp,nb,nquad_bsp
	real*8 tx(l+2*kp-1),charge,rmax
	complex*16 x(1:nb,1:nb-2),nm(1:nb-2),energy,sn,psi,kr
	complex*16,allocatable :: du_over(:,:),temp(:),y(:)
	integer j,i
	real*8,allocatable :: bb1(:)

	if(real(energy).lt.0d0)then
	allocate(du_over(2*kp-1,nb-1))
	call integ_over(l,kp,nb,tx,nquad_bsp,du_over)

	do j=1,nb-2
	allocate(temp(1:nb-1),y(1:nb-1))
	temp(1:nb-1)=x(2:nb,j)
	call zgbmv('N',nb-1,nb-1,kp-1,kp-1,c1,du_over,2*kp-1,temp,1,c0,y,1)
	nm(j)=dot_product(conjg(temp),y)
	nm(j)=1d0/sqrt(nm(j))
	deallocate(temp,y)
	enddo
	deallocate(du_over)
	endif

	if(real(energy).gt.0d0)then
	allocate(bb1(nb))
	call evaluate_bsplines1(l,nb,kp,rmax,tx,bb1)
	do j=1,nb-2
	sn=c0
		do i=1,nb
		sn=sn+x(i,j)*bb1(i)
		enddo
	kr=sqrt(2d0*energy)
	psi=ci*(kr*rmax+charge*log(2d0*kr*rmax)/kr) 
	psi=exp(psi)
	nm(j)=sn/psi
	enddo
	deallocate(bb1)
	endif

	return
	end subroutine 

!---------------------------------------------------------------------------------------
    	subroutine sort(n,y,A)
!*************************************************************************
!	Sorts vector y in increasing absolute value and the corresponding
!	vector A(:,j) the same way. Output in the same y and A.
!*************************************************************************
  	use complex_module
	implicit none
	integer i,n
	integer,allocatable :: indx(:)
	complex*16 y(1:n),A(1:n,1:n)
	real*8,allocatable :: rytemp(:)
	complex*16,allocatable :: Atemp(:,:),cytemp(:)

	allocate(rytemp(1:n))
	do i=1,n
	rytemp(i)=abs(y(i))
	enddo

	allocate(indx(1:n))
	call indexx(n,rytemp,indx)
	deallocate(rytemp)
	allocate(Atemp(1:n,1:n),cytemp(1:n))
	Atemp=A
	do i=1,n
            cytemp(i)=y(indx(i))
            A(:,i)=Atemp(:,indx(i))
	enddo
	y=cytemp
 	deallocate(indx) 
	deallocate(Atemp,cytemp)

	if(allocated(indx)) deallocate(indx)
	if(allocated(Atemp)) deallocate(Atemp)
	if(allocated(rytemp)) deallocate(rytemp)
	if(allocated(cytemp)) deallocate(cytemp)

	return
   	end subroutine 
!-------------------------------------------------------------------------------------------------
     	subroutine integ_kin(l,kp,nb,tx1,nquad_bsp,du)
!*****************************************************************************************
!	Kinetic energy -0.5*<B(jp,kp)|d²/dx²|B(j,kp)>=0.5*<d/dx(B(jp,kp))|d/dx(B(j,kp))>
!*****************************************************************************************
	use complex_module
	use matrix_module
	implicit none
    	integer i,j,iq,s,l,kp,nb,nquad_bsp
	real*8 a0,tx1(l+2*kp-1)
	complex*16 du(2*kp-1,nb-1)

	!Diagonal elements
	du=c0
	do j=2,nb
		do s=0,kp-1
		if(tx1(j+s).ne.tx1(j+s+1))then
		a0=(tx1(j+s+1)-tx1(j+s))/2d0
			do iq=1,nquad_bsp
			du(kp,j-1)=du(kp,j-1)+wd(iq)*dbb(j,(j+s-kp)*nquad_bsp+iq)*dbb(j,(j+s-kp)*nquad_bsp+iq)*a0
			enddo
		endif
		enddo
	enddo
	
	!Upper diagonal elements
	do i=1,kp-1
		do j=i+2,nb
			do s=0,kp-i-1
			if(tx1(j+s).ne.tx1(j+s+1))then
			a0=(tx1(j+s+1)-tx1(j+s))/2d0
				do iq=1,nquad_bsp
				du(kp-i,j-1)=du(kp-i,j-1)+wd(iq)*dbb(j-i,(j+s-kp)*nquad_bsp+iq)*dbb(j,(j+s-kp)*nquad_bsp+iq)*a0
				enddo
			endif
			enddo
		du(kp+i,j-i-1)=du(kp-i,j-1)
		enddo
	enddo

	du=real(du*0.5d0)+ci*0d0

   	return 
	end
!--------------------------------------------------------------------------------------------------------
    	subroutine integ_pot(l,kp,nb,tx1,rmax,ishort,alpha,charge,u0,rc0,delta0,u1,rc1,delta1,u2,rc2,delta2,nquad_bsp,du)
!*****************************************************************************************
!	Short range potential <B(jp,kp)|-v(r)|B(j,kp)>
!*****************************************************************************************
	use complex_module
	use matrix_module
	implicit none
    	integer i,j,iq,s,l,kp,nb,ishort,nquad_bsp
	real*8 a0,vr,tx1(l+2*kp-1),rmax,alpha,charge
	complex*16 du(2*kp-1,nb-1)
	real*8,external :: potential
	real*8 u0,rc0,delta0,u1,rc1,delta1,u2,rc2,delta2

	!Diagonal elements
	du=c0
	do j=2,nb
		do s=0,kp-1
		if(tx1(j+s).ne.tx1(j+s+1))then
		a0=(tx1(j+s+1)-tx1(j+s))/2d0
			do iq=1,nquad_bsp
			vr=potential(xd((j+s-kp)*nquad_bsp+iq),ishort,charge,alpha,u0,rc0,delta0,rmax)
			du(kp,j-1)=du(kp,j-1)+wd(iq)*bb(j,(j+s-kp)*nquad_bsp+iq)*bb(j,(j+s-kp)*nquad_bsp+iq)*vr*a0
			enddo
		endif
		enddo
	enddo
	
	!Upper diagonal elements
	do i=1,kp-1
		do j=i+2,nb
			do s=0,kp-i-1
			if(tx1(j+s).ne.tx1(j+s+1))then
			a0=(tx1(j+s+1)-tx1(j+s))/2d0
				do iq=1,nquad_bsp
			        vr=potential(xd((j+s-kp)*nquad_bsp+iq),ishort,charge,alpha,u0,rc0,delta0,rmax)
				du(kp-i,j-1)=du(kp-i,j-1)+wd(iq)*bb(j-i,(j+s-kp)*nquad_bsp+iq)*bb(j,(j+s-kp)*nquad_bsp+iq)*vr*a0
				enddo
			endif
			enddo
		du(kp+i,j-i-1)=du(kp-i,j-1)
		enddo
	enddo

	du=(-1d0)*du
	
    	return 
    	end
!------------------------------------------------------------------------------------------------------------------------------
   	subroutine integ_over(l,kp,nb,tx1,nquad_bsp,du)
!*****************************************************************************************
!	Overlap <B(jp,kp)|B(j,kp)>
!*****************************************************************************************
	use matrix_module
	use complex_module
	implicit none
    	integer i,j,iq,s,l,kp,nb,nquad_bsp
	real*8 a0,tx1(l+2*kp-1)
	complex*16 du(2*kp-1,nb-1)

	!Diagonal elements
	du=c0
	do j=2,nb
		do s=0,kp-1
		if(tx1(j+s).ne.tx1(j+s+1))then
		a0=(tx1(j+s+1)-tx1(j+s))/2d0
			do iq=1,nquad_bsp
			du(kp,j-1)=du(kp,j-1)+wd(iq)*bb(j,(j+s-kp)*nquad_bsp+iq)*bb(j,(j+s-kp)*nquad_bsp+iq)*a0
			enddo
		endif
		enddo
	enddo
	
	!Upper diagonal elements
	do i=1,kp-1
		do j=i+2,nb
			do s=0,kp-i-1
			if(tx1(j+s).ne.tx1(j+s+1))then
			a0=(tx1(j+s+1)-tx1(j+s))/2d0
				do iq=1,nquad_bsp
				du(kp-i,j-1)=du(kp-i,j-1)+wd(iq)*bb(j-i,(j+s-kp)*nquad_bsp+iq)*bb(j,(j+s-kp)*nquad_bsp+iq)*a0
				enddo
			endif
			enddo
		du(kp+i,j-i-1)=du(kp-i,j-1)
		enddo
	enddo


	du=real(du)+ci*0d0

   	return 
   	end
!------------------------------------------------------------------------------------------------------------------------------
   	subroutine integ_cou(ilong,u0,rc0,delta0,u1,rc1,delta1,u2,rc2,delta2,charge,rmax,l,kp,nb,tx1,nquad_bsp,du)
!*****************************************************************************************
!	Coulomb <B(jp,kp)|-1/r|B(j,kp)>
!*****************************************************************************************
	use matrix_module
	use complex_module
	implicit none
    	integer i,j,iq,s,l,kp,nb,nquad_bsp,ilong
	real*8 a0,tx1(l+2*kp-1),charge,u0,rc0,delta0,u1,rc1,delta1,u2,rc2,delta2,vr,rmax
	complex*16 du(2*kp-1,nb-1)
	real*8,external :: potential

	!Diagonal elements
	du=c0
	do j=2,nb
		do s=0,kp-1
		if(tx1(j+s).ne.tx1(j+s+1))then
		a0=(tx1(j+s+1)-tx1(j+s))/2d0
			do iq=1,nquad_bsp
			vr=potential(xd((j+s-kp)*nquad_bsp+iq),ilong,charge,1d0,u0,rc0,delta0,rmax)
			du(kp,j-1)=du(kp,j-1)+wd(iq)*bb(j,(j+s-kp)*nquad_bsp+iq)*bb(j,(j+s-kp)*nquad_bsp+iq)*vr*a0
			enddo
		endif
		enddo
	enddo
	
	!Upper diagonal elements
	do i=1,kp-1
		do j=i+2,nb
			do s=0,kp-i-1
			if(tx1(j+s).ne.tx1(j+s+1))then
			a0=(tx1(j+s+1)-tx1(j+s))/2d0
				do iq=1,nquad_bsp
			        vr=potential(xd((j+s-kp)*nquad_bsp+iq),ilong,charge,1d0,u0,rc0,delta0,rmax)
				du(kp-i,j-1)=du(kp-i,j-1)+wd(iq)*bb(j-i,(j+s-kp)*nquad_bsp+iq)*bb(j,(j+s-kp)*nquad_bsp+iq)*vr*a0
				enddo
			endif
			enddo
		du(kp+i,j-i-1)=du(kp-i,j-1)
		enddo
	enddo

	du=real(du)+ci*0d0

   	return 
   	end
!------------------------------------------------------------------------------------------------------------------------------
   	subroutine integ_l(ele,l,kp,nb,tx1,nquad_bsp,du)
!*****************************************************************************************
!	Centrifugal potential <B(jp,kp)|ele*(ele+1)/r²|B(j,kp)>
!*****************************************************************************************
	use complex_module
	use matrix_module
	implicit none
    	integer i,j,iq,s,ele,l,kp,nb,nquad_bsp
	real*8 a0,tx1(l+2*kp-1)
	complex*16 du(2*kp-1,nb-1)

	!Diagonal elements
	du=c0
	do j=2,nb
		do s=0,kp-1
		if(tx1(j+s).ne.tx1(j+s+1))then
		a0=(tx1(j+s+1)-tx1(j+s))/2d0
			do iq=1,nquad_bsp
			du(kp,j-1)=du(kp,j-1)+wd(iq)*bb(j,(j+s-kp)*nquad_bsp+iq)*bb(j,(j+s-kp)*nquad_bsp+iq)*&
					a0/xd((j+s-kp)*nquad_bsp+iq)**2
			enddo
		endif
		enddo
	enddo
	
	!Upper diagonal elements
	do i=1,kp-1
		do j=i+2,nb
			do s=0,kp-i-1
			if(tx1(j+s).ne.tx1(j+s+1))then
			a0=(tx1(j+s+1)-tx1(j+s))/2d0
				do iq=1,nquad_bsp
				du(kp-i,j-1)=du(kp-i,j-1)+wd(iq)*bb(j-i,(j+s-kp)*nquad_bsp+iq)*bb(j,(j+s-kp)*nquad_bsp+iq)*&
						a0/xd((j+s-kp)*nquad_bsp+iq)**2
				enddo
			endif
			enddo
		du(kp+i,j-i-1)=du(kp-i,j-1)
		enddo
	enddo

	du=0.5d0*ele*(ele+1d0)*real(du)+ci*0d0

   	return 
   	end
