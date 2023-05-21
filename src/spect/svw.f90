	subroutine svw(en,l,fc,xq2)
	use datainput_module
	implicit none
	integer l,n,ntrial,ncase,i
	real*8 en,xq2(nquad),fc(nquad),tolx,tolf
	real*8 f,fp,g,gp,err
	real*8,allocatable :: x(:)

	tolx=1d-8
	tolf=1d-10
	ntrial=100
	if(ilong.eq.5)then
	n=5
	allocate(x(n))
	x=1d0
	call mnewt(en,l,ntrial,x,n,tolx,tolf)

	do i=1,nquad
		if(xq2(i).le.rc0)then
		call scoul(-charge,en,l,xq2(i),f,fp,g,gp,err)
		fc(i)=x(1)*f
		endif
		if(xq2(i).gt.rc0.and.xq2(i).le.rc0+delta0)then
		call scoul(-charge,en+u0,l,xq2(i),f,fp,g,gp,err)
		fc(i)=x(2)*f+x(3)*g
		endif
		if(xq2(i).gt.rc0+delta0)then
		call scoul(-charge,en,l,xq2(i),f,fp,g,gp,err)
		fc(i)=x(4)*f+x(5)*g
		endif
	enddo

	deallocate(x)
	endif

	if(ilong.eq.6)then
	n=9
	allocate(x(n))
	x=1d0
	call mnewt(en,l,ntrial,x,n,tolx,tolf)

	do i=1,nquad
		if(xq2(i).le.rc0)then
		call scoul(-charge,en,l,xq2(i),f,fp,g,gp,err)
		fc(i)=x(1)*f
		endif
		if(xq2(i).gt.rc0.and.xq2(i).le.rc0+delta0)then
		call scoul(-charge,en+u0,l,xq2(i),f,fp,g,gp,err)
		fc(i)=x(2)*f+x(3)*g
		endif
		if(xq2(i).gt.rc0+delta0.and.xq2(i).le.rc1)then
		call scoul(-charge,en,l,xq2(i),f,fp,g,gp,err)
		fc(i)=x(4)*f+x(5)*g
		endif
		if(xq2(i).gt.rc1.and.xq2(i).le.rc1+delta1)then
		call scoul(-charge,en+u1,l,xq2(i),f,fp,g,gp,err)
		fc(i)=x(6)*f+x(7)*g
		endif
		if(xq2(i).gt.rc1+delta1)then
		call scoul(-charge,en,l,xq2(i),f,fp,g,gp,err)
		fc(i)=x(8)*f+x(9)*g
		endif
	enddo
	deallocate(x)
	endif
	
	return
	end
