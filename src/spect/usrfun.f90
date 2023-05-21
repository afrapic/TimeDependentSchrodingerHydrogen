	subroutine usrfun(en,l,x,n,NP,fvec,fjac)
	use datainput_module
	implicit none
	real*8 x(n),fvec(NP),fjac(NP,NP),en
	integer n,NP,l
	real*8 fa1,fb1,fb2,fb11,fb21,fc1,fc2,g,gp,err
	real*8 dfa1,dfb1,dfb2,dfb11,dfb21,dfc1,dfc2
	real*8 fc11,fc21,fd1,fd2,fd11,fd21,fe1,fe2
	real*8 dfc11,dfc21,dfd1,dfd2,dfd11,dfd21,dfe1,dfe2

	if(ilong.eq.5)then

	call scoul(-charge,en,l,rc0,fa1,dfa1,g,gp,err)
	call scoul(-charge,en+U0,l,rc0,fb1,dfb1,fb2,dfb2,err)
	call scoul(-charge,en+u0,l,rc0+delta0,fb11,dfb11,fb21,dfb21,err)
	call scoul(-charge,en,l,rc0+delta0,fc1,dfc1,fc2,dfc2,err)

	fvec(1)=x(1)*fa1-x(2)*fb1-x(3)*fb2
	fvec(2)=x(1)*dfa1-x(2)*dfb1-x(3)*dfb2 
	fvec(3)=x(2)*fb11+x(3)*fb21-x(4)*fc1-x(5)*fc2 
	fvec(4)=x(2)*dfb11+x(3)*dfb21-x(4)*dfc1-x(5)*dfc2 
	fvec(5)=x(4)**2+x(5)**2-1d0

	fjac=0d0
	fjac(1,1)=fa1
	fjac(1,2)=-fb1
	fjac(1,3)=-fb2
	fjac(2,1)=dfa1
	fjac(2,2)=-dfb1
	fjac(2,3)=-dfb2
	fjac(3,2)=fb11
	fjac(3,3)=fb21
	fjac(3,4)=-fc1
	fjac(3,5)=-fc2
	fjac(4,2)=dfb11
	fjac(4,3)=dfb21
	fjac(4,4)=-dfc1
	fjac(4,5)=-dfc2
	fjac(5,4)=2d0*x(4)
	fjac(5,5)=2d0*x(5)
	endif

	if(ilong.eq.6)then

	call scoul(-charge,en,l,rc0,fa1,dfa1,g,gp,err)
	call scoul(-charge,en+u0,l,rc0,fb1,dfb1,fb2,dfb2,err)
	call scoul(-charge,en+u0,l,rc0+delta0,fb11,dfb11,fb21,dfb21,err)
	call scoul(-charge,en,l,rc0+delta0,fc1,dfc1,fc2,dfc2,err)
	call scoul(-charge,en,l,rc1,fc11,dfc11,fc21,dfc21,err)
	call scoul(-charge,en+u1,l,rc1,fd1,dfd1,fd2,dfd2,err)
	call scoul(-charge,en+u1,l,rc1+delta1,fd11,dfd11,fd21,dfd21,err)
	call scoul(-charge,en,l,rc1+delta1,fe1,dfe1,fe2,dfe2,err)

	fvec(1)=x(1)*fa1-x(2)*fb1-x(3)*fb2
	fvec(2)=x(1)*dfa1-x(2)*dfb1-x(3)*dfb2 
	fvec(3)=x(2)*fb11+x(3)*fb21-x(4)*fc1-x(5)*fc2 
	fvec(4)=x(2)*dfb11+x(3)*dfb21-x(4)*dfc1-x(5)*dfc2 
	fvec(5)=x(4)*fc11+x(5)*fc21-x(6)*fd1-x(7)*fd2 
	fvec(6)=x(4)*dfc11+x(5)*dfc21-x(6)*dfd1-x(7)*dfd2
	fvec(7)=x(6)*fd11+x(7)*fd21-x(8)*fe1-x(9)*fe2 
	fvec(8)=x(6)*dfd11+x(7)*dfd21-x(8)*dfe1-x(9)*dfe2	
	fvec(9)=x(8)**2+x(9)**2-1d0

	fjac=0d0
	fjac(1,1)=fa1
	fjac(1,2)=-fb1
	fjac(1,3)=-fb2
	fjac(2,1)=dfa1
	fjac(2,2)=-dfb1
	fjac(2,3)=-dfb2
	fjac(3,2)=fb11
	fjac(3,3)=fb21
	fjac(3,4)=-fc1
	fjac(3,5)=-fc2
	fjac(4,2)=dfb11
	fjac(4,3)=dfb21
	fjac(4,4)=-dfc1
	fjac(4,5)=-dfc2
	fjac(5,4)=fc11
	fjac(5,5)=fc21
	fjac(5,6)=-fd1
	fjac(5,7)=-fd2
	fjac(6,4)=dfc11
	fjac(6,5)=dfc21
	fjac(6,6)=-dfd1
	fjac(6,7)=-dfd2
	fjac(7,6)=fd11
	fjac(7,7)=fd21
	fjac(7,8)=-fe1
	fjac(7,9)=-fe2
	fjac(8,6)=dfd11
	fjac(8,7)=dfd21
	fjac(8,8)=-dfe1
	fjac(8,9)=-dfe2
	fjac(9,8)=2d0*x(8)
	fjac(9,9)=2d0*x(9)

	endif


	return
	end
