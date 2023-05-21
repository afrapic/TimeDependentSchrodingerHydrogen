	subroutine propag_fatunla(h_atomic,m_over,m_chol,h_int_up,yt)
	use datainput_module
	use file_io_module
	use scratchdir_module
	use complex_module
    	implicit none
	real*8 h_atomic(ntot1,norder)
	real*8 m_over(ntot1,norder)
	real*8 m_chol(ntot1,norder)
	real*8 h_int_up(ntot2,norder)
    	integer it,pulse,isave,ierr
	real*8  errf,nrm,h0
	real*8 dt_save,tsave
	real*8 t,t_new,tfinal,timei,timef,time1
	complex*16 yt(ntot1),y_save(ntot1),error(ntot1),prod
	real*8,external :: dsecnd
	integer,parameter :: maxstep=2000000000

	open(unit=45,file=trim(dir_scratch)//'tsave-'//trim(file_output)//'.dat')

	!Initialize variables
	t=tini_pulse
	tfinal=tfin_pulse
	tfin=tfin_pulse
	call norm2(m_over,yt,nrm)
	yt=yt/nrm
	y_save=yt
	h0=h
	dt_save=(tfin-tini_pulse)/dble(num_save_vec)

	!prints initial vector
	call write_coef(0,tini_pulse,yt)

	!write(6,*) 'Pulse propagation'
	isave=1
	tsave=tini_pulse+isave*dt_save

	timei=dsecnd()
	do it=1,maxstep
	y_save=yt

	!Fatula explicit scheme
	call fatunla(h0,t,h_atomic,m_over,m_chol,h_int_up,yt,error)
	prod=dot_product(error,error)
	errf=real(prod)

	!If errf<acctrunc continues
	if(errf.lt.acctrunc)then
	t_new=t+h0
	y_save=yt
	write(6,*) t_new,h0,errf

		!Exit if tfinal reached
		if(t_new.eq.tfinal)then
		timef=dsecnd()
		call write_coef(-1,tfinal,yt)
		call write_coef(num_save_vec,tfinal,yt)
		exit
		endif

		!Saves vector if reached t_save
		if(isave.ne.num_save_vec)then
		if(t_new.ge.tsave)then
		call write_coef(isave,t_new,yt)
		isave=isave+1
		tsave=tini_pulse+isave*dt_save
		endif
		endif

	
		!Adjust time step at the end of the pulse and at the end of free propagation
     		if(t_new.gt.tfinal) then
		h0=tfinal-t
		t_new=t+h0
        	t=t_new
		yt=y_save
       	else
			!!If errf<1d-15 increases time step
      		if(errf.lt.acctrunc*1d-5) then
         		h0=h0*inch
         		t=t_new
      		else
         		t=t_new
         		endif
      	endif
   	else
	!If errf>1d-10 rejects time step, decreases and goes back.
   	yt=y_save
	h0=h0/inch
		if(h0.lt.minh) then
      	write(6,*) "program stop because step become to small" 		
		stop
      	endif
 	endif
		
	enddo

	if(it.eq.maxstep)then
		if(t.ne.tfinal)then
		write(6,*) 'maxstep number of iterations in propagator and tfinal not reached. Increase maxstep in input'
		stop
		endif
	endif

	call print_wtime(timef-timei)

    return
    end
!*******************************************************************************************************************************
	subroutine fatunla(h0,t,h_atomic,m_over,m_chol,h_int_up,vec,error)
	use datainput_module
	use file_io_module
	use scratchdir_module
	use complex_module
	implicit none
	integer i,j,pulse
	real*8 h_atomic(ntot1,norder)
	real*8 m_over(ntot1,norder)
	real*8 m_chol(ntot1,norder)
	real*8 h_int_up(ntot2,norder)
	real*8 h0,t
	complex*16 vec(ntot1),R(ntot1),S(ntot1),w1(ntot1),w2(ntot1)
	complex*16 fd(-1:4,ntot1),error(ntot1)

	!Fatounla scheme execution
	fd(-1,:)=vec(:)
	call deriv_comput(t,fd,h_atomic,m_chol,h_int_up)
	call stiff_param(h0,fd,w1,w2,R,S)

	!Computes the new vector 
	do j=1,ntot1
		vec(j)=vec(j)+R(j)*fd(0,j)+S(j)*fd(1,j)
	enddo

	!Computes truncation error
	error=c0
	do j=1,ntot1
		error(j)=h0**5*(fd(4,j)+(w2(j)**3-w1(j)*w2(j)**2+w2(j)*w1(j)**2-w1(j)**3)*fd(1,j)&
			   -w1(j)*w2(j)*(w1(j)**2-w1(j)*w2(j)+w2(j)**2)*fd(0,j))/120.0d+0
	enddo

	return
	end
!**********************************************************************************************************************************
	subroutine deriv_comput(t,fd,h_atomic,m_chol,h_int_up)          
	use datainput_module
	use file_io_module
	use scratchdir_module
	use complex_module
	implicit none
	integer i,pulse
	real*8 t
	real*8 h_atomic(ntot1,norder)
	real*8 m_chol(ntot1,norder)
	real*8 h_int_up(ntot2,norder)
 	complex*16 fd(-1:4,ntot1),gt(0:4)
	complex*16 sumf(ntot1),fdn(ntot1)

	!Computes different derivatives of pulse
	call pulse_deriv(t,gt)

	!Evaluates the zero order derivative
	call deriv(t,0,fd,gt,h_atomic,m_chol,h_int_up)
	!Multiplies by inverse 
	fdn(:)=fd(0,:)
	sumf=c0
	call minverse(m_chol,fdn,sumf)
	fd(0,:)=sumf(:)
                
	!Evaluates the first order derivative
	call deriv(t,1,fd,gt,h_atomic,m_chol,h_int_up)
	!Multiplies by inverse 
	fdn(:)=fd(1,:)
	sumf=c0
	call minverse(m_chol,fdn,sumf)
	fd(1,:)=sumf(:)
                 
	!Evaluates the second order derivative
	call deriv(t,2,fd,gt,h_atomic,m_chol,h_int_up)
	!Multiplies by inverse 
	fdn(:)=fd(2,:)
	sumf=c0
	call minverse(m_chol,fdn,sumf)
	fd(2,:)=sumf(:)
        
	!Evaluates the third order derivative
	call deriv(t,3,fd,gt,h_atomic,m_chol,h_int_up)
	!Multiplies by inverse 
	fdn(:)=fd(3,:)
	sumf=c0
	call minverse(m_chol,fdn,sumf)
	fd(3,:)=sumf(:)
              
	!Evaluates the fourth order derivative
	call deriv(t,4,fd,gt,h_atomic,m_chol,h_int_up)
	!Multiplies by inverse 
	fdn(:)=fd(4,:)
	sumf=c0
	call minverse(m_chol,fdn,sumf)
	fd(4,:)=sumf(:)  


	return
	end
! *********************************************************************** 
	subroutine deriv(t,nderiv,fd,gt,h_atomic,m_chol,h_int_up)                  
	use datainput_module
	use file_io_module
	use scratchdir_module
	use complex_module
	implicit none
	integer i,nderiv
	real*8 h_atomic(ntot1,norder)
	real*8 m_chol(ntot1,norder)
	real*8 h_int_up(ntot2,norder)
	real*8 t            
	complex*16  fd(-1:4,ntot1),gt(0:4)
	complex*16,allocatable :: y_a(:),y_i(:),fdp(:)

	!Computes the corresponding vector for the derivatives
	if(nderiv.eq.0)then
    	fd(0,:)=gt(0)*fd(-1,:)
    	endif    
	if(nderiv.eq.1)then
    	fd(1,:)=gt(0)*fd(0,:)+gt(1)*fd(-1,:)
    	endif    
	if(nderiv.eq.2)then
	fd(2,:)=gt(0)*fd(1,:)+2d0*gt(1)*fd(0,:)+gt(2)*fd(-1,:)
	endif
	if(nderiv.eq.3)then
	fd(3,:)=gt(0)*fd(2,:)+3d0*gt(1)*fd(1,:)+3d0*gt(2)*fd(0,:)+gt(3)*fd(-1,:)
	endif
	if(nderiv.eq.4)then
	fd(4,:)=gt(0)*fd(3,:)+4d0*gt(1)*fd(2,:)+6d0*gt(2)*fd(1,:)+4d0*gt(3)*fd(0,:)+gt(4)*fd(-1,:)
	endif
	
	allocate(y_a(ntot1),y_i(ntot1),fdp(ntot1))
	y_a=c0
	y_i=c0
	fdp(:)=fd(nderiv-1,:)
	call hatomic_vec(h_atomic,fdp,y_a)
	fdp=c0   
	fdp(:)=fd(nderiv,:)
	call hint_vec(c1,h_int_up,fdp,y_i)
	deallocate(fdp)

	fd(nderiv,:)=-ci*(y_a(:)+y_i(:))
	deallocate(y_a,y_i)                  

	return
	end
!****************************************************************
	subroutine pulse_deriv(t,gt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!	The time dependency of the interaction with the pulse is written as
!!!	H(int)=g(t)*mat_int*constants
!!!	with g(t)=	A(t) for velocity gauge
!!!			E(t)=-dA(t)/dt for length gauge
!!!	Then constants are added in the multiplication routines.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use datainput_module
	use pi_module
	use complex_module
	implicit none   
	real*8 t,t1,t2
	complex*16 gt(0:4),cte

	gt=c0
	if(t.gt.tfin_pulse)then
	gt=c0
	return
	endif

	if(gauge.eq.1)then
	cte=cmplx(amax)
	endif
	if(gauge.eq.2)then
	cte=-ci*amax
	endif

	if(env.eq.1)then
      	gt(0)=sin(omega*t)*sin(pi*t/tau)**2
		gt(1)=sin(pi*t/tau)*((2d0*pi*Cos((pi*t)/tau)*Sin(omega*t))/tau+omega*Cos(omega*t)*Sin((pi*t)/tau))
		gt(2)=((-(omega**2*tau**2)+(4d0*pi**2+omega**2*tau**2)*Cos((2d0*pi*t)/tau))*Sin(omega*t)+4d0*omega*pi*tau*&
			Cos(omega*t)*Sin((2d0*pi*t)/tau))/(2d0*tau**2)
		gt(3)=((omega*tau*Cos(omega*t)*(-(omega**2*tau**2)+(12d0*pi**2+omega**2*tau**2)*Cos((2d0*pi*t)/tau)))/2d0-&
			pi*(4d0*pi**2+3d0*omega**2*tau**2)*Sin(omega*t)*Sin((2d0*pi*t)/tau))/tau**3
		gt(4)=-((-(omega**4*tau**4)+(16d0*pi**4+24d0*omega**2*pi**2*tau**2+omega**4*tau**4)*Cos((2d0*pi*t)/tau))&
			*Sin(omega*t)+8d0*omega*pi*tau*(4d0*pi**2+omega**2*tau**2)*Cos(omega*t)*Sin((2d0*pi*t)/tau))/&
     			(2d0*tau**4)
 	endif


     	return
     	end
!---------------------------------------------------------------------------------------------------------------------------
   	subroutine stiff_param(h0,fd,w1,w2,R,S)
	use datainput_module
	use file_io_module
	use scratchdir_module
	use complex_module
    	implicit none         
    	integer i,j     
    	real*8 h0
    	complex*16 w1(ntot1),w2(ntot1),R(ntot1),S(ntot1),d,e,denom,fd(-1:4,ntot1)
    	complex*16 tetha(ntot1), lambda(ntot1)
   	logical crit(ntot1)

	!Computes the stiffness parameters *******    
	lambda=c0
	tetha=c0
	crit=.true.

	do j=1,ntot1
		denom=fd(1,j)**2-fd(0,j)*fd(2,j)
			if(abs(denom).lt.accden) then
			R(j)=c0
			S(j)=c0
			w1(j)=c0
			w2(j)=c0
			crit(j)=.false. 
			d=c0
			e=c0  
			else
			d=(fd(0,j)*fd(3,j)-fd(1,j)*fd(2,j))/denom
			e=(fd(1,j)*fd(3,j)-fd(2,j)**2)/denom
			w1(j)=0.5d0*(-d+cdsqrt(d**2+4d0*e))
			w2(j)=w1(j)+d                            
       		endif
	enddo

	do j=1,ntot1
		if(crit(j)) then
			if(dble(w1(j)*h0).gt.maxomh.or.dble(-w2(j)*h0).gt.maxomh.or.(abs(fd(0,j)).lt.acc.and.abs(fd(1,j)).lt.acc))then            
			R(j)=c0
			S(j)=c0
			crit(j)=.false.
			endif
		endif
	enddo

	do j=1,ntot1
		if(crit(j)) then
			if(abs(w1(j)).lt.acctrunc) then
			tetha(j)=h0/w2(j)
			else
			tetha(j)=(cdexp(w1(j)*h0)-1d0)/(w1(j)*(w1(j)+w2(j)))
			endif
			if(abs(w2(j)).lt.acctrunc) then
			lambda(j)=-h0/w1(j)
			else
			lambda(j)=(cdexp(-w2(j)*h0)-1d0)/(w2(j)*(w1(j)+w2(j)))
			endif            
		R(j)=w2(j)*tetha(j)-w1(j)*lambda(j)
		S(j)=tetha(j)+lambda(j)           
		endif
	enddo

	return
	end











