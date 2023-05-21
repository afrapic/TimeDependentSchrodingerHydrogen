	subroutine read_input
!*************************************************************************
!	This subroutine reads the input file passed in execution called
!	sturmian-#.input, where # is the number associated later to the
!	outputs.
!	subroutine by A L Frapiccini
!*************************************************************************
    	use datainput_module
	use file_io_module
	use pi_module
    	implicit none
	integer iflag

    	namelist /orderlist/norder
    	namelist /energylist/energy
    	namelist /potlist/ilong,ishort,alpha,charge,u0,rc0,delta0,u1,rc1,delta1,iboundary
    	namelist /lqlist/lq
    	namelist /bsplinelist/kp,nb,rmax,nquad_bsp
	namelist /quadlist/nquad
	namelist /pulselist/peaki,omega,ncycle,omega1,ncycle1,tsep,phase,env,gauge
	namelist /propaglist/qnumber,h,nkry,num_save_vec,tolrndm
	namelist /fatunlalist/acc,accden,maxomh,acctrunc,minh,inch
	namelist /plotlist/plot_extract,plot_vq,plot_rdens,plot_wavepkt,plot_pop,plot_spect_td,plot_rm,ngrids,ngridsk,emax

	iflag=0
    	open(unit=1000,file=file_input)
	read(1000,orderlist)
	read(1000,energylist)
	read(1000,potlist)
	read(1000,lqlist)
	read(1000,bsplinelist)
	read(1000,quadlist)
	read(1000,pulselist)
	read(1000,propaglist) 
	!read(1000,fatunlalist)
	read(1000,plotlist) 
	close(1000)

	if(gauge.eq.1)then
    	amax=5.3378d-09*dsqrt(peaki)
	endif
	if(gauge.eq.2)then
    	amax=5.3378d-09*dsqrt(peaki)/omega
	endif
	tau=2d0*pi*ncycle/omega
	if(env.eq.1)then
	tini_pulse=0
	tfin_pulse=tau
	endif
	if(env.eq.2)then
	tau=2*pi*ncycle/omega
	tau1=2*pi*ncycle1/omega1
	tini1=0d0
	tfin1=tau
	tini2=tfin1+tsep
	tfin2=tini2+tau1
	tini_pulse=tini1
	tfin_pulse=tfin2
	endif
	if(env.eq.3)then
	tau=2*pi*ncycle/omega
	tau1=2*pi*ncycle1/omega
	tini1=0d0
	tini_pulse=tini1
	tfin_pulse=tau1+2*tau
	endif

	if(iboundary.ne.0)then
	iflag=-1
	write(6,*) 'For TDSE iboundary(1)=0'
	endif

	if(iflag.ne.0)then
	write(6,*) 'ERROR IN INPUT FILE'
	stop
	endif

	call name_file(file_input,file_output)

	write(6,*) 'Reading input file from: ',file_input
	write(6,*)

	write(6,*) 'DATA FOR TDSE BASIS'
    	write(6,*) 'Size of the basis:',norder
    	write(6,*) 'Energy:',energy
    	write(6,*) 'Maximum angular momentum:',lq
	if(ilong.eq.0)then
	write(6,*) 'Long range Coulomb',charge
	endif
	if(ilong.eq.5)then
	write(6,*) 'Long range Fullerene+Coulomb',charge
	write(6,*) 'u0:',u0
	write(6,*) 'rc0:',rc0
	write(6,*) 'delta0:',delta0
	endif
	if(ilong.eq.6)then
	write(6,*) 'Long range Two Fullerene+Coulomb',charge
	write(6,*) 'u0:',u0
	write(6,*) 'rc0:',rc0
	write(6,*) 'delta0:',delta0
	write(6,*) 'u1:',u1
	write(6,*) 'rc1:',rc1
	write(6,*) 'delta1:',delta1
	endif
	if(ishort.eq.0)then
	write(6,*) 'Short range Coulomb box'
	endif
	if(ishort.eq.1)then
	write(6,*) 'Short range Yukawa',alpha
	endif
	if(ishort.eq.2)then
	write(6,*) 'Short range Hulten',alpha
 	endif
	if(ishort.eq.3)then
	write(6,*) 'Short range exponential',alpha
 	endif
	if(ishort.eq.4)then
	write(6,*) 'Short range box'
 	endif
	
    	if(real(energy).gt.0d0)then
		if(iboundary.eq.-1)then
            write(6,*) 'Incoming wave boundary condition'
		end if
		if(iboundary.eq.1)then
            write(6,*) 'Outgoing wave boundary condition'
		end if
		if(iboundary.eq.0)then
            write(6,*) 'Stationary wave boundary condition'
        	end if
    	endif
    	if(real(energy).lt.0d0)then
        	if(abs(iboundary).eq.1)then
            write(6,*) 'Exponential decay boundary condition'
        	end if
        	if(iboundary.eq.0)then
            write(6,*) 'Box boundary condition'
        	end if
    	endif
	write(6,*)
    	write(6,*) 'Bspline order',kp
    	write(6,*) 'Number of Bsplines used',nb
	write(6,*) 'Size of the box:',rmax
   	write(6,*) 'Quadrature points for bsplines:',nquad_bsp
	write(6,*) 

	write(6,*) 'PULSE DATA'
	if(gauge.eq.2)then
	write(6,*) 'Velocity gauge'
	endif
	if(gauge.eq.1)then
	write(6,*) 'Length gauge'
	endif
	if(env.eq.1)then
	write(6,'(a)') 'One frequency Sine square envelope'
	write(6,'(a3,1x,f8.5)') 'a0:',amax
	write(6,'(a4,1x,f10.5)') 'tau:',tau
	write(6,'(a6,1x,f6.3)') 'omega:',omega
	write(6,'(a24,1x,f8.2,1x,a1,1x,f8.2,1x,a1)') 'Interval of the pulse: [',tini_pulse,',',tfin_pulse,']'
	endif
	if(env.eq.2)then
	write(6,'(a)') 'Two frequency Sine square envelope'
	write(6,'(a3,1x,f8.5)') 'a0:',amax
	write(6,'(a4,1x,f10.5)') 'tau:',tau
	write(6,'(a6,1x,f5.2)') 'omega:',omega
	write(6,'(a4,1x,f10.5)') 'tau1:',tau1
	write(6,'(a6,1x,f5.2)') 'omega1:',omega1
	write(6,'(a20,1x,f5.2)') 'Time between pulses:',tsep
	write(6,'(a24,1x,f8.2,1x,a1,1x,f8.2,1x,a1)') 'Interval of the pulse: [',tini_pulse,',',tfin_pulse,']'
	endif
	if(env.eq.3)then
	write(6,'(a)') 'Sine square envelope with a plateau'
	write(6,'(a3,1x,f8.5)') 'a0:',amax
	write(6,'(a4,1x,f10.5)') 'tau:',tau
	write(6,'(a6,1x,f5.2)') 'omega:',omega
	write(6,'(a15,1x,f10.5)') 'tau1 (plateau):',tau1
	write(6,'(a24,1x,f8.2,1x,a1,1x,f8.2,1x,a1)') 'Interval of the pulse: [',tini_pulse,',',tfin_pulse,']'
	endif
	write(6,*)

	write(6,*) 'TDSE data'
	write(6,*) 'Initial state n:',qnumber(1),'l:',qnumber(2)
	write(6,*) 'Stepsize for Arnoldi:',h
	write(6,*) 'Number of Krylov vectors:',nkry

	write(6,*)

	ntot1=(lq+1)*norder
	ntot2=lq*norder 

	return
	end
