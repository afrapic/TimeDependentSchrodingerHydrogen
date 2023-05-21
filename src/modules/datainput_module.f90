	module datainput_module
!*******************************************************
!	Stores the data in the input file and 
!	derivated variables
!	module by A L Frapiccini
!*******************************************************
	implicit none
!sturmians data for time propagation
	integer scatt,inp
	integer norder,ntot1,ntot2,ntot3
	complex*16 energy
	integer ilong
	integer ishort
	real*8 alpha
	real*8 charge 
	real*8 u0,u1,u2
	real*8 rc0,rc1,rc2
	real*8 delta0,delta1,delta2
	integer iboundary
	integer lq
	integer kp
	integer nb
	real*8 rmax
	integer nquad_bsp
!Data to integrate
	integer nquad
!Data for the pulse
	real*8 omega,omega1
	real*8 peaki
	real*8 ncycle,ncycle1
	real*8 phase,tsep
	integer env,gauge
	real*8 period,tini_pulse,tfin_pulse,tau,amax,tfin,tini1,tfin1,tini2,tfin2,tau1
!Data for the propagation
	real*8 h,tolrndm
	integer nkry
	integer num_save_vec
	integer qnumber(2)
	real*8 acc,accden,maxomh,acctrunc,minh,inch
!Data for plotting
	integer plot_rdens,plot_wavepkt,plot_pop,plot_extract,plot_spect_td,plot_rm,ngrids,ngridsk,plot_vq
	real*8 emax

	end module
