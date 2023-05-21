	subroutine print_wtime(time)
!***********************************************************************************
!	Transforms the time given by mpi_wtime to second,minute,hour,day format
!	and prints to screen.
!	subroutine by A L Frapiccini
!***********************************************************************************
	implicit none
	real*8 time
	real*8 days,hours,minutes,seconds

	if(time.lt.60d0)then
	seconds=time
	write(6,'(a5,1x,f6.2,1x,a4)') 'Time:',seconds,' sec.'
	return
	endif

	if(time.ge.60d0)then
		minutes=1d0*int(time/60d0)
		seconds=mod(time,60d0)
		if(minutes.lt.60d0)then
		write(6,'(a5,1x,i2,1x,a4,1x,f6.2,1x,a4)') 'Time:',int(minutes),' min.',seconds,' sec.'
		return
		endif
		if(minutes.ge.60d0)then
		hours=1d0*int(minutes/60d0)
		minutes=mod(minutes,60d0)
			if(hours.lt.24d0)then
			write(6,'(a5,1x,i2,a4,1x,i2,1x,a4,1x,f6.2,1x,a4)') 'Time:',int(hours),' hrs.',int(minutes),' min.',seconds,' sec.'
			return
			endif
			if(hours.ge.24d0)then
			days=1d0*int(hours/24d0)
			hours=mod(hours,24d0)
			write(6,'(a5,1x,i2,a5,1x,i2,a4,1x,i2,1x,a4,1x,f6.2,1x,a4)') 'Time:',int(days),' days',int(hours),&
					' hrs.',int(minutes),' min.',seconds,' sec.'
			return
			endif
		endif
	endif

	return
	end
