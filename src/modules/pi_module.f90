      module pi_module
!**************************************
!  factors involving pi.
!**************************************
         implicit none
         save
      real*8, parameter :: fourthpi =  0.785398163397448d0
      real*8, parameter :: halfpi   = 2.d0*fourthpi
      real*8, parameter :: pi       = 4.d0*fourthpi
      real*8, parameter :: twopi    = 8.d0*fourthpi
      end module pi_module
