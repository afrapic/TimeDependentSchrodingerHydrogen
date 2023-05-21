	module matrix_module
!****************************************************
!	Stores variables related to the Bsplines
!	module by A L Frapiccini 
!****************************************************
	implicit none
	real*8,allocatable :: tx(:)
	real*8,allocatable :: xq(:),wq(:),xs0(:),xs1(:),xs2(:)
	real*8,allocatable :: xd(:),wd(:),bb(:,:),dbb(:,:)
	complex*16,allocatable :: beta(:,:),st(:,:,:),dst(:,:,:),ss0(:,:,:),ss1(:,:,:),ss2(:,:,:),st_bnd(:,:,:)
	complex*16,allocatable :: sa(:,:,:),sb(:,:,:),sc(:,:,:),sd(:,:,:)

	end module
