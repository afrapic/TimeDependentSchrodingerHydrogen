module file_io_module
!********************************************************************************************
!	Stores subroutines related to input/output of data
!	module by F Colavecchia, additions by A L Frapiccini
!********************************************************************************************
!use datainput_module
implicit none

 	character*20 file_output
 	character*20 file_input
	character*11 progname
	character*7,parameter :: dir_output='output/'	

	!Interface to write formatted output matrix or vector (complex*16,real*8 or integer)
 	interface write_matrix
        module procedure write_complex_matrix,write_real8_matrix, &
        write_complex_vector,write_real8_vector,write_integer_matrix,write_integer_vector
    	end interface

	!Interface to write unformatted output matrix or vector (complex*16,real*8 or integer)
    	interface write_matrix_raw
        module procedure write_complex_matrix_raw,write_real8_matrix_raw,  &
        write_complex_vector_raw,write_real8_vector_raw,write_integer_matrix_raw,write_integer_vector_raw, &
	  write_complex_tensor_raw,write_real8_tensor_raw,write_integer_tensor_raw
    	end interface

	!Interface to read formatted output matrix or vector (complex*16,real*8 or integer)
    	interface read_matrix
        module procedure read_complex_matrix,read_real8_matrix,    &
        read_complex_vector,read_real8_vector,read_integer_matrix,read_integer_vector
    	end interface

	!Interface to read unformatted output matrix or vector (complex*16,real*8 or integer)
    	interface read_matrix_raw
        module procedure read_complex_matrix_raw,read_real8_matrix_raw,  &
        read_complex_vector_raw,read_real8_vector_raw,read_integer_matrix_raw,read_integer_vector_raw, &
        read_complex_tensor_raw,read_real8_tensor_raw,read_integer_tensor_raw
    	end interface

	contains 

	subroutine name_file(file_in,file_out)
	!Gets the number of input from file_in and stores in file_out as a character*2
	character*20 file_in,file_out
	integer ii,ifi,s,lq,i

    	ii=1
	ifi=1
   	 s=len_trim(file_in)
    	do i=1,s
        if(file_in(i:i).eq.'-') ii=i+1
        if(file_in(i:i).eq.'.') then
            ifi=i-1
            exit
        endif
    	enddo
    	read(file_in(ii:ifi),*) lq
	if(floor(dble(lq)/100d0).lt.1d0)then
    	write(file_out,'(i2.2)') lq
	else
	write(file_out,'(i3.3)') lq
	endif
	
	end subroutine

    !WRITE
	!
    !   complex*16, formatted
    !
    	subroutine write_complex_matrix(h,filename)
        implicit none
        integer mh,nh,i,j
        complex*16 :: h(:,:)
        character(len=*) filename
        mh = size(h,1)
        nh = size(h,2)
        open(unit=1000,file=filename)
        do i=1,mh
            do j=1,nh
                write(1000,'(2es18.8)') h(i,j)
            end do
        end do
        return
    	end subroutine write_complex_matrix
    	subroutine write_complex_vector(h,filename)
        implicit none
        integer mh,i
        complex*16 :: h(:)
        character(len=*) filename
        mh = size(h,1)
        open(unit=1000,file=filename)
        do i=1,mh
            write(1000,'(2es18.8)') h(i)
        end do
        return
    	end subroutine write_complex_vector
    !
    !   real*8, formatted
    !
    	subroutine write_real8_matrix(h,filename)
        implicit none
        integer mh,nh,i,j
        real*8 :: h(:,:)
        character(len=*) filename
        mh = size(h,1)
        nh = size(h,2)

        open(unit=1000,file=filename)
        do i=1,mh
            do j=1,nh
                write(1000,'(2es18.8)') h(i,j)
            end do
        end do
        return
    	end subroutine write_real8_matrix
    	subroutine write_real8_vector(h,filename)
        implicit none
        integer mh,i
        real*8 :: h(:)
        character(len=*) filename
        mh = size(h,1)
        open(unit=1000,file=filename)
        do i=1,mh
            write(1000,'(2es18.8)') h(i)
        end do
        return
    	end subroutine write_real8_vector
        subroutine write_integer_matrix(h,filename)
        implicit none
        integer :: h(:,:)
        integer i,j,mh,nh
        character(len=*) filename
        mh=size(h,1)
        nh=size(h,2)
        open(unit=1000,file=filename)
        do i=1,mh
        do j=1,nh
        write(1000,*) h(i,j)
        enddo
        enddo
        end subroutine write_integer_matrix
        subroutine write_integer_vector(h,filename)
        implicit none
        integer :: h(:)
        integer i,mh
        character(len=*) filename
	mh=size(h,1)
        open(unit=1000,file=filename)
	do i=1,mh
        write(1000,*) h(i)
	enddo
        end subroutine write_integer_vector
        
    !WRITE
	!
    !   complex*16, unformatted
    !
    	subroutine write_complex_matrix_raw(h,filename)
        implicit none
        integer mh,nh
        complex*16 :: h(:,:)
        character(len=*) filename
        mh = size(h,1)
        nh = size(h,2)
        open(unit=1000,file=filename,form='unformatted')
        write(1000) h
        return
    	end subroutine write_complex_matrix_raw
    	subroutine write_complex_vector_raw(h,filename)
        implicit none
        integer mh
        complex*16 :: h(:)
        character(len=*) filename
        mh = size(h,1)
        open(unit=1000,file=filename,form='unformatted')
        write(1000) h
        return
    	end subroutine write_complex_vector_raw
    !
    !   real*8 unformatted
    !
   	subroutine write_real8_matrix_raw(h,filename)
        implicit none
        integer mh,nh
        real*8 :: h(:,:)
        character(len=*) filename
        mh = size(h,1)
        nh = size(h,2)

        open(unit=1000,file=filename,form='unformatted')
        write(1000) h
        return
    	end subroutine write_real8_matrix_raw
    	subroutine write_real8_vector_raw(h,filename)
        implicit none
        integer mh
	  real*8 :: h(:)
        character(len=*) filename
        mh = size(h,1)
        open(unit=1000,file=filename,form='unformatted')
        write(1000) h
        return
    	end subroutine write_real8_vector_raw

    !   integer, unformatted
    !
    	subroutine write_integer_matrix_raw(h,filename)
        implicit none
        integer h(:,:)
        character(len=*) filename
        open(unit=1000,file=filename,form='unformatted')
        write(1000) h
        return
    	end subroutine write_integer_matrix_raw
    	subroutine write_integer_vector_raw(h,filename)
        implicit none
        integer mh
        integer :: h(:)
        character(len=*) filename
        mh = size(h,1)
        open(unit=1000,file=filename,form='unformatted')
        write(1000) h
        return
    	end subroutine write_integer_vector_raw

   	subroutine write_complex_tensor_raw(h,filename)
        implicit none
        integer mh,nh,jh
        complex*16 :: h(:,:,:)
        character(len=*) filename
        mh = size(h,1)
        nh = size(h,2)
        jh = size(h,3)
        open(unit=1000,file=filename,form='unformatted')
        write(1000) h
        return
    	end subroutine write_complex_tensor_raw
    !   real*8 unformatted
    !
   	subroutine write_real8_tensor_raw(h,filename)
        implicit none
        integer mh,nh,jh
        real*8 :: h(:,:,:)
        character(len=*) filename
        mh = size(h,1)
        nh = size(h,2)
        jh = size(h,3)

        open(unit=1000,file=filename,form='unformatted')
        write(1000) h
        return
    	end subroutine write_real8_tensor_raw
    !   integer, unformatted
    !
    	subroutine write_integer_tensor_raw(h,filename)
        implicit none
        integer mh,nh,jh
        integer h(:,:,:)
        character(len=*) filename
        open(unit=1000,file=filename,form='unformatted')
        write(1000) h
        return
    	end subroutine write_integer_tensor_raw
 
  !READ
	!
    !   complex*16, formatted
    !
	subroutine read_complex_matrix(h,filename)
        implicit none
        integer mh,nh,i,j
        complex*16 :: h(:,:)
        character(len=*) filename
        mh = size(h,1)
        nh = size(h,2)

        open(unit=1000,file=filename)
        do i=1,mh
            do j=1,nh
               read(1000,'(2es18.8)') h(i,j)
            end do
        end do
        return
    	end subroutine read_complex_matrix
    	subroutine read_complex_vector(h,filename)
        implicit none
        integer mh,i
        complex*16 :: h(:)
        character(len=*) filename
        mh = size(h,1)

        open(unit=1000,file=filename)
        do i=1,mh
            read(1000,'(2es18.8)') h(i)
        end do
        return
    	end subroutine read_complex_vector
    !
    !   real*8, formatted
    !
    	subroutine read_real8_matrix(h,filename)
        implicit none
        integer mh,nh,i,j
        real*8 :: h(:,:)
        character(len=*) filename
        mh = size(h,1)
        nh = size(h,2)

        open(unit=1000,file=filename)
        do i=1,mh
            do j=1,nh
                read(1000,'(2es18.8)') h(i,j)
            end do
        end do
        return
    	end subroutine read_real8_matrix
    	subroutine read_real8_vector(h,filename)
        implicit none
        integer mh,i
        real*8 :: h(:)
        character(len=*) filename
        mh = size(h,1)

        open(unit=1000,file=filename)
        do i=1,mh
            read(1000,'(2es18.8)') h(i)
        end do
        return
    	end subroutine read_real8_vector
        subroutine read_integer_matrix(h,filename)
        implicit none
        integer :: h(:,:)
        integer mh,nh,i,j
        character(len=*) filename
        mh=size(h,1)
        nh=size(h,2)
        open(unit=1000,file=filename)
        do i=1,mh
        do j=1,nh
        read(1000,*) h(i,j)
        enddo
        enddo
        end subroutine read_integer_matrix
        subroutine read_integer_vector(h,filename)
        implicit none
        integer :: h(:)
        integer mh,i
        character(len=*) filename
	mh=size(h,1)
        open(unit=1000,file=filename)
	do i=1,mh
        read(1000,*) h(i)
	enddo
        end subroutine read_integer_vector
  !READ
	!
    !   complex*16, unformatted
    !
   	subroutine read_complex_matrix_raw(h,filename)
        implicit none
        integer mh,nh
        complex*16 :: h(:,:)
        character(len=*) filename
        mh = size(h,1)
        nh = size(h,2)

        open(unit=1000,file=filename,form='unformatted')
        read(1000) h
        return
   	end subroutine read_complex_matrix_raw
    	subroutine read_complex_vector_raw(h,filename)
        implicit none
        integer mh
        complex*16 :: h(:)
        character(len=*) filename
        mh = size(h,1)

        open(unit=1000,file=filename,form='unformatted')
        read(1000) h
        return
    	end subroutine read_complex_vector_raw
    !
    !   real*8 unformatted
    !
    	subroutine read_real8_matrix_raw(h,filename)
        implicit none
        integer mh,nh
        real*8 :: h(:,:)
        character(len=*) filename
        mh = size(h,1)
        nh = size(h,2)

        open(unit=1000,file=filename,form='unformatted')
        read(1000) h
        return
    	end subroutine read_real8_matrix_raw
    	subroutine read_real8_vector_raw(h,filename)
        implicit none
        integer mh
        real*8 :: h(:)
        character(len=*) filename
        mh = size(h,1)

        open(unit=1000,file=filename,form='unformatted')
        read(1000) h
        return
    	end subroutine read_real8_vector_raw

!   integer, unformatted
    !
   	subroutine read_integer_matrix_raw(h,filename)
        implicit none
        integer h
        character(len=*) filename

        open(unit=1000,file=filename,form='unformatted')
        read(1000) h
        return
   	end subroutine read_integer_matrix_raw
    	subroutine read_integer_vector_raw(h,filename)
        implicit none
        integer mh
        integer :: h(:)
        character(len=*) filename
        mh = size(h,1)

        open(unit=1000,file=filename,form='unformatted')
        read(1000) h
        return
    	end subroutine read_integer_vector_raw

   	subroutine read_complex_tensor_raw(h,filename)
        implicit none
        integer mh,nh
        complex*16 :: h(:,:,:)
        character(len=*) filename
        mh = size(h,1)
        nh = size(h,2)

        open(unit=1000,file=filename,form='unformatted')
        read(1000) h
        return
   	end subroutine read_complex_tensor_raw
    !
    !   real*8 unformatted
    !
    	subroutine read_real8_tensor_raw(h,filename)
        implicit none
        integer mh,nh
        real*8 :: h(:,:,:)
        character(len=*) filename
        mh = size(h,1)
        nh = size(h,2)

        open(unit=1000,file=filename,form='unformatted')
        read(1000) h
        return
    	end subroutine read_real8_tensor_raw

!   integer, unformatted
    !
   	subroutine read_integer_tensor_raw(h,filename)
        implicit none
        integer h(:,:,:)
        character(len=*) filename

        open(unit=1000,file=filename,form='unformatted')
        read(1000) h
        return
   	end subroutine read_integer_tensor_raw

end module file_io_module
