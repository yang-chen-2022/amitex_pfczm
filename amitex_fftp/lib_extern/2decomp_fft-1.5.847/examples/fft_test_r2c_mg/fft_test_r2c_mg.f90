!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main test program for the FFT r2c/c2r interface
!  also demonstrate the use of the IO library 
!
! Example for using the enhanced API for multigrid purpose
!---------------------------------------------------------
!
!   One call to : decomp_2d_init
!   multiple initialization of grids          : decomp_2d_fft_init(PHYSICAL_IN_X, nx, ny, nz,Igrid,Ngrid)
!                      => for write function  : get decomp. info. from get_decomp_ff_info (new function in decomp_2d_fft)
!   changing from one grid to another for FFT : associate_pointers_decomp_2d_fft(Igrid)
!   to obtain sp. and ph. pencil sizes        : decomp_2d_fft_get_size(fft_start,fft_end,fft_size, xstart0, xend0, xsize0)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program fft_test_r2c_mg

  use decomp_2d 
  use decomp_2d_fft
  use glassman
  use decomp_2d_io
  
  use MPI
  
  implicit none
  !include "fftw3.f"
  
  integer            :: nx=4, ny=4, nz=6
  integer, parameter :: p_row=2, p_col=3
  
  real(mytype), allocatable, dimension(:,:,:) :: in, in2
  complex(mytype), allocatable, dimension(:,:,:) :: out
  
  integer, dimension(3)         :: fft_start, fft_end, fft_size
  integer, dimension(3)         :: xstart0,xend0,xsize0
  type(DECOMP_INFO)             :: decomp_ph, decomp_sp
    
  real(mytype), allocatable, dimension(:,:,:) :: in_global, in_g2, in_g3
  complex(mytype), allocatable, dimension(:,:,:) :: out_global
  
  integer (kind=MPI_OFFSET_KIND) :: filesize, disp
  
  real(mytype) :: err
  !integer*8 :: plan
  integer :: fh, ierror, i,j,k, n,iol
 
  !---------------------------Declarations multi-grilles
  integer :: Igrid, Ngrid
  
  call MPI_INIT(ierror)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  
  call decomp_2d_fft_init(PHYSICAL_IN_X, nx, ny, nz,Igrid=1)
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

  call decomp_2d_fft_init(PHYSICAL_IN_X, 2*nx, 2*ny, 2*nz,Igrid=2)
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  !Ngrid = 1 ! to test check error => OK
  !Igrid = 1 ! 
  !call decomp_2d_fft_init(PHYSICAL_IN_X, 2*nx, 2*ny, 2*nz,Igrid,Ngrid)
  !call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  !Ngrid = 2 ! to test check error => OK
  !Igrid = 3 !
  !call decomp_2d_fft_init(PHYSICAL_IN_X, 2*nx, 2*ny, 2*nz,Igrid,Ngrid)
  !call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  
do Igrid = 1,2
if (nrank==0) print *, "----------------- IGRID = ", Igrid  

  if (Igrid ==2) then
    nx = 2*nx; ny = 2*ny; nz = 2*nz
  end if 
  if (allocated(in_global)) deallocate(in_global)
  if (allocated(in_g2)) deallocate(in_g2)
  if (allocated(in_g3)) deallocate(in_g3)
  if (allocated(out_global)) deallocate(out_global)
  allocate(in_global(nx,ny,nz))
  allocate(in_g2(nx,ny,nz))
  allocate(in_g3(nx,ny,nz))
  allocate(out_global(nx/2+1,ny,nz))

  !-------------------------------------associate pointers and get grid sizes (in physical and spectral space)
  call associate_pointers_decomp_2d_fft(Igrid)
  call decomp_2d_fft_get_size(fft_start,fft_end,fft_size, xstart0, xend0, xsize0)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute a small problem all on rank 0 as reference
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call random_number(in_global)
  
  if (nrank==0) then
     write(*,*) '*** Reference serial computation on rank 0 only'
     write(*,*) ' global real input'
     do i=1,nx
        write(*,20) ((in_global(i,j,k),j=1,ny),k=1,nz)
     end do
     
     ! Using a 3D FFT routine supplied by this library
     call glassman_3d_r2c(in_global,nx,ny,nz,out_global)
     
     ! If using FFTW library:
     !  - uncomment the FFTW include file & plan above
     !  - uncomment the follwing function calls
     !  - change names to dfftw... for double precision 
     !call sfftw_plan_dft_r2c_3d(plan,nx,ny,nz, &
     !     in_global,out_global,FFTW_ESTIMATE)
     !call sfftw_execute_dft_r2c(plan,in_global,out_global)
     
     write(*,*) ' global complex output'
     do i=1,nx/2+1
        write(*,10) ((out_global(i,j,k),j=1,ny),k=1,nz)
     end do
  end if
10 format(1x,6(:,'(',F5.2,',',F5.2,')'))
20 format(1x,6F5.2)
  
  ! File for testing IO
  if (Igrid == 1) then
  call MPI_FILE_OPEN(MPI_COMM_WORLD, 'fftdata1', &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
       fh, ierror)
  elseif (Igrid== 2) then
  call MPI_FILE_OPEN(MPI_COMM_WORLD, 'fftdata2', &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
       fh, ierror)
  end if
  filesize = 0_MPI_OFFSET_KIND
  call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
  disp = 0_MPI_OFFSET_KIND

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test the real-to-complex interface (r2c) 
  
  !  input is X-pencil real    data whose global size is nx*ny*nz
  ! output is Z-pencil complex data whose global size is (nx/2+1)*ny*nz
  allocate (in(xstart0(1):xend0(1),xstart0(2):xend0(2),xstart0(3):xend0(3)))
  
  allocate (out(fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  
  ! each processor gets its local portion of global data
  do k=xstart0(3),xend0(3)
     do j=xstart0(2),xend0(2)
        do i=xstart0(1),xend0(1)
           in(i,j,k) = in_global(i,j,k)
        end do
     end do
  end do
  
  ! write input to file
  call get_decomp_fft_info(decomp_ph,decomp_sp)
  call decomp_2d_write_var(fh,disp,1,in, decomp_ph)
  
  if (nrank==0) then
     write(*,*) ' '
     write(*,*) '*** Distributed computation (X-pencil input)'
     write(*,*) ' real input held by rank 0:'
     write(*,20) in
  end if
  
  ! compute r2c transform 
  call decomp_2d_fft_3d(in,out)
  
  if (nrank==0) then
     write(*,*) ' - after forward transform'
     write(*,*) ' complex output held by rank 0:'
     write(*,10) out
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test the complex-to-real interface (c2r) 

  allocate (in2(xstart0(1):xend0(1),xstart0(2):xend0(2),xstart0(3):xend0(3)))
  
  ! compute c2r transform
  call decomp_2d_fft_3d(out,in2)
  
  ! normalisation
  in2 = in2 / real(nx) / real(ny) / real(nz)
  
  ! write the data recovered by inverse FFT to file
  call get_decomp_fft_info(decomp_ph,decomp_sp)
  call decomp_2d_write_var(fh,disp,1,in2,decomp_ph)
  
  if (nrank==0) then
     write(*,*) ' - after backward transform and normalisation'
     write(*,*) ' real output held by rank 0:'
     write(*,20) in2
  end if
  
  deallocate(in,in2,out)
  !call decomp_2d_fft_finalize
  
  call MPI_FILE_CLOSE(fh,ierror)
  
  ! check on rank 0 if input data is properly recovered
  ! this also tests the IO routines
  if (nrank==0) then
     in_g2(1,1,1) = real(0., mytype)
     inquire(iolength=iol) in_g2(1,1,1)
     if (Igrid == 1) then
        OPEN(10, FILE='fftdata1', FORM='unformatted', &
          ACCESS='DIRECT', RECL=iol)
     elseif (Igrid == 2) then
        OPEN(10, FILE='fftdata2', FORM='unformatted', &
          ACCESS='DIRECT', RECL=iol)
     end if     
     n=1
     do k=1,nz
        do j=1,ny
           do i=1,nx
              read(10,rec=n) in_g2(i,j,k)
              n=n+1
           end do
        end do
     end do
     do k=1,nz
        do j=1,ny
           do i=1,nx
              read(10,rec=n) in_g3(i,j,k)
              n=n+1
           end do
        end do
     end do
     err = 0._mytype
     do k=1,nz
        do j=1,ny
           do i=1,nx
              err = err + (in_g2(i,j,k)-in_g3(i,j,k))**2
           end do
        end do
     end do
     err = err / real(nx,mytype) / real(ny,mytype) / real(nz,mytype)
     write(*,*) ' error / mesh point: ', sqrt(err)
  end if

end do
  
  call decomp_2d_fft_finalize
  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)
  
end program fft_test_r2c_mg
