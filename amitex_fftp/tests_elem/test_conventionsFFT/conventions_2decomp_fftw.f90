!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main test program for the FFT r2c/c2r interface
!  also demonstrate the use of the IO library 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program conventions_2decomp_fftw

#ifdef NEW_2DECOMP
  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
#else
  use decomp_2d  
#endif  
  use decomp_2d_fft
!  use glassman
  use decomp_2d_io
  
  use MPI
  
  implicit none
  include "fftw3.f"
  
!  integer, parameter :: nx=5, ny=5, nz=5
!  integer, parameter :: nx=2, ny=2, nz=2
!  integer, parameter :: nx=97, ny=12, nz=3
  integer, parameter :: nx=4, ny=3, nz=3
!  integer, parameter :: nx=40, ny=30, nz=30
  integer            :: p_row, p_col
  
  real(mytype), allocatable, dimension(:,:,:) :: in, in2
  complex(mytype), allocatable, dimension(:,:,:) :: out
  
  integer, dimension(3) :: fft_start, fft_end, fft_size
  
  real(mytype), dimension(nx,ny,nz) :: in_global, in_g2, in_g3
  complex(mytype), dimension(nx/2+1,ny,nz) :: out_global
  
  integer (kind=MPI_OFFSET_KIND) :: filesize, disp
  
  real(mytype) :: err,tmp,tmp2,puissF
  !integer*8 :: plan
  integer :: fh, ierror, i,j,k, n,iol
  integer*8 :: plan  

  
  call MPI_INIT(ierror)
  p_row=2
  p_col=2
  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call decomp_2d_fft_init
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute a small problem all on rank 0 as reference
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call random_number(in_global)

if(nrank==0) print *, "===========tests convention 2decomp======================"
if(nrank==0) print *, "GLOBAL : moy = ",sum(in_global)/(nx*ny*nz)
if(nrank==0) print *, "GLOBAL somme des carre = ",sum(in_global * in_global)


  if (nrank==0) then
     ! Using a 3D FFT routine supplied by this library
     !call glassman_3d_r2c(in_global,nx,ny,nz,out_global)
     
     ! If using FFTW library:
     !  - uncomment the FFTW include file & plan above
     !  - uncomment the follwing function calls
     !  - change names to dfftw... for double precision 
     call dfftw_plan_dft_r2c_3d(plan,nx,ny,nz, &
          in_global,out_global,FFTW_ESTIMATE)
     call dfftw_execute_dft_r2c(plan,in_global,out_global)

  end if

if(nrank==0) print *, "===========tests convention 2decomp======================"
if(nrank==0) print *,"GLOBAL  : (F(1,1,1))/Ntot = ", real(out_global(1,1,1) / (nx*ny*nz),mytype)
if(nrank==0) print *,"GLOBAL  : (somme des carres F)/Ntot = ", &
         2_mytype * sum(real(out_global*dconjg(out_global),mytype)) / real(nx*ny*nz,mytype) &
        - sum(real(out_global(1,:,:)*dconjg(out_global(1,:,:)),mytype)) / real(nx*ny*nz,mytype)&
        - sum(real(out_global(nx/2+1,:,:)*conjg(out_global(nx/2+1,:,:)),mytype)) / real(nx*ny*nz,mytype)
!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test the real-to-complex interface (r2c) 
  
  !  input is X-pencil real    data whose global size is nx*ny*nz
  ! output is Z-pencil complex data whose global size is (nx/2+1)*ny*nz
  allocate (in(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
  
  call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)
  allocate (out(fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  
  ! each processor gets its local portion of global data
  do k=xstart(3),xend(3)
     do j=xstart(2),xend(2)
        do i=xstart(1),xend(1)
           in(i,j,k) = in_global(i,j,k)
        end do
     end do
  end do
    
  ! compute r2c transform 
  call decomp_2d_fft_3d(in,out)
  

if (fft_start(1)==1 .AND. fft_start(2)==1 .AND. (fft_start(3)==1)) then
print *, "===========tests convention 2decomp======================"
print *,"DISTRIB F(1,1,1) = ", real(out(1,1,1) / (nx*ny*nz))
end if

!2x somme des carres sur le 1/2 espace de Fourier
  puissF=0
  tmp=0
  tmp = 2*sum(real(out*conjg(out),mytype))
  call MPI_AllReduce(tmp,puissF,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierror)
  
! correction des plans (1,:,:) et (nx/2+1,:,:) pris deux fois en compte
  tmp=0
  tmp2=0
  if(fft_start(1)==1) then
  tmp = sum(real(out(1,:,:)*conjg(out(1,:,:)),mytype))
  end if
  call MPI_AllReduce(tmp,tmp2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierror)
  puissF = puissF - tmp2

  tmp=0
  tmp2=0
  if(modulo(nx,2) == 0 .AND. fft_end(1)==nx/2+1) then
  tmp = sum(real(out(nx/2+1,:,:)*conjg(out(nx/2+1,:,:)),mytype))
  end if
  call MPI_AllReduce(tmp,tmp2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierror)
  puissF = puissF -  tmp2



  if(nrank==0) print *, "===========tests convention 2decomp======================"
  if(nrank==0) print *,"DISTRIB  : (somme des carres F)/Ntot =  = ", &
              puissF / real(nx*ny*nz,mytype)

  deallocate(in,out) 
  
  call decomp_2d_fft_finalize
  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)
  
end program conventions_2decomp_fftw
