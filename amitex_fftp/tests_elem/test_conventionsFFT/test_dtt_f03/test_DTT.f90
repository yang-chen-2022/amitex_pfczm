!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main test program  for non-distributed DTT  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_DTT

  use ISO_FORTRAN_ENV
  use, intrinsic :: iso_c_binding

  
  implicit none
    
  !include "fftw3.f"
  include "fftw3.f03"
  
  
  integer, parameter :: nx=40, ny=30, nz=30
    
  
  real(REAL64), dimension(nx,ny,nz) :: in, in2, outr
  
  type(C_PTR) :: plan  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Generate Random Number
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call random_number(in)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DST1/iDST1 multiple 1D-Z omitting First and Last value 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  outr = 0
  in2 = in
  plan = fftw_plan_many_r2r(1, [nz-2], nx*ny, &             !rank, size, howmany 
        in,  [nz], nx*ny, 1, &                              !inembed,istride,idist
        outr,[nz], nx*ny, 1, [FFTW_RODFT00],FFTW_ESTIMATE)  !onembed,ostride,odist
         
  call fftw_execute_r2r(plan, in(1,1,2), outr(1,1,1))

  !plan = fftw_plan_many_r2r(1, [nz-2],nx*ny,&
  !       outr, [nz], nx*ny, 1, &
  !       in2,  [nz], nx*ny, 1, [FFTW_RODFT00], FFTW_ESTIMATE)
         
  call fftw_execute_r2r(plan, outr(1,1,1), in2(1,1,2))

  in2(:,:,2:nz-1) = in2(:,:,2:nz-1) / real(2*nz-2,REAL64)   ! normalization 2(nz-1) - OK vs fftw


  ! Output screen
  print *, "------------ Error DST1/iDST1 multiple 1D-Z omitting First and Last Value"
  print *, sqrt( sum(((in(:,:,2:nz-1) - in2(:,:,2:nz-1))**2)/(nx*ny*(nz-2))) )

  
end program test_DTT
