!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main test program  for non-distributed DTT  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_DTT

  use ISO_FORTRAN_ENV

  implicit none
    
  include "fftw3.f"
  
  
!  integer, parameter :: nx=5, ny=5, nz=5
!  integer, parameter :: nx=2, ny=2, nz=2
!  integer, parameter :: nx=97, ny=12, nz=3
  integer, parameter :: nx=2, ny=2, nz=6
!  integer, parameter :: nx=40, ny=30, nz=30
    
  
  real(REAL64), dimension(nx,ny,nz)        :: in, in2, outr
  real(REAL64), dimension(nx,ny,nz-2)      :: inzred2, outrzred2
  real(REAL64), dimension(nx,ny,nz-1)      :: inzred1, outrzred1
  complex(REAL64), dimension(nx/2+1,ny,nz) :: outx
  complex(REAL64), dimension(nx,ny,nz/2+1) :: outz
  complex(REAL64), dimension(nx,ny,nz/2)         :: outzred2
  complex(REAL64), dimension(nx,ny,(nz-1)/2+1)   :: outzred1
  
  

  integer   :: nx_2,nz_2
  integer*8 :: plan  

  nx_2 = nx / 2 + 1
  nz_2 = nz / 2 + 1
  

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Generate Random Number
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call random_number(in)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! FFT/iFFT 3D
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call dfftw_plan_dft_r2c_3d(plan,nx,ny,nz, &
                             in,outx,FFTW_ESTIMATE)
  call dfftw_execute_dft_r2c(plan,in,outx)

  call dfftw_plan_dft_c2r_3d(plan,nx,ny,nz, &
                             outx,in,FFTW_ESTIMATE)
  call dfftw_execute_dft_c2r(plan,outx,in2)
  
  in2 = in2 / real(nx*ny*nz,REAL64)

print *, "------------ Error FFT/iFFT 3D"
print *, sqrt( sum(((in - in2)**2)/(nx*ny*nz)) )
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! FFT/iFFT multiple 1D-X
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call dfftw_plan_many_dft_r2c(plan, 1, nx,ny*nz, &
         in, nx, 1, nx, &                    !inembed,istride,idist   ! jth element of kth sequence is located at j*idist + k*istride
         outx, nx_2, 1, nx_2, FFTW_ESTIMATE)  !onembed,ostride,odist
         
  call dfftw_execute_dft_r2c(plan, in, outx)
!print *, outx
  call dfftw_plan_many_dft_c2r(plan, 1, nx ,ny*nz, & 
         outx,nx_2, 1, nx_2,&                !inembed,istride,idist   ! jth element of kth sequence is located at j*idist + k*istride
         in2, nx, 1, nx, FFTW_ESTIMATE)      !onembed,ostride,odist

  call dfftw_execute_dft_c2r(plan, outx, in2)

  in2 = in2 / real(nx,REAL64)

print *, "------------ Error FFT/iFFT multiple 1D-X"
print *, sqrt( sum(((in - in2)**2)/(nx*ny*nz)) )


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! FFT/iFFT multiple 1D-Z
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
   call dfftw_plan_many_dft_r2c(plan, 1, nz, nx*ny, &
         in, nz,nx*ny, 1, &                   !inembed,istride,idist ! jth element of kth sequence is located at j*idist + k*istride
         outz, nz_2, nx*ny, 1, FFTW_ESTIMATE) !onembed,ostride,odist
         
   call dfftw_execute_dft_r2c(plan, in, outz)
   
!print *, outz

   call dfftw_plan_many_dft_c2r(plan, 1, nz,nx*ny,&
         outz, nz_2, nx*ny, 1, &
         in2, nz,nx*ny, 1, FFTW_ESTIMATE)
         
  call dfftw_execute_dft_c2r(plan, outz, in2)

  in2 = in2 / real(nz,REAL64)

print *, "------------ Error FFT/iFFT multiple 1D-Z"
print *, sqrt( sum(((in - in2)**2)/(nx*ny*nz)) )


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! FFT/iFFT multiple 1D-Z omitting First and Last value 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  outz = 0
  call dfftw_plan_many_dft_r2c(plan, 1, nz-2, nx*ny, &
        in, nz,nx*ny, 1, &                   !inembed,istride,idist ! jth element of kth sequence is located at j*idist + k*istride
        outz, nz_2, nx*ny, 1, FFTW_ESTIMATE) !onembed,ostride,odist
         
  call dfftw_execute_dft_r2c(plan, in(1,1,2), outz(1,1,1))

  call dfftw_plan_many_dft_c2r(plan, 1, nz-2,nx*ny,&
         outz, nz_2, nx*ny, 1, &
         in2, nz,nx*ny, 1, FFTW_ESTIMATE)
         
  call dfftw_execute_dft_c2r(plan, outz(1,1,1), in2(1,1,2))

  in2(:,:,2:nz-1) = in2(:,:,2:nz-1) / real(nz-2,REAL64)


print *, "------------ Error FFT/iFFT multiple 1D-Z omitting First and Last Value"
print *, sqrt( sum(((in - in2)**2)/(nx*ny*(nz-2))) )


  ! FFT multiple 1D-Z only on the Zone of Interest
  inzred2 = in(:,:,2:nz-1)
  outzred2 = 0.
  call dfftw_plan_many_dft_r2c(plan, 1, nz-2, nx*ny, &
        inzred2, nz-2,nx*ny, 1, &                   !inembed,istride,idist ! jth element of kth sequence is located at j*idist + k*istride
        outzred2, nz_2-1, nx*ny, 1, FFTW_ESTIMATE)  !onembed,ostride,odist
         
  call dfftw_execute_dft_r2c(plan, inzred2, outzred2)

print *, "------------ Error (FFT multiple 1D-Z only on the Zone of Interest)  vs (FFT multiple 1D-Z omitting First and Last Value)"
print *, "                 real and imaginary part "
print *, sqrt( sum(((real(outz(:,:,1:nz_2-1),REAL64) - real(outzred2,REAL64))**2)/(nx*ny*(nz_2-1))) )
print *, sqrt( sum(((dimag(outz(:,:,1:nz_2-1)) - dimag(outzred2))**2)/(nx*ny*(nz_2-1))) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! FFT/iFFT multiple 1D-Z omitting Last value 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  outz = 0
  call dfftw_plan_many_dft_r2c(plan, 1, nz-1, nx*ny, &
        in, nz,nx*ny, 1, &                   !inembed,istride,idist ! jth element of kth sequence is located at j*idist + k*istride
        outz, nz_2, nx*ny, 1, FFTW_ESTIMATE) !onembed,ostride,odist
         
  call dfftw_execute_dft_r2c(plan, in(1,1,1), outz(1,1,1)) 

  call dfftw_plan_many_dft_c2r(plan, 1, nz-1,nx*ny,&
         outz, nz_2, nx*ny, 1, &
         in2, nz,nx*ny, 1, FFTW_ESTIMATE)
         
  call dfftw_execute_dft_c2r(plan, outz(1,1,1), in2(1,1,1))

  in2(:,:,1:nz-1) = in2(:,:,1:nz-1) / real(nz-1,REAL64)


print *, "------------ Error FFT/iFFT multiple 1D-Z omitting Last Value"
print *, sqrt( sum(((in - in2)**2)/(nx*ny*(nz-1))) )

  ! FFT multiple 1D-Z only on the Zone of Interest
  inzred1 = in(:,:,1:nz-1)
  outzred1 = 0.
  call dfftw_plan_many_dft_r2c(plan, 1, nz-1, nx*ny, &
        inzred1, nz-1,nx*ny, 1, &                   !inembed,istride,idist ! jth element of kth sequence is located at j*idist + k*istride
        outzred1, (nz-1)/2+1, nx*ny, 1, FFTW_ESTIMATE)  !onembed,ostride,odist
         
  call dfftw_execute_dft_r2c(plan, inzred1, outzred1)

print *, "------------ Error (FFT multiple 1D-Z only on the Zone of Interest)  vs (FFT multiple 1D-Z omitting Last Value)"
print *, "                 real and imaginary part "
print *, sqrt( sum(((real(outz(:,:,1:(nz-1)/2+1),REAL64) - real(outzred1,REAL64))**2)/(nx*ny*((nz-1)/2+1))) )
print *, sqrt( sum(((dimag(outz(:,:,1:(nz-1)/2+1)) - dimag(outzred1))**2)/(nx*ny*((nz-1)/2+1))) )


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DST1/iDST1 multiple 1D-Z omitting First and Last value 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  outr = 0
  in2 = in
  call dfftw_plan_many_r2r(plan, 1, nz-2, nx*ny, &
        in, nz,nx*ny, 1, &                   !inembed,istride,idist ! jth element of kth sequence is located at j*idist + k*istride
        outr, nz, nx*ny, 1, FFTW_RODFT00,FFTW_ESTIMATE)  !onembed,ostride,odist
         
  call dfftw_execute_r2r(plan, in(1,1,2), outr(1,1,1))

  call dfftw_plan_many_r2r(plan, 1, nz-2,nx*ny,&
         outr, nz, nx*ny, 1, &
         in2, nz,nx*ny, 1, FFTW_RODFT00, FFTW_ESTIMATE)
         
  call dfftw_execute_r2r(plan, outr(1,1,1), in2(1,1,2))

  in2(:,:,2:nz-1) = in2(:,:,2:nz-1) / real(2*nz-2,REAL64)   ! normalization 2(nz-1) - OK vs fftw


print *, "------------ Error DST1/iDST1 multiple 1D-Z omitting First and Last Value"
print *, sqrt( sum(((in - in2)**2)/(nx*ny*(nz-2))) )

  ! FFT multiple 1D-Z only on the Zone of Interest
  inzred2 = in(:,:,2:nz-1)
  outrzred2 = 0.

  call dfftw_plan_many_r2r(plan, 1, nz-2, nx*ny, &
        inzred2, nz-2,nx*ny, 1, &                   !inembed,istride,idist ! jth element of kth sequence is located at j*idist + k*istride
        outrzred2, nz - 2, nx*ny, 1,  FFTW_RODFT00, FFTW_ESTIMATE)  !onembed,ostride,odist
         
  call dfftw_execute_r2r(plan, inzred2, outrzred2)

print *, "--------- Error (DST multiple 1D-Z only on the Zone of Interest)  vs (DST1 multiple 1D-Z omitting First and Last Value)"
print *, sqrt( sum(((outr(:,:,1:nz-2) - outrzred2)**2)/(nx*ny*(nz-2))) )



  
end program test_DTT
