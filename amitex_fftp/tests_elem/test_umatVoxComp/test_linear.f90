program test_linear
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use mpi
#ifdef NEW_2DECOMP
  use decomp_2d_constants, only: mytype
#else
  use decomp_2d, only: mytype
#endif
  use linear_mod, only: LUSolve
  implicit none

  real(mytype) :: a(3,3), b(3), original(3,3), rhs(3)
  real(mytype) :: a6(6,6), b6(6)
  real(mytype) :: small(1,1), value(1)
  integer :: ierr, info, i

  call MPI_Init(ierr)
  a = reshape([1._mytype, 1._mytype, 0._mytype, &
               1._mytype, 1._mytype, 1._mytype, &
               0._mytype, 1._mytype, 1._mytype], [3,3])
  b = [2._mytype, 3._mytype, 2._mytype]
  call LUSolve(a, b, info)
  if (info /= 0 .or. maxval(abs(b - 1._mytype)) > 1.e-12_mytype) error stop 1

  a6 = 0._mytype
  do i = 1, 6
     a6(i,modulo(i,6)+1) = real(i, mytype)
  end do
  b6 = sum(a6, dim=2)
  call LUSolve(a6, b6, info)
  if (info /= 0 .or. maxval(abs(b6 - 1._mytype)) > 1.e-12_mytype) error stop 7

  a = reshape([1._mytype, 1._mytype, 0._mytype, &
               1._mytype, 1._mytype, 0._mytype, &
               0._mytype, 0._mytype, 1._mytype], [3,3])
  b = [1._mytype, 2._mytype, 3._mytype]
  original = a
  rhs = b
  call LUSolve(a, b, info)
  if (info /= 2 .or. any(a /= original) .or. any(b /= rhs)) error stop 2

  a = 0._mytype
  a(1,1) = 1._mytype
  a(2,2) = 1._mytype
  call LUSolve(a, b, info)
  if (info /= 3) error stop 3

  small(1,1) = 1.e-20_mytype
  value = 2.e-20_mytype
  call LUSolve(small, value, info)
  if (info /= 0 .or. abs(value(1) - 2._mytype) > 1.e-12_mytype) error stop 4
  small(1,1) = 0._mytype
  call LUSolve(small, value, info)
  if (info /= 1) error stop 5

  small(1,1) = ieee_value(0._mytype, ieee_quiet_nan)
  call LUSolve(small, value, info)
  if (info /= -1) error stop 6

  print *, "LU regression tests passed"
  call MPI_Finalize(ierr)
end program test_linear
