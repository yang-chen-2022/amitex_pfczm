!===================================================================================================
!
!       MODULE amitex_mpi_mod
! 
!> Managgement of MPI global configuration
!!
!!
!===================================================================================================
module amitex_mpi_mod
  use mpi, only: MPI_COMM_WORLD
  implicit none
  private

  !> Global communicator. To be set before AMITEX initialisation.
  integer :: amitex_mpi_comm = MPI_COMM_WORLD

  public :: amitex_mpi_comm

end module amitex_mpi_mod