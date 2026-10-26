!123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!===================================================================================================
!
!       MODULE RESOLUTION_USER_MOD (minimalist white module in src)
!>
!>      Module for the USER algorithm : implement the algorithm 
!>
!!
!!
module resolution_user_mod

  use ISO_FORTRAN_ENV

#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi,       only : nrank
#else
  use decomp_2d, only : mytype, nrank 
#endif  

  use amitex_mod
  use error_mod

  use amitex_user_mod

  implicit none

  private 

  public :: resolution_user

contains

!===========================================================================================================
!!==========================================================================================================
!> resolution_user : white function here
!!                                              
!!==========================================================================================================
subroutine resolution_user()
  
  implicit none
  
  if (nrank==0) then

     write(Flog,"(A)") "ERROR :"
     write(Flog,"(A)") "      <User Algo=""user""> is defined in xml algorithm file "
     write(Flog,"(A)") "      BUT 'resolution_user' is not implemented "

     write(output_unit,"(A)") "ERROR :"
     write(output_unit,"(A)") "      <User Algo=""user""> is defined in xml algorithm file "
     write(output_unit,"(A)") "      BUT 'resolution_user' is not implemented "
  end if

end subroutine resolution_user


end module resolution_user_mod



