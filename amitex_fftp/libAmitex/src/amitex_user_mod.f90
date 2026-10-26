!===============================================================================
!
!       MODULE AMITEX_USER_MOD (minimalist white module in src)
!>
!>   - Declaration and initialization of USER defined variables 
!>   - User-defined procedures uses in resolution_user_mod or standard_user_mod 
!>
!!
!!
!===============================================================================
module amitex_user_mod

  use ISO_FORTRAN_ENV

#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi,       only : nrank
#else
  use decomp_2d, only : nrank,mytype 
#endif  

  use amitex_mod
  use error_mod

  implicit none

  private

  public :: init_user_variables

contains


!!------------------------------------------------------------------------------
!>                                         ALLOCATE AND INITIALIZE PUBLIC FIELDS 

subroutine init_user_variables()

  implicit none

  if (nrank==0) then
     write(Flog,"(A)") "ERROR :"
     write(Flog,"(A)") "      <User ...> is defined in xml algorithm file "
     write(Flog,"(A)") "      BUT 'init_user_variables' is not implemented "

     write(output_unit,"(A)") "ERROR init_user_variables :"
     write(output_unit,"(A)") "      <User ...> is defined in xml algorithm file "
     write(output_unit,"(A)") "      BUT 'init_user_variables' is not implemented "
  end if 

end subroutine init_user_variables

end module amitex_user_mod
