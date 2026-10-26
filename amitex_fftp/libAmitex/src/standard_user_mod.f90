!123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!===================================================================================================
!
!       MODULE STANDARD_USER_MOD (minimalist white module in src)
!>      
!>
!!
!!
module standard_user_mod

  use ISO_FORTRAN_ENV

#ifdef NEW_2DECOMP
  use decomp_2d_mpi,       only : nrank
#else
  use decomp_2d, only : nrank 
#endif  

  use amitex_mod
  use error_mod

  use amitex_user_mod

  implicit none

  private 

  public :: before_unpas_user, after_unpas_user

contains

!===========================================================================================================
!!==========================================================================================================
!> before_unpas_user : write what must be done just BEFORE unpas (using the standard algorithm)
!!                       
!!                       
!!==========================================================================================================
subroutine before_unpas_user(load_n,load_incr, ind_tps)
  
  implicit none

  integer,intent(in)     :: load_n,load_incr, ind_tps

 ! trick to avoid gcc-warning
 integer                 :: bidon
 bidon = load_n + load_incr + ind_tps


  if (nrank==0) then

     write(Flog,"(A)") "ERROR :"
     write(Flog,"(A)") "      <User Algo=""standard""> is defined in xml algorithm file "
     write(Flog,"(A)") "      BUT 'before_unpas_user' is not implemented "

     write(output_unit,"(A)") "ERROR :"
     write(output_unit,"(A)") "      <User Algo=""standard"">  is defined in xml algorithm file "
     write(output_unit,"(A)") "      BUT 'before_unpas_user' is not implemented "
  end if

end subroutine before_unpas_user


!===========================================================================================================
!!==========================================================================================================
!> after_unpas_user : write what must be done just AFTER unpas (using the standard algorithm)
!!                       
!!                       
!!==========================================================================================================
subroutine after_unpas_user(load_n,load_incr, ind_tps)
  
  implicit none

  integer,intent(in)     :: load_n,load_incr, ind_tps

  ! trick to avoid gcc-warning
  integer                 :: bidon
  bidon = load_n + load_incr + ind_tps

  if (nrank==0) then

     write(Flog,"(A)") "ERROR :"
     write(Flog,"(A)") "      <User Algo=""standard""> is defined in xml algorithm file "
     write(Flog,"(A)") "      BUT 'after_unpas_user' is not implemented "

     write(output_unit,"(A)") "ERROR :"
     write(output_unit,"(A)") "      <User Algo=""standard"">  is defined in xml algorithm file "
     write(output_unit,"(A)") "      BUT 'after_unpas_user' is not implemented "
  end if

end subroutine after_unpas_user


end module standard_user_mod



