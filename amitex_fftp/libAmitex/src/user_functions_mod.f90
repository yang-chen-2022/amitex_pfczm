!===========================================================================================================
!
! MODULE USER_FUNCTIONS_MOD :
!
!> USER-DEFINED FUNCTIONS MODULE
!!
!!
!!  Procedures :
!!  ------------
!!
!!      load_interruption_average            : load interruption on average stress/strain criteria 
!!      load_interruption_internal_variables : load interruption on internal variables
!!
!!
!!  Usage :  
!!  -------  1 - in algorithm.xml : 
!!                      add < user_functions_mod Value="true" /> 
!!
!!           2 - generate the dynamic library libuser_functions.so with :
!!
!!           3 - befor using amitex : prepend the environment variable LD_LIBRARY_PATH
!!                      LD_LIBRARY_PATH=path_to_libuser_functions.so:$LD_LIBRARY_PATH
!!
!===========================================================================================================
module user_functions_mod

  use ISO_FORTRAN_ENV

  private

  public :: user_load_interruption_average, user_load_interruption_internal_variables 

contains
!===========================================================================================================
!> user_load_interruption_average : LOAD INTERRUPTION ON STRESS/STRAIN AVERAGES
!!
!!   \param[out] test_out       logical 
!!   \param[out] alpha          real64 fraction of the increment until interruption
!!   \param[in]  load_n         integer index of the partial loading        
!!   \param[in]  sigMoy0        stress tensor (1:6 or 9) at the begining of the increment
!!   \param[in]  sigMoy         stress tensor (1:6 or 9) at the end
!!   \param[in]  defMoy0        strain tensor
!!   \param[in]  defMoy         strain tensor
!!   \param[in]  local_load_0   local load on the cell at the begining of the increment
!!   \param[in]  local_load_1   local load on the cell at the end of the increment
!!   \param[in]  local_load_t   (-2:0) times associated to the load
!!   \param[in]  user_param     real64 vector of user parameters associated to
!!                                  <User_interruption Index="i" Value="xx"/> 
!!                                  in loading.xml files
!!
!!
!!   NOTATION : 
!!      Small strain assumption : ordering 11,22,33,12,13,23
!!                                Voigt notation : x2 for components ij (i~=j)
!!                                                 of the STRAIN tensor
!!
!!      Finite strains : stress is PK1 with ordering 11 22 33 12 13 23 21 31 32
!!
!!
!! COPY of comments from loading_mod.f90 : 
!! local_load%t1(1:algo_param%nTensDef) :                        pilotage 
!! local_load%t1(algo_param%nTensDef+1:2*algo_param%nTensDef) :  valeurs associées au pilotage
!! local_load%t1(2*algo_param%nTensDef+1) :                      temperature
!! local_load%t1(nb_param next indices) :                        parametres externes (s'ils existent)
!! local_load%t1(27 next indices) :                              gradgradU components (if exists)
!!
!! local_load%t_load(-2:0):                            (-2,-1,0) respectively previous, begin and  end times
!===========================================================================================================
subroutine user_load_interruption_average(test_out,alpha,load_n,sigMoy0,sigMoy,defMoy0,defMoy,&
                                          local_load_0,local_load_1,local_load_t, user_param)
  implicit none
  real(kind=REAL64),dimension(:), intent(in) :: sigMoy0, sigMoy, defMoy0, defMoy
  real(kind=REAL64),dimension(:), intent(in) :: local_load_0, local_load_1, local_load_t
  integer,intent(in)                         :: load_n
  real(kind=REAL64),intent(out)              :: alpha 
  logical, intent(out)                       :: test_out
  real(kind=REAL64),dimension(:),intent(in)  :: user_param

  alpha = load_n + sigMoy0(1) + sigMoy(1) + defMoy0(1) + defMoy(1) + size(user_param) &! trick to avoid gcc-warning (unused arg)
          + local_load_0(1) + local_load_1(1)+ local_load_t(1)
  alpha = 1.

  test_out =.false.

end subroutine user_load_interruption_average

!===========================================================================================================
!> user_load_interruption_internal_variables : LOAD INTERRUPTION ON INTERNAL VARIABLES
!!
!!   \param[out] test_out       logical 
!!   \param[out] alpha          real64 fraction of the increment until interruption
!!   \param[in]  load_n         integer index of the partial loading                   
!!   \param[in]  Varint         real64 array (npoints,nVarint) of internal variables
!!   \param[in]  Varint0        real64 array (npoints,nVarint) of internal variables
!!   \param[in]  numM           integer material number (from 1 to N materials)
!!   \param[in]  user_param     real64 vector of user parameters associated to
!!                                  <User_interruption Index="i" Value="xx"/> 
!!                                  in loading.xml files
!!
!===========================================================================================================
subroutine user_load_interruption_internal_variables(test_out,alpha,load_n,VarInt0,Varint,numM,user_param)
  implicit none
  integer,intent(in)                           :: numM
  real(kind=REAL64),dimension(:,:), intent(in) :: VarInt,Varint0
  integer,intent(in)                           :: load_n
  real(kind=REAL64),intent(out)                :: alpha 
  logical, intent(out)                         :: test_out
  real(kind=REAL64),dimension(:),intent(in)    :: user_param

  alpha = load_n + numM + size(Varint) + size(Varint0) + size(user_param)! trick to avoid gcc-warning (unused arg)
  alpha = 1. 

  test_out =.false.

end subroutine user_load_interruption_internal_variables

end module user_functions_mod
