!==============================================================================
!
!       MODULE ERROR_MOD : 
!> Gestion des erreurs
!!
!! 
!!
!!  
!!  Subroutines
!! - amitex_abort :          affichage du message d'erreur et du niveau d'erreur.
!! - check_amitex_abort :    interrompt l'execution en cas d'erreur detectee precedemment.
!!
!!
!==============================================================================
module error_mod

  use mpi
  use amitex_mpi_mod, only: amitex_mpi_comm

  private 
  
  !> Pointer "publiques" vers composantes de SIMU_AMITEX
  public :: error
  
  !> Fonctions publiques
  public :: amitex_abort, check_amitex_abort

  logical, pointer :: error ! initialisee a .FALSE. dans SIMU%error

contains




!==============================================================================
!
!                          SUBROUTINE AMITEX_ABORT
!> Fonction d'arret du programme en cas d'erreur.\n
!!
!!   Affiche le rang du processus concerne, le type d'erreur et le message d'erreur.\n
!!   Enregistre la presence d'erreur (ne necessitant pas l'arret immediat du programme).\n
!!   Interrompt le programme si necessaire.
!!
!!
!! \param[in]  msg: (chaine de caracteres) message d'erreur
!! \param[in]  level: (integer) niveau d'erreur
!!                   - -2  rien
!!                   - -1  information 
!!                   - 0   avertissement, ne necessite pas l'interruption du programme
!!                   - 1   erreur, necessite l'interruption du programme (peut etre differe)
!!                   - 2   erreur, necessite l'arret immediat du programme
!! \param[in]  r: (entier), optional, rang du processus concerne pour l'ecriture 
!!               des messages d'erreur
!!               si 'r' n'est pas present, tous les processus sont impliques
!!
!!         IMPORTANT
!!             en cas de level 2, meme si r est present, 
!!             TOUS les processus appellent MPI_ABORT 
!!
!!------------------------------------------------------------------------------
!! MODIFICATIONS :
!!----------------
!! 24/03/2016 - LG : 	ecriture des message d'erreur sur stdout 
!!			arret du prog. par check_amitex_abort
!==============================================================================
subroutine amitex_abort( msg, level, r)

use, intrinsic :: iso_fortran_env, only : stdout=>output_unit

  implicit none

  character(len=*),intent(in) :: msg
  integer,intent(in) :: level
  integer,optional :: r

  integer :: rank, r0
  character(len=11) :: err_type
  integer :: ierror
 

  call MPI_COMM_RANK(amitex_mpi_comm,rank,ierror)
  
 r0=rank
 if(present(r)) r0=r

 if(level==-2) then
    ! on ne fait rien 
 elseif(level==-1) then
    ! information (localisation de l'erreur)
    err_type=" "
    if(rank==r0) write(stdout,"(A,I0,2A)") "AMITEX_FFTP INFO : rank ",rank,err_type,trim(msg) 
 else if(level==0) then
    ! simple avertissement
    err_type=" warning : "
    if(rank==r0) write(stdout,"(A,I0,2A)") "AMITEX_FFTP WARNING : rank ",rank,err_type,trim(msg)
 else if(level==1)then
    ! erreur, arret differe
    error=.TRUE.
    err_type=" error : "
    if(rank==r0) write(stdout,"(A,I0,2A)") "AMITEX_FFTP ERROR : rank ",rank,err_type,trim(msg) 
 else if(level==2)then
    ! erreur, arret immediat du programme
    error=.TRUE.
    err_type = " error : "
    if(rank==r0) write(stdout,"(A,I0,2A)") "AMITEX_FFTP ERROR : rank ",rank,err_type,trim(msg)
    call MPI_BARRIER(amitex_mpi_comm,ierror)
    call check_amitex_abort()
 end if
 

end subroutine amitex_abort
!==============================================================================



!==============================================================================
!
!                       SUBROUTINE CHECK_AMITEX_ABORT
!
!>   Interrompt le programme si une erreur est survenue.
!!
!!\param[in]  r: (entier) optionnel, rang du processus concerne
!!               si 'r' n'est pas present, tous les processus sont impliques
!!
!==============================================================================
subroutine check_amitex_abort(r)
  implicit none

  integer,intent(in), optional :: r
  integer :: rank, r0, ierror
  
  call MPI_COMM_RANK(amitex_mpi_comm,rank,ierror)
  
  r0=rank
  if(present(r)) r0=r

  !! En cas d'erreur: arret de l'execution  
  if(error .and. rank==r0 )then
    call MPI_Abort(amitex_mpi_comm, 1, ierror)
    call MPI_FINALIZE(ierror)
  end if

end subroutine check_amitex_abort
!==============================================================================

end module error_mod


