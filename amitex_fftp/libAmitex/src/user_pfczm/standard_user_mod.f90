!123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!===================================================================================================
!
!       MODULE STANDARD_USER_MOD 
!>
!>      Module for the standard algorithm : implement BEFORE and AFTER unpas user defined procedures 
!>
!!
!!
module standard_user_mod


  use ISO_FORTRAN_ENV


! MPI AND 2DECOMP MODULES
!------------------------
  use MPI             
  use decomp_2d
  use decomp_2d_fft

! ALL AMITEX MODULES (complete list of amitex modules)
!------------------------------------------------------
  use algo_functions_mod  
  use amitex_mod
  use amitex_user_mod     
  use error_mod
  use field_mod
  use green_mod
  use io2_amitex_mod
  use io_amitex_mod    
  use linear_mod
  use loading_mod
  use material_mod
  use non_local_mod
  use param_algo_mod
  !!use resolution_mod       this module uses standard_user_mod
  !!use resolution_user_mod  no need
  use sortie_std_mod
  !!use standard_user_mod    this file
  
#ifdef OPENMP
  use omp_lib
#endif

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
  ! ================================================================
  real(mytype) :: timePF, timeM
  real(mytype) :: critPF       !criterion for phase-field problem
  integer      :: nItPF        !number of iterations for phase field problem
  real(mytype)           :: dt, t1
  integer                :: ierror
  ! ================================================================

  !-- Calcul du pas de temps
  if (load_incr==0) then
     dt = load(load_n)%time(load_incr)
  else
     dt = load(load_n)%time(load_incr) -  load(load_n)%time(load_incr-1)
  end if

  !! solve the phase-field problem (evaluate HH and damageVar, global variables)
  call resolPF_visc(ind_tps,dt,nItPF,critPF,timePF)
  call irreversibility_dBulk() !added by YC 2021.06.12

  if(nrank==0) write(Flog,"(A)")          "Phase-field problem ------------------------------------"
  if(nrank==0) write(Flog,"(A,I9)")       "    number of iterations:",nItPF
  if(nrank==0) write(Flog,"(A,E15.8)")    "    convergence residual:",critPF
  if(nrank==0) write(Flog,"(A,E15.8,A,/)")"           time duration:",timePF,"s"
if(nrank==0) write(*,"(A)")          "Phase-field problem ------------------------------------"
if(nrank==0) write(*,"(A,I9)")       "    number of iterations:",nItPF
if(nrank==0) print *, 'd_max',maxval(damageVar)

  !! interphase damage evolution (CZM)
  call resolCZM()

       !!2020.04.24
       !!special treatment for composite voxels
       !! --> if d>kk & dI>kk, then d=kk
       call damageCV_adhoc()


  !! update the damage field and internal variables (of UMAT)
  call updateIntVar()
!  call update_dInterphase()
  call update_dBulk4Interphase()

end subroutine before_unpas_user


!===========================================================================================================
!!==========================================================================================================
!> after_unpas_user : write what must be done just AFTER unpas (using the standard algorithm)
!!                       
!!                       
!!==========================================================================================================
subroutine after_unpas_user(load_n,load_incr,ind_tps)
  
  implicit none

  integer,intent(in)     :: load_n,load_incr, ind_tps

  real(mytype)           :: dt, t1
  integer                :: ierror
  integer                :: i
  character(len=200)     :: tmp_char2,tmp_char1,tmp_char   !< variable caractere temporaire
  real(mytype),dimension(algo_param%nTensDef) :: sigAV

  call mpi_barrier(mpi_comm_world,ierror)

  !-- Calcul du pas de temps
  if (load_incr==0) then
     dt = load(load_n)%time(load_incr)
  else
     dt = load(load_n)%time(load_incr) -  load(load_n)%time(load_incr-1)
  end if


 !debug -------------------------------------------------
 call field_Mean(Sig,grid%ntot,algo_param%nTensDef,sigAV)
    sigAVmax = maxval(sigAV)
 if (nrank==0) print *,'sigAV', sigAV
! if (nrank==0) print *, 'local_loading,local_loading0, dlocal_loading',  local_loading,local_loading0, dlocal_loading
! if (nrank==0) print *, 'load(load_n)%pevolution',load(load_n)%pevolution
 !2020.09.01---------------------------------------------

 
!     !if the stress is lower than a certain level, then print vtk, abort the computation
!    if (sigAVmax<sigAVmaxInHistoEver*0.5) then
!      call writeVTK(ind_tps)
!      if (nrank == 0) write(OUTPUT_UNIT,*) "[PF VTK Write] dt, t = ", dt, load(load_n)%time(load_incr)
!      call amitex_abort("Critical low stress detected, computation aborted.",2,0)
!    end if


    !update the sigAVmaxInHistoEver
    if (sigAVmax>sigAVmaxInHistoEver) then
       sigAVmaxInHistoEver = sigAVmax
    end if



!  !================================================================================
!  !                                                                      SORTIE VTK
!  !             TODO AMITEX - mettre en place une structure pour analyser les temps 
!
!  if(load_incr/=0) then ! Pas de sortie vtk pour les pas de temps fictifs
!    if(extract%tpsVTK(ind_tps))then
!      call mpi_barrier(mpi_comm_world,ierror)
!      t1 = MPI_Wtime()
!
!      ! Sortie vtk pour les variables de champ de phase
!      !---------------------------------------
!      ! a convertir si + de 1 PF
!      !Sortie vtk des champs de phase
!      write(tmp_char1,"(I4)") ind_tps
!      do i=1,NdamageVar
!        write(tmp_char2,"(I4)") i
!        tmp_char="_eta"//trim(adjustl(tmp_char2))//"_"//trim(adjustl(tmp_char1))
!   
!        call print_field_vtk(damageVar,trim(fic_vtk)//trim(tmp_char)//".vtk","eta",6,1,i,1)
!        !call write_header_vtk(trim(fic_vtk)//trim(tmp_char)//".vtk",1000, &
!        !      grid%nx,grid%ny,grid%nz,grid%dx,grid%dy,grid%dz)
!        !TEMPfield32=real(eta0(:,i))
!        !call write_bin(1,TEMPfield32,"eta",trim(fic_vtk)//trim(tmp_char)//".vtk",1000)
!      end do
!          
!      if (nrank == 0) write(OUTPUT_UNIT,*) "[PF VTK Write] dt, t = ", dt, load(load_n)%time(load_incr)
!
!      call mpi_barrier(mpi_comm_world,ierror)
!      if (nrank == 0) write(Flog,"(A,E15.8)")"Temps ecriture VTK PF : ",MPI_Wtime() - t1
!
!      !t9_it = t9_it + MPI_Wtime() - t1 TODO
!    end if
!  end if
!
end subroutine after_unpas_user


end module standard_user_mod



