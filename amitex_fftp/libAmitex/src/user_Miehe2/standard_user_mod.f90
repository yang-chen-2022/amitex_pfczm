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
#ifdef NEW_2DECOMP
  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
#else
  use decomp_2d
#endif             
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
  use amitex_mpi_mod, only: amitex_mpi_comm
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
  real(mytype) :: timePF
  real(mytype) :: critPF       !criterion for phase-field problem
  integer      :: nItPF        !number of iterations for phase field problem
  real(mytype) :: dt

  integer :: bidon
  bidon = ind_tps
  ! ================================================================

  !-- Calcul du pas de temps
  if (load_incr==0) then
     dt = load(load_n)%time(load_incr)
  else
     dt = load(load_n)%time(load_incr) -  load(load_n)%time(load_incr-1)
  end if

  !! solve the phase-field problem (evaluate HH and damageVar, global variables)
  call resolPF_visc(dt,nItPF,critPF,timePF)

  if(nrank==0) write(Flog,"(A)")          "Phase-field problem ------------------------------------"
  if(nrank==0) write(Flog,"(A,I9)")       "    number of iterations:",nItPF
  if(nrank==0) write(Flog,"(A,E15.8)")    "    convergence residual:",critPF
  if(nrank==0) write(Flog,"(A,E15.8,A,/)")"           time duration:",timePF,"s"
if(nrank==0) write(*,"(A)")          "Phase-field problem ------------------------------------"
if(nrank==0) write(*,"(A,I9)")       "    number of iterations:",nItPF

  !! update the damage field and internal variables (of UMAT)
  !! ========================================================
  call updateIntVar()



!!-------------------------------delete me
!print *, 'algo_param%filter_typeD',algo_param%filter_typeD
!print *, 'algo_param%filter_radiusD',algo_param%filter_radiusD
!  call initFreq_YC("no_filter",algo_param%filter_radiusD)
!!----------------------------------------

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
  

  !! update the cauchy stress and internal variables!TODO,this has been done
  !during the process of mechanical solution, no need to repeat it here
  !! ===============================================
  call mpi_barrier(amitex_mpi_comm,ierror)
!  Sig0 = Sig
!  do i=1,size(MattotP)
!     MattotP(i)%VarInt0(1:MattotP(i)%nVarInt,:)=MattotP(i)%VarInt(1:MattotP(i)%nVarInt,:)
!  end do

  !! update the history field
  !! ========================
  call updateHistory()     !VarInt(2) --> HH

  !-- Calcul du pas de temps
  if (load_incr==0) then
     dt = load(load_n)%time(load_incr)
  else
     dt = load(load_n)%time(load_incr) -  load(load_n)%time(load_incr-1)
  end if

  !================================================================================
  !                                                                      SORTIE VTK
  !             TODO AMITEX - mettre en place une structure pour analyser les temps 

  if(load_incr/=0) then ! Pas de sortie vtk pour les pas de temps fictifs
    if(extract%tpsVTK(ind_tps))then
      call mpi_barrier(amitex_mpi_comm,ierror)
      t1 = MPI_Wtime()

      ! Sortie vtk pour les variables de champ de phase
      !---------------------------------------
      ! a convertir si + de 1 PF
      !Sortie vtk des champs de phase
      write(tmp_char1,"(I4)") ind_tps
      do i=1,NdamageVar
        write(tmp_char2,"(I4)") i
        tmp_char="_eta"//trim(adjustl(tmp_char2))//"_"//trim(adjustl(tmp_char1))
   
        !call print_field_vtk(damageVar(:,i),trim(fic_vtk)//trim(tmp_char)//".vtk","eta"//trim(tmp_char2))
        !call write_header_vtk(trim(fic_vtk)//trim(tmp_char)//".vtk",1000, &
        !      grid%nx,grid%ny,grid%nz,grid%dx,grid%dy,grid%dz)
        !TEMPfield32=real(eta0(:,i))
        !call write_bin(1,TEMPfield32,"eta",trim(fic_vtk)//trim(tmp_char)//".vtk",1000)
      end do
          
      if (nrank == 0) write(OUTPUT_UNIT,*) "[PF VTK Write] dt, t = ", dt, load(load_n)%time(load_incr)

      call mpi_barrier(amitex_mpi_comm,ierror)
      if (nrank == 0) write(Flog,"(A,E15.8)")"Temps ecriture VTK PF : ",MPI_Wtime() - t1

      !t9_it = t9_it + MPI_Wtime() - t1 TODO
    end if
  end if

end subroutine after_unpas_user



!===================================================================================================
!>       SUBROUTINE D'AFFECTATION DES FREQUENCES
!!
!!          format (pinceaux-Z, taille globale (nx/2+1)*ny*nz)
!!
!! LIMITATIONS : GRILLE A PIXELS CUBIQUES
!!
!! \param[out] FREQ                 Tableau des frequences
!! \param[in] filter_type           Type de filtre
!! \param[in] filter_radius         Rayon du filtre
!!
!===================================================================================================
  subroutine initFreq_YC(filter_type,filter_radius)

    implicit none
    real(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3) :: FREQ_YC

    integer                               :: i,j,k,nx_2,ny_2,nz_2
    integer,dimension(1:grid%nx/2 + 1)         :: idx
    integer,dimension(1:grid%ny)               :: idy
    integer,dimension(1:grid%nz)               :: idz

    real(mytype)                  :: T1,T2,T3,DF1,DF2,DF3,U1,U2,U3,ros2,h12,h22,h32,cote_1,cote_2,cote_3,ratio,dx1,dx2,dx3
#ifdef DOUBLE_PREC
    real(mytype),parameter        :: PI = 4._mytype*DATAN(1._mytype), sqrt3 = sqrt(3._mytype)
#else
    real(mytype),parameter        :: PI = 4._mytype*ATAN(1._mytype), sqrt3 = sqrt(3._mytype)
#endif
    !> Rayon du filtre
    real(mytype),intent(in)                        :: filter_radius
    !> Type de filtre
    character(len=*),intent(in)                    :: filter_type

    ! FREQUENCES

    ! Ancienne version : on ramenait la taille maximale de la boite a 1
    !T0 = float(maxval((/nx,ny,nz/)))
    !en anisotrope on aurait pu choisir float(maxval((/nx*dx,ny*dy,nz*dz/)))
    !T1 = float(nx)/T0
    !T2 = float(ny)/T0
    !T3 = float(nz)/T0
  
    ! Nouvelle version : necessite d'appliquer un facteur multiplicatif sur les criteres
    !                    permet de prendre en compte une dimension "physique"
    T1 = float(grid%nx)*grid%dx
    T2 = float(grid%ny)*grid%dy
    T3 = float(grid%nz)*grid%dz

    DF1 = 2._mytype * PI / T1
    DF2 = 2._mytype * PI / T2
    DF3 = 2._mytype * PI / T3
    ! indices pour ramener les frequences nulles en (1,1,1)
    ! attention aux cas pairs et impairs
    nx_2 = grid%nx / 2
    idx = (/(nx_2 + i,i=1,nx_2+1) /)

    ny_2 = grid%ny / 2
    if (mod(grid%ny,2) .eq. 0) then
       idy = (/(ny_2 + i,i=1,ny_2),(i,i=1,ny_2)/)
    else
       idy = (/(ny_2 + i,i=1,ny_2+1),(i,i=1,ny_2)/)
print *, 'ny_2', ny_2
    endif

    nz_2 = grid%nz / 2
    if (mod(grid%nz,2) .eq. 0) then
       idz = (/(nz_2 + i,i=1,nz_2),(i,i=1,nz_2)/)
    else
       idz = (/(nz_2 + i,i=1,nz_2+1),(i,i=1,nz_2)/)
    endif
print *, 'nx,ny,nz=', grid%nx,grid%ny,grid%nz
print *, 'idx', idx
print *, 'idy', idy
print *, 'idz', idz
print *, 'fft_stat,fft_end',fft_start,fft_end
    if(filter_type == "no_filter")then
       do k = fft_start(3),fft_end(3)
          do j = fft_start(2),fft_end(2)
             do i = fft_start(1),fft_end(1)
                FREQ_YC(i,j,k,1) = (-nx_2 - 1 + idx(i)) * DF1
                FREQ_YC(i,j,k,2) = (-ny_2 - 1 + idy(j)) * DF2
                FREQ_YC(i,j,k,3) = (-nz_2 - 1 + idz(k)) * DF3
             end do
          end do
       end do
print *, 'FREQ_YC', FREQ_YC
    elseif(filter_type == "hexa")then
       ros2 = filter_radius/2._mytype;
       do k = fft_start(3),fft_end(3)
          do j = fft_start(2),fft_end(2)
             do i = fft_start(1),fft_end(1)
                U1 = (-nx_2 - 1 + idx(i)) * DF1
                U2 = (-ny_2 - 1 + idy(j)) * DF2
                U3 = (-nz_2 - 1 + idz(k)) * DF3
                FREQ_YC(i,j,k,1) = (grid%nx/(ros2*T1))*sin(ros2*U1*T1/grid%nx)*cos(ros2*U2*T2/grid%ny)*cos(ros2*U3*T3/grid%nz)
                FREQ_YC(i,j,k,2) = (grid%ny/(ros2*T2))*cos(ros2*U1*T1/grid%nx)*sin(ros2*U2*T2/grid%ny)*cos(ros2*U3*T3/grid%nz)
                FREQ_YC(i,j,k,3) = (grid%nz/(ros2*T3))*cos(ros2*U1*T1/grid%nx)*cos(ros2*U2*T2/grid%ny)*sin(ros2*U3*T3/grid%nz)
             end do
          end do
       end do
    elseif(filter_type == "octa")then
       ros2 = filter_radius/2._mytype
       dx1 = (T1/grid%nx)
       dx2 = (T2/grid%ny)
       dx3 = (T3/grid%nz)
       h12 = dx1*dx1
       h22 = dx2*dx2
       h32 = dx3*dx3
       cote_1 = sqrt(h22 + h32)
       cote_2 = sqrt(h12 + h32)
       cote_3 = sqrt(h12 + h22)
       ratio =  sqrt((cote_1 + cote_2 + cote_3)*(-cote_1 + cote_2 + cote_3)*&
            (cote_1 - cote_2 + cote_3)*(cote_1 + cote_2 - cote_3)) / (sqrt3*filter_radius * dx1 * dx2 * dx3)
       do k = fft_start(3),fft_end(3)
          do j = fft_start(2),fft_end(2)
             do i = fft_start(1),fft_end(1)
                U1 = (-nx_2 - 1 + idx(i)) * DF1
                U2 = (-ny_2 - 1 + idy(j)) * DF2
                U3 = (-nz_2 - 1 + idz(k)) * DF3
                FREQ_YC(i,j,k,1) = ratio * sin(ros2*U1*T1/grid%nx)
                FREQ_YC(i,j,k,2) = ratio * sin(ros2*U2*T2/grid%ny)
                FREQ_YC(i,j,k,3) = ratio * sin(ros2*U3*T3/grid%nz)
             end do
          end do
       end do
    else
       call amitex_abort("Erreur (InitFreq) : Filtre "//filter_type //" inconnu",2,0)
    end if
  end subroutine initFreq_YC



end module standard_user_mod



