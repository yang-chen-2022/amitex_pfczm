!===============================================================================
!
!       MODULE MIEHE_MOD : 
!>
!>      Module of high level functions for Miehe 'non-local' implementation
!>
!!      Y.Chen & L.Gelebart - 2018-2019
!!
!!
!===============================================================================
module miehe_mod  

! INTRINSIC MODULES
!------------------
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

! AMITEX MODULES 
!----------------
  use amitex_mod
  use material_mod
  use green_mod  !, only : grid
  use field_mod  !, only : fft_start,fft_end,field_fft
  use algo_functions_mod , only : act3
  use param_algo_mod !,only Nloc_param
  use error_mod
  use amitex_mpi_mod, only: amitex_mpi_comm

! IMPLICIT NONE - PRIVATE
!------------------------
  implicit none
  private

! PUBLIC FUNCTIONS
!-----------------
  public :: resolPF_visc, updateHistory

! PUBLIC VARIABLES
!-----------------

  ! Fields
  public :: AA, BB, tau, tau0, ACT3_Rdamage, ACT3_Udamage
                                      ! Fields in Real space
  public :: tauF, FreqLaplacian       ! Fields in Fourier space


!!------------------------------------------------------------------------------
!>                                                                 PUBLIC FIELDS 

  !> real space (ntot,ncomp)
  real(mytype),allocatable,dimension(:,:)        :: tau, tau0, AA, BB
  real(mytype),allocatable,dimension(:,:,:)      :: ACT3_Rdamage, ACT3_Udamage

  !> Fourier space ((nx,ny,nz,ncomp)
  complex(mytype),allocatable,dimension(:,:,:,:) :: tauF
  real(mytype),allocatable,dimension(:,:,:)      :: FreqLaplacian


!!------------------------------------------------------------------------------


contains


!==================================================================================================
!>                 RESOLPF_VISC  solve the phase-field problem
!
!  param[in]
!       -     dt   : time increment between t-1 and t (current)
!       - list_mat : list of materials concerned by the non-local model
!       - modelInd : non-local model number (NLocMod_Num in mat.xml and algo.xml)
!
!  param[out]
!       -  nitPF : number of iterations in the phase-field resolution
!       - critPF : average value of the residual of local phase-field PDE
!       - timePF : time duration used for solving the phase-field problem
!
!  Global variables 
!       -      d  =  Nloc(modelInd)%Var(:,1) : damage variable
!       -      H  =  Nloc(modelInd)%Var(:,2) : history field evaluated from umat
!       -      H+ = GNloc(modelInd)%Var(:,1) : history field max(H+_0,H)
!
!          input  : GNloc (H+)
!          output : Nloc  (d, H is not modified)
!
!==================================================================================================
subroutine resolPF_visc(load_incr,dt,dt_,nitPF,critPF,timePF,list_mat,modelInd)

  implicit none

  integer,intent(in)                  :: load_incr
  real(mytype),intent(in)             :: dt,dt_
  integer, intent(in),dimension(:)    :: list_mat
  integer,intent(in)                  :: modelInd
  real(mytype),intent(out)            :: critPF
  real(mytype),intent(out)            :: timePF
  integer,intent(out)                 :: nItPF  ! nombre d'iterations pour le pas de chargement courant

  integer                             :: ierr
  logical                             :: testCV_PF   !< test of convergence for phase-field
  real(mytype)                        :: AA0,AAmin,AAmax
  integer                             :: iact3           !< compteur pour ACV
  logical                             :: lact3           !< indice d'activation de l'ACV
  logical                             :: acv_test,  acc_CV_PF
                                                                  !< sortie act3 : 
                                                                  !< permet de ne pas stocker un champ non accelere
  integer                             :: niterMaxPF
  real(mytype)                        :: tol_critPF, sumTau,sumDiffTau

  !----------------------------------------------------------------------------------------------
  call mpi_barrier(amitex_mpi_comm,ierr)
  timePF = MPI_WTIME() - timePF

  tol_critPF = Nloc_param(modelInd)%p_real(1)
  niterMaxPF = int(Nloc_param(modelInd)%p_real(2))
  if (trim(Nloc_param(modelInd)%p_string(1))=="true") acc_CV_PF=.true.

  critPF=0.
  nItPF=0
  testCV_PF  = .true.
  timePF = 0.

  iact3=0
  lact3=.false.
  acv_test = .true.

  !! compute the terms associated to history field
  !! =============================================
  call getHistoryTerm_visc(dt,list_mat,modelInd)  !AA=1/lc**2+visc/dt/lc/gc+2*HH/lc/gc
                                                  !BB=visc/dt/lc/gc*eta0+2*HH/lc/gc

  ! "reference phase" AA0 : different choices
  ! -----------------------------------------
  call MPI_AllReduce(minval(AA),AAmin,1,real_type,MPI_MIN,amitex_mpi_comm,ierr) !min of AA
  call MPI_AllReduce(maxval(AA),AAmax,1,real_type,MPI_MAX,amitex_mpi_comm,ierr) !max of AA
  AA0 = (AAmin+AAmax) / 2._mytype  !(min+max) / 2

  !!                 Initialize Damage from previous iteration
  !! =========================================================
  if (trim(Nloc_param(modelInd)%P_string(2))=="true") then
     call initDamage(load_incr,dt,dt_,list_mat,modelInd)
  end if

  !! compute the polarization term with viscous regularization
  !! =========================================================
  !tau0 = BB - (AA-AA0)*damageVar
  tau0(:,1) = BB(:,1) - (AA(:,1)-AA0)*NLOC(modelInd)%Var(:,1)

  !! ------------------------------------------------------ loooooooooooop
  !! ---------------------------------------------------------------------
  do while (testCV_PF)

     !! stock for convergence acceleration
     !! ==================================
     !if (algo_param%acc_CV .and. acv_test) then
     if (acc_CV_PF .and. acv_test) then
         iact3 = iact3+1
         if (iact3==5) iact3 = 1
         !ACT3_Udamage(:,1,iact3) = damageVar(:,1)
         ACT3_Udamage(:,1,iact3) = NLOC(modelInd)%Var(:,1)
     end if

     !! FFT of the polarization term
     !! ============================
     call mpi_barrier(amitex_mpi_comm,ierr) !TODO:verify whether the mpi_barrier is required
     call field_fft(tau0,tauF,1,1)
     tauF = tauF / real(grid%ntot,mytype)

     !! solution in Fourier space
     !! =========================
     call mpi_barrier(amitex_mpi_comm,ierr) !TODO:verify whether the mpi_barrier is required
     call apply_greenPF_visc(AA0)  !damage variable is stocked in ConF

     !! damage variable in real space
     !! =============================
     call mpi_barrier(amitex_mpi_comm,ierr) !TODO:verify whether the mpi_barrier is required
     !call field_ifft(damageVar,tauF,1,1)
     call field_ifft(NLOC(modelInd)%Var(:,1),tauF,1,1)

          !! convergence acceleration
          !! ========================
          ! on stocke le residu
          !if (algo_param%acc_CV .and. acv_test) then 
          if (acc_CV_PF .and. acv_test) then
             !ACT3_Rdamage(:,1,iact3) = damageVar(:,1) - ACT3_Udamage(:,1,iact3)
             ACT3_Rdamage(:,1,iact3) = NLOC(modelInd)%Var(:,1) - ACT3_Udamage(:,1,iact3)
          end if

          ! on rebascule le test ACV a .true.
          if (.not. acv_test) acv_test = .true.
          ! on applique l'ACV
          lact3 = .false.

          !if (algo_param%acc_CV) then 
          if (acc_CV_PF) then
             if (nItPF>=4 .AND. modulo(nItPF-4,3)==0) then 
                lact3 = .true.
!                call act3(ACT3_Udamage,ACT3_Rdamage,damageVar,GradQD,acv_test)
                call act3(ACT3_Udamage,ACT3_Rdamage,NLOC(modelInd)%Var(:,1:1),GradQD,acv_test)
             end if
          end if

     !! compute the polarization at current time, t
     !! ===========================================
     tau(:,1) = BB(:,1) - (AA(:,1)-AA0)*NLOC(modelInd)%Var(:,1)

     !! convergence test: residual of local PDE
     !! =======================================
     critPF = sum( (tau-tau0)**2 )
     call mpi_barrier(amitex_mpi_comm,ierr)
     call MPI_AllReduce(critPF,sumDiffTau,1,real_type,MPI_SUM,amitex_mpi_comm,ierr)
     critPF = sum(tau**2)
     call mpi_barrier(amitex_mpi_comm,ierr)
     call MPI_AllReduce(critPF,sumTau,1,real_type,MPI_SUM,amitex_mpi_comm,ierr)
     if (sumTau==0.) then !at the first loading step, def=0 -> HH=0 -> Con=0
        critPF = 0.
     else
        critPF = sqrt( sumDiffTau ) / sqrt(sumTau)
     end if

     if (critPF<tol_critPF) testCV_PF = .false.

     !! update the polarization term
     !! ============================
     tau0 = tau

     !! update the increment counter
     !! ============================
     nitPF = nitPF+1

     ! exit if no convergence until (nb_it>nb_max)
     ! ===========================================
     if (nitPF>niterMaxPF) then
        testCV_PF = .false. !=> TODO the following....
        !decrease the time increment to a half (t-1,t-1/2,t)
        !re-run resolPF => output only the result of the last time step (t)
        !conditional exit: if (dt<dt_min)
     end if

  end do!while loop

  ! time elapse
  !------------
  call mpi_barrier(amitex_mpi_comm,ierr)
  timePF = MPI_WTIME() - timePF

end subroutine resolPF_visc


!==================================================================================================
!>          GET_HISTORY_TERM_VISC  
!  
!  input fields  : d in Nloc(1) 
!                  H+ in GNloc(1)
!  output fields : AA = 1/lc**2 + visc/dt/lc/gc + 2*HH/lc/gc
!                  2*HH/lc/gc + visc/dt/lc/gc*d
!==================================================================================================
subroutine getHistoryTerm_visc(dt,list_mat,modelInd)

  implicit none

  !> Input 
  real(mytype),intent(in)          :: dt
  integer, intent(in),dimension(:) :: list_mat
  integer,intent(in)               :: modelInd

  !> Local Variables
  real(mytype)                                          :: lc, visc, gc, visco
  integer,dimension(Nloc_models(modelInd)%NCoeff_nloc)  :: Ind_CoeffNloc  ! indices of coefficients in mattotp()%Coeff
  integer                                               :: i0,i      
  integer(kind=INT64)                                   :: j,k,l,m        ! loop indices


  !-------------------------------------------------------------------------------------
  ! Remark : 
  ! For non_local models : replace the simple 'do i=1,size(mattotP)'
  !               by  'do i0=1,size(list_mat) / do i=1,size(mattotP) / if (mattotP(i)%numM == list_mat(i0))' 
  !
  do i0=1,size(list_mat)
  do i=1,size(mattotP)
  if (mattotP(i)%numM == list_mat(i0)) then
     if (mattotP(i)%NPhase >1 .OR. mattotP(i)%Interphase) then
        call amitex_abort("The non-local model is not extended yet to composite voxels",2)
     end if

     l=1
     do j=1,size(MattotP(i)%Zone(:,1)) !loop over the number of zones 
        ! material coefficients
        Ind_CoeffNloc = Nloc_models(modelInd)%Ind_CoeffNloc(i0,:) !> to get the non-local model coefficients from the material coefficients 
        visc = MattotP(i)%Coeff(Ind_CoeffNloc(1),j)
        lc   = MattotP(i)%Coeff(Ind_CoeffNloc(2),j)
        gc   = MattotP(i)%Coeff(Ind_CoeffNloc(3),j)
        visco = visc/dt/lc/gc                      !> viscous term: visc/dt/lc/gc

        do k=l,MattotP(i)%zone(j,1)
           !linear index of voxel position
           m = MattotP(i)%pos(k)
           !compute A(x) = 1/lc/lc + visc/dt/lc/gc + 2*HH/lc/gc
           !AA(m,1) = 1._mytype/lc/lc + visco + 2._mytype*HH(m,1)/lc/gc
           AA(m,1) = 1._mytype/lc/lc + visco + 2._mytype*GNloc(ModelInd)%Var(m,1)/lc/gc

           !compute B(x) = 2*HH/lc/gc + visc/dt/lc/gc*eta0
           !BB(m,1) = 2._mytype*HH(m,1)/lc/gc + visco*damageVar(m,1)
           BB(m,1) = 2._mytype*GNloc(modelInd)%Var(m,1)/lc/gc + visco*Nloc(modelInd)%Var(m,1)
        end do
        l=MattotP(i)%zone(j,1)+1
     end do
  end if
  end do
  end do


end subroutine getHistoryTerm_visc


!==================================================================================================
!>          apply_GreenPF_visc  
!  
!  input fields  : tauF, polarization  
!  output fields : tauF, damage d  with  tauF = tauF / (AA0 + Freq^2)
!
!==================================================================================================

subroutine apply_GreenPF_visc(AA0)
!! solve the phase field problem in Fourier space
!  <-> apply the Green operator
!      using a polarization term without unit
!      using a viscous regularization
!  param [in]
!        -  AA0 : the average of AA, AA = 1/lc/lc + visc/dt/lc/gc + 2*HH/lc/gc
   implicit none
   real(mytype),intent(in)      :: AA0
   integer(kind=8)   :: i,j,h
   do i=fft_start(3),fft_end(3)
      do j=fft_start(2),fft_end(2)
         do h=fft_start(1),fft_end(1)
            tauF(h,j,i,1) = tauF(h,j,i,1) / ( AA0 + (FREQ(h,j,i,1)**2+FREQ(h,j,i,2)**2+FREQ(h,j,i,3)**2) )
            !tauF(h,j,i,1) = tauF(h,j,i,1) / ( AA0 - FreqLaplacian(h,j,i) )
         end do
      end do
   end do
end subroutine apply_GreenPF_visc


!==================================================================================================
!>       UPDATE HISTORY VARIABLE : HH = max(HH,Varint_H) in MattotP(k)
!!                                 if k is in the non-local model
!!
!! \param[in] list_mat            integer array (nmat)
!!                                list of nmat material IDs concerned by the non-local model 
!! \param[in] list_varInt_Nloc    integer array (nmat,Ngnloc)
!!                                for each material i (ID=list_material(i)), list_varInt_Nloc(i,:)
!!                                -> list of indices in MattotP%Varint corresponding to Nloc(:)  
!==================================================================================================
subroutine updateHistory(list_mat,list_varInt_Nloc,modelInd)

  implicit none
  integer,dimension(:), intent(in)     :: list_mat
  integer,dimension(:,:),intent(in)    :: list_varInt_Nloc
  integer, intent(in)                  :: modelInd

  integer                              :: i0,i,indVarNloc
  integer(kind=INT64)                  :: j,m, k,l

  !!update the internal variables and assign the history field into one of the internal variables

  do i0 = 1,size(list_mat)
  do i=1,size(mattotP)
  if (mattotP(i)%numM == list_mat(i0)) then
     indVarNloc =  list_varInt_Nloc(i0,2)
     !indice minimum de zone
     l=1
     do j=1,size(MattotP(i)%zone(:,1))
        do k=l,MattotP(i)%zone(j,1)
           !linear index of voxel position
           m = MattotP(i)%pos(k)
           !update the history field by internal variable (VarInt(2)-->HH)
           if (GNloc(modelInd)%Var(m,1)<MattotP(i)%VarInt(indVarNloc,k)) then
              GNloc(modelInd)%Var(m,1) = MattotP(i)%VarInt(indVarNloc,k)
           end if
        end do
        l=MattotP(i)%zone(j,1)+1
     end do
  end if
  end do
  end do

end subroutine updateHistory


!==================================================================================================
!>     INITDAMAGE : propose an initial damage field from previous solutions
!!
!!     Formally : d = d0 + (d0 - d_)
!!         or     NLOC(modelInd)%Var(:,1) = mattotP(:)%Varint0 
!!                                   + (mattotP(:)%Varint0 - NLOC(modelInd)%Var(:,1)) * (dt / dt_)
!!
!==================================================================================================
subroutine initDamage(load_incr,dt,dt_,list_mat,modelInd)

  implicit none
 
  !> Inputs
  integer, intent(in)              :: load_incr
  real(mytype),intent(in)          :: dt, dt_
  integer,dimension(:), intent(in) :: list_mat
  integer, intent(in)              :: modelInd

  !> Local Variables 
  integer                          :: i0,i,indVarNloc
  integer(kind=INT64)              :: j,m, k,l

  if (load_incr > 2) then

  do i0 = 1,size(list_mat)
  do i=1,size(mattotP)
  if (mattotP(i)%numM == list_mat(i0)) then
     indVarNloc =  Nloc_models(modelInd)%Ind_VarNloc(i0,1)  ! damage field
     !indice minimum de zone
     l=1
     do j=1,size(MattotP(i)%zone(:,1))
        do k=l,MattotP(i)%zone(j,1)
           !linear index of voxel position
           m = MattotP(i)%pos(k)
           !extrapolated damage field
           Nloc(modelInd)%Var(m,1) = MattotP(i)%VarInt0(indVarNloc,k) * (1. + (dt / dt_)) &
                                   - Nloc(modelInd)%Var(m,1) * (dt / dt_)
        end do
        l=MattotP(i)%zone(j,1)+1
     end do
  end if
  end do
  end do

  end if

end subroutine initDamage







end module miehe_mod



