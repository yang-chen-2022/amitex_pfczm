!===============================================================================
!
!       MODULE AMITEX_USER_MOD : 
!>
!>      Module for the USER defined procedure :
!>
!>                  resolution_user() in resolution_user_mod.f90
!>                  OR
!>                  before_unpas_user, after_unpas_user in standard_user_mod.f90
!!
!!
!!     TODO : EXTENSION TO MULTIPLE FIELDS
!!
!===============================================================================
module amitex_user_mod

! INTRINSIC MODULES
!------------------
  use ISO_FORTRAN_ENV

! MPI AND 2DECOMP MODULES
!------------------------
  use MPI             
  use decomp_2d
  use decomp_2d_fft

! ALL AMITEX MODULES
!-----------------------
  use algo_functions_mod  ! complete list of amitex modules
  use amitex_mod
  !!use amitex_user_mod     ! this file
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
  !!use resolution_mod       these modules use amitex_user_mod
  !!use resolution_user_mod  
  use sortie_std_mod
  !!use standard_user_mod  
  
#ifdef OPENMP
  use omp_lib
#endif

  implicit none

  private

! Declare the procedure used in amitex_fftp before the resolution (REQUIRED) 
! Initialize global user-variables
  public :: init_user_variables

! Declare global user-variables as public
  public :: damageVar, AA, BB, tau, tau0, HH, ACT3_Rdamage, ACT3_Udamage, NdamageVar
  public :: tauF, FreqLaplacian
  public :: sigAVmax, sig0AVmax, sigAVmaxInHisto, trigger, sigAVmaxInHistoEver

! Declare global procedures as public (used in resolution_user_mod or standard_user_mod)
  public :: initDamageHH,resolPF_visc,getHistoryTerm_visc,apply_GreenPF_visc
  public :: updateIntVar, initFreqLaplacian_pfczm
  public :: update_dInterphase,resolCZM,update_dBulk4Interphase,damageCV_adhoc
  public :: rs,tred1,tred2,tqlrat,tql2,pythag 
  public :: irreversibility_dBulk


!!------------------------------------------------------------------------------
!>                                                                 PUBLIC FIELDS 

  !> real space (ntot,ncomp)
  real(mytype),allocatable,dimension(:,:)        :: damageVar, tau, tau0, HH, AA, BB
  real(mytype),allocatable,dimension(:,:,:)      :: ACT3_Rdamage, ACT3_Udamage

  !> Fourier space ((nx,ny,nz,ncomp)
  complex(mytype),allocatable,dimension(:,:,:,:) :: tauF
  real(mytype),allocatable,dimension(:,:,:)      :: FreqLaplacian


!!------------------------------------------------------------------------------
!>		          			  parameters for auto-print vtks
  real(mytype)           :: sigAVmax,sig0AVmax,sigAVmaxInHisto,sigAVmaxInHistoEver
  logical                :: trigger

!!------------------------------------------------------------------------------
!>                                                                       indices
  integer :: i, NdamageVar


contains



!!------------------------------------------------------------------------------
!>                                         ALLOCATE AND INITIALIZE PUBLIC FIELDS 

subroutine init_user_variables()

  implicit none
  integer :: alloc_stat

  !!------------------------------------------------------------------------------
  !>                                                     INITIALISATION DES CHAMPS
  !>                                                 POUR LE PHASE-FIELD YANG CHEN

  ! strain energy history: in real space
  allocate(HH(xsize(1)*xsize(2)*xsize(3),1),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_variables 1)",2)
  HH=0.

  ! damage variable
  allocate(damageVar(xsize(1)*xsize(2)*xsize(3),1),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_variables 2)",2)

  ! polarization term for damage problem
  allocate(tau(xsize(1)*xsize(2)*xsize(3),1),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_variables 3)",2)
  tau=0.
  allocate(tau0(xsize(1)*xsize(2)*xsize(3),1),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_variables 3)",2)
  tau0=0.

  ! in Fourier space
  allocate(tauF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1),&
           stat=alloc_stat)   
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_variables 4)",2)
  tauF=0.
 
  ! process variables for computing the polarization term tau
  allocate(AA(xsize(1)*xsize(2)*xsize(3),1),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_variables 5)",2)
  AA=0.
  allocate(BB(xsize(1)*xsize(2)*xsize(3),1),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_variables 5)",2)
  BB=0.
  
  ! frequencies for Laplacian operator
  allocate(FreqLaplacian(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_variables 6)",2)
  FreqLaplacian=0.
 
  ! storage variables for convergence acceleration of phase-field solution
  if (user_param%p_string(1)=="true") then
     allocate(ACT3_Rdamage(xsize(1)*xsize(2)*xsize(3),1,4),stat=alloc_stat)
     allocate(ACT3_Udamage(xsize(1)*xsize(2)*xsize(3),1,4),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_variables 8)",2)
     ACT3_Rdamage=0.
     ACT3_Udamage=0.
  end if

  ! initialisation
  NdamageVar=1

  ! initialize the damage variable(s)
  call initDamageHH()

  ! initialize the lame coeff of interphase
  call initInterphasePROPS2()

  ! initialise the freqency for laplacian operator
  call initFreqLaplacian_pfczm()

  ! parameters for auto-printing vtk
  trigger = .true.
  sigAVmax = 0.
  sig0AVmax = 0.
  sigAVmaxInHisto = 0.
  sigAVmaxInHistoEver = 0.

end subroutine init_user_variables





! ****************************************************************************************************
! ****************************************************************************************************
! ****************************************************************************************************
! ****************************************************************************************************
! ****************************************************************************************************
! ****************************************************************************************************
! ****************************************************************************************************
! ****************************************************************************************************
! ****************************************************************************************************
! ****************************************************************************************************
! ====================================================================================================
!                                 subroutines for resolutionPF_YC
! ====================================================================================================

subroutine initFreqLaplacian_pfczm()
!> initialize the frequency square that corresponds to laplaian operator
!   - inspired from the subroutine "field_second_order_partial_centeredF"
! schema de derivation : DIFFERENCE FINIE CENTREE "approximation 27 voxels"
  implicit none
  integer                     :: i,j,k
#ifdef DOUBLE_PREC
  real(mytype),parameter      :: PI2 = 8._mytype*DATAN(1._mytype) !2*pi
#else
  real(mytype),parameter      :: PI2 = 8._mytype*ATAN(1._mytype)
#endif
  do k = fft_start(3),fft_end(3)
     do j = fft_start(2),fft_end(2)
        do i = fft_start(1),fft_end(1)
           FreqLaplacian(i,j,k) = (2._mytype/(grid%dx*grid%dx)) * (cos(PI2*(i-1)/grid%nx)-1._mytype) &
                                + (2._mytype/(grid%dy*grid%dy)) * (cos(PI2*(j-1)/grid%ny)-1._mytype) &
                                + (2._mytype/(grid%dz*grid%dz)) * (cos(PI2*(k-1)/grid%nz)-1._mytype)
        end do
     end do
  end do
end subroutine initFreqLaplacian_pfczm

subroutine initDamageHH()
!!Initialize the damage and history fields
  implicit none
  integer :: i,j,k,l,m, ii
  integer :: id  !indices
  real(mytype),allocatable,dimension(:)  :: fv
  real(mytype)  :: fvtot
  integer,parameter :: id_d=1, id_h=2 !need to be consistent with "mat_*.xml"
  !initialize the damage field by internal variable
  !------------------------------------------------
  do i=1,size(mattotP)
     !indice minimum de zone
     l=1
     !no initialisation for interphase materials
     if (mattotP(i)%Interphase) cycle !no calculation for interphase material
     !no initialisation for composite voxels
     if (mattotP(i)%Nphase>1) then !for composite voxels-----------------------------------------------------
        allocate(fv(mattotP(i)%Nphase))
        !update the current internal variable
        MattotP(i)%VarInt=MattotP(i)%VarInt0
        !initialise the damageVar
        do j=1,size(MattotP(i)%zone(:,1))
           !volume fraction of each phase
           do ii=1,mattotP(i)%Nphase-1
              id = 1 + mattotP(i)%Nphase + ii !index to fv_i
              fv(ii) = dble(mattotP(i)%Coeff(id,j))
           end do
           fv(mattotP(i)%Nphase) = 1._mytype-sum(fv(1:mattotP(i)%Nphase-1))
           do k=l,MattotP(i)%zone(j,1)
              !linear index of voxel position
              m = MattotP(i)%pos(k)
              !initialize the damageVar and history field (not applied to interphase)
              id = mattotP(i)%Nphase + 9
              if (trim(mattotP(i)%LawName_comp(1))=='czm_alfanosacco_explicit') then
                 damageVar(m,1) = 0
                 HH(m,1) = 0
                 fvtot = 0
              else
                 damageVar(m,1) = dble(mattotP(i)%VarInt0(id+id_d,k))*fv(1)
                 HH(m,1) =  dble(mattotP(i)%VarInt0(id+id_h,k))*fv(1)
                 fvtot = fv(1)
              end if
              do ii=2,mattotP(i)%Nphase
                 id = id + mattotP(i)%VarInt(ii-1,k) + 9
                 if (trim(mattotP(i)%LawName_comp(ii))=='czm_alfanosacco_explicit') cycle
                 damageVar(m,1) = damageVar(m,1) + dble(mattotP(i)%VarInt0(id+id_d,k))*fv(ii)
                 HH(m,1) = HH(m,1) + dble(MattotP(i)%VarInt0(id+id_h,k))*fv(ii)
                 fvtot = fvtot + fv(ii)
              end do
              damageVar(m,1) = damageVar(m,1) / fvtot
              HH(m,1) = HH(m,1) / fvtot
           end do
           l=MattotP(i)%zone(j,1)+1
        end do
        deallocate(fv)
     else                          !for homogeneous voxels---------------------------------------------------
        !update the current internal variable
        MattotP(i)%VarInt=MattotP(i)%VarInt0
        !initialise the damage variable and history field
        do j=1,size(MattotP(i)%zone(:,1))
           do k=l,MattotP(i)%zone(j,1)
              !linear index of voxel position
              m = MattotP(i)%pos(k)
              !initialize the damage field by internal variable
              damageVar(m,1) = dble(MattotP(i)%VarInt0(1,k))
              HH(m,1) = dble(MattotP(i)%VarInt0(2,k))
           end do
           l=MattotP(i)%zone(j,1)+1
        end do
     end if
  end do
end subroutine initDamageHH

subroutine initInterphasePROPS()
!Initialize the PROPS for the interphase umat COEFS
!for interphase only
  implicit none
  integer :: i,j,l, ii
  integer :: id
  integer,parameter :: NPROPS = 7
  real(mytype),dimension(NPROPS) :: PROPS
  real(mytype)             :: Kn,Ks
  real(mytype),parameter   :: pi=3.141592653589793238462643383279502884
  real(mytype),allocatable,dimension(:)  :: fv
  do i=1,size(mattotP)
     !indice minimum de zone
     l=1
     !initialisation for composite voxels only
     if (mattotP(i)%Nphase>1) then 
        allocate(fv(mattotP(i)%Nphase))
        do j=1,size(MattotP(i)%zone(:,1))
           fv(1:mattotP(i)%Nphase-1) = mattotP(i)%Coeff(mattotP(i)%Nphase+2:mattotP(i)%Nphase*2,j)
           fv(mattotP(i)%Nphase) = 1. - sum(fv(1:mattotP(i)%Nphase-1))
           id = 2*mattotP(i)%Nphase + 8 !index to the 1st Coeff_i 
           if (trim(mattotP(i)%LawName_comp(1))=='czm_alfanosacco_explicit') then
              !list of property coefficients
              PROPS = mattotP(i)%Coeff(id+1:id+NPROPS,j)
              !stiffness of CZM
              Kn = PROPS(1)**2 / (1.-PROPS(5)) / 2. / PROPS(2)
              Ks = PROPS(3)**2 / (1.-PROPS(5)) / 2. / PROPS(4)
              !equivalent Lame coefficients (lambda,mu)
              mattotP(i)%Coeff(id-1,j) = (Kn-2.*Ks)*mattotP(i)%Coeff(id+7,j)
              mattotP(i)%Coeff(id,j) = Ks*mattotP(i)%Coeff(id+7,j)
           else
              do ii=2,mattotP(i)%Nphase
                 id = id + mattotP(i)%Coeff(ii,j) + 2    !index to the 1st Coeff_i
                 if (trim(mattotP(i)%LawName_comp(ii))=='czm_alfanosacco_explicit') then
                    !list of property coefficients
                    PROPS = mattotP(i)%Coeff(id+1:id+NPROPS,j)
                    !stiffness of CZM
                    Kn = PROPS(1)**2 / (1.-PROPS(5)) / 2. / PROPS(2)
                    Ks = PROPS(3)**2 / (1.-PROPS(5)) / 2. / PROPS(4)
                    !equivalent Lame coefficients (lambda,mu)
                    mattotP(i)%Coeff(id-1,j) = (Kn-2.*Ks)*mattotP(i)%Coeff(id+7,j)
                                               !NOTE,this could be negative !
                    mattotP(i)%Coeff(id,j) = Ks*mattotP(i)%Coeff(id+7,j)
                 end if
              end do
           end if
           l=MattotP(i)%zone(j,1)+1
        end do
        deallocate(fv)
     end if
  end do
end subroutine initInterphasePROPS

subroutine initInterphasePROPS2()
!Initialize the PROPS for the interphase umat COEFS
!for interphase only
  implicit none
  integer :: i,j,l, ii
  integer :: id
  integer,parameter :: NPROPS = 7
  real(mytype),dimension(NPROPS) :: PROPS
  real(mytype)             :: Kn,Ks, nu
  real(mytype),parameter   :: pi=3.141592653589793238462643383279502884
  real(mytype),allocatable,dimension(:)  :: fv
  real(mytype),dimension(3)  :: e1,e2
  do i=1,size(mattotP)
     !indice minimum de zone
     l=1
     !initialisation for composite voxels only
     if (mattotP(i)%Nphase>1) then 
        allocate(fv(mattotP(i)%Nphase))
        do j=1,size(MattotP(i)%zone(:,1))
           fv(1:mattotP(i)%Nphase-1) = mattotP(i)%Coeff(mattotP(i)%Nphase+2:mattotP(i)%Nphase*2,j)
           fv(mattotP(i)%Nphase) = 1. - sum(fv(1:mattotP(i)%Nphase-1))
           !---------------------------------------------------------------
           !local coefficients (orientation vectors),added by YC 2020.05.04
           e1 = dble(mattotP(i)%Coeff(mattotP(i)%Nphase*2+1:mattotP(i)%Nphase*2+3,j))
           e1 = e1 / norm2(e1)
           e2 = dble(mattotP(i)%Coeff(mattotP(i)%Nphase*2+4:mattotP(i)%Nphase*2+6,j))
           e2 = e2 - (dot_product(e2,e1) * e1)
           e2 = e2 / norm2(e2)
           !---------------------------------------------------------------
           id = 2*mattotP(i)%Nphase + 8 !index to the 1st Coeff_i 
           if (trim(mattotP(i)%LawName_comp(1))=='czm_alfanosacco_explicit') then
              !list of property coefficients
              PROPS = mattotP(i)%Coeff(id+1:id+NPROPS,j)
              !characteristic length (thickness of the interphase in the voxel)
              mattotP(i)%Coeff(id+7,j) = fv(1) * (grid%dx*grid%dy*grid%dz)**(1.0/3.0) !2026.10.23
              !stiffness of CZM
              Kn = PROPS(1)**2 / (1.-PROPS(5)) / 2. / PROPS(2)
              Ks = PROPS(3)**2 / (1.-PROPS(5)) / 2. / PROPS(4)
              !equivalent Lame coefficients (lambda,mu)
              nu = 0.2 !ATTENTION, here we assume the Poisson's ratio is 0.2
              mattotP(i)%Coeff(id-1,j) = Kn*mattotP(i)%Coeff(id+7,j)*nu/(1+nu)/(1-2*nu)
              mattotP(i)%Coeff(id,j) = Ks*mattotP(i)%Coeff(id+7,j)
              !input the orientation vectors,added by YC 2020.05.04
              mattotP(i)%Coeff(id+8:id+10,j) = e1
              mattotP(i)%Coeff(id+11:id+13,j) = e2
           else
              do ii=2,mattotP(i)%Nphase
                 id = id + mattotP(i)%Coeff(ii,j) + 2    !index to the 1st Coeff_i
                 if (trim(mattotP(i)%LawName_comp(ii))=='czm_alfanosacco_explicit') then
                    !list of property coefficients
                    PROPS = mattotP(i)%Coeff(id+1:id+NPROPS,j)
                    !characteristic length (thickness of the interphase in the voxel)
                    mattotP(i)%Coeff(id+7,j) = fv(ii) * (grid%dx*grid%dy*grid%dz)**(1.0/3.0) !2026.10.23
                    !stiffness of CZM
                    Kn = PROPS(1)**2 / (1.-PROPS(5)) / 2. / PROPS(2)
                    Ks = PROPS(3)**2 / (1.-PROPS(5)) / 2. / PROPS(4)
                    !equivalent Lame coefficients (lambda,mu)
                    nu = 0.2 !ATTENTION, here we assume the Poisson's ratio is 0.2
                    mattotP(i)%Coeff(id-1,j) = Kn*mattotP(i)%Coeff(id+7,j)*nu/(1+nu)/(1-2*nu)
                    mattotP(i)%Coeff(id,j) = Ks*mattotP(i)%Coeff(id+7,j)
                    !input the orientation vectors,added by YC 2020.05.04
                    mattotP(i)%Coeff(id+8:id+10,j) = e1
                    mattotP(i)%Coeff(id+11:id+13,j) = e2
                 end if
              end do
           end if
           l=MattotP(i)%zone(j,1)+1
        end do
        deallocate(fv)
     end if
  end do
end subroutine initInterphasePROPS2

subroutine resolPF_visc(ind_tps,dt,nitPF,critPF,timePF)
!! solve the phase-field problem
!  param[in]
!       - ind_tps: index of loading step
!       -     dt : time increment between t-1 and t (current)
!       -   visc : viscous parameter
!       -     lc : characteristic length
!       -niterMaxPF: max number of iteration of phase-field solution
!  param[inout] global variables in amitex_mod
!       -     HH : history field
!       -damageVar : damage variable
!  param[out]
!       -  nitPF : number of iterations in the phase-field resolution
!       - critPF : average value of the residual of local phase-field PDE
!       - timePF : time duration used for solving the phase-field problem
  implicit none
  integer,intent(in)                  :: ind_tps !< indice du pas de temps courant
  real(mytype),intent(in)             :: dt
  integer                             :: niterMaxPF
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
  real(mytype)                        :: tol_critPF, sumTau,sumDiffTau
  real(mytype)                        :: t1
  integer                             :: i,j,k     !< indice de boucle 
  integer                             :: m,l       !index of voxel, minimum number of zone


  critPF=0.
  tol_critPF = user_param%p_real(1)
  nItPF=0
  testCV_PF  = .true.
  timePF = 0.

  iact3=0
  lact3=.false.
  acv_test = .true.

  acc_CV_PF = user_param%p_string(1)=="true"
  niterMaxPF = user_param%p_real(2)


  !! compute the terms associated to history field
  !! =============================================
  call getHistoryTerm_visc(dt)  !AA=1/lc**2+visc/dt/lc/gc+2*HH/lc/gc
                                !BB=visc/dt/lc/gc*eta0+2*HH/lc/gc

  ! "reference phase" AA0 : different choices
  ! -----------------------------------------
  call MPI_AllReduce(minval(AA),AAmin,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierr) !min of AA
  call MPI_AllReduce(maxval(AA),AAmax,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr) !max of AA
  AA0 = (AAmin+AAmax) / 2._mytype  !(min+max) / 2

  !! compute the polarization term with viscous regularization
  !! =========================================================
  tau0 = BB - (AA-AA0)*damageVar

  !! ------------------------------------------------------ loooooooooooop
  !! ---------------------------------------------------------------------
  do while (testCV_PF)
     t1 = MPI_WTIME()

           !! stock for convergence acceleration
           !! ==================================
           !if (algo_param%acc_CV .and. acv_test) then
           if (acc_CV_PF .and. acv_test) then
              iact3 = iact3+1
              if (iact3==5) iact3 = 1
              ACT3_Udamage(:,1,iact3) = damageVar(:,1)
           end if

     !! FFT of the polarization term
     !! ============================
     call mpi_barrier(mpi_comm_world,ierr) !TODO:verify whether the mpi_barrier is required
     call field_fft(tau0,tauF,1,1)
     tauF = tauF / real(grid%ntot,mytype)

     !! solution in Fourier space
     !! =========================
     call mpi_barrier(mpi_comm_world,ierr) !TODO:verify whether the mpi_barrier is required
     if (user_param%p_string(2)=="defaultLap") then
        call apply_greenPF_visc_defaultLap(AA0)  !damage variable is stocked in ConF
     else
        call apply_greenPF_visc(AA0)  !damage variable is stocked in ConF
     end if

     !! damage variable in real space
     !! =============================
     call mpi_barrier(mpi_comm_world,ierr) !TODO:verify whether the mpi_barrier is required
     call field_ifft(damageVar,tauF,1,1)

          !! convergence acceleration
          !! ========================
          ! on stocke le residu
          !if (algo_param%acc_CV .and. acv_test) then 
          if (acc_CV_PF .and. acv_test) then
             ACT3_Rdamage(:,1,iact3) = damageVar(:,1) - ACT3_Udamage(:,1,iact3)
          end if

          ! on rebascule le test ACV a .true.
          if (.not. acv_test) acv_test = .true.
          ! on applique l'ACV
          lact3 = .false.

          !if (algo_param%acc_CV) then 
          if (acc_CV_PF) then
             if (nItPF>=4 .AND. modulo(nItPF-4,3)==0) then 
                lact3 = .true.
                call act3(ACT3_Udamage,ACT3_Rdamage,damageVar,GradQD,acv_test)
             end if
          end if

     !! compute the polarization at current time, t
     !! ===========================================
     tau = BB - (AA-AA0)*damageVar

     !! convergence test: residual of local PDE
     !! =======================================
     critPF = sum( (tau-tau0)**2 )
     call mpi_barrier(mpi_comm_world,ierr)
     call MPI_AllReduce(critPF,sumDiffTau,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
     critPF = sum(tau**2)
     call mpi_barrier(mpi_comm_world,ierr)
     call MPI_AllReduce(critPF,sumTau,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
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

     !! time counter
     !! ============
     call mpi_barrier(mpi_comm_world,ierr)
     timePF = timePF + MPI_WTIME() - t1

     ! exit if no convergence until (nb_it>nb_max)
     ! ===========================================
     if (nitPF>niterMaxPF) then
        testCV_PF = .false. !=> TODO the following....
        !decrease the time increment to a half (t-1,t-1/2,t)
        !re-run resolPF => output only the result of the last time step (t)
        !conditional exit: if (dt<dt_min)
     end if

  end do!while loop
end subroutine resolPF_visc

subroutine getHistoryTerm_visc(dt)
!! compute the right-hand side of the phase-field PDE
!! using a form without unit
!  param[in]
!       -   dt : time increment between t-1 and t
!       -   HH : history field
!       - damageVar : damage field at the previous time step
!  param[out]
!       -   AA : viscous term, AA = 1/lc**2 + visc/dt/lc/gc + 2*HH/lc/gc
!       -   BB : viscous term, BB = 2*HH/lc/gc + visc/dt/lc/gc*damageVar
  implicit none
  real(mytype),intent(in) :: dt
  real(mytype)            :: lc,gc,visc,visco,psiP,fvtot
  integer               :: i,j,k,ii     !< indice de boucle 
  integer               :: m,l       !index of voxel, minimum number of zone
  integer :: id  !indices
  real(mytype),allocatable,dimension(:)  :: fv
  integer,parameter :: id_visc=3, id_lc=4, id_gc=5 !need to be consistent with "mat_*.xml"
  integer,parameter :: id_psiP=2        !need to be consistent with "mat_*.xml"
  !                                                                   Loop of every voxel
  ! -------------------------------------------------------------------------------------
  do i=1,size(mattotP)
     if (mattotP(i)%Interphase) cycle !no action for interphase material
     !indice minimum de zone
     l=1
     if (mattotP(i)%Nphase>1) then   !for composite voxels --------------------------------------------------
        allocate(fv(mattotP(i)%Nphase))
        do j=1,size(MattotP(i)%zone(:,1))
           !volume fraction of each phase
           do ii=1,mattotP(i)%Nphase-1
              id = 1 + mattotP(i)%Nphase + ii !index to fv_i
              fv(ii) = dble(mattotP(i)%Coeff(id,j))
           end do
           fv(mattotP(i)%Nphase) = 1._mytype-sum(fv(1:mattotP(i)%Nphase-1))
           !average values of the coefficients lc, gc
           id = 2*mattotP(i)%Nphase + 8               !index to the 1st Coeff_i
           if (trim(mattotP(i)%LawName_comp(1))/='czm_alfanosacco_explicit') then
              visc = dble(mattotP(i)%Coeff(id+id_visc,j)) * fv(1)
              lc = dble(mattotP(i)%Coeff(id+id_lc,j)) * fv(1)
              gc = dble(mattotP(i)%Coeff(id+id_gc,j)) * fv(1)
              fvtot = fv(1)
           else
              visc = 0
              lc = 0
              gc = 0
              fvtot = 0
           end if
           do ii=2,mattotP(i)%Nphase
              id = id + mattotP(i)%Coeff(ii,j) + 2    !index to the 1st Coeff_i
              if (trim(mattotP(i)%LawName_comp(ii))=='czm_alfanosacco_explicit') cycle
              visc = visc + dble(mattotP(i)%Coeff(id+id_visc,j)) * fv(ii)
              lc = lc + dble(mattotP(i)%Coeff(id+id_lc,j)) * fv(ii)
              gc = gc + dble(mattotP(i)%Coeff(id+id_gc,j)) * fv(ii)
              fvtot = fvtot + fv(ii)
           end do
           visc = visc / fvtot
           lc = lc / fvtot
           gc = gc / fvtot
           !viscous term
           visco = visc/dt/lc/gc
           !total strain energy of the composite voxel
           do k=l,MattotP(i)%zone(j,1)
              !linear index of voxel position
              m = MattotP(i)%pos(k)
              !strain energy
              id = mattotP(i)%Nphase + 9
              if (trim(mattotP(i)%LawName_comp(1))/='czm_alfanosacco_explicit') then
                 psiP = dble(mattotP(i)%VarInt(id+id_psiP,k)) * fv(1)
                 fvtot = fv(1)
              else
                 psiP = 0
                 fvtot = 0
              end if
              do ii=2,mattotP(i)%Nphase
                 id = id + mattotP(i)%VarInt(ii-1,k) + 9 !index to the 1st VarInt_i
                 if (trim(mattotP(i)%LawName_comp(ii))=='czm_alfanosacco_explicit') cycle
                 psiP = psiP + dble(mattotP(i)%VarInt(id+id_psiP,k)) * fv(ii)
                 fvtot = fvtot + fv(ii)
              end do
              psiP = psiP / fvtot
              !history field
              if (psiP>HH(m,1))    HH(m,1)=psiP
              !A(x) = 1/lc/lc + 2*HH/lc/gc
              AA(m,1) = 1._mytype/lc/lc + visco + 2._mytype*HH(m,1)/lc/gc
              !B(x) = 2*HH/lc/gc
              BB(m,1) = 2._mytype*HH(m,1)/lc/gc + visco*damageVar(m,1)
           end do
           l=MattotP(i)%zone(j,1)+1
        end do
        deallocate(fv)
     else                            !for homogeneous voxels ------------------------------------------------
        do j=1,size(MattotP(i)%zone(:,1))
           !fracture energy
           visc = dble(MattotP(i)%Coeff(MattotP(i)%Ncoeff-2,j))
           lc = dble(MattotP(i)%Coeff(MattotP(i)%Ncoeff-1,j))
           gc = dble(MattotP(i)%Coeff(MattotP(i)%Ncoeff,j))
           !viscous term: visc/dt/lc/gc
           visco = visc/dt/lc/gc
           do k=l,MattotP(i)%zone(j,1)
              !linear index of voxel position
              m = MattotP(i)%pos(k)
              !history field
              if (HH(m,1)<MattotP(i)%VarInt(2,k)) then
                 HH(m,1) = MattotP(i)%VarInt(2,k)
              end if
              !compute A(x) = 1/lc/lc + visc/dt/lc/gc + 2*HH/lc/gc
              AA(m,1) = 1._mytype/lc/lc + visco + 2._mytype*HH(m,1)/lc/gc
              !compute B(x) = 2*HH/lc/gc + visc/dt/lc/gc*damageVar0
              BB(m,1) = 2._mytype*HH(m,1)/lc/gc + visco*damageVar(m,1)
           end do
           l=MattotP(i)%zone(j,1)+1
        end do
     end if
  end do
end subroutine getHistoryTerm_visc

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
!            tauF(h,j,i,1) = tauF(h,j,i,1) / ( AA0 + (FREQ(h,j,i,1)**2+FREQ(h,j,i,2)**2+FREQ(h,j,i,3)**2) )
            tauF(h,j,i,1) = tauF(h,j,i,1) / ( AA0 - FreqLaplacian(h,j,i) )
         end do
      end do
   end do
end subroutine apply_GreenPF_visc

subroutine apply_GreenPF_visc_defaultLap(AA0)
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
         end do
      end do
   end do
end subroutine apply_GreenPF_visc_defaultLap

subroutine updateIntVar()
!!update the internal variables and assign the damage field into one of the internal variables
  implicit none
  integer :: i,j,k,l,m,ii,id
  integer,parameter :: id_d=1  !need to be consistent with "mat_*.xml"
  do i=1,size(mattotP)
     !indice minimum de zone
     l=1
     if (mattotP(i)%Interphase) cycle !no calculation for interphase material
     if (mattotP(i)%Nphase>1) then   !for composite voxels --------------------------------------------------
        do j=1,size(MattotP(i)%zone(:,1))
           do k=l,MattotP(i)%zone(j,1)
              !linear index of voxel position
              m = MattotP(i)%pos(k)
              !update internal variable by damage variable (damageVar->VarInt(1))
              id = mattotP(i)%Nphase + 9
              if (trim(mattotP(i)%LawName_comp(1))/='czm_alfanosacco_explicit') then
                 mattotP(i)%VarInt(id+id_d,k) = damageVar(m,1)
              end if
              do ii=2,mattotP(i)%Nphase
                 id = id + mattotP(i)%VarInt(ii-1,k) + 9
                 if (trim(mattotP(i)%LawName_comp(ii))/='czm_alfanosacco_explicit') then
                    MattotP(i)%VarInt(id+id_d,k) = damageVar(m,1)
                 end if
              end do
              !update internal variable at t-1 (VarInt->VarInt0)
              MattotP(i)%VarInt0(:,k) = MattotP(i)%VarInt(:,k)
           end do
           l=MattotP(i)%zone(j,1)+1
        end do
     else                            !for homogeneous voxels ------------------------------------------------
        do j=1,size(MattotP(i)%zone(:,1))
           do k=l,MattotP(i)%zone(j,1)
              !linear index of voxel position
              m = MattotP(i)%pos(k)
              !update internal variable by damage variable (damageVar->VarInt(1))
              MattotP(i)%VarInt(1,k) = damageVar(m,1)
              !update internal variable at t-1 (VarInt->VarInt0)
              MattotP(i)%VarInt0(:,k) = MattotP(i)%VarInt(:,k)
           end do
           l=MattotP(i)%zone(j,1)+1
        end do
     end if 
  end do
end subroutine updateIntVar

subroutine update_dInterphase()
!! update the damage variable of interphase at t-1
  implicit none
  integer               :: i,j,k,ii     !< indice de boucle 
  integer               :: l       !index of voxel, minimum number of zone
  integer               :: m       !index of voxel position
  integer :: id  !indices
  do i=1,size(mattotP)
     !indice minimum de zone
     l=1
     if (mattotP(i)%Nphase>1) then
        do j=1,size(MattotP(i)%zone(:,1))
           do k=l,MattotP(i)%zone(j,1)
              !linear index of voxel position
              m = MattotP(i)%pos(k)
              !update the internal variables
              id = mattotP(i)%Nphase + 9
              if (trim(mattotP(i)%LawName_comp(1))=='czm_alfanosacco_explicit') then
                 MattotP(i)%VarInt(id+5,k) = MattotP(i)%VarInt(id+1,k)
                 MattotP(i)%VarInt(id+6,k) = damageVar(m,1) !damageVar is used in interphase UMAT
              else
                 do ii=2,mattotP(i)%Nphase
                    id = id + mattotP(i)%VarInt(ii-1,k)+9  !index to the 1st VarInt_i
                    if (trim(mattotP(i)%LawName_comp(ii))=='czm_alfanosacco_explicit') then
                       MattotP(i)%VarInt(id+5,k) = MattotP(i)%VarInt(id+1,k) 
                       MattotP(i)%VarInt(id+6,k) = damageVar(m,1) !damageVar is used in interphase UMAT
                       exit
                    end if
                 end do
              end if
           end do
           l=MattotP(i)%zone(j,1)+1
        end do
     end if
  end do
end subroutine update_dInterphase
subroutine update_dBulk4Interphase()
!! update the bulk damage variable for interphase UMAT
  implicit none
  integer               :: i,j,k,ii     !< indice de boucle 
  integer               :: l       !index of voxel, minimum number of zone
  integer               :: m       !index of voxel position
  integer :: id  !indices
  do i=1,size(mattotP)
     !indice minimum de zone
     l=1
     if (mattotP(i)%Nphase>1) then
        do j=1,size(MattotP(i)%zone(:,1))
           do k=l,MattotP(i)%zone(j,1)
              !linear index of voxel position
              m = MattotP(i)%pos(k)
              !update the internal variables
              id = mattotP(i)%Nphase + 9
              if (trim(mattotP(i)%LawName_comp(1))=='czm_alfanosacco_explicit') then
                 MattotP(i)%VarInt0(id+6,k) = damageVar(m,1) !damageVar is used in interphase UMAT
              else
                 do ii=2,mattotP(i)%Nphase
                    id = id + mattotP(i)%VarInt(ii-1,k)+9  !index to the 1st VarInt_i
                    if (trim(mattotP(i)%LawName_comp(ii))=='czm_alfanosacco_explicit') then
                       MattotP(i)%VarInt0(id+6,k) = damageVar(m,1) !damageVar is used in interphase UMAT
                       exit
                    end if
                 end do
              end if
           end do
           l=MattotP(i)%zone(j,1)+1
        end do
     end if
  end do
end subroutine update_dBulk4Interphase

subroutine resolCZM()
!! compute the damage variable of interphase
  implicit none
  real(mytype),dimension(14)   :: PROPS
  real(mytype)                 :: dInterphase,s01,s02,s1,s2,s2_1,s2_2,def1,def4,def5
  real(mytype),parameter       :: pi=3.141592653589793238462643383279502884
  real(mytype)                 :: lcI   !characteristic length for interface fracture
  integer               :: i,j,k,ii     !< indice de boucle 
  integer               :: m,l       !index of voxel, minimum number of zone
  integer :: id  !indices
  !                                                                   Loop of every voxel
  ! -------------------------------------------------------------------------------------
  do i=1,size(mattotP)
     !indice minimum de zone
     l=1
     if (mattotP(i)%Nphase>1) then   !only for composite voxels ---------------------------------------------
        do j=1,size(MattotP(i)%zone(:,1))
           !read the coefficients from UMAT props
           id = 2*mattotP(i)%Nphase + 8 !index to the 1st Coeff_i 
           if (trim(mattotP(i)%LawName_comp(1))=='czm_alfanosacco_explicit') then
              PROPS = dble(mattotP(i)%Coeff(int(id+1):int(id+mattotP(i)%Coeff(1+1,j)),j))
           else
              do ii=2,mattotP(i)%Nphase
                 id = id + mattotP(i)%Coeff(ii,j) + 2    !index to the 1st Coeff_i
                 if (trim(mattotP(i)%LawName_comp(ii))=='czm_alfanosacco_explicit') then
                    PROPS = dble(mattotP(i)%Coeff(int(id+1):int(id+mattotP(i)%Coeff(1+ii,j)),j))
                    exit
                 end if
              end do
           end if
           do k=l,MattotP(i)%zone(j,1)
              !linear index of voxel position
              m = MattotP(i)%pos(k)
              !compute CZM parameters (with damage interaction)
              s01 = (1.-PROPS(5)) * ( 2.*PROPS(2) / ( PROPS(1) ) )
              s02 = (1.-PROPS(5)) * ( 2.*PROPS(4) / ( PROPS(3) ) )
              !read data from internal variables of UMAT
              id = mattotP(i)%Nphase + 9
              if (trim(mattotP(i)%LawName_comp(1))=='czm_alfanosacco_explicit') then
                 !read internal variables
!                 id = id + 9                       !index to the 1st VarInt_i
                 dInterphase = dble(mattotP(i)%VarInt(id+1,k))
                 def1 = dble(mattotP(i)%VarInt(id-8,k))
                 def4 = dble(mattotP(i)%VarInt(id-7,k))
                 def5 = dble(mattotP(i)%VarInt(id-6,k))
                 !strain --> displ. jump
                 s1 = def1 * PROPS(7)
                 s2_1 = def4*PROPS(7)*0.5 
                 s2_2 = def5*PROPS(7)*0.5
                 s2 = sqrt( s2_1**2 + s2_2**2 ) !principal sliding displ
                 !damage variable of the interphase
                 call damageUpdate_(s1,s2,s01,s02,PROPS(5),dInterphase)
                 !call damageUpdate_wS(s1,s2,s01,s02,PROPS(5),PROPS(15),dInterphase)
              else
                 do ii=2,mattotP(i)%Nphase
                    id = id + mattotP(i)%VarInt(ii-1,k)  !index to the 1st stress/strain component
                    if (trim(mattotP(i)%LawName_comp(ii))=='czm_alfanosacco_explicit') then
                       !read internal variables
                       id = id + 9                       !index to the 1st VarInt_i
                       dInterphase = dble(mattotP(i)%VarInt(id+1,k))
                       def1 = dble(mattotP(i)%VarInt(id-8,k))
                       def4 = dble(mattotP(i)%VarInt(id-7,k))
                       def5 = dble(mattotP(i)%VarInt(id-6,k))
                       !strain --> displ. jump
                       s1 = def1 * PROPS(7)
                       s2_1 = def4*PROPS(7)*0.5 
                       s2_2 = def5*PROPS(7)*0.5
                       s2 = sqrt( s2_1**2 + s2_2**2 ) !principal sliding displ
                       !damage variable of the interphase
                       call damageUpdate_(s1,s2,s01,s02,PROPS(5),dInterphase)
                       !call damageUpdate_wS(s1,s2,s01,s02,PROPS(5),PROPS(15),dInterphase)
                       exit
                    end if
!                    id = id + 9                       !index to the 1st VarInt_i
                 end do
             end if
             !update the internal variables
              id = mattotP(i)%Nphase + 9
              if (trim(mattotP(i)%LawName_comp(1))=='czm_alfanosacco_explicit') then
                 !manual control of the irreversibility & viscous regularisation
                 if (dInterphase<MattotP(i)%VarInt(id+1,k)) then
                    dInterphase =  MattotP(i)%VarInt(id+1,k)
                 else
                    dInterphase = 1./(1.+PROPS(14)) * dInterphase + &
                                  PROPS(14)/(1.+PROPS(14)) * MattotP(i)%VarInt(id+1,k)
                 end if
                 !update
                 MattotP(i)%VarInt(id+1,k) = dInterphase
                 MattotP(i)%VarInt(id+5,k) = dInterphase
              else
                 do ii=2,mattotP(i)%Nphase
                    id = id + mattotP(i)%VarInt(ii-1,k)+9  !index to the 1st VarInt_i
                    if (trim(mattotP(i)%LawName_comp(ii))=='czm_alfanosacco_explicit') then
                       !manual control of the irreversibility & viscous regularisation
                       if (dInterphase<MattotP(i)%VarInt(id+1,k)) then
                          dInterphase =  MattotP(i)%VarInt(id+1,k)
                       else
                          dInterphase = 1./(1.+PROPS(14)) * dInterphase + &
                                        PROPS(14)/(1.+PROPS(14)) * MattotP(i)%VarInt(id+1,k)
                       end if
                       !update
                       MattotP(i)%VarInt(id+1,k) = dInterphase
                       MattotP(i)%VarInt(id+5,k) = dInterphase
                       exit
                    end if
                 end do
              end if
           end do
           l=MattotP(i)%zone(j,1)+1
        end do
     end if
  end do
end subroutine resolCZM

subroutine damageUpdate_(s1,s2,s01,s02,eta,d)
  implicit none
  double precision,intent(inout)    :: d
  double precision,intent(in)       :: s1,s2,s01,s02,eta
  double precision                  :: beta,s1P,d1
  double precision,parameter        :: kk=1.
  !positive part of normal displ. <s1>+
  if (s1<0) then
     s1P = 0.
  else
     s1P = s1
  end if
  !mixed-mode damage criterion
  beta = sqrt( (s1P/s01)**2 + (s2/s02)**2 ) - 1.
  !damage evolution
  if (beta<=0) then
     return
  else
     d1 = beta/(1.+beta)/eta
     !if (d1>1) d1 = 1
     if (d1>kk) d1 = kk
     if (d1<d) d1 = d
     d = d1
  end if
end subroutine damageUpdate_

subroutine damageCV_adhoc()
! special treatment for the damage variables of composite voxles
!  --> if both d and dI > kk, then d==kk
  implicit none
  real(mytype) :: kk, dBulk, dInterphase
  integer :: i,j,k,m,id,ii,l
  do i=1,size(mattotP)
     !indice minimum de zone
     l=1
     if (mattotP(i)%Nphase>1) then   !only for composite voxels ---------------------------------------------
        do j=1,size(MattotP(i)%zone(:,1))
           do k=l,MattotP(i)%zone(j,1)
              !linear index of voxel position
              m = MattotP(i)%pos(k)
              !extract the bulk damage variable
              dBulk = damageVar(m,1)
              !extract the interphase damage variable
              id = mattotP(i)%Nphase + 9
              if (trim(mattotP(i)%LawName_comp(1))=='czm_alfanosacco_explicit') then
                 dInterphase = MattotP(i)%VarInt(id+1,k)
              else
                 do ii=2,mattotP(i)%Nphase
                    id = id + mattotP(i)%VarInt(ii-1,k)+9  !index to the 1st VarInt_i
                    if (trim(mattotP(i)%LawName_comp(ii))=='czm_alfanosacco_explicit') then
                       dInterphase = MattotP(i)%VarInt(id+1,k)
                       exit
                    end if
                 end do
              end if
              !check and treat
              if (dBulk>kk .AND. dInterphase>kk) then
                 damageVar(m,1) = kk
              end if
           end do
           l=MattotP(i)%zone(j,1)+1
        end do
     end if
  end do

end subroutine damageCV_adhoc






subroutine damageUpdate_wS(s1,s2,s01,s02,eta,dS,d)
  implicit none
  double precision,intent(inout)    :: d
  double precision,intent(in)       :: s1,s2,s01,s02,eta
  double precision,intent(in)       :: dS
  double precision                  :: beta,s1P,d1
  double precision,parameter        :: kk=1. - 1.e-6
  !positive part of normal displ. <s1>+
  if (s1<0) then
     s1P = 0.
  else
     s1P = s1
  end if
  !mixed-mode damage criterion
  beta = sqrt( (s1P/s01)**2 + (s2/s02)**2 ) - 1.
  beta = beta / dS
  !damage evolution
  if (beta<=0) then
     return
  else
     d1 = beta/(1.+beta)/eta
     if (d1>1) d1 = 1
     if (d1<d) d1 = d
     d = d1
  end if
end subroutine damageUpdate_wS


subroutine ChangeCoord_stress(sigx,e1,e2,e3,sigx_loc)
! transform from global to local coordinates for stress tensor
! ATTENTION: only useful for strain tensor --> (eps4,eps5,eps6)=2*(eps12,eps13,eps23)
!                                              this holds in both local and global coords
!  input:
!              sigx - "tensor" in the global coordinates
!         e1,e2,e3 - 3 axes describing the local coordinates
! output:
!             sigx_loc - "tensor" in the local coordinates
  implicit none
  double precision,dimension(3),intent(in)  :: e1,e2,e3
  double precision,dimension(6),intent(in)  :: sigx
  double precision,dimension(6),intent(out) :: sigx_loc
  sigx_loc(1) = ( sigx(1)*e1(1) + sigx(4)*e1(2) + sigx(5)*e1(3) ) * e1(1) &
              + ( sigx(4)*e1(1) + sigx(2)*e1(2) + sigx(6)*e1(3) ) * e1(2) &
              + ( sigx(5)*e1(1) + sigx(6)*e1(2) + sigx(3)*e1(3) ) * e1(3)
  sigx_loc(2) = ( sigx(1)*e2(1) + sigx(4)*e2(2) + sigx(5)*e2(3) ) * e2(1) &
              + ( sigx(4)*e2(1) + sigx(2)*e2(2) + sigx(6)*e2(3) ) * e2(2) &
              + ( sigx(5)*e2(1) + sigx(6)*e2(2) + sigx(3)*e2(3) ) * e2(3)
  sigx_loc(3) = ( sigx(1)*e3(1) + sigx(4)*e3(2) + sigx(5)*e3(3) ) * e3(1) &
              + ( sigx(4)*e3(1) + sigx(2)*e3(2) + sigx(6)*e3(3) ) * e3(2) &
              + ( sigx(5)*e3(1) + sigx(6)*e3(2) + sigx(3)*e3(3) ) * e3(3)
  sigx_loc(4) = ( sigx(1)*e2(1) + sigx(4)*e2(2) + sigx(5)*e2(3) ) * e1(1) &
              + ( sigx(4)*e2(1) + sigx(2)*e2(2) + sigx(6)*e2(3) ) * e1(2) &
              + ( sigx(5)*e2(1) + sigx(6)*e2(2) + sigx(3)*e2(3) ) * e1(3) 
  sigx_loc(5) = ( sigx(1)*e3(1) + sigx(4)*e3(2) + sigx(5)*e3(3) ) * e1(1) &
              + ( sigx(4)*e3(1) + sigx(2)*e3(2) + sigx(6)*e3(3) ) * e1(2) &
              + ( sigx(5)*e3(1) + sigx(6)*e3(2) + sigx(3)*e3(3) ) * e1(3)
  sigx_loc(6) = ( sigx(1)*e3(1) + sigx(4)*e3(2) + sigx(5)*e3(3) ) * e2(1) &
              + ( sigx(4)*e3(1) + sigx(2)*e3(2) + sigx(6)*e3(3) ) * e2(2) &
              + ( sigx(5)*e3(1) + sigx(6)*e3(2) + sigx(3)*e3(3) ) * e2(3) 
end subroutine ChangeCoord_stress


subroutine irreversibility_dBulk()
!!Enforce the irreversibility of bulk damage variable
  implicit none
  integer :: i,j,k,l,m, ii
  integer,parameter :: id_d=1 !need to be consistent with "mat_*.xml"
  integer :: id  !indices
  do i=1,size(mattotP)
     !indice minimum de zone
     l=1
     !no operation for interphase materials
     if (mattotP(i)%Interphase) cycle !no calculation for interphase material
     !no operation for composite voxels
     if (mattotP(i)%Nphase>1) then !for composite voxels-----------------------------------------------------
        !update the current internal variable
        MattotP(i)%VarInt=MattotP(i)%VarInt0
        do j=1,size(MattotP(i)%zone(:,1))
           !volume fraction of each phase
           do ii=1,mattotP(i)%Nphase-1
              id = 1 + mattotP(i)%Nphase + ii !index to fv_i
           end do
           do k=l,MattotP(i)%zone(j,1)
              !linear index of voxel position
              m = MattotP(i)%pos(k)
              !force the irrversibility (not applied to interphase)
              id = mattotP(i)%Nphase + 9
              if (trim(mattotP(i)%LawName_comp(1))/='czm_alfanosacco_explicit') then
                 damageVar(m,1) = maxval( (/damageVar(m,1), dble(mattotP(i)%VarInt(id+id_d,k))/) )
              end if
              do ii=2,mattotP(i)%Nphase
                 id = id + mattotP(i)%VarInt(ii-1,k) + 9
                 if (trim(mattotP(i)%LawName_comp(ii))=='czm_alfanosacco_explicit') cycle
                 damageVar(m,1) = maxval( (/damageVar(m,1), dble(mattotP(i)%VarInt(id+id_d,k))/) )
              end do
           end do
           l=MattotP(i)%zone(j,1)+1
        end do
     else                          !for homogeneous voxels---------------------------------------------------
        !update the current internal variable
        MattotP(i)%VarInt=MattotP(i)%VarInt0
        do j=1,size(MattotP(i)%zone(:,1))
           do k=l,MattotP(i)%zone(j,1)
              !linear index of voxel position
              m = MattotP(i)%pos(k)
              !force the irreversibility
              damageVar(m,1) = maxval( (/damageVar(m,1), dble(MattotP(i)%VarInt(1,k))/) )
           end do
           l=MattotP(i)%zone(j,1)+1
        end do
     end if
  end do
end subroutine irreversibility_dBulk













































!! net code ****************************************************************************************
!! *************************************************************************************************
subroutine rs ( n, a, w, matz, z, ierr )

!*****************************************************************************80
!
!! RS computes eigenvalues and eigenvectors of real symmetric matrix.
!
!  Discussion:
!
!    RS calls the recommended sequence of EISPACK routines
!    to find the eigenvalues and eigenvectors (if desired)
!    of a real symmetric matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2018
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the real symmetric matrix.
!
!    Input, logical MATZ, is false if only eigenvalues are desired, 
!    and true if both eigenvalues and eigenvectors are desired.
!
!    Output, real ( kind = 8 ) W(N), the eigenvalues in ascending order.
!
!    Output, real ( kind = 8 ) Z(N,N), contains the eigenvectors, if MATZ
!    is true.
!
!    Output, integer ( kind = 4 ) IERR, is set equal to an error
!    completion code described in the documentation for TQLRAT and TQL2.
!    The normal completion code is zero.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) fv1(n)
  real ( kind = 8 ) fv2(n)
  integer ( kind = 4 ) ierr
  logical matz
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) z(n,n)

  if ( .not. matz ) then

    call tred1 ( n, a, w, fv1, fv2 )

    call tqlrat ( n, w, fv2, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'RS - Fatal error!'
      write ( *, '(a)' ) '  Error return from TQLRAT.'
      return
    end if

  else

    call tred2 ( n, a, w, fv1, z )

    call tql2 ( n, w, fv1, z, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'RS - Fatal error!'
      write ( *, '(a)' ) '  Error return from TQL2.'
      return
    end if

  end if

  return
end subroutine rs
!!

subroutine tred1 ( n, a, d, e, e2 )

!*****************************************************************************80
!
!! TRED1 transforms a real symmetric matrix to symmetric tridiagonal form.
!
!  Discussion:
!
!    TRED1 reduces a real symmetric matrix to a symmetric
!    tridiagonal matrix using orthogonal similarity transformations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 March 2018
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Martin, Reinsch, James Wilkinson,
!    TRED1,
!    Numerische Mathematik,
!    Volume 11, pages 181-195, 1968.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input/output, real ( kind = 8 ) A(N,N), on input, contains the real
!    symmetric matrix.  Only the lower triangle of the matrix need be supplied.
!    On output, A contains information about the orthogonal transformations
!    used in the reduction in its strict lower triangle.
!    The full upper triangle of A is unaltered.
!
!    Output, real ( kind = 8 ) D(N), contains the diagonal elements of the
!    tridiagonal matrix.
!
!    Output, real ( kind = 8 ) E(N), contains the subdiagonal elements of the
!    tridiagonal matrix in its last N-1 positions.  E(1) is set to zero.
!
!    Output, real ( kind = 8 ) E2(N), contains the squares of the corresponding
!    elements of E.  E2 may coincide with E if the squares are not needed.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) e2(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) scale

  d(1:n) = a(n,1:n)

  do i = 1, n
    a(n,i) = a(i,i)
  end do

  do i = n, 1, -1

    l = i - 1
    h = 0.0D+00
!
!  Scale row.
!
    scale = sum ( abs ( d(1:l) ) )

    if ( scale == 0.0D+00 ) then

      do j = 1, l
        d(j) = a(l,j)
        a(l,j) = a(i,j)
        a(i,j) = 0.0D+00
      end do

      e(i) = 0.0D+00
      e2(i) = 0.0D+00

      cycle

    end if

    d(1:l) = d(1:l) / scale

    do k = 1, l
      h = h + d(k) * d(k)
    end do

    e2(i) = h * scale * scale
    f = d(l)
    g = - sign ( sqrt ( h ), f )
    e(i) = scale * g
    h = h - f * g
    d(l) = f - g

    if ( 1 <= l ) then
!
!  Form A * U.
!
      e(1:l) = 0.0D+00

      do j = 1, l

        f = d(j)
        g = e(j) + a(j,j) * f
        do k = j + 1, l
          g = g + a(k,j) * d(k)
          e(k) = e(k) + a(k,j) * f
        end do

        e(j) = g

      end do
!
!  Form P.
!
      f = 0.0D+00

      do j = 1, l
        e(j) = e(j) / h
        f = f + e(j) * d(j)
      end do

      h = f / ( h + h )
!
!  Form Q.
!
      e(1:l) = e(1:l) - h * d(1:l)
!
!  Form reduced A.
!
      do j = 1, l

        f = d(j)
        g = e(j)

        a(j:l,j) = a(j:l,j) - f * e(j:l) - g * d(j:l)

      end do

    end if

    do j = 1, l
      f = d(j)
      d(j) = a(l,j)
      a(l,j) = a(i,j)
      a(i,j) = f * scale
    end do

  end do

  return
end subroutine tred1
!!
subroutine tred2 ( n, a, d, e, z )

!*****************************************************************************80
!
!! TRED2 transforms a real symmetric matrix to symmetric tridiagonal form.
!
!  Discussion:
!
!    TRED2 reduces a real symmetric matrix to a
!    symmetric tridiagonal matrix using and accumulating
!    orthogonal similarity transformations.
!
!    A and Z may coincide, in which case a single storage area is used
!    for the input of A and the output of Z.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 March 2018
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Martin, Reinsch, James Wilkinson,
!    TRED2,
!    Numerische Mathematik,
!    Volume 11, pages 181-195, 1968.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the real symmetric input matrix.  Only the
!    lower triangle of the matrix need be supplied.
!
!    Output, real ( kind = 8 ) D(N), the diagonal elements of the tridiagonal
!    matrix.
!
!    Output, real ( kind = 8 ) E(N), contains the subdiagonal elements of the
!    tridiagonal matrix in E(2:N).  E(1) is set to zero.
!
!    Output, real ( kind = 8 ) Z(N,N), the orthogonal transformation matrix
!    produced in the reduction.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  real ( kind = 8 ) hh
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) scale
  real ( kind = 8 ) z(n,n)

  do i = 1, n
    z(i:n,i) = a(i:n,i)
  end do

  d(1:n) = a(n,1:n)

  do i = n, 2, -1

    l = i - 1
    h = 0.0D+00
!
!  Scale row.
!
    scale = sum ( abs ( d(1:l) ) )

    if ( scale == 0.0D+00 ) then

      e(i) = d(l)

      do j = 1, l
        d(j) = z(l,j)
        z(i,j) = 0.0D+00
        z(j,i) = 0.0D+00
      end do

      d(i) = 0.0D+00

      cycle

    end if

    d(1:l) = d(1:l) / scale

    h = h + dot_product ( d(1:l), d(1:l) )

    f = d(l)
    g = - sign ( sqrt ( h ), f )
    e(i) = scale * g
    h = h - f * g
    d(l) = f - g
!
!  Form A*U.
!
    e(1:l) = 0.0D+00

    do j = 1, l

      f = d(j)
      z(j,i) = f
      g = e(j) + z(j,j) * f

      do k = j + 1, l
        g = g + z(k,j) * d(k)
        e(k) = e(k) + z(k,j) * f
      end do

      e(j) = g

    end do
!
!  Form P.
!
    e(1:l) = e(1:l) / h

    f = dot_product ( e(1:l), d(1:l) )

    hh = 0.5D+00 * f / h
!
!  Form Q.
!
    e(1:l) = e(1:l) - hh * d(1:l)
!
!  Form reduced A.
!
    do j = 1, l

      f = d(j)
      g = e(j)

      z(j:l,j) = z(j:l,j) - f * e(j:l) - g * d(j:l)

      d(j) = z(l,j)
      z(i,j) = 0.0D+00

    end do

    d(i) = h

  end do
!
!  Accumulation of transformation matrices.
!
  do i = 2, n

!   l = i - 1
    z(n,i-1) = z(i-1,i-1)
    z(i-1,i-1) = 1.0D+00
    h = d(i)

    if ( h /= 0.0D+00 ) then

      d(1:i-1) = z(1:i-1,i) / h

      do j = 1, i - 1

        g = dot_product ( z(1:i-1,i), z(1:i-1,j) )

        do k = 1, i - 1
          z(k,j) = z(k,j) - g * d(k)
        end do

      end do

    end if

    z(1:i-1,i) = 0.0D+00

  end do

  d(1:n) = z(n,1:n)

  z(n,1:n-1) = 0.0D+00
  z(n,n) = 1.0D+00

  e(1) = 0.0D+00

  return
end subroutine tred2
!!
subroutine tqlrat ( n, d, e2, ierr )

!*****************************************************************************80
!
!! TQLRAT computes all eigenvalues of a real symmetric tridiagonal matrix.
!
!  Discussion:
!
!    TQLRAT finds the eigenvalues of a symmetric
!    tridiagonal matrix by the rational QL method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 March 2018
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    C Reinsch,
!    Algorithm 464, TQLRAT,
!    Communications of the ACM,
!    Volume 16, page 689, 1973.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N).  On input, D contains the diagonal
!    elements of the matrix.  On output, D contains the eigenvalues in ascending
!    order.  If an error exit was made, then the eigenvalues are correct
!    in positions 1 through IERR-1, but may not be the smallest eigenvalues.
!
!    Input/output, real ( kind = 8 ) E2(N), contains in positions 2 through N 
!    the squares of the subdiagonal elements of the matrix.  E2(1) is
!    arbitrary.  On output, E2 has been overwritten by workspace
!    information.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, for no error,
!    J, if the J-th eigenvalue could not be determined after 30 iterations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e2(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) its
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( kind = 8 ) p
  !real ( kind = 8 ) pythag
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t

  ierr = 0

  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e2(i-1) = e2(i)
  end do

  f = 0.0D+00
  t = 0.0D+00
  e2(n) = 0.0D+00

  do l = 1, n

    its = 0
    h = abs ( d(l) ) + sqrt ( e2(l) )

    if ( t <= h ) then

      t = h
      b = abs ( t ) * epsilon ( b )
      c = b * b

    end if
!
!  Look for small squared sub-diagonal element.
!
    do m = l, n
      if ( e2(m) <= c ) then
        exit
      end if
    end do

    if ( m /= l ) then

      do

        if ( 30 <= its ) then
          ierr = l
          return
        end if

        its = its + 1
!
!  Form shift.
!
        l1 = l + 1
        s = sqrt ( e2(l) )
        g = d(l)
        p = ( d(l1) - g ) / ( 2.0D+00 * s )
        call pythag( p, 1.0D+00, r )
        d(l) = s / ( p + sign ( r, p ) )
        h = g - d(l)
        d(l1:n) = d(l1:n) - h
        f = f + h
!
!  Rational QL transformation.
!
        g = d(m)
        if ( g == 0.0D+00 ) then
          g = b
        end if

        h = g
        s = 0.0D+00
        mml = m - l

        do i = m - 1, l, -1
          p = g * h
          r = p + e2(i)
          e2(i+1) = s * r
          s = e2(i) / r
          d(i+1) = h + s * ( h + d(i) )
          g = d(i) - e2(i) / g
          if ( g == 0.0D+00 ) then
            g = b
          end if
          h = g * p / r
        end do

        e2(l) = s * g
        d(l) = h
!
!  Guard against underflow in convergence test.
!
        if ( h == 0.0D+00 ) then
          exit
        end if

        if ( abs ( e2(l) ) <= abs ( c / h ) ) then
          exit
        end if

        e2(l) = h * e2(l)

        if ( e2(l) == 0.0D+00 ) then
          exit
        end if

      end do

    end if

    p = d(l) + f
!
!  Order the eigenvalues.
!
    do i = l, 1, -1
      if ( i == 1 ) then
        d(i) = p
        exit
      else if ( d(i-1) <= p ) then
        d(i) = p
        exit
      end if
      d(i) = d(i-1)
    end do

  end do

  return
end subroutine tqlrat
!!

subroutine tql2 ( n, d, e, z, ierr )

!*****************************************************************************80
!
!! TQL2 computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
!
!  Discussion:
!
!    TQL2 finds the eigenvalues and eigenvectors of a symmetric
!    tridiagonal matrix by the QL method.  The eigenvectors of a full
!    symmetric matrix can also be found if TRED2 has been used to reduce this
!    full matrix to tridiagonal form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 March 2018
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Bowdler, Martin, Reinsch, James Wilkinson,
!    TQL2,
!    Numerische Mathematik,
!    Volume 11, pages 293-306, 1968.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N).  On input, the diagonal elements of
!    the matrix.  On output, the eigenvalues in ascending order.  If an error
!    exit is made, the eigenvalues are correct but unordered for indices
!    1,2,...,IERR-1.
!
!    Input/output, real ( kind = 8 ) E(N).  On input, E(2:N) contains the
!    subdiagonal elements of the input matrix, and E(1) is arbitrary.
!    On output, E has been destroyed.
!
!    Input, real ( kind = 8 ) Z(N,N).  On input, the transformation matrix
!    produced in the reduction by TRED2, if performed.  If the eigenvectors of
!    the tridiagonal matrix are desired, Z must contain the identity matrix.
!    On output, Z contains the orthonormal eigenvectors of the symmetric
!    tridiagonal (or full) matrix.  If an error exit is made, Z contains
!    the eigenvectors associated with the stored eigenvalues.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, normal return,
!    J, if the J-th eigenvalue has not been determined after
!    30 iterations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) dl1
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) el1
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) its
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( kind = 8 ) p
  !real ( kind = 8 ) pythag
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) s2
  real ( kind = 8 ) t
  real ( kind = 8 ) tst1
  real ( kind = 8 ) tst2
  real ( kind = 8 ) z(n,n)

  ierr = 0

  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e(i-1) = e(i)
  end do

  f = 0.0D+00
  tst1 = 0.0D+00
  e(n) = 0.0D+00

  do l = 1, n

    its = 0
    h = abs ( d(l) ) + abs ( e(l) )
    tst1 = max ( tst1, h )
!
!  Look for a small sub-diagonal element.
!
    do m = l, n
      tst2 = tst1 + abs ( e(m) )
      if ( tst2 == tst1 ) then
        exit
      end if
    end do

    if ( m /= l ) then

      do

        if ( 30 <= its ) then
          ierr = l
          return
        end if

        its = its + 1
!
!  Form shift.
!
        l1 = l + 1
        l2 = l1 + 1
        g = d(l)
        p = ( d(l1) - g ) / ( 2.0D+00 * e(l) )
        call pythag( p, 1.0D+00 , r )
        d(l) = e(l) / ( p + sign ( r, p ) )
        d(l1) = e(l) * ( p + sign ( r, p ) )
        dl1 = d(l1)
        h = g - d(l)
        d(l2:n) = d(l2:n) - h
        f = f + h
!
!  QL transformation.
!
        p = d(m)
        c = 1.0D+00
        c2 = c
        el1 = e(l1)
        s = 0.0D+00
        mml = m - l

        do i = m - 1, l, -1

          c3 = c2
          c2 = c
          s2 = s
          g = c * e(i)
          h = c * p
          call pythag( p, e(i) , r )
          e(i+1) = s * r
          s = e(i) / r
          c = p / r
          p = c * d(i) - s * g
          d(i+1) = h + s * ( c * g + s * d(i) )
!
!  Form vector.
!
          do k = 1, n
            h = z(k,i+1)
            z(k,i+1) = s * z(k,i) + c * h
            z(k,i) = c * z(k,i) - s * h
          end do

        end do

        p = - s * s2 * c3 * el1 * e(l) / dl1
        e(l) = s * p
        d(l) = c * p
        tst2 = tst1 + abs ( e(l) )

        if ( tst2 <= tst1 ) then
          exit
        end if

      end do

    end if

    d(l) = d(l) + f

  end do
!
!  Order eigenvalues and eigenvectors.
!
  do i = 1, n - 1

    k = i
    p = d(i)

    do j = i + 1, n

      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if

    end do

    if ( k /= i ) then

      d(k) = d(i)
      d(i) = p

      do j = 1, n
        t      = z(j,i)
        z(j,i) = z(j,k)
        z(j,k) = t
      end do

    end if

  end do

  return
end subroutine tql2
!!

subroutine pythag ( a, b, pythag1 )

!*****************************************************************************80
!
!! PYTHAG computes SQRT ( A * A + B * B ) carefully.
!
!  Discussion:
!
!    The formula
!
!      PYTHAG1 = sqrt ( A * A + B * B )
!
!    is reasonably accurate, but can fail if, for example, A^2 is larger
!    than the machine overflow.  The formula can lose most of its accuracy
!    if the sum of the squares is very large or very small.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the two legs of a right triangle.
!
!    Output, real ( kind = 8 ) PYTHAG1, the length of the hypotenuse.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) p
  real ( kind = 8 ) pythag1
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) u

  p = max ( abs ( a ), abs ( b ) )

  if ( p /= 0.0D+00 ) then

    r = ( min ( abs ( a ), abs ( b ) ) / p )**2

    do

      t = 4.0D+00 + r

      if ( t == 4.0D+00 ) then
        exit
      end if

      s = r / t
      u = 1.0D+00 + 2.0D+00 * s
      p = u * p
      r = ( s / u )**2 * r

    end do

  end if

  pythag1 = p

  return
end subroutine pythag

! ====================================================================================================
!                                fin des subroutines for resolutionPF_YC
! ====================================================================================================
! ****************************************************************************************************
! ****************************************************************************************************
! ****************************************************************************************************
! ****************************************************************************************************
! ****************************************************************************************************
! ****************************************************************************************************
! ****************************************************************************************************
! ****************************************************************************************************
! ****************************************************************************************************
! ****************************************************************************************************
! ****************************************************************************************************





end module amitex_user_mod
