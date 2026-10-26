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
#ifdef NEW_2DECOMP
  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
#else
  use decomp_2d
#endif             
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
  use amitex_mpi_mod, only: amitex_mpi_comm
  
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

! Declare global procedures as public (used in resolution_user_mod or standard_user_mod)
  public :: initDamageHH,resolPF_visc,getHistoryTerm_visc,updateHistory, apply_GreenPF_visc
  public :: updateIntVar
  public :: rs,tred1,tred2,tqlrat,tql2,pythag 


!!------------------------------------------------------------------------------
!>                                                                 PUBLIC FIELDS 

  !> real space (ntot,ncomp)
  real(mytype),allocatable,dimension(:,:)        :: damageVar, tau, tau0, HH, AA, BB
  real(mytype),allocatable,dimension(:,:,:)      :: ACT3_Rdamage, ACT3_Udamage

  !> Fourier space ((nx,ny,nz,ncomp)
  complex(mytype),allocatable,dimension(:,:,:,:) :: tauF
  real(mytype),allocatable,dimension(:,:,:)      :: FreqLaplacian


!!------------------------------------------------------------------------------
!>                                                                       indices
  integer :: NdamageVar


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

  ! initialise the freqency for laplacian operator
  call initFreqLaplacian(FreqLaplacian)

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
subroutine initDamageHH()
!!Initialize the damage and history fields
  implicit none
  integer(INT64) :: i,j,k,l,m
  !initialize the damage field by internal variable
  !------------------------------------------------
  do i=1,size(mattotP)
     !indice minimum de zone
     l=1
     !update the current internal variable
     MattotP(i)%VarInt=MattotP(i)%VarInt0
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
  end do
end subroutine initDamageHH

subroutine resolPF_visc(dt,nitPF,critPF,timePF)
!! solve the phase-field problem
!  param[in]
!       - ind_tps: index of loading step
!       -     dt : time increment between t-1 and t (current)
!       -   visc : viscous parameter
!       -     lc : characteristic length
!       -niterMaxPF: max number of iteration of phase-field solution
!  param[inout] global variables in amitex_mod
!       -     HH : history field
!       -    eta : damage variable
!  param[out]
!       -  nitPF : number of iterations in the phase-field resolution
!       - critPF : average value of the residual of local phase-field PDE
!       - timePF : time duration used for solving the phase-field problem
  implicit none

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
 

  critPF=0.
  tol_critPF = user_param%p_real(1)
  nItPF=0
  testCV_PF  = .true.
  timePF = 0.

  iact3=0
  lact3=.false.
  acv_test = .true.

  acc_CV_PF = user_param%p_string(1)=="true"
  niterMaxPF = int(user_param%p_real(2))


  !! compute the terms associated to history field
  !! =============================================
  call getHistoryTerm_visc(dt)  !AA=1/lc**2+visc/dt/lc/gc+2*HH/lc/gc
                                !BB=visc/dt/lc/gc*eta0+2*HH/lc/gc

  ! "reference phase" AA0 : different choices
  ! -----------------------------------------
  call MPI_AllReduce(minval(AA),AAmin,1,real_type,MPI_MIN,amitex_mpi_comm,ierr) !min of AA
  call MPI_AllReduce(maxval(AA),AAmax,1,real_type,MPI_MAX,amitex_mpi_comm,ierr) !max of AA
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

     !! time counter
     !! ============
     call mpi_barrier(amitex_mpi_comm,ierr)
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
!       -   BB : viscous term, BB = 2*HH/lc/gc + visc/dt/lc/gc*eta0
  implicit none
  real(mytype),intent(in) :: dt
  real(mytype)            :: lc, visc, gc, visco
  integer(INT64)          :: i,j,k     !< indice de boucle 
  integer(INT64)          :: m,l       !index of voxel, minimum number of zone
  !                                                                   Loop of every voxel
  ! -------------------------------------------------------------------------------------
  do i=1,size(mattotP)
     !indice minimum de zone
     l=1
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
           !compute A(x) = 1/lc/lc + visc/dt/lc/gc + 2*HH/lc/gc
           AA(m,1) = 1._mytype/lc/lc + visco + 2._mytype*HH(m,1)/lc/gc
           !compute B(x) = 2*HH/lc/gc + visc/dt/lc/gc*eta0
           BB(m,1) = 2._mytype*HH(m,1)/lc/gc + visco*damageVar(m,1)
        end do
        l=MattotP(i)%zone(j,1)+1
     end do
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
   integer(INT64)   :: i,j,h
   do i=fft_start(3),fft_end(3)
      do j=fft_start(2),fft_end(2)
         do h=fft_start(1),fft_end(1)
            !tauF(h,j,i,1) = tauF(h,j,i,1) / ( AA0 + (FREQ(h,j,i,1)**2+FREQ(h,j,i,2)**2+FREQ(h,j,i,3)**2) )
            tauF(h,j,i,1) = tauF(h,j,i,1) / ( AA0 - FreqLaplacian(h,j,i) )
         end do
      end do
   end do
end subroutine apply_GreenPF_visc

subroutine updateIntVar()
!!update the internal variables and assign the damage field into one of the internal variables
  implicit none
  integer(INT64) :: i,j,k,l,m
  do i=1,size(mattotP)
     !indice minimum de zone
     l=1
     do j=1,size(MattotP(i)%zone(:,1))
        do k=l,MattotP(i)%zone(j,1)
           !linear index of voxel position
           m = MattotP(i)%pos(k)
           !update internal variable by damage variable (eta->VarInt(1))
           MattotP(i)%VarInt(1,k) = damageVar(m,1)
           !update internal variable at t-1 (VarInt->VarInt0)
           MattotP(i)%VarInt0(:,k) = MattotP(i)%VarInt(:,k)
        end do
        l=MattotP(i)%zone(j,1)+1
     end do
  end do
end subroutine updateIntVar

subroutine updateHistory()
!!update the internal variables and assign the history field into one of the internal variables
  implicit none
  integer(INT64) :: i,j,k,l,m
  do i=1,size(mattotP)
     !indice minimum de zone
     l=1
     do j=1,size(MattotP(i)%zone(:,1))
        do k=l,MattotP(i)%zone(j,1)
           !linear index of voxel position
           m = MattotP(i)%pos(k)
           !update the history field by internal variable (VarInt(2)-->HH)
           if (HH(m,1)<MattotP(i)%VarInt(2,k)) then
              HH(m,1) = MattotP(i)%VarInt(2,k)
           end if
        end do
        l=MattotP(i)%zone(j,1)+1
     end do
  end do
end subroutine updateHistory

















































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
  c = 0
  b = 0
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
  s2 = 0
  c3 = 0
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
