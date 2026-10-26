!123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!===================================================================================================
!
!       MODULE RESOLUTION_MOD : 
!>
!> Resolution du probleme avec une methode FFT
!!
!!  
!!  Subroutines
!!
!! - resolution_NL_base : loop over time steps
!!                        call unpas_NL_base
!!                        but also (for user defined 'couplings') : 
!!                               - Nloc_call
!!                               OR
!!                               - before(after)_unpas_user
!!
!!
!!
!
module resolution_mod

  use ISO_FORTRAN_ENV

  use MPI             
#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype  ! Type
  use decomp_2d_mpi,       only : nrank   ! Variable
#else
  use decomp_2d,          only : mytype,& ! Type
                                 nrank    ! Variable
#endif  

  use io_amitex_mod,      only : write_stdout0 ! function
  use io2_amitex_mod,     only : times_io,& ! pointer variable
                                 writeVTK   ! function
  use algo_functions_mod, only : sigToPK1   ! function
  use param_algo_mod,     only : algo_param, user_param ! pointer variable
  use loading_mod,        only : LOCALLOAD,& !Type
                                 local_load, local_loadD, load, extract,&  ! Pointer Variables
                                 get_current_loading, get_current_loadingD ! Functions
  use material_mod,       only : mattotP, times_m, VCinfo ! Pointer Variables
  use green_mod,          only : times_g                  ! Pointer Variables
  use error_mod,          only : amitex_abort, check_amitex_abort ! Functions
  use amitex_mod,         only : Def, Def0, Sig, Sig0, PK1, FluxD, FluxD0, CVFORsauv, Flog ! Pointer Variables
  use field_mod,          only : times_f, grid,& ! Pointer Variables
                                 field_Mean ! Functions
  use non_local_mod,      only : nloc_call ! Functions
  use sortie_std_mod,     only : times_s,&   ! Pointer Variables
                                 sortie_std  ! Functions
  use NL_base_mod,        only : times_b, info_NL_b,& ! Pointer Variables
                                 unpas_NL_base, log_after_unpas, log_final_times ! Functions
  use standard_user_mod,  only : before_unpas_user, after_unpas_user ! Functions
  use user_functions_mod, only : user_load_interruption_average, user_load_interruption_internal_variables ! Functions
  use simu_mod,           only : simu_restart_write, simu_restart_read
  use amitex_mpi_mod, only: amitex_mpi_comm

  private
   
  !> Pointer "publiques" vers composantes de SIMU_AMITEX
  ! nothing

  !> Variables "publiques" utilisees lors de l'initialisation
  ! nothing

  !> Types publiques (pour definition de SIMU_AMITEX)
  ! nothing

  !> Fonctions publiques
  public :: resolution_NL_base

  
contains

!===========================================================================================================
!!==========================================================================================================
!> resolution_NL_base : Schema de base pour resoudre un probleme non lineaire
!!                       en utilisant des transformees de Fourier rapides
!!                       Cette procédure contient l'algorithme de resolution
!!                       la procédure traite la mécanique et la diffusion
!!
!! NOTATIONS (entree xml et champs) Sigma et Def (HPP) : 
!!                 CAST3M (Voigt + ordre 11 22 33 12 13 23)
!!                 PK1 et grad(u) (GDEF) : ordre 11 22 33 12 13 23 21 31 32 sans notation
!!
!!                             Variables dans MATERIAL_MOD
!!--------------------------------------------------------
!!   LambdaMu0   Coefficients [Lambda0,mu0] du milieu de référence
!!   C0          matrice de rigidite du materiau de reference
!!   k0D         coefficient de conductivité (diffusion) 
!!               taille (nVarD)
!!
!!   nmateriaux  nombre de matériaux
!!   
!!   VCinfo      info sur algorithme laminate
!!
!!           CHAMPS paralleles accessibles dans AMITEX_MOD
!!--------------------------------------------------------
!! Descrition des champs dans l'espace de Fourier :
!!             pinceaux Z, (nx/2+1,ny,nz,3,nVarD ou NTensDef) 
!!             indices 3D : (fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),xx)
!!
!! Descrition des champs dans l'espace Reel :
!!             pinceaux X, (nx,ny,nz,3,nVarD) ou (nx,ny,nz,NTensDef) 
!!             indice 1D : (1:amx_size(1)*amx_size(2)*amx_size(3),NTensDef ou 3,nVarD)
!!
!!   Sig      Champ de contrainte de Cauchy
!!   PK1      Champ de contrainte de Piola-Kirchhoff en grandes transformations
!!   Def      Champ de deformation en HPP et grad(u) en GD
!!   Def0     Champ de deformation en HPP et grad(u) en GD a l'instant precedent
!!   Sig0     Champ de contrainte de Cauchy a l'instant precedent
!!   SigF     Contrainte de Cauchy en HPP ou contrainte de Piola-Kirchhoff en GD,
!!            espace Fourier 
!!   DefF     Deformation en HPP ou grad(u) en GD, espace Fourier 
!!            pinceaux-Z (nx/2+1,ny,nz,nTensDef) 
!!   FluxD0   Champ de vecteurs flux (:,3,nVarD) (diffusion)
!!   FluxD    Champ de vecteurs flux (:,3,nVarD) (diffusion) 
!!   GradQD0  Champ de gradient de variable de diffusion (:,3,nVarD) (diffusion)
!!   GradQD   Champ de gradient de variable de diffusion (:,3,nVarD) (diffusion) 
!!   FluxDF   Champ de vecteur flux dans Fourier (diffudion)
!!   GradQDF  Champ de gradient dans Fourier 
!!
!! Tableaux pour l'acceleration de convergence
!!                 pour la mecanique (ntot,Ntens,4)
!!                 pour la diffusion (ntot,3*nVarD,4)
!!   ACT3_R   Residus pour l'acceleration de convergence
!!   ACT3_U   Increment des champs de deformations pour l'acceleration de convergence
!!   ACT3_RD  Idem ACT3_R pour diffusion
!!   ACT3_UD  Idem ACT3_U pour la diffusion
!!
!!                        autres variables dans AMITEX_MOD
!!--------------------------------------------------------
!!   Flog	Unites logiques pour l'ouverture du fichier log
!!   fic_vtk    Racine du nom des fichiers de sortie vtk ou (mz)std
!!
!!                                Variables dans GREEN_MOD
!!--------------------------------------------------------
!!   FREQ     Tableau des fréquences, espace Fourier 
!!            pinceaux-Z (nx/2+1,ny,nz,3)
!!            indices 3D : (fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),3)
!!
!!                              Variables dans LOADING_MOD
!!--------------------------------------------------------
!!   local_load%t1         chargement mecanique impose a chaque iteration
!!   local_loadD%t1        chargement "diffusion" impose a chaque iteration
!!
!!
!!                              Variables dans NL_BASE_MOD
!!--------------------------------------------------------
!!   crit_b                criteria for resolution_NL_base
!!   times_b               elapsed times for resolution_NL_base
!!
!!------------------------------------------------------------
!!                             AUTRES VARIABLES DE TYPE-DERIVE
!!
!!      MattotP         tableau de structure MATERIAL   (material_mod.f90)
!!            -> description du materiau
!!      
!!      load            tableau de structure LOADING    (loading_mod.f90))
!!            -> description du chargement
!!
!!      initValExt      structure INITLOADEXT           (loading_mod.f90)
!!            -> description de l'initialisation du chargement externe
!!
!!      extract         structure PARAMEXTRACT          (loading_mod.f90)
!!            -> description des sorties vtk
!!
!!      algo_param      structure PARAM_ALGO            (param_algo_mod.f90)
!!            -> description des parametres algorithmiques
!!
!!      grid            structure GRID_DIM              (green_mod.f90)
!!            -> description de la grille 
!!       
!!      VCinfo          structure VOXCOMP_INFO          (material_mod.f90)
!!            -> info de suivi des algo voxels composites
!!
!!===============================================================================
subroutine resolution_NL_base()

!!==============================================================================
!!                                                                  DECLARATIONS
!!==============================================================================

  implicit none
!!------------------------------------------------------------------------------
!>                                                     PARAMETRES "PARALLELISME"

  integer :: ierror                                      !< erreur relative au fonction MPI
  
!!------------------------------------------------------------------------------
!>                                                       PARAMETRES "CHARGEMENT"

  integer                                     :: ind_tps !< indice du pas de temps courant
  integer                                     :: load_n, load_incr
                                                         !< numero du chargement partiel 
                                                         !< numero de l'increment du chargement

  logical      :: test_load_interruption  !< test interruption du chargement courant
  real(mytype) :: alpha_interruption      !< parametre pour ajuster le chargement lors d'une interruption
                                     
!!------------------------------------------------------------------------------
!>                                           PARAMETRES DE SUIVI DE L'ALGORITHME

  integer               :: nIt                    !< nombres d'iterations dans un pas de chargement
  integer               :: nIt0                   !< nombre d'iterations dans le pas de temps "fictif"
  integer               :: nitTot                 !< nombres d'iterations total

!!------------------------------------------------------------------------------
!>                                  PARAMETRES DE SUIVI DE L'ALGORITHME LAMINATE

  logical                       :: testLaminate,testLaminateloc !< vrai si présence de voxels composites laminate 
                                                                !< uniquement utile pour sortir es infos laminate

  !< suivi performance algorithme laminate (material_mod)
  !<               VCinfo%nSub_lam         : nbre de recours à la subdivision du pas de temps 
  !<               VCinfo%nIncr_lam        : nbre de sous pas de temps par subdivision 
  !<               VCinfo%nIt_lam          : nbre d'itérations par appel de l'algorithme laminate
  !<               VCinfo%nCall_lam        : nbre d'appel à l'algorithme laminate
  !<               VCinfo%nSub_lam_pas          
  !<               VCinfo%nIncr_lam_pas        
  !<               VCinfo%nIt_lam_pas          
  !<               VCinfo%nCall_lam_pas        


!!------------------------------------------------------------------------------
!>                                            PARAMETRES DE L'ALGO DE RESOLUTION

  real(mytype),dimension(algo_param%nTensDef) :: defMoy, sigMoy      !<  contrainte, def... moyennes fin d'increment
  real(mytype),dimension(algo_param%nTensDef) :: defMoy0, sigMoy0    !<  contrainte, def... moyennes debut d'increment 
                                                                     !   utile dans user_load_interruption_average
  real(mytype),dimension(algo_param%nTensDef) :: defMoy00, sigMoy00  !<  contrainte, def... moyennes debut de chargement partiel

  real(mytype),dimension(3,algo_param%nVarD)  :: gradQDMoy, FluxDMoy     !<  grad et flux moyens fin d'increment
  real(mytype),dimension(3,algo_param%nVarD)  :: gradQDMoy0, FluxDMoy0   !<  grad et flux moyens debut d'increment
  real(mytype),dimension(3,algo_param%nVarD)  :: gradQDMoy00, FluxDMoy00 !<  grad et flux moyens debut de chargement partiel

!!------------------------------------------------------------------------------
!>                                     PARAMETRES DE SUIVI DES TEMPS D'EXECUTION

  !! type double precision au lieu de real_mytype (=> requis par MPI_WTIME)
  double precision :: t1
  double precision :: t_fft0, t_ifft0, t_behavior0, t_green0, t_crit0, t_wvtk0,&
                      t_wstd0, t_unpas0
  !< t1          : date de repere avant les fonctions
  !! t_xxx0      : temps cumules en debut de pas
  
!!------------------------------------------------------------------------------
!>                                                                        DIVERS

  integer               :: i,j                      !< indice de boucle 
  logical               :: tmp_bool
  real(mytype)          :: tmp_real
  integer,dimension(2)  :: tmp_int

!! Variables pour test restart
!  character(len=250)   :: totofile
!  integer              :: totofid

!!==============================================================================
!!                                                               INITIALISATIONS
!!==============================================================================
  call MPI_BARRIER(amitex_mpi_comm,ierror)
  times_b%startalgo = MPI_WTIME()

  nIttot= 0
  nIt0 = 0

  if (algo_param%Mechanics) then
    defMoy=0                 !< contraintes et deformations moyennes
    sigMoy=0
    defMoy0 = 0
    sigMoy0 = 0
    defMoy00 = 0
    sigMoy00 = 0
  end if 
  if (algo_param%Diffusion) then  
    GradQDMoy = 0
    FluxDMoy = 0
    GradQDMoy0 = 0
    FluxDMoy0 = 0
    GradQDMoy00 = 0
    FluxDMoy00 = 0
  end if
  ind_Tps= 1               !< indice du premier pas

  t1 = 0                   !< Temps d'execution

  testlaminateloc = .false.  !< test de presence materiau "laminate"
  do i=1,size(mattotP)
   if (trim(mattotP(i)%lawname) == "laminate") testlaminateloc = .true.
  end do
  call MPI_Allreduce(testlaminateloc, testlaminate, 1, MPI_LOGICAL, MPI_LOR, amitex_mpi_comm, ierror)    

  !***************************************************************************************************
  !                                                          DEBUT BOUCLE SUR LES CHARGEMENTS PARTIELS
  BOUCLE_CHARGEMENTS:&
  do load_n=1,size(load)
     test_load_interruption = .false.
     alpha_interruption = 1.
     !************************************************************************************************
     !                                             DEBUT BOUCLE SUR INCREMENTS D'UN CHARGEMENT PARTIEL
     BOUCLE_INCREMENTS:&
     do load_incr = 0,load(load_n)%NIncr
     
        !================================================================================
        !                                                                 RESTART (BEGIN)

        if (algo_param%Mechanics) then

          ! RESTART : BEFORE STARTING POINT
          if (load(load_n)%skip_before_restart .or. load_incr <= load(load_n)%restart_inc) then
            ! READ averages at the begining of the partial loading used for restart
            if (.not. load(load_n)%skip_before_restart) then
               call simu_restart_read(tmp_int(1),tmp_int(2),0,defMoy00,sigMoy00, gradQDMoy00, FluxDMoy00)
            end if
            ! MAJ CHARGEMENT local_load(D)
            if(load_incr == 0) ind_tps=ind_tps-1  ! on n'incremente pas ind_tps pour le pas de temps fictif
            call update_local_load(load_n,load_incr,defMoy00, sigMoy00, gradQDMoy00, FluxDMoy00)
            if (nrank==0) write(OUTPUT_UNIT,*) "SKIP before RESTART, load_n, load_incr :", load_n, load_incr
            ind_tps=ind_tps+1
            ! CYCLE
            cycle
          end if

          ! RESTART : STARTING POINT
          if (load_incr == load(load_n)%restart_inc+1) then
            call field_Mean(SIG0,grid%ntot,6,sigMoy0)
            call field_Mean(DEF0,grid%ntot,algo_param%nTensDef,DefMoy0)
          end if
        end if
        
        !================================================================================
        !                                                     REFERENCES MESURES DE TEMPS
        !

        t_fft0      = times_f%fft
        t_ifft0     = times_f%ifft
        t_behavior0 = times_m%behavior
        t_green0    = times_g%apply
        t_crit0     = times_g%crit
        t_wvtk0     = times_io%wvtk
        t_wstd0     = times_s%wtot
        t_unpas0    = times_b%unpas
        
        call MPI_BARRIER(amitex_mpi_comm,ierror) 
        times_b%startstep = MPI_WTIME()
 
        !================================================================================
        !                                                    MAJ CHARGEMENT local_load(D)
        !
        !                             actualisation de local_load(D)%t0, t1, dt et t_load
        !

        if(load_incr == 0) ind_tps=ind_tps-1  ! on n'incremente pas ind_tps pour le pas de temps fictif
        
        call update_local_load(load_n,load_incr,defMoy00, sigMoy00, gradQDMoy00, FluxDMoy00)

        !================================================================================
        !                                                             APPEL UNPAS_NL_BASE
        !
        !
 
        if (user_param%test .AND. trim(user_param%algo) == "standard") then
           call before_unpas_user(load_n,load_incr,ind_tps)
        end if

        if (algo_param%Nloc) then !---- Mecanique non locale resolution explicite
           call Nloc_call(load_n,load_incr,ind_tps,"befor",algo_param%nloc_explicit)
        end if
        
!print *,"AV UNPAS_NL_BASE", nrank
        call unpas_NL_base(load_n,load_incr,ind_tps,nIt)
        
        if (algo_param%substepping .and. info_NL_b%reach_nitermax) then 
           call substepping(load_n,load_incr,ind_tps)
        end if

        if (user_param%test .AND. trim(user_param%algo) == "standard") then
           call after_unpas_user(load_n,load_incr,ind_tps)
        end if

        call MPI_BARRIER(amitex_mpi_comm,ierror) !necessaire pour le calcul non local ?

        if (algo_param%Nloc) then !---- Mecanique non locale resolution explicite
!print *,"AV NLOC_CALL", nrank
           call Nloc_call(load_n,load_incr,ind_tps,"after",algo_param%nloc_explicit)
        end if
        
        !================================================================================
        !                                                       ABORT IF nitermax_reached

        if (info_NL_b%reach_nitermax) then
           call amitex_abort("Max number of iteration (Nitermax) reached, stop simulation.",2,0)
        end if

!print *,"AV USER LOAD INTERRUPTION", nrank
        !================================================================================
        !                                                          USER LOAD INTERRUPTION 

        !--TODO : prevoir de recuperer les moyennes depuis unpas_NL_base

        if (algo_param%Mechanics .AND. load(load_n)%user_interruption_flag) then
        if(algo_param%HPP) then
           call field_Mean(SIG,grid%ntot,algo_param%nTensDef,sigMoy)
        else
           call field_Mean(PK1,grid%ntot,algo_param%nTensDef,sigMoy)
        end if
        call field_Mean(Def,grid%ntot,algo_param%nTensDef,defMoy)
         
        !--- Critere d'arret sur contrainte-deformation moyenne
        call user_load_interruption_average(&
                          tmp_bool,tmp_real,&
                          load_n,sigMoy0,sigMoy,defMoy0,defMoy,&
                          local_load%t0,local_load%t1,local_load%t_load,&
                          load(load_n)%user_interruption)
        test_load_interruption = test_load_interruption .OR. tmp_bool
        alpha_interruption = min(alpha_interruption, tmp_real)


        !--- Critere d'arret sur variable interne
        do i=1,size(MattotP)
             call user_load_interruption_internal_variables(&
                          tmp_bool,tmp_real,&
                          load_n,MattotP(i)%VarInt0,MattotP(i)%VarInt,MattotP(i)%numM,&
                          load(load_n)%user_interruption) 
             test_load_interruption = test_load_interruption .OR. tmp_bool
             alpha_interruption = min(alpha_interruption, tmp_real)
        end do
        end if
        
!print *,"AV INTERRUPTION ALPHA_INTERRUPTION", nrank
        !================================================================================
        !                                                 INTERRUPTION ALPHA_INTERRUPTION
        ! 
        !         on ajuste SIG, DEF, mattot%Varint, local_load%t1 et local_load%t_load(0) 
        !
        !                    + UPDATE load()%time : 
        !                    chargement courant : on fige le temps des increments suivants
        !                    chargements suivants : on decale les temps
        !
        ! TODO : a revoir ou a limiter a une interruption sans reprise
 
        if (algo_param%Mechanics .AND. load(load_n)%user_interruption_flag) then       
        if (test_load_interruption .and. load_incr .NE. 0 .and. alpha_interruption .ne. 1.) then
           if (alpha_interruption > 1.) then
             call amitex_abort("Load Interruption (in user_functions_mod) with alpha > 1 ",2)
           else
             Sig  = Sig0 + alpha_interruption * (Sig - Sig0)
             Def  = Def0 + alpha_interruption * (Def - Def0)
             if (.not. algo_param%HPP) call SigToPK1(Sig,Def,PK1)
             do i=1,size(MattotP)
               MattotP(i)%VarInt = MattotP(i)%VarInt0 &
                                 + alpha_interruption*(MattotP(i)%VarInt-MattotP(i)%VarInt0) 
             end do
             local_load%t1(algo_param%nTensDef:) = local_load%t0(algo_param%nTensDef:) +&
               alpha_interruption*(local_load%t1(algo_param%nTensDef:)-local_load%t0(algo_param%nTensDef:))
             local_load%t_load(0) = local_load%t_load(-1) +&
                                    alpha_interruption*(local_load%t_load(0) - local_load%t_load(-1))
             
             ! Decalage du temps des chargements suivants
             do i=load_n+1,size(load)
                load(i)%time=load(i)%time - (load(load_n)%time(load(load_n)%Nincr) - local_load%t_load(0))
             end do 
             ! On fige le reste du chargement courant (time)
             load(load_n)%time(load_incr:) = load(load_n)%time(load_incr)
 
           end if
        end if  
        end if

!print *,"AV FICHIER DE SORTIE STANDARD", nrank
        !================================================================================
        !                                                      FICHIER DE SORTIE STANDARD
        !
        ! contraintes et deformations moyennes, ecarts types
        ! +verif valeurs moyenne sur la cellule
        !
        ! \TODO Dissocier sortie standard et evaluation des moyennes 
        !

        !> pour un premier pas de chargement, on tient compte des itérations du pas de temps fictif 
        if(load_incr .eq. 0) nIt0 = nIt
        if(load_incr .eq. 1) nIt = nIt + nIt0

        if(load_incr/=0) then  ! Pas de sortie std pour les pas de temps fictifs
            if (algo_param%Mechanics) tmp_real = local_load%t_load(0)
            if (algo_param%Diffusion) tmp_real = local_loadD%t_load(0)
            call sortie_std(tmp_real,ind_tps,nIt,&
                            sigMoy,defMoy,FluxDMoy, GradQDMoy)
        end if

        !> on retire nIt0 ajoute plus haut pour un premier pas de chargement
        if(load_incr .eq. 1) nIt = nIt - nIt0 
         
!print *,"AV FICHIER DE SORTIE VTK", nrank
        !================================================================================
        !                                                           FICHIER DE SORTIE VTK
        
        if(load_incr/=0) then ! Pas de sortie vtk pour les pas de temps fictifs
           ! sortie sur demande ou sur interruption
           if(extract%tpsVTK(ind_tps) .OR. test_load_interruption)then 

              call writeVTK(ind_tps)

           end if
        end if
        
!print *,"AV ACTUALISATION DE LA CONTRAINTE", nrank
        !================================================================================
        !                ACTUALISATION DE LA CONTRAINTE DE CAUCHY, DES VARIABLES INTERNES
        !
        
        if (algo_param%Mechanics) then
          do i=1,size(MattotP)
             MattotP(i)%VarInt0=MattotP(i)%VarInt
          end do
          !Sig0 = Sig
          ! modif intel15/oneapi2022 pour cas test fissure_1proc sur marquises
          do i=1,6
          do j= 1,size(Sig0,dim=1)
             Sig0(j,i) = Sig(j,i)
          end do
          end do
          sigMoy0 = sigMoy
          defMoy0 = defMoy
        end if 
        if (algo_param%Diffusion) then
          FluxD0 = FluxD
          GradQDMoy0 = GradQDMoy
          FluxDMoy0  = FluxDMoy 
        end if 

        !================================================================================
        !                               ACTUALISATION MOYENNES EN DEBUT DE CHARG. PARTIEL
        !        
                   
        if (algo_param%Mechanics) then                  
        if(load_incr==load(load_n)%NIncr .OR. test_load_interruption)then  
            sigMoy00 = SigMoy
            defMoy00 = DefMoy  
        end if
        end if
        if (algo_param%Diffusion) then
        if(load_incr==load(load_n)%NIncr)then  
            GradQDMoy00 = GradQDMoy
            FluxDMoy00  = FluxDMoy 
        end if
        end if
        
        !================================================================================
        !                                                                   RESTART (END)
        !                                                                      save files
        !
        ! Not last increment
        if (load_incr /= 0 .and. load(load_n)%restart_save .and. mod(load_incr,load(load_n)%restart_save_period)==0 &
           .and. load_incr /= load(load_n)%Nincr) then
            call simu_restart_write(load_n,load_incr,defMoy00, sigMoy00, gradQDMoy00, FluxDMoy00)
        end if

        ! Last increment
        if (load_incr == load(load_n)%Nincr .and. load(load_n)%restart_save) then  
          if (load_n /=size(load)) then ! not last loading
            call simu_restart_write(load_n+1,0,defMoy00, sigMoy00, gradQDMoy00, FluxDMoy00)
          else                          ! last loading
            call simu_restart_write(load_n,load_incr,defMoy00, sigMoy00, gradQDMoy00, FluxDMoy00)
          end if
        end if

        !================================================================================
        !                                                         UPDATE ITERATION COUNTS

        nitTot=nitTot+nIt                                             ! Incrementation du nombre d'iterations totales
        
        VCinfo%nCall_lam_pas = VCinfo%nCall_lam_pas*(nIt+1)           ! Incrementation du nombre d'appel à umatLaminate pour ce pas
        VCinfo%nCall_lam = VCinfo%nCall_lam + VCinfo%nCall_lam_pas    ! Incrementation du nombre d'appel total à umatLaminate
        VCinfo%nIt_lam   = VCinfo%nIt_lam   + VCinfo%nIt_lam_pas      ! Incrementation du nombre d'iterations totales pour l'algorithme laminate
        VCinfo%nSub_lam  = VCinfo%nSub_lam  + VCinfo%nSub_lam_pas     ! Incrementation du nombre total de fois ou il a été 
                                                                      ! nécessaire de subdiviser le pas de temps de laminate
        VCinfo%nIncr_lam = VCinfo%nIncr_lam + VCinfo%nIncr_lam_pas    ! Incrementation du nombre de sous pas de temps total utilisé pour l'algorithme laminate


        !================================================================================
        !                                                                      LOG OUTPUT

        call log_after_unpas(load_incr,ind_tps,nIt,&
                  t_fft0,t_ifft0,t_behavior0,t_green0,t_crit0,t_wvtk0,t_wstd0,t_unpas0,&
                  testLaminate)
        
        call check_amitex_abort(0) !necessaire???

        !================================================================================
        !                                                                  LOAD INCREMENT
        !                                                          USER LOAD INTERRUPTION 

        ind_tps = ind_tps+1

        if (test_load_interruption .AND. load_incr .NE. 0) exit BOUCLE_INCREMENTS
        
!print *,"FIN BOUCLE INCREMENTS", ind_tps
     end do BOUCLE_INCREMENTS
     !                                               FIN BOUCLE SUR INCREMENTS D'UN CHARGEMENT PARTIEL
     !************************************************************************************************
  end do BOUCLE_CHARGEMENTS
  !                                                            FIN BOUCLE SUR LES CHARGEMENTS PARTIELS
  !***************************************************************************************************

  !================================================================================
  !                                                                FINAL LOG OUTPUT
  
  call log_final_times(nitTot)
  
end subroutine resolution_NL_base


!==============================================================================
!      SUBROUTINE UPDATE_LOCAL_LOAD
!
!>  RECUPERATION DU CHARGEMENT COURANT (load_n,load_incr)
!!                                      -> local_load(D)%t0, t1 and %dt
!!                                      -> MAJ local_load%t_load        
!!  MAJ local_load%t_load :   t(-2)=t(-1), t(-1)=t(0), t(0)=load(load_nb)%time(incr)
!!
!! RAPPEL (voir loading_mod.f90)
!! MECA
!! local_load%t1(1:algo_param%nTensDef) :                       pilotage 
!! local_load%t1(algo_param%nTensDef+1:2*algo_param%nTensDef) : valeurs associées au pilotage
!! local_load%t1(2*algo_param%nTensDef+1) :                     temperature
!! local_load%t1(nb_param next indices) :                       parametres externes (s'ils existent)
!! local_load%t1(27 next indices) :                             gradgradU components (if exists)
!!
!! DIFFUSION
!! local_loadD%t1(1:3*algo_param%nVarD) :                        pilotage 
!! local_loadD%t1(3*algo_param%nVarD+1:6*algo_param%nVarD) :     valeurs associées au pilotage
!! local_loadD%t1(6*algo_param%nVarD+1) :                        temperature
!! local_loadD%t1(6*algo_param%nVarD+2:) :                       parametres externes (s'ils existent)
!!
!!
!!  \param[in]    load_n, numero de chargement partiel
!!  \param[in]    load_incr, indice du pas de calcul
!!  \param[in]    defMoy     | grandeurs moyennes en debut de pas
!!  \param[in]    sigMoy     |    -> utile en cas de "maintien" au changement 
!!  \param[in]    gradQDMoy  |       de chargement partiel
!!  \param[in]    FluxDMoy   |
!!
!------------------------------------------------------------------------------
subroutine update_local_load(load_n,load_incr,defMoy, sigMoy, gradQDMoy, FluxDMoy)

  implicit none

  integer,intent(in)                                     :: load_n, load_incr
                                                                       !< numero du chargement partiel 
                                                                       !< numero de l'increment du chargement 
                                                                       !<(load_incr=0 : increment supplementaire en deb. chargement)
  real(mytype),dimension(algo_param%nTensDef),intent(in)   :: defMoy, sigMoy      !< Contrainte, def... moyennes en deput de chargement partiel
  real(mytype),dimension(3,algo_param%nVarD),intent(in)    :: gradQDMoy, FluxDMoy !< utile en cas de maintien (fluage par ex.)


  if (algo_param%Mechanics) then
     local_load%t0 = local_load%t1
     call get_current_loading(load_n,load_incr,defMoy, sigMoy, local_load%t_load,&
                              local_load%t0,local_load%t1)
     local_load%dt = local_load%t1 - local_load%t0
  end if

  if (algo_param%Diffusion) then
     local_loadD%t0 = local_loadD%t1
     call get_current_loadingD(load_n,load_incr,gradQDMoy, FluxDMoy, local_loadD%t_load,&
                               local_loadD%t0,local_loadD%t1)
     local_loadD%dt = local_loadD%t1 - local_loadD%t0
  end if 

end subroutine update_local_load
!============================================================================== 



!===========================================================================================================
!> SUBSTEPPING
!!
!! Recursive function :
!!  - nothing done if : convergence or maximal depth achieved or first 'fictious' increment
!!  - Nitermax increases with depth according to a geometric progression : 
!!           algo_param%Nitermax2 =   algo_param%Nitermax2  * (algo_param%Geom_ratio)^depth   
!!  - uses the best DEF stored in CVFORsauv%DEF to initialize the first sub-step         
!!  
!! L. GELEBART 01 / 2023
!===========================================================================================================
recursive subroutine substepping(load_n,load_incr,ind_tps)
  implicit none
 
  integer,intent(in)                   :: load_n,load_incr, ind_tps
  
  real(mytype)                               :: dtsub
  integer                                    :: i,j,nit, Ntensdef
  integer,save                               :: current_depth = 0
  character(len=100)                         :: tmp_char
  character(len=16)                          :: init_def
  real(mytype), allocatable, dimension(:)    :: load_t0
  real(mytype), allocatable, dimension(:)    :: load_t1
  real(mytype), allocatable, dimension(:)    :: load_dt
  real(mytype), dimension(-2:0)              :: load_t_load
  type(LOCALLOAD)                            :: local_load_sub

  
  Ntensdef = algo_param%nTensDef
  init_def = algo_param%init_def
  
  !> TEST - SUBSTEPPING 
  if (info_NL_b%reach_nitermax .and. load_incr .ne. 0 .and. current_depth .le. algo_param%depth) then

    info_NL_b%reach_nitermax = .false.
    current_depth = current_depth + 1    
    algo_param%Nitermax2 = algo_param%Nitermax2 * algo_param%Geom_ratio
    
    !> STORE current loading 
    load_t0     = local_load%t0      ! allocations à la volee.
    load_t1     = local_load%t1      
    load_dt     = local_load%dt      
    load_t_load = local_load%t_load      

    !> INITIALIZE subloading
    local_load_sub%t0     = load_t0      ! allocations à la volee.
    local_load_sub%t1     = load_t1  
    local_load_sub%dt     = load_dt 
    local_load_sub%t_load = load_t_load  
    
    local_load_sub%dt(nTensDef+1:2*nTensDef) = load_dt(nTensDef+1:2*nTensDef) / algo_param%Nsub    
    dtsub = (load_t_load(0) - load_t_load(-1)) / algo_param%Nsub

    !>===================================================== LOOP OVER SUB-STEPS
    !>
    do i = 1, algo_param%Nsub    
       !>===============================================================
       !                                           BEFORE SOLVE SUB-STEP

       ! MODIFY INIT_DEF for i = 1
       if (i == 1) then
          !DEF = DEF0 + (DEF - DEF0) / algo_param%Nsub           ! use last unconverged DEF
          DEF = DEF0 + (CVFORsauv%DEF - DEF0) / algo_param%Nsub  ! use best unconverged DEF
          !DEF = DEF0                                            ! use DEF0 = DEF 
          CVFORsauv%critmin = huge(1._mytype)
          algo_param%init_def="none"
       end if
       
       ! ADJUST LOADING
       local_load_sub%t0(nTensDef+1:2*nTensDef) = load_t0(nTensDef+1:2*nTensDef) &
                                               +  (i -1) * local_load_sub%dt(nTensDef+1:2*nTensDef)
       local_load_sub%t1(nTensDef+1:2*nTensDef) = load_t0(nTensDef+1:2*nTensDef) &
                                               +  i * local_load_sub%dt(nTensDef+1:2*nTensDef) 
       local_load_sub%t_load(0) = load_t_load(-1) + i * dtsub

       !>===============================================================
       !                                                  SOLVE SUB-STEP
       
       ! UNPAS (WARNING : local_load = local_load_sub in unpas_NL_base)
       call unpas_NL_base(load_n,load_incr,ind_tps,nIt,local_load_sub)
       write(tmp_char,*) "================SUB-STEPPING - depth ", current_depth," step ", i , "nit", nit
       call write_stdout0(trim(tmp_char))

       ! RECURSIVE SUBSTEPPING
       if (info_NL_b%reach_nitermax .and. current_depth < algo_param%depth) then
          call substepping(load_n,load_incr,ind_tps)
       end if

       ! WRITE VTK BEFORE "QUIT" IF NON-CONVERGENCE
       if (info_NL_b%reach_nitermax .and. current_depth == algo_param%depth) then
          if(nrank == 0) write(Flog,"(A)")"Substepping failure - Write VTK before quit."
          call writeVTK(ind_tps)
          call amitex_abort("Supstepping - maximum depth reached without convergence", 2, 0)
       end if
       
       ! UPDATE FIELDS
       do j=1,size(MattotP)
         MattotP(j)%VarInt0=MattotP(j)%VarInt
       end do
       Sig0 = Sig
       !Def0 = Def  ! DON'T!!! (update done in initDef in unpas_NL_base!)
       
       !sigMoy0 = sigMoy  ! user_load_interruption will not work
       !defMoy0 = defMoy
       
       !>===============================================================
       !                                            AFTER SOLVE SUB-STEP
       
       ! BACK INIT_DEF for i = 1
       if (i == 1) algo_param%init_def = init_def

       ! UPDATE TIMES
       local_load_sub%t_load(-2) = local_load_sub%t_load(-1)
       local_load_sub%t_load(-1) = local_load_sub%t_load(0)
              
    end do
    
    !>================================================= AFTER LOOP OVER SUB-STEPS
    !> IF CONVERGENCE OK : 
    !>       - DECREMENT CURRENT DEPTH
    !>       - DIVIDE 'BACK' NITERMAX2
    !>       - BACK to initial local_load (for use at depth-1)
    !>       - NEW DEF0 (for use at depth-1)
    if (.not. info_NL_b%reach_nitermax) then
       current_depth = current_depth -1
       algo_param%Nitermax2 = algo_param%Nitermax2 / algo_param%Geom_ratio
       
       local_load%t0     = load_t0           
       local_load%t1     = load_t1           
       local_load%dt     = load_dt           
       local_load%t_load = load_t_load   
            
       Def0 = Def - (Def - Def0) * algo_param%Nsub  
    end if

  end if
  !> END SUBSTEPPING 

end subroutine substepping


end module resolution_mod
