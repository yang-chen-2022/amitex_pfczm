!123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!===================================================================================================
!
!       MODULE NL_BASE_MOD 
!>
!> SOLVE 'ONE STEP' ('UN PAS' in french!) WITH THE BASIC SCHEME (WITH CONVERGENCE ACCELARATION)
!!
!!  
!!  Subroutines
!! - unpas_NL_base : solves one step (unpas) with the basic Scheme (Moulinec-Suquet)
!!                   + CAST3M convergence acceleration (ACT3)
!!
!! - log_after_unpas, log_final_times, write_crit_log : 
!!                    Write in "log" file and/or standard output
!
module NL_base_mod

  use ISO_FORTRAN_ENV
  use IEEE_arithmetic

  use MPI             ! on charge les modules afin d'avoir acces aux fonctions
#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi,       only : nrank
#else
  use decomp_2d, only : nrank,mytype 
#endif  

  use io_amitex_mod,  only: write_stdout0 ! function
  use io2_amitex_mod
  use algo_functions_mod
  use param_algo_mod
  use loading_mod
  use material_mod
  use green_mod
  use error_mod
  use amitex_mod
  use field_mod
  use non_local_mod
  use sortie_std_mod
  use amitex_mpi_mod, only: amitex_mpi_comm
!  use user_functions_mod

  private
   
  !> Pointer "publiques" vers composantes de SIMU_AMITEX
  public ::  crit_b, times_b, info_NL_b

  !> Variables "publiques" utilisees lors de l'initialisation
  ! nothing

  !> Types publiques (pour definition de SIMU_AMITEX)
  public ::  CRIT_BASE, TIMES_BASE, INFO_NL_BASE

  !> Fonctions publiques
  public :: unpas_NL_base, log_after_unpas, log_final_times, write_crit_log

  !!---------------------------------------------------------------
  !>                                   TYPES CRIT_BASE & TIMES_BASE
  !>                     gather criteria used by resolution_NL_base
  !> 
  type CRIT_BASE
                                                        !< MECANIQUE
    real(mytype)                 :: eq                  !< critere d'equilibre
    real(mytype)                 :: SigMoy,DefMoy       !< critere sur les moyennes
    real(mytype)                 :: Cptb                !< critere de compatibilite
  
                                                                  !< DIFFUSION
    real(mytype),allocatable,dimension(:)  :: eqD                 !< critere d'equilibre 
    real(mytype),allocatable,dimension(:)  :: FluxDMoy,GradQDMoy  !< critere sur les moyennes
    real(mytype),allocatable,dimension(:)  :: CptbD               !< critere de compatibilite

  end type CRIT_BASE
  
  type TIMES_BASE
     double precision       :: startalgo    = 0  
     double precision       :: startstep    = 0  
     double precision       :: unpas    = 0       
  end type TIMES_BASE
  
  type INFO_NL_BASE
     logical                :: reach_nitermax = .false. 
     !integer                :: nit_step  = 0 (TODO)
     !integer                :: nit_tot   = 0
     integer                :: nit_unpas  = 0
  end type INFO_NL_BASE
  
  type(CRIT_BASE), pointer       :: crit_b 
  type(TIMES_BASE), pointer      :: times_b 
  type(INFO_NL_BASE), pointer    :: info_NL_b
  
contains

!===========================================================================================================
!!==========================================================================================================
!> unpas_NL_base : Schema de base pour resoudre un probleme non lineaire
!!                       en utilisant des transformees de Fourier rapides
!!                       Cette procédure contient l'algorithme de resolution
!!                       la procédure traite la mécanique et la diffusion
!!
!!                 Resoud le probleme sur un pas, est appele par resolution_NL_base 
!!                       qui deroule les pas successifs
!!
!! WARNING 1 : if algo_param%init_def/='none' 
!!             at the begining :
!!               DEF0 (before call) = def. at previous step (n-2) 
!!               DEF  (before call) = def. at last step (n-1)
!!             after iniDef :
!!               DEF0 = def.. at last step (n-1)
!!               DEF  = initialization of def. at current step (n)
!!             at the end :
!!               DEF0 (after call) = def. at last step (n-1)
!!               DEF (after call)  = def. at current step (n) 
!!
!! WARNING 2 : if locload is given as argument, then local_load is overloaded (VERY IMPORTANT!)
!!
!! NOTATIONS (entree xml et champs) Sigma et Def (HPP) : 
!!                 CAST3M (Voigt + ordre 11 22 33 12 13 23)
!!                 PK1 et grad(u) (GDEF) : ordre 11 22 33 12 13 23 21 31 32 sans notation
!!
!!
!!                        LISTE DES VARIABLES MODIFIEES (A finir!)
!!--------------------------------------------------------
!!   VCinfo%nIt_lam_pas, nSub_lam_pas, nIncr_lam_pas, nCall_lam_pas
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
!!   VCinfo      infos sur integration modele composite
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
!!   Def_star Champ de deformation additionnel defini tel que : Def = E + Def_per + Def_star
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
!!   Flog	Unites logiques pour l'ouverture du fichier .log
!!   fic_vtk    Racine du nom des fichiers de sortie vtk ou .(mz)std
!!
!!                                Variables dans GREEN_MOD
!!--------------------------------------------------------
!!   FREQ     Tableau des fréquences, espace Fourier 
!!            pinceaux-Z (nx/2+1,ny,nz,3)
!!            indices 3D : (fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),3)
!!
!!                              Variables dans LOADING_MOD
!!--------------------------------------------------------
!!   local_load%t1 et t0 et dt   chargement mecanique et Temp et Param. Ext. impose a chaque iteration
!!                               %t1 est utilise pour imposer le chargement "méca"
!!                               %t0 et %dt sont utilises dans behavior pour imposer un chargement en T et Param Ext.      
!!
!!   local_loadD%t1 et dt        chargement "diffusion" impose a chaque iteration
!!
!!                              Variables dans NL_BASE_MOD
!!--------------------------------------------------------
!!   crit_b                criteria for resolution_NL_base
!!   times_b               elapsed times for resolution_NL_base
!!
!!------------------------------------------------------------
!!
!!   \param[in]    ind_tps (indice de pas de temps)
!!   \param[in]    load_n,load_incr
!!                 numero de chargement partiel et indice du pas de calcul
!!   \param[out]   nIt         
!!                 nombre d'iterations effectuees
!!   \param[in]    locload  (optional)   Variables de type LOCALLOAD. Si presence :
!!   \param[in]    locloadD (optional)     surcharge local_loading / local_loadingD  
!!                                         apres get_current_loading
!!
!!------------------------------------------------------------------------------
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
!!
!!===============================================================================

subroutine unpas_NL_base(load_n,load_incr,ind_tps,nIt,locload,locloadD)

!!==============================================================================
!!                                                                  DECLARATIONS
!!==============================================================================

  implicit none
!!------------------------------------------------------------------------------
!>                                                     PARAMETRES "PARALLELISME"

  integer                          :: ierror             !< erreur relative au fonction MPI
  
!!------------------------------------------------------------------------------
!>                                          PARAMETRES "DIMENSIONS DE LA GRILLE"

  integer                          :: ni2                !< dimension sur 2


!!------------------------------------------------------------------------------
!>                                            PARAMETRES "MATERIAU DE REFERENCE"

  real(mytype)                     :: lambda, mu

!!------------------------------------------------------------------------------
!>                                                       PARAMETRES "CHARGEMENT"

  integer,intent(in)                                     :: ind_tps    !< indice du pas de temps courant 
  real(mytype),dimension(algo_param%nTensDef)            :: Eimp       !< deformation moyenne imposee
  real(mytype),dimension(3*algo_param%nVarD)             :: EimpD      !< gradQ moyen impose
  integer,intent(in)                                     :: load_n, load_incr
                                                                       !< numero du chargement partiel 
                                                                       !< numero de l'increment du chargement 
                                                                       !<(load_incr=0 : increment supplementaire en deb. chargement)
  integer                                                :: nb_stress_pilot,nb_flux_pilot 
                                                                       !< nombre de directions pilotees en contrainte
                                     
!!------------------------------------------------------------------------------
!>                                              PARAMETRES "CHARGEMENT" OPTIONEL

  type(LOCALLOAD), intent(in), optional                  :: locload, locloadD

!!------------------------------------------------------------------------------
!>                                       ACCELERATION DE CONVERGENCE (SI ACTIVE)
!>                                            tableaux paralleles (ntot,Ntens,4)
!>                                          tableaux paralleles (ntot,3*nVarD,4)
!>  increments de champs de deformation et residus des 4 pas de temps precedents

  integer                           :: iact3           !< compteur pour ACV
  logical                           :: lact3           !< indice d'activation de l'ACV
  logical                           :: acv_test, acv_testD        
                                                                  !< sortie act3 : 
                                                                  !< permet de ne pas stocker un champ non accelere

!!------------------------------------------------------------------------------
!>                                                CONVERGENCE FORCEE (SI ACTIVE)
!>                     
  logical                           :: testCVFOR        !< test pour activer la convergence forcee
  integer                           :: countCVFOR       !< compteur de convergence forcee 

!!------------------------------------------------------------------------------
!>                                                             NB ITERATIONS MAX
!>                     
  integer                            :: Nitermax    ! algo_param%Nitermax or algo_param%Nitermax2
  
!!------------------------------------------------------------------------------
!>                                           PARAMETRES DE SUIVI DE L'ALGORITHME

  integer,intent(out)               :: nIt                    !< nombres d'iterations dans un pas de chargement

!!------------------------------------------------------------------------------
!>                                  PARAMETRES DE SUIVI DE L'ALGORITHME LAMINATE

  !< suivi performance algorithme laminate 
  !<               nSub_laminate         : nbre de recours à la subdivision du pas de temps 
  !<               nIncr_laminate        : nbre de sous pas de temps par subdivision 
  !<               nIt_laminate          : nbre d'itérations par appel de l'algorithme laminate
  !<               nCall_laminate        : nbre d'appel à l'algorithme laminate
  integer(kind=8)               :: nSub_laminate_temp,nIncr_laminate_temp,nIt_laminate_temp !< pour un appel au comportement (behavior)

!!------------------------------------------------------------------------------
!>                                            PARAMETRES DE L'ALGO DE RESOLUTION

  logical                                                  :: test_out            !< test de sortie de boucle while


!!------------------------------------------------------------------------------
!>                                                          COMPATIBILITE ABAQUS
  integer(kind=8)                    :: KSTEP 

!!------------------------------------------------------------------------------
!>                                     PARAMETRES DE SUIVI DES TEMPS D'EXECUTION

  !! type double precision au lieu de real_mytype (=> requis par MPI_WTIME)
  double precision                   :: t_ini,t1,t10

  !< t1          : date de repere avant les fonctions
  !! t10         : temps cumule pris par les 10 dernieres iterations
  
!!------------------------------------------------------------------------------
!>                                                                        DIVERS

  integer               :: i,j,k,l,m !< indice de boucle 
  integer               :: nVarD     !< nombre de variables de diffusion 
  logical               :: testNaN1,testNaN2,testdef

  integer               :: TPRED_deb !< indice in local_load%t1 for T 
  integer               :: TPRED_fin !< indice in local_load%t1 for the last param_ext

!!==============================================================================
!!                                                                    INIT TIMES
!!==============================================================================
  call MPI_BARRIER(amitex_mpi_comm,ierror)
  t_ini = MPI_WTIME() 
  t10 = t_ini  !< date avant debut boucle while pour cumul 10 iterations

!!==============================================================================
!!                                                           AUXILIARY VARIABLES
!!==============================================================================

  TPRED_deb = 2*algo_param%nTensDef+1 
  TPRED_fin = 2*algo_param%nTensDef+1+initValExt%nb_param    

!================================================================================ 
!                                                                 INITIALISATIONS
!

  lambda=0; mu=0
  nVarD = algo_param%nVarD
  if (algo_param%Mechanics) then
    lambda = Matref%LambdaMu0(1)    !< materiau de référence
    mu = Matref%LambdaMu0(2)
    Eimp=0
  end if
  if (algo_param%Diffusion) then  
    EimpD=0
  end if
  nSub_laminate_temp=0
  nIncr_laminate_temp=0
  nIt_laminate_temp=0

  Nitermax = algo_param%Nitermax
  if (algo_param%substepping .and. load_incr /=0) Nitermax = algo_param%Nitermax2

       !================================================================================
       !                                                       BOUCLE CONVERGENCE FORCEE
       !
       !  Forced convergence stops when testCFFOR is set to .false.
       ! 
       testCVFOR  = .true.
       countCVFOR = 0
       do while (testCVFOR)

        VCinfo%nSub_lam_pas          = 0    !< indicateurs performances algorithme 
        VCinfo%nIncr_lam_pas         = 0    !  laminate pour le pas de chargement
        VCinfo%nIt_lam_pas           = 0
        VCinfo%nCall_lam_pas         = 0
        !----------------------------------------------------------------------
        !                                    GESTION DU DEMI-PAS SUPPLEMENTAIRE

        !! Le chargement load_incr=0 correspond a un chargement intermediaire
        !! entre le chargement final impose a load_n-1
        !! et la premiere iteration du chargement load_n 
        !! avec un ratio defini dans loading_mod
        KSTEP = ind_tps       ! attention, du fait du pas de temps fictif, 
                              ! KSTEP est repete deux fois a chaque debut de chargement

        if(nrank==0)then
          if(load_incr == 0)then
              write(OUTPUT_UNIT,"(/,A)")"=================================================="
              write(OUTPUT_UNIT,"(A,I8,A)")"Pas de temps fictif :",ind_tps,".5"
              write(Flog,"(/,A)")"=================================================="
              write(Flog,"(A,I8,A)")"Pas de temps fictif :",ind_tps,".5"
          else
              write(OUTPUT_UNIT,"(/,A)")"=================================================="
              write(OUTPUT_UNIT,"(A,I8)")"Pas de temps :",ind_tps
              write(Flog,"(/,A)")"=================================================="
              write(Flog,"(A,I8)")"Pas de temps :",ind_tps
          end if 
        end if

        !----------------------------------------------------------------------      
        !                                EN PRESENCE D'UN CHARGEMENT LOCAL_LOAD      
        !

        if (algo_param%Mechanics .AND. present(locload)) local_load = locload
        if (algo_param%Diffusion .AND. present(locloadD)) local_loadD = locloadD

        !----------------------------------------------------------------------      
        !                                                    ANALYSE CHARGEMENT      
        !

        if (algo_param%Mechanics) nb_stress_pilot = count(local_load%t1(1:algo_param%nTensDef) == 1)
        if (algo_param%Diffusion) nb_flux_pilot = count(local_loadD%t1(1:3*nVarD) == 1)  

        !----------------------------------------------------------------------
        !                                 PROPOSITION D'UN CHAMP DE DEFORMATION
        !
        ! avec CV forcee : 
        !    1ere initialisation identique au cas sans CV forcee
        !    PUIS
        !    si algo_param%init_cvfor=='default' : 
        !       initialisation identique au cas sans CV forcee (DEF -> DEF0 et initialisation de DEF)
        !    si algo_param%init_cvfor=='best': 
        !       on recupere le champ DEF choisi en fin de CV forcee precedente 
        !    'default' semble donner les meilleurs resultats (en endommagement)
        !    'best' semble donner de bons resultats en plasticite
        ! sans CV forcee :  le test ci-dessous est toujours vrai
        !
        if (algo_param%init_cvfor=='default' .OR. countCVFOR == 0) then
        if (algo_param%Mechanics) then
            if (algo_param%init_def /= 'none') then !-- peut etre mis a 'none' dans un fonction user pour proposer une autre initialisation 
            if (algo_param%HPP_nsym) then
                ! prise en compte du gradient complet du deplacement dans un 
                ! calcul HPP : initialisation de Def_nsym
                call initDef(local_load%t1(1:2*algo_param%nTensDef),&
                         nb_stress_pilot,local_load%t_load,load_incr,load_n,Matref%C0,Def0,Def,&
                         Def_nsym0,Def_nsym); !load_incr,load_n  pour tester cas particulier 0 et 1
            else 
                call initDef(local_load%t1(1:2*algo_param%nTensDef),&
                         nb_stress_pilot,local_load%t_load,load_incr,load_n,Matref%C0,Def0,Def);

            end if
            else
                call write_stdout0("STRAINS (DEF,DEF0) NOT INITIALIZED IN UNPAS_NL_BASE")
            end if
        end if

        if (algo_param%Diffusion) then
            call initDefD(local_loadD%t1(1:6*algo_param%nVarD),&
                         nb_flux_pilot,local_loadD%t_load,load_incr,load_n,Matref%K0D,GradQD0,GradQD);
        end if
        end if

        !----------------------------------------------------------------------
        !                              INITIALISATION DES CRITERES ET COMPTEURS

        if (algo_param%Diffusion) then
            crit_b%eqD=0_mytype
            crit_b%CptbD=0_mytype
            crit_b%FluxDMoy=0_mytype
            crit_b%GradQDMoy=0_mytype  
        end if
        if (algo_param%Mechanics) then
            crit_b%eq=0_mytype
            crit_b%Cptb=0_mytype
            crit_b%DefMoy=0_mytype
            crit_b%SigMoy=0_mytype
        end if

        nIt=-1  ! nombre d'iterations pour le pas de chargement courant
                ! passe a 0 a l'entree de la boucle while
                ! puis l'algo force au moins une iteration 
        iact3=0
        lact3=.false.
        acv_test = .true.
        acv_testD = .true.

        !----------------------------------------------------------------------
        !                                                          DefF=FFT(Def)
        !              format DefF: (pinceaux-Z, taille globale (nx/2+1)*ny*nz)
        !                                             on normalise DefF par ntot

!PATCH TEMPORAIRE !!
!TODO URGENT : patch qui fonctionne pour AME + def_star avec itmax=1
!              A REVOIR (egalement pour gradgradU)
!              Pour calculer DefF :
!                  Def ne doit comporter que la composante périodique!!!
! A REFLECHIR : travailler sur Def = Def_tilde + def_moy
!               ajouter/enlever def_star (et/ou gradgradU) avant/apres le calcul de la contrainte 
if (add_def_star .AND. load_incr/=0 .and. algo_param%init_Def=='default')  Def = Def - Def_star           
        if (algo_param%Mechanics) then
           call field_fft(Def,DefF,algo_param%nTensDef,1)
           DefF=DefF/ real(grid%ntot,mytype)
        end if
if (add_def_star .AND. load_incr/=0 .and. algo_param%init_Def=='default')  Def = Def + Def_star           

        if (algo_param%Diffusion) then
           call field_fft(gradQD,GradQDF,3,algo_param%nvarD)
           GradQDF=GradQDF / real(grid%ntot,mytype)
        end if

       !pour tester calcul critere en deformation sur un champ simple
       !  call eval_criteres(SigF,DefF,FREQ,local_load%t1(1:2*algo_param%nTensDef),&
       !        load(load_n)%DirStress_flag, load(load_n)%DirStress2,& 
       !        crit_b%eq,crit_b%SigMoy,crit_b%DefMoy,crit_b%Cptb,grid%ntot)


        !----------------------------------------------------------------------
        !                           POUR LE 1er PAS DE CHARGEMENT (load_incr=0)
        !         PROPOSITION D'UN CHAMP DE DEFORMATION (SUITE): AJOUT gragradU
        !                                                        AJOUT Def_star
        !
        !                      on doit ajouter gragradU APRES avoir evalue DefF
        !
        if (algo_param%init_cvfor=='default' .OR. countCVFOR == 0) then
        if (algo_param%Mechanics) then
            ! Cas gradgradU impose
            ! Def = def + (gradgradU) (_sym si HPP)
            if (n_gradgradU==27 .AND. load_incr==0) then !test pas tres joli <=> test gradgradU impose
               call add_gradgradU(Def,algo_param%nTensDef,&
                          local_load%t1(TPRED_fin+1:TPRED_fin+27),load(load_n)%gevolution,&
                           grid%dx,grid%dy,grid%dz,grid%nx,grid%ny,grid%nz)
            end if
            if (add_def_star .AND. load_incr==0) then 
               Def = Def + Def_star
            end if            
         end if
         end if

        !----------------------------------------------------------------------
        !                                                          COMPORTEMENT

        if (algo_param%Mechanics) then !----- MECANIQUE
             if (algo_param%HPP_nsym)  then
                ! prise en compte du gradient complet du deplacement dans un 
                ! calcul HPP 
                call behavior(SIG,SIG0,DEF,DEF0,local_load%t0(TPRED_deb:TPRED_fin),&
                           local_load%dt(TPRED_deb:TPRED_fin),local_load%t_load,KSTEP,&
                           nSub_laminate_temp,nIncr_laminate_temp,nIt_laminate_temp,VCinfo%nCall_lam_pas,&
                           Def_nsym=Def_nsym,Def_nsym0=Def_nsym0)
             else
                call behavior(SIG,SIG0,DEF,DEF0,local_load%t0(TPRED_deb:TPRED_fin),&
                           local_load%dt(TPRED_deb:TPRED_fin),local_load%t_load,KSTEP,&
                           nSub_laminate_temp,nIncr_laminate_temp,nIt_laminate_temp,VCinfo%nCall_lam_pas)
 
             end if
           if(.not. algo_param%HPP) then
              call SigToPK1(Sig,Def,PK1)
           end if
           if (algo_param%Nloc .and. .false. ) then !---- Mecanique non locale resolution implicite (A REVOIR)
              !TODO revoir la possibilite d'introduire de l'implicite
              call Nloc_call(load_n,load_incr,ind_tps,"after", algo_param%nloc_implicit)
           end if
        end if

        if (algo_param%Diffusion) then !----- DIFFUSION
           call behaviorD(FluxD,FluxD0,GradQD,GradQD0,local_loadD%t0(6*algo_param%nVarD+1:),&
                          local_loadD%dt(6*algo_param%nVarD+1:),local_loadD%t_load,KSTEP)
        end if

        ! Incrémentation des compteurs après appel au comportement
        !       a ce moment nCall_laminate contient le nombre de voxels composites "laminate"
        !       pour ce pas de chargement (en prévision de cas ou ce nombre changerait 
        !       d'un pas de chargement à l'autre)
        VCinfo%nSub_lam_pas   = VCinfo%nSub_lam_pas  + nSub_laminate_temp
        VCinfo%nIncr_lam_pas  = VCinfo%nIncr_lam_pas + nIncr_laminate_temp
        VCinfo%nIt_lam_pas    = VCinfo%nIt_lam_pas   + nIt_laminate_temp

        !================================================================================
        !                                                                    BOUCLE WHILE
        do
          nIt=nIt+1   ! actualisation du compteur de boucle, 
                      ! au 1er passage avant tests des criteres on a i=0

!interet de cette initialisation des criteres moyens?
!         if (algo_param%Mechanics) then
!              crit_b%DefMoy=0
!              crit_b%SigMoy=0  
!         end if
!         if (algo_param%Diffusion) then
!              crit_b%FluxDMoy=0
!              crit_b%GradQDMoy=0      
!         end if
           !----------------------------------------------------------------------
           !                                FFT DE LA CONTRAINTE : SigF = FFT(Sig)
           !              format SigF: (pinceaux-Z, taille globale (nx/2+1)*ny*nz)
          

           if(algo_param%Mechanics) then  !-------MECANIQUE
           if(algo_param%HPP) then !>--HPP
              call field_fft(Sig,SigF,algo_param%nTensDef,1)
           else                    !>---GRANDES DEF
              call field_fft(PK1,SigF,algo_param%nTensDef,1)
              call mpi_barrier(amitex_mpi_comm,ierror)
              ! correction dimensions paires
              if (algo_param%CorrPair) then
              if(modulo(grid%ny,2) == 0) then
                 ni2 = grid%ny/2 +1
                 if(ni2>= fft_start(2) .AND. ni2<= fft_end(2) ) then 
                    SigF(fft_start(1):fft_end(1),ni2,fft_start(3):fft_end(3),1:9)=&
                         cmplx(0,kind=mytype)
                 end if
              end if
              if(modulo(grid%nz,2) == 0) then
                 ni2 = grid%nz/2 +1
                 if(ni2>= fft_start(3) .AND. ni2<= fft_end(3) ) then 
                    SigF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),ni2,1:9)=&
                         cmplx(0,kind=mytype)
                 end if
              end if
              end if
           end if 
           end if

           if(algo_param%Diffusion) then  !-------DIFFUSION
              call field_fft(FluxD,FluxDF,3,algo_param%nvarD)
           end if

!print *, "AV ACV",nrank,shape(ACT3_U),shape(Def)
           !----------------------------------------------------------------------
           !                                 STOCKAGE ACCCELERATION DE CONVERGENCE
           ! TODO : dissocier les deux acv meca et diffusion

           ! on ne stocke pas s'il n'y a pas eu ACV (acv_test = .false.)
           if(algo_param%acc_CV .and. acv_test .and. acv_testD) then
              iact3 = iact3+1
              if (iact3==5) then 
                 iact3 = 1
              end if
              if (algo_param%Mechanics) then
                if (algo_param%HPP_nsym) then
                ! prise en compte du gradient complet du deplacement dans un 
                ! calcul HPP : stockage de gradient de u pour l'acceleration
                    !ACT3_U(:,:,iact3) = Def_nsym(:,:)
                    ! modif intel15/oneapi2022 pour cas test fissure_1proc sur marquises machine batch 28p (ok sur machine interactive 14p)
                    do i = 1,algo_param%nTensDef
                    do j = 1, size(ACT3_U,dim=1)
                       ACT3_U(j,i,iact3) = Def_nsym(j,i)
                    end do
                    end do

                else
                    !ACT3_U(:,:,iact3) = Def(:,:)
                    ! modif intel15/oneapi2022 pour cas test fissure_1proc sur marquises machine batch 28p (ok sur machine interactive 14p)
                    do i = 1, algo_param%nTensDef
                    do j = 1, size(ACT3_U,dim=1)
                       ACT3_U(j,i,iact3) = Def(j,i)
                    end do
                    end do
                end if
              end if
              if (algo_param%Diffusion) then
                  do j = 1,nVarD
                  do i = 1,3
                  do l = 1, size(ACT3_UD,dim=1)
                     k = i+3*(j-1)
                     ACT3_UD(l,k,iact3) = GradQD(l,i,j)           
                  end do         
                  end do
                  end do
              end if
           end if
!print *, "AP ACV",nrank

           !----------------------------------------------------------------------
           !                                 EVALUATION DES CRITERES EN CONTRAINTE
           !                                                          + VERIF. NaN
           
           if(algo_param%Mechanics) then !----------MECANIQUE
           if(algo_param%HPP) then
              call eval_criteres(SigF,DefF,FREQ,local_load%t1(1:2*algo_param%nTensDef),&
                   load(load_n)%DirStress_flag, load(load_n)%DirStress2,& 
                   crit_b%eq,crit_b%SigMoy,crit_b%DefMoy,crit_b%Cptb,0_INT64)
              !A DECOMMENTER pour critere en def. a chaque pas
              ! call eval_criteres(SigF,DefF,FREQ,local_load%t1(1:2*algo_param%nTensDef),&
              !     load(load_n)%DirStress_flag, load(load_n)%DirStress2,& 
              !     crit_b%eq,crit_b%SigMoy,crit_b%DefMoy,crit_b%Cptb,grid%ntot)
!if (nrank==0) print *, "CRITEQ                                   =   ", crit_b%eq
           else
              call eval_criteres_GD(SigF,DefF,FREQ,local_load%t1(1:2*algo_param%nTensDef),&
                   load(load_n)%DirStress_flag, load(load_n)%DirStress2,& 
                   crit_b%eq,crit_b%SigMoy,crit_b%DefMoy,crit_b%Cptb,0_INT64)
              !
              !block pour tester l'utilisation du criteres HPP
              !------------------------------------------------
              !block
              !complex(mytype),allocatable,dimension(:,:,:,:)  ::SigF2
              !allocate(SigF2(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:6))
              !SIGF2=(SigF(:,:,:,1:6)+ SigF(:,:,:,(/1,2,3,7,8,9/)))/2.
              !call eval_criteres(SigF2,DefF(:,:,:,1:6),&
              !     FREQ,(/ local_load%t1(1:6), local_load%t1(10:15) /),&
              !     load(load_n)%DirStress_flag, load(load_n)%DirStress2,& 
              !     crit_b%eq,crit_b%SigMoy,crit_b%DefMoy,crit_b%Cptb,0_INT64)
              !end block
           end if
           end if

           if(algo_param%Diffusion) then !----------DIFFUSION
              call eval_criteresD(FluxDF,GradQDF,FREQ,local_loadD%t1(1:6*nVarD),&
                   crit_b%eqD,crit_b%FluxDMoy,crit_b%GradQDMoy,crit_b%CptbD,0_INT64)        
           end if 

                                         !----------TESTS CRITERES NAN
           testNaN1=.false.
           testNaN2=.false.
           if (algo_param%Mechanics) then
                  testNaN1=ieee_is_nan(crit_b%eq)
                  testNaN2=ieee_is_nan(crit_b%SigMoy)
           end if
           if (algo_param%Diffusion) then
                do j = 1,nVarD
                        testNaN1 = testNaN1 .or. ieee_is_nan(crit_b%eqD(j))
                        testNaN2 = testNaN2 .or. ieee_is_nan(crit_b%FluxDmoy(j))
                end do
           end if
           if(testNaN1) then
              if(nrank==0) then
                 write (Flog,"(A,I8)")"-------- ERREUR : CRITERE D'EQUILIBRE NaN----- iteration ",nIT
                 if (algo_param%Mechanics) write (Flog,"(A,E15.8)") "crit_eq = ", crit_b%eq
                 if (algo_param%Diffusion) write (Flog,"(A,(E15.8))") "crit_eqD = ", crit_b%eqD
              end if
              call amitex_abort("Critere d'equilibre NaN",2,0)
           end if
           if(testNaN2) then
              if(nrank==0) then
                 write (Flog,"(A,I8)")&
                 "-------- ERREUR : CRITERE SUR LA CONTRAINTE (ou LE FLUX) MOYENNE NaN----- iteration ",nIT
                 if (algo_param%Mechanics) write (Flog,"(A,E15.8)") "critSigMoy = ", crit_b%SigMoy
                 if (algo_param%Diffusion) write (Flog,"(A,(E15.8))") "critFluxDmoy) = ", crit_b%FluxDmoy
              end if
              call amitex_abort("Critere sur la contrainte (ou le flux) moyenne NaN",2,0)
           end if

!print *, "AV STOCKAGE CV FOR",nrank
           !----------------------------------------------------------------------
           !                                           STOCKAGE CONVERGENCE FORCEE
           !                                                       OR SUB-STEPPING
           !
           ! nIT >1 est INDISPENSABLE dans le cas algo_param%init_cvfor='best'
           !        pour eviter de tomber sur des cycles d'iterations de CV forcee identiques
           if ((algo_param%CVfor .OR. algo_param%substepping) .AND. (nIt > 1)) then
           if (crit_b%eq < CVFORsauv%critmin) then
              CVFORsauv%critmin=crit_b%eq
              CVFORsauv%Def=Def
              if (nrank==0) write (Flog,"(A,(E15.8))") "crit_eq NON-CONV= ", crit_b%eq
              if (nrank==0) write (OUTPUT_UNIT,"(A,(E15.8))") "crit_eq NON-CONV= ", crit_b%eq
           end if
           end if

           !----------------------------------------------------------------------
           !                                  SORTIE CONDITIONNELLE BOUCLE INFINIE
           !
           !    criteres Contrainte (Flux) Macro et Equilibre < tol
           !    Si demande (Nitermin=1 et/ou Nitermin_acv=1) :
           !      minimum 1 itération (= application de Green)
           !      minimum 1 itération après accélération de convergence (= appli. Green)
           !
           !     => application de Green assure la compatibilité 
           !    (on verifie néanmoins le critere de compatibilité avant de sortie)
           !
           !                                  AVANT DE SORTIR: 

           !! EVALUATION DU TEST DE SORTIE
           test_out=.true.

           if (algo_param%NiterMin_acv==1) then 
              test_out = (.NOT. lact3)
           end if

           if (algo_param%NiterMin==1) then
              test_out = test_out .AND. (nIt .ge. 1) 
           end if

           if (algo_param%Mechanics) then 
                test_out = test_out .AND. (crit_b%eq .LT. algo_param%tol_criteq) .AND. &
                           (crit_b%SigMoy .LT. algo_param%tol_critsig)
           end if

           if (algo_param%Diffusion) then
                  test_out = test_out .AND. &
                        (max(maxval(crit_b%eqD),maxval(crit_b%FluxDMoy)) .LT. algo_param%tol_criteq)
           end if

           !! SORTIE
           if (test_out) then
              !>-----------EVALUATION DES CRITERES EN DEFORMATION (OU GRADQ)

              if(algo_param%Mechanics) then  !-------MECANIQUE
              if(algo_param%HPP) then
                 call eval_criteres(SigF,DefF,FREQ,local_load%t1(1:2*algo_param%nTensDef),&
                      load(load_n)%DirStress_flag, load(load_n)%DirStress2,& 
                      crit_b%eq,crit_b%SigMoy,crit_b%DefMoy,crit_b%Cptb,grid%ntot)
              else
                 call eval_criteres_GD(SigF,DefF,FREQ,local_load%t1(1:2*algo_param%nTensDef),&
                      load(load_n)%DirStress_flag, load(load_n)%DirStress2,& 
                      crit_b%eq,crit_b%SigMoy,crit_b%DefMoy,crit_b%Cptb,grid%ntot)
              end if
              end if
              if(algo_param%Diffusion) then !-------DIFFUSION
                 call eval_criteresD(FluxDF,GradQDF,FREQ,local_loadD%t1(1:6*nVarD),&
                   crit_b%eqD,crit_b%FluxDMoy,crit_b%GradQDMoy,crit_b%CptbD,grid%ntot)
              end if

!print *, "AP EVAL CRITERES DEF",nrank
              !>-----------TESTS NAN 
              testNaN1=.false.
              testNaN2=.false.
              if (algo_param%Mechanics) then
                     testNaN1=ieee_is_nan(crit_b%Cptb)
                     testNaN2=ieee_is_nan(crit_b%DefMoy)
              end if
              if (algo_param%Diffusion) then
                     do j = 1,nVarD
                        testNaN1 = testNaN1 .or. ieee_is_nan(crit_b%CptbD(j))
                        testNaN2 = testNaN2 .or. ieee_is_nan(crit_b%GradQDmoy(j))
                     end do
              end if
              if(testNaN1) then
                 if(nrank==0) then
                    write (Flog,"(A)")"-------- ERREUR : CRITERE DE COMPATIBILITE NaN-----"
                    if (algo_param%Mechanics) write (Flog,"(A,E15.8)") "critCptb = ", crit_b%eq
                    if (algo_param%Diffusion)write (Flog,"(A,(E15.8))") "critCptbD = ", crit_b%eqD
                 end if
                 call amitex_abort("Critere de compatibilite NaN",2,0)
              end if
              if(testNaN2) then
                 if(nrank==0) then
                    write (Flog,"(A)")&
                         "-------- ERREUR : CRITERE SUR LA DEFORMATION (ou GRADQD) MOYENNE NaN-----"
                    if (algo_param%Mechanics) write (Flog,"(A,E15.8)") "critDefMoy = ", crit_b%DefMoy
                    if (algo_param%Diffusion) write (Flog,"(A,(E15.8))") "critFluxDmoy) = ", crit_b%GradQDmoy
                 end if
                 call amitex_abort("Critere sur la deformation (ou GRADQD) moyenne NaN",2,0)
              end if

!print *, "AV SORTIE SUR CRITERES EN DEF",nrank
              !>-----------SORTIE SUR CRITERES EN DEF (ET GRADQD) (compatibilite & moy) 
              testdef = .true.
              if (algo_param%Mechanics) then
                  testdef = max(crit_b%DefMoy,crit_b%Cptb) .LT. algo_param%tol_critdef
              end if
              if (algo_param%Diffusion) then
                 testdef = testdef & 
                     .and. (max(maxval(crit_b%GradQDMoy),maxval(crit_b%CptbD)) .LT. algo_param%tol_critdef)
              end if
!print *, "testdef",testdef
              if(testdef) then
                 !! si ce critere est respecte, on sort de la boucle
                 !! Avant de sortir, en grandes transformations et si CorrPair=.true.
                 !! et si la resolution en y ou en z est paire, on fait une FFT inverse 
                 !! pour avoir le tenseur de
                 !! Piola-Kirchhoff sans oscillation et on recupere le tenseur des
                 !! contraintes de Cauchy 

                   if (algo_param%Mechanics .AND. algo_param%CorrPair) then
                   if((.NOT. algo_param%HPP) .AND. (modulo(grid%ny,2)==0 .OR. modulo(grid%nz,2)==0)) then
                       call field_ifft(PK1,SigF,algo_param%nTensDef,1)
                       PK1=PK1/ real(grid%ntot,mytype)
                       call PK1ToSig(Sig,Def,PK1)
                   end if 
                   end if
                 testCVFOR=.false.  ! Sortie sur criteres respectes : test convergence forcee => .FALSE.
                 !! on reinitialise le critere min  de la CV forcee (ou sub-stepping) pour une iteration de CVfor  
                 CVFORsauv%critmin = 1.e100_mytype
                 exit
              end if
           end if
           
!print *, "AV EVAL SORTIE CONDITIONNELLE ERREUR",nrank
           !----------------------------------------------------------------------
           !                                          SORTIE CONDITIONNELLE ERREUR 
           !                        nombre d'iterations max atteint sans CV forcee

           ! TODO : calculer detF + sortie si pb, ici ou meme a chaque iteration?...
           
           if(nIt .gt. nIterMax .AND. (.not. algo_param%CVfor) ) then
              !! Avant de quitter le programme :
              !! on sort les images vtk des champs de la derniere iteration effectuee
              if (.not. algo_param%substepping) then ! no write if substepping
              if (.not. algo_param%CVfor  .or. (algo_param%CVfor .and. (nit .eq. algo_param%Nit_cvfor))) then
                 if(nrank == 0) write(Flog,"(A)")"Ecriture des VTK du dernier pas de temps converge."
                 call writeVTK(ind_tps)
              end if
              end if

              if (nrank == 0) then
                 write(OUTPUT_UNIT,"(A)") "Maximum number of iterations reached."
                 write(OUTPUT_UNIT,"(/,A,/)")"-----------------------------"
                 write (Flog,"(A)") "Maximum number of iterations reached."
                 write (Flog,"(A,/)")"-----------------------------"
              end if
              
              !! Si on fait trop d'iterations on arrete le calcul
              !call amitex_abort("Nombre d iterations maximal atteint, le calcul est arrete.",2,0)
              
              !! Si on fait trop d'iterations, on sort de la procedure (pour traitement ulterieur) 
              info_NL_b%reach_nitermax = .true.
              return
 
           end if
           
!print *, "AV EVAL SORTIE CONDITIONNELLE CONVERGENCE FORCEE",nrank
           !----------------------------------------------------------------------
           !                              SORTIE CONDITIONNELLE CONVERGENCE FORCEE 
           !                                           1ERE ITERATION DE CV FORCEE
           !                                              critere nItermax atteint
           !                      RQ: pas de convergence forcee sur le premier pas
           ! 
           if(algo_param%CVfor .AND. (load_incr .ne. 0) .AND. &
              (countCVfor == 0) .AND. (nIt .gt. nIterMax)) then              
              testCVFOR=.true.
              exit
           end if

           !----------------------------------------------------------------------
           !                              SORTIE CONDITIONNELLE CONVERGENCE FORCEE
           !                                        AUTRES ITERATIONS DE CV FORCEE
           !                                             critere Nit_cvfor atteint
           !                     RQ : pas de convergence forcee sur le premier pas
           !  
           if (algo_param%CVfor .and. (nit .eq. algo_param%Nit_cvfor) .and. (load_incr .ne. 0)) then
              if (countCVFOR .eq. (algo_param%Ncvfor-1) ) then
                if(nrank == 0) write(Flog,"(A)")"Ecriture VTK du dernier pas de convergence forcée."
                call writeVTK(ind_tps)
              end if
              testCVFOR=.true.
              exit
           end if
           
!print *, "AV MAJ Eimp",nrank
           !----------------------------------------------------------------------
           !                            MISE A JOUR DE LA DEFORMATION IMPOSEE Eimp

           if(algo_param%Mechanics) then !--MECANIQUE
           if(nb_stress_pilot>0) then 
              call updateEimp(Eimp, local_load%t1(1:2*algo_param%nTensDef),nb_stress_pilot,&
                load(load_n)%DirStress_flag, load(load_n)%DirStress_cauchy, &
                load(load_n)%DirStress,load(load_n)%DirStress2, &
                Matref%C0, SigF, DefF, grid%ntot,fft_start)
           else
              Eimp = local_load%t1(algo_param%nTensDef+1:2*algo_param%nTensDef);
           end if
           end if

           if(algo_param%Diffusion) then !--DIFFUSION
           if(nb_flux_pilot>0) then
              call updateEimpD(EimpD, local_loadD%t1(1:6*algo_param%nVarD),&
                   Matref%K0D, FluxDF, GradQDF, grid%ntot,fft_start)
           else
              EimpD = local_loadD%t1(3*nVarD+1:6*nVarD);
           end if
           end if
           
!print *, "AV GREEN",nrank
           !----------------------------------------------------------------------
           !                                                     APPLICATION GREEN  
           ! format dans Fourier: pinceaux-Z, taille globale ((nx/2)*ny*nz)
           ! ATTENTION : on traite ici les cas des grilles paires
           ! ATTENTION : en simple precision la FFT doit etre realisee sur un champs normalise
           !            (sinon, risque de depassement du reel max) 
           !            sortie de apply_green non normalisee
           ! REMARQUE : SIGF est utilise ici en entree/sortie de appli_Green
           ! 

        if(algo_param%Mechanics) then !----- MECANIQUE
           !! DefF = -\Gamma_0 * (SigF - C_0:DefF)          !modif intel/oneapi (:,:,:,1) -> (i,j,k,1)
!print *, "AV DefF = -\Gamma_0 * (SigF - C_0:DefF)"         !   sinon segfault aleatoires
           do k=fft_start(3),fft_end(3)
           do j=fft_start(2),fft_end(2)
           do i=fft_start(1),fft_end(1)     
           SigF(i,j,k,1) = SigF(i,j,k,1) / real(grid%ntot,mytype) &
                -  lambda*(DefF(i,j,k,1)+DefF(i,j,k,2)+DefF(i,j,k,3))- 2._mytype*mu*DefF(i,j,k,1)
           SigF(i,j,k,2) = SigF(i,j,k,2) / real(grid%ntot,mytype) &
                -  lambda*(DefF(i,j,k,1)+DefF(i,j,k,2)+DefF(i,j,k,3))- 2._mytype*mu*DefF(i,j,k,2)
           SigF(i,j,k,3) = SigF(i,j,k,3) / real(grid%ntot,mytype) &
                -  lambda*(DefF(i,j,k,1)+DefF(i,j,k,2)+DefF(i,j,k,3))- 2._mytype*mu*DefF(i,j,k,3)
           end do
           end do
           end do

           if(algo_param%HPP) then
              do k=fft_start(3),fft_end(3)
              do j=fft_start(2),fft_end(2)
              do i=fft_start(1),fft_end(1)     
              SigF(i,j,k,4:algo_param%nTensDef) = SigF(i,j,k,4:algo_param%nTensDef) &
                   / real(grid%ntot,mytype) - mu*DefF(i,j,k,4:algo_param%nTensDef)
              end do
              end do
              end do
              !! Ici SigF = FFT(sigma - C0:epsilon)
              if (algo_param%HPP_nsym) then
                ! prise en compte du gradient complet du deplacement dans un 
                ! calcul HPP : SigF = sym(DefF_nsym)
                call apply_Green(SigF,FREQ,Matref%LambdaMu0,DefF_nsym)
              else
                call apply_Green(SigF,FREQ,Matref%LambdaMu0)
              end if
           else
              if(algo_param%C0sym) then ! -> version symetrique : C0:gradU_sym
                 do k=fft_start(3),fft_end(3)
                 do j=fft_start(2),fft_end(2)
                 do i=fft_start(1),fft_end(1)     
                 SigF(i,j,k,4:algo_param%nTensDef) = SigF(i,j,k,4:algo_param%nTensDef) &
                    / real(grid%ntot,mytype) &
                    - mu*( DefF(i,j,k,4:algo_param%nTensDef) + DefF(i,j,k,(/7,8,9,4,5,6/)) ) 
                 end do
                 end do
                 end do 
              else                      ! -> version non symetrique : C0:gradU
                 do k=fft_start(3),fft_end(3)
                 do j=fft_start(2),fft_end(2)
                 do i=fft_start(1),fft_end(1)     
                 SigF(i,j,k,4:algo_param%nTensDef) = SigF(i,j,k,4:algo_param%nTensDef) &
                    / real(grid%ntot,mytype) - 2._mytype*mu*DefF(i,j,k,4:algo_param%nTensDef)
                 end do
                 end do
                 end do 
              end if
              !! Ici SigF = FFT(Pi - C0:grad(u))
              call apply_Green_GD(SigF,FREQ,Matref%LambdaMu0)
           end if
           !! En sortie SigF = FFT(-Gamma_0 * Tau)
           !DefF = SigF
           ! modif intel15/oneapi2022 pour cas test fissure_1proc sur marquises machine batch 28p (ok sur machine interactive 14p)
!print *, "AV DefF=SigF"
           do l=1,algo_param%nTensDef
           do k=fft_start(3),fft_end(3)
           do j=fft_start(2),fft_end(2)
           do i=fft_start(1),fft_end(1)     
               DefF(i,j,k,l) = SigF(i,j,k,l)
           end do
           end do
           end do
           end do
!print *, "AP DefF=SigF"
           
           !! On impose la valeur moyenne Eimp a Def en modifiant DefF(0)
           !! si la frequence 0 est dans le processus
           if ((fft_start(1) .eq. 1) .AND. (fft_start(2) .eq. 1) .AND. (fft_start(3) .eq. 1)) then
              DefF(1,1,1,1:algo_param%nTensDef) = cmplx(Eimp,kind=mytype)
           end if
           !! Si HPP et prise en compte complete du gradient du deplacement
           !! On impose la valeur moyenne du gradient du deplacement a Eimp 
           !!  en modifiant DefF_nsym(0) si la frequence 0 est dans le processus
           if (algo_param%HPP_nsym) then
                if ((fft_start(1) .eq. 1) .AND. (fft_start(2) .eq. 1) .AND. (fft_start(3) .eq. 1)) then
                    DefF_nsym(1,1,1,1) = cmplx(Eimp(1),kind=mytype)
                    DefF_nsym(1,1,1,2) = cmplx(Eimp(2),kind=mytype)
                    DefF_nsym(1,1,1,3) = cmplx(Eimp(3),kind=mytype)
                    DefF_nsym(1,1,1,4) = cmplx(Eimp(4),kind=mytype)/2._mytype
                    DefF_nsym(1,1,1,5) = cmplx(Eimp(5),kind=mytype)/2._mytype
                    DefF_nsym(1,1,1,6) = cmplx(Eimp(6),kind=mytype)/2._mytype
                    DefF_nsym(1,1,1,7) = cmplx(Eimp(4),kind=mytype)/2._mytype
                    DefF_nsym(1,1,1,8) = cmplx(Eimp(5),kind=mytype)/2._mytype
                    DefF_nsym(1,1,1,9) = cmplx(Eimp(6),kind=mytype)/2._mytype
                end if
           end if

           !! On annule les fréquences extremes si ny et nz sont pairs
           if(.NOT. algo_param%HPP .AND. algo_param%CorrPair) then 
              if(modulo(grid%ny,2) == 0) then
                 ni2 = grid%ny/2 +1
                 if(ni2>= fft_start(2) .AND. ni2<= fft_end(2) ) then 
                    DefF(fft_start(1):fft_end(1),ni2,fft_start(3):fft_end(3),1:9)=&
                         cmplx(0,kind=mytype)
                 end if
              end if
              if(modulo(grid%nz,2) == 0) then
                 ni2 = grid%nz/2 +1
                 if(ni2>= fft_start(3) .AND. ni2<= fft_end(3) ) then 
                    DefF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),ni2,1:9)=&
                         cmplx(0,kind=mytype)
                 end if
              end if
           end if

        end if !-----FIN MECANIQUE-------------------------------------

        if(algo_param%Diffusion) then !----- DIFFUSION-----------------

           !! DefF = -\Gamma_0 * (SigF - C_0:DefF)
           do i=1,nVarD
           do j=1,3
           do m=fft_start(3),fft_end(3)   ! modif intel15/oneapi2022 (:,:,:,j,i) -> (k,l,m,j,i)
           do l=fft_start(2),fft_end(2)
           do k=fft_start(1),fft_end(1)     

                FluxDF(k,l,m,j,i) =  FluxDF(k,l,m,j,i) / real(grid%ntot,mytype) &
                                 +  Matref%K0D(i)*GradQDF(k,l,m,j,i)
           end do
           end do
           end do
           end do
           end do

           !! Ici SigF = FFT(sigma - C0:epsilon)
           call apply_GreenD(FluxDF, FREQ,Matref%K0D)

           !! En sortie SigF = FFT(-Gamma_0 * Tau)
           GradQDF = FluxDF

           !! On impose la valeur moyenne EimpD a GradQ en modifiant GradQDF(0)
           !! si la frequence 0 est dans le processus
           if ((fft_start(1) .eq. 1) .AND. (fft_start(2) .eq. 1) .AND. (fft_start(3) .eq. 1)) then
              do i = 1,nVarD
                   GradQDF(1,1,1,1:3,i) = cmplx(EimpD(1+3*(i-1):3*i),kind=mytype)
              end do
           end if

        end if !-----FIN DIFFUSION

           !----------------------------------------------------------------------
           !                                                      Def = iFFT(DefF)
           !                                         here Def stores (E + Def_per)
           ! retour espace reel => format (pinceaux-X, taille globale, nx,ny,nz)      

           if(algo_param%Mechanics) then !----- MECANIQUE
               call field_ifft(Def,DefF,algo_param%nTensDef,1)
               if (algo_param%HPP_nsym) then
                ! prise en compte du gradient complet du deplacement dans un 
                ! calcul HPP 
                    call field_ifft(Def_nsym,DefF_nsym,9,1)
               end if
           end if 

           if(algo_param%Diffusion) then !----- DIFFUSION
               call field_ifft(GradQD,GradQDF,3,algo_param%nvarD)
           end if 

           !----------------------------------------------------------------------
           !                           ON AJOUTE LE GRADIENT DE DEFORMATION IMPOSE  
           !                                                           ET DEF_STAR    
           if(algo_param%Mechanics) then !----- MECANIQUE
               !Def = Def + Defimp_hetero
               if (n_gradgradU==27) then
                  call add_gradgradU(Def,algo_param%nTensDef,&
                          local_load%t1(TPRED_fin+1:TPRED_fin+27),load(load_n)%gevolution,&
                           grid%dx,grid%dy,grid%dz,grid%nx,grid%ny,grid%nz)
               end if
               if (add_def_star) Def = Def + Def_star
           end if 

           !----------------------------------------------------------------------
           !                          APPLICATION DE L'ACCELERATION DE CONVERGENCE
           ! lact3 : flag pour signaler que l'ACV vient d'être utilisée 
           ! TODO : dissocier les deux acv

           ! on stocke le residu
           if (algo_param%acc_CV .and. acv_test .and. acv_testD) then 
              if (algo_param%Mechanics) then
                if (algo_param%HPP_nsym) then
                ! prise en compte du gradient complet du deplacement dans un 
                ! calcul HPP : residu calcule sur le gradient de u
                    !ACT3_R(:,:,iact3) = Def_nsym(:,:) - ACT3_U(:,:,iact3)
                    ! modif intel15/oneapi2022 pour cas test fissure_1proc sur marquises machine batch 28p (ok sur machine interactive 14p)
                    do i=1,algo_param%nTensDef
                    do j=1,size(ACT3_R,dim=1)
                        ACT3_R(j,i,iact3) = Def_nsym(j,i) - ACT3_U(j,i,iact3)
                    end do
                    end do
                else 
                    !ACT3_R(:,:,iact3) = Def(:,:) - ACT3_U(:,:,iact3)
                    ! modif intel15/oneapi2022 pour cas test fissure_1proc sur marquises machine batch 28p (ok sur machine interactive 14p)
                    do i=1,algo_param%nTensDef
                    do j=1,size(ACT3_R,dim=1)
                        ACT3_R(j,i,iact3) = Def(j,i) - ACT3_U(j,i,iact3)
                    end do
                    end do
                end if
              end if
              if (algo_param%Diffusion) then
                  do j = 1,nVarD
                  do i = 1,3
                  do l=1,size(ACT3_RD,dim=1)
                     k = i+3*(j-1)
                     ACT3_RD(l,k,iact3) = GradQD(l,i,j) - ACT3_UD(l,k,iact3)                     
                  end do
                  end do
                  end do
              end if
           end if

           ! on rebascule le test ACV a .true.
           if (.not. acv_test) acv_test = .true.
           if (.not. acv_testD) acv_testD = .true.

           ! on applique l'ACV
           lact3 = .false.
           if (algo_param%acc_CV) then 
              if(nIt>=4 .AND. modulo(nIt-4,algo_param%modACV)==0) then 
                 lact3 = .true.
                 if (algo_param%Mechanics) then
                    if (algo_param%HPP_nsym) then     
                         ! prise en compte du gradient complet du deplacement dans un 
                         ! calcul HPP : ACV
                       call act3(ACT3_U,ACT3_R,Def_nsym,GradQD,acv_test)
                       call field_fft(Def_nsym,DefF_nsym,9,1)
                       DefF_nsym=DefF_nsym/real(grid%ntot,mytype)

                       ! Calcul de Def = sym(Def_nsym) correspondant 
                       do j=1,size(Def,dim=1)              ! modif intel/oneapi (:,1) -> (j,1)
                       Def(j,1) = Def_nsym(j,1)
                       Def(j,2) = Def_nsym(j,2)
                       Def(j,3) = Def_nsym(j,3)
                       Def(j,4) = Def_nsym(j,4) + Def_nsym(j,7)
                       Def(j,5) = Def_nsym(j,5) + Def_nsym(j,8)
                       Def(j,6) = Def_nsym(j,6) + Def_nsym(j,9)   
                       end do
                       call field_fft(Def,DefF,algo_param%nTensDef,1)
                       DefF=DefF/ real(grid%ntot,mytype)                                      
                    else
                       call act3(ACT3_U,ACT3_R,Def,GradQD,acv_test)
                       if (n_gradgradU==27) then !si gradgradU, on retire la partie gradient pour le calcul de defF
                          call add_gradgradU(Def,algo_param%nTensDef,&
                           -local_load%t1(TPRED_fin+1:TPRED_fin+27),load(load_n)%gevolution,&
                           grid%dx,grid%dy,grid%dz,grid%nx,grid%ny,grid%nz)
                       end if
                       if (add_def_star) Def = Def - Def_star        
                                      
                       call field_fft(Def,DefF,algo_param%nTensDef,1)
                       
                       DefF=DefF/ real(grid%ntot,mytype) 
                       if (n_gradgradU==27) then
                          call add_gradgradU(Def,algo_param%nTensDef,&
                           local_load%t1(TPRED_fin+1:TPRED_fin+27),load(load_n)%gevolution,&
                           grid%dx,grid%dy,grid%dz,grid%nx,grid%ny,grid%nz)
                       end if       
                       if (add_def_star) Def = Def + Def_star                                                         
                    end if
                 end if
                 if (algo_param%Diffusion) then
                   call act3(ACT3_UD,ACT3_RD,Def,GradQD,acv_testD)
                   call field_fft(GradQD,GradQDF,3,algo_param%nvarD)
                   GradQDF=GradQDF/ real(grid%ntot,mytype)
                 end if
              end if
           end if

           !----------------------------------------------------------------------
           !                                            EVALUATION DU COMPORTEMENT

           if (algo_param%Mechanics) then !-----MECANIQUE
             if (algo_param%HPP_nsym)  then
                ! prise en compte du gradient complet du deplacement dans un 
                ! calcul HPP 
                call behavior(SIG,SIG0,DEF,DEF0,local_load%t0(TPRED_deb:TPRED_fin),&
                           local_load%dt(TPRED_deb:TPRED_fin),local_load%t_load,KSTEP,&
                           nSub_laminate_temp,nIncr_laminate_temp,nIt_laminate_temp,&
                           Def_nsym=Def_nsym,Def_nsym0=Def_nsym0)
             else
                call behavior(SIG,SIG0,DEF,DEF0,local_load%t0(TPRED_deb:TPRED_fin),&
                           local_load%dt(TPRED_deb:TPRED_fin),local_load%t_load,KSTEP,&
                           nSub_laminate_temp,nIncr_laminate_temp,nIt_laminate_temp)
             end if
             if(.not. algo_param%HPP) then
                call SigToPK1(Sig,Def,PK1)
             end if
             if (algo_param%Nloc .and. .false.) then !---- Mecanique non locale resolution implicite (A REVOIR)
                ! TODO revoir possiblite d'introduire de l'implicite
                call Nloc_call(load_n,load_incr,ind_tps,"after",algo_param%nloc_implicit)
             end if
           end if

           if (algo_param%Diffusion) then !-----DIFFUSION
             call behaviorD(FluxD,FluxD0,GradQD,GradQD0,local_loadD%t0(6*algo_param%nVarD+1:),&
                           local_loadD%dt(6*algo_param%nVarD+1:),local_loadD%t_load,KSTEP)
           end if

           VCinfo%nSub_lam_pas  = VCinfo%nSub_lam_pas  + nSub_laminate_temp
           VCinfo%nIncr_lam_pas = VCinfo%nIncr_lam_pas + nIncr_laminate_temp
           VCinfo%nIt_lam_pas   = VCinfo%nIt_lam_pas   + nIt_laminate_temp

           !----------------------------------------------------------------------
           !                          AFFICHAGE fichier.log (MODULO 10 ITERATIONS)

           
           !nIt=nIt+1;      ! incrementation du nombre d'iterations 
             ! A DECOMMENTER pour une sortie des criteres a chaque iteration dans le .log 
             !  if (nrank==0) call write_crit_log(Flog)
             !  if (nrank==0) write (Flog,"(/,A,/)")"-----------------------------"
             ! A DECOMMENTER pour une sortie nIt,crit_b%eq a chaque iteration dans le .log 
             ! if(nrank==0)then
             !     write (Flog,"(I8,E15.8)") nIt,crit_b%eq
             ! end if

!-- sortie energie de deformation
!block
!   real(mytype) :: energy_loc,energy
  
!   energy_loc = sum(SIG*DEF)
!   call MPI_Allreduce(energy_loc,energy,1,MPI_DOUBLE_PRECISION,MPI_SUM,amitex_mpi_comm,ierror)

!   if(nrank==0)then
!     if (nIT .eq. 0 .and.load_incr==0) then
!       open(unit=78941, file=trim(fic_vtk)//"_energy.log",form="formatted",status="replace", action="write")
!     end if
!     write (78941,"(E15.8,I8)") Energy,nIT
!   end if

!end block


           if((nIt .gt. 0)  .and. modulo(nIt,10)==0) then  ! affichage toutes les 10 iterations
              call MPI_BARRIER(amitex_mpi_comm,ierror) !necessite?
              t1 = MPI_WTIME()

              if(nrank==0)then
                 write (Flog,"(/,A,/)")"-----------------------------"
                 if(load_incr==0) then
                    write (Flog,"(A,I8,A)")"pas de chargement fictif ",ind_tps,".5"
                 else
                    write (Flog,"(A,I8)")"pas de chargement ",ind_tps
                 end if
                 write (Flog,"(I8,A)") nIt, " iterations"
                 !sortie fichier log
                 call write_crit_log(Flog)
                 !sortie ecran
                 call write_crit_log(OUTPUT_UNIT)
                 write (Flog,"(A,E15.8)") "Temps total apres initialisation :",t1-times_b%startalgo
                 write (Flog,"(A,E15.8)") "Temps des 10 dernieres iterations :",t1-t10
                 write (Flog, "(A,/,A)") "Temps :",&
                   "FFT_contrainte IFFT_deformation     umat        apply_green   eval_criteres  ecriture"
                 !! A afficher en sortie
                 write(OUTPUT_UNIT,"(/,A,/)")"-----------------------------"
                 if(load_incr==0) then
                    write(OUTPUT_UNIT,"(A,I8,A)")"pas de chargement fictif",ind_tps,".5"
                 else
                    write(OUTPUT_UNIT,"(A,I8)")"pas de chargement ",ind_tps
                 end if
                 write(OUTPUT_UNIT,"(I8,A)") nIt, " iterations"
                 write (Flog,"(6E15.8)") times_f%fft, times_f%ifft, times_m%behavior, times_g%apply, times_g%crit, times_io%wvtk
              end if
              t10 = MPI_WTIME()
           end if

        end do
        !                                                                FIN BOUCLE WHILE
        !================================================================================
        !                                                              CONVERGENCE FORCEE 
        !                                              
        !
        !------------ACTUALISATION FORCEE DES VARIABLES INTERNES
        !
!print *,"AP WHILE - AV CVFOR"
        if (testCVFOR .and. algo_param%CVfor) then

          !--On recalcule le comportement pour evaluer les variables internes 
          !  avec la meilleure solution en def.
          if (algo_param%HPP_nsym)  then
             ! prise en compte du gradient complet du deplacement dans un 
             ! calcul HPP 
             call behavior(SIG,SIG0,CVFORsauv%Def,DEF0,local_load%t0(TPRED_deb:TPRED_fin),&
                        local_load%dt(TPRED_deb:TPRED_fin),local_load%t_load,KSTEP,&
                        nSub_laminate_temp,nIncr_laminate_temp,nIt_laminate_temp,&
                        Def_nsym=Def_nsym,Def_nsym0=Def_nsym0)
          else
             call behavior(SIG,SIG0,CVFORsauv%Def,DEF0,local_load%t0(TPRED_deb:TPRED_fin),&
                        local_load%dt(TPRED_deb:TPRED_fin),local_load%t_load,KSTEP,&
                        nSub_laminate_temp,nIncr_laminate_temp,nIt_laminate_temp)
          end if

          !--On modifie les variables internes et la contrainte en debut de pas de temps
          do i=1,size(MattotP)
             MattotP(i)%VarInt0(1:MattotP(i)%nVarInt,:)=MattotP(i)%VarInt(1:MattotP(i)%nVarInt,:)
          end do  
          Sig0=Sig

         !-- Choix de l'initialisation de la prochaine iteration de CV forcee
         !   pour le cas init_cvfor='best'
         if (algo_param%init_cvfor=='best') then
             !DEF0 = DEF0  not modified
             Def = CVFORsauv%Def
         elseif (algo_param%init_cvfor=='default') then ! done in initDef : DED0 = DEF, DEF = DEF + (DEF-DEF0)*t_ratio) 
         else
             call amitex_abort("UNPAS_NL_BASE : algo_param%init_cvfor different from 'default', 'best'",2,0)
         end if

         !-- Actualisation compteur
         if (nrank==0) write (Flog,"(A,I4,A,E15.8)") &
                       "countCVFOR= ", countCVFOR, " / crit eq = ", CVFORsauv%critmin
         if (nrank==0) write (OUTPUT_UNIT,"(A,I4,A,E15.8)") &
                       "countCVFOR= ", countCVFOR," / crit eq = ", CVFORsauv%critmin
         CVFORsauv%critmin=1.e100_mytype
         countCVFOR = countCVFOR + 1 
        end if

        !
        !------------TEST DE SORTIE DE PROGRAMME SUR LA CONVERGENCE FORCEE
        !
        if (countCVFOR .eq. algo_param%Ncvfor) then 
           if (nrank == 0) then
              write (OUTPUT_UNIT,"(A,I8,A)") "Plus de ",algo_param%Ncvfor," convergences forcees, le calcul est arrete."
              write (OUTPUT_UNIT,"(/,A,/)")"-----------------------------"
              write (Flog,"(A,I8,A)") "Plus de ",algo_param%Ncvfor," convergences forcees, le calcul est arrete."
              write (Flog,"(/,A,/)")"-----------------------------"
           end if
           call amitex_abort("Nombre maximal de convergences forcees atteint",2,0)
        end if

        end do
        !                                                   FIN BOUCLE CONVERGENCE FORCEE
        !================================================================================

  call MPI_BARRIER(amitex_mpi_comm,ierror)
  times_b%unpas = times_b%unpas + MPI_WTIME()- t_ini

  info_NL_b%nit_unpas  = Nit

!print *,"FIN UNPAS_NL_BASE"
end subroutine unpas_NL_base



!==================================================================================
!==================================================================================
!                         SUBROUTINE LOG_AFTER_UNPAS
!
!> Write informations after a loading step in resolution_NL_base
!!
!!
!==================================================================================
subroutine log_after_unpas(load_incr,ind_tps, nIt,&
               t_fft0,t_ifft0,t_behavior0,t_green0,t_crit0,t_wvtk0,t_wstd0,t_unpas0,&
               testLaminate,msg)
  
  implicit none
  
  integer, intent(in)      :: load_incr,ind_tps, nIt
  real(mytype), intent(in) :: t_fft0,t_ifft0,t_behavior0,t_green0,&
                              t_crit0,t_wvtk0,t_wstd0,t_unpas0
  logical, intent(in)      :: testLaminate
  character(len=*),optional,intent(in) :: msg
  integer                  :: ierror  ! mpi error
  double precision         :: t1      ! time measurement


  call MPI_BARRIER(amitex_mpi_comm,ierror)
  t1 = MPI_WTIME()

  if (nrank==0) then
    if (present(msg)) write(OUTPUT_UNIT,"(A)") msg
    write(OUTPUT_UNIT,"(/,A)")     "------------------------------------------ ("//trim(simu_name)//")"
    write(OUTPUT_UNIT,"(A,I8)")"Nombre d'iterations     :",nIt
    write(OUTPUT_UNIT,"(A,/,/)")"=================================================="
    !! Dans le log
    write(Flog,"(/,A)")"--------------------------------------------------"
    if(load_incr==0) then
       write(Flog,"(A,I8,A)")"Fin du pas de chargement fictif :",ind_tps,".5"
    else
       write(Flog,"(A,I8)")"Fin du pas de chargement :",ind_tps
    end if
    write(Flog,"(A,I8)")"Nombre d'iterations     :",nIt 
    call write_crit_log(Flog)

    !! Info umatLaminate
    if (testLaminate) then
       write(Flog,"(A)") "UmatLaminate : "
       write(Flog,"(A,F12.4)")"          - Nombre d'iterations moyen :",&
                   real(VCinfo%nIt_lam_pas)/real(VCinfo%nCall_lam_pas)
       !! Info éventuels pilotages effectués pour umatLaminate
       if (VCinfo%nSub_lam_pas > 0) then
          write(Flog,"(A,I12,A,I12,A,F12.4,A)")"          - Nombre de recours à la subdivision du & 
            & pas de temps/Nombre d'appels  :",VCinfo%nSub_lam_pas,&
              " / ",VCinfo%nCall_lam_pas,"  (",100*real(VCinfo%nSub_lam_pas)/real(VCinfo%nCall_lam_pas),"%)"
          write(Flog,"(A,F12.4)")"          - Nombre moyen de sous pas de temps utilisé pour les &
            & subdivisons      :", real(VCinfo%nIncr_lam_pas)/real(VCinfo%nSub_lam_pas)
       end if
    end if
    
    write(Flog,"(A,E15.8)") "Total time             :",t1-times_b%startalgo
    write(Flog, "(A,/,A)") &
      "FFT             IFFT           umat           apply_green    eval_crit      write VTK      write STD",&
      "Total times :"
    t1 = MPI_Wtime()
    write(Flog,"(7E15.8)") times_f%fft, times_f%ifft, times_m%behavior, times_g%apply,&
                           times_g%crit, times_io%wvtk, times_s%wtot
    write(Flog,"(A,/,7E15.8)")"Loading step times:",times_f%fft - t_fft0, &
         times_f%ifft - t_ifft0, times_m%behavior - t_behavior0, times_g%apply - t_green0,&
         times_g%crit - t_crit0, times_io%wvtk - t_wvtk0, times_s%wtot - t_wstd0
    write(Flog,"(A,E15.8,A,E15.8,A)")"Time in loading step time (and in unpas_NL_base) :",&
         t1-times_b%startstep,"(",times_b%unpas - t_unpas0,")"
    write(Flog,"(A)")"=================================================="
    write(Flog,"(A,/,/)")"=================================================="
    
  end if

end subroutine log_after_unpas

!==================================================================================
!==================================================================================
!                         SUBROUTINE LOG_FINAL_TIMES
!
!> Write final times informations after resolution_NL_base
!!
!!
!==================================================================================
subroutine log_final_times(nitTot)

  implicit none
 
  integer, intent(in)      :: nitTot
  integer                  :: ierror  ! mpi error
  double precision         :: t1      ! time measurement
  

  call MPI_BARRIER(amitex_mpi_comm,ierror)
  t1 = MPI_WTIME()
  
  if (nrank==0) then
     write(Flog,"(/,A)")     "========================================== ("//trim(simu_name)//")"
     write(Flog,"(A)")       "Algorithm times"
     write(Flog,"(A)")       "---------------"    
     write(Flog,"(A,E15.8)") "Cumulated time for  behavior      : ", times_m%behavior
     write(Flog,"(A,E15.8)") "Cumulated time for  eval_crit     : ", times_g%crit
     write(Flog,"(A,E15.8)") "Cumulated time for  apply_green   : ", times_g%apply
     write(Flog,"(A,E15.8)") "Cumulated time for  FFT           : ", times_f%fft
     write(Flog,"(A,E15.8)") "Cumulated time for  iFFT          : ", times_f%ifft
     write(Flog,"(A,E15.8)") "Cumulated time for  write VTK     : ", times_io%wvtk
     write(Flog,"(A,E15.8)") "Cumulated time for  write STD     : ", times_s%wtot
     write(Flog,"(A)")       "-----------------------------------------------------"    
     write(Flog,"(A,E15.8)") "Sum                               : ",&
          times_m%behavior + times_g%crit + times_g%apply + times_f%fft &
          + times_f%ifft + times_io%wvtk + times_s%wtot
     write(Flog,"(A,E15.8)") "Cumulated time in unpas_NL_base   : ", times_b%unpas
     write(Flog,"(A,E15.8)") "Total time in resolution_NL_base  : ", t1-times_b%startalgo

     write(Flog,"(A,/, I12)") "nombre total d'iterations" ,  nitTot ! en francais, utilise dans validation/comparaison.f90

     write(OUTPUT_UNIT,"(/,A)")     "========================================== ("//trim(simu_name)//")"
     write(OUTPUT_UNIT,"(/,A)")     "Algorithm times"
     write(OUTPUT_UNIT,"(A)")       "---------------"    
     write(OUTPUT_UNIT,"(A,E15.8)") "Cumulated time for  behavior      : ", times_m%behavior
     write(OUTPUT_UNIT,"(A,E15.8)") "Cumulated time for  eval_crit     : ", times_g%crit
     write(OUTPUT_UNIT,"(A,E15.8)") "Cumulated time for  apply_green   : ", times_g%apply
     write(OUTPUT_UNIT,"(A,E15.8)") "Cumulated time for  FFT           : ", times_f%fft
     write(OUTPUT_UNIT,"(A,E15.8)") "Cumulated time for  iFFT          : ", times_f%ifft
     write(OUTPUT_UNIT,"(A,E15.8)") "Cumulated time for  write VTK     : ", times_io%wvtk
     write(OUTPUT_UNIT,"(A,E15.8)") "Cumulated time for  write STD     : ", times_s%wtot
     write(OUTPUT_UNIT,"(A)")       "-----------------------------------------------------"    
     write(OUTPUT_UNIT,"(A,E15.8)") "Sum                               : ",&
          times_m%behavior + times_g%crit + times_g%apply + times_f%fft &
          + times_f%ifft + times_io%wvtk + times_s%wtot
     write(OUTPUT_UNIT,"(A,E15.8)") "Cumulated time in unpas_NL_base   : ", times_b%unpas
     write(OUTPUT_UNIT,"(A,E15.8)") "Total time in resolution_NL_base  : ", t1-times_b%startalgo
     
     write(OUTPUT_UNIT,"(/,A,/, I12,/)") "Total number of iterations :" ,  nitTot
  end if

end subroutine log_final_times

!==============================================================================
!      SUBROUTINE WRITE_CRIT_LOG
!
!>  ECRIURE DES CRITERES DANS LE FICHIER LOG
!!      si la diffusion est active, on affiche les valeurs max pour chaque variable
!!      sinon, on affiche 0
!!
!!  crit : type(CRIT_BASE) variable de ce module (NL_base_mod.f90)
!!
!!  \param[in]    FU, unite logique (Flog ou OUTPUT_UNIT)
!!
!------------------------------------------------------------------------------
subroutine write_crit_log(FU)

   implicit none 

   integer, intent(in)                  :: FU

   real(mytype)                         :: A,B,C,D
   character(len=5),parameter           :: FMT_real="E15.8"             !<format d'ecriture d'un reel
!   character(len=6),parameter          :: FMT_real="E25.18"            !<format d'ecriture d'un reel tres precis


   !> Si la diffusion n'est pas activee, les criteres associes sont affiches a 0.
   if (algo_param%Diffusion) then   
      A = maxval(crit_b%eqD)
      B = maxval(crit_b%FluxDMoy)
      C = maxval(crit_b%CptbD) 
      D = maxval(crit_b%GradQDMoy)
   else  
      A = 0.
      B = 0.
      C = 0. 
      D = 0.
   end if

   !> ecriture
   write(FU,"(A)")     "------------------------------------------ ("//trim(simu_name)//")"                               
   write(FU,"(2(A,"//FMT_real//"))") &
                      "criteres d'equilibre                (meca,diff) :", crit_b%eq,' , ', A
   write(FU,"(2(A,"//FMT_real//"))") &
                      "criteres sur la contrainte moyenne  (meca,diff) :", crit_b%SigMoy,' , ', B
   write(FU,"(2(A,"//FMT_real//"))") &
                      "criteres de compatibilite           (meca,diff) :", crit_b%Cptb,' , ', C
   write(FU,"(2(A,"//FMT_real//"))") &
                      "criteres sur la deformation moyenne (meca,diff) :", crit_b%DefMoy,' , ', D

end subroutine write_crit_log
!==============================================================================


end module NL_base_mod
