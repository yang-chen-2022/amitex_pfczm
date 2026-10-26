!===============================================================================
!>                  PROGRAMME PRINCIPAL DU CODE AMITEX_FFTP
!! 
!! LE PROGRAMME PRINCIPAL SE CHARGE :
!!      1/ D'allouer et d'initialiser, a partir des fichiers d'entree,
!!            l'ensemble des donnees necessaires au calcul dans SIMU(1)
!!      2/ D'associer les pointeurs correspondant aux champs de SIMU(1) 
!!      3/ De lancer le calcul (resolution_NL_base par exemple) 
!!            utilisant ces pointeurs
!!      
!! NOTATIONS (valable pour les entrees du code et les champs) :
!!            Sigma et Def (HPP) : CAST3M (Voigt + ordre 11 22 33 12 13 23)
!!            Pi et grad(u) (GDEF) (ordre 11 22 33 12 13 23 21 31 32 sans notation)
!!
!!
!! REGLES DE PROGRAMMATION :utiliser les types definis dans 2decomp 
!!            supprimer les double precision => real(mytype), 
!!                  sauf les double precision sortie de la librairie MPI 
!!                  (ex. : temps sortie de MPI_WTIME())
!!            dans les fonctions MPI, remplacer les types MPI : 
!!                  MPI_DOUBLE_PRECISION=> real_mytype 
!!            creation de complexe : 1 + i = cmplx(1,1,kind=mytype)
!!            attention aux constantes => pi = 3.14.._mytype
!!            float => real(nx,mytype)
!!            Pour passer les transformees de Fourier DefF et SigF aux autres fonctions
!!            ne pas oublier de preciser que les tableaux vont de fft_start à fft_end
!!            sinon Fortran les alloue de 1 à fft_end-fft_start + 1
!!
!! RESTART :  For saving files in a restart folder : use <Restart Every="10" /> 
!!                in partial loadings (see loading.xml)
!!            As soon as the restart folder exists, it is used.
!!            Only available for mechanics.
!!            Only available using the same decomposition (prow,pcol)
!!
!!
!! ATTENTION : en simple precision la FFT doit etre realisee sur un champs normalise
!!            (sinon, risque de depassement du reel max) 
!!
!!==============================================================================
!!
!! UTILISATION DE L'EXECUTABLE :
!!      amitex_fftp -help
!!
!!==============================================================================
program amitex_fftp

!> CHARGEMENT DES DIFFERENTS MODULES

  use ISO_FORTRAN_ENV

  use MPI             
#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi,       only : nrank
  use decomp_2d ,          only : decomp_2d_finalize
  use decomp_2d_fft,       only : decomp_2d_fft_finalize,decomp_2d_fft_set_ngrid
#else
  use decomp_2d,      only : mytype, nrank, decomp_2d_finalize
  use decomp_2d_fft,  only : decomp_2d_fft_finalize  
#endif  

  use io_amitex_mod,  only : write_stdout0 ! function
  use io2_amitex_mod, only : print_coeff_varint_vtk ! function
  use param_algo_mod, only : algo_param, user_param ! pointers
  use material_mod ,  only : VCInfo, testComposite  ! pointers
  use loading_mod,    only : load ! pointer
  use error_mod,      only : amitex_abort !function
  use resolution_mod, only : resolution_NL_base ! function
  use amitex_mod,     only : Sig, Flog, &       ! pointers
                             AMITEXenv, LOGO    ! variables
  use field_mod,      only : field_isNaN_real,& ! function 
                             ntotFP, FREQ, amx_size  ! pointer
  use non_local_mod,  only : init_user_nloc_variables_call ! function
  use simu_mod,       only : SIMU, & ! variable
                             nullify_pointers_simu, init_simu_command_line,&
                             associate_pointers_simu, simu_restart_read !functions
  use amitex_mpi_mod, only: amitex_mpi_comm

 !! USERS
  use amitex_user_mod,     only : init_user_variables  ! function
  use resolution_user_mod, only : resolution_user      ! function

!!==============================================================================
!!                                                                  DECLARATIONS
!!==============================================================================
  implicit none

!!------------------------------------------------------------------------------
!>                                                                    SIMULATION  
  integer                                           :: choix_simu=1 
                                        !< 1 : 1 simu standard 
                                        !! 2 : (dev) 2 simu, on effectue la 2nde

!!------------------------------------------------------------------------------
!>                                                                  PARALLELISME 

  integer :: ierror                     !< erreur relative au fonction MPI


!!------------------------------------------------------------------------------
!>                                                          SUIVI DE L'EXECUTION 

  double precision :: t0,t1,t2          !< double precision au lieu de real_mytype requis par MPI_WTIME
                                        !! t0 : date de depart de l'algorithme
                                        !! t1 : temps pour evaluer la duree de certaines fonctions
                                        !! t2 : date de debut du programme (pour duree de l'initialisation)

!!------------------------------------------------------------------------------
!>                                                                 VERIFICATIONS

  logical :: resolution=.true.          !< lancement de la resolution
  logical :: check                      !< verifications avant la resolution
  integer :: coeff2print =0             !< indice du coeff a sortir (ligne de commande) 
  integer :: varint2print=0             !< indice de varint a sortir (ligne de commande) 
  logical :: help=.false.               !< test pour avoir l'aide et sortir du code

!!------------------------------------------------------------------------------
!>                                                                        DIVERS

  integer            :: alloc_stat      !< erreur lors d'une allocation memoire
  character(len=200) :: version         !< numero de version ecrit dans .log
  logical            :: bool0 

  t0 = 0 ! avoid gcc-warning (unused var)
  
!!==============================================================================
!>                                                            INITIALISATION MPI 
!>                                                        + INITIALISATION nrank
!>                            -> necessaire tant que 2decomp n'est pas inialisee
!!==============================================================================
  call MPI_INIT(ierror)
  t2 = MPI_WTIME()
  call MPI_COMM_RANK(amitex_mpi_comm,nrank,ierror)

  ! write logo on stdout
  call write_stdout0(achar(10)//LOGO//achar(10))

  ! write version on stdout
  version="x.y.z" ! do not modify this line
  call write_stdout0("VERSION "//trim(version))

!!==============================================================================
!!                                                     INITIALISATION DE LA SIMU
!!==============================================================================

  if (choix_simu == 1) allocate(SIMU(1))
  if (choix_simu == 2) allocate(SIMU(2))

#ifdef NEW_2DECOMP  
  call decomp_2d_fft_set_ngrid(size(SIMU))
#endif

  call nullify_pointers_simu()  
  call init_simu_command_line(SIMU(1),version, help,coeff2print,varint2print, Igrid = 1)

  ! fin de programme si "amitex_fftp -help"
  if (help) then
     goto 666
  end if

  if (choix_simu == 2) then
     call nullify_pointers_simu()
     call init_simu_command_line(SIMU(2),version, help,coeff2print,varint2print,&
                 "cmd_files/bilayer_2mat_2zone_xmla.in", Igrid = 1)
  end if 
  
!!==============================================================================
!!                                                                    CHOIX SIMU 
!!==============================================================================
  call nullify_pointers_simu()

  select case(choix_simu)
  case(1)
     call associate_pointers_simu(SIMU,1)
  case(2)
     call associate_pointers_simu(SIMU,2)   
  end select
  
!!==============================================================================
!!                                               CHECK AVANT LANCEMENT DU CALCUL
!!                            TODO : creer une fonction check_simu dans simu_mod
!!==============================================================================
! possibilite de supprimer tous les checks
  check=.true.

  if (check) then

  ! check NaN Freq
  bool0=.false.
  call field_isNaN_real(Freq,3*ntotFP,bool0)
  if (bool0) call amitex_abort("Presence de NaN dans Freq (amitex)",2)

  ! check "Voxels Composites" ET "Modeles Non-Locaux"
  if (testComposite .AND. algo_param%Nloc) then 
      call amitex_abort("Utilisation simultanee de modeles non-locaux &
              & et de voxels composites en cours de developpement (amitex)",0,0)
  end if

  ! check "utilisateur (algo complet ou fonctions)" et "Modeles non-locaux"
  if (user_param%test .AND. algo_param%Nloc) then 
      call amitex_abort("Utilisation simultanee de modeles non-locaux &
              & et d'un noeud <User> dans le fichier algorithme (.xml)  (amitex)",2)
  end if

  end if
  call write_stdout0("After a few check") 
  
!!==============================================================================
!!                                                         SORTIE VTK POUR VERIF 
!!                                   (option -coeff2print i et/ou varin2print i)
!!                              on utilise Sig pour construire le champ de coeff 
!!                                                mais on ne lance pas le calcul
!!==============================================================================
  ! En diffusion pure, Sig n'est pas alloué => on alloue un champ a une composante
  if(.not. allocated(SIMU(1)%Sig)) then    
        allocate(SIMU(1)%Sig(amx_size(1)*amx_size(2)*amx_size(3),1),stat=alloc_stat)
        if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex13)",2)
        Sig => SIMU(1)%Sig
  end if

  call print_coeff_varint_vtk(coeff2print,varint2print)
  
  if ((coeff2print > 0) .OR. (varint2print > 0)) then
    resolution=.false.
  end if 

!!==============================================================================
!!                                                           INITIALIZATION USER
!!==============================================================================

if (user_param%test) then
  call write_stdout0("Before init_user_variable") 
  call init_user_variables()
  call write_stdout0("After init_user_variable") 
end if

!!==============================================================================
!!                                                      INITIALIZATION NLOC USER
!!      allocation and initialization of additional variables for user_nloc1,2,3
!!                                                                  if necessary
!!==============================================================================

if (algo_param%Nloc) then
    call write_stdout0("Before init_user_nloc_variables_call") 
    call init_user_nloc_variables_call()
    call write_stdout0("After init_user_nloc_variables_call") 
end if

!!==============================================================================
!!                                                                  READ RESTART
!!                                      and define : load(i)%skip_before_restart
!!                                                   load(iload)%restart_inc
!!==============================================================================
block
integer :: iload, iincr, i 

  call simu_restart_read(iload,iincr,read_option=0)       ! read restart in infos.txt

  if (iload /=0) then ! no restart if no files or pb reading restart files
    do i = 1, iload-1
      load(i)%skip_before_restart = .true.
    end do    
    if (iincr<=load(iload)%Nincr) then
      load(iload)%restart_inc = iincr
      call simu_restart_read(iload,iincr,read_option=1) ! read and update restart files
    else
      call amitex_abort("Restart issue : restart increment > Nb incr of the partial loading",2,0)
    end if  
  end if

end block 

!!==============================================================================
!!                                                                     ALGO USER
!!==============================================================================

if (resolution .AND. user_param%test .AND. trim(user_param%algo) == "user") then

  !> temps total d'initialisation
  call MPI_BARRIER(amitex_mpi_comm,ierror) 
  t0 = MPI_WTIME()
  if(nrank==0) then
      write(Flog,fmt="(A,E15.8)") "Duration TOTAL INITIALIZATION (s) : ", t0 - t2
      write(OUTPUT_UNIT,fmt="(A,E15.8)") "Duration TOTAL INITIALIZATION (s) : ", t0 - t2
  end if

  !> Lancement de l'algorithme USER
  call resolution_user()

end if


!!==============================================================================
!!                                                                 ALGO STANDARD
!!==============================================================================
!!
!! POUR LE CALCUL: L'ALGORITHME PEUT UTILISER, EN PLUS DES PARAMETRES TRANSMIS,
!!            LES STRUCTURES SUIVANTES AFFECTEES DURANT LES INITIALISATIONS :
!!      
!!      MattotP         tableau de structure MATERIAL   (material_mod.f90)
!!              -> description du materiau
!!
!!      load            tableau de structure LOADING    (loading_mod.f90))
!!              -> description du chargement
!!
!!      initValExt      structure INITLOADEXT           (loading_mod.f90)
!!              -> description de l'initialisation du chargement externe
!!
!!      extract         structure PARAMEXTRACT          (loading_mod.f90)
!!              -> description des sorties vtk
!!
!!      algo_param      structure PARAM_ALGO            (param_algo_mod.f90)
!!              -> description des parametres algorithmiques
!!
!!      grid            structure GRID_DIM              (green_mod.f90)
!!              -> description de la grille 
!!
!!==============================================================================
call write_stdout0("Just before standard algorithm") 
if(resolution .AND. .NOT. (user_param%test .AND. trim(user_param%algo) == "user" )) then

  !> temps total d'initialisation
  call MPI_BARRIER(amitex_mpi_comm,ierror) 
  t0 = MPI_WTIME()
  if(nrank==0) then
      write(Flog,fmt="(A,E15.8)") "Duration TOTAL INITAIALIZATION (s) : ", t0 - t2
      write(OUTPUT_UNIT,fmt="(A,E15.8)") "Duration TOTAL INITAIALIZATION (s) : ", t0 - t2
  end if

  !> Lancement de l'algorithme
  if (algo_param%scheme_type == "basic_scheme") then
     call resolution_NL_base()
  end if

  if(nrank==0) then
     if (testComposite) then
        write(Flog,"(A)") "------------------------------"
        write(Flog,"(A)") "UmatLaminate "
        write(Flog,"(A,F12.4)") "          - nombre d'iterations moyen par voxel algo laminate :" ,&
                    real(VCinfo%nIt_lam)/real(VCinfo%nCall_lam)
        write(Flog,"(A,I12,A,I12,A,F12.4,A)") &
          "          - Nombre de recours à la subdivision du pas de temps/Nombre d'appels  :",VCinfo%nSub_lam,&
          " / ",VCinfo%nCall_lam,"  (",100*real(VCinfo%nSub_lam)/real(VCinfo%nCall_lam),"%)"
        if (VCinfo%nSub_lam >0) write(Flog,"(A,F12.4)") &
          "          - Nombre moyen subdivisions dans PilotageLaminate  :",&
          real(VCinfo%nIncr_lam)/real(VCinfo%nSub_lam)
        write(Flog,"(A,/)") "------------------------------"
     end if
  end if

end if


!!==============================================================================
!!                                            DESALLOCATIONS ET FIN DU PROGRAMME
!!==============================================================================
  call write_stdout0("Deallocations and end of program")

  if(allocated(SIMU)) deallocate(SIMU)

  if(allocated(AMITEXenv%path)) deallocate(AMITEXenv%path)


  call decomp_2d_fft_finalize

  call MPI_BARRIER(amitex_mpi_comm,ierror) 
  t1 = MPI_WTIME()
  if (nrank==0) then
    write(Flog,"(A,E15.8)") 'temps global (hors initialisations) =', t1-t0
    write(OUTPUT_UNIT,"(A,E15.8)") 'temps global (hors initialisations) =', t1-t0
    close(Flog)
  end if
  call decomp_2d_finalize !must be after the last usage of nrank!
666  call MPI_FINALIZE(ierror)
end program amitex_fftp

!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
