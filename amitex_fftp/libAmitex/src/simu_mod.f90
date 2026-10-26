!123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!===================================================================================================
!
!       MODULE SIMU_MOD 
! 
!> MANAGEMENT OF AMITEX SIMULATIONS : INITIALIZING, CHANGING SIMULATIONS
!!
!!      Simulations are stored in a 1D array of SIMU_AMITEX objects
!!      Access to simulations is done via pointers towards SIMU_AMITEX components
!!
!!
!!  Objets publiques : type SIMU_AMITEX
!!
!!  Subroutines (publiques) : 
!!      init_simu_command_line  : affecte un objet SIMU_AMITEX a partir de la ligne de commande 
!!                                                                      (simulation principale)
!!                                      et/ou d'un fichier regroupant ces memes infos (NAMELIST)
!!                                associe les pointeurs correspondants
!!                                desalloue les donnees intermediaires (chargees a la lecture)
!!      nullify_pointers_simu   : mise a l'etat nul de tous les pointeurs
!!      associate_pointers_simu : association de tous les pointeurs a une simulation
!!
!===================================================================================================

module simu_mod

  use ISO_FORTRAN_ENV
  use iso_c_binding

  use mpi
#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype,& ! TYPE
                                  PHYSICAL_IN_X ! variable
  use decomp_2d_mpi,       only : nrank  ! variable
  use decomp_2d,           only : DECOMP_INFO,& !TYPE
                                  !xstart, xend, xsize,&  !variable 
                                  !nx_global,ny_global, nz_global, & ! variables 
                                  decomp_2d_init, decomp_2d_finalize   ! functions
  use decomp_2d_fft,  only : decomp_2d_fft_init,decomp_2d_fft_finalize,decomp_2d_fft_get_size,& !functions
                             decomp_2d_fft_get_sp,decomp_2d_fft_get_ph,decomp_2d_fft_set_ngrid,&
                             decomp_2d_fft_get_ngrid,decomp_2d_fft_use_grid 
#else
  use decomp_2d,      only : mytype, DECOMP_INFO, & !TYPES
                             nrank,&                ! variables
                             !xsize, xstart, xend,&  ! pointer variables (with modified decomp_2d for multigrids)
                             !nx_global, ny_global, nz_global,&  ! pointer variables (with modified decomp_2d for multigrids)
                             decomp_2d_init, decomp_2d_finalize   ! functions
                             
  use decomp_2d_fft,  only : PHYSICAL_IN_X,& ! variables
                             decomp_2d_fft_init,decomp_2d_fft_finalize,decomp_2d_fft_get_size,& !functions
                             associate_pointers_decomp_2d_fft, get_decomp_fft_info
#endif                             
  use amitex_mod,     only : SAUVCVFOR,NLOC_FIELD,GNLOC_FIELD,& ! TYPES
                             fic_vtk0,fic_log0,Flog0,& ! init variables
                             fic_vtk,fic_log,Flog,simu_name,& ! pointers variables
                             Sig,Sig0,Def,Def0,Def_nsym,Def_nsym0,PK1,Def_star,Dep,SigF,DefF,DefF_nsym,depF,& ! pointers variables
                             FluxD,FluxD0, GradQD,gradQD0,FluxDF, GradQDF,& ! pointers variables
                             ACT3_R,ACT3_U,ACT3_RD,ACT3_UD,& ! pointers variables
                             CVFORsauv,FvolN,FvolNF,TEMPfield32,Nloc,GNloc,& ! pointers variables
                             copyXML,lire_commande,init_log, get_AMITEXenv,set_prow_pcol,& ! functions
                             test_FvolN !provisoire
  use param_algo_mod, only : PARAM_ALGO,PARAM_ALGO_USER,PARAM_ALGO_NLOC,& !TYPE
                             algo_param0,user_param0,nloc_param0,& ! init variables
                             algo_param,user_param,nloc_param,& ! pointers
                             read_param, read_param_composite,print_param_algo,& !functions
                             print_param_algo_composite !functions
  use material_mod,   only : MATERIAL, NMATERIAL, MATERIAL_REFERENCE, INFO_VOXCOMP,TIMES_MAT,NL_MODEL,& !TYPES
                             MattotP0,nmateriaux0,matref0,Nloc_models0,MatComposite,Interphases,& !init variables
                             nmateriaux, Matref, MattotP, VCinfo, testComposite,times_m,Nloc_models,& ! pointers
                             test_composite,read_composite,read_mat,eval_nb_materiaux,& !functions
                             print_material_structure,deallocate_MatComposite,& !functions
                             read_mat_composite,get_interphase_list,deallocate_Interphases,& !functions
                             deallocate_Matref0,deallocate_mattotP0 ! functions
  use green_mod,      only : TIMES_GREEN,& !TYPES
                             times_g       ! Pointer variable                       
  use loading_mod,    only : LOADING,INITLOADEXT,LOCALLOAD, PARAM_EXTRACT,& !TYPES
                             initValExt0, local_load0, local_loadd0, n_gradgradU0, add_def_star0, extract0,load0,& !init variables
                             load, initValExt,extract,local_load, local_loadD,n_gradgradU,add_def_star, & !pointers
                             read_load, print_loading,print_initloadext,print_param_extract,&  !functions
                             deallocate_load, deallocate_initvalext, deallocate_local_load, deallocate_param_extract0 !functions
  use sortie_std_mod, only : BUFFER_STD, STRUCT_MAT_SORTIE, TIMES_STD,& !types 
                             buffer0,Matcomp_sortie0,& !init variables
                             buffer,Matcomp_sortie,times_s,& !pointers
                             init_buffers, initsortie_std,init_matcomp_sortie, desallocation_sortie_std !functions
  use field_mod,      only : ph_decomp, sp_decomp, fft_start, fft_end, fft_size, ntotP, ntotFP,times_f,& ! pointer variables
                             amx_start, amx_end, amx_size,&
                             FREQ,FREQ_2,grid,&
                             GRID_DIM, TIMES_FIELD,& !TYPES
                             initGrid, print_grid, initFreq,initFreq_2 ! functions
  use io_amitex_mod,  only : write_stdout0,read_header_vtk,write_file0 !functions
  use io2_amitex_mod, only : times_io,& ! pointer variables
                             TIMES_IO2 !TYPES
  use non_local_mod,  only : init_Nloc_Variables !functions
  use error_mod,      only : error,& ! pointer
                             check_amitex_abort, amitex_abort !functions
  use NL_base_mod,    only : CRIT_BASE,TIMES_BASE,INFO_NL_BASE,& !TYPES
                             crit_b,times_b,info_NL_b !pointers
  use amitex_mpi_mod, only: amitex_mpi_comm

  implicit none

  private
  

  !> Variables publiques
  public :: SIMU, Iactive_simu
  
  !> Pointer "publiques" vers composantes de SIMU_AMITEX
  public :: Igrid_simu
  
  !> Type publique 
  public :: SIMU_AMITEX

  !> Fonctions publiques
  public :: init_simu_command_line, nullify_pointers_simu, associate_pointers_simu,&
            simu_restart_write, simu_restart_read


  !> Private function  : TODO - move elsewhere

  interface
    function mkdir(path,mode) bind(c,name="mkdir")
      use iso_c_binding
      integer(c_int) :: mkdir
      character(kind=c_char,len=1) :: path(*)
      integer(c_int16_t), value :: mode
    end function mkdir
  end interface


!!------------------------------------------------------------------------------
!>                                                                    SIMULATION
  type SIMU_AMITEX
        !---amitex_mod.f90 - fields (voir description in amitex_mod.f90)
        real(mytype),allocatable,dimension(:,:)          :: Sig,Sig0 
        real(mytype),allocatable,dimension(:,:)          :: Def,Def0
        real(mytype),allocatable,dimension(:,:)          :: Def_nsym,Def_nsym0                                                  
        real(mytype),allocatable,dimension(:,:)          :: PK1
        real(mytype),allocatable,dimension(:,:)          :: Def_star
        real(mytype),allocatable,dimension(:,:)          :: Dep
        real(mytype),allocatable,dimension(:,:,:)        :: FluxD,FluxD0  
        real(mytype),allocatable,dimension(:,:,:)        :: GradQD,GradQD0
        complex(mytype),allocatable,dimension(:,:,:,:)   :: SigF 
        complex(mytype),allocatable,dimension(:,:,:,:)   :: DefF
        complex(mytype),allocatable,dimension(:,:,:,:)   :: DefF_nsym
        complex(mytype),allocatable,dimension(:,:,:,:)   :: DepF
        complex(mytype),allocatable,dimension(:,:,:,:,:) :: FluxDF
        complex(mytype),allocatable,dimension(:,:,:,:,:) :: GradQDF
        real(mytype),allocatable,dimension(:,:)          :: FvolN 
        complex(mytype),allocatable,dimension(:,:,:,:)   :: FvolNF
        real(mytype),allocatable,dimension(:,:,:)        :: ACT3_R, ACT3_U
        real(mytype),allocatable,dimension(:,:,:)        :: ACT3_RD, ACT3_UD
        type(sauvCVFOR)                                  :: CVFORsauv
        real(kind=REAL32),allocatable,dimension(:)       :: TEMPfield32   
        type(NLOC_FIELD),allocatable,dimension(:)        :: Nloc
        type(GNLOC_FIELD),allocatable,dimension(:)       :: GNloc
        !---amitex_mod.f90 - output files
        integer                                          :: Flog
        character(len=200)                               :: fic_vtk, fic_log
        character(len=100)                               :: simu_name
        !---param_algo_mod.f90
        type(PARAM_ALGO)                                 :: algo_param
        type(PARAM_ALGO_USER)                            :: user_param
        type(PARAM_ALGO_NLOC),allocatable,dimension(:)   :: Nloc_param
        !---material_mod.f90
        type(MATERIAL),allocatable,dimension(:)          :: MattotP
        type(NMATERIAL)                                  :: nmateriaux
        type(MATERIAL_REFERENCE)                         :: Matref
        type(INFO_VOXCOMP)                               :: VCinfo
        logical                                          :: testComposite
        type(NL_MODEL),allocatable,dimension(:)          :: Nloc_models
        type(TIMES_MAT)                                  :: times_m
        !---green_mod.f90
        type(TIMES_GREEN)                                :: times_g
        !---loading_mod.f90
        type(LOADING),allocatable,dimension(:)           :: load
        type(INITLOADEXT)                                :: initValExt
        type(LOCALLOAD)                                  :: local_load
        type(LOCALLOAD)                                  :: local_loadD
        type(PARAM_EXTRACT)                              :: extract
        integer                                          :: n_gradgradU = 0
        logical                                          :: add_def_star = .false.
        !---sortie_std_mod.f90
        type(BUFFER_STD)                                 :: buffer
        type(STRUCT_MAT_SORTIE),allocatable, dimension(:):: Matcomp_sortie
        type(TIMES_STD)                                  :: times_s
        !---io2_amitex_mod.f90
        type(TIMES_IO2)                                  :: times_io
        !---field_mod.f90
        integer,dimension(3)                             :: fft_start,fft_end,fft_size
        integer,dimension(3)                             :: amx_start,amx_end,amx_size
        integer(kind=8)                                  :: ntotP, ntotFP
#ifdef NEW_2DECOMP
        type(DECOMP_INFO),pointer                        :: ph_decomp   ! usefull for halo-cell
        type(DECOMP_INFO),pointer                        :: sp_decomp      
#else
        type(DECOMP_INFO)                                :: ph_decomp   ! usefull for halo-cell
        type(DECOMP_INFO)                                :: sp_decomp      
#endif  
        type(TIMES_FIELD)                                :: times_f
        real(mytype),allocatable,dimension(:,:,:,:)      :: FREQ
        complex(mytype),allocatable,dimension(:,:,:)     :: FREQ_2
        type(GRID_DIM)                                   :: grid
        !---NL_base_mod.f90
        type(CRIT_BASE)                                  :: crit_b
        type(TIMES_BASE)                                 :: times_b
        type(INFO_NL_BASE)                               :: info_NL_b
        !---error_mod.f90
        logical                                          :: error = .false.
        !---simu_mod
        integer                                          :: Igrid_simu=0     ! index of the grid in FFT_multigrid (from decomp_2d_fft)
                                                                             ! initialized at 0 to check initialization
        contains 
           GENERIC   :: WRITE(UNFORMATTED) => write_simu
           GENERIC   :: READ(UNFORMATTED)  => read_simu
           procedure :: write_simu
           procedure :: read_simu
  end type SIMU_AMITEX                                                       

  type(SIMU_AMITEX),allocatable,target,dimension(:)      :: SIMU
  type(SIMU_AMITEX)                                      :: SIMU00         ! used to refresh initialization variables
  integer                                                :: Iactive_simu   ! Index of the active simulation
  integer, pointer                                       :: Igrid_simu     ! Index on the grid for the active simulation

contains

subroutine write_simu(SIMU, iunit, istat, imsg)

  implicit none
  class(SIMU_AMITEX),intent(in)    :: SIMU  
  integer,intent(in)               :: iunit
  integer,intent(out)              :: istat
  character(len=*), intent(inout)  :: imsg
  
  !if (allocated(SIMU%SIG))       write(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%SIG
  if (allocated(SIMU%SIG0))      write(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%SIG0
  if (allocated(SIMU%DEF))       write(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%DEF
  if (allocated(SIMU%DEF0))      write(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%DEF0
  if (allocated(SIMU%DEF_nsym))  write(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%DEF_nsym
  if (allocated(SIMU%DEF_nsym0)) write(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%DEF_nsym0
  !if (allocated(SIMU%PK1))       write(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%PK1
  !if (allocated(SIMU%DEF_star))  write(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%DEF_star
  !if (allocated(SIMU%Dep))       write(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%Dep
  !if (allocated(SIMU%FluxD))     write(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%FluxD
  if (allocated(SIMU%FluxD0))    write(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%FluxD0
  if (allocated(SIMU%GradQD))    write(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%GradQD
  if (allocated(SIMU%GradQD0))   write(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%GradQD0
  if (allocated(SIMU%MattotP))   write(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%MattotP
  
  if (nrank == 0 ) write(OUTPUT_UNIT,*) "    WRITE_SIMU : iostat =", istat

end subroutine write_simu

subroutine read_simu(SIMU, iunit, istat, imsg)

  implicit none
  class(SIMU_AMITEX),intent(inout) :: SIMU  
  integer,intent(in)               :: iunit
  integer,intent(out)              :: istat
  character(len=*), intent(inout)  :: imsg
  
  !if (allocated(SIMU%SIG))       read(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%SIG
  if (allocated(SIMU%SIG0))      read(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%SIG0
  if (allocated(SIMU%DEF))       read(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%DEF
  if (allocated(SIMU%DEF0))      read(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%DEF0
  if (allocated(SIMU%DEF_nsym))  read(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%DEF_nsym
  if (allocated(SIMU%DEF_nsym0)) read(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%DEF_nsym0
  !if (allocated(SIMU%PK1))       read(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%PK1
  !if (allocated(SIMU%DEF_star))  read(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%DEF_star
  !if (allocated(SIMU%Dep))       read(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%Dep
  !if (allocated(SIMU%FluxD))     read(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%FluxD
  if (allocated(SIMU%FluxD0))    read(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%FluxD0
  if (allocated(SIMU%GradQD))    read(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%GradQD
  if (allocated(SIMU%GradQD0))   read(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%GradQD0
  if (allocated(SIMU%MattotP))   read(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) SIMU%MattotP

  if (nrank == 0 ) write(OUTPUT_UNIT,*) "    READ_SIMU : iostat", istat


end subroutine read_simu

!==================================================================================

!==================================================================================
!                         SUBROUTINE INIT_SIMU_COMMAND_LINE
!
!> Initialize a simulation object from the command line or from an input file
!!                                AND 
!! Associate the corresponding pointers
!!
!!  \param[out] SIMU         SIMU_AMITEX object, filled during the procedure
!!  \param[in]  version      string for writing the version number in log file
!!  \param[in]  file_cmd     (optional) name of the input file
!!  \param[out] help         True if command line with -help
!!  \param[out] coeff2print  index of the coefficient to output (vtk)
!!  \param[out] varint2print index of the internal variable to output (vtk)
!!  \param[in]  file_cmd     (optional) file name if initialisation from a file
!!  \param[in]  Igrid        (optional) Index of the grid on which the simulation is performed
!!                           default : Igrid=1
!!
!! Les pointers courants ci-dessous sont egalement positionnes
!!   fic_vtk,fic_log,Flog (amitex_mod)
!!   fft_start,fft_end,,ntotP,ntotFP (field_mod)
!!   mattotP (material_mod) pour desallouer MattotP0 rapidement
!!
!==================================================================================

subroutine init_simu_command_line(SIMU,version,help, coeff2print,varint2print,file_cmd, Igrid)
  
  implicit none
  
  type(SIMU_AMITEX), target, intent(out)   :: SIMU
  character(len=*), intent(in)             :: version
  character(len=*), optional, intent(in)   :: file_cmd
  integer, optional, intent(in)            :: Igrid
  
!!------------------------------------------------------------------------------
!>                                                     OPTIONS LIGNE DE COMMANDE
  integer,intent(out) :: coeff2print              !< indice du coeff a sortir (ligne de commande) 
  integer,intent(out) :: varint2print             !< indice de varint a sortir (ligne de commande) 
  logical,intent(out) :: help                     !< test pour avoir l'aide et sortir du code
  
!!------------------------------------------------------------------------------
!>                                                       DIMENSIONS DE LA GRILLE

  integer             :: nx, ny, nz                    !< dimensions de la cellule
  real(mytype)        :: dx, dy, dz, x0, y0 , z0       !< spacing and origin (parametre de l'entete vtk)
  !> verification de  la coherence des fichiers vtk (s'il y en a 2 en entree)
  integer             :: nx_z, ny_z, nz_z              !< dimensions de la cellule
  real(mytype)        :: dx_z,dy_z,dz_z,x0_z,y0_z,z0_z !< spacing and origin (parametre de l'entete vtk)

!!------------------------------------------------------------------------------
!>                                                      FICHIERS D'ENTREE/SORTIE

  character(len=200)  :: fic_numM, fic_numZ, fic_mat, fic_algo, fic_char, fic_local
                                                !< nom des fichiers
  integer             :: Fstd, FMstd, FZstd     !< FUnit pour initilisation des fichiers sortie
                                                
  real(mytype)        :: coeff_buffer_zstd      !< coefficient multiplicatif pour regler la taille
                                                !  du buffer_z pour sortie_std : 
                                                !  taille max du buffer = coeff_buffer_zstd x taille d'un champ de reels
                                                
!!------------------------------------------------------------------------------
!>                                                         MATERIAU DE REFERENCE

  real(mytype)        :: lambda, mu, lp2mu,dm   !< lambda,mu, lambda+2.mu,2.mu

!!------------------------------------------------------------------------------
!>                                                                      MATERIAU
!>                                                        SANS VOXELS COMPOSITES

  !> variables temporaires desalloues apres initialisation de Mattot 
  integer(kind=8),allocatable,dimension(:,:,:)  :: numZ         !< Champ des numeros de zone (num_zone.vtk)
  integer,allocatable,dimension(:,:,:)          :: numM         !< Champ des numeros de loi  (num_mat.vtk) 
                                                
!!------------------------------------------------------------------------------
!>                                                                    CHARGEMENT

  integer             :: nb_Tps, nb_param        !< nombre de pas de temps 
                                                  !< nombre de parametres externes

!!------------------------------------------------------------------------------
!>                                                                   MULTIGRILLE

  integer             :: Igrid0=1                   !< index of the grid on 

!!------------------------------------------------------------------------------
!>                                                                  PARALLELISME 
  
  integer             :: p_row=0, p_col, p_row_old  !< nb. de lignes, de colonnes 
                                                    !! pour la decomposition de 2decomp
  integer             :: ierror                   !< erreur relative au fonction MPI
  integer             :: type_mpi_m, type_mpi_z   !< type MPI associe a la lecture des fichiers vtk(si utilise)
                                                  !! MPI_INTEGER(1,2,4 ou 8)
  logical,save        :: test_decomp_init=.false. !< test si 2decomp a deja ete initialise
  logical,dimension(3):: pbc=.true.               !< Periodic BC si utlisation de Halo cells

!!------------------------------------------------------------------------------
!>                                                          SUIVI DE L'EXECUTION 

  double precision     :: t1      !< double precision au lieu de real_mytype requis par MPI_WTIME
                                  !! temps pour evaluer la duree de certaines fonctions

!!------------------------------------------------------------------------------
!>                                                                        DIVERS

  integer :: ndeb,nfin,i          !< indices de boucle et de tableau local_load%t1, t0, dt
  integer :: iload,iincr          !< restart starting point
  character(len=200)::  fic_vtk1,fic_log1 !< to adjust file output when using restart
  integer :: io_stat              !< erreur lors d'une lecture ou ecriture de fichier
  integer :: alloc_stat           !< erreur lors d'une allocation memoire


  allocate(numZ(1,1,1));deallocate(numZ) ! avoid gcc-warning
  allocate(numM(1,1,1));deallocate(numM) ! avoid gcc-warning
  
!!------------------------------------------------------------------------------
!>                                                               CHECK ARGUMENTS 
!!                                      
if (present(Igrid)) then
if (Igrid <= 0) then
      call amitex_abort('init_simu_command_line : Igrid must be > 0', 2,0)
end if
end if

!!------------------------------------------------------------------------------
!>                                               LECTURE DE LA LIGNE DE COMMANDE 
!!                                      ET DE LA VARIABLE D'ENVIRONNEMENT AMITEX
!!                              -> noms et unites logiques des fichiers
!!                              -> nx,ny,nz,dx,dy,dz, si pas de vtk en entree
!!                              -> coeff2print,varint2print,help

  error => SIMU%error

  nx = 0; ny = 0; nz = 0;
  dx = 0.; dy = 0.; dz = 0.;
  x0 = 0.; y0 = 0.; z0 = 0.;

  if (present(file_cmd)) then    
    call lire_commande(fic_mat, fic_algo, fic_char, fic_numM, fic_numZ, fic_vtk0,&
       fic_log0,&
       nx,ny,nz,dx,dy,dz,coeff2print,varint2print,help,file_cmd)
  else 
    call lire_commande(fic_mat, fic_algo, fic_char, fic_numM, fic_numZ, fic_vtk0,&
       fic_log0,&
       nx,ny,nz,dx,dy,dz,coeff2print,varint2print,help)
    if (help) return
  end if

  ! RESTART :  replace file names (fic_vtk1), during initialization
  call simu_restart_read_iload_iincr(fic_vtk0,iload,iincr)
  fic_vtk1 = fic_vtk0
  fic_log1 = fic_log0
  if (iload*iincr/=0) then ! Restart point
    write(fic_vtk1,fmt="(2A,I0,A,I0)") trim(fic_vtk0),"_",iload,"_",iincr
    write(fic_log1,fmt="(2A,I0,A,I0)") trim(fic_log0),"_",iload,"_",iincr
  end if
  
  ! AFFECTATION SIMU ET MAJ POINTEURS DE SIMU POUR SORTIES 
  SIMU%fic_vtk = fic_vtk1
  SIMU%fic_log = fic_log1
  SIMU%Flog = Flog0
  fic_vtk => SIMU%fic_vtk
  fic_log => SIMU%fic_log
  Flog => SIMU%Flog
  
  call check_amitex_abort(0)    !verifie les messages d'erreur d'amitex_abort dans lire_commande 

  call init_log()
  
  ! write version in log file
  if(nrank==0) write(Flog,"(A)") "VERSION "//trim(version)

  call get_AMITEXenv()

  call mpi_barrier(amitex_mpi_comm,ierror)

!!--------------------------------------      call amitex_abort('init_simu_command_line : trying to initialize Simu twice', 2,0)----------------------------------------
!>                                     COPIE DES FICHIERS XML DANS UN REPERTOIRE

  t1 = MPI_WTIME()
  call copyXML(fic_mat, fic_algo, fic_char, fic_vtk,r=0) ! r=0 : rank0 makes the copy
  call mpi_barrier(amitex_mpi_comm,ierror)
  if(nrank==0)write(Flog,"(A,E15.8)") "Duration COPY XML files (s) : ", MPI_WTIME() - t1

!!------------------------------------------------------------------------------
!>                                          LECTURE DES ENTETES DES FICHIERS VTK 
!!                                                -> nx,ny,nz,dx,dy,dz

  t1 = MPI_WTIME()
  !lecture de l'entete des fichiers numM et numZ
  if(trim(fic_numM)=="" .and. trim(fic_numZ)=="" .and. (nx==0 .or. ny==0 .or. nz==0)) then
     call amitex_abort("At least one vtk file must be given, or, provide &
                       &nx,ny and nz for a unique material with on zone per voxel (amitex)",1,0)
  else if(trim(fic_numM)=="" .and. (nx==0 .and. ny==0 .and. nz==0))then
     call read_header_vtk(fic_numZ, nx,ny,nz,dx,dy,dz,x0,y0,z0,type_mpi_z)
  else if(nx==0 .and. ny==0 .and. nz==0) then
     call read_header_vtk(fic_numM, nx,ny,nz,dx,dy,dz,x0,y0,z0,type_mpi_m)
     ! dans le cas ou les deux fichiers sont presents :
     ! verification de la concordance des dimensions
     if(trim(fic_numZ)/="")then
        call read_header_vtk(fic_numZ, nx_z,ny_z,nz_z,dx_z,dy_z,dz_z,x0_z,y0_z,z0_z,type_mpi_z)
        if(nx/=nx_z .or. ny/=ny_z .or. nz/=nz_z)then
           call amitex_abort("Grid dimensions (DIMENSIONS) in vtk files (mate and zone) are different (simu_mod)",1,0)
        end if
        if(dx/=dx_z .or. dy/=dy_z .or. dz/=dz_z)then
           call amitex_abort("Voxel sizes (SPACING) in vtk files (mate and zone) are different (simu_mod)",0,0)
        end if
        if(x0/=x0_z .or. y0/=y0_z .or. z0/=z0_z)then
           call amitex_abort("Origin coordinates (ORIGIN) in vtk files (mate and zone) are different (simu_mod)",0,0)
        end if
     end if
  end if
  ! affectation de la structure grid (type(GRID_DIM) definie dans field_mod)
  ! puis ecriture sur le fichier .log
  call initGrid(nx,ny,nz,dx,dy,dz,x0,y0,z0)
  SIMU%grid = grid
  call print_grid(Flog,0)

  call mpi_barrier(amitex_mpi_comm,ierror)
  if(nrank==0)write(Flog,"(A,E15.8)") "Duration READ HEADINGS vtk files (s) : ", MPI_WTIME() - t1
  call check_amitex_abort()

!!------------------------------------------------------------------------------
!>                                INITIALISATION 2DECOMP, 2DECOMP_FFT, FIELD_MOD 

  t1 = MPI_WTIME()
  p_row_old = p_row
  call set_pRow_pCol(p_row,p_col,nx,ny,nz)                !determination de p_row et p_col
if (.not. test_decomp_init) then                          !TODO : gestion de differentes decompositions 
#ifdef NEW_2DECOMP
  call decomp_2d_init(nx,ny,nz,p_row,p_col,pbc,amitex_mpi_comm) !initialisation decomp_2d
#else
  call decomp_2d_init(nx,ny,nz,p_row,p_col,pbc)           !initialisation decomp_2d
#endif
  test_decomp_init=.true.
else
  if (p_row_old /= p_row) &
   call amitex_abort("If 2DECOMP&FFT fails with 'Invalid 2D processor grid' message :"//ACHAR(10)// &
        "Try option '-2decomp low_p_row' or '-2decomp user' with appropriate decomposition for all simulations",0,0)
end if
  
  if (present(Igrid)) Igrid0 = Igrid

#ifdef NEW_2DECOMP
   if(decomp_2d_fft_get_ngrid() < Igrid0) then
      call decomp_2d_fft_set_ngrid(Igrid0)
   end if
   call decomp_2d_fft_init(PHYSICAL_IN_X, nx, ny, nz, Igrid0)
   call decomp_2d_fft_use_grid(Igrid0)
   call decomp_2d_fft_get_size(SIMU%fft_start,SIMU%fft_end,SIMU%fft_size)
   SIMU%ph_decomp => decomp_2d_fft_get_ph()
   SIMU%sp_decomp => decomp_2d_fft_get_sp()
   SIMU%amx_start = SIMU%ph_decomp%xst
   SIMU%amx_end   = SIMU%ph_decomp%xen
   SIMU%amx_size  = SIMU%ph_decomp%xsz 
#else
  !call decomp_2d_fft_init                                 !initialisation decomp_2d_fft
  call decomp_2d_fft_init(PHYSICAL_IN_X, nx, ny, nz, Igrid0)
  !call decomp_2d_fft_get_size(SIMU%fft_start,SIMU%fft_end,SIMU%fft_size) 
  call decomp_2d_fft_get_size(SIMU%fft_start,SIMU%fft_end,SIMU%fft_size, &
                              SIMU%amx_start,SIMU%amx_end,SIMU%amx_size) 
                                                          !fft_start,end,size 
                                                          !amx_start,end,size   => pointers in field_mod
#endif
  
  
  if (SIMU%Igrid_simu == 0) then  
      SIMU%Igrid_simu = Igrid0
  elseif (SIMU%Igrid_simu /= Igrid0) then
      call amitex_abort('init_simu_command_line : trying to initialize Simu twice', 2,0)
  end if
  
  Igrid_simu => SIMU%Igrid_simu
  ! below already done in decomp_2d_init
#ifdef NEW_2DECOMP
  amx_start     => SIMU%amx_start
  amx_end       => SIMU%amx_end
  amx_size      => SIMU%amx_size
#else  
  amx_start     => SIMU%amx_start
  amx_end       => SIMU%amx_end
  amx_size      => SIMU%amx_size
#endif

  ! get decomp from associated pointers (ph,sp) in decomp_2d_fft
#ifdef NEW_2DECOMP
  !SIMU%ph_decomp = decomp_2d_fft_get_ph()   ! DONE PREVIOUSLY
  !SIMU%sp_decomp = decomp_2d_fft_get_sp()
#else
  call get_decomp_fft_info(SIMU%ph_decomp,SIMU%sp_decomp)
#endif
  ph_decomp => SIMU%ph_decomp
  sp_decomp => SIMU%sp_decomp
  
  SIMU%ntotP = amx_size(1)*amx_size(2)*amx_size(3)
  SIMU%ntotFP = SIMU%fft_size(1) * SIMU%fft_size(2) * SIMU%fft_size(3)  

  fft_start =>   SIMU%fft_start
  fft_end   =>   SIMU%fft_end
  fft_size  =>   SIMU%fft_size
  ntotP     =>   SIMU%ntotP
  ntotFP    =>   SIMU%ntotFP

  call mpi_barrier(amitex_mpi_comm,ierror)

  if(nrank==0)write(Flog,"(3(A,I0))") "Dimensions : ",nx,"x",ny,"x",nz
  if(nrank==0)write(Flog,"(A,E15.8)") "Duration INITIALIZE 2DECOMP&FFT (s): ",MPI_WTIME() - t1

!!------------------------------------------------------------------------------
!>                                    INITIALISATION DES PARAMETRES D'ALGORITHME 

  call write_stdout0("before read_param,read_param_composite")

  t1 = MPI_WTIME()
  if(nrank==0)write(Flog,"(2A)")"Parametres de l'algorithme definis dans le fichier : ", trim(fic_algo)
  
  call read_param(fic_algo)

  !parametres specifiques "voxels composites"
  SIMU%testComposite = test_composite(fic_mat)
  if (SIMU%testComposite) call read_param_composite(fic_algo) 

  call MPI_Barrier(amitex_mpi_comm,ierror)  
  if(nrank==0) write(Flog,"(A,E15.8)") &
               "Duration READING algorithm parameters (s) : ", MPI_WTIME()-t1
  
  call check_amitex_abort()

  call write_stdout0("after read_param,read_param_composite")

  ! ecriture de la variable de type param_algo (pinceau 0) dans .log 
  ! repoussee car param_algo est aussi modifie dans read_mat (en presence de modeles non-locaux)


!!------------------------------------------------------------------------------
!>                                                 LECTURE DES DONNEES COMPOSITE 

  if (SIMU%testComposite) then
     t1 = MPI_WTIME()
     
     ! Lecture des données d'interphase
     call get_interphase_list(fic_mat)

     ! Lecture des données composites, allocation et initialisation de la structure MatComposite
     call read_composite(fic_mat)

    call  MPI_Barrier(amitex_mpi_comm,ierror)   
    if (nrank==0)    write(Flog,"(A,E15.8)") &
          "Duration READING composite voxels data files (s) : ", MPI_WTIME()-t1
  else 
     if (nrank==0) then
        write(Flog,"(A)") "Simulation without composite voxels"
     end if
  end if

!!------------------------------------------------------------------------------
!>                 INITIALISATION PARTIELLE DU TABLEAU DE "MATERIAL" : GEOMETRIE 
!>    (une partie de MattotP0 est affectee dans read_geom ou read_geom_composite)

  t1 = MPI_WTIME()
  ! Allocation des tableaux de numero de zone et de materiaux
  allocate(numM(amx_start(1):amx_end(1), amx_start(2):amx_end(2), amx_start(3):amx_end(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex1)",2)
  numM = 0

  allocate(numZ(amx_start(1):amx_end(1), amx_start(2):amx_end(2), amx_start(3):amx_end(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex1)",2)
  numZ = 0

  ! Creation des champs numM et numZ
  if(nrank==0) then
  write(Flog,"(A)")" "
  write(Flog,"(A)")"LECTURE DES PROPRIETES MATERIAU"
     if(trim(fic_numM)=="")then  
        write(Flog,"(A)")"Pas de numero de materiau defini dans un fichier vtk."
     else
        write(Flog,"(2A)")"Numeros de materiaux definis dans le fichier : "&
             , trim(fic_numM)
     end if
     if(trim(fic_numZ)=="")then  
        write(Flog,"(A)")"Pas de numero de zone defini dans un fichier vtk."
     else
        write(Flog,"(2A)")"Numeros de zones definis dans le fichier : "&
             , trim(fic_numZ)
     end if
  end if
  
  ! Affectation d'une partie de mattotP0
  !-------------------------------------
  if (fic_numM /="" .OR. fic_numZ /="") then
  if (SIMU%testComposite) then
     ! Cas avec voxels composites
     call read_geom_composite(fic_numM, fic_numZ,fic_vtk,type_mpi_m,type_mpi_z, numM, numZ,nx,ny)
  else
     ! Cas sans voxels composites
     call read_geom(fic_numM, fic_numZ, type_mpi_m,type_mpi_z, numM, numZ)
  end if
  end if
  ! Cas particulier d'un materiau continu (un voxel par zone)
  if (fic_numM =="" .and. fic_numZ =="") then
     numM=1
     call read_geom_continuous(nx,ny)
  end if
  
  call check_amitex_abort()

  call MPI_Barrier(amitex_mpi_comm,ierror)  
  if(nrank==0)write(Flog,fmt="(A,E15.8)") &
       "Duration READING the geometry (s): ",MPI_WTIME()-t1

!!------------------------------------------------------------------------------
!>         INITIALISATION PARTIELLE DU TABLEAU DE "MATERIAL" : COEFF, VAR.INT... 
!>       (une partie de MatotP est affectee dans read_mat ou read_mat_composite)

  call write_stdout0(&
       "before read_mat,read_mat_composite, print_param_algo, print_material_structure")

  t1 = MPI_WTIME()
  if(nrank==0)write(Flog,"(2A)")"Materials defined in file : ", trim(fic_mat)

  !affectation de mattotP0 et lecture materiau de reference
  if (algo_param0%Diffusion) allocate(Matref0%K0D(algo_param0%nVarD))
  if (SIMU%testComposite) then
     ! Cas avec voxels composites
     call read_mat_composite(fic_mat, Matref0%lambdamu0, nmateriaux0%n,nmateriaux0%n_composites)
  else
     ! Cas sans voxels composites
     if(.not. allocated(Matref0%K0D))  allocate(Matref0%K0D(1))  ! allocation bidon (TODO mettre ces variables dans un module)
     call read_mat(fic_mat, Matref0%lambdamu0, Matref0%K0D, nmateriaux0%n)
  end if

  ! Decompte des differents nombres de materiaux (dans la cellule, dans le pinceau,
  !                                               en comptant ou non les materiaux composites)
  call eval_nb_materiaux

  call MPI_Barrier(amitex_mpi_comm, ierror)
  if(nrank==0)write(Flog,fmt="(A,E15.8)") &
       "Duration BUILDING 'material' derived type : ",MPI_WTIME()-t1

  ! ecriture de la variable de type param_algo (pinceau 0) dans .log
  call print_param_algo(Flog,0)
  if (SIMU%testComposite) call print_param_algo_composite(Flog,0)

  ! ecriture de la structure materiau sur le pinceau 0 dans .log
  call print_material_structure(Flog,0,nx,ny,nz) 

  call write_stdout0(&
       "after read_mat,read_mat_composite, print_param_algo, print_material_structure")

!!------------------------------------------------------------------------------
!>                                          DESALLOCATION DES DONNES TEMPORAIRES

  ! Desallocation des  tableaux temporaires numM et numZ
  deallocate(numZ,stat=alloc_stat)
  deallocate(numM,stat=alloc_stat)

  ! Desallocation du tableau temporaire MatComposite
  if (allocated(MatComposite)) call deallocate_MatComposite()
  if (allocated(Interphases)) call deallocate_Interphases()

!!------------------------------------------------------------------------------
!>                                               INITIALISATION DE LA MATRICE C0 
!>                         

! MECANIQUE
if (algo_param0%Mechanics) then
  lambda = Matref0%lambdamu0(1)
  mu = Matref0%lambdamu0(2)
  dm = 2._mytype*mu
  lp2mu = lambda + dm
  ! Matrice C0 prenant en compte les notations pour Def et Sig
  allocate(Matref0%C0(algo_param0%nTensDef,algo_param0%nTensDef))
  if(algo_param0%nTensDef == 6) then
     Matref0%C0(1,:) = (/ lp2mu, lambda, lambda ,0._mytype,0._mytype,0._mytype /)
     Matref0%C0(2,:) = (/ lambda, lp2mu, lambda ,0._mytype,0._mytype,0._mytype /)
     Matref0%C0(3,:) = (/ lambda, lambda, lp2mu ,0._mytype,0._mytype,0._mytype /)
     Matref0%C0(4,:) = (/ 0._mytype,0._mytype,0._mytype,mu,0._mytype,0._mytype /)
     Matref0%C0(5,:) = (/ 0._mytype,0._mytype,0._mytype,0._mytype,mu,0._mytype /)
     Matref0%C0(6,:) = (/ 0._mytype,0._mytype,0._mytype,0._mytype,0._mytype,mu /)
  else
     Matref0%C0(1,:) = (/ lp2mu, lambda, lambda ,0._mytype,0._mytype,0._mytype,0._mytype,&
          0._mytype,0._mytype /)
     Matref0%C0(2,:) = (/ lambda, lp2mu, lambda ,0._mytype,0._mytype,0._mytype,0._mytype,&
          0._mytype,0._mytype /)
     Matref0%C0(3,:) = (/ lambda, lambda, lp2mu ,0._mytype,0._mytype,0._mytype,0._mytype,&
          0._mytype,0._mytype /)
     Matref0%C0(4,:) = (/ 0._mytype,0._mytype,0._mytype,dm,0._mytype,0._mytype,0._mytype,&
          0._mytype,0._mytype /)
     Matref0%C0(5,:) = (/ 0._mytype,0._mytype,0._mytype,0._mytype,dm,0._mytype,0._mytype,&
          0._mytype,0._mytype /)
     Matref0%C0(6,:) = (/ 0._mytype,0._mytype,0._mytype,0._mytype,0._mytype,dm,0._mytype,&
          0._mytype,0._mytype /)
     Matref0%C0(7,:) = (/ 0._mytype,0._mytype,0._mytype,0._mytype,0._mytype,0._mytype,dm,&
          0._mytype,0._mytype /)
     Matref0%C0(8,:) = (/ 0._mytype,0._mytype,0._mytype,0._mytype,0._mytype,0._mytype,&
          0._mytype,dm,0._mytype /)
     Matref0%C0(9,:) = (/ 0._mytype,0._mytype,0._mytype,0._mytype,0._mytype,0._mytype,&
          0._mytype,0._mytype,dm /)
  end if

! DIFFUSION - matrice du comportement "isotrope" n'est pas necessaire comme en meca 
end if

!!------------------------------------------------------------------------------
!>                                     INITIALISATION DE LA STRUCTURE CHARGEMENT 
!>                        (structures LOADING (load) et INITLOADEXT (initValExt)
!>                                                  initialisees dans read_load)

  t1 = MPI_WTIME()
  if(nrank==0)write(Flog,"(2A)")"Chargement defini dans le fichier : ", trim(fic_char)
  call read_load(fic_char,nmateriaux0%n,nb_tps)

  call MPI_Barrier(amitex_mpi_comm, ierror)
  if(nrank==0) then
     write(Flog,fmt="(A,E15.8)") &
          "Duree de la construction de la structure chargement : ",MPI_WTIME()-t1
     write(Flog,fmt="(A,I12)") &
          "Nombre de pas de chargement : ",nb_tps
  end if
  
  call check_amitex_abort()
  
  ! ecriture de la structure chargement sur le pinceau 0 dans .log
  call print_loading(Flog,0) 
  call print_initloadext(Flog,0) 
  call print_param_extract(Flog,0,nmateriaux0%n)
  call write_stdout0("End of loading initialization") 

!!------------------------------------------------------------------------------
!>                             INITIALISATION DES TABLEAUX DE CHARGEMENT COURANT 
!>                                local_load0%t0, local_load0%t1, local_load0%dt
!>                          &  local_loadD0%t0, local_loadD0%t1, local_loadD0%dt

  nb_param=initValExt0%nb_param

! MECANIQUE : local_load0
if (algo_param0%Mechanics) then

  !> Allocation et initialisation du chargement debut de pas : local_load0%t0
  allocate(local_load0%t0(2*algo_param0%nTensDef+1+nb_param+n_gradgradU0),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex3)",2)
  !! On initialise le type de pilotage et les valeurs associees a 0 et gradgradU a 0
  local_load0%t0 = 0
  !! On affecte la temperature initiale
  local_load0%t0(2*algo_param0%nTensDef+1) = initValExt0%temp
  !! Si besoin, on affecte les valeurs initiales des parametres externes
  if(nb_param>0) then
    ndeb = 2*algo_param0%nTensDef+1+1
    nfin = 2*algo_param0%nTensDef+1+initValExt0%nb_param
    local_load0%t0(ndeb:nfin) = initValExt0%param_values(:)
  end if

  !> Allocation et initialisation du chargement courant : local_load0%t1
  allocate(local_load0%t1(2*algo_param0%nTensDef+1+nb_param+n_gradgradU0),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex3)",2)
  local_load0%t1 = local_load0%t0

  !> Allocation et initialisation de d'increment de chargement : local_load0%dt
  allocate(local_load0%dt(2*algo_param0%nTensDef+1+nb_param+n_gradgradU0),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex3)",2)
  local_load0%dt = 0


! DIFFUSION : local_loadD0
elseif (algo_param0%Diffusion) then 
  !> Allocation et initialisation du chargement courant : local_loadD0%t0
  allocate(local_loadD0%t0(2*3*algo_param0%nVarD+1+nb_param),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex4)",2)
  !! On initialise le type de pilotage et les valeurs associees a 0
  local_loadD0%t0(1:2*3*algo_param0%nVarD) = 0
  !! On affecte la temperature initiale
  local_loadD0%t0(2*3*algo_param0%nVarD+1) = initValExt0%temp
  !! Si besoin, on affecte les valeurs initiales des parametres externes
  if(nb_param>0) local_loadD0%t0(2*3*algo_param0%nVarD+2:) = initValExt0%param_values(:)

  !> Allocation et initialisation du chargement courant : local_loadD0%t1
  allocate(local_loadD0%t1(2*3*algo_param0%nVarD+1+nb_param),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex4)",2)
  local_loadD0%t1 = local_loadD0%t0

  !> Allocation et initialisation du chargement courant : local_loadD0%dt
  allocate(local_loadD0%dt(2*3*algo_param0%nVarD+1+nb_param),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex4)",2)
  local_loadD0%dt = 0

end if

!!------------------------------------------------------------------------------
!>          INITIALISATION DES FICHIERS DE SORTIE STANDARD (.std, .mstd et .zstd)
!>          TODO : a deplacer dans sortie_std_mod.f90

  if(nrank==0)then
     !! Initialisation du fichier .std
     open(newunit=Fstd, file=trim(fic_vtk)//".std",form="formatted", status="replace", action="write",iostat= io_stat)  
     if ( io_stat /= 0 ) then
        write(OUTPUT_UNIT,"(3A,I0)")" Probleme a l'ouverture (fichier: ",trim(fic_vtk)//".std",") (amitex)",io_stat
        write(Flog,"(3A,I0)")" Probleme a l'ouverture (fichier: ",trim(fic_vtk)//".std",") (amitex)",io_stat
     else
        !! on appelle initSortie_std avec 0 pour preciser que le fichier concerne les moyennes globales
        call initSortie_std(Fstd,0)
        close(Fstd)
     end if
     call write_stdout0("End initialize .std") 

     !! Initialisation du fichier .mstd
     if (nmateriaux0%n_non_composites>1) then
        open(newunit=FMstd, file=trim(fic_vtk)//".mstd",form="formatted", status="replace", &
             action="write",iostat= io_stat)  
        if ( io_stat /= 0 ) then
           write(OUTPUT_UNIT,"(3A,I0)")" Probleme a l'ouverture (fichier: ",fic_vtk,".mstd) (amitex)",&
                io_stat
           write(Flog,"(3A,I0)")" Probleme a l'ouverture (fichier: ",fic_vtk,&
                ".mstd) (amitex)",io_stat
        else
           !! On donne l'oppose du nombre de materiaux a la fonction initSortie_std
           !! Cela permet de savoir que c'est une moyenne par zone et de donner le nombre de materiaux dans le domaine
           !call initSortie_std(FMstd,-nmateriaux0%n) ! Pour utiliser ancienne fonction sortie_std   
           call initSortie_std(FMstd,-nmateriaux0%n_non_composites)
           close(FMstd)
        end if
     end if
     call write_stdout0("End initialize .mstd") 

     !! Si on va sortir des moyennes par zone, on initialise les fichiers .zstd
     if(extract0%printZSTD) then 
        !do i=1,nmateriaux0%n   ! Pour utiliser ancienne fonction sortie_std     
        do i=1,nmateriaux0%n_non_composites
           if(extract0%printMeanZ(i)) then
              !! Initialisation du fichier .zstd
              write(fic_local,fmt="(A,I0,A)") trim(fic_vtk)//"_",i,".zstd"
              open(newunit=FZstd, file=trim(fic_local),form="formatted", status="replace", action="write",iostat= io_stat)  
              if ( io_stat /= 0 ) then
                 write(OUTPUT_UNIT,"(3A,I0)")" Probleme a l'ouverture (fichier: ",fic_local,") (amitex)",io_stat
                 write(Flog,"(3A,I0)")" Probleme a l'ouverture (fichier: ",fic_local,") (amitex)",io_stat
              else
                 !! On donne le numero de zone concerne a la fonction initSortie_std 
                 !! Ce nombre etant > 0 la fonction initialise un fichier de moyenne par zone 
                 call initSortie_std(FZstd,i)
                 close(FZstd)
              end if
           end if
        end do
     end if
     call write_stdout0("End initialize .zstd") 
  end if


  !sortie du programme en cas de probleme a l'ouverture d'un des fichiers
  call MPI_Bcast(io_stat,1,MPI_INTEGER,0,amitex_mpi_comm,ierror)
  if(io_stat/=0)then
     call write_stdout0("Problem when opening file") 
     call decomp_2d_fft_finalize
     call decomp_2d_finalize 
     call MPI_Abort(amitex_mpi_comm,io_stat,ierror)
     stop
  end if

  ! Initialisation et allocation des buffers pour sortie_std si on en a besoin. 
  if (nmateriaux0%n_non_composites > 1 .OR. extract0%printZSTD) then
     coeff_buffer_zstd = 1.
     call init_buffers(coeff_buffer_zstd)
  end if

  ! Initialisation de la structure permettant la prise en compte des voxels composites
  ! dans sortie_std_mod.
  ! modif LG 06-2019: uniquement en presence de vox Composites
  ! ATTENTION - la structure matcomp_sortie peut devenir très "lourde" en mémoir si 
  !             on utilise des voxels composites avec un grand nombre de zones
  if (nmateriaux0%n_composites > 0) then
        call write_stdout0("WARNING : composite voxels with a large number of zone is memory consuming")
        call write_file0(Flog,"WARNING : composite voxels with a large number of zone is memory consuming")
     call init_matcomp_sortie
  end if

  call write_stdout0("End of output initialization") 

!!------------------------------------------------------------------------------
!>                                                           MAJ POINTER MATTOTP 
!>                        + DESALLOCATIONS MATTOTP0 POUR LIMITE L'ESPACE MEMOIRE
!>                                            on duplique MattotP c'est pas top! 
!>                         -> a desallouer avant l'allocation des autres champs!

  allocate(SIMU%MattotP(size(MattotP0)))
  SIMU%MattotP = MattotP0
  MattotP => SIMU%MattotP
  call deallocate_MattotP0()


!!------------------------------------------------------------------------------
!>                         INITIALISATION DES CHAMPS DE FREQUENCE, FLUX ET GRADQ 
!>                                                             POUR LA DIFFUSION

if (algo_param0%Diffusion) then

  ! Flux, GradQ
  allocate(SIMU%FluxD0(amx_size(1)*amx_size(2)*amx_size(3),3,algo_param0%nVarD),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex5)",2)

  allocate(SIMU%GradQD0(amx_size(1)*amx_size(2)*amx_size(3),3,algo_param0%nVarD),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex6)",2)

  allocate(SIMU%FluxD(amx_size(1)*amx_size(2)*amx_size(3),3,algo_param0%nVarD),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex7)",2)

  allocate(SIMU%GradQD(amx_size(1)*amx_size(2)*amx_size(3),3,algo_param0%nVarD),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex8)",2)


  ! FluxF, GradQF espace spectral (pinceaux-Z, taille globale ~(nx/2+1)*ny*nz)
  allocate(SIMU%FluxDF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),&
       3,algo_param0%nVarD),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex9)",2)

  allocate(SIMU%gradQDF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),&
       3,algo_param0%nVarD),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex10)",2)

  allocate(SIMU%FREQ(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex11)",2)
  
  allocate(SIMU%FREQ_2(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex11b)",2)

  ! Stockage pour acceleration de convergence
  if(algo_param0%acc_CV) then
     allocate(SIMU%ACT3_RD(amx_size(1)*amx_size(2)*amx_size(3),3*algo_param0%nVarD,4),stat=alloc_stat)
     allocate(SIMU%ACT3_UD(amx_size(1)*amx_size(2)*amx_size(3),3*algo_param0%nVarD,4),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex12)",2)
  end if

  !initialisations
  SIMU%FluxD0  = 0._mytype
  SIMU%FluxD   = 0._mytype
  SIMU%GradQD0 = 0._mytype
  SIMU%GradQD  = 0._mytype
  SIMU%FluxDF  = 0._mytype
  SIMU%GradQDF = 0._mytype
  if(allocated(SIMU%ACT3_RD)) SIMU%ACT3_RD = 0._mytype
  if(allocated(SIMU%ACT3_UD)) SIMU%ACT3_UD = 0._mytype

  call initFREQ(SIMU%FREQ,'no_filter',0._mytype)
  call initFREQ_2(SIMU%Freq,SIMU%FREQ_2)
  call initFREQ(SIMU%FREQ,algo_param0%filter_typeD,algo_param0%filter_radiusD)
else
  !to avoid : forrtl: severe (408): fort: (7): Attempt to use pointer GRADQD when it is not associated with a target
  !           when gradQD called in ACT3
  allocate(SIMU%GradQD(0,3,algo_param0%nVarD),stat=alloc_stat)
end if

  call write_stdout0("End of allocations (Diffusion)") 

!!------------------------------------------------------------------------------
!>             INITIALISATION DES CHAMPS DE FREQUENCE, CONTRAINTE ET DEFORMATION 
!>                                       VARIABLES NON-LOCALES POUR LA MECANIQUE
if (algo_param0%Mechanics) then

  ! Contrainte, deformation
  allocate(SIMU%Sig(amx_size(1)*amx_size(2)*amx_size(3),1:6),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex13)",2)

  allocate(SIMU%Def(amx_size(1)*amx_size(2)*amx_size(3),1:algo_param0%nTensDef),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex14)",2)

  allocate(SIMU%Sig0(amx_size(1)*amx_size(2)*amx_size(3),1:6),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex15)",2)

  allocate(SIMU%Def0(amx_size(1)*amx_size(2)*amx_size(3),1:algo_param0%nTensDef),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex16)",2)
  if(.not. algo_param0%HPP) then
     allocate(SIMU%PK1(amx_size(1)*amx_size(2)*amx_size(3),1:algo_param0%nTensDef),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex17)",2)
  end if
!print *, "AV Def_nsym"  
  if (algo_param0%HPP_nsym) then
    ! allocation des champs dans le cas ou l'on souhaite prendre en compte le gradient complet
    ! du deplacement dans un calcul HPP
    allocate(SIMU%Def_nsym(amx_size(1)*amx_size(2)*amx_size(3),1:9),stat=alloc_stat)
    if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex18)",2)
    
    allocate(SIMU%Def_nsym0(amx_size(1)*amx_size(2)*amx_size(3),1:9),stat=alloc_stat)
    if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex19)",2)
  end if
     
  ! Contrainte, deformation, espace spectral (pinceaux-Z, taille globale ~(nx/2+1)*ny*nz)
  allocate(SIMU%SigF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),&
       1:algo_param0%nTensDef),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex20)",2)

  allocate(SIMU%DefF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),&
       1:algo_param0%nTensDef),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex21)",2)

!print *, "AV FREQ"  
  if (.NOT. algo_param0%Diffusion) then
     allocate(SIMU%FREQ(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex22)",2)
     allocate(SIMU%FREQ_2(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3)),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex22b)",2)
  end if
 
  if (algo_param0%HPP_nsym) then
    ! allocation des champs dans le cas ou l'on souhaite prendre en compte le gradient complet
    ! du deplacement dans un calcul HPP
    allocate(SIMU%DefF_nsym(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:9),stat=alloc_stat)
    if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex23)",2)
  end if

!print *, "AV ACV"    
  ! Stockage pour acceleration de convergence
  if (algo_param0%acc_CV) then
    if (algo_param0%HPP_nsym) then
    ! allocation des champs dans le cas ou l'on souhaite prendre en compte le gradient complet
    ! du deplacement dans un calcul HPP : acceleration de convergence sur le gradient de u. 
        allocate(SIMU%ACT3_R(amx_size(1)*amx_size(2)*amx_size(3),1:9,4),stat=alloc_stat)
        allocate(SIMU%ACT3_U(amx_size(1)*amx_size(2)*amx_size(3),1:9,4),stat=alloc_stat)
        if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex24)",2)        
    else
        allocate(SIMU%ACT3_R(amx_size(1)*amx_size(2)*amx_size(3),1:algo_param0%nTensDef,4),stat=alloc_stat)
        allocate(SIMU%ACT3_U(amx_size(1)*amx_size(2)*amx_size(3),1:algo_param0%nTensDef,4),stat=alloc_stat)
        if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex25)",2)
    end if
  end if

!print *, "AV CVFOR"    
  ! Stockage pour convergence forcee or substepping
  if (algo_param0%CVfor .OR. algo_param0%substepping) then
    if (algo_param0%HPP_nsym) then
       allocate(SIMU%CVFORsauv%Def(amx_size(1)*amx_size(2)*amx_size(3),1:9),stat=alloc_stat)
       if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex26)",2)
    else
       allocate(SIMU%CVFORsauv%Def(amx_size(1)*amx_size(2)*amx_size(3),1:algo_param0%nTensDef),stat=alloc_stat)
       if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex27)",2)
    end if
    SIMU%CVFORsauv%Def     = 0._mytype
    SIMU%CVFORsauv%critmin = 1.e100_mytype
  end if

!print *, "AV FVOLN"    
  ! Forces Volumiques Nodales
  ! TODO evaluer testFvolN en fonction des entrees
  test_FvolN = .false.
  if (test_FvolN) then
    allocate(SIMU%FvolN(amx_size(1)*amx_size(2)*amx_size(3),1:3),stat=alloc_stat)
    if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (FvolN)",2)
    allocate(SIMU%FvolNF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),&
       1:3),stat=alloc_stat)
    if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (FvolNF)",2)
  end if

  !initialisations
  SIMU%sig =0.
  SIMU%Def0=0.
  SIMU%Def =0.
  SIMU%Sig0=0.
  SIMU%DefF=0.
  SIMU%SigF=0.
  if(allocated(SIMU%PK1)) SIMU%PK1 = 0.
  if(allocated(SIMU%ACT3_R)) SIMU%ACT3_R=0.
  if(allocated(SIMU%ACT3_U)) SIMU%ACT3_U=0.
  if(allocated(SIMU%FvolN)) SIMU%FvolN=0.
  if(allocated(SIMU%FvolNF)) SIMU%FvolNF=0.
  if(allocated(SIMU%Def_nsym)) SIMU%Def_nsym=0.
  if(allocated(SIMU%Def_nsym0)) SIMU%Def_nsym0=0.
  if(allocated(SIMU%DefF_nsym)) SIMU%DefF_nsym=0.

  if (.NOT. algo_param0%Diffusion) then
     call initFREQ(SIMU%FREQ,'no_filter',0._mytype)
     call initFREQ_2(SIMU%FREQ,SIMU%FREQ_2)
     call initFREQ(SIMU%FREQ,algo_param0%filter_type,algo_param0%filter_radius)
  end if 

!print *, "AV init_Nloc"
  ! Allocation et initialisation des variables non locales
  if (algo_param0%Nloc) then
     allocate(SIMU%Nloc_models(size(Nloc_models0)))
     SIMU%Nloc_models = Nloc_models0
     Nloc_models => SIMU%Nloc_models
     call init_Nloc_variables(SIMU%Nloc,SIMU%GNloc)
     deallocate(Nloc_models0)    
  end if
  
else
  !to avoid : forrtl: severe (408): fort: (7): Attempt to use pointer DEF when it is not associated with a target
  !           when DEF called in ACT3
  allocate(SIMU%Def(0,1:algo_param0%nTensDef),stat=alloc_stat)
end if

  call write_stdout0("End of allocations (Mechanics)") 

!!------------------------------------------------------------------------------
!>                           ALLOCATION CRITERES POUR DIFFUSION (NL_base_mod) 
!>
if (algo_param0%Diffusion) then
   allocate(SIMU%crit_b%eqD(algo_param0%nVarD))
   allocate(SIMU%crit_b%FluxDMoy(algo_param0%nVarD))
   allocate(SIMU%crit_b%GradQDMoy(algo_param0%nVarD))
   allocate(SIMU%crit_b%CptbD(algo_param0%nVarD))
   SIMU%crit_b%eqD       = 1.e100_mytype
   SIMU%crit_b%FluxDMoy  = 1.e100_mytype
   SIMU%crit_b%GradQDMoy = 1.e100_mytype
   SIMU%crit_b%CptbD     = 1.e100_mytype
end if
SIMU%crit_b%eq     = 1.e100_mytype
SIMU%crit_b%Cptb   = 1.e100_mytype
SIMU%crit_b%SigMoy = 1.e100_mytype
SIMU%crit_b%DefMoy = 1.e100_mytype

!!------------------------------------------------------------------------------
!>                                        INITIALISATION DE TABLEAUX TEMPORAIRES 
!>

  allocate(SIMU%TEMPfield32(amx_size(1)*amx_size(2)*amx_size(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort(&
                     "espace memoire disponible insuffisant (amitex_fftp TEMPfield32)",2)
  SIMU%TEMPfield32=0.

  call write_stdout0("End of fields allocations and initializations") 

!!==============================================================================
!!                                                                 MAJ POINTEURS 
!!                     + DESALLOCATIONS VARIABLES D'INITIALISATION SI NECESSAIRE  
!!                         necessaire si objet ou composante d'objet allocatable
!!
!!  Remarque : certaines MAJ ont deja ete effectuee precedemment
!!             elles son commentees ci-dessous
!!==============================================================================
! ATTENTION : dans SIMU(1)%Matcomp_sortie = Matcomp_sortie0
!             le " = " realise l'allocation-affectation 
!             => ne passe pas avec intel15 (necessite allocation PUIS affectation)

  !---param_algo_mod
  SIMU%algo_param = algo_param0
  algo_param => SIMU%algo_param
  algo_param0 = SIMU00%algo_param  
  SIMU%user_param = user_param0
  user_param => SIMU%user_param
  user_param0 = SIMU00%user_param
  if(algo_param%Nloc) then 
     allocate(SIMU%Nloc_param(size(Nloc_param0)))
     SIMU%Nloc_param = Nloc_param0
     Nloc_param => SIMU%Nloc_param     
     deallocate(Nloc_param0)
  else
     nullify(Nloc_param)
  end if

  !---material_mod
!  allocate(SIMU%MattotP(size(MattotP0)))
!  SIMU%MattotP = MattotP0
!  MattotP => SIMU%MattotP
!  call deallocate_MattotP0()
  SIMU%nmateriaux=nmateriaux0
  nmateriaux => SIMU%nmateriaux
  nmateriaux0 = SIMU00%nmateriaux
  
  SIMU%Matref=Matref0
  Matref => SIMU%Matref
  call deallocate_Matref0()
  
  VCinfo => SIMU%VCinfo  ! les donnees de VCinfo sont initialisees a la construction (est-ce " norme fortran" compatible?)
                         ! ATTENTION : etaient egalement initialisees dans resolution_mod (supprimé)
  testComposite => SIMU%testComposite
  times_m => SIMU%times_m

!  if(algo_param%Nloc) then 
!     allocate(SIMU%Nloc_models(size(Nloc_models0)))      !deja fait
!     SIMU%Nloc_models = Nloc_models0
!     Nloc_models => SIMU%Nloc_models
!     deallocate(Nloc_models0)   
!  else
!     nullify(Nloc_models)
!  end if
  
  !---green_mod
  times_g => SIMU%times_g

  !---loading_mod
  allocate(SIMU%load(size(load0)))
  SIMU%load = load0
  load => SIMU%load
  !if(allocated(load0)) deallocate(load0) !ne marche pas ?!!!
  if(allocated(load0)) call deallocate_load(load0) 
  
  SIMU%initValExt = initValExt0
  initValExt => SIMU%initValExt
  call deallocate_initValExt(initValExt0)
  
  SIMU%local_load = local_load0
  local_load => SIMU%local_load
  call deallocate_local_load(local_load0)
  
  SIMU%local_loadD = local_loadD0
  local_loadD => SIMU%local_loadD
  call deallocate_local_load(local_loadD0)

  SIMU%extract = extract0
  extract => SIMU%extract
  call deallocate_param_extract0()
  
  SIMU%n_gradgradU = n_gradgradU0
  n_gradgradU => SIMU%n_gradgradU
  n_gradgradU0 = SIMU00%n_gradgradU
  
  SIMU%add_def_star = add_def_star0
  add_def_star => SIMU%add_def_star
  add_def_star0 = SIMU00%add_def_star
  
  !---sortie_std_mod
  if (allocated(buffer0%t)) then
     SIMU%buffer = buffer0
     buffer => SIMU%buffer
     buffer0 = SIMU00%buffer
  else 
     nullify(buffer)
  end if
  if (allocated(Matcomp_sortie0)) then
     allocate(SIMU%Matcomp_sortie(size(Matcomp_sortie0)))
     SIMU%Matcomp_sortie = Matcomp_sortie0 
     Matcomp_sortie => SIMU%Matcomp_sortie
     deallocate(Matcomp_sortie0)
  else
     nullify(Matcomp_sortie)
  end if
  call desallocation_sortie_std()
  times_s => SIMU%times_s

  !---io2_amitex_mod.f90
  times_io => SIMU%times_io
  
  !---amitex_mod - champs
  Sig       => SIMU%Sig
  Sig0      => SIMU%Sig0
  Def       => SIMU%Def
  Def0      => SIMU%Def0
  SigF      => SIMU%SigF
  DefF      => SIMU%DefF
  TEMPfield32 => SIMU%TEMPfield32
  CVFORsauv => SIMU%CVFORsauv 

  if (SIMU%add_def_star) then  
     Def_star  => SIMU%Def_star
  else
     nullify(Def_star)
  end if
  if (allocated(SIMU%PK1)) then  
     PK1       => SIMU%PK1
  else
     nullify(PK1)
  end if
  if (allocated(SIMU%Dep)) then  
     Dep       => SIMU%Dep
  else
     nullify(Dep)
  end if
  
  if (allocated(SIMU%FluxD)) then  
     FluxD     => SIMU%FluxD
     FluxD0    => SIMU%FluxD0
     GradQD    => SIMU%GradQD
     GradQD0   => SIMU%GradQD0
     FluxDF    => SIMU%FluxDF
     GradQDF   => SIMU%GradQDF
  else 
     nullify(FluxD)
     nullify(FluxD0)
     !nullify(GradQD)
     !to avoid : forrtl: severe (408): fort: (7): Attempt to use pointer GRADQD when it is not associated with a target
     !           when GRADQD called in ACT3
     GradQD    => SIMU%GradQD
     nullify(GradQD0)
     nullify(FluxDF)
     nullify(GradQDF)
  end if 
  if (allocated(SIMU%DefF_nsym)) then
     DefF_nsym => SIMU%DefF_nsym
     Def_nsym  => SIMU%Def_nsym
     Def_nsym0 => SIMU%Def_nsym0
  else
     nullify(DefF_nsym)
     nullify(Def_nsym)
     nullify(Def_nsym0)
  end if 
  if (allocated(SIMU%DepF)) then
     DepF => SIMU%DepF
  else
     nullify(DepF)
  end if  
  if (allocated(SIMU%FvolN)) then
     FvolN     => SIMU%FvolN
  else 
     nullify(FvolN)
  end if  
  if (allocated(SIMU%FvolNF)) then
     FvolNF    => SIMU%FvolNF
  else 
     nullify(FvolNF)
  end if  
  if (allocated(SIMU%ACT3_R)) then
     ACT3_R    => SIMU%ACT3_R
     ACT3_U    => SIMU%ACT3_U
  else
     nullify(ACT3_R)
     nullify(ACT3_U)
  end if
  if (allocated(SIMU%ACT3_RD)) then
     ACT3_RD    => SIMU%ACT3_RD
     ACT3_UD    => SIMU%ACT3_UD
  else
     nullify(ACT3_RD)
     nullify(ACT3_UD)
  end if
  if (allocated(SIMU%Nloc)) then
     Nloc  => SIMU%Nloc
     GNloc => SIMU%GNloc
  else
     nullify(Nloc)
     nullify(GNloc)
  end if
  
  !---amitex_mod - sorties fichier
!  fic_vtk => SIMU%fic_vtk
!  fic_log => SIMU%fic_log
!  Flog => SIMU%Flog
  fic_vtk1 = fic_vtk0   ! if RESTART :  replace filenames BACK 
  SIMU%fic_vtk = fic_vtk1
  fic_vtk => SIMU%fic_vtk

  SIMU%simu_name = "MASTER"
  simu_name => SIMU%simu_name
  fic_vtk0 = SIMU00%fic_vtk
  fic_log0 = SIMU00%fic_log
  Flog0    = SIMU00%Flog
  

  !---field_mod 
!  fft_start =>   SIMU%fft_start
!  fft_end   =>   SIMU%fft_end
!  fft_size  =>   SIMU%fft_size
!  amx_start    =>   SIMU%amx_start
!  amx_end      =>   SIMU%amx_end
!  amx_size     =>   SIMU%amx_size
!  ntotP     =>   SIMU%ntotP
!  ntotFP    =>   SIMU%ntotFP
   FREQ      => SIMU%FREQ 
   FREQ_2    => SIMU%FREQ_2 
   grid      => SIMU%grid
   times_f   =>   SIMU%times_f

  !---NL_base_mod 
  crit_b    => SIMU%crit_b
  times_b   => SIMU%times_b
  info_NL_b => SIMU%info_NL_b

  !---error_mod
  !error => SIMU%error   ! done at the begining of the subroutine
  

  !---simu_mod
  !Igrid_simu => SIMU%Igrid_simu ! done at the begining of the subroutine
  
!!==============================================================================
!!                                                  ADJUST STD FILES FOR RESTART 
  if (iload*iincr/=0) then
    call restart_adjust_std_files(fic_vtk1,iload,iincr)
  end if

  
end subroutine init_simu_command_line


!==================================================================================
!                         SUBROUTINE NULLIFY_POINTERS_SIMU
!
!> Nullify all the pointers required by a simulation
!!
!==================================================================================
subroutine nullify_pointers_simu()
  
  implicit none
  
  !---param_algo_mod
  nullify(algo_param)
  nullify(user_param)
  nullify(nloc_param)

  !---material_mod
  nullify(MattotP)
  nullify(nmateriaux)
  nullify(Matref)
  nullify(VCinfo)
  nullify(testComposite)
  nullify(times_m)
  nullify(Nloc_models)
  
  !---green_mod
  nullify(times_g)

  !---loading_mod
  nullify(load)
  nullify(initValExt)
  nullify(local_load)
  nullify(local_loadD)
  nullify(extract)
  nullify(n_gradgradU)
  nullify(add_def_star)
  
  !---sortie_std_mod
  nullify(buffer)
  nullify(Matcomp_sortie)
  nullify(times_s)

  !---io2_amitex_mod.f90
  nullify(times_io)

  !---amitex_mod - champs
  nullify(Sig)
  nullify(Sig0)
  nullify(Def)
  nullify(Def0)
  nullify(SigF)
  nullify(DefF)
  nullify(TEMPfield32)
  nullify(CVFORsauv)
  nullify(PK1)
  nullify(Def_star)
  nullify(Dep)
  nullify(FluxD)
  nullify(FluxD0)
  nullify(GradQD)
  nullify(GradQD0)
  nullify(FluxDF)
  nullify(GradQDF)
  nullify(DefF_nsym)
  nullify(DepF)
  nullify(Def_nsym)
  nullify(Def_nsym0)
  nullify(FvolN)
  nullify(FvolNF)
  nullify(ACT3_R)
  nullify(ACT3_U)
  nullify(ACT3_RD)
  nullify(ACT3_UD)
  nullify(Nloc)
  nullify(GNloc)
  
  !---amitex_mod - sorties fichier
  nullify(fic_vtk)
  nullify(fic_log)
  nullify(Flog)
  nullify(simu_name)

  !---field_mod 
  nullify(sp_decomp)
  nullify(ph_decomp)
  nullify(fft_start)
  nullify(fft_end)
  nullify(fft_size)
  nullify(amx_start)
  nullify(amx_end)
  nullify(amx_size)
  nullify(ntotP)
  nullify(ntotFP)
  nullify(times_f)
  nullify(FREQ)
  nullify(FREQ_2)
  nullify(grid)
    
  !---NL_base_mod - criteres
  nullify(crit_b)
  nullify(times_b)
  nullify(info_NL_b)
  
  !---error_mod
  nullify(error)
    
  !---simu_mod
  nullify(Igrid_simu)

end subroutine nullify_pointers_simu

!==================================================================================
!                         SUBROUTINE ASSOCIATE_POINTERS_SIMU
!
!> Associate all the pointers required by a simulation 
!>         (including the new pointers in decomp_2d_fft modified for multigrids)
!!
!!  \param[in] SIMU        array(:) of SIMU_AMITEX objects
!!  \param[in] I           index of the simulation in SIMU(:)
!!
!==================================================================================
subroutine associate_pointers_simu(SIMU,I)
  
  implicit none
  
  integer, intent(in)                                    :: I
  type(SIMU_AMITEX), dimension(:), target, intent(inout) :: SIMU


  !---param_algo_mod
  algo_param => SIMU(I)%algo_param
  user_param => SIMU(I)%user_param
  if(SIMU(I)%algo_param%Nloc) then
     Nloc_param => SIMU(I)%Nloc_param
  else
     nullify(Nloc_param)
  end if

  !---material_mod
  MattotP    => SIMU(I)%MattotP
  nmateriaux => SIMU(I)%nmateriaux
  Matref     => SIMU(I)%Matref
  VCinfo     => SIMU(I)%VCinfo  ! les donnees de VCinfo sont initialisees a la construction (est-ce " norme fortran" compatible?)
                                ! ATTENTION : etaient egalement initialisees dans resolution_mod (supprimé)
  testComposite => SIMU(I)%testComposite
  times_m    => SIMU(I)%times_m
  if(SIMU(I)%algo_param%Nloc) then
     Nloc_models => SIMU(I)%Nloc_models
  else
     nullify(Nloc_models)
  end if
  
  !---green_mod
  times_g    => SIMU(I)%times_g

  !---loading_mod
  load        => SIMU(I)%load
  initValExt  => SIMU(I)%initValExt
  local_load  => SIMU(I)%local_load
  local_loadD => SIMU(I)%local_loadD
  extract     => SIMU(I)%extract
  n_gradgradU => SIMU(I)%n_gradgradU
  add_def_star => SIMU(I)%add_def_star
  
  !---sortie_std_mod
  if (allocated(SIMU(I)%buffer%t)) then
     buffer   => SIMU(I)%buffer
  else 
     nullify(buffer)
  end if
  if (allocated(SIMU(I)%Matcomp_sortie)) then
     Matcomp_sortie => SIMU(I)%Matcomp_sortie
  else
     nullify(Matcomp_sortie)
  end if
  times_s    => SIMU(I)%times_s

  !---io2_amitex_mod
  times_io    => SIMU(I)%times_io

  !---amitex_mod - champs
  Sig         => SIMU(I)%Sig
  Sig0        => SIMU(I)%Sig0
  Def         => SIMU(I)%Def
  Def0        => SIMU(I)%Def0
  SigF        => SIMU(I)%SigF
  DefF        => SIMU(I)%DefF
  TEMPfield32 => SIMU(I)%TEMPfield32
  CVFORsauv   => SIMU(I)%CVFORsauv 

  if (allocated(SIMU(I)%PK1)) then  
     PK1       => SIMU(I)%PK1
  else
     nullify(PK1)
  end if
  if (add_def_star) then  
     Def_star  => SIMU(I)%Def_star
  else
     nullify(Def_star)
  end if
  if (allocated(SIMU(I)%Dep)) then  
     Dep       => SIMU(I)%Dep
  else
     nullify(Dep)
  end if

  if (allocated(SIMU(I)%FluxD)) then  
     FluxD     => SIMU(I)%FluxD
     FluxD0    => SIMU(I)%FluxD0
     GradQD    => SIMU(I)%GradQD
     GradQD0   => SIMU(I)%GradQD0
     FluxDF    => SIMU(I)%FluxDF
     GradQDF   => SIMU(I)%GradQDF
  else 
     nullify(FluxD)
     nullify(FluxD0)
     !nullify(GradQD)
     !to avoid : forrtl: severe (408): fort: (7): Attempt to use pointer GRADQD when it is not associated with a target
     !           when GRADQD called in ACT3
     GradQD    => SIMU(I)%GradQD
     nullify(GradQD0)
     nullify(FluxDF)
     nullify(GradQDF)
  end if 
  if (allocated(SIMU(I)%DefF_nsym)) then
     DefF_nsym => SIMU(I)%DefF_nsym
     Def_nsym  => SIMU(I)%Def_nsym
     Def_nsym0 => SIMU(I)%Def_nsym0
  else
     nullify(DefF_nsym)
     nullify(Def_nsym)
     nullify(Def_nsym0)
  end if 
  if (allocated(SIMU(I)%DepF)) then  
     DepF       => SIMU(I)%DepF
  else
     nullify(DepF)
  end if
  if (allocated(SIMU(I)%FvolN)) then
     FvolN     => SIMU(I)%FvolN
  else 
     nullify(FvolN)
  end if  
  if (allocated(SIMU(I)%FvolNF)) then
     FvolNF    => SIMU(I)%FvolNF
  else 
     nullify(FvolNF)
  end if  
  if (allocated(SIMU(I)%ACT3_R)) then
     ACT3_R    => SIMU(I)%ACT3_R
     ACT3_U    => SIMU(I)%ACT3_U
  else
     nullify(ACT3_R)
     nullify(ACT3_U)
  end if
  if (allocated(SIMU(I)%ACT3_RD)) then
     ACT3_RD    => SIMU(I)%ACT3_RD
     ACT3_UD    => SIMU(I)%ACT3_UD
  else
     nullify(ACT3_RD)
     nullify(ACT3_UD)
  end if
  if (allocated(SIMU(I)%Nloc)) then
     Nloc  => SIMU(I)%Nloc
     GNloc => SIMU(I)%GNloc
  else
     nullify(Nloc)
     nullify(GNloc)
  end if
  
  !---amitex_mod - sorties fichier
  fic_vtk => SIMU(I)%fic_vtk
  fic_log => SIMU(I)%fic_log
  Flog    => SIMU(I)%Flog
  simu_name => SIMU(I)%simu_name

  !---decomp_2d_fft  
#ifdef NEW_2DECOMP
  call decomp_2d_fft_use_grid(SIMU(I)%Igrid_simu)
  ! In case number of grids is changed, the array fft_engines may be reallocated, which
  ! may change the addresses of fft_engines(Igrid_simu)%ph and fft_engines(Igrid_simu)%sp,
  ! invalidating ph_decomp and sp_decomp.
  SIMU(I)%ph_decomp => decomp_2d_fft_get_ph()
  SIMU(I)%sp_decomp => decomp_2d_fft_get_sp()
#else
  call associate_pointers_decomp_2d_fft(SIMU(I)%Igrid_simu)
#endif

  !---field_mod 
  ph_decomp  =>   SIMU(I)%ph_decomp
  sp_decomp  =>   SIMU(I)%sp_decomp
  fft_start  =>   SIMU(I)%fft_start
  fft_end    =>   SIMU(I)%fft_end
  fft_size   =>   SIMU(I)%fft_size
  amx_start    =>   SIMU(I)%amx_start
  amx_end      =>   SIMU(I)%amx_end
  amx_size     =>   SIMU(I)%amx_size
  ntotP      =>   SIMU(I)%ntotP
  ntotFP     =>   SIMU(I)%ntotFP
  times_f    =>   SIMU(I)%times_f
  FREQ       =>   SIMU(I)%FREQ 
  FREQ_2     =>   SIMU(I)%FREQ_2 
  grid       =>   SIMU(I)%grid
  
  !---NL_base_mod - criteres
  crit_b    => SIMU(I)%crit_b
  times_b   => SIMU(I)%times_b
  info_NL_b => SIMU(I)%info_NL_b

  !---error_mod
  error     => SIMU(I)%error
  
  !---simu_mod
  Igrid_simu =>  SIMU(I)%Igrid_simu

  !---Active simu
  Iactive_simu = I
  
end subroutine associate_pointers_simu

!==================================================================================
!                         SUBROUTINE SIMU_RESTART_WRITE
!
!> WRITE RESTART FILES : infos.txt       for loading informations
!!                       simu_nrank.save for fields and material (ie. int. var.)
!!
!!
!! param[in]  iload    : saved partial loading for the restart
!! param[in]  iincr    : saved increment, in loading iload, for the restart 
!! param[in]  defMoy    : (optional) average stress at the begining of partial loading
!! param[in]  sigMoy    : (optional) average stress at the begining of partial loading
!! param[in]  gradQDMoy : (optional) average grad(T) at the begining of partial loading
!! param[in]  FluxDMoy  : (optional) average flux at the begining of partial loading
!!
!!
!! TODO : use io2_amitex_mod/print_coeff_varint_vtk instead of writing one file/proc
!!        
!==================================================================================
subroutine simu_restart_write(iload,iincr,defMoy00,sigMoy00, gradQDMoy00, FluxDMoy00)
  
  implicit none

  integer, intent(in) :: iload,iincr                                                   ! (inout) prevent nvfortran-warning
  real(mytype),dimension(algo_param%nTensDef),intent(in)   :: defMoy00, sigMoy00       ! 
  real(mytype),dimension(3,algo_param%nVarD),intent(in)    :: gradQDMoy00, FluxDMoy00  !
  
  character(12)       :: extfile
  integer             :: fid,ios,estat,cstat, ierror
  character(len=100)  :: mess,lineOS
  logical             :: exists 
  
  namelist /RESTART_INFOS1/iload,iincr!,prow,pcol
  namelist /RESTART_INFOS2/defMoy00, sigMoy00,gradQDMoy00,FluxDMoy00

  ios=iload;ios=iincr !avoid warning 'unused'
  
  ! CREATE RESTART DIRECTORY (if necessary -> PREFER CREATING the directory BEFORE executing AMITEX!!)
  exists = .false.  
  inquire(FILE=trim(fic_vtk)//"_restart", EXIST=exists)
#ifdef __INTEL_COMPILER
  inquire(DIRECTORY=trim(fic_vtk)//"_restart", EXIST=exists)
#endif

  if (.not. exists) then
    if (nrank == 0) then
      lineOS = "mkdir -p "//trim(fic_vtk)//"_restart"//achar(0) ! adding //achar(0) to circumvent an intel15 bug
      call execute_command_line (trim(lineOS),exitstat=estat,cmdstat=cstat,wait=.true.)
      if (cstat /= 0) then ! 2nd possibility if pb with 1st one 
        ios = mkdir(trim(fic_vtk)//"_restart"//achar(0), int(o'755',c_int16_t))
      end if
    end if

    !----WAIT for creation of directory (if not : problem on marquises when creating simu_xxx.save files on several nodes)
    !    1 - Creating the directory on several nodes is "long" (with both methods, execute_command_line or mkdir interface)!!!
    !    2 - => Create directory before launching amitex_fftp!!
    !    TODO : understand this...
    do while (.not. exists)
      call execute_command_line("sleep 1")
      inquire(FILE=trim(fic_vtk)//"_restart", EXIST=exists)
    end do
    !----END WAIT

    call MPI_BARRIER(amitex_mpi_comm,ierror)
  end if    

  ! SAVE SIMU
  mess = "write simu"
  write(extfile,'(I0)') nrank

  open(newunit=fid,file=trim(fic_vtk)//"_restart/simu_"//trim(adjustl(extfile))//".save",&
       form='UNFORMATTED', action='WRITE',status='REPLACE',position='REWIND',&
       access='SEQUENTIAL')
       
  write(UNIT=fid,IOSTAT=ios,IOMSG=mess) SIMU  !! TODO : i-loop over SIMU(i) ??
  close(fid)

  ! SAVE INFOS (iload,incr) : TODO : prow, pcol for safety check
  open(newunit=fid,file=trim(fic_vtk)//"_restart/infos.txt",&
       form='FORMATTED', action='WRITE',status='REPLACE',position='REWIND',&
       access='SEQUENTIAL')
  write(UNIT=fid,NML=RESTART_INFOS1)
  write(UNIT=fid,NML=RESTART_INFOS2)
  close(fid)
  
end subroutine simu_restart_write

!==================================================================================
!                         SUBROUTINE SIMU_RESTART_READ
!
!> READ RESTART FILES : only infos.txt if read_option == 0       (loading infos)
!!                      also simu_nrank/save if read_option == 1 (field infos)
!!
!! param[in]   read_option : if /=1, read infos.txt (iload,iincr, xxxMoy00)
!!                           if ==1, read also simu_nrank.save files to update SIMU 
!!                                        (SIMU : variable of the module)
!! param[out]  iload0    : last saved loading for the restart
!! param[out]  iincr0    : last saved increment, in loading iload0, for the restart 
!!                  if problem when reading infos.txt or if infos.txt do not exist
!!                  iload0 = incr0 = 0    
!!
!! param[out]  defMoy    : (optional) average stress at the begining of partial loading
!! param[out]  sigMoy    : (optional) average stress at the begining of partial loading
!! param[out]  gradQDMoy : (optional) average grad(T) at the begining of partial loading
!! param[out]  FluxDMoy  : (optional) average flux at the begining of partial loading
!!
!!
!! TODO : use io2_amitex_mod/print_coeff_varint_vtk instead of writing one file/proc
!!        -> write a new subroutine io2_amitex_mod/read_coeff_varint_vtk
!==================================================================================
subroutine simu_restart_read(iload0,iincr0, read_option,defMoy,sigMoy, gradQDMoy, FluxDMoy)
  
  implicit none

  integer, intent(out) :: iload0,iincr0
  integer,intent(in)   :: read_option
  real(mytype),dimension(algo_param%nTensDef),optional, intent(out)   :: defMoy, sigMoy     
  real(mytype),dimension(3,algo_param%nVarD),optional, intent(out)    :: gradQDMoy, FluxDMoy 

  real(mytype),dimension(algo_param%nTensDef)   :: defMoy00, sigMoy00     
  real(mytype),dimension(3,algo_param%nVarD)    :: gradQDMoy00, FluxDMoy00   
  character(12)       :: extfile
  integer             :: fid,ios1,ios2
  integer             :: iload,iincr
  character(len=100)  :: mess
  logical             :: test_file
  
  namelist /RESTART_INFOS1/iload,iincr
  namelist /RESTART_INFOS2/defMoy00, sigMoy00,gradQDMoy00,FluxDMoy00

  ! READ INFOS (iload,incr) : TODO : prow, pcol for safety check
  inquire(file=trim(fic_vtk)//"_restart/infos.txt", exist=test_file)
  
  if (test_file) then
    open(newunit=fid,file=trim(fic_vtk)//"_restart/infos.txt",&
       form='FORMATTED', action='READ',position='REWIND',&
       access='SEQUENTIAL')
    read(UNIT=fid,NML=RESTART_INFOS1,IOSTAT=ios1)
    read(UNIT=fid,NML=RESTART_INFOS2,IOSTAT=ios2)
    if (ios1 /=0 .or. ios2 /=0) then
      if (nrank==0) write(OUTPUT_UNIT,'(A)') "PROBLEM READING RESTART_INFOS file - set iload=iinc=0 for restart" 
      iload=0;iincr=0
    end if
    close(fid)
  else
    if (nrank==0) write(OUTPUT_UNIT,'(2A)') "NO RESTART FILE in ", trim(fic_vtk)//"_restart/infos.txt"      
    iload=0; iincr=0
  end if
  
  ! READ SIMU
  if (read_option==1) then
    mess = "write simu"
    write(extfile,'(I0)') nrank

    open(newunit=fid,file=trim(fic_vtk)//"_restart/simu_"//trim(adjustl(extfile))//".save",&
       form='UNFORMATTED', action='READ',position='REWIND',&
       access='SEQUENTIAL')
       
    read(UNIT=fid,IOSTAT=ios1,IOMSG=mess) SIMU !! TODO : i-loop over SIMU(i) ??
    close(fid)
  end if
  
  iload0 = iload
  iincr0 = iincr
  if (present(defMoy))    defMoy    = defMoy00
  if (present(sigMoy))    sigMoy    = sigMoy00
  if (present(FluxDMoy))  FluxDMoy  = FluxDMoy00
  if (present(gradQDmoy)) gradQDmoy = gradQDmoy00 
  
end subroutine simu_restart_read
!==================================================================================
!                         SUBROUTINE SIMU_RESTART_READ_ILOAD_IINCR
!
!> READ RESTART FILES (only ILOAD, IINCR in infos.txt)
!!
!==================================================================================
subroutine simu_restart_read_iload_iincr(fic_vtk0,iload0,iincr0)
  
  implicit none

  character(len=*),intent(in) :: fic_vtk0
  integer, intent(out)        :: iload0,iincr0

  integer             :: fid,ios
  integer             :: iload,iincr
  character(len=100)  :: mess
  logical             :: test_file
  
  namelist /RESTART_INFOS1/iload,iincr

  ! READ INFOS (iload,incr) : TODO : prow, pcol for safety check
  inquire(file=trim(fic_vtk0)//"_restart/infos.txt", exist=test_file)
  
  if (test_file) then
    open(newunit=fid,file=trim(fic_vtk0)//"_restart/infos.txt",&
       form='FORMATTED', action='READ',position='REWIND',&
       access='SEQUENTIAL')
    read(UNIT=fid,NML=RESTART_INFOS1,IOSTAT=ios)
    if (ios /=0) then
      if (nrank==0) write(OUTPUT_UNIT,'(A,I0)') "PROBLEM READING RESTART_INFOS file - set iload=iinc=0 for restart. ios =",ios
      iload=0;iincr=0
    end if
    close(fid)
  else
    if (nrank==0) write(OUTPUT_UNIT,'(2A)') "NO RESTART FILE in ", trim(fic_vtk0)//"_restart/infos.txt"      
    iload=0; iincr=0
  end if
  
  iload0 = iload
  iincr0 = iincr
  
end subroutine simu_restart_read_iload_iincr

!==================================================================================
!                         SUBROUTINE RESTART_ADJUST_STD_FILES
!
!> Remove last lines after the restart starting point
!!
!!   .std file   : OK
!!   .mstd files : OK
!!   .zstd files : OK
!==================================================================================
subroutine restart_adjust_std_files(fic_vtk0,iload0,iincr0)
  implicit none

  character(len=*),intent(in) :: fic_vtk0
  integer, intent(in)         :: iload0,iincr0
  
  logical                                       :: test_file
  integer                                       :: ios, nblines,nheader, i,j,imat,fid, nremove,nremove0
  integer                                       :: ierror
  integer(kind=8)                               :: nZone
  character(len=1000), dimension(:),allocatable :: lines
  character(len=1000)                           :: line, fic_name
  character(len=8)                              :: tmp_file="tmp_file"
  

if (nrank==0) then  
  nremove=0
  do i = 1, iload0-1
    nremove = nremove + load(i)%Nincr
  end do
  nremove = nremove + iincr0 + 1

  !--------std file
  inquire(file=trim(fic_vtk0)//".std", exist=test_file)
  if (test_file) then
    !open file
    fic_name = trim(fic_vtk0)//".std"
    open(newunit=fid,file=trim(fic_name),&
       form='FORMATTED', action='READ',position='REWIND',&
       access='SEQUENTIAL',iostat=ios)
    ! check ios   
    
    ! count lines
    nblines = 0
    nheader = 0
    do
      read(fid, '(A)', iostat=ios) line
      if (ios /= 0) exit
      nblines = nblines + 1
      if (line(1:1) == '#') nheader = nheader + 1
    end do  
    rewind(fid)
    
    ! store in vector lines
    allocate(lines(nblines))
    do i = 1, nblines
     read(fid, '(A)', iostat=ios) lines(i)
     if (ios /= 0) exit
    end do        
    close(fid)
    
    ! write in file
    open(newunit=fid, file=trim(fic_name), form='FORMATTED',status='replace', action='write',&
        access='SEQUENTIAL',iostat=ios)
    
    nremove0 = nblines - (nremove + nheader)
    do i = 1, nblines - nremove0
       write(fid, '(A)') trim(lines(i))
    end do
    close(fid)
    deallocate(lines)
    
  else
    write(OUTPUT_UNIT,'(A)'), trim(fic_vtk0)//".std doesn't exist (restart)"
  end if
  
  !--------mstd file
  nremove = nremove0 * nmateriaux%n
  fic_name=trim(fic_vtk0)//".mstd"
  inquire(file=trim(fic_name), exist=test_file)
  if (test_file) then
    !open file
    open(newunit=fid,file=trim(fic_vtk0)//".mstd",&
       form='FORMATTED', action='READ',position='REWIND',&
       access='SEQUENTIAL',iostat=ios)
    ! check ios   
    
    ! count lines
    nblines = 0
    do
      read(fid, '(A)', iostat=ios) line
      if (ios /= 0) exit
      nblines = nblines + 1
    end do  
    rewind(fid)
    
    ! store in vector lines
    allocate(lines(nblines))
    do i = 1, nblines
     read(fid, '(A)', iostat=ios) lines(i)
     if (ios /= 0) exit
    end do        
    close(fid)
    
    ! write in tmp_file
    open(newunit=fid, file=trim(fic_name), form='FORMATTED',status='replace', action='write',&
        access='SEQUENTIAL',iostat=ios)
    
    do i = 1, nblines - nremove
       write(fid, '(A)') trim(lines(i))
    end do
    close(fid)
    deallocate(lines)
    
  else
    write(OUTPUT_UNIT,'(A)') trim(fic_vtk0)//".mstd doesn't exist (restart)"
  end if
  
  !--------zstd file
  do imat = 1, nmateriaux%n
  
    nremove = nremove0 * nmateriaux%nz(imat)
    write(fic_name,'(A,I0,A)') trim(fic_vtk0)//"_",imat,".zstd" 
    inquire(file=trim(fic_name), exist=test_file)
    
    if (test_file) then
      !open file
      open(newunit=fid,file=trim(fic_name),&
         form='FORMATTED', action='READ',position='REWIND',&
         access='SEQUENTIAL',iostat=ios)
      ! check ios   
      
      ! count lines
      nblines = 0
      do
        read(fid, '(A)', iostat=ios) line
        if (ios /= 0) exit
        nblines = nblines + 1
      end do  
      rewind(fid)
      
      ! store in vector lines
      allocate(lines(nblines))
      do i = 1, nblines
       read(fid, '(A)', iostat=ios) lines(i)
       if (ios /= 0) exit
      end do        
      close(fid)
      
      ! write in tmp_file
      open(newunit=fid, file=trim(fic_name), form='FORMATTED',status='replace', action='write',&
        access='SEQUENTIAL',iostat=ios)
      
      do i = 1, nblines - nremove
         write(fid, '(A)') trim(lines(i))
      end do
      close(fid)
      deallocate(lines)
      
    else
      write(OUTPUT_UNIT,'(A)') trim(fic_name)//" doesn't exist (restart)"      
    end if
  end do
  
  !----------------------DELETE UNDESIRED FILES
  write(fic_name,'(A,I0,A,I0,A)') trim(fic_vtk0)//"_",iload0,"_",iincr0,".std" 
  call execute_command_line('rm "' // trim(fic_name) // '" > /dev/null 2>&1')
  write(fic_name,'(A,I0,A,I0,A)') trim(fic_vtk0)//"_",iload0,"_",iincr0,".mstd" 
  call execute_command_line('rm "' // trim(fic_name) // '" > /dev/null 2>&1')
  do imat = 1, nmateriaux%n
    write(fic_name,'(A,I0,A,I0,A,I0,A)') trim(fic_vtk0)//"_",iload0,"_",iincr0,"_",imat,".zstd" 
    call execute_command_line('rm "' // trim(fic_name) // '" > /dev/null 2>&1')  
  end do
end if  
  call MPI_Barrier(amitex_mpi_comm, ierror)

end subroutine restart_adjust_std_files


end module simu_mod
