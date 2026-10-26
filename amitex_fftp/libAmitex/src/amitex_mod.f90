!==============================================================================
!
!       MODULE AMITEX_MOD : 
!> Module comprenant les fonctions necessaires a amitex_fftp.f90
!!
!!  
!!  Subroutines
!! - lire_commande :    Lit la ligne de commande
!! - get_file :         Lit un nom de fichier en ligne de commande, retourne 
!!                      le nom du fichier et une unite de fichier s'il existe.
!! - set_pRow_pCol :    Determination du nombre de lignes et de colonnes de la 
!!                       decomposition a partir du nombre de processus utilises
!! - copyXML       :    Copie les fichiers XML
!! - get_AMITEXenv :    Renseigne la structure AMITEXenv (path, stat, length)
!!
!==============================================================================
module amitex_mod

  use ISO_FORTRAN_ENV

  use mpi
#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi,       only : nrank
#else
  use decomp_2d, only : nrank,mytype 
#endif  

  use io_amitex_mod
  use error_mod
  use amitex_mpi_mod, only: amitex_mpi_comm
  implicit none

  private
  
  !> Pointer "publiques" vers composantes de SIMU_AMITEX
  public :: Sig,Sig0,Def,Def0,Def_nsym,Def_nsym0,PK1,Def_star,Dep,SigF,DefF,DefF_nsym,DepF,&
            FluxD,FluxD0, GradQD,gradQD0,FluxDF, GradQDF,&
            ACT3_R,ACT3_U,ACT3_RD,ACT3_UD,&
            Nloc, GNloc, CVFORsauv,&
            FvolN,FvolNF, TEMPfield32,& 
            fic_vtk,fic_log,Flog,&
            simu_name
  
  !> Variables "publiques" utilisees lors de l'initialisation
  public :: fic_vtk0,fic_log0, Flog0

  !> Variables publiques a transformer en pointer (ou variables 'provisoires')
  public :: test_FvolN  ! TODO : variable provisoire...  

  !> Variables publiques
  public :: AMITEXenv, LOGO

  !> Types publiques (pour definition de SIMU_AMITEX)
  public :: SAUVCVFOR, NLOC_FIELD, GNLOC_FIELD

  !> Fonctions publiques
  public :: lire_commande, get_file,set_pRow_pCol, copyXML, get_AMITEXenv,&
            init_log

  !!------------------------------------------------------------------------------
  !>                                                      VARIABLE D'ENVIRONNEMENT 
  !>                                                  
  type envAMITEX
        character(len=:),allocatable         :: path
        integer                              :: path_stat
  end type envAMITEX

  !-----------------------------------------------
  type(envAMITEX)        :: AMITEXenv

  type sauvCVFOR
        real(mytype)                            :: critmin
        real(mytype),allocatable,dimension(:,:) :: Def
  end type sauvCVFOR
  
  !!------------------------------------------------------------------------------
  !>                                                          NOM DE LA SIMULATION
  
  character(len=100),pointer     :: simu_name 

  !!------------------------------------------------------------------------------
  !>                                                            FICHIERS DE SORTIE

  !> racine utile pour construire le nom des fichiers de sortie (valeur par defaut)
  character(len=200),target     ::  fic_vtk0="sortie" 
             
  !> Nom de fichiers de sortie .log (et .std) (valeur par defaut)
  character(len=200),target     ::  fic_log0="sortie.log"
  character(len=200),pointer    ::  fic_vtk, fic_log
 
  !> unite logique des fichiers de sortie  (valeurs par defaut)
  integer, target               ::  Flog0 != 4004 
  integer, pointer              ::  Flog
  
  !!------------------------------------------------------------------------------
  !>                                                                 DECOMPOSITION

  logical                                    :: decomp2D_low_p_row = .false.  ! decomposition with low p_row if .true.
  logical                                    :: decomp2D_user = .false.       ! user decomposition if .true.

  !!------------------------------------------------------------------------------
  !>                                                   CONTRAINTES ET DEFORMATIONS 
  !>                                                           tableaux paralleles

  !> dans l'espace reel : tableaux (ntot,Ntens) - pinceaux X
  !>                      notation 1D (1:amx_size(1)*amx_size(2)*amx_size(3),Ntens ou 4, nVarD)
  real(mytype),pointer,dimension(:,:)          :: Sig,Sig0 !< Contrainte de Cauchy (ntot,6)
  real(mytype),pointer,dimension(:,:)          :: Def,Def0 !< SI HPP Deformation linearisee (ntot,6)
                                                           !! SINON gradient du deplacement (ntot,9)
  real(mytype),pointer,dimension(:,:)          :: Def_nsym,Def_nsym0
                                                           !< Gradient du deplacement (ntot,9)
                                                           !! SI HPP et SI Prise en compte complete du gradient
                                                           !! du deplacement   
  real(mytype),pointer,dimension(:,:)          :: PK1      !< SI NON HPP: Contrainte de Piola Kirchoff 1
                                                           !! (ntot,9)
                                                           !! SI HPP : non alloué
  real(mytype),pointer,dimension(:,:)          :: Def_star !< Additionnal Deformation (or Displacement gradient)                                                    
  real(mytype),pointer,dimension(:,:)          :: Dep      !< Displacement field                                                   
  real(mytype),pointer,dimension(:,:,:)        :: FluxD,FluxD0     !< Flux (ntot,3,nVarD) pour la Diffusion
  real(mytype),pointer,dimension(:,:,:)        :: GradQD,GradQD0   !< Gradient Quantite Diffusion (ntot,3,nVarD)


  !> dans l'espace de Fourier : tableaux (nx/2+1,ny,nz,Ntens), ATTENTION AU SENS DE SigF - pinceaux Z
  !>                      notation 3D (fft_start(1):fft_end(1),...(2),...(3),Ntens ou 4, nVarD)
  complex(mytype),pointer,dimension(:,:,:,:)   :: SigF      !< Contrainte de Cauchy (HPP) ou PK1 (sinon) 
  complex(mytype),pointer,dimension(:,:,:,:)   :: DefF      !< Déformation linearisee (HPP) 
                                                            !! Gradient du deplacement (sinon))
  complex(mytype),pointer,dimension(:,:,:,:)   :: DefF_nsym !< Gradient du deplacement si prise
                                                            !! en compte complete du gradient du deplacement
  complex(mytype),pointer,dimension(:,:,:,:)   :: DepF      !< Deplacement
  complex(mytype),pointer,dimension(:,:,:,:,:) :: FluxDF    !< Flux (Diffusion) (nx/2+1,ny,nz,3,nVarD)
  complex(mytype),pointer,dimension(:,:,:,:,:) :: GradQDF   !< Gradient Quantite Diffusion (nx/2+1,ny,nz,3,nVarD)

  !!------------------------------------------------------------------------------
  !>                                                     FORCES VOLUMIQUES NODALES 
  !>                                                         tableaux paralleles
  real(mytype),pointer,dimension(:,:)          :: FvolN       !< Force Volumique aux Noeuds (ntot,3)
  complex(mytype),pointer,dimension(:,:,:,:)   :: FvolNF      !< Force Volumique aux Noeuds (nx,ny,nz,3)
  logical                                      :: test_FvolN

  !!------------------------------------------------------------------------------
  !>                                        CHAMP DE DEFORMATION IMPOSE HETEROGENE
  !>                                                         tableaux paralleles
  !real(mytype),allocatable,dimension(:,:)        :: Defimp_hetero       !< Champs de deformation impose (ntot,6)
  !complex(mytype),allocatable,dimension(:,:,:,:) :: Defimp_heteroF      !< Champs de deformation impose (nx,ny,nz,3)
  !logical                                        :: test_Defimp_hetero

  !!------------------------------------------------------------------------------
  !>                                         TABLEAUX CHAMPS DES MODELES NON LOCAUX 
  !>                                                           tableaux paralleles
  !>          un indice par modele non-local
  !!          champs dans lequel seront stockees les composantes des champs de varInt
  !!                              Nloc_field  : variables internes a deriver 
  !!                              GNloc_field : variables non locales issus des calculs non locaux
  !!                                       a repasser a MattotP(..)%VarInt(..,..)
  !---------------------------------------------------------------------------------
  type NLOC_FIELD
    !> dans l'espace reel : tableaux (ntot,param_algo%Nloc_models(i)%nNLoc) - pinceaux X
    !!                      notation 1D (1:amx_size(1)*amx_size(2)*amx_size(3),param_algo%Nloc_models(i)%nNLoc)
     real(mytype),allocatable,dimension(:,:)            :: Var
    !> dans l'espace de Fourier : tableaux (nx/2+1,ny,nz,param_algo%Nloc_models(i)%nNLoc)
    !!                      notation 3D (fft_start(1):fft_end(1),...(2),...(3),param_algo%Nloc_models(i)%nNLoc)
     complex(mytype),allocatable,dimension(:,:,:,:)     :: VarF
  end type NLOC_FIELD
  type GNLOC_FIELD
    !> dans l'espace reel : tableaux (ntot,param_algo%Nloc_models(i)%nGNLoc) - pinceaux X
    !>                      notation 1D (1:amx_size(1)*amx_size(2)*amx_size(3),param_algo%Nloc_models(i)%nGNLoc)
     real(mytype),allocatable,dimension(:,:)            :: Var
    !> dans l'espace de Fourier : tableaux (nx/2+1,ny,nz,param_algo%Nloc_models(i)%nGNLoc)
    !>                      notation 3D (fft_start(1):fft_end(1),...(2),...(3),param_algo%Nloc_models(i)%nGNLoc)
     complex(mytype),allocatable,dimension(:,:,:,:)     :: VarF
  end type GNLOC_FIELD


  type(NLOC_FIELD),pointer,dimension(:)             :: Nloc
  type(GNLOC_FIELD),pointer,dimension(:)            :: GNloc

  !!------------------------------------------------------------------------------
  !>                                       ACCELERATION DE CONVERGENCE (SI ACTIVE)
  !>                                         tableaux paralleles 
  !>                                            pour la mecanique (ntot,Ntens,4)
  !>                                            pour la diffusion (ntot,3*nVarD,4)

  real(mytype),pointer,dimension(:,:,:)      :: ACT3_R, ACT3_U
  real(mytype),pointer,dimension(:,:,:)      :: ACT3_RD, ACT3_UD

  !!------------------------------------------------------------------------------
  !>                                                CONVERGENCE FORCEE (SI ACTIVE)
  !>                                                 OR SUBSTEPPING (IF ACTIVATED)
  !>       variable regroupant le meilleur couple (critere , champ de deformation) 
  !>                                lors des convergences forcees / or substepping

  type(SAUVCVFOR),pointer                    :: CVFORsauv


  !!------------------------------------------------------------------------------
  !>                                                          TABLEAUX TEMPORAIRES

  real(kind=REAL32),pointer,dimension(:)     :: TEMPfield32   

  !!------------------------------------------------------------------------------
  !>                                                                          LOGO
  character(*), parameter :: LOGO=&
"                     _ _            ______ ______ _______ _____  "//ACHAR(10)//&
"     /\             (_) |          |  ____|  ____|__   __|  __ \ "//ACHAR(10)//&
"    /  \   _ __ ___  _| |_ _____  _| |__  | |__     | |  | |__) |"//ACHAR(10)//&
"   / /\ \ | '_ ` _ \| | __/ _ \ \/ /  __| |  __|    | |  |  ___/ "//ACHAR(10)//&
"  / ____ \| | | | | | | ||  __/>  <| |    | |       | |  | |     "//ACHAR(10)//&
" /_/    \_\_| |_| |_|_|\__\___/_/\_\_|    |_|       |_|  |_|     "//ACHAR(10)                                                          


contains 

!==============================================================================
!
!>   SUBROUTINE D'ECRITURE DE L'AIDE EN LIGNE DE COMMANDE
!!
!! Ecrit l'aide sur stdout \n
!!
!====================================================================
subroutine help_command_line()
call write_stdout0("")
call write_stdout0("")
call write_stdout0("RESUME THE POSSIBILITIES TO RUN AMITEX")
call write_stdout0("--------------------------------------")
call write_stdout0("for a detailled description : http://www.maisondelasimulation.fr/projects/amitex/html/user.html  ")
call write_stdout0("below: OPENMPI is used for a parallel use  ")
call write_stdout0("")
call write_stdout0("0/ help (this message)")
call write_stdout0(&
"mpirun -np n amitex_fftp -help   OR SIMPLY :  amitex_fftp -help")
call write_stdout0("")
call write_stdout0("")
call write_stdout0("1/ general case")
call write_stdout0(&
"mpirun -np n amitex_fftp  -nm <num_mat.vtk> -nz <num_zone.vtk> -a <algo.xml> -c <load.xml> -m <mat.xml> -s <output>")
call write_stdout0("")
call write_stdout0("")
call write_stdout0("2/ one material (assumes num_mat.vtk full of 1)")
call write_stdout0(&
"mpirun -np n amitex_fftp  -nz <num_zone.vtk> -a <algo.xml> -c <load.xml> -m <mat.xml>  -s <output>")
call write_stdout0("")
call write_stdout0("")
call write_stdout0("3/ one zone (assumes num_zone.vtk full of 1)")
call write_stdout0(&
"mpirun -np n amitex_fftp  -nm <num_mat.vtk>  -a <algo.xml> -c <load.xml> -m <mat.xml>  -s <output>")
call write_stdout0("")
call write_stdout0("")
call write_stdout0("4/ one material with one zone per voxel (assumes num_zone.vtk varying from 1 to the number of voxels) ")
call write_stdout0("   assumes also dx=dy=dz=1. ")
call write_stdout0(&
"mpirun -np n amitex_fftp  -a <algo.xml> -c <load.xml> -m <mat.xml> -NX nx -NY ny -NZ nz -s <output>")
call write_stdout0("")
call write_stdout0("")
call write_stdout0("5/ one material with one zone per voxel (assumes num_zone.vtk varying from 1 to the number of voxels) ")
call write_stdout0(&
"mpirun -np n amitex_fftp  -a <algo.xml> -c <load.xml> -m <mat.xml> -NX nx -NY ny -NZ nz -DX dx -DY dy -DZ dz -s <output>")
call write_stdout0("")
call write_stdout0("")
call write_stdout0("ADDITIONNAL OPTIONS (TO CHECK INPUT DATA)")
call write_stdout0("-----------------------------------------")
call write_stdout0("  -coeff2print  i : generate a vtk file for the ith coefficient")
call write_stdout0("  -varint2print i : generate a vtk file for the ith internal variable")
call write_stdout0("  (1 field per material complemented with 0, no field if material has  less than i coeff. or int. var.)")
call write_stdout0("")
call write_stdout0("ADDITIONNAL OPTIONS (TO CHOOSE 2D PENCIL DECOMPOSITION)")
call write_stdout0("-------------------------------------------------------")
call write_stdout0("  The default tries to find a decomposition as square as possible (p_row as close as possible to p_col")
call write_stdout0("  -2decomp low_p_row : will chose a 2D pencil decomposition with the minimum p_row value")
call write_stdout0("  -2decomp user : a user-defined decomposition  can be used from file user_2decomp.txt ")
call write_stdout0("           containing two lines (p_row and p_col)")
call write_stdout0("")
call write_stdout0("")

end subroutine help_command_line


!==============================================================================
!
!>   SUBROUTINE DE LECTURE DES ARGUMENTS DE LA LIGNE DE COMMANDE.
!!
!! Lit la ligne de commande qui doit etre de la forme\n
!!
!!mpirun -np n ./amitex_fftp  -nm <num_mat.vtk> -nz <num_zone.vtk> -a <algo.xml> -c <charg.xml> \n
!!                            -m <mat.xml> -s <sortie> \n
!!OU
!!
!!mpirun -np n ./amitex_fftp  -nz <num_zone.vtk> -a <algo.xml> -c <charg.xml> \n
!!                            -m <mat.xml>  -s <sortie> \n
!!
!!OU
!!
!!mpirun -np n ./amitex_fftp  -nm <num_mat.vtk>  -a <algo.xml> -c <charg.xml> \n
!!                            -m <mat.xml>  -s <sortie> \n
!!
!!OU
!!
!!mpirun -np n ./amitex_fftp  -a <algo.xml> -c <charg.xml> \n
!!                            -m <mat.xml> -NX nx -NY ny -NZ nz -s <sortie> \n
!!
!!OU
!!
!!mpirun -np n ./amitex_fftp  -a <algo.xml> -c <charg.xml> \n
!!                            -m <mat.xml> -NX nx -NY ny -NZ nz -DX dx -DY dy -DZ dz -s <sortie> \n
!!
!!
!!OU, POUR VERIFIER LA MISE EN DONNEE (le calcul ne se lance pas)
!!
!! ajouter les options  -coeff2print i
!!                      -varint2print i
!!
!! Ces deux options fournissent une sortie vtk du coeff (varint) i de chaque materiau 
!!
!!
!! - n : nombre de processus utilises
!!
!! - -nm <num_mat.vtk>  nom du fichier vtk contentant les numeros de materiaux de chaque voxel
!! - -nz <num_zone.vtk>  nom du fichier vtk contentant les numeros de zones de chaque voxel
!! - -a <algo.xml>  nom du fichier xml contenant les parametres de l'algorithme
!! - -c <charg.xml> nom du fichier xml decrivant le chargement
!! - -m <mat.xml>   nom du fichier xml contenant les informations des materiaux
!! - -NX nx         dimension suivant x (si -nm et -nz non fournis)
!! - -NY ny         dimension suivant y (idem)
!! - -NZ nz         dimension suivant z (idem)
!! - -DX dx         taille de voxel suivant x (si -nm et -nz non fournis, optionnel, 1. par defaut)
!! - -DY dy         taille de voxel suivant y (idem)
!! - -DZ dz         taille de voxel suivant z (idem)
!!
!!optionnel:
!! - -s  <sortie>   racine du nom des fichiers de sortie du programme
!!                Le fichier d'information sera nomme "nomFichier.std"
!!                Les fichiers vtk seront nommes "nomFichier_i.vtk" ou i est
!!                l'indice du pas de temps auquel le fichier est ecrit.
!! 
!! \param[out] fic_mat: (chaine de caracteres) nom du fichier de materiaux
!! \param[out] fic_algo: (chaine de caracteres) nom du fichier des parametres de l'algorithme
!! \param[out] fic_char: (chaine de caracteres) nom du fichier de chargement
!! \param[out] fic_log: (chaine de caracteres) nom du fichier de sortie.log
!! \param[out] fic_vtk: (chaine de caracteres) radical du nom des fichiers vtk
!! \param[out] nx,ny,nz: dimension de la cellule si les fichiers vtk ne sont pas donnes
!! \param[out] dx,dy,dz: tailles des voxels si les fichiers vtk ne sont pas donnes
!! \param[out] coeff2print numero de coeff a sortir (vtk) pour verification
!! \param[out] varint2print numero de varint a sortir (vtk) pour verification
!! \param[out] help     logical, true if "amitex_fftp -help" 
!! \param[in]  file_cmd (optional) file name for input (instead of command line arg)
!!
!====================================================================
subroutine lire_commande(fic_mat, fic_algo, fic_char, fic_numM, fic_numZ, fic_vtk,&
     fic_log,&
     nx,ny,nz,dx,dy,dz,coeff2print,varint2print,help,file_cmd)

  implicit none

  character(len=200),intent(out)  :: fic_mat, fic_algo, fic_char, fic_vtk, fic_log,fic_numM, fic_numZ
           !> erreur compilo intel si len=*, incompatible avec NAMELIST
  character(len=200)            :: fic_std
  integer                       :: Fmat, Falgo, Fchar,Fnmat,Fnzone
  character(len=200)            :: arg !arguments du code fortran lance en ligne de commande
  character(len=200)            :: err_msg
  integer                       :: i, FT,ios
  integer,intent(out)           :: nx,ny,nz
  real(mytype),intent(out)      :: dx,dy,dz
  integer,intent(out)           :: coeff2print,varint2print
  logical                       :: existe, file_exists
  integer,dimension(20)         :: test_arg0, test_arg1
  logical                       :: testPB
  logical,intent(out)           :: help
  character(len=*),optional,intent(in)   :: file_cmd 
  integer                                :: Fcmd
  
  namelist /CMD/ fic_mat, fic_algo, fic_char, fic_vtk, fic_numM, fic_numZ,&
                 nx,ny,nz,dx,dy,dz, coeff2print, varint2print

  !valeurs par defaut
  help = .false.
  coeff2print  = 0
  varint2print = 0

  !test arguments
  test_arg0 = 0
  test_arg1 = 0

  !valeur initiales de nx,ny,nz
  nx=0
  ny=0
  nz=0

  !valeur par defaut de dx,dy,dz
  dx=1._mytype
  dy=1._mytype
  dz=1._mytype

  ! valeurs par defaut des fichier de sortie (passees en variables public du module)
  !fic_vtk="sortie"
  fic_std ="sortie.std"
  !fic_log="sortie.log"

  !! fichier d'entree: nom vides, unites logiques a -1
  ! en fin de subroutine:
  ! unite logique <0 => erreur dans la ligne de commande
  fic_mat=""
  fic_algo=""
  fic_char=""
  fic_numM=""
  fic_numZ=""
  Fmat = -1
  Falgo = -1
  Fchar = -1
  Fnzone = -1
  Fnmat = -1
  
if (.not. present(file_cmd)) then !******* CAS LECTURE LIGNE DE COMMANDE

  ! start file unit after default values 0 (stderr), 5 (stdin) and 6 (stdout)
  FT = Flog0 +11

  do i=1,command_argument_count(),2

     call getarg(i,arg)

     if(arg =="-m")then
        call get_file(i+1, fic_mat,Fmat , FT, existe)
        if(.NOT.existe) call write_stdout0("File "//trim(fic_mat)//" absent")
        test_arg1(1)=1
     else if(arg =="-a")then
        call get_file(i+1, fic_algo,Falgo, FT, existe)
        if(.NOT.existe) call write_stdout0("File "//trim(fic_algo)//" absent")
        test_arg1(2)=1
     else if(arg =="-c")then
        call get_file(i+1, fic_char,Fchar, FT, existe)
        if(.NOT.existe) call write_stdout0("File "//trim(fic_char)//" absent")
        test_arg1(3)=1
     else if(arg =="-nm")then
        call get_file(i+1, fic_numM,Fnmat, FT, existe)
        if(.NOT.existe) call write_stdout0("File "//trim(fic_numM)//" absent")
        test_arg1(4)=1
     else if(arg =="-nz")then
        call get_file(i+1, fic_numZ,Fnzone, FT, existe)
        if(.NOT.existe) call write_stdout0("File "//trim(fic_numZ)//" absent")
        test_arg1(5)=1
     else if(arg=="-NX")then
        call getarg(i+1,arg)
        read(arg,"(I8)") nx
        test_arg1(6)=1
     else if(arg=="-NY")then
        call getarg(i+1,arg)
        read(arg,"(I8)") ny
        test_arg1(7)=1
     else if(arg=="-NZ")then
        call getarg(i+1,arg)
        read(arg,"(I8)") nz
        test_arg1(8)=1
     else if(arg=="-DX")then
        call getarg(i+1,arg)
        read(arg,*) dx
        test_arg1(9)=1
     else if(arg=="-DY")then
        call getarg(i+1,arg)
        read(arg,*) dy
        test_arg1(10)=1
     else if(arg=="-DZ")then
        call getarg(i+1,arg)
        read(arg,*) dz
        test_arg1(11)=1
     else if(arg=="-s")then
        call getarg(i+1,arg)
        read(arg,"(A)") fic_vtk
        fic_std = trim(fic_vtk)//".std"
        fic_log = trim(fic_vtk)//".log"
        test_arg1(12)=1
     else if(arg=="-coeff2print")then
        call getarg(i+1,arg)
        read(arg,"(I8)") coeff2print
        test_arg1(13)=1
     else if(arg=="-varint2print")then
        call getarg(i+1,arg)
        read(arg,"(I8)") varint2print
        test_arg1(14)=1
     else if(arg=="-2decomp")then ! optional argument
        call getarg(i+1,arg)
        if (trim(arg)=="low_p_row") decomp2D_low_p_row = .true.
        if (trim(arg)=="user") decomp2D_user = .true.
     else if(arg=="-help")then
        call help_command_line() 
        help=.true. 
        return
     else
        call write_stdout0("invalid command line")
        call write_stdout0("unknown argument : "//trim(arg))
        call amitex_abort("unknown argument '"//trim(arg)//&
                         &"' :'amitex_fftp -help' for a description of valid command lines ",2,0)
     end if
  end do

  ! TEST les arguments de la ligne de commande
  !-------------------------------------------
  testPB=.true.

  !!amitex_fftp  -nm <num_mat.vtk> -nz <num_zone.vtk> -a <algo.xml> -c <charg.xml> 
  !!                            -m <mat.xml>  
  test_arg0=0
  test_arg0( (/1,2,3,4,5/) )=1
  if (all(test_arg0(1:11) == test_arg1(1:11))) testPB=.false.

  !!amitex_fftp  -nz <num_zone.vtk> -a <algo.xml> -c <charg.xml> 
  !!                            -m <mat.xml>  
  test_arg0=0
  test_arg0( (/1,2,3,5/) )=1
  if (all(test_arg0(1:11) == test_arg1(1:11))) testPB=.false.

  !!amitex_fftp  -nm <num_mat.vtk>  -a <algo.xml> -c <charg.xml> 
  !!                            -m <mat.xml>  
  test_arg0=0
  test_arg0( (/1,2,3,4/) )=1
  if (all(test_arg0(1:11) == test_arg1(1:11))) testPB=.false.

  !!amitex_fftp  -a <algo.xml> -c <charg.xml> 
  !!                            -m <mat.xml> -NX nx -NY ny -NZ nz 
  test_arg0=0
  test_arg0( (/1,2,3,6,7,8/) )=1
  if (all(test_arg0(1:11) == test_arg1(1:11))) testPB=.false.

  !!amitex_fftp  -a <algo.xml> -c <charg.xml> 
  !!                            -m <mat.xml> -NX nx -NY ny -NZ nz -DX dx -DY dy -DZ dz 
  test_arg0=0
  test_arg0( (/1,2,3,6,7,8,9,10,11/) )=1
  if (all(test_arg0(1:11) == test_arg1(1:11))) testPB=.false.


  ! TRAITEMENT DES PBS 
  !-------------------------------------------

  if (testPB) then 
        call amitex_abort(&
        "invalid command line : 'amitex_fftp -help' for a description of valid command lines ",2,0)
  end if

  call check_amitex_abort(0)

else !***********************************************CAS LECTURE FICHIER

  !                 Ouverture du fichier d'entree
  !----------------------------------------------
  INQUIRE(file=trim(file_cmd),exist=file_exists)
  if(.not. file_exists) then
    write(err_msg,fmt="(3A)") "The file : ",&
         trim(file_cmd)," does not exist (lire_commande)"
    call amitex_abort(trim(err_msg),2,0)
  end if
  open(newunit=Fcmd, file=trim(file_cmd),form="formatted", status="old", action="read",iostat= ios)
  if ( ios /= 0 ) then
     write(err_msg,fmt="(3A,I0,A)") "Probleme a l'ouverture du fichier : ",trim(file_cmd),&
          " , ",ios," (lire_commande)"
     call amitex_abort(trim(err_msg),2,0) 
  end if

  !                              Lecture NAMELIST
  !----------------------------------------------
  read(Fcmd, nml=CMD)
  close(Fcmd)
  fic_std = trim(fic_vtk)//".std"
  fic_log = trim(fic_vtk)//".log"
  
  !                        TEST EXISTENCE FICHIERS
  !-----------------------------------------------
  INQUIRE(FILE=fic_mat, EXIST=file_exists)
  if (file_exists) Fmat = 1
  INQUIRE(FILE=fic_algo, EXIST=file_exists)
  if (file_exists) Falgo = 1
  INQUIRE(FILE=fic_char, EXIST=file_exists)
  if (file_exists) Fchar = 1
  INQUIRE(FILE=fic_numM, EXIST=file_exists)
  if (file_exists) Fnmat = 1
  INQUIRE(FILE=fic_numZ, EXIST=file_exists)
  if (file_exists) Fnzone = 1
  
  
end if !**********FIN DES AFFECTATIONS, DEBUTS TESTS 'COMMUNS'

  ! TEST la presence des fichiers xml
  !----------------------------------
  if(Fmat==-1) call amitex_abort("No .xml file for the 'material' properties",1,0)
  if(Falgo==-1) call amitex_abort("No .xml file for the 'algorithm' parameters",1,0)
  if(Fchar==-1) call amitex_abort("No .xml file for the 'loading'",1,0)
  if(fic_numZ /= "" .and. Fnzone==-1) call amitex_abort("No .vtk file for the 'zones' distribution",1,0)
  if(fic_numM /= "" .and. Fnmat==-1) call amitex_abort("No .vtk file for the 'materials' distribution",1,0)

  ! TEST cas sans fichier vtk
  !---------------------------------- 
  if(fic_numM=="" .and. fic_numZ=="")then
     if(nrank==0) then
        write(OUTPUT_UNIT,"(A)") "No vtk file given in the input :"
        write(OUTPUT_UNIT,"(A)") "      - a unique material is considered"
        write(OUTPUT_UNIT,"(A)") "      - one zone per voxel is considered"
        write(OUTPUT_UNIT,"(A,3I8)") "Numbers of voxels in each direction", nx, ny, nz
        write(OUTPUT_UNIT,"(A,3E15.8)") "Voxel size in each direction", dx, dy, dz
     end if
     if(nx==0 .OR. ny==0 .OR. nz==0) then
        call amitex_abort("No vtk file given in the input : &
                          & -NX, -NY, -NZ (-DX, -DY, -DZ) must be given in the command line",1,0)
     end if
  end if

  call check_amitex_abort(0)

end subroutine lire_commande
!==============================================================================

!==============================================================================
!
!>   SUBROUTINE DE COPIE DES FICHIERS XML
!!
!!
!! 
!! \param[in] fic_mat: (chaine de caracteres) nom du fichier de materiaux
!! \param[in] fic_algo: (chaine de caracteres) nom du fichier des parametres de l'algorithme
!! \param[in] fic_char: (chaine de caracteres) nom du fichier de chargement
!! \param[in] Fmat,Falgo,Fchar: (entier) unite logique des fichiers ci-dessus
!! \param[in] fic_vtk: (chaine de caracteres) radical du nom des fichiers vtk
!! \param[in] FWrite: (entier) unite logique utilise pour sauvegarder les fichiers
!! \param[in] Flog: (entier) unite logique du fichier log
!! \param[in] r : (entier) rang du proc realisant la copie
!!
!==============================================================================
subroutine copyXML(fic_mat, fic_algo, fic_char, fic_vtk,r)
  
  implicit none

  integer,intent(in)           :: r
  character(len=*),intent(in)  :: fic_mat, fic_algo, fic_char, fic_vtk
  integer                      :: cstat

  cstat=0

  if (nrank==r) then
    call copyfile(fic_mat ,trim(fic_vtk)//"_mat.xml" ,cstat)
    if (cstat /=0) write(Flog,"(A)") &
       "info:"//trim(fic_vtk)//"_mat.xml : pb with execute_command_line -> fortran copy"
    call copyfile(fic_algo,trim(fic_vtk)//"_algo.xml",cstat)
    if (cstat /=0) write(Flog,"(A)") &
       "info:"//trim(fic_vtk)//"_algo.xml : pb with execute_command_line -> fortran copy"
    call copyfile(fic_char,trim(fic_vtk)//"_char.xml",cstat)
    if (cstat /=0) write(Flog,"(A)") &
       "info:"//trim(fic_vtk)//"_char.xml : pb with execute_command_line -> fortran copy"
  end if

end subroutine copyXML

!===================================================================
subroutine copyfile(fic_in,fic_out,cstat)

  implicit none  

  character(len=*),intent(in)  :: fic_in, fic_out
  integer,intent(out)          :: cstat 
  integer                      :: Fin,Fout
  integer                      :: estat, i, j
  character(1000)              :: line, lineOS

  estat = 0
  cstat = 0

  !! Copy file with OS (execute_command_line)
  lineOS = "cp -f "//fic_in//" "//fic_out//achar(0)  ! adding //achar(0) to circumvent an intel15 bug
  call execute_command_line (trim(lineOS),exitstat=estat,cmdstat=cstat,wait=.true.)

  !! If not supported -> Copy file with read/write fortran
  if (cstat /= 0) then
    open(newunit=Fin, file=fic_in,form="formatted", status="old",action="read",iostat=j)  
    open(newunit=Fout, file=fic_out,form="formatted",status="replace", action="write",iostat=j)
    read(unit=Fin, fmt="(a)", iostat=i ) line
    do while ( i == 0 )
       write(unit=Fout, fmt="(a)",iostat=j) trim(line)
       read (unit=Fin, fmt="(a)", iostat=i) line
    end do
    close( unit=Fin )
    close( unit=Fout )
  end if

end subroutine copyfile


!===================================================================
!
!       SUBROUTINE GET_FILE
!>  Lit un nom de fichier en ligne de commande, retourne le nom du fichier
!!  et une unite de fichier s'il existe.
!! \param[in]       ind_arg: (entier) numero de l'argument
!! \param[in]       cur_file_u: (entier) derniere unite de fichier utilisee
!! \param[out]       cur_file_u: (entier) unite de fichier utilise
!! \param[out]       file_exist: (bool) existance du fichier
!! \param[out]       file_unit: (entier) unite de fichier associee au fichier en argument
!! \param[out]       file_name: (chaine de caracteres) nom du fichier en argument
!
!===================================================================
subroutine get_file( ind_arg, file_name, file_unit, cur_file_u, file_exist )

  implicit none

  integer, intent(in)           :: ind_arg
  integer, intent(inout)        :: cur_file_u
  logical, intent(out)        :: file_exist
  integer, intent(out)          :: file_unit
  character(len=*),intent(out)  :: file_name
  character(len=200)            :: arg !< arguments du code fortran lance en ligne de commande

  call getarg(ind_arg,arg)
  read(arg,"(A)") file_name
  INQUIRE(FILE=file_name, EXIST=file_exist)
  if(file_exist) then
     file_unit=cur_file_u
     cur_file_u=cur_file_u +1
  else
     file_unit=-1
  end if
end subroutine get_file


!===================================================================
!
!               SUBROUTINE SET _PROW_PCOL
!
!> Determination du nombre de lignes et de colonnes de la decomposition
!! a partir du nombre de processus utilises
!!
!! Cherche le couple de diviseur optimal (les plus proches) de nproc
!! pour réaliser la decomposition
!! 2decomp impose p_col < min(nx/2+1,ny), p_row < min(ny,nz)
!! La decomposition optimale calculee peut être non compatible avec 
!! ces conditions. 
!! Dans ce cas de figure, une decomposition 1D en nproc processus
!! sous reserve que nproc < min(nx/2+1,ny) ou < min(ny,nz)
!!
!! \param[out] p_row: (entier) nombre de lignes de la decomposition
!! \param[out] p_col: (entier) nombre de colonnes de la decomposition
!! \param[in]  nx,ny,nz : (entier) discretisation 
!===================================================================
subroutine set_pRow_pCol(p_row,p_col,nx,ny,nz)

  implicit none
  integer,intent(out) :: p_row,p_col
  integer,intent(in)  :: nx,ny,nz
  integer             :: p_row0=0,p_col0=0
  integer             :: n,ierror, d ,k
  integer             :: p_col_max,p_row_max
  character(len=200)  :: err_msg
  logical             :: file_exists=.false., affected=.false.,lopened=.false.
  integer             :: FID=7807 

  ! Initializations
  p_row = 0; p_col = 0
  call MPI_Comm_size(amitex_mpi_comm, n,ierror)  !nombre de processus utilises

  ! Check Max number of process
  p_row_max = min(nx/2+1,ny)
  p_col_max = min(ny,nz)

  if (n > p_row_max*p_col_max) then   
       write(err_msg,fmt="(A,I0,A,I0)") &
            "Number of proc ",n," too high, max possible = ", p_row_max*p_col_max
       call amitex_abort(err_msg,2,0)
  end if

  !-------------------Option -2decompp user : reading p_row,p_col in file optim_2decomp.txt
  if (decomp2D_user) then
    if (nrank==0) then
      inquire(FILE="user_2decomp.txt",EXIST=file_exists)
      if (File_exists) then
        open(unit=FID, file="user_2decomp.txt",form="formatted", action="read")
        read(FID,*) p_row0
        read(FID,*) p_col0
        close(FID)
      else
        call amitex_abort("Problem with option -2decomp user : file user_2decomp.txt is absent",-1,0)
      end if
    end if

    call MPI_Allreduce(p_row0, p_row, 1, MPI_INTEGER, MPI_MAX, amitex_mpi_comm, ierror)
    call MPI_Allreduce(p_col0, p_col, 1, MPI_INTEGER, MPI_MAX, amitex_mpi_comm, ierror)

    if (p_row*p_col .ne. n) then       
       write(err_msg,fmt="(A,I0,A,I0,A,I0)") &
            "(p_row, p_col) (",p_row,",",p_col,"), in user_2decomp.txt are not compatible with nproc ", n
       call amitex_abort(err_msg,-1,0)
    else
       affected = .true.
    end if 
  end if

  !-------------------Option -2decompp low_p_row : p_row as low as possible
  if (decomp2D_low_p_row) then

  ! Recherche dans l'ordre croissant de lignes
  do k=1,p_row_max
     p_row = k
     p_col = n / k
     if ((p_row * p_col == n) .and. (p_col .le. p_col_max)) then
        affected = .true.
        exit
     end if
  end do
    
  ! Si pas de couple trouve : recherche dans l'ordre croissant des colonnes
  if (.not. affected) then
  do k=1,p_col_max
     p_col = k
     p_row = n / k
     if ((p_row * p_col == n) .and. (p_row .le. p_row_max)) then
        affected = .true.
        exit
     end if
  end do
  end if

  ! Information si pas de couple trouve - on utilise option par defaut
  if (.not. affected) then
     write(err_msg,fmt="(A)") &
          "Problem when evaluating (p_row,p_col) with option -low_p_row : no couple found"
     call amitex_abort(err_msg,-1,0)
  end if
  end if

  !-------------------Option by default
  if (.not. affected) then
    p_col = int(sqrt(real(n)))
    p_row=p_col

    d=n-p_row*p_col
    do while( d /= 0)
       if(d>0)then
          p_col=p_col+((d-1)/p_row)+1
       else
          p_row=p_row+((d+1)/p_col)-1
       end if
       d=n-p_row*p_col
    end do

    if ((p_col > p_col_max) .or. (p_row > p_row_max)) then
       !! Decomposition 2D non acceptable par 2dcomp  avec le couple de valeurs 
       !! determinees par l'algorithme
       !! --> decomposition 1D 

       write(err_msg,fmt="(A,I0,A)") "Décomposition 2D de la cellule en ",n," processus incompatible &
            &avec 2decomp --> decomposition 1D"
       call amitex_abort(err_msg,-1,0)

       if (p_col_max >= n) then
          p_col = n
          p_row = 1
       elseif (p_row_max >= n) then
          p_row = n
          p_col = 1
       else
          ! Decomposition 1D impossible 
          write(err_msg,fmt="(A,I0)") "Decomposition 1D impossible, nombre de processus demande trop eleve :&
                &nombre max possible = ", max(p_row_max,p_col_max)
          call amitex_abort(err_msg,2,0)        
       end if
    end if
  end if


  if(nrank==0)then
    inquire(UNIT=Flog, OPENED=lopened)   ! if using set_prow_pcol before initializing Flog
    if (lopened) then
      write(Flog,"(A,I0,A)") "Simulation performed on ",n," process with :"
      write(Flog,"(I0,A)") p_row," lines"
      write(Flog,"(I0,A)") p_col," columns"
    end if
    write(OUTPUT_UNIT,"(A,I0,A)") "Simulation performed on ",n," process with :"
    write(OUTPUT_UNIT,"(I0,A)") p_row," lines"
    write(OUTPUT_UNIT,"(I0,A)") p_col," columns"
  end if

end subroutine set_pRow_pCol

!===================================================================
!
!               SUBROUTINE GET_AMITEXENV
!
!> Lecture de la variable d'environnement AMITEX_PATH 
!!       -> AMITEXenv%path et AMITEXenv%path_stat
!!
!!            
!!
!===================================================================
subroutine get_AMITEXenv()

  implicit none 
  
  integer :: vlength
  integer :: vstat

  vlength = 0
  vstat   = 0

if (.not. allocated(AMITEXenv%path)) then
  !lecture de la variable d'environnement AMITEX
  call GET_ENVIRONMENT_VARIABLE ( NAME ="AMITEX_PATH", LENGTH = vlength,STATUS=vstat)
  allocate ( character ( len = vlength ) :: AMITEXenv%path )

  if (vstat == 0) then
    call GET_ENVIRONMENT_VARIABLE ( NAME ="AMITEX_PATH", VALUE= AMITEXenv%path,STATUS=vstat)
    AMITEXenv%path_stat = vstat
  end if

  ! check error
  if (vstat==1) then
    call amitex_abort(&
    "Undefined environment variable AMITEX_PATH",0,0)
  else if (vstat==2) then 
    call amitex_abort(&
    "Undefined environment variable AMITEX_PATH (processor does not support env. variables)",0,0)
  else if (vstat .ne. 0) then
     call amitex_abort(&
    "Undefined environment variable AMITEX_PATH (unknown status from get_environment_variable)",0,0)
  else 
    call write_stdout0("AMITEX_PATH = "//AMITEXenv%path)
    call write_file0(Flog,"AMITEX_PATH = "//AMITEXenv%path)
  end if
end if  

end subroutine get_AMITEXenv

!===================================================================
!
! SUBROUTINE CHECK_LD_LIBRARY_PATH
!
!> Check LD_LIBRARY_PATH : FONCTION INUTILE CAR
!!
!!   si LD_LIBRARY_PATH n'est pas correct, 
!!         amitex_fftp plante des le lancement       
!!
!! (on garde comme exemple...
!!  peut-etre a modifier pour checker certaines coherences...)
!===================================================================
subroutine check_LD_LIBRARY_PATH()  ! gcc-warning accepte (unsused function)

  implicit none 
  
  integer                      :: vlength
  integer                      :: vstat
  character(len=:),allocatable :: ld_path
  integer                      :: pos

  vlength = 0
  vstat   = 0

  !lecture de la variable d'environnement LD_LIBRARY_PATH
  call GET_ENVIRONMENT_VARIABLE ( NAME ="LD_LIBRARY_PATH", LENGTH = vlength,STATUS=vstat)
  allocate ( character ( len = vlength ) :: ld_path )

  if (vstat == 0) then
    call GET_ENVIRONMENT_VARIABLE ( NAME ="LD_LIBRARY_PATH", VALUE= ld_path,STATUS=vstat)
  end if

  ! check error
  if (vstat==1) then
    call amitex_abort(&
    "Undefined environment variable LD_LIBRARY_PATH",2)
  else if (vstat==2) then 
    call amitex_abort(&
    "Undefined environment variable LD_LIBRARY_PATH (processor does not support env. variables)",2)
  else if (vstat .ne. 0) then
     call amitex_abort(&
    "Undefined environment variable LD_LIBRARY_PATH (unknown status from get_environment_variable)",2)
  end if
  
  ! check for user_functions in libAmitex/lib
  pos = index(ld_path,AMITEXenv%path//"/libAmitex/lib")

  if (pos==0) then
    call write_stdout0("*********** ERROR LD_LIBRARY_PATH, since version 8.17.2 :")
    call write_stdout0(AMITEXenv%path//"/libAmitex/lib must be added to LD_LIBRARY_PATH")
    call write_stdout0("depending on your linux shell, use :")
    call write_stdout0("    export LD_LIBRARY_PATH=your_path:$LD_LIBRARY_PATH  (bash/sh  shell)")
    call write_stdout0("    setenv LD_LIBRARY_PATH your_path:$LD_LIBRARY_PATH  (csh/tcsh shell)")

    call write_file0(Flog,"*********** ERROR LD_LIBRARY_PATH, since version 8.17.2 :")
    call write_file0(Flog,AMITEXenv%path//"/libAmitex/lib must be added to LD_LIBRARY_PATH")
    call write_file0(Flog,"depending on your linux shell, use :")
    call write_file0(Flog,"    export LD_LIBRARY_PATH=your_path:$LD_LIBRARY_PATH  (bash/sh  shell)")
    call write_file0(Flog,"    setenv LD_LIBRARY_PATH your_path:$LD_LIBRARY_PATH  (csh/tcsh shell)")

    call amitex_abort(AMITEXenv%path//"/libAmitex/lib must be added to LD_LIBRARY_PATH",2)

  end if
  
end subroutine check_LD_LIBRARY_PATH


!===========================================================================
!
!               SUBROUTINE INIT_LOG
!
!>    Open the log file from Flog and fic_log public variables
!!
!! \param[in] fic_log0 : (character) optional, log file name if not default
!!
!===========================================================================

subroutine init_log(fic_log0)

  implicit none 

  character(len=*), intent(in),optional :: fic_log0
  integer                               :: io_stat

  if (present(fic_log0)) fic_log=fic_log0

  if(nrank==0) then
     open(newunit=Flog, file=fic_log,form="formatted", status="replace", action="write",iostat= io_stat)  
     if ( io_stat /= 0 ) then
        write(OUTPUT_UNIT,"(3A,I0)")" Probleme a l'ouverture (fichier: ",fic_log,") (amitex)",io_stat
     end if
  end if

end subroutine init_log

end module amitex_mod
