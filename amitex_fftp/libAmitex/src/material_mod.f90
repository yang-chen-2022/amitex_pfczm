!==============================================================================
!       MODULE MATERIAL_MOD : 
!> Definition, initialisation et desallocation de la structure MATERIAL
!!
!! structure \n
!! MATERIAL :  structure definissant la cellule et ses proprietes \n
!!
!! subroutines\n 
!! read_mat :             initialisation des proprietes de la cellule: loi,
!!                        coefficients, variables internes \n
!! behavior :             comportement de la cellule (appel a UMAT) \n
!! initParam_umat :       initialisation des parametres de la fonction UMAT \n
!! deallocate_MattotP0 :  desallocation de MattotP0
!!
!! REMARQUE 1
!!  read_geom a ete sortie du module pour envisager de ne pas donner toutes
!!  les sources
!!  => necessite de mettre taille_mattotP en public
!!
!! REMARQUE 2
!!  pour l'utilisation de umat (compatibilité CAST3M, mfront etc...) :  
!!  on impose une utilisation de variables en double precision (voir behavior)
!! 

module material_mod

  use ISO_FORTRAN_ENV
  use IEEE_arithmetic

  use iso_c_binding

  use linear_mod
#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype, real_type
  use decomp_2d_mpi,       only : nrank, nproc
#else
  use decomp_2d, only : mytype, nrank, nproc, real_type
#endif  

  use mpi
  use Fox_dom
  use io_amitex_mod
  use error_mod
  use amitex_mod, only         : AMITEXenv          ! variable
  use field_mod, only          : amx_size           ! pointer
  use param_algo_mod , only    : algo_param0, &! initialization variable (used in read_mat, read_mat_coposite)
                                 algo_param    ! pointer
  use algo_functions_mod, only : act3local ! function
  use amitex_mpi_mod, only: amitex_mpi_comm
#ifdef OPENMP
  use omp_lib
#endif


  private

  !> Pointer "publiques" vers composantes de SIMU_AMITEX
  public ::     nmateriaux, Matref, MattotP, VCinfo, testComposite, times_m,&
                Nloc_models

  !> Variables "publiques" utilisees lors de l'initialisation
  public ::     nmateriaux0, Matref0, MattotP0,Nloc_models0,&
                MatComposite, Interphases 

  !> Variables publiques a transformer en pointer (ou variables 'provisoires')
  !RAS

  !> Types publiques (pour definition de SIMU_AMITEX)
  public ::     MATERIAL, MATERIAL_REFERENCE, NMATERIAL, INFO_VOXCOMP, TIMES_MAT, NL_MODEL

  !> Fonctions publiques
  public ::     print_material_structure,&
                eval_nb_materiaux,&
                read_mat,&
                behavior,behaviorD,&
                test_composite,&
                get_interphase_list,&
                read_composite,&
                read_mat_composite,&
                assemble_field_from_varint,assemble_field_from_varint0,&
                assign_varint_from_field,assign_varint0_from_field,&
                get_homogeneous_coeff, get_minmax_coeff, get_minmax_coeff_function, &
                multiply_field_coeff,apply_field_coeff_function,&
                CalcRotMatrices,&
                CalcRotMatrices_GD,&
                deallocate_MattotP0,deallocate_Matref0,&
                deallocate_MatComposite,&
                deallocate_Interphases,&
                MIN_MPI_I8,MAX_MPI_I8,MAX_MPI_I,SUM_MPI_I8,SUM_MPI_I8_LARGE,&
                umat_load, umat_call,initParam_umat ! derniere ligne publique pour tests elementaires

!!------------------------------------------------------------------------------
!>                                                            TIMES MEASUREMENTS

  type TIMES_MAT
     double precision :: behavior = 0                     !< times in behavior/behaviorD
  end type TIMES_MAT

!!------------------------------------------------------------------------------
!>                                                    DENOMBREMENT DES MATERIAUX

  type NMATERIAL
     integer :: n                                         !< Nombre de materiaux dans la cellule
     integer(kind=8), allocatable, dimension(:) :: nZ     !< Nombre de zones par materiau
     integer :: n_loc                                     !< Nombre de materiaux sur le pinceau
     integer :: n_non_composites                          !< Nombre de materiaux non composites dans la cellule
     integer :: n_non_composites_loc                      !< Nombre de materiaux non composites sur le pinceau
     integer :: n_composites                              !< Nombre de materiaux composites dans la cellule
     integer :: n_composites_loc                          !< Nombre de materiaux composites sur le pinceau  
  end type NMATERIAL


!!------------------------------------------------------------------------------
!>                                                         MATERIAU DE REFERENCE

  type MATERIAL_REFERENCE
     real(mytype),  dimension(2)   :: LambdaMu0            !< [lambda0,Mu0]
     real(mytype),allocatable,dimension(:,:) :: C0         !< C0 : matrice de rigidite du materiau de reference
                                                           !< (6,6) en HPP ou (9,9) en GD
     real(mytype),allocatable,dimension(:)   :: k0D!,A0D   !< k0D - coefficient de conductivité (diffusion) 
                                                           !< taille (nVarD)
                                                           !< TODO : A0D - coeff. associe a la partie instationnaire
  end type MATERIAL_REFERENCE


!> Structure Materiau valable pour les materiaux "homogenes" et "composites"
!!      un "materiau" composite (mattotP(i)) étant défini par :
!!              - le modèle composite  : défini dans libname, lawname
!!              - chaque phase associée à un materiau homogene : 
!!                      définis dans numM_composite et numZ_composite
!!              - les coeff propres au modèle composite (fractions volumique, normale à l'interface...) 
!!                      + les coeff. de chaque phase
!!              - les variables internes propres au modèle composite (deformation AP, contrainte plane...) 
!!                      + les variables internes propres au modèle
!!
!!             REMARQUE :
!!                      dans le cas d'une utilisation de materiaux "composite" pour la prise en compte
!!                      d'interfaces (interphases) (i.e. voxels composites) on a 1 zone par voxel
!!                      => nombre de zone du materiau = nombre de voxels du matériau
!!
!!             IMPORTANT : TRAITEMENT DES INTERPHASES
!!                      La structure permet la prise en compte de materiaux ou de zones dits "interphase"
!!                      => aucun voxel de la cellule n'est constitue de ces materiaux ou de ces zones 
!!                      => ces materiaux/zones "interphase" sont present comme phases des materiaux composites,
!!                         leurs variables internes/coefficients sont donc contenus dans les tableaux 
!!                         %VarInt et %Coeff des materiaux composites
!!                      => ces materiaux / zones interphase disposent d'emplacement specifiques dans mattotP
!!                         pour y lire leurs variables internes et coefficients depuis le mat.xml
!!                         ATTENTION : ces emplacements sont utilises uniquement pour stocker ces donnees a la 
!!                                     lecture des fichiers xml. 
!!                                     Ces donnes seront par la suite copie au sein d'un materiau composite, 
!!                                     et seront utilises dans les calculs du comportement uniquement depuis
!!                                     ce materiau composite
!!                      => Les donnees specifiques "interphase" initialisees a la lecture des donnees 
!!                         sont stockees dans mattotP :
!!                                pour un materiau "interphase" : un element de mattotP
!!                                pour des zones "interphase" : des lignes ajoutees aux tableaux
!!                                            mattotP%Zone, mattotP%Coeff et mattotP%Varint
!!                      => PAR CONVENTION, les donnees specifiques "interphase" (voir ci-dessus) sont stockees 
!!                         UNIQUEMENT SUR LE PINCEAU "0".
!!                      => Les donnees specifiques "interphase" de mattotP sont utilisees dans le code 
!!                         UNIQUEMENT A L'INITIALISATION (pour definir les materiaux composites)
!!
!!                         A L'ISSUE DE L'INITIALISATION LES DONNEES "INTERPHASE" NE SONT PLUS UTILISEES,
!!                         ELLES ONT ETE REPARTIES AU SEIN DES MATERIAUX COMPOSITES
!!
!!             TODO :  AFIN D'EVITER DES RISQUES DE CONFUSION,
!!                     APRES INITIALISATIONS "call clean_interphase" 
!!                          -> SUPPRIMER DU PINCEAU 0  :
!!                                     LES MATERIAUX D'INTERPHASE de mattotP
!!                                     LES ZONES D'INTERPHASE DE CHAQUE MATERIAU
!!
!!             AMELIORATION POSSIBLE : 
!!                      la structure Materiau definie ci-dessous implique de memoriser 
!!                      les coefficients de chaque zone meme si certains sont constants sur 
!!                      toutes les zones du materiau
!!
!!                      IDEE : définir un tableau de taille variable avec pour chaque coeff "i"
!!                                     coeff(i)%val de dimension soit 1 soit Nzone
!!                             définir un tableau indicateur (dimension Ncoeff) : coeff_constant
!!                                     liste des indices des coefficients constants
!!
!!

  type MATERIAL
    !> Nombre de phase : si Nphase=1 cas "homogene"; si Nphase > 1 :cas "composite"
    integer                                     :: Nphase = 1
    !> Balise indiquant si le matériau est un matériau d'interphase (.true. ou .false.)
    !>        Materiau d'interphase : materiau n'existant pas en tant que voxel "homogene"
    !>        ATTENTION : par convention, les materiaux d'interphase sont sur le noeud 0
    logical                                     :: Interphase = .false.
    !> Liste des zones d'interphase si le matériau en contient
    !>        Zone d'interphase : zone n'existant pas en tant que voxel "homogene"
    integer(kind=8),allocatable,dimension(:)    :: Zones_interphase
    !> numero du materiau
    integer                                     :: numM     
    !> material name (optional, and only defined if Nphase=1)
    character(len=100)                          :: Matename     
    !> LOI MECANIQUE associe au materiau
    !> nom de la librairie (chemin complet, sensible a la casse, max 200 caracteres)
    character(len=200,kind=C_CHAR)              :: Libname
    !> nom de loi (converti en minuscule -> insensible a la casse, max 200 caracteres)
    !>    Cas composite : "voigt" ou "laminate" ou "reuss"
    character(len=200,kind=C_CHAR)              :: Lawname     
    !> LOI DIFFUSION associe au materiau 
    !>               Terme de DIFFUSION K (q=-k grad)
    character(len=200,kind=C_CHAR)              :: LibnameK     
    character(len=200,kind=C_CHAR)              :: LawnameK     
                !TODO : cas terme source + instationnaire
                !    !>               Terme INSTATIONNAIRE A (div q = S-Ad./dt)
                !    character(len=200,kind=C_CHAR)              :: LibnameA     
                !    character(len=200,kind=C_CHAR)              :: LawnameA     
                !    !>               Terme SOURCE S (div q = S-Ad./dt)
                !    character(len=200,kind=C_CHAR)              :: LibnameS     
                !    character(len=200,kind=C_CHAR)              :: LawnameS     

    !> nombre de coefficients    
    integer                                     :: Ncoeff
    integer                                     :: NcoeffK
                !    integer                                     :: NcoeffA
                !    integer                                     :: NcoeffS
                
    !> names of coefficients and internal variables : coeff_name(Ncoeff), varint_name(Nvarint) 
    !        (optional and only defined if Nphase=1)
    !>       (usefull for multi-physic extensions)
    character(len=100), allocatable, dimension(:)  :: Coeffname     
    character(len=100), allocatable, dimension(:)  :: Varintname     

    !> nombre de variable internes 
    integer                                     :: NvarInt  
    !> position linearisee des voxels dans le pinceau, ordonnes par numero 
    !!    de zone croissant
    !! ATTENTION : %pos -  est vide pour un materiau d'interphase
    !!                  -  ne contient pas d'élements pour les zones interphase
    !!                     par definition
    integer(kind=8),allocatable,dimension(:)    :: pos 
    !> Zone(i,1): indice dans pos du dernier element de la zone 'Zone(i,2)'
    !! Zone(i,2): numZ associe (necessaire pour moyenne, ecart types ...)
    !!   Pour une zone d'interphase,  
    !!       uniquement sur le pinceau 0, on ajoute en fin de tableau autant de lignes que de zones interphase
    !!       pour ces lignes : Zone(I_hom+i,1) : Zone(I_hom,1) + i
    !!                                          --> cet element ne correspond pas a un element de %pos 
    !!                                              mais a une colonne de %VarInt sur le pinceau 0 (voir %VarInt)
    !!                                              * Lors de l'initalisation des variables internes, cet element est lu  
    !!                                                et indique l'emplacement de la colonne de %VarInt ou stocker les donnees
    !!                                                dans %VarInt lors de la lecture des variables internes pour le materiau homogene
    !!                                              * Lors de l'initialisation des materiaux composites, cet element est lu 
    !!                                                et indique l'emplacement des variables internes a recuperer dans le tableau 
    !!                                                %VarInt
    !!                         Zone(I_hom+i,2) : numero de zone (egalement contenu dans %Zones_Interphase)
    !!   Pour un materiau "INTERPHASE":
    !!       Zone(i,1) : i 
    !!                   --> cet element ne correspond pas a un element de %pos  mais a une colonne de %VarInt
    !!                                              * Lors de l'initalisation des variables internes, cet element est lu  
    !!                                                et indique l'emplacement de la colonne de %VarInt ou stocker les donnees
    !!                                                dans %VarInt lors de la lecture des variables internes pour le materiau homogene
    !!                                              * Lors de l'initialisation des materiaux composites, cet element est lu 
    !!                                                et indique l'emplacement des variables internes a recuperer dans le tableau 
    !!                                                %VarInt
    !!       Zone(i,2) : numero de zone
    integer(kind=8),allocatable,dimension(:,:)  :: Zone
    !> coefficients, dimensions: (nCoeff, nZone)
    real(mytype),allocatable, dimension(:,:)    :: Coeff
    real(mytype),allocatable, dimension(:,:)    :: CoeffK
                !    real(mytype),allocatable, dimension(:,:)    :: CoeffA
                !    real(mytype),allocatable, dimension(:,:)    :: CoeffS
    !> variables internes :
    !>     dans le cas de modeles non-locaux, les variables "non-locales" (GNloc)
    !>     sont stockees a la suite des variables locales dans l'ordre des indices de modeles non-locaux 
    !>     (NlocMod_num dans param_mat.xml) 
    !>                    dimensions (nVarInt,NPosition) 
    !>                            ou (nVarint,nZones) pour un materiau "INTERPHASE" 
    !>                            ou sur le pincau 0 si presence de zones interphases
    !>                                (nVarint, Nposition+Nzones_interphases)   
    !>                                les indices des Nzones_interphases colonnes du tableau 
    !>                                correspondent dans ce cas aux valeurs dans les Nzones_interphases  
    !>                                lignes de %Zone(:,1)
    real(mytype),allocatable, dimension(:,:)    :: VarInt
    !> variables internes dimensions: (nVarInt,NPosition)
    !>                             ou (nVarint,nZones) pour un materiau "INTERPHASE" 
    !>                             ou sur le pincau 0 si presence de zones interphases
    !>                                (nVarint, Nposition+Nzones_interphases)   
    !>                                les indices des Nzones_interphases colonnes du tableau 
    !>                                correspondent dans ce cas aux valeurs dans les Nzones_interphases  
    !>                                lignes de %Zone(:,1)
    real(mytype),allocatable, dimension(:,:)    :: VarInt0


    !> Entrees supplementaires allouees dans le cas composite (Nphase>1)
    !!------------------------------------------------------------------
    !> liste des numeros de materiau associé à chacune des phases
    !!          dimension (Nphase)
    integer, allocatable,dimension(:)           :: numM_composite
    !> liste des numeros de zone du materiau associé à chacune des phase
    !!          numZ_composite(i,j)=numero de la zone GLOBALE du materiau numM_composite(i) pour la zone "j" LOCALE du matériau composite
    !!          dimension (Nphase,nZone)
    integer(kind=INT64), allocatable,dimension(:,:)          :: numZ_composite
    !> tableau de noms de librairie pour chacune des phases (chemin complet, sensible a la casse, max 200 caracteres) (NPhases)
    character(len=200,kind=C_CHAR), allocatable,dimension(:) :: Libname_comp
    !> tableau de noms de loi pour chacune des phases (chemin complet, sensible a la casse, max 200 caracteres) (NPhases)
    character(len=200,kind=C_CHAR), allocatable,dimension(:) :: Lawname_comp
    !>
    !> Composantes ci-dessous utiles pour l'instant (v6.2.0) dans read_mat_composite 
    !!               pour construction de %Coeff et %Varint
    !!                                                 => interet de les conserver ?
    !!------------------------------------------------------------------------------
    !> coefficients pour loi composite (laminate et reuss pour l'instant) dimensions : (nCoeff,nZone)
    !> la taille 1 doit être cohérente avec la loi composite (2 lignes pour laminate et reuss,0 pour Voigt)
    !> il s'agit des "Coeff_composites" donnes dans les fichiers materiaux .xml
    real(mytype),allocatable, dimension(:,:)    :: Coeff_comp
    !> nombre de coefficients matériau de chacune des phases (dimension Nphase)
    integer,allocatable,dimension(:)            :: NCoeff_comp
    !> nombre de variable internes de chacune des phases (dimension Nphase)
    integer,allocatable,dimension(:)            :: NvarInt_comp

    !> Parametres d'algorithmes pour les materiaux composites, utiles pour unifier le calcul du 
    !  comportement des voxels composites 
    !---------------------------------------------------------------------------------
     !> Acceleration de convergence
     logical                                 :: acc_CV_composite
     !> Tolerance sur le critere d'equilibre (a multiplier par le critere global tol_criteq)
     !>                    (par defaut 0.01)
     real(mytype)                            :: tol_criteq_composite
     !> Type d'initialisation (default,proportional, linear)
     !> "default" : increments de deformation homogene (identique sur toutes les phases)
     !> "proportional" : inrements de def. conservant la proportionalite des deformations entre phases
     !> "linear" : increments de deformation interpoles a partir du pas precedent
     character(len=16)                       :: Init_composite
     !> Paramètres de pilotage
     integer(kind=4)                         :: Nmax_subdivision_composite
     !> Type de matrice tangente utilisé (elastic, mfront_elastic, mfront_tangent, numerical)
     !> "elastic"        : comportement linéaire isotrope par phase donne par Coeff_composite dans materiau.xml
     !> "mfront_elastic" : tenseur d'élasticité renvoyee par Mfront
     !> "mfront_tangent" : matrice tangente cohérente renvoyee par Mfront
     !> "mfront_secant"  : matrice secant (elastique endommagee) renvoyee par Mfront 
     !                     si codee dans la loi
     !> "numerical"     : matrice tangente obtenue par differentiation numerique
     character(len=16)                       :: Jacobian_type_composite
     !> Valeur absolue de la perturbation pour le calcul de la matrice tangente numerique
     !  du modèle laminate
     real(mytype)                            :: Perturbation_composite
     contains 
       GENERIC   :: WRITE(UNFORMATTED) => write_mate
       GENERIC   :: READ(UNFORMATTED)  => read_mate
       procedure :: write_mate
       procedure :: read_mate

  end type MATERIAL

!------------------------------------------------------------------------------
!> Structure composite valable uniquement pour les matériaux composites
!!                  Structure TEMPORAIRE permettant de récupérer les données composites
!!                  facilement à partir des fichiers d'entrée avant de les réintroduire
!!                  dans le structure MATERIAL
!!      
!!    7/12/2020 : on supprime les champs stockes dans %pos_global, Zone, Fv, Normale, 
!!                                                      Direction, Surfaces
!!                 on rajoute %Nzones_global et %dir pour lire ces champs au moment requis
!!
  type COMPOSITE
   !> Tableau contenant le numéro des phases de chacun des matériaux
   !! composant le voxel composite
     integer,allocatable,dimension(:)             ::  num_Phases

   !> nom de loi composite ('voigt','laminate','reuss' pour l'instant uniquement)
     character(len=200,kind=C_CHAR)               ::  Lawname

   !> nom du repertoire de stockage des donnees composites
     character(len=400)                           ::  dir
     
   !> nombre de zone total = nombre total de voxels composites sur la cellule (rappel :1 zone = 1 voxel)
     integer(INT64)                               ::  Nzones_global

  end type COMPOSITE

!------------------------------------------------------------------------------
!> Structure INTERPHASE valable uniquement pour les matériaux composites
!!                  Structure TEMPORAIRE permettant d'identifier les "Interphases"
!!                  avant de les réintroduire dans le structure MATERIAL
!!                  Initialisée par get_interphase_list à partir du fichier mat.xml
!! "Interphase" :   materiau ou zone n'existant pas sous forme de voxel "homogene"

  type INTERPHASE
   !> Tableau contenant le numéro des phases de chacun des matériaux
   !! composant le voxel composite
     integer(kind=8)                              ::  numM    ! numéro du matériau concerné
     logical                                      ::  all     ! vrai si tout le matériau est une interphase
                                                              ! faux si seulement des zones du matériau sont
                                                              ! des interphases
     integer(kind=8),allocatable,dimension(:)     ::  zones   ! alloué si des zones du matériau sont des 
                                                              ! interphases
                                                              ! contient la liste des numéros de zone 
                                                              ! correspondants
     integer(kind=8)                              :: Nzones   ! nombre de zones du matériau d'interphases
                                                              ! non utilise si all = .true.
                                                              ! (uniquement pour les materiaux d'interphase)
  end type INTERPHASE

!------------------------------------------------------------------------------
!> NL_Model : Sous-structure de param_algo de description des modeles 
!>            non locaux 
!>          
!>         Les modeles non locaux sont codes actuellement dans le module
!>         non_local_mod. Le role de ces routines est de calculer 
!>         des variables internes non locales (gradients etc...)
!>         a partir de variables internes locales.
!>
!>         materiaux concernes par un modele = materiaux dont les variables
!>         internes sont utilisees par le modele non local et/ou 
!>         recuperant des variables internes non locales necessaires
!>         au calcul de son comportement mecanique
!------------------------------------------------------------------------------
type NL_MODEL   
     !> Nom du modele non local, passe a NL_call pour 
     !> execution de la routine correspondante
     character(len=200)                         :: Modelname 
     !> Numero d'identification du modele non local
     !integer                                    :: NlocMod_num !=> rendu coherent avec indice du tableau Nloc_models 

     !> Nombre de composantes scalaires des  variables internes locales 
     !> utilisees par le modele non local pour evaluer les variables 
     !> internes non locales
     integer                                    :: Nnloc
     !> Nombre de composantes scalaires des variables internes non
     !> locales evaluees par le modele non local (variables "Gradient")
     integer                                    :: NGnloc
     !> Nombre de coefficients du modele non-local a rechercher dans les coefficients 
     !> des materiaux associes au modelemodele non local (variables "Gradient")
     integer                                    :: Ncoeff_nloc
     !> Numeros d'identification des materiaux concernes par le modele
     !>                       taille (nmat)
     integer,allocatable,dimension(:)           :: numM
     !> Tableau contenant pour chaque materiau les indices des variables 
     !> internes locales necessaires aux calculs du modele non local NlocMod_num
     !>                       taille (nmat,Nnloc)
     !>                       variables a recuperer comme suit :
     !>                       j tel que: mattoP(j)%numM = Nloc_models%numM(i) 
     !>                       mattotP(j)%VarInt(Ind_VarNloc(i,:)) 
     integer,allocatable,dimension(:,:)         :: Ind_VarNloc
     !> Tableau contenant pour chaque materiau les indices des variables 
     !> internes non locales issues des calculs du modele non local NlocMod_num
     !>                       taille (nmat,Nnloc)
     !>                       variables a stocker comme suit : 
     !>                       j tel que: mattoP(j)%numM = Nloc_models%numM(i) 
     !>                       mattotP(j)%VarInt(Ind_VarGNloc(i,:)) 
     integer,allocatable,dimension(:,:)         :: Ind_VarGNloc
     !> Tableau contenant pour chaque materiau les indices des coefficients 
     !> necessaires aux calculs du modele non local NlocMod_num
     !>                       taille (nmat,Ncoeff_nloc)
     !>                       variables a recuperer comme suit : 
     !>                       j tel que: mattoP(j)%numM = Nloc_models%numM(i) 
     !>                       mattotP(j)%Coeff(Ind_CoeffNloc(i,:),:) 
     integer,allocatable,dimension(:,:)         :: Ind_CoeffNloc

end type NL_MODEL

!!------------------------------------------------------------------------------
!>                                    SUIVI DE L'EXECUTION ALGO VOXELS COMPOSITES

type INFO_VOXCOMP
  integer(kind=8)  :: nIt_lam=0      !< nombre total d'iterations de l'algorithme multi-couches
                                        !! (~ nbre de voxel x nbre d'appels par voxel) 
  integer(kind=8)  :: nSub_lam=0     !< nombre total de recours à une subdivision du pas de temps 
                                        !! (~ nbre de voxel x nbre de subdivisions par voxel)
  integer(kind=8)  :: nIncr_lam=0    !< nombre total de sous pas de temps total utilisé par PilotageLaminate
                                        !! (~ nbre de voxel x nbre de sous-pas voxel)
  integer(kind=8)  :: nCall_lam=0    !< nombre total d'appels à umatLaminate au cours du calcul
                                        !! (~ nbre de voxel x nbre d'appel par voxel)
                                        
  integer(kind=8)  :: nIt_lam_pas     !< IDEM mais sur un pas de chargement
  integer(kind=8)  :: nSub_lam_pas
  integer(kind=8)  :: nIncr_lam_pas
  integer(kind=8)  :: nCall_lam_pas         
  
end type INFO_VOXCOMP



!> structures 'NMATERIAL'
  type(NMATERIAL)                           :: Nmateriaux0
  type(NMATERIAL), pointer                  :: Nmateriaux
  
!> structures 'MATERIAL_REFERENCE'
  type(MATERIAL_REFERENCE)                  :: Matref0
  type(MATERIAL_REFERENCE), pointer         :: Matref
  
!> tableau de structures 'MATERIAL'
  type(MATERIAL), allocatable,dimension(:)  :: MattotP0
  type(MATERIAL), pointer,dimension(:)      :: MattotP

!> tableau de structures 'COMPOSITE'
  type(COMPOSITE),allocatable,dimension(:)  :: MatComposite

!> tableau de structures 'INTERPHASE'
  type(INTERPHASE),allocatable,dimension(:) :: Interphases

!> tableau de structures 'NL_Model'
  type(NL_Model),allocatable,dimension(:)   :: Nloc_models0
  type(NL_Model),pointer,dimension(:)       :: Nloc_models

!> structure 'INFO_VOXCOMP'
  type(INFO_VOXCOMP),pointer                :: VCinfo
  
!> structure 'TIMES_MAT'
  type(TIMES_MAT),pointer                   :: times_m
  
!> Test de presence de voxels composites
  logical, pointer                          :: testComposite

! ATTENTION :  MattotP est parallelisé, pour chaque proc,
!              sa taille est égale au nombre de materiaux presents dans le pinceau!
!-------------------------------------------------------------------------------


! parametres UMAT - fixes sur le pas de temps
!--------------------------------------------
  character(len=16)                   :: CMNAME
  integer(kind=INT64)                 :: NDI, NSHR, NOEL, NPT, layer, kspt, KSTEP, NTENS
  double precision,dimension(6)       :: ddsddt, drplde
  double precision,dimension(6,6)     :: ddsdde
  double precision                    :: sse, spd, scd, rpl, drpldt, TEMP,DTEMP, &
                                         pnewdt, celent
  double precision,allocatable, dimension(:) :: PREDEF, DPRED
  double precision, dimension(3)      :: COORDS
  double precision, dimension(3,3)    :: drot


!> Indices d'identification dans le repère d el'interface (modèle Multicouches)
  integer, dimension(3)                                   :: IndAP = (/1,4,5/)
  integer, dimension(3)                                   :: IndPl = (/2,3,6/)


!------------------------------------------------------------------------------
!> Gestion de l'interface umat

!> interface umat_load
  interface
     function umat_load(l,n) bind(c,name="umat_load")
       import :: c_ptr
       import :: c_char
       implicit none
       type(c_ptr) :: umat_load
       character(kind=c_char),  intent(in) :: l(*)
       character(kind=c_char),  intent(in) :: n(*)
     end function umat_load
  end interface

!> interface umat_call
!MODIF LG : on remplace INTEGER*8 par INTEGER(KIND=c_long_long)
!MODIF LG 7/8/2018: on remplace INTEGER*8 par INTEGER(KIND=c_int64_t)
  interface
     subroutine umat_call(f,STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL, &
       DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,   &
       PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,           &
       NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,            &
       NOEL,NPT,LAYER,KSPT,KSTEP,KINC) bind(c,name="umat_call")
       import :: c_ptr
       import :: c_double
       import :: c_char
       import :: c_int
       import :: C_INT64_T
       implicit none
       type(c_ptr), value,    intent(in)  :: f
       REAL(KIND=c_double),   intent(out) :: STRESS(*)
       REAL(KIND=c_double),   intent(out) :: STATEV(NSTATV)
       REAL(KIND=c_double),   intent(in) :: DDSDDE(NTENS,NTENS)
       REAL(KIND=c_double),   intent(in) :: SSE
       REAL(KIND=c_double),   intent(in) :: SPD
       REAL(KIND=c_double),   intent(in) :: SCD
       REAL(KIND=c_double),   intent(in) :: RPL
       REAL(KIND=c_double),   intent(in) :: DDSDDT(NTENS)
       REAL(KIND=c_double),   intent(in) :: DRPLDE(NTENS)
       REAL(KIND=c_double),   intent(in) :: DRPLDT
       REAL(KIND=c_double),   intent(in)  :: STRAN(*)
       REAL(KIND=c_double),   intent(in)  :: DSTRAN(*)
       REAL(KIND=c_double),   intent(in)  :: TIME(*)
       REAL(KIND=c_double),   intent(in)  :: DTIME
       REAL(KIND=c_double),   intent(in)  :: TEMP
       REAL(KIND=c_double),   intent(in)  :: DTEMP
       REAL(KIND=c_double),   intent(in)  :: PREDEF(*)
       REAL(KIND=c_double),   intent(in)  :: DPRED(*)
       CHARACTER(KIND=c_char),intent(in)  :: CMNAME(*)
       INTEGER(KIND=C_INT64_T),   intent(in)  :: NDI
       INTEGER(KIND=C_INT64_T),   intent(in)  :: NSHR
       INTEGER(KIND=C_INT64_T),   intent(in)  :: NTENS
       INTEGER(KIND=C_INT64_T),   intent(in)  :: NSTATV
       REAL(KIND=c_double),   intent(in)  :: PROPS(NPROPS)
       INTEGER(KIND=C_INT64_T),   intent(in)  :: NPROPS
       REAL(KIND=c_double),   intent(in)  :: COORDS(3)
       REAL(KIND=c_double),   intent(in)  :: DROT(3,3)
       REAL(KIND=c_double),   intent(in) :: PNEWDT
       REAL(KIND=c_double),   intent(in) :: CELENT
       REAL(KIND=c_double),   intent(in) :: DFGRD0(3,3)
       REAL(KIND=c_double),   intent(in) :: DFGRD1(3,3)
       INTEGER(KIND=C_INT64_T),   intent(in)  :: NOEL
       INTEGER(KIND=C_INT64_T),   intent(in)  :: NPT
       INTEGER(KIND=C_INT64_T),   intent(in)  :: LAYER
       INTEGER(KIND=C_INT64_T),   intent(in)  :: KSPT
       INTEGER(KIND=C_INT64_T),   intent(in)  :: KSTEP
       INTEGER(KIND=C_INT64_T),   intent(in)  :: KINC
     end subroutine umat_call
  end interface

!> interface umatD_call
  interface
     subroutine umatD_call(f,FLUXD,GRADQD,DGRADQD,TIME,DTIME,TEMP,DTEMP,   &
       PREDEF,DPRED,NVARD,PROPS,NPROPS,KINC) bind(c,name="umatD_call")
       import :: c_ptr
       import :: c_double
       import :: c_char
       import :: c_int
       import :: C_INT64_T
       implicit none
       type(c_ptr), value,    intent(in)  :: f
       REAL(KIND=c_double),   intent(out) :: FLUXD(3,NVARD)
       REAL(KIND=c_double),   intent(in)  :: GRADQD(3,NVARD)
       REAL(KIND=c_double),   intent(in)  :: DGRADQD(3,NVARD)
       REAL(KIND=c_double),   intent(in)  :: TIME(*)
       REAL(KIND=c_double),   intent(in)  :: DTIME
       REAL(KIND=c_double),   intent(in)  :: TEMP
       REAL(KIND=c_double),   intent(in)  :: DTEMP
       REAL(KIND=c_double),   intent(in)  :: PREDEF(*)
       REAL(KIND=c_double),   intent(in)  :: DPRED(*)
       INTEGER(KIND=C_INT64_T),   intent(in)  :: NVARD
       REAL(KIND=c_double),   intent(in)        :: PROPS(NPROPS)
       INTEGER(KIND=C_INT64_T),   intent(in)  :: NPROPS
       INTEGER(KIND=C_INT64_T),   intent(in)  :: KINC
     end subroutine umatD_call
  end interface


!> Function interface for use in get_minmax_coeff_function
!---------------------------------------------------------
 abstract interface 
   function func0(coeff0) result(values0)
      real(kind(0.0D0)), dimension(:,:),intent(in) :: coeff0  ! ATTENTION - gcc n'accepte pas real(mytype) dans l'interface
      real(kind(0.0D0)), dimension(size(coeff0,2)) :: values0 !             on utilise la definition de mytype dans 2decomp (kind(0.0D0))
   end function func0
 end interface
 
!> Function interface for use in apply_field_coeff_function
!
! param[in] : Field0  : dimension(nzone,ncomp,nvar) part (zone) of a field
! param[in] : coeff0  : list of coefficients of the non-local model (homogeneous on the zone)
! 
!---------------------------------------------------------------------
 abstract interface 
   function func1(Field0,coeff0) result(Field1)
      real(kind(0.0D0)), dimension(:,:,:),intent(in) :: Field0  ! ATTENTION - gcc n'accepte pas real(mytype) dans l'interface
      real(kind(0.0D0)), dimension(:),intent(in)     :: coeff0  !             on utilise la definition de mytype dans 2decomp (kind(0.0D0))
      real(kind(0.0D0)), dimension(size(Field0,1),size(Field0,2), size(Field0,3)) :: Field1
   end function func1
 end interface

contains


subroutine write_mate(MATE, iunit, istat, imsg)

  implicit none
  class(MATERIAL),intent(in)       :: MATE  
  integer,intent(in)               :: iunit
  integer,intent(out)              :: istat
  character(len=*), intent(inout)  :: imsg
  
  !if (allocated(MATE%varint))  write(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) MATE%varint
  if (allocated(MATE%varint0)) write(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) MATE%varint0

  if (nrank == 0 ) write(OUTPUT_UNIT,*) "    WRITE_MATE : iostat =", istat

end subroutine write_mate

subroutine read_mate(MATE, iunit, istat, imsg)

  implicit none
  class(MATERIAL),intent(inout)    :: MATE  
  integer,intent(in)               :: iunit
  integer,intent(out)              :: istat
  character(len=*), intent(inout)  :: imsg
  
  !if (allocated(MATE%varint))  read(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) MATE%varint
  if (allocated(MATE%varint0)) read(UNIT=iunit,IOSTAT=istat,IOMSG=imsg) MATE%varint0
  if (nrank == 0 ) write(OUTPUT_UNIT,*) "    READ_MATE : iostat", istat


end subroutine read_mate

!==================================================================================
!==================================================================================
!                         SUBROUTINE PRINT_MATERIAL_STRUCTURE
!
!> Ecriture des composantes de la structure materiau
!!
!!
!!  \param[in] Flog:   (entier) unite logique du fichier de sortie .log
!!  \param[in] nrank0: (entier) pinceau pour lequel on effectue la sortie
!!
!==================================================================================
subroutine print_material_structure(Flog,nrank0,nx,ny,nz)
  
  implicit none
  integer, intent(in)   :: Flog
  integer, intent(in)   :: nrank0
  integer, intent(in)   :: nx, ny, nz                   !< dimensions de la cellule
  integer               :: i,j,k
  integer(KIND=8)       :: ntot_comp
  
  ntot_comp = 0
  if (nrank==nrank0) then
     write(Flog,"(A)") " "
     write(Flog,"(A,I0)") "STRUCTURE MATERIAU, pinceau ", nrank0
     do i=1,size(MattotP0)
        write(Flog,"(A,I0)") "----Indice ", i 
        write(Flog,"(A,I0)") "numM ", MattotP0(i)%numM
        if (MattotP0(i)%Nphase == 1) then
          write(Flog,"(A,A)") "Material name : ", trim(MattotP0(i)%Matename)
        end if
        if (MattotP0(i)%Interphase) then
           write(Flog,"(A)") " Materiau Interphase : present uniquement au sein des voxels composites "
        end if
        if (MattotP0(i)%Nphase > 1) then
           write(Flog,"(A)",advance='no') "  Materiau Composite constitue des materiaux : { "
           do j=1,size(MattotP0(i)%numM_composite)
              write(Flog,"(I0,A)",advance='no') MattotP0(i)%numM_composite(j)," "
           end do
           write(Flog,"(A)")"} "
           ntot_comp = ntot_comp + size(MattotP0(i)%pos,1)
        end if
        write(Flog,"(A,I0)") "taille du vecteur 'position' ", size(MattotP0(i)%pos,1)
        write(Flog,"(A,I0)") "taille du tableau 'zone' ", size(MattotP0(i)%zone,1)
        if (allocated(MattotP0(i)%Zones_interphase)) then
           write(Flog,"(A,I0,A)") "Materiau comprenant ",size(MattotP0(i)%Zones_interphase)," zones Interphases &
                         & presentes uniquement au sein des voxels composites "
           if (size(MattotP0(i)%Zones_interphase) > 10) then
               write(Flog,"(A)",advance='no') "  Liste des zones interphase : {"
               do j=1,5
                     write(Flog,"(I0,A)",advance='no') MattotP0(i)%Zones_interphase(j)
               end do
               write(Flog,"(A)",advance='no') " ... "
               do j=size(MattotP0(i)%Zones_interphase)-4,size(MattotP0(i)%Zones_interphase)
                     write(Flog,"(I0,A)",advance='no') MattotP0(i)%Zones_interphase(j)
               end do
               write(Flog,"(A)")"} "
            else
               write(Flog,"(A)",advance='no') "  Liste des zones interphase : {"
               do j=1,size(MattotP0(i)%Zones_interphase)
                     write(Flog,"(I0,A)",advance='no') MattotP0(i)%Zones_interphase(j)," "
               end do
               write(Flog,"(A)")"} "
            end if
        end if

        if (algo_param0%Mechanics) then
        if (MattotP0(i)%Nphase > 1) then
           write(Flog,"(2A)") "Modele d'homogeneisation : ", trim(MattotP0(i)%Lawname)
        else
           write(Flog,"(2A)") "libname ", trim(MattotP0(i)%Libname)
           write(Flog,"(2A)") "lawname ", trim(MattotP0(i)%Lawname)
        end if
        write(Flog,"(A,I0)") "NCoeff ", MattotP0(i)%NCoeff
        write(Flog,"(A,I0)") "Nvarint ", MattotP0(i)%Nvarint
        write(Flog,"(A,I0,A,I0)") "tailles du tableau 'Coeff' :", size(MattotP0(i)%Coeff,1)," - ",size(MattotP0(i)%Coeff,2)
        
        if (MattotP0(i)%Nphase == 1) then
        do k=1,MattotP0(i)%NCoeff
           write(Flog,"(A,I0,A,A,E15.8)") "Coeff(",k,",1), ", trim(MattotP0(i)%Coeffname(k))//":", MattotP0(i)%Coeff(k,1)
        end do
        
        do k=1,MattotP0(i)%Nvarint
           write(Flog,"(A,I0,A,A,E15.8)") "Varint0(",k,",1), ", trim(MattotP0(i)%Varintname(k))//":", MattotP0(i)%Varint0(k,1)
        end do
        end if 
           
        end if

        if (algo_param0%Diffusion) then
        write(Flog,"(2A)") "libnameK ", trim(MattotP0(i)%LibnameK)
        write(Flog,"(2A)") "lawnameK ", trim(MattotP0(i)%LawnameK)
!        write(Flog,"(2A)") "libnameA ", MattotP0(i)%LibnameA
!        write(Flog,"(2A)") "lawnameA ", MattotP0(i)%LawnameA
!        write(Flog,"(2A)") "libnameS ", MattotP0(i)%LibnameS
!        write(Flog,"(2A)") "lawnameS ", MattotP0(i)%LawnameS
        write(Flog,"(A,I0)") "NCoeffK ", MattotP0(i)%NCoeffK
!        write(Flog,"(A,I0)") "NCoeffA ", MattotP0(i)%NCoeffA
!        write(Flog,"(A,I0)") "NCoeffS ", MattotP0(i)%NCoeffS
        write(Flog,"(A,I0,A,I0)") "tailles du tableau 'CoeffK' :", size(MattotP0(i)%CoeffK,1)," - ",size(MattotP0(i)%CoeffK,2)
        write(Flog,"(A,E15.8)") "CoeffK(1,1) ", MattotP0(i)%CoeffK(1,1)
!        write(Flog,"(A,I0,A,I0)") "tailles du tableau 'CoeffA' :", size(MattotP0(i)%CoeffA,1)," - ",size(MattotP0(i)%CoeffA,2)
!        write(Flog,"(A,E15.8)") "CoeffA(1,1) ", MattotP0(i)%CoeffA(1,1)
!        write(Flog,"(A,I0,A,I0)") "tailles du tableau 'CoeffS' :", size(MattotP0(i)%CoeffS,1)," - ",size(MattotP0(i)%CoeffS,2)
!        write(Flog,"(A,E15.8)") "CoeffS(1,1) ", MattotP0(i)%CoeffS(1,1)
        end if
     end do
     write(Flog,"(A)") " "
     if (ntot_comp > 0) then
        write(Flog,'(A,I12,F12.4,A)') "   Nombre total/Fraction de voxels composites : ",ntot_comp,&
                                         100*real(ntot_comp)/(real(nx)*real(ny)*real(nz))," %"
        write(Flog,"(A)") " "
     end if

     if (algo_param0%Nloc) then
       write(Flog,"(A)") " "
       write(Flog,"(A)") "Specific NON LOCAL model"
       write(Flog,"(A)") "------------------------"
       do i=1,size(Nloc_models0)
           write(Flog,"(A,I0)") &
                 "Non local model ",i
           write(Flog,"(A,A,A)") &
                 "Implementation in subroutine : ",trim(Nloc_models0(i)%Modelname)
           write(Flog,"(A,I5,A)") &
                 "Number of input fields : ",Nloc_models0(i)%Nnloc," component(s)"
           write(Flog,"(A,I5,A)") &
                 "Number of output fields : ",Nloc_models0(i)%Ngnloc," component(s)"
           write(Flog,"(A)",advance='no') &
                 "Materials concerned by the model : "
           do j=1,size(Nloc_models0(i)%numM)
               write(Flog,"(I0,A)",advance='no') Nloc_models0(i)%numM(j)," "
           end do
           do j=1,size(Nloc_models0(i)%numM)
              write(Flog,"(/,A,I0)") " -- material ",Nloc_models0(i)%numM(j)

              write(Flog,"(A)",advance='no')"          int. var. indices used as input fields : "
              do k=1,size(Nloc_models0(i)%Ind_VarNloc(j,:))
                 write(Flog,"(I0,A)",advance='no') Nloc_models0(i)%Ind_VarNloc(j,k)," "
              end do

              write(Flog,"(/,A)",advance='no')"          int. var. indices used as output fields : "
              do k=1,size(Nloc_models0(i)%Ind_VarGNloc(j,:))
                 write(Flog,"(I0,A)",advance='no') Nloc_models0(i)%Ind_VarGNloc(j,k)," "
              end do

              write(Flog,"(/,A)",advance='no')"          coeff. indices used in the non-local model : "
              do k=1,size(Nloc_models0(i)%Ind_CoeffNloc(j,:))
                 write(Flog,"(I0,A)",advance='no') Nloc_models0(i)%Ind_CoeffNloc(j,k)," "
              end do

           end do
           write(Flog,"(/,A)") " "
       end do
     end if

  end if

end subroutine print_material_structure

!==================================================================================
!==================================================================================
!                         SUBROUTINE READ_MAT
!> Definition des materiaux
!!
!! A partir du fichier xml file_mat
!!  lit et affecte:
!!       - les proprietes du materiau de reference (lambda0, mu0)
!!       - nom de librairie et de loi pour chacun des materiaux (MattotP0(i)%Libname et MattotP0%Lawname )
!!       - les coefficients et variables internes ( MattotP0(i)%Coeff, MattotP0(i)%VarInt )
!!  et dans le cas de modeles non-locaux lit et affecte:
!!       - les indices de variables internes du materiau, pour chaque modele NL:
!!                     Nloc_models0(l)%Ind_VarNloc(m,:)
!!       - les indices de variables internes "NON-LOCALES" du materiau, pour chaque modele NL:
!!                     Nloc_models0(l)%Ind_VarGNloc(m,:)
!!
!!  \param[in] file_mat: (chaine de caracteres) nom du fichier xml renseignant les proprietes materiau
!!  \param[in] nmate: (entier) nombre de materiaux
!!
!! \param[out] lambdamu0: (reels) coeff. de Lame materiau de reference\n
!!             K0D      : (reels) prop. du materiau de reference \n
!!
!! Modifie aussi MattotP0
!!
!!
!! CORRECTION LE 4/4/2019 (Aldo Marano):
!!        la lecture parallele binaire (getbincoefffromfile via get_vect_xml)
!!        pouvait restee bloquer dans le cas de materiau absent sur un pinceau
!!        car get_vect_xml n'etait pas appelee
!!        -> necessite de traiter ce cas 
!!      VOIR AUSSI CORRECTIONS get_vect_xml et getbincoeffomfile et getCoeffFromFile
!!
!!
!! \TODO Prendre en compte les cas ou il y a plus de coefficients que de zones :
!!       on pourrait alors s'arreter de lire des que le nombre de zones est atteint
!==================================================================================
!subroutine read_mat(file_mat, lambdamu0, K0D,A0D, nmate)
subroutine read_mat(file_mat, lambdamu0, K0D, nmate)

  implicit none

  character(len=*), intent(in)                         :: file_mat
  real(mytype),dimension(2), intent(out)               :: lambdamu0
  real(mytype),dimension(algo_param0%nVarD),intent(out) :: K0D !,A0D
  integer, intent(in)                                  :: nmate

  type(node), pointer              :: fi, cur_node, cur_node2
  type(nodeList), pointer          :: material, coeff, intVar, coeffK
  integer(kind=4)                  :: i,j, m, nMat, ierror, ind
  character(len=1000)              :: err_msg,libname,lawname, matename,libnameK,lawnameK
  integer                          :: arraySize1, alloc_stat
  integer(kind=8)                  :: nbZone, arraySize2 
  integer,allocatable,dimension(:) :: tmp_array,tmp_arrayc
 
  integer(kind=8), dimension(nmate) :: zones
  integer(kind=8), dimension(nmate) :: zones_tmp

  integer(kind=8)      :: k,l
  logical                                                         :: material_presence
  integer(kind=4)                                                 :: mat_local_Index
  real(mytype),dimension(:,:),allocatable                         :: empty_array_real
  integer(kind=8),dimension(:,:),allocatable                      :: empty_array_zone

 
  allocate(tmp_array(1));allocate(tmp_arrayc(1));deallocate(tmp_array);deallocate(tmp_arrayc) ! avoid gcc-warning
  
  alloc_stat = 0  

  ! allocation des tableaux vides utilises comme entree de get_vect_xml sur 
  ! les pinceaux n'ayant rien a lire (materiau non present sur le pinceau)
  allocate(empty_array_real(0,0))
  allocate(empty_array_zone(0,2)) ! alloue a la taille (0,2) pour coller a la forme
                                 ! du tableau MattotP0()%Zone, ayant 2 colonnes 
                                 ! et autant de lignes que de zones, dans le cas
                                 ! ou le materiau ne possede aucune zone sur le pinceau
                                 
  ! Initialisation
  K0D=0_mytype
  !A0D=0._mytype
  lambdamu0=0._mytype

  ! zones_tmp: nombre maximal de zones pour chaque materiau et chaque proc
  zones_tmp=0
  do i=1,size(MattotP0)
    zones_tmp(MattotP0(i)%numM)=maxval(MattotP0(i)%zone(:,2))
  end do
  ! zones : nombre maximal de zones (dans tous les processus) pour chaque materiau
  call MPI_Allreduce(zones_tmp,zones,nmate, MPI_INTEGER8,MPI_MAX,amitex_mpi_comm,ierror)

  fi => parseFile(trim(file_mat))

  ! lecture des proprietes du materiau de reference
  if (algo_param0%Mechanics .eqv. .true.) then  
  call GetNodeList(fi, material, "Reference_Material", 1,1)  
  if(associated(material))then
    cur_node => item(material,0)
    call get_real_xml(cur_node,"Lambda0",lambdamu0(1),1,0)
    call get_real_xml(cur_node,"Mu0",lambdamu0(2),1,0)
  end if
  end if

  if (algo_param0%Diffusion .eqv. .true.) then  
  call GetNodeList(fi, material, "Reference_MaterialD", 1,1)  
  if(associated(material))then
    cur_node => item(material,0)
    call get_real_xml(cur_node,"K0",K0D(1),1,0)  !TODO extension a plusieurs variables
!    call get_real_xml(cur_node,"A0",A0D,1,0)
  end if
  end if

  !====================================================================== LECTURE MODELES NON-LOCAUX  
  call read_nloc_models(fi)
  
!  !====================================================================== BOUCLE SUR MODELES NON-LOCAUX  
!  call write_stdout0("before reading Non_local_modeling (read_param, material_mod)")
!  node_list => getElementsByTagName(fi,"Non_local_modeling")

!  if(getLength(node_list)>0)then
       
!    allocate(Nloc_models0(getLength(node_list)),stat=alloc_stat)
!    if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_param)",2,0)
       
!    do i=1,getLength(node_list)

!         cur_node => item(node_list,i-1)

!         !! Recuperation du numero d'identification du modele non local
!         string="NLocMod_num"
!         call get_int_xml(cur_node, string, ind0,1)

!         !! Recuperation du numero de modele
!         string="Modelname"
!         call get_str_xml(cur_node, trim(string), tmp_string,1)
!         Nloc_models0(ind0)%Modelname = tmp_string(1:200)

!         !! Recuperation du nombre de composantes du champ de variables
!         !! a deriver au sein du modele
!         string="Nnloc"
!         call get_int_xml(cur_node, string, Nloc_models0(ind0)%Nnloc,1)

!         !! Recuperation du nombre de composantes du champ de variables
!         !! non locales (= "gradients") au sein du modele
!         string="Ngnloc"
!         call get_int_xml(cur_node, string, Nloc_models0(ind0)%Ngnloc,1)

!         !! Recuperation du nombre de coefficients du modeles non local
!         !! a recuperer dans les coefficients materiaux
!         string="Ncoeff_nloc"
!         call get_int_xml(cur_node, string, Nloc_models0(ind0)%Ncoeff_nloc,1)

!         !! Recuperation des numeros des materiaux concernes par le modele
!         call getNodeList(cur_node,sub_list,"numM",-1,1)
!         if(getLength(sub_list)>0)then
!            if(getLength(sub_list)>1) then
!               !! Plusieurs noeuds numM dans le NL_model i
!               !! -> Erreur 
!               call amitex_abort(" More than one node 'numM' node &
!                    &'Non_local_modeling' (read_param)",2,0)
!            end if

!            sub_node => item(sub_list,0)
!            string="Nmat"
!            call get_int_xml(sub_node, string,n,1)

!            allocate(Nloc_models0(ind0)%numM(n),stat=alloc_stat)
!            if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_param)",2,0)
!            allocate(Nloc_models0(ind0)%Ind_VarNloc(n,Nloc_models0(i)%Nnloc),stat=alloc_stat)
!            if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_param)",2,0)
!            allocate(Nloc_models0(ind0)%Ind_VarGNloc(n,Nloc_models0(i)%Ngnloc),stat=alloc_stat)
!            if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_param)",2,0)
!            allocate(Nloc_models0(ind0)%Ind_CoeffNloc(n,Nloc_models0(i)%Ncoeff_nloc),stat=alloc_stat)
!            if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_param)",2,0)

!            allocate(tmp_array(n),stat=alloc_stat)
!            if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_param)",2,0)

!            call extractDataContent(sub_node,numM_list)
!            tmp_array=0
!            n=size(tmp_array)
!            call get_int_vect( numM_list, tmp_array, n, err)
!            Nloc_models0(ind0)%numM = tmp_array
!            deallocate(tmp_array)
!         else
!            call amitex_abort("noeud 'numM' manquant dans un noeud 'Non_local_modeling' (read_param)",2,0)
!         end if
!    end do
        
!  else
!         algo_param0%Nloc = .FALSE.
!  end if
!  call write_stdout0("after reading Non_local_modeling (read_param, material_mod)")
!  call mpi_barrier(amitex_mpi_comm,ierror)
!  !================================================================== FIN BOUCLE SUR MODELES NON-LOCAUX  

  !==========================================================================  BOUCLE SUR LES MATERIAUX

  ! identification des noeuds xml "Material" (compatibilite avec le nombre de matériaux issus des fichiers vtk) 
  material => getElementsByTagName(fi,"Material")
  if(nmate/=getLength(material))then
    write(err_msg,fmt="(A,I0,A,I0)") "Nombre de materiaux differents (read_mat) : vtk :",nmate,"; xml :",getLength(material)
    call amitex_abort(err_msg,2,0)
    call MPI_Barrier(amitex_mpi_comm,ierror)
  end if

  do i=0,nmate-1

    cur_node => item(material,i)
    ! numero de materiau du noeud i
    call get_int_xml(cur_node,"numM",nMat,2,0)

    if (algo_param0%Mechanics .eqv. .true.) then  !!------------------------------------lecture de la partie mecanique

       ! nom de la librairie (chemin complet ou "amitex") associée (pas en minuscule)
       call get_str_xml_nolowercase(cur_node,"Lib",libname,1,0)

       ! Default library if Lib="" in xml file
       if (libname .eq. "") then
       if (AMITEXenv%path_stat .ne. 0) then
          call amitex_abort(&
          "Undefined environment variable AMITEX_FFTP",2,0) ! ERROR -> STOP
       else
          libname = AMITEXenv%path // "/libAmitex/src/materiaux/libUmatAmitex.so"
       end if
       end if

       ! nom de la loi associée 
       call get_str_xml(cur_node,"Law",lawname,1,0)
       
       ! material name
       matename=""
       call get_str_xml(cur_node,"Name",matename,-2)

       ! listes des noeuds de coefficients et variables internes
       coeff => getElementsByTagName(cur_node,"Coeff")
       intVar=> getElementsByTagName(cur_node,"IntVar")

        material_presence = .false. 
        mat_local_Index = -1
        do j=1,size(MattotP0)
            ! recherche du materiau dans MattotP0 pour savoir si il est sur le pinceau
            if(nMat== MattotP0(j)%numM)then
                material_presence = .true. 
                mat_local_Index = j
            end if
        end do

            ! bibliotheque .so et nom de loi   
            ! affectés si le materiau est sur le pinceau
        if (material_presence) MattotP0(mat_local_Index)%Libname = libname(1:200)
        if (material_presence) MattotP0(mat_local_Index)%Lawname = lawname(1:200)
        if (material_presence) MattotP0(mat_local_Index)%Matename = matename(1:100)

!       ! numero de loi
!       MattotP0(j)%NumeLoi = adjustr(NumeLoi)
        nbZone=zones(nMat)  ! nombre maximal de zones du materiau

        !------- lecture des coefficients
        if (material_presence) then
        ! la materiau est present sur le pinceaux : mise en forme des donnes necessaires
        ! a la recuperation des coefficients

            ! mise en forme des entrees et appel a get_vect_xml dans le cas ou le materiau est present
            arraySize1= getLength(coeff)           ! nombre de coefficients
            MattotP0(mat_local_Index)%NCoeff=arraySize1
            arraySize2= size(MattotP0(mat_local_Index)%zone(:,1),kind=8) ! nombre de zones du materiau sur le processus
            ! allocation du tableau coeff
            allocate(MattotP0(mat_local_Index)%coeff( arraySize1, arraySize2),stat=alloc_stat )
            if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant MattotP0%coeff (read_mat)",2)
            ! lecture des valeurs des coefficients
            call get_vect_xml(nMat,coeff,arraySize1,arraySize2,MattotP0(mat_local_Index)%coeff,nbZone,arraySize2,&
                              MattotP0(mat_local_Index)%zone,"Coeff",MattotP0(mat_local_Index)%pos)
            ! allocation du tableau Coeffname
            allocate(MattotP0(mat_local_Index)%Coeffname(arraySize1),stat=alloc_stat )
            if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant MattotP0%Coeffname (read_mat)",2)

            ! Read Coeff names :
            ! REMARK : below, no check on ind, already done in get_vect_xml
            do m=0,arraySize1-1
               cur_node2 => item(coeff,m)
               ind = -1
               matename=""
               call get_str_xml(cur_node2,"Name",matename,-2)
               call get_int_xml(cur_node2,"Index",ind,1)
               MattotP0(mat_local_Index)%Coeffname(ind) = matename(1:100)
            end do
                        
           do l=1,arraySize2
             do k=1,arraySize1
               if(MattotP0(mat_local_Index)%coeff(k,l)/=MattotP0(mat_local_Index)%coeff(k,l))then
                 write(err_msg,fmt="(3(A,I0),A)") "Coefficient:",k, " zone:",MattotP0(j)%zone(l,2),&
                   " materiau:",nMat," sans valeur (read_mat)"
                 call amitex_abort(err_msg,1)
               end if
             end do
           end do
        else
        ! materiau non present sur le pinceau : on envoie des tableaux vides à get_vect_xml
            arraysize1 = 0
            arraysize2 = 0
            ! appel a get_vect_xml pour assurer que les pinceaux ne possedant pas le materiau 
            ! soit associés a la lecture parallèle des fichiers binaires
            !call get_vect_xml(nMat,coeff,arraySize1,arraysize2,empty_array_real,nbZone,&
            !                  arraySize2,empty_array_zone,"Coeff",MattotP0(mat_local_Index)%pos) 
            !attention :  ici mat_local_indice=-1 -> on passe un tableau INT64 de taille nulle!!!
            call get_vect_xml(nMat,coeff,arraySize1,arraysize2,empty_array_real,nbZone,&
                              arraySize2,empty_array_zone,"Coeff",[integer(kind=INT64) ::])

        end if
                   
        !------- lecture des variables internes
        arraySize1=getLength(intVar)        ! nombre de variables internes
        if (material_presence) MattotP0(mat_local_Index)%NVarInt=arraySize1
        if(arraySize1>0)then
        ! variables internes a lire 
            if (material_presence) then
            ! lecture des variables internes dans le cas ou le materiau est sur le pinceau
              arraySize2=size(MattotP0(mat_local_Index)%pos,kind=8)   ! nombre de voxels du materiau dans le processus
              l=size(MattotP0(mat_local_Index)%zone(:,1),kind=8)      ! nombre de zones 
                                                       ! ATTENTION : pas de bug ici en presence de zones interphases
                                                       ! elles sont comprises dans Zone(:,1) et doivent etre comptees
                                                       ! pour pouvoir les lire
              allocate(MattotP0(mat_local_Index)%VarInt(arraySize1,arraySize2),stat=alloc_stat)
              if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant MattotP0%VarInt (read_mat)",2)
              allocate(MattotP0(mat_local_Index)%VarInt0(arraySize1,arraySize2),stat=alloc_stat)
              if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant MattotP0%VarInt0 (read_mat)",2)
              MattotP0(mat_local_Index)%VarInt0=-1
              ! lecture des variables internes
              call get_vect_xml(nMat,intVar,arraySize1,arraySize2,MattotP0(mat_local_Index)%VarInt0,nbZone,&
                                l,MattotP0(mat_local_Index)%zone,"IntVar",MattotP0(mat_local_Index)%pos)


              ! allocation du tableau Varintname
              allocate(MattotP0(mat_local_Index)%Varintname(arraySize1),stat=alloc_stat )
              if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant MattotP0%Varintname (read_mat)",2)

              ! Read Varint names :
              ! REMARK : below, no check on ind, already done in get_vect_xml
              do m=0,arraySize1-1
                 cur_node2 => item(intVar,m)
                 ind = -1
                 matename=""
                 call get_str_xml(cur_node2,"Name",matename,-2)
                 call get_int_xml(cur_node2,"Index",ind,1)
                 MattotP0(mat_local_Index)%Varintname(ind) = matename(1:100)
              end do

              do l=1,size(MattotP0(mat_local_Index)%zone(:,1)) 
                ! ATTENTION
                ! pas de bug ici vis a vis de la presence de zones interphases. Les lignes correspondantes de Zone(:,1)
                ! contiennent les numeros de colonne de %VarInt où stocker les variables internes de ces zones --> OK 
                do k=1,arraySize1
                  if( MattotP0(mat_local_Index)%VarInt0(k,MattotP0(mat_local_Index)%zone(l,1)) /= &
                      MattotP0(mat_local_Index)%VarInt0(k,MattotP0(mat_local_Index)%zone(l,1)))then
                    write(err_msg,fmt="(3(A,I0),A)") "Variables internes:",k, " zone:",MattotP0(mat_local_Index)%zone(l,2),&
                    " materiau:",nMat," sans valeur (read_mat)"
                    call amitex_abort(err_msg,1)
                  end if
                end do
              end do

            else
            ! materiau non present sur le pinceau : on envoie des tableaux vides à get_vect_xml
                arraysize1 = 0
                arraysize2 = 0
                ! appel a get_vect_xml pour assurer que les pinceaux ne possedant pas le materiau 
                ! soit associés a la lecture parallèle des fichiers binaires
                call get_vect_xml(nMat,intVar,arraySize1,arraysize2,empty_array_real,nbZone,arraySize2,&
                                  empty_array_zone,"IntVar",[integer(kind=INT64) ::])    
            end if
        else
        ! pas de variables internes a lire : allocation de tableaux vides
          if (material_presence) then
          allocate(MattotP0(mat_local_Index)%VarInt(0,0),stat=alloc_stat)
          if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant MattotP0%VarInt(0,0) (read_mat)",2)
          allocate(MattotP0(mat_local_Index)%VarInt0(0,0),stat=alloc_stat)
          if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant MattotP0%VarInt0(0,0) (read_mat)",2)
          end if
        end if  
          
          !------------------------------------------------------------------------------------------
          ! MECANIQUE NON LOCALE
          ! lecture des indices des variables internes a deriver et stockage dans 
          ! algo_param0 (plus simple pour l'appel a un modele non local que de passer
          ! par MattotP0 --> informations requises dans chaque pinceau or mattot est local

          if (algo_param0%Nloc) then
              cur_node => item(material,i)
              call read_nloc_indices(cur_node,nMat,intVar)
              
!             !---------------------------------------------------------------IndVarNloc / IndCoeffNloc
!             ! Recuperation des noeuds IndexVarNloc
!             IndVarNloc => getElementsByTagName(cur_node,"IndexVarNloc") 
!             IndCoeffNloc => getElementsByTagName(cur_node,"IndexCoeffNloc") 

!             if (getLength(IndVarNloc) .ne. getLength(IndCoeffNloc)) then
!             if(alloc_stat /=0) call amitex_abort(&
!               "The number of nodes IndVarNloc and IndCoeffNloc are different (read_mat)",2,0)
!             end if 

!             ! boucle sur ces noeuds (boucle sur les modeles non-locaux concernant ce materiau) 
!             !     Rq : tests de verif des entrees realises ulterieurement
!             do k=1,getLength(IndVarNloc)

!                sub_node => item(IndVarNloc,int(k-1,kind=4)) ! recuperation du kieme noeud IndVarNloc
!                sub_nodec => item(IndCoeffNloc,int(k-1,kind=4)) ! recuperation du kieme noeud IndCoeffNloc
!                ! recuperation du numero du modele concerne
!                call get_int_xml(sub_node,"NLocMod_num",ind,2,0)
!                call get_int_xml(sub_nodec,"NLocMod_num",ind0,2,0)
!                test = .false.
!                      ! modele correspondant trouve --> recherche du rang du materiau dans le modele
!                      do m=1,size(Nloc_models0(ind)%numM)
!                      if (Nloc_models0(ind)%numM(m) == nMat) then
!                            ! allocation de tmp_array 
!                            allocate(tmp_array(Nloc_models0(ind)%Nnloc),stat=alloc_stat)
!                            if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant tmp_array (read_mat)",2,0)
!                            tmp_array = 0
!                            allocate(tmp_arrayc(Nloc_models0(ind)%Ncoeff_nloc),stat=alloc_stat)
!                            if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant tmp_arrayc (read_mat)",2,0)
!                            tmp_arrayc = 0
!                            test = .true.
                            
!                            ! Recuperation des indices dans le tableau de variable interne                
!                            call extractDataContent(sub_node,IndNloc_list)
!                            call extractDataContent(sub_nodec,IndNloc_listc)
!                            Nloc_models0(ind)%Ind_VarNloc(m,:) = 0
!                            Nloc_models0(ind)%Ind_CoeffNloc(m,:) = 0
!                            call get_int_vect(IndNloc_list,tmp_array,Nloc_models0(ind)%Nnloc , err)
!                            call get_int_vect(IndNloc_listc,tmp_arrayc,Nloc_models0(ind)%NCoeff_nloc , err)
!                            Nloc_models0(ind)%Ind_VarNloc(m,:) = tmp_array
!                            Nloc_models0(ind)%Ind_CoeffNloc(m,:) = tmp_arrayc
!                            ! test ! tous les elements du tableau initialises ?
!                            if (any(Nloc_models0(ind)%Ind_VarNloc(m,:)==0)) then
!                               !! Plusieurs noeuds numM dans le NL_model i
!                               !! On ecrit un warning 
!                               write(err_msg,fmt="(A,I0,A)") "Noeud 'IndexVarNloc' incomplet pour le materiau ",&
!                                    nMat,", nombre d'indices inferieur a",Nloc_models0(ind)%Nnloc,&
!                                    " pour le modele non local ", ind," (read_mat)"
!                               call amitex_abort(err_msg,1)
!                            end if
!                            if (any(Nloc_models0(ind)%Ind_CoeffNloc(m,:)==0)) then
!                               !! Plusieurs noeuds numM dans le NL_model i
!                               !! On ecrit un warning 
!                               write(err_msg,fmt="(A,I0,A)") "Noeud 'IndexCoeffNloc' incomplet pour le materiau ",&
!                                    nMat,", nombre d'indices inferieur a",Nloc_models0(ind)%NCoeff_nloc,&
!                                    " pour le modele non local ", ind," (read_mat)"
!                               call amitex_abort(err_msg,1)
!                            end if
!                            deallocate(tmp_array)
!                            deallocate(tmp_arrayc)
!                      end if
!                      end do

!                if (.not.test) then
!                   write(err_msg,fmt="(2(A,I0),A)") "Noeud IndVarNloc refere au modele non local",ind,&
!                        "pour le  materiau:",nMat,", mais le modele n'est pas associe a ce materiau &
!                             & dans le fichier param_algo.xml (read_mat)"
!                   call amitex_abort(err_msg,1)
!                end if

!                if (allocated(tmp_array)) deallocate(tmp_array)
!             end do  ! Fin lecture des indices IndVarNloc

!             !---------------------------------------------------------------IndVarGNloc
!             !             les indices GNloc sont les derniers indices du vecteur VarInt
!             !                             dans l'odre des numeros de modeles croissants

!             ! On place le compteur en fin de vecteur Varint
!             Compteur = getLength(intVar)
!             ! boucle sur les modeles par numero decroissant 
!             do k = size(Nloc_models0),1,-1 
!                 ! Affectation de Nloc_models0(k)%Ind_VarGNloc pour le materiau concerne
!                 !  (index=nMat du noeud material dans mat.xml)
!                 do m=1,size(Nloc_models0(k)%numM)
!                 if (Nloc_models0(k)%numM(m) == nMat) then
!                    ! calcul et affectation des indices
!                    Nloc_models0(k)%Ind_VarGNloc(m,:) = &
!                              (/(p, p=(Compteur- Nloc_models0(k)%NGnloc +1),Compteur) /)
!                    ! decalage du compteur
!                    Compteur = Compteur - Nloc_models0(k)%NGnloc
!                 end if
!                 end do
!             end do ! fin calcul des indices IndVarGNloc
          
          end if ! FIN MECANIQUE NON LOCALE

    end if !! fin de la lecture de la partie "mecanique"

    if (algo_param0%Diffusion .eqv. .true.) then  !!------------------------------------lecture de la partie Diffusion

       ! nom de la librairie (chemin complet ou "amitex") associée (pas en minuscule)
       call get_str_xml_nolowercase(cur_node,"LibK",libnameK,1,0)
       !call get_str_xml_nolowercase(cur_node,"LibA",libnameA,1,0)
       !call get_str_xml_nolowercase(cur_node,"LibS",libnameS,1,0)

       ! Default library if LibK="" in xml file
       if (libnameK .eq. "") then
       if (AMITEXenv%path_stat .ne. 0) then
          call amitex_abort(&
          "Undefined environment variable AMITEX_FFTP",2,0) ! ERROR -> STOP
       else
          libnameK = AMITEXenv%path // "/libAmitex/src/materiauxK/libUmatAmitexK.so"
       end if
       end if
       
       ! nom de la loi associée 
       call get_str_xml(cur_node,"LawK",lawnameK,1,0)
       !call get_str_xml(cur_node,"LawA",lawnameA,1,0)
       !call get_str_xml(cur_node,"LawS",lawnameS,1,0)

       ! listes des noeuds de coefficients 
       coeffK => getElementsByTagName(cur_node,"CoeffK")
       !coeffA => getElementsByTagName(cur_node,"CoeffA")
       !coeffS => getElementsByTagName(cur_node,"CoeffS")


        material_presence = .false. 
        mat_local_Index = -1
        do j=1,size(MattotP0)
            ! recherche du materiau dans MattotP0 pour savoir si il est sur le pinceau
            if(nMat== MattotP0(j)%numM)then
                material_presence = .true. 
                mat_local_Index = j
            end if
        end do

       ! bibliotheque .so et nom de loi   
       if (material_presence) MattotP0(mat_local_Index)%LibnameK = libnameK(1:200)
       if (material_presence) MattotP0(mat_local_Index)%LawnameK = lawnameK(1:200)

       nbZone=zones(nMat)  ! nombre maximal de zones du materiau

       ! lecture des coefficients K
       if (material_presence) then
           arraySize1= getLength(coeffK)           ! nombre de coefficients K
           MattotP0(mat_local_Index)%NCoeffK=arraySize1
           arraySize2= size(MattotP0(mat_local_Index)%zone(:,1),kind=8) ! nombre de zones du materiau sur le processus
           allocate(MattotP0(mat_local_Index)%coeffK( arraySize1, arraySize2),stat=alloc_stat )
           if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant MattotP0%coeffK (read_mat)",2)
           ! lecture des valeurs des coefficients
           call get_vect_xml(nMat,coeffK,arraySize1,arraySize2,MattotP0(mat_local_Index)%coeffK,nbZone,arraySize2,&
                             MattotP0(mat_local_Index)%zone,"CoeffK",MattotP0(mat_local_Index)%pos)

           do l=1,arraySize2
             do k=1,arraySize1
               if(MattotP0(mat_local_Index)%coeffK(k,l)/=MattotP0(mat_local_Index)%coeffK(k,l))then
                 write(err_msg,fmt="(3(A,I0),A)") "CoefficientK:",k, " zone:",MattotP0(mat_local_Index)%zone(l,2),&
                   " materiau:",nMat," sans valeur (read_mat)"
                 call amitex_abort(err_msg,1)
               end if
             end do
           end do
        else 
        ! materiau non present sur le pinceau : on envoie des tableaux vides à get_vect_xml
            arraysize1 = 0
            arraysize2 = 0
            ! appel a get_vect_xml pour assurer que les pinceaux ne possedant pas le materiau 
            ! soit associés a la lecture parallèle des fichiers binaires
            call get_vect_xml(nMat,coeffK,arraySize1,arraysize2,empty_array_real,nbZone,&
                              arraySize2,empty_array_zone,"CoeffK",MattotP0(mat_local_Index)%pos)        
        end if
        
           !TODO 
           ! lecture des coefficients A (copie-colle lecture coeff K)
           ! lecture des coefficients S (copie-colle lecture coeff K)

    end if !! fin de la lecture de la partie "Diffusion"


  end do
  !====================================================================== FIN BOUCLE SUR LES MATERIAUX  

  call destroy(fi)

  call initParam_umat()

end subroutine read_mat

!==================================================================================
!                         SUBROUTINE READ_NLOC_MODELS
!
!> Lecture des modeles non-locaux (appele par read_mat et read_mat_composite)
!!
!!    * Lecture du noeud "Non_local_modeling" dans le fichier mate.xml
!!    * Allocation et affectation associees au sein de nloc_modesl0 
!!  
!!
!!  \param[in] fi (fi => parseFile(trim(file_mat))
!!
!! 03/2022 : creation, on isole cette partie commune a read_mat et read_mat_composite
!!           en vue de l'utiliser dans read_mat_composite
!==================================================================================
subroutine read_nloc_models(fi)  

  implicit none
  type(node), pointer, intent(in) :: fi
  
  type(node), pointer              :: cur_node, sub_node
  type(nodeList), pointer          :: node_list, sub_list     ! material, coeff, intVar, coeffK, IndVarNloc, IndCoeffNloc
  integer(kind=4)                  :: ind0, n     !  i,j, nMat, ierror,ind,

  integer               :: alloc_stat,ierror,err
  integer               :: i
  character(len=1000)   :: string, tmp_string
  character(len=100000) :: numM_list
  integer,allocatable,dimension(:) :: tmp_array
  
  !====================================================================== BOUCLE SUR MODELES NON-LOCAUX  
  call write_stdout0("before reading Non_local_modeling (read_param, material_mod)")
  node_list => getElementsByTagName(fi,"Non_local_modeling")

  if(getLength(node_list)>0)then
       
    allocate(Nloc_models0(getLength(node_list)),stat=alloc_stat)
    if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_param)",2,0)
       
    do i=1,getLength(node_list)

         cur_node => item(node_list,i-1)

         !! Recuperation du numero d'identification du modele non local
         string="NLocMod_num"
         call get_int_xml(cur_node, string, ind0,1)

         !! Recuperation du numero de modele
         string="Modelname"
         call get_str_xml(cur_node, trim(string), tmp_string,1)
         Nloc_models0(ind0)%Modelname = tmp_string(1:200)

         !! Recuperation du nombre de composantes du champ de variables
         !! a deriver au sein du modele
         string="Nnloc"
         call get_int_xml(cur_node, string, Nloc_models0(ind0)%Nnloc,1)

         !! Recuperation du nombre de composantes du champ de variables
         !! non locales (= "gradients") au sein du modele
         string="Ngnloc"
         call get_int_xml(cur_node, string, Nloc_models0(ind0)%Ngnloc,1)

         !! Recuperation du nombre de coefficients du modeles non local
         !! a recuperer dans les coefficients materiaux
         string="Ncoeff_nloc"
         call get_int_xml(cur_node, string, Nloc_models0(ind0)%Ncoeff_nloc,1)

         !! Recuperation des numeros des materiaux concernes par le modele
         call getNodeList(cur_node,sub_list,"numM",-1,1)
         if(getLength(sub_list)>0)then
            if(getLength(sub_list)>1) then
               !! Plusieurs noeuds numM dans le NL_model i
               !! -> Erreur 
               call amitex_abort(" More than one node 'numM' node &
                    &'Non_local_modeling' (read_param)",2,0)
            end if

            sub_node => item(sub_list,0)
            string="Nmat"
            call get_int_xml(sub_node, string,n,1)

            allocate(Nloc_models0(ind0)%numM(n),stat=alloc_stat)
            if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_param)",2,0)
            allocate(Nloc_models0(ind0)%Ind_VarNloc(n,Nloc_models0(i)%Nnloc),stat=alloc_stat)
            if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_param)",2,0)
            allocate(Nloc_models0(ind0)%Ind_VarGNloc(n,Nloc_models0(i)%Ngnloc),stat=alloc_stat)
            if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_param)",2,0)
            allocate(Nloc_models0(ind0)%Ind_CoeffNloc(n,Nloc_models0(i)%Ncoeff_nloc),stat=alloc_stat)
            if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_param)",2,0)

            allocate(tmp_array(n),stat=alloc_stat)
            if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_param)",2,0)

            call extractDataContent(sub_node,numM_list)
            tmp_array=0
            n=size(tmp_array)
            call get_int_vect( numM_list, tmp_array, n, err)
            Nloc_models0(ind0)%numM = tmp_array
            deallocate(tmp_array)
         else
            call amitex_abort("noeud 'numM' manquant dans un noeud 'Non_local_modeling' (read_param)",2,0)
         end if
    end do
        
  else
         algo_param0%Nloc = .FALSE.
  end if
  call write_stdout0("after reading Non_local_modeling (read_param, material_mod)")
  call mpi_barrier(amitex_mpi_comm,ierror)
  !================================================================== FIN BOUCLE SUR MODELES NON-LOCAUX  
end subroutine read_nloc_models  

!==================================================================================
!                         SUBROUTINE READ_NLOC_INDICES
!
!> Lecture des indices de coeff et varint pour les modeles non-locaux 
!>         (appele par read_mat et read_mat_composite)
!!
!!  Pour chaque materiau du fichier mate.xml :
!!            lit et affecte IndexVarNloc, IndexCoeffNloc
!!            affecte IndevVarGNloc
!!
!!  \param[in] cur_node  type(node) correspondant au materiau considere
!!  \param[in] nMat      indice nMat du materiau
!!  \param[in] intVar    type(nodeList) correspondant a la liste de variables internes 
!!                       du materiau nMat
!!
!! 03/2022 : creation, on isole cette partie commune a read_mat et read_mat_composite
!!           en vue de l'utiliser dans read_mat_composite
!==================================================================================
  
subroutine read_nloc_indices(cur_node,nMat,intVar)
  
  implicit none
  type(node), pointer,intent(in)        :: cur_node
  type(nodeList), pointer, intent(in)   :: intVar
  integer, intent(in)                   :: nMat
  
  type(node), pointer                   :: sub_node,sub_nodec
  type(nodeList), pointer               :: IndVarNloc, IndCoeffNloc  
  
  integer                               :: alloc_stat, err
  integer                               :: k, m, p, ind, ind0, compteur
  logical                               :: test
  integer,allocatable,dimension(:)      :: tmp_array,tmp_arrayc
  character(len=1000)                   :: err_msg 
  character(len=100000)                 :: IndNloc_list,IndNloc_listc
  
    alloc_stat = 0
    !---------------------------------------------------------------IndVarNloc / IndCoeffNloc
    ! Recuperation des noeuds IndexVarNloc
    IndVarNloc => getElementsByTagName(cur_node,"IndexVarNloc") 
    IndCoeffNloc => getElementsByTagName(cur_node,"IndexCoeffNloc") 

    if (getLength(IndVarNloc) .ne. getLength(IndCoeffNloc)) then
    if(alloc_stat /=0) call amitex_abort(&
      "The number of nodes IndVarNloc and IndCoeffNloc are different (read_mat)",2,0)
    end if 

    ! boucle sur ces noeuds (boucle sur les modeles non-locaux concernant ce materiau) 
    !     Rq : tests de verif des entrees realises ulterieurement
    do k=1,getLength(IndVarNloc)

       sub_node => item(IndVarNloc,int(k-1,kind=4)) ! recuperation du kieme noeud IndVarNloc
       sub_nodec => item(IndCoeffNloc,int(k-1,kind=4)) ! recuperation du kieme noeud IndCoeffNloc
       ! recuperation du numero du modele concerne
       call get_int_xml(sub_node,"NLocMod_num",ind,2,0)
       call get_int_xml(sub_nodec,"NLocMod_num",ind0,2,0)
       test = .false.
    ! modele correspondant trouve --> recherche du rang du materiau dans le modele
    do m=1,size(Nloc_models0(ind)%numM)
    if (Nloc_models0(ind)%numM(m) == nMat) then
          ! allocation de tmp_array 
          allocate(tmp_array(Nloc_models0(ind)%Nnloc),stat=alloc_stat)
          if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant tmp_array (read_mat)",2,0)
          tmp_array = 0
          allocate(tmp_arrayc(Nloc_models0(ind)%Ncoeff_nloc),stat=alloc_stat)
          if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant tmp_arrayc (read_mat)",2,0)
          tmp_arrayc = 0
          test = .true.
          
          ! Recuperation des indices dans le tableau de variable interne       
          call extractDataContent(sub_node,IndNloc_list)
          call extractDataContent(sub_nodec,IndNloc_listc)
          Nloc_models0(ind)%Ind_VarNloc(m,:) = 0
          Nloc_models0(ind)%Ind_CoeffNloc(m,:) = 0
          call get_int_vect(IndNloc_list,tmp_array,Nloc_models0(ind)%Nnloc , err)
          call get_int_vect(IndNloc_listc,tmp_arrayc,Nloc_models0(ind)%NCoeff_nloc , err)
          Nloc_models0(ind)%Ind_VarNloc(m,:) = tmp_array
          Nloc_models0(ind)%Ind_CoeffNloc(m,:) = tmp_arrayc
          ! test ! tous les elements du tableau initialises ?
          if (any(Nloc_models0(ind)%Ind_VarNloc(m,:)==0)) then
             !! Plusieurs noeuds numM dans le NL_model i
             !! On ecrit un warning 
             write(err_msg,fmt="(A,I0,A)") "Noeud 'IndexVarNloc' incomplet pour le materiau ",&
         nMat,", nombre d'indices inferieur a",Nloc_models0(ind)%Nnloc,&
         " pour le modele non local ", ind," (read_mat)"
             call amitex_abort(err_msg,1)
          end if
          if (any(Nloc_models0(ind)%Ind_CoeffNloc(m,:)==0)) then
             !! Plusieurs noeuds numM dans le NL_model i
             !! On ecrit un warning 
             write(err_msg,fmt="(A,I0,A)") "Noeud 'IndexCoeffNloc' incomplet pour le materiau ",&
         nMat,", nombre d'indices inferieur a",Nloc_models0(ind)%NCoeff_nloc,&
         " pour le modele non local ", ind," (read_mat)"
             call amitex_abort(err_msg,1)
          end if
          deallocate(tmp_array)
          deallocate(tmp_arrayc)
    end if
    end do

       if (.not.test) then
          write(err_msg,fmt="(2(A,I0),A)") "Noeud IndVarNloc refere au modele non local",ind,&
      "pour le  materiau:",nMat,", mais le modele n'est pas associe a ce materiau &
           & dans le fichier param_algo.xml (read_mat)"
          call amitex_abort(err_msg,1)
       end if

       if (allocated(tmp_array)) deallocate(tmp_array)
    end do  ! Fin lecture des indices IndVarNloc

    !---------------------------------------------------------------IndVarGNloc
    !    les indices GNloc sont les derniers indices du vecteur VarInt
    !           dans l'odre des numeros de modeles croissants

    ! On place le compteur en fin de vecteur Varint
    Compteur = getLength(intVar)
    ! boucle sur les modeles par numero decroissant 
    do k = size(Nloc_models0),1,-1 
        ! Affectation de Nloc_models0(k)%Ind_VarGNloc pour le materiau concerne
        !  (index=nMat du noeud material dans mat.xml)
        do m=1,size(Nloc_models0(k)%numM)
        if (Nloc_models0(k)%numM(m) == nMat) then
           ! calcul et affectation des indices
           Nloc_models0(k)%Ind_VarGNloc(m,:) = &
            (/(p, p=(Compteur- Nloc_models0(k)%NGnloc +1),Compteur) /)
           ! decalage du compteur
           Compteur = Compteur - Nloc_models0(k)%NGnloc
        end if
        end do
    end do ! fin calcul des indices IndVarGNloc
          
end subroutine read_nloc_indices

!==============================================================================
!
!                      SUBROUTINE BEHAVIOR
!
!> Calcul du comportement du materiau
!!
!!  Appel a la fonction umat pour chacun des voxels de la cellule
!!       l'ordre de parcours respecte la structure MATERIAL.
!!
!!
!!  WARNING : SIG=SIG0 en debut de procedure
!!              Utilise mattotp%varint0 en entree de umat
!!  
!!  \param[out]    SIG: champ des tenseurs de contraintes (a l'instant t)
!!  \param[in]     SIG0: champ des tenseurs de contraintes (au debut du chargement)
!!  \param[in]     DEF: champ des tenseurs de deformations (a l'instant t)
!!  \param[in]     DEF0: champ des tenseurs de deformations (au debut du chargement)
!!  \param[in]     TEMPRED: temperature et parametres externes de la loi de comportement
!!  \param[in]     DTEMPRED: increments de temperature et parametres externes de la loi de comportement
!!  \param[in]     t_load: temps aux instants t, t-1 et t-2
!!  \param[in]     KSTEP indice de chargement 
!!  \param[in]     Def_nsym : gradient du deplacement si pris en compte dans un calcul
!!                            HPP (optionnel) (a l'instant t)
!!  \param[in]     Def_nsym0 : gradient du deplacement si pris en compte dans un calcul
!!                            HPP (optionnel) (au debut du chargement)
!!  \param[in]     ext_law (optionnel), chaine de caracteres 'extension' du nom de "Lawname"
!!                           pour utilisation dans des modeles 'user'
!!  \param[in]     ext_lib (optionnel), idem pour le nom de la librairie
!!                           (Pas encore implemente : TODO si besoin)
!!  \param[in]     use_sig0 (optional), by default uses SIG=SIG0
!!                      if true (default), uses sig=sig0 as umat input, 
!!                      else uses sig (sig0 remains unmodified)
!!  \param[in]     use_varint0 (optional), 
!!                      if true (default), uses mattotP%varint0 as umat input, 
!!                      else uses mattotP%varint
!!  \param[out]    nSub_tot : nombre de recours à la subdivision du pas de temps dans PilotageLaminate
!!  \param[out]    nIncr_tot : nombre total de sous pas de temps utilisés dans PilotageLaminate
!!  \param[out]    nIt_tot : nombre total d'itérations dans umatLaminate
!!  \param[out]    nVoxComp_tot : nombre total de voxels composites "laminate"
!!
!! Modifie egalement les variables internes dans MattotP
!!
!! UMAT 'extension'
!! ----------------
!!    ext_law (+ext_lib envisagee) permet de faire appel a behavior avec des lois
!!           lawname_"ext_law" pour les metriaux homogenes et composites
!!
!!    
!!    UTILISATION : 
!!           modeles non-locaux pour mettre a jour une variable interne Nloc
!!
!!    ATTENTION : les variables internes "mécaniques" ne doivent pas etre modifiees
!!                la contrainte (STRESS) ne doit pas etre modifiee
!!
!!    ATTENTION :  cet appel suppose que les donnees d'entrees 
!!                 (Varint0, sig et def+ddef) sont equilibrees
!!
!!         CECI EST TRES IMPORTANT POUR LES VOXELS COMPOSITES 
!!         car le chargement de chacune des phases est stocke dans le vecteur 
!!         variable interne 0 (voir details dans entetes de umatVoigt, umatReuss...)
!!         il est donc PRIMORDIAL QUE varint0 soit actualise avant cet appel
!!   
!==============================================================================
subroutine behavior(SIG, SIG0, DEF, DEF0, TEMPRED,DTEMPRED, t_load, KSTEP0,&
                    nSub_tot,nIncr_tot,nIt_tot,&
                    nVoxComp_tot_s,Def_nsym,Def_nsym0,ext_law,ext_lib,ext_VC,use_varint0,use_sig0)

  implicit none

  real(mytype), dimension(amx_size(1)*amx_size(2)*amx_size(3),6)  :: SIG, SIG0 
  real(mytype), dimension(amx_size(1)*amx_size(2)*amx_size(3),algo_param%nTensDef) :: DEF, DEF0
  real(mytype), dimension(amx_size(1)*amx_size(2)*amx_size(3),9),optional          :: Def_nsym,Def_nsym0
  real(mytype),dimension(-2:0), intent(in)               :: t_load
  real(mytype),dimension(:), intent(in)                  :: TEMPRED,DTEMPRED
  integer(kind=8),intent(in)                             :: KSTEP0
  integer(kind=8),intent(out)                            :: nSub_tot,nIncr_tot,nIt_tot
  integer(kind=8),intent(out),optional                   :: nVoxComp_tot_s
  character(len=*),intent(in),optional                   :: ext_law,ext_lib
  logical,intent(in),optional                            :: ext_VC      !if true, extended behavior is applied to VC (default) 
  logical,intent(in),optional                            :: use_varint0 !if true (default), uses varint0 as umat input, else uses varint 
  logical,intent(in),optional                            :: use_sig0    !if true (default), uses sig0 as umat input, else uses sig 


  type(c_ptr),allocatable,dimension(:) :: umatptr    ! pointeur de fonction
  integer(kind=INT64)                  :: KINC
  integer                              :: i,h
#ifdef OPENMP
  integer                              :: nbthread
#endif
  integer(kind=INT64)                  :: j,k,l,m,maxphases
  double precision,dimension(2)        :: dt,times         ! dt = (dt_new, dt_old)
  character(len=200)                   :: err_msg
  logical                              :: check_extension, ext_VC0, use_varint00, use_sig00 

  double precision, dimension(6)                   :: sig_tmp
  double precision, dimension(3,3)                 :: DFGRD0,DFGRD1
  double precision, dimension(algo_param%nTensDef) :: def0_tmp, ddef_tmp

  double precision,allocatable, dimension(:)       :: Coeff, VarInt
  integer(kind=INT64)                              :: nVarInt, nCoeff, nVarint_umatcall
                                                   ! if (nVarint = 0, nVarint_umat_call = 1)
                                                   ! else nVarint_umatcall=nVarint (compatibility with MFRONT/UMAT/CAST3M)
  integer(kind=INT64)                              :: nVoxComp,nVoxComp_tot
  integer(kind=INT32)                              :: temp_int1,temp_int2
  integer                                          :: ierror
  double precision                                 :: t1 !mesure de temps 

  integer(kind=INT64)                              :: Nsub_loc,Nincr_loc,Nit_loc ! compteurs pour suivi des performances de Pilotage/umatLaminate


  call MPI_BARRIER(amitex_mpi_comm,ierror)          
  t1 = MPI_WTIME()

  KINC=1

  ! Chose stress used as input for umat 
  use_sig00 = .true.
  if (present(use_sig0)) use_sig00 = use_sig0
  if (use_sig00) SIG=SIG0

  dt(1)=t_load(0)-t_load(-1)
  dt(2)=t_load(-1)-t_load(-2)
  times=dble((/0._mytype,t_load(-1)/))

  ! tableau de coefficients et de variables internes (ici, temporairement, leur valeur max, utilisee pour allouer les tableaux)
  nCoeff=maxVal(MattotP(:)%Ncoeff)
  nVarInt=maxVal(MattotP(:)%NvarInt)
  
  ! Compteurs du nombre de voxels ayant eu recours à la subdivision du pas de temps pour l'integration du modèle  multi-couches 
  Nsub_loc = 0
  nSub_tot = 0
  ! Compteurs du nombre de sous pas de temps utilisés lors subdivision du pas de temps pour l'integration du modèle  multi-couches 
  Nincr_loc = 0
  nIncr_tot = 0
  ! Compteurs d'itérations par appel pour l'intégration du modèle multi-couches 
  Nit_loc      = 0
  nIt_tot      = 0
  ! Compteurs nombre total de voxels composites sur le pinceau
  nVoxComp     = 0
  nVoxComp_tot = 0

  allocate(Coeff(nCoeff));
  allocate(VarInt(max(1,nVarInt))) 
  Coeff=0
  VarInt=0

  ! Traitement du choix des variables internes transmises dans umat
  use_varint00 = .true.
  if (present(use_varint0)) use_varint00 = use_varint0

  ! Traitement du cas 'extension' 
  check_extension = .false.
  if (present(ext_law)) check_extension = .true.
  if (present(ext_lib)) call amitex_abort("Behavior(...) : use of variable ext_lib not implemented",2,0)
  
  ext_VC0=.true.
  if (present(ext_VC)) ext_VC0 = ext_VC

  !! Initialisation des grandeurs UMAT indep du mteriau TEMP, DTEMP, PREDEF, DPRED, KSTEP
  !! Si on entre des parametres externes en donnees on les envoie dans le UMAT
  !! Sinon on construit un tableau avec une valeur valant 0
  
  KSTEP = KSTEP0
  TEMP = dble(TEMPRED(1))
  DTEMP = dble(DTEMPRED(1))
  if (size(TEMPRED) .NE. size(DTEMPRED)) call amitex_abort("TEMPRED and DTEMPRED have diffrent sizes (behavior)",2)

  if(size(TEMPRED)>1) then
     if (.not. allocated(PREDEF)) allocate(PREDEF(size(TEMPRED)-1))
     PREDEF(:)=dble(TEMPRED(2:))
     if (.not. allocated(DPRED)) allocate(DPRED(size(DTEMPRED)-1))
     DPRED(:)=dble(DTEMPRED(2:))
  else
     if (.not. allocated(PREDEF)) allocate(PREDEF(1))
     PREDEF=0
     if (.not. allocated(DPRED)) allocate(DPRED(1))
     DPRED=0
  end if

#ifdef OPENMP
  !$OMP PARALLEL
  nbthread = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
#endif

  ! on récupère le plus grand nombre de phases des matériaux du pinceau
  maxphases = maxval(mattotP(:)%Nphase)

  ! TODO : test a faire sur les entrees
  if ((maxphases > 1) .and. (.not.algo_param%HPP)) then
     write(err_msg,"(A)") "Voxels composites en GD non implémentés : utiliser HPP"
     call amitex_abort(err_msg,1)
  end if


  ! allocation du tableau de pointeurs de fonctions
  allocate(umatptr(maxphases)) 
  
  do i=1,size(mattotP)
     ! En cas de matériau d'interphase : rien a faire, le materiau n'existe que comme constituant de materiaux composites
     if (mattotP(i)%Interphase) then
        cycle
     end if

     !affectation du tableau de pointeurs de fonction (cas "composite" ou "homogene", "classique" et "extension")
     !--------------------------- CAS COMPOSITE
     if (mattotP(i)%NPhase > 1) then  
       ! incrémentation du nombre de voxels composites multi-couches ou reuss sur le pinceau
       if (trim(mattotP(i)%lawname) == "laminate") nVoxComp = nVoxComp + size(MattotP(i)%zone(:,1))
       if (trim(mattotP(i)%lawname) == "reuss") nVoxComp = nVoxComp + size(MattotP(i)%zone(:,1))
       ! pas de bug ici en presence de zones interphases --> les mattot consideres sont composites 
       ! ils n'ont donc pas de zones interphases
       
       !------- cas general
       if (.NOT. check_extension) then 
         do j=1,mattotP(i)%NPhase
           umatptr(j) = umat_load(trim(mattotP(i)%Libname_comp(j))//C_NULL_CHAR,&
                                  trim(MattotP(i)%Lawname_comp(j))//C_NULL_CHAR)
         end do
       !------- cas extension
       else 
         do j=1,mattotP(i)%NPhase
           umatptr(j) = umat_load(trim(mattotP(i)%Libname_comp(j))//C_NULL_CHAR,&
                                  trim(MattotP(i)%Lawname_comp(j))//"_"//trim(ext_law)//C_NULL_CHAR)
         end do
       end if
     !--------------------------- CAS HOMOGENE
     else 
       !------- cas general
       if (.NOT. check_extension) then
         umatptr(1) = umat_load(trim(MattotP(i)%Libname)//C_NULL_CHAR,&
                                trim(MattotP(i)%Lawname)//C_NULL_CHAR)         
       !------- cas extension
       else
         umatptr(1) = umat_load(trim(MattotP(i)%Libname)//C_NULL_CHAR,&
                                trim(MattotP(i)%Lawname)//"_"//trim(ext_law)//C_NULL_CHAR)
       end if
     end if

     ! nombre de coefficient et de variables internes du materiau
     nCoeff=MattotP(i)%Ncoeff
     nVarInt=MattotP(i)%NvarInt
     nVarint_umatcall = nVarint

     ! indice minimum de zone
     l=1
     do j=1,size(MattotP(i)%zone(:,1),kind=INT64)
        
        if (allocated(mattotP(i)%Zones_interphase)) then
           if (any(mattotP(i)%Zones_interphase == mattotP(i)%zone(j,2))) then
              ! si zone interphase : cycle, pas de calculs
              cycle
           end if
        end if

        Coeff(1:nCoeff)=dble(MattotP(i)%Coeff(:,j))
        
#ifdef OPENMP
        !$OMP PARALLEL NUM_THREADS(nbthread) private(m,sig_tmp,ddef_tmp,def0_tmp,VarInt)
        !$OMP DO
#endif
        do k=l,MattotP(i)%zone(j,1)

           m=MattotP(i)%pos(k)
           sig_tmp=dble(sig(m,1:6))
           if(nVarInt>0) then
              if (use_varint00) then
                 VarInt(1:nVarInt)=dble(MattotP(i)%VarInt0(:,k))   ! default case use varint0 as umat input
              else
                 VarInt(1:nVarInt)=dble(MattotP(i)%VarInt(:,k))
              end if 
           else
              VarInt=0
              nVarint_umatcall = 1 ! Pour compatibilite MFRONT (et compatibilite umat/CAST3M) : nVarint >= 1
           end if
           if(algo_param%HPP) then
              ddef_tmp=dble(def(m,1:algo_param%nTensDef)-def0(m,1:algo_param%nTensDef))
              def0_tmp=dble(def0(m,1:algo_param%nTensDef))
              if (algo_param%HPP_nsym) then
                  !! prise en compte complete du tenseur gradient du deplacement 
                  !! dans un calcul HPP 
                  DFGRD0(1,1) = dble(1.+Def_nsym0(m,1))
                  DFGRD0(2,1) = dble(Def_nsym0(m,7))
                  DFGRD0(3,1) = dble(Def_nsym0(m,8))
                  DFGRD0(1,2) = dble(Def_nsym0(m,4))
                  DFGRD0(2,2) = dble(1.+Def_nsym0(m,2))
                  DFGRD0(3,2) = dble(Def_nsym0(m,9))
                  DFGRD0(1,3) = dble(Def_nsym0(m,5))
                  DFGRD0(2,3) = dble(Def_nsym0(m,6))
                  DFGRD0(3,3) = dble(1.+Def_nsym0(m,3))
                  DFGRD1(1,1) = dble(1.+Def_nsym(m,1))
                  DFGRD1(2,1) = dble(Def_nsym(m,7))
                  DFGRD1(3,1) = dble(Def_nsym(m,8))
                  DFGRD1(1,2) = dble(Def_nsym(m,4))
                  DFGRD1(2,2) = dble(1.+Def_nsym(m,2))
                  DFGRD1(3,2) = dble(Def_nsym(m,9))
                  DFGRD1(1,3) = dble(Def_nsym(m,5))
                  DFGRD1(2,3) = dble(Def_nsym(m,6))
                  DFGRD1(3,3) = dble(1.+Def_nsym(m,3))
              else
                  !! Si on veut evaluer une loi GDEF avec une resolution en HPP on passe
                  !! le tenseur gradient de la transformation (symmetrique) ---------- A VOIR !!!!!!
                  DFGRD0(1,1) = dble(1.+def0(m,1))
                  DFGRD0(2,1) = dble(def0(m,4))/2._mytype
                  DFGRD0(3,1) = dble(def0(m,5))/2._mytype
                  DFGRD0(1,2) = dble(def0(m,4))/2._mytype
                  DFGRD0(2,2) = dble(1.+def0(m,2))
                  DFGRD0(3,2) = dble(def0(m,6))/2._mytype
                  DFGRD0(1,3) = dble(def0(m,5))/2._mytype
                  DFGRD0(2,3) = dble(def0(m,6))/2._mytype
                  DFGRD0(3,3) = dble(1.+def0(m,3))
                  DFGRD1(1,1) = dble(1.+def(m,1))
                  DFGRD1(2,1) = dble(def(m,4))/2._mytype
                  DFGRD1(3,1) = dble(def(m,5))/2._mytype
                  DFGRD1(1,2) = dble(def(m,4))/2._mytype
                  DFGRD1(2,2) = dble(1.+def(m,2))
                  DFGRD1(3,2) = dble(def(m,6))/2._mytype
                  DFGRD1(1,3) = dble(def(m,5))/2._mytype
                  DFGRD1(2,3) = dble(def(m,6))/2._mytype
                  DFGRD1(3,3) = dble(1.+def(m,3))
              end if   

              ! APPEL AU COMPORTEMENT (VOXELS COMPOSITES OU HOMOGENES)
              if (mattotP(i)%NPhase > 1 .AND. (.NOT. check_extension)) then
                 !VOXEL COMPOSITE classique (UMAT sans extension)
                 
                 !! TODO : Maintenant qu'on n'appelle plus que pilotage_composite on pourrait s'affranchir
                 !!        de la distinction de cas ci-dessous if(voigt/reuss/laminate)
                 
                 if (trim(mattotP(i)%lawname) == "voigt") then
                    ! Loi d'homogénéisation de VOIGT
                    call PilotageComposite(i,umatptr(1:mattotP(i)%NPhase),sig_tmp,VarInt(1:nVarInt),&
                                   def0_tmp,ddef_tmp,times,dt,&
                                   nVarInt,Coeff(1:nCoeff),nCoeff,DFGRD0,DFGRD1,KINC,temp_int1,temp_int2)
                    if(KINC < 1) then
                       write(err_msg,"(A,I0,A)") "Probleme lors de l execution de la loi UMATVOIGT &
                            & pour le matériau ",mattotP(i)%numM," (behavior)"
                       call amitex_abort(err_msg,1)
                    end if
                 else if  (trim(mattotP(i)%lawname) == "laminate") then
                    ! Loi d'homogénéisation MULTI-COUCHES
                    ! temp_int1 -> nbre d'intération moyen umatLaminate au cours de l'intégration
                    ! temp_int2 -> nbre de sous pas de temps ayant été utilisés pour l'intégration (0 = pas de subdivision du pas de temps)
                    call  PilotageComposite(i,umatptr(1:mattotP(i)%Nphase),sig_tmp,VarInt(1:nVarInt),&
                                          def0_tmp,ddef_tmp,times,dt,&
                                          nVarInt,Coeff(1:nCoeff),nCoeff,DFGRD0, DFGRD1,KINC,temp_int1,temp_int2)

                    ! Incrémentation des compteurs pour chaque voxel de chaque matériau composite pour tout le pinceau
                    Nit_loc  = Nit_loc + temp_int1 !nombre d'itérations algo multicouche
                    if (temp_int2 > 0) Nsub_loc = Nsub_loc + 1 ! Nbre de voxels ayant necessité la subdivision du pas de temps
                    Nincr_loc = Nincr_loc + temp_int2 ! Nbre de sous pas temps total

                    if((KINC < 1) .and. (KINC /=-10)) then
                       write(err_msg,"(A,I0,A)") "Probleme lors de l'execution de la loi UMATLAMINATE &
                            & pour le matériau ",mattotP(i)%numM," (behavior)"
                       call amitex_abort(err_msg,1)
                    end if
                 else if (trim(mattotP(i)%lawname) == "reuss") then
                    ! Loi d'homogénéisation de REUSS
                    ! TODO (?) : prevoir de passer par PilotageVoxcomp (PilotageLaminate legerement modifie)
                     call PilotageComposite(i,umatptr(1:mattotP(i)%NPhase),sig_tmp,VarInt(1:nVarInt),&
                                   def0_tmp,ddef_tmp,times,dt,&
                                   nVarInt,Coeff(1:nCoeff),nCoeff,DFGRD0,DFGRD1,KINC,temp_int1,temp_int2)
                    !TODO : gestion des erreurs
                    if(KINC < 1) then
                       write(err_msg,"(A,I0,A)") "Probleme lors de l execution de la loi UMATREUSS &
                            & pour le matériau ",mattotP(i)%numM," (behavior)"
                       call amitex_abort(err_msg,2)
                    end if
                 end if
              elseif (mattotP(i)%NPhase > 1 .AND. check_extension .AND. ext_VC0) then
                 !VOXEL COMPOSITE : appel du comportement UMAT "extension"
                 call  umatComposite_ext(i,umatptr(1:mattotP(i)%NPhase),sig_tmp,VarInt(1:nVarInt),&
                                   def0_tmp,ddef_tmp,times,dt,&
                                   nVarInt,Coeff(1:nCoeff),nCoeff,DFGRD0,DFGRD1,KINC)
              elseif (mattotP(i)%NPhase == 1 ) then
                 !VOXEL HOMOGENE : appel du comportement UMAT ("classique" ou "extension")
                 call umat_call(umatptr(1),sig_tmp,VarInt(1:nVarInt),ddsdde,&
                      sse,spd,scd, rpl, ddsddt, drplde, drpldt,&
                      def0_tmp , ddef_tmp,times,&
                      dt(1),TEMP,DTEMP,PREDEF,DPRED,&
                      CMNAME,NDI,NSHR,NTENS,nVarInt_umatcall,&
                      Coeff(1:nCoeff),nCoeff,coords,&
                      DROT,pnewdt,celent,DFGRD0, DFGRD1,NOEL, NPT, LAYER, KSPT, KSTEP,kinc)
                 !TODO traitement erreur (KINC) cas homogene
              end if
           else ! HPP else GDEF (PAS ENCORE DE VOXELS COMPOSITES)
              !! Calcul de Id+Grad(u) a partir des Grad(u)
              !! def suit l'ordre 11 22 33 12 13 23 21 31 32
              DFGRD0(1,1) = dble(1.+def0(m,1))
              DFGRD0(2,1) = dble(def0(m,7))
              DFGRD0(3,1) = dble(def0(m,8))
              DFGRD0(1,2) = dble(def0(m,4))
              DFGRD0(2,2) = dble(1.+def0(m,2))
              DFGRD0(3,2) = dble(def0(m,9))
              DFGRD0(1,3) = dble(def0(m,5))
              DFGRD0(2,3) = dble(def0(m,6))
              DFGRD0(3,3) = dble(1.+def0(m,3))
              DFGRD1(1,1) = dble(1.+def(m,1))
              DFGRD1(2,1) = dble(def(m,7))
              DFGRD1(3,1) = dble(def(m,8))
              DFGRD1(1,2) = dble(def(m,4))
              DFGRD1(2,2) = dble(1.+def(m,2))
              DFGRD1(3,2) = dble(def(m,9))
              DFGRD1(1,3) = dble(def(m,5))
              DFGRD1(2,3) = dble(def(m,6))
              DFGRD1(3,3) = dble(1.+def(m,3))
              !! Si on veut evaluer une loi HPP avec une resolution en GDEF on passe
              !! le tenseur de deformation linearise ---------- A VOIR !!!!!!
              def0_tmp(1) = dble(def0(m,1))
              def0_tmp(2) = dble(def0(m,2))
              def0_tmp(3) = dble(def0(m,3))
              def0_tmp(4) = DFGRD0(1,2)+DFGRD0(2,1)
              def0_tmp(5) = DFGRD0(1,3)+DFGRD0(3,1)
              def0_tmp(6) = DFGRD0(2,3)+DFGRD0(3,2)
              ddef_tmp(1) = dble(def(m,1)) - def0_tmp(1)
              ddef_tmp(2) = dble(def(m,2)) - def0_tmp(2) 
              ddef_tmp(3) = dble(def(m,3)) - def0_tmp(3)
              ddef_tmp(4) = DFGRD1(1,2)+DFGRD1(2,1) - def0_tmp(4)
              ddef_tmp(5) = DFGRD1(1,3)+DFGRD1(3,1) - def0_tmp(5)
              ddef_tmp(6) = DFGRD1(2,3)+DFGRD1(3,2) - def0_tmp(6)

              call umat_call(umatptr(1),sig_tmp,VarInt(1:nVarInt),ddsdde,&
                    sse,spd,scd, rpl, ddsddt, drplde, drpldt,&
                    def0_tmp , ddef_tmp,times,&
                    dt(1),TEMP,DTEMP,PREDEF,DPRED,&
                    CMNAME,NDI,NSHR,NTENS,nVarInt_umatcall,&
                    Coeff(1:nCoeff),nCoeff,coords,&
                    DROT,pnewdt,celent,DFGRD0, DFGRD1,NOEL, NPT, LAYER, KSPT, KSTEP,kinc)
              !TODO traitement erreur (KINC) cas homogene
           end if  ! FIN 

           ! AFFECTATION DES DE "CHAMPS" CONTRAINTE ET VARIABLE INTERNE 
           SIG(m,1:6)=real(sig_tmp,mytype)
           if(nVarint>0)then
              MattotP(i)%VarInt(:,k)= real(VarInt(1:nVarInt),mytype)
           end if

           ! TRAITEMENT DES ERREURS
           if(KINC == -10) then
              write(err_msg,"(A)") "Probleme lors de l'execution de la loi UMATLAMINATE,&
                                   & plus de 200 itérations de l'algorithme, calcul arrêté  (behavior)"
              call amitex_abort(err_msg,1)
           end if
           if(KINC == -11) then
              write(err_msg,"(A)") "Probleme lors de l'execution de la loi UMATREUSS,&
                                   & plus de 200 itérations de l'algorithme, calcul arrêté  (behavior)"
              call amitex_abort(err_msg,1)
           end if

           !! test isnan sur sig et varint en sortie de umat :
           do h=1,6
              if (ieee_is_nan(sig_tmp(h))) then
                 call amitex_abort("Contrainte a la sortie de umat NaN (behavior)",2,0)
                 stop
              end if
           end do
           do h=1,int(nVarint,INT32)
              if (ieee_is_nan(VarInt(h))) then
                 write(err_msg,"(I4)") h
                 call amitex_abort("Variable interne "//trim(err_msg)//"a la sortie de umat NaN (behavior)",2,0)
                 stop
              end if
           end do


        end do
#ifdef OPENMP
        !$OMP END DO
        !$OMP END PARALLEL
#endif

        l=MattotP(i)%zone(j,1)+1
     end do

  end do
  ! Incrémentation des compteurs sur tous les voxels composites de la cellule
  ! Nombre total de voxels composite sur la cellule
  call MPI_Allreduce(nVoxComp, nVoxComp_tot, 1, MPI_INTEGER8, MPI_SUM, amitex_mpi_comm, ierror)
  ! Nombre total d'itération de l'algorithme multi-couche
  call MPI_Allreduce(Nit_loc, nIt_tot, 1, MPI_INTEGER8, MPI_SUM, amitex_mpi_comm, ierror)
  ! Nombre total de voxels ayant eu recours à la subdivision du pas de temps pour l'integration du modèle  multi-couches
  call MPI_Allreduce(Nsub_loc, Nsub_tot, 1, MPI_INTEGER8, MPI_SUM, amitex_mpi_comm, ierror)
  ! Nombre total  de sous pas de temps utilisés lors subdivision du pas de temps pour l'integration du modèle  multi-couches
  call MPI_Allreduce(Nincr_loc, nIncr_tot, 1, MPI_INTEGER8, MPI_SUM, amitex_mpi_comm, ierror)

  if (present(nVoxComp_tot_s)) nVoxComp_tot_s = nVoxComp_tot

  if (allocated(Coeff)) deallocate(Coeff)
  if (allocated(VarInt)) deallocate(VarInt)
  if (allocated(PREDEF)) deallocate(PREDEF)
  if (allocated(DPRED)) deallocate(DPRED)
  if (allocated(umatptr)) deallocate(umatptr)
  call check_amitex_abort()
  
  call MPI_BARRIER(amitex_mpi_comm,ierror)  
  times_m%behavior = times_m%behavior + MPI_WTIME() - t1        

end subroutine behavior

!==============================================================================
!
!                      SUBROUTINE BEHAVIORD
!
!> Calcul du comportement du materiau pour la diffusion (GradQD -> FluxD)
!!
!!  Appel a la fonction umat pour chacun des voxels de la cellule
!!       l'ordre de parcours respecte la structure MATERIAL.
!!
!!  
!!  \param[out]    FluxD: champ des tenseurs des flux (a l'instant t,ici fin de pas de temps)
!!  \param[in]     FluxD0: champ des tenseurs flux de diffusion (au debut du chargement)
!!  \param[in]     GradQD: champ des tenseurs gradient des variables de diffusion (a l'instant t)
!!  \param[in]     GradQD0: champ des tenseurs de deformations (au debut du chargement)
!!  \param[in]     TEMPRED: temperature et parametres externes de la loi de comportement
!!  \param[in]     t_load: temps aux instants t, t-1 et t-2
!!  \param[in]     KSTEP indice de chargement 
!!
!! Modifie egalement les variables internes dans MattotP
!!
!==============================================================================
subroutine behaviorD(FluxD, FluxD0, GradQD, GradQD0, TEMPRED,DTEMPRED,t_load,KSTEP0)

  implicit none

  real(mytype), dimension(amx_size(1)*amx_size(2)*amx_size(3),3,algo_param%nVarD),intent(out) :: FluxD 
  real(mytype), dimension(amx_size(1)*amx_size(2)*amx_size(3),3,algo_param%nVarD),intent(in)  :: FluxD0,GradQD, GradQD0
  real(mytype),dimension(-2:0), intent(in)                                           :: t_load
  real(mytype),dimension(:), intent(in)                                              :: TEMPRED,DTEMPRED
  integer(kind=8),intent(in)                                                         :: KSTEP0

  type(c_ptr),allocatable,dimension(:) :: umatptr    ! pointeur de fonction
  integer(kind=8)                      :: KINC
  integer            :: i,h
  integer(kind=INT64):: i0                !INT64 : compatibilite nVarD
#ifdef OPENMP
  integer            :: nbthread
#endif
  integer(kind=8)    :: j,k,l,m
  double precision   :: dt
  double precision,dimension(2)   :: times
  character(len=200) :: err_msg

  double precision, dimension(3,algo_param%nVarD) :: fluxd_tmp
  double precision, dimension(3,algo_param%nVarD) :: gradQD0_tmp, dGradQD_tmp

  double precision,allocatable, dimension(:) :: Coeff!,nVarInt
  integer(kind=INT64) :: nCoeff,nVarD      !INT64 : compatibiite umat

  integer                     :: ierror
  double precision            :: t1
  
  call MPI_BARRIER(amitex_mpi_comm,ierror)          
  t1 = MPI_WTIME()

  KINC = 1

  !Nombre de variables de diffusion (INT64 pour interface umat)
  nVarD = int(algo_param%nVarD,INT64)

  FluxD = FluxD0
  dt=t_load(0)-t_load(-1)
  times=dble((/0._mytype,t_load(-1)/))

  ! tableau de coefficients et de variables internes
  nCoeff=int(maxVal(MattotP(:)%NcoeffK),INT64)


  allocate(Coeff(nCoeff));
  !allocate(VarInt(nVarInt));
  Coeff=0
  !VarInt=0

  !! Si on entre des parametres externes en donnees on les envoie dans le UMAT
  !! Sinon on construit un tableau avec une valeur valant 0
  KSTEP=KSTEP0
  TEMP = dble(TEMPRED(1))
  DTEMP = dble(DTEMPRED(1))
  if (size(TEMPRED) .NE. size(DTEMPRED)) call amitex_abort("TEMPRED and DTEMPRED have diffrent sizes (behavior)",2)

  if(size(TEMPRED)>1) then
     if (.not. allocated(PREDEF)) allocate(PREDEF(size(TEMPRED)-1))
     PREDEF(:)=dble(TEMPRED(2:))
     if (.not. allocated(DPRED)) allocate(DPRED(size(DTEMPRED)-1))
     DPRED(:)=dble(DTEMPRED(2:))
  else
     if (.not. allocated(PREDEF)) allocate(PREDEF(1))
     PREDEF=0
     if (.not. allocated(DPRED)) allocate(DPRED(1))
     DPRED=0
  end if


#ifdef OPENMP
  !$OMP PARALLEL
  nbthread = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
#endif

  ! allocation du tableau de pointeurs de fonctions (pour voxels composites)
  allocate(umatptr(1)) 

  do i=1,size(mattotP)

     !pointeur de fonction
     umatptr(1) = umat_load(trim(MattotP(i)%LibnameK)//C_NULL_CHAR,trim(MattotP(i)%LawnameK)//C_NULL_CHAR)

     ! nombre de coefficient et de variables internes du materiau
     nCoeff=MattotP(i)%NcoeffK
     !nVarInt=MattotP(i)%NvarInt

     ! indice minimum de zone
     l=1
     !MattotP(i)%VarInt=MattotP(i)%VarInt0
     do j=1,size(MattotP(i)%zone(:,1))

        Coeff(1:nCoeff)=dble(MattotP(i)%CoeffK(:,j))

#ifdef OPENMP
        !$OMP PARALLEL NUM_THREADS(nbthread) private(m,FluxD_tmp,dGradQD_tmp,GradQD0_tmp)
        !$OMP DO
#endif
        do k=l,MattotP(i)%zone(j,1)

           m=MattotP(i)%pos(k)
           !sig_tmp=dble(sig(m,1:6))
           do i0=1,nVarD
               FluxD_tmp(:,i0)=dble(FluxD(m,1:3,i0))
           end do
           !if(nVarInt>1) then
           !   VarInt(1:nVarInt)=dble(MattotP(i)%VarInt(:,k))
           !else
           !   VarInt=0
           !end if

           do i0=1,nVarD
              dGradQD_tmp(:,i0)=dble(gradQD(m,1:3,i0)-GradQD0(m,1:3,i0))
              GradQD0_tmp(:,i0)=dble(GradQD0(m,1:3,i0))
           end do

           call umatD_call(umatptr(1),fluxd_tmp,gradQD0_tmp,dgradQD_tmp,times,&
                    dt,TEMP,DTEMP,PREDEF,DPRED,NVARD,Coeff(1:nCoeff),nCoeff,kinc)

           if(KINC < 1) then
              write(err_msg,"(3A,I0,A)") "Probleme lors de l'execution de la loi UMATD ",&
                   CMNAME," KINC = ",KINC," verifier que la loi est bien implementee dans &
                   &ce programme (behaviorD)"
              call amitex_abort(err_msg,1)
           end if
           do i0=1,nVarD
              FluxD(m,:,i0)=real(FluxD_tmp(1:3,i0),mytype)
           end do

           !if(nVarint>1)then
           !   MattotP(i)%VarInt(:,k)= real(VarInt(1:nVarInt),mytype)
           !end if

           !! test isnan sur FluxD en sortie de umat :
           do i0=1,nVarD
           do h=1,3
              if (ieee_is_nan(FluxD_tmp(h,i0))) then
                 call amitex_abort("Flux a la sortie de umatD NaN (behavior)",2,0)
              end if
           end do 
           end do
           !do h=1,nVarint
           !   if (ieee_is_nan(VarInt(h))) then
           !      call amitex_abort("Variable interne a la sortie de umat NaN (behavior)",2,0)
           !   end if
           !end do 

        end do
#ifdef OPENMP
        !$OMP END DO
        !$OMP END PARALLEL
#endif

        l=MattotP(i)%zone(j,1)+1
     end do
  end do

  if (allocated(Coeff)) deallocate(Coeff)
  if (allocated(PREDEF)) deallocate(PREDEF)
  if (allocated(DPRED)) deallocate(DPRED)
  if (allocated(umatptr)) deallocate(umatptr)
  call check_amitex_abort()
  
  call MPI_BARRIER(amitex_mpi_comm,ierror)  
  times_m%behavior = times_m%behavior + MPI_WTIME() - t1        

end subroutine behaviorD

!==============================================================================
!
!               SUBROUTINE INITPARAM_UMAT
!
!>  Initialisation des parametres de la fonction UMAT fixes sur le pas de temps
!
!
!
!==============================================================================
subroutine initParam_umat()

      ! STRESS  
      ! STATEV
      ddsdde = 0._mytype ! matrice jacobienne (ntens,ntens)
      sse = 0._mytype
      spd = 0._mytype
      scd = 0._mytype
      rpl = 0._mytype
      ddsddt = 0._mytype
      drplde = 0._mytype
      drpldt = 0._mytype
      !STRAN  
      !DSTRAN  
      !TIME (2)
      !DTIME
      TEMP = 0._mytype
      DTEMP = 0._mytype
      !PREDEF 
      !DPRED  
      CMNAME="0"
      NDI = 2
      NSHR = 3  ! nbre de composante extra-diagonale du tenseur de contrainte
      NTENS=6   ! pour umat, il s'agit du tenseur des contraintes de Cauchy
      !NSTATV 
      !PROPS  
      !NPROPS     
      COORDS = 0._mytype
      drot = 0._mytype     ! adherance CAST3M - matrice de passage repere local element - repere global
      drot(1,1)=1._mytype  ! adherence ABAQUS - matrice d'increments de rotation
      drot(2,2)=1._mytype  ! adherence AMITEX = adherence CAST3M avec repere local = repere global  
      drot(3,3)=1._mytype  !                   -> drot=Id
      pnewdt = 1e12_mytype  ! 1e12 for identical substepping behavior between MFRONT/AMITEX and MFRONT/CAST3M
      celent = 0._mytype
      !DFGRD0 
      !DFGRD1 
      NOEL = 0
      NPT = 0
      layer = 0
      kspt = 0
      KSTEP = 0
      !KINC 

end subroutine initParam_umat

!==============================================================================
!      SUBROUTINE DEALLOCATE_MATTOTP0
!
!>  DESALLOCATION DU TABLEAU DE MATERIAUX
!!
!!  INUTILE : dallocate(MattotP0) doit fonctionner
!
!------------------------------------------------------------------------------
subroutine deallocate_MattotP0()

  implicit none

!  integer :: i

! INUTILE : deallocate(MattotP0) fait tout ce travail
!  do i = 1, size(MattotP0)
!    if(allocated(MattotP0(i)%pos)) deallocate(MattotP0(i)%pos)
!    if(allocated(MattotP0(i)%Zone)) deallocate(MattotP0(i)%Zone)
!    if(allocated(MattotP0(i)%Coeff)) deallocate(MattotP0(i)%Coeff)
!    if(allocated(MattotP0(i)%varInt)) deallocate(MattotP0(i)%varInt)
!    if(allocated(MattotP0(i)%varInt0)) deallocate(MattotP0(i)%varInt0)
!    if(allocated(MattotP0(i)%numM_composite)) deallocate(MattotP0(i)%numM_composite)
!    if(allocated(MattotP0(i)%numZ_composite)) deallocate(MattotP0(i)%numZ_composite)
!    if(allocated(MattotP0(i)%Coeff_comp)) deallocate(MattotP0(i)%Coeff_comp)
!    if(allocated(MattotP0(i)%NCoeff_comp)) deallocate(MattotP0(i)%NCoeff_comp)
!    if(allocated(MattotP0(i)%NvarInt_comp)) deallocate(MattotP0(i)%NvarInt_comp)
!    if(allocated(MattotP0(i)%Libname_comp)) deallocate(MattotP0(i)%Libname_comp)
!    if(allocated(MattotP0(i)%Lawname_comp)) deallocate(MattotP0(i)%Lawname_comp)
!  end do

  deallocate(MattotP0)
end subroutine deallocate_MattotP0

!==============================================================================
!      SUBROUTINE DEALLOCATE_MATREF0
!
!>  DESALLOCATION DU TABLEAU DE MATERIAUX
!
!------------------------------------------------------------------------------
subroutine deallocate_Matref0()

  implicit none
  
  if(allocated(Matref0%C0)) deallocate(Matref0%C0)
  if(allocated(Matref0%K0D)) deallocate(Matref0%K0D)
  

end subroutine deallocate_Matref0

!==============================================================================
!      SUBROUTINE MIN_MPI_I8
!
!>  APPEL DE MPI_ALL_REDUCE, 1, INTEGER8, MIN 
!!  \param[in]     i : integer value distributed on pencils
!!  \param[out]    j : min(i)
!!  \param[out]    ierror : error output from MPI_all_reduce
!------------------------------------------------------------------------------
subroutine MIN_MPI_I8(i,j,ierror)
    implicit none
    integer(kind=8), intent(in) :: i
    integer(kind=8), intent(out):: j
    integer, intent(out)        :: ierror

    call MPI_Allreduce(i, j, 1, MPI_INTEGER8, MPI_MIN, amitex_mpi_comm, ierror)

end subroutine MIN_MPI_I8

!==============================================================================
!      SUBROUTINE MAX_MPI_I8
!
!>  APPEL DE MPI_ALL_REDUCE, 1, INTEGER8, MAX
!!  \param[in]     i : integer8 value distributed on pencils
!!  \param[out]    j : max(i)
!!  \param[out]    ierror : error output from MPI_all_reduce
!------------------------------------------------------------------------------
subroutine MAX_MPI_I8(i,j,ierror)
    implicit none
    integer(kind=8), intent(in) :: i
    integer(kind=8), intent(out):: j
    integer, intent(out)        :: ierror

    call MPI_Allreduce(i, j, 1, MPI_INTEGER8, MPI_MAX, amitex_mpi_comm, ierror)

end subroutine MAX_MPI_I8

!==============================================================================
!      SUBROUTINE MAX_MPI_I
!
!>  APPEL DE MPI_ALL_REDUCE, 1, INTEGER, MAX 
!!  \param[in]     i : integer value distributed on pencils
!!  \param[out]    j : max(i)
!!  \param[out]    ierror : error output from MPI_all_reduce
!------------------------------------------------------------------------------
subroutine MAX_MPI_I(i,j,ierror)
    implicit none
    integer, intent(in) :: i
    integer, intent(out):: j
    integer, intent(out):: ierror

    call MPI_Allreduce(i, j, 1, MPI_INTEGER, MPI_MAX, amitex_mpi_comm, ierror)

end subroutine MAX_MPI_I

!==============================================================================
!      SUBROUTINE SUM_MPI_I8
!
!>  APPEL DE MPI_ALL_REDUCE, 1, INTEGER8, SUM 
!!  \param[in]     i : integer8 value distributed on pencils
!!  \param[out]    j : sum(i)
!!  \param[out]    ierror : error output from MPI_all_reduce
!------------------------------------------------------------------------------
subroutine SUM_MPI_I8(i,j,ierror)
    implicit none
    integer(kind=8), intent(in) :: i
    integer(kind=8), intent(out):: j
    integer, intent(out)        :: ierror

    call MPI_Allreduce(i, j, 1, MPI_INTEGER8, MPI_SUM, amitex_mpi_comm, ierror)

end subroutine SUM_MPI_I8

!==============================================================================
!      SUBROUTINE SUM_MPI_I8_LARGE
!
!>  APPEL DE MPI_ALL_REDUCE, SIZE_VECT, INTEGER8, SUM 
!!  \param[in]     Vect_int_in : vector of integer8 (size N) )value distributed on pencils
!!  \param[in]     Size_vec    : size of the vector to communicate
!!  \param[out]    Vect_int_out : sum(i)
!!  \param[out]    ierror : error output from MPI_all_reduce
!!
!!  TO DO : implementer la gestion du cas ou la taille du vecteur a sommer
!!          depasse la taille du type int4. Il faut decouper le vecteur_in
!!          en sous blocs et les sommer en parallele dans le bloc correspondant
!!          du fichier out
!------------------------------------------------------------------------------
subroutine SUM_MPI_I8_LARGE(Vect_int_in,Vect_int_out,Size_vect,ierror)
    implicit none
    integer(kind=8),intent(in)                        :: Size_vect
    integer(kind=8),dimension(Size_vect),intent(in)   :: Vect_int_in
    integer(kind=8),dimension(Size_vect),intent(out)  :: Vect_int_out
    integer, intent(out)                              :: ierror
    
    integer(kind=8)                                   :: limit_size
    
    limit_size = 2**30-1
    if (Size_vect >= limit_size) then
        call amitex_abort("La taille du vecteur a sommer depasse la taille pouvant être communique par MPI_Allreduce.&
                           & La somme MPI par buffer n'est pas encore implementee (SUM_MPI_I8_LARGE)",2,0)
    else 
        call MPI_Allreduce(Vect_int_in, Vect_int_out, int(Size_vect,kind=4), MPI_INTEGER8, MPI_SUM, amitex_mpi_comm, ierror)
    end if

end subroutine SUM_MPI_I8_LARGE

!==================================================================================
!==================================================================================
!                         FONCTION TEST_COMPOSITE
!
!> Détecte  si le calcul s'effectue avec ou sans Voxels Composites
!!
!! A partir du fichier xml file_mat
!!       - Repère si la balise Material_composite est présente et renseignée dans le
!!         fichier
!!       - Renvoie la valeur vrai si la balise est présente et renseignée, faux sinon
!!
!!  \param[in]  file_mat: (chaine de caracteres) nom du fichier xml renseignant les proprietes materiau
!!
!!  \param[out] test_composite:	".true." ou ".false"
!!
!==================================================================================
  function test_composite(file_mat)

    implicit none

    character(len=*),intent(in)      :: file_mat
    logical                          :: test_composite

    type(node), pointer              :: fi
    type(nodeList), pointer          :: material

    fi => parseFile(trim(file_mat))

    material => getElementsByTagName(fi,"Material_composite")

    test_composite = .false.
    if (getLength(material) > 0) test_composite = .true.

  end function test_composite

!==================================================================================
!==================================================================================
!                         FONCTION GET_INTERPHASE_LIST
!> Détecte  si le calcul s'effectue avec ou sans Voxels Composites
!! Si le calcul 
!!
!! A partir du fichier xml file_mat
!!       - Repère si la balise Material_composite est présente et renseignée dans le
!!         fichier
!!       - Renvoie la valeur vrai si la balise est présente et renseignée, faux sinon
!!
!!  \param[in]  file_mat: (chaine de caracteres) nom du fichier xml renseignant les proprietes materiau
!!
!!  TEST_COMPOSITE[out] renvoie ".true." ou ".false"
!!
!==================================================================================
subroutine  get_interphase_list(file_mat)

  implicit none

  character(len=*),intent(in)                           :: file_mat
  integer                                               :: ninterphases
  integer                                               :: i,alloc_stat,numM,n,err
  integer,allocatable,dimension(:)                      :: tmp_array
  type(node), pointer                                   :: fi,cur_node
  type(nodeList), pointer                               :: interphase,mat_list,zone_list, sub_list  
  character(len=100000)                                 :: int_zone_list
  character(len=200)                                    :: err_msg

  fi         => parseFile(trim(file_mat))
  interphase => getElementsByTagName(fi,"Interphase")

  if (getLength(interphase) >0 ) then

     cur_node   => item(interphase,0)
     mat_list   => getElementsByTagName(cur_node,"Interphase_material")
     zone_list  => getElementsByTagName(cur_node,"Interphase_zone_list")
     ninterphases = getLength(mat_list) + getLength(zone_list)

     ! allocation de la structure INTERPHASES
     allocate(Interphases(ninterphases),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort("Structure Interphases : espace memoire disponible insuffisant (get_interphase_list)",2)
     
     ! Mise en donnée des matériau d'interphase 
     do i=1,getLength(mat_list)
        numM=-1
        cur_node => item(mat_list,i-1)
        call get_int_xml(cur_node,"numM",numM,1,0)  ! stockage du numéro de matériau concerne
        Interphases(i)%numM = numM
        Interphases(i)%all = .true.     ! balise reperant un matériau d'interphase : aucun voxel homogène constitué de ce matériau

        if (numM == -1) then           
           write(err_msg,fmt="(A,I4,A)") "Numero de materiau invalide noeud Interphase_material",i," (get_interphase_list)"
           call amitex_abort(err_msg,1)
        end if

        ! Récupération du nombre de zones concernées et allocation du tableau
        call get_int_xml(cur_node,"Nzones",n,1,0)
        Interphases(i)%Nzones = n;
     end do

     ! Mise en donnée des zones d'interphase
     do i=1,getLength(zone_list)
        numM=-1
        cur_node => item(zone_list,i-1) 
        ! stockage du numéro de matériau concerne
        call get_int_xml(cur_node,"numM",numM,1,0) 
        Interphases(i+getLength(mat_list))%numM = numM
        Interphases(i+getLength(mat_list))%all = .false.  ! balise reperant un materiau comprenant 
                                                          ! des zones interphases : seulement certaines
                                                          ! zones du matériau numM sont des interphases

        if (numM == -1) then           
           write(err_msg,fmt="(A,I4,A)") "Numero de materiau invalide noeud Interphase_zone_list",i&
                                         ," (get_interphase_list)"
           call amitex_abort(err_msg,1)
        end if

        ! Récupération du nombre de zones concernées et allocation du tableau
        call get_int_xml(cur_node,"Nzones",n,1,0)
        Interphases(i)%Nzones = n;
        allocate(Interphases(i+getLength(mat_list))%zones(n),stat=alloc_stat)
        if(alloc_stat /=0) call amitex_abort("Structure Interphases : espace memoire disponible insuffisant &
                                             &(get_interphase_list)",2)

        ! Récupération des numéros de zone concernés
        call getNodeList(cur_node,sub_list,"ZoneList",-1,1)
        if(getLength(sub_list)>0)then
           allocate(tmp_array(n),stat=alloc_stat)
           if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (get_interphase_list)",2)
           call extractDataContent(item(sub_list,0),int_zone_list)
           call get_int_vect( int_zone_list,tmp_array, n, err)
           Interphases(i+getLength(mat_list))%zones = tmp_array ! gcc-warning accepte
           deallocate(tmp_array)
        end if
        
     end do

  end if

end subroutine get_interphase_list

!==================================================================================
!==================================================================================
!                      FONCTION READ_COMPOSITE 
!> Lecture des fichiers de données composites indiqués dans le fichier materiau
!!        Les données sont stockées dans la structure 'Composite'
!!
!! A partir du fichier xml file_mat
!!   lit et affecte :
!! 
!!      A partir du fichier 'liste_composite_materials.txt'
!!            - le nombre de matériaux composites et leur composition (numéro des 
!!              matériaux homogènes qui les constituent
!!            - la loi d'intégration du comportement composite de chaque matériau
!!        
!!      Dans le répertoire du matériau : rep_i_j_k pour un matériau constitué des 
!!                                       matériaux i, j et k
!!            - la position linéairisée des voxels ayant ce comportement dans le 
!!              pinceau
!!            - les numéros de zones associés à chacune des phases constitutives
!!              par voxel 
!!            - Les fractions volumiques des phases constitutives par voxel
!!            - Les données géométriques éventuelles (vecteurs d'orientation de
!!              l'interface : Normale,Direction tangente, Surface des interfaces)
!!
!!  \param[in]  : file_mat
!!
!!   ATTENTION : la lecture peut-être très lente en fonction du réseau
!!
!!  Out : 
!!      Crée et modifie la structure Composite
!!
!!  MODIF LG 7/12/2020 : structure Composite ne stocke plus les champs pos_globale,
!!                                           Zone, Fv, Normale, Direction, Surface
!!
!==================================================================================

subroutine read_composite(file_mat)

  implicit none

  character(len=*), intent(in)                  :: file_mat

  integer                                       :: i, nligne, alloc_stat,local_offset
  integer(kind=INT64)                           :: nZones
  integer,dimension(10)                         :: numM
  type(node),pointer                            :: fi,cur_node
  type(nodeList),pointer                        :: material, CoeffComp
  character(len=200)                            :: lawname,ligne,ssdir,filename2
  character(len=400)                            :: err_msg,filename,dir
  character(len=200)                            :: entete
  integer(kind=4)                               :: nMat,ierror,Unit1,Unit2,ios
  logical                                       :: file_exists

  fi => parseFile(trim(file_mat))

  material => getElementsByTagName(fi,"Material_composite") ! récupération du noeud Material_composite
  cur_node => item(material,0) 

  CoeffComp => getElementsByTagName(cur_node,"Coeff_composite") ! récupération du noeud Coeff_composite
  cur_node => item(CoeffComp,0)

  call get_str_xml_nolowercase(cur_node,"directory",dir,1,0)  ! récupération du répertoire où sont stockées les fichiers composites
  
  filename = trim(dir) // '/list_composite_materials.txt'  ! reconstruction du nom de fichier liste_materials
  inquire(FILE=filename,EXIST=file_exists)
  if(.not.file_exists)then 
     write(err_msg,fmt="(3A)") "Le fichier: '",trim(filename),"' n'existe pas ,(read_composite )"
     call amitex_abort(err_msg,1,0)
     return
  end if

  ! ouverture du fichier contenant la définition des matériaux composites
  open(newunit=Unit1,file=trim(filename), action="read", form="formatted", access="sequential",&
           & position="rewind",iostat=ierror, status="old")

  ! gestion d'erreur à l'ouverture du fichier
  if (ierror /= 0)  then
     write(err_msg,fmt="(3(A))") "Echec de l'ouverture du fichier : ", trim(filename), "(read_composite)"
     call amitex_abort(err_msg,1)
  else 

     nligne = 0
     ! PREMIERE LECTURE du fichier : détermination du nombre de matériaux composites 
     !------------------------------------------------------------------------------
     do while (ierror==0)
        read(unit=Unit1, iostat=ierror, fmt="(A)",advance='YES') ligne
        if (ierror==0) nligne = nligne + 1
     end do

     ! allocation de la structure composite
     allocate(MatComposite(nligne),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_composite)",2)

     ! Retour au début du fichier
     rewind(unit=Unit1)
     nligne =0

     ! SECONDE LECTURE: liste des matériaux composites
     !------------------------------------------------
     do 
        nligne = nligne +1
        read(unit=Unit1, iostat=ierror, fmt="(A)",advance='YES') ligne
        if (ierror > 0) then
           ! gestion des erreurs de lecture
           write(err_msg,fmt="(3(A))") "Erreur lors de la lecture du fichier : ", trim(filename), "(read_composite)"
           call amitex_abort(err_msg,1)
 
        else if (ierror < 0) then
           ! fin de fichier
           exit

        else 
           ! lecture de la ligne : lecture de la loi composite
           !--------------------------------------------------
           i = verify(ligne,'0123456789 ') ! i : indice du premier caractère non numérique/non espace
           lawname = trim(ligne(i:))
           ! vérification : nom de loi précisé implémenté (voigt ou laminate) A CONSERVER ???
           ! autre check possible : nom de loi nom précisé 
           if (trim(lawname) /= 'voigt' .and. trim(lawname) /= 'laminate' .and. trim(lawname) /= 'reuss') then
              write(err_msg,fmt="(3(A))") "Loi composite ",trim(lawname),&
                                " non implémentée  : utiliser 'laminate', 'reuss' ou 'Voigt' (read_composite)"
              call amitex_abort(err_msg,1)
           end if

           ! vérification/warning : si un ou plusieurs  numéros matériaux indiqués après le nom de la loi utilisée
           ! ils sont ignorés
           if (scan(lawname,'0123456789') /= 0) then
              write(err_msg,fmt="(3(A),I0,(A))") "Numéros de matériaux constitutifs écrits &
                   & après le nom de loi composite dans le fichier  : ", &
                     trim(filename)," à la ligne",nligne, " ignorés (read_composite)"
              call amitex_abort(err_msg,0)
           end if

           ! Affectation de la loi de comportement composite à la structure composite
           MatComposite(nligne)%Lawname = trim(lawname)

           ! on ne garde que les numéros de matériau dans ligne
           ligne = ligne(1:i-1)
           ! lecture de la ligne : lecture des numéros de matériaux constitutifs
           !--------------------------------------------------------------------
           i = 1
           nMat = 1

           ! initialisation du nom de sous-répertoire contenant les données du matériau composite actuel
           ssdir = "rep"

           do while ((i /= -1) .and. (i<len(trim(ligne))+1) .and. (nMat <= 10))
              ! récupération du numéro de phase suivant
              call nextint(ligne,i,numM(nMat))   
              write(ssdir,"(A,A,I1)")trim(ssdir),"_",numM(nMat)
              if (i<len(trim(ligne))+1) nMat = nMat + 1
              
              ! warning : on ne prends pas en compte des composites à plus de 10 phases (?????)
              if (nMat == 10) then
                 write(err_msg,fmt="(3(A),I0,(A))") "Plus de 10 phases renseignée dans le fichier  : ",&
                            filename," à la ligne",nligne,&
                            " les phases suivantes seront ignorées pour ce matériau (read_composite)"
                 call amitex_abort(err_msg,0)
              end if
           end do

           ! allocation du tableau Phases de la structure Composite puis affectation
           allocate(MatComposite(nligne)%num_Phases(nMat),stat=alloc_stat)
           if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_composite)",2)
           MatComposite(nligne)%num_Phases(:) = numM(1:nMat)
           

           ! Lecture du répertoire propre au matériau composite : données géométriques
           !--------------------------------------------------------------------------
           
           ! Fichier contenant le vecteur position (indices) linéarisée des voxels
           filename2 = trim(dir) //'/' // trim(ssdir) // "/pos.bin"
           matComposite(nligne)%dir = trim(dir) //'/' // trim(ssdir) 

           ! Lecture du nombre de zones composite, pour vérification du nombre de valeurs dans 
           ! chaque fichier composite, et compatibilité avec GETBINCOEFFFROMFILE           
           !! On verifie que le fichier existe
           inquire(FILE=filename2,EXIST=file_exists)
           if(.not.file_exists)then 
              write(err_msg,fmt="(3A)") "Le fichier: '",trim(filename2),"' n'existe pas ,(read_composite )"
              call amitex_abort(err_msg,1,0)
              return
           end if
           !! On ouvre le fichier pour lire l'en-tete
           open(newunit=Unit2, file=trim(filename2),form="formatted", status="old", action="read",iostat= ios)
           if ( ios /= 0 ) then
              write(err_msg,fmt="(3A,I0,A)") "Probleme a l'ouverture du fichier : ",trim(filename2),&
                   " , ",ios," (read_composite)"
              call amitex_abort(err_msg,2,0) 
           end if
           !! Lecture et stockage de l'en-tete (les 2 premieres lignes) dans le tableau entete
           read(Unit2,FMT='(A)') entete
           !! Fermeture du fichier d'entree
           close(unit=Unit2)
           !! Recuperation du nombre de zones dans le fichier binaire
           local_offset=1
           call nextInt8(entete,local_offset,nZones)

           matComposite(nligne)%Nzones_global = nZones


        end if

     end do
  end if

  call check_amitex_abort()

end subroutine read_composite

!==================================================================================
!==================================================================================
!                         SUBROUTINE READ_MAT_COMPOSITE
!> Definition des materiaux
!!
!! A partir du fichier xml file_mat
!!  lit et affecte:
!!       - les proprietes du materiau de reference (lambda0, mu0)
!!       - nom de librairie et de loi pour chacun des materiaux (MattotP0(i)%Libname et MattotP0%Lawname )
!!       - les coefficients et variables internes ( MattotP0(i)%Coeff, MattotP0(i)%VarInt )
!!
!! Version de READ_MAT adaptée au traitement des matériaux_composite. Gère notamment :
!!       - l'initialisation des données des matériaux homogènes dans le pinceau 
!!         de façon identique à READ_MAT
!!       - l'initialisation des coefficients et variables internes des matériaux composites
!!         à partir des données initialisées dans les matériaux homogènes et de la structure
!!         MatComposite
!!       - la récupération de données dans des structures matériaux présentent dans d'autres pinceaux
!!         (cas de composites à la frontière du pinceau ou des composites avec interphases)
!!       - Initialisation des matériaux d'interphase
!!
!!  \param[in]  file_mat: (chaine de caracteres) nom du fichier xml renseignant les proprietes materiau
!!  \param[in] nmate: (entier) nombre de materiaux
!!
!! \param[out] lambdamu0: (reels) proprietes du materiau de reference\n
!!
!! Modifie aussi MattotP0
!!
!!
!!
!! CORRECTION LE 4/4/2019 (Aldo Marano):
!!        la lecture parallele binaire (getbincoefffromfile via get_vect_xml)
!!        pouvait restee bloquer dans le cas de materiau absent sur un pinceau
!!        car get_vect_xml n'etait pas appelee
!!        -> necessite de traiter ce cas 
!!      VOIR AUSSI CORRECTIONS get_vect_xml et getbincoeffomfile
!!
!! ATTENTION LE 8/11/2019 (LG) ;
!!        La lecture des variables internes des materiaux 'composites' n'est pas compatible avec
!!        la nouveau Type d'initialisation Type="Variable"
!!        => Pas de plantage attendu MAIS :
!!              on affecte au materiau composite defini par (numM,numZ) 
!!              la variable interne du dernier point de la zone numZ du materiau numM
!!
!! MODIFICATION LE 7/12/2020 (LG) :
!!        Les champs MatComposite(i)%pos_global, zone, fv etc... ne sont plus stockés
!!        On lit les champs "lorsque c'est necessaire"
!!
!! \TODO Prendre en compte les cas ou il y a plus de coefficients que de zones :
!!       on pourrait alors s'arreter de lire des que le nombre de zones est atteint
!==================================================================================
subroutine read_mat_composite(file_mat, lambdamu0, nmate,nmate_composites)

  implicit none

  character(len=*), intent(in)                                    :: file_mat
  real(mytype),dimension(2), intent(out)                          :: lambdamu0
  integer, intent(in)                                             :: nmate, nmate_composites

  type(node), pointer                                             :: fi, cur_node, cur_node2
  type(nodeList), pointer                                         :: material, coeff, intVar, coeff_comp
  integer(kind=4)                                                 :: i,j, m, ind, nMat, ierror
  character(len=200)                                              :: err_msg,libname, matename, lawname
  integer                                                         :: arraySize1, alloc_stat
  integer(kind=8)                                                 :: nbZone, arraySize2 

  integer(kind=8), dimension(nmate)                               :: zones
  integer(kind=8), dimension(nmate)                               :: zones_tmp

  integer(kind=8)                                                 :: k,l,i_64
  logical                                                         :: material_presence
  integer(kind=4)                                                 :: mat_local_Index
  real(mytype),dimension(:,:),allocatable                         :: empty_array_real
  integer(kind=8),dimension(:,:),allocatable                      :: empty_array_zone
  
  real(mytype),dimension(:),allocatable                           :: array_real_1d
  integer(INT64),dimension(:),allocatable                         :: array_int64_1d

  ! Déclaration des variables propres à l'initialisation des composites
  !--------------------------------------------------------------------
  character(len=400)                                             :: file_name
  character(len=100)                                             :: function_name
  character(len=2)                                               :: iph_str
  integer                                                        :: FU
  integer                                                        :: nmate_hom,Ncoeff_lawC,NvarInt_lawC
  integer                                                        :: nZoneR
  integer(kind=INT64)                                            :: nZmax,ZoneCible
  integer,dimension(nmate-nmate_composites)                      :: NCoeff_list,NVarInt_list,NCoeff_comp_list 
  character(len=200),dimension(nmate-nmate_composites)           :: lawname_list,libname_list
       ! variables de boucles pour l'initialisation des composites
  integer                                                        :: icomp,iphi,imS,izS,imR,iph
  integer(INT64)                                                 :: izone
  integer                                                        :: Ind_ini,Ind_end
       ! variable de stockage des indices des données à transmettre/recevoir dans le pinceau
       ! source/cible
  integer                                                        :: imS0,izS0 ! indices dans le pinceau source
  integer                                                        :: imR0      ! indices dans le pinceau cible
       ! variables de rang de processus pour les communications MPI
  integer                                                        :: rksend0,rksend1 ! rang des processus sources (0 local, 1 global)
  integer                                                        :: rkrecv0 ,rkrecv_M, rkrecv ! rang des processus cibles (local)
  integer,dimension(nproc)                                       :: rkrecv1 ! rangs des processus cibles (liste globale)
  integer                                                        :: etiquette
  integer, dimension( MPI_STATUS_SIZE)                           :: statut
       ! variables de tailles dynamiques pour affectations multiples (supprimer si non retenu)
  integer,dimension(:),allocatable                               :: izR0_list
  logical,dimension(:),allocatable                               :: ZoneInit
  real(mytype),dimension(:,:),allocatable                        :: Coeff_tab2
  integer                                                        :: Nphases, test_laminate
  integer(INT64)                                                 :: NzoneP


  ! 
  allocate(Coeff_tab2(1,1)); deallocate(Coeff_tab2)  ! avoid gcc-warning

  ! allocation des tableaux vides utilises comme entree de get_vect_xml sur 
  ! les pinceaux n'ayant rien a lire (materiau non present sur le pinceau)
  allocate(empty_array_real(0,0))
  allocate(empty_array_zone(0,2)) ! alloue a la taille (0,2) pour coller a la forme
                                 ! du tableau MattotP0()%Zone, ayant 2 colonnes 
                                 ! et autant de lignes que de zones, dans le cas
                                 ! ou le materiau ne possede aucune zone sur le pinceau

  !Initialisation
  function_name="read_mat_composite"
  Ncoeff_lawC = 0

  ! Calcul du nombre de matériaux homogènes
  nmate_hom = nmate - nmate_composites

  ! zones_tmp: nombre maximal de zones pour chaque materiau et chaque proc
  zones_tmp=0
  do i=1,size(MattotP0)
    zones_tmp(MattotP0(i)%numM)=maxval(mattotP0(i)%zone(:,2))
  end do
  ! zones : nombre maximal de zones (dans tous les processus) pour chaque materiau
  call MPI_Allreduce(zones_tmp,zones,nmate, MPI_INTEGER8,MPI_MAX,amitex_mpi_comm,ierror)

  fi => parseFile(trim(file_mat))

  call GetNodeList(fi, material, "Reference_Material", 1,1)  
  ! lecture des proprietes du materiau de reference
  if(associated(material))then
    cur_node => item(material,0)
    call get_real_xml(cur_node,"Lambda0",lambdamu0(1),1,0)
    call get_real_xml(cur_node,"Mu0",lambdamu0(2),1,0)
  end if

    material => getElementsByTagName(fi,"Material")
  ! verification du nombre de materiaux
  ! le nombre de noeuds de material donne le nombre de matériaux  homogènes de la cellule
  ! le nombre total de matériaux doit être ce nombre additionné du nombre de matériaux composites
  if( nmate /= ( getLength(material)+nmate_composites) )then
    write(err_msg,fmt="(A,I0,A,I0)") &
            "Nombre de materiaux differents (read_mat_composite) :vtk (composite) :",&
            nmate,"; xml + ncomposites :",getLength(material)
    call amitex_abort(err_msg,2,0)
    call MPI_Barrier(amitex_mpi_comm,ierror)
  end if

  !====================================================================== LECTURE MODELES NON-LOCAUX  
  call read_nloc_models(fi)


  !---------------------------------------------------------------------
  !                                          BOUCLE  MATERIAUX HOMOGENES
  !---------------------------------------------------------------------
  do i=0,nmate-nmate_composites-1

    cur_node => item(material,i)
    ! numero de materiau du noeud i
    call get_int_xml(cur_node,"numM",nMat,2,0)

    ! nom de la librairie (chemin complet ou "amitex") associée (pas en minuscule)
    call get_str_xml_nolowercase(cur_node,"Lib",libname,1,0)

    ! Default library if Lib="" in xml file
    if (libname .eq. "") then
    if (AMITEXenv%path_stat .ne. 0) then
       call amitex_abort(&
       "Undefined environment variable AMITEX_FFTP",2,0) ! ERROR -> STOP
    else
       libname = AMITEXenv%path // "/libAmitex/src/materiaux/libUmatAmitex.so"
    end if
    end if

    ! nom de la loi associée 
     call get_str_xml(cur_node,"Law",lawname,1,0)

    ! material name
    matename=""
    call get_str_xml(cur_node,"Name",matename,-2)

    ! listes des noeuds de coefficients et variables internes
    coeff => getElementsByTagName(cur_node,"Coeff")
    coeff_comp => getElementsByTagName(cur_node,"Coeff_composite")
    intVar=> getElementsByTagName(cur_node,"IntVar")

    ! Mise en mémoire de certaines données pour affectation dans les matériaux composites ultérieurement
    NCoeff_list(i+1)      = getLength(coeff)
    NCoeff_comp_list(i+1) = getlength(coeff_comp)
    NVarInt_list(i+1)     = getLength(intVar)
    libname_list(i+1)     = libname
    lawname_list(i+1)     = lawname
    
    material_presence = .false. 
    mat_local_Index = -1
    do j=1,size(MattotP0)
        ! recherche du materiau dans MattotP0 pour savoir si il est sur le pinceau
        if(nMat== MattotP0(j)%numM)then
            material_presence = .true. 
            mat_local_Index = j
        end if
    end do

        ! bibliotheque .so et nom de loi   
        ! affectés si le materiau est sur le pinceau
    if (material_presence) MattotP0(mat_local_Index)%Libname = libname
    if (material_presence) MattotP0(mat_local_Index)%Lawname = lawname
    if (material_presence) MattotP0(mat_local_Index)%Matename = matename(1:100)
            
    nbZone=zones(nMat)  ! nombre maximal de zones du materiau

    ! lecture des coefficients
    if (material_presence) then
    ! la materiau est present sur le pinceaux : mise en forme des donnes necessaires
    ! a la recuperation des coefficients
    
        ! message d'avertissement si le materiau est une interphase a plusieurs zones
        if (MattotP0(mat_local_Index)%Interphase .and. (nbZone >1)) then
           write(err_msg,fmt="(3(A,I0),A,I4,A)") "Matériau d'interphase ", MattotP0(mat_local_Index)%numM,&
                                                   " découpé en ",nbZone," zones (read_mat_composite)"
           call amitex_abort(err_msg,-1)
        end if

        ! mise en forme des entrees et appel a get_vect_xml dans le cas ou le materiau est present
        arraySize1= getLength(coeff)           ! nombre de coefficients
        MattotP0(mat_local_Index)%NCoeff=arraySize1
        arraySize2= size(MattotP0(mat_local_Index)%zone(:,1),kind=8) ! nombre de zones du materiau sur le processus
        ! allocation du tableau coeff
        allocate(MattotP0(mat_local_Index)%coeff( arraySize1, arraySize2),stat=alloc_stat )
        if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_mat_composite)",2)
        ! lecture des valeurs des coefficients
        call get_vect_xml(nMat,coeff,arraySize1,arraySize2,MattotP0(mat_local_Index)%coeff,nbZone,arraySize2,&
                          MattotP0(mat_local_Index)%zone,"Coeff",MattotP0(mat_local_Index)%pos,test_composite=.true.)    
                          
        ! allocation du tableau Coeffname
        allocate(MattotP0(mat_local_Index)%Coeffname(arraySize1),stat=alloc_stat )
        if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant MattotP0%Coeffname (read_mat)",2)

        ! Read Coeff names :
        ! REMARK : below, no check on ind, already done in get_vect_xml
        do m=0,arraySize1-1
           cur_node2 => item(coeff,m)
           ind = -1
           matename=""
           call get_str_xml(cur_node2,"Name",matename,-2)
           call get_int_xml(cur_node2,"Index",ind,1)
           MattotP0(mat_local_Index)%Coeffname(ind) = matename(1:100)
        end do                     
        
        ! verification de la lecture correcte des coefficients (si le materiau est present sur le pinceau)
        do l=1,arraySize2
          do k=1,arraySize1
            if(MattotP0(mat_local_Index)%coeff(k,l)/=MattotP0(mat_local_Index)%coeff(k,l))then
              write(err_msg,fmt="(3(A,I0),A)") "Coefficient:",k, " zone:",MattotP0(mat_local_Index)%zone(l,2),&
                " materiau:",nMat," sans valeur (read_mat_composite)"
              call amitex_abort(err_msg,1)
            end if
          end do
        end do  
    else
    ! materiau non present sur le pinceau : on envoie des tableaux vides à get_vect_xml
        arraysize1 = 0
        arraysize2 = 0
        ! appel a get_vect_xml pour assurer que les pinceaux ne possedant pas le materiau 
        ! soit associés a la lecture parallèle des fichiers binaires
        !call get_vect_xml(nMat,coeff,arraySize1,arraysize2,empty_array_real,nbZone,&
        !                  arraySize2,empty_array_zone,"Coeff",MattotP0(mat_local_Index)%pos,test_composite=.true.)
        allocate(array_int64_1d(0))
        call get_vect_xml(nMat,coeff,arraySize1,arraysize2,empty_array_real,nbZone,&
                          arraySize2,empty_array_zone,"Coeff",array_int64_1d,test_composite=.true.)
        deallocate(array_int64_1d)
    end if

        
    if (algo_param0%Jacobian_type_laminate == "elastic" .or. algo_param0%Jacobian_type_reuss == "elastic") then
       ! lecture des coefficients composites (utilisés par le modèle multi-couche composite)
       if (material_presence) then
           arraySize1= getLength(coeff_comp)           ! nombre de coefficients
           arraySize2= size(MattotP0(mat_local_Index)%zone(:,1),kind=8) ! nombre de zones du materiau sur le processus
                                                         ! pas de bug en presence de zones interphase 
                                                         ! elles possedent bien des coefficients composites a lire
                                                         ! et stocker
           allocate(MattotP0(mat_local_Index)%Coeff_comp( arraySize1, arraySize2),stat=alloc_stat )
           if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant &
                                                & (read_mat_composite)",2)
           ! lecture des valeurs des coefficients
           call get_vect_xml(nMat,coeff_comp,arraySize1,arraySize2,MattotP0(mat_local_Index)%Coeff_comp,&
                            nbZone,arraySize2,MattotP0(mat_local_Index)%zone,"Coeff_composite",&
                            MattotP0(mat_local_Index)%pos,test_composite=.true.)
        else 
        ! materiau non present sur le pinceau : on envoie des tableaux vides à get_vect_xml
            arraysize1 = 0
            arraysize2 = 0
            ! appel a get_vect_xml pour assurer que les pinceaux ne possedant pas le materiau 
            ! soit associés a la lecture parallèle des fichiers binaires
            !call get_vect_xml(nMat,coeff_comp,arraySize1,arraysize2,empty_array_real,nbZone,arraySize2,&
            !                  empty_array_zone,"Coeff_composite",MattotP0(mat_local_Index)%pos,test_composite=.true.)
            allocate(array_int64_1d(0))          
            call get_vect_xml(nMat,coeff_comp,arraySize1,arraysize2,empty_array_real,nbZone,arraySize2,&
                              empty_array_zone,"Coeff_composite",array_int64_1d,test_composite=.true.)
            deallocate(array_int64_1d)
            
            ! verification de la lecture correcte des coefficients composites (si le materiau est present sur le pinceau)
            do l=1,arraySize2
              do k=1,arraySize1
                if(MattotP0(mat_local_Index)%coeff_comp(k,l)/=MattotP0(mat_local_Index)%coeff_comp(k,l))then
                  write(err_msg,fmt="(3(A,I0),A)") "Coefficient composite :",k, " zone:",MattotP0(mat_local_Index)%zone(l,2),&
                    " materiau:",nMat," sans valeur (read_mat_composite)"
                  call amitex_abort(err_msg,1)
                end if
              end do
            end do
        end if
    else
       ! Si la matrice tangente demandé pour l'inégration du modèle multi-couche n'est pas la matrice d'élasticité équivalente :
       !  les coeffs comp sont initialisés à 0 (pas optimal...) car non utilisés 
       if (material_presence) then
           arraySize1 = 2
           arraySize2= size(MattotP0(mat_local_Index)%zone(:,1),kind=8) ! nombre de zones du materiau sur le processus
                                                         ! pas de bug en presence de zones interphase 
                                                         ! elles possedent bien des coefficients composites a lire
                                                         ! et stocker
           allocate(MattotP0(mat_local_Index)%Coeff_comp( arraySize1, arraySize2),stat=alloc_stat )
           if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant &
                                               & (read_mat_composite)",2)
           MattotP0(mat_local_Index)%Coeff_comp = 0.
       end if
    end if

    
    ! lecture des variables internes
        arraySize1=getLength(intVar)        ! nombre de variables internes
        if (material_presence) MattotP0(mat_local_Index)%NVarInt=arraySize1
        if(arraySize1>0)then
        ! variables internes a lire 
            if (material_presence) then
            ! lecture des variables internes dans le cas ou le materiau est sur le pinceau
              arraySize2=size(MattotP0(mat_local_Index)%pos,kind=8)   ! nombre de voxels du materiau dans le processus
              if (MattotP0(mat_local_Index)%Interphase) then
                 ! Si le matériau est un matériau d'interphase, alors il a un tableau pos vide -> on crée tout de même un tableau pour contenir 
                 ! les VarInt pour l'initialisation des voxels composites
                 ! arraySize2 = 1
                 arraySize2=size(MattotP0(mat_local_Index)%zone(:,1))  
              end if
              if (allocated(MattotP0(mat_local_Index)%Zones_interphase) .and. (nrank == 0)) then
                 ! Si des zones interphases sont presentes et que l'on se trouve sur le pinceau 0
                 ! il faut agrandir le tableau %VarInt pour contenir les variables internes des zones interphase
                 ! pour l'initialisation des voxels composites
                 arraySize2 = arraySize2 + size(MattotP0(mat_local_Index)%Zones_interphase)
              end if
              l=size(MattotP0(mat_local_Index)%zone(:,1),kind=8)      ! nombre de zones 
                                                       ! ATTENTION : pas de bug ici en presence de zones interphases
                                                       ! elles sont comprises dans Zone(:,1) et doivent etre comptees
                                                       ! pour pouvoir les lire
              allocate(MattotP0(mat_local_Index)%VarInt(arraySize1,arraySize2),stat=alloc_stat)
              if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_mat_composite)",2)
              allocate(MattotP0(mat_local_Index)%VarInt0(arraySize1,arraySize2),stat=alloc_stat)
              if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_mat_composite)",2)
              MattotP0(mat_local_Index)%VarInt0=-1
              ! lecture des variables internes
              call get_vect_xml(nMat,intVar,arraySize1,arraySize2,MattotP0(mat_local_Index)%VarInt0,nbZone,l,&
                                MattotP0(mat_local_Index)%zone,"IntVar",MattotP0(mat_local_Index)%pos,test_composite=.true.)

              ! allocation du tableau Varintname
              allocate(MattotP0(mat_local_Index)%Varintname(arraySize1),stat=alloc_stat )
              if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant MattotP0%Varintname (read_mat)",2)

              ! Read Varint names :
              ! REMARK : below, no check on ind, already done in get_vect_xml
              do m=0,arraySize1-1
                 cur_node2 => item(intVar,m)
                 ind = -1
                 matename=""
                 call get_str_xml(cur_node2,"Name",matename,-2)
                 call get_int_xml(cur_node2,"Index",ind,1)
                 MattotP0(mat_local_Index)%Varintname(ind) = matename(1:100)
              end do

              do l=1,size(MattotP0(mat_local_Index)%zone(:,1)) 
                ! ATTENTION
                ! pas de bug ici vis a vis de la presence de zones interphases. Les lignes correspondantes de Zone(:,1)
                ! contiennent les numeros de colonne de %VarInt où stocker les variables internes de ces zones --> OK 
                do k=1,arraySize1
                  if( MattotP0(mat_local_Index)%VarInt0(k,MattotP0(mat_local_Index)%zone(l,1)) /= &
                      MattotP0(mat_local_Index)%VarInt0(k,MattotP0(mat_local_Index)%zone(l,1)))then
                    write(err_msg,fmt="(3(A,I0),A)") "Variables internes:",k, " zone:",MattotP0(mat_local_Index)%zone(l,2),&
                    " materiau:",nMat," sans valeur (read_mat_composite)"
                    call amitex_abort(err_msg,1)
                  end if
                end do
              end do

            else
            ! materiau non present sur le pinceau : on envoie des tableaux vides à get_vect_xml
                arraysize1 = 0
                arraysize2 = 0
                ! appel a get_vect_xml pour assurer que les pinceaux ne possedant pas le materiau 
                ! soit associés a la lecture parallèle des fichiers binaires
                !call get_vect_xml(nMat,intVar,arraySize1,arraysize2,empty_array_real,nbZone,&
                !                  arraySize2,empty_array_zone,"IntVar",MattotP0(mat_local_Index)%pos,test_composite=.true.)  
                allocate(array_int64_1d(0))
                call get_vect_xml(nMat,intVar,arraySize1,arraysize2,empty_array_real,nbZone,&
                                  arraySize2,empty_array_zone,"IntVar",array_int64_1d,test_composite=.true.) 
                deallocate(array_int64_1d) 
            end if
        else
        ! pas de variables internes a lire : allocation de tableaux vides
          if (material_presence) then
          allocate(MattotP0(mat_local_Index)%VarInt(0,0),stat=alloc_stat)
          if(alloc_stat /=0)&
             call amitex_abort("Espace memoire disponible insuffisant MattotP0%VarInt(0,0)(read_mat_composite)",2)
          allocate(MattotP0(mat_local_Index)%VarInt0(0,0),stat=alloc_stat)
          if(alloc_stat /=0) &
             call amitex_abort("Espace memoire disponible insuffisant MattotP0%VarInt0(0,0)(read_mat_composite)",2)
          end if
        end if
        
        !LECTURE DES INDICES DES COEFF ET VAR INT DES MODELES NON NOCAUX  
        if (algo_param0%Nloc) then
              cur_node => item(material,i)
              call read_nloc_indices(cur_node,nMat,intVar)
        end if

 end do

  ! désallocation des tableaux vides utilises comme entree de get_vect_xml sur 
  ! les pinceaux n'ayant rien a lire.
  deallocate(empty_array_real)
  deallocate(empty_array_zone)
  
  call destroy(fi)

  call MPI_BARRIER(amitex_mpi_comm,ierror)
  
  !---------------------------------------------------------------------
  !                                                 MATERIAUX COMPOSITES   
  !---------------------------------------------------------------------
  
  do icomp = 1,nmate_composites

     ! Boucle sur les matériaux composites de la cellule (icomp)
     if (nrank == 0) write(OUTPUT_UNIT,'(A,I4,A)') "Initialisation matériau composite : ",icomp," en cours ........"

     rkrecv_M = -1
     imR0     = -1

     ! Recherche des pinceaux possédant ce matériau
     do imR=1,size(MattotP0)
        ! Boucle sur les matériaux contenus dans le pinceau : recherche du matériau à affecter
        if (MattotP0(imR)%numM == nmate_hom + icomp) then
           rkrecv_M = nrank
           imR0 = imR
           exit ! info trouvée on sort de la boucle
        end if
     end do

     ! on evalue le nbre de phases du materiau composite icomp
     nPhases = size(MatComposite(icomp)%num_Phases)
     
     ! on evalue si le materiau est de type "laminate"
     test_laminate = 0
     if (trim(MatComposite(icomp)%Lawname) == "laminate") then
        test_laminate = 1
     end if 
     

     !----------------------------------------------------------------------------------
     ! rkrecv_M /= -1 sur les pinceaux ou le matériau composite à initialiser est présent
     ! imR0 est l'indice local du matériau à initialiser dans le tableau MattotP0
     !----------------------------------------------------------------------------------

     
     ! si le matériau a été repéré dans le pinceau son indice est imR0 et on procède aux affectations
     if (rkrecv_M >= 0) then

        ! Affectations des champs globaux du matériau composite icomp
        !-------------------------------------------------------------
        !  Allocations

        ! Nombre de coefficients de la loi de comportement de la phase i
        allocate(MattotP0(imR0)%NCoeff_comp(MattotP0(imR0)%Nphase),stat=alloc_stat)
        if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant NCoeff_comp (read_mat_composite)",2)

        ! Nombre de variables internes de la phase i

        allocate(MattotP0(imR0)%NvarInt_comp(MattotP0(imR0)%Nphase),stat=alloc_stat)
        if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant NVarInt_comp (read_mat_composite)",2)

        ! Librairie contenant la loi de comportement de la phase i
        allocate(MattotP0(imR0)%Libname_comp(MattotP0(imR0)%Nphase),stat=alloc_stat)
        if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant Libname_comp (read_mat_composite)",2)

        ! Loi de comportement de la phase i
        allocate(MattotP0(imR0)%Lawname_comp(MattotP0(imR0)%Nphase),stat=alloc_stat)
        if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant Lawname_comp  (read_mat_composite)",2)

        ! Nom de la loi de comportement composite
        MattotP0(imR0)%lawname = MatComposite(icomp)%Lawname

        !   Affectactions des grandeurs par phase
        do i=1,MattotP0(imR0)%Nphase
           ! Nombre de coefficients de la loi de comportement de la phase i
           ! + nombre de coefficients de la phase propres au modèle composite
           MattotP0(imR0)%NCoeff_comp(i) = NCoeff_list(MattotP0(imR0)%numM_composite(i))

           ! Nombre de variables internes de la phase i
           MattotP0(imR0)%NvarInt_comp(i) = NVarInt_list(MattotP0(imR0)%numM_composite(i)) 

           ! Librairie contenant la loi de comportement de la phase i
           MattotP0(imR0)%Libname_comp(i) = libname_list(MattotP0(imR0)%numM_composite(i)) 

           ! Loi de comportement de la phase i
           MattotP0(imR0)%Lawname_comp(i) = lawname_list(MattotP0(imR0)%numM_composite(i)) 
        end do

        ! Nombre de coefficients/varInt + allocations des tableaux locaux
        !--------------------------------------------------------------------
        !    nombre de coefficients/variables internes specifiquement liées à la loi composite
        !          (n'inclue pas les coeff. elastiques equivalent pour l'algo NR modifie)
        ! On affecte aussi les parametres d'algo pour pouvoir unifier les calculs de comportement

        if (trim(MattotP0(imR0)%lawname) == "laminate") then
           ! voir subroutine umatlaminate
           NCoeff_lawC  = 6 + 2*MattotP0(imR0)%Nphase 
           NVarInt_lawC = 10*MattotP0(imR0)%Nphase
           MattotP0(imR0)%acc_CV_composite = algo_param0%acc_CV_laminate
           MattotP0(imR0)%tol_criteq_composite = algo_param0%tol_criteq_laminate
           MattotP0(imR0)%Init_composite = algo_param0%Init_laminate
           MattotP0(imR0)%Nmax_subdivision_composite = algo_param0%Nmax_subdivision_laminate
           MattotP0(imR0)%Jacobian_type_composite = algo_param0%Jacobian_type_laminate
           MattotP0(imR0)%Perturbation_composite = algo_param0%Perturbation_laminate

        else if (trim(MattotP0(imR0)%lawname) == "voigt")  then
           ! voir subroutine umatvoigt
           NCoeff_lawC  = 2*MattotP0(imR0)%Nphase 
           NVarInt_lawC = 7*MattotP0(imR0)%Nphase 
           ! Pour voigt, les valeurs qui suivent ne servent a rien car la convergence est immediate
           ! On remplit des valeurs par defaut pour eviter des bugs 
           MattotP0(imR0)%acc_CV_composite = .false.
           MattotP0(imR0)%tol_criteq_composite = 0.01
           MattotP0(imR0)%Init_composite = 'default'
           MattotP0(imR0)%Nmax_subdivision_composite = 0
           MattotP0(imR0)%Jacobian_type_composite = ""
           MattotP0(imR0)%Perturbation_composite = 0.

        else if (trim(MattotP0(imR0)%lawname) == "reuss")  then
           ! voir subroutine umatreuss
           NCoeff_lawC  = 2*MattotP0(imR0)%Nphase 
           NVarInt_lawC = 13*MattotP0(imR0)%Nphase 
           MattotP0(imR0)%acc_CV_composite = algo_param0%acc_CV_reuss
           MattotP0(imR0)%tol_criteq_composite = algo_param0%tol_criteq_reuss
           MattotP0(imR0)%Init_composite = algo_param0%Init_reuss
           MattotP0(imR0)%Nmax_subdivision_composite = algo_param0%Nmax_subdivision_reuss
           MattotP0(imR0)%Jacobian_type_composite = algo_param0%Jacobian_type_reuss
           MattotP0(imR0)%Perturbation_composite = algo_param0%Perturbation_reuss
        else 
           write(err_msg,fmt="(3(A,I0),A)") "Matériau composite ",icomp,&
 " loi de comportement composite inconnue. Lois implémentées : 'laminate', 'voigt' ou 'reuss' (read_mat_composite)"
           call amitex_abort(err_msg,1)
        end if

        ! Allocation du tableau Coeff 
        arraysize1 = NCoeff_lawC + sum( MattotP0(imR0)%NCoeff_comp)
        if (trim(MattotP0(imR0)%lawname) == "laminate") arraySize1 = arraySize1 + 2*MattotP0(imR0)%Nphase 
        if (trim(MattotP0(imR0)%lawname) == "reuss")    arraySize1 = arraySize1 + 2*MattotP0(imR0)%Nphase 
        arraysize2 = size(MattotP0(imR0)%Zone(:,1),kind=8)
        allocate(MattotP0(imR0)%coeff( arraySize1, arraySize2),stat=alloc_stat )
        if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant &
                                            & coeff mat composite (read_mat_composite)",2)

        MattotP0(imR0)%coeff = 0_mytype

        ! Affectation de NCoeff
        MattotP0(imR0)%NCoeff  = arraysize1

        if (allocated(Coeff_tab2)) deallocate(Coeff_tab2)

        ! Affectation des coefficients communs aux différentes loi composites 
        allocate(Coeff_tab2(MattotP0(imR0)%Nphase + 1,arraySize2),stat = alloc_stat)
        if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant &
                                            & coeff_tab2 Coeff communs/phases (read_mat_composite)",2)
        Coeff_tab2 = 0._mytype

        ! on stocke les coefficient dans un tableau temporaire Coeff_tab2 pour 
        ! pouvoir tout allouer d'un coup (évite une longue boucle)
        Coeff_tab2(1,1) = MattotP0(imR0)%Nphase 
        Coeff_tab2(2:1+MattotP0(imR0)%Nphase,1) = MattotP0(imR0)%NCoeff_comp
        !Coeff_tab2(:,:) = spread(Coeff_tab2(:,1),2,arraySize2)   ! bug intel-oneapi 2021.1 (arraysize2 trop grand?)
        call spread_perso(Coeff_tab2,Coeff_tab2(:,1))  

        ! on affecte les coefficient nPhase ainsi que les nCoeff_i
        MattotP0(imR0)%Coeff(1:1+MattotP0(imR0)%Nphase,:) = Coeff_tab2(:,:)  

        !--------------------------------------------- LECTURES FICHIER DES DONNEES COMPOSITES
        ! On affecte les fractions volumiques        
        ! MattotP0(imR0)%Coeff(2+MattotP0(imR0)%Nphase:2*MattotP0(imR0)%Nphase,:) = &
        !               transpose(MatComposite(icomp)%Fv(MattotP0(imR0)%Zone(:,2),1:MattotP0(imR0)%Nphase-1))
        allocate(array_real_1d(size(MattotP0(imR0)%Zone(:,2))))
        nZoneP = int(size(MattotP0(imR0)%Zone(:,2)),INT64)
        
        do iph = 1, Nphases-1
          write(iph_str,fmt="(I0)") iph 
          file_name=trim(MatComposite(icomp)%dir)//'/fv'// trim(iph_str) //'.bin'
          FU = 724212395+iph
          call check_file(file_name,function_name)
          call  getBinCoeffFromFile(file_name, MattotP0(imR0)%Zone(:,2), nZoneP, array_real_1d, FU )
          MattotP0(imR0)%Coeff(1+Nphases + iph,:) = array_real_1d
        end do

        ! On affecte les coefficients specifiques a la loi laminate
        !if (trim(MattotP0(imR)%lawname) == "laminate") then
        if (test_laminate == 1) then
           ! Composantes de la normale à l'interface
           ! MattotP0(imR)%Coeff(2*MattotP0(imR)%Nphase+1:2*MattotP0(imR)%Nphase+3,:) = &
           !            transpose(MatComposite(icomp)%Normale(MattotP0(imR)%Zone(:,2),:))
           do iph = 1, 3
             write(iph_str,fmt="(I0)") iph 
             if(iph==1) file_name=trim(MatComposite(icomp)%dir)//'/N12x.bin'
             if(iph==2) file_name=trim(MatComposite(icomp)%dir)//'/N12y.bin'
             if(iph==3) file_name=trim(MatComposite(icomp)%dir)//'/N12z.bin'
             FU = 22015352 + iph
             call check_file(file_name,function_name)
             call  getBinCoeffFromFile(file_name, MattotP0(imR0)%Zone(:,2), nZoneP, array_real_1d, FU )
             MattotP0(imR)%Coeff(2*Nphases+iph,:) = array_real_1d
           end do
           
           ! Composantes de la direction tangente à l'interface
           ! MattotP0(imR)%Coeff(2*MattotP0(imR)%Nphase+4:2*MattotP0(imR)%Nphase+6,:) = &
           !            transpose(MatComposite(icomp)%Direction(MattotP0(imR)%Zone(:,2),:))
           do iph = 1, 3
             write(iph_str,fmt="(I0)") iph 
             if(iph==1) file_name=trim(MatComposite(icomp)%dir)//'/Tx.bin'
             if(iph==2) file_name=trim(MatComposite(icomp)%dir)//'/Ty.bin'
             if(iph==3) file_name=trim(MatComposite(icomp)%dir)//'/Tz.bin'
             FU = 22015352 + 3 + iph
             call check_file(file_name,function_name)
             call  getBinCoeffFromFile(file_name, MattotP0(imR0)%Zone(:,2), nZoneP, array_real_1d, FU )
             MattotP0(imR)%Coeff(2*Nphases+3+iph,:) = array_real_1d
           end do                     
        end if
        
        deallocate(array_real_1d)
        !--------------------------------------------- FIN LECTURES FICHIER DES DONNEES COMPOSITES

        if (allocated(Coeff_tab2)) deallocate(Coeff_tab2)

        ! Allocation des tableaux VarInt et VarInt0 
        arraysize1 = NVarInt_lawC + sum( MattotP0(imR0)%NVarInt_comp)
        arraysize2 = size(MattotP0(imR0)%Zone(:,1),kind=8)
        allocate(MattotP0(imR0)%VarInt( arraySize1, arraySize2),stat=alloc_stat )
        write(err_msg,fmt="(A,2(I0),A,I0,A,I0,A)") "Espace memoire disponible insuffisant VarInt, taille ",&
                      arraysize1,arraysize2," NVarInt_lawC = ",NVarInt_lawC,&
                      " rkrecv_M : ",rkrecv_M,"(read_mat_composite)"
        if(alloc_stat /=0) call amitex_abort(err_msg,2)
        allocate(MattotP0(imR0)%VarInt0( arraySize1, arraySize2),stat=alloc_stat )
        if(alloc_stat /=0) &
               call amitex_abort("Espace memoire disponible insuffisant VarInt (read_mat_composite)" ,2)

        ! Affectation de NVarint
        MattotP0(imR0)%NVarint = arraysize1

        MattotP0(imR0)%VarInt  = 0.

        ! Affectation des variables internes communes aux différentes loi composites  (nvarint(i))
        allocate(Coeff_tab2(MattotP0(imR0)%Nphase,arraySize2),stat = alloc_stat)
        if(alloc_stat /=0)  call amitex_abort("Espace memoire disponible insuffisant &
                                 & coeff_tab2 affectation VarInt communes/phases (read_mat_composite)",2)

        ! on stocke les coefficient dans un tableau temporaire Coeff_tab2 
        ! pour pouvoir tout allouer d'un coup (évite une longue boucle)
        Coeff_tab2(:,1) = MattotP0(imR0)%NvarInt_comp           
                                                               
        !Coeff_tab2(:,:) = spread(Coeff_tab2(:,1),2,arraySize2) ! bug intel-oneapi 2021.1 (arraysize2 trop grand?)
        call spread_perso(Coeff_tab2,Coeff_tab2(:,1))  
                                                               
        MattotP0(imR0)%VarInt(1:MattotP0(imR0)%Nphase,:) = Coeff_tab2

        MattotP0(imR0)%VarInt0 = MattotP0(imR0)%VarInt

        if (allocated(Coeff_tab2)) deallocate(Coeff_tab2)

     else      ! Materiau composite absent du pinceau MAIS
               !               TOUS LES PINCEAUX doivent appeler getbinCoeffFromFile
               !----------------------------------------------------------------------------------------
        allocate(array_real_1d(0))
        allocate(array_int64_1d(0))
        nZoneP = 0_INT64
      
        do iph = 1, Nphases-1
          write(iph_str,fmt="(I0)") iph 
          file_name=trim(MatComposite(icomp)%dir)//'/fv'//trim(iph_str)//'.bin'
          FU = 724212395+iph
          call check_file(file_name,function_name)
          call  getBinCoeffFromFile(file_name, array_int64_1d, nZoneP, array_real_1d, FU )
        end do
       
        ! Lecture des donnees 'laminate'
        if (test_laminate>0) then
        !lecture des composantes du vecteur normal
        do iph = 1, 3
          write(iph_str,fmt="(I0)") iph 
          if(iph==1) file_name=trim(MatComposite(icomp)%dir)//'/N12x.bin'
          if(iph==2) file_name=trim(MatComposite(icomp)%dir)//'/N12y.bin'
          if(iph==3) file_name=trim(MatComposite(icomp)%dir)//'/N12z.bin'
          FU = 22015352+iph
          call check_file(file_name,function_name)
          call  getBinCoeffFromFile(file_name, array_int64_1d, nZoneP, array_real_1d, FU )
        end do
        
        !lecture des composantes du vecteur tangentiel
        do iph = 1, 3
          write(iph_str,fmt="(I0)") iph 
          if(iph==1) file_name=trim(MatComposite(icomp)%dir)//'/Tx.bin'
          if(iph==2) file_name=trim(MatComposite(icomp)%dir)//'/Ty.bin'
          if(iph==3) file_name=trim(MatComposite(icomp)%dir)//'/Tz.bin'
          FU = 22015352 + 3 + iph
          call check_file(file_name,function_name)
          call  getBinCoeffFromFile(file_name, array_int64_1d, nZoneP, array_real_1d, FU )
        end do   
        end if                  

        deallocate(array_real_1d)
        deallocate(array_int64_1d)
        
     end if    ! FIN TEST  presence materiaux composite sur pinceau (if (rkrecv_M >= 0))
     
     
     ! Affectations des champs locaux du matériau composite icomp
     !-------------------------------------------------------------
     !do iphi = 1,size(MatComposite(icomp)%num_Phases)
     do iphi = 1,Nphases
        ! Boucle sur les phases du matériaux composite icomp (iphi)

        ! Lecture Fichier des tableaux de zones de chaque phase -> array_real_1d 
        !                             attention, getBinCoeffFromFile renvoie des reels
        !   TODO : peut-on envisager un lecture par blocs pour limiter l'empreinte memoire?  
        nZoneP = MatComposite(icomp)%Nzones_global
        allocate(array_int64_1d(nZoneP))
        allocate(array_real_1d(nZoneP))
        array_int64_1d = (/ (i_64,i_64=1,nZoneP) /)
        write(file_name,"(A,A,I0,A)") trim(MatComposite(icomp)%dir),'/zone',iphi,'.bin'
        call check_file(file_name,function_name)
        call getBinCoeffFromFile(file_name, array_int64_1d, nZoneP, array_real_1d, FU )
        
        ! Allocation d'un tableau ZoneInit permettant de retracer les voxels composites ayant été initialisées
        !nZmax = maxval(MatComposite(icomp)%Zone(:,iphi+1))
        nZmax = int(maxval(array_real_1d),INT64)
        allocate(ZoneInit(nZmax),stat=alloc_stat)
        if(alloc_stat /=0) &
                call amitex_abort("Espace memoire disponible insuffisant ZoneInit (read_mat_composite)",2)
        ZoneInit = .false. ! gcc-warning accepte 

        do izone = 1,MatComposite(icomp)%Nzones_global
!        do izone = 1,size(MatComposite(icomp)%Zone(:,1))
           ! Boucle sur les zones du matériau composite icomp (izone)

           ! Test : si le voxel est déjà initialisé, on passe au suivant
           !ZoneCible = MatComposite(icomp)%Zone(izone,iphi+1)
           ZoneCible = int(array_real_1d(izone),INT64)
           if (ZoneInit(ZoneCible)) cycle

           ZoneInit(ZoneCible) = .true. ! le statut du voxel passe à initialisé 
 
           !! Repérage d'un pinceau possédant la combinaison Matériau/Zone homogène à transmettre
           !--------------------------------------------------------------------------------
           rksend0 = -1 ! mise à -1 de rksend sur tous les pinceaux pour distinguer la non détection du rang 0
           izS0 = -1
           imS0 = -1

           do imS = 1,size(MattotP0)
              ! Boucle sur les matériaux contenus dans le pinceau (imS)
              
              if (MattotP0(imS)%numM == MatComposite(icomp)%num_Phases(iphi)) then
                 ! la phase est présente dans le pinceau : on cherche maintenant la zone

                 do izS =1,size(MattotP0(imS)%Zone(:,1))
                    ! Boucle sur les zones du matériau imS dans le pinceau courant (izS)

!                    if (MattotP0(imS)%Zone(izS,2) == MatComposite(icomp)%Zone(izone,iPhi+1)) then
                    if (MattotP0(imS)%Zone(izS,2) == int(array_real_1d(izone),INT64) ) then
                       ! la zone recherchée du matériau imS dans le pinceau est également présente :
                       ! récupération de la localisation des informations (rang du processus et indices dans les tableaux)
                       imS0 = imS
                       izS0 = izS
                       rksend0 = nrank

                       exit ! info trouvée, sortie de la boucle
                    end if

                    ! fin boucle sur les zones du matériau imS (izS)
                 end do
              end if
              if (rksend0 /= -1) exit ! info trouvée, on sort de la boucle

              ! fin boucle sur les matériaux contenus dans le pinceau (imS)
           end do

           !-------------------------------------------------------------------------------------------------
           ! rksend0 /= 0 sur les pinceaux qui possèdent la combinaison Matériau/Zone Homogène à affecter
           ! imS0         est l'indice local du matériau homogène contenant l'information sur ces pinceaux dans 
           !              le tableau MattotP0
           ! izS0         est l'indice local de la zone du matériau homogène imS0 recherchée dans les tableaux
           !              de la structure MattotP0%Zone/Coeff/VarInt (premier indice pour Zone, second pour 
           !              Coeff, VarInt, VarInt0)
           !-------------------------------------------------------------------------------------------------

           ! Récupération d'un seul rang de processus source pour communication des données ultérieure
           call MPI_ALLreduce(rksend0,rksend1,1,MPI_INTEGER,MPI_MAX,amitex_mpi_comm,ierror)
           
           if (rksend1 == -1) then
!              write(err_msg,fmt="(A,I4,A,I4,A)") &
!                   "Combinaison matériau zone : ",MatComposite(icomp)%num_Phases(iphi),"/"&
!                   ,MatComposite(icomp)%Zone(izone,iPhi+1)," non trouvée pour &
!                   &initialisation du materiau composite (read_mat_composite)"
              write(err_msg,fmt="(A,I4,A,I4,A)") &
                   "Combinaison matériau zone : ",MatComposite(icomp)%num_Phases(iphi),"/"&
                   ,int(array_real_1d(izone),INT64)," non trouvée pour &
                   &initialisation du materiau composite (read_mat_composite)"
              call amitex_abort(err_msg,2)
           end if

           !! Repérage des matériaux possédant la combinaison Matériau/Zone homogène à affecter
           !  à un matériau composite, et n'ayant pas l'info dans leur pinceau
           !-----------------------------------------------------------------------------------     
           rkrecv1(:) = -1
           rkrecv0    = -1
           rkrecv     = -1
           
           if (rkrecv_M >= 0) then
              ! la matériau homogène a été trouvé précédemment dans le pinceau, il faut maintenant chercher si la zone 
              ! du matériau homogène à rechercher s'y trouve également

              ! on compte combien de voxels composites composés de la combinaison
              ! Matériau_homogène/Zone actuellement recherchée sont à initialiser pour cette phase
              nZoneR = count(MattotP0(imR0)%numZ_composite(:,iphi) == ZoneCible) 
              if (nZoneR > 0) then
                 ! dans le cas ou ce nombre est non nul : on conserve les indices locaux de tous les voxels
                 ! à initialiser
                 allocate(izR0_list(nZoneR),stat=alloc_stat)
                 if(alloc_stat /=0)  &
                  call amitex_abort("Espace memoire disponible insuffisant izR0_list (read_mat_composite)",1)
                 izR0_list = pack( (/ (i,i=1,size(MattotP0(imR0)%numZ_composite(:,iphi) ) )  /) ,&
                             MattotP0(imR0)%numZ_composite(:,iphi) == ZoneCible) 
                             ! liste des indices des voxels ayant la phase recherchée (Matériau/Zone)
                 rkrecv = nrank ! le pinceau est alors signalé comme devant recevoir les données locales de cette combinaison (Matériau/Zone) pour le matériau composite icomp
              end if

           end if

           !-------------------------------------------------------------------------------------------------
           ! rkrecv /= -1  sur les pinceaux ou la combinaison Materiau/Zone composite à initialiser est présente (doit potentiellement recevoir l'info)
           ! rkrecv0 /= -1 sur les pinceaux ou la combinaison Materiau/Zone composite à initialiser est présente (doit effectivement recevoir l'info)
           !               et la combinaison Matériau/Zone homogène à affecter absente
           !-------------------------------------------------------------------------------------------------

           if ((rkrecv /= -1) .and. (rksend0 == -1)) rkrecv0 = rkrecv 

           ! Récupération de tous les rangs des processus qui doivent recevoir l'info 
           call MPI_ALLgather(rkrecv0,1,MPI_INTEGER,rkrecv1,1,MPI_INTEGER,amitex_mpi_comm,ierror)  ! rkrecv1(i) = rkrecv0 du pinceau i-1

           !-----------------------------------------------------------------------------------------
           ! AFFECTATIONS LOCALES
           !-----------------------
           !       Si le pinceau possède le matériau à affecter et les données recherchées 
           !       (ie rkrecv0 /= -1 et rksend0 /= -1), pas de communications MPI
           !       Communications point à point MPI pour tous les pinceau ayant le matériau à affecter
           !       mais pas les données recherchées
           !-----------------------------------------------------------------------------------------
           
           if (rkrecv /= -1) then
              ! le matériau et des voxels composites à initialiser sont présents sur le pinceau 

              
              if (rksend0 /= -1) then
                 ! le matériau homogène et la zone à affecter sont présents sur le pinceau

                 if (allocated(Coeff_tab2)) deallocate(Coeff_tab2)

                 !------------------------------------ Coefficients -----------------------------------
                 ! allocation du tableau temporaire Coeff_tab2 qui contiendra l'information à affecter pour les coefficients
                 ! pour 'laminate' et 'reuss' on doit rajouter deux coeff lamda,mu utiles pour l'algo NR modifie
                 arraySize1 = MattotP0(imR0)%NCoeff_comp(iphi)
                 if  (trim(MattotP0(imR0)%lawname) == "laminate") arraySize1 = arraySize1 + 2
                 if  (trim(MattotP0(imR0)%lawname) == "reuss") arraySize1 = arraySize1 + 2
                  arraySize2 = nZoneR
                 allocate(Coeff_tab2(arraySize1,arraySize2),stat = alloc_stat)
                 if(alloc_stat /=0) then
                    write(err_msg,"(A,I4,I6,A)") & 
                    "Espace memoire disponible insuffisant Coeff_tab2 reception locale Coeff ArraySizes : ( ",&
                     arraySize1,arraySize2," )(read_mat_composite)"
                    call amitex_abort(err_msg ,2)
                 end if 
                 Coeff_tab2 = 0._mytype
                 
                 ! affectation des coefficients
                 if (trim(MattotP0(imR0)%lawname) == "laminate") then
                    if (Ncoeff_comp_list(MattotP0(imR0)%numM_composite(iphi)) /= 2 )then
                       write(err_msg,fmt="(A,I0,A)") &
                           "La loi 'laminate' requiert deux coefficients composites par matériau homogène, ",&
                           Ncoeff_comp_list(MattotP0(imR0)%numM_composite(iphi))," reçus (read_mat_composite)"
                       call amitex_abort(err_msg,1)
                    end if
                    ! affectation des coefficients de la phase propres au modèle laminate (lambda,mu equivalents initiaux)
                    Coeff_tab2(1:Ncoeff_comp_list(MattotP0(imR0)%numM_composite(iphi)),1) = &
                              MattotP0(imS0)%Coeff_comp(:,izS0) 
                                                               
                    ! affectation des coefficients homogènes propres au comportement de la phase
                    Coeff_tab2(1+Ncoeff_comp_list(MattotP0(imR0)%numM_composite(iphi)):,1) = &
                              MattotP0(imS0)%Coeff(:,izS0)     
                 else if (trim(MattotP0(imR0)%lawname) == "reuss") then
                    if (Ncoeff_comp_list(MattotP0(imR0)%numM_composite(iphi)) /= 2 )then
                       write(err_msg,fmt="(A,I0,A)") &
                           "La loi 'reuss' requiert deux coefficients composites par matériau homogène, ",&
                           Ncoeff_comp_list(MattotP0(imR0)%numM_composite(iphi))," reçus (read_mat_composite)"
                       call amitex_abort(err_msg,1)
                    end if
                    ! affectation des coefficients de la phase propres au modèle Reuss (lambda,mu equivalents initiaux)
                    Coeff_tab2(1:Ncoeff_comp_list(MattotP0(imR0)%numM_composite(iphi)),1) = &
                              MattotP0(imS0)%Coeff_comp(:,izS0) 
                                                               
                    ! affectation des coefficients homogènes propres au comportement de la phase
                    Coeff_tab2(1+Ncoeff_comp_list(MattotP0(imR0)%numM_composite(iphi)):,1) = &
                              MattotP0(imS0)%Coeff(:,izS0)     
                 else if (trim(MattotP0(imR0)%lawname) == "voigt") then
                    ! affectation des coefficients homogènes propres au comportement de la phase
                    Coeff_tab2(:,1) = MattotP0(imS0)%Coeff(:,izS0) 
                 end if

                 !Coeff_tab2 = spread(Coeff_tab2(:,1),2,arraySize2) ! bug intel-oneapi 2021.1 (arraysize2 trop grand?)
                 call spread_perso(Coeff_tab2,Coeff_tab2(:,1))  

                 Ind_ini = nCoeff_lawC + sum( MattotP0(imR0)%NCoeff_comp(1:iphi-1)) +1
                 if  (trim(MattotP0(imR0)%lawname) == "laminate") Ind_ini = Ind_ini + 2*(iphi -1)
                 if  (trim(MattotP0(imR0)%lawname) == "reuss")    Ind_ini = Ind_ini + 2*(iphi -1)
                 Ind_end = Ind_ini + arraySize1 - 1

                 MattotP0(imR0)%Coeff(Ind_ini:Ind_end,izR0_list) = Coeff_tab2 

                 if (allocated(Coeff_tab2)) deallocate(Coeff_tab2)


                 !------------------------------------ Variables internes -----------------------------------
                 ! allocation du tableau temporaire Coeff_tab2 qui contiendra l'information à affecter pour les variables internes
                 arraySize1 = MattotP0(imR0)%NVarInt_comp(iphi)
                 if (arraySize1 > 0) then

                    arraySize2 = nZoneR
                    allocate(Coeff_tab2(arraySize1,arraySize2),stat = alloc_stat)
                    if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant &
                                                & Coeff_tab2 reception locale VarInt (read_mat_composite)",2)
                    coeff_tab2 = 0._mytype

                    ! affectation des variables internes - ATTENTION : cela limite l'utilisation des materiaux composites
                    !                                                  a des cas ou les variables internes sont initialisees 
                    !                                                  "Constant" ou "Constant_Zone" MAIS PAS "Variable"!!!
                    Coeff_tab2(:,1) = MattotP0(imS0)%VarInt0(:,MattotP0(imS0)%zone(izS0,1))  

                    ! affectation des coefficients homogènes propres au comportement de la phase
                    !Coeff_tab2(:,:) = spread(Coeff_tab2(:,1),2,arraySize2)  ! bug intel-oneapi 2021.1 (arraysize2 trop grand?)
                    call spread_perso(Coeff_tab2,Coeff_tab2(:,1))  


                    if (trim(MattotP0(imR0)%lawname) == "laminate") then
                       Ind_ini    = MattotP0(imR0)%Nphase + 9*(iphi) + &
                                    sum( MattotP0(imR0)%NVarint_comp(1:iphi-1)) + 1 
                       Ind_end    = Ind_ini + arraySize1 - 1
                    else if (trim(MattotP0(imR0)%lawname) == "reuss") then
                       Ind_ini    = MattotP0(imR0)%Nphase + 12*(iphi) + &
                                    sum( MattotP0(imR0)%NVarint_comp(1:iphi-1)) + 1 
                       Ind_end    = Ind_ini + arraySize1 - 1
                    else if (trim(MattotP0(imR0)%lawname) == "voigt") then
                       Ind_ini    = MattotP0(imR0)%Nphase + 6*(iphi) + &
                                    sum( MattotP0(imR0)%NVarint_comp(1:iphi-1)) + 1 
                       Ind_end    = Ind_ini + arraySize1 - 1
                    end if

                    MattotP0(imR0)%VarInt(Ind_ini:Ind_end,izR0_list) = Coeff_tab2 
                    MattotP0(imR0)%VarInt0(Ind_ini:Ind_end,izR0_list) = &
                                    MattotP0(imR0)%VarInt(Ind_ini:Ind_end,izR0_list)

                    if (allocated(Coeff_tab2)) deallocate(Coeff_tab2)
                 end if

              else 
                 ! matériau à initialiser présent dans le pinceau mais matériau à affecter absent : communication nécessaire
                 if (allocated(Coeff_tab2)) deallocate(Coeff_tab2)

                 !------------------------------------ Coefficients -----------------------------------
                 ! allocation du tableau temporaire Coeff_tab2 qui contiendra l'information à affecter
                 arraySize1 = MattotP0(imR0)%NCoeff_comp(iphi)
                 if  (trim(MattotP0(imR0)%lawname) == "laminate") arraySize1 = arraySize1 + 2
                 if  (trim(MattotP0(imR0)%lawname) == "reuss") arraySize1 = arraySize1 + 2
                 arraySize2 = nZoneR
                 allocate(Coeff_tab2(arraySize1,arraySize2),stat = alloc_stat)
                 if(alloc_stat /=0) then
                    write(err_msg,"(A,I4,I6,A)") "Espace memoire disponible insuffisant Coeff_tab2 pour&
                                & reception depuis un pinceau externe des Coeff , ArraySizes : ( ",&
                                  arraySize1,arraySize2," )(read_mat_composite)"
                    call amitex_abort(err_msg ,2)
                 end if
                 Coeff_tab2 = 0._mytype
                 
                 ! affectation des coefficients
                 !------------------------------
                 if (trim(MattotP0(imR0)%lawname) == "laminate" .or. trim(MattotP0(imR0)%lawname) == "reuss") then
                    
                    if (Ncoeff_comp_list(MattotP0(imR0)%numM_composite(iphi)) /= 2 )then
                       write(err_msg,fmt="(A)") &
                           "La loi laminate requiert deux coefficients composites par matériau homogène, ",&
                           Ncoeff_comp_list(MattotP0(imR0)%numM_composite(iphi))," reçus (read_mat_composite)"
                       call amitex_abort(err_msg,2)
                    end if
                    ! Reception coefficients composites
                    etiquette = 100 ! 100 correspond à l'envoi de MattotP0(imS0)%Coeff_composites(:,izS0) depuis le pinceau rksend1
                    call MPI_RECV( Coeff_tab2(1:Ncoeff_comp_list(MattotP0(imR0)%numM_composite(iphi)),1),&
                          Ncoeff_comp_list(MattotP0(imR0)%numM_composite(iphi)),real_type,rksend1,etiquette,&
                          amitex_mpi_comm,statut,ierror)

                    ! Reception coefficients homogène de la combinaison Materiau/Zone à affecter pour les coefficients
                    etiquette = 101 ! 101 correspond à l'envoi de MattotP0(imS0)%Coeff(:,izS0) depuis le pinceau rksend1
                    call MPI_RECV( Coeff_tab2(1+Ncoeff_comp_list(MattotP0(imR0)%numM_composite(iphi)):size(Coeff_tab2,1),1),&
                             Ncoeff_list(MattotP0(imR0)%numM_composite(iphi)),real_type,rksend1,etiquette,&
                             amitex_mpi_comm,statut,ierror) 

                 else if (trim(MattotP0(imR0)%lawname) == "voigt") then

                    ! Reception coefficients homogène de la combinaison Materiau/Zone à affecter pour les coefficients
                    etiquette = 101 ! 101 correspond à l'envoi de MattotP0(imS0)%Coeff(:,izS0) depuis le pinceau rksend1
                    call MPI_RECV( Coeff_tab2(:,1),Ncoeff_list(MattotP0(imR0)%numM_composite(iphi)),&
                                   real_type,rksend1,etiquette,amitex_mpi_comm,statut,ierror) 

                 end if

                 !Coeff_tab2 = spread(Coeff_tab2(:,1),2,arraySize2) ! bug intel-oneapi 2021.1 (arraysize2 trop grand?)
                 call spread_perso(Coeff_tab2,Coeff_tab2(:,1))  

                 Ind_ini = nCoeff_lawC + sum( MattotP0(imR0)%NCoeff_comp(1:iphi-1)) +1
                 if  (trim(MattotP0(imR0)%lawname) == "laminate") Ind_ini = Ind_ini + 2*(iphi -1)
                 if  (trim(MattotP0(imR0)%lawname) == "reuss")    Ind_ini = Ind_ini + 2*(iphi -1)
                 Ind_end = Ind_ini + arraySize1 - 1
                 MattotP0(imR0)%Coeff(Ind_ini:Ind_end,izR0_list) = Coeff_tab2  !! ATTENTION INDICES A GAUCHE
 

                 if (allocated(Coeff_tab2)) deallocate(Coeff_tab2)

                 !------------------------------------ Variables internes -----------------------------------
                 ! allocation du tableau temporaire Coeff_tab2 qui contiendra l'information à affecter pour les variables internes
                 arraySize1 = MattotP0(imR0)%NVarInt_comp(iphi)
                 if (arraySize1 > 0) then
                    arraySize2 = nZoneR
                    allocate(Coeff_tab2(arraySize1,arraySize2),stat = alloc_stat)
                    if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant Coeff_tab2& 
                             & pour reception depuis un pinceau externe des VarInt (read_mat_composite)",2)

                    ! Reception des variables internes homogène de la combinaison Materiau/Zone à affecter pour les coefficients
                    etiquette = 200 ! 200 correspond à l'envoi de MattotP0(imS0)%VarInt(:,izS0) depuis le pinceau rksend1
                    call MPI_RECV( Coeff_tab2(:,1),arraySize1,real_type,rksend1,&
                                   etiquette,amitex_mpi_comm,statut,ierror) 

                    !Coeff_tab2 = spread(Coeff_tab2(:,1),2,arraySize2) ! bug intel-oneapi 2021.1 (arraysize2 trop grand?)
                    call spread_perso(Coeff_tab2,Coeff_tab2(:,1))  

                    if (trim(MattotP0(imR0)%lawname) == "laminate") then
                       Ind_ini    = MattotP0(imR0)%Nphase + 9*(iphi) + &
                                    sum( MattotP0(imR0)%NVarint_comp(1:iphi-1)) + 1 
                       Ind_end    = Ind_ini + arraySize1 - 1
                    else if (trim(MattotP0(imR0)%lawname) == "reuss") then
                       Ind_ini    = MattotP0(imR0)%Nphase + 12*(iphi) + &
                                    sum( MattotP0(imR0)%NVarint_comp(1:iphi-1)) + 1 
                       Ind_end    = Ind_ini + arraySize1 - 1
                    else if (trim(MattotP0(imR0)%lawname) == "voigt") then
                       Ind_ini    = MattotP0(imR0)%Nphase + 6*(iphi) + &
                                    sum( MattotP0(imR0)%NVarint_comp(1:iphi-1)) + 1 
                       Ind_end    = Ind_ini + arraySize1 - 1
                    end if

                    MattotP0(imR0)%VarInt(Ind_ini:Ind_end,izR0_list) = Coeff_tab2 
                    MattotP0(imR0)%VarInt0(Ind_ini:Ind_end,izR0_list) = &
                                    MattotP0(imR0)%VarInt(Ind_ini:Ind_end,izR0_list)

                 end if

              end if
              ! le matériau et des voxels composites à initialiser sont présents sur le pinceau Fin récupération des données
           end if

           ! Envoi de données
           if (nrank == rksend1) then
              ! Envoi des coefficiens Composites

              do i =1,size(rkrecv1) ! Boucle sur les pinceau : envoi uniquement à ceux ayant besoin d'informations indisponibles localement
                 if (rkrecv1(i) /= -1) then

                    if (trim(MatComposite(icomp)%Lawname) == "laminate" .or. trim(MatComposite(icomp)%Lawname) == "reuss") then
                       ! Envoi coefficients composites
                       etiquette = 100 ! 100 correspond à l'envoi de MattotP0(imS0)%Coeff_composites(:,izS0) depuis le pinceau rksend1
                       call MPI_SEND( MattotP0(imS0)%Coeff_comp(:,izS0),&
                               size(MattotP0(imS0)%Coeff_comp(:,izS0)),real_type,rkrecv1(i),etiquette,&
                               amitex_mpi_comm,ierror)
                    end if

                    ! Envoi coefficients homogènes
                    etiquette = 101 ! 101 correspond à l'envoi de MattotP0(imS0)%Coeff(:,izS0) depuis le pinceau rksend1
                    call MPI_SEND(MattotP0(imS0)%Coeff(:,izS0) ,&
                             size(MattotP0(imS0)%Coeff(:,izS0)),real_type,rkrecv1(i),etiquette,&
                             amitex_mpi_comm,ierror) 
                    
                    ! Envoi variables internes homogènes
                    etiquette = 200 ! 200 correspond à l'envoi de MattotP0(imS0)%VarInt(:,izS0) depuis le pinceau rksend1
                    call MPI_SEND( MattotP0(imS0)%VarInt0(:,MattotP0(imS0)%Zone(izS0,1)) ,&
                             size(MattotP0(imS0)%VarInt0(:,MattotP0(imS0)%Zone(izS0,1))),real_type,rkrecv1(i),etiquette,&
                             amitex_mpi_comm,ierror) 

                 end if
              end do

           end if
           
           if (allocated(izR0_list)) then
              deallocate(izR0_list)
           end if

           !--------------------------------------------------
           ! FIN AFFECTATIONS LOCALES
           !--------------------------------------------------

           ! fin boucle sur les zones du matériau composite icomp (izone)
        end do

        if (allocated(ZoneInit)) then
           deallocate(ZoneInit)
        end if
        
        deallocate(array_int64_1d)
        deallocate(array_real_1d)
        ! fin boucle sur les phases du matériau composite icomp (iphi)
     end do
     
     ! Synchronisation des processus après chaque initialisation de matériau
     call MPI_BARRIER(amitex_mpi_comm,ierror)
     if (nrank == 0) write(OUTPUT_UNIT,'(A,I4,A)') &
                    "Initialisation matériau composite : ",icomp," terminée."

     ! fin boucle sur les matériaux composites de la cellule (icomp)
  end do

  call initParam_umat()

end subroutine read_mat_composite


!===================================================================
!
!               SUBROUTINE ASSEMBLE_FIELD_FROM_VARINT
!               SUBROUTINE ASSEMBLE_FIELD_FROM_VARINT0
!
!> Assemble les valeurs de variables internes contenues dans la structure
!! MattotP pour former un champ. 
!! Uniquement pour les champs de reels
!!
!!   \param[in] : list_mat   [vecteur d'entier]  
!!                liste des materiaux possedant les variables internes
!!                constituant le champ
!!   \param[in] : list_varInt
!!                liste des indices des variables internes dont sont 
!!                les differentes composantes du champ a construire
!!                pour chaque materiau precise dans liste_mat
!!                Chaque materiau peut avoir des indices de variables
!!                internes differents, neanmoins le nombre d'indices
!!                reste le meme == nbre de composantes du champ a construire
!!                [tableaux d'entiers kind=8 nmat x nindices]
!!
!!   \param[in] : ncomp
!!                nombre de variables internes a recuperer = nbre de composantes
!!                du champ
!!
!!   \param[out ] : Field
!!                  Champ distribue sur la decomposition 2decomp
!!                  a construire a partir des valeurs de variables
!!                  internes precisees par liste_VarInt et liste_mat        
!!
!!    ATTENTION : NON ADAPTE AU CAS DE MATERIAUX COMPOSITES
!!                ADAPTATION A FAIRE EN UTILISANT LE TRAVAIL DE 
!!                J. DUVERGE :
!!                structure qui recense quel materiau composite 
!!                contient quelle phase/zone etc.. pour faciliter
!!                traitement des sorties --> adapte a la construction
!!                de champ egalement
!!                + NECESSITE UN CHOIX : comment definir le champ sur 
!!                les voxels composites ???        
!!
!===================================================================
subroutine Assemble_field_from_VarInt(list_mat,list_varInt,ncomp,Field)

  implicit none
  ! Entrees :
  integer,dimension(:),intent(in)                                       :: list_mat
  integer,dimension(:,:),intent(in)                                     :: list_varInt
  integer, intent(in)                                                   :: ncomp
  ! Sorties :   
  real(mytype),dimension(amx_size(1)*amx_size(2)*amx_size(3),ncomp),intent(out)  :: Field
  ! Autres variables 
  integer                                                               :: i,j,k
  !integer,dimension(2)                                                  :: order2 = (/ 2, 1 /)
  
  !! Initialisation du champ a 0.
  !Field(:,:) = 0.

  !! Remplissage du champ
  do i=1,size(list_mat)
     !! Materiau par materiau
 
     ! On parcours les materiaux sur le pinceau
     do j=1,size(mattotP)
        ! Champ rempli uniquement si le materiau precise se trouve sur le pinceau
        if (mattotP(j)%numM == list_mat(i)) then

            do k=1,ncomp
                Field(mattotP(j)%pos,k) = mattotP(j)%VarInt(list_varInt(i,k),:)
            end do
        
          ! ANCIENNE VERSION : 
          ! probleme en compilation intel avec reshape pour des gros champs --> segmentation fault
          ! mal compris.....
          ! Les tableaux ont des formes inverses (NposxNvarInt) vs (NvarIntxNpos)
          ! l'option order de reshape permet d'inverser l'ordre d'ordonancement des valeurs
          ! du second tableau pour faire coincider les formes
          !     Field(mattotP(j)%pos,:) = reshape(mattotP(j)%VarInt(list_varInt(i,:),:),&
          !                             (/size(mattotP(j)%pos),size(list_varInt(i,:)) /),order=order2) 
                                    
        end if     
     end do
  end do

end subroutine Assemble_field_from_VarInt
!-----------------------------------------------------------------------
! copy-paste of the previous function - replacing Varint by Varint0 
!-----------------------------------------------------------------------
subroutine Assemble_field_from_VarInt0(list_mat,list_varInt,ncomp,Field)

  implicit none
  ! Entrees :
  integer,dimension(:),intent(in)                                       :: list_mat
  integer,dimension(:,:),intent(in)                                     :: list_varInt
  integer, intent(in)                                                   :: ncomp
  ! Sorties :   
  real(mytype),dimension(amx_size(1)*amx_size(2)*amx_size(3),ncomp),intent(out)  :: Field
  ! Autres variables 
  integer                                                               :: i,j,k
  
  !! Initialisation du champ a 0.
  !Field(:,:) = 0.

  !! Remplissage du champ
  do i=1,size(list_mat)
     !! Materiau par materiau
     ! On parcours les materiaux sur le pinceau
     do j=1,size(mattotP)
        ! Champ rempli uniquement si le materiau precise se trouve sur le pinceau
        if (mattotP(j)%numM == list_mat(i)) then
            do k=1,ncomp
                Field(mattotP(j)%pos,k) = mattotP(j)%VarInt0(list_varInt(i,k),:)
            end do
        end if     
     end do
  end do

end subroutine Assemble_field_from_VarInt0

!===================================================================
!
!               SUBROUTINE ASSIGN_VARINT_FROM_FIELD
!               SUBROUTINE ASSIGN_VARINT0_FROM_FIELD
!
!> Met a jour les valeurs de variables internes contenues dans la structure
!! MattotP a parti d'un champ  
!! Uniquement pour les champs de reels
!!
!!   \param[in] : list_mat   [vecteur d'entier]  
!!                liste des materiaux possedant les variables internes
!!                constituant le champ
!!   \param[in] : list_varInt
!!                liste des indices des variables internes dont sont 
!!                les differentes composantes du champ a construire
!!                pour chaque materiau precise dans liste_mat
!!                Chaque materiau peut avoir des indices de variables
!!                internes differents, neanmoins le nombre d'indices
!!                reste le meme == nbre de composantes du champ a construire
!!                [tableaux d'entiers kind=8 nmat x nindices]
!!
!!   \param[in] : ncomp
!!                nombre de variables internes a recuperer = nbre de composantes
!!                du champ
!!
!!   \param[in] : Field
!!                  Champ distribue sur la decomposition 2decomp
!!                  a construire a partir des valeurs de variables
!!                  internes precisees par liste_VarInt et liste_mat        
!!
!!    ATTENTION : NON ADAPTE AU CAS DE MATERIAUX COMPOSITES
!!                ADAPTATION A FAIRE EN UTILISANT LE TRAVAIL DE 
!!                J. DUVERGE :
!!                structure qui recense quel materiau composite 
!!                contient quelle phase/zone etc.. pour faciliter
!!                traitement des sorties --> adapte a la construction
!!                de champ egalement
!!                + NECESSITE UN CHOIX : comment definir le champ sur 
!!                les voxels composites ???        
!!
!===================================================================
subroutine Assign_VarInt_from_Field(list_mat,list_VarInt,ncomp,Field)

  implicit none
  ! Entrees :
  integer,dimension(:),intent(in)                                        :: list_mat
  integer,dimension(:,:),intent(in)                                      :: list_varInt
  integer,intent(in)                                                     :: ncomp
  real(mytype),dimension(amx_size(1)*amx_size(2)*amx_size(3),ncomp),intent(in)    :: Field
  ! Autres variables 
  integer                                                                :: i,j,k
  !integer,dimension(2)                                                   :: order2 = (/ 2, 1 /)

  do i=1,size(list_mat)
     !! Materiau par materiau
 
     ! On parcours les materiaux sur le pinceau
     do j=1,size(mattotP)
        ! Champ remplit uniquement si le materiau precise se trouve sur le pinceau
        if (mattotP(j)%numM == list_mat(i)) then
            do k=1,ncomp
                mattotP(j)%VarInt(list_varInt(i,k),:) = Field(mattotP(j)%pos,k)
            end do 
            
          ! ANCIENNE VERSION : 
          ! probleme en compilation intel avec reshape pour des gros champs --> segmentation fault
          ! mal compris.....
          ! Les tableaux ont des formes inverses (NposxNvarInt) vs (NvarIntxNpos)
          ! l'option order de reshape permet d'inverser l'ordre d'ordonancement des valeurs
          ! du second tableau pour faire coincider les formes
          !   mattotP(j)%VarInt(list_varInt(i,:),:) = reshape(Field(mattotP(j)%pos,:),&
          !                            (/size(list_varInt(i,:)),size(mattotP(j)%pos) /),order=order2)

        end if     
     end do
  end do

end subroutine Assign_VarInt_from_Field
!-----------------------------------------------------------------------
! copy-paste of the previous function - replacing Varint by Varint0 
!-----------------------------------------------------------------------
subroutine Assign_VarInt0_from_Field(list_mat,list_VarInt,ncomp,Field)

  implicit none
  ! Entrees :
  integer,dimension(:),intent(in)                                        :: list_mat
  integer,dimension(:,:),intent(in)                                      :: list_varInt
  integer,intent(in)                                                     :: ncomp
  real(mytype),dimension(amx_size(1)*amx_size(2)*amx_size(3),ncomp),intent(in)    :: Field
  ! Autres variables 
  integer                                                                :: i,j,k

  do i=1,size(list_mat)
     !! Materiau par materiau
     ! On parcours les materiaux sur le pinceau
     do j=1,size(mattotP)
        ! Champ remplit uniquement si le materiau precise se trouve sur le pinceau
        if (mattotP(j)%numM == list_mat(i)) then
            do k=1,ncomp
                mattotP(j)%VarInt0(list_varInt(i,k),:) = Field(mattotP(j)%pos,k)
            end do             
        end if     
     end do
  end do

end subroutine Assign_VarInt0_from_Field

!==============================================================================
!
!>   SUBROUTINE get_homogeneous_coeff
!!
!!   Get and assign a homogeneous coefficient from a non-local model
!!
!!   1 - Check the coefficient (given through in MattotP%Coeff) is homogeneous
!!       ERROR if not homogeneous
!!   2 - Assign the coefficient to the whole cell (even if the non-local model
!!       is not present on a pencil) 
!!
!! param[out]    coeff,     output value   
!! param[in]     modelInd,  index of the non-local model
!! param[in]     ind_coeff, index of the coefficient in the non-local model
!! param[in]     list_mat(:),  list of materials numM in the non-local model 
!!
!==============================================================================
subroutine get_homogeneous_coeff(coeff, modelInd, ind_coeff, list_mat)

 implicit none
 real(mytype), intent(out)       :: coeff      ! output value 
 integer,intent(in)              :: modelInd   ! index of the non-local model
 integer,intent(in)              :: ind_coeff  ! index of the coefficient in the non-local model
 integer,dimension(:),intent(in) :: list_mat   ! list of materials numM in the non-local model 

 integer                         :: i,j, ierr
 integer                         :: icoeff ! index of the coeff in the mattotP%Coeff(icoeff,izone)
 real(mytype)                    :: maxv,minv, coeff_tmp1, coeff_tmp2
 logical                         :: testP  ! test presence du modele sur le pinceau
 character(len=200)              :: err_msg

 maxv = 0
 minv = 0
 icoeff = 0
 coeff = 0
 coeff_tmp1 = huge(maxv)
 coeff_tmp2 = 0.
 testP = .false. 

 ! error message
 write(err_msg,fmt="(A,I0,A,I0,A)") "non-local coeff ",ind_coeff," for model ",&
                                    modelInd," is not homogeneous (get_homogeneous_coeff)"

 
 !--- Homogeneite par pinceau (si ok -> coeff_tmp1 = coeff_tmp2)
 do i=1,size(list_mat)
    do j=1,size(mattotP)
    if (mattotP(j)%numM == list_mat(i)) then
       testP = .true.
       !------ tests homogeneite entre zones d'un materiau (sur un pinceau)
       icoeff  = Nloc_models(modelInd)%Ind_CoeffNloc(i,ind_coeff) 
       maxv   = maxval(MattotP(j)%Coeff(icoeff,:))      
       minv   = minval(MattotP(j)%Coeff(icoeff,:))      
       if (maxv /= minv) then 
           ! ERROR
           call amitex_abort(trim(err_msg),2)
       else
           coeff_tmp2 = maxv
       end if
       !------ tests homogeneite entre materiaux (sur un pinceau)
       if (coeff_tmp1 /=huge(maxv) .and. coeff_tmp2 /= coeff_tmp1) then 
           ! ERROR
           call amitex_abort(trim(err_msg),2)
       end if
       coeff_tmp1 = coeff_tmp2
    end if
    end do
 end do

 !--- Homogeneite sur les pinceaux et affectation si le model n'est pas present sur le pinceau
 if (.not. testP) coeff_tmp1 = -huge(coeff_tmp1)
 call MPI_allreduce(coeff_tmp1,maxv,1,real_type,MPI_MAX,amitex_mpi_comm,ierr)
 
 if (.not. testP) coeff_tmp1 = huge(coeff_tmp1)  
 call MPI_allreduce(coeff_tmp1,minv,1,real_type,MPI_MIN,amitex_mpi_comm,ierr)
 
 if (maxv==minv) then
    coeff = maxv
 else
    !ERROR
    call amitex_abort(trim(err_msg),2)
 end if

end subroutine get_homogeneous_coeff

!==============================================================================
!
!>   SUBROUTINE get_minmax_coeff
!!
!!   Get and assign the min and max values of a coefficient from a non-local model
!!
!! param[out]    coeffmin,     output value   
!! param[out]    coeffmax,     output value   
!! param[in]     modelInd,  index of the non-local model
!! param[in]     ind_coeff, index of the coefficient in the non-local model
!! param[in]     list_mat(:),  list of materials numM in the non-local model 
!!
!==============================================================================
subroutine get_minmax_coeff(coeffmin,coeffmax, modelInd, ind_coeff, list_mat)

 implicit none
 real(mytype), intent(out)       :: coeffmin, coeffmax  ! output value 
 integer,intent(in)              :: modelInd            ! index of the non-local model
 integer,intent(in)              :: ind_coeff           ! index of the coefficient in the non-local model
 integer,dimension(:),intent(in) :: list_mat            ! list of materials numM in the non-local model 

 real(mytype)                    :: tmp_min, tmp_max
 integer                         :: i,j, icoeff
 integer                         :: ierror
 
 !--- min et max sur les pinceaux + initilisation (utile si coeff non present sur le pinceaux)
 coeffmin = huge(coeffmin)
 coeffmax = -huge(coeffmax)
 do i=1,size(list_mat)
    do j=1,size(mattotP)
    if (mattotP(j)%numM == list_mat(i)) then
       icoeff  = Nloc_models(modelInd)%Ind_CoeffNloc(i,ind_coeff) 
       tmp_max   = maxval(MattotP(j)%Coeff(icoeff,:))      
       tmp_min   = minval(MattotP(j)%Coeff(icoeff,:)) 
       coeffmin = min(coeffmin, tmp_min)
       coeffmax = max(coeffmax, tmp_max)     
    end if
    end do
 end do
 
 call MPI_ALLREDUCE(MPI_IN_PLACE, coeffmin, 1, real_type, MPI_MIN, amitex_mpi_comm, ierror)
 call MPI_ALLREDUCE(MPI_IN_PLACE, coeffmax, 1, real_type, MPI_MAX, amitex_mpi_comm, ierror)

end subroutine get_minmax_coeff

!==============================================================================
!
!>   SUBROUTINE multiply_field_coeff
!!
!!   Multiply a field by a coefficient from a non-local model
!!
!! param[inout]  AA field input / output   
!! param[in]     modelInd,  index of the non-local model
!! param[in]     ind_coeff, index of the coefficient in the non-local model
!! param[in]     list_mat(:),  list of materials numM in the non-local model 
!! param[in]     ncomp,  number of components of the field
!! param[in]     nvar,  number of variable in the field 
!!               
!! Note : AA is reshaped as (amx_size(1)*amx_size(2)*amx_size(3),ncomp,nvar)
!!
!==============================================================================
subroutine multiply_field_coeff(AA, modelInd, ind_coeff,list_mat, ncomp, nvar)

 implicit none
 integer,intent(in)              :: modelInd            ! index of the non-local model
 integer,intent(in)              :: ind_coeff           ! index of the coefficient in the non-local model
 integer,dimension(:),intent(in) :: list_mat            ! list of materials numM in the non-local model 
 integer                         :: ncomp,nvar
 real(mytype), dimension(amx_size(1)*amx_size(2)*amx_size(3),ncomp,nvar),intent(inout)  :: AA

 real(mytype)                    :: coeff
 integer                         :: i,j,k,icoeff
 integer(INT64)                  :: lmin, lmax
 
 do i=1,size(list_mat)
    do j=1,size(mattotP)
    if (mattotP(j)%numM == list_mat(i)) then
       icoeff  = Nloc_models(modelInd)%Ind_CoeffNloc(i,ind_coeff) 
       lmin=1

       do k=1, size(mattotP(j)%zone,1)
         coeff = MattotP(j)%Coeff(icoeff,k)
         
         if (k>1) lmin = MattotP(j)%zone(k-1,1) + 1
         lmax = MattotP(j)%zone(k,1)
         
         AA(MattotP(j)%pos(lmin:lmax),:,:) = AA(MattotP(j)%pos(lmin:lmax),:,:) * coeff
         
       end do
    end if
    end do
 end do
 

end subroutine multiply_field_coeff

!==============================================================================
!
!>   SUBROUTINE get_minmax_coeff_function
!!
!!   Get and assign the min and max values of a function of coefficients
!!             from a non-local model
!!
!!   The function is an external function defined by the user 
!!               (-> see abstract interface in the module)
!!
!!
!!   Example : - if A is the first coeff of the non-local modem, and B the 3rd
!!             - the user wants the min/max of B/A
!!    Then  :  - use (/1,3/) for ind_coeff
!!             - define an external function fun 
!!                             real(mytype),dimension(size(coeff,2)) fun(coeff)
!!                             fun(:) = coeff(2,:) / coeff(1,:)  !  2 corresponds to index 3 (B)
!!                                                               !  1 corresponds to index 1 (A) 
!!
!!
!!
!! param[out]    valmin,     output value   
!! param[out]    valmax,     output value   
!! param[in]     modelInd,  index of the non-local model
!! param[in]     function of coefficient (input coeff(1:n_ind_coeff,Nzones), output : values(Nzones)
!! param[in]     list_mat(:),  list of materials numM in the non-local model 
!!
!==============================================================================
subroutine get_minmax_coeff_function(valmin,valmax, modelInd, fun0, list_mat)
 implicit none
 real(mytype), intent(out)       :: valmin, valmax  ! output value 
 integer,intent(in)              :: modelInd            ! index of the non-local model
 procedure(func0)                :: fun0
 integer,dimension(:),intent(in) :: list_mat            ! list of materials numM in the non-local model 

 

 real(mytype)                            :: tmp_min, tmp_max
 integer                                 :: i,j,k
 integer                                 :: icoeff, ncoeff, nzone
 real(mytype),allocatable,dimension(:,:) :: coeff
 real(mytype),allocatable,dimension(:)   :: values
 integer                                 :: ierror
 
 
 ncoeff = size(Nloc_models(modelInd)%Ind_CoeffNloc,2)
 
 !--- min et max sur les pinceaux + initilisation (utile si coeff non present sur le pinceaux)
 valmin = huge(valmin)
 valmax = -huge(valmax)
 do i=1,size(list_mat)
    do j=1,size(mattotP)
    if (mattotP(j)%numM == list_mat(i)) then
       nzone = size(MattotP(list_mat(i))%Coeff,2)
       allocate(coeff(ncoeff,nzone))
       allocate(values(nzone))
       do k = 1,ncoeff
          icoeff  = Nloc_models(modelInd)%Ind_CoeffNloc(i,k)
          coeff(k,:) = MattotP(list_mat(i))%Coeff(icoeff,:)
       end do
       values = fun0(coeff)   ! user definition
       tmp_max = maxval(values)
       tmp_min = minval(values)
       valmax = max(valmax, tmp_max)
       valmin = min(valmin, tmp_min)
       deallocate(coeff)
       deallocate(values)
    end if  
    end do
 end do
 
 call MPI_ALLREDUCE(MPI_IN_PLACE, valmin, 1, real_type, MPI_MIN, amitex_mpi_comm, ierror)
 call MPI_ALLREDUCE(MPI_IN_PLACE, valmax, 1, real_type, MPI_MAX, amitex_mpi_comm, ierror)


end subroutine get_minmax_coeff_function

!==============================================================================
!
!>   FUNCTION apply_field_coeff_function
!!
!!   Apply a function on a field  with coefficient from a non-local model :
!!
!!   symbolically:   Fout = func1(Fin, coeff)
!!
!! param[out]    Fout 
!! param[in]     Fin  
!! param[in]     modelInd,  index of the non-local model
!! param[in]     list_mat(:),  list of materials numM in the non-local model 
!! param[in]     ncomp,  number of components of the field
!! param[in]     nvar,  number of variable in the field 
!!               
!! Note : Fin is reshaped as (amx_size(1)*amx_size(2)*amx_size(3),ncomp,nvar)
!!
!==============================================================================
subroutine apply_field_coeff_function(Fin, Fout, modelInd,list_mat, ncomp, nvar,fun1)

 implicit none
 integer,intent(in)              :: modelInd            ! index of the non-local model
 integer,dimension(:),intent(in) :: list_mat            ! list of materials numM in the non-local model 
 integer                         :: ncomp,nvar
 real(mytype), dimension(amx_size(1)*amx_size(2)*amx_size(3),ncomp,nvar),intent(in)   :: Fin
 real(mytype), dimension(amx_size(1)*amx_size(2)*amx_size(3),ncomp,nvar)              :: Fout
 procedure(func1)                 :: fun1

 real(mytype),allocatable,dimension(:,:)  :: coeff
 integer                                  :: i,j,k,icoeff, ncoeff
 integer(INT64)                           :: lmin, lmax, nzone
 
 ncoeff = size(Nloc_models(modelInd)%Ind_CoeffNloc,2)

  do i=1,size(list_mat)
    do j=1,size(mattotP)
    if (mattotP(j)%numM == list_mat(i)) then
       nzone = size(MattotP(list_mat(i))%Coeff,2)
       allocate(coeff(ncoeff,nzone))
       
       !--- get table coeff(size(ind_coeff),nzone)
       do k = 1,ncoeff
          icoeff  = Nloc_models(modelInd)%Ind_CoeffNloc(i,k)
          coeff(k,:) = MattotP(list_mat(i))%Coeff(icoeff,:)
       end do
       
       !--- apply function on zones 
       lmin=1
       do k=1, size(mattotP(j)%zone,1) !loop on zones
         
         if (k>1) lmin = MattotP(j)%zone(k-1,1) + 1
         lmax = MattotP(j)%zone(k,1)
         
         Fout(MattotP(j)%pos(lmin:lmax),:,:) = fun1(Fin(MattotP(j)%pos(lmin:lmax),:,:), coeff(:,k))
         
       end do
       deallocate(coeff)
    end if

    end do
  end do

end subroutine apply_field_coeff_function


!=======================================================================
!>           FONCTION UMATVOIGT
!!-----------------------------------------------------------------------
!! integration d'un modele composite de Voigt à N phases (déformation homogène)
!! formalisme HPP
!!
!!-----------------------------------------------------------------------
!!
!!   \param[in]   umatTab Tableau du type dérivé Umatptr
!!                - contient les pointeurs de fonction vers le comportement de chaque phase
!!
!!    ============================================================================
!!      COEFFICIENTS
!!      -----------------
!!
!!                Coeff(ncoeff) tableau de coefficients, ordre defini ci-dessous
!!                    avec ncoeff = 2N+sum(ncoeffi)
!!
!!      1  - nPhase     nombre de phases dans le matériau composite
!!      2  - ncoeff1    nombre de coeff de chaque loi
!!      3  - ncoeff2
!!           .......
!!      1+N - ncoeffN
!!      
!!      2+N - FV 1      fractions volumiques de chaque phase
!!      3+N - FV 2
!!            ....
!!      2N  - FV N-1
!!
!!      1+2N,2N+ncoeff1
!!                 Coeff1            coefficients pour le materiau 1
!!      1+2N+ncoeff1, 2N+ncoeff1+ncoeff2
!!                 Coeff2            coefficients pour le materiau 2
!!            ....
!!            ....
!!      1+2N+ncoeff1+...+ncoeffn-1, 2N+ncoeff1+......+ncoeffN
!!                  CoeffN           coefficients pour le materiau n
!!
!!    ============================================================================
!!
!!    ============================================================================
!!      VARIABLES INTERNES
!!      -----------------------
!!
!!      Varint(nvarint) tableau de variables internes, ordre defini ci-dessous
!!                      avec nvarint= 7N + sum(nVari)
!!
!!      1  - nvar1            nombre de var. internes de chaque loi
!!      2  - nvar2
!!           .....
!!      N  - nvarN
!!
!!           1+N,6+N
!!               Sig1            Contraintes dans la phase 1 (notation Voigt)
!!           7+N,6+N+nVar1
!!               VarInt1         Variables internes de la phase 1
!!               ........
!!               ........
!!           1+6*(i-1)+N+nVar1+...+nVari-1,6*i+N+nVar1+....+nVari-1
!!               Sigi           Contraintes dans la phase i (notation Voigt)
!!           1+6*i+N+nVar1+...+nVari-1,6*i+N+nvar1+...+nVari
!!               VarintI        Variables internes de la phase i
!!               ........
!!               ........
!!           1+6*(N-1)+N+nVar1+...+nVarn-1,6*N+N+nVar1+.....+nVarn-1
!!               SigN       contraintes dans la phase N
!!
!!           1+6*N+N+nVar1+.....+nVarn-1,6*N+N+sum(nVari)
!!               VarintN    Variables internes dans la phase N
!!
!!    ============================================================================
!!
!!
!! \param[out] Sig, Varint
!!
!!-------------------------------------------------------------------------
!! ATTENTION - ATTENTION - ATTENTION : notation de Voigt, 
!!            pour la deformation : x2 sur les termes de cisaillement 
!!                        soit ici NT0,NT1 et T0T1 
!!            pour la contrainte  : pas de facteur 2
!!-------------------------------------------------------------------------
!!
!! Notation CAST3M (Voigt, ordre 11 22 33 12 13 23)
!=======================================================================


!=======================================================================
!>           FONCTION UMATREUSS
!!-----------------------------------------------------------------------
!! integration d'un modele composite de REUSS a N phases
!!          
!! formalisme HPP
!!
!!-----------------------------------------------------------------------
!!
!! \param[in]   dt	Vecteur des pas de temps nouveau et ancien (dt_new,dt_old)
!! \param[in]   umatTab Tableau du type dérivé Umatptr contenant les pointeurs 
!!                      de fonction vers le comportement de chaque phase
!!
!!    ==========================================================================
!!      COEFFICIENTS
!!      ----------------------
!!
!!	Coeff(ncoeff) tableau de coefficients, ordre defini ci-dessous
!!                    avec ncoeff = 4*N+ncoeff1+.............+ncoeffn
!!
!!
!!	1    - nPhase                           nombre de phases dans le matériau composite
!!	2    - ncoeff1                          nombre de coeff de chaque loi
!!	3    - ncoeff2
!!             .......
!!	1+N  - ncoeffN
!!	
!!      2+N  - FV 1	                 	fractions volumiques de chaque phase (n-1 données la nieme vérifiant FV N = 1 - sum(FV I)
!!      3+N  - FV 2
!!             ....
!!      2N   - FV N-1
!!
!!      2N+1,2N+2
!!             (/ Lambdaeq, Mueq /)             coefficients de Lamé équivalents au comportement élastique initial du matériau 1
!!      2N+3,2N+2+ncoeff1
!!             Coeff1                           coefficients pour le materiau 1
!!      2N+3+ncoeff1,2N+4+ncoeff1
!!             (/ Lambdaeq, Mueq /)             coefficients de Lamé équivalents au comportement élastique initial du matériau 2
!!      2N+5+ncoeff1, 2N+4+ncoeff1+ncoeff2
!!             Coeff2                           coefficients pour le materiau 2
!!            ....
!!            ....
!!      2N+2*(N-1)+ncoeff1+.....+ncoeffn-1+1,4*N+ncoeff1+.............+ncoeffn-1
!!             (/ Lambdaeq, Mueq /)             coefficients de Lamé équivalents au comportement élastique initial du matériau N
!!      4*N+ncoeff1+.............+ncoeffn-1+1,4*N+ncoeff1+.............+ncoeffn
!!             CoeffN                           coefficients pour le materiau n
!!    ===========================================================================
!!
!!
!!    ===========================================================================
!!      VARIABLES INTERNES 
!!      -----------------------
!!
!!      Varint(nvarint) tableau de variables internes, ordre defini ci-dessous
!!                      avec nvarint= 13N + sum(nVari)
!!
!!      1  - nvar1            nombre de var. internes de chaque loi
!!      2  - nvar2
!!           .....
!!      N  - nvarN
!!
!!      1+N,6+N
!!           Def1                        composantes de la defomation, phase 1 
!!
!!      7+N,12+N
!!           Def1_old                    composantes de la deformation pas precedent, phase 1
!!
!!      13+N,12+N+nvar1
!!           VarInt1                     variables internes, phase 1
!!
!!      13+N+nvar1,18+N+nvar1
!!           Def2                        composantes de la defomation, phase 2 
!!
!!      19+N+nvar1,24+N+nvar1
!!           Def2_old                    composantes de la deformation pas precedent, phase 2
!!                                      
!!      25+N+nvar1,24+N+nvar1+nvar2
!!           VarInt2                     variables internes, phase2
!!
!!       .........
!!       .........
!!       .........
!!
!!      N+12*(N-1)+1+nvar1+...+nvarn-1,N+12*(N-1)+6+nvar1+...+nvarn-1
!!           Defn                        composantes de la defomation, phase N
!!
!!      N+12*(N-1)+7+nvar1+...+nvarn-1, N+12*(N-1)+12+nvar1+...+nvarn-1
!!           Defn_old                    composantes de la deformation pas precedent, phase N
!!
!!      13*N+nvar1+...+nvarn-1+1,13*N+nvar1+...+nvarn
!!           Varintn                     variables internes, phaseN
!!
!!    ============================================================================
!!
!!   autres paramètres d'entrée : voir formalisme UMAT
!!
!! \param[out] Sig, Varint
!!
!!-------------------------------------------------------------------------
!! ATTENTION - ATTENTION - ATTENTION : notation de Voigt, 
!!            pour la deformation : x2 sur les termes de cisaillement 
!!            pour la contrainte  : pas de facteur 2
!!-------------------------------------------------------------------------
!!
!! Notation CAST3M (Voigt, ordre 11 22 33 12 13 23)
!=======================================================================


!=======================================================================
!>           FONCTION UMATLAMINATE
!!-----------------------------------------------------------------------
!! integration d'un modele composite LAMINATE a N phases
!!           
!! formalisme HPP
!!
!!-----------------------------------------------------------------------
!!
!! \param[in]   dt	Vecteur des pas de temps nouveau et ancien (dt_new,dt_old)
!! \param[in]   umatTab Tableau du type dérivé Umatptr contenant les pointeurs 
!!                      de fonction vers le comportement de chaque phase
!!
!!    ==========================================================================
!!      COEFFICIENTS
!!      ----------------------
!!
!!	Coeff(ncoeff) tableau de coefficients, ordre defini ci-dessous
!!                    avec ncoeff = 4*N+6+ncoeff1+.............+ncoeffn
!!
!!
!!	1    - nPhase                           nombre de phases dans le matériau composite
!!	2    - ncoeff1                          nombre de coeff de chaque loi
!!	3    - ncoeff2
!!             .......
!!	1+N  - ncoeffN
!!	
!!      2+N  - FV 1	                 	fractions volumiques de chaque phase (n-1 données la nieme vérifiant FV N = 1 - sum(FV I)
!!      3+N  - FV 2
!!             ....
!!      2N   - FV N-1
!!
!!      2N+1 - Nx                               Composantes du vecteur normal à l'interface
!!      2N+2 - Ny
!!      2N+3 - Nz
!!      2N+4 - Tx                               Composantes d'une direction tangeant à l'interface
!!      2N+5 - Ty
!!      2N+6 - Tz       
!!
!!      2N+7,2N+8
!!             (/ Lambdaeq, Mueq /)             coefficients de Lamé équivalents au comportement élastique initial du matériau 1
!!      2N+9,2N+8+ncoeff1
!!             Coeff1                           coefficients pour le materiau 1
!!      2N+9+ncoeff1,2N+10+ncoeff1
!!             (/ Lambdaeq, Mueq /)             coefficients de Lamé équivalents au comportement élastique initial du matériau 2
!!      2N+11+ncoeff1, 2N+10+ncoeff1+ncoeff2
!!             Coeff2                           coefficients pour le materiau 2
!!            ....
!!            ....
!!      2N+6+2*(N-1)+ncoeff1+.....+ncoeffn-1+1,4*N+6+ncoeff1+.............+ncoeffn-1
!!             (/ Lambdaeq, Mueq /)             coefficients de Lamé équivalents au comportement élastique initial du matériau N
!!      4*N+6+ncoeff1+.............+ncoeffn-1+1,4*N+6+ncoeff1+.............+ncoeffn
!!             CoeffN                           coefficients pour le materiau n
!!    ===========================================================================
!!
!!
!!    ===========================================================================
!!      VARIABLES INTERNES 
!!      -----------------------
!!
!!      Varint(nvarint) tableau de variables internes, ordre defini ci-dessous
!!                      avec nvarint= 10N + sum(nVari)
!!
!!      1  - nvar1            nombre de var. internes de chaque loi
!!      2  - nvar2
!!           .....
!!      N  - nvarN
!!
!!      1+N,3+N
!!           DefA1                      composantes de la defomation Anti-Plane
!!                                      dans la phase 1 (ordre NN,NT0,NT1, notation Voigt)
!!      4+N,6+N
!!           SigP1                      composantes de la contrainte Plane
!!                                      dans la phase 1 (ordre T0T0,T1T1,T0T1, notation Voigt)
!!      7+N,9+N
!!           DefA1_old                  composantes de la deformation Anti-Plane
!!                                      du pas de temps précédent  
!!      10+N,9+N+nvar1
!!           VarInt1                    variables internes phase 1
!!
!!      10+N+nvar1,12+N+nvar1
!!           DefA2                      composantes de la defomation Anti-Plane
!!                                      dans la phase 2 (ordre NN,NT0,NT1, notation Voigt)
!!      13+N+nvar1,15+N+nvar1
!!           SigP2                      composantes de la contrainte Plane
!!                                      dans la phase 2 (ordre T0T0,T1T1,T0T1, notation Voigt)
!!      16+N+nvar1,18+N+nvar1
!!           DefA2_old                  composantes de la deformation Anti-Plane
!!                                      du pas de temps précédent  
!!      19+N+nvar1,18+N+nvar1+nvar2
!!           VarInt2                    variables internes phase2
!!
!!       .........
!!       .........
!!       .........
!!
!!      N+9*(N-1)+1+nvar1+...+nvarn-1,N+9*(N-1)+3+nvar1+...+nvarn-1
!!           DefAn                      composantes de la defomation Anti-Plane
!!                                      dans la phase N (ordre NN,NT0,NT1, notation Voigt)
!!      N+9*(N-1)+4+nvar1+...+nvarn-1,N+9*(N-1)+6+nvar1+...+nvarn-1
!!           SigPn                      composantes de la contrainte Plane
!!                                      dans la phase N (ordre T0T0,T1T1,T0T1, notation Voigt)
!!      N+9*(N-1)+7+nvar1+...+nvarn-1, N+9*(N-1)+9+nvar1+...+nvarn-1
!!           DefAn_old                  composantes de la deformation Anti-Plane
!!                                      du pas de temps précédent  
!!      10*N+nvar1+...+nvarn-1+1,10*N+nvar1+...+nvarn
!!
!!    ============================================================================
!!
!!   autres paramètres d'entrée : voir formalisme UMAT
!!
!! \param[out] Sig, Varint
!!
!!-------------------------------------------------------------------------
!! ATTENTION - ATTENTION - ATTENTION : notation de Voigt, 
!!            pour la deformation : x2 sur les termes de cisaillement 
!!                                  soit ici NT0,NT1 et T0T1 
!!            pour la contrainte  : pas de facteur 2
!!-------------------------------------------------------------------------
!!
!! Notation CAST3M (Voigt, ordre 11 22 33 12 13 23)
!=======================================================================

!=======================================================================
!>           SUBROUTINE UMATCOMPOSITE
!!-----------------------------------------------------------------------
!! integration du modele composite a N phases
!!           
!!-----------------------------------------------------------------------
!!
!> Se referer aux en-tetes des fonctions umatReuss umatLaminate et umatVoigt
!! conservees ci-dessus pour avoir les descriptions des compositions 
!! des tableaux Coeff et VarInt.
!!
!!-----------------------------------------------------------------------
!!
!! \param[in]   numM Numero local du materiau correspondant au voxel composite en question
!! \param[in]   dt	Vecteur des pas de temps nouveau et ancien (dt_new,dt_old)
!! \param[in]   umatTab Tableau du type dérivé Umatptr contenant les pointeurs 
!!                      de fonction vers le comportement de chaque phase
!!
!!   autres paramètres d'entrée : voir formalisme UMAT
!!
!! \param[out] Sig, Varint
!!
!!-------------------------------------------------------------------------
!! ATTENTION - ATTENTION - ATTENTION : notation de Voigt, 
!!            pour la deformation : x2 sur les termes de cisaillement 
!!                                  soit ici NT0,NT1 et T0T1 
!!            pour la contrainte  : pas de facteur 2
!!-------------------------------------------------------------------------
!!
!! Notation CAST3M (Voigt, ordre 11 22 33 12 13 23)
!=======================================================================
subroutine umatComposite(numM,umatptr,sig,varint,def0,ddef,times,dt,&
                        nvarint,coeff,ncoeff,DFGRD0, DFGRD1,KINC,Cont)

  implicit none
  ! Arguments entrée-sortie 
  integer, intent(in)                                     :: numM
  real(mytype),dimension(ntens),intent(inout)             :: Sig
  real(mytype),dimension(ntens),intent(in)                :: Def0, DDef
  integer(kind=INT64),intent(in)                          :: ncoeff,nvarint
  real(mytype),dimension(nvarint),intent(inout),target    :: Varint
  real(mytype),dimension(ncoeff),intent(in),target        :: Coeff
  real(mytype),dimension(2),intent(in)                    :: times
  real(mytype),dimension(2),intent(in)                    :: dt      !(dt_new,dt_old)
  real(mytype), dimension(3,3),intent(in)                 :: DFGRD0, DFGRD1
  integer(kind=INT64),intent(inout)                       :: KINC
  integer(kind=INT32),intent(out)                         :: Cont  
  !----Pointeur de fonction
  type(c_ptr),dimension(int(Coeff(1))), intent(in)        :: umatptr  
  !----Variables de boucle
  integer                                                 :: i, j
  integer                                                 :: indVarint
  !----Variables erreur
  integer                                                 :: alloc_stat
  character(len=200)                                      :: err_msg
  !----Variables Composites
  !    le recours au pointeurs permet d'éviter les tableaux dynamiques de taille nPhases
  !    plus pénibles à utiliser
  integer(kind=INT32)                                     :: nPhases    
  integer(kind=INT64),dimension(:),allocatable            :: nVarintloc
  integer(kind=INT32)                                     :: taille_residu
  real(mytype),dimension(:),allocatable                   :: Fv,Residu,SecMembre,Residuold
  real(mytype),dimension(nvarint)                         :: VarInt0
  real(mytype),dimension(6)                               :: SigMoy
  real(mytype),dimension(:,:),allocatable                 :: DF 
  real(mytype),dimension(6,6)                             :: RotEps,iRotEps,RotSig,iRotSig  ! Pour le cas laminate
  real(mytype),dimension(6)                               :: DdefIn
  real(mytype),dimension(:,:),allocatable                 :: DdefLoc
  integer                                                 :: Contacv
  real(mytype)                                            :: Normresidu,NormSigma
  real(mytype), dimension(:,:),allocatable                :: ACT3_U,ACT3_R 
  logical                                                 :: erreur
  integer,dimension(:),allocatable                        :: ind_def

!modification_YC//////////////////////////////////////////////////////    
  integer :: stat
  real(mytype) :: criterion
!//////////////////////////////////////////////////////    

  ! Preparation pour le calcul du residu
  call init_taille_residu(numM,taille_residu,ind_def)

  !Initialisations bidons pour supprimer des gcc-warning 
  Normresidu = DFGRD0(1,1); Normresidu = DFGRD1(1,1);Normresidu = 0._mytype  

  ! Stockage des valeurs initiales des variables internes 
  Varint0    = Varint
  ! Initialisation des détecteurs d'erreurs 
  erreur     = .false.

  ! Nombre de Phases dans le matériau
  nPhases = int(Coeff(1),INT32)

  ! Allocation des tableaux Fv  
  allocate(Fv(nPhases),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (umatComposite)",2)
  Fv = 0.

  ! Allocation du tableau NVarInt
  allocate(nVarintloc(nPhases),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (umatComposite)",2)

  ! Nombre de varint de chaque phase
  nVarintloc = int(Varint(1:nPhases),kind=8)

  ! Allocation de la matrice tangente de l'algorithme de Newton-Raphson modifié
  allocate(DF(taille_residu*(nPhases-1),taille_residu*(nPhases-1)))
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (umatComposite)",2)
  DF = 0.

  ! Allocation du résidu et du second membre du système Newton-Raphson
  allocate(Residu(taille_residu*(nPhases-1)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (umatComposite)",2)
  Residu = 0.
  allocate(Residuold(taille_residu*(nPhases-1)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (umatComposite)",2)
  Residuold = 0.  !gcc-warning accepte
  allocate(SecMembre(taille_residu*(nPhases-1)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (umatComposite)",2)
  SecMembre = 0.

  ! Allocation du tableau des incréments de déformation de chaque phase
  allocate(DdefLoc(6,nPhases),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (umatComposite)",2)
  DdefLoc = 0.

  ! Allocation des tableaux mémoire pour l'accélération de convergence si nécessaire
  if (MattotP(numM)%acc_CV_composite) then
     allocate(ACT3_U(taille_residu*(nPhases-1),4),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (umatComposite)",2)
     allocate(ACT3_R(taille_residu*(nPhases-1),4),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (umatComposite)",2)
  end if

  ! Fractions volumiques de chaque phase
  Fv(1:nPhases-1) = Coeff(nPhases+2:2*nPhases) ! Fractions volumiques des n-1 premières phases passées en entrée
  Fv(nPhases)     = 1- sum(Fv(1:nPhases-1))    ! Fraction volumique de la dernière phase déduite des n-1 précédentes

  ! Test fractions volumiques nulles
  ! TODO placer ce test avant cette fonction (pas necesaire de le refaire a chaque passage)
  if (any(FV < 1e-4)) then
     write(err_msg,"(A)") "Une ou plusieurs fractions volumiques <1e-4  (umatComposite)"
     call amitex_abort(err_msg,2)
  end if


  !-------------------------------------------------------------------------------------------------
  !!                I - CALCULS PRELIMINAIRES ET INITIALISATIONS ALGO
  !-------------------------------------------------------------------------------------------------

  call init_newton(numM,taille_residu,nPhases,nVarintloc,Coeff,Varint,DDef,DdefLoc,&
                       dt,Fv,DF,DdefIn,RotEps,iRotEps,RotSig,iRotSig)
  
  !------------------------------------------------------------------------------------------------------
  !!                           II - BOUCLE NEWTON RAPHSON
  !------------------------------------------------------------------------------------------------------
  
  
  ! Compteurs d'itérations (global pour stopper si non convergence et autre  pour l'accélération de convergence acv)
  Cont    = 0
  Contacv = 0

  ! début boucle
  DO
     NormResidu = 0.
     Residu = 0.
     SigMoy(:)  = Sig

     ! CALCUL DU RESIDU
     call CalcResiduComposite(numM,umatptr,SigMoy,Def0,DdefLoc,ncoeff,nvarint,Varint,Varint0, &
                              Coeff,times,dt,Fv,Residu,KINC,erreur,RotEps, &
                              iRotEps,RotSig,iRotSig,DF)     
     ! test erreur
     if ((erreur) .and. (KINC /= -20)) then
        write(err_msg,"(A)") "Erreur lors de l'evaluation du residu (umatComposite)"
        call amitex_abort(err_msg,2)
        exit
     end if
     
     if  (KINC == -20) then
        exit
     end if

     if (trim(MattotP(numM)%Jacobian_type_composite) == 'numerical') then 

        ! CALCUL DE LA MATRICE TANGENTE NUMERIQUE
        !--------------------------------------------------------
!!! TODO : pouvoir appeler MatriceTangenteNumerique en reuss
        call MatriceTangenteNumerique(numM,umatptr,Sig,Def0,DdefLoc,ncoeff,nvarint,Varint0,&
             Coeff,times,dt,Fv,Residu,KINC,&
             RotEps,iRotEps,RotSig,iRotSig,erreur,DF)

        ! test erreur
        if (erreur) then
           write(err_msg,"(A)") &
                         "Erreur lors du calcul de la matrice tangente numerique (umatComposite)"
           call amitex_abort(err_msg,2)
           exit
        end if

     end if


     ! Passage de l'incrément de déformation dans le repère de l'interface pour laminate
     !--------------------------------------------------------------------
     if (MattotP(numM)%lawname=="laminate") DdefLoc = matmul(RotEps,DdefLoc)

     ! Stockage des incréments et résidus pour accélération de convergence
     if (MattotP(numM)%acc_CV_composite) then
        Contacv = Contacv + 1
        if (Contacv .eq. 5) Contacv = 1
        ACT3_R(:,Contacv) = Residu
        ACT3_U(:,Contacv) = reshape( DdefLoc(ind_def,2:nPhases), (/ taille_residu*(nPhases-1) /) )
     end if

     ! Evaluation de la norme de la contrainte moyenne et de celle du résidu
     !----------------------------------------------------------------------
     NormSigma  = sqrt(dot_product(Sigmoy,Sigmoy))
     NormResidu = sqrt(dot_product(Residu,Residu))


     ! TEST CONVERGENCE
     !-----------------
!modification_yangChen//////////////////////////////////////////////////////    
     !IF (Normresidu < MattotP(numM)%tol_criteq_composite*algo_param%tol_criteq*NormSigma) then  !original

     criterion = MattotP(numM)%tol_criteq_composite*algo_param%tol_criteq*NormSigma !YC:2022.03.19
     if (criterion<1.e-4)  criterion = 1.e-4
     IF (Normresidu < criterion) then
!///////////////////////////////////////////////////////////////////////////////////////////
        ! Convergence atteinte 
        ! print *,"#Cont = ", Cont
        ! print *, "normresidu/normSigma, cont", Normresidu/NormSigma, cont, MattotP(numM)%tol_criteq_composite*algo_param%tol_criteq
        EXIT
     END IF

     ! RESOLUTION SYSTEME :
     !-------------------------------------------
     if (mattotP(numM)%acc_CV_composite .AND. (Cont > 3) .AND.  (modulo(Cont-4,3) .EQ. 0) ) then
         ! Calcul de l'incrément suivant par la méthode d'accélération de convergence
        call ACT3LOCAL(ACT3_U,ACT3_R,SecMembre)
        DdefLoc(ind_def,2:nPhases) = reshape(SecMembre, (/ taille_residu, int(nPhases-1) /) )
     else
        ! Résolution du système linéaire Newton-Raphson modifié
        SecMembre = -Residu
        call lusolve(DF,SecMembre)    ! résolution du système
        ! Actualisation de l'incrément de déformation
        !--------------------------------------------
        DdefLoc(ind_def,2:nPhases) = DdefLoc(ind_def,2:nPhases) + &
                                     reshape(SecMembre, (/ taille_residu, int(nPhases-1) /) )
     end if
     ! Calcul de l'incrément de déformation de la phase 1 permettant de vérifier la déformation moyenne imposée
     DdefLoc(ind_def,1) = (1/FV(1)) * DdefIn(ind_def)
     do i = 2,nPhases
        DdefLoc(ind_def,1) = Ddefloc(ind_def,1) - (FV(i)/FV(1))*Ddefloc(ind_def,i)
     end do

     ! Pour laminate, retour dans le repère global de l'incrément de déformation pour calcul du comportement suivant
     !------------------------------------------------------------------------------------------------
     if (MattotP(numM)%lawname=="laminate") DdefLoc = matmul(iRotEps,DdefLoc)
     
     ! Incrémentation du compteur
     Cont = Cont + 1
     
     ! Si on dépasse 200 itérations => NON CONVERGENCE : on sort de la boucle et on arrête le calcul
     ! Marqueur caractéristique KINC = -10
!/////////////////////////////////////////////////////////////////////
!modification YC, 2022.03.29
     !if (Cont >= 200) then
     !   print *, "(umatComposite) Cont. max reached (200), normresidu/normSigma =", Normresidu/normSigma 
     !   KINC = -10
     !   EXIT
     !end if
     if (Cont >= 1000) then
         if (nrank==0) print *,'Cont reaches 1000, --> exit, but no KINC=-10'
         if (nrank==0) print *,'          Normresidu: ', Normresidu
         EXIT
     end if
!/////////////////////////////////////////////////////////////////////
     
  END DO

    
     !!!!! ATTENTION : après convergence, l'incrément de déformation solution Ddefloc est exprimé dans le repère
     !!!!!             locale lié à l'interface (N,T1,T2,NT1,NT2,T1T2) pour LAMINATE


  !------------------------------------------------------------------------------------------------------
  !!                              III - POST TRAITEMENTS
  !------------------------------------------------------------------------------------------------------ 
  
  ! Mise à jour de la contrainte moyenne
  Sig = Sigmoy

     
  ! Mise à jour des deformations de chaque phase (deformation anti-plane pour laminate)
  do i=1,nPhases
     indVarint     =  nPhases+((6+taille_residu)*(i-1))   + int(sum(nVarintloc(1:i-1)),INT32)
     if (indVarint+6+taille_residu+nVarintloc(i)>nvarint) then
        write(err_msg,"(A)") "Indices hors des bornes varint (umatComposite)"
        call amitex_abort(err_msg,2)
     end if
     do j = 1, taille_residu
        Varint(indVarint+j) = Varint0(indVarint+j) + DdefLoc(ind_def(j),i) 
     end do
  end do

       
  ! Desallocation des tableaux NVarIntloc, Fv  
  deallocate(nVarintloc)
  deallocate(Fv)

  ! Desallocation de la matrice tangente de l'algorithme de Newton-Raphson modifié
  deallocate(DF)

  ! Desallocation du résidu et du second membre du système Newton-Raphson
  deallocate(Residu)
  deallocate(SecMembre)

  ! Desallocation du tableau des incréments de déformation de chaque phase
  deallocate(DdefLoc)

  ! Desallocation des tableaux mémoire pour l'accélération de convergence si nécessaire
  if (MattotP(numM)%acc_CV_composite) then
     deallocate(ACT3_U)
     deallocate(ACT3_R)
  end if

end subroutine umatComposite

!=======================================================================
!>           SUBROUTINE UMATCOMPOSITE_EXT
!!-----------------------------------------------------------------------
!! APPLICATION DES LOIS UMAT "EXTENSION" DANS CHACUNE DES PHASES
!!
!!       ATTENTION : suppose que les champs varint0,def+def et sig 
!!                   sont "equilibres" lors de l'appel a behavior
!!
!! UTILISATION VISEE : modeles non-locaux
!! on appelle UMAT "extension" pour evaluer une nouvelle variable interne
!!           
!! ENTREES : voir UMATCOMPOSITE
!!
!! FONCTIONNEMENT : ajout d'arguments optionels (umat_only, DDef)
!!                  a calcResiduComposite (utilise dans umatComposite)                 
!!
!=======================================================================
subroutine umatComposite_ext(numM,umatptr,sig,varint,def0,ddef,times,dt,&
                        nvarint,coeff,ncoeff,DFGRD0, DFGRD1,KINC)

  implicit none
  ! Arguments entrée-sortie 
  integer, intent(in)                                     :: numM
  real(mytype),dimension(ntens),intent(inout)             :: Sig
  real(mytype),dimension(ntens),intent(in)                :: Def0, DDef
  integer(kind=INT64),intent(in)                          :: ncoeff,nvarint
  real(mytype),dimension(nvarint),intent(inout),target    :: Varint
  real(mytype),dimension(ncoeff),intent(in),target        :: Coeff
  real(mytype),dimension(2),intent(in)                    :: times
  real(mytype),dimension(2),intent(in)                    :: dt      !(dt_new,dt_old)
  real(mytype), dimension(3,3),intent(in)                 :: DFGRD0, DFGRD1
  integer(kind=INT64),intent(inout)                       :: KINC
  !----Pointeur de fonction
  type(c_ptr),dimension(int(Coeff(1))), intent(in)        :: umatptr  

 
  logical                                                 :: umat_only = .true.
  integer                                                 :: nPhases
 
  ! ADHERENCE a CalcResiduComposite
  real(mytype),dimension(1,1)                             :: DdefLoc !non utilise ici
  real(mytype),dimension(nvarint)                         :: VarInt0
  real(mytype),dimension(:),allocatable                   :: Fv
  real(mytype),dimension(1)                               :: Residu  !non utilise ici
  logical                                                 :: erreur
  real(mytype),dimension(1,1)                             :: DF      !non utilise ici 
  real(mytype),dimension(6,6)                             :: RotEps,iRotEps,RotSig,iRotSig  ! Pour le cas laminate
  
  integer                                                 :: alloc_stat 
  character(len=200)                                      :: err_msg


  !Initialisations bidons pour supprimer des gcc-warning 
  Residu(1) = DFGRD0(1,1); Residu(1) = DFGRD1(1,1)  

  ! NOMBRE DE PHASES
  nPhases = int(Coeff(1),INT32)

  allocate(Fv(nPhases),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (umatComposite)",2)
  Fv = 0.
  ! Fractions volumiques de chaque phase
  Fv(1:nPhases-1) = Coeff(nPhases+2:2*nPhases) ! Fractions volumiques des n-1 premières phases passées en entrée
  Fv(nPhases)     = 1- sum(Fv(1:nPhases-1))    ! Fraction volumique de la dernière phase déduite des n-1 précédentes

  ! Test fractions volumiques nulles
  ! TODO placer ce test avant cette fonction (pas necesaire de le refaire a chaque passage)
  if (any(FV < 1e-4)) then
     write(err_msg,"(A)") "Une ou plusieurs fractions volumiques <1e-4  (umatComposite)"
     call amitex_abort(err_msg,2)
  end if
  
  ! ADHERENCE a CalcResiduComposite   
  DdefLoc = 0.
  Residu = 0.
  DF = 0.
  
  Varint0    = Varint
  erreur     = .false.

 
  ! Appel a UMAT "extension" via parametres optionnels de CalcResiduComposite
  call CalcResiduComposite(numM,umatptr,Sig,Def0,DdefLoc,ncoeff,nvarint,Varint,Varint0, &
                              Coeff,times,dt,Fv,Residu,KINC,erreur,RotEps, &
                              iRotEps,RotSig,iRotSig,DF,&
                              umat_only,DDef) ! arguments optionels pour UMAT "extension" (pas de calcul de residu)

end subroutine umatComposite_ext
!==============================================================================
!
!                      SUBROUTINE INIT_TAILLES_RESIDU
!
!> Initialise la taille du residu et les indices des composantes de la deformation
!! qui interviennent dans le residu.
!!
!!
!! \param[in]  numM          : le numero (local) du materiau composite en question
!! 
!! \param[out] taille_residu : la taille du residu (6 pour reuss, 3 pour laminate 
!!                             et normalement 0 pour voigt mais on met 1 pour les 
!!                             raisons ecrites dans les commentaires ci-dessous)
!!
!! \param[out] ind_def       : les composantes de la deformation qui apparaissent
!!                             dans le residu
!!
!==============================================================================


subroutine init_taille_residu(numM,taille_residu,ind_def)

implicit none


  integer,intent(in)                                       :: numM
  integer,intent(out)                                      :: taille_residu
  integer,dimension(:),allocatable,intent(out)             :: ind_def
  integer                                                  :: alloc_stat

  if (MattotP(numM)%lawname=="voigt") then              ! En realite, le residu fait une taille 0 pour voigt
     taille_residu = 1                                  ! car les contraintes sont egales par definition du 
     allocate(ind_def(6),stat = alloc_stat)             ! modele. On alloue cependant une taille 1 au residu
     if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (umatComposite)",2)
     ind_def = (/0/)                                    ! On alloue cependant une taille 1 au residu pour 
                                                        ! pouvoir le mettre egal a 0 (des la premiere 
                                                        ! iteration) et quoi qu'il arrive, on sortira au
                                                        ! premier test de convergence.

  elseif (MattotP(numM)%lawname=="reuss") then
     taille_residu = 6
     allocate(ind_def(6),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (umatComposite)",2)
     ind_def = (/1,2,3,4,5,6/)


  elseif (MattotP(numM)%lawname=="laminate") then
     taille_residu = 3
     allocate(ind_def(3),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (umatComposite)",2)
     ind_def = IndAP
  else
     allocate(ind_def(1)); taille_residu = 0 ! supprime un gcc-warning
     call amitex_abort("composite behavior name different from Voigt, Reuss or Laminate",2)
  end if


end subroutine init_taille_residu



!==============================================================================
!
!                      SUBROUTINE INIT_NEWTON
!
!> Initialise l'increment de deformation et calcule la matrice tangente elastique
!! lorsque c'est necessaire.
!!
!!
!!
!==============================================================================


subroutine init_newton(numM,taille_residu,nPhases,nVarintloc,Coeff,Varint,DDef,DdefLoc,&
                       dt,Fv,DF,DdefIn,RotEps,iRotEps,RotSig,iRotSig)

implicit none


  integer,intent(in)                                      :: numM ! numero local du materiau
  integer,intent(inout)                                   :: taille_residu 
  integer,intent(in)                                      :: nPhases
  integer(kind=INT64),dimension(:),intent(in)             :: nVarintloc !nb de varint de chaque phase
  real(mytype),dimension(:),intent(in),target             :: Coeff
  real(mytype),dimension(:),intent(in),target             :: Varint
  real(mytype),dimension(ntens),intent(in)                :: DDef
  real(mytype),dimension(:,:),allocatable                 :: DdefLoc
  real(mytype),dimension(2),intent(in)                    :: dt 
  real(mytype),dimension(:),intent(in)                    :: Fv ! fraction volumique de chaque phase
  real(mytype),dimension(:,:),intent(inout)               :: DF 
  real(mytype),dimension(6),intent(out)                   :: DdefIn
  real(mytype),dimension(6,6),intent(out)                 :: RotEps,iRotEps,RotSig,iRotSig  ! Pour le cas laminate
  integer                                                 :: indCoeff, i


  if (MattotP(numM)%lawname=="laminate") then
     indCoeff = 2*nPhases
     ! Calcul des matrices de changement de bases en notation Voigt 
     call CalcRotMatrices(Coeff(indCoeff+1:indCoeff+3),&
                       Coeff(indCoeff+4:indCoeff+6),RotEps,IRotEps,RotSig,IRotSig) 

     ! Construction de la matrice tangente approchée pour l'algorithme de Newton-Raphson 
     ! (si demande)
     !----------------------------------------------------------------------------------
     if (trim(Mattotp(numM)%Jacobian_type_composite) == 'elastic') then 
        call MatriceTangenteElastiqueLaminate(Coeff,DF,FV)
     end if

     ! Incrément de déformation moyen -- repère interface
     DdefIn        =  matmul(RotEps,ddef)

     ! Incrément de déformation de chaque phase -- DdefLoc dans le repère global
     call InitDeltaEpsLaminate(DdefLoc,Varint,trim(algo_param%Init_laminate),DdefIn,FV,nVarintloc,dt,iRotEps)

  elseif (MattotP(numM)%lawname=="reuss") then
     ! Construction de la matrice tangente approchée pour l'algorithme de Newton-Raphson 
     ! (si demande)
     !----------------------------------------------------------------------------------
     if (trim(Mattotp(numM)%Jacobian_type_composite) == 'elastic') then 
        call MatriceTangenteElastiqueReuss(Coeff,DF,FV)
     end if

     DdefIn = ddef

     ! Increment de deformation de chaque phase --> DdefLoc 
     call InitDeltaEpsReuss(DdefLoc,Varint,algo_param%Init_reuss,DdefIn,FV,nVarintloc,dt)

  elseif (MattotP(numM)%lawname=="voigt") then
     ! Pas de matrice tangente pour voigt car on n'a pas de systeme a resoudre
     DdefIn = ddef
     do i = 1, size(DdefLoc(1,:))
        DdefLoc(:,i)=ddef
     end do
     taille_residu = 0 ! on remet la vraie taille pour la suite, 'taille_residu = 1' servait uniquement 
                       ! a allouer les differents tableaux avec une taille non nulle
  end if     
  

end subroutine init_newton

!==============================================================================
!
!                      SUBROUTINE PILOTAGECOMPOSITE
!
!> Calcul du comportement pour un voxel composite multi-couche
!!
!!  1 - appel à umatComposite
!!  2 - en cas de non convergence de l'algorithme d'homogeneisation/ d'une loi locale à une phase du voxel :
!!            -> division du pas de temps en deux sous pas de temps égaux et calcul du comportement
!!            -> subdivisions successives jusqu'à convergence du calcul 
!!            -> nombre de sudivisions maximum autorisé fixé par l'utilisateur (5 par défaut)  --> algo_param%Nmax_subdivision
!!
!!  ATTENTION : * ne décompose pour le moment que les incréments de temps et de déformation
!!                  ==> les paramètres externes sont supposés constants sur le pas de temps 
!!                      ayant leur valeur à t+dt
!!  
!!  \param[out] Cont : nombre d'itération moyen dans umatComposite
!!  \param[out] Nincrements : nombre de sous pas de temps utilisés pour intégrer le modèle multi-couche

!!
!!  Pour le reste : Même entrées que umatComposite --> cf subroutine umatComposite pour détails
!!
!! Modifie egalement les variables internes dans MattotP
!!
!==============================================================================
subroutine  PilotageComposite(numM,umatptr,sig,varint,def0,ddef,times,dt,&
                        nvarint,coeff,ncoeff,DFGRD0, DFGRD1,KINC,Cont,Nincrements)

  implicit none
  ! Arguments entrée-sortie 
  integer, intent(in)                                    :: numM ! numero local du materiau composite
  real(mytype),dimension(ntens),intent(inout)            :: Sig
  real(mytype),dimension(ntens),intent(in)               :: Def0, DDef
  integer(kind=8),intent(in)                             :: ncoeff,nvarint
  real(mytype),dimension(nvarint),intent(inout)          :: Varint
  real(mytype),dimension(ncoeff),intent(in)              :: Coeff
  real(mytype),dimension(2),intent(in)                   :: times
  real(mytype),dimension(2),intent(in)                   :: dt
  real(mytype), dimension(3,3),intent(in)                :: DFGRD0, DFGRD1
  integer(kind=8),intent(inout)                          :: KINC
  integer(kind=INT32),intent(out)                        :: Cont
  integer(kind=INT32),intent(out)                        :: Nincrements
  !----Pointeur de fonction
  type(c_ptr),dimension(int(Coeff(1))), intent(in)        :: umatptr 
  !----Variables de pilotage
  real(mytype),dimension(ntens)                          :: Sig0
  real(mytype),dimension(nvarint)                        :: Varint0
  real(mytype),dimension(ntens)                          :: def00, dddef
  real(mytype),dimension(2)                              :: ddt, times0
  real(mytype), dimension(3,3)                           :: DDFGRD0, DDFGRD1
  !----Variables de boucles
  integer(kind=4)                                        :: i,j  
  !----Variables erreur
  character(len=200)                                     :: err_msg
  !----Autres
  integer(kind=INT32)                                    :: nIt  

  ! Initialisations
  !----------------
  ddt =  DFGRD0(1,1:2); ddt = DFGRD1(1,1:2); ddt = [0.,0.] ! initialisations bidons pour 
                                                           ! supprimer des gcc-warning
  Cont    = 0
  nIt     = 0
  Sig0    = Sig
  Varint0 = Varint
  Nincrements = 0  
 
  ! Premier appel au comportement
  !-------------------------------------------------------------------
  call  umatComposite(numM,umatptr,sig,varint,def0,ddef,times,dt,&
       nvarint,coeff,ncoeff,DFGRD0, DFGRD1,KINC,Cont)
  
  
  ! En cas de non convergence : subdivision du pas de temps
  !------------------------------------------------------------
  do j=1,MattotP(numM)%Nmax_subdivision_composite

     if (KINC == 1) then
        exit
     end if

     ! Nombre d'incréments de subdivision : on divise par 2 tant que le calcul ne converge pas
     !                                      jusqu'à atteindre le nombre max d'itérations
     Nincrements = 2**j     
     ! Contrainte moyenne et variables internes reinitialisées
     Sig    = Sig0
     Varint = Varint0
     Cont = 0
     ! Indice d'erreur remis réinitialisé
     KINC = 1 
     ! Nouveaux incréments  de temps et de déformation:
     ddt(1) = dt(1)/Nincrements
     dddef  = DDef /Nincrements  
     ! Déformation initiale et instant initial réinitialisés
     def00  = Def0
     times0 = times

     ! Mise à jour des tenseurs gradient de transformation pour passage loi grandes def en HPP
     DDFGRD0 = DFGRD0
     !DDFGRD1
     DDFGRD1(1,1) = dble(1.+def00(1)+dddef(1))
     DDFGRD1(2,1) = dble(def00(4)+dddef(4))/2._mytype
     DDFGRD1(3,1) = dble(def00(5)+dddef(5))/2._mytype
     DDFGRD1(1,2) = dble(def00(4)+dddef(4))/2._mytype
     DDFGRD1(2,2) = dble(1.+def00(2)+dddef(2))
     DDFGRD1(3,2) = dble(def00(6)+dddef(6))/2._mytype
     DDFGRD1(1,3) = dble(def00(5)+dddef(5))/2._mytype
     DDFGRD1(2,3) = dble(def00(6)+dddef(6))/2._mytype
     DDFGRD1(3,3) = dble(1.+def00(3)+dddef(3))

     ! Appels successifs à umatComposite pour chaque sous-pas de temps
     do i = 1,Nincrements

        if (i > 1) then
           ddt(2) = ddt(1)
        end if

        call  umatComposite(numM,umatptr,sig,varint,def00,dddef,times,ddt,&
             nvarint,coeff,ncoeff,DDFGRD0, DDFGRD1,KINC,nIt)

        Cont = Cont + nIt
        nIt  = 0

        if ((KINC == -10) .or. (KINC == -20)) then
           ! sortie de boucle et resubdivision si non convergence
           exit
        end if

        ! Mise à jour valeurs initiales
        def00 = def00 + dddef
        times0(2) = times0(2) + ddt(1)

        !! Mise à jour des tenseurs gradient de transformation initiaux et finaux
        !DFGRD0
        DDFGRD0(1,1) = dble(1.+def00(1))
        DDFGRD0(2,1) = dble(def00(4))/2._mytype
        DDFGRD0(3,1) = dble(def00(5))/2._mytype
        DDFGRD0(1,2) = dble(def00(4))/2._mytype
        DDFGRD0(2,2) = dble(1.+def00(2))
        DDFGRD0(3,2) = dble(def00(6))/2._mytype
        DDFGRD0(1,3) = dble(def00(5))/2._mytype
        DDFGRD0(2,3) = dble(def00(6))/2._mytype
        DDFGRD0(3,3) = dble(1.+def00(3))
        !DDFGRD1
        DDFGRD1(1,1) = dble(1.+def00(1)+dddef(1))
        DDFGRD1(2,1) = dble(def00(4)+dddef(4))/2._mytype
        DDFGRD1(3,1) = dble(def00(5)+dddef(5))/2._mytype
        DDFGRD1(1,2) = dble(def00(4)+dddef(4))/2._mytype
        DDFGRD1(2,2) = dble(1.+def00(2)+dddef(2))
        DDFGRD1(3,2) = dble(def00(6)+dddef(6))/2._mytype
        DDFGRD1(1,3) = dble(def00(5)+dddef(5))/2._mytype
        DDFGRD1(2,3) = dble(def00(6)+dddef(6))/2._mytype
        DDFGRD1(3,3) = dble(1.+def00(3)+dddef(3))


     end do

     ! Calcul du nombre d'itération moyen par appel dans l'algorithme
     Cont = Cont/Nincrements
     
  end do

  if (KINC == -10) then
     write(err_msg,"(A,I4,A)") "Non convergence umatComposite pour une subdivison en ",Nincrements,& 
       " sous pas de temps -> subdivision maximale autorisée : arrêt du calcul (Pilotage umatComposite)"
     call amitex_abort(err_msg,2)
  elseif (KINC == -20) then
     write(err_msg,"(A,I4,A)") "Non convergence umat locale pour une subdivision en ",Nincrements,&
       " sous pas de temps -> subdivision maximale autorisée :  arrêt du calcul (Pilotage umatComposite)"
     call amitex_abort(err_msg,2)
  end if

end subroutine PilotageComposite

!==============================================================================
!      SUBROUTINE DEALLOCATE_MATCOMPOSITE
!
!>  DESALLOCATION DU TABLEAU DE MATERIAUX
!
!------------------------------------------------------------------------------
subroutine deallocate_MatComposite()

  implicit none

  integer :: i

  do i = 1, size(MatComposite)
    if(allocated(MatComposite(i)%num_Phases)) deallocate(MatComposite(i)%num_Phases)
!    if(allocated(MatComposite(i)%Zone)) deallocate(MatComposite(i)%Zone)
!    if(allocated(MatComposite(i)%pos_globale)) deallocate(MatComposite(i)%pos_globale)
!    if(allocated(MatComposite(i)%Fv)) deallocate(MatComposite(i)%Fv)
!    if(allocated(MatComposite(i)%Normale)) deallocate(MatComposite(i)%Normale)
!    if(allocated(MatComposite(i)%Direction)) deallocate(MatComposite(i)%Direction)
!    if(allocated(MatComposite(i)%Surfaces)) deallocate(MatComposite(i)%Surfaces)
  end do

  deallocate(MatComposite)
end subroutine deallocate_MatComposite

!==============================================================================
!      SUBROUTINE DEALLOCATE_INTERPHASES
!
!>  DESALLOCATION DU TABLEAU INTERPHASES
!
!------------------------------------------------------------------------------
subroutine deallocate_Interphases()

  implicit none

  integer :: i

  do i = 1, size(Interphases)
    if(allocated(Interphases(i)%zones)) deallocate(Interphases(i)%zones)
  end do

  deallocate(Interphases)
end subroutine deallocate_Interphases

!======================================================================
!                     Procédures propres umatLaminate
!
!  Routines utilisées pour intégrer le modèle multi-couche 
!  NON PUBLIQUES
!
! Calc_rot_matrices
! InitDeltaEpsLaminate
! CalcResiduLaminate
! MatriceTangenteElastiqueLaminate
! MatriceTangenteNumerique
!======================================================================

!!===================================================================
!!                 CALC_ROT_MATRICES                                !
!!                                                                  !
!! Calcul des matrices de changement de base en notation de Voigt   !
!! pour les tenseurs d'ordre 2 symmetriques                         !
!!                          HPP                                     !
!! --------------------------------------------------------------   !
!! Epsinterface = RotEps*Epsglobal  Epsglobal = iRotEps*Epsinterface!
!! Siginterface = RotSig*Sigglobal  Sigglobal = iRotSig*Siginterface!
!!
!! \param[in] N,T : vecteurs normal/tangent à l'interface entre les !
!!                  couches du laminé                               !
!! \param[ou] RotEps,iRotEps,RotSig,iRotSig : matrices de changement!
!!                                            de base               !
!!                                                                  !
!!===================================================================
subroutine CalcRotMatrices(N,T,RotEps,IRotEps,RotSig,IRotSig)

  implicit none

  real(mytype), dimension(3),intent(in)       :: N,T  ! Coordonnées des vecteurs normaux et tangents à l'interface
  real(mytype), dimension(3,3)                :: R,invR    
  real(mytype), dimension(6,6),intent(out)    :: RotEps,iRotEps,RotSig,iRotSig
  real(mytype),dimension(6,6)                 :: iRotEps2,iRotSig2
  integer                                     :: i,j

  !! R = [N T T2] matrice de passage du repère global au repère de l'interface
  !! T2 = cross(N,T) (seconde direction tangente à l'interface pour compléter le repère)

  R(:,1) = N
  R(:,2) = T
  R(1,3) = (N(2)*T(3)) - (N(3)*T(2))
  R(2,3) = (N(3)*T(1)) - (N(1)*T(3))
  R(3,3) = (N(1)*T(2)) - (N(2)*T(1))

  call InvLu(R,invR)

  do i = 1,3,1

     do j = 1,3,1
        RotEps(i,j) = R(j,i)**2;
     end do

     RotEps(i,4) = R(1,i)*R(2,i); 
     RotEps(i,5) = R(1,i)*R(3,i);
     RotEps(i,6) = R(2,i)*R(3,i);

     RotEps(4,i) = 2*R(i,1)*R(i,2);
     RotEps(5,i) = 2*R(i,1)*R(i,3);
     RotEps(6,i) = 2*R(i,2)*R(i,3);
  end do


  RotEps(4,4) = R(1,1)*R(2,2) + R(2,1)*R(1,2);
  RotEps(4,5) = R(1,1)*R(3,2) + R(3,1)*R(1,2);
  RotEps(4,6) = R(2,1)*R(3,2) + R(3,1)*R(2,2);
  RotEps(5,4) = R(1,1)*R(2,3) + R(2,1)*R(1,3);
  RotEps(5,5) = R(1,1)*R(3,3) + R(3,1)*R(1,3);
  RotEps(5,6) = R(2,1)*R(3,3) + R(3,1)*R(2,3);
  RotEps(6,4) = R(1,2)*R(2,3) + R(2,2)*R(1,3);
  RotEps(6,5) = R(1,2)*R(3,3) + R(3,2)*R(1,3);
  RotEps(6,6) = R(2,2)*R(3,3) + R(3,2)*R(2,3);

   do i = 1,3,1

     do j = 1,3,1
        iRotEps2(i,j) = invR(j,i)**2;
     end do

     iRotEps2(i,4) = invR(1,i)*invR(2,i); 
     iRotEps2(i,5) = invR(1,i)*invR(3,i);
     iRotEps2(i,6) = invR(2,i)*invR(3,i);

     iRotEps2(4,i) = 2*invR(i,1)*invR(i,2);
     iRotEps2(5,i) = 2*invR(i,1)*invR(i,3);
     iRotEps2(6,i) = 2*invR(i,2)*invR(i,3);
  end do


  iRotEps2(4,4) = invR(1,1)*invR(2,2) + invR(2,1)*invR(1,2);
  iRotEps2(4,5) = invR(1,1)*invR(3,2) + invR(3,1)*invR(1,2);
  iRotEps2(4,6) = invR(2,1)*invR(3,2) + invR(3,1)*invR(2,2);
  iRotEps2(5,4) = invR(1,1)*invR(2,3) + invR(2,1)*invR(1,3);
  iRotEps2(5,5) = invR(1,1)*invR(3,3) + invR(3,1)*invR(1,3);
  iRotEps2(5,6) = invR(2,1)*invR(3,3) + invR(3,1)*invR(2,3);
  iRotEps2(6,4) = invR(1,2)*invR(2,3) + invR(2,2)*invR(1,3);
  iRotEps2(6,5) = invR(1,2)*invR(3,3) + invR(3,2)*invR(1,3);
  iRotEps2(6,6) = invR(2,2)*invR(3,3) + invR(3,2)*invR(2,3); 

  RotSig = RotEps
  RotSig(1:3,4:6) = 2*RotSig(1:3,4:6)
  RotSig(4:6,1:3) = 0.5*RotSig(4:6,1:3)

  iRotSig2 = iRotEps2
  iRotSig2(1:3,4:6) = 2*iRotSig2(1:3,4:6)
  iRotSig2(4:6,1:3) = 0.5*iRotSig2(4:6,1:3)

  iRotEps = iRotEps2
  iRotSig = iRotSig2

end subroutine CalcRotMatrices



!!===================================================================
!!                 CALC_ROT_MATRICES_GD                             !
!!                                                                  !
!! Calcul des matrices de changement de base en notation vectorielle!
!! pour les tenseurs d'ordre 2                                      !
!!                          GD                                      !
!! --------------------------------------------------------------   !
!! T : tenseur dans l'ancienne base en notation vectorielle         !
!! T': tenseur dans la nouvelle base definie par V1 et V2 en notation
!!     vectorielle                                                  !
!! T' = Rot*T  et T = iRot*T'                                       !
!!                                                                  !   
!!                                                                  !    
!! \param[in]  V1,V2 : 2 premiers vecteurs de la nouvelle base       !
!! \param[out] Rot,iRot : matrices de changement de base  (9,9)     !
!!                   Notation vectorielle 11,22,33,12,13,23,21,31,32!                                             !
!!===================================================================
subroutine CalcRotMatrices_GD(V1,V2,Rot,iRot)

  implicit none

  real(mytype), dimension(3),intent(in)       :: V1,V2  ! Coordonnées des vecteurs de la nouvelle base dans la base globale
  real(mytype), dimension(3,3)                :: R,invR    
  real(mytype), dimension(9,9),intent(out)    :: Rot,iRot
  integer                                     :: i,j   ! indices notation vectorielle
  integer,dimension(2)                        :: ind1,ind2 ! indices matrices de rotation

  !! R = [V1 V2 V3] matrice de passage du repère global au repère de l'interface
  !! V3 = cross(V1,V2) (seconde direction tangente à l'interface pour compléter le repère)

  R(:,1) = V1
  R(:,2) = V2
  R(1,3) = (V1(2)*V2(3)) - (V1(3)*V2(2))
  R(2,3) = (V1(3)*V2(1)) - (V1(1)*V2(3))
  R(3,3) = (V1(1)*V2(2)) - (V1(2)*V2(1))

  call InvLu(R,invR)

  do i = 1,9
     do j = 1,9
        ind1 = lin_indices_t2(i);
        ind2 = lin_indices_t2(j);
        Rot(i,j)  = R(ind2(1),ind1(1))*R(ind2(2),ind1(2));
        iRot(i,j) = invR(ind2(1),ind1(1))*invR(ind2(2),ind1(2));
     end do
  end do

end subroutine CalcRotMatrices_GD

!==============================================================================
!       FUNCTION LIN_INDICES_T2
!------------------------------------------------------------------------------
!> Renvoie les deux indices de composantes d'un tenseur sous sa forme matricielle
!  à partir des indices sous sa forme vectorielle. Pour les tenseurs non 
!   symmétriques uniquement
!
!  forme matricielle T = T11 T12 T13  =  T1 T4 T5
!                        T21 T22 T23     T7 T2 T6
!                        T31 T32 T33     T8 T9 T3
!
!  forme vectorielle T = T1 T2 T3 T4 T5 T6 T7 T8 T9
!
!  Correspondance : 1=11 2=22 3=33 4=12 5=13 6=23 7=21 8=31 9=32
!
!  \param[in]  :  i  indices de la composante sous forme vectorielle
!
!  \param[out]       indices de la composante sous forme matricielle
!==============================================================================
function lin_indices_t2(i)

  implicit none
  integer,intent(in)     :: i
  integer,dimension(2)   :: lin_indices_t2

  select case(i)
    case(1)
        lin_indices_t2(1) = 1
        lin_indices_t2(2) = 1
    case(2)
        lin_indices_t2(1) = 2
        lin_indices_t2(2) = 2
    case(3)
        lin_indices_t2(1) = 3
        lin_indices_t2(2) = 3
    case(4)
        lin_indices_t2(1) = 1
        lin_indices_t2(2) = 2
    case(5)
        lin_indices_t2(1) = 1
        lin_indices_t2(2) = 3
    case(6)
        lin_indices_t2(1) = 2
        lin_indices_t2(2) = 3
    case(7)
        lin_indices_t2(1) = 2
        lin_indices_t2(2) = 1
    case(8)
        lin_indices_t2(1) = 3
        lin_indices_t2(2) = 1
    case(9)
        lin_indices_t2(1) = 3
        lin_indices_t2(2) = 2
  end select

end function lin_indices_t2

!!===================================================================
!!                      INIT_DELTA_EPS_Laminate                              !
!!                      --------------                              !
!!                                                                  !
!! Initialisation de l'incrément de déformation pour l'algorithme   !
!! de Newton-Rapshon de umatLaminate. Trois types d'initialisations !
!!                                                                  !
!! Initialisation par défaut :                                      !
!! --------------------------                                       !
!! Homogène dans chaque phase égale à l'DEpsMoyen imposé            !
!!                                                                  !
!!                                                                  !
!! Initialisation proportionnelle :                                 ! 
!! -------------------------------                                  !
!! L'incrément de déformation de chaque phase est calculé en        !
!! conservant la proportion de la déformation de cette phase vis à  !
!! vis de la déformation moyenne, par rapport à l'incrément de      !
!! déformation moyen imposé.                                        !
!! Si Eps1 = Alpha*EpsMoyen, alors DEps1 = Alpha*DEpsMoyen          !
!!                                                                  !
!! Initialisation linéaire :                                        ! 
!! -------------------------------                                  !
!! L'incrément de déformation de chaque phase est calculé en        !
!! conservant par interpolation linéaire basée sur le pas de temps  !
!! précédent :                                                      !
!! DEpsi = (dt/dt_old)*(Epsi,t - Epsi,t-dt)                         !
!!                                                                  !
!! Entrées/sorties : voir umatLaminate                              !
!!===================================================================
subroutine InitDeltaEpsLaminate(DdefLoc,Varint,Init,DdefIn,FV,nVarintloc,dt,iRotEps)

  implicit none

  ! Arguments
  real(mytype),dimension(:,:),intent(inout)        :: DdefLoc
  real(mytype),dimension(:),intent(in),target      :: Varint
  real(mytype),dimension(6),intent(in)             :: DDefIn
  real(mytype),dimension(:),intent(in)             :: Fv
  character(len=*),intent(in)                      :: Init
  integer(kind=INT64),dimension(:),intent(in)      :: nVarintloc !INT64 : compatibilite umat
  real(mytype),dimension(2),intent(in)             :: dt
  real(mytype),dimension(6,6)                      :: iRotEps

  ! Variables locales
  integer(kind=INT32)                              :: indVarint,nPhases
  integer(kind=INT32)                              :: i,j,k
  real(mytype),dimension(:),pointer                :: Epsloc,Epsloc_old,Epsloc_Temp
  real(mytype),dimension(3)                        :: EpsRef
  real(mytype)                                     :: Somme = 0_mytype
  character(len=16)                                :: Initloc

  nPhases = size(FV)

  Initloc = Init
  Epsloc     => NULL()
  Epsloc_old => NULL()
  EpsRef     =  DdefIn(IndAP)

  ! Initialisation homogène de l'incrément de déformation (par défaut)
  DdefLoc(:,:) = spread(DdefIn,2,nPhases)

  do i=1,nPhases-1

     indVarint  =  nPhases+(9*(i-1)) + int(sum(nVarintloc(1:i-1)),INT32)
     Epsloc     => Varint(indVarint+1:indVarint+3)
     Epsloc_old => Varint(indVarint+7:indVarint+9)

     if (Initloc == 'proportional' .or. Initloc == 'linear') then

        select case(Initloc)

        case('proportional')
           ! Initloclocialisation Proportionnelle de l'incrément de déformation 
           !   voir annexe stage A.Marano NT - DEN/DANS/DMN/SRMA/LC2M/NT/2016-3641/A
           do j=1,3
              if (Epsloc(j) > 1d-6) then
                 ! initialisation réalisée pour chaque composante ou la déformation n'est pas nulle
                 ! si elle l'est, on garde la valeur de l'incrément de déformation moyen
                 Somme = FV(i)
                 do k = i+1,nPhases
                    indVarint   =  nPhases+(9*(k-1)) + int(sum(nVarintloc(1:k-1)),INT32)
                    Epsloc_temp => Varint(indVarint+1:indVarint+3)
                    Somme = Somme + FV(k)*(Epsloc_temp(j)/Epsloc(j))              
                 end do
                 DdefLoc(indAP(j),i) = EpsRef(indAP(j)) / Somme
                 EpsRef(indAP(j))    = EpsRef(indAP(j)) - FV(i)*DdefLoc(indAP(j),i)
              end if
           end do

        case('linear')
           if ((dt(2) == 0)) then
              Initloc = ' '
           else if(all(Epsloc == Epsloc_old) ) then                   
              ! on garde l'initialisation par la déformation moyenne pour éviter un incrément nul de def
           else 
              ! Initloclocialisation par interpolation linéaire sur chaque phase de l'incrément de déformation
              DdefLoc(indAP,i) = (dt(1)/dt(2))*(Epsloc - Epsloc_old)  
           end if
        end select

     end if

     ! Calcul incrémental de la déformation de la phase N pour respecter la déformation moyenne imposée
     if (i==1) DdefLoc(indAP,nPhases) = DdefIn(indAP) / FV(nPhases) 
     DdefLoc(indAP,nPhases) = DdefLoc(indAP,nPhases) - (FV(i)/FV(nPhases))*DdefLoc(indAP,i)

     ! la déformation initiale est stockée à la place de la déformation du pas de temps
     ! précédent pour le prochain appel à umatLaminate
     Epsloc_old = Epsloc
     Ddefloc(:,i) = matmul(iRotEps,Ddefloc(:,i)) ! retour repère global
     if (i==nPhases-1)  then
        Ddefloc(:,i+1) = matmul(iRotEps,Ddefloc(:,i+1))
        indVarint  =  nPhases+(9*i) + int(sum(nVarintloc(1:i)),INT32)
        Epsloc     => Varint(indVarint+1:indVarint+3)
        Epsloc_old => Varint(indVarint+7:indVarint+9)
        Epsloc_old = Epsloc
     end if

  end do 
  ! supprime les erreurs dues aux rotations 
  where (abs(Ddefloc) < 1e-8) Ddefloc = 0.
end subroutine InitDeltaEpsLaminate

!!===================================================================
!!                      INIT_DELTA_EPS_REUSS                        !
!!                      --------------------                        !
!!                                                                  !
!! Initialisation de l'incrément de déformation pour l'algorithme   !
!! de Newton-Rapshon de umatReuss. Trois types d'initialisations    !
!!                                                                  !
!! Initialisation par défaut :                                      !
!! --------------------------                                       !
!! Homogène dans chaque phase égale à l'DEpsMoyen imposé            !
!!                                                                  !
!!                                                                  !
!! Initialisation proportionnelle :                                 ! 
!! -------------------------------                                  !
!! L'incrément de déformation de chaque phase est calculé en        !
!! conservant la proportion de la déformation de cette phase vis à  !
!! vis de la déformation moyenne, par rapport à l'incrément de      !
!! déformation moyen imposé.                                        !
!! Si Eps1 = Alpha*EpsMoyen, alors DEps1 = Alpha*DEpsMoyen          !
!!                                                                  !
!! Initialisation linéaire :                                        ! 
!! -------------------------------                                  !
!! L'incrément de déformation de chaque phase est calculé en        !
!! conservant par interpolation linéaire basée sur le pas de temps  !
!! précédent :                                                      !
!! DEpsi = (dt/dt_old)*(Epsi,t - Epsi,t-dt)                         !
!!                                                                  !
!! Entrées/sorties : voir umatLaminate                              !
!!===================================================================
subroutine InitDeltaEpsReuss(DdefLoc,Varint,Init,DdefIn,FV,nVarintloc,dt)

  implicit none

  ! Arguments
  real(mytype),dimension(:,:),intent(inout)        :: DdefLoc
  real(mytype),dimension(:),intent(in),target      :: Varint
  real(mytype),dimension(6),intent(in)             :: DDefIn
  real(mytype),dimension(:),intent(in)             :: Fv
  character(len=16),intent(in)                     :: Init
  integer(kind=INT64),dimension(:),intent(in)      :: nVarintloc !INT64 : compatibilite umat
  real(mytype),dimension(2),intent(in)             :: dt

  ! Variables locales
  integer(kind=INT32)                              :: indVarint,nPhases
  integer(kind=INT32)                              :: i,j,k
  real(mytype),dimension(:),pointer                :: Epsloc,Epsloc_old,Epsloc_Temp, Epsloc_old_temp
  real(mytype),dimension(6)                        :: EpsRef
  real(mytype)                                     :: Somme = 0_mytype
  character(len=16)                                :: Initloc
  logical                                          :: testPB

  nPhases = size(FV)
  Initloc = Init
  Epsloc     => NULL()
  Epsloc_old => NULL()
  EpsRef     =  DDefIn

  testPB=.false.

  ! Initialisation homogène de l'incrément de déformation (ce qui est fait par défaut)
  DdefLoc(:,:) = spread(DdefIn,2,nPhases)

  DOLOOP_phases:&
  do i=1,nPhases

     indVarint  =  nPhases+(12*(i-1)) + int(sum(nVarintloc(1:i-1)),INT32)
     Epsloc     => Varint(indVarint+1:indVarint+6)
     Epsloc_old => Varint(indVarint+7:indVarint+12)

     select case(trim(Initloc))

        !Assure la proportionalite entre les composantes des increments de deformation
        !Avantages :
        !      increment de deformation moyenne assuree : Bien pour changement de trajets de chargements
        !      ne necessite pas dt_old
        case('proportional')
           if (minval(abs(Epsloc - Epsloc_old)) < 1e-14) then !cas d'increments de deformation trop faibles 
                    testPB = .true.
                    exit DOLOOP_phases
           end if
           do j=1,6
                 Somme = FV(i)
                 do k = i+1,nPhases
                    indVarint   =  nPhases+(12*(k-1)) + int(sum(nVarintloc(1:k-1)),INT32)
                    Epsloc_temp => Varint(indVarint+1:indVarint+6)
                    Epsloc_old_temp => Varint(indVarint+7:indVarint+12)
                    Somme = Somme + (FV(k) * (Epsloc_temp(j)-Epsloc_old_temp(j)) / (Epsloc(j)-Epsloc_old(j)))
                 end do
                 ! somme = 0 sur la phase 1 si on a un increment moyen impose nul
                 !       dans ce cas, on conserve la valeur par defaut
                 if (abs(Somme) < 1e-14) then
                    testPB = .true.
                    exit DOLOOP_phases
                 end if 
                 DdefLoc(j,i) = EpsRef(j) / Somme
                 EpsRef(j)    = EpsRef(j) - FV(i)*DdefLoc(j,i)
           end do

        case('linear')
           !TODO (amelioration possible) : a voir si le nouveau cas "proportional" pose des pbs
           !ne faire du lineaire que s'il n'y a pas de changement trop marque de trajet de chargement
           if (dt(2) /= 0) then
              DdefLoc(:,i) = (dt(1)/dt(2))*(Epsloc - Epsloc_old)
           end if
     end select


     ! la déformation initiale est stockée à la place de la déformation du pas de temps
     ! précédent pour le prochain appel à umatReuss
     Epsloc_old = Epsloc
  end do DOLOOP_phases 

  ! Au moindre pb (RARE normalement), on reprend l'initilisation homogene
  if (testPB) DdefLoc(:,:) = spread(DdefIn,2,nPhases)
 

end subroutine InitDeltaEpsReuss



!!===================================================================
!!                      CALC_RESIDU_COMPOSITE
!!                      ---------------------
!!  Calcul la fonction Résidu utilisé pour caractériser la solution 
!!  du modèle d'homogeneisation 
!!  - Appel au comportement de chaque phase composant le voxel composite
!!
!!  - Gestion des erreurs spécifique
!!             - distinction entre non convergence d'une loi locale 
!!               et détection de Nan (pour subidivision du pas de temps dans le 1er cas)
!!
!!  - Calcul de la matrice tangente associée à l'incrément de déformation
!!    en entrée (optionnel)
!!
!!  Entrées/sorties : voir UmatLaminate
!!
!! 30/03/2022 : Ajout arguments optionnels pour appel UMAT "extension"
!!             
!!===================================================================


subroutine  CalcResiduComposite(numM,umatptr,SigMoy,Def0,DdefLoc,ncoeff,nvarint,Varint,Varint0,&
                               Coeff,times,dt,Fv,Residu,KINC,&
                               erreur,RotEps,iRotEps,RotSig,iRotSig,DF,umat_only,DDef)
  implicit none
  ! Arguments entrée-sortie 
  integer, intent(in)                                    :: numM ! numero LOCAL du materiau composite
  real(mytype),dimension(ntens),intent(inout)            :: SigMoy
  real(mytype),dimension(ntens),intent(in)               :: Def0
  real(mytype),dimension(:,:),intent(in)                 :: DdefLoc
  real(mytype),dimension(:,:),intent(inout)              :: DF
  integer(kind=8),intent(in)                             :: ncoeff,nvarint
  real(mytype),dimension(nvarint),intent(inout),target   :: Varint
  real(mytype),dimension(nvarint),intent(in)             :: Varint0
  real(mytype),dimension(ncoeff),intent(in),target       :: Coeff
  type(c_ptr),dimension(int(Coeff(1))), intent(in)       :: umatptr   !----Pointeur de fonction
  real(mytype),dimension(2),intent(in)                   :: times
  real(mytype),dimension(2),intent(in)                   :: dt
  real(mytype),dimension(:),intent(in)                   :: Fv
  real(mytype),dimension(:),intent(inout)                :: Residu
  real(mytype),dimension(6,6),intent(in)                 :: RotEps,iRotEps,RotSig,iRotSig
  integer(kind=8),intent(inout)                          :: KINC
  logical,intent(out)                                    :: erreur
  
  logical,intent(in),optional                            :: umat_only  ! ARGUMENTS OPTIONNELS POUR APPEL UMAT "EXTENSION"
  real(mytype),dimension(ntens),intent(in),optional      :: DDef

  !----Variables Sortie 
  real(mytype),dimension(ntens)                          :: Sig
  !----Variables erreur
  integer                                                :: alloc_stat
  character(len=200)                                     :: err_msg
  !----Variables de boucle
  integer(kind=INT32)                                    :: i,h
  integer(kind=INT32)                                    :: indVarint
  ! variables locales
  integer(kind=INT32)                                    :: nPhases
  integer(kind=INT64),dimension(:),allocatable           :: nCoeffloc,nVarintloc
                                                         ! INT64 pour compatibilite interface umat
  integer(kind=INT64)                                    :: nVarintloc_umatcall  
                                                         ! if (nVarintloc(i)=0) nVarintloc_umatcall=1
                                                         ! else nVarintloc_umat=nVarintloc(i)  
                                                         ! => for compatibility (MFRONT/UMAT/CAST3M)
  real(mytype),dimension(6)                              :: Sigloc, Epsloc
  real(mytype),dimension(:),pointer                      :: Varintloc,Coeffloc 
  real(mytype), dimension(3,3)                           :: DFGRD0_loc, DFGRD1_loc
  real(mytype),dimension(3,3)                            :: Cloc,Cloc1 
  real(mytype),dimension(6)                              :: DEpsloc

  
  logical                                                :: umat_only0

  
  !initialisation
  indVarint = 0 ! avoid gcc-warning
  umat_only0 = .false.
  if (present(umat_only)) umat_only0 = umat_only
  
  ! Nombre de Phases dans le matériau
  nPhases = int(Coeff(1),INT32)

  ! Allocation des tableaux NCoeff, NVarInt
  allocate(nVarintloc(nPhases),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort(&
               "Espace memoire disponible insuffisant (CalcResiduComposite)",2)
  allocate(nCoeffloc(nPhases),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort(&
               "Espace memoire disponible insuffisant (CalcResiduComposite)",2)

  ! Mise à l'état NULL des pointeurs
  Varintloc  => NULL()
  Coeffloc   => NULL()

  ! Nombre de varint de chaque phase
  nVarintloc = int(Varint(1:nPhases),kind=INT64)

  ! Nombre de coeff pour chaque phase
  nCoeffloc  = int(Coeff(2:nPhases+1),kind=INT64)

  Sig = Sigmoy  ! sauvegarde l'entree sigMoy
  Sigmoy = 0.
  ddsdde = 0.
  

  ! Boucle de calcul du comportement, évalutation du résidu et de la contrainte/varint pour chaque phase
  do i=1,nPhases
  
     !--------------------------------------- RECUPERATION DES TENSEURS LOCAUX ET DES COEFFICIENTS
     
     if (.not. umat_only0) then
        call get_tenseurs_locaux(numM,i,nPhases,indVarint,nCoeffloc,nVarintloc,Coeff,Varint,Varint0,&
                            Varintloc,Coeffloc,RotEps,iRotEps,RotSig,iRotSig,Def0,Sig,Sigloc,Epsloc)
     else
        call get_tenseurs_locaux(numM,i,nPhases,indVarint,nCoeffloc,nVarintloc,Coeff,Varint,Varint0,&
                            Varintloc,Coeffloc,RotEps,iRotEps,RotSig,iRotSig,Def0,Sig,Sigloc,Epsloc,&
                            DDef,Depsloc)
     end if


     !! Si on veut evaluer une loi GDEF avec une resolution en HPP on passe
     !! le tenseur gradient de la transformation (symmetrique) ---------- A VOIR !!!!!
     ! gradient de def initial
     DFGRD0_loc(1,1) = dble(1.+Epsloc(1))
     DFGRD0_loc(2,1) = dble(Epsloc(4))/2._mytype
     DFGRD0_loc(3,1) = dble(Epsloc(5))/2._mytype
     DFGRD0_loc(1,2) = dble(Epsloc(4))/2._mytype
     DFGRD0_loc(2,2) = dble(1.+Epsloc(2))
     DFGRD0_loc(3,2) = dble(Epsloc(6))/2._mytype
     DFGRD0_loc(1,3) = dble(Epsloc(5))/2._mytype
     DFGRD0_loc(2,3) = dble(Epsloc(6))/2._mytype
     DFGRD0_loc(3,3) = dble(1.+Epsloc(3))
     ! gradient de def t+dt
     if (.not. umat_only0) then
     DFGRD1_loc(1,1) = dble(1.+ Epsloc(1) + DdefLoc(1,i))
     DFGRD1_loc(2,1) = dble(Epsloc(4) + DdefLoc(4,i))/2._mytype
     DFGRD1_loc(3,1) = dble(Epsloc(5) + DdefLoc(5,i))/2._mytype
     DFGRD1_loc(1,2) = dble(Epsloc(4) + DdefLoc(4,i))/2._mytype
     DFGRD1_loc(2,2) = dble(1.+Epsloc(2) + DdefLoc(2,i))
     DFGRD1_loc(3,2) = dble(Epsloc(6) + DdefLoc(6,i))/2._mytype
     DFGRD1_loc(1,3) = dble(Epsloc(5) + DdefLoc(5,i))/2._mytype
     DFGRD1_loc(2,3) = dble(Epsloc(6) + DdefLoc(6,i))/2._mytype
     DFGRD1_loc(3,3) = dble(1.+Epsloc(3) + DdefLoc(3,i)) 
     else
     DFGRD1_loc(1,1) = dble(1.+ Epsloc(1) + Depsloc(1))
     DFGRD1_loc(2,1) = dble(Epsloc(4) + Depsloc(4))/2._mytype
     DFGRD1_loc(3,1) = dble(Epsloc(5) + Depsloc(5))/2._mytype
     DFGRD1_loc(1,2) = dble(Epsloc(4) + Depsloc(4))/2._mytype
     DFGRD1_loc(2,2) = dble(1.+Epsloc(2) + Depsloc(2))
     DFGRD1_loc(3,2) = dble(Epsloc(6) + Depsloc(6))/2._mytype
     DFGRD1_loc(1,3) = dble(Epsloc(5) + Depsloc(5))/2._mytype
     DFGRD1_loc(2,3) = dble(Epsloc(6) + Depsloc(6))/2._mytype
     DFGRD1_loc(3,3) = dble(1.+Epsloc(3) + Depsloc(3))      
     end if


     ! passage à Mfront de l'opérateur tangent à calculer (à implémenter dans
     ! le fichier Mfront : ddsdde(1,1) correspond à la variable mfront "smt"
     ! --> écrire un test dans le fichier mfront pour calculer le bon opérateur
     ! 1 - tenseur d'élasticité
     ! 4 - matrice tangente cohérente
     
     if( (trim(Mattotp(numM)%lawname) == 'laminate') .and. &
                (trim(MattotP(numM)%Jacobian_type_composite) == 'mfront_elastic') ) ddsdde(1,1) = 1
     if( (trim(Mattotp(numM)%lawname) == 'laminate') .and. &
                (trim(MattotP(numM)%Jacobian_type_composite) == 'mfront_secant') )  ddsdde(1,1) = 2
     if( (trim(Mattotp(numM)%lawname) == 'laminate') .and. &
                (trim(MattotP(numM)%Jacobian_type_composite) == 'mfront_tangent') ) ddsdde(1,1) = 4
     

     
     ! appel au comportement phase i (les variables sont dans le repère global !!!)
     !!-------------------------------------------------------------------------

     nVarintloc_umatcall = nVarintloc(i)     
     if (nVarintloc_umatcall == 0) nVarintloc_umatcall = 1 !pour compatibilite MFRONT (et compatibilite umat/CAST3M) : nVarint>=1  
     if (.not. umat_only0) then
     call umat_call(umatptr(i),Sigloc,Varintloc,ddsdde,sse,spd,scd, &
          RPL, DDSDDT, DRPLDE, DRPLDT,& 
          Epsloc,Ddefloc(:,i),times,dt(1),&
          temp,dtemp,PREDEF,DPRED,&
          CMNAME,NDI,nshr,ntens,nVarintloc_umatcall,&    
          Coeffloc,nCoeffloc(i),coords,&        
          DROT,pnewdt,celent,DFGRD0_loc, DFGRD1_loc,&
          NOEL, NPT, LAYER, KSPT, KSTEP,KINC)
     else
     call umat_call(umatptr(i),Sigloc,Varintloc,ddsdde,sse,spd,scd, &
          RPL, DDSDDT, DRPLDE, DRPLDT,& 
          Epsloc,Depsloc,times,dt(1),&
          temp,dtemp,PREDEF,DPRED,&
          CMNAME,NDI,nshr,ntens,nVarintloc_umatcall,&    
          Coeffloc,nCoeffloc(i),coords,&        
          DROT,pnewdt,celent,DFGRD0_loc, DFGRD1_loc,&
          NOEL, NPT, LAYER, KSPT, KSTEP,KINC)
     end if
     
     ! Tests Erreurs
     if (KINC<0 .AND. KINC/=-10) then
        erreur = .true.
        KINC = -20
        call amitex_abort("(CalcResiduComposite)",2)
        exit
     end if
     
     do h=1,6
        if (ieee_is_nan(Sigloc(h))) then
           write(err_msg,'(A,I2,A)') &
                "Contrainte à la sortie de umat NaN pour la phase :",i,"  (CalcResiduComposite)"
           call amitex_abort(err_msg,0)
           erreur = .true.
           KINC = -20
        end if
        exit
     end do

     do h=1,int(nVarintloc(i),INT32)
        if (ieee_is_nan(Varintloc(h))) then
           write(err_msg,'(A,I2,A,I4,A)') &
              "Variable(s) interne(s) à la sortie de umat NaN pour la phase :",i," (CalcResiduComposite)"
           call amitex_abort(err_msg,0)
           erreur = .true.
           KINC = -20
           exit
        end if
     end do

   
     if (.not. umat_only0) then ! partie shuntee pour un appel umat "extension"
     ! Mise à jour de la matrice tangente du système (assemblage)
     !!----------------------------------------------

     if(  (trim(Mattotp(numM)%lawname) == 'laminate') .AND. &
         ((trim(Mattotp(numM)%Jacobian_type_composite) == 'mfront_tangent') .or. &
          (trim(Mattotp(numM)%Jacobian_type_composite) == 'mfront_elastic'))       ) then

        if (i == 1) then
           Cloc1 = matmul(RotSig(IndAP,:),  matmul(ddsdde,iRotEps(:,IndAP)))
        else
           Cloc = matmul(RotSig(IndAP,:),  matmul(ddsdde,iRotEps(:,IndAP)))

           DF(:,(i-2)*3+1:3*(i-1)) = (FV(i)/FV(1)) * &
                      reshape( spread( Cloc1,2,nPhases-1 ), (/ int(3*(nPhases-1)) , 3 /))
           DF((i-2)*3+1:3*(i-1),(i-2)*3+1:3*(i-1)) = DF((i-2)*3+1:3*(i-1),(i-2)*3+1:3*(i-1)) + Cloc
        end if

        ddsdde = 0.

     end if


     ! Calcul du residu (nul pour voigt, ecart des contraintes pour reuss et ecart des contrainte planes pour laminate)
     !-----------------------------------------
     
     call assemblage_residu(numM,i,nPhases,indVarint,Varint,RotSig,Sigloc,Residu)
     
     
     ! Evaluation de la contrainte moyenne
     !------------------------------------
     Sigmoy = Sigmoy + FV(i)*Sigloc
     else   ! fin partie shuntee pour un appel umat "extension"
     Sigmoy = Sig  ! sigmoy ne doit pas etre modifie lors d'un appel a umat-extension
     end if 
     
  end do
  
  
  
  ! Mise à l'état NULL des pointeurs
  Varintloc  => NULL()
  Coeffloc   => NULL()

  ! Desallocation des tableaux NCoeff, NVarInt
  deallocate(nVarintloc)
  deallocate(nCoeffloc)
  
  
end subroutine  CalcResiduComposite



!==============================================================================
!
!                      SUBROUTINE GET_TENSEURS_LOCAUX
!
!> Initialise les coefficients et les variables internes de chaque phase, ainsi que les valeurs
!! des tenseurs de contraintes et de deformations locaux avant le calcul du residu.
!!
!! 30/03/2022 : Ajout arguments optionnels pour utilisation UMAT "extension"
!!              via appel a umatcomposite_ext appelant calcresiducomposite
!==============================================================================

subroutine get_tenseurs_locaux(numM,i,nPhases,indVarint,nCoeffloc,nVarintloc,Coeff,Varint,Varint0,&
                            Varintloc,Coeffloc,RotEps,iRotEps,RotSig,iRotSig,Def0,Sig,Sigloc,Epsloc,DDef,Depsloc)

implicit none

  integer,intent(in)                                :: numM ! numero local du materiau composite
  integer,intent(in)                                :: i    ! numero de phase
  integer,intent(in)                                :: nPhases ! nombre de phases
  integer(kind=INT32),intent(out)                   :: indVarint
  integer(kind=INT64),dimension(:),intent(in)       :: nCoeffloc,nVarintloc ! nombre de coeff/varint par phase
  real(mytype),dimension(:),intent(in),target       :: Coeff
  real(mytype),dimension(:),intent(inout),target    :: Varint
  real(mytype),dimension(:),intent(in)              :: Varint0
  real(mytype),dimension(:),intent(inout),pointer   :: Varintloc,Coeffloc 
  real(mytype),dimension(:),intent(in)              :: Def0
  real(mytype),dimension(:),intent(in)              :: Sig
  real(mytype),dimension(6,6),intent(in)            :: RotEps,iRotEps,RotSig,iRotSig
  real(mytype),dimension(6),intent(out)             :: Sigloc, Epsloc
  integer(kind=INT32)                               :: indCoeff

  real(mytype),dimension(:),intent(in),optional     :: DDef  
  real(mytype),dimension(6),intent(out),optional    :: DEpsloc

    

  ! Cas Laminate
  if(MattotP(numM)%lawname=="laminate") then
       
     indCoeff   =  2*nPhases + 6 + 2*(i) + int(sum(NCoeffloc(1:i-1)),INT32)
     indVarint  =  nPhases+9*(i-1) + int(sum(nVarintloc(1:i-1)),INT32)

     ! Reconstruction de la déformation initiale de la phase i dans le repère de l'interface
     ! puis passage au repère global
     Epsloc        = matmul(RotEps,def0)                      !> déformation plane imposée par la valeur moyenne
     Epsloc(IndAp) = Varint0(indVarint+1:indVarint+3)   !> déformation normale récupérée depuis VarInt
     Epsloc        = matmul(iRotEps,Epsloc)  ! retour repère global

     ! reconstruction de la contrainte initiale de la phase i
!TODO - ERRROR a l'execution avec openapi2022 + DEBUG = -g -check all -fpe0 -warn -traceback -debug extended
!pas compris
!forrtl: severe (408): fort: (3): Subscript #1 of the array SIGLOC has value 0 which is less than the lower bound of 1
!if (nrank==0) print *, indAP, shape(Sigloc)
!print *, "nrank, Sigloc(IndPl) ",nrank, Sigloc(IndPl)
!print *, "nrank, sigloc(indAP) ",nrank, sigloc(indAP)  
     Sigloc(IndAP(:)) = matmul(RotSig(indAP,:),Sig)        !> contrainte normale imposée par la valeur moyenne
     Sigloc(IndPl) = Varint0(indVarint+4:indVarint+6)   !> contrainte plane récupérée depuis VarInt0 (sigma init)
     Sigloc        = matmul(iRotSig,Sigloc)  ! retour repère global
        
     ! récupération des Coeffs/Varint
     Varint( indVarint+10:indVarint+9+nVarintloc(i) ) =  &
                 Varint0( indVarint+10:indVarint+9+nVarintloc(i) )  ! reprise de la valeur initiale des variables internes
     Varintloc => Varint( indVarint+10:indVarint+9+nVarintloc(i) )
     Coeffloc  => Coeff( indCoeff+1:indCoeff+NCoeffloc(i) )
     
     if (present(DDef) .and. present(Depsloc)) then
        Epsloc        = matmul(RotEps,def0)                      !> déformation plane imposée par la valeur moyenne
        Epsloc(IndAp) = Varint0(indVarint+7:indVarint+9)   !> déformation normale récupérée depuis VarInt
        Epsloc        = matmul(iRotEps,Epsloc)  ! retour repère global
        Depsloc        = matmul(RotEps,Ddef)                !> déformation plane imposée par la valeur moyenne
        Depsloc(IndAp) = Varint0(indVarint+1:indVarint+3)-Varint0(indVarint+7:indVarint+9)   !> déformation normale récupérée depuis VarInt
        Depsloc        = matmul(iRotEps,Depsloc)             ! retour repère global
     end if
     
  ! Cas Reuss
  elseif (MattotP(numM)%lawname=="reuss") then

     indCoeff   =  2*nPhases + 2*(i) + int(sum(NCoeffloc(1:i-1)),INT32)
     indVarint  =  nPhases+12*(i-1) + int(sum(nVarintloc(1:i-1)),INT32)

     ! déformation initiale de la phase i récupérée depuis VarInt
     Epsloc = Varint0(indVarint+1:indVarint+6)
        
     ! contrainte initiale de la phase i egale a la contrainte moyenne
     Sigloc = Sig
        
     ! récupération des Coeffs/Varint
     Varint( indVarint+13:indVarint+12+nVarintloc(i) ) =  &
                  Varint0( indVarint+13:indVarint+12+nVarintloc(i) )  
           ! reprise de la valeur initiale des variables internes
     Varintloc => Varint( indVarint+13:indVarint+12+nVarintloc(i) )
     Coeffloc  => Coeff( indCoeff+1:indCoeff+NCoeffloc(i) )
     if (present(DDef) .and. present(Depsloc)) then
        Epsloc  = Varint0(indVarint+7:indVarint+12)
        Depsloc = Varint0(indVarint+1:indVarint+6)-Varint0(indVarint+7:indVarint+12)
     end if
        
  ! Cas Voigt
  elseif (MattotP(numM)%lawname=="voigt") then
     
     indCoeff   =  2*nPhases + int(sum(NCoeffloc(1:i-1)),INT32)
     indVarint  =  nPhases+6*(i-1) + int(sum(nVarintloc(1:i-1)),INT32)

     ! déformation initiale de la phase i egale a la deformation moyenne
     Epsloc = def0
       
     ! contrainte initiale de la phase i recuperee depuis les variables internes
     Sigloc = Varint0(indVarint+1:indVarint+6)
     
     ! récupération des Coeffs/Varint
     Varint( indVarint+7:indVarint+6+nVarintloc(i) ) =  &
                  Varint0( indVarint+7:indVarint+6+nVarintloc(i) )  
           ! reprise de la valeur initiale des variables internes
     Varintloc => Varint( indVarint+7:indVarint+6+nVarintloc(i) )
     Coeffloc  => Coeff( indCoeff+1:indCoeff+NCoeffloc(i) )
     if (present(DDef) .and. present(Depsloc)) then
        Depsloc = DDef
     end if
        
  end if


end subroutine get_tenseurs_locaux



!==============================================================================
!
!                      SUBROUTINE ASSEMBLAGE_RESIDU
!
!> Assemble le residu a partir des calculs des contraintes des differentes phases et
!! range la valeur de la contrainte dans les variables internes du materiau composite
!!
!!
!==============================================================================

subroutine assemblage_residu(numM,i,nPhases,indVarint,Varint,RotSig,Sigloc,Residu)

implicit none

  integer,intent(in)                                :: numM ! numero local du materiau composite
  integer,intent(in)                                :: i    ! numero de phase
  integer,intent(in)                                :: nPhases ! nombre de phases
  integer(kind=INT32),intent(in)                    :: indVarint
  real(mytype),dimension(:),intent(inout),target    :: Varint
  real(mytype),dimension(6,6),intent(in)            :: RotSig
  real(mytype),dimension(6),intent(in)              :: Sigloc
  real(mytype),dimension(:),intent(inout)           :: Residu

  real(mytype),dimension(3,int(nPhases-1,INT32))    :: ltmp
  real(mytype),dimension(6,int(nPhases-1,INT32))    :: rtmp


  ! Cas Laminate
  if(MattotP(numM)%lawname=="laminate") then
     ! Mise à jour de la contrainte plane (valeur conservée en cas de convergence, pas appelée dans le calcul)
     !------------------------------------
     Varint(indVarint+4:indVarint+6) = matmul(RotSig(IndPl,:),Sigloc)
     if (i==1) then
        ! - Sigma_1 sur tous les termes du résidu
        ltmp = spread(matmul(RotSig(IndAP,:),Sigloc),2,nPhases-1)
        Residu(:) = Residu(:) - reshape(ltmp, (/3*(nPhases-1)/) )
     else
        ! + Sigma_i sur le terme correspondant du résidu
        Residu(3*(i-2)+1:3*(i-1)) = Residu(3*(i-2)+1:3*(i-1)) + matmul(RotSig(IndAP,:),Sigloc)
     end if

  ! Cas Reuss
  elseif (MattotP(numM)%lawname=="reuss") then
     if (i==1) then
        ! - Sigma_1 sur tous les termes du résidu
        rtmp = spread(Sigloc,2,nPhases-1)
        Residu(:) = Residu(:) - reshape(rtmp, (/int(6*(nPhases-1)) /) )
     else
        ! + Sigma_i sur le terme correspondant du résidu
        Residu(6*(i-2)+1:6*(i-1)) = Residu(6*(i-2)+1:6*(i-1)) + Sigloc
     end if
     
  ! Cas Voigt
  elseif (MattotP(numM)%lawname=="voigt") then
     Residu = 0
     ! Mise a jour de la contrainte
     Varint(indVarInt+1:indVarInt+6)=SigLoc
     
  end if



end subroutine assemblage_residu



!!===================================================================
!!                 MATRICE_TANGENTE_ELASTIQUE
!!                 --------------------------
!!
!!  Calcule la matrice elastique "equivalente" du
!!  modèle multi-couche pour intégration avec la méthode de 
!!  Newton-Raphson modifiée
!!
!!  Entrées : voir umatLaminate
!!  \param[out] DF0 : matrice Newton-modifié
!!===================================================================
subroutine  MatriceTangenteElastiqueLaminate(Coeff,DF0,FV)

  implicit none
  !----Variables entrees-sorties
  real(mytype),dimension(:,:),intent(inout)               :: DF0 
  real(mytype),dimension(:),intent(in),target             :: Coeff
  real(mytype),dimension(:),intent(in)                    :: Fv
  !----Variables de boucle
  integer                                                 :: i
  integer                                                 :: indCoeff
  !----Variables locales
  integer(kind=INT32)                                     :: nPhases  
  real(mytype),dimension(3,3)                             :: Cloc,Cloc1
  real(mytype),dimension(:),pointer                       :: LambdaMuEq
  integer(kind=INT32),dimension(:),allocatable            :: nCoeffloc
  integer                                                 :: alloc_stat

  ! Mise à l'état NULL du pointeur
  LambdaMuEq => NULL()
  ! Nombre de Phases dans le matériau
  nPhases = int(Coeff(1),INT32)

    ! Allocation des tableaux NCoeff, NVarInt
  allocate(nCoeffloc(nPhases),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (CalcResiduLaminate)",2)
  ! Nombre de coeff pour chaque phase
  nCoeffloc  = int(Coeff(2:nPhases+1),kind=INT32)

  Cloc1(:,:) = 0.
  Cloc(:,:)  = 0.

  ! Matrice tangente de la phase 1 
  indCoeff   = 2*nPhases + 7
  LambdaMuEq => Coeff(indCoeff:indCoeff+1)
  Cloc1(1,1) = 2*LambdaMuEq(2) + LambdaMuEq(1)
  Cloc1(2,2) = LambdaMuEq(2)
  Cloc1(3,3) = LambdaMuEq(2)

  ! Remplissage de la matrice tangente phase par phase colonne par colonne

  do i=2,nPhases
     ! Matrice tangente de la phase i
     indCoeff   = 2*nPhases + 7 + 2*(i-1) +sum(NCoeffloc(1:i-1))
     LambdaMuEq => Coeff(indCoeff:indCoeff+1)
     Cloc(1,1) = 2*LambdaMuEq(2) + LambdaMuEq(1)
     Cloc(2,2) = LambdaMuEq(2)
     Cloc(3,3) = LambdaMuEq(2)

     DF0(:,(i-2)*3+1:3*(i-1)) = (FV(i)/FV(1)) * reshape( spread( Cloc1,2,nPhases-1 ), (/ int(3*(nPhases-1)) , 3 /))

     DF0((i-2)*3+1:3*(i-1),(i-2)*3+1:3*(i-1)) = DF0((i-2)*3+1:3*(i-1),(i-2)*3+1:3*(i-1)) + Cloc
  end do

  LambdaMuEq => NULL()
  deallocate(nCoeffloc)

end subroutine MatriceTangenteElastiqueLaminate 

!!===================================================================
!!                 MATRICE_TANGENTE_ELASTIQUE_REUSS
!!                 ---------------------------------
!!
!!  Calcule la matrice elastique "equivalente" du
!!  modèle de Reuss pour intégration avec la méthode de 
!!  Newton-Raphson modifiée
!!
!!  Entrées : voir umatReuss
!!  \param[out] DF0 : matrice 'tangente' pour Newton-modifié
!!===================================================================

!---------------------------------------------------------------------------------
!TODO : Fusionner les fonctions MatriceTangentElastique des que possible
!---------------------------------------------------------------------------------
subroutine  MatriceTangenteElastiqueReuss(Coeff,DF0,FV)

  implicit none
  !----Variables entrees-sorties
  real(mytype),dimension(:,:),intent(inout)               :: DF0 
  real(mytype),dimension(:),intent(in),target             :: Coeff
  real(mytype),dimension(:),intent(in)                    :: Fv
  !----Variables de boucle
  integer                                                 :: i
  integer                                                 :: indCoeff
  !----Variables locales
  integer(kind=INT32)                                     :: nPhases  
  real(mytype),dimension(6,6)                             :: Cloc,Cloc1
  real(mytype),dimension(:),pointer                       :: LambdaMuEq
  integer(kind=INT32),dimension(:),allocatable            :: nCoeffloc
  integer                                                 :: alloc_stat

  ! Mise à l'état NULL du pointeur
  LambdaMuEq => NULL()
  ! Nombre de Phases dans le matériau
  nPhases = int(Coeff(1),INT32)

  ! Allocation des tableaux NCoeff, NVarInt
  allocate(nCoeffloc(nPhases),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (CalcResiduReuss)",2)
  ! Nombre de coeff pour chaque phase
  nCoeffloc  = int(Coeff(2:nPhases+1),kind=INT32)

  Cloc1(:,:) = 0.
  Cloc(:,:)  = 0.

  ! Matrice tangente de la phase 1 
  indCoeff   = 2*nPhases + 1
  LambdaMuEq => Coeff(indCoeff:indCoeff+1)
  Cloc1(1,1) = 2*LambdaMuEq(2) + LambdaMuEq(1)
  Cloc1(2,2) = Cloc1(1,1)
  Cloc1(3,3) = Cloc1(1,1)
  Cloc1(4,4) = LambdaMuEq(2)
  Cloc1(5,5) = Cloc1(4,4)
  Cloc1(6,6) = Cloc1(4,4)
  Cloc1(1,2) = LambdaMuEq(1)
  Cloc1(2,1) = Cloc1(1,2)
  Cloc1(1,3) = Cloc1(1,2)
  Cloc1(3,1) = Cloc1(1,2)
  Cloc1(2,3) = Cloc1(1,2)
  Cloc1(3,2) = Cloc1(1,2)

  ! Remplissage de la matrice tangente phase par phase colonne par colonne

  do i=2,nPhases
     ! Matrice tangente de la phase i
     indCoeff   = 2*nPhases + 1 + 2*(i-1) +sum(NCoeffloc(1:i-1))
     LambdaMuEq => Coeff(indCoeff:indCoeff+1)
     Cloc(1,1) = 2*LambdaMuEq(2) + LambdaMuEq(1)
     Cloc(2,2) = Cloc(1,1)
     Cloc(3,3) = Cloc(1,1)
     Cloc(4,4) = LambdaMuEq(2)
     Cloc(5,5) = Cloc(4,4)
     Cloc(6,6) = Cloc(4,4)
     Cloc(1,2) = LambdaMuEq(1)
     Cloc(2,1) = Cloc(1,2)
     Cloc(1,3) = Cloc(1,2)
     Cloc(3,1) = Cloc(1,2)
     Cloc(2,3) = Cloc(1,2)
     Cloc(3,2) = Cloc(1,2)

     DF0(:,(i-2)*6+1:6*(i-1)) =  &
        (FV(i)/FV(1)) * reshape( spread( Cloc1,2,nPhases-1 ), (/ int(6*(nPhases-1)) , 6 /))

     DF0((i-2)*6+1:6*(i-1),(i-2)*6+1:6*(i-1)) = DF0((i-2)*6+1:6*(i-1),(i-2)*6+1:6*(i-1)) + Cloc

  end do

  LambdaMuEq => NULL()

end subroutine MatriceTangenteElastiqueReuss



!!===================================================================!
!!                 MATRICE_TANGENTE_NUMERIQUE
!!                 --------------------------
!!
!!  Calcule la matrice tangente du  modèle multi-couche pour intégration 
!!  avec la méthode de Newton-Raphson 
!!  Calcul effectué par différenciation numérique sur chaque composante
!!  du vecteur de variables (DdefLoc = incrément de déformation antiplane)
!!
!!  --> La fonction résidu est évalué pour chaque perturbation et les composantes
!!      de sa matrice jacobienne sont évaluées par la méthode des différences finies 
!!
!!  La valeur de la perturbation affectée à chaque composante est donnée
!!  par algo_param%Perturbation_laminate --> choix utilisateur (par défaut 1e-8)
!!
!!  Entrées : voir umatLaminate
!!  \param[out] DF : matrice tangente
!!===================================================================
subroutine  MatriceTangenteNumerique(numM,umatptr,Sig,Def0,DdefLoc,ncoeff,nvarint,Varint0,&
                               Coeff,times,dt,Fv,Residu,KINC,&
                               RotEps,iRotEps,RotSig,iRotSig,erreur,DF)

  implicit none
  ! Arguments entrée-sortie 
  integer, intent(in)                                    :: numM ! numero local du materiau composite
  real(mytype),dimension(ntens),intent(in)               :: Sig
  real(mytype),dimension(ntens),intent(in)               :: Def0
  real(mytype),dimension(:,:),intent(in)                 :: DdefLoc
  real(mytype),dimension(:,:),intent(inout)              :: DF
  integer(kind=8),intent(in)                             :: ncoeff,nvarint
  real(mytype),dimension(nvarint),intent(in)             :: Varint0
  real(mytype),dimension(ncoeff),intent(in)              :: Coeff
  type(c_ptr),dimension(int(Coeff(1))), intent(in)       :: umatptr   !----Pointeur de fonction
  real(mytype),dimension(2),intent(in)                   :: times
  real(mytype),dimension(2),intent(in)                   :: dt
  real(mytype),dimension(:),intent(in)                   :: Fv
  real(mytype),dimension(:),intent(in)                   :: Residu
  real(mytype),dimension(6,6),intent(in)                 :: RotEps,iRotEps,RotSig,iRotSig
  integer(kind=8),intent(inout)                          :: KINC
  logical,intent(out)                                    :: erreur
  !----Variables erreur
  integer                                                :: alloc_stat
  !----Variables de boucle
  integer(kind=INT32)                                    :: i,j
  ! variables locales
  integer(kind=INT32)                                    :: nPhases,indColonne
  real(mytype),dimension(:),allocatable                  :: Residu2,DRes  
  real(mytype),dimension(:,:),allocatable                :: Ddefloc_perturbe
  real(mytype),dimension(6,3)                            :: DEps_pert
  real(mytype),dimension(nvarint)                        :: VarI
  real(mytype),dimension(ntens)                          :: Sig2
  
  
  ! Nombre de Phases dans le matériau
  nPhases = int(Coeff(1),kind=INT32)

  ! Allocation du résidu
  allocate(Residu2(3*(nPhases-1)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (MatriceTangenteNumerique)",2)
  Residu2 = 0._mytype
  allocate(DRes(3*(nPhases-1)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (MatriceTangenteNumerique)",2)
  DRes = 0._mytype  ! gcc-warning accepte

  ! Allocation du tableau des incréments de déformation perturbé de chaque phase
  allocate(DdefLoc_perturbe(6,nPhases),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (umatLaminate)",2)
  DdefLoc_perturbe = 0._mytype

  DF = 0._mytype

  ! Calcul de l'incrément de déformation élémentaire pour différentiation numérique dans le repère global
  ! pour une perturbation dans le repère de l'interface
  !-------------------------------------------------------------------------------------------------------
  DEps_pert(:,1) = (/1.,0.,0.,0.,0.,0. /)
  DEps_pert(:,2) = (/0.,0.,0.,1.,0.,0. /) 
  DEps_pert(:,3) = (/0.,0.,0.,0.,1.,0. /) 
  DEps_pert = DEps_pert*algo_param%Perturbation_laminate 

  Deps_pert = matmul(iRotEps,DEps_pert)
  
  ! Petite variation de la déformation anti plane  pour chaque phase
  do j=2,nPhases

     ! Petite variation de chacune des trois composantes de la déformation anti plane
     do i = 1,3

        indColonne = i+(j-2)*3

        VarI = Varint0
        Sig2 = Sig

        Ddefloc_perturbe = DdefLoc
        Ddefloc_perturbe(:,j) = Ddefloc_perturbe(:,j) + Deps_pert(:,i)
        Ddefloc_perturbe(:,1) = Ddefloc_perturbe(:,1) - (FV(j)/FV(1))*Deps_pert(:,i)

        call CalcResiduComposite(numM,umatptr,Sig2,Def0,DdefLoc_perturbe,ncoeff,nvarint,VarI,Varint0,&
             Coeff,times,dt,Fv,Residu2,KINC,&
             erreur,RotEps,iRotEps,RotSig,iRotSig,DF)    

        ! test erreur
        if (erreur) then
           exit
        end if

        DRes = Residu2 - Residu
        DRes = DRes /algo_param%Perturbation_laminate
        
        DF(:,indColonne) = DRes

        DRes       = 0.
        Residu2    = 0.

     end do

     ! test erreur
     if (erreur) then
        exit
     end if

  end do

  deallocate(Residu2)
  deallocate(Dres)
  deallocate(Ddefloc_perturbe)

end subroutine MatriceTangenteNumerique

!==================================================================================================
!!                  EVAL_NB_MATERIAUX
!!           ------------------------------
!!
!! Subroutine appelee dans amitex_fftp apres avoir initialise les tableaux MattotP
!!
!!  > Sert a compter les differents nombres de materiaux :
!!
!! deja calcule  nmateriaux:                    (entier) nombre de materiaux dans la cellule
!! deja calcule  nmateriaux_composites:         (entier) nombre de materiaux composites dans la cellule
!! \param[out]   nz(1:nmateriaux)               (entier-kind=8) nombre de zone par materiau
!! \param[out]   nmateriaux_non_composites:     (entier) nombre de materiaux non composites dans la cellule
!! \param[out]   nmateriaux_loc:                (entier) nombre de materiaux dans le pinceau
!! \param[out]   nmateriaux_composites_loc:     (entier) nombre de materiaux composites dans le pinceau
!! \param[out]   nmateriaux_non_composites_loc: (entier) nombre de materiaux non composites dans le pinceau
!!
!--------------------------------------------------------------------------------------------------

subroutine eval_nb_materiaux

  implicit none


  integer             :: i,j                ! indice de boucle

  integer             :: check_nbmat, &     ! nmateriaux et nmateriaux_composites ont deja ete calcules,
                         check_nbmat_comp   ! on verifie qu'ils ont la bonne valeur

  integer             :: ierror             ! erreur MPI
  
  integer(kind=8)     :: nZone


  nmateriaux0%n_loc = 0
  nmateriaux0%n_composites_loc = 0
  nmateriaux0%n_non_composites_loc = 0
  check_nbmat = 0
  check_nbmat_comp = 0
  nmateriaux0%n_non_composites = 0

  ! Evaluation des nombres de materiaux non-composites et totaux sur la cellule
  do i = 1, size(MattotP0)
     nmateriaux0%n_loc = max(nmateriaux0%n_loc,MattotP0(i)%numM)
     if(MattotP0(i)%NPhase==1) nmateriaux0%n_non_composites_loc = &
                              max(nmateriaux0%n_non_composites_loc,MattotP0(i)%numM)
  end do

  call MPI_BARRIER(amitex_mpi_comm,ierror) 
  call MPI_Allreduce(nmateriaux0%n_non_composites_loc,nmateriaux0%n_non_composites,&
                     1,MPI_INTEGER,MPI_MAX,amitex_mpi_comm,ierror)
  call MPI_BARRIER(amitex_mpi_comm,ierror) 
  call MPI_Allreduce(nmateriaux0%n_loc,check_nbmat,1,MPI_INTEGER,MPI_MAX,amitex_mpi_comm,ierror)


  ! Evaluation des nombres de materiaux non-composites et totaux sur le pinceau
  nmateriaux0%n_loc = size(MattotP0)
  nmateriaux0%n_non_composites_loc = 0
  do i = 1, nmateriaux0%n_loc
     if(MattotP0(i)%NPhase==1) nmateriaux0%n_non_composites_loc = nmateriaux0%n_non_composites_loc + 1
  end do

  ! Deduction des nombres locaux et globaux de materiaux composites et check 
  if(nmateriaux0%n/=check_nbmat) call amitex_abort(&
   "mauvais calcul du nombre de materiaux (eval_nb_materiaux / material_mod)",2)
  check_nbmat_comp = nmateriaux0%n - nmateriaux0%n_non_composites
  if(nmateriaux0%n_composites/=check_nbmat_comp) call amitex_abort(&
   "mauvais calcul du nombre de materiaux composites (eval_nb_materiaux / material_mod)",2)
  nmateriaux0%n_composites_loc = nmateriaux0%n_loc - nmateriaux0%n_non_composites_loc
  
  ! Evaluate the number of zone per material 
  allocate(nmateriaux0%nz(nmateriaux0%n))
  do i = 1,nmateriaux0%n
    do j = 1,size(mattotP0)
      nZone=0
      if (mattotP0(j)%numM == i) nZone = maxval(mattotP0(j)%Zone(:,2))
    end do
    call MAX_MPI_I8(nZone, nmateriaux0%nz(i),ierror)
  end do

end subroutine eval_nb_materiaux

!!===================================================================
!!                 SPREAD_PERSO
!!                 --------------------------
!!
!!  Procedure mise en place pour contourner un bug de la fonction spread 
!!       avec intel-oneapi 2021.1 avec les voxels composites
!!
!!
!!  array(:,i) = array1(:)
!!
!!  \param[in]    array1, dimension(:)
!!  \param[inout] array2, dimension(:,:)
!!===================================================================

subroutine spread_perso(array2, array1)

   real(mytype), dimension(:,:), intent(inout) :: array2
   real(mytype), dimension(:), intent(in)      :: array1

   integer :: i
      
   !  check size   
   if (size(array2,1) /= size(array1,1)) then
   call amitex_abort(&
   "function spread_perso in material_mod : size(array2,1) /= size(array1,1)",2)
   end if
   
   do i = 1,size(array2,2)
        array2(:,i) = array1(:)
   end do
   
end subroutine spread_perso

end module material_mod
