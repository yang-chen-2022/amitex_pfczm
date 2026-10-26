!====================================================================
!
!       MODULE PARAM_ALGO_MOD
!
!> Definition, initialisation des parametres de l'algorithme
!!
!!
!!  Structure \n
!! PARAM_ALGO :    parametres de l'algorithme
!!
!!  Subroutine \n
!! read_param :    lecture des parametres de l'algorithme dans le fichier xml
!!
!====================================================================

module param_algo_mod

  use, intrinsic :: IEEE_ARITHMETIC

#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi,       only : nrank
#else
  use decomp_2d, only : nrank,mytype 
#endif  
  use Fox_dom
  use io_amitex_mod
  use error_mod

  private

  !> Pointer "publiques" vers composantes de SIMU_AMITEX
  public :: algo_param, user_param, Nloc_param

  !> Variables "publiques" utilisees lors de l'initialisation
  public :: algo_param0, user_param0, Nloc_param0

  !> Types publiques (pour definition de SIMU_AMITEX)
  public :: PARAM_ALGO, PARAM_ALGO_USER,PARAM_ALGO_NLOC

  !> Fonctions publiques
  public :: print_param_algo,print_param_algo_composite,&
            read_param, read_param_composite

!====================================================================
!> Structure des parametres de l'algorithme
!--------------------------------------------------------------------
  type PARAM_ALGO
     !>--------- COMMUN MECANIQUE - DIFFUSION
     !> Nombre d'iterations max
     integer                                 :: nIterMax=1000
     !> Nombre d'iterations min (0 ou 1) 
     integer                                 :: nIterMin=1
     !> Nombre d'iterations min (0 ou 1) apres une acceleration de convergence
     integer                                 :: nIterMin_acv=1
     !> Acceleration de convergence
     logical                                 :: acc_CV 
     integer                                 :: modACV=3          ! acceleration every modACV iterations (default=3)
     !> Tolerance sur le critere d'equilibre
     real(mytype)                            :: tol_criteq=1e-4           !valeur par defaut (>1e-12 et <1e-3) 
     !> Tolerance pour les criteres en deformation (critere de compatibilite 
     !! et critere sur la deformation moyenne)
     real(mytype)                            :: tol_critdef=1e-10_mytype  !valeur par defaut (>1e-14 et <1e-6)
     !> Tolerance pour le critere en contrainte moyenne
     real(mytype)                            :: tol_critsig=1e-4_mytype   !valeur par defaut (>1e-12 et <0.1)
     !> Type de schema
     character(len=16)                       :: scheme_type 
     !> Initialisation de la deformation en debut de pas
     character(len=16)                       :: init_def = 'default'    ! linear extrapolation (Def0,Def)_old -> (Def0,Def)_new
                                                                        ! si 'previous' (initialisation avec la solution du pas précédent)
                                                                        !                Def0 = Def and Def = Def 
                                                                        ! si 'None', nothing done : allows to arbitrary define Def0 and Def                                                              
                                                                        !             definie dans une fonction utilisateur

     !> Sub-stepping 
     logical                                 :: substepping = .false.
     !> Sub-stepping division
     integer                                 :: Nsub  = 10         ! if Nitermax2 eached, 
     !> Sub-stepping - max number of iteration (used instead of algo_param%NIterMax if load_incr /= 0)
     integer                                 :: nIterMax2  = 10
     real(mytype)                            :: geom_ratio = 2.   ! for depth = i : Nitmax = nItermax2 * (geom_ratio)^i
     !> Sub-Stepping depth
     integer                                 :: depth = 5

     !> Prise en compte de la Mécanique
     logical                                 :: Mechanics
     !> Prise en compte de la Diffusion
     logical                                 :: Diffusion
     !> Prise en compte de la convergence forcee
     logical                                 :: CVfor = .false.
     integer                                 :: Nit_cvfor = 200         ! nbre d'iterations avant de faire une convergence forcee
     integer                                 :: Ncvfor = 100            ! nbre de convergences forcees avant sortie programme
     character(len=16)                       :: init_cvfor = 'default'  ! initialisation de la deformation pour chaque iteration de CVFOR :
                                                                        ! 'default' : initialisation identique au cas sans CVFOR
                                                                        ! 'last' : dernier champ de l'iteration CVFOR precedente
                                                                        ! 'best' : meilleur champ de l'iteration CVFOR precedente

     !>--------- SPECIFIQUE MECANIQUE
     !> Rayon du filtre (operateur de Green)
     real(mytype)                            :: filter_radius
     !> Type de filtre (operateur de Green)
     character(len=16)                       :: filter_type
     !> Hypothese des petites perturbations
     logical                                 :: HPP
     !> Hypothese des petites perturbations non symmetrique : 
     !! active le calcul du tenseur gradient du deplacement 
     !! dans le cas HPP, ainsi que son utilisation au sein 
     !! des lois de comportement UMAT. 
     !! ATTENTION : Le gradient de deplacement moyen est impose symetrique
     logical                                 :: HPP_nsym = .false.
     !> Nombre de composantes pour le tenseur de deformation (6 or 9)
     integer                                 :: nTensDef
     !> Correction pour les dimensions paires (.true. uniquement dans le cas "non_filtre")
     logical                                 :: CorrPair = .false.
     !> Tenseur d'elasticite du materiau de reference en Grandes Transformations
     logical                                 :: C0sym = .true.
     
     !>--------Specifique MECANIQUE NON LOCALE
     !> Modelisation non locale utilisee
     logical                                 :: Nloc
     logical,allocatable,dimension(:)        :: nloc_explicit   ! .true. if Nloc_models(i) is implemented explicitly around 
                                                                !           the local mechanics algorithm
                                                                ! default is .true. 
     logical,allocatable,dimension(:)        :: nloc_implicit   ! .not. nloc_explicit
                                                                !       implicit introduction is still under progress  

     !>--------- SPECIFIQUE DIFFUSION
     !> Rayon du filtre (operateur de Green)
     real(mytype)                            :: filter_radiusD
     !> Type de filtre (operateur de Green)
     character(len=16)                       :: filter_typeD
     !> Type de schema
     character(len=16)                       :: scheme_typeD 
     !> Nombre de variables pour la diffusion (1 pour commencer!)
     integer                                 :: nVarD
     !> Correction pour les dimensions paires (.true. uniquement dans le cas "non_filtre")
     logical                                 :: CorrPairD = .false.
     !> Type de calcul (stationnaire)
     logical                                 :: StationaryD


     !>--------- SPECIFIQUE ALGO COMPOSITE LAMINATE/REUSS
     !> Acceleration de convergence
     logical                                 :: acc_CV_laminate 
     logical                                 :: acc_CV_reuss 
     !> Tolerance sur le critere d'equilibre (a multiplier par le critere global tol_criteq)
     !>                    (par defaut 0.01)
     real(mytype)                            :: tol_criteq_laminate 
     real(mytype)                            :: tol_criteq_reuss 
     !> Type d'initialisation (default,proportional, linear)
     !> "default" : increments de deformation homogene (identique sur toutes les phases)
     !> "proportional" : inrements de def. conservant la proportionalite des deformations entre phases
     !> "linear" : increments de deformation interpoles a partir du pas precedent
     character(len=16)                       :: Init_laminate
     character(len=16)                       :: Init_reuss
     !> Paramètres de pilotage
     integer(kind=4)                         :: Nmax_subdivision_laminate
     integer(kind=4)                         :: Nmax_subdivision_reuss
     !> Type de matrice tangente utilisé (elastic, mfront_elastic, mfront_tangent, numerical)
     !> "elastic"        : comportement linéaire isotrope par phase donne par Coeff_composite dans materiau.xml
     !> "mfront_elastic" : tenseur d'élasticité renvoyee par Mfront
     !> "mfront_tangent" : matrice tangente cohérente renvoyee par Mfront
     !> "mfront_secant"  : matrice secant (elastique endommagee) renvoyee par Mfront 
     !                     si codee dans la loi
     !> "numerical"     : matrice tangente obtenue par differentiation numerique
     character(len=16)                       :: Jacobian_type_laminate
     character(len=16)                       :: Jacobian_type_reuss
     !> Valeur absolue de la perturbation pour le calcul de la matrice tangente numerique
     !  du modèle laminate
     real(mytype)                            :: Perturbation_laminate
     real(mytype)                            :: Perturbation_reuss


  end type PARAM_ALGO
!--------------------------------------------------------------------

!> Structure de parametres USER de l'algorithme
!--------------------------------------------------------------------
  type PARAM_ALGO_USER
    logical                                           :: test=.false. ! .true. if using user defined functions (2 possibilities)   
    character(len=10)                                 :: algo         ! "standard" : user functons in the standard algorithm
                                                                      ! "user" : a complete algorihm
    real(mytype), allocatable, dimension(:)           :: p_real       ! real parameters
    character(len=1000), allocatable, dimension(:)    :: p_string     ! string parameters           
  end type PARAM_ALGO_USER

!> Structure de parametres pour l'utilisation de modeles NLOC
!--------------------------------------------------------------------
  type PARAM_ALGO_NLOC 
    !integer                                           :: NlocMod_num    ! rendu coherent avec l'indice du tableau Nloc_models  
    real(mytype), allocatable, dimension(:)           :: p_real          ! real parameters
    character(len=1000), allocatable, dimension(:)    :: p_string        ! string parameters           
  end type PARAM_ALGO_NLOC

!!------------------------------------------------------------------------------
!>                                                                  DECLARATIONS

!> Derived types defined above, for the STANDARD algorithm
type(PARAM_ALGO)                                  :: algo_param0
type(PARAM_ALGO),pointer                          :: algo_param
type(PARAM_ALGO_NLOC),allocatable,dimension(:)    :: Nloc_param0
type(PARAM_ALGO_NLOC),pointer,dimension(:)        :: Nloc_param

!> Derived types defined above, for a USER defined algorithm
type(PARAM_ALGO_USER)                             :: user_param0
type(PARAM_ALGO_USER),pointer                     :: user_param

contains

!==================================================================================
!==================================================================================
!                         SUBROUTINE PRINT_PARAM_ALGO
!
!> Ecriture des composantes d'une variable de type param_algo
!!
!!
!!  \param[in] Flog:   (entier) unite logique du fichier de sortie .log
!!  \param[in] nrank0: (entier) pinceau pour lequel on effectue la sortie
!!
!==================================================================================
subroutine print_param_algo(Flog,nrank0)
  
  implicit none
  integer, intent(in)   :: Flog
  integer, intent(in)   :: nrank0
  integer               :: i,k


  if(nrank==nrank0)then 
  write(Flog,"(A)") " "
  write(Flog,"(A)") "STRUCTURE PARAM_ALGO "
  write(Flog,"(A)") "====================="
  write(Flog,"(A)") " "

  write(Flog,"(A)") "Commun MECANIQUE-DIFFUSION"
  write(Flog,"(2A)") "Algorithme : ",algo_param0%scheme_type 
  write(Flog,"(A,E15.8)") &
       "Critere d'equilibre : ", algo_param0%tol_criteq
  write(Flog,"(A,E15.8)") &
       "Critere sur la contrainte moyenne : ", algo_param0%tol_critsig
  write(Flog,"(A,E15.8)") &
       "Critere en deformation : ", algo_param0%tol_critdef
       
  if(algo_param0%acc_CV) then
     write(Flog,"(A)") "Simulation with convergence acceleration."
     write(Flog,"(A,I6,A)") "Acceleration every ",algo_param0%modACV," iterations."
  else
     write(Flog,"(A)") "Simulation without convergence acceleration."
  end if
  write(Flog,"(A,I6)") &
       "Max number of iterations : ", algo_param0%nItermax
  write(Flog,"(A,I6)") &
       "Min number of iterations : ", algo_param0%nItermin
  write(Flog,"(A,I6)") &
       "Min number of iterations after an acceleration (0 or 1) : ", algo_param0%nItermin_acv

  if(algo_param0%substepping) then
       write(Flog,"(A)") "Simulation with SUBSTEPPING."
       write(Flog,"(A,I6)") "   Division parameter (Nsub)",algo_param0%Nsub
       write(Flog,"(A,I6)") "   Max. number of iterations before sub-step (Nitermax2)",algo_param0%Nitermax2
       write(Flog,"(A,E15.8)") "   Geometrical ratio to increase Max number of iterations (Nitermax2*(geom_ratio)^depth)" &
                                                                                      ,algo_param0%geom_ratio
       write(Flog,"(A,I6)") "   Depth (Depth) ",algo_param0%depth
  else
       write(Flog,"(A)") "Simulation without SUBSTEPPING."  
  end if

  if(algo_param0%CVfor) then
     write(Flog,"(A)") "Calcul effectue avec convergence forcee."
     write(Flog,"(A,I5)") "Nombre d'iterations avant de forcer : ", algo_param0%Nit_cvfor
     write(Flog,"(A,I5)") "Nombre de convergences forcees avant sortie du programme : ", algo_param0%Ncvfor
     write(Flog,"(A,A)") "Initilisation de la deformation a chaque CV forcee : ", trim(algo_param0%init_cvfor)
  else
     write(Flog,"(A)") "Calcul effectue sans convergence forcee."
  end if 
  
  write(Flog,"(A)") " "
  write(Flog,"(A,A)") &
       "Initilisation ('default' ou 'previous') : ", algo_param0%init_def

  if(algo_param0%Mechanics) then
     write(Flog,"(A)") " "
     write(Flog,"(A)") "Specifique MECANIQUE"
     write(Flog,"(3A,E15.8)") &
          "Filtre de l'operateur de Green : ", algo_param0%filter_type," R = ", algo_param0%filter_radius
     if(algo_param0%HPP) then
        write(Flog,"(A)") "Hypothese des petites perturbations."
        if (algo_param0%HPP_nsym) then
            write(Flog,"(A)") "Prise en compte complete du gradient du deplacement dans le calcul HPP (9 composantes)"
        end if
     else
        write(Flog,"(A)") "Calcul effectue en grandes transformations."
        if (algo_param0%C0sym) then
           write(Flog,"(A)") "   avec un milieu de reference C0 symetrique"
        else
           write(Flog,"(A)") "   avec un milieu de reference C0 non-symetrique"
        end if
     end if
     if(algo_param0%CorrPair) then
        write(Flog,"(A)") "Correction specifique pour dimensions de grille paires"
     end if

     !----------------NLOC
     if (algo_param0%Nloc) then
        write(Flog,"(A)") " "
        write(Flog,"(A)") "Specifique MECANIQUE NON LOCALE"
        do i=1,size(Nloc_param0)
           write(Flog,"(A,I0)") &
                 "Non local model ", i

           ! implementation type
           if(algo_param0%nloc_explicit(i)) then
              write(Flog,"(A)")"   Explicit implementation in the standard alorithm"
           else
              write(Flog,"(A)")"   Implicit implementation in the FFT resolution (NOT YET AVAILABLE)"              
           end if

           ! real parameters
           if(allocated(Nloc_param0(i)%p_real)) then
           do k=1,size(Nloc_param0(i)%p_real)
              write(Flog,"(A,I0,A,E15.8)") "   Real Parameter   ",k, " = ", Nloc_param0(i)%p_real(k)
           end do
           end if

           ! string parameters
           if(allocated(Nloc_param0(i)%p_string)) then
           do k=1,size(Nloc_param0(i)%p_string)
              write(Flog,"(A,I0,A,A)") "   String Parameter ",k, " = ", trim(Nloc_param0(i)%p_string(k))
           end do
           end if

           write(Flog,"(A)")" "
        end do
     end if
     !----------------END NLOC

  end if

   

  if(algo_param0%Diffusion) then
     write(Flog,"(A)") " "
     write(Flog,"(A)") "Specifique DIFFUSION"
     if (algo_param0%StationaryD) then
          write(Flog,"(A)") "Calcul stationnaire"
     else
          write(Flog,"(A)") "Calcul non-stationnaire"
     end if
     write(Flog,"(3A,E15.8)") &
          "Filtre de l'operateur de Green : ", algo_param0%filter_typeD," R = ", algo_param0%filter_radiusD
     write(Flog,"(A,I4)") "Nombre de variables de diffusion : ", algo_param0%nVarD
     if(algo_param0%CorrPairD) then
        write(Flog,"(A)") "Correction specifique pour dimensions de grille paires"
     end if
     write(Flog,"(A)") " "
  end if

  

  ! <User> parameters
  if (user_param0%test) then
     write(Flog,"(A)") " "
     write(Flog,"(A)") "Specific USER"

     ! algorithm
     write(Flog,"(A)") " "
     write(Flog,"(A)") " Algorithm "//trim(user_param0%algo)
     write(Flog,"(A)") " "

     ! real parameters
     if(allocated(user_param0%p_real)) then
        do k=1,size(user_param0%p_real)
           write(Flog,"(A,I0,A,E15.8)") "   Real Parameter   ",k, " = ", user_param0%p_real(k)
        end do
     end if

     ! string parameters
     if(allocated(user_param0%p_string)) then
        do k=1,size(user_param0%p_string)
           write(Flog,"(A,I0,A,A)") "   String Parameter ",k, " = ", trim(user_param0%p_string(k))
        end do
     end if
  end if 

  end if !nrank==0

  !----------------------------------------------------CHECK CONSISTENCY
  
  if(algo_param0%Diffusion .and. algo_param0%Mechanics) then
      if ((algo_param0%filter_typeD .ne. algo_param0%filter_type) &
          .or. (algo_param0%filter_radiusD .ne. algo_param0%filter_radius)) then
          call amitex_abort("parametres de filtres differents (mecanique/diffusion)",2,0)
      end if
  end if 
  
  if (algo_param0%substepping .and. algo_param0%CVfor) then
       call amitex_abort("Forced-Convergence and Sub-Stepping procedures cannot be used together",2,0)    
  end if

end subroutine print_param_algo

!==================================================================================
!==================================================================================
!                         SUBROUTINE PRINT_PARAM_ALGO_COMPOSITE
!
!> Ecriture des composantes specifiques a l'utilisation de voxels composites
!! d'une variable de type param_algo
!!
!!
!!  \param[in] Flog:   (entier) unite logique du fichier de sortie .log
!!  \param[in] nrank0: (entier) pinceau pour lequel on effectue la sortie
!!
!==================================================================================
subroutine print_param_algo_composite(Flog,nrank0)
  
  implicit none
  integer, intent(in)   :: Flog
  integer, intent(in)   :: nrank0

  if(nrank==nrank0)then 
  write(Flog,"(A)") " "
  write(Flog,"(A)") "Specifique VOXELS COMPOSITES LAMINATE"
  if (algo_param0%acc_CV_laminate) then
     write(Flog,"(A)") "   Algorithme : Newton modifie (raideur constante) avec acceleration de convergence"
  else
     write(Flog,"(A)") "   Algorithme : Newton modifié (raideur constante)"
  end if
  write(Flog,"(A,E15.8)") "   Critere de convergence (a multiplier par le critere global) = ", &
                          algo_param0%tol_criteq_laminate
  write(Flog,"(A,L7)")    "   Acceleration de convergence : ", algo_param0%acc_CV_laminate
  write(Flog,"(A,A)")     "   Type d'initialisation : ", algo_param0%Init_laminate
  write(Flog,"(A,I4)")    "   Nombre max de sous-decoupage : ", algo_param0%Nmax_subdivision_laminate
  write(Flog,"(A,A)")     "   Type de matrice tangente utilise : ", algo_param0%Jacobian_type_laminate
  write(Flog,"(A,E15.8)") "   Perturbation si tangente numerique : ", algo_param0%Perturbation_laminate
 
  write(Flog,"(A)") " "
  write(Flog,"(A)") "Specifique VOXELS COMPOSITES REUSS"
  if (algo_param0%acc_CV_reuss) then
     write(Flog,"(A)") "   Algorithme : Newton modifie (raideur constante) avec acceleration de convergence"
  else
     write(Flog,"(A)") "   Algorithme : Newton modifié (raideur constante)"
  end if
  write(Flog,"(A,E15.8)") "   Critere de convergence (a multiplier par le critere global) = ", &
                          algo_param0%tol_criteq_reuss
  write(Flog,"(A,L7)")    "   Acceleration de convergence : ", algo_param0%acc_CV_reuss
  write(Flog,"(A,A)")     "   Type d'initialisation : ", algo_param0%Init_reuss
  write(Flog,"(A,A)")     "   Type de matrice tangente utilise : ", algo_param0%Jacobian_type_reuss

!TODO prendre en compte ces possibilites pour les voxels composites REUSS
!  write(Flog,"(A,I4)")    "   Nombre max de sous-decoupage : ", algo_param0%Nmax_subdivision_reuss

!  write(Flog,"(A,E15.8)") "   Perturbation si tangente numerique : ", algo_param0%Perturbation_reuss


  select case (trim(algo_param0%Init_reuss))
  case("default","linear","proportional")
  case default
      call amitex_abort("algo_param0%Init_reuss different de ""default"",&
               & ""linear"" ou ""proportional"" (print_param_algo_composite)",2,0)
  end select

  if (trim(algo_param0%Jacobian_type_reuss) .ne. "elastic") then
   call amitex_abort("algo_param0%Jacobian_type_reuss different de ""elastic"" (print_param_algo_composite)",2,0)
  end if

  end if

end subroutine print_param_algo_composite

!====================================================================
!====================================================================
!                         SUBROUTINE READ_PARAM
!
!>   Lecture des parametres de l'algorithme dans un fichier xml
!!
!! Plusieurs parametres sont lus :
!!   - La tolerance pour le critere d'equilibre (par defaut 1e-4 et doit
!!     etre superieure a 1e-3).
!!   - L'acceleration de convergence (vaut default, true ou false).
!!     La valeur par defaut est true.
!!   - Le type de filtre (default, no_filter, octa ou hexa).
!!     La valeur par defaut est hexa.
!!   - Le schema (pour l'instant le seul schema implemente est
!!     basic_scheme)
!! \param[in]   file_algo: nom du fichier xml (chaine de caracteres)
!!
!====================================================================
subroutine read_param(file_algo)

  implicit none
  
  character(len=*), intent(in) :: file_algo
  type(node) ,pointer          :: fi, cur_node, cur_node0,cur_elem
  type(nodeList), pointer      :: node_list, node_list0, elem_list, userParam
  character(len=200)           :: string, tmp_string,err_msg
  !> Parametres par defaut (meca)
  real(mytype), parameter      :: max_tol_criteq = (10._mytype)**(-3)
#ifdef DOUBLE_PREC
  real(mytype), parameter      :: min_tol_criteq = (10._mytype)**(-12)
#else
  real(mytype), parameter      :: min_tol_criteq = (10._mytype)**(-5)
#endif
  character(len=25), parameter :: default_filter = "hexa", default_scheme = "basic_scheme"
  character(len=25), parameter :: default_filterD = "hexa"
  integer                      :: i,j,alloc_stat,ind,ind0

  integer                      :: nuser_param
  real(mytype)                 :: crit_default


  fi => parseFile(trim(file_algo))

!!---------------------LECTURE DES PARAMETRES COMMUN "MECANIQUE-DIFFUSION"

  node_list0 => getElementsByTagName(fi,"Algorithm")
  if(getLength(node_list0)==1) then
    cur_node0 => item(node_list0,0)
  else 
    call amitex_abort("Le noeud 'Algorithm' n'est pas renseigne dans le fichier xml (read_param)",1,0)
  end if

  !! On recupere le type de schema utilise
  string = "Type"
  call get_str_xml_default(cur_node0, trim(string), algo_param0%scheme_type,0,"basic_scheme",1)
  if((algo_param0%scheme_type /= "basic_scheme") .AND. (algo_param0%scheme_type /= "user")) then
     call amitex_abort("Type de schema inconnu (read_param)",1,0)
  end if

  !! On recupere le critere de convergence demande (si non renseigne, valeur par defaut)
  string = "Convergence_Criterion"
  elem_list => getElementsByTagName(cur_node0, string)
  if(getLength(elem_list)>0) then
     cur_elem => item(elem_list, 0)       
     string = "Value"
     !ici la valeur par defaut est algo_param0%tol_criteq, initilalisee avec algo_param0
     crit_default = algo_param0%tol_criteq
     call get_real_xml_default(cur_elem, string, algo_param0%tol_criteq,0,crit_default,1)
  end if
  if(algo_param0%tol_criteq > max_tol_criteq) then
     call amitex_abort("Critere d'arret de l'algorithme trop eleve (read_param)",2,0)
  elseif(algo_param0%tol_criteq < min_tol_criteq) then
     call amitex_abort("Critere d'arret de l'algorithme trop faible (read_param)",2,0)
  end if

  !! On recupere le type d'initialisation ('default', par defaut ou 'previous')
  string = "Initialize"
  elem_list => getElementsByTagName(cur_node0, string)
  if(getLength(elem_list)>0)then
     cur_elem => item(elem_list, 0)       
     string = "Value"
     call get_str_xml(cur_elem, trim(string), tmp_string,1)
     algo_param0%init_def = trim(tmp_string)
     if((algo_param0%init_def .NE. "previous") .AND. (algo_param0%init_def .NE. "default")) then
        call amitex_abort(" 'Initialize' must be set to 'previous' or 'default' (read_param)",2,0)
     end if
  end if


  !! On recupere le critere de compatibilite (si non renseigne, valeur par defaut)
  string = "Convergence_Criterion_Compatibility"
  elem_list => getElementsByTagName(cur_node0, string)
  if(getLength(elem_list)>0) then
     cur_elem => item(elem_list, 0)       
     string = "Value"
     !ici la valeur par defaut est algo_param0%tol_critdef, initilalisee avec algo_param0
     crit_default = algo_param0%tol_critdef
     call get_real_xml_default(cur_elem, string, algo_param0%tol_critdef,0,crit_default,1)
  end if
  if(algo_param0%tol_critdef > 1.e-6_mytype) then
     call amitex_abort("STOP : compatibility criterion too high (> 1e-6) (read_param)",2,0)
  elseif(algo_param0%tol_critsig < 1.e-14) then
     call amitex_abort("STOP : compatibility criterion too low (<1e-14) (read_param)",2,0)
  end if


  !! On recupere le critere de convergence demande (si non renseigne, valeur par defaut)
  string = "Convergence_Criterion_Smacro"
  elem_list => getElementsByTagName(cur_node0, string)
  if(getLength(elem_list)>0) then
     cur_elem => item(elem_list, 0)       
     string = "Value"
     !ici la valeur par defaut est algo_param0%tol_critsig, initilalisee avec algo_param0
     crit_default = algo_param0%tol_critsig
     call get_real_xml_default(cur_elem, string, algo_param0%tol_critsig,0,crit_default,1)
  end if
  if(algo_param0%tol_critsig > 0.1_mytype) then
     call amitex_abort("Critere Smacro trop eleve (read_param)",2,0)
  elseif(algo_param0%tol_critsig < min_tol_criteq) then
     call amitex_abort("Critere Smacro trop faible (read_param)",2,0)
  end if

  !! On recupere le nombre d'iteration max (si non renseigne, valeur initiale par defaut 1000)
  string = "Nitermax"
  elem_list => getElementsByTagName(cur_node0, string)
  if(getLength(elem_list)>0) then
     cur_elem => item(elem_list, 0)       
     string = "Value"
     call get_int_xml(cur_elem,string,algo_param0%nItermax,1)
  end if

  !! On recupere le nombre d'iteration min
  string = "Nitermin"
  elem_list => getElementsByTagName(cur_node0, string)
  if(getLength(elem_list)>0) then
     cur_elem => item(elem_list, 0)       
     string = "Value"
     call get_int_xml(cur_elem,string,algo_param0%nItermin,1)
  end if
  if ((algo_param0%nItermin .ne. 0) .AND. (algo_param0%nItermin .ne. 1)) then
        call amitex_abort("Nitermin doit valoir 0 ou 1 (read_param)",1,0)
  end if

  !! On recupere le nombre d'iteration min
  string = "Nitermin_acv"
  elem_list => getElementsByTagName(cur_node0, string)
  if(getLength(elem_list)>0) then
     cur_elem => item(elem_list, 0)       
     string = "Value"
     call get_int_xml(cur_elem,string,algo_param0%nItermin_acv,1)
  end if
  if ((algo_param0%nItermin_acv .ne. 0) .AND. (algo_param0%nItermin_acv .ne. 1)) then
        call amitex_abort("Nitermin_acv doit valoir 0 ou 1 (read_param)",1,0)
  end if

  !! On verifie si l'utilisateur veut utiliser l'acceleration de convergence
  string = "Convergence_Acceleration"
  elem_list => getElementsByTagName(cur_node0, string)
  if(getLength(elem_list)>0)then
     !Value
     cur_elem => item(elem_list, 0)       
     string = "Value"
     call get_str_xml(cur_elem, trim(string), tmp_string,1)
     if(tmp_string == "true") then
        algo_param0%acc_CV = .TRUE.
     elseif(tmp_string == "false") then
        algo_param0%acc_CV = .FALSE.           
     else
        call amitex_abort("Acceleration de convergence doit etre true ou false (read_param)",1,0)
     end if
     
     !modACV
     string = "modACV"
     call get_int_xml(cur_elem, string, algo_param0%modACV,-2)

  else 
     call amitex_abort("Acceleration de convergence non renseignee (read_param)",1,0)
  end if

  !! On verifie si l'utilisateur veut utiliser sub-stepping
  string = "Substepping"
  elem_list => getElementsByTagName(cur_node0, string)
  if(getLength(elem_list)>0)then

     algo_param0%substepping = .true.
     cur_elem => item(elem_list, 0)       
   
     !Nitermax2
     string = "Nitermax2"
     call get_int_xml(cur_elem, string, algo_param0%Nitermax2,-2)
     
     !Geom_ratio
     string = "Geom_ratio"
     !call get_int_xml(cur_elem, string, algo_param0%geom_ratio,-2)
     call get_real_xml(cur_elem, string,algo_param0%geom_ratio,-2)
     
     !Nsub
     string = "Nsub"
     call get_int_xml(cur_elem, string, algo_param0%Nsub,-2)
     
     !Depth
     string = "Depth"
     call get_int_xml(cur_elem, string, algo_param0%Depth,-2)

  else 
     call amitex_abort("No SUBSTEPPING (read_param)",-2,0)
  end if

  !! On verifie si l'utilisateur veut utiliser la convergence forcee
  string = "Convergence_Forced"
  elem_list => getElementsByTagName(cur_node0, string)
  if(getLength(elem_list)>0)then
     !Value
     cur_elem => item(elem_list, 0)       
     string = "Value"
     call get_str_xml(cur_elem, trim(string), tmp_string,1)
     if(tmp_string == "true") then
        algo_param0%CVfor = .TRUE.
     elseif(tmp_string == "false") then
     else
        call amitex_abort("Convergence_Forced doit etre true ou false (read_param)",1,0)
     end if

     !Nit_cvfor
     string = "Nit_cvfor"
     call get_int_xml(cur_elem, string, algo_param0%Nit_cvfor,1)

     !Ncvfor
     string = "Ncvfor"
     call get_int_xml(cur_elem, string, algo_param0%Ncvfor,1)

     !Init_cvfor
     string = "Init_cvfor"
     call get_str_xml(cur_elem, trim(string), algo_param0%init_cvfor,-1,0)
     if((algo_param0%init_cvfor .NE. "default") .AND. (algo_param0%init_cvfor .NE. "last") &
        .AND. (algo_param0%init_cvfor .NE. "best")) then
        call amitex_abort("Init_cvfor mal renseigne : 'default', 'best' ou 'last' (read_param)",1,0)
     end if

  end if



!!---------------------LECTURE DES PARAMETRES POUR L'ALGO "MECANIQUE"

  node_list0 => getElementsByTagName(fi,"Mechanics")
  if(getLength(node_list0)==1) then
    algo_param0%Mechanics=.true.
    cur_node0 => item(node_list0,0)
  elseif(getLength(node_list0)==0) then
    algo_param0%Mechanics=.false.
  elseif(getLength(node_list0)>1) then
    algo_param0%Mechanics=.false.
    call amitex_abort("Le noeud 'Mechanics' est present plus d'une fois (read_param)",1,0)
  end if

 if (algo_param0%Mechanics .eqv. .true.) then

    !! On recupere le filtre et on fixe le rayon par default
    node_list => getElementsByTagName(cur_node0,"Filter")
    if(getLength(node_list)==1)then
       cur_node => item(node_list,0)
       !! lecture du type de filtre (et affectation a defaut_filter si on lit "default", cf get_str_xml_default)
       string = "Type"
       call get_str_xml_default(cur_node,trim(string),algo_param0%filter_type,0,default_filter,1)     
       !! rayon du filtre en fonction du type
       if(algo_param0%filter_type == "hexa")then
          algo_param0%filter_radius = 1
       elseif(algo_param0%filter_type == "octa")then
          algo_param0%filter_radius = 2
       elseif(algo_param0%filter_type == "no_filter") then
          algo_param0%filter_radius = 0
          algo_param0%CorrPair = .true.
       else
          call amitex_abort("Filtre inconnu (read_param)",1,0)
       end if
    else
       call amitex_abort("Filtre non renseigne ou renseigne plusieurs fois (read_param)",1,0)
    end if

    !! On recupere le type de calcul : HPP ou grandes transformations
    node_list => getElementsByTagName(cur_node0,"Small_Perturbations")
    if(getLength(node_list)==1)then
       cur_node => item(node_list,0)
       !! On recupere la valeur du booleen
       string = "Value"
       call get_str_xml(cur_node, trim(string), tmp_string,1)
       if(tmp_string == "true") then
          algo_param0%HPP = .TRUE.
          algo_param0%nTensDef = 6
       elseif(tmp_string == "false") then
          algo_param0%HPP = .FALSE.           
          algo_param0%nTensDef = 9
       else
          call amitex_abort("In Small_Perturbation, Value must be &
               &'true' or 'false' (read_param)",1,0)
       end if
       
       !! On recupere la valeur du booleen
       string = "Displacement_Gradient"
       tmp_string="sym"
       call get_str_xml(cur_node, trim(string), tmp_string,-2,0)
       if(tmp_string == "sym") then
              algo_param0%HPP_nsym = .FALSE.
       elseif(tmp_string == "nsym") then
              algo_param0%HPP_nsym = .TRUE.
       else
          call amitex_abort("In Small_Perturbation, Displacement_Gradient must be &
               &'sym' or 'nsym' (read_param)",1,0)
       end if

    else
       call amitex_abort("Hypothese des petites perturbations non renseignee ou renseignee plusieurs&
            & fois. Les calculs seront faits par la suite en faisant cette hypothese (read_param)"&
            ,-1,0)
       algo_param0%HPP = .TRUE.
       algo_param0%nTensDef = 6
    end if    

    !! On recupere le choix d'utiliser un C0 symetrique ou non (en grande transformation)
    node_list => getElementsByTagName(cur_node0,"C0sym")
    if(getLength(node_list)==1)then
       cur_node => item(node_list,0)
       !! On recupere la valeur du booleen
       string = "Value"
       call get_str_xml(cur_node, trim(string), tmp_string,1)
       if(tmp_string == "true") then
          algo_param0%C0sym = .TRUE.
       elseif(tmp_string == "false") then
          algo_param0%C0sym = .FALSE.           
       else
          call amitex_abort("C0sym doit etre &
               &'true' ou 'false' (read_param)",1,0)
       end if
       if (algo_param0%HPP) then
          call amitex_abort("Calcul en petites perturbations, &
                         & renseigner C0sym est inutile (read_param)",-1,0)
       end if
    else
       if (.not. algo_param0%HPP) then
        call amitex_abort("Calcul en grandes transformations avec C0 symetrique (read_param)",-1,0)
       end if
    end if

 end if !fin de la condition sur algo_param0%Mechanics

  !!---------------------LECTURE DES PARAMETRES POUR LA "MECANIQUE NON LOCALE" 

  !! Les indices NLocMod_num correspondent aux indices de Nloc_param0 

  node_list => getElementsByTagName(fi,"Non_local_algorithm")

  if(getLength(node_list)>0)then
     algo_param0%Nloc = .TRUE.

     allocate(Nloc_param0(getLength(node_list)),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort(&
                      "Espace memoire disponible insuffisant (read_param, Non_local_modeling)",2,0)
     allocate(algo_param0%nloc_explicit(getLength(node_list)),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort(&
                      "Espace memoire disponible insuffisant (read_param, Non_local_modeling)",2,0)
     allocate(algo_param0%nloc_implicit(getLength(node_list)),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort(&
                      "Espace memoire disponible insuffisant (read_param, Non_local_modeling)",2,0)
     ! default values
     algo_param0%nloc_explicit =.true.
     algo_param0%nloc_implicit =.false.

     do i=1,getLength(node_list)

          cur_node => item(node_list,i-1)

          !! Recuperation du numero d'identification du modele non local
          string="NLocMod_num"
          call get_int_xml(cur_node, string, ind0,1)

          !! Recuperation du type d'integration - explicite par defaut
          string="Algo"
          call get_str_xml(cur_node, trim(string), tmp_string,1)

          if (trim(tmp_string) == "implicit") then       ! nloc_implicit est traite apres avec .not. nloc_explicit
             algo_param0%nloc_explicit(ind0)=.false.
             call amitex_abort("Non_local_algorithm -> Algo : 'implicit' is not yet available",1,0)  
          elseif (trim(tmp_string) == "explicit") then
          else
             call amitex_abort("Non_local_algorithm -> Algo : Possible values are 'implicit' &
                 & or 'explicit'. Defaul value ('explicit') is used  (read_param)",-1,0)  
          end if

          !! Recuperation de parametres utilisateurs (copier-coller-adapter de READING PARAMETERS FOR "USER" ALGORITHM)

          !> --- READING REAL PARAMETERS
          userParam=> getElementsByTagName(cur_node,"P_real")
          nuser_param=getLength(userParam)    

          if (nuser_param>0) then
            !> Allocation
            allocate(Nloc_param0(ind0)%p_real(nuser_param),stat=alloc_stat)
            if(alloc_stat /=0) call amitex_abort(&
               "Espace memoire disponible insuffisant (read_param) pour alloue user_param%p_real",2,0)
            !> lecture de toutes les valeurs de parametres (avec traitement des erreurs)
            do j=0,nuser_param-1
              cur_node0 => item(userParam,j)
              call get_int_xml(cur_node0,"Index",ind,1)
              if (ind>nuser_param .OR. ind<1) then
                 write(err_msg,fmt="(A)") &
                 " indice negatif ou nul, ou superieur au nombre de parametres P_real <user> (read_param)"
                 call amitex_abort(err_msg, 2,0)
              end if
              call get_real_xml(cur_node0,"Value",Nloc_param0(ind0)%p_real(ind),1)
            end do
          end if ! test nuser_param

          !> --- READING STRING PARAMETERS
          userParam=> getElementsByTagName(cur_node,"P_string")
          nuser_param=getLength(userParam)

          if (nuser_param>0) then
            !> Allocation
            allocate(Nloc_param0(ind0)%p_string(nuser_param),stat=alloc_stat)
            if(alloc_stat /=0) call amitex_abort(&
                  "Espace memoire disponible insuffisant (read_param) pour alloue user_param%p_string",2,0)
            !> lecture de toutes les valeurs de parametres (avec traitement des erreurs)
            do j=0,nuser_param-1
              cur_node0 => item(userParam,j)
              call get_int_xml(cur_node0,"Index",ind,0)
              if (ind>nuser_param .OR. ind<1) then
                 write(err_msg,fmt="(A)")&
                 " indice negatif ou nul, ou superieur au nombre de parametres P_string <user> (read_param)"
                 call amitex_abort(err_msg, 2,0)
              end if
              call get_str_xml_nolowercase(cur_node0,"Value",Nloc_param0(ind0)%p_string(ind),0)
            end do
          end if ! test nuser_param
     
     end do !fin liste sur les noeuds "Non_local_modeling"

     algo_param0%nloc_implicit = .not. algo_param0%nloc_explicit

  else
     algo_param0%Nloc = .FALSE.
  end if

 
!!---------------------LECTURE DES PARAMETRES POUR L'ALGO "DIFFUSION"

!TODO reflexion sur le parametres par defaut
  algo_param0%nVarD=0  

  node_list0 => getElementsByTagName(fi,"Diffusion")
  if(getLength(node_list0)==1) then
    algo_param0%Diffusion=.true.
    cur_node0 => item(node_list0,0)
  elseif(getLength(node_list0)==0) then
    algo_param0%Diffusion=.false.
  elseif(getLength(node_list0)>1) then
    algo_param0%Diffusion=.false.
    call amitex_abort("Le noeud 'Diffusion' est present plus d'une fois (read_param)",1,0)
  end if

  if (algo_param0%Diffusion .eqv. .true.) then

    !! On recupere le filtre et on fixe le rayon par default
    node_list => getElementsByTagName(cur_node0,"Filter")
    if(getLength(node_list)==1)then
       cur_node => item(node_list,0)
       !! lecture du type de filtre (et affectation a defaut_filter si on lit "default", cf get_str_xml_default)
       string = "Type"
       call get_str_xml_default(cur_node,trim(string),algo_param0%filter_typeD,0,default_filter,1)     
       !! rayon du filtre en fonction du type
       if(algo_param0%filter_typeD == "hexa")then
          algo_param0%filter_radiusD = 1
       elseif(algo_param0%filter_typeD == "octa")then
          algo_param0%filter_radiusD = 2
       elseif(algo_param0%filter_typeD == "no_filter") then
          algo_param0%filter_radiusD = 0
          algo_param0%CorrPairD = .true.
       else
          call amitex_abort("Filtre inconnu (read_param)",1,0)
       end if
    else
       call amitex_abort("Filtre non renseigne ou renseigne plusieurs fois (read_param)",1,0)
    end if

    !! On recupere le type de calcul "stationnaire" 
    node_list => getElementsByTagName(cur_node0,"Stationary")
    if(getLength(node_list)==1)then
       cur_node => item(node_list,0)
       !! On recupere le type de filtre
       string = "Value"
       call get_str_xml(cur_node, trim(string), tmp_string,1)
       if(tmp_string == "true") then
          algo_param0%StationaryD = .TRUE.
       elseif(tmp_string == "false") then
          algo_param0%StationaryD = .FALSE.           
       else
          call amitex_abort("Stationary doit etre &
               &'true' ou 'false' (read_param)",1,0)
       end if
    else
       call amitex_abort("Stationary non renseigne ou renseigne plusieurs fois (read_param)",1,0)
    end if

    !! On fixe pour l'instant le nombre de vriable a 1 
    algo_param0%nVarD=1
 
  end if !fin de la condition sur algo_param0%Diffusion


!!--------------------- READING PARAMETERS FOR "USER" ALGORITHM

  node_list0 => getElementsByTagName(fi,"User")
  if(getLength(node_list0)== 1) then  ! Reading "User" parameters
     cur_node0 => item(node_list0,0)
     user_param0%test=.true.

     !!-- READING ALGO : "standard" or "user"
     string = "Algo"
     call get_str_xml(cur_node0, trim(string), user_param0%algo,1)
     if((user_param0%algo /= "standard") .AND. (user_param0%algo /= "user")) then
        call amitex_abort("Unknown Algo in node <User> : standard or user (read_param)",1,0)
     end if

     !> --- READING REAL PARAMETERS
     cur_node => item(node_list0,0)
     userParam=> getElementsByTagName(cur_node,"P_real")
     nuser_param=getLength(userParam)    
     if (nuser_param>0) then
     !> Allocation
     allocate(user_param0%p_real(nuser_param),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_param) pour alloue user_param%p_real",2,0)
     !> lecture de toutes les valeurs de parametres (avec traitement des erreurs)
     do i=0,nuser_param-1
         cur_node => item(userParam,i)
         call get_int_xml(cur_node,"Index",ind,1)
         if (ind>nuser_param .OR. ind<1) then
            write(err_msg,fmt="(A)") " indice negatif ou nul, ou superieur au nombre de parametres P_real <user> (read_param)"
            call amitex_abort(err_msg, 2,0)
         end if
         call get_real_xml(cur_node,"Value",user_param0%p_real(ind),1)
     end do
     end if ! test nuser_param

     !> --- READING STRING PARAMETERS
     cur_node => item(node_list0,0)
     userParam=> getElementsByTagName(cur_node,"P_string")
     nuser_param=getLength(userParam)
     if (nuser_param>0) then
     !> Allocation
     allocate(user_param0%p_string(nuser_param),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_param) pour alloue user_param%p_string",2,0)
     !> lecture de toutes les valeurs de parametres (avec traitement des erreurs)
     do i=0,nuser_param-1
         cur_node => item(userParam,i)
         call get_int_xml(cur_node,"Index",ind,0)
         if (ind>nuser_param .OR. ind<1) then
            write(err_msg,fmt="(A)") " indice negatif ou nul, ou superieur au nombre de parametres P_string <user> (read_param)"
            call amitex_abort(err_msg, 2,0)
         end if
         call get_str_xml_nolowercase(cur_node,"Value",user_param0%p_string(ind),0)
     end do
     end if ! test nuser_param

  elseif (getLength(node_list0) > 1) then
    call amitex_abort("Le noeud 'User' est present plus d'une fois (read_param)",2,0)
  end if

!!---------------------FIN DE LECTURE DES PARAMETRES
 call destroy(fi)

end subroutine read_param


!====================================================================
!====================================================================
!                         SUBROUTINE READ_PARAM_COMPOSITES
!
!>   Lecture des parametres des algorithmes composites dans un fichier xml
!!
!! Plusieurs parametres sont lus, à destination de l'algorithme d'intégration 
!! du modèle Multi-couche (Reuss plus tard égalemennt ??) :
!!   - La tolerance pour le critere d'equilibre umatLaminate
!!     (par defaut 1e-2 et doit etre superieure a 1e-1).
!!   - L'acceleration de convergence (vaut default, true ou false) pour umatLaminate
!!     La valeur par defaut est true.
!!   - Le type d'initialisation choisie pour l'algorithme umatLaminate
!!   - Le type de pilotage retenu pour l'algorithme umatLaminate
!!      * Nombre de sous pas temps par pas de temps : par défaut 1
!!
!! \param[in]   file_algo: nom du fichier xml (chaine de caracteres)
!!
!====================================================================
subroutine read_param_composite(file_algo)

  implicit none

  character(len=*), intent(in) :: file_algo
  type(node) ,pointer          :: fi, cur_node, cur_elem
  type(nodeList), pointer      :: node_list, elem_list
  character(len=50)            :: string, tmp_string
  !> Parametres par defaut
  real(mytype), parameter      :: default_tol_criteq = (10._mytype)**(-4)
  real(mytype), parameter      :: default_tol_criteq_voxcomp = (10._mytype)**(-2)
  real(mytype), parameter      :: max_tol_criteq = (10._mytype)**(-3)
  real(mytype), parameter      :: max_tol_criteq_voxcomp = (10._mytype)**(-1)
#ifdef DOUBLE_PREC
  real(mytype), parameter      :: min_tol_criteq = (10._mytype)**(-12)
#else
  real(mytype), parameter      :: min_tol_criteq = (10._mytype)**(-5)
#endif
  real(mytype), parameter      :: min_tol_criteq_voxcomp = (10._mytype)**(-15) !very low value (~ discard its usage) 
  real(mytype), parameter      :: default_perturbation = (10._mytype)**(-8)
  character(len=25)            :: Node_name,nomvc
  integer                      :: i

  fi => parseFile(trim(file_algo))

 !Boucle sur les 2 choix possibles (1 : laminate, 2 : reuss)
 do i = 1,2
 if (i == 1) then 
   Node_name = "Algorithm_laminate"
   nomvc="laminate"
 end if
 if (i == 2) then
   Node_name = "Algorithm_reuss"
   nomvc= "reuss"
 end if

! Paramètres pour les modèles composites Laminate / reuss
 node_list => getElementsByTagName(fi,trim(Node_name))

 if(getLength(node_list)>0)then
    ! Noeud présent dans le fichier xml algorithme
    cur_node => item(node_list,0)
    !! On recupere le critere de convergence demande
    string = "Convergence_Criterion"
    elem_list => getElementsByTagName(cur_node, string)
    if(getLength(elem_list)>0) then
       cur_elem => item(elem_list, 0)       
       string = "Value"
       if (i==1) call get_real_xml_default(cur_elem, string, algo_param0%tol_criteq_laminate,0,default_tol_criteq_voxcomp,1)
       if (i==2) call get_real_xml_default(cur_elem, string, algo_param0%tol_criteq_reuss,0,default_tol_criteq_voxcomp,1)
    else
       if (i==1) algo_param0%tol_criteq_laminate = default_tol_criteq_voxcomp
       if (i==2) algo_param0%tol_criteq_reuss = default_tol_criteq_voxcomp
    end if
    if (i==1) then
      if(algo_param0%tol_criteq_laminate > max_tol_criteq_voxcomp) then
        call amitex_abort("Critere d'arret de l'algorithme laminate trop eleve (read_param_composite)",2,0)
      elseif(algo_param0%tol_criteq_laminate < min_tol_criteq_voxcomp) then
        call amitex_abort("Critere d'arret de l'algorithme laminate trop faible (read_param_composite)",2,0)
      end if
    end if 
    if (i==2) then
      if(algo_param0%tol_criteq_reuss > max_tol_criteq_voxcomp) then
        call amitex_abort("Critere d'arret de l'algorithme reuss trop eleve (read_param_composite)",2,0)
      elseif(algo_param0%tol_criteq_reuss < min_tol_criteq_voxcomp) then
        call amitex_abort("Critere d'arret de l'algorithme reuss trop faible (read_param_composite)",2,0)
      end if
    end if 

    !! On verifie si l'utilisateur veut utiliser 
    !! l'acceleration de convergence
    string = "Convergence_Acceleration"
    elem_list => getElementsByTagName(cur_node, string)
    if(getLength(elem_list)>0)then
       cur_elem => item(elem_list, 0)       
       string = "Value"
       call get_str_xml(cur_elem, trim(string), tmp_string,1)
       if(tmp_string == "true") then
          if (i==1) algo_param0%acc_CV_laminate = .TRUE.
          if (i==2) algo_param0%acc_CV_reuss = .TRUE.
       elseif(tmp_string == "false") then
          if (i==1) algo_param0%acc_CV_laminate = .FALSE.           
          if (i==2) algo_param0%acc_CV_reuss= .FALSE.           
       else
          if (i==1) algo_param0%acc_CV_laminate = .true.  
          if (i==2) algo_param0%acc_CV_reuss = .true.  
          call amitex_abort("Acceleration de convergence pour le modele "//trim(nomvc)//" &
                         & mal renseignee :activee par défaut (read_param_composite)",-1,0)
       end if
    else 
          if (i==1) algo_param0%acc_CV_laminate = .true.  
          if (i==2) algo_param0%acc_CV_reuss = .true.  
          call amitex_abort("Acceleration de convergence pour le modele "//trim(nomvc)//" &
                         & mal renseignee :activee par défaut (read_param_composite)",-1,0)
    end if

    !! On recupere le type d'initialisation choisie par l'utilisateur
    string = "Initialisation_type"
    elem_list => getElementsByTagName(cur_node, string)
    if(getLength(elem_list)>0)then
       cur_elem => item(elem_list, 0)       
       string = "Value"
       call get_str_xml(cur_elem, trim(string), tmp_string,1)
       if((trim(tmp_string) == "proportional") .or. (trim(tmp_string) == "linear") .or. &
          (trim(tmp_string) == "default") ) then
          if (i==1) algo_param0%Init_laminate = trim(tmp_string)
          if (i==2) algo_param0%Init_reuss = trim(tmp_string)
       else
          if (i==1) algo_param0%Init_laminate = "default"
          if (i==2) algo_param0%Init_reuss = "default"
          call amitex_abort("Type d'initialisation "//trim(nomvc)//" mal renseignee : &
                            &valeur par defaut (read_param_composite)",-1,0)
       end if
    end if

    !! On récupère le nombre de sous pas de temps d'intégration de loi laminate souhaité par l'utilisateur
    string = "Nmax_subdivision"
    elem_list => getElementsByTagName(cur_node, string)
    if(getLength(elem_list)>0)then
       cur_elem => item(elem_list, 0)       
       string = "Value"
       if (i==1) call get_int_xml(cur_elem, string,algo_param0%Nmax_subdivision_laminate ,1)
       if (i==2) call get_int_xml(cur_elem, string,algo_param0%Nmax_subdivision_reuss ,1)
    else
       if (i==1) algo_param0%Nmax_subdivision_laminate = 5
       if (i==2) algo_param0%Nmax_subdivision_reuss = 5
       call amitex_abort("Nombre de subdivision maximale du pas de temps pour intégration du modèle &
            & "//trim(nomvc)//" : valeur par defaut : 5 <=> 32 sous incréments (read_param_composite)",-1,0)
    end if

    !! On recupere le type de matrice tangente choisie par l'utilisateur
    string = "Jacobian_matrix"
    elem_list => getElementsByTagName(cur_node, string)
    if(getLength(elem_list)>0)then
       cur_elem => item(elem_list, 0)       
       string = "Value"
       call get_str_xml(cur_elem, trim(string), tmp_string,1)
       if ( (trim(tmp_string) == "elastic") .or.  ( trim(tmp_string) == "mfront_elastic" ) .or. &
          ( trim(tmp_string) == "mfront_tangent" ) .or. ( trim(tmp_string) == "numerical" ) .or. &
           ( trim(tmp_string) == "mfront_secant") ) then
          if (i==1) algo_param0%Jacobian_type_laminate = trim(tmp_string)        
          if (i==2) algo_param0%Jacobian_type_reuss = trim(tmp_string)        
       else
          if (i==1) algo_param0%Jacobian_type_laminate = "numerical"
          if (i==2) algo_param0%Jacobian_type_reuss = "elastic"
          call amitex_abort("Type de matrice tangente "//trim(nomvc)//"  mal renseignée :&
                            & jacobienne calculée numériquement par defaut (read_param_composite)",-1,0)
       end if
    else
          if (i==1) algo_param0%Jacobian_type_laminate = "numerical"
          if (i==2) algo_param0%Jacobian_type_reuss = "elastic"
          call amitex_abort("Type de matrice tangente "//trim(nomvc)//" valeur par défaut :&
                            & jacobienne calculée numériquement (read_param_composite)",-1,0)       
    end if

    if (i==1 .and. trim(algo_param0%Jacobian_type_laminate) == "numerical") then
       !! On recupere la valeur de la perturbation pour le calcul numérique de la matrice tangente 
       !  du modèle multi-couches
       string = "Numerical_Jacobian_perturbation"
       elem_list => getElementsByTagName(cur_node, string)
       if(getLength(elem_list)>0)then
          cur_elem => item(elem_list, 0)       
          string = "Value"
          call get_real_xml_default(cur_elem, string, algo_param0%Perturbation_laminate,0,default_perturbation,1)
       else
          algo_param0%Perturbation_laminate = default_perturbation
          call amitex_abort("Valeur perturbation pour évaluation de la matrice tangente numerique du modele&
                           & Multi-couche : valeur par defaut -> 1e-6 (read_param_composite)",-1,0)
       end if
    end if
    if (i==2 .and. trim(algo_param0%Jacobian_type_reuss) == "numerical") then
       !! On recupere la valeur de la perturbation pour le calcul numérique de la matrice tangente 
       !  du modèle de reuss
       string = "Numerical_Jacobian_perturbation"
       elem_list => getElementsByTagName(cur_node, string)
       if(getLength(elem_list)>0)then
          cur_elem => item(elem_list, 0)       
          string = "Value"
          call get_real_xml_default(cur_elem, string, algo_param0%Perturbation_reuss,0,default_perturbation,1)
       else
          algo_param0%Perturbation_reuss = default_perturbation
          call amitex_abort("Valeur perturbation pour évaluation de la matrice tangente numerique du modele&
                           & de Reuss : valeur par defaut -> 1e-6 (read_param_composite)",-1,0)
       end if
    end if


 else
    ! Noeud non présent dans le fichier xml algorithme
    if (i==1) then
      call amitex_abort("Pas de paramètres pour l'algorithmes laminate :&
                        & parametres algorithmiques par défaut (read_param_composite)",-1,0)   
      algo_param0%tol_criteq_laminate = default_tol_criteq_voxcomp
      algo_param0%acc_CV_laminate = .true.
      algo_param0%Init_laminate = "default"
      algo_param0%Nmax_subdivision_laminate  = 5
      algo_param0%Jacobian_type_laminate = "numerical"
      algo_param0%Perturbation_laminate = default_perturbation
    end if
    if (i==2) then
      call amitex_abort("Pas de paramètres pour l'algorithmes reuss :&
                        & parametres algorithmiques par défaut (read_param_composite)",-1,0)   

      algo_param0%tol_criteq_reuss = default_tol_criteq_voxcomp
      algo_param0%acc_CV_reuss = .true.
      algo_param0%Init_reuss = "default"
      algo_param0%Nmax_subdivision_reuss  = 0
      algo_param0%Jacobian_type_reuss = "elastic"
      algo_param0%Perturbation_reuss = default_perturbation
    end if 
 end if

 end do ! fin boucle sur les deux choix possibles

 call destroy(fi)

end subroutine read_param_composite

end module param_algo_mod
