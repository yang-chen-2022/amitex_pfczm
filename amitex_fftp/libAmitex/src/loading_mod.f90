!123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!===================================================================================================
!
!    MODULE LOADING_MOD
!
!>  Definition, initialisation et desallocation du chargement et des parametres d'extraction
!!
!!  Structures :
!! - LOADING       : structure du chargement
!! - PARAM_EXTRACT : structure des parametres d'extraction
!! - INITLOADEXT   : structure d'initialisation des parametres exterieurs (dont temperature)
!!
!! Subroutines :
!! - read_load :     lecture du chargement et des parametres d'extraction a partir du
!!               fichier xml
!! - get_loadList_xml : lecture d'un pas de chargement
!! - getLoadTime :  extrait du fichier xml les informations des temps de chargement
!! - getLoadParam :  extrait du fichier xml, pour un pas de chargement,les informations 
!!                   relatives a la temperature et aux parametres externes
!! - get_load :      extrait du fichier xml, pour un pas de chargement et une
!!               composante les informations relatives au  pilotage
!! - repartition_discretization : associe un entier au type de discretisation du temps
!!               souhaite
!! - repartition_evolution : associe un entier au type d'evolution du chargement
!!               souhaite 
!! - set_time:              calcule le temps de chaque pas de chargement
!! - get_current_nb_param:  renvoie le nombre de parametres externes du chargement courant
!! - get_current_loading:    calcule le chargement courant
!! - get_loading_value:     calcul d'une composante du chargement courant (a
!!                 l'execution)
!! - deallocate_load: desallocation du chargement
!! - deallocate_param_extract: desallocation des parametres d'extraction
!!
!!
!!
!===================================================================================================

module loading_mod

#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi,       only : nrank
#else
  use decomp_2d, only : nrank,mytype 
#endif  
  use MPI
  use Fox_dom
  use io_amitex_mod
  use material_mod
  use error_mod
  use param_algo_mod
  use amitex_mpi_mod, only: amitex_mpi_comm

  private

  !> Pointer "publiques" vers composantes de SIMU_AMITEX
  public ::  load, initValExt,extract,&
             local_load, local_loadD,&
             n_gradgradU, add_def_star 

  !> Variables "publiques" utilisees lors de l'initialisation
  public ::  load0, initValExt0, extract0,&
             local_load0, local_loadD0, &
             n_gradgradU0, add_def_star0

  !> Types publiques (pour definition de SIMU_AMITEX)
  public ::  LOADING, INITLOADEXT,LOCALLOAD, PARAM_EXTRACT

  !> Fonctions publiques
  public ::  read_load,print_loading,print_initLoadExt,print_param_extract,&
             get_current_loading,get_current_loadingD,get_dt_from_load,&
             deallocate_load, deallocate_param_extract0,deallocate_initValExt,deallocate_local_load
 
!===================================================================================================
!> Chargement a un instant : tableau de reel defini a partir du tableau de type LOADING 
!!                           (voir fonction get_current_loading)
!!

  type LOCALLOAD
     real(mytype),allocatable,dimension(:) :: t0,t1,dt
                                                        !< MECANIQUE
                                                        !< chargement impose a l'instant t 
                                                        !< type pilotage, valeurs imposees, temperature, param. ext
                                                        !< tableau (Ntens+Ntens+1+nb. param.ext.+ (27 gradgradU))
                                                        !<
                                                        !< DIFFUSION
                                                        !< chargement impose a l'instant t 
                                                        !< type pilotage, valeurs imposees, temperature, param. ext
                                                        !< tableau (3nVarD+3nVarD+1+nb. param.ext.)
                                                        
     real(mytype), dimension(-2:0) :: t_load =0         !<temps de chargement aux instant t-2, t-1, t
  end type LOCALLOAD

!! MECA
!! local_load%t1(1:algo_param%nTensDef) :                        pilotage 
!! local_load%t1(algo_param%nTensDef+1:2*algo_param%nTensDef) :  valeurs associées au pilotage
!! local_load%t1(2*algo_param%nTensDef+1) :                      temperature
!! local_load%t1(nb_param next indices) :                        parametres externes (s'ils existent)
!! local_load%t1(27 next indices) :                              gradgradU components (if exists)
!
!! DIFFUSION
!! local_loadD%t1(1:3*algo_param%nVarD) :                        pilotage 
!! local_loadD%t1(3*algo_param%nVarD+1:6*algo_param%nVarD) :     valeurs associées au pilotage
!! local_loadD%t1(6*algo_param%nVarD+1) :                        temperature
!! local_loadD%t1(6*algo_param%nVarD+2:) :                       parametres externes (s'ils existent)


!===================================================================================================
!> Structure d'initialisation des parametres externes
!!
!! is_temp_initialized : booleen permettant de savoir si la temperature est initialisee 
!!                       (vaut .FALSE. par defaut)
!! nb_param : nombre de parametres externes initialises (vaut 0 par defaut)
!!
!! Si la temperature est initialisee, il faut qu'elle soit precisee a chaque chargement
!! Si la temperature est precisee a un chargement il faut qu'elle ait ete initialisee
!! De meme si nb_param>0, il faut que chaque parametre soit precise a chaque chargement 
!! Si un parametre est precise alors que nb_param=0 => erreur
  type INITLOADEXT
     logical                               :: is_temp_initialized
     integer                               :: nb_param     
     real(mytype)                          :: temp
     real(mytype),allocatable,dimension(:) :: param_values
  end type INITLOADEXT
 
!===================================================================================================
!> Structure de chargement
!!
!!
!! Notation de Voigt (en HPP), ordre (11, 22, 33, 12, 13, 23)
!!
!! NIncr : nombre d'increment du chargement
!!
!! time  : (tableau) temps a chaque increment du chargement
!!
!! discretization : discretisation du temps
!!               - 0 "User"      discretisation definie par l'utilisateur (liste dans le fichier xml)
!!               - 1 "Linear"    discretization lineaire
!!                          -# borne inferieure : 0 pour le premier chargement
!!                                                temps a l'increment precedent sinon
!!                          -# borne superieure : "Tfinal" defini dans le fichier xml
!!
!! driving:
!!               - 0 pilotage en deformation
!!               - 1 pilotage en contrainte
!!
!! value:        - valeur final imposee (en contrainte ou deformation)
!!
!! evolution:
!!               - "Constant" conserve la contrainte/deformation moyenne a l'instant precedent pour les increments
!!               - "Linear"   evolution lineaire de la valeur moyenne a l'instant precedent jusqu'a la valeur: "value"
!!
!------------------------------------------------------------------------------
  type LOADING


    !> Nombre d'increment de temps
    integer                                     :: NIncr 
    !> Temps a chaque increment, d'indice 0:Nincr avec time(0)=temps du demi-pas initial
    !> ATTENTION : ce vecteur est modifie en cas d'interruption detectee : 
    !>     chargement courant : (time(i>i_interruption) = t_interruption)
    !>     chargements suivants : decalage du temps induit par l'interruption du chargement courant
    real(mytype),allocatable,dimension(:)       :: time  
    !> Type de discretisation du temps, 0=user, 1=linear
    integer                                     :: Discretization 
    !> Type de pilotage 0: deformation,  1: contrainte (vecteur nTens)
    integer, allocatable,dimension(:)           :: driving   
    !> Type de pilotage 0: deformation,  1: contrainte (vecteur 3*nVarD, nVarD=1 pour commencer!)
    integer, allocatable,dimension(:)           :: drivingD   
    !> Valeurs imposees en fin de chargement mecanique, par composante (vecteur nTens)
    real(mytype), allocatable,dimension(:)      :: mvalue
    !> Valeurs imposees du gradient de gardient de U, en fin de chargement mecanique, tableau (3,3,3)
    real(mytype), dimension(3,3,3)              :: gvalue=0.
    !> Valeurs imposees en fin de chargement diffusion, par composante (vecteur 3*nVarD, nVarD=1 pour commencer!)
    real(mytype), allocatable,dimension(:)      :: dvalue
    !> Valeurs imposees en fin de chargement pour la temperature et les parametres externes
    real(mytype)                                :: tvalue
    real(mytype),allocatable,dimension(:)       :: pvalue
    !> Type d'evolution du chargement mécanique, par composante (vecteur nTens), 0 : constant, 1=linear
    integer, allocatable,dimension(:)           :: mevolution 
    !> Type d'evolution du chargement mécanique, par composante (3,3,3), 0 = non accepte, 1=linear, -1=non pilotee
    integer, dimension(3,3,3)                   :: gevolution=-1 
    !> Type d'evolution du chargement diffusion, par composante (vecteur 3*nVarD), 0 : constant, 1=linear
    integer, allocatable,dimension(:)           :: devolution 
    !> Type d'evolution de la temperature et des parametres externes, 0 : constant, 1=linear
    integer                                     :: tevolution 
    integer,allocatable,dimension(:)            :: pevolution
    !> Chargement a direction de contrainte imposee, DirStress (vecteur Ntens)
    logical                                     :: DirStress_flag = .false.
    logical                                     :: DirStress_cauchy 
    real(mytype), allocatable,dimension(:)      :: DirStress
    !> DirStress2 egal a DirStress en general.
    !> Si DirStress_cauchy : conversion Cauchy -> PK1
    !>                       reevalue a chaque pas, utilise pour la convergence Macro 
    real(mytype), allocatable,dimension(:)      :: DirStress2
    !> Introduction de parametres pour une interruption 'User'
    logical                                     :: User_interruption_flag = .false.
    real(mytype), allocatable,dimension(:)      :: User_interruption
    !> Restart 
    logical                                     :: restart_save        = .false.  ! .true. if saving files for restart is required in the loading
    logical                                     :: skip_before_restart = .false.  ! .true. if the loading if before a restart starting point
    integer                                     :: restart_inc         = -1       ! starting increment (read in saved restart files)
                                                                                  ! default = -1 (no restarting increment)
    integer                                     :: restart_save_period = huge(1)  ! save restart files every restart_save_period  increments

  end type LOADING
!------------------------------------------------------------------------------

!===================================================================================================
!> Tableau d'entiers codes sur 1 octet
  type TAB_LOGICAL
    logical,allocatable,dimension(:) :: val
  end type TAB_LOGICAL
!------------------------------------------------------------------------------

!===================================================================================================
!> Structure des parametres d'extraction des donnees
!!
!! \todo Changer TABLEAU_INT1 en TABLEAU_LOGICAL
!------------------------------------------------------------------------------
  type PARAM_EXTRACT
    !> Deformations a extraire dans le fichier VTK 
    logical                                     :: defVTK=.false.
    !> Contraintes a extraire dans le fichier VTK 
    logical                                     :: sigVTK=.false.
    !> Flux a extraire dans le fichier VTK 
    logical                                     :: FluxDVTK=.false.
    !> Flux a extraire dans le fichier VTK 
    logical                                     :: GradDVTK=.false.
    !> Variables internes a extraire dans le fichier VTK, pour chaque materiau \n  
    !! dimension : nombre de materiaux  
    !! varVTK%val : tableau de logique precisant si on sort la VI (dimension nombre de VI)
    type(TAB_LOGICAL),allocatable,dimension(:)  :: varVTK
    !> Test indiquant les temps auxquels on extrait les donnees VTK \n
    !! dimension: nombre de pas de temps 
    logical,allocatable,dimension(:)            :: tpsVTK   

    !> Test indiquant, pour chaque materiau, de sortir des valeurs moyennes par zone \n
    !! dimension: nombre de materiaux    
    logical, allocatable, dimension(:)          :: printMeanZ
    !> Variables internes dont on souhaite une sortie par zone, pour chaque materiau \n
    !! dimension: nombre de matériaux 
    !! varZSTD%val : tableau de logique precisant si on sort la VI (dimension nombre de VI)
    type(TAB_LOGICAL),allocatable,dimension(:)  :: varZSTD
    !> Test indiquant les temps auxquels on extrait les donnees par zone dans le .zstd \n
    !! dimension: nombre de pas de temps
    logical,allocatable,dimension(:)            :: tpsZSTD     

    
    !> Nombre de pas de temps qu'il faut sauvegarder par zone dans chaque chargement \n
    !! dimension : nombre de chargements élémentaires
    integer, allocatable, dimension(:)          :: numberZstd   
    !> Variable pour savoir si au moins une valeur par zone doit etre sortie
    logical                                     :: printZSTD

    !> Tests indiquant les temps auxquels on extrait les donnees par cellule et par materiau (.std et .mstd)\n
    !! dimension: nombre de pas de temps
    logical,allocatable,dimension(:)            :: tpsSTD     

    !> Nombre de pas de temps qu'il faut sauvegarder (/cellule et /materiau) dans chaque chargement \n
    !! dimension : nombre de chargements élémentaires
    integer, allocatable, dimension(:)          :: numberStd   


  end type PARAM_EXTRACT
!------------------------------------------------------------------------------

!> Le chargement
  type(LOADING), allocatable,dimension(:) :: load0
  type(LOADING), pointer,dimension(:)     :: load
  type(INITLOADEXT)                       :: initValExt0
  type(INITLOADEXT),pointer               :: initValExt
  type(LOCALLOAD)                         :: local_load0
  type(LOCALLOAD),pointer                 :: local_load
  type(LOCALLOAD)                         :: local_loadD0
  type(LOCALLOAD),pointer                 :: local_loadD

!> Additionnal gradgradU loading
  integer                                 :: n_gradgradU0 = 0 !< number of components to impose for gradgradU (0 or 27) 
  integer,pointer                         :: n_gradgradU

!> Additionnal "star" (compatible) strain loading
  logical                                 :: add_def_star0 = .false.
  logical, pointer                        :: add_def_star
  
!> Les parametres d'extraction
  type(PARAM_EXTRACT)                     :: extract0
  type(PARAM_EXTRACT),pointer             :: extract
  
!> Ratio pour le calcul du pas de temps 0 en debut de chargement
  real(mytype), parameter :: ratio=0.4_mytype
  

contains


!===================================================================================================
!===================================================================================================
!                         SUBROUTINE CHECK_LOADING
!
!> Verifie la compatibilite de donnees de la variable de type loading
!!
!! TODO : developper les verifications
!!
!!  \param[in] Flog:   (entier) unite logique du fichier de sortie .log
!!  \param[in] nrank0: (entier) pinceau pour lequel on effectue la sortie
!!
!===================================================================================================
subroutine check_loading(Flog,nrank0)

  implicit none
  integer, intent(in)   :: Flog
  integer, intent(in)   :: nrank0
  integer               :: i,k,l
  logical               :: test
  character(len=200)    :: err_msg


  if (nrank==nrank0) then
  do i = 1,size(load0) !boucle sur les chargements
  
  !! SYMETRIE de DirStress si DirStress_cauchy
  if (load0(i)%DirStress_cauchy .AND. (.not. algo_param0%HPP) .and. algo_param0%Mechanics) then
     test = (load0(i)%DirStress(4)==load0(i)%DirStress(7)) &
       .and.(load0(i)%DirStress(5)==load0(i)%DirStress(8)) &
       .and.(load0(i)%DirStress(6)==load0(i)%DirStress(9)) 
     if (.not. test) then
        write(err_msg,fmt="(A,I0,A)") &
           "tenseur 'DirStress', noeud 'Loading' ",i,&
           " non symetrique (incompatible avec <DirStress Type=""cauchy"" /> "
        call amitex_abort(err_msg,1,0)
        write(Flog,"(A)") err_msg
     end if
  end if

  !! PLUS d'une composante de deformation impose
  if (load0(i)%DirStress_flag .and. algo_param0%Mechanics) then
     l=0
     do k=1,algo_param0%Ntensdef
        if (load0(i)%driving(k) == 0) l = l+1
     end do   
     if (l==0) then
        write(err_msg,fmt="(A,I0,A)") &
        "pas de composante en deformation imposee, noeud 'Loading' ",i,&
        " (incompatible avec ""DirStress"" impose) "
        call amitex_abort(err_msg,1,0)
        write(Flog,"(A)") err_msg
     elseif (l .gt. 1) then
        write(err_msg,fmt="(A,I0,A)") &
           "plus d'une composante en deformation imposee, noeud 'Loading' ",i,&
           " (pas forcement compatible avec ""DirStress"" impose) "
        call amitex_abort(err_msg,0,0)
        write(Flog,"(A)") err_msg
     end if
  end if
  
  !! RESTART
  if (load0(i)%restart_save .and. (load0(i)%restart_save_period==huge(1).or. load0(i)%restart_save_period<=0)) then
        write(err_msg,fmt="(A)") &
           "Restart is activated but the restart_save_period is huge or <=0  "
        call amitex_abort(err_msg,1,0)
        write(Flog,"(A)") err_msg
  end if  
  end do ! fin boucle sur les chargements
  end if ! fin si nrank==0

end subroutine check_loading


!===================================================================================================
!===================================================================================================
!                         SUBROUTINE PRINT_LOADING
!
!> Ecriture des composantes d'une variable de type loading
!!
!!
!!  \param[in] Flog:   (entier) unite logique du fichier de sortie .log
!!  \param[in] nrank0: (entier) pinceau pour lequel on effectue la sortie
!!
!===================================================================================================
subroutine print_loading(Flog,nrank0)
  
  implicit none
  integer, intent(in)   :: Flog
  integer, intent(in)   :: nrank0
  integer               :: i,k,l

  call check_loading(Flog,nrank0)

  if (nrank==nrank0) then
     write(Flog,"(A)") " "
     write(Flog,"(A,I0)") "STRUCTURE LOADING, pinceau ", nrank0
     do i=1,size(load0)
        write(Flog,"(A,I0)") "----Indice ", i 
        write(Flog,"(A,I0)") "Nincr ", load0(i)%Nincr
        write(Flog,"(A,E15.8)") "Temps final ", load0(i)%time(load0(i)%Nincr)
        write(Flog,"(A,I0)") "Discretisation (0=user, 1=lin.) ", load0(i)%Discretization

        if (algo_param0%Mechanics) then
        do k=1,algo_param0%nTensDef
           write(Flog,"(A,I0)") "Pilotage (0=Strain, 1=Stress) ", load0(i)%driving(k)
           write(Flog,"(A,I0)") "Evolution, Meca, (0=const., 1=lin.) ", load0(i)%mevolution(k)
           if (load0(i)%mevolution(k).eq.1) &
              write(Flog,"(A,E15.8)") "Valeur imposee, Meca, fin de charg. ", load0(i)%mvalue(k)
        end do
        if (load0(i)%dirStress_flag) then
           do k=1,algo_param0%nTensDef
              write(Flog,"(A,I0,A,E15.8)") "Direction ",k," de contrainte fixe ", load0(i)%dirStress(k)
           end do
           if     ((.not. algo_param0%HPP) .and. load0(i)%DirStress_cauchy) then
              write(Flog,"(A)") "Type de contrainte pour la direction de contrainte imposee : Cauchy"
           elseif ((.not. algo_param0%HPP) .and. (.not. load0(i)%DirStress_cauchy)) then
              write(Flog,"(A)") "Type de contrainte pour la direction de contrainte imposee : PK1"
           end if
        end if
        if (load0(i)%User_interruption_flag) then
           do k=1,size(load0(i)%User_interruption)
              write(Flog,"(A,I0,A,E15.8)") "User_interruption ",k," : ", load0(i)%User_interruption(k)
           end do
        end if
        
        end if

        if (algo_param0%Diffusion) then
        do k=1,algo_param0%nVarD
        do l=1,3
          write(Flog,"(A,I0)") "Pilotage (0=GradD, 1=FluxD) ", load0(i)%drivingD((k-1)+l)
          write(Flog,"(A,I0)") "Evolution, Diff, (0=const., 1=lin.) ", load0(i)%devolution((k-1)+l)
          if (load0(i)%devolution((k-1)+l).eq.1) &
             write(Flog,"(A,E15.8)") "Valeur imposee, Diffusion, fin de charg. ", load0(i)%dvalue((k-1)+l)
        end do
        end do
        end if

        write(Flog,"(A,I0)") "Evolution, Temp., (0=const., 1=lin.) ", load0(i)%tevolution
        write(Flog,"(A,E15.8)") "Temperature, fin de charg. ", load0(i)%tvalue
        if (initValExt0%nb_param>0) then   
        do k=1,initValExt0%nb_param  
            write(Flog,"(A,I0)") "Evolution, par. ext., (0=const., 1=lin.) ", load0(i)%pevolution(k)
            if (load0(i)%pevolution(k).eq.1) &
               write(Flog,"(A,E15.8)") "Parametres ext., fin de charg. ", load0(i)%pvalue(k)
        end do
        end if  

        if (load0(i)%restart_save) then
          write(Flog,"(A)") "RESTART INFOS"
          write(Flog,"(A,L)") "Save restart activated : ", load0(i)%restart_save
          write(Flog,"(A,I0,A)") "Save restart-files every : ", load0(i)%restart_save_period, " increments"
          write(Flog,"(A,L)") "Skip partial loading (before restart point): ", load0(i)%skip_before_restart      
          write(Flog,"(A,I0)") "Restarting increment (restart ppoint): ", load0(i)%restart_inc         
        end if

     end do
     write(Flog,"(A)") " "
  end if

end subroutine print_loading

!===================================================================================================
!                         SUBROUTINE PRINT_INITLOADEXT
!
!> Ecriture des composantes de d'une variable de type loading
!!
!!
!!  \param[in] Flog:   (entier) unite logique du fichier de sortie .log
!!  \param[in] nrank0: (entier) pinceau pour lequel on effectue la sortie
!!
!===================================================================================================
subroutine print_initloadext(Flog,nrank0)
  
  implicit none
  integer, intent(in)   :: Flog
  integer, intent(in)   :: nrank0
  integer               :: k

  if (nrank==nrank0) then
     write(Flog,"(A)") " "
     write(Flog,"(A,I0)") "STRUCTURE INITLOADEXT, pinceau ", nrank0
     if (initValExt0%is_temp_initialized) then
        write(Flog,"(A)") "Temperature initialisee"
     else
        write(Flog,"(A)") "Temperature non initialisee par l'utilisateur "
     end if
     write(Flog,"(A,I0)") "Nombre de parametres exterieurs ", initValExt0%nb_param
     write(Flog,"(A,E15.8)") "Temperature initiale ", initValExt0%temp
     if (initValExt0%nb_param>0) then   
     do k=1,initValExt0%nb_param  
        write(Flog,"(A,E15.8)") "Parametre exterieur initial ", initValExt0%param_values(k)
     end do
     end if  
     write(Flog,"(A)") " "
  end if

end subroutine print_initloadext

!===================================================================================================
!                         SUBROUTINE PRINT_PARAM_EXTRACT
!
!> Ecriture des composantes d'une variable de type param_extract
!!
!!
!!  \param[in] Flog:   (entier) unite logique du fichier de sortie .log
!!  \param[in] nrank0: (entier) pinceau pour lequel on effectue la sortie
!!  \param[in] nmate: (entier) nombre de materiaux
!!
!===================================================================================================
subroutine print_param_extract(Flog,nrank0,nmate)
  
  implicit none
  integer, intent(in)   :: Flog
  integer, intent(in)   :: nrank0,nmate
  integer               :: k, sumvtk,sumzstd
  
  sumvtk = count(extract0%tpsVTK)
  sumzstd = count(extract0%tpsZSTD)

  if (nrank==nrank0) then
     write(Flog,"(A)") " "
     write(Flog,"(A,I0)") "STRUCTURE PARAM_EXTRACT, pinceau ", nrank0
     write(Flog,"(A,I0)") "Nombre de pas de temps ou extraire des fichiers vtk  ", sumvtk
     write(Flog,"(A,I0)") "Nombre de pas de temps ou extraire des fichiers zstd ", sumzstd
     if (algo_param0%Mechanics) then
         write(Flog,"(A,L2)") "Sortie vtk des contraintes ", extract0%sigVTK
         write(Flog,"(A,L2)") "Sortie vtk des deformation ", extract0%defVTK
         do k = 1,size(MattotP0)
            write(Flog,"(A,I0)") "Nombre de variables internes a sortir en vtk ", size(extract0%varVTK(k)%val)
            write(Flog,"(A,I0)") "Nombre de variables internes a sortir par zone ", size(extract0%varZSTD(k)%val)
         end do
     end if 
     if (algo_param0%Diffusion) then
         write(Flog,"(A,L2)") "Sortie vtk des FluxD ", extract0%FluxDVTK
         write(Flog,"(A,L2)") "Sortie vtk des GradD ", extract0%GradDVTK
     end if 

     do k = 1,nmate
       write(Flog,"(A,I0,A,L2)") "Sortie des moyennes par zone pour le materiau ", k," : ", extract0%printMeanZ(k)
     end do

     do k = 1,size(load0)
         write(Flog,"(A,I0,A,I0)") "Nombre de sortie .zstd pour le chargement ", k," : ", extract0%numberZstd(k)
     end do
     write(Flog,"(A)") " "

     do k = 1,size(load0)
         write(Flog,"(A,I0,A,I0)") "Nombre de sortie .std pour le chargement ", k," : ", extract0%numberStd(k)
     end do
     write(Flog,"(A)") " "
  end if

end subroutine print_param_extract


!===========================================================================================================
!> SUBROUTINE get_dt_from_load(dt,dt_,load_n, load_incr)
!!
!! \param[out]     dt (real)            delta(t) between increments load_incr and load_incr-1
!! \param[out]     dt_ (real)           delta(t) between increments load_incr-1 and load_incr-2
!! \param[in]      load_n (integer)     index of the loading
!! \param[in]      load_incr (integer)  increment in the loading load_n   
!! 
!! Remark : special cases for load_incr=0 or 1
!===========================================================================================================
subroutine get_dt_from_load(dt,dt_,load_n, load_incr)

 implicit none
 
 real(mytype), intent(out)         :: dt, dt_
 integer, intent(in)               :: load_n, load_incr

 integer,dimension(1)              :: ub
 
 
 if (load_n > 1) ub = ubound(load(load_n-1)%time)

 if (load_incr==0 .and. load_n==1) then
    dt = load(load_n)%time(load_incr)
    dt_= 0.
 elseif (load_incr==1 .and. load_n==1) then
    dt = load(load_n)%time(load_incr) - load(load_n)%time(load_incr-1)
    dt_= load(load_n)%time(load_incr-1)
 elseif (load_incr==0 .and. load_n > 1) then
    dt  = load(load_n)%time(load_incr)      - load(load_n-1)%time(ub(1))
    dt_ = load(load_n-1)%time(ub(1)) - load(load_n-1)%time(ub(1)-1)
 elseif (load_incr==1 .and. load_n > 1) then
    dt  = load(load_n)%time(load_incr)    -  load(load_n)%time(load_incr-1)
    dt_ = load(load_n)%time(load_incr-1)  -  load(load_n-1)%time(ub(1)) 
 else
    dt  = load(load_n)%time(load_incr)   -  load(load_n)%time(load_incr-1)
    dt_ = load(load_n)%time(load_incr-1) -  load(load_n)%time(load_incr-2)
 end if
 
end subroutine get_dt_from_load



!===================================================================================================
!
!                     SUBROUTINE READ_LOAD
!
!> Lecture du chargement
!!
!! A partir d'un fichier xml, lit pour chaque chargement
!!     - le type de pilotage, le type d'evolution et la valeur finale imposee.
!!     - le temps en fin du chargement, le nombre d'increment de temps et
!!       le type de discretisation du temps.
!!     - les donnees a extraires (contraintes, deformations, variables internes)
!!              et les increments de temps auxquels extraire ces donnees.
!! Les resultats sont stockes dans chargement (type LOADING) 
!! et extract (type PARAM_EXTRACT)
!!
!! \param[in]   file_char: (chaine de caracteres) nom du fichier xml contenant
!!                  les informations de chargement
!! \param[in]   nmate: (entier) nombre de materiaux
!!
!! \param[out]   nb_tps (entier) nombre de pas de chargement
!!
!!
!===================================================================================================
subroutine read_load(file_char,nmate,nb_tps)

  implicit none
  
  character(len=*), intent(in)  :: file_char
  integer,intent(in)            :: nmate
  integer, intent(out)          :: nb_tps
 
  type(node), pointer        :: fi, cur_node
  type(nodeList), pointer    :: load_list, sub_list
  
  integer :: i, j, k, n, nb_load, err 
  integer,allocatable,dimension(:) :: tmp_array
  integer :: tag
  character(len=300)    :: err_msg
  character(len=100000) :: vtk_list
  logical :: ok_i
  
  integer :: alloc_stat

  allocate(tmp_array(1));deallocate(tmp_array) !avoid gcc-warning
  
  err=0
  fi =>parsefile(trim(file_char))

  !!!--------------------------
  !! donnees vtk a extraire
  !! composantes %defVTK,sigVTK,...  .false. par defaut, mises a .true. si renseignee dans Output
  !! si le pointeur du noeud "output" n'est pas associé => pas de sortie vtk 
  !! en pratique ce mecanisme ne marche pas (plus? a-t-il marché?) : load_list sort de GetNodeList à l'état associated=.true.
  call GetNodeList(fi, load_list, "Output", 1,1)
  !!call GetNodeList(fi, load_list, "Output", -1,1) !! tentative pour corriger
  !!if(nrank==0)print*,"associated(load_list) = ", associated(load_list)
  if(associated(load_list) .and. (algo_param0%Mechanics .eqv. .true.))then 
    cur_node => item(load_list,0)
    call GetNodeList(cur_node , sub_list, "vtk_StressStrain", 1,1)
    cur_node => item(sub_list,0)
    call get_int_xml(cur_node,"Strain",i,1,0)
    if (i==1) extract0%defVTK=.true.
    call get_int_xml(cur_node,"Stress",i,1,0)
    if (i==1) extract0%sigVTK=.true.
  end if
  if(associated(load_list) .and. algo_param0%Diffusion)then 
    cur_node => item(load_list,0)
    call GetNodeList(cur_node , sub_list, "vtk_FluxDGradD", 1,1)
    cur_node => item(sub_list,0)
    call get_int_xml(cur_node,"FluxD",i,1,0)
    if (i==1) extract0%FluxDVTK=.true.
    call get_int_xml(cur_node,"GradD",i,1,0)
    if (i==1) extract0%GradDVTK=.true.
  end if

    
  allocate(tmp_array(nmate),stat=alloc_stat)       ! nombre de variable internes par materiau (parallelise)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_load1)",2)
  tmp_array = 0 
  if (algo_param0%Mechanics) then
    allocate(extract0%VarVTK(nmate),stat=alloc_stat)  ! nombre de variable internes par materiau (global)
    allocate(extract0%varZSTD(nmate),stat=alloc_stat)  ! nombre de variable internes par materiau (global)
  end if
  ! booleen pour savoir si on imprime les moyennes par zone dans chaque materiau
  allocate(extract0%printMeanZ(nmate),stat=alloc_stat) 
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_load2)",2)
  extract0%printMeanZ=.false.
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_load3)",2)

  if (algo_param0%Mechanics) then
  do i=1,size(MattotP0)
    !! nombre de variables internes pour chaque materiau,
    !! permet la correspondance : numero de materiau global <-> parallelise
    tmp_array(MattotP0(i)%numM) = MattotP0(i)%NVarInt
    ! les materiaux non presents dans le pinceau ont "0" variables internes
  end do
  do i=1,nmate
    ! remplace les "0" dus a l'abscence de materiau par le nombre de variables
    ! internes
    call MPI_Allreduce(tmp_array(i),n,1,MPI_INTEGER,MPI_MAX,amitex_mpi_comm,err)
    ! allocation des tableaux extract0%varVTK
    allocate(extract0%varVTK(i)%val(n),stat=alloc_stat)
    allocate(extract0%varZSTD(i)%val(n),stat=alloc_stat)
    if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_load4)",2)
    extract0%varVTK(i)%val= .false. 
    extract0%varZSTD(i)%val= .false. 
  end do
  end if
  deallocate(tmp_array)

  !! liste des noeuds "vtk_IntVar" : liste des variables internes a extraire
  !! si le pointeur du noeud "output" n'est pas associé => pas de sortie spécifique (zstd) 
  if(associated(load_list))then 
    if (algo_param0%Mechanics) then
    sub_list => getElementsByTagName(item(load_list,0),"vtk_IntVarList")  
    !! Mise a jour des variables dans le chargement pour savoir quelles variables internes
    !! doivent etre ecrites
    if(getLength(sub_list)>0) call set_intvar_to_print(sub_list,nmate)
    end if
    !! Liste des noeuds "zone" : liste des materiaux ou des sorties par zone seront ecrites
    sub_list => getElementsByTagName(item(load_list,0),"Zone")  
    !! Mise a jour des variables dans le chargement pour savoir quelles variables internes
    !! doivent etre ecrites
    if(getLength(sub_list)>0) then
       extract0%printZSTD = .TRUE.
       call set_zone_to_print(sub_list,nmate)
    else
       extract0%printZSTD = .FALSE.
    end if
  else
     extract0%printZSTD = .FALSE.
  end if
  !! Initialisation des parametres externes et de la temperature
  call getInitLoadExt(fi)

  !------------------------------------------------------------------------------------------
  !!                                                                  LECTURE DES CHARGEMENTS 

  !! Liste des noeuds de chargement
  load_list => getElementsByTagName(fi,"Loading")
  nb_load = getlength(load_list)
  !! Allocation de l'ensemble des chargements
  allocate(load0(nb_load),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_load5)",2)
  load0(:)%NIncr=-1
  if(getLength(load_list)==0)then
     write(err_msg,fmt="(A)") "Aucun chargement (read_load6)"
     call amitex_abort(err_msg,1,0)
  end if
  !! Pour tous les chargements
  do i=0,getLength(load_list)-1
     ok_i=.TRUE.
     tag=-1
     !! On definit le numero du chargement
     call get_int_xml(item(load_list,i),"Tag",tag,1)
     if(tag==-1)then
        ok_i=.FALSE.
     end if
     !! verification de la validite du numero
     if(tag>getLength(load_list) .or. tag<1 )then
        write(err_msg,fmt="(A,I0,A)") "Numero de chargement (Tag) invalide", tag," (read_load7)"
        call amitex_abort(err_msg,1,0)
        ok_i=.FALSE.
     else if(load0(tag)%NIncr>-1.)then
        write(err_msg,fmt="(A,I0,A)") "Numero de chargement (Tag) deja traite", tag, " (read_load8)"
        call amitex_abort(err_msg,1,0)
        ok_i=.FALSE.
     else
        !! Allocation des tableaux pour ce chargement
        if(algo_param0%Mechanics) then
        allocate(load0(tag)%driving(algo_param0%nTensDef),stat=alloc_stat)
        allocate(load0(tag)%mvalue(algo_param0%nTensDef),stat=alloc_stat)
        allocate(load0(tag)%mevolution(algo_param0%nTensDef),stat=alloc_stat)
        allocate(load0(tag)%DirStress(algo_param0%nTensDef),stat=alloc_stat)
        allocate(load0(tag)%DirStress2(algo_param0%nTensDef),stat=alloc_stat)
        end if
        if(algo_param0%Diffusion) then
        allocate(load0(tag)%drivingD(3*algo_param0%nVarD),stat=alloc_stat)
        allocate(load0(tag)%dvalue(3*algo_param0%nVarD),stat=alloc_stat)
        allocate(load0(tag)%devolution(3*algo_param0%nVarD),stat=alloc_stat)
        end if

        if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_load9)",2)
        !! Lecture du chargement
        call get_loadList_xml(item(load_list,i),tag)
     end if
     if(.not. ok_i) then
        write(err_msg,fmt="(A,I0,A)") "Noeud 'Loading'",i+1, " (read_load10)"
        call amitex_abort(err_msg,1,0)
     end if
  end do
  
  call check_amitex_abort()
  nb_tps=sum(load0(:)%NIncr)
  n= maxval(load0(:)%NIncr)
  !! Temps des chargements
  call set_time()
  !-------------------
  !! Temps d'extraction
  
  !! Allocation du tableau des temps d'extraction
  allocate(extract0%tpsVTK(nb_tps),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_load11)",2)
  allocate(extract0%tpsZSTD(nb_tps),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_load12)",2)
  allocate(extract0%tpsSTD(nb_tps),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_load12b)",2)
  allocate(tmp_array(n),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_load13)",2)
  extract0%tpsVTK = .FALSE.
  extract0%tpsZSTD = .FALSE.
  extract0%tpsSTD = .FALSE.

  !  k=0
  allocate(extract0%numberZstd(getLength(load_list)))
  allocate(extract0%numberStd(getLength(load_list)))
  extract0%numberZstd = 0
  do i=0,getLength(load_list)-1
     call get_int_xml(item(load_list,i),"Tag",tag,1)
     if(tag < getLength(load_list)+1 .and. tag > 0 ) then

        !! Lecture de la liste des temps d'ecriture des moyennes par zone (.zstd)
        !! pour le chargement "tag"
        call getNodeList(item(load_list,i),sub_list,"Output_zone",-1,1)
        if(getLength(sub_list)>0) then
           if(getLength(sub_list)>1) then
              !! Plusieurs noeuds Output_zone dans le chargement i
              !! On ecrit un warning 
              write(err_msg,fmt="(2(A,I0),A)") " Plusieurs noeuds 'Output_zone' dans le &
                   &noeud de chargement ", i," (",getLength(sub_list),"). On ne considere &
                   &que le premier (read_load)"
              call amitex_abort(err_msg,-1,0)
           end if
           call get_int_xml(item(sub_list,0), "Number", extract0%numberZstd(tag), err)
           ! k+1: premier increment du chargement 
           k=0
           if(tag>1) k=sum(load0(1:tag-1)%NIncr)
           if(extract0%numberZstd(tag)<1) then
              !! On veut ecrire un nombre de pas de temps inferieur ou egal a 0
              !! On ecrit un warning 
              write(err_msg,fmt="(A,I0,A)") "Nombre de pas de temps devant etre ecrits &
                   &dans le .zstd (",extract0%numberZstd(tag),") inferieur ou egal a 0 : &
                   &on n'ecrira pas de moyenne par zone dans ce chargement (read_load)"
              call amitex_abort(err_msg,-1,0)
              extract0%tpsZSTD(k+1:k+load0(tag)%NIncr) = .FALSE.
           elseif(extract0%numberZstd(tag)>load0(tag)%NIncr) then
              !! On veut ecrire plus de pas de temps qu'il n'y en a dans ce chargement
              !! On ecrit un warning 
              write(err_msg,fmt="(A,I0,A)") "Nombre de pas de temps devant etre ecrits &
                   &dans le .zstd (",extract0%numberZstd(tag),") superieur au nombre de &
                   &pas dans ce chargement (read_load)"
              call amitex_abort(err_msg,-1,0)
              extract0%tpsZSTD(k+1:k+load0(tag)%NIncr) = .TRUE.
           else
              !! Sinon on choisit les pas ou on calcule les 
              !! moyennes par zone de la facon suivante :
              !! On choisit d'abord de sortir la moyenne par zone a la derniere iteration 
              !! du chargement. Les autres pas sont obtenus en decoupant 
              !! de maniere reguliere l'intervalle entre le debut et la fin.
              do j=1,extract0%numberZstd(tag)
                 extract0%tpsZSTD(k+int(j*load0(tag)%NIncr/extract0%numberZstd(tag)) )=.TRUE.
              end do
           end if
        end if

        !! Lecture de la liste des temps d'ecriture des moyennes par cellule et materiau (.std et .mstd)
        !! pour le chargement "tag"
        extract0%numberStd(tag)=load0(tag)%NIncr   !par defaut on impose Nincr (sortie a chaque increment)
        call getNodeList(item(load_list,i),sub_list,"Output_cell",-1,1)
        if(getLength(sub_list)>0) then
           if(getLength(sub_list)>1) then
              !! Plusieurs noeuds Output_zone dans le chargement i
              !! On ecrit un warning 
              write(err_msg,fmt="(2(A,I0),A)") " Plusieurs noeuds 'Output_cell' dans le &
                   &noeud de chargement ", i," (",getLength(sub_list),"). On ne considere &
                   &que le premier (read_load)"
              call amitex_abort(err_msg,-1,0)
           end if
           call get_int_xml(item(sub_list,0), "Number", extract0%numberStd(tag), err)
           if(extract0%numberStd(tag)<1) then
              !! On veut ecrire un nombre de pas de temps inferieur ou egal a 0
              !! On ecrit un warning 
              write(err_msg,fmt="(A,I0,A)") "Nombre de pas de temps devant etre ecrits &
                   &dans le .std (",extract0%numberStd(tag),") inferieur ou egal a 0 : &
                   &on n'ecrira pas de moyenne du dernier chargement (read_load)"
              call amitex_abort(err_msg,-1,0)
              extract0%numberStd(tag)=1
           elseif(extract0%numberStd(tag)>load0(tag)%NIncr) then
              !! On veut ecrire plus de pas de temps qu'il n'y en a dans ce chargement
              !! On ecrit un warning 
              write(err_msg,fmt="(A,I0,A)") "Nombre de pas de temps devant etre ecrits &
                   &dans .std et .mstd (",extract0%numberStd(tag),") superieur au nombre de &
                   &pas dans ce chargement (read_load)"
              call amitex_abort(err_msg,-1,0)
              extract0%numberStd(tag)=load0(tag)%NIncr
           end if
        end if
        !! On choisit d'abord de sortir la moyenne par zone a la derniere iteration 
        !! du chargement. Les autres pas sont obtenus en decoupant 
        !! de maniere reguliere l'intervalle entre le debut et la fin.
        k=0
        if(tag>1) k=sum(load0(1:tag-1)%NIncr)
        do j=1,extract0%numberStd(tag)
           extract0%tpsSTD(k+int(j*int(load0(tag)%NIncr/extract0%numberStd(tag))) )=.TRUE.
        end do

        !! Lecture de la liste des temps d'ecriture des VTK
        !! pour le chargement "tag"
        call getNodeList(item(load_list,i),sub_list,"Output_vtkList",-1,1)
        if(getLength(sub_list)>0)then
           if(getLength(sub_list)>1) then
              !! Plusieurs noeuds Output_vtkList dans le chargement i
              !! On ecrit un warning 
              write(err_msg,fmt="(2(A,I0),A)") " Plusieurs noeuds 'Output_vtkList' dans le &
                   &noeud de chargement ", i," (",getLength(sub_list),"). On ne considere &
                   &que le premier (read_load)"
              call amitex_abort(err_msg,-1,0)
           end if
           call extractDataContent(item(sub_list,0),vtk_list)
           tmp_array=0
           n=size(tmp_array)
           call get_int_vect( vtk_list, tmp_array, n, err)
           ! k+1: premier increment du chargement 
           k=0
           if(tag>1) k=sum(load0(1:tag-1)%NIncr)
           do j=1,n
              if( tmp_array(j)>load0(tag)%NIncr .or. tmp_array(j)<1)then
                 write(err_msg,fmt="(A,I0,A)") "Indice de temps d'extraction invalide ",tmp_array(j)," (read_load)"
                 call amitex_abort(err_msg,1,0)
              else
                 !! On ajoute les pas de temps concernes a la liste des pas de temps
                 !! pour lesquels on sort des VTK.
                 extract0%tpsVTK( k+tmp_array(j) )= .true.
              end if
           end do
        end if
        !      k=k+load0(tag)%NIncr
     end if
  end do
  if(allocated(tmp_array)) deallocate(tmp_array)
  call destroy(fi)
end subroutine read_load
!===================================================================================================


!===================================================================================================
!
!               SUBROUTINE SET_ZONE_TO_PRINT
!> Modification de la structure extract pour savoir quelles valeurs par zone
!! seront imprimees
!! \param[in]  zone_list: (nodeList) liste de noeuds donnant les materiaux
!!                         pour lesquels on doit extraire les moyennes par zone
!!                         et les variables internes pour lesquels il faut calculer les
!!                         moyennes par zone
!! \param[in]   nmate: (entier) nombre de materiaux
!!
!===================================================================================================
subroutine set_zone_to_print(zone_list,nmate)
  implicit none
  
  type(nodeList), pointer,intent(in)    :: zone_list
  type(nodeList), pointer               :: sub_list
  integer,intent(in)                    :: nmate
  type(node), pointer                   :: cur_node
  integer                               :: n,i,j,numM,alloc_stat,err
  logical                               :: ok_i
  character(len=200)                    :: err_msg
  integer,allocatable,dimension(:)      :: tmp_array 
  character(len=100000) :: var_int_list

  allocate(tmp_array(1));deallocate(tmp_array) ! avoid gcc-warning
  
  extract0%printMeanZ = .false.
  !! zone_list : liste des noeuds "Zone" (dans output)
  do i=0,getLength(zone_list)-1
    ok_i=.TRUE.
    numM=-1
    cur_node => item(zone_list,i)
    call get_int_xml(cur_node,"numM",numM,1)
    
    if(numM==-1)then  ! numero de materiau non renseigne
      ok_i=.FALSE.
      
    else if(numM>nmate .or. numM<1)then  ! numero de materiau invalide
      ok_i=.FALSE.
      write(err_msg,fmt="(A)") "Numero de materiau invalide (read_load)"
      call amitex_abort(err_msg,1)
      
    else !numero de materiau valide
      if(extract0%printMeanZ(numM)) then
        ! le numero de materiau a deja ete traite (l'initialisation est faite a 0)
        ok_i=.FALSE.
        write(err_msg,fmt="(A)") "Numero de materiau deja traite ", numM
      end if
      extract0%printMeanZ(numM) = .true.
    end if

    if(.not. ok_i)then
      ! s'il y a eu erreur, on precise le noeud
      write(err_msg, fmt="(A,I0,A)") "Noeud 'Zone' ", i+1," (read_load)"
      call amitex_abort(err_msg,1,0)
    end if

    !! n = nombre de variables internes du materiau
    !!     Rq : ci-dessous, extract0%varZSTD(numM)%val est un tableau de logical
    !!          de la taille du nombre de Var. int.
    if (algo_param0%Mechanics) then
    n=size(extract0%varZSTD(numM)%val)
    if( n >0) then
       allocate(tmp_array(n),stat=alloc_stat)
       if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_load)",2)
       call getNodeList(item(zone_list,i),sub_list,"VarIntList",-1,1)
       if(getLength(sub_list)>0)then
          call extractDataContent(item(sub_list,0),var_int_list)
          tmp_array=0
          n=size(tmp_array) 
          call get_int_vect( var_int_list, tmp_array, n, err)
          do j=1,n
             extract0%varZSTD(numM)%val( tmp_array(j) ) = .true.
          end do
       end if
    end if
    if(allocated(tmp_array))deallocate(tmp_array)
    end if
  end do
end subroutine set_zone_to_print
!===================================================================================================

!===================================================================================================
!
!               SUBROUTINE SET_INTVAR_TO_PRINT
!> Modification de la structure extract pour savoir quelles variables internes
!! seront imprimees
!! \param[in] sub_list : liste de noeuds ou sont precises les variables internes a ecrire
!! \param[in] nmate : nombre de materiaux dans le modele actuel
!!
!===================================================================================================
subroutine set_intvar_to_print(sub_list,nmate)
  implicit none
  
  type(nodeList), pointer,intent(in)    :: sub_list
  integer,intent(in)                    :: nmate
  type(node), pointer                   :: cur_node
  integer                               :: n,k,i,numM,alloc_stat,err
  logical                               :: ok_i
  character(len=200)                    :: err_msg,string
  integer,allocatable,dimension(:)      :: tmp_array 

  allocate(tmp_array(1)); deallocate(tmp_array) ! avoid gcc-warning
  
  !! sub_list : liste des noeuds "IntVar" (dans output)
  do i=0,getLength(sub_list)-1
    ok_i=.TRUE.
    numM=-1
    cur_node => item(sub_list,i)
    call get_int_xml(cur_node,"numM",numM,1)
    
    if(numM==-1)then  ! numero de materiau non renseigne
      ok_i=.FALSE.      
      write(err_msg,fmt="(A)") "Numero de materiau non renseigne (read_load)"
      call amitex_abort(err_msg,1)
    else if(numM>nmate .or. numM<1)then  ! numero de materiau invalide
      ok_i=.FALSE.
      write(err_msg,fmt="(A)") "Numero de materiau invalide (read_load)"
      call amitex_abort(err_msg,1)
    else !numero de materiau valide
      !if(.not. all(extract0%varVTK(numM)%val .eq. -1)) then
      if(.not. all(.not. extract0%varVTK(numM)%val) ) then
        ! le numero de materiau a deja ete traite (l'initialisation est faite a .false.)
        ok_i=.FALSE.
        write(err_msg,fmt="(A,I0)") "Numero de materiau deja traite ", numM
      end if
      !! n = nombre de variables internes du materiau
      n=size(extract0%VarVTK(numM)%val)
      if( n >0) then
        allocate(tmp_array(n),stat=alloc_stat)
        if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_load)",2)
        tmp_array=0 
        string = "false"
        call extractDataContent(cur_node,string)
        if(string .eq. "false")then  ! la liste est absente
           ok_i=.FALSE.
        else
           err=0
           ! on lit un tableau d'entiers de taille maximale 'n'
           call get_int_vect(string, tmp_array, n , err)
           ! n: taille de "List"
           if(err==1)then
              write(err_msg,fmt="(A)")"Trop d'arguments dans 'List' (read_load)"
              call amitex_abort(err_msg,1,0)
              ok_i=.FALSE.
           else
              do k=1,n
                 !aucun element de la liste ne doit avoir une valeur superieure au
                 !nombre de variables internes
                 if(tmp_array(k)> size(tmp_array) )then
                    write(err_msg,fmt="(A)")"Indice de variable interne non valide (read_load)"
                    call amitex_abort(err_msg,1,0)
                    ok_i=.FALSE.
                 else
                    ! la variable tmp_array(k) sera extraite
                    extract0%varVTK(numM)%val( tmp_array(k) ) = .true.
                 end if
              end do
           end if
        end if
        deallocate(tmp_array)
      end if
    end if
    if(.not. ok_i)then
      ! si il y a eu erreur, on precise le noeud
      write(err_msg, fmt="(A,I0,A)") "Noeud 'IntVar' ", i+1," (read_load)"
      call amitex_abort(err_msg,1,0)
    end if
  end do
end subroutine set_intvar_to_print
!===================================================================================================

!===================================================================================================
!
!               SUBROUTINE GETINITLOADEXT
!> Recupere les donnees pour initialiser la temperature et les parametres externes
!! seulement si la balise InitLoadExt est presente
!! \param[in]  file_node: (node) noeud representant le fichier de chargement
!!
!===================================================================================================
subroutine getInitLoadExt(file_node)

  implicit none
  
  type(node),pointer, intent(in)   :: file_node  
  type(nodeList), pointer          :: init_list, node_list
  type(node),pointer               :: cur_node      ! noeud courant
  character(len=300)               :: err_msg
  integer                          :: param_i, param_index
  logical,allocatable,dimension(:) :: is_param_initialized

  call GetNodeList(file_node, init_list, "InitLoadExt",-1,0)

  if(getLength(init_list)>1) then
     !! Si la balise InitLoadExt est trouve plusieurs fois on affiche une erreur
     write(err_msg,fmt="(A)") "Parametres externes initialises deux fois (getInitLoadExt)."
     call amitex_abort(err_msg,1)
  elseif(getLength(init_list)==0) then
     !! Si les parametres externes et la temperature ne sont pas initialises, on 
     !! met le booleen a .FALSE. et on precise que le nb de parametre externe est nul
     initValExt0%is_temp_initialized = .FALSE.
     initValExt0%nb_param = 0
     initValExt0%temp = 0
     !! On affiche un message pour signaler a l'utilisateur que ca n'a pas ete initialise
     write(err_msg,fmt="(A)") " Temperature et parametres externes non initialises. &
          &On suppose que la temperature est egale a 0 et qu'il y a un &
          &parametre externe valant 0 (getInitLoadExt)"
     call amitex_abort(err_msg,-1,0)
  else
     !! recherche la composante 'T' (temperature) dans la balise InitLoadExt
     node_list => getElementsByTagName(item(init_list,0),"T")
     !! Si la temperature n'est pas precisee on met le booleen a .FALSE.
     if(getLength(node_list)==0) then
        initValExt0%is_temp_initialized = .FALSE.
        !! On suppose que la temperature reste constante egale a 0 et on affiche un message
        write(err_msg,fmt="(A)") " Temperature non initialisee. &
             &On suppose qu'elle est egale a 0 (getInitLoadExt)"
        call amitex_abort(err_msg,-1,0)
        initValExt0%temp = 0
     elseif(getLength(node_list)>1) then
        !! Si il y a plusieurs temperatures precisees, erreur
        write(err_msg,fmt="(A)") "Plusieurs temperatures precisees (getInitLoadExt)"
        call amitex_abort(err_msg,1,0)
     else
        !! Sinon on met le booleen a .TRUE. et on recupere la valeur
        initValExt0%is_temp_initialized = .TRUE.
        cur_node => item(node_list,0)
        call get_real_xml(cur_node,"Value", initValExt0%temp,0,1)
     end if
     !! recherche la composante 'Param' (parametre externe) dans la balise InitLoadExt
     node_list => getElementsByTagName(item(init_list,0),"Param")
     initValExt0%nb_param = getLength(node_list)
     if(initValExt0%nb_param>0) then 
        !! Si le nombre de parametres initialises est non nul
        allocate(is_param_initialized(initValExt0%nb_param))
        is_param_initialized = .FALSE.
        !! On definit un tableau pour savoir si chaque parametre a ete initialise
        allocate(initValExt0%param_values(initValExt0%nb_param))
        do param_i=1,initValExt0%nb_param
           cur_node => item(node_list,param_i-1)
           !! On recupere l'indice du parametre
           call get_int_xml(cur_node,"Index",param_index,0,1)
           if(is_param_initialized(initValExt0%nb_param)) then
              !! Si le parametre a deja ete initialise erreur
              write(err_msg,fmt="(A,I0,A)") "Parametre deja initialise indice ",param_index," (getInitLoadExt)"
              call amitex_abort(err_msg,2,0)
           else
              is_param_initialized(param_index) = .TRUE.
              !! On verifie que cet indice est entre les bornes
              if(param_index<1 .OR. param_index > initValExt0%nb_param) then
                 write(err_msg,fmt="(A,I0,A)") "Indice non valide ",param_index," (getInitLoadExt)"
                 call amitex_abort(err_msg,2,0)
              end if
              call get_real_xml(cur_node,"Value", initValExt0%param_values(param_index),0,1)
           end if
        end do
        deallocate(is_param_initialized)
     else
        !! Dans la fonction behavior, on suppose qu'il y a un parametre
        !! externe valant 0. On affiche un message pour informer l'utilisateur
        write(err_msg,fmt="(A)") " Parametres externes non initialises. &
             &On suppose qu'il y a un parametre externe valant 0 (getInitLoadExt)"
        call amitex_abort(err_msg,-1,0)
     end if
  end if
end subroutine getInitLoadExt

!===================================================================================================
!
!               SUBROUTINE GET_LOADLIST_XML
!> Definition du chargement pour chaque composante
!! Lit et extrait pour chaque composante (xx, yy, zz, xy, xz, yz [et yx, zx, zy 
!!                                        en grandes transformations]):
!!      - le type de pilotage (contrainte ou deformation)
!!      - le type d'evolution (constante, lineaire)
!!      - la valeur finale (en contrainte ou deformation, selon le pilotage)
!!      - Les valeurs et les types d'evolution pour la temperature et les
!!        parametres externes
!! Et aussi les parametres associes a :
!!      - chargement en direction de contrainte imposee
!!      - chargement en gradient de gradient de deplacement
!!      - interruption de chargement implementee par l'utilisateur
!!      
!! \param[in]  loading_node: (node) noeud de chargement "Loading"
!! \param[in]  tag: (entier) numero du chargement
!!
!===================================================================================================
subroutine get_loadList_xml(loading_node, tag)

  implicit none
  
  type(node),pointer, intent(in)     :: loading_node  ! loading node
  integer,intent(in)                 :: tag           ! loading tag
  type(nodeList),pointer             :: node_list,elem_list,userParam     ! liste de sous noeud
  type(node),pointer                 :: cur_node, cur_elem      ! noeud courant
  character(len=15)                  :: string        
  character(len=200)                 :: err_msg
  character(len=1),dimension(3)      :: comp=["x","y","z"]
  integer,dimension(3)               :: dir
  integer                            :: alloc_stat,i,j,k, ind, nuser_param

  call getLoadTime(tag, loading_node)
  call getLoadParam(tag,loading_node)

  ! cas general mecanique
  if (algo_param0%Mechanics) then
  call get_load(tag,loading_node,"xx",1)
  call get_load(tag,loading_node,"yy",2)
  call get_load(tag,loading_node,"zz",3)
  call get_load(tag,loading_node,"xy",4)
  call get_load(tag,loading_node,"xz",5)
  call get_load(tag,loading_node,"yz",6)
  if(algo_param0%HPP) then
     call getNodeList(loading_node,node_list,"yx",-1,1)
     if(getLength(node_list)>0)then
        write(err_msg,fmt="(A,I0,A)") "Composante ""yx"" presente en HPP au noeud &
             &'Loading' ",tag,". Elle ne sera pas lue (get_loadList_xml)"
        call amitex_abort(err_msg,-1,0)
     end if
     call getNodeList(loading_node,node_list,"zy",-1,1)
     if(getLength(node_list)>0)then
        write(err_msg,fmt="(A,I0,A)") "Composante ""zy"" presente en HPP au noeud &
             &'Loading' ",tag,". Elle ne sera pas lue (get_loadList_xml)"
        call amitex_abort(err_msg,-1,0)
     end if
     call getNodeList(loading_node,node_list,"zx",-1,1)
     if(getLength(node_list)>0)then
        write(err_msg,fmt="(A,I0,A)") "Composante ""zx"" presente en HPP au noeud &
             &'Loading' ",tag,". Elle ne sera pas lue (get_loadList_xml)"
        call amitex_abort(err_msg,-1,0)
     end if
  else
     call get_load(tag,loading_node,"yx",7)
     call get_load(tag,loading_node,"zx",8)
     call get_load(tag,loading_node,"zy",9)
  end if
  end if

  ! Recherche si DirStress defini sur CAUCHY ou PK1
  if (algo_param0%Mechanics .AND. (.not. algo_param0%HPP) .AND. load0(tag)%DirStress_flag) then
  call getNodeList(loading_node,node_list,"DirStress",1,1)
  if(associated(node_list))then
     cur_node => item(node_list,0)
     string=""
     call get_str_xml(cur_node,"Type", string,1,0)
     if(string=="cauchy") then
        load0(tag)%DirStress_cauchy=.true.
     elseif(string=="pk1") then
        load0(tag)%DirStress_cauchy=.false.
     else
        write(err_msg,fmt="(3A,I0)") &
        "Type de 'DirStress' inconnu : '",string,"' (get_loadList_xml) noeud 'Loading' ",tag
        call amitex_abort(err_msg,1,0)
     end if
  else
     write(err_msg,fmt="(A,I0,A)") &
           "Composante 'DirStress', noeud 'Loading' ",tag," absente (get_loadList_xml)"
     call amitex_abort(err_msg,1,0)
  end if
  end if ! fin condition Meca + Gdef + DirStress

  ! On initialise DirStress2
  if (algo_param0%Mechanics .AND. load0(tag)%DirStress_flag) then
     load0(tag)%DirStress2 = load0(tag)%DirStress
  end if 
  
  ! Cas des Gradients de deformation imposes
  if (algo_param0%Mechanics .AND. algo_param0%HPP) then
     do k=1,3
     do j=1,3
     do i=1,3
        dir = [i,j,k]
        call get_load_graddef(tag,loading_node,  &
            "gradgradU_"//comp(dir(1))//comp(dir(2))//"_"//comp(dir(3)),dir)
     end do
     end do
     end do   
  end if
  
  ! Cas User_interruption (ADAPTATION lecture parametres user P_real dans param_algo_mod.f90) 
  if (algo_param0%Mechanics) then
     userParam=> getElementsByTagName(loading_node,"User_interruption")
     nuser_param=getLength(userParam)    
     if (nuser_param>0) then
     load0(tag)%User_interruption_flag=.true.
     !> Allocation
     allocate(load0(tag)%User_interruption(nuser_param),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort(&
        "Espace memoire disponible insuffisant (get_LoadList_xml) pour allouer load(tag)%User_interruption",2,0)
     !> lecture de toutes les valeurs de parametres (avec traitement des erreurs)
     do i=0,nuser_param-1
         cur_node => item(userParam,i)
         call get_int_xml(cur_node,"Index",ind,1)
         if (ind>nuser_param .OR. ind<1) then
            write(err_msg,fmt="(A)") &
            " Index <=0, or higher than the number of parameters <User_interruption... /> in .xml load file (get_LoadList_xml)"
            call amitex_abort(err_msg, 2,0)
         end if
         call get_real_xml(cur_node,"Value",load0(tag)%User_interruption(ind),1)
     end do
     end if
  end if
  
  ! cas general diffusion
  if (algo_param0%Diffusion) then
  call get_load(tag,loading_node,"x0",1)
  call get_load(tag,loading_node,"y0",2)
  call get_load(tag,loading_node,"z0",3)
  end if 

  ! Restart
  elem_list => getElementsByTagName(loading_node, "Restart")
  if(getLength(elem_list)==1) then
     cur_elem => item(elem_list, 0)       
     string = "Every"
     call get_int_xml(cur_elem,string,load0(tag)%restart_save_period,1)
     load0(tag)%restart_save=.true.
  end if

end subroutine  get_loadList_xml


!===================================================================================================
!
!                    SUBROUTINE GETLOADTIME
!
!> Lecture des parametres de temps pour le chargement
!!  Lit et extrait pour le chargement donne:
!!       - le nombre d'increment de temps
!!       - le type de discretisation ( defini par l'utilisateur, lineaire)
!!       - la liste de temps si defini par l'utilisateur, le temps en fin
!!           de chargement sinon
!!
!! \param[in]  load_nb: (entier) numero du chargement
!! \param[in]  loading_node: (type node) noeud de chargement "Loading"
!!
!!
!===================================================================================================
subroutine getLoadTime(load_nb,loading_node)

  implicit none
  
  integer, intent(in)                   :: load_nb  ! numero du chargement
  type(node),intent(in),pointer         :: loading_node  ! noeud de chargement
  real(mytype)                          :: last_time, current_time
  type(nodeList),pointer                :: time_nodeList
  type(node),pointer         :: time_node ! noeud d'ou on extrait les temps

  integer                               :: j  ! numero de pas de temps local
  character(len=100000)                 :: time_list
  character(len=300)                    :: string,err_msg
  integer                               :: n, alloc_stat,err
  logical                               :: error
  
  error=.FALSE.
  last_time = 0._mytype
  if(load_nb .GE. 2)  last_time = load0(load_nb-1)%time(load0(load_nb-1)%NIncr)
  
  call getNodeList(loading_node,time_nodeList,"Time_Discretization",1,1)
  if(associated(time_nodeList)) then
     time_node => item(time_nodeList,0)
     n=-1
     ! nombre d'increment du chargement
     call get_int_xml(time_node,"Nincr",n,1)
     if(n.LE.0) then
        error=.TRUE.
     else
        ! allocation du tableau de temps
        allocate(load0(load_nb)%time(0:n),stat=alloc_stat)
        if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (getLoadTime)",2)
        load0(load_nb)%NIncr=n
     end if
     ! type de discretisation
     string=""
     call get_str_xml(time_node ,"Discretization",string,1)

     if(string=="user")then
        ! discretisation definie par l'utilisateur
        call getNodeList(loading_node,time_nodeList,"Time_List",1,1)
        ! liste de temps de l'utilisateur
        call extractDataContent(item(time_nodeList,0),time_list)
        ! on lit un tableau d'entiers de taille maximale 'n'
        call get_real_vect(time_list, load0(load_nb)%time(1:n), n , err)
        ! n: taille de "List"
        !! On verifie que les pas de temps sont positifs
        do j=1,load0(load_nb)%NIncr
           current_time = load0(load_nb)%time(j)
           if(current_time .LE. last_time) then
              write(err_msg,fmt="(2(A,E15.8),A)") "Temps suivant inferieur au temps courant : ",&
                   current_time," < ",last_time," (getLoadTime)"
              call amitex_abort(trim(err_msg),2,0)
           end if
           last_time=current_time
        end do
     else if(string=="linear")then
        ! Discretisation lineaire
        ! Temps en fin de  chargement
        call get_real_xml(time_node, "Tfinal", load0(load_nb)%time(n),1)
        !! On verifie que les pas de temps sont positifs
        if(load0(load_nb)%time(n) .LE. last_time) then
           write(err_msg,fmt="(2(A,E15.8),A)") "Temps final inferieur ou egal au temps courant : ",&
                load0(load_nb)%time(n)," <= ",last_time," (getLoadTime)"
           call amitex_abort(trim(err_msg),2,0)
        end if
        last_time = load0(load_nb)%time(n)
     else if(string=="")then
        error=.TRUE.
     end if
     if(error) then
        write(string,fmt="(A,I0,A)") "Erreur noeud 'Loading' ",load_nb," (getLoadTime)"
        call amitex_abort(string,1)
     else
        call repartition_discretization(string, load0(load_nb)%discretization)
        if(load0(load_nb)%discretization==-1) then
           write(string,fmt="(A,I0,A)") "Erreur noeud 'Loading' ",load_nb," (getLoadTime)"
           call amitex_abort(string,1,0)
        end if
     end if
  else
     write(string,fmt="(A,I0,A)") "Erreur noeud 'Loading' ",load_nb," (getLoadTime)"
     call amitex_abort(string,1,0)
  end if
end subroutine getLoadTime

!===================================================================================================
!
!                     SUBROUTINE GETLOADPARAM
!
!> Lecture des parametres externes
!! Lit dans le fichier xml, pour un noeud de chargement et une composante donnes,
!!       - La temperature imposee et le type d'evolution souhaite
!!       - Les parametres externes imposes et le type d'evolution souhaite
!! Si la temperature n'est pas precisee on suppose qu'elle est egale a 0
!! Si aucun parametre n'est precise on suppose que c'est un tableau de taille 1 
!! avec la valeur 0
!! \param[in]  load_nb: (entier)  numero du chargement
!! \param[in]  loading_node: (type node) noeud de chargement (fichier xml) "Loading"
!!
!!
!===================================================================================================

subroutine getLoadParam(load_nb, loading_node)
  
  implicit none

  integer, intent(in)            :: load_nb      !numero de chargement, numero de composante
  type(node),pointer, intent(in) :: loading_node ! noeud de chargement

  type(node),pointer            :: cur_node      ! noeud courant
  type(nodeList),pointer        :: node_list     ! liste de sous noeud
  character(len=15)             :: string        ! pilotage
  character(len=200)             :: err_msg       ! message d'erreur
  integer                       :: param_i, param_index, alloc_stat
  logical :: error

  string=""
  error=.FALSE.

  ! recherche la composante 'T' (temperature) dans le chargement
  node_list => getElementsByTagName(loading_node,"T")
  if(getLength(node_list)==0 .AND. initValExt0%is_temp_initialized) then
     error=.TRUE.
     write(err_msg,fmt="(A)") "Temperature initialisee non definie (getLoadParam)"
     call amitex_abort(err_msg,1,0)
  elseif(getLength(node_list)>0 .AND. .NOT. initValExt0%is_temp_initialized) then
     error=.TRUE.
     write(err_msg,fmt="(A)") "Temperature definie et non initialisee (getLoadParam)"
     call amitex_abort(err_msg,1,0)
  end if
  if(getLength(node_list)==0) then
    load0(load_nb)%tvalue=0
    load0(load_nb)%tevolution=0
  elseif(getLength(node_list)>1) then
     write(err_msg,fmt="(A)") "Plusieurs temperatures precisees (getLoadParam)"
     call amitex_abort(err_msg,1,0)
     error=.TRUE.
  else
     cur_node => item(node_list,0)
     call get_str_xml(cur_node,"Evolution",string,1,0)
     if(string=="linear")then
        !! On initialise tvalue a huge(1._mytype)
        load0(load_nb)%tvalue=huge(1._mytype)
        !! cas lineaire: recherche la valeur imposee (en fin de chargement)
        call get_real_xml(cur_node,"Value", load0(load_nb)%tvalue,1,0)
        if(load0(load_nb)%tvalue==huge(1._mytype))then
           error=.TRUE.
        else
           call repartition_evolution(string,load0(load_nb)%tevolution)
        end if 
     else if(string=="constant")then
        call repartition_evolution(string,load0(load_nb)%tevolution)
     else if(string=="")then
        error=.TRUE.
        write(err_msg,fmt="(A)") "Type d'evolution non specifie (getLoadParam)"
        call amitex_abort(err_msg,1,0)
     else
        error=.TRUE.
        write(err_msg,fmt="(3A)") "Evolution non implementee : ",string," (getLoadParam)"
        call amitex_abort(err_msg,1,0)
     end if
  endif 
  if(error)then
     write(err_msg,fmt="(A,I0,A)") "Temperature noeud 'Loading' ",load_nb," (getLoadParam)"
     call amitex_abort(err_msg,-1,0)
  end if
  error=.FALSE.
  !! Fin du calcul de la temperature
  ! recherche la composante 'Param' (parametres externes) dans le chargement
  node_list => getElementsByTagName(loading_node,"Param")
  if(initValExt0%nb_param/=getLength(node_list)) then
     error=.TRUE.
     write(err_msg,fmt="(A)") "Nombre de parametres externes incoherent avec l'initialisation (getLoadParam)"
     call amitex_abort(err_msg,1,0)
  else
     allocate(load0(load_nb)%pvalue(initValExt0%nb_param),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (readLoadParam)",2)
     allocate(load0(load_nb)%pevolution(initValExt0%nb_param),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (readLoadParam)",2)
     load0(load_nb)%pevolution = -2
     do param_i=1,initValExt0%nb_param
        cur_node => item(node_list,param_i-1)
        call get_int_xml(cur_node,"Index",param_index,0,1)
        if(param_index<1 .OR. param_index > initValExt0%nb_param) then
           error=.TRUE.
           write(err_msg,fmt="(2(A,I0),A)") "Indice non valide ",param_index," noeud 'Loading' ",load_nb," (getLoadParam)"
           call amitex_abort(err_msg,2,0)
        elseif(load0(load_nb)%pevolution(param_index) /= -2) then
           error=.TRUE.
           write(err_msg,fmt="(2(A,I0),A)") "Indice deja entre ",param_index," noeud 'Loading' ",load_nb," (getLoadParam)"
           call amitex_abort(err_msg,1,0)
        end if
        call get_str_xml(cur_node,"Evolution",string,1,0)
        if(string=="linear")then
           !! On initialise tvalue a huge(1._mytype)
           load0(load_nb)%pvalue(param_index)=huge(1._mytype)
           !! cas lineaire: rec2herche la valeur imposee (en fin de chargement)
           call get_real_xml(cur_node,"Value", load0(load_nb)%pvalue(param_index),1,0)
           if(load0(load_nb)%tvalue==huge(1._mytype))then
              error=.TRUE.
           else
              call repartition_evolution(string,load0(load_nb)%pevolution(param_index))
           end if
        else if(string=="constant")then
           call repartition_evolution(string,load0(load_nb)%pevolution(param_index))
        else if(string=="")then
           error=.TRUE.
           write(err_msg,fmt="(A)") "Type d'evolution non specifie (getLoadParam)"
           call amitex_abort(err_msg,1,0)
        else
           error=.TRUE.
           write(err_msg,fmt="(3A)") "Evolution non implementee : ",string," (getLoadParam)"
           call amitex_abort(err_msg,1,0)
        end if
     end do
  endif
  if(error)then
     write(err_msg,fmt="(A,I0,A)") "Parametre externe noeud 'Loading' ",load_nb," (getLoadParam)"
     call amitex_abort(err_msg,-1,0)
  end if



end subroutine getLoadParam



!===================================================================================================
!
!                     SUBROUTINE GET_LOAD
!
!> Lecture du chargement
!! Lit dans le fichier xml, pour un noeud de chargement et une composante donnes,
!!       - le type de pilotage
!!       - le type d'evolution (de la valeur imposee)
!!       - la valeur imposee (si necessaire) 
!!
!! \param[in]  direction: (entier)  numero du chargement
!! \param[in]  load_nb: (entier)  numero de l'increment
!! \param[in]  loading_node: (type node) noeud de chargement (fichier xml) "Loading"
!! \param[in]  comp: (chaine de caracteres) composante traitee
!!
!!
!!pilotage: contrainte: 1
!!          deformation: 0
!!
!!ordre des composantes: 11,22,33,12,13,23,(+Voigt en HPP, 21,31,32 en GD)
!===================================================================================================

subroutine get_load(load_nb, loading_node, comp,direction)
  
  implicit none

  integer, intent(in)            :: load_nb, direction   !numero de chargement, numero de composante
  type(node),pointer, intent(in) :: loading_node ! noeud de chargement
  character(len=2),intent(in)    :: comp         ! nom de la composante

  type(node),pointer            :: cur_node      ! noeud courant
  type(nodeList),pointer        :: node_list     ! list de sous noeud
  character(len=15)             :: string        ! pilotage
  character(len=50)             :: err_msg       ! message d'erreur

  logical :: error, diffusion,test1,test2,test3

  string=""
  error=.FALSE.
  
  ! teste si le chargement correspond a un chargement de diffusion
  diffusion=.false.
  if (comp=="x0" .or. comp=="y0" .or. comp=="z0") diffusion=.true.

  ! recherche la composante 'comp' dans le chargement
  call getNodeList(loading_node,node_list,comp,1,1)

  if(associated(node_list))then
     cur_node => item(node_list,0)

     ! identification du type de pilotage 
     string=""
     call get_str_xml(cur_node,"Driving", string,1,0)
     if(string=="strain" .and. (diffusion .eqv. .false.)) then
        load0(load_nb)%driving(direction)=0
     else if(string=="stress".and. (diffusion .eqv. .false.)) then
        load0(load_nb)%driving(direction)=1
     else if(string=="gradd" .and. (diffusion .eqv. .true.)) then
        load0(load_nb)%drivingD(direction)=0
        diffusion=.true.
     else if(string=="fluxd" .and. (diffusion .eqv. .true.)) then
        load0(load_nb)%drivingD(direction)=1
        diffusion=.true.
     else if(string=="")then
        error=.TRUE.
     else
        write(err_msg,fmt="(3A)") "Pilotage inconnu : '",string,"' (get_load)"
        call amitex_abort(err_msg,1,0)
        error=.TRUE.
     end if

     ! identification d'un chargement a multiaxialite de contrainte imposee dirStress
     if(algo_param0%Mechanics) then
        load0(load_nb)%DirStress(direction) = huge(1._mytype)
        call get_real_xml(cur_node,"DirStress", load0(load_nb)%DirStress(direction),-2,0) 
        if (load0(load_nb)%DirStress(direction) .NE. huge(1._mytype)) load0(load_nb)%DirStress_flag=.true. 
        test1 = .not. load0(load_nb)%DirStress_flag
        test2 = load0(load_nb)%DirStress_flag .and. load0(load_nb)%driving(direction)==0
        test3 = load0(load_nb)%DirStress_flag .and. load0(load_nb)%driving(direction)==1
     end if

     !CAS MECANIQUE
     ! cas "usuel" : DirStress_flag = .false. 
     !           OU  (DirStress_flag = .true. ET composante imposee en deformation)
     if (algo_param0%Mechanics .and. (test1 .or. test2)) then
     string=""
     call get_str_xml(cur_node,"Evolution",string,1,0)
     if(string=="linear".and. (diffusion .eqv. .false.))then
        !! On initialise mvalue a huge(1._mytype)
        load0(load_nb)%mvalue(direction)=huge(1._mytype)
        ! cas lineaire: recherche la valeur imposee (en fin de chargement)
        call get_real_xml(cur_node,"Value", load0(load_nb)%mvalue(direction),1,0)
        if(load0(load_nb)%mvalue(direction)==huge(1._mytype))then
           error=.TRUE.
        else
           call repartition_evolution(string,load0(load_nb)%mevolution(direction))
        end if
     else if(string=="constant" .and. (diffusion .eqv. .false.)) then
        call repartition_evolution(string,load0(load_nb)%mevolution(direction))
     else if(string=="")then
        error=.TRUE.
        write(err_msg,fmt="(A)") "Type d'evolution non specifie (get_load)"
        call amitex_abort(err_msg,1,0)
     else
        error=.TRUE.
        write(err_msg,fmt="(3A)") "Evolution non implementee : ",string," (get_load)"
        call amitex_abort(err_msg,1,0)
     end if
     end if

     ! cas "particulier" : DirStress_flag = .true. ET composante imposee en contrainte
     !   => on impose une evolution bidon" (constante nulle)
     if(algo_param0%Mechanics .and. test3) then
        call repartition_evolution("constant",load0(load_nb)%mevolution(direction))
        load0(load_nb)%mvalue(direction) = 0.
     end if

     !CAS DIFFUSION
     if (algo_param0%Diffusion) then
     string=""
     call get_str_xml(cur_node,"Evolution",string,1,0)
     if(string=="linear" .and. (diffusion .eqv. .true.)) then
        !! On initialise dvalue a huge(1._mytype)
        load0(load_nb)%dvalue(direction)=huge(1._mytype)
        ! cas lineaire: recherche la valeur imposee (en fin de chargement)
        call get_real_xml(cur_node,"Value", load0(load_nb)%dvalue(direction),1,0)
        if(load0(load_nb)%dvalue(direction)==huge(1._mytype))then
           error=.TRUE.
        else
           call repartition_evolution(string,load0(load_nb)%devolution(direction))
        end if
     else if(string=="constant" .and. (diffusion .eqv. .true.)) then
        call repartition_evolution(string,load0(load_nb)%devolution(direction))
     else if(string=="")then
        error=.TRUE.
        write(err_msg,fmt="(A)") "Type d'evolution non specifie (get_load)"
        call amitex_abort(err_msg,1,0)
     else
        error=.TRUE.
        write(err_msg,fmt="(3A)") "Evolution non implementee : ",string," (get_load)"
        call amitex_abort(err_msg,1,0)
     end if
     end if

     ! GESTION ERREUR
     if(error)then
        write(err_msg,fmt="(3A,I0,A)") "Component ", comp, " node 'Loading' ",load_nb," (get_load)"
        call amitex_abort(err_msg,-1,0)
     end if
  else
     write(err_msg,fmt="(3A,I0,A)") "Component ", comp,", node 'Loading' ",load_nb," absent (get_load)"
     call amitex_abort(err_msg,-1,0)
  end if

end subroutine get_load

!------------- SAME FUNCTION FOR DEFORMATION GRADIENTS

subroutine get_load_graddef(load_nb, loading_node, comp,dir)
  
  implicit none

  integer, intent(in)            :: load_nb      ! numero de chargement
  integer,dimension(3)           :: dir          ! indice de composante
  type(node),pointer, intent(in) :: loading_node ! noeud de chargement
  character(len=14),intent(in)   :: comp         ! nom de la composante

  type(node),pointer             :: cur_node     ! noeud courant
  type(nodeList),pointer         :: node_list    ! list de sous noeud
  character(len=15)              :: string       ! pilotage
  character(len=100)             :: err_msg      ! message d'erreur
  real(mytype)                   :: val0          ! temporary value


  ! recherche la composante 'comp' dans le chargement
  ! on ne specifie pas le nombre de noeud (-1)
  call getNodeList(loading_node,node_list,comp,-1,-2)

  if (getLength(node_list) > 1 ) then
     write(err_msg,fmt="(3A)") "Error : component ", comp," appears more than once in load.xml (get_load)"
     call amitex_abort(trim(err_msg),1,0)
  end if
  
  if(getLength(node_list) == 1) then
     n_gradgradU0 = 27
     cur_node => item(node_list,0)

     ! read Evolution and value
     string=""
     call get_str_xml(cur_node,"Evolution",string,1,0)

     if(string=="linear") then
        ! Initialize the component of gvalue to huge(1._mytype)
        val0=huge(1._mytype)
        ! cas lineaire: recherche la valeur imposee (en fin de chargement)
        call get_real_xml(cur_node,"Value", val0,1,0)
        if(val0==huge(1._mytype)) then
           write(err_msg,fmt="(A)") "Error : Value not affected for a deformation gradient component (get_load0)"
           call amitex_abort(trim(err_msg),1,0)
        else
           load0(load_nb)%gvalue(dir(1),dir(2),dir(3)) = val0 + load0(load_nb)%gvalue(dir(1),dir(2),dir(3))
           call repartition_evolution(string,load0(load_nb)%gevolution(dir(1),dir(2),dir(3)))
        end if
     else if(string=="constant") then
        call repartition_evolution(string,load(load_nb)%gevolution(dir(1),dir(2),dir(3)))
     else if(string=="")then
        write(err_msg,fmt="(A)") "'Evolution' to be specified (get_loadgraddef)"
        call amitex_abort(trim(err_msg),1,0)
     else
        write(err_msg,fmt="(3A)") "Unrecognized 'Evolution' : ",string," (get_loadgraddef)"
        call amitex_abort(trim(err_msg),1,0)
     end if
  end if

end subroutine get_load_graddef

!===================================================================================================
!
!              SUBROUTINE REPARTITION_DISCRETIZATION
!
!  
!>  renvoie un entier en fonction du type de discretisation du temps de chargement souhaite
!!
!! \param[in] discretization: (chaine de caracteres) nom de la discretisation
!!
!! \param[out]   n: (entier) code de la discretisation
!!                 - -1 discretisation inconnue ou non implementee
!!                 - 0 discretisation definie par l'utilisateur
!!                 - 1 discretisation lineaire du temps
!!
!===================================================================================================
subroutine repartition_discretization(discretization, n)

  implicit none
  
  character(len=*), intent(in)  :: discretization
  integer, intent(out)          :: n
  character(len=200)            :: err_msg
  
  n=-1
  if(trim(discretization)=="user")then
    n=0
  else if(trim(discretization)=="linear") then
    n=1
  else
    write(err_msg,fmt="(3A)") "Loi de discretisation non implementee : ",&
         trim(discretization)," (repartition_discretization)"
    call amitex_abort(err_msg,1,0)
  end if
  

end subroutine repartition_discretization
!===================================================================================================


!===================================================================================================
!
!              SUBROUTINE REPARTITION_EVOLUTION
!
!  
!>  renvoie un entier en fonction du type d'evolution du chargement souhaite 
!!
!! \param[in]   evol: (chaine de caracteres) nom de l'evolution
!!
!! \param[out]   n: (entier) code de la discretisation
!!                 - -1 evolution inconnue ou non implementee
!!                 - 0  constante
!!                 - 1  lineaire
!!
!===================================================================================================
subroutine repartition_evolution(evol, n)

  implicit none
  
  character(len=*), intent(in)  :: evol
  integer, intent(out)          :: n
  character(len=200)               :: err_msg
  
  if(trim(evol)=="constant") then
    n=0
  else if(trim(evol)=="linear") then
    n=1
  else
    n=-1
    write(err_msg,fmt="(3A)") "Unimplemented 'Evolution' : ",&
         trim(evol),"(repartition_evolution)"
    call amitex_abort(err_msg,1,0)
  end if
end subroutine repartition_evolution


!===================================================================================================
!
!                     SUBROUTINE GET_CURRENT_NB_PARAM
!
!
!>  renvoie le nombre de parametres externes au chargement courant
!!
!! \param[in]   load_nb  (entier) numero du chargement
!!
!! \param[out]  nb_param  (entier) nombre de parametres externes
!!
!!
!===================================================================================================
subroutine get_current_nb_param(load_nb, nb_param) ! gcc-warning accepte (unused)
  
  implicit none
  
  integer, intent(in)                         :: load_nb ! numero du chargement
  integer,intent(out)     :: nb_param

  nb_param = size(load0(load_nb)%pevolution)

end subroutine get_current_nb_param

!===================================================================================================
!
!                     SUBROUTINE GET_CURRENT_LOADING
!
!
!>  renvoie le chargement a l'instant courant
!!
!! \param[in]   load_nb  (entier) numero du chargement
!! \param[in]   incr  (entier) numero de l'increment de temps au sein du chargement
!! \param[in]   defMoy  (reel) deformation moyenne a l'instant precedent 
!!                      (dimension algo_param%nTensDef)
!! \param[in]   sigMoy  (reel) contrainte moyenne a l'instant precedent 
!!                      (dimension algo_param%nTensDef)
!! \param[in]   local_load0  (reel) chargement local a l'instant precedent
!!                  - local_load(1,algo_param%nTensDef) : 
!!                      type de pilotage (0 deformation, 1 contrainte)
!!                  - local_load(algo_param%nTensDef+1,2*algo_param%nTensDef) : 
!!                      valeur imposee pour la contrainte ou la deformation
!!                  - local_load(2*algo_param%nTensDef+1) : temperature imposee
!!                  - local_load(next nb_param indices) : parametres externes imposes
!!                  - local_load(next 27 indices) : gradgradU components
!! \param[in]   t  (reel) tableau des temps aux instants t-3, t-2, t-1
!!
!!
!!
!!
!! \param[out]  local_load  (reel) chargement local actuel
!!                - local_load(1,algo_param%nTensDef) : 
!!                    type de pilotage (0 deformation, 1 contrainte)
!!                - local_load(algo_param%nTensDef+1,2*algo_param%nTensDef) : valeur imposee
!!                - local_load(2*algo_param%nTensDef+1) : temperature imposee
!!                - local_load(next nb_param indices) : parametres externes imposes
!!                - local_load(next 27 indices) : gradgradU components
!! \param[out]  t  (reel) tableau des temps aux instants t-2, t-1, t
!!
!!
!===================================================================================================
subroutine get_current_loading(load_nb, incr, defMoy, sigMoy, t , local_load0, local_load)
  
  implicit none
  
  integer, intent(in)                            :: load_nb, incr  ! numero du chargement et increment 
  real(mytype),dimension(algo_param%nTensDef), intent(in)       :: defMoy,sigMoy
  real(mytype),dimension(:), intent(in)          :: local_load0
  real(mytype),dimension(:), intent(out)         :: local_load
  real(mytype),dimension(-2:0),intent(inout)     :: t
  integer                                        :: ndeb,nfin

  integer :: j

  ! on initialise local_load a local_load0 (si chargement Constant, on ne fait rien)
  local_load = local_load0

  t(-2)=t(-1)         !t-2
  t(-1)=t(0)          !t-1
  t(0)=load(load_nb)%time(incr)
  local_load(1:algo_param%nTensDef)=load(load_nb)%driving
  !valeur du chargement local
  do j=1,algo_param%nTensDef
    if(load(load_nb)%driving(j)==1)then
      !pilotage en contrainte
      call get_loading_value(load_nb, incr, j, sigMoy(j), local_load(j+algo_param%nTensDef))
    else
      !pilotage en deformation
      call get_loading_value(load_nb, incr, j, defMoy(j), local_load(j+algo_param%nTensDef))
    end if
  end do

  ! temperature
  call get_temp_value(load_nb, incr, local_load(2*algo_param%nTensDef+1))

  ! external parameters
  if(initValExt%nb_param>0) then
    ndeb = 2*algo_param%nTensDef+1+1
    nfin = 2*algo_param%nTensDef+1+initValExt%nb_param
    call get_PREDEF_value(load_nb, incr,local_load(ndeb:nfin))
  end if

  ! gragradU components
  if(n_gradgradU==27) then !test pas tres joli <=> test gradgradU impose
    ndeb = 2*algo_param%nTensDef+1+initValExt%nb_param+1
    nfin = 2*algo_param%nTensDef+1+initValExt%nb_param+27
    call get_gradgradU_value(load_nb, incr,local_load(ndeb:nfin))
  end if

end subroutine get_current_loading
!------------------------------------------------------------------------------
!Adaptation diffusion
subroutine get_current_loadingD(load_nb, incr, gradQDMoy, FluxDMoy, t , local_loadD0, local_loadD)
  
  implicit none
  
  integer, intent(in)                         :: load_nb, incr  ! numero du chargement et increment 
  real(mytype),dimension(3,algo_param%nVarD), intent(in)       :: gradQDMoy,FluxDMoy
  real(mytype),dimension(:), intent(in)       :: local_loadD0
  real(mytype),dimension(:), intent(out)      :: local_loadD
  real(mytype),dimension(-2:0),intent(inout)  :: t

  integer :: j,ivar,direction

  ! on initialise local_load a local_load0 (si chargement Constant, on ne fait rien)
  local_loadD = local_loadD0

  t(-2)=t(-1)         !t-2
  t(-1)=t(0)          !t-1
  t(0)=load(load_nb)%time(incr)
  local_loadD(1:3*algo_param%nVarD)=load(load_nb)%drivingD

  !valeur du chargement local
  do ivar=1,algo_param%nVarD
  do j=1,3
    direction = 3*(ivar-1)+j
    if(load(load_nb)%drivingD(j)==1)then
      !pilotage en FluxD
      call get_loading_value(load_nb, incr, direction, FluxDMoy(j,ivar), local_loadD(3*algo_param%nVarD+direction))
    else

      !pilotage en GradQ
      call get_loading_value(load_nb, incr, direction, gradQDMoy(j,ivar), &
                             local_loadD(3*algo_param%nVarD+direction))
    end if
  end do
  end do

  call get_temp_value(load_nb, incr, local_loadD(6*algo_param%nVarD+1))

  if(initValExt%nb_param>0) call get_PREDEF_value(load_nb, incr, local_loadD(6*algo_param%nVarD+2:))
end subroutine get_current_loadingD
!===================================================================================================




!===================================================================================================
!
!                           SUBROUTINE SET_TIME
!
!
!>  Calcule le temps de chaque pas de chargement (stockes dans chargement%time).
!!  
!!  la loi de discretisation est donnee par le champ "chargement%Discretization"
!!       - 0 : discretisation definie par l'utilisateur (dans le fichier xml de chargement)
!!       - 1 : discretisation lineaire
!!             -#  borne inferieure : pour  le premier chargement :0 
!!                                 sinon le temps final du chargement precedent
!!             -#  borne superieure: "Tfinal" du fichier xml stocke dans "chargement%time(nincr)"
!!  Dans les deux cas, on introduit un pas de temps intermediaire t(i)
!!  au debut de chaque nouveau chargement tel que
!!  t(i) = t(i-1) + ratio*(t(i+1)- t(i-1))
!!  ou les valeurs t(i-1) et t(i+1) sont donnes dans la definition du chargement.
!===================================================================================================
subroutine set_time()

  implicit none

  integer i, j, nincr
  real(mytype) :: t, dt
  character(len=200) :: err_msg
  
  t= 0._mytype
  
  do i=1,size(load0)
     
     nincr=load0(i)%NIncr
     if(load0(i)%Discretization==1)then
        !discretisation lineaire du temps
        dt= (load0(i)%time(nincr)-t )/nincr
        do j=1,nincr-1
           load0(i)%time(j)= t+j*dt
        end do
        t= load0(i)%time(nincr)
     else if(load0(i)%Discretization/=0)then
        ! 0 : discretisation definie par l'utilisateur
        write(err_msg,fmt="(A,I0)") "Loi de discretisation du temps non implementee : ",&
             load0(i)%Discretization," (set_time)"
        call amitex_abort(err_msg,1,0)
     end if
     !! Calcul du temps 0 : t(i-1) + ratio*(t(i+1)- t(i-1))
     if(i==1) then
        !! t(0) = 0
        load0(i)%time(0) = load0(i)%time(1)*ratio
     else
        load0(i)%time(0) = (1._mytype - ratio) * load0(i-1)%time(load0(i-1)%Nincr)&
             + load0(i)%time(1)*ratio
     end if
     t=load0(i)%time(nincr)     
  end do

end subroutine set_time

!===================================================================================================
!
!                       SUBROUTINE GET_LOADING_VALUE
!
!> Calcul du chargement courant
!! Permet de calculer une composante du chargement courant
!! en fonction du chargement precedent et des moyennes des contraintes et deformations.
!! On introduit egalement un pas de temps fictif au debut de chargement.
!! Pour ce pas de temps, l'ecart de temps avec le temps precedent est ratio*dt_n
!! ou dt_n est le pas de temps initialement impose.
!!
!! \param[in]  load_nb  (entier) numero du chargement
!! \param[in]  incr  (entier) numero de l'increment de temps (au sein du chargement)
!! \param[in]  direction  (entier) indice de composante du tenseur
!! \param[in]  moy  (reel) indice 'direction' de valeur moyenne de la contrainte ou 
!!                         de la deformation  (selon le pilotage) a l'instant precedent
!!
!! \param[in,out]   local_load  (reel) valeur de la composante direction du chargement
!!                       - a l'instant precedent (input)
!!                       - a l'instant courant (output)
!!
!===================================================================================================
subroutine get_loading_value(load_nb, incr, direction, moy, local_load)

  implicit none

  integer,intent(in) :: load_nb    ! numero du chargement
  integer,intent(in) :: incr ! numero de l'increment (dans le chargement load_nb)
  integer,intent(in) :: direction    ! indice du tenseur
  real(mytype), intent(in) :: moy           !valeur moyenne (sig/def) en fin du chargement precedent
  real(mytype), intent(inout):: local_load   !valeur du chargement a l'increment precedent
  real(mytype)               :: dt_Deltat !< ratio de pas de temps : 
  !! (t_{load_nb}^incr - t_{load_nb}^{incr+1})/(t_{load_nb}^0 - t_{load_nb}^{n})
  real(mytype)         :: value     !! valeur a imposer en fin de chargement
  integer              :: evolution !! type d'evolution
  integer              :: n, n_old  !! Nombre de pas de chargement pour ce chargement et le precedent
  character(len=200)   :: err_msg

  value = 0._mytype
  n_old = 0
  n = load(load_nb)%NIncr
  if(load_nb/=1)then
     n_old = load(load_nb-1)%NIncr
  end if
  if (algo_param%Mechanics) then
    evolution = load(load_nb)%mevolution(direction)
    value = load(load_nb)%mvalue(direction)  
  end if 
  if (algo_param%Diffusion) then
    evolution = load(load_nb)%devolution(direction)
    value = load(load_nb)%dvalue(direction)  
  end if 

  if(evolution==0) then
  !evolution constante
    local_load = moy
  else if(evolution==1) then
     !evolution lineaire
     !! Pour savoir avec quel ratio on incremente le chargement,
     !! on calcule le rapport dt/Deltat avec :
     !! - dt : pas de temps entre l'increment incr du pas de chargement
     !!        courant (load_nb) et l'increment precedent
     !! - Deltat : ecart entre le temps final du chargement courant (load_nb)
     !!            et le temps final du chargement precedent
     !! On ecrit d'abord dt_Deltat = dt
     if(incr==0)then
        !! Cas du pas 0 on va chercher le temps au pas de chargement 0
        !! du chargement courant et la valeur du temps a la fin du chargement precedent
        if(load_nb==1)then
           ! t_init = 0
           dt_Deltat = load(load_nb)%time(incr)
        else
           dt_Deltat = load(load_nb)%time(incr) - load(load_nb-1)%time(n_old)
        end if
     else        
        dt_Deltat = load(load_nb)%time(incr) - load(load_nb)%time(incr-1)
     end if
     !! Puis dt_Deltat = dt_Deltat/Deltat
     if(load_nb==1)then
        !! t_init = 0
        dt_Deltat = dt_Deltat / load(load_nb)%time(n)
     else
        dt_Deltat = dt_Deltat / (load(load_nb)%time(n) - load(load_nb-1)%time(n_old))
     end if
     if(incr==0)then
        !! Chargement particulier pour le pas 0 (depend de la valeur finale du
        !! chargement precedent)
        local_load = moy + (value-moy)*dt_Deltat
     else
        !! Sinon on evolue lineairement jusqu'a la valeur finale
        local_load = local_load + (value-moy)*dt_Deltat
        ! -> suppose que moy est initialise pour load_nb=1
     end if
  else
     write(err_msg,fmt="(A,I0,A)") "Loi d'evolution du chargement non implementee : ",&
          evolution," (get_loading_value)"
     call amitex_abort(err_msg,1,0)
  end if
  
end subroutine get_loading_value

!===================================================================================================
!
!                       SUBROUTINE GET_TEMP_VALUE
!
!> Calcul des valeurs imposees pour la temperature
!! Permet de calculer la temperature du chargement courant
!! en fonction du chargement precedent.
!! On introduit egalement un pas de temps fictif au debut de chargement.
!! Pour ce pas de temps, l'ecart de temps avec le temps precedent est ratio*dt_n
!! ou dt_n est le pas de temps initialement impose.
!!
!! \param[in]  load_nb  (entier) numero du chargement
!! \param[in]  incr  (entier) numero de l'increment de temps (au sein du chargement)
!! \param[in,out]   temp  (reel) valeur de la temperature
!!                       - a l'instant precedent (input)
!!                       - a l'instant courant (output)
!!
!===================================================================================================
subroutine get_temp_value(load_nb, incr, temp)

  implicit none

  integer,intent(in) :: load_nb    ! numero du chargement
  integer,intent(in) :: incr ! numero de l'increment (dans le chargement load_nb)
  real(mytype), intent(inout):: temp   !valeur du chargement a l'increment precedent
  real(mytype)               :: dt_Deltat !< ratio de pas de temps : 
  !! (t_{load_nb}^incr - t_{load_nb}^{incr+1})/(t_{load_nb}^0 - t_{load_nb}^{n})
  integer :: n, n_old !! Nombre de pas de chargement pour ce chargement et le precedent
  character(len=200)   :: err_msg

  n = load(load_nb)%NIncr
  n_old = 0
  if(load_nb/=1)then
     n_old = load(load_nb-1)%NIncr
  end if
  !! Si load(load_nb)%tevolution==0 on a rien a faire temp est deja a la bonne valeur
  if(load(load_nb)%tevolution==1) then
     !evolution lineaire
     !! Pour savoir avec quel ratio on incremente le chargement,
     !! on calcule le rapport dt/Deltat avec :
     !! - dt : pas de temps entre l'increment incr du pas de chargement
     !!        courant (load_nb) et l'increment precedent
     !! - Deltat : ecart entre le temps final du chargement courant (load_nb)
     !!            et le temps final du chargement precedent
     !! On ecrit d'abord dt_Deltat = dt
     if(incr==0)then
        !! Cas du pas 0 on va chercher le temps au pas de chargement 0
        !! du chargement courant et la valeur du temps a la fin du chargement precedent
        if(load_nb==1)then
           ! t_init = 0
           dt_Deltat = load(load_nb)%time(incr)
        else
           dt_Deltat = load(load_nb)%time(incr) - load(load_nb-1)%time(n_old)
        end if
     else        
        dt_Deltat = load(load_nb)%time(incr) - load(load_nb)%time(incr-1)
     end if
     !! Puis dt_Deltat = dt_Deltat/Deltat
     if(load_nb==1)then
        if(incr == 0) then
           !! t_init = 0
           dt_Deltat = dt_Deltat / load(load_nb)%time(n)
        else
           dt_Deltat = dt_Deltat / (load(load_nb)%time(n) -load(load_nb)%time(incr-1)) 
        end if
     else
        if(incr == 0) then
           dt_Deltat = dt_Deltat / (load(load_nb)%time(n) - load(load_nb-1)%time(n_old))
        else
           dt_Deltat = dt_Deltat / (load(load_nb)%time(n) -load(load_nb)%time(incr-1)) 
        end if
     end if
     temp = temp + (load(load_nb)%tvalue-temp)*dt_Deltat
  elseif(load(load_nb)%tevolution/=0) then
     write(err_msg,fmt="(A,I0,A)") "Loi d'evolution du chargement non implementee : ",&
          load(load_nb)%tevolution," (get_temp_value)"
     call amitex_abort(err_msg,1,0)
  end if
end subroutine get_temp_value

!===================================================================================================
!
!                       SUBROUTINE GET_PREDEF_VALUE
!
!> Calcul des valeurs imposees pour les parametres externes
!! Permet de calculer les parametres externes du chargement courant
!! en fonction du chargement precedent.
!! On introduit egalement un pas de temps fictif au debut de chargement.
!! Pour ce pas de temps, l'ecart de temps avec le temps precedent est ratio*dt_n
!! ou dt_n est le pas de temps initialement impose.
!!
!! \param[in]  load_nb  (entier) numero du chargement
!! \param[in]  incr  (entier) numero de l'increment de temps (au sein du chargement)
!! \param[in,out]   PREDEF (reel) valeur des parametres externes
!!                       - a l'instant precedent (input)
!!                       - a l'instant courant (output)
!!
!------------------------------------------------------------------------------
subroutine get_predef_value(load_nb, incr, PREDEF)

  implicit none

  integer,intent(in) :: load_nb    ! numero du chargement
  integer,intent(in) :: incr ! numero de l'increment (dans le chargement load_nb)
  real(mytype),dimension(:),intent(inout):: PREDEF !valeur des parametres externes
                                                   !! a l'increment precedent
  real(mytype)               :: dt_Deltat !< ratio de pas de temps : 
  !! (t_{load_nb}^incr - t_{load_nb}^{incr+1})/(t_{load_nb}^0 - t_{load_nb}^{n})
  integer :: n, n_old, i, nb_param !! Nombre de pas de chargement pour ce chargement et le precedent
                         !! i : indice pour iterer sur les parametres externes
                         !! nb_param : nombre de parametres externes
  character(len=200)   :: err_msg

  n = load(load_nb)%NIncr
  n_old = 0
  if(load_nb/=1)then
     n_old = load(load_nb-1)%NIncr
  end if
  nb_param = size(predef)

  if(nb_param /= size(load(load_nb)%pevolution)) then
     write(err_msg,fmt="(A,I0,A)") "Nombre de parametres externes mal initialise : ",&
          nb_param," (get_predef_value)"
     call amitex_abort(err_msg,2,0)
  end if

  do i=1,nb_param
     !! Si load(load_nb)%pevolution(i)==0 PREDEF(i) est deja a la bonne valeur
     if(load(load_nb)%pevolution(i)==1) then
        !evolution lineaire
        !! Pour savoir avec quel ratio on incremente le chargement,
        !! on calcule le rapport dt/Deltat avec :
        !! - dt : pas de temps entre l'increment incr du pas de chargement
        !!        courant (load_nb) et l'increment precedent
        !! - Deltat : ecart entre le temps final du chargement courant (load_nb)
        !!            et le temps final du chargement precedent
        !! On ecrit d'abord dt_Deltat = dt
        if(incr==0)then
           !! Cas du pas 0 on va chercher le temps au pas de chargement 0
           !! du chargement courant et la valeur du temps a la fin du chargement precedent
           if(load_nb==1)then
              ! t_init = 0
              dt_Deltat = load(load_nb)%time(incr)
           else
              dt_Deltat = load(load_nb)%time(incr) - load(load_nb-1)%time(n_old)
           end if
        else        
           dt_Deltat = load(load_nb)%time(incr) - load(load_nb)%time(incr-1)
        end if
        !! Puis dt_Deltat = dt_Deltat/Deltat
        if(load_nb==1)then
           if(incr == 0) then
              !! t_init = 0
              dt_Deltat = dt_Deltat / load(load_nb)%time(n)
           else
              dt_Deltat = dt_Deltat / (load(load_nb)%time(n) -load(load_nb)%time(incr-1)) 
           end if
        else
           if(incr == 0) then
              dt_Deltat = dt_Deltat / (load(load_nb)%time(n) - load(load_nb-1)%time(n_old))
           else
              dt_Deltat = dt_Deltat / (load(load_nb)%time(n) -load(load_nb)%time(incr-1)) 
           end if
        end if
        predef(i) = predef(i) + (load(load_nb)%pvalue(i)-predef(i))*dt_Deltat
     elseif(load(load_nb)%pevolution(i)/=0) then 
        write(err_msg,fmt="(A,I0,A)") "Loi d'evolution du chargement non implementee : ",&
             load(load_nb)%pevolution," (get_predef_value)"
        call amitex_abort(err_msg,1,0)
     end if
  end do
end subroutine get_predef_value
!===================================================================================================
!
!                       SUBROUTINE GET_GRADGRADU_VALUE
!
!> Calcul des valeurs imposees pour le gradgradU impose
!! Permet de calculer les composantes de gradgradU du chargement courant
!! en fonction du chargement precedent.
!! On introduit egalement un pas de temps fictif au debut de chargement.
!! Pour ce pas de temps, l'ecart de temps avec le temps precedent est ratio*dt_n
!! ou dt_n est le pas de temps initialement impose.
!!
!! \param[in]  load_nb  (entier) numero du chargement
!! \param[in]  incr  (entier) numero de l'increment de temps (au sein du chargement)
!! \param[in,out]   GRADGRADU (reel) valeur des parametres externes
!!                       - a l'instant precedent (input)
!!                       - a l'instant courant (output)
!!
!------------------------------------------------------------------------------
subroutine get_gradgradU_value(load_nb, incr, gradgradU)

  implicit none

  integer,intent(in) :: load_nb    ! numero du chargement
  integer,intent(in) :: incr ! numero de l'increment (dans le chargement load_nb)
  real(mytype),dimension(3,3,3),intent(inout):: gradgradU !valeur de gradgradU impose
                                                   !! a l'increment precedent
  real(mytype)       :: dt_Deltat !< ratio de pas de temps : 
  !! (t_{load_nb}^incr - t_{load_nb}^{incr+1})/(t_{load_nb}^0 - t_{load_nb}^{n})
  integer            :: n, n_old, i,j,k


  n = load(load_nb)%NIncr
  n_old = 0
  if(load_nb/=1)then
     n_old = load(load_nb-1)%NIncr
  end if

  do k=1,3
  do j=1,3
  do i=1,3
     !! Si load(load_nb)%pevolution(i)==0 gradgradU(i) est deja a la bonne valeur
     if(load(load_nb)%gevolution(i,j,k)==1) then
        !evolution lineaire
        !! Pour savoir avec quel ratio on incremente le chargement,
        !! on calcule le rapport dt/Deltat avec :
        !! - dt : pas de temps entre l'increment incr du pas de chargement
        !!        courant (load_nb) et l'increment precedent
        !! - Deltat : ecart entre le temps final du chargement courant (load_nb)
        !!            et le temps final du chargement precedent
        !! On ecrit d'abord dt_Deltat = dt
        if(incr==0)then
           !! Cas du pas 0 on va chercher le temps au pas de chargement 0
           !! du chargement courant et la valeur du temps a la fin du chargement precedent
           if(load_nb==1)then
              ! t_init = 0
              dt_Deltat = load(load_nb)%time(incr)
           else
              dt_Deltat = load(load_nb)%time(incr) - load(load_nb-1)%time(n_old)
           end if
        else        
           dt_Deltat = load(load_nb)%time(incr) - load(load_nb)%time(incr-1)
        end if
        !! Puis dt_Deltat = dt_Deltat/Deltat
        if(load_nb==1)then
           if(incr == 0) then
              !! t_init = 0
              dt_Deltat = dt_Deltat / load(load_nb)%time(n)
           else
              dt_Deltat = dt_Deltat / (load(load_nb)%time(n) -load(load_nb)%time(incr-1)) 
           end if
        else
           if(incr == 0) then
              dt_Deltat = dt_Deltat / (load(load_nb)%time(n) - load(load_nb-1)%time(n_old))
           else
              dt_Deltat = dt_Deltat / (load(load_nb)%time(n) -load(load_nb)%time(incr-1)) 
           end if
        end if

        gradgradU(i,j,k) = gradgradU(i,j,k) + (load(load_nb)%gvalue(i,j,k)-gradgradU(i,j,k))*dt_Deltat
     end if
     
  end do
  end do
  end do
  
end subroutine get_gradgradU_value

!===================================================================================================
!
!         SUBROUTINE DEALLOCATE_CHARGEMENT
!
!>DESALLOCATION DU TABLEAU DE CHARGEMEMENT
!
!===================================================================================================
subroutine deallocate_load(load)

  implicit none

  type(LOADING),allocatable,dimension(:), intent(inout) :: load

!  integer::i
!INUTILE
!  do i=1,size(load1)
!    if (allocated(load1(i)%time)) deallocate(load1(i)%time)
!    if (allocated(load0(i)%driving)) deallocate(load0(i)%driving)
!    if (allocated(load0(i)%drivingD)) deallocate(load0(i)%drivingD)
!    if (allocated(load0(i)%mvalue)) deallocate(load0(i)%mvalue)
!    if (allocated(load0(i)%mevolution)) deallocate(load0(i)%mevolution)
!    if (allocated(load0(i)%pvalue)) deallocate(load0(i)%pvalue)
!    if (allocated(load0(i)%pevolution)) deallocate(load0(i)%pevolution)
!    if (allocated(load0(i)%dvalue)) deallocate(load0(i)%dvalue)
!    if (allocated(load0(i)%devolution)) deallocate(load0(i)%devolution)
!  end do
  deallocate(load)

end subroutine deallocate_load

subroutine deallocate_initValExt(initValExt)

  implicit none
  type(INITLOADEXT),intent(inout) :: initValExt

  if (allocated(initValExt%param_values)) deallocate(initValExt%param_values)

end subroutine deallocate_initValExt

subroutine deallocate_local_load(local_load)

  implicit none
  type(LOCALLOAD),intent(inout) :: local_load

  if (allocated(local_load%t0)) deallocate(local_load%t0)
  if (allocated(local_load%t1)) deallocate(local_load%t1)
  if (allocated(local_load%dt)) deallocate(local_load%dt)

end subroutine deallocate_local_load


!===================================================================================================
!
!        SUBROUTINE DEALLOCATE_PARAM_EXTRACT
!
!>DESALLOCATION DU TABLEAU DES PARAMETRES D'EXTRACTION
!
!===================================================================================================
subroutine deallocate_param_extract0()
  implicit none

  integer::i
  if (algo_param%Mechanics)then
  do i=1,size(extract0%varZSTD)
    if(allocated(extract0%varZSTD(i)%val)) deallocate(extract0%varZSTD(i)%val)
  end do
  if(allocated(extract0%varZSTD)) deallocate(extract0%varZSTD)
  do i=1,size(extract0%VarVTK)
    if(allocated(extract0%varVTK(i)%val)) deallocate(extract0%varVTK(i)%val)
  end do
  if(allocated(extract0%varVTK)) deallocate(extract0%VarVTK)
  end if
  if(allocated(extract0%tpsVTK)) deallocate(extract0%tpsVTK)
  if(allocated(extract0%tpsZSTD)) deallocate(extract0%tpsZSTD)
  if(allocated(extract0%tpsSTD)) deallocate(extract0%tpsSTD)
  if(allocated(extract0%printMeanZ)) deallocate(extract0%printMeanZ)
  if(allocated(extract0%numberZstd)) deallocate(extract0%numberZstd)
  if(allocated(extract0%numberStd)) deallocate(extract0%numberStd)
end subroutine deallocate_param_extract0


end module loading_mod
