!123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!===================================================================================================
!!
!>  INITIALISATION DES FICHIERS DE SORTIE:  INITSORTIE_STD, SORTIE_STD
!!
!===================================================================================================
!!                                                   DESCRIPTION GENERALE
!!-----------------------------------------------------------------------
!!
!! LES FONCTIONS: initSortie_std, sortie_std, calcul_somme_zones, init_buffers, remplir_buffer,
!!                print_mzstd et sortie_std_macro
!!
!! Le fichier de sortie standard contient des donnees choisies par l'utilisateur
!!   l'entete est ecrite dans la fonction "initSortie_std" de ce fichier
!!   les donnees sont ecrites a chaque pas de temps et definies dans la
!!   fonction "sortie_std" de ce fichier
!!
!!
!! L'evaluation des grandeurs par zone et materiau se fait par l'intermediaire d'un buffer
!!   dont l'organisation est decrite dans  : remplir_buffer
!!
!!
!! REMARQUE :  On force les differentes tailles de buffer a etre < 2^31 car on est limite par la 
!!  |INT32|    communication MPI (voir MPI_reduce dans print_mzstd). Il s'ensuit que toutes les variables entieres relatives a ces
!!  |  vs |    tailles (comme le nombre de zones presentes dans le buffer, la taille d'une zone
!!  |INT64|    ou encore les indices de remplissage) sont de type INT32. 
!!             Si un jour on peut communiquer via MPI des tableaux de taille > 2^31 alors on
!!             pourra repasser tous ces entiers en INT64 et ne pas contraindre les tailles.
!!
!!
!===================================================================================================

module sortie_std_mod

  use ISO_FORTRAN_ENV

  use MPI
#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype, real_type
  use decomp_2d_mpi,       only : nrank
#else
  use decomp_2d,      only : mytype,real_type, nrank
#endif  
  use material_mod,   only : nmateriaux,nmateriaux0,&
                             MattotP,MattotP0, &
                             CalcRotMatrices
  use loading_mod,    only : extract,extract0
  use param_algo_mod, only : algo_param,algo_param0
  use error_mod,      only : amitex_abort
  use field_mod,      only : field_mean, var_to_std, grid, amx_size
  use amitex_mod,     only : Sig, Def, PK1, FluxD, GradQD, &
                             fic_vtk
  use amitex_mpi_mod, only: amitex_mpi_comm

  implicit none

  private 

  private :: calcul_somme_zone, calcul_somme_composite, remplir_buffer_mecanique,&
             sortie_std_macro, remplir_buffer_diffusion, print_mzstd, detF

  private :: eps

  !> Pointer "publiques" vers composantes de SIMU_AMITEX
  public :: buffer,  Matcomp_sortie,times_s

  !> Variables "publiques" utilisees lors de l'initialisation
  public :: buffer0, Matcomp_sortie0

  !> Variables publiques a transformer en pointer (ou variables 'provisoires')
  ! RAS

  !> Types publiques (pour definition de SIMU_AMITEX)
  public :: BUFFER_STD, STRUCT_MAT_SORTIE,STRUCT_ZONE_SORTIE,TIMES_STD

  !> Fonctions publiques
  public  :: initSortie_std, sortie_std, init_buffers, init_matcomp_sortie, &
             desallocation_sortie_std
             
             



  ! besoin de ce parametre pour mettre a zero les variances trop proches du zero machine
  ! et qui soulevent des erreurs dans certains cas (homogene par exemple)
  double precision,parameter                     :: eps=1e-6_8 

  ! buffer pour limiter le nombre de communications a faire pour pouvoir reunir les donnees 
  ! a imprimer dans les fichiers .std (tres utile si le nombre de zones est grand)
  !  double precision, allocatable, dimension(:)    :: buffer_z, buffer_z_tmp
  !  double precision, allocatable, dimension(:)    :: buffer_m, buffer_m_tmp
  !  double precision, allocatable, dimension(:)    :: buffer, buffer_tmp
  ! nouvelle structure BUFFER_STD (buffer et buffer_tmp sont dans buffer%t et buffer%t_tmp
  type BUFFER_STD
     double precision, allocatable, dimension(:)    :: z, z_tmp
     double precision, allocatable, dimension(:)    :: m, m_tmp
     double precision, allocatable, dimension(:)    :: t, t_tmp
  end type BUFFER_STD
  type(BUFFER_STD)          :: buffer0  ! utilisee pour lors des initialisations
  type(BUFFER_STD), pointer :: buffer
  !------------------------------------------------------------------------------------------------
  ! On cree un type qui va nous permettre de stocker les positions des voxels composites dans
  ! lesquels sont presentes chaque zone de chaque materiau. Cela nous permet de prendre en compte
  ! les contributions des voxels composites dans les differents calculs de moyennes et ecarts-types
  ! le prototype est d'avoir :
  !     struct_mat_sortie(materiau i)%zone(zone j)%pos 
  !     = l'ensemble des voxels composites du pinceau dans lesquels est presente la zone j 
  !       du materiau i, ou i et j sont les numeros globaux.
  !
  ! On cree un type derive d'un type derive car les tailles des tableaux etant variables, on ne 
  ! peut pas en faire des matrices.
  type STRUCT_ZONE_SORTIE
    ! Nombre de voxels composites contenant la zone j du materiau i 
    integer                                                   :: nombre_occurences = 0
    ! Tableau contenant les numeros (locaux) des materiaux composites associes a ces occurences,
    ! de taille nombre_occurences
    integer,allocatable,dimension(:)                          :: numM_comp
    ! Tableau contenant les numeros (locaux) des zones associees a ces occurences,
    ! de taille nombre_occurences
    integer(kind=8),allocatable,dimension(:)                  :: numZ_comp
    ! Tableau contenant les numeros (locaux) des phases dans le materiau composite,
    ! de taille nombre_occurences
    integer,allocatable,dimension(:)                          :: numPhase_comp
    ! Indice servant a remplir les tableaux numM_comp et numZ_comp
    integer                                                   :: indice_remplissage = 1
  end type STRUCT_ZONE_SORTIE

  type STRUCT_MAT_SORTIE
    ! Un tableau dont la taille est le nombre total de zones du materiau i
    ! on aura autant de 'struct_mat_sortie' que de materiaux, sur le meme modele que MattotP
    type(Struct_zone_sortie),allocatable,dimension(:)         :: zone 
  end type STRUCT_MAT_SORTIE

  type TIMES_STD
     double precision       :: wtot    = 0  ! temps total cumule ecritures std
  end type TIMES_STD


  type(Struct_mat_sortie), allocatable, dimension(:)          :: Matcomp_sortie0
  type(Struct_mat_sortie), pointer, dimension(:)              :: Matcomp_sortie
  type(TIMES_STD),pointer                                     :: times_s
  !------------------------------------------------------------------------------------------------

contains


!=======================================================================
!           FONCTION INITSORTIE_STD
!-----------------------------------------------------------------------
!
!>Permet d'ecrire l'en-tete des fichiers de sortie standards
!!
!!           Cet en-tete est a ecrire en utilisant la fonction " write(Fstd,"FMT") "
!!           ou "FMT" est le format des donnees a ecrire.
!! \param[in]       Fstd: (entier)  unite logique du fichier sortie standard (fichier d'information)
!! \param[in]       numM: (entier) si > 0 : numero de materiau correspondant au fichier (moyenne par zone)
!!                                 si < 0 : - nombre de materiaux dans le domaine (moyenne par materiau)
!!                                 si = 0 : on ecrit les moyennes globales
!!
!!
!!
!!-----------------------------------------------------------------------
!!                               IMPORTANT \n 
!!Lors de l'appel de "write" il faut utiliser l'unite logique "Fstd". \n 
!!DANS LE CAS CONTRAIRE: on ne peut pas garantir que l'entete soit ecrit
!!                       dans le bon fichier 
!=======================================================================
subroutine initSortie_std(Fstd,numM)


  implicit none
  !-----------------------------------------------------------------------
  !variables non modifiables
  integer, intent(in)           :: Fstd, numM
  !-----------------------------------------------------------------------

  integer(kind = 8)             :: p,p1,p2,p3,p4,nb_colonnes,i,nb_lignes,ivar
                                        ! nb_colonnes : permet de suivre la description des colonnes dans l'entete
                                        ! nb_lignes : nombre de lignes de 0 a mettre apres l'entete
  real(mytype),dimension(max(6*algo_param0%nVarD+1,31))   :: tab_null  
                                        ! tableau de valeur nulle pour ecrire la ligne de 0
  character(len=8)              :: p1_char,p2_char,p3_char,p4_char 
                                        ! nombre de colonnes dans les fichiers => utile pour la construction des "format" d'ecriture
                                        ! mecanique : p1 et p2 (valeurs moyennes et ecarts types)
                                        ! diffusion : p3 et p4 (valeurs moyennes et ecarts types)   
                                        ! astuce : on ne rajoute pas d'espace (1X) pour l'ecriture des ecarts-types (positifs) 
                                     
  character(len=100)            :: format


  !! tableau de valeurs nulles
  tab_null = 0._mytype 

  !! identification des nombres de colonnes
     if (algo_param0%HPP) then
       p1 = 13
       p2 = 12
     else
       p1 = 31
       p2 = 30
     end if
     p = 6*algo_param0%nVarD
     p3 = p+1
     p4 = p
     write(p1_char,"(I4)") p1 
     write(p2_char,"(I4)") p2 
     write(p3_char,"(I4)") p3 
     write(p4_char,"(I4)") p4 


  if(algo_param0%Mechanics .and. algo_param0%Diffusion) then  !!-----MECANIQUE ET DIFFUSION A TRAITER
     write(Fstd,"(A)")"# sortie_STD : MECANIQUE ET DIFFUSION NON ENCORE PRIS EN COMPTE"
  end if

  
  if(numM>0) then
     !entete du fichier contenant les moyennes par zones pour le materiau numM
     write(Fstd,"(A)")"# Le fichier est ecrit avec les zones a la suite des autres"
     write(Fstd,"(A)")"# Exemple en gnuplot pour tracer la contrainte moyenne xx en fonction&
          & de la deformation moyenne xx dans la zone i"
     if(algo_param0%HPP) then
        write(Fstd,"(A)")"# plot ""fichier.zstd"" every nbZone::i-1 u 8:2"
     else
        write(Fstd,"(A)")"# plot ""fichier.zstd"" every nbZone::i-1 u 17:2"
     end if
     write(Fstd,"(A)")"# ou nbZone est le nombre de zones dans le materiau"
     !! Si on calcule les moyennes par zone on n'ecrit pas de lignes de 0 a cause des variables internes qui ne valent pas toujours 0
     nb_lignes = 0
  elseif(numM<0) then
     !entete du fichier contenant les moyennes par materiau
     !! -numM est le nombre de materiaux et on ecrit une ligne de 0 par materiau
     nb_lignes = -numM 
     write(Fstd,"(A)")"# Le fichier est ecrit avec les materiaux a la suite des autres"
     write(Fstd,"(A)")"# Exemple en gnuplot pour tracer la contrainte moyenne xx en fonction&
          & de la deformation moyenne xx dans le materiau i"
     if(algo_param0%HPP) then
        write(Fstd,"(A)")"# plot ""fichier.mstd"" every nbMat::i-1 u 8:2"
     else
        write(Fstd,"(A)")"# plot ""fichier.mstd"" every nbMat::i-1 u 17:2"
     end if
     write(Fstd,"(A)")"# ou nbMat est le nombre de materiaux dans le domaine"
  else
     !! Si numM = 0, on ecrit les moyennes globales et une ligne de 0
     nb_lignes = 1
  end if


  !! Debut de la description des colonnes
  write(Fstd,"(A)")   "#   1ere colonne : le temps"
  nb_colonnes=1
  if(algo_param0%Mechanics) then  !!------------------MECANIQUE
!  write(Fstd,"(A)")   "#   Quantites MECANIQUE"
    if(algo_param0%HPP) then
       write(Fstd,"(A)")"#   2-7e colonne : la contrainte moyenne en xx,yy,zz,xy,xz,yz"
       write(Fstd,"(A)")"#  8-13e colonne : la deformation moyenne en xx,yy,zz,xy,xz,yz"
       write(Fstd,"(A)")"# 14-19e colonne : l'ecart type de la contrainte moyenne en &
          &xx,yy,zz,xy,xz,yz"
       write(Fstd,"(A)")"# 20-25e colonne : l'ecart type de la deformation moyenne en &
          &xx,yy,zz,xy,xz,yz"
       nb_colonnes = 26
    else
       write(Fstd,"(A)")"#   2-7e colonne : la contrainte de Cauchy moyenne en &
          &xx,yy,zz,xy,xz,yz"
       write(Fstd,"(A)")"#  8-16e colonne : la contrainte de Piola-Kirchhoff moyenne en &
          &xx,yy,zz,xy,xz,yz,yx,zx,zy"
       write(Fstd,"(A)")"# 17-22e colonne : la deformation de Green-Lagrange moyenne en &
          &xx,yy,zz,xy,xz,yz"
       write(Fstd,"(A)")"# 23-31e colonne : le gradient de u moyen en &
          &xx,yy,zz,xy,xz,yz,yx,zx,zy"
       write(Fstd,"(A)")"# 32-37e colonne : l'ecart type de la contrainte de Cauchy  &
          &xx,yy,zz,xy,xz,yz"
       write(Fstd,"(A)")"# 38-46e colonne : l'ecart type de la contrainte de Piola-Kirchhoff  &
          &xx,yy,zz,xy,xz,yz,yx,zx,zy"
       write(Fstd,"(A)")"# 47-52e colonne : l'ecart type de la deformation de Green-Lagrange  &
          &xx,yy,zz,xy,xz,yz"
       write(Fstd,"(A)")"# 53-61e colonne : l'ecart type du gradient de u  &
          &xx,yy,zz,xy,xz,yz,yx,zx,zy"
       nb_colonnes = 62
    end if
  

    !! Pour l'instant on ecrit les moyennes des variable internes uniquement par zone
    if(numM>0) then
       !! On parcourt les variables internes        
       do p=1,size(extract0%varZSTD(numM)%val,kind=8)
          ! si on doit extraire la variable interne
          if(extract0%varZSTD(numM)%val(p)) then
             write(Fstd,"(A,I0,A,I0)")"# ",nb_colonnes,"e colonne : la valeur moyenne de la variable interne ",p
             write(Fstd,"(A,I0,A,I0)")"# ",nb_colonnes+1,"e colonne : l'ecart type de la variable interne ",p
             nb_colonnes = nb_colonnes +2
          end if
       end do
    end if

  else if (algo_param0%Diffusion) then   !!------------DIFFUSION
  write(Fstd,"(A)")   "#   Quantites DIFFUSION"
     do ivar=1,algo_param0%nVarD 
       p=nb_colonnes+(ivar-1)*12 + 1
       write(Fstd,"(A,I0,A,I0,A,I0,A)")"#   ",p,"-",p+2,"e colonne : le flux moyen (",ivar,") en x,y,z"
       write(Fstd,"(A,I0,A,I0,A,I0,A)")"#   ",p+3,"-",p+5,"e colonne : le gradient moyen (",ivar,") en x,y,z"
       write(Fstd,"(A,I0,A,I0,A,I0,A)")"#   ",p+6,"-",p+8,"e colonne : l'ecart-type du flux (",ivar,") en x,y,z"
       write(Fstd,"(A,I0,A,I0,A,I0,A)")"#   ",p+9,"-",p+11,"e colonne : l'ecart-type du gradient (",ivar,") en x,y,z"
     end do
     nb_colonnes=nb_colonnes+12*(algo_param0%nVarD)
  end if
  
  ! on ajoute le nombre d'iteration au fichier des moyennes globales (.std)
  if (numM==0) then
     write(Fstd,"(A,I0,A)")"#    ",nb_colonnes,"e colonne : nombre d'iterations necessaires pour ce pas de &
          &chargement"
  end if

  !! Ecritures des lignes de zeros initiales (si besoin : pas pour les fichiers zstd)
  !!        signification de 1X : on ajoute un espace, necessaire pour les quantites potistives ou negatives
  !!                              pas necessaire pour les ecarts-types (>0) => espace = la place du signe -
  do i=1,nb_lignes
  if (algo_param0%Mechanics .and. (.not. algo_param0%Diffusion)) then  !-----MECA PURE
     if (numM==0) then
       format="("//trim(p1_char)//"(E15.8,1X),"//trim(p2_char)//"E15.8,I4)"
       write(Fstd,fmt=trim(format)) tab_null(1:p1),tab_null(1:p2),0
     else 
       format="("//trim(p1_char)//"(E15.8,1X),"//trim(p2_char)//"E15.8)"
       write(Fstd,fmt=trim(format)) tab_null(1:p1),tab_null(1:p2)
     end if
  end if
  if (algo_param0%Diffusion .and. (.not. algo_param0%Mechanics)) then  !-----DIFFUSION PURE
     if (numM==0) then
       format="("//trim(p3_char)//"(E15.8,1X),"//trim(p4_char)//"E15.8,I4)"
       write(Fstd,fmt=trim(format)) tab_null(1:p3),tab_null(1:p4),0
     else 
       format="("//trim(p3_char)//"(E15.8,1X),"//trim(p4_char)//"E15.8)"
       write(Fstd,fmt=trim(format)) tab_null(1:p3),tab_null(1:p4)
     end if
  end if
  if (algo_param0%Diffusion .and. algo_param0%Mechanics) then          !-----MECANIQUE & DIFFUSION
     if (numM==0) then
       format="("//trim(p1_char)//"(E15.8,1X),"//trim(p2_char)//"E15.8,"//trim(p3_char)//"(E15.8,1X),"//trim(p4_char)//"E15.8,I4)"
       write(Fstd,fmt=trim(format)) tab_null(1:p1),tab_null(1:p2),tab_null(1:p3),tab_null(1:p4),0
     else 
       format="("//trim(p1_char)//"(E15.8,1X),"//trim(p2_char)//"E15.8,"//trim(p3_char)//"(E15.8,1X),"//trim(p4_char)//"E15.8)"
       write(Fstd,fmt=trim(format)) tab_null(1:p1),tab_null(1:p2),tab_null(1:p3),tab_null(1:p4)
     end if
  end if
  end do !fin boucle sur nombre de lignes


end subroutine initSortie_std



!=============================================================================
!      SUBROUTINE INIT_MATCOMP_SORTIE
!-----------------------------------------------------------------------------
!> Initialise et remplis la structure Matcomp_sortie0 qui permet de prendre en compte la
!! contribution des voxels composites lors de l'ecriture des fichiers .std.
!! Est appelee lors de la phase d'initialisation du code AMITEX et Matcomp_sortie sera utilisee
!! dans la fonction sortie_std.
!!
!! Cela nous permet d'eviter de parcourir tous les voxels composites pour chaque zone de chaque
!! materiau pour le calcul des moyennes par zone.
!=============================================================================

subroutine init_matcomp_sortie

  implicit none

                                      !! Declarations !!

  ! pour les erreurs d'allocations, erreur MPI
  integer                                                  :: alloc_stat=0, ierror

  ! nombre et numero des phases constituant le materiau composite en question 
  integer                                                  :: nbPhases, phase
  ! nombre de voxels (de zones) qui constituent ledit materiau composite
  integer(kind=8)                                          :: nb_vox_comp

  ! Pour determiner le nombre total de zones des materiaux
  integer(kind=8),allocatable,dimension(:)                 :: nZones_global, nZones
  integer(kind=8)                                          :: zone_globale
  integer                                                  :: mat_global

  ! indices de boucle
  integer                                                  :: i
  integer(kind=8)                                          :: j

!-------------------------------------------------------------------------------------
  allocate(nZones_global(1));allocate(nZones(1));deallocate(nZones_global);deallocate(nZones) !avoid gcc-warning

  ! On commence par allouer les tableaux de structure struct_mat_sortie et struct_zone_sortie
  ! Pour cela on n'a besoin que des nombres de materiaux et de zones par materiaux

  allocate(Matcomp_sortie0(nmateriaux0%n_non_composites), stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort(&
                     "Espace memoire disponible insuffisant (Matcomp_sortie0, sortie_std_mod)",2)
  allocate(nZones(nmateriaux0%n_non_composites), stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort(&
                     "Espace memoire disponible insuffisant (nZones, sortie_std_mod)",2)
  nZones = 0
  allocate(nZones_global(nmateriaux0%n_non_composites), stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort(&
                     "Espace memoire disponible insuffisant (nZones_tmp, sortie_std_mod)",2)
  nZones_global = 0 

  ! Recherche du nombre total de zones

  do i=1,nmateriaux0%n_non_composites_loc
     nZones(MattotP0(i)%numM) = size(MattotP0(i)%zone(:,1),kind=8) 
     nZones(MattotP0(i)%numM) = MattotP0(i)%zone(nZones(MattotP0(i)%numM),2)
  end do
  call MPI_BARRIER(amitex_mpi_comm,ierror) 
  call MPI_Allreduce(nZones,nZones_global,nmateriaux0%n_non_composites,MPI_INTEGER8,MPI_MAX,amitex_mpi_comm,ierror)
  ! Ici, nZones(i) est le nombre total de zones du materiau i dans la cellule
  ! Ceci est vrai egalement pour les materiau interphase, et prend en compte egalement 
  ! les zones interphase. Ces deux derniers cas sont en effet couvert par construction du tableau %Zone
  ! en presence d'un materiau interphase ou de zones interphase

  do i=1,nmateriaux0%n_non_composites
     allocate(Matcomp_sortie0(i)%zone(nZones_global(i)), stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort(&
                  "Espace memoire disponible insuffisant (Matcomp_sortie0%zone, sortie_std_mod)",2)
  end do

  ! Maintenant on va compter les occurences de chaque zone de chaque materiau dans les voxels
  ! composites pour allouer le tableau pos.
  ! On parcours les voxels composites et si on trouve le materiau i, zone j, on ajoute 1 a 
  ! struct_mat_sortie(materiau i)%zone(zone j)%nombre_occurences
  do i = nmateriaux0%n_non_composites_loc+1 , nmateriaux0%n_loc
     nbPhases = MattotP0(i)%NPhase
     ! nombre de zones du materiau composite i dans le pinceau
     nb_vox_comp = size(MattotP0(i)%numZ_composite(:,1),kind=8) 
     ! On parcourt le tableau numZ_composite
     do phase = 1, nbPhases
        do j = 1, nb_vox_comp
           zone_globale = MattotP0(i)%numZ_composite(j,phase)
           mat_global = MattotP0(i)%numM_composite(phase)
           Matcomp_sortie0(mat_global)%zone(zone_globale)%nombre_occurences = &
           Matcomp_sortie0(mat_global)%zone(zone_globale)%nombre_occurences + 1
        end do
     end do
  end do

  ! A ce stade, on a les bons nombres d'occurences, on peut allouer les tableaux numX_comp et les remplir
  do i = 1 , nmateriaux0%n_non_composites
     do j = 1, nZones_global(i)
        allocate(Matcomp_sortie0(i)%zone(j)%numM_comp(Matcomp_sortie0(i)%zone(j)%nombre_occurences)&
                                                                             ,stat=alloc_stat)
        if(alloc_stat /=0) call amitex_abort(&
                  "Espace memoire disponible insuffisant (Matcomp_sortie0%pos, sortie_std_mod)",2)
        allocate(Matcomp_sortie0(i)%zone(j)%numZ_comp(Matcomp_sortie0(i)%zone(j)%nombre_occurences)&
                                                                             ,stat=alloc_stat)
        if(alloc_stat /=0) call amitex_abort(&
                  "Espace memoire disponible insuffisant (Matcomp_sortie0%pos, sortie_std_mod)",2)
        allocate(Matcomp_sortie0(i)%zone(j)%numPhase_comp(Matcomp_sortie0(i)%zone(j)%nombre_occurences)&
                                                                             ,stat=alloc_stat)
        if(alloc_stat /=0) call amitex_abort(&
                  "Espace memoire disponible insuffisant (Matcomp_sortie0%pos, sortie_std_mod)",2)
     end do
  end do

  ! Remplissage, on re-parcours les voxels composites
  do i = nmateriaux0%n_non_composites_loc+1 , nmateriaux0%n_loc
     nbPhases = MattotP0(i)%NPhase
     nb_vox_comp = size(MattotP0(i)%numZ_composite(:,1),kind=8)
     do phase = 1, nbPhases
        do j = 1, nb_vox_comp
           zone_globale = MattotP0(i)%numZ_composite(j,phase)
           mat_global = MattotP0(i)%numM_composite(phase)
           Matcomp_sortie0(mat_global)%zone(zone_globale)%numZ_comp &
(Matcomp_sortie0(mat_global)%zone(zone_globale)%indice_remplissage) = j
           Matcomp_sortie0(mat_global)%zone(zone_globale)%numM_comp &
(Matcomp_sortie0(mat_global)%zone(zone_globale)%indice_remplissage) = i
           Matcomp_sortie0(mat_global)%zone(zone_globale)%numPhase_comp &
(Matcomp_sortie0(mat_global)%zone(zone_globale)%indice_remplissage) = phase
           Matcomp_sortie0(mat_global)%zone(zone_globale)%indice_remplissage = &
           Matcomp_sortie0(mat_global)%zone(zone_globale)%indice_remplissage + 1
        end do
     end do
  end do

  deallocate(nZones)
  deallocate(nZones_global)

end subroutine init_matcomp_sortie


!!===============================================================================================
!! ALLOCATION DES BUFFERS POUR SORTIE_STD
!------------------------------------------------------------------------------------------------
!! On fixe une taille a ne pas depasser, qu'on appellera taille_memoire_max, qu'on fixe relativement
!! a la taille d'un champ sur un pinceau (amx_size(1)*amx_size(2)*amx_size(3))  
!!
!! On regarde ensuite la plus grande taille de buffer dont on aura besoin
!! (qui correspond au plus grand nombre de zones que multiplie le nombre d'informations par zone) et
!! on alloue le buffer a la plus petite de ces deux tailles pour tout le calcul.
!!
!! On limite egalement la taille max du buffer%z a 2^31 - 1 
!!    -> assure que la taille du buffer est bien un integer (kind=4) pour l'utilisation de MPI_reduce 
!!       dans print_zstd
!!
!! Pour buffer et buffer%m, on suppose qu'on aura pas de probleme de taille et on les alloue
!! directement a la taille necessaire.
!!
!! En sortie, les buffers buffer(%t,%z,%m) et buffer(%t,%z,%m)_tmp sont alloues  
!!
!================================================================================================


subroutine init_buffers(coeff_buffer_zstd)

  implicit none


  integer                       :: ierror       ! erreur MPI
  integer                       :: alloc_stat=0 ! erreur alloc

  ! Donnees relatives a la construction du buffer
 ! taille max que l'on se fixe pour le buffer
  integer(kind=8)               :: taille_memoire_max
 ! la taille qu'on alloue finalement au buffer
  integer                       :: taille_du_buffer
 ! la taille d'une zone dans le buffer
  integer                       :: taille_zone
 ! le plus grand nombre de zones des materiaux
  integer(kind=8)               :: nZones_max, nZones_tmp
  integer(kind=8)               :: cell_size ! taille d'un pinceau
 ! Variables internes
  integer                       :: nVI_to_print, nVI_to_print_max, p, i

   !< coefficient multiplicatif pour regler la taille
   !  du buffer%z pour sortie_std
  real(mytype), intent(in)              :: coeff_buffer_zstd


  ! les champs sont locaux (aux pinceaux)
  ! On fixe la taille a ne pas depasser comme etant la taille du plus grand champ,
  ! ie nombre de voxels total du plus gros pinceau
  cell_size=int(amx_size(1),kind=8)*int(amx_size(2),kind=8)*int(amx_size(3),kind=8)
  call MPI_BARRIER(amitex_mpi_comm,ierror) 
  call MPI_Allreduce(cell_size,taille_memoire_max,1,MPI_INTEGER8,MPI_MAX,amitex_mpi_comm,ierror)
  taille_memoire_max = int(taille_memoire_max * coeff_buffer_zstd)


  nZones_max = 0
  nVI_to_print_max = 0
  ! On determine quel est le plus grand nombre de zones et de variables internes
  ! pour avoir la plus petite taille suffisante au buffer pour tout ranger dedans
  ! en un seul coup
  do i = 1, nmateriaux0%n_non_composites_loc
     nZones_tmp = size(MattotP0(i)%zone(:,1),kind=8) 
     nZones_tmp = MattotP0(i)%zone(nZones_tmp,2)
     nZones_max = max(nZones_max,nZones_tmp)
     if (extract0%printMeanZ(i) .AND. algo_param0%mechanics) then
        nVI_to_print = 0
        do p = 1, size(extract0%varZSTD(i)%val)
           if (extract0%varZSTD(i)%val(p)) nVI_to_print = nVI_to_print + 1
        end do
        nVI_to_print_max = max(nVI_to_print,nVI_to_print_max)
     end if
  end do

  call MPI_BARRIER(amitex_mpi_comm,ierror) 
  call MPI_Allreduce(nZones_max,nZones_tmp,1,MPI_INTEGER8,MPI_MAX,amitex_mpi_comm,ierror)
  nZones_max = nZones_tmp
  ! maintenant, nZones_max contient le plus grand nombre de zones presentes dans un meme materiau
  ! y compris les zones interphase. Valable egalement pour les materiaux interphase
  ! nVI_to_print_max le plus grand nombre de variables internes a sortir dans fichier zstd

  ! On calcule la taille que prend une zone (un materiau / la cellule) dans le buffer associe
  ! (voir organisationdans remplir_buffer)
  taille_zone = 1
  if (algo_param0%Mechanics) then
     if (algo_param0%hpp) then
        taille_zone = taille_zone + 24
     else
        taille_zone = taille_zone + 61
     end if
  end if
  if (algo_param0%Diffusion) then
     taille_zone = taille_zone + 12*algo_param0%nVarD
  end if

  ! taille_zone est pour l'instant la place qu'occupe la cellule dans buffer, c'est aussi la 
  ! taille qu'occupe un materiau dans buffer%m
  ! On alloue alors buffer et buffer%m 
  if (.NOT. allocated(buffer0%m)) allocate(buffer0%m(taille_zone*nmateriaux0%n_non_composites),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort(&
              "init_buffers (sortie_std_mod) : Espace memoire disponible insuffisant (buffer0%m)",2)
  if (.NOT. allocated(buffer0%t)) allocate(buffer0%t(taille_zone),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort(&
                "init_buffers (sortie_std_mod) : Espace memoire disponible insuffisant (buffer0)",2)
  !! une copie de ces buffers pour le mpi_reduce
  if (.NOT. allocated(buffer0%m_tmp)) allocate(buffer0%m_tmp(taille_zone*nmateriaux0%n_non_composites),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort(&
              "init_buffers (sortie_std_mod) : Espace memoire disponible insuffisant (buffer0%m)",2)
  if (.NOT. allocated(buffer0%t_tmp)) allocate(buffer0%t_tmp(taille_zone),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort(&
                "init_buffers (sortie_std_mod) : Espace memoire disponible insuffisant (buffer0)",2)

  ! On ajoute les variables internes pour avoir la place que prend une zone dans buffer0%z
  taille_zone = taille_zone + 2*nVI_to_print_max

  
  ! On verifie que la taille max permet de rentrer 1000 zones dans le buffer0. Un tel buffer pese
  ! au max 2,4Mo (en considerant un calcul GD et ~ 250 variables internes) 
  ! ce qui n'est pas grand chose, et evite de devoir boucler plein de fois pour remplir
  ! le buffer dans le cas ou on a une petite cellule (et donc la taille d'un champ est faible).
  ! on impose une taille maxi inferieure a 2^31
  taille_memoire_max = max(1000*taille_zone, taille_memoire_max)
  taille_memoire_max = min((2_8**31) - 1,taille_memoire_max)


  ! On peut allouer le buffer a la bonne taille, en deux exemplaires pour le mpi_reduce.
  ! Soit le materiau le plus gros (en nombre de zones) rentre en entier dans le buffer, 
  ! auquel cas on alloue la place suffisante pour ce materiau, soit il ne rentre pas et 
  ! on alloue la taille max.
  taille_du_buffer = int(min(taille_memoire_max,taille_zone*nZones_max),kind=INT32)

  if (.NOT. allocated(buffer0%z)) allocate(buffer0%z(taille_du_buffer),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort(&
              "init_buffers (sortie_std_mod) : Espace memoire disponible insuffisant (buffer%z)",2)
  if (.NOT. allocated(buffer0%z_tmp)) allocate(buffer0%z_tmp(taille_du_buffer),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort(&
          "init_buffers (sortie_std_mod) : Espace memoire disponible insuffisant (buffer%z_tmp)",2)  


end subroutine init_buffers



!==============================================================================
!       SUBROUTINE SORTIE_STD
!------------------------------------------------------------------------------
!> Ecriture des donnees sur le fichier de sortie standard.
!! Fonction  appelee a chaque fin de pas de calcul. \n
!! Permet de calculer et d'ecrire des grandeurs d'interet dans le fichier de sortie standard
!! (par exemple : contrainte moyenne, moyennes par phase etc...).
!!
!!
!! \param[in]      t: (reel) temps du chargement
!! \param[in]      ind_tps: (entier) pas de temps courant
!! \param[in]      nb_it: (entier) nombre d'iterations
!! \param[out]     MoyStress, MoyDef,FluxDMoy,gradQDMoy :
!!                       tableaux de reels (dimensions 6 ou 9), contraintes 
!!                       et deformations moyennes utilisees comme valeurs initiales
!!                       pour le pas de chargement suivant. En HPP, cela correspond a la
!!                       contrainte de Cauchy et a la deformation. 
!!                       En grandes transformations, cela correspond a la contrainte
!!                       de Piola-Kirchhoff et au gradient du deplacement
!!
!! REMARQUE :
!!
!! Dans les buffers les "moyennes" sont les sommes des valeurs sur les voxels
!! et les "ecarts-types" sont les sommes des carres des valeurs. 
!!
!! On ne les transformera en moyenne/ecart-type a proprement
!! parler qu'en fin de fonction lors de l'ecriture dans les fichiers.
!!
!!
!! REMARQUE BUFFERS :
!!
!!    buffer (pour sortie macro) et buffer%m (pour sortie par materiau)
!!    ont les memes caracteristiques
!!    buffer%z (pour sortie par zone) contient en plus, si demande, 
!!    des variables internes 
!!
!!------------------------------------------------------------------------------
!!                 09/2016 CORRECTION BUG (PLANTAGES ALEATOIRES AVEC openMPI1.8)
!!                                      AJOUT DE MPI_BARRIER DERRIERE MPI_REDUCE
!!                                                  (mal compris mais ça marche)
!!
!!------------------------------------------------------------------------------
!!                                                                     IMPORTANT \n
!! Lors de l'appel de "write" il faut utiliser l'unite logique "Fstd". \n
!! DANS LE CAS CONTRAIRE: on ne peut pas garantir que les donnees soient ecrites
!!                       dans le bon fichier.\n
!!
!!!
!! Les differents tableaux (pinceaux en X) sont definies de "amx_start(i)" a "amx_end(i)" 
!! ou i correspond a l'axe: 
!! - 'x' si i=1
!! - 'y' si i=2
!! - 'z' si i=3
!! Les valeurs de amx_start et amx_end varient pour chaque processus.
!!Donc pour atteindre une position specifique (par exemple Sig(2,5,1,:) )
!! on prendra soin de verifier que ces coordonnees existent pour le processus
!! (un processus est defini par son rang).\n
!!
!!------------------------------------------------------------------------------
!!                                                                     IMPORTANT\n
!!Les contraintes et deformations moyennes doivent etres calculees puisqu'elles
!!sont utilisees comme valeurs initiales certains chargement
!!(le premier chargement d'un ensemble defini dans le fichier xml de chargement).
!!
!!
!!------------------------------------------------------------------------------
!!                                                                     IMPORTANT\n
!!Pour des raisons de precision sur les moyennes et les ecarts types, on utilise
!!des variables temporaires de type double precision.
!!
!! Si tel n'est pas le cas, en simple precision, on obtient une valeur de
!! E[X^2]-E[X]^2 negative et donc une valeur d'ecart type de "NaN"
!!
!==============================================================================
  

subroutine sortie_std(t,ind_tps,nb_it,MoyStress,MoyDef,FluxDMoy,gradQDMoy)
  !
  !On considere ici que les phases sont definies par les numeros de materiau
  !


  implicit none
  !------------------------------------------------------------------------------

  ! Variables d'entree
  real(mytype), intent(in)   :: t
  integer, intent(in)        :: ind_tps, nb_it

  ! Variables liees aux fichiers de sortie
  character(len=200)         :: fic_std, fic_local
  Integer                    :: FZstd,Fstd

  ! Sorties de la routine
  ! MoyStress est le tenseur de cauchy en HPP et de Piola-Kirchhoff en grandes deformations
  ! MoyDef est le tenseur des deformations en HPP et le gradient de u en grandes deformations
  real(mytype),dimension(algo_param%nTensDef),intent(out) :: MoyStress, MoyDef
  real(mytype),dimension(3,algo_param%nVarD),intent(out)  :: FluxDMoy, gradQDMoy

!------------------------------------------------------------------------------------------

  ! Test pour sortir (ou non) les resultats dans les fichiers .std et .mstd
  logical                                        :: testPrintStd  
  ! variables relatives aux buffers
  ! tailleBuff(_m / _z) represente la taille que prend la cellule (un materiau / une zone)
  ! dans le buffer (en nombre de cases du tableau)
  integer                                        :: tailleBuff_z, tailleBuff_m
  ! le nombre de fois qu'on doit remplir le buffer
  integer                                        :: nb_tranche
  ! nBuff_m(_z) est le nombre de materiaux (zones) qu'on a dans le buffer
  integer                                        :: nBuff_z, nBuff_m
  ! le plus grand nombre de zones qui rentrent dans le buffer%z
  integer                                        :: nZonesMax

  ! indices relatifs au remplissage et a la lecture des buffers
  integer                                        :: pos_remplissage_buff_m
  integer                                        :: pos_remplissage_buff_z
  integer                                        :: pos_remplissage_buff
  integer                                        :: ind_tranche, indice_diff, indice
  integer                                        :: debut_mat_i_dans_buff_m

  integer                                        :: ierror     ! erreur MPI

  !> moyenne et ecart type des contraintes et deformations par zone
 ! nombre de voxels dans la zone
  double precision                                        :: nVox_z
 ! somme des determinants de F
  double precision                                        :: Sdet_z
 ! contrainte de Cauchy (moyenne et ecart-type)
  double precision,dimension(6)                           :: MS_z, EctS_z
 ! deformation de Green-Lagrange (moyenne et ecart-type)
  double precision,dimension(6)                           :: ME_z, EctE_z
 ! flux et gradient des variables de diffusion (moyenne et ecart-type)
  double precision,dimension(3,algo_param%nVarD)          :: MS_zD,ME_zD, EctS_zD, EctE_zD
 ! contrainte de Piola-Kirchhoff et gradient de u (moyenne et ecart-type)
  double precision,dimension(merge(0,1,algo_param%HPP)*9) :: MP_z,MG_z, EctP_z,EctG_z
 ! variables internes (moyenne et ecart-type)
  double precision,allocatable,dimension(:)               :: varInt_z,EctVI_z
  integer                                                 :: nVI_to_print, p
 ! variable de test 
  logical                                                 :: test_calc
  ! indices de materiaux et de zones
  integer(kind=8)                                         :: numZ, j, nbZone
  integer                                                 :: i, numM, zone
  double precision                                        :: t1 ! time measurement


!------------------------------------------------------------------------------------------
  call MPI_BARRIER(amitex_mpi_comm,ierror)          
  t1 = MPI_WTIME()

  Sdet_z = 0. !avoid gcc-warning

  !> On recupere le test pour sortir les moyennes sur la cellule et par materiau (.std et .mstd)
  testPrintStd = extract%tpsSTD(ind_tps)

  !! Ouverture de fichier std et (ouverture,flush,fermeture) du fichier .zstd (pourquoi?) : INTERET???
  fic_std = trim(fic_vtk)//".std"
  if(nrank==0)then
     !fichier .std
     open(newunit=Fstd, file=fic_std,form="formatted", status="old",position="append", action="write")
     if(extract%printZSTD .AND. extract%tpsZSTD(ind_tps)) then
        do i=1,nmateriaux%n_non_composites
           if(extract%printMeanZ(i)) then
              write(fic_local,fmt="(A,I0,A)") trim(fic_vtk)//"_",i,".zstd"
              open(newunit=FZstd, file=trim(fic_local),form="formatted", status="old",&
                   position="append", action="write")    
              flush(FZstd)
              close(FZstd)
           end if
        end do
     end if
     close(Fstd)
  end if
  
  !! Dans le cas ou on n'a que le .std a sortir, on passe par un calcul vectoriel qui calcule plus
  !  efficacement les moyennes macroscopiques.
  if (nmateriaux%n_non_composites==1 .AND. &
          (.NOT. extract%printMeanZ(1) .OR. .NOT.extract%tpsZSTD(ind_tps))) then
     call sortie_std_macro(t,nb_it,MoyStress,MoyDef,FluxDMoy,gradQDMoy,testPrintStd)
     call MPI_BARRIER(amitex_mpi_comm,ierror)          
     times_s%wtot = times_s%wtot + MPI_WTIME() - t1
     return
  end if

  !! Decompte des tailles des buffers et initialisation des differents indices
  tailleBuff_m = 1   ! on aura toujours le nombre de voxels total composant le materiau
  tailleBuff_z = 1   ! on aura toujours le nombre de voxels total composant la zone
  nVI_to_print = 0
  if (algo_param%Mechanics) then !-------MECANIQUE
    if (algo_param%HPP) then
      ! on compte moyenne et ecart-type pour chaque composante
      tailleBuff_m = tailleBuff_m + 2*6 + 2*6
      tailleBuff_z = tailleBuff_z + 2*6 + 2*6
    else !contraintes PK (P), gradient de U (G) en Grandes Transfo. et le déterminant de la transformation
       tailleBuff_m = tailleBuff_m + 2*6 + 2*6 + 2*9 + 2*9 + 1
       tailleBuff_z = tailleBuff_z + 2*6 + 2*6 + 2*9 + 2*9 + 1
    end if
  end if
  if (algo_param%Diffusion) then !-------DIFFUSION
    tailleBuff_m = tailleBuff_m + 2*2*3*algo_param%nVarD
    tailleBuff_z = tailleBuff_z + 2*2*3*algo_param%nVarD
    ! on a la moyenne et l'ecart-type de deux vecteurs a trois composantes (2*2*3)
  end if


  pos_remplissage_buff_z = 1
  buffer%z = 0
  pos_remplissage_buff_m = 1
  buffer%m = 0
  nBuff_m = nmateriaux%n_non_composites
  buffer%t = 0
  i = 1        !! i represente l'indice local du materiau

  do numM=1,nmateriaux%n_non_composites  !================================== BOUCLE SUR LES MATERIAUX

     pos_remplissage_buff_m = 1 + (numM-1)*tailleBuff_m

     !! Recherche du nombre total de zones du materiau numM :
     !! Si numM est dans le processus alors numM = MattotP(i)%numM
     !! car i evolue de maniere croissante et suit l'indice global numM
     !! Les materiaux "Interphase", qui ne sont presents que sous forme de constituants de
     !! "voxels composites", sont maintenant inclus, ils seront calcules dans 
     !! calcul_somme_composite avec tous les voxels composites.
     numZ=0
     if (i <= nmateriaux%n_non_composites_loc) then
     if(MattotP(i)%numM==numM) then
        !! nbZone = Nombre de zones sur le pinceau
        nbZone=size(MattotP(i)%zone(:,1),kind=8) 
        !! numZ = Indice maximal des zones sur le pinceau
        numZ=MattotP(i)%zone(nbZone,2)
     else
        numZ=0
     end if
     end if
     !! nbZone = Nombre de zones total
     call MPI_BARRIER(amitex_mpi_comm,ierror) 
     call MPI_Allreduce(numZ,nbZone,1,MPI_INTEGER8,MPI_MAX,amitex_mpi_comm,ierror)

     !! j est l'indice de zone local
     j=1

     ! On compte les variables internes a afficher
     nVI_to_print = 0
     if (algo_param%Mechanics) then
       if(extract%printMeanZ(numM) .AND. extract%tpsZSTD(ind_tps)) then
          do p = 1, size(extract%varZSTD(numM)%val)
             if (extract%varZSTD(numM)%val(p)) nVI_to_print = nVI_to_print + 1
          end do
       end if
     end if
     allocate(varInt_z(nVI_to_print))
     allocate(EctVI_z(nVI_to_print))
     tailleBuff_z = tailleBuff_z + 2*nVI_to_print

     ! On a le nombre de zones, on vérifie que ça ne fait pas trop pour le buffer%z,
     ! auquel cas on partitionne les zones et on les traite "tranche par tranche"
     ! Ex : On ne peut rentrer que 100 zones dans le buffer et on en a 450, alors on 
     !      les traite 100 par 100, on a donc 5 tranches de zones (4 de 100 et 1 de 50)
     ! La condition de taille s'exprime par 
     ! 64bits*nb_Zones*taille_zone_dans_buffer  <  taille max du buffer (fixée à 500Mo ??)


     !nombre maximal de zones qu'on peut rentrer dans buffer%z
     nZonesMax = int(size(buffer%z)/tailleBuff_z)

     !nombre de fois qu'on aura a remplir le buffer
     nb_tranche = int((nbZone-1)/nZonesMax)+1


     do ind_tranche = 1, nb_tranche  !-----------------------BOUCLE SUR LES TRANCHES DE ZONES
        !!                                  a chaque tour, on remplit, ecrit et vide buffer%z
        !!-----------------------------------------------------------------------------------
        buffer%z = 0 
        pos_remplissage_buff_z = 1 !curseur de position dans buffer%z

        !! On regarde combien on met de zone dans le buffer, soit c'est la dernière
        !! tranche et on met ce qu'il reste, soit on en met le nombre max
        if(ind_tranche==nb_tranche) then 
          nBuff_z = int(mod(nbZone-1,nZonesMax),kind=INT32)+1
        else
          nBuff_z = nZonesMax
        end if

        do zone=1,nBuff_z   !------------------------------------------------BOUCLE SUR LES ZONES
           !!                     EVALUATION DES MOY. ET ECART-TYPES PAR ZONE
           !!---------------------------------------------------------------- 
           !! INITIALISATION a 0
           !! Moyennes et ecart-types des contraintes (S) et des deformations (E) 
           !! dans tout la zone numZ du materiau numM
           pos_remplissage_buff_z = 1 + (zone-1)*tailleBuff_z
           numZ = zone + (ind_tranche-1)*nZonesMax ! indice global de zone
           nVox_z = 0
           test_calc = .true.   !  par defaut on calcule toute les zones. On teste plus bas si la 
                                !  zone numZ est une zone interphase
           if(algo_param%Mechanics) then
             MS_z=0
             ME_z=0
             EctS_z=0
             EctE_z=0
             if(.NOT.algo_param%HPP) then
               Sdet_z=0
               MP_z=0
               MG_z=0
               EctP_z=0
               EctG_z=0
             end if
             if(nVI_to_print>0) then
               varInt_z=0
               EctVI_z=0
             end if
           end if
           if (algo_param%Diffusion) then
             MS_zD=0
             ME_zD=0
             EctS_zD=0
             EctE_zD=0
           end if

           !! Si le materiau numM et la zone numZ sont presents dans le pinceau
           !! Si la zone numZ est dans le processus alors numZ = MattotP(i)%zone(j,2)
           !! car j evolue de maniere croissante et suit l'indice global numZ

           if (i <= nmateriaux%n_non_composites_loc) then
           if (.NOT.MattotP(i)%Interphase) then
           if (j <= size(MattotP(i)%zone(:,1),kind=8)) then 
           if(MattotP(i)%numM==numM .AND. MattotP(i)%zone(j,2)==numZ ) then
              ! C'est ici qu'on calcule les moyennes (sommes en realite) par zone
              ! pour les voxels homogenes (de materiaux non composites)
              if (allocated(mattotP(i)%Zones_interphase)) then
                 ! le materiau contient des zones Interphases : test pour savoir si 
                 ! la zone designee par numZ en est une
                 if (any(mattotP(i)%Zones_interphase == numZ)) then
                    ! si zone interphase : test_calc = .false., pas de calculs
                    ! les sommes restent a 0
                    test_calc = .false.
                 end if
              end if
           if (test_calc) then   
              call calcul_somme_zone(i,numM,j,Sdet_z,MS_z,ME_z,EctS_z,EctE_z,MS_zD,&
                   ME_zD,EctS_zD,EctE_zD,MP_z,MG_z,EctP_z,EctG_z,nVox_z,nVI_to_print,varInt_z,EctVI_z)
           end if
              j = j + 1
           end if
           end if
           end if
           end if

           ! On calcule la contribution apportee par les portions du materiau numM, zone numZ situees
           ! dans les voxels composites du pinceau :
           ! sommes sur les voxels composites (valeurs ponderees par les fractions volumiques)
           ! modif LG 06-2019: uniquement si presence de Voxel composite
           if (nmateriaux%n_composites > 0) then
              call calcul_somme_composite(numM, numZ, nVI_to_print, MS_z, ME_z, EctS_z, EctE_z,&
                                                              nVox_z, varInt_z, EctVI_z)
           end if

           !------------ REMPLISSAGE DES BUFFER
           pos_remplissage_buff = 1
           debut_mat_i_dans_buff_m = pos_remplissage_buff_m

           !-- remplissage du nombre de voxels et +1 sur les curseurs
           !     initialement dans remplir_buffer, deplace ici apres split de la fonction 
           !     (suite a bug avec intel19 pas compris)
           buffer%t(pos_remplissage_buff) = buffer%t(pos_remplissage_buff) + nVox_z
           pos_remplissage_buff = pos_remplissage_buff +1
           buffer%m(pos_remplissage_buff_m) = buffer%m(pos_remplissage_buff_m) + nVox_z
           pos_remplissage_buff_m = pos_remplissage_buff_m +1
           buffer%z(pos_remplissage_buff_z) = buffer%z(pos_remplissage_buff_z) + nVox_z
           pos_remplissage_buff_z = pos_remplissage_buff_z +1

           !-- fonction remplir_buffer, splittee suite a bug inte19 (pas compris)
           call remplir_buffer_mecanique(pos_remplissage_buff_z,&
     pos_remplissage_buff_m, pos_remplissage_buff, Sdet_z,MS_z,ME_z,EctS_z,EctE_z,&
     MP_z,MG_z,EctP_z,EctG_z,nVI_to_print,varInt_z,EctVI_z)

           call remplir_buffer_diffusion(pos_remplissage_buff_z,&
     pos_remplissage_buff_m, pos_remplissage_buff,MS_zD,ME_zD,EctS_zD,EctE_zD)

           !-- a la fin de chaque zone, on ramene le curseur materiau a sa position 
           !   initialement dans remplir_buffer, deplace ici apres split de la fonction 
           !   (suite a bug avec intel19 pas compris)
           pos_remplissage_buff_m = debut_mat_i_dans_buff_m

        end do     !------------------------------------------------FIN BOUCLE SUR LES ZONES

        !! Le buffer%z est plein et est pret, on l'imprime dans le fichier .zstd
        !! lorsqu'il est demandé par l'utilisateur

        if (extract%printMeanZ(numM) .AND. extract%tpsZSTD(ind_tps)) then
           fic_local = ""
           write(fic_local,fmt="(A,I0,A)") trim(fic_vtk)//"_",numM,".zstd"
           call print_mzstd(buffer%z,buffer%z_tmp,nBuff_z,tailleBuff_z,nb_it,t,&
                                                    nVI_to_print,fic_local,.true.)
        end if
     end do !---------------------------------------------FIN BOUCLE SUR LES TRANCHES DE ZONES

     ! Le nombre de variables internes a afficher est propre a chaque materiau
     if (allocated(varInt_z)) deallocate(varInt_z)
     if (allocated(EctVI_z)) deallocate(EctVI_z)
     tailleBuff_z = tailleBuff_z - 2*nVI_to_print

     ! on incremente le numero de materiau local (si le matériau numM etait present dans le pinceau)
     if (i <= nmateriaux%n_non_composites_loc) then
        if(MattotP(i)%numM==numM) i=i+1         
     end if
  end do  !--------------------------------------------------FIN BOUCLE SUR LES MATERIAUX


  !! On imprime le fichier .mstd quand on a plus d'un materiau
  if (nmateriaux%n_non_composites>1) then
     fic_local = ""
     write(fic_local,fmt="(A)") trim(fic_vtk)//".mstd"
     call print_mzstd(buffer%m,buffer%m_tmp,nBuff_m,tailleBuff_m,nb_it,t,0,fic_local,testPrintStd)
  end if

  !! print du fichier .std
  fic_local = ""
  write(fic_local,fmt="(A)") trim(fic_vtk)//".std"
  call print_mzstd(buffer%t,buffer%t_tmp,1,tailleBuff_m,nb_it,t,0,fic_local,testPrintStd)

  !! on affecte les grandeurs moyennes de sortie de la procédure
  if (nrank==0) then
     indice = 2
     if(algo_param%Mechanics) then
       if(algo_param%HPP) then
          ! cas ou MoyStress est le champ des contraintes de Cauchy
          MoyStress = real(buffer%t(indice:indice+5),mytype)
          indice = indice + 12
          ! cas ou MoyDef est la deformation de Green-Lagrange
          MoyDef = real(buffer%t(indice:indice+algo_param%nTensDef-1),mytype)
          indice = indice + 12
       else
          ! cas ou MoyStress est le champ des contraintes de Piola-Kirchhoff
          MoyStress = real(buffer%t(indice+25:indice+33),mytype)
          ! cas ou MoyDef est le gradient de u
          MoyDef = real(buffer%t(indice+43:indice+51),mytype)
          indice = indice + 61
       end if
     end if
     if(algo_param%Diffusion) then
        do indice_diff = 1,algo_param%nVarD
          FluxDMoy(:,indice_diff) = real(buffer%t(indice:indice+2),mytype)
          indice = indice + 6
          GradQDMoy(:,indice_diff) = real(buffer%t(indice:indice+2),mytype)
          indice = indice + 6
        end do
     end if 
  end if

  !! On envoie a tous les processeurs les valeurs moyennes (seul 0 a les valeurs globales)
  call MPI_Bcast(MoyDef,algo_param%nTensDef,real_type,0,amitex_mpi_comm,ierror)
  call MPI_Bcast(MoyStress,algo_param%nTensDef,real_type,0,amitex_mpi_comm,ierror)
  call MPI_Bcast(GradQDMoy,3*algo_param%nVarD,real_type,0,amitex_mpi_comm,ierror)
  call MPI_Bcast(FluxDMoy,3*algo_param%nVarD,real_type,0,amitex_mpi_comm,ierror)

  call MPI_BARRIER(amitex_mpi_comm,ierror)          
  times_s%wtot = times_s%wtot + MPI_WTIME() - t1
  
end subroutine sortie_std


!==============================================================================
!       SUBROUTINE CALCUL_SOMME_ZONE
!------------------------------------------------------------------------------
!!
!> Calcul des valeurs moyennes et des ecarts-types par zone des quantites demandees
!! par l'utilisateur (contrainte, deformation, variables internes).
!! Pour une zone d'un materiau demmandee, parcours tous les voxels homogenes et
!! ajoute les contributions du voxel sur les moyennes et les ecarts-types
!!
!! Permet d'alleger la routine sortie_std
!!
!!
!! REMARQUE :
!!
!! A ce stade, les "moyennes" sont les sommes des valeurs sur les voxels et les
!! "ecarts-types" sont les sommes des carres des valeurs. On ne les transformera
!! en moyenne/ecart-type a proprement parler que lors de l'ecriture dans les fichiers.
!! 
!==============================================================================

subroutine calcul_somme_zone(numM,i,numZ,Sdet_z,MS_z,ME_z,EctS_z,EctE_z,MS_zD,&
     ME_zD,EctS_zD,EctE_zD,MP_z,MG_z,EctP_z,EctG_z,nVox_z,nVI_to_print,varInt_z,EctVI_z)


  implicit none


  !------------------------------------------------------------------------------
  ! variables non modifiables

  integer(kind=8),intent(in)              :: numZ
  integer,intent(in)                      :: numM, nVI_to_print, i 

  !------------------------------------------------------------------------------
  ! variables modifiables

  ! moyenne et ecart type des contraintes et deformation par zone :

 ! nombre de voxels dans la zone
  double precision,intent(out)                               :: nVox_z
 ! somme des determinants de F
  double precision,intent(out)                               :: Sdet_z
 ! contrainte de Cauchy
  double precision,dimension(6),intent(out)                  :: MS_z, EctS_z
 ! deformation de Green-Lagrange
  double precision,dimension(6), intent(out)                 :: ME_z, EctE_z
 ! flux et gradient des variables de diffusion
  double precision,dimension(3,algo_param%nVarD),intent(out) :: MS_zD,ME_zD, EctS_zD, EctE_zD
 ! contrainte de Piola-Kirchhoff et gradient de u
  double precision,dimension(merge(0,1,algo_param%HPP)*9),intent(out) :: MP_z,MG_z, EctP_z,EctG_z
 ! variables internes
  double precision,dimension(nVI_to_print),intent(out)       :: varInt_z,EctVI_z
  integer                                                    :: indice_VI
  integer                                                    :: p

  ! valeurs par voxel des quantites d'interet, temporaires
  double precision,dimension(6)                                :: sig_tmp
  double precision,dimension(6)                                :: def_tmp
  double precision,dimension(3,algo_param%nVarD)               :: FluxD_tmp,gradQD_tmp
  double precision,dimension(merge(0,1,algo_param%HPP)*9)      :: gradu_tmp, PK1_tmp
  double precision                                             :: det_tmp
  double precision                                             :: vi_tmp

  ! indices relatifs aux voxels de la zone
  integer(kind=8)                         :: k, l, indZone_min, indZone_max
  !-----------------------------------------------------------------------------


  !! Les voxels sont ordonnes par zone dans Mattotp(numM)%pos(k)
  !! L'indice du dernier voxel de la zone numZ est donne par 
  !! MattotP(numM)%zone(numZ,1)

  if (numZ==1) indZone_min = 1
  if (numZ >1) indZone_min = MattotP(numM)%zone(numZ-1,1)+1
  indZone_max = MattotP(numM)%zone(numZ,1)
  nVox_z = dble(indZone_max-indZone_min+1)
  MS_z=0
  ME_z=0
  EctS_z=0
  EctE_z=0
  Sdet_z=0
  MP_z=0
  MG_z=0
  EctP_z=0
  EctG_z=0
  varInt_z=0
  EctVI_z=0
  MS_zD=0
  ME_zD=0
  EctS_zD=0
  EctE_zD=0

  do k=indZone_min,indZone_max
     l=MattotP(numM)%pos(k)
     if(algo_param%Mechanics) then   !------- MECANIQUE 
        sig_tmp=dble(sig(l,:))
        if(algo_param%HPP) then
           def_tmp=dble(def(l,:))
           !! Somme des contraintes et deformations sur la zone
           MS_z=MS_z+sig_tmp
           ME_z=ME_z+def_tmp
           !! Somme sur la zone des carres des contraintes, deformations
           EctS_z=EctS_z+sig_tmp*sig_tmp
           EctE_z=EctE_z+def_tmp*def_tmp
        else
           gradu_tmp=dble(Def(l,:))
           call detF(gradu_tmp,det_tmp)
           !det_tmp=1._mytype ! permet de shunter l'effet de la variation de volume sur 
                              ! moy et ecart-types de la contrainte de Cauchy
           PK1_tmp=dble(PK1(l,:))
           !! On calcule le tenseur de Green-Lagrange en ce point
           def_tmp(1) = gradu_tmp(1) + 0.5_8 * (gradu_tmp(1)**2 &
                + gradu_tmp(7)**2 + gradu_tmp(8)**2)
           def_tmp(2) = gradu_tmp(2) + 0.5_8 * (gradu_tmp(2)**2 &
                + gradu_tmp(4)**2 + gradu_tmp(9)**2)
           def_tmp(3) = gradu_tmp(3) + 0.5_8 * (gradu_tmp(3)**2 &
                + gradu_tmp(5)**2 + gradu_tmp(6)**2)
           def_tmp(4) = 0.5_8 * (gradu_tmp(4) + gradu_tmp(7) + &
                gradu_tmp(1)*gradu_tmp(4)+ gradu_tmp(7)*gradu_tmp(2) &
                + gradu_tmp(8)*gradu_tmp(9))
           def_tmp(5) = 0.5_8 * (gradu_tmp(5) + gradu_tmp(8) + &
                gradu_tmp(1)*gradu_tmp(5) + gradu_tmp(7)*gradu_tmp(6) + &
                gradu_tmp(8)*gradu_tmp(3))
           def_tmp(6) = 0.5_8 * (gradu_tmp(6) + gradu_tmp(9) + &
                gradu_tmp(4)*gradu_tmp(5) + gradu_tmp(2)*gradu_tmp(6) + &
                gradu_tmp(9)*gradu_tmp(3))
           !! Somme sur la zone des contraintes, deformations,...
           !! pour Cauchy, on prend en compre la variation de volume
           Sdet_z = Sdet_z + det_tmp
           MS_z = MS_z + sig_tmp * det_tmp
           ME_z = ME_z + def_tmp
           MP_z = MP_z + PK1_tmp
           MG_z = MG_z + gradu_tmp
           !! Somme sur la zone des carres des contraintes, deformations,...
           EctS_z = EctS_z + sig_tmp*sig_tmp * det_tmp
           EctE_z = EctE_z + def_tmp*def_tmp
           EctP_z = EctP_z + PK1_tmp*PK1_tmp
           EctG_z = EctG_z + gradu_tmp*gradu_tmp
        end if                                                            ! end if HPP/GT
        !! On parcourt les variables internes
        indice_VI = 1
        do p=1,size(extract%varZSTD(i)%val)
           ! si on doit extraire la variable interne
           if(extract%varZSTD(i)%val(p) .AND. nVI_to_print >0) then  
              vi_tmp=dble(MattotP(numM)%VarInt(p,k))
              ! Somme sur la zone
              varInt_z(indice_VI)=varInt_z(indice_VI)+vi_tmp
              !! Somme sur la zone du carre de la variable interne
              EctVI_z(indice_VI)=EctVI_z(indice_VI)+vi_tmp*vi_tmp
              indice_VI = indice_VI + 1
           end if
        end do
     end if                          !-------FIN MECANIQUE
     if (algo_param%Diffusion) then  !-------DIFFUSION
        FluxD_tmp=dble(FluxD(l,:,:))
        gradQD_tmp=dble(gradQD(l,:,:))
        !! Somme sur la zone des contraintes, deformations
        MS_zD=MS_zD+FluxD_tmp
        ME_zD=ME_zD+GradQD_tmp
        !! Somme sur la zone des carres des contraintes, deformations
        EctS_zD=EctS_zD+FluxD_tmp*FluxD_tmp
        EctE_zD=EctE_zD+gradQD_tmp*gradQD_tmp
     end if                          !-------FIN DIFFUSION
  end do


end subroutine calcul_somme_zone



!===============================================================================
!      SUBROUTINE CALCUL_SOMME_COMPOSITE
!-------------------------------------------------------------------------------
!> Subroutine qui calcule les contributions des grandeurs d'interet dans les voxels composites.
!! Prend en entree un numero (global) de materiau et de zone et ressort la somme des grandeurs
!! d'interet sur tous les voxels composites du pinceau composes de ladite zone. Ces valeurs seront 
!! ensuite ajoutees dans les buffers pour l'ecriture des fichiers .std.
!! Fonction appelee dans sortie_std pour chaque zone de chaque materiau.
!!
!!
!! \param[in]        numM            (entier)             Numero (global) du materiau
!! \param[in]        numZ            (entier)             Numero (global) de la zone
!! \param[in]        nVI_to_print    (entier)             Nombre de variables internes a afficher
!!
!! \param[inout]       nVox_z          (double precision)
!! \param[inout]       MS_z,  EctS_z   (double precision)   Moyenne/Ecart-type contrainte de Cauchy
!! \param[inout]       ME_z,  EctE_z   (double precision)   Moyenne/Ecart-type deformation Green-Lagrange
!! \param[inout]       Sdet_z          (double precision)   Somme des determinants de la transformation
!! \param[inout]       MP_z,  EctP_z   (double precision)   Moyenne/Ecart-type Piola-Kirchhoff
!! \param[inout]       MG_z,  EctG_z   (double precision)   Moyenne/Ecart-type Gradient de u
!! \param[inout]       MS_zD, EctS_zD  (double precision)   Moyenne/Ecart-type Flux diffusion
!! \param[inout]       ME_zD, EctE_zD  (double precision)   Moyenne/Ecart-type Gradient diffusion
!! \param[inout]       varInt_z, EctVI_z (double precision) Moyenne/Ecart-type Variables internes
!!
!!
!! REMARQUE :
!!
!! A ce stade, les "moyennes" sont les sommes des valeurs sur les voxels et les "ecarts-types" sont
!! les sommes des carres des valeurs. On ne les transformera en moyenne/ecart-type a proprement
!! parler que lors de l'ecriture dans les fichiers.
!!
!===============================================================================

subroutine calcul_somme_composite(numM, numZ, nVI_to_print, MS_z, ME_z, EctS_z, EctE_z,&
                                                               nVox_z, varInt_z, EctVI_z)

  implicit none

  ! parametres d'entree
  integer,intent(in)                                                  :: numM
  integer(kind=8),intent(in)                                          :: numZ
  integer,intent(in)                                                  :: nVI_to_print

  !moyenne et ecart type des contraintes et deformations par zone, et variables internes
  double precision,dimension(6),intent(inout)                         :: MS_z, EctS_z, ME_z, EctE_z
  double precision,allocatable,dimension(:),intent(inout)             :: varInt_z, EctVI_z
  double precision,intent(inout)                                      :: nVox_z

  ! le nombre de voxels composites contenant numM/numZ 
  integer                                                             :: nVox_comp
  integer                                                             :: voxel ! indice de boucle

  ! numero de materiau / de zone / position des voxels composites contenant numM/numZ 
  integer                                                             :: mat_comp
  integer(kind=8)                                                     :: zone_comp
  integer(kind=8)                                                     :: pos_voxel_comp

  ! nombre de phases du materiau composite
  integer                                                             :: nPhase 
  ! indice de la phase que constitie numM/numZ dans le materiau composite
  integer                                                             :: indice_phase
  ! indice pour trouver les donnees relatives a la phase dans le tableau des coefficients
  ! du voxel composite
  integer                                                             :: indice_coeff_phase

  ! matrices de contrainte et deformation locales (par phase) et globales (sur le voxel composite)
  real(mytype), dimension(6)                               :: Sig_loc, Def_loc, Sig_glob, Def_glob
  ! variables internes
  real(mytype), allocatable, dimension(:)                  :: varInt_tmp
  integer                                                  :: indiceVI, p
  ! fraction volumique de numM/numZ dans le voxel
  double precision                                         :: frac_vol 

  ! pour le modele laminate :
  ! vecteurs normal et tangentiel de la surface
  real(mytype), dimension(3)                                          :: n, t0, t1
  ! matrices de passage du repere global (x,y,z) au repere local (n,t0,t1)
  real(mytype), dimension(6,6)                             :: RotEps,iRotEps,RotSig,iRotSig
  ! Indices d'identification dans le repère de l'interface (modele Multicouches)
  integer, dimension(3)                                    :: IndAP = (/1,4,5/)
  integer, dimension(3)                                    :: IndPl = (/2,3,6/)


  allocate(varInt_tmp(nVI_to_print))
  varInt_tmp = 0

  Sig_loc = 0
  Def_loc = 0

  nVox_comp = Matcomp_sortie(numM)%zone(numZ)%nombre_occurences

  do voxel = 1, nVox_comp
     ! On lit directement le numero de materiau du voxel composite,
     ! ainsi que le numero de phase du materiau qui nous interesse dans le voxel coposite
     mat_comp = Matcomp_sortie(numM)%zone(numZ)%numM_comp(voxel)
     indice_phase = Matcomp_sortie(numM)%zone(numZ)%numPhase_comp(voxel)
     ! numZ_comp contient le numero local de zone du materiau composite en question.
     ! c'est donc l'indice dans le tableau zone(:,2) du vrai numero de zone, ou encore
     ! l'indice dans le tableau pos du voxel en question (car on a ici 1zone = 1voxel)
     zone_comp = Matcomp_sortie(numM)%zone(numZ)%numZ_comp(voxel)
     pos_voxel_comp = MattotP(mat_comp)%pos(zone_comp)
     ! zone_comp est desormais le 'vrai' numero de zone du voxel composite

     nPhase = Mattotp(mat_comp)%nPhase
     ! la fraction volumique de la derniere phase du materiau composite n'apparait pas dans
     ! le tableau des coefficients 
     if (indice_phase == nPhase) then 
        frac_vol = 1. - sum(MattotP(mat_comp)%Coeff(2+nPhase:2*nPhase,zone_comp))
     else
        frac_vol = MattotP(mat_comp)%Coeff(1+nPhase+indice_phase,zone_comp)
     end if
        


     ! cas 'voigt'
     if (Mattotp(mat_comp)%Lawname=='voigt') then

        ! On trouve la bonne phase dans les coefficients du voxel composite
        indice_coeff_phase = int(sum(MattotP(mat_comp)%Varint(1:indice_phase-1,zone_comp)))
        indice_coeff_phase = 1 + indice_coeff_phase + nPhase + 6*(indice_phase-1)

        ! On recupere la contrainte locale
        Sig_loc = MattotP(mat_comp)%VarInt(indice_coeff_phase:indice_coeff_phase+5,zone_comp)

        ! La deformation locale vaut la deformation globale (voigt)
        Def_loc = def(pos_voxel_comp,:)

        ! Variables internes
        indiceVI = 1
        do p = 1, size(extract%varZSTD(numM)%val)
           if (extract%varZSTD(numM)%val(p) .AND. nVI_to_print >0) then
              varInt_tmp(indiceVI) = MattotP(mat_comp)%VarInt(indice_coeff_phase+5+p,zone_comp)
              indiceVI = indiceVI + 1
           end if
        end do
     ! FIN Voigt


     ! cas 'reuss'
     else if (Mattotp(mat_comp)%Lawname=='reuss') then

        ! On trouve la bonne phase dans les coefficients du voxel composite
        indice_coeff_phase = int(sum(MattotP(mat_comp)%Varint(1:indice_phase-1,zone_comp)))
        indice_coeff_phase = 1 + indice_coeff_phase + nPhase + 12*(indice_phase-1)

        ! On recupere la deformation locale
        Def_loc = MattotP(mat_comp)%VarInt(indice_coeff_phase:indice_coeff_phase+5,zone_comp)

        ! La contrainte locale vaut la contrainte globale (reuss)
        Sig_loc = sig(pos_voxel_comp,:)

        ! Variables internes
        indiceVI = 1
        do p = 1, size(extract%varZSTD(numM)%val .AND. nVI_to_print >0 )
           if (extract%varZSTD(numM)%val(p)) then
              varInt_tmp(indiceVI) = MattotP(mat_comp)%VarInt(indice_coeff_phase+11+p,zone_comp)
              indiceVI = indiceVI + 1
           end if
        end do
     ! FIN Reuss


     ! cas 'laminate'
     else if (Mattotp(mat_comp)%Lawname=='laminate') then

        ! vecteurs de la base locale
        n  = MattotP(mat_comp)%Coeff(2*nPhase+1:2*nPhase+3,zone_comp)
        t0 = MattotP(mat_comp)%Coeff(2*nPhase+4:2*nPhase+6,zone_comp)
        t1(2) = n(3)*t0(1)-t0(3)*n(1)
        t1(3) = n(1)*t0(2)-t0(1)*n(2)

        ! On trouve la bonne phase dans les coefficients du voxel composite
        indice_coeff_phase = int(sum(MattotP(mat_comp)%Varint(1:indice_phase-1,zone_comp)))
        indice_coeff_phase = 1 + indice_coeff_phase + nPhase + 9*(indice_phase-1)

        ! On prend la contrainte et la deformation globale (moyenne sur le voxel), on la bascule dans le repere local,
        ! on remplace les composantes planes ou anti-planes par les valeurs locales stockees
        ! dans les variables internes du voxel composite pour avoir les valeurs locales et on 
        ! rebascule dans le repere global.

        ! Donnees globales dans le repere global
        Def_glob = def(pos_voxel_comp,:)
        Sig_glob = sig(pos_voxel_comp,:)
        ! Donnees globales dans le repere local
        call CalcRotMatrices(n,t0,RotEps,IRotEps,RotSig,IRotSig)
        Def_glob = matmul(RotEps,Def_glob)
        Sig_glob = matmul(RotSig,Sig_glob)
        ! recuperation de la deformation et contrainte locales
        Def_loc(IndAP) = MattotP(mat_comp)%VarInt(indice_coeff_phase:indice_coeff_phase+2,zone_comp)
        Def_loc(Indpl) = Def_glob(IndPl)
        Sig_loc(IndPl) = MattotP(mat_comp)%VarInt(indice_coeff_phase+3:indice_coeff_phase+5,zone_comp)
        Sig_loc(IndAP) = Sig_glob(IndAP)
        ! Passage dans le repere global
        Def_loc = matmul(iRotEps,Def_loc)
        Sig_loc = matmul(iRotSig,Sig_loc)

        ! Variables internes
        indiceVI = 1
        do p = 1, size(extract%varZSTD(numM)%val .AND. nVI_to_print >0)
           if (extract%varZSTD(numM)%val(p)) then
             VarInt_tmp(indiceVI) = MattotP(mat_comp)%VarInt(indice_coeff_phase+8+p,zone_comp)
             indiceVI = indiceVI + 1
           end if
        end do
     end if
     ! FIN Laminate

     ! On rajoute la contribution du voxel composite
     nVox_z   = nVox_z   + dble(frac_vol)
     MS_z     = MS_z     + dble(frac_vol*Sig_loc)
     ME_z     = ME_z     + dble(frac_vol*Def_loc)
     EctS_z   = EctS_z   + dble(frac_vol*Sig_loc*Sig_loc)
     EctE_z   = EctE_z   + dble(frac_vol*Def_loc*Def_loc)
     if (nVI_to_print>0) then
       VarInt_z = VarInt_z + dble(frac_vol*VarInt_tmp)
       EctVI_z  = EctVI_z  + dble(frac_vol*VarInt_tmp*VarInt_tmp)
     end if
  end do
  deallocate(varInt_tmp)

end subroutine calcul_somme_composite



!==============================================================================
!       SUBROUTINE REMPLIR_BUFFER
!
! ATTENTION : suite a un bug avec compilateur intel19,remplir_buffer 
!             a ete splitee en deux parties (meca et diffusion)
!
!------------------------------------------------------------------------------
!> Prend les donnees calculees dans calcul_somme_zones et remplis le buffer%z 
!! en mettant a jour les indices qu'il faut, et rajoute les contributions de chaque
!! zone dans buffer%m.
!!
!!
!! \param[inout]           buffer%z, buffer%m, buffer%t
!!
!! \param[inout]           pos_remplissage_buffer_z, pos_remplissage_buffer_m
!!
!! \param[in]              nVox_z               nombre de voxels constituant la zone
!!
!! \param[in]              numM, numZ           numeros de materiau et de zone consideres
!!
!! \param[in]              ind_tps              pas de temps courant
!!
!!
!! \param[in]    prefixe :   M       pour Moyenne
!!                           Ect     pour Ecart-type
!!
!!               radical :   S       pour la contrainte de Cauchy
!!                           E       pour la deformation de Green-Lagrange
!!                           P       pour le tenseur de Piola-Kirshhoff
!!                           G       pour le gradient de u
!!                           VI      pour les variables internes
!!
!!               suffixe :   _z      signifie que c'est une somme sur les voxels de la zone
!!                           _m      signifie que c'est une somme sur les voxels du materiau
!!                           _zD, _mD            si c'est une variable de diffusion
!!
!!
!! Organisation du buffer : On range les quantites par zone les unes a la suite des autres.
!!                          Pour chaque zone l'organisation est la suivante.
!!
!! La première case contient le nombre de voxels composant la zone / le materiau / la cellule.
!! On range ensuite dans l'ordre la mecanique puis la diffusion (si les deux sont presents).
!! Chaque grandeur prend deux fois sa place dans le buffer car on y stoque sa valeur moyenne 
!! et son ecart-type.
!!
!! Si on a de la mecanique :
!! En HPP on a le tenseur des contraintes de Cauchy de taille 6, puis le tenseur
!! des deformations de taille 6.
!! En grandes deformations on a le tenseur de contraintes de Cauchy de taille 6, la deformation
!! de Green-Lagrange de taille 6, le determinant de la transformation (lui ne prend qu'une case,
!! pas d'ecart-type), les contraintes de Piola-Kirchhoff de taille 9 et le gradient de u de taille 9.
!! Pour le buffer de zones on a enfin les variables internes en moyenne et ecart-type les unes
!! apres les autres.
!!
!! Si on a de la diffusion :
!! On traite chaque variable l'une apres l'autre. On a le flux moyen, son ecart-type, le gradient
!! moyen et son ecart-type, qui ont chacun une taille 3.
!------------------------------------------------------------------------------

subroutine remplir_buffer_mecanique(pos_remplissage_buff_z,&
      pos_remplissage_buff_m,pos_remplissage_buff,Sdet_z,MS_z,ME_z,EctS_z,EctE_z,&
      MP_z,MG_z,EctP_z,EctG_z,nVI_to_print,varInt_z,EctVI_z)

  implicit none

  ! indices relatifs au remplissage des buffers (ce sont des curseurs)
  integer,intent(inout)                                     :: pos_remplissage_buff_z
  integer,intent(inout)                                     :: pos_remplissage_buff_m
  integer,intent(inout)                                     :: pos_remplissage_buff


 ! somme des determinants de F
  double precision,intent(in)                               :: Sdet_z
 ! contrainte de Cauchy
  double precision,dimension(6),intent(in)                  :: MS_z, EctS_z
 ! deformation de Green-Lagrange
  double precision,dimension(6), intent(in)                 :: ME_z, EctE_z
 ! contrainte de Piola-Kirchhoff et gradient de u
  double precision,dimension(merge(0,1,algo_param%HPP)*9),intent(in) :: MP_z,MG_z, EctP_z,EctG_z
 ! variables internes
  integer,intent(in)                                         :: nVI_to_print
  double precision,dimension(nVI_to_print),intent(in)        :: varInt_z, EctVI_z

  if(algo_param%Mechanics) then                                    !------ MECANIQUE

    ! remplissage de la moyenne de la contrainte de Cauchy
    buffer%z(pos_remplissage_buff_z:pos_remplissage_buff_z+5)=MS_z
    pos_remplissage_buff_z = pos_remplissage_buff_z +6
    buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+5) = &
    buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+5) + MS_z
    pos_remplissage_buff_m = pos_remplissage_buff_m +6
    buffer%t(pos_remplissage_buff:pos_remplissage_buff+5) = &
    buffer%t(pos_remplissage_buff:pos_remplissage_buff+5) + MS_z
    pos_remplissage_buff = pos_remplissage_buff +6

    ! remplissage de l'ecart-type de la contrainte de Cauchy
    buffer%z(pos_remplissage_buff_z:pos_remplissage_buff_z+5) = EctS_z
    pos_remplissage_buff_z = pos_remplissage_buff_z +6
    buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+5) = &
    buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+5) + EctS_z
    pos_remplissage_buff_m = pos_remplissage_buff_m +6
    buffer%t(pos_remplissage_buff:pos_remplissage_buff+5) = &
    buffer%t(pos_remplissage_buff:pos_remplissage_buff+5) + EctS_z
    pos_remplissage_buff = pos_remplissage_buff +6

    ! remplissage de la moyenne de la deformation de Green-Lagrange
    buffer%z(pos_remplissage_buff_z:pos_remplissage_buff_z+5) = ME_z
    pos_remplissage_buff_z = pos_remplissage_buff_z +6
    buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+5) = &
    buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+5) + ME_z
    pos_remplissage_buff_m = pos_remplissage_buff_m +6
    buffer%t(pos_remplissage_buff:pos_remplissage_buff+5) = &
    buffer%t(pos_remplissage_buff:pos_remplissage_buff+5) + ME_z
    pos_remplissage_buff = pos_remplissage_buff +6

    ! remplissage de l'ecart-type de la deformation de Green-Lagrange
    buffer%z(pos_remplissage_buff_z:pos_remplissage_buff_z+5) = EctE_z
    pos_remplissage_buff_z = pos_remplissage_buff_z +6
    buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+5) = &
    buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+5) + EctE_z
    pos_remplissage_buff_m = pos_remplissage_buff_m +6
    buffer%t(pos_remplissage_buff:pos_remplissage_buff+5) = &
    buffer%t(pos_remplissage_buff:pos_remplissage_buff+5) + EctE_z
    pos_remplissage_buff = pos_remplissage_buff +6

    if(.NOT.algo_param%HPP) then

      ! remplissage du determinant de F
      buffer%z(pos_remplissage_buff_z)=Sdet_z
      pos_remplissage_buff_z = pos_remplissage_buff_z +1
      buffer%m(pos_remplissage_buff_m) = buffer%m(pos_remplissage_buff_m) + Sdet_z
      pos_remplissage_buff_m = pos_remplissage_buff_m +1
      buffer%t(pos_remplissage_buff) = buffer%t(pos_remplissage_buff) + Sdet_z
      pos_remplissage_buff = pos_remplissage_buff +1

      ! remplissage de la moyenne de la contrainte de Piola-Kirchhoff
      buffer%z(pos_remplissage_buff_z:pos_remplissage_buff_z+8) = MP_z
      pos_remplissage_buff_z = pos_remplissage_buff_z +9
      buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+8) = &
      buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+8) + MP_z
      pos_remplissage_buff_m = pos_remplissage_buff_m +9
      buffer%t(pos_remplissage_buff:pos_remplissage_buff+8) = &
      buffer%t(pos_remplissage_buff:pos_remplissage_buff+8) + MP_z
      pos_remplissage_buff = pos_remplissage_buff +9

      ! remplissage de l'ecart-type de la contrainte de Piola-Kirchhoff
      buffer%z(pos_remplissage_buff_z:pos_remplissage_buff_z+8) = EctP_z
      pos_remplissage_buff_z = pos_remplissage_buff_z +9
      buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+8) = &
      buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+8) + EctP_z
      pos_remplissage_buff_m = pos_remplissage_buff_m +9
      buffer%t(pos_remplissage_buff:pos_remplissage_buff+8) = &
      buffer%t(pos_remplissage_buff:pos_remplissage_buff+8) + EctP_z
      pos_remplissage_buff = pos_remplissage_buff +9

      ! remplissage de la moyenne du gradient de u
      buffer%z(pos_remplissage_buff_z:pos_remplissage_buff_z+8) = MG_z
      pos_remplissage_buff_z = pos_remplissage_buff_z +9
      buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+8) = &
      buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+8) + MG_z
      pos_remplissage_buff_m = pos_remplissage_buff_m +9
      buffer%t(pos_remplissage_buff:pos_remplissage_buff+8) = &
      buffer%t(pos_remplissage_buff:pos_remplissage_buff+8) + MG_z
      pos_remplissage_buff = pos_remplissage_buff +9

      ! remplissage de l'ecart-type du gradient de u
      buffer%z(pos_remplissage_buff_z:pos_remplissage_buff_z+8) = EctG_z
      pos_remplissage_buff_z = pos_remplissage_buff_z +9
      buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+8) = &
      buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+8) + EctG_z
      pos_remplissage_buff_m = pos_remplissage_buff_m +9
      buffer%t(pos_remplissage_buff:pos_remplissage_buff+8) = &
      buffer%t(pos_remplissage_buff:pos_remplissage_buff+8) + EctG_z
      pos_remplissage_buff = pos_remplissage_buff +9

    end if

    ! remplissage des variables internes pour le buffer des zones
    if (nVI_to_print>0) then
      buffer%z(pos_remplissage_buff_z:pos_remplissage_buff_z+2*nVI_to_print-1:2)=varInt_z
      buffer%z(pos_remplissage_buff_z+1:pos_remplissage_buff_z+2*nVI_to_print:2)=EctVI_z
      pos_remplissage_buff_z = pos_remplissage_buff_z + 2*nVI_to_print
    end if
  end if                                                        !------ FIN MECANIQUE

end subroutine remplir_buffer_mecanique

subroutine remplir_buffer_diffusion(pos_remplissage_buff_z,&
      pos_remplissage_buff_m,pos_remplissage_buff,MS_zD,ME_zD,EctS_zD,EctE_zD)

  implicit none

  ! indices relatifs au remplissage des buffers (ce sont des curseurs)
  integer,intent(inout)                                     :: pos_remplissage_buff_z
  integer,intent(inout)                                     :: pos_remplissage_buff_m
  integer,intent(inout)                                     :: pos_remplissage_buff
  integer                                                   :: indice_diff


 ! flux et gradient des variables de diffusion
  double precision,dimension(3,algo_param%nVarD),intent(in) :: MS_zD, ME_zD, EctS_zD, EctE_zD
 

  if (algo_param%Diffusion) then                                !------ DIFFUSION

  ! une variable de diffusion apres l'autre
  do indice_diff=1,algo_param%nVarD
    ! la moyenne du flux
    buffer%z(pos_remplissage_buff_z:pos_remplissage_buff_z+2) = MS_zD(:,indice_diff)
    pos_remplissage_buff_z = pos_remplissage_buff_z +3
    buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+2) = &
    buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+2) + MS_zD(:,indice_diff)
    pos_remplissage_buff_m = pos_remplissage_buff_m +3
    buffer%t(pos_remplissage_buff:pos_remplissage_buff+2) = &
    buffer%t(pos_remplissage_buff:pos_remplissage_buff+2) + MS_zD(:,indice_diff)
    pos_remplissage_buff = pos_remplissage_buff +3

    ! l'ecart-type du flux
    buffer%z(pos_remplissage_buff_z:pos_remplissage_buff_z+2) = EctS_zD(:,indice_diff)
    pos_remplissage_buff_z = pos_remplissage_buff_z +3
    buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+2) = &
    buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+2) + EctS_zD(:,indice_diff)
    pos_remplissage_buff_m = pos_remplissage_buff_m +3
    buffer%t(pos_remplissage_buff:pos_remplissage_buff+2) = &
    buffer%t(pos_remplissage_buff:pos_remplissage_buff+2) + EctS_zD(:,indice_diff)
    pos_remplissage_buff = pos_remplissage_buff +3

    ! la moyenne du gradient
    buffer%z(pos_remplissage_buff_z:pos_remplissage_buff_z+2) = ME_zD(:,indice_diff)
    pos_remplissage_buff_z = pos_remplissage_buff_z +3
    buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+2) = &
    buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+2) + ME_zD(:,indice_diff)
    pos_remplissage_buff_m = pos_remplissage_buff_m +3
    buffer%t(pos_remplissage_buff:pos_remplissage_buff+2) = &
    buffer%t(pos_remplissage_buff:pos_remplissage_buff+2) + ME_zD(:,indice_diff)
    pos_remplissage_buff = pos_remplissage_buff +3

    ! l'ecart-type du gradient
    buffer%z(pos_remplissage_buff_z:pos_remplissage_buff_z+2) = EctE_zD(:,indice_diff)
    pos_remplissage_buff_z = pos_remplissage_buff_z +3
    buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+2) = &
    buffer%m(pos_remplissage_buff_m:pos_remplissage_buff_m+2) + EctE_zD(:,indice_diff)
    pos_remplissage_buff_m = pos_remplissage_buff_m +3
    buffer%t(pos_remplissage_buff:pos_remplissage_buff+2) = &
    buffer%t(pos_remplissage_buff:pos_remplissage_buff+2) + EctE_zD(:,indice_diff)
    pos_remplissage_buff = pos_remplissage_buff +3

  end do
  end if                                                        !------ FIN DIFFUSION

end subroutine remplir_buffer_diffusion


!=================================================================================================
!       SUBROUTINE PRINT_MZSTD
!-------------------------------------------------------------------------------------------------
!> Rassemble et traite toutes les donnees presentes dans le buffer, puis les ecrit dans le 
!! fichier .std, .mstd ou .zstd. 
!!
!> Est appelee quand buffer_z est plein, et quand buffer et buffer_m ont fini d'etre remplis.
!!
!> On utilise une seule et meme fonction pour traiter les 3 buffers. On peut identifier 
!! quel buffer est en train d'etre traite via la variable fic_unit qui vaut Fstd, FMstd 
!! ou FZstd. Les trois buffers sont traites de facon quasi similaire.
!! Les seules differences sont :
!!                   - On imprime le nombre d'iteration a la fin de chaque ligne dans le .std
!!                   - Le nom du fichier zstd doit faire apparaitre le numero de materiau 
!!                     correspondant (essai.std / essai.mstd / essai_1.zstd / essai_6.zstd)
!!                   - Dans le fichier zstd sont imprimees les variables internes demandees
!!
!!
!> \param[inout]   buffer_mz     :le buffer qui contient toutes les donnees de zones/materiaux/cellule
!! \param[inout]   buffer_mz_tmp :variable temporaire tampon (valeur in et out indifferente)
!! \param[in]   nBuff_mz      :le nombre de zones/materiaux presents dans le buffer
!!                             vaut 1 pour le buffer macro (une seule cellule dans la cellule !)
!! \param[in]   tailleBuff_mz :la taille que prend une zone/materiau/cellule dans le buffer
!! \param[in]   nbIT          :(entier) le nombre d'iterations
!! \param[in]   t             :(reel) temps du chargement
!! \param[in]   nVI_to_print  :(entier) nbre de variables internes a afficher (pour le zstd, 0 sinon)
!! \param[in]   fic_local     :(chaine de caracteres) nom du fichier dans lequel on ecrit
!! \param[in]   testPrint     :(logical) pour imprimer ou non dna le fichier
!=================================================================================================


subroutine print_mzstd(buffer_mz,buffer_mz_tmp,nBuff_mz,tailleBuff_mz,nbIt,&
                                          t,nVI_to_print,fic_local,testPrint)

  ! le buffer a afficher
  double precision,dimension(:),allocatable,intent(inout)              :: buffer_mz
  double precision,dimension(:),allocatable,intent(inout)              :: buffer_mz_tmp
  ! le nombre de cases qu'occupe la cellule / une zone / un materiau dans le buffer
  integer,intent(in)                                       :: tailleBuff_mz
  ! le nombre de materiaux ou de zones presents dans le buffer
  integer,intent(in)                                       :: nBuff_mz

  integer,intent(in)                                       :: nbIt         ! nombre d'iterations
  real(mytype), intent(in)                                 :: t            ! temps
  integer,intent(in)                                       :: nVI_to_print ! nombre de variables internes
                                                                           ! a afficher (pour le zstd,
                                                                           ! (vaut 0 pour les autres)
  logical, intent(in)                                      :: testPrint

  ! indices servant de curseur lors de l'ecriture
  integer                                                  :: zone, indice
  integer                                                  :: indice_diff, indice_VI

  double precision                                         :: nVox    ! nombre de voxels
  double precision                                         :: Sdet_mz ! determinant de F

  character(len=200),intent(in)                            :: fic_local ! nom du fichier
  character(len=200)                                       :: format    ! format d'ecriture
  integer                                                  :: fic_unit  ! unite du fichier
  integer                                                  :: ierror    ! erreur MPI

  !-----------------------------------------------------------------------------------------------

  ! On somme les valeurs de chaque processeur pour avoir les valeurs globales par zone/materiau/cellule

  call MPI_BARRIER(amitex_mpi_comm,ierror)
  call MPI_Reduce(buffer_mz(1:nBuff_mz*tailleBuff_mz),buffer_mz_tmp(1:nBuff_mz*tailleBuff_mz),&
                  nBuff_mz*tailleBuff_mz,MPI_DOUBLE_PRECISION,MPI_SUM,0,amitex_mpi_comm,ierror)
  call MPI_BARRIER(amitex_mpi_comm,ierror)
  buffer_mz = buffer_mz_tmp


  !-----------------------------------------------------------------------------------------------
  ! On divise par le nombre de voxels pour transformer les sommes en moyennes et on transforme
  ! les sommes de carres en ecarts-types, en verifiant qu'elles soient bien positives
  if(nrank==0) then
    do zone=1,nBuff_mz
      indice = 2 + (zone-1)*tailleBuff_mz
      nVox = buffer_mz(indice-1)
      if (algo_param%Mechanics) then                                 !------- MECANIQUE
        if (algo_param%HPP) then                                     !------- HPP
          ! on divise par le nombre de voxels pour avoir les moyennes de Cauchy et Green-Lagrange
          buffer_mz(indice:indice+23) = &
          buffer_mz(indice:indice+23)/nVox
          ! variance Cauchy
          buffer_mz(indice+6:indice+11) = buffer_mz(indice+6:indice+11) - &
          buffer_mz(indice:indice+5)*buffer_mz(indice:indice+5)
          indice = indice + 12
          ! variance Green-Lagrange
          buffer_mz(indice+6:indice+11) = buffer_mz(indice+6:indice+11) - &
          buffer_mz(indice:indice+5)*buffer_mz(indice:indice+5)
          indice = indice - 12
          ! calcul ces ecarts-types

          call var_to_std(buffer_mz(indice:indice+5),buffer_mz(indice+6:indice+11),6,eps)
          indice = indice + 12
          call var_to_std(buffer_mz(indice:indice+5),buffer_mz(indice+6:indice+11),6,eps)
          indice = indice + 12

        else                                                          !------- GRANDES DEF
          Sdet_mz = buffer_mz(indice+24)
          buffer_mz(indice:indice+11)=buffer_mz(indice:indice+11)/Sdet_mz     !moyenne Cauchy
          buffer_mz(indice+12:indice+23)=buffer_mz(indice+12:indice+23)/nVox !moyenne Green-Lagrange
          buffer_mz(indice+25:indice+60)=buffer_mz(indice+25:indice+60)/nVox !moyenne Piola-Kirchhoff
                                                                                      ! et grad u
          buffer_mz(indice+6:indice+11) = buffer_mz(indice+6:indice+11) - &
          buffer_mz(indice:indice+5)*buffer_mz(indice:indice+5)           ! variance Cauchy

          buffer_mz(indice+18:indice+23) = buffer_mz(indice+18:indice+23) - &
          buffer_mz(indice+12:indice+17)*buffer_mz(indice+12:indice+17)   ! variance Green-Lagrange

          buffer_mz(indice+34:indice+42) = buffer_mz(indice+34:indice+42) - &
          buffer_mz(indice+25:indice+33)*buffer_mz(indice+25:indice+33)   ! variance Piola-Kirchhoff

          buffer_mz(indice+52:indice+60) = buffer_mz(indice+52:indice+60) - &
          buffer_mz(indice+43:indice+51)*buffer_mz(indice+43:indice+51)   ! variance grad u

          ! calcul ces ecarts-types
          call var_to_std(buffer_mz(indice:indice+5),buffer_mz(indice+6:indice+11),6,eps)
          call var_to_std(buffer_mz(indice+12:indice+17),buffer_mz(indice+18:indice+23),6,eps)
          call var_to_std(buffer_mz(indice+25:indice+33),buffer_mz(indice+34:indice+42),9,eps)
          call var_to_std(buffer_mz(indice+43:indice+51),buffer_mz(indice+52:indice+60),9,eps)
          indice = indice + 61
        end if

        ! les variables internes
        buffer_mz(indice:indice+2*nVI_to_print-1)=buffer_mz(indice:indice+2*nVI_to_print-1)/nVox !moyenne
        do indice_VI = 1 , nVI_to_print
          buffer_mz(indice+1) = buffer_mz(indice+1) - buffer_mz(indice) * buffer_mz(indice) !variance
          call var_to_std(buffer_mz(indice),buffer_mz(indice+1),1,eps) !ecart-type
          indice = indice + 2
        end do
      end if                                                  !------- FIN MECANIQUE


      if (algo_param%Diffusion) then                          !------- DIFFUSION
          buffer_mz(indice:indice+12*algo_param%nVarD-1) = &
          buffer_mz(indice:indice+12*algo_param%nVarD-1) / nVox
          do indice_diff = 1,2*algo_param%nVarD
            buffer_mz(indice+3:indice+5) = buffer_mz(indice+3:indice+5) - &
            buffer_mz(indice:indice+2)*buffer_mz(indice:indice+2)
            call var_to_std(buffer_mz(indice:indice+2),buffer_mz(indice+3:indice+5),3,eps)
            indice = indice + 6
          end do
      end if                                                  !------- FIN DIFFUSION
    end do

    ! fin du traitement des donnees
    !-----------------------------------------------------------------------------------
    ! debut du print dans le fichier .(m/z)std
    if (testPrint) then
    open(newunit=fic_unit, file=trim(fic_local),form="formatted", status="old",&
                                            position="append", action="write")

    do zone=1,nBuff_mz
      write(fic_unit, fmt="(E15.8,1X)", advance="no") t
      indice = 2 + (zone-1)*tailleBuff_mz
      if (algo_param%Mechanics) then                                 !------- MECANIQUE
        if (algo_param%HPP) then
          write(fic_unit, fmt="(12(E15.8,1X),12E15.8)", advance="no")&  !astuce : ecart-type positifs, on ne rajoute pas de 1X (interet?)
              buffer_mz(indice:indice+5),&     ! ecriture moyenne Cauchy
              buffer_mz(indice+12:indice+17),& ! ecriture moyenne Green-Lagrange
              buffer_mz(indice+6:indice+11),&  ! ecriture ecart-type Cauchy
              buffer_mz(indice+18:indice+23)   ! ecriture ecart-type Green-Lagrange
          indice = indice + 24
        else
          ! ecriture des moyennes de piola-K et gauss, puis de leurs ecarts-types
          write(fic_unit, fmt="(30(E15.8,1X),30E15.8)", advance="no")&  !astuce : ecart-type positifs, on ne rajoute pas de 1X (interet?)
              buffer_mz(indice   :indice+ 5),&   ! ecriture moyenne Cauchy
              buffer_mz(indice+25:indice+33),&   ! ecriture moyenne Piola-Kirchhoff
              buffer_mz(indice+12:indice+17),&   ! ecriture moyenne Green-Lagrange
              buffer_mz(indice+43:indice+51),&   ! ecriture moyenne grad u
              buffer_mz(indice+ 6:indice+11),&   ! ecriture ecart-type Cauchy
              buffer_mz(indice+34:indice+42),&   ! ecriture ecart-type Piola-Kirchhoff
              buffer_mz(indice+18:indice+23),&   ! ecriture ecart-type Green-Lagrange
              buffer_mz(indice+52:indice+60)     ! ecriture ecart-type grad u
          indice = indice + 61
        end if
        ! ecriture des variables internes pour le buffer zones
        if (nVI_to_print>0) then
          format = ""
          write(format,"(I4)") nVI_to_print
          format = "("//trim(format)//"(1X,E15.8,E15.8))"  !astuce : ecart-type positifs, on ne rajoute pas de 1X (interet?)
          write(fic_unit, fmt=trim(format), advance="no") buffer_mz(indice:indice-1+2*nVI_to_print)
          indice = indice + 2*nVI_to_print
        end if
      end if                                                      !------- FIN MECANIQUE
      if (algo_param%Diffusion) then                                 !------- DIFFUSION
        format = ""
        write(format,"(2(A,I4),A)") "(",6*algo_param%nVarD,"(E15.8,1X),",6*algo_param%nVarD,"E15.8)"
        write(fic_unit, fmt=trim(format), advance="no") buffer_mz(indice:indice+3*algo_param%nVarD-1), &
                                     buffer_mz(indice+6*algo_param%nVarD:indice+9*algo_param%nVarD-1), &
                                     buffer_mz(indice+3*algo_param%nVarD:indice+6*algo_param%nVarD-1), &
                                     buffer_mz(indice+9*algo_param%nVarD:indice+12*algo_param%nVarD-1)
        indice = indice + 12*algo_param%nVarD
      end if                                                      !------- FIN DIFFUSION
      !if (fic_unit==Fstd) write(fic_unit, fmt = "(I4)", advance = "no") nbIt ! le nombre d'iterations 
      !                                                                       ! a afficher dans le .std
      if (fic_local(len_trim(fic_local)-3:len_trim(fic_local)) == ".std") then
          write(fic_unit, fmt = "(I6)", advance = "no") nbIt
      end if 
      write(fic_unit, fmt="()", advance="yes") !retour a la ligne
    end do
    flush(fic_unit)
    close(fic_unit)
    end if !fin testPrint
    ! fin du print .(m/z)std
    !-----------------------------------------------------------------------------------

  end if !nrank == 0

end subroutine print_mzstd 

       
       
!==============================================================================
!       SUBROUTINE SORTIE_STD_MACRO
!------------------------------------------------------------------------------
!> Ecriture des donnees sur le fichier de sortie standard (sit testPrint)
!! Fonction  appelee a chaque fin de pas de calcul. \n
!! Permet de calculer et d'ecrire des grandeurs d'interet dans le fichier de sortie standard
!! Uniquement pour la cellule complete :
!!               -> moyenne et écart type de la déformation sur toute la cellule
!!               -> moyenne et écart type de la contrainte sur toute la cellule
!!
!!
!! \param[in]      t: (reel) temps du chargement
!! \param[in]      nb_it: (entier) nombre d'iterations
!! \param[out]     MoyStress, MoyDef: tableaux de reels (dimensions 6 ou 9), contraintes 
!!                       et deformations moyennes utilisees comme valeurs initiales
!!                       pour le pas de chargement suivant. En HPP, cela correspond a la
!!                       contrainte de Cauchy et a la deformation. 
!!                       En grandes transformations, cela correspond a la contrainte
!!                       de Piola-Kirchhoff et au gradient du deplacement
!!                 FluxDMoy, gradQDMoy : tableaux reels (dimensions (3,nVar))
!! \param[in]      testPrint :(logical) pour afficher la sortie dans le .std
!!
!!
!!
!!REMARQUE \n
!!       Le rang du processus est defini par la variable 'nrank'
!!
!!                                                                     IMPORTANT \n
!! Lors de l'appel de "write" il faut utiliser l'unite logique "Fstd". \n
!! DANS LE CAS CONTRAIRE: on ne peut pas garantir que les donnees soient ecrites
!!                       dans le bon fichier.\n
!!
!!!
!! Les differents tableaux (pinceaux en X) sont definies de "amx_start(i)" a "amx_end(i)" 
!! ou i correspond a l'axe: 
!! - 'x' si i=1
!! - 'y' si i=2
!! - 'z' si i=3
!! Les valeurs de amx_start et amx_end varient pour chaque processus.
!!Donc pour atteindre une position specifique (par exemple Sig(2,5,1,:) )
!! on prendra soin de verifier que ces coordonnees existent pour le processus
!! (un processus est defini par son rang).\n
!!
!!------------------------------------------------------------------------------
!!                                                                     IMPORTANT\n
!!Les contraintes et deformations moyennes doivent etres calculees puisqu'elles
!!sont utilisees comme valeurs initiales certains chargement
!!(le premier chargement d'un ensemble defini dans le fichier xml de chargement).
!!
!!
!!------------------------------------------------------------------------------
!!                                                                     IMPORTANT\n
!!Pour des raisons de precision sur les moyennes et les ecarts types, on utilise
!!des variables temporaires de type double precision.
!!
!! Si tel n'est pas le cas, en simple precision, on obtient une valeur de
!! E[X^2]-E[X]^2 negative et donc une valeur d'ecart type de "NaN"
!!
!==============================================================================
subroutine sortie_std_macro(t, nb_it, MoyStress, MoyDef, FluxDMoy, gradQDMoy,testPrint)

  !
  !On considere ici que les phases sont definies par les numeros de materiau
  !


  implicit none
  !------------------------------------------------------------------------------
  !variables non modifiables  
  real(mytype), intent(in)   :: t
  real(mytype),dimension(algo_param%nTensDef),intent(out)  :: MoyDef,MoyStress
  real(mytype),dimension(3,algo_param%nVarD),intent(out)   :: FluxDMoy,gradQDMoy

  integer, intent(in)        :: nb_it
  logical, intent(in)        :: testPrint

  !format d'affichage des reels
  character(len=5),parameter :: FMT_real="E15.8"

  character(len=200)         :: fic_std
  integer                    :: Fstd

  !------------------------------------------------------------------------------
  !variables modifiables
  ! moyennes, moyennes par phase et ecart types

  !moyenne et ecart type des contraintes et deformation sur la cellule
  double precision,dimension(6)                           ::  MS, ME,EctS, EctE
  double precision,dimension(6)                           ::  MS_tmp, EctS_tmp, ME_tmp, EctE_tmp
  double precision,dimension(3,algo_param%nVarD)          ::  MSD, MED,EctSD, EctED
  double precision,dimension(merge(0,1,algo_param%HPP)*9) ::  MP, EctP,MG,EctG
  double precision                                        ::  somme_det, somme_det_tmp
  double precision,dimension(6)                           ::  def_tmp
  double precision                                        ::  det
  double precision, dimension(9)                          ::  gradu_tmp

  ! variables pour la definition des formats
  integer(kind=8)               :: p,p1,p2,p3,p4
  character(len=8)              :: p1_char,p2_char,p3_char,p4_char 
					! nombre de colonnes dans les fichiers => utile pour la construction des "format" d'ecriture
					! mecanique : p1 et p2 (valeurs moyennes et ecarts types)
					! diffusion : p3 et p4 (valeurs moyennes et ecarts types)   
                                        ! astuce : on ne rajoute pas d'espace (1X) pour l'ecriture des ecarts-types (positifs)  
  character(len=200)            :: format
  integer                       :: ierror, i


  !! identification des nombres de colonnes a ecrire
     if (algo_param%HPP) then
       p1 = 13
       p2 = 12
     else
       p1 = 31
       p2 = 30
     end if
     p = 6*algo_param%nVarD
     p3 = p+1
     p4 = p
     write(p1_char,"(I4)") p1 
     write(p2_char,"(I4)") p2 
     write(p3_char,"(I4)") p3 
     write(p4_char,"(I4)") p4 


  !! Ouverture de fichier std
  !fic_std = trim(fic_vtk)//".std"
  !if(nrank==0)then
  !   !fichier .std
  !   open(newunit=Fstd, file=fic_std,form="formatted", status="old",position="append", action="write")
  !end if

  !! INITIALISATION a 0 
  !! Moyennes et ecart-types des contraintes (S) et des deformations (E) 
  !! dans tout le domaine

  if (algo_param%Mechanics) then !-------MECANIQUE
    MS=0       !Moyennes
    ME=0
    EctS=0     !Ecart-types
    EctE=0

    if(.NOT.algo_param%HPP) then !contraintes PK (P) et gradient de U (G) en Grandes Transfo.
       MP=0
       MG=0
       EctP=0
       EctG=0
    end if
  end if
  if (algo_param%Diffusion) then !-------DIFFUSION
    MSD=0
    MED=0
    EctSD=0
    EctED=0
  end if

  !! CALCUL DES MOYENNES
  
  if (algo_param%Mechanics) then !-------MECANIQUE

     if(algo_param%HPP) then

        call field_Mean(Sig,grid%ntot,6,MS,EctS) ! Moyenne et Variance  contrainte de Cauchy
        call field_Mean(Def,grid%ntot,6,ME,EctE) ! Moyenne déformation 

     else !contraintes PK (P) et gradient de U (G) en Grandes Transfo.

        somme_det = 0.
        do i = 1, size(sig(:,1))
        
           gradu_tmp = Def(i,:)
           call detF(gradu_tmp,det)
           somme_det = somme_det + det
           
             !! On calcule le tenseur de Green-Lagrange
           def_tmp(1) = dble(Def(i,1) + 0.5_8 * (Def(i,1)**2 &
                       + Def(i,7)**2 + Def(i,8)**2))
           def_tmp(2) = dble(Def(i,2) + 0.5_8 * (Def(i,2)**2 &
                       + Def(i,4)**2 + Def(i,9)**2))
           def_tmp(3) = dble(Def(i,3) + 0.5_8 * (Def(i,3)**2 &
                       + Def(i,5)**2 + Def(i,6)**2))
           def_tmp(4) = dble(0.5_8 * (Def(i,4) + Def(i,7) + &
                       Def(i,1)*Def(i,4)+ Def(i,7)*Def(i,2) &
                       + Def(i,8)*Def(i,9)))
           def_tmp(5) = dble(0.5_8 * (Def(i,5) + Def(i,8) + &
                       Def(i,1)*Def(i,5) + Def(i,7)*Def(i,6) + &
                       Def(i,8)*Def(i,3)))
           def_tmp(6) = dble(0.5_8 * (Def(i,6) + Def(i,9) + &
                       Def(i,4)*Def(i,5) + Def(i,2)*Def(i,6) + &
                       Def(i,9)*Def(i,3)))
                       
           MS = MS + sig(i,:)*det
           ME = ME + def_tmp
           EctS = EctS + sig(i,:)*sig(i,:)*det
           EctE = EctE + def_tmp*def_tmp
           
        end do
     
        call MPI_Allreduce(somme_det,somme_det_tmp,1,MPI_DOUBLE_PRECISION,&
                                                             MPI_SUM,amitex_mpi_comm,ierror)
        call MPI_Allreduce(EctS,EctS_tmp,6,MPI_DOUBLE_PRECISION,MPI_SUM,amitex_mpi_comm,ierror)
        call MPI_Allreduce(MS,MS_tmp,6,MPI_DOUBLE_PRECISION,MPI_SUM,amitex_mpi_comm,ierror)
        MS = MS_tmp/somme_det_tmp
        EctS = EctS_tmp/somme_det_tmp - MS*MS
        call var_to_std(MS,EctS,6,eps)
        
        call MPI_Allreduce(EctE,EctE_tmp,6,MPI_DOUBLE_PRECISION,MPI_SUM,amitex_mpi_comm,ierror)
        call MPI_Allreduce(ME,ME_tmp,6,MPI_DOUBLE_PRECISION,MPI_SUM,amitex_mpi_comm,ierror)
        ME = ME_tmp / grid%ntot
        EctE = EctE_tmp / grid%ntot - ME*ME
        call var_to_std(ME,EctE,6,eps)

        call field_Mean(Def,grid%ntot,9,MG,EctG)      ! Moyenne et Variance du gradient de U
        call field_Mean(PK1,grid%ntot,9,MP,EctP)      ! Moyenne et Variance contrainte Piola Kirchoff 1

     end if
  end if                  !---------FIN MECANIQUE
  
  if (algo_param%Diffusion) then !---------DIFFUSION

     call field_Mean(FluxD,grid%ntot,3*algo_param%nVarD,MSD,EctSD) ! Moyenne et Variance  contrainte de Cauchy
     call field_Mean(GradQD,grid%ntot,3*algo_param%nVarD,MED,EctED) ! Moyenne déformation 

  end if                  !---------FIN DIFFUSION 

  !! ECRITURE Fichier .std
  if (nrank==0) then
  if(testPrint) then
     fic_std = trim(fic_vtk)//".std"
     open(newunit=Fstd, file=fic_std,form="formatted", status="old",position="append", action="write")
     if (algo_param%Mechanics .and. .not.algo_param%Diffusion) then !-Mecanique pure
        format="("//trim(p1_char)//"(E15.8,1X),"//trim(p2_char)//"E15.8,I4)"
        if(algo_param%HPP) then
           write(Fstd,fmt=trim(format)) t,MS,ME,EctS,EctE,nb_it
        else 
           write(Fstd,fmt=trim(format)) t,MS,MP,ME,MG,EctS,EctP,EctE,EctG,nb_it
        end if
     end if
     if (algo_param%Diffusion .and. .not.algo_param%Mechanics) then  !-Diffusion pure
        format="("//trim(p3_char)//"(E15.8,1X),"//trim(p4_char)//"E15.8,I4)"
        write(Fstd,fmt=trim(format)) t,MSD(:,1),MED(:,1),EctSD(:,1),EctED(:,1),nb_it
     end if
     flush(Fstd)
     close(Fstd)
  end if
  end if 
  

  !! on affecte les grandeurs moyennes de sortie de la procédure
  if(algo_param%Mechanics) then
       if(algo_param%HPP) then
          MoyDef = real(ME,mytype)
          MoyStress = real(MS,mytype)
       else
          MoyDef = real(MG,mytype)
          MoyStress = real(MP,mytype)
       end if
  end if
  if(algo_param%Diffusion) then
          GradQDMoy = real(MED,mytype)
          FluxDMoy = real(MSD,mytype)
  end if 

end subroutine sortie_std_macro


!==============================================================================
!       SUBROUTINE DETF
!------------------------------------------------------------------------------
!> Calcul du determinant de (1 + G) 
!!      ou G est une matrice 3x3 donnee sous forme de vecteur
!! 
!! \param[in]     G sous forme de vecteur selon l'ordre 
!!                11 22 33 12 13 23 21 31 32
!!
!! \param[out]    res, determinant de (1 + G) 
!!
!==============================================================================
subroutine detF(G,res)

  implicit none
  real(mytype),dimension(9), intent(in)  :: G
  real(mytype), intent(out)              :: res  


! copie de det.eso de cast3m
!        DET= A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
!           +A(2,1)*(A(3,2)*A(1,3)-A(1,2)*A(3,3)) &
!           +A(3,1)*(A(1,2)*A(2,3)-A(2,2)*A(1,3)) &

  !calcul
  res = (G(1)+1._mytype)*((G(2)+1._mytype)*(G(3)+1._mytype)-G(9)*G(6)) &
      + G(7)*(G(9)*G(5)-G(4)*(G(3)+1._mytype)) &
      + G(8)*(G(4)*G(6)-(G(2)+1._mytype)*G(5))



end subroutine detF



!===============================================================================
!    SUBROUTINE DESALLOCATION_SORTIE_STD
!-------------------------------------------------------------------------------
!> Desaloue les buffers et la structure Matcomp_sortie
!!
!! Interet : Garder les variables 'Matcomp_sortie', 'buffer', 'buffer_z'...
!!           internes au module sortie_std tout en les desallouant a la fin du
!!           programme principal amitex_fftp.
!!
!-------------------------------------------------------------------------------
subroutine desallocation_sortie_std

  implicit none

integer    :: i, j

  ! Desallocation des buffers
  if (allocated(buffer0%t)) deallocate(buffer0%t)
  if (allocated(buffer0%t_tmp)) deallocate(buffer0%t_tmp)
  if (allocated(buffer0%m)) deallocate(buffer0%m)
  if (allocated(buffer0%m_tmp)) deallocate(buffer0%m_tmp)
  if (allocated(buffer0%z)) deallocate(buffer0%z)
  if (allocated(buffer0%z_tmp)) deallocate(buffer0%z_tmp)

  ! Desallocation de Matcomp_sortie
  if(allocated(Matcomp_sortie0)) then
  do i = 1, size(Matcomp_sortie0)
     do j = 1, size(Matcomp_sortie0(i)%zone)
        if (allocated(Matcomp_sortie0(i)%zone(j)%numM_comp)) then
            deallocate(Matcomp_sortie0(i)%zone(j)%numM_comp)
        end if
        if (allocated(Matcomp_sortie0(i)%zone(j)%numZ_comp)) then
            deallocate(Matcomp_sortie0(i)%zone(j)%numZ_comp)
        end if
        if (allocated(Matcomp_sortie0(i)%zone(j)%numPhase_comp)) then
            deallocate(Matcomp_sortie0(i)%zone(j)%numPhase_comp)
        end if
     end do
     if(allocated(Matcomp_sortie0(i)%zone))deallocate(Matcomp_sortie0(i)%zone)
  end do
  deallocate(Matcomp_sortie0)
  end if


end subroutine desallocation_sortie_std



end module sortie_std_mod
