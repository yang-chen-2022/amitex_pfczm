
!====================================================================
!
!       MODULE NON_LOCAL_MOD :
!
!> "White" module dedicated to non-local modelling (or more generally to coupling) 
!! In practice the non-local modelling is developped in a user-defined module : 
!!                         non_local_user_mod.f90    
!! 
!!
!!
!!
!====================================================================

module non_local_mod
  
  use ISO_FORTRAN_ENV

  use MPI
#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi,       only : nrank
#else
  use decomp_2d , only : mytype, nrank
#endif  

  use error_mod

  ! Champ de variables non locales defini dans amitex_mod
  use amitex_mod, only : Nloc, GNloc,& ! pointers
                         NLOC_FIELD,GNLOC_FIELD !Types
  
  use param_algo_mod , only : algo_param

  use material_mod, only : Assign_VarInt_from_Field,Assemble_field_from_VarInt,& ! functions
                           mattotP, Nloc_models !pointers
  use field_mod
  use green_mod

  use non_local_user_mod

  implicit none

  private

  !> variables publiques
  !public :: 

  !> routines publiques
  public :: init_Nloc_variables, Nloc_call, deallocate_Nloc_Variables,init_user_Nloc_variables_call

contains
  
  !=====================================================================================
  !                             SUBROUTINE INIT_NLOC_VARIABLES
  !
  !>  Allocation et initialisation a 0 des champs de variables non locales stockees dans 
  !>                  les variables GNloc(i) et Nloc(i) (i indice du modele non-local)
  !!
  !!
  !=====================================================================================
  subroutine init_Nloc_Variables(Nloc,GNloc)

    implicit none
    
    type(NLOC_FIELD), allocatable, dimension(:), intent(inout)  :: Nloc
    type(GNLOC_FIELD), allocatable, dimension(:), intent(inout) :: GNloc
    
    integer                                :: i,alloc_stat      
       

    !> Allocation du tableau Nloc
    allocate(Nloc(size(Nloc_models)),stat=alloc_stat)
    if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (initNloc)",2) 

    !> Allocation du tableaux GNloc
    allocate(GNloc(size(Nloc_models)),stat=alloc_stat)
    if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (initNloc)",2) 

    do i=1,size(Nloc_models)
       !> Allocation du tableau Nloc%Var pour chaque modele (espace reel)
       allocate(Nloc(i)%Var(amx_size(1)*amx_size(2)*amx_size(3),1:Nloc_models(i)%Nnloc),stat=alloc_stat)
       if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (initNloc)",2)  

       !> Allocation du tableau Nloc%VarF pour chaque modele (espace spectral)
       allocate(Nloc(i)%VarF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),&
            1:Nloc_models(i)%Nnloc),stat=alloc_stat)
       if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (initNloc)",2)


       !> Initialisation des tableaux à 0
       Nloc(i)%Var = 0.
       Nloc(i)%VarF = 0.

       !> Allocation du tableau GNloc%Var pour chaque modele (espace reel)
       allocate(GNloc(i)%Var(amx_size(1)*amx_size(2)*amx_size(3),1:Nloc_models(i)%NGnloc),stat=alloc_stat)
       if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (initNloc)",2)  

       !> Allocation du tableau GNloc%VarF pour chaque modele (espace spectral)
       allocate(GNloc(i)%VarF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),&
            1:Nloc_models(i)%NGnloc),stat=alloc_stat)
       if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (initNloc)",2)


       !> Initialisation des tableaux à 0
       GNloc(i)%Var = 0.
       GNloc(i)%VarF = 0.

    end do

  end subroutine init_Nloc_Variables

!==================================================================================
!                             SUBROUTINE DEALLOCATE_NLOC
!
!> Desalloue les champs de variables non locales stockees dans 
!>                les variables GNloc(i) et Nloc(i) (i indice du modele non-local)
!!
!==================================================================================
subroutine deallocate_Nloc_variables(Nloc,GNloc)

  implicit none
  
  type(NLOC_FIELD), allocatable, dimension(:), intent(inout)  :: Nloc
  type(GNLOC_FIELD), allocatable, dimension(:), intent(inout) :: GNloc

  integer  :: i   

  if (allocated(Nloc)) then
     do i =1,size(Nloc)
       if (allocated(Nloc(i)%Var)) deallocate(Nloc(i)%Var)
       if (allocated(Nloc(i)%Var)) deallocate(Nloc(i)%VarF)
     end do
     deallocate(Nloc)
  end if

  if (allocated(GNloc)) then
     do i =1,size(GNloc)
       if (allocated(GNloc(i)%Var)) deallocate(GNloc(i)%Var)
       if (allocated(GNloc(i)%Var)) deallocate(GNloc(i)%VarF)
     end do
     deallocate(GNloc)
  end if 
end subroutine deallocate_Nloc_variables


!==================================================================================
!            SUBROUTINE INIT_USER_NLOC_VARIABLES_CALL
!
!> Appel des modeles non locaux 'utilisateur' (user_nloc1,user_nloc2,user_nloc3)
!!
!!
!!  ATTENTION : La lecture des .xml par Fox ne prend pas en compte la casse
!!              ==> ECRIRE TOUS LES TESTS EN MINUSCULES !!!!
!!
!==================================================================================
subroutine init_user_Nloc_variables_call()

    implicit none
    integer                                :: i
    integer                                :: alloc_stat
    integer                                :: Nnloc, NGnloc
    integer,allocatable,dimension(:)       :: list_mat          
    integer,allocatable,dimension(:,:)     :: list_varInt_Nloc  
    integer,allocatable,dimension(:,:)     :: list_varInt_GNloc   

    do i=1,size(Nloc_models)

       ! On recupere les modeles 'locaux' impliques dans le modele non-local
       ! et les listes de variables internes associes (une liste par modele 'local')
       !!----------------------------------------------------------------------------
       Nnloc  = Nloc_models(i)%Nnloc
       NGnloc = Nloc_models(i)%NGnloc

       allocate(list_mat(size(Nloc_models(i)%numM)),stat=alloc_stat) 
       if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant 1 (Nloc_call)",2)
       list_mat = 0

       allocate(list_varInt_Nloc(size(Nloc_models(i)%numM),Nnloc),stat=alloc_stat)
       if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant 2 (Nloc_call)",2)
       list_varInt_Nloc = 0   

       allocate(list_varInt_GNloc(size(Nloc_models(i)%numM),NGnloc),stat=alloc_stat)
       if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant 3 (Nloc_call)",2)
       list_varInt_GNloc = 0   

       list_mat = Nloc_models(i)%numM
       list_varInt_Nloc  = Nloc_models(i)%Ind_VarNloc
       list_varInt_GNloc = Nloc_models(i)%Ind_VarGNloc

       ! Appel du modele correspondant  : CHAINES TESTS EN MINISCULE (voir plus haut)
       !!---------------------------------

        select case(trim(Nloc_models(i)%Modelname))
        ! pas d'initialisation specifique pour ces trois modeles
        case("test_nloc")

        ! possibilite de faire co-exister 3 modeles "non-locaux" (en general on en a un seul) 
        case("user_nloc1")
           call init_user_nloc_variables1(Nnloc,NGnloc,list_mat,list_VarInt_Nloc,list_VarInt_GNloc,i)
        case("user_nloc2")
           call init_user_nloc_variables2(Nnloc,NGnloc,list_mat,list_VarInt_Nloc,list_VarInt_GNloc,i)
        case("user_nloc3")
           call init_user_nloc_variables3(Nnloc,NGnloc,list_mat,list_VarInt_Nloc,list_VarInt_GNloc,i)
        case default
           call amitex_abort("Modele non local renseigne non reconnu (Nloc_call)",2,0)
        end select

        deallocate(list_varInt_Nloc)
        deallocate(list_varInt_GNloc)
        deallocate(list_mat)
    end do

end subroutine init_user_Nloc_variables_call


!==================================================================================
!                             SUBROUTINE NLOC_CALL
!
!> Appel des modeles non locaux demandes en entree par l'utilisateur
!!
!!  \param[in] before_after:   character(5) "befor" 
!!                                       or "after"
!!
!!
!!  Remarque : lawname  (chaine caracteres) nom du modele non local a appeler
!!
!!  ATTENTION : La lecture des .xml par Fox ne prend pas en compte la casse
!!              ==> ECRIRE TOUS LES TESTS EN MINUSCULES !!!!
!!
!==================================================================================
subroutine Nloc_call(load_n,load_incr,ind_tps,before_after,bool_tab)

    implicit none
    character(len=5),intent(in)                     :: before_after
    integer,intent(in)                              :: load_n,load_incr, ind_tps
    logical,dimension(size(Nloc_models)),intent(in) :: bool_tab

    integer                                :: i
    integer                                :: alloc_stat
    integer                                :: Nnloc, NGnloc
    integer,allocatable,dimension(:)       :: list_mat  
    integer,allocatable,dimension(:,:)     :: list_varInt_Nloc  
    integer,allocatable,dimension(:,:)     :: list_varInt_GNloc   

    do i=1,size(Nloc_models)
    if (bool_tab(i)) then

       ! On recupere les modeles 'locaux' impliques dans le modele non-local
       ! et les listes de variables internes associes (une liste par modele 'local')
       !!----------------------------------------------------------------------------
       Nnloc  = Nloc_models(i)%Nnloc
       NGnloc = Nloc_models(i)%NGnloc

       allocate(list_mat(size(Nloc_models(i)%numM)),stat=alloc_stat) 
       if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant 1 (Nloc_call)",2)
       list_mat = 0

       allocate(list_varInt_Nloc(size(Nloc_models(i)%numM),Nnloc),stat=alloc_stat)
       if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant 2 (Nloc_call)",2)
       list_varInt_Nloc = 0   

       allocate(list_varInt_GNloc(size(Nloc_models(i)%numM),NGnloc),stat=alloc_stat)
       if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant 3 (Nloc_call)",2)
       list_varInt_GNloc = 0   

       list_mat = Nloc_models(i)%numM
       list_varInt_Nloc  = Nloc_models(i)%Ind_VarNloc
       list_varInt_GNloc = Nloc_models(i)%Ind_VarGNloc


       ! Construction du champ de variables locales a deriver
       !!-----------------------------------------------------
       !call Assemble_field_from_VarInt(list_mat,list_varInt_Nloc,Nnloc,Nloc(i)%Var)


       ! Appel du modele correspondant  : CHAINES TESTS EN MINISCULE (voir plus haut)
       !!---------------------------------
       if (before_after .EQ. "after") then
        select case(trim(Nloc_models(i)%Modelname))
        case("test_nloc")
           call Test_nloc(Nnloc,NGnloc,list_mat,list_VarInt_Nloc,list_VarInt_GNloc,i,&
                                 load_n,load_incr, ind_tps)
        ! possibilite de faire co-exister 3 modeles "non-locaux" (en general on en a un seul) 
        case("user_nloc1")
           call user_nloc_after1(Nnloc,NGnloc,list_mat,list_VarInt_Nloc,list_VarInt_GNloc,i,&
                                 load_n,load_incr, ind_tps)
        case("user_nloc2")
           call user_nloc_after2(Nnloc,NGnloc,list_mat,list_VarInt_Nloc,list_VarInt_GNloc,i,&
                                 load_n,load_incr, ind_tps)
        case("user_nloc3")
           call user_nloc_after3(Nnloc,NGnloc,list_mat,list_VarInt_Nloc,list_VarInt_GNloc,i,&
                                 load_n,load_incr, ind_tps)
        case default
           call amitex_abort("Modele non local renseigne non reconnu (Nloc_call)",2,0)
        end select
       else if (before_after .EQ. "befor") then
        select case(trim(Nloc_models(i)%Modelname))
        ! pas de before pour ces trois modeles
        case("test_nloc")
        ! possibilite de faire co-exister 3 modeles "non-locaux" (en general on en a un seul) 
        case("user_nloc1")
           call user_nloc_before1(Nnloc,NGnloc,list_mat,list_VarInt_Nloc,list_VarInt_GNloc,i,&
                                  load_n,load_incr, ind_tps)
        case("user_nloc2")
           call user_nloc_before2(Nnloc,NGnloc,list_mat,list_VarInt_Nloc,list_VarInt_GNloc,i,&
                                  load_n,load_incr, ind_tps)
        case("user_nloc3")
           call user_nloc_before3(Nnloc,NGnloc,list_mat,list_VarInt_Nloc,list_VarInt_GNloc,i,&
                                  load_n,load_incr, ind_tps)
        case default
           call amitex_abort("Modele non local renseigne non reconnu (Nloc_call)",2,0)
        end select
       else
        call amitex_abort(&
              "Variable before_after is not ""after"" neither ""befor"" (Nloc_call)",2,0)
       end if

       ! Repartition des variables internes non locales dans la structure mattotP
       !!----------------------------------------------------------
       !call Assign_VarInt_from_Field(list_mat,list_varInt_GNloc,NGnloc,GNloc(i)%Var)

       deallocate(list_varInt_Nloc)
       deallocate(list_varInt_GNloc)
       deallocate(list_mat)
    end if
    end do

end subroutine Nloc_call

!==================================================================================
!==================================================================================
!                             MODELES NON LOCAUX 
!
!> Les routines qui suivent constituent l'implementation de modeles non locaux specifiques
!  comprenant la specification des variables internes concernees par les calculs non locaux
!  et le calcul proprement dit des variables non locales 
!
! ATTENTION : les calculs font appel a MattotP()%VarInt et MattotP()%Coeff
!             le modele doit donc etre en coherence avec les lois de comportement
!             utilisees dans l'ordonnancement des variables internes et des 
!             coefficients
!
!       FORMAT ENTREES POUR LES MODELES NON LOCAUX
!      ---------------------------------------------------
!             Rq : Ce formalisme permet d'appeler tous les modeles depuis Nloc_call 
!                  sur le modele des appels UMAT
!                  Nloc_call met en forme les entrees des modeles non locaux 
!                  au format presente ci-dessous, a partir des informations
!                  fournies dans les fichiers .xml. Aucun travail de mise en forme
!                  des entrees/sorties n'est necessaire pour le developpement d'un modele
!                  non local
!                  Chaque utilisateur dispose donc des variables d'entree suivantes 
!                  pour developper son modele non local. 
!
!      MEMO : Lors du dvp d'un modele non local, penser a implementer son appel dans Nloc_call 
!     -----   en l'associant a un nom precis
!
!      Entrees des modeles non locaux mises en forme par Nloc_call
!      -----------------------------------------------------------
!               - Nnloc  : nombre de composantes du champ NLOC a construire
!                 ------   c'est a dire le nombre de champs scalaires qui sont 
!                          utilises par le modele pour construire les champs non locaux
!                          (entier x1) 
!
!               - NGnloc : nombre de composantes du champ GNLOC a calculer
!                 ------   c'est a dire le nombre de champs scalaires qui sont 
!                          calcules par le modele 
!                          (entier x1) 
!
!               - list_mat : tableaux d'entiers de taille (nmat)
!                 --------   Contient la liste des numM des materiaux concernes par le modele non local
!                            i.e. les materiaux dont les %VarInt servent a construire NLOC
!                            et/ou contiendront les valeurs du champ GNLOC
!
!               - list_varInt_Nloc : tableaux d'entiers de taille (nmat,Nnloc)
!                 -----------------  Contient ligne par ligne la liste des indices 
!                                    des variables internes de chaque materiau a recuperer
!                                    pour construire le champ NLOC
!                    list_varInt_Nloc(i,:) = ensemble des indices du materiau list_mat(i)
!                                            a recuperer pour construire NLOC
!
!               - list_varInt_GNloc : tableaux d'entiers de taille (nmat,Nnloc)
!                 -----------------   Contient ligne par ligne la liste des indices 
!                                     des variables internes non locales de chaque materiau
!                                     a recuperer depuis la champ GNLOC
!                    list_varInt_GNloc(i,:) = ensemble des indices du materiau list_mat(i)
!                                             ou stocker les valeurs du champ GNLOC
!
!                - modelInd : entier 
!                 ----------  indice du modele non local integre permettant de selectionner
!                             les bons elements des structures NLOC et GNLOC
!                             Dans le modele non local, chaque appel a ces champs doit 
!                             etre sous la forme |  NLOC(modelInd)%... 
!                                                | GNLOC(modelInd)%...
!                             La coherence des tailles et des donnes est assure par 
!                             Nloc_call.
!
!      Principe general de developpement  (voir user_GURTINSSnloc3)
!      ---------------------------------------------------
!               1- Transfert : MattotP%Varint -> champ de variables locales (espace reel)
!                          ==>   call Assemble_field_from_VarInt(list_mat,list_varInt_Nloc,Nnloc,Nloc(modelInd)%Var)
!
!               2- NLOC : si necessaire, changement de base repere local-> global 
! 
!               3- On passe le champ NLOC dans l'espace de Fourier
!                          ==>   call field_fft(Nloc(modelInd)%Var,Nloc(modelInd)%VarF,Nnloc,1)
!
!               4- On calcule dans l'espace de Fourier le champ de variables GNLOC
!                  non locales, notamment a l'aide des operateurs de derivation
!                           DEVELOPPEMENT PROPREMENT DIT DES CALCULS DU MODELE NON LOCAL ICI
!                          ==> on construit ainsi GNLOC(modelInd)%VarF
!
!               5- On repasse  GNLOC dans l'espace reel 
!                          ==> call field_ifft(GNloc(modelInd)%Var(...),GNloc(modelInd)%VarF(...),NGnloc,1)  
!                              NE PAS OUBLIER LA RENORMALISATION :   
!                              GNloc(modelInd)%VarF = GNloc(modelInd)%VarF / real(grid%ntot,mytype
!
!               6- GNLOC : si necessaire, changement de base repere local-> global
!
!               7- Transfert : champs de variables non-locales -> MattotP%Varint (espace reel)
!                          En general:
!                          ==>  call Assign_VarInt_from_Field(list_mat,list_varInt_GNloc,NGnloc,GNloc(modelInd)%Var)
!
!      Remarque : Des calculs intermediaires pourraient etre realises entre les differentes etapes selon les modeles 
!
!      Remarque : operateurs de derivation disponibles dans le field_mod
!
!      Remarque : si la fonction est très simple, il peut etre inutile de proceder aux transferts 
!                 de Varint vers NLOC puis de GNloc vers Varint 
!
!==================================================================================
!==================================================================================
 
!==================================================================================
!                             SUBROUTINE Test_nloc
!
!> Test algorithme non local : utilise pour tester et valider les operateurs de 
!                              derivation
!!
!==================================================================================
subroutine Test_nloc(Nnloc,NGnloc,list_mat,list_varInt_Nloc,list_varInt_GNloc,modelInd,&
                                 load_n,load_incr, ind_tps)

 implicit none
 !> Entrees 
 integer,intent(in)                   :: Nnloc, NGnloc
 integer,dimension(:), intent(in)     :: list_mat
 integer,dimension(:,:),intent(in)    :: list_varInt_Nloc
 integer,dimension(:,:),intent(in)    :: list_varInt_GNloc
 integer,intent(in)                   :: modelInd
 integer,intent(in)                   :: load_n,load_incr, ind_tps

 ! Autre variables
 integer,dimension(3)                 :: N ! dimensions de la cellule 
 real(mytype),dimension(3)            :: d ! dimensions des voxels 

 ! trick to avoid gcc-warning
 integer                              :: bidon
 bidon = Nnloc + NGnloc + list_mat(1) + list_varInt_Nloc(1,1) + list_varInt_GNloc(1,1) + modelInd &
       + load_n + load_incr + ind_tps

 N(1) = grid%nx
 N(2) = grid%ny
 N(3) = grid%nz
 d(1) = grid%dx
 d(2) = grid%dy
 d(3) = grid%dz

 ! Construction du champ de variables locales a deriver 
 ! choix ici : valeur nulle là ou le modele non-local n'est pas defini
 !!-------------------------------------------------------------------
 Nloc(modelInd)%Var=0.
 call Assemble_field_from_VarInt(list_mat,list_varInt_Nloc,Nnloc,Nloc(modelInd)%Var)

 ! II - Passage dans l'espace de Fourier 
 !!----------------------------------------------------------
 call field_fft(Nloc(modelInd)%Var,Nloc(modelInd)%VarF,Nnloc,1)

 ! III - Derivation
 !!----------------------------------------------------------
 !call laplace_scalar(Nloc(modelInd)%VarF(:,:,:,1),GNloc(modelInd)%VarF(:,:,:,1),N,d)
 !call laplace_scalar(Nloc(modelInd)%VarF(:,:,:,2),GNloc(modelInd)%VarF(:,:,:,2),N,d,FREQ)

 ! calcul de gamma,2 --> gamma = Nloc(modelInd)%VarF(:,:,:,6)
 call field_first_order_partial_centeredF(GNloc(modelInd)%VarF(:,:,:,10),1._mytype,Nloc(modelInd)%VarF(:,:,:,6),2,1,1,N(2),d(2))
 ! calcul du rot(Hp)
 call rotF(Nloc(modelInd)%VarF,GNloc(modelInd)%VarF(:,:,:,1:9),N,d) 

 ! IV - Retour dans l'espace Reel pour le champ "gradient"
 !!----------------------------------------------------------
 GNloc(modelInd)%VarF = GNloc(modelInd)%VarF / real(grid%ntot,mytype) ! renormalisation
 call field_ifft(GNloc(modelInd)%Var,GNloc(modelInd)%VarF,NGnloc,1)
 ! calcul de la norme 
 call field_normcomp(GNloc(modelInd)%Var(:,1:9),GNloc(modelInd)%Var(:,11),9,1)
 ! calcul de la norme 
 call field_normcomp(Nloc(modelInd)%Var(:,1:9),GNloc(modelInd)%Var(:,12),9,1)

 ! V - Repartition des variables internes non locales dans la structure mattotP
 !!----------------------------------------------------------
 call Assign_VarInt_from_Field(list_mat,list_varInt_GNloc,NGnloc,GNloc(modelInd)%Var)

end subroutine Test_nloc

end module non_local_mod
