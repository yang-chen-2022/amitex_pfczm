!===============================================================================
!
! MODULE NON_LOCAL_USER_MOD : 
!> User-defined non-local models (white module - to be implemented by the user)
!!
!!           GINZBURG-LANDAU PHASE FIELD FOR GRAIN-GROWTH
!! 
!===============================================================================
module non_local_user_mod

  use ISO_FORTRAN_ENV

! MPI AND 2DECOMP MODULES
!------------------------
  use MPI
#ifdef NEW_2DECOMP
  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
#else
  use decomp_2d 
#endif  
  use decomp_2d_fft

! AMITEX MODULES 
!----------------
  use error_mod
  use amitex_mod       !, only : Nloc,GNloc  ! non-local fields, access to other fields
  use param_algo_mod   !, only : algo_param
  use material_mod     !, only : Assign_VarInt_from_Field,Assemble_field_from_VarInt,&
                       !         get_homogeneous_coeff,&
                       !         mattotP,calcrotmatrices_GD, Nloc_models
  use field_mod        !, only : fft_start, fft_end, field_fft, field_ifft 
  use green_mod        !, only : FREQ, grid
  use loading_mod      !, only : load,extract
  use io2_amitex_mod   !, only : print_field_vtk

  implicit none

  private

  !> public subroutines (one per non-local model) - DO NOT TOUCH
  public :: user_nloc_after1,user_nloc_after2,user_nloc_after3,&
            user_nloc_before1, user_nloc_before2, user_nloc_before3,&
            init_user_nloc_variables1,init_user_nloc_variables2,init_user_nloc_variables3

  !> USER variables  
  real(mytype),allocatable,dimension(:)     :: alpha, kinetic   ! coefficients (homogenes) des champ de phase (one/phase field)
  real(mytype),allocatable,dimension(:,:,:) :: FREQ2            ! Frequences au carre
  real(mytype),allocatable,dimension(:)     :: Etasum           ! sum of phase fields weighted by indices (for output)
  
  real(mytype),pointer,dimension(:,:)         :: Eta            ! ALIASES
  real(mytype),pointer,dimension(:,:)         :: df_dEta        ! Associated at the begining of user_nloc_after1
  complex(mytype),pointer,dimension(:,:,:,:)  :: EtaF           ! (last dimension is the phase index)
  

contains

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
!      Principe de developpement  (voir GURTIN_SS)
!      --------------------------
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
!                          ==>  call Assign_VarInt_from_Field(list_mat,list_varInt_GNloc,NGnloc,GNloc(modelInd)%Var)
!
!      Remarque : Des calculs intermediaires pourraient etre realises entre les differentes etapes selon les modeles 
!
!      Remarque : operateurs de derivation disponibles dans le field_mod
!
!      Remarque : si la fonction est très simple, il peut etre inutile de proceder aux transferts 
!                 de Varint vers NLOC puis de GNloc vers Varint 
!
!
!==================================================================================

!==================================================================================
!      GRAINGROWTH3 : INIT_USER_NLOC_VARIABLES1 &  USER_NLOC_AFTER1
!
!> Implementation of a Gindzburg-Laudau phase-field for grain-growth
!! 
!!            with the behavior graingrowth_10grains_3.f90
!!
!!                 See EXAMPLES in cas_tests_PF 
!!
!!
!!     COEFF : 
!!          alpha(1)   : 1st non-local coeff
!!          kinetic(1) : 2nd non-local coeff
!!          alpha(2)   : 3rd non-local coeff
!!          kinetic(2) : 4th non-local coeff
!!          ...
!!
!!     NLOC : 
!!          df_dEta : Nloc(modelInd)%Var(:,:), last index for phase field number
!!
!!     GNLOC : 
!!          Eta     : GNloc(modelInd)%Var(:,:), last index for phase field number
!!          EtaF    : GNloc(modelInd)%VarF(:,:,:,:), last index for phase field number
!!
!! LIMITATIONS : HOMOGENEOUS PHASE-FIELD COEFF (alpha and kinetic)
!!
!==================================================================================

subroutine init_user_nloc_variables1(Nnloc,NGnloc,list_mat,list_varInt_Nloc,list_varInt_GNloc,modelInd)

 implicit none
 !> Entrees (NE PAS TOUCHER)
 integer,intent(in)                   :: Nnloc, NGnloc
 integer,dimension(:), intent(in)     :: list_mat
 integer,dimension(:,:),intent(in)    :: list_varInt_Nloc
 integer,dimension(:,:),intent(in)    :: list_varInt_GNloc
 integer,intent(in)                   :: modelInd
 
 !> autres variables
 integer                              :: alloc_stat
 integer                              :: i
 integer                              :: Nphases     ! nb de champ de phases

 ! trick to avoid gcc-warning
 integer                              :: bidon
 bidon = Nnloc + NGnloc + list_mat(1) + list_varInt_Nloc(1,1) + list_varInt_GNloc(1,1) + modelInd 

 !================================================================================
 !                                                  INITIALIZE FIELDS FROM VARINT0
 !                                               (initial values in mate.xml file)

 call Assemble_field_from_VarInt0(list_mat,list_varInt_GNloc,NGnloc,GNloc(modelInd)%Var) !-> Eta
 call Assemble_field_from_VarInt0(list_mat,list_varInt_Nloc,Nnloc,Nloc(modelInd)%Var)    !-> df_dEta

 !================================================================================
 !                                                    INITIALIZE HOMOGENEOUS COEFF
 !                                                       (values in mate.xml file)
 !  alpha(i)   : 1st non-local coeff for phase i
 !  kinetic(i) : 2nd non-local coeff for phase i
 !  get_homogeneous_coeff verifie si les coeff entres sont bien homogenes

 Nphases = NGnloc 

 allocate(alpha(Nphases))  ! nb coeff = nb champ de phases
 allocate(kinetic(Nphases))
 
 do i = 1, Nphases
    call get_homogeneous_coeff(alpha(i),  modelInd,  2*i-1, list_mat)
    call get_homogeneous_coeff(kinetic(i), modelInd, 2*i  , list_mat)
 end do

 !================================================================================
 !                                                 OTHER ALLOCATION/INITIALIZATION
 !                                               

 allocate(FREQ2(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3)),stat=alloc_stat)
 if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_nloc_variables1 1)",2)

 FREQ2 = FREQ(:,:,:,1)**2 + FREQ(:,:,:,2)**2 + FREQ(:,:,:,3)**2


 allocate(Etasum(lbound(GNloc(modelInd)%Var,1):ubound(GNloc(modelInd)%Var,1)),stat=alloc_stat)
 if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_nloc_variables1 1)",2)

 Etasum = 0.

end subroutine init_user_nloc_variables1


!==================================================================================
!      GRAINGROWTH3 : INIT_USER_NLOC_VARIABLES1 & USER_NLOC_AFTER1
!
!> Implementation of a Gindzburg-Laudau phase-field for grain-growth
!! 
!!            with the behavior graingrowth_10grains_3.f90
!!
!!                 See EXAMPLES in cas_tests_PF 
!!
!!
!!     COEFF : 
!!          alpha(1)   : 1st non-local coeff
!!          kinetic(1) : 2nd non-local coeff
!!          alpha(2)   : 3rd non-local coeff
!!          kinetic(2) : 4th non-local coeff
!!          ...
!!
!!     NLOC : 
!!          df_dEta : Nloc(modelInd)%Var(:,:), last index for phase field number
!!
!!     GNLOC : 
!!          Eta     : GNloc(modelInd)%Var(:,:), last index for phase field number
!!          EtaF    : GNloc(modelInd)%VarF(:,:,:,:), last index for phase field number
!!
!! LIMITATIONS : HOMOGENEOUS PHASE-FIELD COEFF (alpha and kinetic)
!!
!==================================================================================

subroutine user_nloc_after1(Nnloc,NGnloc,list_mat,list_varInt_Nloc,list_varInt_GNloc,modelInd,&
                            load_n,load_incr, ind_tps)

 implicit none
 !> Entrees (DO NOT TOUCH)
 integer,intent(in)                   :: Nnloc, NGnloc
 integer,dimension(:), intent(in)     :: list_mat
 integer,dimension(:,:),intent(in)    :: list_varInt_Nloc
 integer,dimension(:,:),intent(in)    :: list_varInt_GNloc
 integer,intent(in)                   :: modelInd
 integer,intent(in)                   :: load_n,load_incr, ind_tps

 ! Autres variables (USER variables)
 real(mytype)                         :: dt
 integer                              :: i
 integer                              :: Nphases
 character(len=200)                   :: tmp_char1,tmp_char   !< variable caractere temporaire


 ! trick to avoid gcc-warning
 integer                              :: bidon
 bidon = Nnloc + NGnloc + list_mat(1) + list_varInt_Nloc(1,1) + list_varInt_GNloc(1,1) + modelInd &
       + load_n + load_incr + ind_tps
       
 !================================================================================
 !                                                   POINTER ASSOCIATION (ALIASES)
 !
 df_dEta => Nloc(modelInd)%Var
 Eta     => GNloc(modelInd)%Var
 EtaF    => GNloc(modelInd)%VarF
 
 !================================================================================
 !                                                                 INITIALISATIONS
 !
 !---- pas de temps
 if (load_incr==0) then
    dt = load(load_n)%time(load_incr)
 else
    dt = load(load_n)%time(load_incr) -  load(load_n)%time(load_incr-1)
 end if
 
 !---- Number of phase fields
 Nphases = NGnloc 
 
 !================================================================================
 !                                                   EVOLUTION DES CHAMPS DE PHASE

 ! GET df/dEta from MattotP%Varint (via Nloc%Var)
 ! GET Eta from MattotP%Varint (via GNloc%Var)
 !--------------------------------------------------
 call Assemble_field_from_VarInt(list_mat,list_varInt_Nloc,Nnloc,Nloc(modelInd)%Var)     !-> df_deta (pointer on Nloc(modelInd)%Var)
 call Assemble_field_from_VarInt(list_mat,list_varInt_GNloc,NGnloc,GNloc(modelInd)%Var)  !-> Eta  (pointer on GNloc(modelInd)%Var)

 do i=1,Nphases
  
    ! (Eta - L.dt.df_dEta) -> Eta
    !--------------------------------------------------
    Eta(:,i) = Eta(:,i) - Kinetic(i)*dt*df_dEta(:,i)

    ! FOURIER TRANSFORM -> EtaF
    !--------------------------------------------------
    call field_fft(Eta(:,i),EtaF(:,:,:,i),1,1) 
    EtaF(:,:,:,i)=EtaF(:,:,:,i)/ real(grid%ntot,mytype)

    ! SEMI-IMPLICIT SCHEME
    !--------------------------------------------------
    EtaF(:,:,:,i) = EtaF(:,:,:,i)/(1._mytype + FREQ2(:,:,:) * Kinetic(i) * alpha(i) * dt)

    ! INVERSE FOURIER TRANSFORM                  -> Eta
    !--------------------------------------------------
    call field_ifft(Eta(:,i),EtaF(:,:,:,i),1,1)

 end do

 ! STORE "FINAL" Eta in MattotP%Varint (via GNloc%Var)
 !--------------------------------------------------
 call Assign_VarInt_from_Field(list_mat,list_varInt_GNloc,NGnloc,GNloc(modelInd)%Var) !<-Eta (pointer on GNloc(modelInd)%Var(:,1))


 !================================================================================
 !                                                                      SORTIE VTK
 !            ci-dessous : ind_tps optional -> ecriture si extract%tpsVTK(ind_tps)
  
 if(load_incr/=0) then ! Pas de sortie vtk pour les pas de temps fictifs
 
     Etasum = 0.
     do i = 1,Nphases
        Etasum = Etasum + Eta(:,i)*real(i,mytype)
     end  do

     write(tmp_char1,"(I6)") ind_tps
     tmp_char="_etasum_"//trim(adjustl(tmp_char1))   
     call print_field_vtk(Etasum,trim(fic_vtk)//trim(tmp_char)//".vtk","etasum", ind_tps) 

 end if


end  subroutine user_nloc_after1


!##################################################################################
!##################################################################################
!#########################################################          white functions
!##################################################################################
!##################################################################################
!##########################################          do not modify if not necessary
!##################################################################################
!##################################################################################
subroutine user_nloc_after2(Nnloc,NGnloc,list_mat,list_varInt_Nloc,list_varInt_GNloc,modelInd,&
                            load_n,load_incr, ind_tps)

 implicit none
 !> Entrees 
 integer,intent(in)                   :: Nnloc, NGnloc
 integer,dimension(:), intent(in)     :: list_mat
 integer,dimension(:,:),intent(in)    :: list_varInt_Nloc
 integer,dimension(:,:),intent(in)    :: list_varInt_GNloc
 integer,intent(in)                   :: modelInd
 integer,intent(in)                   :: load_n,load_incr, ind_tps

 ! trick to avoid gcc-warning
 integer                              :: bidon
 bidon = Nnloc + NGnloc + list_mat(1) + list_varInt_Nloc(1,1) + list_varInt_GNloc(1,1) + modelInd &
       + load_n + load_incr + ind_tps

end subroutine user_nloc_after2

!----------------------------------------------------------------------------------
subroutine user_nloc_after3(Nnloc,NGnloc,list_mat,list_varInt_Nloc,list_varInt_GNloc,modelInd,&
                            load_n,load_incr, ind_tps)

 implicit none
 !> Entrees 
 integer,intent(in)                   :: Nnloc, NGnloc
 integer,dimension(:), intent(in)     :: list_mat
 integer,dimension(:,:),intent(in)    :: list_varInt_Nloc
 integer,dimension(:,:),intent(in)    :: list_varInt_GNloc
 integer,intent(in)                   :: modelInd
 integer,intent(in)                   :: load_n,load_incr, ind_tps

 ! trick to avoid gcc-warning
 integer                              :: bidon
 bidon = Nnloc + NGnloc + list_mat(1) + list_varInt_Nloc(1,1) + list_varInt_GNloc(1,1) + modelInd &
       + load_n + load_incr + ind_tps

end subroutine user_nloc_after3

!----------------------------------------------------------------------------------
subroutine user_nloc_before1(Nnloc,NGnloc,list_mat,list_varInt_Nloc,list_varInt_GNloc,modelInd,&
                            load_n,load_incr, ind_tps)

 implicit none
 !> Entrees 
 integer,intent(in)                   :: Nnloc, NGnloc
 integer,dimension(:), intent(in)     :: list_mat
 integer,dimension(:,:),intent(in)    :: list_varInt_Nloc
 integer,dimension(:,:),intent(in)    :: list_varInt_GNloc
 integer,intent(in)                   :: modelInd
 integer,intent(in)                   :: load_n,load_incr, ind_tps

 ! trick to avoid gcc-warning
 integer                              :: bidon
 bidon = Nnloc + NGnloc + list_mat(1) + list_varInt_Nloc(1,1) + list_varInt_GNloc(1,1) + modelInd &
       + load_n + load_incr + ind_tps

end subroutine user_nloc_before1

!----------------------------------------------------------------------------------
subroutine user_nloc_before2(Nnloc,NGnloc,list_mat,list_varInt_Nloc,list_varInt_GNloc,modelInd,&
                            load_n,load_incr, ind_tps)

 implicit none
 !> Entrees 
 integer,intent(in)                   :: Nnloc, NGnloc
 integer,dimension(:), intent(in)     :: list_mat
 integer,dimension(:,:),intent(in)    :: list_varInt_Nloc
 integer,dimension(:,:),intent(in)    :: list_varInt_GNloc
 integer,intent(in)                   :: modelInd
 integer,intent(in)                   :: load_n,load_incr, ind_tps

 ! trick to avoid gcc-warning
 integer                              :: bidon
 bidon = Nnloc + NGnloc + list_mat(1) + list_varInt_Nloc(1,1) + list_varInt_GNloc(1,1) + modelInd &
       + load_n + load_incr + ind_tps

end subroutine user_nloc_before2

!----------------------------------------------------------------------------------
subroutine user_nloc_before3(Nnloc,NGnloc,list_mat,list_varInt_Nloc,list_varInt_GNloc,modelInd,&
                            load_n,load_incr, ind_tps)

 implicit none
 !> Entrees 
 integer,intent(in)                   :: Nnloc, NGnloc
 integer,dimension(:), intent(in)     :: list_mat
 integer,dimension(:,:),intent(in)    :: list_varInt_Nloc
 integer,dimension(:,:),intent(in)    :: list_varInt_GNloc
 integer,intent(in)                   :: modelInd
 integer,intent(in)                   :: load_n,load_incr, ind_tps

 ! trick to avoid gcc-warning
 integer                              :: bidon
 bidon = Nnloc + NGnloc + list_mat(1) + list_varInt_Nloc(1,1) + list_varInt_GNloc(1,1) + modelInd &
       + load_n + load_incr + ind_tps

end subroutine user_nloc_before3

!----------------------------------------------------------------------------------
subroutine init_user_nloc_variables2(Nnloc,NGnloc,list_mat,list_varInt_Nloc,list_varInt_GNloc,modelInd)

 implicit none
 !> Entrees 
 integer,intent(in)                   :: Nnloc, NGnloc
 integer,dimension(:), intent(in)     :: list_mat
 integer,dimension(:,:),intent(in)    :: list_varInt_Nloc
 integer,dimension(:,:),intent(in)    :: list_varInt_GNloc
 integer,intent(in)                   :: modelInd

 ! trick to avoid gcc-warning
 integer                              :: bidon
 bidon = Nnloc + NGnloc + list_mat(1) + list_varInt_Nloc(1,1) + list_varInt_GNloc(1,1) + modelInd 

end subroutine init_user_nloc_variables2

!----------------------------------------------------------------------------------
subroutine init_user_nloc_variables3(Nnloc,NGnloc,list_mat,list_varInt_Nloc,list_varInt_GNloc,modelInd)

 implicit none
 !> Entrees 
 integer,intent(in)                   :: Nnloc, NGnloc
 integer,dimension(:), intent(in)     :: list_mat
 integer,dimension(:,:),intent(in)    :: list_varInt_Nloc
 integer,dimension(:,:),intent(in)    :: list_varInt_GNloc
 integer,intent(in)                   :: modelInd

 ! trick to avoid gcc-warning
 integer                              :: bidon
 bidon = Nnloc + NGnloc + list_mat(1) + list_varInt_Nloc(1,1) + list_varInt_GNloc(1,1) + modelInd 

end subroutine init_user_nloc_variables3



end module non_local_user_mod
