!===============================================================================
!
!> MODULE NON_LOCAL_USER_MOD : IMPLEMENTATION OF MIEHE PHASE FIELD IN THE AMITEX
!>                             NON LOCAL FRAMEWORK
!! 
!!
!!  Y.Chen & L.Gelebart - 2018-2019
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

! AMITEX MODULES (all available)
!-------------------------------
  use error_mod
  use amitex_mod      
  use amitex_mod      !, only : Nloc,GNloc  non local fields
  use param_algo_mod  !, only : Nloc_models
  use material_mod    !, only : Assemble_field_from_VarInt0, mattotP
                           
  use field_mod
  use green_mod
  use loading_mod     !, only : load

! LOCAL MODULES
!--------------
  use Miehe_mod

  implicit none

  private

  !> public subroutines (one per non-local model)
  public :: user_nloc_after1, user_nloc_after2, user_nloc_after3,&
            user_nloc_before1, user_nloc_before2, user_nloc_before3,&
            init_user_nloc_variables1,init_user_nloc_variables2,init_user_nloc_variables3

contains

!===================================================================================================
!> init_user_nloc_variables1 - For Phase Field Miehe
!!     initialize NLOC variables (ie d and HH) with mattotP%Varint0
!!     initialize square frequencies with ... "TO VERIFY"
!!     allocate additional fields for resolPF_visc and initilalize them to 0
!!
!===================================================================================================

subroutine init_user_nloc_variables1(Nnloc,NGnloc,list_mat,list_varInt_Nloc,list_varInt_GNloc,modelInd)

  implicit none

  !> Input 
  integer,intent(in)                   :: Nnloc, NGnloc
  integer,dimension(:), intent(in)     :: list_mat
  integer,dimension(:,:),intent(in)    :: list_varInt_Nloc
  integer,dimension(:,:),intent(in)    :: list_varInt_GNloc
  integer,intent(in)                   :: modelInd

  !> Local Variables
  integer                              :: alloc_stat

  ! trick to avoid gcc-warning
  integer                              :: bidon
  bidon = Nnloc + NGnloc + list_mat(1) + list_varInt_Nloc(1,1) + list_varInt_GNloc(1,1) + modelInd 

  !!------------------------------------------------------------------------------
  !>                           INITIALIZATION OF NLOC FIELDS BY INTERNAL VARIABLES
  !>                                           Nloc variables 1 and 2 <-> d and HH

  ! Initialize Nloc fields with corresponding internal variables
  call Assemble_field_from_VarInt0(list_mat,list_varInt_Nloc,Nnloc,Nloc(modelInd)%Var)

  ! polarization term for damage problem
  allocate(tau(xsize(1)*xsize(2)*xsize(3),1),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_nloc_variables1 3)",2)
  tau=0.
  allocate(tau0(xsize(1)*xsize(2)*xsize(3),1),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_nloc_variables1 3)",2)
  tau0=0.

  ! in Fourier space
  allocate(tauF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1),&
           stat=alloc_stat)   
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_nloc_variables1 4)",2)
  tauF=0.
 
  ! process variables for computing the polarization term tau
  allocate(AA(xsize(1)*xsize(2)*xsize(3),1),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_nloc_variables1 5)",2)
  AA=0.
  allocate(BB(xsize(1)*xsize(2)*xsize(3),1),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_nloc_variables1 5)",2)
  BB=0.
  
  ! storage variables for convergence acceleration of phase-field solution
     if (nloc_param(modelInd)%p_string(1)=="true") then
     allocate(ACT3_Rdamage(xsize(1)*xsize(2)*xsize(3),1,4),stat=alloc_stat)
     allocate(ACT3_Udamage(xsize(1)*xsize(2)*xsize(3),1,4),stat=alloc_stat)
     if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_nloc_variables1 8)",2)
     ACT3_Rdamage=0.
     ACT3_Udamage=0.
  end if

  ! frequencies for Laplacian operator
if (1==2) then !desactivation
  allocate(FreqLaplacian(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_nloc_variables1 6)",2)
  FreqLaplacian=0.
   call initFreqLaplacian(FreqLaplacian)
end if

end subroutine init_user_nloc_variables1


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
!               => FAIT DANS NLOC_CALL (a partir du 01/07/2019)
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
!               => FAIT DANS NLOC_CALL(a partir du 01/07/2019)
!
!      Remarque : Des calculs intermediaires pourraient etre realises entre les differentes etapes selon les modeles 
!
!      Remarque : operateurs de derivation disponibles dans le field_mod
!
!==================================================================================


!==================================================================================
!> user_nloc_before1 - For Phase Field Miehe
!>                     Evaluate Damage Phase Field
!
!==================================================================================
subroutine user_nloc_before1(Nnloc,NGnloc,list_mat,list_varInt_Nloc,list_varInt_GNloc,modelInd,&
                             load_n,load_incr, ind_tps)

 implicit none

 !> Input 
 integer,intent(in)                   :: Nnloc, NGnloc
 integer,dimension(:), intent(in)     :: list_mat
 integer,dimension(:,:),intent(in)    :: list_varInt_Nloc
 integer,dimension(:,:),intent(in)    :: list_varInt_GNloc
 integer,intent(in)                   :: modelInd
 integer,intent(in)                   :: load_n,load_incr, ind_tps

 !> Local variables
 real(mytype)                         :: dt, dt_
 real(mytype)                         :: timePF
 real(mytype)                         :: critPF       !criterion for phase-field problem
 integer                              :: nItPF        !number of iterations for phase field problem
 integer,save                         :: nITPFcum

 ! trick to avoid gcc-warning
 integer                              :: bidon
 bidon = Nnloc + NGnloc + list_mat(1) + list_varInt_Nloc(1,1) + list_varInt_GNloc(1,1) + modelInd &
       + load_n + load_incr + ind_tps

 !-----------------------------------------------------------------------------

 !-- Evaluate time steps dt and dt_
 if (load_incr==0) then
    dt = load(load_n)%time(load_incr)
    dt_= 0.
 elseif (load_incr==1) then
    dt  = load(load_n)%time(load_incr)   -  load(load_n)%time(load_incr-1)
    dt_ = load(load_n)%time(load_incr)    
 else
    dt  = load(load_n)%time(load_incr)   -  load(load_n)%time(load_incr-1)
    dt_ = load(load_n)%time(load_incr-1) -  load(load_n)%time(load_incr-2)
 end if


 !! solve the phase-field problem (evaluate HH and damageVar, global variables)
 !-----------------------------------------------------------------------------
 call resolPF_visc(load_incr,dt,dt_,nItPF,critPF,timePF,list_mat,modelInd)


 !!                    update the damage field and internal variables (of UMAT)
 !-----------------------------------------------------------------------------
 call Assign_VarInt0_from_Field(list_mat,list_VarInt_Nloc,Nnloc,Nloc(modelInd)%Var)   


 !!                       save the initial fields (d and H), in Varint, in Nloc
 !                               for use in the next step as the (i-1) solution 
 !                        to predictict an initial damage field in resolPF_visc  
 !-----------------------------------------------------------------------------
 if (trim(nloc_param(modelInd)%P_string(2))=="true") then
 if (load_incr > 1) then
   call Assemble_Field_From_Varint(list_mat,list_VarInt_Nloc,Nnloc,Nloc(modelInd)%Var)
 end if
 end if

 !-----------------------------------------------------------------------------
 nItPFcum = nItPFcum + nItPF
 if (nrank==0) write(OUTPUT_UNIT,*) "NItPFcum = ",NItPFcum
 if (nrank==0) write(Flog,*) "NItPFcum = ",NItPFcum

end  subroutine user_nloc_before1


!==================================================================================
!> user_nloc_after1 - For Phase Field Miehe
!>                    Evaluate History Field
!
!==================================================================================
subroutine user_nloc_after1(Nnloc,NGnloc,list_mat,list_varInt_Nloc,list_varInt_GNloc,modelInd,&
                            load_n,load_incr, ind_tps)

 implicit none

 !> Input 
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

 !-----------------------------------------------------------------------------

 !!                                                        Update History Field
 !                                  GNloc = max(GNloc, H dans mattotP()%Varint)
 !-----------------------------------------------------------------------------

 call updateHistory(list_mat,list_varInt_Nloc,modelInd)  

end  subroutine user_nloc_after1



!##################################################################################
!##################################################################################
!#########################################################          white functions
!##################################################################################
!##################################################################################
!##########################################          do not modify if not necessary
!##################################################################################
!##################################################################################
!----------------------------------------------------------------------------------
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

 ! trick to avaoid gcc-warning
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

!==================================================================================
!
! copy-paste of previous functions, replacing user_nloc_before by init_user_nloc_variables
!
!==================================================================================

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
