!===============================================================================
!
! MODULE NON_LOCAL_USER_MOD : 
!> User-defined non-local models (white module - to be implemented by the user)
!! 
!===============================================================================
module non_local_user_mod

  use ISO_FORTRAN_ENV

! MPI AND 2DECOMP MODULES
!------------------------
  use MPI
#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi,       only : nrank
#else
  use decomp_2d, only : nrank,mytype 
#endif  

! AMITEX MODULES (all available)
!-------------------------------
  use error_mod
  use amitex_mod      
  
  implicit none

  private

  !> public subroutines (one per non-local model)
  public :: user_nloc_after1, user_nloc_after2, user_nloc_after3,&
            user_nloc_before1, user_nloc_before2, user_nloc_before3,&
            init_user_nloc_variables1,init_user_nloc_variables2,init_user_nloc_variables3

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

subroutine user_nloc_after1(Nnloc,NGnloc,list_mat,list_varInt_Nloc,list_varInt_GNloc,modelInd,&
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

  if (nrank==0) then
     write(Flog,"(A)") "ERROR :"
     write(Flog,"(A)") "      <Non_local_modeling Modelname= ""user_nloc1"" ....> is defined in xml algorithm file "
     write(Flog,"(A)") "      BUT 'user_nloc_after1' is not implemented "

     write(output_unit,"(A)") "ERROR :"
     write(output_unit,"(A)") "      <Non_local_modeling Modelname= ""user_nloc1"" ....> is defined in xml algorithm file "
     write(output_unit,"(A)") "      BUT 'user_nloc_after1' is not implemented "
  end if

end subroutine user_nloc_after1

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

 ! trick to avoid gcc-warning
 integer                              :: bidon
 bidon = Nnloc + NGnloc + list_mat(1) + list_varInt_Nloc(1,1) + list_varInt_GNloc(1,1) + modelInd &
       + load_n + load_incr + ind_tps

  if (nrank==0) then
     write(Flog,"(A)") "ERROR :"
     write(Flog,"(A)") "      <Non_local_modeling Modelname= ""user_nloc2"" ....> is defined in xml algorithm file "
     write(Flog,"(A)") "      BUT 'user_nloc_after2' is not implemented "

     write(output_unit,"(A)") "ERROR :"
     write(output_unit,"(A)") "      <Non_local_modeling Modelname= ""user_nloc2"" ....> is defined in xml algorithm file "
     write(output_unit,"(A)") "      BUT 'user_nloc_after2' is not implemented "
  end if

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

  if (nrank==0) then
     write(Flog,"(A)") "ERROR :"
     write(Flog,"(A)") "      <Non_local_modeling Modelname= ""user_nloc3"" ....> is defined in xml algorithm file "
     write(Flog,"(A)") "      BUT 'user_nloc_after3' is not implemented "

     write(output_unit,"(A)") "ERROR :"
     write(output_unit,"(A)") "      <Non_local_modeling Modelname= ""user_nloc3"" ....> is defined in xml algorithm file "
     write(output_unit,"(A)") "      BUT 'user_nloc_after3' is not implemented "
  end if

end subroutine user_nloc_after3


!==================================================================================
!
! copy-paste of previous functions, replacing after by before
!
!==================================================================================

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

  if (nrank==0) then
     write(Flog,"(A)") "ERROR :"
     write(Flog,"(A)") "      <Non_local_modeling Modelname= ""user_nloc1"" ....> is defined in xml algorithm file "
     write(Flog,"(A)") "      BUT 'user_nloc_before1' is not implemented "

     write(output_unit,"(A)") "ERROR :"
     write(output_unit,"(A)") "      <Non_local_modeling Modelname= ""user_nloc1"" ....> is defined in xml algorithm file "
     write(output_unit,"(A)") "      BUT 'user_nloc_before1' is not implemented "
  end if

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

  if (nrank==0) then
     write(Flog,"(A)") "ERROR :"
     write(Flog,"(A)") "      <Non_local_modeling Modelname= ""user_nloc2"" ....> is defined in xml algorithm file "
     write(Flog,"(A)") "      BUT 'user_nloc_before2' is not implemented "

     write(output_unit,"(A)") "ERROR :"
     write(output_unit,"(A)") "      <Non_local_modeling Modelname= ""user_nloc2"" ....> is defined in xml algorithm file "
     write(output_unit,"(A)") "      BUT 'user_nloc_before2' is not implemented "
  end if

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

  if (nrank==0) then
     write(Flog,"(A)") "ERROR :"
     write(Flog,"(A)") "      <Non_local_modeling Modelname= ""user_nloc3"" ....> is defined in xml algorithm file "
     write(Flog,"(A)") "      BUT 'user_nloc_before3' is not implemented "

     write(output_unit,"(A)") "ERROR :"
     write(output_unit,"(A)") "      <Non_local_modeling Modelname= ""user_nloc3"" ....> is defined in xml algorithm file "
     write(output_unit,"(A)") "      BUT 'user_nloc_before3' is not implemented "
  end if

end subroutine user_nloc_before3

!==================================================================================
!
! copy-paste of previous functions, replacing user_nloc_before by init_user_nloc_variables
!
!==================================================================================

subroutine init_user_nloc_variables1(Nnloc,NGnloc,list_mat,list_varInt_Nloc,list_varInt_GNloc,modelInd)

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

  if (nrank==0) then
     write(Flog,"(A)") "ERROR :"
     write(Flog,"(A)") "      <Non_local_modeling Modelname= ""user_nloc1"" ....> is defined in xml algorithm file "
     write(Flog,"(A)") "      BUT 'init_user_nloc_variables1' is not implemented "

     write(output_unit,"(A)") "ERROR :"
     write(output_unit,"(A)") "      <Non_local_modeling Modelname= ""user_nloc1"" ....> is defined in xml algorithm file "
     write(output_unit,"(A)") "      BUT 'init_user_nloc_variables1' is not implemented "
  end if

end subroutine init_user_nloc_variables1

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

  if (nrank==0) then
     write(Flog,"(A)") "ERROR :"
     write(Flog,"(A)") "      <Non_local_modeling Modelname= ""user_nloc2"" ....> is defined in xml algorithm file "
     write(Flog,"(A)") "      BUT 'init_user_nloc_variables2' is not implemented "

     write(output_unit,"(A)") "ERROR :"
     write(output_unit,"(A)") "      <Non_local_modeling Modelname= ""user_nloc2"" ....> is defined in xml algorithm file "
     write(output_unit,"(A)") "      BUT 'init_user_nloc_variables2' is not implemented "
  end if

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

  if (nrank==0) then
     write(Flog,"(A)") "ERROR :"
     write(Flog,"(A)") "      <Non_local_modeling Modelname= ""user_nloc3"" ....> is defined in xml algorithm file "
     write(Flog,"(A)") "      BUT 'init_user_nloc_variables3' is not implemented "

     write(output_unit,"(A)") "ERROR :"
     write(output_unit,"(A)") "      <Non_local_modeling Modelname= ""user_nloc3"" ....> is defined in xml algorithm file "
     write(output_unit,"(A)") "      BUT 'init_user_nloc_variables3' is not implemented "
  end if

end subroutine init_user_nloc_variables3


end module non_local_user_mod
