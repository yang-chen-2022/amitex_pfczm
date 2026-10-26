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
  use decomp_2d, only : mytype, nrank
#endif  
  use decomp_2d_fft

! AMITEX MODULES (all available)
!-------------------------------
  use error_mod
  use amitex_mod       !, only : Nloc,GNloc  ! non-local fields, access to other fields
  use param_algo_mod   ! , only : algo_param
  use material_mod     !, only : Assign_VarInt_from_Field,Assemble_field_from_VarInt,&
                       !    mattotP,calcrotmatrices_GD, Nloc_models
  use field_mod
  use green_mod


  implicit none

  private

  !> public subroutines (one per non-local model)
  public :: user_nloc_after1,user_nloc_after2,user_nloc_after3,&
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
!                             SUBROUTINE GURTIN_SS
!
!> Implementation du modele de Gurtin formule en petites deformations
!
!     Contrainte interne ajoutee a la loi d'écoulement valant -lambda*rot(rot(Hp)):mu
!
!     NLOC : 9 composantes --> Hp = Fp - Id
!
!     GNLOC : 10 composantes 
!               * 1-9 --> -rot(rot(Hp)) = S
!               *  10 --> norm(rot(Hp)) = norm(Nye)
!
! Hp est toujour dans le repere local du crystal (convention de l'UMAT/MFRONT) au sein du 
! tableau VarInt. C'est simplement lors de son stockage temporaire dans
! Nloc%Var qu'il est tourne dans le repere global. C'est indispensable pour le calcul
! de ses derivees car les operateurs de derivation sont definis dans le repere global
! de même que la FFT. 
!
! De même S est toujours dans le repere local au sein du tableau Varint (convention UMAT/MFRONT)
! Il est obtenu par calcul de derivees dans cette fonction dans le repere global. Il faut 
! donc le tourner dans le repere local avant de le renvoyer dans le tableau
! Varint des materiaux.
!
! OPTIONS "utilisateur"
!----------------------
! nloc_param(modelInd)%P_real(1) : 3 options pour les derivation discrete 
!  si = 1. : rot = DF centree	                     rot(rot()) = operateur rotrot DF
!  si = 2. : rot = DF centree	                     rot(rot()) = rot DF centree applique 2 fois 
!  si = 3. : rot = DF de la meca (hexa par defaut)   rot(rot()) = rot DF de la meca applique 2 fois
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

 !-- check user defined parameters

 if ((nloc_param(modelInd)%P_real(1) .ne. 1) .and.  &
     (nloc_param(modelInd)%P_real(1) .ne. 2) .and.  &
     (nloc_param(modelInd)%P_real(1) .ne. 3)) then
   call amitex_abort("Parameter P_real(1) for gurtin_ss different from 1, 2 or 3",2,0)
 end if 

end subroutine init_user_nloc_variables1



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

 ! Autre variables
 integer,dimension(3)                         :: N          ! dimensions de la cellule 
 real(mytype),dimension(3)                    :: d          ! dimensions des voxels 
 integer                                      :: i,j        ! variables de boucles (sur materiau etc...)
 integer(kind=INT64)                          :: k,l,m      ! variables de boucles (sur positioin, zone etc..)
 integer(kind=INT64)                          :: indice_pos  
 real(mytype), dimension(3)                   :: V1,V2      ! Coordonnées des vecteurs directeurs du repere local du crystal
 real(mytype), dimension(9,9)                 :: RotMat,iRotMat   ! Matrices de changement de base entre repere local et repere du crystal
 real(mytype), dimension(9)                   :: Champ_tmp 
 integer, dimension(3)                        :: ind_coeff_v1, ind_coeff_v2

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
 

 ! Passage du champ Hp du repère local (sortie d'UMAT) au repère global
 !--------------------------------------------------------------------
  do i=1,size(list_mat)
     !! On parcout les materiaux concernes par le modele de Gurtin HPP
 
     ! On parcours les materiaux sur le pinceau
     do j=1,size(mattotP)
        ! changement de base uniquement si le materiau precise se trouve sur le pinceau
        if (mattotP(j)%numM == list_mat(i)) then
!~             write(OUTPUT_UNIT,*) "Rang : ",nrank," materiau : ",list_mat(i)," nombre voxels : ",size(mattotP(j)%pos)
!~             write(OUTPUT_UNIT,*) "Rang : ",nrank," materiau : ",mattotP(j)%numM," tableau zone : ",mattotP(j)%zone
!~             write(OUTPUT_UNIT,*) "Rang : ",nrank," materiau : ",mattotP(i)%numM," tableau zone : ",mattotP(i)%zone
            
            ! indice minimum de zone (local)
            m=1            
            ! boucle sur les zones locales du materiau : 
            ! => les differents grains du crystal
            do k=1,size(mattotP(j)%zone(:,1))
            
                ! recuperation des vecteurs directeurs du repere local
                !V1 = mattotP(j)%Coeff(10:12,k)
                !V2 = mattotP(j)%Coeff(13:15,k)
                ind_coeff_v1 = Nloc_models(modelInd)%Ind_CoeffNloc(i,1:3)
                ind_coeff_v2 = Nloc_models(modelInd)%Ind_CoeffNloc(i,4:6)
                V1 = mattotP(j)%Coeff(ind_coeff_v1,k)
                V2 = mattotP(j)%Coeff(ind_coeff_v2,k)

                
                ! Calcul des matrices de changement de base
                call CalcRotMatrices_GD(V1,V2,RotMat,iRotMat)
            
                ! boucle sur les voxels de la zone 
                do l=m,MattotP(j)%zone(k,1)
                    indice_pos = mattotP(j)%pos(l) ! indice lineaire du voxel dans le champ
                    
                    ! recuperation du tenseur Hp sur le voxel 
                    Champ_tmp = Nloc(modelInd)%Var(indice_pos,:)
!~                     write(OUTPUT_UNIT,*) "Rang : ",nrank,"Pos : ", indice_pos, "Champ avant rot : ",Nloc(modelInd)%Var(indice_pos,:)
                    ! rotation du champ 
                    Nloc(modelInd)%Var(indice_pos,:) = matmul(iRotMat,Champ_tmp)
!~                     write(OUTPUT_UNIT,*) "Pos : ", indice_pos, "Champ apres rot : ", Nloc(modelInd)%Var(indice_pos,:)
                end do 
                m=MattotP(j)%zone(k,1)+1
            end do
        end if     
     end do
  end do

 ! Passage dans l'espace de Fourier (repere global)
 !!----------------------------------------------------------
 call field_fft(Nloc(modelInd)%Var,Nloc(modelInd)%VarF,Nnloc,1)


 ! ----> a ce moment  Nloc%Var = Hp

 ! Derivations dans l'espace de Fourier : calcul des champs (repere global)
 !!----------------------------------------------------------

 !     norm(rot(Hp)) Hp = NLOC  --> GNLOC(10)  
 !!----------------------------------------------------------
   ! calcul du rot(Hp)
   if (nloc_param(modelInd)%P_real(1) .ne. 3) then                            ! cas 1 et 2 : rot difference finies centree
      call rotF(Nloc(modelInd)%VarF,GNloc(modelInd)%VarF(:,:,:,1:9),N,d)      !
   else 
      call rotF(Nloc(modelInd)%VarF,GNloc(modelInd)%VarF(:,:,:,1:9),N,d,FREQ) ! cas 3 : rot a partir des frequences du pb meca (hexa par defaut)
   end if                                                                              

   ! renormalisation
   GNloc(modelInd)%VarF(:,:,:,1:9) = GNloc(modelInd)%VarF(:,:,:,1:9) / real(grid%ntot,mytype) 
   ! retour dans l'espace reel
   call field_ifft(GNloc(modelInd)%Var(:,1:9),GNloc(modelInd)%VarF(:,:,:,1:9),9,1)  
   ! calcul de la norme 
   call field_normcomp(GNloc(modelInd)%Var(:,1:9),GNloc(modelInd)%Var(:,10),9,1)

 !      -rot(rot(Hp)) Hp = NLOC --> GNLOC(1:9)
 !!----------------------------------------------------------
   !! calcul du rotrot(Hp) (ecrase le rot)
   if     (nloc_param(modelInd)%P_real(1)==1) then
     call rotrotF(Nloc(modelInd)%VarF,GNloc(modelInd)%VarF(:,:,:,1:9),N,d)          ! cas 1 : rot_rot() differences finies 
   elseif (nloc_param(modelInd)%P_real(1)==2) then
     Nloc(modelInd)%VarF = GNloc(modelInd)%VarF(:,:,:,1:9) * real(grid%ntot,mytype) ! Nloc = rot(Hp) 
     call rotF(Nloc(modelInd)%VarF,GNloc(modelInd)%VarF(:,:,:,1:9),N,d)             ! cas 2 : rot(rot()) differences finies centrees
   elseif (nloc_param(modelInd)%P_real(1)==3) then
     call rotrotF(Nloc(modelInd)%VarF,GNloc(modelInd)%VarF(:,:,:,1:9),N,d,FREQ)     ! cas 3 : rot(rot()) avecfrequences du pb meca
   end if

   ! renormalisation
   GNloc(modelInd)%VarF(:,:,:,1:9) = GNloc(modelInd)%VarF(:,:,:,1:9) / real(grid%ntot,mytype) 
   ! retour dans l'espace reel
   call field_ifft(GNloc(modelInd)%Var(:,1:9),GNloc(modelInd)%VarF(:,:,:,1:9),9,1)  
   ! a ce moment ----> GNloc()%Var(:,1,9) = rot(rot(Hp))
   ! s = - A * rot(rot(Fp))   --> la multiplication par A est faite dans l'UMAT
   GNloc(modelInd)%Var(:,1:9) = -1*GNloc(modelInd)%Var(:,1:9)

 ! Passage du champ S du repère global au repère local
 !--------------------------------------------------------------------
  do i=1,size(list_mat)
     !! On parcout les materiaux concernes par le modele de Gurtin HPP
 
     ! On parcours les materiaux sur le pinceau
     do j=1,size(mattotP)
        ! changement de base uniquement si le materiau precise se trouve sur le pinceau
        if (mattotP(j)%numM == list_mat(i)) then
            
            ! indice minimum de zone (local)
            m=1            
            ! boucle sur les zones locales du materiau : 
            ! => les differents grains du crystal
            do k=1,size(mattotP(j)%zone(:,1))
            
                ! recuperation des vecteurs directeurs du repere local
                !V1 = mattotP(j)%Coeff(10:12,k)
                !V2 = mattotP(j)%Coeff(13:15,k)
                ind_coeff_v1 = Nloc_models(modelInd)%Ind_CoeffNloc(i,1:3)
                ind_coeff_v2 = Nloc_models(modelInd)%Ind_CoeffNloc(i,4:6)
                V1 = mattotP(j)%Coeff(ind_coeff_v1,k)
                V2 = mattotP(j)%Coeff(ind_coeff_v2,k)

                ! Calcul des matrices de changement de base
                call CalcRotMatrices_GD(V1,V2,RotMat,iRotMat)
            
                ! boucle sur les voxels de la zone 
                do l=m,MattotP(j)%zone(k,1)
                    indice_pos = mattotP(j)%pos(l) ! indice lineaire du voxel dans le champ
                    
                    ! recuperation du tenseur S
                    Champ_tmp = GNloc(modelInd)%Var(indice_pos,1:9)
                    
                    ! rotation du repere global au local
                    GNloc(modelInd)%Var(indice_pos,1:9) = matmul(RotMat,Champ_tmp)
                    
                end do 
                m=MattotP(j)%zone(k,1)+1
            end do
        end if     
     end do
  end do

 ! Repartition des variables internes non locales dans la structure mattotP
 !!----------------------------------------------------------
 call Assign_VarInt_from_Field(list_mat,list_varInt_GNloc,NGnloc,GNloc(modelInd)%Var)

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
