!===================================================================
!
!       MODULE ALGO_FUNCTIONS
!> Module gerant les calculs lies a l'algorithme iteratif
!!
!!
!!  Subroutines
!! - initDef : calcul du champ de deformation a la premiere iteration
!!             de l'algorithme principal
!! - initDefD 
!! - updateEimp : Mise a jour du champ de deformation a imposer
!! - acti3 : calcul de la solution obtenue en appliquant 
!!           l'acceleration de convergence
!! - SigToPK1 : Subroutine de calcul du champ de contrainte de Piola-Kirchhoff
!! - PK1ToSig : Subroutine de calcul du champ de contrainte de Cauchy
!!
!===================================================================
module algo_functions_mod

#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype, real_type
  use decomp_2d_mpi,       only : nrank
#else
  use decomp_2d, only : nrank,mytype,real_type 
#endif  
  use mpi
  use error_mod
  use linear_mod
  use param_algo_mod
  use amitex_mpi_mod, only: amitex_mpi_comm
#ifdef OPENMP
  use omp_lib
#endif
  private

  public:: initDef, initDefD, updateEimp,updateEimpD, &
           acti3, ACT3,ACT3_Moy,PK1ToSig,SigToPK1,&
           ACT3local, optimize_scalar

contains


!===================================================================
!
!> \brief     SUBROUTINE D'INITIALISATION DU CHAMP DE DEFORMATION
!!
!! Calcul de la deformation initiale pour la methode iterative
!!
!! Au premier increment du temps (load_n = 1, load_incr = 0),
!! on choisit \f$ \epsilon_0 = E^{imp}_n \f$
!! ou \f$ E^{imp}_{n-1} \f$ est la deformation imposee
!! dans le cas ou uniquement les deformations sont imposees.\n
!! Si certaines contraintes sont imposees on calcule 
!! la deformation Eimp correspondante en supposant que 
!! le materiau est lineaire elastique de tenseur C0.
!! On applique ensuite la meme formule.
!!
!! Pour les premiers increments d'un nouveau chargement (load_n > 1, load_incr = 0),
!! on choisit \f$ \epsilon_n = \epsilon_{n-1} \f$
!!
!! Dans les autres cas, on calcule
!! \f[ \epsilon_n = \epsilon_{n-1} + \frac{t_n - t_{n-1}}{t_{n-1} - t_{n-2}}(\epsilon_{n-1} - \epsilon_{n-2}).\f]
!!
!! \param[in] local_loading: chargement impose
!!                -  local_loading(1:algo_param%nTensDef) : pilotage
!!                -  local_loading(algo_param%nTensDef+1:2*algo_param%nTensDef) : valeur du pilotage
!! \param[in] nb_stress_pilot: nombre de composantes ou le pilotage est en contrainte
!! \param[in] t_load:  Temps de chargement aux instant t, t-1, t-2
!! \param[in] load_incr:  increment du chargement
!! \param[in] load_n:  numero du pas de chargement
!! \param[in] C0:  Tenseur decrivant le materiau de reference
!! \param[in] Def0:  déformation au pas de temps t-2
!! \param[in] Def:  déformation au pas de temps précédent
!!
!! \param[out] Def0:  déformation au pas de temps precedent
!! \param[out] Def:  proposition de déformation initiale
!-----------------------------------------------------------------
  subroutine initDefD(local_loadingD,nb_flux_pilot,t_load,load_incr,load_n,K0D,GradQD0,GradQD)

    implicit none

    real(mytype),dimension(6*algo_param%nVarD), intent(in)     :: local_loadingD
    Integer, intent(in)                        :: load_incr, load_n
    real(mytype),dimension(-2:0),intent(in)    :: t_load
    real(mytype),dimension(:,:,:), intent(inout) :: GradQD0, GradQD
    Integer                                    :: i,m,j,l,ivar
    integer,intent(in)                         :: nb_flux_pilot !< nombre de composantes ou le pilotage est en contrainte
    real(mytype),dimension(nb_flux_pilot)      :: EimpS
    real(mytype),dimension(3*algo_param%nVarD-nb_flux_pilot)    :: EimpE
    !! Contraintes et deformations moyennes
    real(mytype),dimension(algo_param%nVarD), intent(in) :: K0D
    !! Uniquement au premier chargement :
    !! On calcule le champs de deformation impose correspondant
    !! au pilotage impose
    if(load_n==1 .AND. load_incr==0) then
       GradQD0 = GradQD
       if(nb_flux_pilot>0) then
          ! CALCUL EimpE et EimpS
          j=1
          l=1
          do ivar=1,algo_param%nVarD
          do i = 1,3
             m = 3*(ivar-1) + i
             if(local_loadingD(m) == 1) then
                EimpS(j) = - local_loadingD(m+3*algo_param%nVarD) / K0D(ivar)
                j=j+1
             else
                EimpE(l) = local_loadingD(m+3*algo_param%nVarD)
                l=l+1
             end if
          end do
          end do
          ! Affectation de GradQD          
          j=1
          l=1
          do ivar=1,algo_param%nVarD
          do i = 1,3
             m = 3*(ivar-1) + i
             if(local_loadingD(m) == 1) then
                GradQD(:,i,ivar) = EimpS(j)
                j=j+1
             else
                GradQD(:,i,ivar) = EimpE(l)
                l=l+1
             end if
          end do
          end do
       else
          do ivar=1,algo_param%nVarD
          do i = 1,3
             GradQD(:,i,ivar) = local_loadingD(3*algo_param%nVarD + 3*(ivar-1) + i)
          end do
          end do
       end if      
    elseif(load_incr==0) then
       !! Cas du pas de temps 0 quand on change de chargement,
       !! On fait un demi pas de temps en prenant comme condition initiale
       !! la solution obtenue au pas de temps precedent
       GradQD0 = GradQD
    else
       ! Def0 = Def
       ! DDef = Def - Def0(iteration precedente)
       GradQD = GradQD + (t_load(0)-t_load(-1))/(t_load(-1)-t_load(-2))*(GradQD-GradQD0)
       ! Calcul arithmetique evitant un stockage temporaire
       ! Equivalent a Def0 = Def avant le calcul precedent
       GradQD0= (GradQD*(t_load(-1)-t_load(-2))+GradQD0*(t_load(0)-t_load(-1)))&
            /(t_load(0)-t_load(-2))
    end if
  end subroutine initDefD
!====================================================================

!====================================================================
!>  \brief SUBROUTINE DE MISE A JOUR DU CHAMP DE DEFORMATION IMPOSE
!!
!! Calcule le champ de deformation moyen qu'il faut imposer.\n
!! Doit prendre en compte le cas ou une contrainte est imposee
!!
!! \param[out] Eimp: Champ de deformation impose a l'instant t. Calcule
!! a partir des deformations et des contraintes imposees 
!! \param[in] local_loading: Chargement impose
!!               -  local_loading(1:3*algo_param%nVarD) : pilotage
!!               -  local_loading(3*algo_param%nVarD+1:6*algo_param%nVarD) : valeur du pilotage
!! \param[in] C0: tenseur decrivant le materiau de reference
!! \param[in] SigF: Contrainte a l'instant t-1 dans l'espace de Fourier
!! \param[in] DefF: Deformation a l'instant t-1 dans l'espace de Fourier
!! \param[in] ntot: Nombre de voxels dans la cellule
!! \param[in] fft_start: indices de la grille de fréquences (pour chaque pinceau)
!====================================================================
  subroutine updateEimpD(EimpD, local_loadingD,K0D, FluxDF, GradQDF, ntot,fft_start)


    implicit none

    integer(KIND=8),intent(in)                            :: ntot
    !> transformee de Fourier discrete des contraintes et des deformations
    complex(mytype),dimension(:,:,:,:,:),intent(in)       :: FluxDF,GradQDF 
    !> chargement en un temps donne
    real(mytype),dimension(6*algo_param%nVarD),intent(in) :: local_loadingD 
    real(mytype),dimension(3*algo_param%nVarD),intent(out):: EimpD !< deformation imposee
    !> Matrice de raideur du materiau de reference
    real(mytype),dimension(algo_param%nVarD), intent(in)  :: k0D 

    !> Valeur moyenne de la polarisation, Eimp temporaire
    real(mytype),dimension(3*algo_param%nVarD)            :: Tau, EimpD0 

    integer                           :: i,m,ivar !<indices
    integer, dimension(3), intent(in) :: fft_start
    integer                           :: ierr
 
   
    EimpD = 0._mytype    
    EimpD0 = 0._mytype    
    Tau=0._mytype

    if((fft_start(1) .eq. 1) .AND. (fft_start(2) .eq. 1) .AND. (fft_start(3) .eq. 1)) then
    !! On recupere les valeurs imposees
    EimpD0 = local_loadingD(3*algo_param%nVarD+1:6*algo_param%nVarD)
   
       do ivar=1,algo_param%nVarD
       !! On calcule \tau = \FluxD_{moy} + K0 : \GradQD_{moy}
       !! La renormalisation par ntot est deja faite pour DefF
       tau(3*(ivar-1)+1:3*ivar) = real(FluxDF(1,1,1,:,ivar)/real(ntot, mytype),mytype) + &
            K0D(ivar)*real(GradQDF(1,1,1,:,ivar),mytype)
       do i = 1,3
          m = 3*(ivar-1) + i
          !! On calcule et on affecte a EimpD0 les composantes imposees en Flux
          if(local_loadingD(m) == 1) then
             EimpD0(m) = (tau(m) - local_loadingD(m+3*algo_param%nVarD)) / K0D(ivar)
          end if
       end do
       end do

    end if

    !! Astuce pour distribuer Eimp sur tous les procs sans connaitre 
    !! le proc. contenant la fréquence (1,1,1) :
    !! On initialise a 0 au depart et on fait la somme a la fin du calcul
    call MPI_AllReduce(EimpD0,EimpD,3*algo_param%nVarD,real_type,MPI_SUM,amitex_mpi_comm,ierr)

  
  end subroutine updateEimpD
!====================================================================








































!===================================================================
!
!> \brief     SUBROUTINE D'INITIALISATION DU CHAMP DE DEFORMATION
!!
!! Calcul de la deformation initiale pour la methode iterative
!!
!! Au premier increment du temps (load_n = 1, load_incr = 0),
!! on choisit \f$ \epsilon_0 = E^{imp}_n \f$
!! ou \f$ E^{imp}_{n-1} \f$ est la deformation imposee
!! dans le cas ou uniquement les deformations sont imposees.\n
!! Si certaines contraintes sont imposees on calcule 
!! la deformation Eimp correspondante en supposant que 
!! le materiau est lineaire elastique de tenseur C0.
!! On applique ensuite la meme formule.
!!
!! Pour les premiers increments d'un nouveau chargement (load_n > 1, load_incr = 0),
!! on choisit \f$ \epsilon_n = \epsilon_{n-1} \f$
!!
!! Dans les autres cas, on calcule
!! \f[ \epsilon_n = \epsilon_{n-1} + \frac{t_n - t_{n-1}}{t_{n-1} - t_{n-2}}(\epsilon_{n-1} - \epsilon_{n-2}).\f]
!!
!! \param[in] local_loading: chargement impose
!!                -  local_loading(1:algo_param%nTensDef) : pilotage
!!                -  local_loading(algo_param%nTensDef+1:2*algo_param%nTensDef) : valeur du pilotage
!! \param[in] nb_stress_pilot: nombre de composantes ou le pilotage est en contrainte
!! \param[in] t_load:  Temps de chargement aux instant t, t-1, t-2
!! \param[in] load_incr:  increment du chargement
!! \param[in] load_n:  numero du pas de chargement
!! \param[in] C0:  Tenseur decrivant le materiau de reference
!! \param[in] Def0:  déformation au pas de temps t-2
!! \param[in] Def:  déformation au pas de temps précédent
!! \param[in] Def_nsym : gradient du deplacement au pas de temps precedent (optionnel)
!! \param[in] Def_nsym0 : gradient du deplacement au pas de temps t-2 (optionnel)
!!
!! \param[out] Def0:  déformation au pas de temps precedent
!! \param[out] Def:  proposition de déformation initiale
!! \param[out] Def_nsym0 : gradient du deplacement au pas de temps precedent
!! \param[out] Def_nsym : proposition de gradient du deplacement initial
!-----------------------------------------------------------------
  subroutine initDef(local_loading,nb_stress_pilot,t_load,load_incr,load_n,C0,Def0,Def,Def_nsym0,Def_nsym)

    implicit none

    real(mytype),dimension(2*algo_param%nTensDef), intent(in)     :: local_loading
    Integer, intent(in)                        :: load_incr, load_n
    real(mytype),dimension(-2:0),intent(in)    :: t_load
    real(mytype),dimension(:,:), intent(inout) :: Def0, Def
    Integer                                    :: i,j,k,l,m
    integer,intent(in)                    :: nb_stress_pilot !< nombre de composantes ou le pilotage est en contrainte
    real(mytype),dimension(nb_stress_pilot)      :: Simp, EimpS
    real(mytype),dimension(algo_param%nTensDef-nb_stress_pilot)    :: EimpE
    !! Contraintes et deformations moyennes
    real(mytype),dimension(nb_stress_pilot,algo_param%nTensDef-nb_stress_pilot) :: CSE
    real(mytype),dimension(nb_stress_pilot,nb_stress_pilot) :: CSS !< variables intermediaires
    real(mytype),dimension(algo_param%nTensDef,algo_param%nTensDef), intent(in) :: C0
    real(mytype),dimension(:,:), intent(inout),optional       :: Def_nsym0,Def_nsym


    ! CAS PARTICULIER 1er increment du 1er chargement
    !! On calcule le champs de deformation impose correspondant au pilotage impose
    if(load_n==1 .AND. load_incr==0) then
       Def0 = Def
       if(nb_stress_pilot>0) then
          j=1
          k=1
          l=1
          do i = 1,algo_param%nTensDef
             if(local_loading(i) == 1) then
                Simp(j) = local_loading(i+algo_param%nTensDef)
                k=1
                do m = 1,algo_param%nTensDef
                   if(local_loading(m) == 1) then
                      CSS(k,j) = C0(m,i)
                      k = k+1
                   end if
                end do
                j=j+1
             else
                EimpE(l) = local_loading(i+algo_param%nTensDef)
                k=1
                do m = 1,algo_param%nTensDef
                   if(local_loading(m) == 1) then
                      CSE(k,l) = C0(m,i)
                      k = k+1
                   end if
                end do
                l=l+1
             end if
          end do
          EimpS = Simp - MATMUL(CSE,EimpE)
          call LUSolve(CSS, EimpS)
          j=1
          k=1
          do i = 1,algo_param%nTensDef
             if(local_loading(i) == 1) then
                Def(:,i) = EimpS(j)
                j=j+1
             else
                Def(:,i) = EimpE(k)
                k=k+1
             end if
          end do
       else
          do i = 1,algo_param%nTensDef
             Def(:,i) = local_loading(i+algo_param%nTensDef)
          end do
       end if      
            if (algo_param%HPP_nsym) then
                ! Dans le cas de la prise en compte complete du gradient du 
                ! deplacement au sein d'un calcul HPP : le nouveau gradient
                ! de deplacement est impose symetrique et egal a Def. 
                Def_nsym0 = Def_nsym
                Def_nsym(:,1) = Def(:,1) 
                Def_nsym(:,2) = Def(:,2) 
                Def_nsym(:,3) = Def(:,3) 
                Def_nsym(:,4) = Def(:,4)/2 
                Def_nsym(:,5) = Def(:,5)/2 
                Def_nsym(:,6) = Def(:,6)/2 
                Def_nsym(:,7) = Def(:,4)/2 
                Def_nsym(:,8) = Def(:,5)/2 
                Def_nsym(:,9) = Def(:,6)/2 
            end if

    ! CAS PARTICULIER 1er increment de chargement
    elseif(load_incr==0) then
       !! Cas du pas de temps 0 quand on change de chargement,
       !! On fait un demi pas de temps en prenant comme condition initiale
       !! la solution obtenue au pas de temps precedent
       Def0 = Def
       if (algo_param%HPP_nsym) then
            ! Dans le cas de la prise en compte complete du gradient du 
            ! deplacement au sein d'un calcul HPP : 
            Def_nsym0 = Def_nsym
       end if

    ! CAS GENERAL 'default'
    elseif(algo_param%init_def=="default") then
       ! Def0 = Def
       ! DDef =Def - Def0(iteration precedente)
       Def = Def + (t_load(0)-t_load(-1))/(t_load(-1)-t_load(-2))*(Def-Def0)
       ! Calcul arithmetique evitant un stockage temporaire
       ! Equivalent a Def0 = Def avant le calcul precedent
       
       Def0= (Def*(t_load(-1)-t_load(-2))+Def0*(t_load(0)-t_load(-1)))&
            /(t_load(0)-t_load(-2))
        if (algo_param%HPP_nsym) then
            ! Dans le cas de la prise en compte complete du gradient du 
            ! deplacement au sein d'un calcul HPP : meme type d'initialisation
            Def_nsym = Def_nsym + (t_load(0)-t_load(-1))/(t_load(-1)-t_load(-2))&
                       *(Def_nsym-Def_nsym0)
            Def_nsym0= (Def_nsym*(t_load(-1)-t_load(-2))+Def_nsym0*(t_load(0)-t_load(-1)))&
                       /(t_load(0)-t_load(-2))
        end if

    ! CAS GENERAL 'previous'
    elseif(algo_param%init_def=="previous") then
        Def0=Def

    ! CASE 'None'    
    elseif(algo_param%init_def=="none") then ! case if algo_param%init_def='None' or something else
        !nothing to do  : Def0 = Def0 and Def = Def !
    else
        call amitex_abort("initDef : algo_param%init_def different from 'default', 'previous','none'",2,0)
    end if


  end subroutine initDef
!====================================================================



!====================================================================
!>  \brief SUBROUTINE DE MISE A JOUR DE LA DEFORMATION MOYENNE IMPOSE
!!
!! Calcule la deformation moyenne qu'il faut imposer.\n
!! Doit prendre en compte le cas ou une contrainte est imposee
!!
!! \param[out] Eimp: Deformation imposee a l'instant t. Calcule
!! a partir des deformations et des contraintes imposees 
!! \param[in] local_loading: Chargement impose
!!               -  local_loading(1:algo_param%nTensDef) : pilotage
!!               -  local_loading(algo_param%nTensDef+1:2*algo_param%nTensDef) : valeur du pilotage
!! \param[in] nb_stress_pilot nombre de composantes ou le pilotage est en contrainte
!! \param[in] C0: tenseur decrivant le materiau de reference
!! \param[in] SigF: Contrainte a l'instant t-1 dans l'espace de Fourier
!! \param[in] DefF: Deformation a l'instant t-1 dans l'espace de Fourier
!! \param[in] ntot: Nombre de voxels dans la cellule
!! \param[in] fft_start: dimensions de la grille de fréquences
!! \param[in] DirStress_flag : .true. si direction de contrainte multiaxiale impose
!!                et pilotage deformation
!! \param[in] DirStress_cauchy : .true. si contrainte de Cauchy, sinon PK1
!! \param[in] DirStress : (vecteur Ntens) direction de contrainte imposee
!! \param[in] DirStress2 : (vecteur Ntens) direction de contrainte imposee PK1
!!               reevalue si GDEF et DirStress_cauchy
!====================================================================
  subroutine updateEimp(Eimp, local_loading, nb_stress_pilot, &
            DirStress_flag, DirStress_cauchy,DirStress,DirStress2,C0, SigF, DefF, ntot,fft_start)

    implicit none
    !> nombre de voxels dans la cellule
    integer(KIND=8),intent(in)                               :: ntot
    !> transformee de Fourier discrete des contraintes et des deformations
    complex(mytype),dimension(:,:,:,:),intent(in)            :: SigF,DefF 
    !> chargement en un temps donne
    real(mytype),dimension(2*algo_param%nTensDef),intent(in) :: local_loading 
    !> nombre de composantes ou le pilotage est en contrainte
    integer,intent(in)                    :: nb_stress_pilot 
    !> test si chargement en direction de contrainte imposee (et pilotage deformation)
    !>      + test du type de contrainte en Grandes def. : cauchy ou pk1
    logical, intent(in)                                      :: DirStress_flag
    logical, intent(in)                                      :: DirStress_cauchy
    !> Direction de contrainte imposee (et pilotage deformation)
    !>      + direction pk1 si contrainte de Cauchy imposee 
    real(mytype),dimension(algo_param%nTensDef), intent(in)  :: DirStress
    real(mytype),dimension(algo_param%nTensDef), intent(out) :: DirStress2
    real(mytype),dimension(algo_param%nTensDef),intent(out)  :: Eimp !< deformation imposee
    !> Matrice de raideur du materiau de reference
    real(mytype),dimension(algo_param%nTensDef,algo_param%nTensDef), intent(in) :: C0 
    !> Valeur moyenne de la polarisation, Eimp temporaire, Direction de contrainte imposee
    real(mytype),dimension(algo_param%nTensDef)              :: Tau, Eimp0, Simp_DirStress,&
                                                                Smoy, Emoy,DirStress20
    !< Restriction de C0 aux composantes en contrainte imposees
    real(mytype),dimension(nb_stress_pilot,algo_param%nTensDef-nb_stress_pilot) :: CSE
    real(mytype),dimension(nb_stress_pilot,nb_stress_pilot)  ::CSS 
    real(mytype),dimension(nb_stress_pilot)                  :: A_C0, TauS, Simp, EimpS
    real(mytype),dimension(algo_param%nTensDef-nb_stress_pilot)   :: EimpE 
    real(mytype) :: num,denom
    !< variables intermediaires : 
    !!  - TauS est la polarisation moyenne sur les composantes ou la contrainte est imposee
    !!  - Simp(i) sont les contraintes imposees
    !!  - EimpS(i) sont les deformations a imposer correspondant aux Simp
    !!  - EimpE(i) sont les deformations imposees
    !!  - A_C0 = CSE * EimpE
    integer               :: i,j,k,l,m !<indices
    integer, dimension(3), intent(in) :: fft_start
    integer               :: ierr
 
    Eimp = 0._mytype    
    Eimp0 = 0._mytype    
    DirStress2 = 0._mytype
    DirStress20 = 0._mytype


    !! ASTUCE :  on effectue le calcul sur le pinceau contenant la Freq nulle
    !!           PUIS on redistribue (MPI_ALLREDUCE, SUM)
    if((fft_start(1) .eq. 1) .AND. (fft_start(2) .eq. 1) .AND. (fft_start(3) .eq. 1)) then
    
    !! La renormalisation par ntot est deja faite pour DefF
    Smoy = real(SigF(1,1,1,:)/real(ntot, mytype),mytype)
    Emoy = real(DefF(1,1,1,:),mytype)

    ! Chargement en direction de contrainte imposee (et pilotage deformation)
    ! => evaluation de Simp_dirstress 
    if(DirStress_flag) then
       DirStress20 = DirStress
       !! evaluation de DirStress PK1 si GDEF et contrainte de Cauchy imposee
       if((.not. algo_param%HPP) .AND. DirStress_cauchy) then
           call DirStress_SigToPK1(DirStress,DirStress20,Emoy)
       end if       
       !! alpha = (Smoy.DirStress)/(DirStress.DirStress)
       num = dot_product(Smoy,DirStress20) 
       if (algo_param%HPP) num = num + dot_product(Smoy(4:6),DirStress20(4:6))
       denom = dot_product(DirStress20,DirStress20) 
       if (algo_param%HPP) denom = denom + dot_product(DirStress20(4:6),DirStress20(4:6))
       Simp_dirstress = DirStress20 * num / denom
    end if
     
    Eimp0 = local_loading(algo_param%nTensDef+1:2*algo_param%nTensDef)
    !! On calcule \tau = \Sigma_{moy} - C_0 : E_{moy}
    tau = Smoy - MATMUL(C0,Emoy)

    j=1
    k=1
    l=1
    do i = 1,algo_param%nTensDef
       if(local_loading(i) == 1) then
          k=1
          do m = 1,algo_param%nTensDef
             if(local_loading(m) == 1) then
                CSS(k,j) = C0(m,i)
                k = k+1
             end if
          end do
          TauS(j) = Tau(i)
          if (DirStress_flag) then
              Simp(j) = Simp_DirStress(i)
          else 
              Simp(j) = local_loading(i+algo_param%nTensDef);
          end if
          j=j+1
       else
          EimpE(l) = local_loading(i+algo_param%nTensDef);

          k=1
          do m = 1,algo_param%nTensDef
             if(local_loading(m) == 1) then
                CSE(k,l) = C0(m,i)
                k = k+1
             end if
          end do
          l=l+1
       end if
    end do

    A_C0 = MATMUL(CSE,EimpE)

    EimpS = Simp - A_C0 - TauS
    call LUSolve(CSS,EimpS)

    j=1
    do i = 1,algo_param%nTensDef
       if(local_loading(i) == 1) then
          Eimp0(i) = EimpS(j)
          j=j+1
       end if
    end do

    end if !! end if test sur fft_start

    !! Astuce pour distribuer Eimp sur tous les procs sans connaitre 
    !! le proc. contenant la fréquence (1,1,1) :
    !! On initialise a 0 au depart et on fait la somme a la fin du calcul
    call MPI_AllReduce(Eimp0,Eimp,algo_param%nTensDef,real_type,MPI_SUM,amitex_mpi_comm,ierr)
    call MPI_AllReduce(DirStress20,DirStress2,algo_param%nTensDef,real_type,MPI_SUM,amitex_mpi_comm,ierr)

  end subroutine updateEimp
!====================================================================


!===================================================================
!
!        FONCTION ACTI3
!
!> \brief Propose une solution au solveur iteratif en utilisant 
!! les quatre calculs precedents.\n
!!
!! Extrait notice CAST3M:\n
!! A partir de 4 champs de deplacements et des 4 champs d'increments
!! de forces correctrices associes, ACTI3 evalue le champs d'increments
!! de deplacement permettant de resoudre le probleme sur la restriction
!! de l'application tangente au sous-espace de dimension 3 defini par
!! les champs en entree
!!
!! \param[in]     ACT3_U: champs de déformation
!!                 tableaux de réels (amx_size(1)*amx_size(2)*amx_size(3),algo_param%nTensDef,4)
!! \param[in]     ACT3_R: residus correspondants
!! \param[in,out] Def: champ de deformation propose 
!!
!-----------------------------------------------------------------
  subroutine acti3(ACT3_U,ACT3_R,Def)

    implicit none
   
    !> CONTRAINTES ET DEFORMATIONS
    !! dans l'espace reel, coherent avec la structure "MATERIAL"
    real(mytype),dimension(:,:,:), intent(in)  :: ACT3_R, ACT3_U
    real(mytype),dimension(:,:),intent(inout)  :: Def
    real(mytype)                               :: A2,B2,J2,D2,E2,K2,A3,J3,L,J4,M,J5,N,maxR,tmp
    real(mytype),dimension(9)                  :: A,B
    integer :: ierr
    !! Tolerance epsilon differente en simple et double precision
#ifdef DOUBLE_PREC
    real(mytype),parameter                     :: epsilon = 1e-30
#else
    real(mytype),parameter                     :: epsilon = 1e-20
#endif

    L = 0._mytype
    M = 0._mytype
    N = 0._mytype
    !! On calcule maxR : la valeur maximale sur tous les residus stockes
    tmp = maxval(abs(ACT3_R(:,:,:)))
    call MPI_allreduce( tmp,maxR,1,real_type,MPI_MAX,amitex_mpi_comm,ierr)
    
    !! On renormalise toutes les sommes en divisant par maxR pour eviter d'atteindre 
    !! la valeur maximale d'un reel 
    A(1) = SUM(((ACT3_R(:,:,2) -ACT3_R(:,:,1))/maxR)**2)
    A(2) = SUM(((ACT3_R(:,:,3) -ACT3_R(:,:,1))/maxR)**2)
    A(3) = SUM(((ACT3_R(:,:,4) -ACT3_R(:,:,1))/maxR)**2)
    A(4) = SUM((ACT3_R(:,:,2)-ACT3_R(:,:,1))/maxR*(ACT3_R(:,:,3) -ACT3_R(:,:,1))/maxR)
    A(5) = SUM((ACT3_R(:,:,2)-ACT3_R(:,:,1))/maxR*(ACT3_R(:,:,4) -ACT3_R(:,:,1))/maxR)
    A(6) = SUM((ACT3_R(:,:,3)-ACT3_R(:,:,1))/maxR*(ACT3_R(:,:,4) -ACT3_R(:,:,1))/maxR)
    A(7) = -SUM(ACT3_R(:,:,1)/maxR*(ACT3_R(:,:,2)-ACT3_R(:,:,1))/maxR)
    A(8) = -SUM(ACT3_R(:,:,1)/maxR*(ACT3_R(:,:,3)-ACT3_R(:,:,1))/maxR)
    A(9) = -SUM(ACT3_R(:,:,1)/maxR*(ACT3_R(:,:,4)-ACT3_R(:,:,1))/maxR)
    call MPI_allreduce( A,B,9,real_type,MPI_SUM,amitex_mpi_comm,ierr)
    if (ABS(B(3)) > epsilon) then
       A2=B(1)-(B(5)**2/B(3)); B2=B(4)-(B(6)*B(5)/B(3));   J2=B(7)-(B(9)*B(5)/B(3));
       D2=B(4)-(B(5)*B(6)/B(3)); E2=B(2)-(B(6)**2/B(3));   K2=B(8)-(B(9)*B(6)/B(3));
       if (abs(E2) > epsilon) then
          A3=A2-(D2*B2/E2)
          J3=J2-(K2*B2/E2)
          if (abs(A3) > epsilon) then
             L=J3/A3
             J4=J2-(A2*L)             
             if (abs(B2) > epsilon) then
                M=J4/B2
                J5=B(7)-(B(1)*L)-(B(4)*M)
                if (abs(B(5)) > epsilon) then
                   N=J5/B(5);
                end if
             end if
          end if
       end if
    end if
    Def = L * (ACT3_U(:,:,2) - ACT3_U(:,:,1)) + M *(ACT3_U(:,:,3) - ACT3_U(:,:,1)) + N  * (ACT3_U(:,:,4) - ACT3_U(:,:,1)) 
    Def = Def + ACT3_U(:,:,1);
    
  end subroutine acti3
!====================================================================
!===================================================================
!
!        FONCTION ACT3
!
!> \brief Propose une solution au solveur iteratif en utilisant 
!! les quatre calculs precedents.\n
!!
!! Extrait notice CAST3M:\n
!! REND UN CHAMP PRIMAL EXTRAPOLE QUI MINIMISE LE CHAMP DUAL.
!! METHODE : ON ESTIME L'APPLICATION TANGENTE ET ON L'UTILISE
!! POUR RESOUDRE LE PROBLEME
!!
!! \param[in] ACT3_U: champs de déformation
!!                 tableaux de réels (amx_size(1)*amx_size(2)*amx_size(3),algo_param%nTensDef,4)
!! \param[in] ACT3_R: residus correspondants
!! \param[out] Def: champ de deformation propose 
!! \param[out] GradQ : champ de Gradient propose
!! \param[out] acv_test : .false. si pas d'acceleration
!!
!-----------------------------------------------------------------
  subroutine ACT3(ACT3_U,ACT3_R,Def,GradQD,acv_test)

    implicit none
   
    !> CONTRAINTES ET DEFORMATIONS
    !! dans l'espace reel, coherent avec la structure "MATERIAL"
    real(mytype),dimension(:,:,:), intent(in)  :: ACT3_R, ACT3_U
    real(mytype),dimension(:,:),intent(out)    :: Def
    real(mytype),dimension(:,:,:),intent(out)  :: GradQD
    logical                                    :: acv_test
    real(mytype)                               :: L,M,N,maxR,tmp
    real(mytype)                               :: F12,F22,F32,F12p,F22p,F32p,F1F2,F1F3,&
         F2F3,F0F1,F0F2,F0F3,F1F2p,F1F3p,F2F3p,F0F1p,F0F2p,F0F3p,&
         A1,B1,C1,J1,D1,E1,F1,K1,G1,H1,I1,L1,A2,B2,J2,D2,E2,K2,A3,B3,J3,B4,J4,C5,J5,&
         ref
    integer :: ierr
    !! Tolerance epsilon differente en simple et double precision
#ifdef DOUBLE_PREC
    real(mytype),parameter                     :: XPETIT = 1e-30_mytype
#else
    real(mytype),parameter                     :: XPETIT = 1e-20_mytype
#endif
    integer                                    :: i,j,k
    
    !! INITIALISATIONS (potentiellement non-initialises d'apres gcc -Wall)
    K1 = 0._mytype
    J1 = 0._mytype
    E1 = 0._mytype
    D1 = 0._mytype
    B1 = 0._mytype
    A1 = 0._mytype


    acv_test = .true.

    !!CAS  A 4 POINTS :
    !!-----------------
    !!famille des 3 vecteurs variations du champ dual
    L = 0._mytype
    M = 0._mytype
    N = 0._mytype
    !! On calcule maxR : la valeur maximale sur tous les residus stockes
    tmp = maxval(abs(ACT3_R(:,:,:)))
    call MPI_allreduce( tmp,maxR,1,real_type,MPI_MAX,amitex_mpi_comm,ierr)
    !! On renormalise toutes les sommes en divisant par maxR pour eviter d'atteindre 
    !! la valeur maximale d'un reel 

    F12p = SUM(((ACT3_R(:,:,2)-ACT3_R(:,:,1))/maxR)**2)
    F22p = SUM(((ACT3_R(:,:,3)-ACT3_R(:,:,1))/maxR)**2)
    F32p = SUM(((ACT3_R(:,:,4)-ACT3_R(:,:,1))/maxR)**2)
    call MPI_allreduce(F12p,F12,1,real_type,MPI_SUM,amitex_mpi_comm,ierr)
    call MPI_allreduce(F22p,F22,1,real_type,MPI_SUM,amitex_mpi_comm,ierr)
    call MPI_allreduce(F32p,F32,1,real_type,MPI_SUM,amitex_mpi_comm,ierr)
    if(f12.LT.XPETIT*100._mytype) f12=XPETIT*100._mytype
    if(f22.LT.XPETIT*100._mytype) f22=XPETIT*100._mytype
    if(f32.LT.XPETIT*100._mytype) f32=XPETIT*100._mytype

    !!Détection de flip-flop ?
    !!(1 et 3 sont trop proches, malgré 2 différent de 1) 
    if (f12/f22.gt.1.2_mytype .and. f12/f32.gt.0.707_mytype .and.&
         f12/f32.lt.1.414_mytype) then
       call amitex_abort("act3 : plan B",-1,0)
       ! CALL ERREUR(858)
       L = (sqrt(f12)+sqrt(f22))/(2*sqrt(f12))
       M = 0._mytype
       N = 0._mytype
       GOTO 1010
    endif
    !!Construction d'une famille libre orthogonale
    F1F2p = SUM((ACT3_R(:,:,2)-ACT3_R(:,:,1))/maxR*(ACT3_R(:,:,3) -ACT3_R(:,:,1))/maxR)
    F1F3p = SUM((ACT3_R(:,:,2)-ACT3_R(:,:,1))/maxR*(ACT3_R(:,:,4) -ACT3_R(:,:,1))/maxR)
    F2F3p = SUM((ACT3_R(:,:,3)-ACT3_R(:,:,1))/maxR*(ACT3_R(:,:,4) -ACT3_R(:,:,1))/maxR)
    F0F1p = SUM(ACT3_R(:,:,1)/maxR*(ACT3_R(:,:,2)-ACT3_R(:,:,1))/maxR)
    F0F2p = SUM(ACT3_R(:,:,1)/maxR*(ACT3_R(:,:,3)-ACT3_R(:,:,1))/maxR)
    F0F3p = SUM(ACT3_R(:,:,1)/maxR*(ACT3_R(:,:,4)-ACT3_R(:,:,1))/maxR)
    call MPI_allreduce( F1F2p,F1F2,1,real_type,MPI_SUM,amitex_mpi_comm,ierr)
    call MPI_allreduce( F1F3p,F1F3,1,real_type,MPI_SUM,amitex_mpi_comm,ierr)
    call MPI_allreduce( F2F3p,F2F3,1,real_type,MPI_SUM,amitex_mpi_comm,ierr)
    call MPI_allreduce( F0F1p,F0F1,1,real_type,MPI_SUM,amitex_mpi_comm,ierr)
    call MPI_allreduce( F0F2p,F0F2,1,real_type,MPI_SUM,amitex_mpi_comm,ierr)
    call MPI_allreduce( F0F3p,F0F3,1,real_type,MPI_SUM,amitex_mpi_comm,ierr)


    A1  =  F12
    B1  =  F1F2
    C1  =  F1F3
    J1 = -F0F1
    D1  =  F1F2
    E1  =  F22
    F1  =  F2F3
    K1 = -F0F2
    G1  =  F1F3
    H1  =  F2F3
    I1 =  F32
    L1 = -F0F3
    !!write (6,*) ' f12 f22 f32 ',f12,f22,f32
    
    !!Calcul des coordonnées de IFO0 dans la famille libre
    !!et élimination des vecteurs liés
    IF (ABS(I1) .LT. xpetit) GOTO 2000

    A2  = A1  - ( G1*C1/I1)
    B2  = B1  - ( H1*C1/I1)
    J2 = J1 - (L1*C1/I1)
    D2  = D1  - ( G1*F1/I1)
    E2  = E1  - ( H1*F1/I1)
    K2 = K1 - (L1*F1/I1)

    IF (ABS(E2) .LE. ABS(E1)*1E-6_mytype) THEN
       !!      write (6,*) ' act3 e2 e1 '
       GOTO 2000
    ENDIF

    A3  = A2 -(D2*B2/E2)
    J3 = J2-(K2*B2/E2)

    IF (ABS(A3) .LE. ABS(A1)*1E-6_mytype) THEN
       !!      write (6,*) ' act3 a3 a2 '
       GOTO 2000
    ENDIF
    IF (ABS(A2) .LE. ABS(A1)*1E-6_mytype) THEN
       !!      write (6,*) ' act3 a2 a1 '
       GOTO 2000
    ENDIF

    L=J3/A3

    B4  = B2
    J4 = J2 - (A2*L)

    IF (ABS(B4) .LE. ABS(B1)*1E-6_mytype) THEN
       !!      write (6,*) ' act3 b4 b1 '
       GOTO 2000
    ENDIF

    M  = J4/B4

    C5  = C1
    J5 = J1 - (A1*L) - (B1*M)

    IF (ABS(C5) .LT. xpetit) GOTO 2000

    N=J5/C5

    ref=sqrt(a1**2+b1**2+c1**2)
    !!write (6,*) ' act3 ',ref,XI1,A3,b4,c5

    GOTO 1010

2000 CONTINUE
    !!le systeme est singulier on essaye avec un vecteur de moins

    !!write (6,*) ' act3 ',xi1,e2,e1,a3,a1,b4,b1,c5
    IF (ABS(E1) .LT. xpetit) GOTO 2100

    A2  = A1  - (D1 *B1/E1)
    J2 = J1 - (K1*B1/E1)

    IF (ABS(A2) .LE. ABS(A1)*1E-6_mytype) THEN
       !!   write (6,*) ' act3 a2 a1 '
       GOTO 2100
    ENDIF

    L  = J2/A2

    B3  = B1
    J3 = J1 - (A1*L)

    IF (ABS(B3) .LT. xpetit) GOTO 2100
    !!
    M = J3/B3

    N = 0._mytype

    !!act3 : reduction a 2 dimensions
    call amitex_abort("act3 : reduction a 2 dimensions",-1,0)

    GOTO 1010

2100 CONTINUE
    !!le systeme est singulier on essaye avec un vecteur de moins
    IF (ABS(A1) .LT. xpetit ) GOTO 1000

    L  = J1/A1

    M  = 0._mytype
    N  = 0._mytype

    !!act3 : reduction a 1 dimensions
    call amitex_abort("act3 : reduction a 1 dimension",-1,0)
    GOTO 1010

1000 CONTINUE
    !!act3 : pas d accélération
    call amitex_abort("act3 : pas d'acceleration",-1,0)
    L = 0._mytype
    M = 0._mytype
    N = 0._mytype
    acv_test = .false.

1010 CONTINUE


    !!Vérification qu'on ne retombe pas sur l'un des duaux entrés
    !!Sinon, on réessaye en tenant compte de moins de vecteurs
    if (abs(l-1._mytype).lt.1e-1_mytype .and. abs(m).lt.1e-1_mytype&
         .and. abs(n).lt.1e-1_mytype) goto 1000
    if (abs(m-1._mytype).lt.1e-1_mytype .and. abs(n).lt.1e-1_mytype&
         .and. abs(l).lt.1e-1_mytype) goto 2100
    if (abs(n-1._mytype).lt.1e-1_mytype .and. abs(l).lt.1e-1_mytype&
         .and. abs(m).lt.1e-1_mytype) goto 2000

    if (algo_param%Mechanics) then
       Def = L * (ACT3_U(:,:,2) - ACT3_U(:,:,1)) + M *(ACT3_U(:,:,3) - ACT3_U(:,:,1)) + N  * (ACT3_U(:,:,4) - ACT3_U(:,:,1)) 
       Def = Def + ACT3_U(:,:,1)
    end if
    if (algo_param%Diffusion) then  
       do j=1,algo_param%nVarD
       do i = 1,3
          k = i+3*(j-1)
          GradQD(:,i,j) = L * (ACT3_U(:,k,2) - ACT3_U(:,k,1)) + M *(ACT3_U(:,k,3) - ACT3_U(:,k,1)) &
                        + N  * (ACT3_U(:,k,4) - ACT3_U(:,k,1)) 
          GradQD(:,i,j) = GradQD(:,i,j) + ACT3_U(:,k,1)
       end do
       end do 
    end if 
    RETURN
    
  end subroutine ACT3
!====================================================================

!====================================================================
!===================================================================
!
!        FONCTION ACT3LOCAL
!
!> \brief Propose une solution au solveur iteratif en utilisant 
!! les quatre calculs precedents.\n
!!  Version non parallèle de ACT3
!!
!! Extrait notice CAST3M:\n
!! REND UN CHAMP PRIMAL EXTRAPOLE QUI MINIMISE LE CHAMP DUAL.
!! METHODE : ON ESTIME L'APPLICATION TANGENTE ET ON L'UTILISE
!! POUR RESOUDRE LE PROBLEME
!!
!! \param[in] ACT3_U: champs de déformation
!!                 tableaux de réels (amx_size(1)*amx_size(2)*amx_size(3),algo_param%nTensDef*ntot,4)
!! \param[in] ACT3_R: residus correspondants
!! \param[in,out] Def: champ de deformation propose 
!!
!-----------------------------------------------------------------
 subroutine ACT3LOCAL(ACT3_U,ACT3_R,Def)

    implicit none

    !> CONTRAINTES ET DEFORMATIONS
    !! dans l'espace reel, coherent avec la structure "MATERIAL"
    real(mytype),dimension(:,:), intent(in)       :: ACT3_R, ACT3_U
    real(mytype),dimension(:),intent(inout)       :: Def
    real(mytype)                                  :: L,M,N,maxR
    real(mytype)                                  :: F12,F22,F32,F1F2,F1F3,&
         F2F3,F0F1,F0F2,F0F3,&
         A1,B1,C1,J1,D1,E1,F1,K1,G1,H1,I1,L1,A2,B2,J2,D2,E2,K2,A3,B3,J3,B4,J4,C5,J5,&
         ref
    !! Tolerance epsilon differente en simple et double precision

    real(mytype),parameter                     :: XPETIT = 1e-30

    !! INITIALISATIONS (potentiellement non-initialises d'apres gcc -Wall)
    K1 = 0._mytype
    J1 = 0._mytype
    E1 = 0._mytype
    D1 = 0._mytype
    B1 = 0._mytype
    A1 = 0._mytype


    !!CAS  A 4 POINTS :
    !!-----------------
    !!famille des 3 vecteurs variations du champ dual
    L = 0._mytype
    M = 0._mytype
    N = 0._mytype
    !! On calcule maxR : la valeur maximale sur tous les residus stockes
    maxR = maxval(abs(ACT3_R(:,:)))
    !! On renormalise toutes les sommes en divisant par maxR pour eviter d'atteindre 
    !! la valeur maximale d'un reel 

    f12 = SUM(((ACT3_R(:,2)-ACT3_R(:,1))/maxR)**2)
    f22 = SUM(((ACT3_R(:,3)-ACT3_R(:,1))/maxR)**2)
    f32 = SUM(((ACT3_R(:,4)-ACT3_R(:,1))/maxR)**2)
    if(f12.LT.XPETIT*100._mytype) f12=XPETIT*100._mytype
    if(f22.LT.XPETIT*100._mytype) f22=XPETIT*100._mytype
    if(f32.LT.XPETIT*100._mytype) f32=XPETIT*100._mytype

    !!Détection de flip-flop ?
    !!(1 et 3 sont trop proches, malgré 2 différent de 1) 
    if (f12/f22.gt.1.2_mytype .and. f12/f32.gt.0.707_mytype .and.&
         f12/f32.lt.1.414_mytype) then

       !       call amitex_abort("act3 : plan B",-1,0)   !!! A Remettre lors de l'intégration à AMITEX !!!
       !       ! CALL ERREUR(858)
       L = (sqrt(f12)+sqrt(f22))/(2*sqrt(f12))
       M = 0._mytype
       N = 0._mytype
       GOTO 1010
    endif
    !!Construction d'une famille libre orthogonale
    F1F2 = SUM((ACT3_R(:,2)-ACT3_R(:,1))/maxR*(ACT3_R(:,3) -ACT3_R(:,1))/maxR)
    F1F3 = SUM((ACT3_R(:,2)-ACT3_R(:,1))/maxR*(ACT3_R(:,4) -ACT3_R(:,1))/maxR)
    F2F3 = SUM((ACT3_R(:,3)-ACT3_R(:,1))/maxR*(ACT3_R(:,4) -ACT3_R(:,1))/maxR)
    F0F1 = SUM(ACT3_R(:,1)/maxR*(ACT3_R(:,2)-ACT3_R(:,1))/maxR)
    F0F2 = SUM(ACT3_R(:,1)/maxR*(ACT3_R(:,3)-ACT3_R(:,1))/maxR)
    F0F3 = SUM(ACT3_R(:,1)/maxR*(ACT3_R(:,4)-ACT3_R(:,1))/maxR)


    A1  =  F12
    B1  =  F1F2
    C1  =  F1F3
    J1 = -F0F1
    D1  =  F1F2
    E1  =  F22
    F1  =  F2F3
    K1 = -F0F2
    G1  =  F1F3
    H1  =  F2F3
    I1 =  F32
    L1 = -F0F3
    !!write (6,*) ' f12 f22 f32 ',f12,f22,f32

    !!Calcul des coordonnées de IFO0 dans la famille libre
    !!et élimination des vecteurs liés
    IF (ABS(I1) .LT. xpetit) GOTO 2000

    A2  = A1  - ( G1*C1/I1)
    B2  = B1  - ( H1*C1/I1)
    J2 = J1 - (L1*C1/I1)
    D2  = D1  - ( G1*F1/I1)
    E2  = E1  - ( H1*F1/I1)
    K2 = K1 - (L1*F1/I1)

    IF (ABS(E2) .LE. ABS(E1)*1E-6_mytype) THEN
       !!      write (6,*) ' act3 e2 e1 '
       GOTO 2000
    ENDIF

    A3  = A2 -(D2*B2/E2)
    J3 = J2-(K2*B2/E2)

    IF (ABS(A3) .LE. ABS(A1)*1E-6_mytype) THEN
       !!      write (6,*) ' act3 a3 a2 '
       GOTO 2000
    ENDIF
    IF (ABS(A2) .LE. ABS(A1)*1E-6_mytype) THEN
       !!      write (6,*) ' act3 a2 a1 '
       GOTO 2000
    ENDIF

    L=J3/A3

    B4  = B2
    J4 = J2 - (A2*L)

    IF (ABS(B4) .LE. ABS(B1)*1E-6_mytype) THEN
       !!      write (6,*) ' act3 b4 b1 '
       GOTO 2000
    ENDIF

    M  = J4/B4

    C5  = C1
    J5 = J1 - (A1*L) - (B1*M)

    IF (ABS(C5) .LT. xpetit) GOTO 2000

    N=J5/C5

    ref=sqrt(a1**2+b1**2+c1**2)
    !!write (6,*) ' act3 ',ref,XI1,A3,b4,c5

    GOTO 1010

2000 CONTINUE
    !!le systeme est singulier on essaye avec un vecteur de moins

    !!write (6,*) ' act3 ',xi1,e2,e1,a3,a1,b4,b1,c5
    IF (ABS(E1) .LT. xpetit) GOTO 2100

    A2  = A1  - (D1 *B1/E1)
    J2 = J1 - (K1*B1/E1)

    IF (ABS(A2) .LE. ABS(A1)*1E-6_mytype) THEN
       !!   write (6,*) ' act3 a2 a1 '
       GOTO 2100
    ENDIF

    L  = J2/A2

    B3  = B1
    J3 = J1 - (A1*L)

    IF (ABS(B3) .LT. xpetit) GOTO 2100
    !!
    M = J3/B3

    N = 0._mytype

    !!act3 : reduction a 2 dimensions
    !    call amitex_abort("act3 : reduction a 2 dimensions",-1,0)      !!! A Remettre lors de l'intégration à AMITEX !!!

    GOTO 1010

2100 CONTINUE
    !!le systeme est singulier on essaye avec un vecteur de moins
    IF (ABS(A1) .LT. xpetit ) GOTO 1000

    L  = J1/A1

    M  = 0._mytype
    N  = 0._mytype

    !!act3 : reduction a 1 dimensions
    !    call amitex_abort("act3 : reduction a 1 dimension",-1,0)  !!! A Remettre lors de l'intégration à AMITEX !!!
    GOTO 1010

1000 CONTINUE
    !!act3 : pas d accélération
    !    call amitex_abort("act3 : pas d'acceleration",-1,0)    !!! A Remettre lors de l'intégration à AMITEX !!!
    L = 0._mytype
    M = 0._mytype
    N = 0._mytype

1010 CONTINUE


    !!Vérification qu'on ne retombe pas sur l'un des duaux entrés
    !!Sinon, on réessaye en tenant compte de moins de vecteurs
    if (abs(l-1._mytype).lt.1e-1_mytype .and. abs(m).lt.1e-1_mytype&
         .and. abs(n).lt.1e-1_mytype) goto 1000
    if (abs(m-1._mytype).lt.1e-1_mytype .and. abs(n).lt.1e-1_mytype&
         .and. abs(l).lt.1e-1_mytype) goto 2100
    if (abs(n-1._mytype).lt.1e-1_mytype .and. abs(l).lt.1e-1_mytype&
         .and. abs(m).lt.1e-1_mytype) goto 2000

    Def = L * (ACT3_U(:,2) - ACT3_U(:,1)) + M *(ACT3_U(:,3) - ACT3_U(:,1)) + N  * (ACT3_U(:,4) - ACT3_U(:,1)) ;
    Def = Def + ACT3_U(:,1);
    RETURN

  end subroutine ACT3LOCAL
!====================================================================
  
!====================================================================
!===================================================================
!
!        FONCTION ACT3_Moy
!
!> \brief Propose une solution au solveur iteratif en utilisant 
!! les quatre calculs precedents.\n
!!
!! Extrait notice CAST3M:\n
!! REND UN CHAMP PRIMAL EXTRAPOLE QUI MINIMISE LE CHAMP DUAL.
!! METHODE : ON ESTIME L'APPLICATION TANGENTE ET ON L'UTILISE
!! POUR RESOUDRE LE PROBLEME
!!
!! \param[in] ACT3_U: Deformation moyenne imposee
!!                 tableaux de réels (nb_stress_pilot,4)
!! \param[in] ACT3_R: Contrainte moyenne qu'on voudrait imposee - contrainte moyenne obtenue
!! \param[in,out] Eimp: Deformation moyenne proposee 
!!
!-----------------------------------------------------------------
  subroutine ACT3_Moy(ACT3_U,ACT3_R,Eimp)

    implicit none
   
    !> CONTRAINTES ET DEFORMATIONS
    !! dans l'espace reel, coherent avec la structure "MATERIAL"
    real(mytype),dimension(:,:), intent(in)  :: ACT3_R, ACT3_U
    real(mytype),dimension(:),intent(inout)  :: Eimp
    real(mytype)                               :: L,M,N
    real(mytype)                               :: F12,F22,F32,F1F2,F1F3,&
         F2F3,F0F1,F0F2,F0F3,&
         A1,B1,C1,J1,D1,E1,F1,K1,G1,H1,I1,L1,A2,B2,J2,D2,E2,K2,A3,B3,J3,B4,J4,C5,J5,&
         ref
    !! Tolerance epsilon differente en simple et double precision
#ifdef DOUBLE_PREC
    real(mytype),parameter                     :: XPETIT = 1e-30
#else
    real(mytype),parameter                     :: XPETIT = 1e-20
#endif

    !! INITIALISATIONS (potentiellement non-initialises d'apres gcc -Wall)
    K1 = 0._mytype
    J1 = 0._mytype
    E1 = 0._mytype
    D1 = 0._mytype
    B1 = 0._mytype
    A1 = 0._mytype
    
    !!CAS  A 4 POINTS :
    !!-----------------
    !!famille des 3 vecteurs variations du champ dual
    L = 0._mytype
    M = 0._mytype
    N = 0._mytype

    F12 = SUM(((ACT3_R(:,2)-ACT3_R(:,1)))**2)
    F22 = SUM(((ACT3_R(:,3)-ACT3_R(:,1)))**2)
    F32 = SUM(((ACT3_R(:,4)-ACT3_R(:,1)))**2)
    if(f12.LT.XPETIT*100._mytype) f12=XPETIT*100._mytype
    if(f22.LT.XPETIT*100._mytype) f22=XPETIT*100._mytype
    if(f32.LT.XPETIT*100._mytype) f32=XPETIT*100._mytype

    !!Détection de flip-flop ?
    !!(1 et 3 sont trop proches, malgré 2 différent de 1) 
    if (f12/f22.gt.1.2_mytype .and. f12/f32.gt.0.707_mytype .and.&
         f12/f32.lt.1.414_mytype) then
       call amitex_abort("act3 : plan B",-1,0)
       ! CALL ERREUR(858)
       L = (sqrt(f12)+sqrt(f22))/(2*sqrt(f12))
       M = 0._mytype
       N = 0._mytype
       GOTO 1010
    endif
    !!Construction d'une famille libre orthogonale
    F1F2 = SUM((ACT3_R(:,2)-ACT3_R(:,1))*(ACT3_R(:,3) -ACT3_R(:,1)))
    F1F3 = SUM((ACT3_R(:,2)-ACT3_R(:,1))*(ACT3_R(:,4) -ACT3_R(:,1)))
    F2F3 = SUM((ACT3_R(:,3)-ACT3_R(:,1))*(ACT3_R(:,4) -ACT3_R(:,1)))
    F0F1 = SUM(ACT3_R(:,1)*(ACT3_R(:,2)-ACT3_R(:,1)))
    F0F2 = SUM(ACT3_R(:,1)*(ACT3_R(:,3)-ACT3_R(:,1)))
    F0F3 = SUM(ACT3_R(:,1)*(ACT3_R(:,4)-ACT3_R(:,1)))

    A1  =  F12
    B1  =  F1F2
    C1  =  F1F3
    J1 = -F0F1
    D1  =  F1F2
    E1  =  F22
    F1  =  F2F3
    K1 = -F0F2
    G1  =  F1F3
    H1  =  F2F3
    I1 =  F32
    L1 = -F0F3
    !!write (6,*) ' f12 f22 f32 ',f12,f22,f32
    
    !!Calcul des coordonnées de IFO0 dans la famille libre
    !!et élimination des vecteurs liés
    IF (ABS(I1) .LT. xpetit) GOTO 2000

    A2  = A1  - ( G1*C1/I1)
    B2  = B1  - ( H1*C1/I1)
    J2 = J1 - (L1*C1/I1)
    D2  = D1  - ( G1*F1/I1)
    E2  = E1  - ( H1*F1/I1)
    K2 = K1 - (L1*F1/I1)

    IF (ABS(E2) .LE. ABS(E1)*1E-6_mytype) THEN
       !!      write (6,*) ' act3 e2 e1 '
       GOTO 2000
    ENDIF

    A3  = A2 -(D2*B2/E2)
    J3 = J2-(K2*B2/E2)

    IF (ABS(A3) .LE. ABS(A1)*1E-6_mytype) THEN
       !!      write (6,*) ' act3 a3 a2 '
       GOTO 2000
    ENDIF
    IF (ABS(A2) .LE. ABS(A1)*1E-6_mytype) THEN
       !!      write (6,*) ' act3 a2 a1 '
       GOTO 2000
    ENDIF

    L=J3/A3

    B4  = B2
    J4 = J2 - (A2*L)

    IF (ABS(B4) .LE. ABS(B1)*1E-6_mytype) THEN
       !!      write (6,*) ' act3 b4 b1 '
       GOTO 2000
    ENDIF

    M  = J4/B4

    C5  = C1
    J5 = J1 - (A1*L) - (B1*M)

    IF (ABS(C5) .LT. xpetit) GOTO 2000

    N=J5/C5

    ref=sqrt(a1**2+b1**2+c1**2)
    !!write (6,*) ' act3 ',ref,XI1,A3,b4,c5

    GOTO 1010

2000 CONTINUE
    !!le systeme est singulier on essaye avec un vecteur de moins

    !!write (6,*) ' act3 ',xi1,e2,e1,a3,a1,b4,b1,c5
    IF (ABS(E1) .LT. xpetit) GOTO 2100

    A2  = A1  - (D1 *B1/E1)
    J2 = J1 - (K1*B1/E1)

    IF (ABS(A2) .LE. ABS(A1)*1E-6_mytype) THEN
       !!   write (6,*) ' act3 a2 a1 '
       GOTO 2100
    ENDIF

    L  = J2/A2

    B3  = B1
    J3 = J1 - (A1*L)

    IF (ABS(B3) .LT. xpetit) GOTO 2100
    !!
    M = J3/B3

    N = 0._mytype

    !!act3 : reduction a 2 dimensions
    call amitex_abort("act3 : reduction a 2 dimensions",-1,0)

    GOTO 1010

2100 CONTINUE
    !!le systeme est singulier on essaye avec un vecteur de moins
    IF (ABS(A1) .LT. xpetit ) GOTO 1000

    L  = J1/A1

    M  = 0._mytype
    N  = 0._mytype

    !!act3 : reduction a 1 dimensions
    call amitex_abort("act3 : reduction a 1 dimension",-1,0)
    GOTO 1010

1000 CONTINUE
    !!act3 : pas d accélération
    call amitex_abort("act3 : pas d'acceleration",-1,0)
    L = 0._mytype
    M = 0._mytype
    N = 0._mytype

1010 CONTINUE


    !!Vérification qu'on ne retombe pas sur l'un des duaux entrés
    !!Sinon, on réessaye en tenant compte de moins de vecteurs
    if (abs(l-1._mytype).lt.1e-1_mytype .and. abs(m).lt.1e-1_mytype&
         .and. abs(n).lt.1e-1_mytype) goto 1000
    if (abs(m-1._mytype).lt.1e-1_mytype .and. abs(n).lt.1e-1_mytype&
         .and. abs(l).lt.1e-1_mytype) goto 2100
    if (abs(n-1._mytype).lt.1e-1_mytype .and. abs(l).lt.1e-1_mytype&
         .and. abs(m).lt.1e-1_mytype) goto 2000

    Eimp = L * (ACT3_U(:,2) - ACT3_U(:,1)) + M *(ACT3_U(:,3) - ACT3_U(:,1)) + N  * (ACT3_U(:,4) - ACT3_U(:,1))
    Eimp = Eimp + ACT3_U(:,1)
    RETURN
    
  end subroutine ACT3_Moy
!====================================================================





!====================================================================

!====================================================================
!>       SUBROUTINE DE CALCUL DU CHAMP DE CONTRAINTE DE PIOLA-KIRCHHOFF
!!       A PARTIR DU CHAMP DE CONTRAINTE DE CAUCHY ET DU GRADIENT DES DEFORMATIONS
!!
!! NOTATIONS (entree et champs) Sigma : CAST3M (ordre 11 22 33 12 13 23)
!!            PK1 et grad(u) : ordre 11 22 33 12 13 23 21 31 32 sans notation
!!
!! \param[out] PK1   Champ des contraintes de Piola-Kirchhoff
!! \param[in] Sig    Champ des contraintes de Cauchy
!! \param[in] gradU  Champ de gradient des deplacements
!!
!====================================================================
subroutine SigToPK1(Sig,gradU,PK1)

  implicit none
  real(mytype),dimension(:,:),intent(in)      :: gradU,Sig
  real(mytype),dimension(:,:),intent(inout)   :: PK1
  real(mytype),dimension(3,3)                 :: MatSig,MatStress,COFACTOR
  integer                                     :: i
  !-----------PARALLELISATION OMP
#ifdef OPENMP
  integer                                     :: nbthread
!$OMP PARALLEL
  nbthread = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL

!$OMP PARALLEL NUM_THREADS(nbthread) private(COFACTOR,MatSig,MatStress)
!$OMP DO
#endif
  do i = 1,size(gradU,1)
     !! On calcule la matrice des cofacteurs de Id + gradU
     COFACTOR(1,1) = +((1+ gradU(i,2))*(1+ gradU(i,3))-gradU(i,6) &
          *gradU(i,9))
     COFACTOR(1,2) = -(gradU(i,7)*(1+ gradU(i,3))-gradU(i,6)*gradU(i,8))
     COFACTOR(1,3) = +(gradU(i,7)*gradU(i,9)-(1+ gradU(i,2))*gradU(i,8))
     COFACTOR(2,1) = -(gradU(i,4)*(1+ gradU(i,3))-gradU(i,5)*gradU(i,9)) 
     COFACTOR(2,2) = +((1+ gradU(i,1))*(1+ gradU(i,3))-gradU(i,5)&
          *gradU(i,8))
     COFACTOR(2,3) = -((1+ gradU(i,1))*gradU(i,9)-gradU(i,4)*gradU(i,8))
     COFACTOR(3,1) = +(gradU(i,4)*gradU(i,6)-gradU(i,5)*(1+ gradU(i,2)))
     COFACTOR(3,2) = -((1+ gradU(i,1))*gradU(i,6)-gradU(i,5)*gradU(i,7))
     COFACTOR(3,3) = +((1+ gradU(i,1))*(1+ gradU(i,2))-gradU(i,4)&
          *gradU(i,7))

     !---- On met sigma sous forme de matrice symetrique dans MatSig
     !---- en prenant en compte la notation
     MatSig(1,1) = Sig(i,1)
     MatSig(2,2) = Sig(i,2)
     MatSig(3,3) = Sig(i,3)
     MatSig(1,2) = Sig(i,4)
     MatSig(1,3) = Sig(i,5)
     MatSig(2,3) = Sig(i,6)
     MatSig(2,1) = Sig(i,4)
     MatSig(3,1) = Sig(i,5)
     MatSig(3,2) = Sig(i,6)
     !---- PK1 = Det(F) * sig * TRANSPOSE(F^{-1}) = sig * cofactor
     MatStress=MATMUL(MatSig,COFACTOR)
     !---- On met PK1 sous forme de vecteur avec la notation 11 22 33 12 13 23 21 31 32
     PK1(i,1) = MatStress(1,1)
     PK1(i,2) = MatStress(2,2)
     PK1(i,3) = MatStress(3,3)
     PK1(i,4) = MatStress(1,2)
     PK1(i,5) = MatStress(1,3)
     PK1(i,6) = MatStress(2,3)
     PK1(i,7) = MatStress(2,1)
     PK1(i,8) = MatStress(3,1)
     PK1(i,9) = MatStress(3,2)
  end do
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif

end subroutine SigToPK1
!====================================================================

!====================================================================
!>       SUBROUTINE DE CALCUL DU CHAMP DE CONTRAINTE DE CAUCHY
!!       A PARTIR DU CHAMP DE CONTRAINTE DE PIOLA-KIRCHHOFF ET DU GRADIENT DES DEFORMATIONS
!!
!! NOTATIONS (entree et champs) Sigma : CAST3M (ordre 11 22 33 12 13 23)
!!            PK1 et grad(u) : ordre 11 22 33 12 13 23 21 31 32 sans notation
!!
!! \param[in] PK1   Champ des contraintes de Piola-Kirchhoff
!! \param[out] Sig    Champ des contraintes de Cauchy
!! \param[in] gradU  Champ de gradient des deplacements
!!
!====================================================================
subroutine PK1ToSig(Sig,gradU,PK1)

  implicit none
  real(mytype),dimension(:,:),intent(in)      :: gradU,PK1
  real(mytype),dimension(:,:),intent(inout)   :: Sig
  real(mytype),dimension(3,3)                 :: MatSig,MatStress,tMatF
  real(mytype)                                :: invDET,W1,W2,W3
  integer                                     :: i
  !-----------PARALLELISATION OMP
#ifdef OPENMP
  integer                                     :: nbthread
!$OMP PARALLEL
  nbthread = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL

!$OMP PARALLEL NUM_THREADS(nbthread) private(tMatF,W1,W2,W3,invDET,MatSig,MatStress)
!$OMP DO
#endif
  do i = 1,size(gradU,1)
     !! On calcule la matrice transposee de F = Id + grad(u)
     tMatF(1,1) = 1._mytype + gradU(i,1)
     tMatF(2,2) = 1._mytype + gradU(i,2)
     tMatF(3,3) = 1._mytype + gradU(i,3)
     tMatF(2,1) = gradU(i,4)
     tMatF(3,1) = gradU(i,5)
     tMatF(3,2) = gradU(i,6)
     tMatF(1,2) = gradU(i,7)
     tMatF(1,3) = gradU(i,8)
     tMatF(2,3) = gradU(i,9)
     !!   on calcule l'inverse du determinant de F egal a celui de sa transposee
     W1 = tMatF(2,2)*tMatF(3,3) - tMatF(3,2)*tMatF(2,3)
     W2 = tMatF(3,1)*tMatF(2,3) - tMatF(2,1)*tMatF(3,3)
     W3 = tMatF(2,1)*tMatF(3,2) - tMatF(3,1)*tMatF(2,2)
     invDET = 1. / (tMatF(1,1)*W1 + tMatF(1,2)*W2 + tMatF(1,3)*W3)

     !---- On met PK1 sous forme de matrice symetrique dans MatStress
     !---- en prenant en compte la notation 11 22 33 12 13 23 21 31 32
     MatStress(1,1) = PK1(i,1)
     MatStress(2,2) = PK1(i,2)
     MatStress(3,3) = PK1(i,3)
     MatStress(1,2) = PK1(i,4)
     MatStress(1,3) = PK1(i,5)
     MatStress(2,3) = PK1(i,6)
     MatStress(2,1) = PK1(i,7)
     MatStress(3,1) = PK1(i,8)
     MatStress(3,2) = PK1(i,9)
     !---- sig = 1/Det(F) * PK1 * TRANSPOSE(F)
     MatSig=invDET*MATMUL(MatStress,tMatF)
     !---- On met sigma sous forme de vecteur avec la notation 11 22 33 12 13 23
     Sig(i,1) = MatSig(1,1)
     Sig(i,2) = MatSig(2,2)
     Sig(i,3) = MatSig(3,3)
     Sig(i,4) = MatSig(1,2)
     Sig(i,5) = MatSig(1,3)
     Sig(i,6) = MatSig(2,3)
  end do
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif

end subroutine PK1ToSig
!====================================================================

!====================================================================
!>       EVALUATION DE LA DIRECTION DE CONTRAINTE PK1
!!       A PARTIR DE LA DIRECTION DE CONTRAINTE DE CAUCHY
!!
!! NOTATIONS DirStress, DirStress2 et gradU : 
!!           ordre 11 22 33 12 13 23 21 31 32 sans notation
!!
!! \param[out] DirStress2 Direction de chargement PK1 (tableau 9)
!! \param[in]  DirStress  Direction de chargement Cauchy (tableau 9)
!! \param[in]  gradU      Gradient de deplacement (tableau 9)
!!
!! ATTENTION :  dans ce contexte, la contrainte de Cauchy est exprime
!!              avec ses 9 composantes (on aura pris soin au prealable 
!!              de verifier que les cisaillements sont symetriques)
!====================================================================
subroutine DirStress_SigToPK1(DirStress,DirStress2,gradU)

  implicit none
  real(mytype),dimension(9),intent(in)      :: gradU,DirStress
  real(mytype),dimension(9),intent(out)     :: DirStress2
  real(mytype),dimension(3,3)               :: MatSig,MatStress,COFACTOR

     !! On calcule la matrice des cofacteurs de Id + gradU
     COFACTOR(1,1) = +((1+ gradU(2))*(1+ gradU(3))-gradU(6) &
          *gradU(9))
     COFACTOR(1,2) = -(gradU(7)*(1+ gradU(3))-gradU(6)*gradU(8))
     COFACTOR(1,3) = +(gradU(7)*gradU(9)-(1+ gradU(2))*gradU(8))
     COFACTOR(2,1) = -(gradU(4)*(1+ gradU(3))-gradU(5)*gradU(9)) 
     COFACTOR(2,2) = +((1+ gradU(1))*(1+ gradU(3))-gradU(5)&
          *gradU(8))
     COFACTOR(2,3) = -((1+ gradU(1))*gradU(9)-gradU(4)*gradU(8))
     COFACTOR(3,1) = +(gradU(4)*gradU(6)-gradU(5)*(1+ gradU(2)))
     COFACTOR(3,2) = -((1+ gradU(1))*gradU(6)-gradU(5)*gradU(7))
     COFACTOR(3,3) = +((1+ gradU(1))*(1+ gradU(2))-gradU(4)&
          *gradU(7))

     !---- On met sigma sous forme de matrice symetrique dans MatSig
     !---- en prenant en compte la notation
     MatSig(1,1) = DirStress(1)
     MatSig(2,2) = DirStress(2)
     MatSig(3,3) = DirStress(3)
     MatSig(1,2) = DirStress(4)
     MatSig(1,3) = DirStress(5)
     MatSig(2,3) = DirStress(6)
     MatSig(2,1) = DirStress(7)
     MatSig(3,1) = DirStress(8)
     MatSig(3,2) = DirStress(9)
     !---- PK1 = Det(F) * sig * TRANSPOSE(F^{-1}) = sig * cofactor
     MatStress=MATMUL(MatSig,COFACTOR)
     !---- On met PK1 sous forme de vecteur avec la notation 11 22 33 12 13 23 21 31 32
     DirStress2(1) = MatStress(1,1)
     DirStress2(2) = MatStress(2,2)
     DirStress2(3) = MatStress(3,3)
     DirStress2(4) = MatStress(1,2)
     DirStress2(5) = MatStress(1,3)
     DirStress2(6) = MatStress(2,3)
     DirStress2(7) = MatStress(2,1)
     DirStress2(8) = MatStress(3,1)
     DirStress2(9) = MatStress(3,2)


end subroutine DirStress_SigToPK1
!====================================================================





!================================================================================
!>      optimize_scalar(i_BC, Tc_in, Tc_res, tol_T_BC,Tc_norm)
!! 
!!      optimize Tc_in(0) to vanish Tc_res(0)
!!      for an input TC_in(0), Tc_res(0) output a new proposal for Tc_in(0)
!!      uses previous values in Tc_in(1:2), Tc_res(1:2)
!!          
!!  \param[in]    i_BC        iteration index
!!  \param[inout] Tc_in       solution for Tc_in applied to the field at t
!!                            dimension(0:2)
!!  \param[inout] Tc_res      residual (Tc-Tc_after) associated to Tc_in
!!                            dimension(0:2)
!!  \param[in]    tol_T_BC    tolerance on normalized error (Tc_res / Tc_norm)
!!  \patam[in]    Tc_norm     normalization (generally norm(DT))
!!
!! ALGO :
!!--------- 
!!   for first iterations : Tc_in(0) = Tc_in(0) + 1*(Tc_res(0))
!!   for next iterations : 
!!       optimize Tc_in(0) from Tc_res(0) and previous couples (solution, residual)  
!!                                                   i.e. (Tc_in(1:2), Tc_res(1:2)) 
!!       if (Tc_in(2)-Tc_in(1) too low, uses :  Tc_in(0) = Tc_in(0) + 1*(Tc_res(0))
!!
!================================================================================
subroutine optimize_scalar(i_BC, Tc_in, Tc_res, tol_T_BC,Tc_norm)

  implicit none
  
  integer, intent(in)                        :: i_BC
  real(mytype), dimension(0:2),intent(inout) :: Tc_in, Tc_res
  real(mytype),intent(in)                    :: tol_T_BC, Tc_norm
  
  
  logical       :: acv2
  real(mytype)  :: p0, m0
 
  ! init local variables 
  acv2 = .false.

  if (i_BC<=0) then
    call amitex_abort("optimize_BC_translation : iteration index i_BC must be > 0 ",2,0)
  end if
  
  ! ALGO 
  !   for first iterations : Tc_in(0) = Tc_in(0) + 1*(Tc_res(0))
  !   for next iterations : 
  !       optimize Tc_in(0) from Tc_res(0) and previous (solution, residual)  (Tc_in(1:2), Tc_res(1:2)) 
  !       if (Tc_in(2)-Tc_in(1) too low, uses :  Tc_in(0) = Tc_in(0) + 1*(Tc_res(0))
  !
     if (i_BC == 1) then
        Tc_in(1)  = Tc_in(0)
        Tc_res(1) = Tc_res(0)
        Tc_in(0)  = Tc_in(0) + 1*(Tc_res(0)) ! ce facteur "1" semble pouvoir etre optimise mais depend sans doute de bcp de choses...
                                             ! valeur "1" est sans doute + sure
     elseif (i_BC == 2) then
        Tc_in(2)  = Tc_in(0)
        Tc_res(2) = Tc_res(0)
        Tc_in(0)  = Tc_in(0) + 1*(Tc_res(0)) ! ce facteur "1" semble pouvoir etre optimise mais depend sans doute de bcp de choses...
                                             ! valeur "1" est sans doute + sure
     else
        m0 = (Tc_res(2) - Tc_res(1))/ (Tc_in(2) - Tc_in(1))
        p0 = Tc_res(2) - m0 * Tc_in(2)
        if ((abs(Tc_in(2) - Tc_in(1))/Tc_norm)> 2.*tol_T_BC) acv2 = .true.
        if (mod(i_BC,2)/=0) then
          Tc_in(1)  = Tc_in(0)
          Tc_res(1) = Tc_res(0)
        else
          Tc_in(2)  = Tc_in(0)
          Tc_res(2) = Tc_res(0)
        end if
        ! if (Tc_in2 - Tc_in1) is too low , use the classical fix-point
        if (acv2) then 
           Tc_in(0) = - p0 / m0
        else
           Tc_in(0) = Tc_in(0) + 1*(Tc_res(0))! ce facteur "1" semble pouvoir etre optimise mais depend sans doute de bcp de choses...
                                              ! valeur "1" est sans doute + sure
        end if
     end if     
     
end subroutine optimize_scalar


end module algo_functions_mod
