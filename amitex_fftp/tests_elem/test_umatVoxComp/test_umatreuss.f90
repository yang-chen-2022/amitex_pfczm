program  test_umatreuss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> CHARGEMENT DES MODULES
  use ISO_C_BINDING
  use material_mod
  use param_algo_mod

!!==============================================================================
!!						      		    DECLARATIONS
!!==============================================================================
  implicit none
!-----------------------------------------------------------------------------------------------
!         PENSER A AJUSTER LA DECLARATION varInt(nvarint) et coeff(ncoeff) en fonction de la loi
!-----------------------------------------------------------------------------------------------
!----------------------------------------------- PARAMETRES DE UMAT
  integer,parameter                :: nphases =2        ! A AJUSTER
  double precision, dimension(6)   :: sig,def,ddef,def0
  double precision                 :: temp,dtemp
  double precision, dimension(3,3) :: DFGRD0, DFGRD1
  double precision, dimension(6*nphases) ::  coeff      !umatReuss : 6*nb_phases en elasticite
  double precision, dimension(1*nphases + 12*nphases)   ::  varInt   !umatReuss : (1+12)*nb_phases en elasticite
  double precision, dimension(1)   :: predef
  integer(kind=8)                  :: kinc, ncoeff,nvarint, KSTEP
  double precision                 :: t,dt 
  double precision, dimension(2)   :: dt_tab
!--------------------------------------------- Tableau de pointeurs de fonctions
  type(c_ptr),dimension(nphases)   :: umatptr    ! 2 phases
!---------------------------------------------- AUTRES PARAMETRES

  character(len=200,kind=C_CHAR)   :: lawname
  character(len=200,kind=C_CHAR)   :: libname

  real     :: t1,t2
  integer  :: i,Niter, cas



!==========================================================================

  t1=0.
  t2=0.
  call cpu_time(t1)


!                                 Choix de la loi
!------------------------------------------------
!

 libname='/home/gelebart/amitex_fftp/libAmitex/src/materiaux/libUmatAmitex.so'
 lawname = 'elasiso'
 print *, "libname, lawname", trim(libname),' ', trim(lawname)

 umatptr(1) = umat_load(trim(libname)//C_NULL_CHAR,trim(lawname)//C_NULL_CHAR)
 umatptr(2) = umat_load(trim(libname)//C_NULL_CHAR,trim(lawname)//C_NULL_CHAR)


! Definition du chargement en deformation imposee
!------------------------------------------------
  t=0.           ! temps initial
  Niter = 1      ! nombre d'increments de temps
  dt=100._8/Niter  ! increment de temps 

  temp = 293.15_8  ! temperature initiale
  dtemp = 0.     ! increment de temperature

  def0 =0.       ! deformation initiale

  ddef(1)=-0.0025_8/Niter ! increments de deformation (6 composantes, notation de Voigt)
  ddef(2)=-0.0025_8/Niter
  ddef(3)=0.01_8/Niter
  ddef(4)=0.0
  ddef(5)=0.0
  ddef(6)=0.0
  
!           Initialisation des variables internes
!------------------------------------------------
   nvarint=1*nphases + 12*nphases
   varint = 0.
   varint(1) = 0	! nb varint de la phase 1
   varint(2) = 0        ! nb varint de la phase 2

!                                 Coeff materiaux
!------------------------------------------------
   ncoeff = 6*nphases !pour 2 phases elastiques

   cas = 2

   ! cas homogene
   if (cas .eq. 1) then
   coeff(1) =  2 ! nb_phases
   coeff(2) =  2 ! nb coeff phase 1
   coeff(3) =  2 ! nb coeff phase 2
   coeff(4) =  0.5_8    ! fv phase 1
   coeff(5) =   100000. ! lambdaeq phase 1
   coeff(6) =   50000.  ! Mueq phase 1
   coeff(7) =   100000. ! lambda phase 1
   coeff(8) =   50000.  ! Mu phase 1
   coeff(9) =   100000. ! lambdaeq phase 2
   coeff(10) =  50000.  ! Mueq phase 2
   coeff(11) =  100000. ! lambda phase 2
   coeff(12) =  50000.  ! Mu phase 2

   ! cas heterogene
   else if (cas .eq. 2) then
   coeff(1) =  2 ! nb_phases
   coeff(2) =  2 ! nb coeff phase 1
   coeff(3) =  2 ! nb coeff phase 2
   coeff(4) =  0.1_8     ! fv phase 1
   coeff(5) =   1000000. ! lambdaeq phase 1
   coeff(6) =   500000.  ! Mueq phase 1
   coeff(7) =   1000000. ! lambda phase 1
   coeff(8) =   500000.  ! Mu phase 1
   coeff(9) =   100000.  ! lambdaeq phase 2
   coeff(10) =  50000.   ! Mueq phase 2
   coeff(11) =  100000.  ! lambda phase 2
   coeff(12) =  50000.  ! Mu phase 2
   end if
print *,"coeff = ", coeff

! INITIALISATION POUR UTILISATION UMAT      
! HYPER IMPORTANT (sinon NTENS n'est pas connu)
!                                 (non utilisees)
!------------------------------------------------
  call initParam_umat() ! TRES IMPORTANT
  KINC = 1
  KSTEP = 0
  DFGRD0 = 0.
  DFGRD1 = 0.
  PREDEF = 0.

!
! CHOIX DES PARAMETRES D'INTEGRATION DU MODELE DE REUSS
!--------------------------------------------------------------
  algo_param%tol_criteq = 1e-4_8

  algo_param%tol_criteq_reuss = 1e-2_8
  algo_param%acc_CV_reuss = .false.
  algo_param%Init_reuss = "TODO"
  algo_param%Jacobian_type_reuss="elastic"

!-------------------------------- APPEL INCREMENTAL A LA LOI
  sig = 0.
  def0 = 0.

  write (*,"(A)"), "# colonne 1    : temps"
  write (*,"(A)"), "# colonne 2-7  : deformation"
  write (*,"(A)"), "# colonne 8-13 : contrainte"
  write (*,"(13(1X,E15.8))"), t, def0,sig
  do i=1,Niter
    t=t+dt
    def=def0+ddef
    dt_tab = 0
    dt_tab(1) = dt
    call umatReuss(umatptr,sig,varint,def0,ddef,(/0._8,t-dt/),dt_tab,temp,PREDEF,&
                     nvarint,coeff,ncoeff,DFGRD0, DFGRD1,KSTEP,KINC)

    def0=def
    write (*,"(13(1X,E15.8))"), t, def0,sig
  end do
!--------------------------------

  call cpu_time(t2)
  write (*,"(A,E15.8)"),  "#temps ecoule", t2-t1




!==============================================================================
!  VERIFICATION MODELE DE REUSS
!  compatibilite
print *,"-----------"
print *, "Def_impose = ", def0
print *, "Defmoy = ", 0.1_8*varint(3:8) + 0.9_8*varint(15:20)
print *,"-----------"
print *, "Sigmoy = ", sig
print *,"-----------"
print * ,'sig1 = ', coeff(7)*sum(varint(3:5))+2*coeff(8)*varint(3)
print * ,'sig2 = ', coeff(7)*sum(varint(3:5))+2*coeff(8)*varint(4)
print * ,'sig3 = ', coeff(7)*sum(varint(3:5))+2*coeff(8)*varint(5)
print *,"-----------"
print * ,'sig1 = ', coeff(11)*sum(varint(15:17))+2*coeff(12)*varint(15)
print * ,'sig2 = ', coeff(11)*sum(varint(15:17))+2*coeff(12)*varint(16)
print * ,'sig3 = ', coeff(11)*sum(varint(15:17))+2*coeff(12)*varint(17)


end program test_umatreuss
