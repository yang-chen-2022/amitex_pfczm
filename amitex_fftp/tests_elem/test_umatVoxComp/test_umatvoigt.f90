program  test_umatvoigt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> CHARGEMENT DES MODULES
  use ISO_C_BINDING
  use material_mod

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
  double precision, dimension(4*nphases) ::  coeff      !umatVoigt : 4*nb_phases en elasticite
  double precision, dimension(1*nphases + 6*nphases)   ::  varInt   !umatVoigt : (1+6)*nb_phases en elasticite
  double precision, dimension(1)   :: predef
  integer(kind=8)                  :: kinc, ncoeff,nvarint, KSTEP
  double precision :: t,dt
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
   nvarint=nphases + 6*nphases
   varint = 0.
   varint(1) = 0	! nb varint de la phase 1
   varint(2) = 0        ! nb varint de la phase 2

!                                 Coeff materiaux
!------------------------------------------------
   ncoeff = 4*nphases !pour 2 phases elastiques

   cas = 2

   ! cas homogene
   if (cas .eq. 1) then
   coeff(1) =  2 ! nb_phases
   coeff(2) =  2 ! nb coeff phase 1
   coeff(3) =  2 ! nb coeff phase 2
   coeff(4) =  0.5_8 ! fv phase 1
   coeff(5) =  100000. ! lambda phase 1
   coeff(6) =  50000. ! Mu phase 1
   coeff(7) =  100000. ! lambda phase 2
   coeff(8) =  50000. ! Mu phase 2

   ! cas heterogene
   else if (cas .eq. 2) then
   coeff(1) =  2 ! nb_phases
   coeff(2) =  2 ! nb coeff phase 1
   coeff(3) =  2 ! nb coeff phase 2
   coeff(4) =  0.1_8 ! fv phase 1
   coeff(5) =  1000000. ! lambda phase 1
   coeff(6) =  500000. ! Mu phase 1
   coeff(7) =  100000. ! lambda phase 2
   coeff(8) =  50000. ! Mu phase 2
   end if
print *,"coeff = ", coeff
! pour verifification - cas homogene (OK) 
!       Tr(def)=0.005
!       l = 100e3
!       m = 50e3
!       sig1 = l*0.005 + 2*m*-0.0025  
!       sig2 = sig1
!       sig3 = l*0.005 + 2*m*0.01
!       sig1 = sig2 = 250
!       sig3 = 1500

! pour verifification - cas heterogene (OK) 
!       Tr(def)=0.005
!       l = 0.1*100e4 + 0.9*100e3
!       m = 0.1*50e4 + 0.9*50e3
!       l = 100e3
!       m = 50e3
!       sig1 = l*0.005 + 2*m*-0.0025  
!       sig2 = sig1
!       sig3 = l*0.005 + 2*m*0.01
!       sig1 = sig2 = 475
!       sig3 = 2850

! INITIALISATION POUR UTILISATION UMAT      
! HYPER IMPORTANT (sinon NTENS n'est pas connu)
!                                 (non utilisees)
!------------------------------------------------
  call initParam_umat()
  KINC = 1
  KSTEP = 0
  DFGRD0 = 0.
  DFGRD1 = 0.
  PREDEF = 0.


!-------------------------------- APPEL INCREMENTAL A LA LOI
  sig =0.
  def0 =0.
  write (*,"(A)"), "# colonne 1    : temps"
  write (*,"(A)"), "# colonne 2-7  : deformation"
  write (*,"(A)"), "# colonne 8-13 : contrainte"
  write (*,"(13(1X,E15.8))"), t, def0,sig
  do i=1,Niter
    t=t+dt
    def=def0+ddef
    call umatVoigt(umatptr,sig,varint,def0,ddef,(/0._8,t-dt/),dt,temp,PREDEF,&
                     nvarint,coeff,ncoeff,DFGRD0, DFGRD1,KSTEP,KINC)
    def0=def
  write (*,"(13(1X,E15.8))"), t, def0,sig
  end do
!--------------------------------

  call cpu_time(t2)
  write (*,"(A,E15.8)"),  "#temps ecoule", t2-t1



print *, "varint = ", varint
print *, "sig = ", sig


end program test_umatvoigt
