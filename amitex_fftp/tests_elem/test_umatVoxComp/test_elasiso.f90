program  test_elasiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> CHARGEMENT DES DIFFERENTS MODULES
   use ISO_C_BINDING
!  use MPI             
!  use decomp_2d
!  use decomp_2d_fft

!  use io_amitex_mod
!  use green_mod
!  use param_algo_mod
!  use loading_mod
  use material_mod
!  use error_mod
!  use resolution_mod
!  use amitex_mod
!  use field_mod

!!==============================================================================
!!						      		    DECLARATIONS
!!==============================================================================
  implicit none
!-----------------------------------------------------------------------------
!         PENSER A AJUSTER LA DECLARATION varInt(nvarint) et coeff(ncoeff)
!-----------------------------------------------------------------------------
!----------------------------------------------- PARAMETRES DE UMAT
  double precision, dimension(6)   :: sig,def,ddef,def0,ddsddt,drplde
  double precision                 :: temp,dtemp,SSE, SPD, SCD,RPL,drpldt,pnewdt,celent
  double precision, dimension(6,6) :: ddsdde
  double precision, dimension(3,3) :: DFGRD0, DFGRD1, DROT
  double precision, dimension(3)   :: coords
  double precision, dimension(2)   ::  coeff    !A AJUSTER A LA LOI
  double precision, dimension(1)   ::  varInt   !A AJUSTER A LA LOI
  double precision, dimension(1)   :: predef,dpred
  character(len=16)                :: law
  integer(kind=8)                  :: kinc,ndi,ntens,ncoeff,nvarint,nshr,NOEL, NPT, LAYER, KSPT, KSTEP
  double precision :: t,dt
!---------------------------------------------- AUTRES PARAMETRES

  character(len=200,kind=C_CHAR)   :: lawname
  character(len=200,kind=C_CHAR)   :: libname

  real     :: t1,t2
  integer  :: i,Niter

  type(c_ptr) :: umatptr

!==========================================================================
  law="               1"   !ne sert plus dans le cas d'un appel par umat_call

  t1=0.
  t2=0.
  call cpu_time(t1)


!                                 Choix de la loi
!------------------------------------------------
!

 libname='/home/gelebart/amitex_fftp/libAmitex/src/materiaux/libUmatAmitex.so'
 lawname = 'elasiso'
 print *, "libname, lawname", trim(libname), trim(lawname)

 umatptr = umat_load(trim(libname)//C_NULL_CHAR,trim(lawname)//C_NULL_CHAR)

! Definition du chargement en deformation imposee
!------------------------------------------------
  t=0._8           ! temps initial
  Niter = 10      ! nombre d'increments de temps
  dt=100._8/Niter  ! increment de temps 

  temp = 293.15_8  ! temperature initiale
  dtemp = 0._8     ! increment de temperature

  def0 =0._8       ! deformation initiale

  ddef(1)=-0.0025_8/Niter ! increments de deformation (6 composantes, notation de Voigt)
  ddef(2)=-0.0025_8/Niter
  ddef(3)=0.01_8/Niter
  ddef(4)=0.0_8
  ddef(5)=0.0_8
  ddef(6)=0.0_8
  
!           Initialisation des variables internes
!------------------------------------------------
   nvarint=0
   varint=dble(0)

!                                 Coeff materiaux
!------------------------------------------------
   ncoeff = 2
   coeff(1) =   100000_8
   coeff(2)=    50000_8

! pour verifification (OK) 
!       Tr(def)=0.005
!       sig(1) = 100e3*0.005 + 2*50e3*-0.0025 = 100e3(0.0025) = 250 
!       sig(2) = 250
!       sig(3) = 100e3(0.015)=1500

!      autres init. pour adherence au format umat
!                                 (non utilisees)
!------------------------------------------------
  sig =0._8
  def0 =0._8
  ndi = 2
  ntens=6
  ddsdde=0._8
  sse = 0._8
  spd=0._8
  scd=0._8
  rpl=0._8
  ddsddt=0._8
  drplde=0._8
  drpldt=0._8
  predef=293.15_8
  dpred=0._8
  nshr=0
  coords=0._8
  drot=0._8
  drot(1,1)=1._8
  drot(2,2)=1._8
  drot(3,3)=1._8
  pnewdt=0._8
  celent=0._8
  NOEL= 0
  NPT=0 
  LAYER=0 
  KSPT=0 
  KSTEP=0
  kinc=1

!-------------------------------- APPEL INCREMENTAL A LA LOI
  write (*,"(A)"), "# colonne 1    : temps"
  write (*,"(A)"), "# colonne 2-7  : deformation"
  write (*,"(A)"), "# colonne 8-13 : contrainte"
  write (*,"(13(1X,E15.8))"), t, def0,sig
  do i=1,Niter
  t=t+dt
  def=def0+ddef

  call umat_call(umatptr,sig,varint,ddsdde,sse,spd,scd, &
                RPL, DDSDDT, DRPLDE, DRPLDT,&
                def0,ddef,(/0._8,t-dt/),dt,&
                temp,dtemp,PREDEF,DPRED,&
                law,ndi,nshr,ntens,nvarint,&
                coeff,ncoeff,coords,&
                DROT,pnewdt,celent,DFGRD0, DFGRD1,&
                NOEL, NPT, LAYER, KSPT, KSTEP,kinc)

  def0=def
  write (*,"(13(1X,E15.8))"), t, def0,sig
  end do
!--------------------------------

  call cpu_time(t2)
  write (*,"(A,E15.8)"),  "#temps ecoule", t2-t1

end program test_elasiso
