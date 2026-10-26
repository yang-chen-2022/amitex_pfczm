program  test_CC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> CHARGEMENT DES DIFFERENTS MODULES
  use ISO_C_BINDING
  use material_mod

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
  double precision, dimension(27)   ::  coeff    !A AJUSTER A LA LOI
  double precision, dimension(49)   ::  varInt   !A AJUSTER A LA LOI
  double precision, dimension(1)   :: predef,dpred
  character(len=16)                :: law
  integer(kind=8)                  :: kinc,ndi,ntens,ncoeff,nvarint,nshr,NOEL, NPT, LAYER, KSPT, KSTEP
  double precision :: t,dt
!--------------------------------------------- Pointeur de fonctions
  type(c_ptr) :: umatptr  
!---------------------------------------------- AUTRES PARAMETRES

  character(len=200,kind=C_CHAR)   :: lawname
  character(len=200,kind=C_CHAR)   :: libname

  real     :: t1,t2
  integer  :: i,Niter


!==========================================================================

  t1=0.
  t2=0.
  call cpu_time(t1)


!                                 Choix de la loi
!------------------------------------------------
!

 libname='/home/gelebart/amitex_fftp/cas_tests/comportements/polyxCC/comportement_umat/libUmatAmitex.so'
 lawname = 'umatbcchpp'
 print *, "libname, lawname", trim(libname), trim(lawname)

 umatptr = umat_load(trim(libname)//C_NULL_CHAR,trim(lawname)//C_NULL_CHAR)

! Definition du chargement en deformation imposee
!------------------------------------------------
  t=0._8           ! temps initial
  Niter = 100      ! nombre d'increments de temps
  dt=100._8/Niter  ! increment de temps 

  temp = 293.15_8  ! temperature initiale
  dtemp = 0._8     ! increment de temperature

  def0 =0._8       ! deformation initiale

  ddef(1)=-0.005_8/Niter ! increments de deformation (6 composantes, notation de Voigt)
  ddef(2)=-0.005_8/Niter
  ddef(3)=0.01_8/Niter
  ddef(4)=0.0_8
  ddef(5)=0.0_8
  ddef(6)=0.0_8
  
!           Initialisation des variables internes
!------------------------------------------------
   nvarint=49
   varint=dble(0)

!                                 Coeff materiaux (SOTERIA Fer Pur BCC 12s  Perform60)
!------------------------------------------------
   ncoeff = 27
   coeff(1) =  236.412E3
   coeff(2) =  0.35
   coeff(3) =  275.2E3
   coeff(4) =  112.4E3
   coeff(5) =  87.56E3
   coeff(6) =  363.
   coeff(7) =  0.
   coeff(8) =  0.25
   coeff(9) =  10e-6
   coeff(10) =  1e11
   coeff(11) =  2.481e-7
   coeff(12) =  100
   coeff(13) =  3.
   coeff(14) =  2.e-6
   coeff(15) =  0.84
   coeff(16) =  6e5
   coeff(17) =  6e5
   coeff(18) =  1e-6
   coeff(19) =  50
   coeff(20) =  0
   coeff(21) =  1.
   coeff(22) =  0.
   coeff(23) =  5.e-4  !gpmoyen??
   coeff(24) =  temp
   coeff(25) =  0
   coeff(26) =  0
   coeff(27) =  0





! INITIALISATION POUR UTILISATION UMAT      
! HYPER IMPORTANT (sinon NTENS n'est pas connu)
!                                 (non utilisees)
!------------------------------------------------
!  call initParam_umat()
      Law="0"
      KINC = 1
      NSHR = 0
      NOEL = 0
      NPT = 0
      layer = 0
      kspt = 0
      !KSTEP = 0
      NTENS=6 !pour umat, il s'agit du tenseur des contraintes de Cauchy
      DTEMP = 0.
      DPRED = 0.
      ddsdde = 0.
      sse = 0.
      spd = 0.
      scd = 0.
      rpl = 0.
      ddsddt = 0.
      drplde = 0.
      drpldt = 0.
      COORDS = 0.
      drot = 0.
      drot(1,1)=1.
      drot(2,2)=1.
      drot(3,3)=1.
      pnewdt = 0.
      celent = 0.
      DFGRD0 = 0.
      DFGRD1 = 0.


!-------------------------------- APPEL INCREMENTAL A LA LOI
  sig =0.
  def0 =0.
  write (*,"(A)"), "# colonne 1    : temps"
  write (*,"(A)"), "# colonne 2-7  : deformation"
  write (*,"(A)"), "# colonne 8-13 : contrainte"
  write (*,"(13(1X,E15.8))"), t, def0,sig
!---
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
!---
  ddef = -ddef !on inverse le sens de sollicitation
  do i=1,Niter/2
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

end program test_CC
