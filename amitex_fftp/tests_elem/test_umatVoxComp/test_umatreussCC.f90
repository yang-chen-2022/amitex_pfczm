program  test_umatreussCC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> CHARGEMENT DES MODULES
  use ISO_C_BINDING
  use material_mod
  use param_algo_mod
  use mpi 

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
  double precision, dimension(31*nphases)   ::  coeff      !umatReuss pour umatBCCHPP
  double precision, dimension(62*nphases)   ::  varInt   !umatReuss pour umatBCCHPP
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
  integer  :: i,Niter, cas,ierror

  call MPI_INIT(ierror) ! necessaire si l'on souhaite avoir les message d'erreur amitex_abort

!==========================================================================

  t1=0.
  t2=0.
  call cpu_time(t1)


!                                 Choix de la loi
!------------------------------------------------
!

 libname='/home/gelebart/amitex_fftp/cas_tests/comportements/polyxCC/comportement_umat/libUmatAmitex.so'
 lawname = 'umatbcchpp'
 print *, "libname, lawname", trim(libname),'   ', trim(lawname)

 umatptr(1) = umat_load(trim(libname)//C_NULL_CHAR,trim(lawname)//C_NULL_CHAR)
 umatptr(2) = umat_load(trim(libname)//C_NULL_CHAR,trim(lawname)//C_NULL_CHAR)


! Definition du chargement en deformation imposee
!------------------------------------------------
  t=0.           ! temps initial
  Niter = 100      ! nombre d'increments de temps
  dt=100._8/Niter  ! increment de temps 

  temp = 293.15_8  ! temperature initiale
  dtemp = 0.     ! increment de temperature

  def0 =0.       ! deformation initiale

  ddef(1)=-0.005_8/Niter ! increments de deformation (6 composantes, notation de Voigt)
  ddef(2)=-0.005_8/Niter
  ddef(3)=0.01_8/Niter
  ddef(4)=0.0
  ddef(5)=0.0
  ddef(6)=0.0
  
!           Initialisation des variables internes
!------------------------------------------------
   nvarint=1*nphases + 12*nphases + 49*nphases
   varint = 0.
   varint(1) = 49	! nb varint de la phase 1
   varint(2) = 49        ! nb varint de la phase 2

!                                 Coeff materiaux
!------------------------------------------------
   ncoeff = 2*nphases +nphases*(27+2)

   cas = 2

   ! cas homogene
   if (cas .eq. 1) then
   coeff(1) =  2 ! nb_phases
   coeff(2) =  27 ! nb coeff phase 1
   coeff(3) =  27 ! nb coeff phase 2
   coeff(4) =  0.1_8    ! fv phase 1

   coeff(5) =   2.0431e5 ! lambdaeq phase 1
   coeff(6) =   0.8756e5  ! Mueq phase 1
   coeff(7) =  236.412E3
   coeff(8) =  0.35
   coeff(9) =  275.2E3
   coeff(10) =  112.4E3
   coeff(11) =  87.56E3
   coeff(12) =  363.
   coeff(13) =  0.
   coeff(14) =  0.25
   coeff(15) =  10e-6
   coeff(16) =  1e11
   coeff(17) =  2.481e-7
   coeff(18) =  100
   coeff(19) =  3.
   coeff(20) =  2.e-6
   coeff(21) =  0.84
   coeff(22) =  6e5
   coeff(23) =  6e5
   coeff(24) =  1e-6
   coeff(25) =  50
   coeff(26) =  0
   coeff(27) =  1.
   coeff(28) =  0.
   coeff(29) =  5.e-4 !gpmoyen??
   coeff(30) =  temp 
   coeff(31) =  0
   coeff(32) =  0
   coeff(33) =  0

   coeff(34) =   2.0431e5 ! lambdaeq phase 1
   coeff(35) =   0.8756e5  ! Mueq phase 1
   coeff(36:36+26) = coeff(7:33)

   ! cas heterogene
   else if (cas .eq. 2) then
   coeff(1) =  2 ! nb_phases
   coeff(2) =  27 ! nb coeff phase 1
   coeff(3) =  27 ! nb coeff phase 2
   coeff(4) =  0.01_8    ! fv phase 1

   coeff(5) =   2.0431e5 ! lambdaeq phase 1
   coeff(6) =   0.8756e5  ! Mueq phase 1
   coeff(7) =  236.412E3
   coeff(8) =  0.35
   coeff(9) =  275.2E3
   coeff(10) =  112.4E3
   coeff(11) =  87.56E3
   coeff(12) =  363.
   coeff(13) =  0.
   coeff(14) =  0.25
   coeff(15) =  10e-6
   coeff(16) =  1e11
   coeff(17) =  2.481e-7
   coeff(18) =  100
   coeff(19) =  3.
   coeff(20) =  2.e-6
   coeff(21) =  0.84
   coeff(22) =  6e5
   coeff(23) =  6e5
   coeff(24) =  1e-6
   coeff(25) =  50
   coeff(26) =  0
   coeff(27) =  1.
   coeff(28) =  0.
   coeff(29) =  5.e-4 !gpmoyen??
   coeff(30) =  temp 
   coeff(31) =  0
   coeff(32) =  0
   coeff(33) =  0

   coeff(34) =   2.0431e5  ! lambdaeq phase 1
   coeff(35) =   0.8756e5  ! Mueq phase 1
   coeff(36:36+26) = coeff(7:33)
   coeff(60) =  1.	! changement de rotation
   coeff(61) =  2.
   coeff(62) =  3.
  

   end if
!print *,"coeff = ", coeff

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
!  algo_param%acc_CV_reuss = .false.
  algo_param%acc_CV_reuss = .true.
  algo_param%Init_reuss = "default"
!  algo_param%Init_reuss = "linear"
!  algo_param%Init_reuss = "proportional"
  algo_param%Jacobian_type_reuss="elastic"

!-------------------------------- APPEL INCREMENTAL A LA LOI
  sig = 0.
  def0 = 0.
  dt_tab = 0

  write (*,"(A)"), "# colonne 1    : temps"
  write (*,"(A)"), "# colonne 2-7  : deformation"
  write (*,"(A)"), "# colonne 8-13 : contrainte"
  write (*,"(13(1X,E15.8))"), t, def0,sig
!Niter=10
  ddef(4)=0.000001_8/Niter
  ddef(5)=0.000001_8/Niter
  ddef(6)=0.000001_8/Niter

!---
  do i=1,Niter
    t=t+dt
    def=def0+ddef
    dt_tab(1) = dt

    call umatReuss(umatptr,sig,varint,def0,ddef,(/0._8,t-dt/),dt_tab,temp,PREDEF,&
                     nvarint,coeff,ncoeff,DFGRD0, DFGRD1,KSTEP,KINC)

    def0=def
    dt_tab(2)=dt_tab(1)
    write (*,"(13(1X,E15.8))"), t, def0,sig
  end do
!---
!Niter = 2
!  ddef(4)=1e-13
!  ddef(5)=1e-13
!  ddef(6)=0
!--
  ddef = -ddef !on inverse le sens de sollicitation
  do i=1,Niter/2
    t=t+dt
    def=def0+ddef
    dt_tab(1) = dt

    call umatReuss(umatptr,sig,varint,def0,ddef,(/0._8,t-dt/),dt_tab,temp,PREDEF,&
                     nvarint,coeff,ncoeff,DFGRD0, DFGRD1,KSTEP,KINC)

    def0=def
    dt_tab(2)=dt_tab(1)
    write (*,"(13(1X,E15.8))"), t, def0,sig
  end do
!--------------------------------

  call cpu_time(t2)
  write (*,"(A,E15.8)"),  "#temps ecoule", t2-t1




!==============================================================================
!  VERIFICATION MODELE DE REUSS
!  compatibilite
print *, "#-----------"
print *, "#Def_impose = ", def0
print *, "#Defmoy = ", coeff(4)*varint(3:8) + (1-coeff(4))*varint(2+12+49+1:2+12+49+6)
print *, "#-----------"
print *, "#Def1 = ", varint(3:8) 
print *, "#-----------"
print *, "#Def2 = ", varint(2+12+49+1:2+12+49+6)
print *, "#-----------"
print *, "#Sigmoy = ", sig
print *, "#-----------"

  call MPI_FINALIZE(ierror)
end program test_umatreussCC
