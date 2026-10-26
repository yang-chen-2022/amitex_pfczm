program  test_umatvoigtCC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> CHARGEMENT DES MODULES
  use ISO_C_BINDING
  use material_mod
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
  double precision, dimension(2*nphases+2*27) ::  coeff      !umatVoigt avec umatBCCHPP
  double precision, dimension(1*nphases + 6*nphases + 2* 49)   ::  varInt   !umatVoigt avec umatBCCHPP elasticite
  double precision, dimension(1)   :: predef
  integer(kind=8)                  :: kinc, ncoeff,nvarint, KSTEP
  double precision :: t,dt
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
  t=0.             ! temps initial
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
   nvarint=1*nphases + 6*nphases + 2*49
   varint = 0.
   varint(1) = 49	! nb varint de la phase 1
   varint(2) = 49        ! nb varint de la phase 2

!                                 Coeff materiaux (SOTERIA Fer Pur BCC 12s  Perform60)
!------------------------------------------------
   ncoeff = 2*nphases + 2*27 !pour 2 phases elastiques
   cas = 2
   ! cas homogene
   if (cas .eq. 1) then
   coeff(1) =  2 ! nb_phases
   coeff(2) =  27 ! nb coeff phase 1
   coeff(3) =  27 ! nb coeff phase 2
   coeff(4) =  0.5_8 ! fv phase 1

   coeff(5) =  236.412E3
   coeff(6) =  0.35
   coeff(7) =  275.2E3
   coeff(8) =  112.4E3
   coeff(9) =  87.56E3
   coeff(10) =  363.
   coeff(11) =  0.
   coeff(12) =  0.25
   coeff(13) =  10e-6
   coeff(14) =  1e11
   coeff(15) =  2.481e-7
   coeff(16) =  100
   coeff(17) =  3.
   coeff(18) =  2.e-6
   coeff(19) =  0.84
   coeff(20) =  6e5
   coeff(21) =  6e5
   coeff(22) =  1e-6
   coeff(23) =  50
   coeff(24) =  0
   coeff(25) =  1.
   coeff(26) =  0.
   coeff(27) =  5.e-4 !gpmoyen??
   coeff(28) =  temp 
   coeff(29) =  0
   coeff(30) =  0
   coeff(31) =  0
   
   coeff(32:32+26) = coeff(5:31)


   ! cas heterogene
   else if (cas .eq. 2) then
   coeff(1) =  2 ! nb_phases
   coeff(2) =  27 ! nb coeff phase 1
   coeff(3) =  27 ! nb coeff phase 2
   coeff(4) =  0.5_8 ! fv phase 1

   !Phase 1
   coeff(5) =  236.412E3
   coeff(6) =  0.35
   coeff(7) =  275.2E3
   coeff(8) =  112.4E3
   coeff(9) =  87.56E3
   coeff(10) =  363.
   coeff(11) =  0.
   coeff(12) =  0.25
   coeff(13) =  10e-6
   coeff(14) =  1e11
   coeff(15) =  2.481e-7
   coeff(16) =  100
   coeff(17) =  3.
   coeff(18) =  2.e-6
   coeff(19) =  0.84
   coeff(20) =  6e5
   coeff(21) =  6e5
   coeff(22) =  1e-6
   coeff(23) =  50
   coeff(24) =  0
   coeff(25) =  1.
   coeff(26) =  0.
   coeff(27) =  5.e-4 !gpmoyen??
   coeff(28) =  temp 
   coeff(29) =  0
   coeff(30) =  0
   coeff(31) =  0
   
   !Phase 2
   coeff(32:32+26) = coeff(5:31)
   coeff(56) =  1. ! changement de rotation
   coeff(57) =  2.
   coeff(58) =  3.

   end if
   !print *,"coeff = ", coeff

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
!---
  do i=1,Niter
    t=t+dt
    def=def0+ddef
    call umatVoigt(umatptr,sig,varint,def0,ddef,(/0._8,t-dt/),dt,temp,PREDEF,&
                     nvarint,coeff,ncoeff,DFGRD0, DFGRD1,KSTEP,KINC)
    def0=def
  write (*,"(13(1X,E15.8))"), t, def0,sig
  end do
!---
  ddef = -ddef !on inverse le sens de sollicitation
  do i=1,Niter/2
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

!==============================================================================
!  VERIFICATION MODELE DE VOIGT
!  compatibilite

print *, "#-----------"
print *, "#Sig = ", sig
print *, "#-----------"
print * ,'#moy(sigi) = ', coeff(4)*varint(nphases+1:nphases+6)+(1-coeff(4))*varint(nphases+6+49+1:nphases+6+49+6)



  call MPI_FINALIZE(ierror)
end program test_umatvoigtCC
