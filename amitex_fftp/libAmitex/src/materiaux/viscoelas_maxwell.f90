!=======================================================================
!           FONCTION VISCOELAS_MAXWELL
!-----------------------------------------------------------------------
!
! modele viscoelastique lineaire maxwell generalise pour k et mu
!
!-----------------------------------------------------------------------
!input:
!      NTENS     : (entier) dimension des tenseurs (contraintes et deformations)
!      NPROPS    : (entier)
!      NSTATV    : (entier) nombre de variables internes
!      DTIME     : (reel) pas de temps
!      DSTRAN    : (tableau reel(6)) tenseurs des increments de deformations
!      PROPS     : (tableau reel(NPROPS)) proprietes du materiaux
!
!
!output:
!      STRESS    : (tableau reel(6)) tenseur des contraintes
!      STATEV    : (tableau reel(6)) variables internes
!
! Notation CAST3M (Voigt, ordre 11 22 33 12 13 23)
!=======================================================================

SUBROUTINE viscoelas_maxwell( STRESS, STATEV, DDSDDE, SSE, SPD, SCD,&
                       RPL, DDSDDT, DRPLDE, DRPLDT,&
                       STRAN, DSTRAN, TIME, DTIME,&
                       TEMP, DTEMP, PREDEF, DPRED,&
                       CMNAME, NDI, NSHR, NTENS, NSTATV,&
                       PROPS, NPROPS, COORDS,&
                       DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,&
                       NOEL, NPT, LAYER, KSPT, KSTEP, KINC ) &
                       BIND(C,name='viscoelas_maxwell') ! added for compatibility with nvfortran

  implicit none

!     Arguments de l'interface umat
!--------------------------------------------------------
      !CHARACTER*16  CMNAME ! modified for compatibility with nvfortran
      CHARACTER*1  CMNAME
!
      INTEGER*8    NDI, NSHR, NTENS, NSTATV, NPROPS,&
                   NOEL, NPT, LAYER, KSPT, KSTEP, KINC
!
      REAL*8       STRESS(NTENS), STATEV(NSTATV),&
                   DDSDDE(NTENS,NTENS),&
                   SSE, SPD, SCD,&
                   RPL, DDSDDT(NTENS), DRPLDE(NTENS), DRPLDT,&
                   STRAN(NTENS), DSTRAN(NTENS),&
                   TIME(2), DTIME,&
                   TEMP, DTEMP, PREDEF(*), DPRED(*),&
                   PROPS(NPROPS),&
                   COORDS(3),&
                   DROT(3,3),&
                   PNEWDT,&
                   CELENT,&
                   DFGRD0(3,3), DFGRD1(3,3)
!--------------------------------------------------------

  integer                       :: NCHAIN, M, i
  double precision              :: depsiv, xk, xmu, Tmu, CONS, CONS2, CONS3, FLUA
  double precision, dimension(6):: de



  depsiv=DSTRAN(1)+DSTRAN(2)+DSTRAN(3)

  do i=1,NTENS
    de(i)=DSTRAN(i)-depsiv/3._8
    if (i.GT.3) then
      de(i)=DSTRAN(i)/2._8
    end if
    STRESS(i) = 0._8
  end do

!  xk0 = PROPS(1)
!  xmu0 = PROPS(2)

! nombre total de chaines
         NCHAIN = int(PROPS(3))

! recuperation de k, mu et tau:
!
  DO M=0,NCHAIN-1
! DO 10 M=0,0
    if (M.GT.0) then
      xk=PROPS(3*M+3)
      xmu=PROPS(3*M+4)
      Tmu=PROPS(3*M+5)
      FLUA=EXP(-DTIME/Tmu)
      CONS=2._8*xmu*Tmu/DTIME*(1._8-FLUA)
      CONS2=xk*Tmu/DTIME*(1._8-FLUA)

    else
      xk=PROPS(4)
      xmu=PROPS(5)
      CONS=2._8*xmu
      CONS2=xk
      FLUA = 1._8
    end if
!
    cons3 = cons2*depsiv
!           
! calcul des contraintes dans chaque chaine
! STATEV(M+1)=sigk1;  STATEV(M+2)=sigk2; STATEV(M+3)=sigk3;
! STATEV(M+4)=sigmu1;  STATEV(M+5)=sigmu2; STATEV(M+6)=sigmu3;
! STATEV(M+7)=sigmu4;  STATEV(M+8)=sigmu5; STATEV(M+9)=sigmu6;
    do i=1,NTENS
      STATEV(i+(NTENS+3)*M+3) = STATEV(i+(NTENS+3)*M+3)*FLUA + cons*de(i)
      if (i.LT.4)then
        STATEV(i+(NTENS+3)*M)=STATEV(i+(NTENS+3)*M)*FLUA+cons3
        STRESS(i) = STRESS(i)+STATEV(i+(NTENS+3)*M)
      endif
      STRESS(i)=STRESS(i)+STATEV(i+(NTENS+3)*M+3)
    end do
  end do
  !continue 10
end subroutine viscoelas_maxwell
