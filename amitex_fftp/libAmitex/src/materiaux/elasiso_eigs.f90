!=======================================================================
!>           FONCTION ELASISO_EIGS
!!-----------------------------------------------------------------------
!! modele elastique lineaire isotrope avec deformation libre (variable interne)
!!
!!-----------------------------------------------------------------------
!!
!! \param[in]  STRAN   (tableau reel(6)) tenseur des deformations
!! \param[in]  DSTRAN  (tableau reel(6)) increment du tenseur des deformations
!! \param[in]  PROPS   (tableau reel(2)) coefficients lambda, mu,
!!             lambda et mu sont les coefficients de Lame
!! \param[in]  STATEV  (tableau reel(6)) tenseur des deformations libres
!!
!! \param[out] STRESS  (tableau reel(6)) tenseur des contraintes
!!
!! Notation CAST3M (Voigt, ordre 11 22 33 12 13 23)
!=======================================================================
SUBROUTINE elasiso_eigs( STRESS, STATEV, DDSDDE, SSE, SPD, SCD,&
                       RPL, DDSDDT, DRPLDE, DRPLDT,&
                       STRAN, DSTRAN, TIME, DTIME,&
                       TEMP, DTEMP, PREDEF, DPRED,&
                       CMNAME, NDI, NSHR, NTENS, NSTATV,&
                       PROPS, NPROPS, COORDS,&
                       DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,&
                       NOEL, NPT, LAYER, KSPT, KSTEP, KINC ) &
                       BIND(C,name='elasiso_eigs')! added for compatibility with nvfortran


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

  double precision,dimension(6)      :: Def
  double precision                   :: Tr_Def,Mu2,Mu0,Lambda0
  integer                            :: i

  Def = STRAN+DSTRAN
  Lambda0 = PROPS(1)
  Mu0 = PROPS(2)
  Mu2 = 2._8 * Mu0
  Tr_Def = Def(1) + Def(2) + Def(3) - STATEV(1) - STATEV(2) - STATEV(3)
  do i = 1,3
     STRESS(i) = Lambda0 * Tr_Def + Mu2 * (Def(i) - STATEV(i))
  end do
  do i = 4,6
     !! Pour les coefficients extra-diagonaux (2 \mu \epsilon) on n'a pas besoin
     !! de prendre en compte la dilatation thermique et la notation de Voigt permet
     !! de remplacer 2 \mu par \mu
     STRESS(i) = Mu0 * (Def(i) -STATEV(i))
  end do

end subroutine elasiso_eigs

