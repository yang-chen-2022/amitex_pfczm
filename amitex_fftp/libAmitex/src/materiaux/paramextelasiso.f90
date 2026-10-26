!=======================================================================
!>           FONCTION PARAMEXTELASISO
!!-----------------------------------------------------------------------
!! modele thermoelastique lineaire isotrope - ou la 'temperature' 
!! est lue dans PREDEF,DPRED
!!-----------------------------------------------------------------------
!!
!! \param[in]  STRAN   (tableau reel(6)) tenseur des deformations
!! \param[in]  DSTRAN  (tableau reel(6)) increment du tenseur des deformations
!! \param[in]  PREDEF  (tableau reel(1)) 'temperature' en debut de pas
!! \param[in]  DPRED   (tableau reel(1)) increment de 'temperature'
!! \param[in]  PROPS   (tableau reel(4)) coefficients lambda, mu,
!!                                       alpha et T0
!!             lambda et mu sont les coefficients de Lame
!!             alpha est le coefficient de dilatation thermique
!!             T0 est la temperature de reference
!! \param[out] STRESS     (tableau reel(6)) tenseur des contraintes
!!
!! Notation CAST3M (Voigt, ordre 11 22 33 12 13 23)
!=======================================================================
SUBROUTINE paramextelasiso( STRESS, STATEV, DDSDDE, SSE, SPD, SCD,&
                       RPL, DDSDDT, DRPLDE, DRPLDT,&
                       STRAN, DSTRAN, TIME, DTIME,&
                       TEMP, DTEMP, PREDEF, DPRED,&
                       CMNAME, NDI, NSHR, NTENS, NSTATV,&
                       PROPS, NPROPS, COORDS,&
                       DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,&
                       NOEL, NPT, LAYER, KSPT, KSTEP, KINC )&
                       BIND(C,name='paramextelasiso') ! added for compatibility with nvfortran

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
  !! On remplace Def par Def + Coeff(3) * (TEMPF - COEFF(4)) * Id
  !! Il faut donc rajouter 3 * Coeff(3) * (TEMPF - COEFF(4)) a la trace
  Tr_Def = Def(1) + Def(2) + Def(3) - 3.0_8 * PROPS(3) * (PREDEF(1) + DPRED(1) - PROPS(4))
  do i = 1,3
     !! De meme on rajoute Coeff(3) * (TEMPF - COEFF(4)) pour les coefficients diagonaux
     STRESS(i) = Lambda0 * Tr_Def + Mu2 * (Def(i) - PROPS(3) * (PREDEF(1) + DPRED(1) - PROPS(4)))
  end do
  do i = 4,6
     !! Pour les coefficients extra-diagonaux (2 \mu \epsilon) on n'a pas besoin
     !! de prendre en compte la dilatation thermique et la notation de Voigt permet
     !! de remplacer 2 \mu par \mu
     STRESS(i) = Mu0 * Def(i)
  end do

end subroutine paramextelasiso
