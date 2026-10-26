!=======================================================================
!>           FONCTION ELASISO for Phase-Field
!!-----------------------------------------------------------------------
!!
!!-----------------------------------------------------------------------
!!
!! \param[in]  STRAN   (tableau reel(6)) tenseur des deformations
!! \param[in]  DSTRAN  (tableau reel(6)) inccremenet tenseur des deformations
!! \param[in]  PROPS   (tableau reel(2)) coefficients (lambda, mu)
!!
!! \param[inout]  STATEV (tableau reel(2)) damage and history fields
!!
!! \param[out] STRESS    (tableau reel(6)) tenseur des contraintes
!!
!! Notation CAST3M (Voigt, ordre 11 22 33 12 13 23)
!! NOTE: this subroutine requires another subroutine "rs" for computing the eigen*
!=================================================================================
SUBROUTINE elasisoSW( STRESS, STATEV, DDSDDE, SSE, SPD, SCD,&
                       RPL, DDSDDT, DRPLDE, DRPLDT,&
                       STRAN, DSTRAN, TIME, DTIME,&
                       TEMP, DTEMP, PREDEF, DPRED,&
                       CMNAME, NDI, NSHR, NTENS, NSTATV,&
                       PROPS, NPROPS, COORDS,&
                       DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,&
                       NOEL, NPT, LAYER, KSPT, KSTEP, KINC )

  use decomp_2d

  implicit none

!     Arguments de l'interface umat
!--------------------------------------------------------
      CHARACTER*16  CMNAME
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

  double precision,dimension(6)    :: Def
  double precision,dimension(3)    :: Def_sw
  double precision                 :: Tr_Def
  !get the current strain
  Def = STRAN+DSTRAN

  !swelling strain
  Def_sw(1:3) = PROPS(3) * (PREDEF(1)+DPRED(1)) / 3
  Def(1:3) = Def(1:3) - Def_sw
   
  !Hooke's law
  Tr_Def = Def(1) + Def(2) + Def(3)
  STRESS(1:3) = PROPS(1)*Tr_Def + 2.*PROPS(2)*Def(1:3)
  STRESS(4:6) = PROPS(2)*Def(4:6)
 
end subroutine elasisoSW

