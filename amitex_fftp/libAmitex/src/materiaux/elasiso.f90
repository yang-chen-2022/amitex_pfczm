!=======================================================================
!>           FONCTION ELASISO
!!-----------------------------------------------------------------------
!! modele elastique lineaire isotrope 
!!
!!-----------------------------------------------------------------------
!!
!! \param[in]  STRAN   (tableau reel(6)) tenseur des deformations
!! \param[in]  DSTRAN  (tableau reel(6)) inccremenet tenseur des deformations
!! \param[in]  PROPS   (tableau reel(2)) coefficients lambda et mu
!!
!! \param[out] STRESS  (tableau reel(6)) tenseur des contraintes
!!
!! Notation CAST3M (Voigt, ordre 11 22 33 12 13 23)
!=======================================================================
SUBROUTINE elasiso( STRESS, STATEV, DDSDDE, SSE, SPD, SCD,&
                       RPL, DDSDDT, DRPLDE, DRPLDT,&
                       STRAN, DSTRAN, TIME, DTIME,&
                       TEMP, DTEMP, PREDEF, DPRED,&
                       CMNAME, NDI, NSHR, NTENS, NSTATV,&
                       PROPS, NPROPS, COORDS,&
                       DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,&
                       NOEL, NPT, LAYER, KSPT, KSTEP, KINC ) &
                       BIND(C,name='elasiso')  ! added for compatibility with nvfortran                  

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

  double precision,dimension(6)                 :: Def
  double precision                              :: Tr_Def,Mu2
  integer                                       :: i

  Def = STRAN+DSTRAN
  Mu2 = 2._8 * PROPS(2)
  Tr_Def = Def(1) + Def(2) + Def(3)
  do i = 1,3
     STRESS(i) = (PROPS(1)) * Tr_Def + Mu2 * Def(i)
  end do
  do i = 4,6
     STRESS(i) = PROPS(2) * Def(i)
  end do

end subroutine elasiso
