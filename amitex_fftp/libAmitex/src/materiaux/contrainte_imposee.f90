!=======================================================================
!>           FONCTION CONTRAINTE_IMPOSEE
!!-----------------------------------------------------------------------
!! return a constant stress tensor with  6 components given in PROPS
!!
!!-----------------------------------------------------------------------
!!
!! \param[in]  PROPS   (tableau reel(6)) components of the stress to be applied
!!
!! \param[out] STRESS  (tableau reel(6)) stress tensor STRESS=PROPS
!!
!! Notation CAST3M (Voigt, ordre 11 22 33 12 13 23)
!=======================================================================
SUBROUTINE contrainte_imposee( STRESS, STATEV, DDSDDE, SSE, SPD, SCD,&
                       RPL, DDSDDT, DRPLDE, DRPLDT,&
                       STRAN, DSTRAN, TIME, DTIME,&
                       TEMP, DTEMP, PREDEF, DPRED,&
                       CMNAME, NDI, NSHR, NTENS, NSTATV,&
                       PROPS, NPROPS, COORDS,&
                       DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,&
                       NOEL, NPT, LAYER, KSPT, KSTEP, KINC ) &
                       BIND(C,name='contrainte_imposee') ! added for compatibility with nvfortran


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

  integer         :: i


  do i = 1,6
     STRESS(i) = PROPS(i) 
  end do

end subroutine contrainte_imposee
