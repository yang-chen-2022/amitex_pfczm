!=======================================================================
!>           FUNCTION FOURIER_ISO
!!-----------------------------------------------------------------------
!! 
!!     Isotrope Fourier law with polarization, q = - k. grad(T) 
!!
!!-----------------------------------------------------------------------
!!
!! \param[in]  GRADQD  double(3), temperature gradient
!! \param[in]  DGRADQD double(3), increment of the temperature gradient
!! \param[in]  PROPS   double(4) 
!!                         PROPS(1) = k
!!                         PROPS(2:4) = tau
!!
!! \param[out] FLUXD   double(3), output flux
!!
!!
!=======================================================================
SUBROUTINE Fourier_iso(FLUXD,GRADQD,DGRADQD,TIME,DTIME,TEMP,DTEMP,   &
       PREDEF,DPRED,NVARD,PROPS,NPROPS,KINC) &
       BIND(C,name='fourier_iso')       

  implicit none

!     Arguments de l'interface umat
!--------------------------------------------------------
!
      INTEGER*8    NVARD, NPROPS,KINC
!
      REAL*8       FLUXD(3,NVARD), &
                   GRADQD(3,NVARD),&
                   DGRADQD(3,NVARD),&
                   TIME(2), DTIME,&
                   TEMP, DTEMP, PREDEF(*), DPRED(*),&
                   PROPS(NPROPS)
!--------------------------------------------------------


  FluxD = - PROPS(1) * (GRADQD + DGRADQD)

end subroutine Fourier_iso
