!=======================================================================
!>           FUNCTION FOURIER_ISO_POLARIZATION
!!-----------------------------------------------------------------------
!! 
!!     Isotrope Fourier law with polarization, q = - k. grad(T) + tau 
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
!! Warning : Only available for NVARD=1
!=======================================================================
SUBROUTINE Fourier_iso_polarization(FLUXD,GRADQD,DGRADQD,TIME,DTIME,TEMP,DTEMP,&
       PREDEF,DPRED,NVARD,PROPS,NPROPS,KINC) &
       BIND(C,name='fourier_iso_polarization')
       
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
  FluxD(:,1) = FluxD(:,1) + PROPS(2:4)


end subroutine Fourier_iso_polarization
