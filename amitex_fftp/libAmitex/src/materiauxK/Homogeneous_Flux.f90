!=======================================================================
!>           FUNCTION HOMOGENEOUS_FLUX
!!-----------------------------------------------------------------------
!! 
!!     Output a Homogeneous Flux given from coefficients (PROPS(1:3)) 
!!
!!-----------------------------------------------------------------------
!!
!! \param[in]  PROPS   double(3), components of the output flux
!!
!! \param[out] FluxD   double(3), output flux
!!
!! Warning : Only available for NVARD=1
!=======================================================================
SUBROUTINE Homogeneous_Flux(FLUXD,GRADQD,DGRADQD,TIME,DTIME,TEMP,DTEMP, &
       PREDEF,DPRED,NVARD,PROPS,NPROPS,KINC) &
       BIND(C,name='homogeneous_flux')
       
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

  FluxD(:,1) = PROPS(1:3)

end subroutine Homogeneous_Flux
