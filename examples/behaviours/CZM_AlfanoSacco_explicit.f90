!=======================================================================
!>           FONCTION Cohesive zone for interphase
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
!=================================================================================SUBROUTINE interphase( STRESS, STATEV, DDSDDE, SSE, SPD, SCD,&
SUBROUTINE CZM_AlfanoSacco_explicit( STRESS, STATEV, DDSDDE, SSE, SPD, SCD,&
                       RPL, DDSDDT, DRPLDE, DRPLDT,&
                       STRAN, DSTRAN, TIME, DTIME,&
                       TEMP, DTEMP, PREDEF, DPRED,&
                       CMNAME, NDI, NSHR, NTENS, NSTATV,&
                       PROPS, NPROPS, COORDS,&
                       DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,&
                       NOEL, NPT, LAYER, KSPT, KSTEP, KINC )
  use decomp_2d
  implicit none
!  Arguments de l'interface umat
!--------------------------------------------------------
   CHARACTER*16  CMNAME
   INTEGER*8    NDI, NSHR, NTENS, NSTATV, NPROPS,&
                NOEL, NPT, LAYER, KSPT, KSTEP, KINC
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

   !local variables    
   double precision,dimension(6) :: def
   double precision,dimension(3) :: e1,e2,e3
   double precision :: K1,K2,s1,s2,s01,s02,s2_1,s2_2,ang,sig,tau,sig0,tau0,&
                       GcI,GcII,s2di0,s2di,sigd_el,taud_el,mu_fric,eta,d,&
                       fricCrit,d_lambda,taud,s2_sign
   double precision :: G
   double precision,dimension(2) :: n2
   double precision,parameter   :: pi=3.141592653589793238462643383279502884
   double precision,parameter   :: kk=1.-1e-3 
   double precision,parameter   :: eps_machine=epsilon(1.d0)
   !read the material properties
   sig0 = PROPS(1)
   GcI = PROPS(2)
   tau0 = PROPS(3)
   GcII = PROPS(4)
   eta = PROPS(5)
   mu_fric = PROPS(6)  !friction coefficient
   CELENT = PROPS(7)   !characteristic length
   !read the local axes
   e1 = PROPS(8:10)/(norm2(PROPS(8:10)))
   e2 = PROPS(11:13) - (dot_product(PROPS(11:13),e1) * e1)
   e2 = e2 / norm2(e2)
   e3(1) = e1(2)*e2(3) - e1(3)*e2(2)
   e3(2) = e1(3)*e2(1) - e1(1)*e2(3)
   e3(3) = e1(1)*e2(2) - e1(2)*e2(1)
   !compute other CZM parameters
   s01 = (1.-eta) * (2.*GcI/sig0)
   s02 = (1.-eta) * (2.*GcII/tau0)
   K1 = sig0**2 / (1.-eta) / 2. / GcI
   K2 = tau0**2 / (1.-eta) / 2. / GcII
   !strain intrement 
   def = STRAN + DSTRAN
   !change of coordinates (global --> local)
   call ChangeCoord(def,e1,e2,e3,def)
   !strain --> displacement jump
   s1 = def(1) * CELENT
   s2_1 = def(4)*CELENT*0.5 
   s2_2 = def(5)*CELENT*0.5
   s2 = sqrt( s2_1**2 + s2_2**2 ) !principal sliding displ
   !direction vector of the shear displ.
   if (s2/=0) then
      n2 = (/ s2_1/s2, s2_2/s2 /)
   else
      n2 = STATEV(3:4)!if s2==0, use the direction of the previous step
   end if
   !sign of the shear displacement
   if ( (n2(1)>eps_machine) .OR. ((abs(n2(1))<=eps_machine).AND.(n2(2)>0)) ) then
         !Use machine precision instead of zero to reduce numerical noise.
      s2_sign = 1.
   else
      s2_sign = -1.
   end if
   !update the shear displ. with the sign
   s2 = s2 * s2_sign
   !damage variable
   d = STATEV(1)
   if (d>kk) d = kk !2020.04.24
   !"elastic" trial values
   s2di = STATEV(2) !residual inelastic shear displ.
   if (s1>=0) then
      sigd_el = 0.
   else
      sigd_el = K1 * s1
   end if
   taud_el = K2 * (s2 - s2di)
   sig = (1.-d) * K1 * s1 + d * sigd_el
   !shear stress with friction law
   fricCrit = mu_fric*sigd_el + abs(taud_el)
   if ( fricCrit<=0 ) then !no shear displ.
      tau = (1.-d) * K2 * s2 + d * taud_el
   else !additional shear displ.
      !shear stress with friction law
      d_lambda = fricCrit / K2
      taud = -mu_fric*sigd_el * taud_el/abs(taud_el)
      !total stress
      tau = (1.-d) * K2 * s2 + d * taud
      s2di = s2di + d_lambda * taud_el/abs(taud_el)
   end if
   !take into account the sign of the shear displ.
   tau = tau * s2_sign
   !stress tensor (anti-plane components) : tension
   STRESS(1) = sig
   !stress tensor (anti-plane compoenents) : shear
   STRESS(4) = tau * n2(1) !projection along the two axes
   STRESS(5) = tau * n2(2)
   !stress tensor (in-plane components)
   STRESS(2) = 0 !NOTE:this exacerbates convergence issues when crack grows from interface to bulk
   STRESS(3) = 0
   STRESS(6) = 0
   !change of coordinates: local --> global
   call ChangeCoord2(STRESS,e1,e2,e3,STRESS)
   !update the state variables
!   STATEV(1) = d
   STATEV(2) = s2di
   STATEV(3) = n2(1)
   STATEV(4) = n2(2)

end subroutine CZM_AlfanoSacco_explicit


!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************
!******************************************** Subroutines ******************************************
!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************
!123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!===================================================================================================

subroutine damageUpdate_(s1,s2,s01,s02,eta,d)
  implicit none
  double precision,intent(out)   :: d
  double precision,intent(in)    :: s1,s2,s01,s02,eta
  double precision               :: beta,s1P
  !positive part of normal displ. <s1>+
  if (s1<0) then
     s1P = 0.
  else
     s1P = s1
  end if
  !mixed-mode damage criterion
  beta = sqrt( (s1P/s01)**2 + (s2/s02)**2 ) - 1.
  !damage evolution
  if (beta<=0) then
     return
  else
     d = beta/(1.+beta)/eta
  end if
end subroutine damageUpdate_


subroutine ChangeCoord(defx,e1,e2,e3,defx1)
! transform from global to local coordinates
! ATTENTION: only useful for strain tensor --> (eps4,eps5,eps6)=2*(eps12,eps13,eps23)
!                                              this holds in both local and global coords
!  input:
!              defx - "tensor" in the global coordinates
!         e1,e2,e3 - 3 axes describing the local coordinates
! output:
!             defx1 - "tensor" in the local coordinates
  implicit none
  double precision,dimension(3),intent(in)  :: e1,e2,e3
  double precision,dimension(6),intent(in)  :: defx
  double precision,dimension(6),intent(out) :: defx1
  defx1(1) = (      defx(1)*e1(1) + 0.5*defx(4)*e1(2) + 0.5*defx(5)*e1(3))* e1(1) &
             + (0.5*defx(4)*e1(1) +     defx(2)*e1(2) + 0.5*defx(6)*e1(3))* e1(2) &
             + (0.5*defx(5)*e1(1) + 0.5*defx(6)*e1(2) +     defx(3)*e1(3))* e1(3)
  defx1(2) = (      defx(1)*e2(1) + 0.5*defx(4)*e2(2) + 0.5*defx(5)*e2(3))* e2(1) &
             + (0.5*defx(4)*e2(1) +     defx(2)*e2(2) + 0.5*defx(6)*e2(3))* e2(2) &
             + (0.5*defx(5)*e2(1) + 0.5*defx(6)*e2(2) +     defx(3)*e2(3))* e2(3)
  defx1(3) = (      defx(1)*e3(1) + 0.5*defx(4)*e3(2) + 0.5*defx(5)*e3(3))* e3(1) &
             + (0.5*defx(4)*e3(1) +     defx(2)*e3(2) + 0.5*defx(6)*e3(3))* e3(2) &
             + (0.5*defx(5)*e3(1) + 0.5*defx(6)*e3(2) +     defx(3)*e3(3))* e3(3)
  defx1(4) = ( (    defx(1)*e2(1) + 0.5*defx(4)*e2(2) + 0.5*defx(5)*e2(3))* e1(1) &
             + (0.5*defx(4)*e2(1) +     defx(2)*e2(2) + 0.5*defx(6)*e2(3))* e1(2) &
             + (0.5*defx(5)*e2(1) + 0.5*defx(6)*e2(2) +     defx(3)*e2(3))* e1(3) ) * 2.
  defx1(5) = ( (    defx(1)*e3(1) + 0.5*defx(4)*e3(2) + 0.5*defx(5)*e3(3))* e1(1) &
             + (0.5*defx(4)*e3(1) +     defx(2)*e3(2) + 0.5*defx(6)*e3(3))* e1(2) &
             + (0.5*defx(5)*e3(1) + 0.5*defx(6)*e3(2) +     defx(3)*e3(3))* e1(3) ) * 2.
  defx1(6) = ( (    defx(1)*e3(1) + 0.5*defx(4)*e3(2) + 0.5*defx(5)*e3(3))* e2(1) &
             + (0.5*defx(4)*e3(1) +     defx(2)*e3(2) + 0.5*defx(6)*e3(3))* e2(2) &
             + (0.5*defx(5)*e3(1) + 0.5*defx(6)*e3(2) +     defx(3)*e3(3))* e2(3) ) * 2.
end subroutine ChangeCoord

subroutine ChangeCoord2(sig1,e1,e2,e3,sig)
!transform from local to global coordinates
!ATTENTION: only useful for stress --> [sig4,sig5,sig6]=[sig12,sig13,sig23]
!                                      this holds for both local and global coords
!!  input:
!!             sig1 - "tensor" in the local coordinates
!!         e1,e2,e3 - 3 axes describing the local coordinates
!! output:
!!              sig - "tensor" in the global coordinates
  double precision,dimension(3)::e1,e2,e3
  double precision,dimension(6)::sig,sig1
  sig(1) =  (sig1(1)*e1(1) + sig1(4)*e2(1) + sig1(5)*e3(1))* e1(1) &
          + (sig1(4)*e1(1) + sig1(2)*e2(1) + sig1(6)*e3(1))* e2(1) &
          + (sig1(5)*e1(1) + sig1(6)*e2(1) + sig1(3)*e3(1))* e3(1)
  sig(2) =  (sig1(1)*e1(2) + sig1(4)*e2(2) + sig1(5)*e3(2))* e1(2) &
          + (sig1(4)*e1(2) + sig1(2)*e2(2) + sig1(6)*e3(2))* e2(2) &
          + (sig1(5)*e1(2) + sig1(6)*e2(2) + sig1(3)*e3(2))* e3(2)
  sig(3) =  (sig1(1)*e1(3) + sig1(4)*e2(3) + sig1(5)*e3(3))* e1(3) &
          + (sig1(4)*e1(3) + sig1(2)*e2(3) + sig1(6)*e3(3))* e2(3) &
          + (sig1(5)*e1(3) + sig1(6)*e2(3) + sig1(3)*e3(3))* e3(3)
  sig(4) =  (sig1(1)*e1(2) + sig1(4)*e2(2) + sig1(5)*e3(2))* e1(1) &
          + (sig1(4)*e1(2) + sig1(2)*e2(2) + sig1(6)*e3(2))* e2(1) &
          + (sig1(5)*e1(2) + sig1(6)*e2(2) + sig1(3)*e3(2))* e3(1)
  sig(5) =  (sig1(1)*e1(3) + sig1(4)*e2(3) + sig1(5)*e3(3))* e1(1) &
          + (sig1(4)*e1(3) + sig1(2)*e2(3) + sig1(6)*e3(3))* e2(1) &
          + (sig1(5)*e1(3) + sig1(6)*e2(3) + sig1(3)*e3(3))* e3(1)
  sig(6) =  (sig1(1)*e1(3) + sig1(4)*e2(3) + sig1(5)*e3(3))* e1(2) &
          + (sig1(4)*e1(3) + sig1(2)*e2(3) + sig1(6)*e3(3))* e2(2) &
          + (sig1(5)*e1(3) + sig1(6)*e2(3) + sig1(3)*e3(3))* e3(2)
end subroutine ChangeCoord2
