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
SUBROUTINE elasisoPF( STRESS, STATEV, DDSDDE, SSE, SPD, SCD,&
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

  double precision,dimension(6)    :: Def, DefP, DefN
  double precision,dimension(3)    :: Def_sw
  double precision, dimension(3)   :: lbda  !eigenvalues in ascending order
  double precision, dimension(3,3) :: vec   !eigenvectors: vec(:,1), vec(:,2), vec(:,3)
  double precision,dimension(6,6)  :: PP0, PN0 !projection tensor, DefP=PP*Def
  double precision                 :: Tr_Def, Tr_Def_P,Tr_Def_N, psiP,gd
  integer                          :: i
  double precision, parameter      :: kk = 1e-6 !YC:2022.03.19
  integer ( kind = 4 )             :: ierr

  !get the current strain
  Def = STRAN+DSTRAN

  !split the strain into postive and negative parts
  call splitStrain( Def, DefP, DefN, PP0,PN0 )
     !-------------------------------------------------------------------------
     ! print *,'in umat,def=',def
     ! print *,'in umat,old method,defP=',defP
     ! !eigenvalues and eigenvectors of the strain tensor
     ! call rs( 3, reshape( (/def(1), def(4)*0.5, def(5)*0.5,&
     !                        def(4)*0.5, def(2), def(6)*0.5,&
     !                        def(5)*0.5, def(6)*0.5, def(3)/) ,(/3,3/)),&
     !          lbda, .true., vec, ierr )
     ! call projTensor(lbda,vec,PP,PN) !cf. [Cervera Tesei, IJSS 2018]
     ! call tensorMultipl4ordX2ord(PP,def,defP)
     ! defP(4:6) = defP(4:6)*2.
     ! defN = def - defP
     ! print *,'in umat,new method,defP=',defP
     !NOTE-1:This way to split the strain tensor gives exaxtly the same result 
     !       as what we have proposed in [Chen et al.,CMAME 2019, Appendix A]
     !NOTE-2:Both methods have the same issue of convergence when the unloading
     !       reaches zero (tested for shearXY load on One Element)
     !       Another test without strain spliting did not encouter the
     !       We also tested an unloading that avoided zero (-1e-3)
     !       --> so this issue should be related to the strain spliting
     !       Alternative solution: instead of zero, we set the unloading "node"
     !                             as a SMALL VALUE with the sign changed
     !--------------------------------------------------------------------------

  !compute the stress
  Tr_Def = Def(1) + Def(2) + Def(3)
  Tr_Def_P = ( Tr_Def + abs(Tr_Def) ) / 2.
  Tr_Def_N = ( Tr_Def - abs(Tr_Def) ) / 2.
  gd = (1.-STATEV(1))**2 + kk
  STRESS(1:3) = gd * ( PROPS(1)*Tr_Def_P + 2.*PROPS(2)*DefP(1:3) ) &
                   + ( PROPS(1)*Tr_Def_N + 2.*PROPS(2)*DefN(1:3) )
  STRESS(4:6) = gd * ( PROPS(2)*DefP(4:6) ) &
                   + ( PROPS(2)*DefN(4:6) )
  !compute the history field
  psiP = 0.5*PROPS(1) * ( Tr_Def_P )**2 &
           + PROPS(2) * ( sum(DefP(1:3)**2) + sum(DefP(4:6)**2)/2. )
  !update the history field
  STATEV(2) = psiP
 
end subroutine elasisoPF

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

subroutine splitStrain( def0, defP, defN, PP,PN )
    !split strain tensor into positive / negative parts according to their eigenvalues
    !---------------------------------------------------------------------------------
    ! param[in]
    !      > def0 : strain tensor (6x1 array)
    !
    ! param[out]
    !      > defP : positive part of strain (6x1 array)
    !      > defN : negative part of strain (6x1 array)
    !      >   PP : (optional) positive projection tensor, defP=PP:def
    !      >   PN : (optional) negative projection tensor, defN=PN:def
    !
    implicit none
!    double precision, dimension(3,3) :: eps0, epsP, epsN
    double precision, dimension(6)   :: defP, defN, def0
    double precision, dimension(6,6) :: PP,PN
    double precision, dimension(3)   :: lbda  !eigenvalues in ascending order
    double precision, dimension(3,3) :: vec   !eigenvectors: vec(:,1), vec(:,2), vec(:,3)
    integer ( kind = 4 ) :: ierr
    integer(kind=8)   :: i
    !eigenvalues and eigenvectors of the strain tensor
    call rs( 3, reshape( (/def0(1), def0(4)*0.5, def0(5)*0.5,&
                           def0(4)*0.5, def0(2), def0(6)*0.5,&
                           def0(5)*0.5, def0(6)*0.5, def0(3)/) ,(/3,3/)),&
             lbda, .true., vec, ierr )
    !construct the projection tensor PP
    do i=1,3
       if (lbda(i)>=0) then
          lbda(i) = 1
       else
          lbda(i) = 0
       end if
    end do
    PP(1,1) = lbda(1)*vec(1,1)*vec(1,1) + lbda(2)*vec(1,2)*vec(1,2) + lbda(3)*vec(1,3)*vec(1,3)
    PP(2,2) = lbda(1)*vec(2,1)*vec(2,1) + lbda(2)*vec(2,2)*vec(2,2) + lbda(3)*vec(2,3)*vec(2,3)
    PP(3,3) = lbda(1)*vec(3,1)*vec(3,1) + lbda(2)*vec(3,2)*vec(3,2) + lbda(3)*vec(3,3)*vec(3,3)
    PP(4,1) = lbda(1)*vec(1,1)*vec(2,1) + lbda(2)*vec(1,2)*vec(2,2) + lbda(3)*vec(1,3)*vec(2,3)
    PP(4,2) = lbda(1)*vec(1,1)*vec(2,1) + lbda(2)*vec(1,2)*vec(2,2) + lbda(3)*vec(1,3)*vec(2,3)
    PP(5,1) = lbda(1)*vec(1,1)*vec(3,1) + lbda(2)*vec(1,2)*vec(3,2) + lbda(3)*vec(1,3)*vec(3,3)
    PP(5,3) = lbda(1)*vec(1,1)*vec(3,1) + lbda(2)*vec(1,2)*vec(3,2) + lbda(3)*vec(1,3)*vec(3,3)
    PP(6,2) = lbda(1)*vec(2,1)*vec(3,1) + lbda(2)*vec(2,2)*vec(3,2) + lbda(3)*vec(2,3)*vec(3,3)
    PP(6,3) = lbda(1)*vec(2,1)*vec(3,1) + lbda(2)*vec(2,2)*vec(3,2) + lbda(3)*vec(2,3)*vec(3,3)
    PP(1,4) = 0.5*PP(4,2)
    PP(1,5) = 0.5*PP(5,3)
    PP(2,4) = 0.5*PP(4,1)
    PP(2,6) = 0.5*PP(6,3)
    PP(3,5) = 0.5*PP(5,1)
    PP(3,6) = 0.5*PP(6,2)
    PP(4,4) = 0.5*( PP(1,1) + PP(2,2) )
    PP(5,5) = 0.5*( PP(1,1) + PP(3,3) )
    PP(6,6) = 0.5*( PP(2,2) + PP(3,3) )
    PP(4,5) = PP(2,6)
    PP(4,6) = PP(1,5)
    PP(5,4) = PP(3,6)
    PP(5,6) = PP(1,4)
    PP(6,4) = PP(3,5)
    PP(6,5) = PP(2,4)
    PP(1,2) = 0.
    PP(1,3) = 0.
    PP(1,6) = 0.
    PP(2,1) = 0.
    PP(2,3) = 0.
    PP(2,5) = 0.
    PP(3,1) = 0.
    PP(3,2) = 0.
    PP(3,4) = 0.
    PP(4,3) = 0.
    PP(5,2) = 0.
    PP(6,1) = 0.
    !split the strain by multiplying the projection tensor
    do i=1,6
       defP(i) = sum(PP(i,:)*def0)
    end do
    PN = reshape( (/1.,0.,0.,0.,0.,0.,&
                    0.,1.,0.,0.,0.,0.,&
                    0.,0.,1.,0.,0.,0.,&
                    0.,0.,0.,1.,0.,0.,&
                    0.,0.,0.,0.,1.,0.,&
                    0.,0.,0.,0.,0.,1./), (/6,6/) )-PP
    do i=1,6
       defN(i) = sum(PN(i,:)*def0)
    end do
end subroutine splitStrain

