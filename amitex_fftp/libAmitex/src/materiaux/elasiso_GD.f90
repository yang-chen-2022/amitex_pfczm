
!=======================================================================
!> MATTOVEC(DEFOR,STRAN) Transforme une matrice sym 3x3 en vecteur à 6 composantes 
!=======================================================================
SUBROUTINE  MATTOVEC(DEFOR,STRAN)

  implicit none

  REAL*8 STRAN(6),DEFOR(3,3)
  
  STRAN(1)=DEFOR(1,1)
  STRAN(2)=DEFOR(2,2)
  STRAN(3)=DEFOR(3,3)
  STRAN(4)=DEFOR(2,1)
  STRAN(5)=DEFOR(3,1)
  STRAN(6)=DEFOR(3,2)
END SUBROUTINE MATTOVEC
!=======================================================================


!=======================================================================
!> Calcul du VRAI determinant (VDET) de la matrice A donnee en lignes   
!=======================================================================
subroutine VDET(A,DET)
  
  implicit none
  
  double precision,DIMENSION(3,3),intent(in) :: A
  double precision,intent(out) :: DET
  double precision :: A11,A21,A31,A12,A22,A32,A33,A13,A23,W1,W2,W3

  A11 = A(1,1)
  A22 = A(2,2)
  A33 = A(3,3)
  A12 = A(1,2)
  A13 = A(1,3)
  A23 = A(2,3)
  A21 = A(2,1)
  A31 = A(3,1)
  A32 = A(3,2)
  W1 = A22*A33 - A32*A23
  W2 = A31*A23 - A21*A33
  W3 = A21*A32 - A31*A22
  DET = A11*W1 + A12*W2 + A13*W3
  
end subroutine VDET

!=======================================================================
!>           FONCTION ELASISO_GD
!!-----------------------------------------------------------------------
!! modele elastique lineaire isotrope 
!!
!!-----------------------------------------------------------------------
!!
!! \param[in]  DFGRD1  (tableau reel(9)) tenseur 1 + grad(u) 
!! \param[in]  PROPS   (tableau reel(2)) coefficients lambda et mu
!! \param[out] STRESS  (tableau reel(6)) tenseur des contraintes de Cauchy
!!
!! Notation CAST3M pour Sigma ordre 11 22 33 12 13 23
!! Notation pour DFGRD1 ordre 11 21 31 12 22 32 13 23 33
!=======================================================================
SUBROUTINE elasiso_GD( STRESS, STATEV, DDSDDE, SSE, SPD, SCD,&
                       RPL, DDSDDT, DRPLDE, DRPLDT,&
                       STRAN, DSTRAN, TIME, DTIME,&
                       TEMP, DTEMP, PREDEF, DPRED,&
                       CMNAME, NDI, NSHR, NTENS, NSTATV,&
                       PROPS, NPROPS, COORDS,&
                       DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,&
                       NOEL, NPT, LAYER, KSPT, KSTEP, KINC ) &
                       BIND(C,name='elasiso_gd')! added for compatibility with nvfortran


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

  double precision,dimension(3,3)  :: STRESSMAT,EDEF
  !! todo : changer en logical ou integer
  double precision, dimension(3,3), parameter :: ID=reshape((/1,0,0,0,1,0,0,0,1/), (/3,3/))
  double precision                            ::    DETDFGRD,MU2,TREDEF

  !---- Loi elastique isotrope en grandes deformations : 
  !----     e = 1/2*(transpose(F)*F - ID)  
  EDEF   = 0.5_8 * (MATMUL(TRANSPOSE(DFGRD1),DFGRD1)- ID)
  !----     MU2 = 2 * mu
  MU2    = 2._8 * PROPS(2)
  !----     TREDEF = Trace(e)
  TREDEF = EDEF(1,1) + EDEF(2,2) + EDEF(3,3)
  call VDET(DFGRD1,DETDFGRD)
  !----     Sig = 1/det(F)  F *(lambda*Tr(e)ID +2 mu e)*transpose(F)
  STRESSMAT = 1._8/DETDFGRD * MATMUL(DFGRD1,MATMUL(PROPS(1)*TREDEF*ID+&
       MU2*EDEF,TRANSPOSE(DFGRD1)))
  CALL MATTOVEC(STRESSMAT,STRESS)
  !---- On met Sig sous forme d'un vecteur de taille 6 (matrice symetrique)
  RETURN
end subroutine elasiso_GD
!=======================================================================
