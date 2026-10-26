!=======================================================================
!>           FONCTION ELASANISO
!!-----------------------------------------------------------------------
!! modele elastique lineaire anisotrope 
!!
!!-----------------------------------------------------------------------
!!
!! \param[in]  STRAN   (tableau reel(6)) tenseur des deformations
!! \param[in]  DSTRAN  (tableau reel(6)) inccremenet tenseur des deformations
!! \param[in]  PROPS   (tableau reel(9+3+3+6)) 
!!                          9 coefficients rangés par ligne
!!                          1 2 3 0 0 0
!!                          2 4 5 0 0 0
!!                          3 5 6 0 0 0
!!                          0 0 0 7 0 0
!!                          0 0 0 0 8 0
!!                          0 0 0 0 0 9
!!
!!                          vecteurs e1 et e2 de la base locale
!!
!!                          6 composantes de deformation libre 
!!                             (notation de Voigt, ordre CAST3M)
!!
!! \param[out] STRESS     (tableau reel(6)) tenseur des contraintes
!!
!! Notation CAST3M (Voigt, ordre 11 22 33 12 13 23)
!=======================================================================
SUBROUTINE elasaniso_eigs( STRESS, STATEV, DDSDDE, SSE, SPD, SCD,&
                       RPL, DDSDDT, DRPLDE, DRPLDT,&
                       STRAN, DSTRAN, TIME, DTIME,&
                       TEMP, DTEMP, PREDEF, DPRED,&
                       CMNAME, NDI, NSHR, NTENS, NSTATV,&
                       PROPS, NPROPS, COORDS,&
                       DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,&
                       NOEL, NPT, LAYER, KSPT, KSTEP, KINC ) &
                       BIND(C,name='elasaniso_eigs') ! added for compatibility with nvfortran


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
  double precision,dimension(6)      :: sig1,def1, defL
  double precision,dimension(3)      :: e1,e2,e3,v1,v2,v3

  Def = STRAN+DSTRAN

! Calcul de e1 (e1 normalisé) 
  e1 = PROPS(10:12)/(norm2(PROPS(10:12)))

! Calcul de e2 (deuxième vecteur local) e2 = e2 - (e2.e1)e1, normalisé
  e2 = PROPS(13:15) - (dot_product(PROPS(13:15),e1) * e1) 
  e2 = e2 / norm2(e2)

! Calcul de e3 (troisieme vecteur local)
  e3(1) = e1(2)*e2(3) - e1(3)*e2(2)
  e3(2) = e1(3)*e2(1) - e1(1)*e2(3) 
  e3(3) = e1(1)*e2(2) - e1(2)*e2(1) 

! Deformation Libre
  defL = PROPS(16:21)

! Passage de la deformation dans la base locale (def1)
  def1(1) = (           Def(1)*e1(1) + 0.5_8*Def(4)*e1(2) + 0.5_8*Def(5)*e1(3))* e1(1) &
          + (0.5_8*Def(4)*e1(1) +            Def(2)*e1(2) + 0.5_8*Def(6)*e1(3))* e1(2) & 
          + (0.5_8*Def(5)*e1(1) + 0.5_8*Def(6)*e1(2) +            Def(3)*e1(3))* e1(3)

  def1(2) = (           Def(1)*e2(1) + 0.5_8*Def(4)*e2(2) + 0.5_8*Def(5)*e2(3))* e2(1) &
          + (0.5_8*Def(4)*e2(1) +            Def(2)*e2(2) + 0.5_8*Def(6)*e2(3))* e2(2) & 
          + (0.5_8*Def(5)*e2(1) + 0.5_8*Def(6)*e2(2) +            Def(3)*e2(3))* e2(3)

  def1(3) = (           Def(1)*e3(1) + 0.5_8*Def(4)*e3(2) + 0.5_8*Def(5)*e3(3))* e3(1) &
          + (0.5_8*Def(4)*e3(1) +            Def(2)*e3(2) + 0.5_8*Def(6)*e3(3))* e3(2) & 
          + (0.5_8*Def(5)*e3(1) + 0.5_8*Def(6)*e3(2) +            Def(3)*e3(3))* e3(3)

  def1(4) = (           Def(1)*e1(1) + 0.5_8*Def(4)*e1(2) + 0.5_8*Def(5)*e1(3))* e2(1) &
          + (0.5_8*Def(4)*e1(1) +            Def(2)*e1(2) + 0.5_8*Def(6)*e1(3))* e2(2) & 
          + (0.5_8*Def(5)*e1(1) + 0.5_8*Def(6)*e1(2) +            Def(3)*e1(3))* e2(3) &
          + (           Def(1)*e2(1) + 0.5_8*Def(4)*e2(2) + 0.5_8*Def(5)*e2(3))* e1(1) &
          + (0.5_8*Def(4)*e2(1) +            Def(2)*e2(2) + 0.5_8*Def(6)*e2(3))* e1(2) & 
          + (0.5_8*Def(5)*e2(1) + 0.5_8*Def(6)*e2(2) +            Def(3)*e2(3))* e1(3)

  def1(5) = (           Def(1)*e1(1) + 0.5_8*Def(4)*e1(2) + 0.5_8*Def(5)*e1(3))* e3(1) &
          + (0.5_8*Def(4)*e1(1) +            Def(2)*e1(2) + 0.5_8*Def(6)*e1(3))* e3(2) & 
          + (0.5_8*Def(5)*e1(1) + 0.5_8*Def(6)*e1(2) +            Def(3)*e1(3))* e3(3) &
          + (           Def(1)*e3(1) + 0.5_8*Def(4)*e3(2) + 0.5_8*Def(5)*e3(3))* e1(1) &
          + (0.5_8*Def(4)*e3(1) +            Def(2)*e3(2) + 0.5_8*Def(6)*e3(3))* e1(2) & 
          + (0.5_8*Def(5)*e3(1) + 0.5_8*Def(6)*e3(2) +            Def(3)*e3(3))* e1(3)

  def1(6) = (           Def(1)*e2(1) + 0.5_8*Def(4)*e2(2) + 0.5_8*Def(5)*e2(3))* e3(1) &
          + (0.5_8*Def(4)*e2(1) +            Def(2)*e2(2) + 0.5_8*Def(6)*e2(3))* e3(2) & 
          + (0.5_8*Def(5)*e2(1) + 0.5_8*Def(6)*e2(2) +            Def(3)*e2(3))* e3(3) &
          + (           Def(1)*e3(1) + 0.5_8*Def(4)*e3(2) + 0.5_8*Def(5)*e3(3))* e2(1) &
          + (0.5_8*Def(4)*e3(1) +            Def(2)*e3(2) + 0.5_8*Def(6)*e3(3))* e2(2) & 
          + (0.5_8*Def(5)*e3(1) + 0.5_8*Def(6)*e3(2) +            Def(3)*e3(3))* e2(3)


! Calcul de la contrainte dans la base locale
  def1 = def1 - defL
  sig1(1) = PROPS(1)*def1(1) + PROPS(2)*def1(2) + PROPS(3)*def1(3) 
  sig1(2) = PROPS(2)*def1(1) + PROPS(4)*def1(2) + PROPS(5)*def1(3) 
  sig1(3) = PROPS(3)*def1(1) + PROPS(5)*def1(2) + PROPS(6)*def1(3) 
  sig1(4) = PROPS(7)*def1(4) 
  sig1(5) = PROPS(8)*def1(5) 
  sig1(6) = PROPS(9)*def1(6) 

! Retour de la contrainte dans la base macrosopique
  
  v1(1) = e1(1)
  v1(2) = e2(1)
  v1(3) = e3(1)

  v2(1) = e1(2)
  v2(2) = e2(2)
  v2(3) = e3(2)

  v3(1) = e1(3)
  v3(2) = e2(3)
  v3(3) = e3(3)

  STRESS(1) =  (sig1(1)*v1(1) + sig1(4)*v1(2) + sig1(5)*v1(3))* v1(1) &
          + (sig1(4)*v1(1) + sig1(2)*v1(2) + sig1(6)*v1(3))* v1(2) & 
          + (sig1(5)*v1(1) + sig1(6)*v1(2) + sig1(3)*v1(3))* v1(3)

  STRESS(2) =  (sig1(1)*v2(1) + sig1(4)*v2(2) + sig1(5)*v2(3))* v2(1) &
          + (sig1(4)*v2(1) + sig1(2)*v2(2) + sig1(6)*v2(3))* v2(2) & 
          + (sig1(5)*v2(1) + sig1(6)*v2(2) + sig1(3)*v2(3))* v2(3)

  STRESS(3) =  (sig1(1)*v3(1) + sig1(4)*v3(2) + sig1(5)*v3(3))* v3(1) &
          + (sig1(4)*v3(1) + sig1(2)*v3(2) + sig1(6)*v3(3))* v3(2) & 
          + (sig1(5)*v3(1) + sig1(6)*v3(2) + sig1(3)*v3(3))* v3(3)

  STRESS(4) =  (sig1(1)*v1(1) + sig1(4)*v1(2) + sig1(5)*v1(3))* v2(1) &
          + (sig1(4)*v1(1) + sig1(2)*v1(2) + sig1(6)*v1(3))* v2(2) & 
          + (sig1(5)*v1(1) + sig1(6)*v1(2) + sig1(3)*v1(3))* v2(3) &
          + (sig1(1)*v2(1) + sig1(4)*v2(2) + sig1(5)*v2(3))* v1(1) &
          + (sig1(4)*v2(1) + sig1(2)*v2(2) + sig1(6)*v2(3))* v1(2) & 
          + (sig1(5)*v2(1) + sig1(6)*v2(2) + sig1(3)*v2(3))* v1(3)

  STRESS(5) =  (sig1(1)*v1(1) + sig1(4)*v1(2) + sig1(5)*v1(3))* v3(1) &
          + (sig1(4)*v1(1) + sig1(2)*v1(2) + sig1(6)*v1(3))* v3(2) & 
          + (sig1(5)*v1(1) + sig1(6)*v1(2) + sig1(3)*v1(3))* v3(3) &
          + (sig1(1)*v3(1) + sig1(4)*v3(2) + sig1(5)*v3(3))* v1(1) &
          + (sig1(4)*v3(1) + sig1(2)*v3(2) + sig1(6)*v3(3))* v1(2) & 
          + (sig1(5)*v3(1) + sig1(6)*v3(2) + sig1(3)*v3(3))* v1(3)
 
  STRESS(6) =  (sig1(1)*v2(1) + sig1(4)*v2(2) + sig1(5)*v2(3))* v3(1) &
          + (sig1(4)*v2(1) + sig1(2)*v2(2) + sig1(6)*v2(3))* v3(2) & 
          + (sig1(5)*v2(1) + sig1(6)*v2(2) + sig1(3)*v2(3))* v3(3) & 
          + (sig1(1)*v3(1) + sig1(4)*v3(2) + sig1(5)*v3(3))* v2(1) &
          + (sig1(4)*v3(1) + sig1(2)*v3(2) + sig1(6)*v3(3))* v2(2) & 
          + (sig1(5)*v3(1) + sig1(6)*v3(2) + sig1(3)*v3(3))* v2(3)

  STRESS(4)=STRESS(4)/2._8
  STRESS(5)=STRESS(5)/2._8
  STRESS(6)=STRESS(6)/2._8

end subroutine elasaniso_eigs
