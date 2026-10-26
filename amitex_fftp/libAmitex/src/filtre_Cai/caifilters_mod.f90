module caifilters_mod

  use iso_c_binding
  use ISO_FORTRAN_ENV
  use mpi
  use decomp_2d, only : mytype
  use field_mod, only : fft_start,fft_end
  use error_mod, only : amitex_abort

  private

  !> Fonctions publiques
  public :: field_filter_caiF, wk, wtildek
  
  interface

     real(c_double) function wk(a,k) bind(c)
       use iso_c_binding
       real(c_double), intent(in) :: a
       real(c_double), intent(in) :: k
     end function wk

     real(c_double) function wtildek(a,k) bind(c)
       use iso_c_binding
       real(c_double), intent(in) :: a
       real(c_double), intent(in) :: k
     end function wtildek

  end interface


contains

!==============================================================================
!       SUBROUTINE FIELD_FILTER_CAIF
!------------------------------------------------------------------------------
!> Applique un Filtre a un champ dans l'espace de Fourier 
!!
!!              field_filterF(fieldF,FILTER,ncomp,nvar,N,d)
!!
!!              FieldF = Filtre * fieldF
!!
!! FILTER : 1   : Filtre "Moyenne" 8 pts (avec changement de support noeuds->centres)
!!          -1  : Filtre "Moyenne" 8 pts (avec changement de support centres->noeuds)
!!          2   : Filtre "Moyenne" 8 pts applique 2 fois (retour support initial)
!!
!!          100 : Filtre de Gauss isotrope, param(1) = sigma 
!!                   g=Aexp(-X^2 / (2.sigma^2) (A = facteur de normalisation (integrale = 1)
!!
!!          200 : Filtre de Cai (WK)
!!          201 : Filtre de Cai (WKtilde)
!!
!!   \param[inout] fieldF      champ entree/sortie, (nxP/2,nyP,nzP,ncomp,nvar)    (profil indifferent)
!!   \param[in]  ncomp         nombre de composantes du champ 
!!   \param[in]  nvar          nombre de champs (1 en mecanique)
!!
!!   \param[in]  N(3)          dimension de la cellule dans la direction dir
!!   \param[in]  d(3)          taille des voxels dans la direction dir
!!
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (fieldF et D_fieldF)
!!
!! REMARQUE : les filtres 1, -1, 2 ne dependant pas de la dimension d des voxels  
!!
!==============================================================================
subroutine field_filter_caiF(fieldF,FILTER,ncomp,nvar,N,d,param)

  implicit none
  
  !Entrees/Sorties :
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),ncomp,nvar),&
                intent(inout)          :: fieldF 
  integer,intent(in)                   :: ncomp,nvar,FILTER
  integer,dimension(3),intent(in)      :: N
  real(mytype),dimension(3),intent(in) :: d
  real(mytype),dimension(:),intent(in),optional :: param
  
  ! Autre variables
  integer                     :: i,j,k,l,m ! variables de boucles
#ifdef DOUBLE_PREC
  real(mytype),parameter      :: PI = 4._mytype*DATAN(1._mytype)
#else
  real(mytype),parameter      :: PI = 4._mytype*ATAN(1._mytype)
#endif
  real(mytype),dimension(3)   :: F          ! Frequence
  complex(mytype)             :: imP = cmplx(0,1,kind=mytype) ! le nombre complexe i
  
  integer,dimension(1:N(1)/2 + 1) :: idx    !> Utiles pour definition des frequences (no_filter)
  integer,dimension(1:N(2))       :: idy
  integer,dimension(1:N(3))       :: idz
  integer                         :: nx_2, ny_2, nz_2
  real(mytype)                    :: DF1,DF2,DF3
  real(mytype)                    :: fact, FiltreF
  complex(mytype),dimension(ncomp,nvar) :: fieldF0


  !--- Definition de variables pour l'application de Filtre Gaussien et autres
  fact = 0.  
  DF1 = 2._mytype * PI / (float(N(1))*d(1))
  DF2 = 2._mytype * PI / (float(N(2))*d(2))
  DF3 = 2._mytype * PI / (float(N(3))*d(3))
  
  ! indices pour ramener les frequences nulles en (1,1,1)
  ! attention aux cas pairs et impairs
  nx_2 = N(1) / 2
  idx = (/(nx_2 + i,i=1,nx_2+1) /)

  ny_2 = N(2) / 2
  if (mod(N(2),2) .eq. 0) then
     idy = (/(ny_2 + i,i=1,ny_2),(i,i=1,ny_2)/)
  else
     idy = (/(ny_2 + i,i=1,ny_2+1),(i,i=1,ny_2)/)
  endif

  nz_2 = N(3) / 2
  if (mod(N(3),2) .eq. 0) then
     idz = (/(nz_2 + i,i=1,nz_2),(i,i=1,nz_2)/)
  else
     idz = (/(nz_2 + i,i=1,nz_2+1),(i,i=1,nz_2)/)
  endif

  !--- Filtre 8 voisins applique UNE fois (Centres => Noeuds)
  if (FILTER == -1) then 
  do m=1,nvar
  do l=1,ncomp
  do k=fft_start(3),fft_end(3)
  do j=fft_start(2),fft_end(2)
  do i=fft_start(1),fft_end(1)
    F(1) = PI*(i-1)/N(1)   ! attention ici F= F/2.
    F(2) = PI*(j-1)/N(2)
    F(3) = PI*(k-1)/N(3)
    fieldF(i,j,k,l,m) = fieldF(i,j,k,l,m) * &
                        (cos(F(1)) * cos(F(2)) * cos(F(3))) * &! attention ici F= F/2.
                        exp(-imP*(F(1)+F(2)+F(3))) 
  end do
  end do
  end do
  end do
  end do
  
  !--- Filtre 8 voisins applique UNE fois (Noeuds => Centres)
  elseif (FILTER == 1) then 
  do m=1,nvar
  do l=1,ncomp
  do k=fft_start(3),fft_end(3)
  do j=fft_start(2),fft_end(2)
  do i=fft_start(1),fft_end(1)
    F(1) = PI*(i-1)/N(1)   ! attention ici F= F/2.
    F(2) = PI*(j-1)/N(2)
    F(3) = PI*(k-1)/N(3)
    fieldF(i,j,k,l,m) = fieldF(i,j,k,l,m) * &
                        (cos(F(1)) * cos(F(2)) * cos(F(3))) * & ! attention ici F= F/2.
                        exp(imP*(F(1)+F(2)+F(3))) 
  end do
  end do
  end do
  end do
  end do
  
  !--- Filtre 8 voisins applique DEUX fois (support inchange : C => N => C ou N => C => N)
  elseif (FILTER == 2) then
  do m=1,nvar
  do l=1,ncomp
  do k=fft_start(3),fft_end(3)
  do j=fft_start(2),fft_end(2)
  do i=fft_start(1),fft_end(1)
    F(1) = PI*(i-1)/N(1)   ! attention ici F= F/2.
    F(2) = PI*(j-1)/N(2)
    F(3) = PI*(k-1)/N(3)
    fieldF(i,j,k,l,m) = fieldF(i,j,k,l,m) * &
                         (cos(F(1)) * cos(F(2)) * cos(F(3))) ** 2 ! attention ici F= F/2.
  end do
  end do
  end do
  end do
  end do
  
  !--- Filtre de GAUSS
  !    
  !    g = Aexp(-alpha.X^2)   avec A = (alpha/pi)^3/2 (normalisation)
  !
  !      -> gF = exp(-freq^2/(4.alpha))
  ! 
  !      alpha = 1/(2.sigma^2)  -> 1/(4.alpha) = sigma^2 / 2
  !
  elseif (FILTER == 100) then
  
    fact = param(1)**2 / 2._mytype
    do m=1,nvar
    do l=1,ncomp
    do k=fft_start(3),fft_end(3)
    do j=fft_start(2),fft_end(2)
    do i=fft_start(1),fft_end(1)
      F(1) = (-nx_2 - 1 + idx(i)) * DF1
      F(2) = (-ny_2 - 1 + idy(j)) * DF2
      F(3) = (-nz_2 - 1 + idz(k)) * DF3
      FiltreF = exp(-sum(F*F) * fact)
      fieldF(i,j,k,l,m) = fieldF(i,j,k,l,m) * FiltreF
    end do
    end do
    end do
    end do
    end do
    
  !--- Filtre de CAI wk
  !    
  elseif (FILTER == 200) then
    ! store field for nul Frequency
    if(fft_start(1) == 1 .and. fft_start(2) == 1 .and. fft_start(3)==1) then
    do m=1,nvar
    do l=1,ncomp
      FieldF0(l,m)=fieldF(1,1,1,l,m)
    end do
    end do
    end if

    do m=1,nvar
    do l=1,ncomp
    do k=fft_start(3),fft_end(3)
    do j=fft_start(2),fft_end(2)
    do i=fft_start(1),fft_end(1)
      F(1) = (-nx_2 - 1 + idx(i)) * DF1
      F(2) = (-ny_2 - 1 + idy(j)) * DF2
      F(3) = (-nz_2 - 1 + idz(k)) * DF3
      FiltreF =  wk(param(1),sqrt(sum(F*F)))
      fieldF(i,j,k,l,m) = fieldF(i,j,k,l,m) * FiltreF
    end do
    end do
    end do
    end do
    end do
    
    ! apply field for nul Frequency
    if(fft_start(1) == 1 .and. fft_start(2) == 1 .and. fft_start(3)==1) then
    do m=1,nvar
    do l=1,ncomp
      fieldF(1,1,1,l,m) = FieldF0(l,m)
    end do
    end do
    end if
    
  !--- Filtre de CAI wtildek
  !    
  elseif (FILTER == 201) then
  
    do m=1,nvar
    do l=1,ncomp
    do k=fft_start(3),fft_end(3)
    do j=fft_start(2),fft_end(2)
    do i=fft_start(1),fft_end(1)
      F(1) = (-nx_2 - 1 + idx(i)) * DF1
      F(2) = (-ny_2 - 1 + idy(j)) * DF2
      F(3) = (-nz_2 - 1 + idz(k)) * DF3
      FiltreF = wtildek(param(1),sqrt(sum(F*F)))
      fieldF(i,j,k,l,m) = fieldF(i,j,k,l,m) * FiltreF
    end do
    end do
    end do
    end do
    end do
  
  else
     call amitex_abort("field_filterF : FILTER different from -1,1,2,100,200,201 - not yet implemented " ,2,0)  
  end if

end subroutine field_filter_caiF

end module caifilters_mod
