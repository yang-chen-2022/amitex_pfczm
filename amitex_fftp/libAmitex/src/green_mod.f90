!123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!===================================================================================================
!
!       MODULE GREEN_MOD
!> Module gerant les calculs lies au tenseur de Green
!!
!!
!!  Subroutines
!! - initGrid :      initialisation d'une variable de type GRID_DIM
!! - print_grid :    ecriture d'une variable de type GRID_DIM
!! - initFreq :      initialisation du tableau de frequences
!! - apply_green :   application du tenseur de Green
!! - apply_green_GD :   application du tenseur de Green en grandes transformations
!! - eval_criteres :  calcul des criteres (equilibre, compatiblite, moyenne) en Mecanique
!! - eval_criteres_GD : calcul des criteres (equilibre, compatiblite, moyenne) en Mecanique
!!                    Grandes Transformations
!! - eval_criteresD : calcul des criteres (equilibre, compatiblite, moyenne) en Diffusion
!!
!! \todo Reflechir a l'utilisation d'apply_green pour evaluer le critere
!===================================================================================================
module green_mod

  use ISO_FORTRAN_ENV

#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype, real_type
  use decomp_2d_mpi,       only : nrank
#else
  use decomp_2d,      only : mytype, nrank, real_type
#endif  

  use mpi
  use amitex_mod
  use param_algo_mod
  use error_mod
  use field_mod      , only : fft_start, fft_end, field_sum_squareF, grid, FREQ_2
  use amitex_mpi_mod, only: amitex_mpi_comm
#ifdef OPENMP
  use omp_lib 
#endif

  private

  !> Pointer "publiques" vers composantes de SIMU_AMITEX
  public :: times_g 

  !> Variables "publiques" utilisees lors de l'initialisation
  !public :: grid0

  !> Variables publiques a transformer en pointer (ou variables 'provisoires')
  ! RAS

  !> Types publiques (pour definition de SIMU_AMITEX)
  public :: TIMES_GREEN

  !> Fonctions publiques
  public :: apply_green, apply_green_GD, apply_greenD,&
            eval_criteres, eval_criteres_GD, eval_criteresD

!------------------------------------------------------------------------------
!> Structure TIMES_GREEN

  type TIMES_GREEN
     double precision       :: apply  = 0  ! temps cumule dans apply_green
     double precision       :: crit   = 0  ! temps cumule dans eval_crit
  end type TIMES_GREEN

!------------------------------------------------------------------------------
!> Declaration des variables
  type(TIMES_GREEN),pointer :: times_g

contains

!===================================================================================================
!> Subroutine d'application du tenseur de Green
!! sans passer par le stockage de GREEN0
!! 
!! Notation CAST3M (Voigt+ordre 11 22 33 12 13 23)
!!
!!   \param[in] SigF   contrainte (notation CAST3M), espace Fourier 
!!                   pinceaux-Z (nx/2+1,ny,nz,6)
!!   \param[in] FREQ  tableau de fréquence, espace Fourier 
!!                   pinceaux-Z (nx/2+1,ny,nz,3)
!!   \param[in] LambdaMu0   Coefficients [Lambda0,mu0] du milieu de référence
!!   \param[out] SigF vaut - gamma0*SigF en sortie
!!                    avec les notations CAST3M pour la deformation
!!   \param[out] DefF_nsym gradient complet du deplacement (9 composantes)
!!               (optionnel) 
!!
!===================================================================================================
  subroutine apply_green(SigF, FREQ,LambdaMu0,DefF_nsym)

    implicit none
    real(mytype), dimension(2), intent(in) :: LambdaMu0
    complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:6), &
         intent(inout)  :: SigF
    real(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3), &
         intent(in)  :: FREQ
    complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:9),&
         intent(out),optional :: DefF_nsym
    ! variables temporaires pour la correction Paire (desuet)
    complex(mytype), allocatable, dimension(:,:,:) :: SigFx, SigFy, SigFz !gcc-warning accepte 
                                                                          !     (allocate dans if)

    integer                       :: nx2,ny2,nz2
    real(mytype)                  :: norme2, A, C, ik0, imu0      ! variables temporaires
    complex(mytype)               :: B                            ! variable temporaire
    complex(mytype)               :: trsig                        ! variable temporaire
    complex(mytype), dimension(3) :: u,fi                         ! variables temporaires
    complex(mytype)               :: imP = cmplx(0,1,kind=mytype) ! le nombre complexe i
    integer                       :: i,j,h                        ! variables temporaires
#ifdef OPENMP
    integer                       :: nbthread 
#endif
    integer                       :: ierror
    double precision              :: t1
  
    call MPI_BARRIER(amitex_mpi_comm,ierror)          
    t1 = MPI_WTIME()

#ifdef OPENMP
!$OMP PARALLEL
  nbthread = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
#endif


    ! On stocke dans les SigFi les valeurs de SigF 
    ! reutilisees dans le calcul pour le cas pair
    if(algo_param%CorrPair) then
    if(modulo(grid%nx,2) == 0) then
       nx2 = grid%nx/2 +1
       if(nx2>=fft_start(1) .AND. nx2<=fft_end(1)) then
          allocate(SigFx(fft_start(2):fft_end(2),fft_start(3):fft_end(3),6))
          do i=1,6
             do j=fft_start(3),fft_end(3)
                do h=fft_start(2),fft_end(2)
                   SigFx(h,j,i) = SigF(nx2,h,j,i)
                end do
             end do
          end do
       end if
    end if
    if(modulo(grid%ny,2) == 0) then
       ny2 = grid%ny/2 +1
       if(ny2>=fft_start(2) .AND. ny2<=fft_end(2)) then
          allocate(SigFy(fft_start(1):fft_end(1),fft_start(3):fft_end(3),6))
          do i=1,6
             do j=fft_start(3),fft_end(3)
                do h=fft_start(1),fft_end(1)
                   SigFy(h,j,i) = SigF(h,ny2,j,i)
                end do
             end do
          end do
       end if
    end if
    if(modulo(grid%nz,2) == 0) then
       nz2 = grid%nz/2 +1
       if(nz2>=fft_start(3) .AND. nz2<=fft_end(3)) then
          allocate(SigFz(fft_start(1):fft_end(1),fft_start(2):fft_end(2),6))
          do i=1,6
             do j=fft_start(2),fft_end(2)
                do h=fft_start(1),fft_end(1)
                   SigFz(h,j,i) = SigF(h,j,nz2,i)
                end do
             end do
          end do
       end if
    end if
    end if
!!!!!!!!!!!!!!!!!!
    ! Introduction de la constante C=(3K+G)/(3K+4G)
    !avec: G=mu et 3K=(3lambda+2mu)
    C=(LambdaMu0(1)+ LambdaMu0(2))/(LambdaMu0(1)+ 2._mytype*LambdaMu0(2))

#ifdef OPENMP
!$OMP PARALLEL NUM_THREADS(nbthread) private(fi,norme2,A,B,u)
!$OMP DO
#endif
    do i=fft_start(3),fft_end(3)
       do j=fft_start(2),fft_end(2)
          do h=fft_start(1),fft_end(1)
             ! fi = div(SigF)
             fi(1)= imP * (  FREQ(h,j,i,1)* SigF(h,j,i,1) +FREQ(h,j,i,2)* SigF(h,j,i,4) + FREQ(h,j,i,3)* SigF(h,j,i,5) )
             fi(2)= imP * (  FREQ(h,j,i,1)* SigF(h,j,i,4) +FREQ(h,j,i,2)* SigF(h,j,i,2) + FREQ(h,j,i,3)* SigF(h,j,i,6) )
             fi(3)= imP * (  FREQ(h,j,i,1)* SigF(h,j,i,5) +FREQ(h,j,i,2)* SigF(h,j,i,6) + FREQ(h,j,i,3)* SigF(h,j,i,3) )

             ! forces volumiques nodales : on normalise ici par ntot
             if (associated(FvolNF)) then
                fi = fi + FvolNF(h,j,i,:)*FREQ_2(h,j,i)/real(grid%ntot,mytype) 
             end if

             norme2 = FREQ(h,j,i,1) **2 +  FREQ(h,j,i,2) **2 +  FREQ(h,j,i,3) **2
             if ((h .eq. 1) .AND. (i .eq. 1) .AND. (j .eq. 1)) then
                norme2=1
             end if
             !! On fixe differentes constantes utiles au calcul
             ! A= 1/(G |FREQ|**2)
             ! B= (3K+G)/(3K+4G) / (|FREQ|**2) * prod_scal(FREQ, fi)  =>  B = C/(|FREQ|**2) prod_scal(FREQ, fi)
             ! avec:
             ! G=mu et 3K=(3lambda+2mu)
             A = 1._mytype / ( LambdaMu0(2) * norme2)
             B = ( C/ norme2) * (FREQ(h,j,i,1)*fi(1)+ FREQ(h,j,i,2)*fi(2)+ FREQ(h,j,i,3)*fi(3) )  !

             u = A*( fi - B * FREQ(h,j,i,:) )

             SigF(h,j,i,1) = imP * FREQ(h,j,i,1)*u(1)
             SigF(h,j,i,2) = imP * FREQ(h,j,i,2)*u(2)
             SigF(h,j,i,3) = imP * FREQ(h,j,i,3)*u(3)
             SigF(h,j,i,4) = imP * ( FREQ(h,j,i,1)*u(2) + FREQ(h,j,i,2)*u(1) )
             SigF(h,j,i,5) = imP * ( FREQ(h,j,i,1)*u(3) + FREQ(h,j,i,3)*u(1) )
             SigF(h,j,i,6) = imP * ( FREQ(h,j,i,2)*u(3) + FREQ(h,j,i,3)*u(2) )
             
             if (present(DefF_nsym)) then
                ! Calcul du champ complet gradient du deplacement
                DefF_nsym(h,j,i,1) = imP * FREQ(h,j,i,1)*u(1) ! grad(u)(1,1)
                DefF_nsym(h,j,i,2) = imP * FREQ(h,j,i,2)*u(2) ! grad(u)(2,2)
                DefF_nsym(h,j,i,3) = imP * FREQ(h,j,i,3)*u(3) ! grad(u)(3,3)
                DefF_nsym(h,j,i,4) = imP * FREQ(h,j,i,2)*u(1) ! grad(u)(1,2)
                DefF_nsym(h,j,i,5) = imP * FREQ(h,j,i,3)*u(1) ! grad(u)(1,3)
                DefF_nsym(h,j,i,6) = imP * FREQ(h,j,i,3)*u(2) ! grad(u)(2,3)
                DefF_nsym(h,j,i,7) = imP * FREQ(h,j,i,1)*u(2) ! grad(u)(2,1)
                DefF_nsym(h,j,i,8) = imP * FREQ(h,j,i,1)*u(3) ! grad(u)(3,1)
                DefF_nsym(h,j,i,9) = imP * FREQ(h,j,i,2)*u(3) ! grad(u)(3,2)
             end if 
          end do
       end do
    end do
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif


!!! Traitement des cas pairs
  if (algo_param%CorrPair) then 
    ! On fixe des constantes
    ik0 = 1._mytype/(LambdaMu0(1)+2._mytype*LambdaMu0(2)/3._mytype)
    imu0 = 1._mytype/(LambdaMu0(2))
    ! On teste la parite sur chaque composante
    if(modulo(grid%nx,2) == 0) then
       nx2= grid%nx/2+1
       if(nx2>=fft_start(1) .AND. nx2<=fft_end(1)) then
          !! Traitement reserve au cas pair (voir Moulinec Suquet '98)
          do j=fft_start(3),fft_end(3)
             do h=fft_start(2),fft_end(2)
                trsig = -(ik0/9._mytype - imu0/6._mytype)*(SigFx(h,j,1)+SigFx(h,j,2)+SigFx(h,j,3))
                SigF(nx2,h,j,1) = trsig - SigFx(h,j,1)*imu0/2._mytype
                SigF(nx2,h,j,2) = trsig - SigFx(h,j,2)*imu0/2._mytype
                SigF(nx2,h,j,3) = trsig - SigFx(h,j,3)*imu0/2._mytype
                SigF(nx2,h,j,4) = -SigFx(h,j,4)*imu0
                SigF(nx2,h,j,5) = -SigFx(h,j,5)*imu0
                SigF(nx2,h,j,6) = -SigFx(h,j,6)*imu0
             end do
          end do
          deallocate(SigFx)
       end if
    end if
    !! Idem pour la direction y
    if(modulo(grid%ny,2) == 0) then
       ny2= grid%ny/2+1
       if(ny2>=fft_start(2) .AND. ny2<=fft_end(2)) then
          do j=fft_start(3),fft_end(3)
             do h=fft_start(1),fft_end(1)
                trsig = -(ik0/9._mytype - imu0/6._mytype)*(SigFy(h,j,1)+SigFy(h,j,2)+SigFy(h,j,3))
                SigF(h,ny2,j,1) = trsig - SigFy(h,j,1)*imu0/2._mytype
                SigF(h,ny2,j,2) = trsig - SigFy(h,j,2)*imu0/2._mytype
                SigF(h,ny2,j,3) = trsig - SigFy(h,j,3)*imu0/2._mytype
                SigF(h,ny2,j,4) = -SigFy(h,j,4)*imu0
                SigF(h,ny2,j,5) = -SigFy(h,j,5)*imu0
                SigF(h,ny2,j,6) = -SigFy(h,j,6)*imu0
             end do
          end do
          deallocate(SigFy)
       end if
    end if
    !! Idem pour la direction z
    if(modulo(grid%nz,2) == 0) then
       nz2= grid%nz/2+1
       if(nz2>=fft_start(3) .AND. nz2<=fft_end(3)) then
          do j=fft_start(2),fft_end(2)
             do h=fft_start(1),fft_end(1)
                trsig = -(ik0/9._mytype - imu0/6._mytype)*(SigFz(h,j,1)+SigFz(h,j,2)+SigFz(h,j,3))
                SigF(h,j,nz2,1) = trsig - SigFz(h,j,1)*imu0/2._mytype
                SigF(h,j,nz2,2) = trsig - SigFz(h,j,2)*imu0/2._mytype
                SigF(h,j,nz2,3) = trsig - SigFz(h,j,3)*imu0/2._mytype
                SigF(h,j,nz2,4) = -SigFz(h,j,4)*imu0
                SigF(h,j,nz2,5) = -SigFz(h,j,5)*imu0
                SigF(h,j,nz2,6) = -SigFz(h,j,6)*imu0
             end do
          end do
          deallocate(SigFz)
       end if
    end if
   end if ! fin if cas_pair pris en compte
   
  call MPI_BARRIER(amitex_mpi_comm,ierror)          
  times_g%apply = times_g%apply + MPI_WTIME() - t1

  end subroutine apply_green
!==================================================================


!====================================================================
!> Subroutine d'application du tenseur de Green en grandes transformations
!! sans passer par le stockage de GREEN0
!! 
!! Notation : ordre 11 22 33 12 13 23 21 31 32
!!
!!   \param[in] PiF   contrainte (notation CAST3M), espace Fourier 
!!                   pinceaux-Z (nx/2+1,ny,nz,9)
!!   \param[in] FREQ  tableau de fréquence, espace Fourier 
!!                   pinceaux-Z (nx/2+1,ny,nz,3)
!!   \param[in] LambdaMu0   Coefficients [Lambda0,mu0] du milieu de référence
!!   \param[out] PiF vaut - gamma0*PiF en sortie
!!                    avec les notations CAST3M pour la deformation
!!
!!
!! REMARQUE : 12/12/2018, prise en compte option d'un C0 symetrique
!====================================================================
  subroutine apply_green_GD(PiF, FREQ,LambdaMu0)

    implicit none

    real(mytype), dimension(2), intent(in) :: LambdaMu0
    complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:9), &
         intent(inout)  :: PiF
    real(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3), &
         intent(in)  :: FREQ
    real(mytype) :: norme2, A, C ! variables temporaires
    complex(mytype) :: B ! variable temporaire
    complex(mytype), dimension(3) :: u,fi ! variables temporaires
    complex(mytype) :: imP = cmplx(0,1,kind=mytype) ! le nombre complexe i
    integer :: i,j,h ! variables temporaires
#ifdef OPENMP
    integer                       :: nbthread
#endif
    integer                       :: ierror
    double precision              :: t1
  
    call MPI_BARRIER(amitex_mpi_comm,ierror)          
    t1 = MPI_WTIME()

#ifdef OPENMP
!$OMP PARALLEL
  nbthread = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
#endif


    ! Introduction de la constante C
    if(algo_param%C0sym) then    !-> version symetrique : C0:gradU_sym
        C=(LambdaMu0(1)+LambdaMu0(2))/(LambdaMu0(1)+ 2._mytype*LambdaMu0(2))
    else                         !-> version non symetrique : C0:gradU
        C=LambdaMu0(1)/(LambdaMu0(1)+ 2._mytype*LambdaMu0(2)) 
    end if
#ifdef OPENMP
!$OMP PARALLEL NUM_THREADS(nbthread) private(fi,norme2,A,B,u)
!$OMP DO
#endif
    do i=fft_start(3),fft_end(3)
       do j=fft_start(2),fft_end(2)
          do h=fft_start(1),fft_end(1)
             !!!!!!!!!!!!!!!!!!!!!!!!!!! A REFAIRE 
             !  fi = div(PiF)
             !! fi = i Freq_i*\hat{\Pi}_i = i(Freq_i(1)*pi_i(1,1) + Freq_i(2)*pi_i(1,2)
             !!        + Freq_i(3)*pi_i(1,3)) e_1 +...

             fi(1)= imP * (  FREQ(h,j,i,1)*PiF(h,j,i,1)+FREQ(h,j,i,2)*PiF(h,j,i,4)+FREQ(h,j,i,3)*PiF(h,j,i,5) )
             fi(2)= imP * (  FREQ(h,j,i,1)*PiF(h,j,i,7)+FREQ(h,j,i,2)*PiF(h,j,i,2)+FREQ(h,j,i,3)*PiF(h,j,i,6) )
             fi(3)= imP * (  FREQ(h,j,i,1)*PiF(h,j,i,8)+FREQ(h,j,i,2)*PiF(h,j,i,9)+FREQ(h,j,i,3)*PiF(h,j,i,3) )

             ! forces volumiques nodales : on normalise ici par ntot
             if (associated(FvolNF)) then
                fi = fi + FvolNF(h,j,i,:)*FREQ_2(h,j,i)/real(grid%ntot,mytype) 
             end if

             norme2 = FREQ(h,j,i,1) **2 +  FREQ(h,j,i,2) **2 +  FREQ(h,j,i,3) **2
             if ((h .eq. 1) .AND. (i .eq. 1) .AND. (j .eq. 1)) then
                norme2=1
             end if
             !! On fixe differentes constantes utiles au calcul
             ! A= 1/(G |FREQ|**2)
             ! B= (lambda)/(lambda+2mu) / (|FREQ|**2) * prod_scal(FREQ, fi) =>  B = C/(|FREQ|**2) prod_scal(FREQ, fi)
             ! G=2mu et C=lambda/(lambda+2mu)
             if(algo_param%C0sym) then   !-> version symetrique : C0:gradU_sym
                A = 1._mytype / (LambdaMu0(2) * norme2) 
             else                        !-> version non symetrique : C0:gradU
                A = 1._mytype / (2._mytype* LambdaMu0(2) * norme2) 
             end if

             B = ( C/ norme2) * (FREQ(h,j,i,1)*fi(1)+ FREQ(h,j,i,2)*fi(2)+ FREQ(h,j,i,3)*fi(3) )  !

             u = A*( fi - B * FREQ(h,j,i,:) )

             PiF(h,j,i,1) = imP * FREQ(h,j,i,1)*u(1) ! grad(u)(1,1)
             PiF(h,j,i,2) = imP * FREQ(h,j,i,2)*u(2) ! grad(u)(2,2)
             PiF(h,j,i,3) = imP * FREQ(h,j,i,3)*u(3) ! grad(u)(3,3)
             PiF(h,j,i,4) = imP * FREQ(h,j,i,2)*u(1) ! grad(u)(1,2)
             PiF(h,j,i,5) = imP * FREQ(h,j,i,3)*u(1) ! grad(u)(1,3)
             PiF(h,j,i,6) = imP * FREQ(h,j,i,3)*u(2) ! grad(u)(2,3)
             PiF(h,j,i,7) = imP * FREQ(h,j,i,1)*u(2) ! grad(u)(2,1)
             PiF(h,j,i,8) = imP * FREQ(h,j,i,1)*u(3) ! grad(u)(3,1)
             PiF(h,j,i,9) = imP * FREQ(h,j,i,2)*u(3) ! grad(u)(3,2)
          end do
       end do
    end do
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif

  call MPI_BARRIER(amitex_mpi_comm,ierror)          
  times_g%apply = times_g%apply + MPI_WTIME() - t1

  end subroutine apply_green_GD
!==================================================================



!====================================================================
!> SUBROUTINE D'EVALUATION DU CRITERE D'EQUILIBRE ET 
!!            DU CRITERE SUR LES VALEURS MOYENNES EN PETITES PERTURBATIONS
!!
!!          format (pinceaux-Z, taille globale (nx/2+1)*ny*nz)
!!
!! Notation CAST3M (Voigt + ordre 11,22,33,12,13,23)
!!
!! 
!! \param[in]      SigF: (complexe) champs des contraintes dans l'espace de Fourier
!! \param[in]      DefF: (complexe) champs des deformations dans l'espace de Fourier
!!                                  DeF est ici prealablement normalise par ntot
!! \param[in]      FREQ: (reel) champs des frequences
!! \param[in]      local_loading   chargement impose
!!                              -  local_loading(1:6) : pilotage
!!                              -  local_loading(7:12) : valeur du pilotage
!! \param[in]      ntot: (entier) nombre de cellules dans le domaine si on veut calculer le critere sur la deformation
!!                                0 si on veut calculer les criteres sur la contrainte
!! \param[in]      DirStress_flag   .true. si contrainte multiaxiale impose avec pilotage deformation
!! \param[in]      DirStress        (vecteur 6) direction de tenseur de contrainte correspondant
!! \param[out]     criteq: (reel) critere d'equilibre
!! \param[out]     critSigMoy: (reel) ecart entre les contraintes moyennes obtenues 
!!                                 et les contraintes moyennes imposees
!! \param[out]     critDefMoy: (reel) ecart entre les deformations moyennes obtenues 
!!                                 et les deformations moyennes imposees
!! \param[out]     critCptb: (reel) critere de compatibilite
!!
!====================================================================
  subroutine eval_criteres(SigF, DefF, FREQ, local_loading,&
       DirStress_flag, DirStress,& 
       criteq,critSigMoy,critDefMoy,critCptb,ntot)

    implicit none

    complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:6),&
         intent(inout)  :: SigF
    complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:6),&
         intent(in)  :: DefF
    real(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3),&
         intent(in)  :: FREQ
    real(mytype), dimension(12), intent(in)     :: local_loading
    logical, intent(in)                         :: DirStress_flag
    real(mytype), dimension(6), intent(in)      :: DirStress
    !> nombre de voxels de la cellule ou 0 si on ne souhaite pas calculer le critere en compatbilite
    integer(kind=8), intent(in)                 :: ntot
    real(mytype),dimension(6)     :: alpha,tmp6
    real(mytype),dimension(1)     :: tmp1
    real(mytype)                  :: alpha_max
    complex(mytype), dimension(3) :: fi
    integer                       :: i,j,h, ierr, nb_stress_pilot,n0x,n0y,n0z
    integer(kind=8)               :: cell_size
#ifdef OPENMP
    integer                       :: nbthread
#endif
    real(mytype),intent(out)      :: criteq, critDefMoy, critSigMoy,critCptb
    real(mytype)                  :: sum2Sig, sum2forces,sum2div,tmp,normeDDef,normeDSig,tmpi,&
                                     sum2Def,sum2Icomp,dR
#ifdef DOUBLE_PREC
    real(mytype),parameter        :: epsilon = 10._mytype**(-300)
#else
    real(mytype),parameter        :: epsilon = 10._mytype**(-30)
#endif
    complex(mytype),parameter     :: imP = cmplx(0,1,kind=mytype) ! le nombre complexe i
    integer                       :: ierror
    double precision              :: t1
  
    call MPI_BARRIER(amitex_mpi_comm,ierror)          
    t1 = MPI_WTIME()

#ifdef OPENMP
!$OMP PARALLEL
  nbthread = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
#endif


    cell_size= grid%ntot
    nb_stress_pilot = count(local_loading(1:6) == 1)
    if(ntot>0) then
       !! Si on rentre dans cette fonction avec computeDef=1 on calcule
       !!========================================================================================
       !!                                                                 CRITERES EN DEFORMATION
       !!
       !! On veut calculer ici 2 criteres differents :
       !! - le critere sur la def moyenne
       !! - critCptb est le critere de compatibilite
       !!
       !! ATTNETION : en entree DefF est divise par ntot 
       !!            -> pas besoin de diviser a nouveau pour obtenir la deformation moyenne 
       !!            -> corriger la somme des carres
       !!========================================================================================
       !!--------------------------------------------------------------------------
       !!                                        SOMME DES CARRES DE LA DEFORMATION 
       !!                                        sum2Def =  \sum_{n=0}^{N} |x[n]|^2
       !!
       !! ATTENTION : DefF est normalise par Ntot -> penser a corriger 
       !!
       !! Dans le cas ou on impose un chargement thermique avec deformation ou 
       !! contrainte moyenne nulle, la premiere valeur du champ de deformation proposee
       !! par l'algorithme est identiquement nulle. Ainsi, on ne peut pas normaliser
       !! par cette valeur => seuil
       tmp6=0
       call field_sum_squareF(DefF,tmp6,grid%nx,grid%ny,grid%nz,6,1)
       tmp6=tmp6 * grid%ntot * grid%ntot !--- correction car DefF est ici normalise par ntot
       sum2Def = sum(tmp6(1:3)) + sum(tmp6(4:6)) / 2._mytype
       if(sum2Def<epsilon) sum2Def = epsilon

       !!--------------------------------------------------------------------------
       !!                                        CRITERE SUR LA DEFORMATION MOYENNE
       !!                         critDefMoy = norme(Defmoy) / sqrt(sum2Def * dV/V)
       !!                                    = norme(Defmoy) / sqrt(sum2Def / ntot)
       !!      
       !!             Si aucune deformation moyenne n'est imposee le critere vaut 0
       critDefMoy = 0
       if(nb_stress_pilot<6) then
          !! la deformation moyenne est imposee dans au moins une direction.
          !! On veut calculer l'ecart relatif entre la deformation moyenne
          !! et celle que l'on veut imposer (pour certaines directions)
          normeDDef = 0 
          tmp = 0
          if((fft_start(1) .eq. 1) .AND. (fft_start(2) .eq. 1) .AND. (fft_start(3) .eq. 1)) then
             !! normeDDef est la norme 2 de l'ecart entre la deformation moyenne
             !! et la deformation imposee (sur les composantes concernees)
             !! Prise en compte du fait que DefF est divise par ntot
             do i=1,3
                if(local_loading(i) == 0) then
                   normeDDef = normeDDef + (real(DefF(1,1,1,i),mytype) - &
                        local_loading(i+6))**2
                end if
             end do
             do i=4,6
                if(local_loading(i) == 0) then
                   normeDDef = normeDDef + 0.5_mytype*(real(DefF(1,1,1,i),mytype) - &
                        local_loading(i+6))**2
                end if
             end do
             tmp = sqrt(normeDDef/(sum2Def/grid%ntot))
             !! On normalise par la norme L2 de la deformation
          end if
          !! On recupere le max sur tous les processus pour que tous les processus
          !! aient la bonne valeur de critDefMoy
          call MPI_AllReduce(tmp,critDefMoy,1,real_type,MPI_MAX,amitex_mpi_comm,ierr)
       end if


       !!--------------------------------------------------------------------------
       !!                                        CALCUL DU CRITERE DE COMPATIBILITE         
       !!                  
       !!            critcptb = sqrt(somme des carres du tenseur d'incompatibilite )
       !!                   / sqrt(somme des carres de la deformation)
       !!   
       !!            crtcptb = cricptb*Tmax*Tmax
       !! 
       !! Faute d'argument tel que pour le critere d'equilibre, 
       !! on assure l'independance a la dimension en multipliant pat Tmax^2
       !! on aurait pu choisir dx^2. Le choix realise a le merite d'etre conservatif.  
       !!
       !! Traitement eventuel des cas pairs (a supprimer, plus utilise) :
       !!
       !! On elimine les frequences introduites dans le cas des resolutions paires
       !! En fait, on definit n0i qui vaut 
       !! - n0i = -1 si le nombre de voxels dans la direction i est impair
       !! - n0i = ni/2+1 si le nombre de voxels dans la direction i est pair
       !! Si on est sur un des indices n0i on ne fait pas le calcul
       !! dR=dR+LeviCivita(r,p,l)*LeviCivita(k,q,m)*Def(l,m)*FREQ(h,j,i,p)*FREQ(h,j,i,q)  
       !! ce qui est equivalent a prendre des frequences nulles sur ces points
       !! Ce cas n'arrive pas si n0i = -1
       n0x=-1
       n0y=-1
       n0z=-1

       ! traitemen dimensions impaires
       if (algo_param%CorrPair) then      
       if (mod(grid%nx,2) .eq. 0) then
          n0x = grid%nx/2+1
       end if
       if (mod(grid%ny,2) .eq. 0) then
          n0y = grid%ny/2+1
       end if
       if (mod(grid%nz,2) .eq. 0) then
          n0z = grid%nz/2+1
       end if
       end if

       tmp=0
#ifdef OPENMP
!$OMP PARALLEL NUM_THREADS(nbthread) private(dR)
!$OMP DO REDUCTION(+:tmp)
#endif
       do i=fft_start(3),fft_end(3)
          if(i .EQ. n0z) cycle
          do j=fft_start(2),fft_end(2)
             if(j .EQ. n0y) cycle
             do h=fft_start(1),fft_end(1)
                if(h .EQ. n0x) cycle
!                tmp=0

                tmpi = 2._mytype   !tmpi permet de prendre en compte la somme totale a partir d'un espace (1:nx/2+1,1:ny,1:nz)
                if (h == 1) tmpi=1._mytype
                if (modulo(grid%nx,2) == 0 .AND. h == grid%nx / 2 + 1) tmpi = 1._mytype

                dR = real(DefF(h,j,i,3),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,2)&
                     + real(DefF(h,j,i,2),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,3)&
                     - real(DefF(h,j,i,6),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,3)
                tmp = tmp + (tmpi*dR**2)/grid%ntot
                dR = real(DefF(h,j,i,1),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,3) +&
                     real(DefF(h,j,i,3),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,1) -&
                     real(DefF(h,j,i,5),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,1)
                tmp = tmp + (tmpi*dR**2)/grid%ntot
                dR = real(DefF(h,j,i,1),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,2) +&
                     real(DefF(h,j,i,2),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,1) -&
                     real(DefF(h,j,i,4),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,2)
                tmp = tmp + (tmpi*dR**2)/grid%ntot
                dR = -real(DefF(h,j,i,1),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,2) +&
                     0.5_mytype*real(DefF(h,j,i,5),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,1) -&
                     0.5_mytype*real(DefF(h,j,i,6),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,1) +&
                     0.5_mytype*real(DefF(h,j,i,4),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,1)
                tmp = tmp + (tmpi*2._mytype*dR**2)/grid%ntot
                dR = -real(DefF(h,j,i,2),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,3) +&
                     0.5_mytype*real(DefF(h,j,i,6),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,1) -&
                     0.5_mytype*real(DefF(h,j,i,5),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,2) +&
                     0.5_mytype*real(DefF(h,j,i,4),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,3)
                tmp = tmp + (tmpi*2._mytype*dR**2)/grid%ntot
                dR = -real(DefF(h,j,i,3),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,2) +&
                     0.5_mytype*real(DefF(h,j,i,5),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,2) -&
                     0.5_mytype*real(DefF(h,j,i,4),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,3) +&
                     0.5_mytype*real(DefF(h,j,i,6),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,3)
                tmp = tmp + (tmpi*2._mytype*dR**2)/grid%ntot
                !!Parties imaginaires
                dR = real(aimag(DefF(h,j,i,3)),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,2) +&
                     real(aimag(DefF(h,j,i,2)),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,3) -&
                     real(aimag(DefF(h,j,i,6)),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,3)
                tmp = tmp + (tmpi*dR**2)/grid%ntot
                dR = real(aimag(DefF(h,j,i,1)),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,3) +&
                     real(aimag(DefF(h,j,i,3)),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,1) -&
                     real(aimag(DefF(h,j,i,5)),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,1)
                tmp = tmp + (tmpi*dR**2)/grid%ntot
                dR = real(aimag(DefF(h,j,i,1)),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,2) +&
                     real(aimag(DefF(h,j,i,2)),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,1) -&
                     real(aimag(DefF(h,j,i,4)),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,2)
                tmp = tmp + (tmpi*dR**2)/grid%ntot
                dR = -real(aimag(DefF(h,j,i,1)),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,2) +&
                     0.5_mytype*real(aimag(DefF(h,j,i,5)),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,1) -&
                     0.5_mytype*real(aimag(DefF(h,j,i,6)),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,1) +&
                     0.5_mytype*real(aimag(DefF(h,j,i,4)),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,1)
                tmp = tmp + (tmpi*2._mytype*dR**2)/grid%ntot
                dR = -real(aimag(DefF(h,j,i,2)),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,3) +&
                     0.5_mytype*real(aimag(DefF(h,j,i,6)),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,1) -&
                     0.5_mytype*real(aimag(DefF(h,j,i,5)),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,2) +&
                     0.5_mytype*real(aimag(DefF(h,j,i,4)),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,3)
                tmp = tmp + (tmpi*2._mytype*dR**2)/grid%ntot
                dR = -real(aimag(DefF(h,j,i,3)),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,2) +&
                     0.5_mytype*real(aimag(DefF(h,j,i,5)),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,2) -&
                     0.5_mytype*real(aimag(DefF(h,j,i,4)),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,3) +&
                     0.5_mytype*real(aimag(DefF(h,j,i,6)),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,3)
                tmp = tmp + (tmpi*2._mytype*dR**2)/grid%ntot
             end do
          end do
       end do
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif

       call MPI_ALLREDUCE(tmp,sum2Icomp,1,real_type,MPI_SUM,amitex_mpi_comm,ierr)
       sum2Icomp = sum2Icomp * grid%ntot * grid%ntot !--- Correction car defF normalisee (voir correction sum2Def)
       critCptb = sqrt(sum2Icomp/sum2Def)

       ! NORMALISATION du critere 
       ! => assurer l'independance au choix de dx,dy,dz
       ! => introduction possible d'une dimension "physique"
       critCptb = critCptb * grid%Tmax * grid%Tmax

    else
       !!========================================================================================
       !!                                                                  CRITERES EN CONTRAINTE
       !!
       !! On veut calculer ici 2 criteres differents :
       !! - criteq est le critere d'equilibre : norm(div(sigma).dV)_L^2/norm(sigma.n.dS)_L^2
       !! - critSigMoy est le critere sur la valeur moyenne de la contrainte
       !!========================================================================================
       !!--------------------------------------------------------------------------
       !!                                         SOMME DES CARRES DE LA CONTRAINTE 
       !!                                        sum2Sig =  \sum_{n=0}^{N} |x[n]|^2
       !!
       tmp6=0
       call field_sum_squareF(SigF,tmp6,grid%nx,grid%ny,grid%nz,6,1)
       sum2Sig = sum(tmp6) + sum(tmp6(4:6))

       !!--------------------------------------------------------------------------
       !!                  CALCUL DU CRITERE SUR LA VALEUR MOYENNE DE LA CONTRAINTE
       !!
       !!                         critSigMoy = norme(Sigmoy) / sqrt(sum2Sig * dV/V)
       !!                                    = norme(Sigmoy) / sqrt(sum2Sig / ntot)
       !!
       !!    remarque : ce critere est calcule sur le pinceau contenant (1,1,1)
       !!               il faut redistribuer

       critSigMoy = 0
       normeDSig = 0
       tmp = 0
       !! Si aucune contrainte n'est imposee ce critere vaut 0
       !! Cas d'un chargement par composante
       if (.not. DirStress_flag) then
       if((nb_stress_pilot>0) .AND. (fft_start(1) .eq. 1) .AND. (fft_start(2) .eq. 1) .AND. (fft_start(3) .eq. 1)) then
          !! normeDSig est la norme 2 de l'ecart entre la contrainte moyenne
          !! et la contrainte imposee (sur les composantes concernees)
          do i=1,3
             if(local_loading(i) == 1) then
                normeDSig = normeDSig + (real(SigF(1,1,1,i),mytype)/grid%ntot - &
                     local_loading(i+6))**2
             end if
          end do
          do i=4,6
             if(local_loading(i) == 1) then
                normeDSig = normeDSig + 2._mytype*(real(SigF(1,1,1,i),mytype)/grid%ntot - &
                     local_loading(i+6))**2
             end if
          end do
          tmp = sqrt(normeDSig/(sum2Sig/grid%ntot))
       end if
       !! Cas d'un chargement multiaxial en contrainte pilote en deformation
       elseif ((fft_start(1) .eq. 1) .AND. (fft_start(2) .eq. 1) .AND. (fft_start(3) .eq. 1)) then 
          alpha = 0.
          where (DirStress .ne. 0)
             alpha = (real(SigF(1,1,1,:),mytype) / grid%ntot) / DirStress
          end where
          if (abs(maxval(alpha)) >= abs(minval(alpha))) then
             alpha_max = maxval(alpha)
          else
             alpha_max = minval(alpha)
          end if
          do i=1,3
                normeDSig = normeDSig + (real(SigF(1,1,1,i),mytype)/ grid%ntot - &
                     alpha_max*DirStress(i))**2
          end do
          do i=4,6
                normeDSig = normeDSig + 2._mytype*(real(SigF(1,1,1,i),mytype)/ grid%ntot - &
                     alpha_max*DirStress(i))**2
          end do
          tmp = sqrt(normeDSig/(sum2Sig/grid%ntot))
       end if
       !! On recupere le max sur tous les processus pour que tous les processus
       !! aient la bonne valeur de critSigMoy
       call MPI_AllReduce(tmp,critSigMoy,1,real_type,MPI_MAX,amitex_mpi_comm,ierr)   

       !!--------------------------------------------------------------------------
       !!                      CALCUL DU CRITERE D'EQUILIBRE         
       !!                  
       !!            criteq = sqrt(somme des carres de la force feq)
       !!                   / sqrt(somme des carres des forces f1,f2,f3)
       !!
       !!            feq = div(sig)*dV
       !!            f1 = Sig.e1.dx2.dx3
       !!            f2 = Sig.e2.dx1.dx3
       !!            f3 = Sig.e3.dx1.dx2
       !!
       !!            ci-dessous : criteq = sqrt(sum2div)*dV / sqrt(sum2forces)
       !!

       ! NUMERATEUR : somme sur les points des carres de div(sigma) (a multiplier par dV -> etape finale)
       criteq=0 !ici role de variable temporaire
#ifdef OPENMP
!$OMP PARALLEL NUM_THREADS(nbthread) private(tmp,fi)
!$OMP DO REDUCTION(+:criteq)
#endif
       do i=fft_start(3),fft_end(3)
          tmp=0
         do j=fft_start(2),fft_end(2)
             do h=fft_start(1),fft_end(1)
                ! fi = div(SigF)
                fi(1)=  FREQ(h,j,i,1)* SigF(h,j,i,1) +FREQ(h,j,i,2)*SigF(h,j,i,4) + FREQ(h,j,i,3)* SigF(h,j,i,5)
                fi(2)=  FREQ(h,j,i,1)* SigF(h,j,i,4) +FREQ(h,j,i,2)*SigF(h,j,i,2) + FREQ(h,j,i,3)* SigF(h,j,i,6)
                fi(3)=  FREQ(h,j,i,1)* SigF(h,j,i,5) +FREQ(h,j,i,2)*SigF(h,j,i,6) + FREQ(h,j,i,3)* SigF(h,j,i,3)

                ! prise en compte des forces volumiques nodales (ici on a i.fi, ce qui ne change pas la norme)
                if (associated(FvolNF)) fi = -fi + imP*FvolNF(h,j,i,:) * FREQ_2(h,j,i)          

                tmpi = 2._mytype
                if (h == 1) tmpi=1._mytype
                if (modulo(grid%nx,2) == 0 .AND. h == grid%nx / 2 + 1) tmpi = 1._mytype
                tmp=tmp+tmpi * real( fi(1)*conjg(fi(1))+fi(2)*conjg(fi(2))+fi(3)*conjg(fi(3)),mytype) / grid%ntot
             end do
          end do
          criteq=criteq+tmp 
       end do
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
       call MPI_AllReduce(criteq,sum2div,1,real_type,MPI_SUM,amitex_mpi_comm,ierr)

       ! DENOMINATEUR : somme des carres des forces Sig.e1.dx2.dx3, Sig.e2.dx1.dx3 et Sig.e3.dx1.dx2
       sum2forces = 0._mytype
       tmp1 = 0._mytype
       call field_sum_squareF(SigF(:,:,:,1),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dy*grid%dz)**2
       call field_sum_squareF(SigF(:,:,:,4),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dy*grid%dz)**2
       call field_sum_squareF(SigF(:,:,:,5),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dy*grid%dz)**2

       call field_sum_squareF(SigF(:,:,:,4),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dx*grid%dz)**2
       call field_sum_squareF(SigF(:,:,:,2),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dx*grid%dz)**2
       call field_sum_squareF(SigF(:,:,:,6),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dx*grid%dz)**2

       call field_sum_squareF(SigF(:,:,:,5),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dx*grid%dy)**2
       call field_sum_squareF(SigF(:,:,:,6),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dx*grid%dy)**2
       call field_sum_squareF(SigF(:,:,:,3),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dx*grid%dy)**2

       ! etape finale 
       if (sum2forces == 0.) sum2forces = 1.     ! evite la division par zero
       !criteq=sqrt(sum2div/sum2Sig) * grid%dx   !valable pour voxels cubiques
       criteq=sqrt(sum2div/sum2forces) *grid%dV  !voxels anisotropes (valide voxels composite ci-dessus)

    end if

  call MPI_BARRIER(amitex_mpi_comm,ierror)          
  times_g%crit = times_g%crit + MPI_WTIME() - t1

  end subroutine eval_criteres
!====================================================================


!====================================================================
!====================================================================
!> SUBROUTINE D'EVALUATION DU CRITERE D'EQUILIBRE ET 
!!            DU CRITERE SUR LES VALEURS MOYENNES EN GRANDES TRANSFORMATIONS
!!
!!          format (pinceaux-Z, taille globale (nx/2+1)*ny*nz)
!!
!! Ordre 11 22 33 12 13 23 21 31 32
!!
!! 
!! \param[in]      PiF: (complexe) champs des contraintes Piola Kirchhoff dans l'espace de Fourier
!! \param[in]      DefF: (complexe) champs des deformations dans l'espace de Fourier
!! \param[in]      FREQ: (reel) champs des frequences
!! \param[in]      local_loading   chargement impose
!!                              -  local_loading(1:6) : pilotage
!!                              -  local_loading(7:12) : valeur du pilotage
!! \param[in]      ntot: (entier) nombre de cellules dans le domaine si on veut calculer 
!!                       le critere sur la deformation,
!!                       0 si on veut calculer les criteres sur la contrainte
!! \param[in]      DirStress_flag   .true. si contrainte multiaxiale impose avec pilotage deformation
!! \param[in]      DirStress        (vecteur 9) direction de tenseur de contrainte correspondant
!!
!! \param[out]     criteq: (reel) critere d'equilibre
!! \param[out]     critPiMoy: (reel) ecart entre les contraintes moyennes obtenues 
!!                                 et les contraintes moyennes imposees
!! \param[out]     critDefMoy: (reel) ecart entre les deformations moyennes obtenues 
!!                                 et les deformations moyennes imposees
!! \param[out]     critCptb: (reel) critere de compatibilite
!!
!====================================================================
  subroutine eval_criteres_GD(PiF, DefF, FREQ, local_loading,&
       DirStress_flag, DirStress,& 
       criteq,critPiMoy,critDefMoy,critCptb,ntot)

    implicit none

    complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:9),&
         intent(inout)  :: PiF
    complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:9),&
         intent(in)  :: DefF
    real(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3),&
         intent(in)  :: FREQ
    real(mytype),dimension(18), intent(in)  :: local_loading
    !> nombre de voxels de la cellule ou 0 si on ne souhaite pas calculer le critere en compatbilite
    integer(kind=8), intent(in)                     :: ntot
    logical, intent(in)                             :: DirStress_flag
    real(mytype), dimension(9),intent(in)           :: DirStress
   
    complex(mytype), dimension(3) :: fi
    integer                       :: i,j,h,k,l,m,p,q,r,ierr, nb_stress_pilot,n0x,n0y,n0z
    integer(kind=8)               :: cell_size
#ifdef OPENMP
    integer                       :: nbthread
#endif
    real(mytype),intent(out)      :: criteq, critDefMoy, critPiMoy,critCptb
    real(mytype)                  :: sum2Sig,sum2div,sum2forces,tmp,tmpi,w,&
                                     sum2def,sum2Icomp,normeDPi,normeDDef
    real(mytype),dimension(9)     :: alpha,tmp9
    real(mytype),dimension(1)     :: tmp1
    real(mytype)                  :: alpha_max
    complex(mytype)               :: imP = cmplx(0,1,kind=mytype) ! le nombre complexe i
    integer(kind=1),dimension(3,3,3) :: MatLeviC
    integer(kind=1)                  :: prodMatLC
    real(mytype),dimension(3,3)      :: MatDefR, MatDefI,dR,dI
#ifdef DOUBLE_PREC
    real(mytype),parameter                     :: epsilon = 10._mytype**(-300)
#else
    real(mytype),parameter                     :: epsilon = 10._mytype**(-30)
#endif
    integer                       :: ierror
    double precision              :: t1
  
    call MPI_BARRIER(amitex_mpi_comm,ierror)          
    t1 = MPI_WTIME()

#ifdef OPENMP
    !$OMP PARALLEL
    nbthread = OMP_GET_NUM_THREADS()
    !$OMP END PARALLEL
#endif

    !! Initialisation de la matrice de Levi Civita
    MatLeviC = 0
    MatLeviC(1,2,3)=1
    MatLeviC(2,3,1)=1
    MatLeviC(3,1,2)=1
    MatLeviC(3,2,1)=-1
    MatLeviC(1,3,2)=-1
    MatLeviC(2,1,3)=-1

    cell_size= grid%ntot
    nb_stress_pilot = count(local_loading(1:9) == 1)
    if(ntot>0) then
       !! Si on rentre dans cette fonction avec computeDef=1 on calcule
       !!========================================================================================
       !!                                                                 CRITERES EN DEFORMATION
       !!
       !! On veut calculer ici 2 criteres differents :
       !! - le critere sur la def moyenne
       !! - critCptb est le critere de compatibilite
       !!
       !! ATTENTION : en entree DefF est divise par ntot 
       !!            -> pas besoin de diviser a nouveau pour obtenir la deformation moyenne 
       !!            -> corriger la somme des carres
       !!========================================================================================
       !!--------------------------------------------------------------------------
       !!                                        SOMME DES CARRES DE LA DEFORMATION 
       !!                                        sum2Def =  \sum_{n=0}^{N} |x[n]|^2
       !!
       !! ATTENTION : DefF est normalise par Ntot-> a corriger
       !!
       !! Dans le cas ou on impose un chargement thermique avec deformation ou 
       !! contrainte moyenne nulle, la premiere valeur du champ de deformation proposee
       !! par l'algorithme est identiquement nulle. Ainsi, on ne peut pas normaliser
       !! par cette valeur => seuil
       tmp9=0
       call field_sum_squareF(DefF,tmp9,grid%nx,grid%ny,grid%nz,9,1)
       tmp9=tmp9 * grid%ntot * grid%ntot !--- correction car DefF est ici normalise par ntot
       sum2Def = sum(tmp9)
       if(sum2Def<epsilon) sum2Def = epsilon

       !!--------------------------------------------------------------------------
       !!                                        CRITERE SUR LA DEFORMATION MOYENNE
       !!                         critDefMoy = norme(Defmoy) / sqrt(sum2Def * dV/V)
       !!                                    = norme(Defmoy) / sqrt(sum2Def / ntot)
       !!      
       !!             Si aucune deformation moyenne n'est imposee le critere vaut 0
       critDefMoy = 0
       if(nb_stress_pilot<9) then
          !! la deformation moyenne est imposee dans au moins une direction.
          !! On veut calculer l'ecart relatif entre la deformation moyenne
          !! et celle que l'on veut imposer (pour certaines directions)
          normeDDef = 0 
          tmp = 0
          if((fft_start(1) .eq. 1) .AND. (fft_start(2) .eq. 1) .AND. (fft_start(3) .eq. 1)) then
             !! normeDDef est la norme 2 de l'ecart entre la deformation moyenne
             !! et la deformation imposee (sur les composantes concernees)
             !! Prise en compte du fait que DefF est divise par ntot 
             do i=1,9
                if(local_loading(i) == 0) then
                   normeDDef = normeDDef + (real(DefF(1,1,1,i),mytype) - &
                        local_loading(i+9))**2
                end if
             end do
             tmp = sqrt(normeDDef/(sum2Def/grid%ntot))
             !! On normalise par la norme L2 de la deformation
          end if
          !! On recupere le max sur tous les processus pour que tous les processus
          !! aient la bonne valeur de critDefMoy
          call MPI_AllReduce(tmp,critDefMoy,1,real_type,MPI_MAX,amitex_mpi_comm,ierr)
       end if

       !!--------------------------------------------------------------------------
       !!                                        CALCUL DU CRITERE DE COMPATIBILITE         
       !!                  
       !!            critcptb = sqrt(somme des carres du tenseur d'incompatibilite )
       !!                   / sqrt(somme des carres de la deformation)
       !!   
       !!            crtcptb = cricptb*Tmax*Tmax
       !! 
       !! Faute d'argument tel que pour le critere d'equilibre, 
       !! on assure l'independance a la dimension en multipliant pat Tmax^2
       !! on aurait pu choisir dx^2. Le choix realise a le merite d'etre conservatif.  
       !!
       !! Traitement eventuel des cas pairs (a supprimer, plus utilise) :
       !!
       !! On elimine les frequences introduites dans le cas des resolutions paires
       !! En fait, on definit n0i qui vaut 
       !! - n0i = -1 si le nombre de voxels dans la direction i est impair
       !! - n0i = ni/2+1 si le nombre de voxels dans la direction i est pair
       !! Si on est sur un des indices n0i on ne fait pas le calcul
       !! dR=dR+LeviCivita(r,p,l)*LeviCivita(k,q,m)*Def(l,m)*FREQ(h,j,i,p)*FREQ(h,j,i,q)  
       !! ce qui est equivalent a prendre des frequences nulles sur ces points
       !! Ce cas n'arrive pas si n0i = -1
       n0x=-1
       n0y=-1
       n0z=-1   
       ! traitement dimensions impaires
       if (algo_param%CorrPair) then    
       if (mod(grid%nx,2) .eq. 0) then
          n0x = grid%nx/2+1
       end if
       if (mod(grid%ny,2) .eq. 0) then
          n0y = grid%ny/2+1
       end if
       if (mod(grid%nz,2) .eq. 0) then
          n0z = grid%nz/2+1
       end if
       end if
       tmp=0

#ifdef OPENMP
       !$OMP PARALLEL NUM_THREADS(nbthread) private(MatDefR,MatDefI,dR,dI,w)
       !$OMP DO REDUCTION(+:tmp)
#endif
       do i=fft_start(3),fft_end(3)
          if(i .EQ. n0z) cycle
          do j=fft_start(2),fft_end(2)
             if(j .EQ. n0y) cycle
             do h=fft_start(1),fft_end(1)
                if(h .EQ. n0x) cycle
                tmpi = 2._mytype   !tmpi permet de prendre en compte la somme totale a partir 
                                   !d'un espace (1:nx/2+1,1:ny,1:nz)
                if (h == 1) tmpi=1._mytype
                if (modulo(grid%nx,2) == 0 .AND. h == grid%nx / 2 + 1) tmpi = 1._mytype

                w=0
                MatDefR(1,1)= real(DefF(h,j,i,1),mytype)
                MatDefR(2,2)= real(DefF(h,j,i,2),mytype)
                MatDefR(3,3)= real(DefF(h,j,i,3),mytype)
                MatDefR(1,2)= real(DefF(h,j,i,4),mytype)
                MatDefR(1,3)= real(DefF(h,j,i,5),mytype)
                MatDefR(2,3)= real(DefF(h,j,i,6),mytype)
                MatDefR(2,1)= real(DefF(h,j,i,7),mytype)
                MatDefR(3,1)= real(DefF(h,j,i,8),mytype)
                MatDefR(3,2)= real(DefF(h,j,i,9),mytype)

                MatDefI(1,1)= real(aimag(DefF(h,j,i,1)),mytype)
                MatDefI(2,2)= real(aimag(DefF(h,j,i,2)),mytype)
                MatDefI(3,3)= real(aimag(DefF(h,j,i,3)),mytype)
                MatDefI(1,2)= real(aimag(DefF(h,j,i,4)),mytype)
                MatDefI(1,3)= real(aimag(DefF(h,j,i,5)),mytype)
                MatDefI(2,3)= real(aimag(DefF(h,j,i,6)),mytype)
                MatDefI(2,1)= real(aimag(DefF(h,j,i,7)),mytype)
                MatDefI(3,1)= real(aimag(DefF(h,j,i,8)),mytype)
                MatDefI(3,2)= real(aimag(DefF(h,j,i,9)),mytype)
                dR = 0
                dI = 0
                do r=1,3
                   do k=1,3
                      do l=1,3
                         do m=1,3
                            do p=1,3
                               do q=1,3
                                  prodMatLC = MatLeviC(r,p,l)*MatLeviC(k,q,m)
                                  if(prodMatLC .NE. 0) then
                                     dR(r,k)=dR(r,k)+prodMatLC*&
                                          MatDefR(l,m)*FREQ(h,j,i,p)*FREQ(h,j,i,q)
                                     dI(r,k)=dI(r,k)+prodMatLC*&
                                          MatDefI(l,m)*FREQ(h,j,i,p)*FREQ(h,j,i,q)
                                  end if
                               end do
                            end do
                         end do
                      end do
                      w = w+ dR(r,k)**2 + dI(r,k)**2
                   end do
                end do
                tmp = tmp + (w * tmpi) / grid%ntot
             end do
          end do
       end do
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif

       call MPI_ALLREDUCE(tmp,sum2Icomp,1,real_type,MPI_SUM,amitex_mpi_comm,ierr)
       sum2Icomp = sum2Icomp * grid%ntot * grid%ntot !--- Correction car defF normalisee (voir correction sum2Def)
       critCptb = sqrt(sum2Icomp/sum2Def)

       ! NORMALISATION du critere  
       ! => assurer l'independance au choix de dx,dy,dz
       ! => introduction possible d'une dimension "physique"
       critCptb = critCptb * grid%Tmax * grid%Tmax

    else
       !!========================================================================================
       !!                                                                  CRITERES EN CONTRAINTE
       !!
       !! On veut calculer ici 2 criteres differents :
       !! - criteq est le critere d'equilibre : norm(div(sigma).dV)_L^2/norm(sigma.n.dS)_L^2
       !! - critSigMoy est le critere sur la valeur moyenne de la contrainte
       !!========================================================================================
       !!--------------------------------------------------------------------------
       !!                                         SOMME DES CARRES DE LA CONTRAINTE 
       !!                                        sum2Sig =  \sum_{n=0}^{N} |x[n]|^2
       !!
       tmp9=0
       call field_sum_squareF(PiF,tmp9,grid%nx,grid%ny,grid%nz,9,1)
       sum2Sig = sum(tmp9)

       !!--------------------------------------------------------------------------
       !!                  CALCUL DU CRITERE SUR LA VALEUR MOYENNE DE LA CONTRAINTE
       !!
       !!                         critSigMoy = norme(Sigmoy) / sqrt(sum2Sig * dV/V)
       !!                                    = norme(Sigmoy) / sqrt(sum2Sig / ntot)
       !!
       !!    remarque : ce critere est calcule sur le pinceau contenant (1,1,1)
       !!               il faut redistribuer

       critPiMoy = 0
       normeDPi = 0
       tmp = 0
       !! Si aucune contrainte n'est imposee ce critere vaut 0      
       !! Cas d'un chargement par composante
       if (.not. DirStress_flag) then
       if((nb_stress_pilot>0) &
         .AND. (fft_start(1) .eq. 1) .AND. (fft_start(2) .eq. 1) .AND. (fft_start(3) .eq. 1)) then
          !! normeDPi est la norme L2 de l'ecart entre la contrainte moyenne
          !! et la contrainte imposee (sur les composantes concernees)
          do i=1,9
             if(local_loading(i) == 1) then
                normeDPi = normeDPi + (real(PiF(1,1,1,i),mytype)/grid%ntot - &
                     local_loading(i+9))**2
             end if
          end do
          tmp = sqrt(normeDPi/(sum2Sig/grid%ntot))
       end if
       !! Cas d'un chargement multiaxial en contrainte pilote en deformation
       elseif ((fft_start(1) .eq. 1) .AND. (fft_start(2) .eq. 1) .AND. (fft_start(3) .eq. 1)) then 
          alpha = 0.
          where (abs(DirStress) .gt. 1.e-6_mytype)
             alpha = real(PiF(1,1,1,:)/grid%ntot,mytype) / DirStress
          end where
          if (abs(maxval(alpha)) >= abs(minval(alpha))) then
             alpha_max = maxval(alpha)
          else
             alpha_max = minval(alpha)
          end if
          do i=1,9
                normeDPi = normeDPi + (real(PiF(1,1,1,i),mytype)/grid%ntot - &
                     alpha_max*DirStress(i))**2
          end do
          tmp = sqrt(normeDPi/(sum2Sig/grid%ntot))
       end if

       !! On recupere le max sur tous les processus pour que tous les processus
       !! aient la bonne valeur de critPiMoy
       call MPI_AllReduce(tmp,critPiMoy,1,real_type,MPI_MAX,amitex_mpi_comm,ierr)   
       
       
       !!--------------------------------------------------------------------------
       !!                      CALCUL DU CRITERE D'EQUILIBRE         
       !!                  
       !!            criteq = sqrt(somme des carres de la force feq)
       !!                   / sqrt(somme des carres des forces f1,f2,f3)
       !!
       !!            feq = div(pi)*dV
       !!            f1 = Pi.e1.dx2.dx3
       !!            f2 = Pi.e2.dx1.dx3
       !!            f3 = Pi.e3.dx1.dx2
       !!
       !!            ci-dessous : criteq = sqrt(sum2div)*dV / sqrt(sum2forces)
       !!
       criteq=0 !ici role de variable temporaire
#ifdef OPENMP
!$OMP PARALLEL NUM_THREADS(nbthread) private(tmp,fi)
!$OMP DO REDUCTION(+:criteq)
#endif
       do i=fft_start(3),fft_end(3)
          tmp=0
          do j=fft_start(2),fft_end(2)
             do h=fft_start(1),fft_end(1)
                ! fi = div(PiF)
                fi(1)=  (FREQ(h,j,i,1)* PiF(h,j,i,1) +FREQ(h,j,i,2)*PiF(h,j,i,4) + FREQ(h,j,i,3)* PiF(h,j,i,5) )
                fi(2)=  (FREQ(h,j,i,1)* PiF(h,j,i,7) +FREQ(h,j,i,2)*PiF(h,j,i,2) + FREQ(h,j,i,3)* PiF(h,j,i,6) )
                fi(3)=  (FREQ(h,j,i,1)* PiF(h,j,i,8) +FREQ(h,j,i,2)*PiF(h,j,i,9) + FREQ(h,j,i,3)* PiF(h,j,i,3) )

                ! prise en compte des forces volumiques nodales (ici on a i.fi, ce qui ne change pas la norme)
                if (associated(FvolNF)) fi = -fi + imP*FvolNF(h,j,i,:) * FREQ_2(h,j,i)          

                tmpi = 2._mytype
                if (h == 1) tmpi=1._mytype
                if (modulo(grid%nx,2) == 0 .AND. h == grid%nx / 2 + 1) tmpi = 1._mytype
                tmp=tmp+ (tmpi*real( fi(1)*conjg(fi(1))+fi(2)*conjg(fi(2))+fi(3)*conjg(fi(3))  ,mytype) &
                         / grid%ntot)
             end do
          end do
          criteq=criteq+tmp
       end do
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
       call MPI_AllReduce(criteq,sum2div,1,real_type,MPI_SUM,amitex_mpi_comm,ierr)

       ! DENOMINATEUR : somme des carres des forces Pi.e1.dx2.dx3, Pi.e2.dx1.dx3 et Pi.e3.dx1.dx2
       ! Ordre 11 22 33 12 13 23 21 31 32
       sum2forces = 0
       call field_sum_squareF(PiF(:,:,:,1),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dy*grid%dz)**2
       call field_sum_squareF(PiF(:,:,:,7),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dy*grid%dz)**2
       call field_sum_squareF(PiF(:,:,:,8),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dy*grid%dz)**2

       call field_sum_squareF(PiF(:,:,:,4),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dx*grid%dz)**2
       call field_sum_squareF(PiF(:,:,:,2),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dx*grid%dz)**2
       call field_sum_squareF(PiF(:,:,:,9),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dx*grid%dz)**2

       call field_sum_squareF(PiF(:,:,:,5),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dx*grid%dy)**2
       call field_sum_squareF(PiF(:,:,:,6),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dx*grid%dy)**2
       call field_sum_squareF(PiF(:,:,:,3),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dx*grid%dy)**2

       ! etape finale 
       if (sum2forces == 0.) sum2forces = 1.     ! evite la division par zero
       !criteq=sqrt(sum2div/sum2Sig) * grid%dx   !valable pour voxels cubiques
       criteq=sqrt(sum2div/sum2forces) *grid%dV  !voxels anisotropes (valide voxels composite ci-dessus)


    end if

  call MPI_BARRIER(amitex_mpi_comm,ierror)          
  times_g%crit = times_g%crit + MPI_WTIME() - t1

  end subroutine eval_criteres_GD
!====================================================================




!====================================================================
!> SUBROUTINE D'EVALUATION DU CRITERE D'EQUILIBRE ET 
!!            DU CRITERE SUR LES VALEURS MOYENNES EN DIFFUSION
!!
!!          format (pinceaux-Z, taille globale (nx/2+1)*ny*nz)
!!
!!
!! 
!! \param[in]      SigF: (complexe) champs des contraintes dans l'espace de Fourier
!! \param[in]      DefF: (complexe) champs des deformations dans l'espace de Fourier
!!                                  DeF est ici prealablement normalise par ntot
!! \param[in]      FREQ: (reel) champs des frequences
!! \param[in]      local_loading   chargement impose
!!                              -  local_loading(1:6) : pilotage
!!                              -  local_loading(7:12) : valeur du pilotage
!! \param[in]      fft_start, fft_end: (reel) bornes des tableaux SIGF et FREQ
!! \param[in]      nx,ny,nz:    Dimensions dans chaque direction
!! \param[in]      ntot: (entier) nombre de cellules dans le domaine si on veut calculer le critere sur la deformation
!!                                0 si on veut calculer les criteres sur la contrainte
!! \param[out]     criteq: (reel) critere d'equilibre
!! \param[out]     critSigMoy: (reel) ecart entre les contraintes moyennes obtenues 
!!                                 et les contraintes moyennes imposees
!! \param[out]     critDefMoy: (reel) ecart entre les deformations moyennes obtenues 
!!                                 et les deformations moyennes imposees
!! \param[out]     critCptb: (reel) critere de compatibilite
!!
!====================================================================
  subroutine eval_criteresD(FluxDF, GradQDF, FREQ, local_loadingD,&
       criteqD,critFluxDMoy,critGradQDMoy,critCptbD,ntot)

    implicit none

    complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3,algo_param%nVarD),&
         intent(in)  :: FluxDF
    complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3,algo_param%nVarD),&
         intent(in)  :: GradQDF
    real(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3),&
         intent(in)  :: FREQ
    real(mytype),dimension(6*algo_param%nVArD), intent(in)  :: local_loadingD
    !> nombre de voxels de la cellule ou 0 si on ne souhaite pas calculer le critere en compatbilite
    integer(kind=8), intent(in)             :: ntot
    complex(mytype)                         :: fi
    integer                                 :: i,j,h,ivar, ierr, nb_flux_pilot,n0x,n0y,n0z,nVarD
    integer(kind=8)                         :: cell_size
#ifdef OPENMP
    integer                                 :: nbthread
#endif
    real(mytype),intent(out),dimension(:)   :: criteqD, critGradQDMoy, critFluxDMoy,critCptbD
    real(mytype)                            :: tmp,normeDFluxD,normeDGradQD,dR
    real(mytype)                            :: sum2FluxD, sum2forces,sum2div,tmpi,&
                                               sum2gradQD,sum2Icomp
    real(mytype),dimension(3)               :: tmp3
    real(mytype),dimension(1)               :: tmp1

#ifdef DOUBLE_PREC
    real(mytype),parameter                     :: epsilon = 10._mytype**(-300)
#else
    real(mytype),parameter                     :: epsilon = 10._mytype**(-30)
#endif
    integer                       :: ierror
    double precision              :: t1
  
    call MPI_BARRIER(amitex_mpi_comm,ierror)          
    t1 = MPI_WTIME()
    
#ifdef OPENMP
!$OMP PARALLEL
  nbthread = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
#endif
   
    nVarD = algo_param%nVarD
    cell_size= grid%ntot

!---------- BOUCLE SUR LES VARIABLES DE DIFFUSION
do ivar=1,nVarD
    nb_flux_pilot = count(local_loadingD(1+3*(ivar-1):3*ivar) == 1)

    if(ntot>0) then
       !! Si on rentre dans cette fonction avec computeDef=1 on calcule
       !! - le critere sur la def moyenne
       !! - critCptb est le critere de compatibilite
       !!========================================================================================
       !!                                                                 CRITERES EN DEFORMATION
       !!
       !! On veut calculer ici 2 criteres differents :
       !! - le critere sur la def moyenne
       !! - critCptb est le critere de compatibilite
       !!
       !! ATTNETION : en entree DefF est divise par ntot 
       !!            -> pas besoin de diviser a nouveau pour obtenir la deformation moyenne 
       !!            -> corriger la somme des carres
       !!========================================================================================
       !!--------------------------------------------------------------------------
       !!                                        SOMME DES CARRES DE LA DEFORMATION 
       !!                                        sum2Def =  \sum_{n=0}^{N} |x[n]|^2
       !!
       !! ATTENTION : DefF est normalise par Ntot -> penser a corriger 
       !!
       !! Dans le cas ou on impose un chargement thermique avec deformation ou 
       !! contrainte moyenne nulle, la premiere valeur du champ de deformation proposee
       !! par l'algorithme est identiquement nulle. Ainsi, on ne peut pas normaliser
       !! par cette valeur => seuil

       tmp3=0
       call field_sum_squareF(GradQDF(:,:,:,:,ivar),tmp3,grid%nx,grid%ny,grid%nz,3,1)
       tmp3=tmp3 * grid%ntot * grid%ntot !--- correction car DefF est ici normalise par ntot
       sum2gradQD = sum(tmp3(1:3))
       if(sum2gradQD<epsilon) sum2gradQD = epsilon

       !!--------------------------------------------------------------------------
       !!                                        CRITERE SUR LA DEFORMATION MOYENNE
       !!                         critDefMoy = norme(Defmoy) / sqrt(sum2Def * dV/V)
       !!                                    = norme(Defmoy) / sqrt(sum2Def / ntot)
       !!      
       !!             Si aucune deformation moyenne n'est imposee le critere vaut 0
       critGradQDMoy(ivar) = 0
       if(nb_flux_pilot<3) then
          !! la deformation moyenne est imposee dans au moins une direction.
          !! On veut calculer l'ecart relatif entre la deformation moyenne
          !! et celle que l'on veut imposer (pour certaines directions)
          normeDGradQD = 0._mytype 
          tmp = 0._mytype
          if((fft_start(1) .eq. 1) .AND. (fft_start(2) .eq. 1) .AND. (fft_start(3) .eq. 1)) then
             !! normeDDef est la norme L2 de l'ecart entre la deformation moyenne
             !! et la deformation imposee (sur les composantes concernees)
             !! Prise en compte du fait que DefF est divise par ntot en simple precision
             do i=1+3*(ivar-1),3*ivar
                if(local_loadingD(i) == 0) then
                   normeDGradQD = normeDGradQD + (real(GradQDF(1,1,1,i,ivar),mytype) - &
                        local_loadingD(i+3*nVarD))**2
                end if
             end do
             !! On calcule tmp le critere sur les valeurs moyennes de la deformation.
             !! Le but est de verifier que les deformations moyennes correspondent bien a
             !! celles que l'on veut imposer
             tmp = sqrt(normeDGradQD/(sum2gradQD/grid%ntot))
             !! On normalise par la norme L2 de la deformation
          end if
          !! On recupere le max sur tous les processus pour que tous les processus
          !! aient la bonne valeur de critDefMoy
          call MPI_AllReduce(tmp,critGradQDMoy(ivar),1,real_type,MPI_MAX,amitex_mpi_comm,ierr)
       end if

       !!--------------------------------------------------------------------------
       !!                                        CALCUL DU CRITERE DE COMPATIBILITE         
       !!                  
       !!            critcptb = sqrt(somme des carres du tenseur d'incompatibilite )
       !!                   / sqrt(somme des carres de la deformation)
       !!   
       !!            crtcptb = cricptb*Tmax*Tmax
       !! 
       !! Faute d'argument tel que pour le critere d'equilibre, 
       !! on assure l'independance a la dimension en multipliant pat Tmax^2
       !! on aurait pu choisir dx^2. Le choix realise a le merite d'etre conservatif.  
       !!
       !! Traitement eventuel des cas pairs (a supprimer, plus utilise) :
       !!
       !! On elimine les frequences introduites dans le cas des resolutions paires
       !! En fait, on definit n0i qui vaut 
       !! - n0i = -1 si le nombre de voxels dans la direction i est impair
       !! - n0i = ni/2+1 si le nombre de voxels dans la direction i est pair
       !! Si on est sur un des indices n0i on ne fait pas le calcul
       !! dR=dR+LeviCivita(r,p,l)*LeviCivita(k,q,m)*Def(l,m)*FREQ(h,j,i,p)*FREQ(h,j,i,q)  
       !! ce qui est equivalent a prendre des frequences nulles sur ces points
       !! Ce cas n'arrive pas si n0i = -1
       n0x=-1
       n0y=-1
       n0z=-1      
       if (algo_param%CorrPairD) then 
       if (mod(grid%nx,2) .eq. 0) then
          n0x = grid%nx/2+1
       end if
       if (mod(grid%ny,2) .eq. 0) then
          n0y = grid%ny/2+1
       end if
       if (mod(grid%nz,2) .eq. 0) then
          n0z = grid%nz/2+1
       end if
       end if
       tmp=0

#ifdef OPENMP
!$OMP PARALLEL NUM_THREADS(nbthread) private(dR)
!$OMP DO REDUCTION(+:tmp)
#endif
       do i=fft_start(3),fft_end(3)
          if(i .EQ. n0z) cycle
          do j=fft_start(2),fft_end(2)
             if(j .EQ. n0y) cycle
             do h=fft_start(1),fft_end(1)
                if(h .EQ. n0x) cycle
                tmp=0

                tmpi = 2._mytype   !tmpi permet de prendre en compte la somme totale a partir d'un espace (1:nx/2+1,1:ny,1:nz)
                if (h == 1) tmpi=1._mytype
                if (modulo(grid%nx,2) == 0 .AND. h == grid%nx / 2 + 1) tmpi = 1._mytype

                dR =   real(GradQDF(h,j,i,2,ivar),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,2)&
                     - real(GradQDF(h,j,i,1,ivar),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,2)&
                     - real(GradQDF(h,j,i,1,ivar),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,3)&
                     + real(GradQDF(h,j,i,3,ivar),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,3)
                tmp = tmp + (tmpi*dR**2)/grid%ntot
                dR =   real(GradQDF(h,j,i,3,ivar),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,3)&
                     - real(GradQDF(h,j,i,2,ivar),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,3)&
                     - real(GradQDF(h,j,i,2,ivar),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,1)&
                     + real(GradQDF(h,j,i,1,ivar),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,1)
                tmp = tmp + (tmpi*dR**2)/grid%ntot
                dR =   real(GradQDF(h,j,i,1,ivar),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,1)&
                     - real(GradQDF(h,j,i,3,ivar),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,1)&
                     - real(GradQDF(h,j,i,3,ivar),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,2)&
                     + real(GradQDF(h,j,i,2,ivar),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,2)
                tmp = tmp + (tmpi*dR**2)/grid%ntot
                !!Parties imaginaires
                dR =   real(aimag(GradQDF(h,j,i,2,ivar)),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,2)&
                     - real(aimag(GradQDF(h,j,i,1,ivar)),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,2)&
                     - real(aimag(GradQDF(h,j,i,1,ivar)),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,3)&
                     + real(aimag(GradQDF(h,j,i,3,ivar)),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,3)
                tmp = tmp + (tmpi*dR**2)/grid%ntot
                dR =   real(aimag(GradQDF(h,j,i,3,ivar)),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,3)&
                     - real(aimag(GradQDF(h,j,i,2,ivar)),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,3)&
                     - real(aimag(GradQDF(h,j,i,2,ivar)),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,1)&
                     + real(aimag(GradQDF(h,j,i,1,ivar)),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,1)
                tmp = tmp + (tmpi*dR**2)/grid%ntot
                dR =   real(aimag(GradQDF(h,j,i,1,ivar)),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,1)&
                     - real(aimag(GradQDF(h,j,i,3,ivar)),mytype)*FREQ(h,j,i,1)*FREQ(h,j,i,1)&
                     - real(aimag(GradQDF(h,j,i,3,ivar)),mytype)*FREQ(h,j,i,2)*FREQ(h,j,i,2)&
                     + real(aimag(GradQDF(h,j,i,2,ivar)),mytype)*FREQ(h,j,i,3)*FREQ(h,j,i,2)
                tmp = tmp + (tmpi*dR**2)/grid%ntot
             end do
          end do
       end do
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif

      call MPI_ALLREDUCE(tmp,sum2Icomp,1,real_type,MPI_SUM,amitex_mpi_comm,ierr)
       sum2Icomp = sum2Icomp * grid%ntot * grid%ntot !--- Correction car defF normalisee (voir correction sum2Def)
       critCptbD(ivar) = sqrt(sum2Icomp/sum2gradQD)

       ! NORMALISATION du critere 
       ! => assurer l'independance au choix de dx,dy,dz
       ! => introduction possible d'une dimension "physique"
       critCptbD(ivar) = critCptbD(ivar) * grid%Tmax * grid%Tmax

    else
       !!========================================================================================
       !!                                                                  CRITERES EN CONTRAINTE
       !!
       !! On veut calculer ici 2 criteres differents :
       !! - criteq est le critere d'equilibre : norm(div(sigma).dV)_L^2/norm(sigma.n.dS)_L^2
       !! - critSigMoy est le critere sur la valeur moyenne de la contrainte
       !!========================================================================================
       !!--------------------------------------------------------------------------
       !!                                         SOMME DES CARRES DE LA CONTRAINTE 
       !!                                        sum2Sig =  \sum_{n=0}^{N} |x[n]|^2
       !!
       tmp3=0
       call field_sum_squareF(FluxDF(:,:,:,:,ivar),tmp3,grid%nx,grid%ny,grid%nz,3,1)
       sum2FluxD = sum(tmp3)

      !!--------------------------------------------------------------------------
       !!                  CALCUL DU CRITERE SUR LA VALEUR MOYENNE DE LA CONTRAINTE
       !!
       !!                         critSigMoy = norme(Sigmoy) / sqrt(sum2Sig * dV/V)
       !!                                    = norme(Sigmoy) / sqrt(sum2Sig / ntot)
       !!
       !!    remarque : ce critere est calcule sur le pinceau contenant (1,1,1)
       !!               il faut redistribuer

       critFluxDMoy(ivar) = 0
       normeDFluxD = 0
       tmp = 0
       !! Si aucune contrainte n'est imposee ce critere vaut 0
       if((nb_flux_pilot>0) .AND. (fft_start(1) .eq. 1) .AND. (fft_start(2) .eq. 1) .AND. (fft_start(3) .eq. 1)) then
          !! normeDSig est la norme 2 de l'ecart entre la contrainte moyenne
          !! et la contrainte imposee (sur les composantes concernees)
          do i=1+3*(ivar-1),3*ivar
             if(local_loadingD(i) == 1) then
                normeDFluxD = normeDFluxD + (real(FluxDF(1,1,1,i,ivar),mytype)/grid%ntot - &
                     local_loadingD(i+3*nVarD))**2
             end if
          end do

          !! On met dans tmp le critere sur les valeurs moyennes de la contrainte
          !! Le but est de verifier que les contraintes
          !! moyennes correspondent bien a celles que l'on veut imposer
          tmp = sqrt(normeDFluxD/(sum2FluxD/grid%ntot))
       end if
       !! On recupere le max sur tous les processus pour que tous les processus
       !! aient la bonne valeur de critSigMoy
       call MPI_AllReduce(tmp,critFluxDMoy(ivar),1,real_type,MPI_MAX,amitex_mpi_comm,ierr)   

       !!--------------------------------------------------------------------------
       !!                      CALCUL DU CRITERE D'EQUILIBRE         
       !!                  
       !!            criteq = sqrt(somme des carres de la force feq)
       !!                   / sqrt(somme des carres des forces f1,f2,f3)
       !!
       !!            feq = div(sig)*dV
       !!            f1 = Sig.e1.dx2.dx3
       !!            f2 = Sig.e2.dx1.dx3
       !!            f3 = Sig.e3.dx1.dx2
       !!
       !!            ci-dessous : criteq = sqrt(sum2div)*dV / sqrt(sum2forces)
       !!

       ! NUMERATEUR : somme sur les points des carres de div(sigma) (a multiplier par dV -> etape finale)
       criteqD(ivar)=0 !ici role de variable temporaire
#ifdef OPENMP
!$OMP PARALLEL NUM_THREADS(nbthread) private(tmp,fi)
!$OMP DO REDUCTION(+:criteqD)
#endif
       do i=fft_start(3),fft_end(3)
          tmp=0
          do j=fft_start(2),fft_end(2)
             do h=fft_start(1),fft_end(1)
                ! fi = div(SigF)
                fi=  FREQ(h,j,i,1)* FluxDF(h,j,i,1,ivar) +FREQ(h,j,i,2)*FluxDF(h,j,i,2,ivar) + FREQ(h,j,i,3)* FluxDF(h,j,i,3,ivar)
                tmpi = 2._mytype
                if (h == 1) tmpi=1._mytype
                if (modulo(grid%nx,2) == 0 .AND. h == grid%nx / 2 + 1) tmpi = 1._mytype
                tmp=tmp+tmpi * real( fi*conjg(fi), mytype) / grid%ntot
             end do
          end do
          criteqD(ivar)=criteqD(ivar)+tmp
       end do
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
       call MPI_AllReduce(criteqD(ivar),sum2div,1,real_type,MPI_SUM,amitex_mpi_comm,ierr)

       ! DENOMINATEUR : somme des carres des forces Sig.e1.dx2.dx3, Sig.e2.dx1.dx3 et Sig.e3.dx1.dx2
       sum2forces = 0._mytype
       tmp1 = 0._mytype
       call field_sum_squareF(FluxDF(:,:,:,1,ivar),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dy*grid%dz)**2

       call field_sum_squareF(FluxDF(:,:,:,2,ivar),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dx*grid%dz)**2

       call field_sum_squareF(FluxDF(:,:,:,3,ivar),tmp1,grid%nx,grid%ny,grid%nz,1,1)
       sum2forces = sum2forces + tmp1(1)*(grid%dx*grid%dy)**2

       ! etape finale 
       if (sum2forces == 0.) sum2forces = 1.     ! evite la division par zero
       !criteq=sqrt(sum2div/sum2Sig) * grid%dx   !valable pour voxels cubiques
       criteqD(ivar)=sqrt(sum2div/sum2forces) *grid%dV  !voxels anisotropes (valide voxels composite ci-dessus)

    end if

end do
!---------- BOUCLE SUR LES VARIABLES DE DIFFUSION

  call MPI_BARRIER(amitex_mpi_comm,ierror)          
  times_g%crit = times_g%crit + MPI_WTIME() - t1

  end subroutine eval_criteresD
!====================================================================



!===================================================================================================
!> Subroutine d'application du tenseur de Green
!! sans passer par le stockage de GREEN0
!! 
!! Notation CAST3M (Voigt+ordre 11 22 33 12 13 23)
!!
!!   \param[in] SigF   contrainte (notation CAST3M), espace Fourier 
!!                   pinceaux-Z (nx/2+1,ny,nz,6)
!!   \param[in] FREQ  tableau de fréquence, espace Fourier 
!!                   pinceaux-Z (nx/2+1,ny,nz,3)
!!   \param[in] LambdaMu0   Coefficients [Lambda0,mu0] du milieu de référence
!!   \param[in] fft_start, fft_end    limites des tableaux ou seront stockees
!!                                    les FFT   
!!   \param[out] SigF vaut - gamma0*SigF en sortie
!!                    avec les notations CAST3M pour la deformation
!!
!===================================================================================================
  subroutine apply_greenD(FluxDF, FREQ, K0D)

    implicit none

    real(mytype), dimension(algo_param%nVarD), intent(in) :: K0D
    complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),3,algo_param%nVarD), &
         intent(inout)                                     :: FluxDF
    real(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3), &
         intent(in)                                        :: FREQ
    ! variables temporaires pour la correction Paire (desuet)
    complex(mytype), allocatable, dimension(:,:,:,:) :: FluxDFx,FluxDFy,FluxDFz !gcc-warning accepte 
                                                                                ! (allocate dans if)
                                              
    integer                    :: nx2,ny2,nz2
    real(mytype)    :: norme2                            ! variable temporaire (norme2 de la frequence)
    complex(mytype), dimension(algo_param%nVarD) :: Q,fi ! variables temporaires (quantite et divergence)
    complex(mytype) :: imP = cmplx(0,1,kind=mytype) ! le nombre complexe i
    integer :: i,j,h,ivar,nVarD ! variables temporaires
#ifdef OPENMP
    integer  :: nbthread 
#endif
    integer                       :: ierror
    double precision              :: t1
  
    call MPI_BARRIER(amitex_mpi_comm,ierror)          
    t1 = MPI_WTIME()

#ifdef OPENMP
!$OMP PARALLEL
  nbthread = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
#endif

    nVarD = algo_param%nVarD

    ! On stocke dans les SigFi les valeurs de SigF 
    ! reutilisees dans le calcul pour le cas pair
    if (algo_param%CorrPairD) then
    if(modulo(grid%nx,2) == 0) then
       nx2 = grid%nx/2 +1
       if(nx2>=fft_start(1) .AND. nx2<=fft_end(1)) then
          allocate(FluxDFx(fft_start(2):fft_end(2),fft_start(3):fft_end(3),3,nVarD))
          do ivar=1,nVarD
          do i=1,3
             do j=fft_start(3),fft_end(3)
                do h=fft_start(2),fft_end(2)
                   FluxDFx(h,j,i,ivar) = FluxDF(nx2,h,j,i,ivar)
                end do
             end do
          end do
          end do
       end if
    end if
    if(modulo(grid%ny,2) == 0) then
       ny2 = grid%ny/2 +1
       if(ny2>=fft_start(2) .AND. ny2<=fft_end(2)) then
          allocate(FluxDFy(fft_start(1):fft_end(1),fft_start(3):fft_end(3),3,nVarD))
          do ivar=1,nVarD
          do i=1,3
             do j=fft_start(3),fft_end(3)
                do h=fft_start(1),fft_end(1)
                   FluxDFy(h,j,i,ivar) = FluxDF(h,ny2,j,i,ivar)
                end do
             end do
          end do
          end do
       end if
    end if
    if(modulo(grid%nz,2) == 0) then
       nz2 = grid%nz/2 +1
       if(nz2>=fft_start(3) .AND. nz2<=fft_end(3)) then
          allocate(FluxDFz(fft_start(1):fft_end(1),fft_start(2):fft_end(2),3,nVarD))
          do ivar=1,nVarD
          do i=1,3
             do j=fft_start(2),fft_end(2)
                do h=fft_start(1),fft_end(1)
                   FluxDFz(h,j,i,ivar) = FluxDF(h,j,nz2,i,ivar)
                end do
             end do
          end do
          end do
       end if
    end if
    end if
!!!!!!!!!!!!!!!!!!

#ifdef OPENMP
!$OMP PARALLEL NUM_THREADS(nbthread) private(fi,norme2,Q)
!$OMP DO
#endif
    do ivar =1,nVarD
    do i=fft_start(3),fft_end(3)
       do j=fft_start(2),fft_end(2)
          do h=fft_start(1),fft_end(1)
             ! fi = div(SigF)
             fi(ivar)= imP * (  FREQ(h,j,i,1)* FluxDF(h,j,i,1,ivar) &
                              + FREQ(h,j,i,2)* FluxDF(h,j,i,2,ivar) &
                              + FREQ(h,j,i,3)* FluxDF(h,j,i,3,ivar))

             norme2 = FREQ(h,j,i,1) **2 +  FREQ(h,j,i,2) **2 +  FREQ(h,j,i,3) **2
             if ((h .eq. 1) .AND. (i .eq. 1) .AND. (j .eq. 1)) then
                norme2=1
             end if
             !! Calcul de Q = -(i.tau.xsi)/(k0.xsi2)
             Q(ivar) = -fi(ivar) / (K0D(ivar)*norme2)
             
             !! calcul du gradent de Q mis dans FluxDF 
             FluxDF(h,j,i,1,ivar) = imP * FREQ(h,j,i,1)*Q(ivar)
             FluxDF(h,j,i,2,ivar) = imP * FREQ(h,j,i,2)*Q(ivar)
             FluxDF(h,j,i,3,ivar) = imP * FREQ(h,j,i,3)*Q(ivar)
          end do
       end do
    end do
    end do
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif


!!! Traitement des cas pairs
    if (algo_param%CorrPairD) then
    ! On teste la parite sur chaque composante
    if(modulo(grid%nx,2) == 0) then
       nx2= grid%nx/2+1
       if(nx2>=fft_start(1) .AND. nx2<=fft_end(1)) then
          !! Traitement reserve au cas pair (voir Moulinec Suquet '98)
          do ivar=1,nVarD
          do i=1,3
          do j=fft_start(3),fft_end(3)
             do h=fft_start(2),fft_end(2)
!                SigF(nx2,h,j,6) = -SigFx(h,j,6)*imu0
                FluxDF(nx2,h,j,i,ivar) = FluxDFx(h,j,i,ivar) / K0D(ivar)
             end do
          end do
          end do
          end do
          deallocate(FluxDFx)
       end if
    end if
    !! Idem pour la direction y
    if(modulo(grid%ny,2) == 0) then
       ny2= grid%ny/2+1
       if(ny2>=fft_start(2) .AND. ny2<=fft_end(2)) then
          do ivar=1,nVarD
          do i=1,3
          do j=fft_start(3),fft_end(3)
             do h=fft_start(1),fft_end(1)
!                SigF(h,ny2,j,6) = -SigFy(h,j,6)*imu0
                FluxDF(h,ny2,j,i,ivar) = FluxDFy(h,j,i,ivar) / K0D(ivar)
             end do
          end do
          end do
          end do
          deallocate(FluxDFy)
       end if
    end if
    !! Idem pour la direction z
    if(modulo(grid%nz,2) == 0) then
       nz2= grid%nz/2+1
       if(nz2>=fft_start(3) .AND. nz2<=fft_end(3)) then
          do ivar=1,nVarD
          do i=1,3
          do j=fft_start(2),fft_end(2)
             do h=fft_start(1),fft_end(1)
!                SigF(h,j,nz2,6) = -SigFz(h,j,6)*imu0
                FluxDF(h,j,nz2,i,ivar) = FluxDFz(h,j,i,ivar) / K0D(ivar)
             end do
          end do
          end do
          end do
          deallocate(FluxDFz)
       end if
    end if
    end if ! fin test CorrPairD

  call MPI_BARRIER(amitex_mpi_comm,ierror)          
  times_g%apply = times_g%apply + MPI_WTIME() - t1

  end subroutine apply_greenD
!==================================================================


end module green_mod
