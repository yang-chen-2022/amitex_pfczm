program test_eval_criteres

!> CHARGEMENT DES DIFFERENTS MODULES
  use MPI             
#ifdef NEW_2DECOMP
  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
#else
  use decomp_2d  
#endif  
  use decomp_2d_fft

  use green_mod
  use error_mod
  use field_mod
  use amitex_lib
  use simu_mod
  
!!==============================================================================
!!						      		    DECLARATIONS
!!==============================================================================
  implicit none
!!------------------------------------------------------------------------------
!> 							 DIMENSIONS DE LA GRILLE

!  integer,parameter               :: nx=200, ny=50, nz=100         !< dimensions de la cellule
!  real(mytype),parameter          :: Lx=200, Ly=50, Lz=100         !< dimensions de la cellule
  integer,parameter               :: nx=200, ny=200, nz=200         !< dimensions de la cellule
  real(mytype),parameter          :: Lx=10, Ly=10, Lz=10         !< dimensions de la cellule
  real(mytype)                    :: dx, dy, dz,x0=0.,y0=0.,z0=0.     
               

!!------------------------------------------------------------------------------
!> 						   		    PARALLELISME 

  integer :: p_row=2, p_col=2           !< nb. de lignes, de colonnes pour la decomposition de 2decomp
  integer :: ierror                     !< erreur relative au fonction MPI

!!------------------------------------------------------------------------------
!>                                  		     CONTRAINTES ET DEFORMATIONS 
!>					    		tableaux paralleles

  !> dans l'espace reel : tableaux (ntot,Ntens)
  real(mytype),allocatable,dimension(:,:,:,:)    :: Sig,Def !< Contrainte de Cauchy (ntot,6)
  real(mytype),allocatable,dimension(:,:,:,:)    :: divSig !< Contrainte de Cauchy (ntot,6)
  real(mytype),allocatable,dimension(:,:,:,:)    :: PK1     


  !> dans l'espace de Fourier : tableaux (nx,ny,nz,Ntens), ATTENTION AU SENS DE SigF
  complex(mytype),allocatable,dimension(:,:,:,:)        :: SigF,DefF    !< Contrainte de Cauchy 
  complex(mytype),allocatable,dimension(:,:,:,:)        :: PK1F    !< PK1 (sinon) 


!!------------------------------------------------------------------------------
!>                                  		     COORDONNES GRILLE + divergence

  real(mytype),allocatable,dimension(:,:,:)        :: X,Y,Z


!!------------------------------------------------------------------------------
!> 							 		  DIVERS

  integer :: i,j,k                      !< indice de boucle
  integer :: alloc_stat                 !< erreur lors d'une allocation memoire

  real(mytype),parameter        :: PI = 4._mytype*DATAN(1._mytype), sqrt3 = sqrt(3._mytype)
  real(mytype)                  :: tmp, ttt, crit_eq,critSigMoy,critDefMoy,critCptb,sum2sig,&
                                   critSigMoy0,crit_eq0,normDsig,sum2div
  real(mytype),dimension(12)    :: loading
  real(mytype),dimension(18)    :: loadingGD
  real(mytype),dimension(9)     :: dzeros9
  real(mytype),dimension(6)     :: dzeros6,dsigmoy,somme
  integer(kind=8)               :: i0

!!==============================================================================
!!						  ALLOCATIONS ET INITIALISATIONS
!!
!!              ATTENTION : L'ORDRE DES INITIALISATIONS CI-DESSOUS EST IMPORTANT
!!
!!==============================================================================

!!------------------------------------------------------------------------------
!>                               			      INITIALISATION MPI 

  call MPI_INIT(ierror)
  !t2 = MPI_WTIME()

!!------------------------------------------------------------------------------
!>                               			 PRE-INITIALISATION AMITEX_LIB 

  allocate(SIMU(1))
  call init_amitex_lib_0()

!!------------------------------------------------------------------------------
!>                          		     INITIALISATION 2DECOMP, 2DECOMP_FFT 

  call decomp_2d_init(nx,ny,nz,p_row,p_col)     !initialisation decomp_2d
  call decomp_2d_fft_init(PHYSICAL_IN_X, nx, ny, nz)                      !initialisation decomp_2d_fft
  call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)

!!------------------------------------------------------------------------------
!>                          		          FIN INITIALISATION AMITEX_LIB 

  dx = Lx/nx
  dy = Ly/ny
  dz = Lz/nz
  call init_amitex_lib(nx,ny,nz,dx,dy,dz,x0,y0,z0,ficlog="test_criteres.log")!,filter_type="no_filter",filter_radius=0._8)

  call mpi_barrier(mpi_comm_world,ierror)
  !t1 = MPI_WTIME() - t2

  if(nrank==0)print *, "Dimensions : ",nx,"x",ny,"x",nz

!!------------------------------------------------------------------------------
!>	       INITIALISATION DES CHAMPS DE FREQUENCE, CONTRAINTE ET DEFORMATION 
!>                                                             POUR LA MECANIQUE

  ! Contrainte, deformation
  allocate(Sig(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),1:6),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex13)",2)

  allocate(Def(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),1:6),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex13)",2)

  allocate(PK1(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),1:9),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex13)",2)
     
  ! Contrainte, deformation, espace spectral (pinceaux-Z, taille globale ~(nx/2+1)*ny*nz)
  allocate(SigF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),&
       1:6),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex18)",2)

  allocate(DefF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),&
       1:6),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex18)",2)

  allocate(PK1F(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),&
       1:9),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex18)",2)

  ! Coordonnees grille
  allocate(X(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex13)",2)
  allocate(Y(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex13)",2)
  allocate(Z(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex13)",2)

  allocate(divSig(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),3),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex13)",2)
  
  ! Initialisations a 0
  sig =0.
  SigF=0.
  DefF=0
  X=0.;Y=0.;Z=0.

!!==============================================================================
!!					      INITIALISATION CHAMP DE CONTRAINTE
!!                                                            ON IMPOS PK1 = SIG
!!                                                  CALCUL ANALYTIQUE DIVERGENCE
!!==============================================================================
  do k=xstart(3),xend(3)
  do j=xstart(2),xend(2)
  do i=xstart(1),xend(1)
      X(i,j,k) = (i - 0.5_mytype) * dx
      Y(i,j,k) = (j - 0.5_mytype) * dy
      Z(i,j,k) = (k - 0.5_mytype) * dz
      Sig(i,j,k,1) = 1._mytype 
      Sig(i,j,k,2) = 2._mytype
      Sig(i,j,k,3) = 3._mytype
      Sig(i,j,k,4) = 4._mytype
      Sig(i,j,k,5) = 5._mytype
      Sig(i,j,k,6) = 6._mytype
  end do
  end do
  end do
  call mpi_barrier(mpi_comm_world,ierror)

  Sig(:,:,:,1) = Sig(:,:,:,1) + sin(10_mytype*pi*X/(nx*dx))
  Sig(:,:,:,2) = Sig(:,:,:,2) + cos(10._mytype*pi*Y/(ny*dy))
  Sig(:,:,:,3) = Sig(:,:,:,3) + sin(124._mytype + 10._mytype*pi*Z/(nz*dZ))
  Sig(:,:,:,4) = Sig(:,:,:,4) + sin(10._mytype*pi*X/(nx*dx)) + sin(10._mytype*pi*Y/(ny*dy)) 
  Sig(:,:,:,5) = Sig(:,:,:,5) + sin(10._mytype*pi*X/(nx*dx)) + sin(10._mytype*pi*Z/(nz*dz)) 
  Sig(:,:,:,6) = Sig(:,:,:,6) + sin(10._mytype*pi*Z/(nz*dz)) + sin(10._mytype*pi*Y/(ny*dy)) 

  PK1(:,:,:,1) = Sig(:,:,:,1)
  PK1(:,:,:,2) = Sig(:,:,:,2)
  PK1(:,:,:,3) = Sig(:,:,:,3)
  PK1(:,:,:,4) = Sig(:,:,:,4)
  PK1(:,:,:,5) = Sig(:,:,:,5)
  PK1(:,:,:,6) = Sig(:,:,:,6)
  PK1(:,:,:,7) = Sig(:,:,:,4)
  PK1(:,:,:,8) = Sig(:,:,:,5)
  PK1(:,:,:,9) = Sig(:,:,:,6)



!! CALCUL ANALYTIQUE DE DIV 
  divSig=0
  divSig(:,:,:,1) = (10._mytype*pi/(dx*nx)) * cos (10._mytype*pi*X / (nx*dx)) &
                  + (10_mytype*pi/(dy*ny)) * cos(10._mytype*pi*Y/(ny*dy)) &
                  + (10_mytype*pi/(dz*nz)) * cos(10._mytype*pi*Z/(nz*dz)) 

  divSig(:,:,:,2) =-(10_mytype*pi/(dy*ny)) * sin (10._mytype*pi*Y / (ny*dy)) &
                  + (10_mytype*pi/(dx*nx)) * cos(10._mytype*pi*X/(nx*dx)) &
                  + (10_mytype*pi/(dz*nz)) * cos(10._mytype*pi*Z/(nz*dz)) 

  divSig(:,:,:,3) =(10_mytype*pi/(dz*nz)) * cos (124._mytype + 10_mytype*pi*Z / (nz*dz)) &
                  + (10_mytype*pi/(dx*nx)) * cos(10._mytype*pi*X/(nx*dx)) &
                  + (10_mytype*pi/(dy*ny)) * cos(10._mytype*pi*Y/(ny*dy)) 


!!==============================================================================
!!					    TRANSFO FOURIER + VERIF DES MOYENNES
!!                                                       
!!==============================================================================

  call field_fft(Sig,SigF,6,1)
  call field_fft(PK1,PK1F,9,1)

  if (nrank==0) print*, "======================================================"
  if (nrank==0) print*, "                                    VERIF DES MOYENNES"
  if (nrank==0) print*, "======================================================"
 
  ! ESPACE REEL
  do i=1,6
  tmp = sum(Sig(:,:,:,i))
  call MPI_AllReduce(tmp,somme(i),1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
  if (nrank==0) print*, " somme(Sig(:,:,:,1) / ntot = ", somme(i) / real(nx*ny*nz,mytype)
  end do

  ! ESPACE FOURIER
  do i = 1,6
  if (nrank==0) print*, " SigF(1,1,1,:) / ntot = ", real(SigF(1,1,1,i) / real(nx*ny*nz,mytype))
  end do

!!==============================================================================
!!					                        CHARGEMENT MACRO
!!                                                            CONTRAINTE IMPOSEE
!!==============================================================================
   loading=[1,1,1,1,1,1,0,2,3,4,5,6]
   loadingGD=[1,1,1,1,1,1,1,1,1,0,2,3,4,5,6,4,5,6]

!!==============================================================================
!!					       EVALUATION DU CRITERE ESPACE REEL
!!                                                       
!!==============================================================================
!! CRITERE MACRO
  ttt = 1._mytype
  tmp=0
  normDsig=0.
  do i = 1,6
     if (i>3) ttt=2_mytype
     tmp = tmp + sum(ttt*Sig(:,:,:,i)*Sig(:,:,:,i))
     normDsig = normDsig + (somme(i)/(nx*ny*nz)-loading(6+i))*(somme(i)/(nx*ny*nz)-loading(6+i)) * ttt
  end do

  call MPI_AllReduce(tmp,sum2Sig,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
 
  normDsig = sqrt(normDsig)
  critSigmoy0 = normDsig / sqrt( sum2Sig/(nx*ny*nz) )



!! CALCUL ANALYTIQUE DE DIV 

  tmp = sum(divSig(:,:,:,1)*divSig(:,:,:,1) + divSig(:,:,:,2)*divSig(:,:,:,2) + divSig(:,:,:,3)*divSig(:,:,:,3))
  call MPI_AllReduce(tmp,sum2div,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)



  crit_eq0 = dx*sqrt(sum2div / sum2Sig)
  if (nrank==0) print*, "======================================================"
  if (nrank==0) print*, "EVALUATION DES CRITERES DEPUIS L'ESPACE REEL (CAS HPP)"
  if (nrank==0) print*, "======================================================"
  if (nrank==0) print*, " crit_eq0 (analytique) = ", crit_eq0 
  if (nrank==0) print*, " critSigMoy0 = ", critSigMoy0


!!==============================================================================
!!					    EVALUATION DU CRITERE ESPACE FOURIER
!!                                                       
!!==============================================================================

   dzeros6=0_mytype
   i0 = 0   ! on calcule les criteres d'equilibre
   !i0 = 1  ! on calcule les criteres en deformation
   call eval_criteres(SigF,DefF,FREQ,loading,&
                   .false.,dzeros6,& 
                   crit_eq,critSigMoy,critDefMoy,critCptb,i0)

  if (nrank==0) print*, "======================================================"
  if (nrank==0) print*, "                     SORTIE DE EVAL_CRITERES (CAS HPP)"
  if (nrank==0) print*, "======================================================"

  if (nrank==0) print*, " crit_eq = ", crit_eq 
  if (nrank==0) print*, " critSigMoy = ", critSigMoy


   call eval_criteres_GD(PK1F,DefF,FREQ,loadingGD,&
                   .false.,dzeros9,& 
                   crit_eq,critSigMoy,critDefMoy,critCptb,i0)

  if (nrank==0) print*, "======================================================"
  if (nrank==0) print*, "                  SORTIE DE EVAL_CRITERES_GD (CAS GD) "
  if (nrank==0) print*, "======================================================"

  if (nrank==0) print*, " crit_eq = ", crit_eq 
  if (nrank==0) print*, " critSigMoy = ", critSigMoy


!!==============================================================================
!!					      DESALLOCATIONS ET FIN DU PROGRAMME
!!==============================================================================
  call decomp_2d_fft_finalize
  call decomp_2d_finalize 
  call MPI_FINALIZE(ierror)

end program test_eval_criteres

