program test_field_filterF
!> CHARGEMENT DES DIFFERENTS MODULES

  use ISO_FORTRAN_ENV

  use MPI
#ifdef NEW_2DECOMP
  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
#else
  use decomp_2d
#endif  
  use decomp_2d_fft
  use io2_amitex_mod
  use amitex_mod
  use error_mod
  use green_mod
  use field_mod
  use amitex_lib
  use simu_mod


!!==============================================================================
!!                                                                  DECLARATIONS
!!==============================================================================
  implicit none
!!------------------------------------------------------------------------------
!>                                                       DIMENSIONS DE LA GRILLE   : ATTENTION, LECTURE D'UN VTK DE 65x65x65
!>                                                                       => ON DOIT UTILISER nx=ny=nz=65

  integer,parameter               :: nx=64, ny=64, nz=64         !< dimensions de la cellule
  real(mytype),parameter          :: Lx=10, Ly=10, Lz=10            !< dimensions de la cellule
  real(mytype)                    :: dx, dy, dz,x0=0.,y0=0.,z0=0.     


!!------------------------------------------------------------------------------
!>                                                                  PARALLELISME

  integer                      :: p_row=1, p_col=1           !< nb. de lignes, de colonnes pour la decomposition de 2decomp
  logical,dimension(3)         :: periodic_bc=.true.
  integer                      :: ierror                     !< erreur relative au fonction MPI

!!------------------------------------------------------------------------------
!>                                                           tableaux paralleles

  !> dans l'espace reel : tableaux (ntot)
  real(mytype),allocatable,dimension(:)                    :: V0, V1, Vm1, V2 
  real(mytype),allocatable,dimension(:,:,:)                :: V0c ! copy 3D and halo of V0
  real(mytype),allocatable,dimension(:,:,:)                :: V1r, Vm1r, V2r, X,Y,Z
  
  complex(mytype),allocatable,dimension(:,:,:)             :: VF
  real(mytype)                                             :: max_err,sumV,sigma
  integer                                                  :: Nsigma
  double precision                                         :: t1
  real(mytype),parameter                                   :: PI = 4._mytype*DATAN(1._mytype)
  
  
!!------------------------------------------------------------------------------
!>                                                           divers

  integer :: alloc_stat                 !< erreur lors d'une allocation memoire

!!------------------------------------------------------------------------------
!>                                                            INITIALISATION MPI

  call MPI_INIT(ierror)
  
!!------------------------------------------------------------------------------
!>                               			 PRE-INITIALISATION AMITEX_LIB 

  allocate(SIMU(1))
  call init_amitex_lib_0()

!!------------------------------------------------------------------------------
!>                                           INITIALISATION 2DECOMP, 2DECOMP_FFT

  call decomp_2d_init(nx,ny,nz,p_row,p_col,periodic_bc)     !initialisation decomp_2d
  call decomp_2d_fft_init(PHYSICAL_IN_X, nx, ny, nz)                      !initialisation decomp_2d_fft
  call decomp_2d_fft_get_size(fft_start,fft_end,fft_size) 

!!------------------------------------------------------------------------------
!>                          		          FIN INITIALISATION AMITEX_LIB 

  dx = Lx/nx
  dy = Ly/ny
  dz = Lz/nz
  call init_amitex_lib(nx,ny,nz,dx,dy,dz,x0,y0,z0,ficlog="test_criteres.log")!,filter_type="no_filter",filter_radius=0._8)

  call mpi_barrier(mpi_comm_world,ierror)
  
  
!!------------------------------------------------------------------------------
!>                                        INITIALISATION DE TABLEAUX TEMPORAIRES
!>

  allocate(TEMPfield32(xsize(1)*xsize(2)*xsize(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort(&
                     "espace memoire disponible insuffisant (amitex_fftp TEMPfield32)",2)
  TEMPfield32=0.


!!------------------------------------------------------------------------------
!>                                            ALLOCATION - INITIALISATION CHAMPS
!>                                                             

  allocate(V0(xsize(1)*xsize(2)*xsize(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(V1(xsize(1)*xsize(2)*xsize(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(Vm1(xsize(1)*xsize(2)*xsize(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(V2(xsize(1)*xsize(2)*xsize(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)

  allocate(V0c(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(V1r(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(Vm1r(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(V2r(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)

  allocate(X(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(Y(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(Z(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)

  allocate(VF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)

! definition des coordonnes de la grille (pout alidation filtre Gauss)
  call gene_cartesian_coord(X,Y,Z,dx,dy,dz)
  X = X - dx / 2._mytype
  Y = Y - dx / 2._mytype
  Z = Z - dz / 2._mytype

! On assigne V0 a la valeur du champ  
  !V0  = real(nrank,mytype)
  !V0c = real(nrank,mytype)
! Random value
  block
    integer :: size0
    integer, allocatable, dimension(:) :: seed
    call random_seed(SIZE=size0)
    allocate(seed(size0))
    seed=1
    call random_seed(PUT=seed)
    call random_number(V0)
    call random_seed(PUT=seed)
    call random_number(V0c)
  end block  
 
  V1 = 0.
  Vm1 = 0.
  V2 = 0.
  VF = 0.

!!------------------------------------------------------------------------------
!>                                          APPLICATION DES FILTRES DANS FOURIER
!>                                                             

   t1 = MPI_WTIME()
      
   ! Filtre 8 voisin Noeuds -> Centres
   call field_fft(V0c,VF,1,1)
   call field_filterF(VF,FILTER=1,ncomp=1,nvar=1,&
                   N=(/grid%nx,grid%ny,grid%nz /),&
                   d=(/grid%dx,grid%dy,grid%dz /))
   call field_ifft(V1,VF,1,1)
   V1 = V1 / grid%ntot
!   call print_field_vtk(V1,"V1.vtk","V1") 
      
   ! Filtre 8 voisin Centres -> Noeuds
   call field_fft(V0,VF,1,1)
   call field_filterF(VF,FILTER=-1,ncomp=1,nvar=1,&
                   N=(/grid%nx,grid%ny,grid%nz /),&
                   d=(/grid%dx,grid%dy,grid%dz /))
   call field_ifft(Vm1,VF,1,1)
   Vm1 = Vm1 / grid%ntot
!   call print_field_vtk(Vm1,"Vm1.vtk","Vm1") 
 
   ! Filtre 8 voisin "double" Centres -> Noeuds -> Centres
   call field_fft(V0,VF,1,1)
   call field_filterF(VF,FILTER=2,ncomp=1,nvar=1,&
                   N=(/grid%nx,grid%ny,grid%nz /),&
                   d=(/grid%dx,grid%dy,grid%dz /))
   call field_ifft(V2,VF,1,1)
   V2 = V2 / grid%ntot
!   call print_field_vtk(V2,"V2.vtk","V2")
 
   if(nrank==0) print*, "Times for Fourier Filters (s) = ",MPI_WTIME() - t1
   
!!------------------------------------------------------------------------------
!>                                      APPLICATION DES FILTRES DANS ESPACE REEL
!>                                                             

   t1 = MPI_WTIME()
               
  !-----------------------------------Filtre 8 voisin Noeuds -> Centres

  V1r = V0c
  call field_filter(V1r,FILTER=1,ncomp=1,nvar=1)
  !  call print_field_vtk(reshape(V1r, (/ xsize(1)*xsize(2)*xsize(3) /) ),"V1r.vtk","V1r") 

  !-----------------------------------Filtre 8 voisin Centres -> Noeuds

  Vm1r = V0c
  call field_filter(Vm1r,FILTER=-1,ncomp=1,nvar=1)
  !  call print_field_vtk(reshape(Vm1r, (/ xsize(1)*xsize(2)*xsize(3) /) ),"Vm1r.vtk","Vm1r")   

  !-----------------------------------Filtre 8 applique deux fois  
   
  V2r = V0c
  call field_filter(V2r,FILTER=1,ncomp=1,nvar=1)
  call field_filter(V2r,FILTER=-1,ncomp=1,nvar=1)
  !  call print_field_vtk(reshape(V2r, (/ xsize(1)*xsize(2)*xsize(3) /) ),"V2r.vtk","V2r")  
  
  if(nrank==0) print*, "Times for Real Filters (s) = ",MPI_WTIME() - t1

!!------------------------------------------------------------------------------
!>                                           VALIDATION : CALCUL DES DIFFERENCES
!>                                                             

  
  max_err = huge(max_err)
  call diff_fields(V1,V1r,max_err,typediff=1)
  if (nrank==0) print * , "max(abs(V1r - V1)) = ", max_err

  max_err = huge(max_err)
  call diff_fields(Vm1,Vm1r,max_err,typediff=1)
  if (nrank==0) print * , "max(abs(Vm1r - Vm1)) = ", max_err

  max_err = huge(max_err)
  call diff_fields(V2,V2r,max_err,typediff=1)
  if (nrank==0) print * , "max(abs(V2r - V2)) = ", max_err


!!------------------------------------------------------------------------------
!>                                                          TEST FILTRE GAUSSIEN
!>                                                             
   Nsigma = 5
   sigma = Nsigma * grid%dx

   !----------Initialisation : Dirac en (1,1,1)
   V0 = 0_mytype
   if (xstart(1) == 1 .and. xstart(2) == 1 .and. xstart(3) == 1) V0(1) = 1.
   
   !----------Application filtre Gaussien dans l'ESPACE DE FOURIER   
   call field_fft(V0,VF,1,1)
   call field_filterF(VF,FILTER=100,ncomp=1,nvar=1,&
                   N=(/grid%nx,grid%ny,grid%nz /),&
                   d=(/grid%dx,grid%dy,grid%dz /),param=(/ sigma /))
   call field_ifft(V1,VF,1,1)
   V1 = V1 / grid%ntot
   !call print_field_vtk(reshape(V1, (/ xsize(1)*xsize(2)*xsize(3) /) ),"GaussF.vtk","V1")  
   
   !----------Verification de la normalisation
   sumV = sum(V1)  
   call MPI_Allreduce(MPI_IN_PLACE,SumV,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
   if (nrank==0) print * , "GAUSS : sum(V1) = ", sumV

   !----------Application filtre Gaussien dans l'ESPACE REEL
   V2 = reshape((2*PI*sigma*sigma) ** (-1.5_mytype) * exp(-(X*X+Y*Y+Z*Z)/(2._mytype * sigma*sigma)),&
                 (/ xsize(1)*xsize(2)*xsize(3) /))
   V2 = V2 * dx*dy*dz
   !call print_field_vtk(reshape(V2, (/ xsize(1)*xsize(2)*xsize(3) /) ),"GaussR.vtk","V2") 
   
   !-----------Verification du filtre Fourier vs Reel sur une petite zone
   V1r = reshape(V1, (/ xsize(1),xsize(2),xsize(3) /))
   V2r = reshape(V2, (/ xsize(1),xsize(2),xsize(3) /))
   
   if (xstart(1) == 1 .and. xstart(2) == 1 .and. xstart(3) == 1) then
      max_err = maxval(abs((V1r(1:Nsigma,1:Nsigma,1:Nsigma) - V2r(1:Nsigma,1:Nsigma,1:Nsigma))/V1r(1:Nsigma,1:Nsigma,1:Nsigma)))
      print * , "GAUSS : max_err_rel (maxval(abs(Delta_V / V))) sur domaine (1:Nsigma)^3 = ", max_err / V1r(1,1,1)
   end if

!!------------------------------------------------------------------------------
!>                       TEST FILTRE REEL "SOMME / 6 voisins pour Champ de phase 
!>                                                             

   !----------Initialisation : Dirac en (1,1,1)
   V0 = 0_mytype
   if (xstart(1) == 1 .and. xstart(2) == 1 .and. xstart(3) == 1) V0(1) = 1.

   V1 = V0
   call field_filter(V1,1000,1,1,(/grid%nx,grid%ny,grid%nz/),(/grid%dx,grid%dy,grid%dz /))

   ! validation visuelle sur le champ : 6 voxels allumés 3 dans le coin (1,1,1), 
   !                                    puis 1 par coin le plus proche
   !call print_field_vtk(reshape(V1, (/ xsize(1)*xsize(2)*xsize(3) /) ),"V1.vtk","V1")   
   

  call decomp_2d_fft_finalize
  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)

!=======================================================================
CONTAINS

subroutine diff_fields(field1,field2,Dmax,typediff)
    implicit none
    real(mytype),dimension(xsize(1)*xsize(2)*xsize(3)),intent(in) :: field1,field2
    integer,intent(in)       :: typediff
    real(mytype),intent(out) :: Dmax
    
    integer :: nv    
    integer :: ierror
    
    nv = xsize(1)*xsize(2)*xsize(3)
    Dmax = 0.
    
    if (typediff == 1) then
      Dmax = maxval(abs(field1 - field2))
      call MPI_Allreduce(MPI_IN_PLACE,Dmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
    else
      call amitex_abort("diff_fields : typediff different from 1 - not yet implemented " ,2,0)  
    end if

end subroutine diff_fields

end program test_field_filterF


