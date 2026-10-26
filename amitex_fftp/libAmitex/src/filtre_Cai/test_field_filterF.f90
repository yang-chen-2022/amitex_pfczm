program test_field_filterF
!> CHARGEMENT DES DIFFERENTS MODULES

  use ISO_FORTRAN_ENV

  use MPI
  use decomp_2d
  use decomp_2d_fft
  use io2_amitex_mod
  use amitex_mod
  use error_mod
  use green_mod
  use field_mod
  use caifilters_mod



!!==============================================================================
!!                                                                  DECLARATIONS
!!==============================================================================
  implicit none
!!------------------------------------------------------------------------------
!>                                                       DIMENSIONS DE LA GRILLE   : ATTENTION, LECTURE D'UN VTK DE 65x65x65
!>                                                                       => ON DOIT UTILISER nx=ny=nz=65

  integer,parameter               :: nx=65, ny=65, nz=65            !< dimensions de la cellule
  real(mytype),parameter          :: Lx=10, Ly=10, Lz=10            !< dimensions de la cellule
  real(mytype)                    :: dx, dy, dz


!!------------------------------------------------------------------------------
!>                                                                  PARALLELISME

  integer :: p_row=1, p_col=1           !< nb. de lignes, de colonnes pour la decomposition de 2decomp
  logical,dimension(3)         :: periodic_bc=.true.
  integer :: ierror                     !< erreur relative au fonction MPI
  integer, dimension(3),target :: fft_start0, fft_end0, fft_size0
  type(TIMES_FIELD),target     :: times_f0       
  type(TIMES_IO2),target       :: times_io0


!!------------------------------------------------------------------------------
!>                                                           tableaux paralleles

  !> dans l'espace reel : tableaux (ntot)
  real(mytype),allocatable,dimension(:)                    :: V0, V1, Vm1, V2 
  real(mytype),allocatable,dimension(:,:,:)                :: V0c ! copy 3D and halo of V0
  real(mytype),allocatable,dimension(:,:,:)                :: V1r, Vm1r, V2r, X,Y,Z, Filtre
  
  complex(mytype),allocatable,dimension(:,:,:)             :: VF,FiltreF
  real(mytype)                                             :: max_err,sumV,sigma,acai
  integer                                                  :: i,j,k,Nsigma,Ncai
  double precision                                         :: t0,t1
  real(mytype),parameter                                   :: PI = 4._mytype*DATAN(1._mytype)
  
!!------------------------------------------------------------------------------
!>                                                           divers

  integer :: alloc_stat                 !< erreur lors d'une allocation memoire
  integer :: io_stat                    !< erreur lors d'une lecture ou ecriture de fichier

!!------------------------------------------------------------------------------
!>                                                            INITIALISATION MPI

  call MPI_INIT(ierror)

!!------------------------------------------------------------------------------
!>                                                     INITIALISATION FLOG, FVTK
  Flog0 = 11487
  Flog => Flog0
  fic_log0 = "titi"
  fic_vtk0 = "toto"
  fic_log => fic_log0
  fic_vtk => fic_vtk0
  times_io => times_io0

  if(nrank==0) then
     open(unit=Flog, file="sortie.log",form="formatted", status="replace", action="write",iostat= io_stat)
     if ( io_stat /= 0 ) then
        write(*,"(A,I0)")" Probleme a l'ouverture (fichier: sortie.log) (amitex)",io_stat
     end if
  end if

!!------------------------------------------------------------------------------
!>                                           INITIALISATION 2DECOMP, 2DECOMP_FFT

  call decomp_2d_init(nx,ny,nz,p_row,p_col,periodic_bc)     !initialisation decomp_2d
  call decomp_2d_fft_init 

  call decomp_2d_fft_get_size(fft_start0,fft_end0,fft_size0) 
  fft_start => fft_start0
  fft_end   => fft_end0
  fft_size  => fft_size0
  
  times_f => times_f0
  
  if(nrank==0)print *, "Dimensions : ",nx,"x",ny,"x",nz, p_row,p_col

!!------------------------------------------------------------------------------
!>                                                   INITIALISATION OBJET GRILLE

  dx = Lx/nx
  dy = Ly/ny
  dz = Lz/nz
  call initGrid(nx,ny,nz,dx,dy,dz,0._mytype,0._mytype,0._mytype)
  grid => grid0 
  
  
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
  allocate(Filtre(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)

  allocate(VF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(FiltreF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)

! definition des coordonnes de la grille (pout alidation filtre Gauss)
  call gene_cartesian_coord(X,Y,Z,dx,dy,dz)
  X = X - dx / 2._mytype
  Y = Y - dx / 2._mytype
  Z = Z - dz / 2._mytype

! On assigne V0 a la valeur du champ  
  V0  = real(nrank,mytype)
  V0c = real(nrank,mytype)
!  call print_field_vtk(V0,"V0.vtk","V0") 
 
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


!!==============================================================================
!>                                                          TEST FILTRE GAUSSIEN
!>                                                             
   Nsigma = 3
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
   call MPI_Allreduce(MPI_IN_PLACE,SumV,1,MPI_DOUBLE_PRECISION,MPI_SUM,amitex_mpi_comm,ierror)
   if (nrank==0) print * , "GAUSS : sum(V1) = ", sumV

   !----------Application filtre Gaussien dans l'ESPACE REEL
   V2 = reshape((2._mytype*PI*sigma*sigma) ** (-1.5_mytype) * exp(-(X*X+Y*Y+Z*Z)/(2._mytype * sigma*sigma)),&
                 (/ xsize(1)*xsize(2)*xsize(3) /))
   V2 = V2 * dx*dy*dz
   !call print_field_vtk(reshape(V2, (/ xsize(1)*xsize(2)*xsize(3) /) ),"GaussR.vtk","V2") 
   
   !----------Verification de la normalisation
   sumV = sum(V2)  
   call MPI_Allreduce(MPI_IN_PLACE,SumV,1,MPI_DOUBLE_PRECISION,MPI_SUM,amitex_mpi_comm,ierror)
   if (nrank==0) print * , "GAUSS : sum(V2) = ", sumV
   
   !-----------Verification du filtre Fourier vs Reel sur une petite zone
   V1r = reshape(V1, (/ xsize(1),xsize(2),xsize(3) /))
   V2r = reshape(V2, (/ xsize(1),xsize(2),xsize(3) /))
   
   if (xstart(1) == 1 .and. xstart(2) == 1 .and. xstart(3) == 1) then
!print *, V1(1:2*Nsigma)
!print *, V2(1:2*Nsigma)
      max_err = maxval(abs((V1r(1:Nsigma,1:Nsigma,1:Nsigma) - V2r(1:Nsigma,1:Nsigma,1:Nsigma))/V1r(1:Nsigma,1:Nsigma,1:Nsigma)))
      print * , "GAUSS : max_err_rel (maxval(abs(Delta_V)) / V(1,1,1)) sur domaine (1:Nsigma)^3 = ", max_err / V1r(1,1,1)
   end if



!!==============================================================================
!>                                                               TEST FILTRE CAI
!>                                                             
   Ncai = 2    
   acai = Ncai * grid%dx

   !----------Initialisation : Dirac en (1,1,1)
   V0 = 0_mytype
   if (xstart(1) == 1 .and. xstart(2) == 1 .and. xstart(3) == 1) V0(1) = 1.
   
   !----------Application filtre de Cai dans l'ESPACE DE FOURIER (METHODE 1)
   call field_fft(V0,VF,1,1)
   call field_filter_caiF(VF,FILTER=200,ncomp=1,nvar=1,&
                   N=(/grid%nx,grid%ny,grid%nz /),&
                   d=(/grid%dx,grid%dy,grid%dz /),param=(/ acai /))
   call field_ifft(V1,VF,1,1)
   V1 = V1 / grid%ntot
   call print_field_vtk(reshape(V1, (/ xsize(1)*xsize(2)*xsize(3) /) ),"CaiF.vtk","V1")  
   
   !----------Verification de la normalisation
   sumV = sum(V1)  
   call MPI_Allreduce(MPI_IN_PLACE,SumV,1,MPI_DOUBLE_PRECISION,MPI_SUM,amitex_mpi_comm,ierror)
   if (nrank==0) print * , "CAI : sum(V1) = ", sumV

   
   !----------Application filtre de CAI dans l'ESPACE REEL
   V2 = reshape( (15._mytype/(8._mytype*PI*acai**3)) &
                *(1._mytype +  (X*X+Y*Y+Z*Z) / (acai**2)) ** (-7._mytype / 2._mytype),&
                 (/ xsize(1)*xsize(2)*xsize(3) /))
   V2 = V2 * dx*dy*dz
   call print_field_vtk(reshape(V2, (/ xsize(1)*xsize(2)*xsize(3) /) ),"CaiR.vtk","V2") 

   !----------Verification de la normalisation
   sumV = sum(V2)  
   call MPI_Allreduce(MPI_IN_PLACE,SumV,1,MPI_DOUBLE_PRECISION,MPI_SUM,amitex_mpi_comm,ierror)
   if (nrank==0) print * , "CAI : sum(V2) = ", sumV
   
   !-----------Verification du filtre Fourier vs Reel sur une petite zone
   V1r = reshape(V1, (/ xsize(1),xsize(2),xsize(3) /))
   V2r = reshape(V2, (/ xsize(1),xsize(2),xsize(3) /))
   
   if (xstart(1) == 1 .and. xstart(2) == 1 .and. xstart(3) == 1) then
      max_err = maxval(abs((V1r(1:Ncai,1:Ncai,1:Ncai) - V2r(1:Ncai,1:Ncai,1:Ncai))))
      print * , "CAI : max_err_rel (maxval(abs(Delta_V)) / V(1,1,1)))) sur domaine (1:Ncai)^3 = ", max_err / V1r(1,1,1)
      print * , "CAI : err en 1,1,1 (abs(Delta_V(1,1,1))) / V(1,1,1)))) sur domaine (1:Ncai)^3 = ", abs(V1(1)-V2(1)) / V1r(1,1,1)
   end if

   !----------Application filtre de Cai dans l'ESPACE DE FOURIER (METHODE 2)
   Filtre=0
   do i = 0,1
   do j = 0,1
   do k = 0,1
      where((X-i*Lx) <Lx/2. .and. (X-i*Lx) >= -Lx/2. .and. & 
            (Y-j*Ly) <Ly/2. .and. (Y-j*Ly) >= -Ly/2. .and. & 
            (Z-k*Lz) <Lz/2. .and. (Z-k*Lz) >= -Lz/2.) 
      Filtre = Filtre + (15._mytype/(8._mytype*PI*acai**3)) &
           *(1._mytype +  ((X-i*Lx)**2+(Y-j*Ly)**2+(Z-k*Lz)**2) / (acai**2)) ** (-7._mytype / 2._mytype)
      end where
   end do
   end do
   end do
   Filtre = Filtre *  dx*dy*dz

   call field_fft(Filtre,FiltreF,1,1)
   call field_fft(V0,VF,1,1)
   VF = FiltreF*VF
   call field_ifft(V1,VF,1,1)
   V1 = V1 / grid%ntot
   call print_field_vtk(reshape(V1, (/ xsize(1)*xsize(2)*xsize(3) /) ),"CaiF2.vtk","V1")  

   !----------Verification de la normalisation
   sumV = sum(V1)  
   call MPI_Allreduce(MPI_IN_PLACE,SumV,1,MPI_DOUBLE_PRECISION,MPI_SUM,amitex_mpi_comm,ierror)
   if (nrank==0) print * , "CAI2 : sum(V1) = ", sumV
   
   
   !-----------Verification du filtre Fourier vs Reel sur une petite zone
   V1r = reshape(V1, (/ xsize(1),xsize(2),xsize(3) /))
   V2r = reshape(V2, (/ xsize(1),xsize(2),xsize(3) /))
   
   if (xstart(1) == 1 .and. xstart(2) == 1 .and. xstart(3) == 1) then
      max_err = maxval(abs((V1r(1:Ncai,1:Ncai,1:Ncai) - V2r(1:Ncai,1:Ncai,1:Ncai))))
      print * , "CAI2 : max_err_rel (maxval(abs(Delta_V)) / V(1,1,1)))) sur domaine (1:Ncai)^3 = ", max_err / V1r(1,1,1)
      print * , "CAI2 : err en 1,1,1 (abs(Delta_V(1,1,1))) / V(1,1,1)))) sur domaine (1:Ncai)^3 = ", abs(V1(1)-V2(1)) / V1r(1,1,1)
   end if

  
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
      call MPI_Allreduce(MPI_IN_PLACE,Dmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,amitex_mpi_comm,ierror)
    else
      call amitex_abort("diff_fields : typediff different from 1 - not yet implemented " ,2,0)  
    end if

end subroutine diff_fields

end program test_field_filterF


