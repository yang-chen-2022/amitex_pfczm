!===================================================================
!> Code pour lire des fichiers VTK en grandes transformations et les reecrire sur 
!! le maillage deforme
!!
!! ATTENTION : Dep est le seul tableau de rang (ncomp,:) et non (:,ncomp) 
!!		(pour Dep, ncomp=3)
!!
!!==============================================================================
program deformedShape

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
  use amitex_lib

  use deformed_shape_mod

  use error_mod
  use io_amitex_mod
  use amitex_mod   , only : Flog, &
                            set_prow_pcol ! function
  use green_mod
  use field_mod
  
  use simu_mod

  implicit none


  character(len=200)                      :: fic_in, fic_out
  real(mytype),dimension(:,:),allocatable :: dep, field
  character(len=20)                       :: numS                      !< pas de chargement a traiter 
  integer                                 :: nx, ny, nz                !< dimensions de la cellule
  integer(KIND=8)                         :: ntot                      !< nb. de voxels (nx*ny*nz)
  real(mytype)                            :: dx,dy,dz                  !< spacing (parametre de l'entete vtk)
  real(mytype)                            :: x0,y0,z0                  !< origin (parametre de l'entete vtk)
  integer                                 :: p_row, p_col              !< nb. de lignes, de colonnes pour la decomposition de 2decomp
                                                                       !< Description des pinceaux de l'espace spectral (FFT)
  integer                                 :: ierror                    !< erreur relative au fonction MPI
  integer                                 :: alloc_stat                !< erreur lors d'une allocation memoire
  integer                                 :: type_mpi                  !< type MPI lu par read_header_vtk
  complex(mytype),allocatable,dimension(:,:,:,:)    :: GradUF          !< transformee de Fourier du gradient du deplacement
  complex(mytype),allocatable,dimension(:,:,:,:)    :: DepF            !< transformee de Fourier du deplacement

  !> Rayon du filtre
  real(mytype)                            		:: filter_radius
  !> Type de filtre
  character(len=16)                       		:: filter_type

  !> Booleen pour savoir si les donnees d'entree sont dans un fichier ou plusieurs
  logical                                 		:: outprog
  character(len=16)                                      :: choice    !< choix de la ligne de commande -def, -sig ou -pk1

  !! Initialisation de MPI
  call MPI_INIT(ierror)
  call MPI_COMM_RANK(mpi_comm_world,nrank,ierror) !necessaire tant que 2decomp n'est pas initialise

  !! Initialization of amitex_lib : 1st part, just after MPI_INIT
  !!                                2nd part, after 2decomp/fft initialiszation
  allocate(SIMU(1))
  call init_amitex_lib_0()

  !! Recuperation du pas de calcul, de la racine des fichiers VTK a lire, du type de tenseur a utiliser
  call lire_commande_post(numS, fic_in, fic_out,choice,outprog)
  if(outprog) goto 666

  !! Lecture de l'en-tete du fichier de deformation pour connaitre les dimensions de la grille
  call read_header_vtk(trim(fic_in)//"_def1_"//trim(ADJUSTL(numS))//".vtk",nx,ny,nz,dx,dy,dz,x0,y0,z0,type_mpi)

  !! Initialisation MPI 2decomp
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  call set_pRow_pCol(p_row,p_col,nx,ny,nz)      !determination de p_row et p_col
  call decomp_2d_init(nx,ny,nz,p_row,p_col)     !initialisation decomp_2d
  call decomp_2d_fft_init                       !initialisation decomp_2d_fft
  call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)

  !! Initialization of amitex_lib (second part, juste after 2decomp_fft init_get_size)
  call init_amitex_lib(nx,ny,nz,dx,dy,dz,x0,y0,z0,ficlog="deformedShape.log",filter_type="no_filter",filter_radius=0._8)
  !call init_amitex_lib(nx,ny,nz,dx,dy,dz,ficlog="deformedShape.log",filter_type="hexa",filter_radius=1._8)

  if (nrank==0) write(Flog,"(A)") achar(10)//"log file for deformedShape"//achar(10)//&
                                             "--------------------------"//achar(10)
  call print_Grid(Flog,0)


  !!------------------------------------------------------------------------------
  !!                  Initialisation des tableaux dans l'espace reel
  !!
  !! ATTENTION : Dep est le seul tableau de rang (ncomp,:) et non (:,ncomp) 
  !!	=> necessaire pour l'ecriture VTK au format structures_grid
  !!		(Dep represente la position des points)
  !!------------------------------------------------------------------------------
  allocate(Dep(3,xsize(1)*xsize(2)*xsize(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (deformedShape)",2)
  allocate(field(xsize(1)*xsize(2)*xsize(3),9),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (deformedShape)",2)
  !!------------------------------------------------------------------------------
  !!                  Initialisation des tableaux dans l'espace de Fourier
  !!------------------------------------------------------------------------------
  allocate(GradUF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),&
       9),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex)",2)
  allocate(DepF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),&
       3),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex)",2)



  !! Possibilite de modifier les frequences....
  !! Choix du filtre
  !! Pour l'instant, uniquement le cas sans filtre fonctionne (a approfondir)
  ! filter_type ="hexa"
  filter_type ="no_filter"
  !! Choix du rayon du filtre
  if(filter_type == "hexa")then
     filter_radius = 1._mytype
  elseif(filter_type == "octa")then
     filter_radius = 2._mytype
  elseif(filter_type == "no_filter") then
     filter_radius = 0._mytype
  else
     call amitex_abort("Filtre inconnu (read_param)",1,0)
  end if
  !call initFreq(FREQ,trim(filter_type),filter_radius)
  ntot=nx*ny
  ntot=ntot*nz
  call computeDep(numS, fic_in, Dep, field,DepF,GradUF,nx,ny,nz,dx,dy,dz,ntot)
  call printVTK(numS, fic_in,fic_out,9,real(Dep,kind=REAL32),field,nx,ny,nz,choice)

  deallocate(dep,field,graduf,depf)
666  call MPI_FINALIZE(ierror)

end program deformedShape

!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
