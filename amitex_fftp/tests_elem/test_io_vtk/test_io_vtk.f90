program test_io_vtk
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

  integer,parameter               :: nx=65, ny=65, nz=65         !< dimensions de la cellule
  real(mytype),parameter          :: Lx=10, Ly=10, Lz=10            !< dimensions de la cellule
  real(mytype)                    :: dx, dy, dz,x0=0.,y0=0.,z0=0.     


!!------------------------------------------------------------------------------
!>                                                                  PARALLELISME

  integer :: p_row=2, p_col=7           !< nb. de lignes, de colonnes pour la decomposition de 2decomp
  integer :: ierror                     !< erreur relative au fonction MPI
!  integer, dimension(3) :: fft_start, fft_end, fft_size
                                        !< Description des pinceaux de l'espace spectral (FFT)
                                        !<    => defini dans green_mod

!!------------------------------------------------------------------------------
!>                                                           tableaux paralleles

  !> dans l'espace reel : tableaux (ntot)
  real(mytype),allocatable,dimension(:)                    :: Var,Var2 
  real(mytype)                                             :: max_err,min_err
  integer(kind=INT64),allocatable,dimension(:)             :: IVar
  

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

  call decomp_2d_init(nx,ny,nz,p_row,p_col)                !initialisation decomp_2d
  call decomp_2d_fft_init(PHYSICAL_IN_X, nx, ny, nz)       !initialisation decomp_2d_fft
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

  allocate(Var(xsize(1)*xsize(2)*xsize(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(Var2(xsize(1)*xsize(2)*xsize(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)

  allocate(IVar(xsize(1)*xsize(2)*xsize(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)

  Var(:)=real(nrank,mytype)
  Var2=0.
  Ivar = 0

!!==============================================================================
!!------------------------------------------------------------------------------
!>                                       ECRITURE (64bits -> 32 BITS par défaut)
!>                                                             

  if (nrank==0) then
    print *, "" 
    print *, " Write out_real32.vtk (64bits amitex-> 32 bits vtk) "
    print *, ""
  end if 

  call print_field_vtk(Var,"out_real32.vtk","tutu")

!!------------------------------------------------------------------------------
!>                                                     LECTURE 32 bits -> 64 bits
!>                                                             

  if (nrank==0) then
    print *, "" 
    print *, " Read out_real32.vtk (32bits vtk-> 64 bits amitex) "
    print *, ""
  end if 

  call read_field_vtk(Var2,"out_real32.vtk")

!!------------------------------------------------------------------------------
!>                                                                         CHECK

  call MPI_Allreduce(maxval(Var2 - Var), max_err, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
  call MPI_Allreduce(minval(Var2 - Var), min_err, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)

  if (nrank==0) then
    print *, "" 
    print *, "check (write - read) : maxval, minval", max_err, min_err 
    print *, ""
    if (max_err==0 .and. min_err==0) then
       print *, ""
       print *, "CHECK(assign field -> print_field_vtk -> read_field_vtk) : OK" 
       print *, "=============================================================" 
    else
       print *, ""
       print *, "PROBLEM"
       print *, "======="
    end if
  end if 

!!==============================================================================
!!------------------------------------------------------------------------------
!> LECTURE CHAMP ENTIER 32bits vtk (unsigned_short) -> entier 64bits (Ivar en 64b)
!>                                                             

  if (nrank==0) then
    print *, "" 
    print *, " Read mat_al2p5_65.vtk (32bits vtk-> 64 bits amitex) "
    print *, ""
  end if 

  call read_field_vtk(IVar,"../../cas_tests/microstructures/al2p5_65/mat_al2p5_65.vtk")
 

!!------------------------------------------------------------------------------
!>                ECRITURE CHAMP ENTIER (64 bits amitex -> flottant 32 bits vtk)
!>                                                             

  if (nrank==0) then
    print *, "" 
    print *, " Write out_real32.vtk (64bits amitex-> 32 bits vtk) "
    print *, ""
  end if 

  Var = Ivar ! entier 64 -> flottant 64 bits
  call print_field_vtk(Var,"out_real32.vtk","tutu")

!!------------------------------------------------------------------------------
!>                                          LECTURE 32 bits vtk-> 64 bits amitex
!>                                                             

  if (nrank==0) then
    print *, "" 
    print *, " Read out_real32.vtk (64bits vtk-> 8 bits amitex) "
    print *, ""
  end if 

  call read_field_vtk(Var2,"out_real32.vtk")

 !!------------------------------------------------------------------------------
 !>                                                                         CHECK

  call MPI_Allreduce(maxval(Var2 - Var), max_err, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
  call MPI_Allreduce(minval(Var2 - Var), min_err, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)

  if (nrank==0) then
    print *, "" 
    print *, "check (write - read) : maxval, minval", max_err, min_err 
    print *, ""
    if (max_err==0 .and. min_err==0) then
       print *, ""
       print *, "CHECK(read_field_vtk (int) -> print_field_vtk -> read_field_vtk) : OK" 
       print *, "=====================================================================" 
    else
       print *, ""
       print *, "PROBLEM"
       print *, "=======" 
    end if

  end if 


  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)

end program test_io_vtk
