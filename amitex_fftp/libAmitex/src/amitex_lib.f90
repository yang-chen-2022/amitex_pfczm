module amitex_lib

  use ISO_FORTRAN_ENV

  use MPI
#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi,       only : nrank
  use decomp_2d,           only : xstart,xend,xsize
  use decomp_2d_fft,       only : decomp_2d_fft_get_sp,decomp_2d_fft_get_ph
#else
  use decomp_2d,      only : mytype, nrank, xstart,xend,xsize
  use decomp_2d_fft,  only : get_decomp_fft_info
#endif  

  use amitex_mod
  use green_mod
  use field_mod
  use error_mod
  use io2_amitex_mod
  
  use simu_mod

  private

  public :: init_amitex_lib_0, init_amitex_lib

contains
!===================================================================
!===================================================================
!TODO : REGROUPER init_amitex_lib_0 et init_amitex_lib
!              et ajouter les initialisation 2decomp (voir simu_mod)
!
!===================================================================
!===================================================================
!
!               SUBROUTINE INIT_AMITEX_LIB_0
!
!>   First set of initializations for using AMITEX as a library
!!
!!             TO BE USED "JUST AFTER MPI_INIT"
!!
!!           Initilization - Association of :
!!              
!!                  Flog      
!!                  fic_log   
!!                  fic_vtk   (maybe not necessary here)
!!                  fft_start 
!!                  fft_end   
!!                  fft_size  
!!                  amx_start 
!!                  amx_end   
!!                  amx_size  
!!                  times_f   
!!                  error      
!!
!!
!===================================================================
subroutine init_amitex_lib_0()

  implicit none
!
! COPY FROM simu_mod, TYPE SIMU_AMITEX, restriction to usefull variables
!  
!        !---amitex_mod.f90 - fields (voir description in amitex_mod.f90)
!       x real(kind=REAL32),allocatable,dimension(:)       :: TEMPfield32   
!        !---amitex_mod.f90 - output files
!       * integer                                          :: Flog
!       * character(len=200)                               :: fic_vtk, fic_log
!        !---green_mod.f90
!       x real(mytype),allocatable,dimension(:,:,:,:)      :: FREQ
!       x type(GRID_DIM)                                   :: grid
!        !---field_mod.f90
!       * integer,dimension(3)                             :: fft_start,fft_end,fft_size
!       * type(TIMES_FIELD)                                :: times_f
!        !---error_mod.f90
!       * logical                                          :: error=.false.

  
  !> dummy initialization, modified after
  SIMU(1)%Flog = 8617          
  SIMU(1)%fic_log = "flogname" 
  SIMU(1)%fic_vtk = "fvtkname" 
  SIMU(1)%fft_start = 0
  SIMU(1)%fft_end = 0
  SIMU(1)%fft_size = 0
  SIMU(1)%amx_start = 0
  SIMU(1)%amx_end = 0
  SIMU(1)%amx_size = 0
  !SIMU(1)%times_f  ! objet initialisé dans le type
  !SIMU(1)%error    ! initialise
  
  Flog      => SIMU(1)%Flog
  fic_log   => SIMU(1)%fic_log
  fic_vtk   => SIMU(1)%fic_vtk
  fft_start => SIMU(1)%fft_start
  fft_end   => SIMU(1)%fft_end
  fft_size  => SIMU(1)%fft_size
  amx_start => SIMU(1)%amx_start
  amx_end   => SIMU(1)%amx_end
  amx_size  => SIMU(1)%amx_size
  times_io  => SIMU(1)%times_io
  times_f   => SIMU(1)%times_f
  times_g   => SIMU(1)%times_g
  error     => SIMU(1)%error 

end subroutine init_amitex_lib_0

!===================================================================
!
!               SUBROUTINE INIT_AMITEX_LIB
!
!>    Initializations for using AMITEX as a library
!!
!!         TO BE USED "JUST AFTER INITIALIZING 2DECOMP/FFT"
!!
!!         - open the log file
!!         - Allocation - Initilization - Association of
!!                  frequencies FREQ
!!                  TEMPfields32, for writing vtk fields 
!!
!!
!!
!===================================================================
subroutine init_amitex_lib(nx,ny,nz,dx,dy,dz,x0,y0,z0,ficlog,filter_type,filter_radius)

  implicit none 

  character(len=*), intent(in),optional          :: ficlog
  character(len=*),intent(in),optional           :: filter_type
  real(mytype),intent(in),optional               :: filter_radius   !no_filter, hexa, octa
  integer, intent(in)                            :: nx,ny,nz
  real(mytype),intent(in)                        :: dx,dy,dz
  real(mytype),intent(in)                        :: x0,y0,z0
  integer                                        :: alloc_stat
  character(200)                                 :: filter_type0
  real(mytype)                                   :: filter_radius0

  !> Ouverture du fichier log
  if (present(ficlog)) then
     call init_log(ficlog)
  else
     call init_log()
  end if
  
  !> Recuperation x(start,end,size)  
  amx_start = xstart
  amx_end   = xend
  amx_size  = xsize

  !> Initialisation de 'grid' (module green_mod) : grid
  call initGrid(nx,ny,nz,dx,dy,dz,x0,y0,z0)  ! initialize grid0
  SIMU(1)%grid = grid
  grid => SIMU(1)%grid  

  !> Allocation - Initilisation - Association des frequences FREQ
  filter_type0 = "hexa"
  filter_radius0 = 1

  if (present(filter_type))   filter_type0   = filter_type
  if (present(filter_radius)) filter_radius0 = filter_radius

  allocate(SIMU(1)%FREQ(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_amitex_lib)",2,0)
  
  call initFREQ(SIMU(1)%FREQ,filter_type0,filter_radius0) ! affecte grid0

  FREQ => SIMU(1)%FREQ

  !> Allocation - Initilisation - Association de TEMPfield32 
  allocate(SIMU(1)%TEMPfield32(amx_size(1)*amx_size(2)*amx_size(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort(&
                     "espace memoire disponible insuffisant (init_amitex_lib)",2,0)
                     
  SIMU(1)%TEMPfield32=0.
 
  TEMPfield32 => SIMU(1)%TEMPfield32

  !> 
#ifdef NEW_2DECOMP   
  SIMU(1)%ph_decomp => decomp_2d_fft_get_ph()
  SIMU(1)%sp_decomp => decomp_2d_fft_get_sp()
  ph_decomp => SIMU(1)%ph_decomp
  sp_decomp => SIMU(1)%sp_decomp
#else
  call get_decomp_fft_info(SIMU(1)%ph_decomp,SIMU(1)%sp_decomp)
#endif


end subroutine init_amitex_lib


end module amitex_lib
