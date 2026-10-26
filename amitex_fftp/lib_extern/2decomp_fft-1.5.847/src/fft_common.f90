!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This file contains common code shared by all FFT engines

  integer, parameter, public :: DECOMP_2D_FFT_FORWARD = -1
  integer, parameter, public :: DECOMP_2D_FFT_BACKWARD = 1
  
  ! Physical space data can be stored in either X-pencil or Z-pencil
  integer, parameter, public :: PHYSICAL_IN_X = 1
  integer, parameter, public :: PHYSICAL_IN_Z = 3 

  integer, pointer, save :: format                 ! input X-pencil or Z-pencil
  
  ! The libary can only be initialised once
  integer, save :: initialised = 0 

  ! Global size of the FFT
  integer, pointer, save :: nx_fft, ny_fft, nz_fft

  ! 2D processor grid (common to the different grids)
  integer, save, dimension(2) :: dims

  ! Decomposition objects (internal use)
  TYPE(DECOMP_INFO), pointer, save :: ph  ! physical space
  TYPE(DECOMP_INFO), pointer, save :: sp  ! spectral space

  ! Workspace to store the intermediate Y-pencil data
  ! *** TODO: investigate how to use only one workspace array
  complex(mytype), pointer, dimension(:,:,:) :: wk2_c2c, wk2_r2c
  complex(mytype), pointer, dimension(:,:,:) :: wk13

  ! Derived-type gathering informations for 2decomp_2D_fft for multigrid purpose
  type DECOMP_FFT_MULTIGRID
      integer                  :: format
      integer                  :: nx_fft=0, ny_fft=0, nz_fft=0 
      TYPE(DECOMP_INFO)        :: ph  ! physical space
      TYPE(DECOMP_INFO)        :: sp  ! spectral space
      complex(mytype), allocatable, dimension(:,:,:) :: wk2_c2c, wk2_r2c
      complex(mytype), allocatable, dimension(:,:,:) :: wk13
      !plan?
  end type DECOMP_FFT_MULTIGRID
  
  type(DECOMP_FFT_MULTIGRID), allocatable, dimension(:), target :: FFT_multigrid
  type(DECOMP_FFT_MULTIGRID), allocatable, dimension(:)         :: FFT_multigrid_tmp ! for reallocation purpose (move_alloc)

  public :: decomp_2d_fft_init, decomp_2d_fft_3d, &
       decomp_2d_fft_finalize, decomp_2d_fft_get_size,&
       associate_pointers_decomp_2d_fft, get_decomp_fft_info !,FFT_multigrid

  
  ! Declare generic interfaces to handle different inputs
  
  interface decomp_2d_fft_init
     module procedure fft_init_noarg
     module procedure fft_init_arg
     module procedure fft_init_general
     module procedure fft_init_general_multigrid
  end interface
  
  interface decomp_2d_fft_3d
     module procedure fft_3d_c2c
     module procedure fft_3d_r2c
     module procedure fft_3d_c2r
  end interface

  
contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Return the decomposition objects of 'active' pointers
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine get_decomp_fft_info(ph_decomp, sp_decomp)

    implicit none

    TYPE(DECOMP_INFO), intent(OUT) :: ph_decomp, sp_decomp

    integer               :: errorcode

    if (.not. associated(ph) .OR. .not. associated(sp)) then
      errorcode = 4
      call decomp_2d_abort(errorcode, &
        'decomp_2d_fft, get_decomp_fft_info : sp or ph not associated')    
    end if
    
    ph_decomp = ph
    sp_decomp = sp
        
    return
  end subroutine get_decomp_fft_info
    

    
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialise the FFT module
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_init_noarg
    
    implicit none
    
    call fft_init_arg(PHYSICAL_IN_X)  ! default input is X-pencil data
    
    return
  end subroutine fft_init_noarg

  subroutine fft_init_arg(pencil)     ! allow to handle Z-pencil input

    implicit none

    integer, intent(IN) :: pencil

    call fft_init_general(pencil, nx_global, ny_global, nz_global)

    return
  end subroutine fft_init_arg

  subroutine fft_init_general(pencil, nx, ny, nz)
    implicit none

    integer, intent(IN) :: pencil
    integer, intent(IN) :: nx, ny, nz
  
    integer :: Igrid
    
    Igrid = 1
  
    call fft_init_general_multigrid(pencil, nx, ny, nz, Igrid)
  
  end subroutine fft_init_general

  ! Initialise the FFT library to perform arbitrary size transforms on various grids
  !       Increase size(FFT_multigrid) if necessary
  !       Initialize FFT_multigrid(Igrid) (if not yet initialized)
  !       Associate pointers
  subroutine fft_init_general_multigrid(pencil, nx, ny, nz, Igrid)

    implicit none

    integer, intent(IN) :: pencil
    integer, intent(IN) :: nx, ny, nz
    integer, intent(IN) :: Igrid
    
    logical, dimension(2) :: dummy_periods
    integer, dimension(2) :: dummy_coords
    integer               :: status, errorcode, ierror
    integer               :: Ngrid0

    ! allocate FFT_multigrid if first call : size = Igrid
    if (initialised == 0) then
      allocate(FFT_multigrid(Igrid))
    end if

    Ngrid0 = size(FFT_multigrid)

    ! reallocate FFT_multigrid if Igrid > size(FFT_multigrid) : size = Igrid
    if (Igrid > Ngrid0) then
       allocate(FFT_multigrid_tmp(Igrid))
       FFT_multigrid_tmp(1:Ngrid0) = FFT_multigrid
       call move_alloc(FROM=FFT_multigrid_tmp,TO=FFT_multigrid)
    end if
    
    ! check if Igrid has been initialized before with different properties
    if (FFT_multigrid(Igrid)%nx_fft /= 0 .AND. FFT_multigrid(Igrid)%nx_fft /= nx &
                                         .AND. FFT_multigrid(Igrid)%nx_fft /= ny &
                                         .AND. FFT_multigrid(Igrid)%nx_fft /= nz &
                                         .AND. FFT_multigrid(Igrid)%format /= pencil) then
       errorcode = 4
       call decomp_2d_abort(errorcode, &
                'decomp_2d_fft, fft_init_general_multigrid : FFT_multigrid(Igrid) already initialized')    
    end if
    
    !
    !-------------------(Initialize AND associate) pointers if NOT initialized, 
    !                            OR
    !                   associate pointers if initialized
    !
    if (FFT_multigrid(Igrid)%nx_fft == 0) then ! initialize if not initialized before 
    FFT_multigrid(Igrid)%format = pencil
    FFT_multigrid(Igrid)%nx_fft = nx
    FFT_multigrid(Igrid)%ny_fft = ny
    FFT_multigrid(Igrid)%nz_fft = nz
    format => FFT_multigrid(Igrid)%format
    nx_fft => FFT_multigrid(Igrid)%nx_fft
    ny_fft => FFT_multigrid(Igrid)%ny_fft
    nz_fft => FFT_multigrid(Igrid)%nz_fft
 
    ! determine the processor grid in use
    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
         dims, dummy_periods, dummy_coords, ierror)

    ! for c2r/r2c interface:
    ! if in physical space, a real array is of size: nx*ny*nz
    ! in spectral space, the complex array is of size:
    !         (nx/2+1)*ny*nz, if PHYSICAL_IN_X
    !      or nx*ny*(nz/2+1), if PHYSICAL_IN_Z

    call decomp_info_init(nx, ny, nz, FFT_multigrid(Igrid)%ph)
    if (format==PHYSICAL_IN_X) then
       call decomp_info_init(nx/2+1, ny, nz, FFT_multigrid(Igrid)%sp)
    else if (format==PHYSICAL_IN_Z) then
       call decomp_info_init(nx, ny, nz/2+1, FFT_multigrid(Igrid)%sp)
    end if

    ph => FFT_multigrid(Igrid)%ph
    sp => FFT_multigrid(Igrid)%sp

    allocate(FFT_multigrid(Igrid)%wk2_c2c(ph%ysz(1),ph%ysz(2),ph%ysz(3)), STAT=status)
    allocate(FFT_multigrid(Igrid)%wk2_r2c(sp%ysz(1),sp%ysz(2),sp%ysz(3)), STAT=status)
    wk2_c2c => FFT_multigrid(Igrid)%wk2_c2c
    wk2_r2c => FFT_multigrid(Igrid)%wk2_r2c
    
    if (format==PHYSICAL_IN_X) then
       allocate(FFT_multigrid(Igrid)%wk13(sp%xsz(1),sp%xsz(2),sp%xsz(3)), STAT=status)
    else if (format==PHYSICAL_IN_Z) then
       allocate(FFT_multigrid(Igrid)%wk13(sp%zsz(1),sp%zsz(2),sp%zsz(3)), STAT=status)
    end if
    
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising FFT')
    end if
    wk13 => FFT_multigrid(Igrid)%wk13

    call init_fft_engine(Igrid)
    
    else
          call associate_pointers_decomp_2d_fft(Igrid) ! associate pointers if already initialized
          call associate_pointers_decomp_2d_fft_spec(Igrid)
    end if
    initialised = initialised + 1
    return
  end subroutine fft_init_general_multigrid
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! associate pointers for decomp_2d_fft
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine associate_pointers_decomp_2d_fft(Igrid)

  implicit none
  integer,intent(in) :: Igrid
  
  integer            :: errorcode

  if (Igrid > size(FFT_multigrid)) then
       errorcode = 4
       call decomp_2d_abort(errorcode, &
            'decomp_2d_fft, associate_pointers_decomp_2d_fft : Igrid > size(FFT_multigrid)')
  end if
  
  format    => FFT_multigrid(Igrid)%format
  nx_fft    => FFT_multigrid(Igrid)%nx_fft
  ny_fft    => FFT_multigrid(Igrid)%ny_fft
  nz_fft    => FFT_multigrid(Igrid)%nz_fft
  ph        => FFT_multigrid(Igrid)%ph  ! physical space
  sp        => FFT_multigrid(Igrid)%sp  ! spectral space
  wk2_c2c   => FFT_multigrid(Igrid)%wk2_c2c
  wk2_r2c   => FFT_multigrid(Igrid)%wk2_r2c
  wk13      => FFT_multigrid(Igrid)%wk13

  call associate_pointers_decomp_2d_fft_spec(Igrid)

  end subroutine associate_pointers_decomp_2d_fft

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Final clean up
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_fft_finalize
    
    implicit none
    
    integer :: Igrid
    
    do Igrid = 1,size(FFT_multigrid)

      call decomp_info_finalize(FFT_multigrid(Igrid)%ph)
      call decomp_info_finalize(FFT_multigrid(Igrid)%sp)

      deallocate(FFT_multigrid(Igrid)%wk2_c2c,&
                 FFT_multigrid(Igrid)%wk2_r2c,&
                 FFT_multigrid(Igrid)%wk13)

    end do
    
    call finalize_fft_engine
    
    deallocate(FFT_multigrid)
    
    initialised = 0

    return
  end subroutine decomp_2d_fft_finalize


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Return the size, starting/ending index of the distributed array 
  !  whose global size is (nx/2+1)*ny*nz, for defining data structures
  !  in r2c and c2r interfaces
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_fft_get_size(istart, iend, isize,jstart, jend, jsize)
    
    implicit none
    integer, dimension(3), intent(OUT) :: istart, iend, isize
    integer, dimension(3), optional, intent(OUT) :: jstart, jend, jsize
        
    if (format==PHYSICAL_IN_X) then
       istart = sp%zst
       iend   = sp%zen
       isize  = sp%zsz
    else if (format==PHYSICAL_IN_Z) then
       istart = sp%xst
       iend   = sp%xen
       isize  = sp%xsz
    end if
    
    if (present(jstart) .and. present(jend) .and. present(jsize)) then
    if (format==PHYSICAL_IN_X) then
       jstart = ph%xst
       jend   = ph%xen
       jsize  = ph%xsz
    else if (format==PHYSICAL_IN_Z) then
       jstart = ph%zst
       jend   = ph%zen
       jsize  = ph%zsz
    end if
    end if
    
    return
  end subroutine decomp_2d_fft_get_size
