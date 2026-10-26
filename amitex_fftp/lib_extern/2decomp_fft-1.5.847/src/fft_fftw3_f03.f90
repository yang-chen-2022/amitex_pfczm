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

! This is the FFTW implementation of the FFT library using 
! the Fortran 2003 interface introduced in FFTW 3.3-beta1

module decomp_2d_fft
  
  use decomp_2d, only : mytype, DECOMP_INFO, &         ! TYPES
                        nrank, DECOMP_2D_COMM_CART_X,& ! VARIABLES common to all "grids"
                        nx_global, ny_global, nz_global, & ! variables used ONLY with decomp_2d_fft_init and NO ARGUMENT
                        decomp_info_init, decomp_info_finalize, decomp_2d_abort, & 
                        transpose_y_to_z, transpose_y_to_x, transpose_x_to_y, transpose_z_to_y  
  
  use, intrinsic :: iso_c_binding
  
  implicit none

  include "fftw3.f03"
  
  private        ! Make everything private unless declared public

  ! engine-specific global variables
  integer, save :: plan_type = FFTW_MEASURE

  ! FFTW plans
  ! j=1,2,3 corresponds to the 1D FFTs in X,Y,Z direction, respectively
  ! For c2c transforms: 
  !     use plan(-1,j) for  forward transform; 
  !     use plan( 1,j) for backward transform;
  ! For r2c/c2r transforms:
  !     use plan(0,j) for r2c transforms;
  !     use plan(2,j) for c2r transforms;
  type(C_PTR), pointer, dimension(:,:), save :: plan

  integer, parameter, public :: DECOMP_2D_FFT_FORWARD = -1
  integer, parameter, public :: DECOMP_2D_FFT_BACKWARD = 1
  
  ! Physical space data can be stored in either X-pencil or Z-pencil
  integer, parameter, public :: PHYSICAL_IN_X = 1
  integer, parameter, public :: PHYSICAL_IN_Z = 3 

  integer, pointer, save :: format          ! input X-pencil or Z-pencil

  ! Number of initialisation
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
  complex(mytype), pointer :: wk2_c2c(:,:,:), wk2_r2c(:,:,:), wk13(:,:,:)
  type(C_PTR), pointer     :: wk2_c2c_p, wk2_r2c_p, wk13_p

  ! Derived-type gathering informations for 2decomp_2D_fft for multigrid purpose
  type DECOMP_FFT_MULTIGRID
      integer                  :: format
      integer                  :: nx_fft=0, ny_fft=0, nz_fft=0 
      TYPE(DECOMP_INFO)        :: ph  ! physical space
      TYPE(DECOMP_INFO)        :: sp  ! spectral space
      complex(mytype), pointer :: wk2_c2c(:,:,:), wk2_r2c(:,:,:), wk13(:,:,:)
      type(C_PTR)              :: wk2_c2c_p, wk2_r2c_p, wk13_p   
      type(C_PTR)              :: plan(-1:2,3)   
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
    integer               :: status, errorcode
    integer(C_SIZE_T)     :: sz, ierror
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

    sz = ph%ysz(1)*ph%ysz(2)*ph%ysz(3)
    FFT_multigrid(Igrid)%wk2_c2c_p = fftw_alloc_complex(sz)
    call c_f_pointer(FFT_multigrid(Igrid)%wk2_c2c_p,FFT_multigrid(Igrid)%wk2_c2c,[ph%ysz(1),ph%ysz(2),ph%ysz(3)])
    wk2_c2c_p => FFT_multigrid(Igrid)%wk2_c2c_p
    wk2_c2c   => FFT_multigrid(Igrid)%wk2_c2c

    sz = sp%ysz(1)*sp%ysz(2)*sp%ysz(3)
    FFT_multigrid(Igrid)%wk2_r2c_p = fftw_alloc_complex(sz)
    call c_f_pointer(FFT_multigrid(Igrid)%wk2_r2c_p,FFT_multigrid(Igrid)%wk2_r2c,[sp%ysz(1),sp%ysz(2),sp%ysz(3)])
    wk2_r2c_p => FFT_multigrid(Igrid)%wk2_r2c_p
    wk2_r2c => FFT_multigrid(Igrid)%wk2_r2c

    if (format==PHYSICAL_IN_X) then
       sz = sp%xsz(1)*sp%xsz(2)*sp%xsz(3)
       FFT_multigrid(Igrid)%wk13_p = fftw_alloc_complex(sz)
       call c_f_pointer(FFT_multigrid(Igrid)%wk13_p,FFT_multigrid(Igrid)%wk13,[sp%xsz(1),sp%xsz(2),sp%xsz(3)])
    else if (format==PHYSICAL_IN_Z) then
       sz = sp%zsz(1)*sp%zsz(2)*sp%zsz(3)
       FFT_multigrid(Igrid)%wk13_p = fftw_alloc_complex(sz)
       call c_f_pointer(FFT_multigrid(Igrid)%wk13_p,FFT_multigrid(Igrid)%wk13,[sp%zsz(1),sp%zsz(2),sp%zsz(3)])
    end if
    wk13_p => FFT_multigrid(Igrid)%wk13_p
    wk13   => FFT_multigrid(Igrid)%wk13
    
    call init_fft_engine(Igrid)
    plan => FFT_multigrid(Igrid)%plan
    
    else 
      call associate_pointers_decomp_2d_fft(Igrid) ! associate pointers if already initialized
    end if ! initialize if not initialized before 
    
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
  wk2_c2c_p => FFT_multigrid(Igrid)%wk2_c2c_p
  wk2_r2c_p => FFT_multigrid(Igrid)%wk2_r2c_p
  wk13_p    => FFT_multigrid(Igrid)%wk13_p   
  plan      => FFT_multigrid(Igrid)%plan  

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

      call fftw_free(FFT_multigrid(Igrid)%wk2_c2c_p)
      call fftw_free(FFT_multigrid(Igrid)%wk2_r2c_p)
      call fftw_free(FFT_multigrid(Igrid)%wk13_p)
      
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


  ! Return a FFTW3 plan for multiple 1D c2c FFTs in X direction
  subroutine c2c_1m_x_plan(plan1, decomp, isign)

    implicit none

    type(C_PTR) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: isign

    complex(mytype), pointer :: a1(:,:,:)
    type(C_PTR) :: a1_p
    integer(C_SIZE_T) :: sz

    sz = decomp%xsz(1)*decomp%xsz(2)*decomp%xsz(3)
    a1_p = fftw_alloc_complex(sz)
    call c_f_pointer(a1_p,a1,[decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)])

#ifdef DOUBLE_PREC
    plan1 =  fftw_plan_many_dft(1, decomp%xsz(1), &
         decomp%xsz(2)*decomp%xsz(3), a1, decomp%xsz(1), 1, &
         decomp%xsz(1), a1, decomp%xsz(1), 1, decomp%xsz(1), &
         isign, plan_type)
#else
    plan1 = fftwf_plan_many_dft(1, decomp%xsz(1), &
         decomp%xsz(2)*decomp%xsz(3), a1, decomp%xsz(1), 1, &
         decomp%xsz(1), a1, decomp%xsz(1), 1, decomp%xsz(1), &
         isign, plan_type)
#endif

    call fftw_free(a1_p)

    return
  end subroutine c2c_1m_x_plan

  ! Return a FFTW3 plan for multiple 1D c2c FFTs in Y direction
  subroutine c2c_1m_y_plan(plan1, decomp, isign)

    implicit none

    type(C_PTR) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: isign

    complex(mytype), pointer :: a1(:,:)
    type(C_PTR) :: a1_p
    integer(C_SIZE_T) :: sz

    ! Due to memory pattern of 3D arrays, 1D FFTs along Y have to be
    ! done one Z-plane at a time. So plan for 2D data sets here.
    sz = decomp%ysz(1)*decomp%ysz(2)
    a1_p = fftw_alloc_complex(sz)
    call c_f_pointer(a1_p,a1,[decomp%ysz(1),decomp%ysz(2)])

#ifdef DOUBLE_PREC
    plan1 =  fftw_plan_many_dft(1, decomp%ysz(2), decomp%ysz(1), &
         a1, decomp%ysz(2), decomp%ysz(1), 1, a1, decomp%ysz(2), &
         decomp%ysz(1), 1, isign, plan_type)
#else
    plan1 = fftwf_plan_many_dft(1, decomp%ysz(2), decomp%ysz(1), &
         a1, decomp%ysz(2), decomp%ysz(1), 1, a1, decomp%ysz(2), &
         decomp%ysz(1), 1, isign, plan_type)
#endif

    call fftw_free(a1_p)

    return
  end subroutine c2c_1m_y_plan

  ! Return a FFTW3 plan for multiple 1D c2c FFTs in Z direction
  subroutine c2c_1m_z_plan(plan1, decomp, isign)

    implicit none

    type(C_PTR) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: isign

    complex(mytype), pointer :: a1(:,:,:)
    type(C_PTR) :: a1_p
    integer(C_SIZE_T) :: sz

    sz = decomp%zsz(1)*decomp%zsz(2)*decomp%zsz(3)
    a1_p = fftw_alloc_complex(sz)
    call c_f_pointer(a1_p,a1,[decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)])

#ifdef DOUBLE_PREC
    plan1 =  fftw_plan_many_dft(1, decomp%zsz(3), &
         decomp%zsz(1)*decomp%zsz(2), a1, decomp%zsz(3), &
         decomp%zsz(1)*decomp%zsz(2), 1, a1, decomp%zsz(3), &
         decomp%zsz(1)*decomp%zsz(2), 1, isign, plan_type)
#else
    plan1 = fftwf_plan_many_dft(1, decomp%zsz(3), &
         decomp%zsz(1)*decomp%zsz(2), a1, decomp%zsz(3), &
         decomp%zsz(1)*decomp%zsz(2), 1, a1, decomp%zsz(3), &
         decomp%zsz(1)*decomp%zsz(2), 1, isign, plan_type)
#endif

    call fftw_free(a1_p)

    return
  end subroutine c2c_1m_z_plan

  ! Return a FFTW3 plan for multiple 1D r2c FFTs in X direction
  subroutine r2c_1m_x_plan(plan1, decomp_ph, decomp_sp)

    implicit none

    type(C_PTR) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp

    real(mytype), pointer :: a1(:,:,:)
    complex(mytype), pointer :: a2(:,:,:)
    type(C_PTR) :: a1_p, a2_p
    integer(C_SIZE_T) :: sz

    sz = decomp_ph%xsz(1)*decomp_ph%xsz(2)*decomp_ph%xsz(3)
    a1_p = fftw_alloc_real(sz)
    call c_f_pointer(a1_p,a1, &
         [decomp_ph%xsz(1),decomp_ph%xsz(2),decomp_ph%xsz(3)])
    sz = decomp_sp%xsz(1)*decomp_sp%xsz(2)*decomp_sp%xsz(3)
    a2_p = fftw_alloc_complex(sz)
    call c_f_pointer(a2_p,a2, &
         [decomp_sp%xsz(1),decomp_sp%xsz(2),decomp_sp%xsz(3)])

#ifdef DOUBLE_PREC
    plan1 =  fftw_plan_many_dft_r2c(1, decomp_ph%xsz(1), &
         decomp_ph%xsz(2)*decomp_ph%xsz(3), a1, decomp_ph%xsz(1), 1, &
         decomp_ph%xsz(1), a2, decomp_sp%xsz(1), 1, decomp_sp%xsz(1), &
         plan_type)
#else
    plan1 = fftwf_plan_many_dft_r2c(1, decomp_ph%xsz(1), &
         decomp_ph%xsz(2)*decomp_ph%xsz(3), a1, decomp_ph%xsz(1), 1, &
         decomp_ph%xsz(1), a2, decomp_sp%xsz(1), 1, decomp_sp%xsz(1), &
         plan_type)
#endif

    call fftw_free(a1_p)
    call fftw_free(a2_p)

    return
  end subroutine r2c_1m_x_plan

  ! Return a FFTW3 plan for multiple 1D c2r FFTs in X direction
  subroutine c2r_1m_x_plan(plan1, decomp_sp, decomp_ph)

    implicit none

    type(C_PTR) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph

    complex(mytype), pointer :: a1(:,:,:)
    real(mytype), pointer :: a2(:,:,:)
    type(C_PTR) :: a1_p, a2_p
    integer(C_SIZE_T) :: sz

    sz = decomp_sp%xsz(1)*decomp_sp%xsz(2)*decomp_sp%xsz(3)
    a1_p = fftw_alloc_complex(sz)
    call c_f_pointer(a1_p,a1, &
         [decomp_sp%xsz(1),decomp_sp%xsz(2),decomp_sp%xsz(3)])
    sz = decomp_ph%xsz(1)*decomp_ph%xsz(2)*decomp_ph%xsz(3)
    a2_p = fftw_alloc_real(sz)
    call c_f_pointer(a2_p,a2, &
         [decomp_ph%xsz(1),decomp_ph%xsz(2),decomp_ph%xsz(3)])

#ifdef DOUBLE_PREC
    plan1 =  fftw_plan_many_dft_c2r(1, decomp_ph%xsz(1), &
         decomp_ph%xsz(2)*decomp_ph%xsz(3), a1, decomp_sp%xsz(1), 1, &
         decomp_sp%xsz(1), a2, decomp_ph%xsz(1), 1, decomp_ph%xsz(1), &
         plan_type)
#else
    plan1 = fftwf_plan_many_dft_c2r(1, decomp_ph%xsz(1), &
         decomp_ph%xsz(2)*decomp_ph%xsz(3), a1, decomp_sp%xsz(1), 1, &
         decomp_sp%xsz(1), a2, decomp_ph%xsz(1), 1, decomp_ph%xsz(1), &
         plan_type)
#endif

    call fftw_free(a1_p)
    call fftw_free(a2_p)

    return
  end subroutine c2r_1m_x_plan

  ! Return a FFTW3 plan for multiple 1D r2c FFTs in Z direction
  subroutine r2c_1m_z_plan(plan1, decomp_ph, decomp_sp)

    implicit none

    type(C_PTR) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp

    real(mytype), pointer :: a1(:,:,:)
    complex(mytype), pointer :: a2(:,:,:)
    type(C_PTR) :: a1_p, a2_p
    integer(C_SIZE_T) :: sz

    sz = decomp_ph%zsz(1)*decomp_ph%zsz(2)*decomp_ph%zsz(3)
    a1_p = fftw_alloc_real(sz)
    call c_f_pointer(a1_p,a1, &
         [decomp_ph%zsz(1),decomp_ph%zsz(2),decomp_ph%zsz(3)])
    sz = decomp_sp%zsz(1)*decomp_sp%zsz(2)*decomp_sp%zsz(3)
    a2_p = fftw_alloc_complex(sz)
    call c_f_pointer(a2_p,a2, &
         [decomp_sp%zsz(1),decomp_sp%zsz(2),decomp_sp%zsz(3)])

#ifdef DOUBLE_PREC
    plan1 =  fftw_plan_many_dft_r2c(1, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), a1, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), 1, a2, decomp_sp%zsz(3), &
         decomp_sp%zsz(1)*decomp_sp%zsz(2), 1, plan_type)
#else
    plan1 = fftwf_plan_many_dft_r2c(1, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), a1, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), 1, a2, decomp_sp%zsz(3), &
         decomp_sp%zsz(1)*decomp_sp%zsz(2), 1, plan_type)
#endif

    call fftw_free(a1_p)
    call fftw_free(a2_p)

    return
  end subroutine r2c_1m_z_plan

  ! Return a FFTW3 plan for multiple 1D c2r FFTs in Z direction
  subroutine c2r_1m_z_plan(plan1, decomp_sp, decomp_ph)

    implicit none

    type(C_PTR) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph

    complex(mytype), pointer :: a1(:,:,:)
    real(mytype), pointer :: a2(:,:,:)
    type(C_PTR) :: a1_p, a2_p
    integer(C_SIZE_T) :: sz

    sz = decomp_sp%zsz(1)*decomp_sp%zsz(2)*decomp_sp%zsz(3)
    a1_p = fftw_alloc_complex(sz)
    call c_f_pointer(a1_p,a1, &
         [decomp_sp%zsz(1),decomp_sp%zsz(2),decomp_sp%zsz(3)])
    sz = decomp_ph%zsz(1)*decomp_ph%zsz(2)*decomp_ph%zsz(3)
    a2_p = fftw_alloc_real(sz)
    call c_f_pointer(a2_p,a2, &
         [decomp_ph%zsz(1),decomp_ph%zsz(2),decomp_ph%zsz(3)])

#ifdef DOUBLE_PREC
    plan1 =  fftw_plan_many_dft_c2r(1, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), a1, decomp_sp%zsz(3), &
         decomp_sp%zsz(1)*decomp_sp%zsz(2), 1, a2, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), 1, plan_type)
#else
    plan1 = fftwf_plan_many_dft_c2r(1, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), a1, decomp_sp%zsz(3), &
         decomp_sp%zsz(1)*decomp_sp%zsz(2), 1, a2, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), 1, plan_type)
#endif

    call fftw_free(a1_p)
    call fftw_free(a2_p)

    return
  end subroutine c2r_1m_z_plan


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time initialisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_fft_engine(Igrid)

    implicit none

    integer, intent(IN) :: Igrid

    if (nrank==0) then
       write(*,*) ' '
       write(*,*) '***** Using the FFTW (F2003 interface) engine *****'
       write(*,*) ' '
    end if

    if (format == PHYSICAL_IN_X) then

       ! For C2C transforms
       call c2c_1m_x_plan(FFT_multigrid(Igrid)%plan(-1,1), ph, FFTW_FORWARD )
       call c2c_1m_y_plan(FFT_multigrid(Igrid)%plan(-1,2), ph, FFTW_FORWARD )
       call c2c_1m_z_plan(FFT_multigrid(Igrid)%plan(-1,3), ph, FFTW_FORWARD )
       call c2c_1m_z_plan(FFT_multigrid(Igrid)%plan( 1,3), ph, FFTW_BACKWARD)
       call c2c_1m_y_plan(FFT_multigrid(Igrid)%plan( 1,2), ph, FFTW_BACKWARD)
       call c2c_1m_x_plan(FFT_multigrid(Igrid)%plan( 1,1), ph, FFTW_BACKWARD)
       
       ! For R2C/C2R tranforms
       call r2c_1m_x_plan(FFT_multigrid(Igrid)%plan(0,1), ph, sp)
       call c2c_1m_y_plan(FFT_multigrid(Igrid)%plan(0,2), sp, FFTW_FORWARD )
       call c2c_1m_z_plan(FFT_multigrid(Igrid)%plan(0,3), sp, FFTW_FORWARD )
       call c2c_1m_z_plan(FFT_multigrid(Igrid)%plan(2,3), sp, FFTW_BACKWARD)
       call c2c_1m_y_plan(FFT_multigrid(Igrid)%plan(2,2), sp, FFTW_BACKWARD)
       call c2r_1m_x_plan(FFT_multigrid(Igrid)%plan(2,1), sp, ph)

    else if (format == PHYSICAL_IN_Z) then

       ! For C2C transforms
       call c2c_1m_z_plan(FFT_multigrid(Igrid)%plan(-1,3), ph, FFTW_FORWARD )
       call c2c_1m_y_plan(FFT_multigrid(Igrid)%plan(-1,2), ph, FFTW_FORWARD ) 
       call c2c_1m_x_plan(FFT_multigrid(Igrid)%plan(-1,1), ph, FFTW_FORWARD )
       call c2c_1m_x_plan(FFT_multigrid(Igrid)%plan( 1,1), ph, FFTW_BACKWARD)
       call c2c_1m_y_plan(FFT_multigrid(Igrid)%plan( 1,2), ph, FFTW_BACKWARD)
       call c2c_1m_z_plan(FFT_multigrid(Igrid)%plan( 1,3), ph, FFTW_BACKWARD)
       
       ! For R2C/C2R tranforms
       call r2c_1m_z_plan(FFT_multigrid(Igrid)%plan(0,3), ph, sp)
       call c2c_1m_y_plan(FFT_multigrid(Igrid)%plan(0,2), sp, FFTW_FORWARD )
       call c2c_1m_x_plan(FFT_multigrid(Igrid)%plan(0,1), sp, FFTW_FORWARD )
       call c2c_1m_x_plan(FFT_multigrid(Igrid)%plan(2,1), sp, FFTW_BACKWARD)
       call c2c_1m_y_plan(FFT_multigrid(Igrid)%plan(2,2), sp, FFTW_BACKWARD)
       call c2r_1m_z_plan(FFT_multigrid(Igrid)%plan(2,3), sp, ph)
       
    end if

    return
  end subroutine init_fft_engine


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time finalisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine finalize_fft_engine

    implicit none

    integer :: Igrid,i,j
    
    do Igrid = 1, size(FFT_multigrid)
    do j=1,3
       do i=-1,2
#ifdef DOUBLE_PREC
          call fftw_destroy_plan(FFT_multigrid(Igrid)%plan(i,j))
#else
          call fftwf_destroy_plan(FFT_multigrid(Igrid)%plan(i,j))
#endif
       end do
    end do
    end do

    return
  end subroutine finalize_fft_engine


  ! Following routines calculate multiple one-dimensional FFTs to form 
  ! the basis of three-dimensional FFTs.

  ! c2c transform, multiple 1D FFTs in x direction
  subroutine c2c_1m_x(inout, isign, plan1)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    type(C_PTR) :: plan1

#ifdef DOUBLE_PREC
    call fftw_execute_dft(plan1, inout, inout)
#else
    call fftwf_execute_dft(plan1, inout, inout)
#endif

    return
  end subroutine c2c_1m_x


  ! c2c transform, multiple 1D FFTs in y direction
  subroutine c2c_1m_y(inout, isign, plan1)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    type(C_PTR) :: plan1

    integer :: k, s3

    s3 = size(inout,3)

    do k=1,s3  ! transform on one Z-plane at a time
#ifdef DOUBLE_PREC
       call fftw_execute_dft(plan1, inout(:,:,k), inout(:,:,k))
#else
       call fftwf_execute_dft(plan1, inout(:,:,k), inout(:,:,k))
#endif
    end do

    return
  end subroutine c2c_1m_y

  ! c2c transform, multiple 1D FFTs in z direction
  subroutine c2c_1m_z(inout, isign, plan1)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    type(C_PTR) :: plan1

#ifdef DOUBLE_PREC
       call fftw_execute_dft(plan1, inout, inout)
#else
       call fftwf_execute_dft(plan1, inout, inout)
#endif

    return
  end subroutine c2c_1m_z

  ! r2c transform, multiple 1D FFTs in x direction
  subroutine r2c_1m_x(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(INOUT)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output

#ifdef DOUBLE_PREC
    call fftw_execute_dft_r2c(plan(0,1), input, output)
#else
    call fftwf_execute_dft_r2c(plan(0,1), input, output)
#endif    

    return

  end subroutine r2c_1m_x

  ! r2c transform, multiple 1D FFTs in z direction
  subroutine r2c_1m_z(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(INOUT)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output

#ifdef DOUBLE_PREC
    call fftw_execute_dft_r2c(plan(0,3), input, output)
#else
    call fftwf_execute_dft_r2c(plan(0,3), input, output)
#endif

    return

  end subroutine r2c_1m_z

  ! c2r transform, multiple 1D FFTs in x direction
  subroutine c2r_1m_x(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT)  ::  input
    real(mytype), dimension(:,:,:), intent(OUT) :: output

#ifdef DOUBLE_PREC
    call fftw_execute_dft_c2r(plan(2,1), input, output)
#else
    call fftwf_execute_dft_c2r(plan(2,1), input, output)
#endif

    return

  end subroutine c2r_1m_x

  ! c2r transform, multiple 1D FFTs in z direction
  subroutine c2r_1m_z(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: input
    real(mytype), dimension(:,:,:), intent(OUT) :: output

#ifdef DOUBLE_PREC
    call fftw_execute_dft_c2r(plan(2,3), input, output)
#else
    call fftwf_execute_dft_c2r(plan(2,3), input, output)
#endif    

    return

  end subroutine c2r_1m_z



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D FFT - complex to complex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_c2c(in, out, isign)
    
    implicit none
    
    complex(mytype), dimension(:,:,:), intent(INOUT) :: in
    complex(mytype), dimension(:,:,:), intent(OUT) :: out
    integer, intent(IN) :: isign

#ifndef OVERWRITE
    complex(mytype), pointer :: wk1(:,:,:)
    integer(C_SIZE_T) :: sz
    type(C_PTR) :: wk1_p
#endif

    if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_FORWARD .OR.  &
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_BACKWARD) then
       
       ! ===== 1D FFTs in X =====
#ifdef OVERWRITE
       call c2c_1m_x(in,isign,plan(isign,1))
#else
       sz = ph%xsz(1)*ph%xsz(2)*ph%xsz(3)
       wk1_p = fftw_alloc_complex(sz)
       call c_f_pointer(wk1_p, wk1, [ph%xsz(1),ph%xsz(2),ph%xsz(3)])
       wk1 = in
       call c2c_1m_x(wk1,isign,plan(isign,1))
#endif

       ! ===== Swap X --> Y; 1D FFTs in Y =====

       if (dims(1)>1) then
#ifdef OVERWRITE
          call transpose_x_to_y(in,wk2_c2c,ph)
#else
          call transpose_x_to_y(wk1,wk2_c2c,ph)
#endif
          call c2c_1m_y(wk2_c2c,isign,plan(isign,2))
       else
#ifdef OVERWRITE
          call c2c_1m_y(in,isign,plan(isign,2))
#else
          call c2c_1m_y(wk1,isign,plan(isign,2))
#endif
       end if

       ! ===== Swap Y --> Z; 1D FFTs in Z =====
       if (dims(1)>1) then
          call transpose_y_to_z(wk2_c2c,out,ph)
       else
#ifdef OVERWRITE
          call transpose_y_to_z(in,out,ph)
#else
          call transpose_y_to_z(wk1,out,ph)
#endif
       end if
       call c2c_1m_z(out,isign,plan(isign,3))

    else if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_BACKWARD &
         .OR. & 
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_FORWARD) then

       ! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
       call c2c_1m_z(in,isign,plan(isign,3))
#else
       sz = ph%zsz(1)*ph%zsz(2)*ph%zsz(3)
       wk1_p = fftw_alloc_complex(sz)
       call c_f_pointer(wk1_p, wk1, [ph%zsz(1),ph%zsz(2),ph%zsz(3)])
       wk1 = in
       call c2c_1m_z(wk1,isign,plan(isign,3))
#endif

       ! ===== Swap Z --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
#ifdef OVERWRITE
          call transpose_z_to_y(in,wk2_c2c,ph)
#else
          call transpose_z_to_y(wk1,wk2_c2c,ph)
#endif
          call c2c_1m_y(wk2_c2c,isign,plan(isign,2))
       else  ! out==wk2_c2c if 1D decomposition
#ifdef OVERWRITE
          call transpose_z_to_y(in,out,ph)
#else
          call transpose_z_to_y(wk1,out,ph)
#endif
          call c2c_1m_y(out,isign,plan(isign,2))
       end if

       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          call transpose_y_to_x(wk2_c2c,out,ph)
       end if
       call c2c_1m_x(out,isign,plan(isign,1))
       
    end if

#ifndef OVERWRITE
    call fftw_free(wk1_p)
#endif

    return
  end subroutine fft_3d_c2c

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D forward FFT - real to complex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_r2c(in_r, out_c)
    
    implicit none
    
    real(mytype), dimension(:,:,:), intent(INOUT) :: in_r
    complex(mytype), dimension(:,:,:), intent(OUT) :: out_c

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in X =====
       call r2c_1m_x(in_r,wk13)

       ! ===== Swap X --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
          call transpose_x_to_y(wk13,wk2_r2c,sp)
          call c2c_1m_y(wk2_r2c,-1,plan(0,2))
       else
          call c2c_1m_y(wk13,-1,plan(0,2))
       end if

       ! ===== Swap Y --> Z; 1D FFTs in Z =====
       if (dims(1)>1) then
          call transpose_y_to_z(wk2_r2c,out_c,sp)
       else
          call transpose_y_to_z(wk13,out_c,sp)
       end if
       call c2c_1m_z(out_c,-1,plan(0,3))
                
    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in Z =====
       call r2c_1m_z(in_r,wk13)

       ! ===== Swap Z --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
          call transpose_z_to_y(wk13,wk2_r2c,sp)
          call c2c_1m_y(wk2_r2c,-1,plan(0,2))
       else  ! out_c==wk2_r2c if 1D decomposition
          call transpose_z_to_y(wk13,out_c,sp)
          call c2c_1m_y(out_c,-1,plan(0,2))
       end if

       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          call transpose_y_to_x(wk2_r2c,out_c,sp)
       end if
       call c2c_1m_x(out_c,-1,plan(0,1))

    end if
    
    return
  end subroutine fft_3d_r2c
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D inverse FFT - complex to real
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_c2r(in_c, out_r)
    
    implicit none
    
    complex(mytype), dimension(:,:,:), intent(INOUT) :: in_c
    real(mytype), dimension(:,:,:), intent(OUT) :: out_r

#ifndef OVERWRITE
    complex(mytype), pointer :: wk1(:,:,:)
    integer(C_SIZE_T) :: sz
    type(C_PTR) :: wk1_p
#endif

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
       call c2c_1m_z(in_c,1,plan(2,3))       
#else
       sz = sp%zsz(1)*sp%zsz(2)*sp%zsz(3)
       wk1_p = fftw_alloc_complex(sz)
       call c_f_pointer(wk1_p, wk1, [sp%zsz(1),sp%zsz(2),sp%zsz(3)])
       wk1 = in_c
       call c2c_1m_z(wk1,1,plan(2,3))
#endif

       ! ===== Swap Z --> Y; 1D FFTs in Y =====
#ifdef OVERWRITE
       call transpose_z_to_y(in_c,wk2_r2c,sp)
#else
       call transpose_z_to_y(wk1,wk2_r2c,sp)
#endif
       call c2c_1m_y(wk2_r2c,1,plan(2,2))

       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          call transpose_y_to_x(wk2_r2c,wk13,sp)
          call c2r_1m_x(wk13,out_r)
       else
          call c2r_1m_x(wk2_r2c,out_r)
       end if

    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in X =====
#ifdef OVERWRITE
       call c2c_1m_x(in_c,1,plan(2,1))
#else
       sz = sp%xsz(1)*sp%xsz(2)*sp%xsz(3)
       wk1_p = fftw_alloc_complex(sz)
       call c_f_pointer(wk1_p, wk1, [sp%xsz(1),sp%xsz(2),sp%xsz(3)])
       wk1 = in_c
       call c2c_1m_x(wk1,1,plan(2,1))
#endif

       ! ===== Swap X --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
#ifdef OVERWRITE
          call transpose_x_to_y(in_c,wk2_r2c,sp)
#else
          call transpose_x_to_y(wk1,wk2_r2c,sp)
#endif
          call c2c_1m_y(wk2_r2c,1,plan(2,2))
       else  ! in_c==wk2_r2c if 1D decomposition
#ifdef OVERWRITE
          call c2c_1m_y(in_c,1,plan(2,2))
#else
          call c2c_1m_y(wk1,1,plan(2,2))
#endif
       end if

       ! ===== Swap Y --> Z; 1D FFTs in Z =====
       if (dims(1)>1) then
          call transpose_y_to_z(wk2_r2c,wk13,sp)
       else
#ifdef OVERWRITE
          call transpose_y_to_z(in_c,wk13,sp)
#else
          call transpose_y_to_z(wk1,wk13,sp)
#endif
       end if
       call c2r_1m_z(wk13,out_r)

    end if

#ifndef OVERWRITE
    call fftw_free(wk1_p)
#endif

    return
  end subroutine fft_3d_c2r

  
end module decomp_2d_fft
