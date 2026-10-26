program test_field_divgrad
!> CHARGEMENT DES DIFFERENTS MODULES

  use ISO_FORTRAN_ENV

  use MPI
#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype,PHYSICAL_IN_X
  use decomp_2d_mpi,       only : nrank
  use decomp_2d,           only : xsize, xstart,xend,&
                                  decomp_2d_init, decomp_2d_finalize
#else
  use decomp_2d,      only : mytype,&                           ! TYPE
                             nrank, &                           ! Variable
                             xsize, xstart,xend,&               ! Pointer Variable
                             decomp_2d_init, decomp_2d_finalize ! Functions
#endif  

  
  use decomp_2d_fft,  only : decomp_2d_fft_init, decomp_2d_fft_get_size,decomp_2d_fft_finalize ! Functions
  use amitex_mod,     only : TEMPfield32 ! Pointer Variables
  use error_mod,      only : amitex_abort ! Function
  use field_mod,      only : TIMES_FIELD, & ! Type
                             grid, fft_size, fft_start, fft_end,&   ! Pointer Variables
                             field_div, field_grad, field_divF, field_gradF,& ! Functions
                             field_ifft, field_fft, &
                             initGrid, initFreq, initFreq_2, gene_cartesian_coord 
  use amitex_lib
  use simu_mod

!!==============================================================================
!!                                                                  DECLARATIONS
!!==============================================================================
  implicit none
!!------------------------------------------------------------------------------
!>                                                       DIMENSIONS DE LA GRILLE   
!>                                                                    

  integer,parameter               :: nx=63, ny=18, nz=127        !< dimensions de la cellule (voxel)
  real(mytype),parameter          :: Lx=10, Ly=20, Lz=30         !< dimensions de la cellule (taille)
  real(mytype)                    :: dx, dy, dz,x0=0.,y0=0.,z0=0.     


!!------------------------------------------------------------------------------
!>                                                                  PARALLELISME
!>                                     + pointer target for using amitex modules

  integer                      :: p_row=1, p_col=14           !< nb. de lignes, de colonnes pour la decomposition de 2decomp
  logical,dimension(3)         :: periodic_bc=.true.
  integer                      :: ierror                     !< erreur relative au fonction MPI


!!------------------------------------------------------------------------------
!>                                                           tableaux paralleles

  !> dans l'espace reel : tableaux (ntot)
  real(mytype),allocatable,dimension(:,:)                  :: V01, V03, V01f, V01r, V06, V03r, V03f,V09, V06r, V06f, V09r,V09f
  real(mytype),allocatable,dimension(:,:,:)                :: X,Y,Z
  
  real(mytype),allocatable,dimension(:,:,:,:)              :: FREQ_
  complex(mytype),allocatable,dimension(:,:,:,:)           :: V1F,V3F,V6F,V9F,FREQ_2_
  real(mytype)                                             :: max_err
  integer                                                  :: i,Icase
  integer, dimension(3)                                    :: BC0
  double precision                                         :: t0,t1
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

  if(nrank==0)print *, "Dimensions : ",nx,"x",ny,"x",nz, p_row,p_col

!!------------------------------------------------------------------------------
!>                          		          FIN INITIALISATION AMITEX_LIB 

  dx = Lx/nx
  dy = Ly/ny
  dz = Lz/nz
  call init_amitex_lib(nx,ny,nz,dx,dy,dz,x0,y0,z0,ficlog="test_field_divgrad.log")!,filter_type="no_filter",filter_radius=0._8)


!!------------------------------------------------------------------------------
!>                                  INITIALISATION FREQUENCES (LOCALES AU PROGR)
!>                                           FREQUENCES MODIFIEES EF-IR (hexa,1)

  allocate(FREQ_2_(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)  
  allocate(FREQ_(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),3),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  
  call initFREQ(Freq_,'no_filter',0._mytype)
  call initFREQ_2(Freq_,FREQ_2_)
  call initFREQ(Freq_,'hexa',1._mytype)
 
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


  allocate(V09(xsize(1)*xsize(2)*xsize(3),9),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(V09r(xsize(1)*xsize(2)*xsize(3),9),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(V09f(xsize(1)*xsize(2)*xsize(3),9),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  
  allocate(V06(xsize(1)*xsize(2)*xsize(3),6),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(V06r(xsize(1)*xsize(2)*xsize(3),6),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(V06f(xsize(1)*xsize(2)*xsize(3),6),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)


  allocate(V03(xsize(1)*xsize(2)*xsize(3),3),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(V03r(xsize(1)*xsize(2)*xsize(3),3),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(V03f(xsize(1)*xsize(2)*xsize(3),3),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  
  allocate(V01(xsize(1)*xsize(2)*xsize(3),1),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(V01f(xsize(1)*xsize(2)*xsize(3),1),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(V01r(xsize(1)*xsize(2)*xsize(3),1),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  
  allocate(X(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(Y(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  allocate(Z(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)


  allocate(V1F(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  
  allocate(V3F(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),3),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  
  allocate(V6F(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),6),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)
  
  allocate(V9F(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),9),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant ",2)


! definition des coordonnes de la grille (pour validation filtre Gauss)
  call gene_cartesian_coord(X,Y,Z,dx,dy,dz)
  X = X - dx / 2._mytype
  Y = Y - dx / 2._mytype
  Z = Z - dz / 2._mytype

DO Icase = 1,2

  select case(Icase)
  ! On assigne V0i un champ aleatoire
  case(1)
    if (nrank==0) print * , "******************************************* PBC + RANDOM FIELD"
    call random_number(V01)
    call random_number(V03)
    call random_number(V06)
    call random_number(V09)
    BC0 = (/ 0,0,0/)
  ! Champ linear
  case(2)
    if (nrank==0) print * , "******************************************* LIN. EXTRAP. BC + LINEAR FIELD"
    V01(:,1) = 1._mytype*reshape(X,(/xsize(1)*xsize(2)*xsize(3)/)) +&
             2._mytype*reshape(Y,(/xsize(1)*xsize(2)*xsize(3)/)) + &
             3._mytype*reshape(Z,(/xsize(1)*xsize(2)*xsize(3)/)) +7._mytype
    do i=1,3; V03(:,i) = V01(:,1); end do
    do i=1,6; V06(:,i) = V01(:,1); end do
    do i=1,9; V09(:,i) = V01(:,1); end do
    BC0 = (/ 1,1,1/)
  
  end select 
!  call print_field_vtk(V0,"V0.vtk","V0") 

!!*************************************************************************  
!!************************************************************************* DIVERGENCE
!!************************************************************************* 
  if (nrank==0) print * , "********************************* DIVERGENCE"

!!------------------------------------------------------------------------------
!>                                                         Div(V0) - CAS VECTEUR
!>                                                             

  t1 = MPI_WTIME()
  
  !--------------------------------------------------CENTRES -> NOEUDS
  !--------------------------- ESPACE REEL  
  call field_div(V03,V01r,3,1,d=(/grid%dx,grid%dy,grid%dz /),N=(/grid%nx,grid%ny,grid%nz /),dtype=-1,BC=BC0)

  !--------------------------- ESPACE FOURIER
  call field_fft(V03,V3F,3,1)
  call field_divF(V3F,V1F,3,1,FREQ_,FREQ_2_,dtype=-1)
  call field_ifft(V01f,V1F,1,1)
  V01f = V01f / grid%ntot
  if (Icase==2) then
     V01f(:,1)=6._mytype
  end if  
 
  max_err = huge(max_err)
  call diff_fields(V01r,V01f,max_err,2,1,grid%ntot)
  if (nrank==0) print * , "----- VECTEUR Centres->Noeuds : sum(abs(V01r - V01f))/ntot = ", max_err
!if (Icase==2) then
!  print *, "diff div_1 : ", nrank, maxval(V01r(:,1)) - minval(V01r(:,1)), minval(V01r(:,1)), maxval(V01r(:,1)) 
!end if  
  
  !--------------------------------------------------NOEUDS -> CENTRES
  !--------------------------- ESPACE REEL  
  call field_div(V03,V01r,3,1,d=(/grid%dx,grid%dy,grid%dz /),N=(/grid%nx,grid%ny,grid%nz /),dtype=1,BC=BC0)

  !--------------------------- ESPACE FOURIER
  call field_fft(V03,V3F,3,1)
  call field_divF(V3F,V1F,3,1,FREQ_,FREQ_2_,dtype=1)
  call field_ifft(V01f,V1F,1,1)
  V01f = V01f / grid%ntot
  if (Icase==2) then
     V01f(:,1)=6._mytype
  end if  
 
  max_err = huge(max_err)
  call diff_fields(V01r,V01f,max_err,2,1,grid%ntot)
  if (nrank==0) print * , "----- VECTEUR Noeuds->Centres : sum(abs(V01r - V01f))/ntot = ", max_err
 

!!------------------------------------------------------------------------------
!>                                           Div(V0) - CAS TENSEUR 2 SYMMETRIQUE
!>                                                             

  !--------------------------------------------------CENTRES -> NOEUDS
  !--------------------------- ESPACE REEL  
  call field_div(V06,V03r,6,3,d=(/grid%dx,grid%dy,grid%dz /),N=(/grid%nx,grid%ny,grid%nz /),&
                 dtype=-1,voigt_convention="strain",BC=BC0)

  !--------------------------- ESPACE FOURIER
  call field_fft(V06,V6F,6,1)
  call field_divF(V6F,V3F,6,3,FREQ_,FREQ_2_,dtype=-1,voigt_convention="strain")
  call field_ifft(V03f,V3F,3,1)
  V03f = V03f / grid%ntot
  if (Icase==2) then
     V03f(:,1)=3.5_mytype;V03f(:,2)=4._mytype;V03f(:,3)=4.5_mytype;
  end if  
  
  max_err = huge(max_err)
  call diff_fields(V03r,V03f,max_err,2,3,grid%ntot)
  if (nrank==0) print * , "----- TENSEUR2 SYM Centres->Noeuds : sum(abs(V03r - V03f))/ntot = ", max_err
!if (Icase==2) then
!  print *, "diff div_1 : ", nrank, maxval(V03r(:,1)) - minval(V03r(:,1)), maxval(V03r(:,1)) 
!  print *, "diff div_2 : ",nrank, maxval(V03r(:,2)) - minval(V03r(:,2)), maxval(V03r(:,2)) 
!  print *, "diff div_3 : ",nrank, maxval(V03r(:,3)) - minval(V03r(:,3)), maxval(V03r(:,3)) 
!end if  
  
  !--------------------------------------------------NOEUDS -> CENTRES
  V03r=0
  V03f=0
  !--------------------------- ESPACE REEL  
  call field_div(V06,V03r,6,3,d=(/grid%dx,grid%dy,grid%dz /),N=(/grid%nx,grid%ny,grid%nz /),&
                 dtype=1,voigt_convention="strain",BC=BC0)

  !--------------------------- ESPACE FOURIER
  call field_fft(V06,V6F,6,1)
  call field_divF(V6F,V3F,6,3,FREQ_,FREQ_2_,dtype=1,voigt_convention="strain")
  call field_ifft(V03f,V3F,3,1)
  V03f = V03f / grid%ntot
  if (Icase==2) then
     V03f(:,1)=3.5_mytype;V03f(:,2)=4._mytype;V03f(:,3)=4.5_mytype;
  end if  
  
  max_err = huge(max_err)
  call diff_fields(V03r,V03f,max_err,2,3,grid%ntot)
  if (nrank==0) print * , "----- TENSEUR2 SYM Noeuds->Centres : sum(abs(V03r - V03f))/ntot = ", max_err

!!------------------------------------------------------------------------------
!>                                                       Div(V0) - CAS TENSEUR 2
!>                                                             

  !--------------------------------------------------CENTRES -> NOEUDS
  !--------------------------- ESPACE REEL  
  call field_div(V09,V03r,9,3,d=(/grid%dx,grid%dy,grid%dz /),N=(/grid%nx,grid%ny,grid%nz /),dtype=-1,BC=BC0)

  !--------------------------- ESPACE FOURIER
  call field_fft(V09,V9F,9,1)
  call field_divF(V9F,V3F,9,3,FREQ_,FREQ_2_,dtype=-1)
  call field_ifft(V03f,V3F,3,1)
  V03f = V03f / grid%ntot
  if (Icase==2) then
     V03f(:,1)=6._mytype;V03f(:,2)=6._mytype;V03f(:,3)=6._mytype;
  end if  
  
  max_err = huge(max_err)
  call diff_fields(V03r,V03f,max_err,2,3,grid%ntot)
  if (nrank==0) print * , "----- TENSEUR 2 Centres->Noeuds : sum(abs(V03r - V03f))/ntot = ", max_err
  
  !--------------------------------------------------NOEUDS -> CENTRES
  !--------------------------- ESPACE REEL  
  call field_div(V09,V03r,9,3,d=(/grid%dx,grid%dy,grid%dz /),N=(/grid%nx,grid%ny,grid%nz /),dtype=1,BC=BC0)

  !--------------------------- ESPACE FOURIER
  call field_fft(V09,V9F,9,1)
  call field_divF(V9F,V3F,9,3,FREQ_,FREQ_2_,dtype=1)
  call field_ifft(V03f,V3F,3,1)
  V03f = V03f / grid%ntot
  if (Icase==2) then
     V03f(:,1)=6._mytype;V03f(:,2)=6._mytype;V03f(:,3)=6._mytype;
  end if  
  
  max_err = huge(max_err)
  call diff_fields(V03r,V03f,max_err,2,3,grid%ntot)
  if (nrank==0) print * , "----- TENSEUR 2 Noeuds->Centres : sum(abs(V03r - V03f))/ntot = ", max_err

!print *, "V03r(1:3,1:3)-V03f(1:5,1) =", V03r(1:4,1:3)-V03f(1:4,1:3)

!!*************************************************************************  
!!************************************************************************* GRADIENT
!!************************************************************************* 
  if (nrank==0) print * , "********************************* GRADIENT"

!!------------------------------------------------------------------------------
!>                                                       Grad(V0) - CAS SCALAIRE
!>                                                             

  !--------------------------------------------------CENTRES -> NOEUDS
  !--------------------------- ESPACE REEL  
  call field_grad(V01,V03r,1,3,(/grid%dx,grid%dy,grid%dz /),(/grid%nx,grid%ny,grid%nz /),dtype=-1,BC=BC0)
  
  !--------------------------- ESPACE FOURIER
  call field_fft(V01,V1F,1,1)
  call field_gradF(V1F,V3F,1,3,FREQ_,FREQ_2_,-1)
  call field_ifft(V03f,V3F,3,1)
  V03f = V03f / grid%ntot
  if (Icase==2) then
     V03f(:,1)=1._mytype;V03f(:,2)=2._mytype;V03f(:,3)=3._mytype 
  end if  
  max_err = huge(max_err)
  call diff_fields(V03r,V03f,max_err,2,3,grid%ntot)
  if (nrank==0) print * , "----- SCALAIRE Centres->Noeuds : sum(abs(V03r - V03f))/ntot = ", max_err

! Should be 0 if VOr is a linear field  
!if (Icase==2) then
!  print *, "diff grad_1 : ", nrank, maxval(V03r(:,1)) - minval(V03r(:,1)), maxval(V03r(:,1)) 
!  print *, "diff grad_2 : ",nrank, maxval(V03r(:,2)) - minval(V03r(:,2)), maxval(V03r(:,2)) 
!  print *, "diff grad_3 : ",nrank, maxval(V03r(:,3)) - minval(V03r(:,3)), maxval(V03r(:,3)) 
!end if  
  !--------------------------------------------------NOEUDS -> CENTRES
  !--------------------------- ESPACE REEL  
  call field_grad(V01,V03r,1,3,(/grid%dx,grid%dy,grid%dz /),(/grid%nx,grid%ny,grid%nz /),dtype=1,BC=BC0)
  
  !--------------------------- ESPACE FOURIER
  call field_fft(V01,V1F,1,1)
  call field_gradF(V1F,V3F,1,3,FREQ_,FREQ_2_,1)
  call field_ifft(V03f,V3F,3,1)
  V03f = V03f / grid%ntot
  if (Icase==2) then
     V03f(:,1)=1._mytype;V03f(:,2)=2._mytype;V03f(:,3)=3._mytype 
  end if  
    
  max_err = huge(max_err)
  call diff_fields(V03r,V03f,max_err,2,3,grid%ntot)
  if (nrank==0) print * , "----- SCALAIRE Noeuds->Centres : sum(abs(V03r - V03f))/ntot = ", max_err

!!------------------------------------------------------------------------------
!>                                                       Grad(V0) - CAS VECTEUR
!>                                                             OUT TENSEUR 2 SYM

  !--------------------------------------------------CENTRES -> NOEUDS
  !--------------------------- ESPACE REEL  
  call field_grad(V03,V06r,3,6,(/grid%dx,grid%dy,grid%dz /),(/grid%nx,grid%ny,grid%nz /),&
                  dtype=-1,voigt_convention="strain",BC=BC0)
  
  !--------------------------- ESPACE FOURIER
  call field_fft(V03,V3F,3,1)
  call field_gradF(V3F,V6F,3,6,FREQ_,FREQ_2_,-1,"strain")
  call field_ifft(V06f,V6F,6,1)
  V06f = V06f / grid%ntot
  if (Icase==2) then
     V06f(:,1)=1._mytype;V06f(:,2)=2._mytype;V06f(:,3)=3._mytype 
     V06f(:,4)=3._mytype;V06f(:,5)=4._mytype;V06f(:,6)=5._mytype 
  end if  
    
  max_err = huge(max_err)
  call diff_fields(V06r,V06f,max_err,2,6,grid%ntot)
  if (nrank==0) print * , "----- VECTEUR (out SYM) Centres->Noeuds : sum(abs(V06r - V06f))/ntot = ", max_err

!if (Icase==2) then
!  print *, "diff grad_1 : ", nrank, maxval(V06r(:,1)) - minval(V06r(:,1)), maxval(V06r(:,1)) 
!  print *, "diff grad_2 : ", nrank, maxval(V06r(:,2)) - minval(V06r(:,2)), maxval(V06r(:,2)) 
!  print *, "diff grad_3 : ", nrank, maxval(V06r(:,3)) - minval(V06r(:,3)), maxval(V06r(:,3)) 
!  print *, "diff grad_4 : ", nrank, maxval(V06r(:,4)) - minval(V06r(:,4)), maxval(V06r(:,4)) 
!  print *, "diff grad_5 : ", nrank, maxval(V06r(:,5)) - minval(V06r(:,5)), maxval(V06r(:,5)) 
!  print *, "diff grad_6 : ", nrank, maxval(V06r(:,6)) - minval(V06r(:,6)), maxval(V06r(:,6)) 
!end if  

  !--------------------------------------------------NOEUDS -> CENTRES
  !--------------------------- ESPACE REEL  
  call field_grad(V03,V06r,3,6,(/grid%dx,grid%dy,grid%dz /),(/grid%nx,grid%ny,grid%nz /),&
                  dtype=1,voigt_convention="strain",BC=BC0)
  
  !--------------------------- ESPACE FOURIER
  call field_fft(V03,V3F,3,1)
  call field_gradF(V3F,V6F,3,6,FREQ_,FREQ_2_,1,"strain")
  call field_ifft(V06f,V6F,6,1)
  V06f = V06f / grid%ntot
  if (Icase==2) then
     V06f(:,1)=1._mytype;V06f(:,2)=2._mytype;V06f(:,3)=3._mytype 
     V06f(:,4)=3._mytype;V06f(:,5)=4._mytype;V06f(:,6)=5._mytype 
  end if  
    
  max_err = huge(max_err)
  call diff_fields(V06r,V06f,max_err,2,6,grid%ntot)
  if (nrank==0) print * , "----- VECTEUR (out SYM) Noeuds->Centres : sum(abs(V06r - V06f))/ntot = ", max_err


!!------------------------------------------------------------------------------
!>                                                       Grad(V0) - CAS VECTEUR
!>                                                        OUT TENSEUR 2 NON SYM

  !--------------------------------------------------CENTRES -> NOEUDS
  !--------------------------- ESPACE REEL  
  call field_grad(V03,V09r,3,9,(/grid%dx,grid%dy,grid%dz /),(/grid%nx,grid%ny,grid%nz /),dtype=-1,BC=BC0)
  
  !--------------------------- ESPACE FOURIER
  call field_fft(V03,V3F,3,1)
  call field_gradF(V3F,V9F,3,9,FREQ_,FREQ_2_,-1)
  call field_ifft(V09f,V9F,9,1)
  V09f = V09f / grid%ntot
  if (Icase==2) then
     V09f(:,1)=1._mytype;V09f(:,2)=2._mytype;V09f(:,3)=3._mytype 
     V09f(:,4)=2._mytype;V09f(:,5)=3._mytype;V09f(:,6)=3._mytype 
     V09f(:,7)=1._mytype;V09f(:,8)=1._mytype;V09f(:,9)=2._mytype 
  end if  
    
  max_err = huge(max_err)
  call diff_fields(V09r,V09f,max_err,2,9,grid%ntot)
  if (nrank==0) print * , "----- VECTEUR Centres->Noeuds : sum(abs(V09r - V09f))/ntot = ", max_err

  !--------------------------------------------------NOEUDS CENTRES
  !--------------------------- ESPACE REEL  
  call field_grad(V03,V09r,3,9,(/grid%dx,grid%dy,grid%dz /),(/grid%nx,grid%ny,grid%nz /),dtype=1,BC=BC0)
  
  !--------------------------- ESPACE FOURIER
  call field_fft(V03,V3F,3,1)
  call field_gradF(V3F,V9F,3,9,FREQ_,FREQ_2_,1)
  call field_ifft(V09f,V9F,9,1)
  V09f = V09f / grid%ntot
  if (Icase==2) then
     V09f(:,1)=1._mytype;V09f(:,2)=2._mytype;V09f(:,3)=3._mytype 
     V09f(:,4)=2._mytype;V09f(:,5)=3._mytype;V09f(:,6)=3._mytype 
     V09f(:,7)=1._mytype;V09f(:,8)=1._mytype;V09f(:,9)=2._mytype 
  end if  
    
  max_err = huge(max_err)
  call diff_fields(V09r,V09f,max_err,2,9,grid%ntot)
  if (nrank==0) print * , "----- VECTEUR Noeuds->Centres : sum(abs(V09r - V09f))/ntot = ", max_err


END DO ! Tests case : random with BC=000 / linear with BC=111
 
  call decomp_2d_fft_finalize
  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)

!=======================================================================
CONTAINS

subroutine diff_fields(field1,field2,Dmax,typediff,ncomp,ntot)
    implicit none
    integer,intent(in)       :: ncomp
    real(mytype),dimension(xsize(1)*xsize(2)*xsize(3),ncomp),intent(in) :: field1,field2
    integer,intent(in)       :: typediff
    integer(kind=INT64),intent(in)       :: ntot
    real(mytype),intent(out) :: Dmax
    
    integer :: nv    
    integer :: ierror
    
    nv = xsize(1)*xsize(2)*xsize(3)
    Dmax = 0.
    
    if (typediff == 1) then
      Dmax = maxval(abs(field1 - field2))
      call MPI_Allreduce(MPI_IN_PLACE,Dmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
    elseif (typediff == 2) then
      Dmax = sum(abs(field1 - field2))
      call MPI_Allreduce(MPI_IN_PLACE,Dmax,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
      Dmax = Dmax / ntot
    else
      write(*,*) "ERROR : diff_fields : typediff different from 1 - not yet implemented "  
    end if

end subroutine diff_fields

end program test_field_divgrad


