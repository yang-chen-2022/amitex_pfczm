!123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!===================================================================================================
!
! MODULE IO2_AMITEX_MOD : 
!> Regroupement de subroutines AVANCEES liees aux entrees et aux sorties
!! Utilise les subroutine ELEMENTAIRES de io_amitex_mod 
!!
!!  
!!  Subroutines
!! - writeVTK :         Ecriture des fichiers VTK demandes en sortie
!!
!===================================================================================================
module io2_amitex_mod

  use ISO_FORTRAN_ENV

  use mpi
  use material_mod
  use loading_mod
#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi,       only : nrank
#else
  use decomp_2d, only : nrank,mytype 
#endif  
  use io_amitex_mod
  use error_mod
  use param_algo_mod
  use field_mod, only : grid
  use amitex_mod 
  use amitex_mpi_mod, only: amitex_mpi_comm
  implicit none

  private 
  
  
  !> Pointer "publiques" vers composantes de SIMU_AMITEX
  public :: times_io 

  !> Types publiques (pour definition de SIMU_AMITEX)
  public :: TIMES_IO2
  
  !> Fonctions publiques
  public :: writeVTK, print_field_vtk, read_field_vtk, print_coeff_varint_vtk

  interface read_field_vtk
    module procedure  read_field_vtk_real64
    module procedure  read_field_vtk_int32
    module procedure  read_field_vtk_int64
  end interface read_field_vtk

  interface print_field_vtk
    module procedure  print_fieldreal32_vtk
    module procedure  print_fieldreal64_vtk
  end interface print_field_vtk
  
!------------------------------------------------------------------------------
!> Structure TIMES_IO2

  type TIMES_IO2
     double precision       :: wvtk    = 0  ! temps cumule ecritures vtk
  end type TIMES_IO2

  type(TIMES_IO2),pointer   :: times_io
  
contains 

!===================================================================================================
!
!               SUBROUTINE writeVTK
!
!> Ecriture des fichiers VTK demandes en sortie
!! 
!! \param[in] ind_tps    (entier) pas de temps courant
!!
!! ATTENTION : POUR ASSURER LA COMPATIBILITE AVEC LE PROGRAMME DE POST-TRAITEMENT
!!             NE PAS MODIFIER LE NOM DES VARIABLES DANS L'EN-TETE DES VTK : SIG, DEF ET PK1
!!
!! ATTENTION : le nombre de quantites diffusives est ici limite a 9 
!!
!!
!===================================================================================================
  subroutine writeVTK(ind_tps)

    implicit none

    integer, intent(in) :: ind_tps
    integer             :: i,j,k
    integer(kind=8)     :: p
    character(len=200)  :: numS, suffixe, nomVarInt, numVar
    character(len=1)    :: tmp_i,tmp_j

    numS=""
    write(numS,"(I0)") ind_tps

    !============================================================================= MECANIQUE
    if(algo_param%Mechanics)then

    !-------------------------------------------Ecriture des contraintes
    if(extract%sigVTK)then
       ! contraintes
          do i=1,6
             write(tmp_i,"(I1)") i
             suffixe = "_sig"//trim(tmp_i)//"_"//trim(ADJUSTL(numS))//".vtk"
             call print_field_vtk(Sig(:,i),trim(fic_vtk)//trim(suffixe),"Sig_"//tmp_i)
          end do

       ! contraintes de Piola-Kirchhoff
          if(.not. algo_param%HPP) then
            do i=1,9
               write(tmp_i,"(I1)") i
               suffixe = "_pi"//trim(tmp_i)//"_"//trim(ADJUSTL(numS))//".vtk"
               call print_field_vtk(PK1(:,i),trim(fic_vtk)//trim(suffixe),"PK1_"//tmp_i)
            end do
          end if
    end if

    !-------------------------------------------Ecriture des deformations
    if(extract%defVTK)then
          !! Deformations
          !! En petites perturbations, cela correspond a 
          !! \epsilon = 1/2(\nabla u + {}^t\nabla u)
          !! En grandes transformations, cela correspond a \nabla u
          do i=1,algo_param%nTensDef
             write(tmp_i,"(I1)") i
             suffixe = "_def"//trim(tmp_i)//"_"//trim(ADJUSTL(numS))//".vtk"
             call print_field_vtk(Def(:,i),trim(fic_vtk)//trim(suffixe),"Def_"//tmp_i)
          end do
    end if

    !-------------------------------------------Ecriture des Variables Internes
    !suffixe du nom de fichier des variables internes
    suffixe = "_"//trim(ADJUSTL(numS))//".vtk"
    k=1
    do i=1,nmateriaux%n
       do j=1,size(extract%VarVTK(i)%val)
          ! si on doit extraire la variable interne
          if(extract%varVTK(i)%val(j)) then
             TEMPfield32=0.

             ! si le materiau 'i' se trouve dans le pinceau
             ! on stocke dans TEMPfield32 la valeur de la variable interne si elle existe (0 sinon)
             if (k<=size(MattotP)) then
             if(mattotP(k)%numM == i) then
                do p=1,size(MattotP(k)%pos,kind=8)
                   TEMPfield32(MattotP(k)%pos(p))=real(MattotP(k)%VarInt(j,p),kind=REAL32)
                end do
             end if
             end if

             ! Nom de fichier : nom_base_M'numM'_varInt'numVar'_'step'.vtk
             write(numVar,*) j   
             nomVarInt = "_varInt"
             nomVarInt = trim(nomVarInt)//trim(AdjustL(numVar))

             write(numVar,*) i
             nomVarInt = "M"//trim(AdjustL(numVar))//trim(nomVarint)

             ! Ecriture u fichier vtk (fonction generique) 
             call print_field_vtk(TEMPfield32,trim(fic_vtk)//"_"//trim(nomVarint)//trim(suffixe),trim(nomVarInt))
          end if
       end do
       if (k<=size(MattotP)) then
       if(MattotP(k)%numM==i)then
          k=k+1
       end if
       end if
    end do

    end if ! fin partie MECANIQUE

    !============================================================================= DIFFUSION
    if(algo_param%Diffusion)then

    !-------------------------------------------Ecriture des gradients de Q (Diffusion)
    if(extract%gradDVTK)then
       do j=1,algo_param%nVarD
         write(tmp_j,"(I1)") j
         if (j>9) call amitex_abort("Nombre de quantites diffusive > 9 (writeVTK)",2) !sinon coder j sur 2 caracteres et utiliser trim(adjustl())
            do i=1,3
               write(tmp_i,"(I1)") i
               suffixe = "_"//trim(ADJUSTL(numS))//"_gradQ"//trim(tmp_j)//"_"//trim(tmp_i)//".vtk"
               call print_field_vtk(GradQD(:,i,j),trim(fic_vtk)//trim(suffixe),"gradQ"//trim(tmp_j)//"_"//trim(tmp_i))
            end do
       end do
    end if

    !-------------------------------------------Ecriture des Flux de Q (Diffusion)
    if(extract%fluxDVTK)then
       do j=1,algo_param%nVarD
         write(tmp_j,"(I1)") j
         if (j>9) call amitex_abort("Nombre de quantites diffusive > 9 (writeVTK)",2) !sinon coder j sur 2 caracteres et utiliser trim(adjustl())
            do i=1,3
               write(tmp_i,"(I1)") i
               suffixe = "_"//trim(ADJUSTL(numS))//"_Flux"//trim(tmp_j)//"_"//trim(tmp_i)//".vtk"
               call print_field_vtk(FluxD(:,i,j),trim(fic_vtk)//trim(suffixe),"Flux"//trim(tmp_j)//"_"//trim(tmp_i))
            end do
       end do
    end if
    end if ! fin partie DIFFUSION

  end subroutine writeVTK
!==============================================================================

!===================================================================================================
!
!               SUBROUTINE print_field_vtk
!
!> Ecriture d'un champ Field de real(mytype) format en 32 bits 
!!
!! si ind_tps present : ecriture si extract%tpsVTK(ind_tps) = .true.
!!
!!
!!  \param[in]   Field real(mytype), tableau (nxP,nyP,nzP) (profil indifferent) 
!!  \param[in]   Filename, name of the vtk file
!!  \param[in]   Compname, name of the component writen in the vtk file
!!  \param[in]   (optional)ind_tps, index in the current loading step
!!  
!!
!===================================================================================================
  subroutine print_fieldreal64_vtk(Field,Filename,Compname,ind_tps)

    implicit none

    real(kind=REAL64),dimension(:),intent(in) :: Field
    character(len=*),intent(in)               :: Filename,Compname
    integer, intent(in),optional              :: ind_tps
    integer                                   :: ierror
    double precision                          :: t1
    logical                                   :: test
  
    call MPI_BARRIER(amitex_mpi_comm,ierror)          
    t1 = MPI_WTIME()
  
    test = .true.
    if (present(ind_tps)) then
       if(.not. extract%tpsVTK(ind_tps)) test = .false. 
    end if
    
    if (test) then
       if(nrank==0) call write_header_vtk(trim(Filename),&
                        grid%nx,grid%ny,grid%nz,&
                        grid%dx,grid%dy,grid%dz,&
                        grid%x0,grid%y0,grid%z0)
       call convert_field_for_write_bin32(TEMPfield32,Field) ! conversion bigendian
       call write_bin32(1,TEMPfield32,Compname,Filename)
    end if
    
    call MPI_BARRIER(amitex_mpi_comm,ierror)          
    times_io%wvtk = times_io%wvtk + MPI_WTIME() - t1

  end subroutine print_fieldreal64_vtk
!------------------------------------------------------  
  subroutine print_fieldreal32_vtk(Field,Filename,Compname,ind_tps)

    implicit none

    real(kind=REAL32),dimension(:),intent(in) :: Field
    character(len=*),intent(in)               :: Filename,Compname
    integer                                   :: ierror
    integer, intent(in),optional              :: ind_tps
    double precision                          :: t1
    logical                                   :: test
 
    call MPI_BARRIER(amitex_mpi_comm,ierror)          
    t1 = MPI_WTIME()
    
    test = .true.
    if (present(ind_tps)) then
       if(.not. extract%tpsVTK(ind_tps)) test = .false. 
    end if

    if (test) then
       if(nrank==0) call write_header_vtk(trim(Filename),&
                        grid%nx,grid%ny,grid%nz,&
                        grid%dx,grid%dy,grid%dz,&
                        grid%x0,grid%y0,grid%z0)
       call convert_field_for_write_bin32(TEMPfield32,Field) ! conversion bigendian
       call write_bin32(1,TEMPfield32,Compname,Filename)
    end if
    
    call MPI_BARRIER(amitex_mpi_comm,ierror)          
    times_io%wvtk = times_io%wvtk + MPI_WTIME() - t1

  end subroutine print_fieldreal32_vtk


!===================================================================================================
!
!               SUBROUTINE read_field_vtk
!
!> Lecture d'un champ Field de reels (une seule composante, decomposition en pinceaux X) 
!>                                   
!!
!!
!!  \param[out]  Field real(mytype), tableau (nxP x nyP x nzP) (profil 1D implicite) 
!!  \param[in]   Filename, name of the vtk file
!!  
!!
!===================================================================================================
!------------------------------------------------------
  subroutine read_field_vtk_real64(Field,Filename)

    implicit none

    real(mytype), dimension(:),intent(out)  :: Field
    character(len=*),intent(in)             :: Filename
    integer                                 :: nx,ny,nz, type_mpi
    real(mytype)                            :: dx,dy,dz,x0,y0,z0
    integer,parameter                       :: ipencil=1 ! decompopsition suivant X

    ! read_header_vtk 
    call read_header_vtk(Filename, nx,ny,nz,dx,dy,dz,x0,y0,z0, type_mpi)

    ! check the consistency between nx,ny,nz and the grid used in the code
    call check_grid(nx,ny,nz)

    ! read binary data
    call read_bin(Filename,ipencil,type_mpi,Field)

  end subroutine read_field_vtk_real64
!------------------------------------------------------
  subroutine read_field_vtk_int32(Field,Filename)

    implicit none

    integer(kind=INT32), dimension(:),intent(out)  :: Field
    character(len=*),intent(in)             :: Filename
    integer                                 :: nx,ny,nz, type_mpi
    real(mytype)                            :: dx,dy,dz,x0,y0,z0
    integer,parameter                       :: ipencil=1 ! decompopsition suivant X

    ! read_header_vtk 
    call read_header_vtk(trim(Filename), nx,ny,nz,dx,dy,dz,x0,y0,z0, type_mpi)

    ! check the consistency between nx,ny,nz and the grid used in the code
    call check_grid(nx,ny,nz)

    ! read binary data
    call read_bin(trim(Filename),ipencil,type_mpi,Field)

  end subroutine read_field_vtk_int32
!------------------------------------------------------
  subroutine read_field_vtk_int64(Field,Filename)

    implicit none

    integer(kind=INT64), dimension(:),intent(out)  :: Field
    character(len=*),intent(in)             :: Filename
    integer                                 :: nx,ny,nz, type_mpi
    real(mytype)                            :: dx,dy,dz,x0,y0,z0
    integer,parameter                       :: ipencil=1 ! decompopsition suivant X

    ! read_header_vtk 
    call read_header_vtk(trim(Filename),nx,ny,nz,dx,dy,dz,x0,y0,z0, type_mpi)

    ! check the consistency between nx,ny,nz and the grid used in the code
    call check_grid(nx,ny,nz)

    ! read binary data
    call read_bin(trim(Filename),ipencil,type_mpi,Field)

  end subroutine read_field_vtk_int64


  !------------------------------------------------------
  subroutine check_grid(nx,ny,nz)

    implicit none 

    integer,intent(in) :: nx,ny,nz

    if (nx .ne. grid%nx) call amitex_abort("Error in read_field_vtk (nx in vtk different from grid%nx) " ,2)
    if (ny .ne. grid%ny) call amitex_abort("Error in read_field_vtk (ny in vtk different from grid%ny) " ,2)
    if (nz .ne. grid%nz) call amitex_abort("Error in read_field_vtk (nz in vtk different from grid%nz) " ,2)
  
  end subroutine check_grid

!===================================================================================================
!
!               SUBROUTINE print_coeff_varint_vtk
!
!> vtk output for coefficients or internal variables in MattotP
!>                                   
!!
!!
!!  \param[in]  coeff2print, index of the coefficient to output 
!!  \param[in]  varint2print, index of the internal variable to output 
!!  
!!              No output if the coeff2print or varint2print = 0 (or <0)
!!
!===================================================================================================
  subroutine print_coeff_varint_vtk(coeff2print,varint2print)
  
    implicit none 

    integer,intent(in) :: coeff2print, varint2print
    
    integer            :: ierror
    logical            :: bool0,bool1
    integer            :: i,j,k 
    character(len=200) :: tmp_char1,tmp_char
    real(mytype)       :: t1

  if ((coeff2print > 0) .OR. (varint2print > 0)) then
  
  do i=1,nmateriaux%n
     ! recontruction du champ de coeff pour le materiau i (range dans sig(:,1))
     !     uniquement si le coefficient existe dans ce materiau
     ! idem pour variables internes dans sig(:,2)
     bool0 = .false.
     bool1 = .false.
     sig(:,1) = 0.
     sig(:,2) = 0.

     !--- CREATION DU CHAMP DE COEFF/VAR. INT. A SORTIR    
     do j=1,size(mattotP)
     if ((mattotP(j)%numM==i) .and. (.not. mattotP(j)%Interphase)) then

       !--- Coefficients
       if(coeff2print .le. size(mattotP(j)%Coeff,1).AND. (coeff2print > 0)) then
         sig(mattotP(j)%pos(1:mattotP(j)%Zone(1,1)),1) = mattotP(j)%Coeff(coeff2print,1)
         do k=2,size(mattotP(j)%Zone,1)
            if (allocated(mattotP(j)%Zones_interphase)) then
               if (any(mattotP(j)%Zones_interphase == mattotP(j)%zone(k,2))) then
                  ! si zone interphase : cycle, pas de calculs
                  cycle
               end if
            end if
            sig(mattotP(j)%pos(mattotP(j)%Zone(k-1,1)+1:mattotP(j)%Zone(k,1)),1)=mattotP(j)%Coeff(coeff2print,k)
         end do
         bool0=.true.
       end if
       
       !--- Variables internes
       if(varint2print .le. mattotP(j)%Nvarint .AND. (varint2print > 0)) then
         sig(mattotP(j)%pos,2) = mattotP(j)%Varint0(varint2print,:)
         bool1=.true.
       end if
  
     end if
     end do
 
     CALL MPI_ALLREDUCE(MPI_IN_PLACE, bool0, 1, MPI_LOGICAL, MPI_LOR, amitex_mpi_comm, ierror)
     CALL MPI_ALLREDUCE(MPI_IN_PLACE, bool1, 1, MPI_LOGICAL, MPI_LOR, amitex_mpi_comm, ierror)

     !--- SORTIE VTK DU CHAMP DE COEFF/VAR. INT.    
     !--- Coefficients
     if (bool0 .AND. (coeff2print > 0)) then
       t1 = MPI_WTIME()
       !
       write(tmp_char1,"(I4)") coeff2print
       tmp_char="_coeff_"//trim(adjustl(tmp_char1))
       write(tmp_char1,"(I4)") i
       tmp_char=trim(tmp_char)//"_mat_"//trim(adjustl(tmp_char1))
       call print_field_vtk(Sig(:,1),trim(fic_vtk)//trim(tmp_char)//".vtk",trim(tmp_char))
       !
       call MPI_Barrier(amitex_mpi_comm, ierror)
       if(nrank==0) then
       write(Flog,fmt="(A,E15.8)") "VTK writing time on proc 0 (s) (option -coeff2print) : ",MPI_WTIME()-t1
       write(OUTPUT_UNIT,fmt="(A,E15.8)") "VTK writing time on proc 0 (s) (option -coeff2print) : ",MPI_WTIME()-t1
       end if
     end if

     !--- Variables internes
     if (bool1 .AND. (varint2print > 0)) then
       t1 = MPI_WTIME()
       !
       write(tmp_char1,"(I4)") varint2print
       tmp_char="_varint_"//trim(adjustl(tmp_char1))
       write(tmp_char1,"(I4)") i
       tmp_char=trim(tmp_char)//"_mat_"//trim(adjustl(tmp_char1))
       call print_field_vtk(Sig(:,2),trim(fic_vtk)//trim(tmp_char)//".vtk",trim(tmp_char))
       !
       call MPI_Barrier(amitex_mpi_comm, ierror)
       if(nrank==0) then
       write(Flog,fmt="(A,E15.8)") "VTK writing time on proc 0 (s) (option -varint2print) : ",MPI_WTIME()-t1
       write(OUTPUT_UNIT,fmt="(A,E15.8)") "VTK writing time on proc 0 (s) (option -varint2print) : ",MPI_WTIME()-t1
       end if
     end if

  end do

  end if

  
  end subroutine print_coeff_varint_vtk

end module io2_amitex_mod


