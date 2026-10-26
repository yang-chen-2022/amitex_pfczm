!===============================================================================
!
! MODULE DEFORMED_SHAPE : 
!> Regroupement de subroutines pour le programme deformedShape
!! 
!!
!!  Subroutines
!! - lire_commande_post :    Lit la ligne de commande
!! - computeDep :            Calcule le champ de deplacement
!! - read_header_vtk_post :  Lecture de l'en-tete du fichier vtk principal
!! - read_bin_post :      Lecture de la partie binaire d'un fichier vtk
!! - print_vtk :          Ecriture du fichier VTK
!! - write_bin_dep :     Ecriture du champ de deplacement en binaire
!! - write_header_vtk_dep :  Ecriture de l'en-tete du fichier VTK (configuration actuelle)
!!
!!
!!
!===============================================================================
module deformed_shape_mod
  use ISO_FORTRAN_ENV
  use mpi
#ifdef NEW_2DECOMP
  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
#else
  use decomp_2d
#endif             
  use decomp_2d_io
  use decomp_2d_fft
  use error_mod
  use io_amitex_mod
  use amitex_mod , only : TEMPfield32
  use green_mod
  use field_mod


  private

  ! liste des subroutines publiques  
  public :: lire_commande_post,computeDep, printVTK


contains


!==============================================================================
!
!>   Subroutine de lecture des arguments de la ligne de commande.
!!
!! Lit la ligne de commande qui doit etre de la forme\n
!! mpirun ./deformedShape loading_step <pas_de_calcul> <fic_vtk_d_entree> [<fic_vtk_de_sortie>] -sig (ou -def ou -pk1)
!!
!! - <pas_de_calcul>      indice du pas de calcul
!! - <fic_vtk_d_entree>   racine des fichiers vtk sortis par le programme amitex_fftp
!! - <fic_vtk_de_sortie>  racine des fichiers sur lesquels seront ecrits les resultats.
!!                        Si cette valeur n'est pas precise on ecrit dans les fichiers "sortie"
!!
!! 
!! \param[out] numS: (chaine de caracteres) loading step
!! \param[out] fic_in: (chaine de caracteres) racine des fichiers d'entree
!! \param[out] fic_out: (chaine de caracteres) racine des fichiers de sortie
!! \param[out] choice: (logical) character(4) -sig or -def or -pk1
!! \param[out] outprog: (logical) .true. to leave the program
!!
!! WARNING : since version 6.2.0, the name of output vtk file is 
!!           "root_name_in"_"field"_"step".vtk
!! WARNING : since version 8.14.1, amitex systematically outputs 
!!           one vtk file / component
!!            
!==============================================================================
subroutine lire_commande_post(numS, fic_in, fic_out,choice,outprog)


  implicit none

  character(len=*),intent(inout)  :: fic_in, fic_out
  character(len=20),intent(out)   :: numS
  character(len=16),intent(out)   :: choice 
  logical,intent(out)             :: outprog
  logical                         :: existe
  integer                         :: i,nargs
  character(len=200)              :: fic_in_i
  character(len=200)              :: err_msg,char_tmp
  character(len=*), parameter     :: help_msg = achar(10)//&
         "Help deformedShape : mpirun deformedShape loading_step <root_filename_in> &
                   &[<root_filename_out>] -def (or -sig or -pk1 -or -M(numM)_varInt(numVI))"//achar(10)//&
         "--------------------        default :<root_filename_out> = <root_filename_in>"//achar(10)//&
         achar(10)//&
         "      Output the Displacement gradient (-def), or Cauchy stress (-sig) or 1st Piola-Kirchoff stress (-pk1), field"&
         //achar(10)//&
         "      in the deformed configuration from their quantities on the initial configuration (Amitex_FFTP output)."&
         //achar(10)//&
         "      For internal variable, the example -M2_varInt4 is for the Internal Variable 4 in Material 2."& 
         //achar(10)//&
         "      The user must provide : the root_filename_in and the loading_step (see AMITEX_FFT output)."&
         //achar(10)//&
         "      If not given, the root_filename_out will be the root_filename_in"

  !! Initialize
  outprog = .false.
  fic_out = ""

  nargs = command_argument_count()
  !! Help info (if command line with no entry or -help)
  if(nargs == 1 .OR. nargs == 0) then
  if (nargs == 1) then
     call getarg(1,char_tmp)
  else
     char_tmp="-help"
  end if
  if (char_tmp=="-help") then
     call write_stdout0(help_msg)
     outprog=.true.
     goto 666
  end if
  end if

  !! On verifie qu'il y a assez d'arguments
  if(nargs < 3) then
     call amitex_abort(help_msg,2,0)
  end if

  !! On recupere le pas de chargement des vtks que l'on souhaite traiter
  call getarg(1,numS)
  if (nrank==0) call write_stdout0("     The loading_step chosen is : "// numS)

  !! On recupere la racine des fichiers de deformations
  call getarg(2,fic_in)
  if (nrank==0) call write_stdout0("     The root_filename_in is    : "// fic_in)

  INQUIRE(FILE=trim(fic_in)//"_def_"//trim(ADJUSTL(numS))//".vtk", EXIST=existe)
  !! Cette option est desormais supprimee => on utilise un champ / composante
  if(existe) then 
     call amitex_abort(achar(10)//"SINCE V8.14.1 : "//achar(10)//&
                  "Using a single vtk file for GRAD_U is no more possible."//achar(10)//&
                  "      => USE ONE VTK FILE PER COMPONENT OF GRAD_U",2,0)
  else
  !! Si le fichier def.vtk n'existe pas, on verifie que les donnees ne sont pas en
  !! plusieurs fichiers
     do i=1,9
        write (fic_in_i, "(A,I0,A)") trim(fic_in)//"_def",i,"_"//trim(ADJUSTL(numS))//".vtk"
        INQUIRE(FILE=fic_in_i, EXIST=existe)
        !! Si un des fichiers de deformation defi.vtk n'existe pas 
        !! on affiche une erreur.
        if((.NOT.existe)) then 
           write (err_msg,"(A,A,A)")"Le fichier VTK de deformation ",&
                trim(fic_in_i)," n'existe pas"
           call amitex_abort(err_msg,2,0)
        end if
     end do
  end if

  !! lecture et test du dernier argument (-def or -sig or -pk1 or -M'numM'_varInt'numVint')
  call getarg(nargs,choice)
  if (nrank==0) call write_stdout0("     The option is              : "// trim(choice))

  if (choice(1:4) /= "-def" .AND. choice(1:4) /= "-sig" .AND. choice(1:4) /= "-pk1" .AND. &
      choice(1:2) /= "-M") then
     call amitex_abort(help_msg,2,0)
  end if 
  if (choice(1:2) == "-M" .AND. index(choice,"varInt") ==0 ) then
       call amitex_abort(help_msg,2,0)
  end if 

  !! S'il y a un 2e argument, celui-ci sera la racine du fichier vtk de sortie
  if(nargs == 4 ) then
     call getarg(3,fic_out)
  else
     !! Sinon, par defaut, la racine s'appelle "sortie"
     fic_out = fic_in
  end if
  if (nrank==0) call write_stdout0("     The root_filename_out is   : "// fic_in)


666  fic_out = trim(fic_out)

end subroutine lire_commande_post
!-------------------------------------------------------------------



!==================================================================================
!                           Subroutine computeDep
!> Lit le fichier contenant le champ de deformation et calcule le champ de deplacement
!!
!!
!! \param[in]  fic_in: (chaine de caracteres) fichier vtk ou sont ecrits les deformations
!! \param[inout]  Dep: champ de deplacement
!! \param[inout]  DepF: champ de deplacement dans l'espace de Fourier
!! \param[inout]  GradU: gradient du deplacement
!! \param[inout]  GradUF: gradient du deplacement dans l'espace de Fourier
!! \param[in] fft_start, fft_end:    limites des tableaux ou seront stockees les FFT   
!! \param[in] ntot:    Nombre de voxels global
!! \param[in] dx,dy,dz:    Spacing
!!
!==================================================================================
subroutine computeDep(numS, fic_in, dep, gradu,depf,graduf,nx,ny,nz,dx,dy,dz,ntot)

  implicit none

  character(len=20),intent(in)                 :: numS
  real(mytype),allocatable,dimension(:,:,:,:)  :: Freq_0
  integer, intent(in)                          :: nx,ny,nz
  real(mytype),dimension(3,xsize(1)*xsize(2)*xsize(3))                            :: Dep
  real(mytype),dimension(3,xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: Dep1
  real(mytype),dimension(:,:),intent(inout)    :: GradU
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:9),intent(inout) :: GradUF
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3),intent(inout) :: DepF
  real(mytype), intent(in)        :: dx,dy,dz
  real(mytype),dimension(9)       :: F,F0
  character(len=*), intent(in)    :: fic_in
  character(len=200)              :: fic_def_i
  integer                         :: ierr
  integer(kind=4)                 :: i,j,h
  complex(mytype)                 :: imp = cmplx(0,1,kind=mytype) ! le nombre complexe i
  complex(mytype)                 :: decal
  real(mytype)                    :: norme2
  integer(kind=8), intent(in)     :: ntot

  !! on recupere le champ de deformation
  do i=1,9
     write (fic_def_i, "(A,I0,A)") trim(fic_in)//"_def",i,"_"//trim(ADJUSTL(numS))//".vtk"
     call read_bin(fic_def_i,1,mpi_real,GradU(:,i))
  end do

  !! on calcule la transformee de fourier de la deformation
  call field_fft(GradU,GradUF,9,1)
  GradUF = GradUF/real(ntot,mytype)
  F0 = 0
  if(fft_start(1) == 1 .AND. fft_start(2) == 1 .AND. fft_start(3) == 1) then
     do i = 1,3
        ! on rajoute 1 sur la diagonale pour calculer x = id + moy(grad(u)) 
        F0(i)= 1. + GradUF(1,1,1,i)
     end do
     do i = 4,9
        F0(i)= GradUF(1,1,1,i)
     end do
  end if
  !! Astuce pour distribuer F sur tous les procs sans connaitre 
  !! le proc. contenant la fréquence (1,1,1) :
  !! On initialise a 0 au depart et on fait la somme a la fin du calcul
  call MPI_AllReduce(F0,F,9,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
  
  !pour tester utilisation des frequences filtrees : pour l'instant sans succes
  !allocate(Freq_0(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3))
  !call initFreq(Freq_0,"no_filter",0._mytype)
  !call initFreq(Freq,"hexa",1._mytype)

  !! on calcule la transformee de fourier de u tilde par la formule 
  !! \f[ \hat{u} = - i\frac{graduf*freq}{\abs{freq}^2} \f]
  do i=fft_start(3),fft_end(3)
     do j=fft_start(2),fft_end(2)
        do h=fft_start(1),fft_end(1)
           norme2 = freq(h,j,i,1) **2 +  freq(h,j,i,2) **2 +  freq(h,j,i,3) **2
	   ! Pour tester utilisation des freq filtree
           ! decal = exp(imp * (Freq_0(h,j,i,1) * dx/2._mytype + Freq_0(h,j,i,2) * dy/2._mytype +Freq_0(h,j,i,3) * dz/2._mytype))
	   ! 
           if ((h .eq. 1) .and. (i .eq. 1) .and. (j .eq. 1)) then
              norme2= 1._mytype
              decal=1._mytype
           end if
           !! notation : ordre 11 22 33 12 13 23 21 31 32
           !!                   1  2  3  4  5  6  7  8  9
           decal = 1._mytype
           DepF(h,j,i,1) =  - imp*( GradUF(h,j,i,1) * Freq(h,j,i,1) + &
                GradUF(h,j,i,4) * Freq(h,j,i,2) +  GradUF(h,j,i,5) * Freq(h,j,i,3))*decal/norme2
           DepF(h,j,i,2) =  - imp*( GradUF(h,j,i,7) * Freq(h,j,i,1) + &
                GradUF(h,j,i,2) * Freq(h,j,i,2) +  GradUF(h,j,i,6) * Freq(h,j,i,3))*decal/norme2
           DepF(h,j,i,3) =  - imp*( GradUF(h,j,i,8) * Freq(h,j,i,1) + &
                GradUF(h,j,i,9) * Freq(h,j,i,2) +  GradUF(h,j,i,3) * Freq(h,j,i,3))*decal/norme2
        end do
     end do
  end do

  !! On effectue la transformee de Fourier inverse pour retrouver \tilde{u}
  do i = 1,3
     call field_ifft(Dep1(i,:,:,:),depf(:,:,:,i),1,1)
  end do

  do i=xstart(3),xend(3)
     do j=xstart(2),xend(2)
        do h=xstart(1),xend(1)
           ! index = i+(j-1)*xsize(3)+(h-1)*xsize(3)*xsize(2)
           !! u = \tilde{u} + F * X
           dep1(1,h,j,i) = dep1(1,h,j,i) + (-dx/2._mytype + dx * h)*F(1)&
           + (-dy/2._mytype + dy * j)*F(6) + (-dz/2._mytype + dz * i)*F(5)
           dep1(2,h,j,i) = dep1(2,h,j,i) + (-dx/2._mytype + dx * h)*F(9)&
           + (-dy/2._mytype + dy * j)*F(2) + (-dz/2._mytype + dz * i)*F(4)
           dep1(3,h,j,i) = dep1(3,h,j,i) + (-dx/2._mytype + dx * h)*F(8)&
           + (-dy/2._mytype + dy * j)*F(7) + (-dz/2._mytype + dz * i)*F(3)
           ! dep1(1,h,j,i) = dep1(1,h,j,i) + dx * (h-1)*F(1)&
           ! + dy * (j-1)*F(6) + dz * (i-1)*F(5)
           ! dep1(2,h,j,i) = dep1(2,h,j,i) + dx * (h-1)*F(9)&
           ! + dy * (j-1)*F(2) + dz * (i-1)*F(4)
           ! dep1(3,h,j,i) = dep1(3,h,j,i) + dx * (h-1)*F(8)&
           ! + dy * (j-1)*F(7) + dz * (i-1)*F(3)
        end do
     end do
  end do
  dep = reshape(dep1,(/3,xsize(1)*xsize(2)*xsize(3)/))
  if(allocated(Freq_0)) deallocate(Freq_0)
     
end subroutine computedep
!===============================================================================

!===============================================================================
!
!            subroutine read_bin_post
!
!> lecture de la partie binaire d'un fichier vtk (sortie d'amitex_fftp)
!!
!!  ATTENTION :
!!  Suppose que la taille en octets de deux sauts de ligne et des lignes
!!       SCALARS XXX_i float
!!       LOOKUP_TABLE default
!!  => cela represente 43 octets 
!!
!! \param[in] ipencil: orientation des pinceaux du tableau var
!! \param[inout] var: tableau a lire 
!! \param[in] nomvar: nom de la variable a extraire du vtk
!! \param[in] nomfic: nom du fichier
!! \param[in] type_mpi: type des donnees
!! \param[in] index: composante a lire
!! \param[in] size: nombre de composantes pour cette variable
!!
!!
!-------------------------------------------------------------------------------
  subroutine read_bin_post(nomfic,ipencil,type_mpi,var,index,size)
    
    use mpi

    implicit none
    
    character(len=*), intent(in)                     :: nomfic
    integer, intent(in)                              :: ipencil,type_mpi
    integer, intent(in), optional                    :: index,size
    real(kind=REAL32), dimension(:), intent(inout)   :: var

    type(decomp_info)                                :: decomp
    integer(kind=mpi_offset_kind)                    :: disp, filesize
    integer, dimension(3)                            :: sizes, subsizes, starts
    integer                                          :: ierror, newtype, fh, noctets,i
    integer(kind=1)                                  :: debut

    !variable du traitement de l'endianness
    integer(kind=1), dimension(2)::testendian
    integer (kind=2)             ::testendian2
    equivalence(testendian,testendian2)

    !variable pour definir un type de donnees adapte
    integer :: type_end
    integer,parameter ::nbbl=4
    integer,dimension(:),allocatable ::blocklen
    integer(kind=mpi_address_kind),dimension(:), allocatable :: dplt

    !! nOctets = nombre d'octets par nombre en fonction du type de nombre
    call MPI_TYPE_SIZE(type_mpi,nOctets,ierror)

    testEndian=(/int(1,kind=1),int(0,kind=1)/)

    call get_decomp_info(decomp)

    !! tailles, debut et etendues des sous-tableaux
    !! starts ramene a 0 (pour MPI) 
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)

    !! taille et bornes des tableaux selon l'orientation:
    !! pinceaux selon x, y ou z   
    if (ipencil == 1) then
       subsizes(1) = decomp%xsz(1)
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = decomp%xst(1)-1
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = decomp%ysz(2)
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = decomp%yst(2)-1
       starts(3) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = decomp%zsz(3)
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = decomp%zst(3)-1
    endif

    !! Creation du type utilise par mpi en fonction de l'endianness
    if(testEndian2==1)then

       allocate(blockLen(nOctets),dplt(nOctets))
       blockLen=1
       do i=1,nOctets
          dplt(i)=nOctets-i
       end do
       if(nOctets>1)then
          !creation de la structure d'entier 
          call mpi_type_create_hindexed(int(nOctets,kind=4),blockLen,dplt,MPI_integer1,type_end,ierror)
          call MPI_Type_commit(type_end,ierror)
       else 
          type_end = type_mpi
       end if
    else
       type_end = type_mpi
    end if


    !creation de la structure de donnee (tableau d'entier defini ci-dessus)
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, type_mpi, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)

    call MPI_FILE_OPEN(MPI_COMM_WORLD, nomFic, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, &
         fh, ierror)

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    !recherche du debut de l'enregistrement des donnees
    !estimation du debut de l'enregistrement
    call MPI_File_get_size(fh,filesize,ierror)

    disp = sizes(1)
    disp = disp*sizes(2)
    disp = disp*sizes(3)
    if(Present(index)) then
       if(.not. present(size)) then
          call amitex_abort( "Erreur si on precise un indice on doit preciser le nombre &
               &de composantes. (read_bin_post)" ,2)
       end if
       disp = disp*int(nOctets,kind=8)*(size+1-index)
       !! On rajoute la taille en octets de deux sauts de ligne et des lignes
       !!       SCALARS XXX_i float
       !!       LOOKUP_TABLE default
       !! Cela represente 43 octets -> 43 dans la ligne ci-dessous
       disp = disp + 43*(size-index)
    else
       disp = disp*int(nOctets,kind=8)
    end if
    disp = (filesize) - disp
    if(disp <1)then  
       call amitex_abort( "Erreur lors de la lecture du fichier vtk " ,2)
    end if

    call MPI_File_seek(fh,disp,MPI_SEEK_SET ,ierror)
    call MPI_File_read(fh, debut, 1, MPI_INTEGER1, MPI_STATUS_IGNORE, ierror)

    !! Recherche du premier saut de ligne avant le debut estime des donnees
    !! (10 en ASCII correspond a un saut de ligne)
    i=0
    debut=1
    do while( (debut/=10) )
       i=i+int(1,kind=1)
       disp = disp -1
       call MPI_File_seek(fh,disp,MPI_SEEK_SET ,ierror)
       call MPI_File_read(fh, debut, 1, MPI_INTEGER1, MPI_STATUS_IGNORE, ierror)
    end do
    !On place le debut apres le saut de ligne trouve (ci-dessus)
    disp = disp +1

    call MPI_FILE_SET_VIEW(fh,disp,type_end, &
         newtype,'native',MPI_INFO_NULL,ierror)

    !lecture des donnees
    call MPI_FILE_READ_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         type_end, MPI_STATUS_IGNORE, ierror)
    !fermeture du fichier
    call MPI_FILE_CLOSE(fh,ierror)
    !destruction des type et tableau
    call MPI_TYPE_FREE(newtype,ierror)
    if(testEndian2==1 .AND. nOctets/=1) call MPI_TYPE_FREE(type_end,ierror)
    if(allocated(blockLen)) deallocate(blockLen)
    if(allocated(dplt)) deallocate(dplt)

    return

  end subroutine  read_bin_post
!===============================================================================

!===============================================================================
!
!            SUBROUTINE PRINT_VTK
!> Ecriture du fichier VTK contenant les differents champs dans le maillage deforme
!!
!!
!!   
!!
!! \param[in] fic_in: racine des fichiers de donnees
!! \param[in] fic_out: nom du fichier de sortie
!! \param[in] Fvtk: unite du fichier de sortie
!! \param[in] Dep: Deplacement dans l'espace reel
!! \param[in] field: vecteur contenant au depart le gradient de u
!! \param[in] nx,ny,nz:      (entiers) dimensions de la cellule
!!
!!
!-------------------------------------------------------------------------------
subroutine printVTK(numS, fic_in,fic_out,Fvtk,Dep,field,nx,ny,nz,choice)

  implicit none

  character(len=20),intent(in)      :: numS
  integer                           :: newtype, fh, type_reel
  integer,dimension(4)              ::blockLen
  integer(kind=mpi_offset_kind)     :: disp, filesize
  character(len=*),intent(in)       :: fic_in, fic_out
  character(len=*),intent(in)       :: choice
  character(len=40)                 :: nomVar
  character(len=200)                :: fic_in_i,fic_out_i
  real(kind=REAL32),dimension(:,:),intent(in)                :: Dep
  real(kind=REAL32),dimension(xsize(1)*xsize(2)*xsize(3)*3)  :: Dep1
  real(mytype),dimension(xsize(1)*xsize(2)*xsize(3))    :: temp
  real(mytype),dimension(:,:),intent(inout)             :: field
  integer, intent(in)               :: Fvtk, nx,ny,nz
  integer(kind=8)                   :: ntot,index
  integer                           :: i,j,h,ios,ierror,rank
  integer(kind=MPI_ADDRESS_KIND),dimension(4):: dplt
  logical                           :: existe


  !! Un fichier par type de tenseur (Def, Sig ou Pi)
  if     (choice(1:4) == "-def") then
    fic_out_i = trim(fic_out)//"_Def_"//trim(ADJUSTL(numS))//".vtk"
  elseif (choice(1:4) == "-sig") then
    fic_out_i = trim(fic_out)//"_Sig_"//trim(ADJUSTL(numS))//".vtk"
  elseif (choice(1:4) == "-pk1") then
    fic_out_i = trim(fic_out)//"_Pi_"//trim(ADJUSTL(numS))//".vtk"
  elseif (choice(1:2) == "-M") then
    fic_out_i = trim(fic_out)//"_V"//trim(choice(2:))//"_"//trim(ADJUSTL(numS))//".vtk"
  end if

  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
  ntot = nx
  ntot = ntot*ny
  ntot = ntot*nz
  !! Ecriture de l'en-tete
  call write_header_vtk_dep(fic_out_i, Fvtk,nx,ny,nz)

  !! Ecriture en binaire de la liste des points
  call write_bin_dep(1,Dep,fic_out_i)
  !! Ecriture du champ de deformation
  call MPI_Barrier(MPI_COMM_WORLD,ierror) 
  if(nrank==0) then
     open(unit=Fvtk, file= fic_out_i,form="formatted", &
          status="old", action="write", position="APPEND",iostat= ios)
     if ( ios /= 0 ) then
        print *," Probleme a l'ouverture (fichier: ",fic_out_i,")",ios," (printVTK)"
        stop 
     end if

     write(Fvtk,"(A,I0)")"POINT_DATA   ", ntot
     close(Fvtk)
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierror) 

  !! Gradient de u
  if (choice(1:4)=="-def") then
  do i=1,9
     write (nomVar, "(A,I0)") "GradU_",i
     !TEMPfield32=real(field(:,i))
     !call write_bin(1,TEMPfield32,trim(nomVar),fic_out_i,Fvtk)
     call convert_field_for_write_bin32(TEMPfield32,Field(:,i)) ! conversion bigendian
     call write_bin32(1,TEMPfield32,trim(nomVar),fic_out_i)
  end do
  end if

  !! Contrainte de Cauchy
  if (choice(1:4)=="-sig") then
  do i=1,6
     write (fic_in_i, "(A,I0,A)") trim(fic_in)//"_sig",i,"_"//trim(ADJUSTL(numS))//".vtk"
     INQUIRE(FILE=fic_in_i, EXIST=existe)
     !! permet d'ignorer les donnees de contrainte absentes
     if(.NOT. existe) cycle
     write (nomVar, "(A,I0)") "Sig_",i
     !! on recupere le champ de contrainte
     !call read_bin_post(fic_in_i,1,mpi_real,temp)
     call read_bin(fic_in_i,1,MPI_REAL,temp)
     !TEMPfield32=real(temp)
     !call write_bin(1,TEMPfield32,trim(nomVar),fic_out_i,Fvtk)
     call convert_field_for_write_bin32(TEMPfield32,temp) ! conversion bigendian
     call write_bin32(1,TEMPfield32,trim(nomVar),fic_out_i)
  end do
  end if

  !! Contrainte de Piola-Kirchhoff
  if (choice(1:4)=="-pk1") then
  do i=1,9
     write (fic_in_i, "(A,I0,A)") trim(fic_in)//"_pi",i,"_"//trim(ADJUSTL(numS))//".vtk"
     INQUIRE(FILE=fic_in_i, EXIST=existe)
     !! permet d'ignorer les donnees de contrainte absentes
     if(.NOT. existe) cycle
     write (nomVar, "(A,I0)") "PK1_",i
     !! on recupere le champ de contrainte
     !call read_bin_post(fic_in_i,1,mpi_real,temp)
     call read_bin(fic_in_i,1,MPI_REAL,temp)
     !TEMPfield32=real(temp)
     !call write_bin(1,TEMPfield32,trim(nomVar),fic_out_i,Fvtk)
     call convert_field_for_write_bin32(TEMPfield32,temp) ! conversion bigendian
     call write_bin32(1,TEMPfield32,trim(nomVar),fic_out_i)
  end do
  end if

  !! Variable interne
  if (choice(1:2)=="-M") then
     write (fic_in_i, "(A)") trim(fic_in)//"_"//trim(choice(2:))//"_"//trim(ADJUSTL(numS))//".vtk"     
     write (nomVar, "(A)") trim(choice(2:))

     INQUIRE(FILE=fic_in_i, EXIST=existe)
     if(.NOT. existe) then
        call amitex_abort("Internal variable file not found"//achar(10),2,0)
     end if

     !! on recupere le champ de donnees binaire
     call read_bin(fic_in_i,1,MPI_REAL,temp)
     call convert_field_for_write_bin32(TEMPfield32,temp) ! conversion bigendian
     call write_bin32(1,TEMPfield32,trim(nomVar),fic_out_i)
  end if
 
end subroutine printVTK
!===============================================================================

!===============================================================================
!
!            SUBROUTINE WRITE_BIN_DEP
!> Ecriture du champ de deplacement en binaire
!!
!!
!!   
!!
!! \param[in] ipencil: orientation des pinceaux du tableau var
!! \param[in] var: tableau a ecrire 
!! \param[in] nomFic: nom du fichier
!!
!!
!-------------------------------------------------------------------------------
  subroutine write_bin_dep(ipencil,var0,nomFic)
    
    implicit none
    
    integer, intent(in)                 :: ipencil
    real(kind=REAL32),dimension(1:3,xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in)::var0
    real(kind=REAL32),dimension(1:3,xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) ::var
    character(len=*),intent(in)         :: nomFic
    
    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(4) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh, ios, type_reel,i
    
    !variable du traitement de l'endianness
    integer(kind=1), dimension(2)::testEndian
    integer (kind=2)             ::testEndian2
    equivalence(testEndian,testEndian2)
    
    !variable pour definir un type de donnees adapte
    !On ecrit des reels en simple precision (4 octets)
    integer,dimension(4) ::blockLen
    integer(kind=MPI_ADDRESS_KIND),dimension(4):: dplt

    testEndian=(/int(1,kind=1),int(0,kind=1)/)
    var = var0  
    if(testEndian2==1)call SWAP_F4(var)
    
    call MPI_Barrier(MPI_COMM_WORLD,ierror) 
    call get_decomp_info(decomp)
    
!tailles, debut et etendues des sous tableau 
!starts ramene a 0 (pour MPI) 
    sizes(1) = 3
    sizes(2) = decomp%xsz(1)
    sizes(3) = decomp%ysz(2)
    sizes(4) = decomp%zsz(3)
    
!taille et bornes des tableau selon l'orientation:
!pinceaux selon x, y ou z   
    if (ipencil == 1) then
       subsizes(1) = 3
       subsizes(2) = decomp%xsz(1)
       subsizes(3) = decomp%xsz(2)
       subsizes(4) = decomp%xsz(3)
       starts(1) = 0
       starts(2) = decomp%xst(1)-1
       starts(3) = decomp%xst(2)-1
       starts(4) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = 3
       subsizes(2) = decomp%ysz(1)
       subsizes(3) = decomp%ysz(2)
       subsizes(4) = decomp%ysz(3)
       starts(1) = 0
       starts(2) = decomp%yst(1)-1
       starts(3) = decomp%yst(2)-1
       starts(4) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = 3
       subsizes(2) = decomp%zsz(1)
       subsizes(3) = decomp%zsz(2)
       subsizes(4) = decomp%zsz(3)
       starts(1) = 0
       starts(2) = decomp%zst(1)-1
       starts(3) = decomp%zst(2)-1
       starts(4) = decomp%zst(3)-1
    endif
    
!creation du type utilise par mpi en fonction de l'endianness
!  if(testEndian2==1)then
!    
!    blockLen=1
!    do i=1,4
!      dplt(i)=4-i
!    end do
!    
!    call mpi_type_create_hindexed(4,blockLen,dplt,MPI_integer1,type_reel,ierror)
!    call MPI_Type_commit(type_reel,ierror)
!    
!  else
!    type_reel = MPI_REAL
!  end if
  type_reel = MPI_REAL
  
!creation de la structure de donnees
    call MPI_TYPE_CREATE_SUBARRAY(4, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, MPI_REAL, newtype, ierror)    
    call MPI_TYPE_COMMIT(newtype,ierror)


!ouverture du fichier    
    call MPI_FILE_OPEN(MPI_COMM_WORLD, nomFic, MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                    MPI_INFO_NULL, fh, ierror)


! deplacement ( permet de ne pas ecraser l'en-tete)   
    call MPI_File_get_size(fh,filesize,ierror)
    disp = filesize

!ecriture des donnees

   call MPI_FILE_SET_VIEW(fh,disp,type_reel, &
         newtype,'native',MPI_INFO_NULL,ierror)
   call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3)*subsizes(4), &
         type_reel, MPI_STATUS_IGNORE, ierror)

!fermeture du fichier
    call MPI_FILE_CLOSE(fh,ierror)
!    if(testEndian2==1)call MPI_TYPE_FREE(type_reel,ierror)
    call MPI_TYPE_FREE(newtype,ierror)

  end subroutine  write_bin_dep
!-------------------------------------------------------------------------------

!===============================================================================
!
!                   SUBROUTINE WRITE_HEADER_VTK_DEP : 
!> Ecriture de l'en-tete du fichier vtk dans la configuration actuelle
!!
!! \param[in]   nomFic: (chaine de characteres) nom du fichier d'ecriture
!! \param[in]         Fvtk: (entier) unite loique du fichier "nomFic"
!! \param[in]         nx,ny,nz: (entiers) dimensions de la cellule
!!
!-------------------------------------------------------------------------------
subroutine write_header_vtk_dep(nomFic, Fvtk, nx, ny, nz)

  implicit none
  character(len=*), intent(in)  :: nomFic
  integer, intent(in)           :: Fvtk, nx,ny,nz
  integer :: ios 
  integer(kind=8)               :: ntot

  ntot = nx
  ntot = ntot*ny
  ntot = ntot*nz

!ouverture du fichier de sortie
  open(unit=Fvtk, file=nomFic,form="formatted", &
       status="replace", action="write",iostat= ios)
  if ( ios /= 0 ) then
     print *," Probleme a l'ouverture (fichier: ",nomFic,")",ios
     stop 
  end if

  write(Fvtk,"(A)")"# vtk DataFile Version 4.5"
  write(Fvtk,"(A)")"Materiau"
  write(Fvtk,"(A)")"BINARY"
  write(Fvtk,"(A)")"DATASET STRUCTURED_GRID"
  write(Fvtk,"(A,I5,I5,I5)")"DIMENSIONS  ",nx,ny,nz
  write(Fvtk, "(A,I0,A)")"POINTS    ",ntot,"   float"

  close(Fvtk)


end subroutine write_header_vtk_dep

! now (03/2025) available in io_amitex_mod
!!***********************************************************************************************************************************
!!  SWAP_F4
!!
!!  Swap bytes for a four-byte REAL.
!!  After calling this subroutine, the input number will be replaced by the output number.
!!***********************************************************************************************************************************

!      ELEMENTAL SUBROUTINE SWAP_F4 (FLOAT4)

!      IMPLICIT NONE

!      REAL(KIND=4), INTENT(IN OUT) :: FLOAT4

!      INTEGER(KIND=1), DIMENSION(4) :: BYTE_ARR, BYTE_ARR_TMP
!      INTEGER :: I


!      BYTE_ARR = TRANSFER (FLOAT4, BYTE_ARR)

!      BYTE_ARR_TMP = BYTE_ARR

!      DO I = 1, 4
!         BYTE_ARR(I) = BYTE_ARR_TMP(5-I)
!      END DO

!      FLOAT4 = TRANSFER (BYTE_ARR, FLOAT4)

!      RETURN

!      END SUBROUTINE SWAP_F4



end module deformed_shape_mod
