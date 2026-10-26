!===============================================================================
!
! MODULE IO_AMITEX_MOD : 
!> Regroupement de subroutines ELEMENTAIRES liees aux entrees et aux sorties
!! 	les autres modules (material_mod, param_algo_mod, loading_mod) utilisent
!!	ces subroutines ELEMENTAIRES
!! 
!!
!!  Subroutines
!! - write_file0 :    ecriture dans un fichier par le noeud nrank=0
!! - write_stdout0 :  ecriture sur stdout par le noeud nrank=0 
!! - read_header_vtk :       ecriture de l'en-tete d'un fichier vtk
!! - write_header_vtk :      ecriture de l'en-tete d'un fichier vtk
!! - read_bin :      lecture de la partie binaire d'un fichier vtk
!! - write_bin32 :   ecriture de la partie binaire d'un fichier vtk
!! - nextInt :       recherche de l'entier suivant dans une chaine de caracteres
!! - nextReal :      recherche du reel suivant dans une chaine de caracteres
!! - get_int_vect :  lecture d'un tableau d'entier dans une chaine de caracteres
!! - getNodeList :   rechercher une liste de noeud dans un fichier xml
!! - isPresent :     verifier la presence d'un attribut dans un fichier xml
!! - get_int_xml :   lecture d'un entier dans un fichier xml
!! - get_real_xml :  lecture d'un reel dans un fichier xml
!! - get_real_xml_default :  lecture d'un reel dans un fichier xml avec 
!!                    une valeur par defaut
!! - get_str_xml :   lecture d'une chaine de caracteres dans un fichier xml
!! - get_str_xml_nolowercase : idem sans modifier la casse
!! - get_str_xml_default :  lecture d'une chaine de caracteres dans un fichier xml avec 
!!                    une valeur par defaut
!! - get_vect_xml :  lecture des coefficients ou variables internes du materiau
!!                dans un fichier xml
!! - getCoeffFromFile :     lecture de la liste des valeurs d'un coefficents 
!!               (ou d'une variables internes) dans un fichier ASCII (cas constant par
!!               zone) pour un materiau
!! - getBinCoeffFromFile :  lecture de la liste des valeurs d'un coefficents 
!!               (ou d'une variables internes) dans un fichier binaire (cas constant par
!!               zone) pour un materiau
!!
!!
!! ATTENTION : il s'agit ici d'un module de subroutines ELEMENTAIRES
!!	=> ne pas definir de subroutine plus evoluees necessitant 
!!	l'ajout de nouveau module
!!
!===============================================================================
module io_amitex_mod
  use ISO_FORTRAN_ENV
  use, intrinsic :: IEEE_ARITHMETIC
  use mpi
#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype, real_type
  use decomp_2d_mpi,       only : nrank
  use decomp_2d,           only : DECOMP_INFO
  use decomp_2d_fft,       only : decomp_2d_fft_get_sp,decomp_2d_fft_get_ph
#else
  use decomp_2d,     only : mytype, nrank, real_type, DECOMP_INFO
  use decomp_2d_fft, only : get_decomp_fft_info
#endif  
  use error_mod
  use field_mod,     only : amx_start, amx_end, amx_size
  use amitex_mpi_mod, only: amitex_mpi_comm
  use Fox_dom

  private

  ! liste des subroutines publiques  
  public :: check_file, write_file0, write_stdout0, read_bin, read_bin4,read_bin8,read_header_vtk, & 
       convert_field_for_write_bin32, write_bin32, write_header_vtk, &
       GetNodeList, get_int_xml, get_real_xml, get_real_xml_default, get_str_xml,&
       get_str_xml_default,get_str_xml_nolowercase, get_int_vect,get_real_vect, get_vect_xml, stringToLower,&
       nextReal, nextInt8, nextInt, type2mpi_type,&
       getBinCoeffFromFile, SWAP_F4, SWAP_F8, SWAP_I2, SWAP_I4, SWAP_I8

  !> Lecture de la partie binaire d'un fichier vtk 
  !!
  !! Deux implementations en fonction de la taille des entiers
  !! dans le tableau a lire.
  interface read_bin
     module procedure read_bin4
     module procedure read_bin8
     module procedure readreal_bin8
  end interface read_bin

  !> Convert a field (32 or 64 bits) in a new field (32bits bigendian) to be written on disk (vtk)
  interface convert_field_for_write_bin32
     module procedure convert_field32_for_write_bin32
     module procedure convert_field64_for_write_bin32
  end interface convert_field_for_write_bin32
contains

!===============================================================================
!
!                   SUBROUTINE write_stdout0
!
!> Ecrit une chaine de caractere sur stdout par le noeud nrank=0
!!
!! \param[in]   string: chaine de caractere a ecrire
!!
!-------------------------------------------------------------------------------
subroutine write_stdout0(string)
   
   implicit none

   character(len=*),intent(in)    :: string

   if(nrank==0) write(OUTPUT_UNIT,"(A)") trim(string)

end subroutine write_stdout0
!===============================================================================
!
!                   SUBROUTINE check_file
!
!> verifie que le nom de fichier file_name existe sinon amitex_abort(err_msg,2,0)
!!
!! \param[in]   file_name     : nom du fichier a tester
!! \param[in]   function_name : nom de la fonction appelante 
!!
!!
!! ATTENTION : interruption immediate, seul le process 0 ecrit le message d'erreur
!!             cette fonction doit etre appele a minima par le process 0
!-------------------------------------------------------------------------------
subroutine check_file(file_name,function_name)
   
   implicit none

   character(len=*),intent(in)    :: file_name
   character(len=*),intent(in)    :: function_name
   character(len=500)             :: err_msg
   logical                        :: file_exists

   inquire(FILE=file_name,EXIST=file_exists)
   if(.not.file_exists)then 
        write(err_msg,fmt="(5A)") "Le fichier: '",trim(file_name),"' n'existe pas (",trim(function_name),")"
        call amitex_abort(trim(err_msg),2,0)
   end if

end subroutine check_file
!===============================================================================
!
!                   SUBROUTINE write_file0
!
!> Ecrit une chaine de caractere sur un fichier par le noeud nrank=0
!!
!! \param[in]   FU    : file unit
!! \param[in]   string: chaine de caractere a ecrire
!!
!-------------------------------------------------------------------------------
subroutine write_file0(FU,string)
   
   implicit none

   integer,intent(in)             :: FU
   character(len=*),intent(in)    :: string

   if(nrank==0) write(FU,"(A)") trim(string)

end subroutine write_file0

!===============================================================================
!
!                   SUBROUTINE TYPE2MPI_TYPE
!> Traduit un nom de type en type MPI
!!
!! \param[in]   type_name: (chaine de characteres) nom du type de donnees
!! \param[out]  mpi_type: type de donnees au format MPI_TYPE
!!
!-------------------------------------------------------------------------------
subroutine type2mpi_type(type_name,mpi_type) 
  
  implicit none
  
  character(len=*), intent(in)  :: type_name
  integer, intent(out)          :: mpi_type
  character(len=150)            :: err_msg
  
  if(index(type_name,"char")>0)then
     mpi_type = MPI_INTEGER1
  else if(index(type_name,"unsigned_char")>0) then
     mpi_type = MPI_INTEGER1
     write(err_msg, fmt="(A)") "Attention !!! unsigned_char &
          & Les donnees du fichier doivent etre comprises entre 0 et 127"
     call amitex_abort(err_msg,0,0)
  else if(index(type_name,"short")>0) then
     mpi_type = MPI_INTEGER2
  else if(index(type_name,"unsigned_short")>0) then
     mpi_type = MPI_INTEGER2
     write(err_msg, fmt="(A)") "Attention !!! unsigned_short &
          & Les donnees du fichier doivent etre comprises entre 0 et 2^15 -1"
     call amitex_abort(err_msg,0,0)
  else if(index(type_name,"int")>0) then
     mpi_type = MPI_INTEGER4
  else if(index(type_name,"unsigned_int")>0) then
     mpi_type = MPI_INTEGER4
     write(err_msg, fmt="(A)") "Attention !!! unsigned_int &
          & Les donnees du fichier doivent etre comprises entre 0 et 2^31 -1"
     call amitex_abort(err_msg,0,0)
  else if(index(type_name,"long")>0) then
     mpi_type = MPI_INTEGER8
  else if(index(type_name,"unsigned_long")>0) then
     mpi_type = MPI_INTEGER8
     write(err_msg, fmt="(A)") "Attention !!! unsigned_long &
          & Les donnees du fichier doivent etre comprises entre 0 et 2^63 -1"
     call amitex_abort(err_msg,0,0)
  else if(index(type_name,"float")>0) then
     mpi_type = MPI_REAL
  else if(index(type_name,"double")>0) then
     mpi_type = MPI_DOUBLE_PRECISION
  else
     write(err_msg, fmt="(A)") "Type de donnees non conforme. (type2mpi_type)"
     call amitex_abort(err_msg,1,0)
  end if
end subroutine type2mpi_type
!===============================================================================
!
!        SUBROUTINE READ_HEADER_VTK : 
!> Lecture de l'en-tete du fichier vtk
!!
!! Arrete le programme si l'en-tete n'est pas conforme ou si 
!! les donnees ne sont pas des entiers (sur 2 ou 4 octets).
!!
!! \param[in]   nomFic        (chaine de characteres)nom du fichier a lire
!! \param[out]  nx,ny,nz      (entiers) DIMENSIONS in vtk file
!! \param[out]  dx,dy,dz      SPACING in vtk file
!! \param[out]  x0,y0,z0      ORIGIN in vtk file
!! 
!! \param[out]  type_mpi      type de donnees contenu dans le fichier \n
!!                            entier sur 1, 2 ou 4 octets :\n
!!                            MPI_INTEGER1, MPI_INTEGER2 ou MPI_INTEGER4
!!
!!------------------------------------------------------------------------------
!! MODIFICATIONS :
!!----------------
!! 21/01/2016 - LG : analyse ligne par ligne de l'entete
!!                   passage de tous les amitex_abort(err_msg, 1, 0)
!!                   en amitex_abort(err_msg,1) sinon les messages n'apparaisent 
!!                   pas (mal compris)
!!
!! 21/06/2019 - LG : relaxe l'ordre des lignes pour pouvoir directement lire
!!                   les fichiers vtk issus de python (+ import vtk)
!!
!! 16/04/2021 - LG : ajout des sorties x0,y0,z0 
!!
!! 30/06/2021 - LG : correction bug (M.Josien : pb lecture '5e-7', '5.e-7' ok)
!!                   suppression utilisation nextReal
!!                   lecture des triplets par "read(chaine,*,IOSTAT=ios) data"
!!                           -> utilisation format libre
!!                              + tests ios et affectation des datas
!-------------------------------------------------------------------------------
subroutine read_header_vtk(nomFic, nx,ny,nz,dx,dy,dz,x0,y0,z0, type_mpi)

  implicit none

  character(len=*), intent(in)   :: nomFic
  integer                        :: FU

  integer, intent(out)           :: nx,ny,nz, type_mpi
  real(mytype), intent(out)      :: dx,dy,dz
  real(mytype), intent(out)      :: x0,y0,z0
  
  integer                        :: i, p, ios, n
  character(len=150)             :: err_msg
  integer(kind=8)                :: cell_data, cn
  logical                        :: file_exists
  logical,dimension(8)           :: tab_valid 
  character(len=200), dimension(10)::entete
 
  n=0
  cell_data=0
  tab_valid = .false.
  dx = huge(dx); dy=dx; dz=dx
  nx = huge(nx); ny=nx; nz=nx
  x0 = huge(x0); y0=x0; z0=x0

  !! Ouverture du fichier d'entree
  INQUIRE(file=trim(nomFic),exist=file_exists)
  if(.not.file_exists) then
    write(err_msg,fmt="(3A)") "Le fichier fichier : ",&
         trim(nomFic)," n'existe pas (read_header_vtk)"
    call amitex_abort(err_msg,2,0)
  end if
  open(newunit=FU, file=trim(nomFic),form="formatted", status="old", action="read",iostat= ios)
  if ( ios /= 0 ) then
     write(err_msg,fmt="(3A,I0,A)") "Probleme a l'ouverture du fichier : ",trim(nomFic),&
          " , ",ios," (read_header_vtk)"
     call amitex_abort(err_msg,2,0) 
  end if

  !! Lecture et stockage de l'en-tete (les 10 premieres lignes) dans le tableau entete
  do i=1,10
      read(FU,FMT='(A)') entete(i)
  end do
  !! Fermeture du fichier d'entree
  close(unit=FU)  


  !!-------------Boucle sur les lignes de l'entete
  do i=1,10

  !1 lecture du mot BINARY
  if(index(entete(i),"BINARY")>0 .AND. i==3)then
      n=n+1
      tab_valid(1)=.true.
  end if

  !2 lecture du mot DATASET STRUCTURED_POINTS
  if(index(entete(i),"DATASET STRUCTURED_POINTS")>0 .AND. i==4)then
      n=n+1
      tab_valid(2)=.true.
  end if

  !3 lecture des DIMENSIONS
  if (index(entete(i),"DIMENSIONS")>0) then
      p=10
      !! Lecture des dimensions de la cellule
      !! Subroutine nextInt: erreur de lecture => p=-1
      !! En cas d'erreur: n=-1
      !call nextInt(entete(i),p,nx)
      !nx=nx-1
      !if(p==-1) n=-1
      !call nextInt(entete(i),p,ny)
      !ny=ny-1
      !if(p==-1) n=-1
      !call nextInt(entete(i),p,nz)
      !nz=nz-1
      !if(p==-1) n=-1
      !if(n == -1) then
      !  write(err_msg, fmt="(A)") "Invalid or missing DIMENSIONS (read_header_vtk)"
      !  call amitex_abort(err_msg,1)
      !! else if((nx/2)*2==nx .or. (ny/2)*2==ny .or.(nz/2)*2==nz) then
      !!   write(err_msg, fmt="(A,3I5)") "les dimensions paires ne sont pas encore traitees. (read_header_vtk)"
      !!   call amitex_abort(err_msg,1,0)
      !end if
      
      ! utilisation de la fonction de lecture fortran!!!
      read(entete(i),*,IOSTAT=ios) err_msg,nx,ny,nz
      if (ios > 0) call amitex_abort("Problem in (read_header_vtk), line DIMENSIONS",2)
      nx = nx - 1
      ny = ny - 1
      nz = nz - 1
      if(nx<=0 .or. ny<=0 .or. nz<=0) then
        write(err_msg, fmt="(A)") "DIMENSION nul or negative (read_header_vtk)"
        call amitex_abort(err_msg,1)
      else if(nx==huge(nx) .or. ny==huge(nx) .or. nz==huge(nx)) then
        write(err_msg, fmt="(A)") "Unaffected DIMENSION (read_header_vtk)"
        call amitex_abort(err_msg,1)      
      end if

      n=n+1
      tab_valid(3)=.true.
  end if

  !4 lecture des SPACING
  if (index(entete(i),"SPACING")>0) then
      p=7
      ! lecture des dimensions des voxels
      ! en cas d'erreur: n=-1
      !call nextReal(entete(i),p,dx)
      !if(p==-1) n=-1
      !call nextReal(entete(i),p,dy)
      !if(p==-1) n=-1
      !call nextReal(entete(i),p,dz)
      !if(p==-1) n=-1
      !if(n == -1) then
      !  write(err_msg, fmt="(A)") "Invalid or missing SPACING (read_header_vtk)"
      !  call amitex_abort(err_msg,1)
      !else if(dx==0 .or. dy==0 .or. dz==0) then
      !  write(err_msg, fmt="(A)") "Spacing nul (read_header_vtk)"
      !  call amitex_abort(err_msg,1)
      !end if
      
      ! utilisation de la fonction de lecture fortran!!!
      read(entete(i),*,IOSTAT=ios) err_msg,dx,dy,dz
      if (ios > 0) call amitex_abort("Problem in (read_header_vtk), line SPACING",2)

      if(dx<=0 .or. dy<=0 .or. dz<=0) then
        write(err_msg, fmt="(A)") "SPACING nul or negative (read_header_vtk)"
        call amitex_abort(err_msg,1)
      else if(dx==huge(dx) .or. dy==huge(dx) .or. dz==huge(dx)) then
        write(err_msg, fmt="(A)") "Unaffected SPACING (read_header_vtk)"
        call amitex_abort(err_msg,1)      
      end if
      
      n=n+1
      tab_valid(4)=.true.
  end if

  !5 lecture du nombre de cellules CELL_DATA    
  if (index(entete(i),"CELL_DATA")>0) then
      p=9
      call nextInt8(entete(i),p,cell_data)
      n=n+1
      tab_valid(5)=.true.
  end if
  
  !6 lecture du type de donnees SCALARS
  if(index(entete(i),"SCALARS")>0 .AND. i==9) then
    !! Recuperation du type de donnees mpi correspondant au type de donnees vtk
    call type2mpi_type(trim(entete(i)),type_mpi)
    !! On verifie que l'on recupere bien des entiers codes sur moins de 8 octets
    if(type_mpi /= MPI_INTEGER1 .AND. type_mpi /= MPI_INTEGER2 &
         .AND. type_mpi /= MPI_INTEGER4 .AND. type_mpi /= MPI_INTEGER8 &
         .AND. type_mpi /= MPI_REAL .AND. type_mpi /= MPI_DOUBLE_PRECISION) then
      write(err_msg, fmt="(A)") "Non conform integer type (read_header_vtk)"
      call amitex_abort(err_msg,1,0)
    end if
    n=n+1
    tab_valid(6)=.true.
  end if
  
  !7 lecture du mot LOOKUP_TABLE
  if (index(entete(i),"LOOKUP_TABLE")>0 .AND. i==10)then
      n=n+1
      tab_valid(7)=.true.
  end if
  
  !8 lecture des coordonnees ORIGIN
  if (index(entete(i),"ORIGIN")>0) then
      p=6
      ! lecture des coordonnees de l'origine de la grille (coin mini)
      ! en cas d'erreur: n=-1
      !call nextReal(entete(i),p,x0)
      !if(p==-1) n=-1
      !call nextReal(entete(i),p,y0)
      !if(p==-1) n=-1
      !call nextReal(entete(i),p,z0)
      !if(p==-1) n=-1
      !if(n == -1) then
      !  write(err_msg, fmt="(A)") "Invalid or missing ORIGIN (read_header_vtk)"
      !  call amitex_abort(err_msg,1)
      !end if
      
      ! utilisation de la fonction de lecture fortran!!!
      read(entete(i),*,IOSTAT=ios) err_msg,x0,y0,z0
      if (ios > 0) call amitex_abort("Problem in (read_header_vtk), line ORIGIN",2)

      if(x0==huge(x0) .or. y0==huge(x0) .or. z0==huge(x0)) then
        write(err_msg, fmt="(A)") "Unaffected ORIGIN (read_header_vtk)"
        call amitex_abort(err_msg,1)      
      end if
  
      n=n+1
      tab_valid(8)=.true.
  end if



  end do !-- fin boucle sur les lignes de l'entete



  !verification de la coherence : dimensions-cell_data
  cn=nx*ny
  cn=cn*nz
  if (cn/= cell_data) then
    write(err_msg, fmt="(A,I0,A,3(I0,A))") "La taille de la cellule (CELL_DATA=",cell_data,&
         ") et ses dimensions ",nx,"x", ny,"x",nz, " ne sont pas coherentes (read_header_vtk)"
    call amitex_abort(err_msg,1)
  end if

  !------------------------- VERIFICATIONS
  if (n/=8) then
    write(err_msg,fmt="(3A)") "En-tete du fichier vtk '",trim(nomFic),"' mal forme (read_header_vtk)"
    call amitex_abort(err_msg,1)
  end if

  if (.not. tab_valid(1)) then
      write(err_msg, fmt="(A)") "BINARY absent (read_header_vtk)"
      call amitex_abort(err_msg,1)
  end if 
  if (.not. tab_valid(2)) then 
      write(err_msg, fmt="(A)") "DATASET STRUCTURED_POINTS absent (read_header_vtk)"
      call amitex_abort(err_msg,1)
  end if 
  if (.not. tab_valid(3)) then
      write(err_msg, fmt="(A)") "DIMENSIONS absent (read_header_vtk)"
      call amitex_abort(err_msg,1)
  end if 
  if (.not. tab_valid(4)) then
      write(err_msg, fmt="(A)") "SPACING non présent (read_header_vtk)"
      call amitex_abort(err_msg,1)
  end if 
  if (.not. tab_valid(5)) then
      write(err_msg, fmt="(A)") "CELL_DATA absent (read_header_vtk)"
      call amitex_abort(err_msg,1)
  end if 
  if (.not. tab_valid(6)) then
      write(err_msg, fmt="(A)") "SCALARS absent (read_header_vtk)"
      call amitex_abort(err_msg,1)
  end if 
  if (.not. tab_valid(7)) then 
      write(err_msg, fmt="(A)") "LOOKUP_TABLE absent (read_header_vtk)"
      call amitex_abort(err_msg,1)
  end if 
  if (.not. tab_valid(8)) then 
      write(err_msg, fmt="(A)") "ORIGIN absent (read_header_vtk)"
      call amitex_abort(err_msg,1)
  end if 


  call check_amitex_abort()

end subroutine read_header_vtk
!-------------------------------------------------------------------------------


!===============================================================================
!
!                   SUBROUTINE WRITE_HEADER_VTK : 
!> Ecriture de l'entete du fichier vtk par le noeud 0
!!
!! \param[in]         nomFic: (chaine de characteres) nom du fichier d'ecriture
!! \param[in]         nx,ny,nz: (entiers) dimensions de la cellule
!! \param[in]         dx,dy,dz: (reels) espacement dans la cellule
!! \param[in]         x0,y0,z0: (reels) origine de la cellule
!!
!-------------------------------------------------------------------------------
subroutine write_header_vtk(nomFic,nx,ny,nz,dx,dy,dz,x0,y0,z0)

  implicit none
  character(len=*), intent(in)  :: nomFic
  integer, intent(in)           :: nx,ny,nz
  integer                       :: Fvtk 
  real(mytype), intent(in)      :: dx,dy,dz
  real(mytype), intent(in)      :: x0,y0,z0
  integer :: ios 
  integer(kind=8)               :: ntot

  ntot = nx
  ntot = ntot*ny
  ntot = ntot*nz

  !ecriture sur le fichier de sortie (proc 0)
  if (nrank == 0) then
    open(newunit=Fvtk, file=nomFic,form="formatted", &
                              status="replace", action="write",iostat= ios)
      if ( ios /= 0 ) then
        write(OUTPUT_UNIT,"(3A,I4)") " Problem when opening (file : ",nomFic,") - ios = ",ios
        stop 
      end if
    
    write(Fvtk,"(A)")"# vtk DataFile Version 4.5"
    write(Fvtk,"(A)")"Materiau"
    write(Fvtk,"(A)")"BINARY"
    write(Fvtk,"(A)")"DATASET STRUCTURED_POINTS"
    write(Fvtk,"(A,I6,I6,I6)")"DIMENSIONS  ",nx+1,ny+1,nz+1
    write(Fvtk,"(A,3E15.7)")"ORIGIN     ",x0,y0,z0
    write(Fvtk,"(A,3E15.7)")"SPACING    ",dx,dy,dz
    write(Fvtk,"(A,I12)")"CELL_DATA    ",ntot
  close(Fvtk)
 end if

end subroutine write_header_vtk
!-------------------------------------------------------------------------------




!===============================================================================
!===============================================================================
!
!
! FONCTIONS DE LECTURE ET D'ECRITURE DE LA PARTIE BINAIRE DE FICHIER VTK
!
!
!> Les donnees lues et ecrites sont des tableaux a trois dimensions 
!! (repartis avec 2decomp). \n
!!
!! L'ecriture se fait a la fin du fichier specifie. \n
!! La lecture se fait sur les nx.ny.nz.2(ou 4) derniers octets du fichier
!! a partir du dernier saut de ligne detecte. \n
!! Ce saut de ligne est repere par l'octet 00001010 (valeurs decimale: 10)
!!  et ou: 
!!      - nx,ny,nz sont les dimensions du tableau
!!      - 2 (ou 4) est le nombre d'octets sur lesquels sont codes les donnees
!!     (il s'agit de integer(kind=2) ou integer(kind=4))
!! La fonction read_bin4 est utilisee pour lire les numeros de materiaux (kind = 4)
!! et la fonction read_bin8 est utilisee pour lire les numeros de zones (kind = 8)
!!
!!!! CHANGEMENT 14/12/2018 : inversion ordre de byte APRES lecture MPI 
!!             => generalisation de la modif faite a getbinCoeffFromFile
!!                meme si ici, l'effet est moindre (mais positif)
!!
!===============================================================================
!===============================================================================

!===============================================================================
!
!      SUBROUTINE READ_BIN : 
!> Lecture de la partie binaire d'un fichier vtk 
!!
!!
!! \param[in]       nomFic: nom du fichier
!! \param[in]       ipencil: orientation des pinceaux du tableau var
!! \param[in]       type_mpi: type des donnees lues (MPI_INTEGER2 ou MPI_INTEGER4)
!! \param[out]      var: tableau a lire 
!-------------------------------------------------------------------------------
  subroutine read_bin4(nomFic,ipencil,type_mpi,var)
    
    implicit none

    character(len=*), intent(in)                     :: nomFic
    integer, intent(in)                              :: ipencil,type_mpi
    integer(kind=4), dimension(1:amx_size(1)*amx_size(2)*amx_size(3)), intent(out)   :: var

    integer(kind=2), allocatable, dimension(:)       :: tmp_int2
    integer(kind=1), allocatable, dimension(:)       :: tmp_int1
#ifdef NEW_2DECOMP
    TYPE(DECOMP_INFO),pointer                        :: decomp, sp_decomp
#else
    TYPE(DECOMP_INFO)                                :: decomp, sp_decomp
#endif
    integer(kind=MPI_OFFSET_KIND)                    :: disp, filesize
    integer, dimension(3)                            :: sizes, subsizes, starts
    integer                                          :: ierror, newtype, fh, nOctets,i
    integer(kind=1)                                  :: debut


    !variable du traitement de l'endianness
    integer(kind=1), dimension(2)::testEndian
    integer (kind=2)             ::testEndian2
    equivalence(testEndian,testEndian2)

    !variable pour definir un type de donnees adapte
    integer :: type_end
!    integer,dimension(:),allocatable ::blockLen
!    integer(kind=MPI_ADDRESS_KIND),dimension(:), allocatable :: dplt

    testEndian=(/int(1,kind=1),int(0,kind=1)/)
#ifdef NEW_2DECOMP
    decomp    => decomp_2d_fft_get_ph()
    sp_decomp => decomp_2d_fft_get_sp()
#else
    call get_decomp_fft_info(decomp,sp_decomp)
#endif
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

    !! nOctets = nombre d'octets par entier en fonction du type d'entier
    call MPI_TYPE_SIZE(type_mpi,nOctets,ierror)
    if(type_mpi /= MPI_INTEGER1 .AND. type_mpi /= MPI_INTEGER2 &
         .AND. type_mpi /= MPI_INTEGER4) then
      call amitex_abort("Type de donnees non conforme. (read_bin4)",1,0)
    end if

    !! Creation du type utilise par mpi en fonction de l'endianness (modif le 14/12/2018)
!    if(testEndian2==1)then
!       allocate(blockLen(nOctets),dplt(nOctets))
!       blockLen=1
!       do i=1,nOctets
!          dplt(i)=nOctets-i
!       end do
!       !creation de la structure d'entier 
!       call mpi_type_create_hindexed(int(nOctets,kind=4),blockLen,dplt,MPI_integer1,type_end,ierror)
!       call MPI_Type_commit(type_end,ierror)
!    else
!       type_end = type_mpi
!    end if
    type_end=type_mpi

    !creation de la structure de donnee (tableau d'entier defini ci-dessus)
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, type_mpi, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)

    call MPI_FILE_OPEN(amitex_mpi_comm, nomFic, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, &
         fh, ierror)

    call MPI_BARRIER(amitex_mpi_comm,ierror)
    !recherche du debut de l'enregistrement des donnees
    !estimation du debut de l'enregistrement
    call MPI_File_get_size(fh,filesize,ierror)

    disp = sizes(1)
    disp = disp*sizes(2)
    disp = disp*sizes(3)
    disp = disp*int(nOctets,kind=8)
    disp = (filesize) - disp
    if(disp <1)then  
       call amitex_abort( "Erreur lors de la lecture du fichier vtk (read_bin4) " ,2)
    end if

    call MPI_File_seek(fh,disp,MPI_SEEK_SET ,ierror)
    call MPI_File_read(fh, debut, 1, MPI_INTEGER1, MPI_STATUS_IGNORE, ierror)

    !! Recherche du premier saut de ligne avant le debut estime des donnees
    !! (10 en ASCII correspond a un saut de ligne)
    i=0
    do while( (debut/=10) )
       i=i+int(1,kind=1)
       disp = disp -1
       call MPI_File_seek(fh,disp,MPI_SEEK_SET ,ierror)
       call MPI_File_read(fh, debut, 1, MPI_INTEGER1, MPI_STATUS_IGNORE, ierror)
    end do
    !si les caracteres precedent sont aussi 10, on remonte au 1er
    ! ajout LG : 07/12/2020, correction lecture si 1er entier des donnees est un 10
    do while( (debut==10) )
       disp = disp -1
       call MPI_File_seek(fh,disp,MPI_SEEK_SET ,ierror)
       call MPI_File_read(fh, debut, 1, MPI_INTEGER1, MPI_STATUS_IGNORE, ierror)
    end do
    disp = disp + 1

    !On place le debut apres le saut de ligne trouve (ci-dessus)
    disp = disp +1

    call MPI_FILE_SET_VIEW(fh,disp,type_end, &
         newtype,'native',MPI_INFO_NULL,ierror)

    !lecture des donnees (selon le type)
    if (nOctets==1) then
       !lecture des donnees (integer(kind=1))
       allocate(tmp_int1(1:amx_size(1)*amx_size(2)*amx_size(3)))
       call MPI_FILE_READ_ALL(fh, tmp_int1, &
            subsizes(1)*subsizes(2)*subsizes(3), &
            type_end, MPI_STATUS_IGNORE, ierror)
       var=int(tmp_int1,4)
       deallocate(tmp_int1)
    else if (nOctets==2) then
       !lecture des donnees (integer(kind=2))
       allocate(tmp_int2(1:amx_size(1)*amx_size(2)*amx_size(3)))
       call MPI_FILE_READ_ALL(fh, tmp_int2, &
            subsizes(1)*subsizes(2)*subsizes(3), &
            type_end, MPI_STATUS_IGNORE, ierror)
       if(testEndian2==1)call SWAP_I2(tmp_int2)
       var=int(tmp_int2,4)
       deallocate(tmp_int2)
    elseif (nOctets==4) then
       !lecture des donnees (integer(kind=4))
       call MPI_FILE_READ_ALL(fh, var, &
            subsizes(1)*subsizes(2)*subsizes(3), &
            type_end, MPI_STATUS_IGNORE, ierror)
       if(testEndian2==1)call SWAP_I4(var)
    end if

    !fermeture du fichier
    call MPI_FILE_CLOSE(fh,ierror)
    !destruction des type et tableau
    call MPI_TYPE_FREE(newtype,ierror)
!    if(testEndian2==1) call MPI_TYPE_FREE(type_end,ierror)
!    if(allocated(blockLen)) deallocate(blockLen)
!    if(allocated(dplt)) deallocate(dplt)

    return
  end subroutine read_bin4
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  subroutine read_bin8(nomFic,ipencil,type_mpi,var)

    implicit none

    character(len=*), intent(in)                     :: nomFic
    integer, intent(in)                              :: ipencil,type_mpi
    integer(kind=8), dimension(1:amx_size(1)*amx_size(2)*amx_size(3)), intent(out) :: var

    integer(kind=1), allocatable, dimension(:)   :: tmp_int1
    integer(kind=2), allocatable, dimension(:)   :: tmp_int2
    integer(kind=4), allocatable, dimension(:)   :: tmp_int4
#ifdef NEW_2DECOMP
    TYPE(DECOMP_INFO),pointer                        :: decomp, sp_decomp
#else
    TYPE(DECOMP_INFO)                                :: decomp, sp_decomp
#endif
    integer(kind=MPI_OFFSET_KIND)                    :: disp, filesize
    integer, dimension(3)                            :: sizes, subsizes, starts
    integer                                          :: ierror, newtype, fh, nOctets,i
    integer(kind=1)                                  :: debut


    !variable du traitement de l'endianness
    integer(kind=1), dimension(2)::testEndian
    integer (kind=2)             ::testEndian2
    equivalence(testEndian,testEndian2)

    !variable pour definir un type de donnees adapte
    integer :: type_end
!    integer,dimension(:),allocatable ::blockLen
!    integer(kind=MPI_ADDRESS_KIND),dimension(:), allocatable :: dplt

    testEndian=(/int(1,kind=1),int(0,kind=1)/)

#ifdef NEW_2DECOMP
    decomp    => decomp_2d_fft_get_ph()
    sp_decomp => decomp_2d_fft_get_sp()
#else
    call get_decomp_fft_info(decomp,sp_decomp)
#endif

    !tailles, debut et etendues des sous tableau
    !starts ramene a 0 (pour MPI) 
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)

    !taille et bornes des tableau selon l'orientation:
    !pinceaux selon x, y ou z   
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

    !! nOctets = nombre d'octets par entier en fonction du type d'entier
    call MPI_TYPE_SIZE(type_mpi,nOctets,ierror)
    if(type_mpi /= MPI_INTEGER1 .AND. type_mpi /= MPI_INTEGER2 &
         .AND. type_mpi /= MPI_INTEGER4 .AND. type_mpi /= MPI_INTEGER8) then
      call amitex_abort("Type de donnees non conforme. (read_bin8)",1,0)
    end if

    !! Creation du type utilise par mpi en fonction de l'endianness (modif le 14/12/2018)
!    if(testEndian2==1)then
!       allocate(blockLen(nOctets),dplt(nOctets))
!       blockLen=1
!       do i=1,nOctets
!          dplt(i)=nOctets-i
!       end do
!       call mpi_type_create_hindexed(int(nOctets,kind=4),blockLen,dplt,MPI_integer1,type_end,ierror)
!       call MPI_Type_commit(type_end,ierror)
!    else
!       type_end = type_mpi
!    end if
     type_end = type_mpi

    !! Creation de la structure de donnees (tableau d'entiers defini ci-dessus)

    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, type_mpi, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(amitex_mpi_comm, nomFic, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, &
         fh, ierror)
    call MPI_BARRIER(amitex_mpi_comm,ierror)

    !recherche du debut de l'enregistrement des donnees
    !estimation du debut de l'enregistrement
    call MPI_File_get_size(fh,filesize,ierror)

    disp = sizes(1)
    disp = disp*sizes(2)
    disp = disp*sizes(3)
    disp = disp*int(nOctets,kind=8)
    disp = (filesize) - disp
    if(disp <1)then  
       call amitex_abort("Erreur lors de la lecture du fichier vtk (read_bin8) ",2)
    end if

    call MPI_File_seek(fh,disp,MPI_SEEK_SET ,ierror)
    call MPI_File_read(fh, debut, 1, MPI_INTEGER1, MPI_STATUS_IGNORE, ierror)
    !! Recherche du premier saut de ligne avant le debut estime des donnees
    !! (10 en ASCII correspond a un saut de ligne)
    i=0
    do while( (debut/=10) )
       i=i+int(1,kind=1)
       disp = disp -1
       call MPI_File_seek(fh,disp,MPI_SEEK_SET ,ierror)
       call MPI_File_read(fh, debut, 1, MPI_INTEGER1, MPI_STATUS_IGNORE, ierror)
    end do
    !si les caracteres precedent sont aussi 10, on remonte au 1er
    ! ajout LG : 07/12/2020, correction lecture si 1er entier des donnees est un 10
    do while( (debut==10) )
       disp = disp -1
       call MPI_File_seek(fh,disp,MPI_SEEK_SET ,ierror)
       call MPI_File_read(fh, debut, 1, MPI_INTEGER1, MPI_STATUS_IGNORE, ierror)
    end do
    disp = disp + 1
    
    !On place le debut apres le saut de ligne trouve (ci-dessus)
    disp = disp +1
    call MPI_FILE_SET_VIEW(fh,disp,type_end, &
         newtype,'native',MPI_INFO_NULL,ierror)
    !lecture des donnees (selon le type)
    if (nOctets==1) then
       !lecture des donnees (integer(kind=1))
       allocate(tmp_int1(1:amx_size(1)*amx_size(2)*amx_size(3)))
       call MPI_FILE_READ_ALL(fh, tmp_int1, &
            subsizes(1)*subsizes(2)*subsizes(3), &
            type_end, MPI_STATUS_IGNORE, ierror)
       var=int(tmp_int1,8)
       deallocate(tmp_int1)
    else if (nOctets==2) then
       !lecture des donnees (integer(kind=2))
       allocate(tmp_int2(1:amx_size(1)*amx_size(2)*amx_size(3)))
       call MPI_FILE_READ_ALL(fh, tmp_int2, &
            subsizes(1)*subsizes(2)*subsizes(3), &
            type_end, MPI_STATUS_IGNORE, ierror)
       if(testEndian2==1)call SWAP_I2(tmp_int2)
       var=int(tmp_int2,8)
       deallocate(tmp_int2)
    else if (nOctets==4) then
       !lecture des donnees (integer(kind=4))
       allocate(tmp_int4(1:amx_size(1)*amx_size(2)*amx_size(3)))
       call MPI_FILE_READ_ALL(fh, tmp_int4, &
            subsizes(1)*subsizes(2)*subsizes(3), &
            type_end, MPI_STATUS_IGNORE, ierror)
       if(testEndian2==1)call SWAP_I4(tmp_int4)
       var=int(tmp_int4,8)
       deallocate(tmp_int4)
    else if(nOctets==8) then
       !lecture des donnees (integer(kind=8))
       call MPI_FILE_READ_ALL(fh, var, &
            subsizes(1)*subsizes(2)*subsizes(3), &
            type_end, MPI_STATUS_IGNORE, ierror)
       if(testEndian2==1)call SWAP_I8(var)
    end if
    !fermeture du fichier
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)

!    if(testEndian2==1) call MPI_TYPE_FREE(type_end,ierror)
!    if(allocated(blockLen)) deallocate(blockLen)
!    if(allocated(dplt)) deallocate(dplt)

    return
  end subroutine read_bin8 
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine readreal_bin8(nomFic,ipencil,type_mpi,var)

    implicit none

    character(len=*), intent(in)                     :: nomFic
    integer, intent(in)                              :: ipencil,type_mpi
    real(kind=8), dimension(1:amx_size(1)*amx_size(2)*amx_size(3)), intent(out)      :: var

    real(kind=4), allocatable, dimension(:)          :: tmp_real4
#ifdef NEW_2DECOMP
    TYPE(DECOMP_INFO),pointer                        :: decomp, sp_decomp
#else
    TYPE(DECOMP_INFO)                                :: decomp, sp_decomp
#endif
    integer(kind=MPI_OFFSET_KIND)                    :: disp, filesize
    integer, dimension(3)                            :: sizes, subsizes, starts
    integer                                          :: ierror, newtype, fh, nOctets,i
    integer(kind=1)                                  :: debut


    !variable du traitement de l'endianness
    integer(kind=1), dimension(2)::testEndian
    integer (kind=2)             ::testEndian2
    equivalence(testEndian,testEndian2)

    !variable pour definir un type de donnees adapte
    integer :: type_end
!    integer,dimension(:),allocatable ::blockLen
!    integer(kind=MPI_ADDRESS_KIND),dimension(:), allocatable :: dplt

    testEndian=(/int(1,kind=1),int(0,kind=1)/)

#ifdef NEW_2DECOMP
    decomp    => decomp_2d_fft_get_ph()
    sp_decomp => decomp_2d_fft_get_sp()
#else
    call get_decomp_fft_info(decomp,sp_decomp)
#endif

    !tailles, debut et etendues des sous tableau
    !starts ramene a 0 (pour MPI) 
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)

    !taille et bornes des tableau selon l'orientation:
    !pinceaux selon x, y ou z   
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

    !! nOctets = nombre d'octets par reel en fonction du type d'entier
    call MPI_TYPE_SIZE(type_mpi,nOctets,ierror)
    if(type_mpi /= MPI_REAL .AND. type_mpi /= MPI_DOUBLE_PRECISION) then
      call amitex_abort("Type de donnees non conforme. (readreal_bin8)",1,0)
    end if

    !! Creation du type utilise par mpi en fonction de l'endianness (modif le 14/12/2018)
!    if(testEndian2==1)then
!       allocate(blockLen(nOctets),dplt(nOctets))
!       blockLen=1
!       do i=1,nOctets
!          dplt(i)=nOctets-i
!       end do
!       call mpi_type_create_hindexed(int(nOctets,kind=4),blockLen,dplt,MPI_integer1,type_end,ierror)
!       call MPI_Type_commit(type_end,ierror)
!    else
!       type_end = type_mpi
!    end if
     type_end = type_mpi

    !! Creation de la structure de donnees (tableau d'entiers defini ci-dessus)

    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, type_mpi, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(amitex_mpi_comm, nomFic, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, &
         fh, ierror)
    call MPI_BARRIER(amitex_mpi_comm,ierror)

    !recherche du debut de l'enregistrement des donnees
    !estimation du debut de l'enregistrement
    call MPI_File_get_size(fh,filesize,ierror)

    disp = sizes(1)
    disp = disp*sizes(2)
    disp = disp*sizes(3)
    disp = disp*int(nOctets,kind=8)
    disp = (filesize) - disp
    if(disp <1)then  
       call amitex_abort("Erreur lors de la lecture du fichier vtk (readreal_bin8) ",2)
    end if

    call MPI_File_seek(fh,disp,MPI_SEEK_SET ,ierror)
    call MPI_File_read(fh, debut, 1, MPI_INTEGER1, MPI_STATUS_IGNORE, ierror)
    !! Recherche du premier saut de ligne avant le debut estime des donnees
    !! (10 en ASCII correspond a un saut de ligne)
    i=0
    do while( (debut/=10) )
       i=i+int(1,kind=1)
       disp = disp -1
       call MPI_File_seek(fh,disp,MPI_SEEK_SET ,ierror)
       call MPI_File_read(fh, debut, 1, MPI_INTEGER1, MPI_STATUS_IGNORE, ierror)
    end do
    !si les caracteres precedent sont aussi 10, on remonte au 1er
    ! ajout LG : 07/12/2020, correction lecture si 1er entier des donnees est un 10
    do while( (debut==10) )
       disp = disp -1
       call MPI_File_seek(fh,disp,MPI_SEEK_SET ,ierror)
       call MPI_File_read(fh, debut, 1, MPI_INTEGER1, MPI_STATUS_IGNORE, ierror)
    end do
    disp = disp + 1

    !On place le debut apres le saut de ligne trouve (ci-dessus)
    disp = disp +1
    call MPI_FILE_SET_VIEW(fh,disp,type_end, &
         newtype,'native',MPI_INFO_NULL,ierror)

    !lecture des donnees (selon le type)
    if (type_mpi == MPI_REAL) then
       allocate(tmp_real4(1:amx_size(1)*amx_size(2)*amx_size(3)))
       !lecture des donnees (real(kind=4))
       call MPI_FILE_READ_ALL(fh, tmp_real4, &
            subsizes(1)*subsizes(2)*subsizes(3), &
            type_end, MPI_STATUS_IGNORE, ierror)
       if(testEndian2==1)call SWAP_F4(tmp_real4)
       var=real(tmp_real4,8)
       deallocate(tmp_real4)
    elseif (type_mpi == MPI_DOUBLE_PRECISION) then
       !lecture des donnees (real(kind=8))
       call MPI_FILE_READ_ALL(fh, var, &
            subsizes(1)*subsizes(2)*subsizes(3), &
            type_end, MPI_STATUS_IGNORE, ierror)
       if(testEndian2==1)call SWAP_F8(var)
    end if

    !fermeture du fichier
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)

!    if(testEndian2==1) call MPI_TYPE_FREE(type_end,ierror)
!    if(allocated(blockLen)) deallocate(blockLen)
!    if(allocated(dplt)) deallocate(dplt)

    return
  end subroutine readreal_bin8 
!-------------------------------------------------------------------------------



!===============================================================================
!
!            SUBROUTINE WRITE_BIN32
!> Ecriture de la partie binaire d'un fichier vtk 
!>        Par defaut avec 2 lignes ASCII : LOOKUP_TABLE et nom de la composante
!>        OU sans ces lignes ascii
!!
!! ATTENTION : pour var, 
!!      on n'utilise volontairement pas le profil implicite pour reshaper la variable
!!      => prévoir de passer la taille du tableau pour verif
!!
!! \param[in] ipencil: orientation des pinceaux du tableau var
!! \param[in] var: tableau a ecrire (de type REAL32)
!! \param[in] nomVar: nom de la variable a mettre dans le vtk
!! \param[in] nomFic: nom du fichier
!! \param[in] ascii_comp: variable optionelle (.true. ou .false.) pour
!!            .true. si non present
!!            ajoute deux lignes ASCII (LOOKUP_TABLE et nom de la composante
!!                   pour les fichiers vtk)
!!
!!Les donnees du fichier sont ecrites sur 4 octets (reel en simple precision)
!!
!! 28/10/2019 : traitement endianness dans convert_field_for_write_bin32
!!
!-------------------------------------------------------------------------------
  subroutine write_bin32(ipencil,var,nomVar,nomFic,ascii_comp)
    
    implicit none
    
    integer, intent(in)                 :: ipencil
    integer                             :: Fvtk
    real(kind=REAL32),dimension(amx_start(1):amx_end(1),amx_start(2):amx_end(2),amx_start(3):amx_end(3)), intent(in) &
                                        :: var
    character(len=*),intent(in)         :: nomFic, nomVar
    logical , optional , intent(in)     :: ascii_comp
    
    logical                             :: ascii_comp_loc
#ifdef NEW_2DECOMP
    TYPE(DECOMP_INFO),pointer                        :: decomp, sp_decomp
#else
    TYPE(DECOMP_INFO)                                :: decomp, sp_decomp
#endif
    integer(kind=MPI_OFFSET_KIND)       :: filesize, disp
    integer, dimension(3)               :: sizes, subsizes, starts
    integer                             :: ierror, newtype, fh, ios, type_reel,rank
    real(kind=REAL64)                   :: t1
    
   !initialisation temps
   t1 = MPI_WTIME()

    !variable pour definir un type de donnees adapte
    !On ecrit des reels en simple precision (4 octets)
!    integer,dimension(4) ::blockLen
!    integer(kind=MPI_ADDRESS_KIND),dimension(4):: dplt

    
   if (      present(ascii_comp)) ascii_comp_loc=ascii_comp
   if (.not. present(ascii_comp)) ascii_comp_loc=.true.
    
   call MPI_Barrier(amitex_mpi_comm,ierror) 
   call MPI_COMM_RANK(amitex_mpi_comm, rank, ierror) !rang du processus actuel
   
   ! ecriture des lignes ASCII par le proc 0 
   if(rank==0) then
      open(newunit=Fvtk, file= nomFic,form="formatted", &
                           status="old", action="write", position="APPEND",iostat= ios)
        if ( ios /= 0 ) then
          write(OUTPUT_UNIT,"(3A,I4,A)") " Probleme a l'ouverture (fichier: ",nomFic,") - ios = ",ios," (write_bin32)"
          stop 
        end if
        if (ascii_comp_loc .eqv. .true.) then
          !write(Fvtk,"(/,/,3A)")"SCALARS ",nomVar," float"   !interet des /,/, ???
          write(Fvtk,"(3A)")"SCALARS ",nomVar," float"
          write(Fvtk,"(A)")"LOOKUP_TABLE default"
        end if      

      close(Fvtk)
   end if
 
   call MPI_Barrier(amitex_mpi_comm,ierror) 
#ifdef NEW_2DECOMP
    decomp    => decomp_2d_fft_get_ph()
    sp_decomp => decomp_2d_fft_get_sp()
#else
    call get_decomp_fft_info(decomp,sp_decomp)
#endif
    
!tailles, debut et etendues des sous tableau 
!starts ramene a 0 (pour MPI) 
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)
    
!taille et bornes des tableau selon l'orientation:
!pinceaux selon x, y ou z   
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
    

!creation du type utilise par mpi en fonction de l'endianness
!  if(testEndian2==1)then  
!    blockLen=1
!    do i=1,4
!      dplt(i)=4-i
!    end do 
!    call mpi_type_create_hindexed(4,blockLen,dplt,MPI_integer1,type_reel,ierror)
!    call MPI_Type_commit(type_reel,ierror)  
!  else
!    type_reel = MPI_REAL
!  end if
   type_reel = MPI_REAL

!creation de la structure de donnees
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, MPI_REAL, newtype, ierror)    
    call MPI_TYPE_COMMIT(newtype,ierror)


!ouverture du fichier    
    call MPI_FILE_OPEN(amitex_mpi_comm, nomFic, MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                    MPI_INFO_NULL, fh, ierror)


! deplacement ( permet de ne pas ecraser l'en-tete)   
    call MPI_File_get_size(fh,filesize,ierror)
    disp = filesize

!ecriture des donnees

   call MPI_FILE_SET_VIEW(fh,disp,type_reel, &
         newtype,'native',MPI_INFO_NULL,ierror)
   call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         type_reel, MPI_STATUS_IGNORE, ierror)

!fermeture du fichier
    call MPI_FILE_CLOSE(fh,ierror)
!    if(testEndian2==1)call MPI_TYPE_FREE(type_reel,ierror)
    call MPI_TYPE_FREE(newtype,ierror)

!ecriture temps
    call MPI_Barrier(amitex_mpi_comm, ierror)
    if(nrank==0) then
    write(OUTPUT_UNIT,fmt="(A,E15.8)") "write_bin32 time on proc 0 (s) : ",MPI_WTIME()-t1
    end if

  end subroutine  write_bin32

!===============================================================================
!
!    SUBROUTINE CONVERT_FIELDxx_FOR_WRITE_BIN32
!
!>        CONVERT    (64bits -> 32bits) if necessary 
!>    and SWAP BYTES (from little to big endian if necessary 
!!
!! \param[in] ipencil: orientation des pinceaux du tableau var
!! \param[out] var: tableau a ecrire (de type REAL32)
!!
!!
!!
!-------------------------------------------------------------------------------
  subroutine convert_field32_for_write_bin32(tmp32big, Var)

    implicit none  
  
    real(kind=REAL32), dimension(:),intent(in)  :: Var
    real(kind=REAL32), dimension(:),intent(out) :: tmp32big

    !variables pour le traitement de l'endianness
    integer(kind=1), dimension(2)::testEndian
    integer(kind=2)              ::testEndian2
    equivalence(testEndian,testEndian2) ! testEndian2 partage la même zone mémoire que testEndian
    
    !> Conversion 32 bits
    tmp32big = Var

    !> Swap big endian, if necessary
    testEndian=(/int(1,kind=1),int(0,kind=1)/)
    if(testEndian2==1)call SWAP_F4(tmp32big) !traitement endianness

  end subroutine convert_field32_for_write_bin32

!-------------------------------------------------------------------------------
  subroutine convert_field64_for_write_bin32(tmp32big, Var)

    implicit none  
  
    real(kind=REAL64), dimension(:),intent(in)  :: Var
    real(kind=REAL32), dimension(:),intent(out) :: tmp32big

    !variables pour le traitement de l'endianness
    integer(kind=1), dimension(2)  ::testEndian
    integer(kind=2)                ::testEndian2
    equivalence(testEndian,testEndian2) ! testEndian2 partage la même zone mémoire que testEndian
    
    !> Conversion 32 bits
    tmp32big = real(Var,kind=REAL32)

    !> Swap big endian, if necessary
    testEndian=(/int(1,kind=1),int(0,kind=1)/)
    if(testEndian2==1)call SWAP_F4(tmp32big) !traitement endianness

  end subroutine convert_field64_for_write_bin32

!===============================================================================
!
!Fonctions d'initialisation du materiau avec fichier vtk
!
!===============================================================================


!===============================================================================
!
!  SUBROUTINE INITVTK1 : Initialisation du nombre de materiaux et des dimensions de la cellule enregistrement de l'en-tete
!
! \param[in] ficVTK: nom du fichier d'entree (format vtk)
! \param[in]  FU: unite logique du fichier ficVTK
!
! \param[out] entete: en-tete de 10 lignes du fichier vtk
! \param[out]        nmateriau: nombre de materiaux dans la cellule
! \param[out]        nx,ny,nz: dimensions de la cellule
! \param[out]        dx,dy,dz: spacing (present dans l'entete vtk)
! \param[out]        type_mpi: type des donnees du fichier vtk (entier sur 2 ou 4 octets :
!                        MPI_INTEGER2 ou MPI_INTEGER4)
!-------------------------------------------------------------------------------
!subroutine initVTK1(ficVTK, FU, entete,nx,ny,nz,dx,dy,dz,type_mpi)
! 
!    use MPI
!
!    implicit none
!    
!    character(len=*), intent(in)                :: ficVTK
!    integer,intent(in)                          :: FU
!    character(len=*), dimension(11), intent(out):: entete
!    integer, intent(out)                        :: nx,ny,nz,type_mpi
!    real(mytype), intent(out)                           :: dx,dy,dz
!    
!    integer                                     :: rank, ierror
!
!    call MPI_Comm_rank(amitex_mpi_comm, rank, ierror)
!
!    if(rank==0) call read_header_vtk(trim(ficVTK), FU, entete, nx, ny, nz, dx, dy, dz, type_mpi)
!
!    call MPI_Bcast(nx,1,MPI_INTEGER,0,amitex_mpi_comm,ierror)
!    call MPI_Bcast(ny,1,MPI_INTEGER,0,amitex_mpi_comm,ierror)
!    call MPI_Bcast(nz,1,MPI_INTEGER,0,amitex_mpi_comm,ierror)
!    call MPI_Bcast(type_mpi,1,MPI_INTEGER,0,amitex_mpi_comm,ierror)
!    call MPI_Barrier(amitex_mpi_comm,ierror)
!
!
!  end subroutine initVTK1
!-------------------------------------------------------------------------------

  
  
  
!===============================================================================
!
!               FONCTION INV_BORDER_INT4
!
!inverse l'ordre des octets d'un entier code sur 4 octets(kind=4)
!   conversion: big endian <-> little endian
!
!input:
!   tab     tableau d'entier dont on veut inverser les octets
!   starts  tableau de 3 entiers: indices de debut du tableau (trois dimensions)
!   end     tableau de 3 entiers: indices de fin du tableau (trois dimensions)
!output:
!   tab     tableau d'entier dont les octets ont ete inverses
!
!===============================================================================
!    subroutine inv_bOrdre_int4(tab,starts,ends)
!        
!        integer,dimension(3), intent(in)                :: starts, ends
!        integer(kind=4),dimension(starts(1):ends(1),starts(2):ends(2),starts(3):ends(3)), intent(inout) :: tab
!        integer    (kind=4) :: tmp
!        real i,j,k
!        
!        
!        print *, "bytes swap"
!        
!        do i=starts(1),ends(1)
!            do j=starts(2),ends(2)
!                do k=starts(3),ends(3)
!                    call MVBITS(tab(i,j,k),0,8,tmp,24)
!                    call MVBITS(tab(i,j,k),8,8,tmp,16)
!                    call MVBITS(tab(i,j,k),16,8,tmp,8)
!                    call MVBITS(tab(i,j,k),24,8,tmp,0)
!                    tab(i,j,k)=tmp
!                end do
!            end do
!        end do
!      
!    end subroutine inv_bOrdre_int4
!-------------------------------------------------------------------------------




!===============================================================================
!===============================================================================
!fonction de lecture (parsing) 
!
!recherche d'un entier ou d'un reel
!
!===============================================================================
!===============================================================================


!===============================================================================
!
!                  SUBROUTINE NEXTINT
!> Lecture du premier entier dans une chaine de caracteres
!!
!! \param[in] ligne: chaine de caracteres ou l'on recherche un entier
!! \param[inout]     i: indice du caractere a partir duquel on recherche l'entier
!!                   i<=1 : debut de chaine
!! \param[out]    n: premier entier trouve
!! \param[out]    i: indice du caractere apres le premier entier trouve
!!
!! s'il n'y a pas d'entier dans la chaine, i vaut -1
!-------------------------------------------------------------------------------
subroutine nextInt(ligne,i,n)

  implicit none

  character(len=*),intent(in)    :: ligne
  integer, intent(inout)         :: i
  integer, intent(out)           :: n
  integer                        :: tmp,i_init,k

  n=0
  k=0
  i_init = i
  !! On cherche ou est le prochain entier
  !
  !ci-dessous : cree une erreur avec -g -fcheck=all
  ! do while(index('0123456789', ligne(i:i) ) == 0 .AND. i < len(ligne)+1)
  !    i=i+1
  ! end do
  do while (i < len(ligne)+1)
     if (index('0123456789', ligne(i:i) ) /= 0) exit
     i=i+1
  end do

  if(i==len(ligne)+1) then
    i=-1
  else
     !! Tant que les caracteres sont des entiers on continue de lire
     do while(index('0123456789', ligne(i:i)  ) /= 0)
        read(ligne(i:i),*)tmp
        k=k+1
        i=i+1
        !n = 10*n+tmp
        if (k<9) then    ! entier sur 32 bits = 2^(32-1)~1e9 => choix 8 chiffres max
           n = 10*n+tmp
        else
           call amitex_abort("ERREUR (nextInt):Lecture d'un entier sur plus de 8 caracteres", 2, 0) 
        end if
     end do
  end if
end subroutine nextInt
!-------------------------------------------------------------------------------


subroutine nextInt8(ligne,i,n)

  implicit none

  character(len=*),intent(in)    :: ligne
  integer, intent(inout)         :: i
  integer(kind=8), intent(out)   :: n
  integer                        :: tmp,k

  n=0
  k=0
  !! On cherche ou est le prochain entier
  !
  !ci-dessous : cree une erreur avec -g -fcheck=all
  !do while(index('0123456789', ligne(i:i) ) == 0 .AND. i<len(ligne)+1)
  !  i=i+1
  !end do
  do while (i < len(ligne)+1)
     if (index('0123456789', ligne(i:i) ) /= 0) exit
     i=i+1
  end do

  if(i==len(ligne)+1) then
    i=-1
  else
    do while(index('0123456789', ligne(i:i)  ) /= 0 .AND. i<len(ligne)+1)
      read(ligne(i:i),*)tmp
      k=k+1
      i=i+1
      !n = 10*n+tmp
      if (k<18) then    ! entier sur 64 bits = 2^(64-1)~1e18 => choix 17 chiffres max
         n = 10*n+tmp
      else
         call amitex_abort("ERREUR (nextInt8):Lecture d'un entier, &
              & ou d'une mantisse ou partie decimale sur plus de 17 caracteres", 2, 0)
      end if
    end do
  end if

end subroutine nextInt8


!===============================================================================
!
!     SUBROUTINE NEXTREAL : 
!> Lecture du premier reel dans une chaine de caracteres
!!
!!
!!
!!
!! \param[in] ligne: chaine de caracteres ou on recherche un reel
!! \param[in]  i: indice du caractere a partir duquel on recherche un reel
!!
!! \param[out]  n: premier reel trouve
!! \param[out]  i: indice du caractere apres le premier reel trouve
!!
!! s'il n'y a pas de reel dans la chaine, n vaut 0
!!
!!
!!------------------------------------------------------------------------------
!! CORRECTIONS : 
!!     21/01/2016 : lecture des réels "1. " ou "1.e2" (il fallait ("1.0" "1.0e1") 
!!     16/01/2018 : lecture des reels dont la partie entiere ou decimale
!!         est un entier superieur a la limite des 32 bits : nextint -> nextint8
!-------------------------------------------------------------------------------
subroutine nextReal(ligne,i,n)

  implicit none

  character(len=*), intent(in)  :: ligne
  integer, intent(inout)        :: i
  real(mytype), intent(out)     :: n
  integer(INT64)                :: tmp
  integer                       :: deci
  real(mytype)                  :: r

  n=0.
  deci = 1



    tmp=0
    call nextint8(ligne,i,tmp)
    n=real(tmp)
    if (index('.,',ligne(i:i))/=0 ) then
      if (index(' eE',ligne(i+1:i+1))==0 ) then !lecture "1. " ou "1.e2"..
         deci=i
         tmp=0
         call nextInt8(ligne,i,tmp)
         i=i-1
         r=10._mytype**(real(deci-i,mytype))
         n = n+real(tmp)*r
      end if
    end if
    i=i+1
    if (index('eE',ligne(i:i))/=0 ) then
      i=i+1
      if (index('-',ligne(i:i))/=0 )  then
        deci=-1
      else
        deci =1
        i=i-1
      end if
      tmp=0
      call nextInt8(ligne,i,tmp)
      n=n*(10.**dble(tmp*deci))
    end if
end subroutine nextReal
!---------------------------------------------------------------

!===============================================================================
!
!       SUBROUTINE GET_REAL_VECT : 
!> Lecture d'un vecteur de reels dans une chaine de caracteres
!!
!! \param[in]  line: (chaine de caracteres) ligne a parser
!! \param[in]       n: (entier) nombre maximal de reels lus
!! \param[out]    array: (entier) tableau de reels lus a partir de 'line'
!! \param[out]       n: (entier) nombre de reels lus
!! \param[out]       error: (entier) 1 s'il y a trop de reels 0 sinon
!!
!!
!! 30/06/21 : LG - remplace letcure par nextReal par "read(chaine,*,IOSTAT=ios) data"
!!                 
!-------------------------------------------------------------------------------
subroutine get_real_vect( line, array, n, error)

  implicit none

  integer, intent(inout)               :: n
  real(mytype),dimension(n),intent(out):: array
  character(len=*), intent(in)         :: line
  integer, intent(out)                 :: error

  integer                              :: i,j,ios
  real(mytype),dimension(n+1)          :: array2

  
  !i=1
  !j=0
  !error = 0
  !val = 0
  !do while(i>-1 .AND. j < n)
  !  j=j+1
  !  call nextReal(line,i,val)
  !  if(i>-1)then
  !    if(j > n)then
  !      error=1
  !      return
  !    else
  !      array(j)=val
  !    end if
  !  end if
  !end do

  !n=j-1
  
  !-------------IMPLEMENTATION with read instead of nextreal
  error = 0
  array2 = huge(array)
  read(line,*,IOSTAT=ios) array2
  
  if(ios>0) call amitex_abort("problem reading a list of real in xml file, probably loading.xml",2)
  
  j=n+1
  do i = n+1,1,-1
     if (array2(n)==huge(array)) j=i
  end do  
  
  if (j == n + 1) error = 1
  
  n = j - 1
  array = array2(1:n)
  
end subroutine get_real_vect

!===============================================================================
!
!       SUBROUTINE GET_INT_VECT : 
!> Lecture d'un vecteur d'entiers dans une chaine de caracteres
!!
!! \param[in]  line: (chaine de caracteres) ligne a parser
!! \param[in]       n: (entier) nombre maximal d'entiers lus
!! \param[out]    array: (entier) tableau d'entiers lus a partir de 'line'
!! \param[out]       n: (entier) nombre d'entiers lus
!! \param[out]       error: (entier) 1 s'il y a trop d'entiers 0 sinon
!-------------------------------------------------------------------------------
subroutine get_int_vect( line, array, n, error)

  implicit none

  integer, intent(inout)             :: n
  integer, dimension(n),intent(out)  :: array
  character(len=*), intent(in)     :: line
  integer, intent(out)               :: error

  integer          :: i,j, val
  
  i=1
  j=0
  error = 0
  val = 0
  
  do while(i>-1)
    j=j+1
    call nextInt(line,i,val)
    if(i>-1)then
      if(j>n)then
        error=1
        return
      else
        array(j)=val
      end if
    end if
  end do

  n=j-1

end subroutine get_int_vect





!==============================================================================
!
!    SUBROUTINE GET_VECT_XML : 
!
!> Permet d'initialiser les tableaux "coeff" et "VarInt" de mattot a partir de fichier(s) vtk
!!
!! \param[in]   nMat: (entier) numero du materiau
!! \param[in]   node_list: liste des noeuds contenant les valeurs a extraire
!! \param[in]   m,n: (entier) dimension du tableau extrait
!! \param[in]   nZone: (entier) nombre de zone du materiau nMat 
!! \param[in]   nZoneP: (entier) nombre de zone du materiau nMat pour le processus
!! \param[in]   ind_zone: (entier) tableau MattotP0%zone du materiau
!! \param[in]   vect_type: (chaine de caracteres) type de donnees (Coeff ou IntVar)
!! \param[in]   pos: (entier) vecteur MattotP0%pos du materiau
!! \param[in]   test_composite: (logical), optional, .true. if using composite voxels
!!
!! \param[out]   array: (reel) tableau des coefficients/ variables internes
!!
!!
!!
!! ATTENTION: \n
!!       Les composantes du tableau non renseignees par le fichier xml
!!       ont pour valeur NaN
!!
!! ATTENTION : \n
!!       Cette fonction doit etre appelee et doit passer par les fonctions de lecture de fichier
!!       binaire meme si le materiau n'est pas present sur le pinceau (nZoneP=0)
!!       -> cette lecture parallele doit etre faite par TOUS les proc. (meme ceux qui ne lisent rien!) 
!!
!! CHANGEMENT 04/04/2019 : corrections de la fonction pour les cas ou nZoneP=0 et array vide
!!                         voir aussi correction de getBinCoeffFromFile
!!
!==============================================================================
subroutine get_vect_xml( nMat, node_list, m,n, array, nZone, nZoneP, ind_zone, vect_type,pos,test_composite )

  implicit none
  
  integer,intent(in)                            :: nMat      ! material 
  type(nodeList),pointer, intent(inout)         :: node_list ! list of node containing array's values
  integer, intent(in)                           :: m         ! m,n :array size,  m: number of coeff / int Var
  integer(kind=8), intent(in)                   :: n         ! n: number of zone/voxels in current process
  integer(kind=8),intent(in)                    :: nZone,nZoneP     ! nZone number of zone and zone in process
  real(mytype),dimension(m,n), intent(out)      :: array     ! array to extract from xml file
  integer(kind=8),dimension(:), intent(in)      :: pos       ! position vector for the material concerned MattotP0%pos
  integer(kind=8),dimension(nZoneP,2),intent(in):: ind_zone  ! array zone from MattotP0 structure
  character(len=*), intent(in)                  :: vect_type ! type of the array to extract (Coeff or IntVar)
  logical,intent(in),optional                   :: test_composite ! test if composite voxels are used
  
  
  type(node),pointer        :: cur_node,attribut
  type(nodeList),pointer    :: node_subList
  integer                   :: i, j, ind, numZ
  integer                   :: alloc_stat
  integer(kind=8)           :: k, l, numZ_found
  character(len=50)         :: coeff_type, fmt_file
  real(mytype)              :: value0
  real(mytype),dimension(:), allocatable :: tmp_array
  character(len=200)        :: file_name, err_msg
  logical                   :: ok, ok_i, ok_j
  logical                   :: empty_rank
  real(mytype), allocatable, dimension(:) :: TEMP64      ! For reading 'VTK' file
  integer                   :: type_mpi,nx,ny,nz,ierror  !
  real(mytype)              :: dx,dy,dz,x0,y0,z0         !
  integer(kind=8)           :: n0_loc,n0
  logical                   :: test_composite1 = .false.

  !! Detection d'une utilisation de voxels composites
  if (present(test_composite)) test_composite1 = test_composite
  

  !! Determination prealable de la presence ou non de donnees a lire sur le pinceau
  !  leur absence est signalee par des tableaux vides en entree et nZoneP = 0
  empty_rank = .false.
  if (nZoneP == 0) empty_rank = .true.

  ok=.TRUE.
  value0=-1._mytype
  !array=sqrt(value)
  array=IEEE_VALUE(value0,IEEE_QUIET_NAN)
  !array = huge(value0) ! compatible with gcc -ffast-math

  !chaque noeud de node_list contient un coefficient (une variable interne) du
  !materiau nMat
  do i=0,getLength(node_list)-1
     ! noeud: coeff / IntVar
     cur_node => item(node_list,i)
     ind = -1
     coeff_type=""

     !noeud i valide
     ok_i=.TRUE.
     !indice du coefficient (variable interne)
     call get_int_xml(cur_node,"Index",ind,1)
     if (.not. empty_rank) then
         if(ind > getLength(node_list))then
            ok_i=.FALSE.
            write(err_msg,fmt="(A)") " indice superieur au nombre d'elements (get_vect_xml)"
            call amitex_abort(err_msg, 0,0)
         ! TODO : polyxCC_def_imp fails with nvfortran if uncommented (comment with nvfortran or find solution!!)
         else if( array(ind,1) == array(ind,1))then !remplacer par test IEEE_is_nan
         !else if( ieee_is_nan(array(ind,1)) )then ! pas mieux pour nvfortran 
         !else if( array(ind,1) /= huge(value0))then ! compatible with gcc -ffast-math option (NaN incompatible)
            write(err_msg,fmt="(A)") " indice present deux fois (get_vect_xml)"
            call amitex_abort(err_msg, 1,0)
            ok_i=.FALSE.
         end if
     end if
     call get_str_xml(cur_node, "Type", coeff_type,1)
     if(ind==-1 .or. coeff_type=="")then
        ok_i = .FALSE.
     end if

     if(.not. ok_i)then
        ! erreur, passer au suivant

     else if(trim(coeff_type)=="constant") then
        value0=0.2384316844789e-21 ! valeur quelconque
        call get_real_xml(cur_node,"Value",value0,1)
        if(value0==0.2384316844789e-21)then
           ok_i=.FALSE. !attribut value absent ...
        else
           if (.not. empty_rank)  array(ind,:) = value0  ! <=> MattotP0(j)%coeff(ind,:)=value
        end if

     else if(trim(coeff_type) == "constant_zone")then
        !! Ce coefficient/variable interne est constant par zone
        attribut => getAttributeNode(cur_node,"File")
        if(associated(attribut))then
           !! Si les valeurs par zone sont donnees dans un fichier, 
           !! on recupere le nom du fichier
           file_name = getAttribute(cur_node, "File")
           allocate(tmp_array(nZoneP))
           value0=-1._mytype
           tmp_array=IEEE_VALUE(value0,IEEE_QUIET_NAN)
           !tmp_array=sqrt(value)  !astuce pour fabriquer un NaN
           !! On cherche sous quel format sont les donnees
           attribut => getAttributeNode(cur_node,"Format")
           if(associated(attribut))then
              fmt_file = getAttribute(cur_node,"Format")
              call stringToLower(fmt_file)
              if(trim(fmt_file) == "binary") then 
                 !! Si le format est binaire, on lit le fichier binaire 
                 !! (l'unite du fichier pour la lecture est le numero d'indice + 20)
                 call getBinCoeffFromFile(file_name, ind_zone(:,2), nZoneP, tmp_array, ind+20)
              elseif(trim(fmt_file) == "ascii") then
                 !! Si le format est ASCII, on lit le fichier avec la fonction getCoeffFromFile
                 !! (l'unite du fichier pour la lecture est le numero d'indice + 20)
                 call getCoeffFromFile(file_name, ind_zone(:,2), nZoneP, tmp_array, ind+20)
              else
                 !! Si c'est un autre format on affiche une erreur
                 write(err_msg,fmt="(3A)") "Format de fichier inconnu : ",&
                      fmt_file," (get_vect_xml)"
                 call amitex_abort(err_msg,2,0)
              end if
           else 
              !! Si le format n'est pas precise, on suppose que c'est de l'ASCII
              call getCoeffFromFile(file_name, ind_zone(:,2), nZoneP, tmp_array, ind+20)
           end if
           !! On met les valeurs dans le tableau global regroupant toutes les zones
           if(vect_type== "Coeff" .OR. vect_type== "CoeffK" .OR. vect_type== "Coeff_composite" )then
              if (.not. empty_rank) array(ind,:)=tmp_array(:)
           else if(vect_type=="IntVar")then
              l=1
              if (.not. empty_rank) then
                  do k=1,nZoneP
                     array(ind,l:ind_zone(k,1))=tmp_array(k)
                     l=ind_zone(k,1)+1
                  end do
              end if
           end if
           if(allocated(tmp_array)) deallocate(tmp_array)
        elseif(.not. empty_rank) then !On shunt cette partie si materiau absent (empty_rank=0)
           !! Si le coefficient/variable interne est constant par zone et qu'aucun
           !! fichier n'est donne en entree, les donnees doivent se trouver dans le 
           !! fichier XML dans des noeuds Zone.
           call getNodeList(cur_node, node_subList, "Zone", -1,1)
           !! Les noeuds de node_sublist contiennent le numero de chaque zone
           !! et la valeur du coefficient (de la variable interne) associe
           !! dont l'indice est ind
           if(.not. associated(node_subList)) then
              !! Si il n'y a pas de noeud Zone => erreur
              ok_i=.FALSE.
           else
              !! Sinon, on parcourt ces noeuds pour recuperer les valeurs par zone
              k=1   ! k sera l'indice de zone courant
              do j=0,getLength(node_subList)-1
                 ok_j=.TRUE.
                 numZ=-1
                 !! On recupere le numero de zone (numZ) du noeud "Zone" courant
                 call get_int_xml(item(node_subList,j),"numZ",numZ,1)
                 if(numZ .LE. 0 .OR. numZ>nZone) then
                    !! Si le numero de zone est negatif ou superieur au nombre de zones 
                    !! dans le materiau => warning
                    ok_j = .FALSE.
                    if(numZ == -1) then
                       write(err_msg,fmt="(A,I0,A)") "numZ =-1 ou introuvable sur le noeud 'Zone' ",&
                            j+1," (get_vect_xml)"
                    elseif(numZ>nZone)then
                       write(err_msg,fmt="(A,I0,A)") "numero de zone invalide, noeud 'Zone' ",&
                            j+1," (get_vect_xml)"
                    else
                       write(err_msg,fmt="(2(A,I0),A)") "numZ <=0 (",&
                            numZ,") sur le noeud 'Zone' ",j+1," (get_vect_xml)"
                    end if
                    call amitex_abort(err_msg,-1,0)
                 else 
                    !! On cherche le numero local correspondant a la zone numZ
                    !! En partant de l'indice de zone local k correspondant a la zone
                    !! suivant celle correspondante au noeud "Zone" precedent
                    numZ_found=k
                    do while(numZ_found>0) 
                       !! recherche de numZ dans ind_zone
                       !! on garde la variable numZ_found qui est notre point de depart
                       if(numZ==ind_zone(k,2))then
                          numZ_found=0
                          !! si la zone numZ correspond a la zone k de ce processus
                          !! on recupere la valeur
                          value0=-1
                          call get_real_xml(item(node_subList,j), "Value",value0,1)
                          if(value0==-1) then
                             ok_j = .FALSE.
                             write(err_msg, fmt="(A,I0,A)") "noeud 'Zone' ",j+1," (get_vect_xml)"
                             call amitex_abort(err_msg,-1,0)
                          end if
                          !! On determine si c'est un coefficient ou une variable interne
                          if(vect_type=="Coeff" .OR. vect_type=="CoeffK" .OR. vect_type=="Coeff_composite")then
                             if (.not. empty_rank)  array(ind,k)=value0 !    MattotP0(j)%coeff( numZ , ind )=value
                          else if(vect_type=="IntVar")then
                             !! l est l'indice local du premier element de la zone k
                             !! dans le processus
                             if (.not. empty_rank) then 
                                 if(k==1)then
                                    l=1         ! position du premier element de la zone 1
                                 else
                                    l=ind_zone(k-1,1)+1 ! position du premier element de la zone k
                                 end if
                                 
                                 array(ind,l:ind_zone(k,1))=value0 
                                 !! ind_zone(k,1) est l'indice local du dernier element de la zone k
                                 !! tous les elements de la zone k sont initialises a value
                             end if
                          else
                             ok_j=.FALSE.
                             write(err_msg,fmt="(3A)") " type attendu: Coeff ou CoeffK ou IntVar, type present: '",&
                                  vect_Type,"' (get_vect_xml)"
                             call amitex_abort(err_msg,1)
                          end if
                       end if
                       k=k+1
                       if(k>nZoneP)then
                          !! Si on a depasse l'indice local max on revient on debut
                          k=1
                       else if(ind_zone(k,2)==0)then
                          !! Sinon, si le numero de zone global est 0 on passe au k suivant 
                          k=k+1
                       end if
                       !! Si k=numZ_found cela signifie que l'on a fait le tour des zones 
                       !! locales et que la zone numZ n'a pas ete trouvee lors
                       !! du parcours de ind_zone => cette zone n'est pas dans ce processus
                       if(k==numZ_found) then
                          numZ_found=0
                       end if
                    end do

                    if(numZ_found==-1)then
                       ok_j=.FALSE.
                    end if

                 end if
                 ok_i = ok_j .and. ok_i
              end do ! fin boucle sur les noeuds 'Zone'
           end if ! fin definition des coeff/ varInt dans les noeuds 'Zone'
        end if ! fin definition des coeff/ varInt

     !------------------------------- VARIABLES INTERNES INITIALES "VARIABLES" : LECTURE FICHIER VTK
     else if(trim(coeff_type) == "variable")then

        !! variable interne introduite sous forme de champ (vtk)
        attribut => getAttributeNode(cur_node,"File")
        if(associated(attribut))then
           !! on recupere le nom du fichier
           file_name = getAttribute(cur_node, "File")
           !! On cherche sous quel format sont les donnees
           attribut => getAttributeNode(cur_node,"Format")
           if(associated(attribut))then
              fmt_file = getAttribute(cur_node,"Format")
              call stringToLower(fmt_file)
              if(trim(fmt_file) == "vtk") then 
                 !! Si le format est vtk, on lit le fichier dans un tableau temporaire
                 !! puis on affecte les variables internes a partir des positions
                 if (.not. allocated(TEMP64)) then 
                    allocate(TEMP64(amx_size(1)*amx_size(2)*amx_size(3)),stat=alloc_stat)
                    if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (get_vect_xml, TEMP64)",2)
                 end if
                 call read_header_vtk(file_name, nx,ny,nz,dx,dy,dz,x0,y0,z0, type_mpi)
                 !Check consistency : nx *ny * nz = somme_procs(n0)!!!
                 n0_loc = amx_size(1)*amx_size(2)*amx_size(3)
                 call MPI_Allreduce(n0_loc,n0,1,MPI_INTEGER8,MPI_SUM,amitex_mpi_comm,ierror)
                 if (n0 /= nx*ny*nz) then
                    write(err_msg,fmt="(A)") "Size of the vtk file given for internal variable in 'material'.xml&
                                          & is not consistent with the grid size"
                    call amitex_abort(err_msg,2,0)
                 end if 
                 ! check mpi_type
                 if(type_mpi /= MPI_REAL .AND. type_mpi /= MPI_DOUBLE_PRECISION) then
                    call amitex_abort("Data-type in vtk file given in 'material'.xml for 'Variable' IntVar"// &
                                     achar(10)//"      => must be 'float' or 'double'",2,0)
                 end if
                 call read_bin(file_name,1,type_mpi,TEMP64)
                 if (.not. empty_rank) array(ind,:)=TEMP64(pos)
                 !! WARNING message when using composite voxels
                 if (test_composite1) call amitex_abort(achar(10)//achar(10)//&
         "WARNING - WARNING - WARNING "// achar(10)//&
         "=========================== "// achar(10)//&
         "   Using Composite Voxels with Type=Variable for initializing internal variables is DANGEROUS."// achar(10)//&
         "   Up to now, a phase in a composite voxel is defined by (numM, numZ),"// achar(10)//&
         "   its internal variable will correspond to the last voxel of the zone numZ in material numM."//achar(10)//&
         achar(10)//&
         "   => PREFER USING ONE ZONE PER VOXEL IN THAT CASE."//achar(10)//&
         achar(10),0,0) ! warning (0) on node0 (0)

              else
                 !! Si c'est un autre format on affiche une erreur
                 write(err_msg,fmt="(3A)") "Unknown 'Format' for 'Variable' internal variable &
                      & in 'material'.xml : ",fmt_file,". Format available : vtk (get_vect_xml)"
                 call amitex_abort(err_msg,2,0)
              end if
           else
              write(err_msg,fmt="(A)") "'Format' not precised for 'Variable' internal variable &
                                       &in 'material'.xml (get_vect_xml)"
              call amitex_abort(err_msg,2,0)
           end if
         else
           write(err_msg,fmt="(A)") "'File' not precised for 'Variable' internal variable&
                                    & in 'material'.xml (get_vect_xml)"
           call amitex_abort(err_msg,2,0)
         end if

     else
        ok_i=.FALSE.
        write(err_msg,fmt="(3A)") " type de coefficient non reconnu: ",coeff_type," (get_vect_xml)"
        call amitex_abort(err_msg, 1,0)
     end if

     if(.not. ok_i)then
        ok=.FALSE.
        write(err_msg, fmt="(A,I0, 3A,I0,A)") "materiau numM=",nMat,", noeud '",&
             getNodeName(cur_node),"' ",i+1, "(get_vect_xml)"
        call amitex_abort(err_msg,-1,0)
     end if

  end do

end subroutine get_vect_xml
!==============================================================================


!==============================================================================
!
!                   SUBROUTINE GETBINCOEFFFROMFILE : 
!
!> Lecture dans un fichier binaire des valeurs d'un coefficient (variables internes) 
!> pour un materiau selon la zone (cas constant par zone) 
!!
!! \param[in]   file_name: (chaine de caracteres) nom du fichier a lire
!! \param[in]   ind_zone: (entier) tableau des indices des zones presentes dans le materiau
!!                        et dans ce processus
!! \param[in]   nZoneP: (entier) nombre de zones presentes dans le materiau dans ce processus
!! \param[in]   FU: (entier) unite du fichier ou sont lus les coefficients
!! \param[out]  array: (reel) tableau des coefficients/variables internes locaux
!!
!!
!! FONCTIONNEMENT \n
!! On charge sur chaque proc les indices compris entre ind_zone_min et ind_zone_max
!!
!! ATTENTION \n
!! On peut charger au maximum 2^31-1 zones par pinceaux (ind_zone_max-ind_zone_min+1)
!!      => PB associe au nbre de donnees a transmettre dans les fonctions MPI : kind=4
!!      => cette erreur est detectee a la compilation pour certains couples "compilo/MPI"
!!
!!
!! INCONVENIENT POTENTIEL \n
!! Encombrement memoire! 
!!
!! CHANGEMENT 14/12/2018 : inversion ordre de byte APRES lecture MPI 
!!             => pb de lenteur excessive a lecture de donnees composites
!!
!! CHANGEMENT 04/04/2019 : correction d'une erreur si nZoneP=0 et array vide  
!!            getBinCoeffFromFile (appelee par get_vect_xml appelee par read_mat)
!!            pouvait ne pas etre appelee par tous les 
!!            pinceaux (i.e. ceux n'ayant pas de coefficients a lire dans le fichier binaire),
!!            ce qui occasionnait des bugs après la tentative de lecture
!!            du fichier de coefficients (processus coinces dans la procedure)
!!            Pour résoudre ce problème, des tableaux vides peuvent maintenant etre passes
!!            a la procedure pour qu'elle puisse etre appellee par les processus
!!            n'ayant pas de coefficients a lire depuis le fichier
!!       VOIR CORRECTIONS get_vect_xml et read_mat
!==============================================================================
subroutine getBinCoeffFromFile(file_name, ind_zone, nZoneP, array, FU )

  implicit none

  character(len=*),intent(in)                   :: file_name
  integer(kind=8),intent(in)                    :: nZoneP
  integer(kind=8),dimension(nZoneP),intent(in)  :: ind_zone
  real(mytype),dimension(nZoneP), intent(out)   :: array
  !! integer,intent(in)                            :: FU
  !> Ne compile pas en OpenMPI si intent(in). Raison inconnue
  integer                                       :: FU
  character(len=200), dimension(2)              :: entete
  integer(kind=8)                               :: i, nZoneMax, nMax,nCoeffFile,nSum
  integer(kind=8)                               :: ind_zone_max,ind_zone_min
  integer(kind=4)                               :: nCoeffProc 
         !kind=8 ne passe pas avec les fonctions MPI_IO pour certains couples "compilateur/MPI"
  real(mytype)                                  :: valueSum, tmpValue
  character(len=200)                            :: err_msg
  integer                                       :: ierror, ios, alloc_stat
  integer                                       :: local_offset
  logical                                       :: file_exists
  integer                                       :: type_mpi,nOctets
  integer(kind=MPI_OFFSET_KIND)                 :: disp, filesize,offset
  integer(kind=1)                               :: debut
    
  !variables du traitement de l'endianness
  integer(kind=1), dimension(2)::testEndian
  integer (kind=2)             ::testEndian2
  equivalence(testEndian,testEndian2)
  
  ! flag de traitement du cas on aucun coefficient n'est a lire sur le pinceau
  logical :: empty_rank

  !variables pour definir un type MPI de donnees adapte
  integer :: type_end
!  integer,dimension(:),allocatable ::blockLen
!  integer(kind=MPI_ADDRESS_KIND),dimension(:), allocatable :: dplt
  
  !> Variables temporaires qui seront utilisees pour la lecture
  integer(kind=INT8),dimension(:),allocatable      :: tmp_int1
  integer(kind=INT16),dimension(:),allocatable     :: tmp_int2
  integer(kind=INT32),dimension(:),allocatable     :: tmp_int4
  integer(kind=INT64),dimension(:),allocatable     :: tmp_int8
  real(kind=REAL32),dimension(:),allocatable       :: tmp_real
  real(kind=REAL64),dimension(:),allocatable       :: tmp_double
  
  ! avoid gcc-warning
  allocate(tmp_int1(1)); deallocate(tmp_int1)
  allocate(tmp_int2(1)); deallocate(tmp_int2)
  allocate(tmp_int4(1)); deallocate(tmp_int4)
  allocate(tmp_int8(1)); deallocate(tmp_int8)
  allocate(tmp_real(1)); deallocate(tmp_real)
  allocate(tmp_double(1)); deallocate(tmp_double)

  !! Determination prealable de la presence ou non de donnes a lire sur le pinceau
  !  leur absence est signalee par des tableaux vides en entree et nZoneP = 0
  empty_rank = .false.
  if (nZoneP == 0) empty_rank = .true.

  !!================================================= PRELIMINAIRES
  !!                                                 lecture entete
  !!                              definition du type mpi (type_end)
  !!                        recherche du debut des donnees binaires

  testEndian=(/int(1,kind=1),int(0,kind=1)/)

  !! On verifie que le fichier existe
  inquire(FILE=file_name,EXIST=file_exists)
  if(.not.file_exists)then 
    write(err_msg,fmt="(3A)") "Le fichier: '",trim(file_name),"' n'existe pas ,(getBinCoeffFromFile )"
    call amitex_abort(err_msg,1,0)
    return
  end if

  !! Lecture et stockage de l'en-tete (les 2 premieres lignes) dans le tableau entete
  open(unit=FU, file=trim(file_name),form="formatted", status="old", action="read",iostat= ios)
  if ( ios /= 0 ) then
     write(err_msg,fmt="(3A,I0,A)") "Probleme a l'ouverture du fichier : ",trim(file_name),&
          " , ",ios," (getBinCoeffFromFile)"
     call amitex_abort(err_msg,2,0) 
  end if
  do i=1,2
      read(FU,FMT='(A)') entete(i)
  end do
  close(unit=FU)

  !! Recuperation du nombre de zones dans le fichier binaire 
  !! local_offset= 1: on cherche un entier a partir du debut de la ligne 
  local_offset=1
  call nextInt8(entete(1),local_offset,nCoeffFile)
  
  !! Recuperation du format des donnees ecrites en binaire
  call type2mpi_type(trim(entete(2)),type_mpi)
  
  !! nOctets = nombre d'octets par entier en fonction du type d'entier
  call MPI_TYPE_SIZE(type_mpi,nOctets,ierror)

  !! Creation du type utilise par mpi en fonction de l'endianness : CHANGEMENT
  !!                                     14/12/2018 -> on inverse les bytes APRES la lecture
!  if(testEndian2==1)then
!     allocate(blockLen(nOctets),dplt(nOctets))
!     blockLen=1
!     do i=1,nOctets
!        dplt(i)=nOctets-i
!     end do
!     if(nOctets>1)then 
!        !creation de la structure
!        call mpi_type_create_hindexed(int(nOctets,kind=4),blockLen,dplt,MPI_integer1,type_end,ierror)
!        call MPI_Type_commit(type_end,ierror)
!     else
!        type_end = type_mpi
!     end if
!  else
!     type_end = type_mpi
!  end if
!  call MPI_TYPE_COMMIT(type_end,ierror)
  type_end=type_mpi


  !! nZone est le numero de zone maximal de ce materiau dans ce processus 
  !! => nMax est le numero de zone maximal de ce materiau i.e.le nombre de zones dans le materiau
  if (empty_rank) then
     ! pas de donnes a recuperer pour le pinceau : on utilise les indices pour allouer
    ! des tableau vides allocate(0:0)
      ind_zone_max = 0
      ind_zone_min = 0
  else 
      ind_zone_max = maxval(ind_zone)
      ind_zone_min = minval(ind_zone)
  end if
  call MPI_Allreduce(ind_zone_max,nMax,1,MPI_INTEGER8,MPI_MAX,amitex_mpi_comm,ierror)

  !! Si le nombre de zones est trop eleve => erreur
  nZoneMax=int(2,kind=8)**31-1  
  if(nMax>nZoneMax)then
     write(err_msg,fmt="(A)") "Nombre de zones trop important &
          & (max=2^31 -1) (getBinCoeffFromFile)"
     call amitex_abort(err_msg,1,0)
     return
  end if

  !! On verifie que le nombre de zones dans le materiau est superieur a celui entre
  !! dans le fichier binaire
  if(nCoeffFile < nMax)then
     write(err_msg,fmt="(3A,I8,A,I8)") "Nombre de zones trop faible dans le fichier de coefficient "&
          ,trim(file_name)," par rapport au fichier VTK (getBinCoeffFromFile) : nCoeffFile="&
          ,nCoeffFile," nMax= ", nMax
     call amitex_abort(err_msg,2,0)
     return
  end if

  !! Ouverture du fichier en parallele
  call MPI_FILE_OPEN(amitex_mpi_comm, file_name, &
       MPI_MODE_RDONLY, MPI_INFO_NULL, &
       FU, ierror)
  call MPI_BARRIER(amitex_mpi_comm,ierror)
  
  !! On recupere la taille du fichier
  call MPI_File_get_size(FU,filesize,ierror)

  disp = nCoeffFile*int(nOctets,kind=8)
  !! nCoeffFile*nOctets doit etre la taille des donnees entrees
  disp = (filesize) - disp
  if(disp <1)then  
     !! Si le fichier est plus petit que le nombre de donnees qu'il doit contenir => erreur
     write(err_msg,fmt="(3A)") "Erreur lors de la lecture du fichier de coefficient ", &
          file_name," (getBinCoeffFromFile)"
     call amitex_abort(err_msg,2)
  end if
  call MPI_File_seek(FU,disp,MPI_SEEK_SET ,ierror)
  call MPI_File_read(FU, debut, 1, MPI_INTEGER1, MPI_STATUS_IGNORE, ierror)
  !! Recherche du premier saut de ligne avant le debut estime des donnees
  !! (10 en ASCII correspond a un saut de ligne)
  i=0
  do while( (debut/=10) )
     i=i+int(1,kind=1)
     disp = disp -1
     call MPI_File_seek(FU,disp,MPI_SEEK_SET,ierror)
     call MPI_File_read(FU,debut,1,MPI_INTEGER1,MPI_STATUS_IGNORE,ierror)
  end do

  !On place le debut apres le saut de ligne trouve ci-dessus
  disp = disp +1

  !!======================================= LECTURE PARALLELE DES DONNEES BINAIRES
  !!          chaque processeur lit les donnees entre ind_zone_min et ind_zone_max                        
  !!                                                  dans tmp_array(1:nCoeffProc)
  !!
  !! Pour chaque type : Allocation / Lecture / Affectation / Conversion real(mytype) 
  !!
  if (empty_rank) then
    ! cas d'absence de donnes a lire sur le pinceau
    ! 0 coefficients a lire
      nCoeffProc = 0
      offset = disp  
  else
      nCoeffProc = int(ind_zone_max - ind_zone_min,kind=INT32) + 1
      offset = disp + (ind_zone_min - 1) * nOctets
  end if

  ! Allocation / Affectation / Conversion pour chaque type_mpi
  if(type_mpi == MPI_REAL) then
        allocate(tmp_real(ind_zone_min:ind_zone_max),stat=alloc_stat); tmp_real=0 
        call MPI_FILE_READ_AT(FU, offset, tmp_real(ind_zone_min:ind_zone_max), nCoeffProc, type_mpi, MPI_STATUS_IGNORE, ierror) !gcc-warning accepte
        if(testEndian2==1)call SWAP_F4(tmp_real)
        array = real(tmp_real(ind_zone),mytype) 

  elseif(type_mpi == MPI_DOUBLE_PRECISION) then
        allocate(tmp_double(ind_zone_min:ind_zone_max),stat=alloc_stat); tmp_double=0 
        call MPI_FILE_READ_AT(FU, offset, tmp_double(ind_zone_min:ind_zone_max), nCoeffProc, type_mpi, MPI_STATUS_IGNORE, ierror) 
        if(testEndian2==1)call SWAP_F8(tmp_double)
        array = real(tmp_double(ind_zone),mytype) 

  elseif(type_mpi == MPI_INTEGER1) then
        allocate(tmp_int1(ind_zone_min:ind_zone_max),stat=alloc_stat); tmp_int1=0 
        call MPI_FILE_READ_AT(FU, offset, tmp_int1(ind_zone_min:ind_zone_max), nCoeffProc, type_mpi, MPI_STATUS_IGNORE, ierror) 
        array = real(tmp_int1(ind_zone),mytype) 

  elseif(type_mpi == MPI_INTEGER2) then
        allocate(tmp_int2(ind_zone_min:ind_zone_max),stat=alloc_stat);tmp_int2=0 
        call MPI_FILE_READ_AT(FU, offset, tmp_int2(ind_zone_min:ind_zone_max), nCoeffProc, type_mpi, MPI_STATUS_IGNORE, ierror) 
        if(testEndian2==1)call SWAP_I2(tmp_int2)
        array = real(tmp_int2(ind_zone),mytype) 

  elseif(type_mpi == MPI_INTEGER4) then
        allocate(tmp_int4(ind_zone_min:ind_zone_max),stat=alloc_stat); tmp_int4=0 
        call MPI_FILE_READ_AT(FU, offset, tmp_int4(ind_zone_min:ind_zone_max), nCoeffProc, type_mpi, MPI_STATUS_IGNORE, ierror) 
        if(testEndian2==1)call SWAP_I4(tmp_int4)
        array = real(tmp_int4(ind_zone),mytype) 

  elseif(type_mpi == MPI_INTEGER8) then
        allocate(tmp_int8(ind_zone_min:ind_zone_max),stat=alloc_stat);tmp_int8=0 
        call MPI_FILE_READ_AT(FU, offset, tmp_int8(ind_zone_min:ind_zone_max), nCoeffProc, type_mpi, MPI_STATUS_IGNORE, ierror) 
        if(testEndian2==1)call SWAP_I8(tmp_int8)
        array = real(tmp_int8(ind_zone),mytype) 

  else
     write(err_msg,fmt="(A)") "Format de nombre non reconnu"
     call amitex_abort(err_msg,1,0)
     return
  end if

  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_composite)",2)

!  ! Lecture / Affectation / Conversion
!  if(type_mpi == MPI_REAL) then
!       call MPI_FILE_READ_AT(FU, offset, tmp_real(ind_zone_min:ind_zone_max), nCoeffProc, type_mpi, MPI_STATUS_IGNORE, ierror) !gcc-warning accepte
!       if(testEndian2==1)call SWAP_F4(tmp_real)
!       array = real(tmp_real(ind_zone),mytype) !gcc-warning accepte
!  elseif(type_mpi == MPI_DOUBLE_PRECISION) then
!       call MPI_FILE_READ_AT(FU, offset, tmp_double(ind_zone_min:ind_zone_max), nCoeffProc, type_mpi, MPI_STATUS_IGNORE, ierror) !gcc-warning accepte
!       if(testEndian2==1)call SWAP_F8(tmp_double)
!       array = real(tmp_double(ind_zone),mytype) !gcc-warning accepte
!  elseif(type_mpi == MPI_INTEGER1) then
!       call MPI_FILE_READ_AT(FU, offset, tmp_int1(ind_zone_min:ind_zone_max), nCoeffProc, type_mpi, MPI_STATUS_IGNORE, ierror) !gcc-warning accepte
!       array = real(tmp_int1(ind_zone),mytype) !gcc-warning accepte
!  elseif(type_mpi == MPI_INTEGER2) then
!       call MPI_FILE_READ_AT(FU, offset, tmp_int2(ind_zone_min:ind_zone_max), nCoeffProc, type_mpi, MPI_STATUS_IGNORE, ierror) !gcc-warning accepte
!       if(testEndian2==1)call SWAP_I2(tmp_int2)
!       array = real(tmp_int2(ind_zone),mytype) !gcc-warning accepte
!  elseif(type_mpi == MPI_INTEGER4) then
!       call MPI_FILE_READ_AT(FU, offset, tmp_int4(ind_zone_min:ind_zone_max), nCoeffProc, type_mpi, MPI_STATUS_IGNORE, ierror) !gcc-warning accepte
!       if(testEndian2==1)call SWAP_I4(tmp_int4)
!       array = real(tmp_int4(ind_zone),mytype) !gcc-warning accepte
!  elseif(type_mpi == MPI_INTEGER8) then
!       call MPI_FILE_READ_AT(FU, offset, tmp_int8(ind_zone_min:ind_zone_max), nCoeffProc, type_mpi, MPI_STATUS_IGNORE, ierror) !gcc-warning accepte
!       if(testEndian2==1)call SWAP_I8(tmp_int8)
!       array = real(tmp_int8(ind_zone),mytype) !gcc-warning accepte
!  end if

  !!================================================ AFFICHAGE DE QUELQUES VALEURS
  !! Affichage des valeurs du coefficient dans la premiere et la derniere zone :
  !! On recupere la valeur dans la premiere (1) et la derniere zone (nMax) dans chaque processus
  !! ou on les trouve et on obtient valeur en divisant la somme des valeurs 
  !! trouvees dans chaque processus par le nombre de processus concernes

  tmpValue = 0
  i=0
  if(ind_zone_min==1) then
     tmpValue = array(minloc(ind_zone,1))
     i=1
  end if
  call MPI_Reduce(tmpValue,valueSum,1,real_type,MPI_SUM,0,amitex_mpi_comm,ierror)
  call MPI_Reduce(i,nSum,1,MPI_INTEGER8,MPI_SUM,0,amitex_mpi_comm,ierror)
  if (nrank==0) write(OUTPUT_UNIT,"(3A,E15.8)") "Valeur coefficient dans le fichier ",&
                      trim(file_name)," pour la zone 1 = ", valueSum/nSum

  tmpValue = 0
  i=0
  if(ind_zone_max == nMax) then
     tmpValue = array(maxloc(ind_zone,1))
     i=1
  end if
  call MPI_Reduce(tmpValue,valueSum,1,real_type,MPI_SUM,0,amitex_mpi_comm,ierror)
  call MPI_Reduce(i,nSum,1,MPI_INTEGER8,MPI_SUM,0,amitex_mpi_comm,ierror)
  if (nrank==0) write(OUTPUT_UNIT,"(3A,I8,A,E15.8)") "Valeur coefficient dans le fichier ",&
                      trim(file_name)," pour la derniere zone ",nMax ," = ", valueSum/nSum


  !! Fermeture du fichier
  call MPI_FILE_CLOSE(FU,ierror)
!  if(testEndian2==1 .AND. nOctets/=1) call MPI_TYPE_FREE(type_end,ierror)

  !! Desallocations 
  if (allocated(tmp_double)) deallocate(tmp_double)  
  if (allocated(tmp_real)) deallocate(tmp_real)  
  if (allocated(tmp_int1)) deallocate(tmp_int1)  
  if (allocated(tmp_int2)) deallocate(tmp_int2)  
  if (allocated(tmp_int4)) deallocate(tmp_int4)  
  if (allocated(tmp_int8)) deallocate(tmp_int8)  
!  if (allocated(blocklen)) deallocate(blocklen)  
!  if (allocated(dplt)) deallocate(dplt)  

end subroutine getBinCoeffFromFile


!==============================================================================
!
!                   SUBROUTINE GETCOEFFFROMFILE : 
!> Lecture dans un fichier ASCII les valeurs d'un coefficient (variables internes)
!> pour un materiau selon la zone (cas constant par zone) 
!!
!! \param[in]   file_name: (chaine de caracteres) nom du fichier a lire
!! \param[in]   ind_zone: (entier) tableau des indices des zones presentes dans le materiau
!! \param[in]   nZoneP: (entier) nombre de zones presentes dans le materiau
!! \param[in]   FU: (entier) unite du fichier ou sont lus les coefficients
!! \param[out]  array: (reel) tableau des coefficients/ variables internes
!!
!!
!! ATTENTION \n
!! On peut definir au maximum 2^31-1 zones
!!
!==============================================================================
subroutine getCoeffFromFile(file_name, ind_zone, nZoneP, array, FU )

  implicit none

  character(len=*),intent(in)                   :: file_name
  integer(kind=8),intent(in)                    :: nZoneP
  integer(kind=8),dimension(nZoneP),intent(in)  :: ind_zone
  !! Inverse du tableau ind_zone
  integer(kind=8),dimension(:),allocatable      :: inv_ind_zone  
  real(mytype),dimension(nZoneP), intent(out)   :: array

  integer, intent(in)                   :: FU
  integer(kind=8)                       :: i,j, nZone, nMax
  character(len=200)                    :: err_msg
  integer                               :: ierror, ios
  logical                               :: file_exists,empty_rank

  inquire(FILE=file_name,EXIST=file_exists)
  if(.not.file_exists)then 
    write(err_msg,fmt="(3A)") "Le fichier: '",trim(file_name),"' n'existe pas ,(getCoeffFromFile )"
    call amitex_abort(err_msg,1,0)
    return
  end if
  
  !! Determination prealable de la presence ou non de donnes a lire sur le pinceau
  !  leur absence est signalee par des tableaux vides en entree et nZoneP = 0
  empty_rank = .false.
  if (nZoneP == 0) empty_rank = .true.
  
  !! nZone est le numero de zone maximal de ce materiau dans ce processus
  if(.not. empty_rank) then
    nZone = ind_zone(nZoneP)
  else
    nZone=0
  end if

  call MPI_Allreduce(nZone,nMax,1,MPI_INTEGER8,MPI_MAX,amitex_mpi_comm,ierror)
  !! nMax est le numero de zone maximal de ce materiau 
  !! ie. le nombre de zones dans le materiau
  nZone=int(2,kind=8)**31-1
  if(nMax>nZone)then
     write(err_msg,fmt="(A)") "Nombre de zones trop important &
          & (max=2^31 -1) (getCoeffFromFile)"
     call amitex_abort(err_msg,1,0)
     return
  end if
  !! On inverse le tableau ind_zone
  !! inv_ind_zone permet de trouver le numero de zone local au processeur a partir
  !! du numero de zone global
  allocate(inv_ind_zone(nMax))
  inv_ind_zone = -1
  if (.not. empty_rank) then
     do i=1,nZoneP
       inv_ind_zone(ind_zone(i))=i
     end do
  end if
  
  open(unit=FU, file=trim(file_name)&
       ,action="read")
  !! On boucle sur les zones du materiau
  do i=1,nMax
     j = inv_ind_zone(i)
     if(j>0) then
        !! Si la zone est presente dans ce processus on recuper les donnees de la ligne
        read(FU,"(F17.0)",iostat=ios) array(j)
        if(ios/=0) then
           write(err_msg,fmt="(A,I0)") &
                "Erreur lecture coefficient. ios = ", ios
           call amitex_abort(err_msg,1,0)
           return
        end if
     else
        !! Si cette zone n'est pas dans le processus on saute la ligne
        read(FU,*,iostat=ios)
        if(ios/=0) then
           write(err_msg,fmt="(A,I0)") &
                "Erreur lecture coefficient. ios = ", ios
           call amitex_abort(err_msg,1,0)
           return
        end if
     end if
  end do
  close(FU)
  deallocate(inv_ind_zone)

end subroutine getCoeffFromFile
!==============================================================================


!==============================================================================
!
!  SUBROUTINE GET_REAL_XML : 
!> Permet de lire la valeur d'un attribut (reel) dans un fichier xml,
!> ne modifie pas la valeur si le nom de l'attribut n'est pas correct
!!
!!
!! \param[in]   cur_node: (type node) noeud dans lequel se trouve l'attribut
!! \param[in]   attribute_name: (chaine de caracteres) nom de l'attribut
!! \param[in]   err_code: (entier) niveau d'erreur en cas d'absence de la valeur
!! \param[in]   r: (entier) optionnel, rang du processus implique
!!                       (si r est absent, tout les processus sont impliques)
!!
!! \param[out] value: (reel) valeur de l'attribut
!
!==============================================================================
subroutine get_real_xml(cur_node,attribute_name,value,err_code,r)

  implicit none
  
  type(node), pointer,intent(in)  :: cur_node
  character(len=*), intent(in)    :: attribute_name
  real(mytype),intent(out)        :: value
  integer,intent(in)              :: err_code
  integer,intent(in),optional     :: r

    if(isPresent(cur_node,attribute_name,err_code,r)) then
      call extractDataAttribute(cur_node,attribute_name,value)
    end if

end subroutine get_real_xml

!==============================================================================
!
!  SUBROUTINE GET_REAL_XML_DEFAULT : 
!> Permet de lire la valeur d'un attribut (reel) dans un fichier xml
!! en attribuant une valeur par defaut
!!
!! \param[in]   cur_node: (type node) noeud dans lequel se trouve l'attribut
!! \param[in]   attribute_name: (chaine de caracteres) nom de l'attribut
!! \param[in]   err_code: (entier) niveau d'erreur en cas d'absence de la valeur
!! \param[in]   default : (reel) valeur par defaut
!! \param[in]   r: (entier) optionnel, rang du processus implique
!!                       (si r est absent, tout les processus sont impliques)
!! \param[out]  value: (reel) valeur obtenue
!
!==============================================================================
subroutine get_real_xml_default(cur_node,attribute_name,value,err_code,default,r)

  implicit none

  integer,intent(in),optional     :: r  
  type(node), pointer,intent(in)  :: cur_node
  character(len=*), intent(in)    :: attribute_name
  real(mytype),intent(in)         :: default
  real(mytype),intent(out)        :: value
  integer,intent(in)              :: err_code
  character(len=50)               :: attribute

  if(isPresent(cur_node,attribute_name,err_code,r)) then
     attribute = getAttribute(cur_node,attribute_name)
     call stringToLower(attribute)
     if(attribute == "default") then
        value = default        
     else
        call extractDataAttribute(cur_node,attribute_name,value)
     end if
  end if

end subroutine get_real_xml_default

!===================================================================
!
!                           SUBROUTINE GET_INT_XML
!
!>  permet de lire la valeur d'un attribut (entier) dans un fichier xml
!!
!! \param[in]  cur_node: (type node) noeud dans lequel se trouve l'attribut
!! \param[in]  attribute_name: (chaine de caracteres) nom de l'attribut
!! \param[in]  err_code: (entier) niveau d'erreur en cas d'absence de la valeur
!! \param[in]  r: (entier) optionnel, rang du processus implique
!!                (si r est absent, tous les processus sont impliques)
!!
!! \param[out]   value: (entier) valeur de l'attribut
!!
!===================================================================
subroutine get_int_xml(cur_node,attribute_name,value,err_code,r)

  implicit none
  
  type(node), pointer, intent(in) :: cur_node
  character(len=*), intent(in)    :: attribute_name
  integer, intent(out)            :: value
  integer,intent(in)              :: err_code
  integer,intent(in),optional     :: r

    if(isPresent(cur_node,attribute_name,err_code,r)) then
       call extractDataAttribute(cur_node,attribute_name,value)
    end if

end subroutine get_int_xml

!===================================================================
!
!                   SUBROUTINE GET_STR_XML
!
!>  permet de lire la valeur d'un attribut (chaine de caracteres) dans un fichier xml
!!  la chaine de caractere est transformée en minuscules
!!
!! \param[in]  cur_node: (type node) noeud dans lequel se trouve l'attribut
!! \param[in]   attribute_name: (chaine de caracteres) nom de l'attribut
!! \param[in]   err_code: (entier) niveau d'erreur en cas d'absence de la valeur
!! \param[in]   r: (entier) optionnel, rang du processus implique
!!              (si r est absent, tous les processus sont impliques)
!! \param[out]   value: (chaine de caracteres) valeur de l'attribut
!
!===================================================================
subroutine get_str_xml(cur_node,attribute_name,value,err_code,r)

  implicit none
  
  type(node), pointer, intent(in) :: cur_node
  character(len=*), intent(in)    :: attribute_name
  character(len=*), intent(out)   :: value
  integer,intent(in)              :: err_code
  integer,intent(in),optional     :: r

  if(isPresent(cur_node,attribute_name,err_code,r)) then
     value =getAttribute(cur_node,attribute_name)
  end if
   
  call stringToLower(value)
end subroutine get_str_xml

!===================================================================
!
!                   SUBROUTINE GET_STR_XML_NOLOWERCASE
!
!>  permet de lire la valeur d'un attribut (chaine de caracteres) dans un fichier xml
!!  la chaine de caractere n'est pas transformée en minuscules (No Lower Case)
!!
!! \param[in]  cur_node: (type node) noeud dans lequel se trouve l'attribut
!! \param[in]   attribute_name: (chaine de caracteres) nom de l'attribut
!! \param[in]   err_code: (entier) niveau d'erreur en cas d'absence de la valeur
!! \param[in]   r: (entier) optionnel, rang du processus implique
!!              (si r est absent, tous les processus sont impliques)
!! \param[out]   value: (chaine de caracteres) valeur de l'attribut
!
!===================================================================
subroutine get_str_xml_nolowercase(cur_node,attribute_name,value,err_code,r)

  implicit none
  
  type(node), pointer, intent(in) :: cur_node
  character(len=*), intent(in)    :: attribute_name
  character(len=*), intent(out)   :: value
  integer,intent(in)              :: err_code
  integer,intent(in),optional     :: r

  if(isPresent(cur_node,attribute_name,err_code,r)) then
     value =getAttribute(cur_node,attribute_name)
  end if
   
end subroutine get_str_xml_nolowercase


!===================================================================
!
!                   SUBROUTINE GET_STR_XML_DEFAULT
!
!> Permet de lire la valeur d'un attribut (chaine de caracteres) dans un fichier xml
!! en attribuant une valeur par defaut
!!
!! \param[in]   cur_node: (type node) noeud dans lequel se trouve l'attribut
!! \param[in]   attribute_name: (chaine de caracteres) nom de l'attribut
!! \param[in]   err_code: (entier) niveau d'erreur en cas d'absence de la valeur
!! \param[in]   default : (reel) valeur par defaut
!! \param[in]   r: (entier) optionnel, rang du processus implique
!!              (si r est absent, tous les processus sont impliques)
!! \param[out]   value: (chaine de caracteres) valeur de l'attribut
!
!===================================================================
subroutine get_str_xml_default(cur_node,attribute_name,value,err_code,default,r)

  implicit none
  
  type(node), pointer, intent(in) :: cur_node
  character(len=*), intent(in)    :: attribute_name
  character(len=*), intent(out)   :: value
  integer,intent(in)              :: err_code
  character(len=*),intent(in)     :: default
  integer,intent(in),optional     :: r

  if(isPresent(cur_node,attribute_name,err_code,r)) then
     value = getAttribute(cur_node,attribute_name)
     call stringToLower(value)
     if(value=="default") value = default
  end if

end subroutine get_str_xml_default

!==============================================================================
!
!                       FUNCTION ISPRESENT
!
!>  Permet de savoir si un attribut est present dans un noeud (xml) renvoie "TRUE" si l'attribut est present, "FALSE" sinon
!!
!! \param[in] node_tested: (type node) noeud (xml) dans lequel on cherche l'attribut
!! \param[in]   attribute_tested: (chaine de caracteres) nom de l'attribut recherche
!! \param[in]   err_code: (entier) niveau d'erreur en cas d'absence de l'attribut
!! \param[in]   r: (entier) optionnel, rang du processus implique pour l'erreur
!
!==============================================================================
function isPresent(node_tested, attribute_tested, err_code,r)

  implicit none

  type(node), pointer, intent(in) :: node_tested
  character(len=*),intent(in)     :: attribute_tested
  integer,intent(in)              :: err_code
  integer,intent(in),optional     :: r

  type(node), pointer             :: att
  logical                         :: isPresent
  character(len=len(attribute_tested)+100)              :: err_msg


  att => getAttributeNode(node_tested,trim(attribute_tested))
  if(.not.associated(att))then
    isPresent = .False.
    write(err_msg,fmt="(3A)") "attribut: '", trim(attribute_tested), "' absent (isPresent)"
    call amitex_abort(trim(err_msg), err_code,r)
  else
    isPresent = .True.
  end if  

end function isPresent




!====================================================================
!
!                       SUBROUTINE  GETNODELIST 
!> Renvoie la liste des noeuds, portant un nom donne, contenus dans un noeud "parent" 
!!
!!
!! \param[in] parent: (type node) noeud dans lequel on recherche la liste
!! \param[in]   list_name: (chaine de caracteres) nom des noeuds recherches
!! \param[in]   n: (entier) taille de la liste
!!                   ( si ,n=-1, aucune taille n'est specifie)
!! \param[in]   err_code: (entier) en cas d'erreur, niveaude l'erreur (warnig, error...)
!!
!! \param[out]   list: (type nodeList) liste des noeuds du noeud "parent" portant le nom
!!               list_name
!!
!====================================================================
subroutine getNodeList(parent, list, list_name, n, err_code)

  implicit none

  type(node),pointer, intent(in)      :: parent
  type(nodeList),pointer, intent(out) :: list
  character(len=*), intent(in)        :: list_name
  integer, intent(in)                 :: n, err_code
  
  character(len=150) :: err_msg
  
  list => getElementsByTagName(parent,list_name)
  if(getLength(list)/=n .and. n>-1) then
    write(err_msg,fmt="(A,I0,3A)") "exactement ",n," noeud(s) '",list_name,"' devrai(en)t etre present(s) (GetNodeList)"
    call amitex_abort(err_msg,err_code)
  end if

  !if(getLength(list) .eq. 0) then  !! tentative pour rendre le noeud output optionnel 
  !  list => null()
  !end if

  
end subroutine GetNodeList

!====================================================================
!
!> Convertit une chaine de caracteres en minuscules
!!
!!
!! \param[in] string: Chaine de caracteres initiale
!! \param[out] string: Chaine de caracteres uniquement avec des minuscules
!!
!====================================================================
subroutine stringToLower( string )

  implicit none
    character(len=*), intent(inout) :: string 
    integer                          :: i, k

    do i = 1,len(string) 
        k = iachar(string(i:i)) 
        if ( k >= iachar('A') .and. k <= iachar('Z') ) then 
            k = k + iachar('a') - iachar('A') 
            string(i:i) = achar(k) 
        endif 
    enddo     
end subroutine stringToLower 




!***********************************************************************************************************************************
!
!                                                        B Y T E S W A P
!
!  Module:       BYTESWAP
!
!  Programmer:   Dr. David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!
!  Date:         July 23, 2005
!
!  Language:     Fortran-90
!
!  Description:  This module includes several subroutines to reverse the byte order of INTEGER and REAL variables.
!                This is useful when reading binary data from a file intended for use on a compute whose byte
!                order is opposite that of the computer on which the Fortran program is to be run.
!
!                   SWAP_I2    Swap bytes of a 2-byte INTEGER
!                   SWAP_I4    Swap bytes of a 4-byte INTEGER
!                   SWAP_I8    Swap bytes of a 4-byte INTEGER (ajout LG)
!                   SWAP_F4    Swap bytes of a 4-byte REAL
!                   SWAP_F8    Swap bytes of a 8-byte REAL
!
!***********************************************************************************************************************************

!***********************************************************************************************************************************
!  SWAP_I2
!
!  Swap bytes for a two-byte INTEGER.
!  After calling this subroutine, the input integer will be replaced by the output integer.
!***********************************************************************************************************************************

      ELEMENTAL SUBROUTINE SWAP_I2 (BYTE2)

      IMPLICIT NONE

      INTEGER(KIND=2), INTENT(IN OUT) :: BYTE2

      INTEGER(KIND=1), DIMENSION(2) :: BYTE_ARR, BYTE_ARR_TMP
      INTEGER :: I


      BYTE_ARR = TRANSFER (BYTE2, BYTE_ARR)

      BYTE_ARR_TMP = BYTE_ARR

      DO I = 1, 2
         BYTE_ARR(I) = BYTE_ARR_TMP(3-I)
      END DO

      BYTE2 = TRANSFER (BYTE_ARR, BYTE2)

      RETURN

      END SUBROUTINE SWAP_I2



!***********************************************************************************************************************************
!  SWAP_I4
!
!  Swap bytes for a four-byte INTEGER.
!  After calling this subroutine, the input integer will be replaced by the output integer.
!***********************************************************************************************************************************

      ELEMENTAL SUBROUTINE SWAP_I4 (BYTE4)

      IMPLICIT NONE

      INTEGER(KIND=4), INTENT(IN OUT) :: BYTE4

      INTEGER(KIND=1), DIMENSION(4) :: BYTE_ARR, BYTE_ARR_TMP
      INTEGER :: I


      BYTE_ARR = TRANSFER (BYTE4, BYTE_ARR)

      BYTE_ARR_TMP = BYTE_ARR

      DO I = 1, 4
         BYTE_ARR(I) = BYTE_ARR_TMP(5-I)
      END DO

      BYTE4 = TRANSFER (BYTE_ARR, BYTE4)

      RETURN

      END SUBROUTINE SWAP_I4




!***********************************************************************************************************************************
!  SWAP_I8
!
!  Swap bytes for a eight-byte INTEGER.
!  After calling this subroutine, the input integer will be replaced by the output integer.
!***********************************************************************************************************************************

      ELEMENTAL SUBROUTINE SWAP_I8 (BYTE8)

      IMPLICIT NONE

      INTEGER(KIND=8), INTENT(IN OUT) :: BYTE8

      INTEGER(KIND=1), DIMENSION(8) :: BYTE_ARR, BYTE_ARR_TMP
      INTEGER :: I


      BYTE_ARR = TRANSFER (BYTE8, BYTE_ARR)

      BYTE_ARR_TMP = BYTE_ARR

      DO I = 1, 8
         BYTE_ARR(I) = BYTE_ARR_TMP(9-I)
      END DO

      BYTE8 = TRANSFER (BYTE_ARR, BYTE8)

      RETURN

      END SUBROUTINE SWAP_I8
!***********************************************************************************************************************************
!  SWAP_F4
!
!  Swap bytes for a four-byte REAL.
!  After calling this subroutine, the input number will be replaced by the output number.
!***********************************************************************************************************************************

      ELEMENTAL SUBROUTINE SWAP_F4 (FLOAT4)

      IMPLICIT NONE

      REAL(KIND=4), INTENT(IN OUT) :: FLOAT4

      INTEGER(KIND=1), DIMENSION(4) :: BYTE_ARR, BYTE_ARR_TMP
      INTEGER :: I


      BYTE_ARR = TRANSFER (FLOAT4, BYTE_ARR)

      BYTE_ARR_TMP = BYTE_ARR

      DO I = 1, 4
         BYTE_ARR(I) = BYTE_ARR_TMP(5-I)
      END DO

      FLOAT4 = TRANSFER (BYTE_ARR, FLOAT4)

      RETURN

      END SUBROUTINE SWAP_F4



!***********************************************************************************************************************************
!  SWAP_F8
!
!  Swap bytes for an eight-byte REAL.
!  After calling this subroutine, the input number will be replaced by the output number.
!***********************************************************************************************************************************

      ELEMENTAL SUBROUTINE SWAP_F8 (FLOAT8)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN OUT) :: FLOAT8

      INTEGER(KIND=1), DIMENSION(8) :: BYTE_ARR, BYTE_ARR_TMP
      INTEGER :: I


      BYTE_ARR = TRANSFER (FLOAT8, BYTE_ARR)

      BYTE_ARR_TMP = BYTE_ARR

      DO I = 1, 8
         BYTE_ARR(I) = BYTE_ARR_TMP(9-I)
      END DO

      FLOAT8 = TRANSFER (BYTE_ARR, FLOAT8)

      RETURN

      END SUBROUTINE SWAP_F8


end module io_amitex_mod
