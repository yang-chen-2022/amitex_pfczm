!123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!===================================================================================================
!                         SUBROUTINE READ GEOM CONTINUOUS
!
!> Construction partielle de MattotP0 : un seul materiau et une zone par voxel 
!!
!!  creer tableau MattotP0 (directement reparti sur les processus)
!!       - numero de materiau                   (MattotP0%numM )
!!       - tableau de positions triees par zone  (MattotP0%pos )
!!       - numero de zone present dans le pinceau et indice max de la zone
!!                                              (MattotP0%zone )
!! ICI:
!!      numM=1
!!      pos(:)   : indice de position local : 1:nbvoxP
!!      zone(:,1): indice de position local : 1:nbvoxP
!!      zone(:,2): numeros de zone globale  
!!                 (defini selon le parcours des tableaux en fortran) 
!!      (nbvoxP = nombre de voxels dans le pinceau)
!!
!!
!! \param[in] nx,ny (nz n'est pas necessaire, on utilise le decoupage 2decomp)
!!
!===================================================================================================
subroutine read_geom_continuous(nx,ny)

  use ISO_FORTRAN_ENV

  use error_mod
  use material_mod
#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype
#else
  use decomp_2d , only : mytype, nrank
#endif  

  use field_mod, only : amx_start, amx_end, amx_size

  implicit none

  integer,intent(in)    :: nx,ny
  integer(kind=8)       :: nbvoxP,posG, posP
  integer               :: alloc_stat,i,j,k

  nbvoxP = amx_size(1)*amx_size(2)*amx_size(3)

  !> Allocation de MattotP0 (un seul materiau)
  allocate(MattotP0(1),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_geom_continuous)",2)

  !> Allocation du nombre de "position" (tableau pos = nombre de voxel dans le pinceau)
  allocate(MattotP0(1)%pos(nbvoxP),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_geom_continuous)",2)

  !> Allocation du nombre de "zone" (tableau zone = nombre de voxel dans le pinceau)
  allocate(MattotP0(1)%zone(nbvoxP,2),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_geom_continuous)",2)

  !> Affectation du materiau
  MattotP0(1)%numM = 1

  !> Affectation des positions et des zones 
  !!    zone(i,1) : indice max de la zone zone(i,2)
  !!    zone(i,2) : numero de la zone (numZ)
  !! ICI:
  !!    zone(i,1) : indice de position (1:nbvoxP)
  !!    zone(i,2) : numero de zone globale  (defini selon le parcours des tableaux en fortran) 

  posP = 0
  do k=amx_start(3),amx_end(3)         
  do j=amx_start(2),amx_end(2)         
  do i=amx_start(1),amx_end(1)    
     posP = PosP + 1     
     posG = i + ((j-1)*nx) + ((k-1)*nx*ny)
     MattotP0(1)%pos(posP) = posP
     MattotP0(1)%zone(posP,1) = posP
     MattotP0(1)%zone(posP,2) = posG
  end do 
  end do 
  end do 

  nmateriaux0%n = 1 
  
end subroutine read_geom_continuous
!==================================================================================
!==================================================================================
!                         SUBROUTINE READ GEOM  
!
!> Lecture de la geometrie et construction partielle de MattotP0
!!
!!  A partir de la lecture des fichers vtk (fic_numM et fic_numZ),
!!  creer tableau MattotP0 (directement reparti sur les processus)
!!       - numero de materiau                   (MattotP0%numM )
!!       - tableau de positions triee par zone  (MattotP0%pos )
!!       - numero de zone present dans le pinceau et indice max de la zone
!!                                              (MattotP0%zone )
!!
!! \param[in] fic_numM: (chaine de caracteres) nom du fichier vtk 
!!                      decrivant la geometrie des numeros de materiau
!! \param[in] fic_numZ: (chaine de caracteres) nom du fichier vtk 
!!                       decrivant la geometrie des numeros de zone
!! \param[in] type_mpi_m,type_mpi_z: (entier) type de donnees lues mpi_integer2 ou mpi_integer4
!!                          dans les fichiers (fic_numM, fic_numZ)
!!                          depend de l'entete du fichier : type short ou int
!!
!! \param[out]   numM, numZ: (entiers) champs de numeros de materiau et de zone
!! 
!! Modifie aussi MattotP0, evalue nmateriaux0 (module material_mod)
!!
!==================================================================================
subroutine read_geom(fic_numM, fic_numZ, type_mpi_m, type_mpi_z, numM, numZ)

  use ISO_FORTRAN_ENV

  use material_mod
#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi,       only : nrank
#else
  use decomp_2d , only : mytype, nrank
#endif 
  use field_mod, only : amx_start, amx_end 
  use io_amitex_mod
  use error_mod

  implicit none

  character(len=*), intent(in)  :: fic_numM, fic_numZ
  integer, intent(in)           :: type_mpi_m, type_mpi_z
  integer,dimension(amx_start(1):amx_end(1),amx_start(2):amx_end(2),amx_start(3):amx_end(3)),&
                    intent(out) ::numM
  integer(kind=8),dimension(amx_start(1):amx_end(1),amx_start(2):amx_end(2),amx_start(3):amx_end(3)),&
                    intent(out) ::numZ

!> tableau temporaire pour la creation de MattotP0
  integer(kind=8),allocatable, dimension(:,:)::  tailles_mattotP

  integer ::  nMat, ierror
  integer(kind=8) :: nZone, i,j,k,l,m, nTot
  integer(kind=8),allocatable,dimension(:) :: ind_zones,temp_tailles

  character(len=200) :: err_msg
  integer :: alloc_stat

  numM=-1
  numZ=-1

  ! lecture des fichiers (s'ils existent)
  ! fichier materiaux
  if(fic_numM/="") then 
    call read_bin4(trim(fic_numM), 1,type_mpi_m, numM) !numM n'est pas au format "vecteur" (numM(:))
                                                       !on ne peut utiliser l'interface generique read_bin
    i=minval(numM)
!    call MPI_Allreduce(i, j, 1, MPI_INTEGER8, MPI_MIN, amitex_mpi_comm, ierror)
    call MIN_MPI_I8(i,j,ierror)
    ! indice minimal de materiau: doit etre "0" ou "1"
    if(j==0)then
      ! numM reindexe a partir de "1"
      numM=numM+1
    else if(j<0.or.j>1) then
      write(err_msg,fmt="(A,I0,A)") "Numero de materiau minimal invalide: ",j," (attendu: 0 ou 1), (read_geom)"
      call amitex_abort(err_msg,1)
    end if
  else 
    numM=1
    nMat=1
  end if

  ! fichier zone
  if(fic_numZ/="") then 
    call read_bin8(trim(fic_numZ), 1,type_mpi_z, numZ) !numZ n'est pas au format "vecteur" (numZ(:))
                                                       !on ne peut utiliser l'interface generique read_bin
    i=minval(numZ)
!    call MPI_Allreduce(i, j, 1, MPI_INTEGER8, MPI_MIN, amitex_mpi_comm, ierror)
    call MIN_MPI_I8(i,j,ierror)
    ! indice minimal de zone: doit etre "0" ou "1"
    if(j==0)then
      ! numZ reindexe a partir de "1"
      numZ=numZ+1
    else if(j<0.or.j>1)then
      write(err_msg,fmt="(A,I4,A)") "Numero de zone minimal invalide: ",j," (attendu: 0 ou 1) (read_geom)"
      call amitex_abort(err_msg,1)
    end if
  else
    numZ=1
    nZone=1
  end if

  call check_amitex_abort()

! nombre de materiaux et de zones  
  nMat=maxval(numM)
  i = maxval(numZ)
!  call MPI_Allreduce(nMat, nmateriaux0%n, 1, MPI_INTEGER, MPI_MAX, amitex_mpi_comm, ierror)
!  call MPI_Allreduce(i, nZone, 1, MPI_INTEGER8, MPI_MAX, amitex_mpi_comm, ierror)
  call MAX_MPI_I(nMat, nmateriaux0%n,ierror)
  call MAX_MPI_I8(i, nZone,ierror)

! tailles_mattotp2(i,1) : nombre de voxels du materiau i dans le pinceau
! tailles_mattotp2(i,2) : nombre de zones du materiau  i dans le pinceau
! tailles_mattotp2(i,j+2) : nombre de voxels pour la zone j dans le pinceau
  allocate(tailles_mattotP(nmateriaux0%n,2+nZone),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_geom)",2)
  tailles_mattotP=0

! On compte le nombre de voxels par materiau et par zone
  do k= amx_start(3),amx_end(3)
    do j= amx_start(2),amx_end(2)
      do i= amx_start(1),amx_end(1)
        ! 1ere colonne : compte du nombre de voxels pour chaque materiau
        tailles_mattotP(numM(i,j,k), 1) = tailles_mattotP(numM(i,j,k), 1)+ 1
        ! colonnes suivantes : compte du nb de voxels pour chaque zone du materiau
        tailles_mattotP(numM(i,j,k),numZ(i,j,k)+2) = tailles_mattotP(numM(i,j,k),numZ(i,j,k)+2)+1
      end do
    end do
  end do

  j=0
  l=0
  do i=1,nmateriaux0%n
    !  j: nombre de materiau dans le pinceau
    if(tailles_mattotP(i,1)>0)j=j+1
    do k=1,nZone
    !  nombre de zone du materiau i dans le pinceau (stocke dans tailles_mattotP(i,2))
      if(tailles_mattotP(i,k+2)>0) tailles_mattotP(i,2)=tailles_mattotP(i,2)+1
    end do
  end do

! recapitulatif du comptage ci-dessus  
! tailles_mattotp2(i,1) : nombre de voxels du materiau i dans le pinceau
! tailles_mattotp2(i,2) : nombre de zones du materiau  i dans le pinceau
! tailles_mattotp2(i,j+2) : nombre de voxels pour la zone j dans le pinceau
  

  ! taille de mattotp: nombre de materiau dans le pinceau
  allocate(MattotP0(j),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_geom)",2)
    nMat = 0 ! numero de materiau dans le pinceau

    ! determination de la presence du materiau dans le pinceau, du nombre de
    ! voxels associes, construction du tableau de position 'pos'  (ordonne par zone)
    do i=1,nmateriaux0%n
      !si le materiau i existe dans le pinceau:
      if(tailles_mattotp(i,1)>0 ) then
      
        ! numero du materiau (dans le pinceau) correspondant au materiau i (global)
        nMat= nMat+1    
        MattotP0(nMat)%numM=int(i,kind=4)
        ! allocation du nombre de "position" (tableau pos)
        allocate(MattotP0(nMat)%pos(tailles_mattotP(i,1 )),stat=alloc_stat)
        if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_geom)",2)

        ! allocation du tableau d'indice des zones
        allocate(MattotP0(nMat)%zone(tailles_mattotp(i,2) ,2),stat=alloc_stat)
        if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_geom)",2)
        ! mattotp2(nmat)%zone(i,1):nombre de voxels de la zone
        ! mattotp2(nmat)%zone(i,2):numero de zone (non parallelise)

        ! affectation des indices max des zones
        l=0
        m=0
        do k=3,size(tailles_mattotp(i,:),kind=8)  ! boucle sur les positions
          if(tailles_mattotp(i,k)>0)then          ! presence de la zone k-2 (materiau i) dans le pinceau
            l=l+1                                 ! numero de la zone parallelise (pour le pinceau)
            m= m + tailles_mattotp(i,k)           ! nombre de voxels cumules (de toutes les zones deja parcourues)
            MattotP0(nMat)%zone(l,1)= m            ! indice maximal de la zone dans le tableau pos
            MattotP0(nMat)%zone(l,2)= k-2          ! numero de zone (non parallelise)
          end if
        end do

        allocate(ind_zones(tailles_mattotp(i,2)),stat=alloc_stat) ! ind_zone: indice courant de chaque zone
        if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_geom)",2)
        ind_zones(1)=1
        if(tailles_mattotp(i,2)>1)then
          !ind_zones( zone )= indice max (zone-1) +1 (premier indice de zone)
          ind_zones(2:tailles_mattotp(i,2)) = MattotP0(nMat)%zone(1:tailles_mattotp(i,2)-1,1)+1
        end if

        ! numero de l'element a placer dans le tableau 'pos'
        nTot = 1
        ! affectation des positions
        do j=amx_start(3),amx_end(3)
          do k=amx_start(2),amx_end(2)
            do l=amx_start(1),amx_end(1)
              if(numM(l,k,j)==i)then
                nzone = numZ(l,k,j) ! numero de zone a trouver
                m=1  ! indice de la zone actuelle
                ! recherche de la zone correspondante (parallelise)
                do while(MattotP0(nMat)%zone(m,2)/=numZ(l,k,j) )
                  m=m+1 
                end do
                ! position courante
                MattotP0(nMat)%pos(ind_zones(m))=nTot
                ! nouvel indice courant de zone
                ind_zones(m)=ind_zones(m)+1
              end if
              nTot=nTot+1      
            end do
          end do
        end do
        if(allocated(ind_zones)) deallocate(ind_zones)
         
      end if
   end do

! verification de la continuite des nombres de materiaux/zones
    nzone = size(tailles_mattotp(1,:))-2
    allocate(ind_zones(nzone),stat=alloc_stat)  ! ind_zone(i): nombre de voxels pour la zone i, NON PARALLELISE (MPI_Reduce, MPI_SUM) 
    if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_geom_composite)",2)
    allocate(temp_tailles(nzone),stat=alloc_stat)
    if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_geom_composite)",2)
    
  do i=1,nmateriaux0%n
    ! SOLUTION INITIALE
    ! nzone = size(tailles_mattotp(1,:))-2
    ! allocate(ind_zones(nzone),stat=alloc_stat)  ! ind_zone(i): nombre de voxels pour la zone i, NON PARALLELISE (MPI_Reduce, MPI_SUM) 
    ! call MPI_Allreduce(tailles_mattotp(i,3:nZone+2),ind_zones(:), nzone, MPI_INTEGER8, MPI_SUM, amitex_mpi_comm,ierror)
    !
    ! 2 PROBLEMES ICI : 1/ nzone en kind=8 (incompatible openMPI -> kind=4)
    !                   2/ le passage de 'tailles_mattotp(i,3:nZone+2)' en argument créé un problème
    !                      origine mal comprise...  possibilité : donnees non-contigues...
    !                      fonctionne en creant un tableau temporaire 1D (de donnees contigues...)
    ! 
    ! SOLUTION ACTUELLE pour compter les voxels en une seule somme parallele par materiau : 
    ! on passe par la fonction SUM_MPI_I8_LARGE qui permet de realiser des sommes paralleles avec un nombre 
    ! d'elements sueprieur a la valeur max des int_kind=4 (2^31 -1) 
    ! ATTENTION : pour les vecteurs depassant effectivement cette taille, la fonction n'est pas encore
    !            implementee (solution de secours pour le moment)
    ! ATTENTION : peut être passer l'allocation et la definition de temp_tailles dans SUM_MPI_I8_LARGE

    temp_tailles(:) = tailles_mattotp(i,3:nZone+2)
    call SUM_MPI_I8_LARGE(temp_tailles,ind_zones(:),nzone,ierror)    

    if(nrank==0) then   
      if(all(ind_zones .eq. 0))then
        write(err_msg,fmt="(A,I0,A)") "Materiau ",i," absent (lecture fichier vtk) (read_geom)"
        call amitex_abort(err_msg,1,0) 
      end if
      k=0
      do j=1,nzone ! size(ind_zones)
        if(ind_zones(j)==0 .and. k==0)then
          ! la zone portant le numero j est absente
          ! si une zone de numero > j est presente il s'agit d'une erreur
          k=j
        else if(ind_zones(j)>0 .and. k>0) then
          ! cas d'une zone absente
          if(j-1>k)then
            write (err_msg,fmt="(I0,A,I0)") k," a ",j-1
          else
            write (err_msg,fmt="(I0)") k
          end if
          write(err_msg,fmt="(A,I0,3A)") "Materiau ",i,", zone(s) ", trim(err_msg),&
                                         " absente(s) (lecture fichier vtk)  (read_geom)"
          call amitex_abort(err_msg,1,0)
          k=0
        else if(ind_zones(j)<0) then
          write(err_msg, fmt="(A,I0,A)") "Subroutine read_geom, ind_zone(j) negatif,j= ",j, " (read_geom)"
          call amitex_abort(err_msg,1,0)
        end if
      end do
    end if
  end do

  
  if(allocated(ind_zones)) deallocate(ind_zones)

  if(allocated(tailles_MattotP)) deallocate(tailles_MattotP)

  if (allocated(temp_tailles)) deallocate(temp_tailles)

  ! arret si erreur (ci-dessus)
  call check_amitex_abort()

end subroutine read_geom
!==================================================================================


!==================================================================================
!==================================================================================
!                         SUBROUTINE READ GEOM  COMPOSITE
!
!> Lecture de la geometrie et construction partielle de MattotP0
!!
!!  A partir de la lecture des fichers vtk (fic_numM et fic_numZ),
!!    et des variables de type COMPOSITE et INTERPHASE 
!!  creer tableau MattotP0 (directement reparti sur les processus)
!!       - numero de materiau                   (MattotP0%numM )
!!       - tableau de positions triee par zone  (MattotP0%pos )
!!       - numero de zone present dans le pinceau etnombre de position par zone
!!                                              (MattotP0%zone )
!!       Pour reperer les "interphases" (materiau ou zone n'existant pas sous forme 
!!                                      de voxel "homogene")
!!       - booleen si le materiau est une interphase     (MattotP0%Interphase)
!!       - liste des zones d'interphases (s'il y en a)   (MattotP0%Zones_Interphase)
!!       Pour definir les voxels composites :
!!       - nombre de phase pour les materiaux composites (MattotP0%Nphase)
!!       - numeriau de materiau et de zone associe a chaque phase du materiau composite
!!                                    (MattotP0%numM_composite & MattotP0%numZ_composite)
!!
!!
!! \param[in] fic_numM: (chaine de caracteres) nom du fichier vtk 
!!                      decrivant la geometrie des numeros de materiau
!! \param[in] fic_numZ: (chaine de caracteres) nom du fichier vtk 
!!                       decrivant la geometrie des numeros de zone
!! \param[in] type_mpi_m,type_mpi_z: (entier) type de donnees lues mpi_integer2 ou mpi_integer4
!!                          dans les fichiers (fic_numM, fic_numZ)
!!                          depend de l'entete du fichier : type short ou int
!!
!! \param[out]  numM, numZ: (entiers) champs de numeros de materiau et de zone
!! 
!! Modifie aussi MattotP0 (evalue aussi nmateriaux0 variable de material_mod)
!!
!!
!! 7/12/2020 MODIFICATION LG - Les donnees lourdes ne sont plus stockees dans MatComposite
!!                             on les lit "a la demande"
!!
!==================================================================================
 subroutine read_geom_composite(fic_numM,fic_numZ,fic_vtk,type_mpi_m,type_mpi_z,numM,numZ,nx,ny)

  use ISO_FORTRAN_ENV
  use mpi
  use material_mod
#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi,       only : nrank
#else
  use decomp_2d , only : mytype, nrank
#endif  
  use field_mod, only : amx_start, amx_end
  use io_amitex_mod
  use error_mod

  implicit none

  character(len=*), intent(in)                         :: fic_numM, fic_numZ, fic_vtk
  integer, intent(in)                                  :: type_mpi_m, type_mpi_z
  integer, intent(in)                                  :: nx,ny
  integer,dimension(amx_start(1):amx_end(1),amx_start(2):amx_end(2),amx_start(3):amx_end(3)),&
                    intent(out)                        :: numM
  integer(kind=8),dimension(amx_start(1):amx_end(1),amx_start(2):amx_end(2),amx_start(3):amx_end(3)),&
                    intent(out)                        :: numZ
!> tableau temporaire pour la creation de MattotP0
  integer(kind=8),allocatable, dimension(:,:)          :: tailles_mattotP

  integer(kind=INT32)                                  :: nMat, ierror, nmateriauxhom, nphases, i,j,k
  integer(kind=INT32)                                  :: indice_iph
  integer(kind=INT64)                                  :: nZone, i_64, j_64,k_64,l,m, nTot,indice_pos
  integer(kind=INT64)                                  :: nzone_iph  
  integer(kind=8),allocatable,dimension(:)             :: ind_zones,temp_tailles

  character(len=200)                                   :: err_msg
  character(len=400)                                   :: file_name
  character(len=100)                                   :: function_name
  integer                                              :: FU
  integer                                              :: alloc_stat
  logical                                              :: Iph,Iph_z
  
  real(mytype),allocatable,dimension(:)                :: array_real_1d,pos_globalr
  integer(INT64),allocatable,dimension(:)              :: array_int64_1d
  integer(kind=INT64)                                  :: NzoneP
  integer                                              :: iphi

  numM=-1
  numZ=-1
  function_name="read_geom_composite"

  ! lecture des fichiers (s'ils existent)
  ! fichier materiaux
  if(fic_numM/="") then 
    call read_bin4(trim(fic_numM), 1,type_mpi_m, numM) !numM n'est pas au format "vecteur" (numM(:))
                                                       !on ne peut utiliser l'interface generique read_bin
    i_64=minval(numM)
!    call MPI_Allreduce(i, j, 1, MPI_INTEGER8, MPI_MIN, amitex_mpi_comm, ierror)
    call MIN_MPI_I8(i_64,j_64,ierror)
    ! indice minimal de materiau: doit etre "0" ou "1"
    if(j_64==0)then
      ! numM reindexe a partir de "1"
      numM=numM+1
    else if(j_64<0.or.j_64>1) then
      write(err_msg,fmt="(A,I0,A)") "Numero de materiau minimal invalide: ",j_64,&
                                    " (attendu: 0 ou 1), (read_geom_composite)"
      call amitex_abort(err_msg,1)
    end if
  else 
    numM=1
    nMat=1
  end if

  ! fichier zone
  if(fic_numZ/="") then 
    call read_bin8(trim(fic_numZ), 1,type_mpi_z, numZ) !numM n'est pas au format "vecteur" (numM(:))
                                                       !on ne peut utiliser l'interface generique read_bin
    i_64=minval(numZ)
!    call MPI_Allreduce(i, j, 1, MPI_INTEGER8, MPI_MIN, amitex_mpi_comm, ierror)
    call MIN_MPI_I8(i_64,j_64,ierror)
    ! indice minimal de zone: doit etre "0" ou "1"
    if(j_64==0)then
      ! numZ reindexe a partir de "1"
      numZ=numZ+1
    else if(j_64<0.or.j_64>1)then
      write(err_msg,fmt="(A,I4,A)") "Numero de zone minimal invalide: ",j_64,&
                                    " (attendu: 0 ou 1) (read_geom_composite)"
      call amitex_abort(err_msg,1)
    end if
  else
    numZ=1
    nZone=1
  end if

  call check_amitex_abort()

! nombre de materiaux et de zones  
  nMat=maxval(numM)
  i_64 = maxval(numZ)
!  call MPI_Allreduce(nMat, nmateriaux0%n, 1, MPI_INTEGER, MPI_MAX, amitex_mpi_comm, ierror)
!  call MPI_Allreduce(i, nZone, 1, MPI_INTEGER8, MPI_MAX, amitex_mpi_comm, ierror)
  call MAX_MPI_I(nMat, nmateriaux0%n,ierror)
  call MAX_MPI_I8(i_64, nZone,ierror)
  
! Traitement des numeros de materiaux d'interphases non pris en compte : 
! nmateriaux vaut la valeur max du champ numM en entrée. Si un matériau d'interface
! est present, il se peut que son numM soit supérieur à nmateriaux 
! => on teste le numM le plus elevee present dans la structure Interphases
!    pour corriger si besoin
  if (allocated(Interphases)) then
     if (nmateriaux0%n < maxval(Interphases(:)%numM)) then
        nmateriaux0%n = maxval( (/ nmateriaux0%n,int(maxval(Interphases(:)%numM)) /) )
     end if
  end if


  ! Ajout du nombre de matériaux composites pour obtenir le total exact de matériaux dans la cellule
  nmateriauxhom    = nmateriaux0%n
  nmateriaux0%n = nmateriaux0%n + size(MatComposite)
  nmateriaux0%n_composites = nmateriaux0%n - nmateriauxhom
  
  ! Zones composites non présentes dans le champ numZ en entrée : mise à jour du nombre maximum de zones pour un 
  ! matériau dans la cellule
  do i=1,size(MatComposite)
!     if (size(MatComposite(i)%Zone(:,1)) > nZone) nZone = size(MatComposite(i)%Zone(:,1))
     if (MatComposite(i)%Nzones_global > nZone) nZone = MatComposite(i)%Nzones_global
  end do

! tailles_mattotp2(i,1) : nombre de voxels du materiau i dans le pinceau
! tailles_mattotp2(i,2) : nombre de zones du materiau  i dans le pinceau
! tailles_mattotp2(i,j+2) : nombre de voxels pour la zone j dans le pinceau
  allocate(tailles_mattotP(nmateriaux0%n,2+nZone),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_geom_composite)",2)
  tailles_mattotP=0

!  modification des champs numZ et numM : incorporation des numM et numZ des voxels composites
! On boucle sur les matériaux composites
  do i=1,nmateriaux0%n_composites
  
     ! Lecture Fichier position -> pos_globalr 
     !                             attention, getBinCoeffFromFile renvoie des reels
     !
     !  TODO : peut-on envisager un lecture par blocs pour limiter l'empreinte memoire?  
     !
     nZoneP = MatComposite(i)%Nzones_global
     allocate(array_int64_1d(nZoneP))
     allocate(pos_globalr(nZoneP))
     array_int64_1d = (/ (i_64,i_64=1,nZoneP) /)
     write(file_name,"(A,A)") trim(MatComposite(i)%dir),'/pos.bin'
     FU = 94480158
     call check_file(file_name,function_name)
     call getBinCoeffFromFile(file_name, array_int64_1d, nZoneP, pos_globalr, FU )
  
     ! On boucle sur les voxels composites du matériaux i
     !do j_64 = 1,size(MatComposite(i)%pos_globale)
     do j_64 = 1,size(pos_globalr)
        ! récupération de l'indice linéaire global du voxel composite
        !indice_pos = MatComposite(i)%pos_globale(j_64)
        indice_pos = int(pos_globalr(j_64),INT64)

        ! détermination des indices/coordonnées correspondants (divisions euclidiennes car type=int) (repère global X,Y,Z)

        ! Indice en Z 
        m = indice_pos/(nx*ny) ! m est le nombre de "tranches" de taille (nx*ny) sous le voxel considéré
                               ! le voxel se trouve donc dans la tranche m+1 (coordonnée Z)

        if (indice_pos - (m*nx*ny) == 0) then
           ! On se trouve sur le dernier voxel d'une tranche : le résultat de la division euclidienne donne cette fois 
           ! la tranche réelle du voxel, il faut donc lui retrancher 1 pour que m+1 corresponde à son indice en Z réel
           m = m-1
        end if

        ! on vérifie que le voxel se trouve dans le pinceau (coordonnée Z)
        if ( (m + 1 >= amx_start(3)) .and. (m + 1<= amx_end(3))) then 
           
           ! Indice en Y
           indice_pos = indice_pos - (m*nx*ny) ! premier reste
           l = indice_pos/nx ! l est le nombre de rangées situées avant celle du voxel dans la tranche m+1 
                             ! l'indice de la rangée où il se trouve est donc l+1 (coordonnée Y)
           
           if (indice_pos - (l*nx) == 0) then
              l = l-1 ! On se trouve sur le dernier voxel d'une rangée : le résultat de la division euclidienne donne cette fois 
                      ! la trancherangée réelle du voxel, il faut donc lui retrancher 1 pour que m+1 corresponde à son indice en Y réel
           end if
           
           ! on vérifie que le voxel se trouve dans le pinceau (coordonnée Y)
           if ( (l+1 >= amx_start(2)) .and. (l+1 <= amx_end(2))) then 

              k_64 = indice_pos - (l*nx) ! Indice en X : reste de la division euclidienne précédente

              ! on vérifie que le voxel se trouve dans le pinceau (coordonnée X)

              if ( (k_64 >= amx_start(1)) .and. (k_64 <= amx_end(1))) then 
                 ! VOXEL COMPOSITE DANS LE PINCEAU -> mise à jour de numM et numZ
                 numM(k_64,l+1,m+1) = i + nmateriauxhom ! numérotation des matériaux composite commence après le dernier matériau homogène (1+nmateriauxhom)
                 !numZ(k_64,l+1,m+1) = MatComposite(i)%Zone(j_64,1) ! indice de zone du voxel composite
                 numZ(k_64,l+1,m+1) = j_64                          ! indice de zone du voxel composite

              end if
           end if
        end if
     end do
     deallocate(array_int64_1d)
     deallocate(pos_globalr)
  end do

!==========================================================================================================================================
!
! SORTIES TEST : Reconstruction correcte des champs numM, numZ avec prise en compte des voxels composites
!                sortie dans le répertoire de sortie spécifié en ligne de commande pour vérification de la reconstruction correcte de 
!                la micro-structure

! SORTIE VTK TEST numM
! call write_header_vtk((trim(fic_vtk)//'_numM_reconstruit.vtk'),1000, nx,ny,nz,1._mytype,1._mytype,1._mytype)
! TEMPfield32=real(numM(amx_start(1):amx_end(1),amx_start(2):amx_end(2),amx_start(3):amx_end(3)))
! call write_bin(1,TEMPfield32,"numM_",(trim(fic_vtk)//'_numM_reconstruit.vtk'),1000)

! SORTIE VTK TEST numZ
! call write_header_vtk((trim(fic_vtk)//'_numZ_reconstruit.vtk'),1000, nx,ny,nz,1._mytype,1._mytype,1._mytype)
! TEMPfield32=real(numZ(amx_start(1):amx_end(1),amx_start(2):amx_end(2),amx_start(3):amx_end(3)))
! call write_bin(1,TEMPfield32,"numZ_",(trim(fic_vtk)//'_numZ_reconstruit.vtk'),1000)
!
  err_msg=trim(fic_vtk)   ! bidon - evite gcc-warning tant que la partie ci-dessus reste commentee
!===========================================================================================================================================


! On compte le nombre de voxels par materiau et par zone
  do k= amx_start(3),amx_end(3)
    do j= amx_start(2),amx_end(2)
      do i= amx_start(1),amx_end(1)
        ! 1ere colonne : compte du nombre de voxels pour chaque materiau
        tailles_mattotP(numM(i,j,k), 1) = tailles_mattotP(numM(i,j,k), 1)+ 1
        ! colonnes suivantes : compte du nb de voxels pour chaque zone du materiau
        tailles_mattotP(numM(i,j,k),numZ(i,j,k)+2) = tailles_mattotP(numM(i,j,k),numZ(i,j,k)+2)+1
      end do
    end do
  end do

  ! Comptage du nombre de matériaux et du nombre de zones de chaque matériau dans le pinceau
  j=0
  l=0
  do i=1,nmateriaux0%n

     Iph = .false.
     if (allocated(Interphases)) then 
        do k=1,size(Interphases)
           if ((Interphases(k)%numM == i) ) then
              Iph = .true.
              exit
           end if
        end do
     end if
     ! Iph repère si le matériau est un matériau
     ! d'interphase

     !  j: nombre de materiau dans le pinceau
     if (allocated(Interphases)) then
        ! Matériaux d'interphase présents
        if ((nrank == 0) .and. Iph) then
           ! si on est sur le processus 0, on rajoute les matériaux d'interphase dans le décompte
           ! (convention pour récupérer les données d'interphase lors
           ! l'initialisation des données composites)
           if ((tailles_mattotP(i,1) > 0) .and. Interphases(k)%all) then
              ! erreur =>  présence de voxels homogènes affectés du matériau d'interphase i
              write(err_msg,fmt="(A,I0,A,I0,A)") "Le materiau d'interphase ",i,&
                              "est affecté à ",tailles_mattotP(i,1),&
                              " voxels : le retirer de la liste d'interphases (read_geom_composite)"
              call amitex_abort(err_msg,2,0) 
           elseif ((tailles_mattotP(i,1) > 0) .and. ( .not. Interphases(k)%all)) then
                ! Le materiau est un materiau comprenant des zones interphases, et a des voxels homogenes
                ! presents sur le pinceau : on le prend en compte normalement
                j = j+1;
           else
                ! le materiau est un materiau d'interphase et n'a aucun voxel homogene present sur le pinceau
                ! on le prend en compte malgre tout sur le pinceau 0
                j = j+1;
                tailles_mattotP(i,1) = 1 ! mise à un pour initialisation du matériau dans la structure MattotP0 du pinceau 0 
                                        ! (voir ci-dessous)
                                        ! son nombre de zones est laisse a 0 pour le moment
           end if

        else if (tailles_mattotP(i,1) > 0) then
           ! si le matériau n'est pas un matériau d'interphase, on le prend en compte
               if (Iph) then
                  if (Interphases(k)%all) then
                  ! erreur => présence de voxels homogènes affectés du matériau d'interphase i
                  write(err_msg,fmt="(A,I0,A,I0,A)") "Le materiau d'interphase ",i,&
                                  "est affecté à ",tailles_mattotP(i,1),&
                                  " voxels : le retirer de la liste d'interphases (read_geom_composite)"
                  call amitex_abort(err_msg,2,0) 
                  end if
               end if
           j = j +1
        end if
     else if(tailles_mattotP(i,1)>0 ) then
        ! si le matériau est présent dans le pinceau, on le prend en compte
        j = j+1
     end if

     do k_64=1,nZone
        !  nombre de zone du materiau i dans le pinceau (stocke dans tailles_mattotP(i,2))
        if(tailles_mattotP(i,k_64+2)>0) tailles_mattotP(i,2)=tailles_mattotP(i,2)+1
     end do
  end do

! recapitulatif du comptage ci-dessus  
! tailles_mattotp(i,1) : nombre de voxels du materiau i dans le pinceau
! tailles_mattotp(i,2) : nombre de zones du materiau  i dans le pinceau
! tailles_mattotp(i,j+2) : nombre de voxels pour la zone j dans le pinceau

  ! taille de MattotP0: nombre de materiau dans le pinceau
  allocate(MattotP0(j),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_geom_composite)",2)
    nMat = 0 ! numero de materiau dans le pinceau

    ! determination de la presence du materiau dans le pinceau, du nombre de
    ! voxels associes, construction du tableau de position 'pos'  (ordonne par zone)
    do i=1,nmateriaux0%n

     Iph = .false.
     Iph_z = .false.
     if (allocated(Interphases)) then 
        do k=1,size(Interphases)
           if ((Interphases(k)%numM == i) .and. Interphases(k)%all) then
              Iph = .true. ! matériau d'interphase
              indice_iph = k ! indice du matériau dans la structure Interphases
              exit
           elseif ((Interphases(k)%numM == i) .and. .not.(Interphases(k)%all)) then
              Iph_z = .true. ! présence de zones d'interphase pour ce matériau
              indice_iph = k ! indice du matériau dans la structure Interphases
              exit
           end if
        end do
     end if  
 
      if(tailles_mattotp(i,1)>0 ) then
        ! le materiau i existe dans le pinceau:
      
        ! numero du materiau nMat (dans le pinceau) correspondant au materiau i (i numéro global)
        nMat= nMat+1    
        MattotP0(nMat)%numM=int(i,kind=4)
       
        if ((tailles_mattotp(i,2) == 0) .and. (Iph .or. Iph_z )) then
        ! matériau d'interphase (nombre de zones nul)  
        ! OU materiau ayant des zones interphases, et sans voxels sur le pinceau 0
        ! son tableau pos n'est pas alloué, Zones est alloué selon le nombre de zones demandees
        ! en entree
        ! et le champ Interphase de la structure MattotP0 est mis à ".true." (.false. par défaut)

           if (nrank /= 0) then ! vérification du respect de la convention : matériaux d'interphases dans le pinceau 0 uniquement
              write(err_msg,fmt="(3(A,I0),A)") "Matériau d'interphase ", MattotP0(j)%numM,&
                              " présent dans un pinceau de rang différent de 0 (read_geom_composite)"
              call amitex_abort(err_msg,1)
           end if

           MattotP0(nMat)%Interphase = Interphases(indice_iph)%all ! affectation du statut du materiau :
                                                         ! faux si le materiau a des zones interphases et des voxels homogenes
                                                         ! vrai si il est vrai materiau d'interphase (aucun voxel homogene dans la cellule) 
           
           ! allocation du tableau d'indice des zones
           ! Le nombre de zones de ce materiau est specifie dans la structure Interphases 
           allocate(MattotP0(nMat)%zone(Interphases(indice_iph)%Nzones,2),stat=alloc_stat)
           if(alloc_stat /=0) &
                        call amitex_abort("Espace memoire disponible insuffisant (read_geom_composite)",2)

           if (Iph_z) then
               ! materiau avec zones interphases
               do k_64 = 1,Interphases(indice_iph)%Nzones
                     ! Toutes les zones du materiau d'interphase seront stockees sur le pinceau 0. 
                     ! les deux indices de zone locaux et globaux sont donc les meme pour 
                     ! son tableau zone
                     ! ATTENTION : si le nombre de zones interphase est tres important, peu poser
                     !             pb d'avoir toutes les donnes sur un seul pinceau
                     MattotP0(nMat)%zone(k_64,1) = k_64
                     MattotP0(nMat)%zone(k_64,2) = Interphases(indice_iph)%zones(k_64)
               end do
              allocate(MattotP0(nMat)%Zones_interphase(Interphases(indice_iph)%Nzones),stat=alloc_stat)
              if(alloc_stat /=0) &
                     call amitex_abort("Espace memoire disponible insuffisant (read_geom_composite)",2)  
              MattotP0(nMat)%Zones_interphase = Interphases(indice_iph)%zones
               
           else
               ! materiau d'interphase      
               do k_64 = 1,Interphases(indice_iph)%Nzones
                     ! Toutes les zones du materiau d'interphase seront stockees sur le pinceau 0. 
                     ! les deux indices de zone locaux et globaux sont donc les meme pour 
                     ! son tableau zone
                     ! ATTENTION : si le nombre de zones interphase est tres important, peu poser
                     !             pb d'avoir toutes les donnes sur un seul pinceau
                     MattotP0(nMat)%zone(k_64,1) = k_64
                     MattotP0(nMat)%zone(k_64,2) = k_64
               end do
           end if

        else if (tailles_mattotp(i,2) > 0) then
  
           ! allocation du nombre de "position" (tableau pos)
           allocate(MattotP0(nMat)%pos(tailles_mattotP(i,1 )),stat=alloc_stat)
           if(alloc_stat /=0) &
                         call amitex_abort("Espace memoire disponible insuffisant (read_geom_composite)",2)

           ! allocation du tableau d'indice des zones
           if ((nrank == 0) .and. Iph_z) then
              ! sur le pinceau 0 on rajoute les zones interphases qui ne sont affectées à 
              ! aucun voxel dans le tableau Zone. Elles sont ajoutées en fin de tableau
              ! afin de pouvoir lire les coefficients matériau les caractérisant dans 
              ! read_mat_composite
              ! ATTENTION : si toutes les zones du matériau sont des interphases, 
              ! le test ne le détectera pas 
              ! solution : ajouter +1 à tailles_mattotP(i,2) ou passer le mat en 
              ! matériau d'interphase. (Suppose de savoir de combien de zones il se 
              ! constitue. 
              nzone_iph = size(Interphases(indice_iph)%zones)
              allocate(MattotP0(nMat)%zone(tailles_mattotp(i,2)+nzone_iph ,2),stat=alloc_stat)
              if(alloc_stat /=0) &
                   call amitex_abort("Espace memoire disponible insuffisant (read_geom_composite)",2)
           else
              allocate(MattotP0(nMat)%zone(tailles_mattotp(i,2) ,2),stat=alloc_stat)
              if(alloc_stat /=0) &
                   call amitex_abort("Espace memoire disponible insuffisant (read_geom_composite)",2)  
           end if
           ! mattotp2(nmat)%zone(i,1):nombre de voxels de la zone
           ! mattotp2(nmat)%zone(i,2):numero de zone (non parallelise)

           ! affectation des indices max des zones
           l=0
           m=0
           do k_64=3,size(tailles_mattotp(i,:),kind=8)  ! boucle sur les positions
              if(tailles_mattotp(i,k_64)>0)then          ! presence de la zone k-2 (materiau i) dans le pinceau
                 l=l+1                                 ! numero de la zone parallelise (pour le pinceau)
                 m= m + tailles_mattotp(i,k_64)           ! nombre de voxels cumules (de toutes les zones deja parcourues)
                 MattotP0(nMat)%zone(l,1)= m            ! indice maximal de la zone dans le tableau pos
                 MattotP0(nMat)%zone(l,2)= k_64-2          ! numero de zone (non parallelise)
              end if
           end do
           
           if  ((nrank == 0) .and. Iph_z) then
              do k_64=1,nzone_iph
                 l = l+1
                 MattotP0(nMat)%zone(l,1)= m+k_64       ! on incremente le dernier indice de %pos 
                                                    ! --> ces indices ne correspondent plus a des elements de %pos
                                                    ! --> ces indices correspondent a des colonnes des tableaux 
                                                    !     %VarInt et %VarInt0 sur le pinceau 0 destinees a stocker
                                                    !     les variables internes des zones interphases en attendant 
                                                    !     d'etre recuperees au sein des voxels composites
                                                    ! --> les zones d'interphases ne sont pas prises
                                                    ! en compte lors du calcul du comportement
                 MattotP0(nMat)%zone(l,2)= Interphases(indice_iph)%zones(k_64)
              end do
              allocate(MattotP0(nMat)%Zones_interphase(nzone_iph),stat=alloc_stat)
              if(alloc_stat /=0) &
                   call amitex_abort("Espace memoire disponible insuffisant (read_geom_composite)",2)  
              MattotP0(nMat)%Zones_interphase = Interphases(indice_iph)%zones
           end if

           allocate(ind_zones(tailles_mattotp(i,2)),stat=alloc_stat) ! ind_zone: indice courant de chaque zone
           if(alloc_stat /=0) &
                         call amitex_abort("Espace memoire disponible insuffisant (read_geom_composite)",2)
           ind_zones(1)=1
           if(tailles_mattotp(i,2)>1)then
              !ind_zones( zone )= indice max (zone-1) +1 (premier indice de zone)
              ind_zones(2:tailles_mattotp(i,2)) = MattotP0(nMat)%zone(1:tailles_mattotp(i,2)-1,1)+1
           end if

           ! numero de l'element a placer dans le tableau 'pos'
           nTot = 1
           ! affectation des positions
           ! ATTENTION CETTE BOUCLE EST TRES CHRONOPHAGE (si bcp de voxels composites)  
           ! C'est un copier-coller de read_geom  
           do j=amx_start(3),amx_end(3)
              do k=amx_start(2),amx_end(2)
                 do l=amx_start(1),amx_end(1)
                    if(numM(l,k,j)==i)then
                       nzone = numZ(l,k,j) ! numero de zone a trouver
                       m=1  ! indice de la zone actuelle
                       ! recherche de la zone correspondante (parallelise)
                       do while(MattotP0(nMat)%zone(m,2)/=numZ(l,k,j) )
                          m=m+1 
                       end do
                       ! position courante
                       MattotP0(nMat)%pos(ind_zones(m))=nTot
                       ! nouvel indice courant de zone
                       ind_zones(m)=ind_zones(m)+1
                    end if
                    nTot=nTot+1      
                 end do
              end do
           end do
           if(allocated(ind_zones)) deallocate(ind_zones)
        end if
     end if
  end do

  ! Affectations spécifiques : informations relatives à la géométrie des voxels composites
  !                            et à leur composition (dans le pinceau) à initialiser dans 
  !                            la structure MattotP0 à partir de la structure temporaire 
  !                            MatComposite (voir material_mod)

  nMat = 1
  do i=nmateriauxhom+1,nmateriaux0%n
  
     ! on evalue le nbre de phase du materiau
     nPhases = size(MatComposite(i-nmateriauxhom)%num_Phases)
     
     !si le materiau i existe dans le pinceau:
     if(tailles_mattotp(i,1)>0 ) then
        ! recherche du matériau (composite) i dans le tableau MattotP0 dans le pinceau
        do while(MattotP0(nMat)%numM /= i)
           nMat = nMat +1
        end do

        nPhases = size(MatComposite(i-nmateriauxhom)%num_Phases)

        ! affectation du nombre de phases du matériau
        MattotP0(nMat)%Nphase = nPhases

        ! allocation du tableau numM_composite (nPhases)
        allocate(MattotP0(nMat)%numM_composite(nPhases),stat=alloc_stat)
        if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_geom_composite)",2)

        ! affectation du tableau numM_composite
        MattotP0(nMat)%numM_composite(:) = MatComposite(i-nmateriauxhom)%num_Phases

        ! allocation du tableau numZ_composite (nZones,nPhases)
        allocate(MattotP0(nMat)%numZ_composite(tailles_mattotp(i,2),nPhases),stat=alloc_stat)
        if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_geom_composite)",2)

        ! affection du tableau numZ_composite
        ! MattotP0(nMat)%numZ_composite(:,:) = &
        !             MatComposite(i-nmateriauxhom)%Zone(MattotP0(nMat)%Zone(:,2),2:nPhases+1)
        NzoneP = size(MattotP0(nMat)%Zone(:,2))
        allocate(array_real_1d(NzoneP))
        do iphi = 1, Nphases
          write(file_name,"(A,A,I0,A)") trim(MatComposite(i-nmateriauxhom)%dir),'/zone',iphi,'.bin'
          FU = 35010797+iphi
          call check_file(file_name,function_name)
          call  getBinCoeffFromFile(file_name, MattotP0(nMat)%Zone(:,2), nZoneP, array_real_1d, FU )
          MattotP0(nMat)%numZ_composite(:,iphi) = int(array_real_1d(:),INT64)
        end do
        deallocate(array_real_1d)  
     else 
        !
        ! Appel a getBinCoeffFromFile avec des tableaux vides
        !
        allocate(array_int64_1d(0))
        allocate(array_real_1d(0))
        nZoneP=0_INT64
        do iphi = 1, Nphases
          write(file_name,"(A,A,I0,A)") trim(MatComposite(i-nmateriauxhom)%dir),'/zone',iphi,'.bin'
          FU = 35010797+iphi
          call check_file(file_name,function_name)
          call  getBinCoeffFromFile(file_name, array_int64_1d, nZoneP, array_real_1d, FU )
        end do
        deallocate(array_int64_1d)
        deallocate(array_real_1d)
     end if          
  end do

! verification de la continuite des nombres de materiaux/zones
    nzone = size(tailles_mattotp(1,:))-2
    allocate(ind_zones(nzone),stat=alloc_stat)  ! ind_zone(i): nombre de voxels pour la zone i, NON PARALLELISE (MPI_Reduce, MPI_SUM) 
    if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_geom_composite)",2)
    allocate(temp_tailles(nzone),stat=alloc_stat)
    if(alloc_stat /=0) call amitex_abort("Espace memoire disponible insuffisant (read_geom_composite)",2)
    
  do i=1,nmateriaux0%n

     Iph = .false.
     Iph_z = .false.
     if (allocated(Interphases)) then 
        do k=1,size(Interphases)
           if ((Interphases(k)%numM == i) .and. Interphases(k)%all) then
              Iph = .true.
              exit
           elseif ((Interphases(k)%numM == i) .and. .not.(Interphases(k)%all)) then
              Iph_z = .true.
              indice_iph = k ! présence de zones d'interphase pour ce matériau
              exit
           end if
        end do
     end if

    ! SOLUTION INITIALE
    ! nzone = size(tailles_mattotp(1,:))-2
    ! allocate(ind_zones(nzone),stat=alloc_stat)  ! ind_zone(i): nombre de voxels pour la zone i, NON PARALLELISE (MPI_Reduce, MPI_SUM) 
    ! call MPI_Allreduce(tailles_mattotp(i,3:nZone+2),ind_zones(:), nzone, MPI_INTEGER8, MPI_SUM, amitex_mpi_comm,ierror)
    !
    ! 2 PROBLEMES ICI : 1/ nzone en kind=8 (incompatible openMPI -> kind=4)
    !                   2/ le passage de 'tailles_mattotp(i,3:nZone+2)' en argument créé un problème
    !                      origine mal comprise...  possibilité : donnees non-contigues...
    !                      fonctionne en creant un tableau temporaire 1D (de donnees contigues...)
    ! 
    ! SOLUTION ACTUELLE pour compter les voxels en une seule somme parallele par materiau : 
    ! on passe par la fonction SUM_MPI_I8_LARGE qui permet de realiser des sommes paralleles avec un nombre 
    ! d'elements sueprieur a la valeur max des int_kind=4 (2^31 -1) 
    ! ATTENTION : pour les vecteurs depassant effectivement cette taille, la fonction n'est pas encore
    !            implementee (solution de secours pour le moment)
    ! ATTENTION : peut être passer l'allocation et la definition de temp_tailles dans SUM_MPI_I8_LARGE

    temp_tailles(:) = tailles_mattotp(i,3:nZone+2)
    call SUM_MPI_I8_LARGE(temp_tailles,ind_zones(:),nzone,ierror)    

    if(nrank==0) then   
       if(all(ind_zones .eq. 0) .and. (.not.Iph) ) then
          ! Aucun voxel détecté pour ce materiau et le matériau n'est pas un materiau d'interphase
          ! Il est possible que le materiau ne contienne que des zones interphases, cas
          ! non traite --> il faut passer le materiau en materiau d'interphase contenant plusieurs
          ! zones 
          write(err_msg,fmt="(A,I0,A)") "Materiau ",i," absent (lecture fichier vtk) (read_geom_composite)"
          call amitex_abort(err_msg,1,0) 
       else if (all(ind_zones .eq. 0) .and. Iph ) then
          ! confirmation qu'aucun voxel n'as été affecté au matériau d'interphase -> Info et non erreur
          write(err_msg,fmt="(A,I0,A)") "Aucun voxel affecté au matériau d'interphase (matériau ",i,&
               ") (lecture fichier vtk) (read_geom_composite)"
          call amitex_abort(err_msg,-1,0) 
       end if
       k_64=0
       do j_64=1,nzone ! size(ind_zones)
          if(ind_zones(j_64)==0 .and. k_64==0)then
             ! la zone portant le numero j est absente
             ! si une zone de numero > j est presente il s'agit d'une erreur
             k_64=j_64
             if (Iph_z) then
                if (any(Interphases(indice_iph)%zones == k_64)) then
                   ! zone k absente mais zone d'interphase : pas d'erreur
                   ! remise a 0 de k et reprise du test 
                   k_64 = 0
                end if
             end if
          else if(ind_zones(j_64)>0 .and. k_64>0 .and. (.not.Iph) ) then
        ! cas d'une zone absente
           ! zone et absente et pas une zone interphase : erreur
          if(j_64-1>k_64)then
            write (err_msg,fmt="(I0,A,I0)") k_64," a ",j_64-1
          else
            write (err_msg,fmt="(I0)") k_64
          end if
          write(err_msg,fmt="(A,I0,3A)") "Materiau ",i,", zone(s) ", trim(err_msg),&
                                         " absente(s) (lecture fichier vtk)  (read_geom_composite)"
          call amitex_abort(err_msg,1,0)
          k_64=0
        else if(ind_zones(j_64)<0) then
          write(err_msg, fmt="(A,I0,A)") "Subroutine read_geom, ind_zone(j) negatif,j= ",j_64,&
                                         " (read_geom_composite)"
          call amitex_abort(err_msg,1,0)
       end if
      end do
   end if
  end do

  
  if(allocated(ind_zones)) deallocate(ind_zones)

  if(allocated(tailles_MattotP)) deallocate(tailles_MattotP)

  if (allocated(temp_tailles)) deallocate(temp_tailles)

  ! arret si erreur (ci-dessus)
  call check_amitex_abort()

end subroutine read_geom_composite
!==================================================================================


