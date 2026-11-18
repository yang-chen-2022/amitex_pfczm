!
!> MODULE to solve linear systems of equations
!!
!! Copyright (C) 2006  Alberto Ramos <alberto@martin.ft.uam.es>
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!! USA
!! 
!!
!! $ v. 1.0; Released: 15/09/2006; $

! ***********************************************************
! *
MODULE linear_mod
! *
! ***********************************************************

  use ISO_FORTRAN_ENV
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

#ifdef NEW_2DECOMP
  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi,       only : nrank
#else
  use decomp_2d, only : nrank,mytype 
#endif  
  use error_mod

  implicit none
 
  private
  
  public:: InvLU, LUSolve


CONTAINS

!  *********************************************
!>   Pivote M to arrange the elemnts in the
!!   diag. big.
!  *                                           *
  Subroutine Pivoting(M,Ipiv,Idet)
!  *                                           *
!  *********************************************

!  *********************************************

    Real (mytype), Intent (inout) :: M(:,:)
    Integer, Intent (out) :: Ipiv(:), Idet

    Real (mytype) :: Rval(Size(M,2))
    Integer :: Ipos(1), Kval, I

    ! Set  initial Idet, Ipiv
    Idet = 1
    forall (I=1:Size(Ipiv)) Ipiv(I) = I

    Do I = 1, Size(M,1) - 1
       Ipos = MaxLoc(Abs(M(I:,I))) + (I-1)

       Rval(:) = M(I,:)
       M(I,:) = M(Ipos(1),:)
       M(Ipos(1),:) = Rval(:)

       Kval = Ipiv(I)
       Ipiv(I) = Ipiv(Ipos(1))
       Ipiv(Ipos(1)) = Kval

       If (Ipiv(I) .ne. I) Idet = -Idet
    End Do

    Return
  End Subroutine Pivoting

!  *********************************************
!>   Makes LU decomposition of 
!! matrix M, with pivoting Ipiv. It returns  
!! in Idet if the number of premutations is 
!! even or odd.
!  *                                           *
  Subroutine lu(M,Ipiv,Idet,info)
!  *                                           *
!  *********************************************
!  *********************************************

    Real (mytype), Intent (inout) :: M(:,:)
    Integer, Intent (out) :: Ipiv(:), Idet, info

    Real (mytype) :: row(Size(M,2))
    Integer :: I, J, Idim, pivot, tmp

    Idim = Size(M,1)
    Idet = 1
    info = 0
    do I = 1, Idim
       Ipiv(I) = I
    end do
    if (.not. all(ieee_is_finite(M))) then
       info = -1
       return
    end if

    do I = 1, Idim
       pivot = maxloc(abs(M(I:Idim,I)), dim=1) + I - 1
       if (M(pivot,I) == 0._mytype) then
          info = I
          return
       end if
       if (pivot /= I) then
          row = M(I,:)
          M(I,:) = M(pivot,:)
          M(pivot,:) = row
          tmp = Ipiv(I)
          Ipiv(I) = Ipiv(pivot)
          Ipiv(pivot) = tmp
          Idet = -Idet
       end if
       do J = I+1, Idim
          M(J,I) = M(J,I) / M(I,I)
          M(J,I+1:Idim) = M(J,I+1:Idim) - M(J,I)*M(I,I+1:Idim)
       end do
       if (.not. all(ieee_is_finite(M))) then
          info = -1
          return
       end if
    end do
  End Subroutine lu

!  *********************************************
!> Solve a linear set of equations using LU 
!! decomposition. Both M and b are overwritten on success.
!! Optional info returns zero on success, the zero pivot index on a singular
!! matrix, or -1 for non-finite data. On failure M and b are restored.
!  *                                           *
  Subroutine LUSolve(M, b, info)
!  *                                           *
!  *********************************************
!  *********************************************

    Real (mytype), Intent (inout) :: M(:,:), b(:)
    Integer, Intent (out), optional :: info

    Real (mytype) :: bcp(Size(M,1)), original(Size(M,1),Size(M,2))
    Integer :: Ipiv(Size(M,1)), Id, Idim, I, status

    Idim = Size(M,1)
    if(Idim < 1) then
       !print*, "Matrice vide. Impossible d'inverser le systeme."
       !call abort
       call amitex_abort("Matrice vide. Impossible d'inverser le systeme.",2)
    end if
    if(size(M,2) /= Idim .or. size(b) /= Idim) &
       call amitex_abort("LUSolve: incompatible matrix or vector dimensions",2)
    original = M
    bcp = b
    status = 0
    if (.not. all(ieee_is_finite(b))) then
       status = -1
    else
       CALL LU(M, Ipiv, Id, status)
    end if
    if (status == 0) then
       ! Now we Permutate b
       Do I = 1, Idim
          b(I) = bcp(Ipiv(I))
       End Do
       ! First solve Lx = b
       Do I = 2, Idim
          b(I) = b(I) - Dot_Product(M(I,1:I-1),b(1:I-1))
       End Do

       ! Now solve Ux = b
       b(Idim) = b(Idim) / M(Idim, Idim)
       Do I = Idim - 1, 1, -1
          b(I) = b(I) - Dot_Product(M(I, I+1:Idim),b(I+1:Idim))
          b(I) = b(I) / M(I, I)
       End Do
       if (.not. all(ieee_is_finite(b))) status = -1
    end if

    if (status /= 0) then
       M = original
       b = bcp
       if (.not. present(info)) then
          write(OUTPUT_UNIT,'(A,I0,A,I0)') "LU failure: rank ", nrank, ", pivot/status ", status
          write(OUTPUT_UNIT,'(A)') "Original matrix (rows):"
          do I = 1, Idim
             write(OUTPUT_UNIT,'(*(ES24.16E3,1X))') original(I,:)
          end do
          write(OUTPUT_UNIT,'(A)') "Right-hand side:"
          write(OUTPUT_UNIT,'(*(ES24.16E3,1X))') bcp
          flush(OUTPUT_UNIT)
          call amitex_abort("Erreur LU",2)
       end if
    end if
    if (present(info)) info = status
  End Subroutine LUSolve

!  *********************************************
!> Compute the determinant of the matrix M.
!  *                                           *
  Real (mytype) Function Det(M)
!  *                                           *
!  *********************************************
!  *********************************************

    Real (mytype), Intent (in) :: M(:,:)
    
    Real (mytype) :: Mcp(Size(M,1),Size(M,1))
    Integer :: Ipiv(Size(M,1)), Id, I, J, Idim
 
    Idim = Size(M,1)
    Mcp = M
    ! First make pivoting
    CALL Pivoting(Mcp, Ipiv, Id)

    
    ! NOW: LU Decomposition
    ! Separamos el paso I = 1
    Do J = 2, Idim
       Mcp(J, 1) = Mcp(J, 1) / Mcp(1, 1)
    End Do

    Det = Mcp(1,1)
    Do I = 2, Idim
       Mcp(I, I) = Mcp(I, I) - &
            & Dot_Product(Mcp(I,1:I-1), Mcp(1:I-1, I)) 

       Do J = i+1, Idim
          Mcp(I, J) = Mcp(I, J) - Dot_Product(Mcp(I, 1:I-1),&
               & Mcp(1:I-1, J))
          Mcp(J, I) = Mcp(J, I) - Dot_Product(Mcp(J, 1:I-1)&
               &,Mcp(1:I-1, I))

          If (Abs(Mcp(I,I)) < Epsilon(1.0_mytype)) Then
             Det = 0.0_mytype
             Return
          End If
          Mcp(J, I) = Mcp(J, I) / Mcp(I, I)
       End Do
       Det = Det * Mcp(I,I)
    End Do

    Det = Real(Id,mytype)*Det

    Return
  End Function Det

!  *********************************************
!> Compute the inverse of the matrix M 
!! using LU decomposition
!  *                                           *
  Subroutine InvLU(M,Inv)
!  *                                           *
!  *********************************************

!  *********************************************

    Real (mytype), Intent (in) :: M(:,:)
    Real (mytype), Intent (out) :: Inv(size(M,1),size(M,1))
    
    Real (mytype) :: one(size(M,1)), Mcp(size(M,1),size(M,1))
    Integer :: Idim, i
    Inv = M

    Idim = size(M,1)
    if(Idim /= size(M,2)) then
       !print*, "Matrice non carree. Impossible de l'inverser."
       !call abort
       call amitex_abort("Matrice non carree. Impossible de l'inverser.",2)
    end if
    if(Idim < 1) then
       !print*, "Matrice vide. Impossible de l'inverser."
       !call abort
       call amitex_abort("Matrice vide. Impossible de l'inverser.",2)
    end if
    do i = 1,Idim
       one=0._mytype
       one(i) = 1._mytype
       Mcp = M
       call LUSolve(Mcp, one)
       Inv(:,i)=one
    end do
    
  End Subroutine InvLU

End MODULE linear_mod

