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
  Subroutine lu(M,Ipiv,Idet)
!  *                                           *
!  *********************************************
!  *********************************************

    Real (mytype), Intent (inout) :: M(:,:)
    Integer, Intent (out) :: Ipiv(:), Idet

    Integer :: I, J, Idim

    Idim = Size(M,1)
    ! First make pivoting
    CALL Pivoting(M, Ipiv, Idet)

    ! NOW: LU Decomposition
    ! The first step is done apart
    Do J = 2, Idim
       M(J, 1) = M(J, 1) / M(1, 1)
    End Do

    Do I = 2, Idim
       M(I, I) = M(I, I) - &
            & Dot_Product(M(I,1:I-1),M(1:I-1, I)) 

       Do J = i+1, Idim
          M(I, J) = M(I, J) - Dot_Product(M(I, 1:I-1),&
               & M(1:I-1, J))
          M(J, I) = M(J, I) - Dot_Product(M(J, 1:I-1)&
               &,M(1:I-1, I))

          If (Abs(M(I,I)) < Epsilon(1.0_mytype)) then
             if(nrank==0) write(OUTPUT_UNIT,*) "Erreur LU : M"
             if(nrank==0) write(OUTPUT_UNIT,'(6(F10.7))') M
             call amitex_abort("Erreur LU",2)
             !CALL Abort
            end if            
          M(J, I) = M(J, I) / M(I, I)
       End Do
    End Do

    Return
  End Subroutine lu

!  *********************************************
!> Solve a linear set of equations using LU 
!! decomposition. Both M and b are overwritten
!  *                                           *
  Subroutine LUSolve(M, b)
!  *                                           *
!  *********************************************
!  *********************************************

    Real (mytype), Intent (inout) :: M(:,:), b(:)

    Real (mytype) :: bcp(Size(M,1))
    Integer :: Ipiv(Size(M,1)), Id, Idim, I

    Idim = Size(M,1)
    if(Idim < 1) then
       !print*, "Matrice vide. Impossible d'inverser le systeme."
       !call abort
       call amitex_abort("Matrice vide. Impossible d'inverser le systeme.",2)
    end if
    CALL LU(M, Ipiv, Id)
    bcp = b
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
    Return
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

