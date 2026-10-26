!! (YC) codes copied from internet (see copyright below), for computing the eigen*

!! net code ****************************************************************************************
!! *************************************************************************************************
subroutine rs ( n, a, w, matz, z, ierr )

!*****************************************************************************80
!
!! RS computes eigenvalues and eigenvectors of real symmetric matrix.
!
!  Discussion:
!
!    RS calls the recommended sequence of EISPACK routines
!    to find the eigenvalues and eigenvectors (if desired)
!    of a real symmetric matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2018
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the real symmetric matrix.
!
!    Input, logical MATZ, is false if only eigenvalues are desired, 
!    and true if both eigenvalues and eigenvectors are desired.
!
!    Output, real ( kind = 8 ) W(N), the eigenvalues in ascending order.
!
!    Output, real ( kind = 8 ) Z(N,N), contains the eigenvectors, if MATZ
!    is true.
!
!    Output, integer ( kind = 4 ) IERR, is set equal to an error
!    completion code described in the documentation for TQLRAT and TQL2.
!    The normal completion code is zero.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) fv1(n)
  real ( kind = 8 ) fv2(n)
  integer ( kind = 4 ) ierr
  logical matz
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) z(n,n)

  if ( .not. matz ) then

    call tred1 ( n, a, w, fv1, fv2 )

    call tqlrat ( n, w, fv2, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'RS - Fatal error!'
      write ( *, '(a)' ) '  Error return from TQLRAT.'
      return
    end if

  else

    call tred2 ( n, a, w, fv1, z )

    call tql2 ( n, w, fv1, z, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'RS - Fatal error!'
      write ( *, '(a)' ) '  Error return from TQL2.'
      return
    end if

  end if

  return
end subroutine rs
!!

subroutine tred1 ( n, a, d, e, e2 )

!*****************************************************************************80
!
!! TRED1 transforms a real symmetric matrix to symmetric tridiagonal form.
!
!  Discussion:
!
!    TRED1 reduces a real symmetric matrix to a symmetric
!    tridiagonal matrix using orthogonal similarity transformations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 March 2018
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Martin, Reinsch, James Wilkinson,
!    TRED1,
!    Numerische Mathematik,
!    Volume 11, pages 181-195, 1968.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input/output, real ( kind = 8 ) A(N,N), on input, contains the real
!    symmetric matrix.  Only the lower triangle of the matrix need be supplied.
!    On output, A contains information about the orthogonal transformations
!    used in the reduction in its strict lower triangle.
!    The full upper triangle of A is unaltered.
!
!    Output, real ( kind = 8 ) D(N), contains the diagonal elements of the
!    tridiagonal matrix.
!
!    Output, real ( kind = 8 ) E(N), contains the subdiagonal elements of the
!    tridiagonal matrix in its last N-1 positions.  E(1) is set to zero.
!
!    Output, real ( kind = 8 ) E2(N), contains the squares of the corresponding
!    elements of E.  E2 may coincide with E if the squares are not needed.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) e2(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) scale

  d(1:n) = a(n,1:n)

  do i = 1, n
    a(n,i) = a(i,i)
  end do

  do i = n, 1, -1

    l = i - 1
    h = 0.0D+00
!
!  Scale row.
!
    scale = sum ( abs ( d(1:l) ) )

    if ( scale == 0.0D+00 ) then

      do j = 1, l
        d(j) = a(l,j)
        a(l,j) = a(i,j)
        a(i,j) = 0.0D+00
      end do

      e(i) = 0.0D+00
      e2(i) = 0.0D+00

      cycle

    end if

    d(1:l) = d(1:l) / scale

    do k = 1, l
      h = h + d(k) * d(k)
    end do

    e2(i) = h * scale * scale
    f = d(l)
    g = - sign ( sqrt ( h ), f )
    e(i) = scale * g
    h = h - f * g
    d(l) = f - g

    if ( 1 <= l ) then
!
!  Form A * U.
!
      e(1:l) = 0.0D+00

      do j = 1, l

        f = d(j)
        g = e(j) + a(j,j) * f
        do k = j + 1, l
          g = g + a(k,j) * d(k)
          e(k) = e(k) + a(k,j) * f
        end do

        e(j) = g

      end do
!
!  Form P.
!
      f = 0.0D+00

      do j = 1, l
        e(j) = e(j) / h
        f = f + e(j) * d(j)
      end do

      h = f / ( h + h )
!
!  Form Q.
!
      e(1:l) = e(1:l) - h * d(1:l)
!
!  Form reduced A.
!
      do j = 1, l

        f = d(j)
        g = e(j)

        a(j:l,j) = a(j:l,j) - f * e(j:l) - g * d(j:l)

      end do

    end if

    do j = 1, l
      f = d(j)
      d(j) = a(l,j)
      a(l,j) = a(i,j)
      a(i,j) = f * scale
    end do

  end do

  return
end subroutine tred1
!!
subroutine tred2 ( n, a, d, e, z )

!*****************************************************************************80
!
!! TRED2 transforms a real symmetric matrix to symmetric tridiagonal form.
!
!  Discussion:
!
!    TRED2 reduces a real symmetric matrix to a
!    symmetric tridiagonal matrix using and accumulating
!    orthogonal similarity transformations.
!
!    A and Z may coincide, in which case a single storage area is used
!    for the input of A and the output of Z.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 March 2018
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Martin, Reinsch, James Wilkinson,
!    TRED2,
!    Numerische Mathematik,
!    Volume 11, pages 181-195, 1968.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the real symmetric input matrix.  Only the
!    lower triangle of the matrix need be supplied.
!
!    Output, real ( kind = 8 ) D(N), the diagonal elements of the tridiagonal
!    matrix.
!
!    Output, real ( kind = 8 ) E(N), contains the subdiagonal elements of the
!    tridiagonal matrix in E(2:N).  E(1) is set to zero.
!
!    Output, real ( kind = 8 ) Z(N,N), the orthogonal transformation matrix
!    produced in the reduction.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  real ( kind = 8 ) hh
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) scale
  real ( kind = 8 ) z(n,n)

  do i = 1, n
    z(i:n,i) = a(i:n,i)
  end do

  d(1:n) = a(n,1:n)

  do i = n, 2, -1

    l = i - 1
    h = 0.0D+00
!
!  Scale row.
!
    scale = sum ( abs ( d(1:l) ) )

    if ( scale == 0.0D+00 ) then

      e(i) = d(l)

      do j = 1, l
        d(j) = z(l,j)
        z(i,j) = 0.0D+00
        z(j,i) = 0.0D+00
      end do

      d(i) = 0.0D+00

      cycle

    end if

    d(1:l) = d(1:l) / scale

    h = h + dot_product ( d(1:l), d(1:l) )

    f = d(l)
    g = - sign ( sqrt ( h ), f )
    e(i) = scale * g
    h = h - f * g
    d(l) = f - g
!
!  Form A*U.
!
    e(1:l) = 0.0D+00

    do j = 1, l

      f = d(j)
      z(j,i) = f
      g = e(j) + z(j,j) * f

      do k = j + 1, l
        g = g + z(k,j) * d(k)
        e(k) = e(k) + z(k,j) * f
      end do

      e(j) = g

    end do
!
!  Form P.
!
    e(1:l) = e(1:l) / h

    f = dot_product ( e(1:l), d(1:l) )

    hh = 0.5D+00 * f / h
!
!  Form Q.
!
    e(1:l) = e(1:l) - hh * d(1:l)
!
!  Form reduced A.
!
    do j = 1, l

      f = d(j)
      g = e(j)

      z(j:l,j) = z(j:l,j) - f * e(j:l) - g * d(j:l)

      d(j) = z(l,j)
      z(i,j) = 0.0D+00

    end do

    d(i) = h

  end do
!
!  Accumulation of transformation matrices.
!
  do i = 2, n

!   l = i - 1
    z(n,i-1) = z(i-1,i-1)
    z(i-1,i-1) = 1.0D+00
    h = d(i)

    if ( h /= 0.0D+00 ) then

      d(1:i-1) = z(1:i-1,i) / h

      do j = 1, i - 1

        g = dot_product ( z(1:i-1,i), z(1:i-1,j) )

        do k = 1, i - 1
          z(k,j) = z(k,j) - g * d(k)
        end do

      end do

    end if

    z(1:i-1,i) = 0.0D+00

  end do

  d(1:n) = z(n,1:n)

  z(n,1:n-1) = 0.0D+00
  z(n,n) = 1.0D+00

  e(1) = 0.0D+00

  return
end subroutine tred2
!!
subroutine tqlrat ( n, d, e2, ierr )

!*****************************************************************************80
!
!! TQLRAT computes all eigenvalues of a real symmetric tridiagonal matrix.
!
!  Discussion:
!
!    TQLRAT finds the eigenvalues of a symmetric
!    tridiagonal matrix by the rational QL method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 March 2018
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    C Reinsch,
!    Algorithm 464, TQLRAT,
!    Communications of the ACM,
!    Volume 16, page 689, 1973.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N).  On input, D contains the diagonal
!    elements of the matrix.  On output, D contains the eigenvalues in ascending
!    order.  If an error exit was made, then the eigenvalues are correct
!    in positions 1 through IERR-1, but may not be the smallest eigenvalues.
!
!    Input/output, real ( kind = 8 ) E2(N), contains in positions 2 through N 
!    the squares of the subdiagonal elements of the matrix.  E2(1) is
!    arbitrary.  On output, E2 has been overwritten by workspace
!    information.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, for no error,
!    J, if the J-th eigenvalue could not be determined after 30 iterations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e2(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) its
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( kind = 8 ) p
  !real ( kind = 8 ) pythag
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t

  ierr = 0

  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e2(i-1) = e2(i)
  end do

  f = 0.0D+00
  t = 0.0D+00
  e2(n) = 0.0D+00

  do l = 1, n

    its = 0
    h = abs ( d(l) ) + sqrt ( e2(l) )

    if ( t <= h ) then

      t = h
      b = abs ( t ) * epsilon ( b )
      c = b * b

    end if
!
!  Look for small squared sub-diagonal element.
!
    do m = l, n
      if ( e2(m) <= c ) then
        exit
      end if
    end do

    if ( m /= l ) then

      do

        if ( 30 <= its ) then
          ierr = l
          return
        end if

        its = its + 1
!
!  Form shift.
!
        l1 = l + 1
        s = sqrt ( e2(l) )
        g = d(l)
        p = ( d(l1) - g ) / ( 2.0D+00 * s )
        call pythag( p, 1.0D+00, r )
        d(l) = s / ( p + sign ( r, p ) )
        h = g - d(l)
        d(l1:n) = d(l1:n) - h
        f = f + h
!
!  Rational QL transformation.
!
        g = d(m)
        if ( g == 0.0D+00 ) then
          g = b
        end if

        h = g
        s = 0.0D+00
        mml = m - l

        do i = m - 1, l, -1
          p = g * h
          r = p + e2(i)
          e2(i+1) = s * r
          s = e2(i) / r
          d(i+1) = h + s * ( h + d(i) )
          g = d(i) - e2(i) / g
          if ( g == 0.0D+00 ) then
            g = b
          end if
          h = g * p / r
        end do

        e2(l) = s * g
        d(l) = h
!
!  Guard against underflow in convergence test.
!
        if ( h == 0.0D+00 ) then
          exit
        end if

        if ( abs ( e2(l) ) <= abs ( c / h ) ) then
          exit
        end if

        e2(l) = h * e2(l)

        if ( e2(l) == 0.0D+00 ) then
          exit
        end if

      end do

    end if

    p = d(l) + f
!
!  Order the eigenvalues.
!
    do i = l, 1, -1
      if ( i == 1 ) then
        d(i) = p
        exit
      else if ( d(i-1) <= p ) then
        d(i) = p
        exit
      end if
      d(i) = d(i-1)
    end do

  end do

  return
end subroutine tqlrat
!!

subroutine tql2 ( n, d, e, z, ierr )

!*****************************************************************************80
!
!! TQL2 computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
!
!  Discussion:
!
!    TQL2 finds the eigenvalues and eigenvectors of a symmetric
!    tridiagonal matrix by the QL method.  The eigenvectors of a full
!    symmetric matrix can also be found if TRED2 has been used to reduce this
!    full matrix to tridiagonal form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 March 2018
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Bowdler, Martin, Reinsch, James Wilkinson,
!    TQL2,
!    Numerische Mathematik,
!    Volume 11, pages 293-306, 1968.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N).  On input, the diagonal elements of
!    the matrix.  On output, the eigenvalues in ascending order.  If an error
!    exit is made, the eigenvalues are correct but unordered for indices
!    1,2,...,IERR-1.
!
!    Input/output, real ( kind = 8 ) E(N).  On input, E(2:N) contains the
!    subdiagonal elements of the input matrix, and E(1) is arbitrary.
!    On output, E has been destroyed.
!
!    Input, real ( kind = 8 ) Z(N,N).  On input, the transformation matrix
!    produced in the reduction by TRED2, if performed.  If the eigenvectors of
!    the tridiagonal matrix are desired, Z must contain the identity matrix.
!    On output, Z contains the orthonormal eigenvectors of the symmetric
!    tridiagonal (or full) matrix.  If an error exit is made, Z contains
!    the eigenvectors associated with the stored eigenvalues.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, normal return,
!    J, if the J-th eigenvalue has not been determined after
!    30 iterations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) dl1
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) el1
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) its
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( kind = 8 ) p
  !real ( kind = 8 ) pythag
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) s2
  real ( kind = 8 ) t
  real ( kind = 8 ) tst1
  real ( kind = 8 ) tst2
  real ( kind = 8 ) z(n,n)

  ierr = 0

  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e(i-1) = e(i)
  end do

  f = 0.0D+00
  tst1 = 0.0D+00
  e(n) = 0.0D+00

  do l = 1, n

    its = 0
    h = abs ( d(l) ) + abs ( e(l) )
    tst1 = max ( tst1, h )
!
!  Look for a small sub-diagonal element.
!
    do m = l, n
      tst2 = tst1 + abs ( e(m) )
      if ( tst2 == tst1 ) then
        exit
      end if
    end do

    if ( m /= l ) then

      do

        if ( 30 <= its ) then
          ierr = l
          return
        end if

        its = its + 1
!
!  Form shift.
!
        l1 = l + 1
        l2 = l1 + 1
        g = d(l)
        p = ( d(l1) - g ) / ( 2.0D+00 * e(l) )
        call pythag( p, 1.0D+00 , r )
        d(l) = e(l) / ( p + sign ( r, p ) )
        d(l1) = e(l) * ( p + sign ( r, p ) )
        dl1 = d(l1)
        h = g - d(l)
        d(l2:n) = d(l2:n) - h
        f = f + h
!
!  QL transformation.
!
        p = d(m)
        c = 1.0D+00
        c2 = c
        el1 = e(l1)
        s = 0.0D+00
        mml = m - l

        do i = m - 1, l, -1

          c3 = c2
          c2 = c
          s2 = s
          g = c * e(i)
          h = c * p
          call pythag( p, e(i) , r )
          e(i+1) = s * r
          s = e(i) / r
          c = p / r
          p = c * d(i) - s * g
          d(i+1) = h + s * ( c * g + s * d(i) )
!
!  Form vector.
!
          do k = 1, n
            h = z(k,i+1)
            z(k,i+1) = s * z(k,i) + c * h
            z(k,i) = c * z(k,i) - s * h
          end do

        end do

        p = - s * s2 * c3 * el1 * e(l) / dl1
        e(l) = s * p
        d(l) = c * p
        tst2 = tst1 + abs ( e(l) )

        if ( tst2 <= tst1 ) then
          exit
        end if

      end do

    end if

    d(l) = d(l) + f

  end do
!
!  Order eigenvalues and eigenvectors.
!
  do i = 1, n - 1

    k = i
    p = d(i)

    do j = i + 1, n

      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if

    end do

    if ( k /= i ) then

      d(k) = d(i)
      d(i) = p

      do j = 1, n
        t      = z(j,i)
        z(j,i) = z(j,k)
        z(j,k) = t
      end do

    end if

  end do

  return
end subroutine tql2
!!

subroutine pythag ( a, b, pythag1 )

!*****************************************************************************80
!
!! PYTHAG computes SQRT ( A * A + B * B ) carefully.
!
!  Discussion:
!
!    The formula
!
!      PYTHAG1 = sqrt ( A * A + B * B )
!
!    is reasonably accurate, but can fail if, for example, A^2 is larger
!    than the machine overflow.  The formula can lose most of its accuracy
!    if the sum of the squares is very large or very small.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the two legs of a right triangle.
!
!    Output, real ( kind = 8 ) PYTHAG1, the length of the hypotenuse.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) p
  real ( kind = 8 ) pythag1
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) u

  p = max ( abs ( a ), abs ( b ) )

  if ( p /= 0.0D+00 ) then

    r = ( min ( abs ( a ), abs ( b ) ) / p )**2

    do

      t = 4.0D+00 + r

      if ( t == 4.0D+00 ) then
        exit
      end if

      s = r / t
      u = 1.0D+00 + 2.0D+00 * s
      p = u * p
      r = ( s / u )**2 * r

    end do

  end if

  pythag1 = p

  return
end subroutine pythag
