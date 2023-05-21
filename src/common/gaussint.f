      subroutine gaussq(kind, n, alpha, beta, kpts, endpts, b, t, w)
c  *********************************************************************
c  *                                                                   *
c  *  This set of routines computes the nodes t(j) and weights w(j)    *
c  *  for gaussian-type quadrature rules with pre-assigned nodes.      *
c  *  These are used when one wishes to approximate:                   *
c  *                                                                   *
c  *              integral (from a to b)  f(x) w(x) dx                 *
c  *                                                                   *
c  *                           n                                       *
c  *     by                   sum w  f(t )                             *
c  *                          j=1  j    j                              *
c  *                                                                   *
c  *  Note that w(x) and w(j) have no connection with each other.      *
c  *  here, w(x) is one of six possible non-negative weight functions  *
c  *  (listed below), and f(x) is the function to be integrated.       *
c  *  Gaussian quadrature is particularly useful on infinite intervals *
c  *  (with appropriate weight functions), since then other techniques *
c  *  often fail.                                                      *
c  *  Associated with each weight function w(x) is a set of orthogonal *
c  *  polynomials. The nodes t(j) are just the zeroes of the proper    *
c  *  n-th degree polynomial. Underflows may sometimes occur, but they *
c  *  are harmless.                                                    *
c  *                                                                   *
c  *  References                                                       *
c  *  ==========                                                       *
c  *                                                                   *
c  *   1.  golub, g. h., and welsch, j. h., "calculation of gaussian   *
c  *       quadrature rules," mathematics of computation 23 (april,    *
c  *       1969), pp. 221-230.                                         *
c  *   2.  golub, g. h., "some modified matrix eigenvalue problems,"   *
c  *       siam review 15 (april, 1973), pp. 318-334 (section 7).      *
c  *   3.  stroud and secrest, gaussian quadrature formulas, prentice- *
c  *       hall, englewood cliffs, n.j., 1966.                         *
c  *                                                                   *
c  *  Arguments                                                        *
c  *  =========                                                        *
c  *                                                                   *
c  *  kind      integer scalar;                                        *
c  *            on entry, a control parameter which takes values       *
c  *            between 1 and 6 and gives the type of quadrature       *
c  *            rule. Unchanged on exit.                               *
c  *                                                                   *
c  *     kind = 1:  legendre quadrature, w(x) = 1 on (-1, 1)           *
c  *     kind = 2:  chebyshev quadrature of the first kind             *
c  *                w(x) = 1/sqrt(1 - x*x) on (-1, +1)                 *
c  *     kind = 3:  chebyshev quadrature of the second kind            *
c  *                w(x) = sqrt(1 - x*x) on (-1, 1)                    *
c  *     kind = 4:  hermite quadrature, w(x) = exp(-x*x) on            *
c  *                (-infinity, +infinity)                             *
c  *     kind = 5:  jacobi quadrature, w(x) = (1-x)**alpha * (1+x)**   *
c  *                beta on (-1, 1), alpha, beta .gt. -1.              *
c  *                note: kind=2 and 3 are a special case of this.     *
c  *     kind = 6:  generalized laguerre quadrature, w(x) = exp(-x)*   *
c  *                x**alpha on (0, +infinity), alpha .gt. -1          *
c  *                                                                   *
c  *  n         integer scalar;                                        *
c  *            on entry, gives the number of points used for the      *
c  *            quadrature rule. Unchanged on exit.                    *
c  *                                                                   *
c  *  alpha     real*8 scalar;                                         *
c  *            on entry, parameter used only for gauss-jacobi and     *
c  *            gauss-laguerre quadrature (otherwise use 0.0d0). alpha *
c  *            enters the definition of the weight function (see      *
c  *            above). Unchanged on exit.                             *
c  *                                                                   *
c  *  beta      real*8 scalar;                                         *
c  *            on entry, parameter used only for gauss-jacobi         *
c  *            quadrature (otherwise use 0.0d0). beta enters the      *
c  *            the definition of the weight function (see above).     *
c  *            Unchanged on exit.                                     *
c  *                                                                   *
c  *  kpts      integer scalar;                                        *
c  *            on entry, parameter equal to 0 unless left or right    *
c  *            end-points (or both) of the interval are required to   *
c  *            be a node (this is called gauss-radau or gauss-lobatto *
c  *            quadrature). Then kpts is the number of fixed end-     *
c  *            points (1 or 2). Unchanged on exit.                    *
c  *                                                                   *
c  *  endpts    real*8 (2 x 1) array;                                  *
c  *            on entry, contains the values of any fixed endpoints   *
c  *            if kpts=1 or 2. Unchanged on exit.                     *
c  *                                                                   *
c  *  b         real*8 (1000 x 1) array;                               *
c  *            scratch array.                                         *
c  *                                                                   *
c  *  t         real*8 (1000 x 1) array;                               *
c  *            on exit, contains the desired nodes t(j).              *
c  *                                                                   *
c  *  w         real*8 (1000 x 1) array;                               *
c  *            on exit, contains the desired weights w(j).            *
c  *                                                                   *
c  *  Auxiliary routines                                               *
c  *  ==================                                               *
c  *                                                                   *
c  *  class         subroutine                                         *
c  *                                                                   *
c  *  solve         subroutine                                         *
c  *                                                                   *
c  *********************************************************************
      implicit none

c  ------ SCALAR ARGUMENTS ------  

      integer        kind, n, kpts 

      real*8         alpha, beta 

c  ------ ARRAY ARGUMENTS ------

      real*8         t(1000), w(1000), b(1000), endpts(2)

c  ------ LOCAL SCALARS ------

      integer        i, ierr      

      real*8         t1, gam, muzero    

c  ------ AUXILIARY FUNCTIONS ------

      real*8         solve

 
      call class (kind, n, alpha, beta, b, t, muzero)
 
c *** The matrix of coefficients is assumed to be symmetric.
c *** the array t contains the diagonal elements, the array
c *** b the off-diagonal elements.
c *** Makes appropriate changes in the lower right 2 by 2
c *** submatrix.
c
      if (kpts.eq.0)  go to 100
      if (kpts.eq.2)  go to  50

c *** If kpts=1, only t(n) must be changed.
 
      t(n) = solve(endpts(1), n, t, b)*b(n-1)**2 + endpts(1)
      go to 100
 
c *** If kpts=2, t(n) and b(n-1) must be recomputed
 
   50 gam = solve(endpts(1), n, t, b)
      t1 = ((endpts(1) - endpts(2))/(solve(endpts(2), n, t, b) - gam))
      b(n-1) = dsqrt(t1)
      t(n) = endpts(1) + gam*t1
 
c *** Note that the indices of the elements of b run from 1 to n-1
c *** and thus the value of b(n) is arbitrary.
c *** Now computes the eigenvalues of the symmetric tridiagonal
c *** matrix, which has been modified as necessary.
c *** The method used is a ql-type method with origin shifting
 
  100 w(1) = 1.0d0
      do 105 i = 2, n
         w(i) = 0.0d0
  105 continue  
      call gausq2 (n, t, b, w, ierr)
      do 110 i = 1, n
         w(i) = muzero * w(i) * w(i)
  110 continue  
 
      return
      end




      real*8 function solve(shift, n, a, b)

c  *********************************************************************
c  *                                                                   *
c  *  This procedure performs elimination to solve for the n-th        *
c  *  component of the solution delta to the equation                  *
c  *                                                                   *
c  *            (jn - shift*identity) * delta  = en,                   *
c  *                                                                   *
c  *  where en is the vector of all zeroes except for 1 in the n-th    *
c  *  position. The matrix jn is symmetric tridiagonal, with diagonal  *
c  *  elements a(i), off-diagonal elements b(i).  this equation must   *
c  *  be solved to obtain the appropriate changes in the lower 2 by 2  *
c  *  submatrix of coefficients for orthogonal polynomials.            *
c  *                                                                   *
c  *  Arguments                                                        *
c  *  =========                                                        *
c  *                                                                   *
c  *  shift     real*8 scalar;                                         *
c  *            on entry, the coefficient shift in the above equation. *
c  *            unchanged on exit.                                     *
c  *                                                                   *
c  *  n         integer scalar;                                        *
c  *            on entry, actual dimension of vectors a and b.         *
c  *            unchanged on exit.                                     *
c  *                                                                   *
c  *  a         real*8 (1 x 1000) array;                               *
c  *            on entry, contains the diagonal elements of matrix jn. *
c  *            Unchanged on exit.                                     *
c  *                                                                   *
c  *  b         real*8 (1 x 1000) array;                               *
c  *            on entry, contains the offdiagonal elements of matrix  *
c  *            jn. Unchanged on exit.                                 *
c  *                                                                   *
c  *********************************************************************
 
      implicit none

c  ------ SCALAR ARGUMENTS ------

      integer        n

      real*8         shift

c  ------ ARRAY ARGUMENTS ------

      real*8         a(1000), b(1000)

c  ------ LOCAL SCALARS ------

      integer        i, nm1

      real*8         alpha


      alpha = a(1) - shift
      nm1 = n - 1
      do 10 i = 2, nm1
         alpha = a(i) - shift - b(i-1)**2/alpha
   10 continue 
      solve = 1.0d0/alpha
      return
      end




      subroutine class(kind, n, alpha, beta, b, a, muzero)
c
c           this procedure supplies the coefficients a(j), b(j) of the
c        recurrence relation
c
c             b p (x) = (x - a ) p   (x) - b   p   (x)
c              j j            j   j-1       j-1 j-2
c
c        for the various classical (normalized) orthogonal polynomials,
c        and the zero-th moment
c
c             muzero = integral w(x) dx
c
c        of the given polynomial's weight function w(x).  since the
c        polynomials are orthonormalized, the tridiagonal matrix is
c        guaranteed to be symmetric.
c
c           the input parameter alpha is used only for laguerre and
c        jacobi polynomials, and the parameter beta is used only for
c        jacobi polynomials.  the laguerre and jacobi polynomials
c        require the gamma function.
c
c        WARNING: the gamma function provided here might be no accurate 
c                 enough

      double precision a(1000), b(1000), muzero, alpha, beta
      double precision abi, a2b2, pi, dsqrt, ab
c      double precision gamma  
c
      pi = 4.0d0 * datan(1.0d0)
      nm1 = n - 1
      go to (10, 20, 30, 40, 50, 60), kind
c
c              kind = 1:  legendre polynomials p(x)
c              on (-1, +1), w(x) = 1.
c
   10 muzero = 2.0d0
      do 11 i = 1, nm1
         a(i) = 0.0d0
         abi = i
   11    b(i) = abi/dsqrt(4*abi*abi - 1.0d0)
      a(n) = 0.0d0
      return
c
c              kind = 2:  chebyshev polynomials of the first kind t(x)
c              on (-1, +1), w(x) = 1 / sqrt(1 - x*x)
c
   20 muzero = pi
      do 21 i = 1, nm1
         a(i) = 0.0d0
   21    b(i) = 0.5d0
      b(1) = dsqrt(0.5d0)
      a(n) = 0.0d0
      return
c
c              kind = 3:  chebyshev polynomials of the second kind u(x)
c              on (-1, +1), w(x) = sqrt(1 - x*x)
c
   30 muzero = pi/2.0d0
      do 31 i = 1, nm1
         a(i) = 0.0d0
   31    b(i) = 0.5d0
      a(n) = 0.0d0
      return
c
c              kind = 4:  hermite polynomials h(x) on (-infinity,
c              +infinity), w(x) = exp(-x**2)
c
   40 muzero = dsqrt(pi)
      do 41 i = 1, nm1
         a(i) = 0.0d0
   41    b(i) = dsqrt(i/2.0d0)
      a(n) = 0.0d0
      return
c
c              kind = 5:  jacobi polynomials p(alpha, beta)(x) on
c              (-1, +1), w(x) = (1-x)**alpha + (1+x)**beta, alpha and
c              beta greater than -1
c
   50 ab = alpha + beta
      abi = 2.0d0 + ab
c      muzero = 2.0d0 ** (ab + 1.0d0) * gamma(alpha + 1.0d0) * gamma(
c     x beta + 1.0d0) / gamma(abi)
      a(1) = (beta - alpha)/abi
      b(1) = dsqrt(4.0d0*(1.0d0 + alpha)*(1.0d0 + beta)/((abi + 1.0d0)*
     1  abi*abi))
      a2b2 = beta*beta - alpha*alpha
      do 51 i = 2, nm1
         abi = 2.0d0*i + ab
         a(i) = a2b2/((abi - 2.0d0)*abi)
   51    b(i) = dsqrt (4.0d0*i*(i + alpha)*(i + beta)*(i + ab)/
     1   ((abi*abi - 1)*abi*abi))
      abi = 2.0d0*n + ab
      a(n) = a2b2/((abi - 2.0d0)*abi)
      return
c
c              kind = 6:  laguerre polynomials l(alpha)(x) on
c              (0, +infinity), w(x) = exp(-x) * x**alpha, alpha greater
c              than -1.
c           WARNING: we assumes here that alpha=0 --> otherwise, the 
c                    next line should be reactivated.

   60 muzero = 0d0 !gamma(alpha + 1.0d0)
      do 61 i = 1, nm1
         a(i) = 2.0d0*i - 1.0d0 + alpha
   61    b(i) = dsqrt(i*(i + alpha))
      a(n) = 2.0d0*n - 1 + alpha
      return
      end





      subroutine gausq2(n, d, e, z, ierr)
c
c     this subroutine is a translation of an algol procedure,
c     num. math. 12, 377-383(1968) by martin and wilkinson,
c     as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c     this is a modified version of the 'eispack' routine imtql2.
c
c     this subroutine finds the eigenvalues and first components of the
c     eigenvectors of a symmetric tridiagonal matrix by the implicit ql
c     method.
c
c     on input:
c
c        n is the order of the matrix;
c
c        d contains the diagonal elements of the input matrix;
c
c        e contains the subdiagonal elements of the input matrix
c          in its first n-1 positions.  e(n) is arbitrary;
c
c        z contains the first row of the identity matrix.
c
c      on output:
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1, 2, ..., ierr-1;
c
c        e has been destroyed;
c
c        z contains the first components of the orthonormal eigenvectors
c          of the symmetric tridiagonal matrix.  if an error exit is
c          made, z contains the eigenvectors associated with the stored
c          eigenvalues;
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     ------------------------------------------------------------------
c
      integer i, j, k, l, m, n, ii, mml, ierr
      real*8 d(1000), e(1000), z(1000), b, c, f, g, p, r, s, machep
      real*8 dsqrt, dabs, dsign
c
      machep=1.0d-15
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      e(n) = 0.0d0
      do 240 l = 1, n
         j = 0
c     :::::::::: look for small sub-diagonal element ::::::::::
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            if (dabs(e(m)) .le. machep * (dabs(d(m)) + dabs(d(m+1))))
     x         go to 120
  110    continue
c
  120    p = d(l)
         if (m .eq. l) go to 240
         if (j .eq. 30) go to 1000
         j = j + 1
c     :::::::::: form shift ::::::::::
         g = (d(l+1) - p) / (2.0d0 * e(l))
         r = dsqrt(g*g+1.0d0)
         g = d(m) - p + e(l) / (g + dsign(r, g))
         s = 1.0d0
         c = 1.0d0
         p = 0.0d0
         mml = m - l
c
c     :::::::::: for i=m-1 step -1 until l do -- ::::::::::
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            if (dabs(f) .lt. dabs(g)) go to 150
            c = g / f
            r = dsqrt(c*c+1.0d0)
            e(i+1) = f * r
            s = 1.0d0 / r
            c = c * s
            go to 160
  150       s = f / g
            r = dsqrt(s*s+1.0d0)
            e(i+1) = g * r
            c = 1.0d0 / r
            s = s * c
  160       g = d(i+1) - p
            r = (d(i) - g) * s + 2.0d0 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
c     :::::::::: form first component of vector ::::::::::
            f = z(i+1)
            z(i+1) = s * z(i) + c * f
  200       z(i) = c * z(i) - s * f
c
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0d0
         go to 105
  240 continue
c
c     :::::::::: order eigenvalues and eigenvectors ::::::::::
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
         p = z(i)
         z(i) = z(k)
         z(k) = p
  300 continue
c
      go to 1001
c     :::::::::: set error -- no convergence to an
c                eigenvalue after 30 iterations ::::::::::
 1000 ierr = l
 1001 return
c     :::::::::: last card of gausq2 ::::::::::
      end
c
C-----------------------------------------------------------------------
      SUBROUTINE MATPRT (AA,NL,NC,NLMX,ITYPE,IFORM,ANAME)
C#######################################################################
C#    MATPRT IS THE DRIVER FOR LPRINT                                  #
C#    PRINTS OUT A MATRIX OF INTEGER,REAL OR COMPLEX NUMBERS  A(M,N).  #
C#---------------------------------------------------------------------#
C#    AA     : MATRIX TO BE PRINTED                                    #
C#    NL     : NUMBER OF LINES TO BE PRINTED                           #
C#    NC     : NUMBER OF COLUMNS TO BE PRINTED                         #
C#    NLMX   : ROW DIMENSION OF AA IN THE CALLING PROGRAM              #
C#    ITYPE  : DESCRIBES THE TYPE OF THE NUMBERS:                      #
C#             4 FOR SINGLE PRECISION OR INTEGERS                      #
C#             8 FOR DOUBLE PRECISION OR INTEGERS                      #
C#             24 (28) FOR SINGLE (DOUBLE) PRECISION COMPLEX           #
C#             ( A COMPLEX MATRIX IS SEPARATED INTO TWO REAL MATRICES  #
C#               OF REAL AND IMAGINARY PARTS)                          #
C#    IFORM  : DESCRIBES THE FORMAT CODE;                              #
C#             1 FOR F FORMAT CODE, 2 FOR E FORMAT CODE,               #
C#             3 FOR G FORMAT CODE, 4 FOR I FORMAT CODE,               #
C#             ADD 10 TO OBTAIN HIGH ACCURACY PRINTOUT                 #
C#    NAME   : 8 CHARACTERS TO SPECIFY THE NAME OF THE MATRIX ON LIST. #
C#---------------------------------------------------------------------#
C#    J.M.L. (1980); MAJOR REVISION 10/1982 (CRAY VERSION)             #
C#######################################################################
      CHARACTER*8 ANAME,ANAMP,BLANK
      CHARACTER*30 ANC(2)
      LOGICAL COMPLX,CRAY
      DOUBLE PRECISION AA8
      DIMENSION AA(*),AA4(20),AA8(10)
      EQUIVALENCE (AA4(1),AA8(1))
      DATA BLANK / '        '/
     &    ,ANC   / ' REAL PART OF THE MATRIX      '
     &           , ' IMAGINARY PART OF THE MATRIX ' /
C
      ISKIP = 0
      GO TO 10
C
C     ---------------------------------------------------
      ENTRY MATPR1 (AA,NL,NC,NLMX,ITYPE,IFORM,ANAME,ISKP)
C     ---------------------------------------------------
      ISKIP = ISKP
C
   10 IF (NC .LE. 0) THEN
         PRINT 9002,ANAME,NC
         RETURN
      ENDIF
C
      IF (NL*NC .GT. 10000) GO TO 8000
      CALL MCHINE (CRAY)
C
      IF (ISKIP.GE.1) WRITE (6,9200)
C
      LWORD = ITYPE
      IF (IFORM .EQ. 4 .OR. IFORM .EQ. 14) GO TO 4000
C
      IF (IFORM .LT. 10) NCMX = 10
      IF (IFORM .GT. 10) NCMX =  5
      COMPLX = ITYPE .EQ. 24 .OR. ITYPE .EQ. 28
      IF (COMPLX) GO TO 200
      IF (CRAY) GO TO 1000
      IF (LWORD .EQ. 4) GO TO 1000
      IF (LWORD .EQ. 8) GO TO 1100
C
 200  LWORD = ITYPE-20
      IF (CRAY) GO TO 2000
      IF (LWORD .EQ. 4) GO TO 2000
      IF (LWORD .EQ. 8) GO TO 2100
      PRINT 9000,ITYPE
      RETURN
C
C --- REAL NUMBERS (SINGLE PRECISION)
C
 1000 JS = -1
      NCLEFT = NC
 1010 NCP = MIN0 (NCLEFT,NCMX)
      ANAMP = ANAME
      IF (NL .GT. 20) PRINT 9001
      DO 1011 I = 1,NL
         DO 1012 J = 1,NCP
            IJ = I+(JS+J)*NLMX
            AA8(J) = AA(IJ)
 1012    CONTINUE
         CALL LPRINT (AA8,AA4,I,NCP,IFORM,ANAMP)
         ANAMP = BLANK
 1011 CONTINUE
      NCLEFT = NCLEFT-NCMX
      JS = JS+NCMX
      IF (NCLEFT .LE. 0) RETURN
      GO TO 1010
C
C --- REAL NUMBERS (DOUBLE PRECISION)
C
 1100 JS = -1
      NCLEFT = NC
 1110 NCP = MIN0 (NCLEFT,NCMX)
      ANAMP = ANAME
      IF (NL .GT. 20) PRINT 9001
      DO 1111 I = 1,NL
         DO 1112 J = 1,NCP
            IJ = I+(JS+J)*NLMX
            AA4(2*J-1) = AA(2*IJ-1)
            AA4(2*J)   = AA(2*IJ  )
 1112    CONTINUE
         CALL LPRINT (AA8,AA4,I,NCP,IFORM,ANAMP)
         ANAMP = BLANK
 1111 CONTINUE
      NCLEFT = NCLEFT-NCMX
      JS = JS+NCMX
      IF (NCLEFT .LE. 0) RETURN
      GO TO 1110
C
C --- COMPLEX NUMBERS (SINGLE PRECISION)
C
 2000 DO 2010 IP = 1,2
         PRINT 9010,ANC(IP),ANAME
         JS = -1
         NCLEFT = NC
 2020    NCP = MIN0 (NCLEFT,NCMX)
         IF (NL .GT. 20) PRINT 9001
         DO 2011 I = 1,NL
            DO 2012 J = 1,NCP
               IJ = I+(JS+J)*NLMX
               AA8(J) = AA(2*IJ+IP-2)
 2012       CONTINUE
            CALL LPRINT (AA8,AA4,I,NCP,IFORM,BLANK)
 2011    CONTINUE
         NCLEFT = NCLEFT-NCMX
         JS = JS+NCMX
         IF (NCLEFT .LE. 0) GO TO 2010
         GO TO 2020
 2010 CONTINUE
      RETURN
C
C --- COMPLEX NUMBERS (DOUBLE PRECISION)
C
 2100 DO 2110 IP = 1,2
         PRINT 9010,ANC(IP),ANAME
         JS = -1
         NCLEFT = NC
 2120    NCP = MIN0 (NCLEFT,NCMX)
         IF (NL .GT. 20) PRINT 9001
         DO 2111 I = 1,NL
            DO 2112 J = 1,NCP
               IJ = I+(JS+J)*NLMX
               AA4(2*J-1) = AA(2*(2*IJ+IP-2)-1)
               AA4(2*J)   = AA(2*(2*IJ+IP-2)  )
 2112       CONTINUE
            CALL LPRINT (AA8,AA4,I,NCP,IFORM,BLANK)
 2111    CONTINUE
         NCLEFT = NCLEFT-NCMX
         JS = JS+NCMX
         IF (NCLEFT .LE. 0) GO TO 2110
         GO TO 2120
 2110 CONTINUE
      RETURN
C
C --- INTEGER NUMBERS
C
 4000 IF (IFORM .EQ.  4) NCMX = 20
      IF (IFORM .EQ. 14) NCMX = 10
      JS = -1
      NCLEFT = NC
 4010 NCP = MIN0 (NCLEFT,NCMX)
      ANAMP = ANAME
      IF (NL .GT. 20) PRINT 9001
      DO 4011 I = 1,NL
         DO 4012 J = 1,NCP
            IJ = I+(JS+J)*NLMX
            AA4(J) = AA(IJ)
 4012    CONTINUE
         CALL LPRINT (AA8,AA4,I,NCP,IFORM,ANAMP)
         ANAMP = BLANK
 4011 CONTINUE
      NCLEFT = NCLEFT-NCMX
      JS = JS+NCMX
      IF (NCLEFT .LE. 0) RETURN
      GO TO 4010
 8000 PRINT 9100,NL,NC,NLMX
      RETURN
C
 9000 FORMAT (' ITYPE = ',I10,' SHOULD BE EQUAL TO 4,8,24 OR 28 ;'
     &       ,' MATPRT SUBROUTINE HAS NO ACTION')
 9001 FORMAT ('      ')
 9002 FORMAT (' WARNING; MATPRT DID NOT WRITE MATRIX NAMED ',A8
     &       ,' BECAUSE NUMBER OF COLUMNS IS ',I3)
 9010 FORMAT (A30,A8)
 9100 FORMAT (' WARNING ---- MATPRT ROUTINE; NL,NC,NLMX = ',3I8
     &       ,' ; TOO MANY LINES TO PRINT')
 9200 FORMAT (/1X)
      END



      SUBROUTINE MCHINE (CRAY)
      LOGICAL CRAY
      CRAY = .FALSE.
      RETURN
      END

 
 
      SUBROUTINE LPRINT (AA,IA,IL,NC,IFORM,ANAME)
C#######################################################################
C#    PRINTS OUT A LINE OF INTEGER OR REAL NUMBERS                     #
C#---------------------------------------------------------------------#
C#    AA    : LINE TO BE PRINTED  (REAL REPRESENTATION).               #
C#    IA    : LINE TO BE PRINTED  (INTEGER REPRESENTATION).            #
C#    NC    : NUMBER OF ELEMENTS IN THE LINE.                          #
C#    IFORM : DESCRIBES THE FORMAT CODE                                #
C#            1 FOR F FORMAT CODE, 2 FOR E FORMAT CODE,                #
C#            3 FOR G FORMAT CODE, 4 FOR I FORMAT CODE,                #
C#            ADD 10 TO GET HIGH NUMBER OF DIGITS PRINTOUT.            #
C#---------------------------------------------------------------------#
C#    J.M.L. : 10/1982                                                 #
C#######################################################################
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 ANAME
      DIMENSION AA(*),IA(*)
C
      IGO = IFORM
      IF (IGO .GT. 10) IGO = IGO-6
      GO TO (100,200,300,400,1100,1200,1300,1400),IGO
 100  PRINT 9010, ANAME,IL,(AA(J),J = 1,NC)
      RETURN
 200  PRINT 9020, ANAME,IL,(AA(J),J = 1,NC)
      RETURN
 300  PRINT 9030, ANAME,IL,(AA(J),J = 1,NC)
      RETURN
 400  PRINT 9040, ANAME,IL,(IA(J),J = 1,NC)
      RETURN
 1100 PRINT 9110, ANAME,IL,(AA(J),J = 1,NC)
      RETURN
 1200 PRINT 9120, ANAME,IL,(AA(J),J = 1,NC)
      RETURN
 1300 PRINT 9130, ANAME,IL,(AA(J),J = 1,NC)
      RETURN
 1400 PRINT 9140, ANAME,IL,(IA(J),J = 1,NC)
      RETURN
C
 9010 FORMAT (' ',A8,I3,10F12.6)
 9020 FORMAT (' ',A8,I3,1P,10E12.4)
 9030 FORMAT (' ',A8,I3,1P,10G12.5)
 9040 FORMAT (' ',A8,I3,20I6)
 9110 FORMAT (' ',A8,I3,5F24.16)
 9120 FORMAT (' ',A8,I3,1P,5E24.16)
 9130 FORMAT (' ',A8,I3,1P,5G24.16)
 9140 FORMAT (' ',A8,I3,10I12)
      END



      SUBROUTINE MTCOPY (AA,BB,II,JJ,NA,NB)
      DOUBLE PRECISION AA(NA,*),BB(NB,*)
      REAL   A4(NA,*),B4(NB,*)
C
      DO 10 J = 1,JJ
      DO 10 I = 1,II
      AA(I,J) = BB(I,J)
 10   CONTINUE
      RETURN
C
      ENTRY MTCOP8 (AA,BB,II,JJ,NA,NB,IKEY)
      IF (IKEY.GT.0) GO TO 30
      DO 20 J=1,JJ
      DO 20 I=1,II
      AA(I,J) = BB(I,J)
 20   CONTINUE
      RETURN
C
 30   DO 40 J=1,JJ
      DO 40 I=1,II
      AA(J,I) = BB(I,J)
 40   CONTINUE
      RETURN
C
      ENTRY MTCOP4 (A4,B4,II,JJ,NA,NB,IKEY)
      IF (IKEY.GT.0) GO TO 60
      DO 50 J=1,JJ
      DO 50 I=1,II
      A4(I,J) = B4(I,J)
 50   CONTINUE
      RETURN
C
 60   DO 70 J=1,JJ
      DO 70 I=1,II
      A4(J,I) = B4(I,J)
 70   CONTINUE
      RETURN
      END



      SUBROUTINE MTINIT (A8,II,JJ,NA)
      DOUBLE PRECISION A8,DFACT
      DIMENSION A8(NA,*),A4(NA,*)
C
      DO 5 J = 1,JJ
         DO 6 I = 1,II
            A8(I,J) = 0.D0
 6       CONTINUE
 5    CONTINUE
      RETURN
C
      ENTRY MTINI8 (A8,II,JJ,NA,IKEY)
      DFACT   = DFLOAT(IKEY)
      DO 30 J = 1,JJ
      DO 20 I = 1,II
      A8(I,J) = 0.D0
 20   CONTINUE
      IF (J.GT.II) GO TO 30
      A8(J,J) = DFACT
 30   CONTINUE
      RETURN
C
      ENTRY MTINI4 (A4,II,JJ,NA,IKEY)
      FACT = FLOAT(IKEY)
      DO 50 J=1,JJ
      DO 40 I=1,II
      A4(I,J) = 0.E0
 40   CONTINUE
      IF (J.GT.II) GO TO 50
      A4(J,J) = FACT
 50   CONTINUE
      RETURN
C
      ENTRY M4INIT (A4,II,JJ,NA)
      DO 60 J=1,JJ
      DO 61 I=1,II
      A4(I,J) = 0.E0
 61   CONTINUE
 60   CONTINUE
      RETURN
      END



      SUBROUTINE MTMULT (A,B,C, L,M,N, NA,NB,NC)
C#######################################################################
C#    PERFORMS MATRIX MULTIPLICATION                                   #
C#         A      = B        * C         : ENTRY MTMULT                #
C#          (L,N)    (L,M)      (M,N)                                  #
C#                        (T)                                          #
C#         A      = B        * C         : ENTRY MTMULL                #
C#          (L,N)    (M,L)      (M,N)                                  #
C#                                   (T)                               #
C#         A      = B        * C         : ENTRY MTMULR                #
C#          (L,N)    (L,M)      (N,M)                                  #
C#---------------------------------------------------------------------#
C#    NA,NB,NC : ROW DIMENSIONS OF A,B,C IN THE CALLING PROGRAM        #
C#######################################################################
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NA,*),B(NB,*),C(NC,*)
C
      DO 100 I = 1,L
         DO 101 J = 1,N
            SUM = 0.D0
            DO 102 K = 1,M
               SUM = SUM+B(I,K)*C(K,J)
 102        CONTINUE
            A(I,J) = SUM
 101     CONTINUE
 100  CONTINUE
      RETURN
C
      ENTRY MTMULL (A,B,C, L,M,N, NA,NB,NC)
      DO 200 I = 1,L
         DO 201 J = 1,N
            SUM = 0.D0
            DO 202 K = 1,M
               SUM = SUM+B(K,I)*C(K,J)
 202        CONTINUE
            A(I,J) = SUM
 201     CONTINUE
 200  CONTINUE
      RETURN
C
      ENTRY MTMULR (A,B,C, L,M,N, NA,NB,NC)
      DO 300 I = 1,L
         DO 301 J = 1,N
            SUM = 0.D0
            DO 302 K = 1,M
               SUM = SUM+B(I,K)*C(J,K)
 302        CONTINUE
            A(I,J) = SUM
 301     CONTINUE
 300  CONTINUE
      RETURN
      END




