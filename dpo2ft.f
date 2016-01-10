      SUBROUTINE DPO2FT( UPLO, N, A, LDA, C, LDC, INFO )
*
*  -- Iain Murray, January 2016. Derived in part from DPOTF2, which is a  --
*  -- LAPACK routine (version 3.3.1)                                      --
*  -- LAPACK is a software package provided by Univ. of Tennessee,        --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2011                                                          --
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDC, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
      DOUBLE PRECISION   C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DPO2FT computes the reverse mode sensitivity of DPOTF2, the Cholesky
*  factorization A of a real symmetric positive definite matrix X.
*
*  The original factorization had the form
*     X = A**T * A ,  if UPLO = 'U', or
*     X = A  * A**T,  if UPLO = 'L',
*  where A is an upper or lower triangular matrix respectively.
*
*  The triangular matrix C provides sensitivities of a scalar cost
*  function with respect to X, the output of the Cholesky algorithm A.
*  This routine transforms C into a triangular matrix of sensitivities
*  with respect to the input matrix X.
*
*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          matrices A and C are to be used.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrices A and C.  N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          A contains a Cholesky decomposition of a matrix, usually from
*          a previous call to DPOTF2. If UPLO = 'U', the leading n by n
*          upper triangular part of A is used, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of A is used, and the
*          strictly upper triangular part of A is not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*
*          On exit, if INFO = 0, the sensitivity of the cost function
*          with respect to the full matrix X = A**T *A  or X = A*A**T.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C.  LDC >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*               or if INFO = -8, we hit a not-implemented error, NI_ERR
*          > 0: if INFO = k, A( k, k ) was found to be -ve, and so
*               cannot have been returned from DPOTF2.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER            ONE_I, NI_ERR
      PARAMETER          ( ONE_I = 1, NI_ERR = -8 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT, DISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DSCAL, XERBLA, DGER
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( ONE_I, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPO2FT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Assume A is from a factorization X = A**T *A.
*
         DO 10 J = N, 1, -1
*           NOT IMPLEMENTED YET.
            CALL XERBLA( 'DPO2FT', NI_ERR )
            RETURN
   10    CONTINUE
      ELSE
*
*        Assume A is from a factorization X = A*A**T.
*
         DO 20 J = N, 1, -1
            AJJ = A( J, J )
            IF( AJJ.LE.ZERO ) THEN
               GO TO 30
            END IF
*           C(J,J) = C(J,J) - A(J+1:N,J)'*C(J+1:N,J) / A(J,J);
            C( J, J ) = C( J, J ) - DDOT(
     $            N-J, A( J+1, J ), ONE_I, C( J+1, J ), ONE_I ) / AJJ
*           C(J:N,J) = C(J:N,J) / A(J,J);
            CALL DSCAL( N-J+1, ONE / AJJ, C( J, J ), ONE_I )
*           C(J,1:J-1) = C(J,1:J-1) - C(J:N,J)'*A(J:N,1:J-1);
            CALL DGEMV( 'Transpose', N-J+1, J-1, -ONE, A( J, 1 ), LDA,
     $            C( J, J ), ONE_I, ONE, C( J, 1 ), LDC )
*           C(J+1:N,1:J-1) = C(J+1:N,1:J-1) - C(J+1:N,J)*A(J,1:J-1);
            CALL DGER( N-J, J-1, -ONE, C( J+1, J ), ONE_I, A( J, 1 ),
     $            LDA, C( J+1, 1 ), LDC)
*           %C(J,J) = 0.5 * C(J,J); % can take out of loop if like.
            C( J, J ) = 0.5D+0 * C( J, J )
   20    CONTINUE
      END IF
      GO TO 40
*
   30 CONTINUE
      INFO = J
*
   40 CONTINUE
      RETURN
*
*     End of DPO2FT
*
      END
