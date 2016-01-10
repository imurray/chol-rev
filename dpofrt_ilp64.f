      SUBROUTINE DPOFRT( UPLO, N, A, LDA, C, LDC, INFO )
*
*  -- Iain Murray, January 2016. Derived in part from DPOTRF, which is a  --
*  -- LAPACK routine (version 3.3.1)                                      --
*  -- LAPACK is a software package provided by Univ. of Tennessee,        --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2011                                                          --
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER*8          INFO, LDA, LDC, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
      DOUBLE PRECISION   C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DPOFRT computes the reverse mode sensitivity of DPOTRF, the Cholesky
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
*  This is the block version of the algorithm, calling Level 3 BLAS.
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
*          a previous call to DPOTRF. If UPLO = 'U', the leading n by n
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
*          WARNING: In the current version, some of the elements in the
*          unreferenced diagonal will be overwritten with zeros.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C.  LDC >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*               or if INFO = -8, we hit a not-implemented error, NI_ERR
*          > 0: if INFO = k, A( k, k ) was found to be -ve, and so
*               cannot have been returned from DPOTRF.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      INTEGER*8          ONE_I, NI_ERR
      PARAMETER          ( ONE_I = 1, NI_ERR = -8 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER*8          J, JB, JL, NB, K, L
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER*8          ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DSCAL, DGEMM, DPO2FT, DSYRK, DTRSM
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
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
         CALL XERBLA( 'DPOFRT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( ONE_I, 'DPOTRF', UPLO, N, -ONE_I, -ONE_I, -ONE_I )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
         CALL DPO2FT( UPLO, N, A, LDA, C, LDC, INFO )
      ELSE
*
*        Use blocked code.
*
         IF( UPPER ) THEN
*
*           Assume A is from a factorization X = A**T *A.
*
            DO 10 JL = N-NB+1, -NB, -NB
*              NOT IMPLEMENTED YET.
               CALL XERBLA( 'DPOFRT', NI_ERR )
               RETURN
   10       CONTINUE
*
         ELSE
*
*           Assume A is from a factorization X = A*A**T.
*
            DO 20 JL = N-NB+1, 1-NB, -NB
               J = MAX( ONE_I, JL )
               JB = NB - (J - JL)
               IF( J+JB.LE.N ) THEN
*                 % Lb_B = Lb_B/L_T
*                 L_bar(J+JB:N, J:J+JB-1) = L_bar(J+JB:N, J:J+JB-1)/L(J:J+JB-1, J:J+JB-1);
                  CALL DTRSM( 'Right', 'Lower', 'No transpose',
     $                        'Non-unit', N-J-JB+1, JB, ONE, A( J, J ),
     $                         LDA, C( J+JB, J ), LDC )
*                 % Lb_T -= tril(Lb_B'*L_B)
*                 L_bar(J:J+JB-1, J:J+JB-1) = L_bar(J:J+JB-1, J:J+JB-1) - tril(L_bar(J+JB:N, J:J+JB-1)'*L(J+JB:N, J:J+JB-1));
*                 % But omitting tril(.) and overwriting with stuff and then zeros later.
*                 WASTEFUL: only need bottom left triangle of answer.
                  CALL DGEMM( 'Transpose', 'No transpose', JB, JB,
     $                        N-J-JB+1, -ONE, C( J+JB, J ), LDC,
     $                        A( J+JB, J ), LDA, ONE, C( J, J ), LDC )
               END IF
*              % Lb_T = dpotf2_rev(L_T, Lb_T)
*              L_bar(J:J+JB-1, J:J+JB-1) = dpotf2_rev_inplace(L(J:J+JB-1, J:J+JB-1), L_bar(J:J+JB-1, J:J+JB-1));
               CALL DPO2FT( UPLO, JB, A( J, J ), LDA, C( J, J),
     $                         LDC, INFO )
               IF( INFO.NE.0 )
     $            GO TO 30
*              Copy lower triangle into upper triangle for update of block C below,
*              and temporarily double diagonal element
*              There's probably a better way to do it...
               DO 25 K = J+1, J+JB-1
                  CALL DCOPY ( JB-(K-J), C( K, K-1 ), ONE_I,
     $                          C( K-1, K ), LDC )
   25          CONTINUE
               IF( J.GT.1 ) THEN
                  CALL DSCAL( JB, 2.0D+0, C( J, J ), LDC+1 )
*                 % Lb_D -= Lb_B*L_C
*                 L_bar(J+JB:N, 1:J-1) = L_bar(J+JB:N, 1:J-1) - L_bar(J+JB:N, J:J+JB-1)*L(J:J+JB-1, 1:J-1);
                  CALL DGEMM( 'No transpose', 'No transpose', N-J-JB+1,
     $                        J-1, JB, -ONE, C( J+JB, J ), LDC,
     $                        A( J, 1 ), LDA, ONE, C( J+JB, 1 ), LDC )
*                 % Lb_C -= Lb_B'*L_D
*                 L_bar(J:J+JB-1, 1:J-1) = L_bar(J:J+JB-1, 1:J-1) - L_bar(J+JB:N, J:J+JB-1)'*L(J+JB:N, 1:J-1);
                  CALL DGEMM( 'Transpose', 'No transpose', JB, J-1,
     $                        N-J-JB+1, -ONE, C( J+JB, J ), LDC,
     $                        A( J+JB, 1 ), LDA, ONE, C( J, 1 ), LDC )
*                 % Lb_C -= (Lb_T + Lb_T')*C
*                 L_bar(J:J+JB-1, 1:J-1) = L_bar(J:J+JB-1, 1:J-1) - (L_bar(J:J+JB-1, J:J+JB-1) + L_bar(J:J+JB-1, J:J+JB-1)')*L(J:J+JB-1, 1:J-1);
                  CALL DGEMM( 'No transpose', 'No transpose', JB, J-1,
     $                        JB, -ONE, C( J, J ), LDC,
     $                        A( J, 1 ), LDA, ONE, C( J, 1 ), LDC )
*                 Rescale diagonal back again
                  CALL DSCAL( JB, 0.5D+0, C( J, J ), LDC+1 )
               END IF
*              Zero unwanted bits above diagonal
               DO 26 K = J+1, J+JB-1
                  DO 27 L = 1, K-1
                     C( L, K ) = 0.0D+0
   27             CONTINUE
   26          CONTINUE
*
   20       CONTINUE
         END IF
      END IF
      GO TO 40
*
   30 CONTINUE
      INFO = INFO + J - 1
*
   40 CONTINUE
      RETURN
*
*     End of DPOFRT
*
      END
