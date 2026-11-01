

#ifndef __GMS_SLICOT_C_IFACE_H__
#define __GMS_SLICOT_C_IFACE_H__ 011020260721


/******************************************************
  C - interface to SLICOT F77 Control Library.
  Wrapper routines declarations.
*******************************************************/
#include <complex>

#if defined (__cplusplus)

extern "C" {

		// Subroutine extended description.
	/*
	  C     SLICOT RELEASE 5.5.
C
C     Copyright (c) 2002-2012 NICONET e.V.
C
C     PURPOSE
C
C     To find a controllable realization for the linear time-invariant
C     single-input system
C
C             dX/dt = A * X + B * U,
C
C     where A is an N-by-N matrix and B is an N element vector which
C     are reduced by this routine to orthogonal canonical form using
C     (and optionally accumulating) orthogonal similarity
C     transformations.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBZ    CHARACTER*1
C             Indicates whether the user wishes to accumulate in a
C             matrix Z the orthogonal similarity transformations for
C             reducing the system, as follows:
C             = 'N':  Do not form Z and do not store the orthogonal
C                     transformations;
C             = 'F':  Do not form Z, but store the orthogonal
C                     transformations in the factored form;
C             = 'I':  Z is initialized to the unit matrix and the
C                     orthogonal transformation matrix Z is returned.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the original state-space representation,
C             i.e. the order of the matrix A.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the original state dynamics matrix A.
C             On exit, the leading NCONT-by-NCONT upper Hessenberg
C             part of this array contains the canonical form of the
C             state dynamics matrix, given by Z' * A * Z, of a
C             controllable realization for the original system. The
C             elements below the first subdiagonal are set to zero.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, the original input/state vector B.
C             On exit, the leading NCONT elements of this array contain
C             canonical form of the input/state vector, given by Z' * B,
C             with all elements but B(1) set to zero.
C
C     NCONT   (output) INTEGER
C             The order of the controllable state-space representation.
C
C     Z       (output) DOUBLE PRECISION array, dimension (LDZ,N)
C             If JOBZ = 'I', then the leading N-by-N part of this array
C             contains the matrix of accumulated orthogonal similarity
C             transformations which reduces the given system to
C             orthogonal canonical form.
C             If JOBZ = 'F', the elements below the diagonal, with the
C             array TAU, represent the orthogonal transformation matrix
C             as a product of elementary reflectors. The transformation
C             matrix can then be obtained by calling the LAPACK Library
C             routine DORGQR.
C             If JOBZ = 'N', the array Z is not referenced and can be
C             supplied as a dummy array (i.e. set parameter LDZ = 1 and
C             declare this array to be Z(1,1) in the calling program).
C
C     LDZ     INTEGER
C             The leading dimension of array Z. If JOBZ = 'I' or
C             JOBZ = 'F', LDZ >= MAX(1,N); if JOBZ = 'N', LDZ >= 1.
C
C     TAU     (output) DOUBLE PRECISION array, dimension (N)
C             The elements of TAU contain the scalar factors of the
C             elementary reflectors used in the reduction of B and A.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used in determining the
C             controllability of (A,B). If the user sets TOL > 0, then
C             the given value of TOL is used as an absolute tolerance;
C             elements with absolute value less than TOL are considered
C             neglijible. If the user sets TOL <= 0, then an implicitly
C             computed, default tolerance, defined by
C             TOLDEF = N*EPS*MAX( NORM(A), NORM(B) ) is used instead,
C             where EPS is the machine precision (see LAPACK Library
C             routine DLAMCH).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK. LDWORK >= MAX(1,N).
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     The Householder matrix which reduces all but the first element
C     of vector B to zero is found and this orthogonal similarity
C     transformation is applied to the matrix A. The resulting A is then
C     reduced to upper Hessenberg form by a sequence of Householder
C     transformations. Finally, the order of the controllable state-
C     space representation (NCONT) is determined by finding the position
C     of the first sub-diagonal element of A which is below an
C     appropriate zero threshold, either TOL or TOLDEF (see parameter
C     TOL); if NORM(B) is smaller than this threshold, NCONT is set to
C     zero, and no computations for reducing the system to orthogonal
C     canonical form are performed.
C
C     REFERENCES
C
C     [1] Konstantinov, M.M., Petkov, P.Hr. and Christov, N.D.
C         Orthogonal Invariants and Canonical Forms for Linear
C         Controllable Systems.
C         Proc. 8th IFAC World Congress, Kyoto, 1, pp. 49-54, 1981.
C
C     [2] Hammarling, S.J.
C         Notes on the use of orthogonal similarity transformations in
C         control.
C         NPL Report DITC 8/82, August 1982.
C
C     [3] Paige, C.C
C         Properties of numerical algorithms related to computing
C         controllability.
C         IEEE Trans. Auto. Contr., AC-26, pp. 130-138, 1981.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations and is backward stable.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Oct. 1996.
C     Supersedes Release 2.0 routine AB01AD by T.W.C. Williams,
C     Kingston Polytechnic, United Kingdom, October 1982.
C
C     REVISIONS
C
C     V. Sima, February 16, 1998, October 19, 2001, February 2, 2005.
C
C     KEYWORDS
C
C     Controllability, minimal realization, orthogonal canonical form,
C     orthogonal transformation.
C
C     ******************************************************************
C
	*/

		void AB01MD(unsigned char *, int *, double *, int *, double *, int *,
		            double *, int *, double *, double *, double *, int *, int *);


		/*
		 C     SLICOT RELEASE 5.5.
C
C     Copyright (c) 2002-2012 NICONET e.V.
C
C     PURPOSE
C
C     To find a controllable realization for the linear time-invariant
C     multi-input system
C
C             dX/dt = A * X + B * U,
C
C     where A and B are N-by-N and N-by-M matrices, respectively,
C     which are reduced by this routine to orthogonal canonical form
C     using (and optionally accumulating) orthogonal similarity
C     transformations.  Specifically, the pair (A, B) is reduced to
C     the pair (Ac, Bc),  Ac = Z' * A * Z,  Bc = Z' * B,  given by
C
C             [ Acont     *    ]         [ Bcont ]
C        Ac = [                ],   Bc = [       ],
C             [   0    Auncont ]         [   0   ]
C
C        and
C
C                [ A11 A12  . . .  A1,p-1 A1p ]         [ B1 ]
C                [ A21 A22  . . .  A2,p-1 A2p ]         [ 0  ]
C                [  0  A32  . . .  A3,p-1 A3p ]         [ 0  ]
C        Acont = [  .   .   . . .    .     .  ],   Bc = [ .  ],
C                [  .   .     . .    .     .  ]         [ .  ]
C                [  .   .       .    .     .  ]         [ .  ]
C                [  0   0   . . .  Ap,p-1 App ]         [ 0  ]
C
C     where the blocks  B1, A21, ..., Ap,p-1  have full row ranks and
C     p is the controllability index of the pair.  The size of the
C     block  Auncont is equal to the dimension of the uncontrollable
C     subspace of the pair (A, B).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBZ    CHARACTER*1
C             Indicates whether the user wishes to accumulate in a
C             matrix Z the orthogonal similarity transformations for
C             reducing the system, as follows:
C             = 'N':  Do not form Z and do not store the orthogonal
C                     transformations;
C             = 'F':  Do not form Z, but store the orthogonal
C                     transformations in the factored form;
C             = 'I':  Z is initialized to the unit matrix and the
C                     orthogonal transformation matrix Z is returned.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the original state-space representation,
C             i.e. the order of the matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs, or of columns of B.  M >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the original state dynamics matrix A.
C             On exit, the leading NCONT-by-NCONT part contains the
C             upper block Hessenberg state dynamics matrix Acont in Ac,
C             given by Z' * A * Z, of a controllable realization for
C             the original system. The elements below the first block-
C             subdiagonal are set to zero.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the input matrix B.
C             On exit, the leading NCONT-by-M part of this array
C             contains the transformed input matrix Bcont in Bc, given
C             by Z' * B, with all elements but the first block set to
C             zero.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     NCONT   (output) INTEGER
C             The order of the controllable state-space representation.
C
C     INDCON  (output) INTEGER
C             The controllability index of the controllable part of the
C             system representation.
C
C     NBLK    (output) INTEGER array, dimension (N)
C             The leading INDCON elements of this array contain the
C             the orders of the diagonal blocks of Acont.
C
C     Z       (output) DOUBLE PRECISION array, dimension (LDZ,N)
C             If JOBZ = 'I', then the leading N-by-N part of this
C             array contains the matrix of accumulated orthogonal
C             similarity transformations which reduces the given system
C             to orthogonal canonical form.
C             If JOBZ = 'F', the elements below the diagonal, with the
C             array TAU, represent the orthogonal transformation matrix
C             as a product of elementary reflectors. The transformation
C             matrix can then be obtained by calling the LAPACK Library
C             routine DORGQR.
C             If JOBZ = 'N', the array Z is not referenced and can be
C             supplied as a dummy array (i.e. set parameter LDZ = 1 and
C             declare this array to be Z(1,1) in the calling program).
C
C     LDZ     INTEGER
C             The leading dimension of array Z. If JOBZ = 'I' or
C             JOBZ = 'F', LDZ >= MAX(1,N); if JOBZ = 'N', LDZ >= 1.
C
C     TAU     (output) DOUBLE PRECISION array, dimension (N)
C             The elements of TAU contain the scalar factors of the
C             elementary reflectors used in the reduction of B and A.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used in rank determination when
C             transforming (A, B). If the user sets TOL > 0, then
C             the given value of TOL is used as a lower bound for the
C             reciprocal condition number (see the description of the
C             argument RCOND in the SLICOT routine MB03OD);  a
C             (sub)matrix whose estimated condition number is less than
C             1/TOL is considered to be of full rank.  If the user sets
C             TOL <= 0, then an implicitly computed, default tolerance,
C             defined by  TOLDEF = N*N*EPS,  is used instead, where EPS
C             is the machine precision (see LAPACK Library routine
C             DLAMCH).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (M)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1, N, 3*M).
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     Matrix B is first QR-decomposed and the appropriate orthogonal
C     similarity transformation applied to the matrix A. Leaving the
C     first rank(B) states unchanged, the remaining lower left block
C     of A is then QR-decomposed and the new orthogonal matrix, Q1,
C     is also applied to the right of A to complete the similarity
C     transformation. By continuing in this manner, a completely
C     controllable state-space pair (Acont, Bcont) is found for the
C     given (A, B), where Acont is upper block Hessenberg with each
C     subdiagonal block of full row rank, and Bcont is zero apart from
C     its (independent) first rank(B) rows.
C     NOTE that the system controllability indices are easily
C     calculated from the dimensions of the blocks of Acont.
C
C     REFERENCES
C
C     [1] Konstantinov, M.M., Petkov, P.Hr. and Christov, N.D.
C         Orthogonal Invariants and Canonical Forms for Linear
C         Controllable Systems.
C         Proc. 8th IFAC World Congress, Kyoto, 1, pp. 49-54, 1981.
C
C     [2] Paige, C.C.
C         Properties of numerical algorithms related to computing
C         controllablity.
C         IEEE Trans. Auto. Contr., AC-26, pp. 130-138, 1981.
C
C     [3] Petkov, P.Hr., Konstantinov, M.M., Gu, D.W. and
C         Postlethwaite, I.
C         Optimal Pole Assignment Design of Linear Multi-Input Systems.
C         Leicester University, Report 99-11, May 1996.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations and is backward stable.
C
C     FURTHER COMMENTS
C
C     If the system matrices A and B are badly scaled, it would be
C     useful to scale them with SLICOT routine TB01ID, before calling
C     the routine.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996.
C     Supersedes Release 2.0 routine AB01BD by P.Hr. Petkov.
C
C     REVISIONS
C
C     January 14, 1997, June 4, 1997, February 13, 1998,
C     September 22, 2003, February 29, 2004.
C
C     KEYWORDS
C
C     Controllability, minimal realization, orthogonal canonical form,
C     orthogonal transformation.
C
		*/

		void AB01ND(unsigned char *, int *, int *, double *, int *, double *, int *,
		            int *, int *, int *, double *, int *, double *, double *, int *,
					double *, int *, int *);

		/*
		 SLICOT RELEASE 5.5.
C
C     Copyright (c) 2002-2012 NICONET e.V.
C
C     PURPOSE
C
C     To reduce the matrices A and B using (and optionally accumulating)
C     state-space and input-space transformations U and V respectively,
C     such that the pair of matrices
C
C        Ac = U' * A * U,    Bc = U' * B * V
C
C     are in upper "staircase" form. Specifically,
C
C             [ Acont     *    ]         [ Bcont ]
C        Ac = [                ],   Bc = [       ],
C             [   0    Auncont ]         [   0   ]
C
C        and
C
C                [ A11 A12  . . .  A1,p-1 A1p ]         [ B1 ]
C                [ A21 A22  . . .  A2,p-1 A2p ]         [ 0  ]
C                [  0  A32  . . .  A3,p-1 A3p ]         [ 0  ]
C        Acont = [  .   .   . . .    .     .  ],   Bc = [ .  ],
C                [  .   .     . .    .     .  ]         [ .  ]
C                [  .   .       .    .     .  ]         [ .  ]
C                [  0   0   . . .  Ap,p-1 App ]         [ 0  ]
C
C     where the blocks  B1, A21, ..., Ap,p-1  have full row ranks and
C     p is the controllability index of the pair.  The size of the
C     block Auncont is equal to the dimension of the uncontrollable
C     subspace of the pair (A, B).  The first stage of the reduction,
C     the "forward" stage, accomplishes the reduction to the orthogonal
C     canonical form (see SLICOT library routine AB01ND). The blocks
C     B1, A21, ..., Ap,p-1 are further reduced in a second, "backward"
C     stage to upper triangular form using RQ factorization. Each of
C     these stages is optional.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     STAGES  CHARACTER*1
C             Specifies the reduction stages to be performed as follows:
C             = 'F':  Perform the forward stage only;
C             = 'B':  Perform the backward stage only;
C             = 'A':  Perform both (all) stages.
C
C     JOBU    CHARACTER*1
C             Indicates whether the user wishes to accumulate in a
C             matrix U the state-space transformations as follows:
C             = 'N':  Do not form U;
C             = 'I':  U is internally initialized to the unit matrix (if
C                     STAGES <> 'B'), or updated (if STAGES = 'B'), and
C                     the orthogonal transformation matrix U is
C                     returned.
C
C     JOBV    CHARACTER*1
C             Indicates whether the user wishes to accumulate in a
C             matrix V the input-space transformations as follows:
C             = 'N':  Do not form V;
C             = 'I':  V is initialized to the unit matrix and the
C                     orthogonal transformation matrix V is returned.
C             JOBV is not referenced if STAGES = 'F'.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The actual state dimension, i.e. the order of the
C             matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The actual input dimension.  M >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state transition matrix A to be transformed.
C             If STAGES = 'B', A should be in the orthogonal canonical
C             form, as returned by SLICOT library routine AB01ND.
C             On exit, the leading N-by-N part of this array contains
C             the transformed state transition matrix U' * A * U.
C             The leading NCONT-by-NCONT part contains the upper block
C             Hessenberg state matrix Acont in Ac, given by U' * A * U,
C             of a controllable realization for the original system.
C             The elements below the first block-subdiagonal are set to
C             zero.  If STAGES <> 'F', the subdiagonal blocks of A are
C             triangularized by RQ factorization, and the annihilated
C             elements are explicitly zeroed.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the input matrix B to be transformed.
C             If STAGES = 'B', B should be in the orthogonal canonical
C             form, as returned by SLICOT library routine AB01ND.
C             On exit with STAGES = 'F', the leading N-by-M part of
C             this array contains the transformed input matrix U' * B,
C             with all elements but the first block set to zero.
C             On exit with STAGES <> 'F', the leading N-by-M part of
C             this array contains the transformed input matrix
C             U' * B * V, with all elements but the first block set to
C             zero and the first block in upper triangular form.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     U       (input/output) DOUBLE PRECISION array, dimension (LDU,N)
C             If STAGES <> 'B' or JOBU = 'N', then U need not be set
C             on entry.
C             If STAGES = 'B' and JOBU = 'I', then, on entry, the
C             leading N-by-N part of this array must contain the
C             transformation matrix U that reduced the pair to the
C             orthogonal canonical form.
C             On exit, if JOBU = 'I', the leading N-by-N part of this
C             array contains the transformation matrix U that performed
C             the specified reduction.
C             If JOBU = 'N', the array U is not referenced and can be
C             supplied as a dummy array (i.e. set parameter LDU = 1 and
C             declare this array to be U(1,1) in the calling program).
C
C     LDU     INTEGER
C             The leading dimension of array U.
C             If JOBU = 'I', LDU >= MAX(1,N);  if JOBU = 'N', LDU >= 1.
C
C     V       (output) DOUBLE PRECISION array, dimension (LDV,M)
C             If JOBV = 'I', then the leading M-by-M part of this array
C             contains the transformation matrix V.
C             If STAGES = 'F', or JOBV = 'N', the array V is not
C             referenced and can be supplied as a dummy array (i.e. set
C             parameter  LDV = 1 and declare this array to be V(1,1) in
C             the calling program).
C
C     LDV     INTEGER
C             The leading dimension of array V.
C             If STAGES <> 'F' and JOBV = 'I', LDV >= MAX(1,M);
C             if STAGES = 'F' or JOBV = 'N', LDV >= 1.
C
C     NCONT   (input/output) INTEGER
C             The order of the controllable state-space representation.
C             NCONT is input only if STAGES = 'B'.
C
C     INDCON  (input/output) INTEGER
C             The number of stairs in the staircase form (also, the
C             controllability index of the controllable part of the
C             system representation).
C             INDCON is input only if STAGES = 'B'.
C
C     KSTAIR  (input/output) INTEGER array, dimension (N)
C             The leading INDCON elements of this array contain the
C             dimensions of the stairs, or, also, the orders of the
C             diagonal blocks of Acont.
C             KSTAIR is input if STAGES = 'B', and output otherwise.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used in rank determination when
C             transforming (A, B). If the user sets TOL > 0, then
C             the given value of TOL is used as a lower bound for the
C             reciprocal condition number (see the description of the
C             argument RCOND in the SLICOT routine MB03OD);  a
C             (sub)matrix whose estimated condition number is less than
C             1/TOL is considered to be of full rank.  If the user sets
C             TOL <= 0, then an implicitly computed, default tolerance,
C             defined by  TOLDEF = N*N*EPS,  is used instead, where EPS
C             is the machine precision (see LAPACK Library routine
C             DLAMCH).
C             TOL is not referenced if STAGES = 'B'.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (M)
C             IWORK is not referenced if STAGES = 'B'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             If STAGES <> 'B', LDWORK >= MAX(1, N + MAX(N,3*M));
C             If STAGES =  'B', LDWORK >= MAX(1, M + MAX(N,M)).
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     Staircase reduction of the pencil [B|sI - A] is used. Orthogonal
C     transformations U and V are constructed such that
C
C
C                        |B |sI-A      *  . . .  *      *       |
C                        | 1|    11       .      .      .       |
C                        |  |  A    sI-A    .    .      .       |
C                        |  |   21      22    .  .      .       |
C                        |  |        .     .     *      *       |
C     [U'BV|sI - U'AU] = |0 |     0    .     .                  |
C                        |  |            A     sI-A     *       |
C                        |  |             p,p-1    pp           |
C                        |  |                                   |
C                        |0 |         0          0   sI-A       |
C                        |  |                            p+1,p+1|
C
C
C     where the i-th diagonal block of U'AU has dimension KSTAIR(i),
C     for i = 1,...,p. The value of p is returned in INDCON. The last
C     block contains the uncontrollable modes of the (A,B)-pair which
C     are also the generalized eigenvalues of the above pencil.
C
C     The complete reduction is performed in two stages. The first,
C     forward stage accomplishes the reduction to the orthogonal
C     canonical form. The second, backward stage consists in further
C     reduction to triangular form by applying left and right orthogonal
C     transformations.
C
C     REFERENCES
C
C     [1] Van Dooren, P.
C         The generalized eigenvalue problem in linear system theory.
C         IEEE Trans. Auto. Contr., AC-26, pp. 111-129, 1981.
C
C     [2] Miminis, G. and Paige, C.
C         An algorithm for pole assignment of time-invariant multi-input
C         linear systems.
C         Proc. 21st IEEE CDC, Orlando, Florida, 1, pp. 62-67, 1982.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires O((N + M) x N**2) operations and is
C     backward stable (see [1]).
C
C     FURTHER COMMENTS
C
C     If the system matrices A and B are badly scaled, it would be
C     useful to scale them with SLICOT routine TB01ID, before calling
C     the routine.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996.
C     Supersedes Release 2.0 routine AB01CD by M. Vanbegin, and
C     P. Van Dooren, Philips Research Laboratory, Brussels, Belgium.
C
C     REVISIONS
C
C     January 14, 1997, February 12, 1998, September 22, 2003.
C
C     KEYWORDS
C
C     Controllability, generalized eigenvalue problem, orthogonal
C     transformation, staircase form.
C
		*/

		void AB01OD(unsigned char *, unsigned char *, unsigned char *, int *, int *, double *,
		           int *, double *, int *, double *, int *, double *, int *, int *, int *,
				   int *, double *, int *, double *, int *, int *);


		/*
		SLICOT RELEASE 5.5.
		C
		C     Copyright (c) 2002-2012 NICONET e.V.
		C
		C     PURPOSE
		C
		C     To perform a transformation on the parameters (A,B,C,D) of a
		C     system, which is equivalent to a bilinear transformation of the
		C     corresponding transfer function matrix.
		C
		C     ARGUMENTS
		C
		C     Mode Parameters
		C
		C     TYPE    CHARACTER*1
		C             Indicates the type of the original system and the
		C             transformation to be performed as follows:
		C             = 'D':  discrete-time   -> continuous-time;
		C             = 'C':  continuous-time -> discrete-time.
		C
		C     Input/Output Parameters
		C
		C     N       (input) INTEGER
		C             The order of the state matrix A.  N >= 0.
		C
		C     M       (input) INTEGER
		C             The number of system inputs.  M >= 0.
		C
		C     P       (input) INTEGER
		C             The number of system outputs.  P >= 0.
		C
		C     ALPHA,  (input) DOUBLE PRECISION
		C     BETA    Parameters specifying the bilinear transformation.
		C             Recommended values for stable systems: ALPHA = 1,
		C             BETA = 1.  ALPHA <> 0, BETA <> 0.
		C
		C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
		C             On entry, the leading N-by-N part of this array must
		C             contain the state matrix A of the original system.
		C             On exit, the leading N-by-N part of this array contains
		C                              _
		C             the state matrix A of the transformed system.
		C
		C     LDA     INTEGER
		C             The leading dimension of array A.  LDA >= MAX(1,N).
		C
		C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
		C             On entry, the leading N-by-M part of this array must
		C             contain the input matrix B of the original system.
		C             On exit, the leading N-by-M part of this array contains
		C                              _
		C             the input matrix B of the transformed system.
		C
		C     LDB     INTEGER
		C             The leading dimension of array B.  LDB >= MAX(1,N).
		C
		C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
		C             On entry, the leading P-by-N part of this array must
		C             contain the output matrix C of the original system.
		C             On exit, the leading P-by-N part of this array contains
		C                               _
		C             the output matrix C of the transformed system.
		C
		C     LDC     INTEGER
		C             The leading dimension of array C.  LDC >= MAX(1,P).
		C
		C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
		C             On entry, the leading P-by-M part of this array must
		C             contain the input/output matrix D for the original system.
		C             On exit, the leading P-by-M part of this array contains
		C                                     _
		C             the input/output matrix D of the transformed system.
		C
		C     LDD     INTEGER
		C             The leading dimension of array D.  LDD >= MAX(1,P).
		C
		C     Workspace
		C
		C     IWORK   INTEGER array, dimension (N)
		C
		C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
		C             On exit, if INFO = 0, DWORK(1) returns the optimal value
		C             of LDWORK.
		C
		C     LDWORK  INTEGER
		C             The length of the array DWORK.  LDWORK >= MAX(1,N).
		C             For optimum performance LDWORK >= MAX(1,N*NB), where NB
		C             is the optimal blocksize.
		C
		C     Error Indicator
		C
		C     INFO    INTEGER
		C             = 0:  successful exit;
		C             < 0:  if INFO = -i, the i-th argument had an illegal
		C                   value;
		C             = 1:  if the matrix (ALPHA*I + A) is exactly singular;
		C             = 2:  if the matrix  (BETA*I - A) is exactly singular.
		C
		C     METHOD
		C
		C     The parameters of the discrete-time system are transformed into
		C     the parameters of the continuous-time system (TYPE = 'D'), or
		C     vice-versa (TYPE = 'C') by the transformation:
		C
		C     1.  Discrete -> continuous
		C         _                     -1
		C         A = beta*(alpha*I + A)  * (A - alpha*I)
		C         _                                     -1
		C         B = sqrt(2*alpha*beta) * (alpha*I + A)  * B
		C         _                                         -1
		C         C = sqrt(2*alpha*beta) * C * (alpha*I + A)
		C         _                        -1
		C         D = D - C * (alpha*I + A)  * B
		C
		C     which is equivalent to the bilinear transformation
		C
		C                       z - alpha
		C         z -> s = beta ---------  .
		C                       z + alpha
		C
		C     of one transfer matrix onto the other.
		C
		C     2.  Continuous -> discrete
		C         _                     -1
		C         A = alpha*(beta*I - A)  * (beta*I + A)
		C         _                                    -1
		C         B = sqrt(2*alpha*beta) * (beta*I - A)  * B
		C         _                                        -1
		C         C = sqrt(2*alpha*beta) * C * (beta*I - A)
		C         _                       -1
		C         D = D + C * (beta*I - A)  * B
		C
		C     which is equivalent to the bilinear transformation
		C
		C                      beta + s
		C       s -> z = alpha -------- .
		C                      beta - s
		C
		C     of one transfer matrix onto the other.
		C
		C     REFERENCES
		C
		C     [1] Al-Saggaf, U.M. and Franklin, G.F.
		C         Model reduction via balanced realizations: a extension and
		C         frequency weighting techniques.
		C         IEEE Trans. Autom. Contr., AC-33, pp. 687-692, 1988.
		C
		C     NUMERICAL ASPECTS
		C                                                      3
		C     The time taken is approximately proportional to N .
		C     The accuracy depends mainly on the condition number of the matrix
		C     to be inverted.
		C
		C     CONTRIBUTORS
		C
		C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, and
		C                  A. Varga, German Aerospace Research Establishment,
		C                  Oberpfaffenhofen, Germany, Nov. 1996.
		C     Supersedes Release 2.0 routine AB04AD by W. van der Linden, and
		C     A.J. Geurts, Technische Hogeschool Eindhoven, Holland.
		C
		C     REVISIONS
		C
		C     -
		C
		C     KEYWORDS
		C
		C     Bilinear transformation, continuous-time system, discrete-time
		C     system, state-space model.
		C
		*/

		void AB04MD(unsigned char *, int *, int *, int *, double *, double *, double *,
		            int *, double *, int *, double *, int *, double *, int *, int *,
					double *, int *, int *);

		/*
		SLICOT RELEASE 5.5.
		C
		C     Copyright (c) 2002-2012 NICONET e.V.
		C
		C     PURPOSE
		C
		C     To obtain the state-space model (A,B,C,D) for the cascaded
		C     inter-connection of two systems, each given in state-space form.
		C
		C     ARGUMENTS
		C
		C     Mode Parameters
		C
		C     UPLO    CHARACTER*1
		C             Indicates whether the user wishes to obtain the matrix A
		C             in the upper or lower block diagonal form, as follows:
		C             = 'U':  Obtain A in the upper block diagonal form;
		C             = 'L':  Obtain A in the lower block diagonal form.
		C
		C     OVER    CHARACTER*1
		C             Indicates whether the user wishes to overlap pairs of
		C             arrays, as follows:
		C             = 'N':  Do not overlap;
		C             = 'O':  Overlap pairs of arrays: A1 and A, B1 and B,
		C                     C1 and C, and D1 and D (for UPLO = 'L'), or A2
		C                     and A, B2 and B, C2 and C, and D2 and D (for
		C                     UPLO = 'U'), i.e. the same name is effectively
		C                     used for each pair (for all pairs) in the routine
		C                     call.  In this case, setting LDA1 = LDA,
		C                     LDB1 = LDB, LDC1 = LDC, and LDD1 = LDD, or
		C                     LDA2 = LDA, LDB2 = LDB, LDC2 = LDC, and LDD2 = LDD
		C                     will give maximum efficiency.
		C
		C     Input/Output Parameters
		C
		C     N1      (input) INTEGER
		C             The number of state variables in the first system, i.e.
		C             the order of the matrix A1.  N1 >= 0.
		C
		C     M1      (input) INTEGER
		C             The number of input variables for the first system.
		C             M1 >= 0.
		C
		C     P1      (input) INTEGER
		C             The number of output variables from the first system and
		C             the number of input variables for the second system.
		C             P1 >= 0.
		C
		C     N2      (input) INTEGER
		C             The number of state variables in the second system, i.e.
		C             the order of the matrix A2.  N2 >= 0.
		C
		C     P2      (input) INTEGER
		C             The number of output variables from the second system.
		C             P2 >= 0.
		C
		C     A1      (input) DOUBLE PRECISION array, dimension (LDA1,N1)
		C             The leading N1-by-N1 part of this array must contain the
		C             state transition matrix A1 for the first system.
		C
		C     LDA1    INTEGER
		C             The leading dimension of array A1.  LDA1 >= MAX(1,N1).
		C
		C     B1      (input) DOUBLE PRECISION array, dimension (LDB1,M1)
		C             The leading N1-by-M1 part of this array must contain the
		C             input/state matrix B1 for the first system.
		C
		C     LDB1    INTEGER
		C             The leading dimension of array B1.  LDB1 >= MAX(1,N1).
		C
		C     C1      (input) DOUBLE PRECISION array, dimension (LDC1,N1)
		C             The leading P1-by-N1 part of this array must contain the
		C             state/output matrix C1 for the first system.
		C
		C     LDC1    INTEGER
		C             The leading dimension of array C1.
		C             LDC1 >= MAX(1,P1) if N1 > 0.
		C             LDC1 >= 1 if N1 = 0.
		C
		C     D1      (input) DOUBLE PRECISION array, dimension (LDD1,M1)
		C             The leading P1-by-M1 part of this array must contain the
		C             input/output matrix D1 for the first system.
		C
		C     LDD1    INTEGER
		C             The leading dimension of array D1.  LDD1 >= MAX(1,P1).
		C
		C     A2      (input) DOUBLE PRECISION array, dimension (LDA2,N2)
		C             The leading N2-by-N2 part of this array must contain the
		C             state transition matrix A2 for the second system.
		C
		C     LDA2    INTEGER
		C             The leading dimension of array A2.  LDA2 >= MAX(1,N2).
		C
		C     B2      (input) DOUBLE PRECISION array, dimension (LDB2,P1)
		C             The leading N2-by-P1 part of this array must contain the
		C             input/state matrix B2 for the second system.
		C
		C     LDB2    INTEGER
		C             The leading dimension of array B2.  LDB2 >= MAX(1,N2).
		C
		C     C2      (input) DOUBLE PRECISION array, dimension (LDC2,N2)
		C             The leading P2-by-N2 part of this array must contain the
		C             state/output matrix C2 for the second system.
		C
		C     LDC2    INTEGER
		C             The leading dimension of array C2.
		C             LDC2 >= MAX(1,P2) if N2 > 0.
		C             LDC2 >= 1 if N2 = 0.
		C
		C     D2      (input) DOUBLE PRECISION array, dimension (LDD2,P1)
		C             The leading P2-by-P1 part of this array must contain the
		C             input/output matrix D2 for the second system.
		C
		C     LDD2    INTEGER
		C             The leading dimension of array D2.  LDD2 >= MAX(1,P2).
		C
		C     N       (output) INTEGER
		C             The number of state variables (N1 + N2) in the resulting
		C             system, i.e. the order of the matrix A, the number of rows
		C             of B and the number of columns of C.
		C
		C     A       (output) DOUBLE PRECISION array, dimension (LDA,N1+N2)
		C             The leading N-by-N part of this array contains the state
		C             transition matrix A for the cascaded system.
		C             If OVER = 'O', the array A can overlap A1, if UPLO = 'L',
		C             or A2, if UPLO = 'U'.
		C
		C     LDA     INTEGER
		C             The leading dimension of array A.  LDA >= MAX(1,N1+N2).
		C
		C     B       (output) DOUBLE PRECISION array, dimension (LDB,M1)
		C             The leading N-by-M1 part of this array contains the
		C             input/state matrix B for the cascaded system.
		C             If OVER = 'O', the array B can overlap B1, if UPLO = 'L',
		C             or B2, if UPLO = 'U'.
		C
		C     LDB     INTEGER
		C             The leading dimension of array B.  LDB >= MAX(1,N1+N2).
		C
		C     C       (output) DOUBLE PRECISION array, dimension (LDC,N1+N2)
		C             The leading P2-by-N part of this array contains the
		C             state/output matrix C for the cascaded system.
		C             If OVER = 'O', the array C can overlap C1, if UPLO = 'L',
		C             or C2, if UPLO = 'U'.
		C
		C     LDC     INTEGER
		C             The leading dimension of array C.
		C             LDC >= MAX(1,P2) if N1+N2 > 0.
		C             LDC >= 1 if N1+N2 = 0.
		C
		C     D       (output) DOUBLE PRECISION array, dimension (LDD,M1)
		C             The leading P2-by-M1 part of this array contains the
		C             input/output matrix D for the cascaded system.
		C             If OVER = 'O', the array D can overlap D1, if UPLO = 'L',
		C             or D2, if UPLO = 'U'.
		C
		C     LDD     INTEGER
		C             The leading dimension of array D.  LDD >= MAX(1,P2).
		C
		C     Workspace
		C
		C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
		C             The array DWORK is not referenced if OVER = 'N'.
		C
		C     LDWORK  INTEGER
		C             The length of the array DWORK.
		C             LDWORK >= MAX( 1, P1*MAX(N1, M1, N2, P2) ) if OVER = 'O'.
		C             LDWORK >= 1 if OVER = 'N'.
		C
		C     Error Indicator
		C
		C     INFO    INTEGER
		C             = 0:  successful exit;
		C             < 0:  if INFO = -i, the i-th argument had an illegal
		C                   value.
		C
		C     METHOD
		C
		C     After cascaded inter-connection of the two systems
		C
		C     X1'     = A1*X1 + B1*U
		C     V       = C1*X1 + D1*U
		C
		C     X2'     = A2*X2 + B2*V
		C     Y       = C2*X2 + D2*V
		C
		C     (where  '  denotes differentiation with respect to time)
		C
		C     the following state-space model will be obtained:
		C
		C     X'      = A*X + B*U
		C     Y       = C*X + D*U
		C
		C     where matrix  A  has the form   ( A1     0 ),
		C                                     ( B2*C1  A2)
		C
		C           matrix  B  has the form  (  B1   ),
		C                                    ( B2*D1 )
		C
		C           matrix  C  has the form  ( D2*C1  C2 ) and
		C
		C           matrix  D  has the form  ( D2*D1 ).
		C
		C     This form is returned by the routine when UPLO = 'L'.  Note that
		C     when A1 and A2 are block lower triangular, the resulting state
		C     matrix is also block lower triangular.
		C
		C     By applying a similarity transformation to the system above,
		C     using the matrix  ( 0  I ),  where  I  is the identity matrix of
		C                       ( J  0 )
		C     order  N2,  and  J  is the identity matrix of order  N1,  the
		C     system matrices become
		C
		C           A = ( A2  B2*C1 ),
		C               ( 0     A1  )
		C
		C           B = ( B2*D1 ),
		C               (  B1   )
		C
		C           C = ( C2  D2*C1 ) and
		C
		C           D = ( D2*D1 ).
		C
		C     This form is returned by the routine when UPLO = 'U'.  Note that
		C     when A1 and A2 are block upper triangular (for instance, in the
		C     real Schur form), the resulting state matrix is also block upper
		C     triangular.
		C
		C     REFERENCES
		C
		C     None
		C
		C     NUMERICAL ASPECTS
		C
		C     The algorithm requires P1*(N1+M1)*(N2+P2) operations.
		C
		C     CONTRIBUTORS
		C
		C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, and
		C                  A. Varga, German Aerospace Research Establishment,
		C                  Oberpfaffenhofen, Germany, Nov. 1996.
		C     Supersedes Release 2.0 routine AB05AD by C.J.Benson, Kingston
		C     Polytechnic, United Kingdom, January 1982.
		C
		C     REVISIONS
		C
		C     V. Sima, Research Institute for Informatics, Bucharest, July 2003,
		C     Feb. 2004.
		C
		C     KEYWORDS
		C
		C     Cascade control, continuous-time system, multivariable
		C     system, state-space model, state-space representation.
		*/

		void AB05MD(unsigned char *, unsigned char *, int *, int *, int *, int *, int *, double *,
		            int *, double *, int *, double *, int *, double *, int *, double *, int *, double *,
					int *, double *, int *, double *, int *, int *, double *, int *, double *, int *,
					double *, int *, double *, int *, double *, int *, int *);

		/*
		SLICOT RELEASE 5.5.
		C
		C     Copyright (c) 2002-2012 NICONET e.V.
		C
		C     PURPOSE
		C
		C     To obtain the state-space model (A,B,C,D) for the feedback
		C     inter-connection of two systems, each given in state-space form.
		C
		C     ARGUMENTS
		C
		C     Mode Parameters
		C
		C     OVER    CHARACTER*1
		C             Indicates whether the user wishes to overlap pairs of
		C             arrays, as follows:
		C             = 'N':  Do not overlap;
		C             = 'O':  Overlap pairs of arrays: A1 and A, B1 and B,
		C                     C1 and C, and D1 and D, i.e. the same name is
		C                     effectively used for each pair (for all pairs)
		C                     in the routine call.  In this case, setting
		C                     LDA1 = LDA, LDB1 = LDB, LDC1 = LDC, and LDD1 = LDD
		C                     will give maximum efficiency.
		C
		C     Input/Output Parameters
		C
		C     N1      (input) INTEGER
		C             The number of state variables in the first system, i.e.
		C             the order of the matrix A1.  N1 >= 0.
		C
		C     M1      (input) INTEGER
		C             The number of input variables for the first system and the
		C             number of output variables from the second system.
		C             M1 >= 0.
		C
		C     P1      (input) INTEGER
		C             The number of output variables from the first system and
		C             the number of input variables for the second system.
		C             P1 >= 0.
		C
		C     N2      (input) INTEGER
		C             The number of state variables in the second system, i.e.
		C             the order of the matrix A2.  N2 >= 0.
		C
		C     ALPHA   (input) DOUBLE PRECISION
		C             A coefficient multiplying the transfer-function matrix
		C             (or the output equation) of the second system.
		C             ALPHA = +1 corresponds to positive feedback, and
		C             ALPHA = -1 corresponds to negative feedback.
		C
		C     A1      (input) DOUBLE PRECISION array, dimension (LDA1,N1)
		C             The leading N1-by-N1 part of this array must contain the
		C             state transition matrix A1 for the first system.
		C
		C     LDA1    INTEGER
		C             The leading dimension of array A1.  LDA1 >= MAX(1,N1).
		C
		C     B1      (input) DOUBLE PRECISION array, dimension (LDB1,M1)
		C             The leading N1-by-M1 part of this array must contain the
		C             input/state matrix B1 for the first system.
		C
		C     LDB1    INTEGER
		C             The leading dimension of array B1.  LDB1 >= MAX(1,N1).
		C
		C     C1      (input) DOUBLE PRECISION array, dimension (LDC1,N1)
		C             The leading P1-by-N1 part of this array must contain the
		C             state/output matrix C1 for the first system.
		C
		C     LDC1    INTEGER
		C             The leading dimension of array C1.
		C             LDC1 >= MAX(1,P1) if N1 > 0.
		C             LDC1 >= 1 if N1 = 0.
		C
		C     D1      (input) DOUBLE PRECISION array, dimension (LDD1,M1)
		C             The leading P1-by-M1 part of this array must contain the
		C             input/output matrix D1 for the first system.
		C
		C     LDD1    INTEGER
		C             The leading dimension of array D1.  LDD1 >= MAX(1,P1).
		C
		C     A2      (input) DOUBLE PRECISION array, dimension (LDA2,N2)
		C             The leading N2-by-N2 part of this array must contain the
		C             state transition matrix A2 for the second system.
		C
		C     LDA2    INTEGER
		C             The leading dimension of array A2.  LDA2 >= MAX(1,N2).
		C
		C     B2      (input) DOUBLE PRECISION array, dimension (LDB2,P1)
		C             The leading N2-by-P1 part of this array must contain the
		C             input/state matrix B2 for the second system.
		C
		C     LDB2    INTEGER
		C             The leading dimension of array B2.  LDB2 >= MAX(1,N2).
		C
		C     C2      (input) DOUBLE PRECISION array, dimension (LDC2,N2)
		C             The leading M1-by-N2 part of this array must contain the
		C             state/output matrix C2 for the second system.
		C
		C     LDC2    INTEGER
		C             The leading dimension of array C2.
		C             LDC2 >= MAX(1,M1) if N2 > 0.
		C             LDC2 >= 1 if N2 = 0.
		C
		C     D2      (input) DOUBLE PRECISION array, dimension (LDD2,P1)
		C             The leading M1-by-P1 part of this array must contain the
		C             input/output matrix D2 for the second system.
		C
		C     LDD2    INTEGER
		C             The leading dimension of array D2.  LDD2 >= MAX(1,M1).
		C
		C     N       (output) INTEGER
		C             The number of state variables (N1 + N2) in the connected
		C             system, i.e. the order of the matrix A, the number of rows
		C             of B and the number of columns of C.
		C
		C     A       (output) DOUBLE PRECISION array, dimension (LDA,N1+N2)
		C             The leading N-by-N part of this array contains the state
		C             transition matrix A for the connected system.
		C             The array A can overlap A1 if OVER = 'O'.
		C
		C     LDA     INTEGER
		C             The leading dimension of array A.  LDA >= MAX(1,N1+N2).
		C
		C     B       (output) DOUBLE PRECISION array, dimension (LDB,M1)
		C             The leading N-by-M1 part of this array contains the
		C             input/state matrix B for the connected system.
		C             The array B can overlap B1 if OVER = 'O'.
		C
		C     LDB     INTEGER
		C             The leading dimension of array B.  LDB >= MAX(1,N1+N2).
		C
		C     C       (output) DOUBLE PRECISION array, dimension (LDC,N1+N2)
		C             The leading P1-by-N part of this array contains the
		C             state/output matrix C for the connected system.
		C             The array C can overlap C1 if OVER = 'O'.
		C
		C     LDC     INTEGER
		C             The leading dimension of array C.
		C             LDC >= MAX(1,P1) if N1+N2 > 0.
		C             LDC >= 1 if N1+N2 = 0.
		C
		C     D       (output) DOUBLE PRECISION array, dimension (LDD,M1)
		C             The leading P1-by-M1 part of this array contains the
		C             input/output matrix D for the connected system.
		C             The array D can overlap D1 if OVER = 'O'.
		C
		C     LDD     INTEGER
		C             The leading dimension of array D.  LDD >= MAX(1,P1).
		C
		C     Workspace
		C
		C     IWORK   INTEGER array, dimension (P1)
		C
		C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
		C
		C     LDWORK  INTEGER
		C             The length of the array DWORK.             If OVER = 'N',
		C             LDWORK >= MAX(1, P1*P1, M1*M1, N1*P1), and if OVER = 'O',
		C             LDWORK >= MAX(1, N1*P1 + MAX( P1*P1, M1*M1, N1*P1) ),
		C                                                        if M1 <= N*N2;
		C             LDWORK >= MAX(1, N1*P1 + MAX( P1*P1, M1*(M1+1), N1*P1) ),
		C                                                        if M1 >  N*N2.
		C
		C     Error Indicator
		C
		C     INFO    INTEGER
		C             = 0:  successful exit;
		C             < 0:  if INFO = -i, the i-th argument had an illegal
		C                   value.
		C             > 0:  if INFO = i, 1 <= i <= P1, the system is not
		C                   completely controllable. That is, the matrix
		C                   (I + ALPHA*D1*D2) is exactly singular (the element
		C                   U(i,i) of the upper triangular factor of LU
		C                   factorization is exactly zero), possibly due to
		C                   rounding errors.
		C
		C     METHOD
		C
		C     After feedback inter-connection of the two systems,
		C
		C     X1'     = A1*X1 + B1*U1
		C     Y1      = C1*X1 + D1*U1
		C
		C     X2'     = A2*X2 + B2*U2
		C     Y2      = C2*X2 + D2*U2
		C
		C     (where  '  denotes differentiation with respect to time)
		C
		C     the following state-space model will be obtained:
		C
		C     X'      = A*X  +  B*U
		C     Y       = C*X  +  D*U
		C
		C     where       U = U1 + alpha*Y2,    X  =  ( X1 ),
		C                 Y = Y1 = U2,                ( X2 )
		C
		C     matrix  A  has the form
		C
		C     ( A1  -  alpha*B1*E12*D2*C1       -  alpha*B1*E12*C2    ),
		C     (        B2*E21*C1            A2  -  alpha*B2*E21*D1*C2 )
		C
		C     matrix  B  has the form
		C
		C     (  B1*E12    ),
		C     (  B2*E21*D1 )
		C
		C     matrix  C  has the form
		C
		C     (  E21*C1     -  alpha*E21*D1*C2 ),
		C
		C     matrix D  has the form
		C
		C     (  E21*D1 ),
		C
		C     E21  =  ( I + alpha*D1*D2 )-INVERSE and
		C     E12  =  ( I + alpha*D2*D1 )-INVERSE = I - alpha*D2*E21*D1.
		C
		C     Taking N1 = 0 and/or N2 = 0 on the routine call will solve the
		C     constant plant and/or constant feedback cases.
		C
		C     REFERENCES
		C
		C     None
		C
		C     NUMERICAL ASPECTS
		C
		C     None
		C
		C     CONTRIBUTOR
		C
		C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996.
		C     Supersedes Release 2.0 routine AB05BD by C.J.Benson, Kingston
		C     Polytechnic, United Kingdom, January 1982.
		C
		C     REVISIONS
		C
		C     V. Sima, Research Institute for Informatics, Bucharest, July 2003,
		C     Feb. 2004.
		C
		C     KEYWORDS
		C
		C     Continuous-time system, multivariable system, state-space model,
		C     state-space representation.
		*/

		void AB05ND(unsigned char *, int *, int *, int *, int *, double *, double *, int *,
		           double *, int *, double *, int *, double *, int *, double *, int *,
				   double *, int *, double *, int *, double *, int *, int *, double *,
				   int *, double *, int *, double *, int *, double *, int *, int *, double *, int *, int *);

		/*
		C     SLICOT RELEASE 5.5.
		C
		C     Copyright (c) 2002-2012 NICONET e.V.
		C
		C     PURPOSE
		C
		C     To obtain the state-space model (A,B,C,D) for rowwise
		C     concatenation (parallel inter-connection on outputs, with separate
		C     inputs) of two systems, each given in state-space form.
		C
		C     ARGUMENTS
		C
		C     Mode Parameters
		C
		C     OVER    CHARACTER*1
		C             Indicates whether the user wishes to overlap pairs of
		C             arrays, as follows:
		C             = 'N':  Do not overlap;
		C             = 'O':  Overlap pairs of arrays: A1 and A, B1 and B,
		C                     C1 and C, and D1 and D, i.e. the same name is
		C                     effectively used for each pair (for all pairs)
		C                     in the routine call.  In this case, setting
		C                     LDA1 = LDA, LDB1 = LDB, LDC1 = LDC, and LDD1 = LDD
		C                     will give maximum efficiency.
		C
		C     Input/Output Parameters
		C
		C     N1      (input) INTEGER
		C             The number of state variables in the first system, i.e.
		C             the order of the matrix A1.  N1 >= 0.
		C
		C     M1      (input) INTEGER
		C             The number of input variables for the first system.
		C             M1 >= 0.
		C
		C     P1      (input) INTEGER
		C             The number of output variables from each system.  P1 >= 0.
		C
		C     N2      (input) INTEGER
		C             The number of state variables in the second system, i.e.
		C             the order of the matrix A2.  N2 >= 0.
		C
		C     M2      (input) INTEGER
		C             The number of input variables for the second system.
		C             M2 >= 0.
		C
		C     ALPHA   (input) DOUBLE PRECISION
		C             A coefficient multiplying the transfer-function matrix
		C             (or the output equation) of the second system.
		C
		C     A1      (input) DOUBLE PRECISION array, dimension (LDA1,N1)
		C             The leading N1-by-N1 part of this array must contain the
		C             state transition matrix A1 for the first system.
		C
		C     LDA1    INTEGER
		C             The leading dimension of array A1.  LDA1 >= MAX(1,N1).
		C
		C     B1      (input) DOUBLE PRECISION array, dimension (LDB1,M1)
		C             The leading N1-by-M1 part of this array must contain the
		C             input/state matrix B1 for the first system.
		C
		C     LDB1    INTEGER
		C             The leading dimension of array B1.  LDB1 >= MAX(1,N1).
		C
		C     C1      (input) DOUBLE PRECISION array, dimension (LDC1,N1)
		C             The leading P1-by-N1 part of this array must contain the
		C             state/output matrix C1 for the first system.
		C
		C     LDC1    INTEGER
		C             The leading dimension of array C1.
		C             LDC1 >= MAX(1,P1) if N1 > 0.
		C             LDC1 >= 1 if N1 = 0.
		C
		C     D1      (input) DOUBLE PRECISION array, dimension (LDD1,M1)
		C             The leading P1-by-M1 part of this array must contain the
		C             input/output matrix D1 for the first system.
		C
		C     LDD1    INTEGER
		C             The leading dimension of array D1.  LDD1 >= MAX(1,P1).
		C
		C     A2      (input) DOUBLE PRECISION array, dimension (LDA2,N2)
		C             The leading N2-by-N2 part of this array must contain the
		C             state transition matrix A2 for the second system.
		C
		C     LDA2    INTEGER
		C             The leading dimension of array A2.  LDA2 >= MAX(1,N2).
		C
		C     B2      (input) DOUBLE PRECISION array, dimension (LDB2,M2)
		C             The leading N2-by-M2 part of this array must contain the
		C             input/state matrix B2 for the second system.
		C
		C     LDB2    INTEGER
		C             The leading dimension of array B2.  LDB2 >= MAX(1,N2).
		C
		C     C2      (input) DOUBLE PRECISION array, dimension (LDC2,N2)
		C             The leading P1-by-N2 part of this array must contain the
		C             state/output matrix C2 for the second system.
		C
		C     LDC2    INTEGER
		C             The leading dimension of array C2.
		C             LDC2 >= MAX(1,P1) if N2 > 0.
		C             LDC2 >= 1 if N2 = 0.
		C
		C     D2      (input) DOUBLE PRECISION array, dimension (LDD2,M2)
		C             The leading P1-by-M2 part of this array must contain the
		C             input/output matrix D2 for the second system.
		C
		C     LDD2    INTEGER
		C             The leading dimension of array D2.  LDD2 >= MAX(1,P1).
		C
		C     N       (output) INTEGER
		C             The number of state variables (N1 + N2) in the connected
		C             system, i.e. the order of the matrix A, the number of rows
		C             of B and the number of columns of C.
		C
		C     M       (output) INTEGER
		C             The number of input variables (M1 + M2) for the connected
		C             system, i.e. the number of columns of B and D.
		C
		C     A       (output) DOUBLE PRECISION array, dimension (LDA,N1+N2)
		C             The leading N-by-N part of this array contains the state
		C             transition matrix A for the connected system.
		C             The array A can overlap A1 if OVER = 'O'.
		C
		C     LDA     INTEGER
		C             The leading dimension of array A.  LDA >= MAX(1,N1+N2).
		C
		C     B       (output) DOUBLE PRECISION array, dimension (LDB,M1+M2)
		C             The leading N-by-M part of this array contains the
		C             input/state matrix B for the connected system.
		C             The array B can overlap B1 if OVER = 'O'.
		C
		C     LDB     INTEGER
		C             The leading dimension of array B.  LDB >= MAX(1,N1+N2).
		C
		C     C       (output) DOUBLE PRECISION array, dimension (LDC,N1+N2)
		C             The leading P1-by-N part of this array contains the
		C             state/output matrix C for the connected system.
		C             The array C can overlap C1 if OVER = 'O'.
		C
		C     LDC     INTEGER
		C             The leading dimension of array C.
		C             LDC >= MAX(1,P1) if N1+N2 > 0.
		C             LDC >= 1 if N1+N2 = 0.
		C
		C     D       (output) DOUBLE PRECISION array, dimension (LDD,M1+M2)
		C             The leading P1-by-M part of this array contains the
		C             input/output matrix D for the connected system.
		C             The array D can overlap D1 if OVER = 'O'.
		C
		C     LDD     INTEGER
		C             The leading dimension of array D.  LDD >= MAX(1,P1).
		C
		C     Error Indicator
		C
		C     INFO    INTEGER
		C             = 0:  successful exit;
		C             < 0:  if INFO = -i, the i-th argument had an illegal
		C                   value.
		C
		C     METHOD
		C
		C     After rowwise concatenation (parallel inter-connection with
		C     separate inputs) of the two systems,
		C
		C     X1'     = A1*X1 + B1*U
		C     Y1      = C1*X1 + D1*U
		C
		C     X2'     = A2*X2 + B2*V
		C     Y2      = C2*X2 + D2*V
		C
		C     (where  '  denotes differentiation with respect to time),
		C
		C     with the output equation for the second system multiplied by a
		C     scalar alpha, the following state-space model will be obtained:
		C
		C     X'      = A*X + B*(U)
		C                       (V)
		C
		C     Y       = C*X + D*(U)
		C                       (V)
		C
		C     where matrix  A  has the form    ( A1   0  ),
		C                                      ( 0    A2 )
		C
		C           matrix  B  has the form    ( B1   0  ),
		C                                      ( 0    B2 )
		C
		C           matrix  C  has the form    ( C1   alpha*C2 ) and
		C
		C           matrix  D  has the form    ( D1   alpha*D2 ).
		C
		C     REFERENCES
		C
		C     None
		C
		C     NUMERICAL ASPECTS
		C
		C     None
		C
		C     CONTRIBUTOR
		C
		C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Oct. 1996.
		C     Supersedes Release 2.0 routine AB05CD by C.J.Benson, Kingston
		C     Polytechnic, United Kingdom, January 1982.
		C
		C     REVISIONS
		C
		C     V. Sima, Research Institute for Informatics, Bucharest, July 2003,
		C     Feb. 2004.
		C
		C     KEYWORDS
		C
		C     Continuous-time system, multivariable system, state-space model,
		C     state-space representation.
		*/

		void AB05OD(unsigned char *, int *, int *, int *, int *, int *, double *, double *, int *,
		            double *, int *, double *, int *, double *, int *, double *, int *, double *,
					int *, double *, int *, double *, int *, int *, int *, double *, int *, double *,
					int *, double *, int *, double *, int *, int *);

		/*
		SLICOT RELEASE 5.5.
		C
		C     Copyright (c) 2002-2012 NICONET e.V.
		C
		C     PURPOSE
		C
		C     To compute the state-space model G = (A,B,C,D) corresponding to
		C     the sum G = G1 + alpha*G2, where G1 = (A1,B1,C1,D1) and
		C     G2 = (A2,B2,C2,D2).  G, G1, and G2 are the transfer-function
		C     matrices of the corresponding state-space models.
		C
		C     ARGUMENTS
		C
		C     Mode Parameters
		C
		C     OVER    CHARACTER*1
		C             Indicates whether the user wishes to overlap pairs of
		C             arrays, as follows:
		C             = 'N':  Do not overlap;
		C             = 'O':  Overlap pairs of arrays: A1 and A, B1 and B,
		C                     C1 and C, and D1 and D, i.e. the same name is
		C                     effectively used for each pair (for all pairs)
		C                     in the routine call.  In this case, setting
		C                     LDA1 = LDA, LDB1 = LDB, LDC1 = LDC, and LDD1 = LDD
		C                     will give maximum efficiency.
		C
		C     Input/Output Parameters
		C
		C     N1      (input) INTEGER
		C             The number of state variables in the first system, i.e.
		C             the order of the matrix A1, the number of rows of B1 and
		C             the number of columns of C1.  N1 >= 0.
		C
		C     M       (input) INTEGER
		C             The number of input variables of the two systems, i.e. the
		C             number of columns of matrices B1, D1, B2 and D2.  M >= 0.
		C
		C     P       (input) INTEGER
		C             The number of output variables of the two systems, i.e.
		C             the number of rows of matrices C1, D1, C2 and D2.  P >= 0.
		C
		C     N2      (input) INTEGER
		C             The number of state variables in the second system, i.e.
		C             the order of the matrix A2, the number of rows of B2 and
		C             the number of columns of C2.  N2 >= 0.
		C
		C     ALPHA   (input) DOUBLE PRECISION
		C             The coefficient multiplying G2.
		C
		C     A1      (input) DOUBLE PRECISION array, dimension (LDA1,N1)
		C             The leading N1-by-N1 part of this array must contain the
		C             state transition matrix A1 for the first system.
		C
		C     LDA1    INTEGER
		C             The leading dimension of array A1.  LDA1 >= MAX(1,N1).
		C
		C     B1      (input) DOUBLE PRECISION array, dimension (LDB1,M)
		C             The leading N1-by-M part of this array must contain the
		C             input/state matrix B1 for the first system.
		C
		C     LDB1    INTEGER
		C             The leading dimension of array B1.  LDB1 >= MAX(1,N1).
		C
		C     C1      (input) DOUBLE PRECISION array, dimension (LDC1,N1)
		C             The leading P-by-N1 part of this array must contain the
		C             state/output matrix C1 for the first system.
		C
		C     LDC1    INTEGER
		C             The leading dimension of array C1.
		C             LDC1 >= MAX(1,P) if N1 > 0.
		C             LDC1 >= 1 if N1 = 0.
		C
		C     D1      (input) DOUBLE PRECISION array, dimension (LDD1,M)
		C             The leading P-by-M part of this array must contain the
		C             input/output matrix D1 for the first system.
		C
		C     LDD1    INTEGER
		C             The leading dimension of array D1.  LDD1 >= MAX(1,P).
		C
		C     A2      (input) DOUBLE PRECISION array, dimension (LDA2,N2)
		C             The leading N2-by-N2 part of this array must contain the
		C             state transition matrix A2 for the second system.
		C
		C     LDA2    INTEGER
		C             The leading dimension of array A2.  LDA2 >= MAX(1,N2).
		C
		C     B2      (input) DOUBLE PRECISION array, dimension (LDB2,M)
		C             The leading N2-by-M part of this array must contain the
		C             input/state matrix B2 for the second system.
		C
		C     LDB2    INTEGER
		C             The leading dimension of array B2.  LDB2 >= MAX(1,N2).
		C
		C     C2      (input) DOUBLE PRECISION array, dimension (LDC2,N2)
		C             The leading P-by-N2 part of this array must contain the
		C             state/output matrix C2 for the second system.
		C
		C     LDC2    INTEGER
		C             The leading dimension of array C2.
		C             LDC2 >= MAX(1,P) if N2 > 0.
		C             LDC2 >= 1 if N2 = 0.
		C
		C     D2      (input) DOUBLE PRECISION array, dimension (LDD2,M)
		C             The leading P-by-M part of this array must contain the
		C             input/output matrix D2 for the second system.
		C
		C     LDD2    INTEGER
		C             The leading dimension of array D2.  LDD2 >= MAX(1,P).
		C
		C     N       (output) INTEGER
		C             The number of state variables (N1 + N2) in the resulting
		C             system, i.e. the order of the matrix A, the number of rows
		C             of B and the number of columns of C.
		C
		C     A       (output) DOUBLE PRECISION array, dimension (LDA,N1+N2)
		C             The leading N-by-N part of this array contains the state
		C             transition matrix A for the resulting system.
		C             The array A can overlap A1 if OVER = 'O'.
		C
		C     LDA     INTEGER
		C             The leading dimension of array A.  LDA >= MAX(1,N1+N2).
		C
		C     B       (output) DOUBLE PRECISION array, dimension (LDB,M)
		C             The leading N-by-M part of this array contains the
		C             input/state matrix B for the resulting system.
		C             The array B can overlap B1 if OVER = 'O'.
		C
		C     LDB     INTEGER
		C             The leading dimension of array B.  LDB >= MAX(1,N1+N2).
		C
		C     C       (output) DOUBLE PRECISION array, dimension (LDC,N1+N2)
		C             The leading P-by-N part of this array contains the
		C             state/output matrix C for the resulting system.
		C             The array C can overlap C1 if OVER = 'O'.
		C
		C     LDC     INTEGER
		C             The leading dimension of array C.
		C             LDC >= MAX(1,P) if N1+N2 > 0.
		C             LDC >= 1 if N1+N2 = 0.
		C
		C     D       (output) DOUBLE PRECISION array, dimension (LDD,M)
		C             The leading P-by-M part of this array contains the
		C             input/output matrix D for the resulting system.
		C             The array D can overlap D1 if OVER = 'O'.
		C
		C     LDD     INTEGER
		C             The leading dimension of array D.  LDD >= MAX(1,P).
		C
		C     Error Indicator
		C
		C     INFO    INTEGER
		C             = 0:  successful exit;
		C             < 0:  if INFO = -i, the i-th argument had an illegal
		C                   value.
		C
		C     METHOD
		C
		C     The matrices of the resulting systems are determined as:
		C
		C           ( A1   0  )             ( B1 )
		C       A = (         ) ,       B = (    ) ,
		C           ( 0    A2 )             ( B2 )
		C
		C       C = ( C1  alpha*C2 ) ,  D = D1 + alpha*D2 .
		C
		C     REFERENCES
		C
		C     None
		C
		C     NUMERICAL ASPECTS
		C
		C     None
		C
		C     CONTRIBUTORS
		C
		C     A. Varga, German Aerospace Research Establishment,
		C     Oberpfaffenhofen, Germany, and V. Sima, Katholieke Univ. Leuven,
		C     Belgium, Nov. 1996.
		C
		C     REVISIONS
		C
		C     V. Sima, Research Institute for Informatics, Bucharest, July 2003,
		C     Feb. 2004.
		C
		C     KEYWORDS
		C
		C     Multivariable system, state-space model, state-space
		C     representation.
		*/

		void AB05PD(unsigned char *, int *, int *, int *, int *, double *, double *, int *,
		            double *, int *, double *, int *, double *, int *, double *, int *, 
					double *, int *, double *, int *, double *, int *, int *, double *,
					int *, double *, int *, double *, int *, double *, int *, int *);

		/*
		SLICOT RELEASE 5.5.
		C
		C     Copyright (c) 2002-2012 NICONET e.V.
		C
		C     PURPOSE
		C
		C     To append two systems G1 and G2 in state-space form together.
		C     If G1 = (A1,B1,C1,D1) and G2 = (A2,B2,C2,D2) are the state-space
		C     models of the given two systems having the transfer-function
		C     matrices G1 and G2, respectively, this subroutine constructs the
		C     state-space model G = (A,B,C,D) which corresponds to the
		C     transfer-function matrix
		C
		C                           ( G1 0  )
		C                       G = (       )
		C                           ( 0  G2 )
		C
		C     ARGUMENTS
		C
		C     Mode Parameters
		C
		C     OVER    CHARACTER*1
		C             Indicates whether the user wishes to overlap pairs of
		C             arrays, as follows:
		C             = 'N':  Do not overlap;
		C             = 'O':  Overlap pairs of arrays: A1 and A, B1 and B,
		C                     C1 and C, and D1 and D, i.e. the same name is
		C                     effectively used for each pair (for all pairs)
		C                     in the routine call.  In this case, setting
		C                     LDA1 = LDA, LDB1 = LDB, LDC1 = LDC, and LDD1 = LDD
		C                     will give maximum efficiency.
		C
		C     Input/Output Parameters
		C
		C     N1      (input) INTEGER
		C             The number of state variables in the first system, i.e.
		C             the order of the matrix A1, the number of rows of B1 and
		C             the number of columns of C1.  N1 >= 0.
		C
		C     M1      (input) INTEGER
		C             The number of input variables in the first system, i.e.
		C             the number of columns of matrices B1 and D1.  M1 >= 0.
		C
		C     P1      (input) INTEGER
		C             The number of output variables in the first system, i.e.
		C             the number of rows of matrices C1 and D1.  P1 >= 0.
		C
		C     N2      (input) INTEGER
		C             The number of state variables in the second system, i.e.
		C             the order of the matrix A2, the number of rows of B2 and
		C             the number of columns of C2.  N2 >= 0.
		C
		C     M2      (input) INTEGER
		C             The number of input variables in the second system, i.e.
		C             the number of columns of matrices B2 and D2.  M2 >= 0.
		C
		C     P2      (input) INTEGER
		C             The number of output variables in the second system, i.e.
		C             the number of rows of matrices C2 and D2.  P2 >= 0.
		C
		C     A1      (input) DOUBLE PRECISION array, dimension (LDA1,N1)
		C             The leading N1-by-N1 part of this array must contain the
		C             state transition matrix A1 for the first system.
		C
		C     LDA1    INTEGER
		C             The leading dimension of array A1.  LDA1 >= MAX(1,N1).
		C
		C     B1      (input) DOUBLE PRECISION array, dimension (LDB1,M1)
		C             The leading N1-by-M1 part of this array must contain the
		C             input/state matrix B1 for the first system.
		C
		C     LDB1    INTEGER
		C             The leading dimension of array B1.  LDB1 >= MAX(1,N1).
		C
		C     C1      (input) DOUBLE PRECISION array, dimension (LDC1,N1)
		C             The leading P1-by-N1 part of this array must contain the
		C             state/output matrix C1 for the first system.
		C
		C     LDC1    INTEGER
		C             The leading dimension of array C1.
		C             LDC1 >= MAX(1,P1) if N1 > 0.
		C             LDC1 >= 1 if N1 = 0.
		C
		C     D1      (input) DOUBLE PRECISION array, dimension (LDD1,M1)
		C             The leading P1-by-M1 part of this array must contain the
		C             input/output matrix D1 for the first system.
		C
		C     LDD1    INTEGER
		C             The leading dimension of array D1.  LDD1 >= MAX(1,P1).
		C
		C     A2      (input) DOUBLE PRECISION array, dimension (LDA2,N2)
		C             The leading N2-by-N2 part of this array must contain the
		C             state transition matrix A2 for the second system.
		C
		C     LDA2    INTEGER
		C             The leading dimension of array A2.  LDA2 >= MAX(1,N2).
		C
		C     B2      (input) DOUBLE PRECISION array, dimension (LDB2,M2)
		C             The leading N2-by-M2 part of this array must contain the
		C             input/state matrix B2 for the second system.
		C
		C     LDB2    INTEGER
		C             The leading dimension of array B2.  LDB2 >= MAX(1,N2).
		C
		C     C2      (input) DOUBLE PRECISION array, dimension (LDC2,N2)
		C             The leading P2-by-N2 part of this array must contain the
		C             state/output matrix C2 for the second system.
		C
		C     LDC2    INTEGER
		C             The leading dimension of array C2.
		C             LDC2 >= MAX(1,P2) if N2 > 0.
		C             LDC2 >= 1 if N2 = 0.
		C
		C     D2      (input) DOUBLE PRECISION array, dimension (LDD2,M2)
		C             The leading P2-by-M2 part of this array must contain the
		C             input/output matrix D2 for the second system.
		C
		C     LDD2    INTEGER
		C             The leading dimension of array D2.  LDD2 >= MAX(1,P2).
		C
		C     N       (output) INTEGER
		C             The number of state variables (N1 + N2) in the resulting
		C             system, i.e. the order of the matrix A, the number of rows
		C             of B and the number of columns of C.
		C
		C     M       (output) INTEGER
		C             The number of input variables (M1 + M2) in the resulting
		C             system, i.e. the number of columns of B and D.
		C
		C     P       (output) INTEGER
		C             The number of output variables (P1 + P2) of the resulting
		C             system, i.e. the number of rows of C and D.
		C
		C     A       (output) DOUBLE PRECISION array, dimension (LDA,N1+N2)
		C             The leading N-by-N part of this array contains the state
		C             transition matrix A for the resulting system.
		C             The array A can overlap A1 if OVER = 'O'.
		C
		C     LDA     INTEGER
		C             The leading dimension of array A.  LDA >= MAX(1,N1+N2).
		C
		C     B       (output) DOUBLE PRECISION array, dimension (LDB,M1+M2)
		C             The leading N-by-M part of this array contains the
		C             input/state matrix B for the resulting system.
		C             The array B can overlap B1 if OVER = 'O'.
		C
		C     LDB     INTEGER
		C             The leading dimension of array B.  LDB >= MAX(1,N1+N2).
		C
		C     C       (output) DOUBLE PRECISION array, dimension (LDC,N1+N2)
		C             The leading P-by-N part of this array contains the
		C             state/output matrix C for the resulting system.
		C             The array C can overlap C1 if OVER = 'O'.
		C
		C     LDC     INTEGER
		C             The leading dimension of array C.
		C             LDC >= MAX(1,P1+P2) if N1+N2 > 0.
		C             LDC >= 1 if N1+N2 = 0.
		C
		C     D       (output) DOUBLE PRECISION array, dimension (LDD,M1+M2)
		C             The leading P-by-M part of this array contains the
		C             input/output matrix D for the resulting system.
		C             The array D can overlap D1 if OVER = 'O'.
		C
		C     LDD     INTEGER
		C             The leading dimension of array D.  LDD >= MAX(1,P1+P2).
		C
		C     Error Indicator
		C
		C     INFO    INTEGER
		C             = 0:  successful exit;
		C             < 0:  if INFO = -i, the i-th argument had an illegal
		C                   value.
		C
		C     METHOD
		C
		C     The matrices of the resulting systems are determined as:
		C
		C           ( A1   0  )         ( B1  0  )
		C       A = (         ) ,   B = (        ) ,
		C           ( 0    A2 )         ( 0   B2 )
		C
		C           ( C1   0  )         ( D1  0  )
		C       C = (         ) ,   D = (        ) .
		C           ( 0    C2 )         ( 0   D2 )
		C
		C     REFERENCES
		C
		C     None
		C
		C     CONTRIBUTORS
		C
		C     A. Varga, German Aerospace Research Establishment,
		C     Oberpfaffenhofen, Germany, and V. Sima, Katholieke Univ. Leuven,
		C     Belgium, Nov. 1996.
		C
		C     REVISIONS
		C
		C     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2004.
		C
		C     KEYWORDS
		C
		C     Multivariable system, state-space model, state-space
		C     representation.
		*/

		void AB05QD(unsigned char *, int *, int *, int *, int *, int *, int *, double *, int *,
		            double *, int *, double *, int *, double *, int *, double *, int *, double *,
					int *, double *, int *, double *, int *, int *, int *, int *, double *, int *,
					double *, int *, double *, int *, double *, int *, int *);

		/*
		SLICOT RELEASE 5.5.
		C
		C     Copyright (c) 2002-2012 NICONET e.V.
		C
		C     PURPOSE
		C
		C     To construct for a given state space system (A,B,C,D) the closed-
		C     loop system (Ac,Bc,Cc,Dc) corresponding to the mixed output and
		C     state feedback control law
		C
		C          u = alpha*F*y + beta*K*x + G*v
		C          z = H*y.
		C
		C     ARGUMENTS
		C
		C     Mode Parameters
		C
		C     FBTYPE  CHARACTER*1
		C             Specifies the type of the feedback law as follows:
		C             = 'I':  Unitary output feedback (F = I);
		C             = 'O':  General output feedback.
		C
		C     JOBD    CHARACTER*1
		C             Specifies whether or not a non-zero matrix D appears
		C             in the given state space model:
		C             = 'D':  D is present;
		C             = 'Z':  D is assumed a zero matrix.
		C
		C     Input/Output Parameters
		C
		C     N       (input) INTEGER
		C             The dimension of state vector x, i.e. the order of the
		C             matrix A, the number of rows of B and the number of
		C             columns of C.  N >= 0.
		C
		C     M       (input) INTEGER
		C             The dimension of input vector u, i.e. the number of
		C             columns of matrices B and D, and the number of rows of F.
		C             M >= 0.
		C
		C     P       (input) INTEGER
		C             The dimension of output vector y, i.e. the number of rows
		C             of matrices C and D, and the number of columns of F.
		C             P >= 0 and P = M if FBTYPE = 'I'.
		C
		C     MV      (input) INTEGER
		C             The dimension of the new input vector v, i.e. the number
		C             of columns of matrix G.  MV >= 0.
		C
		C     PZ      (input) INTEGER.
		C             The dimension of the new output vector z, i.e. the number
		C             of rows of matrix H.  PZ >= 0.
		C
		C     ALPHA   (input) DOUBLE PRECISION
		C             The coefficient alpha in the output feedback law.
		C
		C     BETA    (input) DOUBLE PRECISION.
		C             The coefficient beta in the state feedback law.
		C
		C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
		C             On entry, the leading N-by-N part of this array must
		C             contain the system state transition matrix A.
		C             On exit, the leading N-by-N part of this array contains
		C             the state matrix Ac of the closed-loop system.
		C
		C     LDA     INTEGER
		C             The leading dimension of array A.  LDA >= MAX(1,N).
		C
		C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
		C             On entry, the leading N-by-M part of this array must
		C             contain the system input matrix B.
		C             On exit, the leading N-by-M part of this array contains
		C             the intermediary input matrix B1 (see METHOD).
		C
		C     LDB     INTEGER
		C             The leading dimension of array B.  LDB >= MAX(1,N).
		C
		C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
		C             On entry, the leading P-by-N part of this array must
		C             contain the system output matrix C.
		C             On exit, the leading P-by-N part of this array contains
		C             the intermediary output matrix C1+BETA*D1*K (see METHOD).
		C
		C     LDC     INTEGER
		C             The leading dimension of array C.
		C             LDC >= MAX(1,P) if N > 0.
		C             LDC >= 1 if N = 0.
		C
		C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
		C             On entry, if JOBD = 'D', the leading P-by-M part of this
		C             array must contain the system direct input/output
		C             transmission matrix D.
		C             On exit, the leading P-by-M part of this array contains
		C             the intermediary direct input/output transmission matrix
		C             D1 (see METHOD).
		C             The array D is not referenced if JOBD = 'Z'.
		C
		C     LDD     INTEGER
		C             The leading dimension of array D.
		C             LDD >= MAX(1,P) if JOBD = 'D'.
		C             LDD >= 1 if JOBD = 'Z'.
		C
		C     F       (input) DOUBLE PRECISION array, dimension (LDF,P)
		C             If FBTYPE = 'O', the leading M-by-P part of this array
		C             must contain the output feedback matrix F.
		C             If FBTYPE = 'I', then the feedback matrix is assumed to be
		C             an M x M order identity matrix.
		C             The array F is not referenced if FBTYPE = 'I'  or
		C             ALPHA = 0.
		C
		C     LDF     INTEGER
		C             The leading dimension of array F.
		C             LDF >= MAX(1,M) if FBTYPE = 'O' and ALPHA <> 0.
		C             LDF >= 1 if FBTYPE = 'I' or ALPHA = 0.
		C
		C     K       (input) DOUBLE PRECISION array, dimension (LDK,N)
		C             The leading M-by-N part of this array must contain the
		C             state feedback matrix K.
		C             The array K is not referenced if BETA = 0.
		C
		C     LDK     INTEGER
		C             The leading dimension of the array K.
		C             LDK >= MAX(1,M) if BETA <> 0.
		C             LDK >= 1 if BETA = 0.
		C
		C     G       (input) DOUBLE PRECISION array, dimension (LDG,MV)
		C             The leading M-by-MV part of this array must contain the
		C             system input scaling matrix G.
		C
		C     LDG     INTEGER
		C             The leading dimension of the array G.  LDG >= MAX(1,M).
		C
		C     H       (input) DOUBLE PRECISION array, dimension (LDH,P)
		C             The leading PZ-by-P part of this array must contain the
		C             system output scaling matrix H.
		C
		C     LDH     INTEGER
		C             The leading dimension of the array H.  LDH >= MAX(1,PZ).
		C
		C     RCOND   (output) DOUBLE PRECISION
		C             The reciprocal condition number of the matrix
		C             I - alpha*D*F.
		C
		C     BC      (output) DOUBLE PRECISION array, dimension (LDBC,MV)
		C             The leading N-by-MV part of this array contains the input
		C             matrix Bc of the closed-loop system.
		C
		C     LDBC    INTEGER
		C             The leading dimension of array BC.  LDBC >= MAX(1,N).
		C
		C     CC      (output) DOUBLE PRECISION array, dimension (LDCC,N)
		C             The leading PZ-by-N part of this array contains the
		C             system output matrix Cc of the closed-loop system.
		C
		C     LDCC    INTEGER
		C             The leading dimension of array CC.
		C             LDCC >= MAX(1,PZ) if N > 0.
		C             LDCC >= 1 if N = 0.
		C
		C     DC      (output) DOUBLE PRECISION array, dimension (LDDC,MV)
		C             If JOBD = 'D', the leading PZ-by-MV part of this array
		C             contains the direct input/output transmission matrix Dc
		C             of the closed-loop system.
		C             The array DC is not referenced if JOBD = 'Z'.
		C
		C     LDDC    INTEGER
		C             The leading dimension of array DC.
		C             LDDC >= MAX(1,PZ) if JOBD = 'D'.
		C             LDDC >= 1 if JOBD = 'Z'.
		C
		C     Workspace
		C
		C     IWORK   INTEGER array, dimension (LIWORK)
		C             LIWORK >= MAX(1,2*P) if JOBD = 'D'.
		C             LIWORK >= 1 if JOBD = 'Z'.
		C             IWORK is not referenced if JOBD = 'Z'.
		C
		C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
		C
		C     LDWORK  INTEGER
		C             The length of the array DWORK.
		C             LDWORK >= wspace, where
		C                   wspace = MAX( 1, M, P*MV, P*P + 4*P ) if JOBD = 'D',
		C                   wspace = MAX( 1, M ) if JOBD = 'Z'.
		C             For best performance, LDWORK >= MAX( wspace, N*M, N*P ).
		C
		C     Error Indicator
		C
		C     INFO    INTEGER
		C             = 0:  successful exit;
		C             < 0:  if INFO = -i, the i-th argument had an illegal
		C                   value;
		C             = 1:  if the matrix I - alpha*D*F is numerically singular.
		C
		C     METHOD
		C
		C     The matrices of the closed-loop system have the expressions:
		C
		C     Ac = A1 + beta*B1*K,      Bc = B1*G,
		C     Cc = H*(C1 + beta*D1*K),  Dc = H*D1*G,
		C
		C     where
		C
		C     A1 = A + alpha*B*F*E*C,   B1 = B + alpha*B*F*E*D,
		C     C1 = E*C,                 D1 = E*D,
		C
		C     with E = (I - alpha*D*F)**-1.
		C
		C     NUMERICAL ASPECTS
		C
		C     The accuracy of computations basically depends on the conditioning
		C     of the matrix I - alpha*D*F. If RCOND is very small, it is likely
		C     that the computed results are inaccurate.
		C
		C     CONTRIBUTORS
		C
		C     A. Varga, German Aerospace Research Establishment,
		C     Oberpfaffenhofen, Germany, and V. Sima, Katholieke Univ. Leuven,
		C     Belgium, Nov. 1996.
		C
		C     REVISIONS
		C
		C     January 14, 1997, February 18, 1998.
		C     V. Sima, Research Institute for Informatics, Bucharest, July 2003,
		C     Jan. 2005.
		C
		C     KEYWORDS
		C
		C     Multivariable system, state-space model, state-space
		C     representation.
		*/

		void AB05RD(unsigned char *, unsigned char *, int *, int *, int *, int *, int *, double *,
		            double *, double *, int *, double *, int *, double *, int *, double *, int *,
					double *, int *, double *, int *, double *, int *, double *, int *, double *,
					double *, int *, double *, int *, double *, int *, int *, double *, int *, int *);

		/*
		SLICOT RELEASE 5.5.
		C
		C     Copyright (c) 2002-2012 NICONET e.V.
		C
		C     PURPOSE
		C
		C     To construct for a given state space system (A,B,C,D) the closed-
		C     loop system (Ac,Bc,Cc,Dc) corresponding to the output feedback
		C     control law
		C
		C          u = alpha*F*y + v.
		C
		C     ARGUMENTS
		C
		C     Mode Parameters
		C
		C     FBTYPE  CHARACTER*1
		C             Specifies the type of the feedback law as follows:
		C             = 'I':  Unitary output feedback (F = I);
		C             = 'O':  General output feedback.
		C
		C     JOBD    CHARACTER*1
		C             Specifies whether or not a non-zero matrix D appears in
		C             the given state space model:
		C             = 'D':  D is present;
		C             = 'Z':  D is assumed a zero matrix.
		C
		C     Input/Output Parameters
		C
		C     N       (input) INTEGER
		C             The number of state variables, i.e. the order of the
		C             matrix A, the number of rows of B and the number of
		C             columns of C.  N >= 0.
		C
		C     M       (input) INTEGER
		C             The number of input variables, i.e. the number of columns
		C             of matrices B and D, and the number of rows of F.  M >= 0.
		C
		C     P       (input) INTEGER
		C             The number of output variables, i.e. the number of rows of
		C             matrices C and D, and the number of columns of F.  P >= 0
		C             and P = M if FBTYPE = 'I'.
		C
		C     ALPHA   (input) DOUBLE PRECISION
		C             The coefficient alpha in the output feedback law.
		C
		C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
		C             On entry, the leading N-by-N part of this array must
		C             contain the system state transition matrix A.
		C             On exit, the leading N-by-N part of this array contains
		C             the state matrix Ac of the closed-loop system.
		C
		C     LDA     INTEGER
		C             The leading dimension of array A.  LDA >= MAX(1,N).
		C
		C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
		C             On entry, the leading N-by-M part of this array must
		C             contain the system input matrix B.
		C             On exit, the leading N-by-M part of this array contains
		C             the input matrix Bc of the closed-loop system.
		C
		C     LDB     INTEGER
		C             The leading dimension of array B.  LDB >= MAX(1,N).
		C
		C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
		C             On entry, the leading P-by-N part of this array must
		C             contain the system output matrix C.
		C             On exit, the leading P-by-N part of this array contains
		C             the output matrix Cc of the closed-loop system.
		C
		C     LDC     INTEGER
		C             The leading dimension of array C.
		C             LDC >= MAX(1,P) if N > 0.
		C             LDC >= 1 if N = 0.
		C
		C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
		C             On entry, the leading P-by-M part of this array must
		C             contain the system direct input/output transmission
		C             matrix D.
		C             On exit, if JOBD = 'D', the leading P-by-M part of this
		C             array contains the direct input/output transmission
		C             matrix Dc of the closed-loop system.
		C             The array D is not referenced if JOBD = 'Z'.
		C
		C     LDD     INTEGER
		C             The leading dimension of array D.
		C             LDD >= MAX(1,P) if JOBD = 'D'.
		C             LDD >= 1 if JOBD = 'Z'.
		C
		C     F       (input) DOUBLE PRECISION array, dimension (LDF,P)
		C             If FBTYPE = 'O', the leading M-by-P part of this array
		C             must contain the output feedback matrix F.
		C             If FBTYPE = 'I', then the feedback matrix is assumed to be
		C             an M x M order identity matrix.
		C             The array F is not referenced if FBTYPE = 'I' or
		C             ALPHA = 0.
		C
		C     LDF     INTEGER
		C             The leading dimension of array F.
		C             LDF >= MAX(1,M) if FBTYPE = 'O' and ALPHA <> 0.
		C             LDF >= 1 if FBTYPE = 'I' or ALPHA = 0.
		C
		C     RCOND   (output) DOUBLE PRECISION
		C             The reciprocal condition number of the matrix
		C             I - alpha*D*F.
		C
		C     Workspace
		C
		C     IWORK   INTEGER array, dimension (LIWORK)
		C             LIWORK >= MAX(1,2*P) if JOBD = 'D'.
		C             LIWORK >= 1 if JOBD = 'Z'.
		C             IWORK is not referenced if JOBD = 'Z'.
		C
		C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
		C
		C     LDWORK  INTEGER
		C             The length of the array DWORK.
		C             LDWORK >= wspace, where
		C                       wspace = MAX( 1, M, P*P + 4*P ) if JOBD = 'D',
		C                       wspace = MAX( 1, M ) if JOBD = 'Z'.
		C             For best performance, LDWORK >= MAX( wspace, N*M, N*P ).
		C
		C     Error Indicator
		C
		C     INFO    INTEGER
		C             = 0:  successful exit;
		C             < 0:  if INFO = -i, the i-th argument had an illegal
		C                   value;
		C             = 1:  if the matrix I - alpha*D*F is numerically singular.
		C
		C     METHOD
		C
		C     The matrices of the closed-loop system have the expressions:
		C
		C     Ac = A + alpha*B*F*E*C,  Bc = B + alpha*B*F*E*D,
		C     Cc = E*C,                Dc = E*D,
		C
		C     where E = (I - alpha*D*F)**-1.
		C
		C     NUMERICAL ASPECTS
		C
		C     The accuracy of computations basically depends on the conditioning
		C     of the matrix I - alpha*D*F.  If RCOND is very small, it is likely
		C     that the computed results are inaccurate.
		C
		C     CONTRIBUTORS
		C
		C     A. Varga, German Aerospace Research Establishment,
		C     Oberpfaffenhofen, Germany, and V. Sima, Katholieke Univ. Leuven,
		C     Belgium, Nov. 1996.
		C
		C     REVISIONS
		C
		C     January 14, 1997.
		C     V. Sima, Research Institute for Informatics, Bucharest, July 2003.
		C
		C     KEYWORDS
		C
		C     Multivariable system, state-space model, state-space
		C     representation.
		*/

		void  AB05SD(unsigned char *, unsigned char *, int *, int *, int *, double *, double *,
		             int *, double *, int *, double *, int *, double *, int *, double *, int *,
					 double *, int *, double *, int *, int *);

		/*
		SLICOT RELEASE 5.5.
		C
		C     Copyright (c) 2002-2012 NICONET e.V.
		C
		C     PURPOSE
		C
		C     To find the dual of a given state-space representation.
		C
		C     ARGUMENTS
		C
		C     Mode Parameters
		C
		C     JOBD    CHARACTER*1
		C             Specifies whether or not a non-zero matrix D appears in
		C             the given state space model:
		C             = 'D':  D is present;
		C             = 'Z':  D is assumed a zero matrix.
		C
		C     Input/Output Parameters
		C
		C     N       (input) INTEGER
		C             The order of the state-space representation.  N >= 0.
		C
		C     M       (input) INTEGER
		C             The number of system inputs.  M >= 0.
		C
		C     P       (input) INTEGER
		C             The number of system outputs.  P >= 0.
		C
		C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
		C             On entry, the leading N-by-N part of this array must
		C             contain the original state dynamics matrix A.
		C             On exit, the leading N-by-N part of this array contains
		C             the dual state dynamics matrix A'.
		C
		C     LDA     INTEGER
		C             The leading dimension of array A.  LDA >= MAX(1,N).
		C
		C     B       (input/output) DOUBLE PRECISION array, dimension
		C             (LDB,MAX(M,P))
		C             On entry, the leading N-by-M part of this array must
		C             contain the original input/state matrix B.
		C             On exit, the leading N-by-P part of this array contains
		C             the dual input/state matrix C'.
		C
		C     LDB     INTEGER
		C             The leading dimension of array B.  LDB >= MAX(1,N).
		C
		C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
		C             On entry, the leading P-by-N part of this array must
		C             contain the original state/output matrix C.
		C             On exit, the leading M-by-N part of this array contains
		C             the dual state/output matrix B'.
		C
		C     LDC     INTEGER
		C             The leading dimension of array C.
		C             LDC >= MAX(1,M,P) if N > 0.
		C             LDC >= 1 if N = 0.
		C
		C     D       (input/output) DOUBLE PRECISION array, dimension
		C             (LDD,MAX(M,P))
		C             On entry, if JOBD = 'D', the leading P-by-M part of this
		C             array must contain the original direct transmission
		C             matrix D.
		C             On exit, if JOBD = 'D', the leading M-by-P part of this
		C             array contains the dual direct transmission matrix D'.
		C             The array D is not referenced if JOBD = 'Z'.
		C
		C     LDD     INTEGER
		C             The leading dimension of array D.
		C             LDD >= MAX(1,M,P) if JOBD = 'D'.
		C             LDD >= 1 if JOBD = 'Z'.
		C
		C     Error Indicator
		C
		C     INFO    INTEGER
		C             = 0:  successful exit;
		C             < 0:  if INFO = -i, the i-th argument had an illegal
		C                   value.
		C
		C     METHOD
		C
		C     If the given state-space representation is the M-input/P-output
		C     (A,B,C,D), its dual is simply the P-input/M-output (A',C',B',D').
		C
		C     REFERENCES
		C
		C     None
		C
		C     NUMERICAL ASPECTS
		C
		C     None
		C
		C     CONTRIBUTOR
		C
		C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996.
		C     Supersedes Release 2.0 routine AB07AD by T.W.C.Williams, Kingston
		C     Polytechnic, United Kingdom, March 1982.
		C
		C     REVISIONS
		C
		C     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2004.
		C
		C     KEYWORDS
		C
		C     Dual system, state-space model, state-space representation.
		C
		*/

		void AB07MD(unsigned char *, int *, int *, int *, double *, int *, double *, int *,
		            double *, int *, double *, int *, int *);

		/*
		SLICOT RELEASE 5.5.
		C
		C     Copyright (c) 2002-2012 NICONET e.V.
		C
		C     PURPOSE
		C
		C     To compute the inverse (Ai,Bi,Ci,Di) of a given system (A,B,C,D).
		C
		C     ARGUMENTS
		C
		C     Input/Output Parameters
		C
		C     N       (input) INTEGER
		C             The order of the state matrix A.  N >= 0.
		C
		C     M       (input) INTEGER
		C             The number of system inputs and outputs.  M >= 0.
		C
		C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
		C             On entry, the leading N-by-N part of this array must
		C             contain the state matrix A of the original system.
		C             On exit, the leading N-by-N part of this array contains
		C             the state matrix Ai of the inverse system.
		C
		C     LDA     INTEGER
		C             The leading dimension of the array A.  LDA >= MAX(1,N).
		C
		C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
		C             On entry, the leading N-by-M part of this array must
		C             contain the input matrix B of the original system.
		C             On exit, the leading N-by-M part of this array contains
		C             the input matrix Bi of the inverse system.
		C
		C     LDB     INTEGER
		C             The leading dimension of the array B.  LDB >= MAX(1,N).
		C
		C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
		C             On entry, the leading M-by-N part of this array must
		C             contain the output matrix C of the original system.
		C             On exit, the leading M-by-N part of this array contains
		C             the output matrix Ci of the inverse system.
		C
		C     LDC     INTEGER
		C             The leading dimension of the array C.  LDC >= MAX(1,M).
		C
		C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
		C             On entry, the leading M-by-M part of this array must
		C             contain the feedthrough matrix D of the original system.
		C             On exit, the leading M-by-M part of this array contains
		C             the feedthrough matrix Di of the inverse system.
		C
		C     LDD     INTEGER
		C             The leading dimension of the array D.  LDD >= MAX(1,M).
		C
		C     RCOND   (output) DOUBLE PRECISION
		C             The estimated reciprocal condition number of the
		C             feedthrough matrix D of the original system.
		C
		C     Workspace
		C
		C     IWORK   INTEGER array, dimension (2*M)
		C
		C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
		C             On exit, if INFO = 0 or M+1, DWORK(1) returns the optimal
		C             value of LDWORK.
		C
		C     LDWORK  INTEGER
		C             The length of the array DWORK.  LDWORK >= MAX(1,4*M).
		C             For good performance, LDWORK should be larger.
		C
		C             If LDWORK = -1, then a workspace query is assumed;
		C             the routine only calculates the optimal size of the
		C             DWORK array, returns this value as the first entry of
		C             the DWORK array, and no error message related to LDWORK
		C             is issued by XERBLA.
		C
		C     Error Indicator
		C
		C     INFO    INTEGER
		C             = 0:  successful exit;
		C             < 0:  if INFO = -i, the i-th argument had an illegal
		C                   value;
		C             = i:  the matrix D is exactly singular; the (i,i) diagonal
		C                   element is zero, i <= M; RCOND was set to zero;
		C             = M+1:  the matrix D is numerically singular, i.e., RCOND
		C                   is less than the relative machine precision, EPS
		C                   (see LAPACK Library routine DLAMCH). The
		C                   calculations have been completed, but the results
		C                   could be very inaccurate.
		C
		C     METHOD
		C
		C     The matrices of the inverse system are computed with the formulas:
		C                   -1              -1         -1           -1
		C       Ai = A - B*D  *C,  Bi = -B*D  ,  Ci = D  *C,  Di = D  .
		C
		C     NUMERICAL ASPECTS
		C
		C     The accuracy depends mainly on the condition number of the matrix
		C     D to be inverted. The estimated reciprocal condition number is
		C     returned in RCOND.
		C
		C     CONTRIBUTORS
		C
		C     A. Varga, German Aerospace Center, Oberpfaffenhofen, March 2000.
		C     D. Sima, University of Bucharest, April 2000.
		C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000.
		C     Based on the routine SYSINV, A. Varga, 1992.
		C
		C     REVISIONS
		C
		C     A. Varga, German Aerospace Center, Oberpfaffenhofen, July 2000.
		C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2011.
		C
		C     KEYWORDS
		C
		C     Inverse system, state-space model, state-space representation.
		*/

		void AB07ND(int *, int *, double *, int *, double *, int *, double *, int *, 
		            double *, int *, double *, int *, double *, int *, int *);

		/*
		SLICOT RELEASE 5.5.
		C
		C     Copyright (c) 2002-2012 NICONET e.V.
		C
		C     PURPOSE
		C
		C     To compute the normal rank of the transfer-function matrix of a
		C     state-space model (A,B,C,D).
		C
		C     ARGUMENTS
		C
		C     Mode Parameters
		C
		C     EQUIL   CHARACTER*1
		C             Specifies whether the user wishes to balance the compound
		C             matrix (see METHOD) as follows:
		C             = 'S':  Perform balancing (scaling);
		C             = 'N':  Do not perform balancing.
		C
		C     Input/Output Parameters
		C
		C     N       (input) INTEGER
		C             The number of state variables, i.e., the order of the
		C             matrix A.  N >= 0.
		C
		C     M       (input) INTEGER
		C             The number of system inputs.  M >= 0.
		C
		C     P       (input) INTEGER
		C             The number of system outputs.  P >= 0.
		C
		C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
		C             The leading N-by-N part of this array must contain the
		C             state dynamics matrix A of the system.
		C
		C     LDA     INTEGER
		C             The leading dimension of array A.  LDA >= MAX(1,N).
		C
		C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
		C             The leading N-by-M part of this array must contain the
		C             input/state matrix B of the system.
		C
		C     LDB     INTEGER
		C             The leading dimension of array B.  LDB >= MAX(1,N).
		C
		C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
		C             The leading P-by-N part of this array must contain the
		C             state/output matrix C of the system.
		C
		C     LDC     INTEGER
		C             The leading dimension of array C.  LDC >= MAX(1,P).
		C
		C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
		C             The leading P-by-M part of this array must contain the
		C             direct transmission matrix D of the system.
		C
		C     LDD     INTEGER
		C             The leading dimension of array D.  LDD >= MAX(1,P).
		C
		C     RANK    (output) INTEGER
		C             The normal rank of the transfer-function matrix.
		C
		C     Tolerances
		C
		C     TOL     DOUBLE PRECISION
		C             A tolerance used in rank decisions to determine the
		C             effective rank, which is defined as the order of the
		C             largest leading (or trailing) triangular submatrix in the
		C             QR (or RQ) factorization with column (or row) pivoting
		C             whose estimated condition number is less than 1/TOL.
		C             If the user sets TOL to be less than SQRT((N+P)*(N+M))*EPS
		C             then the tolerance is taken as SQRT((N+P)*(N+M))*EPS,
		C             where EPS is the machine precision (see LAPACK Library
		C             Routine DLAMCH).
		C
		C     Workspace
		C
		C     IWORK   INTEGER array, dimension (2*N+MAX(M,P)+1)
		C
		C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
		C             On exit, if INFO = 0, DWORK(1) returns the optimal value
		C             of LDWORK.
		C
		C     LDWORK  INTEGER
		C             The length of the array DWORK.
		C             LDWORK >= (N+P)*(N+M) +
		C                       MAX( MIN(P,M) + MAX(3*M-1,N), 1,
		C                            MIN(P,N) + MAX(3*P-1,N+P,N+M) )
		C             For optimum performance LDWORK should be larger.
		C
		C             If LDWORK = -1, then a workspace query is assumed;
		C             the routine only calculates the optimal size of the
		C             DWORK array, returns this value as the first entry of
		C             the DWORK array, and no error message related to LDWORK
		C             is issued by XERBLA.
		C
		C     Error Indicator
		C
		C     INFO    INTEGER
		C             = 0:  successful exit;
		C             < 0:  if INFO = -i, the i-th argument had an illegal
		C                   value.
		C
		C     METHOD
		C
		C     The routine reduces the (N+P)-by-(M+N) compound matrix (B  A)
		C                                                            (D  C)
		C
		C     to one with the same invariant zeros and with D of full row rank.
		C     The normal rank of the transfer-function matrix is the rank of D.
		C
		C     REFERENCES
		C
		C     [1] Svaricek, F.
		C         Computation of the Structural Invariants of Linear
		C         Multivariable Systems with an Extended Version of
		C         the Program ZEROS.
		C         System & Control Letters, 6, pp. 261-266, 1985.
		C
		C     [2] Emami-Naeini, A. and Van Dooren, P.
		C         Computation of Zeros of Linear Multivariable Systems.
		C         Automatica, 18, pp. 415-430, 1982.
		C
		C     NUMERICAL ASPECTS
		C
		C     The algorithm is backward stable (see [2] and [1]).
		C
		C     CONTRIBUTOR
		C
		C     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2001.
		C
		C     REVISIONS
		C
		C     V. Sima, Research Institute for Informatics, Bucharest, June 2001,
		C     Dec. 2003, Jan. 2009, Mar. 2009, Apr. 2009.
		C
		C     KEYWORDS
		C
		C     Multivariable system, orthogonal transformation,
		C     structural invariant.
		*/

		void AB08MD(unsigned char *, int *, int *, int *, double *, int *, double *, int *,
		            double *, int *, double *, int *, int *, double *, int *, double *, int *, int *);

		/*
		SLICOT RELEASE 5.5.
		C
		C     Copyright (c) 2002-2012 NICONET e.V.
		C
		C     PURPOSE
		C
		C     To compute the normal rank of the transfer-function matrix of a
		C     state-space model (A,B,C,D).
		C
		C     ARGUMENTS
		C
		C     Mode Parameters
		C
		C     EQUIL   CHARACTER*1
		C             Specifies whether the user wishes to balance the compound
		C             matrix (see METHOD) as follows:
		C             = 'S':  Perform balancing (scaling);
		C             = 'N':  Do not perform balancing.
		C
		C     Input/Output Parameters
		C
		C     N       (input) INTEGER
		C             The number of state variables, i.e., the order of the
		C             matrix A.  N >= 0.
		C
		C     M       (input) INTEGER
		C             The number of system inputs.  M >= 0.
		C
		C     P       (input) INTEGER
		C             The number of system outputs.  P >= 0.
		C
		C     A       (input) COMPLEX*16 array, dimension (LDA,N)
		C             The leading N-by-N part of this array must contain the
		C             state dynamics matrix A of the system.
		C
		C     LDA     INTEGER
		C             The leading dimension of array A.  LDA >= MAX(1,N).
		C
		C     B       (input) COMPLEX*16 array, dimension (LDB,M)
		C             The leading N-by-M part of this array must contain the
		C             input/state matrix B of the system.
		C
		C     LDB     INTEGER
		C             The leading dimension of array B.  LDB >= MAX(1,N).
		C
		C     C       (input) COMPLEX*16 array, dimension (LDC,N)
		C             The leading P-by-N part of this array must contain the
		C             state/output matrix C of the system.
		C
		C     LDC     INTEGER
		C             The leading dimension of array C.  LDC >= MAX(1,P).
		C
		C     D       (input) COMPLEX*16 array, dimension (LDD,M)
		C             The leading P-by-M part of this array must contain the
		C             direct transmission matrix D of the system.
		C
		C     LDD     INTEGER
		C             The leading dimension of array D.  LDD >= MAX(1,P).
		C
		C     RANK    (output) INTEGER
		C             The normal rank of the transfer-function matrix.
		C
		C     Tolerances
		C
		C     TOL     DOUBLE PRECISION
		C             A tolerance used in rank decisions to determine the
		C             effective rank, which is defined as the order of the
		C             largest leading (or trailing) triangular submatrix in the
		C             QR (or RQ) factorization with column (or row) pivoting
		C             whose estimated condition number is less than 1/TOL.
		C             If the user sets TOL to be less than SQRT((N+P)*(N+M))*EPS
		C             then the tolerance is taken as SQRT((N+P)*(N+M))*EPS,
		C             where EPS is the machine precision (see LAPACK Library
		C             Routine DLAMCH).
		C
		C     Workspace
		C
		C     IWORK   INTEGER array, dimension (2*N+MAX(M,P)+1)
		C
		C     DWORK   DOUBLE PRECISION array, dimension (2*MAX(M,P))
		C
		C     ZWORK   COMPLEX*16 array, dimension (LZWORK)
		C             On exit, if INFO = 0, ZWORK(1) returns the optimal value
		C             of LZWORK.
		C
		C     LZWORK  INTEGER
		C             The length of the array ZWORK.
		C             LZWORK >= (N+P)*(N+M) + MAX(MIN(P,M) + MAX(3*M-1,N), 1,
		C                                         MIN(P,N) + MAX(3*P-1,N+P,N+M))
		C             For optimum performance LZWORK should be larger.
		C
		C             If LZWORK = -1, then a workspace query is assumed;
		C             the routine only calculates the optimal size of the
		C             ZWORK array, returns this value as the first entry of
		C             the ZWORK array, and no error message related to LZWORK
		C             is issued by XERBLA.
		C
		C     Error Indicator
		C
		C     INFO    INTEGER
		C             = 0:  successful exit;
		C             < 0:  if INFO = -i, the i-th argument had an illegal
		C                   value.
		C
		C     METHOD
		C
		C     The routine reduces the (N+P)-by-(M+N) compound matrix (B  A)
		C                                                            (D  C)
		C
		C     to one with the same invariant zeros and with D of full row rank.
		C     The normal rank of the transfer-function matrix is the rank of D.
		C
		C     REFERENCES
		C
		C     [1] Svaricek, F.
		C         Computation of the Structural Invariants of Linear
		C         Multivariable Systems with an Extended Version of
		C         the Program ZEROS.
		C         System & Control Letters, 6, pp. 261-266, 1985.
		C
		C     [2] Emami-Naeini, A. and Van Dooren, P.
		C         Computation of Zeros of Linear Multivariable Systems.
		C         Automatica, 18, pp. 415-430, 1982.
		C
		C     NUMERICAL ASPECTS
		C
		C     The algorithm is backward stable (see [2] and [1]).
		C
		C     CONTRIBUTOR
		C
		C     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2001.
		C     Complex version: V. Sima, Research Institute for Informatics,
		C     Bucharest, Dec. 2008.
		C
		C     REVISIONS
		C
		C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2009,
		C     Apr. 2009.
		C
		C     KEYWORDS
		C
		C     Multivariable system, unitary transformation,
		C     structural invariant.
		*/

		void AB08MZ(unsigned char *, int *, int *, int *, std::complex<double> *, int *, std::complex<double> *,
		            int *, std::complex<double> *, int *, std::complex<double> *, int *, int *, double *, int *,
					double *, std::complex<double> *, int *, int *);

		/*
		SLICOT RELEASE 5.5.
		C
		C     Copyright (c) 2002-2012 NICONET e.V.
		C
		C     PURPOSE
		C
		C     To construct for a linear multivariable system described by a
		C     state-space model (A,B,C,D) a regular pencil (A - lambda*B ) which
		C                                                    f          f
		C     has the invariant zeros of the system as generalized eigenvalues.
		C     The routine also computes the orders of the infinite zeros and the
		C     right and left Kronecker indices of the system (A,B,C,D).
		C
		C     ARGUMENTS
		C
		C     Mode Parameters
		C
		C     EQUIL   CHARACTER*1
		C             Specifies whether the user wishes to balance the compound
		C             matrix (see METHOD) as follows:
		C             = 'S':  Perform balancing (scaling);
		C             = 'N':  Do not perform balancing.
		C
		C     Input/Output Parameters
		C
		C     N       (input) INTEGER
		C             The number of state variables, i.e., the order of the
		C             matrix A.  N >= 0.
		C
		C     M       (input) INTEGER
		C             The number of system inputs.  M >= 0.
		C
		C     P       (input) INTEGER
		C             The number of system outputs.  P >= 0.
		C
		C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
		C             The leading N-by-N part of this array must contain the
		C             state dynamics matrix A of the system.
		C
		C     LDA     INTEGER
		C             The leading dimension of array A.  LDA >= MAX(1,N).
		C
		C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
		C             The leading N-by-M part of this array must contain the
		C             input/state matrix B of the system.
		C
		C     LDB     INTEGER
		C             The leading dimension of array B.  LDB >= MAX(1,N).
		C
		C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
		C             The leading P-by-N part of this array must contain the
		C             state/output matrix C of the system.
		C
		C     LDC     INTEGER
		C             The leading dimension of array C.  LDC >= MAX(1,P).
		C
		C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
		C             The leading P-by-M part of this array must contain the
		C             direct transmission matrix D of the system.
		C
		C     LDD     INTEGER
		C             The leading dimension of array D.  LDD >= MAX(1,P).
		C
		C     NU      (output) INTEGER
		C             The number of (finite) invariant zeros.
		C
		C     RANK    (output) INTEGER
		C             The normal rank of the transfer function matrix.
		C
		C     DINFZ   (output) INTEGER
		C             The maximum degree of infinite elementary divisors.
		C
		C     NKROR   (output) INTEGER
		C             The number of right Kronecker indices.
		C
		C     NKROL   (output) INTEGER
		C             The number of left Kronecker indices.
		C
		C     INFZ    (output) INTEGER array, dimension (N)
		C             The leading DINFZ elements of INFZ contain information
		C             on the infinite elementary divisors as follows:
		C             the system has INFZ(i) infinite elementary divisors
		C             of degree i, where i = 1,2,...,DINFZ.
		C
		C     KRONR   (output) INTEGER array, dimension (MAX(N,M)+1)
		C             The leading NKROR elements of this array contain the
		C             right Kronecker (column) indices.
		C
		C     KRONL   (output) INTEGER array, dimension (MAX(N,P)+1)
		C             The leading NKROL elements of this array contain the
		C             left Kronecker (row) indices.
		C
		C     AF      (output) DOUBLE PRECISION array, dimension
		C             (LDAF,N+MIN(P,M))
		C             The leading NU-by-NU part of this array contains the
		C             coefficient matrix A  of the reduced pencil. The remainder
		C                                 f
		C             of the leading (N+M)-by-(N+MIN(P,M)) part is used as
		C             internal workspace.
		C
		C     LDAF    INTEGER
		C             The leading dimension of array AF.  LDAF >= MAX(1,N+M).
		C
		C     BF      (output) DOUBLE PRECISION array, dimension (LDBF,N+M)
		C             The leading NU-by-NU part of this array contains the
		C             coefficient matrix B  of the reduced pencil. The
		C                                 f
		C             remainder of the leading (N+P)-by-(N+M) part is used as
		C             internal workspace.
		C
		C     LDBF    INTEGER
		C             The leading dimension of array BF.  LDBF >= MAX(1,N+P).
		C
		C     Tolerances
		C
		C     TOL     DOUBLE PRECISION
		C             A tolerance used in rank decisions to determine the
		C             effective rank, which is defined as the order of the
		C             largest leading (or trailing) triangular submatrix in the
		C             QR (or RQ) factorization with column (or row) pivoting
		C             whose estimated condition number is less than 1/TOL.
		C             If the user sets TOL to be less than SQRT((N+P)*(N+M))*EPS
		C             then the tolerance is taken as SQRT((N+P)*(N+M))*EPS,
		C             where EPS is the machine precision (see LAPACK Library
		C             Routine DLAMCH).
		C
		C     Workspace
		C
		C     IWORK   INTEGER array, dimension (MAX(M,P))
		C
		C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
		C             On exit, if INFO = 0, DWORK(1) returns the optimal value
		C             of LDWORK.
		C
		C     LDWORK  INTEGER
		C             The length of the array DWORK.
		C             LDWORK >= MAX( 1, MIN(P,M) + MAX(3*M-1,N),
		C                               MIN(P,N) + MAX(3*P-1,N+P,N+M),
		C                               MIN(M,N) + MAX(3*M-1,N+M) ).
		C             An upper bound is MAX(s,N) + MAX(3*s-1,N+s), with
		C             s = MAX(M,P).
		C             For optimum performance LDWORK should be larger.
		C
		C             If LDWORK = -1, then a workspace query is assumed;
		C             the routine only calculates the optimal size of the
		C             DWORK array, returns this value as the first entry of
		C             the DWORK array, and no error message related to LDWORK
		C             is issued by XERBLA.
		C
		C     Error Indicator
		C
		C     INFO    INTEGER
		C             = 0:  successful exit;
		C             < 0:  if INFO = -i, the i-th argument had an illegal
		C                   value.
		C
		C     METHOD
		C
		C     The routine extracts from the system matrix of a state-space
		C     system (A,B,C,D) a regular pencil A - lambda*B  which has the
		C                                        f          f
		C     invariant zeros of the system as generalized eigenvalues as
		C     follows:
		C
		C        (a) construct the (N+P)-by-(N+M) compound matrix (B  A);
		C                                                         (D  C)
		C
		C        (b) reduce the above system to one with the same invariant
		C            zeros and with D of full row rank;
		C
		C        (c) pertranspose the system;
		C
		C        (d) reduce the system to one with the same invariant zeros and
		C            with D square invertible;
		C
		C        (e) perform a unitary transformation on the columns of
		C            (A - lambda*I  B) in order to reduce it to
		C            (      C       D)
		C
		C            (A  - lambda*B   X)
		C            ( f           f   ), with Y and B  square invertible;
		C            (      0         Y)              f
		C
		C        (f) compute the right and left Kronecker indices of the system
		C            (A,B,C,D), which together with the orders of the infinite
		C            zeros (determined by steps (a) - (e)) constitute the
		C            complete set of structural invariants under strict
		C            equivalence transformations of a linear system.
		C
		C     REFERENCES
		C
		C     [1] Svaricek, F.
		C         Computation of the Structural Invariants of Linear
		C         Multivariable Systems with an Extended Version of
		C         the Program ZEROS.
		C         System & Control Letters, 6, pp. 261-266, 1985.
		C
		C     [2] Emami-Naeini, A. and Van Dooren, P.
		C         Computation of Zeros of Linear Multivariable Systems.
		C         Automatica, 18, pp. 415-430, 1982.
		C
		C     NUMERICAL ASPECTS
		C
		C     The algorithm is backward stable (see [2] and [1]).
		C
		C     FURTHER COMMENTS
		C
		C     In order to compute the invariant zeros of the system explicitly,
		C     a call to this routine may be followed by a call to the LAPACK
		C     Library routine DGGEV with A = A , B = B  and N = NU.
		C                                     f       f
		C     If RANK = 0, the routine DGEEV can be used (since B = I).
		C                                                        f
		C     CONTRIBUTOR
		C
		C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996.
		C     Supersedes Release 2.0 routine AB08BD by F. Svaricek.
		C
		C     REVISIONS
		C
		C     Oct. 1997, Feb. 1998, Dec. 2003, March 2004, Jan. 2009, Mar. 2009,
		C     Apr. 2009, Apr. 2011.
		C
		C     KEYWORDS
		C
		C     Generalized eigenvalue problem, Kronecker indices, multivariable
		C     system, orthogonal transformation, structural invariant.
		*/

		void AB08ND(unsigned char *, int *, int *, int *, double *, int *, double *, int *,
		           double *, int *, double *, int *, int *, int *, int *, int *, int *, 
				   int *, int *, int *, double *, int *, double *, int *, double *, int *,
				   double *, int *, int *);

		/*
		SLICOT RELEASE 5.5.
		C
		C     Copyright (c) 2002-2012 NICONET e.V.
		C
		C     PURPOSE
		C
		C     To extract from the (N+P)-by-(M+N) system
		C                  ( B  A )
		C                  ( D  C )
		C     an (NU+MU)-by-(M+NU) "reduced" system
		C                  ( B' A')
		C                  ( D' C')
		C     having the same transmission zeros but with D' of full row rank.
		C
		C     ARGUMENTS
		C
		C     Input/Output Parameters
		C
		C     N       (input) INTEGER
		C             The number of state variables.  N >= 0.
		C
		C     M       (input) INTEGER
		C             The number of system inputs.  M >= 0.
		C
		C     P       (input) INTEGER
		C             The number of system outputs.  P >= 0.
		C
		C     RO      (input/output) INTEGER
		C             On entry,
		C             = P     for the original system;
		C             = MAX(P-M, 0) for the pertransposed system.
		C             On exit, RO contains the last computed rank.
		C
		C     SIGMA   (input/output) INTEGER
		C             On entry,
		C             = 0  for the original system;
		C             = M  for the pertransposed system.
		C             On exit, SIGMA contains the last computed value sigma in
		C             the algorithm.
		C
		C     SVLMAX  (input) DOUBLE PRECISION
		C             During each reduction step, the rank-revealing QR
		C             factorization of a matrix stops when the estimated minimum
		C             singular value is smaller than TOL * MAX(SVLMAX,EMSV),
		C             where EMSV is the estimated maximum singular value.
		C             SVLMAX >= 0.
		C
		C     ABCD    (input/output) DOUBLE PRECISION array, dimension
		C             (LDABCD,M+N)
		C             On entry, the leading (N+P)-by-(M+N) part of this array
		C             must contain the compound input matrix of the system.
		C             On exit, the leading (NU+MU)-by-(M+NU) part of this array
		C             contains the reduced compound input matrix of the system.
		C
		C     LDABCD  INTEGER
		C             The leading dimension of array ABCD.
		C             LDABCD >= MAX(1,N+P).
		C
		C     NINFZ   (input/output) INTEGER
		C             On entry, the currently computed number of infinite zeros.
		C             It should be initialized to zero on the first call.
		C             NINFZ >= 0.
		C             On exit, the number of infinite zeros.
		C
		C     INFZ    (input/output) INTEGER array, dimension (N)
		C             On entry, INFZ(i) must contain the current number of
		C             infinite zeros of degree i, where i = 1,2,...,N, found in
		C             the previous call(s) of the routine. It should be
		C             initialized to zero on the first call.
		C             On exit, INFZ(i) contains the number of infinite zeros of
		C             degree i, where i = 1,2,...,N.
		C
		C     KRONL   (input/output) INTEGER array, dimension (N+1)
		C             On entry, this array must contain the currently computed
		C             left Kronecker (row) indices found in the previous call(s)
		C             of the routine. It should be initialized to zero on the
		C             first call.
		C             On exit, the leading NKROL elements of this array contain
		C             the left Kronecker (row) indices.
		C
		C     MU      (output) INTEGER
		C             The normal rank of the transfer function matrix of the
		C             original system.
		C
		C     NU      (output) INTEGER
		C             The dimension of the reduced system matrix and the number
		C             of (finite) invariant zeros if D' is invertible.
		C
		C     NKROL   (output) INTEGER
		C             The number of left Kronecker indices.
		C
		C     Tolerances
		C
		C     TOL     DOUBLE PRECISION
		C             A tolerance used in rank decisions to determine the
		C             effective rank, which is defined as the order of the
		C             largest leading (or trailing) triangular submatrix in the
		C             QR (or RQ) factorization with column (or row) pivoting
		C             whose estimated condition number is less than 1/TOL.
		C             NOTE that when SVLMAX > 0, the estimated ranks could be
		C             less than those defined above (see SVLMAX).
		C
		C     Workspace
		C
		C     IWORK   INTEGER array, dimension (MAX(M,P))
		C
		C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
		C             On exit, if INFO = 0, DWORK(1) returns the optimal value
		C             of LDWORK.
		C
		C     LDWORK  INTEGER
		C             The length of the array DWORK.
		C             LDWORK >= MAX( 1, MIN(P,M) + MAX(3*M-1,N),
		C                               MIN(P,N) + MAX(3*P-1,N+P,N+M) ).
		C             For optimum performance LDWORK should be larger.
		C
		C             If LDWORK = -1, then a workspace query is assumed;
		C             the routine only calculates the optimal size of the
		C             DWORK array, returns this value as the first entry of
		C             the DWORK array, and no error message related to LDWORK
		C             is issued by XERBLA.
		C
		C     Error Indicator
		C
		C     INFO    INTEGER
		C             = 0:  successful exit;
		C             < 0:  if INFO = -i, the i-th argument had an illegal
		C                   value.
		C
		C     REFERENCES
		C
		C     [1] Svaricek, F.
		C         Computation of the Structural Invariants of Linear
		C         Multivariable Systems with an Extended Version of
		C         the Program ZEROS.
		C         System & Control Letters, 6, pp. 261-266, 1985.
		C
		C     [2] Emami-Naeini, A. and Van Dooren, P.
		C         Computation of Zeros of Linear Multivariable Systems.
		C         Automatica, 18, pp. 415-430, 1982.
		C
		C     NUMERICAL ASPECTS
		C
		C     The algorithm is backward stable.
		C
		C     CONTRIBUTOR
		C
		C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996.
		C     Supersedes Release 2.0 routine AB08BZ by F. Svaricek.
		C
		C     REVISIONS
		C
		C     V. Sima, Oct. 1997; Feb. 1998, Jan. 2009, Apr. 2009, Apr. 2011.
		C     A. Varga, May 1999; May 2001.
		C
		C     KEYWORDS
		C
		C     Generalized eigenvalue problem, Kronecker indices, multivariable
		C     system, orthogonal transformation, structural invariant.
		C
		*/

		void AB08NX(int *, int *, int *, int *, int *, double *, double *, int *, int *, int *,
		            int *, int *, int *, int *, double *, int *, double *, int *, int *);

		/*
		SLICOT RELEASE 5.5.
		C
		C     Copyright (c) 2002-2012 NICONET e.V.
		C
		C     PURPOSE
		C
		C     To construct for a linear multivariable system described by a
		C     state-space model (A,B,C,D) a regular pencil (A - lambda*B ) which
		C                                                    f          f
		C     has the invariant zeros of the system as generalized eigenvalues.
		C     The routine also computes the orders of the infinite zeros and the
		C     right and left Kronecker indices of the system (A,B,C,D).
		C
		C     ARGUMENTS
		C
		C     Mode Parameters
		C
		C     EQUIL   CHARACTER*1
		C             Specifies whether the user wishes to balance the compound
		C             matrix (see METHOD) as follows:
		C             = 'S':  Perform balancing (scaling);
		C             = 'N':  Do not perform balancing.
		C
		C     Input/Output Parameters
		C
		C     N       (input) INTEGER
		C             The number of state variables, i.e., the order of the
		C             matrix A.  N >= 0.
		C
		C     M       (input) INTEGER
		C             The number of system inputs.  M >= 0.
		C
		C     P       (input) INTEGER
		C             The number of system outputs.  P >= 0.
		C
		C     A       (input) COMPLEX*16 array, dimension (LDA,N)
		C             The leading N-by-N part of this array must contain the
		C             state dynamics matrix A of the system.
		C
		C     LDA     INTEGER
		C             The leading dimension of array A.  LDA >= MAX(1,N).
		C
		C     B       (input) COMPLEX*16 array, dimension (LDB,M)
		C             The leading N-by-M part of this array must contain the
		C             input/state matrix B of the system.
		C
		C     LDB     INTEGER
		C             The leading dimension of array B.  LDB >= MAX(1,N).
		C
		C     C       (input) COMPLEX*16 array, dimension (LDC,N)
		C             The leading P-by-N part of this array must contain the
		C             state/output matrix C of the system.
		C
		C     LDC     INTEGER
		C             The leading dimension of array C.  LDC >= MAX(1,P).
		C
		C     D       (input) COMPLEX*16 array, dimension (LDD,M)
		C             The leading P-by-M part of this array must contain the
		C             direct transmission matrix D of the system.
		C
		C     LDD     INTEGER
		C             The leading dimension of array D.  LDD >= MAX(1,P).
		C
		C     NU      (output) INTEGER
		C             The number of (finite) invariant zeros.
		C
		C     RANK    (output) INTEGER
		C             The normal rank of the transfer function matrix.
		C
		C     DINFZ   (output) INTEGER
		C             The maximum degree of infinite elementary divisors.
		C
		C     NKROR   (output) INTEGER
		C             The number of right Kronecker indices.
		C
		C     NKROL   (output) INTEGER
		C             The number of left Kronecker indices.
		C
		C     INFZ    (output) INTEGER array, dimension (N)
		C             The leading DINFZ elements of INFZ contain information
		C             on the infinite elementary divisors as follows:
		C             the system has INFZ(i) infinite elementary divisors
		C             of degree i, where i = 1,2,...,DINFZ.
		C
		C     KRONR   (output) INTEGER array, dimension (MAX(N,M)+1)
		C             The leading NKROR elements of this array contain the
		C             right Kronecker (column) indices.
		C
		C     KRONL   (output) INTEGER array, dimension (MAX(N,P)+1)
		C             The leading NKROL elements of this array contain the
		C             left Kronecker (row) indices.
		C
		C     AF      (output) COMPLEX*16 array, dimension (LDAF,N+MIN(P,M))
		C             The leading NU-by-NU part of this array contains the
		C             coefficient matrix A  of the reduced pencil. The remainder
		C                                 f
		C             of the leading (N+M)-by-(N+MIN(P,M)) part is used as
		C             internal workspace.
		C
		C     LDAF    INTEGER
		C             The leading dimension of array AF.  LDAF >= MAX(1,N+M).
		C
		C     BF      (output) COMPLEX*16 array, dimension (LDBF,N+M)
		C             The leading NU-by-NU part of this array contains the
		C             coefficient matrix B  of the reduced pencil. The
		C                                 f
		C             remainder of the leading (N+P)-by-(N+M) part is used as
		C             internal workspace.
		C
		C     LDBF    INTEGER
		C             The leading dimension of array BF.  LDBF >= MAX(1,N+P).
		C
		C     Tolerances
		C
		C     TOL     DOUBLE PRECISION
		C             A tolerance used in rank decisions to determine the
		C             effective rank, which is defined as the order of the
		C             largest leading (or trailing) triangular submatrix in the
		C             QR (or RQ) factorization with column (or row) pivoting
		C             whose estimated condition number is less than 1/TOL.
		C             If the user sets TOL to be less than SQRT((N+P)*(N+M))*EPS
		C             then the tolerance is taken as SQRT((N+P)*(N+M))*EPS,
		C             where EPS is the machine precision (see LAPACK Library
		C             Routine DLAMCH).
		C
		C     Workspace
		C
		C     IWORK   INTEGER array, dimension (MAX(M,P))
		C
		C     DWORK   DOUBLE PRECISION array, dimension (MAX(N,2*MAX(P,M)))
		C
		C     ZWORK   DOUBLE PRECISION array, dimension (LZWORK)
		C             On exit, if INFO = 0, ZWORK(1) returns the optimal value
		C             of LZWORK.
		C
		C     LZWORK  INTEGER
		C             The length of the array ZWORK.
		C             LZWORK >= MAX( 1, MIN(P,M) + MAX(3*M-1,N),
		C                               MIN(P,N) + MAX(3*P-1,N+P,N+M),
		C                               MIN(M,N) + MAX(3*M-1,N+M) ).
		C             An upper bound is MAX(s,N) + MAX(3*s-1,N+s), with
		C             s = MAX(M,P).
		C             For optimum performance LZWORK should be larger.
		C
		C             If LZWORK = -1, then a workspace query is assumed;
		C             the routine only calculates the optimal size of the
		C             ZWORK array, returns this value as the first entry of
		C             the ZWORK array, and no error message related to LZWORK
		C             is issued by XERBLA.
		C
		C     Error Indicator
		C
		C     INFO    INTEGER
		C             = 0:  successful exit;
		C             < 0:  if INFO = -i, the i-th argument had an illegal
		C                   value.
		C
		C     METHOD
		C
		C     The routine extracts from the system matrix of a state-space
		C     system (A,B,C,D) a regular pencil A - lambda*B  which has the
		C                                        f          f
		C     invariant zeros of the system as generalized eigenvalues as
		C     follows:
		C
		C        (a) construct the (N+P)-by-(N+M) compound matrix (B  A);
		C                                                         (D  C)
		C
		C        (b) reduce the above system to one with the same invariant
		C            zeros and with D of full row rank;
		C
		C        (c) pertranspose the system;
		C
		C        (d) reduce the system to one with the same invariant zeros and
		C            with D square invertible;
		C
		C        (e) perform a unitary transformation on the columns of
		C            (A - lambda*I  B) in order to reduce it to
		C            (      C       D)
		C
		C            (A  - lambda*B   X)
		C            ( f           f   ), with Y and B  square invertible;
		C            (     0          Y)              f
		C
		C        (f) compute the right and left Kronecker indices of the system
		C            (A,B,C,D), which together with the orders of the infinite
		C            zeros (determined by steps (a) - (e)) constitute the
		C            complete set of structural invariants under strict
		C            equivalence transformations of a linear system.
		C
		C     REFERENCES
		C
		C     [1] Svaricek, F.
		C         Computation of the Structural Invariants of Linear
		C         Multivariable Systems with an Extended Version of
		C         the Program ZEROS.
		C         System & Control Letters, 6, pp. 261-266, 1985.
		C
		C     [2] Emami-Naeini, A. and Van Dooren, P.
		C         Computation of Zeros of Linear Multivariable Systems.
		C         Automatica, 18, pp. 415-430, 1982.
		C
		C     NUMERICAL ASPECTS
		C
		C     The algorithm is backward stable (see [2] and [1]).
		C
		C     FURTHER COMMENTS
		C
		C     In order to compute the invariant zeros of the system explicitly,
		C     a call to this routine may be followed by a call to the LAPACK
		C     Library routine ZGGEV with A = A , B = B  and N = NU.
		C                                     f       f
		C     If RANK = 0, the routine ZGEEV can be used (since B = I).
		C                                                        f
		C     CONTRIBUTOR
		C
		C     V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996.
		C     Complex version: V. Sima, Research Institute for Informatics,
		C     Bucharest, Nov. 2008.
		C
		C     REVISIONS
		C
		C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2009,
		C     Apr. 2009, Apr. 2011.
		C
		C     KEYWORDS
		C
		C     Generalized eigenvalue problem, Kronecker indices, multivariable
		C     system, unitary transformation, structural invariant.
		C
		*/


		void AB08NZ(unsigned char *, int *, int *, int *, std::complex<double> *, int *, std::complex<double> *, int *,
		            std::complex<double> *, int *, std::complex<double> *, int *, int *, int *, int *, int *, int *,
					int *, int *, int *, std::complex<double> *, int *, std::complex<double> *, int *, double *, int *,
					double *, std::complex<double> *, int *, int *);

		/*
		SLICOT RELEASE 5.5.
		C
		C     Copyright (c) 2002-2012 NICONET e.V.
		C
		C     PURPOSE
		C
		C     To compute a reduced order model (Ar,Br,Cr) for a stable original
		C     state-space representation (A,B,C) by using either the square-root
		C     or the balancing-free square-root Balance & Truncate (B & T)
		C     model reduction method.
		C
		C     ARGUMENTS
		C
		C     Mode Parameters
		C
		C     DICO    CHARACTER*1
		C             Specifies the type of the original system as follows:
		C             = 'C':  continuous-time system;
		C             = 'D':  discrete-time system.
		C
		C     JOB     CHARACTER*1
		C             Specifies the model reduction approach to be used
		C             as follows:
		C             = 'B':  use the square-root Balance & Truncate method;
		C             = 'N':  use the balancing-free square-root
		C                     Balance & Truncate method.
		C
		C     EQUIL   CHARACTER*1
		C             Specifies whether the user wishes to preliminarily
		C             equilibrate the triplet (A,B,C) as follows:
		C             = 'S':  perform equilibration (scaling);
		C             = 'N':  do not perform equilibration.
		C
		C     ORDSEL  CHARACTER*1
		C             Specifies the order selection method as follows:
		C             = 'F':  the resulting order NR is fixed;
		C             = 'A':  the resulting order NR is automatically determined
		C                     on basis of the given tolerance TOL.
		C
		C     Input/Output Parameters
		C
		C     N       (input) INTEGER
		C             The order of the original state-space representation, i.e.
		C             the order of the matrix A.  N >= 0.
		C
		C     M       (input) INTEGER
		C             The number of system inputs.  M >= 0.
		C
		C     P       (input) INTEGER
		C             The number of system outputs.  P >= 0.
		C
		C     NR      (input/output) INTEGER
		C             On entry with ORDSEL = 'F', NR is the desired order of the
		C             resulting reduced order system.  0 <= NR <= N.
		C             On exit, if INFO = 0, NR is the order of the resulting
		C             reduced order model. NR is set as follows:
		C             if ORDSEL = 'F', NR is equal to MIN(NR,NMIN), where NR
		C             is the desired order on entry and NMIN is the order of a
		C             minimal realization of the given system; NMIN is
		C             determined as the number of Hankel singular values greater
		C             than N*EPS*HNORM(A,B,C), where EPS is the machine
		C             precision (see LAPACK Library Routine DLAMCH) and
		C             HNORM(A,B,C) is the Hankel norm of the system (computed
		C             in HSV(1));
		C             if ORDSEL = 'A', NR is equal to the number of Hankel
		C             singular values greater than MAX(TOL,N*EPS*HNORM(A,B,C)).
		C
		C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
		C             On entry, the leading N-by-N part of this array must
		C             contain the state dynamics matrix A.
		C             On exit, if INFO = 0, the leading NR-by-NR part of this
		C             array contains the state dynamics matrix Ar of the reduced
		C             order system.
		C
		C     LDA     INTEGER
		C             The leading dimension of array A.  LDA >= MAX(1,N).
		C
		C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
		C             On entry, the leading N-by-M part of this array must
		C             contain the original input/state matrix B.
		C             On exit, if INFO = 0, the leading NR-by-M part of this
		C             array contains the input/state matrix Br of the reduced
		C             order system.
		C
		C     LDB     INTEGER
		C             The leading dimension of array B.  LDB >= MAX(1,N).
		C
		C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
		C             On entry, the leading P-by-N part of this array must
		C             contain the original state/output matrix C.
		C             On exit, if INFO = 0, the leading P-by-NR part of this
		C             array contains the state/output matrix Cr of the reduced
		C             order system.
		C
		C     LDC     INTEGER
		C             The leading dimension of array C.  LDC >= MAX(1,P).
		C
		C     HSV     (output) DOUBLE PRECISION array, dimension (N)
		C             If INFO = 0, it contains the Hankel singular values of
		C             the original system ordered decreasingly. HSV(1) is the
		C             Hankel norm of the system.
		C
		C     Tolerances
		C
		C     TOL     DOUBLE PRECISION
		C             If ORDSEL = 'A', TOL contains the tolerance for
		C             determining the order of reduced system.
		C             For model reduction, the recommended value is
		C             TOL = c*HNORM(A,B,C), where c is a constant in the
		C             interval [0.00001,0.001], and HNORM(A,B,C) is the
		C             Hankel-norm of the given system (computed in HSV(1)).
		C             For computing a minimal realization, the recommended
		C             value is TOL = N*EPS*HNORM(A,B,C), where EPS is the
		C             machine precision (see LAPACK Library Routine DLAMCH).
		C             This value is used by default if TOL <= 0 on entry.
		C             If ORDSEL = 'F', the value of TOL is ignored.
		C
		C     Workspace
		C
		C     IWORK   INTEGER array, dimension (LIWORK)
		C             LIWORK = 0, if JOB = 'B';
		C             LIWORK = N, if JOB = 'N'.
		C
		C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
		C             On exit, if INFO = 0, DWORK(1) returns the optimal value
		C             of LDWORK.
		C
		C     LDWORK  INTEGER
		C             The length of the array DWORK.
		C             LDWORK >= MAX(1,N*(2*N+MAX(N,M,P)+5)+N*(N+1)/2).
		C             For optimum performance LDWORK should be larger.
		C
		C     Warning Indicator
		C
		C     IWARN   INTEGER
		C             = 0:  no warning;
		C             = 1:  with ORDSEL = 'F', the selected order NR is greater
		C                   than the order of a minimal realization of the
		C                   given system. In this case, the resulting NR is
		C                   set automatically to a value corresponding to the
		C                   order of a minimal realization of the system.
		C
		C     Error Indicator
		C
		C     INFO    INTEGER
		C             = 0:  successful exit;
		C             < 0:  if INFO = -i, the i-th argument had an illegal
		C                   value;
		C             = 1:  the reduction of A to the real Schur form failed;
		C             = 2:  the state matrix A is not stable (if DICO = 'C')
		C                   or not convergent (if DICO = 'D');
		C             = 3:  the computation of Hankel singular values failed.
		C
		C     METHOD
		C
		C     Let be the stable linear system
		C
		C          d[x(t)] = Ax(t) + Bu(t)
		C          y(t)    = Cx(t)                               (1)
		C
		C     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1)
		C     for a discrete-time system. The subroutine AB09AD determines for
		C     the given system (1), the matrices of a reduced order system
		C
		C          d[z(t)] = Ar*z(t) + Br*u(t)
		C          yr(t)   = Cr*z(t)                             (2)
		C
		C     such that
		C
		C           HSV(NR) <= INFNORM(G-Gr) <= 2*[HSV(NR+1) + ... + HSV(N)],
		C
		C     where G and Gr are transfer-function matrices of the systems
		C     (A,B,C) and (Ar,Br,Cr), respectively, and INFNORM(G) is the
		C     infinity-norm of G.
		C
		C     If JOB = 'B', the square-root Balance & Truncate method of [1]
		C     is used and, for DICO = 'C', the resulting model is balanced.
		C     By setting TOL <= 0, the routine can be used to compute balanced
		C     minimal state-space realizations of stable systems.
		C
		C     If JOB = 'N', the balancing-free square-root version of the
		C     Balance & Truncate method [2] is used.
		C     By setting TOL <= 0, the routine can be used to compute minimal
		C     state-space realizations of stable systems.
		C
		C     REFERENCES
		C
		C     [1] Tombs M.S. and Postlethwaite I.
		C         Truncated balanced realization of stable, non-minimal
		C         state-space systems.
		C         Int. J. Control, Vol. 46, pp. 1319-1330, 1987.
		C
		C     [2] Varga A.
		C         Efficient minimal realization procedure based on balancing.
		C         Proc. of IMACS/IFAC Symp. MCTS, Lille, France, May 1991,
		C         A. El Moudui, P. Borne, S. G. Tzafestas (Eds.),
		C         Vol. 2, pp. 42-46.
		C
		C     NUMERICAL ASPECTS
		C
		C     The implemented methods rely on accuracy enhancing square-root or
		C     balancing-free square-root techniques.
		C                                         3
		C     The algorithms require less than 30N  floating point operations.
		C
		C     CONTRIBUTOR
		C
		C     C. Oara and A. Varga, German Aerospace Center,
		C     DLR Oberpfaffenhofen, March 1998.
		C     Based on the RASP routines SRBT and SRBFT.
		C
		C     REVISIONS
		C
		C     May 2, 1998.
		C     November 11, 1998, V. Sima, Research Institute for Informatics,
		C     Bucharest.
		C
		C     KEYWORDS
		C
		C     Balancing, minimal state-space representation, model reduction,
		C     multivariable system, state-space model.
		*/


		void AB09AD(unsigned  char *, unsigned char *, unsigned char *, unsigned char *,
		            int *, int *, int *, int *, double *, int *, double *, int *, 
					double *, int *, double *, double *, int *, double *, int *, int *, int *);

		/*
		SLICOT RELEASE 5.5.
		C
		C     Copyright (c) 2002-2012 NICONET e.V.
		C
		C     PURPOSE
		C
		C     To compute a reduced order model (Ar,Br,Cr) for a stable original
		C     state-space representation (A,B,C) by using either the square-root
		C     or the balancing-free square-root Balance & Truncate model
		C     reduction method. The state dynamics matrix A of the original
		C     system is an upper quasi-triangular matrix in real Schur canonical
		C     form. The matrices of the reduced order system are computed using
		C     the truncation formulas:
		C
		C          Ar = TI * A * T ,  Br = TI * B ,  Cr = C * T .
		C
		C     ARGUMENTS
		C
		C     Mode Parameters
		C
		C     DICO    CHARACTER*1
		C             Specifies the type of the original system as follows:
		C             = 'C':  continuous-time system;
		C             = 'D':  discrete-time system.
		C
		C     JOB     CHARACTER*1
		C             Specifies the model reduction approach to be used
		C             as follows:
		C             = 'B':  use the square-root Balance & Truncate method;
		C             = 'N':  use the balancing-free square-root
		C                     Balance & Truncate method.
		C
		C     ORDSEL  CHARACTER*1
		C             Specifies the order selection method as follows:
		C             = 'F':  the resulting order NR is fixed;
		C             = 'A':  the resulting order NR is automatically determined
		C                     on basis of the given tolerance TOL.
		C
		C     Input/Output Parameters
		C
		C     N       (input) INTEGER
		C             The order of the original state-space representation, i.e.
		C             the order of the matrix A.  N >= 0.
		C
		C     M       (input) INTEGER
		C             The number of system inputs.  M >= 0.
		C
		C     P       (input) INTEGER
		C             The number of system outputs.  P >= 0.
		C
		C     NR      (input/output) INTEGER
		C             On entry with ORDSEL = 'F', NR is the desired order of the
		C             resulting reduced order system.  0 <= NR <= N.
		C             On exit, if INFO = 0, NR is the order of the resulting
		C             reduced order model. NR is set as follows:
		C             if ORDSEL = 'F', NR is equal to MIN(NR,NMIN), where NR
		C             is the desired order on entry and NMIN is the order of a
		C             minimal realization of the given system; NMIN is
		C             determined as the number of Hankel singular values greater
		C             than N*EPS*HNORM(A,B,C), where EPS is the machine
		C             precision (see LAPACK Library Routine DLAMCH) and
		C             HNORM(A,B,C) is the Hankel norm of the system (computed
		C             in HSV(1));
		C             if ORDSEL = 'A', NR is equal to the number of Hankel
		C             singular values greater than MAX(TOL,N*EPS*HNORM(A,B,C)).
		C
		C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
		C             On entry, the leading N-by-N part of this array must
		C             contain the state dynamics matrix A in a real Schur
		C             canonical form.
		C             On exit, if INFO = 0, the leading NR-by-NR part of this
		C             array contains the state dynamics matrix Ar of the
		C             reduced order system.
		C
		C     LDA     INTEGER
		C             The leading dimension of array A.  LDA >= MAX(1,N).
		C
		C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
		C             On entry, the leading N-by-M part of this array must
		C             contain the original input/state matrix B.
		C             On exit, if INFO = 0, the leading NR-by-M part of this
		C             array contains the input/state matrix Br of the reduced
		C             order system.
		C
		C     LDB     INTEGER
		C             The leading dimension of array B.  LDB >= MAX(1,N).
		C
		C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
		C             On entry, the leading P-by-N part of this array must
		C             contain the original state/output matrix C.
		C             On exit, if INFO = 0, the leading P-by-NR part of this
		C             array contains the state/output matrix Cr of the reduced
		C             order system.
		C
		C     LDC     INTEGER
		C             The leading dimension of array C.  LDC >= MAX(1,P).
		C
		C     HSV     (output) DOUBLE PRECISION array, dimension (N)
		C             If INFO = 0, it contains the Hankel singular values of
		C             the original system ordered decreasingly. HSV(1) is the
		C             Hankel norm of the system.
		C
		C     T       (output) DOUBLE PRECISION array, dimension (LDT,N)
		C             If INFO = 0 and NR > 0, the leading N-by-NR part of this
		C             array contains the right truncation matrix T.
		C
		C     LDT     INTEGER
		C             The leading dimension of array T.  LDT >= MAX(1,N).
		C
		C     TI      (output) DOUBLE PRECISION array, dimension (LDTI,N)
		C             If INFO = 0 and NR > 0, the leading NR-by-N part of this
		C             array contains the left truncation matrix TI.
		C
		C     LDTI    INTEGER
		C             The leading dimension of array TI.  LDTI >= MAX(1,N).
		C
		C     Tolerances
		C
		C     TOL     DOUBLE PRECISION
		C             If ORDSEL = 'A', TOL contains the tolerance for
		C             determining the order of reduced system.
		C             For model reduction, the recommended value is
		C             TOL = c*HNORM(A,B,C), where c is a constant in the
		C             interval [0.00001,0.001], and HNORM(A,B,C) is the
		C             Hankel-norm of the given system (computed in HSV(1)).
		C             For computing a minimal realization, the recommended
		C             value is TOL = N*EPS*HNORM(A,B,C), where EPS is the
		C             machine precision (see LAPACK Library Routine DLAMCH).
		C             This value is used by default if TOL <= 0 on entry.
		C             If ORDSEL = 'F', the value of TOL is ignored.
		C
		C     Workspace
		C
		C     IWORK   INTEGER array, dimension (LIWORK)
		C             LIWORK = 0, if JOB = 'B', or
		C             LIWORK = N, if JOB = 'N'.
		C
		C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
		C             On exit, if INFO = 0, DWORK(1) returns the optimal value
		C             of LDWORK.
		C
		C     LDWORK  INTEGER
		C             The length of the array DWORK.
		C             LDWORK >= MAX(1,N*(MAX(N,M,P)+5) + N*(N+1)/2).
		C             For optimum performance LDWORK should be larger.
		C
		C     Warning Indicator
		C
		C     IWARN   INTEGER
		C             = 0:  no warning;
		C             = 1:  with ORDSEL = 'F', the selected order NR is greater
		C                   than the order of a minimal realization of the
		C                   given system. In this case, the resulting NR is
		C                   set automatically to a value corresponding to the
		C                   order of a minimal realization of the system.
		C
		C     Error Indicator
		C
		C     INFO    INTEGER
		C             = 0:  successful exit;
		C             < 0:  if INFO = -i, the i-th argument had an illegal
		C                   value;
		C             = 1:  the state matrix A is not stable (if DICO = 'C')
		C                   or not convergent (if DICO = 'D');
		C             = 2:  the computation of Hankel singular values failed.
		C
		C     METHOD
		C
		C     Let be the stable linear system
		C
		C          d[x(t)] = Ax(t) + Bu(t)
		C          y(t)    = Cx(t)                               (1)
		C
		C     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1)
		C     for a discrete-time system. The subroutine AB09AX determines for
		C     the given system (1), the matrices of a reduced NR order system
		C
		C          d[z(t)] = Ar*z(t) + Br*u(t)
		C          yr(t)   = Cr*z(t)                             (2)
		C
		C     such that
		C
		C           HSV(NR) <= INFNORM(G-Gr) <= 2*[HSV(NR+1) + ... + HSV(N)],
		C
		C     where G and Gr are transfer-function matrices of the systems
		C     (A,B,C) and (Ar,Br,Cr), respectively, and INFNORM(G) is the
		C     infinity-norm of G.
		C
		C     If JOB = 'B', the square-root Balance & Truncate method of [1]
		C     is used and, for DICO = 'C', the resulting model is balanced.
		C     By setting TOL <= 0, the routine can be used to compute balanced
		C     minimal state-space realizations of stable systems.
		C
		C     If JOB = 'N', the balancing-free square-root version of the
		C     Balance & Truncate method [2] is used.
		C     By setting TOL <= 0, the routine can be used to compute minimal
		C     state-space realizations of stable systems.
		C
		C     REFERENCES
		C
		C     [1] Tombs M.S. and Postlethwaite I.
		C         Truncated balanced realization of stable, non-minimal
		C         state-space systems.
		C         Int. J. Control, Vol. 46, pp. 1319-1330, 1987.
		C
		C     [2] Varga A.
		C         Efficient minimal realization procedure based on balancing.
		C         Proc. of IMACS/IFAC Symp. MCTS, Lille, France, May 1991,
		C         A. El Moudui, P. Borne, S. G. Tzafestas (Eds.),
		C         Vol. 2, pp. 42-46.
		C
		C     NUMERICAL ASPECTS
		C
		C     The implemented methods rely on accuracy enhancing square-root or
		C     balancing-free square-root techniques.
		C                                         3
		C     The algorithms require less than 30N  floating point operations.
		C
		C     CONTRIBUTOR
		C
		C     A. Varga, German Aerospace Center,
		C     DLR Oberpfaffenhofen, March 1998.
		C     Based on the RASP routines SRBT1 and SRBFT1.
		C
		C     REVISIONS
		C
		C     May 2, 1998.
		C     November 11, 1998, V. Sima, Research Institute for Informatics,
		C     Bucharest.
		C     December 1998, V. Sima, Katholieke Univ. Leuven, Leuven.
		C     February 14, 1999, A. Varga, German Aerospace Center.
		C     February 22, 1999, V. Sima, Research Institute for Informatics.
		C     February 27, 2000, V. Sima, Research Institute for Informatics.
		C
		C     KEYWORDS
		C
		C     Balancing, minimal state-space representation, model reduction,
		C     multivariable system, state-space model.
		*/

		void AB09AX(unsigned char *, unsigned char *, unsigned char *, int *, int *, int *, int *,
		            double *, int *, double *, int *, double *, int *, double *, double *, int *, 
					double *, int *, double *, int *, double *, int *, int *, int *);

		/*
		SLICOT RELEASE 5.5.
		C
		C     Copyright (c) 2002-2012 NICONET e.V.
		C
		C     PURPOSE
		C
		C     To compute a reduced order model (Ar,Br,Cr,Dr) for a stable
		C     original state-space representation (A,B,C,D) by using either the
		C     square-root or the balancing-free square-root Singular
		C     Perturbation Approximation (SPA) model reduction method.
		C
		C     ARGUMENTS
		C
		C     Mode Parameters
		C
		C     DICO    CHARACTER*1
		C             Specifies the type of the original system as follows:
		C             = 'C':  continuous-time system;
		C             = 'D':  discrete-time system.
		C
		C     JOB     CHARACTER*1
		C             Specifies the model reduction approach to be used
		C             as follows:
		C             = 'B':  use the square-root SPA method;
		C             = 'N':  use the balancing-free square-root SPA method.
		C
		C     EQUIL   CHARACTER*1
		C             Specifies whether the user wishes to preliminarily
		C             equilibrate the triplet (A,B,C) as follows:
		C             = 'S':  perform equilibration (scaling);
		C             = 'N':  do not perform equilibration.
		C
		C     ORDSEL  CHARACTER*1
		C             Specifies the order selection method as follows:
		C             = 'F':  the resulting order NR is fixed;
		C             = 'A':  the resulting order NR is automatically determined
		C                     on basis of the given tolerance TOL1.
		C
		C     Input/Output Parameters
		C
		C     N       (input) INTEGER
		C             The order of the original state-space representation, i.e.
		C             the order of the matrix A.  N >= 0.
		C
		C     M       (input) INTEGER
		C             The number of system inputs.  M >= 0.
		C
		C     P       (input) INTEGER
		C             The number of system outputs.  P >= 0.
		C
		C     NR      (input/output) INTEGER
		C             On entry with ORDSEL = 'F', NR is the desired order of
		C             the resulting reduced order system.  0 <= NR <= N.
		C             On exit, if INFO = 0, NR is the order of the resulting
		C             reduced order model. NR is set as follows:
		C             if ORDSEL = 'F', NR is equal to MIN(NR,NMIN), where NR
		C             is the desired order on entry and NMIN is the order of a
		C             minimal realization of the given system; NMIN is
		C             determined as the number of Hankel singular values greater
		C             than N*EPS*HNORM(A,B,C), where EPS is the machine
		C             precision (see LAPACK Library Routine DLAMCH) and
		C             HNORM(A,B,C) is the Hankel norm of the system (computed
		C             in HSV(1));
		C             if ORDSEL = 'A', NR is equal to the number of Hankel
		C             singular values greater than MAX(TOL1,N*EPS*HNORM(A,B,C)).
		C
		C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
		C             On entry, the leading N-by-N part of this array must
		C             contain the state dynamics matrix A.
		C             On exit, if INFO = 0, the leading NR-by-NR part of this
		C             array contains the state dynamics matrix Ar of the
		C             reduced order system.
		C
		C     LDA     INTEGER
		C             The leading dimension of array A.  LDA >= MAX(1,N).
		C
		C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
		C             On entry, the leading N-by-M part of this array must
		C             contain the original input/state matrix B.
		C             On exit, if INFO = 0, the leading NR-by-M part of this
		C             array contains the input/state matrix Br of the reduced
		C             order system.
		C
		C     LDB     INTEGER
		C             The leading dimension of array B.  LDB >= MAX(1,N).
		C
		C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
		C             On entry, the leading P-by-N part of this array must
		C             contain the original state/output matrix C.
		C             On exit, if INFO = 0, the leading P-by-NR part of this
		C             array contains the state/output matrix Cr of the reduced
		C             order system.
		C
		C     LDC     INTEGER
		C             The leading dimension of array C.  LDC >= MAX(1,P).
		C
		C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
		C             On entry, the leading P-by-M part of this array must
		C             contain the original input/output matrix D.
		C             On exit, if INFO = 0, the leading P-by-M part of this
		C             array contains the input/output matrix Dr of the reduced
		C             order system.
		C
		C     LDD     INTEGER
		C             The leading dimension of array D.  LDD >= MAX(1,P).
		C
		C     HSV     (output) DOUBLE PRECISION array, dimension (N)
		C             If INFO = 0, it contains the Hankel singular values of
		C             the original system ordered decreasingly. HSV(1) is the
		C             Hankel norm of the system.
		C
		C     Tolerances
		C
		C     TOL1    DOUBLE PRECISION
		C             If ORDSEL = 'A', TOL1 contains the tolerance for
		C             determining the order of reduced system.
		C             For model reduction, the recommended value is
		C             TOL1 = c*HNORM(A,B,C), where c is a constant in the
		C             interval [0.00001,0.001], and HNORM(A,B,C) is the
		C             Hankel-norm of the given system (computed in HSV(1)).
		C             For computing a minimal realization, the recommended
		C             value is TOL1 = N*EPS*HNORM(A,B,C), where EPS is the
		C             machine precision (see LAPACK Library Routine DLAMCH).
		C             This value is used by default if TOL1 <= 0 on entry.
		C             If ORDSEL = 'F', the value of TOL1 is ignored.
		C
		C     TOL2    DOUBLE PRECISION
		C             The tolerance for determining the order of a minimal
		C             realization of the given system. The recommended value is
		C             TOL2 = N*EPS*HNORM(A,B,C). This value is used by default
		C             if TOL2 <= 0 on entry.
		C             If TOL2 > 0, then TOL2 <= TOL1.
		C
		C     Workspace
		C
		C     IWORK   INTEGER array, dimension MAX(1,2*N)
		C             On exit with INFO = 0, IWORK(1) contains the order of the
		C             minimal realization of the system.
		C
		C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
		C             On exit, if INFO = 0, DWORK(1) returns the optimal value
		C             of LDWORK.
		C
		C     LDWORK  INTEGER
		C             The length of the array DWORK.
		C             LDWORK >= MAX(1,N*(2*N+MAX(N,M,P)+5)+N*(N+1)/2).
		C             For optimum performance LDWORK should be larger.
		C
		C     Warning Indicator
		C
		C     IWARN   INTEGER
		C             = 0:  no warning;
		C             = 1:  with ORDSEL = 'F', the selected order NR is greater
		C                   than the order of a minimal realization of the
		C                   given system. In this case, the resulting NR is
		C                   set automatically to a value corresponding to the
		C                   order of a minimal realization of the system.
		C
		C     Error Indicator
		C
		C     INFO    INTEGER
		C             = 0:  successful exit;
		C             < 0:  if INFO = -i, the i-th argument had an illegal
		C                   value;
		C             = 1:  the reduction of A to the real Schur form failed;
		C             = 2:  the state matrix A is not stable (if DICO = 'C')
		C                   or not convergent (if DICO = 'D');
		C             = 3:  the computation of Hankel singular values failed.
		C
		C     METHOD
		C
		C     Let be the stable linear system
		C
		C          d[x(t)] = Ax(t) + Bu(t)
		C          y(t)    = Cx(t) + Du(t)                           (1)
		C
		C     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1)
		C     for a discrete-time system. The subroutine AB09BD determines for
		C     the given system (1), the matrices of a reduced order system
		C
		C          d[z(t)] = Ar*z(t) + Br*u(t)
		C          yr(t)   = Cr*z(t) + Dr*u(t)                       (2)
		C
		C     such that
		C
		C           HSV(NR) <= INFNORM(G-Gr) <= 2*[HSV(NR+1) + ... + HSV(N)],
		C
		C     where G and Gr are transfer-function matrices of the systems
		C     (A,B,C,D) and (Ar,Br,Cr,Dr), respectively, and INFNORM(G) is the
		C     infinity-norm of G.
		C
		C     If JOB = 'B', the balancing-based square-root SPA method of [1]
		C     is used and the resulting model is balanced.
		C
		C     If JOB = 'N', the balancing-free square-root SPA method of [2]
		C     is used.
		C     By setting TOL1 = TOL2, the routine can be used to compute
		C     Balance & Truncate approximations.
		C
		C     REFERENCES
		C
		C     [1] Liu Y. and Anderson B.D.O.
		C         Singular Perturbation Approximation of Balanced Systems,
		C         Int. J. Control, Vol. 50, pp. 1379-1405, 1989.
		C
		C     [2] Varga A.
		C         Balancing-free square-root algorithm for computing singular
		C         perturbation approximations.
		C         Proc. 30-th IEEE CDC,  Brighton, Dec. 11-13, 1991,
		C         Vol. 2, pp. 1062-1065.
		C
		C     NUMERICAL ASPECTS
		C
		C     The implemented methods rely on accuracy enhancing square-root or
		C     balancing-free square-root techniques.
		C                                         3
		C     The algorithms require less than 30N  floating point operations.
		C
		C     CONTRIBUTOR
		C
		C     C. Oara and A. Varga, German Aerospace Center,
		C     DLR Oberpfaffenhofen, March 1998.
		C     Based on the RASP routine SRBFSP.
		C
		C     REVISIONS
		C
		C     May 2, 1998.
		C     November 11, 1998, V. Sima, Research Institute for Informatics,
		C     Bucharest.
		C
		C     KEYWORDS
		C
		C     Balancing, minimal state-space representation, model reduction,
		C     multivariable system, singular perturbation approximation,
		C     state-space model.
		*/

		void AB09BD(unsigned char *, unsigned char *, unsigned char *, unsigned char *, int *,
		            int *, int *, int *, double *, int *, double *, int *, double *, int *, 
					double *, int *, double *, double *, double *, int *, double *, int *, int *, int *);


/*
SUBROUTINE MB3OYZ( M, N, A, LDA, RCOND, SVLMAX, RANK, SVAL, JPVT, &
     TAU, DWORK, ZWORK, INFO )
    INTEGER            INFO, LDA, M, N, RANK
      DOUBLE PRECISION   RCOND, SVLMAX
!C     .. Array Arguments ..

#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER            JPVT( * )
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: JPVT
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      COMPLEX*16         A( LDA, * ), TAU( * ), ZWORK( * )
#else
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: A
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: TAU
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: ZWORK
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION   DWORK( * ), SVAL( 3 )
#else
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(3) :: SVAL
#endif
*/
void MB30YZ(int *, int *, std::complex<double> __restrict *,
            int *, double *, double *, int *,
	    double __restrict *, int __restrict *,
	    std::complex<double> __restrict *, double __restrict *,
	    std::complex<double> __restrict *, int *);

/*
SUBROUTINE MB3PYZ( M, N, A, LDA, RCOND, SVLMAX, RANK, SVAL, JPVT,
  TAU, DWORK, ZWORK, INFO )
    INTEGER            INFO, LDA, M, N, RANK
      DOUBLE PRECISION   RCOND, SVLMAX
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER            JPVT( * )
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: JPVT
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1      
      COMPLEX*16         A( LDA, * ), TAU( * ), ZWORK( * )
#else
       COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: A
       COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: TAU
       COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: ZWORK
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
       DOUBLE PRECISION   DWORK( * ), SVAL( 3 )
#else      
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
       DOUBLE PRECISION, DIMENSION(3) :: SVAL
#endif
*/
void MB30YZ(int *,
            int *,
	    std::complex<double> __restrict *,
	    int *,
	    double *,
	    double *,
	    int *,
	    double __restrict *,
	    int __restrict *,
	    std::coomplex<double> __restrict *,
	    double __restrict *,
	    std::complex<double> __restrict *,
	    int *);

/*
SUBROUTINE FB01QD( JOBK, MULTBQ, N, M, P, S, LDS, A, LDA, B,  &
                        LDB, Q, LDQ, C, LDC, R, LDR, K, LDK, TOL, &
                        IWORK, DWORK, LDWORK, INFO )
    CHARACTER         JOBK, MULTBQ
      INTEGER           INFO, LDA, LDB, LDC, LDK, LDQ, LDR, LDS, LDWORK, &
                        M, N, P
      DOUBLE PRECISION  TOL
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER           IWORK(*)
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), &
           K(LDK,*), Q(LDQ,*), R(LDR,*), S(LDS,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: K
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Q
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: S
*/
void FB01QD(const char *,
            const char *,
	    int *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
            int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE FB01SD( JOBX, MULTAB, MULTRC, N, M, P, SINV, LDSINV, &
          AINV, LDAINV, B, LDB, RINV, LDRINV, C, LDC,           &
          QINV, LDQINV, X, RINVY, Z, E, TOL, IWORK,             &
          DWORK, LDWORK, INFO )
   CHARACTER         JOBX, MULTAB, MULTRC
      INTEGER           INFO, LDAINV, LDB, LDC, LDQINV, LDRINV, LDSINV, &
                        LDWORK, M, N, P
      DOUBLE PRECISION  TOL
!C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER           IWORK(*)
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  AINV(LDAINV,*), B(LDB,*), C(LDC,*), DWORK(*), &
                        E(*), QINV(LDQINV,*), RINV(LDRINV,*), RINVY(*), &
                        SINV(LDSINV,*), X(*), Z(*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AINV
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: E
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: QINV
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: RINV
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: RINVY
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: SINV
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: X
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: Z
#endif
*/
void FB01SD(const char *,
            const char *,
	    const char *,
	    int *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE FB01VD( N, M, L, P, LDP, A, LDA, B, LDB, C, LDC, Q, &
  LDQ, R, LDR, K, LDK, TOL, IWORK, DWORK, LDWORK,INFO)
    .. Scalar Arguments ..
      INTEGER           INFO, L, LDA, LDB, LDC, LDK, LDP, LDQ, LDR, &
                        LDWORK, M, N
      DOUBLE PRECISION  TOL
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER           IWORK(*)
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1      
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), &
           K(LDK,*), P(LDP,*), Q(LDQ,*), R(LDR,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: K
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: P
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Q
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
#endif
*/
void FB01VD(int *,
            int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
 SUBROUTINE MB02OD( SIDE, UPLO, TRANS, DIAG, NORM, M, N, ALPHA, A,
   LDA, B, LDB, RCOND, TOL, IWORK, DWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER          DIAG, NORM, SIDE, TRANS, UPLO
      INTEGER            INFO, LDA, LDB, M, N
      DOUBLE PRECISION   ALPHA, RCOND, TOL
!C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER            IWORK(*)
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION   A(LDA,*), B(LDB,*), DWORK(*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
#endif
*/
void MB02OD(const char *,
            const char *,
	    const char *,
            const char *,
            const char *,
	    int *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *);

/*
 SUBROUTINE FB01TD( JOBX, MULTRC, N, M, P, SINV, LDSINV, AINV, &
        LDAINV, AINVB, LDAINB, RINV, LDRINV, C, LDC,          &
        QINV, LDQINV, X, RINVY, Z, E, TOL, IWORK,             &
        DWORK, LDWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         JOBX, MULTRC
      INTEGER           INFO, LDAINB, LDAINV, LDC, LDQINV, LDRINV, &
                        LDSINV, LDWORK, M, N, P
      DOUBLE PRECISION  TOL
!C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER           IWORK(*)
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  AINV(LDAINV,*), AINVB(LDAINB,*), C(LDC,*), &
                       DWORK(*), E(*), QINV(LDQINV,*), RINV(LDRINV,*), &
                       RINVY(*), SINV(LDSINV,*), X(*), Z(*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AINV
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AINVB
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: E
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: QINV
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: RINV
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: RINVY
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: SINV
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: X
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: Z
#endif
*/
void FB01TD(const char *,
            const char *,
	    int *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);



/*
SUBROUTINE MB04KD( UPLO, N, M, P, R, LDR, A, LDA, B, LDB, C, LDC, &
     TAU, DWORK )
!C     .. Scalar Arguments ..
      CHARACTER         UPLO
      INTEGER           LDA, LDB, LDC, LDR, M, N, P
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), &
           R(LDR,*), TAU(*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
#endif
*/

void MB04KD(const char *,
            int *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    double __restrict *);

/*
 SUBROUTINE MB04ID(N, M, P, L, A, LDA, B, LDB, TAU, DWORK, LDWORK,INFO)
!C     .. Scalar Arguments ..
      INTEGER           INFO, L, LDA, LDB, LDWORK, M, N, P
!C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), DWORK(*), TAU(*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
#endif
*/
void MB04ID(int *,
            int *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE MB04LD( UPLO, N, M, P, L, LDL, A, LDA, B, LDB, C, LDC &
                  TAU, DWORK )
!C     .. Scalar Arguments ..
      CHARACTER         UPLO
      INTEGER           LDA, LDB, LDC, LDL, M, N, P
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), &
           L(LDL,*), TAU(*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: L
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
#endif
*/
void MB04LD(const char *,
            int *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    double __restrict *);

/*
SUBROUTINE SB01BD( DICO, N, M, NP, ALPHA, A, LDA, B, LDB, WR, WI, &
   NFP, NAP, NUP, F, LDF, Z, LDZ, TOL, DWORK, &
   LDWORK, IWARN, INFO )
!C     .. Scalar Arguments ..
      CHARACTER        DICO
      INTEGER          INFO, IWARN, LDA, LDB, LDF, LDWORK, LDZ, M, N, &
                       NAP, NFP, NP, NUP
      DOUBLE PRECISION ALPHA, TOL
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
     DOUBLE PRECISION A(LDA,*), B(LDB,*), DWORK(*), F(LDF,*), &
          WI(*), WR(*), Z(LDZ,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: F
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WI
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WR
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Z
#endif
*/
void SB01BD(const char *,
            int *,
	    int *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    double __restrict *,
	    int *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    int *,
	    int *);

/*
 SUBROUTINE SB01BX(REIG,N,XR,XI,WR,WI,S,P)
 !C     .. Scalar Arguments ..
      LOGICAL          REIG
      INTEGER          N
      DOUBLE PRECISION P, S, XI ,XR
      !C     .. Array Arguments ..

      DOUBLE PRECISION WI(*), WR(*)
*/
void SB01BD(int *,
            int *,
	    double *,
	    double *,
	    double __restrict *,
	    double __restrict *,
	    double *,
	    double *);

/*
 SUBROUTINE SB01BY(N,M,S,P,A B,F,TOL,DWORK,INFO)
!C     .. Scalar Arguments ..
      INTEGER           INFO, M, N
      DOUBLE PRECISION  P, S, TOL
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      !DOUBLE PRECISION  A(N,*), B(N,*), DWORK(*), F(M,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: F
#endif
*/
void SB01BY(int *,
            int *,
	    double *,
	    double *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double *,
	    double __restrict *,
	    int *);

/*
  SUBROUTINE MB03QD( DICO, STDOM, JOBU, N, NLOW, NSUP, ALPHA, &
       A, LDA, U, LDU, NDIM, DWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER        DICO, JOBU, STDOM
      INTEGER          INFO, LDA, LDU, N, NDIM, NLOW, NSUP
      DOUBLE PRECISION ALPHA
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION A(LDA,*), DWORK(*), U(LDU,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
#endif
*/
void MB03QD(const char *,
            const char *,
	    const char *,
	    int *,
	    int *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int *,
	    double __restrict *,
	    int *);

/*
 SUBROUTINE MB03QY(N,L,A,LDA,U,LDU,E1,E2,INFO)
!C     .. Scalar Arguments ..
      INTEGER          INFO, L, LDA, LDU, N
      DOUBLE PRECISION E1, E2
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION A(LDA,*), U(LDU,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
#endif      
*/
void MB03QY(int *,
            int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double *,
	    int *);

/*
SUBROUTINE TB04CD( JOBD, EQUIL, N, M, P, NPZ, A, LDA, B, LDB, C,  &
LDC, D, LDD, NZ, LDNZ, NP, LDNP, ZEROSR,                          &
ZEROSI, POLESR, POLESI, GAINS, LDGAIN, TOL,                       &
IWORK, DWORK, LDWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER          EQUIL, JOBD
      DOUBLE PRECISION   TOL
      INTEGER            INFO, LDA, LDB, LDC, LDD, LDGAIN, LDNP, LDNZ, &
                         LDWORK, M, N, NPZ, P
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION   A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*), &
                         DWORK(*), GAINS(LDGAIN,*), POLESI(*),   &
                         POLESR(*), ZEROSI(*), ZEROSR(*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: GAINS
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: POLESI
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: POLESR
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ZEROSI
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ZEROSR
#endif    
*/
void TB04CD(const char *,
            const char *,
	    int *,
	    int *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int __restrict *,
	    int *,
	    int __restrict *,
	    int *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    int *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
 SUBROUTINE TB01ID( JOB, N, M, P, MAXRED, A, LDA, B, LDB, C, LDC,  &
      SCALE, INFO)
!C     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            INFO, LDA, LDB, LDC, M, N, P
      DOUBLE PRECISION   MAXRED
!C     ..
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ), &
           SCALE( * )
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: SCALE
#endif
*/
void TB01ID(const char *,
            int *,
	    int *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *);

/*
SUBROUTINE TB01ZD( JOBZ, N, P, A, LDA, B, C, LDC, NCONT, Z, LDZ,
    TAU, TOL, DWORK, LDWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         JOBZ
      INTEGER           INFO, LDA, LDC, LDWORK, LDZ, N, NCONT, P
      DOUBLE PRECISION  TOL
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), B(*), C(LDC,*), DWORK(*), TAU(*), &
           Z(LDZ,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TAU
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Z
#endif      
*/
void TB01ZD(const char *,
            int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    double __restrict *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    double *,
	    double __restrict *,
	    int *,
	    int *);

/*
 SUBROUTINE MB01PD( SCUN, TYPE, M, N, KL, KU, ANRM, NBL, NROWS, A, &
   LDA, INFO )
!C     .. Scalar Arguments ..
      CHARACTER          SCUN, TYPE
      INTEGER            INFO, KL, KU, LDA, M, MN, N, NBL
      DOUBLE PRECISION   ANRM
      !C     .. Array Arguments ..
   
      INTEGER            NROWS ( * )
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION   A( LDA, * )
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
#endif
*/
void MB01PD(const char *,
            const char *,
	    int *,
	    int *,
	    int *,
	    int *,
	    double *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    int * );

/*
 SUBROUTINE MB01QD( TYPE, M, N, KL, KU, CFROM, CTO, NBL, NROWS, A,
   LDA, INFO )
!C     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N, NBL
      DOUBLE PRECISION   CFROM, CTO
!C     ..
      !C     .. Array Arguments ..

      INTEGER            NROWS ( * )
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1      
      DOUBLE PRECISION   A( LDA, * )
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
#endif
*/
void MB01QD(const char *,
            int *,
	    int *,
	    int *,
	    int *,
	    double *,
	    double *,
	    int *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE TB05AD( BALEIG, INITA, N, M, P, FREQ, A, LDA, B, LDB, &
  C, LDC, RCOND, G, LDG, EVRE, EVIM, HINVB,                      &
  LDHINV, IWORK, DWORK, LDWORK, ZWORK, LZWORK,INFO)
!C     .. Scalar Arguments ..
      CHARACTER         BALEIG, INITA
      INTEGER           INFO, LDA, LDB, LDC, LDG, LDHINV, LDWORK, &
                        LZWORK, M, N, P
      DOUBLE PRECISION  RCOND
      COMPLEX*16        FREQ
      !C     .. Array Arguments ..

      INTEGER           IWORK(*)
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), EVIM(*), &
           EVRE(*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: EVIM
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: EVRE
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      COMPLEX*16        ZWORK(*), G(LDG,*), HINVB(LDHINV,*)
#else
      COMPLEX(16), DIMENSION(:), ALLOCATABLE :: ZWORK
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: G
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: HINVB
#endif
*/
void TB05AD(const char *,
            const char *,
	    int *,
	    int *,
	    int *,
	    std::complex<double> *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    std::complex<double> __restrict *,
	    int *,
	    double __restrict *,
	    double __restrict *,
	    std::complex<double> __restrict *,
	    int *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    std::complex<double> __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE MB02RZ(TRANS,N,NRHS,H,LDH,IPIV,B,LDB,INFO)
!C     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDB, LDH, N, NRHS
!C     ..
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER            IPIV( * )
      COMPLEX*16         B( LDB, * ), H( LDH, * )
#else
      INTEGER,     DIMENSION(:), ALLOCATABLE :: IPIV
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: B,H
#endif
*/
void MB02RZ(const char *,
            int *,
	    int *,
	    std::complex<double> __restrict *,
	    int *,
	    int __restrict *,
	    std::complex<double> __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE SB02MD( DICO, HINV, UPLO, SCAL, SORT, N, A, LDA, G, &
    LDG, Q, LDQ, RCOND, WR, WI, S, LDS, U, LDU,                &
    IWORK, DWORK, LDWORK, BWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         DICO, HINV, SCAL, SORT, UPLO
      INTEGER           INFO, LDA, LDG, LDQ, LDS, LDU, LDWORK, N
      DOUBLE PRECISION  RCOND
      !C     .. Array Arguments ..

      LOGICAL           BWORK(*)
      INTEGER           IWORK(*)
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), DWORK(*), G(LDG,*), Q(LDQ,*), &
           S(LDS,*), U(LDU,*), WR(*), WI(*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: G
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Q
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: S
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WR
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: WI
#endif

*/
void SB02MD(const char *,
            const char *,
	    const char *,
	    const char *,
	    const char *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    int *);

/*
SUBROUTINE SB02MU(DICO, HINV, UPLO, N, A, LDA, G, LDG, Q, LDQ, S, &
      LDS, IWORK, DWORK, LDWORK, INFO)
!C     .. Scalar Arguments ..
      CHARACTER         DICO, HINV, UPLO
      INTEGER           INFO, LDA, LDG, LDQ, LDS, LDWORK, N
      !C     .. Array Arguments ..
 
      INTEGER           IWORK(*)
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1      
      DOUBLE PRECISION  A(LDA,*), DWORK(*), G(LDG,*), Q(LDQ,*), &
           S(LDS,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: G
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Q
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: S
#endif
*/
void SB02MU(const char *,
            const char *,
	    const char *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE SB02ND( DICO, FACT, UPLO, JOBL, N, M, P, A, LDA, B, &
    LDB, R, LDR, IPIV, L, LDL, X, LDX, RNORM, F,               &
    LDF, OUFACT, IWORK, DWORK, LDWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         DICO, FACT, JOBL, UPLO
      INTEGER           INFO, LDA, LDB, LDF, LDL, LDR, LDWORK, LDX, M, &
                        N, P
      DOUBLE PRECISION  RNORM
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1 
      INTEGER           IPIV(*), IWORK(*), OUFACT(2)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), DWORK(*), F(LDF,*), &
           L(LDL,*), R(LDR,*), X(LDX,*)
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      INTEGER, DIMENSION(2) :: OUFACT
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: F
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: L
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X
#endif      
*/
void SB02ND(const char *,
            const char *,
	    const char *,
	    const char *,
	    int *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);


/*
SUBROUTINE SB02MT( JOBG, JOBL, FACT, UPLO, N, M, A, LDA, B, LDB, &
Q, LDQ, R, LDR, L, LDL, IPIV, OUFACT, G, LDG,    &
IWORK, DWORK, LDWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         FACT, JOBG, JOBL, UPLO
      INTEGER           INFO, LDA, LDB, LDG, LDL, LDQ, LDR, LDWORK, M, &
                        N, OUFACT
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      !INTEGER           IPIV(*), IWORK(*)
      !DOUBLE PRECISION  A(LDA,*), B(LDB,*), DWORK(*), G(LDG,*), &
      !     L(LDL,*), Q(LDQ,*), R(LDR,*)
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: G
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: L
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Q
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
#endif
*/
void SB02MT(const char *,
            const char *,
	    const char *,
	    const char *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE MB01RB( SIDE, UPLO, TRANS, M, N, ALPHA, BETA, R, LDR, &
      A, LDA, B, LDB, INFO)
!C     .. Scalar Arguments ..
      CHARACTER         SIDE, TRANS, UPLO
      INTEGER           INFO, LDA, LDB, LDR, M, N
      DOUBLE PRECISION  ALPHA, BETA
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), R(LDR,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
#endif      
*/
void MB01RB(const char *,
            const char *,
	    const char *,
	    int *,
	    int *,
	    double *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE MB01RX( SIDE, UPLO, TRANS, M, N, ALPHA, BETA, R, LDR, &
      A, LDA, B, LDB, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         SIDE, TRANS, UPLO
      INTEGER           INFO, LDA, LDB, LDR, M, N
      DOUBLE PRECISION  ALPHA, BETA
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1  
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), R(LDR,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
#endif   
*/
void MB01RX(const char *,
            const char *,
	    const char *,
	    int *,
	    int *,
	    double *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE SB10ZD( N, M, NP, A, LDA, B, LDB, C, LDC, D, LDD, &
  FACTOR, AK, LDAK, BK, LDBK, CK, LDCK, DK, &
LDDK, RCOND, TOL, IWORK, DWORK, LDWORK, BWORK, &
INFO )
   !C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDAK, LDB, LDBK, LDC, LDCK, LDD, &
                         LDDK, LDWORK, M, N, NP
      DOUBLE PRECISION   FACTOR, TOL
!C     ..
!C     .. Array Arguments ..
       
      INTEGER            IWORK( * )
      LOGICAL            BWORK( * )
      DOUBLE PRECISION, DIMENSION(6) :: RCOND
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1  
      DOUBLE PRECISION   A ( LDA,  * ), AK( LDAK, * ), B ( LDB,  * ), &
                         BK( LDBK, * ), C ( LDC,  * ), CK( LDCK, * ), &
                         D ( LDD,  * ), DK( LDDK, * ), DWORK( * )
                         
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
     
#endif
*/
void SB10ZD(int *,
            int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    int *);

/*
SUBROUTINE SB02OD( DICO, JOBB, FACT, UPLO, JOBL, SORT, N, M, P, A, &
LDA, B, LDB, Q, LDQ, R, LDR, L, LDL, RCOND, X,                     &
LDX, ALFAR, ALFAI, BETA, S, LDS, T, LDT, U,                        &
LDU, TOL, IWORK, DWORK, LDWORK, BWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         DICO, FACT, JOBB, JOBL, SORT, UPLO
      INTEGER           INFO, LDA, LDB, LDL, LDQ, LDR, LDS, LDT, LDU, &
                        LDWORK, LDX, M, N, P
      DOUBLE PRECISION  RCOND, TOL
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1 
      LOGICAL           BWORK(*)
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), ALFAI(*), ALFAR(*), B(LDB,*), BETA(*), &
                        DWORK(*), L(LDL,*), Q(LDQ,*), R(LDR,*),          &
                        S(LDS,*), T(LDT,*), U(LDU,*), X(LDX,*)
#else
      LOGICAL, DIMENSION(:), ALLOCATABLE :: BWORK
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: ALFAI
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: ALFAR
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: BETA
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: L
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Q
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: S
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: T
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X
#endif     
*/
void SB02OD(const char *,
            const char *,
	    const char *,
	    const char *,
	    const char *,
	    const char *,
	    int *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    int *);

/*
SUBROUTINE SB02OY( TYPE, DICO, JOBB, FACT, UPLO, JOBL, JOBE, N, M, &
  P, A, LDA, B, LDB, Q, LDQ, R, LDR, L, LDL, E,                    &
  LDE, AF, LDAF, BF, LDBF, TOL, IWORK, DWORK,                      &
  LDWORK, INFO)
!C     .. Scalar Arguments ..
      CHARACTER         DICO, FACT, JOBB, JOBE, JOBL, TYPE, UPLO
      INTEGER           INFO, LDA, LDAF, LDB, LDBF, LDE, LDL, LDQ, LDR, &
                        LDWORK, M, N, P
      DOUBLE PRECISION  TOL
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), AF(LDAF,*), B(LDB,*), BF(LDBF,*),  &
           DWORK(*), E(LDE,*), L(LDL,*), Q(LDQ,*), R(LDR,*)
#else
      INTEGER, DIMENSION(:),ALLOCATABLE :: IWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AF
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BF
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: E
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: L
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Q
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
#endif
*/
void SB02OY(const char *,
            const char *,
	    const char *,
	    const char *,
	    const char *,
	    const char *,
	    const char *,
	    int *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE SB03TD( JOB, FACT, TRANA, UPLO, LYAPUN, N, SCALE, A,
LDA, T, LDT, U, LDU, C, LDC, X, LDX, SEP,
RCOND, FERR, WR, WI, IWORK, DWORK, LDWORK,INFO)
!C     .. Scalar Arguments ..
      CHARACTER          FACT, JOB, LYAPUN, TRANA, UPLO
      INTEGER            INFO, LDA, LDC, LDT, LDU, LDWORK, LDX, N
      DOUBLE PRECISION   FERR, RCOND, SCALE, SEP
!C     ..
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), DWORK( * ), &
                         T( LDT, * ), U( LDU, * ), WI( * ), WR( * ), &
                         X( LDX, * )
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: T
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WI
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WR
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X
#endif
*/
void SB03TD(const char *,
            const char *,
	    const char *,
	    const char *,
	    const char *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double *,
	    double *,
	    double __restrict *,
	    double __restrict *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE MB01RU( UPLO, TRANS, M, N, ALPHA, BETA, R, LDR, A, LDA,
  X, LDX, DWORK, LDWORK, INFO)
!C     .. Scalar Arguments ..
      CHARACTER         TRANS, UPLO
      INTEGER           INFO, LDA, LDR, LDWORK, LDX, M, N
      DOUBLE PRECISION  ALPHA, BETA
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      !DOUBLE PRECISION  A(LDA,*), DWORK(*), R(LDR,*), X(LDX,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X
#endif
*/
void MB01RU(const char *,
            const char *,
	    int *,
	    int *,
	    double *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE SB03MY( TRANA, N, A, LDA, C, LDC, SCALE, INFO)
!C     .. Scalar Arguments ..
      CHARACTER          TRANA
      INTEGER            INFO, LDA, LDC, N
      DOUBLE PRECISION   SCALE
!C     ..
!C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * )
!C     ..
*/
void SB03MY(const char *,
            int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    int *);

/*
SUBROUTINE SB03MW( LTRAN, LUPPER, T, LDT, B, LDB, SCALE, X, LDX,
  XNORM, INFO )
!C     .. Scalar Arguments ..
      LOGICAL            LTRAN, LUPPER
      INTEGER            INFO, LDB, LDT, LDX
      DOUBLE PRECISION   SCALE, XNORM
!C     ..
      !C     .. Array Arguments ..
  
      DOUBLE PRECISION   B( LDB, * ), T( LDT, * ), X( LDX, * )
*/
void SB03MW(int *,
            int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double *,
	    int *);

/*
SUBROUTINE SB03QD( JOB, FACT, TRANA, UPLO, LYAPUN, N, SCALE, A, &
LDA, T, LDT, U, LDU, C, LDC, X, LDX, SEP,  &
RCOND, FERR, IWORK, DWORK, LDWORK, INFO)
!C     .. Scalar Arguments ..
      CHARACTER          FACT, JOB, LYAPUN, TRANA, UPLO
      INTEGER            INFO, LDA, LDC, LDT, LDU, LDWORK, LDX, N
      DOUBLE PRECISION   FERR, RCOND, SCALE, SEP
!C     ..
      !C!     .. Array Arguments ..
  
      INTEGER            IWORK( * )
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), DWORK( * ), &
           T( LDT, * ), U( LDU, * ), X( LDX, * )
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: T
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X
#endif
*/
void SB03QD(const char *,
            const char *,
	    const char *,
	    const char *,
	    const char *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE MB01UW( SIDE, TRANS, M, N, ALPHA, H, LDH, A, LDA, &
     DWORK, LDWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         SIDE, TRANS
      INTEGER           INFO, LDA, LDH, LDWORK, M, N
      DOUBLE PRECISION  ALPHA
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), DWORK(*), H(LDH,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: H
#endif
*/
void MB01UW(const char *,
            const char *,
	    int *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int *);

/*
 SUBROUTINE SB03QX( TRANA, UPLO, LYAPUN, N, XANORM, T, LDT, U, LDU, &
       R, LDR, FERR, IWORK, DWORK, LDWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER          LYAPUN, TRANA, UPLO
      INTEGER            INFO, LDR, LDT, LDU, LDWORK, N
      DOUBLE PRECISION   FERR, XANORM
!C     ..
      !C     .. Array Arguments ..

      INTEGER            IWORK( * )
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1      
     DOUBLE PRECISION   DWORK( * ), R( LDR, * ), T( LDT, * ), &
     U( LDU, * )
#else 
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: T
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
#endif
*/
void SB03QX(const char *,
            const char *,
	    const char *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE SB03QY( JOB, TRANA, LYAPUN, N, T, LDT, U, LDU, X, LDX, &
     SEP, THNORM, IWORK, DWORK, LDWORK, INFO)
!C     .. Scalar Arguments ..
      CHARACTER          JOB, LYAPUN, TRANA
      INTEGER            INFO, LDT, LDU, LDWORK, LDX, N
      DOUBLE PRECISION   SEP, THNORM
!C     ..
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER            IWORK( * )
      DOUBLE PRECISION   DWORK( * ), T( LDT, * ), U( LDU, * ), &
           X( LDX, * )
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: T
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X
#endif
*/
void SB03QY(const char *,
            const char *,
	    const char *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE FD01AD( JP, L, LAMBDA, XIN, YIN, EFOR, XF, EPSBCK, &
CTETA, STETA, YQ, EPOS, EOUT, SALPH, IWARN, &
INFO )
!C     .. Scalar Arguments ..
      CHARACTER         JP
      INTEGER           INFO, IWARN, L
      DOUBLE PRECISION  EFOR, EOUT, EPOS, LAMBDA, XIN, YIN
!C     .. Array Arguments ..
      !DOUBLE PRECISION  CTETA(*), EPSBCK(*), SALPH(*), STETA(*), XF(*), &
      !                  YQ(*)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CTETA
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: EPSBCK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: SALPH
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: STETA
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: XF
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: YQ

*/
void FD01AD(const char *,
            int *,
	    double *,
	    double *,
	    double *,
	    double *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double *,
	    double *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE TD05AD( UNITF, OUTPUT, NP1, MP1, W, A, B, VALR, VALI, &
     INFO )
!C     .. Scalar Arguments ..
      CHARACTER         OUTPUT, UNITF
      INTEGER           INFO, MP1, NP1
      DOUBLE PRECISION  VALI, VALR, W
!C     .. Array Arguments ..
      !DOUBLE PRECISION  A(*), B(*)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: B
*/
void TD05AD(const char *,
            const char *,
	    int *,
	    int *,
	    double *,
	    double __restrict *,
	    double __restrict *,
	    double *,
	    double *,
	    int *);

/*
SUBROUTINE TF01QD( NC, NB, N, IORD, AR, MA, H, LDH, INFO )
!C     .. Scalar Arguments ..
      INTEGER           INFO, LDH, N, NB, NC
!C     .. Array Arguments ..
      INTEGER           IORD(*)
      DOUBLE PRECISION  AR(*), H(LDH,*), MA(*)
*/
void TF01QD(int *,
            int *,
	    int *,
	    int __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE TF01RD( NA, NB, NC, N, A, LDA, B, LDB, C, LDC, H, LDH, &
     DWORK, LDWORK, INFO )
!C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDB, LDC, LDH, LDWORK, N, NA, NB, NC
!C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), H(LDH,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: H
#endif  
*/
void TF01RD(int *,
            int *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE TC05AD( LERI, M, P, SVAL, INDEX, PCOEFF, LDPCO1, &
LDPCO2, QCOEFF, LDQCO1, LDQCO2, RCOND, CFREQR,&
LDCFRE, IWORK, DWORK, ZWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         LERI
      INTEGER           INFO, LDCFRE, LDPCO1, LDPCO2, LDQCO1, LDQCO2, M, P
    
      DOUBLE PRECISION  RCOND
      COMPLEX*16        SVAL
!C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER           INDEX(*), IWORK(*)
      DOUBLE PRECISION  DWORK(*), PCOEFF(LDPCO1,LDPCO2,*), &
                        QCOEFF(LDQCO1,LDQCO2,*)
      COMPLEX*16        CFREQR(LDCFRE,*), ZWORK(*)
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: INDEX, IWORK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: PCOEFF
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: QCOEFF
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: CFREQR
      COMPLEX(16), DIMENSION(:), ALLOCATABLE :: ZWORK
#endif
*/
void TC05AD(const char *,
            int *,
	    int *,
	    std::complex<double> *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    int *,
	    double *,
	    std::complex<double> __restrict *,
	    int *,
	    int __restrict *,
	    double __restrict *,
	    std::complex<double> __restrict *,
	    int *);

/*
SUBROUTINE SB10KD( N, M, NP, A, LDA, B, LDB, C, LDC, FACTOR, &
AK, LDAK, BK, LDBK, CK, LDCK, DK, LDDK, RCOND, &
IWORK, DWORK, LDWORK, BWORK, INFO )
!C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDAK, LDB, LDBK, LDC, LDCK, LDDK, &
                         LDWORK, M, N, NP
      DOUBLE PRECISION   FACTOR
!C     ..
!C     .. Array Arguments ..
       DOUBLE PRECISION, DIMENSION(4) :: RCOND
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
     INTEGER            IWORK( * )
     LOGICAL            BWORK( * )
      DOUBLE PRECISION   A( LDA, * ), AK( LDAK, * ), B( LDB, * ), &
                       BK( LDBK, * ), C( LDC, * ), CK( LDCK, * ), &
                         DK( LDDK, * ), DWORK( * )
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      LOGICAL, DIMENSION(:), ALLOCATABLE :: BWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
#endif     
*/
void SB10KD(int *,
            int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    int *);

/*
SUBROUTINE SB10AD( JOB, N, M, NP, NCON, NMEAS, GAMMA, A, LDA, &
B, LDB, C, LDC, D, LDD, AK, LDAK, BK, LDBK, CK, &
LDCK, DK, LDDK, AC, LDAC, BC, LDBC, CC, LDCC, &
DC, LDDC, RCOND, GTOL, ACTOL, IWORK, LIWORK, &
DWORK, LDWORK, BWORK, LBWORK, INFO )
!C     .. Scalar Arguments ..
      INTEGER            INFO, JOB, LBWORK, LDA, LDAC, LDAK, LDB, LDBC, &
                         LDBK, LDC, LDCC, LDCK, LDD, LDDC, LDDK, LDWORK, &
                         LIWORK, M, N, NCON, NMEAS, NP
      DOUBLE PRECISION   ACTOL, GAMMA, GTOL
!C     ..
!C     .. Array Arguments ..
      LOGICAL            BWORK( * )
      INTEGER            IWORK( * )
      DOUBLE PRECISION, DIMENSION(4) :: RCOND
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION   A( LDA, * ), AC( LDAC, * ), AK( LDAK, * ), &
                        B( LDB, * ), BC( LDBC, * ), BK( LDBK, * ), &
                        C( LDC, * ), CC( LDCC, * ), CK( LDCK, * ), &
                        D( LDD, * ), DC( LDDC, * ), DK( LDDK, * ), &
                        DWORK( * )
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DK
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
#endif
*/
void SB10AD(int *,
            int *,
	    int *,
	    int *,
	    int *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
            double __restrict *,
	    double *,
	    double *,
	    int __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int __restrict *,
            int *,
	    int *);
	    
/*
SUBROUTINE SB10LD( N, M, NP, NCON, NMEAS, A, LDA, B, LDB, C, LDC, &
D, LDD, AK, LDAK, BK, LDBK, CK, LDCK, DK, LDDK, &
AC, LDAC, BC, LDBC, CC, LDCC, DC, LDDC, IWORK, &
DWORK, LDWORK, INFO )
!C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDAC, LDAK, LDB, LDBC, LDBK, LDC, &
                        LDCC, LDCK, LDD, LDDC, LDDK, LDWORK, M, N, &
                        NCON, NMEAS, NP
!C     ..
!C     .. Array Arguments ..
      INTEGER            IWORK( * )
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION   A( LDA, * ), AC( LDAC, * ), AK( LDAK, * ),
                      B( LDB, * ), BC( LDBC, * ), BK( LDBK, * ),
                      C( LDC, * ), CC( LDCC, * ), CK( LDCK, * ),
                      D( LDD, * ), DC( LDDC, * ), DK( LDDK, * ),
                      DWORK( * )
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DK
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
#endif
*/
void SB10LD(int *,
            int *,
	    int *,
	    int *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
            int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE SB10PD( N, M, NP, NCON, NMEAS, A, LDA, B, LDB, C, LDC, &
D, LDD, TU, LDTU, TY, LDTY, RCOND, TOL, DWORK, &
LDWORK, INFO )
!C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LDC, LDD, LDTU, LDTY, LDWORK, &
                         M, N, NCON, NMEAS, NP
      DOUBLE PRECISION   TOL
!C     ..
!C     .. Array Arguments ..
       DOUBLE PRECISION, DIMENSION(2)  :: RCOND
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ), &
                       D( LDD, * ), DWORK( * ), &
                        TU( LDTU, * ), TY( LDTY, * )
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: TU
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: TY
#endif      
*/
void SB10PD(int *,
            int *,
	    int *,
	    int *,
	    int *,
	    int *,
            double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    double *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE SB10QD( N, M, NP, NCON, NMEAS, GAMMA, A, LDA, B, LDB,
C, LDC, D, LDD, F, LDF, H, LDH, X, LDX, Y, LDY,
XYCOND, IWORK, DWORK, LDWORK, BWORK, INFO )
!C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LDC, LDD, LDF, LDH, LDWORK, &
                         LDX, LDY, M, N, NCON, NMEAS, NP
      DOUBLE PRECISION   GAMMA
!C     ..
!C     .. Array Arguments ..
      INTEGER            IWORK( * )
      LOGICAL            BWORK( * )
      DOUBLE PRECISION,DIMENSION(2) :: XYCOND
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
     DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ), &
                        D( LDD, * ), DWORK( * ),  F( LDF, * ), &
                        H( LDH, * ), X( LDX, * ), &
                        Y( LDY, * )
#else     
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: D
      DOUBLE PRECISION,DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: F
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: H
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: X
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: Y
#endif
*/
void SB10QD((int *,
            int *,
	    int *,
	    int *,
	    int *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
            double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    int *);

/*

SUBROUTINE SB02RD( JOB, DICO, HINV, TRANA, UPLO, SCAL, SORT, FACT, &
LYAPUN, N, A, LDA, T, LDT, V, LDV, G, LDG, Q, &
LDQ, X, LDX, SEP, RCOND, FERR, WR, WI, S, LDS, &
IWORK, DWORK, LDWORK, BWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         DICO, FACT, HINV, JOB, LYAPUN, SCAL, SORT, &
                        TRANA, UPLO
      INTEGER           INFO, LDA, LDG, LDQ, LDS, LDT, LDV, LDWORK, LDX, &
                       N
      DOUBLE PRECISION  FERR, RCOND, SEP
!C     .. Array Arguments ..
      LOGICAL           BWORK(*)
      INTEGER           IWORK(*)
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), DWORK(*), G(LDG,*), Q(LDQ,*),
                      S(LDS,*), T(LDT,*), V(LDV,*), WI(*), WR(*),
                       X(LDX,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: G
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Q
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: S
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: T
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: V
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WI
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WR
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X
#endif   
*/
void SB02RD(const char *,
            const char *,
	    const char *,
	    const char *,
	    const char *,
	    const char *,
	    const char *,
	    const char *,
	    const char *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double *,
	    double *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    int *);

/*
SUBROUTINE MB02PD( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, &
EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, &
IWORK, DWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         EQUED, FACT, TRANS
      INTEGER           INFO, LDA, LDAF, LDB, LDX, N, NRHS
      DOUBLE PRECISION  RCOND
!C     ..
!C     .. Array Arguments ..
      INTEGER           IPIV( * ), IWORK( * )
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A( LDA, * ), AF( LDAF, * ), B( LDB, * ), &
                      BERR( * ), C( * ), DWORK( * ), FERR( * ), &
                        R( * ), X( LDX, * )
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AF
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: BERR
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FERR
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: R
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X
#endif
*/
void MB02PD(const char *,
            const char *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    const char *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double __restrict *,
	    double __restrict *,
	    int __restrict *,
	    double __restrict *,
	    int *);

/*
*/



}

#endif

#endif /*__GMS_SLICOT_C_IFACE_H__*/
