/**
 *
 * @file zhemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhemm wrappers
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Lionel Eyraud-Dubois
 * @author Pierre Esterie
 * @date 2025-12-19
 * @precisions normal z -> c
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 * @brief Allocate the required workspaces for asynchronous hemm
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether the hermitian matrix A appears on the
 *          left or right in the operation as follows:
 *          = ChamLeft:      \f[ C = \alpha \times A \times B + \beta \times C \f]
 *          = ChamRight:     \f[ C = \alpha \times B \times A + \beta \times C \f]
 *
 * @param[in] uplo
 *          Specifies whether the upper or lower triangular part of
 *          the hermitian matrix A is to be referenced as follows:
 *          = ChamLower: Only the lower triangular part of the
 *                hermitian matrix A is to be referenced.
 *          = ChamUpper: Only the upper triangular part of the
 *                hermitian matrix A is to be referenced.
 *
 * @param[in] A
 *          The descriptor of the matrix A.
 *
 * @param[in] B
 *          The descriptor of the matrix B.
 *
 * @param[in] C
 *          The descriptor of the matrix C.
 *
 *******************************************************************************
 *
 * @retval An allocated opaque pointer to use in CHAMELEON_zhemm_Tile_Async()
 * and to free with CHAMELEON_zhemm_WS_Free().
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zhemm_Tile_Async
 * @sa CHAMELEON_zhemm_WS_Free
 *
 */
void *CHAMELEON_zhemm_WS_Alloc( cham_side_t        side __attribute__((unused)),
                                cham_uplo_t        uplo __attribute__((unused)),
                                const CHAM_desc_t *A,
                                const CHAM_desc_t *B,
                                const CHAM_desc_t *C )
{
    CHAM_context_t            *chamctxt;
    struct chameleon_pzgemm_s *options;
    int                        P, Q;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        return NULL;
    }

    options = calloc( 1, sizeof(struct chameleon_pzgemm_s) );
    options->alg = ChamGemmAlgAuto;

    /*
     * If only one process, or if generic has been globally enforced, we switch
     * to generic immediately.
     */
    P = chameleon_desc_datadist_get_iparam(C, 0);
    Q = chameleon_desc_datadist_get_iparam(C, 1);
    if ( ((P == 1) && (Q == 1)) ||
         (chamctxt->generic_enabled == CHAMELEON_TRUE) )
    {
        options->alg = ChamGemmAlgGeneric;
    }

    /* Look at environment variable is something enforces the variant. */
    if ( options->alg == ChamGemmAlgAuto )
    {
        char *algostr = chameleon_getenv( "CHAMELEON_GEMM_ALGO" );

        if ( algostr ) {
            if ( strcasecmp( algostr, "summa_c" ) == 0 ) {
                options->alg = ChamGemmAlgSummaC;
            }
            else if ( strcasecmp( algostr, "summa_a" ) == 0  ) {
                options->alg = ChamGemmAlgSummaA;
            }
            else if ( strcasecmp( algostr, "summa_b" ) == 0  ) {
                options->alg = ChamGemmAlgSummaB;
            }
            else if ( strcasecmp( algostr, "generic" ) == 0  ) {
                options->alg = ChamGemmAlgGeneric;
            }
            else if ( strcasecmp( algostr, "auto" ) == 0  ) {
                options->alg = ChamGemmAlgAuto;
            }
            else {
                fprintf( stderr, "ERROR: CHAMELEON_GEMM_ALGO is not one of AUTO, SUMMA_A, SUMMA_B, SUMMA_C, GENERIC => Switch back to Automatic switch\n" );
            }
        }
        chameleon_cleanenv( algostr );
    }

    /* Perform automatic choice if not already enforced. */
    if ( options->alg == ChamGemmAlgAuto )
    {
        double sizeA, sizeB, sizeC;
        double ratio = 1.5; /* Arbitrary ratio to give more weight to writes wrt reads. */

        /* Compute the average array per node for each matrix */
        sizeA = ((double)A->m * (double)A->n) / (double)(chameleon_desc_datadist_get_iparam(A, 0) * chameleon_desc_datadist_get_iparam(A, 1));
        sizeB = ((double)B->m * (double)B->n) / (double)(chameleon_desc_datadist_get_iparam(B, 0) * chameleon_desc_datadist_get_iparam(B, 1));
        sizeC = ((double)C->m * (double)C->n) / (double)(P * Q) * ratio;

        if ( (sizeC > sizeA) && (sizeC > sizeB) ) {
            options->alg = ChamGemmAlgSummaC;
        }
        else {
            if ( sizeA > sizeB ) {
                options->alg = ChamGemmAlgSummaA;
            }
            else {
                options->alg = ChamGemmAlgSummaB;
            }
        }
    }

    assert( options->alg != ChamGemmAlgAuto );

    /* Switch back to generic if descriptors are not 2dbc */
    if ( options->alg == ChamGemmAlgSummaC )
    {
        if ( (A->get_rankof_init != chameleon_getrankof_2d) ||
             (B->get_rankof_init != chameleon_getrankof_2d) ||
             (C->get_rankof_init != chameleon_getrankof_2d) )
        {
            options->alg = ChamGemmAlgGeneric;
        }
    }

    /* Now that we have decided which algorithm, let's allocate the required data structures. */
    if ( (options->alg == ChamGemmAlgSummaC ) &&
         (C->get_rankof_init == chameleon_getrankof_2d ) )
    {
        int lookahead = chamctxt->lookahead;

        chameleon_desc_init_2dtile( chamctxt, &(options->WA), "HEMM_WA", ChamComplexDouble,
                                    A->mb, A->nb, A->mb * C->mt, A->nb * Q * lookahead, P, Q );
        chameleon_desc_init_2dtile( chamctxt, &(options->WB), "HEMM_WB", ChamComplexDouble,
                                    B->mb, B->nb, B->mb * P * lookahead, B->nb * C->nt, P, Q );
    }

    return (void*)options;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 * @brief Free the allocated workspaces for asynchronous hemm
 *
 *******************************************************************************
 *
 * @param[in,out] user_ws
 *          On entry, the opaque pointer allocated by CHAMELEON_zhemm_WS_Alloc()
 *          On exit, all data are freed.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zhemm_Tile_Async
 * @sa CHAMELEON_zhemm_WS_Alloc
 *
 */
void CHAMELEON_zhemm_WS_Free( void *user_ws )
{
    struct chameleon_pzgemm_s *ws = (struct chameleon_pzgemm_s*)user_ws;

    if ( ws->alg == ChamGemmAlgSummaC ) {
        chameleon_desc_destroy( &(ws->WA) );
        chameleon_desc_destroy( &(ws->WB) );
    }
    free( ws );
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zhemm - Performs one of the matrix-matrix operations
 *
 *     \f[ C = \alpha \times A \times B + \beta \times C \f]
 *
 *  or
 *
 *     \f[ C = \alpha \times B \times A + \beta \times C \f]
 *
 *  where alpha and beta are scalars, A is an hermitian matrix and  B and
 *  C are m by n matrices.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether the hermitian matrix A appears on the
 *          left or right in the operation as follows:
 *          = ChamLeft:      \f[ C = \alpha \times A \times B + \beta \times C \f]
 *          = ChamRight:     \f[ C = \alpha \times B \times A + \beta \times C \f]
 *
 * @param[in] uplo
 *          Specifies whether the upper or lower triangular part of
 *          the hermitian matrix A is to be referenced as follows:
 *          = ChamLower:     Only the lower triangular part of the
 *                             hermitian matrix A is to be referenced.
 *          = ChamUpper:     Only the upper triangular part of the
 *                             hermitian matrix A is to be referenced.
 *
 * @param[in] M
 *          Specifies the number of rows of the matrix C. M >= 0.
 *
 * @param[in] N
 *          Specifies the number of columns of the matrix C. N >= 0.
 *
 * @param[in] alpha
 *          Specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is M when side = ChamLeft,
 *          and is N otherwise. Only the uplo triangular part is referenced.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,ka).
 *
 * @param[in] B
 *          B is a LDB-by-N matrix, where the leading M-by-N part of
 *          the array B must contain the matrix B.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,M).
 *
 * @param[in] beta
 *          Specifies the scalar beta.
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array is overwritten by the M by N updated matrix.
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max(1,M).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zhemm_Tile
 * @sa CHAMELEON_chemm
 * @sa CHAMELEON_dhemm
 * @sa CHAMELEON_shemm
 *
 */
int CHAMELEON_zhemm( cham_side_t side, cham_uplo_t uplo, int M, int N,
                     CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA,
                                                  CHAMELEON_Complex64_t *B, int LDB,
                     CHAMELEON_Complex64_t beta,  CHAMELEON_Complex64_t *C, int LDC )
{
    int NB;
    int Am;
    int status;
    CHAM_desc_t descAl, descAt;
    CHAM_desc_t descBl, descBt;
    CHAM_desc_t descCl, descCt;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    void *ws;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zhemm", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ( (side != ChamLeft) && (side != ChamRight) ){
        chameleon_error("CHAMELEON_zhemm", "illegal value of side");
        return -1;
    }
    if ((uplo != ChamLower) && (uplo != ChamUpper)) {
        chameleon_error("CHAMELEON_zhemm", "illegal value of uplo");
        return -2;
    }
    Am = ( side == ChamLeft ) ? M : N;
    if (M < 0) {
        chameleon_error("CHAMELEON_zhemm", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zhemm", "illegal value of N");
        return -4;
    }
    if (LDA < chameleon_max(1, Am)) {
        chameleon_error("CHAMELEON_zhemm", "illegal value of LDA");
        return -7;
    }
    if (LDB < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zhemm", "illegal value of LDB");
        return -9;
    }
    if (LDC < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zhemm", "illegal value of LDC");
        return -12;
    }

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == (CHAMELEON_Complex64_t)0.0) && beta == (CHAMELEON_Complex64_t)1.0))
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZHEMM, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zhemm", "chameleon_tune() failed");
        return status;
    }

    /* Set MT & NT & KT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInput, uplo,
                         A, NB, NB, LDA, Am, Am, Am, sequence, &request );
    chameleon_zlap2tile( chamctxt, "B", &descBl, &descBt, ChamDescInput, ChamUpperLower,
                         B, NB, NB, LDB, N, M,  N, sequence, &request );
    chameleon_zlap2tile( chamctxt, "C", &descCl, &descCt, ChamDescInout, ChamUpperLower,
                         C, NB, NB, LDC, N, M,  N, sequence, &request );

    /* Call the tile interface */
    ws = CHAMELEON_zhemm_WS_Alloc( side, uplo, &descAt, &descBt, &descCt );
    CHAMELEON_zhemm_Tile_Async( side, uplo, alpha, &descAt, &descBt, beta, &descCt, ws, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                         ChamDescInput, uplo, sequence, &request );
    chameleon_ztile2lap( chamctxt, &descBl, &descBt,
                         ChamDescInput, ChamUpperLower, sequence, &request );
    chameleon_ztile2lap( chamctxt, &descCl, &descCt,
                         ChamDescInout, ChamUpperLower, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    CHAMELEON_zhemm_WS_Free( ws );
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );
    chameleon_ztile2lap_cleanup( chamctxt, &descBl, &descBt );
    chameleon_ztile2lap_cleanup( chamctxt, &descCl, &descCt );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zhemm_Tile - Performs Hermitian matrix multiplication.
 *  Tile equivalent of CHAMELEON_zhemm().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether the hermitian matrix A appears on the
 *          left or right in the operation as follows:
 *          = ChamLeft:      \f[ C = \alpha \times A \times B + \beta \times C \f]
 *          = ChamRight:     \f[ C = \alpha \times B \times A + \beta \times C \f]
 *
 * @param[in] uplo
 *          Specifies whether the upper or lower triangular part of
 *          the hermitian matrix A is to be referenced as follows:
 *          = ChamLower:     Only the lower triangular part of the
 *                             hermitian matrix A is to be referenced.
 *          = ChamUpper:     Only the upper triangular part of the
 *                             hermitian matrix A is to be referenced.
 *
 * @param[in] alpha
 *          Specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is M when side = ChamLeft,
 *          and is N otherwise. Only the uplo triangular part is referenced.
 *
 * @param[in] B
 *          B is a LDB-by-N matrix, where the leading M-by-N part of
 *          the array B must contain the matrix B.
 *
 * @param[in] beta
 *          Specifies the scalar beta.
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array is overwritten by the M by N updated matrix.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zhemm
 * @sa CHAMELEON_zhemm_Tile_Async
 * @sa CHAMELEON_chemm_Tile
 * @sa CHAMELEON_dhemm_Tile
 * @sa CHAMELEON_shemm_Tile
 *
 */
int CHAMELEON_zhemm_Tile( cham_side_t side, cham_uplo_t uplo,
                          CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                          CHAMELEON_Complex64_t beta,  CHAM_desc_t *C )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;
    void *ws;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zhemm_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    ws = CHAMELEON_zhemm_WS_Alloc( side, uplo, A, B, C );
    CHAMELEON_zhemm_Tile_Async( side, uplo, alpha, A, B, beta, C, ws, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( B, sequence );
    CHAMELEON_Desc_Flush( C, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    CHAMELEON_zhemm_WS_Free( ws );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile_Async
 *
 *  CHAMELEON_zhemm_Tile_Async - Performs Hermitian matrix multiplication.
 *  Non-blocking equivalent of CHAMELEON_zhemm_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zhemm
 * @sa CHAMELEON_zhemm_Tile
 * @sa CHAMELEON_chemm_Tile_Async
 * @sa CHAMELEON_dhemm_Tile_Async
 * @sa CHAMELEON_shemm_Tile_Async
 *
 */
int CHAMELEON_zhemm_Tile_Async( cham_side_t side, cham_uplo_t uplo,
                                CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                                CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                                void *user_ws,
                                RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    struct chameleon_pzgemm_s *ws;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zhemm_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zhemm_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zhemm_Tile_Async", "NULL request");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == CHAMELEON_SUCCESS) {
        request->status = CHAMELEON_SUCCESS;
    }
    else {
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_SEQUENCE_FLUSHED);
    }

    /* Check descriptors for correctness */
    if (chameleon_desc_check(A) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zhemm_Tile_Async", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(B) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zhemm_Tile_Async", "invalid second descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(C) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zhemm_Tile_Async", "invalid third descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ( (side != ChamLeft) && (side != ChamRight) ){
        chameleon_error("CHAMELEON_zhemm_Tile_Async", "illegal value of side");
        return chameleon_request_fail(sequence, request, -1);
    }
    if ((uplo != ChamLower) && (uplo != ChamUpper)) {
        chameleon_error("CHAMELEON_zhemm_Tile_Async", "illegal value of uplo");
        return chameleon_request_fail(sequence, request, -2);
    }

    /* Check matrices sizes */
    if ( (B->m != C->m) || (B->n != C->n) ) {
        chameleon_error("CHAMELEON_zhemm_Tile_Async", "B and C must have the same size");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ( (A->m != A->n) ||
         ( (side == ChamLeft)  && (A->m != B->m ) ) ||
         ( (side == ChamRight) && (A->m != B->n ) ) ) {
        chameleon_error("CHAMELEON_zhemm_Tile_Async", "Matrix A must be square of size M or N regarding side.");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Check tiles sizes */
    if ( (B->mb != C->mb) || (B->nb != C->nb) ) {
        chameleon_error("CHAMELEON_zhemm_Tile_Async", "B and C must have the same tile sizes");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ( (A->mb != A->nb) ||
         ( (side == ChamLeft)  && (A->mb != B->mb ) ) ||
         ( (side == ChamRight) && (A->mb != B->nb ) ) ) {
        chameleon_error("CHAMELEON_zhemm_Tile_Async", "Matrix A must be square with square tiles wich fits the reagding tile size of B and C");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Check submatrix starting point */
    if( (A->i != 0) || (A->j != 0) ||
        (B->i != 0) || (B->j != 0) ||
        (C->i != 0) || (C->j != 0) ) {
        chameleon_error("CHAMELEON_zhemm_Tile_Async", "Submatrices are not supported for now");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if ( (C->m == 0) || (C->n == 0) ||
         ( (alpha == (CHAMELEON_Complex64_t)0.0) && (beta == (CHAMELEON_Complex64_t)1.0) ) )
    {
        return CHAMELEON_SUCCESS;
    }

    if ( user_ws == NULL ) {
        ws = CHAMELEON_zhemm_WS_Alloc( side, uplo, A, B, C );
    }
    else {
        ws = user_ws;
    }

    chameleon_pzhemm( ws, side, uplo, alpha, A, B, beta, C, sequence, request );

    if ( user_ws == NULL ) {
        CHAMELEON_Desc_Flush( A, sequence );
        CHAMELEON_Desc_Flush( B, sequence );
        CHAMELEON_Desc_Flush( C, sequence );
        chameleon_sequence_wait( chamctxt, sequence );
        CHAMELEON_zhemm_WS_Free( ws );
    }
    return CHAMELEON_SUCCESS;
}
