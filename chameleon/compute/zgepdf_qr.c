/**
 *
 * @file zgepdf_qr.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2016-2020 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgepdf_qr wrapper.
 *
 * This wrapper exists only for testing purpose to ease the creation of
 * timer/testing fot this subroutine of the QDWH/Zolo algorithms. That is why
 * only the Tile version is available.
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 ******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 * @brief Perfom a QR factorization of a matrix A and generate the associated Q.
 *
 * A is of the form [ A1 ] and Q is generated in a similar form [ Q1 ].
 *                  [ A2 ]                                      [ Q2 ]
 *
 * @warning This function is NOT a USER routine and should not be called directly.
 *
 *******************************************************************************
 *
 * @param[in] doqr
 *            Specify if A1 is already factorized or if it needs to be done within this call.
 *
 * @param[in] optid
 *            Specify if A2 is an identity structure matrix and optimization is performed or not.
 *
 * @param[in] qrtreeT
 *            Describe the reduction tree for the factorization of A1.
 *
 * @param[in] qrtreeB
 *            Describe the reduction tree for the factorization of A2.
 *
 * @param[in] A1
 *
 * @param[in] TS1
 *
 * @param[in] TT1
 *
 * @param[in] Q1
 *
 * @param[in] A2
 *
 * @param[in] TS2
 *
 * @param[in] TT2
 *
 * @param[in] Q2
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_cgepdf_qr_Tile
 * @sa CHAMELEON_dgepdf_qr_Tile
 * @sa CHAMELEON_sgepdf_qr_Tile
 *
 */
int CHAMELEON_zgepdf_qr_Tile( int doqr, int optid,
                              const libhqr_tree_t *qrtreeT, const libhqr_tree_t *qrtreeB,
                              CHAM_desc_t *A1, CHAM_desc_t *TS1, CHAM_desc_t *TT1, CHAM_desc_t *Q1,
                              CHAM_desc_t *A2, CHAM_desc_t *TS2, CHAM_desc_t *TT2, CHAM_desc_t *Q2 )
{
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request = RUNTIME_REQUEST_INITIALIZER;
#if defined(CHAMELEON_COPY_DIAG)
    CHAM_desc_t         D1, D2;
#endif
    CHAM_desc_t        *D1ptr = NULL;
    CHAM_desc_t        *D2ptr = NULL;
    int                 status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgepdf_qr_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

#if defined(CHAMELEON_COPY_DIAG)
    {
        int n = A1->n;
        chameleon_zdesc_copy_and_restrict( chamctxt, A1, &D1, A1->m, n );
        D1ptr = &D1;
        chameleon_zdesc_copy_and_restrict( chamctxt, A2, &D2, A2->m, n );
        D2ptr = &D2;
    }
#endif

    chameleon_pzgepdf_qr( 1, doqr, optid, qrtreeT, qrtreeB,
                          A1, TS1, TT1, D1ptr, Q1,
                          A2, TS2, TT2, D2ptr, Q2,
                          sequence, &request );

    CHAMELEON_Desc_Flush( Q1, sequence );
    CHAMELEON_Desc_Flush( Q2, sequence );

    if ( D1ptr != NULL ) {
        CHAMELEON_Desc_Flush( D1ptr, sequence );
        CHAMELEON_Desc_Flush( D2ptr, sequence );
    }
    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );

    if ( D1ptr != NULL ) {
        chameleon_desc_destroy( D1ptr );
        chameleon_desc_destroy( D2ptr );
    }

    return status;
}
