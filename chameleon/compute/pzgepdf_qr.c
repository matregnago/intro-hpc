/**
 *
 * @file pzgepdf_qr.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2016-2020 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgepdf_qr parallel algorithm
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

void
chameleon_pzgepdf_qr( int genD, int doqr, int optid,
                      const libhqr_tree_t *qrtreeT, const libhqr_tree_t *qrtreeB,
                      CHAM_desc_t *A1, CHAM_desc_t *TS1, CHAM_desc_t *TT1, CHAM_desc_t *D1, CHAM_desc_t *Q1,
                      CHAM_desc_t *A2, CHAM_desc_t *TS2, CHAM_desc_t *TT2, CHAM_desc_t *D2, CHAM_desc_t *Q2,
                      RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    int L, KT, Dtop;

    if ( A1->m > A1->n ) {
        KT = A1->nt;
    } else {
        KT = A1->mt;
    }
    L = optid ? A2->m : 0;

    Dtop = genD;

    /* Factorize A1 if needed */
    if ( doqr ) {
        if ( qrtreeT ) {
            chameleon_pzgeqrf_param( genD, KT, qrtreeT, A1,
                                     TS1, TT1, D1, sequence, request );
        }
        else {
            chameleon_pzgeqrf( genD, A1, TS1, D1, sequence, request );
        }
        Dtop = 0;
    }

    /* Factorize A2 */
    if ( qrtreeT ) {
        chameleon_pztpqrt_param( genD, (optid ? ChamUpper : ChamUpperLower), KT, qrtreeB,
                                 A1, A2, TS2, TT2, D2, sequence, request );
    }
    else {
        chameleon_pztpqrt( L, A1, A2, TS2, sequence, request );
    }

    /* Initialize Q1 and Q2 */
    chameleon_pzlaset( ChamUpperLower, 0.0, 1.0, Q1, sequence, request );
    chameleon_pzlaset( ChamUpperLower, 0.0, 0.0, Q2, sequence, request );

    /* Generate Q1 and Q2 */
    if ( qrtreeT ) {
        chameleon_pztpgqrt_param( 0, (optid ? ChamUpper : ChamUpperLower), KT, qrtreeB,
                                  A2, Q1, Q2, TS2, TT2, D2, sequence, request );
        chameleon_pzungqr_param( Dtop, KT, qrtreeT, A1, Q1, TS1, TT1, D1, sequence, request );
    }
    else {
        chameleon_pztpgqrt( KT, L, A2, TS2, Q1, Q2, sequence, request );
        chameleon_pzungqr( Dtop, A1, Q1, TS1, D1, sequence, request );
    }

    /* Flush all read data */
    CHAMELEON_Desc_Flush( A1,  sequence );
    CHAMELEON_Desc_Flush( TS1, sequence );

    CHAMELEON_Desc_Flush( A2,  sequence );
    CHAMELEON_Desc_Flush( TS2, sequence );

    if ( qrtreeT ) {
        CHAMELEON_Desc_Flush( TT1, sequence );
        CHAMELEON_Desc_Flush( TT2, sequence );
    }

    if ( D1 != NULL ) {
        CHAMELEON_Desc_Flush( D1, sequence );
    }
    if ( D2 != NULL ) {
        CHAMELEON_Desc_Flush( D2, sequence );
    }

    return;
}
