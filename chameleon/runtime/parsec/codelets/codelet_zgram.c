/**
 *
 * @file parsec/codelet_zgram.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgram PaRSEC codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zgegram_parsec( parsec_execution_stream_t *context,
                     parsec_task_t             *this_task )
{
    cham_uplo_t uplo;
    int m, n, mt, nt;
    double *Di;
    int lddi;
    double *Dj;
    int lddj;
    double *D;
    double *A;
    int lda;

    parsec_dtd_unpack_args(
        this_task, &uplo, &m, &n, &mt, &nt, &Di, &lddi, &Dj, &lddj, &D, &A, &lda );

    CORE_zgram( uplo,
                m, n, mt, nt,
                Di, lddi,
                Dj, lddj,
                D,
                A, lda);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

static inline int
CORE_zsygram_parsec( parsec_execution_stream_t *context,
                     parsec_task_t             *this_task )
{
    cham_uplo_t uplo;
    int m, n, mt, nt;
    double *Di;
    int lddi;
    double *D;
    double *A;
    int lda;

    parsec_dtd_unpack_args(
        this_task, &uplo, &m, &n, &mt, &nt, &Di, &lddi, &D, &A, &lda );

    CORE_zgram( uplo,
                m, n, mt, nt,
                Di, lddi,
                Di, lddi,
                D,
                A, lda);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zgram( const RUNTIME_option_t *options,
                        cham_uplo_t uplo,
                        int m, int n, int mt, int nt,
                        const CHAM_desc_t *Di, int Dim, int Din,
                        const CHAM_desc_t *Dj, int Djm, int Djn,
                        const CHAM_desc_t *D, int Dm, int Dn,
                        CHAM_desc_t *A, int Am, int An)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileDi = Di->get_blktile( Di, Dim, Din );
    CHAM_tile_t *tileDj = Dj->get_blktile( Dj, Djm, Djn );
    CHAM_tile_t *tileA  = A->get_blktile( A, Am, An );
    double *ptrDi, *ptrDj;

    /*
     * Test if Di is Dj, when we are on the diagonal.
     * This to avoid having the same data twice in inputs (not handled in parsec).
     */
    ptrDi = (double *)(RTBLKADDR( Di, double, Dim, Din ));
    ptrDj = (double *)(RTBLKADDR( Dj, double, Djm, Djn ));
    if (ptrDi == ptrDj) {
        parsec_dtd_taskpool_insert_task(
            PARSEC_dtd_taskpool, CORE_zsygram_parsec, options->priority, "sygram",
            sizeof(cham_uplo_t),   &uplo, VALUE,
            sizeof(int),   &m,    VALUE,
            sizeof(int),   &n,    VALUE,
            sizeof(int),   &mt,   VALUE,
            sizeof(int),   &nt,   VALUE,
            PASSED_BY_REF, RTBLKADDR( Di, double, Dim, Din ), chameleon_parsec_get_arena_index( Di ) | INPUT,
            sizeof(int), &(tileDi->ld), VALUE,
            PASSED_BY_REF, RTBLKADDR( D, double, Dm, Dn ), chameleon_parsec_get_arena_index( D ) | INPUT,
            PASSED_BY_REF, RTBLKADDR( A, double, Am, An ), chameleon_parsec_get_arena_index( A ) | INOUT | AFFINITY,
            sizeof(int), &(tileA->ld), VALUE,
            PARSEC_DTD_ARG_END );
    } else {
        parsec_dtd_taskpool_insert_task(
            PARSEC_dtd_taskpool, CORE_zgegram_parsec, options->priority, "gegram",
            sizeof(cham_uplo_t),   &uplo, VALUE,
            sizeof(int),   &m,    VALUE,
            sizeof(int),   &n,    VALUE,
            sizeof(int),   &mt,   VALUE,
            sizeof(int),   &nt,   VALUE,
            PASSED_BY_REF, RTBLKADDR( Di, double, Dim, Din ), chameleon_parsec_get_arena_index( Di ) | INPUT,
            sizeof(int), &(tileDi->ld), VALUE,
            PASSED_BY_REF, RTBLKADDR( Dj, double, Djm, Djn ), chameleon_parsec_get_arena_index( Dj ) | INPUT,
            sizeof(int), &(tileDj->ld), VALUE,
            PASSED_BY_REF, RTBLKADDR( D, double, Dm, Dn ), chameleon_parsec_get_arena_index( D ) | INPUT,
            PASSED_BY_REF, RTBLKADDR( A, double, Am, An ), chameleon_parsec_get_arena_index( A ) | INOUT | AFFINITY,
            sizeof(int), &(tileA->ld), VALUE,
            PARSEC_DTD_ARG_END );
    }

}
