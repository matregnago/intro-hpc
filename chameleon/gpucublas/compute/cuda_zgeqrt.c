/**
 *
 * @file cuda_zgeqrt.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zgeqrt GPU kernel
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2025-01-29
 * @precisions normal z -> c d s
 *
 */
#include "gpucublas.h"

#if defined(CHAMELEON_USE_MAGMA)
int
CUDA_zgeqrt(
    magma_int_t m, magma_int_t n, magma_int_t nb,
    magmaDoubleComplex *da, magma_int_t ldda,
    magmaDoubleComplex *v,  magma_int_t ldv,
    magmaDoubleComplex *dt, magma_int_t lddt,
    magmaDoubleComplex *t,  magma_int_t ldt,
    magmaDoubleComplex *dd,
    magmaDoubleComplex *d,  magma_int_t ldd,
    magmaDoubleComplex *tau,
    magmaDoubleComplex *hwork,
    magmaDoubleComplex *dwork,
    CUstream stream )
{
#define da_ref(a_1,a_2) ( da+(size_t)(a_2)*(ldda) + (a_1))
#define v_ref( a_1,a_2) ( v +(size_t)(a_2)*(ldv)  + (a_1))
#define dt_ref(a_1,a_2) ( dt+(size_t)(a_2)*(lddt) + (a_1))
#define t_ref( a_1,a_2) ( t +(size_t)(a_2)*(ldt)  + (a_1))

    int i, k, ib, old_i, old_ib, rows, cols;
    double _Complex one=1.;

    if (m < 0) {
        return -1;
    } else if (n < 0) {
        return -2;
    } else if (ldda < chameleon_max(1,m)) {
        return -4;
    }

    k = chameleon_min(m,n);
    if (k == 0) {
        hwork[0] = *((magmaDoubleComplex*) &one);
        return MAGMA_SUCCESS;
    }

    /* lower parts of little T must be zero: memset to 0 for simplicity */
    memset( t_ref(0,0), 0, sizeof(magmaDoubleComplex) * nb * n );
    cudaMemsetAsync( dt_ref(0,0), 0, sizeof(magmaDoubleComplex) * nb * n, stream );

    if ( (nb > 1) && (nb < k) ) {
        /* Use blocked code initially */
        old_i = 0; old_ib = nb;
        for (i = 0; i < k-nb; i += nb) {

            ib = chameleon_min(k-i, nb);
            rows = m -i;
            magma_zgetmatrix_async( rows, ib,
                                    da_ref(i,i), ldda,
                                    v_ref(i,0), ldv, stream );

            if (i>0){
                /* Apply H' to A(i:m,i+2*ib:n) from the left */
                cols = n-old_i-2*old_ib;
                magma_zlarfb_gpu( MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
                                  m-old_i, cols, old_ib,
                                  da_ref(old_i, old_i), ldda, dt_ref(0,old_i), lddt,
                                  da_ref(old_i, old_i+2*old_ib), ldda,
                                  dwork, cols);

                /* store the diagonal */
                magma_zsetmatrix_async( old_ib, old_ib,
                                        d,                    old_ib,
                                        da_ref(old_i, old_i), ldda, stream );
            }

            magma_queue_sync( stream );
            /* Form the triangular factor of the block reflector on the host
             H = H(i) H(i+1) . . . H(i+ib-1) */
            CORE_zgeqrt(rows, ib, ib,
                        (double _Complex*) v_ref(i, 0), ldv,
                        (double _Complex*) t_ref(0, 0), ib,
                        (double _Complex*) tau+i,
                        (double _Complex*) hwork);

            /* Put 0s in the upper triangular part of a panel (and 1s on the
             diagonal); copy the upper triangular in d. */
            CORE_zgesplit(ChamLeft, ChamUnit, chameleon_min(rows,ib), ib,
                          (double _Complex*) v_ref(i, 0), ldv,
                          (double _Complex*) d, ib);

            /* send the custom panel to the GPU */
            magma_zsetmatrix( rows, ib,
                              v_ref(i, 0), ldv,
                              da_ref(i, i), ldda );

            if ( i + ib < n ){
                /* Send the triangular factor T to the GPU */
                magma_zsetmatrix( ib, ib,
                                  t_ref(0, 0), ib,
                                  dt_ref(0, i), lddt );

                if (i+nb < k-nb) {
                    /* Apply H' to A(i:m,i+ib:i+2*ib) from the left */
                    magma_zlarfb_gpu( MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
                                      rows, ib, ib,
                                      da_ref(i,i),    ldda, dt_ref(0,i), lddt,
                                      da_ref(i,i+ib), ldda, dwork, ib);
                }
                else {
                    cols = n-i-ib;
                    magma_zlarfb_gpu( MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
                                      rows, cols, ib,
                                      da_ref(i,i),    ldda, dt_ref(0,i), lddt,
                                      da_ref(i,i+ib), ldda, dwork, cols);
                    cudaThreadSynchronize();
                    /* Fix the diagonal block */
                    magma_zsetmatrix_async( ib, ib,
                                            d,            ib,
                                            da_ref(i, i), ldda,
                                            stream );
                }
                old_i  = i;
                old_ib = ib;
            }
        }
    } else {
        i = 0;
    }

    /* Use unblocked code to factor the last or only block. */
    if (i < k) {
        ib   = n-i;
        rows = m-i;
        magma_zgetmatrix( rows, ib,
                          da_ref(i,i), ldda,
                          v_ref(i,0), ldv );
        CORE_zgeqrt(rows, ib, ib,
                    (double _Complex*) v_ref(i, 0), ldv,
                    (double _Complex*) t_ref(0, 0), ib,
                    (double _Complex*) tau+i,
                    (double _Complex*) hwork);
        /* send the last factorized panel to the GPU */
        magma_zsetmatrix( rows, ib,
                          v_ref(i, 0), ldv,
                          da_ref(i, i), ldda );
        /* Send the triangular factor T to the GPU */
        magma_zsetmatrix( ib, ib,
                          t_ref(0, 0), ib,
                          dt_ref(0, i), lddt );
    }

#undef da_ref
#undef v_ref
#undef dt_ref
#undef t_ref

    return CHAMELEON_SUCCESS;
}
#endif
