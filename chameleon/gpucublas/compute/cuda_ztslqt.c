/**
 *
 * @file cuda_ztslqt.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_ztslqt GPU kernel
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2025-01-29
 * @precisions normal z -> c d s
 *
 */
#include "gpucublas.h"

#if defined(CHAMELEON_USE_MAGMA) && 0
int CUDA_ztslqt(
        magma_int_t m, magma_int_t n, magma_int_t nb,
        magmaDoubleComplex *da1, magma_int_t ldda1,
        magmaDoubleComplex *da2, magma_int_t ldda2,
        magmaDoubleComplex *a2,  magma_int_t lda2,
        magmaDoubleComplex *dt,  magma_int_t lddt,
        magmaDoubleComplex *t,  magma_int_t ldt,
        magmaDoubleComplex *dd,
        magmaDoubleComplex *d,  magma_int_t ldd,
        magmaDoubleComplex *tau,
        magmaDoubleComplex *hwork,
        magmaDoubleComplex *dwork,
        CUstream stream)
{
#define da1_ref(a_1,a_2) ( da1+(size_t)(a_2)*ldda1 + (a_1))
#define da2_ref(a_1,a_2) ( da2+(size_t)(a_2)*ldda2 + (a_1))
#define a2_ref( a_1,a_2) ( a2 +(size_t)(a_2)*lda2  + (a_1))
#define t_ref(  a_1,a_2) ( t  +(size_t)(a_2)*ldt   + (a_1))
#define dt_ref( a_1,a_2) ( dt +(size_t)(a_2)*lddt  + (a_1))
#define d_ref(  a_1,a_2) ( d  +(size_t)(a_2)*ldd   + (a_1))

    int i, k, lddwork, old_i, old_ib, rows, cols;
    int ib;
    double _Complex one=1.;

    if (m < 0) {
        return -1;
    } else if (n < 0) {
        return -2;
    } else if (ldda2 < chameleon_max(1,m)) {
        return -4;
    }

    k = chameleon_min(m,n);
    if (k == 0) {
        hwork[0] = *((magmaDoubleComplex*) &one);
        return MAGMA_SUCCESS;
    }

    lddwork= m;

    /* lower parts of little T must be zero: memset all to 0 for simplicity */
    memset( t, 0, sizeof(magmaDoubleComplex) * nb * n );
    cudaMemset( dt, 0, sizeof(magmaDoubleComplex) * nb * n );

    //k = chameleon_min(m, nb); // m can be lower than IB
    /* copy the first diag tile of A1 from device to host: da1 -> d */
    cublasGetMatrix(nb, nb, sizeof(magmaDoubleComplex),
                    da1_ref(0, 0), ldda1,
                    d, nb);

    /* copy first panel of A2 from device to host: da2 -> a2 */
    cublasGetMatrix(nb, n, sizeof(magmaDoubleComplex),
                    da2_ref(0, 0), ldda2,
                    a2, nb);

    /* This is only blocked code for now */
    for (i = 0; i < m; i += nb) {

        ib = chameleon_min(m-i, nb);
        cols = n;

        /* Send the next panel (diagonal block of A1 & block column of A2)
           to the CPU (in work_a1 and work_a2) and apply tsmqr update on the
           remaining non updated panels */
        if (i>0) {

            /* copy the diag tile of A1 from device to host: da1 -> d */
            cublasGetMatrix(ib, ib, sizeof(magmaDoubleComplex),
                            da1_ref(i, i), ldda1,
                            d, ib);

            /* copy panel of A2 from device to host: da2 -> a2 */
            cublasGetMatrix(ib, cols, sizeof(magmaDoubleComplex),
                            da2_ref(i, 0), ldda2,
                            a2, ib);

            /* Apply H' to A(i+2*ib:m,i:n) from the left */
            rows = m-old_i-2*old_ib;
            if (rows > 0){
                CUDA_ztsmlq(
                        MagmaRight, MagmaConjTrans,
                        rows, old_ib, rows, cols, old_ib, old_ib,
                        da1_ref(old_i+2*old_ib, old_i), ldda1,
                        da2_ref(old_i+2*old_ib, 0), ldda2,
                        da2_ref(old_i, 0), ldda2,
                        dt_ref(0, old_i), lddt,
                        dwork, lddwork,
                        dwork + nb * lddwork, nb,
                        stream );
            }

        }

        /* compute LQ factorization of the panel of A2 ib x cols */
        CORE_ztslqt(ib, cols, ib,
                    (double _Complex*) d, ib,
                    (double _Complex*) a2, ib,
                    (double _Complex*) t, ib,
                    (double _Complex*) tau,
                    (double _Complex*) hwork);

        /* Send the panel from A2 back to the GPU */
        cublasSetMatrix(ib, cols, sizeof(magmaDoubleComplex),
                        a2, ib,
                        da2_ref(i, 0), ldda2);

        /* Send the triangular factor T from hwork to the GPU */
        cublasSetMatrix(ib, cols, sizeof(magmaDoubleComplex),
                        t, ib,
                        dt_ref(0, i), lddt);

        /* get back the diag tile in A1 from host to device: d -> da1 */
        cublasSetMatrix(ib, ib, sizeof(magmaDoubleComplex),
                        d, ib,
                        da1_ref(i, i), ldda1);

        /* tsmlq update on one panel forward (look ahead 1) */
        if (i + ib < m) {

            if (i+2*ib < m){
                rows = ib;
            }
            else{
                rows = m-i-ib;
            }

            /* Apply H' to A(i+ib:i+2*ib,i:n) from the right */
            CUDA_ztsmlq(
                    MagmaRight, MagmaConjTrans,
                    rows, ib, rows, cols, ib, ib,
                    da1_ref(i+ib, i), ldda1,
                    da2_ref(i+ib, 0), ldda2,
                    da2_ref(i, 0), ldda2,
                    dt_ref(0, i), lddt,
                    dwork, lddwork,
                    dwork + nb * lddwork, nb,
                    stream );
            cudaThreadSynchronize();
            old_i = i;
            old_ib = ib;
        }
    }

#undef da1_ref
#undef da2_ref
#undef a2_ref
#undef t_ref
#undef dt_ref
#undef d_ref

    return CHAMELEON_SUCCESS;

}
#endif
