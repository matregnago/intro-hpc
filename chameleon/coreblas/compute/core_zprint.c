/**
 *
 * @file core_zprint.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zprint CPU kernel
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "coreblas/lapacke.h"
#include "coreblas.h"

#define HEADER "--- %s / %p, %d:\n"
#define FOOTER "-------------------------\n"

void
CORE_zprint( FILE *file, const char *header,
             cham_uplo_t uplo, int M, int N,
             const CHAMELEON_Complex64_t *A, int lda )
{
    FILE   *output = (file == NULL) ? stdout : file;
    size_t  len = 0;
    int     i, j;

    len  = 4 + strlen(header) + 3 + 20 + 2 + 10 + 2;
#if defined(PRECISION_z) || defined(PRECISION_c)
    len += ( 2 + 31 * N + 1 ) * M;
#else
    len += ( 2 + 14 * N + 1 ) * M;
#endif
    len += strlen(FOOTER);

    if ( len <= 4096 ) {
        size_t size;
        char  *s = malloc( len );

        size = snprintf( s, len, HEADER, header, (void*)A, lda );
        for(i=0; i<M; i++) {
            size += snprintf( s + size, len - size, "  " );
            for(j=0; j<N; j++) {
#if defined(PRECISION_z) || defined(PRECISION_c)
                size += snprintf( s + size, len - size, " (% e, % e)",
                                  creal( A[j*lda + i] ),
                                  cimag( A[j*lda + i] ) );
#else
                size += snprintf( s + size, len - size, " % e",
                                  A[j*lda + i] );
#endif
            }
            size += snprintf( s + size, len - size, "\n" );
        }
        size += snprintf( s + size, len - size, FOOTER );

        fflush( output );
        fwrite( s, 1, size, output );
        fflush( output );
        free( s );
    }
    else {
        fprintf( output, HEADER, header, (void*)A, lda );
        for(i=0; i<M; i++) {
            fprintf( output, "  " );
            for(j=0; j<N; j++) {
#if defined(PRECISION_z) || defined(PRECISION_c)
                fprintf( output, " (% e, % e)",
                         creal( A[j*lda + i] ),
                         cimag( A[j*lda + i] ));
#else
                fprintf( output, " % e", A[j*lda + i] );
#endif
            }
            fprintf( output, "\n" );
        }
        fprintf( output, FOOTER );
        fflush( output );
    }

    (void)uplo;
}
