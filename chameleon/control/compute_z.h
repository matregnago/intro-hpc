/**
 *
 * @file compute_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon computational functions header
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Alycia Lisito
 * @author Matthieu Kuhn
 * @author Lionel Eyraud-Dubois
 * @author Ana Hourcau
 * @author Pierre Esterie
 * @author Matteo Marcos
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#ifndef _compute_z_h_
#define _compute_z_h_

/**
 * @brief Data structure to handle the GEMM workspaces
 */
struct chameleon_pzgemm_s {
    cham_gemm_t alg;
    CHAM_desc_t WA;
    CHAM_desc_t WB;
};

/**
 * @brief Data structure to handle the LASWP workspaces
 */
struct chameleon_pzlaswp_s {
    CHAM_desc_t   *Wu;                  /**< Workspace used for the row/column permutation        */
    CHAM_perm_t    ws;                  /**< Workspace used for the row/column permutation        */
    CHAM_reduce_t  reduce;              /**< Structure for reduction operations                   */
    int            batch_size_swap;     /**< Batch size for the permutation                       */
    int            batch_size_swap_get; /**< Batch size for the laswp get                         */
    int            batch_size_swap_set; /**< Batch size for the laswp set                         */
    int            allreduce;           /**< Specifies whether the reduction is replicated or not */
};

struct chameleon_pzgetrf_s;

typedef void (*chameleon_zgetrf_update_fct_t) (
    struct chameleon_pzgetrf_s *, CHAM_desc_t *, CHAM_ipiv_t *, int, int, RUNTIME_option_t *, RUNTIME_sequence_t * );

/**
 * @brief Data structure to handle the GETRF workspaces with partial pivoting
 */
struct chameleon_pzgetrf_s {
    struct chameleon_pzlaswp_s    *laswp;                /**< Structure containing the permutation workspace and the reduce data   */
    CHAM_desc_pivot_t              pivot;                /**< Structure containing the workspace used for the panel factorisation  */
    cham_getrf_t                   alg;                  /**< Define the algorithm used to compute the getrf                       */
    cham_bool_t                    backperm_enabled;     /**< Define if the backward perumtation is enabled or not                 */
    int                            ib;                   /**< Internal blocking parameter                                          */
    int                            batch_adaptive;       /**< Whether to use adaptative batch or not                               */
    int                            batch_size;           /**< Batch size                                                           */
    int                            batch_size_blas2;     /**< Batch size for the blas 2 operations of the panel factorization      */
    int                            batch_size_blas3;     /**< Batch size for the blas 3 operations of the panel factorization      */
    int                            ringswitch;           /**< Define when to switch to ring bcast                                  */
    cham_fixdbl_t                  flops_min;            /**< Define size of batched task in MFlops                                */
    chameleon_zgetrf_update_fct_t  panel_permute_update; /**< Function to apply the panel permute and panel update                 */
    CHAM_desc_t                   *Up;                   /**< Workspace used for the panel factorization                           */
    CHAM_desc_t                   *Wl;                   /**< Workspace used for the update                                        */
};

/**
 * @brief   Data structure to handle the GETRF temporary workspaces
 *          for MPI transfers.
 *
 * @comment The idea is to manage explicitely temporary
 *          blocks arising from MPI transfers automatically
 *          inferred by StarPU, hence limiting the total number
 *          of temporary data allocated for these blocks.
 *
 *          The blocks to be sent/received on the network are
 *          copied into those buffers. These copies are
 *          then used by the algorithm in place of the regular
 *          blocks of the problem matrix.
 *
 *          For WL (resp. WU), the number of allocated blocks
 *          corresponds to the number of blocks on the column
 *          (resp. on the line) multiplied by lookahead number
 *          from the current chameleon context.
 *
 *          Then, depending on the block panel index, we access
 *          one of the temporary column blocks of WL and row blocks
 *          of WU int a circular way.
 *
 *          For instance, for the block panel index k, the block
 *          A(m,k) produced by the TRSM(A(k,k),A(m,k)) is stored
 *          into temporary buffer WL(m,k%chamctxt->lookahead).
 *          Similarly, the block A(k,n) is stored into the temporary
 *          block WU(k%chamctxt->lookahead, n).
 *
 *          Notice that, by doing so, the notion of look ahead is
 *          reintroduced : artificial dependencies are implied by
 *          the circular usage of WL and WU temporary workspaces.
 *
 */
struct chameleon_pzgetrf_nopiv_s {
    int use_workspace;

    CHAM_desc_t WL; /**> Workspace to store temporary blocks of the
                         diagonal and the lower part of the problem matrix */
    CHAM_desc_t WU; /**> Workspace to store temporary blocks of the
                         upper part of the problem matrix                  */
};

/**
 * @brief Data structure to handle the Centering-Scaled workspaces
 */
struct chameleon_pzcesca_s {
    CHAM_desc_t Wgcol;
    CHAM_desc_t Wgrow;
    CHAM_desc_t Wgelt;
    CHAM_desc_t Wdcol;
    CHAM_desc_t Wdrow;
};

/**
 * @brief Data structure to handle the GRAM workspaces
 */
struct chameleon_pzgram_s {
    CHAM_desc_t Wcol;
    CHAM_desc_t Welt;
};

/**
 *  Declarations of internal sequential functions
 */
int chameleon_zshift( CHAM_context_t *chamctxt, int m, int n, CHAMELEON_Complex64_t *A,
                      int nprob, int me, int ne, int L,
                      RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );

/**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 */
#if defined(PRECISION_z) || defined(PRECISION_d)
void chameleon_pzgered( cham_uplo_t uplo, double prec, CHAM_desc_t *A,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzhered( cham_trans_t trans, cham_uplo_t uplo, double prec, CHAM_desc_t *A,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzgerst( cham_uplo_t uplo, CHAM_desc_t *A,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
#endif
int chameleon_pzgebrd( int genD, cham_job_t jobu, cham_job_t jobvt,
                       CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                       CHAMELEON_Complex64_t *U, int LDU, CHAMELEON_Complex64_t *VT, int LDVT,
                       double *E, double *S,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzgelqf( int genD, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgelqfrh( int genD, int BS, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                          RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgenm2( double tol, CHAM_desc_t *A, double *result,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzgemm( struct chameleon_pzgemm_s *options, cham_trans_t transA, cham_trans_t transB, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgepdf_qdwh( cham_mtxtype_t trans, CHAM_desc_t *descU, CHAM_desc_t *descH, gepdf_info_t *info,
                             RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzgepdf_qr( int genD, int doqr, int optid, const libhqr_tree_t *qrtreeT, const libhqr_tree_t *qrtreeB, CHAM_desc_t *A1, CHAM_desc_t *TS1, CHAM_desc_t *TT1, CHAM_desc_t *D1, CHAM_desc_t *Q1, CHAM_desc_t *A2, CHAM_desc_t *TS2, CHAM_desc_t *TT2, CHAM_desc_t *D2, CHAM_desc_t *Q2,
                           RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzgeqrf( int genD, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgeqrfrh( int genD, int BS, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                          RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgetrf( struct chameleon_pzgetrf_s *ws, CHAM_desc_t *A, CHAM_ipiv_t *IPIV,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzgetrf_incpiv(CHAM_desc_t *A, CHAM_desc_t *L, CHAM_desc_t *D, int *IPIV,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgetrf_nopiv(struct chameleon_pzgetrf_nopiv_s *ws, CHAM_desc_t *A,
                             RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgetrf_reclap(CHAM_desc_t *A, int *IPIV,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgetrf_rectil(CHAM_desc_t *A, int *IPIV,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgetrf_panel_permute_update( struct chameleon_pzgetrf_s *ws,
                                             CHAM_desc_t                *A,
                                             CHAM_ipiv_t                *ipiv,
                                             int                         k,
                                             int                         n,
                                             RUNTIME_option_t           *options,
                                             RUNTIME_sequence_t         *sequence );
void chameleon_pzgetrf_panel_update( struct chameleon_pzgetrf_s *ws,
                                     CHAM_desc_t                *A,
                                     CHAM_ipiv_t                *ipiv,
                                     int                         k,
                                     int                         n,
                                     RUNTIME_option_t           *options,
                                     RUNTIME_sequence_t         *sequence );
void chameleon_pzhegst(int itype, cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzhemm( struct chameleon_pzgemm_s *ws,cham_side_t side, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzherk(cham_uplo_t uplo, cham_trans_t trans, double alpha, CHAM_desc_t *A, double beta, CHAM_desc_t *C,
                      RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzher2k(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, double beta, CHAM_desc_t *C,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzhetrd_he2hb(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *E,
                             RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzlacpy(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzlag2c(CHAM_desc_t *A, CHAM_desc_t *SB,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzlange_generic( cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A,
                                double *result,
                                RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzlansy_generic( cham_normtype_t norm, cham_uplo_t uplo, cham_trans_t trans,
                                CHAM_desc_t *A, double *result,
                                RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzlascal(cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzlaset( cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta, CHAM_desc_t *A,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzlaset2(cham_uplo_t uplo, CHAMELEON_Complex64_t alpha,                          CHAM_desc_t *A,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzlaswp( struct chameleon_pzlaswp_s *ws, cham_dir_t dir, CHAM_desc_t *A, CHAM_ipiv_t *IPIV,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzlaswp_panel( struct chameleon_pzlaswp_s *ws, cham_bool_t inplace, cham_dir_t dir, CHAM_desc_t *A, CHAM_ipiv_t *ipiv,
                              int k, int n, RUNTIME_option_t *options, RUNTIME_sequence_t *sequence );
void chameleon_pzlaswpc( struct chameleon_pzlaswp_s *ws, cham_dir_t dir, CHAM_desc_t *A, CHAM_ipiv_t *IPIV,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzlatms( cham_dist_t idist, unsigned long long int seed, cham_sym_t sym, double *D, int mode, double cond, double dmax, CHAM_desc_t *A,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzlatro( cham_uplo_t uplo, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *B,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzlauum(cham_uplo_t uplo, CHAM_desc_t *A,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzplghe(double bump, cham_uplo_t uplo, CHAM_desc_t *A, int bigM, int m0, int n0, unsigned long long int seed,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzplgsy(CHAMELEON_Complex64_t bump, cham_uplo_t uplo, CHAM_desc_t *A, int bigM, int m0, int n0, unsigned long long int seed,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzplrnt(CHAM_desc_t *A, int bigM, int m0, int n0, unsigned long long int seed,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzplrnk(int K, CHAM_desc_t *C, unsigned long long int seedA, unsigned long long int seedB,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzpotrf(cham_uplo_t uplo, CHAM_desc_t *A,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzpotrimm(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *C,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzshift(int, int, int, CHAMELEON_Complex64_t *, int *, int, int, int,
                       RUNTIME_sequence_t*, RUNTIME_request_t*);
void chameleon_pzsymm( struct chameleon_pzgemm_s *ws,cham_side_t side, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzsyrk(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                      RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzsyr2k(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzsytrf(cham_uplo_t uplo, CHAM_desc_t *A,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pztile2band(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *descAB,
                           RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pztpgqrt( int kt, int L, CHAM_desc_t *V2, CHAM_desc_t *T2, CHAM_desc_t *Q1, CHAM_desc_t *Q2,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pztpqrt( int L, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pztradd(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAMELEON_Complex64_t beta, CHAM_desc_t *B,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pztrmm(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                      RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pztrsm(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                      RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pztrsmpl(CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *L, int *IPIV,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pztrsmrv(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *W,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pztrtri(cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzungbr(cham_side_t side, CHAM_desc_t *A, CHAM_desc_t *O, CHAM_desc_t *T,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzungbrrh(cham_side_t side, CHAM_desc_t *A, CHAM_desc_t *O, CHAM_desc_t *T,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzungqr( int genD, CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, CHAM_desc_t *D,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzungqrrh( int genD, int BS, CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, CHAM_desc_t *D,
                          RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzunglq( int genD, CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, CHAM_desc_t *D,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzunglqrh( int genD, int BS, CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, CHAM_desc_t *D,
                          RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzungtr(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzunmqr( int genD, cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, CHAM_desc_t *D,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzunmqrrh( int genD, int BS, cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, CHAM_desc_t *D,
                          RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzunmlq( int genD, cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, CHAM_desc_t *D,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzunmlqrh( int genD, int BS, cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, CHAM_desc_t *D,
                          RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int chameleon_pzgelqf_step( int genD, int k, int ib,
                            CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                            RUNTIME_option_t *options,
                            RUNTIME_sequence_t *sequence );
int chameleon_pzgeqrf_step( int genD, int k, int ib,
                            CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                            RUNTIME_option_t *options,
                            RUNTIME_sequence_t *sequence );
int  chameleon_pzgelqf_param_step( int genD, cham_uplo_t uplo, int k, int ib,
                                   const libhqr_tree_t *qrtree, int *tiles,
                                   CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                                   RUNTIME_option_t *options,
                                   RUNTIME_sequence_t *sequence );
int  chameleon_pzgeqrf_param_step( int genD, cham_uplo_t uplo, int k, int ib,
                                   const libhqr_tree_t *qrtree, int *tiles,
                                   CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                                   RUNTIME_option_t *options,
                                   RUNTIME_sequence_t *sequence );
void chameleon_pzungqr_param_step( int genD, cham_uplo_t uplo, int k, int ib,
                                   const libhqr_tree_t *qrtree, int nbtiles, int *tiles,
                                   CHAM_desc_t *A, CHAM_desc_t *Q,
                                   CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                                   RUNTIME_option_t *options,
                                   RUNTIME_sequence_t *sequence );
void chameleon_pzgelqf_param( int genD, int K, const libhqr_tree_t *qrtree,
                              CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgeqrf_param( int genD, int K, const libhqr_tree_t *qrtree,
                              CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzunmlq_param( int genD, const libhqr_tree_t *qrtree, cham_side_t side, cham_trans_t trans,
                              CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzunmqr_param( int genD, const libhqr_tree_t *qrtree, cham_side_t side, cham_trans_t trans,
                              CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzunglq_param( int genD, const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *Q,
                              CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzungqr_param( int genD, int K, const libhqr_tree_t *qrtree,
                              CHAM_desc_t *A, CHAM_desc_t *Q,
                              CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pztpgqrt_param( int genD, cham_uplo_t uplo, int kt, const libhqr_tree_t *qrtree,
                               CHAM_desc_t *V2, CHAM_desc_t *Q1, CHAM_desc_t *Q2,
                               CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *DD,
                               RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pztpqrt_param( int genD, cham_uplo_t uplo, int K,
                              const libhqr_tree_t *qrtree,
                              CHAM_desc_t *ATop, CHAM_desc_t *A,
                              CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );

/**
 * Centered-Scaled function prototypes
 */
void chameleon_pzcesca( struct chameleon_pzcesca_s *ws, int center, int scale, cham_store_t axis, CHAM_desc_t *A,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
/**
 * Gram function prototypes
 */
void chameleon_pzgram( struct chameleon_pzgram_s *ws, cham_uplo_t uplo, CHAM_desc_t *A,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );

/**
 * Specific functions called when option is already initialized
 */
void chameleon_pzbcast_tile( cham_store_t dir, cham_bcast_t algo, const CHAM_desc_t *A, int Am, int An, const CHAM_desc_t *W, int Wm, int Wn, RUNTIME_option_t *options );

/**
 *  Macro for matrix conversion / Lapack interface
 */
static inline int
chameleon_zdesc_alloc_diag( const CHAM_context_t *chamctxt,
                            CHAM_desc_t *descA, const char *name,
                            int nb, int m, int n, int p, int q )
{
    int   rc, diag_m = chameleon_min( m, n );
    char *subname;
    chameleon_asprintf( &subname, "%s_diag", name );
    rc = chameleon_desc_init( chamctxt, descA, subname, CHAMELEON_MAT_ALLOC_TILE,
                                ChamComplexDouble, nb, nb, diag_m, nb, diag_m, nb, p, q,
                                chameleon_getaddr_diag,
                                chameleon_getblkldd_ccrb,
                                chameleon_getrankof_2d_diag, NULL );
    free( subname );
    return rc;
}

#define chameleon_zdesc_alloc( descA, name, mb, nb, m, n, free)         \
    {                                                                   \
        int rc;                                                         \
        rc = chameleon_desc_init_local( chamctxt, &(descA), name,       \
                                        ChamComplexDouble,              \
                                        (mb), (nb),                     \
                                        (m), (n) );                     \
        if ( rc != CHAMELEON_SUCCESS ) {                                \
            {free;}                                                     \
            return rc;                                                  \
        }                                                               \
    }

/**
 * @brief Create a copy of a descriptor restricted to a smaller size.
 * @param[in]  descIn  The input descriptor from which the structure should be copied.
 * @param[out] descOut The output descriptor that is a copy of the input one with allocation on the fly.
 * @param[in]  m       The number of rows of the output descriptor.
 * @param[in]  n       The number of columns of the output descriptor.
 * @return CHAMELEON_SUCCESS on success, the associated error on failure.
 */
static inline int
chameleon_zdesc_copy_and_restrict( const CHAM_context_t *chamctxt,
                                   const CHAM_desc_t *descIn,
                                   CHAM_desc_t *descOut,
                                   int m, int n )
{
    int rc;
    char *subname;
    chameleon_asprintf( &subname, "%s_restricted", descIn->name );
    rc = chameleon_desc_init( chamctxt, descOut, subname, CHAMELEON_MAT_ALLOC_TILE,
                              ChamComplexDouble, descIn->mb, descIn->nb,
                              m, n, m, n,
                              chameleon_desc_datadist_get_iparam(descIn, 0),
                              chameleon_desc_datadist_get_iparam(descIn, 1),
                              NULL, NULL,
                              descIn->get_rankof_init,
                              descIn->get_rankof_init_arg );
    free( subname );
    return rc;
}

/**
 * @brief Internal function to convert the lapack format to tile format in
 * LAPACK interface calls
 */
static inline int
chameleon_zlap2tile( const CHAM_context_t *chamctxt, const char *name,
                     CHAM_desc_t *descAl, CHAM_desc_t *descAt,
                     cham_desc_access_t mode, cham_uplo_t uplo,
                     CHAMELEON_Complex64_t *A, int mb, int nb, int lm, int ln, int m, int n,
                     RUNTIME_sequence_t *seq, RUNTIME_request_t *req )
{
    char *fullname;
    if ( CHAMELEON_TRANSLATION == ChamOutOfPlace ) {
        /* Initialize the Lapack descriptor */
        chameleon_asprintf( &fullname, "%sl", name );
        chameleon_desc_init( chamctxt, descAl, fullname, A, ChamComplexDouble, mb, nb,
                             lm, ln, m, n, 1, 1,
                             chameleon_getaddr_cm, chameleon_getblkldd_cm, NULL, NULL );
        descAl->styp = ChamCM;
        free( fullname );

        /* Initialize the tile descriptor */
        chameleon_asprintf( &fullname, "%st", name );
        chameleon_desc_init( chamctxt, descAt, fullname, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, mb, nb,
                             lm, ln, m, n, 1, 1,
                             chameleon_getaddr_ccrb, chameleon_getblkldd_ccrb, NULL, NULL );
        free( fullname );

        if ( mode & ChamDescInput ) {
            chameleon_pzlacpy( uplo, descAl, descAt, seq, req );
        }
    }
    else {
        /* Initialize the tile descriptor */
        chameleon_asprintf( &fullname, "%st", name );
        chameleon_desc_init( chamctxt, descAt, fullname, A, ChamComplexDouble, mb, nb,
                             lm, ln, m, n, 1, 1,
                             chameleon_getaddr_cm, chameleon_getblkldd_cm, NULL, NULL );
        free( fullname );
    }
    return CHAMELEON_SUCCESS;
}

/**
 * @brief Internal function to convert back the tile format to the lapack format
 * in LAPACK interface calls
 */
static inline int
chameleon_ztile2lap( const CHAM_context_t *chamctxt, CHAM_desc_t *descAl, CHAM_desc_t *descAt,
                     cham_desc_access_t mode, cham_uplo_t uplo,
                     RUNTIME_sequence_t *seq, RUNTIME_request_t *req )
{
    if ( CHAMELEON_TRANSLATION == ChamOutOfPlace ) {
        if ( mode & ChamDescOutput ) {
            chameleon_pzlacpy( uplo, descAt, descAl, seq, req );
        }
        RUNTIME_desc_flush( descAl, seq );
    }
    RUNTIME_desc_flush( descAt, seq );

    return CHAMELEON_SUCCESS;
}

/**
 * @brief Internal function to cleanup the temporary data from the layout
 * conversions in LAPACK interface calls
 */
static inline void
chameleon_ztile2lap_cleanup( const CHAM_context_t *chamctxt, CHAM_desc_t *descAl, CHAM_desc_t *descAt )
{
    (void)chamctxt;
    if ( CHAMELEON_TRANSLATION == ChamOutOfPlace ) {
        chameleon_desc_destroy( descAl );
    }
    chameleon_desc_destroy( descAt );
}

#endif /* _compute_z_h_ */
