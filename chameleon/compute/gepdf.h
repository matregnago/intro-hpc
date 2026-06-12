/**
 *
 * @file gepdf.h
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon GEPDF auxiliary function.
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-02-18
 *
 */
#ifndef _chameleon_gepdf_h_
#define _chameleon_gepdf_h_

/**
 * @brief Compute the parameters of the QDWH iteration.
 *
 * This function is kept outside the generation process to be kept in double
 * precision. Single precision may induce incorrect computations of the
 * parameters and the number of iterations.
 *
 * @param[in] Li
 *        The estimated condition number on input.
 *
 * @param[out] params
 *        Array of size 3 that holds the parameters from the coming iteration
 *        computed from Li on exit.
 *
 * @retval The new estimated condition number after the iteration
 *
 */
static inline double
chameleon_gepdf_parameters( double Li, double *params )
{
    double L2, dd, sqd, a1;
    double a, b, c;

    L2  = Li * Li;
    /* Add this check because I add L2 that was slowly moving to 1.0 + eps, and generating NaN in the pow call */
    L2  = (L2 < 1.) ? L2 : 1.;
    dd  = pow( 4. * (1. - L2 ) / (L2 * L2), 1./3. );
    sqd = sqrt(1. + dd);
    a1  = sqd + sqrt( 8. - 4. * dd + 8. * (2. - L2) / (L2 * sqd) ) / 2.;
    a  = a1;
    b  = (a - 1.) * (a - 1.) / 4.;
    c  = a + b - 1.;

    params[0] = a;
    params[1] = b;
    params[2] = c;

    /* Return the new Li */
    return Li * (a + b * L2) / (1. + c * L2);
}

#endif /* _chameleon_gepdf_h_ */
