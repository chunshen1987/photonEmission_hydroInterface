#include "gauss_quadrature.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

int gauss_quadrature(
    int order, int kind, double alpha, double beta, double a, double b,
    double x[], double w[]) {
    double eps = 1e-100;  // avoid nan
    // Construct the rule.
    switch (kind) {
        case 1:  // Legendre,  1.0, integration interval: (a, b)
            alpha = 0.0e0;
            beta = 0.0e0;
            break;
        case 2:  // Chebyshev Type 1, ((b-x)*(x-a))^-0.5), integration interval:
                 // (a, b)
            alpha = 0.0e0;
            beta = 0.0e0;
            break;
        case 3:  // Gegenbauer, ((b-x)*(x-a))^alpha  integration interval: (a,
                 // b)
            // alpha is a real parameter greater than -1
            beta = 0.0e0;
            break;
        case 4:  // Jacobi, (b-x)^alpha*(x-a)^beta integration interval: (a, b)
            // alpha is the exponent of (1-x), which must be greater than -1.
            // beta is the exponent of (1+x), which must be greater than -1.
            break;
        case 5:  // Generalized Laguerre, (x-a)^alpha*exp(-b*(x-a)) integration
                 // interval: (a,inf)
            // b is the scale factor in the exponential, and is typically 1.
            // alpha is a real parameter greater than -1
            alpha = 0.0e0;  // The generalized version does not work currently
                            // (Sept. 19, 2012)
            beta = 0.0e0;
            break;
        case 6:  // Generalized Hermite, |x-a|^alpha*exp(-b*(x-a)^2) integration
                 // interval: (-inf, inf)
            // alpha is the parameter for the generalized Gauss-Hermite
            // quadrature rule. The value of alpha may be any real value greater
            // than -1.0. Specifying alpha=0.0 results in the basic
            // (non-generalized) rule.
            beta = 0.0e0;
            a = 0.0;  // the center point
            // b is the scale factor for the gaussian width
            break;
        case 7:  // Exponential, |x-(a+b)/2.0|^alpha, integration interval:
                 // (a,b)
            // alpha is a real parameter greater than -1
            beta = 0.0e0;
            break;
        case 8:  // Rational (x-a)^alpha*(x+b)^beta, integration interval: (a,
                 // inf)
            cout << " mode 8 is not support yet!" << endl;
            exit(1);
            break;
        case 9:  // Chebyshev Type 2, ((b-x)*(x-a))^(+0.5), integration interval
                 // (a, b)
            alpha = 0.0e0;
            beta = 0.0e0;
            break;
        default:
            cout << " there is no rule of type: " << kind << endl;
            exit(1);
    }

    // calculate the gauss-quadrature points and weights
    cgqf(order, kind, alpha, beta, a, b, x, w);

    // scale the weights by their correspond functional form
    for (int i = 0; i < order; i++) {
        switch (kind) {
            case 1:  // Legendre,  1.0, integration interval: (a, b)
                break;
            case 2:  // Chebyshev Type 1, ((b-x)*(x-a))^-0.5), integration
                     // interval: (a, b)
                w[i] = w[i] * sqrt((b - x[i]) * (x[i] - a));
                break;
            case 3:  // Gegenbauer, ((b-x)*(x-a))^alpha  integration interval:
                     // (a, b)
                w[i] = w[i] / (pow((b - x[i]) * (x[i] - a), alpha) + eps);
                break;
            case 4:  // Jacobi, (b-x)^alpha*(x-a)^beta integration interval: (a,
                     // b)
                w[i] = w[i] / (pow((b - x[i]), alpha) + eps)
                       / (pow((x[i] - a), beta) + eps);
                break;
            case 5:  // Generalized Laguerre, (x-a)^alpha*exp(-b*(x-a))
                     // integration interval: (a,inf)
                w[i] = w[i] / (pow(fabs(x[i] - a), alpha) + eps)
                       / (exp(-b * (x[i] - a)) + eps);
                break;
            case 6:  // Generalized Hermite, |x-a|^alpha*exp(-b*(x-a)^2)
                     // integration interval: (-inf, inf)
                w[i] = w[i] / (pow(fabs(x[i] - a), alpha) + eps)
                       / (exp(-b * (x[i] - a) * (x[i] - a)) + eps);
                break;
            case 7:  // Exponential, |x-(a+b)/2.0|^alpha, integration interval:
                     // (a,b)
                w[i] = w[i] / (pow(fabs(x[i] - (a + b) / 2.0), alpha) + eps);
                break;
            case 8:  // Rational (x-a)^alpha*(x+b)^beta, integration interval:
                     // (a, inf)
                w[i] = w[i] / (pow((x[i] - a), alpha) + eps)
                       / (pow((x[i] + b), beta) + eps);
                break;
            case 9:  // Chebyshev Type 2, ((b-x)*(x-a))^(+0.5), integration
                     // interval (a, b)
                w[i] = w[i] / (sqrt((b - x[i]) * (x[i] - a)) + eps);
                break;
            default:
                cout << " there is no rule of type: " << kind << endl;
                exit(1);
        }
    }
    return 0;
}

int gauss_quadrature_standard(
    int order, int kind, double alpha, double beta, double a, double b,
    double x[], double w[]) {
    //  Compute the Gauss quadrature formula for default values of A and B.
    //
    cdgqf(order, kind, alpha, beta, x, w);

    return 0;
}

int scale_gausspoints(
    int order, int kind, double alpha, double beta, double a, double b,
    double x[], double w[]) {
    double eps = 1e-100;  // avoid nan
    // Construct the rule.
    switch (kind) {
        case 1:  // Legendre,  1.0, integration interval: (a, b)
            alpha = 0.0e0;
            beta = 0.0e0;
            break;
        case 2:  // Chebyshev Type 1, ((b-x)*(x-a))^-0.5), integration interval:
                 // (a, b)
            alpha = 0.0e0;
            beta = 0.0e0;
            break;
        case 3:  // Gegenbauer, ((b-x)*(x-a))^alpha  integration interval: (a,
                 // b)
            // alpha is a real parameter greater than -1
            beta = 0.0e0;
            break;
        case 4:  // Jacobi, (b-x)^alpha*(x-a)^beta integration interval: (a, b)
            // alpha is the exponent of (1-x), which must be greater than -1.
            // beta is the exponent of (1+x), which must be greater than -1.
            break;
        case 5:  // Generalized Laguerre, (x-a)^alpha*exp(-b*(x-a)) integration
                 // interval: (a,inf)
            // b is the scale factor in the exponential, and is typically 1.
            // alpha is a real parameter greater than -1
            beta = 0.0e0;
            break;
        case 6:  // Generalized Hermite, |x-a|^alpha*exp(-b*(x-a)^2) integration
                 // interval: (-inf, inf)
            // alpha is the parameter for the generalized Gauss-Hermite
            // quadrature rule. The value of alpha may be any real value greater
            // than -1.0. Specifying alpha=0.0 results in the basic
            // (non-generalized) rule.
            beta = 0.0e0;
            a = 0.0;  // the center point
            // b is the scale factor for the gaussian width
            break;
        case 7:  // Exponential, |x-(a+b)/2.0|^alpha, integration interval:
                 // (a,b)
            // alpha is a real parameter greater than -1
            beta = 0.0e0;
            break;
        case 8:  // Rational (x-a)^alpha*(x+b)^beta, integration interval: (a,
                 // inf)
            cout << " mode 8 is not support yet!" << endl;
            exit(1);
            break;
        case 9:  // Chebyshev Type 2, ((b-x)*(x-a))^(+0.5), integration interval
                 // (a, b)
            alpha = 0.0e0;
            beta = 0.0e0;
            break;
        default:
            cout << " there is no rule of type: " << kind << endl;
            exit(1);
    }

    // calculate the gauss-quadrature points and weights
    cgqf_onlyscale(order, kind, alpha, beta, a, b, x, w);

    // scale the weights by their correspond functional form
    for (int i = 0; i < order; i++) {
        switch (kind) {
            case 1:  // Legendre,  1.0, integration interval: (a, b)
                break;
            case 2:  // Chebyshev Type 1, ((b-x)*(x-a))^-0.5), integration
                     // interval: (a, b)
                w[i] = w[i] * sqrt((b - x[i]) * (x[i] - a));
                break;
            case 3:  // Gegenbauer, ((b-x)*(x-a))^alpha  integration interval:
                     // (a, b)
                w[i] = w[i] / (pow((b - x[i]) * (x[i] - a), alpha) + eps);
                break;
            case 4:  // Jacobi, (b-x)^alpha*(x-a)^beta integration interval: (a,
                     // b)
                w[i] = w[i] / (pow((b - x[i]), alpha) + eps)
                       / (pow((x[i] - a), beta) + eps);
                break;
            case 5:  // Generalized Laguerre, (x-a)^alpha*exp(-b*(x-a))
                     // integration interval: (a,inf)
                w[i] = w[i] / (pow(fabs(x[i] - a), alpha) + eps)
                       / (exp(-b * (x[i] - a)) + eps);
                break;
            case 6:  // Generalized Hermite, |x-a|^alpha*exp(-b*(x-a)^2)
                     // integration interval: (-inf, inf)
                w[i] = w[i] / (pow(fabs(x[i] - a), alpha) + eps)
                       / (exp(-b * (x[i] - a) * (x[i] - a)) + eps);
                break;
            case 7:  // Exponential, |x-(a+b)/2.0|^alpha, integration interval:
                     // (a,b)
                w[i] = w[i] / (pow(fabs(x[i] - (a + b) / 2.0), alpha) + eps);
                break;
            case 8:  // Rational (x-a)^alpha*(x+b)^beta, integration interval:
                     // (a, inf)
                w[i] = w[i] / (pow((x[i] - a), alpha) + eps)
                       / (pow((x[i] + b), beta) + eps);
                break;
            case 9:  // Chebyshev Type 2, ((b-x)*(x-a))^(+0.5), integration
                     // interval (a, b)
                w[i] = w[i] / (sqrt((b - x[i]) * (x[i] - a)) + eps);
                break;
            default:
                cout << " there is no rule of type: " << kind << endl;
                exit(1);
        }
    }
    return 0;
}

//****************************************************************************80

void cdgqf(
    int nt, int kind, double alpha, double beta, double t[], double wts[])

//****************************************************************************80
//
//  Purpose:
//
//    CDGQF computes a Gauss quadrature formula with default A, B and simple
//    knots.
//
//  Discussion:
//
//    This routine computes all the knots and weights of a Gauss quadrature
//    formula with a classical weight function with default values for A and B,
//    and only simple knots.
//
//    There are no moments checks and no printing is done.
//
//    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Output, double T[NT], the knots.
//
//    Output, double WTS[NT], the weights.
//
{
    double *aj;
    double *bj;
    double zemu;

    parchk(kind, 2 * nt, alpha, beta);
    //
    //  Get the Jacobi matrix and zero-th moment.
    //
    aj = new double[nt];
    bj = new double[nt];

    zemu = class_matrix(kind, nt, alpha, beta, aj, bj);
    //
    //  Compute the knots and weights.
    //
    sgqf(nt, aj, bj, zemu, t, wts);

    delete[] aj;
    delete[] bj;

    return;
}
//****************************************************************************80

void cgqf(
    int nt, int kind, double alpha, double beta, double a, double b, double t[],
    double wts[])

//****************************************************************************80
//
//  Purpose:
//
//    CGQF computes knots and weights of a Gauss quadrature formula.
//
//  Discussion:
//
//    The user may specify the interval (A,B).
//
//    Only simple knots are produced.
//
//    Use routine EIQFS to evaluate this quadrature formula.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Input, double A, B, the interval endpoints, or
//    other parameters.
//
//    Output, double T[NT], the knots.
//
//    Output, double WTS[NT], the weights.
//
{
    int i;
    int *mlt;
    int *ndx;
    //
    //  Compute the Gauss quadrature formula for default values of A and B.
    //
    cdgqf(nt, kind, alpha, beta, t, wts);
    //
    //  Prepare to scale the quadrature formula to other weight function with
    //  valid A and B.
    //
    mlt = new int[nt];
    for (i = 0; i < nt; i++) {
        mlt[i] = 1;
    }
    ndx = new int[nt];
    for (i = 0; i < nt; i++) {
        ndx[i] = i + 1;
    }
    scqf(nt, t, mlt, wts, nt, ndx, wts, t, kind, alpha, beta, a, b);

    delete[] mlt;
    delete[] ndx;

    return;
}

//****************************************************************************80

void cgqf_onlyscale(
    int nt, int kind, double alpha, double beta, double a, double b, double t[],
    double wts[]) {
    int i;
    int *mlt;
    int *ndx;

    //  Prepare to scale the quadrature formula to other weight function with
    //  valid A and B.
    //
    mlt = new int[nt];
    for (i = 0; i < nt; i++) {
        mlt[i] = 1;
    }
    ndx = new int[nt];
    for (i = 0; i < nt; i++) {
        ndx[i] = i + 1;
    }
    scqf(nt, t, mlt, wts, nt, ndx, wts, t, kind, alpha, beta, a, b);

    delete[] mlt;
    delete[] ndx;

    return;
}

//****************************************************************************80

double class_matrix(
    int kind, int m, double alpha, double beta, double aj[], double bj[])

//****************************************************************************80
//
//  Purpose:
//
//    CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
//
//  Discussion:
//
//    This routine computes the diagonal AJ and sub-diagonal BJ
//    elements of the order M tridiagonal symmetric Jacobi matrix
//    associated with the polynomials orthogonal with respect to
//    the weight function specified by KIND.
//
//    For weight functions 1-7, M elements are defined in BJ even
//    though only M-1 are needed.  For weight function 8, BJ(M) is
//    set to zero.
//
//    The zero-th moment of the weight function is returned in ZEMU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
//
//    Input, int M, the order of the Jacobi matrix.
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Output, double AJ[M], BJ[M], the diagonal and subdiagonal
//    of the Jacobi matrix.
//
//    Output, double CLASS_MATRIX, the zero-th moment.
//
{
    double a2b2;
    double ab;
    double aba;
    double abi;
    double abj;
    double abti;
    double apone;
    int i;
    double pi = 3.14159265358979323846264338327950;
    double temp;
    double temp2;
    double zemu = 0.;

    temp = r8_epsilon();

    parchk(kind, 2 * m - 1, alpha, beta);

    temp2 = 0.5;

    if (500.0 * temp < r8_abs(pow(tgamma(temp2), 2) - pi)) {
        cout << "\n";
        cout << "CLASS_MATRIX - Fatal error!\n";
        cout << "  Gamma function does not match machine parameters.\n";
        exit(1);
    }

    if (kind == 1) {
        ab = 0.0;

        zemu = 2.0 / (ab + 1.0);

        for (i = 0; i < m; i++) {
            aj[i] = 0.0;
        }

        for (i = 1; i <= m; i++) {
            abi = i + ab * (i % 2);
            abj = 2 * i + ab;
            bj[i - 1] = sqrt(abi * abi / (abj * abj - 1.0));
        }
    } else if (kind == 2) {
        zemu = pi;

        for (i = 0; i < m; i++) {
            aj[i] = 0.0;
        }

        bj[0] = sqrt(0.5);
        for (i = 1; i < m; i++) {
            bj[i] = 0.5;
        }
    } else if (kind == 3) {
        ab = alpha * 2.0;
        zemu =
            pow(2.0, ab + 1.0) * pow(tgamma(alpha + 1.0), 2) / tgamma(ab + 2.0);

        for (i = 0; i < m; i++) {
            aj[i] = 0.0;
        }

        bj[0] = sqrt(1.0 / (2.0 * alpha + 3.0));
        for (i = 2; i <= m; i++) {
            bj[i - 1] = sqrt(i * (i + ab) / (4.0 * pow(i + alpha, 2) - 1.0));
        }
    } else if (kind == 4) {
        ab = alpha + beta;
        abi = 2.0 + ab;
        zemu = pow(2.0, ab + 1.0) * tgamma(alpha + 1.0) * tgamma(beta + 1.0)
               / tgamma(abi);
        aj[0] = (beta - alpha) / abi;
        bj[0] = sqrt(
            4.0 * (1.0 + alpha) * (1.0 + beta) / ((abi + 1.0) * abi * abi));
        a2b2 = beta * beta - alpha * alpha;

        for (i = 2; i <= m; i++) {
            abi = 2.0 * i + ab;
            aj[i - 1] = a2b2 / ((abi - 2.0) * abi);
            abi = abi * abi;
            bj[i - 1] = sqrt(
                4.0 * i * (i + alpha) * (i + beta) * (i + ab)
                / ((abi - 1.0) * abi));
        }
    } else if (kind == 5) {
        zemu = tgamma(alpha + 1.0);

        for (i = 1; i <= m; i++) {
            aj[i - 1] = 2.0 * i - 1.0 + alpha;
            bj[i - 1] = sqrt(i * (i + alpha));
        }
    } else if (kind == 6) {
        zemu = tgamma((alpha + 1.0) / 2.0);

        for (i = 0; i < m; i++) {
            aj[i] = 0.0;
        }

        for (i = 1; i <= m; i++) {
            bj[i - 1] = sqrt((i + alpha * (i % 2)) / 2.0);
        }
    } else if (kind == 7) {
        ab = alpha;
        zemu = 2.0 / (ab + 1.0);

        for (i = 0; i < m; i++) {
            aj[i] = 0.0;
        }

        for (i = 1; i <= m; i++) {
            abi = i + ab * (i % 2);
            abj = 2 * i + ab;
            bj[i - 1] = sqrt(abi * abi / (abj * abj - 1.0));
        }
    } else if (kind == 8) {
        ab = alpha + beta;
        zemu = tgamma(alpha + 1.0) * tgamma(-(ab + 1.0)) / tgamma(-beta);
        apone = alpha + 1.0;
        aba = ab * apone;
        aj[0] = -apone / (ab + 2.0);
        bj[0] = -aj[0] * (beta + 1.0) / (ab + 2.0) / (ab + 3.0);
        for (i = 2; i <= m; i++) {
            abti = ab + 2.0 * i;
            aj[i - 1] = aba + 2.0 * (ab + i) * (i - 1);
            aj[i - 1] = -aj[i - 1] / abti / (abti - 2.0);
        }

        for (i = 2; i <= m - 1; i++) {
            abti = ab + 2.0 * i;
            bj[i - 1] = i * (alpha + i) / (abti - 1.0) * (beta + i)
                        / (abti * abti) * (ab + i) / (abti + 1.0);
        }
        bj[m - 1] = 0.0;
        for (i = 0; i < m; i++) {
            bj[i] = sqrt(bj[i]);
        }
    } else if (kind == 9) {
        zemu = pi / 2.0;

        for (i = 0; i < m; i++) {
            aj[i] = 0.0;
        }

        for (i = 0; i < m; i++) {
            bj[i] = 0.5;
        }
    }

    return zemu;
}
//****************************************************************************80

void imtqlx(int n, double d[], double e[], double z[])

//****************************************************************************80
//
//  Purpose:
//
//    IMTQLX diagonalizes a symmetric tridiagonal matrix.
//
//  Discussion:
//
//    This routine is a slightly modified version of the EISPACK routine to
//    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
//
//    The authors thank the authors of EISPACK for permission to use this
//    routine.
//
//    It has been modified to produce the product Q' * Z, where Z is an input
//    vector and Q is the orthogonal matrix diagonalizing the input matrix.
//    The changes consist (essentialy) of applying the orthogonal
//    transformations directly to Z as they are generated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//    Roger Martin, James Wilkinson,
//    The Implicit QL Algorithm,
//    Numerische Mathematik,
//    Volume 12, Number 5, December 1968, pages 377-383.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double D(N), the diagonal entries of the matrix.
//    On output, the information in D has been overwritten.
//
//    Input/output, double E(N), the subdiagonal entries of the
//    matrix, in entries E(1) through E(N-1).  On output, the information in
//    E has been overwritten.
//
//    Input/output, double Z(N).  On input, a vector.  On output,
//    the value of Q' * Z, where Q is the matrix that diagonalizes the
//    input symmetric tridiagonal matrix.
//
{
    double b;
    double c;
    double f;
    double g;
    int i;
    int ii;
    int itn = 30;
    int j;
    int k;
    int l;
    int m;
    int mml;
    double p;
    double prec;
    double r;
    double s;

    prec = r8_epsilon();

    if (n == 1) {
        return;
    }

    e[n - 1] = 0.0;

    for (l = 1; l <= n; l++) {
        j = 0;
        for (;;) {
            for (m = l; m <= n; m++) {
                if (m == n) {
                    break;
                }

                if (r8_abs(e[m - 1])
                    <= prec * (r8_abs(d[m - 1]) + r8_abs(d[m]))) {
                    break;
                }
            }
            p = d[l - 1];
            if (m == l) {
                break;
            }
            if (itn <= j) {
                cout << "\n";
                cout << "IMTQLX - Fatal error!\n";
                cout << "  Iteration limit exceeded\n";
                exit(1);
            }
            j = j + 1;
            g = (d[l] - p) / (2.0 * e[l - 1]);
            r = sqrt(g * g + 1.0);
            g = d[m - 1] - p + e[l - 1] / (g + r8_abs(r) * r8_sign(g));
            s = 1.0;
            c = 1.0;
            p = 0.0;
            mml = m - l;

            for (ii = 1; ii <= mml; ii++) {
                i = m - ii;
                f = s * e[i - 1];
                b = c * e[i - 1];

                if (r8_abs(g) <= r8_abs(f)) {
                    c = g / f;
                    r = sqrt(c * c + 1.0);
                    e[i] = f * r;
                    s = 1.0 / r;
                    c = c * s;
                } else {
                    s = f / g;
                    r = sqrt(s * s + 1.0);
                    e[i] = g * r;
                    c = 1.0 / r;
                    s = s * c;
                }
                g = d[i] - p;
                r = (d[i - 1] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i] = g + p;
                g = c * r - b;
                f = z[i];
                z[i] = s * z[i - 1] + c * f;
                z[i - 1] = c * z[i - 1] - s * f;
            }
            d[l - 1] = d[l - 1] - p;
            e[l - 1] = g;
            e[m - 1] = 0.0;
        }
    }
    //
    //  Sorting.
    //
    for (ii = 2; ii <= m; ii++) {
        i = ii - 1;
        k = i;
        p = d[i - 1];

        for (j = ii; j <= n; j++) {
            if (d[j - 1] < p) {
                k = j;
                p = d[j - 1];
            }
        }

        if (k != i) {
            d[k - 1] = d[i - 1];
            d[i - 1] = p;
            p = z[i - 1];
            z[i - 1] = z[k - 1];
            z[k - 1] = p;
        }
    }
    return;
}
//****************************************************************************80

void parchk(int kind, int m, double alpha, double beta)

//****************************************************************************80
//
//  Purpose:
//
//    PARCHK checks parameters ALPHA and BETA for classical weight functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
//
//    Input, int M, the order of the highest moment to
//    be calculated.  This value is only needed when KIND = 8.
//
//    Input, double ALPHA, BETA, the parameters, if required
//    by the value of KIND.
//
{
    double tmp;

    if (kind <= 0) {
        cout << "\n";
        cout << "PARCHK - Fatal error!\n";
        cout << "  KIND <= 0.\n";
        exit(1);
    }
    //
    //  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
    //
    if (3 <= kind && kind <= 8 && alpha <= -1.0) {
        cout << "\n";
        cout << "PARCHK - Fatal error!\n";
        cout << "  3 <= KIND and ALPHA <= -1.\n";
        exit(1);
    }
    //
    //  Check BETA for Jacobi.
    //
    if (kind == 4 && beta <= -1.0) {
        cout << "\n";
        cout << "PARCHK - Fatal error!\n";
        cout << "  KIND == 4 and BETA <= -1.0.\n";
        exit(1);
    }
    //
    //  Check ALPHA and BETA for rational.
    //
    if (kind == 8) {
        tmp = alpha + beta + m + 1.0;
        if (0.0 <= tmp || tmp <= beta) {
            cout << "\n";
            cout << "PARCHK - Fatal error!\n";
            cout << "  KIND == 8 but condition on ALPHA and BETA fails.\n";
            exit(1);
        }
    }
    return;
}
//****************************************************************************80

double r8_abs(double x)

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
    double value;

    if (0.0 <= x) {
        value = x;
    } else {
        value = -x;
    }
    return value;
}
//****************************************************************************80

double r8_epsilon()

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
    double value;

    value = 1.0;

    while (1.0 < (double)(1.0 + value)) {
        value = value / 2.0;
    }

    value = 2.0 * value;

    return value;
}
//****************************************************************************80

double r8_sign(double x)

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
    double value;

    if (x < 0.0) {
        value = -1.0;
    } else {
        value = 1.0;
    }
    return value;
}
//****************************************************************************80

void r8mat_write(string output_filename, int m, int n, double table[])

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file with no header.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the table data.
//
{
    int i;
    int j;
    ofstream output;
    //
    //  Open the file.
    //
    output.open(output_filename.c_str());

    if (!output) {
        cerr << "\n";
        cerr << "R8MAT_WRITE - Fatal error!\n";
        cerr << "  Could not open the output file.\n";
        return;
    }
    //
    //  Write the data.
    //
    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            output << "  " << setw(24) << setprecision(16) << table[i + j * m];
        }
        output << "\n";
    }
    //
    //  Close the file.
    //
    output.close();

    return;
}
//****************************************************************************80

void scqf(
    int nt, double t[], int mlt[], double wts[], int nwts, int ndx[],
    double swts[], double st[], int kind, double alpha, double beta, double a,
    double b)

//****************************************************************************80
//
//  Purpose:
//
//    SCQF scales a quadrature formula to a nonstandard interval.
//
//  Discussion:
//
//    The arrays WTS and SWTS may coincide.
//
//    The arrays T and ST may coincide.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, double T[NT], the original knots.
//
//    Input, int MLT[NT], the multiplicity of the knots.
//
//    Input, double WTS[NWTS], the weights.
//
//    Input, int NWTS, the number of weights.
//
//    Input, int NDX[NT], used to index the array WTS.
//    For more details see the comments in CAWIQ.
//
//    Output, double SWTS[NWTS], the scaled weights.
//
//    Output, double ST[NT], the scaled knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Input, double A, B, the interval endpoints.
//
{
    double al = 0.;
    double be = 0.;
    int i;
    int k;
    int l;
    double p;
    double shft = 0;
    double slp = 0;
    double temp;
    double tmp;

    temp = r8_epsilon();

    parchk(kind, 1, alpha, beta);

    if (kind == 1) {
        al = 0.0;
        be = 0.0;
        if (r8_abs(b - a) <= temp) {
            cout << "\n";
            cout << "SCQF - Fatal error!\n";
            cout << "  |B - A| too small.\n";
            cout << "kind =" << kind << endl;
            cout << "A =" << a << ", B =" << b << endl;
            exit(1);
        }
        shft = (a + b) / 2.0;
        slp = (b - a) / 2.0;
    } else if (kind == 2) {
        al = -0.5;
        be = -0.5;
        if (r8_abs(b - a) <= temp) {
            cout << "\n";
            cout << "SCQF - Fatal error!\n";
            cout << "  |B - A| too small.\n";
            cout << "kind =" << kind << endl;
            cout << "A =" << a << ", B =" << b << endl;
            exit(1);
        }
        shft = (a + b) / 2.0;
        slp = (b - a) / 2.0;
    } else if (kind == 3) {
        al = alpha;
        be = alpha;
        if (r8_abs(b - a) <= temp) {
            cout << "\n";
            cout << "SCQF - Fatal error!\n";
            cout << "  |B - A| too small.\n";
            cout << "kind =" << kind << endl;
            cout << "A =" << a << ", B =" << b << endl;
            exit(1);
        }
        shft = (a + b) / 2.0;
        slp = (b - a) / 2.0;
    } else if (kind == 4) {
        al = alpha;
        be = beta;

        if (r8_abs(b - a) <= temp) {
            cout << "\n";
            cout << "SCQF - Fatal error!\n";
            cout << "  |B - A| too small.\n";
            cout << "kind =" << kind << endl;
            cout << "A =" << a << ", B =" << b << endl;
            exit(1);
        }
        shft = (a + b) / 2.0;
        slp = (b - a) / 2.0;
    } else if (kind == 5) {
        if (b <= 0.0) {
            cout << "\n";
            cout << "SCQF - Fatal error!\n";
            cout << "  B <= 0\n";
            cout << "kind =" << kind << endl;
            exit(1);
        }
        shft = a;
        slp = 1.0 / b;
        al = alpha;
        be = 0.0;
    } else if (kind == 6) {
        if (b <= 0.0) {
            cout << "\n";
            cout << "SCQF - Fatal error!\n";
            cout << "  B <= 0.\n";
            cout << "kind =" << kind << endl;
            exit(1);
        }
        shft = a;
        slp = 1.0 / sqrt(b);
        al = alpha;
        be = 0.0;
    } else if (kind == 7) {
        al = alpha;
        be = 0.0;
        if (r8_abs(b - a) <= temp) {
            cout << "\n";
            cout << "SCQF - Fatal error!\n";
            cout << "  |B - A| too small.\n";
            cout << "kind =" << kind << endl;
            cout << "A =" << a << ", B =" << b << endl;
            exit(1);
        }
        shft = (a + b) / 2.0;
        slp = (b - a) / 2.0;
    } else if (kind == 8) {
        if (a + b <= 0.0) {
            cout << "\n";
            cout << "SCQF - Fatal error!\n";
            cout << "kind =" << kind << endl;
            cout << "  A + B <= 0.\n";
            exit(1);
        }
        shft = a;
        slp = a + b;
        al = alpha;
        be = beta;
    } else if (kind == 9) {
        al = 0.5;
        be = 0.5;
        if (r8_abs(b - a) <= temp) {
            cout << "\n";
            cout << "SCQF - Fatal error!\n";
            cout << "  |B - A| too small.\n";
            cout << "kind =" << kind << endl;
            cout << "A =" << a << ", B =" << b << endl;
            exit(1);
        }
        shft = (a + b) / 2.0;
        slp = (b - a) / 2.0;
    }

    p = pow(slp, al + be + 1.0);

    for (k = 0; k < nt; k++) {
        st[k] = shft + slp * t[k];
        l = abs(ndx[k]);

        if (l != 0) {
            tmp = p;
            for (i = l - 1; i <= l - 1 + mlt[k] - 1; i++) {
                swts[i] = wts[i] * tmp;
                tmp = tmp * slp;
            }
        }
    }
    return;
}
//****************************************************************************80

void sgqf(
    int nt, double aj[], double bj[], double zemu, double t[], double wts[])

//****************************************************************************80
//
//  Purpose:
//
//    SGQF computes knots and weights of a Gauss Quadrature formula.
//
//  Discussion:
//
//    This routine computes all the knots and weights of a Gauss quadrature
//    formula with simple knots from the Jacobi matrix and the zero-th
//    moment of the weight function, using the Golub-Welsch technique.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, double AJ[NT], the diagonal of the Jacobi matrix.
//
//    Input/output, double BJ[NT], the subdiagonal of the Jacobi
//    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.
//
//    Input, double ZEMU, the zero-th moment of the weight function.
//
//    Output, double T[NT], the knots.
//
//    Output, double WTS[NT], the weights.
//
{
    int i;
    //
    //  Exit if the zero-th moment is not positive.
    //
    if (zemu <= 0.0) {
        cout << "\n";
        cout << "SGQF - Fatal error!\n";
        cout << "  ZEMU <= 0.\n";
        exit(1);
    }
    //
    //  Set up vectors for IMTQLX.
    //
    for (i = 0; i < nt; i++) {
        t[i] = aj[i];
    }
    wts[0] = sqrt(zemu);
    for (i = 1; i < nt; i++) {
        wts[i] = 0.0;
    }
    //
    //  Diagonalize the Jacobi matrix.
    //
    imtqlx(nt, t, bj, wts);

    for (i = 0; i < nt; i++) {
        wts[i] = wts[i] * wts[i];
    }

    return;
}
