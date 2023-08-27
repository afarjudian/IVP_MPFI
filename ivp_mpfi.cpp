#include <cstdlib>

#include "gmpxx.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "ivp_mpfi.h"



#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/multiprecision/mpfi.hpp>


using namespace boost::multiprecision;
using namespace boost::numeric::ublas;

mpfi_float evalQuadMPFIAt(quadMPFI qm, mpfi_float x) {
    return (qm.a * x * x) + (qm.b * x) + qm.c;
}

mpfi_float evalQEnclosureMPFI(quadEnclosureMPFI qe, mpfi_float x) {
    return hull(evalQuadMPFIAt(qe.lowerEnclosure, x),
            evalQuadMPFIAt(qe.upperEnclosure, x));
}

quadMPFI subtract(quadMPFI qm1, quadMPFI qm2) {
    quadMPFI result;
    result.a = qm1.a - qm2.a;
    result.b = qm1.b - qm2.b;
    result.c = qm1.c - qm2.c;

    return result;
}

mpfi_float evalCubicMPFIAt(cubicMPFI cm, mpfi_float x) {
    return (cm.a * x * x * x) + (cm.b * x * x) + (cm.c * x) + cm.d;
}

mpfi_float evalCEnclosureMPFI(cubicEnclosureMPFI ce, mpfi_float x) {
    return hull(evalCubicMPFIAt(ce.lowerEnclosure, x),
            evalCubicMPFIAt(ce.upperEnclosure, x));
}

cubicMPFI subtract(cubicMPFI cm1, cubicMPFI cm2) {
    cubicMPFI result;
    result.a = cm1.a - cm2.a;
    result.b = cm1.b - cm2.b;
    result.c = cm1.c - cm2.c;
    result.c = cm1.d - cm2.d;

    return result;
}

mpfi_float symExpand(mpfi_float a, mpfi_float r) {

    mpfi_float aPr = a + r;
    mpfi_float aMr = a - r;
    return hull(aPr, aMr);
}

vector<mpfi_float> symExpand(vector<mpfi_float> v, mpfi_float r) {

    size_t n = v.size();
    vector< mpfi_float> result(n);

    for (size_t i = 0; i < n; i++) {
        result(i) = symExpand(v(i), r);
    }
    return result;
}

matrix< mpfi_float> integrate(matrix< mpfi_float> y,
        mpfi_float domain,
        vector<mpfi_float> y0) {


    size_t n = y.size1(); // dimension
    size_t numOfIntervals = y.size2();

    if (n != y0.size()) {
        std::cerr << "y0 must have " << n << " elements." << std::endl;
        exit(EXIT_FAILURE);
    }

    // step size
    mpfi_float h(width(domain) / numOfIntervals);

    //    std::cout << "step size: " << h << std::endl;


    // the Piecewise linear integral has one element more than the number of 
    // intervals
    matrix< mpfi_float> pwl(n, numOfIntervals + 1);

    for (size_t i = 0; i < n; i++) {

        pwl(i, 0) = y0(i);

        for (size_t j = 1; j < numOfIntervals + 1; j++) {
            pwl(i, j) = pwl(i, j - 1) + (h * y(i, j - 1));
        }
    }


    // Now, the piecewise constant enclosure of the pwlIntegral
    return pwl_to_pwc(pwl);
}

matrix< mpfi_float> applyField(
        vector<mpfi_float> (*f)(vector<mpfi_float>),
        matrix< mpfi_float> y) {

    size_t n = y.size1();
    size_t numOfIntervals = y.size2();

    matrix< mpfi_float> fy(n, numOfIntervals);

    for (int j = 0; j < numOfIntervals; j++) {

        matrix_column< matrix < mpfi_float >> col(y, j);

        vector< mpfi_float> fyj = f(col);

        for (int i = 0; i < n; i++) {
            fy(i, j) = fyj(i);
        }
    }

    return fy;
}

matrix< mpfi_float> bisect(matrix< mpfi_float> y) {

    size_t n = y.size1();
    size_t numOfIntervals = y.size2();

    matrix< mpfi_float> result(n, numOfIntervals * 2);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < numOfIntervals; j++) {
            result(i, 2 * j) = y(i, j);
            result(i, 2 * j + 1) = y(i, j);
        }
    }

    return result;

}

void setInitialBound(
        matrix< mpfi_float> * y,
        vector< mpfi_float> b) {

    size_t n = (*y).size1();
    size_t numOfIntervals = (*y).size2();

    if (n != b.size()) {
        std::cerr << "b must have " << n << " elements." << std::endl;
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < numOfIntervals; j++) {
            (*y)(i, j) = b(i);
        }
    }



}

mpfr_float maxWidthLastColumn(matrix< mpfi_float> y) {

    size_t n = y.size1();
    size_t numOfIntervals = y.size2();

    mpfr_float w = 0;

    for (size_t i = 0; i < n; i++) {
        w = max(width(y(i, numOfIntervals - 1)), w);
    }

    return w;
}

mpfr_float maxWidth(matrix< mpfi_float> y) {

    size_t n = y.size1();
    size_t numOfIntervals = y.size2();

    mpfr_float w = 0;

    // When converting non-autonomous IVPs to autonomous ones,
    // the first component is usually time, which may be ignored
    //    for (size_t i = 1; i < n; i++) { 
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < numOfIntervals; j++) {
            w = max(width(y(i, j)), w);
        }
    }

    return w;
}

matrix< mpfi_float> enforceBounds(matrix< mpfi_float> y,
        vector< mpfi_float> b) {


    size_t n = y.size1();
    size_t numOfIntervals = y.size2();

    if (n != b.size()) {
        std::cerr << "b must have " << n << " elements." << std::endl;
        exit(EXIT_FAILURE);
    }


    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < numOfIntervals; j++) {
            y(i, j) = boost::multiprecision::intersect(y(i, j), b(i));
        }
    }
    return y;
}

matrix< mpfi_float> pwl_to_pwc(matrix< mpfi_float> pwl) {

    size_t n = pwl.size1();
    size_t numOfIntervals = pwl.size2() - 1;

    matrix< mpfi_float> pwc(n, numOfIntervals);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < numOfIntervals; j++) {
            pwc(i, j) = hull(pwl(i, j), pwl(i, j + 1));
        }
    }

    return pwc;
}

matrix< mpfi_float> pwq_to_pwc(matrix< quadEnclosureMPFI> pwq,
        mpfi_float domain) {

    size_t n = pwq.size1();
    size_t numOfIntervals = pwq.size2() - 1;

    matrix< mpfi_float> pwc(n, numOfIntervals);


    // the left end point of the domain
    mpfr_float a = boost::multiprecision::lower(domain);

    mpfi_float aMPFI(a);

    // step size
    mpfi_float h(width(domain) / numOfIntervals);


    // variables to be used in the loop    
    mpfi_float subInterval, left, right;

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < numOfIntervals; j++) {
            left = aMPFI + (j * h);
            right = aMPFI + ((j + 1) * h);

            subInterval = hull(left, right);

            pwc(i, j) = evalQEnclosureMPFI(pwq(i, j), subInterval);

        }
    }

    return pwc;
}

matrix< mpfi_float> pwCubic_to_pwc(matrix< cubicEnclosureMPFI> pwCubic,
        mpfi_float domain) {

    size_t n = pwCubic.size1();
    size_t numOfIntervals = pwCubic.size2() - 1;

    matrix< mpfi_float> pwc(n, numOfIntervals);


    // the left end point of the domain
    mpfr_float a = boost::multiprecision::lower(domain);

    mpfi_float aMPFI(a);

    // step size
    mpfi_float h(width(domain) / numOfIntervals);


    // variables to be used in the loop    
    mpfi_float subInterval, left, right;

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < numOfIntervals; j++) {
            left = aMPFI + (j * h);
            right = aMPFI + ((j + 1) * h);

            subInterval = hull(left, right);

            pwc(i, j) = evalCEnclosureMPFI(pwCubic(i, j), subInterval);

        }
    }

    return pwc;
}

matrix< mpfi_float> euler_FOC_MPFI(
        vector< mpfi_float> (*f)(vector< mpfi_float> y),
        vector< mpfi_float> y0,
        mpfi_float domain,
        mpfi_float M, // field bound
        size_t depth,
        size_t prec// precision        
        ) {

    mpfi_float::default_precision(prec);

    size_t n = y0.size();
    size_t numOfIntervals = 1 << depth;

    // first obtain a piecewise linear enclosure
    // which must have one element more than numOfIntervals in each row
    matrix< mpfi_float> pwl(n, numOfIntervals + 1);

    // assign y0 to the first column
    for (size_t i = 0; i < n; i++) {
        pwl(i, 0) = y0(i);
    }

    // step size
    mpfi_float h(width(domain) / numOfIntervals);

    // This is the amount of the symmetric expansion within the Euler
    // operator. In the paper it is written as \Delta q_i M
    mpfi_float hM = h * M;

    // variables to be used in the loop
    vector< mpfi_float> prev(n);
    vector< mpfi_float> prevExpanded(n);
    vector< mpfi_float> fPrevExpanded(n);

    // here is the Euler operator. We need to proceed column-by-column
    for (size_t j = 1; j < numOfIntervals + 1; j++) {
        // obtain the result up to the previous interval
        for (size_t i = 0; i < n; i++) {
            prev(i) = pwl(i, j - 1);
        }

        prevExpanded = symExpand(prev, hM);
        fPrevExpanded = f(prevExpanded);

        for (size_t i = 0; i < n; i++) {
            pwl(i, j) = prev(i) + h * fPrevExpanded(i);
        }
    }

    return pwl_to_pwc(pwl);
}


// Note that we only consider the scalar case for Runge-Kutta Euler
// But to keep it in line with the other implementations, the algorithm
// deals with vectors and matrices.

matrix< mpfi_float> euler_RK_MPFI_Scalar(
        vector< mpfi_float> (*f)(vector< mpfi_float> y),
        matrix< mpfi_float> (*df)(vector< mpfi_float> y),
        vector< mpfi_float> y0,
        mpfi_float domain,
        mpfi_float M,
        mpfi_float M1,
        mpfi_float M2,
        size_t depth,
        size_t prec
        ) {


    mpfi_float::default_precision(prec);

    size_t n = y0.size(); // this should really be 1 for now (scalar)
    size_t numOfIntervals = 1 << depth;

    // first obtain a piecewise cubic (pwCubic) enclosure
    // which must have one element more than numOfIntervals in each row
    matrix< cubicEnclosureMPFI> pwCubic(n, numOfIntervals + 1);

    // assign y0 to the first column
    for (size_t j = 0; j < n; j++) {

        mpfr_float upperY0 = boost::multiprecision::upper(y0(j));
        mpfr_float lowerY0 = boost::multiprecision::lower(y0(j));

        mpfi_float upperY0MPFI(upperY0);
        mpfi_float lowerY0MPFI(lowerY0);


        pwCubic(j, 0).upperEnclosure.a = 0;
        pwCubic(j, 0).upperEnclosure.b = 0;
        pwCubic(j, 0).upperEnclosure.c = 0;
        pwCubic(j, 0).upperEnclosure.d = upperY0MPFI;
        //        pwCubic(j, 0).upperEnclosure.d = y0(j);

        pwCubic(j, 0).lowerEnclosure.a = 0;
        pwCubic(j, 0).lowerEnclosure.b = 0;
        pwCubic(j, 0).lowerEnclosure.c = 0;
        pwCubic(j, 0).lowerEnclosure.d = lowerY0MPFI;
        //        pwCubic(j, 0).lowerEnclosure.d = y0(j);
    }


    // the left end point of the domain
    mpfr_float a = boost::multiprecision::lower(domain);
    mpfi_float aMPFI(a);


    // step size
    mpfi_float h(width(domain) / numOfIntervals);

    /* The implementation follows the same ideas and notations as that of 
     * the temporal discretization paper
     */

    // some variables to be used in the loop
    mpfi_float leftEndPoint;

    vector< mpfi_float>
            yLeftEndPoint(n),
            fyLeftEndPoint(n),
            fdfyLeftEndPoint(n);

    matrix <mpfi_float> dfyLeftEndPoint(n, n);

    // The following is the alpha value which appears as the coefficient 
    // of h^3 in the error term. i.e. [-alpha, alpha] *h^3
    mpfi_float alpha = (M2 * M + M1 * M1) * M / 6;


    /*
     * following ENTCS notation of Prop 5.9, the following are rho, nu and
     * sigma, up and down, respectively.
     */

    mpfr_float ru, rd,
            nu, nd,
            su, sd;

    /*
     * MPFI variants
     */
    mpfi_float ruM, rdM,
            nuM, ndM,
            suM, sdM;


    for (size_t i = 1; i < numOfIntervals + 1; i++) {

        // this is q_{i-1}, ENTCS paper notation
        leftEndPoint = aMPFI + h * (i - 1);

        // y(q_{i-1}), ENTCS paper notation
        for (size_t j = 0; j < n; j++) {
            yLeftEndPoint(j) =
                    evalCEnclosureMPFI(pwCubic(j, i - 1), leftEndPoint);
        }

        // u (y(q_{i-1})), , ENTCS paper notation
        fyLeftEndPoint = f(yLeftEndPoint);
        // u' (y(q_{i-1})), , ENTCS paper notation
        dfyLeftEndPoint = df(yLeftEndPoint);

        fdfyLeftEndPoint = prod(dfyLeftEndPoint, fyLeftEndPoint);


        for (size_t j = 0; j < n; j++) {
            rd = boost::multiprecision::lower(yLeftEndPoint(j));
            rdM = mpfi_float(rd);
            ru = boost::multiprecision::upper(yLeftEndPoint(j));
            ruM = mpfi_float(ru);

            nd = boost::multiprecision::lower(fyLeftEndPoint(j));
            ndM = mpfi_float(nd);
            nu = boost::multiprecision::upper(fyLeftEndPoint(j));
            nuM = mpfi_float(nu);

            sd = boost::multiprecision::lower(fdfyLeftEndPoint(j));
            sdM = mpfi_float(sd);
            su = boost::multiprecision::upper(fdfyLeftEndPoint(j));
            suM = mpfi_float(su);


            // TODO: +++ Amin Up to here*** We need cubic enclosures


            // finally, setting the coefficients. 
            // Note that, the expressions must be first converted to 
            // ax^3 + bx^2 + cx + d format for the 
            // following to make sense


            pwCubic(j, i).lowerEnclosure.a = -alpha;

            pwCubic(j, i).lowerEnclosure.b =
                    sdM / 2
                    + 3 * alpha * leftEndPoint;

            pwCubic(j, i).lowerEnclosure.c =
                    ndM
                    - sdM * leftEndPoint
                    - 3 * alpha * leftEndPoint*leftEndPoint;

            pwCubic(j, i).lowerEnclosure.d =
                    rdM
                    - ndM * leftEndPoint
                    + sdM * leftEndPoint * leftEndPoint / 2
                    + alpha * leftEndPoint * leftEndPoint * leftEndPoint;


            pwCubic(j, i).upperEnclosure.a = alpha;

            pwCubic(j, i).upperEnclosure.b =
                    suM / 2
                    - 3 * alpha * leftEndPoint;

            pwCubic(j, i).upperEnclosure.c =
                    nuM
                    - suM * leftEndPoint
                    + 3 * alpha * leftEndPoint*leftEndPoint;

            pwCubic(j, i).upperEnclosure.d =
                    ruM
                    - nuM * leftEndPoint
                    + suM * leftEndPoint * leftEndPoint / 2
                    - alpha * leftEndPoint * leftEndPoint * leftEndPoint;



        } // end for j
    } // end for i


    return pwCubic_to_pwc(pwCubic, domain);
} // end of euler_RK_MPFI

matrix< mpfi_float> euler_SOC_MPFI(
        vector< mpfi_float> (*f)(vector< mpfi_float> y),
        matrix< mpfi_float> (*df)(vector< mpfi_float> y),
        vector< mpfi_float> y0,
        mpfi_float domain,
        mpfi_float M,
        size_t depth,
        size_t prec
        ) {


    mpfi_float::default_precision(prec);

    size_t n = y0.size();
    size_t numOfIntervals = 1 << depth;

    // first obtain a piecewise quadratic (pwq) enclosure
    // which must have one element more than numOfIntervals in each row
    matrix< quadEnclosureMPFI> pwq(n, numOfIntervals + 1);

    // assign y0 to the first column
    for (size_t j = 0; j < n; j++) {

        mpfr_float upperY0 = boost::multiprecision::upper(y0(j));
        mpfr_float lowerY0 = boost::multiprecision::lower(y0(j));

        mpfi_float upperY0MPFI(upperY0);
        mpfi_float lowerY0MPFI(lowerY0);


        pwq(j, 0).upperEnclosure.a = 0;
        pwq(j, 0).upperEnclosure.b = 0;
        pwq(j, 0).upperEnclosure.c = upperY0MPFI;

        pwq(j, 0).lowerEnclosure.a = 0;
        pwq(j, 0).lowerEnclosure.b = 0;
        pwq(j, 0).lowerEnclosure.c = lowerY0MPFI;
    }


    // the left end point of the domain
    mpfr_float a = boost::multiprecision::lower(domain);
    mpfi_float aMPFI(a);


    // step size
    mpfi_float h(width(domain) / numOfIntervals);

    // This is the amount of the symmetric expansion within the Euler
    // operator. In the paper it is written as \Delta q_i M
    mpfi_float hM = h * M;


    /* The implementation follows the same ideas and notations as that of 
     * Prop 5.9 of the ENTCS Paper
     */

    // some variables to be used in the loop
    mpfi_float leftEndPoint;

    vector< mpfi_float>
            yLeftEndPoint(n),
            fyLeftEndPoint(n),
            yLeftEndPointSymExpand(n),
            fSymExpand(n),
            fdfSymExpand(n);

    matrix <mpfi_float> dfSymExpand(n, n);

    /*
     * following ENTCS notation of Prop 5.9, the following are rho, nu and
     * sigma, up and down, respectively.
     */

    mpfr_float ru, rd,
            nu, nd,
            su, sd;

    /*
     * MPFI variants
     */
    mpfi_float ruM, rdM,
            nuM, ndM,
            suM, sdM;


    for (size_t i = 1; i < numOfIntervals + 1; i++) {

        // this is q_{i-1}, ENTCS paper notation
        leftEndPoint = aMPFI + h * (i - 1);

        // y(q_{i-1}), ENTCS paper notation
        for (size_t j = 0; j < n; j++) {
            yLeftEndPoint(j) = evalQEnclosureMPFI(pwq(j, i - 1), leftEndPoint);
        }

        // u (y(q_{i-1})), , ENTCS paper notation
        fyLeftEndPoint = f(yLeftEndPoint);


        yLeftEndPointSymExpand = symExpand(yLeftEndPoint, hM);
        // the following should be obvious (follow ENTCS)        
        fSymExpand = f(yLeftEndPointSymExpand);
        dfSymExpand = df(yLeftEndPointSymExpand);

        fdfSymExpand = prod(dfSymExpand, fSymExpand);


        for (size_t j = 0; j < n; j++) {
            rd = boost::multiprecision::lower(yLeftEndPoint(j));
            rdM = mpfi_float(rd);
            ru = boost::multiprecision::upper(yLeftEndPoint(j));
            ruM = mpfi_float(ru);

            nd = boost::multiprecision::lower(fyLeftEndPoint(j));
            ndM = mpfi_float(nd);
            nu = boost::multiprecision::upper(fyLeftEndPoint(j));
            nuM = mpfi_float(nu);

            sd = boost::multiprecision::lower(fdfSymExpand(j));
            sdM = mpfi_float(sd);
            su = boost::multiprecision::upper(fdfSymExpand(j));
            suM = mpfi_float(su);


            // finally, setting the coefficients. 
            // Note that, the expressions of the ENTCS paper (Prop 5.9)
            // must be first converted to ax^2 + bx + c format for the 
            // following to make sense

            pwq(j, i).lowerEnclosure.a = sdM / 2;
            pwq(j, i).lowerEnclosure.b = ndM - (sdM * leftEndPoint);
            pwq(j, i).lowerEnclosure.c =
                    (sdM * leftEndPoint * leftEndPoint) / 2
                    - (ndM * leftEndPoint)
                    + rdM;

            pwq(j, i).upperEnclosure.a = suM / 2;
            pwq(j, i).upperEnclosure.b = nuM - (suM * leftEndPoint);
            pwq(j, i).upperEnclosure.c =
                    (suM * leftEndPoint * leftEndPoint) / 2
                    - (nuM * leftEndPoint)
                    + ruM;

        } // end for j
    } // end for i


    return pwq_to_pwc(pwq, domain);
} // end of euler_RK_MPFI_Scalar



vector< mpfi_float> evalAt(matrix< mpfi_float> y,
        mpfi_float domain,
        mpfi_float x) {

    size_t n = y.size1();
    size_t numOfIntervals = y.size2();

    // the two end points of the domain
    mpfr_float a = boost::multiprecision::lower(domain);
    mpfr_float b = boost::multiprecision::upper(domain);

    // endpoints of x
    mpfr_float xL = boost::multiprecision::lower(x);
    mpfr_float xU = boost::multiprecision::upper(x);

    // making sure we stay within bound
    mpfr_float leftEnd = max(a, xL);
    mpfr_float rightEnd = min(b, xU);


    // width of each subinterval
    mpfr_float h = (b - a) / numOfIntervals;

    mpfr_float d1 = (leftEnd - a) / h;
    mpfr_float d2 = (rightEnd - a) / h;

    size_t leftIndex = mpfr_get_ui(d1.backend().data(), MPFR_RNDD);
    leftIndex = leftIndex > 0 ? leftIndex : 0;

    size_t rightIndex = mpfr_get_ui(d2.backend().data(), MPFR_RNDU);
    rightIndex = rightIndex < numOfIntervals - 1 ?
            rightIndex : numOfIntervals - 1;


    // create a vector of size n
    vector< mpfi_float> result(n);

    // first assign the interval at the leftIndex
    for (size_t i = 0; i < n; i++) {
        result(i) = y(i, leftIndex);
    }

    // now form the hull with the remaining, iteratively
    for (size_t j = leftIndex + 1; j <= rightIndex; j++) {
        for (size_t i = 0; i < n; i++) {
            result(i) = hull(result(i), y(i, j));
        }
    }

    return result;
}
