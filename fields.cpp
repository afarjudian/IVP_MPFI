
#include <boost/multiprecision/mpfi.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>



using namespace boost::multiprecision;
using namespace boost::numeric::ublas;

vector< mpfi_float> identity(vector< mpfi_float> y) {
    return y;
}

matrix< mpfi_float> d_identity(vector< mpfi_float> y) {
    matrix< mpfi_float> result(1, 1);
    result(0, 0) = 1;
    return result;
}

vector< mpfi_float> cos(vector< mpfi_float> y) {
    vector< mpfi_float> result(1);
    result(0) = cos(y(0));
    return result;
}

matrix< mpfi_float> d_cos(vector< mpfi_float> y) {
    matrix< mpfi_float> result(1, 1);
    result(0, 0) = -sin(y(0));
    return result;
}

vector< mpfi_float> sinEq(vector< mpfi_float> y) {

    // commenting out for efficiency
    //    if (y.size() != 2) {
    //        std::cerr << "y must have 2 elements." << std::endl;
    //        exit(EXIT_FAILURE);
    //    }

    vector< mpfi_float> result(2);
    result(0) = y(1);
    result(1) = -y(0);

    return result;
}

matrix< mpfi_float> d_sinEq(vector< mpfi_float> y) {
    matrix< mpfi_float> result(2, 2);

    result(0, 0) = 0;
    result(0, 1) = 1;
    result(1, 0) = -1;
    result(1, 1) = 0;

    return result;
}

/*
 * the IVP y'(t) =  cos( t) y(t), y(0) = 1 has the solution
 * y(t) = exp( sin(  t))
 * 
 * We assign y1(t) = t, y2( t):= exp(sin(t))
 */

vector< mpfi_float> expSin(vector< mpfi_float> y) {


    // commenting out for efficiency
    //    if (y.size() != 2) {
    //        std::cerr << "y must have 2 elements." << std::endl;
    //        exit(EXIT_FAILURE);
    //    }

    mpfi_float mpfiTEN(10);

    vector< mpfi_float> result(2);
    result(0) = 1;
    result(1) = mpfiTEN * cos(mpfiTEN * y(0)) * y(1);

    return result;

}

matrix< mpfi_float> d_expSin(vector< mpfi_float> y) {


    // commenting out for efficiency
    //    if (y.size() != 2) {
    //        std::cerr << "y must have 2 elements." << std::endl;
    //        exit(EXIT_FAILURE);
    //    }

    mpfi_float mpfiTEN(10);
    mpfi_float mpfiMinusHundred(-100);

    matrix< mpfi_float> result(2, 2);
    result(0,0) = 0;
    result(0,1) = 0;
    result(1,0) = mpfiMinusHundred * y( 1) * sin( mpfiTEN * y( 0));
    result(1,1) = mpfiTEN * cos( mpfiTEN * y( 0));

    return result;

}


vector< mpfi_float> absSin(vector< mpfi_float> y) {


    // commenting out for efficiency
    //    if (y.size() != 2) {
    //        std::cerr << "y must have 2 elements." << std::endl;
    //        exit(EXIT_FAILURE);
    //    }

    vector< mpfi_float> result(2);
    result(0) = 1;
    result(1) = abs(sin(y(0) + y(1)));

    return result;

}

matrix< mpfi_float> d_absSin(vector< mpfi_float> y) {


    // commenting out for efficiency
    //    if (y.size() != 2) {
    //        std::cerr << "y must have 2 elements." << std::endl;
    //        exit(EXIT_FAILURE);
    //    }

    /*
     * This is crucial, the \hat{L}-derivative of abs function 
     * relies on the following coefficient, which is either
     * 1, -1, or [-1,1]
     */
    mpfi_float coeff;

    /*
     * we should check whether sin( x+y) is positive, negative, or it
     * contains 0
     */
    mpfi_float sinSum = sin(y(0) + y(1));

    mpfr_float a = boost::multiprecision::lower(sinSum);
    mpfr_float b = boost::multiprecision::upper(sinSum);

    if (0 < a) {
        coeff = 1;
    } else if (b < 0) {
        coeff = -1;
    } else {
        coeff = mpfi_float(-1, 1);
    }


    matrix< mpfi_float> result(2, 2);
    result(0, 0) = 0;
    result(0, 1) = 0;
    result(1, 0) = coeff * cos(y(0) + y(1));
    result(1, 1) = coeff * cos(y(0) + y(1));

    return result;

}
