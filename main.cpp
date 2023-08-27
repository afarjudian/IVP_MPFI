/* 
 * File:   main.cpp
 * Author: amin
 *
 * Created on August 26, 2023, 6:17 PM
 */

#include <cstdlib>
#include <ctime>
#include <chrono>

#include <algorithm>

#include "gmpxx.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"



#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/multiprecision/mpfi.hpp>

#include "ivp_mpfi.h"
#include "fields.h"

#include <iostream>

// mpfi_float::default_precision seems to set the decimal precision
// rather than binary
#define PREC 50

// prec for std::cout
#define STD_PREC 7

//using namespace std;

using namespace boost::numeric::ublas;
using namespace boost::multiprecision;
using namespace boost::math::constants;




////////////////////////////////////
// global variables


// for time keeping
std::clock_t clock_start, clock_end;
auto chrono_start = std::chrono::high_resolution_clock::now();
auto chrono_end = std::chrono::high_resolution_clock::now();
double cpuTime;

// run the experiments from a minimum depth of (usually) 1 to the maxDepth
// int maxDepth = 16;
int maxDepth = 10;
int maxPicardIter = 15;



////////////////////////////////

/*
 * y' = y
 * y(0) = 1
 * solution y(t) = exp(t)
 */

void identityEulerFirstOrder(void) {

    /*
     * setting things up for a simple y'=y equation over [0,1], with y(0) = 1
     */

    vector< mpfi_float> y0(1);
    y0(0) = 1;

    vector< mpfi_float> b(1);
    mpfi_t bound;
    mpfi_init(bound);
    mpfi_interv_ui(bound, 0, 4);

    b(0) = bound;

    mpfi_float domain(0, 1);


    mpfi_float x(0.4, 0.7);

    //     M has to be chosen carefully
    mpfi_float M_exponential(3);
    mpfi_float M1_exponential(1);
    mpfi_float M2_exponential(0);


    ////////////////////////////////
    // Euler First Order


    std::cout << "==================================== " <<
            "\n y' = y, y(0) = 1" <<
            "\n solution y(t) = exp(t)" <<
            "\n \n First-order Euler" << std::endl;



    // The first column shows the depth, the second the maxWidth, and the third
    // contains the time (CPU time in ms)
    matrix< double> results(maxDepth, 3);

    matrix< mpfi_float> resultEulerFO;

    for (int depth = 1; depth <= maxDepth; depth++) {

        results(depth - 1, 0) = depth;

        clock_start = std::clock();

        resultEulerFO = euler_FOC_MPFI(
                identity,
                y0,
                domain,
                M_exponential,
                depth
                );

        clock_end = std::clock();
        cpuTime = 1000.0 * (clock_end - clock_start) / CLOCKS_PER_SEC;

        results(depth - 1, 2) = cpuTime;
        results(depth - 1, 1) = mpfr_get_d(maxWidth(resultEulerFO).backend().data(), MPFR_RNDU);
    }


    std::cout << std::fixed << std::setprecision(STD_PREC);
    std::cout << "==========================\n";
    std::cout << "depth\t maxWidth\t cpuTime\n\n";

    for (int depth = 1; depth <= maxDepth; depth++) {
        std::cout << results(depth - 1, 0) << "\t";
        std::cout << results(depth - 1, 1) << "\t";
        std::cout << results(depth - 1, 2) << "\n";
    }

    std::cout << "==========================\n";


    size_t nEulerFO = resultEulerFO.size1();
    size_t numOfIntervals_EulerFO = resultEulerFO.size2();


    std::cout << "y(End) = " <<
            resultEulerFO(0, numOfIntervals_EulerFO - 1) << std::endl;
    std::cout << "maxWidth = " << maxWidth(resultEulerFO) << std::endl;

    std::cout << "nEulerFO = " << nEulerFO << std::endl;
    std::cout << "numOfIntervals_EulerFO = " << numOfIntervals_EulerFO << std::endl;


    vector< mpfi_float> yEulerFO_x = evalAt(resultEulerFO, domain, x);
    std::cout << "x = " << x << std::endl;
    std::cout << "yEulerFO_x = " << yEulerFO_x << std::endl;
    std::cout << "========================\n";

}

void identityEulerSecondOrder(void) {
    /*
     * setting things up for a simple y'=y equation over [0,1], with y(0) = 1
     */

    vector< mpfi_float> y0(1);
    y0(0) = 1;

    vector< mpfi_float> b(1);
    mpfi_t bound;
    mpfi_init(bound);
    mpfi_interv_ui(bound, 0, 4);

    b(0) = bound;

    mpfi_float domain(0, 1);


    mpfi_float x(0.4, 0.7);

    //     M has to be chosen carefully
    mpfi_float M_exponential(3);
    mpfi_float M1_exponential(1);
    mpfi_float M2_exponential(0);

    ///////////////////////////////////////////////////////
    // second order Euler    
    std::cout << "==================================== " <<
            "\n y' = y, y(0) = 1" <<
            "\n solution y(t) = exp(t)" <<
            "\n \n Second-order Euler" << std::endl;



    // The first column shows the depth, the second the maxWidth, and the third
    // contains the time (CPU time in ms)
    matrix< double> results(maxDepth, 3);

    size_t nEulerSO;
    size_t numOfIntervals_EulerSO;

    matrix< mpfi_float> resultEulerSO;


    for (int depth = 1; depth <= maxDepth; depth++) {

        results(depth - 1, 0) = depth;

        clock_start = std::clock();

        resultEulerSO = euler_SOC_MPFI(
                identity,
                d_identity,
                y0,
                domain,
                M_exponential,
                depth
                );

        clock_end = std::clock();
        cpuTime = 1000.0 * (clock_end - clock_start) / CLOCKS_PER_SEC;

        results(depth - 1, 2) = cpuTime;
        results(depth - 1, 1) = mpfr_get_d(maxWidth(resultEulerSO).backend().data(), MPFR_RNDU);
    }

    std::cout << std::fixed << std::setprecision(STD_PREC);
    std::cout << "==========================\n";
    std::cout << "depth\t maxWidth\t cpuTime\n\n";

    for (int depth = 1; depth <= maxDepth; depth++) {
        std::cout << results(depth - 1, 0) << "\t";
        std::cout << results(depth - 1, 1) << "\t";
        std::cout << results(depth - 1, 2) << "\n";
    }


    std::cout << "==========================\n";


    nEulerSO = resultEulerSO.size1();
    numOfIntervals_EulerSO = resultEulerSO.size2();


    std::cout << "y(End) = " <<
            resultEulerSO(0, numOfIntervals_EulerSO - 1) << std::endl;
    std::cout << "maxWidth = " << maxWidth(resultEulerSO) << std::endl;

    std::cout << "nEulerSO = " << nEulerSO << std::endl;
    std::cout << "numOfIntervals_EulerSO = " << numOfIntervals_EulerSO << std::endl;


    vector< mpfi_float> yEulerSO_x = evalAt(resultEulerSO, domain, x);
    std::cout << "x = " << x << std::endl;
    std::cout << "yEulerSO_x = " << yEulerSO_x << std::endl;
    std::cout << "========================\n";


}

void identityEulerRungeKutta(void) {
    /*
     * setting things up for a simple y'=y equation over [0,1], with y(0) = 1
     */

    vector< mpfi_float> y0(1);
    y0(0) = 1;

    vector< mpfi_float> b(1);
    mpfi_t bound;
    mpfi_init(bound);
    mpfi_interv_ui(bound, 0, 4);

    b(0) = bound;

    mpfi_float domain(0, 1);


    mpfi_float x(0.4, 0.7);

    //     M has to be chosen carefully
    mpfi_float M_exponential(3);
    mpfi_float M1_exponential(1);
    mpfi_float M2_exponential(0);



    ////////////////////////////////
    // Euler Runge-Kutta


    std::cout << "==================================== " <<
            "\n y' = y, y(0) = 1" <<
            "\n solution y(t) = exp(t)" <<
            "\n \n Euler Runge-Kutta" << std::endl;



    // The first column shows the depth, the second the maxWidth, and the third
    // contains the time (CPU time in ms)
    matrix< double> results(maxDepth, 3);

    matrix< mpfi_float> resultEulerRK;

    for (int depth = 1; depth <= maxDepth; depth++) {

        results(depth - 1, 0) = depth;

        clock_start = std::clock();

        resultEulerRK = euler_RK_MPFI_Scalar(
                identity,
                d_identity,
                y0,
                domain,
                M_exponential,
                M1_exponential,
                M2_exponential,
                depth
                );


        clock_end = std::clock();
        cpuTime = 1000.0 * (clock_end - clock_start) / CLOCKS_PER_SEC;

        results(depth - 1, 2) = cpuTime;
        results(depth - 1, 1) = mpfr_get_d(maxWidth(resultEulerRK).backend().data(), MPFR_RNDU);
    }

    std::cout << std::fixed << std::setprecision(STD_PREC);
    std::cout << "==========================\n";
    std::cout << "depth\t maxWidth\t cpuTime\n\n";

    for (int depth = 1; depth <= maxDepth; depth++) {
        std::cout << results(depth - 1, 0) << "\t";
        std::cout << results(depth - 1, 1) << "\t";
        std::cout << results(depth - 1, 2) << "\n";


    }

    std::cout << "==========================\n";


    size_t nEulerRK = resultEulerRK.size1();
    size_t numOfIntervals_EulerRK = resultEulerRK.size2();


    std::cout << "y(End) = " <<
            resultEulerRK(0, numOfIntervals_EulerRK - 1) << std::endl;
    std::cout << "maxWidth = " << maxWidth(resultEulerRK) << std::endl;

    std::cout << "nEulerRK = " << nEulerRK << std::endl;
    std::cout << "numOfIntervals_EulerRK = " << numOfIntervals_EulerRK << std::endl;


    vector< mpfi_float> yEulerRK_x = evalAt(resultEulerRK, domain, x);
    std::cout << "x = " << x << std::endl;
    std::cout << "yEulerRK_x = " << yEulerRK_x << std::endl;
    std::cout << "========================\n";

}



///////////////////////////////

/*
 * the equation y'(t) = 10 cos(10 t) y(t), y(0) = 1
 * with the solution y(t) = exp( sin( 10 t))
 * We take  y1(t) = t, y2( t):= exp(sin(t))
 * thus [y1(0),y2(0)] = [0,1]
 * 
 * domain = [0,0.1]
 * initial bounds are [0,0.1]x[0,3]
 * 
 */



void tenCosTenTEulerFirstOrder(void) {

    vector< mpfi_float> y0(2);
    y0(0) = 0;
    y0(1) = 1;

    vector< mpfi_float> b(2);
    mpfi_t bound;
    mpfi_init(bound);
    mpfi_interv_d(bound, 0, 0.1);

    b(0) = bound;

    mpfi_interv_ui(bound, 0, 3);
    b(1) = bound;



    mpfi_float domain(0, 0.1);
    //
    //    mpfi_float x;
    //
    //
    mpfi_float M_expSin(30);


    //////////////////////////////////////////////////////
    // First-Order Euler


    std::cout << "==================================== " <<
            "\n y'(t) = 10 cos(10 t) y(t), y(0) = 1" <<
            "\n solution y(t) = exp( sin( 10 t))" <<
            "\n \n First-order Euler" << std::endl;


    //    The first column shows the depth, the second the maxWidth, and the third
    //    contains the time(CPU time in ms)
    matrix< double> results(maxDepth, 3);

    matrix< mpfi_float> resultEulerFO;


    for (int depth = 1; depth <= maxDepth; depth++) {

        results(depth - 1, 0) = depth;

        clock_start = std::clock();

        resultEulerFO = euler_FOC_MPFI(
                expSin,
                y0,
                domain,
                M_expSin,
                depth
                );

        clock_end = std::clock();
        cpuTime = 1000.0 * (clock_end - clock_start) / CLOCKS_PER_SEC;

        results(depth - 1, 2) = cpuTime;
        results(depth - 1, 1) = mpfr_get_d(maxWidth(resultEulerFO).backend().data(), MPFR_RNDU);
    }


    std::cout << std::fixed << std::setprecision(STD_PREC);
    std::cout << "==========================\n";
    std::cout << "depth\t maxWidth\t cpuTime\n\n";

    for (int depth = 1; depth <= maxDepth; depth++) {
        std::cout << results(depth - 1, 0) << "\t";
        std::cout << results(depth - 1, 1) << "\t";
        std::cout << results(depth - 1, 2) << "\n";
    }

    std::cout << "==========================\n";

}

void tenCosTenTEulerSecondOrder(void) {

    vector< mpfi_float> y0(2);
    y0(0) = 0;
    y0(1) = 1;

    vector< mpfi_float> b(2);
    mpfi_t bound;
    mpfi_init(bound);
    mpfi_interv_d(bound, 0, 0.1);

    b(0) = bound;

    mpfi_interv_ui(bound, 0, 3);
    b(1) = bound;



    mpfi_float domain(0, 0.1);
    //
    //    mpfi_float x;
    //
    //
    mpfi_float M_expSin(30);


    ////////////////////
    // Second-order Euler for y'(t) = 10 cos(10 t) y(t), y(0) = 1


    std::cout << "==================================== " <<
            "\n y'(t) = 10 cos(10 t) y(t), y(0) = 1" <<
            "\n solution y(t) = exp( sin( 10 t))" <<
            "\n \n Second-order Euler" << std::endl;



    // The first column shows the depth, the second the maxWidth, and the third
    // contains the time (CPU time in ms)
    matrix< double> results(maxDepth, 3);

    //        size_t nEulerSO;
    //        size_t numOfIntervals_EulerSO;


    matrix< mpfi_float> resultEulerSO;

    for (int depth = 1; depth <= maxDepth; depth++) {

        results(depth - 1, 0) = depth;

        clock_start = std::clock();

        resultEulerSO = euler_SOC_MPFI(
                expSin,
                d_expSin,
                y0,
                domain,
                M_expSin,
                depth);

        clock_end = std::clock();
        cpuTime = 1000.0 * (clock_end - clock_start) / CLOCKS_PER_SEC;

        results(depth - 1, 2) = cpuTime;
        results(depth - 1, 1) = mpfr_get_d(maxWidth(resultEulerSO).backend().data(), MPFR_RNDU);
    }


    std::cout << std::fixed << std::setprecision(STD_PREC);
    std::cout << "==========================\n";
    std::cout << "depth\t maxWidth\t cpuTime\n\n";

    for (int depth = 1; depth <= maxDepth; depth++) {
        std::cout << results(depth - 1, 0) << "\t";
        std::cout << results(depth - 1, 1) << "\t";
        std::cout << results(depth - 1, 2) << "\n";
    }

    mpfi_float x(0.05);

    vector< mpfi_float> yEulerSO_x = evalAt(resultEulerSO, domain, x);
    std::cout << "x = " << x << std::endl;
    std::cout << "yEulerSO_x = " << yEulerSO_x << std::endl;

    std::cout << "==========================\n";

}



//////////////////////////////////////////
// y'(t) =  | sin( t+y(t))|, y(0) = 1 

void absSinEulerFirstOrder(void) {
    vector< mpfi_float> y0(2);
    y0(0) = 0;
    y0(1) = 1;



    vector< mpfi_float> b(2);
    mpfi_t bound;
    mpfi_init(bound);
    mpfi_interv_ui(bound, 0, 5);

    b(0) = bound;


    mpfi_interv_ui(bound, -5, 5);
    b(1) = bound;



    mpfi_float domain(0, 5);

    mpfi_float x;


    mpfi_float M_absSin(2);



    ////////////////////////////////
    // Euler First Order


    std::cout << "==================================== " <<
            "\n y'(t) =  | sin( t+y(t))|, y(0) = 1 " <<
            "\n \n First-order Euler" << std::endl;


    // The first column shows the depth, the second the maxWidth, and the third
    // contains the time (CPU time in ms)
    matrix< double> results(maxDepth, 3);

    matrix< mpfi_float> resultEulerFO;


    for (int depth = 1; depth <= maxDepth; depth++) {

        results(depth - 1, 0) = depth;

        clock_start = std::clock();

        resultEulerFO = euler_FOC_MPFI(
                absSin,
                y0,
                domain,
                M_absSin,
                depth
                );



        clock_end = std::clock();
        cpuTime = 1000.0 * (clock_end - clock_start) / CLOCKS_PER_SEC;

        results(depth - 1, 2) = cpuTime;
        results(depth - 1, 1) = mpfr_get_d(maxWidth(resultEulerFO).backend().data(), MPFR_RNDU);
    }

    std::cout << std::fixed << std::setprecision(STD_PREC);
    std::cout << "==========================\n";
    std::cout << "depth\t maxWidth\t cpuTime\n\n";

    for (int depth = 1; depth <= maxDepth; depth++) {
        std::cout << results(depth - 1, 0) << "\t";
        std::cout << results(depth - 1, 1) << "\t";
        std::cout << results(depth - 1, 2) << "\n";
    }

    std::cout << "==========================\n";


}

void absSinEulerSecondOrder(void) {
    vector< mpfi_float> y0(2);
    y0(0) = 0;
    y0(1) = 1;



    vector< mpfi_float> b(2);
    mpfi_t bound;
    mpfi_init(bound);
    mpfi_interv_ui(bound, 0, 5);

    b(0) = bound;


    mpfi_interv_ui(bound, -5, 5);
    b(1) = bound;



    mpfi_float domain(0, 5);

    mpfi_float x;


    mpfi_float M_absSin(2);


    //    Second - order Euler for y'(t) =  | sin( t+y(t))|, y(0) = 1 
    std::cout << "==================================== " <<
            "\n y'(t) =  | sin( t+y(t))|, y(0) = 1 " <<
            "\n \n Second-order Euler" << std::endl;



    // The first column shows the depth, the second the maxWidth, and the third
    // contains the time (CPU time in ms)
    matrix< double> results(maxDepth, 3);

    size_t nEulerSO;
    size_t numOfIntervals_EulerSO;


    for (int depth = 1; depth <= maxDepth; depth++) {

        results(depth - 1, 0) = depth;

        clock_start = std::clock();

        matrix< mpfi_float> resultEulerSO = euler_SOC_MPFI(
                absSin,
                d_absSin,
                y0,
                domain,
                M_absSin,
                depth);

        clock_end = std::clock();
        cpuTime = 1000.0 * (clock_end - clock_start) / CLOCKS_PER_SEC;

        results(depth - 1, 2) = cpuTime;
        results(depth - 1, 1) = mpfr_get_d(maxWidth(resultEulerSO).backend().data(), MPFR_RNDU);
    }

    std::cout << std::fixed << std::setprecision(STD_PREC);
    std::cout << "==========================\n";
    std::cout << "depth\t maxWidth\t cpuTime\n\n";

    for (int depth = 1; depth <= maxDepth; depth++) {
        std::cout << results(depth - 1, 0) << "\t";
        std::cout << results(depth - 1, 1) << "\t";
        std::cout << results(depth - 1, 2) << "\n";
    }

    std::cout << "==========================\n";
}

/*
 * y'(t)=cos( y(t)) 
 * y(0) = 0
 * solution y(t) = 2 atan(  tanh(t/2 )) 
 */


void cosEulerRungeKutta(void) {

    vector< mpfi_float> y0(1);
    y0(0) = 0;

    vector< mpfi_float> b(1);
    mpfi_t bound;
    mpfi_init(bound);
    mpfi_interv_ui(bound, 0, 2);

    b(0) = bound;

    mpfi_float domain(0, 5);


    mpfi_float x(0.4, 0.7);

    //     M has to be chosen carefully

    // All of the following should really be 1
    // But we leave some room for inaccuracies
    mpfi_float M_cos(1);
    mpfi_float M1_cos(1);
    mpfi_float M2_cos(1);




    std::cout << "==================================== " <<
            "\n y'(t)=cos( y(t)), y(0) = 0" <<
            "\n solution y(t) = 2 atan(  tanh(t/2 )) " <<
            "\n \n Euler Runge-Kutta" << std::endl;



    // The first column shows the depth, the second the maxWidth, and the third
    // contains the time (CPU time in ms)
    matrix< double> results(maxDepth, 3);

    size_t nEulerRK;
    size_t numOfIntervals_EulerRK;

    matrix< mpfi_float> resultEulerRK;


    for (int depth = 1; depth <= maxDepth; depth++) {

        results(depth - 1, 0) = depth;

        clock_start = std::clock();

        resultEulerRK = euler_RK_MPFI_Scalar(
                cos,
                d_cos,
                y0,
                domain,
                M_cos,
                M1_cos,
                M2_cos,
                depth
                );

        clock_end = std::clock();
        cpuTime = 1000.0 * (clock_end - clock_start) / CLOCKS_PER_SEC;

        results(depth - 1, 2) = cpuTime;
        results(depth - 1, 1) = mpfr_get_d(maxWidth(resultEulerRK).backend().data(), MPFR_RNDU);
    }


    std::cout << std::fixed << std::setprecision(STD_PREC);
    std::cout << "==========================\n";
    std::cout << "depth\t maxWidth\t cpuTime\n\n";

    for (int depth = 1; depth <= maxDepth; depth++) {
        std::cout << results(depth - 1, 0) << "\t";
        std::cout << results(depth - 1, 1) << "\t";
        std::cout << results(depth - 1, 2) << "\n";
    }

    std::cout << "==========================\n";


    nEulerRK = resultEulerRK.size1();
    numOfIntervals_EulerRK = resultEulerRK.size2();


    std::cout << "y(End) = " <<
            resultEulerRK(0, numOfIntervals_EulerRK - 1) << std::endl;
    std::cout << "maxWidth = " << maxWidth(resultEulerRK) << std::endl;

    std::cout << "nEulerRK = " << nEulerRK << std::endl;
    std::cout << "numOfIntervals_EulerRK = " << numOfIntervals_EulerRK << std::endl;


    vector< mpfi_float> yEulerRK_x = evalAt(resultEulerRK, domain, x);
    std::cout << "x = " << x << std::endl;
    std::cout << "yEulerRK_x = " << yEulerRK_x << std::endl;
    std::cout << "========================\n";

}

void cosEulerFisrtOrder(void) {

    vector< mpfi_float> y0(1);
    y0(0) = 0;

    vector< mpfi_float> b(1);
    mpfi_t bound;
    mpfi_init(bound);
    mpfi_interv_ui(bound, 0, 2);

    b(0) = bound;

    mpfi_float domain(0, 5);


    mpfi_float x(0.4, 0.7);

    //     M has to be chosen carefully

    // All of the following should really be 1
    // But we leave some room for inaccuracies
    mpfi_float M_cos(1);
    mpfi_float M1_cos(1);
    mpfi_float M2_cos(1);

    std::cout << "==================================== " <<
            "\n y'(t)=cos( y(t)), y(0) = 0" <<
            "\n solution y(t) = 2 atan(  tanh(t/2 )) " <<
            "\n \n First-order Euler" << std::endl;



    // The first column shows the depth, the second the maxWidth, and the third
    // contains the time (CPU time in ms)
    matrix< double> results(maxDepth, 3);

    size_t nEulerFO;
    size_t numOfIntervals_EulerFO;

    matrix< mpfi_float> resultEulerFO;


    for (int depth = 1; depth <= maxDepth; depth++) {

        results(depth - 1, 0) = depth;

        clock_start = std::clock();

        resultEulerFO = euler_FOC_MPFI(
                cos,
                y0,
                domain,
                M_cos,
                depth
                );

        clock_end = std::clock();
        cpuTime = 1000.0 * (clock_end - clock_start) / CLOCKS_PER_SEC;

        results(depth - 1, 2) = cpuTime;
        results(depth - 1, 1) = mpfr_get_d(maxWidth(resultEulerFO).backend().data(), MPFR_RNDU);
    }

    std::cout << std::fixed << std::setprecision(STD_PREC);
    std::cout << "==========================\n";
    std::cout << "depth\t maxWidth\t cpuTime\n\n";

    for (int depth = 1; depth <= maxDepth; depth++) {
        std::cout << results(depth - 1, 0) << "\t";
        std::cout << results(depth - 1, 1) << "\t";
        std::cout << results(depth - 1, 2) << "\n";

    }

    std::cout << "==========================\n";


    nEulerFO = resultEulerFO.size1();
    numOfIntervals_EulerFO = resultEulerFO.size2();


    std::cout << "y(End) = " <<
            resultEulerFO(0, numOfIntervals_EulerFO - 1) << std::endl;
    std::cout << "maxWidth = " << maxWidth(resultEulerFO) << std::endl;

    std::cout << "nEulerFO = " << nEulerFO << std::endl;
    std::cout << "numOfIntervals_EulerFO = " << numOfIntervals_EulerFO << std::endl;


    vector< mpfi_float> yEulerFO_x = evalAt(resultEulerFO, domain, x);
    std::cout << "x = " << x << std::endl;
    std::cout << "yEulerFO_x = " << yEulerFO_x << std::endl;
    std::cout << "========================\n";

}

void cosEulerSecondOrder(void) {

    vector< mpfi_float> y0(1);
    y0(0) = 0;

    vector< mpfi_float> b(1);
    mpfi_t bound;
    mpfi_init(bound);
    mpfi_interv_ui(bound, 0, 2);

    b(0) = bound;

    mpfi_float domain(0, 5);


    mpfi_float x(0.4, 0.7);

    //     M has to be chosen carefully

    // All of the following should really be 1
    // But we leave some room for inaccuracies
    mpfi_float M_cos(1);
    mpfi_float M1_cos(1);
    mpfi_float M2_cos(1);


    std::cout << "==================================== " <<
            "\n y'(t)=cos( y(t)), y(0) = 0" <<
            "\n solution y(t) = 2 atan(  tanh(t/2 )) " <<
            "\n \n Second-order Euler" << std::endl;



    // The first column shows the depth, the second the maxWidth, and the third
    // contains the time (CPU time in ms)
    matrix< double> results(maxDepth, 3);

    size_t nEulerSO;
    size_t numOfIntervals_EulerSO;

    matrix< mpfi_float> resultEulerSO;


    for (int depth = 1; depth <= maxDepth; depth++) {

        results(depth - 1, 0) = depth;

        clock_start = std::clock();

        resultEulerSO = euler_SOC_MPFI(
                cos,
                d_cos,
                y0,
                domain,
                M_cos,
                depth
                );

        clock_end = std::clock();
        cpuTime = 1000.0 * (clock_end - clock_start) / CLOCKS_PER_SEC;

        results(depth - 1, 2) = cpuTime;
        results(depth - 1, 1) = mpfr_get_d(maxWidth(resultEulerSO).backend().data(), MPFR_RNDU);
    }

    std::cout << std::fixed << std::setprecision(STD_PREC);
    std::cout << "==========================\n";
    std::cout << "depth\t maxWidth\t cpuTime\n\n";

    for (int depth = 1; depth <= maxDepth; depth++) {
        std::cout << results(depth - 1, 0) << "\t";
        std::cout << results(depth - 1, 1) << "\t";
        std::cout << results(depth - 1, 2) << "\n";
    }

    std::cout << "==========================\n";
    nEulerSO = resultEulerSO.size1();
    numOfIntervals_EulerSO = resultEulerSO.size2();


    std::cout << "y(End) = " <<
            resultEulerSO(0, numOfIntervals_EulerSO - 1) << std::endl;
    std::cout << "maxWidth = " << maxWidth(resultEulerSO) << std::endl;

    std::cout << "nEulerSO = " << nEulerSO << std::endl;
    std::cout << "numOfIntervals_EulerSO = " << numOfIntervals_EulerSO << std::endl;


    vector< mpfi_float> yEulerSO_x = evalAt(resultEulerSO, domain, x);
    std::cout << "x = " << x << std::endl;
    std::cout << "yEulerSO_x = " << yEulerSO_x << std::endl;
    std::cout << "========================\n";

}

/**
 * Just the main function for testing and running everything!
 */
int main(int argc, char** argv) {



    std::cout << std::setprecision(STD_PREC);

    mpfi_float::default_precision(PREC);



    ////////////////////////////////

    /*
     * y' = y
     * y(0) = 1
     * solution y(t) = exp(t)
     * 
     * 
     * Uncomment any of the following to run the corresponding experiment
     */

    //    identityEulerRungeKutta();
    //    identityEulerFirstOrder();
    //    identityEulerSecondOrder();



    /*
     * y'(t)=cos( y(t)) 
     * y(0) = 0
     * solution y(t) = 2 atan(  tanh(t/2 )) 
     * 
     *  Uncomment any of the following to run the corresponding experiment
     */

    //    cosEulerRungeKutta();
    //    cosEulerFisrtOrder();
    cosEulerSecondOrder();



    /*
     * the equation y'(t) = 10 cos(10 t) y(t), y(0) = 1
     * with the solution y(t) = exp( sin( 10 t))
     * We take  y1(t) = t, y2( t):= exp(sin(t))
     * thus [y1(0),y2(0)] = [0,1]
     * 
     * domain = [0,0.1]
     * initial bounds are [0,0.1]x[0,3]
     * 
     * Uncomment any of the following to run the corresponding experiment
     */

    //    tenCosTenTEulerFirstOrder();
    //    tenCosTenTEulerSecondOrder();



    /*
     * y'(t) =  | sin( t+y(t))|, y(0) = 1 
     * 
     * Uncomment any of the following to run the corresponding experiment
     */

    //    absSinEulerFirstOrder();
    //    absSinEulerSecondOrder();


    return EXIT_SUCCESS;
}




