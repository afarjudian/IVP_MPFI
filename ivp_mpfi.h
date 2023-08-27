/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ivp_mpfi.h
 * Author: amin
 *
 * Created on September 6, 2020, 1:36 PM
 */

#ifndef IVP_MPFI_H
#define IVP_MPFI_H

#include <boost/multiprecision/mpfi.hpp>
#include <boost/numeric/ublas/vector.hpp>



using namespace boost::multiprecision;
using namespace boost::numeric::ublas;

/**
 * a structure representing the quadratic polynomial:
 *    a x^2 + b x + c
 */
typedef struct {
    mpfi_float a;
    mpfi_float b;
    mpfi_float c;
} quadMPFI;

/**
 * A structure representing a quadratic enclosure over an interval. 
 * The number of enclosures is (normally) 1+numOfIntervals
 * The reason is that the first enclosure only tells what the value is
 * at the left end-point of the domain.
 * 
 */
typedef struct {
    quadMPFI upperEnclosure;
    quadMPFI lowerEnclosure;

} quadEnclosureMPFI;

/**
 * Evaluates the quadratic polynomial
 *  a x^2 + b x + c
 * at interval x
 * 
 * @param qm quadratic polynomial
 * @param x value to be evaluated at
 * @return quadratic polynomial evaluated at x
 */
mpfi_float evalQuadMPFIAt(quadMPFI qm, mpfi_float x);

/**
 * Evaluates a given quadratic enclosure over an interval x
 * @param qe quadratic enclosure
 * @param x interval
 * @return qe( x)
 */
mpfi_float evalQEnclosureMPFI(quadEnclosureMPFI qe, mpfi_float x);

/**
 * given quadratic polynomials ( a1, b1, c1) and ( a2, b2, c2)
 * returns ( a1 - a2, b1 - b2, c1 - c2)
 * 
 * @param qm1 first argument
 * @param qm2 second argument
 * @return the result of subtraction
 */
quadMPFI subtract(quadMPFI qm1, quadMPFI qm2);

/**
 * a structure representing the cubic polynomial:
 *    a x^3 + b x^2 + c x + d
 */
typedef struct {
    mpfi_float a;
    mpfi_float b;
    mpfi_float c;
    mpfi_float d;
} cubicMPFI;

/**
 * A structure representing a cubic enclosure over an interval. 
 * The number of enclosures is (normally) 1+numOfIntervals
 * The reason is that the first enclosure only tells what the value is
 * at the left end-point of the domain.
 * 
 */
typedef struct {
    cubicMPFI upperEnclosure;
    cubicMPFI lowerEnclosure;

} cubicEnclosureMPFI;

/**
 * Evaluates the cubic polynomial
 *  a x^3 + b x^2 + c x + d
 * at interval x
 * 
 * @param cm cubic polynomial
 * @param x value to be evaluated at
 * @return cubic polynomial evaluated at x
 */
mpfi_float evalCubicMPFIAt(cubicMPFI cm, mpfi_float x);

/**
 * Evaluates a given cubic enclosure over an interval x
 * @param ce cubic enclosure
 * @param x interval
 * @return ce( x)
 */
mpfi_float evalCEnclosureMPFI(cubicEnclosureMPFI ce, mpfi_float x);

/**
 * given cubic polynomials ( a1, b1, c1, d1) and ( a2, b2, c2, d2)
 * returns ( a1 - a2, b1 - b2, c1 - c2, d1 - d2)
 * 
 * @param cm1 first argument
 * @param cm2 second argument
 * @return the result of subtraction
 */
cubicMPFI subtract(cubicMPFI cm1, cubicMPFI cm2);

/**
 * Symmetric expansion of an interval [x,y] with r >= 0 is
 * defined as [x-r, y+r]
 * @param a interval to be expanded
 * @param r amount of expansion 
 * @return SymExpand( a, r)
 */
mpfi_float symExpand(mpfi_float a, mpfi_float r);


/**
 * The vector version of symExpand
 * @param v vector to be symmetrically expanded
 * @param r amount of expansion, r >= 0
 * @return SymExpand( v, r)
 */
vector<mpfi_float> symExpand(vector<mpfi_float> v, mpfi_float r);

/**
 * @param domain = [a,b], and interval
 * @param y: an interval enclosure of a function [a,b] --> R^n
 * y has n rows, and the number of columns determines the number of intervals
 * 
 * @param y0 an n-tuple of intervals, acting as the additive constant
 * 
 * @return y0 + \int_[a,b] y dt
 * 
 */
matrix< mpfi_float> integrate(matrix< mpfi_float> y,
        mpfi_float domain,
        vector<mpfi_float> y0);



/**
 * applies the given field over the n-tuple of enclosures
 * @param f field
 * @param y n-tuple of enclosures
 * @return f(y)
 */

matrix< mpfi_float> applyField(
        vector<mpfi_float> (*f)(vector<mpfi_float>),
        matrix< mpfi_float> y);




/**
 * Given a set of intervals, bisects them.
 * @param y n-tuple of intervals
 * @return bisection of y
 */
matrix< mpfi_float> bisect(matrix< mpfi_float> y);




/**
 * Initializes y to the given vector b
 * @param y pointer to matrix of enclosures
 * @param b vector of initial bound
 */

void setInitialBound(
        matrix< mpfi_float> * y,
        vector< mpfi_float> b);



/**
 * returns the maximum width of the last column of y
 * @param y input matrix
 * @return maximum width of the last column of y
 */

mpfr_float maxWidthLastColumn(matrix< mpfi_float> y);

/**
 * returns the maximum width of y
 * @param y input matrix
 * @return maximum width of y
 */

mpfr_float maxWidth(matrix< mpfi_float> y);



/**
 * Returns a piecewise constant enclosure of a piecewise linear enclosure
 * @param pwl piecewise linear enclosure
 * @return pwc enclosure of pwl
 */
matrix< mpfi_float> pwl_to_pwc(matrix< mpfi_float> pwl);


/**
 * Given a piece-wise quadratic enclosure, returns a pwc enclosure 
 * 
 * @param pwq a matrix where (i,j)-th element is the quadratic enclosure of the 
 * i-th coordinate over the j-th interval
 * @param domain the domain as in interval
 * @return pwc piece-wise constant enclosure of the given pwq enclosure
 */
matrix< mpfi_float> pwq_to_pwc(matrix< quadEnclosureMPFI> pwq, mpfi_float domain);

/**
 * Given a piece-wise cubic enclosure, returns a piece-wise constant enclosure 
 * 
 * @param pwCubic a matrix where (i,j)-th element is the cubic enclosure of the 
 * i-th coordinate over the j-th interval
 * @param domain the domain as in interval
 * @return pwc piece-wise constant enclosure of the given pwCubic enclosure
 */
matrix< mpfi_float> pwCubic_to_pwc(matrix< cubicEnclosureMPFI> pwCubic, mpfi_float domain);

/**
 * Solves the IVP y'(t) = f(y(t)) using first-order interval Euler method
 * with constant expansion,
 * over domain = [a,b] subject to y(a)=y0 
 * 
 * field f: [-K,K]^n -> [-M,M]^n
 * solution y: [a,b] -> [-K,K]^n
 * with (b-a) <= (K / M) 
 * 
 * 
 * @param f field
 * @param y0 initial value
 * @param domain the interval [a,b] over which to solve the IVP
 * @param M bound on the range of the field f (should be chosen very carefully)
 * @param depth The number of intervals is 2^depth
 * @param prec precision of representation
 * @return a piecewise constant enclosure of the solution
 */
matrix< mpfi_float> euler_FOC_MPFI(
        vector< mpfi_float> (*f)(vector< mpfi_float> y),
        vector< mpfi_float> y0,
        mpfi_float domain,
        mpfi_float M, // field bound
        size_t depth = 10,
        size_t prec = 500 // precision        
        );


/**
 * Solves the IVP y'(t) = f(y(t)) using second-order domain theoretic 
 * Euler method, over domain = [a,b] subject to y(a)=y0 
 * 
 * field f: [-K,K]^n -> [-M,M]^n
 * \hat{L}-derivative of the field: df: [-K,K]^n -> [-M',M']^(nxn)
 * solution y: [a,b] -> [-K,K]^n
 * with (b-a) <= (K / (M(1+nM'))) 
 * 
 * 
 * @param f field
 * @param df \hat{L}-derivative of the field
 * @param y0 initial value
 * @param domain the interval [a,b] over which to solve the IVP
 * @param M bound on the range of the field f (should be chosen very carefully)
 * @param depth The number of intervals is 2^depth
 * @param prec precision of representation
 * @return a piecewise constant enclosure of the solution
 */
matrix< mpfi_float> euler_SOC_MPFI(
        vector< mpfi_float> (*f)(vector< mpfi_float> y),
        matrix< mpfi_float> (*df)(vector< mpfi_float> y),
        vector< mpfi_float> y0,
        mpfi_float domain,
        mpfi_float M, // field bound
        size_t depth = 10,
        size_t prec = 500// precision        
        );




// Note that we only consider the scalar case for Runge-Kutta Euler
// But to keep it in line with the other implementations, the algorithm
// deals with vectors and matrices.

/**
 * 
 * @param f field
 * @param df \hat{L}-derivative of the field
 * @param y0 initial value
 * @param domain the interval [a,b] over which to solve the IVP
 * @param M bound on the range of the field f (should be chosen very carefully)
 * @param M1 bound on the range of the derivative of field f (should be chosen very carefully)
 * @param M2 bound on the range of the second derivative of field f (should be chosen very carefully)
 * @param depth The number of intervals is 2^depth
 * @param prec precision of representation
 * @return a piecewise constant enclosure of the solution

 */

matrix< mpfi_float> euler_RK_MPFI_Scalar(
        vector< mpfi_float> (*f)(vector< mpfi_float> y),
        matrix< mpfi_float> (*df)(vector< mpfi_float> y),
        vector< mpfi_float> y0,
        mpfi_float domain,
        mpfi_float M,
        mpfi_float M1,
        mpfi_float M2,
        size_t depth = 10,
        size_t prec = 500// precision        
        );


/**
 * Evaluates a matrix of enclosures at a given interval
 * @param y matrix of enclosures
 * @param domain the interval [a,b] over which y is defined
 * @param x interval over which y is evaluates
 * @return a vector of n intervals, where n = y.size1()
 */
vector< mpfi_float> evalAt(matrix< mpfi_float> y,
        mpfi_float domain,
        mpfi_float x);


/**
 * For a given enclosure matrix y, intersects each element
 * with the appropriate component of the bound b
 * @param y enclosure matrix
 * @param b bound
 * @return intersection of y and b
 */
matrix< mpfi_float> enforceBounds(matrix< mpfi_float> y,
        vector< mpfi_float> b);

#endif /* IVP_MPFI_H */

