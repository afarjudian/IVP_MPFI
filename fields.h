/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   fields.h
 * Author: amin
 *
 * Created on September 6, 2020, 3:33 PM
 */

#ifndef FIELDS_H
#define FIELDS_H


#include <boost/multiprecision/mpfi.hpp>



using namespace boost::multiprecision;
using namespace boost::numeric::ublas;


/**
 * field for y'=y
 * @param y 1D vector
 * @return y
 */
vector< mpfi_float> identity( vector< mpfi_float> y);

/**
 * \hat{L}-derivative for the identity field, arising from the IVP y'=y
 * @param y 1D vector
 * @return 1x1 matrix with the only element being 1, i.e., [1]
 */
matrix< mpfi_float> d_identity(vector< mpfi_float> y);


/**
 * field for y'(t)=cos(y(t))
 * @param y 1D vector
 * @return cos(y)
 */
vector< mpfi_float> cos( vector< mpfi_float> y);

/**
 * \hat{L}-derivative for the cos field, arising from the IVP y'(t)=cos(y(t))
 * @param y 1D vector
 * @return 1x1 matrix with the only element being cos(y(0))
 */
matrix< mpfi_float> d_cos(vector< mpfi_float> y);



/**
 * field for y''=y. We assign y2:=y1'
 * 
 * @param y 2D vector [y1(t), y2(t)]
 * @return [y2(t), -y1(t)]
 */
vector< mpfi_float> sinEq( vector< mpfi_float> y);


/**
 * \hat{L}-derivative for the field of the IVP y''=y.
 * We assign y2:=y1' 
 * The field is f(x,y)=(y, -x). Hence df(x,y) is the matrix
 * [ 0 1;
 *  -1 0]
 * 
 * @param y 2D vector [y1(t), y2(t)]
 * @return the matrix [ 0 1; -1 0]
 */
matrix< mpfi_float> d_sinEq( vector< mpfi_float> y);



/**
 * the IVP y'(t) =  10 cos( 10 t) y(t), y(0) = 1 has the solution
 * y(t) = exp( sin(10 t))
 * 
 * We assign y1(t) = t, y2( t):= exp( sin( 10 t))
 * @param y 2D vector [y1(t), y2(t)]
 * @return  [ 1, 10 cos( 10 y1(t)) y2(t)]
 */

vector< mpfi_float> expSin(vector< mpfi_float> y);

/**
 * \hat{L}-derivative of the field for the IVP
 * y'(t) =  10 cos( 10 t) y(t), y(0) = 1, which has the solution
 * y(t) = exp( sin(10 t))
 * 
 * We assign y1(t) = t, y2( t):= exp( sin( 10 t))
 * 
 * 
 * @param y 2D vector [y1(t), y2(t)]
 * @return the matrix [ 0 0; r s] with r = -100 y2(t) sin( 10 y1(t)) and s = 10 cos( 10 y1(t))
 */
matrix< mpfi_float> d_expSin(vector< mpfi_float> y);



/**
 * the IVP y'(t) =  | sin( t+y(t))|, y(0) = 1 
 * with non-differentiable field
 * the solution is checked against an octave implementation, as it 
 * does not have a closed-form solution.
 * 
 * The Octave implementation is called absSin.m in 
 * Academic/Programming/Octave/Diff_EQ/Non_Differentiable_Field
 * 
 * We assign y1(t) = t, y2( t):= y(t)
 * 
 * @param y 2D vector [y1(t), y2(t)]
 * @return [1, |sin(y1(t) + y2(t))|]
 */
vector< mpfi_float> absSin(vector< mpfi_float> y);


/**
 * \hat{L}-derivative of the field for the IVP
 * y'(t) =  | sin( t+y(t))|, y(0) = 1 
 * with non-differentiable field
 * the solution is checked against an octave implementation, as it 
 * does not have a closed-form solution
 * The Octave implementation is called absSin.m in 
 * Academic/Programming/Octave/Diff_EQ/Non_Differentiable_Field
 * 
 * We assign y1(t) = t, y2( t):= y(t)
 * 
 * @param y 2D vector [y1(t), y2(t)]
 * @return the matrix [ 0 0; r s] with r = s = coeff * cos( y1(t) + y2(t)) and 
 * coeff is 1, -1, or [-1,1] depending on whether sin( y1(t) + y2(t)) is greater
 * than 0, less than 0, or contains 0 (respectively).
 */
matrix< mpfi_float> d_absSin(vector< mpfi_float> y);

#endif /* FIELDS_H */

