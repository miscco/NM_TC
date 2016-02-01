/*
 *	Copyright (c) 2015 University of Lübeck
 *
 *	Permission is hereby granted, free of charge, to any person obtaining a copy
 *	of this software and associated documentation files (the "Software"), to deal
 *	in the Software without restriction, including without limitation the rights
 *	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *	copies of the Software, and to permit persons to whom the Software is
 *	furnished to do so, subject to the following conditions:
 *
 *	The above copyright notice and this permission notice shall be included in
 *	all copies or substantial portions of the Software.
 *
 *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *	THE SOFTWARE.
 *
 *	AUTHORS:	Michael Schellenberger Costa: mschellenbergercosta@gmail.com
 *              Stefanie Gareis: gareis@inb.uni-luebeck.de
 */

/****************************************************************************************************/
/*                                       Random number streams                                      */
/****************************************************************************************************/
#pragma once
#include <random>

/****************************************************************************************************/
/*									Struct for normal distribution                                  */
/****************************************************************************************************/
struct random_stream_normal
{
    /* Random number engine: Mersenne-Twister */
    std::mt19937_64                     mt;
    /* Random number generator: Normal-distribution */
    std::normal_distribution<double>    norm_dist;

    /* Constructors */
    random_stream_normal(){}
    random_stream_normal(double mean, double stddev)
    : mt(rand()) , norm_dist(mean, stddev)
    {}

    /* Overwrites the function-call operator "( )" */
    double operator( )(void) {
        return norm_dist(mt);
    }
};
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/

/****************************************************************************************************/
/*									Struct for uniform int distribution                             */
/****************************************************************************************************/
struct random_stream_uniform_int
{
    /* Random number engine: Mersenne-Twister */
    std::mt19937_64                     mt;
    /* Random number generator: Uniform integer-distribution */
    std::uniform_int_distribution<>     uniform_dist;

    /* Constructors */
    random_stream_uniform_int(){}
    random_stream_uniform_int(double lower_bound, double upper_bound)
    : mt(rand()) , uniform_dist(lower_bound, upper_bound)
    {}

    /* Overwrites the function-call operator "( )" */
    double operator( )(void) {
        return uniform_dist(mt);
    }
};
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/
