/*
*	Copyright (c) 2014 Michael Schellenberger Costa
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
*/

/****************************************************************************************************/
/*									Definition of all macros used									*/
/****************************************************************************************************/
#pragma once
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/facilities/empty.hpp>


/****************************************************************************************************/
/*									Macro for vector initialization									*/
/****************************************************************************************************/
#define _INIT(x)	{x, 0.0, 0.0, 0.0, 0.0}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									Macros for calculation of nth RK term							*/
/****************************************************************************************************/
#define _GET_VARNAMES(r, data, i, elem) BOOST_PP_COMMA_IF( i ) BOOST_PP_CAT(data, elem)
#define _SET_RK1(r, data, elem)			BOOST_PP_CAT(data, elem) = (elem[0]);
#define _SET_RK2(r, data, elem)			BOOST_PP_CAT(data, elem) = (elem[0] + elem[1] * 0.5);
#define _SET_RK3(r, data, elem)			BOOST_PP_CAT(data, elem) = (elem[0] + elem[2] * 0.5);
#define _SET_RK4(r, data, elem)			BOOST_PP_CAT(data, elem) = (elem[0] + elem[3]);

#define _SWITCH(...) 																		\
		double BOOST_PP_SEQ_FOR_EACH_I(_GET_VARNAMES, var_, __VA_ARGS__);					\
		switch(N) {																			\
			default:																		\
				BOOST_PP_SEQ_FOR_EACH(_SET_RK1, var_, __VA_ARGS__)							\
				break;																		\
			case 2:																			\
				BOOST_PP_SEQ_FOR_EACH(_SET_RK2, var_, __VA_ARGS__)							\
				break;																		\
			case 3:																			\
				BOOST_PP_SEQ_FOR_EACH(_SET_RK3, var_, __VA_ARGS__)							\
				break;																		\
			case 4:																			\
				BOOST_PP_SEQ_FOR_EACH(_SET_RK4, var_, __VA_ARGS__)							\
				break;																		\
		}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									Macros for repeated strings										*/
/****************************************************************************************************/
#define _REPEAT1(z, num, elem) BOOST_PP_COMMA_IF(num) elem
#define _REPEAT(elem, num) 	   BOOST_PP_REPEAT(num , _REPEAT1 , elem)
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
