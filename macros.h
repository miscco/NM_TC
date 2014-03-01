/****************************************************************************************************/
/*						header that defines the macros used in the simulation						*/
/****************************************************************************************************/
#pragma once
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/facilities/empty.hpp>

// macro for the initialization of the vectors
#define _INIT(x)	{x, 0.0, 0.0, 0.0, 0.0}

// macros to get the respective variables for RK terms
#define _RK1(x) (x[0])
#define _RK2(x) (x[0] + x[1] * 0.5)
#define _RK3(x) (x[0] + x[2] * 0.5)
#define _RK4(x) (x[0] + x[3])

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

// macros to add the RK terms
#define _ADD_RK(elem) elem [0] += (elem [1] + elem [2] *2 + elem [3] *2 + elem [4])/6;

// macro for repeated entry
#define _REPEAT1(z, num, elem) BOOST_PP_COMMA_IF(num) elem
#define _REPEAT(elem, num) 	 BOOST_PP_REPEAT(num , _REPEAT1 , elem)

/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/
