/****************************************************************************************************/
/*									functions for data storage										*/
/****************************************************************************************************/
#pragma once
#include <iostream>
#include <vector>
#include "Cortical_Column.h"
using std::vector;

// saving the fluctuations of the populations
inline void get_data(int counter, Cortical_Column& C, double* Ve) {
		Ve 	[counter] = Col.Ve	[0];
}
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/


/****************************************************************************************************/
/*							function that creates a matlab data container							*/
/****************************************************************************************************/
// function to create a MATLAB data container
mxArray* SetMexArray(int N, int M) {
	mxArray* Array	= mxCreateDoubleMatrix(0, 0, mxREAL);
    mxSetM(Array, N);
    mxSetN(Array, M);
    mxSetData(Array, mxMalloc(sizeof(double)*M*N));
    return Array;
}
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/
