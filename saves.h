/****************************************************************************************************/
/*									functions for data storage										*/
/****************************************************************************************************/
#pragma once
#include <iostream>
#include "Cortical_Column.h"
#include "Thalamic_Column.h"

// saving the fluctuations of the populations
inline void get_data(int counter, Cortical_Column& Cortex, Thalamic_Column& Thalamus, double* Ve, double* Vt) {
	Ve 	[counter] = Cortex.Ve	[0];
	Vt 	[counter] = Thalamus.Vt	[0];
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
