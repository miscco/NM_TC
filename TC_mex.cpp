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
 *
 *	Based on:	A thalamocortical neural mass model of the EEG during NREM sleep and its response
 *				to auditory stimulation.
 *				M Schellenberger Costa, A Weigenand, H-VV Ngo, L Marshall, J Born, T Martinetz,
 *				JC Claussen.
 *				PLoS Computational Biology http://dx.doi.org/10.1371/journal.pcbi.1005022
 */

/******************************************************************************/
/* Implementation of the simulation as MATLAB routine (mex compiler)		  */
/* mex command is given by:													  */
/* mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" TC_mex.cpp Cortical_Column.cpp    */
/*                                          Thalamic_Column.cpp               */
/******************************************************************************/
#include "mex.h"
#include "matrix.h"

#include <iterator>
#include <vector>

#include "Cortical_Column.h"
#include "Data_Storage.h"
#include "ODE.h"
#include "Stimulation.h"
#include "Thalamic_Column.h"
mxArray* GetMexArray(int N, int M);
mxArray* get_marker(Stim &stim);

/******************************************************************************/
/*                          Fixed simulation settings						  */
/******************************************************************************/
extern const int onset	= 20;		/* Time until data is stored in  s		  */
extern const int res 	= 1E4;		/* Number of iteration steps per s		  */
extern const int red 	= 1E2;		/* Number of iterations steps not saved	  */
extern const double dt 	= 1E3/res;	/* Duration of a time step in ms		  */
extern const double h	= sqrt(dt); /* Square root of dt for SRK iteration	  */

/******************************************************************************/
/*                              Simulation routine	 						  */
/*								lhs defines outputs							  */
/*								rhs defines inputs							  */
/******************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Initialize the seeder */
    srand(time(NULL));

    /* Fetch inputs */
    const int T				= (int) (mxGetScalar(prhs[0]));	/* Duration of simulation in s			*/
    const int Time 			= (T+onset)*res;				/* Total number of iteration steps		*/
    double* Param_Cortex	= mxGetPr (prhs[1]);			/* Parameters of cortical module		*/
    double* Param_Thalamus	= mxGetPr (prhs[2]);			/* Parameters of thalamic module		*/
    double* Connections		= mxGetPr (prhs[3]);			/* Connectivity values C <-> T			*/
    double* var_stim	 	= mxGetPr (prhs[4]);			/* Parameters of stimulation protocol	*/

    /* Initialize the populations */
    Cortical_Column Cortex	 = Cortical_Column(Param_Cortex,   Connections);
    Thalamic_Column Thalamus = Thalamic_Column(Param_Thalamus, Connections);

    /* Link both modules */
    Cortex.get_Thalamus(Thalamus);
    Thalamus.get_Cortex(Cortex);

    /* Initialize the stimulation protocol */
    Stim Stimulation(Cortex, Thalamus, var_stim);

    /* Create data containers */
    std::vector<mxArray*> dataArray;
    dataArray.reserve(4);
    dataArray.push_back(GetMexArray(1, T*res/red));	// Vt
    dataArray.push_back(GetMexArray(1, T*res/red));	// Vr
    dataArray.push_back(GetMexArray(1, T*res/red));	// Ca
    dataArray.push_back(GetMexArray(1, T*res/red));	// act_h

    /* Pointer to the data blocks */
    std::vector<double*> dataPointer;
    dataPointer.reserve(dataArray.size());
    for (auto &dataptr : dataArray) {
        dataPointer.push_back(mxGetPr(dataptr));
    }

    /* Simulation */
    int count = 0;
    for (unsigned t=0; t < Time; ++t) {
        ODE (Cortex, Thalamus);
        Stimulation.check_stim(t);
        if(t >= onset*res && t%red == 0){
            get_data(count, Cortex, Thalamus, dataPointer);
            ++count;
        }
    }

    /* Return the data containers */
    nlhs = dataArray.size()+1;
    for (auto &dataptr : dataArray) {
        plhs[std::distance(&dataptr, dataArray.data())] = dataptr;
    }
    plhs[dataArray.size()] = get_marker(Stimulation);

    return;
}

/******************************************************************************/
/*                          Create MATLAB data containers					  */
/******************************************************************************/
mxArray* GetMexArray(int N, int M) {
    mxArray* Array	= mxCreateDoubleMatrix(0, 0, mxREAL);
    mxSetM(Array, N);
    mxSetN(Array, M);
    {mxSetData(Array, mxMalloc(sizeof(double)*M*N));}
    return Array;
}

mxArray* get_marker(Stim &stim) {
    extern const int red;
    mxArray* marker	= mxCreateDoubleMatrix(0, 0, mxREAL);
    mxSetM(marker, 1);
    mxSetN(marker, stim.marker_stimulation.size());
    mxSetData(marker, mxMalloc(sizeof(double)*stim.marker_stimulation.size()));
    double* Pr_Marker = mxGetPr(marker);
    unsigned counter  = 0;
    /* Division by res transforms marker time from dt to sampling rate */
    for(auto & elem : stim.marker_stimulation) {
        Pr_Marker[counter++] = elem/red;
    }
    return marker;
}
