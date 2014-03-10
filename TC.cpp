/********************************************************************************************************/
/* 	Implementation of the simulation as MATLAB routine (mex compiler)									*/
/* 	mex command is given by:																			*/
/* 	mex CXXFLAGS="\$CXXFLAGS -std=gnu++0x -fpermissive" TC.cpp Cortical_Column.cpp Thalamic_Column.cpp	*/
/*	The Simulation requires the following boost libraries:	Preprocessor								*/
/*															Random										*/
/********************************************************************************************************/
#include <ctime>
#include "mex.h"
#include "matrix.h"
#include "Cortical_Column.h"
#include "Thalamic_Column.h"
#include "Stimulation.h"
#include "saves.h"
#include "ODE.h"
using std::vector;

/****************************************************************************************************/
/*									fixed simulation parameters										*/
/****************************************************************************************************/
extern const int onset	= 5;
extern const int res 	= 1E4;
extern const int red 	= res/100;
extern const double dt 	= 1E3/res;
extern const double h	= sqrt(dt);
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										Simulation routine											*/
/****************************************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	// Initializing the seeder.
	srand(time(0));

	// inputs
	const int T				= (int)(mxGetScalar(prhs[0]));
	double* Param_Cortex	= mxGetPr (prhs[1]);
	double* Param_Thalamus	= mxGetPr (prhs[2]);
	double* var_stim	 	= mxGetPr (prhs[3]);
	const int Time 			= (T+onset)*res;

	// Initializing the populations;
	Cortical_Column Cortex(Param_Cortex);
	Thalamic_Column Thalamus(Param_Thalamus);

	// Linking both modules
	Cortex.get_Thalamus(Thalamus);
	Thalamus.get_Cortex(Cortex);

	// Initialize the stimulation protocol
	Stim	Stimulation(Thalamus, var_stim);

	// setting up the data containers
	mxArray* Ve		= SetMexArray(1, T*res/red);
	mxArray* Vt		= SetMexArray(1, T*res/red);

	// Pointer to the actual data block
	double* Pr_Ve	= mxGetPr(Ve);
	double* Pr_Vt	= mxGetPr(Vt);

	int count = 0;
	for (int t=0; t<Time; ++t) {
		ODE (Cortex, Thalamus);
		Stimulation.check_stim(t);
		if(t>=onset*res && t%red==0){
			get_data(count, Cortex, Thalamus, Pr_Ve, Pr_Vt);
			++count;
		}
	}

	// output of the simulation
	plhs[0] = Ve;
	plhs[1] = Vt;
return;
}
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/
