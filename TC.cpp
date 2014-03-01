/****************************************************************************************************/
/*								main file for generation of mex function							*/
/****************************************************************************************************/
#include <vector>
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
extern const int 	onset 	= 5;
extern const int 	res 	= 1E4;
extern const double dt 		= 1E3/res;
extern const double h		= sqrt(dt);
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
	double* var_stim	 	= mxGetPr (prhs[1]);
	const int Time 			= (T+onset)*res;

	// Initializing the populations;
	Cortical_Column Cortex;
	Thalamic_Column Thalamus;

	// Initialize the stimulation protocoll
	Stim	Stimulation(Var_Stim, Thalamus);

	// setting up the data containers
	mxArray* Ve		= SetMexArray(1, T*red);

	int count = 0;
	for (int t=0; t<Time; ++t) {
		ODE (Cortex, Thalamus);

		if(t>=onset*res){
			get_data(count, Cortex, Ve);
			++count;
		}
	}

	// output of the simulation
	plhs[0] = getMexArray(Ve);
return;
}
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/
