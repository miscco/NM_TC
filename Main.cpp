/****************************************************************************************************/
/*		Main file for compilation tests																*/
/*		The Simulation requires the following boost libraries:	Preprocessor						*/
/*																Random								*/
/****************************************************************************************************/
#include <iostream>
#include <ctime>
#include "Cortical_Column.h"
#include "Thalamic_Column.h"
#include "ODE.h"

/****************************************************************************************************/
/*										Fixed simulation settings									*/
/****************************************************************************************************/
extern const int T 		= 60;
extern const int res 	= 1E4;
extern const double dt 	= 1E3/res;
extern const double h	= sqrt(dt);
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										Main simulation routine										*/
/****************************************************************************************************/
int main(void) {
	// Initializing the seeder.
	srand(time(0));

	// Initializing the populations;
	Cortical_Column Cortex;
	Thalamic_Column Thalamus;

	// Connect both modules
	Cortex.get_Thalamus(Thalamus);
	Thalamus.get_Cortex(Cortex);

	// takes the time of the simulation
	time_t start,end;
	time (&start);
	// simulation
	for (int t=0; t< T*res; ++t) {
		ODE (Cortex, Thalamus);
	}

	time (&end);
	// time consumed by the simulation
	double dif = difftime(end,start);
	std::cout << "simulation done!\n";
	std::cout << "took " << dif 	<< " seconds" << "\n";
	std::cout << "end\n";
}
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/
