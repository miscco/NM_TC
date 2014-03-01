/****************************************************************************************************/
/*									C++ main file for code verification								*/
/****************************************************************************************************/
#include <iostream>
#include <ctime>
#include <vector>

#include "Cortical_Column.h"
#include "Thalamic_Column.h"
#include "Stimulation.h"
#include "ODE.h"

using std::vector;

extern const int 	T 		= 100;
extern const int 	onset 	= 5;
extern const int 	res 	= 1E4;
extern const double dt 		= 1E3/res;
extern const double h		= sqrt(dt);

// simulation of the thalamo-cortical model
int main(void) {
	// Initializing the seeder.
	srand(time(0));

	// Initializing the populations;
	Cortical_Column Cortex;
	Thalamic_Column Thalamus;

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
