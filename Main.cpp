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
extern const int T		= 30;								/* Simulation length s					*/
extern const int res 	= 1E4;								/* Number of iteration steps per s		*/
extern const int red 	= 1E2;                              /* Fraction of iterations that is saved	*/
extern const double dt 	= 1E3/res;							/* Duration of a iteration step in ms	*/
extern const double h	= sqrt(dt);							/* Square root of dt for SRK iteration	*/
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										Main simulation routine										*/
/****************************************************************************************************/
int main(void) {
	/* Initializing the seeder */
	srand(time(0));

	/* Initialize the populations */
	Cortical_Column Cortex;
	Thalamic_Column Thalamus;

	/* Connect both modules */
	Cortex.get_Thalamus(Thalamus);
	Thalamus.get_Cortex(Cortex);

	/* Take the time of the simulation */
	time_t start,end;
	time (&start);
	/* Simulation */
	for (int t=0; t< T*res; ++t) {
		ODE (Cortex, Thalamus);
	}

	time (&end);
	/* Time consumed by the simulation */
	double dif = difftime(end,start);
	std::cout << "simulation done!\n";
	std::cout << "took " << dif 	<< " seconds" << "\n";
	std::cout << "end\n";
}
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/
