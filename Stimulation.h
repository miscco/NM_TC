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
/*								Implementation of the stimulation protocol							*/
/****************************************************************************************************/
#pragma once
#include "Cortical_Column.h"
#include "Thalamic_Column.h"

/****************************************************************************************************/
/*											Stimulation object										*/
/****************************************************************************************************/
class Stim {
public:
	/* Empty constructor for compiling */
	Stim(void);

	Stim(Cortical_Column& C, Thalamic_Column& T, double* var)
	{ Cortex   = &C;
	  Thalamus = &T;
	  setup(var);}

	/* Setup with respect to stimulation mode */
	void setup		(double* var_stim) {
		extern const int onset;
		extern const int res;
		extern const int dt;

		/* Mode of stimulation */
		mode		= (int) var_stim[0];

		/* Scale the stimulation strength from s^-1 (Hz) to ms^-1 */
		strength 	= 		var_stim[1] / 1000;

		/* Scale duration from ms to dt */
		duration 	= (int) var_stim[2] * res / 1000;

		/* Scale the ISI from s to ms */
		ISI 		= (int) var_stim[3] * res;

		/* Scale time to stimulus from ms to dt */
		time_to_stim= (int) var_stim[4] * res / 1000;

		/* Set the onset correction for the marker */
		correction = onset * res;

		if(mode==1) {
			time_to_stim = (int) (onset+1) * res;
		}

		correction = onset * res;

	}

	void check_stim	(int time) {

		/* Check if stimulation should start */
		switch (mode) {

		/* No stimulation */
		default:
			break;

		/* Periodic stimulation */
		case 1:
			/* Check if stimulation time is reached */
			if(time == time_to_stim) {
				/* Switch stimulation on */
				stimulation_started 	= true;
				Thalamus->set_input(strength);

				/* Update the timer */
				time_to_stim += ISI;

				/* Add marker */
				marker_threshold.push_back(0);
				marker_minimum.push_back(0);
				marker_stimulation.push_back(time - correction);
			}
			break;

		/* Phase dependent up state stimulation */
		case 2:
			/* Search for threshold */
			if(!stimulation_started && !minimum_found && !threshold_crossed && time>correction) {
				if(Cortex->Ve[0]<=threshold) {
					threshold_crossed 	= true;
					marker_threshold.push_back(time - correction);
				}
			}

			/* Search for minimum */
			if(threshold_crossed) {
				if(Cortex->Ve[0]>Ve_old) {
					threshold_crossed 	= false;
					minimum_found 		= true;
					marker_minimum.push_back(time - correction);
					Ve_old = 0;
				} else {
					Ve_old = Cortex->Ve[0];
				}
			}

			/* Wait until the stimulation should start */
			if(minimum_found) {
				count_to_start++;


				/* Start stimulation after time_to_stim has passed */
				if(count_to_start==time_to_stim) {
					minimum_found 			= false;
					stimulation_started 	= true;
					count_to_start 			= 0;
					marker_stimulation.push_back(time - correction);
					Thalamus->set_input(strength);
				}
			}
			break;

		/* Phase dependent down state stimulation */
		case 3:
				/* Search for threshold */
				if(!stimulation_started && !minimum_found && !threshold_crossed && time>correction) {
					if(Cortex->Ve[0]<=threshold) {
						threshold_crossed 		= true;
						marker_threshold.push_back(time - correction);
					}
				}

				/* Search for minimum */
				if(threshold_crossed) {
					if(Cortex->Ve[0]>Ve_old) {
						threshold_crossed 		= false;
						minimum_found 			= true;
						marker_minimum.push_back(time - correction);
						Ve_old = 0;
					} else {
						Ve_old = Cortex->Ve[0];
					}
				}

				/* Start the stimulation */
				if(minimum_found) {
					minimum_found 			= false;
					stimulation_started 	= true;
					marker_stimulation.push_back(time - correction);
					Thalamus->set_input(strength);
				}
				break;
		}

		/* Wait to switch the stimulation off */
		if(stimulation_started) {
			count_duration++;

			/* Switch stimulation off */
			if(count_duration==duration) {
				stimulation_started 	= false;
				count_duration			= 0;
				Thalamus->set_input(0.0);
			}
		}
	}

	/* Create MATLAB container for marker storage */
	mxArray* get_marker(void) {
		extern const int res;
		mxArray* Marker	= mxCreateDoubleMatrix(0, 0, mxREAL);
	    mxSetM(Marker, 3);
	    mxSetN(Marker, marker_stimulation.size());
	    mxSetData(Marker, mxMalloc(sizeof(double)*3*marker_stimulation.size()));
		double* Pr_Marker	= mxGetPr(Marker);
		for(unsigned i=0; i<marker_stimulation.size(); ++i) {
			/* Division by res transforms marker time from sample rate to sec */
			Pr_Marker[0+i*3] = marker_threshold[i]/res;
			Pr_Marker[1+i*3] = marker_minimum[i]/res;
			Pr_Marker[2+i*3] = marker_stimulation[i]/res;
		}
		return Marker;
	}

private:
	/* Stimulation strength */
	double strength 	= 0.0;

	/* Duration of the stimulation */
	int duration 		= 500;

	/* Inter-stimulus-interval in case of periodic stimulation */
	int ISI				= 5E4;

	/* Threshold for phase dependent stimulation */
	double threshold	= -70;

	/* Time until stimulus after minimum was found */
	int	time_to_stim	= 900;

	/* Mode of stimulation 				*/
	/* 0 == none 						*/
	/* 1 == periodic 					*/
	/* 2 == phase dependent up state 	*/
	/* 3 == phase dependent down state 	*/
	int mode			= 0;

	/* Internal variables */
	/* Simulation on for TRUE and off for FALSE */
	bool stimulation_started= false;

	/* Threshold has been reached */
	bool threshold_crossed	= false;

	/* Minimum found */
	bool minimum_found		= false;

	/* Onset in time steps to correct the given time of the markers */
	int correction			= 10000;

	/* Counter for stimulation duration */
	int count_duration		= 0;

	/* Counter after minimum */
	int count_to_start 		= 0;

	/* Old voltage value for minimum detection */
	double Ve_old			= 0;

	/* Pointer to cortical column */
	Cortical_Column* Cortex;

	/* Pointer to thalamic column */
	Thalamic_Column* Thalamus;

	/* Data containers */
	std::vector<double>	marker_threshold;
	std::vector<double>	marker_minimum;
	std::vector<double>	marker_stimulation;
};
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/
