/****************************************************************************************************/
/*								Implementation of the stimulation protocol							*/
/****************************************************************************************************/
#pragma once
#include "Thalamic_Column.h"

class Cortical_Column;

class Stim {
public:
	// empty constructor for compiling
	Stim(void);

	Stim(double* var)
	{setup(var);}

	// scaling from SI to simulation variables s -> ms
	void setup		(double* var_stim) {
		extern const int onset;
		extern const int res;
		extern const int dt;

		// scale the stimulation strength with respect to ms^-1
		strength = var_stim[0] / 1000;

		// scale the stimulation variables with respect to simulation resolution
		ISI 	 = (int) var_stim[1] * res;

		// stimulation starts after the onset
		start 	 = (int)(var_stim[2] + onset) *res;

		// rescale duration with respect to dt
		duration = (int) var_stim[3]/dt;
	}

	void check_stim	(int time) {
		// check whether a stimulation should start
		// the duration of the stimulation is ignored
		if(time==(start + count_stim*ISI)){
			// turn the stimulation on
			mode = 1;
			Thalamus->set_input(strength);

		}

		// check whether a stimulation should end
		if(mode ==1 && count_dur ==duration) {
			// turn off the stimulation
			mode = 0;
			Thalamus->set_input(0.0);

			// add counter for stimulation occurence
			count_stim++ ;

			// reset the stimulation counter
			count_dur = 0;
		}

		// if stimulation is on track its duration
		if(mode==1){
			count_dur++;
		}
	}

private:

	// stimulation strength
	double strength = 0.0;

	// inter stimulus intervall
	int ISI = 0;

	// onset until stimulation starts
	int start = 0;

	// duiration of the stimulation
	int duration = 0;

	// counter for stimulation events
	int count_stim = 0;

	// counter for stimulation length
	int count_dur  = 0;

	// Simulation on for TRUE and off for FALSE
	bool mode = 0;

	// Thalamic column
	Thalamic_Column* Thalamus;
};
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/
