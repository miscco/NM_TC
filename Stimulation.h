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
	/* empty constructor for compiling */
	Stim(void);

	Stim(Cortical_Column& C, Thalamic_Column& T, double* var)
	{ Cortex   = &C;
	  Thalamus = &T;
	  setup(var);}

	/* setup with respect to stimulation mode */
	void setup		(double* var_stim) {
		extern const int onset;
		extern const int res;
		extern const int dt;

		/* mode of stimulation */
		mode		= (int) var_stim[0];

		/* scale the stimulation strength from s^-1 to ms^-1 */
		strength 	= 		var_stim[1] / 1000;

		/* scale duration from ms to dt */
		duration 	= (int) var_stim[2] * res / 1000;

		/* scale the ISI from s to ms^-1 */
		ISI 		= (int) var_stim[3] * res;

		/* scale time to stimulus from ms to dt */
		time_to_stim= (int) var_stim[4] * res / 1000;

		if(mode==1) {
			time_to_stim = (onset+1) * res;
		}

		correction = onset * res;

	}

	void check_stim	(int time) {

		/* check if stimulation should start */
		switch (mode) {

		/* no stimulation */
		default:
			break;

		/* periodic stimulation */
		case 1:
			/* check if time is reached */
			if(time == time_to_stim) {
				/* switch stimulation on */
				stimulation_started 	= true;
				Thalamus->set_input(strength);

				/* update the timer */
				time_to_stim += ISI;


				marker_threshold.push_back(0);
				marker_minimum.push_back(0);
				marker_stimulation.push_back(time - correction);
			}
			break;

		/* phase dependent up state stimulation */
		case 2:
			/* search for threshold */
			if(!stimulation_started && !minimum_found && !threshold_crossed && time>correction) {
				if(Cortex->Ve[0]<=threshold) {
					threshold_crossed 	= true;
					marker_threshold.push_back(time - correction);
				}
			}

			/* search for minimum */
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

			/* wait until the stimulation should start */
			if(minimum_found) {
				count_to_start++;


				/* start stimulation after time_to_stim has passed */
				if(count_to_start==time_to_stim) {
					minimum_found 			= false;
					stimulation_started 	= true;
					count_to_start 			= 0;
					marker_stimulation.push_back(time - correction);
					Thalamus->set_input(strength);
				}
			}
			break;

		/* phase dependent down state stimulation */
		case 3:
				/* search for threshold */
				if(!stimulation_started && !minimum_found && !threshold_crossed && time>correction) {
					if(Cortex->Ve[0]<=threshold) {
						threshold_crossed 		= true;
						marker_threshold.push_back(time - correction);
					}
				}

				/* search for minimum */
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

				/* start the stimulation */
				if(minimum_found) {
					minimum_found 			= false;
					stimulation_started 	= true;
					marker_stimulation.push_back(time - correction);
					Thalamus->set_input(strength);
				}
				break;
		}

		/* wait to switch the stimulation off */
		if(stimulation_started) {
			count_duration++;

			/* switch stimulation off */
			if(count_duration==duration) {
				stimulation_started 	= false;
				count_duration			= 0;
				Thalamus->set_input(0.0);
			}
		}
	}

	mxArray* get_marker(void) {
		mxArray* Marker	= mxCreateDoubleMatrix(0, 0, mxREAL);
	    mxSetM(Marker, 3);
	    mxSetN(Marker, marker_stimulation.size());
	    mxSetData(Marker, mxMalloc(sizeof(double)*3*marker_stimulation.size()));
		double* Pr_Marker	= mxGetPr(Marker);
		for(unsigned i=0; i<marker_stimulation.size(); ++i) {
			Pr_Marker[0+i*3] = marker_threshold[i];
			Pr_Marker[1+i*3] = marker_minimum[i];
			Pr_Marker[2+i*3] = marker_stimulation[i];
		}
		return Marker;
	}

private:

	/* Stimulation parameters */
	/* stimulation strength */
	double strength 	= 0.0;

	/* duration of the stimulation */
	int duration 		= 500;

	/* inter stimulus intervall in case of periodic stimulation */
	int ISI				= 5E4;

	/* threshold for phase dependent stimulation */
	double threshold	= -80;

	/* time until stimulus after minimum was found */
	int	time_to_stim	= 5500;

	/* mode of stimulation 				*/
	/* 0 == none 						*/
	/* 1 == periodic 					*/
	/* 2 == phase dependent up state 	*/
	/* 3 == phase dependent down state 	*/
	int mode			= 0;

	/* Internal variables */
	/* Simulation on for TRUE and off for FALSE */
	bool stimulation_started= false;

	/* threshold has been reached */
	bool threshold_crossed	= false;

	/* minimum found */
	bool minimum_found		= false;

	/* onset in timesteps to correct the given time of the markers */
	int correction			= 10000;

	/* counter for stimulation duration */
	int count_duration		= 0;

	/* counter after minimum */
	int count_to_start 		= 0;

	/* old voltage value for minimum detection */
	double Ve_old			= 0;

	/* pointer to cortical column */
	Cortical_Column* Cortex;

	/* pointer to thalamic column */
	Thalamic_Column* Thalamus;

	/* Data containers */
	std::vector<int>	marker_threshold;
	std::vector<int>	marker_minimum;
	std::vector<int>	marker_stimulation;
};
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/
