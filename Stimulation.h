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
/*					Implementation of the stimulation protocol				  */
/******************************************************************************/
#pragma once
#include <vector>

#include "Cortical_Column.h"
#include "Random_Stream.h"
#include "Thalamic_Column.h"

/******************************************************************************/
/*								Stimulation object							  */
/******************************************************************************/
class Stim {
public:
    /* Constructor with references and stimulation variables */
    Stim(Cortical_Column& C, Thalamic_Column& T, double* var)
    { Cortex   = &C;
      Thalamus = &T;
      setup(var);}

    /* Initialize stimulation class with respect to stimulation mode */
    void setup		(double* var_stim);

    /* Check whether stimulation should be started/stopped */
    void check_stim	(int time);
private:
    /* Mode of stimulation 	*/
    /* 0 == none 			*/
    /* 1 == semi-periodic	*/
    /* 2 == phase dependent */
    int     mode                    = 0;

    /* Default values already in dt: E1==ms,  E4==s	*/
    /* Stimulation strength 						*/
    double 	strength 				= 0.0;

    /* Duration of the stimulation */
    int 	duration 				= 120E1;

    /* Interval between different Stimulus events */
    int 	ISI						= 5E4;

    /* Range of Inter Stimulus Interval */
    int 	ISI_range 				= 1E4;

    /* Number of stimuli in case of multiple stimuli per event */
    int		number_of_stimuli 		= 1;

    /* Time until next stimulus */
    /* Function varies between different stimulation modes */
    int		time_to_stimuli			= 350E1;

    /* Time between stimuli in case of multiple stimuli per event */
    int 	time_between_stimuli 	= 1050E1;

    /* Threshold for phase dependent stimulation */
    double 	threshold				= -68;

    /* Internal variables */
    /* Simulation on for TRUE and off for FALSE */
    bool 	stimulation_started 	= false;

    /* If threshold has been reached */
    bool 	threshold_crossed		= false;

    /* If minimum was found */
    bool 	minimum_found			= false;

    /* If a stimulation event has occurred and there is a minimal time (pause) until the next one */
    bool 	stimulation_paused 		= false;

    /* If burst mode is enabled */
    bool    burst_enabled            = false;

    /* In case of bursted stimulation start the bursts */
    bool 	burst_started 			= true;

    /* Length of a burst stimulus */
    int 	burst_length 			= 10;

    /* Length of silence between burst stimuli */
    int		burst_ISI 				= 10;

    /* Onset in time steps where no data is recorded, so stimulation has to be delayed */
    int 	onset_correction		= 10E4;

    /* Counter for number of stimuli that occurred within a stimulation event */
    int 	count_stimuli 			= 1;

    /* Counter for number of bursts per stimuli */
    int 	count_bursts 			= 0;

    /* Counter for stimulation duration since started*/
    int 	count_duration			= 0;

    /* Counter after minimum was found */
    int 	count_to_start 			= 0;

    /* Counter for time between two stimulation events (with multiple tones) */
    int 	count_pause 			= 0;

    /* Old voltage value for minimum detection */
    double 	Vp_old					= 0.0;

    /* Pointer to columns */
    Cortical_Column* Cortex;
    Thalamic_Column* Thalamus;

    /* Data containers */
    std::vector<int>		marker_stimulation;

    /* Random number generator in case of semi-periodic stimulation */
    randomStreamUniformInt Uniform_Distribution = randomStreamUniformInt(0, 0);

    /* Create MATLAB container for marker storage */
    friend mxArray* get_marker(Stim &stim);
};

/******************************************************************************/
/*							Function definitions							  */
/******************************************************************************/
void Stim::setup (double* var_stim) {
    extern const int onset;
    extern const int res;

    /* Set the onset onset_correction for the marker */
    onset_correction 		= onset * res;

    /* Mode of stimulation */
    mode					= (int) var_stim[0];

    /* Scale the stimulation strength from s^-1 (Hz) to ms^-1 */
    strength 				= 		var_stim[1] / 1000;

    /* Scale duration from ms to dt */
    duration 				= (int) var_stim[2] * res / 1000;

    /* Scale the inter stimulus event interval from s to dt */
    ISI 					= (int) var_stim[3] * res;

    /* Scale inter stimulus event interval range from s to dt */
    ISI_range 				= (int) var_stim[4] * res;

    /* Number of stimuli per Stimulus event */
    number_of_stimuli		= (int) var_stim[5];

    /* Scale time_between_stimuli from ms to dt */
    time_between_stimuli 	= (int) var_stim[6] * res / 1000;

    /* Scale the length of burst_length and burst_ISI from ms to dt*/
    burst_length 			= (int) 2  * res / 1000;
    burst_ISI 				= (int) 28 * res / 1000;

    /* If ISI is fixed do not create RNG */
    if (mode == 1) {
        /* Set first time_to_stimuli to 1 sec after onset */
        time_to_stimuli = (int) (onset+1) * res;

        /* If ISI is random create RNG */
        if (ISI_range != 0){
            /* Generate uniform distribution */
            Uniform_Distribution = randomStreamUniformInt(ISI-ISI_range, ISI+ISI_range);
        }
    } else {
        /* In case of phase dependent stimulation, time_to_stim is the time from minimum detection to start of stimulation */
        /* Scale time_to_stimuli from ms to dt */
        time_to_stimuli = (int) var_stim[7] * res / 1000;
    }
}

void Stim::check_stim	(int time) {
    /* Check if stimulation should start */
    switch (mode) {

    /* No stimulation */
    default:
        break;

    /* Semi-periodic stimulation */
    case 1:
        /* Check if stimulation time is reached */
        if(time == time_to_stimuli) {
            /* Switch stimulation on */
            stimulation_started 	= true;
            Thalamus->set_input(strength);

            /* Add marker for the first stimuli in the event */
            if(count_stimuli == 1) {
                marker_stimulation.push_back(time - onset_correction);
            }

            /* Check if multiple stimuli should be applied */
            if (count_stimuli < number_of_stimuli) {
                /* Update the timer with respect to time between stimuli */
                time_to_stimuli += time_between_stimuli;
                count_stimuli++;
            }
            /* After last stimulus in event update the timer with respect to (random) ISI*/
            else {
                time_to_stimuli += (ISI_range==0)? ISI : Uniform_Distribution();

                /* Reset the stimulus counter for next stimulation event */
                count_stimuli = 1;
            }
        }
        break;

    /* Phase dependent stimulation */
    case 2:
        /* Search for threshold */
        if(!stimulation_started && !minimum_found && !threshold_crossed && time>onset_correction && !stimulation_paused) {
            if(Cortex->Vp[0]<=threshold) {
                threshold_crossed 	= true;
            }
        }

        /* Search for minimum */
        if(threshold_crossed) {
            if(Cortex->Vp[0]>Vp_old) {
                threshold_crossed 	= false;
                minimum_found 		= true;
                Vp_old = 0;
            } else {
                Vp_old = Cortex->Vp[0];
            }
        }

        /* Wait until the stimulation should start */
        if(minimum_found) {
            /* Start stimulation after time_to_stimuli has passed */
            if(count_to_start==time_to_stimuli + (count_stimuli-1) * time_between_stimuli) {
                stimulation_started 	= true;
                Thalamus->set_input(strength);

                /* Add marker for the first stimuli in the event */
                if(count_stimuli == 1) {
                    marker_stimulation.push_back(time - onset_correction);
                }

                /* Check if multiple stimuli should be applied */
                if (count_stimuli < number_of_stimuli) {
                    /* Update the number of stimuli */
                    count_stimuli++;
                } else {
                    /* After last stimulus in event pause the stimulation */
                    minimum_found 			= false;
                    stimulation_paused 		= true;
                    count_to_start 			= 0;

                    /* Reset the stimulus counter for next stimulation event */
                    count_stimuli = 1;
                }
            }
            /* Update counter */
            count_to_start++;
        }
        break;
    }

    /* Actual stimulation protocols */
    if(stimulation_started) {
        /* Wait to switch the stimulation off */
        if(count_duration==duration) {
            stimulation_started 	= false;
            burst_started 			= true;
            count_duration			= 0;
            count_bursts			= 0;
            Thalamus->set_input(0.0);
        }

        count_duration++;
        count_bursts++;

        /* Switch stimulation on and off wrt burst parameters */
        if(burst_enabled) {
            if(burst_started) {
                if(count_bursts%burst_length==0) {
                    count_bursts 	= 0;
                    burst_started 	= false;
                    Thalamus->set_input(0.0);
                }
            } else {
                if(count_bursts%burst_ISI==0) {
                    count_bursts 	= 0;
                    burst_started	= true;
                    Thalamus->set_input(strength);
                }
            }
        }
    }

    /* Wait if there is a pause between stimulation events */
    if(stimulation_paused) {
        if(count_pause == ISI) {
            stimulation_paused	= false;
            count_pause 		= 0;
        }
        count_pause++;
    }
}
