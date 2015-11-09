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
 *				PLoS Computational Biology In Review (in review).
 */

/****************************************************************************************************/
/*									Functions of the cortical module								*/
/****************************************************************************************************/
#include "Cortical_Column.h"

/****************************************************************************************************/
/*										 Initialization of RNG 										*/
/****************************************************************************************************/
void Cortical_Column::set_RNG(void) {
	extern const double dt;
	/* Number of independent random variables */
	int N = 2;

	/* Create RNG for each stream */
	for (int i=0; i<N; ++i){
		/* Add the RNG for I_{l}*/
		MTRands.push_back({ENG(rand()), DIST (0.0, dphi*dt)});

		/* Add the RNG for I_{l,0} */
		MTRands.push_back({ENG(rand()), DIST (0.0, dt)});

		/* Get the random number for the first iteration */
		Rand_vars.push_back(MTRands[2*i]());
		Rand_vars.push_back(MTRands[2*i+1]());
	}
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Firing Rate functions 										*/
/****************************************************************************************************/
/* Pyramidal firing rate */
double Cortical_Column::get_Qe	(int N) const{
	double q = Qe_max / (1 + exp(-C1 * (Ve[N] - theta_e) / sigma_e));
	return q;
}

/* Inhibitory firing rate */
double Cortical_Column::get_Qi	(int N) const{
	double q = Qi_max / (1 + exp(-C1 * (Vi[N] - theta_i) / sigma_i));
	return q;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										Synaptic currents											*/
/****************************************************************************************************/
/* Excitatory input to pyramidal population */
double Cortical_Column::I_ee	(int N) const{
	double I = y_ee[N] * (Ve[N] - E_AMPA);
	return I;
}

/* Inhibitory input to pyramidal population */
double Cortical_Column::I_ie	(int N) const{
	double I = y_ie[N] * (Ve[N] - E_GABA);
	return I;
}
/* Excitatory input to inhibitory population */
double Cortical_Column::I_ei	(int N) const{
	double I = y_ei[N] * (Vi[N] - E_AMPA);
	return I;
}

/* Inhibitory input to inhibitory population */
double Cortical_Column::I_ii	(int N) const{
	double I = y_ii[N] * (Vi[N] - E_GABA);
	return I;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Current functions 											*/
/****************************************************************************************************/
/* Leak current of pyramidal population */
double Cortical_Column::I_L_e	(int N) const{
	double I = g_L * (Ve[N] - E_L_e);
	return I;
}

/* Leak current of inhibitory population */
double Cortical_Column::I_L_i	(int N) const{
	double I = g_L * (Vi[N] - E_L_i);
	return I;
}

/* Sodium dependent potassium current */
double Cortical_Column::I_KNa		(int N)  const{
	double w_KNa  = 0.37/(1+pow(38.7/Na[N], 3.5));
	double I_KNa  = g_KNa * w_KNa * (Ve[N] - E_K);
	return I_KNa;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									 		Potassium pump	 										*/
/****************************************************************************************************/
double Cortical_Column::Na_pump		(int N) const{
	double Na_pump = R_pump*(Na[N]*Na[N]*Na[N]/(Na[N]*Na[N]*Na[N]+3375) - Na_eq*Na_eq*Na_eq/(Na_eq*Na_eq*Na_eq+3375));
	return Na_pump;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 RK noise scaling 											*/
/****************************************************************************************************/
double Cortical_Column::noise_xRK(int N, int M) const{
	return gamma_e * gamma_e * (Rand_vars[2*M] + Rand_vars[2*M+1]/std::sqrt(3))*B[N];
}

double Cortical_Column::noise_aRK(int M) const{
	return gamma_e * gamma_e * (Rand_vars[2*M] - Rand_vars[2*M+1]*std::sqrt(3))/4;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										Calculate the Nth SRK term									*/
/****************************************************************************************************/
void Cortical_Column::set_RK (int N) {
	extern const double dt;
	Ve	[N+1] = Ve  [0] + A[N] * dt*(-(I_L_e(N) + I_ee(N) + I_ie(N))/tau_e - I_KNa(N));
	Vi	[N+1] = Vi  [0] + A[N] * dt*(-(I_L_i(N) + I_ei(N) + I_ii(N))/tau_i);
	Na	[N+1] = Na  [0] + A[N] * dt*(alpha_Na * get_Qe(N) - Na_pump(N))/tau_Na;
	y_ee[N+1] = y_ee[0] + A[N] * dt*(x_ee[N]);
	y_ei[N+1] = y_ei[0] + A[N] * dt*(x_ei[N]);
	y_ie[N+1] = y_ie[0] + A[N] * dt*(x_ie[N]);
	y_ii[N+1] = y_ii[0] + A[N] * dt*(x_ii[N]);
	y	[N+1] = y	[0] + A[N] * dt*(x	 [N]);
	x_ee[N+1] = x_ee[0] + A[N] * dt*(pow(gamma_e, 2) * (N_ee * get_Qe(N) + N_te * Thalamus->y[N] - y_ee[N]) - 2 * gamma_e * x_ee[N]) + noise_xRK(N, 0);
	x_ei[N+1] = x_ei[0] + A[N] * dt*(pow(gamma_e, 2) * (N_ei * get_Qe(N) + N_ti * Thalamus->y[N] - y_ei[N]) - 2 * gamma_e * x_ei[N]) + noise_xRK(N, 1)	;
	x_ie[N+1] = x_ie[0] + A[N] * dt*(pow(gamma_i, 2) * (N_ie * get_Qi(N)						 - y_ie[N]) - 2 * gamma_i * x_ie[N]);
	x_ii[N+1] = x_ii[0] + A[N] * dt*(pow(gamma_i, 2) * (N_ii * get_Qi(N)						 - y_ii[N]) - 2 * gamma_i * x_ii[N]);
	x	[N+1] = x	[0] + A[N] * dt*(pow(nu, 	  2) * (	   get_Qe(N)						 - y   [N])	- 2 * nu	  * x   [N]);
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									Function that adds all SRK terms								*/
/****************************************************************************************************/
void Cortical_Column::add_RK(void) {
	Ve	[0] = (-3*Ve  [0] + 2*Ve  [1] + 4*Ve  [2] + 2*Ve  [3] + Ve	[4])/6;
	Vi	[0] = (-3*Vi  [0] + 2*Vi  [1] + 4*Vi  [2] + 2*Vi  [3] + Vi	[4])/6;
	Na	[0] = (-3*Na  [0] + 2*Na  [1] + 4*Na  [2] + 2*Na  [3] + Na	[4])/6;
	y_ee[0] = (-3*y_ee[0] + 2*y_ee[1] + 4*y_ee[2] + 2*y_ee[3] + y_ee[4])/6;
	y_ei[0] = (-3*y_ei[0] + 2*y_ei[1] + 4*y_ei[2] + 2*y_ei[3] + y_ei[4])/6;
	y_ie[0] = (-3*y_ie[0] + 2*y_ie[1] + 4*y_ie[2] + 2*y_ie[3] + y_ie[4])/6;
	y_ii[0] = (-3*y_ii[0] + 2*y_ii[1] + 4*y_ii[2] + 2*y_ii[3] + y_ii[4])/6;
	y	[0] = (-3*y	  [0] + 2*y	  [1] + 4*y   [2] + 2*y	  [3] + y	[4])/6;
	x_ee[0] = (-3*x_ee[0] + 2*x_ee[1] + 4*x_ee[2] + 2*x_ee[3] + x_ee[4])/6 + noise_aRK(0);
	x_ei[0] = (-3*x_ei[0] + 2*x_ei[1] + 4*x_ei[2] + 2*x_ei[3] + x_ei[4])/6 + noise_aRK(1);
	x_ie[0] = (-3*x_ie[0] + 2*x_ie[1] + 4*x_ie[2] + 2*x_ie[3] + x_ie[4])/6;
	x_ii[0] = (-3*x_ii[0] + 2*x_ii[1] + 4*x_ii[2] + 2*x_ii[3] + x_ii[4])/6;
	x	[0] = (-3*x	  [0] + 2*x	  [1] + 4*x   [2] + 2*x	  [3] + x	[4])/6;

	/* Generate noise for the next iteration */
	for (unsigned i=0; i<Rand_vars.size(); ++i) {
		Rand_vars[i] = MTRands[i]() + input;
	}
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
