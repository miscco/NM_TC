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
 *	The aboVp copyright notice and this permission notice shall be included in
 *	all copies or substantial portions of the Software.
 *
 *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVpNT SHALL THE
 *	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *	THE SOFTWARE.
 *
 *	AUTHORS:	Michael Schellenberger Costa: mschellenbergercosta@gmail.com
 *
 *	Based on:	Characterization of K-Complexes and Slow Wave Activity in a Neural Mass Model
 *				A Weigenand, M Schellenberger Costa, H-VV Ngo, JC Claussen, T Martinetz
 *				PLoS Computational Biology. 2014;10:e1003923
 *
 *				A thalamocortical neural mass model of the EEG during NREM sleep and its response
 *				to auditory stimulation.
 *				M Schellenberger Costa, A Weigenand, H-VV Ngo, L Marshall, J Born, T Martinetz,
 *				JC Claussen.
 *				PLoS Computational Biology (in review).
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
		MTRands.push_back(random_stream_normal(0.0, dphi*dt));

		/* Add the RNG for I_{l,0} */
		MTRands.push_back(random_stream_normal(0.0, dt));

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
double Cortical_Column::get_Qp	(int N) const{
	double q = Qp_max / (1 + exp(-C1 * (Vp[N] - theta_p) / sigma_p));
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
double Cortical_Column::I_ep	(int N) const{
	double I = g_AMPA * s_ep[N] * (Vp[N] - E_AMPA);
	return I;
}

/* Inhibitory input to pyramidal population */
double Cortical_Column::I_gp	(int N) const{
	double I = g_GABA * s_gp[N] * (Vp[N] - E_GABA);
	return I;
}
/* Excitatory input to inhibitory population */
double Cortical_Column::I_ei	(int N) const{
	double I = g_AMPA * s_ei[N] * (Vi[N] - E_AMPA);
	return I;
}

/* Inhibitory input to inhibitory population */
double Cortical_Column::I_gi	(int N) const{
	double I = g_GABA * s_gi[N] * (Vi[N] - E_GABA);
	return I;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Current functions 											*/
/****************************************************************************************************/
/* Leak current of pyramidal population */
double Cortical_Column::I_L_p	(int N) const{
	double I = g_L * (Vp[N] - E_L_p);
	return I;
}

/* Leak current of inhibitory population */
double Cortical_Column::I_L_g	(int N) const{
	double I = g_L * (Vi[N] - E_L_g);
	return I;
}

/* Sodium dependent potassium current */
double Cortical_Column::I_KNa		(int N)  const{
	double w_KNa  = 0.37/(1+pow(38.7/Na[N], 3.5));
	double I_KNa  = g_KNa * w_KNa * (Vp[N] - E_K);
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
	Vp	[N+1] = Vp  [0] + A[N] * dt*(-(I_L_p(N) + I_ep(N) + I_gp(N))/tau_p - I_KNa(N));
	Vi	[N+1] = Vi  [0] + A[N] * dt*(-(I_L_g(N) + I_ei(N) + I_gi(N))/tau_i);
	Na	[N+1] = Na  [0] + A[N] * dt*(alpha_Na * get_Qp(N) - Na_pump(N))/tau_Na;
	s_ep[N+1] = s_ep[0] + A[N] * dt*(x_ep[N]);
	s_ei[N+1] = s_ei[0] + A[N] * dt*(x_ei[N]);
	s_gp[N+1] = s_gp[0] + A[N] * dt*(x_gp[N]);
	s_gi[N+1] = s_gi[0] + A[N] * dt*(x_gi[N]);
	y	[N+1] = y	[0] + A[N] * dt*(x	 [N]);
	x_ep[N+1] = x_ep[0] + A[N] * dt*(pow(gamma_e, 2) * (N_ep * get_Qp(N) + N_tp * Thalamus->y[N] - s_ep[N]) - 2 * gamma_e * x_ep[N]) + noise_xRK(N, 0);
	x_ei[N+1] = x_ei[0] + A[N] * dt*(pow(gamma_e, 2) * (N_ei * get_Qp(N) + N_ti * Thalamus->y[N] - s_ei[N]) - 2 * gamma_e * x_ei[N]) + noise_xRK(N, 1)	;
	x_gp[N+1] = x_gp[0] + A[N] * dt*(pow(gamma_g, 2) * (N_ip * get_Qi(N)						 - s_gp[N]) - 2 * gamma_g * x_gp[N]);
	x_gi[N+1] = x_gi[0] + A[N] * dt*(pow(gamma_g, 2) * (N_ii * get_Qi(N)						 - s_gi[N]) - 2 * gamma_g * x_gi[N]);
	x	[N+1] = x	[0] + A[N] * dt*(pow(nu, 	  2) * (	   get_Qp(N)						 - y   [N])	- 2 * nu	  * x   [N]);
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									Function that adds all SRK terms								*/
/****************************************************************************************************/
void Cortical_Column::add_RK(void) {
	Vp	[0] = (-3*Vp  [0] + 2*Vp  [1] + 4*Vp  [2] + 2*Vp  [3] + Vp	[4])/6;
	Vi	[0] = (-3*Vi  [0] + 2*Vi  [1] + 4*Vi  [2] + 2*Vi  [3] + Vi	[4])/6;
	Na	[0] = (-3*Na  [0] + 2*Na  [1] + 4*Na  [2] + 2*Na  [3] + Na	[4])/6;
	s_ep[0] = (-3*s_ep[0] + 2*s_ep[1] + 4*s_ep[2] + 2*s_ep[3] + s_ep[4])/6;
	s_ei[0] = (-3*s_ei[0] + 2*s_ei[1] + 4*s_ei[2] + 2*s_ei[3] + s_ei[4])/6;
	s_gp[0] = (-3*s_gp[0] + 2*s_gp[1] + 4*s_gp[2] + 2*s_gp[3] + s_gp[4])/6;
	s_gi[0] = (-3*s_gi[0] + 2*s_gi[1] + 4*s_gi[2] + 2*s_gi[3] + s_gi[4])/6;
	y	[0] = (-3*y	  [0] + 2*y	  [1] + 4*y   [2] + 2*y	  [3] + y	[4])/6;
	x_ep[0] = (-3*x_ep[0] + 2*x_ep[1] + 4*x_ep[2] + 2*x_ep[3] + x_ep[4])/6 + noise_aRK(0);
	x_ei[0] = (-3*x_ei[0] + 2*x_ei[1] + 4*x_ei[2] + 2*x_ei[3] + x_ei[4])/6 + noise_aRK(1);
	x_gp[0] = (-3*x_gp[0] + 2*x_gp[1] + 4*x_gp[2] + 2*x_gp[3] + x_gp[4])/6;
	x_gi[0] = (-3*x_gi[0] + 2*x_gi[1] + 4*x_gi[2] + 2*x_gi[3] + x_gi[4])/6;
	x	[0] = (-3*x	  [0] + 2*x	  [1] + 4*x   [2] + 2*x	  [3] + x	[4])/6;

	/* Generate noise for the next iteration */
	for (unsigned i=0; i<Rand_vars.size(); ++i) {
		Rand_vars[i] = MTRands[i]() + input;
	}
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
