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
#include "Thalamic_Column.h"

/****************************************************************************************************/
/*										 Initialization of RNG 										*/
/****************************************************************************************************/
void Thalamic_Column::set_RNG(void) {
	extern const double dt;
	/* Number of independent random variables */
	int N = 1;

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
/* Thalamic relay (TC) firing rate */
double Thalamic_Column::get_Qt	(int N) const{
	double Q = Qt_max/ (1 + exp(-C1 * (Vt[N] - theta_t) / sigma_t));
	return Q;
}

/* Thalamic reticular (RE) firing rate */
double Thalamic_Column::get_Qr	(int N) const{
	double Q = Qr_max / (1 + exp(-C1 * (Vr[N]- theta_r) / sigma_r));
	return Q;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										Synaptic currents											*/
/****************************************************************************************************/
/* Excitatory input to TC population */
double Thalamic_Column::I_et	(int N) const{
	double I = g_AMPA * y_et[N]* (Vt[N]- E_AMPA);
	return I;
}

/* Inhibitory input to TC population */
double Thalamic_Column::I_gt	(int N) const{
	double I = g_AMPA * y_rt[N]* (Vt[N]- E_GABA);
	return I;
}
/* Excitatory input to RE population */
double Thalamic_Column::I_er	(int N) const{
	double I = g_GABA * y_er[N]* (Vr[N]- E_AMPA);
	return I;
}

/* Inhibitory input to RE population */
double Thalamic_Column::I_gr	(int N) const{
	double I = g_GABA * y_rr[N]* (Vr[N]- E_GABA);
	return I;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 	I_T gating	 											*/
/****************************************************************************************************/
/* Activation in TC population after Destexhe 1996 */
double Thalamic_Column::m_inf_T_t	(int N) const{
	double m = 1/(1+exp(-(Vt[N]+59)/6.2));
	return m;
}

/* Activation in RE population after Destexhe 1996 */
double Thalamic_Column::m_inf_T_r	(int N) const{
	double m = 1/(1+exp(-(Vr[N]+52)/7.4));
	return m;
}

/* Deactivation in TC population after Destexhe 1996 */
double Thalamic_Column::h_inf_T_t	(int N) const{
	double h = 1/(1+exp( (Vt[N]+81)/4));
	return h;
}

/* Deactivation in RE population after Destexhe 1996 */
double Thalamic_Column::h_inf_T_r	(int N) const{
	double h = 1/(1+exp( (Vr[N]+80)/5));
	return h;
}

/* deactivation time in RE population after Destexhe 1996 */
double Thalamic_Column::tau_h_T_t	(int N) const{
	double tau =  (30.8 + (211.4 + exp((Vt[N]+115.2)/5))/(1 + exp((Vt[N]+86)/3.2)))/pow(3, 1.2);
	return tau;
}

/* Deactivation time in RE population after Destexhe 1996 */
double Thalamic_Column::tau_h_T_r	(int N) const{
	double tau =  (85 + 1/(exp((Vr[N]+48)/4) + exp(-(Vr[N]+407)/50)))/pow(3, 1.2);
	return tau;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 	I_h gating 												*/
/****************************************************************************************************/
/* Activation in TC population after Destexhe 1993 */
double Thalamic_Column::m_inf_h	(int N) const{
	double h = 1/(1+exp( (Vt[N]+75)/5.5));
	return h;
}

/* Activation time for slow components in TC population after Chen2012 */
double Thalamic_Column::tau_m_h	(int N) const{
	double tau = (20 + 1000/(exp((Vt[N]+ 71.5)/14.2) + exp(-(Vt[N]+ 89)/11.6)));
	return tau;
}

/* Instantaneous calcium binding onto messenger protein after Chen2012 */
double Thalamic_Column::P_h	(int N) const{
	//double P_h = k1 * pow(Ca[N], n_P)/(k1*pow(Ca[N], n_P)+k2);
	double P_h = k1 * Ca[N]* Ca[N]* Ca[N]/(k1* Ca[N]* Ca[N]* Ca[N]+k2);
	return P_h;
}

/* Return I_h activation */
double Thalamic_Column::act_h	(void) const{
	double activation = m_h[0] + g_inc * m_h2[0];
	return activation;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Current functions 											*/
/****************************************************************************************************/
/* Leak current of TC population */
double Thalamic_Column::I_L_t	(int N) const{
	double I = g_L * (Vt[N]- E_L_t);
	return I;
}

/* Potassium leak current of TC population */
double Thalamic_Column::I_LK_t	(int N) const{
	double I = g_LK * (Vt[N]- E_K);
	return I;
}

/* Leak current of RE population */
double Thalamic_Column::I_L_r	(int N) const{
	double I = g_L * (Vr[N]- E_L_r);
	return I;
}

/* Potassium leak current of RE population */
double Thalamic_Column::I_LK_r	(int N) const{
	double I = g_LK	* (Vr[N]- E_K);
	return I;
}

/* T-type current of TC population */
double Thalamic_Column::I_T_t	(int N) const{
	double I = g_T_t * m_inf_T_t(N) * m_inf_T_t(N) * h_T_t[N] * (Vt[N]- E_Ca);
	return I;
}

/* T-type current of RE population */
double Thalamic_Column::I_T_r	(int N) const{
	double I = g_T_r * m_inf_T_r(N) * m_inf_T_r(N) * h_T_r[N] * (Vr[N]- E_Ca);
	return I;
}

/* h-type current of TC population */
double Thalamic_Column::I_h		(int N) const{
	double I = g_h * (m_h[N] + g_inc * m_h2[N]) * (Vt[N]- E_h);
	return I;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 RK noise scaling 											*/
/****************************************************************************************************/
double Thalamic_Column::noise_xRK(int N, int M) const{
	return gamma_e * gamma_e * (Rand_vars[2*M] + Rand_vars[2*M+1]/std::sqrt(3))*B[N];
}

double Thalamic_Column::noise_aRK(int M) const{
	return gamma_e * gamma_e * (Rand_vars[2*M] - Rand_vars[2*M+1]*std::sqrt(3))/4;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										Calculate the Nth SRK term									*/
/****************************************************************************************************/
void Thalamic_Column::set_RK (int N) {
	extern const double dt;
	Vt	  	[N+1] = Vt   [0] + A[N]*dt*(-(I_L_t(N) + I_et(N) + I_gt(N))/tau_t - (I_LK_t(N) + I_T_t(N) + I_h(N)));
	Vr	  	[N+1] = Vr   [0] + A[N]*dt*(-(I_L_r(N) + I_er(N) + I_gr(N))/tau_r - (I_LK_r(N) + I_T_r(N)));
	Ca      [N+1] = Ca   [0] + A[N]*dt*(alpha_Ca * I_T_t(N) - (Ca[N] - Ca_0)/tau_Ca);
	y_et	[N+1] = y_et [0] + A[N]*dt*(x_et[N]);
	y_er	[N+1] = y_er [0] + A[N]*dt*(x_er[N]);
	y_rt	[N+1] = y_rt [0] + A[N]*dt*(x_rt[N]);
	y_rr	[N+1] = y_rr [0] + A[N]*dt*(x_rr[N]);
	y		[N+1] = y	 [0] + A[N]*dt*(x	[N]);
	x_et  	[N+1] = x_et [0] + A[N]*dt*(pow(gamma_e, 2) * (                 + N_pt * Cortex->y[N] - y_et[N]) - 2 * gamma_e * x_et[N]) + noise_xRK(N,0);
	x_er  	[N+1] = x_er [0] + A[N]*dt*(pow(gamma_e, 2) * (N_tr * get_Qt(N)	+ N_pr * Cortex->y[N] - y_er[N]) - 2 * gamma_e * x_er[N]);
	x_rt  	[N+1] = x_rt [0] + A[N]*dt*(pow(gamma_g, 2) * (N_rt * get_Qr(N)						  - y_rt[N]) - 2 * gamma_g * x_rt[N]);
	x_rr  	[N+1] = x_rr [0] + A[N]*dt*(pow(gamma_g, 2) * (N_rr * get_Qr(N)						  - y_rr[N]) - 2 * gamma_g * x_rr[N]);
	x	  	[N+1] = x	 [0] + A[N]*dt*(pow(nu,		 2) * (		  get_Qt(N)						  - y   [N]) - 2 * nu	   * x   [N]);
	h_T_t   [N+1] = h_T_t[0] + A[N]*dt*(h_inf_T_t(N) - h_T_t[N])/tau_h_T_t(N);
	h_T_r 	[N+1] = h_T_r[0] + A[N]*dt*(h_inf_T_r(N) - h_T_r[N])/tau_h_T_r(N);
	m_h 	[N+1] = m_h  [0] + A[N]*dt*((m_inf_h(N) * (1 - m_h2[N]) - m_h[N])/tau_m_h(N) - k3 * P_h(N) * m_h[N] + k4 * m_h2[N]);
	m_h2 	[N+1] = m_h2 [0] + A[N]*dt*(k3 * P_h(N) * m_h[N] - k4 * m_h2[N]);
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									Function that adds all SRK terms								*/
/****************************************************************************************************/
void Thalamic_Column::add_RK(void) {
	Vt	  	[0] =(-3*Vt   [0] + 2*Vt   [1] + 4*Vt	[2] + 2* Vt     [3] + Vt	[4])/6;
	Vr	  	[0] =(-3*Vr   [0] + 2*Vr   [1] + 4*Vr	[2] + 2* Vr     [3] + Vr	[4])/6;
	Ca	  	[0] =(-3*Ca   [0] + 2*Ca   [1] + 4*Ca	[2] + 2* Ca     [3] + Ca	[4])/6;
	y_et	[0] =(-3*y_et [0] + 2*y_et [1] + 4*y_et	[2] + 2* y_et	[3] + y_et	[4])/6;
	y_er	[0] =(-3*y_er [0] + 2*y_er [1] + 4*y_er	[2] + 2* y_er	[3] + y_er	[4])/6;
	y_rt	[0] =(-3*y_rt [0] + 2*y_rt [1] + 4*y_rt	[2] + 2* y_rt	[3] + y_rt	[4])/6;
	y_rr	[0] =(-3*y_rr [0] + 2*y_rr [1] + 4*y_rr	[2] + 2* y_rr	[3] + y_rr	[4])/6;
	y	  	[0] =(-3*y	  [0] + 2*y	   [1] + 4*y	[2] + 2* y		[3] + y		[4])/6;
	x_et  	[0] =(-3*x_et [0] + 2*x_et [1] + 4*x_et	[2] + 2* x_et	[3] + x_et	[4])/6 + noise_aRK(0);
	x_er  	[0] =(-3*x_er [0] + 2*x_er [1] + 4*x_er	[2] + 2* x_er	[3] + x_er	[4])/6;
	x_rt  	[0] =(-3*x_rt [0] + 2*x_rt [1] + 4*x_rt	[2] + 2* x_rt	[3] + x_rt	[4])/6;
	x_rr  	[0] =(-3*x_rr [0] + 2*x_rr [1] + 4*x_rr	[2] + 2* x_rr	[3] + x_rr	[4])/6;
	x	  	[0] =(-3*x	  [0] + 2*x	   [1] + 4*x	[2] + 2* x		[3] + x		[4])/6;
	h_T_t 	[0] =(-3*h_T_t[0] + 2*h_T_t[1] + 4*h_T_t[2] + 2* h_T_t  [3] + h_T_t	[4])/6;
	h_T_r 	[0] =(-3*h_T_r[0] + 2*h_T_r[1] + 4*h_T_r[2] + 2* h_T_r  [3] + h_T_r	[4])/6;
	m_h		[0] =(-3*m_h  [0] + 2*m_h  [1] + 4*m_h	[2] + 2* m_h    [3] + m_h	[4])/6;
	m_h2	[0] =(-3*m_h2 [0] + 2*m_h2 [1] + 4*m_h2	[2] + 2* m_h2	[3] + m_h2	[4])/6;
	/* Generate noise for the next iteration */
	for (unsigned i=0; i<Rand_vars.size(); ++i) {
		Rand_vars[i] = MTRands[i]() + input;
	}
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
