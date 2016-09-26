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
/*							Functions of the thalamic module				  */
/******************************************************************************/
#include "Thalamic_Column.h"

/******************************************************************************/
/*							Initialization of RNG 							  */
/******************************************************************************/
void Thalamic_Column::set_RNG(void) {
    extern const double dt;
    unsigned numRandomVariables = 1;

    MTRands.reserve(2*numRandomVariables);
    Rand_vars.reserve(2*numRandomVariables);
    for (unsigned i=0; i < numRandomVariables; ++i){
        /* Add the RNG for I_{l}*/
        MTRands.push_back(randomStreamNormal(0.0, dphi*dt));

        /* Add the RNG for I_{l,0} */
        MTRands.push_back(randomStreamNormal(0.0, dt));

        /* Get the random number for the first iteration */
        Rand_vars.push_back(MTRands[2*i]());
        Rand_vars.push_back(MTRands[2*i+1]());
    }
}
/******************************************************************************/
/*                          RK noise scaling 								  */
/******************************************************************************/
double Thalamic_Column::noise_xRK(int N, int M) const{
    return gamma_e * gamma_e * (Rand_vars[2*M] + Rand_vars[2*M+1]/std::sqrt(3))*B[N];
}

double Thalamic_Column::noise_aRK(int M) const{
    return gamma_e * gamma_e * (Rand_vars[2*M] - Rand_vars[2*M+1]*std::sqrt(3))/4;
}

/******************************************************************************/
/*                          Firing Rate functions 							  */
/******************************************************************************/
double Thalamic_Column::get_Qt	(int N) const{
    return Qt_max/ (1 + exp(-C1 * (Vt[N] - theta_t) / sigma_t));
}

double Thalamic_Column::get_Qr	(int N) const{
    return Qr_max / (1 + exp(-C1 * (Vr[N]- theta_r) / sigma_r));
}

/******************************************************************************/
/*							Synaptic currents								  */
/******************************************************************************/
/* Excitatory input to TC population */
double Thalamic_Column::I_et	(int N) const{
    return g_AMPA * s_et[N]* (Vt[N]- E_AMPA);
}

/* Inhibitory input to TC population */
double Thalamic_Column::I_gt	(int N) const{
    return g_GABA * s_gt[N]* (Vt[N]- E_GABA);
}
/* Excitatory input to RE population */
double Thalamic_Column::I_er	(int N) const{
    return g_AMPA * s_er[N]* (Vr[N]- E_AMPA);
}

/* Inhibitory input to RE population */
double Thalamic_Column::I_gr	(int N) const{
    return g_GABA * s_gr[N]* (Vr[N]- E_GABA);
}

/******************************************************************************/
/*                          I_T gating functions 							  */
/******************************************************************************/
/* Activation in TC population after Destexhe 1996 */
double Thalamic_Column::m_inf_T_t	(int N) const{
    return 1/(1+exp(-(Vt[N]+59)/6.2));
}

/* Activation in RE population after Destexhe 1996 */
double Thalamic_Column::m_inf_T_r	(int N) const{
    return 1/(1+exp(-(Vr[N]+52)/7.4));
}

/* Deactivation in TC population after Destexhe 1996 */
double Thalamic_Column::h_inf_T_t	(int N) const{
    return 1/(1+exp( (Vt[N]+81)/4));
}

/* Deactivation in RE population after Destexhe 1996 */
double Thalamic_Column::h_inf_T_r	(int N) const{
    return 1/(1+exp( (Vr[N]+80)/5));
}

/* Deactivation time in RE population after Destexhe 1996 */
double Thalamic_Column::tau_h_T_t	(int N) const{
    return (30.8 + (211.4 + exp((Vt[N]+115.2)/5))/(1 + exp((Vt[N]+86)/3.2)))/3.7371928;
}

/* Deactivation time in RE population after Destexhe 1996 */
double Thalamic_Column::tau_h_T_r	(int N) const{
    return (85 + 1/(exp((Vr[N]+48)/4) + exp(-(Vr[N]+407)/50)))/3.7371928;
}

/******************************************************************************/
/*                          I_h gating functions 							  */
/******************************************************************************/
/* Activation in TC population after Destexhe 1993 */
double Thalamic_Column::m_inf_h	(int N) const{
    return 1/(1+exp( (Vt[N]+75)/5.5));
}

/* Activation time for slow components in TC population after Chen2012 */
double Thalamic_Column::tau_m_h	(int N) const{
    return (20 + 1000/(exp((Vt[N]+ 71.5)/14.2) + exp(-(Vt[N]+ 89)/11.6)));
}

/* Instantaneous calcium binding onto messenger protein after Chen2012 */
double Thalamic_Column::P_h	(int N) const{
    //return k1 * pow(Ca[N], n_P)/(k1*pow(Ca[N], n_P)+k2);
    return k1 * Ca[N] * Ca[N] * Ca[N] * Ca[N]/(k1 * Ca[N] * Ca[N] * Ca[N] * Ca[N]+k2);
}

/* Return I_h activation */
double Thalamic_Column::act_h	(void) const{
    return m_h[0] + g_inc * m_h2[0];
}

/******************************************************************************/
/*							Intrinsic currents                                */
/******************************************************************************/
/* Leak current of TC population */
double Thalamic_Column::I_L_t	(int N) const{
    return g_L * (Vt[N]- E_L_t);

}

/* Potassium leak current of TC population */
double Thalamic_Column::I_LK_t	(int N) const{
    return g_LK * (Vt[N]- E_K);
}

/* Leak current of RE population */
double Thalamic_Column::I_L_r	(int N) const{
    return g_L * (Vr[N]- E_L_r);
}

/* Potassium leak current of RE population */
double Thalamic_Column::I_LK_r	(int N) const{
    return g_LK	* (Vr[N]- E_K);
}

/* T-type current of TC population */
double Thalamic_Column::I_T_t	(int N) const{
    return g_T_t * m_inf_T_t(N) * m_inf_T_t(N) * h_T_t[N] * (Vt[N]- E_Ca);
}

/* T-type current of RE population */
double Thalamic_Column::I_T_r	(int N) const{
    return g_T_r * m_inf_T_r(N) * m_inf_T_r(N) * h_T_r[N] * (Vr[N]- E_Ca);
}

/* h-type current of TC population */
double Thalamic_Column::I_h		(int N) const{
    return g_h * (m_h[N] + g_inc * m_h2[N]) * (Vt[N]- E_h);
}

/******************************************************************************/
/*                              SRK iteration                                 */
/******************************************************************************/
void Thalamic_Column::set_RK (int N) {
    extern const double dt;
    Vt	  	[N+1] = Vt   [0] + A[N]*dt*(-(I_L_t(N) + I_et(N) + I_gt(N))/tau_t - C_m * (I_LK_t(N) + I_T_t(N) + I_h(N)));
    Vr	  	[N+1] = Vr   [0] + A[N]*dt*(-(I_L_r(N) + I_er(N) + I_gr(N))/tau_r - C_m * (I_LK_r(N) + I_T_r(N)));
    Ca      [N+1] = Ca   [0] + A[N]*dt*(alpha_Ca * I_T_t(N) - (Ca[N] - Ca_0)/tau_Ca);
    h_T_t   [N+1] = h_T_t[0] + A[N]*dt*(h_inf_T_t(N) - h_T_t[N])/tau_h_T_t(N);
    h_T_r 	[N+1] = h_T_r[0] + A[N]*dt*(h_inf_T_r(N) - h_T_r[N])/tau_h_T_r(N);
    m_h 	[N+1] = m_h  [0] + A[N]*dt*((m_inf_h(N) * (1 - m_h2[N]) - m_h[N])/tau_m_h(N) - k3 * P_h(N) * m_h[N] + k4 * m_h2[N]);
    m_h2 	[N+1] = m_h2 [0] + A[N]*dt*(k3 * P_h(N) * m_h[N] - k4 * m_h2[N]);
    s_et	[N+1] = s_et [0] + A[N]*dt*(x_et[N]);
    s_er	[N+1] = s_er [0] + A[N]*dt*(x_er[N]);
    s_gt	[N+1] = s_gt [0] + A[N]*dt*(x_gt[N]);
    s_gr	[N+1] = s_gr [0] + A[N]*dt*(x_gr[N]);
    y		[N+1] = y	 [0] + A[N]*dt*(x	[N]);
    x_et  	[N+1] = x_et [0] + A[N]*dt*(gamma_e*gamma_e * (                 + N_tp * Cortex->y[N] - s_et[N]) - 2 * gamma_e * x_et[N]) + noise_xRK(N,0);
    x_er  	[N+1] = x_er [0] + A[N]*dt*(gamma_e*gamma_e * (N_rt * get_Qt(N)	+ N_rp * Cortex->y[N] - s_er[N]) - 2 * gamma_e * x_er[N]);
    x_gt  	[N+1] = x_gt [0] + A[N]*dt*(gamma_g*gamma_g * (N_tr * get_Qr(N)						  - s_gt[N]) - 2 * gamma_g * x_gt[N]);
    x_gr  	[N+1] = x_gr [0] + A[N]*dt*(gamma_g*gamma_g * (N_rr * get_Qr(N)						  - s_gr[N]) - 2 * gamma_g * x_gr[N]);
    x	  	[N+1] = x	 [0] + A[N]*dt*(nu * nu         * (		  get_Qt(N)						  - y   [N]) - 2 * nu	   * x   [N]);
}

void Thalamic_Column::add_RK(void) {
    add_RK(Vt);
    add_RK(Vr);
    add_RK(Ca);
    add_RK(s_et);
    add_RK(s_er);
    add_RK(s_gt);
    add_RK(s_gr);
    add_RK(y);
    add_RK_noise(x_et, 0);
    add_RK(x_er);
    add_RK(x_gt);
    add_RK(x_gr);
    add_RK(x);
    add_RK(h_T_t);
    add_RK(h_T_r);
    add_RK(m_h);
    add_RK(m_h2);

    /* Generate noise for the next iteration */
    for (unsigned i=0; i<Rand_vars.size(); ++i) {
        Rand_vars[i] = MTRands[i]() + input;
    }
}
