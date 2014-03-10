/****************************************************************************************************/
/*									Functions of the cortical module								*/
/****************************************************************************************************/
#include "Thalamic_Column.h"

/****************************************************************************************************/
/*										 Initialization of RNG 										*/
/****************************************************************************************************/
void Thalamic_Column::set_RNG(void) {
	// the number of independent streams
	int N = 2;

	// get the RNG
	for (int i=0; i<N; ++i){
		// add the RNG
		MTRands.push_back({ENG(rand()), DIST (mphi_t, dphi_t)});

		// get the random number for the first iteration
		Rand_vars.push_back(MTRands[i]());
	}
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Firing Rate functions 										*/
/****************************************************************************************************/
// thalamic relay (TC) firing rate
double Thalamic_Column::get_Qt	(int N) const{
	_SWITCH((Vt))
	double q = Qt_max/ (1 + exp(-C1 * (var_Vt - theta_t) / sigma_t));
	return q;
}

// thalamic reticular (RE) firing rate
double Thalamic_Column::get_Qr	(int N) const{
	_SWITCH((Vr))
	double q = Qr_max / (1 + exp(-C1 * (var_Vr - theta_r) / sigma_r));
	return q;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										Synaptic currents											*/
/****************************************************************************************************/
// excitatory input to TC population
double Thalamic_Column::I_et	(int N) const{
	_SWITCH((Vt)(Phi_tt))
	double psi = g_AMPA * var_Phi_tt * (var_Vt - E_AMPA);
	return psi;
}

// inhibitory input to TC population
double Thalamic_Column::I_it	(int N) const{
	_SWITCH((Vt)(Phi_rt))
	double psi = g_GABA * var_Phi_rt * (var_Vt - E_GABA);
	return psi;
}
// excitatory input to RE population
double Thalamic_Column::I_er	(int N) const{
	_SWITCH((Vr)(Phi_tr))
	double psi = g_AMPA * var_Phi_tr * (var_Vr - E_AMPA);
	return psi;
}

// inhibitory input to RE population
double Thalamic_Column::I_ir	(int N) const{
	_SWITCH((Vr)(Phi_rr))
	double psi = g_GABA * var_Phi_rr * (var_Vr - E_GABA);
	return psi;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 	I_T gating	 											*/
/****************************************************************************************************/
// instant activation in TC population
// after Destexhe 1996
double Thalamic_Column::m_inf_T_t	(int N) const{
	_SWITCH((Vt))
	double m = 1/(1+exp(-(var_Vt+59)/6.2));
	return m;
}

// instant activation in RE population
// after Destexhe 1996
double Thalamic_Column::m_inf_T_r	(int N) const{
	_SWITCH((Vr))
	double m = 1/(1+exp(-(var_Vr+52)/7.2));
	return m;
}

// instant deactivation in TC population
// after Destexhe 1996
double Thalamic_Column::h_inf_T_t	(int N) const{
	_SWITCH((Vt))
	double h = 1/(1+exp( (var_Vt+81)/4));
	return h;
}

// instant deactivation in RE population
// after Destexhe 1996
double Thalamic_Column::h_inf_T_r	(int N) const{
	_SWITCH((Vr))
	double h = 1/(1+exp( (var_Vr+80)/5));
	return h;
}

// activation time in TC population
// after Bazhenov 1998
double Thalamic_Column::tau_m_T_t	(int N) const{
	_SWITCH((Vt))
	double tau = (0.612 + 1 /(exp(-(var_Vt + 132) / 16.7) + exp((var_Vt + 16.8)/18.2)))/pow(3,1.2);
	return tau;
}

// activation time in RE population
// after Destexhe 1996
double Thalamic_Column::tau_m_T_r	(int N) const{
	_SWITCH((Vr))
	double tau = (1 + 0.33/( exp((var_Vr+27)/10.0) + exp(-(var_Vr+102)/15.0)))/pow(2.5, 1.2);
	return tau;
}

// deactivation time in RE population
// after Destexhe 1996
double Thalamic_Column::tau_h_T_t	(int N) const{
	_SWITCH((Vt))
	double tau =  (30.8 + (211.4 + exp((var_Vt+115.2)/5))/(1 + exp((var_Vt+86)/3.2)))/pow(3, 1.2);
	return tau;
}

// deactivation time in RE population
// after Destexhe 1996
double Thalamic_Column::tau_h_T_r	(int N) const{
	_SWITCH((Vr))
	double tau =  (85 + 1/(exp((var_Vr+48)/4) + exp(-(var_Vr+407)/50)))/pow(3, 1.2);
	return tau;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 	I_h gating 												*/
/****************************************************************************************************/
// instant activation in TC population
// after Destexhe 1993
double Thalamic_Column::m_inf_h	(int N) const{
	_SWITCH((Vt))
	double h = 1/(1+exp( (var_Vt+75)/5.5));
	return h;
}

// activation time for slow components in TC population
// after Destexhe 1993
double Thalamic_Column::tau_m_h	(int N) const{
	_SWITCH((Vt))
	double tau = 1. / (exp(-14.59 - 0.086 * var_Vt) + exp(-1.87 + 0.07 * var_Vt));
	return tau;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Current functions 											*/
/****************************************************************************************************/
// Leak current of TC population
double Thalamic_Column::I_L_t	(int N) const{
	_SWITCH((Vt))
	double I = gL_t * (var_Vt - E_L_t);
	return I;
}

// Potassium leak current of TC population
double Thalamic_Column::I_LK_t	(int N) const{
	_SWITCH((Vt))
	double I = gLK_t * (var_Vt - E_LK_t);
	return I;
}

// Leak current of RE population
double Thalamic_Column::I_L_r	(int N) const{
	_SWITCH((Vr))
	double I = gL_r * (var_Vr - E_L_r);
	return I;
}

// Potassium leak current of RE population
double Thalamic_Column::I_LK_r	(int N) const{
	_SWITCH((Vr))
	double I = gLK_r * (var_Vr - E_LK_r);
	return I;
}

// T-type current of TC population
double Thalamic_Column::I_T_t	(int N) const{
	_SWITCH((Vt)(h_T_t))
	//double I = gTt * pow(var_m_T_t, 2) * var_h_T_t * (var_Vt - E_Ca);
	double I = gTt * pow(m_inf_T_t(N), 2) * var_h_T_t * (var_Vt - E_Ca);
	return I;
}

// T-type current of RE population
double Thalamic_Column::I_T_r	(int N) const{
	_SWITCH((Vr)(h_T_r)(m_T_r))
	double I = gTr * pow(var_m_T_r, 2) * var_h_T_r * (var_Vr - E_Ca);
	//double I = gTr * pow(m_inf_T_r(N), 2) * var_h_T_r * (var_Vr - E_Ca);
	return I;
}

// h-type current of TC population
double Thalamic_Column::I_h		(int N) const{
	_SWITCH((Vt)(m_h)(m_h2))
	double I = gh * (var_m_h + g_inc * var_m_h2) * (var_Vt - E_h);
	return I;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 RK noise scaling 											*/
/****************************************************************************************************/
double Thalamic_Column::noise_xRK(int N) const{
	extern const double h;
	extern const vector<double> B1, B2;
	double n = 1  / h * (B1[N-1] * Rand_vars[0] + B2[N-1] * Rand_vars[1]);
	return n;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										Calculate the Nth SRK term									*/
/****************************************************************************************************/
void Thalamic_Column::set_RK (int N) {
	extern const double dt;
	_SWITCH((Ca)
			(Phi_tt)(Phi_tr)(Phi_rt)(Phi_rr)(phi_t)
			(x_tt)	(x_tr)	(x_rt)	(x_rr)	(y_t)
			(m_T_t)	(m_T_r)	(h_T_t)	(h_T_r)
			(m_h)	(m_h2)	(P_h))
	Vt	  	[N] = dt*(-(I_L_t(N) + I_et(N) + I_it(N))/tau_t - (I_LK_t(N) + I_T_t(N) + I_h(N)));
	Vr	  	[N] = dt*(-(I_L_r(N) + I_er(N) + I_ir(N))/tau_r - (I_LK_r(N) + I_T_r(N)));
	Ca		[N] = dt*(alpha_Ca * I_T_t(N) - (var_Ca - Ca_0)/tau_Ca);
	m_T_t 	[N] = dt*(m_inf_T_t(N) - var_m_T_t)/tau_m_T_t(N);
	m_T_r 	[N] = dt*(m_inf_T_r(N) - var_m_T_r)/tau_m_T_r(N);
	h_T_t 	[N] = dt*(h_inf_T_t(N) - var_h_T_t)/tau_h_T_t(N);
	h_T_r 	[N] = dt*(h_inf_T_r(N) - var_h_T_r)/tau_h_T_r(N);
	m_h 	[N] = dt*((m_inf_h(N) * (1 - var_m_h2) - var_m_h)/tau_m_h(N) - k3 * var_P_h * var_m_h + k4 * var_m_h2);
	m_h2 	[N] = dt*(k3   * var_P_h			 * 	    var_m_h  	- k4   * var_m_h2);
	P_h		[N] = dt*(k1   * pow(var_Ca, n_P) 	 * (1 - var_P_h) 	- k2   * var_P_h);
	Phi_tt	[N] = dt*(var_x_tt);
	Phi_tr	[N] = dt*(var_x_tr);
	Phi_rt	[N] = dt*(var_x_rt);
	Phi_rr	[N] = dt*(var_x_rr);
	phi_t	[N] = dt*(var_y_t);
	x_tt  	[N] = dt*(pow(gamma_e, 2) * (noise_xRK(N) 		+ N_et * Cortex->get_phi(N) 	- var_Phi_tt) - 2 * gamma_e * var_x_tt);
	x_tr  	[N] = dt*(pow(gamma_e, 2) * (N_tr * get_Qt(N)	+ N_er * Cortex->get_phi(N)	- var_Phi_tr) - 2 * gamma_e * var_x_tr);
	x_rt  	[N] = dt*(pow(gamma_i, 2) * (N_rt * get_Qr(N) 									- var_Phi_rt) - 2 * gamma_i * var_x_rt);
	x_rr  	[N] = dt*(pow(gamma_i, 2) * (N_rr * get_Qr(N)									- var_Phi_rr) - 2 * gamma_i * var_x_rr);
	y_t  	[N] = dt*(pow(nu, 	   2) * (		get_Qt(N)									- var_phi_t)  - 2 * nu	 	* var_y_t);
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									Function that adds all SRK terms								*/
/****************************************************************************************************/
void Thalamic_Column::add_RK(void) {
	extern const double h;
	Vt	  	[0] += (Vt		[1] + Vt		[2] * 2 + Vt		[3] * 2 + Vt		[4])/6;
	Vr	  	[0] += (Vr		[1] + Vr		[2] * 2 + Vr		[3] * 2 + Vr		[4])/6;
	Ca	  	[0] += (Ca		[1] + Ca		[2] * 2 + Ca		[3] * 2 + Ca		[4])/6;
	Phi_tt	[0] += (Phi_tt	[1] + Phi_tt	[2] * 2 + Phi_tt	[3] * 2 + Phi_tt	[4])/6;
	Phi_tr	[0] += (Phi_tr	[1] + Phi_tr	[2] * 2 + Phi_tr	[3] * 2 + Phi_tr	[4])/6;
	Phi_rt	[0] += (Phi_rt	[1] + Phi_rt	[2] * 2 + Phi_rt	[3] * 2 + Phi_rt	[4])/6;
	Phi_rr	[0] += (Phi_rr	[1] + Phi_rr	[2] * 2 + Phi_rr	[3] * 2 + Phi_rr	[4])/6;
	phi_t	[0] += (phi_t	[1] + phi_t		[2] * 2 + phi_t		[3] * 2 + phi_t		[4])/6;
	x_tt  	[0] += (x_tt	[1] + x_tt		[2] * 2 + x_tt		[3] * 2 + x_tt		[4])/6 + pow(gamma_e, 2) * h * Rand_vars[0];
	x_tr  	[0] += (x_tr	[1] + x_tr		[2] * 2 + x_tr		[3] * 2 + x_tr		[4])/6;
	x_rt  	[0] += (x_rt	[1] + x_rt		[2] * 2 + x_rt		[3] * 2 + x_rt		[4])/6;
	x_rr  	[0] += (x_rr	[1] + x_rr		[2] * 2 + x_rr		[3] * 2 + x_rr		[4])/6;
	y_t  	[0] += (y_t		[1] + y_t		[2] * 2 + y_t		[3] * 2 + y_t		[4])/6;
	m_T_t 	[0] += (m_T_t	[1] + m_T_t		[2] * 2 + m_T_t		[3] * 2 + m_T_t		[4])/6;
	m_T_r 	[0] += (m_T_r	[1] + m_T_r		[2] * 2 + m_T_r		[3] * 2 + m_T_r		[4])/6;
	h_T_t 	[0] += (h_T_t	[1] + h_T_t		[2] * 2 + h_T_t		[3] * 2 + h_T_t		[4])/6;
	h_T_r 	[0] += (h_T_r	[1] + h_T_r		[2] * 2 + h_T_r		[3] * 2 + h_T_r		[4])/6;
	m_h		[0] += (m_h		[1] + m_h		[2] * 2 + m_h		[3] * 2 + m_h		[4])/6;
	m_h2	[0] += (m_h2	[1] + m_h2		[2] * 2 + m_h2		[3] * 2 + m_h2		[4])/6;
	P_h	 	[0] += (P_h		[1] + P_h		[2] * 2 + P_h		[3] * 2 + P_h		[4])/6;

	// generating the noise for the next iteration
	for (unsigned i=0; i<Rand_vars.size(); ++i) {
		Rand_vars[i] = MTRands[i]() + input;
	}
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
