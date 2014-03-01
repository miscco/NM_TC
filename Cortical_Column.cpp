#include "Cortical_Column.h"
/****************************************************************************************************/
/*										Firing Rate functions										*/
/****************************************************************************************************/
// pyramidal firing rate
double Cortical_Column::get_Qe	(int N) const{
	_SWITCH((Ve))
	double q = Qe_max / (1 + exp(-Scale * (var_Ve - theta_e) / sigma_e));
	return q;
}

// cortical inhibitory firing rate
double Cortical_Column::get_Qi	(int N) const{
	_SWITCH((Vi))
	double q = Qi_max / (1 + exp(-Scale * (var_Vi - theta_i) / sigma_i));
	return q;
}
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/


/****************************************************************************************************/
/*											Synaptic input											*/
/****************************************************************************************************/
// cortical axonal flux
double Cortical_Column::get_phi	(int N) const{
	_SWITCH((phi_e))
	return var_phi_e;
}
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 synaptic currents											*/
/****************************************************************************************************/
// excitatory input to pyramidal population
double Cortical_Column::I_ee	(int N) const{
	_SWITCH((Ve)(Phi_ee))
	double I = var_Phi_ee * (var_Ve - E_AMPA);
	return I;
}

// inhibitory input to pyramidal population
double Cortical_Column::I_ie	(int N) const{
	_SWITCH((Ve)(Phi_ie))
	double I = var_Phi_ie * (var_Ve - E_GABA);
	return I;
}
// excitatory input to inhibitory population
double Cortical_Column::I_ei	(int N) const{
	_SWITCH((Vi)(Phi_ei))
	double I = var_Phi_ei * (var_Vi - E_AMPA);
	return I;
}

// inhibitory input to inhibitory population
double Cortical_Column::I_ii	(int N) const{
	_SWITCH((Vi)(Phi_ii))
	double I = var_Phi_ii * (var_Vi - E_GABA);
	return I;
}
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/



/****************************************************************************************************/
/*										 Current functions											*/
/****************************************************************************************************/
// Leak current of pyramidal population
double Cortical_Column::I_L_e	(int N) const{
	_SWITCH((Ve))
	double I = gL_e * (var_Ve - E_L_e);
	return I;
}

// Leak current of inhibitory population
double Cortical_Column::I_L_i	(int N) const{
	_SWITCH((Vi))
	double I = gL_i * (var_Vi - E_L_i);
	return I;
}

// sodium dependent potassium current
double Cortical_Column::I_KNa		(int N)  const{
	_SWITCH((Ve)(Na))
	double w_KNa  = 0.37/(1+pow(38.7/var_Na, 3.5));
	double I_KNa  = g_KNa * w_KNa * (var_Ve - E_K);
	return I_KNa;
}
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 potassium concentration									*/
/****************************************************************************************************/
double Cortical_Column::Na_pump		(int N) const{
	_SWITCH((Na))
	double Na_pump = R_pump*( pow(var_Na, 3)/(pow(var_Na, 3)+3375)  -  pow(Na_eq, 3)/(pow(Na_eq, 3)+3375));
	return Na_pump;
}
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/



/****************************************************************************************************/
/*										 RK noise scaling 											*/
/****************************************************************************************************/
// function that returns the noise to exitatory population for stochastic RK4
double Cortical_Column::noise_xRK(int N, int M) const{
	extern const double h;
	extern const vector<double> B1, B2;
	double n = 1  / h * (B1[N-1] * Rand_vars[2*M] + B2[N-1] * Rand_vars[2*M+1]);
	return n;
}
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 	ODE functions											*/
/****************************************************************************************************/
// function that calculates the Nth RK term
void Cortical_Column::set_RK		(int N) {
	extern const double dt;
	_SWITCH((Phi_ee)(Phi_ei)(Phi_ie)(Phi_ii)(phi_e)
			(x_ee) 	(x_ei)	(x_ie)	(x_ii)	(y_e))
	Ve	  	[N] = dt*(-(I_L_e(N) + I_ee(N) + I_ie(N))/tau_e - I_KNa(N));
	Vi	  	[N] = dt*(-(I_L_i(N) + I_ei(N) + I_ii(N))/tau_i);
	Na		[N] = dt*(alpha_Na*get_Qe(N) - Na_pump(N))/tau_Na;
	Phi_ee	[N] = dt*(var_x_ee);
	Phi_ei	[N] = dt*(var_x_ei);
	Phi_ie	[N] = dt*(var_x_ie);
	Phi_ii	[N] = dt*(var_x_ii);
	phi_e	[N] = dt*(var_y_e);
	x_ee  	[N] = dt*(pow(gamma_e, 2) * (N_ee * get_Qe(N) + N_te * Thalamus->get_phi(N)	+ noise_xRK(N, 0)	- var_Phi_ee) - 2 * gamma_e * var_x_ee);
	x_ei  	[N] = dt*(pow(gamma_e, 2) * (N_ei * get_Qe(N) + N_ti * Thalamus->get_phi(N)	+ noise_xRK(N, 1)	- var_Phi_ei) - 2 * gamma_e * var_x_ei);
	x_ie  	[N] = dt*(pow(gamma_i, 2) * (N_ie * get_Qi(N) 			  	   									- var_Phi_ie) - 2 * gamma_i * var_x_ie);
	x_ii  	[N] = dt*(pow(gamma_i, 2) * (N_ii * get_Qi(N)		 	  	  									- var_Phi_ii) - 2 * gamma_i * var_x_ii);
	y_e	 	[N] = dt*(pow(nu, 2) 	  * (	    get_Qe(N)			  										- var_phi_e)  - 2 * nu 	  	* var_y_e);
}

// function that ads all the RK terms together
void Cortical_Column::add_RK(void) {
	extern const double h;
	Ve	  	[0] += (Ve		[1] + Ve	[2] * 2 + Ve	[3] * 2 + Ve	[4])/6;
	Vi	  	[0] += (Vi		[1] + Vi	[2] * 2 + Vi	[3] * 2 + Vi	[4])/6;
	Na		[0] += (Na		[1] + Na	[2] * 2 + Na	[3] * 2 + Na	[4])/6;
	Phi_ee	[0] += (Phi_ee	[1] + Phi_ee[2] * 2 + Phi_ee[3] * 2 + Phi_ee[4])/6;
	Phi_ei	[0] += (Phi_ei	[1] + Phi_ei[2] * 2 + Phi_ei[3] * 2 + Phi_ei[4])/6;
	Phi_ie	[0] += (Phi_ie	[1] + Phi_ie[2] * 2 + Phi_ie[3] * 2 + Phi_ie[4])/6;
	Phi_ii	[0] += (Phi_ii	[1] + Phi_ii[2] * 2 + Phi_ii[3] * 2 + Phi_ii[4])/6;
	phi_e	[0] += (phi_e	[1] + phi_e	[2] * 2 + phi_e	[3] * 2 + phi_e	[4])/6;
	x_ee  	[0] += (x_ee	[1] + x_ee	[2] * 2 + x_ee	[3] * 2 + x_ee	[4])/6 + pow(gamma_e, 2) * h * Rand_vars[0];
	x_ei  	[0] += (x_ei	[1] + x_ei	[2] * 2 + x_ei	[3] * 2 + x_ei	[4])/6 + pow(gamma_e, 2) * h * Rand_vars[2];
	x_ie  	[0] += (x_ie	[1] + x_ie	[2] * 2 + x_ie	[3] * 2 + x_ie	[4])/6;
	x_ii  	[0] += (x_ii	[1] + x_ii	[2] * 2 + x_ii	[3] * 2 + x_ii	[4])/6;
	y_e   	[0] += (y_e		[1] + y_e	[2] * 2 + y_e	[3] * 2 + y_e	[4])/6;
	// generating the noise for the next iteration
	for (unsigned i=0; i<Rand_vars.size(); ++i) {
		Rand_vars[i] = MTRands[i]();
	}
}
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/
