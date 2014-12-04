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
/*									Functions of the cortical module								*/
/****************************************************************************************************/
#include "Cortical_Column.h"

/****************************************************************************************************/
/*										 Initialization of RNG 										*/
/****************************************************************************************************/
void Cortical_Column::set_RNG(void) {
	/* Number of independent streams */
	int N = 4;

	/* Create RNG for each stream */
	for (int i=0; i<N; ++i){
		/* Add the RNG */
		MTRands.push_back({ENG(rand()), DIST (mphi, sqrt(dphi))});

		/* Get the random number for the first iteration */
		Rand_vars.push_back(MTRands[i]());
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
	_SWITCH((Ve))
    double q = Qe_max / (1 + exp(-C1 * (var_Ve - theta_e) / sigma_e));
	return q;
}

/* Inhibitory firing rate */
double Cortical_Column::get_Qi	(int N) const{
	_SWITCH((Vi))
    double q = Qi_max / (1 + exp(-C1 * (var_Vi - theta_i) / sigma_i));
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
	_SWITCH((Ve)(Phi_ee))
	double I = var_Phi_ee * (var_Ve - E_AMPA);
	return I;
}

/* Inhibitory input to pyramidal population */
double Cortical_Column::I_ie	(int N) const{
	_SWITCH((Ve)(Phi_ie))
	double I = var_Phi_ie * (var_Ve - E_GABA);
	return I;
}
/* Excitatory input to inhibitory population */
double Cortical_Column::I_ei	(int N) const{
	_SWITCH((Vi)(Phi_ei))
	double I = var_Phi_ei * (var_Vi - E_AMPA);
	return I;
}

/* Inhibitory input to inhibitory population */
double Cortical_Column::I_ii	(int N) const{
	_SWITCH((Vi)(Phi_ii))
	double I = var_Phi_ii * (var_Vi - E_GABA);
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
	_SWITCH((Ve))
	double I = g_L * (var_Ve - E_L_e);
	return I;
}

/* Leak current of inhibitory population */
double Cortical_Column::I_L_i	(int N) const{
	_SWITCH((Vi))
	double I = g_L * (var_Vi - E_L_i);
	return I;
}

/* Sodium dependent potassium current */
double Cortical_Column::I_KNa		(int N)  const{
	_SWITCH((Ve)(Na))
    double w_KNa  = 0.37/(1+pow(38.7/var_Na, 3.5));
	double I_KNa  = g_KNa * w_KNa * (var_Ve - E_K);
	return I_KNa;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									 		Potassium pump	 										*/
/****************************************************************************************************/
double Cortical_Column::Na_pump		(int N) const{
	_SWITCH((Na))
	double Na_pump = R_pump*( pow(var_Na, 3)/(pow(var_Na, 3)+3375)  -  pow(Na_eq, 3)/(pow(Na_eq, 3)+3375));
	return Na_pump;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 RK noise scaling 											*/
/****************************************************************************************************/
double Cortical_Column::noise_xRK(int N, int M) const{
	extern const double h;
	extern const vector<double> B1, B2;
	double n = 1  / h * (B1[N-1] * Rand_vars[2*M] + B2[N-1] * Rand_vars[2*M+1]);
	return n;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										Calculate the Nth SRK term									*/
/****************************************************************************************************/
void Cortical_Column::set_RK (int N) {
	extern const double dt;
	_SWITCH((Phi_ee)(Phi_ei)(Phi_ie)(Phi_ii)(phi)
			(x_ee) 	(x_ei)	(x_ie)	(x_ii)	(y))
	Ve	  	[N] = dt*(-(I_L_e(N) + I_ee(N) + I_ie(N))/tau_e - I_KNa(N));
	Vi	  	[N] = dt*(-(I_L_i(N) + I_ei(N) + I_ii(N))/tau_i);
	Na		[N] = dt*(alpha_Na * get_Qe(N) - Na_pump(N))/tau_Na;
	Phi_ee	[N] = dt*(var_x_ee);
	Phi_ei	[N] = dt*(var_x_ei);
	Phi_ie	[N] = dt*(var_x_ie);
	Phi_ii	[N] = dt*(var_x_ii);
	phi		[N] = dt*(var_y);
	x_ee  	[N] = dt*(pow(gamma_e, 2) * (N_ee * get_Qe(N) + noise_xRK(N, 0)	+ N_te * Thalamus->get_phi(N)	- var_Phi_ee) - 2 * gamma_e * var_x_ee);
	x_ei  	[N] = dt*(pow(gamma_e, 2) * (N_ei * get_Qe(N) + noise_xRK(N, 1)	+ N_ti * Thalamus->get_phi(N)	- var_Phi_ei) - 2 * gamma_e * var_x_ei);
	x_ie  	[N] = dt*(pow(gamma_i, 2) * (N_ie * get_Qi(N) 			  										- var_Phi_ie) - 2 * gamma_i * var_x_ie);
	x_ii  	[N] = dt*(pow(gamma_i, 2) * (N_ii * get_Qi(N)		 	  										- var_Phi_ii) - 2 * gamma_i * var_x_ii);
	y	  	[N] = dt*(pow(nu, 	   2) * (		get_Qe(N)													- var_phi)	  - 2 * nu	 	* var_y);
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									Function that adds all SRK terms								*/
/****************************************************************************************************/
void Cortical_Column::add_RK(void) {
	extern const double h;
	Ve	  	[0] += (Ve		[1] + Ve	[2] * 2 + Ve	[3] * 2 + Ve	[4])/6;
	Vi	  	[0] += (Vi		[1] + Vi	[2] * 2 + Vi	[3] * 2 + Vi	[4])/6;
	Na		[0] += (Na		[1] + Na	[2] * 2 + Na	[3] * 2 + Na	[4])/6;
	Phi_ee	[0] += (Phi_ee	[1] + Phi_ee[2] * 2 + Phi_ee[3] * 2 + Phi_ee[4])/6;
	Phi_ei	[0] += (Phi_ei	[1] + Phi_ei[2] * 2 + Phi_ei[3] * 2 + Phi_ei[4])/6;
	Phi_ie	[0] += (Phi_ie	[1] + Phi_ie[2] * 2 + Phi_ie[3] * 2 + Phi_ie[4])/6;
	Phi_ii	[0] += (Phi_ii	[1] + Phi_ii[2] * 2 + Phi_ii[3] * 2 + Phi_ii[4])/6;
	phi		[0] += (phi		[1] + phi	[2] * 2 + phi	[3] * 2 + phi	[4])/6;
	x_ee  	[0] += (x_ee	[1] + x_ee	[2] * 2 + x_ee	[3] * 2 + x_ee	[4])/6 + pow(gamma_e, 2) * h * Rand_vars[0];
	x_ei  	[0] += (x_ei	[1] + x_ei	[2] * 2 + x_ei	[3] * 2 + x_ei	[4])/6 + pow(gamma_e, 2) * h * Rand_vars[2];
	x_ie  	[0] += (x_ie	[1] + x_ie	[2] * 2 + x_ie	[3] * 2 + x_ie	[4])/6;
	x_ii  	[0] += (x_ii	[1] + x_ii	[2] * 2 + x_ii	[3] * 2 + x_ii	[4])/6;
	y	 	[0] += (y		[1] + y		[2] * 2 + y		[3] * 2 + y		[4])/6;

	/* Generat noise for the next iteration */
	for (unsigned i=0; i<Rand_vars.size(); ++i) {
		Rand_vars[i] = MTRands[i]();
	}
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
