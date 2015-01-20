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

/************************************************************************************************/
/*								Header file of a cortical module								*/
/************************************************************************************************/
#pragma once
#include <cmath>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include "Thalamic_Column.h"
#include "macros.h"
using std::vector;

class Thalamic_Column;
/****************************************************************************************************/
/*										Typedefs for RNG											*/
/****************************************************************************************************/
typedef boost::mt11213b                    	ENG;    /* Mersenne Twister		*/
typedef boost::normal_distribution<double>	DIST;   /* Normal Distribution	*/
typedef boost::variate_generator<ENG,DIST> 	GEN;    /* Variate generator	*/
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*								Implementation of the cortical module 								*/
/****************************************************************************************************/
class Cortical_Column {
public:
	/* Constructors */
	Cortical_Column(void)
	{set_RNG();}

	Cortical_Column(double* Param, double* Con)
	 :sigma_e 	(Param[0]),	g_KNa	(Param[1]), 	  dphi	(Param[2]),
	  N_te		(Con[2]),	N_ti	(Con[3])
	{set_RNG();}

	/* Connect to the thalamic module */
	void	get_Thalamus(Thalamic_Column& T) {Thalamus = &T;}

	/* Return axonal flux */
	double 	get_phi		(int N) const {_SWITCH((phi)); return var_phi;}

	/* Data storage access */
	friend void get_data (int, Cortical_Column&, Thalamic_Column&, _REPEAT(double*, 4));

    /* ODE functions */
    void 	set_RK		(int);
    void 	add_RK	 	(void);

	/* Stimulation protocol access */
	friend class 	Stim;

private:
    /* Initialize the RNGs */
    void 	set_RNG		(void);

    /* Firing rates */
    double 	get_Qe		(int) const;
    double 	get_Qi		(int) const;

    /* Currents */
    double 	I_ee		(int) const;
    double 	I_ei		(int) const;
    double 	I_ie		(int) const;
    double 	I_ii		(int) const;
    double 	I_L_e		(int) const;
    double 	I_L_i		(int) const;
    double 	I_KNa		(int) const;

    /* Potassium pump */
    double 	Na_pump		(int) const;

    /* Noise function */
    double 	noise_xRK 	(int, int) const;

    /* Population variables */
	vector<double> 	Ve		= _INIT(E_L_e),		/* excitatory membrane voltage						*/
					Vi		= _INIT(E_L_i),		/* inhibitory membrane voltage						*/
					Na		= _INIT(Na_eq),		/* Na concentration									*/
					Phi_ee	= _INIT(0.0),		/* PostSP from excitatory to excitatory population	*/
					Phi_ei	= _INIT(0.0),		/* PostSP from excitatory to inhibitory population	*/
					Phi_ie	= _INIT(0.0),		/* PostSP from inhibitory to excitatory population	*/
					Phi_ii	= _INIT(0.0),		/* PostSP from inhibitory to inhibitory population	*/
					phi		= _INIT(0.0),		/* axonal flux										*/
					x_ee	= _INIT(0.0),		/* derivative of Phi_ee								*/
					x_ei	= _INIT(0.0),		/* derivative of Phi_ei								*/
					x_ie	= _INIT(0.0),		/* derivative of Phi_ie				 				*/
					x_ii	= _INIT(0.0),		/* derivative of Phi_ii 							*/
					y		= _INIT(0.0);		/* derivative of phi 								*/

	/* Random number generators */
	vector<GEN>		MTRands;

	/* Container for noise */
	vector<double>	Rand_vars;

	/* Declaration and Initialization of parameters */
	/* Membrane time in ms */
	const double 	tau_e 		= 30;
	const double 	tau_i 		= 30;

	/* Maximum firing rate in ms^-1 */
	const double 	Qe_max		= 30.E-3;
	const double 	Qi_max		= 60.E-3;

	/* Sigmoid threshold in mV */
	const double 	theta_e		= -58.5;
	const double 	theta_i		= -58.5;

	/* Sigmoid gain in mV */
	const double 	sigma_e		= 4;
	const double 	sigma_i		= 6;

	/* Scaling parameter for sigmoidal mapping (dimensionless) */
	const double 	C1          = (3.14159265/sqrt(3));

	/* parameters of the firing adaption */
	const double 	alpha_Na	= 2;			/* Sodium influx per spike			in mM ms 	*/
    const double 	tau_Na		= 1.7;			/* Sodium time constant 			in ms 		*/

    const double 	R_pump   	= 0.09;        	/* Na-K pump  constant              in mM/ms 	*/
	const double 	Na_eq    	= 9.5;         	/* Na-eq concentration              in mM 		*/

	/* PSP rise time in ms^-1 */
	const double 	gamma_e		= 70E-3;
	const double 	gamma_i		= 58.6E-3;

	/* Axonal flux time constant */
    const double 	nu			= 60E-3;

	/* Conductivities in mS/cm^-2 */
	/* Leak */
	const double 	g_L    		= 1;

	/* KNa */
	const double	g_KNa		= 1.33;

	/* Reversal potentials in mV */
	/* Synaptic */
	const double 	E_AMPA  	= 0;
	const double 	E_GABA  	= -70;

	/* Leak */
    const double 	E_L_e 		= -64;
	const double 	E_L_i 		= -64;

	/* Potassium */
	const double 	E_K    		= -100;

	/* Noise parameters in ms^-1 */
	const double 	mphi		= 0E-3;
    const double	dphi		= 60E-3;
	double			input		= 0.0;

	/* Connectivities (dimensionless) */
    const double 	N_ee		= 115;
	const double 	N_ei		= 72;
    const double 	N_ie		= 90;
    const double 	N_ii		= 90;
	const double 	N_te		= 0;
	const double 	N_ti		= 0;

	/* Pointer to thalamic column */
	Thalamic_Column* Thalamus;
};
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/

