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

/************************************************************************************************/
/*								Header file of a thalamic module								*/
/************************************************************************************************/
#pragma once
#include <cmath>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include "Cortical_Column.h"
using std::vector;

class Cortical_Column;

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
/*									Macro for vector initialization									*/
/****************************************************************************************************/
#ifndef _INIT
#define _INIT(x)	{x, 0.0, 0.0, 0.0, 0.0}
#endif
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*								Implementation of the thalamic module 								*/
/****************************************************************************************************/
class Thalamic_Column {
public:
	/* Constructors for compile check */
	Thalamic_Column(void)
	{set_RNG();}

	/* Constructor for simulation */
	Thalamic_Column(double* Param, double* Con)
	: g_LK	(Param[1]),	g_h (Param[0]),
	  N_et 		(Con[0]),	N_er	(Con[1])
	{set_RNG();}

	/* Get the pointer to the cortical module */
	void	get_Cortex	(Cortical_Column& C) {Cortex = &C;}

	/* ODE functions */
	void 	set_RK		(int);
	void 	add_RK	 	(void);

	/* Data storage  access */
	friend void get_data (int, Cortical_Column&, Thalamic_Column&, vector<double*>);
	friend class Cortical_Column;

	/* Set strength of input */
	void	set_input	(double I) {input = I;}

private:
	/* Initialize the RNGs */
	void 	set_RNG		(void);

	/* Firing rates */
	double 	get_Qt		(int) const;
	double 	get_Qr		(int) const;

	/* Synaptic currents */
	double 	I_et		(int) const;
	double 	I_it		(int) const;
	double 	I_er		(int) const;
	double 	I_ir		(int) const;

	/* Activation functions */
	double  m_inf_T_t	(int) const;
	double  m_inf_T_r	(int) const;
	double  m_inf_h		(int) const;
	double  tau_m_h		(int) const;
	double  P_h			(int) const;
	double  act_h		(void)const;

	double  m_inf_hs	(int) const;
	double  tau_m_hs	(int) const;
	double  tau_m_hf	(int) const;

	/* Deactivation functions */
	double  h_inf_T_t	(int) const;
	double  h_inf_T_r	(int) const;
	double  tau_h_T_t	(int) const;
	double  tau_h_T_r	(int) const;

	/* Currents */
	double 	I_L_t		(int) const;
	double 	I_L_r		(int) const;
	double 	I_LK_t		(int) const;
	double 	I_LK_r		(int) const;
	double 	I_T_t		(int) const;
	double 	I_T_r		(int) const;
	double 	I_h			(int) const;

	/* Noise functions */
	double 	noise_xRK 	(int,int) const;
	double 	noise_aRK 	(int) const;

	/* Population variables */
	vector<double> 	Vt		= _INIT(E_L_t),		/* TC membrane voltage								*/
					Vr		= _INIT(E_L_r),		/* RE membrane voltage								*/
					Ca		= _INIT(Ca_0),		/* Calcium concentration of TC population			*/
					y_tt	= _INIT(0.0),		/* PostSP from TC population to TC population		*/
					y_tr	= _INIT(0.0),		/* PostSP from TC population to RE population		*/
					y_rt	= _INIT(0.0),		/* PostSP from RE population to TC population		*/
					y_rr	= _INIT(0.0),		/* PostSP from RE population to RE population		*/
					y		= _INIT(0.0),		/* axonal flux										*/
					x_tt	= _INIT(0.0),		/* derivative of y_tt								*/
					x_tr	= _INIT(0.0),		/* derivative of y_tr								*/
					x_rt	= _INIT(0.0),		/* derivative of y_rt								*/
					x_rr	= _INIT(0.0),		/* derivative of y_rr								*/
					x		= _INIT(0.0),		/* derivative of y									*/
					h_T_t	= _INIT(0.0),		/* inactivation of T channel						*/
					h_T_r	= _INIT(0.0),		/* inactivation of T channel						*/
					m_h		= _INIT(0.0),		/* activation 	of h   channel						*/
					m_h2	= _INIT(0.0);		/* activation 	of h   channel bound with protein 	*/

	/* Random number generators */
	vector<GEN>		MTRands;

	/* Container for noise */
	vector<double>	Rand_vars;

	/* Declaration and Initialization of parameters */
	/* Membrane time in ms */
	const double 	tau_t 		= 20;
	const double 	tau_r 		= 20;

	/* Maximum firing rate in ms^-1 */
	const double 	Qt_max		= 400.E-3;
	const double 	Qr_max		= 400.E-3;

	/* Sigmoid threshold in mV */
	const double 	theta_t		= -58.5;
	const double 	theta_r		= -58.5;

	/* Sigmoid gain in mV */
	const double 	sigma_t		= 6;
	const double 	sigma_r		= 6;

	/* Scaling parameter for sigmoidal mapping (dimensionless) */
	const double 	C1          = (3.14159265/sqrt(3));

	/* PSP rise time in ms^-1 */
	const double 	gamma_e		= 70E-3;
	const double 	gamma_i		= 100E-3;

	/* axonal flux time constant */
	const double 	nu			= 120E-3;

	/* Conductivities in mS/cm^-2 */
	/* Leak current */
	const double 	g_L  		= 1;

	/* Potassium leak current */
	const double 	g_LK 		= 0.02;

	/* T current */
	const double	g_T_t		= 3;
	const double	g_T_r		= 2.3;

	/* h current */
	const double	g_h			= 0.051;

	/* Reversal potentials in mV */
	/* Synaptic */
	const double 	E_AMPA  	= 0;
	const double 	E_GABA  	= -70;

	/* Leak */
	const double 	E_L_t 		= -70;
	const double 	E_L_r 		= -70;

	/* Potassium */
	const double 	E_K    		= -100;

	/* I_T current */
	const double 	E_Ca    	= 120;

	/* I_h current */
	const double 	E_h    		= -40;

	/* Calcium parameters */
	const double	alpha_Ca	= -51.8E-6;			/* influx per spike in nmol		*/
	const double	tau_Ca		= 10;				/* calcium time constant in ms	*/
	const double	Ca_0		= 2.4E-4;			/* resting concentration 		*/

	/* I_h activation parameters */
	const double 	k1			= 2.5E7;
	const double 	k2			= 4E-4;
	const double 	k3			= 1E-1;
	const double 	k4			= 1E-3;
	const double 	n_P			= 4;
	const double 	g_inc		= 2;

	/* Noise parameters in ms^-1 */
	const double 	mphi		= 0E-3;
	const double	dphi		= 20E-1;
	double			input		= 0.0;

	/* Connectivities (dimensionless) */
	const double 	N_tr		= 3;
	const double 	N_rt		= 5;
	const double 	N_rr		= 19;
	const double 	N_et		= 2.6;
	const double 	N_er		= 2.6;

	/* Pointer to cortical column */
	Cortical_Column* Cortex;

	/* SRK integration parameters */
	const vector<double> A = {0.5, 0.5, 1.0, 1.0};
	const vector<double> B = {0.75, 0.75, 0.0, 0.0};
};
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
