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

/************************************************************************************************/
/*								Header file of a cortical module								*/
/************************************************************************************************/
#pragma once
#include <cmath>
#include <vector>
#include "Random_Stream.h"
#include "Thalamic_Column.h"
using std::vector;
class Thalamic_Column;

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
/*								Implementation of the cortical module 								*/
/****************************************************************************************************/
class Cortical_Column {
public:
	/* Constructors */
	Cortical_Column(void)
	{set_RNG();}

	Cortical_Column(double* Param, double* Con)
		:sigma_p 	(Param[0]),	g_KNa	(Param[1]), 	  dphi	(Param[2]),
		  N_pt		(Con[2]),	N_it	(Con[3])
	{set_RNG();}

	/* Connect to the thalamic module */
	void	get_Thalamus(Thalamic_Column& T) {Thalamus = &T;}

	/* Data storage access */
	friend void get_data (int, Cortical_Column&, Thalamic_Column&, vector<double*>);

	/* ODE functions */
	void 	set_RK		(int);
	void 	add_RK	 	(void);

	/* Stimulation protocol access */
	friend class Stim;
	friend class Thalamic_Column;

private:
	/* Declaration of private functions */
	/* Initialize the RNGs */
	void 	set_RNG		(void);

	/* Firing rates */
	double 	get_Qp		(int) const;
	double 	get_Qi		(int) const;

	/* Currents */
	double 	I_ep		(int) const;
	double 	I_ei		(int) const;
	double 	I_gp		(int) const;
	double 	I_gi		(int) const;
	double 	I_L_p		(int) const;
	double 	I_L_g		(int) const;
	double 	I_KNa		(int) const;

	/* Potassium pump */
	double 	Na_pump		(int) const;

	/* Noise functions */
	double 	noise_xRK 	(int,int) const;
	double 	noise_aRK 	(int) const;

	/* Declaration and Initialization of parameters */
	/* Membrane time in ms */
	const double 	tau_p 		= 30;
	const double 	tau_i 		= 30;

	/* Maximum firing rate in ms^-1 */
	const double 	Qp_max		= 30.E-3;
	const double 	Qi_max		= 60.E-3;

	/* Sigmoid threshold in mV */
	const double 	theta_p		= -58.5;
	const double 	theta_i		= -58.5;

	/* Sigmoid gain in mV */
	const double 	sigma_p		= 4;
	const double 	sigma_i		= 6;

	/* Scaling parameter for sigmoidal mapping (dimensionless) */
	const double 	C1          = (M_PI/sqrt(3));

	/* parameters of the firing adaption */
	const double 	alpha_Na	= 2;			/* Sodium influx per spike			in mM ms 	*/
	const double 	tau_Na		= 1.7;			/* Sodium time constant 			in ms 		*/

	const double 	R_pump   	= 0.09;        	/* Na-K pump  constant              in mM/ms 	*/
	const double 	Na_eq    	= 9.5;         	/* Na-eq concentration              in mM 		*/

	/* PSP rise time in ms^-1 */
	const double 	gamma_e		= 70E-3;
	const double 	gamma_g		= 58.6E-3;

	/* Axonal flux time constant */
	const double 	nu			= 120E-3;

	/* Leak weight in aU*/
	const double 	g_L    		= 1.;

	/* Synaptic weight in ms */
	const double 	g_AMPA 		= 1.;
	const double 	g_GABA 		= 1.;

	/* Conductivity */
	/* KNa in mS/cm^2 */
	const double	g_KNa		= 1.33;

	/* Reversal potentials in mV */
	/* Synaptic */
	const double 	E_AMPA  	= 0;
	const double 	E_GABA  	= -70;

	/* Leak */
	const double 	E_L_p 		= -64;
	const double 	E_L_g 		= -64;

	/* Potassium */
	const double 	E_K    		= -100;

	/* Noise parameters in ms^-1 */
	const double 	mphi		= 0E-3;
	const double	dphi		= 20E-1;
	double			input		= 0.0;

	/* Connectivities (dimensionless) */
	const double 	N_pp		= 115;
	const double 	N_ip		= 72;
	const double 	N_pi		= 90;
	const double 	N_ii		= 90;
	const double 	N_pt		= 2.5;
	const double 	N_it		= 2.5;

	/* Pointer to thalamic column */
	Thalamic_Column* Thalamus;

	/* Parameters for SRK4 iteration */
	const vector<double> A = {0.5,  0.5,  1.0, 1.0};
	const vector<double> B = {0.75, 0.75, 0.0, 0.0};

	/* Random number generators */
	vector<random_stream_normal> MTRands;

	/* Container for noise */
	vector<double>	Rand_vars;

	/* Population variables																			*/
	vector<double> 	Vp		= _INIT(E_L_p),		/* excitatory membrane voltage						*/
	Vi		= _INIT(E_L_g),		/* inhibitory membrane voltage						*/
	Na		= _INIT(Na_eq),		/* Na concentration									*/
	s_ep	= _INIT(0.0),		/* PostSP from excitatory to excitatory population	*/
	s_ei	= _INIT(0.0),		/* PostSP from excitatory to inhibitory population	*/
	s_gp	= _INIT(0.0),		/* PostSP from inhibitory to excitatory population	*/
	s_gi	= _INIT(0.0),		/* PostSP from inhibitory to inhibitory population	*/
	y		= _INIT(0.0),		/* axonal flux										*/
	x_ep	= _INIT(0.0),		/* derivative of s_ep								*/
	x_ei	= _INIT(0.0),		/* derivative of s_ei								*/
	x_gp	= _INIT(0.0),		/* derivative of s_gp				 				*/
	x_gi	= _INIT(0.0),		/* derivative of s_gi								*/
	x		= _INIT(0.0);		/* derivative of y									*/
};
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/

