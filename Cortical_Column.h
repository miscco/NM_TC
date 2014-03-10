/************************************************************************************************/
/*								Header file of a cortical module								*/
/************************************************************************************************/
#pragma once
#include <cmath>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include "macros.h"
#include "parameters.h"
#include "Thalamic_Column.h"
using std::vector;
class Thalamic_Column;

/****************************************************************************************************/
/*										Typedefs for RNG											*/
/****************************************************************************************************/
typedef boost::mt11213b                    	ENG;    // Mersenne Twister
typedef boost::normal_distribution<double>	DIST;   // Normal Distribution
typedef boost::variate_generator<ENG,DIST> 	GEN;    // Variate generator
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*								Implementation of the cortical module 								*/
/****************************************************************************************************/
class Cortical_Column {
public:
	// Constructors
	Cortical_Column(void)
	: Ve		(_INIT(E_L_e)),	Vi 	   	(_INIT(E_L_i)),	Na	 	(_INIT(Na_eq)),
	  Phi_ee	(_INIT(0.0)), 	Phi_ei 	(_INIT(0.0)), 	Phi_ie 	(_INIT(0.0)), 	Phi_ii	(_INIT(0.0)), 	phi_e	(_INIT(0.0)),
	  x_ee 		(_INIT(0.0)), 	x_ei   	(_INIT(0.0)),	x_ie   	(_INIT(0.0)), 	x_ii	(_INIT(0.0)), 	y_e		(_INIT(0.0)),
	  alpha_Na 	(0), 			tau_Na	(0),			g_KNa	(0),	  		theta_e	(0),
	  sigma_e 	(0), 			dphi_c	(0)
	{set_RNG();}

	Cortical_Column(double* Par)
	: Ve		(_INIT(E_L_e)),	Vi 	   	(_INIT(E_L_i)),	Na	 	(_INIT(Na_eq)),
	  Phi_ee	(_INIT(0.0)), 	Phi_ei 	(_INIT(0.0)), 	Phi_ie 	(_INIT(0.0)), 	Phi_ii	(_INIT(0.0)), 	phi_e	(_INIT(0.0)),
	  x_ee 		(_INIT(0.0)), 	x_ei   	(_INIT(0.0)),	x_ie   	(_INIT(0.0)), 	x_ii	(_INIT(0.0)), 	y_e		(_INIT(0.0)),
	  alpha_Na 	(Par[0]), 		tau_Na	(Par[1]),		g_KNa	(Par[2]),	  	theta_e	(Par[3]),
	  sigma_e 	(Par[4]), 		dphi_c	(Par[5])
	{set_RNG();}

	// get the pointer to the cortical module
	void	get_Thalamus(Thalamic_Column& T) {Thalamus = &T;}

	// Initialize the RNGs
	void 	set_RNG		(void);

	// Firing rates
	double 	get_Qe		(int) const;
	double 	get_Qi		(int) const;

	// Currents
	double 	I_ee		(int) const;
	double 	I_ei		(int) const;
	double 	I_ie		(int) const;
	double 	I_ii		(int) const;
	double 	I_L_e		(int) const;
	double 	I_L_i		(int) const;
	double 	I_KNa		(int) const;

	// Potassium pump
	double 	Na_pump		(int) const;

	// Noise function
	double 	noise_xRK 	(int, int) const;

	// ODE functions
	void 	set_RK		(int);
	void 	add_RK	 	(void);

	// Data storage
	friend void get_data (int, Cortical_Column&, Thalamic_Column&, _REPEAT(double*, 2));

private:
	// Population variables
	vector<double> 	Ve,			// excitatory membrane voltage
					Vi,			// inhibitory membrane voltage
					Na,			// Na concentration
					Phi_ee,		// PostSP from excitatory to excitatory population
					Phi_ei,		// PostSP from excitatory to inhibitory population
					Phi_ie,		// PostSP from inhibitory to excitatory population
					Phi_ii,		// PostSP from inhibitory to inhibitory population
					phi_e,		// axonal flux
					x_ee,		// derivative of Phi_ee
					x_ei,		// derivative of Phi_ei
					x_ie,		// derivative of Phi_ie
					x_ii,		// derivative of Phi_ii
					y_e;		// derivative of phi_e

	// Adaption parameters
	double			alpha_Na,	// Sodium influx per spike
					tau_Na,		// Sodium time constant
					g_KNa;		// KNa conductance

	// Firing rate parameters
	double			theta_e,	// pyramidal firing threshold
					sigma_e;	// pyramidal gain

	// Noise parameters
	double 			dphi_c;

	// pointer to the thalamic module
	Thalamic_Column*	Thalamus;

	// Random number generators
	vector<GEN>		MTRands;

	// Container for noise
	vector<double>	Rand_vars;
};
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/

