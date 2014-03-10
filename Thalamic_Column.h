/************************************************************************************************/
/*								Header file of a thalamic module								*/
/************************************************************************************************/
#pragma once
#include <cmath>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include "macros.h"
#include "parameters.h"
#include "Cortical_Column.h"
using std::vector;
class Cortical_Column;

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
/*								Implementation of the thalamic module 								*/
/****************************************************************************************************/
class Thalamic_Column {
public:
	// Constructors
	Thalamic_Column(void)
	: Vt 	(_INIT(E_L_t)), Vr 	  	(_INIT(E_L_r)),	Ca	  	(_INIT(Ca_0)),
	  Phi_tt(_INIT(0.0)), 	Phi_tr 	(_INIT(0.0)), 	Phi_rt 	(_INIT(0.0)), 	Phi_rr 	(_INIT(0.0)),  	phi_t	(_INIT(0.0)),
	  x_tt 	(_INIT(0.0)), 	x_tr   	(_INIT(0.0)), 	x_rt   	(_INIT(0.0)),  	x_rr	(_INIT(0.0)),  	y_t	(_INIT(0.0)),
	  h_T_t	(_INIT(0.0)),	h_T_r  	(_INIT(0.0)),	m_T_t  	(_INIT(0.0)),	m_T_r	(_INIT(0.0)),
	  m_h	(_INIT(0.0)),	m_h2	(_INIT(0.0)),	P_h		(_INIT(0.0)),
	  N_tr 	(210), 		   	N_rt   	(410), 			N_rr 	(800),			input 	(0.0)
	{set_RNG();}

	// Constructors
	Thalamic_Column(double* Par)
	: Vt 	(_INIT(0)), 	Vr 	  	(_INIT(0)), 	Ca	  	(_INIT(Ca_0)),
	  Phi_tt(_INIT(0.0)), 	Phi_tr 	(_INIT(0.0)), 	Phi_rt 	(_INIT(0.0)), 	Phi_rr 	(_INIT(0.0)),  	phi_t	(_INIT(0.0)),
	  x_tt 	(_INIT(0.0)), 	x_tr   	(_INIT(0.0)), 	x_rt   	(_INIT(0.0)),  	x_rr	(_INIT(0.0)),  	y_t	(_INIT(0.0)),
	  h_T_t	(_INIT(0.0)),	h_T_r  	(_INIT(0.0)),	m_T_t  	(_INIT(0.0)),	m_T_r	(_INIT(0.0)),
	  m_h	(_INIT(0.0)),	m_h2	(_INIT(0.0)),	P_h		(_INIT(0.0)),
	  N_tr 	(Par[0]), 	   	N_rt	(Par[1]), 		N_rr	(Par[2]),		input 	(0.0)
	{set_RNG();}

	// get the pointer to the cortical module
	void	get_Cortex	(Cortical_Column& C) {Cortex = &C;}

	// change the strength of input
	void	set_input	(double I) {input = I;}

	// Initialize the RNGs
	void 	set_RNG		(void);

	// Firing rates
	double 	get_Qt		(int) const;
	double 	get_Qr		(int) const;

	// Synaptic currents
	double 	I_et		(int) const;
	double 	I_it		(int) const;
	double 	I_er		(int) const;
	double 	I_ir		(int) const;

	// Activation functions
	double  m_inf_T_t	(int) const;
	double  m_inf_T_r	(int) const;
	double  tau_m_T_t	(int) const;
	double  tau_m_T_r	(int) const;
	double  m_inf_h		(int) const;
	double  tau_m_h		(int) const;

	// Deactivation functions
	double  h_inf_T_t	(int) const;
	double  h_inf_T_r	(int) const;
	double  tau_h_T_t	(int) const;
	double  tau_h_T_r	(int) const;

	// Currents
	double 	I_L_t		(int) const;
	double 	I_L_r		(int) const;
	double 	I_LK_t		(int) const;
	double 	I_LK_r		(int) const;
	double 	I_T_t		(int) const;
	double 	I_T_r		(int) const;
	double 	I_h			(int) const;

	// Noise functions
	double 	noise_xRK 	(int) const;

	// ODE functions
	void 	set_RK		(int);
	void 	add_RK	 	(void);

	// Data storage
	friend void get_data (int, Cortical_Column&, Thalamic_Column&, _REPEAT(double*, 2));

private:
	// Population variables
	vector<double> 	Vt,			// TC membrane voltage
					Vr,			// RE membrane voltage
					Ca,			// Calcium concentration of TC population
					Phi_tt,		// PostSP from TC population to TC population
					Phi_tr,		// PostSP from TC population to RE population
					Phi_rt,		// PostSP from RE population to TC population
					Phi_rr,		// PostSP from RE population to RE population
					phi_t,		// axonal flux
					x_tt,		// derivative of Phi_tt
					x_tr,		// derivative of Phi_tr
					x_rt,		// derivative of Phi_rt
					x_rr,		// derivative of Phi_rr
					y_t,		// derivative of phi_t
					h_T_t,		// inactivation of T channel
					h_T_r,		// inactivation of T channel
					m_T_t,		// activation 	of T channel
					m_T_r,		// activation 	of T channel
					m_h,		// activation 	of h   channel
					m_h2,		// activation 	of h   channel bound with protein
					P_h;		// fraction of protein bound with calcium

	// Connectivities
	double			N_tr,		// TC to RE loop
					N_rt,		// RE to TC loop
					N_rr;		// RE self  loop

	// Noise parameters
	double 			input;

	// pointer to the cortical module
	Cortical_Column* 	Cortex;

	// Random number generators
	vector<GEN>		MTRands;

	// Container for noise
	vector<double>	Rand_vars;
};
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
