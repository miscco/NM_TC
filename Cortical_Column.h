/****************************************************************************************************/
/*									header file of a cortical module								*/
/****************************************************************************************************/
#pragma once
#include <vector>
#include <cmath>
#include "macros.h"
#include "parameters.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include "Thalamic_Column.h"
using std::vector;
class Thalamic_Column;

/****************************************************************************************************/
/*										 typedefs for the RNG										*/
/****************************************************************************************************/
typedef boost::mt11213b                    	ENG;    // Mersenne Twister
typedef boost::normal_distribution<double>	DIST;   // Normal Distribution
typedef boost::variate_generator<ENG,DIST> 	GEN;    // Variate generator
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/


/****************************************************************************************************/
/*							implementation of the cortical module									*/
/****************************************************************************************************/
class Cortical_Column {
public:
	// Constructors
	Cortical_Column(void)
	: Ve	(_INIT(E_L_e)),	Vi 	   	(_INIT(E_L_i)),	Na	 	(_INIT(Na_eq)),
	  Phi_ee(_INIT(0.0)), 	Phi_ei 	(_INIT(0.0)), 	Phi_ie 	(_INIT(0.0)), 	Phi_ii	(_INIT(0.0)),	phi_e	(_INIT(0.0)),
	  x_ee 	(_INIT(0.0)), 	x_ei   	(_INIT(0.0)),	x_ie   	(_INIT(0.0)), 	x_ii	(_INIT(0.0)),	y_e  	(_INIT(0.0))
	{MTRands   = {{ENG(rand()), DIST (mphi_sc, dphi_sc)}, {ENG(rand()), DIST (mphi_sc, dphi_sc)}, {ENG(rand()), DIST (mphi_sc, dphi_sc)}, {ENG(rand()), DIST (mphi_sc, dphi_sc)}};
	 Rand_vars = {MTRands[0](), MTRands[1](), MTRands[2](), MTRands[3]()};}


	Cortical_Column(int N_Cols)
	: Ve	(_INIT(E_L_e)),	Vi 	   	(_INIT(E_L_i)),	Na	 	(_INIT(Na_eq)),
	  Phi_ee(_INIT(0.0)), 	Phi_ei 	(_INIT(0.0)), 	Phi_ie 	(_INIT(0.0)), 	Phi_ii	(_INIT(0.0)),  phi_e	(_INIT(0.0)),
	  x_ee 	(_INIT(0.0)), 	x_ei   	(_INIT(0.0)),	x_ie   	(_INIT(0.0)), 	x_ii	(_INIT(0.0)),  y_e  	(_INIT(0.0))
	{MTRands   = {{ENG(rand()), DIST (mphi_sc, dphi_sc)}, {ENG(rand()), DIST (mphi_sc, dphi_sc)}, {ENG(rand()), DIST (mphi_sc, dphi_sc)}, {ENG(rand()), DIST (mphi_sc, dphi_sc)}};
	 Rand_vars = {MTRands[0](), MTRands[1](), MTRands[2](), MTRands[3]()};}

	// get the pointer to the thalamic column
	void	get_Thalamus(Thalamic_Column& T) {Thalamus = &T;}

	// firing rate functions
	double 	get_Qe		(int) const;
	double 	get_Qi		(int) const;
	double 	get_phi		(int) const;

	// current functions
	double 	I_ee		(int) const;
	double 	I_ei		(int) const;
	double 	I_ie		(int) const;
	double 	I_ii		(int) const;
	double 	I_L_e		(int) const;
	double 	I_L_i		(int) const;
	double 	I_KNa		(int) const;

	// potassium concentration
	double 	Na_pump		(int) const;

	// external input functions
	double 	get_inp_e	(int) const;
	double 	get_inp_i	(int) const;

	// noise functions
	double 	noise_xRK 	(int, int) const;

	// ODE functions
	void 	set_RK		(int);
	void 	add_RK	 	(void);

	friend void get_data (int, Cortical_Column&, _REPEAT(double*, 1));

private:
	// population variables
	vector<double> 	Ve,			// excitatory membrane voltage
					Vi,			// inhibitory membrane voltage
					Na,			// Na concentration
					Phi_ee,		// PostSP from excitatory to excitatory population
					Phi_ei,		// PostSP from excitatory to inhibitory population
					Phi_ie,		// PostSP from inhibitory to excitatory population
					Phi_ii,		// PostSP from inhibitory to inhibitory population
					phi_e,		// axonal flux from pyramidal population
					x_ee,		// derivative of Phi_ee
					x_ei,		// derivative of Phi_ei
					x_ie,		// derivative of Phi_ie
					x_ii,		// derivative of Phi_ii
					y_e;		// derivative of phi_e

	// pointer to other thalamic columns
	Thalamic_Column*Thalamus;

	// random number generators
	vector<GEN>		MTRands;

	// container for random numbers
	vector<double>	Rand_vars;
};
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/
