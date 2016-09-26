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
/*						Implementation of a thalamic module					  */
/******************************************************************************/
#pragma once
#include <cmath>
#include <vector>

#include "Cortical_Column.h"
#include "Random_Stream.h"
class Cortical_Column;

class Thalamic_Column {
public:
    /* Constructor for simulation */
    Thalamic_Column(double* Param, double* Con)
        : g_LK		(Param[0]),	g_h 	(Param[1]),
          N_tp 		(Con[0]),	N_rp	(Con[1])
    {set_RNG();}

    /* Get the pointer to the cortical module */
    void	get_Cortex	(Cortical_Column& C) {Cortex = &C;}

    /* ODE functions */
    void 	set_RK		(int);
    void 	add_RK	 	(void);

    /* Set strength of external input */
    void	set_input	(double I) {input = I;}

private:
    /* Declaration of private functions */
    /* Initialize the RNGs */
    void 	set_RNG		(void);

    /* Firing rates */
    double 	get_Qt		(int) const;
    double 	get_Qr		(int) const;

    /* Synaptic currents */
    double 	I_et		(int) const;
    double 	I_gt		(int) const;
    double 	I_er		(int) const;
    double 	I_gr		(int) const;

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

    /* Helper functions */
    inline std::vector<double> init (double value)
    {return {value, 0.0, 0.0, 0.0, 0.0};}

    inline void add_RK (std::vector<double>& var)
    {var[0] = (-3*var[0] + 2*var[1] + 4*var[2] + 2*var[3] + var[4])/6;}

    inline void add_RK_noise (std::vector<double>& var, unsigned noise)
    {var[0] = (-3*var[0] + 2*var[1] + 4*var[2] + 2*var[3] + var[4])/6 + noise_aRK(noise);}

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
    const double 	sigma_t		= 6.;
    const double 	sigma_r		= 6.;

    /* Scaling parameter for sigmoidal mapping (dimensionless) */
    const double 	C1          = (M_PI/sqrt(3));

    /* PSP rise time in ms^-1 */
    const double 	gamma_e		= 70E-3;
    const double 	gamma_g		= 100E-3;

    /* Axonal flux time constant in ms^-1*/
    const double 	nu			= 120E-3;

    /* Membrane capacitance in muF/cm^2 */
    const double	C_m			= 1.;

    /* Leak weight in aU */
    const double 	g_L    		= 1.;

    /* Synaptic weights in ms */
    const double 	g_AMPA 		= 1.;
    const double 	g_GABA 		= 1.;

    /* Conductivities */
    /* Potassium leak current in mS/m^2 */
    const double 	g_LK 		= 0.02;

    /* T current in mS/m^2 */
    const double	g_T_t		= 3;
    const double	g_T_r		= 2.3;

    /* h current in mS/m^2 */
    const double	g_h			= 0.06;

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
    const double	dphi		= 20E-3;
    double			input		= 0.0;

    /* Connectivities (dimensionless) */
    const double 	N_rt		= 3.;
    const double 	N_tr		= 5.;
    const double 	N_rr		= 25.;

    /* Connectivities from cortex (dimensionless) */
    const double 	N_tp		= 2.6;
    const double 	N_rp		= 2.6;

    /* Pointer to cortical column */
    Cortical_Column* Cortex;

    /* Parameters for SRK4 iteration */
    const std::vector<double> A = {0.5,  0.5,  1.0, 1.0};
    const std::vector<double> B = {0.75, 0.75, 0.0, 0.0};

    /* Random number generators */
    std::vector<randomStreamNormal> MTRands;

    /* Container for noise */
    std::vector<double>	Rand_vars;

    /* Population variables																			*/
    std::vector<double> Vt		= init(E_L_t),		/* TC membrane voltage								*/
                        Vr		= init(E_L_r),		/* RE membrane voltage								*/
                        Ca		= init(Ca_0),		/* Calcium concentration of TC population			*/
                        s_et	= init(0.0),		/* PostSP from TC population to TC population		*/
                        s_er	= init(0.0),		/* PostSP from TC population to RE population		*/
                        s_gt	= init(0.0),		/* PostSP from RE population to TC population		*/
                        s_gr	= init(0.0),		/* PostSP from RE population to RE population		*/
                        y		= init(0.0),		/* axonal flux										*/
                        x_et	= init(0.0),		/* derivative of s_et								*/
                        x_er	= init(0.0),		/* derivative of s_er								*/
                        x_gt	= init(0.0),		/* derivative of s_gt								*/
                        x_gr	= init(0.0),		/* derivative of s_gr								*/
                        x		= init(0.0),		/* derivative of y									*/
                        h_T_t	= init(0.0),		/* inactivation of T channel						*/
                        h_T_r	= init(0.0),		/* inactivation of T channel						*/
                        m_h		= init(0.0),		/* activation 	of h   channel						*/
                        m_h2	= init(0.0);		/* activation 	of h   channel bound with protein 	*/

    /* Data storage  access */
    friend void get_data (int, Cortical_Column&, Thalamic_Column&, std::vector<double*>);
    friend class Cortical_Column;
};
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
