#ifndef __INTEGRATION_H__
#define __INTEGRATION_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>

#include "mex.h"
#include "macros.h"

/**
*	@brief Function to define the ODEs for ship roll model

* 	\param 	t 			Current time of integration
*	\param 	posIn[2] 	Vector of initial condition
*	\param  VelOut[2]   Right hand side of the system of ODEs
*  	\param  params      Parameters for the odes

*  	\return Returns nothing
*/
int ship_roll_model(double t, const double posIn[2], double velOut[2], void *params);

/** 
* 	@brief  			Function to evolve an initial condition for 2D vector fields

* 	\param  curr_t      Current time of integration, not my current variable naming style. 
* 	\param  tau         Length of integration i.e., final time = curr_t + tau
* 	\param  curr_pos[]  Initial condition for integration 
* 	\param  *iter_pos   Final result of integration
*	\param 	params 		Obsolete input argument, as the parameters are now being set as 
						global variables

* 	\return Returns nothing
*/

int ship_roll_model_stoc_forcing(double t, const double posIn[2], double velOut[2], void *params);

int evolve_pt(double curr_t, double tau, double curr_pos[], double *iter_pos, double params);

/**
* 	@brief  			Function to evolve an initial condition with check for some event in time

* 	\param  curr_t      Current time of integration, not my current variable naming style. 
*	\param  tau         Length of integration i.e., final time = curr_t + tau
*	\param  curr_pos[]  Initial condition for integration 
*	\param  *iter_pos   Final result of integration

*	\return Returns nothing
*/

int evolve_pt_event(double curr_t, double tau, double curr_pos[], double *iter_pos, double params);

double event_fun(const double Y[], double valueFun);

/**
*	@brief 				Function to rotate trajectories for twisted pipe; duplicate function exists under perron.c

**/

double eventfun(const double Y[],int which); 

#endif
