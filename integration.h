#ifndef __INTEGRATION_H__
#define __INTEGRATION_H__

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_roots.h>
#include <math.h>

#include "macros.h"
// #include "parameters.h"
// #include "vector_field.h"

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
int evolve_pt(double curr_t, double tau, double curr_pos[], double *iter_pos, double params);

/**
* 	@brief  			Function to evolve an initial condition with check for some event in time

* 	\param  curr_t      Current time of integration, not my current variable naming style. 
*	\param  tau         Length of integration i.e., final time = curr_t + tau
*	\param  curr_pos[]  Initial condition for integration 
*	\param  *iter_pos   Final result of integration

*	\return Returns nothing
*/
int evolve_pt_time_event(double curr_t, double tau, double curr_pos[], double *iter_pos);

/**
*	@brief 				Function to integrate initial conditions for the time range and terminate
						integration when the event occurs. 
**/
int iterate_event_check(double timeSpan[], double curr[], double *output, double *eventTime, int *eventFlag);

/**
*	@brief 				Function to define the event function, which is being checked to terminate
						integration
**/
double eventfun(const double Y[], int which);

/**
*	@brief 				Function to rotate trajectories for twisted pipe; duplicate function exists under perron.c

**/
int rotate_pos(double loc_out[], double *posInSecondCell, double angleOfRotation);

#endif
