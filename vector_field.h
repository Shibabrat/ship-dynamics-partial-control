#ifndef __VECTOR_FIELD_H__
#define __VECTOR_FIELD_H__

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_roots.h>
#include <math.h>

#include "parameters.h"

/**
*	@brief		Function to control the vector fields using the SYSTEM name defined in parameters.c
**/
int velocity_field(double t, const double posIn[], double velOut[], void *params);

/**
*	@brief		Function to define double gyre velocity field
**/
int ship_roll_model(double t, const double posIn[], double velOut[], void *params);

/**
*	@brief		Function to define ODES for ship roll model 
**/
int double_gyre(double t, const double posIn[], double velOut[], void *params);

/**
*	@brief		Function to define cell slip velocity field
**/
int cell_slip(double t, const double posIn[], double velOut[], void *params);

/**
*	@brief		Function to define velocity field for non-autonomous blinking flow with no slip 
**/
int blinking_no_slip(double t, const double posIn[], double velOut[], void *params);

/**
*	@brief		Function to define velocity field for autonomous blinking flow with no slip 
**/
int no_blinking_no_slip(double t, const double posIn[], double velOut[], void *params);

/**
*	@brief		Function to define velocity field for twisted pipe with theta as independent variable (2D odes)
**/
int twisted_pipe_theta(double t, const double posIn[], double velOut[], void *params);

/**
*	@brief		Function to define velocity field for twisted pipe with time as independent variable (3D odes)
**/
int twisted_pipe_time(double t, const double posIn[], double velOut[], void *params);

/**
*	@brief		Function to define time periodic sine map
**/
int time_periodic_sine_map(double t, const double posIn[], double velOut[], void *params);

/**
*	@brief		Function to obtain modulus of real numbers
**/
double mod_double(double dividend, double divisor);

#endif
