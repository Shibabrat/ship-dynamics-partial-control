#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "macros.h"

extern	char	SYSTEM[SHORTSTRING];		///< Global variable for the name of system
extern 	char 	GEOMETRY[SHORTSTRING];  	///< Global variable for the geometry of domain 
extern 	char 	PATH_OUTPUT[LONGSTRING];
extern 	char 	PMAP_OUTPUT_PATH[LONGSTRING];	
extern 	char 	DATA_OUTFILE[SHORTSTRING];
extern 	char 	BASE_FILENAME[SHORTSTRING];

extern 	double 	PARAMETER_01;
extern	double	PARAMETER_02;
extern 	double	PARAMETER_03;

extern 	double 	T_BEGIN;
extern 	double	T_END;
extern  double  T_DELTA;
extern 	double 	T_STAR;
extern 	double	BOXES;
extern 	double	TESTPOINTS;

extern double DIM_01_MIN;
extern double DIM_01_MAX;
extern double DIM_02_MIN;
extern double DIM_02_MAX;
extern long int DIM_01_BOXES;
extern long int DIM_02_BOXES;

extern 	int 	NDIM;
extern 	int 	DIRECTION; 

//Function prototypes
/**
*	@brief	Function to set parameters are at compile time

* 	@return Returns nothing
**/
void set_system_parameters();

/**
*	@brief	Function to receive command line input, this version (20th Dec 2013) accepts 
			time step for transfer operator computation or the time instant to observe 
			Poincare sections.

* 	@return Returns nothing
**/
void get_system_parameters(char *argv[]);

/**
*	@brief	Function to set paramters that are derived from input parameters or paths

*	@return Returns nothing
**/
void set_derived_parameters();

/**
*	@brief	Function to copy data to output directory using indexed file names

*	@return Returns nothing
**/
void move_data_to_output_dir(long int datafirstframe, long int datalastframe);

/**
*	@brief 	Function to move file to outpur directory
	
* 	@return Returns nothing
**/
void move_file_to_output_dir(char fileName[]);

#endif