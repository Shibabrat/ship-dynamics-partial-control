/**

*	@file 	parameters.c

*	@brief 	Functions for setting parameters for transfer operator and Poincare section computations.

*   @bug	All the variables are declared in parameters.h and defined in globals.h    

*   @copyright  MERL[2013]
**/
#include "parameters.h"

void set_system_parameters()

{
	/*Double gyre flow*/
//  	 strcpy(SYSTEM,"double_gyre");
//  	 strcpy(GEOMETRY, "rectangle");
 	 //strcpy(PATH_OUTPUT, "/homes/naik/data/double-gyre/transfer-op-output");		//Local directory for files
 	 //strcpy(PMAP_OUTPUT_PATH, "/homes/naik/data/double-gyre/poincare-map-output");
//  	 strcpy(PATH_OUTPUT, "/home/uwhpsc/share-windows/mechatronics-MERL/transfer-op-output");		//Local directory for files
//  	 strcpy(PMAP_OUTPUT_PATH, "/home/uwhpsc/share-windows/mechatronics-MERL/poincare-map-output");		//Local directory for files
 	 // strcpy(PATH_OUTPUT, "/data/naik/no-blinking-no-slip/boxes-4096/case01");		//Server directory for files
//  	 strcpy(BASE_FILENAME, "double_gyre");
     
     /* Ship roll model */
//     strcpy(SYSTEM,"ship_roll_model");
//     strcpy(GEOMETRY,"rectangle");
//     strcpy(PATH_OUTPUT, "/Users/shibabrat/Documents/Work/Research-Computations/partial-control-sculpting");		//Local directory for files
//     strcpy(PMAP_OUTPUT_PATH, "/Users/shibabrat/Documents/Work/Research-Computations/partial-control-sculpting");		//Local directory for files
//  	strcpy(PATH_OUTPUT, "/data/naik/no-blinking-no-slip/boxes-4096/case01");		//Server directory for files
//  	strcpy(BASE_FILENAME, "ship_roll_model");
     
 	 // /*Parameters for the vector field, see example function in vector_field.c*/
//  	 PARAMETER_01 = 4.94;
//  	 PARAMETER_02 = 90*(pi/180);
//  	 PARAMETER_03 = 2*pi;
 
 	 // /*Time variables*/
//  	 T_BEGIN = 0;						//Time to begin the computation of transfer operator
//  	 T_END = 2.0;							//Time to end the computation of transfer operator
//  	 TESTPOINTS = 100;
 
 	 // /*Problem domain*/
//  	 DIM_01_MIN = 0;
// //  	 DIM_01_MAX = 2;
//  	 DIM_02_MIN = 0;
//  	 DIM_02_MAX = 1;
//  	 DIM_01_BOXES = 10;
//  	 DIM_02_BOXES = 10;
 
 	 // /*Dimension of the vector field and the direction of integration*/
 	 NDIM = 2;
//  	 DIRECTION = 1;					//Backward time integration: -1

	/*Blinking/no blinking flow with no slip*/
//	strcpy(SYSTEM,"no_blinking_no_slip");
//	strcpy(GEOMETRY, "rectangle");
//	strcpy(PATH_OUTPUT, "/homes/naik/data/no-blinking-no-slip/transfer-op-output");		//Local directory for files
//	// strcpy(PATH_OUTPUT, "/data/naik/no-blinking-no-slip/boxes-4096/case01");		//Server directory for files 
//	strcpy(PMAP_OUTPUT_PATH, "/homes/naik/data/no-blinking-no-slip/poincare-map-output");
//	strcpy(BASE_FILENAME, "no_blinking_no_slip");
//
//	//Parameters for the vector field, see example function in vector_field.c
//	PARAMETER_01 = 0.0;
//	PARAMETER_02 = 1.2;
//	PARAMETER_03 = 0.0;
//
//	//Time variables
//	T_BEGIN = 0;						//Time to begin the computation of transfer operator
//	T_END = 0.1;							//Time to end the computation of transfer operator
//	//BOXES = 400;
//	TESTPOINTS = 100;
//
//	//Problem domain
//	DIM_01_MIN = -1;
//	DIM_01_MAX = 2;
//	DIM_02_MIN = -1;
//	DIM_02_MAX = 1;
//	DIM_01_BOXES = 40;
//	DIM_02_BOXES = 40;
//
//	//Dimension of the vector field and the direction of integration
// 	NDIM = 2;
//	DIRECTION = 1;					//Backward time integration: -1

	/*Twisted pipe 2D model:Jones,Aref-1989*/
	// strcpy(SYSTEM,"twisted_pipe_theta");
	// strcpy(GEOMETRY, "circle");
	// strcpy(PMAP_OUTPUT_PATH,"/homes/naik/data/twisted-pipe/poincare-map-output");
	// strcpy(PATH_OUTPUT, "/homes/naik/data/twisted-pipe/transfer-op-output");	
	// strcpy(BASE_FILENAME, "twisted_pipe_theta");

	// //Parameters for the vector field, see example function in vector_field.c
	// PARAMETER_01 = 400;
	// PARAMETER_02 = 4;
	// PARAMETER_03 = pi/16;

	// //Time variables
	// T_BEGIN = 0;						//Start angle the computation of transfer operator
	// T_END = 250*2*pi;					//Angle of pipe to end the computation of transfer operator
	// TESTPOINTS = 100;

	// //Problem domain
	// DIM_01_MIN = 0;
	// DIM_01_MAX = 1;
	// DIM_02_MIN = 0;
	// DIM_02_MAX = 2*pi;
	// DIM_01_BOXES = 1;
	// DIM_02_BOXES = 1;

	// //Dimension of the vector field and the direction of integration
	// NDIM = 2;
	// DIRECTION = 1;					//Backward time integration: -1

}

// void get_system_parameters(char *argv[])
// 
// {
// 	T_DELTA = atof(argv[1]);		//Time step for transfer operator computation, aka operator splitting time step
// 										//Plane to observe Poincare section 
// }

// void set_derived_parameters()
// 
// {
// 	// Check for correct input path name otherwise do the necessary
// 	if( PATH_OUTPUT[strlen(PATH_OUTPUT)-1] != '/' )
// 		sprintf(PATH_OUTPUT,"%s/",PATH_OUTPUT);
// 
// }

// void move_data_to_output_dir(long int datafirstframe, long int datalastframe)
// 
// {
// 
// 	char syscommand[LONGSTRING];
//   	long int f0, f1, ff;
// 
//   	// set_derived_parameters();
// 
// 	printf("Moving data to output directory ... \n");
// 	if(datafirstframe < datalastframe)
// 	{
// 		f0 = datafirstframe;
// 		f1 = datalastframe;
// 	}
// 	else
// 	{
// 		f0 = datalastframe;
// 		f1 = datafirstframe;
// 	}
// 
// 	// printf("%s\n", DATA_OUTPUT);
// 	for(ff=f0;ff<=f1;++ff)
// 	{
// 		sprintf(syscommand,"mv -u %s%03ld.txt %s", DATA_OUTFILE,ff,PATH_OUTPUT);
// 		system(syscommand);
// 	}

// }
// 
// void move_file_to_output_dir(char fileName[])
// {
// 
// 	char syscommand[LONGSTRING];
// 
// 	sprintf(syscommand,"mv -u %s %s",fileName,PATH_OUTPUT);
// 	system(syscommand);
// 
// }
