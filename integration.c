/** 
*   @file       integration.c

*   @brief      Functions to call GSL library options for integrating ODEs.

*   @bug        Event check integration may need more tests

*   @copyright  Shibabrat Naik[2015]
**/
#include "integration.h"

int ship_roll_model(double t, const double posIn[2], double velOut[2], void *params){

	int dirtn = 1;
    
    double Mt, phi, pPhi;
	/* Roll angle of vanishing stability */
    double phiCritical = 0.88;
    double H = 4.94;
    double chi = 90*(pi/180);
    
    /* Parameters for Edith Terkol */
    double b1 = 0.0043, b2 = 0.0225, c1 = 0.384, c2 = 0.1296, 
        c3 = 1.0368, c4 = -4.059, c5 = 2.4052, I = 1174, wM = 0;
    
    /* Parameters for Edith Terkol, without damping */
//     double b1 = 0.0, b2 = 0.0, c1 = 0.384, c2 = 0.1296, 
//         c3 = 1.0368, c4 = -4.059, c5 = 2.4052, I = 1174, wM = 0;
    
    /* Disturbance due to regular seas */
    double alpha0 = 0.73, omegaN = 0.62, omegaE = 0.527, lambda = 221.94;
    
    Mt = I*alpha0*pow(omegaN,2)*pi*(H/lambda)*(sin(chi))*sin(omegaE*t);
	// mexPrintf("%lf\n",Mt);

    phi = posIn[0]; 
    pPhi = posIn[1];
    
    velOut[0] = dirtn*pPhi;
    velOut[1] = dirtn*(-b1*pPhi - b2*fabs(pPhi)*pPhi - c1*phi - 
            c2*fabs(phi)*phi - c3*pow(phi,3) -c4*fabs(phi)*pow(phi,3) - 
            c5*pow(phi,5) + Mt/I);
    
    return GSL_SUCCESS;
}

int ship_roll_model_stoc_forcing(double t, const double posIn[2], double velOut[2], void *params){

	int dirtn = -1, i, N;

	double omegaI, omegaIE, tempRatio, waveEnergySpect, fPhi = 0, 
		phi, pPhi;
	/* Roll angle of vanishing stability */
    double phiCritical = 0.88;
    double H = 4.94;
    double chi = 90*(pi/180);

	/* Parameters for Edith Terkol */
    double b1 = 0.0043, b2 = 0.0225, c1 = 0.384, c2 = 0.1296, 
        c3 = 1.0368, c4 = -4.059, c5 = 2.4052, I = 1174, wM = 0;

    /* Disturbance due to regular seas */
    double alpha0 = 0.73, dOmega = 0.01, finOmega = 2, initOmega = 0.01, 
    	omegaNPhi = 0.62, omegaZ = 0.527;
	double U = 4*0.514444444; // m/s vessel speed, 1 knot = 0.514444444 m / s

    double epsilonI;

    N = ( finOmega - initOmega )/dOmega + 1;
    // printf("%d\n", N);

	// const gsl_rng_type *T;
	// gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
	// gsl_rng_env_setup();
	// T = gsl_rng_default;
	// r = gsl_rng_alloc (T);


    gsl_rng *r;
    if((r = gsl_rng_alloc(gsl_rng_mt19937)) == NULL) {
    printf("ERROR: Could not create random number generator\n");
    exit(1);
  	}
  	gsl_rng_set(r, 1167);

	char fileName[] = "test_forcing.txt";
	FILE *writeFile;
	writeFile = fopen(fileName,"w");
    for (i = 0; i < N; ++i)
    {
		epsilonI = 2*pi*gsl_rng_uniform (r);
		// mexPrintf ("%lf\n", epsilonI);
		fprintf(writeFile, "%lf\n", epsilonI);

		omegaI = initOmega + i*dOmega;
		omegaIE = omegaI - ((pow(omegaI,2)*U)/gAdg)*cos(chi);
    	tempRatio = pow((omegaZ/omegaI),4);
    	waveEnergySpect = 0.11*pow(H,2)*(tempRatio/omegaI)*exp(-0.44*tempRatio);
		fPhi = fPhi + pow(omegaI,2)*sqrt(waveEnergySpect)*sin(omegaIE*t + epsilonI);	
	}	
	fPhi = pow(omegaNPhi,2)*sin(chi)*alpha0*(sqrt(2*dOmega)/gAdg)*fPhi;
	
	
	
	fclose(writeFile);

	// mexPrintf("%lf\n",fPhi);
	gsl_rng_free (r);

	phi = posIn[0]; 
    pPhi = posIn[1];
    
    velOut[0] = dirtn*pPhi;
    velOut[1] = dirtn*(-b1*pPhi - b2*fabs(pPhi)*pPhi - c1*phi - 
            c2*fabs(phi)*phi - c3*pow(phi,3) -c4*fabs(phi)*pow(phi,3) - 
            c5*pow(phi,5) + fPhi);

    return GSL_SUCCESS;

}

int evolve_pt(double curr_t, double tau, double curr_pos[], double *iter_pos, double params)

{
	// int NDIM = 2;
    //Setting up the GSL integrator routines
    const gsl_odeiv_step_type *T = gsl_odeiv_step_rkf45;
    gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, NDIM);
    
    gsl_odeiv_control *c = gsl_odeiv_control_y_new(1e-12, 1e-12);
    gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(NDIM);

    gsl_odeiv_system sys = {ship_roll_model, NULL, NDIM, &params};
//     gsl_odeiv_system sys = {ship_roll_model_stoc_forcing, NULL, NDIM, &params};

    int i;
    double t, tNext;
    double tMin, tMax, deltaT;
    double h = 1e-10; 
    double y[NDIM];
    
    tMin = curr_t;
    tMax = curr_t + tau;
    deltaT = tau/1;

    //Initial conditions for integrator
    for (i = 0; i < NDIM; i++)
    {
        y[i] = curr_pos[i];        
    }
    
    t = tMin;    //Initialization of time variable
    for (tNext = tMin + deltaT; tNext <= tMax; tNext += deltaT)
    {
        while (t < tNext)
        {
            int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, tNext, &h, y);
        
            if (status != GSL_SUCCESS)
            {
                printf ("error, return value=%d\n", status);
                break;
            }
        }
        // printf("%lf %lf %lf\n", t, y[0],y[1]);
    }

    for (i = 0; i < NDIM; i++)
    {
        iter_pos[i] = y[i]; 
    }

    // printf("%lf %lf %lf %lf\n",curr_pos[0], curr_pos[1], iter_pos[0], iter_pos[1]);
    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);
    
    return EXIT_SUCCESS;
}

int evolve_pt_event(double curr_t, double tau, double curr_pos[], double *iter_pos, double params)

/* Fixed time stepping to locate event during integration
*
*/


{
	int i, j, k;			
	double deltaT = 1e-2;		
	int isterminal = 1;
	double eventTime, eventFlag = 1;
	double phiCritical = 0.88;
	double new_curr_pos[NDIM], new_iter_pos[NDIM];

	double ti = curr_t, tf = curr_t + tau;
	long int numSteps = 0;

	for (i = 0; i < NDIM; ++i){		
		new_curr_pos[i] = curr_pos[i];
		new_iter_pos[i] = iter_pos[i];
	}

	/* Evolving the point*/
	while (ti < tf){
		evolve_pt(ti, deltaT, new_curr_pos, new_iter_pos, params);
		numSteps = numSteps + 1;	
		// printf("Time instant %e\n",ti);

		/* Trying to catch the event, if it occured */
		if (event_fun(new_iter_pos,phiCritical) > 0 || event_fun(new_iter_pos,-phiCritical) < 0){
			if (isterminal == 1){
				eventTime = ti;
				eventFlag = 0;
				
				for (i = 0; i < NDIM; ++i)	
					iter_pos[i] = new_iter_pos[i];

				break;
			}
			else
				continue;
			// printf("Capsize, stop integration and locate the zero of the event\n");
		}
		else if (event_fun(new_iter_pos,phiCritical) < 0 || event_fun(new_iter_pos,-phiCritical) > 0){
			// printf("Safe, so go to the next time step\n");
			eventTime = ti;
			eventFlag = 1;
		}

		ti = ti + deltaT;
		for (i = 0; i < NDIM; ++i)
			new_curr_pos[i] = new_iter_pos[i];	
	}

	char fileName[] = "tests.txt";
	FILE *writeFile;
	writeFile = fopen(fileName,"a");
	if (ti >= tf)
	{
		for (k = 0; k < NDIM - 1; ++k)
			fprintf(writeFile, "%lf\t", curr_pos[k]);

		fprintf(writeFile, "%lf\n", curr_pos[NDIM-1]);
		for (i = 0; i < NDIM; ++i)	
			iter_pos[i] = new_iter_pos[i];
	}
	
	fclose(writeFile);

	return GSL_SUCCESS;	
}

double event_fun(const double Y[], double valueFun) 

/**
*   @brief      This is the event function (0-crossings to be detected), "which" provides the direction
                Use "which" to -1 for detecting positive->negative 0-crossing, +1 for detecting negative->positive 0-crossing
**/
{ 
    double theta = Y[0];
    double result;
    result = (theta - valueFun);

    return result; 
}

double eventfun(const double Y[],int which) 
//This is the event function (0-crossings to be detected), "which" provides the direction
//Use "which" to -1 for detecting positive->negative 0-crossing, +1 for detecting negative->positive 0-crossing
{ 
	int d;
	double x = Y[0];
	double y = Y[1];
	double theta = Y[2];
	double result;
	// result = which*(theta - pi);
	result = which*(x - 0.88);

	return result; 
}
 
