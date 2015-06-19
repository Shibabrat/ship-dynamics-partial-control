/**
*   @file   vector_field.c

*   @brief  Module for defining the ODEs for a dynamical system.

*   @bug    Parameters are not sent as input to the velocity field rather set as global variables in parameters.c

*   @copyright MERL[2013]
**/

#include "vector_field.h"

int velocity_field(double t, const double posIn[], double velOut[], void *params)
{

    //  This is where a new vector field can be introduced
	if ( strcmp(SYSTEM,"ship_roll_model") == 0 ){
		ship_roll_model(t, (double *)posIn, (double *)velOut, params);
	}
	else if ( strcmp(SYSTEM,"no_blinking_no_slip") == 0 )
    {
		no_blinking_no_slip(t, (double *)posIn, (double *)velOut, params);
    }
    else if ( strcmp(SYSTEM,"blinking_no_slip") == 0 )
    {
        blinking_no_slip(t, (double *)posIn, (double *)velOut, params);
    }
    else if ( strcmp(SYSTEM,"cell_slip") == 0 )
    {
		cell_slip(t, (double *)posIn, (double *)velOut, params);
    }
    else if ( strcmp(SYSTEM,"double_gyre") == 0 )
    {
		double_gyre(t, (double *)posIn, (double *)velOut, params);
    }
    else if ( strcmp(SYSTEM,"twisted_pipe_theta") == 0 )
    {
		twisted_pipe_theta(t, (double *)posIn, (double *)velOut, params);
    }
    else if ( strcmp(SYSTEM,"twisted_pipe_time") == 0)
    {
        twisted_pipe_time(t, (double *)posIn, (double *)velOut, params);
    }
    else if ( strcmp(SYSTEM,"time_periodic_sine_map") == 0 )
    {
		time_periodic_sine_map(t, (double *)posIn, (double *)velOut, params);
    }
    else
    {
        printf("New SYSTEM not defined, aborting!!!\n");
        exit(EXIT_FAILURE);
    }

	return GSL_SUCCESS;
}

	/* Vector field for ship roll model */
int ship_roll_model(double t, const double posIn[2], double velOut[2], void *params){

    int dirtn = DIRECTION;
    
    double Mt, phi, pPhi;
	/* Roll angle of vanishing stability */
    double phiCritical = 0.88;
    double H = PARAMETER_01;
    double chi = PARAMETER_02;
    
    /* Parameters for Edith Terkol */
    double b1 = 0.0043, b2 = 0.0225, c1 = 0.384, c2 = 0.1296, 
        c3 = 1.0368, c4 = -4.059, c5 = 2.4052, I = 1174, wM = 0;

    /* Disturbance due to regular seas */
    double alpha0 = 0.73, omegaN = 0.62, omegaE = 0.527, lambda = 221.94;
    
    Mt = I*alpha0*pow(omegaN,2)*pi*(H/lambda)*(sin(chi))*sin(omegaE*t);
    
    phi = posIn[0]; 
    pPhi = posIn[1];
    
    velOut[0] = pPhi;
    velOut[1] = -b1*pPhi - b2*fabs(pPhi)*pPhi - c1*phi - 
            c2*fabs(phi)*phi - c3*pow(phi,3) -c4*fabs(phi)*pow(phi,3) - 
            c5*pow(phi,5) + Mt/I;
    
    return GSL_SUCCESS;
}


int blinking_no_slip(double t, const double posIn[2], double velOut[2], void *params)
    
    /*Vector field for non-autonomous blinking flow with no slip*/
    
{

    int dirtn = DIRECTION;

    double tSwitch = PARAMETER_01, epsilon1 = PARAMETER_02, epsilon2 = PARAMETER_03;
    double epsilon;
    // double *params_ptr = (double *)params;
    double x,y, x1, x2;
    
    double c1,b1,c2,b2;

    x = posIn[0];
    y = posIn[1];

    if (fmod(t,2.0*tSwitch) < tSwitch)
    {
        epsilon = epsilon1;
        c1 = 0.75 - 0.5*epsilon;
        c2 = 0.75 + 0.5*epsilon;
        b1 = c1 - 1;
        b2 = 2 - c2;
        x1 = (x - b1)/c1;
        x2 = (x - b2)/c2;

        if (x < (0.5 - epsilon))
        {
            velOut[0] = -dirtn*0.5*pi*cos(0.5*pi*x1)*sin(0.5*pi*y)*sin(cos(0.5*pi*x1)*cos(0.5*pi*y));
            velOut[1] = dirtn*(0.5*pi/c1)*sin(0.5*pi*x1)*cos(0.5*pi*y)*sin(cos(0.5*pi*x1)*cos(0.5*pi*y));
        }
        else if (x >= (0.5 - epsilon))
        {
            velOut[0] = dirtn*0.5*pi*cos(0.5*pi*x2)*sin(0.5*pi*y)*sin(cos(0.5*pi*x2)*cos(0.5*pi*y));
            velOut[1] = -dirtn*(0.5*pi/c2)*sin(0.5*pi*x2)*cos(0.5*pi*y)*sin(cos(0.5*pi*x2)*cos(0.5*pi*y));
        }
    }    
    else if (fmod(t,2.0*tSwitch) > tSwitch)
    {
        epsilon = -epsilon2;
        c1 = 0.75 - 0.5*epsilon;
        c2 = 0.75 + 0.5*epsilon;
        b1 = c1 - 1;
        b2 = 2 - c2;
        x1 = (x - b1)/c1;
        x2 = (x - b2)/c2;

        if (x < (0.5 - epsilon))
        {
            velOut[0] = -dirtn*0.5*pi*cos(0.5*pi*x1)*sin(0.5*pi*y)*sin(cos(0.5*pi*x1)*cos(0.5*pi*y));
            velOut[1] = dirtn*(0.5*pi/c1)*sin(0.5*pi*x1)*cos(0.5*pi*y)*sin(cos(0.5*pi*x1)*cos(0.5*pi*y));
        }
        else if (x >= (0.5 - epsilon))
        {
            velOut[0] = dirtn*0.5*pi*cos(0.5*pi*x2)*sin(0.5*pi*y)*sin(cos(0.5*pi*x2)*cos(0.5*pi*y));
            velOut[1] = -dirtn*(0.5*pi/c2)*sin(0.5*pi*x2)*cos(0.5*pi*y)*sin(cos(0.5*pi*x2)*cos(0.5*pi*y));
        }
    }

    return GSL_SUCCESS;
}

int no_blinking_no_slip(double t, const double posIn[], double velOut[], void *params)
      
{

    int dirtn = DIRECTION;

    // double *params_ptr = (double *)params;
    double x,y, x1, x2;
    double epsilon;
    double c1,b1,c2,b2;
    
    epsilon = PARAMETER_02;
    // epsilon = *params_ptr;
    // printf("%lf\n", epsilon);

    x = posIn[0];
    y = posIn[1];
    velOut[0] = 0.0;
    velOut[1] = 0.0;   

    c1 = 0.75 - 0.5*epsilon;
    c2 = 0.75 + 0.5*epsilon;
    b1 = c1 - 1;
    b2 = 2 - c2;
    x1 = (x - b1)/c1;
    x2 = (x - b2)/c2;

    if (x < (0.5 - epsilon))
    {
        velOut[0] = -dirtn*0.5*pi*cos(0.5*pi*x1)*sin(0.5*pi*y)*sin(cos(0.5*pi*x1)*cos(0.5*pi*y));
        velOut[1] = dirtn*(0.5*pi/c1)*sin(0.5*pi*x1)*cos(0.5*pi*y)*sin(cos(0.5*pi*x1)*cos(0.5*pi*y));
    }
    else if (x >= (0.5 - epsilon))
    {
        velOut[0] = dirtn*0.5*pi*cos(0.5*pi*x2)*sin(0.5*pi*y)*sin(cos(0.5*pi*x2)*cos(0.5*pi*y));
        velOut[1] = -dirtn*(0.5*pi/c2)*sin(0.5*pi*x2)*cos(0.5*pi*y)*sin(cos(0.5*pi*x2)*cos(0.5*pi*y));
    }

    // printf("%lf %lf %lf %lf\n", posIn[0], posIn[1], velOut[0], velOut[1]);

    return GSL_SUCCESS;
}


int double_gyre(double t, const double posIn[], double velOut[], void *params)
{
    int dirtn = DIRECTION;
    double *params_ptr = (double *)params;
    double x, y, f, dfdx;
    double A = PARAMETER_01, epsilon = *params_ptr, omega = PARAMETER_03;
// 	printf("%lf, %lf, %lf\n", A, epsilon, omega);

    x = posIn[0];
    y = posIn[1];
    velOut[0] = 0.0;
    velOut[1] = 0.0;    

//Not having if-statement, improves time by 0.21 secs for 1600 integrations

    f = epsilon*sin(omega*t)*pow(x,2) + (1-2*epsilon*sin(omega*t))*x;
    dfdx = 2*epsilon*sin(omega*t)*x + (1-2*epsilon*sin(omega*t));
    velOut[0] = -dirtn*pi*A*sin(pi*f)*cos(pi*y);
    velOut[1] = dirtn*pi*A*cos(pi*f)*sin(pi*y)*dfdx;

    return GSL_SUCCESS;
}

int cell_slip(double t, const double posIn[], double velOut[], void *params)
{
    double x,y;
    int dirtn = DIRECTION;
    x = posIn[0];
    y = posIn[1];

    velOut[0] = dirtn*(-pi/2)*cos(pi*x/2)*sin(pi*y/2);
    velOut[1] = dirtn*(pi/2)*sin(pi*x/2)*cos(pi*y/2);

    return GSL_SUCCESS;
}

int twisted_pipe_theta(double t, const double posIn[],double velOut[], void *params)
{
    int dirtn = DIRECTION;
    double alpha = PARAMETER_01; 

    // double *params_ptr = (double *)params;
    double x,y;
    x = posIn[0];
    y = posIn[1];
    velOut[0] = 0.0;
    velOut[1] = 0.0;    

    // printf("%lf %lf %lf %lf\n", posIn[0], posIn[1], velOut[0], velOut[1]);

    // velOut[0] = dirtn*(alpha/288.0)*(4 - 5*pow(x,2) - 23*pow(y,2) + pow(x,4) + 8*pow(x,2)*pow(y,2) + 7*pow(y,4));
    // velOut[1] = dirtn*(alpha/48)*x*y*( 3 - pow(x,2) - pow(y,2) );

    double beta = PARAMETER_02;
    velOut[0] = dirtn*((alpha/beta)/1152)*( 4 - 5*pow(x,2) - 23*pow(y,2) + pow(x,4) + 8*pow(x,2)*pow(y,2) + 7*pow(y,4) );
    velOut[1] = dirtn*((alpha/beta)/192)*x*y*( 3 - pow(x,2) - pow(y,2) );
    
    // printf("%lf %lf %lf %lf\n", posIn[0], posIn[1], velOut[0], velOut[1]);

    return GSL_SUCCESS;
}

int twisted_pipe_time(double t, const double Y[], double velOut[], void *params)

//THIS PROVIDES THE R.H.S OF THE SYSTEM OF ODES

{

    int dirtn = DIRECTION;
    // double *params_ptr = (double *)params;
    double x, y;
    double r, h, hDash;
    double alpha = PARAMETER_01, beta = PARAMETER_02;
    x = Y[0];
    y = Y[1];

    r = sqrt(pow(x,2) + pow(y,2));
    h = 0.25*(4 - pow(r,2) )*( pow((1 - pow(r,2)),2) );
    hDash = ( r*(pow(r,2)-1)*(4-pow(r,2)) - r*0.5*( pow((1-pow(r,2)),2) ) );

    velOut[0] = dirtn*(alpha/1152.0)*(h + (pow(y,2)/r)*hDash );
    velOut[1] = -dirtn*(alpha*x*y*hDash)/(1152*r);
    velOut[2] = dirtn*0.25*beta*(1 - pow(r,2));

    return GSL_SUCCESS;
}

int time_periodic_sine_map(double t, const double posIn[], double velOut[], void *params)
{
    double timePeriod = *(double *)params;
    double x,y;
    double phase = mod_double(t, timePeriod);

    x = posIn[0];
    y = posIn[1];

    /*Accounting for the switching of flow*/
    if ( (phase < timePeriod/2) & (phase >= 0) )
    {
    /*Accounting for the periodic boundary condition*/
        if ( (x < 0) || (x > 1) )
            x = mod_double(x,1);
        if ( (y < 0) || (y > 1) )
            y = mod_double(y,1);

        velOut[0] = sin(2*pi*y);
        velOut[1] = 0;
    }
    else if ( (phase < timePeriod) & (phase >= timePeriod/2) )
    {
    /*Accounting for the periodic boundary condition*/
        if ( (x < 0) || (x > 1) )
            x = mod_double(x,1);
        if ( (y < 0) || (y > 1) )
            y = mod_double(y,1);

        velOut[0] = 0;
        velOut[1] = sin(2*pi*x);
    }

    return GSL_SUCCESS;
}

double mod_double(double dividend, double divisor)
{
    int quotient = floor(dividend/divisor);
    double remainder = dividend - quotient*divisor;

    return remainder;
}
