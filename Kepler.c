#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>
#include <libgen.h>
#include <math.h>

#define MG -1.32749e11 
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))

//----------------------------------------
// PROTOTYPES
//----------------------------------------
double   *Runge_Kutta(double, double, double, double, double, double *);
double   Leap_Frog_position(double, double);
double   Leap_Frog_xvelocity(double, double, double);
double   Leap_Frog_yvelocity(double, double, double);


//----------------------------------------
// *************MAIN**********************
//----------------------------------------
void main()
{
//printf("%f",Function_dvx(1.02554, 6.548) );
  int     i;
  double  x, y, aux_x, aux_y;
  double  vx, vy, aux_vx, aux_vy, vx_i, vy_i;
  double  delta_t;
  int     i_max;
  FILE   *fp_tmp;
 
double *result;
  result = (double *) calloc(4, sizeof(double));
//RUNGE-KUTTA
  //Boundary conditions
  x       = 0.0;
  y       = -1.0;
  vx      = 29.8e3/(1.496e11*sqrt(6.67e-11*1.99e+30/pow3(1.496e11)));    
  vy      = 0.0;
//  printf("Give the step in time: delta_t \n");
//  scanf("%f", &(delta_t));
  delta_t = 0.656;
  printf("Give the number of steps in time: i\n");
  scanf("%i", &(i_max));

  fp_tmp = fopen("rungekutta.dat", "w");

  for (i=0; i<i_max; i++)
    {
    result = Runge_Kutta( x, y, vx, vy, delta_t, result);


    x  = x + result[0];
    y  = y + result[1];
    vx = vx + result[2];
    vy = vy + result[3];

    fprintf(fp_tmp, "%f		%f\n", x,y);
    }
  fclose(fp_tmp);

//LEAP-FROG
  //Boundary conditions
  x       = 0.0;
  y       = -1.0;
  vx      = 29.8e3/(1.496e11*sqrt(6.67e-11*1.99e+30/pow3(1.496e11)));    //(km/s)
  vy      = 0.0;

  fp_tmp = fopen("leapfrog.dat", "w");

  vx = vx + delta_t/2.0*(-x/pow3(sqrt(pow2(x) + pow2(y)))); //vx_1/2
  vy = vy + delta_t/2.0*(-y/pow3(sqrt(pow2(x) + pow2(y)))); //vy_1/2

  for (i=1; i<(i_max-1); i++)
    {
    x  = x + Leap_Frog_position(vx, delta_t);
    y  = y + Leap_Frog_position(vy, delta_t);
    vx = vx + Leap_Frog_xvelocity(x, y, delta_t);
    vy = vy + Leap_Frog_yvelocity(x, y, delta_t);

    fprintf(fp_tmp, "%f		%f\n", x,y);
    }
  vx = vx + delta_t/2.0*(-x/pow3(sqrt(pow2(x) + pow2(y)))); //vx_N
  vy = vy + delta_t/2.0*(-y/pow3(sqrt(pow2(x) + pow2(y)))); //vy_N



  fclose(fp_tmp);

}



//**************************************************************************************
//----------------------------------------
//---------------Runge_Kutta()------------
//----------------------------------------
double   *Runge_Kutta(double x, double y, double vx, double vy, double delta_t, double *result)
{
  double kvx1, kvx2, kvx3, kvx4, result_vx;
  double kvy1, kvy2, kvy3, kvy4, result_vy;
  double kx1, kx2, kx3, kx4, result_x;
  double ky1, ky2, ky3, ky4, result_y;
  kx1     = vx;
  ky1     = vy;
  kvx1     = -x/(pow3(sqrt(pow2(x) + pow2(y))));
  kvy1     = -y/(pow3(sqrt(pow2(x) + pow2(y))));

  kx2     = vx + ( -x/(pow3(sqrt(pow2(x) + pow2(y)))) )*(delta_t/2.0);
  ky2     = vy + ( -y/(pow3(sqrt(pow2(x) + pow2(y)))) )*(delta_t/2.0);
  kvx2     = -(x+kx1*delta_t/2.0)/(pow3(sqrt(pow2(x+kx1*delta_t/2.0) + pow2(y+ky1*delta_t/2.0))));
  kvy2     = -(y+ky1*delta_t/2.0)/(pow3(sqrt(pow2(x+kx1*delta_t/2.0) + pow2(y+ky1*delta_t/2.0))));

  kx3     = vx + ( -(x+kx1*delta_t/2.0)/(pow3(sqrt(pow2(x+kx1*delta_t/2.0) + pow2(y+ky1*delta_t/2.0)))) )*(delta_t/2.0);
  ky3     = vy + ( -(y+ky1*delta_t/2.0)/(pow3(sqrt(pow2(x+kx1*delta_t/2.0) + pow2(y+ky1*delta_t/2.0)))) )*(delta_t/2.0);
  kvx3     = -(x+kx2*delta_t/2.0)/(pow3(sqrt(pow2(x+kx2*delta_t/2.0) + pow2(y+ky2*delta_t/2.0))));
  kvy3     = -(y+ky2*delta_t/2.0)/(pow3(sqrt(pow2(x+kx2*delta_t/2.0) + pow2(y+ky2*delta_t/2.0))));

  kx4     = vx + ( -(x+kx2*delta_t/2.0)/(pow3(sqrt(pow2(x+kx2*delta_t/2.0) + pow2(y+ky2*delta_t/2.0)))) )*delta_t;
  ky4     = vy + ( -(y+ky2*delta_t/2.0)/(pow3(sqrt(pow2(x+kx2*delta_t/2.0) + pow2(y+ky2*delta_t/2.0)))) )*delta_t;
  kvx4     = -(x+kx3*delta_t)/(pow3(sqrt(pow2(x+kx3*delta_t) + pow2(y+ky3*delta_t))));
  kvy4     = -(y+ky3*delta_t)/(pow3(sqrt(pow2(x+kx3*delta_t) + pow2(y+ky3*delta_t))));

  result_x = delta_t/6.0*(kx1 + 2.0*kx2 + 2.0*kx3 +kx4);
  result_y = delta_t/6.0*(ky1 + 2.0*ky2 + 2.0*ky3 +ky4);
  result_vx = delta_t/6.0*(kvx1 + 2.0*kvx2 + 2.0*kvx3 +kvx4);
  result_vy = delta_t/6.0*(kvy1 + 2.0*kvy2 + 2.0*kvy3 +kvy4);
  result[0] = result_x;
  result[1] = result_y;
  result[2] = result_vx;
  result[3] = result_vy;
  return result;
}


//----------------------------------------
// Leap_Frog_position()
//----------------------------------------
double   Leap_Frog_position(double v, double delta_t)
{
  return delta_t*v;
}


//----------------------------------------
// Leap_Frog_velocity()
//----------------------------------------
double   Leap_Frog_xvelocity(double x, double y, double delta_t)
{
  return delta_t*(-x/pow3(sqrt(pow2(x) + pow2(y))));
}

double   Leap_Frog_yvelocity(double x, double y, double delta_t)
{
  return delta_t*(-y/pow3(sqrt(pow2(x) + pow2(y))));
}







