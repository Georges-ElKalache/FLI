#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#define pi 3.14159265358979323846

int differential_equation(double t, const double z[], double dzdt[], void* params)
{ 
    //double 1 = 1.00;
    double A = -1 -2*1*z[1];
    double B = -2*1*z[0];
    double C = -1 -2*1*z[1];

    dzdt[0] = z[2];
    dzdt[1] = z[3];
    dzdt[2] = -z[0] -2*1*z[0]*z[1];
    dzdt[3] = -z[1] -1*(z[0]*z[0] - z[1]*z[1]);
    dzdt[4] = z[6];
    dzdt[5] = z[7];
    dzdt[6] = A*z[4] + B*z[5];
    dzdt[7] = B*z[4] + C*z[5];

    return GSL_SUCCESS;
}

float entree() {
    float x;
    printf("Donner un nombre: ");
    scanf("%f",&x);
return x;
    }

int main(){

    clock_t start_time, end_time;
    double execution_time;
    start_time = clock();
    FILE* file =  fopen("file1", "w");

    double x = entree();
    double y = entree();
    double z[8] = {x, y, 0, 0, 0.5, 0.5, 0.5, 0.5};

    gsl_odeiv2_system sys = {differential_equation, NULL, 8, NULL};
    gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);
    
    double t = 0;
    int t_end = 400;
    for (int s = 1; s <= t_end; s++){
        int t1 = 100;
        double ti = s * t1 / 100.0;
        int status = gsl_odeiv2_driver_apply(driver, &t, ti, z);
        if (status != GSL_SUCCESS)
        {
            printf("error, return value=%d\n", status);
            break;
        }

        fprintf(file, "%f %f %f %f %f %f %f %f %f\n", t, z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7]);
    }
    gsl_odeiv2_driver_free (driver);
    fclose(file); 

    end_time = clock();
    execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Execution time: %f seconds\n", execution_time);
}

