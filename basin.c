#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_odeiv2.h>
/*
struct coord{
    float x;
    float y;
};

// Define the differential equation function
int differential_equation(double t, const double y[], double dydt[], void* params)
{
    double alpha = 0.1;       // coefficient b
    double X[] = {1.0, 2.0, 3.0};   // example values for Xn
    double pn[] = {1.0, 1.0, 1.0};  // example values for pn: attractive +1 or repulsive -1 magnet 
    double h = 0.1;      // example value for h
    int n = 3;            // number of terms in the summation

    dydt[0] = y[1];  // dx/dt = v

    double sum_term = 0.0;
    for (int i = 0; i < n; i++)
    {
        double distance = X[i] - y[0];  // Xn - x(t)
        double distance_sq = distance * distance;
        double denominator = distance_sq + h * h;
        double term = pn[i] * distance / (denominator * denominator * sqrt(denominator));
        sum_term += term;
    }

    dydt[1] = -alpha * y[1] - y[0] + sum_term;

    return 0;
}

int main()
{   
    FILE* file =  fopen("result", "w");

    // Define the ODE solver type
    gsl_odeiv2_system sys = {differential_equation, NULL, 2, NULL};

    // Create the ODE solver workspace
    gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);

    // Set initial conditions
    double t = 0;
    double y[2] = { 0, 2 };

    // Solve the differential equation
    double t_end = 20;
    while (t < t_end)
    {
        int status = gsl_odeiv2_driver_apply(driver, &t, t_end, y);

        if (status != 0)
        {
            printf("Error: ODE solver failed with status %d\n", status);
            break;
        }
        fprintf(file, "t = %.2f, x = %.6f, v = %.6f\n", t, y[0], y[1]);

        // Process the solution at each time step
        // Access the solution values in the y array
        // ...

        // Update the time and continue the loop
    }

    // Free memory
    
    gsl_odeiv2_driver_free(driver);
    fclose(file); 

    return 0;
}

*/

