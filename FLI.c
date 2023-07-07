#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#define pi 3.14159265358979323846

// Define the differential equation function
int differential_equation(double t, const double z[], double dzdt[], void* params)
{   
    (void)(t); /* avoid unused parameter warning */
    int N = 4;
    int h = 0.5;
    double a=0;
    double b=0;
    double A=-1;
    double B=0;
    double C=-1;
    double R[N][2]; // Cast the parameters to a 2D array pointer
    double d = sqrt(2);

    for(int i = 0; i<N; i++){
        R[i][0]=cos(i*1.6)*d;
        R[i][1]=sin(i*1.6)*d;
    }

    for(int i=0; i < N; i++){

        double dn = pow((R[i][0]-z[0]),2) + pow((R[i][1]-z[1]),2) + 0.025;
        
        
        a = a + ((R[i][0]-z[0])/pow(dn,2.5));
        b = b + ((R[i][1]-z[1])/pow(dn,2.5));
        A = A -  (pow(h,2) - 4*pow((R[i][0]-z[0]),2) + (pow((R[i][1]-z[1]),2)/pow(dn,3.5)));
        B = B + 5*((R[i][1]-z[1])*(R[i][1]-z[1]))/pow(dn,3.5);
        C = C -  (pow(h,2) + pow((R[i][0]-z[0]),2) -4*pow((R[i][1]-z[1]),2))/pow(dn,3.5);
    }

    dzdt[0] = z[2];
    
    dzdt[1] = z[3];
    dzdt[2] = -z[0] + a;
    dzdt[3] = -z[1] + b;
    dzdt[4] = z[6];
    dzdt[5] = z[7];
    dzdt[6] = A * z[3] + B * z[4];
    dzdt[7] = B * z[3] + C * z[4];

    return GSL_SUCCESS;
}

int main(){
    

    clock_t start_time, end_time;
    double execution_time;
    start_time = clock();
    
    FILE* file =  fopen("result_1", "w");


    int N = 4 ;     // Nombre d'aimants
    double h = 0.5;   // rayon des aimants 
    double d = sqrt(2);      // Distance des aimants à l'origine
    double R[N][2]; // Position des aimants
    int T = 20; //Temps d'arrêt

    for(int i = 0; i<N; i++){
        R[i][0]=cos(i*1.6)*d;
        R[i][1]=sin(i*1.6)*d;
    }

    int K  = 10; // nombre de pas de discrétisation de la carte
    int x_min  = 0;
    int y_min  = 0;
    int x_max  = 2;
    int y_max  = 2;
    
    double delta_y = 0.2;
    double delta_x = 0.2;
    
    double CarteFLI[N+1][N+1];

    double Val_y[K]; // vecteur des valeurs de y
    double Val_x[K]; // vecteur des valeurs de x
    
    for (int i = 0; i < K; i++) {
        Val_y[i] = y_min + delta_y * i;
        Val_x[i] = x_min + delta_x * i ;
    }  

    
    for (int i=0; i<K; i++){
        double y = Val_y[i];
        
        for (int j=0; j<K; j++){
            double x = Val_x[j];
            double z[8] = {x, y, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5};
            double t = 0;
            double t_end = 40;

            gsl_odeiv2_system sys = {differential_equation, NULL, 8, NULL};
            gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-3, 1e-3, 0.0);
            
            for (int s = 1; s <= 40; s++){

                int status = gsl_odeiv2_driver_apply (driver, &t, t_end, z);

                if (status != GSL_SUCCESS){
                    printf ("error, return value=%d\n", status);
                    break;
                    }

                fprintf (file, "%.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", t, z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7]);
            }
            gsl_odeiv2_driver_free (driver);
        }
    }

    fclose(file); 

    end_time = clock();
    execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Execution time: %f seconds\n", execution_time);
    
}

