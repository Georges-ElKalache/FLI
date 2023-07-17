#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#define pi 3.14159265358979323846

float entree() {
    float x;
    printf("Donner un nombre: ");
    scanf("%f",&x);
    return x;
}



// Definition de l'equation differentielle
int differential_equation(double t, const double z[], double dzdt[], void* params)
{   
    (void)(t); 
    double lambda = 1.00;

    double A = -1 -2*lambda*z[1];
    //printf("A = %f\n", A);
    double B = -2*lambda*z[0];
    double C = -1 -2*lambda*z[1];
    //printf("B = %f\n", B);
    //printf("C = %f\n", C);

    dzdt[0] = z[2];
    dzdt[1] = z[3];
    dzdt[2] = -z[0] - lambda*z[0]*z[1];
    dzdt[3] = -z[1] - lambda*(z[0]*z[0] - z[1]*z[1]);
    dzdt[4] = z[6];
    dzdt[5] = z[7];
    dzdt[6] = A * z[3] + B * z[4];
    dzdt[7] = B * z[3] + C * z[4];

    return GSL_SUCCESS;
}


int main(int argc, char *argv[]){

    
    clock_t start_time, end_time;
    double execution_time;
    start_time = clock();
    
    FILE* file =  fopen("result_1", "w");
    /*
    int K  = 10; // nombre de pas de discrétisation de la carte
    int x_min  = 0;
    int y_min  = 0;
    
    double delta_y = (double) 2/K;
    double delta_x = (double) 2/K;
    double Val_y[K]; // vecteur des valeurs de y
    double Val_x[K]; // vecteur des valeurs de x
    
    for (int i = 0; i < K; i++) {
        Val_y[i] = y_min + delta_y * i;
        Val_x[i] = x_minc + delta_x * i ;
    }  

    //int f = 0;
    
    for (int i=0; i<K; i++){
        double y = Val_y[i];
        
        for (int j=0; j<K; j++){
            double x = Val_x[j];
            double z[8] = {x, y, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5};
            double t = 0;
            double t_end = 20;

            gsl_odeiv2_system sys = {differential_equation, NULL, 8, NULL};
            gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);
            
            for (int s = 1; s <= t_end; s++){
                int t1 = 100;
                double ti = s * t1 / 100.0;
                int status = gsl_odeiv2_driver_apply (driver, &t, ti, z);
               
                if (status != GSL_SUCCESS){
                    printf ("error, return value=%d\n", status);
                    break;
                    }
                fprintf (file," %f %f %f %f %f %f %f %f\n", z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7]);
            }
            f++;
            gsl_odeiv2_driver_free (driver);
        }
    }
    */
    double x = entree();
    double y = entree();
    double z[8] = {x, y, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5};
    double t = 0;
    double t_end = 40;

    gsl_odeiv2_system sys = {differential_equation, NULL, 8, NULL};
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);

    for (int s = 1; s <= t_end; s++)
    {
        int t1 = 100;
        double ti = s * t1 / 100.0;
        int status = gsl_odeiv2_driver_apply(driver, &t, ti, z);

        if (status != GSL_SUCCESS)
        {
            printf("error, return value=%d\n", status);
            break;
        }
        fprintf(file, " %f  %f %f %f %f %f %f %f %f\n", t, z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7]);
    }
    gsl_odeiv2_driver_free (driver);
    fclose(file); 
    end_time = clock();
    execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Execution time: %f seconds\n", execution_time);
}

/*
    for(int q = 1; q<8; q++){
        Z[s-1 + f*s][q-1]= z[q];
    }

double NW[4], FLI;
for (int k = 0; k < t_end; k++) {
    for (int l = 0; l < 4; l++)
        NW[l] = log(fabs(Z[k][l + 4]));
    FLI = NW[0];
    for (int l = 1; l < 4; l++)
        FLI = fmax(FLI, NW[l]);
    CarteFLI[j][i] = FLI;

}
    */