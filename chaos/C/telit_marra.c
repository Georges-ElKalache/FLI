#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <time.h>

#define N 4
#define T 20
#define K 10

double h, R[N][2], eps[N];

int fz(double t, const double z[], double dzdt[], void *params) {
    int i;
    double a = 0, b = 0, A = -1, B = 0, C = -1;
    
    for (i = 0; i < N; i++) {
        double dn = pow(R[i][0] - z[0], 2) + pow(R[i][1] - z[1], 2) + pow(h, 2);
        a += eps[i] * (R[i][0] - z[0]) / pow(dn, 2.5);
        b += eps[i] * (R[i][1] - z[1]) / pow(dn, 2.5);
        A -= eps[i] * (pow(h, 2) - 4 * pow(R[i][0] - z[0], 2) + pow(R[i][1] - z[1], 2)) / pow(dn, 3.5);
        B += 5 * eps[i] * (R[i][0] - z[0]) * (R[i][1] - z[1]) / pow(dn, 3.5);
        C -= eps[i] * (pow(h, 2) + pow(R[i][0] - z[0], 2) - 4 * pow(R[i][1] - z[1], 2)) / pow(dn, 3.5);
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

int main() {

    FILE* file =  fopen("result_chat_1", "w");

    clock_t start_time, end_time;
    double execution_time;

    // Start the clock
    start_time = clock();

    int i, j;
    double W0[4] = {0.5, 0.5, 0.5, 0.5};
    double x_min = 0, y_min = 0, x_max = 2, y_max = 2;
    double delta_y = (y_max - y_min) / K;
    double delta_x = (x_max - x_min) / K;
    double CarteFLI[K+1][K+1];
    double Val_y[K+1], Val_x[K+1];

    h = 0.5;
    double d = sqrt(2);
    for (i = 0; i < N; i++) {
        R[i][0] = cos(i * 2 * M_PI / N + M_PI / 4) * d;
        R[i][1] = sin(i * 2 * M_PI / N + M_PI / 4) * d;
        eps[i] = 1;
    }

    for (i = 0; i <= K; i++) {
        Val_y[i] = y_min + i * delta_y;
        for (j = 0; j <= K; j++) {
            double x = x_min + j * delta_x;
            double U0[4] = {x, Val_y[i], 0, 0};
            double Z0[8];

            for (int k = 0; k < 4; k++)
                Z0[k] = U0[k];
            for (int k = 0; k < 4; k++)
                Z0[k + 4] = W0[k];

            // GSL ODE Solver Setup
            const gsl_odeiv2_step_type *step_type = gsl_odeiv2_step_rkf45;
            gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(step_type, 8);
            gsl_odeiv2_control *control = gsl_odeiv2_control_standard_new(1e-3, 1e-3, 1.0, 0.2);
            gsl_odeiv2_evolve *evolve = gsl_odeiv2_evolve_alloc(8);
            gsl_odeiv2_system sys = {fz, NULL, 8, NULL};

            double t, t1 = T;
            double Z[8];
            for (int k = 0; k < 8; k++)
                Z[k] = Z0[k];

            // Solve ODE
            for (t = 0.0; t < t1; t += 0.1) {
                int status = gsl_odeiv2_evolve_apply(evolve, control, step, &sys, &t, t1, &t, Z);

                if (status != GSL_SUCCESS) {
                    printf("ODE solver error, status = %d\n", status);
                    break;
                }
            }

            double NW[4], FLI;
            for (int k = 0; k < T; k++) {
                for (int l = 0; l < 4; l++)
                    NW[l] = log(fabs(Z[l + 4]));
                FLI = NW[0];
                for (int l = 1; l < 4; l++)
                    FLI = fmax(FLI, NW[l]);
                CarteFLI[j][i] = FLI;
            }

            fprintf(file, " %f \n", FLI);
            // Free GSL ODE Solver Resources
            gsl_odeiv2_evolve_free(evolve);
            gsl_odeiv2_control_free(control);
            gsl_odeiv2_step_free(step);
        }
    }

    // Code for creating the figure goes here

  end_time = clock();

    // Calculate the execution time in seconds
execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
printf("Execution time: %f seconds\n", execution_time);
    
}