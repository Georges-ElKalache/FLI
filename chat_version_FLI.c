#include <stdio.h>
#include <math.h>
#include <time.h>

#define N 4
#define T 20
#define K 15

double h, R[N][2], eps[N];

void fz(double t, double z[], double dzdt[]) {
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
}

void ode23s(void (*f)(double, double*, double*), double tspan[], double z0[], int n, double t[], double z[][8]) {
    int i;
    double h = (tspan[1] - tspan[0]) / (n - 1);
    
    for (i = 0; i < n; i++) {
        t[i] = tspan[0] + i * h;
        f(t[i], z[i], z[i]);
        
        if (i < n - 1) {
            double k1[8], k2[8], k3[8], k4[8], temp[8];
            int j;
            
            for (j = 0; j < 8; j++)
                k1[j] = h * z[i][j];
            
            for (j = 0; j < 8; j++)
                temp[j] = z[i][j] + 0.5 * k1[j];
            f(t[i] + 0.5 * h, temp, k2);
            
            for (j = 0; j < 8; j++)
                temp[j] = z[i][j] + 0.5 * k2[j];
            f(t[i] + 0.5 * h, temp, k3);
            
            for (j = 0; j < 8; j++)
                temp[j] = z[i][j] + k3[j];
            f(t[i] + h, temp, k4);
            
            for (j = 0; j < 8; j++)
                z[i + 1][j] = z[i][j] + (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6.0;
        }
    }
}

int main() {

    clock_t start_time, end_time;
    double execution_time;

    // Start the clock
    start_time = clock();

    FILE* file =  fopen("result_chat", "w");

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

            double tspan[2] = {0, T};
            double t[T], Z[T][8];
            ode23s(fz, tspan, Z0, T, t, Z);

            double NW[4], FLI;
            for (int k = 0; k < T; k++) {
                for (int l = 0; l < 4; l++)
                    NW[l] = log(fabs(Z[k][l + 4]));
                FLI = NW[0];
                for (int l = 1; l < 4; l++)
                    FLI = fmax(FLI, NW[l]);
                CarteFLI[j][i] = FLI;
            
            }
            //fprintf(file, " pour x = %f, y = %f => FLI = %f \n", x, Val_y[i], FLI);
            fprintf(file, " %f \n", FLI);
        }
    }

    // Code for creating the figure goes here

  end_time = clock();

    // Calculate the execution time in seconds
execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
printf("Execution time: %f seconds\n", execution_time);
    
}
