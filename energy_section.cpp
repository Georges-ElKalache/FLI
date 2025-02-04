#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <ctime>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;

typedef std::vector<double> state_type;

void FLIsection(int H, int h, int N, int K, int T, double py_min, double py_max, double y_min, double y_max, double px, double py, double x) {
    std::ofstream file("section.csv");
    double delta_py = (py_max - py_min) / K;
    double delta_y = (y_max - y_min) / K;
    double CarteFLI[K][K];
    double start_time = clock();

    for (int i = 0; i <= K; ++i) {
        double V[K];
        double y = y_min + i * delta_y;
        for (int j = 0; j <= K; ++j) {
            double py = py_min + j * delta_py;
            double Q = 0;
            for (int l = 0; l < N; l++) {
                Q += (double)1 / (pow(pow((x - R[l][0]), 2) + pow((y - R[l][1]), 2) + pow(h, 2), (3 / 2)));
            }
            double px2 = (double)2 * H - (pow(py, 2) + pow(x, 2) + pow(y, 2)) + (0.667 * Q);
            double FLI = 0;
            if (px2 >= 0) {
                double px = sqrt(px2);
                state_type U0 = {x, y, px, py, W0[0], W0[1], W0[2], W0[3]};
                state_type Z0 = U0;
                double t = 0;
                double ti = 1;
                while (ti <= T) {
                    integrate_adaptive(stepper, fz, Z0, t, ti, 0.5);
                    ti = ti + 1;
                    t = t + 1;
                    std::vector<double> NW = {
                        std::log(norm(Z0[4])),
                        std::log(norm(Z0[5])),
                        std::log(norm(Z0[6])),
                        std::log(norm(Z0[7]))};
                    FLI = NW[0];
                    for (double value : NW) {
                        if (value > FLI) {
                            FLI = value;
                        }
                    }
                }
            } else {
                FLI = 0;
            }
            CarteFLI[j][i] = FLI;
            V[j] = FLI;
            file << V[j] << " ";
        }
        file << "\n";
    }
    file.close();
    double end_time = clock();
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Execution time: %f seconds\n", execution_time);
}

int main() {
    FLIsection(5, 0.5, 4, 500, 100, -4, 4, -4, 4, 0, 0, 0);
    }

