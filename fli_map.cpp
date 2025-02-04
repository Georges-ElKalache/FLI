#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <matplotlibcpp.h>

double norm(const std::vector<double>& vector) {
    double sumSquares = 0.0;
    for (double value : vector) {
        sumSquares += value * value;
    }
    return std::sqrt(sumSquares);
}

const double pi = 3.14159265358979323846;
using namespace std;
using namespace boost::numeric::odeint;
typedef boost::array<double, 8> state_type;
const double h = 0.5;
const double d = sqrt(2.0);
const int N = 4;

struct fz_functor {
    template <typename T>
    void operator()(const T &z, T &dzdt, double t) const {
        double a = 0.0, b = 0.0, A = -1.0, B = 0.0, C = -1.0;
        vector<vector<double>> R(N, vector<double>(2));
        R[0][0] = 1; R[1][0] = 0; R[2][0] = 0; R[0][1] = 0; R[1][1] = 1; R[2][1] = -1; R[3][1] = 0; R[3][0] = -1;
        for (int i = 0; i <= N-1; ++i) {
            double dn = pow(R[i][0] - z[0], 2) + pow(R[i][1] - z[1], 2) + pow(h, 2);
            a += (R[i][0] - z[0]) / pow(dn, 2.5);
            b += (R[i][1] - z[1]) / pow(dn, 2.5);
            A -= 1 + (h * h - 4.0 * pow(R[i][0] - z[0], 2) + pow(R[i][1] - z[1], 2)) / pow(dn, 3.5);
            B += 5.0 * (R[i][0] - z[0]) * (R[i][1] - z[1]) / pow(dn, 3.5);
            C -= 1 + (h * h + pow(R[i][0] - z[0], 2) - 4.0 * pow(R[i][1] - z[1], 2)) / pow(dn, 3.5);
        }
        dzdt[0] = z[2]; dzdt[1] = z[3]; dzdt[2] = -z[0] + a; dzdt[3] = -z[1] + b; dzdt[4] = z[6];
        dzdt[5] = z[7]; dzdt[6] = A * z[4] + B * z[5]; dzdt[7] = B * z[4] + C * z[5];
    }
};

void FLImap(int K, int T, double x_min, double x_max, double y_min, double y_max, double px, double py) {
    clock_t start_time, end_time;
    double execution_time;
    start_time = clock();
    std::ofstream file("map.csv");
    double W0[4] = {0.5, 0.5, 0.5, 0.5};
    double delta_x = (x_max - x_min) / K;
    double delta_y = (y_max - y_min) / K;
    boost::multi_array<double, 2> CarteFLI(boost::extents[K + 1][K + 1]);
    typedef runge_kutta_cash_karp54<state_type> error_stepper_type;
    typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;
    fz_functor fz;
    controlled_stepper_type stepper;
    double FLI;
    for (int i = 0; i <= K; ++i) {
        double y = y_min + i * delta_y;
        double V[K];
        for (int j = 0; j <= K; ++j) {
            double x = x_min + j * delta_x;
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
                    std::log(norm(Z0[7]))
                };
                FLI = NW[0];
                for (double value : NW) {
                    if (value > FLI) {
                        FLI = value;
                    }
                }
                CarteFLI[j][i] = FLI;
            }
            V[j] = FLI;
            file << V[j] << " ";
        }
        file << "\n";
    }
    file.close();
    end_time = clock();
    execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Execution time: %f seconds\n", execution_time);
}


int main(){
    FLImap(150, 150, -1.7, 1.7, -1.7, 1.7, 0, 0);
}

