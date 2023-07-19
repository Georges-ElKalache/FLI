#include <iostream>
#include <cmath>
#include <fstream>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::array<double, 8> state_type;

const double h = 0.5;
const double d = sqrt(2.0);
const int N = 4;

struct fz_functor
{
    template<typename T>
    void operator()(const T& z, T& dzdt, double t) const
    {
        double a = 0.0, b = 0.0, A = -1.0, B = 0.0, C = -1.0;
        
        for (int i = 0; i < N; ++i)
        {
            T R;
            R[0] = cos(i * 2.0 * M_PI / N + M_PI / 4.0) * d;
            R[1] = sin(i * 2.0 * M_PI / N + M_PI / 4.0) * d;
            double dn = pow(R[0] - z[0], 2) + pow(R[1] - z[1], 2) + pow(h, 2);
            a += (R[0] - z[0]) / pow(dn, 2.5);
            b += (R[1] - z[1]) / pow(dn, 2.5);
            A -= (h * h - 4.0 * pow(R[0] - z[0], 2) + pow(R[1] - z[1], 2)) / pow(dn, 3.5);
            B += 5.0 * (R[0] - z[0]) * (R[1] - z[1]) / pow(dn, 3.5);
            C -= (h * h + pow(R[0] - z[0], 2) - 4.0 * pow(R[1] - z[1], 2)) / pow(dn, 3.5);
        }
        
        dzdt[0] = z[2];
        dzdt[1] = z[3];
        dzdt[2] = -z[0] + a;
        dzdt[3] = -z[1] + b;
        dzdt[4] = z[6];
        dzdt[5] = z[7];
        dzdt[6] = A * z[4] + B * z[5];
        dzdt[7] = B * z[4] + C * z[5];
    }
};

int main()
{

    std::ofstream file("file1");
    double T = 20.0;
    double W0[4] = { 0.5, 0.5, 0.5, 0.5 };
    int K = 10;
    double x_min = 0.0, y_min = 0.0, x_max = 2.0, y_max = 2.0;
    double delta_x = (x_max - x_min) / K;
    double delta_y = (y_max - y_min) / K;
    boost::multi_array<double, 2> CarteFLI(boost::extents[K + 1][K + 1]);
    
    typedef runge_kutta_cash_karp54<state_type> error_stepper_type;
    typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;
    
    fz_functor fz;
    controlled_stepper_type stepper;
    
    for (int i = 0; i <= K; ++i)
    {
        double y = y_min + i * delta_y;
        
        for (int j = 0; j <= K; ++j)
        {
            double x = x_min + j * delta_x;
            state_type U0 = { x, y, 0.0, 0.0, W0[0], W0[1], W0[2], W0[3] };
            state_type Z0 = U0;
            
            integrate_adaptive(stepper, fz, Z0, 0.0, T, 0.01);
            file << t << " " << z[0] << " " << z[1] << " " << z[2] << " " << z[3] << " " << z[4] << " " << z[5] << " " << z[6] << " " << z[7] << "\n";

        }
    }
    file.close();
    return 0;
}