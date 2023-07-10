#include <iostream>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

// Define the differential equation function
void differential_equation(const vector<double>& z, vector<double>& dzdt, double t)
{
    // Define the ODE system
    // dzdt[0] = dz/dt
    // dzdt[1] = d^2z/dt^2
    dzdt[0] = z[1];
    dzdt[1] = -z[0];
}

int main()
{
    // Define the initial conditions
    vector<double> z = { 1.0, 0.0 }; // Initial position and velocity

    // Define the time span
    double t_start = 0.0;
    double t_end = 10.0;
    double dt = 0.1;

    // Define the stepper
    runge_kutta4<vector<double>> stepper;

    // Integrate the ODE using the stepper
    for (double t = t_start; t < t_end; t += dt)
    {
        stepper.do_step(differential_equation, z, t, dt);
        cout << "t: " << t << ", z: " << z[0] << endl;
    }

    return 0;
}
