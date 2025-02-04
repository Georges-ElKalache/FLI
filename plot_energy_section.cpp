#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <matplotlibcpp.h>

namespace plt = matplotlibcpp;
void plotFLI_section(const std::string &filename, int N, int H, int K, int T, double x_min, double x_max, double y_min, double y_max, double x, double h, double time) {
    std::vector<std::vector<double>> CarteFLI(K + 1, std::vector<double>(K + 1));
    std::ifstream file(filename);
    for (int i = 0; i <= K; ++i) {
        for (int j = 0; j <= K; ++j) {
            file >> CarteFLI[i][j];
        }
    }
    file.close();
    std::vector<double> Val_x(K + 1), Val_y(K + 1);
    double delta_x = (x_max - x_min) / K;
    double delta_y = (y_max - y_min) / K;
    for (int i = 0; i <= K; ++i) {
        Val_x[i] = x_min + i * delta_x;
        Val_y[i] = y_min + i * delta_y;
    }
    plt::figure_size(800, 600);
    plt::imshow(CarteFLI, {{"extent", {x_min, x_max, y_min, y_max}}, {"origin", "lower"}, {"aspect", "auto"}});
    plt::colorbar();
    plt::xlabel("$py$", {{"fontsize", 14}, {"fontweight", "bold"}, {"color", "k"}});
    plt::ylabel("$y$", {{"fontsize", 14}, {"fontweight", "bold"}, {"color", "k"}});
    plt::title("H=" + std::to_string(H) + ", K=" + std::to_string(K) + " T=" + std::to_string(T) + " N=" + std::to_string(N) + " x=" + std::to_string(x) + " h=" + std::to_string(h) + " time=" + std::to_string(time));
    plt::show();
}
int main() {
    plotFLI_section("section.csv", 4, 5, 500, 100, -4, 4, -4, 4, 0, 0.5, 0);
    return 0;
}