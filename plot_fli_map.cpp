#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <matplotlibcpp.h>

namespace plt = matplotlibcpp;

void plotFLI_map(const std::string& filename, int N, int K, int T, double x_min, double x_max, double y_min, double y_max, double px, double py, double h, double time) {
    std::vector<std::vector<double>> CarteFLI(K + 1, std::vector<double>(K + 1));
    double delta_y = (y_max - y_min) / K;
    double delta_x = (x_max - x_min) / K;
    std::vector<double> Val_y, Val_x;

    for (int i = 0; i <= K; ++i) {
        Val_y.push_back(y_min + i * delta_y);
        Val_x.push_back(x_min + i * delta_x);
    }

    std::ifstream file(filename);
    std::string line;
    for (int i = 0; i <= K; ++i) {
        std::getline(file, line);
        std::istringstream iss(line);
        for (int j = 0; j <= K; ++j) {
            iss >> CarteFLI[i][j];
        }
    }
    file.close();

    plt::imshow(CarteFLI, {{"extent", {x_min, x_max, y_min, y_max}}, {"origin", "lower"}, {"aspect", "auto"}});
    plt::colorbar();
    plt::xlabel("$x$");
    plt::ylabel("$y$");
    plt::title("K=" + std::to_string(K) + " T=" + std::to_string(T) + " N=" + std::to_string(N) + " px=" + std::to_string(px) + " py=" + std::to_string(py) + " h=" + std::to_string(h) + " time=" + std::to_string(time));
    plt::show();
}
int main(){
    plotFLI_map("map.csv", 4, 150, 150, -1.7, 1.7, -1.7, 1.7, 0, 0, 0.5, 1231);
}