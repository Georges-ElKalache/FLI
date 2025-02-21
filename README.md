# Fast Lyapunov Indicator for chaos quantative analysis in dynamical systems 

## Overview
**FLI** is a C++ program that models the dynamics of the Magnetic Pendulum, considering gravitational, magnetic, and frictional forces. It computes the **Fast Lyapunov Indicator** on the system in an attempt to quantify sensitivity to initial conditions and chaotic behavior. FLI is a numerical technique used to study chaos in dynamical systems.

The code is based on the following [paper]([https://github.com/Georges-ElKalache/chaos-quantative-analysis/blob/main/chaos%20in%20the%20magnetic%20pendulum_compressed.pdf])
## Features
- Implements the Runge-Kutta Cash-Karp method for numerical integration  
- Uses Boost ODEInt library for solving differential equations  
- Outputs simulation results to a CSV file for further analysis  
- Efficient execution in C++ with performance tracking  

## Installation

### Prerequisites
Ensure you have the following installed on your system:
- C++ compiler (GCC, Clang, or MSVC)
- Boost C++ Libraries (ODEInt)
- Matplotlibcpp

### Clone the Repository
```sh
git clone https://github.com/Georges-ElKalache/FLI.git
cd FLI
```

## Compilation & Execution

### Using g++
```sh
g++ -std=c++11 -O2 -I /path/to/boost fli_map.cpp -o fli_map
./fli_map

g++ -std=c++11 -O2 -I /path/to/boost plot_fli_map.cpp -o plot_fli_map
./plot_fli_map

g++ -std=c++11 -O2 -I /path/to/boost energy_section.cpp -o energy_section
./energy_section

g++ -std=c++11 -O2 -I /path/to/boost plot_energy_section.cpp -o plot_energy_section
./plot_energy_section
```
Replace `/path/to/boost` with the actual path where Boost is installed.

## File Descriptions

### 1. `fli_map.cpp`
This file generates a **Fast Lyapunov Indicator (FLI) Map**, which solves the system's differential equations and stores results in `map.csv`.

#### **Function: `FLImap`**
**Parameters:**
- `K` - Grid resolution for initial conditions
- `T` - Time span of simulation
- `x_min, x_max, y_min, y_max` - Range for initial conditions
- `px, py` - Initial momenta values

### 2. `plot_fli_map.cpp`
This file reads `map.csv` and visualizes the **FLI map** using `matplotlibcpp`.

#### **Function: `plotFLI_map`**
**Parameters:**
- `filename` - The CSV file with computed FLIs
- `N` - Number of magnets
- `K, T` - Grid size and time span
- `x_min, x_max, y_min, y_max` - Plotting ranges
- `px, py, h, time` - Initial conditions and parameters

### 3. `energy_section.cpp`
This file computes the **FLI energy section**, a different perspective on chaos by analyzing energy variations.

#### **Function: `FLIsection`**
**Parameters:**
- `H` - Energy level
- `h` - Regulatory parameter
- `N` - Number of magnets
- `K, T` - Grid resolution and simulation time
- `py_min, py_max, y_min, y_max` - Ranges for phase space
- `px, py, x` - Initial conditions

### 4. `plot_energy_section.cpp`
This file reads `section.csv` and visualizes the **FLI energy section** using `matplotlibcpp`.

#### **Function: `plotFLI_section`**
**Parameters:**
- `filename` - The CSV file with computed FLIs
- `N, H, K, T` - Number of magnets, energy level, grid size, and time span
- `x_min, x_max, y_min, y_max` - Plotting ranges
- `x, h, time` - Initial conditions

### Output
- The program generates a CSV file (`file2.csv`) containing computed values.
- Execution time is displayed in the terminal.


## Usage
Modify `main.cpp` to set initial conditions, parameters, and integration time as needed. The program allows for:
- Adjusting the number and positions of magnets
- Changing simulation duration and step size
- Testing different values of friction and perturbations

## Results & Applications
The generated FLI maps confirm the chaotic nature of the system, highlighting sensitive and stable regions. This methodology has broader applications, including:
- **Astrodynamics** (stability of satellite orbits)
- **Nonlinear finance models**
- **Neuroscience** (complex neuronal networks)


## References
For more theoretical background, see the paper *"Chaos in the Magnetic Pendulum"* by **Jérôme Perez & Georges ElKalache** (ENSTA Paris, 2024).

