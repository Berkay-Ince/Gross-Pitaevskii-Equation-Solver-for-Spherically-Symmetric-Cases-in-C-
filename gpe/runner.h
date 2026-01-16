#pragma once
#include <vector>
#include <string>
#include <chrono>

struct GPEResult {
    std::string method;        // method used
    double g;              // interaction strength
    double dt;             // time step
    int N;         // number of radial grid points
    double R;              // radial domain size
    int step;          // actual steps taken
    double E_final;             // final expectation energy
    std::chrono::duration<double> wall_time;        // wall time taken
    std::vector<double> r;          // radial grid points
    std::vector<double> psi_final;  // final wavefunction
};

struct variables{
    std::vector<double>& psi;
    const std::vector<double>& V;
    const std::vector<double>& r;
    double dr; double dt; double g;
    const std::string& method; int max_iter;
    const int renorm_every;
    const int report_every;
    const double tol;
    const double E_aim;
    int N; 
    double R;
};

namespace cGPE {

    GPEResult run_simulation(variables var);

} 
