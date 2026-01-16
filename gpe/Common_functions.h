#pragma once
#include <vector>

namespace cGPE 
{
    inline constexpr double pi = 3.141592653589793238462643383279502884;
    double Inf_spat_step(double R, int N);
    std::vector<double> grid(double R, double dr, int N);
    std::vector<double> Initial_wavefunction(const std::vector<double> &r, double R);
    std::vector<double> Potential_function(const std::vector<double> &r);
    void normalize_symmetric(std::vector<double> &psi, double dr);
    double energy_expectation(const std::vector<double> &psi,const std::vector<double> &V,const std::vector<double> &r,double g, double dr);
    std::vector<double> laplacian(const std::vector<double> &psi, double dr);
}