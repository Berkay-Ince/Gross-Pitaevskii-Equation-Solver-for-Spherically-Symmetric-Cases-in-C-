
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>   // for std::cerr
#include "Common_functions.h"
namespace cGPE {

    // Infinitesimal spatial step size
double Inf_spat_step(double R, int N)
{
        return 2.0 * R / N;
}

// Grid: from approximately -R to +R (cell-centered)
std::vector<double> grid(double R, double dr, int N)
{
        std::vector<double> r(N);
        for (int j = 0; j < N; ++j) {
            r[j] = (j - N/2) * dr;
        }
        return r;
}

    // Initial wavefunction psi(r) = r * exp(-r^2)
std::vector<double> Initial_wavefunction(const std::vector<double> &r, double R)
{
        std::vector<double> psi(r.size());
        for (std::size_t i = 0; i < r.size(); ++i) {
            psi[i] = r[i] * std::exp(-r[i] * r[i]);
        }
        return psi;
}

    // Potential V(r) = 0.5 r^2
std::vector<double> Potential_function(const std::vector<double> &r)
{
        std::vector<double> V(r.size());
        for (std::size_t i = 0; i < r.size(); ++i) {
            V[i] = 0.5 * r[i] * r[i];
        }
        return V;
}

    // Normalization for 3D radial wavefunction (psi = r * phi)
void normalize_symmetric(std::vector<double> &psi, double dr)
{
    std::size_t N = psi.size();
    std::size_t mid = N/2; 
    
    double norm = 0.0;

    for (std::size_t i = mid+1; i < N; ++i) {
        norm += psi[i] * psi[i] * dr;
    }
    norm = std::sqrt(4.0 * cGPE::pi * norm);
    for (std::size_t i = 0; i < psi.size(); ++i) {
        psi[i] /= (norm+1e-30);
    }
}

    // Expectation value of energy
double energy_expectation(const std::vector<double> &psi,
                                     const std::vector<double> &V, 
                                     const std::vector<double> &r,
                                     double g, double dr)
{
        double E_kin = 0.0;
        double E_pot = 0.0;
        double E_int = 0.0;

        std::size_t N = psi.size();
        std::size_t mid = N / 2;
        for (std::size_t i =mid+1; i < N-1; ++i) {
            double dpsi = (psi[i+1] - psi[i-1]) / (2.0 * dr);
            E_kin += 0.5 * 4.0 * cGPE::pi * dpsi * dpsi * dr;
            E_pot += 4.0 * cGPE::pi * V[i] *psi[i] * psi[i] * dr;
            if (r[i] != 0.0) {
                E_int += 2.0 * cGPE::pi * g * std::pow(psi[i], 4)
                         / (r[i] * r[i]) * dr;
            }
        }
        return E_kin + E_pot + E_int;
}

std::vector<double> laplacian(const std::vector<double> &psi, double dr)
    {
        std::size_t N = psi.size();
        std::vector<double> kinetic_term(N);

        // interior points
        for (std::size_t i = 1; i < N-1; ++i) {
            kinetic_term[i] = (psi[i+1] - 2.0*psi[i] + psi[i-1]) / (dr*dr);
        }

        // boundaries
        kinetic_term[0]   = (psi[1]     - 2.0*psi[0])     / (dr*dr);
        kinetic_term[N-1] = (-2.0*psi[N-1] + psi[N-2])    / (dr*dr);

        return kinetic_term;
    } 
}
