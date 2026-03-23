#include <vector>
#include <cmath>
#include "parameters.h"
#include "Common_functions.h"
#include <fftw3.h>
#include "thomas_algorithm.h"
#include "fft_part.h"

namespace cGPE {

    // Forward Euler Time Stepping Method for Radially Symmetric GPE

void FE_Method(std::vector<double> &psi,
                                        const std::vector<double> &V,
                                        const std::vector<double> &r, 
                                        double dr, double dt, double g)
    {
        std::size_t N = psi.size();
        std::vector<double> laplacian_psi = cGPE::laplacian(psi, dr);
        for (int i = 0; i <N; i++){
            psi[i] -= dt * (-0.5 * laplacian_psi[i] + V[i] * psi[i] + g * psi[i] *psi[i] * psi[i] / (r[i]*r[i]+ 1e-30));
        }
        psi [0] = 0.0; // Boundary condition at r=-R
        psi [N-1] = 0.0; // Boundary condition at r=R
    }

    // Time-Splitting Spectral Method for Radially Symmetric GPE 

void TSSM_Step(std::vector<double> &psi,
                                    const std::vector<double> &V,
                                    const std::vector<double> &r,
                                    double dr, double dt, double g
                                    ,FFTWCache& c)
{
    const std::size_t N = psi.size();

    // First step
    for(int i = 0; i<N; i++){
        double damping_part = std::exp(-dt * (V[i] + g * psi[i] * psi[i] / (r[i]*r[i]+ 1e-30))*0.5);
        psi[i] *= damping_part;
    }

    // Kinetic evolution
    cGPE::kinetic_step(psi,dt, c);

    // Second step
    for(int i = 0; i<N; i++){
        double damping_part = std::exp(-dt * (V[i] + g * psi[i] * psi[i] / (r[i]*r[i]+ 1e-30))*0.5);
        psi[i] *= damping_part;
    }
}

void CN_Step(std::vector<double> &psi,
                    const std::vector<double> &V,
                    const std::vector<double> &r,
                    double dr, double dt, double g)
{
    const std::size_t N = psi.size();
    double alpha = dt / (4 * dr * dr);
    std::vector<double> diagonal(N);
    for (std::size_t i = 0; i < N; ++i) {
        diagonal[i] =1.0 + 2.0 * alpha + 0.5 * dt * (V[i] + g * psi[i] * psi[i] / (r[i]*r[i]+ 1e-30));
    }
    std::vector<double> lower_diag(N - 1, -alpha);
    std::vector<double> upper_diag(N - 1, -alpha);
    std::vector<double> rhs(N);
    for (std::size_t i = 1; i < N - 1; ++i) {
        rhs[i] = alpha * psi[i - 1] + alpha * psi[i + 1]+ (2.0 -diagonal[i]) * psi[i];}
    std::vector<double> psi_new = thomas::thomas_algorithm(lower_diag, diagonal, upper_diag, rhs);
    psi = psi_new;
}
// RK4 Time Stepping Method for Radially Symmetric GPE

// Hamiltonian for the GPE
std::vector<double> hamiltonian(const std::vector<double> &psi,
                                const std::vector<double> &V,
                                const std::vector<double> &r,
                                double dr, double g)
{
    std::vector<double> H_psi(psi.size());
    std::vector<double> lap = cGPE::laplacian(psi, dr);
    for (std::size_t i = 0; i < psi.size(); ++i) {
        H_psi[i] = 0.5 * lap[i] - V[i] * psi[i] - g * psi[i] * psi[i] * psi[i]/ (r[i]*r[i]+ 1e-30);
    }
    return H_psi;
}
void RK4_Step(std::vector<double> &psi,
                    const std::vector<double> &V,
                    const std::vector<double> &r,
                    double dr, double dt, double g)
{
        std::vector<double> k1(psi.size()), k2(psi.size()), k3(psi.size()), k4(psi.size());
        std::vector<double> derv_temp(psi.size()), derv_temp2(psi.size()), derv_temp3(psi.size()), derv_temp4(psi.size());
        int n = psi.size();
        std::vector<double> psi_temp(n);
        derv_temp = hamiltonian(psi, V, r, dr, g);
        for(int j=0;j < n; ++j)
        {
            k1[j] = dt * derv_temp[j];
            psi_temp[j] = psi[j] + k1[j] / 2.0;
        }
        psi_temp[0] = 0.0; psi_temp[n-1] = 0.0;
        derv_temp2 = hamiltonian(psi_temp, V, r, dr, g);
        for(int j=0;j < n; ++j)
        {
            k2[j] = dt * derv_temp2[j];
            psi_temp[j] = psi[j] + k2[j] / 2.0;
        }
        psi_temp[0] = 0.0; psi_temp[n-1] = 0.0;
        derv_temp3 = hamiltonian(psi_temp, V, r, dr, g);
        for(int j=0;j < n; ++j)
        {
            k3[j] = dt * derv_temp3[j];
            psi_temp[j] = psi[j] + k3[j];
        }
        psi_temp[0] = 0.0; psi_temp[n-1] = 0.0;
        derv_temp4 = hamiltonian(psi_temp, V, r, dr, g);
        for(int j=0;j < n; ++j)
        {
            k4[j] =dt * derv_temp4[j];
        }
        for(int j=0;j < n; ++j)
        {
            psi[j] += (1.0 / 6.0) * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
        }
        psi[0] = 0.0;
        psi[n-1] = 0.0;
}
}
