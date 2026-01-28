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
}
