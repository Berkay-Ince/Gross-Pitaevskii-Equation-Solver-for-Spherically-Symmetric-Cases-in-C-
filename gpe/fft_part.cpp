#include <cmath>
#include "parameters.h"
#include "Common_functions.h"
#include "fft_part.h"

namespace cGPE {

FFTWCache make_cache(std::size_t N, double dr, int flags) {
    FFTWCache c;
    c.N = N; c.dr = dr;
    c.in = (double*)fftw_malloc(sizeof(double)*N);
    c.psi_k = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(N/2+1));

    double dk = 2.0*cGPE::pi/(dr*(double)N);
    c.k2.resize(N/2+1);
    for (std::size_t m=0; m<N/2+1; ++m) {
        double k = (double)m*dk;
        c.k2[m] = k*k;
    }

    c.fwd = fftw_plan_dft_r2c_1d((int)N, c.in, c.psi_k, flags);
    c.inv = fftw_plan_dft_c2r_1d((int)N, c.psi_k, c.in, flags);
    return c;
}

void destroy_cache(FFTWCache& c) {
    fftw_destroy_plan(c.fwd);
    fftw_destroy_plan(c.inv);
    fftw_free(c.in);
    fftw_free(c.psi_k);
}

void kinetic_step(std::vector<double>& psi, double dt, FFTWCache& c) {
    for (std::size_t i=0; i<c.N; ++i) c.in[i] = psi[i];
    fftw_execute(c.fwd);

    for (std::size_t m=0; m<c.N/2+1; ++m) {
        double damp = std::exp(-0.5 * c.k2[m] * dt);
        c.psi_k[m][0] *= damp;
        c.psi_k[m][1] *= damp;
    }

    fftw_execute(c.inv);

    double invN = 1.0/(double)c.N;
    for (std::size_t i=0; i<c.N; ++i) psi[i] = c.in[i]*invN;
}

}