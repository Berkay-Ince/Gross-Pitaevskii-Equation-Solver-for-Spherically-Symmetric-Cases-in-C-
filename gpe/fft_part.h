#pragma once

#include <vector>
#include <cstddef>   // std::size_t
#include <fftw3.h>

struct FFTWCache {
    std::size_t N;
    double dr;
    double* in;
    fftw_complex* psi_k;
    fftw_plan fwd, inv;
    std::vector<double> k2;
};

namespace cGPE {

    FFTWCache make_cache(std::size_t N, double dr, int flags=FFTW_MEASURE);

    void destroy_cache(FFTWCache& c);

    void kinetic_step(std::vector<double>& psi, double dt, FFTWCache& c);

}