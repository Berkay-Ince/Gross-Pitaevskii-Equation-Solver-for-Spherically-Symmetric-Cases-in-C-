#include <vector>
#include <string>
#include <iostream>
#include <chrono>
#include "parameters.h"
#include "Methods.h"
#include "fft_part.h"
#include "Common_functions.h"
#include "runner.h"

namespace cGPE {

    GPEResult run_simulation(variables var)
    {
        // Initialize FFTW cache for TSSP method
        int breakin_step = 0;
        //  Start Timers
        auto t_start = std::chrono::steady_clock::now();
        FFTWCache fftw_cache;
        if (var.method == "TSSM") { // TSSP method
            fftw_cache = make_cache(var.N, var.dr);
        }
        double E_prev = cGPE::energy_expectation(var.psi, var.V, var.r, var.g,var.dr);
        for (int step = 0; step < var.max_iter; ++step) {
            breakin_step = step+1;
            if (var.method == "FE") {
                cGPE::FE_Method(var.psi, var.V, var.r, var.dr, var.dt, var.g);
            } else if (var.method == "CN") {
                cGPE::CN_Step(var.psi, var.V, var.r, var.dr, var.dt, var.g);
            } else if (var.method == "TSSM") {
                cGPE::TSSM_Step(var.psi, var.V, var.r, var.dr, var.dt, var.g, fftw_cache);
            } else {
                throw std::runtime_error("Unknown method selected.");
            }
            if(step % var.renorm_every == 0) {
                cGPE::normalize_symmetric(var.psi,var.dr);
            }
            if(step % var.report_every == 0) {
                double E_curr = cGPE::energy_expectation(var.psi, var.V, var.r, var.g,var.dr);
                double rel = std::abs((E_curr - E_prev)/(E_prev + 1e-30));
                if (rel < var.tol) {
                    break;
                }
                if (E_curr < var.E_aim) {
                    break;
                }
                E_prev = E_curr;
            }
            
        }
        // Clean up FFTW cache
        if (var.method == "TSSM") {
            destroy_cache(fftw_cache);
        }
        auto t_end = std::chrono::steady_clock::now();
        std::chrono::duration<double> wall_time = t_end - t_start;
        double E_final = cGPE::energy_expectation(var.psi, var.V, var.r, var.g,var.dr);
        std::vector<double> psi_final = var.psi;
        return {var.method,
            var.g,
            var.dt,
            var.N,
            var.R,
            breakin_step,
            E_final,
            wall_time,
            var.r,
            psi_final};
    }
}
