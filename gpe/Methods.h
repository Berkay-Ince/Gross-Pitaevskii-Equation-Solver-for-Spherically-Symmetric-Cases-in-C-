#pragma once
#include <vector>
#include "fft_part.h"

namespace cGPE {

void FE_Method(std::vector<double> &psi,
                                    const std::vector<double> &V,
                                    const std::vector<double> &r, 
                                    double dr, double dt, double g);
    
void TSSM_Step(std::vector<double> &psi,
                                const std::vector<double> &V,
                                const std::vector<double> &r,
                                double dr, double dt, double g
                                ,FFTWCache& c);
    
void CN_Step(std::vector<double> &psi,
                    const std::vector<double> &V,
                    const std::vector<double> &r,
                    double dr, double dt, double g);

void RK4_Step(std::vector<double> &psi,
                    const std::vector<double> &V,
                    const std::vector<double> &r,
                    double dr, double dt, double g);       
}
