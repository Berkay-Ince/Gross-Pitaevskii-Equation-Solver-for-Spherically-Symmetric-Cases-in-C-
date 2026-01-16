#pragma once
#include <vector>
#include <string>
#include "gpe/runner.h"
struct run_sweeps
{
    std::vector<std::string> methods;
    std::vector<double> g_values;
    std::vector<std::string> dt_list_str;
    std::vector<double> dt_list_num;
    std::vector<double> R_values;
    
};

namespace cGPE
{
    std::vector<GPEResult> run_sweep(const run_sweeps& sp,
                                int N,
                                double tol,
                                int max_iter,
                                int renorm_every,
                                int report_every,
                                double E_aim,
                                const std::string& out_dir);
}
