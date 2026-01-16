#include "gpe/parameters.h"
#include "gpe/runner.h"
#include "gpe/Common_functions.h"
#include "gpe/for_npz_str_adj.h"
#include <vector>
#include <string>
#include <cstddef>
#include "run_sweeps.h"

namespace cGPE
{

std::vector<GPEResult> run_sweep(const run_sweeps& sp,
                                int N,
                                double tol,
                                int max_iter,
                                int renorm_every,
                                int report_every,
                                double E_aim,
                                const std::string& out_dir)
{
    std::vector<GPEResult> out;

    out.reserve(sp.methods.size() * sp.g_values.size() * sp.R_values.size() * (sp.dt_list_str.size() + sp.dt_list_num.size()));
    
    for (double R : sp.R_values) {
        
        double dr = cGPE::Inf_spat_step(R, N);
        std::vector<double> grid_sweep =cGPE::grid(R, dr, N);
        std::vector<double> dt_conv = cGPE::get_dt_str_to_num(sp.dt_list_str, dr);
        std::vector<double> pot_sweep = cGPE::Potential_function(grid_sweep);

        std::vector<std::pair<std::string,double>> dt_pairs;
        dt_pairs.reserve(sp.dt_list_str.size() + sp.dt_list_num.size());
        
        for (size_t i = 0; i < sp.dt_list_str.size(); ++i) {
            dt_pairs.push_back({sp.dt_list_str[i], dt_conv[i]});
        }
        for (double dt_num : sp.dt_list_num) {
            dt_pairs.push_back({cGPE::to_compact(dt_num), dt_num}); // to_compact from your filename helper
        }
        for (const auto& method : sp.methods) {
            
            for (double g : sp.g_values) {
   
                
                    for (const auto& [dt_input, dt_value] : dt_pairs) {
                        std::vector<double> psi_sweep = cGPE::Initial_wavefunction(grid_sweep, R);
                        
                        variables var
                        {
                            psi_sweep,
                            pot_sweep,
                            grid_sweep,
                            dr,
                            dt_value,
                            g,
                            method,
                            max_iter,
                            renorm_every,
                            report_every,
                            tol,
                            E_aim,
                            N,
                            R
                        };
                        GPEResult res = cGPE::run_simulation(var);
                        out.push_back(res);
                        cGPE::save_result_npz(out_dir, res, dt_input, dt_value, N);
                }
            }
        }
    }
    return out;
}
}
