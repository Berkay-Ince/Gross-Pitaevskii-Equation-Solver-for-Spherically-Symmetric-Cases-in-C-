#include "gpe/parameters.h"
#include "run_sweeps.h"
#include <iostream>
int main() {
    json params = cGPE::load_parameters("input.json");

    // Extract parameters

    // Constant parameters
    int N = params.at("numerics").at("N").get<int>();
    double tol = params.at("numerics").at("tol").get<double>();
    int max_iter = params.at("numerics").at("max_iter").get<int>();
    int renorm_every = params.at("numerics").at("renorm_every").get<int>();
    int report_every = params.at("numerics").at("report_every").get<int>();
    double E_aim = params.at("numerics").at("E_aim").get<double>();
    std::string out_dir = params.at("output").at("out_dir").get<std::string>();

    // Sweep spec
    run_sweeps sp;
    sp.g_values = params.at("g_list").get<std::vector<double>>();
    sp.methods = params.at("methods").get<std::vector<std::string>>();
    sp.R_values = params.at("R_list").get<std::vector<double>>();
    sp.dt_list_str = params.at("dt_list_str").get<std::vector<std::string>>();
    sp.dt_list_num = params.at("dt_list_num").get<std::vector<double>>();

    auto results = cGPE::run_sweep(sp, N, tol, max_iter, renorm_every, report_every, E_aim, out_dir);

    std::cout << "Total runs: " << results.size() << "\n";
    return 0;
}
