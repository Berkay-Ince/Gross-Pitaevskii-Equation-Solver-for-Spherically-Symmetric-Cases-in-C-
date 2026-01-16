#pragma once
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

using json = nlohmann::json;

namespace cGPE
{
    json load_parameters(const std::string& filename);
    double parse_dt_expr(const std::string& expr, double dr);
    std::vector<double> get_dt_str_to_num(const std::vector<std::string>& dt_list_str, double dr);
}
