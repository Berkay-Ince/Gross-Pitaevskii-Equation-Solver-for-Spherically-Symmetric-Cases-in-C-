#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
using json = nlohmann::json;

namespace cGPE 
{
    json load_parameters(const std::string& filename)
    {
        std::ifstream file(filename);
        if (!file) {
            throw std::runtime_error("Could not open file: " + filename);
        }

        json parameters;
        file >> parameters;
        return parameters;
    };
    double parse_dt_expr(const std::string& expr, double dr)
    {
        if (expr == "dr*dr") {
            return dr * dr;
        }

        const std::string suffix = "*dr*dr";

        if (expr.size() > suffix.size() &&
            expr.substr(expr.size() - suffix.size()) == suffix)
        {
            std::string coef_str = expr.substr(0, expr.size() - suffix.size());
            double coef = std::stod(coef_str);
            return coef * dr * dr;
        }

        throw std::runtime_error("Invalid dt expression: " + expr);
    };
    std::vector<double> get_dt_str_to_num(const std::vector<std::string>& dt_list_str, double dr) 
    {
        std::vector<double> dt_list_num;
        for (const auto& expr : dt_list_str) {
            double dt = parse_dt_expr(expr, dr);
            dt_list_num.push_back(dt);
        }
        return dt_list_num;
    }

}