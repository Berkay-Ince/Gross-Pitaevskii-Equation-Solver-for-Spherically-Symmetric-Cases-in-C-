#pragma once

#include <string>
#include <vector>
#include <cstdint>

struct GPEResult;   // forward declaration

namespace cGPE{

std::string sanitize_token(std::string s);
std::string to_compact(double x);

void save_result_npz(const std::string& out_dir,
                     const GPEResult& res,
                     const std::string& dt_input,
                     double dt_value,
                     int N);
}
