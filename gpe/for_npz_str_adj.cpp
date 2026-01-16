#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "external/cnpy/cnpy.h"
#include "runner.h"   // for GPEResult definition

namespace fs = std::filesystem;

namespace cGPE {

static std::string replace_all(std::string s, char from, char to)
{
    std::replace(s.begin(), s.end(), from, to);
    return s;
}

static void npz_save_string(const std::string& npz_path,
                            const std::string& key,
                            const std::string& value,
                            const std::string& mode)
{
    std::vector<uint8_t> bytes(value.begin(), value.end());
    cnpy::npz_save(npz_path, key, bytes.data(), {bytes.size()}, mode);
}

std::string sanitize_token(std::string s)
{
    s = replace_all(s, '.', '_');
    s = replace_all(s, '*', '_');
    s = replace_all(s, ' ', '_');
    s = replace_all(s, '/', '_');
    return s;
}

std::string to_compact(double x)
{
    std::ostringstream oss;
    oss << std::setprecision(16) << x;
    return oss.str();
}

void save_result_npz(const std::string& out_dir,
                     const GPEResult& res,
                     const std::string& dt_input,
                     double dt_value,
                     int N)
{
    fs::create_directories(out_dir);

    std::string fname =
        sanitize_token(res.method)
        + "_g" + sanitize_token(to_compact(res.g))
        + "_R" + sanitize_token(to_compact(res.R))
        + "_"  + sanitize_token(dt_input)
        + ".npz";

    std::string path = (fs::path(out_dir) / fname).string();

    // arrays
    cnpy::npz_save(path, "r",
                   res.r.data(), {res.r.size()}, "w");

    cnpy::npz_save(path, "psi_final",
                   res.psi_final.data(), {res.psi_final.size()}, "a");

    // scalars
    double energy = res.E_final;
    double wall_time_s = res.wall_time.count();
    int steps = res.step;

    cnpy::npz_save(path, "energy", &energy, {1}, "a");
    cnpy::npz_save(path, "steps", &steps, {1}, "a");
    cnpy::npz_save(path, "wall_time_s", &wall_time_s, {1}, "a");

    double g = res.g, R = res.R, dtv = dt_value;
    cnpy::npz_save(path, "g", &g, {1}, "a");
    cnpy::npz_save(path, "R", &R, {1}, "a");
    cnpy::npz_save(path, "dt_value", &dtv, {1}, "a");
    cnpy::npz_save(path, "N", &N, {1}, "a");

    // strings
    npz_save_string(path, "method", res.method, "a");
    npz_save_string(path, "dt_input", dt_input, "a");

    std::cout << "Saved: " << path << "\n";
}

}