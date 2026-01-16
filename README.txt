parameters.h includes   N,   R,  g,  tol,   max_iter,   pi
Common_function.h includes  Inf_spat_step,  Inf_time_step,  grid,   Initial_wavefunction,   Potential_function,    normalize_symmetric,    energy_expectation

g++ -std=c++17 \
  main.cpp run_sweeps.cpp \
  gpe/Common_functions.cpp gpe/fft_part.cpp gpe/for_npz_str_adj.cpp gpe/Methods.cpp gpe/parameters.cpp gpe/runner.cpp gpe/external/cnpy/cnpy.cpp \
  -Igpe -Igpe/external/cnpy \
  -lfftw3 -lz -lm -o run