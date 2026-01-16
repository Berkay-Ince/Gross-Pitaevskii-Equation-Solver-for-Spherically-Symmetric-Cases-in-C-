#pragma once
#include <vector>

// Thomas algorithm for solving tridiagonal systems
namespace thomas{

    inline std::vector<double> thomas_algorithm(const std::vector<double> &a,
                                                const std::vector<double> &b,
                                                const std::vector<double> &c,
                                                const std::vector<double> &r )
    {
        int n = r.size();
        std::vector<double> x(n); //what we want to find

        std::vector<double> gamma(n-1);
        gamma[0]=c[0]/b[0];

        std::vector<double> rho(n);
        rho[0]=r[0]/b[0];

        // Decomposition and forward substitution
        for (int i=1; i<n; ++i)
        {
            double denom = b[i]-a[i-1]*gamma[i-1];
            rho[i]=(r[i]-a[i-1]*rho[i-1])/denom;
            if(i<n-1){
                gamma[i]=c[i]/denom;}
        }

        // Back substitution
        x[n-1]=rho[n-1];
        for (int i=1; i<n; ++i)
        {
            x[n-1-i]=rho[n-1-i]-gamma[n-1-i]*x[n-i];
        }
        return x;
}
}