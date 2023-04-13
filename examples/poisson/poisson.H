#ifndef POISSON_H
#define POISSON_H

#include <vector>
#include <cassert>
#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>
#include <functional>

///
/// Solve the 1-d Poisson problem on a node-centered finite-difference
/// grid with Dirichlet BCs using Gauss-Seidel smoothing
///
class Poisson {

private:

    double xmin;
    double xmax;

    int N;

    std::vector<double> x;
    std::vector<double> phi;
    std::vector<double> f;

    double dx;

public:

    Poisson(double xmin_in, double xmax_in, int N_in)
        : xmin{xmin_in}, xmax{xmax_in}, N{N_in},
          x(N_in, 0), phi(N_in, 0), f(N_in, 0)
    {

        // initialize the coordinates
        assert (xmax > xmin);

        dx = (xmax - xmin) / static_cast<double>(N-1);

        for (int i = 0; i < N; ++i) {
            x[i] = xmin + static_cast<double>(i) * dx;
        }

    }

    ///
    /// set the left Dirichlet boundary condition
    ///
    void set_left_bc(double val) {phi[0] = val;}

    ///
    /// set the right Dirichlet boundary condition
    ///
    void set_right_bc(double val) {phi[N-1] = val;}

    ///
    /// set the source term, f, in L phi = f
    ///
    void set_source(std::function<double(double)> func) {

        for (int i = 0; i < N; ++i) {
            f[i] = func(x[i]);
        }
    }

    ///
    /// return the number of points
    ///
    int npts() {return N;}

    ///
    /// return the coordinate vector
    ///
    const std::vector<double>& xc() {return x;}

    ///
    /// return the source vector
    ///
    std::vector<double>& get_source() {
        return f;
    }

    ///
    /// return the solution vector
    ///
    std::vector<double>& get_phi() {
        return phi;
    }

    ///
    /// do Gauss-Seidel smoothing for n_smooth iterations
    ///
    void smooth(int n_smooth) {

        // perform Gauss-Seidel smoothing

        // we only operate on the interior nodes

        for (int i = 0; i < n_smooth; ++i) {

            for (int j = 1; j < N-1; ++j) {
                phi[j] = 0.5 * (phi[j-1] + phi[j+1] - dx * dx * f[j]);
            }
        }
    }

    ///
    /// solve the Poisson problem via relaxation until the residual
    /// norm is tol compared to the source norm
    ///
    void solve(double tol) {

        double err = std::numeric_limits<double>::max();

        double source_norm = norm(f);

        while (err > tol) {

            smooth(10);

            double r_norm = norm(residual());

            if (source_norm != 0.0) {
                err = r_norm / source_norm;
            } else {
                err = r_norm;
            }
        }

    }

    ///
    /// compute the residual
    ///
    std::vector<double> residual() {
        std::vector<double> r(N, 0);

        for (int i = 1; i < N-1; ++i) {
            r[i] = f[i] - (phi[i+1] - 2.0 * phi[i] + phi[i-1]) / (dx * dx);
        }

        return r;

    }

    ///
    /// given a vector e on our grid, return the L2 norm
    ///
    double norm(const std::vector<double>& e) {
        double l2{0.0};

        for (int i = 0; i < N; ++i) {
            l2 += e[i] * e[i];
        }

        l2 = std::sqrt(dx * l2);

        return l2;
    }

    ///
    /// output the solution to file fname
    ///
    void write_solution(const std::string& fname) {
        auto of = std::ofstream(fname);

        for (int i = 0; i < N; ++i) {
            of << std::setw(20) << x[i] << " " << std::setw(20) << phi[i] << std::setw(20) << f[i] << std::endl;
        }

    }


};
#endif
