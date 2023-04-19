#include <cassert>
#include <vector>

std::vector<double>
tridiag_solve(const std::vector<double>& a,
              const std::vector<double>& b,
              const std::vector<double>& c,
              const std::vector<double>& d) {

    // solve the tridiagonal system of the form Ax = d:
    //
    //   a_i x_{i-1} + b_i x_i + c_i x_{i+1} = d_i
    //
    //     for i = 0, ... , N-1 with a_0 = 0 and c_{N-1} = 0
    //
    //   In matrix form:
    //       b is the main diagonal
    //       a is the subdiagonal
    //       c is the superdiagonal

    int N = static_cast<int>(b.size());

    assert(static_cast<int>(a.size()) == N);
    assert(static_cast<int>(c.size()) == N);
    assert(static_cast<int>(d.size()) == N);

    std::vector<double> cprime(N, 0.0);
    std::vector<double> dprime(N, 0.0);

    std::vector<double> x(N, 0.0);

    cprime[0] = c[0] / b[0];
    dprime[0] = d[0] / b[0];

    for (int i = 0; i < N-1; ++i) {
        cprime[i] = c[i] / (b[i] - cprime[i-1] * a[i]);
        dprime[i] = (d[i] - dprime[i-1] * a[i]) /
            (b[i] - cprime[i-1] * a[i]);
    }

    // there is no cprime for N-1

    dprime[N-1] = (d[N-1] - dprime[N-2] * a[N-1]) /
        (b[N-1] - cprime[N-2] * a[N-1]);

    // back substitution

    x[N-1] = dprime[N-1];
    for (int i = N - 2; i >= 0; --i) {
        x[i] = dprime[i] - cprime[i] * x[i+1];
    }

    return x;
}


std::vector<double>
tridiag_Ax(const std::vector<double>& a,
           const std::vector<double>& b,
           const std::vector<double>& c,
           const std::vector<double>& x) {

    // multiply the tridiagonal matrix A by vector x and return the
    // product vector d */

    int N = static_cast<int>(b.size());

    assert(static_cast<int>(a.size()) == N);
    assert(static_cast<int>(c.size()) == N);
    assert(static_cast<int>(x.size()) == N);

    std::vector<double> d(N, 0.0);

    d[0] = b[0] * x[0] + c[0] * x[1];
    for (int i = 1; i < N-1; ++i) {
        d[i] = a[i] * x[i-1] + b[i] * x[i] + c[i] * x[i+1];
    }
    d[N-1] = a[N-1] * x[N-2] + b[N-1] * x[N-1];

    return d;

}

