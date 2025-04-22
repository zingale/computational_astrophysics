#include <cmath>
#include <iostream>
#include <vector>

#include "tridiag.H"

int main() {

  // test things out

  int N = 10;

  std::vector<double> a(N, 0.0);
  std::vector<double> b(N, 0.0);
  std::vector<double> c(N, 0.0);
  std::vector<double> x(N, 0.0);

  // we'll create a test matrix with b = 4, a = c = 1
  for (int i = 0; i < N; ++i) {
      a[i] = 1.0;
      b[i] = 4.0;
      c[i] = 1.0;
  }

  // we'll make x just be a sequence of ints
  for (int i = 0; i < N; ++i) {
      x[i] = i;
  }

  // get d from our matrix multiply
  auto d = tridiag_Ax(a, b, c, x);

  // now recover x from our tridiag solve
  auto x_solve = tridiag_solve(a, b, c, d);

  for (int i = 0; i < N; ++i) {
      std::cout << x[i] << " " << std::abs(x[i] - x_solve[i]) << std::endl;
  }

}

