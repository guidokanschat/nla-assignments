// Compile with C++17:
// g++ -std=c++17 -o gauss-seidel gauss-seidel.cc

#include "laplace.h"

int main()
{
  const unsigned int d = 2;
  const unsigned int n = 100;
  const unsigned int n_dofs = (d==1) ? n : (d==2) ? n*n : 0;

  if (n_dofs==0) assert(false);

  Laplace<d> A(n_dofs);
  if (n_dofs<8) A.print();

  const std::vector<double> b(n_dofs, 1.);
  std::vector<double> x(n_dofs, 0.);

  const unsigned int n_steps = 50;
  gauss_seidel(A, b, x, n_steps);
}
