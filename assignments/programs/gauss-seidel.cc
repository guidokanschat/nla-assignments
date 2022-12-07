// Compile with C++17:
// g++ -std=c++17 -o gauss-seidel gauss-seidel.cc

#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

class MatrixBase
{
public:
  virtual ~MatrixBase() = default;
  virtual double operator()(const unsigned int, const unsigned int) const = 0;
  virtual std::vector<double> vmult(const std::vector<double>&) const = 0;
};

double norm(const std::vector<double>& v)
{
  double y = 0.;
  for (const auto& x : v)
    y += x*x;
  return sqrt(y);
}

void print(const std::vector<double>& v)
{
  for (const auto& x : v)
    std::cout << x << " ";
  std::cout << "\n";
}

std::vector<double> operator-(const std::vector<double>& v,
			      const std::vector<double>& w)
{
  const unsigned int n = w.size();
  assert(v.size()==n);
  std::vector<double> z(n);
  for (unsigned int i=0; i<n; ++i)
    z[i] = v[i]-w[i];
  return z;
}


template<int d>
class Laplace : public MatrixBase
{
public:
  Laplace(const unsigned int n);
  double operator()(const unsigned int i, const unsigned int j) const override;
  std::vector<double> vmult(const std::vector<double>& v) const override;
  void print() const;
  int size() const;
private:
  const unsigned int n;
};

template<int d>
Laplace<d>::Laplace(const unsigned int n)
  : n(n)
{}

template<int d>
double
Laplace<d>::operator()(const unsigned int i, const unsigned int j) const
{
  if constexpr (d==1){
    assert(i<n && j<n);
    return (i==j) ? 2. : (i==(j+1) || j==(i+1)) ? -1. : 0.;
    }
  else if constexpr (d==2){
    assert(i<(n*n) && j<(n*n));
    return (i==j) ? 4. :
      (i==(j+1) || j==(i+1)) ? -1. :
      ((i==(j+n) || j==(i+n)) ? -1. : 0.);
    }
  else
    assert(false);
}

template<int d>
std::vector<double>
Laplace<d>::vmult(const std::vector<double>& v) const
{
  const unsigned int N = this->size();
  assert(v.size()==N);

  std::vector<double> w(N);
  for (unsigned int i=0; i<N; ++i)
    {
      w[i] = 0.;
      for (unsigned int j=0; j<N; ++j)
	w[i] += (*this)(i,j)*v[j];
    }
  return w;
}

template<int d>
void
Laplace<d>::print() const
{
  for (unsigned int i=0; i<this->size(); ++i){
    for (unsigned int j=0; j<this->size(); ++j)
      printf("%2.f ", (*this)(i,j));
    std::cout << "\n";
  }
}

template<int d>
int
Laplace<d>::size() const
{
  return n;
}


void gauss_seidel(const MatrixBase& A,
		  const std::vector<double>& b,
		  std::vector<double>& x,
		  const unsigned int n_steps = 10)
{
  // std::cout << "x = ";
  // print(x);
  const double r = norm(A.vmult(x) - b);
  std::cout << "Starting residual: r = " << r << " \n";

  const unsigned int n = x.size();
  std::vector<double> y(n, 0.);
  for (unsigned int m = 0; m < n_steps; m++)
    {
      for (unsigned int i = 0; i < n; i++)
	{
	  y[i] = b[i] / A(i,i);
	  for (unsigned int j = 0; j < n; j++)
	    {
	      if (j == i)
	  	continue;
	      y[i] = y[i] - ((A(i,j) / A(i,i)) * x[j]);
	      x[i] = y[i];
	    }
	}
      // std::cout << "x = ";
      // print(x);
      const double r = norm(A.vmult(x) - b);
      std::cout << "step = " << m << "    r = " << r << " \n";
    }
  std::cout << "\n";
}


int main()
{
  const unsigned int d = 2;
  const unsigned int n = 20;
  const unsigned int n_dofs = (d==1) ? n : (d==2) ? n*n : 0;

  if (n_dofs==0) assert(false);

  Laplace<d> A(n_dofs);
  if (n_dofs<8) A.print();

  const std::vector<double> b(n_dofs, 1.);
  std::vector<double> x(n_dofs, 0.);

  const unsigned int n_steps = 50;
  gauss_seidel(A, b, x, n_steps);
}
