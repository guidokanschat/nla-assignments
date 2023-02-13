#pragma once

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

class MatrixBase
{
public:
  virtual ~MatrixBase() = default;
  virtual double operator()(const unsigned int, const unsigned int) const = 0;
  virtual std::vector<double> vmult(const std::vector<double>&) const = 0;
};

template<typename Vector>
void print(const Vector& v)
{
  for (const auto& x : v)
    std::cout << x << " ";
  std::cout << "\n";
}

template<typename Vector>
Vector operator+(const Vector& v,
		 const Vector& w)
{
  const unsigned int n = w.size();
  assert(v.size()==n);
  Vector z(n);
  for (unsigned int i=0; i<n; ++i)
    z[i] = v[i]+w[i];
  return z;
}

template<typename Vector>
Vector operator-(const Vector& v,
		 const Vector& w)
{
  const unsigned int n = w.size();
  assert(v.size()==n);
  Vector z(n);
  for (unsigned int i=0; i<n; ++i)
    z[i] = v[i]-w[i];
  return z;
}

template<typename Vector>
Vector& operator+=(Vector& v,
		   const Vector& w)
{
  const unsigned int n = v.size();
  assert(w.size()==n);
  for (unsigned int i=0; i<n; ++i)
    v[i] += w[i];
  return v;
}

template<typename Vector>
Vector& operator-=(Vector& v,
		   const Vector& w)
{
  const unsigned int n = v.size();
  assert(w.size()==n);
  for (unsigned int i=0; i<n; ++i)
    v[i] -= w[i];
  return v;
}

template<typename Vector>
double operator*(const Vector& v,
		 const Vector& w)
{
  const unsigned int n = w.size();
  assert(v.size()==n);
  double x = 0.;
  for (unsigned int i=0; i<n; ++i)
    x += v[i]*w[i];
  return x;
}

template<typename Vector>
Vector operator*(const double& x,
		 const Vector& v)
{
  const unsigned int n = v.size();
  std::vector<double> z(n);
  for (unsigned int i=0; i<n; ++i)
    z[i] = x*v[i];
  return z;
}

template<typename Vector>
double norm(const Vector& v)
{
  double y = 0.;
  for (const auto& x : v)
    y += x*x;
  return sqrt(y);
}

void gauss_seidel(const MatrixBase& A,
		  const std::vector<double>& b,
		  std::vector<double>& x,
		  const unsigned int n_steps)
{
  const double r = norm(A.vmult(x) - b);
  std::cout << "Starting residual: r = " << r << " \n";
  const unsigned int n = x.size();
  std::vector<double> y(n, 0.);
  for (unsigned int m = 0; m < n_steps; m++){
      for (unsigned int i = 0; i < n; i++){
        y[i] = b[i] / A(i,i);
        for (unsigned int j = 0; j < n; j++)
          {
            if (j == i)
              continue;
            y[i] = y[i] - ((A(i,j) / A(i,i)) * x[j]);
            x[i] = y[i];
          }
      }
      const double r = norm(A.vmult(x) - b);
      std::cout << "step = " << m << "    r = " << r << " \n";
    }
  std::cout << "\n";
}


void steepest_decent(const MatrixBase& A,
		     const std::vector<double>& b,
		     std::vector<double>& x,
		     const double tol,
		     const unsigned int n_steps)
{
  std::vector<double> r = b;
  r -= A.vmult(x);
  const unsigned int n = x.size();
  std::vector<double> p(n);
  for (unsigned int m = 0; m < n_steps; m++)
    {
      const double r2 = r*r;
      std::cout << "step = " << m << "    r = " << sqrt(r2) << " \n";
      if ( r2 < tol*tol)
	return;
      p = A.vmult(r);
      const double alpha = r2/(p*r);
      x += alpha*r;
      r -= alpha*p;
    }
}


void minres(const MatrixBase& A,
	    const std::vector<double>& b,
	    std::vector<double>& x,
	    const double tol,
	    const unsigned int n_steps)
{
  std::vector<double> r = b;
  r -= A.vmult(x);
  const unsigned int n = x.size();
  std::vector<double> p(n);
  for (unsigned int m = 0; m < n_steps; m++)
    {
      const double r2 = r*r;
      std::cout << "step = " << m << "    r = " << sqrt(r2) << " \n";
      if ( r2 < tol*tol)
	return;
      p = A.vmult(r);
      const double alpha = r2/(p*p);
      x += alpha*r;
      r -= alpha*p;
    }
}
