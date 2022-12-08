#pragma once

#include "algorithms.h"

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
  const double h;
};

template<int d>
Laplace<d>::Laplace(const unsigned int n)
  : n(n)
  , h(1./(n+1.))
{}

template<int d>
double
Laplace<d>::operator()(const unsigned int i, const unsigned int j) const
{
  double aij = 0.;

  if constexpr (d==1){
    assert(i<n && j<n);
    aij = (i==j) ? 2. : (i==(j+1) || j==(i+1)) ? -1. : 0.;
    }
  else if constexpr (d==2){
    assert(i<(n*n) && j<(n*n));
    aij = (i==j) ? 4. :
      (i==(j+1) || j==(i+1)) ? -1. :
      ((i==(j+n) || j==(i+n)) ? -1. : 0.);
    }
  else
    assert(false);

  return aij/(h*h);
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
