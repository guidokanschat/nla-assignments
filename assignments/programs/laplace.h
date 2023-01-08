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
  void print_sparsity() const;
  int size() const; // Returns the number of unknowns n_dofs
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
  assert(i<this->size() && j<this->size());

  if constexpr (d==1){
    aij = (i==j) ? 2. : (i==(j+1) || j==(i+1)) ? -1. : 0.;
    }
  else if constexpr (d==2){
    aij = (i==j) ? 4. :
      ( (i==(j+1) && (j%n+1)!=n) || (j==(i+1) && (i%n+1)!=n) ) ? -1. :
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
void
Laplace<d>::print_sparsity() const
{
  for (unsigned int i=0; i<this->size(); ++i){
    for (unsigned int j=0; j<this->size(); ++j)
      printf("%1s ", ((*this)(i,j) != 0. ? "x" : ""));
    std::cout << "\n";
  }
}

template<int d>
int
Laplace<d>::size() const
{
  return (d==1) ? n : (d==2) ? n*n : 0;
}
