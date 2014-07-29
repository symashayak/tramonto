
#ifndef POTENTIALFUNCTIONQUAD_H
#define	POTENTIALFUNCTIONQUAD_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <math.h>
#include "potentialfunction.hpp"

using namespace std;

class PotentialFunctionQUAD : public PotentialFunction {
public:
  PotentialFunctionQUAD(const double min_=0.0, const double max_=10.0);
  ~PotentialFunctionQUAD(){}

  int ContOptParamSize() const { return 2;}

  // calculate function value for given r
  double CalculateF (const double r) const;

  // calculate first derivative w.r.t. ith parameter
  double CalculateDF(const int i, const double r) const;

  // calculate second derivative w.r.t. ith parameter
  double CalculateD2F(const int i, const int j, const double r) const;

  // EQT specific functions

  // calculate integrant r.f(r)
  double CalculateIntRF(const double r1, const double r2) const;

  /// \brief calculate integration r^2 * f(r) dr from r1 to r2
  double CalculateIntR2F(const double r1, const double r2) const;

  double CalculateIntRdF(const double r1, const double r2) const;

  // calculate integrant r.df/dlamdai(r)
  double CalculateIntRdF(const int i, const double r1,
                         const double r2) const;

  /// \brief calculate r^2.df/dlamdai(r)
  double CalculateIntR2dF(const int i, const double r1,
                          const double r2) const ;

  /// \brief calculate r^2.d2f/dlamdai*dlamdaj(r)
  double CalculateIntR2d2F(const int i, const int j,
                           const double r1, const double r2) const;

  // calculate integrant r.d2f/dlamdai*dlamdaj(r)
  double CalculateIntRd2F(const int i, const int j,
                          const double r1, const double r2) const;

  // makes this potential and u C0 continuous at r
  void CutoffContinuous(const double uval);

  // first derivative w.r.t. r
  double CalculateDF(const double r) const;

  // double derivative w.r.t. r
  double CalculateD2F(const double r) const;

protected:

  double IntRF(const double r) const;
  double IntR2F(const double r) const;
  double IntRdF(const double r) const;
  double IntRdF(const int i, const double r) const;
  double IntRd2F(const int i, const int j, const double r) const;
  double IntR2dF(const int i, const double r) const;
  double IntR2d2F(const int i, const int j, const double r) const;

};

#endif
