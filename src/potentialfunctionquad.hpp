
#ifndef POTENTIALFUNCTIONQUAD_H
#define	POTENTIALFUNCTIONQUAD_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <math.h>
#include "potentialfunction.hpp"

using namespace std;

class PotentialFunctionQUAD : public PotentialFunction {
public:
  PotentialFunctionQUAD(const double min_=0.0, const double max_=0.0);
  ~PotentialFunctionQUAD(){}

  void setParam(string filename);
  // calculate function value for given r
  double CalculateF (const double r) const;
  // first derivative w.r.t. r
  double CalculateDF(const double r) const;
  /// \brief calculate integrant r^2 * f(r) dr
  double CalculateIntR2F(const double r) const;

protected:
};

#endif
