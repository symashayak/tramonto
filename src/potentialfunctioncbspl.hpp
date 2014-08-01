
#ifndef POTENTIALFUNCTIONCBSPL_H
#define	POTENTIALFUNCTIONCBSPL_H
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <math.h>
#include "potentialfunction.hpp"

using namespace std;
namespace ub = boost::numeric::ublas;

class PotentialFunctionCBSPL : public PotentialFunction {
 public:
  PotentialFunctionCBSPL(const int nlam_=0,
                         const double min_=0.0, const double max_=0.0);

  ~PotentialFunctionCBSPL(){}

  void setParam(string filename);
  void SaveParam(const string& filename);

  // calculate function value for given r
  double CalculateF (const double r) const;

  // first derivative w.r.t. r
  double CalculateDF(const double r) const;

  int getIndx( const double r) const;

  double getRbreak(const int i ) const { return _rbreak(i); }

 protected:

  double _dr;
  int _nexcl,_ncut;
  int _nbreak;
  double _rexcl;
  string _extrapol_type;

  ub::matrix<double> _M;
  ub::vector<double> _rbreak;
};

#endif
