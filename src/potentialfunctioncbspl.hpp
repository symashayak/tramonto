
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
  PotentialFunctionCBSPL(const int nlam_, const double rexcl_,
                         const string& extrapol_type_,const double min_=0.0, const double max_=10.0);

  ~PotentialFunctionCBSPL(){}

  void setParam(string filename);
  void SaveParam(const string& filename);
  void setOptParam(const int i, const double val);
  double getOptParam(const int i) const;

  // calculate function value for given r
  double CalculateF (const double r) const;

  // calculate first derivative w.r.t. ith parameter
  double CalculateDF(const int i, const double r) const;

  // calculate second derivative w.r.t. ith parameter
  double CalculateD2F(const int i, const int j, const double r) const;

  // EQT specific functions

  int ContOptParamSize() const { return _lam.size() - _ncut - _nexcl; }

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
                           const double r1, const double r2) const {
    return 0.0;
  }

  // calculate integrant r.d2f/dlamdai*dlamdaj(r)
  double CalculateIntRd2F(const int i, const int j,
                          const double r1, const double r2) const {return 0.0;}

  // makes this potential and u C0 continuous at r
  void CutoffContinuous(const double uval);

  // extrapolate excluded parameters in the r<rmin region
  void extrapolExclParam();

  // first derivative w.r.t. r
  double CalculateDF(const double r) const;

  // double derivative w.r.t. r
  double CalculateD2F(const double r) const;

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

  double IntRF(const double r) const {}
  double IntR2F(const double r) const {}
  double IntRdF(const double r) const {}
  double IntRF(const int indx, const double r) const;
  double IntR2F(const int indx, const double r) const;
  double IntRdRF(const int indx, const double r) const;
  double IntRdF(const int i, const double r) const {}
  double IntRdF(const int i, const int indx, const double r) const;
  double IntRd2F(const int i, const int j, const double r) const {return 0;}
  double IntR2dF(const int i, const double r) const {}
  double IntR2dF(const int i, const int indx, const double r) const;
  double IntR2d2F(const int i, const int j, const double r) const {return 0;}

};

#endif
