#include "potentialfunctionlj.hpp"
#include <stdexcept>

PotentialFunctionLJ::PotentialFunctionLJ(const double min_,
                                         const double max_):
  PotentialFunction("LJ",2,min_,max_) {
}

void PotentialFunctionLJ::setParam(string filename) {
  /* input parameter table should contain 2 entries
   * c12 and c6
   */
  votca::tools::Table param;
  param.Load(filename);

  if( param.size() != 4) {
    throw std::runtime_error("Potential parameters size mismatch!\n"
                             "Check input parameter file \""
                             + filename + "\" \nThere should be "
                             + boost::lexical_cast<string>(4) + " parameters");
  } else {
    _min = param.y(0);
    _cut_off = param.y(1);
    _lam(0) = param.y(2);
    _lam(1) = param.y(3);
  }
}

double PotentialFunctionLJ::CalculateF (const double r) const {
  if( r >= _min && r <= _cut_off)
    return _lam(0)/pow(r,12) - _lam(1)/pow(r,6);
  else
    return 0.0;
}

double PotentialFunctionLJ::CalculateDF (const double r) const {
  if( r >= _min && r <= _cut_off)
    {
      double du = -12.0*_lam(0)/pow(r,13) + 6.0*_lam(1)/pow(r,7);
      return du;
    } else
    return 0.0;
}

double PotentialFunctionLJ::CalculateIntR2F (const double r) const {
  if (r >= _min && r <= _cut_off)
      return -_lam(0)/(9.0*pow(r,9)) + _lam(1)/(3.0*pow(r,3));
  else return 0.0;
}
