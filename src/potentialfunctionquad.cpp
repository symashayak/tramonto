
#include "potentialfunctionquad.hpp"

PotentialFunctionQUAD::PotentialFunctionQUAD(const double min_,
                                             const double max_):
  PotentialFunction("QUAD",3,min_,max_)
{
}

void PotentialFunctionQUAD::setParam(string filename) {
  /* input parameter table should contain 5 entries
   * rmin,rcut,a0,a1,a2
   */
  votca::tools::Table param;
  param.Load(filename);

  if( param.size() != 5) {
    throw std::runtime_error("Potential parameters size mismatch!\n"
                             "Check input parameter file \""
                             + filename + "\" \nThere should be "
                             + boost::lexical_cast<string>(5) + " parameters");
  } else {
    _min = param.y(0);
    _cut_off = param.y(1);
    _lam(0) = param.y(2);
    _lam(1) = param.y(3);
    _lam(2) = param.y(4);
  }
}

double PotentialFunctionQUAD::CalculateF(const double r) const {

  if ( r >= _min && r < _cut_off ) {

    return  _lam(0)*(r - _cut_off)*(r - _cut_off)
      + _lam(1)*(r - _cut_off)
      + _lam(2);

  } else {

    return 0.0;

  }

}

double PotentialFunctionQUAD::CalculateDF(const double r) const {

  if ( r >= _min && r < _cut_off ) {

    return  2.0*_lam(0)*(r - _cut_off)
      + _lam(1);

  } else {

    return 0.0;

  }

}

double PotentialFunctionQUAD::CalculateIntR2F(const double r) const {

  if (r >= _min && r <= _cut_off) {
    return (1.0/30.0)*_lam(0)*r*r*r*(6.0*r*r - 15.0*_cut_off*r +
                                     10.0*_cut_off*_cut_off)
      + (1.0/12.0)*_lam(1)*r*r*r*( 3.0*r - 4.0*_cut_off )
      + (1.0/3.0)*_lam(2)*r*r*r;
  } else {
    return 0.0;
  }

}
