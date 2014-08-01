
#include "potentialfunctionquad.hpp"

PotentialFunctionQUAD::PotentialFunctionQUAD(const double min_,
                                             const double max_):
  PotentialFunction("QUAD",3,min_,max_)
{
}

void PotentialFunctionQUAD::setParam(string filename) {

  // Table param;
  // param.Load(filename);

  // if( param.size() != _lam.size()) {

  //     throw std::runtime_error("Potential parameters size mismatch!\n"
  //             "Check input parameter file \""
  //             + filename + "\" \nThere should be "
  //             + boost::lexical_cast<string>( _lam.size() ) + " parameters");
  // } else {
  //     for( int i = 0; i < _lam.size(); i++)
  //         _lam(i) = param.y(i);

  // }
  // _min = param.y(0);
  // _max = param.y(1);
  // _lam(0) = param.y(2);
  // _lam(1) = param.y(3);
  // _lam(2) = param.y(4);

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
