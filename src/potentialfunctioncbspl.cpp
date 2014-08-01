#include "potentialfunctioncbspl.hpp"
#include <stdexcept>

PotentialFunctionCBSPL::PotentialFunctionCBSPL(const int nlam_,
                                               const double min_,
                                               const double max_):
  PotentialFunction("CBSPL",nlam_,min_,max_) {
  _M.resize(4,4,false);
  _M.clear();
  _M(0,0) =  1.0; _M(0,1) =  4.0; _M(0,2) =  1.0; _M(0,3) = 0.0;
  _M(1,0) = -3.0; _M(1,1) =  0.0; _M(1,2) =  3.0; _M(1,3) = 0.0;
  _M(2,0) =  3.0; _M(2,1) = -6.0; _M(2,2) =  3.0; _M(2,3) = 0.0;
  _M(3,0) = -1.0; _M(3,1) =  3.0; _M(3,2) = -3.0; _M(3,3) = 1.0;
  _M /= 6.0;

}

void PotentialFunctionCBSPL::setParam(string filename) {

  // Table param;
  // param.Load(filename);

  // _lam.clear();

  // if( param.size() != _lam.size()) {

  //   throw std::runtime_error("Potential parameters size mismatch!\n"
  //                            "Check input parameter file \""
  //                            + filename + "\" \nThere should be "
  //                            + boost::lexical_cast<string>( _lam.size() ) + " parameters");
  // } else {

  //   for( int i = 0; i < _lam.size(); i++){

  //     _rbreak(i) = param.x(i);
  //     _lam(i) = param.y(i);

  //   }

  // }
  int nlam;// = param.size();
  _lam.resize(nlam);
  _lam.clear();
  _nbreak = nlam - 2;
  _dr = (_cut_off )/( double (nlam - 3) );

  // break point locations
  // since ncoeff = nbreak +2 , r values for last two coefficients are also
  // computed
  _rbreak.resize(nlam, false);
  _rbreak.clear();

  for (int i = 0; i < nlam; i++)
    _rbreak(i) = i * _dr;
}

void PotentialFunctionCBSPL::SaveParam(const string& filename){

  // Table param;
  // param.SetHasYErr(false);
  // param.resize(_lam.size(), false);

  // for (int i = 0; i < _lam.size(); i++)
  //   param.set(i, _rbreak(i), _lam(i), 'i');

  // param.Save(filename);

}

int PotentialFunctionCBSPL::getIndx(const double r) const {

  return min( int( r /_dr ) , int(_lam.size()) - 4 );

}

double PotentialFunctionCBSPL::CalculateF (const double r) const {

  if ( r >= _min && r < _cut_off) {

    ub::vector<double> R;
    ub::vector<double> B;

    int indx = getIndx(r);
    double rk = indx*_dr;
    double t = ( r - rk)/_dr;

    R.resize(4,false); R.clear();
    R(0) = 1.0; R(1) = t; R(2) = t*t; R(3) = t*t*t;

    ub::vector<double> RM = ub::prod(R,_M);

    B.resize(4,false); B.clear();

    B(0) = _lam(indx); B(1) = _lam(indx+1); B(2) = _lam(indx+2);
    B(3) = _lam(indx+3);

    double u = ub::inner_prod(B,RM);

    return u;

  } else
    return 0.0;

}

double PotentialFunctionCBSPL::CalculateDF (const double r) const {

  if ( r >= _min && r <= _cut_off) {

    ub::vector<double> R;
    ub::vector<double> B;

    int indx = getIndx(r);

    double rk = indx*_dr;

    double t = ( r - rk)/_dr;

    R.resize(4,false); R.clear();

    R(0) = 0.0; R(1) = 1.0/_dr; R(2) = 2.0*t/_dr; R(3) = 3.0*t*t/_dr;

    ub::vector<double> RM = ub::prod(R,_M);

    B.resize(4,false); B.clear();

    B(0) = _lam(indx); B(1) = _lam(indx+1); B(2) = _lam(indx+2);
    B(3) = _lam(indx+3);

    double u = ub::inner_prod(B,RM);

    return u;

  } else {

    return 0.0;
  }
}
