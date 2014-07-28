#include "potentialfunctioncbspl.hpp"
#include <stdexcept>

PotentialFunctionCBSPL::PotentialFunctionCBSPL(const int nlam_, const double rexcl_,
                                               const string& extrapol_type_,const double min_, const double max_) :
  PotentialFunction("CBSPL",nlam_,min_,max_) {

  _nbreak = _lam.size() - 2;
  _dr = (_cut_off )/( double (_lam.size() - 3) );

  // break point locations
  // since ncoeff = nbreak +2 , r values for last two coefficients are also
  // computed
  _rbreak.resize(_lam.size(), false);
  _rbreak.clear();

  for (int i = 0; i < _lam.size(); i++)
    _rbreak(i) = i * _dr;

  _rexcl = rexcl_;
  _extrapol_type = extrapol_type_;

  // exclude knots corresponding to r < _min
  _nexcl = min( int( ( _rexcl )/_dr ), _nbreak - 2 );

  // account for finite numerical division of _min/_dr
  // e.g. 0.24/0.02 may result in 11.99999999999999
  //if( _rbreak(_nexcl) == _rexcl ) _nexcl++;

  // 4 is an appropriate number hence hard-coded
  _ncut = 4;

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

}

void PotentialFunctionCBSPL::SaveParam(const string& filename){

  // Table param;
  // param.SetHasYErr(false);
  // param.resize(_lam.size(), false);

  // for (int i = 0; i < _lam.size(); i++)
  //   param.set(i, _rbreak(i), _lam(i), 'i');

  // param.Save(filename);

}

void PotentialFunctionCBSPL::extrapolExclParam(){

  if( _extrapol_type == "linear") {

    // extrapolate first _nexcl knot values using linear extrapolation
    // u(r) = ar + b
    // a = m
    // b = - m*r0 + u0
    // m = (u1-u0)/(r1-r0)
    double u0 = _lam(_nexcl);
    double r0 = _rbreak(_nexcl);
    double m = (_lam(_nexcl + 1) - _lam(_nexcl)) /
      (_rbreak(_nexcl + 1) - _rbreak(_nexcl));
    double a = m;
    double b = -1.0*m*r0 + u0;
    for (int i = 0; i < _nexcl; i++)
      _lam(i) = a*_rbreak(i) + b;

  } else if (_extrapol_type=="exponential"){

    // extrapolate first _nexcl knot values using exponential extrapolation
    // u(r) = a * exp( b * r)
    // a = u0 * exp ( - m * r0/u0 )
    // b = m/u0
    // m = (u1-u0)/(r1-r0)
    double u0 = _lam(_nexcl);
    double r0 = _rbreak(_nexcl);
    double m = (_lam(_nexcl + 1) - _lam(_nexcl)) /
      (_rbreak(_nexcl + 1) - _rbreak(_nexcl));
    double a = u0 * exp(-m * r0 / u0);
    double b = m / u0;
    for (int i = 0; i < _nexcl; i++)
      _lam(i)  = a * exp(b * _rbreak(i));

  } else
    throw std::runtime_error("Extrapolation method  \""
                             + _extrapol_type + "\" is not available yet.\n"
                             + "Please specify either \"linear or exponential\" "
                             + " in options file.");
}

void PotentialFunctionCBSPL::setOptParam(const int i, const double val){

  _lam( i + _nexcl ) = val;

}

double PotentialFunctionCBSPL::getOptParam(const int i) const{

  return _lam( i + _nexcl );

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

// calculate first derivative w.r.t. ith parameter
double PotentialFunctionCBSPL::CalculateDF(const int i, const double r) const{

  if ( r >= _min && r <= _cut_off ) {

    int indx = getIndx(r);

    if ( i >= indx && i <= indx+3 ){

      ub::vector<double> R;

      double rk = indx*_dr;

      double t = ( r - rk)/_dr;

      R.resize(4,false); R.clear();

      R(0) = 1.0; R(1) = t; R(2) = t*t; R(3) = t*t*t;

      ub::vector<double> RM = ub::prod(R,_M);

      return RM(i-indx);

    }else
      return 0.0;

  } else
    return 0;

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

// calculate second derivative w.r.t. ith parameter
double PotentialFunctionCBSPL::CalculateD2F(const int i, const int j, const double r) const {

  // for cubic B-SPlines D2F is zero for all lamdas
  return 0.0;

}

double PotentialFunctionCBSPL::CalculateD2F (const double r) const {

  if ( r >= _min && r <= _cut_off) {

    ub::vector<double> R;
    ub::vector<double> B;

    int indx = getIndx(r);

    double rk = indx*_dr;

    double t = ( r - rk)/_dr;

    R.resize(4,false); R.clear();

    R(0) = 0.0; R(1) = 0.0; R(2) = 2.0/(_dr*_dr); R(3) = 6.0*t/(_dr*_dr);

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

double PotentialFunctionCBSPL::CalculateIntRF(const double r1, const double r2) const {

  double r1_,r2_;

  SetInterval(r1,r2,r1_,r2_);

  // interval in which r1_ lies
  int indx1 = getIndx(r1_);
  // interval in which r2_ lies
  int indx2 = getIndx(r2_);

  double uint = 0.0;
  double b,a;

  if (indx1 < indx2) {
    // integral over interval indx1
    a = r1_;
    b = (indx1 + 1) * _dr;

    uint += IntRF(indx1, b) - IntRF(indx1, a);

    // loop over intervals between indx1 and indx2

    for (int indx = indx1 + 1; indx < indx2; indx++) {

      a = indx*_dr;
      b = a + _dr;

      uint += IntRF(indx, b) - IntRF(indx, a);

    }

    // integral over interval indx2
    a = (indx2) * _dr;
    b = r2_;

    uint += IntRF(indx2, b) - IntRF(indx2, a);

  } else if ( indx1 == indx2 ){

    a = r1_;
    b = r2_;

    uint = IntRF(indx1, b) - IntRF(indx1, a);

  }

  return uint;
}

// compute r*u integration at given interval and r
double PotentialFunctionCBSPL::IntRF(const int indx, const double r) const{

  double rk =  indx*_dr;

  double t = (r - rk) / _dr;

  ub::vector<double> R,B;

  R.resize(4, false);
  R.clear();

  for( int m = 0; m < 4; m++){

    //R(m) = _dr*_dr*pow(t,m+2)/(double(m+2))
    //     + rk*_dr*pow(t,m+1)/(double(m+1));
    R(m) = _dr*( rk + m*r + r)*pow(t,m+1)/(m*m + 3.0*m + 2.0);
  }

  ub::vector<double> RM = ub::prod(R, _M);

  B.resize(4, false);
  B.clear();

  B(0) = _lam(indx);
  B(1) = _lam(indx + 1);
  B(2) = _lam(indx + 2);
  B(3) = _lam(indx + 3);

  return ub::inner_prod(B, RM);

}

double PotentialFunctionCBSPL::CalculateIntRdF(const double r1, const double r2) const {

  double r1_,r2_;

  SetInterval(r1,r2,r1_,r2_);

  // interval in which r1_ lies
  int indx1 = getIndx(r1_);
  // interval in which r2_ lies
  int indx2 = getIndx(r2_);

  double uint = 0.0;
  double b,a;

  if (indx1 < indx2) {
    // integral over interval indx1
    a = r1_;
    b = (indx1 + 1) * _dr;

    uint += IntRdF(indx1, b) - IntRdF(indx1, a);

    // loop over intervals between indx1 and indx2

    for (int indx = indx1 + 1; indx < indx2; indx++) {

      a = indx*_dr;
      b = a + _dr;

      uint += IntRdF(indx, b) - IntRdF(indx, a);

    }

    // integral over interval indx2
    a = (indx2) * _dr;
    b = r2_;

    uint += IntRdF(indx2, b) - IntRdF(indx2, a);

  } else if ( indx1 == indx2 ){

    a = r1_;
    b = r2_;

    uint = IntRdF(indx1, b) - IntRdF(indx1, a);

  }

  return uint;
}

// compute r*u integration at given interval and r
double PotentialFunctionCBSPL::IntRdRF(const int indx, const double r) const{

  double rk =  indx*_dr;

  double t = (r - rk) / _dr;

  ub::vector<double> R,B;

  R.resize(4, false);
  R.clear();

  for( int m = 0; m < 4; m++)
    R(m) = ( rk + m*r)*pow(t,m)/(m + 1.0);

  ub::vector<double> RM = ub::prod(R, _M);

  B.resize(4, false);
  B.clear();

  B(0) = _lam(indx);
  B(1) = _lam(indx + 1);
  B(2) = _lam(indx + 2);
  B(3) = _lam(indx + 3);

  return ub::inner_prod(B, RM);

}

double PotentialFunctionCBSPL::CalculateIntR2F(const double r1, const double r2) const {

  double r1_,r2_;

  SetInterval(r1,r2,r1_,r2_);

  // interval in which r1_ lies
  int indx1 = getIndx(r1_);

  // interval in which r2_ lies
  int indx2 = getIndx(r2_);

  double uint = 0.0;
  double b,a;

  if (indx1 < indx2) {
    // integral over interval indx1
    a = r1_;
    b = (indx1 + 1) * _dr;

    uint += IntR2F(indx1, b) - IntR2F(indx1, a);

    // loop over intervals between indx1 and indx2

    for (int indx = indx1 + 1; indx < indx2; indx++) {

      a = indx*_dr;
      b = a + _dr;

      uint += IntR2F(indx, b) - IntR2F(indx, a);

    }

    // integral over interval indx2
    a = (indx2) * _dr;
    b = r2_;

    uint += IntR2F(indx2, b) - IntR2F(indx2, a);

  } else if ( indx1 == indx2 ){

    a = r1_;
    b = r2_;

    uint = IntR2F(indx1, b) - IntR2F(indx1, a);

  }

  return uint;

}

// compute r^2 * u integration at given interval and r
double PotentialFunctionCBSPL::IntR2F(const int indx, const double r) const{

  double rk =  indx*_dr;
  double t = (r - rk) / _dr;

  ub::vector<double> R,B;

  R.resize(4, false);
  R.clear();

  for( int m = 0; m < 4; m++){

    R(m) = _dr * ( _dr*_dr*pow(t,m+3)/(m+3.0) + 2.0*rk*_dr*pow(t,m+2)/(m+2.0) +
                   rk*rk*pow(t,m+1)/(m+1.0) ) ;

  }

  ub::vector<double> RM = ub::prod(R, _M);

  B.resize(4, false);
  B.clear();

  B(0) = _lam(indx);
  B(1) = _lam(indx + 1);
  B(2) = _lam(indx + 2);
  B(3) = _lam(indx + 3);

  return ub::inner_prod(B, RM);

}

double PotentialFunctionCBSPL::CalculateIntRdF(const int i, const double r1,
                                               const double r2) const {

  double r1_, r2_;

  SetInterval(r1, r2, r1_, r2_);

  // interval in which r1_ lies
  int indx1 = getIndx(r1_);

  // interval in which r2_ lies
  int indx2 = getIndx(r2_);

  double uint = 0.0;
  double b, a;

  if (indx1 < indx2) {
    // integral over interval indx1
    a = r1_;
    b = (indx1 + 1) * _dr;

    uint += IntRdF(i, indx1, b) - IntRdF(i, indx1, a);

    // loop over intervals between indx1 and indx2

    for (int indx = indx1 + 1; indx < indx2; indx++) {

      a = indx*_dr;
      b = a + _dr;

      uint += IntRdF(i, indx, b) - IntRdF(i, indx, a);

    }

    // integral over interval indx2
    a = (indx2) * _dr;
    b = r2_;

    uint += IntRdF(i, indx2, b) - IntRdF(i, indx2, a);

  } else if ( indx1 == indx2) {

    a = r1_;
    b = r2_;

    uint = IntRdF(i, indx1, b) - IntRdF(i, indx1, a);

  }

  return uint;

}

double PotentialFunctionCBSPL::IntRdF(const int i, const int indx, const double r)const {

  if ( i + _nexcl >= indx && i + _nexcl <= indx+3 ){

    double rk = indx*_dr;

    double t = (r - rk) / _dr;

    ub::vector<double> R;

    R.resize(4, false);
    R.clear();

    for (int m = 0; m < 4; m++) {

      //R(m) = _dr * _dr * pow(t, m + 2) / (double(m + 2))
      //        + rk * _dr * pow(t, m + 1) / (double(m + 1));
      R(m) = _dr*( rk + m*r + r)*pow(t,m+1)/(m*m + 3.0*m + 2.0);
    }

    ub::vector<double> RM = ub::prod(R, _M);

    return RM(i + _nexcl - indx);

  } else {

    return 0.0;

  }

}

double PotentialFunctionCBSPL::CalculateIntR2dF(const int i, const double r1,
                                                const double r2) const {

  double r1_, r2_;

  SetInterval(r1, r2, r1_, r2_);

  // interval in which r1_ lies
  int indx1 = getIndx(r1_);

  // interval in which r2_ lies
  int indx2 = getIndx(r2_);

  double uint = 0.0;
  double b, a;

  if (indx1 < indx2) {
    // integral over interval indx1
    a = r1_;
    b = (indx1 + 1) * _dr;

    uint += IntR2dF(i, indx1, b) - IntR2dF(i, indx1, a);

    // loop over intervals between indx1 and indx2

    for (int indx = indx1 + 1; indx < indx2; indx++) {

      a = indx*_dr;
      b = a + _dr;

      uint += IntR2dF(i, indx, b) - IntR2dF(i, indx, a);

    }

    // integral over interval indx2
    a = (indx2) * _dr;
    b = r2_;

    uint += IntR2dF(i, indx2, b) - IntR2dF(i, indx2, a);

  } else if ( indx1 == indx2) {

    a = r1_;
    b = r2_;

    uint = IntR2dF(i, indx1, b) - IntR2dF(i, indx1, a);

  }

  return uint;

}

double PotentialFunctionCBSPL::IntR2dF(const int i, const int indx, const double r)const {

  if ( i + _nexcl >= indx && i + _nexcl <= indx+3 ){

    double rk = indx*_dr;

    double t = (r - rk) / _dr;

    ub::vector<double> R;

    R.resize(4, false);
    R.clear();

    for (int m = 0; m < 4; m++) {

      R(m) = _dr * ( _dr*_dr*pow(t,m+3)/(m+3.0) + 2.0*rk*_dr*pow(t,m+2)/(m+2.0) +
                     rk*rk*pow(t,m+1)/(m+1.0) ) ;
    }

    ub::vector<double> RM = ub::prod(R, _M);

    return RM(i + _nexcl - indx);

  } else {

    return 0.0;

  }

}

void PotentialFunctionCBSPL::CutoffContinuous(const double uval){

  int indx = _lam.size() - 1 ;

  for( int i = 0; i < _ncut; i++){
    _lam(indx-i) = uval;
  }

  extrapolExclParam();

}
