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
  votca::tools::Table param;
  param.Load(filename);
  int nknots = param.size();

  _nbreak = nknots -2;
  _lam.resize(nknots,false);
  _lam.clear();
  _rbreak.resize(nknots,false);
  _rbreak.clear();

  _min = param.x(0);
  for( int i = 0; i < _lam.size(); i++){
    _rbreak(i) = param.x(i);
    _lam(i) = param.y(i);
  }
  // cut-off is the last break point, i.e., nknots-2 th point
  _cut_off = param.x(nknots-3);
  _dr = _rbreak(1) - _rbreak(0);

  SavePotTab("test.dat",0.01,_min,_cut_off);
}

void PotentialFunctionCBSPL::SaveParam(const string& filename){

  votca::tools::Table param;
  param.SetHasYErr(false);
  param.resize(_lam.size(), false);

  for (int i = 0; i < _lam.size(); i++)
    param.set(i, _rbreak(i), _lam(i), 'i');

  param.Save(filename);

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

double PotentialFunctionCBSPL::CalculateIntR2F (const double r) const {

  /* In Tramonto, we need to provide intgration stencil at given r.
   * However, B-splines are piecewise and therefore, we must integrate them piecewise.
   * For piecewise integration, both the upper and lower r limits should be given.
   * Therefore, here, piecewise integration value from min to r is returned.
   * This approach still gives accurate values for integration from r1 to r2, because
   * int(r2,r1) = int(r2,_min)-int(r1,_min)
   */
  cout << _min << "\t" << r << endl;
  double r1_ = _min;
  double r2_ = r;
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
  double rk = indx*_dr;
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
