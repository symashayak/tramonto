#include "potentialfunction.hpp"

PotentialFunction::PotentialFunction(const string& type, const int nlam_,
                                     const double min_,const double max_){
  _type = type;
  _lam.resize(nlam_);
  _lam.clear();
  _min = min_;
  _cut_off = max_;

}


void PotentialFunction::setParam(string filename) {

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

}

void PotentialFunction::SaveParam(const string& filename){

  // Table param;
  // param.SetHasYErr(false);
  // param.resize(_lam.size(), false);

  // for (int i = 0; i < _lam.size(); i++){
  //     param.set(i, i, _lam(i), 'i');
  // }
  // param.Save(filename);

}
void PotentialFunction::SavePotTab(const string& filename,
                                   const double step,
                                   const double rmin, const double rcut){

  int ngrid = (int) ((rcut - rmin) / step + 1.00000001);

  votca::tools::Table pot_tab;
  pot_tab.SetHasYErr(false);
  pot_tab.resize(ngrid, false);

  double r_init;
  int i;
  char flag;
  for (r_init = rmin, i = 0; i < ngrid - 1; r_init += step) {

    // put point, result, flag at point out_x into the table

    if( r_init < _min || r_init > _cut_off )
      flag = 'o';
    else
      flag = 'i';

    pot_tab.set(i++, r_init, CalculateF(r_init), flag);

  }

  if( rcut < _min || rcut > _cut_off )
    flag = 'o';
  else
    flag = 'i';

  pot_tab.set(i, rcut, CalculateF(rcut), flag);

  // save table in the file
  pot_tab.Save(filename);

}

void PotentialFunction::SetInterval(const double r1, const double r2,
                                    double& r1_, double& r2_) const {

  if (r2 < r1) {
    throw std::runtime_error("integration interval error\n");
  }


  // check if r2 falls out of _min and _cut_off range
  if (r2 > _cut_off) {

    r2_ = _cut_off;

  } else if (r2 < _min) {

    r2_ = _min;

  } else {

    r2_ = r2;

  }

  if (r1 < _min) {

    r1_ = _min;

  } else if (r1 > _cut_off) {

    r1_ = _cut_off;

  } else {

    r1_ = r1;

  }


}

