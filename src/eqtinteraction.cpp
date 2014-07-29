/*
* EQTInteraction class
*
* by Sikandar Mashayak , Wednesday, Dec., 14, 2011
*/

#include "eqtinteraction.hpp"

EQTInteraction::EQTInteraction(int id){

  _id = id;
  _type1;
  _type2;
  _name = _type1 + "-" + _type2;

  string filename;
  double rmin,rcut;
  double rexcl = 0;
  int nlam = 0;
  string extrapol_type="";

  rmin;
  rcut;

  nlam;
  rexcl;
  extrapol_type;

  filename;
  _uatm;

  _ucont;

  _rcut;
  _rmin;

}

void EQTInteraction::SavePotTab(const string& filename,
                                const double step,
                                const double rmin, const double rcut) const {

  int ngrid = (int) ((rcut - rmin) / step + 1.00000001);

  double r_init;
  int i;
  char flag;

}

double EQTInteraction::ComputeU(double r)const{

  return _uatm->CalculateF(r) + _ucont->CalculateF(r);

}

double EQTInteraction::ComputeF(double r)const{

  return -1.0*(_uatm->CalculateDF(r) + _ucont->CalculateDF(r));

}

double EQTInteraction::ComputeIntRU(double r1, double r2){

  return _uatm->CalculateIntRF(r1,r2) + _ucont->CalculateIntRF(r1,r2);

}

double EQTInteraction::ComputeAtmIntRU(double r1, double r2){

  return _uatm->CalculateIntRF(r1,r2);

}

double EQTInteraction::ComputeContIntRU(double r1, double r2){

  return _ucont->CalculateIntRF(r1,r2);

}
