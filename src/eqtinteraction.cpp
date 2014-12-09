/*
* EQTInteraction class
*
* by Sikandar Mashayak , Wednesday, Dec., 14, 2011
*/

#include "eqtinteraction.hpp"

EQTInteraction::EQTInteraction(const string& type1,
                               const string& type2){
  _type1 = type1;
  _type2 = type2;
  _name = _type1 + "_" + _type2;

  string filename;

  filename = _name + "_quad.dat";
  _uquad = new PotentialFunctionQUAD();
  _uquad->setParam(filename);

  filename = _name + "_spl.dat";
  _uspl = new PotentialFunctionCBSPL();
  _uspl->setParam(filename);

  // filename = _name + "_spl.dat";
  // _uspl = new PotentialFunctionLJ();
  // _uspl->setParam(filename);
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
  return _uquad->CalculateF(r) + _uspl->CalculateF(r);
  // return _uspl->CalculateF(r);
}

double EQTInteraction::ComputeDU(double r)const{
  return (_uquad->CalculateDF(r) + _uspl->CalculateDF(r));
  // return _uspl->CalculateDF(r);
}

double EQTInteraction::ComputeIntR2U(double r)const{
  return (_uquad->CalculateIntR2F(r) + _uspl->CalculateIntR2F(r));
  // return _uspl->CalculateIntR2F(r);
}
