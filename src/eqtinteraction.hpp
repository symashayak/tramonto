/*
* EQTInteraction class
*
* by Sikandar Mashayak , Wednesday, Dec., 14, 2011
*/

#ifndef EQTINTERACTION_H
#define	EQTINTERACTION_H

#include <string>
#include "potentialfunctioncbspl.hpp"
#include "potentialfunctionquad.hpp"

using namespace std;

class EQTInteraction {

public:

  EQTInteraction() {}
  EQTInteraction(const string& type1,const string& type2);
  ~EQTInteraction() {}
  const string &getName() const { return _name; }
  void setName(const string &name) {  _name=name; }
  const string &getType1() const { return _type1; }
  const string &getType2() const { return _type2; }
  void setType1(const string &type1) { _type1 = type1; }
  void setType2(const string &type2) { _type2 = type2; }
  double ComputeU(double r)const;
  double ComputeDU(double r)const;
  void SavePotTab(const string& filename, const double step,
                  const double rmin, const double rcut) const;
private:
  string _name;
  string _type1,_type2;
  PotentialFunctionQUAD* _uquad;
  PotentialFunctionCBSPL* _uspl;
};

#endif // EQTINTERACTION_H
