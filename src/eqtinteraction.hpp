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

  EQTInteraction(int id);
  ~EQTInteraction() {}

  int getId() const { return _id; }
  void setId(int id) { _id = id; }
  const string &getName() const { return _name; }
  void setName(const string &name) {  _name=name; }
  const string &getType1() const { return _type1; }
  const string &getType2() const { return _type2; }
  const double &Rcut() const { return _rcut; }
  const double &Rmin() const { return _rmin; }
  void setType1(const string &type1) { _type1 = type1; }
  void setType2(const string &type2) { _type2 = type2; }
  double ComputeU(double r)const;
  double ComputeF(double r)const;
  double ComputeIntRU(double r1, double r2);
  double ComputeAtmIntRU(double r1, double r2);
  double ComputeContIntRU(double r1, double r2);
  void SavePotTab(const string& filename, const double step,
                  const double rmin, const double rcut) const;
private:

  int _id;
  string _name;
  string _type1,_type2;

  double _rcut; // larger among uatm and ucont
  double _rmin; // smaller among uatm and ucont

  PotentialFunctionQUAD* _ucont;
  PotentialFunctionCBSPL* _uatm;

};

#endif // EQTINTERACTION_H
