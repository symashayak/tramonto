/*
 *
 * Author: sikandar
 *
 * Created on November 8, 2011, 11:52 PM
 */

#ifndef POTENTIALFUNCTION_H
#define	POTENTIALFUNCTION_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/lexical_cast.hpp>
#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>

using namespace std;
namespace ub = boost::numeric::ublas;

class PotentialFunction {
 public:

  virtual ~PotentialFunction() {}

  // return type of function
  string Type() const { return _type;}

  // read parameters from the input file
  virtual void setParam(string filename);

  // save parameters to the file
  virtual void SaveParam(const string& filename);

  // write potential table for specified interval
  virtual void SavePotTab(const string& filename, const double step,
                          const double rmin, const double rcut);
  // set all parameters
  void setParam(const ub::vector<double> param){ _lam = param; }

  // set ith parameter
  void setParam(const int i, const double val) { _lam(i) = val; }

  // set ith parameter among those to be optimized
  virtual void setOptParam(const int i, const double val) {
    setParam(i,val);
  }

  // set minimum r value to avoid large values
  void setMinDist(const double min) { _min = min; }

  // set cut-off value
  void setCutOffDist(const double cutoff) { _cut_off = cutoff; }

  // calculate function
  virtual double CalculateF (const double r) const = 0;

  // return parameter
  ub::vector<double>& Params() { return _lam; }

  // return ith parameter
  double getParam(const int i) const { return _lam(i); }

  // return size of parameters
  int getParamSize() const { return _lam.size(); }

  // return cut-off value
  double getCutOff() const { return _cut_off; }

  double getMinDist() const { return _min; }

  // first derivative w.r.t. r
  virtual double CalculateDF( const double r) const = 0;

 protected:

  PotentialFunction(const string& type,
                    const int nlam_,const double min_,const double max_);

  ub::vector<double> _lam;
  double _cut_off;
  double _min;
  string _type;

  virtual void SetInterval(const double r1, const double r2,
                           double& r1_, double& r2_) const;
};

#endif
