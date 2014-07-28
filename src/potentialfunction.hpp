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

  // calculate first derivative w.r.t. ith parameter
  virtual double CalculateDF(const int i, const double r) const = 0;

  // calculate second derivative w.r.t. ith parameter
  virtual double CalculateD2F(const int i, const int j, const double r) const = 0;

  // return parameter
  ub::vector<double>& Params() { return _lam; }

  // return ith parameter
  double getParam(const int i) const { return _lam(i); }

  // return ith parameter among those to be optimized
  virtual double getOptParam(const int i) const {
    return getParam(i);
  }

  // return size of parameters
  int getParamSize() const { return _lam.size(); }

  // return number of parameters to be optimized for continuum potential
  // this is different than _nlam since some paramters are fixed by
  // continuity conditions at rmin
  virtual int ContOptParamSize() const { return _lam.size();}

  // return cut-off value
  double getCutOff() const { return _cut_off; }

  double getMinDist() const { return _min; }

  // EQT specific functions

  // calculate integrant r.f(r)
  virtual double CalculateIntRF(const double r1, const double r2) const = 0;

  /// \brief calculate integration r^2 * f(r) dr from r1 to r2
  virtual double CalculateIntR2F(const double r1, const double r2) const = 0;

  /// \brief caculates integration r*(df/dr) from r1 to r2
  virtual double CalculateIntRdF(const double r1, const double r2) const = 0;

  // calculate integrant r.df/dlamdai(r)
  virtual double CalculateIntRdF(const int i, const double r1,
                                 const double r2) const = 0;

  /// \brief calculate r^2.df/dlamdai(r)
  virtual double CalculateIntR2dF(const int i, const double r1,
                                  const double r2) const = 0;

  // calculate integrant r.d2f/dlamdai*dlamdaj(r)
  virtual double CalculateIntRd2F(const int i, const int j,
                                  const double r1, const double r2) const = 0;

  /// \brief calculate r^2.d2f/dlamdai*dlamdaj(r)
  virtual double CalculateIntR2d2F(const int i, const int j,
                                   const double r1, const double r2) const = 0;

  // makes this potential and u continuous at r
  virtual void CutoffContinuous(const double uval) = 0;

  // first derivative w.r.t. r
  virtual double CalculateDF( const double r) const = 0;

  // double derivative w.r.t. r
  virtual double CalculateD2F( const double r) const = 0;

 protected:

  PotentialFunction(const string& type,
                    const int nlam_,const double min_,const double max_);

  ub::vector<double> _lam;
  double _cut_off;
  double _min;
  string _type;

  virtual void SetInterval(const double r1, const double r2,
                           double& r1_, double& r2_) const;

  virtual double IntRF(const double r) const = 0;
  virtual double IntR2F(const double r) const = 0;
  virtual double IntRdF(const double r) const = 0;
  virtual double IntRdF(const int i, const double r) const = 0;
  virtual double IntR2dF(const int i, const double r) const = 0;
  virtual double IntRd2F(const int i, const int j, const double r) const = 0;
  virtual double IntR2d2F(const int i, const int j, const double r) const = 0;

};

#endif
