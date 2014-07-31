/*
* dft_pairPot_eqt class
*
* by Sikandar Mashayak , Tuesday, July, 29, 2014
*/

#include <boost/numeric/ublas/matrix.hpp>
#include "eqtinteraction.hpp"
#include "dft_pairPot_EQT.h"

namespace ub = boost::numeric::ublas;
using namespace std;

// class c_dft_pairPot_EQT {

// public:

//   c_dft_pairPot_EQT():
//     u_ff(NCOMP_MAX,NCOMP_MAX),
//     u_wf(NCOMP_MAX,NWALL_MAX_TYPE),
//     u_ww(NWALL_MAX_TYPE,NWALL_MAX_TYPE)
//   {}
//   ~c_dft_pairPot_EQT() {}

// protected:

//   ub::matrix<EQTInteraction> u_ff;
//   ub::matrix<EQTInteraction> u_wf;
//   ub::matrix<EQTInteraction> u_ww;

// };

ub::matrix<EQTInteraction> u_ff(NCOMP_MAX,NCOMP_MAX);
ub::matrix<EQTInteraction> u_wfu_wf(NCOMP_MAX,NWALL_MAX_TYPE);
ub::matrix<EQTInteraction> u_ww(NWALL_MAX_TYPE,NWALL_MAX_TYPE);

//c_dft_pairPot_EQT pairPot_EQT;

void uEQT_setparams(int context,int i,int j,double *param1,
                    double *param2,double *param3){
  switch (context){
  case FLUID_FLUID:
    *param1 = 0;
    *param2 = i;
    *param3 = j;
    //uff(i,j).setparams();
    break;
  case WALL_FLUID:
    *param1 = 1;
    *param2 = i;
    *param3 = j;
    //uwf(i,j).setparams();
    break;
  case WALL_WALL:
    *param1 = 2;
    *param2 = i;
    *param3 = j;
    //uww(i,j).setparams();
    break;
  default:
    if (Iwrite_screen != NONE) printf("problem with potential context uEQT_setparams\n");
    exit(-1);
  }
  return;
}

double uEQT(double r,int flag,int i,int j){
  if(flag == 0){ // fluid-fluid
    // return uff(i,j).calculateU(r);
  } else if (flag == 1) { // wall-fluid
    // return uwf(i,j).calculateU(r);
  } else if (flag == 2) { // wall-wall
    // return uww(i,j).calculateU(r);
  }
}

double uEQT_DERIV1D(double r,double x,int flag,int i,int j){
  if(flag == 0){ // fluid-fluid
    // return x*uff(i,j).calculateDU(r);
  } else if (flag == 1) { // wall-fluid
    // return x*uwf(i,j).calculateDU(r);
  } else if (flag == 2) { // wall-wall
    // return x*uww(i,j).calculateDU(r);
  }
}

double uEQT_Integral(double r,int i,int j){
}

double uEQT_ATT_CS(double r,int i,int j) {
}


double uEQT_ATT_noCS(double r,int i,int j){
}

void uEQT_InnerCore(int i,int j,double *rCore_left,
                    double *rCore_right,double *epsCore){
}
