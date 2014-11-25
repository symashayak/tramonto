/*
* dft_pairPot_eqt class
*
* by Sikandar Mashayak , Tuesday, July, 29, 2014
*/

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/lexical_cast.hpp>
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

ub::matrix<EQTInteraction*> u_ff(NCOMP_MAX,NCOMP_MAX);
ub::matrix<EQTInteraction*> u_wf(NCOMP_MAX,NWALL_MAX_TYPE);
ub::matrix<EQTInteraction*> u_ww(NWALL_MAX_TYPE,NWALL_MAX_TYPE);

//c_dft_pairPot_EQT pairPot_EQT;

/***************************************************************************
 * uEQT_setparams : set parameters for i,j interation of type context
 *                  and set respective params
 ***************************************************************************/
void uEQT_setparams(int context,int i,int j,double *param1,
                    double *param2,double *param3){
  switch (context){
  case FLUID_FLUID:
    *param1 = 0;
    *param2 = i;
    *param3 = j;
    u_ff(i,j) = new EQTInteraction("ff_"+
                                   boost::lexical_cast<std::string>(i),
                                   boost::lexical_cast<std::string>(j));
    break;
  case WALL_FLUID:
    *param1 = 1;
    *param2 = i;
    *param3 = j;
    u_wf(i,j) = new EQTInteraction("wf_"+
                                   boost::lexical_cast<std::string>(i),
                                   boost::lexical_cast<std::string>(j));
    break;
  case WALL_WALL:
    *param1 = 2;
    *param2 = i;
    *param3 = j;
    u_ww(i,j) = new EQTInteraction("ww_"+
                                   boost::lexical_cast<std::string>(i),
                                   boost::lexical_cast<std::string>(j));
    break;
  default:
    if (Iwrite_screen != NONE) printf("problem with potential context uEQT_setparams\n");
    exit(-1);
  }
  return;
}

/***************************************************************************
 * uEQT : compute potential of a given i,j pair of type flag at r
 ***************************************************************************/
double uEQT(double r,int flag,int i,int j){
  if(flag == 0){ // fluid-fluid
    return u_ff(i,j)->ComputeU(r);
  } else if (flag == 1) { // wall-fluid
    return u_wf(i,j)->ComputeU(r);
  } else if (flag == 2) { // wall-wall
    return u_ww(i,j)->ComputeU(r);
  }
}

/***************************************************************************
 * uEQT_DERIV1D : compute first derivative of potential w.r.t. r
 ***************************************************************************/
double uEQT_DERIV1D(double r,double x,int flag,int i,int j){
  if(flag == 0){ // fluid-fluid
    return (x/r)*u_ff(i,j)->ComputeDU(r);
  } else if (flag == 1) { // wall-fluid
    return (x/r)*u_wf(i,j)->ComputeDU(r);
  } else if (flag == 2) { // wall-wall
    return (x/r)*u_ww(i,j)->ComputeDU(r);
  }
}

/***************************************************************************
 * uEQT_InnerCore : define the properties of the inner core of the potential
 * based on input parameters
 ***************************************************************************/
void uEQT_InnerCore(int i,int j,double *rCore_left,
                    double *rCore_right,double *epsCore){
   switch(Type_CoreATT_R){
   case ATTCORE_SIGMA:
     *rCore_right=Sigma_ff[i][j];
     *rCore_left=0.0;
     break;
   case ATTCORE_UMIN:
     *rCore_right=Rmin_ff[i][j];
     *rCore_left=0.0;
     break;
   case ATTCORE_UCSZERO:
     *rCore_right=Rzero_ff[i][j];
     *rCore_left=0.0;
     break;
   case ATTCORE_SIGTOUMIN:
     *rCore_right=Rmin_ff[i][j];
     *rCore_left=Sigma_ff[i][j];
     break;
   default:
     if (Iwrite_screen != NONE)
       printf("Problem with Type_CoreATT_R - set to %d\n",Type_CoreATT_R);
     exit(-1);
   }
   switch(Type_CoreATT_CONST){
   case CORECONST_UCONST:
     *epsCore=uEQT_ATT_noCS(*rCore_right,i,j);
     break;
   case CORECONST_ZERO:
     *epsCore=0.0;
     break;
   default:
     if (Iwrite_screen != NONE)
       printf("Problem with Type_CoreATT_CONST - set to %d\n",Type_CoreATT_CONST);
     exit(-1);
   }
   return;
}

/***************************************************************************
 * uEQT_ATT_CS : (cutoff-shifted)
 *               the pair potential (based on a 12-6 LJ fluid) that will be
 *               used as the attractive perturbation to a hard sphere
 *               reference fluid in strict mean field DFT calculations
 ***************************************************************************/
double uEQT_ATT_CS(double r,int i,int j) {
  double uatt,r_min;

  switch(Type_CoreATT_R){
  case ATTCORE_SIGMA:
    r_min=Sigma_ff[i][j];
    break;
  case ATTCORE_SIGTOUMIN:
  case ATTCORE_UMIN:
    r_min=Rmin_ff[i][j];
    break;
  case ATTCORE_UCSZERO:
    r_min=Rzero_ff[i][j];
    break;
  }

  if ((r<r_min && Type_CoreATT_CONST==CORECONST_ZERO) ||
      (r<Sigma_ff[i][j] && Type_CoreATT_R==ATTCORE_SIGTOUMIN))
    uatt=0.0;
  else{
     if (r <= Cut_ff[i][j]) {

        if (r < r_min) r = r_min;

        uatt = u_ff(i,j)->ComputeU(r);
     }
     else uatt = 0.0;

  }
  return uatt;
}

/***************************************************************************
 * uEQT_ATT_noCS : (NO cutoff-shifted)
 *               the pair potential (based on a 12-6 LJ fluid) that will be
 *               used as the attractive perturbation to a hard sphere
 *               reference fluid in strict mean field DFT calculations
 ***************************************************************************/
double uEQT_ATT_noCS(double r,int i,int j){
  double uatt,r_min;

  switch(Type_CoreATT_R){
  case ATTCORE_SIGMA:
    r_min=Sigma_ff[i][j];
    break;
  case ATTCORE_SIGTOUMIN:
  case ATTCORE_UMIN:
    r_min=Rmin_ff[i][j];
    break;
  case ATTCORE_UCSZERO:
    r_min=Rzero_ff[i][j];
    break;
  }

  if ((r<r_min && Type_CoreATT_CONST==CORECONST_ZERO) ||
      (r<Sigma_ff[i][j] && Type_CoreATT_R==ATTCORE_SIGTOUMIN))
    uatt=0.0;
  else{
     if (r <= Cut_ff[i][j]) {

        if (r < r_min) r = r_min;

        uatt = u_ff(i,j)->ComputeU(r);
     }
     else uatt = 0.0;
  }
  return uatt;
}

/***************************************************************************
 * uEQT_Integral : the integral of the EQT potential that is used to define
 *                  the magnitude of the DFTMFT UATTRACT integration stencil
 ***************************************************************************/
double uEQT_Integral(double r,int i,int j){
  double uatt_int;
  uatt_int = 4.0 * PI * u_ff(i,j)->ComputeIntR2U(r);
  return uatt_int;
}






