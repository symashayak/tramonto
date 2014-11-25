/* This file was automatically generated.  Do not edit! */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(HAS_VALUES_H)
#include <values.h>
#include <unistd.h>
#include <string.h>
#endif
#include "mpi.h"
#include "az_aztec.h"
#include "az_aztec_defs.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
#define PI    3.141592653589793238462643383279502884197169399375
#define CORECONST_ZERO      1
#define CORECONST_UCONST    0
#define ATTCORE_SIGTOUMIN   3
#define NCOMP_MAX 5
#define ATTCORE_UCSZERO     2
#define ATTCORE_UMIN        1
#define ATTCORE_SIGMA       0
#define NONE       -1
#define NWALL_MAX_TYPE 20
#define WALL_WALL   2
#define NWALL_MAX 600
#define WALL_FLUID  1
#define FLUID_FLUID 0
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
extern int Type_CoreATT_CONST;
extern double Rzero_ff[NCOMP_MAX][NCOMP_MAX];
extern double Rmin_ff[NCOMP_MAX][NCOMP_MAX];
extern int Type_CoreATT_R;
extern int Iwrite_screen;
extern int WallType[NWALL_MAX];
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern double Cut_ff[NCOMP_MAX][NCOMP_MAX];

#ifdef __cplusplus
extern "C" {
#endif
  void uEQT_setparams(int context,int i,int j,double *param1,double *param2,double *param3);
  double uEQT(double r,int flag,int i,int j);
  double uEQT_DERIV1D(double r,double x,int flag,int i,int j);
  double uEQT_Integral(double r,int i,int j);
  double uEQT_ATT_CS(double r,int i,int j);
  double uEQT_ATT_noCS(double r,int i,int j);
  void uEQT_InnerCore(int i,int j,double *rCore_left,double *rCore_right,double *epsCore);
#ifdef __cplusplus
}
#endif
