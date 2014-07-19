/*
//@HEADER
// ******************************************************************** 
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1
// of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

/*
 *  FILE: dft_pairPot_BSpline.c
 *
 *  This file contains routines specific to a strict mean field implementation of
 *  a cubic B-splines based potential.
// \todo implement the split of the potential into the infinit hard core.
 * with a Weeks-Chandler-Anderson split of the LJ potential
 *  into an infinite hard core and an attractive piece.  Note that the hard core
 *  diameters can be modified by the Barker-Henderson approach if desired.
 */

#include "dft_pairPot_BSpline.h"

/******************************************************************************/
/* uBSpline: The cubic B-Splines based potential                         */

double uBSpline(double r,double rcut)
{
  double u;

  if (r <= rcut) {
    // \todo implement u computation
  }
  else u = 0.0;
  return u;
}
/*******************************************************************************/
/* uBSpline_CS_setparams: The parameters for the cut and shifted BSpline potential */
void uBSpline_setparams(int context, int i, int j, double *param1)
{
  switch (context){
     case FLUID_FLUID:
        *param1 = Cut_ff[i][j];
        break;
     case WALL_FLUID:
        *param1 = Cut_wf[i][WallType[j]];
        break;
     case WALL_WALL:
        *param1 = Cut_ww[WallType[i]][WallType[j]];
        break;
     default:
        if (Iwrite_screen != NONE) printf("problem with potential context uBSpline_setparams\n");
        exit(-1);
   }
   return;
}
/*******************************************************************************/
/* uBSpline_DERIV1D: The derivative of a cubic B-splines based potential in the x (or y or z) direction */

double uBSpline_DERIV1D(double r, double x, double rcut)
{
  double uderiv;

  if (r <= rcut) {
    // \todo implement derivative of u computation
  }
  else uderiv = 0.0;
  return uderiv;
}
/******************************************************************************/
/* uBSpline_InnerCore : define the properties of the inner core of the potential based on
                  input parameters */
void uBSpline_InnerCore(int i, int j,double *rCore_left, double *rCore_right)
{
  // \todo implement the InnerCore values
   switch(Type_CoreATT_R){
      case ATTCORE_SIGMA:
        //*rCore_right=
        //*rCore_left=0.0;
        break;
      case ATTCORE_UMIN:
        //*rCore_right=
        //*rCore_left=0.0;
        break;
      case ATTCORE_UCSZERO:
        //*rCore_right
        //*rCore_left=0.0;
        break;
      case ATTCORE_SIGTOUMIN:
        //*rCore_right
        //*rCore_left=Sigma_ff[i][j];
        break;
      default:
        if (Iwrite_screen != NONE) printf("Problem with Type_CoreATT_R - set to %d\n",Type_CoreATT_R);
        exit(-1);
   }
   switch(Type_CoreATT_CONST){
      case CORECONST_UCONST:
        //*epsCore=uBSpline_ATT_noCS(*rCore_right,i,j);
        break;
      case CORECONST_ZERO:
        //*epsCore=0.0;
        break;
      default:
        if (Iwrite_screen != NONE) printf("Problem with Type_CoreATT_CONST - set to %d\n",Type_CoreATT_CONST);
        exit(-1);
   }
   return;
}
/******************************************************************************/
/* uBSpline_ATT_CS: the pair potential (based on a 12-6 LJ fluid) that will be used 
                  as the attractive perturbation to a hard sphere reference fluid
                  in strict mean field DFT calculations */
double uBSpline_ATT_CS(double r,int i, int j)
{
  double uatt,r_min,sigma2,sigma6;
  double r_inv,r2_inv,r6_inv,r12_inv;
  double rc_inv,rc2_inv,rc6_inv,rc12_inv;

  sigma2 = Sigma_ff[i][j]*Sigma_ff[i][j];
  sigma6 = sigma2*sigma2*sigma2;

  rc_inv   = 1.0/Cut_ff[i][j];
  rc2_inv  = rc_inv*rc_inv;
  rc6_inv  = rc2_inv*rc2_inv*rc2_inv;
  rc12_inv = rc6_inv*rc6_inv;

  switch(Type_CoreATT_R){
     case ATTCORE_SIGMA:      r_min=Sigma_ff[i][j]; break;
     case ATTCORE_SIGTOUMIN:  
     case ATTCORE_UMIN:       r_min=Rmin_ff[i][j]; break; /* should be Sigma_ff[i][j]*pow(2.0,1.0/6.0) */
     case ATTCORE_UCSZERO:    r_min=Rzero_ff[i][j]; break;
  }


  if ((r<r_min && Type_CoreATT_CONST==CORECONST_ZERO) ||
      (r<Sigma_ff[i][j] && Type_CoreATT_R==ATTCORE_SIGTOUMIN)) uatt=0.0;  
  else{
     if (r <= Cut_ff[i][j]) {
 
        if (r < r_min) r = r_min; 

        r_inv = 1.0/r;
        r2_inv  = r_inv*r_inv;
        r6_inv  = r2_inv*r2_inv*r2_inv;
        r12_inv = r6_inv*r6_inv;

        uatt = 4.0 * Eps_ff[i][j]* sigma6 * (
               sigma6*(r12_inv - rc12_inv)
                    - (r6_inv  - rc6_inv ) );
     }
     else uatt = 0.0;

  }
  return uatt;
}
/****************************************************************************/
/* uBSpline_ATT_noCS:  calculate the attractive part of
                  12-6 LJ potential at the minimum. */

double uBSpline_ATT_noCS(double r,int i, int j)
{
  double uatt,r_min,sigma2,sigma6;
  double r_inv,r2_inv,r6_inv,r12_inv;

  sigma2 = Sigma_ff[i][j]*Sigma_ff[i][j];
  sigma6 = sigma2*sigma2*sigma2;

  switch(Type_CoreATT_R){
     case ATTCORE_SIGMA:      r_min=Sigma_ff[i][j]; break;
     case ATTCORE_SIGTOUMIN:
     case ATTCORE_UMIN:       r_min=Rmin_ff[i][j]; break; /* should be Sigma_ff[i][j]*pow(2.0,1.0/6.0) */
     case ATTCORE_UCSZERO:    r_min=Rzero_ff[i][j]; break;
  }

  if ((r < r_min && Type_CoreATT_CONST==CORECONST_ZERO) ||
     (r<Sigma_ff[i][j] && Type_CoreATT_R==ATTCORE_SIGTOUMIN))  uatt=0.0;
  else{

      if (r < r_min) r = r_min; 

      r_inv = 1.0/r;

      r2_inv  = r_inv*r_inv;
      r6_inv  = r2_inv*r2_inv*r2_inv;
      r12_inv = r6_inv*r6_inv;

      uatt = 4.0 * Eps_ff[i][j]* sigma6 * ( sigma6*r12_inv  - r6_inv);
      }

  return uatt;
}
/****************************************************************************/
/* uBSpline_IntStencil:  the integral of the 12-6 potential that is used
                        to define the magnitude of the DFTMFT UATTRACT
                        integration stencil. */

double uBSpline_Integral(double r,int i, int j)
{
  double uatt_int, sigma2,sigma6;
  double r_inv,r3_inv,r9_inv;

  sigma2 = Sigma_ff[i][j]*Sigma_ff[i][j];
  sigma6 = sigma2*sigma2*sigma2;

  r_inv = 1.0/r;

  r3_inv  = r_inv*r_inv*r_inv;
  r9_inv  = r3_inv*r3_inv*r3_inv;

  uatt_int = 16 * PI * Eps_ff[i][j]* sigma6 * ( - sigma6*r9_inv/9.0  + r3_inv/3.0 );

  return uatt_int;
}
/****************************************************************************/

