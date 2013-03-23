#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>

#include "EOS.h"
using namespace std;

EOS::EOS()
{
   EOS_epsT = new Table("EOS_tables/EOS_PST.dat"); 
   EOS_mu = new Table("EOS_tables/EOS_Mu.dat"); 
   EOS_coeffs = new Table("EOS_tables/coeff.dat"); 
}

EOS::~EOS()
{
   delete EOS_epsT;
   delete EOS_mu;
   delete EOS_coeffs;
}

double EOS::getOneFromtheOther(int xColnum, int yColnum, double xVal)
{
   double yVal;
   double xMin = EOS_epsT->getFirst(xColnum);
   double xMax = EOS_epsT->getLast(xColnum);
   double yMin = EOS_epsT->getFirst(yColnum);
  
   double extrapCoeff1 = EOS_coeffs->get(1, yColnum-1);
   double extrapCoeff2 = EOS_coeffs->get(2, yColnum-1);
   double e;
   double eMin = EOS_epsT->getFirst(1);
   if(xColnum == 1)
      e = xVal;
   else
      e = EOS_epsT->interp(1, xColnum, xVal, 2);

   string function_name = "EOS::getOneFromtheOther";
   if(xVal < 0)
   {
      int errLevel = 3;
      outputFunctionerror(function_name, "xVal < 0!", xVal, errLevel);
   }
   else if(xVal < xMin)
   {
      int errLevel = 1;
      outputFunctionerror(function_name, "xVal < xMin, using linear extrapolation instead.", xVal, errLevel);
      yVal = e*yMin/eMin;
   }
   else if(xVal > xMax)
   {
      int errLevel = 1;
      outputFunctionerror(function_name, "xVal > xMax, using polynomial extrapolation instead.", xVal, errLevel);
      yVal = extrapCoeff1*pow(e, extrapCoeff2);
   }
   else
      yVal = EOS_epsT->interp(xColnum, yColnum, xVal, 2);
   return(yVal);
}
