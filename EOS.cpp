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

double EOS::getEntropydensityFromtemperature(double T)
{
   double s;
   int entropyColnum = 3;
   int temperatureColnum = 4;
   s = getOneFromtheOther(temperatureColnum, entropyColnum, T);
   return(s);
}

double EOS::getEntropydensityFromenergyDensity(double e)
{
   double s;
   int energyColnum = 1;
   int entropyColnum = 3;
   s = getOneFromtheOther(energyColnum, entropyColnum, e);
   return(s);
}

double EOS::getEnergydensityFromtemperature(double T)
{
   double e;
   int energyColnum = 1;
   int temperatureColnum = 4;
   e = getOneFromtheOther(temperatureColnum, energyColnum, T);
   return(e);
}

double EOS::getEnergydensityFromentropyDensity(double s)
{
   double e;
   int energyColnum = 1;
   int entropyColnum = 3;
   e = getOneFromtheOther(entropyColnum, energyColnum, s);
   return(e);
}

double EOS::getTemperaturefromEntropydensity(double s)
{
   double T;
   int entropyColnum = 3;
   int temperatureColnum = 4;
   T = getOneFromtheOther(entropyColnum, temperatureColnum, s);
   return(T);
}

double EOS::getTemperaturefromEnergydensity(double e)
{
   double T;
   int energyColnum = 1;
   int temperatureColnum = 4;
   T = getOneFromtheOther(energyColnum, temperatureColnum, e);
   return(T);
}

double EOS::getPressurefromEntropydensity(double s)
{
   double p;
   int entropyColnum = 3;
   int pressureColnum = 2;
   p = getOneFromtheOther(entropyColnum, pressureColnum, s);
   return(p);
}

double EOS::getPressurefromTemperature(double T)
{
   double p;
   int temperatureColnum = 4;
   int pressureColnum = 2;
   p = getOneFromtheOther(temperatureColnum, pressureColnum, T);
   return(p);
}

double EOS::getOneFromtheOther(int xColnum, int yColnum, double xVal)
{
   double yVal;
   double xMin = EOS_epsT->getFirst(xColnum);
   double xMax = EOS_epsT->getLast(xColnum);
   double yMin = EOS_epsT->getFirst(yColnum);
  
   double e;

   string function_name = "EOS::getOneFromtheOther";
   if(xVal < 0)
   {
      int errLevel = 3;
      outputFunctionerror(function_name, "xVal < 0!", xVal, errLevel);
   }
   else if(xVal <= xMin)
   {
      int errLevel = 0;
      outputFunctionerror(function_name, "xVal < xMin, using linear extrapolation instead.", xVal, errLevel);
      yVal = xVal*yMin/xMin;
   }
   else if(xVal > xMax)
   {
      int errLevel = 0;
      outputFunctionerror(function_name, "xVal > xMax, using polynomial extrapolation instead.", xVal, errLevel);
      if(xColnum == 1)
         e = xVal;
      else
      {
         double extrapCoeff1_x = EOS_coeffs->get(1, xColnum-1);
         double extrapCoeff2_x = EOS_coeffs->get(2, xColnum-1);
         e = pow(xVal/extrapCoeff1_x, 1./extrapCoeff2_x);
      }
      if(yColnum == 1)
         yVal = e;
      else
      {
         double extrapCoeff1 = EOS_coeffs->get(1, yColnum-1);
         double extrapCoeff2 = EOS_coeffs->get(2, yColnum-1);
         yVal = extrapCoeff1*pow(e, extrapCoeff2);
      }
   }
   else
      yVal = EOS_epsT->interp(xColnum, yColnum, xVal, 2);
   return(yVal);
}
