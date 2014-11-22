#ifndef EOS_H
#define EOS_H

#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>

#include "Arsenal.h"
#include "Table.h"

using namespace std;

class EOS
{
   private:
      Table* EOS_epsT;
      Table* EOS_mu;
      Table* EOS_coeffs;

   public:
      EOS();
      ~EOS();
      
      double getEnergydensityFromtemperature(double T);
      double getEnergydensityFromentropyDensity(double s);
      double getEntropydensityFromtemperature(double T);
      double getEntropydensityFromenergyDensity(double e);
      double getTemperaturefromEntropydensity(double s);
      double getTemperaturefromEnergydensity(double e);
      double getPressurefromEntropydensity(double s);
      double getPressurefromTemperature(double T);
      double getOneFromtheOther(int xColnum, int yColnum, double xVal);
};

#endif
