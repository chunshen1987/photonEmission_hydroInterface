#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>

#include "BjorkenExpansion.h"
using namespace std;

BjorkenExpansion::BjorkenExpansion()
{
   EOS_ptr = new EOS(); 
}

BjorkenExpansion::~BjorkenExpansion()
{
   delete EOS_ptr;
}

void BjorkenExpansion::backtrace_Entropy_Bjorken_1Dlongitudinalexpansion(readindata* frameptr, double tau_1, double tau_2)
{
    double ratio = tau_1/tau_2;
    double tempMin = 0.01;
    for(int i=0; i<nx; i++)
    {
       for(int j=0; j<ny; j++)
       {
          if(frameptr->temp[i][j] > tempMin)
          {
             double sd1 = EOS_ptr->getEntropydensityFromtemperature(frameptr->temp[i][j]);
             double sd2 = sd1*ratio;
             frameptr->temp[i][j] = EOS_ptr->getTemperaturefromEntropydensity(sd2);
             frameptr->ed[i][j] = EOS_ptr->getEnergydensityFromentropyDensity(sd2);
             frameptr->pl[i][j] = EOS_ptr->getPressurefromEntropydensity(sd2);
          }
       }
    }
    return;
}

void BjorkenExpansion::backtrace_Temperature_Bjorken_1Dlongitudinalexpansion(readindata* frameptr, double tau_1, double tau_2)
{
    double ratio = tau_1/tau_2;
    double expon =  1./3.;
    for(int i=0; i<nx; i++)
    {
       for(int j=0; j<ny; j++)
           frameptr->temp[i][j] = frameptr->temp[i][j] * pow(ratio, expon);
    }
    return;
}
