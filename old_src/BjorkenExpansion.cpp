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

void BjorkenExpansion::backtrace_Entropy_Bjorken_1Dlongitudinalexpansion(HydroinfoH5* hydroinfo_ptr, double tau_1, double tau_2)
{
    double ratio = tau_1/tau_2;
    double tempMin = 0.01;
    for(int i=0; i<hydroinfo_ptr->getHydrogridNX(); i++)
    {
       for(int j=0; j<hydroinfo_ptr->getHydrogridNY(); j++)
       {
          double temp_local;
          fluidCell* fluidCellptr = new fluidCell() ;
          hydroinfo_ptr->getHydroinfoOnlattice(0, i, j, fluidCellptr);
          temp_local = fluidCellptr->temperature;

          if(temp_local > tempMin)
          {
             double sd1 = fluidCellptr->sd;
             double sd2 = sd1*ratio;
             fluidCellptr->temperature = EOS_ptr->getTemperaturefromEntropydensity(sd2);
             fluidCellptr->ed = EOS_ptr->getEnergydensityFromentropyDensity(sd2);
             fluidCellptr->pressure = EOS_ptr->getPressurefromEntropydensity(sd2);
          }

          delete fluidCellptr;
       }
    }
    return;
}

void BjorkenExpansion::backtrace_Temperature_Bjorken_1Dlongitudinalexpansion(HydroinfoH5* hydroinfo_ptr, double tau_1, double tau_2)
{
    double ratio = tau_1/tau_2;
    double expon =  1./3.;
    for(int i=0; i<hydroinfo_ptr->getHydrogridNX(); i++)
    {
       for(int j=0; j<hydroinfo_ptr->getHydrogridNY(); j++)
       {
          double temp_local;
          fluidCell* fluidCellptr = new fluidCell() ;
          hydroinfo_ptr->getHydroinfoOnlattice(0, i, j, fluidCellptr);
          temp_local = fluidCellptr->temperature;
          fluidCellptr->temperature = temp_local * pow(ratio, expon);
          delete fluidCellptr;
       }
    }
    return;
}
