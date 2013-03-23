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
    for(int i=0; i<nx; i++)
    {
       for(int j=0; j<ny; j++)
       {
          double sd1 = EOS_ptr->getOneFromtheOther(4, 3, frameptr->temp[i][j]);
          double sd2 = sd1*ratio;
          frameptr->temp[i][j] = EOS_ptr->getOneFromtheOther(3, 4, sd2);
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
