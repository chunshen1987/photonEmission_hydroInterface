#ifndef BJORKENEXPANSION_H
#define BJORKENEXPANSION_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>

#include "EOS.h"
#include "parameter.h"
#include "OSCARreader.h"

using namespace std;

class BjorkenExpansion
{
   private:
      EOS* EOS_ptr;

   public:
      BjorkenExpansion();
      ~BjorkenExpansion();

      void backtrace_Entropy_Bjorken_1Dlongitudinalexpansion(readindata* frameptr, double tau_1, double tau_2);

      void backtrace_Temperature_Bjorken_1Dlongitudinalexpansion(readindata* frameptr, double tau_1, double tau_2);

};

#endif
