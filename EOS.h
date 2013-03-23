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
#include "parameter.h"

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
      
      double getOneFromtheOther(int xColnum, int yColnum, double xVal);
};

#endif
