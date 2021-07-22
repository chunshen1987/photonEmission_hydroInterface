#ifndef SRC_QGP2TO2TOTAL
#define SRC_QGP2TO2TOTAL

#include <vector>
#include <memory>

#include "ThermalPhoton.h"
#include "ParameterReader.h"

class QGP2to2Total : public ThermalPhoton {
 public:
    QGP2to2Total(std::shared_ptr<ParameterReader> paraRdr_in,
                        std::string emissionProcess);
    ~QGP2to2Total() {}
    void analyticRates(double T, std::vector<double> &Eq,
                       std::vector<double> &eqrate_ptr);
};

#endif
