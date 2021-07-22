#ifndef SRC_QGPAMYCOLLINEAR
#define SRC_QGPAMYCOLLINEAR

#include <vector>
#include <memory>

#include "ThermalPhoton.h"
#include "ParameterReader.h"

class QGPAMYCollinear : public ThermalPhoton {
 public:
    QGPAMYCollinear(std::shared_ptr<ParameterReader> paraRdr_in,
                                std::string emissionProcess);
    ~QGPAMYCollinear() {}
    void analyticRates(double T, std::vector<double> &Eq,
                        std::vector<double> &eqrate_ptr);
};

#endif
