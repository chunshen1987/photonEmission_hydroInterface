#ifndef SRC_THERMALDILEPTON_H_
#define SRC_THERMALDILEPTON_H_

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "Arsenal.h"
#include "ParameterReader.h"
#include "Table2D.h"

class ThermalDilepton {
  private:
    std::shared_ptr<ParameterReader> paraRdr;

    int np, nphi, nrapidity, nM;
    int norder;
    int neta;
    std::string rate_path_;

    double dy;

    bool bRateTable_;

    double ***Emission_eqrateTb_ptr;        // muB/T, k/T, M/T

    std::vector<double> EmissionrateTb_Yidxptr;
    double EmissionrateTb_Xmin;
    double EmissionrateTb_Ymin;
    int EmissionrateTb_sizeX;
    int EmissionrateTb_sizeY;
    double EmissionrateTb_dX;
    double EmissionrateTb_dY;

    // photon spectra parameters
    std::string emissionProcess_name;
    double *p, *p_weight;
    double *phi, *phi_weight;
    std::vector<double> y;
    std::vector<double> theta;
    std::vector<double> Minv_;

    double ****dNd2pTdphidydM_eq;

    double **vnpT_cos_eq, **vnpT_sin_eq;
    double **vnM_cos_eq, **vnM_sin_eq;

    double ***vnypT_cos_eq, ***vnypT_sin_eq;

    std::vector<double> vn_cos_eq;
    std::vector<double> vn_sin_eq;


  public:
    ThermalDilepton(
        std::shared_ptr<ParameterReader> paraRdr_in,
        std::string emissionProcess);

    virtual ~ThermalDilepton();

    void setupEmissionrateFromFile(
        double Xmin, double dX, double Ymin, double dY, bool bShearVisCorr,
        bool bBulkVisCorr);
    void readEmissionrate(std::string);

    void setupEmissionrateFromParametrization(
        double Xmin, double dX, int nX, double Ymin, double dY, int nY);

    double get_dy() { return (dy); }

    double getPhotonp(int i) { return (p[i]); }
    double getPhoton_pweight(int i) { return (p_weight[i]); }
    double getPhotonphi(int i) { return (phi[i]); }
    double getPhoton_phiweight(int i) { return (phi_weight[i]); }
    double getPhotontheta(int i) { return (theta[i]); }
    double getPhotonrapidity(int i) { return (y[i]); }

    virtual void analyticRates(
        double T, std::vector<double> &Eq, std::vector<double> &eqrate_ptr);
    virtual void NetBaryonCorrection(
        double T, double muB, std::vector<double> &Eq,
        std::vector<double> &eqrate_ptr) {}

    void checkAnalyticRates();

    void getEmissionRate(
        std::vector<double> &Eq,
        const double Minv, const double T, const double muB,
        std::vector<double> &eqrate_ptr);
    void calThermalPhotonemission(
        std::vector<double> &Eq, int Tb_length, const double Minv, double T,
        std::vector<double> &volume, double fraction);
    void calThermalPhotonemission_3d(
        std::vector<double> &Eq,
        const double Minv, double T, double muB, double volume,
        double fraction);

    void calPhoton_SpvnpT(
        double ***dNd2pTdphipy, double ***vnypT_cos, double ***vnypT_sin,
        double **vnpT_cos, double **vnpT_sin, std::vector<double> &vn_cos,
        std::vector<double> &vn_sin);
    void calPhoton_SpvnpT_shell();
    void outputPhoton_SpvnpT(
        std::string path, std::string type_str, double ***dNd2pTdphidy,
        double ***vnypT_cos, double ***vnypT_sin, double **vnpT_cos,
        double **vnpT_sin, std::vector<double> &vn_cos,
        std::vector<double> &vn_sin);
    void outputPhoton_SpvnpT_shell(std::string path);
    void interpolation2D_bilinear(
        double varX, std::vector<double> &varY, double **Table2D_ptr,
        std::vector<double> &results);
};

#endif  // SRC_THERMALDILEPTON_H
