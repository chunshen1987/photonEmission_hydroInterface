#ifndef SRC_PHOTONEMISSION_H_
#define SRC_PHOTONEMISSION_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>
#include <memory>
#include <vector>

#include "Hydroinfo_h5.h"
#include "Hydroinfo_MUSIC.h"
#include "ThermalPhoton.h"
#include "ParameterReader.h"

using namespace std;

class PhotonEmission {
 private:
    std::shared_ptr<ParameterReader> paraRdr;
    string output_path;

    int neta;
    int np, nphi, nrapidity;
    int norder;

    double gridDx, gridDy, gridDtau;
    double gridX0, gridY0, gridTau0, gridTauf;
    int gridNx, gridNy;

    double T_dec, T_sw_high, T_sw_low;
    double T_cuthigh, T_cutlow;

    int hydro_flag;
    int differential_flag;
    int turn_off_transverse_flow;
    int calHGIdFlag;

    double** lambda; // Lorentz boost transverse only
    double *Eq_localrest_Tb;
    double *pi_photon_Tb;
    double *bulkPi_Tb;
    
    double *dNd2pT_eq;
    double *dNd2pT;
    double ***dNd2pTdphidy_eq;
    double **vnpT_cos_eq, **vnpT_sin_eq;
    double ***dNd2pTdphidy;
    double **vnpT_cos, **vnpT_sin;

    double dNdy_eq, dNdy_tot;
    std::vector<double> vn_sin_eq;
    std::vector<double> vn_cos_eq;
    std::vector<double> vn_cos_tot;
    std::vector<double> vn_sin_tot;

    //photon production processes
    std::unique_ptr<ThermalPhoton> photon_QGP_2_to_2;
    std::unique_ptr<ThermalPhoton> photon_QGP_collinear;
    std::unique_ptr<ThermalPhoton> photon_HG_meson;
    std::unique_ptr<ThermalPhoton> photon_HG_omega;
    std::unique_ptr<ThermalPhoton> photon_HG_rho_spectralfun;
    std::unique_ptr<ThermalPhoton> photon_HG_pipiBremsstrahlung;

    std::unique_ptr<ThermalPhoton> photon_pirho;
    std::unique_ptr<ThermalPhoton> photon_KstarK;
    std::unique_ptr<ThermalPhoton> photon_piK;
    std::unique_ptr<ThermalPhoton> photon_piKstar;
    std::unique_ptr<ThermalPhoton> photon_pipi;
    std::unique_ptr<ThermalPhoton> photon_rhoK;
    std::unique_ptr<ThermalPhoton> photon_rho;
    std::unique_ptr<ThermalPhoton> photon_pirho_omegat;


 public:
    PhotonEmission(std::shared_ptr<ParameterReader> paraRdr_in);
    ~PhotonEmission();

    void set_hydroGridinfo();
    void print_hydroGridinfo();
    void InitializePhotonEmissionRateTables();
    void calPhotonemission(void *hydroinfo_ptr_in, double* eta_ptr, 
                           double* etaweight_ptr);
    void calPhotonemission_3d(void *hydroinfo_ptr_in);
    void calPhoton_total_SpMatrix();
    void calPhoton_SpvnpT_individualchannel();
    void calPhoton_total_Spvn();
    void outputPhoton_total_SpvnpT(string );
    void outputPhotonSpvn();
};

#endif   // SRC_PHOTONEMISSION_H_

