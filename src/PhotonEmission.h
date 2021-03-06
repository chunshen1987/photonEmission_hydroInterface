#ifndef SRC_PHOTONEMISSION_H_
#define SRC_PHOTONEMISSION_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>

#include "Hydroinfo_h5.h"
#include "Hydroinfo_MUSIC.h"
#include "ThermalPhoton.h"
#include "ParameterReader.h"

using namespace std;

class PhotonEmission {
 private:
    ParameterReader *paraRdr;
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
    double* Eq_localrest_Tb;
    double* pi_photon_Tb;
    double* bulkPi_Tb;

    double ***dNd2pTdphidy_eq, *dNd2pT_eq;
    double **vnpT_cos_eq, **vnpT_sin_eq;
    double ***dNd2pTdphidy;
    double *dNd2pT;
    double **vnpT_cos, **vnpT_sin;

    double dNdy_eq, dNdy_tot;
    double *vn_cos_eq, *vn_sin_eq;
    double *vn_cos_tot, *vn_sin_tot;

    //photon production processes
    ThermalPhoton* photon_QGP_2_to_2;
    ThermalPhoton* photon_QGP_collinear;
    ThermalPhoton* photon_HG_meson;
    ThermalPhoton* photon_HG_omega;
    ThermalPhoton* photon_HG_rho_spectralfun;
    ThermalPhoton* photon_HG_pipiBremsstrahlung;

    ThermalPhoton* photon_pirho;
    ThermalPhoton* photon_KstarK;
    ThermalPhoton* photon_piK;
    ThermalPhoton* photon_piKstar;
    ThermalPhoton* photon_pipi;
    ThermalPhoton* photon_rhoK;
    ThermalPhoton* photon_rho;
    ThermalPhoton* photon_pirho_omegat;

 public:
    PhotonEmission(ParameterReader* paraRdr_in);
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

