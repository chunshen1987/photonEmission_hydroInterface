#ifndef PARAMETER_H
#define PARAMETER_H

#include <string>
using namespace std;

//path for rate files
const string rate_path = "ph_rates/";

//hydro grid information
const int nx = 261;
const int ny = 261;
const int neta = 10;
const double eta_i = 0.0;
const double eta_f = 3.0;

//number of points in photon spectra
const int np = 20;
const int nphi = 40;
const int nrapidity = 1;
const int norder = 10;  //calculate photon vn to norder

//parameter set for photon momentum
const double photon_q_i = 0.2;
const double photon_q_f = 4.0;
const double photon_phi_q_i = 0.0;
const double photon_phi_q_f = 2*M_PI;
const double photon_y_i = 0.0;
const double photon_y_f = 0.0;

//p-dependence of the delta f corrections to external particles
const double deltaf_alpha = 2.0;

//transition and freeze out temperature
const double T_dec = 0.150;
//const double T_sw = 0.180;
const double T_sw_high = 0.220;
const double T_sw_low = 0.184;

//photon emission rate table information
const double PhotonemRatetableInfo_Emin = 0.05e0;
const double PhotonemRatetableInfo_Tmin = 0.10e0;
const double PhotonemRatetableInfo_dE = 0.05e0;
const double PhotonemRatetableInfo_dT = 0.002e0;

//starting time, using Bjorken longitudinal 1D expansion from tau_start to tau_0 for hydro evolution
const double tau_start = 0.2;

#endif
