#ifndef PHOTONEMISSION_H
#define PHOTONEMISSION_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>

#include "Hydroinfo_h5.h"
#include "ThermalPhoton.h"
#include "parameter.h"

using namespace std;

class PhotonEmission
{
   private:
      double** lambda; // Lorentz boost transverse only
      double* Eq_localrest_Tb;
      double* pi_photon_Tb;

      double dNd2pTdphidy_hadrontot_eq[np][nphi][nrapidity];
      double dNd2pT_hadrontot_eq[np];
      double vnpT_hadrontot_cos_eq[norder][np];
      double vnpT_hadrontot_sin_eq[norder][np];
      double dNd2pTdphidy_hadrontot[np][nphi][nrapidity];
      double dNd2pT_hadrontot[np];
      double vnpT_hadrontot_cos[norder][np];
      double vnpT_hadrontot_sin[norder][np];

      double dNd2pTdphidy_eq[np][nphi][nrapidity];
      double dNd2pT_eq[np];
      double vnpT_cos_eq[norder][np];
      double vnpT_sin_eq[norder][np];
      double dNd2pTdphidy[np][nphi][nrapidity];
      double dNd2pT[np];
      double vnpT_cos[norder][np];
      double vnpT_sin[norder][np];

      //photon production processes
      ThermalPhoton photon_QGP;
      ThermalPhoton photon_HG;

      ThermalPhoton photon_pirho;
      ThermalPhoton photon_KstarK;
      ThermalPhoton photon_piK;
      ThermalPhoton photon_piKstar;
      ThermalPhoton photon_pipi;
      ThermalPhoton photon_rhoK;
      ThermalPhoton photon_rho;
      ThermalPhoton photon_pirho_omegat;

   public:
      PhotonEmission();
      ~PhotonEmission();

      void print_info();
      void InitializePhotonEmissionRateTables();
      void calPhotonemission(HydroinfoH5* hydroinfo_ptr, double* eta_ptr, double* dvolume);
      void calPhoton_total_SpMatrix();
      void calPhoton_SpvnpT_individualchannel();
      void calPhoton_total_Spvn();
      void outputPhoton_total_SpvnpT(string );
      void outputPhotonSpvn();
};

#endif
