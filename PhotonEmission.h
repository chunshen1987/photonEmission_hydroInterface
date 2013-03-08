#ifndef PHOTONEMISSION_H
#define PHOTONEMISSION_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>

#include "OSCARreader.h"
#include "ThermalPhoton.h"
#include "parameter.h"

using namespace std;

class PhotonEmission
{
   private:
      double** lambda;  // Lorentz boost transverse+longitudinal
      double** lambda_transverse; // Lorentz boost transverse only
      double** RotationM; // Rotation matrix
      double* Rotation_Rzi;
      double* Eq_localrest_Tb;
      double* pi_zz_photon_Tb;

      double dNd2pTdphidy_hadrontot_eq[np][nphi][nrapidity];
      double dNd2pT_hadrontot_eq[np];
      double vnpT_hadrontot_eq[norder][np];
      double dNd2pTdphidy_hadrontot[np][nphi][nrapidity];
      double dNd2pT_hadrontot[np];
      double vnpT_hadrontot[norder][np];

      double dNd2pTdphidy_eq[np][nphi][nrapidity];
      double dNd2pT_eq[np];
      double vnpT_eq[norder][np];
      double dNd2pTdphidy[np][nphi][nrapidity];
      double dNd2pT[np];
      double vnpT[norder][np];

      //photon production processes
      ThermalPhoton photon_QGP;

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
      void calPhotonemission(readindata* frameptr,double* tanheta_ptr, double* volume);
      void calPhoton_total_SpMatrix();
      void calPhoton_SpvnpT_individualchannel();
      void calPhoton_total_Spvn();
      void outputPhoton_total_SpvnpT(string );
      void outputPhotonSpvn();
};

#endif
