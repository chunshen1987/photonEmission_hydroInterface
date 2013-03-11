#ifndef THERMALPHOTON_H
#define THERMALPHOTON_H

#include <string>
#include<fstream>

#include "OSCARreader.h"
#include "Arsenal.h"
#include "parameter.h"

using namespace std;


class ThermalPhoton
{
   private:
      //photon emission rate
      Table2D* Photonemission_eqrateTable_ptr;
      Table2D* Photonemission_viscous_rateTable_ptr;
      
      double** Emission_eqrateTb_ptr;
      double** Emission_viscous_rateTb_ptr;
      double* EmissionrateTb_Yidxptr;
      double EmissionrateTb_Xmin;
      double EmissionrateTb_Ymin;
      int EmissionrateTb_sizeX;
      int EmissionrateTb_sizeY;
      double EmissionrateTb_dX;
      double EmissionrateTb_dY;


      //photon spectra parameters
      string emissionProcess_name;
      double p[np], p_weight[np];
      double phi[nphi], phi_weight[nphi];
      double y[nrapidity];
      double theta[nrapidity];
      double dNd2pTdphidy_eq[np][nphi][nrapidity];
      double dNd2pTdphidy_vis[np][nphi][nrapidity];
      double dNd2pTdphidy_tot[np][nphi][nrapidity];
      double dNd2pT_eq[np], vnpT_eq[norder][np];
      double dNd2pT_vis[np], vnpT_vis[norder][np];
      double dNd2pT_tot[np], vnpT_tot[norder][np];

   public:
      ThermalPhoton();
      ~ThermalPhoton();

      void readEmissionrate(string);

      Table2D* get_eqRatetableptr() {return(Photonemission_eqrateTable_ptr);};
      Table2D* get_visRatetableptr() {return(Photonemission_viscous_rateTable_ptr);};
      double getPhotonp(int i) {return(p[i]);};
      double getPhotonphi(int i) {return(phi[i]);};
      double getPhoton_phiweight(int i) {return(phi_weight[i]);};
      double getPhotontheta(int i) {return(theta[i]);};
      double getPhotonrapidity(int i) {return(y[i]);};
      double getPhotonSpMatrix_eq(int i, int j, int k) {return(dNd2pTdphidy_eq[i][j][k]);};
      double getPhotonSpMatrix_tot(int i, int j, int k) {return(dNd2pTdphidy_tot[i][j][k]);};

      void setRateTableVarXmin(double Xmin) {EmissionrateTb_Xmin = Xmin;};
      void setRateTableVarYmin(double Ymin) {EmissionrateTb_Ymin = Ymin;};
      void setRateTabledvarX(double dX) {EmissionrateTb_dX = dX;};
      void setRateTabledvarY(double dY) {EmissionrateTb_dY = dY;};

      void getPhotonemissionRate(double* Eq, double* pi_zz, int Eq_length, double T, double* eqrate_ptr, double* visrate_ptr);

      //void calPhotonemission(double Eq, double T, double volume, int i, int j, int k);
      
      void calPhotonemission(double* Eq, double* pi_zz, int Tb_length, double T, double* volume, double fraction);
      void calPhoton_SpvnpT();
      void outputPhoton_SpvnpT();
      void interpolation2D_bilinear(double varX, double* varY, int Y_length, double** Table2D_ptr, double* results);

};
#endif
