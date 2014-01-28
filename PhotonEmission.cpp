#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>
#include<vector>
#include<random>

#include "Hydroinfo_h5.h"
#include "ThermalPhoton.h"
#include "tensor_trans.h"
#include "PhotonEmission.h"
#include "ParameterReader.h"

using namespace std;

PhotonEmission::PhotonEmission(ParameterReader* paraRdr_in)
{  
   paraRdr = paraRdr_in;
   output_path = "results/";
   
   set_hydroGridinfo();
   print_hydroGridinfo();

   // read the photon emission rate tables
   InitializePhotonEmissionRateTables();

   lambda = new double* [4];
   for(int i=0; i<4; i++)
   {
       lambda[i] = new double [4];
   }
   for(int i=0;i<4;i++)    //initial by 0
       for(int j=0;j<4;j++)
       {
          lambda[i][j] = 0.0e0;
       }
   
   int Eqtb_length = neta*nrapidity*np*nphi;
   Eq_localrest_Tb = new double [Eqtb_length];
   pi_photon_Tb = new double [Eqtb_length];
   for(int i=0;i<Eqtb_length;i++)
   {
     Eq_localrest_Tb[i] = 0.0;
     pi_photon_Tb[i] = 0.0;
   }

   dNd2pTdphidy_eq = new double** [np];
   dNd2pTdphidy = new double** [np];
   dNd2pT_eq = new double [np];
   dNd2pT = new double [np];
   for(int i=0; i<np; i++)
   {
      dNd2pT_eq[i] = 0.0e0;
      dNd2pT[i] = 0.0e0;
      dNd2pTdphidy_eq[i] = new double* [nphi];
      dNd2pTdphidy[i] = new double* [nphi];
      for(int j=0; j<nphi; j++)
      {
         dNd2pTdphidy_eq[i][j] = new double [nrapidity];
         dNd2pTdphidy[i][j] = new double [nrapidity];
         for(int k=0; k<nrapidity; k++)
         {
            dNd2pTdphidy_eq[i][j][k] = 0.0e0;
            dNd2pTdphidy[i][j][k] = 0.0e0;
         }
      }
   }

   vnpT_cos_eq = new double* [norder];
   vnpT_sin_eq = new double* [norder];
   vnpT_cos = new double* [norder];
   vnpT_sin = new double* [norder];
   for(int order=0; order<norder; order++)
   {
      vnpT_cos_eq[order] = new double [np];
      vnpT_sin_eq[order] = new double [np];
      vnpT_cos[order] = new double [np];
      vnpT_sin[order] = new double [np];
      for(int i = 0; i < np; i++)
      {
         vnpT_cos_eq[order][i] = 0.0e0;
         vnpT_cos[order][i] = 0.0e0;
         vnpT_sin_eq[order][i] = 0.0e0;
         vnpT_sin[order][i] = 0.0e0;
      }
   }

   dNdy_eq = 0.0;
   dNdy_tot = 0.0;
   vn_cos_eq = new double [norder];
   vn_sin_eq = new double [norder];
   vn_cos_tot = new double [norder];
   vn_sin_tot = new double [norder];
   for(int order = 0; order < norder; order++)
   {
      vn_cos_eq[order] = 0.0;
      vn_sin_eq[order] = 0.0;
      vn_cos_tot[order] = 0.0;
      vn_sin_tot[order] = 0.0;
   }

   return;
}

PhotonEmission::~PhotonEmission()
{
   for(int i = 0; i < np; i++)
   {
      for(int j = 0; j < nphi; j++)
      {
         delete[] dNd2pTdphidy_eq[i][j] ;
         delete[] dNd2pTdphidy[i][j];
      }
      delete[] dNd2pTdphidy_eq[i];
      delete[] dNd2pTdphidy[i];
   }
   delete[] dNd2pTdphidy_eq;
   delete[] dNd2pTdphidy;
   delete[] dNd2pT_eq;
   delete[] dNd2pT;

   for(int i =0; i < norder; i++)
   {
      delete[] vnpT_cos_eq[i];
      delete[] vnpT_sin_eq[i];
      delete[] vnpT_cos[i];
      delete[] vnpT_sin[i];
   }
   delete[] vnpT_cos_eq;
   delete[] vnpT_sin_eq;
   delete[] vnpT_cos;
   delete[] vnpT_sin;

   delete [] vn_cos_eq;
   delete [] vn_sin_eq;
   delete [] vn_cos_tot;
   delete [] vn_sin_tot;

   for(int i=0; i<4; i++)
   {
      delete [] lambda[i];
   }
   delete [] lambda;
   delete [] Eq_localrest_Tb;
   delete [] pi_photon_Tb;

   delete photon_QGP;
   delete photon_HG;
   delete photon_pirho;
   delete photon_pirho_omegat;
   delete photon_rho;
   delete photon_rhoK;
   delete photon_pipi;
   delete photon_piKstar;
   delete photon_piK;
   delete photon_KstarK;
   return;
}

void PhotonEmission::set_hydroGridinfo()
{
   gridX0 = paraRdr->getVal("Xmin");
   gridY0 = paraRdr->getVal("Ymin");
   gridDx = paraRdr->getVal("dx");
   gridDy = paraRdr->getVal("dy");
   gridTau0 = paraRdr->getVal("tau_start");
   gridTauf = paraRdr->getVal("tau_end");
   gridDtau = paraRdr->getVal("dTau");

   gridNx = 2*fabs(gridX0)/gridDx + 1;
   gridNy = 2*fabs(gridY0)/gridDy + 1;

   neta = paraRdr->getVal("neta");
   np = paraRdr->getVal("np");
   nphi = paraRdr->getVal("nphi");
   nrapidity = paraRdr->getVal("nrapidity");
   norder = paraRdr->getVal("norder");

   T_dec = paraRdr->getVal("T_dec");
   T_sw_high = paraRdr->getVal("T_sw_high");
   T_sw_low = paraRdr->getVal("T_sw_low");

   T_cuthigh = paraRdr->getVal("T_cuthigh");
   T_cutlow = paraRdr->getVal("T_cutlow");
   
   calHGIdFlag = paraRdr->getVal("CalHGIdFlag");
}

void PhotonEmission::print_hydroGridinfo()
{
   cout << "----------------------------------------" << endl;
   cout << "-- Parameters list for photon emission:" << endl;
   cout << "----------------------------------------" << endl;
   cout << "tau_start =" << paraRdr->getVal("tau_start") << " fm/c." << endl;
   cout << "tau_end =" << paraRdr->getVal("tau_end") << " fm/c." << endl;
   cout << "dTau = " << gridDtau << " fm/c" << endl;
   cout << "X_min = " << gridX0 << " fm/c" << endl;
   cout << "dx = " << gridDx << " fm/c" << endl;
   cout << "Y_min = " << gridY0 << " fm/c" << endl;
   cout << "dy = " << gridDy << " fm/c" << endl;
   cout << endl;

   cout << "T_dec = " << T_dec << " GeV." << endl;
   cout << "T_sw = " << T_sw_low <<  " to " << T_sw_high << " GeV."<< endl;
   cout << "deltaf_alpha = " << paraRdr->getVal("deltaf_alpha") << endl;
   cout << endl;

   cout << "Photon momentum: " << paraRdr->getVal("photon_q_i") 
        << " to " << paraRdr->getVal("photon_q_f") << " GeV, " 
        << "n_q =" << np << endl;
   cout << "Photon momentum angles: " << paraRdr->getVal("photon_phi_q_i") 
        << " to " << paraRdr->getVal("photon_phi_q_f") 
        << ", n_phi=" << nphi << endl;
   cout << "Photon momentum rapidity: " << paraRdr->getVal("photon_y_i") 
        << " to " << paraRdr->getVal("photon_y_f") 
        << ", n_y =" << nrapidity << endl;
   
   cout << "Calculate individual channels in Hadron Resonance Gas phase: ";
   if(calHGIdFlag == 0)
      cout << " No! " << endl;
   else
      cout << " Yes!" << endl;

   cout << "******************************************" << endl;
   return;
}

void PhotonEmission::InitializePhotonEmissionRateTables()
{
   double photonrate_tb_Emin = paraRdr->getVal("PhotonemRatetableInfo_Emin");
   double photonrate_tb_Tmin = paraRdr->getVal("PhotonemRatetableInfo_Tmin");
   double photonrate_tb_dE = paraRdr->getVal("PhotonemRatetableInfo_dE");
   double photonrate_tb_dT = paraRdr->getVal("PhotonemRatetableInfo_dT");
   
   photon_QGP = new ThermalPhoton(paraRdr);
   photon_QGP->setupEmissionrate("QGP_2to2_total", photonrate_tb_Tmin, photonrate_tb_dT, photonrate_tb_Emin, photonrate_tb_dE);
   photon_HG = new ThermalPhoton(paraRdr);
   photon_HG->setupEmissionrate("HG_2to2_total", photonrate_tb_Tmin, photonrate_tb_dT, photonrate_tb_Emin, photonrate_tb_dE);

   photon_pirho = new ThermalPhoton(paraRdr);
   photon_pirho->setupEmissionrate("pion_rho_to_pion_gamma", photonrate_tb_Tmin, photonrate_tb_dT, photonrate_tb_Emin, photonrate_tb_dE);
   photon_KstarK = new ThermalPhoton(paraRdr);
   photon_KstarK->setupEmissionrate("K_Kstar_to_pion_gamma", photonrate_tb_Tmin, photonrate_tb_dT, photonrate_tb_Emin, photonrate_tb_dE);
   photon_piK = new ThermalPhoton(paraRdr);
   photon_piK->setupEmissionrate("pion_K_to_Kstar_gamma", photonrate_tb_Tmin, photonrate_tb_dT, photonrate_tb_Emin, photonrate_tb_dE);
   photon_piKstar = new ThermalPhoton(paraRdr);
   photon_piKstar->setupEmissionrate("pion_Kstar_to_K_gamma",photonrate_tb_Tmin, photonrate_tb_dT, photonrate_tb_Emin, photonrate_tb_dE);
   photon_pipi = new ThermalPhoton(paraRdr);
   photon_pipi->setupEmissionrate("pion_pion_to_rho_gamma", photonrate_tb_Tmin, photonrate_tb_dT, photonrate_tb_Emin, photonrate_tb_dE);
   photon_rhoK = new ThermalPhoton(paraRdr);
   photon_rhoK->setupEmissionrate("rho_K_to_K_gamma", photonrate_tb_Tmin, photonrate_tb_dT, photonrate_tb_Emin, photonrate_tb_dE);
   photon_rho = new ThermalPhoton(paraRdr);
   photon_rho->setupEmissionrate("rho_to_pion_pion_gamma", photonrate_tb_Tmin, photonrate_tb_dT, photonrate_tb_Emin, photonrate_tb_dE);
   photon_pirho_omegat = new ThermalPhoton(paraRdr);
   photon_pirho_omegat->setupEmissionrate("pion_rho_to_omega_to_pion_gamma",photonrate_tb_Tmin, photonrate_tb_dT, photonrate_tb_Emin, photonrate_tb_dE);

   return;
}

void PhotonEmission::calPhotonemission(HydroinfoH5* hydroinfo_ptr, double* eta_ptr, double* etaweight_ptr)
{
  //photon momentum in the lab frame
  double p_q[np], phi_q[nphi], y_q[nrapidity];
  double sin_phiq[nphi], cos_phiq[nphi];
  double p_lab_local[4], p_lab_lowmu[4];
  double flow_u_mu_low[4];
  for(int k=0;k<nrapidity;k++)
     y_q[k] = photon_QGP->getPhotonrapidity(k);
  for(int l=0;l<np;l++) p_q[l] = photon_QGP->getPhotonp(l);
  for(int m=0;m<nphi;m++)
  {
     phi_q[m] = photon_QGP->getPhotonphi(m);
     sin_phiq[m] = sin(phi_q[m]);
     cos_phiq[m] = cos(phi_q[m]);
  }

  double e_local, p_local, temp_local, vx_local, vy_local;
  double tau_local;
  double eta_local;
  double* volume = new double [neta];
  double** pi_tensor_lab = new double* [4];
  for(int i=0; i<4; i++)
    pi_tensor_lab[i] = new double [4];

  fluidCell* fluidCellptr = new fluidCell();

  //main loop begins ...
  //loop over time frame
  if(gridTauf > hydroinfo_ptr->getHydrogridTaumax()) 
     gridTauf = hydroinfo_ptr->getHydrogridTaumax();
  int nFrame = (int)((gridTauf - gridTau0)/gridDtau) + 1;
  for(int frameId = 0; frameId < nFrame; frameId++)
  {
     tau_local = gridTau0 + frameId*gridDtau;
     for(int k=0; k<neta; k++)
        volume[k] = 2 * tau_local * gridDx * gridDy * gridDtau * etaweight_ptr[k]; //volume element: tau*dtau*dx*dy*deta, 2 for symmetry along longitudinal direction
  //loops over the transverse plane
     for(int i=0; i < gridNx; i++)
     {
       double x_local = gridX0 + i*gridDx;
       for(int j=0; j < gridNy; j++)
       {
         double y_local = gridY0 + j*gridDy;

         int idx_Tb = 0;
         hydroinfo_ptr->getHydroinfo(tau_local, x_local, y_local, fluidCellptr);
         temp_local = fluidCellptr->temperature;
         if(temp_local > T_dec && temp_local < T_cuthigh && temp_local > T_cutlow)
         {
           e_local = fluidCellptr->ed;
           p_local = fluidCellptr->pressure;
           vx_local = fluidCellptr->vx;
           vy_local = fluidCellptr->vy;
           for(int mu = 0; mu < 4; mu++)
              for(int nu = 0; nu < 4; nu++)
                 pi_tensor_lab[mu][nu] = fluidCellptr->pi[mu][nu];

           getTransverseflow_u_mu_low(flow_u_mu_low, vx_local, vy_local);
           double prefactor_pimunu = 1./(2.*(e_local + p_local));
           for(int jj=0; jj<neta; jj++)
           {
             eta_local = eta_ptr[jj];
             
             //photon momentum loops
             for(int k=0;k<nrapidity;k++) 
             {
                double cosh_y_minus_eta = cosh(y_q[k] - eta_local);
                double sinh_y_minus_eta = sinh(y_q[k] - eta_local);
             for(int m=0;m<nphi;m++)
             {
             for(int l=0;l<np;l++)
             { 
               p_lab_local[0] = p_q[l]*cosh_y_minus_eta;
               p_lab_local[1] = p_q[l]*cos_phiq[m];
               p_lab_local[2] = p_q[l]*sin_phiq[m];
               p_lab_local[3] = p_q[l]*sinh_y_minus_eta;
               p_lab_lowmu[0] = p_lab_local[0];
               for(int local_i = 1; local_i < 4; local_i++)
                  p_lab_lowmu[local_i] = - p_lab_local[local_i];

               double Eq_localrest_temp = 0.0e0;
               double pi_photon = 0.0e0;
               for(int local_i = 0; local_i < 4; local_i++)
                  Eq_localrest_temp += flow_u_mu_low[local_i]*p_lab_local[local_i];
              
               for(int local_i = 0; local_i < 4; local_i++)
                  for(int local_j = 0; local_j < 4; local_j++)
                     pi_photon += pi_tensor_lab[local_i][local_j]*p_lab_lowmu[local_i]*p_lab_lowmu[local_j];

               Eq_localrest_Tb[idx_Tb] = Eq_localrest_temp;
               pi_photon_Tb[idx_Tb] = pi_photon*prefactor_pimunu;
               idx_Tb++;
             }
             }
             }
           }
           if(temp_local > T_sw_high)
           {
             double QGP_fraction = 1.0;
             photon_QGP->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local, volume, QGP_fraction);
           }
           else if(temp_local > T_sw_low)
           {
             double QGP_fraction = (temp_local - T_sw_low)/(T_sw_high - T_sw_low);
             double HG_fraction = 1 - QGP_fraction;
             photon_QGP->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local, volume, QGP_fraction);
             photon_HG->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local, volume, HG_fraction);
             if(calHGIdFlag == 1)
             {
                photon_pirho->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
                photon_KstarK->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
                photon_piK->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
                photon_piKstar->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
                photon_pipi->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
                photon_rhoK->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
                photon_rho->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
                photon_pirho_omegat->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
             }
           }
           else
           {
             double HG_fraction = 1.0;
             photon_HG->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local, volume, HG_fraction);
             if(calHGIdFlag == 1)
             {
                photon_pirho->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
                photon_KstarK->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
                photon_piK->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
                photon_piKstar->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
                photon_pipi->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
                photon_rhoK->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
                photon_rho->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
                photon_pirho_omegat->calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
             }
           }
         }
       }
     }
     cout<<"frame "<< frameId << " : ";
     cout<<" tau = " << setw(4) << setprecision(3) << tau_local 
         <<" fm/c done!" << endl;
  }

  delete [] volume;
  delete fluidCellptr;
  for(int i=0; i<4; i++)
  {
     delete [] pi_tensor_lab[i];
  }
  delete [] pi_tensor_lab;

  return;
}

void PhotonEmission::calphotonemissionFunction(HydroinfoH5* hydroinfo_ptr, double* eta_ptr, double* etaweight_ptr)
{
  int Tblength = np*nrapidity*neta;
  double* Eq_localrest_table = new double [Tblength];
  double* pi_photon_table = new double [Tblength];
  vector<emissionFunctionmovie> emissionArray;

  //photon momentum in the lab frame
  double p_q[np], phi_q[nphi], y_q[nrapidity];
  double phi_q_weight[nphi];
  double sin_phiq[nphi], cos_phiq[nphi];
  double p_lab_local[4], p_lab_lowmu[4];
  double flow_u_mu_low[4];
  for(int k=0;k<nrapidity;k++)
     y_q[k] = photon_QGP->getPhotonrapidity(k);
  for(int l=0;l<np;l++) p_q[l] = photon_QGP->getPhotonp(l);
  for(int m=0;m<nphi;m++)
  {
     phi_q[m] = photon_QGP->getPhotonphi(m);
     phi_q_weight[m] = photon_QGP->getPhoton_phiweight(m);
     sin_phiq[m] = sin(phi_q[m]);
     cos_phiq[m] = cos(phi_q[m]);
  }

  double e_local, p_local, temp_local, vx_local, vy_local;
  double tau_local;
  double eta_local;
  double* volume = new double [neta];
  double** pi_tensor_lab = new double* [4];
  for(int i=0; i<4; i++)
    pi_tensor_lab[i] = new double [4];

  fluidCell* fluidCellptr = new fluidCell();

  //main loop begins ...
  //loop over time frame
  if(gridTauf > hydroinfo_ptr->getHydrogridTaumax()) 
     gridTauf = hydroinfo_ptr->getHydrogridTaumax();
  int nFrame = (int)((gridTauf - gridTau0)/gridDtau) + 1;
  for(int frameId = 0; frameId < nFrame; frameId++)
  {
     tau_local = gridTau0 + frameId*gridDtau;
     for(int k=0; k<neta; k++)
        volume[k] = 2 * tau_local * gridDx * gridDy * gridDtau * etaweight_ptr[k]; //volume element: tau*dtau*dx*dy*deta, 2 for symmetry along longitudinal direction
  //loops over the transverse plane
     for(int i=0; i < gridNx; i++)
     {
       double x_local = gridX0 + i*gridDx;
       for(int j=0; j < gridNy; j++)
       {
         double y_local = gridY0 + j*gridDy;

         hydroinfo_ptr->getHydroinfo(tau_local, x_local, y_local, fluidCellptr);
         temp_local = fluidCellptr->temperature;
         if(temp_local > T_dec && temp_local < T_cuthigh && temp_local > T_cutlow)
         {
           e_local = fluidCellptr->ed;
           p_local = fluidCellptr->pressure;
           vx_local = fluidCellptr->vx;
           vy_local = fluidCellptr->vy;
           for(int mu = 0; mu < 4; mu++)
              for(int nu = 0; nu < 4; nu++)
                 pi_tensor_lab[mu][nu] = fluidCellptr->pi[mu][nu];

           getTransverseflow_u_mu_low(flow_u_mu_low, vx_local, vy_local);
           double prefactor_pimunu = 1./(2.*(e_local + p_local));
           for(int m=0;m<nphi;m++)
           {
              int idx_Tb = 0;
              for(int jj=0; jj<neta; jj++)
              {
                eta_local = eta_ptr[jj];
                //photon momentum loops
                for(int k=0;k<nrapidity;k++) 
                {
                   double cosh_y_minus_eta = cosh(y_q[k] - eta_local);
                   double sinh_y_minus_eta = sinh(y_q[k] - eta_local);
                   for(int l=0;l<np;l++)
                   { 
                     p_lab_local[0] = p_q[l]*cosh_y_minus_eta;
                     p_lab_local[1] = p_q[l]*cos_phiq[m];
                     p_lab_local[2] = p_q[l]*sin_phiq[m];
                     p_lab_local[3] = p_q[l]*sinh_y_minus_eta;
                     p_lab_lowmu[0] = p_lab_local[0];
                     for(int local_i = 1; local_i < 4; local_i++)
                        p_lab_lowmu[local_i] = - p_lab_local[local_i];

                     double Eq_localrest_temp = 0.0e0;
                     double pi_photon = 0.0e0;
                     for(int local_i = 0; local_i < 4; local_i++)
                        Eq_localrest_temp += flow_u_mu_low[local_i]*p_lab_local[local_i];
                    
                     for(int local_i = 0; local_i < 4; local_i++)
                        for(int local_j = 0; local_j < 4; local_j++)
                           pi_photon += pi_tensor_lab[local_i][local_j]*p_lab_lowmu[local_i]*p_lab_lowmu[local_j];

                     Eq_localrest_table[idx_Tb] = Eq_localrest_temp;
                     pi_photon_table[idx_Tb] = pi_photon*prefactor_pimunu;
                     idx_Tb++;
                   }
                }
              }
              double result = 0.0;
              if(temp_local > T_sw_high)
              {
                double QGP_fraction = 1.0;
                result = photon_QGP->getThermalPhotonemissionFunctioninte(Eq_localrest_table, pi_photon_table, idx_Tb, temp_local, volume, QGP_fraction);
              }
              else if(temp_local > T_sw_low)
              {
                double QGP_fraction = (temp_local - T_sw_low)/(T_sw_high - T_sw_low);
                double HG_fraction = 1 - QGP_fraction;
                result = photon_QGP->getThermalPhotonemissionFunctioninte(Eq_localrest_table, pi_photon_table, idx_Tb, temp_local, volume, QGP_fraction);
                result += photon_HG->getThermalPhotonemissionFunctioninte(Eq_localrest_table, pi_photon_table, idx_Tb, temp_local, volume, HG_fraction);
              }
              else
              {
                double HG_fraction = 1.0;
                result = photon_HG->getThermalPhotonemissionFunctioninte(Eq_localrest_table, pi_photon_table, idx_Tb, temp_local, volume, HG_fraction);
              }
              emissionFunctionmovie emissionFunctiontemp;
              emissionFunctiontemp.tau = tau_local;
              emissionFunctiontemp.x = x_local;
              emissionFunctiontemp.y = y_local;
              emissionFunctiontemp.phi_q = phi_q[m];
              emissionFunctiontemp.S_q = result*phi_q_weight[m];
              emissionArray.push_back(emissionFunctiontemp);
            }
         }
       }
     }
     cout<<"frame "<< frameId << " : ";
     cout<<" tau = " << setw(4) << setprecision(3) << tau_local 
         <<" fm/c done!" << endl;
  }

  cout << "start Monte-Carlo sampling photons ..." << endl;
  int oversample = 10;
  double dNgamma = 0.0;
  int emissionArrayLength = emissionArray.size();
  //building inverse CDF for the emission function
  vector<double> invCDF(emissionArrayLength);
  for(int i = 0; i < emissionArrayLength; i++)
  {
     dNgamma += emissionArray[i].S_q;
     invCDF[i] = dNgamma;
  }
  for(int i = 0; i < emissionArrayLength; i++)
     invCDF[i] = invCDF[i]/dNgamma;
  cout << "total number of photon is: " << dNgamma << " oversample = " << oversample << endl;

  int dNsample = (int)(dNgamma*oversample);

  //generate uniform distributed random numbers
  default_random_engine generator;
  uniform_real_distribution<double> randomDistribution(0.0,1.0);
 
  ofstream outputSample("results/Samples.dat");

  //start Monte-Carlo sampling
  double Eq = 1.0;
  for(int i = 0; i < dNsample; i++)
  {
     double randNumber = randomDistribution(generator);
     int idx = binarySearch(&invCDF, randNumber);
     outputSample << scientific << setw(18) << setprecision(8) 
        << emissionArray[idx].tau << "   " << emissionArray[idx].x << "   " 
        << emissionArray[idx].y << "   " << Eq << "   " << Eq << "   "
        << emissionArray[idx].phi_q << endl;
  }
  outputSample.close();

  delete [] volume;
  delete [] Eq_localrest_table;
  delete [] pi_photon_table;
  delete fluidCellptr;
  for(int i=0; i<4; i++)
  {
     delete [] pi_tensor_lab[i];
  }
  delete [] pi_tensor_lab;

  return;
}


void PhotonEmission::calPhoton_total_SpMatrix()
{
     for(int k=0;k<nrapidity;k++) 
     {
       for(int l=0;l<np;l++)
       { 
         for(int m=0;m<nphi;m++)
         {
           dNd2pTdphidy_eq[l][m][k] =  photon_QGP->getPhotonSpMatrix_eq(l, m, k) + photon_HG->getPhotonSpMatrix_eq(l, m, k);
           dNd2pTdphidy[l][m][k] =  photon_QGP->getPhotonSpMatrix_tot(l, m, k) + photon_HG->getPhotonSpMatrix_tot(l, m, k);
         }
       }
     }
     return;
}

void PhotonEmission::calPhoton_SpvnpT_individualchannel()
{
    photon_QGP->calPhoton_SpvnpT();
    photon_HG->calPhoton_SpvnpT();
    if(calHGIdFlag == 1)
    {
       photon_pirho->calPhoton_SpvnpT();
       photon_KstarK->calPhoton_SpvnpT();
       photon_piK->calPhoton_SpvnpT();
       photon_piKstar->calPhoton_SpvnpT();
       photon_pipi->calPhoton_SpvnpT();
       photon_rhoK->calPhoton_SpvnpT();
       photon_rho->calPhoton_SpvnpT();
       photon_pirho_omegat->calPhoton_SpvnpT();
    }

    return;
}

void PhotonEmission::outputPhotonSpvn()
{
    photon_QGP->outputPhoton_SpvnpT(output_path);
    photon_HG->outputPhoton_SpvnpT(output_path);
    if(calHGIdFlag == 1)
    {
       photon_pirho->outputPhoton_SpvnpT(output_path);
       photon_KstarK->outputPhoton_SpvnpT(output_path);
       photon_piK->outputPhoton_SpvnpT(output_path);
       photon_piKstar->outputPhoton_SpvnpT(output_path);
       photon_pipi->outputPhoton_SpvnpT(output_path);
       photon_rhoK->outputPhoton_SpvnpT(output_path);
       photon_rho->outputPhoton_SpvnpT(output_path);
       photon_pirho_omegat->outputPhoton_SpvnpT(output_path);
    }

    outputPhoton_total_SpvnpT("photon_total");
    return;
}

void PhotonEmission::calPhoton_total_Spvn()
{
   int k = 0;
   for(int i=0;i<np;i++)
   {
       double p = photon_QGP->getPhotonp(i);
       double pweight = photon_QGP->getPhoton_pweight(i);
       for(int j=0;j<nphi;j++)
       {
         double phi = photon_QGP->getPhotonphi(j);
         double phiweight = photon_QGP->getPhoton_phiweight(j);
         dNd2pT_eq[i] += dNd2pTdphidy_eq[i][j][k]*phiweight;
         dNd2pT[i] += dNd2pTdphidy[i][j][k]*phiweight;
         for(int order = 0; order < norder; order++)
         {
            vnpT_cos_eq[order][i] += dNd2pTdphidy_eq[i][j][k]*cos(order*phi)*phiweight;
            vnpT_cos[order][i] += dNd2pTdphidy[i][j][k]*cos(order*phi)*phiweight;
            vnpT_sin_eq[order][i] += dNd2pTdphidy_eq[i][j][k]*sin(order*phi)*phiweight;
            vnpT_sin[order][i] += dNd2pTdphidy[i][j][k]*sin(order*phi)*phiweight;
         }
       }
       dNdy_eq += dNd2pT_eq[i]*p*pweight;
       dNdy_tot += dNd2pT[i]*p*pweight;

       for(int order=0; order<norder; order++)
       {
          vn_cos_eq[order] += vnpT_cos_eq[order][i]*p*pweight;
          vn_sin_eq[order] += vnpT_sin_eq[order][i]*p*pweight;
          vn_cos_tot[order] += vnpT_cos[order][i]*p*pweight;
          vn_sin_tot[order] += vnpT_sin[order][i]*p*pweight;

          vnpT_cos_eq[order][i] = vnpT_cos_eq[order][i]/dNd2pT_eq[i];
          vnpT_cos[order][i] = vnpT_cos[order][i]/dNd2pT[i];
          vnpT_sin_eq[order][i] = vnpT_sin_eq[order][i]/dNd2pT_eq[i];
          vnpT_sin[order][i] = vnpT_sin[order][i]/dNd2pT[i];
       }
       dNd2pT_eq[i] = dNd2pT_eq[i]/(2*M_PI);
       dNd2pT[i] = dNd2pT[i]/(2*M_PI);
   }
   for(int order = 1; order < norder ; order++)
   {
       vn_cos_eq[order] = vn_cos_eq[order]/dNdy_eq;
       vn_sin_eq[order] = vn_sin_eq[order]/dNdy_eq;
       vn_cos_tot[order] = vn_cos_tot[order]/dNdy_tot;
       vn_sin_tot[order] = vn_sin_tot[order]/dNdy_tot;
   }
   return;
}

void PhotonEmission::outputPhoton_total_SpvnpT(string filename)
{
    ostringstream filename_stream_eq_SpMatrix;
    ostringstream filename_stream_eq_Spvn;
    ostringstream filename_stream_SpMatrix;
    ostringstream filename_stream_Spvn;
    ostringstream filename_stream_inte_eq_Spvn;
    ostringstream filename_stream_inte_Spvn;

    filename_stream_eq_SpMatrix << output_path << filename << "_eq_SpMatrix.dat";
    filename_stream_eq_Spvn << output_path << filename << "_eq_Spvn.dat";
    filename_stream_SpMatrix << output_path << filename << "_SpMatrix.dat";
    filename_stream_Spvn << output_path << filename << "_Spvn.dat";
    filename_stream_inte_eq_Spvn << output_path << filename << "_eq_Spvn_inte.dat";
    filename_stream_inte_Spvn << output_path << filename << "_Spvn_inte.dat";

    ofstream fphoton_eq_SpMatrix(filename_stream_eq_SpMatrix.str().c_str());
    ofstream fphoton_eq_Spvn(filename_stream_eq_Spvn.str().c_str());
    ofstream fphotonSpMatrix(filename_stream_SpMatrix.str().c_str());
    ofstream fphotonSpvn(filename_stream_Spvn.str().c_str());
    ofstream fphotoninte_eq_Spvn(filename_stream_inte_eq_Spvn.str().c_str());
    ofstream fphotoninteSpvn(filename_stream_inte_Spvn.str().c_str());

    for(int i=0;i<nphi;i++)
    {
      double phi = photon_QGP->getPhotonphi(i);
      fphoton_eq_SpMatrix << phi << "  ";
      fphotonSpMatrix << phi << "  ";
      for(int j=0;j<np;j++)
        for(int k=0;k<nrapidity;k++)
        {
           fphoton_eq_SpMatrix << scientific << setprecision(6) << setw(16) 
                           << dNd2pTdphidy_eq[j][i][k] << "  ";
           fphotonSpMatrix << scientific << setprecision(6) << setw(16) 
                           << dNd2pTdphidy[j][i][k] << "  ";
        }
      fphoton_eq_SpMatrix << endl;
      fphotonSpMatrix << endl;
    }
    for(int i=0;i<np;i++)
    {
      double pT = photon_QGP->getPhotonp(i);
      fphoton_eq_Spvn << scientific << setprecision(6) << setw(16) 
                  << pT << "  " << dNd2pT_eq[i] << "  " ;
      fphotonSpvn << scientific << setprecision(6) << setw(16) 
                  << pT << "  " << dNd2pT[i] << "  " ;
      for(int order=1; order<norder; order++)
      {
         fphoton_eq_Spvn << scientific << setprecision(6) << setw(16) 
                     << vnpT_cos_eq[order][i] << "  "
                     << vnpT_sin_eq[order][i] << "  "
                     << sqrt(pow(vnpT_cos_eq[order][i], 2) + pow(vnpT_sin_eq[order][i], 2)) << "  ";
         fphotonSpvn << scientific << setprecision(6) << setw(16) 
                     << vnpT_cos[order][i] << "  "
                     << vnpT_sin[order][i] << "  "
                     << sqrt(pow(vnpT_cos[order][i], 2) + pow(vnpT_sin[order][i], 2)) << "  ";
      }
      fphoton_eq_Spvn << endl;
      fphotonSpvn << endl;
    }

    for(int order = 0; order < norder; order++)
    {
       fphotoninte_eq_Spvn << scientific << setprecision(6) << setw(16)
                           << order << "   " << vn_cos_eq[order] << "   " << vn_sin_eq[order] << "   " 
                           << sqrt(pow(vn_cos_eq[order], 2) + pow(vn_sin_eq[order], 2)) << endl;
       fphotoninteSpvn << scientific << setprecision(6) << setw(16)
                       << order << "   " << vn_cos_tot[order] << "   " << vn_sin_tot[order] << "   " 
                       << sqrt(pow(vn_cos_tot[order], 2) + pow(vn_sin_tot[order], 2)) << endl;
    }
    return;
}
