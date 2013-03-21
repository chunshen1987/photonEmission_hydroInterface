#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>

#include "OSCARreader.h"
#include "ThermalPhoton.h"
#include "tensor_trans.h"
#include "PhotonEmission.h"
#include "parameter.h"

using namespace std;

PhotonEmission::PhotonEmission()
{  
   print_info();

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

   for(int i=0; i<np; i++)
   {
      dNd2pT_hadrontot_eq[i] = 0.0e0;
      dNd2pT_hadrontot[i] = 0.0e0;
      dNd2pT_eq[i] = 0.0e0;
      dNd2pT[i] = 0.0e0;
      for(int order=0; order<norder; order++)
      {
         vnpT_cos_eq[order][i] = 0.0e0;
         vnpT_hadrontot_cos[order][i] = 0.0e0;
         vnpT_hadrontot_cos_eq[order][i] = 0.0e0;
         vnpT_cos[order][i] = 0.0e0;
         vnpT_sin_eq[order][i] = 0.0e0;
         vnpT_hadrontot_sin[order][i] = 0.0e0;
         vnpT_hadrontot_sin_eq[order][i] = 0.0e0;
         vnpT_sin[order][i] = 0.0e0;
      }
      for(int j=0; j<nphi; j++)
         for(int k=0; k<nrapidity; k++)
         {
            dNd2pTdphidy_hadrontot_eq[i][j][k] = 0.0e0;
            dNd2pTdphidy_eq[i][j][k] = 0.0e0;
            dNd2pTdphidy_hadrontot[i][j][k] = 0.0e0;
            dNd2pTdphidy[i][j][k] = 0.0e0;
         }
   }

   return;
}

PhotonEmission::~PhotonEmission()
{
   for(int i=0; i<4; i++)
   {
      delete [] lambda[i];
   }
   delete [] lambda;
   delete [] Eq_localrest_Tb;
   delete [] pi_photon_Tb;
   return;
}

void PhotonEmission::print_info()
{
   cout << "----------------------------------------" << endl;
   cout << "-- Parameters list for photon emission:" << endl;
   cout << "----------------------------------------" << endl;
   cout << "tau_start =" << tau_start << " fm/c." << endl;
   cout << "T_dec = " << T_dec << " GeV." << endl;
   cout << "T_sw = " << T_sw_low <<  " to " << T_sw_high << " GeV."<< endl;
   cout << "deltaf_alpha = " << deltaf_alpha << endl;
   cout << endl;
   cout << "Photon momentum: " << photon_q_i << " to " << photon_q_f << " GeV, " << "n_q =" << np << endl;
   cout << "Photon momentum angles: " << photon_phi_q_i << " to " << photon_phi_q_f << ", n_phi=" << nphi << endl;
   cout << "Photon momentum rapidity: " << photon_y_i << " to " << photon_y_f << ", n_y =" << nrapidity << endl;

   cout << "******************************************" << endl;
   return;
}

void PhotonEmission::InitializePhotonEmissionRateTables()
{
   double photonrate_tb_Emin = PhotonemRatetableInfo_Emin;
   double photonrate_tb_Tmin = PhotonemRatetableInfo_Tmin;
   double photonrate_tb_dE = PhotonemRatetableInfo_dE;
   double photonrate_tb_dT = PhotonemRatetableInfo_dT;
   
   photon_QGP.setRateTableVarXmin(photonrate_tb_Tmin);
   photon_QGP.setRateTableVarYmin(photonrate_tb_Emin);
   photon_QGP.setRateTabledvarX(photonrate_tb_dT);
   photon_QGP.setRateTabledvarY(photonrate_tb_dE);
   photon_QGP.readEmissionrate("QGP_2to2_total");

   photon_pirho.setRateTableVarXmin(photonrate_tb_Tmin);
   photon_pirho.setRateTableVarYmin(photonrate_tb_Emin);
   photon_pirho.setRateTabledvarX(photonrate_tb_dT);
   photon_pirho.setRateTabledvarY(photonrate_tb_dE);
   photon_pirho.readEmissionrate("pion_rho_to_pion_gamma");

   photon_KstarK.setRateTableVarXmin(photonrate_tb_Tmin);
   photon_KstarK.setRateTableVarYmin(photonrate_tb_Emin);
   photon_KstarK.setRateTabledvarX(photonrate_tb_dT);
   photon_KstarK.setRateTabledvarY(photonrate_tb_dE);
   photon_KstarK.readEmissionrate("K_Kstar_to_pion_gamma");
   
   photon_piK.setRateTableVarXmin(photonrate_tb_Tmin);
   photon_piK.setRateTableVarYmin(photonrate_tb_Emin);
   photon_piK.setRateTabledvarX(photonrate_tb_dT);
   photon_piK.setRateTabledvarY(photonrate_tb_dE);
   photon_piK.readEmissionrate("pion_K_to_Kstar_gamma");

   photon_piKstar.setRateTableVarXmin(photonrate_tb_Tmin);
   photon_piKstar.setRateTableVarYmin(photonrate_tb_Emin);
   photon_piKstar.setRateTabledvarX(photonrate_tb_dT);
   photon_piKstar.setRateTabledvarY(photonrate_tb_dE);
   photon_piKstar.readEmissionrate("pion_Kstar_to_K_gamma");
  
   photon_pipi.setRateTableVarXmin(photonrate_tb_Tmin);
   photon_pipi.setRateTableVarYmin(photonrate_tb_Emin);
   photon_pipi.setRateTabledvarX(photonrate_tb_dT);
   photon_pipi.setRateTabledvarY(photonrate_tb_dE);
   photon_pipi.readEmissionrate("pion_pion_to_rho_gamma");

   photon_rhoK.setRateTableVarXmin(photonrate_tb_Tmin);
   photon_rhoK.setRateTableVarYmin(photonrate_tb_Emin);
   photon_rhoK.setRateTabledvarX(photonrate_tb_dT);
   photon_rhoK.setRateTabledvarY(photonrate_tb_dE);
   photon_rhoK.readEmissionrate("rho_K_to_K_gamma");

   photon_rho.setRateTableVarXmin(photonrate_tb_Tmin);
   photon_rho.setRateTableVarYmin(photonrate_tb_Emin);
   photon_rho.setRateTabledvarX(photonrate_tb_dT);
   photon_rho.setRateTabledvarY(photonrate_tb_dE);
   photon_rho.readEmissionrate("rho_to_pion_pion_gamma");
   
   photon_pirho_omegat.setRateTableVarXmin(photonrate_tb_Tmin);
   photon_pirho_omegat.setRateTableVarYmin(photonrate_tb_Emin);
   photon_pirho_omegat.setRateTabledvarX(photonrate_tb_dT);
   photon_pirho_omegat.setRateTabledvarY(photonrate_tb_dE);
   photon_pirho_omegat.readEmissionrate("pion_rho_to_omega_to_pion_gamma");
   
   return;
}

void PhotonEmission::calPhotonemission(readindata* frameptr, double* eta_ptr, double* volume)
{
  //photon momentum in the lab frame
  double p_q[np], phi_q[nphi], y_q[nrapidity];
  double sin_phiq[nphi], cos_phiq[nphi];
  double p_lab_local[4], p_lab_lowmu[4];
  double flow_u_mu_low[4];
  for(int k=0;k<nrapidity;k++)
     y_q[k] = photon_pirho.getPhotonrapidity(k);
  for(int l=0;l<np;l++) p_q[l] = photon_pirho.getPhotonp(l);
  for(int m=0;m<nphi;m++)
  {
     phi_q[m] = photon_pirho.getPhotonphi(m);
     sin_phiq[m] = sin(phi_q[m]);
     cos_phiq[m] = cos(phi_q[m]);
  }

  double e_local, p_local, temp_local, vx_local, vy_local;
  double eta_local;
  double** pi_tensor_lab = new double* [4];
  for(int i=0; i<4; i++)
    pi_tensor_lab[i] = new double [4];

  //loops over the transverse plane
  for(int i=0;i<nx;i++)
  {
    for(int j=0;j<ny;j++)
    {
      int idx_Tb = 0;
      temp_local = frameptr->temp[i][j];
      if(temp_local > T_dec)
      {
        e_local = frameptr->ed[i][j];
        p_local = frameptr->pl[i][j];
        vx_local = frameptr->vx[i][j];
        vy_local = frameptr->vy[i][j];
        pi_tensor_lab[0][0] = frameptr->pi00[i][j];
        pi_tensor_lab[0][1] = frameptr->pi01[i][j];
        pi_tensor_lab[0][2] = frameptr->pi02[i][j];
        pi_tensor_lab[0][3] = frameptr->pi03[i][j];
        pi_tensor_lab[1][0] = pi_tensor_lab[0][1];
        pi_tensor_lab[1][1] = frameptr->pi11[i][j];
        pi_tensor_lab[1][2] = frameptr->pi12[i][j];
        pi_tensor_lab[1][3] = frameptr->pi13[i][j];
        pi_tensor_lab[2][0] = pi_tensor_lab[0][2];
        pi_tensor_lab[2][1] = pi_tensor_lab[1][2];
        pi_tensor_lab[2][2] = frameptr->pi22[i][j];
        pi_tensor_lab[2][3] = frameptr->pi23[i][j];
        pi_tensor_lab[3][0] = pi_tensor_lab[0][3];
        pi_tensor_lab[3][1] = pi_tensor_lab[1][3];
        pi_tensor_lab[3][2] = pi_tensor_lab[2][3];
        pi_tensor_lab[3][3] = frameptr->pi33[i][j];

        getTransverseflow_u_mu_low(flow_u_mu_low, vx_local, vy_local);
        double prefactor_pimunu = 1./(2.*(e_local + p_local));
        for(int jj=0; jj<neta; jj++)
        {
          eta_local = eta_ptr[jj];
          
          //photon momentum loop
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
          photon_QGP.calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local, volume, QGP_fraction);
        }
        else if(temp_local > T_sw_low)
        {
          double QGP_fraction = (temp_local - T_sw_low)/(T_sw_high - T_sw_low);
          double HG_fraction = 1 - QGP_fraction;
          photon_QGP.calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local, volume, QGP_fraction);

          photon_pirho.calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
          photon_KstarK.calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
          photon_piK.calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
          photon_piKstar.calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
          photon_pipi.calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
          photon_rhoK.calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
          photon_rho.calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
          photon_pirho_omegat.calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
        }
        else
        {
          double HG_fraction = 1.0;
          photon_pirho.calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
          photon_KstarK.calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
          photon_piK.calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
          photon_piKstar.calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
          photon_pipi.calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
          photon_rhoK.calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
          photon_rho.calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
          photon_pirho_omegat.calThermalPhotonemission(Eq_localrest_Tb, pi_photon_Tb, idx_Tb, temp_local,  volume, HG_fraction);
        }
      }
    }
  }
  for(int i=0; i<4; i++)
  {
     delete [] pi_tensor_lab[i];
  }
  delete [] pi_tensor_lab;

  return;
}

void PhotonEmission::calPhoton_total_SpMatrix()
{
     cout << " Mark" << endl;
     for(int k=0;k<nrapidity;k++) 
     {
       for(int l=0;l<np;l++)
       { 
         for(int m=0;m<nphi;m++)
         {
           dNd2pTdphidy_hadrontot_eq[l][m][k] =  photon_pirho.getPhotonSpMatrix_eq(l, m, k) 
                         + photon_KstarK.getPhotonSpMatrix_eq(l, m, k)
                         + photon_piK.getPhotonSpMatrix_eq(l, m, k)
                         + photon_piKstar.getPhotonSpMatrix_eq(l, m, k)
                         + photon_pipi.getPhotonSpMatrix_eq(l, m, k)
                         + photon_rhoK.getPhotonSpMatrix_eq(l, m, k)
                         + photon_rho.getPhotonSpMatrix_eq(l, m, k);
                         + photon_pirho_omegat.getPhotonSpMatrix_eq(l, m, k);
           dNd2pTdphidy_hadrontot[l][m][k] =  photon_pirho.getPhotonSpMatrix_tot(l, m, k) 
                         + photon_KstarK.getPhotonSpMatrix_tot(l, m, k)
                         + photon_piK.getPhotonSpMatrix_tot(l, m, k)
                         + photon_piKstar.getPhotonSpMatrix_tot(l, m, k)
                         + photon_pipi.getPhotonSpMatrix_tot(l, m, k)
                         + photon_rhoK.getPhotonSpMatrix_tot(l, m, k)
                         + photon_rho.getPhotonSpMatrix_tot(l, m, k);
                         + photon_pirho_omegat.getPhotonSpMatrix_tot(l, m, k);

           dNd2pTdphidy_eq[l][m][k] =  photon_QGP.getPhotonSpMatrix_eq(l, m, k) + dNd2pTdphidy_hadrontot_eq[l][m][k];
           dNd2pTdphidy[l][m][k] =  photon_QGP.getPhotonSpMatrix_tot(l, m, k) + dNd2pTdphidy_hadrontot[l][m][k];

         }
       }
     }
     return;
}



void PhotonEmission::calPhoton_SpvnpT_individualchannel()
{
    photon_QGP.calPhoton_SpvnpT();

    photon_pirho.calPhoton_SpvnpT();
    photon_KstarK.calPhoton_SpvnpT();
    photon_piK.calPhoton_SpvnpT();
    photon_piKstar.calPhoton_SpvnpT();
    photon_pipi.calPhoton_SpvnpT();
    photon_rhoK.calPhoton_SpvnpT();
    photon_rho.calPhoton_SpvnpT();
    photon_pirho_omegat.calPhoton_SpvnpT();

    return;
}

void PhotonEmission::outputPhotonSpvn()
{
    photon_QGP.outputPhoton_SpvnpT();

    photon_pirho.outputPhoton_SpvnpT();
    photon_KstarK.outputPhoton_SpvnpT();
    photon_piK.outputPhoton_SpvnpT();
    photon_piKstar.outputPhoton_SpvnpT();
    photon_pipi.outputPhoton_SpvnpT();
    photon_rhoK.outputPhoton_SpvnpT();
    photon_rho.outputPhoton_SpvnpT();
    photon_pirho_omegat.outputPhoton_SpvnpT();

    outputPhoton_total_SpvnpT("photon_total");
    return;
}

void PhotonEmission::calPhoton_total_Spvn()
{
   int k = 0;
   for(int i=0;i<np;i++)
   {
       for(int j=0;j<nphi;j++)
       {
         double phi = photon_pirho.getPhotonphi(j);
         double phiweight = photon_pirho.getPhoton_phiweight(j);
         dNd2pT_hadrontot_eq[i] += dNd2pTdphidy_hadrontot_eq[i][j][k]*phiweight;
         dNd2pT_hadrontot[i] += dNd2pTdphidy_hadrontot[i][j][k]*phiweight;
         dNd2pT_eq[i] += dNd2pTdphidy_eq[i][j][k]*phiweight;
         dNd2pT[i] += dNd2pTdphidy[i][j][k]*phiweight;
         for(int order=1; order<norder; order++)
         {
            vnpT_hadrontot_cos_eq[order][i] += dNd2pTdphidy_hadrontot_eq[i][j][k]*cos(order*phi)*phiweight;
            vnpT_hadrontot_cos[order][i] += dNd2pTdphidy_hadrontot[i][j][k]*cos(order*phi)*phiweight;
            vnpT_cos_eq[order][i] += dNd2pTdphidy_eq[i][j][k]*cos(order*phi)*phiweight;
            vnpT_cos[order][i] += dNd2pTdphidy[i][j][k]*cos(order*phi)*phiweight;
            vnpT_hadrontot_sin_eq[order][i] += dNd2pTdphidy_hadrontot_eq[i][j][k]*sin(order*phi)*phiweight;
            vnpT_hadrontot_sin[order][i] += dNd2pTdphidy_hadrontot[i][j][k]*sin(order*phi)*phiweight;
            vnpT_sin_eq[order][i] += dNd2pTdphidy_eq[i][j][k]*sin(order*phi)*phiweight;
            vnpT_sin[order][i] += dNd2pTdphidy[i][j][k]*sin(order*phi)*phiweight;
         }
       }
       for(int order=1; order<norder; order++)
       {
          vnpT_hadrontot_cos_eq[order][i] = vnpT_hadrontot_cos_eq[order][i]/dNd2pT_hadrontot_eq[i];
          vnpT_hadrontot_cos[order][i] = vnpT_hadrontot_cos[order][i]/dNd2pT_hadrontot[i];
          vnpT_cos_eq[order][i] = vnpT_cos_eq[order][i]/dNd2pT_eq[i];
          vnpT_cos[order][i] = vnpT_cos[order][i]/dNd2pT[i];
          vnpT_hadrontot_sin_eq[order][i] = vnpT_hadrontot_sin_eq[order][i]/dNd2pT_hadrontot_eq[i];
          vnpT_hadrontot_sin[order][i] = vnpT_hadrontot_sin[order][i]/dNd2pT_hadrontot[i];
          vnpT_sin_eq[order][i] = vnpT_sin_eq[order][i]/dNd2pT_eq[i];
          vnpT_sin[order][i] = vnpT_sin[order][i]/dNd2pT[i];
       }
       
       dNd2pT_hadrontot_eq[i] = dNd2pT_hadrontot_eq[i]/(2*M_PI);
       dNd2pT_hadrontot[i] = dNd2pT_hadrontot[i]/(2*M_PI);
       dNd2pT_eq[i] = dNd2pT_eq[i]/(2*M_PI);
       dNd2pT[i] = dNd2pT[i]/(2*M_PI);
   }
   return;
}

void PhotonEmission::outputPhoton_total_SpvnpT(string filename)
{
    ostringstream filename_stream_Hadrontot_eq_SpMatrix;
    ostringstream filename_stream_Hadrontot_eq_Spvn;
    ostringstream filename_stream_Hadrontot_SpMatrix;
    ostringstream filename_stream_Hadrontot_Spvn;
    ostringstream filename_stream_eq_SpMatrix;
    ostringstream filename_stream_eq_Spvn;
    ostringstream filename_stream_SpMatrix;
    ostringstream filename_stream_Spvn;
    filename_stream_Hadrontot_eq_SpMatrix << filename << "_Hadrontot_eq_SpMatrix.dat";
    filename_stream_Hadrontot_eq_Spvn << filename << "_Hadrontot_eq_Spvn.dat";
    filename_stream_Hadrontot_SpMatrix << filename << "_Hadrontot_SpMatrix.dat";
    filename_stream_Hadrontot_Spvn << filename << "_Hadrontot_Spvn.dat";
    filename_stream_eq_SpMatrix << filename << "_eq_SpMatrix.dat";
    filename_stream_eq_Spvn << filename << "_eq_Spvn.dat";
    filename_stream_SpMatrix << filename << "_SpMatrix.dat";
    filename_stream_Spvn << filename << "_Spvn.dat";

    ofstream fphoton_Hdtot_eq_SpMatrix(filename_stream_Hadrontot_eq_SpMatrix.str().c_str());
    ofstream fphoton_Hdtot_eq_Spvn(filename_stream_Hadrontot_eq_Spvn.str().c_str());
    ofstream fphoton_Hdtot_SpMatrix(filename_stream_Hadrontot_SpMatrix.str().c_str());
    ofstream fphoton_Hdtot_Spvn(filename_stream_Hadrontot_Spvn.str().c_str());
    ofstream fphoton_eq_SpMatrix(filename_stream_eq_SpMatrix.str().c_str());
    ofstream fphoton_eq_Spvn(filename_stream_eq_Spvn.str().c_str());
    ofstream fphotonSpMatrix(filename_stream_SpMatrix.str().c_str());
    ofstream fphotonSpvn(filename_stream_Spvn.str().c_str());

    for(int i=0;i<nphi;i++)
    {
      double phi = photon_pirho.getPhotonphi(i);
      fphoton_eq_SpMatrix << phi << "  ";
      fphoton_Hdtot_eq_SpMatrix << phi << "  ";
      fphotonSpMatrix << phi << "  ";
      fphoton_Hdtot_SpMatrix << phi << "  ";
      for(int j=0;j<np;j++)
        for(int k=0;k<nrapidity;k++)
        {
           fphoton_Hdtot_eq_SpMatrix << scientific << setprecision(6) << setw(16) 
                           << dNd2pTdphidy_hadrontot_eq[j][i][k] << "  ";
           fphoton_eq_SpMatrix << scientific << setprecision(6) << setw(16) 
                           << dNd2pTdphidy_eq[j][i][k] << "  ";
           fphoton_Hdtot_SpMatrix << scientific << setprecision(6) << setw(16) 
                           << dNd2pTdphidy_hadrontot[j][i][k] << "  ";
           fphotonSpMatrix << scientific << setprecision(6) << setw(16) 
                           << dNd2pTdphidy[j][i][k] << "  ";
        }
      fphoton_Hdtot_eq_SpMatrix << endl;
      fphoton_eq_SpMatrix << endl;
      fphoton_Hdtot_SpMatrix << endl;
      fphotonSpMatrix << endl;
    }
    for(int i=0;i<np;i++)
    {
      double pT = photon_pirho.getPhotonp(i);
      fphoton_Hdtot_eq_Spvn << scientific << setprecision(6) << setw(16) 
                  << pT << "  " << dNd2pT_hadrontot_eq[i] << "  " ;
      fphoton_Hdtot_Spvn << scientific << setprecision(6) << setw(16) 
                  << pT << "  " << dNd2pT_hadrontot[i] << "  " ;
      fphoton_eq_Spvn << scientific << setprecision(6) << setw(16) 
                  << pT << "  " << dNd2pT_eq[i] << "  " ;
      fphotonSpvn << scientific << setprecision(6) << setw(16) 
                  << pT << "  " << dNd2pT[i] << "  " ;
      for(int order=1; order<norder; order++)
      {
         fphoton_Hdtot_eq_Spvn << scientific << setprecision(6) << setw(16) 
                     << vnpT_hadrontot_cos_eq[order][i] << "  "
                     << vnpT_hadrontot_sin_eq[order][i] << "  ";
         fphoton_Hdtot_Spvn << scientific << setprecision(6) << setw(16) 
                     << vnpT_hadrontot_cos[order][i] << "  "
                     << vnpT_hadrontot_sin[order][i] << "  ";
         fphoton_eq_Spvn << scientific << setprecision(6) << setw(16) 
                     << vnpT_cos_eq[order][i] << "  "
                     << vnpT_sin_eq[order][i] << "  ";
         fphotonSpvn << scientific << setprecision(6) << setw(16) 
                     << vnpT_cos[order][i] << "  "
                     << vnpT_sin[order][i] << "  ";
      }
      fphoton_Hdtot_eq_Spvn << endl;
      fphoton_Hdtot_Spvn << endl;
      fphoton_eq_Spvn << endl;
      fphotonSpvn << endl;
    }
    
    return;
}
