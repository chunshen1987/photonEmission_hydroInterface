/////////////////////////////////////////////////////////////////////////
//                      hydrodynamics analysis
//                          photon emission 
//
//              author: Chun Shen <shen@mps.ohio-state.edu>
//
//  This program calculates the photon emission from the relativistic
//  heavy ion collision. It reads in viscous hydrodynamics results in 
//  OSCAR format.
//  
//
//  To do in the future:
/////////////////////////////////////////////////////////////////////////

#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>

#include "PhotonEmission.h"
#include "OSCARreader.h"
#include "Stopwatch.h"
#include "Arsenal.h"
#include "parameter.h"
#include "BjorkenExpansion.h"

using namespace std;

int main()
{
  Stopwatch sw;

  sw.tic();
  readindata* frameptr = new readindata ; //hydro data file pionter

  OSCARreader OSCARinputfile("OSCAR2008H.dat");

  //get grid information from OSCAR header
  double grid_t0 = OSCARinputfile.getGridt0();
  double grid_dt, grid_dx, grid_dy;
  grid_dt = OSCARinputfile.getGriddt();
  grid_dx = OSCARinputfile.getGriddx();
  grid_dy = OSCARinputfile.getGriddy();

  int ntime = OSCARinputfile.getTimestep() + (int)((grid_t0 - tau_start)/grid_dt);
  ntime = 2;

  double* eta_ptr = new double [neta];
  double* etaweight_ptr = new double [neta];
  double* tanheta_ptr = new double [neta];
  double* dvolume = new double [neta];
  gauss(neta, 0, eta_i, eta_f, eta_ptr, etaweight_ptr);
  for(int i=0; i<neta ; i++)
  {
     tanheta_ptr[i] = tanh(eta_ptr[i]);
     dvolume[i] = grid_dt*grid_dx*grid_dy*etaweight_ptr[i]; //dtau*dx*dy*deta
  }
  
  OSCARinputfile.readframe_2Dboostinvariant(frameptr);

  PhotonEmission testphoton;

  double tau_local;
  double* volume_element = new double [neta];

  double tau_1, tau_2;
  tau_1 = grid_t0;
  BjorkenExpansion BjorkenExpan;

  for(int itime=0;itime<ntime;itime++) //loop over time evolution
  {
    //OSCARinputfile.readframe(frameptr);
    //tau_local = grid_t0 + itime*grid_dt;
    tau_local = tau_start + itime*grid_dt;
    tau_2 = tau_local;

    if(tau_local <= grid_t0)
    {
       //backtrace_Temperature_Bjorken_1Dlongitudinalexpansion(frameptr, tau_1, tau_2);
       BjorkenExpan.backtrace_Entropy_Bjorken_1Dlongitudinalexpansion(frameptr, tau_1, tau_2);
       tau_1 = tau_2;
    }

    for(int k=0; k<neta; k++)
       volume_element[k] = 2 * tau_local * dvolume[k]; //volume element: tau*dtau*dx*dy*deta, 2 for symmetry along longitudinal direction

    testphoton.calPhotonemission(frameptr, eta_ptr, tanheta_ptr, volume_element);

    cout<<"frame "<< itime << " : ";
    cout<<" tau = " << setw(4) << setprecision(3) << tau_local 
        <<" fm/c done!" << endl;
    
    if(tau_local > grid_t0)
       OSCARinputfile.readframe_2Dboostinvariant(frameptr);

  }

  testphoton.calPhoton_SpvnpT_individualchannel();
  testphoton.calPhoton_total_SpMatrix();
  testphoton.calPhoton_total_Spvn();

  testphoton.outputPhotonSpvn();

  sw.toc();
  cout << "totally takes : " << sw.takeTime() << " seconds." << endl;

//test part

//  double test;
//  test=ph_em_rate_pirho(0.4,0.12,0.18,phrate_tbptr->pirho_mu_ptr);
//  test=interpolation2D_bilinear(phrate_tbptr->pirho_mu_ptr, 0.2,4);
//  cout<<test<< endl;
//  photon_emission(dNd2pTdphi, 0.12, 0.18, 0.2, 0.4, 0, volume, chempot_tbptr); 

/*
  double p[4], k[4], lambda[4][4];
  double vx,vy,vz;

  vx=0.2;
  vy=0.4;
  vz=0;

  p[0] = 1;
  p[1] = sin(pi/6);
  p[2] = cos(pi/6);
  p[3] = 0;

  boost_matrix(lambda,vx,vy,vz);
  boost_trans(p,k,lambda);
  boost_matrix(lambda,-vx,-vy,-vz);
  boost_trans(k,p,lambda);

  cout<< p[0] << "  " << p[1] << "  " << p[2] << "  " << p[3] << endl; 
  cout<< k[0] << "  " << k[1] << "  " << k[2] << "  " << k[3] << endl;
  for(int i=0;i<4;i++)
  {
    for(int j=0;j<4;j++)
       cout << lambda[i][j] << "  ";
    cout<<endl;
  }
*/
/*
  ofstream testresult("test.dat");
  double test1,test2,test3,test4;
  double ph_e0 = 0.2;
  double dph_e = 0.028;
  double ph_e = ph_e0;
  for(int i=0;i<=100;i++)
  {
    test1 = ph_em_rate_pipi(ph_e, 0.12, 0.18, chempot_tbptr);
    test2 = ph_em_rate_pirho(ph_e, 0.12, 0.18, chempot_tbptr);
    test3 = ph_em_rate_pipi_raw(ph_e, 0.12, 0.18, chempot_tbptr);
    test4 = ph_em_rate_pirho_raw(ph_e, 0.12, 0.18, chempot_tbptr);
    testresult << ph_e << "  " <<  test1 << "  " << test3 << "  " 
               << test2 << "  " << test4 << endl;
    ph_e += dph_e;
  }
*/
/*
  double emax = findemax(oldframeptr);  //find the maximum energy
  double einitial = 0.95*emax;          //set highest energy density surface to analysis
  double efreeze = 0.234;                //set the lowest energy density surface to analysis
  double de = pow(efreeze/einitial, 1./(npoint-1));   //set step length for energy density
  cout << "the maxinum energy density in the center is:" 
       << scientific << setprecision(8) << emax << endl;
*/
  return(0);
}

