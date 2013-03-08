#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip> 

#include "OSCARreader.h"
#include "parameter.h"

using namespace std;

OSCARreader::OSCARreader(string infilename)
{
   ostringstream filename_stream;
   filename_stream << infilename;
   fp = new ifstream ;
   fp->open(filename_stream.str().c_str());
   if (fp->is_open()==false)
   {
      cout << "OSCARreader::OSCARreader error: the data file cannot be opened." << endl;
      exit(-1);
   }

   readheader();
   printinfo();
}

OSCARreader::~OSCARreader()
{
   fp->close();
   delete fp;
   delete headerptr;
}

void OSCARreader::readheader()
{
   headerptr= new OSCARheader;
   // read the header
   char dummy[80];
   *(fp) >> headerptr->filename >> headerptr->hydrotype >> headerptr->datainf
         >> dummy >> headerptr->initial >>  dummy >> headerptr->nucleartype
         >> dummy >> headerptr->impactb >> headerptr->maxentropy 
         >> headerptr->thermaltime >> headerptr->etas >> dummy >> headerptr->EOS
         >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy 
         >> headerptr->geometry >> dummy >> dummy
         >> headerptr->timestep >> headerptr->n_x >> headerptr->n_y >> dummy 
         >> dummy >> dummy >> dummy >> headerptr->ti >> headerptr->tf
         >> headerptr->xi 
         >> headerptr->xf 
         >> headerptr->yi 
         >> headerptr->yf;

   for(int i=0;i<13;i++)
       *(fp) >> dummy;
}

void OSCARreader::printinfo()
{
   //print out the header information
   cout << "read in file name: " << headerptr->filename << endl;
   cout <<  headerptr->hydrotype << " " << headerptr->geometry 
        << " hydrodynamic simulation" << endl;
   cout <<  headerptr->initial << " initialization for " 
        << headerptr->nucleartype << " collision" << endl;
   cout << "impact parameter " << headerptr->impactb << " fm, " 
        << "cental entropy density " << headerptr->maxentropy << " fm^-3." 
        << endl;
   cout << "thermalization time " << headerptr->thermaltime << " fm/c, "
        << "specific shear viscosity " << headerptr->etas << endl;
   cout << "equation of state: " << headerptr->EOS << endl;
}

void OSCARreader::readframe_2Dboostinvariant(readindata* dataptr)
{
   int ix, iy;
   //read in data for one frame
   for(int i=0;i<nx;i++)
   {
      for(int j=0;j<ny;j++) 
      {
         *(fp) >> dataptr->itime >> ix >> iy >> dataptr->ed[i][j]
               >> dataptr->pl[i][j] >> dataptr->temp[i][j]
               >> dataptr->r_qgp[i][j]
               >> dataptr->vx[i][j] >> dataptr->vy[i][j]
               >> dataptr->pi00[i][j] >> dataptr->pi01[i][j]
               >> dataptr->pi02[i][j] >> dataptr->pi11[i][j] 
               >> dataptr->pi12[i][j] >> dataptr->pi22[i][j] 
               >> dataptr->pi33[i][j]
               >> dataptr->PPi[i][j] >> dataptr->eta_s[i][j] 
               >> dataptr->taupi[i][j] >> dataptr->tauPib[i][j];
         dataptr->pi03[i][j] = 0.0e0;
         dataptr->pi13[i][j] = 0.0e0;
         dataptr->pi23[i][j] = 0.0e0;
      }
   }
}

void OSCARreader::printframe(readindata* dataptr, string filename)
{
   ostringstream filename_stream;
   filename_stream << filename << ".dat";
   ofstream output(filename_stream.str().c_str());
   //read in data for one frame
   for(int i=0;i<nx;i++)
   {
      for(int j=0;j<ny;j++) 
      {
        output << scientific << setw(16) <<  setprecision(6) 
               << dataptr->itime << "  " <<  i << "  " <<  j << "  " <<  dataptr->ed[i][j]
               << "  " <<  dataptr->pl[i][j] << "  " <<  dataptr->temp[i][j]
               << "  " <<  dataptr->r_qgp[i][j]
               << "  " <<  dataptr->vx[i][j] << "  " <<  dataptr->vy[i][j]
               << "  " <<  dataptr->pi00[i][j] << "  " <<  dataptr->pi01[i][j]
               << "  " <<  dataptr->pi02[i][j] << "  " <<  dataptr->pi11[i][j] 
               << "  " <<  dataptr->pi12[i][j] << "  " <<  dataptr->pi22[i][j] 
               << "  " <<  dataptr->pi33[i][j]
               << "  " <<  dataptr->PPi[i][j] << "  " <<  dataptr->eta_s[i][j] 
               << "  " <<  dataptr->taupi[i][j] << "  " <<  dataptr->tauPib[i][j] 
               << endl;
      }
   }
}


void OSCARreader::readnextframe(readindata** oldpp, readindata** newpp)
{
   //read the next frame into the programe
   readframe_2Dboostinvariant(*(oldpp));
   readindata* tmpptr = *(newpp);
   *(newpp) = *(oldpp);
   *(oldpp) = tmpptr;
}

