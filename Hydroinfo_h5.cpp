#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>

#include "hdf5.h"
#include "Arsenal.h"
#include "Hydroinfo_h5.h"
#include "ParameterReader.h"

using namespace std;

HydroinfoH5::HydroinfoH5(string filename, ParameterReader* paraRdr_in)
{
   paraRdr = paraRdr_in;

   // flag to determine whether to read evolutions for viscous variables
   Visflag = paraRdr->getVal("HydroinfoVisflag");

   herr_t status;
   const char *fileptr = (char*) filename.c_str();
   H5file_id = H5Fopen(fileptr, H5F_ACC_RDWR, H5P_DEFAULT);
   H5groupEventid = H5Gopen(H5file_id, "/Event", H5P_DEFAULT);

   readHydrogridInfo();
   printHydrogridInfo();

   Buffersize = paraRdr->getVal("HydroinfoBuffersize");
   dimensionX = grid_XH - grid_XL + 1;
   dimensionY = grid_YH - grid_YL + 1;

   if(Buffersize < grid_Framenum)
   {
      cout << "Buffersize is too small, increase it to at lease to " << grid_Framenum << endl;
      exit(1);
   }
   //initialize all matrices
   ed = new double** [Buffersize];
   sd = new double** [Buffersize];
   vx = new double** [Buffersize];
   vy = new double** [Buffersize];
   Temperature = new double** [Buffersize];
   Pressure = new double** [Buffersize];
   pi00 = new double** [Buffersize];
   pi01 = new double** [Buffersize];
   pi02 = new double** [Buffersize];
   pi03 = new double** [Buffersize];
   pi11 = new double** [Buffersize];
   pi12 = new double** [Buffersize];
   pi13 = new double** [Buffersize];
   pi22 = new double** [Buffersize];
   pi23 = new double** [Buffersize];
   pi33 = new double** [Buffersize];
   BulkPi = new double** [Buffersize];
   for(int i=0; i<Buffersize; i++)
   {
      ed[i] = new double* [dimensionX];
      sd[i] = new double* [dimensionX];
      vx[i] = new double* [dimensionX];
      vy[i] = new double* [dimensionX];
      Temperature[i] = new double* [dimensionX];
      Pressure[i] = new double* [dimensionX];
      pi00[i] = new double* [dimensionX];
      pi01[i] = new double* [dimensionX];
      pi02[i] = new double* [dimensionX];
      pi03[i] = new double* [dimensionX];
      pi11[i] = new double* [dimensionX];
      pi12[i] = new double* [dimensionX];
      pi13[i] = new double* [dimensionX];
      pi22[i] = new double* [dimensionX];
      pi23[i] = new double* [dimensionX];
      pi33[i] = new double* [dimensionX];
      BulkPi[i] = new double* [dimensionX];
      for(int j=0; j<dimensionX; j++)
      {
         ed[i][j] = new double [dimensionY];
         sd[i][j] = new double [dimensionY];
         vx[i][j] = new double [dimensionY];
         vy[i][j] = new double [dimensionY];
         Temperature[i][j] = new double [dimensionY];
         Pressure[i][j] = new double [dimensionY];
         pi00[i][j] = new double [dimensionY];
         pi01[i][j] = new double [dimensionY];
         pi02[i][j] = new double [dimensionY];
         pi03[i][j] = new double [dimensionY];
         pi11[i][j] = new double [dimensionY];
         pi12[i][j] = new double [dimensionY];
         pi13[i][j] = new double [dimensionY];
         pi22[i][j] = new double [dimensionY];
         pi23[i][j] = new double [dimensionY];
         pi33[i][j] = new double [dimensionY];
         BulkPi[i][j] = new double [dimensionY];
      }
   }
  
   readHydroinfoBuffered_total(); 

   status = H5Gclose(H5groupEventid);
   status = H5Fclose(H5file_id);
}

HydroinfoH5::~HydroinfoH5()
{
   for(int i=0; i<Buffersize; i++)
   {
      for(int j=0; j<dimensionX; j++)
      {
         delete[] ed[i][j];
         delete[] sd[i][j];
         delete[] vx[i][j];
         delete[] vy[i][j];
         delete[] Temperature[i][j];
         delete[] Pressure[i][j];
         delete[] pi00[i][j];
         delete[] pi01[i][j];
         delete[] pi02[i][j];
         delete[] pi03[i][j];
         delete[] pi11[i][j];
         delete[] pi12[i][j];
         delete[] pi13[i][j];
         delete[] pi22[i][j];
         delete[] pi23[i][j];
         delete[] pi33[i][j];
         delete[] BulkPi[i][j];
      }
      delete[] ed[i];
      delete[] sd[i];
      delete[] vx[i];
      delete[] vy[i];
      delete[] Temperature[i];
      delete[] Pressure[i];
      delete[] pi00[i];
      delete[] pi01[i];
      delete[] pi02[i];
      delete[] pi03[i];
      delete[] pi11[i];
      delete[] pi12[i];
      delete[] pi13[i];
      delete[] pi22[i];
      delete[] pi23[i];
      delete[] pi33[i];
      delete[] BulkPi[i];
   }
   delete[] ed;
   delete[] sd;
   delete[] vx;
   delete[] vy;
   delete[] Temperature;
   delete[] Pressure;
   delete[] pi00;
   delete[] pi01;
   delete[] pi02;
   delete[] pi03;
   delete[] pi11;
   delete[] pi12;
   delete[] pi13;
   delete[] pi22;
   delete[] pi23;
   delete[] pi33;
   delete[] BulkPi;
}

void HydroinfoH5::readHydrogridInfo()
{
   herr_t status;

   grid_XL = readH5Attribute_int(H5groupEventid, "XL");
   grid_XH = readH5Attribute_int(H5groupEventid, "XH");
   grid_YL = readH5Attribute_int(H5groupEventid, "YL");
   grid_YH = readH5Attribute_int(H5groupEventid, "YH");
   grid_Tau0 = readH5Attribute_double(H5groupEventid, "Tau0");
   grid_dTau = readH5Attribute_double(H5groupEventid, "dTau");
   grid_dx = readH5Attribute_double(H5groupEventid, "DX");
   grid_dy = readH5Attribute_double(H5groupEventid, "DY");

   grid_X0 = grid_XL * grid_dx;
   grid_Y0 = grid_YL * grid_dy;
   grid_Xmax = grid_XH * grid_dx;
   grid_Ymax = grid_YH * grid_dy;
   
   hsize_t tempFramenum;
   status = H5Gget_num_objs(H5groupEventid, &tempFramenum);
   grid_Framenum = (int) tempFramenum;
   grid_Taumax = grid_Tau0 + (grid_Framenum - 1)*grid_dTau;

   int tempflag = readH5Attribute_int(H5groupEventid, "OutputViscousFlag");
   Visflag = Visflag*tempflag;
}

void HydroinfoH5::printHydrogridInfo()
{
   cout << "-----------------------------------------" << endl;
   cout << "-----------hydro grid info---------------" << endl;
   cout << "-----------------------------------------" << endl;
   cout << "XL = " << grid_XL << endl;
   cout << "XH = " << grid_XH << endl;
   cout << "DX = " << grid_dx << " fm" << endl;
   cout << "YL = " << grid_YL << endl;
   cout << "YH = " << grid_YH << endl;
   cout << "DY = " << grid_dy << " fm" << endl;
   cout << "Tau0 = " << grid_Tau0 << " fm" << endl;
   cout << "dTau = " << grid_dTau << " fm" << endl;
   cout << "Number of Frames: " << grid_Framenum << endl;
   cout << "Taumax = " << grid_Taumax << " fm" << endl;
   cout << "Read in viscous information? ";
   if(Visflag == 1)
      cout << " Yes!" << endl;
   else
      cout << " No!" << endl;
   cout << "-----------------------------------------" << endl;
}

int HydroinfoH5::readH5Attribute_int(hid_t id, string attributeName)
{
   int attributeValue;
   hid_t attr;
   herr_t ret;
   attr = H5Aopen_name(id, attributeName.c_str());
   ret  = H5Aread(attr, H5T_NATIVE_INT, &attributeValue);
   ret =  H5Aclose(attr);
   return(attributeValue);
}

double HydroinfoH5::readH5Attribute_double(hid_t id, string attributeName)
{
   double attributeValue;
   hid_t attr;
   herr_t ret;
   attr = H5Aopen_name(id, attributeName.c_str());
   ret  = H5Aread(attr, H5T_NATIVE_DOUBLE, &attributeValue);
   ret =  H5Aclose(attr);
   return(attributeValue);
}

void HydroinfoH5::readHydroinfoBuffered_total()
{
   hid_t group_id;
   herr_t status;
   
   int frameIdx;
   for(int i=0; i<Buffersize; i++)
   {
      frameIdx = i;
      if(frameIdx < (int) grid_Framenum)
      {
         stringstream frameName;
         frameName << "Frame_" <<  setw(4) << setfill('0') << frameIdx;
         group_id = H5Gopen(H5groupEventid, frameName.str().c_str(), H5P_DEFAULT);
      
         readH5Dataset_double(group_id, "e", ed[i]);
         readH5Dataset_double(group_id, "s", sd[i]);
         readH5Dataset_double(group_id, "Vx", vx[i]);
         readH5Dataset_double(group_id, "Vy", vy[i]);
         readH5Dataset_double(group_id, "Temp", Temperature[i]);
         readH5Dataset_double(group_id, "P", Pressure[i]);
         if(Visflag == 1)
         {
            readH5Dataset_double(group_id, "Pi00", pi00[i]);
            readH5Dataset_double(group_id, "Pi01", pi01[i]);
            readH5Dataset_double(group_id, "Pi02", pi02[i]);
            readH5Dataset_double(group_id, "Pi03", pi03[i]);
            readH5Dataset_double(group_id, "Pi11", pi11[i]);
            readH5Dataset_double(group_id, "Pi12", pi12[i]);
            readH5Dataset_double(group_id, "Pi13", pi13[i]);
            readH5Dataset_double(group_id, "Pi22", pi22[i]);
            readH5Dataset_double(group_id, "Pi23", pi23[i]);
            readH5Dataset_double(group_id, "Pi33", pi33[i]);
            readH5Dataset_double(group_id, "BulkPi", BulkPi[i]);
         }
         status = H5Gclose(group_id);
      }
      else
         break;
   }
}

void HydroinfoH5::readHydroinfoSingleframe(int frameIdx)
{
   hid_t group_id;
   herr_t status;

   if(frameIdx < (int) grid_Framenum)
   {
      stringstream frameName;
      frameName << "Frame_" <<  setw(4) << setfill('0') << frameIdx;
      group_id = H5Gopen(H5groupEventid, frameName.str().c_str(), H5P_DEFAULT);
      
      int Idx = frameIdx;

      readH5Dataset_double(group_id, "e", ed[Idx]);
      readH5Dataset_double(group_id, "s", sd[Idx]);
      readH5Dataset_double(group_id, "Vx", vx[Idx]);
      readH5Dataset_double(group_id, "Vy", vy[Idx]);
      readH5Dataset_double(group_id, "Temp", Temperature[Idx]);
      readH5Dataset_double(group_id, "P", Pressure[Idx]);
      if(Visflag == 1)
      {
         readH5Dataset_double(group_id, "Pi00", pi00[Idx]);
         readH5Dataset_double(group_id, "Pi01", pi01[Idx]);
         readH5Dataset_double(group_id, "Pi02", pi02[Idx]);
         readH5Dataset_double(group_id, "Pi03", pi03[Idx]);
         readH5Dataset_double(group_id, "Pi11", pi11[Idx]);
         readH5Dataset_double(group_id, "Pi12", pi12[Idx]);
         readH5Dataset_double(group_id, "Pi13", pi13[Idx]);
         readH5Dataset_double(group_id, "Pi22", pi22[Idx]);
         readH5Dataset_double(group_id, "Pi23", pi23[Idx]);
         readH5Dataset_double(group_id, "Pi33", pi33[Idx]);
         readH5Dataset_double(group_id, "BulkPi", BulkPi[Idx]);
      }
      status = H5Gclose(group_id);
   }
   else
   {
      cout << "Error: readHydroinfoSingleframe :: frameIdx exceed maximum frame number from hydro" << endl;
      cout << "frameIdx = " << frameIdx << endl;
      exit(1);
   }
}

void HydroinfoH5::readH5Dataset_double(hid_t id, string datasetName, double** dset_data)
{
   herr_t status;
   hid_t dataset_id;
   
   double temp_data[dimensionX][dimensionY];
   dataset_id = H5Dopen(id, datasetName.c_str(), H5P_DEFAULT);
   status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_data);
   for(int i=0; i<dimensionX; i++)
      for(int j=0; j<dimensionY; j++)
         dset_data[i][j] = temp_data[i][j];
   status = H5Dclose(dataset_id);
}

void HydroinfoH5::getHydroinfoOnlattice(int frameIdx, int xIdx, int yIdx, fluidCell* fluidCellptr)
{
   if(frameIdx < 0 || frameIdx > grid_Framenum || xIdx < 0 || xIdx > (grid_XH - grid_XL) || yIdx < 0 || yIdx > (grid_YH - grid_YL))
   {
      cout << "Error: getHydroinfoOnlattice:: Index is wrong" << endl;
      cout << "frameIdx = " << frameIdx << " xIdx = " << xIdx 
           << "yIdx = " << yIdx << endl;
      exit(1);
   }
   fluidCellptr->ed = ed[frameIdx][xIdx][yIdx];
   fluidCellptr->sd = sd[frameIdx][xIdx][yIdx];
   fluidCellptr->vx = vx[frameIdx][xIdx][yIdx];
   fluidCellptr->vy = vy[frameIdx][xIdx][yIdx];
   fluidCellptr->temperature = Temperature[frameIdx][xIdx][yIdx];
   fluidCellptr->pressure = Pressure[frameIdx][xIdx][yIdx];
   fluidCellptr->pi[0][0] = pi00[frameIdx][xIdx][yIdx];
   fluidCellptr->pi[0][1] = pi01[frameIdx][xIdx][yIdx];
   fluidCellptr->pi[0][2] = pi02[frameIdx][xIdx][yIdx];
   fluidCellptr->pi[0][3] = pi03[frameIdx][xIdx][yIdx];
   fluidCellptr->pi[1][0] = fluidCellptr->pi[0][1];
   fluidCellptr->pi[1][1] = pi11[frameIdx][xIdx][yIdx];
   fluidCellptr->pi[1][2] = pi12[frameIdx][xIdx][yIdx];
   fluidCellptr->pi[1][3] = pi13[frameIdx][xIdx][yIdx];
   fluidCellptr->pi[2][0] = fluidCellptr->pi[0][2];
   fluidCellptr->pi[2][1] = fluidCellptr->pi[1][2];
   fluidCellptr->pi[2][2] = pi22[frameIdx][xIdx][yIdx];
   fluidCellptr->pi[2][3] = pi23[frameIdx][xIdx][yIdx];
   fluidCellptr->pi[3][0] = fluidCellptr->pi[0][3];
   fluidCellptr->pi[3][1] = fluidCellptr->pi[1][3];
   fluidCellptr->pi[3][2] = fluidCellptr->pi[2][3];
   fluidCellptr->pi[3][3] = pi33[frameIdx][xIdx][yIdx];
   fluidCellptr->bulkPi = BulkPi[frameIdx][xIdx][yIdx];
}


void HydroinfoH5::getHydroinfo(double tau, double x, double y, fluidCell* fluidCellptr)
{
   if(tau < grid_Tau0 || tau >= grid_Taumax || x < grid_X0 || x >= grid_Xmax || y < grid_Y0 || y >= grid_Ymax)
   {
      setZero_fluidCell(fluidCellptr);
      return;
   }
   int frameIdx, xIdx, yIdx;
   double tauInc, xInc, yInc;
   double temp;

   temp = (tau - grid_Tau0)/grid_dTau;
   frameIdx = (int) floor(temp);
   tauInc = temp - frameIdx;
   
   temp = (x - grid_X0)/grid_dx;
   xIdx = (int) floor(temp);
   xInc = temp - xIdx;

   temp = (y - grid_Y0)/grid_dy;
   yIdx = (int) floor(temp);
   yInc = temp - yIdx;

   fluidCellptr->ed = cubeInterp(tauInc, xInc, yInc, 
      ed[frameIdx][xIdx][yIdx], ed[frameIdx+1][xIdx][yIdx], ed[frameIdx][xIdx+1][yIdx], ed[frameIdx+1][xIdx+1][yIdx], 
      ed[frameIdx][xIdx][yIdx+1], ed[frameIdx+1][xIdx][yIdx+1], ed[frameIdx][xIdx+1][yIdx+1], ed[frameIdx+1][xIdx+1][yIdx+1]);
   fluidCellptr->sd = cubeInterp(tauInc, xInc, yInc, 
      sd[frameIdx][xIdx][yIdx], sd[frameIdx+1][xIdx][yIdx], sd[frameIdx][xIdx+1][yIdx], sd[frameIdx+1][xIdx+1][yIdx], 
      sd[frameIdx][xIdx][yIdx+1], sd[frameIdx+1][xIdx][yIdx+1], sd[frameIdx][xIdx+1][yIdx+1], sd[frameIdx+1][xIdx+1][yIdx+1]);
   fluidCellptr->vx = cubeInterp(tauInc, xInc, yInc, 
      vx[frameIdx][xIdx][yIdx], vx[frameIdx+1][xIdx][yIdx], vx[frameIdx][xIdx+1][yIdx], vx[frameIdx+1][xIdx+1][yIdx], 
      vx[frameIdx][xIdx][yIdx+1], vx[frameIdx+1][xIdx][yIdx+1], vx[frameIdx][xIdx+1][yIdx+1], vx[frameIdx+1][xIdx+1][yIdx+1]);
   fluidCellptr->vy = cubeInterp(tauInc, xInc, yInc, 
      vy[frameIdx][xIdx][yIdx], vy[frameIdx+1][xIdx][yIdx], vy[frameIdx][xIdx+1][yIdx], vy[frameIdx+1][xIdx+1][yIdx], 
      vy[frameIdx][xIdx][yIdx+1], vy[frameIdx+1][xIdx][yIdx+1], vy[frameIdx][xIdx+1][yIdx+1], vy[frameIdx+1][xIdx+1][yIdx+1]);
   fluidCellptr->temperature = cubeInterp(tauInc, xInc, yInc, 
      Temperature[frameIdx][xIdx][yIdx], Temperature[frameIdx+1][xIdx][yIdx], Temperature[frameIdx][xIdx+1][yIdx], Temperature[frameIdx+1][xIdx+1][yIdx], 
      Temperature[frameIdx][xIdx][yIdx+1], Temperature[frameIdx+1][xIdx][yIdx+1], Temperature[frameIdx][xIdx+1][yIdx+1], Temperature[frameIdx+1][xIdx+1][yIdx+1]);
   fluidCellptr->pressure = cubeInterp(tauInc, xInc, yInc, 
      Pressure[frameIdx][xIdx][yIdx], Pressure[frameIdx+1][xIdx][yIdx], Pressure[frameIdx][xIdx+1][yIdx], Pressure[frameIdx+1][xIdx+1][yIdx], 
      Pressure[frameIdx][xIdx][yIdx+1], Pressure[frameIdx+1][xIdx][yIdx+1], Pressure[frameIdx][xIdx+1][yIdx+1], Pressure[frameIdx+1][xIdx+1][yIdx+1]);
   if(Visflag == 1)
   {
      fluidCellptr->pi[0][0] = cubeInterp(tauInc, xInc, yInc, 
         pi00[frameIdx][xIdx][yIdx], pi00[frameIdx+1][xIdx][yIdx], pi00[frameIdx][xIdx+1][yIdx], pi00[frameIdx+1][xIdx+1][yIdx], 
         pi00[frameIdx][xIdx][yIdx+1], pi00[frameIdx+1][xIdx][yIdx+1], pi00[frameIdx][xIdx+1][yIdx+1], pi00[frameIdx+1][xIdx+1][yIdx+1]);
      fluidCellptr->pi[0][1] = cubeInterp(tauInc, xInc, yInc, 
         pi01[frameIdx][xIdx][yIdx], pi01[frameIdx+1][xIdx][yIdx], pi01[frameIdx][xIdx+1][yIdx], pi01[frameIdx+1][xIdx+1][yIdx], 
         pi01[frameIdx][xIdx][yIdx+1], pi01[frameIdx+1][xIdx][yIdx+1], pi01[frameIdx][xIdx+1][yIdx+1], pi01[frameIdx+1][xIdx+1][yIdx+1]);
      fluidCellptr->pi[0][2] = cubeInterp(tauInc, xInc, yInc, 
         pi02[frameIdx][xIdx][yIdx], pi02[frameIdx+1][xIdx][yIdx], pi02[frameIdx][xIdx+1][yIdx], pi02[frameIdx+1][xIdx+1][yIdx], 
         pi02[frameIdx][xIdx][yIdx+1], pi02[frameIdx+1][xIdx][yIdx+1], pi02[frameIdx][xIdx+1][yIdx+1], pi02[frameIdx+1][xIdx+1][yIdx+1]);
      fluidCellptr->pi[0][3] = cubeInterp(tauInc, xInc, yInc, 
         pi03[frameIdx][xIdx][yIdx], pi03[frameIdx+1][xIdx][yIdx], pi03[frameIdx][xIdx+1][yIdx], pi03[frameIdx+1][xIdx+1][yIdx], 
         pi03[frameIdx][xIdx][yIdx+1], pi03[frameIdx+1][xIdx][yIdx+1], pi03[frameIdx][xIdx+1][yIdx+1], pi03[frameIdx+1][xIdx+1][yIdx+1]);
      fluidCellptr->pi[1][0] = fluidCellptr->pi[0][1];
      fluidCellptr->pi[1][1] = cubeInterp(tauInc, xInc, yInc, 
         pi11[frameIdx][xIdx][yIdx], pi11[frameIdx+1][xIdx][yIdx], pi11[frameIdx][xIdx+1][yIdx], pi11[frameIdx+1][xIdx+1][yIdx], 
         pi11[frameIdx][xIdx][yIdx+1], pi11[frameIdx+1][xIdx][yIdx+1], pi11[frameIdx][xIdx+1][yIdx+1], pi11[frameIdx+1][xIdx+1][yIdx+1]);
      fluidCellptr->pi[1][2] = cubeInterp(tauInc, xInc, yInc, 
         pi12[frameIdx][xIdx][yIdx], pi12[frameIdx+1][xIdx][yIdx], pi12[frameIdx][xIdx+1][yIdx], pi12[frameIdx+1][xIdx+1][yIdx], 
         pi12[frameIdx][xIdx][yIdx+1], pi12[frameIdx+1][xIdx][yIdx+1], pi12[frameIdx][xIdx+1][yIdx+1], pi12[frameIdx+1][xIdx+1][yIdx+1]);
      fluidCellptr->pi[1][3] = cubeInterp(tauInc, xInc, yInc, 
         pi13[frameIdx][xIdx][yIdx], pi13[frameIdx+1][xIdx][yIdx], pi13[frameIdx][xIdx+1][yIdx], pi13[frameIdx+1][xIdx+1][yIdx], 
         pi13[frameIdx][xIdx][yIdx+1], pi13[frameIdx+1][xIdx][yIdx+1], pi13[frameIdx][xIdx+1][yIdx+1], pi13[frameIdx+1][xIdx+1][yIdx+1]);
      fluidCellptr->pi[2][0] = fluidCellptr->pi[0][2];
      fluidCellptr->pi[2][1] = fluidCellptr->pi[1][2];
      fluidCellptr->pi[2][2] = cubeInterp(tauInc, xInc, yInc, 
         pi22[frameIdx][xIdx][yIdx], pi22[frameIdx+1][xIdx][yIdx], pi22[frameIdx][xIdx+1][yIdx], pi22[frameIdx+1][xIdx+1][yIdx], 
         pi22[frameIdx][xIdx][yIdx+1], pi22[frameIdx+1][xIdx][yIdx+1], pi22[frameIdx][xIdx+1][yIdx+1], pi22[frameIdx+1][xIdx+1][yIdx+1]);
      fluidCellptr->pi[2][3] = cubeInterp(tauInc, xInc, yInc, 
         pi23[frameIdx][xIdx][yIdx], pi23[frameIdx+1][xIdx][yIdx], pi23[frameIdx][xIdx+1][yIdx], pi23[frameIdx+1][xIdx+1][yIdx], 
         pi23[frameIdx][xIdx][yIdx+1], pi23[frameIdx+1][xIdx][yIdx+1], pi23[frameIdx][xIdx+1][yIdx+1], pi23[frameIdx+1][xIdx+1][yIdx+1]);
      fluidCellptr->pi[3][0] = fluidCellptr->pi[0][3];
      fluidCellptr->pi[3][1] = fluidCellptr->pi[1][3];
      fluidCellptr->pi[3][2] = fluidCellptr->pi[2][3];
      fluidCellptr->pi[3][3] = cubeInterp(tauInc, xInc, yInc, 
         pi33[frameIdx][xIdx][yIdx], pi33[frameIdx+1][xIdx][yIdx], pi33[frameIdx][xIdx+1][yIdx], pi33[frameIdx+1][xIdx+1][yIdx], 
         pi33[frameIdx][xIdx][yIdx+1], pi33[frameIdx+1][xIdx][yIdx+1], pi33[frameIdx][xIdx+1][yIdx+1], pi33[frameIdx+1][xIdx+1][yIdx+1]);
      fluidCellptr->bulkPi = cubeInterp(tauInc, xInc, yInc, 
         BulkPi[frameIdx][xIdx][yIdx], BulkPi[frameIdx+1][xIdx][yIdx], BulkPi[frameIdx][xIdx+1][yIdx], BulkPi[frameIdx+1][xIdx+1][yIdx], 
         BulkPi[frameIdx][xIdx][yIdx+1], BulkPi[frameIdx+1][xIdx][yIdx+1], BulkPi[frameIdx][xIdx+1][yIdx+1], BulkPi[frameIdx+1][xIdx+1][yIdx+1]);
   }
   else
   {
      for(int i=0; i<4; i++)
         for(int j=0; j<4; j++)
            fluidCellptr->pi[i][j] = 0.0e0;
      fluidCellptr->bulkPi = 0.0e0;
   }
}

void HydroinfoH5::setZero_fluidCell(fluidCell* fluidCellptr)
{
   fluidCellptr->ed = 0.0e0;
   fluidCellptr->sd = 0.0e0;
   fluidCellptr->vx = 0.0e0;
   fluidCellptr->vy = 0.0e0;
   fluidCellptr->temperature = 0.0e0;
   fluidCellptr->pressure = 0.0e0;
   for(int i=0; i<4; i++)
      for(int j=0; j<4; j++)
         fluidCellptr->pi[i][j] = 0.0e0;
   fluidCellptr->bulkPi = 0.0e0;
}
