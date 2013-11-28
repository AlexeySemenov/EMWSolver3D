#define _USE_MATH_DEFINES
#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include "netcdf.h"

#define PI 3.1415926535897932384626433832795028841971
#define CC 2.99792458e+8 //скорость света в вакууме
#define CCC 2.99792458e+10 //скорость света в вакууме
#define MU_Z 0.0000012566370614359173//постоянная магнитной пр-ти в вакууме  4.0*M_PI*1.0e-7 
#define EPS_Z 0.0000000000088541878176203892 //диэлектрическая проницаемость 1.0/(CC*CC*MU_Z)
#define FREQ 5.0e+9 // частота источника
#define LAMBDA (CC/FREQ) // чтоб лишние разы не делить при фикс частоте можно забить вручную
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}
const int gridX = 200;//1000;
const int gridY = 600;//300;

double ez[gridX][gridY];

double ex[gridX][gridY];
double ey[gridX][gridY];
double hz[gridX][gridY];

double expr[gridX][gridY];
double eypr[gridX][gridY];
double hzpr[gridX][gridY];

double xDim[gridX];
double yDim[gridY];
double* data;
using namespace std;
void readFile(int timestp, string path)
{   
 //   string fileName;
	//stringstream s;
	//s << "K:\\Ez\\Ez[" << timestp << "].dat";
	//fileName = s.str();
	// ifstream file(fileName.c_str(), ios::in | ios::binary);

	 string fileName2;
	stringstream s2;
	 s2 << path<<"\\Ez[" << timestp << "].dat";
	fileName2 = s2.str();

	 ifstream file2(fileName2.c_str(), ios::in | ios::binary);

	// string fileName3;
	//stringstream s3;
	// s3 << path<<"\\Hx[" << timestp << "].dat";
	//fileName3 = s3.str();
	// ifstream file3(fileName3.c_str(), ios::in | ios::binary);

	//string fileName4;
	//stringstream s4;
	// s4 << path<<"\\Hy[" << timestp << "].dat";
	//fileName4 = s4.str();
	// ifstream file4(fileName4.c_str(), ios::in | ios::binary);

	 double tmp1,tmp2, tmp3, tmp4;
	 //if(!file.is_open())
		//cout << "ERROR::Cannot open file for read! " << fileName << endl ;
	if(!file2.is_open())
		cout << "ERROR::Cannot open file for read! " << fileName2 << endl ;
	//if(!file3.is_open())
	//	cout << "ERROR::Cannot open file for read! " << fileName3 << endl ;
	//if(!file4.is_open())
	//	cout << "ERROR::Cannot open file for read! " << fileName4 << endl ;

	cout << "Reading file... " << timestp <<"...";
	int k = 0;
	for(int i = 0; i < gridX; i++)
		for(int j = 0; j < gridY; j++)
		{

			//file.read((char *) &tmp1, sizeof(double));
			file2.read((char *) &tmp2, sizeof(double));
			//file3.read((char *) &tmp3, sizeof(double));
			//file4.read((char *) &tmp4, sizeof(double));

			//hzcl[i][j] = tmp1;
			ez[i][j] = tmp2;
			ex[i][j] = tmp3;
			ey[i][j] =tmp4;
			//ez[i][j] = tmp1;
			//data[k] = tmp2;
			 k++;
			//hz[i][j] = tmp2 - tmp1;
   }
	cout << gridX*gridY*sizeof(double) << " BYTE red succesfully" << endl;
	//file2.close();
	//file3.close();
	//file4.close();
	     
};

int writeNCF(int from, int to)
{
	 /* IDs for the netCDF file, dimensions, and variables. */
   int ncid, lon_dimid, y_dimid, x_dimid, t_dimid;
   int x_varid, y_varid, hz_varid, ez_varid, ex_varid, ey_varid, t_varid;
   int dimids[3];

   /* The start and count arrays will tell the netCDF library where to
      write our data. */
   size_t start[3], count[3];


   /* Loop indexes. */
   int lvl, lat, lon, rec, i = 0;
   
   /* Error handling. */
   int retval;

   string fileName;
     stringstream s;
	 s << "K:\\ncf\\resTest.nc";
	 fileName = s.str();

   /* Create the file. */
	 if ((retval = nc_create(fileName.c_str(), NC_CLOBBER, &ncid)))
      ERR(retval);

   /* Define the dimensions. The record dimension is defined to have
    * unlimited length - it can grow as needed. In this example it is
    * the time dimension.*/
   if ((retval = nc_def_dim(ncid, "xdim", gridX, &x_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, "ydim", gridY, &y_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, "time", NC_UNLIMITED, &t_dimid)))
      ERR(retval);

   /* Define the coordinate variables. We will only define coordinate
      variables for lat and lon.  Ordinarily we would need to provide
      an array of dimension IDs for each variable's dimensions, but
      since coordinate variables only have one dimension, we can
      simply provide the address of that dimension ID (&lat_dimid) and
      similarly for (&lon_dimid). */
   if ((retval = nc_def_var(ncid, "xdim", NC_DOUBLE, 1, &x_dimid, 
			    &x_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "ydim", NC_DOUBLE, 1, &y_dimid, 
			    &y_varid)))
      ERR(retval);

   if ((retval = nc_def_var(ncid, "time", NC_DOUBLE, 1, &t_dimid, 
			    &t_varid)))
      ERR(retval);


   /* The dimids array is used to pass the dimids of the dimensions of
      the netCDF variables. Both of the netCDF variables we are
      creating share the same four dimensions. In C, the
      unlimited dimension must come first on the list of dimids. */
   dimids[0] = t_dimid;
   dimids[1] = x_dimid;
   dimids[2] = y_dimid;

   /* Define the netCDF variables for the pressure and temperature
    * data. */
   if ((retval = nc_def_var(ncid, "Hz", NC_DOUBLE, 3, 
			    dimids, &hz_varid)))
      ERR(retval);
   //if ((retval = nc_def_var(ncid, "Hx", NC_DOUBLE, 3, 
			//    dimids, &ex_varid)))
   //   ERR(retval);
   //if ((retval = nc_def_var(ncid, "Hy", NC_DOUBLE, 3, 
			//    dimids, &ey_varid)))
   //   ERR(retval);
   //if ((retval = nc_def_var(ncid, "Ez", NC_DOUBLE, 3, 
			//    dimids, &ez_varid)))
   //   ERR(retval);



   /* End define mode. */
   if ((retval = nc_enddef(ncid)))
      ERR(retval);

   /* Write the coordinate variable data. This will put the latitudes
      and longitudes of our data grid into the netCDF file. */
   if ((retval = nc_put_var_double(ncid, x_varid, &xDim[0])))
      ERR(retval);
   if ((retval = nc_put_var_double(ncid, y_varid, &yDim[0])))
      ERR(retval);

   /* These settings tell netcdf to write one timestep of data. (The
     setting of start[0] inside the loop below tells netCDF which
     timestep to write.) */
   count[0] = 1;
   count[1] = gridX;
   count[2] = gridY;
   start[0] = 0;
   start[1] = 0;
   start[2] = 0;

   /* Write the pretend data. This will write our surface pressure and
      surface temperature data. The arrays only hold one timestep worth
      of data. We will just rewrite the same data for each timestep. In
      a real application, the data would change between timesteps. */
   int nmax = to - from;
   for (rec = 0; rec < nmax; rec++)
   {
      start[0] = rec;
	  readFile(rec+from, "K:\\Ez\\Ezres");
      if ((retval = nc_put_vara_double(ncid, hz_varid, start, count, 
				      &ez[0][0])))
	 ERR(retval);
	 // if ((retval = nc_put_vara_double(ncid, ex_varid, start, count, 
		//		      &ex[0][0])))
	 //ERR(retval);
	 // if ((retval = nc_put_vara_double(ncid, ey_varid, start, count, 
		//		      &ey[0][0])))
	 //ERR(retval);

	 // if ((retval = nc_put_vara_double(ncid, ez_varid, start, count, 
		//		      &ez[0][0])))
	 //ERR(retval);

   }

   /* Close the file. */
   if ((retval = nc_close(ncid)))
      ERR(retval);
   
   printf("*** SUCCESS writing  file %s!\n", "Hz.nc");
   return 0;
}

int CountQ(int nmax)
{
	 /* IDs for the netCDF file, dimensions, and variables. */
   int ncid, lon_dimid, y_dimid, x_dimid, t_dimid;
   int x_varid, y_varid, hz_varid, ez_varid, ex_varid, ey_varid, t_varid;
   int dimids[1];
   int l1_varid,l2_varid,l3_varid,l4_varid,l5_varid,l6_varid,l7_varid;
   int l11_varid,l22_varid,l33_varid,l44_varid,l55_varid,l66_varid,l77_varid;

   /* The start and count arrays will tell the netCDF library where to
      write our data. */
   size_t start[1], count[1];


   /* Loop indexes. */
   int lvl, lat, lon, rec, i = 0;
   
   /* Error handling. */
   int retval;

   string fileName;
     stringstream s;
	 s << "K:\\ncf\\TEQAll.nc";
	 fileName = s.str();

   /* Create the file. */
	 if ((retval = nc_create(fileName.c_str(), NC_CLOBBER, &ncid)))
      ERR(retval);

   /* Define the dimensions. The record dimension is defined to have
    * unlimited length - it can grow as needed. In this example it is
    * the time dimension.*/

   if ((retval = nc_def_dim(ncid, "time", NC_UNLIMITED, &t_dimid)))
      ERR(retval);

   /* The dimids array is used to pass the dimids of the dimensions of
      the netCDF variables. Both of the netCDF variables we are
      creating share the same four dimensions. In C, the
      unlimited dimension must come first on the list of dimids. */
   dimids[0] = t_dimid;

   /* Define the netCDF variables for the pressure and temperature
    * data. */
   if ((retval = nc_def_var(ncid, "Lambda 1", NC_DOUBLE, 1, 
			    dimids, &l1_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "Lambda 2", NC_DOUBLE, 1, 
			    dimids, &l2_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "Lambda 3", NC_DOUBLE, 1, 
			    dimids, &l3_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "Lambda 4", NC_DOUBLE, 1, 
			    dimids, &l4_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "Lambda 5", NC_DOUBLE, 1, 
			    dimids, &l5_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "Lambda 6", NC_DOUBLE, 1, 
			    dimids, &l6_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "Lambda 7", NC_DOUBLE, 1, 
			    dimids, &l7_varid)))
      ERR(retval);

   if ((retval = nc_def_var(ncid, "Lambda 11", NC_DOUBLE, 1, 
			    dimids, &l11_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "Lambda 22", NC_DOUBLE, 1, 
			    dimids, &l22_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "Lambda 33", NC_DOUBLE, 1, 
			    dimids, &l33_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "Lambda 44", NC_DOUBLE, 1, 
			    dimids, &l44_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "Lambda 55", NC_DOUBLE, 1, 
			    dimids, &l55_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "Lambda 66", NC_DOUBLE, 1, 
			    dimids, &l66_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "Lambda 77", NC_DOUBLE, 1, 
			    dimids, &l77_varid)))
      ERR(retval);

   if ((retval = nc_def_var(ncid, "Time", NC_DOUBLE, 1, 
			    dimids, &t_varid)))
      ERR(retval);



   /* End define mode. */
   if ((retval = nc_enddef(ncid)))
      ERR(retval);

   count[0] = 1;
   start[1] = 0;


   /* Write the pretend data. This will write our surface pressure and
      surface temperature data. The arrays only hold one timestep worth
      of data. We will just rewrite the same data for each timestep. In
      a real application, the data would change between timesteps. */
   double Qpr=0;
   double Qprpr=0;
   for(rec = 15; rec < nmax-1; rec+=24)
   {
		readFile(rec,"K:\\Hz\\120");
		for(int i = 0; i < gridX; i++)
			for(int j = 0; j < gridY; j++)
			{
				expr[i][j] = ex[i][j];
				eypr[i][j] = ey[i][j];
				hzpr[i][j] = hz[i][j];
			}

		readFile(rec+1, "K:\\Hz\\120");

		double energy = 0;
		double p =0;
		double dx = LAMBDA/120;
		double dt = dx / (2* CC);
		for(int j = 15; j < 375; j++)
			for(int i = 0; i < gridX; i++)
			{
				//double en = (ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j])/2.0 + hz[i][j]*hz[i][j]/2.0;
				double e1 = (expr[i][j]*expr[i][j] + eypr[i][j]*eypr[i][j])/2.0 + hzpr[i][j]*hzpr[i][j]/2.0;
				double e2 = (ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j])/2.0 + hz[i][j]*hz[i][j]/2.0;
				//double en = hz[i][j]*hz[i][j]/2.0;
				//double e1 = hzpr[i][j];
				//double e2 = hz[i][j];
				energy+= e2;
				p+= fabs(e2-e1)/(5.0*dt);
			}
		//cout << "Energy: " << energy << endl;
		//cout << "P: " << p << endl;
			p = abs(p)/((double)(gridX*360));
		double Q;
		double rQ;
		if(p!=0)
		Q = 2*PI*FREQ*energy/p;
		else
			Q = 0;
		cout << "Q=" <<Q<<endl;
		Q=(Q+Qpr+Qprpr)/3.0;
		rQ = 1.0/Q;
		Qprpr = Qpr;
		Qpr = Q;



      start[0] = (rec-15)/24;
	  double tm = start[0]/2;
      if ((retval = nc_put_vara_double(ncid, l1_varid, start, count, 
				      &Q)))
	 ERR(retval);

	  if ((retval = nc_put_vara_double(ncid, l11_varid, start, count, 
				      &rQ)))
	 ERR(retval);

	   if ((retval = nc_put_vara_double(ncid, t_varid, start, count, 
				      &tm)))
	 ERR(retval);

   }

   Qpr=0;
   Qprpr=0;
   for(rec = 15; rec < nmax-1; rec+=24)
   {
		readFile(rec,"K:\\Hz\\180");
		for(int i = 0; i < gridX; i++)
			for(int j = 0; j < gridY; j++)
			{
				expr[i][j] = ex[i][j];
				eypr[i][j] = ey[i][j];
				hzpr[i][j] = hz[i][j];
			}

		readFile(rec+1, "K:\\Hz\\180");

		double energy = 0;
		double p =0;
		double dx = LAMBDA/120;
		double dt = dx / (2* CC);
		for(int j = 15; j < 375; j++)
			for(int i = 0; i < gridX; i++)
			{
				//double en = (ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j])/2.0 + hz[i][j]*hz[i][j]/2.0;
				double e1 = (expr[i][j]*expr[i][j] + eypr[i][j]*eypr[i][j])/2.0 + hzpr[i][j]*hzpr[i][j]/2.0;
				double e2 = (ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j])/2.0 + hz[i][j]*hz[i][j]/2.0;
				//double en = hz[i][j]*hz[i][j]/2.0;
				//double e1 = hzpr[i][j];
				//double e2 = hz[i][j];
				energy+= e2;
				p+= fabs(e2-e1)/(5.0*dt);
			}
		//cout << "Energy: " << energy << endl;
		//cout << "P: " << p << endl;
			p = abs(p)/((double)(gridX*360));
		double Q;
		double rQ;
		if(p!=0)
		Q = 2*PI*FREQ*energy/p;
		else
			Q = 0;
		cout << "Q=" <<Q<<endl;
		Q=(Q+Qpr+Qprpr)/3.0;
		rQ = 1.0/Q;
		Qprpr = Qpr;
		Qpr = Q;



      start[0] = (rec-15)/24;
	  double tm = start[0]/2;
      if ((retval = nc_put_vara_double(ncid, l2_varid, start, count, 
				      &Q)))
	 ERR(retval);

	  if ((retval = nc_put_vara_double(ncid, l22_varid, start, count, 
				      &rQ)))
	 ERR(retval);

	   if ((retval = nc_put_vara_double(ncid, t_varid, start, count, 
				      &tm)))
	 ERR(retval);

   }

   Qpr=0;
   Qprpr=0;
   for(rec = 15; rec < nmax-1; rec+=24)
   {
		readFile(rec,"K:\\Hz\\240");
		for(int i = 0; i < gridX; i++)
			for(int j = 0; j < gridY; j++)
			{
				expr[i][j] = ex[i][j];
				eypr[i][j] = ey[i][j];
				hzpr[i][j] = hz[i][j];
			}

		readFile(rec+1, "K:\\Hz\\240");

		double energy = 0;
		double p =0;
		double dx = LAMBDA/120;
		double dt = dx / (2* CC);
		for(int j = 15; j < 375; j++)
			for(int i = 0; i < gridX; i++)
			{
				//double en = (ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j])/2.0 + hz[i][j]*hz[i][j]/2.0;
				double e1 = (expr[i][j]*expr[i][j] + eypr[i][j]*eypr[i][j])/2.0 + hzpr[i][j]*hzpr[i][j]/2.0;
				double e2 = (ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j])/2.0 + hz[i][j]*hz[i][j]/2.0;
				//double en = hz[i][j]*hz[i][j]/2.0;
				//double e1 = hzpr[i][j];
				//double e2 = hz[i][j];
				energy+= e2;
				p+= fabs(e2-e1)/(5.0*dt);
			}
		//cout << "Energy: " << energy << endl;
		//cout << "P: " << p << endl;
			p = abs(p)/((double)(gridX*360));
		double Q;
		double rQ;
		if(p!=0)
		Q = 2*PI*FREQ*energy/p;
		else
			Q = 0;
		cout << "Q=" <<Q<<endl;
		Q=(Q+Qpr+Qprpr)/3.0;
		rQ = 1.0/Q;
		Qprpr = Qpr;
		Qpr = Q;



      start[0] = (rec-15)/24;
	  double tm = start[0]/2;
      if ((retval = nc_put_vara_double(ncid, l3_varid, start, count, 
				      &Q)))
	 ERR(retval);

	  if ((retval = nc_put_vara_double(ncid, l33_varid, start, count, 
				      &rQ)))
	 ERR(retval);

	   if ((retval = nc_put_vara_double(ncid, t_varid, start, count, 
				      &tm)))
	 ERR(retval);

   }

   Qpr=0;
   Qprpr=0;
   for(rec = 15; rec < nmax-1; rec+=24)
   {
		readFile(rec,"K:\\Hz\\300");
		for(int i = 0; i < gridX; i++)
			for(int j = 0; j < gridY; j++)
			{
				expr[i][j] = ex[i][j];
				eypr[i][j] = ey[i][j];
				hzpr[i][j] = hz[i][j];
			}

		readFile(rec+1, "K:\\Hz\\300");

		double energy = 0;
		double p =0;
		double dx = LAMBDA/120;
		double dt = dx / (2* CC);
		for(int j = 15; j < 375; j++)
			for(int i = 0; i < gridX; i++)
			{
				//double en = (ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j])/2.0 + hz[i][j]*hz[i][j]/2.0;
				double e1 = (expr[i][j]*expr[i][j] + eypr[i][j]*eypr[i][j])/2.0 + hzpr[i][j]*hzpr[i][j]/2.0;
				double e2 = (ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j])/2.0 + hz[i][j]*hz[i][j]/2.0;
				//double en = hz[i][j]*hz[i][j]/2.0;
				//double e1 = hzpr[i][j];
				//double e2 = hz[i][j];
				energy+= e2;
				p+= fabs(e2-e1)/(5.0*dt);
			}
		//cout << "Energy: " << energy << endl;
		//cout << "P: " << p << endl;
			p = abs(p)/((double)(gridX*360));
		double Q;
		double rQ;
		if(p!=0)
		Q = 2*PI*FREQ*energy/p;
		else
			Q = 0;
		cout << "Q=" <<Q<<endl;
		Q=(Q+Qpr+Qprpr)/3.0;
		rQ = 1.0/Q;
		Qprpr = Qpr;
		Qpr = Q;



      start[0] = (rec-15)/24;
	  double tm = start[0]/2;
      if ((retval = nc_put_vara_double(ncid, l4_varid, start, count, 
				      &Q)))
	 ERR(retval);

	  if ((retval = nc_put_vara_double(ncid, l44_varid, start, count, 
				      &rQ)))
	 ERR(retval);

	   if ((retval = nc_put_vara_double(ncid, t_varid, start, count, 
				      &tm)))
	 ERR(retval);

   }
   Qpr=0;
   Qprpr=0;
   for(rec = 15; rec < nmax-1; rec+=24)
   {
		readFile(rec,"K:\\Hz\\360");
		for(int i = 0; i < gridX; i++)
			for(int j = 0; j < gridY; j++)
			{
				expr[i][j] = ex[i][j];
				eypr[i][j] = ey[i][j];
				hzpr[i][j] = hz[i][j];
			}

		readFile(rec+1, "K:\\Hz\\360");

		double energy = 0;
		double p =0;
		double dx = LAMBDA/120;
		double dt = dx / (2* CC);
		for(int j = 15; j < 375; j++)
			for(int i = 0; i < gridX; i++)
			{
				//double en = (ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j])/2.0 + hz[i][j]*hz[i][j]/2.0;
				double e1 = (expr[i][j]*expr[i][j] + eypr[i][j]*eypr[i][j])/2.0 + hzpr[i][j]*hzpr[i][j]/2.0;
				double e2 = (ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j])/2.0 + hz[i][j]*hz[i][j]/2.0;
				//double en = hz[i][j]*hz[i][j]/2.0;
				//double e1 = hzpr[i][j];
				//double e2 = hz[i][j];
				energy+= e2;
				p+= fabs(e2-e1)/(5.0*dt);
			}
		//cout << "Energy: " << energy << endl;
		//cout << "P: " << p << endl;
			p = abs(p)/((double)(gridX*360));
		double Q;
		double rQ;
		if(p!=0)
		Q = 2*PI*FREQ*energy/p;
		else
			Q = 0;
		cout << "Q=" <<Q<<endl;
		Q=(Q+Qpr+Qprpr)/3.0;
		rQ = 1.0/Q;
		Qprpr = Qpr;
		Qpr = Q;



      start[0] = (rec-15)/24;
	  double tm = start[0]/2;
      if ((retval = nc_put_vara_double(ncid, l5_varid, start, count, 
				      &Q)))
	 ERR(retval);

	  if ((retval = nc_put_vara_double(ncid, l55_varid, start, count, 
				      &rQ)))
	 ERR(retval);

	   if ((retval = nc_put_vara_double(ncid, t_varid, start, count, 
				      &tm)))
	 ERR(retval);

   }

   Qpr=0;
   Qprpr=0;
   for(rec = 15; rec < nmax-1; rec+=24)
   {
		readFile(rec,"K:\\Hz\\420");
		for(int i = 0; i < gridX; i++)
			for(int j = 0; j < gridY; j++)
			{
				expr[i][j] = ex[i][j];
				eypr[i][j] = ey[i][j];
				hzpr[i][j] = hz[i][j];
			}

		readFile(rec+1, "K:\\Hz\\420");

		double energy = 0;
		double p =0;
		double dx = LAMBDA/120;
		double dt = dx / (2* CC);
		for(int j = 15; j < 375; j++)
			for(int i = 0; i < gridX; i++)
			{
				//double en = (ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j])/2.0 + hz[i][j]*hz[i][j]/2.0;
				double e1 = (expr[i][j]*expr[i][j] + eypr[i][j]*eypr[i][j])/2.0 + hzpr[i][j]*hzpr[i][j]/2.0;
				double e2 = (ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j])/2.0 + hz[i][j]*hz[i][j]/2.0;
				//double en = hz[i][j]*hz[i][j]/2.0;
				//double e1 = hzpr[i][j];
				//double e2 = hz[i][j];
				energy+= e2;
				p+= fabs(e2-e1)/(5.0*dt);
			}
		//cout << "Energy: " << energy << endl;
		//cout << "P: " << p << endl;
			p = abs(p)/((double)(gridX*360));
		double Q;
		double rQ;
		if(p!=0)
		Q = 2*PI*FREQ*energy/p;
		else
			Q = 0;
		cout << "Q=" <<Q<<endl;
		Q=(Q+Qpr+Qprpr)/3.0;
		rQ = 1.0/Q;
		Qprpr = Qpr;
		Qpr = Q;



      start[0] = (rec-15)/24;
	  double tm = start[0]/2;
      if ((retval = nc_put_vara_double(ncid, l6_varid, start, count, 
				      &Q)))
	 ERR(retval);

	  if ((retval = nc_put_vara_double(ncid, l66_varid, start, count, 
				      &rQ)))
	 ERR(retval);

	   if ((retval = nc_put_vara_double(ncid, t_varid, start, count, 
				      &tm)))
	 ERR(retval);

   }

   Qpr=0;
   Qprpr=0;
   for(rec = 15; rec < nmax-1; rec+=24)
   {
		readFile(rec,"K:\\Hz\\480");
		for(int i = 0; i < gridX; i++)
			for(int j = 0; j < gridY; j++)
			{
				expr[i][j] = ex[i][j];
				eypr[i][j] = ey[i][j];
				hzpr[i][j] = hz[i][j];
			}

		readFile(rec+1, "K:\\Hz\\480");

		double energy = 0;
		double p =0;
		double dx = LAMBDA/120;
		double dt = dx / (2* CC);
		for(int j = 15; j < 375; j++)
			for(int i = 0; i < gridX; i++)
			{
				//double en = (ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j])/2.0 + hz[i][j]*hz[i][j]/2.0;
				double e1 = (expr[i][j]*expr[i][j] + eypr[i][j]*eypr[i][j])/2.0 + hzpr[i][j]*hzpr[i][j]/2.0;
				double e2 = (ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j])/2.0 + hz[i][j]*hz[i][j]/2.0;
				//double en = hz[i][j]*hz[i][j]/2.0;
				//double e1 = hzpr[i][j];
				//double e2 = hz[i][j];
				energy+= e2;
				p+= fabs(e2-e1)/(5.0*dt);
			}
		//cout << "Energy: " << energy << endl;
		//cout << "P: " << p << endl;
			p = abs(p)/((double)(gridX*360));
		double Q;
		double rQ;
		if(p!=0)
		Q = 2*PI*FREQ*energy/p;
		else
			Q = 0;
		cout << "Q=" <<Q<<endl;
		Q=(Q+Qpr+Qprpr)/3.0;
		rQ = 1.0/Q;
		Qprpr = Qpr;
		Qpr = Q;



      start[0] = (rec-15)/24;
	  double tm = start[0]/2;
      if ((retval = nc_put_vara_double(ncid, l7_varid, start, count, 
				      &Q)))
	 ERR(retval);

	  if ((retval = nc_put_vara_double(ncid, l77_varid, start, count, 
				      &rQ)))
	 ERR(retval);

	   if ((retval = nc_put_vara_double(ncid, t_varid, start, count, 
				      &tm)))
	 ERR(retval);

   }
      if ((retval = nc_close(ncid)))
      ERR(retval);
   
   printf("*** SUCCESS writing  file %s!\n", "Hz.nc");
}

int CalcKotiqFunc(int n)
{
	 int ncid, lambda_dimid, y_dimid, x_dimid, t_dimid;
   int lambda_varid, theta_varid, n_varid, F_varid;
   int dimids[1];

   int* Theta_varids;

   /* The start and count arrays will tell the netCDF library where to
      write our data. */
   size_t start[1], count[1];


   /* Loop indexes. */
   int lvl, lat, lon, rec, i = 0;
   
   /* Error handling. */
   int retval;

   Theta_varids = new int[n];

   string fileName;
     stringstream s;
	 s << "K:\\ncf\\KotoFunc3.nc";
	 fileName = s.str();

   /* Create the file. */
	 if ((retval = nc_create(fileName.c_str(), NC_CLOBBER, &ncid)))
      ERR(retval);

   /* Define the dimensions. The record dimension is defined to have
    * unlimited length - it can grow as needed. In this example it is
    * the time dimension.*/
   if ((retval = nc_def_dim(ncid, "lambdadim", gridX, &lambda_dimid)))
      ERR(retval);

   /* Define the coordinate variables. We will only define coordinate
      variables for lat and lon.  Ordinarily we would need to provide
      an array of dimension IDs for each variable's dimensions, but
      since coordinate variables only have one dimension, we can
      simply provide the address of that dimension ID (&lat_dimid) and
      similarly for (&lon_dimid). */
   if ((retval = nc_def_var(ncid, "lambda", NC_DOUBLE, 1, &lambda_dimid, 
			    &lambda_varid)))
      ERR(retval);




   /* The dimids array is used to pass the dimids of the dimensions of
      the netCDF variables. Both of the netCDF variables we are
      creating share the same four dimensions. In C, the
      unlimited dimension must come first on the list of dimids. */
   dimids[0] = lambda_dimid;

   /* Define the netCDF variables for the pressure and temperature
    * data. */
   for(int i = 1; i <= n; i++)
   {
	    string varname;
     stringstream varnamestream;
	 varnamestream << "Theta_" << i;
	 varname = varnamestream.str();
	 if ((retval = nc_def_var(ncid, varname.c_str(), NC_DOUBLE, 1, 
					dimids, &Theta_varids[i-1])))
		  ERR(retval);
   }

   //if ((retval = nc_def_var(ncid, "Ez", NC_DOUBLE, 3, 
			//    dimids, &ez_varid)))
   //   ERR(retval);



   /* End define mode. */
   if ((retval = nc_enddef(ncid)))
      ERR(retval);

    double lambdamin = -10;
   double lambdamax = 10;
   double D= 0.01;
   double T = 2.328;
   double step = (lambdamax-lambdamin)/gridX;
   for(int j = 0; j < gridX; j++)
	   {
		   xDim[j] = lambdamin + step*j;
   }
   /* Write the coordinate variable data. This will put the latitudes
      and longitudes of our data grid into the netCDF file. */
   if ((retval = nc_put_var_double(ncid, lambda_varid, &xDim[0])))
      ERR(retval);

   /* These settings tell netcdf to write one timestep of data. (The
     setting of start[0] inside the loop below tells netCDF which
     timestep to write.) */
   count[0] = gridX;
   start[0] = 0;

   /* Write the pretend data. This will write our surface pressure and
      surface temperature data. The arrays only hold one timestep worth
      of data. We will just rewrite the same data for each timestep. In
      a real application, the data would change between timesteps. */

   double* ThetaF;
   ThetaF = new double[gridX];
  
   
   for (i = 1; i <= n; i++)
   {
	   for(int j = 0; j < gridX; j++)
	   {
		   xDim[j] = lambdamin + step*j;
		   ThetaF[j] = (-T/i)*sqrt(xDim[j]*xDim[j]-(D*i*i+1)*(D*i*i+1)) + (-acos((D*i*i+1)/xDim[j]))/i + 2.0*M_PI/i;
	   }
      if ((retval = nc_put_var_double(ncid, Theta_varids[i-1],
				      &ThetaF[0])))
		ERR(retval);

   }

   /* Close the file. */
   if ((retval = nc_close(ncid)))
      ERR(retval);
   
   cout << "Done!";

}

int main(int argc, char** argv)
{
	for(int i = 0; i < gridX; i++)
		xDim[i] = i/	40.0;

	for(int i = 0; i < gridY; i++)
		yDim[i] = i/40.0;

	writeNCF(0, 2499);

	//CalcKotiqFunc(10); 

	//CountQ(2000);
	//double maxQ = 0;
	//double avrQ = 0;
	//for(int t = 960; t < 999; t++)
	//{
	//readFile(t);
	//for(int i = 0; i < gridX; i++)
	//	for(int j = 0; j < gridY; j++)
	//	{
	//		expr[i][j] = ex[i][j];
	//		eypr[i][j] = ey[i][j];
	//		hzpr[i][j] = hz[i][j];
	//	}

	//readFile(t+1);

	//double energy = 0;
	//double p =0;
	//double dx = LAMBDA/120;
	//double dt = dx / (2* CC);
	//for(int j = 15; j < 375; j++)
	//	for(int i = 0; i < gridX; i++)
	//	{
	//		double en = (ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j])/2.0 + hz[i][j]*hz[i][j]/2.0;
	//		double e1 = expr[i][j]*expr[i][j] + eypr[i][j]*eypr[i][j];
	//		double e2 = ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j];
	//		//double en = hz[i][j]*hz[i][j]/2.0;
	//		//double e1 = hzpr[i][j];
	//		//double e2 = hz[i][j];
	//		energy+= en;
	//		p+= (e2-e1)/(5.0*dt);
	//	}
	////cout << "Energy: " << energy << endl;
	////cout << "P: " << p << endl;
	//	p = abs(p);
	//double Q = 2.0*PI*FREQ*energy/p;
	//cout << "Q-factor: " << Q << endl;
	//avrQ+=Q;
	//if (maxQ < Q) maxQ = Q;
	//}
	//avrQ = avrQ/1000.0;
	//cout << "MAX Q-factor: " << maxQ << endl;
	//cout << "Aver Q-factor: " << avrQ << endl;
	return 1;
}