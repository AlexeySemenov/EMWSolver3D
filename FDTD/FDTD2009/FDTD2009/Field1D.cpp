#include "Field1D.h"
#include <netcdf.h>

using namespace std;

CField1D::CField1D(void)
{
	cout << "Creating NULL 1D Field...";
	Ez = NULL;
	Hy = NULL;

	gridX = 0;
	cout << "done." << endl;
}

CField1D::CField1D(int sizeX)
{
	cout << "Creating 1D Field...";

	Ez = new double[sizeX];
	Hy = new double[sizeX];

	if(!Ez || !Hy) 
	{ 
		cout << endl << "ERROR: Memory allocation! Phase1: E,H-Field";
		system("PAUSE");
		terminate(); 
	}

	for(int i = 0; i < sizeX; i++)
		{
			Ez[i] = 0.0;
			Hy[i] = 0.0;
		}
	gridX = sizeX;

}
int CField1D::InitNetCDF(int timestp, double dx)
{
	int ncid, x_dimid,  varidEz,varidHy, varidX;
    int dimids[1];
	int x, retval;
	string fileName;
	stringstream s;
	s << "C:\\maxwell\\Maxwell1D_x[" << timestp << "].nc";
	fileName = s.str();

	cout << "Writting to File...";

	retval = nc_create(fileName.c_str(), NC_CLOBBER, &ncid);

	
	retval = nc_def_dim(ncid, "x", gridX, &x_dimid);

	dimids[0] = x_dimid;

	retval = nc_def_var(ncid, "Ez-Field", NC_DOUBLE, 1, 
			    dimids, &varidEz);
	retval = nc_def_var(ncid, "Hy-Field", NC_DOUBLE, 1, 
			    dimids, &varidHy);

	retval = nc_def_var(ncid, "X Dim", NC_DOUBLE, 1, &x_dimid,&varidX);

	retval = nc_enddef(ncid);

	retval = nc_put_var_double(ncid, varidEz, &Ez[0]);
	retval = nc_put_var_double(ncid, varidHy, &Hy[0]);

	double* xx;
	xx = new double[gridX];
	for(int i = 0; i < gridX; i++)
	{
		xx[i] = i * dx;
	}

	retval = nc_put_var_double(ncid, varidX, &xx[0]);

	retval = nc_close(ncid);
	cout << "done." << endl;

	return 1;
}

CField1D::~CField1D(void)
{
}
