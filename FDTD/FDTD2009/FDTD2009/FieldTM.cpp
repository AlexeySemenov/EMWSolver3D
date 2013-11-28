#include "FieldTM.h"
#include <netcdf.h>

using namespace std;

CFieldTM::CFieldTM(void)
{
	cout << "Creating NULL 2D TM Field...";
	Hx = NULL;
	Hy = NULL;

	Ez = NULL;

	gridX = 0;
	gridY = 0;
	cout << "done." << endl;
}

CFieldTM::CFieldTM(int sizeX, int sizeY)
{
	cout << "Creating 2D TM Field...";

	Hx = new double*[sizeX];
	Hy = new double*[sizeX];

	if(!Hx || !Hy) 
	{ 
		cout << endl << "ERROR: Memory allocation! Phase1: H-Field";
		system("PAUSE");
		terminate(); 
	}

	Ez = new double*[sizeX];

	if(!Ez) 
	{ 
		cout << endl << "ERROR: Memory allocation! Phase1: E-Field";
		system("PAUSE");
		terminate(); 
	}

	for(int i = 0; i < sizeX; i++)
	{
		Hx[i] = new double[sizeY];
		Hy[i] = new double[sizeY];
		if(!Hx[i] || !Hy[i]) 
		{ 
			cout << endl << "ERROR: Memory allocation! Phase2: H-Field"; 
			system("PAUSE"); 
			terminate(); 
		}

		Ez[i] = new double[sizeY];
		if(!Ez[i]) 
		{ 
			cout << endl << "ERROR: Memory allocation! Phase2: E-Field"; 
			system("PAUSE"); 
			terminate(); 
		}
	}



	for(int i = 0; i < sizeX; i++)
		for(int j = 0; j < sizeY; j++)
		{
			Hx[i][j] = 0.0;
			Hy[i][j] = 0.0;

			Ez[i][j] = 0.0;
		}

	cout << "done." << endl;

	gridX = sizeX;
	gridY = sizeY;

	cout << gridX << " x " << gridY <<  " lattice created." << endl;
}

//int CFieldTM::InitNetCDF(int timestp, double dx)
//{
//	int ncid, x_dimid, y_dimid, varidEx,varidEy,varidHz,varidTMST,
//		varidX, varidY;
//    int dimids[2];
//	int x, y, retval;
//	string fileName;
//	stringstream s;
//	s << "C:\\maxwell\\Maxwell_xy[" << timestp << "].nc";
//	fileName = s.str();
//
//	cout << "Writting to File...";
//
//	retval = nc_create(fileName.c_str(), NC_CLOBBER, &ncid);
//
//	//NcFile dataFile("Maxwell_xy.nc", NcFile::Replace);
//	
//	retval = nc_def_dim(ncid, "x", gridX, &x_dimid);
//	//NcDim* xDim = dataFile.add_dim("x", gridX);
//	retval = nc_def_dim(ncid, "y", gridY, &y_dimid);
//    //NcDim* yDim = dataFile.add_dim("y", gridY);
//
//	dimids[0] = x_dimid;
//    dimids[1] = y_dimid;
//
//	retval = nc_def_var(ncid, "Ex-Field", NC_DOUBLE, 2, 
//			    dimids, &varidEx);
//	retval = nc_def_var(ncid, "Ey-Field", NC_DOUBLE, 2, 
//			    dimids, &varidEy);
//	retval = nc_def_var(ncid, "Hz-Field", NC_DOUBLE, 2, 
//			    dimids, &varidHz);
//	retval = nc_def_var(ncid, "1D Test", NC_DOUBLE, 1, &x_dimid,&varidTMST);
//
//	retval = nc_def_var(ncid, "X Dim", NC_DOUBLE, 1, &x_dimid,&varidX);
//	retval = nc_def_var(ncid, "Y Dim", NC_DOUBLE, 1, &y_dimid,&varidY);
//
//	retval = nc_enddef(ncid);
//
//	retval = nc_put_var_double(ncid, varidEx, &Ex[0][0]);
//	retval = nc_put_var_double(ncid, varidEy, &Ey[0][0]);
//	retval = nc_put_var_double(ncid, varidHz, &Hz[0][0]);
//
//	double* TMST;
//	double* xx;
//	double* yy;
//	TMST = new double[gridX];
//	xx = new double[gridX];
//	yy = new double[gridY];
//	for(int i = 0; i < gridX; i++)
//	{
//		TMST[i] = Hz[i][gridY/2 +1];
//		xx[i] = i * dx;
//		yy[i] = i * dx;
//
//	}
//
//	retval = nc_put_var_double(ncid, varidTMST, &TMST[0]);
//	retval = nc_put_var_double(ncid, varidX, &xx[0]);
//	retval = nc_put_var_double(ncid, varidY, &yy[0]);
//
//	retval = nc_close(ncid);
//	cout << "done." << endl;
//
//	//NcVar *dataEx = dataFile.add_var("Ex-Field", ncDouble, xDim, yDim);
//	//NcVar *dataEy = dataFile.add_var("Ey-Field", ncDouble, xDim, yDim);
//	//NcVar *dataHz = dataFile.add_var("Hz-Field", ncDouble, xDim, yDim);
//	//dataEx->put(&Ex[0][0], gridX, gridY);
//	//dataEy->put(&Ey[0][0], gridX, gridY);
//	//dataHz->put(&Hz[0][0], gridX, gridY);
//
//	return 1;
//}

void CFieldTM::WriteToNetCDF(int timestp, int l, int r, int u, int d)
{
	string fileName;
	stringstream s;
	s << "C:\\hz\\hzpml\\Hz[" << timestp << "].dat";
	fileName = s.str();
	ofstream file(fileName.c_str(), ios::out | ios::binary);
	double tmp;

	if(!file.is_open())
		cout << "ERROR::Cannot open file for write!";

	cout << "Writing to file... ";
	for(int i = l; i < gridX-r; i++)
		for(int j = d; j < gridY-u; j++)
		{
			tmp = Ez[i][j];
			file.write((char *) &tmp, sizeof(double));
		}
	cout << gridX*gridY*sizeof(double) << " BYTM written succesfully" << endl;
	file.close();

}

void CFieldTM::ReadFromFile(int timestp, int l, int r, int u, int d)
{
	string fileName;
	stringstream s;
	s << "C:\\hz\\hzpml\\Hz[" << timestp << "].dat";
	fileName = s.str();
	ifstream file(fileName.c_str(), ios::out);
	double tmp;

	if(!file.is_open())
		cout << "ERROR::Cannot open file for write!";

	cout << "Reading From file... ";
	for(int i = l; i < gridX-r; i++)
		for(int j = d; j < gridY-u; j++)
		{
			file.read((char *) &tmp, sizeof(double));
			Ez[i][j] = tmp;
			
		}
	cout << gridX*gridY*sizeof(double) << " BYTE read succesfully" << endl;
	file.close();

}
CFieldTM::~CFieldTM(void)
{
}
