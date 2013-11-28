#include "FieldTE.h"

using namespace std;

CFieldTE::CFieldTE(void)
{
	cout << "Creating NULL 2D TE Field...";
	Ex = NULL;
	Ey = NULL;

	Hz = NULL;

	gridX = 0;
	gridY = 0;
	cout << "done." << endl;
}

CFieldTE::CFieldTE(int sizeX, int sizeY)
{
	cout << "Creating 2D TE Field...";

	Ex = new double*[sizeX];
	Ey = new double*[sizeX];

	if(!Ex || !Ey) 
	{ 
		cout << endl << "ERROR: Memory allocation! Phase1: E-Field";
		system("PAUSE");
		terminate(); 
	}

	Hz = new double*[sizeX];

	if(!Hz) 
	{ 
		cout << endl << "ERROR: Memory allocation! Phase1: H-Field";
		system("PAUSE");
		terminate(); 
	}

	for(int i = 0; i < sizeX; i++)
	{
		Ex[i] = new double[sizeY];
		Ey[i] = new double[sizeY];
		if(!Ex[i] || !Ey[i]) 
		{ 
			cout << endl << "ERROR: Memory allocation! Phase2: E-Field"; 
			system("PAUSE"); 
			terminate(); 
		}

		Hz[i] = new double[sizeY];
		if(!Hz[i]) 
		{ 
			cout << endl << "ERROR: Memory allocation! Phase2: H-Field"; 
			system("PAUSE"); 
			terminate(); 
		}
	}



	for(int i = 0; i < sizeX; i++)
		for(int j = 0; j < sizeY; j++)
		{
			Ex[i][j] = 0.0;
			Ey[i][j] = 0.0;

			Hz[i][j] = 0.0;
		}

	cout << "done." << endl;

	gridX = sizeX;
	gridY = sizeY;

	cout << gridX << " x " << gridY <<  " lattice created." << endl;
}

//int CFieldTE::InitNetCDF(int timestp, double dx)
//{
//	int ncid, x_dimid, y_dimid, varidEx,varidEy,varidHz,varidTEST,
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
//	retval = nc_def_var(ncid, "1D Test", NC_DOUBLE, 1, &x_dimid,&varidTEST);
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
//	double* TEST;
//	double* xx;
//	double* yy;
//	TEST = new double[gridX];
//	xx = new double[gridX];
//	yy = new double[gridY];
//	for(int i = 0; i < gridX; i++)
//	{
//		TEST[i] = Hz[i][gridY/2 +1];
//		xx[i] = i * dx;
//		yy[i] = i * dx;
//
//	}
//
//	retval = nc_put_var_double(ncid, varidTEST, &TEST[0]);
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

void CFieldTE::WriteToNetCDF(int timestp)
{
	string fileName;
	stringstream s;
	s << "C:\\hz\\hzcl\\Hz[" << timestp << "].dat";
	fileName = s.str();
	ofstream file(fileName.c_str(), ios::out | ios::binary);
	double tmp;

	if(!file.is_open())
		cout << "ERROR::Cannot open file for write!";

	cout << "Writing to file... ";
	for(int i = 0; i < gridX; i++)
		for(int j = 0; j < gridY; j++)
		{
			tmp = Hz[i][j];
			file.write((char *) &tmp, sizeof(double));
		}
	cout << gridX*gridY*sizeof(double) << " BYTE written succesfully" << endl;
	file.close();

}

CFieldTE::~CFieldTE(void)
{
}
