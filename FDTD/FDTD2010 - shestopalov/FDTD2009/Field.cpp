#include "Field.h"

using namespace std;

CField::CField(void)
{
	cout << "Creating NULL 3D Field...";
	Ex = NULL;
	Ey = NULL;
	Ez = NULL;

	Hx = NULL;
	Hy = NULL;
	Hz = NULL;

	gridX = 0;
	gridY = 0;
	gridZ = 0;
	cout << "done." << endl;
}

CField::CField(int sizeX, int sizeY, int sizeZ)
{
	cout << "Creating 3D Field...";

	Ex = new double**[sizeX];
	Ey = new double**[sizeX];
	Ez = new double**[sizeX];

	if(!Ex || !Ey || !Ez) 
	{ 
		cout << endl << "ERROR: Memory allocation! Phase1: E-Field";
		system("PAUSE");
		terminate(); 
	}

	Hx = new double**[sizeX];
	Hy = new double**[sizeX];
	Hz = new double**[sizeX];

	if(!Hx || !Hy || !Hz) 
	{ 
		cout << endl << "ERROR: Memory allocation! Phase1: H-Field";
		system("PAUSE");
		terminate(); 
	}

	for(int i = 0; i < sizeX; i++)
	{
		Ex[i] = new double*[sizeY];
		Ey[i] = new double*[sizeY];
		Ez[i] = new double*[sizeY];
		if(!Ex[i] || !Ey[i] || !Ez[i]) 
		{ 
			cout << endl << "ERROR: Memory allocation! Phase2: E-Field"; 
			system("PAUSE"); 
			terminate(); 
		}

		Hx[i] = new double*[sizeY];
		Hy[i] = new double*[sizeY];
		Hz[i] = new double*[sizeY];
		if(!Hx[i] || !Hy[i] || !Hz[i]) 
		{ 
			cout << endl << "ERROR: Memory allocation! Phase2: H-Field"; 
			system("PAUSE"); 
			terminate(); 
		}
	}

	for(int i = 0; i < sizeX; i++)
		for(int j = 0; j < sizeY; j++)
		{
			Ex[i][j] = new double[sizeZ];
			Ey[i][j] = new double[sizeZ];
			Ez[i][j] = new double[sizeZ];
			if(!Ex[i][j] || !Ey[i][j] || !Ez[i][j]) 
			{ 
				cout << endl << "ERROR: Memory allocation! Phase3: E-Field"; 
				system("PAUSE");
				terminate(); 
			}

			Hx[i][j] = new double[sizeZ];
			Hy[i][j] = new double[sizeZ];
			Hz[i][j] = new double[sizeZ];
			if(!Hx[i][j] || !Hy[i][j] || !Hz[i][j]) 
			{ 
				cout << endl << "ERROR: Memory allocation! Phase3: H-Field"; 
				system("PAUSE");
				terminate(); 
			}
		}


	for(int i = 0; i < sizeX; i++)
		for(int j = 0; j < sizeY; j++)
			for(int k = 0; k < sizeZ; k++)
			{
				Ex[i][j][k] = 0.0;
				Ey[i][j][k] = 0.0;
				Ez[i][j][k] = 0.0;

				Hx[i][j][k] = 0.0;
				Hy[i][j][k] = 0.0;
				Hz[i][j][k] = 0.0;
			}

	cout << "done." << endl;

	gridX = sizeX;
	gridY = sizeY;
	gridZ = sizeZ;

	cout << gridX << " x " << gridY << " x " << gridZ << " lattice created." << endl;


}

void CField::WriteToBinAmp(int timestp, int l, int d, int r, int u, int t, int b)
{
	string fileName;
	stringstream s;
	s << "K:\\Ez\\Amp[" << timestp << "].dat";
	fileName = s.str();
	ofstream file(fileName.c_str(), ios::out | ios::binary);
	double tmp;

	if(!file.is_open())
		cout << "ERROR::Cannot open file for write!";

	cout << "Writing to file... ";
	for(int i = l; i < gridX-r; i++)
		for(int j = d; j < gridY-u; j++)
			for(int k = b; k < gridZ-t; k++)
			{
				if(i==gridX/2)
				{
					tmp = Ey[i][j][k]*Ey[i][j][k];
					file.write((char *) &tmp, sizeof(double));
				}
			}
	cout << gridX*gridY*gridZ*sizeof(double) << " BYTE written succesfully" << endl;
	file.close();

}

void CField::WriteToBinEz(int timestp, int l, int d, int r, int u, int t, int b)
{
	string fileName;
	stringstream s;
	s << "C:\\Ez\\Ez[" << timestp << "].dat";
	fileName = s.str();
	ofstream file(fileName.c_str(), ios::out | ios::binary);
	double tmp;

	if(!file.is_open())
		cout << "ERROR::Cannot open file for write!";

	cout << "Writing to file... ";
	for(int i = l; i < gridX-r; i++)
		for(int j = d; j < gridY-u; j++)
			for(int k = b; k < gridZ-t; k++)
			{
				tmp = Ez[i][j][k];
				file.write((char *) &tmp, sizeof(double));
			}
	cout << gridX*gridY*gridZ*sizeof(double) << " BYTE written succesfully" << endl;
	file.close();

}
void CField::WriteToBinEy(int id, int timestp, int l, int d, int r, int u, int t, int b, double*** addValue)
{
	string fileName;
	stringstream s;
	s << "C:\\Ez\\Ey[" << timestp << "]" << id << ".dat";
	fileName = s.str();
	ofstream file(fileName.c_str(), ios::out | ios::binary);
	double tmp;

	if(!file.is_open())
		cout << "ERROR::Cannot open file for write!";

	cout << "Writing to file... ";
	for(int i = l; i < gridX-r; i++)
		for(int j = d; j < gridY-u; j++)
			for(int k = b; k < gridZ-t; k++)
			{
				tmp = Ey[i][j][k] + addValue[i][j][k];
				file.write((char *) &tmp, sizeof(double));
			}
	cout << gridX*gridY*gridZ*sizeof(double) << " BYTE written succesfully" << endl;
	file.close();

}
void CField::WriteToBinEx(int timestp, int l, int d, int r, int u, int t, int b)
{
	string fileName;
	stringstream s;
	s << "C:\\Ez\\Ex[" << timestp << "].dat";
	fileName = s.str();
	ofstream file(fileName.c_str(), ios::out | ios::binary);
	double tmp;

	if(!file.is_open())
		cout << "ERROR::Cannot open file for write!";

	cout << "Writing to file... ";
	for(int i = l; i < gridX-r; i++)
		for(int j = d; j < gridY-u; j++)
			for(int k = b; k < gridZ-t; k++)
			{
				tmp = Ex[i][j][k];
				file.write((char *) &tmp, sizeof(double));
			}
	cout << gridX*gridY*gridZ*sizeof(double) << " BYTE written succesfully" << endl;
	file.close();

}

void CField::WriteToBinHz(int timestp, int l, int d, int r, int u, int t, int b)
{
	string fileName;
	stringstream s;
	s << "C:\\Ez\\Hz[" << timestp << "].dat";
	fileName = s.str();
	ofstream file(fileName.c_str(), ios::out | ios::binary);
	double tmp;

	if(!file.is_open())
		cout << "ERROR::Cannot open file for write!";

	cout << "Writing to file... ";
	for(int i = l; i < gridX-r; i++)
		for(int j = d; j < gridY-u; j++)
			for(int k = b; k < gridZ-t; k++)
			{
				tmp = Hz[i][j][k];
				file.write((char *) &tmp, sizeof(double));
			}
	cout << gridX*gridY*gridZ*sizeof(double) << " BYTE written succesfully" << endl;
	file.close();

}
void CField::WriteToBinHy(int timestp, int l, int d, int r, int u, int t, int b)
{
	string fileName;
	stringstream s;
	s << "C:\\Ez\\Hy[" << timestp << "].dat";
	fileName = s.str();
	ofstream file(fileName.c_str(), ios::out | ios::binary);
	double tmp;

	if(!file.is_open())
		cout << "ERROR::Cannot open file for write!";

	cout << "Writing to file... ";
	for(int i = l; i < gridX-r; i++)
		for(int j = d; j < gridY-u; j++)
			for(int k = b; k < gridZ-t; k++)
			{
				tmp = Hy[i][j][k];
				file.write((char *) &tmp, sizeof(double));
			}
	cout << gridX*gridY*gridZ*sizeof(double) << " BYTE written succesfully" << endl;
	file.close();

}
void CField::WriteToBinHx(int timestp, int l, int d, int r, int u, int t, int b)
{
	string fileName;
	stringstream s;
	s << "C:\\Ez\\Hx[" << timestp << "].dat";
	fileName = s.str();
	ofstream file(fileName.c_str(), ios::out | ios::binary);
	double tmp;

	if(!file.is_open())
		cout << "ERROR::Cannot open file for write!";

	cout << "Writing to file... ";
	for(int i = l; i < gridX-r; i++)
		for(int j = d; j < gridY-u; j++)
			for(int k = b; k < gridZ-t; k++)
			{
				tmp = Hx[i][j][k];
				file.write((char *) &tmp, sizeof(double));
			}
	cout << gridX*gridY*gridZ*sizeof(double) << " BYTE written succesfully" << endl;
	file.close();

}
CField::~CField(void)
{
	cout << "Creating 3D Field...";
	for(int i = 0; i < gridX; i++)
		for(int j = 0; j < gridY; j++)
		{
			delete[] Ex[i][j];
			delete[] Ey[i][j];
			delete[] Ez[i][j];
			delete[] Hx[i][j];
			delete[] Hy[i][j];
			delete[] Hz[i][j];
		}

	for(int i = 0; i < gridX; i++)
	{
		delete[] Ex[i];
		delete[] Ey[i];
		delete[] Ez[i];
		delete[] Hx[i];
		delete[] Hy[i];
		delete[] Hz[i];
	}

	delete[] Ex;
	delete[] Ey;
	delete[] Ez;

	delete[] Hx;
	delete[] Hy;
	delete[] Hz;

}
