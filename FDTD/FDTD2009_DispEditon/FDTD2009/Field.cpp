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

CField::~CField(void)
{
}
