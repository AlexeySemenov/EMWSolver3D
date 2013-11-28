#include "Field1D.h"

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
void CField1D::WriteToBinEz(int timestp, int l, int r)
{
	string fileName;
	stringstream s;
	s << "C:\\Ez\\Ez1D[" << timestp << "].dat";
	fileName = s.str();
	ofstream file(fileName.c_str(), ios::out | ios::binary);
	double tmp;

	if(!file.is_open())
		cout << "ERROR::Cannot open file for write!";

	cout << "Writing to file... ";
	for(int i = l; i < gridX-r; i++)
		{
			tmp = Ez[i];
			file.write((char *) &tmp, sizeof(double));
		}
	cout << gridX*sizeof(double) << " BYTE written succesfully" << endl;
	file.close();

}

void CField1D::WriteToBinHy(int timestp, int l, int r)
{
	string fileName;
	stringstream s;
	s << "C:\\Ez\\Hy1D[" << timestp << "].dat";
	fileName = s.str();
	ofstream file(fileName.c_str(), ios::out | ios::binary);
	double tmp;

	if(!file.is_open())
		cout << "ERROR::Cannot open file for write!";

	cout << "Writing to file... ";
	for(int i = l; i < gridX-r; i++)
		{
			tmp = Hy[i];
			file.write((char *) &tmp, sizeof(double));
		}
	cout << gridX*sizeof(double) << " BYTE written succesfully" << endl;
	file.close();

}

CField1D::~CField1D(void)
{
}
