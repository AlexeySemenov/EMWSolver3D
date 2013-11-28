#pragma once
#include "stdhead.h"
using namespace std;
class CFieldTE
{
public:
	double** Ex;
	double** Ey;

	double** Hz;

	int gridX;
	int gridY;

	CFieldTE(void);
	CFieldTE(int sizeX, int sizeY);
	int InitNetCDF(int timestp, double dx);
	void WriteToBinHz(int timestp, int l, int d, int r, int u, string path);
	void WriteToBinEx(int timestp, int l, int d, int r, int u, string path);
	void WriteToBinEy(int timestp, int l, int d, int r, int u, string path);

	~CFieldTE(void);
};
