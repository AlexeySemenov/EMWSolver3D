#pragma once
#include "stdhead.h"

class CFieldTM
{
public:
	double** Hx;
	double** Hy;

	double** Ez;

	int gridX;
	int gridY;

	CFieldTM(void);
	CFieldTM(int sizeX, int sizeY);
	//int InitNetCDF(int timestp, double dx);
	void WriteToNetCDF(int timestp, int l, int r, int u, int d);
	void ReadFromFile(int timestp, int l, int r, int u, int d);
	void WriteToBinEz(int timestp, int l, int d, int r, int u);
	void WriteToBinHx(int timestp, int l, int d, int r, int u);
	void WriteToBinHy(int timestp, int l, int d, int r, int u);
	~CFieldTM(void);
};
