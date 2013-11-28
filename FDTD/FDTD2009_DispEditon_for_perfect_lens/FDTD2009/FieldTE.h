#pragma once
#include "stdhead.h"

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
	void WriteToNetCDF(int timestp);
	~CFieldTE(void);
};
