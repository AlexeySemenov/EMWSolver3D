#pragma once
#include "stdhead.h"

class CField1D
{
public:
	double* Ez;
	double* Hy;
	int gridX;

	CField1D(void);
	CField1D(int sizeX);
	int InitNetCDF(int timestp, double dx);
	~CField1D(void);
};
