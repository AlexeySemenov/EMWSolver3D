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
	void WriteToBinEz(int timestp, int l, int r);
	void WriteToBinHy(int timestp, int l, int r);
	~CField1D(void);
};
