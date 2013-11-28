#pragma once
#include "stdhead.h"

class CField
{
public:
	double*** Ex;
	double*** Ey;
	double*** Ez;

	double*** Hx;
	double*** Hy;
	double*** Hz;

	int gridX;
	int gridY;
	int gridZ;

	CField(void);
	CField(int sizeX, int sizeY, int sizeZ);
	~CField(void);
};
