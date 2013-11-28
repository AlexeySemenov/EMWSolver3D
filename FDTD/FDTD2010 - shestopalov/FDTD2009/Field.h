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
	void WriteToBinAmp(int timestp, int l, int d, int r, int u, int t, int b);
	void WriteToBinEz(int timestp, int l, int d, int r, int u, int t, int b);
	void WriteToBinEx(int timestp, int l, int d, int r, int u, int t, int b);
	void WriteToBinEy(int id, int timestp, int l, int d, int r, int u, int t, int b, double*** addValue);

	void WriteToBinHz(int timestp, int l, int d, int r, int u, int t, int b);
	void WriteToBinHx(int timestp, int l, int d, int r, int u, int t, int b);
	void WriteToBinHy(int timestp, int l, int d, int r, int u, int t, int b);
	~CField(void);
};
