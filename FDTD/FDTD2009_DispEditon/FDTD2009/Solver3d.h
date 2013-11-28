#pragma once
#include "stdhead.h"
#include "Field.h"

class CSolver3d
{
public:
	CField* data;

	double*** Ca;
	double*** Cb;
	double*** Da;
	double*** Db;

	double*** Eps;
	double*** Sig;
	double*** SigS;
	double*** Mju;

	double* Hz_inc;
	double* Ey_inc;
	double* Ex_inc;
//===================================================


	double* H_inc;
	double* E_inc;

	double ey_low_m1;
    double ey_low_m2;
    double ey_high_m1;
    double ey_high_m2;

	double ex_low_m1;
    double ex_low_m2;
    double ex_high_m1;
    double ex_high_m2;

	int lTF;
	int rTF;
	int uTF;
	int bTF;

	double S; //Courant number;
	double S2D; //square root S
	double S2DS1D;
	int Nx;
	int Ny;
	int Nz;
	double dx1;
	double dlx;
	double dly;
	double dxy;
	int source_angle;
	int m0;
	int Oja;
	int Oia;
	double CosSource_Angle;
	double SinSource_Angle;
	double theta;
	double A,B,C,D;
	double ConstIncE;
    double ConstIncH;
	double Offset;
	int Ioffset;

//===================================================
	double dt;
	double dx;
	double dy;
	double dz;
//===================================================
	double*** CaFz;
	double*** CbFz;
	double*** CaGz;
	double*** CbGz;
	double*** CaEz;
	double*** CbEz;
	double*** CcEz;

	double*** CaFx;
	double*** CbFx;
	double*** CaGx;
	double*** CbGx;
	double*** CaEx;
	double*** CbEx;
	double*** CcEx;

	double*** CaFy;
	double*** CbFy;
	double*** CaGy;
	double*** CbGy;
	double*** CaEy;
	double*** CbEy;
	double*** CcEy;

	double*** Fz;
	double Fz_pr;
	double*** Gz;
	double Gz_pr;

	double*** Fx;
	double Fx_pr;
	double*** Gx;
	double Gx_pr;

	double*** Fy;
	double Fy_pr;
	double*** Gy;
	double Gy_pr;

	double Bx_pr;
	double By_pr;
	double Bz_pr;

	double*** DaBx;
	double*** DbBx;
	double*** DaHx;
	double*** DbHx;
	double*** DcHx;

	double*** DaBy;
	double*** DbBy;
	double*** DaHy;
	double*** DbHy;
	double*** DcHy;

	double*** DaBz;
	double*** DbBz;
	double*** DaHz;
	double*** DbHz;
	double*** DcHz;

	double*** Bx;
	double*** By;
	double*** Bz;

	double pmlLR;
	double pmlUD;
	double ka_max;
	int  exponent;
	double R_err;
	double boundaryFactor;
	double sigma_max_1;
	double sigma_max_2;
	double t0;
	double tw;
	double rmax;// = 0.00001;
	int orderbc;// = 2;
	int emiTime;

//===================================================


	double Msource;

	void Init3DFDTD_UPML (void);
	void FDTD3DSolver_UPML(int timestep);

	void TimeStepLoop(int N);
	void Debugout(void);

	double source(double t); 
	double gammafunc(double gamma);
	double PlaneWaveLookUpd(int i, int j);
	double interp1(double*y, double xi);

	double Sig_x(double x,bool offset);
	double Sig_y(double y,bool offset);
	double Sig_z(double z,bool offset);

	double K_x(double x);
	double K_y(double y);
	double K_z(double z);

	CSolver3d(int, int, int);
	CSolver3d(void);
	~CSolver3d(void);
};
