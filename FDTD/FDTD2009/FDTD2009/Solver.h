#pragma once
#include "stdhead.h"
#include "Field.h"
#include "FieldTE.h"
#include "FieldTM.h"
#include "Field1D.h"

class CSolver
{
public:
	CField* data;
	CFieldTE* dataTE;
	CFieldTM* dataTM;
	CField1D* data1D;

	double** Ca;
	double** Cb;
	double** Da;
	double** Db;

	double** Eps;
	double** Sig;
	double** SigS;
	double** Mju;

	double* Hz_inc;
	double* Ey_inc;
	double* Ex_inc;
//===================================================
	double* Ca1D;
	double* Cb1D;
	double* Da1D;
	double* Db1D;

	double* Eps1D;
	double* Sig1D;
	double* SigS1D;
	double* Mju1D;

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
	double dx1;
	double dlx;
	double dly;
	double dxy;
	double MurConstY;
	double MurConstX;
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
	double** CaFz;
	double** CbFz;
	double** CaGz;
	double** CbGz;
	double** CaEz;
	double** CbEz;
	double** CcEz;

	double** Fz;
	double Fz_pr;
	double** Gz;
	double Gz_pr;
	double Bx_pr;
	double By_pr;

	double** DaBx;
	double** DbBx;
	double** DaHx;
	double** DbHx;
	double** DcHx;

	double** DaBy;
	double** DbBy;
	double** DaHy;
	double** DbHy;
	double** DcHy;

	double** Bx;
	double** By;

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

	void InitSimpleFDTD (void);
	void SimpleFDTDSolver(int timestep);

	void InitSimpleFDTD_TF (void);
	void SimpleFDTDSolver_TF(int timestep);

	void InitSimpleFDTD_TM (void);
	void SimpleFDTDSolver_TM(int timestep);

	void InitSimpleFDTD_TM_UPML (void);
	void SimpleFDTDSolver_TM_UPML(int timestep);

	void Init3DFDTD_TM_UPML (void);
	void FDTD3DSolver_TM_UPML(int timestep);

	void Init1DSimpleFDTD (void);
	void Simple1DFDTDSolver(int timestep);

	void PointMIncidence(int x, int y, int time);
	void HardSourceMIncidence(int x, int y, int time);

	void FDTD1DIncidence(int x, int time);

	void TimeStepLoop(int N);
	void ContinueTimeStepLoop(int start, int N);
	void Debugout(void);

	double source(double t); 
	double gammafunc(double gamma);
	double PlaneWaveLookUpd(int i, int j);
	double interp1(double*y, double xi);

	double Sig_x(double x,bool offset);
	double Sig_y(double y,bool offset);
	double Sig_z(double z);

	double K_x(double x);
	double K_y(double y);
	double K_z(double z);
	CSolver(void);
	CSolver(const char* par, int, int, int);
	~CSolver(void);
};
