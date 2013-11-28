#include "stdhead.h"
#include "Solver.h"

using namespace std;

CSolver::CSolver(void)
{
	data = new CField(100, 100, 100);

}

CSolver::CSolver(const char* par, int SizeX, int SizeY, int SizeZ)
{
	gridX = SizeX;
	gridY = SizeY;
	if (par == "3D") 
	{
	}

	if (par == "TE")
	{
		dataTE = new CFieldTE(SizeX, SizeY);
		Ca = new double*[dataTE->gridX];
		Cb = new double*[dataTE->gridX];
		Da = new double*[dataTE->gridX];
		Db = new double*[dataTE->gridX];

		Eps = new double*[dataTE->gridX];
		Mju = new double*[dataTE->gridX];
		Sig = new double*[dataTE->gridX];
		SigS = new double*[dataTE->gridX];

		Fz = new double*[dataTE->gridX];
		Gz = new double*[dataTE->gridX];
		Fzeven = new double*[dataTE->gridX];
		Gzeven = new double*[dataTE->gridX];
		Fzodd = new double*[dataTE->gridX];
		Gzodd = new double*[dataTE->gridX];

		Bx = new double*[dataTE->gridX];
		By = new double*[dataTE->gridX];
		Tx = new double*[dataTE->gridX];
		Ty = new double*[dataTE->gridX];

		Bxeven = new double*[dataTE->gridX];
		Byeven = new double*[dataTE->gridX];
		Txeven = new double*[dataTE->gridX];
		Tyeven = new double*[dataTE->gridX];

		Bxodd = new double*[dataTE->gridX];
		Byodd = new double*[dataTE->gridX];
		Txodd = new double*[dataTE->gridX];
		Tyodd = new double*[dataTE->gridX];

		Cb = new double*[dataTE->gridX];
		Db = new double*[dataTE->gridX];

		CaFz = new double*[dataTE->gridX];
		CbFz = new double*[dataTE->gridX];
		CaGz = new double*[dataTE->gridX];
		CbGz = new double*[dataTE->gridX];
		CaEz = new double*[dataTE->gridX];
		CbEz = new double*[dataTE->gridX];
		CcEz = new double*[dataTE->gridX];
		DaBx = new double*[dataTE->gridX];
		DbBx = new double*[dataTE->gridX];
		DaBy = new double*[dataTE->gridX];
		DbBy = new double*[dataTE->gridX];
		DaHx = new double*[dataTE->gridX];
		DbHx = new double*[dataTE->gridX];
		DcHx = new double*[dataTE->gridX];
		DaHy = new double*[dataTE->gridX];
		DbHy = new double*[dataTE->gridX];
		DcHy = new double*[dataTE->gridX];

		CaFx = new double*[dataTE->gridX];
		CbFx = new double*[dataTE->gridX];
		CcFx = new double*[dataTE->gridX];
		CdFx = new double*[dataTE->gridX];
		CeFx = new double*[dataTE->gridX];
		CfFx = new double*[dataTE->gridX];

		CaFy = new double*[dataTE->gridX];
		CbFy = new double*[dataTE->gridX];
		CcFy = new double*[dataTE->gridX];
		CdFy = new double*[dataTE->gridX];
		CeFy = new double*[dataTE->gridX];
		CfFy = new double*[dataTE->gridX];

		CaT = new double*[dataTE->gridX];
		CbT = new double*[dataTE->gridX];
		CcT = new double*[dataTE->gridX];
		CdT = new double*[dataTE->gridX];
		CeT = new double*[dataTE->gridX];
		CfT = new double*[dataTE->gridX];

		gammaE = new double*[dataTE->gridX];
		gammaM = new double*[dataTE->gridX];
		omegaPE = new double*[dataTE->gridX];
		omegaPM = new double*[dataTE->gridX];

		for(int i = 0; i < dataTE->gridX; i++)
		{
			Fz[i] = new double[dataTE->gridY];
			Gz[i] = new double[dataTE->gridY];
			Fzeven[i] = new double[dataTE->gridY];
			Gzeven[i] = new double[dataTE->gridY];
			Fzodd[i] = new double[dataTE->gridY];
			Gzodd[i] = new double[dataTE->gridY];

			Bx[i] = new double[dataTE->gridY];
			By[i] = new double[dataTE->gridY];
			Tx[i] = new double[dataTE->gridY];
			Ty[i] = new double[dataTE->gridY];

			Bxeven[i] = new double[dataTE->gridY];
			Byeven[i] = new double[dataTE->gridY];
			Txeven[i] = new double[dataTE->gridY];
			Tyeven[i] = new double[dataTE->gridY];

			Bxodd[i] = new double[dataTE->gridY];
			Byodd[i] = new double[dataTE->gridY];
			Txodd[i] = new double[dataTE->gridY];
			Tyodd[i] = new double[dataTE->gridY];

			Cb[i] = new double[dataTE->gridY];
			Db[i] = new double[dataTE->gridY];


			CaFz[i] = new double[dataTE->gridY];
			CbFz[i] = new double[dataTE->gridY];
			CaGz[i] = new double[dataTE->gridY];
			CbGz[i] = new double[dataTE->gridY];
			CaEz[i] = new double[dataTE->gridY];
			CbEz[i] = new double[dataTE->gridY];
			CcEz[i] = new double[dataTE->gridY];
			DaBx[i] = new double[dataTE->gridY];
			DbBx[i] = new double[dataTE->gridY];
			DaBy[i] = new double[dataTE->gridY];
			DbBy[i] = new double[dataTE->gridY];
			DaHx[i] = new double[dataTE->gridY];
			DbHx[i] = new double[dataTE->gridY];
			DcHx[i] = new double[dataTE->gridY];
			DaHy[i] = new double[dataTE->gridY];
			DbHy[i] = new double[dataTE->gridY];
			DcHy[i] = new double[dataTE->gridY];

			CaFx[i] = new double[dataTE->gridY];
			CbFx[i] = new double[dataTE->gridY];
			CcFx[i] = new double[dataTE->gridY];
			CdFx[i] = new double[dataTE->gridY];
			CeFx[i] = new double[dataTE->gridY];
			CfFx[i] = new double[dataTE->gridY];

			CaFy[i] = new double[dataTE->gridY];
			CbFy[i] = new double[dataTE->gridY];
			CcFy[i] = new double[dataTE->gridY];
			CdFy[i] = new double[dataTE->gridY];
			CeFy[i] = new double[dataTE->gridY];
			CfFy[i] = new double[dataTE->gridY];

			CaT[i] = new double[dataTE->gridY];
			CbT[i] = new double[dataTE->gridY];
			CcT[i] = new double[dataTE->gridY];
			CdT[i] = new double[dataTE->gridY];
			CeT[i] = new double[dataTE->gridY];
			CfT[i] = new double[dataTE->gridY];

			gammaE[i] = new double[dataTE->gridY];
			gammaM[i] = new double[dataTE->gridY];
			omegaPE[i] = new double[dataTE->gridY];
			omegaPM[i] = new double[dataTE->gridY];

			Eps[i] = new double[dataTE->gridY];
			Mju[i] = new double[dataTE->gridY];
			Sig[i] = new double[dataTE->gridY];
			SigS[i] = new double[dataTE->gridY];
		}


		for(int i = 0; i < dataTE->gridX; i++)
		{
			Ca[i] = new double[dataTE->gridY];
			Cb[i] = new double[dataTE->gridY];
			Da[i] = new double[dataTE->gridY];
			Db[i] = new double[dataTE->gridY];

		}

	}
	if (par == "TM")
	{
		dataTM = new CFieldTM(1000, 1000);
		Ca = new double*[dataTM->gridX];
		Cb = new double*[dataTM->gridX];
		Da = new double*[dataTM->gridX];
		Db = new double*[dataTM->gridX];

		Eps = new double*[dataTM->gridX];
		Mju = new double*[dataTM->gridX];
		Sig = new double*[dataTM->gridX];
		SigS = new double*[dataTM->gridX];

		for(int i = 0; i < dataTM->gridX; i++)
		{
			Ca[i] = new double[dataTM->gridY];
			Cb[i] = new double[dataTM->gridY];
			Da[i] = new double[dataTM->gridY];
			Db[i] = new double[dataTM->gridY];

			Eps[i] = new double[dataTM->gridY];
			Mju[i] = new double[dataTM->gridY];
			Sig[i] = new double[dataTM->gridY];
			SigS[i] = new double[dataTM->gridY];
		}

	}
	if (par == "TM_UPML")
	{
		dataTM = new CFieldTM(SizeX, SizeY);
		Fz = new double*[dataTM->gridX];
		Gz = new double*[dataTM->gridX];
		Fzeven = new double*[dataTM->gridX];
		Gzeven = new double*[dataTM->gridX];
		Fzodd = new double*[dataTM->gridX];
		Gzodd = new double*[dataTM->gridX];

		Bx = new double*[dataTM->gridX];
		By = new double*[dataTM->gridX];
		Tx = new double*[dataTM->gridX];
		Ty = new double*[dataTM->gridX];

		Bxeven = new double*[dataTM->gridX];
		Byeven = new double*[dataTM->gridX];
		Txeven = new double*[dataTM->gridX];
		Tyeven = new double*[dataTM->gridX];

		Bxodd = new double*[dataTM->gridX];
		Byodd = new double*[dataTM->gridX];
		Txodd = new double*[dataTM->gridX];
		Tyodd = new double*[dataTM->gridX];

		Cb = new double*[dataTM->gridX];
		Db = new double*[dataTM->gridX];

		CaFz = new double*[dataTM->gridX];
		CbFz = new double*[dataTM->gridX];
		CaGz = new double*[dataTM->gridX];
		CbGz = new double*[dataTM->gridX];
		CaEz = new double*[dataTM->gridX];
		CbEz = new double*[dataTM->gridX];
		CcEz = new double*[dataTM->gridX];
		DaBx = new double*[dataTM->gridX];
		DbBx = new double*[dataTM->gridX];
		DaBy = new double*[dataTM->gridX];
		DbBy = new double*[dataTM->gridX];
		DaHx = new double*[dataTM->gridX];
		DbHx = new double*[dataTM->gridX];
		DcHx = new double*[dataTM->gridX];
		DaHy = new double*[dataTM->gridX];
		DbHy = new double*[dataTM->gridX];
		DcHy = new double*[dataTM->gridX];

		CaFx = new double*[dataTM->gridX];
		CbFx = new double*[dataTM->gridX];
		CcFx = new double*[dataTM->gridX];
		CdFx = new double*[dataTM->gridX];
		CeFx = new double*[dataTM->gridX];
		CfFx = new double*[dataTM->gridX];

		CaFy = new double*[dataTM->gridX];
		CbFy = new double*[dataTM->gridX];
		CcFy = new double*[dataTM->gridX];
		CdFy = new double*[dataTM->gridX];
		CeFy = new double*[dataTM->gridX];
		CfFy = new double*[dataTM->gridX];

		CaT = new double*[dataTM->gridX];
		CbT = new double*[dataTM->gridX];
		CcT = new double*[dataTM->gridX];
		CdT = new double*[dataTM->gridX];
		CeT = new double*[dataTM->gridX];
		CfT = new double*[dataTM->gridX];

		Eps = new double*[dataTM->gridX];
		Mju = new double*[dataTM->gridX];
		Sig = new double*[dataTM->gridX];
		SigS = new double*[dataTM->gridX];

		gammaE = new double*[dataTM->gridX];
		gammaM = new double*[dataTM->gridX];
		omegaPE = new double*[dataTM->gridX];
		omegaPM = new double*[dataTM->gridX];

		for(int i = 0; i < dataTM->gridX; i++)
		{
			Fz[i] = new double[dataTM->gridY];
			Gz[i] = new double[dataTM->gridY];
			Fzeven[i] = new double[dataTM->gridY];
			Gzeven[i] = new double[dataTM->gridY];
			Fzodd[i] = new double[dataTM->gridY];
			Gzodd[i] = new double[dataTM->gridY];

			Bx[i] = new double[dataTM->gridY];
			By[i] = new double[dataTM->gridY];
			Tx[i] = new double[dataTM->gridY];
			Ty[i] = new double[dataTM->gridY];

			Bxeven[i] = new double[dataTM->gridY];
			Byeven[i] = new double[dataTM->gridY];
			Txeven[i] = new double[dataTM->gridY];
			Tyeven[i] = new double[dataTM->gridY];

			Bxodd[i] = new double[dataTM->gridY];
			Byodd[i] = new double[dataTM->gridY];
			Txodd[i] = new double[dataTM->gridY];
			Tyodd[i] = new double[dataTM->gridY];

			Cb[i] = new double[dataTM->gridY];
			Db[i] = new double[dataTM->gridY];


			CaFz[i] = new double[dataTM->gridY];
			CbFz[i] = new double[dataTM->gridY];
			CaGz[i] = new double[dataTM->gridY];
			CbGz[i] = new double[dataTM->gridY];
			CaEz[i] = new double[dataTM->gridY];
			CbEz[i] = new double[dataTM->gridY];
			CcEz[i] = new double[dataTM->gridY];
			DaBx[i] = new double[dataTM->gridY];
			DbBx[i] = new double[dataTM->gridY];
			DaBy[i] = new double[dataTM->gridY];
			DbBy[i] = new double[dataTM->gridY];
			DaHx[i] = new double[dataTM->gridY];
			DbHx[i] = new double[dataTM->gridY];
			DcHx[i] = new double[dataTM->gridY];
			DaHy[i] = new double[dataTM->gridY];
			DbHy[i] = new double[dataTM->gridY];
			DcHy[i] = new double[dataTM->gridY];

			CaFx[i] = new double[dataTM->gridY];
			CbFx[i] = new double[dataTM->gridY];
			CcFx[i] = new double[dataTM->gridY];
			CdFx[i] = new double[dataTM->gridY];
			CeFx[i] = new double[dataTM->gridY];
			CfFx[i] = new double[dataTM->gridY];

			CaFy[i] = new double[dataTM->gridY];
			CbFy[i] = new double[dataTM->gridY];
			CcFy[i] = new double[dataTM->gridY];
			CdFy[i] = new double[dataTM->gridY];
			CeFy[i] = new double[dataTM->gridY];
			CfFy[i] = new double[dataTM->gridY];

			CaT[i] = new double[dataTM->gridY];
			CbT[i] = new double[dataTM->gridY];
			CcT[i] = new double[dataTM->gridY];
			CdT[i] = new double[dataTM->gridY];
			CeT[i] = new double[dataTM->gridY];
			CfT[i] = new double[dataTM->gridY];

			gammaE[i] = new double[dataTM->gridY];
			gammaM[i] = new double[dataTM->gridY];
			omegaPE[i] = new double[dataTM->gridY];
			omegaPM[i] = new double[dataTM->gridY];

			Eps[i] = new double[dataTM->gridY];
			Mju[i] = new double[dataTM->gridY];
			Sig[i] = new double[dataTM->gridY];
			SigS[i] = new double[dataTM->gridY];
		}

	}
	if (par == "1D")
	{
		data1D = new CField1D(100);
		Ca1D = new double[data1D->gridX];
		Cb1D = new double[data1D->gridX];
		Da1D = new double[data1D->gridX];
		Db1D = new double[data1D->gridX];

		Eps1D = new double[data1D->gridX];
		Mju1D = new double[data1D->gridX];
		Sig1D = new double[data1D->gridX];
		SigS1D = new double[data1D->gridX];
	}
	

}

void CSolver::InitSimpleFDTD_TM_UPML(void)
{
	//Init timestep parameters
	S = 1;
	S2D = S / sqrt(2.0);
	Nx = dataTM->gridX;
	Ny = dataTM->gridY;
	//dx = 1.0e-10;
	//dy = 1.0e-10;
	dz = 0.0;
	dx = LAMBDA/120;
	dy = LAMBDA/120;
	dxy = sqrt(dx*dx+ dy*dy);
	dt = S* dx / (2* CC);

	if (dt > 1.0/(2.0*FREQ))
	{cout << "Warning::Computed dt is greater than Nyquist: " << dt << endl; system("PAUSE");}

	//First order MUr Update Constants
	MurConstY = (CC*dt-dy)/(CC*dt+dy);
	MurConstX = (CC*dt-dx)/(CC*dt+dx);

	//Init  source
	source_angle = 0;
	m0 = 29;
	lTF =30;
	rTF =390;
	bTF =35;
	uTF = 205;
	int ia;
	int ja;
	int ib;
	int jb;
	int mod;
	ia = lTF;
	ja = bTF;
	ib = rTF;
	jb = uTF;

	if ((source_angle >= 0) && (source_angle <= 90))
	{
		Oja = ja;
		Oia = ia;
	}
	else if ((source_angle >=90) && (source_angle <=180))
	{
		Oja = ja;
		Oia = ib;
	}
	else if ((source_angle >=180) && (source_angle <=270))
	{
		Oja = jb;
		Oia = ib;
	}
	else
	{
		Oja = jb;
		Oia = ia;
	}

	// Make Sure everythings adds up nicley
	mod=source_angle % 90;
	if (mod==0){
		if (source_angle == 0)
		{
			CosSource_Angle = 1;
			SinSource_Angle = 0;
		}
		else if (source_angle == 90)
		{
			CosSource_Angle = 0;
			SinSource_Angle = 1;
		}
		else if (source_angle == 180)
		{
			CosSource_Angle = -1;
			SinSource_Angle = 0;
		}
		else if (source_angle == 270)
		{
			CosSource_Angle = 0;
			SinSource_Angle = -1;
		}
	}
	else
	{
		CosSource_Angle = cos(source_angle*M_PI/180.0);
		SinSource_Angle = sin(source_angle*M_PI/180.0);
	}

	//theta = 90*M_PI/180;
 //   A =M_PI/Nx*CosSource_Angle*sin(theta);
 //   B =M_PI/Ny*SinSource_Angle*sin(theta);
 //   C =M_PI/Ny*SinSource_Angle*cos(theta);
 //   D =sin(M_PI*S2D/Nx)*sin(M_PI*S2D/Nx);
	S2DS1D = sqrt((CosSource_Angle*CosSource_Angle*CosSource_Angle*CosSource_Angle
		+SinSource_Angle*SinSource_Angle*SinSource_Angle*SinSource_Angle));

	dx1=S2DS1D*dx;
	ConstIncE = dt/dx1/EPS_Z;
	ConstIncH = dt/dx1/MU_Z;
	Offset = sqrt((ia*CosSource_Angle)*(ia*CosSource_Angle)+
		(ja*SinSource_Angle)*(ja*SinSource_Angle));


	emiTime = -10000000;

	//Init MEDIA
	
	double rcyl;
	double MinEps = 1.0;
	double omegasqr;
	double gammasqr;
	double zeps;
	acyl = 40;
	bcyl = 80;
	double LBslab = 520;
	double RBslab = 1019;
	for(int i = 0; i < dataTM->gridX; i++)
		for(int j = 0; j < dataTM->gridY; j++)
		{



			gammaE[i][j] = 0.0;
			gammaM[i][j] = 0.0;
			omegaPE[i][j] = 0.0;
			omegaPM[i][j] = 0.0;

			Eps[i][j] = 1;
			Mju[i][j] = 1;
			Sig[i][j] = 0;
			SigS[i][j] = 0;


			if((i>= LBslab) && (i<=RBslab) && (j>20) && (j<dataTM->gridY-20))
			{
				omegaPE[i][j] = sqrt(2.0)*OMEGA;
				omegaPM[i][j] = sqrt(2.0)*OMEGA;
				gammaE[i][j] = 0;//0.0005*OMEGA;
				gammaM[i][j] = 0;//0.0005*OMEGA;
			}



			Gz[i][j] = 0.0;
			Fz[i][j] = 0.0;
			Gzeven[i][j] = 0.0;
			Fzeven[i][j] = 0.0;
			Gzodd[i][j] = 0.0;
			Fzodd[i][j] = 0.0;

			Bx[i][j] = 0.0;
			By[i][j] = 0.0;
			Tx[i][j] = 0.0;
			Ty[i][j] = 0.0;

			Bxeven[i][j] = 0.0;
			Byeven[i][j] = 0.0;
			Txeven[i][j] = 0.0;
			Tyeven[i][j] = 0.0;

			Bxodd[i][j] = 0.0;
			Byodd[i][j] = 0.0;
			Txodd[i][j] = 0.0;
			Tyodd[i][j] = 0.0;
		}


	//Init UPML
		double mu_r_x_1, mu_r_x_2, eps_r_x_1,eps_r_x_2;
		
		double boundaryWidth, boundaryHeight;

		mu_r_x_1 = 1;
		mu_r_x_2 = 1;
		eps_r_x_1 = 1;
		eps_r_x_2 = 1;
		pmlUD = 20;
		pmlLR = 20;

		boundaryWidth = (double)pmlLR*dx;
		boundaryHeight = (double)pmlUD*dy;
		

		ka_max = 1.0;
		exponent = 6;
		R_err = 1e-16;
		double eta_1 = sqrt(MU_Z*mu_r_x_1/EPS_Z/eps_r_x_1);
		double eta_2 = sqrt(MU_Z*mu_r_x_2/EPS_Z/eps_r_x_2);
		sigma_max_1= -log(R_err) * (exponent + 1.0) / (2.0 * IMP_Z * boundaryWidth);
		sigma_max_2= -log(R_err) * (exponent + 1.0) / (2.0 * IMP_Z * boundaryHeight);;
		boundaryFactor = sigma_max_1 / ( dx * (pow(boundaryWidth,exponent)) * (exponent + 1));

		tw = 26.53e-12;
		t0 = 4*tw;

	//Init FDTDSolver Consts
		double z = 0;
		


	for(int i = 0; i < dataTM->gridX; i++)
		for(int j = 0; j < dataTM->gridY; j++)
		{
			tmpI = i;
			tmpJ = j;
			Cb[i][j] = (dt / (Eps[i][j]* EPS_Z)) / (1.0 + (Sig[i][j] * dt)/(2 * Eps[i][j]* EPS_Z));

			Db[i][j] = (dt / (Mju[i][j]* MU_Z)) / (1.0 + (SigS[i][j] * dt)/(2 * Mju[i][j]* MU_Z));


			CaFz[i][j] = (2*EPS_Z * Eps[i][j] - Sig[i][j]*dt)/
				(2*EPS_Z * Eps[i][j] + Sig[i][j]*dt);
			CbFz[i][j] = 2 * dt / (2*EPS_Z * Eps[i][j] + Sig[i][j]*dt);

			CaGz[i][j] = (2*EPS_Z * K_x(i) - Sig_x(i,true)*dt)/
				(2*EPS_Z * K_x(i) + Sig_x(i,true)*dt);
			CbGz[i][j] = 2*EPS_Z*dt/(2*EPS_Z * K_x(i) + Sig_x(i,true)*dt);

			CaEz[i][j] = (2*EPS_Z*K_y(j) - Sig_y(j,true)*dt)/
				(2*EPS_Z*K_y(j) + Sig_y(j,true)*dt);
			CbEz[i][j] = (2*EPS_Z * K_z(z) + Sig_z(z)*dt)/
				((2*EPS_Z* K_y(j) + Sig_y(j,true)*dt)*EPS_Z);
			CcEz[i][j] = (2*EPS_Z*K_z(z)-Sig_z(z)*dt)/
				((2*EPS_Z*K_y(j) + Sig_y(j,true)*dt)*EPS_Z);

			DaBx[i][j] = (2*EPS_Z*K_y(j) - Sig_y(j,false)*dt)/
				(2*EPS_Z*K_y(j) + Sig_y(j,false)*dt);
			DbBx[i][j] = 2*EPS_Z*dt/(2*EPS_Z*K_y(j) + Sig_y(j,false)*dt);
			DaHx[i][j] = (2*EPS_Z*K_z(z) - Sig_z(z)*dt)/
				(2*EPS_Z*K_z(z) + Sig_z(z)*dt);
			DbHx[i][j] = (2*EPS_Z*K_x(i) + Sig_x(i,true)*dt)/
				((2*EPS_Z*K_z(z) + Sig_z(z)*dt)*MU_Z);
			DcHx[i][j] = (2*EPS_Z*K_x(i) - Sig_x(i,true)*dt)/
				((2*EPS_Z*K_z(z) + Sig_z(z)*dt)*MU_Z);


			DaBy[i][j] = (2*EPS_Z*K_z(z) - Sig_z(z)*dt)/
				(2*EPS_Z*K_z(z) - Sig_z(z)*dt);
			DbBy[i][j] = 2*EPS_Z*dt/(2*EPS_Z*K_z(z) + Sig_z(z)*dt);
			DaHy[i][j] = (2*EPS_Z*K_z(i) - Sig_x(i,false)*dt)/
				(2*EPS_Z*K_x(i) + Sig_x(i,false)*dt);
			DbHy[i][j] = (2*EPS_Z*K_y(j) + Sig_y(j,true)*dt)/
				((2*EPS_Z*K_x(i) + Sig_x(i,false)*dt)*MU_Z);
			DcHy[i][j] = (2*EPS_Z*K_y(j) - Sig_y(j,true)*dt)/
				((2*EPS_Z*K_x(i) + Sig_x(i,false)*dt)*MU_Z);


			//==========================================================


			double plug = 1.0;
			///*if((i != 0) && (j != 0))
			//{
			//	CaFy[i][j] = (1.0/(plug)) + (dt*((gammaM[i][j-1]+gammaM[i-1][j-1])/2.0)/(2.0*plug));
			//	CbFy[i][j] = 2.0/( plug);
			//	CcFy[i][j] = ((2.0/1.0) - dt*dt*(pow(omegaPM[i][j-1],2)+pow(omegaPM[i-1][j-1], 2))/4.0);
			//	CdFy[i][j] = ((1.0/1.0) - (dt*((gammaM[i][j-1]+gammaM[i-1][j-1])/2.0)/(2.0)) + dt*dt*(pow(omegaPM[i][j-1],2)+pow(omegaPM[i-1][j-1], 2))/8.0);
			//	CeFy[i][j] = (1.0/( plug)) - (dt*((gammaM[i][j-1]+gammaM[i-1][j-1])/2.0)/(2.0*plug));
			//	CfFy[i][j] = ((1.0) + (dt*((gammaM[i][j-1]+gammaM[i-1][j-1])/2.0)/(2.0)) + dt*dt*(pow(omegaPM[i][j-1],2)+pow(omegaPM[i-1][j-1], 2))/8.0);
			//}
			//else
			//{
			//	CaFy[i][j] = (1.0/(plug)) + (dt*((gammaM[i][j])/1.0)/(2.0*plug));
			//	CbFy[i][j] = 2.0/( plug);
			//	CcFy[i][j] = ((2.0/1.0) - dt*dt*pow((omegaPM[i][j])/1.0, 2)/2.0);
			//	CdFy[i][j] = ((1.0/1.0) - (dt*((gammaM[i][j])/1.0)/(2.0)) + dt*dt*pow((omegaPM[i][j])/1.0, 2)/4.0);
			//	CeFy[i][j] = (1.0/( plug)) - (dt*((gammaM[i][j])/1.0)/(2.0*plug));
			//	CfFy[i][j] = ((1.0) + (dt*((gammaM[i][j])/1.0)/(2.0)) + dt*dt*pow((omegaPM[i][j])/1.0, 2)/4.0);			
			//}

			//if((j >= 1))
			//{
			//	CaFx[i][j] = (1.0/(plug)) + (dt*((gammaM[i][j-1]+gammaM[i][j])/2.0)/(2.0*plug));
			//	CbFx[i][j] = 2.0/( plug);
			//	CcFx[i][j] = ((2.0/1.0) - dt*dt*(pow(omegaPM[i][j-1],2)+pow(omegaPM[i][j], 2))/4.0);
			//	CdFx[i][j] = ((1.0/1.0) - (dt*((gammaM[i][j-1]+gammaM[i][j])/2.0)/(2.0)) + dt*dt*(pow(omegaPM[i][j-1],2)+pow(omegaPM[i][j], 2))/8.0);
			//	CeFx[i][j] = (1.0/( plug)) - (dt*((gammaM[i][j-1]+gammaM[i][j])/2.0)/(2.0*plug));
			//	CfFx[i][j] = ((1.0) + (dt*((gammaM[i][j-1]+gammaM[i][j])/2.0)/(2.0)) + dt*dt*(pow(omegaPM[i][j-1],2)+pow(omegaPM[i][j], 2))/8.0);
			//}
			//else
			//{
			//	CaFx[i][j] = (1.0/(plug)) + (dt*((gammaM[i][j])/1.0)/(2.0*plug));
			//	CbFx[i][j] = 2.0/( plug);
			//	CcFx[i][j] = ((2.0/1.0) - dt*dt*pow((omegaPM[i][j])/1.0, 2)/2.0);
			//	CdFx[i][j] = ((1.0/1.0) - (dt*((gammaM[i][j])/1.0)/(2.0)) + dt*dt*pow((omegaPM[i][j])/1.0, 2)/4.0);
			//	CeFx[i][j] = (1.0/( plug)) - (dt*((gammaM[i][j])/1.0)/(2.0*plug));
			//	CfFx[i][j] = ((1.0) + (dt*((gammaM[i][j])/1.0)/(2.0)) + dt*dt*pow((omegaPM[i][j])/1.0, 2)/4.0);			
			//}

			//if((j >= 1))
			//{
			//	CaT[i][j] = (1.0/(plug)) + (dt*gammaE[i][j-1]/(2.0*plug));
			//	CbT[i][j] = 2.0/( plug);
			//	CcT[i][j] = ((2.0/1.0) - dt*dt*omegaPE[i][j-1]*omegaPE[i][j-1]/2.0);
			//	CdT[i][j] = ((1.0/1.0) - (dt*gammaE[i][j-1]/(2.0)) + dt*dt*omegaPE[i][j-1]*omegaPE[i][j-1]/4.0);
			//	CeT[i][j] = (1.0/( plug)) - (dt*gammaE[i][j-1]/(2.0*plug));
			//	CfT[i][j] = ((1.0) + (dt*gammaE[i][j-1]/(2.0)) + dt*dt*omegaPE[i][j-1]*omegaPE[i][j-1]/4.0);
			//}
			//else
			//{
			//	CaT[i][j] = (1.0/(plug)) + (dt*gammaE[i][j]/(2.0*plug));
			//	CbT[i][j] = 2.0/( plug);
			//	CcT[i][j] = ((2.0/1.0) - dt*dt*omegaPE[i][j]*omegaPE[i][j]/2.0);
			//	CdT[i][j] = ((1.0/1.0) - (dt*gammaE[i][j]/(2.0)) + dt*dt*omegaPE[i][j]*omegaPE[i][j]/4.0);
			//	CeT[i][j] = (1.0/( plug)) - (dt*gammaE[i][j]/(2.0*plug));
			//	CfT[i][j] = ((1.0) + (dt*gammaE[i][j]/(2.0)) + dt*dt*omegaPE[i][j]*omegaPE[i][j]/4.0);
			//}		*/

			if((i != gridX-1))
			{
				CaFy[i][j] = (1.0/(plug)) + (dt*((gammaM[i][j]+gammaM[i+1][j])/2.0)/(2.0*plug));
				CbFy[i][j] = 2.0/( plug);
				CcFy[i][j] = ((2.0/1.0) - dt*dt*(pow(omegaPM[i][j],2)+pow(omegaPM[i+1][j], 2))/4.0);
				CdFy[i][j] = ((1.0/1.0) - (dt*((gammaM[i][j]+gammaM[i+1][j])/2.0)/(2.0)) + dt*dt*(pow(omegaPM[i][j],2)+pow(omegaPM[i+1][j], 2))/8.0);
				CeFy[i][j] = (1.0/( plug)) - (dt*((gammaM[i][j]+gammaM[i+1][j])/2.0)/(2.0*plug));
				CfFy[i][j] = ((1.0) + (dt*((gammaM[i][j]+gammaM[i+1][j])/2.0)/(2.0)) + dt*dt*(pow(omegaPM[i][j],2)+pow(omegaPM[i+1][j], 2))/8.0);
			}
			else
			{
				CaFy[i][j] = (1.0/(plug)) + (dt*((gammaM[i][j])/1.0)/(2.0*plug));
				CbFy[i][j] = 2.0/( plug);
				CcFy[i][j] = ((2.0/1.0) - dt*dt*pow((omegaPM[i][j])/1.0, 2)/2.0);
				CdFy[i][j] = ((1.0/1.0) - (dt*((gammaM[i][j])/1.0)/(2.0)) + dt*dt*pow((omegaPM[i][j])/1.0, 2)/4.0);
				CeFy[i][j] = (1.0/( plug)) - (dt*((gammaM[i][j])/1.0)/(2.0*plug));
				CfFy[i][j] = ((1.0) + (dt*((gammaM[i][j])/1.0)/(2.0)) + dt*dt*pow((omegaPM[i][j])/1.0, 2)/4.0);			
			}

			if((j != gridY-1))
			{
				CaFx[i][j] = (1.0/(plug)) + (dt*((gammaM[i][j+1]+gammaM[i][j])/2.0)/(2.0*plug));
				CbFx[i][j] = 2.0/( plug);
				CcFx[i][j] = ((2.0/1.0) - dt*dt*(pow(omegaPM[i][j+1],2)+pow(omegaPM[i][j], 2))/4.0);
				CdFx[i][j] = ((1.0/1.0) - (dt*((gammaM[i][j+1]+gammaM[i][j])/2.0)/(2.0)) + dt*dt*(pow(omegaPM[i][j+1],2)+pow(omegaPM[i][j], 2))/8.0);
				CeFx[i][j] = (1.0/( plug)) - (dt*((gammaM[i][j+1]+gammaM[i][j])/2.0)/(2.0*plug));
				CfFx[i][j] = ((1.0) + (dt*((gammaM[i][j+1]+gammaM[i][j])/2.0)/(2.0)) + dt*dt*(pow(omegaPM[i][j+1],2)+pow(omegaPM[i][j], 2))/8.0);
			}
			else
			{
				CaFx[i][j] = (1.0/(plug)) + (dt*((gammaM[i][j])/1.0)/(2.0*plug));
				CbFx[i][j] = 2.0/( plug);
				CcFx[i][j] = ((2.0/1.0) - dt*dt*pow((omegaPM[i][j])/1.0, 2)/2.0);
				CdFx[i][j] = ((1.0/1.0) - (dt*((gammaM[i][j])/1.0)/(2.0)) + dt*dt*pow((omegaPM[i][j])/1.0, 2)/4.0);
				CeFx[i][j] = (1.0/( plug)) - (dt*((gammaM[i][j])/1.0)/(2.0*plug));
				CfFx[i][j] = ((1.0) + (dt*((gammaM[i][j])/1.0)/(2.0)) + dt*dt*pow((omegaPM[i][j])/1.0, 2)/4.0);			
			}

			CaT[i][j] = (1.0/(plug)) + (dt*gammaE[i][j]/(2.0*plug));
			CbT[i][j] = 2.0/( plug);
			CcT[i][j] = ((2.0/1.0) - dt*dt*omegaPE[i][j]*omegaPE[i][j]/2.0);
			CdT[i][j] = ((1.0/1.0) - (dt*gammaE[i][j]/(2.0)) + dt*dt*omegaPE[i][j]*omegaPE[i][j]/4.0);
			CeT[i][j] = (1.0/( plug)) - (dt*gammaE[i][j]/(2.0*plug));
			CfT[i][j] = ((1.0) + (dt*gammaE[i][j]/(2.0)) + dt*dt*omegaPE[i][j]*omegaPE[i][j]/4.0);	

		}

		
		//	

		//}
	//TF/SF init
		//Ez_inc = new double[dataTM->gridX];
		//Hy_inc = new double[dataTM->gridX];
		//Hx_inc = new double[dataTM->gridX];
		H_inc = new double[2*(dataTM->gridX)];
		E_inc = new double[2*(dataTM->gridX)];

		//for(int i = 0; i < dataTM->gridX; i++)
		//{
		//	cout <<i <<": " << Eps[i][Ny/2] << endl;
		//	if(i%100 == 0) system("PAUSE");
		//}
		
	for(int i = 0; i < 2*(dataTM->gridX); i++)
	{
		//Ez_inc[i] = 0.0;
		//Hy_inc[i] = 0.0;
		//Hx_inc[i] = 0.0;
		H_inc[i] = 0.0;
		E_inc[i] = 0.0;
	}

		


}


void CSolver::InitSimpleFDTD_TE_UPML(void)
{
	//Init timestep parameters
	S = 1;
	S2D = S / sqrt(2.0);
	Nx = dataTE->gridX;
	Ny = dataTE->gridY;
	//dx = 1.0e-10;
	//dy = 1.0e-10;
	dz = 0.0;
	stepLambda = 240.0;
	dx = LAMBDA/stepLambda;
	dy = LAMBDA/stepLambda;
	dxy = sqrt(dx*dx+ dy*dy);
	dt = S* dx / (2* CC);

	if (dt > 1.0/(2.0*FREQ))
	{cout << "Warning::Computed dt is greater than Nyquist: " << dt << endl; system("PAUSE");}

	//First order MUr Update Constants
	MurConstY = (CC*dt-dy)/(CC*dt+dy);
	MurConstX = (CC*dt-dx)/(CC*dt+dx);

	//Init  source
	source_angle = 0;
	m0 = 29;
	lTF =30;
	rTF =390;
	bTF =35;
	uTF = 205;
	int ia;
	int ja;
	int ib;
	int jb;
	int mod;
	ia = lTF;
	ja = bTF;
	ib = rTF;
	jb = uTF;

	if ((source_angle >= 0) && (source_angle <= 90))
	{
		Oja = ja;
		Oia = ia;
	}
	else if ((source_angle >=90) && (source_angle <=180))
	{
		Oja = ja;
		Oia = ib;
	}
	else if ((source_angle >=180) && (source_angle <=270))
	{
		Oja = jb;
		Oia = ib;
	}
	else
	{
		Oja = jb;
		Oia = ia;
	}

	// Make Sure everythings adds up nicley
	mod=source_angle % 90;
	if (mod==0){
		if (source_angle == 0)
		{
			CosSource_Angle = 1;
			SinSource_Angle = 0;
		}
		else if (source_angle == 90)
		{
			CosSource_Angle = 0;
			SinSource_Angle = 1;
		}
		else if (source_angle == 180)
		{
			CosSource_Angle = -1;
			SinSource_Angle = 0;
		}
		else if (source_angle == 270)
		{
			CosSource_Angle = 0;
			SinSource_Angle = -1;
		}
	}
	else
	{
		CosSource_Angle = cos(source_angle*M_PI/180.0);
		SinSource_Angle = sin(source_angle*M_PI/180.0);
	}

	//theta = 90*M_PI/180;
 //   A =M_PI/Nx*CosSource_Angle*sin(theta);
 //   B =M_PI/Ny*SinSource_Angle*sin(theta);
 //   C =M_PI/Ny*SinSource_Angle*cos(theta);
 //   D =sin(M_PI*S2D/Nx)*sin(M_PI*S2D/Nx);
	S2DS1D = sqrt((CosSource_Angle*CosSource_Angle*CosSource_Angle*CosSource_Angle
		+SinSource_Angle*SinSource_Angle*SinSource_Angle*SinSource_Angle));

	dx1=S2DS1D*dx;
	ConstIncE = dt/dx1/EPS_Z;
	ConstIncH = dt/dx1/MU_Z;
	Offset = sqrt((ia*CosSource_Angle)*(ia*CosSource_Angle)+
		(ja*SinSource_Angle)*(ja*SinSource_Angle));


	emiTime = -10000000;

	//Init MEDIA
	
	double rcyl;
	double MinEps = 1.0;
	double omegasqr;
	double gammasqr;
	double zeps;
	
	acyl = 54;
	bcyl = 102;

	LBslab = 200;
	RBslab = 320;
	for(int i = 0; i < dataTE->gridX; i++)
		for(int j = 0; j < dataTE->gridY; j++)
		{



			gammaE[i][j] = 0.0;
			gammaM[i][j] = 0.0;
			omegaPE[i][j] = 0.0;
			omegaPM[i][j] = 0.0;

			Eps[i][j] = 1;
			Mju[i][j] = 1;
			Sig[i][j] = 0;
			SigS[i][j] = 0;


			if( (j>= LBslab) && (j<=RBslab) && (i>=140) && (i<=dataTE->gridX-140))
			{
				omegaPE[i][j] = sqrt(2.0)*OMEGA;
				omegaPM[i][j] = sqrt(2.0)*OMEGA;
				gammaE[i][j] = 0;//0.0005*OMEGA;
				gammaM[i][j] = 0;//0.0005*OMEGA;
			}



			Gz[i][j] = 0.0;
			Fz[i][j] = 0.0;
			Gzeven[i][j] = 0.0;
			Fzeven[i][j] = 0.0;
			Gzodd[i][j] = 0.0;
			Fzodd[i][j] = 0.0;

			Bx[i][j] = 0.0;
			By[i][j] = 0.0;
			Tx[i][j] = 0.0;
			Ty[i][j] = 0.0;

			Bxeven[i][j] = 0.0;
			Byeven[i][j] = 0.0;
			Txeven[i][j] = 0.0;
			Tyeven[i][j] = 0.0;

			Bxodd[i][j] = 0.0;
			Byodd[i][j] = 0.0;
			Txodd[i][j] = 0.0;
			Tyodd[i][j] = 0.0;
		}

	//Init UPML
		double mu_r_x_1, mu_r_x_2, eps_r_x_1,eps_r_x_2;
		
		double boundaryWidth, boundaryHeight;

		mu_r_x_1 = 1;
		mu_r_x_2 = 1;
		eps_r_x_1 = 1;
		eps_r_x_2 = 1;
		pmlUD = 20;
		pmlLR = 20;

		boundaryWidth = (double)pmlLR*dx;
		boundaryHeight = (double)pmlUD*dy;
		

		ka_max = 1.0;
		exponent = 6;
		R_err = 1e-16;
		double eta_1 = sqrt(MU_Z*mu_r_x_1/EPS_Z/eps_r_x_1);
		double eta_2 = sqrt(MU_Z*mu_r_x_2/EPS_Z/eps_r_x_2);
		sigma_max_1= -log(R_err) * (exponent + 1.0) / (2.0 * IMP_Z * boundaryWidth);
		sigma_max_2= -log(R_err) * (exponent + 1.0) / (2.0 * IMP_Z * boundaryHeight);;
		boundaryFactor = sigma_max_1 / ( dx * (pow(boundaryWidth,exponent)) * (exponent + 1));

		tw = 26.53e-12;
		t0 = 4*tw;

	//Init FDTDSolver Consts
		double z = 0;
		


	for(int i = 0; i < dataTE->gridX; i++)
		for(int j = 0; j < dataTE->gridY; j++)
		{
			tmpI = i;
			tmpJ = j;
			Cb[i][j] = (dt / (Eps[i][j]* EPS_Z)) / (1.0 + (Sig[i][j] * dt)/(2 * Eps[i][j]* EPS_Z));

			Db[i][j] = (dt / (Mju[i][j]* MU_Z)) / (1.0 + (SigS[i][j] * dt)/(2 * Mju[i][j]* MU_Z));



			CaGz[i][j] = (2*EPS_Z - Sig_x(i,false)*dt)/
				(2*EPS_Z + Sig_x(i,false)*dt);
			CbGz[i][j] = 2*EPS_Z*dt/(2*EPS_Z + Sig_x(i,false)*dt);

			CaEz[i][j] = (2*EPS_Z - Sig_y(j,false)*dt)/
				(2*EPS_Z + Sig_y(j,false)*dt);
			CbEz[i][j] = (2*EPS_Z * K_z(z) + Sig_z(z)*dt)/
				((2*EPS_Z + Sig_y(j,false)*dt)*MU_Z);
			CcEz[i][j] = (2*EPS_Z*K_z(z)-Sig_z(z)*dt)/
				((2*EPS_Z + Sig_y(j,false)*dt)*MU_Z);

			DaBx[i][j] = (2*EPS_Z - Sig_y(j,true)*dt)/
				(2*EPS_Z + Sig_y(j,true)*dt);
			DbBx[i][j] = 2*EPS_Z*dt/(2*EPS_Z + Sig_y(j,true)*dt);
			DaHx[i][j] = (2*EPS_Z*K_z(z) - Sig_z(z)*dt)/
				(2*EPS_Z*K_z(z) + Sig_z(z)*dt);
			DbHx[i][j] = (2*EPS_Z + Sig_x(i,false)*dt)/
				((2*EPS_Z*K_z(z) + Sig_z(z)*dt)*EPS_Z);
			DcHx[i][j] = (2*EPS_Z - Sig_x(i,false)*dt)/
				((2*EPS_Z*K_z(z) + Sig_z(z)*dt)*EPS_Z);


			DaBy[i][j] = (2*EPS_Z*K_z(z) - Sig_z(z)*dt)/
				(2*EPS_Z*K_z(z) - Sig_z(z)*dt);
			DbBy[i][j] = 2*EPS_Z*dt/(2*EPS_Z*K_z(z) + Sig_z(z)*dt);
			DaHy[i][j] = (2*EPS_Z*K_z(i) - Sig_x(i,true)*dt)/
				(2*EPS_Z*K_x(i) + Sig_x(i,true)*dt);
			DbHy[i][j] = (2*EPS_Z*K_y(j) + Sig_y(j,false)*dt)/
				((2*EPS_Z*K_x(i) + Sig_x(i,true)*dt)*EPS_Z);
			DcHy[i][j] = (2*EPS_Z*K_y(j) - Sig_y(j,false)*dt)/
				((2*EPS_Z*K_x(i) + Sig_x(i,true)*dt)*EPS_Z);


			//==========================================================

			//double plug = 1.0;
			//CaF[i][j] = (1.0/( dt*dt*plug)) + (gammaE[i][j]/(2.0* dt*plug));
			//CbF[i][j] = 2.0/( dt*dt*plug);
			//CcF[i][j] = (2.0/(dt*dt)) - omegaPE[i][j]*omegaPE[i][j]/2.0;
			//CdF[i][j] = (1.0/(dt*dt)) - (gammaE[i][j]/(2*dt)) + omegaPE[i][j]*omegaPE[i][j]/4.0;
			//CeF[i][j] = (1.0/( dt*dt*plug)) - (gammaE[i][j]/(2.0* dt*plug));
			//CfF[i][j] = (1.0/(dt*dt)) + (gammaE[i][j]/(2*dt)) + omegaPE[i][j]*omegaPE[i][j]/4.0;

			//CaT[i][j] = (1.0/( dt*dt*plug)) + (gammaM[i][j]/(2.0* dt*plug));
			//CbT[i][j] = 2.0/( dt*dt*plug);
			//CcT[i][j] = (2.0/(dt*dt)) - omegaPM[i][j]*omegaPM[i][j]/2.0;
			//CdT[i][j] = (1.0/(dt*dt)) - (gammaM[i][j]/(2*dt)) + omegaPM[i][j]*omegaPM[i][j]/4.0;
			//CeT[i][j] = (1.0/( dt*dt*plug)) - (gammaM[i][j]/(2.0* dt*plug));
			//CfT[i][j] = (1.0/(dt*dt)) + (gammaM[i][j]/(2*dt)) + omegaPM[i][j]*omegaPM[i][j]/4.0;
			double plug = 1.0;

			if((i != 0) && (j != 0))
			{
				CaFy[i][j] = (1.0/(plug)) + (dt*((gammaE[i][j-1]+gammaE[i-1][j-1])/2.0)/(2.0*plug));
				CbFy[i][j] = 2.0/( plug);
				CcFy[i][j] = ((2.0/1.0) - dt*dt*(pow(omegaPE[i][j-1],2)+pow(omegaPE[i-1][j-1], 2))/4.0);
				CdFy[i][j] = ((1.0/1.0) - (dt*((gammaE[i][j-1]+gammaE[i-1][j-1])/2.0)/(2.0)) + dt*dt*(pow(omegaPE[i][j-1],2)+pow(omegaPE[i-1][j-1], 2))/8.0);
				CeFy[i][j] = (1.0/( plug)) - (dt*((gammaE[i][j-1]+gammaE[i-1][j-1])/2.0)/(2.0*plug));
				CfFy[i][j] = ((1.0) + (dt*((gammaE[i][j-1]+gammaE[i-1][j-1])/2.0)/(2.0)) + dt*dt*(pow(omegaPE[i][j-1],2)+pow(omegaPE[i-1][j-1], 2))/8.0);
			}
			else
			{
				CaFy[i][j] = (1.0/(plug)) + (dt*((gammaE[i][j])/1.0)/(2.0*plug));
				CbFy[i][j] = 2.0/( plug);
				CcFy[i][j] = ((2.0/1.0) - dt*dt*pow((omegaPE[i][j])/1.0, 2)/2.0);
				CdFy[i][j] = ((1.0/1.0) - (dt*((gammaE[i][j])/1.0)/(2.0)) + dt*dt*pow((omegaPE[i][j])/1.0, 2)/4.0);
				CeFy[i][j] = (1.0/( plug)) - (dt*((gammaE[i][j])/1.0)/(2.0*plug));
				CfFy[i][j] = ((1.0) + (dt*((gammaE[i][j])/1.0)/(2.0)) + dt*dt*pow((omegaPE[i][j])/1.0, 2)/4.0);			
			}

			if((j >= 1))
			{
				CaFx[i][j] = (1.0/(plug)) + (dt*((gammaE[i][j-1]+gammaE[i][j])/2.0)/(2.0*plug));
				CbFx[i][j] = 2.0/( plug);
				CcFx[i][j] = ((2.0/1.0) - dt*dt*(pow(omegaPE[i][j-1],2)+pow(omegaPE[i][j], 2))/4.0);
				CdFx[i][j] = ((1.0/1.0) - (dt*((gammaE[i][j-1]+gammaE[i][j])/2.0)/(2.0)) + dt*dt*(pow(omegaPE[i][j-1],2)+pow(omegaPE[i][j], 2))/8.0);
				CeFx[i][j] = (1.0/( plug)) - (dt*((gammaE[i][j-1]+gammaE[i][j])/2.0)/(2.0*plug));
				CfFx[i][j] = ((1.0) + (dt*((gammaE[i][j-1]+gammaE[i][j])/2.0)/(2.0)) + dt*dt*(pow(omegaPE[i][j-1],2)+pow(omegaPE[i][j], 2))/8.0);
			}
			else
			{
				CaFx[i][j] = (1.0/(plug)) + (dt*((gammaE[i][j])/1.0)/(2.0*plug));
				CbFx[i][j] = 2.0/( plug);
				CcFx[i][j] = ((2.0/1.0) - dt*dt*pow((omegaPE[i][j])/1.0, 2)/2.0);
				CdFx[i][j] = ((1.0/1.0) - (dt*((gammaE[i][j])/1.0)/(2.0)) + dt*dt*pow((omegaPE[i][j])/1.0, 2)/4.0);
				CeFx[i][j] = (1.0/( plug)) - (dt*((gammaE[i][j])/1.0)/(2.0*plug));
				CfFx[i][j] = ((1.0) + (dt*((gammaE[i][j])/1.0)/(2.0)) + dt*dt*pow((omegaPE[i][j])/1.0, 2)/4.0);			
			}

			if((j >= 1))
			{
				CaT[i][j] = (1.0/(plug)) + (dt*gammaM[i][j-1]/(2.0*plug));
				CbT[i][j] = 2.0/( plug);
				CcT[i][j] = ((2.0/1.0) - dt*dt*omegaPM[i][j-1]*omegaPM[i][j-1]/2.0);
				CdT[i][j] = ((1.0/1.0) - (dt*gammaM[i][j-1]/(2.0)) + dt*dt*omegaPM[i][j-1]*omegaPM[i][j-1]/4.0);
				CeT[i][j] = (1.0/( plug)) - (dt*gammaM[i][j-1]/(2.0*plug));
				CfT[i][j] = ((1.0) + (dt*gammaM[i][j-1]/(2.0)) + dt*dt*omegaPM[i][j-1]*omegaPM[i][j-1]/4.0);
			}
			else
			{
				CaT[i][j] = (1.0/(plug)) + (dt*gammaM[i][j]/(2.0*plug));
				CbT[i][j] = 2.0/( plug);
				CcT[i][j] = ((2.0/1.0) - dt*dt*omegaPM[i][j]*omegaPM[i][j]/2.0);
				CdT[i][j] = ((1.0/1.0) - (dt*gammaM[i][j]/(2.0)) + dt*dt*omegaPM[i][j]*omegaPM[i][j]/4.0);
				CeT[i][j] = (1.0/( plug)) - (dt*gammaM[i][j]/(2.0*plug));
				CfT[i][j] = ((1.0) + (dt*gammaM[i][j]/(2.0)) + dt*dt*omegaPM[i][j]*omegaPM[i][j]/4.0);
			}

			//double plug = 1.0;
			//CaF[i][j] = 1;//(1.0/(plug)) + (dt*gammaE[i][j]/(2.0*plug));
			//CbF[i][j] = 2.0/( plug);
			//CcF[i][j] = 2;//((2.0/1.0) - dt*dt*omegaPE[i][j]*omegaPE[i][j]/2.0);
			//CdF[i][j] = 1;//((1.0/1.0) - (dt*gammaE[i][j]/(2.0)) + dt*dt*omegaPE[i][j]*omegaPE[i][j]/4.0);
			//CeF[i][j] = 1;//(1.0/( plug)) - (dt*gammaE[i][j]/(2.0*plug));
			//CfF[i][j] = 1;//((1.0) + (dt*gammaE[i][j]/(2.0)) + dt*dt*omegaPE[i][j]*omegaPE[i][j]/4.0);

			//CaT[i][j] = (1.0/(plug)) + (dt*gammaM[i][j]/(2.0*plug));
			//CbT[i][j] = 2.0/( plug);
			//CcT[i][j] = ((2.0/1.0) - dt*dt*omegaPM[i][j]*omegaPM[i][j]/2.0);
			//CdT[i][j] = ((1.0/1.0) - (dt*gammaM[i][j]/(2.0)) + dt*dt*omegaPM[i][j]*omegaPM[i][j]/4.0);
			//CeT[i][j] = (1.0/( plug)) - (dt*gammaM[i][j]/(2.0*plug));
			//CfT[i][j] = ((1.0) + (dt*gammaM[i][j]/(2.0)) + dt*dt*omegaPM[i][j]*omegaPM[i][j]/4.0);
			

		}

		//for(int i = dataTM->gridY/2; i < dataTM->gridX; i++)
		//for(int j = 0; j < dataTM->gridY/2; j++)
		//{


		//	CaFz[i][j] = CaFz[dataTM->gridX-1-i][dataTM->gridY-1-j];
		//	CbFz[i][j] =CbFz[dataTM->gridX-1-i][dataTM->gridY-1-j];

		//	CaGz[i][j] = CaGz[dataTM->gridX-1-i][dataTM->gridY-1-j];
		//	CbGz[i][j] = CbGz[dataTM->gridX-1-i][dataTM->gridY-1-j];

		//	CaEz[i][j] = CaEz[dataTM->gridX-1-i][dataTM->gridY-1-j];
		//	CbEz[i][j] = CbEz[dataTM->gridX-1-i][dataTM->gridY-1-j];
		//	CcEz[i][j] = CcEz[dataTM->gridX-1-i][dataTM->gridY-1-j];

		//	DaBx[i][j] = DaBx[dataTM->gridX-1-i][dataTM->gridY-1-j];
		//	DbBx[i][j] = DbBx[dataTM->gridX-1-i][dataTM->gridY-1-j];
		//	DaHx[i][j] = DaHx[dataTM->gridX-1-i][dataTM->gridY-1-j];
		//	DbHx[i][j] = DbHx[dataTM->gridX-1-i][dataTM->gridY-1-j];
		//	DcHx[i][j] = DcHx[dataTM->gridX-1-i][dataTM->gridY-1-j];


		//	DaBy[i][j] = DaBy[dataTM->gridX-1-i][dataTM->gridY-1-j];
		//	DbBy[i][j] = DbBy[dataTM->gridX-1-i][dataTM->gridY-1-j];
		//	DaHy[i][j] = DaHy[dataTM->gridX-1-i][dataTM->gridY-1-j];
		//	DbHy[i][j] = DbHy[dataTM->gridX-1-i][dataTM->gridY-1-j];
		//	DcHy[i][j] = DcHy[dataTM->gridX-1-i][dataTM->gridY-1-j];
		//	

		//}
	//TF/SF init
		//Ez_inc = new double[dataTM->gridX];
		//Hy_inc = new double[dataTM->gridX];
		//Hx_inc = new double[dataTM->gridX];
		H_inc = new double[2*(dataTE->gridX)];
		E_inc = new double[2*(dataTE->gridX)];

		//for(int i = 0; i < dataTM->gridX; i++)
		//{
		//	cout <<i <<": " << Eps[i][Ny/2] << endl;
		//	if(i%100 == 0) system("PAUSE");
		//}
		
	for(int i = 0; i < 2*(dataTE->gridX); i++)
	{
		//Ez_inc[i] = 0.0;
		//Hy_inc[i] = 0.0;
		//Hx_inc[i] = 0.0;
		H_inc[i] = 0.0;
		E_inc[i] = 0.0;
	}

		


}

void CSolver::SimpleFDTDSolver_TM_UPML(int timestep)
{
	double rbc, rbc1, tmp1, tmp2,Fztmp,Gztmp,Tx_pr,Ty_pr,Txtmp,Tytmp,Bxtmp,Bytmp;
	//TFSF part
	
	if(timestep<emiTime)
	{

		rbc1 = E_inc[1];
		rbc = E_inc[2*(dataTM->gridX)-2];
		for(int i = 1; i <2*(dataTM->gridX)-1; i++)
		{
			//Ey_inc[i] = (Ey_inc[i] - (dt/(EPS_Z*dx)) * (Hz_inc[i] - Hz_inc[i-1]))*cos(M_PI/4.0);
			//Ex_inc[i] = (Ex_inc[i] - (dt/(EPS_Z*dx)) * (Hz_inc[i] - Hz_inc[i-1]))*sin(M_PI/4.0);
			E_inc[i] = E_inc[i] + ConstIncE * (H_inc[i-1] - H_inc[i]);

		}
		E_inc[0] = rbc1;
		E_inc[2*(dataTM->gridX)-1]=rbc;
		

		 E_inc[m0]=cos(OMEGA*dt*timestep);
	}
	//if(timestep==10)dataTM->Ez[60][60]=100;//sin(OMEGA*dt*timestep);
	//Update E - Fields

	for(int i = 1; i < (dataTM->gridX - 1); i++)
		for(int j = 1; j < (dataTM->gridY - 1); j++)
		{
			
			
			Fz_pr = Fz[i][j];
			Gz_pr = Gz[i][j];
			if(timestep%2 == 0) 
			{
				Fztmp = Fzodd[i][j];
				Gztmp = Gzodd[i][j];
				Fzeven[i][j] = Fz[i][j];
				Gzeven[i][j] = Gz[i][j];

			}
			else
			{
				Fztmp = Fzeven[i][j];
				Gztmp = Gzeven[i][j];
				Fzodd[i][j] = Fz[i][j];
				Gzodd[i][j] = Gz[i][j];
				
			}
			Fz[i][j] =  CaGz[i][j]*Fz[i][j] + (CbGz[i][j])*((dataTM->Hy[i][j] - dataTM->Hy[i-1][j])/dx + 
					(dataTM->Hx[i][j-1] - dataTM->Hx[i][j])/dy);

			Gz[i][j] = (CaT[i][j]*Fz[i][j] - CbT[i][j]*Fz_pr + CcT[i][j]*Gz_pr -
				CdT[i][j]*Gztmp + CeT[i][j]*Fztmp)/CfT[i][j];
			//if((i==113)&&(j==120))cout<<endl << (2*Gz_pr-Gztmp+Fz[i][j] -2*Fz_pr+Fztmp) << endl;

			dataTM->Ez[i][j] = CaEz[i][j]*dataTM->Ez[i][j] + CbEz[i][j]*Gz[i][j] -
				CcEz[i][j]*Gz_pr;
			if((fabs(dataTM->Ez[i][j])>10000))
			{
				cout << i << "/" << j << endl;
				system("PAUSE");
			}


			

	/*dataTM->Ez[i][j] = Ca[i][j] * dataTM->Ez[i][j] + Cb[i][j] * ((dataTM->Hy[i][j] - dataTM->Hy[i-1][j])/dx + 
		(dataTM->Hx[i][j-1] - dataTM->Hx[i][j])/dy);*/
	
		}
	//Correct TF
		if(timestep<=emiTime)
		for(int j = bTF; j <= uTF; j++)
		{
			tmp1 = interp1(H_inc,PlaneWaveLookUpd(lTF-1,j));
			dataTM->Ez[lTF][j] = dataTM->Ez[lTF][j] + (Cb[lTF][j]/dy) * interp1(H_inc,PlaneWaveLookUpd(lTF-1,j))*CosSource_Angle;
			dataTM->Ez[rTF][j] = dataTM->Ez[rTF][j] - (Cb[rTF][j]/dy) * interp1(H_inc,PlaneWaveLookUpd(rTF,j))*CosSource_Angle;
			//tmp1 = interp1(H_inc,PlaneWaveLookUpd(lTF-1,j));
			//Fz[lTF][j] = Fz[lTF][j] + (CbFz[lTF][j]/dy) * interp1(H_inc,PlaneWaveLookUpd(lTF-1,j))*CosSource_Angle;
			//Fz[rTF][j] = Fz[rTF][j] - (CbFz[rTF][j]/dy) * interp1(H_inc,PlaneWaveLookUpd(rTF,j))*CosSource_Angle;
			//Gz[lTF][j] = Gz[lTF][j] + (CbGz[lTF][j]/dy) * interp1(H_inc,PlaneWaveLookUpd(lTF-1,j))*CosSource_Angle;
			//Gz[rTF][j] = Gz[rTF][j] - (CbGz[rTF][j]/dy) * interp1(H_inc,PlaneWaveLookUpd(rTF,j))*CosSource_Angle;
		}
		if(timestep<=emiTime)
		for(int i = lTF; i <= rTF; i++)
		{
			tmp1 = interp1(H_inc,PlaneWaveLookUpd(i,bTF-1));
	        dataTM->Ez[i][bTF] = dataTM->Ez[i][bTF] + (Cb[i][bTF]/dx)* interp1(H_inc,PlaneWaveLookUpd(i,bTF-1))*SinSource_Angle;
			dataTM->Ez[i][uTF] = dataTM->Ez[i][uTF] - (Cb[i][uTF]/dx)* interp1(H_inc,PlaneWaveLookUpd(i,uTF))*SinSource_Angle;
			//tmp1 = interp1(H_inc,PlaneWaveLookUpd(i,bTF-1));
	  //      Fz[i][bTF] = Fz[i][bTF] + (CbFz[i][bTF]/dx)* interp1(H_inc,PlaneWaveLookUpd(i,bTF-1))*SinSource_Angle;
			//Fz[i][uTF] = Fz[i][uTF] - (CbFz[i][uTF]/dx)* interp1(H_inc,PlaneWaveLookUpd(i,uTF))*SinSource_Angle;
			//Gz[i][bTF] = Gz[i][bTF] + (CbGz[i][bTF]/dx)* interp1(H_inc,PlaneWaveLookUpd(i,bTF-1))*SinSource_Angle;
			//Gz[i][uTF] = Gz[i][uTF] - (CbGz[i][uTF]/dx)* interp1(H_inc,PlaneWaveLookUpd(i,uTF))*SinSource_Angle;
		}

	//=======
	if(timestep<=emiTime)
	for(int i = 0; i < (2*(dataTM->gridX) - 1); i++)
		H_inc[i] = H_inc[i] + (ConstIncH) * (E_inc[i] - E_inc[i+1]);
//======================================================================================
    int center = dataTM->gridY/2;
	int sourceWidth = 120;

	int t0=2500;
	double A = 1;
	double B = 1;
	double lineK =tan(23*M_PI/24);
	int x_st = 30;
	int x_en = 500;
	int y_st = 100;
	double lineB = y_st - lineK*x_st;

	for(int i = x_st; i <= x_en; i++)
	{
		double y_en = i*lineK + lineB;
		int j = static_cast< int >(y_en);
		//cout << j << endl;
		dataTM->Ez[i][j] = sin(OMEGA*dt*timestep);
	}
	//dataTM->Ez[268][center+7] = sin(OMEGA*dt*timestep);
	//dataTM->Ez[268][center-7] = sin(OMEGA*dt*timestep);

	//for(int j = center-sourceWidth+1; j < center+sourceWidth; j++)
	//{
	//	A = 25.0*exp(-((j-center)/(((double)sourceWidth/50.0)*sqrt((double)sourceWidth)))*((j-center)/(((double)sourceWidth/50.0)*sqrt((double)sourceWidth))));
 //       B = exp(-(double)(timestep-t0)*(timestep-t0)/((double)(t0*t0)/10.0));
	//	dataTM->Ez[25][j] = A*sin(OMEGA*dt*timestep) + 0.1*dataTM->Ez[25][j];
	//}
	//
	//
	//for(int i = 25; i <=35; i++)
	//{
	//	dataTM->Ez[i][center-sourceWidth] = 0;
	//	dataTM->Ez[i][center+sourceWidth] = 0;
	//}

	//for(int j = 20; j <=center-sourceWidth; j++)
	//{
	//	dataTM->Ez[35][j] = 0;
	//	dataTM->Ez[35][j] = 0;
	//}

	//for(int j = center+sourceWidth; j <=dataTM->gridY-20; j++)
	//{
	//	dataTM->Ez[35][j] = 0;
	//	dataTM->Ez[35][j] = 0;
	//}

	//for(int j = 20; j <=dataTM->gridY-20; j++)
	//{
	//	dataTM->Ez[395][j] = 0;
	//	dataTM->Ez[395][j] = 0;
	//}
	//===================================================================================================================
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	double todel;
	//if(timestep>9)
	//{
	//	//for(int matI = 25; matI <= 45; matI++)
	//	//	dataTM->Ez[matI][70-matI] = sin(OMEGA*dt*timestep);

	//	for(int matI = 30; matI <= 200; matI++)
	//	{
	//		todel = exp(-(double)(matI-100)*(matI-100));
	//		dataTM->Ez[25][matI] = exp(-(double)(matI-100)*(matI-100)*dy)*sin(OMEGA*dt*timestep);
	//	}

	//}
	//Update H - Field
	
	for(int i = 1; i < (dataTM->gridX - 1); i++)
		for(int j = 1; j < (dataTM->gridY - 1); j++)
		{
			
			Bx_pr = Bx[i][j];
			By_pr = By[i][j];

			Tx_pr = Tx[i][j];
			Ty_pr = Ty[i][j];

			if(timestep%2 == 0) 
			{
				Bxtmp = Bxodd[i][j];
				Bytmp = Byodd[i][j];
				Txtmp = Txodd[i][j];
				Tytmp = Tyodd[i][j];
				Bxeven[i][j] = Bx[i][j];
				Byeven[i][j] = By[i][j];

				Txeven[i][j] = Tx[i][j];
				Tyeven[i][j] = Ty[i][j];

			}
			else
			{
				Bxtmp = Bxeven[i][j];
				Bytmp = Byeven[i][j];
				Txtmp = Txeven[i][j];
				Tytmp = Tyeven[i][j];
				Bxodd[i][j] = Bx[i][j];
				Byodd[i][j] = By[i][j];

				Txodd[i][j] = Tx[i][j];
				Tyodd[i][j] = Ty[i][j];
				
			}
			Bx[i][j] = DaBx[i][j]* Bx[i][j] + DbBx[i][j]*
				(((dataTM->Ez[i][j] - dataTM->Ez[i][j+1])/dy));

			Tx[i][j] = (CaFx[i][j]*Bx[i][j] - CbFx[i][j]*Bx_pr + CcFx[i][j]*Tx_pr -
				CdFx[i][j]*Txtmp + CeFx[i][j]*Bxtmp)/CfFx[i][j];


			dataTM->Hx[i][j] = DaHx[i][j]*dataTM->Hx[i][j]+ DbHx[i][j]*Tx[i][j]-
				DcHx[i][j]*Tx_pr;

			
			By[i][j] = DaBy[i][j]* By[i][j] + DbBy[i][j]*
				(((dataTM->Ez[i+1][j] - dataTM->Ez[i][j])/dx));

			Ty[i][j] = (CaFy[i][j]*By[i][j] - CbFy[i][j]*By_pr + CcFy[i][j]*Ty_pr -
				CdFy[i][j]*Tytmp + CeFy[i][j]*Bytmp)/CfFy[i][j];

			dataTM->Hy[i][j] = DaHy[i][j]*dataTM->Hy[i][j]+ DbHy[i][j]*Ty[i][j]-
				DcHy[i][j]*Ty_pr;

			//if((fabs(dataTM->Hx[i][j])>1000))
			//{
			//	system("PAUSE");
			//}

	/*dataTM->Hx[i][j] = Da[i][j] * dataTM->Hx[i][j] + Db[i][j] * 
				(((dataTM->Ez[i][j] - dataTM->Ez[i][j+1])/dy));
	dataTM->Hy[i][j] = Da[i][j] * dataTM->Hy[i][j] + Db[i][j] * 
				(((dataTM->Ez[i+1][j] - dataTM->Ez[i][j])/dx));*/
		}

	//Correct TF
		if(timestep<=emiTime)
		for(int i = lTF; i <= rTF; i++)
		{
			dataTM->Hx[i][bTF-1] = dataTM->Hx[i][bTF-1] +  (Db[i][bTF-1]/dx) * interp1(E_inc,PlaneWaveLookUpd(i,bTF));
			dataTM->Hx[i][uTF] = dataTM->Hx[i][uTF] -  (Db[i][uTF]/dx) * interp1(E_inc,PlaneWaveLookUpd(i,uTF));
			//Bx[i][bTF-1] = Bx[i][bTF-1] +  (DbBx[i][bTF-1]/dx) * interp1(E_inc,PlaneWaveLookUpd(i,bTF));
			//Bx[i][uTF] = Bx[i][uTF] -  (DbBx[i][uTF]/dx) * interp1(E_inc,PlaneWaveLookUpd(i,uTF));
		}
		if(timestep<=emiTime)
		for(int j = bTF; j <= uTF; j++)
		{
			dataTM->Hy[lTF-1][j] = dataTM->Hy[lTF-1][j] -  (Db[lTF-1][j]/dy) *interp1(E_inc,PlaneWaveLookUpd(lTF,j));
			dataTM->Hy[rTF][j] = dataTM->Hy[rTF][j] +  (Db[rTF][j]/dy) *interp1(E_inc,PlaneWaveLookUpd(rTF,j));
			//By[lTF-1][j] = By[lTF-1][j] -  (DbBy[lTF-1][j]/dy) *interp1(E_inc,PlaneWaveLookUpd(lTF,j));
			//By[rTF][j] = By[rTF][j] +  (DbBy[rTF][j]/dy) *interp1(E_inc,PlaneWaveLookUpd(rTF,j));
		}
	
}

void CSolver::SimpleFDTDSolver_TE_UPML(int timestep)
{
	double rbc, rbc1, tmp1, tmp2,Fztmp,Gztmp,Tx_pr,Ty_pr,Txtmp,Tytmp,Bxtmp,Bytmp;
	//TFSF part
	for(int i = 1; i < (dataTE->gridX - 1); i++)
		for(int j = 1; j < (dataTE->gridY - 1); j++)
		{

			Bx_pr = Bx[i][j];
			By_pr = By[i][j];

			Tx_pr = Tx[i][j];
			Ty_pr = Ty[i][j];

			if(timestep%2 == 0) 
			{
				Bxtmp = Bxodd[i][j];
				Bytmp = Byodd[i][j];
				Txtmp = Txodd[i][j];
				Tytmp = Tyodd[i][j];
				Bxeven[i][j] = Bx[i][j];
				Byeven[i][j] = By[i][j];

				Txeven[i][j] = Tx[i][j];
				Tyeven[i][j] = Ty[i][j];

			}
			else
			{
				Bxtmp = Bxeven[i][j];
				Bytmp = Byeven[i][j];
				Txtmp = Txeven[i][j];
				Tytmp = Tyeven[i][j];
				Bxodd[i][j] = Bx[i][j];
				Byodd[i][j] = By[i][j];

				Txodd[i][j] = Tx[i][j];
				Tyodd[i][j] = Ty[i][j];
				
			}
			Bx[i][j] = DaBx[i][j]* Bx[i][j] + DbBx[i][j]*
				(((dataTE->Hz[i][j] - dataTE->Hz[i][j-1])/dy));

			Tx[i][j] = (CaFx[i][j]*Bx[i][j] - CbFx[i][j]*Bx_pr + CcFx[i][j]*Tx_pr -
				CdFx[i][j]*Txtmp + CeFx[i][j]*Bxtmp)/CfFx[i][j];


			dataTE->Ex[i][j] = DaHx[i][j]*dataTE->Ex[i][j]+ DbHx[i][j]*Tx[i][j]-
				DcHx[i][j]*Tx_pr;


				

			By[i][j] = DaBy[i][j]* By[i][j] + DbBy[i][j]*
				(((dataTE->Hz[i-1][j] - dataTE->Hz[i][j])/dx));


			Ty[i][j] = (CaFy[i][j]*By[i][j] - CbFy[i][j]*By_pr + CcFy[i][j]*Ty_pr -
				CdFy[i][j]*Tytmp + CeFy[i][j]*Bytmp)/CfFy[i][j];

			dataTE->Ey[i][j] = DaHy[i][j]*dataTE->Ey[i][j]+ DbHy[i][j]*Ty[i][j]-
				DcHy[i][j]*Ty_pr;

		}

	//if(timestep==10)dataTM->Ez[60][60]=100;//sin(OMEGA*dt*timestep);
	//Update H - Fields
		

	for(int i = 1; i < (dataTE->gridX - 1); i++)
		for(int j = 1; j < (dataTE->gridY - 1); j++)
		{
			
			Fz_pr = Fz[i][j];
			Gz_pr = Gz[i][j];
			if(timestep%2 == 0) 
			{
				Fztmp = Fzodd[i][j];
				Gztmp = Gzodd[i][j];
				Fzeven[i][j] = Fz[i][j];
				Gzeven[i][j] = Gz[i][j];

			}
			else
			{
				Fztmp = Fzeven[i][j];
				Gztmp = Gzeven[i][j];
				Fzodd[i][j] = Fz[i][j];
				Gzodd[i][j] = Gz[i][j];
				
			}
			Fz[i][j] =  CaGz[i][j]*Fz[i][j] + (CbGz[i][j])*((dataTE->Ey[i][j] - dataTE->Ey[i+1][j])/dx + 
					(dataTE->Ex[i][j+1] - dataTE->Ex[i][j])/dy);

			Gz[i][j] = (CaT[i][j]*Fz[i][j] - CbT[i][j]*Fz_pr + CcT[i][j]*Gz_pr -
				CdT[i][j]*Gztmp + CeT[i][j]*Fztmp)/CfT[i][j];
			//if((i==113)&&(j==120))cout<<endl << (2*Gz_pr-Gztmp+Fz[i][j] -2*Fz_pr+Fztmp) << endl;

			dataTE->Hz[i][j] = CaEz[i][j]*dataTE->Hz[i][j] + CbEz[i][j]*Gz[i][j] -
				CcEz[i][j]*Gz_pr;
			if((fabs(dataTE->Hz[i][j])>100))
			{
				system("PAUSE");
			}


			

	/*dataTM->Ez[i][j] = Ca[i][j] * dataTM->Ez[i][j] + Cb[i][j] * ((dataTM->Hy[i][j] - dataTM->Hy[i-1][j])/dx + 
		(dataTM->Hx[i][j-1] - dataTM->Hx[i][j])/dy);*/
	
		}
		//dataTE->Hz[355][dataTE->gridY/2+15] = sin(OMEGA*dt*timestep);
		//dataTE->Hz[355][dataTE->gridY/2-15] = sin(OMEGA*dt*timestep);
	int center = dataTE->gridX/2;
	int sourceWidth = 360;

	int t0=48000;
	double A = 1;
	double B = 1;
	double Kx = 5.0 * 2.0*M_PI/LAMBDA;
	double T0 = 48000;
	double m = 100;

	for(int j = 0; j <4; j++)
	{
		A = exp(-dx*(j+1)*Kx);
        if(timestep <=t0)
		{
			double xON = 1.0 - (m*T0 - timestep)/(m*T0);
			B = 10.0*pow(xON, 3) - 15.0*pow(xON, 4) +6.0*pow(xON, 5);
		}
		B = 1;	
		//dataTE->Hz[center][dataTE->gridY/2] = A*B*sin(OMEGA*dt*timestep) + dataTE->Hz[center][dataTE->gridY/2];
		dataTE->Hz[center+15][140+j] = A*B*sin(OMEGA*dt*timestep) + dataTE->Hz[center+15][140+j];
		dataTE->Hz[center-15][140+j] = A*B*sin(OMEGA*dt*timestep) + dataTE->Hz[center-15][140+j];
	}
	
	/*
	for(int i = 25; i <=35; i++)
	{
		dataTE->Hz[i][center-sourceWidth] = 0;
		dataTE->Hz[i][center+sourceWidth] = 0;
	}

	for(int j = 20; j <=center-sourceWidth; j++)
	{
		dataTE->Hz[35][j] = 0;
		dataTE->Hz[35][j] = 0;
	}

	for(int j = center+sourceWidth; j <=dataTE->gridY-20; j++)
	{
		dataTE->Hz[35][j] = 0;
		dataTE->Hz[35][j] = 0;
	}

	for(int j = 20; j <=dataTE->gridY-20; j++)
	{
		dataTE->Hz[395][j] = 0;
		dataTE->Hz[395][j] = 0;
	}*/

	
}


void CSolver::TimeStepLoopTE(int N)
{
	for(int t = 0; t < N; t++)
	{
		cout << "Solving Step: " << t+1 << " of " << N; 
		//SimpleFDTDSolver_TM_UPML(t);
		SimpleFDTDSolver_TE_UPML(t);
		//HardSourceMIncidence(50,50, t);
		//Simple1DFDTDSolver(t);
		//FDTD1DIncidence(0, t);
		cout << " ...done." << endl;
		//data1D->InitNetCDF(t, dx);
	   // dataTE->InitNetCDF(t, dx);
		//if(t>=250000)
		// dataTM->WriteToNetCDF(t,20,20,20,20);
		//if(t>= 2000)
		if(t>=48000)
		{
			if(t%5 == 0) 
				dataTE->WriteToBinHz(t/5,pmlLR,pmlUD,pmlLR,pmlUD);
			if(t%5 == 0) 
				dataTE->WriteToBinEx(t/5,pmlLR,pmlUD,pmlLR,pmlUD);
			if(t%5 == 0) 
				dataTE->WriteToBinEy(t/5,pmlLR,pmlUD,pmlLR,pmlUD);
		}

		//if(t==200) Debugout();
	}

}

void CSolver::TimeStepLoopTM(int N)
{
	for(int t = 0; t < N; t++)
	{
		cout << "Solving Step: " << t+1 << " of " << N; 
		//SimpleFDTDSolver_TM_UPML(t);
		SimpleFDTDSolver_TM_UPML(t);
		//HardSourceMIncidence(50,50, t);
		//Simple1DFDTDSolver(t);
		//FDTD1DIncidence(0, t);
		cout << " ...done." << endl;
		//data1D->InitNetCDF(t, dx);
	   // dataTE->InitNetCDF(t, dx);
		//if(t>=250000)
		// dataTM->WriteToNetCDF(t,20,20,20,20);
		//if(t>= 2000)
		//if(t>3000)
		if(t>=0)
		{
			if(t%5 == 0) dataTM->WriteToBinEz(t/5,20,20,20,20);
			if(t%5 == 0) dataTM->WriteToBinHx(t/5,20,20,20,20);
			if(t%5 == 0) dataTM->WriteToBinHy(t/5,20,20,20,20);
		}


		//if(t==200) Debugout();
	}

}

void CSolver::ContinueTimeStepLoop(int start, int N)
{
	dataTM->ReadFromFile(start/1000, 20,20,20,20);
	for(int t = start+1; t < N; t++)
	{
		if(t%100 ==0)cout << "Solving Step: " << t+1 << " of " << N; 
		SimpleFDTDSolver_TM_UPML(t);
		//HardSourceMIncidence(50,50, t);
		//Simple1DFDTDSolver(t);
		//FDTD1DIncidence(0, t);
		if(t%100 ==0)cout << " ...done." << endl;
		//data1D->InitNetCDF(t, dx);
	   // dataTE->InitNetCDF(t, dx);
		//if(t>=250000)
		if(t%1000 ==0) dataTM->WriteToNetCDF(t/1000,20,20,20,20);

		//if(t%2 ==0)dataTM->WriteToNetCDF(t/2,20,20,20,20);


		if(t==200) Debugout();
	}

}

void CSolver::Debugout(void)
{
	string fileName;
	stringstream s;
	s << "C:\\hz\\hzpml\\Hz[" << "Debug" << "].txt";
	fileName = s.str();
	ofstream file(fileName.c_str(), ios::out);
	double tmp;

	if(!file.is_open())
		cout << "ERROR::Cannot open file for write!";

	cout << "Writing to file... ";
	for(int i = 0; i < dataTM->gridX; i++){
		for(int j = 0; j < dataTM->gridY; j++)
		{
			tmp = Eps[i][j];
			file << tmp <<"	, ";
		}
		file << endl;
	}
	cout << dataTM->gridX*dataTM->gridY*sizeof(double) << " BYTM written succesfully" << endl;
	file.close();
}

double CSolver::source(double t)
{
	const double E0 = 1.0;
	return E0*sin(OMEGA*dt*t);
}

double CSolver::PlaneWaveLookUpd(int i, int j)
{
	return (1.0/S2DS1D)*(CosSource_Angle*(i-Oia)+SinSource_Angle*(j-Oja))+Offset;
}

double CSolver::interp1(double* y, double xi)
{
	double xi2,temp ;
	xi2 = xi - (int)xi;
	temp = (1-xi2)*y[int(xi)] + xi2*y[int(xi)+1];
	return (1-xi2)*y[int(xi)] + xi2*y[int(xi)+1]	;
}

double CSolver::Sig_x(double x,bool offset)
{
	double tmp;
	double x1;
	double x2;
	double dist;
	dist = pmlLR-x;
	//offset = !offset;

	if((x > pmlLR-1)&&(x < (gridX - pmlLR))) tmp = 0;
	else 
	{
		if(x == pmlLR-1)
		{
			dist = pmlLR-1-x;
			if(offset)
			{
				x1 = (dist + 1.0) * dx;
				x2 = (dist + 0.0) * dx;
				tmp = boundaryFactor * (pow(x1,(exponent+1)) - pow(x2,(exponent+1)));
			}
			else 
			{
				x1 = (dist + 0.5) * dx;
				tmp = boundaryFactor * pow(x1,(exponent+1))  ;
			}
		}
		if(x == (gridX - pmlLR))
		{
			dist = x - (gridX - pmlLR);
			if(!offset)
			{
				x1 = (dist + 1.0) * dx;
				x2 = (dist + 0.0) * dx;
				tmp = boundaryFactor * (pow(x1,(exponent+1)) - pow(x2,(exponent+1)));
			}
			else 
			{
				x1 = (dist + 0.5) * dx;
				tmp = boundaryFactor * pow(x1,(exponent+1))  ;
			}
		}
		if(x < pmlLR-1)
		{
			dist = pmlLR-1-x;
			if(offset)
			{
				x1 = (dist + 1.0) * dx;       // upper bounds for point i
				x2 = (dist + 0.0) * dx;       // lower bounds for point i
			}
			else
			{
				x1 = (dist + 0.5) * dx;       // upper bounds for point i
				x2 = (dist - 0.5) * dx;       // lower bounds for point i
			}
			tmp = boundaryFactor * (pow(x1,(exponent+1)) - pow(x2,(exponent+1)));   //   polynomial grading
		}
		if(x > (gridX - pmlLR)) 
		{
			dist = x - (gridX - pmlLR);
			if(!offset)
			{
				x1 = (dist + 1.0) * dx;       // upper bounds for point i
				x2 = (dist + 0.0) * dx;       // lower bounds for point i
			}
			else
			{
				x1 = (dist + 0.5) * dx;       // upper bounds for point i
				x2 = (dist - 0.5) * dx;       // lower bounds for point i
			}
			tmp = boundaryFactor * (pow(x1,(exponent+1)) - pow(x2,(exponent+1)) );   //   polynomial grading
		}
	}
		//gradientK = 1.0 + (kmax - 1.0)* (pow(x1,(gradingOrder+1)) - pow(x2,(gradingOrder+1)) );


	return tmp;
}
double CSolver::Sig_y(double y, bool offset)
{
	double tmp;
	double x1;
	double x2;
	double dist;
	dist = pmlUD-1-y;
	//offset = !offset;
	if((y > pmlUD-1)&&(y < (gridY - pmlUD))) tmp = 0;
	else 
	{
		if(y == pmlUD-1)
		{
			dist = pmlUD-1-y;
			if(offset)
			{
				x1 = (dist + 1.0) * dy;
				x2 = (dist + 0.0) * dy;
				tmp = boundaryFactor * (pow(x1,(exponent+1)) - pow(x2,(exponent+1)));
			}
			else 
			{
				x1 = (dist + 0.5) * dy;
				tmp = boundaryFactor * pow(x1,(exponent+1))  ;
			}
		}
		if(y == (gridY - pmlUD))
		{
			dist = y - (gridY - pmlUD);
			if(!offset)
			{
				x1 = (dist + 1.0) * dy;
				x2 = (dist + 0.0) * dy;
				tmp = boundaryFactor * (pow(x1,(exponent+1)) - pow(x2,(exponent+1)));
			}
			else 
			{
				x1 = (dist + 0.5) * dy;
				tmp = boundaryFactor * pow(x1,(exponent+1))  ;
			}
		}
		if(y < pmlUD-1)
		{
			dist = pmlUD-1-y;
			if(offset)
			{
				x1 = (dist + 1.0) * dy;       // upper bounds for point i
				x2 = (dist + 0.0) * dy;       // lower bounds for point i
			}
			else
			{
				x1 = (dist + 0.5) * dy;       // upper bounds for point i
				x2 = (dist - 0.5) * dy;       // lower bounds for point i
			}
			tmp = boundaryFactor * (pow(x1,(exponent+1)) - pow(x2,(exponent+1)) );   //   polynomial grading
		}
		if(y > (gridY - pmlUD)) 
		{
			dist = y - (gridY - pmlUD);
			if(!offset)
			{
				x1 = (dist + 1.0) * dy;       // upper bounds for point i
				x2 = (dist + 0.0) * dy;       // lower bounds for point i
			}
			else
			{
				x1 = (dist + 0.5) * dy;       // upper bounds for point i
				x2 = (dist - 0.5) * dy;       // lower bounds for point i
			}
			tmp = boundaryFactor * (pow(x1,(exponent+1)) - pow(x2,(exponent+1)) );   //   polynomial grading
		}
	}
		//gradientK = 1.0 + (kmax - 1.0)* (pow(x1,(gradingOrder+1)) - pow(x2,(gradingOrder+1)) );
	return tmp;
}


double CSolver::Sig_z(double z)
{
	return 0;
}

double CSolver::K_x(double x)
{
	double tmp;
	if( (x > dx*pmlLR)&&(x < dx*(gridX - pmlLR)))tmp = 1.0;
	else if(x <= dx*pmlLR) tmp = 1.0 + (ka_max - 1.0)*pow( (dx*pmlLR - x)/((double) pmlLR*dx) ,exponent);
	else if(x >= dx*(gridX - pmlLR)) tmp = 1.0 + (ka_max - 1.0)*pow( (dx*pmlLR + x - dx*gridX)/((double) pmlLR*dx) ,exponent);
	if (tmp!=1)cout << "K_X = " << tmp <<endl;
	return 1;
}
double CSolver::K_y(double y)
{
	double tmp;
	if( (y > dy*pmlUD)&&(y < dy*(gridY - pmlUD)))tmp = 1.0;
	else if(y <= dy*pmlUD) tmp = 1.0 + (ka_max - 1.0)*pow( (dy*pmlUD - y)/((double) pmlUD*dy) ,exponent);
	else if(y >= dy*(gridY - pmlUD)) tmp = (1.0 + (ka_max - 1.0)*pow( (dy*pmlUD + y - dy*gridY)/((double) pmlUD*dy) ,exponent));
	if (tmp!=1)cout << "K_Y = " << tmp <<endl;
	return 1;
}
double CSolver::K_z(double z)
{
	return 1;
}
CSolver::~CSolver(void)
{
}
