#include "stdhead.h"
#include "Solver3d.h"

using namespace std;

CSolver3d::CSolver3d(void)
{
	data = new CField(100, 100, 100);
}
CSolver3d::CSolver3d(int SizeX, int SizeY, int SizeZ)
{
	data = new CField(SizeX, SizeY, SizeZ);
		Fz = new double**[data->gridX];
		Gz = new double**[data->gridX];

		Fx = new double**[data->gridX];
		Gx = new double**[data->gridX];

		Fy = new double**[data->gridX];
		Gy = new double**[data->gridX];

		Bx = new double**[data->gridX];
		By = new double**[data->gridX];
		Bz = new double**[data->gridX];

		Cb = new double**[data->gridX];
		Db = new double**[data->gridX];

		CaFz = new double**[data->gridX];
		CbFz = new double**[data->gridX];
		CaGz = new double**[data->gridX];
		CbGz = new double**[data->gridX];
		CaEz = new double**[data->gridX];
		CbEz = new double**[data->gridX];
		CcEz = new double**[data->gridX];

		CaFx = new double**[data->gridX];
		CbFx = new double**[data->gridX];
		CaGx = new double**[data->gridX];
		CbGx = new double**[data->gridX];
		CaEx = new double**[data->gridX];
		CbEx = new double**[data->gridX];
		CcEx = new double**[data->gridX];

		CaFy = new double**[data->gridX];
		CbFy = new double**[data->gridX];
		CaGy = new double**[data->gridX];
		CbGy = new double**[data->gridX];
		CaEy = new double**[data->gridX];
		CbEy = new double**[data->gridX];
		CcEy = new double**[data->gridX];

		DaBx = new double**[data->gridX];
		DbBx = new double**[data->gridX];
		DaBy = new double**[data->gridX];
		DbBy = new double**[data->gridX];
		DaBz = new double**[data->gridX];
		DbBz = new double**[data->gridX];
		DaHx = new double**[data->gridX];
		DbHx = new double**[data->gridX];
		DcHx = new double**[data->gridX];
		DaHy = new double**[data->gridX];
		DbHy = new double**[data->gridX];
		DcHy = new double**[data->gridX];
		DaHz = new double**[data->gridX];
		DbHz = new double**[data->gridX];
		DcHz = new double**[data->gridX];

		Eps = new double**[data->gridX];
		Mju = new double**[data->gridX];
		Sig = new double**[data->gridX];
		SigS = new double**[data->gridX];

		for(int i = 0; i < data->gridX; i++)
		{
			Fz[i] = new double*[data->gridY];
			Gz[i] = new double*[data->gridY];
			Fx[i] = new double*[data->gridY];
			Gx[i] = new double*[data->gridY];
			Fy[i] = new double*[data->gridY];
			Gy[i] = new double*[data->gridY];

			Bx[i] = new double*[data->gridY];
			By[i] = new double*[data->gridY];
			Bz[i] = new double*[data->gridY];

			Cb[i] = new double*[data->gridY];
			Db[i] = new double*[data->gridY];


			CaFz[i] = new double*[data->gridY];
			CbFz[i] = new double*[data->gridY];
			CaGz[i] = new double*[data->gridY];
			CbGz[i] = new double*[data->gridY];
			CaEz[i] = new double*[data->gridY];
			CbEz[i] = new double*[data->gridY];
			CcEz[i] = new double*[data->gridY];

			CaFx[i] = new double*[data->gridY];
			CbFx[i] = new double*[data->gridY];
			CaGx[i] = new double*[data->gridY];
			CbGx[i] = new double*[data->gridY];
			CaEx[i] = new double*[data->gridY];
			CbEx[i] = new double*[data->gridY];
			CcEx[i] = new double*[data->gridY];

			CaFy[i] = new double*[data->gridY];
			CbFy[i] = new double*[data->gridY];
			CaGy[i] = new double*[data->gridY];
			CbGy[i] = new double*[data->gridY];
			CaEy[i] = new double*[data->gridY];
			CbEy[i] = new double*[data->gridY];
			CcEy[i] = new double*[data->gridY];

			DaBx[i] = new double*[data->gridY];
			DbBx[i] = new double*[data->gridY];
			DaBy[i] = new double*[data->gridY];
			DbBy[i] = new double*[data->gridY];
			DaBz[i] = new double*[data->gridY];
			DbBz[i] = new double*[data->gridY];
			DaHx[i] = new double*[data->gridY];
			DbHx[i] = new double*[data->gridY];
			DcHx[i] = new double*[data->gridY];
			DaHy[i] = new double*[data->gridY];
			DbHy[i] = new double*[data->gridY];
			DcHy[i] = new double*[data->gridY];
			DaHz[i] = new double*[data->gridY];
			DbHz[i] = new double*[data->gridY];
			DcHz[i] = new double*[data->gridY];

			Eps[i] = new double*[data->gridY];
			Mju[i] = new double*[data->gridY];
			Sig[i] = new double*[data->gridY];
			SigS[i] = new double*[data->gridY];
		}
		for(int i = 0; i < data->gridX; i++)
			for(int j = 0; j < data->gridY; j++)
			{
				Fz[i][j] = new double[data->gridZ];
				Gz[i][j] = new double[data->gridZ];
				Fx[i][j] = new double[data->gridZ];
				Gx[i][j] = new double[data->gridZ];
				Fy[i][j] = new double[data->gridZ];
				Gy[i][j] = new double[data->gridZ];

				Bx[i][j] = new double[data->gridZ];
				By[i][j] = new double[data->gridZ];
				Bz[i][j] = new double[data->gridZ];

				Cb[i][j] = new double[data->gridZ];
				Db[i][j] = new double[data->gridZ];


				CaFz[i][j] = new double[data->gridZ];
				CbFz[i][j] = new double[data->gridZ];
				CaGz[i][j] = new double[data->gridZ];
				CbGz[i][j] = new double[data->gridZ];
				CaEz[i][j] = new double[data->gridZ];
				CbEz[i][j] = new double[data->gridZ];
				CcEz[i][j] = new double[data->gridZ];

				CaFx[i][j] = new double[data->gridZ];
				CbFx[i][j] = new double[data->gridZ];
				CaGx[i][j] = new double[data->gridZ];
				CbGx[i][j] = new double[data->gridZ];
				CaEx[i][j] = new double[data->gridZ];
				CbEx[i][j] = new double[data->gridZ];
				CcEx[i][j] = new double[data->gridZ];

				CaFy[i][j] = new double[data->gridZ];
				CbFy[i][j] = new double[data->gridZ];
				CaGy[i][j] = new double[data->gridZ];
				CbGy[i][j] = new double[data->gridZ];
				CaEy[i][j] = new double[data->gridZ];
				CbEy[i][j] = new double[data->gridZ];
				CcEy[i][j] = new double[data->gridZ];

				DaBx[i][j] = new double[data->gridZ];
				DbBx[i][j] = new double[data->gridZ];
				DaBy[i][j] = new double[data->gridZ];
				DbBy[i][j] = new double[data->gridZ];
				DaBz[i][j] = new double[data->gridZ];
				DbBz[i][j] = new double[data->gridZ];
				DaHx[i][j] = new double[data->gridZ];
				DbHx[i][j] = new double[data->gridZ];
				DcHx[i][j] = new double[data->gridZ];
				DaHy[i][j] = new double[data->gridZ];
				DbHy[i][j] = new double[data->gridZ];
				DcHy[i][j] = new double[data->gridZ];
				DaHz[i][j] = new double[data->gridZ];
				DbHz[i][j] = new double[data->gridZ];
				DcHz[i][j] = new double[data->gridZ];

				Eps[i][j] = new double[data->gridZ];
				Mju[i][j] = new double[data->gridZ];
				Sig[i][j] = new double[data->gridZ];
				SigS[i][j] = new double[data->gridZ];

			}

}





void CSolver3d::Init3DFDTD_UPML(void)
{
	//Init timestep parameters
	S = 1;
	S2D = S / sqrt(2.0);
	Nx = data->gridX;
	Ny = data->gridY;
	Nz = data->gridY;
	//dx = 0.0001;
	//dy = 0.0001;
	dz = LAMBDA/100;
	dx = LAMBDA/100;
	dy = LAMBDA/100;
	dxy = sqrt(dx*dx+ dy*dy + dz*dz);
	dt = S* dxy / (2* CC);

	if (dt > 1.0/(2.0*FREQ))
	{cout << "Warning::Computed dt is greater than Nyquist: " << dt << endl; system("PAUSE");}


	//Init  source
	source_angle = 0;
	m0 = 29;
	lTF =30;
	rTF =1010;
	bTF =30;
	uTF = 1010;
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


	emiTime = -1;

	//Init MEDIA
	double acyl;
	double bcyl;
	double rcyl;
	double MinEps = 1.0;
	acyl = 90;
	bcyl = 360;
	for(int i = 0; i < data->gridX; i++)
		for(int j = 0; j < data->gridY; j++)
			for(int k = 0; k < data->gridZ; k++)
		{



			Eps[i][j][k] = 1.0;
			Mju[i][j][k] = 1.0;
			Sig[i][j][k] = 0.0;
			SigS[i][j][k] = 0.0;


			Gz[i][j][k] = 0.0;
			Fz[i][j][k] = 0.0;
			Gx[i][j][k] = 0.0;
			Fx[i][j][k] = 0.0;
			Gy[i][j][k] = 0.0;
			Fy[i][j][k] = 0.0;
			Bx[i][j][k] = 0.0;
			By[i][j][k] = 0.0;
			Bz[i][j][k] = 0.0;
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
		


	for(int i = 0; i < data->gridX; i++)
		for(int j = 0; j < data->gridY; j++)
			for(int k = 0; k < data->gridZ; k++)
		{
			Cb[i][j][k] = (dt / (Eps[i][j][k]* EPS_Z)) / (1.0 + (Sig[i][j][k] * dt)/(2 * Eps[i][j][k]* EPS_Z));

			Db[i][j][k] = (dt / (Mju[i][j][k]* MU_Z)) / (1.0 + (SigS[i][j][k] * dt)/(2 * Mju[i][j][k]* MU_Z));


			CaFz[i][j][k] = (2*EPS_Z * Eps[i][j][k] - Sig[i][j][k]*dt)/
				(2*EPS_Z * Eps[i][j][k] + Sig[i][j][k]*dt);
			CbFz[i][j][k] = 2 * dt / (2*EPS_Z * Eps[i][j][k] + Sig[i][j][k]*dt);

			CaGz[i][j][k] = (2*EPS_Z * K_x(i) - Sig_x(i,true)*dt)/
				(2*EPS_Z * K_x(i) + Sig_x(i,true)*dt);
			CbGz[i][j][k] = 2*EPS_Z/(2*EPS_Z * K_x(i) + Sig_x(i,true)*dt);

			CaEz[i][j][k] = (2*EPS_Z*K_y(j) - Sig_y(j,true)*dt)/
				(2*EPS_Z*K_y(j) + Sig_y(j,true)*dt);
			CbEz[i][j][k] = (2*EPS_Z * K_z(k) + Sig_z(k,true)*dt)/
				(2*EPS_Z* K_y(j) + Sig_y(j,true)*dt);
			CcEz[i][j][k] = (2*EPS_Z*K_z(z)-Sig_z(k,true)*dt)/
				(2*EPS_Z*K_y(j) + Sig_y(j,true)*dt);

			CaFx[i][j][k] = (2*EPS_Z * Eps[i][j][k] - Sig[i][j][k]*dt)/
				(2*EPS_Z * Eps[i][j][k] + Sig[i][j][k]*dt);
			CbFx[i][j][k] = 2 * dt / (2*EPS_Z * Eps[i][j][k] + Sig[i][j][k]*dt);

			CaGx[i][j][k] = (2*EPS_Z * K_y(j) - Sig_y(j,true)*dt)/
				(2*EPS_Z * K_y(i) + Sig_y(j,true)*dt);
			CbGx[i][j][k] = 2*EPS_Z/(2*EPS_Z * K_y(j) + Sig_y(j,true)*dt);

			CaEx[i][j][k] = (2*EPS_Z*K_z(k) - Sig_z(k,true)*dt)/
				(2*EPS_Z*K_z(k) + Sig_z(k,true)*dt);
			CbEx[i][j][k] = (2*EPS_Z * K_x(i) + Sig_x(i,true)*dt)/
				(2*EPS_Z* K_z(k) + Sig_z(k,true)*dt);
			CcEx[i][j][k] = (2*EPS_Z*K_x(i)-Sig_x(i,true)*dt)/
				(2*EPS_Z*K_z(k) + Sig_z(k,true)*dt);

			CaFy[i][j][k] = (2*EPS_Z * Eps[i][j][k] - Sig[i][j][k]*dt)/
				(2*EPS_Z * Eps[i][j][k] + Sig[i][j][k]*dt);
			CbFy[i][j][k] = 2 * dt / (2*EPS_Z * Eps[i][j][k] + Sig[i][j][k]*dt);

			CaGy[i][j][k] = (2*EPS_Z * K_z(k) - Sig_z(k,true)*dt)/
				(2*EPS_Z * K_z(k) + Sig_z(k,true)*dt);
			CbGy[i][j][k] = 2*EPS_Z/(2*EPS_Z * K_z(k) + Sig_z(k,true)*dt);

			CaEy[i][j][k] = (2*EPS_Z*K_x(i) - Sig_x(i,true)*dt)/
				(2*EPS_Z*K_x(i) + Sig_x(i,true)*dt);
			CbEy[i][j][k] = (2*EPS_Z * K_y(j) + Sig_y(j,true)*dt)/
				(2*EPS_Z* K_x(i) + Sig_x(i,true)*dt);
			CcEy[i][j][k] = (2*EPS_Z*K_y(j)-Sig_y(j,true)*dt)/
				(2*EPS_Z*K_x(i) + Sig_x(i,true)*dt);

			DaBx[i][j][k] = (2*EPS_Z*K_y(j) - Sig_y(j,false)*dt)/
				(2*EPS_Z*K_y(j) + Sig_y(j,false)*dt);
			DbBx[i][j][k] = 2*EPS_Z*dt/(2*EPS_Z*K_y(j) + Sig_y(j,false)*dt);
			DaHx[i][j][k] = (2*EPS_Z*K_z(z) - Sig_z(z,true)*dt)/
				(2*EPS_Z*K_z(z) + Sig_z(z,true)*dt);
			DbHx[i][j][k] = (2*EPS_Z*K_x(i) + Sig_x(i,true)*dt)/
				((2*EPS_Z*K_z(z) + Sig_z(z,true)*dt)*MU_Z*Mju[i][j][k]);
			DcHx[i][j][k] = (2*EPS_Z*K_x(i) - Sig_x(i,true)*dt)/
				((2*EPS_Z*K_z(z) + Sig_z(z,true)*dt)*MU_Z*Mju[i][j][k]);


			DaBy[i][j][k] = (2*EPS_Z*K_z(z) - Sig_z(z,true)*dt)/
				(2*EPS_Z*K_z(z) - Sig_z(z,true)*dt);
			DbBy[i][j][k] = 2*EPS_Z*dt/(2*EPS_Z*K_z(z) + Sig_z(z,true)*dt);
			DaHy[i][j][k] = (2*EPS_Z*K_z(i) - Sig_x(i,false)*dt)/
				(2*EPS_Z*K_x(i) + Sig_x(i,false)*dt);
			DbHy[i][j][k] = (2*EPS_Z*K_y(j) + Sig_y(j,true)*dt)/
				((2*EPS_Z*K_x(i) + Sig_x(i,false)*dt)*MU_Z*Mju[i][j][k]);
			DcHy[i][j][k] = (2*EPS_Z*K_y(j) - Sig_y(j,true)*dt)/
				((2*EPS_Z*K_x(i) + Sig_x(i,false)*dt)*MU_Z*Mju[i][j][k]);

			DaBz[i][j][k] = (2*EPS_Z*K_x(i) - Sig_x(i,false)*dt)/
				(2*EPS_Z*K_x(i) + Sig_x(i,false)*dt);
			DbBz[i][j][k] = 2*EPS_Z*dt/(2*EPS_Z*K_x(i) + Sig_x(i,false)*dt);
			DaHz[i][j][k] = (2*EPS_Z*K_y(j) - Sig_y(j,true)*dt)/
				(2*EPS_Z*K_y(j) + Sig_y(j,true)*dt);
			DbHz[i][j][k] = (2*EPS_Z*K_z(k) + Sig_z(k,true)*dt)/
				((2*EPS_Z*K_y(j) + Sig_y(j,true)*dt)*MU_Z*Mju[i][j][k]);
			DcHz[i][j][k] = (2*EPS_Z*K_z(k) - Sig_z(k,true)*dt)/
				((2*EPS_Z*K_y(j) + Sig_y(j,true)*dt)*MU_Z*Mju[i][j][k]);
			

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
		H_inc = new double[2*(data->gridX)];
		E_inc = new double[2*(data->gridX)];

		//for(int i = 0; i < dataTM->gridX; i++)
		//{
		//	cout <<i <<": " << Eps[i][Ny/2] << endl;
		//	if(i%100 == 0) system("PAUSE");
		//}
		
	for(int i = 0; i < 2*(data->gridX); i++)
	{
		//Ez_inc[i] = 0.0;
		//Hy_inc[i] = 0.0;
		//Hx_inc[i] = 0.0;
		H_inc[i] = 0.0;
		E_inc[i] = 0.0;
	}

		


}

void CSolver3d::FDTD3DSolver_UPML(int timestep)
{
	double rbc, rbc1, tmp1, tmp2;
	//TFSF part
	
	if(timestep<emiTime)
	{

		rbc1 = E_inc[1];
		rbc = E_inc[2*(data->gridX)-2];
		for(int i = 1; i <2*(data->gridX)-1; i++)
		{
			//Ey_inc[i] = (Ey_inc[i] - (dt/(EPS_Z*dx)) * (Hz_inc[i] - Hz_inc[i-1]))*cos(M_PI/4.0);
			//Ex_inc[i] = (Ex_inc[i] - (dt/(EPS_Z*dx)) * (Hz_inc[i] - Hz_inc[i-1]))*sin(M_PI/4.0);
			E_inc[i] = E_inc[i] + ConstIncE * (H_inc[i-1] - H_inc[i]);

		}
		E_inc[0] = rbc1;
		E_inc[2*(data->gridX)-1]=rbc;
		

		 E_inc[m0]=cos(OMEGA*dt*timestep);
	}
	//if(timestep==10)dataTM->Ez[60][60]=100;//sin(OMEGA*dt*timestep);
	//Update E - Fields
	for(int i = 1; i < (data->gridX - 1); i++)
		for(int j = 1; j < (data->gridY - 1); j++)
			for(int k = 1; k < (data->gridZ - 1); k++)
		{

			Fx_pr = Fx[i][j][k];
			Gx_pr = Gx[i][j][k];
			Fx[i][j][k] =  CaFx[i][j][k]*Fx[i][j][k] + (CbFx[i][j][k])*((data->Hz[i][j][k] - data->Hz[i][j-1][k])/dy - 
					(data->Hy[i][j][k] - data->Hy[i][j][k-1])/dz);
			Gx[i][j][k] = CaGx[i][j][k]*Gx[i][j][k] + CbGx[i][j][k]*(Fx[i][j][k]-Fx_pr);
			data->Ex[i][j][k] = CaEx[i][j][k]*data->Ex[i][j][k] + CbEx[i][j][k]*Gx[i][j][k] -
				CcEx[i][j][k]*Gx_pr;

			Fy_pr = Fy[i][j][k];
			Gy_pr = Gy[i][j][k];
			Fy[i][j][k] =  CaFy[i][j][k]*Fy[i][j][k] + (CbFy[i][j][k])*((data->Hx[i][j][k] - data->Hx[i-1][j][k-1])/dz - 
					(data->Hz[i][j][k] - data->Hz[i-1][j][k])/dx);
			Gy[i][j][k] = CaGy[i][j][k]*Gy[i][j][k] + CbGy[i][j][k]*(Fy[i][j][k]-Fy_pr);
			data->Ey[i][j][k] = CaEy[i][j][k]*data->Ey[i][j][k] + CbEy[i][j][k]*Gy[i][j][k] -
				CcEy[i][j][k]*Gy_pr;
			
			Fz_pr = Fz[i][j][k];
			Gz_pr = Gz[i][j][k];
			Fz[i][j][k] =  CaFz[i][j][k]*Fz[i][j][k] + (CbFz[i][j][k]/dx)*((data->Hy[i][j][k] - data->Hy[i-1][j][k]) + 
					(data->Hx[i][j-1][k] - data->Hx[i][j][k]));
			Gz[i][j][k] = CaGz[i][j][k]*Gz[i][j][k] + CbGz[i][j][k]*(Fz[i][j][k]-Fz_pr);
			data->Ez[i][j][k] = CaEz[i][j][k]*data->Ez[i][j][k] + CbEz[i][j][k]*Gz[i][j][k] -
				CcEz[i][j][k]*Gz_pr;
			//if((fabs(dataTM->Ez[i][j])>10000))
			//{
			//	system("PAUSE");
			//}


	/*dataTM->Ez[i][j] = Ca[i][j] * dataTM->Ez[i][j] + Cb[i][j] * ((dataTM->Hy[i][j] - dataTM->Hy[i-1][j])/dx + 
		(dataTM->Hx[i][j-1] - dataTM->Hx[i][j])/dy);*/
	
		}
	//Correct TF
		//if(timestep<=emiTime)
		//for(int j = bTF; j <= uTF; j++)
		//{
		//	tmp1 = interp1(H_inc,PlaneWaveLookUpd(lTF-1,j));
		//	dataTM->Ez[lTF][j] = dataTM->Ez[lTF][j] + (Cb[lTF][j]/dy) * interp1(H_inc,PlaneWaveLookUpd(lTF-1,j))*CosSource_Angle;
		//	dataTM->Ez[rTF][j] = dataTM->Ez[rTF][j] - (Cb[rTF][j]/dy) * interp1(H_inc,PlaneWaveLookUpd(rTF,j))*CosSource_Angle;
		//	//tmp1 = interp1(H_inc,PlaneWaveLookUpd(lTF-1,j));
		//	//Fz[lTF][j] = Fz[lTF][j] + (CbFz[lTF][j]/dy) * interp1(H_inc,PlaneWaveLookUpd(lTF-1,j))*CosSource_Angle;
		//	//Fz[rTF][j] = Fz[rTF][j] - (CbFz[rTF][j]/dy) * interp1(H_inc,PlaneWaveLookUpd(rTF,j))*CosSource_Angle;
		//	//Gz[lTF][j] = Gz[lTF][j] + (CbGz[lTF][j]/dy) * interp1(H_inc,PlaneWaveLookUpd(lTF-1,j))*CosSource_Angle;
		//	//Gz[rTF][j] = Gz[rTF][j] - (CbGz[rTF][j]/dy) * interp1(H_inc,PlaneWaveLookUpd(rTF,j))*CosSource_Angle;
		//}
		//if(timestep<=emiTime)
		//for(int i = lTF; i <= rTF; i++)
		//{
		//	tmp1 = interp1(H_inc,PlaneWaveLookUpd(i,bTF-1));
	 //       dataTM->Ez[i][bTF] = dataTM->Ez[i][bTF] + (Cb[i][bTF]/dx)* interp1(H_inc,PlaneWaveLookUpd(i,bTF-1))*SinSource_Angle;
		//	dataTM->Ez[i][uTF] = dataTM->Ez[i][uTF] - (Cb[i][uTF]/dx)* interp1(H_inc,PlaneWaveLookUpd(i,uTF))*SinSource_Angle;
		//	//tmp1 = interp1(H_inc,PlaneWaveLookUpd(i,bTF-1));
	 // //      Fz[i][bTF] = Fz[i][bTF] + (CbFz[i][bTF]/dx)* interp1(H_inc,PlaneWaveLookUpd(i,bTF-1))*SinSource_Angle;
		//	//Fz[i][uTF] = Fz[i][uTF] - (CbFz[i][uTF]/dx)* interp1(H_inc,PlaneWaveLookUpd(i,uTF))*SinSource_Angle;
		//	//Gz[i][bTF] = Gz[i][bTF] + (CbGz[i][bTF]/dx)* interp1(H_inc,PlaneWaveLookUpd(i,bTF-1))*SinSource_Angle;
		//	//Gz[i][uTF] = Gz[i][uTF] - (CbGz[i][uTF]/dx)* interp1(H_inc,PlaneWaveLookUpd(i,uTF))*SinSource_Angle;
		//}

	//=======
	if(timestep<=emiTime)
	for(int i = 0; i < (2*(data->gridX) - 1); i++)
		H_inc[i] = H_inc[i] + (ConstIncH) * (E_inc[i] - E_inc[i+1]);

	data->Ez[50][50][50] = 10*sin(OMEGA*dt*timestep);
	//Update H - Field
	for(int i = 1; i < (data->gridX - 1); i++)
		for(int j = 1; j < (data->gridY - 1); j++)
			for(int k = 1; k < (data->gridZ - 1); k++)
		{
			Bx_pr = Bx[i][j][k];
			Bx[i][j][k] = DaBx[i][j][k]* Bx[i][j][k] + DbBx[i][j][k]*
				((data->Ez[i][j][k] - data->Ez[i][j+1][k])/dy + (data->Ey[i][j][k+1] - data->Ey[i][j][k])/dz);

			data->Hx[i][j][k] = DaHx[i][j][k]*data->Hx[i][j][k]+ DbHx[i][j][k]*Bx[i][j][k]-
				DcHx[i][j][k]*Bx_pr;

			By_pr = By[i][j][k];
			By[i][j][k] = DaBy[i][j][k]* By[i][j][k] + DbBy[i][j][k]*
				(((data->Ez[i+1][j][k] - data->Ez[i][j][k])/dx)+(data->Ex[i][j][k] - data->Ex[i][j][k+1])/dz);

			data->Hy[i][j][k] = DaHy[i][j][k]*data->Hy[i][j][k]+ DbHy[i][j][k]*By[i][j][k]-
				DcHy[i][j][k]*By_pr;

			Bz_pr = Bz[i][j][k];
			Bz[i][j][k] = DaBz[i][j][k]* Bz[i][j][k] + DbBz[i][j][k]*
				(((data->Ey[i][j][k] - data->Ey[i+1][j][k])/dx)+(data->Ex[i][j+1][k] - data->Ex[i][j][k])/dy);

			data->Hz[i][j][k] = DaHz[i][j][k]*data->Hz[i][j][k]+ DbHz[i][j][k]*By[i][j][k]-
				DcHz[i][j][k]*Bz_pr;



	/*dataTM->Hx[i][j] = Da[i][j] * dataTM->Hx[i][j] + Db[i][j] * 
				(((dataTM->Ez[i][j] - dataTM->Ez[i][j+1])/dy));
	dataTM->Hy[i][j] = Da[i][j] * dataTM->Hy[i][j] + Db[i][j] * 
				(((dataTM->Ez[i+1][j] - dataTM->Ez[i][j])/dx));*/
		}

	//Correct TF
		//if(timestep<=emiTime)
		//for(int i = lTF; i <= rTF; i++)
		//{
		//	dataTM->Hx[i][bTF-1] = dataTM->Hx[i][bTF-1] +  (Db[i][bTF-1]/dx) * interp1(E_inc,PlaneWaveLookUpd(i,bTF));
		//	dataTM->Hx[i][uTF] = dataTM->Hx[i][uTF] -  (Db[i][uTF]/dx) * interp1(E_inc,PlaneWaveLookUpd(i,uTF));
		//	//Bx[i][bTF-1] = Bx[i][bTF-1] +  (DbBx[i][bTF-1]/dx) * interp1(E_inc,PlaneWaveLookUpd(i,bTF));
		//	//Bx[i][uTF] = Bx[i][uTF] -  (DbBx[i][uTF]/dx) * interp1(E_inc,PlaneWaveLookUpd(i,uTF));
		//}
		//if(timestep<=emiTime)
		//for(int j = bTF; j <= uTF; j++)
		//{
		//	dataTM->Hy[lTF-1][j] = dataTM->Hy[lTF-1][j] -  (Db[lTF-1][j]/dy) *interp1(E_inc,PlaneWaveLookUpd(lTF,j));
		//	dataTM->Hy[rTF][j] = dataTM->Hy[rTF][j] +  (Db[rTF][j]/dy) *interp1(E_inc,PlaneWaveLookUpd(rTF,j));
		//	//By[lTF-1][j] = By[lTF-1][j] -  (DbBy[lTF-1][j]/dy) *interp1(E_inc,PlaneWaveLookUpd(lTF,j));
		//	//By[rTF][j] = By[rTF][j] +  (DbBy[rTF][j]/dy) *interp1(E_inc,PlaneWaveLookUpd(rTF,j));
		//}
	
}

void CSolver3d::TimeStepLoop(int N)
{
	for(int t = 0; t < N; t++)
	{
		cout << "Solving Step: " << t+1 << " of " << N; 
		FDTD3DSolver_UPML(t);
		//HardSourceMIncidence(50,50, t);
		//Simple1DFDTDSolver(t);
		//FDTD1DIncidence(0, t);
		cout << " ...done." << endl;
		//data1D->InitNetCDF(t, dx);
	   // dataTE->InitNetCDF(t, dx);
		//if(t>=250000)
		//if(t%400 ==0) dataTM->WriteToNetCDF(t/100,20,20,20,20);

		//if(t%2 ==0)dataTM->WriteToNetCDF(t/2,20,20,20,20);


		if(t==200) Debugout();
	}

}

void CSolver3d::Debugout(void)
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
	for(int i = 0; i < data->gridX; i++){
		for(int j = 0; j < data->gridY; j++)
			for(int k = 0; k < data->gridZ; k++)
		{
			tmp = Eps[i][j][k];
			file << tmp <<"	, ";
		}
		file << endl;
	}
	cout << data->gridX*data->gridY*sizeof(double) << " BYTM written succesfully" << endl;
	file.close();
}

double CSolver3d::source(double t)
{
	const double E0 = 1.0;
	return E0*sin(OMEGA*dt*t);
}

double CSolver3d::PlaneWaveLookUpd(int i, int j)
{
	return (1.0/S2DS1D)*(CosSource_Angle*(i-Oia)+SinSource_Angle*(j-Oja))+Offset;
}

double CSolver3d::interp1(double* y, double xi)
{
	double xi2,temp ;
	xi2 = xi - (int)xi;
	temp = (1-xi2)*y[int(xi)] + xi2*y[int(xi)+1];
	return (1-xi2)*y[int(xi)] + xi2*y[int(xi)+1]	;
}

double CSolver3d::Sig_x(double x,bool offset)
{
	double tmp;
	double x1;
	double x2;
	double dist;
	dist = pmlLR-x;
	//offset = !offset;

	if((x > pmlLR-1)&&(x < (data->gridX - pmlLR))) tmp = 0;
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
		if(x == (data->gridX - pmlLR))
		{
			dist = x - (data->gridX - pmlLR);
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
		if(x > (data->gridX - pmlLR)) 
		{
			dist = x - (data->gridX - pmlLR);
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
double CSolver3d::Sig_y(double y, bool offset)
{
	double tmp;
	double x1;
	double x2;
	double dist;
	dist = pmlUD-1-y;
	//offset = !offset;
	if((y > pmlUD-1)&&(y < (data->gridY - pmlUD))) tmp = 0;
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
		if(y == (data->gridY - pmlUD))
		{
			dist = y - (data->gridY - pmlUD);
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
		if(y > (data->gridY - pmlUD)) 
		{
			dist = y - (data->gridY - pmlUD);
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


double CSolver3d::Sig_z(double x, bool offset)
{
		double tmp;
	double x1;
	double x2;
	double dist;
	dist = pmlLR-x;
	//offset = !offset;

	if((x > pmlLR-1)&&(x < (data->gridX - pmlLR))) tmp = 0;
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
		if(x == (data->gridX - pmlLR))
		{
			dist = x - (data->gridX - pmlLR);
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
		if(x > (data->gridX - pmlLR)) 
		{
			dist = x - (data->gridX - pmlLR);
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

double CSolver3d::K_x(double x)
{
	double tmp;
	if( (x > dx*pmlLR)&&(x < dx*(data->gridX - pmlLR)))tmp = 1.0;
	else if(x <= dx*pmlLR) tmp = 1.0 + (ka_max - 1.0)*pow( (dx*pmlLR - x)/((double) pmlLR*dx) ,exponent);
	else if(x >= dx*(data->gridX - pmlLR)) tmp = 1.0 + (ka_max - 1.0)*pow( (dx*pmlLR + x - dx*data->gridX)/((double) pmlLR*dx) ,exponent);
	if (tmp!=1)cout << "K_X = " << tmp <<endl;
	return 1;
}
double CSolver3d::K_y(double y)
{
	double tmp;
	if( (y > dy*pmlUD)&&(y < dy*(data->gridY - pmlUD)))tmp = 1.0;
	else if(y <= dy*pmlUD) tmp = 1.0 + (ka_max - 1.0)*pow( (dy*pmlUD - y)/((double) pmlUD*dy) ,exponent);
	else if(y >= dy*(data->gridY - pmlUD)) tmp = (1.0 + (ka_max - 1.0)*pow( (dy*pmlUD + y - dy*data->gridY)/((double) pmlUD*dy) ,exponent));
	if (tmp!=1)cout << "K_Y = " << tmp <<endl;
	return 1;
}
double CSolver3d::K_z(double z)
{
	return 1;
}
CSolver3d::~CSolver3d(void)
{
}
